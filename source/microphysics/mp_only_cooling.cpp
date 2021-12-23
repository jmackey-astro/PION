///
/// \file mp_only_cooling.cc
/// \author Jonathan Mackey
/// \date 14.01.2011
///
/// This is a simple microphysics class which only has heating and
/// cooling processes (i.e. no species are tracked explicitly).
///
/// Modifications:\n
/// - 2011.01.17 JM: Fixed bugs, also wrote testing program
///   (commented out -- see the bottom of this file).
/// - 2011.04.11 JM: Refining the model a bit.  Needs more work.
/// - 2011.04.12 JM: Added WSS09 cooling function.  Cooling_flag==8 is the
///    recommended cooling function, since it self-consistently models the
///    hydrogen ioniation-heating and recombination-cooling.  The normalisation
///    is lower than the SD93 curves, possibly because the Oxygen abundance is
///    now lower than it was 15 years ago (Lodders, 2003, ApJ).
///    Removed the DO_HEATING ifdef -- now we have 5 different cooling
///    functions.
///
/// - 2011.05.10 JM: Output cooling rates only if myrank==0 for parallel (so
/// processes
///    don't fight over the file and slow down the code (by a lot!)).
/// - 2015.01.13 JM: Added some comments.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2016.06.21 JM: Temperature() threadsafe.
/// - 2018.01.25 JM: added functions to request n(H+),n(H0),n(e-)

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include <sim_params.h>
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG

#include "microphysics/mp_only_cooling.h"

#include <fstream>
#include <sstream>
using namespace std;

#define KI02 2  ///< Koyama & Inutsuka (2002) function
#define DMcC 3  ///< Dalgarno & McCray (1972) function
#define SD93_CIE 4
#define SD93_PLUS_HEATING 5
#define WSS09_CIE_PLUS_HEATING 6
#define WSS09_CIE_ONLY_COOLING 7
#define WSS09_CIE_LINE_HEAT_COOL 8

// ##################################################################
// ##################################################################

mp_only_cooling::mp_only_cooling(
    const int nv,
    const int ntr,          ///< Number of tracer variables in state vector
    const std::string *tr,  ///< List of tracer variable names.
    struct which_physics *ephys,  ///< pointer to extra physics flags.
    struct rad_sources *rsrcs     ///< radiation sources.
    ) :
    microphysics_base(nv, ntr, tr, ephys, rsrcs),
    cooling_function_SD93CIE(), Hummer94_Hrecomb(), CoolingFn(ephys->cooling),
    nv_prim(nv)
{
  if (!EP->update_erg) {
    spdlog::error(
        "{}: {}", "requested cooling microphysics but no energy update",
        DONT_CALL_ME);
  }
  Integrator_Base::Set_Nvar(1);

  //
  // Mean masses per species/atom: we assume cosmic abundances which gives
  // rho \simeq 1.40 m_p n_H.  For fully ionised gas we have
  //  n_i \simeq 1.1 n_H (from helium).
  //  n_e \simeq 1.2 n_H (from doubly ionised helium)
  //  n_tot \simeq 2.3 n_H (summing the previous lines).
  // So we define muH=1.4mp, mue=1.167mp, mui=1.273mp, mutot=0.609mp
  //
  Mu             = 1.40 * pconst.m_p();
  Mu_tot         = 0.609 * pconst.m_p();
  Mu_tot_over_kB = Mu_tot / pconst.kB();
  Mu_elec        = 1.167 * pconst.m_p();
  Mu_ion         = 1.273 * pconst.m_p();
  inv_Mu2        = 1.0 / (Mu * Mu);
  inv_Mu2_elec_H = 1.0 / (Mu_elec * Mu);

  //
  // Next select which cooling function to set up.
  //
  cooling_flag = EP->cooling;

  switch (cooling_flag) {
    case KI02:
    case DMcC:
    case SD93_CIE:
    case SD93_PLUS_HEATING:
      setup_SD93_cie();
      break;
    case WSS09_CIE_PLUS_HEATING:
    case WSS09_CIE_ONLY_COOLING:
      setup_WSS09_CIE();
      break;
    case WSS09_CIE_LINE_HEAT_COOL:
      spdlog::debug("\tRequested fully ionized gas with WSS09 cooling at high"
                    " temperatures,\n\tand photoionized gas at nebular"
                    " temperatures, with T_eq approx 8000 K.");
      setup_WSS09_CIE_OnlyMetals();
      break;
    default:
      spdlog::error(
          "{}: {}", "Invalid cooling flag in mp_only_cooling", cooling_flag);
      break;
  }

  // generate lookup tables for WSS09_CIE_LINE_HEAT_COOL function.
  gen_mpoc_lookup_tables();

#ifndef NDEBUG
  ostringstream opfile;
  opfile << "coolingNOCHEM_" << cooling_flag << ".txt";
  ofstream outf(opfile.str().c_str());
  if (!outf.is_open()) spdlog::error("{}: {}", "couldn't open outfile", 1);
  outf << "Cooling Curve Data: Temperature(K) Rate(erg/cm^3/s) (n=1 per cc)\n";
  outf.setf(ios_base::scientific);
  outf.precision(6);
  double t = 1.0e0;
  do {
    outf << t << "\t" << Edot(2.34e-24, t) << "\n";
    t *= 1.05;
  } while (t < 1.0e10);
  outf.close();
#endif

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  MaxT_allowed = EP->MaxTemperature;
  MinT_allowed = EP->MinTemperature;
  if (MinT_allowed < 1.0 || MinT_allowed > 1.0e6) MinT_allowed = 1.0;
  if (MaxT_allowed < 1.0e2 || MaxT_allowed > 3.0e10) MaxT_allowed = 1.0e8;
    // cout <<"\t\tAllowed Temperature range: T_min="<<MinT_allowed<<"
    // T_max="; cout <<MaxT_allowed<<" (MAKE SURE THIS IS SET
    // APPROPRIATELY!\n";
#else
  spdlog::error(
      "{}: {}",
      "Please set SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE in "
      "defines/functionality_flags.h if using cooling",
      123);
#endif  // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return;
}

// ##################################################################
// ##################################################################

mp_only_cooling::~mp_only_cooling() {}

// ##################################################################
// ##################################################################

int mp_only_cooling::TimeUpdateMP(
    const pion_flt *p_in,  ///< Primitive Vector to be updated
    pion_flt *p_out,       ///< Destination Vector for updated values.
    const double dt,       ///< Time Step to advance by.
    const double g,        ///< EOS gamma.
    const int,             ///< Switch for what type of integration to use
    double *Tf             ///< final temperature.
)
{
  mp_only_cooling::rho   = p_in[RO];
  mp_only_cooling::gamma = g;
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];

  for (int v = 0; v < nv_prim; v++) {
    if (!isfinite(p_out[v])) {
      spdlog::debug("pout : {}", std::vector<double>(p_in, p_in + nv_prim));
      spdlog::error("{}: {}", "p_out[v]", v);
    }
  }

  double Eint0 = p_in[PG] / (gamma - 1.0);
  double T     = p_out[PG] * Mu_tot_over_kB / p_out[RO];
  if (T < MinT_allowed) {
    Set_Temp(p_out, MinT_allowed, gamma);
    T = MinT_allowed;
  }
  if (T > MaxT_allowed) {
    Set_Temp(p_out, MaxT_allowed, gamma);
    T = MaxT_allowed;
  }
  double Eint = Eint0;
  double tout = 0.0;

  int err = Int_Adaptive_RKCK(1, &Eint, 0.0, dt, 1.0e-2, &Eint, &tout);
  if (err) {
    spdlog::debug("p_in : {}", std::vector<double>(p_in, p_in + nv_prim));
    spdlog::debug(
        "Ein = {}, T={}, Eout={}, dt={}, tout={}", Eint0, T, Eint, dt, tout);
    spdlog::error("{}: {}", "mp_only_cooling integration failed.", err);
  }
  p_out[PG] = Eint * (gamma - 1);
  *Tf       = p_out[PG] * Mu_tot_over_kB / p_out[RO];
  if (*Tf > MaxT_allowed) {
    p_out[PG] *= MaxT_allowed / (*Tf);
    *Tf = MaxT_allowed;
  }
  else if (*Tf < MinT_allowed) {
    p_out[PG] *= MinT_allowed / (*Tf);
    *Tf = MinT_allowed;
  }
  return 0;
}

// ##################################################################
// ##################################################################

int mp_only_cooling::dPdt(
    const int nv,     ///< number of variables we are expecting.
    const double *P,  ///< Current state vector.
    double *R         ///< Rate Vector to write to.
)
{
  R[0] = Edot(rho, P[0] * (gamma - 1.0) * Mu_tot_over_kB / rho);
  return 0;
}

// ##################################################################
// ##################################################################

int mp_only_cooling::Set_Temp(
    pion_flt *p_in,   ///< primitive vector.
    const double T,   ///< temperature requested.
    const double gam  ///< eos gamma.
)
{
#ifndef NDEBUG
  //
  // check for bad input
  //
  if (T <= 0.0) {
    spdlog::error(
        "mp_only_cooling::Set_Temp() Requested negative temperature! T={}", T);
    return 1;
  }
#endif  // NDEBUG
  // cout <<"\n\t\t***** Set_Temp(): p="<<p_in[PG]<<", rho="<< p_in[RO];
  // cout <<", T="<<T<<", new p="<<p_in[RO]*T/Mu_tot_over_kB<<"\n";
  p_in[PG] = p_in[RO] * T / Mu_tot_over_kB;

  return 0;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::Temperature(
    const pion_flt *p_in,  ///< primitive vector
    const double           ///< eos gamma
)
{
  return p_in[PG] * Mu_tot_over_kB / p_in[RO];
}

// ##################################################################
// ##################################################################

double mp_only_cooling::get_n_elec(
    const pion_flt *p_in  ///< primitive state vector.
)
{
  // H and He are fully ionized
  return p_in[RO] / Mu_elec;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::get_n_Hplus(
    const pion_flt *p_in  ///< primitive state vector.
)
{
  // fully ionized hydrogen
  return p_in[RO] / Mu;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::get_n_Hneutral(
    const pion_flt *p_in  ///< primitive state vector.
)
{
  // Assume neutral fraction of 1e-12 (arbitrary)
  return 1.0e-12 * p_in[RO] / Mu;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::timescales(
    const pion_flt *p_in,  ///< Current cell.
    const double gam,      ///< EOS gamma.
    const bool tc,         ///< set to true if including cooling time.
    const bool,            ///< set to true if including recombination time.
    const bool             ///< set to true if including photo-ionsation time.
)
{
  if (!tc) return 1.0e99;

  //
  // First get temperature and internal energy density.
  //
  double Eint    = p_in[PG] / (gam - 1.0);
  double T       = p_in[PG] * Mu_tot_over_kB / p_in[RO];
  double mintime = 1.0e99;
  double rate    = 0.0;

  //
  // we skip the cooling time if the temperature is already near the
  // minimum allowed.
  //
  // We get the cooling rate per unit volume, first for the actual
  // temperature, and then for T/2, and we take the max. of this.
  //
  if (T >= 1.1 * MinT_allowed) {
    rate =
        max(fabs(Edot(p_in[RO], T)),
            fabs(Edot(p_in[RO], max(MinT_allowed, 0.5 * T))));
    //
    // Cooling time is then the energy per unit volume divided by the rate.
    //
    mintime = min(mintime, Eint / rate);
  }

  return mintime;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::Edot(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  double rate = 0.0;
  switch (cooling_flag) {
    case KI02:
      // this function returns rate as positive if cooling, so change
      // sign.
      rate = -CoolingRate(T, 0.0, rho / Mu, 0.0, 0.0);
      break;

    case SD93_CIE:
      rate = Edot_SD93CIE_cool(rho, T);
      break;

    case SD93_PLUS_HEATING:
      rate = Edot_SD93CIE_heat_cool(rho, T);
      break;

    case WSS09_CIE_ONLY_COOLING:
      rate = Edot_WSS09CIE_cool(rho, T);
      break;

    case WSS09_CIE_PLUS_HEATING:
      rate = Edot_WSS09CIE_heat_cool(rho, T);
      break;

    case WSS09_CIE_LINE_HEAT_COOL:
      rate = Edot_WSS09CIE_heat_cool_metallines(rho, T);
      break;

    default:
      spdlog::error(
          "{}: {}", "bad cooling flag in mp_only_cooling::Edot", cooling_flag);
      break;
  }
  return rate;
}

// ##################################################################
// ##################################################################

double mp_only_cooling::Edot_SD93CIE_cool(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  return -(rho * rho / Mu_elec / Mu_ion) * cooling_rate_SD93CIE(T);
}

// ##################################################################
// ##################################################################

//
// cooling from SD93-CIE plus heating assuming full ionisation of H, where
// the heating rate equals recombination rate times 5eV/ionisation.
//
double mp_only_cooling::Edot_SD93CIE_heat_cool(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  return (rho * rho)
         * (2.733e-21 * exp(-0.782991 * log(T)) / Mu_elec / Mu
            - cooling_rate_SD93CIE(T) / Mu_elec / Mu_ion);
}

// ##################################################################
// ##################################################################

// only cooling, uses WSS09-CIE cooling function.
double mp_only_cooling::Edot_WSS09CIE_cool(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  return 2e-26 * rho / Mu - (rho * rho / Mu / Mu) * cooling_rate_SD93CIE(T);
}

// ##################################################################
// ##################################################################

//
// cooling from WSS09-CIE plus heating assuming full ionisation of H, where
// the heating rate equals recombination rate times 5eV/ionisation.
//
double mp_only_cooling::Edot_WSS09CIE_heat_cool(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  return (rho * rho)
         * (2.733e-21 * exp(-0.782991 * log(T)) / Mu_elec / Mu
            - cooling_rate_SD93CIE(T) / Mu / Mu);
}

// ##################################################################
// ##################################################################

double mp_only_cooling::Edot_WSS09CIE_heat_cool_metallines(
    const double rho,  ///< mass density (g/cm3)
    const double T     ///< Temperature (K)
)
{
  // Get T-vector index
  int iT     = 0;
  double dT  = 0.0;
  size_t ihi = lt.NT - 1, ilo = 0, imid = 0;
  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (lt.T[imid] < T)
      ilo = imid;
    else
      ihi = imid;
  } while (ihi - ilo > 1);
  iT          = ilo;
  dT          = T - lt.T[iT];
  double rho2 = rho * rho, rate = 0.0;
  // Henney et al. (2009) forbidden lines (eq. A9)
  rate = -(lt.C_fbdn[iT] + dT * lt.s_C_fbdn[iT]) * rho2 * inv_Mu2_elec_H;
  // take stronger rate between Henney A9 and the WSS09 CIE rate.
  rate = min(rate, -(lt.C_cie[iT] + dT * lt.s_C_cie[iT]) * rho2 * inv_Mu2);
  // Hydrogen cooling due to recombinations and Bremsstrahlung.
  rate -= (lt.C_rrh[iT] + dT * lt.s_C_rrh[iT]) * rho2 * inv_Mu2_elec_H;
  // Helium Bremsstrahlung (using eq.5.15b from Rybicki &
  // Lightman (1978), p.162, with Z=2, n(He)/n(H)=0.1, g_B=1.2
  rate -= (lt.C_ffhe[iT] + dT * lt.s_C_ffhe[iT]) * rho2 * inv_Mu2_elec_H;
  // heating rate = 5eV*RRR*n(e-)*n(H+) = 8.01e-12*Hii_rrr(T)*n^2
  // assuming all hydrogen is ionised.
  rate += 8.01e-12 * (lt.rrhp[iT] + dT * lt.s_rrhp[iT]) * rho2 * inv_Mu2_elec_H;
  return rate;
}

// ##################################################################
// ##################################################################

void mp_only_cooling::gen_mpoc_lookup_tables()
{
  //  Start by generating the logarithmic temperature scale:
  lt.NT        = 200;
  double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
  lt.dlogT =
      (log10(EP->MaxTemperature) - log10(EP->MinTemperature)) / (lt.NT - 1);
  lt.T.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++) {
    lt.T[i] = pow(10.0, log10(EP->MinTemperature) + i * lt.dlogT);
  }

  lt.rrhp.resize(lt.NT);
  lt.C_rrh.resize(lt.NT);
  lt.C_ffhe.resize(lt.NT);
  lt.C_fbdn.resize(lt.NT);
  lt.C_cie.resize(lt.NT);
  // lt.C_dust.resize(lt.NT);

  for (size_t i = 0; i < lt.NT; i++) {
    lt.rrhp[i]   = Hii_rad_recomb_rate(lt.T[i]);
    lt.C_rrh[i]  = Hii_total_cooling(lt.T[i]);
    lt.C_ffhe[i] = 6.72e-28 * sqrt(lt.T[i]);
    lt.C_fbdn[i] =
        1.20e-22
        * exp(-33610.0 / lt.T[i] - (2180.0 * 2180.0 / lt.T[i] / lt.T[i]))
        * exp(-lt.T[i] * lt.T[i] / 5.0e10);
    lt.C_cie[i] = cooling_rate_SD93CIE(lt.T[i]);
    // lt.C_dust[i] = 1.0e-17 * exp(1.5*log(lt.T[i]/2.5e8));
  }

  lt.s_rrhp.resize(lt.NT);
  lt.s_C_rrh.resize(lt.NT);
  lt.s_C_ffhe.resize(lt.NT);
  lt.s_C_fbdn.resize(lt.NT);
  lt.s_C_cie.resize(lt.NT);
  // lt.s_C_dust.resize(lt.NT);
  for (size_t i = 0; i < lt.NT - 1; i++) {
    lt.s_rrhp[i]  = (lt.rrhp[i + 1] - lt.rrhp[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_rrh[i] = (lt.C_rrh[i + 1] - lt.C_rrh[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_ffhe[i] =
        (lt.C_ffhe[i + 1] - lt.C_ffhe[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_fbdn[i] =
        (lt.C_fbdn[i + 1] - lt.C_fbdn[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cie[i] = (lt.C_cie[i + 1] - lt.C_cie[i]) / (lt.T[i + 1] - lt.T[i]);
    // lt.s_C_dust[i] = (lt.C_dust[i+1]-lt.C_dust[i]) / (lt.T[i+1]-lt.T[i]);
  }
  lt.s_rrhp[lt.NT - 1]   = 0.0;
  lt.s_C_rrh[lt.NT - 1]  = 0.0;
  lt.s_C_ffhe[lt.NT - 1] = 0.0;
  lt.s_C_fbdn[lt.NT - 1] = 0.0;
  lt.s_C_cie[lt.NT - 1]  = 0.0;
  // lt.s_C_dust[lt.NT-1] = 0.0;
  return;
}

// ##################################################################
// ##################################################################

// int main()
// {
//  //
//  // test the only_cooling class:
//  // compile with
//  // g++ ../mp_only_cooling.cc ../cooling_SD93_cie.cc ../../global.cc
//  ../../cell_interface.cc -o a.out -lreadline
//  //
//  struct which_physics EP;
//  EP.dynamics=0;
//  EP.raytracing=0;
//  EP.cooling=4;
//  EP.chemistry=0;
//  EP.coll_ionisation=0;
//  EP.rad_recombination=0;
//  EP.phot_ionisation=0;
//  EP.update_erg=1;
//  EP.MP_timestep_limit=1;
//#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
//  EP.MinTemperature = 10.0; ///< Minimum temperature to allow in the
//  simulation. EP.MaxTemperature = 2.0e10; ///< Maximum temperature to allow in
//  the simulation.
//#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
//
//
// class mp_only_cooling MP (5,&EP);

//  double p[5];
//  p[RO] = 1.169e-25;
//  p[VX]=p[VY]=p[VZ]=0.0;
//  p[PG]=0.1;
//  double g=5.0/3.0;
//  MP.Set_Temp(p,1.0e5,g);
//  //for (int i=0;i<10;i++) {
//  //  cout <<"dens="<<p[RO]<<", cooling
//  time="<<MP.timescales(p,g,1,1,1)<<"\n";
//  //  p[RO]*=sqrt(10.0);
//  //  MP.Set_Temp(p,1.0e5,g);
//  //}

//  p[RO] = 1.169e-24;
//  MP.Set_Temp(p,1.0e10,g);
//  double t=0.0, h=1.0e7, T=MP.Temperature(p,g), p_f[5];
//  int err=0;
//  do {
//    //cout <<"t="<<t<<"\tT="<<T<<"\tt_cool=";
//    //cout <<(h=0.3*MP.timescales(p,g,1,1,1))<<"\n";
//    cout <<t<<"\t"<<T<<"\t";
//    cout <<MP.timescales(p,g,1,1,1)<<"\n";
//    h=0.3*MP.timescales(p,g,1,1,1);
//    err += MP.TimeUpdateMP(p,p,h,g,0,&T);
//    if (err) {
//      rep.printVec("p",p,2);
//      spdlog::error("{}: {}", "integration error",err);
//    }
//    t+=h;
//  } while (MP.Temperature(p,g)>1.0e4);

//  return 0;
// }
