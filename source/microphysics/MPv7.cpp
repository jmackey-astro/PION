///
/// \file MPv7.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running calculations with a simple
/// two-temperature isothermal equation of state, where T is T_low
/// when gas is neutral, and T_high when ionised, and interpolated
/// for partial ionisation.
///
/// The file inherits from mpv3, and instead of the Wolfire et al.
/// (2003) and Henney et al. (2009) heating/cooling rates, it sets
/// the temperature based on Hydrogen ion fraction.
///
/// Even if H is molecular, we assume that there is a PDR ahead of
/// the ionisation front, so the ionising photons always hit atomic
/// hydrogen.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15
/// - 2013.03.10 JM: Changed ions/electrons so He is always neutral.
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2013.06.27 JM: changed T(xp) function in get_temperature.
/// - 2013.07.20 JM: Tested a bunch of things; moved PI rate from a
///    function call to a (faster) in-place evaluation.
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2016.03.16 JM: Changed Min/Max neutral fraction, and the method
///    of calculating temperature, to try to improve code.
/// - 2018.03.20 JM: Renamed file.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "microphysics/MPv7.h"

using namespace std;

// ##################################################################
// ##################################################################

MPv7::MPv7(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr,   ///< Number of tracer variables in state vector.
    const std::string* tracers,   ///< List of what the tracer variables mean.
    struct which_physics* ephys,  ///< extra physics stuff.
    struct rad_sources* rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    MPv3(nd, csys, nv, ntr, tracers, ephys, rsrcs, g)
{
#ifdef TESTING
  cout << "MPv7 constructor setting up.\n";
#endif
  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  //
  double X        = 1.0 - EP->Helium_MassFrac;
  mean_mass_per_H = pconst.m_H() / X;  // this is mass per H nucleon.
  //
  // Number of ionised particles per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NION = 1.0;  // +0.25*EP->Helium_MassFrac/X;

  //
  // Number of electrons per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NELEC = 1.0;  // +0.25*EP->Helium_MassFrac/X;

  //
  // The metal MassFrac is unused, so we use that to decide if the
  // neutral medium is neutral or molecular
  // If Metal_MassFrac >0.5 then neutral H is molecular (larger mu).
  // Otherwise, neutral H is atomic (smaller mu).
  //
  if (EP->Metal_MassFrac > 0.5) {
    TTI_Mol = 0.5;
  }
  else {
    TTI_Mol = 1.0;
  }
  TTI_Nnt = TTI_Mol + 0.25 * EP->Helium_MassFrac / X;

  TTI_Thi = EP->MaxTemperature;
  TTI_Tlo = EP->MinTemperature;

  //
  // Set these to be as close to 0 and 1 as possible, so that the
  // temperature is set correctly!
  //
  Min_NeutralFrac = 1.0e-15;
  Max_NeutralFrac = 1.0 - 1.0e-15;
#ifdef TESTING
  cout << "MPv7: Y=" << EP->Helium_MassFrac;
  cout << ", Z=" << EP->Metal_MassFrac << ", mmpH=" << mean_mass_per_H;
  cout << ", NION=" << JM_NION << ", NELEC=" << JM_NELEC << "\n";
#endif  // TESTING
  return;
}

// ##################################################################
// ##################################################################

MPv7::~MPv7()
{
#ifdef TESTING
  cout << "MPv7 destructor.\n";
#endif
  return;
}

// ##################################################################
// ##################################################################

int MPv7::convert_prim2local(
    const double* p_in,  ///< primitive vector from grid cell (length nv_prim)
    double* p_local)
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  mpv_nH = p_in[RO] / mean_mass_per_H;

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0] = 1.0 - p_in[pv_Hp];
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  //
  // internal energy is not used, so this is not really needed.
  //
  p_local[lv_eint] = p_in[PG] / (gamma_minus_one);

#ifdef TESTING
  //
  // Check for NAN/INF
  //
  for (int v = 0; v < 2; v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics", p_local[v]);
  }
  if (mpv_nH < 0.0 || !isfinite(mpv_nH))
    rep.error("Bad density input to MPv7::convert_prim2local", mpv_nH);
#endif  // TESTING

  return 0;
}

// ##################################################################
// ##################################################################

int MPv7::convert_local2prim(
    const double* p_local,
    const double*
        p_in,      ///< input primitive vector from grid cell (length nv_prim)
    double* p_out  ///< updated primitive vector for grid cell (length nv_prim)
)
{
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];

  //
  // Set output H+ fraction
  //
  p_out[pv_Hp] = 1.0 - p_local[lv_H0];
  p_out[pv_Hp] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_out[pv_Hp]));

  //
  // Output pressure is determined by ion fraction and number dens.:
  // - H+ fraction is y,
  // - Total number density, n_tot = nH*((1-y)*Nnt +y*(Nelec+Nion)),
  // - Temperature is T=Tlo + y*(Thi-Tlo) (from function),
  // so p = n_tot*k*T.
  //
  p_out[PG] = get_ntot(mpv_nH, p_out[pv_Hp]) * k_B
              * get_temperature(mpv_nH, 0.0, p_out[pv_Hp]);

#ifdef TESTING
  // cout <<"nH="<< mpv_nH <<", xp="<< p_out[pv_Hp] <<", ntot=";
  // cout <<get_ntot(mpv_nH,p_out[pv_Hp])<<"\n";
#endif

#ifdef TESTING
  if (p_out[pv_Hp] < 0.0 || p_out[pv_Hp] > 1.0 * (1.0 + JM_RELTOL)
      || !isfinite(p_out[pv_Hp]))
    rep.error(
        "Bad output H+ value in MPv7::convert_local2prim", p_out[pv_Hp] - 1.0);
  if (p_out[PG] < 0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in MPv7::convert_local2prim", p_out[PG]);
#endif  // TESTING

  return 0;
}

// ##################################################################
// ##################################################################

double MPv7::get_temperature(
    const double,    ///< nH (per c.c.): UNUSED
    const double,    ///< E_int (per unit volume): UNUSED
    const double xp  ///< x(H+)
)
{
  //
  // returns gas temperature according to T=(2yT_hi+(1-y)*T_lo)/(1+y)
  // which is RJR William's estimate of the correct temperature to
  // get the mixed cell physics right.
  //
  // double frac = ONE_MINUS_EPS;
  // if (xp>frac) return TTI_Thi;
  // else return (xp/frac*(2.0*TTI_Thi-TTI_Tlo) +TTI_Tlo)/(1.0+xp/frac);

  // Dodgy logarithmic interpolation.
  // if (xp<1.0e-4) return TTI_Tlo;
  // else return TTI_Tlo +(TTI_Thi-TTI_Tlo)*(log10(xp)+4.0)/4.0;

  //
  // This is my estimate for gas temperature to try to account for
  // the mixed-cell physics in unresolved ionisation fronts.
  // It smoothes out the pressure, but in the end provides a worse
  // solution than other forms.
  //
  // return TTI_Tlo*TTI_Thi/(TTI_Thi-xp*(TTI_Thi-TTI_Tlo));

  //
  // returns gas temperature according to T=Tlo + y*(Thi-Tlo),
  //
  return (TTI_Tlo + xp * (TTI_Thi - TTI_Tlo));
}

// ##################################################################
// ##################################################################

double MPv7::get_ntot(
    const double nH,  ///< nH
    const double xp  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
)
{
  //
  // This allows for molecular H neutral gas, with TTI_Mol, which is
  // 0.5 if molecular.  This is the (H0/H2) + (He) + (elect.+ions).
  //
  return ((1.0 - xp) * TTI_Mol + (TTI_Nnt - TTI_Mol)
          + xp * (JM_NELEC + JM_NION))
         * nH;
}

// ##################################################################
// ##################################################################

int MPv7::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double*          ///< extra user-data vector (UNUSED)
)
{

#ifdef TESTING
  // cout <<"MPv7::ydot(): Y="<< EP->Helium_MassFrac;
  // cout <<", Z="<< EP->Metal_MassFrac <<", mmpH="<<mean_mass_per_H;
  // cout <<", NION="<< JM_NION <<", NELEC="<< JM_NELEC;
  // cout <<", Nnt="<<TTI_Nnt<<"\n";
#endif  // TESTING

  //
  // Set initial values
  //
  double OneMinusX = max(NV_Ith_S(y_now, lv_H0), Min_NeutralFrac);
  double x_in      = 1.0 - OneMinusX;
  double ne        = JM_NELEC * x_in * mpv_nH;
  // double T = get_temperature(mpv_nH, 0.0, x_in);

  double temp1 = 0.0;
  // double temp2=0.0;

  //
  // oneminusx_dot is in units of 1/s
  //
  double oneminusx_dot = 0.0;
  //
  // Edot is irrelevant here, so just leave it at zero so it never
  // affects the timesteps.  This makes it a bit inefficient, but
  // it's never going to be used in detailed simulations.
  //
  double Edot = 0.0;

  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0.
  //
  // Hi_coll_ion_rates(T, &temp1, &temp2);
  // oneminusx_dot -= temp1*ne*OneMinusX; // the nH is divided out on both
  // sides.

  //
  // photo-ionisation of H:
  // photoionisation rate uses equation 18 in Mellema et al. 2006 (C2-ray
  // paper), noting that their Gamma is the rate per neutral H, so we multiply
  // by 1-x, as in their equation 11.
  //
  if (N_ion_srcs) {
    //
    // set current cell dTau0
    //

    temp1 = mpv_nH * mpv_delta_S * OneMinusX
            * Hi_monochromatic_photo_ion_xsection(JUST_IONISED);

    switch (ion_src_type) {
      case RT_EFFECT_MFION:
        //
        // Rather than divide the discretised rate by n(H0) and then
        // multiply by (1-x) to get oneminusx_dot, we simply divide by
        // n(H) since this is more numerically stable.  To do this,
        // n(H) is passed to the rate function instead of n(H0).
        //
        // Also, instead of using mpv_dTau0 for the cell optical depth,
        // we use the current optical depth (nH*mpv_delta_S*OneMinusX)
        // so that we allow the photoionisation rate to decrease as the
        // number of neutral atoms decreases during the timestep.  This
        // is more stable.
        //
        oneminusx_dot -= Hi_discrete_multifreq_photoion_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_delta_S, mpv_Vshell);
        break;

      case RT_EFFECT_PION_MONO:
        //
        // hardcoded for hv=(1+epsilon)*13.6 eV  monochromatic source.
        //
        if (temp1 < 0.0001) {
          oneminusx_dot -=
              mpv_NIdot * exp(-mpv_Tau0) * temp1 / mpv_Vshell / mpv_nH;
        }
        else {
          oneminusx_dot -= mpv_NIdot * exp(-mpv_Tau0) * (1.0 - exp(-temp1))
                           / mpv_Vshell / mpv_nH;
        }
        break;

      default:
        rep.error("Bad ion_src_type in dYdt()", ion_src_type);
        break;
    }  // switch
  }

  //
  // radiative recombination of H+
  //
  // oneminusx_dot += Hii_rad_recomb_rate(T) *x_in*ne;
  oneminusx_dot += 2.7e-13 * x_in * ne;

  //
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar
  // neighbourhood.
  //
  // oneminusx_dot -= 1.8e-17*OneMinusX;

  //
  // Set Edot to be given by n*k*Tdot/(g-1) with
  //  Tdot=dot(1-x)*(Tlo-Thi)
  //
  // Edot =
  // get_ntot(mpv_nH,x_in)*k_B*oneminusx_dot*(TTI_Tlo-TTI_Thi)/gamma_minus_one;

  NV_Ith_S(y_dot, lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot, lv_eint) = Edot;
  return 0;
}

// ##################################################################
// ##################################################################

double MPv7::get_recombination_rate(
    const int,           ///< ion index in tracer array (optional).
    const double* p_in,  ///< input state vector (primitive).
    const double g       ///< EOS gamma (optional)
)
{
#ifdef FUNCTION_ID
  cout << "MPv7::get_recombination_rate()\n";
#endif  // FUNCTION_ID
  double rate = 0.0;
  double P[nvl];
  // First convert to local variables.
  convert_prim2local(p_in, P);
  // Now get rate
  rate = 2.7e-13 * mpv_nH * mpv_nH * (1.0 - P[lv_H0]) * (1.0 - P[lv_H0])
         * JM_NELEC;

#ifdef FUNCTION_ID
  cout << "MPv7::get_recombination_rate()\n";
#endif  // FUNCTION_ID
  return rate;
}

// ##################################################################
// ##################################################################
