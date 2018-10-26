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
///    Removed the DO_HEATING ifdef -- now we have 5 different cooling functions.
///
/// - 2011.05.10 JM: Output cooling rates only if myrank==0 for parallel (so processes
///    don't fight over the file and slow down the code (by a lot!)).
/// - 2015.01.13 JM: Added some comments.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2016.06.21 JM: Temperature() threadsafe.
/// - 2018.01.25 JM: added functions to request n(H+),n(H0),n(e-)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/mp_only_cooling.h"

#include <sstream>
using namespace std;

#define SD93_CIE 4
#define SD93_PLUS_HEATING 5
#define WSS09_CIE_PLUS_HEATING 6
#define WSS09_CIE_ONLY_COOLING 7
#define WSS09_CIE_LINE_HEAT_COOL 8

// ##################################################################
// ##################################################################


mp_only_cooling::mp_only_cooling(
      const int nv,
      struct which_physics *ephys, ///< pointer to extra physics flags.
      struct rad_sources *rsrcs    ///< radiation sources.
      )
: microphysics_base(ephys,rsrcs),
  cooling_function_SD93CIE(),
  Hummer94_Hrecomb(),
  nv_prim(nv)
{

  //
  // First check that we are updating energy
  //
  if (!EP->update_erg) {
    rep.error("requested cooling microphysics but no energy update",
	      DONT_CALL_ME);
  }

  //
  // Next check that we are limiting timestep by the cooling time,
  // and if not, then set it!
  //
  if (EP->MP_timestep_limit != 1) {
    cout <<"\t\tmp_only_cooling: timestep limiting not set correctl";
    cout <<"y.  Changing from "<<EP->MP_timestep_limit<<" to ";
    EP->MP_timestep_limit = 1;
    cout << EP->MP_timestep_limit <<"\n";
  }

  //
  // Mean masses per species/atom: we assume cosmic abundances which gives
  // rho \simeq 1.40 m_p n_H.  For fully ionised gas we have 
  //  n_i \simeq 1.1 n_H (from helium).
  //  n_e \simeq 1.2 n_H (from doubly ionised helium)
  //  n_tot \simeq 2.3 n_H (summing the previous lines).
  // So we define muH=1.4mp, mue=1.167mp, mui=1.273mp, mutot=0.609mp
  //
  Mu =1.40*pconst.m_p();
  Mu_tot = 0.609*pconst.m_p();
  Mu_tot_over_kB = Mu_tot/pconst.kB();
  Mu_elec = 1.167*pconst.m_p();
  Mu_ion  = 1.273*pconst.m_p();
  inv_Mu2 = 1.0/(Mu*Mu);
  inv_Mu2_elec_H = 1.0/(Mu_elec*Mu);


  //
  // Next select which cooling function to set up.
  //
  cooling_flag = EP->cooling;

  switch (cooling_flag) {
  case SD93_CIE:
  case SD93_PLUS_HEATING:
    setup_SD93_cie();
    break;
  case WSS09_CIE_PLUS_HEATING:
  case WSS09_CIE_ONLY_COOLING:
    setup_WSS09_CIE();
    break;
  case WSS09_CIE_LINE_HEAT_COOL:
    cout <<"\tRequested fully ionized gas with WSS09 cooling at high";
    cout <<" temperatures,\n\tand photoionized gas at nebular";
    cout <<" temperatures, with T_eq approx 8000 K.\n";
    setup_WSS09_CIE_OnlyMetals();
    break;
  default:
    rep.error("Invalid cooling flag in mp_only_cooling",
	      cooling_flag);
    break;
  }

#ifdef TESTING
  ostringstream opfile; opfile << "coolingNOCHEM_" << cooling_flag << ".txt";
  ofstream outf(opfile.str().c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"Cooling Curve Data: Temperature(K) Rate(erg/cm^3/s) (n=1 per cc)\n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=1.0e0;
  do {
    outf << t <<"\t"<< Edot(2.34e-24,t) <<"\n";
    t *=1.05;
  } while (t<1.0e10);
  outf.close();
#endif 

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  MaxT_allowed = EP->MaxTemperature;
  MinT_allowed = EP->MinTemperature;
  if (MinT_allowed <1.0   || MinT_allowed>1.0e6 ) MinT_allowed=1.0;
  if (MaxT_allowed <1.0e2 || MaxT_allowed>3.0e10) MaxT_allowed=1.0e8; 
  cout <<"\t\tAllowed Temperature range: T_min="<<MinT_allowed<<"  T_max=";
  cout <<MaxT_allowed<<" (MAKE SURE THIS IS SET APPROPRIATELY!\n";
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return;
}


// ##################################################################
// ##################################################################




mp_only_cooling::~mp_only_cooling()
{
}



// ##################################################################
// ##################################################################



//
// update internal energy over full timestep.
//
int mp_only_cooling::TimeUpdateMP(
      const pion_flt *p_in, ///< Primitive Vector to be updated
      pion_flt *p_out, ///< Destination Vector for updated values.
      const double dt, ///< Time Step to advance by.
      const double gamma, ///< EOS gamma.
      const int, ///< Switch for what type of integration to use (not used here!)
      double *Tf ///< final temperature.
      )
{

//#ifdef TEST_INF
  for (int v=0;v<nv_prim;v++) {
    if (!isfinite(p_in[v])) {
      cout <<"NAN/INF input to  mp_only_cooling::TimeUpdateMP: ";
      cout <<" v="<<v<<", val="<<p_in[v]<<"\n";
      return 1;
    }
  }
//#endif

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // First get gas temperature, and if it outside the allowed ranges,
  // then set it to the relevant limiting min/max value.
  //
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];
  double T = p_out[PG]*Mu_tot_over_kB/p_out[RO];
  if (T<MinT_allowed) {
    Set_Temp(p_out,MinT_allowed,gamma);
    T = MinT_allowed;
  }
  if (T>MaxT_allowed) {
    Set_Temp(p_out,MaxT_allowed,gamma);
    T = MaxT_allowed;
  }
  double Eint0 = p_in[PG]/(gamma-1.0);
#else
  double T = p_in[PG]*Mu_tot_over_kB/p_in[RO];
  double Eint0 = p_in[PG]/(gamma-1.0);
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  //
  // Now we do a few RK4 steps
  //
  size_t nstep=1;
  double step_dt;
  double Etemp, E0_step, k1, k2, k3, k4;
  do {
    step_dt = dt/nstep;
    E0_step = Eint0;
    T = E0_step*(gamma-1.0)*Mu_tot_over_kB/p_in[RO];

    if (nstep>16) {
      cout <<"Warning: cooling nsteps="<<nstep<<"\n";
      rep.printVec("p_in",p_in,nv_prim);
      cout <<"E0_step="<<E0_step<<", Ein="<<Eint0;
      cout <<", T="<<T<<", Edot="<<Edot(p_in[RO],T)<<"\n";
    }

    for (size_t v=0;v<nstep;v++) {

      k1 = step_dt*Edot(p_in[RO],T);
      Etemp = E0_step+0.5*k1;

      T = Etemp*(gamma-1.0)*Mu_tot_over_kB/p_in[RO];
      k2 = step_dt*Edot(p_in[RO],T);
      Etemp = E0_step+0.5*k2;
      
      T = Etemp*(gamma-1.0)*Mu_tot_over_kB/p_in[RO];
      k3 = step_dt*Edot(p_in[RO],T);
      Etemp = E0_step+k3;

      T = Etemp*(gamma-1.0)*Mu_tot_over_kB/p_in[RO];
      k4 = step_dt*Edot(p_in[RO],T);

      E0_step = E0_step +(k1 +2.0*k2 +2.0*k3 +k4)/6.0;
      T = E0_step*(gamma-1.0)*Mu_tot_over_kB/p_in[RO];

      if (nstep>16) {
        cout <<"MP_only_cooling integration: step="<<v<<", T=";
        cout <<T<<", E0_step="<<E0_step<<"\n";
      }
    }
    
    if (nstep>16) {
      //cout <<"Warning: cooling nsteps="<<nstep<<"\n";
      //rep.printVec("p_in",p_in,nv_prim);
      cout <<"END: T="<<T<<", E0_step="<<E0_step<<", Ein="<<Eint0<<"\n";
    }
    nstep *=2;
  } while (T<MinT_allowed || T>MaxT_allowed || !isfinite(E0_step));

  p_out[PG] = E0_step*(gamma-1);
  *Tf = p_out[PG]*Mu_tot_over_kB/p_out[RO];

  return 0;
}


// ##################################################################
// ##################################################################




//
// Reset pressure so it corresponds to requested temperature.
//
int mp_only_cooling::Set_Temp(
      pion_flt *p_in, ///< primitive vector.
      const double T, ///< temperature requested.
      const double gam ///< eos gamma.
      )
{
#ifdef TESTING
  //
  // check for bad input
  //
  if (T<=0.0) {
    cout <<"mp_only_cooling::Set_Temp() Requested negative ";
    cout <<"temperature! T="<<T<<"\n";
    return 1;
  }
#endif // TESTING
  // cout <<"\n\t\t***** Set_Temp(): p="<<p_in[PG]<<", rho="<< p_in[RO];
  // cout <<", T="<<T<<", new p="<<p_in[RO]*T/Mu_tot_over_kB<<"\n";
  p_in[PG] = p_in[RO]*T/Mu_tot_over_kB;

  return 0;
}



// ##################################################################
// ##################################################################



//
// Returns the gas temperature.  Assumes primitive vector is in
// cgs units and ionised gas with mu=0.7m_p.
//
double mp_only_cooling::Temperature(
    const pion_flt *p_in, ///< primitive vector
    const double   ///< eos gamma
    )
{
  return p_in[PG]*Mu_tot_over_kB/p_in[RO];
}


// ##################################################################
// ##################################################################

//
// Get electron number density (cm^{-3})
//
double mp_only_cooling::get_n_elec(
      const pion_flt *p_in ///< primitive state vector.
      )
{
  // H and He are fully ionized
  return p_in[RO]/Mu_elec;
}


// ##################################################################
// ##################################################################


//
// Get H+ number density (cm^{-3})
//
double mp_only_cooling::get_n_Hplus(
      const pion_flt *p_in ///< primitive state vector.
      )
{
  // fully ionized hydrogen
  return p_in[RO]/Mu;
}


// ##################################################################
// ##################################################################



//
// Get neutral H number density (cm^{-3})
//
double mp_only_cooling::get_n_Hneutral(
      const pion_flt *p_in ///< primitive state vector.
      )
{
  // Assume neutral fraction of 1e-12 (arbitrary)
  return 1.0e-12*p_in[RO]/Mu;
}


// ##################################################################
// ##################################################################



//
// This returns the minimum timescale of the times flagged in the
// arguments.  Time is returned in seconds.  Only cooling flag has
// any effect here.
//
double mp_only_cooling::timescales(
    const pion_flt *p_in, ///< Current cell.
    const double gam,   ///< EOS gamma.
    const bool, ///< set to true if including cooling time.
    const bool, ///< set to true if including recombination time.
    const bool  ///< set to true if including photo-ionsation time.
    )
{
  //
  // First get temperature and internal energy density.
  //
  double Eint = p_in[PG]/(gam-1.0);
  double T    = p_in[PG]*Mu_tot_over_kB/p_in[RO];
  double mintime=1.0e99;
  double rate=0.0;

  //
  // we skip the cooling time if the temperature is already very low;
  // a reasonable temperature could be 10K; there is unlikely to be
  // significant cooling below this temperature, and even if there is
  // it won't be by a large fraction.
  //
  // We get the cooling rate per unit volume, first for the actual
  // temperature, and then for T/2, and we take the max. of this.
  // This is to try to ensure we will get an accurate one-step 
  // integration.
  //
  if (T>=10.0) {
    rate = max(fabs(Edot(p_in[RO],T)), fabs(Edot(p_in[RO],0.5*T)));
    //
    // Cooling time is then the energy per unit volume divided by the rate.
    // Bear in mind we could have net heating, so take fabs(rate)
    //
    mintime = min(mintime, Eint/fabs(rate));
  }
  
  return mintime;
}


// ##################################################################
// ##################################################################




double mp_only_cooling::Edot(
            const double rho, ///< mass density (g/cm3)
            const double T    ///< Temperature (K)
            )
{
  double rate=0.0;
  switch (cooling_flag) {
    case SD93_CIE:
    rate = Edot_SD93CIE_cool(rho,T);
    break;

    case SD93_PLUS_HEATING:
    rate = Edot_SD93CIE_heat_cool(rho,T);
    break;

    case WSS09_CIE_ONLY_COOLING:
    rate = Edot_WSS09CIE_cool(rho,T);
    break;

    case WSS09_CIE_PLUS_HEATING:
    rate = Edot_WSS09CIE_heat_cool(rho,T);
    break;

    case WSS09_CIE_LINE_HEAT_COOL:
    rate = Edot_WSS09CIE_heat_cool_metallines(rho,T);
    break;

    default:
    rep.error("bad cooling flag in mp_only_cooling::Edot",cooling_flag);
    break;
  }
  return rate;
}


// ##################################################################
// ##################################################################




/// only cooling, uses SD93-CIE curve only.
double mp_only_cooling::Edot_SD93CIE_cool(
            const double rho, ///< mass density (g/cm3)
            const double T    ///< Temperature (K)
            )
{
  return -(rho*rho/Mu_elec/Mu_ion)*cooling_rate_SD93CIE(T);
}


// ##################################################################
// ##################################################################


//
// cooling from SD93-CIE plus heating assuming full ionisation of H, where 
// the heating rate equals recombination rate times 5eV/ionisation.
//
double mp_only_cooling::Edot_SD93CIE_heat_cool(
            const double rho, ///< mass density (g/cm3)
            const double T    ///< Temperature (K)
            )
{
  return (rho*rho)*(2.733e-21*exp(-0.782991*log(T))/Mu_elec/Mu 
                    -cooling_rate_SD93CIE(T)/Mu_elec/Mu_ion);
}


// ##################################################################
// ##################################################################




// only cooling, uses WSS09-CIE cooling function.
double mp_only_cooling::Edot_WSS09CIE_cool(
              const double rho, ///< mass density (g/cm3)
              const double T    ///< Temperature (K)
              )
{
  return -(rho*rho/Mu/Mu)*cooling_rate_SD93CIE(T);
}


// ##################################################################
// ##################################################################



// 
// cooling from WSS09-CIE plus heating assuming full ionisation of H, where 
// the heating rate equals recombination rate times 5eV/ionisation.
//
double mp_only_cooling::Edot_WSS09CIE_heat_cool(
              const double rho, ///< mass density (g/cm3)
              const double T    ///< Temperature (K)
              )
{
  return (rho*rho)*(2.733e-21*exp(-0.782991*log(T))/Mu_elec/Mu -cooling_rate_SD93CIE(T)/Mu/Mu);
}


// ##################################################################
// ##################################################################




double mp_only_cooling::Edot_WSS09CIE_heat_cool_metallines(
            const double rho, ///< mass density (g/cm3)
            const double T    ///< Temperature (K)
            )
{
  //
  // COOLING:
  // First forbidden line cooling of e.g. OII,OIII, dominant in HII regions.
  // This is collisionally excited lines of photoionised metals.
  // (Henney et al. 2009,MN,398,157,eq.A9)
  // I exponentially damped this rate at high temperatures because physically
  // oxygen will move to higher ionisation states and that is not accounted for
  // in this simple formula.
  //
  // TODO: Check Oxygen abundance is ok. It is set to 5.81e-4 from Lodders 
  //       (2003,ApJ,591,1220,Tab.4) which is the 'proto-solar nebula' value.
  //       The photospheric value is lower 4.9e-4, and that is used by Wiersma
  //       et al. (2009,MN,393,99).
  //
  // The CIE cooling is from WSS09, which is normalised by [n(H)]^2, not n(e-)n(H+),
  // hence the differing multipliers on the two functions.
  //
  double rho2 = rho*rho;
  double rate = -1.69e-22 *exp(-33610.0/T -(2180.0*2180.0/T/T)) *rho2*inv_Mu2_elec_H *exp(-T*T/5.0e10);
  //cout <<"M="<<rate;
  rate = min(rate, -cooling_rate_SD93CIE(T)*rho2*inv_Mu2);
  //cout <<", CIE="<<-cooling_rate_SD93CIE(T)*rho2*inv_Mu2;
  //
  // Now Hydrogen cooling due to recombinations and Bremsstrahlung.
  //
  rate -= Hii_total_cooling(T) *rho2*inv_Mu2_elec_H;
  //cout <<", H="<<-Hii_total_cooling(T) *rho2*inv_Mu2_elec_H;
  //
  // Now need to add Helium Bremsstrahlung (using eq.5.15b from Rybicki &
  // Lightman (1978), p.162, with Z=2, n(He)/n(H)=0.1, g_B=1.2 (Gaunt-factor))
  //
  rate -= 6.72e-28*sqrt(T)*rho2*inv_Mu2_elec_H;
  //
  // heating rate = 5eV*RRR*n(e-)*n(H+) = 8.01e-12*Hii_rad_recomb_rate(T)*rho*rho/Mu_elec/Mu
  // assuming all hydrogen is ionised.
  //
  rate += 8.01e-12*Hii_rad_recomb_rate(T)*rho2*inv_Mu2_elec_H;
  //cout <<", Heat="<<8.01e-12*Hii_rad_recomb_rate(T)*rho2*inv_Mu2_elec_H<<"\n";
  return rate;
}

// ##################################################################
// ##################################################################



// int main()
// {
//  //
//  // test the only_cooling class:
//  // compile with 
//  // g++ ../mp_only_cooling.cc ../cooling_SD93_cie.cc ../../global.cc ../../cell_interface.cc -o a.out -lreadline
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
//  EP.MinTemperature = 10.0; ///< Minimum temperature to allow in the simulation.
//  EP.MaxTemperature = 2.0e10; ///< Maximum temperature to allow in the simulation.
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
 //  //  cout <<"dens="<<p[RO]<<", cooling time="<<MP.timescales(p,g,1,1,1)<<"\n";
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
 //      rep.error("integration error",err);
 //    }
 //    t+=h;
 //  } while (MP.Temperature(p,g)>1.0e4);

 //  return 0;
 // }

