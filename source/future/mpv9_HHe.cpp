///
/// \file mpv9_HHe.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running calculations with a simple 
/// two-temperature isothermal equation of state, where T is T_low
/// when gas is neutral, and T_high when ionised, and linearly
/// interpolated for partial ionisation.
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
/// - 2013.08.23 JM: Debugging.
/// - 2013.09.16 JM: Debugging (not finished yet!).



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifndef EXCLUDE_MPV9

#include "mpv9_HHe.h"
#include "global.h"
using namespace std;



// ##################################################################
// ##################################################################


void mpv9_HHe::get_error_tolerances(
                double *reltol, ///< relative error tolerance.
                double *atol ///< absolute error tolerances
                )
{
  *reltol = 1.0e-4;
  atol[lv_H0]   = 1.0e-15; ///< minimum neutral fraction I care about.
  atol[lv_He0]  = 1.0e-16; ///< minimum neutral fraction I care about.
  atol[lv_He1]  = 1.0e-16; ///< minimum neutral fraction I care about.
  atol[lv_E] = 1.0e-18; ///< for n=1.0, T=1.0e2, ==> E=1.5e-14
  return;
}



// ##################################################################
// ##################################################################



void mpv9_HHe::get_problem_size(
                  int *ne, ///< number of equations
                  int *np  ///< number of parameters in user_data vector.
                  )
{
  *ne = N_equations;
  *np = N_extradata;
  return;
}



// ##################################################################
// ##################################################################



mpv9_HHe::mpv9_HHe(
          const int nv,                ///< Total number of variables in state vector
          const int ntracer,           ///< Number of tracer variables in state vector.
          const std::string &trtype,   ///< List of what the tracer variables mean.
          struct which_physics *ephys, ///< extra physics stuff.
          const double g               ///< EOS gamma
	  )
  :
  gamma(g), gamma_minus_one(g-1.0), nv_prim(nv)
{
#ifdef TESTING
  cout <<"mpv9_HHe constructor setting up.\n";
#endif

  EP = ephys;

  // ----------------------------------------------------------------
  //
  // Local vector: E_int, y(H0), y(He0), y(He+)
  //
  nvl     = 4;    // Number of local variables.
  lv_E = 0;    // E_{int}.
  lv_H0   = 1;    // x(H0).
  lv_He0  = 2;    // x(He0).
  lv_He1  = 3;    // x(He1).
  nH = 0.0;
  tau =0; tau  = mem.myalloc(tau,4);
  dtau=0; dtau = mem.myalloc(dtau,4);

  N_extradata = 1;
  N_equations = nvl;
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  //
  mean_mass_per_H = m_p()/EP->H_MassFrac; // this is mass per H nucleon.

  //
  // number fraction of H, He = n(He)/n(H)
  //
  X_H  = 1.0;
  X_HE = 0.25*EP->Helium_MassFrac/EP->H_MassFrac;

  //
  // The metal MassFrac is not used for the number of particles, but
  // it is used for heating/cooling via the "metallicity" parameter.
  //
  metallicity = EP->Metal_MassFrac/0.0142;

  //
  // max and min values we allow:
  //
  Min_Nfrac     = 0.0;
  Max_Nfrac     = 1.0;
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // ------- Set up tracer variables (H0, He0, He1). --------
  // ----------------------------------------------------------------
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  //
  // first 6 chars are the type, then list of tracers, each 6 chars long.
  //
  int len = (trtype.length() +5)/6 -1;
  if (len!=ntracer) {
    cout <<"warning: string doesn't match ntracer.  ";
    cout <<"make sure this looks ok: "<<trtype<<"\n";
  }
  //
  // Find H+ fraction in tracer variable list.
  //
  int ftr = nv_prim -ntracer; // first tracer variable.
  string s;
  pv_H0  = 10000;
  pv_He0 = 10000;
  pv_He1 = 10000;
  for (int i=0;i<len;i++) {
    s = trtype.substr(6*(i+1),6); // Get 'i'th tracer variable.
    if      (s=="H0____" || s=="HI____") {
      pv_H0 = ftr+i;
      cout <<"\t\tGot H0 as the "<<pv_H0<<"th el. of P[].\n";
    }
    else if (s=="He0___" || s=="HeI___") {
      pv_He0 = ftr+i;
      cout <<"\t\tGot He0 as the "<<pv_He0<<"th el. of P[].\n";
    }
    else if (s=="He1___" || s=="HeII__") {
      pv_He1 = ftr+i;
      cout <<"\t\tGot He+ as the "<<pv_He1<<"th el. of P[].\n";
    }
  }
  if (pv_H0>nv_prim)
    rep.error("No H0 fraction found in tracer list",trtype);
  if (pv_He0>nv_prim)
    rep.error("No He0 fraction found in tracer list",trtype);
  if (pv_He1>nv_prim)
    rep.error("No He+ fraction found in tracer list",trtype);
  // ================================================================
  // ================================================================

  //
  // ----------------------------------------------------------------
  // ------------------------ RADIATION SOURCE ----------------------
  //
  bool have_src=false;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    //
    // Need to set up the multifrequency tables.  If a HHe-multifreq
    // source is not found then bug out.
    //
    for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
      if (  (SimPM.RS.sources[isrc].type  ==RT_SRC_SINGLE) &&
            (SimPM.RS.sources[isrc].effect==RT_EFFECT_HHE_MFQ) ) {
        //
        // All parameters here are unused except Teff,Rstar, both of
        // which are in cgs units.
        //
        Setup_photoionisation_rate_table(
                SimPM.RS.sources[isrc].Tstar,
                SimPM.RS.sources[isrc].Rstar*Rsun(), 1,1,1,1,1,1);
        have_src=true;
      }
    }
  }
  if (!have_src) rep.error("No radiation source!",SimPM.RS.Nsources);
  // ================================================================
  // ================================================================
  //

  // ----------------------------------------------------------------
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  setup_cvode_solver_without_Jacobian();
  // ================================================================
  // ================================================================

  return;
}



// ##################################################################
// ##################################################################



mpv9_HHe::~mpv9_HHe()
{
#ifdef TESTING
  cout <<"mpv9_HHe destructor.\n";
#endif
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
  tau  = mem.myfree(tau);
  dtau = mem.myfree(dtau);
  return;
}



// ##################################################################
// ##################################################################



int mpv9_HHe::TimeUpdateMP_RTnew(
          const double *p_in,
          const int N_heat,
          const std::vector<struct rt_source_data> &heat_src,
          const int N_ion,
          const std::vector<struct rt_source_data> &ion_src,
          double *p_out,
          const double dt,
          const double,    ///< EOS gamma (UNUSED)
          const int, ///< Switch for integration scheme (UNUSED)
          double *random_stuff ///< final temperature (UNUSED)
          )
{
#ifdef FUNCTION_ID
  cout <<"mpv9_HHe::TimeUpdateMP_RTnew()\n";
#endif // FUNCTION_ID

  if (N_ion!=1) rep.error("No ionising sources!",N_ion);

  //
  // First set local variables for state vector and radiation source
  // properties.
  //
  int err=0;
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to mpv9_HHe::TimeUpdateMP_RTnew()",err);
  }
  for (size_t v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];
#ifdef MPV9_DEBUG
  rep.printVec("update: P",P,nvl);
  rep.printVec("update: Y",NV_DATA_S(y_in),nvl);
#endif // MPV9_DEBUG

  interpret_radiation_data(N_heat,heat_src,N_ion,ion_src);

  //
  // Calculate y-dot[] to see if anything is changing significantly over dt
  //
  double maxdelta=0.0;
  err = ydot(0, y_in, y_out, 0);
  if (err) 
    rep.error("dYdt() returned an error in mpv9_HHe::TimeUpdateMP_RTnew()",err);
  for (size_t v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out,v)*dt/NV_Ith_S(y_in,v)));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < 0.01) {
    for (size_t v=0;v<nvl;v++) {
      NV_Ith_S(y_out,v) = NV_Ith_S(y_in,v) + dt*NV_Ith_S(y_out,v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      rep.error("integration failed: mpv9_HHe::TimeUpdateMP_RTnew()",err);
    }
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (size_t v=0;v<nvl;v++) P[v] = NV_Ith_S(y_out,v);
  //
  // That's it! Convert back to primitive variable in p_out[]. 
  //
  err = convert_local2prim(P,p_in,p_out);

#ifdef FUNCTION_ID
  cout <<"mpv9_HHe::TimeUpdateMP_RTnew()\n";
#endif // FUNCTION_ID
  return err;
}



// ##################################################################
// ##################################################################



double mpv9_HHe::get_recombination_rate(
          const int id,       ///< ion index in tracer array (optional).
          const double *p_in, ///< input state vector (primitive).
          const double g      ///< EOS gamma (optional)
          )
{
#ifdef FUNCTION_ID
  cout <<"mpv9_HHe::get_recombination_rate()\n";
#endif // FUNCTION_ID
  double rate=0.0;
  double P[nvl];
  size_t i;
  //
  // First convert to local variables.
  //
  convert_prim2local(p_in,P);
  //
  // Now get rate
  //
  switch (id) {
    case ION_H_N:
    i = lv_H0;
    rate = 2.7e-13*P[lv_nH]*P[lv_nH]*(1.0-P[i])*(1.0-P[i]);
    break;

    case ION_HE_N:
    i = lv_He0;
    rate = 2.7e-13*P[lv_nH]*P[lv_nH]*(1.0-P[i])*(1.0-P[i]);
    break;
    
    case ION_HE_P:
    i = lv_He1;
    rate = 2.7e-13*P[lv_nH]*P[lv_nH]*(1.0-P[i])*(1.0-P[i]);
    break;
  }
  //cout <<"rate="<<rate<<"\n";

#ifdef FUNCTION_ID
  cout <<"mpv9_HHe::get_recombination_rate()\n";
#endif // FUNCTION_ID
  return rate;
}



// ##################################################################
// ##################################################################



double mpv9_HHe::get_n_el(
        const double *pv, ///< primitive state vector.
        const int id      ///< integer identifier for the element.
        )
{
  double P[nvl];
  double ans=0.0;
  convert_prim2local(pv,P);
  switch (id) {
    case EL_H:
    ans = nH;
    break;

    case EL_HE:
    ans = nH*X_HE;
    break;
    
    default:
    cerr <<" mpv9_HHe::get_n_el() unknown element "<<id<<"\n";
    ans = -1.0e99;
    break;
  }
  return ans;
}




// ##################################################################
// ##################################################################



double mpv9_HHe::Temperature(
            const double *pv, ///< primitive vector
            const double      ///< eos gamma
            )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    cout <<"mpv9_HHe::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];
  convert_prim2local(pv,P);
  return get_temperature(P);
}



// ##################################################################
// ##################################################################




int mpv9_HHe::Set_Temp(
          double *p_pv,   ///< primitive vector.
          const double T, ///< temperature
          const double g  ///< eos gamma.
          )
{
  //
  // Check for T out of range
  //
  if      (T>EP->MaxTemperature) {
    return (Set_Temp(p_pv,EP->MaxTemperature,g));  
  }
  else if (T<EP->MinTemperature) {
    return (Set_Temp(p_pv,EP->MinTemperature,g));  
  }

  //
  // Get local vector, change internal energy to the appropriate
  // value for the requested temperature, and then convert back to
  // primitive variables and return.
  //
  double P[nvl];
  int err = convert_prim2local(p_pv,P);
  P[lv_E] = get_ntot(P)*kB()*T/(gamma_minus_one);
  err += convert_local2prim(P, p_pv, p_pv);
  return err;
}


// ##################################################################
// ##################################################################


///
/// This returns the minimum timescale of all microphysical processes, including
/// reaction times for each species and the total heating/cooling time for the gas.
/// It requires the radiation field as an input, so it has substantially greater
/// capability than the other timescales function.
/// Default setting is DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double mpv9_HHe::timescales_RT(
                    const double *p_in, ///< Current cell state vector.
                    const int N_heat,      ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &heat_src,
                    ///< list of UV-heating column densities and source properties.
                    const int N_ion,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &ion_src,
                    ///< list of ionising src column densities and source properties.
                    const double   ///< EOS gamma.
                    )
{
  int err=0;
  //
  // First convert to local variables.
  //
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to mpv9_HHe::timescales_RT()",err);
  }
  for (size_t v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];
#ifdef MPV9_DEBUG
  rep.printVec("update: P",P,nvl);
  rep.printVec("update: Y",NV_DATA_S(y_in),nvl);
#endif // MPV9_DEBUG

  interpret_radiation_data(N_heat,heat_src,N_ion,ion_src);

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in mpv9_HHe::timescales_RT()",err);
  }

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t=HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the energy change 
  // (heating timescale) is not required for accurately tracking ionisation fronts 
  // (although it may be needed for cooling!).
  // For testing purposes there are ifdefs to allow the code to use a relative 
  // change in neutral fraction and/or the relative change in energy as the 
  // timestep criterion, rather than the default of absolute change in neutral
  // fraction.
  //
  //#ifdef USE_RELATIVE_NEUFRAC_DTLIMIT
  //t = min(t,DTFRAC*max(5.0e-2,NV_Ith_S(y_in, lv_H0))/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"using neutral fraction.\n";
  //#else
  t = min(t,0.25/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  t = min(t,0.5*X_HE/(fabs(NV_Ith_S(y_out, lv_He0))+TINYVALUE));
  t = min(t,0.5*X_HE/(fabs(NV_Ith_S(y_out, lv_He1))+TINYVALUE));
  //t = min(t,(NV_Ith_S(y_in, lv_H0)+TINYVALUE)/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"limit by dx: dt="<<0.25/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\n";
  //#endif


#ifdef MPV9_DEBUG
  if (t<3.16e9) {
  cout <<"MP timescales: xdot="<<NV_Ith_S(y_out, lv_H0);
  cout <<", Edot="<<NV_Ith_S(y_out, lv_E)<<" t_x="<<t;
  }
#endif // MPV9_DEBUG

  //#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  ////
  //// Now cooling/heating time to limit to X% change in energy).
  ////
  //t = min(t,DTFRAC*P[lv_E]/(fabs(NV_Ith_S(y_out, lv_E))+TINYVALUE));
  ////cout <<"using fractional energy change.\n";
  ////cout <<"limit by dE: dt="<<DTFRAC*P[lv_E]/(fabs(NV_Ith_S(y_out, lv_E))+TINYVALUE)<<"\n";
  //#endif

#ifdef MPV9_DEBUG
  if (t<3.16e9) {
  cout <<" and min(t_x,t_e)="<<t<<",  "; rep.printVec("P[1-x,E]",P,nvl);
  }
#endif // MPV9_DEBUG
  return t;
}



// ##################################################################
// ##################################################################



int mpv9_HHe::Tr(const string t)
{
  if      (t=="H0____" || t=="HI____") {
    return lv_H0;
  }
  else if (t=="He0___" || t=="HeI___") {
    return lv_He0;
  }
  else if (t=="He1____" || t=="HeII__") {
    return lv_He1;
  }
  else {
    rep.error("bad input string",t);
  }
  return 100;
}



// ##################################################################
// ##################################################################



// ------------------------------------------------------------------
// ----------------- protected functions follow ---------------------
// ------------------------------------------------------------------



// ##################################################################
// ##################################################################



int mpv9_HHe::convert_prim2local(
          const double *p_in, ///< primitive vector from grid cell (length nv_prim)
          double *p_local
          )
{
  //
  // Set internal energy density and number density of H.
  //
  p_local[lv_E] = p_in[PG]/(gamma_minus_one);
  nH = p_in[RO]/mean_mass_per_H;

  //
  // Set y(i) to be within the required range: [0,1].
  //
  p_local[lv_H0]   = p_in[pv_H0];
  p_local[lv_He0]  = p_in[pv_He0];
  p_local[lv_He1]  = p_in[pv_He1];
#ifdef MPV9_DEBUG
  if (p_local[lv_H0]<-0.01 ||  p_local[lv_H0]>1.01) {
    cout <<"H0 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_H0] <<"\n";
  }
  if (p_local[lv_He0]<-0.01 ||  p_local[lv_He0]>1.01) {
    cout <<"He0 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_He0] <<"\n";
  }
  if (p_local[lv_He1]<-0.01 ||  p_local[lv_He1]>1.01) {
    cout <<"He1 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_He1] <<"\n";
  }
#endif
  p_local[lv_H0]  = max(Min_Nfrac, min(Max_Nfrac, p_local[lv_H0]));
  p_local[lv_He0] = max(Min_Nfrac, min(Max_Nfrac, p_local[lv_He0]));
  p_local[lv_He1] = max(Min_Nfrac, min(Max_Nfrac, p_local[lv_He1]));

  if (p_local[lv_E]<=0.0) {
#ifdef MPV9_DEBUG
    cout <<"Negative pressure input to mpv9_HHe: ";
    cout <<p_local[lv_E]<<"\n";
#endif
    p_local[lv_E] = get_ntot(p_local)*kB()*EP->MinTemperature/(gamma_minus_one);
  }

#ifdef TESTING
  //
  // Check for NAN/INF
  //
  for (size_t v=0;v<nvl;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to mpv9_HHe",v);
  }
  if (nH<0.0 || !isfinite(nH))
    rep.error("Bad density input mpv9_HHe::convert_prim2local",nH);
#endif // TESTING
  
  return 0;
}


// ##################################################################
// ##################################################################


int mpv9_HHe::convert_local2prim(
            const double *p_local,
            const double *p_in, ///< input primitive vector [nv_prim]
            double *p_out       ///< updated primitive vector
            )
{
  for (size_t v=0;v<nv_prim;v++) p_out[v] = p_in[v];

  //
  // Set output ion fractions
  //
#ifdef MPV9_DEBUG
  if (p_local[lv_H0]<-0.01 ||  p_local[lv_H0]>1.01) {
    cout <<"H0 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_H0] <<"\n";
  }
  if (p_local[lv_He0]<-0.01 ||  p_local[lv_He0]>1.01) {
    cout <<"He0 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_He0] <<"\n";
  }
  if (p_local[lv_He1]<-0.01 ||  p_local[lv_He1]>1.01) {
    cout <<"He1 fraction out of bounds to mpv9_HHe: ";
    cout << p_local[lv_He1] <<"\n";
  }
#endif
  p_out[pv_H0]  = p_local[lv_H0];
  p_out[pv_H0]  = max(Min_Nfrac, min(Max_Nfrac, p_out[pv_H0]));
  p_out[pv_He0] = p_local[lv_He0];
  p_out[pv_He0] = max(Min_Nfrac, min(Max_Nfrac*X_HE, p_out[pv_He0]));
  p_out[pv_He1] = p_local[lv_He1];
  p_out[pv_He1] = max(Min_Nfrac, min(Max_Nfrac*X_HE, p_out[pv_He1]));

  //
  // output pressure (cvode should ensure that we don't cool/heat
  // outside the allowed temperature ranges, so we only check if we
  // are debugging.
  //
  p_out[PG]    = p_local[lv_E]*(gamma_minus_one);

#ifdef MPV3_DEBUG
  double T = get_temperature(p_local);
  if (T>1.0001*EP->MaxTemperature) {
    Set_Temp(p_out,EP->MaxTemperature,0);
    cout <<"mpv9_HHe::convert_local2prim() HIGH T. ";
    cout <<"T="<<T<<", obtained from nH="<<nH<<", eint=";
    cout <<p_local[lv_E]<<", x="<<p_out[pv_Hp]<<"... ";
    cout <<" limiting to T="<<EP->MaxTemperature<<"\n";
  }
  if (T<0.9999*EP->MinTemperature) {
    Set_Temp(p_out,EP->MinTemperature,0);
    cout <<"mpv9_HHe::convert_local2prim() LOW  T. ";
    cout <<"T="<<T<<", obtained from nH="<<nH<<", eint=";
    cout <<p_local[lv_E]<<", x="<<p_out[pv_Hp]<<"... ";
    cout <<" limiting to T="<<EP->MaxTemperature<<"\n";
  }
#endif // MPV3_DEBUG

#ifdef TESTING
//cout <<"nH="<< p_local[lv_nH] <<", xp="<< p_out[pv_Hp] <<", ntot=";
//cout <<get_ntot(p_local[lv_nH],p_out[pv_Hp])<<"\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



void mpv9_HHe::interpret_radiation_data(
          const int N_heat,
          const std::vector<struct rt_source_data> &heat_src,
          const int N_ion,
          const std::vector<struct rt_source_data> &ion_src
          )
{
  //
  // For the first ionising source we have:
  //
  // Tau[0] = Int nH*y(H0)*sigma0(H0) dr
  // Tau[1] = Int nH*y(He0)*sigma0(He0) dr
  // Tau[2] = Int nH*y(He1)*sigma0(He1) dr
  // Tau[3] = Int nH*sigma(dust) dr
  // Similarly for dTau[iT]
  //
  // So that is the only source we parse, because that is all the
  // information we need.
  //

  if (N_ion<1) {
    cerr <<"No rad\'n source for mpv9_HHe::interpret_radiation_data";
    cerr <<".  setting all optical depths to zero.\n";
    for (short unsigned int iT=0; iT<MAX_TAU; iT++) {
      tau[iT] = 0.0;
      dtau[iT]= 0.0;
    }
    Vshell  = 1.0e100;
    return;
  }

  for (short unsigned int iT=0; iT<ion_src[0].NTau; iT++) {
    tau[iT] = ion_src[0].Column[iT];
    dtau[iT]= 0.0; //ion_src[0].DelCol[iT]; dtau is unused.
  }
  Vshell  = ion_src[0].Vshell;
  dS      = ion_src[0].dS;
}



// ##################################################################
// ##################################################################



double mpv9_HHe::get_temperature(
    const double *P   ///< local vector
    )
{
  //
  // E = nkT/(g-1), so T=(g-1)E/nk
  return gamma_minus_one*P[lv_E]/(get_ntot(P)*kB());
}



// ##################################################################
// ##################################################################



double mpv9_HHe::get_ntot(
    const double *P   ///< local vector
    )
{
  return nH*(X_H + X_HE) + get_ne(P);
}



// ##################################################################
// ##################################################################



double mpv9_HHe::get_ne(
    const double *P   ///< local vector
    )
{
  return nH*(X_H-P[lv_H0] + 2.0*(X_HE-P[lv_He0]) - P[lv_He1]);
}



// ##################################################################
// ##################################################################




// ##################################################################
// ##################################################################


int mpv9_HHe::ydot(
          double,               ///< current time (UNUSED)
          const N_Vector y_now, ///< current Y-value
          N_Vector y_dot,       ///< vector for Y-dot values
          const double *        ///< extra user-data vector (UNUSED)
          )
{

  double pir[3];
  double phr[3];
  double temp;

  //
  // Set pointers to y, ydot arrays and get electron number density
  // and gas temperature.
  //
  double *y = NV_DATA_S(y_now);
  double *ydot = NV_DATA_S(y_dot);
  double ne = get_ne(y);
  double T = get_temperature(y);
  double y_He2 = std::max(0.0, std::min(X_HE, X_HE -y[lv_He0] -y[lv_He1]));
#ifdef MPV9_DEBUG
  cout <<"YDOT: "; rep.printVec("Ynow",y,4);
  cout <<"YDOT: T="<<T<<", y(HE++)="<<X_HE -y[lv_He0] -y[lv_He1]<<"\n";
#endif // MPV9_DEBUG

  //for (short unsigned int ie=0; ie<N_equations; ie++) ydot[ie]=0.0;

  //
  // We set a minimum electron density based on the idea that Carbon
  // is singly ionised in low density gas.  y(C)=1.5e-4 in the gas
  // phase (by number) (Sofia,1997), approximately, so I add this to
  // the electron number density with an exponential cutoff at high
  // densities.
  //
  ne += nH*1.5e-4*metallicity*exp(-nH*1.0e-4);

  //
  // use current optical depths through cell to get accurate ph-ion
  // rates.
  //
  dtau[0] = nH*y[lv_H0 ]*dS*get_th_xsection(ION_H_N);
  dtau[1] = nH*y[lv_He0]*dS*get_th_xsection(ION_HE_N);
  dtau[2] = nH*y[lv_He1]*dS*get_th_xsection(ION_HE_P);
  dtau[3] = nH*     dS*get_th_xsection(ION_DUST);
#ifdef MPV9_DEBUG
  cout <<"YDOT: "; rep.printVec("dTau",dtau,4);
#endif // MPV9_DEBUG

  // ----------------------------------------------------------------
  //
  // photo-ionisation rates.  All rates are per H nucleon, which
  // needs no modification because the species fractions are per
  // H nucleon also (unlike Friedrich+,2012).
  //
  HHe_photoion_rate(tau[0],dtau[0],tau[1],dtau[1],tau[2],dtau[2],
                    nH, Vshell,
                    pir, phr);
  //temp = exp(-tau[3]); // attenuation by grey dust.
  ydot[lv_H0]  = -pir[0]; // *temp;
  ydot[lv_He0] = -pir[1]; // *temp;
  ydot[lv_He1] = (pir[1]-pir[2]); // *temp;
  ydot[lv_E]   = phr[0] +phr[1] +phr[2]; // *temp;
#ifdef MPV9_DEBUG
  cout <<"** YDOT: nH="<<nH<<", Vshell="<<Vshell<<", ds="<<dS<<"\n"; 
  cout <<"** YDOT: ionising: H0="<<-pir[0]<<", He0="<<-pir[1]<<", He1="<<-pir[2]<<"\n";
  cout <<"** YDOT: heating:  H0="<<phr[0]<<", He0="<<phr[1]<<", He1="<<phr[2]<<"\n";
  //
  // Check for NAN/INF
  //
  for (size_t v=0;v<nvl;v++) {
    if (!isfinite(ydot[v]))
      rep.error("INF/NAN in mpv9_HHe ydot()",v);
  }
#endif // MPV9_DEBUG
  //
  // ----------------------------------------------------------------

#ifdef MPV9_REC
  // ----------------------------------------------------------------
  //
  // Recombination rates.  Here we need only alpha_B for H0, but more
  // for He0 (alpha_A, alpha_1, alpha_B) and He+ (alpha_A, alpha_B,
  // alpha_1, alpha_2).  He+ has dielectronic recombination too, so
  // we can't just use the radiative rates.
  //
  // This scheme is the coupled on-the-spot approximation of 
  // Friedrich+,2012.
  //
  double a1=0.0, a2=0.0, aB=0.0, aA=0.0; // recombination coeffs.
  //
  // First H+ recombinations.
  //
  aB = H1_caseB;
  ydot[lv_H0] += aB*ne*(1.0-y[lv_H0]);

  //
  // Now He+ recombinations, affecting the H0 and He0 fractions.
  // - f0 = fraction of He+ case-1 photons that ionise H0.
  // - f1 = fraction of He+ case-1 photons that ionise He0.
  // -  p = fraction of He+ case-B photons that can ionise H0.
  //
  // Uncoupled OTS has f1=1, f0=0, p=0.
  //
  a1 = He1_case1;
  aB = He1_caseB;
  aA = He1_caseA;
  tau_frac(REGION_B,dtau[0],dtau[1],dtau[2],pir); // use pir[] for f
  phr[0] = p_factor;  // use phr[0] for "p".
  ydot[lv_H0]  -= ne*y[lv_He1]*(a1_He1*pir[0] +aB_He1*phr[0]);
  ydot[lv_He0] += (aA_He1-a1_He1*pir[1])*ne*y_He2;

  //
  // Now He++ recombinations, affecting H0, He0, and He1 fractions.
  // - g0/g1/g2 = fraction of case-1 photons ionising H0/He0/He1.
  // - v = fraction of case-B recombs that have 2-photon emission.
  // - l = fraction of 2-photon photons that can ionise H0 (1.425)
  // - m = fraction of 2-photon photons that can ionise He0 (0.737)
  // - 1-v = fraction of case-B recombs that lead to He-Ly-alpha.
  // - fe = fraction of He-Ly-alpha photons that are absorbed locally.
  // - z/(1-z) = frac. of aborbed He-Ly-alpha photons ionising H0/He0
  // For a single bin in Region B, z=f0, 1-z=f1
  //
  // Uncoupled OTS has g0=0, g1=0, g2=1, v=1, l=0, m=0
  //
  tau_frac(REGION_C,dtau[0],dtau[1],dtau[2],phr);  // use phr[] for g
#define FE 1.0
#define EL 1.425
#define EM 0.737
  //
  // temp=v. Fit to table V in Hummer & Seaton (1964,MN,127,217)
  //
  temp = 0.1764 +2.56074e-3*pow(log10(T),2.70387);
  a1 = He2_case1;
  a2 = He2_case2;
  aA = He2_caseA;
  aB = He2_caseB;
  //
  // H0 loses to recombs to ground state, recombs to n=2, and also
  // from 2-photon case-B recombinations.
  //
  ydot[lv_H0]  -= ne*y_He2*(a1*phr[0] +a2
                  +aB*(temp*(EL-EM+EM*pir[0])+(1.0-temp)*FE*pir[0]));
  //
  // He0 loses to recombs to the ground state and 2-photon case-B
  // recombinations.  Store in temp.var b/c there is a symmetric
  // term in ydot(He1).
  //
  temp = a1*phr[1] +aB*(temp*EM*pir[1]+(1.0-temp)*FE*pir[1]);
  ydot[lv_He0] -= ne*y_He2*(temp);
  //
  // He1 gains from caseA recombinations, loses from a fraction g2 of
  // recombs to the ground state causing re-ionisation, and gains
  // from ionisation of He0 from recombination photons.
  //
  ydot[lv_He1] += ne*y_He2*(aA -a1*phr[2] +temp);


  //
  // ----------------------------------------------------------------
#endif // MPV9_REC

  //
  // now multiply Edot by nH to get units of energy loss/gain per
  // unit volume per second.
  //
  ydot[lv_E] *= nH;

  //
  // We want to limit cooling as we approach the minimum temperature,
  // so we scale the rate to linearly approach zero as we reach Tmin.
  //
  if (ydot[lv_E]<0.0 && T<2.0*EP->MinTemperature) {
    ydot[lv_E] = min(0.0, (ydot[lv_E])*(T-EP->MinTemperature)
                                    /SimPM.EP.MinTemperature);
  }
  return 0;
}


#endif //  if not EXCLUDE_MPV9


