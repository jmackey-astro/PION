///
/// \file MPv10.cpp
/// \author Maggie Goulden
/// \date 2018.10
///
/// Description:
/// - multi-species ionization/recombination non-equilibrium
///   chemistry solver.
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from 
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.
///
/// Modifications:
/// - 2018.10.09 JM: edited header
///
///

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ================================================================
// ================================================================



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


// ##################################################################
// ##################################################################


#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include <set>
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/MPv10.h"

using namespace std;

//#define MPv10_DEBUG


//
// Timestep-limiting is important for making chemistry consistent
// with hydrodynamics.
// This is a good value for MPv3 (see Mackey,2012,A&A,539,A147)
// Will have to do some tests for MPv10...
//
#define DTFRAC 0.25


// ##################################################################
// ##################################################################


void MPv10::get_error_tolerances(
      double *reltol, ///< relative error tolerance.
      double atol[] ///< absolute error tolerances
      )
{
  *reltol = JM_RELTOL;
  atol[lv_H0]   = JM_MINNEU; ///< minimum neutral fraction I care about.
  atol[lv_eint] = JM_MINERG; ///< E_{int}: for n=1.0, T=1.0e4, ==> E=2.07e-12, so say 1e-17?
  return;
}


// ##################################################################
// ##################################################################


void MPv10::get_problem_size(
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


MPv10::MPv10(
      const int nd,   ///< grid dimensions
      const int csys,   ///< Coordinate System flag
      const int nv,   ///< Total number of variables in state vector
      const int ntracer,  ///< Number of tracer variables in state vector.
      const std::string *tracers,  ///< List of what the tracer variables mean.
      struct which_physics *ephys, ///< pointer to extra physics flags.
      struct rad_sources *rsrcs,   ///< radiation sources.
      const double g  ///< EOS Gamma
      )
: microphysics_base(ephys,rsrcs),
  ndim(nd), nv_prim(nv), eos_gamma(g), coord_sys(csys)
{
  
  k_B = pconst.kB();  // Boltzmann constant.
  m_p = pconst.m_p(); // Proton mass.
  m_H = pconst.m_H(); // Hydrogen mass.
  m_He = pconst.m_He(); // Helium mass.

  cout <<"\n---------------------------------------------------------------------\n";
  cout <<"MPv10: a microphysics class.\n";

  // ----------------------------------------------------------------
  // ------- Set up tracer variables:                       ---------
  // ------- 1) Identify elements present in tracer list    ---------
  // ------- 2) Record X_elem_index vector, N_elem          ---------
  // ------- 3) Record: y_ion_index (index in primitive vector),
  //            y_elem_mass_frac (mass frac of y_ion_index spec elem),
  //            y_ion_num_elec (# electrons of y_ion_index species),
  //            N_species_by_elem (arranged to match X_elem_index)
  //            
  // ----------------------------------------------------------------
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  int len = ntracer;
  int ftr = nv_prim -ntracer; // first tracer variable.
  lv_y0_offset = ftr;
  string s; pv_H1p=-1;
  
  // 1) Identify elements present in tracer list. Set keeps only unique objects; no doubling up.
  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    if (s.substr(0,2)=="He"){
      set_elem.insert("He");
    }
    else if (s[0]=='H'){
      set_elem.insert("H");
    }
    else if (s[0]=='C'){
      set_elem.insert("C");
    }
  }
  
  set<string>::iterator it; /// < iterator for set_elem
  // 2) Record X_elem_index vector, N_elem
  N_elem = 0;
  for (it = set_elem.begin(); it != set_elem.end(); ++it) {
      X_elem_index.push_back(ftr + N_elem); ///<record primitive vector index of each element
      N_species_by_elem.push_back(0);
      N_elem++;
  }
  
  // 3) Establish N_species_by_elem, N_species, y_ion_index, y_ion_num_elec, y_elem_mass_frac
  N_species=0;
  int elem_counter=0;
  float mass;
  for (it = set_elem.begin(); it != set_elem.end(); ++it) {
    // Define mass fraction / atomic mass (grams) for this element first.
    if((*it)=="H"){
      //mass_frac = EP->H_MassFrac;
      mass = m_H;
      X_atom_mass.push_back(mass);
      H_ion_index.push_back(0);
    }
    else if ((*it)=="He"){
      //mass_frac = EP->Helium_MassFrac;
      mass = m_He;
      X_atom_mass.push_back(mass);
      He_ion_index.push_back(0); He_ion_index.push_back(0);
    }
          
    // Loop over every tracer and assign species index / mass fraction / num electrons to vectors if that tracer corresponds to the current element.
    for (int i=0;i<len;i++) {
      s = tracers[i];
      // He compared separately as it is the first two characters.
      if (s.substr(0,2)=="He" & (*it)=="He"){
        N_species_by_elem[elem_counter]++;
        y_ion_index.push_back(ftr + N_elem + N_species);
        y_ion_mass_frac_index.push_back(X_elem_index[elem_counter]);
        y_elem_atom_mass.push_back(mass);
        int num_elec_int;
        stringstream ss; ss << s.substr(2,1); ss >> num_elec_int; //Use stringstream to convert string to int.
        y_ion_num_elec.push_back(num_elec_int);
        cout << s.substr(2,1) << ", " << num_elec_int <<"\n\n\n\n\n\n";
        He_ion_index[num_elec_int-1] = ftr + N_elem + N_species;
        N_species++;
      }
      else if (s.substr(0,1)==(*it) & s.substr(0,2)!="He" & s.substr(0,1)=="H"){
        N_species_by_elem[elem_counter]++;
        y_ion_index.push_back(ftr + N_elem + N_species);
        y_ion_mass_frac_index.push_back(X_elem_index[elem_counter]);
        y_elem_atom_mass.push_back(mass);
        int num_elec_int;
        stringstream ss; ss << s.substr(1,1); ss >> num_elec_int; //Use stringstream to convert string to int.
        y_ion_num_elec.push_back(num_elec_int);
        H_ion_index[num_elec_int-1] = ftr + N_elem + N_species;
        N_species++;
      }
    }
    elem_counter++;
  }
  cout << "N_species=" << N_species << ", N_elements=" << N_elem << "\n";  
  
  
  // ================================================================
  // ================================================================
#ifdef TESTING
  cout <<"MPv10:: EP and RS: "<<EP<<"\t"<<RS<<endl;
#endif

  // ----------------------------------------------------------------
  // --- Set up local variables: ion fraction and internal energy density.
  // ----------------------------------------------------------------

  // Get the mean mass per H atom from the He and Z mass fractions.
  // NOTE \Maggie{or let's not assume that...} Assume metal content is low enough to ignore it.
  double X = 1.0-EP->Helium_MassFrac;
  mean_mass_per_H = m_p/X;
  METALLICITY = EP->Metal_MassFrac/0.0142; // in units of solar.
  cout <<"Metallicity = "<<METALLICITY<<" of solar (0.0142)\n";

  setup_local_vectors();
  gamma_minus_one = eos_gamma -1.0;
  Min_NeutralFrac     = JM_MINNEU;
  Max_NeutralFrac     = 1.0-JM_MINNEU;

  //
  // initialise all the radiation variables to values that limit their heating
  // and cooling abilities.
  //
  mpv_nH   = 1.0;  // Hydrogen number density (density of H+ and H)
  mpv_Vshell  = 1.0e54;
  mpv_Tau0     = 1.0e6;
  mpv_dTau0    = 1.0;
  mpv_G0_UV   = 0.0;
  mpv_G0_IR   = 0.0;
  mpv_delta_S = 0.0;
  mpv_NIdot   = 0.0;
  ion_src_type = -1;
  N_ion_srcs  = 0;
  N_diff_srcs = 0;
  // ================================================================
  // ================================================================
  
  // ----------------------------------------------------------------
  // We want to set up the CIE cooling function for metals-only from WSS09
  // (i.e. with H+He cooling subtracted out).
  // ----------------------------------------------------------------
  setup_WSS09_CIE_OnlyMetals();
  // ================================================================
  // ================================================================


  // ----------------------------------------------------------------
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  setup_cvode_solver_without_Jacobian();
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // ---------- output cooling rates for various temperatures -------
  // ----------------------------------------------------------------
  
  //double p[nv_prim];
  //p[RO]=2.338e-24; p[PG]=1.0e-12;
  for (int i=0;i<N_species;i++){
    float n_X_y = 1.0;//p[RO]*( p[ y_ion_mass_frac_index[i]] / y_elem_atom_mass[i]);
    y_ion_number_density.push_back(n_X_y);  //X number density at current cell, organised to match y_ion.
  }
  for (int i=0;i<N_elem;i++){
    float n_X = 1.0;//p[RO]*( p[ X_elem_index[i]] / X_atom_mass[i]);
    X_elem_number_density.push_back(n_X);  //X number density at current cell, organised to match y_ion.
  }

  // ================================================================
  // ================================================================

  cout <<"MPv10: Constructor finished and returning.\n";
  cout <<"---------------------------------------------------------------------\n\n";
  return;
}



// ##################################################################
// ##################################################################


void MPv10::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  nvl     = N_species +1;    // two local variables to integrate
  N_extradata = 0;
  N_equations = N_elem;
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
  lv_H0   = 0;    // x(H0) is the first element in the array NOTE \Maggie{ LEGACY CODE; REMOVE LATER.}
  lv_eint = 1;    // E_{int} is the second element.
  lv_eint = N_species;
  //cout<<"!!!!!!!!!!!!!!!!!! nvl="<<nvl<<"\n";
  return;
}



// ##################################################################
// ##################################################################


MPv10::~MPv10()
{
  //
  // Free vector memory
  //
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}


// ##################################################################
// ##################################################################



int MPv10::Tr(const string s)
{
  if (s.substr(0,2)=="He"){
    int num_elec_int; 
    stringstream ss; ss << s.substr(2,1); ss >> num_elec_int; //Use stringstream to convert string to int.
    int ion_index = He_ion_index[num_elec_int-1];
    return ion_index;
  }
  else if (s.substr(0,1)=="H"){
    int num_elec_int;
    stringstream ss; ss << s.substr(1,1); ss >> num_elec_int; //Use stringstream to convert string to int.
    int ion_index = H_ion_index[num_elec_int-1];
    return ion_index;
  }
  else { return-1;}
    
}




// ##################################################################
// ##################################################################


double MPv10::get_temperature(
      double *y_ion_frac, //< y_ion_fraction (by y_ion_index)
      vector<pion_flt>& y_ion_number_density,//const double nH, ///< nH
      const double E ///< E_int (per unit volume)
      )
{
  //
  // returns gas temperature according to E=nkT/(g-1), => T = E*(g-1)/(K*n)
 return gamma_minus_one*E/(k_B*get_ntot(y_ion_frac,y_ion_number_density));
}

// ##################################################################
// ##################################################################


double MPv10::get_ntot(
      double *y_ion_frac,//<y_ion_fraction (by y_ion_index)
      vector<pion_flt>& y_ion_number_density//const double nH, ///< nH
      )
{
  int species_counter = 0;
  pion_flt n_tot=0;
  
  for (int i=0;i<N_elem;i++){//loop over every element
    int N_elem_species=N_species_by_elem[i];
    pion_flt neutral_frac = 1; //neutral frac, got by subtracting off the ion fractions in the next loop.
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      pion_flt number_density = y_ion_number_density[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      
      n_tot += (1+num_elec)*y_ion_frac[species_counter]*number_density; //add on number of particles got by electrons + ionised atom
      
      neutral_frac -= y_ion_frac[species_counter];
      species_counter ++;
    }
    n_tot += neutral_frac*y_ion_number_density[species_counter-1];
  }
  return n_tot;
}





// ##################################################################
// ##################################################################


int MPv10::convert_prim2local(
      const pion_flt *p_in, ///< primitive vector from grid cell (length nv_prim)
      double *p_local
      )
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  // NOTE \Maggie{ Should also include number densities of all other elements, I think}
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);
  //p_local[lv_H0]   = 1.0-p_in[pv_H1p];
  // loop over N_species instead of just having lv_H0, to identify y(species). Want first element to occur at p_local[0].
  int local_index=-1;
  for (int i=0;i<N_species;i++){
    float n_X_y = p_in[RO]*( p_in[ y_ion_mass_frac_index[i]] / y_elem_atom_mass[i]);
    y_ion_number_density[i] = n_X_y;
    local_index = y_ion_index[i] - lv_y0_offset;
    p_local[ local_index] = p_in[ y_ion_index[i]];
  }
  //mpv_nH = p_in[RO]/mean_mass_per_H;


#ifdef MPv10_DEBUG
  //
  // Check for negative ion fraction, and set to a minimum value if found.
  //
  if (p_local[lv_H0]>1.0) {
    cout <<"MPv10::convert_prim2local: negative ion fraction input: ";
    cout <<"x(H0)="<<p_local[lv_H0] <<", resetting to "<<Max_NeutralFrac<<"\n";
    p_local[lv_H0] = Max_NeutralFrac;
  }
  //
  // Check for bad values:
  //
  if (p_local[lv_H0]>1.01 || p_local[lv_H0]<-0.01) {
    cout <<"MPv10::convert_prim2local: bad ion fraction: ";
    cout <<"x(H0)="<<p_local[lv_H0] <<", resetting to [0,1]\n";
  }
#endif

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  //
  // Check for negative pressure (note this shouldn't happen, so we output a
  // warning) and set to 10K if we find it.
  //
  if (p_local[lv_eint]<=0.0) {
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_prim2local: negative pressure input: p=";
    cout <<p_local[lv_eint]<<", setting to "<<EP->MinTemperature<<"K.\n";
#endif
    pion_flt y_ion_frac[N_species];
    for (int s=0;s<N_species;s++){
      y_ion_frac[s] = p_in[y_ion_index[s]];
    }
    p_local[lv_eint] = get_ntot(y_ion_frac,y_ion_number_density)*k_B*EP->MinTemperature/(gamma_minus_one);
  }


#ifdef MPv10_DEBUG
  //
  // Check for NAN/INF
  //
  for (int v=0;v<2;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics",p_local[v]);
  }
  if (mpv_nH<0.0 || !isfinite(mpv_nH))
    rep.error("Bad density input to MPv10::convert_prim2local",mpv_nH);
#endif // MPv10_DEBUG
  
  return 0;
}


// ##################################################################
// ##################################################################


int MPv10::convert_local2prim(
      const double *p_local,
      const pion_flt *p_in, ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *p_out       ///< updated primitive vector for grid cell (length nv_prim)
      )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];

  p_out[PG]    = p_local[lv_eint]*(gamma_minus_one);
  //p_out[pv_H1p] = 1.0-p_local[lv_H0];

#ifdef MPv10_DEBUG
  if (p_out[pv_H1p]<-10.0*JM_RELTOL || p_out[pv_H1p]>1.0*(1.0+JM_RELTOL) || !isfinite(p_out[pv_H1p])) {
    rep.printVec("p_in",p_in, nv_prim);
    rep.printVec("p_out",p_out, nv_prim);
    rep.printVec("p_local",p_local, nvl);
    rep.error("Bad output H+ value in MPv10::convert_local2prim",p_out[pv_H1p]-1.0); 
  }
  if (p_out[PG]<0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in MPv10::convert_local2prim",p_out[PG]);
#endif // MPv10_DEBUG

  //
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  //p_out[pv_H1p] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[pv_H1p])));

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected x(H+) from p_out[]).
  //
  pion_flt y_ion_frac[N_species];
  for (int s=0;s<N_species;s++){
    y_ion_frac[s] = p_out[y_ion_index[s]];
  }
  double T = get_temperature(y_ion_frac, y_ion_number_density, p_local[lv_eint]);
  if (T>1.01*EP->MaxTemperature) {
    Set_Temp(p_out,EP->MaxTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() HIGH temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MaxTemperature<<"\n";
#endif // MPv10_DEBUG
  }
  if (T<0.99*EP->MinTemperature) {
    Set_Temp(p_out,EP->MinTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() LOW  temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MinTemperature<<"\n";
#endif // MPv10_DEBUG
  }

  return 0;
}



// ##################################################################
// ##################################################################



double MPv10::Temperature(
      const pion_flt *pv, ///< primitive vector
      const double      ///< eos gamma
      )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  // NOTE \Maggie{ This needs changing, too, as the vector only contains hydrogen.}
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    //cout <<"MPv10::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];
  convert_prim2local(pv,P);
  pion_flt y_ion_frac[N_species];
  for (int s=0;s<N_species;s++){
    y_ion_frac[s] = pv[y_ion_index[s]];
  }
  return (get_temperature(y_ion_frac, y_ion_number_density, P[lv_eint]));
}



// ##################################################################
// ##################################################################



int MPv10::Set_Temp(
      pion_flt *p_pv,   ///< primitive vector.
      const double T, ///< temperature
      const double    ///< eos gamma.
      )
{
  //
  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.
  //
  if (p_pv[PG]<=0.0) {
    //cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p_pv[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }
  double P[nvl];
  int err = convert_prim2local(p_pv,P);
  pion_flt y_ion_frac[N_species];
  for (int s=0;s<N_species;s++){
    y_ion_frac[s] = p_pv[y_ion_index[s]];
  }
  P[lv_eint] = get_ntot(y_ion_frac,y_ion_number_density)*k_B*T/(gamma_minus_one);
  err += convert_local2prim(P, p_pv, p_pv);
  return err;
}


// ##################################################################
// ##################################################################



int MPv10::TimeUpdateMP(
      const pion_flt *p_in,   ///< Primitive Vector to be updated.
      pion_flt *p_out,        ///< Destination Vector for updated values.
      const double dt,      ///< Time Step to advance by.
      const double,         ///< EOS gamma.
      const int,            ///< Switch for what type of integration to use.
      double *random_stuff  ///< Vector of extra data (column densities, etc.).
      )
{
  int err=0;
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to MPv10::TimeUpdateMP()",err);
  }
  // Populates CVODE vector with initial conditions (input)
  for (int v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];

  //
  // Calculate y-dot[] to see if anything is changing significantly over dt
  //
  double maxdelta=0.0;
  err = ydot(0, y_in, y_out, 0);
  if (err) 
    rep.error("dYdt() returned an error in MPv10::TimeUpdateMP_RTnew()",err);
  for (int v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out,v)*dt/NV_Ith_S(y_in,v)));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < EULER_CUTOFF) {
    for (int v=0;v<nvl;v++) {
      NV_Ith_S(y_out,v) = NV_Ith_S(y_in,v) + dt*NV_Ith_S(y_out,v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      rep.error("integration failed: MPv10::TimeUpdateMP_RTnew()",err);
    }
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v=0;v<nvl;v++) P[v] = NV_Ith_S(y_out,v);
  err = convert_local2prim(P,p_in,p_out);

  return err;
}



// ##################################################################
// ##################################################################



double MPv10::timescales(
      const pion_flt *p_in, ///< Current cell state vector.
      const double,   ///< EOS gamma.
      const bool, ///< set to 'true' if including cooling time.
      const bool, ///< set to 'true' if including recombination time.
      const bool  ///< set to 'true' if including photo-ionsation time.
      )
{
#ifdef MPv10_DEBUG
  if (RS->Nsources!=0) {
    cout <<"WARNING: MPv10::timescales() using non-RT version!\n";
  }
#endif // MPv10_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin= timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
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
double MPv10::timescales_RT(
      const pion_flt *p_in, ///< Current cell state vector.
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
    rep.error("Bad input state to MPv10::timescales_RT()",err);
  }
  NV_Ith_S(y_in,lv_H0  ) = P[lv_H0];
  NV_Ith_S(y_in,lv_eint) = P[lv_eint];

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in MPv10::timescales_RT()",err);
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
#ifdef USE_RELATIVE_NEUFRAC_DTLIMIT
  t = min(t,DTFRAC*max(5.0e-2,NV_Ith_S(y_in, lv_H0))/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"using neutral fraction.\n";
#else
  t = min(t,DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"limit by dx: dt="<<DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\n";
#endif


#ifdef MPv10_DEBUG
  if (t<3.16e9) {
  cout <<"MP timescales: xdot="<<NV_Ith_S(y_out, lv_H0);
  cout <<", Edot="<<NV_Ith_S(y_out, lv_eint)<<" t_x="<<t;
  }
#endif // MPv10_DEBUG

#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  //
  // Now cooling/heating time to limit to X% change in energy).
  //
  t = min(t,DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE));
  //cout <<"using fractional energy change.\n";
  //cout <<"limit by dE: dt="<<DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE)<<"\n";
#endif

#ifdef MPv10_DEBUG
  if (t<3.16e9) {
  cout <<" and min(t_x,t_e)="<<t<<",  "; rep.printVec("P[1-x,E]",P,nvl);
  }
#endif // MPv10_DEBUG
  return t;
}




// ##################################################################
// ##################################################################



int MPv10::ydot(
      double,               ///< current time (UNUSED)
      const N_Vector y_now, ///< current Y-value
      N_Vector y_dot,       ///< vector for Y-dot values
      const double *        ///< extra user-data vector (UNUSED)
      )
{
  //
  // fixes min-neutral-fraction to Min_NeutralFrac
  // NOTE \Maggie { I thiiiink quite a few things need changing here...}
  double OneMinusX = max(NV_Ith_S(y_now,lv_H0),Min_NeutralFrac); //y0
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double x_in      = 1.0-OneMinusX;
    
  // NOTE \Maggie{this should be updated so as to add up all the electrons freed from ions}
  //double ne        = JM_NELEC*x_in*mpv_nH;
  //
  // First get the temperature.  We assume the total particle number density
  // is given by 1.1*nH*(1+x_in), appropriate for a gas with 10% Helium by 
  // number, and if He is singly ionised whenever H is.
  //
  double ne=0;
  double y_ion_frac[N_species];
  
  for (int s=0;s<N_species;s++){
    int this_ion_index = y_ion_index[s] - lv_y0_offset; //basically, lv_H1 (equivalent to 1 - y[lv_H0])
    int this_elem_mass_frac_index = y_ion_mass_frac_index[s] - lv_y0_offset;
    y_ion_frac[s] = NV_Ith_S(y_now,this_ion_index);
    double X_elem_mass_frac = NV_Ith_S(y_now,this_elem_mass_frac_index);
    ne += y_ion_frac[s]*y_ion_num_elec[s]*X_elem_mass_frac;
  }
  double T = get_temperature(y_ion_frac, y_ion_number_density, E_in);
  cout << "Temperature="<< T<<"\n"<<"Electron density="<<ne<<"\n\n\n";
  //double T = get_temperature_local(mpv_nH, E_in, x_in);


  double temp1=0.0, temp2=0.0;
  double oneminusx_dot=0.0; // oneminusx_dot is in units of 1/s
  double Edot=0.0;
  // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH
  // at the end.

  //
  // We set a minimum electron density based on the idea that Carbon is singly
  // ionised in low density gas.  y(C)=1.5e-4 in the gas phase (by number)
  // (Sofia,1997), approximately, so I add this to the electron number density
  // with an exponential cutoff at high densities.
  // NOTE \Maggie{ this definitely needs a modification}
  //ne += mpv_nH*1.5e-4*METALLICITY*exp(-mpv_nH/1.0e4);


  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0
  // NOTE \Maggie {again, lookup table instead of calculation}
  Hi_coll_ion_rates(T, &temp1, &temp2);
  oneminusx_dot -= temp1*ne*OneMinusX; // the nH is divided out on both sides.
  Edot -= temp2*ne*OneMinusX;
  //cout <<"CI-CR="<< temp2*ne*OneMinusX<<"\n";

  //
  // radiative recombination of H+
  //
  oneminusx_dot += Hii_rad_recomb_rate(T) *x_in*ne;
  //
  // Total H+ cooling: recombination plus free-free
  //
  Edot -= Hii_total_cooling(T) *x_in*ne;
  //cout <<"HII-TC="<<Hii_total_cooling(T) *x_in*ne<<"\n";

  //
  // Add Helium free-free (Z^2*n(He)/n(H) = 0.25*X(He)/X(H) of the H+ free-free
  // rate) The normalisation is scaled so that I multiply by ne*nHp to get the
  // correct cooling rate (i.e. the abundance of He is included in the
  // prefactor).
  //
#ifndef HE_INERT
  // Only if He is ionised, otherwise it has no free-free.
  //NOTE: \Maggie{ JM_NION needs to go as it's just a hack. Come back to this.}
  //Edot -= 1.68e-27*(JM_NION-1.0)*sqrt(T)*x_in*ne;
#endif // HE_INERT

  //
  // collisional excitation cooling of H0 Aggarwal (1983) and Raga+(1997,ApJS).
  //
  Edot -= Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10);
  //cout <<"CE-CR="<<Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10)<<"\n";


  //
  // Cosmic ray heating (HAdCM09 eq.A7).
  //
  Edot += 5.0e-28*OneMinusX;
  //cout <<"CR-HR="<<5.0e-28*OneMinusX<<"\n";

  //
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar neighbourhood.
  //
  oneminusx_dot -= 1.8e-17*OneMinusX;

  //
  // Now COOLING: First forbidden line cooling of e.g. OII,OIII, dominant in
  // HII regions.  This is collisionally excited lines of photoionised metals.
  // (HAdCM09 eq.A9) I have exponentially damped this at high temperatures!
  // This was important!  Oxygen abundance set to 5.37e-4 from
  // Asplund+(2009,ARAA), times 0.77 to account for 23% of O in solid dust.
  //
  temp1 = 1.20e-22*METALLICITY *exp(-33610.0/T -(2180.0*2180.0/T/T))
                               *x_in*ne*exp(-T*T/5.0e10);
  //
  // Fit to Raga, Mellema, Lundqvist (1997) rates for CNO if all are only
  // singly ionised, and for gas phase abundances of CNO of - n(C)/nH =
  // 2.95e-4*0.508 = 1.5e-4 (Sofia+,1997) - n(N)/nH = 7.41e-5 - n(O)/nH =
  // 5.37e-4*0.77 (0.23 goes in solids) These are taken from Asplund et al.
  // 2009.
  //
  //temp1 = 3.0e-22*METALLICITY*exp(-pow(1.4e5/T,0.6))
  //        *exp(-sqrt(mpv_nH/3.0e4)) *x_in*ne*exp(-T*T/5.0e10);


  //
  // Now the Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve.  We
  // take the actual cooling rate to be the max of SD93-CIE and the previous
  // two terms.
  //
  temp2 = cooling_rate_SD93CIE(T) *x_in*x_in*mpv_nH*METALLICITY;
  Edot -= max(temp1,temp2);

  NV_Ith_S(y_dot,lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot,lv_eint) = Edot;
  return 0;
}


