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
  // NOTE \Maggie{hardcoded this for now, but must change atol after.}
  atol[0] = JM_MINNEU; ///< minimum neutral fraction I care about.
  atol[1] = JM_MINERG; ///< E_int: for n=1.0, T=1.0e4, ==> E=2.07e-12, so say 1e-17?
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
  ndim(nd), nv_prim(nv), eos_gamma(g), coord_sys(csys),
  T_min(1e0), T_max(1e9), Num_temps(1e2), delta_log_temp((log10f(T_max) - log10f(T_min))/(Num_temps-1))
{ 
  /// ===================================================================
  ///  Initialise Temperature, Recombination (radiative + dielectronic),
  ///       Ionisation, and (recomb and ionisation) Slopes Tables
  ///  NOTE: SHOULD PROBABLY READ THESE FROM FILES EVENTUALLY... BUT I DON'T FEEL LIKE IT RIGHT NOW, AND I JUST WANT THEM TO WORK FIRST.
  /// ===================================================================

  float Temp_table[Num_temps] = { 1, 1.23285, 1.51991, 1.87382, 2.31013, 2.84804, 3.51119, 4.32876, 5.3367, 6.57933, 8.11131, 10, 12.3285, 15.1991, 18.7382, 23.1013, 28.4804, 35.1119, 43.2876, 53.367, 65.7933, 81.1131, 100, 123.285, 151.991, 187.382, 231.013, 284.804, 351.119, 432.876, 533.67, 657.933, 811.131, 1000, 1232.85, 1519.91, 1873.82, 2310.13, 2848.04, 3511.19, 4328.76, 5336.7, 6579.33, 8111.31, 10000, 12328.5, 15199.1, 18738.2, 23101.3, 28480.4, 35111.9, 43287.6, 53367, 65793.4, 81113.1, 100000, 123285, 151991, 187382, 231013, 284804, 351119, 432876, 533670, 657934, 811131, 1e+06, 1.23285e+06, 1.51991e+06, 1.87382e+06, 2.31013e+06, 2.84804e+06, 3.5112e+06, 4.32876e+06, 5.3367e+06, 6.57934e+06, 8.11131e+06, 1e+07, 1.23285e+07, 1.51991e+07, 1.87382e+07, 2.31013e+07, 2.84804e+07, 3.5112e+07, 4.32876e+07, 5.3367e+07, 6.57933e+07, 8.11131e+07, 1e+08, 1.23285e+08, 1.51991e+08, 1.87382e+08, 2.31013e+08, 2.84804e+08, 3.5112e+08, 4.32876e+08, 5.3367e+08, 6.57934e+08, 8.11131e+08, 1e+09};
  
  float recomb_rate_table[Num_temps][5] = {
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 3.41202e-10, 2.89621e-10, 2.45838e-10, 2.08674e-10, 1.77128e-10, 1.50351e-10, 1.27622e-10, 1.08329e-10, 9.19525e-11, 7.80517e-11, 6.62524e-11, 5.62368e-11, 4.77353e-11, 4.0519e-11, 3.43936e-11, 2.91942e-11, 2.47808e-11, 2.10346e-11, 1.78547e-11, 1.51556e-11, 1.28645e-11, 1.09197e-11, 9.26893e-12, 7.86771e-12, 6.67832e-12, 5.66874e-12, 4.81178e-12, 4.08436e-12, 3.46692e-12, 2.94281e-12, 2.49794e-12, 2.12032e-12, 1.79978e-12, 1.5277e-12, 1.29675e-12, 1.10072e-12, 9.34319e-13, 7.93075e-13, 6.73184e-13, 5.71416e-13, 4.85033e-13, 4.11709e-13, 3.4947e-13, 2.96639e-13, 2.51795e-13, 2.13731e-13, 1.8142e-13, 1.53994e-13, 1.30714e-13, 1.10954e-13, 9.41806e-14, 7.9943e-14, 6.78577e-14, 5.75995e-14, 4.8892e-14, 4.15008e-14, 3.5227e-14, 2.99016e-14, 2.53813e-14, 2.15443e-14, 1.82874e-14, 1.55228e-14, 1.31762e-14, 1.11843e-14, 9.49352e-15, 8.05836e-15, 6.84015e-15, 5.8061e-15, 4.92837e-15, 4.18333e-15, 3.55092e-15, 3.01412e-15, 2.55846e-15, 2.17169e-15, 1.84339e-15, 1.56472e-15, 1.32818e-15, 1.12739e-15, 9.56959e-16, 8.12292e-16, 6.89495e-16, 5.85262e-16, 4.96786e-16, 4.21685e-16, 3.57938e-16, 3.03827e-16, 2.57896e-16, 2.18909e-16, 1.85816e-16, 1.57726e-16, 1.33882e-16, 1.13642e-16, 9.64626e-17, 8.18801e-17, 6.9502e-17, 5.89951e-17, 5.00767e-17, 4.25064e-17, 3.60806e-17},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ 1.33093e-10, 1.17666e-10, 1.03994e-10, 9.18837e-11, 8.11612e-11, 7.16715e-11, 6.32763e-11, 5.58519e-11, 4.92884e-11, 4.34877e-11, 3.83627e-11, 3.38359e-11, 2.98384e-11, 2.63091e-11, 2.3194e-11, 2.04448e-11, 1.8019e-11, 1.5879e-11, 1.39914e-11, 1.23265e-11, 1.08584e-11, 9.56392e-12, 8.42264e-12, 7.41655e-12, 6.52971e-12, 5.74807e-12, 5.0592e-12, 4.45214e-12, 3.91723e-12, 3.44593e-12, 3.0307e-12, 2.6649e-12, 2.34268e-12, 2.05887e-12, 1.80891e-12, 1.58879e-12, 1.39495e-12, 1.2243e-12, 1.07406e-12, 9.41811e-13, 8.25422e-13, 7.23005e-13, 6.32898e-13, 5.5364e-13, 4.8394e-13, 4.22663e-13, 3.68807e-13, 3.2149e-13, 2.79936e-13, 2.43486e-13, 2.11883e-13, 1.87118e-13, 1.791e-13, 2.12928e-13, 3.25733e-13, 5.44478e-13, 8.59022e-13, 1.21617e-12, 1.54153e-12, 1.77147e-12, 1.87401e-12, 1.85128e-12, 1.72897e-12, 1.54243e-12, 1.32565e-12, 1.10534e-12, 8.9925e-13, 7.17141e-13, 5.62742e-13, 4.35843e-13, 3.34002e-13, 2.53772e-13, 1.91481e-13, 1.43671e-13, 1.0731e-13, 7.98579e-14, 5.92515e-14, 4.38562e-14, 3.23971e-14, 2.38937e-14, 1.7599e-14, 1.29486e-14, 9.51838e-15, 6.9916e-15, 5.13234e-15, 3.76548e-15, 2.76137e-15, 2.02421e-15, 1.48333e-15, 1.08664e-15, 7.95822e-16, 5.82695e-16, 4.26555e-16, 3.12194e-16, 2.28454e-16, 1.67149e-16, 1.22278e-16, 8.94407e-17, 6.54141e-17},
{ 5.39147e-10, 4.82303e-10, 4.31218e-10, 3.85324e-10, 3.4411e-10, 3.07114e-10, 2.73921e-10, 2.44153e-10, 2.17473e-10, 1.93575e-10, 1.72182e-10, 1.53045e-10, 1.35938e-10, 1.20656e-10, 1.07016e-10, 9.48513e-11, 8.40103e-11, 7.43571e-11, 6.57686e-11, 5.81336e-11, 5.13516e-11, 4.53321e-11, 3.99936e-11, 3.52626e-11, 3.10731e-11, 2.73656e-11, 2.4087e-11, 2.11895e-11, 1.86304e-11, 1.63716e-11, 1.43789e-11, 1.2622e-11, 1.10736e-11, 9.70986e-12, 8.50918e-12, 7.45259e-12, 6.5232e-12, 5.70605e-12, 4.98789e-12, 4.35698e-12, 3.80298e-12, 3.3167e-12, 2.89006e-12, 2.51593e-12, 2.18799e-12, 1.9007e-12, 1.64915e-12, 1.42904e-12, 1.23657e-12, 1.06839e-12, 9.21565e-13, 7.93493e-13, 6.81901e-13, 5.84779e-13, 5.00362e-13, 4.27095e-13, 3.63609e-13, 3.087e-13, 2.61305e-13, 2.20488e-13, 1.85424e-13, 1.55385e-13, 1.29728e-13, 1.07884e-13, 8.93538e-14, 7.36932e-14, 6.05118e-14, 4.94646e-14, 4.02482e-14, 3.25956e-14, 2.62727e-14, 2.10751e-14, 1.68247e-14, 1.33675e-14, 1.05704e-14, 8.31963e-15, 6.51816e-15, 5.08399e-15, 3.9482e-15, 3.05334e-15, 2.3518e-15, 1.80449e-15, 1.37947e-15, 1.0509e-15, 7.97955e-16, 6.04022e-16, 4.55896e-16, 3.43163e-16, 2.57653e-16, 1.92997e-16, 1.44251e-16, 1.076e-16, 8.01122e-17, 5.95447e-17, 4.41881e-17, 3.27448e-17, 2.4233e-17, 1.79121e-17, 1.32255e-17}
};
  
  
  k_B = pconst.kB();  // Boltzmann constant.
  m_p = pconst.m_p(); // Proton mass.
  m_H = pconst.m_H(); // Hydrogen mass.
  m_He = pconst.m_He(); // Helium mass.

  cout <<"\n---------------------------------------------------------------------\n";
  cout <<"MPv10: a microphysics class.\n";

  // -----------------------------------------------------------------------------
  // --------- Set up tracer variables:                                  ---------
  // --------- (i) Identify elements present in tracer list              ---------
  // --------- (ii) Record X_mass_frac_index vector, N_elem              ---------
  // --------- (iii) Record: y_ion_index_prim (index in primitive vector),    ---------
  // ---------       y_ion_num_elec (# electrons of y_ion_index_prim species),---------
  // ---------       N_species_by_elem (ordered like X_mass_frac_index)  ---------
  // ---------                                                           ---------
  // -----------------------------------------------------------------------------
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  
  int len = ntracer;
  int ftr = nv_prim -ntracer; // first tracer variable.
    
  string s; //pv_H1p=-1;
  
  N_elem = 0; N_species=0;
  int H_index; int He_index; // used to update N_species_by_elem
 
  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    
    // ================================================================
    //          (i) Identify elements present in tracer list
    //          (ii) Record X_mass_frac_index vector, N_elem,
    //               X_elem_atomic_mass, X_elem_number_density
    // ================================================================
    if (s[0]=='X'){
      X_mass_frac_index.push_back(ftr + N_elem); ///<record primitive vector index of each element
      N_species_by_elem.push_back(0);
      
      //=======Helium========
      if (s.substr(2,2)=="He"){
        He_index = N_elem;
        X_elem_atomic_mass.push_back(m_He);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
        He_ion_index.push_back(0); He_ion_index.push_back(0); //initialise the length of H_ion_index
      }
      //=======Hydrogen======
      else if (s[2]=='H'){
        H_index = N_elem;
        X_elem_atomic_mass.push_back(m_H);
        X_elem_number_density.push_back(0); //just to initialise the length of X_elem_number_density
        H_ion_index.push_back(0); //initialise the length of H_ion_index
      }
      N_elem++;
    } 
  }

  //
  // ================================================================
  // (iii) Record:   N_species,    N_species_by_elem,  y_ion_num_elec
  //                 y_ion_index_prim,  y_ip1_index,        y_im1_index
  // ================================================================
  //
  
  lv_y_ion_index_offset = ftr + N_elem; // gives the index at which ions first occur in primitive vector, maps to first index of local vector

  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    //=======Helium========
    if (s.substr(0,2)=="He"){
      cout << "\n\nTesting " << s << "\n";
      species_tracer_initialise(tracers, i, s, "He", 2, He_index, len);
      if (s[2]=='1'){
        cout << s <<"\n";
        y_ip1_index_tables.push_back(4); //index of He2+ in tables
        y_ion_index_tables.push_back(3); //index of He1+ in tables
        y_im1_index_tables.push_back(2); //index of He0 in tables
      }
      else if (s[2]=='2'){
        y_ip1_index_tables.push_back(-1); //doesn't exist in tables
        y_ion_index_tables.push_back(4); //index of He2+ in tables
        y_im1_index_tables.push_back(3); //index of He1+ in tables
      }
    }
    //=======Hydrogen========
    else if (s[0] =='H'){
      cout << "\n\nTesting " << s << "\n";
      species_tracer_initialise(tracers, i, s, "H", 1, H_index, len);
      y_ip1_index_tables.push_back(-1); //doesn't exist in tables
      y_ion_index_tables.push_back(1); //index of H1+ in tables
      y_im1_index_tables.push_back(0); //index of H0 in tables
    }
  }

  cout << "\nAfter reading in tracers, N_species=" << N_species << ", N_elements=" << N_elem << "\n\n";  
  
  
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
  /*mpv_nH   = 1.0;  // Hydrogen number density (density of H+ and H) -- change that where y_now is, and
  mpv_Vshell  = 1.0e54;
  mpv_Tau0     = 1.0e6;  
  mpv_dTau0    = 1.0;
  mpv_G0_UV   = 0.0;
  mpv_G0_IR   = 0.0;
  mpv_delta_S = 0.0;
  mpv_NIdot   = 0.0;
  ion_src_type = -1;
  N_ion_srcs  = 0;
  N_diff_srcs = 0;*/
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
  // ================================================================
  // ================================================================

  cout <<"MPv10: Constructor finished and returning.\n";
  cout <<"---------------------------------------------------------------------\n\n";
  return;
}


// ##################################################################
// ##################################################################

void MPv10::species_tracer_initialise(
    const std::string *tracers,  ///< List of what the tracer variables mean.
    int i, ///<index of current tracer in for loop
    string s, /// < current tracer in for loop
    string el_symbol, ///< element symbol, e.g. "He", "H"
    int el_symbol_length, ///< e.g. "H" is of length 1, "He" of length 2.
    int el_index, /// element index in N_species_by_elem, used in for loops for densities etc
    int length /// < length of tracers vector
    ){
  N_species_by_elem[el_index]++;
  y_ion_index_prim.push_back(lv_y_ion_index_offset + N_species);
  y_ion_index_local.push_back(N_species);
     
  /// Use stringstream to convert electron number string to int.
  int num_elec_int;
  stringstream ss_e; 
  ss_e << s.substr(el_symbol_length,1); 
  ss_e >> num_elec_int; 
  /// Record the electron number and He_ion_index vector (used to identify tracer index from string)
  y_ion_num_elec.push_back(num_elec_int);
  He_ion_index[num_elec_int-1] = lv_y_ion_index_offset + N_species;
        
  /// Define the tracer for the next ion up / down in tracers list
  stringstream ss_ip1;        
  stringstream ss_im1;
  stringstream ss_neutral;
  ss_ip1 << el_symbol << num_elec_int +1 << "+";  
  ss_im1 << el_symbol << num_elec_int -1 << "+";
  ss_neutral << el_symbol << 0 << "+";
  string neutral = ss_neutral.str();
  string ip1 = ss_ip1.str();
  string im1 = ss_im1.str();
      
  /// Check if the next ion up / down exists in tracer list (or if it's neutral) to define y_ip1_index and y_im1_index.
  ///< due to ordering in tracer list, if the next ion exists, its index will be 1 above the current index.
  if (i+1 < length){
    if ( ip1==tracers[i+1] ){
      cout << "Next ion up = " << ip1<<"\n";
      y_ip1_index_local.push_back(N_species + 1); 
    }
    else{
      cout <<"Next ion up doesn't exist\n";
      y_ip1_index_local.push_back(-1); // -1 flag => species doesn't exist 
    }
  }
  else{
    cout <<"Next ion up doesn't exist\n";
    y_ip1_index_local.push_back(-1); // -1 flag => species doesn't exist 
  }
  ///< due to ordering in tracer list, if the previous ion exists, its index will be 1 below the current index.
  if ( im1 == neutral){
    cout <<"Next ion down is neutral \n";
    y_im1_index_local.push_back(-2); // -2 flag => species is neutral
  }
  else if (im1==tracers[i-1] ){//check if the next lowest ion is the neutral species
    cout << "Next ion down = " << im1<<"\n";
    y_im1_index_local.push_back(lv_y_ion_index_offset + N_species - 1); 
  }
  else{
    y_im1_index_local.push_back(-1); // -1 flag => species doesn't exist
    cout <<"Next ion down doesn't exist.\n";
  }
  N_species++;
  
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

//NOTE This is currently incorrect, but I'm not sure if it's ever called anywhere, so I'll leave it for now...
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
      double *y_ion_frac, //< y_ion_fraction (by y_ion_index_prim)
      vector<pion_flt>& X_number_density,//const double nH, ///< nH
      const double E ///< E_int (per unit volume)
      )
{
  //
  // returns gas temperature according to E=nkT/(g-1), => T = E*(g-1)/(K*n)
  return gamma_minus_one*E/(k_B*get_ntot(y_ion_frac,X_number_density));
}

// ##################################################################
// ##################################################################

double MPv10::get_ntot(
      double *y_ion_frac,//<y_ion_fraction (by y_ion_index_prim)
      vector<pion_flt>& X_number_density//const double nH, ///< nH
      )
{
  int species_counter = 0;
  pion_flt n_tot=0;
  
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    pion_flt neutral_frac = 1; //neutral frac, got by subtracting off the ion fractions in the next loop.
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      pion_flt number_density = X_number_density[e]*y_ion_frac[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      
      n_tot += (1+num_elec)*number_density; //add on number of particles got by electrons + ionised atom
      
      neutral_frac -= y_ion_frac[species_counter];
      species_counter ++;
    }
    n_tot += neutral_frac*X_number_density[e];
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
  // ==============================================================
  //                 Set NUMBER DENSITY VECTOR.
  // ==============================================================
  //
  for (int i=0;i<N_elem;i++){//loop over every element
    double n_X = p_in[RO]*( p_in[ X_mass_frac_index[i]] / X_elem_atomic_mass[i]);
    X_elem_number_density[i] = n_X;
  }
  //
  // ==============================================================
  //                 Set INTERNAL ENERGY in local vector.
  // ==============================================================
  //
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);
  //
  // ==============================================================
  //            Set ION FRACTIONS OF ELEMENT in LOCAL vector,
  //            converted from MASS FRACTIONS in PRIM vector.
  // ==============================================================
  //  
  //loop over every species and set the local vector, noting that loc_vec[i]=prim_vec[i]/X_e
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      p_local[ y_ion_index_local[species_counter]] = p_in[ y_ion_index_prim[species_counter]]/p_in[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  
  //for (int v=0;v<nvl;++v) cout << "p_local[ " << v << "] = " << p_local[v]<<"\n";
  //for (int v=0;v<N_elem+N_species;++v) cout << "p_prim[ " << v+lv_y_ion_index_offset-N_elem << "] = " << p_in[v+lv_y_ion_index_offset-N_elem]<<"\n";
  cout << JM_MINNEU <<"\n\n";

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
    //Define y_ion_frac
    pion_flt y_ion_frac[N_species];
    for (int s=0;s<N_species;s++){
      y_ion_frac[ s] = p_local[ s];
    }
    //reset the internal energy (requires using y_ion_frac in get_ntot)
    p_local[lv_eint] = get_ntot(y_ion_frac,X_elem_number_density)*k_B*EP->MinTemperature/(gamma_minus_one);
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
  //
  //=================================================================
  //     Update output primitive vector according to local vector
  //     Also, define y_ion_frac while you're at it.
  //=================================================================
  //
  p_out[PG]    = p_local[lv_eint]*(gamma_minus_one);
  pion_flt y_ion_frac[N_species];

  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      p_out[ y_ion_index_prim[ species_counter]] = p_local[ y_ion_index_local[species_counter]] * p_in[ X_mass_frac_index[e]];
      y_ion_frac[ species_counter] = p_local[ y_ion_index_local[species_counter]];
      species_counter ++;
    }
  }

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

  // NOTE \Maggie{Must come up with alternative to this for every element
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  //p_out[pv_H1p] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[pv_H1p])));

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected x(H+) from p_out[]).
  //
  
  double T = get_temperature(y_ion_frac, X_elem_number_density, p_local[lv_eint]);
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
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      y_ion_frac[ species_counter] = pv[ y_ion_index_prim[ species_counter]] / pv[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  return (get_temperature(y_ion_frac, X_elem_number_density, P[lv_eint]));
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
  
  //Determine y_ion_frac from the primitive vector
  pion_flt y_ion_frac[N_species];
  int species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      y_ion_frac[ species_counter] = p_pv[ y_ion_index_prim[ species_counter]] / p_pv[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  
  //Determine internal energy using get_ntot
  P[lv_eint] = get_ntot(y_ion_frac,X_elem_number_density)*k_B*T/(gamma_minus_one);
  
  //Call convert_local2prim with the new local vector; this will generate a new temperature value;
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

  //
  // Calculate y-dot[] to see if anything is changing significantly over dt
  //
  double maxdelta=0.0;
  err = ydot(0, y_in, y_out, 0);
  cout <<"H+-ion-frac="<< NV_Ith_S(y_in,0) <<"\n";
  cout <<"He+-ion-frac="<< NV_Ith_S(y_in,1) <<"\n";
  cout <<"He2+-ion-frac="<< NV_Ith_S(y_in,2) <<"\n";
  cout <<"E_int="<< NV_Ith_S(y_in,3) <<"\n";

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
  for (int v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];
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
  double OneMinusX = max(NV_Ith_S(y_now,lv_H0),Min_NeutralFrac); //NOTE \Maggie{LEGACY CODE; REMOVE}
  double x_in      = 1.0-OneMinusX; //NOTE \Maggie{LEGACY CODE; REMOVE}

  double ne=0;
  double y_ion_frac[N_species];
  
  int species_counter=0;
  //
  //  ========================================================
  //        Determine ne, y_ion_frac, and number density
  //  ========================================================
  //
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //add to y_ion_frac
      //int this_ion_index = y_ion_index[species_counter] - lv_y_ion_index_offset;
      y_ion_frac[species_counter] = NV_Ith_S(y_now,y_ion_index_local[species_counter]);
      
      //add to ne based on the number of electrons liberated to obtain this ion
      pion_flt number_density = X_elem_number_density[elem]*y_ion_frac[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      ne += num_elec*number_density;
      
      //cout <<"y-ion-frac="<<y_ion_frac[species_counter] << ", y-ion number density="<<number_density <<"\n";
      species_counter ++;
    }
  }
  // Find internal energy
  double E_in      = NV_Ith_S(y_now,lv_eint);

  // Use E_in, y_ion_frac and number density to determine temperature
  double T = get_temperature(y_ion_frac, X_elem_number_density, E_in);
  cout << "Temperature="<< T<<"\n"<<"Electron density="<<ne<<"\n\n\n";

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
  // ne += mpv_nH*1.5e-4*METALLICITY*exp(-mpv_nH/1.0e4);


  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0
  // NOTE \Maggie {again, lookup table instead of calculation}
  
  //
  //  ========================================================
  //          Get Rate of Change of Each Species
  //  ========================================================
  //
  /*species_counter = 0;
  pion_flt coll_ion_out;
  pion_flt coll_ion_in;
  pion_flt recomb_out;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      ///  =============Collisional Ionisation Terms=================
      if (y_ip1_index[species_counter] != -1){ //<<< no species more ionised
         //coll_ion_out = ;
        
      
      species_counter ++;
    }
  }*/

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

  // =========================================================================
  //
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



