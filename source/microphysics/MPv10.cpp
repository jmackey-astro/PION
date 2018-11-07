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
  atol[0] = MPv10_ABSTOL; ///< minimum neutral fraction I care about.
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
  T_min(1e0), T_max(1e9), Num_temps(1e2)
  {
  int max_num_species = 5; //number of species recorded in each table, i.e. the number of possible species recognised by the module
  
  double erg_per_eV = 1.60218e-12;
  /// ===================================================================
  /// Initialise ionisation potentials vector
  /// ===================================================================
  double ionisation_pot_arr[max_num_species] = {13.59844, -1.0e99,
                                                24.58741, 54.41778, -1.0e99}; //energy (eV) required to raise ion from species i to species i+1
  for (int i=0; i<max_num_species; i++) ionisation_pot_arr[i]*=erg_per_eV; //convert eV to erg
  ionisation_potentials.insert(ionisation_potentials.end(), &ionisation_pot_arr[0], &ionisation_pot_arr[5]);

  /// ===================================================================
  ///  Initialise Temperature, Recombination (radiative + dielectronic),
  ///       Ionisation, and (recomb and ionisation) Slopes Tables
  ///  NOTE: SHOULD SET THIS TO READ FROM FILE EVENTUALLY;
  /// ===================================================================
  delta_log_temp=(log10f(T_max) - log10f(T_min))/(Num_temps-1);
  //int Num_temps = 1e2;
  // NOTE 2D arrays were acting up with initialising in the header file, so instead I've flattened into 1D arrays, where
  //  1D_array[ x + width*y] = 2D_array[x,y]
  double Temp_arr[Num_temps] = { 1, 1.23285, 1.51991, 1.87382, 2.31013, 2.84804, 3.51119, 4.32876, 5.3367, 6.57933, 8.11131, 10, 12.3285, 15.1991, 18.7382, 23.1013, 28.4804, 35.1119, 43.2876, 53.367, 65.7933, 81.1131, 100, 123.285, 151.991, 187.382, 231.013, 284.804, 351.119, 432.876, 533.67, 657.933, 811.131, 1000, 1232.85, 1519.91, 1873.82, 2310.13, 2848.04, 3511.19, 4328.76, 5336.7, 6579.33, 8111.31, 10000, 12328.5, 15199.1, 18738.2, 23101.3, 28480.4, 35111.9, 43287.6, 53367, 65793.4, 81113.1, 100000, 123285, 151991, 187382, 231013, 284804, 351119, 432876, 533670, 657934, 811131, 1e+06, 1.23285e+06, 1.51991e+06, 1.87382e+06, 2.31013e+06, 2.84804e+06, 3.5112e+06, 4.32876e+06, 5.3367e+06, 6.57934e+06, 8.11131e+06, 1e+07, 1.23285e+07, 1.51991e+07, 1.87382e+07, 2.31013e+07, 2.84804e+07, 3.5112e+07, 4.32876e+07, 5.3367e+07, 6.57933e+07, 8.11131e+07, 1e+08, 1.23285e+08, 1.51991e+08, 1.87382e+08, 2.31013e+08, 2.84804e+08, 3.5112e+08, 4.32876e+08, 5.3367e+08, 6.57934e+08, 8.11131e+08, 1e+09};
  
  Temp_Table.insert(Temp_Table.end(), &Temp_arr[0], &Temp_arr[Num_temps]);

  
  //***********************************************************
  
  double recomb_arr[max_num_species][Num_temps] = {{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{3.41202e-10, 2.89621e-10, 2.45838e-10, 2.08674e-10, 1.77128e-10, 1.50351e-10, 1.27622e-10, 1.08329e-10, 9.19525e-11, 7.80517e-11, 6.62524e-11, 5.62368e-11, 4.77353e-11, 4.0519e-11, 3.43936e-11, 2.91942e-11, 2.47808e-11, 2.10346e-11, 1.78547e-11, 1.51556e-11, 1.28645e-11, 1.09197e-11, 9.26893e-12, 7.86771e-12, 6.67832e-12, 5.66874e-12, 4.81178e-12, 4.08436e-12, 3.46692e-12, 2.94281e-12, 2.49794e-12, 2.12032e-12, 1.79978e-12, 1.5277e-12, 1.29675e-12, 1.10072e-12, 9.34319e-13, 7.93075e-13, 6.73184e-13, 5.71416e-13, 4.85033e-13, 4.11709e-13, 3.4947e-13, 2.96639e-13, 2.51795e-13, 2.13731e-13, 1.8142e-13, 1.53994e-13, 1.30714e-13, 1.10954e-13, 9.41806e-14, 7.9943e-14, 6.78577e-14, 5.75995e-14, 4.8892e-14, 4.15008e-14, 3.5227e-14, 2.99016e-14, 2.53813e-14, 2.15443e-14, 1.82874e-14, 1.55228e-14, 1.31762e-14, 1.11843e-14, 9.49352e-15, 8.05836e-15, 6.84015e-15, 5.8061e-15, 4.92837e-15, 4.18333e-15, 3.55092e-15, 3.01412e-15, 2.55846e-15, 2.17169e-15, 1.84339e-15, 1.56472e-15, 1.32818e-15, 1.12739e-15, 9.56959e-16, 8.12292e-16, 6.89495e-16, 5.85262e-16, 4.96786e-16, 4.21685e-16, 3.57938e-16, 3.03827e-16, 2.57896e-16, 2.18909e-16, 1.85816e-16, 1.57726e-16, 1.33882e-16, 1.13642e-16, 9.64626e-17, 8.18801e-17, 6.9502e-17, 5.89951e-17, 5.00767e-17, 4.25064e-17, 3.60806e-17, 3.06262e-17},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1.33093e-10, 1.17666e-10, 1.03994e-10, 9.18837e-11, 8.11612e-11, 7.16715e-11, 6.32763e-11, 5.58519e-11, 4.92884e-11, 4.34877e-11, 3.83627e-11, 3.38359e-11, 2.98384e-11, 2.63091e-11, 2.3194e-11, 2.04448e-11, 1.8019e-11, 1.5879e-11, 1.39914e-11, 1.23265e-11, 1.08584e-11, 9.56392e-12, 8.42264e-12, 7.41655e-12, 6.52971e-12, 5.74807e-12, 5.0592e-12, 4.45214e-12, 3.91723e-12, 3.44593e-12, 3.0307e-12, 2.6649e-12, 2.34268e-12, 2.05887e-12, 1.80891e-12, 1.58879e-12, 1.39495e-12, 1.2243e-12, 1.07406e-12, 9.41811e-13, 8.25422e-13, 7.23005e-13, 6.32898e-13, 5.5364e-13, 4.8394e-13, 4.22663e-13, 3.68807e-13, 3.2149e-13, 2.79936e-13, 2.43486e-13, 2.11883e-13, 1.87118e-13, 1.791e-13, 2.12928e-13, 3.25733e-13, 5.44478e-13, 8.59022e-13, 1.21617e-12, 1.54153e-12, 1.77147e-12, 1.87401e-12, 1.85128e-12, 1.72897e-12, 1.54243e-12, 1.32565e-12, 1.10534e-12, 8.9925e-13, 7.17141e-13, 5.62742e-13, 4.35843e-13, 3.34002e-13, 2.53772e-13, 1.91481e-13, 1.43671e-13, 1.0731e-13, 7.98579e-14, 5.92515e-14, 4.38562e-14, 3.23971e-14, 2.38937e-14, 1.7599e-14, 1.29486e-14, 9.51838e-15, 6.9916e-15, 5.13234e-15, 3.76548e-15, 2.76137e-15, 2.02421e-15, 1.48333e-15, 1.08664e-15, 7.95822e-16, 5.82695e-16, 4.26555e-16, 3.12194e-16, 2.28454e-16, 1.67149e-16, 1.22278e-16, 8.94407e-17, 6.54141e-17, 4.78367e-17},
{5.39147e-10, 4.82303e-10, 4.31218e-10, 3.85324e-10, 3.4411e-10, 3.07114e-10, 2.73921e-10, 2.44153e-10, 2.17473e-10, 1.93575e-10, 1.72182e-10, 1.53045e-10, 1.35938e-10, 1.20656e-10, 1.07016e-10, 9.48513e-11, 8.40103e-11, 7.43571e-11, 6.57686e-11, 5.81336e-11, 5.13516e-11, 4.53321e-11, 3.99936e-11, 3.52626e-11, 3.10731e-11, 2.73656e-11, 2.4087e-11, 2.11895e-11, 1.86304e-11, 1.63716e-11, 1.43789e-11, 1.2622e-11, 1.10736e-11, 9.70986e-12, 8.50918e-12, 7.45259e-12, 6.5232e-12, 5.70605e-12, 4.98789e-12, 4.35698e-12, 3.80298e-12, 3.3167e-12, 2.89006e-12, 2.51593e-12, 2.18799e-12, 1.9007e-12, 1.64915e-12, 1.42904e-12, 1.23657e-12, 1.06839e-12, 9.21565e-13, 7.93493e-13, 6.81901e-13, 5.84779e-13, 5.00362e-13, 4.27095e-13, 3.63609e-13, 3.087e-13, 2.61305e-13, 2.20488e-13, 1.85424e-13, 1.55385e-13, 1.29728e-13, 1.07884e-13, 8.93538e-14, 7.36932e-14, 6.05118e-14, 4.94646e-14, 4.02482e-14, 3.25956e-14, 2.62727e-14, 2.10751e-14, 1.68247e-14, 1.33675e-14, 1.05704e-14, 8.31963e-15, 6.51816e-15, 5.08399e-15, 3.9482e-15, 3.05334e-15, 2.3518e-15, 1.80449e-15, 1.37947e-15, 1.0509e-15, 7.97955e-16, 6.04022e-16, 4.55896e-16, 3.43163e-16, 2.57653e-16, 1.92997e-16, 1.44251e-16, 1.076e-16, 8.01122e-17, 5.95447e-17, 4.41881e-17, 3.27448e-17, 2.4233e-17, 1.79121e-17, 1.32255e-17, 9.75524e-18}};
  recomb_rate_table.resize( max_num_species);
  for (int i = 0; i < max_num_species; ++i) recomb_rate_table[i].resize(Num_temps);
  for (int row=0;row < max_num_species; row++) {
    for (int col=0; col<Num_temps; col++){
      recomb_rate_table[row][col] = recomb_arr[row][col];
    }
  }
  
  
  double ionise_arr[max_num_species][Num_temps] = {{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.28882e-22, 1.59626e-19, 1.67741e-17, 7.48925e-16, 1.66962e-14, 2.11854e-13, 1.70127e-12, 9.42312e-12, 3.85978e-11, 1.23704e-10, 3.24709e-10, 7.24291e-10, 1.41442e-09, 2.4769e-09, 3.96526e-09, 5.89348e-09, 8.23294e-09, 1.09161e-08, 1.38445e-08, 1.68983e-08, 1.99467e-08, 2.28586e-08, 2.55124e-08, 2.78051e-08, 2.966e-08, 3.10306e-08, 3.19021e-08, 3.22886e-08, 3.22282e-08, 3.1776e-08, 3.09966e-08, 2.99578e-08, 2.8725e-08, 2.73575e-08, 2.59071e-08, 2.44166e-08, 2.29202e-08, 2.14443e-08, 2.00086e-08, 1.8627e-08, 1.73085e-08, 1.6059e-08, 1.48811e-08, 1.37756e-08, 1.27416e-08, 1.17773e-08, 1.088e-08, 1.00465e-08, 9.27355e-09, 8.55749e-09, 7.89482e-09, 7.28205e-09, 6.71577e-09, 6.19273e-09, 5.70983e-09, 5.26414e-09, 4.8529e-09, 4.47354e-09, 4.12365e-09},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.98674e-33, 2.21776e-28, 9.14578e-25, 8.02956e-22, 2.00991e-19, 1.81788e-17, 7.19909e-16, 1.45878e-14, 1.71622e-13, 1.29871e-12, 6.86858e-12, 2.71553e-11, 8.47503e-11, 2.18214e-10, 4.80361e-10, 9.30467e-10, 1.62316e-09, 2.59807e-09, 3.87311e-09, 5.44211e-09, 7.27585e-09, 9.32517e-09, 1.15249e-08, 1.37981e-08, 1.60613e-08, 1.82293e-08, 2.02211e-08, 2.19656e-08, 2.34072e-08, 2.45094e-08, 2.52564e-08, 2.56527e-08, 2.57199e-08, 2.54927e-08, 2.5014e-08, 2.43301e-08, 2.34871e-08, 2.25279e-08, 2.14905e-08, 2.04072e-08, 1.93045e-08, 1.82034e-08, 1.71199e-08, 1.60659e-08, 1.50499e-08, 1.40776e-08, 1.31522e-08, 1.22756e-08, 1.14483e-08, 1.06698e-08, 9.9389e-09, 9.25409e-09, 8.61343e-09, 8.01483e-09, 7.45609e-09, 6.935e-09, 6.44934e-09, 5.99694e-09, 5.57571e-09},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.09045e-37, 4.97378e-32, 8.42609e-28, 2.29911e-24, 1.42686e-21, 2.66468e-19, 1.87911e-17, 6.01558e-16, 1.01498e-14, 1.01881e-13, 6.71149e-13, 3.14187e-12, 1.1148e-11, 3.1587e-11, 7.45508e-11, 1.51648e-10, 2.73274e-10, 4.4601e-10, 6.71085e-10, 9.44173e-10, 1.25623e-09, 1.59485e-09, 1.94574e-09, 2.29406e-09, 2.62559e-09, 2.92773e-09, 3.19025e-09, 3.40582e-09, 3.57033e-09, 3.68281e-09, 3.74511e-09, 3.76131e-09, 3.73704e-09, 3.67877e-09, 3.59317e-09, 3.48662e-09, 3.36494e-09, 3.23311e-09, 3.09533e-09, 2.95493e-09, 2.81453e-09, 2.67609e-09, 2.54102e-09, 2.41033e-09, 2.28466e-09, 2.1644e-09, 2.04973e-09, 1.94069e-09, 1.83723e-09, 1.7392e-09, 1.64641e-09, 1.55867e-09, 1.47572e-09, 1.39734e-09, 1.32328e-09, 1.2533e-09},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  
  ionise_rate_table.resize( max_num_species);
  for (int i = 0; i < max_num_species; ++i) ionise_rate_table[i].resize(Num_temps);
  for (int row=0;row < max_num_species; row++) {
    for (int col=0; col<Num_temps; col++){
      ionise_rate_table[row][col] = ionise_arr[row][col];
   }
  }
  
  double recomb_slope_arr[max_num_species][Num_temps] = {{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-2.21522e-10, -1.5252e-10, -1.05011e-10, -7.23013e-11, -4.97802e-11, -3.42741e-11, -2.35981e-11, -1.62475e-11, -1.11865e-11, -7.70204e-12, -5.30293e-12, -3.65112e-12, -2.51383e-12, -1.73079e-12, -1.19167e-12, -8.20475e-13, -5.64905e-13, -3.88942e-13, -2.6779e-13, -1.84376e-13, -1.26945e-13, -8.74027e-14, -6.01776e-14, -4.14328e-14, -2.85269e-14, -1.9641e-14, -1.3523e-14, -9.31074e-15, -6.41053e-15, -4.41371e-15, -3.03888e-15, -2.0923e-15, -1.44057e-15, -9.91844e-16, -6.82894e-16, -4.70179e-16, -3.23723e-16, -2.22886e-16, -1.53459e-16, -1.05658e-16, -7.27467e-17, -5.00867e-17, -3.44852e-17, -2.37434e-17, -1.63476e-17, -1.12554e-17, -7.74947e-18, -5.33558e-18, -3.6736e-18, -2.52931e-18, -1.74145e-18, -1.19901e-18, -8.25527e-19, -5.68384e-19, -3.91338e-19, -2.6944e-19, -1.85512e-19, -1.27726e-19, -8.7941e-20, -6.05482e-20, -4.1688e-20, -2.87026e-20, -1.9762e-20, -1.36063e-20, -9.36808e-21, -6.45002e-21, -4.44089e-21, -3.0576e-21, -2.10518e-21, -1.44944e-21, -9.97953e-22, -6.871e-22, -4.73074e-22, -3.25716e-22, -2.24259e-22, -1.54404e-22, -1.06309e-22, -7.31946e-23, -5.03952e-23, -3.46976e-23, -2.38896e-23, -1.64482e-23, -1.13248e-23, -7.7972e-24, -5.36844e-24, -3.69623e-24, -2.54489e-24, -1.75218e-24, -1.20639e-24, -8.30613e-25, -5.71884e-25, -3.93747e-25, -2.71099e-25, -1.86654e-25, -1.28513e-25, -8.84826e-26, -6.09211e-26, -4.19447e-26, -2.88794e-26, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{-6.62552e-11, -4.76266e-11, -3.42193e-11, -2.45753e-11, -1.76418e-11, -1.26596e-11, -9.08102e-12, -6.51184e-12, -4.66806e-12, -3.34536e-12, -2.39681e-12, -1.71679e-12, -1.22942e-12, -8.80227e-13, -6.30093e-13, -4.5096e-13, -3.22701e-13, -2.30886e-13, -1.65172e-13, -1.18146e-13, -8.44984e-14, -6.0427e-14, -4.32085e-14, -3.08933e-14, -2.20862e-14, -1.57884e-14, -1.12855e-14, -8.06615e-15, -5.76471e-15, -4.11959e-15, -2.94371e-15, -2.1033e-15, -1.50269e-15, -1.07349e-15, -7.66812e-16, -5.47689e-16, -3.91139e-16, -2.79304e-16, -1.99418e-16, -1.4236e-16, -1.01611e-16, -7.25124e-17, -5.1736e-17, -3.69037e-17, -2.63166e-17, -1.87609e-17, -1.33698e-17, -9.52385e-18, -6.77637e-18, -4.76551e-18, -3.02904e-18, -7.95572e-19, 2.72232e-18, 7.3634e-18, 1.15818e-17, 1.35086e-17, 1.24414e-17, 9.19337e-18, 5.26998e-18, 1.9064e-18, -3.42799e-19, -1.49599e-18, -1.85071e-18, -1.74451e-18, -1.43811e-18, -1.09117e-18, -7.82098e-19, -5.37852e-19, -3.58567e-19, -2.33414e-19, -1.49152e-19, -9.39312e-20, -5.84776e-20, -3.60743e-20, -2.20922e-20, -1.34508e-20, -8.15135e-21, -4.92128e-21, -2.96219e-21, -1.77863e-21, -1.06586e-21, -6.37693e-22, -3.81022e-22, -2.27413e-22, -1.3561e-22, -8.08049e-23, -4.81181e-23, -2.86382e-23, -1.70365e-23, -1.01307e-23, -6.02207e-24, -3.57863e-24, -2.12604e-24, -1.26275e-24, -7.49841e-25, -4.45179e-25, -2.64255e-25, -1.56834e-25, -9.30667e-26, 0},
{-2.44128e-10, -1.77957e-10, -1.29677e-10, -9.44596e-11, -6.87773e-11, -5.00543e-11, -3.64096e-11, -2.64697e-11, -1.92319e-11, -1.39643e-11, -1.01326e-11, -7.34706e-12, -5.32334e-12, -3.85409e-12, -2.78815e-12, -2.0154e-12, -1.45564e-12, -1.05049e-12, -7.57493e-13, -5.45775e-13, -3.92921e-13, -2.82656e-13, -2.03181e-13, -1.45945e-13, -1.04758e-13, -7.51438e-14, -5.3866e-14, -3.85892e-14, -2.76286e-14, -1.977e-14, -1.4139e-14, -1.01067e-14, -7.22077e-15, -5.15649e-15, -3.68068e-15, -2.62609e-15, -1.87286e-15, -1.33511e-15, -9.51362e-16, -6.77627e-16, -4.82447e-16, -3.43332e-16, -2.44219e-16, -1.73632e-16, -1.23383e-16, -8.76266e-17, -6.2194e-17, -4.41131e-17, -3.12652e-17, -2.21409e-17, -1.56648e-17, -1.10714e-17, -7.81582e-18, -5.51035e-18, -3.87925e-18, -2.7265e-18, -1.91279e-18, -1.3392e-18, -9.35492e-19, -6.51857e-19, -4.52973e-19, -3.13825e-19, -2.16712e-19, -1.49123e-19, -1.02225e-19, -6.97915e-20, -4.74437e-20, -3.21057e-20, -2.16233e-20, -1.44916e-20, -9.66267e-21, -6.40925e-21, -4.22869e-21, -2.77505e-21, -1.8113e-21, -1.17591e-21, -7.59347e-22, -4.87782e-22, -3.1173e-22, -1.98226e-22, -1.25441e-22, -7.90129e-23, -4.95472e-23, -3.09381e-23, -1.92406e-23, -1.19203e-23, -7.35869e-24, -4.52745e-24, -2.77678e-24, -1.69809e-24, -1.03561e-24, -6.29999e-25, -3.82362e-25, -2.31568e-25, -1.39968e-25, -8.4448e-26, -5.08665e-26, -3.05925e-26, -1.83736e-26, 0}};
  
  recomb_slope_table.resize( max_num_species);
  for (int i = 0; i < max_num_species; ++i) recomb_slope_table[i].resize(Num_temps);
  for (int row=0;row < max_num_species; row++) {
    for (int col=0; col<Num_temps; col++){
      recomb_slope_table[row][col] = recomb_slope_arr[row][col];
   }
  }
  
  
  double ionise_slope_arr[max_num_species][Num_temps] = {{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.24717e-25, 1.28032e-22, 1.08451e-20, 3.8765e-19, 6.84884e-18, 6.79836e-17, 4.20852e-16, 1.7698e-15, 5.42375e-15, 1.28336e-14, 2.45855e-14, 3.96436e-14, 5.55377e-14, 6.93535e-14, 7.88041e-14, 8.28107e-14, 8.14957e-14, 7.58164e-14, 6.71169e-14, 5.67714e-14, 4.59684e-14, 3.56163e-14, 2.63286e-14, 1.84509e-14, 1.21079e-14, 7.25695e-15, 3.74255e-15, 1.34631e-15, -1.70616e-16, -1.03641e-15, -1.44883e-15, -1.56645e-15, -1.50797e-15, -1.35665e-15, -1.1672e-15, -9.72954e-16, -7.92298e-16, -6.33829e-16, -5.00134e-16, -3.9041e-16, -3.02173e-16, -2.32299e-16, -1.77617e-16, -1.35218e-16, -1.02581e-16, -7.76026e-17, -5.85731e-17, -4.41291e-17, -3.31975e-17, -2.49441e-17, -1.87243e-17, -1.40444e-17, -1.05274e-17, -7.88713e-18, -5.90655e-18, -4.42183e-18, -3.30941e-18, -2.47629e-18, -1.85257e-18, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.92385e-36, 1.78466e-31, 5.96848e-28, 4.24655e-25, 8.59742e-23, 6.26262e-21, 1.98281e-19, 3.17844e-18, 2.91936e-17, 1.69958e-16, 6.81269e-16, 2.01269e-15, 4.63492e-15, 8.7119e-15, 1.38798e-14, 1.93306e-14, 2.41302e-14, 2.75471e-14, 2.92231e-14, 2.91686e-14, 2.76518e-14, 2.50659e-14, 2.18236e-14, 1.8294e-14, 1.47732e-14, 1.14787e-14, 8.55379e-15, 6.07704e-15, 4.07346e-15, 2.52618e-15, 1.38882e-15, 5.97583e-16, 8.21965e-17, -2.25401e-16, -3.85252e-16, -4.46407e-16, -4.46342e-16, -4.11952e-16, -3.61388e-16, -3.06098e-16, -2.52734e-16, -2.04707e-16, -1.63382e-16, -1.28911e-16, -1.00798e-16, -7.82519e-17, -6.04026e-17, -4.64119e-17, -3.55309e-17, -2.71205e-17, -2.06516e-17, -1.56955e-17, -1.19103e-17, -9.02653e-18, -6.83407e-18, -5.16988e-18, -3.90834e-18, -2.95304e-18, -2.23026e-18, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.6363e-40, 2.13606e-35, 2.93508e-31, 6.49401e-28, 3.26499e-25, 4.92727e-23, 2.79341e-21, 7.12801e-20, 9.47308e-19, 7.38202e-18, 3.71591e-17, 1.30817e-16, 3.43836e-16, 7.11999e-16, 1.21399e-15, 1.76701e-15, 2.26111e-15, 2.60476e-15, 2.75296e-15, 2.70937e-15, 2.51126e-15, 2.21036e-15, 1.85783e-15, 1.4959e-15, 1.15491e-15, 8.53734e-16, 6.01663e-16, 4.00759e-16, 2.48071e-16, 1.3758e-16, 6.18082e-17, 1.30386e-17, -1.58399e-17, -3.0853e-17, -3.67645e-17, -3.71147e-17, -3.43841e-17, -3.02125e-17, -2.56149e-17, -2.11709e-17, -1.7173e-17, -1.37355e-17, -1.08691e-17, -8.53086e-18, -6.65381e-18, -5.16494e-18, -3.99461e-18, -3.08094e-18, -2.37136e-18, -1.82246e-18, -1.39909e-18, -1.07327e-18, -8.22921e-19, -6.30786e-19, -4.83443e-19, -3.7051e-19, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  
  ionise_slope_table.resize( max_num_species);
  for (int i = 0; i < max_num_species; ++i) ionise_slope_table[i].resize(Num_temps);
  for (int row=0;row < max_num_species; row++) {
    for (int col=0; col<Num_temps; col++){
      ionise_slope_table[row][col] = ionise_slope_arr[row][col];
   }
  }
  
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
      cout << "i= " << i << ", s=" << s <<", H_index = " << H_index <<"\n";
      species_tracer_initialise(tracers, i, s, "H", 1, H_index, len);
      y_ip1_index_tables.push_back(-1); //doesn't exist in tables
      y_ion_index_tables.push_back(1); //index of H1+ in tables
      y_im1_index_tables.push_back(0); //index of H0 in tables
    }
  }

  cout << "\nAfter reading in tracers, N_species=" << N_species;
  cout << ", N_elements=" << N_elem << "\n\n";  
  
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
  Min_NeutralFrac     = MPv10_ABSTOL;
  Max_NeutralFrac     = 1.0-MPv10_ABSTOL;

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
//   if (s
//   He_ion_index[num_elec_int-1] = lv_y_ion_index_offset + N_species;
        
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
    y_im1_index_local.push_back(N_species - 1); 
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
  N_equations = nvl;
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
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
//   if (s.substr(0,2)=="He"){
//     int num_elec_int; 
//     stringstream ss; ss << s.substr(2,1); ss >> num_elec_int; //Use stringstream to convert string to int.
//     int ion_index = He_ion_index[num_elec_int-1];
//     return ion_index;
//   }
//   else if (s.substr(0,1)=="H"){
//     int num_elec_int;
//     stringstream ss; ss << s.substr(1,1); ss >> num_elec_int; //Use stringstream to convert string to int.
//     int ion_index = H_ion_index[num_elec_int-1];
//     return ion_index;
//   }
//   else { return-1;}
    
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
  /*cout << "temp=" << gamma_minus_one*E/(k_B*get_ntot(y_ion_frac,X_number_density));
  cout << ", ntot = " <<get_ntot(y_ion_frac,X_number_density) ;
  cout << ", E = " << E << "\n";*/
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
  //rep.printSTLVec("n_x",X_elem_number_density);
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
      p_local[ y_ion_index_local[species_counter]] = 
        p_in[y_ion_index_prim[species_counter]]/p_in[ X_mass_frac_index[e]];
      species_counter ++;
    }
  }
  //rep.printVec("local",p_local,nvl);
  //rep.printVec("prim ",p_in,nv_prim);
  
  //for (int v=0;v<nvl;++v) cout << "p_local[ " << v << "] = " << p_local[v]<<"\n";
  //for (int v=0;v<N_elem+N_species;++v) cout << "p_prim[ " << v+lv_y_ion_index_offset-N_elem << "] = " << p_in[v+lv_y_ion_index_offset-N_elem]<<"\n";

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
#ifdef MPv10_DEBUG
      cout <<"e="<<e<<", s="<<s;
      cout <<", species_counter="<<species_counter;
      cout <<", ion index=";
      cout <<y_ion_index_local[species_counter]<<"\n";
      cout <<"Min_NeutralFrac="<<Min_NeutralFrac;
      cout <<", 1-Max_NeutralFrac="<<1.0-Max_NeutralFrac<<"\n";
      if (p_local[y_ion_index_local[species_counter]]>1.01 || 
          p_local[y_ion_index_local[species_counter]]<-0.01) {
        cout <<"MPv10::convert_prim2local: bad ion fraction: ";
        cout <<"x(H0)="<<p_local[y_ion_index_local[species_counter]];
        cout <<", resetting to [0,1]\n";
      }
#endif
      //p_local[y_ion_index_local[species_counter]] =
      //      max(Min_NeutralFrac, min(Max_NeutralFrac,
      //                               p_local[y_ion_index_local[species_counter]]));
      species_counter ++;
    }
  }

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

  //rep.printVec("prim2local local",p_local,nvl);
  //rep.printVec("prim2local prim ",p_in,nv_prim);

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

  // Set mass fraction tracers to be within the required range (not too close to zero or 1)
  
  species_counter=0;
  for (int e=0;e<N_elem;e++){//loop over every element
    int N_elem_species=N_species_by_elem[e];
    p_out[ X_mass_frac_index[ e]] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[ X_mass_frac_index[ e]])));
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      p_out[ y_ion_index_prim[ species_counter]] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[ y_ion_index_prim[ species_counter]])));
      species_counter ++;
    }
  }
  
  

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected x(H+) from p_out[]).
  //
  
  double T = get_temperature(y_ion_frac, X_elem_number_density, p_local[lv_eint]);
  //cout <<"current temp = " << T<<"\n";
  if (T>1.01*EP->MaxTemperature) {
    Set_Temp(p_out,EP->MaxTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() HIGH temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MaxTemperature<<"\n";
#endif // MPv10_DEBUG
  }
  if (T<0.99*EP->MinTemperature) {
    cout << "Low temperature!!!\n";
    Set_Temp(p_out,EP->MinTemperature,0);
#ifdef MPv10_DEBUG
    cout <<"MPv10::convert_local2prim() LOW  temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_H1p]<<"...  limiting to T="<<EP->MinTemperature<<"\n";
#endif // MPv10_DEBUG
  }
  p_out[PG] = p_local[lv_eint] *(gamma_minus_one);
  cout << "pressure=" <<p_out[PG]<<"\n";
  //rep.printVec("local2prim local",p_local,nvl);
  //rep.printVec("local2prim prim ",p_out,nv_prim);
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
  //cout << "set temperature\n";
  if (p_pv[PG]<=0.0) {
    cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
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
  rep.printVec("p2l start prim ",p_in,nv_prim);
  rep.printVec("p2l start local",P,nvl);
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
  /*cout <<"H+-ion-frac="<< NV_Ith_S(y_in,0) <<"\n";
  cout <<"He+-ion-frac="<< NV_Ith_S(y_in,1) <<"\n";
  cout <<"He2+-ion-frac="<< NV_Ith_S(y_in,2) <<"\n";
  cout <<"E_int="<< NV_Ith_S(y_in,3) <<"\n";*/

  if (err) 
    rep.error("dYdt() returned an error in MPv10::TimeUpdateMP_RTnew()",err);
  for (int v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out,v)*dt/NV_Ith_S(y_in,v)));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < EULER_CUTOFF) {
    cout <<"maxdelta="<<maxdelta<<", Doing euler integration...\n";
    for (int v=0;v<nvl;v++) {
      NV_Ith_S(y_out,v) = NV_Ith_S(y_in,v) + dt*NV_Ith_S(y_out,v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    cout <<"maxdelta="<<maxdelta<<", Doing CVODE integration...\n";
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
  //rep.printVec("l2p end prim ",p_out,nv_prim);
  //rep.printVec("l2p end local",P,nvl);

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
  //LOOP OVER ALL THE IONS HERE, NOT JUST LV_H0, AND THAT SHOULD FIX YOUR SEGFAULT PROBLEM.
  for (int v=0;v<N_equations;v++) t = min(t,DTFRAC/(fabs(NV_Ith_S(y_out, v))+TINYVALUE));
  //cout <<"limit by dx: dt="<<DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\n";


#ifdef MPv10_DEBUG
  if (t<3.16e9) {
  //cout <<"MP timescales: xdot="<<NV_Ith_S(y_out, lv_H0);
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
  //  ========================================================
  //        Determine ne, y_ion_frac, and number density
  //        Also determine neutral fraction for later
  //  ========================================================
  //
  
  double ne=0;
  double y_ion_frac[N_species];
  double X_neutral_frac[N_elem];
  
  int species_counter=0;
  for (int elem=0;elem<N_elem;elem++) {//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++) {//loop over every species in THIS element
      //add to y_ion_frac
      y_ion_frac[species_counter] = NV_Ith_S(y_now,y_ion_index_local[species_counter]);
      
      //add to ne based on the number of electrons liberated to obtain this ion
      pion_flt number_density = X_elem_number_density[elem]*y_ion_frac[species_counter];
      int num_elec = y_ion_num_elec[species_counter];
      ne += num_elec*number_density;
      species_counter ++;
    }
  }

  // Find internal energy
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double dydt[N_equations];
//   for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
//   rep.printVec("ydot going in",dydt,N_equations);
    //initialise ydot to zero.
  for (int v=0;v<N_equations;v++) NV_Ith_S(y_dot, v) =0;

  // Use E_in, y_ion_frac and number density to determine temperature
  double T = get_temperature(y_ion_frac, X_elem_number_density, E_in);
  double oneminusx_dot=0.0; // oneminusx_dot is in units of 1/s
  double Edot=0.0;
  // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH
  // at the end.

  // get neutral fractions
  species_counter=0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    X_neutral_frac[elem] =1.0;
    //cout << "\n neutral_frac=" << neutral_frac << "\n";
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      X_neutral_frac[elem] -=  y_ion_frac[species_counter];
      species_counter ++;
    }
  }
  //rep.printVec("ynow X_neutral_frac",X_neutral_frac,N_elem);
  //rep.printVec("ynow y_ion_frac",y_ion_frac,N_species);
  double yyy[N_equations];
  //for (int v=0;v<N_equations;v++) yyy[v] = NV_Ith_S(y_now,v);
  //rep.printVec("ynow",yyy,N_equations);

  //
  //  ========================================================
  //          Get Rate of Change of Each Species
  //  ========================================================
  //
  /// Start by getting the relevant temperature index:
  if (T > T_max){
    //NOTE RESET ENERGY HERE TOO!!!
    T = T_max;
  }
  else if ( T<T_min){
    T = T_min;
  }
  int temp_index = int (( log10f(T) - log10f(T_min) ) / delta_log_temp );
  double dT = T - Temp_Table[temp_index];


  
  /// ====== Collisional ionisation INTO this species, OUT of previous species ========
  /// this_y_dot(ion) += ionise_rate(im1)*n_e*y(im1) <<< add ionisation from less ionised species to current species
  /// =================================================================================
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    double neutral_frac = X_neutral_frac[elem];   
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists in tracer list
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_ion_rate = ionise_rate_table[y_im1_index_tables[species_counter]] [ temp_index];
        double upper_ion_rate_contrib = dT * ionise_slope_table[ y_im1_index_tables[species_counter]] [temp_index];
        double ionise_rate_im1 = lower_ion_rate + upper_ion_rate_contrib;
        
        //if the less ionised species ISN'T neutral 
        if (y_im1_index_local[species_counter] !=-2){
          double this_y_dot = ionise_rate_im1 *  NV_Ith_S(y_now, y_im1_index_local[species_counter]) *ne;   
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
          NV_Ith_S(y_dot, y_im1_index_local[species_counter]) -= this_y_dot;
          
          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot;
        }
        
        //if the less ionised species IS neutral
        else{
          double this_y_dot = ionise_rate_im1 * neutral_frac * ne;
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot;
        }
      }
      species_counter ++;
    }
  }
  /// ============== Radiative recombination OUT OF this species =====================
  /// y_dot(ion) -= recomb_rate(ion)*n_e*y(ion) <<< subtract recombination to less ionised species
  /// ================================================================================  
  /*species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_recomb_rate = recomb_rate_table [y_ion_index_tables[species_counter] ][temp_index];
        double upper_recomb_rate_contrib = dT * ionise_slope_table[y_ion_index_tables[species_counter] ][temp_index];
        double recomb_rate_ion = lower_recomb_rate + upper_recomb_rate_contrib;

        double this_y_dot = recomb_rate_ion * NV_Ith_S(y_now, y_ion_index_local[species_counter])*ne;
        //Subtract this rate of change from the current species
        NV_Ith_S(y_dot, y_ion_index_local[species_counter]) -= this_y_dot;
        //add this rate of change to the recombination species (provided it isn't neutral)
        /*if (y_im1_index_local[species_counter] !=-2){
          NV_Ith_S(y_dot, y_im1_index_local[species_counter]) += this_y_dot;
        }*/
//         cout << "y_dot [" << species_counter << "] -=" << this_y_dot <<"\n";
//         for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
//         rep.printVec("ydot",dydt,N_equations);
        /// =========  HEATING DUE TO RECOMBINATION OUT OF THIS SPECIES ===========
 /*       double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
        Edot += ion_pot * this_y_dot;
      }
    species_counter ++;
    }
  }*/
  
  
  
  
  

  
 
  /*
  /// ================= Collisional ionisation INTO this species ======================
  /// this_y_dot(ion) += ionise_rate(im1)*n_e*y(im1) <<< add ionisation from less ionised species to current species
  /// =================================================================================
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    double neutral_frac = X_neutral_frac[elem];   
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists in tracer list
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_ion_rate = ionise_rate_table[y_im1_index_tables[species_counter]] [ temp_index];
        double upper_ion_rate_contrib = dT * ionise_slope_table[ y_im1_index_tables[species_counter]] [temp_index];
        double ionise_rate_im1 = lower_ion_rate + upper_ion_rate_contrib;
        
        //if the less ionised species ISN'T neutral 
        if (y_im1_index_local[species_counter] !=-2){
          double this_y_dot = ionise_rate_im1 *  NV_Ith_S(y_now, y_im1_index_local[species_counter]) *ne;   
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
          
          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot;
        }
        
        //if the less ionised species IS neutral
        else{
          double this_y_dot = ionise_rate_im1 * neutral_frac * ne;
          NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;

          /// =========  COOLING DUE TO IONISATION INTO THIS SPECIES ===========
          double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
          Edot -= ion_pot * this_y_dot;
        }
      }
      species_counter ++;
    }
  }
  
  /// ============== Collisional ionisation OUT of this species ======================
  /// y_dot(ion) -= ionise_rate(ion)*n_e*y(ion) <<< subtract ionisation from current species to more ionised species
  /// ================================================================================
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
     
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      // if there does exist a species more ionised in tracer table
      if (y_ip1_index_local[species_counter] != -1){ 
        double lower_ion_rate = ionise_rate_table [y_ion_index_tables[species_counter] ] [temp_index]; //rate at lower bound of temperature
        double upper_ion_rate_contrib = dT * ionise_slope_table [y_ion_index_tables[species_counter] ] [temp_index]; //contribution from upper bound of temperature
        double ionise_rate_ion = lower_ion_rate + upper_ion_rate_contrib; //interpolated rate for exact temperature
        
        double this_y_dot = - ionise_rate_ion * ne * NV_Ith_S(y_now, y_ion_index_local[species_counter]) ;
        NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
        
        /// =========  COOLING DUE TO IONISATION OUT OF THIS SPECIES ===========
        double ion_pot = ionisation_potentials[ y_ion_index_tables[species_counter]];
        Edot -= ion_pot * this_y_dot;
      }
      species_counter ++;
    }
  }
  
  
  /*
  /// ============== Radiative recombination IN to this species ======================
  /// y_dot(ion) += recomb_rate(ip1)*n_e*y(ip1) <<< add recombination from more ionised species to current species
  /// ================================================================================  
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //If there does exist a species more ionised
      if (y_ip1_index_local[species_counter] != -1){ 
        double lower_recomb_rate = recomb_rate_table [y_ip1_index_tables[species_counter]] [temp_index];
        double upper_recomb_rate_contrib = dT * ionise_slope_table [y_ip1_index_tables[species_counter]] [temp_index];
        double recomb_rate_ip1 = lower_recomb_rate + upper_recomb_rate_contrib;
        double this_y_dot = recomb_rate_ip1 * ne * NV_Ith_S(y_now, y_ip1_index_local[species_counter]) ;
        NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
        
        /// =========  HEATING DUE TO RECOMBINATION INTO THIS SPECIES ===========
        double ion_pot = ionisation_potentials[ y_ion_index_tables[species_counter]];
        Edot += ion_pot * this_y_dot;
      }
      species_counter ++;
    }
  }

  /// ============== Radiative recombination OUT OF this species =====================
  /// y_dot(ion) -= recomb_rate(ion)*n_e*y(ion) <<< subtract recombination to less ionised species
  /// ================================================================================  
  species_counter = 0;
  for (int elem=0;elem<N_elem;elem++){//loop over every element
    int N_elem_species=N_species_by_elem[elem];
    
    for (int s=0;s<N_elem_species;s++){//loop over every species in THIS element
      //if the less ionised species exists
      if (y_im1_index_local[species_counter] != -1){ 
        double lower_recomb_rate = recomb_rate_table [y_ion_index_tables[species_counter] ][temp_index];
        double upper_recomb_rate_contrib = dT * ionise_slope_table[y_ion_index_tables[species_counter] ][temp_index];
        double recomb_rate_ion = lower_recomb_rate + upper_recomb_rate_contrib;

        double this_y_dot = - recomb_rate_ion * NV_Ith_S(y_now, y_ion_index_local[species_counter])*ne;
        NV_Ith_S(y_dot, y_ion_index_local[species_counter]) += this_y_dot;
        
        /// =========  HEATING DUE TO RECOMBINATION OUT OF THIS SPECIES ===========
        double ion_pot = ionisation_potentials[ y_im1_index_tables[species_counter]];
        Edot += ion_pot * this_y_dot;
      }
    species_counter ++;
    }
  }*/
  
  //
  // The Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve.
  //
  Edot -= cooling_rate_SD93CIE(T) *ne;
  Edot = 0;
  NV_Ith_S(y_dot,lv_eint) = Edot;

  
  for (int v=0;v<N_equations;v++) dydt[v] = NV_Ith_S(y_dot,v);
  //rep.printVec("ydot returned ",dydt,N_equations);
  
  return 0;
}





