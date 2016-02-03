///
/// \file cooling_SD93_cie.cc
/// \author Jonathan Mackey
///
///
/// Class which calculates the Sutherland & Dopita (1993) cooling
/// rates for collisional ionisation equilibrium.  This is a minimal
/// class which only sets up a spline interpolation and calculates
/// interpolated values.
///
/// Reference: Sutherland & Dopita (1993, ApJS, 88, 253)
///
/// Data file: http://www.mso.anu.edu.au/~ralph/data/cool/m-00.cie
///
/// I am using the normalised curve, so that it can be multiplied
/// by ne*ni to give the appropriate cooling rate.
///
/// I have converted all of the log values to actual values so that I
/// can calculate the spline in linear space.  This is not as good as
/// log-space, but it is sufficiently accurate and runs faster because
/// there is no requirement to convert with log10() and exp() calls
/// at runtime.
///
/// Modifications:\n
/// - 2011.01.14 JM: Written and debugged.
/// - 2011.03.04 JM: Added metal-only cooling for CIE.
/// - 2011.04.12 JM: Changed to log10 interpolation for all cooling curves.
///    Added Wiersma et al (2009) CIE functions, full and metals-only.
/// - 2011.06.20 JM: Got rid of non-ANSI-C exp10 functions
/// - 2012.01.26 JM: replaced 2.303 with GS.ln10()
/// - 2013.08.19 JM: changed all low-T extrapolations to have power
///    law slope of 4.0;
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/interpolate.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/cooling_SD93_cie.h"

using namespace std;

cooling_function_SD93CIE::cooling_function_SD93CIE()
:
  Nspl(91)
{
  Tarray = 0;
  Larray = 0;
  L2array= 0;
  have_set_cooling = false;
}

cooling_function_SD93CIE::~cooling_function_SD93CIE()
{
  Tarray  = mem.myfree(Tarray);
  Larray  = mem.myfree(Larray);
  L2array = mem.myfree(L2array);
}

//
// Set up the spline interpolation arrays, for solar metallicity curve.
//
void cooling_function_SD93CIE::setup_SD93_cie()
{
  if (have_set_cooling) {
    rep.error("setup_SD93_cie() only one cooling function in cooling_function_SD93CIE",1);
  }
  cout <<"\t\t----------------------------------------------------\n";
  cout <<"\t\tSetting up Sutherland & Dopita (1993, ApJS, 88, 253) CIE cooling function\n";
  cout <<"\t\tfrom file m-00.cie from http://www.mso.anu.edu.au/~ralph/data/cool/\n";

  Tarray  = mem.myalloc(Tarray, Nspl);
  Larray  = mem.myalloc(Larray, Nspl);
  L2array = mem.myalloc(L2array,Nspl);
  
  double t1[91] = 
    {1.000000e+04, 1.122018e+04, 1.258925e+04, 1.412538e+04,
     1.584893e+04, 1.778279e+04, 1.995262e+04, 2.238721e+04,
     2.511886e+04, 2.818383e+04, 3.162278e+04, 3.548134e+04,
     3.981072e+04, 4.466836e+04, 5.011872e+04, 5.623413e+04,
     6.309573e+04, 7.079458e+04, 7.943282e+04, 8.912509e+04,
     1.000000e+05, 1.122018e+05, 1.258925e+05, 1.412538e+05,
     1.584893e+05, 1.778279e+05, 1.995262e+05, 2.238721e+05,
     2.511886e+05, 2.818383e+05, 3.162278e+05, 3.548134e+05,
     3.981072e+05, 4.466836e+05, 5.011872e+05, 5.623413e+05,
     6.309573e+05, 7.079458e+05, 7.943282e+05, 8.912509e+05,
     1.000000e+06, 1.122018e+06, 1.258925e+06, 1.412538e+06,
     1.584893e+06, 1.778279e+06, 1.995262e+06, 2.238721e+06,
     2.511886e+06, 2.818383e+06, 3.162278e+06, 3.548134e+06,
     3.981072e+06, 4.466836e+06, 5.011872e+06, 5.623413e+06,
     6.309573e+06, 7.079458e+06, 7.943282e+06, 8.912509e+06,
     1.000000e+07, 1.122018e+07, 1.258925e+07, 1.412538e+07,
     1.584893e+07, 1.778279e+07, 1.995262e+07, 2.238721e+07,
     2.511886e+07, 2.818383e+07, 3.162278e+07, 3.548134e+07,
     3.981072e+07, 4.466836e+07, 5.011872e+07, 5.623413e+07,
     6.309573e+07, 7.079458e+07, 7.943282e+07, 8.912509e+07,
     1.000000e+08, 1.122018e+08, 1.258925e+08, 1.412538e+08,
     1.584893e+08, 1.778279e+08, 1.995262e+08, 2.238721e+08,
     2.511886e+08, 2.818383e+08, 3.162278e+08};
  for (int v=0;v<Nspl;v++)
    Tarray[v] = log10(t1[v]);

  double t2[91] = 
    {8.709636e-24, 3.467369e-23, 6.760830e-23, 1.202264e-22,
     1.621810e-22, 1.584893e-22, 1.380384e-22, 1.258925e-22,
     1.318257e-22, 1.513561e-22, 1.862087e-22, 2.344229e-22,
     2.951209e-22, 3.801894e-22, 4.786301e-22, 6.025596e-22,
     7.244360e-22, 8.511380e-22, 9.772372e-22, 1.047129e-21,
     1.023293e-21, 9.549926e-22, 9.332543e-22, 9.772372e-22,
     1.047129e-21, 1.071519e-21, 1.096478e-21, 1.096478e-21,
     1.023293e-21, 7.413102e-22, 4.466836e-22, 2.818383e-22,
     2.187762e-22, 1.949845e-22, 1.949845e-22, 1.949845e-22,
     1.737801e-22, 1.380384e-22, 1.174898e-22, 1.122018e-22,
     1.096478e-22, 1.096478e-22, 1.096478e-22, 1.122018e-22,
     1.148154e-22, 1.071519e-22, 8.511380e-23, 6.309573e-23,
     4.897788e-23, 4.073803e-23, 3.630781e-23, 3.311311e-23,
     3.162278e-23, 2.951209e-23, 2.754229e-23, 2.570396e-23,
     2.511886e-23, 2.511886e-23, 2.570396e-23, 2.691535e-23,
     2.691535e-23, 2.570396e-23, 2.398833e-23, 2.238721e-23,
     2.089296e-23, 1.995262e-23, 1.905461e-23, 1.862087e-23,
     1.862087e-23, 1.862087e-23, 1.862087e-23, 1.905461e-23,
     1.949845e-23, 1.995262e-23, 2.089296e-23, 2.137962e-23,
     2.238721e-23, 2.290868e-23, 2.398833e-23, 2.511886e-23,
     2.630268e-23, 2.754229e-23, 2.884032e-23, 2.951209e-23,
     3.090295e-23, 3.235937e-23, 3.388442e-23, 3.548134e-23,
     3.715352e-23, 3.981072e-23, 4.168694e-23};
  for (int v=0;v<Nspl;v++)
    Larray[v] = log10(t2[v]);
 
  MinTemp = Tarray[0];
  MaxTemp = Tarray[Nspl-1];
  //
  // logarithmic slope for extrapolation to lower temperatures.
  //
  MinSlope = 8.0 ; // was too shallow! (Larray[1]-Larray[0])/(Tarray[1]-Tarray[0]);
  //
  // logarthmic slope for extrapolation to higher temperatures.
  //
  MaxSlope = (Larray[Nspl-1]-Larray[Nspl-2])/(Tarray[Nspl-1]-Tarray[Nspl-2]);
  //MaxSlope = 0.0;

  cout << "\t\t min-slope="<<MinSlope<<" max-slope="<<MaxSlope<<"\n";

  interpolate.spline(Tarray, Larray, Nspl, 1.e99, 1.e99, L2array);

  have_set_cooling = true;
  cout <<"\t\t----------------------------------------------------\n";
#ifdef TESTING
  ofstream outf("cooling_SD93_cie_solar.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"SD93 Cooling: Temperature(K) Rate(cm^3/s) \n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< cooling_rate_SD93CIE(t)<<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING
  return;
}

//
// Set up the spline interpolation arrays, for Metals--only curve.
//
void cooling_function_SD93CIE::setup_SD93_cie_OnlyMetals()
{
  if (have_set_cooling) {
    rep.error("setup_SD93_OnlyMetals() only one cooling function in cooling_function_SD93CIE",1);
  }
  cout <<"\t\t----------------------------------------------------\n";
  cout <<"\t\tSetting up Sutherland & Dopita (1993, ApJS, 88, 253) CIE cooling function\n";
  cout <<"\t\tfrom file m-00.cie from http://www.mso.anu.edu.au/~ralph/data/cool/\n";
  cout <<"\t\tand file mzero.cie from http://www.mso.anu.edu.au/~ralph/data/cool/\n";

  Tarray  = mem.myalloc(Tarray, Nspl);
  Larray  = mem.myalloc(Larray, Nspl);
  L2array = mem.myalloc(L2array,Nspl);
  
  double t1[91] = 
    {1.000000e+04, 1.122018e+04, 1.258925e+04, 1.412538e+04, \
     1.584893e+04, 1.778279e+04, 1.995262e+04, 2.238721e+04, \
     2.511886e+04, 2.818383e+04, 3.162278e+04, 3.548134e+04, \
     3.981072e+04, 4.466836e+04, 5.011872e+04, 5.623413e+04, \
     6.309573e+04, 7.079458e+04, 7.943282e+04, 8.912509e+04, \
     1.000000e+05, 1.122018e+05, 1.258925e+05, 1.412538e+05, \
     1.584893e+05, 1.778279e+05, 1.995262e+05, 2.238721e+05, \
     2.511886e+05, 2.818383e+05, 3.162278e+05, 3.548134e+05, \
     3.981072e+05, 4.466836e+05, 5.011872e+05, 5.623413e+05, \
     6.309573e+05, 7.079458e+05, 7.943282e+05, 8.912509e+05, \
     1.000000e+06, 1.122018e+06, 1.258925e+06, 1.412538e+06, \
     1.584893e+06, 1.778279e+06, 1.995262e+06, 2.238721e+06, \
     2.511886e+06, 2.818383e+06, 3.162278e+06, 3.548134e+06, \
     3.981072e+06, 4.466836e+06, 5.011872e+06, 5.623413e+06, \
     6.309573e+06, 7.079458e+06, 7.943282e+06, 8.912509e+06, \
     1.000000e+07, 1.122018e+07, 1.258925e+07, 1.412538e+07, \
     1.584893e+07, 1.778279e+07, 1.995262e+07, 2.238721e+07, \
     2.511886e+07, 2.818383e+07, 3.162278e+07, 3.548134e+07, \
     3.981072e+07, 4.466836e+07, 5.011872e+07, 5.623413e+07, \
     6.309573e+07, 7.079458e+07, 7.943282e+07, 8.912509e+07, \
     1.000000e+08, 1.122018e+08, 1.258925e+08, 1.412538e+08, \
     1.584893e+08, 1.778279e+08, 1.995262e+08, 2.238721e+08, \
     2.511886e+08, 2.818383e+08, 3.162278e+08};
    
  for (int v=0;v<Nspl;v++)
    Tarray[v] = log10(t1[v]);

  double t2[91] =
    {4.819184e-24, 2.118406e-23, 2.779758e-23, 3.313008e-23, \
     3.628847e-23, 3.826288e-23, 4.253917e-23, 5.509796e-23, \
     7.687159e-23, 1.077045e-22, 1.507274e-22, 2.055826e-22, \
     2.716786e-22, 3.602368e-22, 4.591316e-22, 5.791173e-22, \
     6.905515e-22, 8.010193e-22, 9.183529e-22, 9.908944e-22, \
     9.765195e-22, 9.186848e-22, 9.044140e-22, 9.543285e-22, \
     1.028932e-21, 1.056384e-21, 1.083889e-21, 1.085763e-21, \
     1.013960e-21, 7.329926e-22, 4.392705e-22, 2.750775e-22, \
     2.126102e-22, 1.892301e-22, 1.894891e-22, 1.897364e-22, \
     1.687682e-22, 1.331406e-22, 1.127035e-22, 1.075245e-22, \
     1.049705e-22, 1.049705e-22, 1.049705e-22, 1.075245e-22, \
     1.100291e-22, 1.023656e-22, 8.021602e-23, 5.819795e-23, \
     4.396601e-23, 3.560941e-23, 3.105973e-23, 2.774279e-23, \
     2.599936e-23, 2.375769e-23, 2.151669e-23, 1.953801e-23, \
     1.866232e-23, 1.835803e-23, 1.878565e-23, 1.967099e-23, \
     1.932957e-23, 1.794149e-23, 1.586002e-23, 1.387583e-23, \
     1.198045e-23, 1.062008e-23, 9.282235e-24, 8.387941e-24, \
     7.656089e-24, 7.139335e-24, 6.598227e-24, 6.465353e-24, \
     6.315879e-24, 6.148781e-24, 6.438564e-24, 5.891455e-24, \
     6.169110e-24, 5.926240e-24, 6.205535e-24, 6.497993e-24, \
     6.350057e-24, 6.649326e-24, 6.962699e-24, 6.069804e-24, \
     6.355865e-24, 6.655408e-24, 6.969068e-24, 6.641024e-24, \
     6.954006e-24, 8.187940e-24, 7.802523e-24};

  for (int v=0;v<Nspl;v++)
    Larray[v] = log10(t2[v]);
 
  MinTemp = Tarray[0];
  MaxTemp = Tarray[Nspl-1];
  //
  // logarithmic slope for extrapolation to lower temperatures.
  //
  MinSlope = 8.0; // was too shallow! (Larray[1]-Larray[0])/(Tarray[1]-Tarray[0]);
  //
  // logarthmic slope for extrapolation to higher temperatures.
  // This is set to zero, because the curve is noisy at high T due
  // to subtracting two very similar numbers.
  //
  MaxSlope = 0.0;

  cout << "\t\t min-slope="<<MinSlope<<" max-slope="<<MaxSlope<<"\n";

  interpolate.spline(Tarray, Larray, Nspl, 1.e99, 1.e99, L2array);

  have_set_cooling = true;
#ifdef TESTING
  ofstream outf("cooling_SD93_cie_metalsonly.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"SD93 Cooling: Temperature(K) Rate(cm^3/s) \n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< cooling_rate_SD93CIE(t)<<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING
  cout <<"\t\t----------------------------------------------------\n";
  return;
}

//
// Set up the spline interpolation arrays, for metal-free curve.
//
void cooling_function_SD93CIE::setup_SD93_cie_MetalFree()
{
  if (have_set_cooling) {
    rep.error("setup_SD93_MetalFree() only one cooling function in cooling_function_SD93CIE",1);
  }
  cout <<"\t\t----------------------------------------------------\n";
  cout <<"\t\tSetting up Sutherland&Dopita(1993,ApJS,88,253) METAL-FREE-CIE cooling function\n";
  cout <<"\t\tfrom file m-zero.cie from http://www.mso.anu.edu.au/~ralph/data/cool/\n";

  Tarray  = mem.myalloc(Tarray, Nspl);
  Larray  = mem.myalloc(Larray, Nspl);
  L2array = mem.myalloc(L2array,Nspl);
  
  double t1[91] =
    {1.000000e+04, 1.122018e+04, 1.258925e+04, 1.412538e+04, \
     1.584893e+04, 1.778279e+04, 1.995262e+04, 2.238721e+04, \
     2.511886e+04, 2.818383e+04, 3.162278e+04, 3.548134e+04, \
     3.981072e+04, 4.466836e+04, 5.011872e+04, 5.623413e+04, \
     6.309573e+04, 7.079458e+04, 7.943282e+04, 8.912509e+04, \
     1.000000e+05, 1.122018e+05, 1.258925e+05, 1.412538e+05, \
     1.584893e+05, 1.778279e+05, 1.995262e+05, 2.238721e+05, \
     2.511886e+05, 2.818383e+05, 3.162278e+05, 3.548134e+05, \
     3.981072e+05, 4.466836e+05, 5.011872e+05, 5.623413e+05, \
     6.309573e+05, 7.079458e+05, 7.943282e+05, 8.912509e+05, \
     1.000000e+06, 1.122018e+06, 1.258925e+06, 1.412538e+06, \
     1.584893e+06, 1.778279e+06, 1.995262e+06, 2.238721e+06, \
     2.511886e+06, 2.818383e+06, 3.162278e+06, 3.548134e+06, \
     3.981072e+06, 4.466836e+06, 5.011872e+06, 5.623413e+06, \
     6.309573e+06, 7.079458e+06, 7.943282e+06, 8.912509e+06, \
     1.000000e+07, 1.122018e+07, 1.258925e+07, 1.412538e+07, \
     1.584893e+07, 1.778279e+07, 1.995262e+07, 2.238721e+07, \
     2.511886e+07, 2.818383e+07, 3.162278e+07, 3.548134e+07, \
     3.981072e+07, 4.466836e+07, 5.011872e+07, 5.623413e+07, \
     6.309573e+07, 7.079458e+07, 7.943282e+07, 8.912509e+07, \
     1.000000e+08, 1.122018e+08, 1.258925e+08, 1.412538e+08, \
     1.584893e+08, 1.778279e+08, 1.995262e+08, 2.238721e+08, \
     2.511886e+08, 2.818383e+08, 3.162278e+08};
    
  for (int v=0;v<Nspl;v++)
    Tarray[v] = log10(t1[v]);

  double t2[91] = 
    {3.890451e-24, 1.348963e-23, 3.981072e-23, 8.709636e-23, \
     1.258925e-22, 1.202264e-22, 9.549926e-23, 7.079458e-23, \
     5.495409e-23, 4.365158e-23, 3.548134e-23, 2.884032e-23, \
     2.344229e-23, 1.995262e-23, 1.949845e-23, 2.344229e-23, \
     3.388442e-23, 5.011872e-23, 5.888437e-23, 5.623413e-23, \
     4.677351e-23, 3.630781e-23, 2.884032e-23, 2.290868e-23, \
     1.819701e-23, 1.513561e-23, 1.258925e-23, 1.071519e-23, \
     9.332543e-24, 8.317638e-24, 7.413102e-24, 6.760830e-24, \
     6.165950e-24, 5.754399e-24, 5.495409e-24, 5.248075e-24, \
     5.011872e-24, 4.897788e-24, 4.786301e-24, 4.677351e-24, \
     4.677351e-24, 4.677351e-24, 4.677351e-24, 4.677351e-24, \
     4.786301e-24, 4.786301e-24, 4.897788e-24, 4.897788e-24, \
     5.011872e-24, 5.128614e-24, 5.248075e-24, 5.370318e-24, \
     5.623413e-24, 5.754399e-24, 6.025596e-24, 6.165950e-24, \
     6.456542e-24, 6.760830e-24, 6.918310e-24, 7.244360e-24, \
     7.585776e-24, 7.762471e-24, 8.128305e-24, 8.511380e-24, \
     8.912509e-24, 9.332543e-24, 9.772372e-24, 1.023293e-23, \
     1.096478e-23, 1.148154e-23, 1.202264e-23, 1.258925e-23, \
     1.318257e-23, 1.380384e-23, 1.445440e-23, 1.548817e-23, \
     1.621810e-23, 1.698244e-23, 1.778279e-23, 1.862087e-23, \
     1.995262e-23, 2.089296e-23, 2.187762e-23, 2.344229e-23, \
     2.454709e-23, 2.570396e-23, 2.691535e-23, 2.884032e-23, \
     3.019952e-23, 3.162278e-23, 3.388442e-23};
  
  for (int v=0;v<Nspl;v++)
    Larray[v] = log10(t2[v]);
 
  MinTemp = Tarray[0];
  MaxTemp = Tarray[Nspl-1];
  //
  // logarithmic slope for extrapolation to lower temperatures.
  //
  MinSlope = 8.0; // was too shallow! (Larray[1]-Larray[0])/(Tarray[1]-Tarray[0]);
  //
  // logarthmic slope for extrapolation to higher temperatures.
  //
  MaxSlope = (Larray[Nspl-1]-Larray[Nspl-2])/(Tarray[Nspl-1]-Tarray[Nspl-2]);
  //MaxSlope = 0.0;

  cout << "\t\t min-slope="<<MinSlope<<" max-slope="<<MaxSlope<<"\n";

  interpolate.spline(Tarray, Larray, Nspl, 1.e99, 1.e99, L2array);

  have_set_cooling = true;
  cout <<"\t\t----------------------------------------------------\n";
#ifdef TESTING
  ofstream outf("cooling_SD93_cie_metalfree.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"SD93 Cooling: Temperature(K) Rate(cm^3/s) \n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< cooling_rate_SD93CIE(t)<<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING
  return;
}

//
// Sets up spline interpolation for metals-only cooling from
// Wiersma, R.P.C., Schaye, J., and Smith, B.D. (2009,MNRAS,393,99)
// (arxiv:0807.3748).  This CIE function is from 
// http://www.strw.leidenuniv.nl/WSS08/z_collis.txt
// and is a factor of a few lower than SD93-CIE.
//
void cooling_function_SD93CIE::setup_WSS09_CIE_OnlyMetals()
{
  if (have_set_cooling) {
    rep.error("setup_WSS09_CIE_OnlyMetals() only one cooling function in cooling_function_SD93CIE",1);
  }
  cout <<"\t\t----------------------------------------------------\n";
  cout <<"\t\tSetting up  Wiersma et al. (2009,MNRAS,393,99) CIE cooling function\n";
  cout <<"\t\tMETALS-ONLY CURVE, resampled onto a 91 point cubic spline interpolation.\n";
  cout <<"\t\tFile from file http://www.strw.leidenuniv.nl/WSS08/z_collis.txt\n";

  Tarray  = mem.myalloc(Tarray, Nspl);
  Larray  = mem.myalloc(Larray, Nspl);
  L2array = mem.myalloc(L2array,Nspl);
  
  double t1[91] = {
    2.00000000e+00, 2.07757611e+00, 2.15515223e+00, 2.23272834e+00, 2.31030446e+00, \
    2.38788057e+00, 2.46545669e+00, 2.54303280e+00, 2.62060892e+00, 2.69818503e+00, \
    2.77576115e+00, 2.85333726e+00, 2.93091337e+00, 3.00848949e+00, 3.08606560e+00, \
    3.16364172e+00, 3.24121783e+00, 3.31879395e+00, 3.39637006e+00, 3.47394618e+00, \
    3.55152229e+00, 3.62909840e+00, 3.70667452e+00, 3.78425063e+00, 3.86182675e+00, \
    3.93940286e+00, 4.01697898e+00, 4.09455509e+00, 4.17213121e+00, 4.24970732e+00, \
    4.32728344e+00, 4.40485955e+00, 4.48243566e+00, 4.56001178e+00, 4.63758789e+00, \
    4.71516401e+00, 4.79274012e+00, 4.87031624e+00, 4.94789235e+00, 5.02546847e+00, \
    5.10304458e+00, 5.18062070e+00, 5.25819681e+00, 5.33577292e+00, 5.41334904e+00, \
    5.49092515e+00, 5.56850127e+00, 5.64607738e+00, 5.72365350e+00, 5.80122961e+00, \
    5.87880573e+00, 5.95638184e+00, 6.03395795e+00, 6.11153407e+00, 6.18911018e+00, \
    6.26668630e+00, 6.34426241e+00, 6.42183853e+00, 6.49941464e+00, 6.57699076e+00, \
    6.65456687e+00, 6.73214299e+00, 6.80971910e+00, 6.88729521e+00, 6.96487133e+00, \
    7.04244744e+00, 7.12002356e+00, 7.19759967e+00, 7.27517579e+00, 7.35275190e+00, \
    7.43032802e+00, 7.50790413e+00, 7.58548025e+00, 7.66305636e+00, 7.74063247e+00, \
    7.81820859e+00, 7.89578470e+00, 7.97336082e+00, 8.05093693e+00, 8.12851305e+00, \
    8.20608916e+00, 8.28366528e+00, 8.36124139e+00, 8.43881750e+00, 8.51639362e+00, \
    8.59396973e+00, 8.67154585e+00, 8.74912196e+00, 8.82669808e+00, 8.90427419e+00, \
    8.98185031e+00};
  for (int v=0;v<Nspl;v++)
    Tarray[v] = t1[v];

  double t2[91] = {
    -2.69042032e+01, -2.68339466e+01, -2.67628015e+01, -2.66852365e+01, -2.66026698e+01, \
    -2.65218150e+01, -2.64469693e+01, -2.63761355e+01, -2.63097777e+01, -2.62474256e+01, \
    -2.61886746e+01, -2.61332877e+01, -2.60808330e+01, -2.60309113e+01, -2.59830826e+01, \
    -2.59369007e+01, -2.58919300e+01, -2.58476214e+01, -2.58031708e+01, -2.57581287e+01, \
    -2.57139260e+01, -2.56680924e+01, -2.56216866e+01, -2.55784123e+01, -2.55358056e+01, \
    -2.54579940e+01, -2.52789911e+01, -2.42634880e+01, -2.31979645e+01, -2.27183209e+01, \
    -2.25726495e+01, -2.24284223e+01, -2.22590643e+01, -2.20877851e+01, -2.19241810e+01, \
    -2.17723986e+01, -2.16330514e+01, -2.15062964e+01, -2.14071669e+01, -2.13475926e+01, \
    -2.13492162e+01, -2.13325337e+01, -2.13034976e+01, -2.12874309e+01, -2.13074247e+01, \
    -2.14856951e+01, -2.16658156e+01, -2.17176117e+01, -2.17351658e+01, -2.17860161e+01, \
    -2.18142313e+01, -2.18029824e+01, -2.18098104e+01, -2.18455343e+01, -2.19092400e+01, \
    -2.20294769e+01, -2.21901200e+01, -2.23345038e+01, -2.24678858e+01, -2.25823022e+01, \
    -2.26539966e+01, -2.26847250e+01, -2.26876913e+01, -2.26767177e+01, -2.26732880e+01, \
    -2.26964528e+01, -2.27613667e+01, -2.28719040e+01, -2.30037799e+01, -2.31212437e+01, \
    -2.32122653e+01, -2.32778695e+01, -2.33214754e+01, -2.33486237e+01, -2.33630273e+01, \
    -2.33677304e+01, -2.33656718e+01, -2.33604758e+01, -2.33515894e+01, -2.33410786e+01, \
    -2.33304239e+01, -2.33191682e+01, -2.33067658e+01, -2.32928461e+01, -2.32761560e+01, \
    -2.32529092e+01, -2.32280201e+01, -2.32018214e+01, -2.31746034e+01, -2.31467139e+01, \
    -2.31183757e+01};
  for (int v=0;v<Nspl;v++)
    Larray[v] = t2[v];
 
  MinTemp = Tarray[0];
  MaxTemp = Tarray[Nspl-1];
  //
  // logarithmic slope for extrapolation to lower temperatures.
  //
  MinSlope = 8.0; // was too shallow! (Larray[1]-Larray[0])/(Tarray[1]-Tarray[0]);
  //
  // logarthmic slope for extrapolation to higher temperatures.
  //
  MaxSlope = (Larray[Nspl-1]-Larray[Nspl-2])/(Tarray[Nspl-1]-Tarray[Nspl-2]);
  //MaxSlope = 0.0;

  cout << "\t\t min-slope="<<MinSlope<<" max-slope="<<MaxSlope<<"\n";

  interpolate.spline(Tarray, Larray, Nspl, 1.e99, 1.e99, L2array);

  have_set_cooling = true;
#ifdef TESTING
  ofstream outf("cooling_WSS09_CIE_metalsonly.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"WSS09 metals-only Cooling: Temperature(K) Rate(cm^3/s) \n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< cooling_rate_SD93CIE(t)<<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING
  cout <<"\t\t----------------------------------------------------\n";
  return;
}

//
// Sets up spline interpolation for total cooling from
// Wiersma, R.P.C., Schaye, J., and Smith, B.D. (2009,MNRAS,393,99)
// (arxiv:0807.3748).  This CIE function is from 
// http://www.strw.leidenuniv.nl/WSS08/z_collis.txt
// and is a factor of a few lower than SD93-CIE.
//
void cooling_function_SD93CIE::setup_WSS09_CIE()
{
  if (have_set_cooling) {
    rep.error("setup_WSS09_CIE() only one cooling function in cooling_function_SD93CIE",1);
  }
  cout <<"\t\t----------------------------------------------------\n";
  cout <<"\t\tSetting up  Wiersma et al. (2009,MNRAS,393,99) CIE cooling function\n";
  cout <<"\t\tTOTAL COOLING CURVE, resampled onto a 91 point cubic spline interpolation.\n";
  cout <<"\t\tFile from file http://www.strw.leidenuniv.nl/WSS08/z_collis.txt\n";

  Tarray  = mem.myalloc(Tarray, Nspl);
  Larray  = mem.myalloc(Larray, Nspl);
  L2array = mem.myalloc(L2array,Nspl);
  
  double t1[91] = {
    3.01443653e+00, 3.08074113e+00, 3.14704572e+00, 3.21335032e+00, 3.27965492e+00, \
    3.34595952e+00, 3.41226411e+00, 3.47856871e+00, 3.54487331e+00, 3.61117791e+00, \
    3.67748250e+00, 3.74378710e+00, 3.81009170e+00, 3.87639630e+00, 3.94270089e+00, \
    4.00900549e+00, 4.07531009e+00, 4.14161469e+00, 4.20791928e+00, 4.27422388e+00, \
    4.34052848e+00, 4.40683308e+00, 4.47313767e+00, 4.53944227e+00, 4.60574687e+00, \
    4.67205147e+00, 4.73835606e+00, 4.80466066e+00, 4.87096526e+00, 4.93726986e+00, \
    5.00357445e+00, 5.06987905e+00, 5.13618365e+00, 5.20248825e+00, 5.26879284e+00, \
    5.33509744e+00, 5.40140204e+00, 5.46770664e+00, 5.53401123e+00, 5.60031583e+00, \
    5.66662043e+00, 5.73292503e+00, 5.79922962e+00, 5.86553422e+00, 5.93183882e+00, \
    5.99814342e+00, 6.06444801e+00, 6.13075261e+00, 6.19705721e+00, 6.26336181e+00, \
    6.32966640e+00, 6.39597100e+00, 6.46227560e+00, 6.52858020e+00, 6.59488479e+00, \
    6.66118939e+00, 6.72749399e+00, 6.79379859e+00, 6.86010319e+00, 6.92640778e+00, \
    6.99271238e+00, 7.05901698e+00, 7.12532158e+00, 7.19162617e+00, 7.25793077e+00, \
    7.32423537e+00, 7.39053997e+00, 7.45684456e+00, 7.52314916e+00, 7.58945376e+00, \
    7.65575836e+00, 7.72206295e+00, 7.78836755e+00, 7.85467215e+00, 7.92097675e+00, \
    7.98728134e+00, 8.05358594e+00, 8.11989054e+00, 8.18619514e+00, 8.25249973e+00, \
    8.31880433e+00, 8.38510893e+00, 8.45141353e+00, 8.51771812e+00, 8.58402272e+00, \
    8.65032732e+00, 8.71663192e+00, 8.78293651e+00, 8.84924111e+00, 8.91554571e+00, \
    8.98185031e+00};
  for (int v=0;v<Nspl;v++)
    Tarray[v] = t1[v];

  double t2[91] = {
    -2.60271176e+01, -2.59862393e+01, -2.59465765e+01, -2.59078199e+01, -2.58697674e+01, \
    -2.58318686e+01, -2.57937411e+01, -2.57551907e+01, -2.57174320e+01, -2.56784011e+01, \
    -2.56384384e+01, -2.55993821e+01, -2.55639829e+01, -2.55237510e+01, -2.54465489e+01, \
    -2.52104553e+01, -2.40439145e+01, -2.27106891e+01, -2.19654772e+01, -2.18194017e+01, \
    -2.18879583e+01, -2.19531192e+01, -2.19590805e+01, -2.19275436e+01, -2.18649150e+01, \
    -2.17724395e+01, -2.16558347e+01, -2.15189076e+01, -2.13864462e+01, -2.13078930e+01, \
    -2.12764474e+01, -2.12848711e+01, -2.13027334e+01, -2.12907589e+01, -2.12777196e+01, \
    -2.12696250e+01, -2.12814920e+01, -2.14013058e+01, -2.15835253e+01, -2.16770509e+01, \
    -2.17018250e+01, -2.17207498e+01, -2.17660049e+01, -2.17952703e+01, -2.17895001e+01, \
    -2.17857233e+01, -2.18042281e+01, -2.18395524e+01, -2.18972244e+01, -2.19956993e+01, \
    -2.21221663e+01, -2.22372740e+01, -2.23363626e+01, -2.24257086e+01, -2.24910629e+01, \
    -2.25272544e+01, -2.25396826e+01, -2.25365125e+01, -2.25236642e+01, -2.25106960e+01, \
    -2.25064580e+01, -2.25154638e+01, -2.25421674e+01, -2.25822164e+01, -2.26216950e+01, \
    -2.26477043e+01, -2.26587787e+01, -2.26584328e+01, -2.26488507e+01, -2.26308448e+01, \
    -2.26093730e+01, -2.25856417e+01, -2.25605022e+01, -2.25345761e+01, -2.25083734e+01, \
    -2.24823140e+01, -2.24523277e+01, -2.24215616e+01, -2.23910551e+01, -2.23608462e+01, \
    -2.23309000e+01, -2.23012518e+01, -2.22719029e+01, -2.22415250e+01, -2.22073473e+01, \
    -2.21733037e+01, -2.21393357e+01, -2.21054096e+01, -2.20714890e+01, -2.20374854e+01, \
    -2.20032642e+01};
  for (int v=0;v<Nspl;v++)
    Larray[v] = t2[v];
 
  MinTemp = Tarray[0];
  MaxTemp = Tarray[Nspl-1];
  //
  // logarithmic slope for extrapolation to lower temperatures.
  //
  MinSlope = 8.0; // was too shallow! (Larray[1]-Larray[0])/(Tarray[1]-Tarray[0]);
  //
  // logarthmic slope for extrapolation to higher temperatures.
  //
  MaxSlope = (Larray[Nspl-1]-Larray[Nspl-2])/(Tarray[Nspl-1]-Tarray[Nspl-2]);
  //MaxSlope = 0.0;

  cout << "\t\t min-slope="<<MinSlope<<" max-slope="<<MaxSlope<<"\n";

  interpolate.spline(Tarray, Larray, Nspl, 1.e99, 1.e99, L2array);

  have_set_cooling = true;
#ifdef TESTING
  ofstream outf("cooling_WSS09_CIE_total.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"WSS09 metals-only Cooling: Temperature(K) Rate(cm^3/s) \n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< cooling_rate_SD93CIE(t) <<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING
  cout <<"\t\t----------------------------------------------------\n";
  return;
}



//
// Calculate rate for a given temperature.  Result is returned in
// units of erg.cm^{3}.s^{-1}
//
double cooling_function_SD93CIE::cooling_rate_SD93CIE(double T
			      ///< Input Temperature.
			      )
{
  //
  // Since adaptive integrators can overshoot, we need to be able to
  // silently deal with negative or infinite temperatures:
  //
  if (T<0.0 || !isfinite(T))
    return HUGE_VAL;

  double rate=0.0;
  //
  // reset temperature to log-space for lookups in tables.
  //
  T = log10(T);
  //
  // Three cases: T>Tmax, T<Tmin, or T in between.
  // For in-between we can lookup the value with splint().
  //
  // We extrapolate according to L(T) = L(Tmax)*(T/Tmax)^alpha
  // in Log space this is log(L) = log(Lmax) +alpha*(log(T)-log(Tmax))
  //
  // At low temperatures we decrease cooling as a power law.
  //
  if (T>MaxTemp) {
    //rate = Larray[Nspl-1] *pow(T/MaxTemp, MaxSlope);
    rate = Larray[Nspl-1] +MaxSlope*(T-MaxTemp);
  }
  else if (T<MinTemp) {
    //rate = Larray[0] *pow(T/MinTemp, MinSlope);
    rate = Larray[0] +MinSlope*(T-MinTemp);
  }
  else {
    interpolate.splint(Tarray, Larray, L2array, Nspl, T, &rate);
  }

  return exp(pconst.ln10()*rate);
}

