/// \file hydrogen_recomb_Hummer94.cc
/// \author Jonathan Mackey
///
/// Set up Spline interpolation table for Hummer (1994) recombination rate
/// and recombination cooling rate for Hydrogen (with and without free-free).
///
/// Modifications:
/// - 2011.03.04 JM: Wrote class.
/// - 2011.03.14 JM: Renamed file.
/// - 2011.04.12 JM: put the main() function in an ifdef for testing.
/// - 2011.04.18 JM: Added an ifdef for RT_TEST_PROBS to return 2.59e-13 as the RRR.
/// - 2011.05.10 JM: Output cooling rates only if myrank==0 for parallel (so processes
///    don't fight over the file and slow down the code (by a lot!)).
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/interpolate.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

//#define TEST_HUMMER94_COOLING_FUNCTION
#include "global.h"
#include "microphysics/hydrogen_recomb_Hummer94.h"
using namespace std;


Hummer94_Hrecomb::Hummer94_Hrecomb()
: kB(1.381e-16)
{
  //
  // Use Hummer (1994) recombination rates and recombination energy loss rates.  These
  // are given as tables, so we need to fit cubic splines to the data.
  //
  hr_t=hr_alpha=hr_alpha2=hr_beta=hr_beta2=0;
  hr_Nspl = 31;

  hr_t      = mem.myalloc(hr_t      , hr_Nspl);
  hr_alpha  = mem.myalloc(hr_alpha  , hr_Nspl);
  hr_alpha2 = mem.myalloc(hr_alpha2 , hr_Nspl);
  hr_beta   = mem.myalloc(hr_beta   , hr_Nspl);
  hr_beta2  = mem.myalloc(hr_beta2  , hr_Nspl);
  hr_btot   = mem.myalloc(hr_btot   , hr_Nspl);
  hr_btot2  = mem.myalloc(hr_btot2  , hr_Nspl);

  double caseb[31] = {9.283e-11, 8.823e-11, 8.361e-11, 7.898e-11, 7.435e-11,
		      6.973e-11, 6.512e-11, 6.054e-11, 5.599e-11, 5.147e-11,
		      4.700e-11, 4.258e-11, 3.823e-11, 3.397e-11, 2.983e-11,
		      2.584e-11, 2.204e-11, 1.847e-11, 1.520e-11, 1.226e-11,
		      9.696e-12, 7.514e-12, 5.710e-12, 4.257e-12, 3.117e-12,
		      2.244e-12, 1.590e-12, 1.110e-12, 7.642e-13, 5.199e-13,
		      3.498e-13};
  double coolb[31] = {8.287e-11, 7.821e-11, 7.356e-11, 6.892e-11, 6.430e-11,
		      5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, 4.170e-11,
		      3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11,
		      1.736e-11, 1.402e-11, 1.103e-11, 8.442e-12, 6.279e-12,
		      4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13,
		      6.023e-13, 3.738e-13, 2.268e-13, 1.348e-13, 7.859e-14,
		      4.499e-14};
  double CoolTot[31]={9.348e-11, 8.889e-11, 8.432e-11, 7.977e-11, 7.525e-11, 7.077e-11,
                      6.633e-11, 6.194e-11, 5.758e-11, 5.332e-11, 4.915e-11, 4.508e-11,
                      4.112e-11, 3.733e-11, 3.373e-11, 3.039e-11, 2.737e-11, 2.472e-11,
                      2.247e-11, 2.062e-11, 1.914e-11, 1.797e-11, 1.704e-11, 1.628e-11,
                      1.563e-11, 1.505e-11, 1.451e-11, 1.402e-11, 1.358e-11, 1.318e-11,
                      1.285e-11};
  //
  // this is to do the spline in linear/linear space... interpolation would
  // be better in log space, but then i'd have to take log(T) and exp(Rate)
  // for every call, and they are expensive...
  //
  for (int i=0; i<hr_Nspl; i++) {
    hr_t[i] = exp(log(10.0)*(1.0 +0.2*static_cast<double>(i)));
    hr_alpha[i] = caseb[i]/sqrt(hr_t[i]);
    hr_beta[i]  = coolb[i]/sqrt(hr_t[i]);
    hr_btot[i]  = CoolTot[i]/sqrt(hr_t[i]);
    hr_alpha2[i]= 0.0;
    hr_beta2[i] = 0.0;
    hr_btot2[i]  = 0.0;
  }
  interpolate.spline(hr_t, hr_alpha, hr_Nspl, 1.e99, 1.e99, hr_alpha2);
  interpolate.spline(hr_t, hr_beta,  hr_Nspl, 1.e99, 1.e99, hr_beta2 );
  interpolate.spline(hr_t, hr_btot,  hr_Nspl, 1.e99, 1.e99, hr_btot2 );

  MinTemp = hr_t[0];
  MaxTemp = hr_t[hr_Nspl-1];
  //
  // logarithmic slope for extrapolation of recombinatation rate alpha_b
  //
  MinSlope_alpha = (log10(hr_alpha[1])-log10(hr_alpha[0]))/
    (log10(hr_t[1])-log10(hr_t[0]));
  MaxSlope_alpha = (log10(hr_alpha[hr_Nspl-1])-log10(hr_alpha[hr_Nspl-2]))/
    (log10(hr_t[hr_Nspl-1])-log10(hr_t[hr_Nspl-2]));
  //cout << "\t\tAlpha min-slope="<<MinSlope_alpha<<" max-slope="<<MaxSlope_alpha<<"\n";
  //
  // logarthmic slope for extrapolation for beta-recomb
  //
  MinSlope_beta  = (log10(hr_beta[1]) -log10(hr_beta[0]))/
    (log10(hr_t[1])-log10(hr_t[0]));
  MaxSlope_beta  = (log10(hr_beta[hr_Nspl-1])-log10(hr_beta[hr_Nspl-2]))/
    (log10(hr_t[hr_Nspl-1])-log10(hr_t[hr_Nspl-2]));
  //cout << "\t\tBeta  min-slope="<<MinSlope_beta<<" max-slope="<<MaxSlope_beta<<"\n";
  //
  // logarthmic slope for extrapolation for beta-total
  //
  MinSlope_btot  = (log10(hr_btot[1]) -log10(hr_btot[0]))/
    (log10(hr_t[1])-log10(hr_t[0]));
  MaxSlope_btot  = (log10(hr_btot[hr_Nspl-1])-log10(hr_btot[hr_Nspl-2]))/
    (log10(hr_t[hr_Nspl-1])-log10(hr_t[hr_Nspl-2]));
  //cout << "\t\tB-tot min-slope="<<MinSlope_btot<<" max-slope="<<MaxSlope_btot<<"\n";

#ifdef TESTING
  ofstream outf("hummer_recomb.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"Hummer Recombination and Cooling Curve Data: Temperature(K) ";
  outf <<"Rate(cm^3/s) CoolB(erg.cm^3/s) Total-Cool[B+FF](erg.cm^3/s)\n";
  outf <<"# Cols: T  alpha_B*sqrt(T)  beta_B*sqrt(T)  beta_B^tot*sqrt(T)";
  outf <<"  alpha_B  beta_B  beta_B^tot.\n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< Hii_rad_recomb_rate(t)*sqrt(t) <<"\t";
    outf << Hii_rad_recomb_cooling(t)*sqrt(t) <<"\t";
    outf << Hii_total_cooling(t)*sqrt(t) <<"\t";
    outf << Hii_rad_recomb_rate(t) <<"\t";
    outf << Hii_rad_recomb_cooling(t) <<"\t";
    outf << Hii_total_cooling(t) <<"\n";
    t *=1.1;
  } while (t<1.0e9);
  outf.close();
#endif //TESTING

  return;
}

Hummer94_Hrecomb::~Hummer94_Hrecomb()
{
  hr_t      = mem.myfree(hr_t);
  hr_alpha  = mem.myfree(hr_alpha);
  hr_alpha2 = mem.myfree(hr_alpha2);
  hr_beta   = mem.myfree(hr_beta);
  hr_beta2  = mem.myfree(hr_beta2);
  hr_btot   = mem.myfree(hr_btot);
  hr_btot2  = mem.myfree(hr_btot2);
  return;
}


double Hummer94_Hrecomb::Hii_rad_recomb_rate(const double T)
{
#ifdef RT_TEST_PROBS
  return 2.59e-13; // set to const so can compare to analytic result
#endif //RT_TEST_PROBS

  //
  // Since adaptive integrators can overshoot, we need to be able to
  // silently deal with negative or infinite temperatures:
  //
  if (T<0.0 || !isfinite(T))
    return HUGE_VAL;

  //
  // Temperature, T, is in Kelvin.  We want to avoid running off the end of the 
  // spline fit, so we limit T to be between [10,1e7].  Outside these ranges we
  // extrapolate the rates based on the logarithmic slope of the function at the
  // relevant end of the range.
  //
  double rate=0.0;
  if      (T>MaxTemp) {
    rate = hr_alpha[hr_Nspl-1] *pow(T/MaxTemp, MaxSlope_alpha);
  }
  else if (T<MinTemp) {
    rate = hr_alpha[0] *pow(T/MinTemp, MinSlope_alpha);
  }
  else {
    interpolate.splint(hr_t, hr_alpha, hr_alpha2, hr_Nspl, T, &rate);
  }

  return rate;
}

double Hummer94_Hrecomb::Hii_rad_recomb_cooling(const double T)
{
  ///
  /// This is the Case B energy loss rate in erg*cm^3/s. Fitting function from
  /// Hummer, 1994, MNRAS, 268, 109.
  ///
  //
  // Since adaptive integrators can overshoot, we need to be able to
  // silently deal with negative or infinite temperatures:
  //
  if (T<0.0 || !isfinite(T))
    return HUGE_VAL;

  //
  // Temperature, T, is in Kelvin.  We want to avoid running off the end of the 
  // spline fit, so we limit T to be between [10,1e7].  Outside these ranges we
  // extrapolate the rates based on the logarithmic slope of the function at the
  // relevant end of the range.
  //
  double rate=0.0;
  if      (T>MaxTemp) {
    rate = hr_beta[hr_Nspl-1] *pow(T/MaxTemp, MaxSlope_beta);
  }
  else if (T<MinTemp) {
    rate = hr_beta[0] *pow(T/MinTemp, MinSlope_beta);
  }
  else {
    interpolate.splint(hr_t, hr_beta, hr_beta2, hr_Nspl, T, &rate);
  }

  return rate*kB*T;
}

double Hummer94_Hrecomb::Hii_total_cooling(const double T)
{
  ///
  /// This is the total energy loss rate in erg*cm^3/s (Case B + Free-free).
  /// Fitting function from Hummer, 1994, MNRAS, 268, 109.
  ///
  //
  // Since adaptive integrators can overshoot, we need to be able to
  // silently deal with negative or infinite temperatures:
  //
  if (T<0.0 || !isfinite(T))
    return HUGE_VAL;

  //
  // Temperature, T, is in Kelvin.  We want to avoid running off the end of the 
  // spline fit, so we limit T to be between [10,1e7].  Outside these ranges we
  // extrapolate the rates based on the logarithmic slope of the function at the
  // relevant end of the range.
  //
  double rate=0.0;
  if      (T>MaxTemp) {
    rate = hr_btot[hr_Nspl-1] *pow(T/MaxTemp, MaxSlope_btot);
  }
  else if (T<MinTemp) {
    rate = hr_btot[0] *pow(T/MinTemp, MinSlope_btot);
  }
  else {
    interpolate.splint(hr_t, hr_btot, hr_btot2, hr_Nspl, T, &rate);
  }
  //cout <<"TOTAL COOLING: T="<<T<<" rate="<<rate<<" kB="<<kB<<" T="<<T<<"\n";
  return rate*kB*T;
}

#ifdef TEST_HUMMER94_COOLING_FUNCTION
#include "cooling_SD93_cie.h"
// g++ -Wall -g Hummer94_Hrecomb.cc cooling_SD93_cie.cc ../global.cc ../cell_interface.cc -lreadline -o temp/testmp
int main()
{

  class Hummer94_Hrecomb rec;
  
  cout <<"rate="<< rec.Hii_rad_recomb_cooling(1.0e7);
  
  class cooling_function_SD93CIE c1;
  c1.setup_SD93_cie_MetalFree();
  class cooling_function_SD93CIE c2;
  c2.setup_SD93_cie_OnlyMetals();
  class cooling_function_SD93CIE c3;
  c3.setup_SD93_cie();

  return 0;
}
#endif // TEST_HUMMER94_COOLING_FUNCTION
