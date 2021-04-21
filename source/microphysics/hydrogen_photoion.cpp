/// \file hydrogen_photoion.cc
/// \author Jonathan Mackey
/// \date 2011.03.14
///
/// \description
/// This contains the class definitions for hydrogen_photoion, which calculates
/// the multifrequency photoionisation rate as a function of optical depth,
/// including spectral hardening due to preferential absorption of lower energy
/// photons.  The algorithm is the discretised photon-conserving formulation of
/// Mellema et al (2006,NewA,11,374).
///
/// Modifications:
/// - 2011.04.17 JM: Added discretised photoionisation rate for monochromatic
///    radiation source.
/// - 2011.04.18 JM: Fixed bugs.
/// - 2011.05.02 JM: Added MinCol, MaxCol to check for spline overruns on the
///    multifreq photoionisation and photoheating rates.
/// - 2011.05.04 JM: Added a discretised multifreq photoionisation rate with
///    an approximation for dtau<<1.  Fixed bugs, simplified code.
/// - 2011.06.20 JM: Got rid of non-ANSI-C exp10 functions
/// - 2011.07.03 JM: Added more digits to ln(10) constant.
/// - 2011.10.08 JM: Added switch to use interpolate.spline/splint instead of
/// the local
///    STL vector one, because the vector functions are slower by about 2.5%.
/// - 2011.11.01 JM: changed to using Tau instead of NH0 (hopefully more
///    efficient by allowing smaller spline array.
/// - 2012.12.26 JM: Added some hacks to study the photoionisation
///    cross section and Blackbody spectrum effects on results.
/// - 2014.03.27 JM: fixed bug in discrete monochromatic PI rate.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2019.09.04 JM: removed vector interpolation because not clear wehre
///    code comes from.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

//#define HACK_MODIFY_BB ///< scale the high energy BB emission
#define OSTERBROCK_XSEC  ///< use osterbrock photoionisation x-section.

#include "microphysics/hydrogen_photoion.h"

using namespace std;

#define LOGTEN 2.302585093

// ##################################################################
// ##################################################################

hydrogen_photoion::hydrogen_photoion()
{
  PI_Tau    = 0;
  PIrate    = 0;
  PIheat    = 0;
  LTPIrate  = 0;
  LTPIheat  = 0;
  PIrt_id   = 0;
  PIht_id   = 0;
  LTPIrt_id = 0;
  LTPIht_id = 0;
  return;
}

// ##################################################################
// ##################################################################

hydrogen_photoion::~hydrogen_photoion()
{
  if (PI_Tau) {
    PI_Tau   = mem.myfree(PI_Tau);
    PIrate   = mem.myfree(PIrate);
    PIheat   = mem.myfree(PIheat);
    LTPIrate = mem.myfree(LTPIrate);
    LTPIheat = mem.myfree(LTPIheat);
  }

  return;
}

// ##################################################################
// ##################################################################

// discretised multifrequency photoionisation rate
double hydrogen_photoion::Hi_discrete_multifreq_photoion_rate(
    const double
        Tau0,  ///< Optical depth of H0 (at 13.6eV) to front edge of cell.
    const double dTau0,  ///< Optical depth of H0 (at 13.6eV) through cell.
    const double nH,     ///< Local number density of H (per cm3).
    const double ds,     ///< path length through cell.
    const double Vshell  ///< Shell volume (cm3).
)
{
  //
  // C2Ray paper (Mellema et al. 2006,NewA,11,374).  equation 6.
  //
  // Logic of this function:
  // if (dtau<<1), use LowTau integration discretised.
  //
  // else use normal discretised integration.
  //
  double dtau = dTau0;
  double ans = 0.0, ans2 = 0.0;

  if (dtau < MinTau) {
    // splint returns value to be multiplied by dNH0/(n(H0)*Vshell)
    ans = max(MinTau, min(MaxTau, Tau0));
    interpolate.splint(PI_Tau, LTPIrate, LTPIrt_id, PI_Nspl, log10(ans), &ans);
    ans = exp(LOGTEN * ans) * dTau0
          / (Hi_monochromatic_photo_ion_xsection(JUST_IONISED) * nH * Vshell);
  }
  else if (dtau < 3.0 * MinTau) {
    // Weighted average of the optically thin and thick calculation.
    ans = max(MinTau, min(MaxTau, Tau0));
    interpolate.splint(PI_Tau, LTPIrate, LTPIrt_id, PI_Nspl, log10(ans), &ans);
    ans = exp(LOGTEN * ans) * dTau0
          / (Hi_monochromatic_photo_ion_xsection(JUST_IONISED) * nH * Vshell);
    ans2 = Hi_multifreq_photoionisation_rate(Tau0, nH, Vshell)
           - Hi_multifreq_photoionisation_rate(Tau0 + dTau0, nH, Vshell);
    ans = ans * 0.5 * (3.0 * MinTau - dtau) / MinTau
          + ans2 * 0.5 * (dtau - MinTau) / MinTau;
  }
  else {
    //
    // just take the difference between the two to get the rate
    //
    ans = Hi_multifreq_photoionisation_rate(Tau0, nH, Vshell)
          - Hi_multifreq_photoionisation_rate(Tau0 + dTau0, nH, Vshell);
  }

  // cout <<"VShell="<<Vshell<<" ds="<<ds<<" pir(near)=";
  // cout <<Hi_multifreq_photoionisation_rate(NH0,     nH,Vshell);
  // cout <<",
  // pir(far)="<<Hi_multifreq_photoionisation_rate(NH0+dNH0,nH,Vshell); cout
  // <<" PIR="<<ans<<"\n";

  return ans;
}

// ##################################################################
// ##################################################################

// discretised multifrequency photoionisation rate
double hydrogen_photoion::Hi_discrete_multifreq_photoheating_rate(
    const double
        Tau0,  ///< Optical depth of H0 (at 13.6eV) to front edge of cell.
    const double dTau0,  ///< Optical depth of H0 (at 13.6eV) through cell.
    const double nH,     ///< Local number density of H (per cm3).
    const double ds,     ///< path length through cell.
    const double Vshell  ///< Shell volume (cm3).
)
{
  //
  // see Hi_discrete_multifreq_photoion_rate() for details.
  //
  double dtau = dTau0;
  double ans = 0.0, ans2 = 0.0;

  if (dtau < MinTau) {
    // splint returns value to be multiplied by dNH0/(n(H0)*Vshell)
    ans = max(MinTau, min(MaxTau, Tau0));
    interpolate.splint(PI_Tau, LTPIheat, LTPIht_id, PI_Nspl, log10(ans), &ans);
    ans = exp(LOGTEN * ans) * dTau0
          / (Hi_monochromatic_photo_ion_xsection(JUST_IONISED) * nH * Vshell);
    //      Hi_multifreq_photoionisation_heating_rate(NH0,     nH,Vshell) -
    //      Hi_multifreq_photoionisation_heating_rate(NH0+dNH0,nH,Vshell)<<"\n";
  }
  else if (dtau < 3.0 * MinTau) {
    // Weighted average of the optically thin and thick calculation.
    ans = max(MinTau, min(MaxTau, Tau0));
    interpolate.splint(PI_Tau, LTPIheat, LTPIht_id, PI_Nspl, log10(ans), &ans);
    ans = exp(LOGTEN * ans) * dTau0
          / (Hi_monochromatic_photo_ion_xsection(JUST_IONISED) * nH * Vshell);
    ans2 =
        Hi_multifreq_photoionisation_heating_rate(Tau0, nH, Vshell)
        - Hi_multifreq_photoionisation_heating_rate(Tau0 + dTau0, nH, Vshell);
    ans = ans * 0.5 * (3.0 * MinTau - dtau) / MinTau
          + ans2 * 0.5 * (dtau - MinTau) / MinTau;
  }
  else {
    //
    // just take the difference between the two to get the rate
    //
    ans = Hi_multifreq_photoionisation_heating_rate(Tau0, nH, Vshell)
          - Hi_multifreq_photoionisation_heating_rate(Tau0 + dTau0, nH, Vshell);
  }
  return ans;
}

// ##################################################################
// ##################################################################

// Get multifrequency photoionisation rate for a given column density of H0,
// local number density of H0, and Shell volume ~4.Pi.[(R+)^3-(R-)^3]/3
// PHYSICALLY THIS IS MEANINGLESS, BECAUSE THERE IS NO POINT IN THE
// VSHELL PARAMETER WITHOUT HAVING A DTAU THROUGH THE SHELL. THIS FUNCTION
// IS A NUMERICAL CONVENIENCE.
//
double hydrogen_photoion::Hi_multifreq_photoionisation_rate(
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV).
    const double nH0,    ///< Local number density of H0 (per cm3).
    const double Vshell  ///< Shell volume (cm3).
)
{
  double ans = 0.0;
  //
  // splint gives back log10(PI-rate).  Check for running off the start or
  // end of the spline with Min/Max-Col values.  Assumes MinTau is small
  // enough that it is effectively zero, and that MaxTau is so large that it
  // gives a negligible photoionisation rate for any Vshell value.
  //
  ans = max(MinTau, min(MaxTau, Tau0));
  // ans = Tau0;
  interpolate.splint(PI_Tau, PIrate, PIrt_id, PI_Nspl, log10(ans), &ans);
  ans = exp(LOGTEN * ans) / (nH0 * Vshell);
  return ans;
}

// ##################################################################
// ##################################################################

//
// Get multifrequency photoionisation heating rate for a given column density
// of H0, local number density of H0, and Shell volume ~4.Pi.[(R+)^3-(R-)^3]/3
//
double hydrogen_photoion::Hi_multifreq_photoionisation_heating_rate(
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV).
    const double nH0,    ///< Local number density of H0 (per cm3).
    const double Vshell  ///< Shell volume (cm3).
)
{
  double ans = 0.0;
  //
  // splint gives back log10(PI-heating-rate).  Check for running off the
  // start or end of the spline with Min/Max-Col values.  Assumes MinTau is
  // small enough that it is effectively zero, and that MaxTau is so large
  // that it gives a negligible photoionisation rate for any Vshell value.
  //
  ans = max(MinTau, min(MaxTau, Tau0));
  interpolate.splint(PI_Tau, PIheat, PIht_id, PI_Nspl, log10(ans), &ans);
  ans = exp(LOGTEN * ans) / (nH0 * Vshell);
  return ans;
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::Hi_monochromatic_photo_ion_xsection_fractional(
    const double E)
{
  //
  // Sutherland & Dopita (2003,Textbook,eq.5.32) give this formula.
  // OSTERBROCK_XSEC replaces this with Osterbrock's (1989) eq. 2.31.
  //
  // N.B. 'E' is assumed to be in cgs units (ergs).
  // This function returns the ratio of the cross-section at energy E to the
  // cross-section at the ionisation edge E0=13.6eV

  if (E < 2.178720e-11)
    return 0.0;
  else {
#ifdef RT_TEST_PROBS
    return 1.0;
#else
#ifdef OSTERBROCK_XSEC
    return 1.34 * exp(-2.99 * log(E / 2.1788e-11))
           - 0.34 * exp(-3.99 * log(E / 2.1788e-11));
#else
    return exp(-3.5 * log(E / 2.18e-11));
#endif
#endif
  }
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::Hi_monochromatic_photo_ion_xsection(const double E)
{
  //
  // Sutherland & Dopita (2003,Textbook,eq.5.32) give this formula, which I
  // think is the same as that given in Osterbrock (1989).
  //
  // N.B. 'E' is assumed to be in cgs units (ergs).
  //
  if (E < 2.178720e-11)
    return 0.0;
  else {
#ifdef RT_TEST_PROBS
    return 6.3042e-18;
#else
#ifdef OSTERBROCK_XSEC
    return 6.3042e-18
           * (1.34 * exp(-2.99 * log(E / 2.1788e-11))
              - 0.34 * exp(-3.99 * log(E / 2.1788e-11)));
#else
    return 6.3042e-18 * exp(-3.5 * log(E / 2.18e-11));
#endif
#endif
  }
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::Hi_discrete_mono_photoion_rate(
    const double
        Tau0,  ///< Optical depth of H0 (at 13.6eV) to front edge of cell.
    const double dTau0,  ///< Optical depth of H0 (at 13.6eV) through cell.
    const double nH,     ///< Local number density of H (per cm3).
    const double Ndot,   ///< ionising photon luminosity of source
    const double E,      ///< Energy of photons (erg)
    const double ds,     ///< path length through cell.
    const double Vshell  ///< Shell volume (cm3).
)
{
  //
  // C2Ray paper (Mellema et al. 2006,NewA,11,374).  equation 6.
  // As for the multi-frequency case, this returns the PI rate per H
  // atom (so no need to multiply by the neutral fraction).  It does
  // this by dividing by the H number density instead of the H0
  // number density, to avoid multiplying and dividing by two small
  // numbers.
  //
  // First scale the cross-section from the ionsation edge to energy E:
  //
  double dtau = dTau0 * Hi_monochromatic_photo_ion_xsection_fractional(E);

  double rate = Ndot
                * exp(-Tau0 * Hi_monochromatic_photo_ion_xsection_fractional(E))
                / Vshell;
  //
  // if dtau<<1 then it is more numerically stable to approximate the
  // (1-exp(-dtau))/nH term by dtau/nH.  Otherwise evaluate the full
  // expression.
  //
  if (dtau < MinTau) {
    rate *= dtau / nH;  // Hi_monochromatic_photo_ion_xsection(E)*ds;
  }
  else {
    rate *= (1.0 - exp(-dtau)) / nH;
  }
  return rate;
}

// ##################################################################
// ##################################################################

void hydrogen_photoion::Setup_photoionisation_rate_table(
    const double Tstar,    ///< BB temperature (K)
    const double Rstar,    ///< Radius of star (cm)
    const double Lstar,    ///< Stellar luminosity (erg/s) (overrides Rstar!).
    const double Tau0min,  ///< Min Optical depth of H0 (at 13.6eV)
    const double Tau0max,  ///< Max Optical depth of H0 (at 13.6eV)
    const double Emax,     ///< Max energy to integrate to.
    const int Nsub,        ///< Number of sub-points in integration.
    const int Nspl         ///< Number of spline points.
)
{
  MinTau = Tau0min;
  MaxTau = Tau0max;
  //
  // First check that the total luminosity is close to the luminosity
  // calculated from 4Pi*R^2 sigma*T^4
  //
  // cout <<Tstar<<"  "<<Rstar<<"  "<<Lstar<<"\n";
  if (Tstar < 1000.0) rep.error("T<1000K.  Object is not a star!", Tstar);
  double L = pconst.StefanBoltzmannConst() * pow(Tstar, 4.0) * 4.0 * M_PI
             * Rstar * Rstar;
  if (fabs(1.0 - L / Lstar) > 0.05 && !pconst.equalD(Lstar, 1.0)) {
    cout << "\tLuminosities don't match! Lstar=" << Lstar;
    cout << " and Stefan-Bolt. law gives L=" << L;
    cout << "\t:Setup_photoionisation_rate_table: Ignoring Lstar and using "
            "4.Pi.R^2.sigma.T^4.\n";
    // cerr <<": NOTE Rstar given will be scaled to give the requested
    // Lstar\n";
    //
    // If L>Lstar, then we want to multiply the radius by sqrt(Lstar/L)<1
    // to decrease the radius and give the requested Luminosity.
    // Maybe we should just bug out?
    // rep.error("Lstar,Rstar,Tstar are inconsistent",Lstar/L);
    // L=sqrt(Lstar/L); // Radius can be multiplied by this below.
  }
  // else {
  //  L=1.0; // So that the radius is scaled by 1.0 below
  //}

#ifdef TESTING
  if (Nspl < 75) {
    cerr << "Setup_photoionisation_rate_table() using Nspl<75 is less "
            "accurate\n";
  }
  if (Nsub < 800) {
    cerr << "Setup_photoionisation_rate_table() using Nsub<800 is less "
            "accurate\n";
  }
#endif  // TESTING

  PI_Nspl = Nspl;

  if (PI_Tau) {
    PI_Tau   = mem.myfree(PI_Tau);
    PIrate   = mem.myfree(PIrate);
    PIheat   = mem.myfree(PIheat);
    LTPIrate = mem.myfree(LTPIrate);
    LTPIheat = mem.myfree(LTPIheat);
  }
  PI_Tau   = mem.myalloc(PI_Tau, Nspl);
  PIrate   = mem.myalloc(PIrate, Nspl);
  PIheat   = mem.myalloc(PIheat, Nspl);
  LTPIrate = mem.myalloc(LTPIrate, Nspl);
  LTPIheat = mem.myalloc(LTPIheat, Nspl);
  //
  // Calculate log10 of photoionisation and photoheating rates for Nspl
  // logarithmically spaced values of Tau0.
  //
  double lTmax = log10(MaxTau);
  double lTmin = log10(MinTau);
  double hh    = (lTmax - lTmin) / (Nspl - 1);
  double Tau0;
  for (int v = 0; v < Nspl; v++) {
    PI_Tau[v] = lTmin + v * hh;
    Tau0      = exp(LOGTEN * PI_Tau[v]);
    // cout <<"\t-- T*="<<Tstar<<", R*="<<Rstar<<", Tau0="<<Tau0<<",
    // Emax="<<Emax<<", Nsub="<<Nsub<<"\n";
    PIrate[v] =
        log10(photoion_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub));
    PIheat[v] = log10(
        photoheating_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub));
    //
    // LOW-DTAU APPROX INTEGRAL ---------------
    //
    LTPIrate[v] =
        log10(PI_LowTau_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub));
    LTPIheat[v] =
        log10(PH_LowTau_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub));
    // ----------------------------------------
    // cout <<"Tau0="<<Tau0<<", pir="<<exp(LOGTEN*PIrate[v])<<",
    // phr="<<exp(LOGTEN*PIheat[v])<<"\n";
  }

  // GSL:
  // Fit a cubic spline to the Photoionisation and Photoheating rates.
  // Large values in the 4th,5th args tell it to use natural boundary
  // conditions, which means set the second derivative to zero at the
  // endpoints. A small value (<1.0e30) indicates that this is the actual
  // value of the first derivative at the boundary values (4th is lower limit,
  // 5th is upper limit).
  // BOOST:
  // 4th and 5th args are the first derivatives at the endpoints.  As
  // far as I can tell the makima() method is continuous and
  // piecewise-constant 1st derivatives (i.e. not differentiable).
  //
  interpolate.spline(PI_Tau, PIrate, PI_Nspl, 0.0, 0.0, PIrt_id);
  interpolate.spline(PI_Tau, PIheat, PI_Nspl, 0.0, 0.0, PIht_id);

  // LOW-DTAU APPROX INTEGRAL ---------------
  interpolate.spline(PI_Tau, LTPIrate, PI_Nspl, 0.0, 0.0, LTPIrt_id);
  interpolate.spline(PI_Tau, LTPIheat, PI_Nspl, 0.0, 0.0, LTPIht_id);
  // ----------------------------------------

  return;
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::photoion_rate_source_integrand(
    const double E,      ///< ergs.
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0    ///< Optical depth of H0 (at 13.6eV)
)
{
  if (E < 2.18e-11) return 0.0;

  //
  // First the E^2 term;
  //
  double ans = E * E;
  //
  // Then exp(-tau), where tau = tau0*sigma(E)/sigma(E0).
  //
  ans *= exp(-Tau0 * Hi_monochromatic_photo_ion_xsection_fractional(E));
  //
  // Then 1/(exp(E/kT)-1)
  //
  ans /= (exp(E / (1.38e-16 * Tstar)) - 1.0);
  //#ifdef HACK_MODIFY_BB
  //  // multiply high energy emission by exp(-(0.03(E/eV))^5)
  //  ans *= exp(-pow(1.87e10*E,5));
  //#endif
  //
  // Then the prefactor: 8.pi^2.Rstar^2/(c^2.h^3)
  //
  ans *= 3.020e59 * Rstar * Rstar;

  return ans;
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::photoheating_rate_source_integrand(
    const double E,      ///< ergs.
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0    ///< Optical depth of H0 (at 13.6eV)
)
{
  return photoion_rate_source_integrand(E, Tstar, Rstar, Tau0) * (E - 2.18e-11);
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::photoion_rate_source_integral(
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV)
    const double Emax,   ///< Max energy to integrate to (ergs).
    const int Nsub       ///< Number of sub-points in integration.
)
{
  double Api = 0.0;

  //
  // Integration in log space: \int ( d(lnE) E*f(E) )
  // New variable is x=ln(E), so \int (dx exp(x)*f(exp(x)))
  //
  double Xmin = log(13.6 * 1.602e-12);  // ionisation threshold.
  double Xmax = log(Emax);
  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (Xmax - Xmin) / Nsub;
  //
  // f(0) lower limit
  //
  Api +=
      exp(Xmin) * photoion_rate_source_integrand(exp(Xmin), Tstar, Rstar, Tau0);
  //
  // f(N) upper limit
  //
  Api +=
      exp(Xmax) * photoion_rate_source_integrand(exp(Xmax), Tstar, Rstar, Tau0);
  //
  // Intermediate points.
  //
  int wt   = 4;
  double E = 0.0, X = 0.0;
  for (int i = 1; i < Nsub; i++) {
    X = Xmin + i * hh;
    E = exp(X);
    Api += wt * E * photoion_rate_source_integrand(E, Tstar, Rstar, Tau0);
    wt = 6 - wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Api *= hh / 3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Api < VERY_TINY_VALUE) {
    Api = VERY_TINY_VALUE;
  }
  return Api;
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::photoheating_rate_source_integral(
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV)
    const double Emax,   ///< Max energy to integrate to.
    const int Nsub       ///< Number of sub-points in integration.
)
{
  double Hpi = 0.0;

  //
  // Integration in log space: \int ( d(lnE) E*f(E) )
  // New variable is x=ln(E), so \int (dx exp(x)*f(exp(x)))
  //
  double Xmin = log(13.6 * 1.602e-12);  // ionisation threshold.
  double Xmax = log(Emax);
  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (Xmax - Xmin) / Nsub;
  //
  // f(0) lower limit
  //
  Hpi += exp(Xmin)
         * photoheating_rate_source_integrand(exp(Xmin), Tstar, Rstar, Tau0);
  //
  // f(N) upper limit
  //
  Hpi += exp(Xmax)
         * photoheating_rate_source_integrand(exp(Xmax), Tstar, Rstar, Tau0);
  //
  // Intermediate points.
  //
  int wt   = 4;
  double E = 0.0, X = 0.0;
  for (int i = 1; i < Nsub; i++) {
    X = Xmin + i * hh;
    E = exp(X);
    Hpi += wt * E * photoheating_rate_source_integrand(E, Tstar, Rstar, Tau0);
    wt = 6 - wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Hpi *= hh / 3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Hpi < VERY_TINY_VALUE) {
    Hpi = VERY_TINY_VALUE;
  }
  return Hpi;
}

// ##################################################################
// ##################################################################

//////////////////////////////////////////
// APPROXIMATE LOW-DTAU FUNCTIONS       //
//////////////////////////////////////////

double hydrogen_photoion::PI_LowTau_rate_source_integrand(
    const double E,      ///< ergs.
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0    ///< Optical depth of H0 (at 13.6eV)
)
{
  return photoion_rate_source_integrand(E, Tstar, Rstar, Tau0)
         * Hi_monochromatic_photo_ion_xsection(E);
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::PH_LowTau_rate_source_integrand(
    const double E,      ///< ergs.
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0    ///< Optical depth of H0 (at 13.6eV)
)
{
  return PI_LowTau_rate_source_integrand(E, Tstar, Rstar, Tau0)
         * (E - 2.18e-11);
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::PI_LowTau_rate_source_integral(
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV)
    const double Emax,   ///< Max energy to integrate to (ergs).
    const int Nsub       ///< Number of sub-points in integration.
)
{
  double Api = 0.0;

  //
  // Integration in log space: \int ( d(lnE) E*f(E) )
  // New variable is x=ln(E), so \int (dx exp(x)*f(exp(x)))
  //
  double Xmin = log(13.6 * 1.602e-12);  // ionisation threshold.
  double Xmax = log(Emax);
  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (Xmax - Xmin) / Nsub;
  //
  // f(0) lower limit
  //
  Api += exp(Xmin)
         * PI_LowTau_rate_source_integrand(exp(Xmin), Tstar, Rstar, Tau0);
  //
  // f(N) upper limit
  //
  Api += exp(Xmax)
         * PI_LowTau_rate_source_integrand(exp(Xmax), Tstar, Rstar, Tau0);
  //
  // Intermediate points.
  //
  int wt   = 4;
  double E = 0.0, X = 0.0;
  for (int i = 1; i < Nsub; i++) {
    X = Xmin + i * hh;
    E = exp(X);
    Api += wt * E * PI_LowTau_rate_source_integrand(E, Tstar, Rstar, Tau0);
    wt = 6 - wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Api *= hh / 3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Api < VERY_TINY_VALUE) {
    Api = VERY_TINY_VALUE;
  }
  return Api;
}

// ##################################################################
// ##################################################################

double hydrogen_photoion::PH_LowTau_rate_source_integral(
    const double Tstar,  ///< BB temperature (K)
    const double Rstar,  ///< Radius of star (cm)
    const double Tau0,   ///< Optical depth of H0 (at 13.6eV)
    const double Emax,   ///< Max energy to integrate to.
    const int Nsub       ///< Number of sub-points in integration.
)
{
  double Hpi = 0.0;

  //
  // Integration in log space: \int ( d(lnE) E*f(E) )
  // New variable is x=ln(E), so \int (dx exp(x)*f(exp(x)))
  //
  double Xmin = log(13.6 * 1.602e-12);  // ionisation threshold.
  double Xmax = log(Emax);
  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (Xmax - Xmin) / Nsub;
  //
  // f(0) lower limit
  //
  Hpi += exp(Xmin)
         * PH_LowTau_rate_source_integrand(exp(Xmin), Tstar, Rstar, Tau0);
  //
  // f(N) upper limit
  //
  Hpi += exp(Xmax)
         * PH_LowTau_rate_source_integrand(exp(Xmax), Tstar, Rstar, Tau0);
  //
  // Intermediate points.
  //
  int wt   = 4;
  double E = 0.0, X = 0.0;
  for (int i = 1; i < Nsub; i++) {
    X = Xmin + i * hh;
    E = exp(X);
    Hpi += wt * E * PH_LowTau_rate_source_integrand(E, Tstar, Rstar, Tau0);
    wt = 6 - wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Hpi *= hh / 3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Hpi < VERY_TINY_VALUE) {
    Hpi = VERY_TINY_VALUE;
  }
  return Hpi;
}

// ##################################################################
// ##################################################################

// TESTING CODE
/*

int main(int argc, char **argv)
{
  if (argc !=6) {
    cerr <<" Bad args...\n";
    exit(1);
  }
  //
  // Stellar temperature and radius:
  //
  double Tstar = atof(argv[1]);
  double Rstar = atof(argv[2])*6.96e10;
  int    Nsub  = atoi(argv[3]);
  double Emax  = atof(argv[4])*1.602e-12;
  int    Nspl  = atoi(argv[5]);

  //double c = 2.9979e10;
  //double h = 6.626e-27;
  //double kB= 1.381e-16;

  double Tau0=6.3042;
  double Api=0.0, Hpi=0.0;
  double One_Over_NH_Vsh = 1.0/(4*M_PI*3.086e18*3.086e18*3.086e16*100);
  class hydrogen_photoion MFQ;

  while (Tau0<8.2e20) {
    Api = MFQ.photoion_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub);
    Hpi = MFQ.photoheating_rate_source_integral(Tstar, Rstar, Tau0, Emax, Nsub);
    cout.setf( ios_base::scientific );
    cout.precision(6);
    //cout << NH0 <<"\t"<< Api <<"\t"<< Hpi <<"\t"<< Hpi/Api/1.602e-12 <<"\t";
    //cout << Api*One_Over_NH_Vsh <<"\t"<< Hpi*One_Over_NH_Vsh <<"\n";
    Tau0 *=3.16227766; // half a decade
    //NH0 *=1.584893192; // one fifth of a decade
  }

  MFQ.Setup_photoionisation_rate_table(Tstar,Rstar,1.0e39,1.0e-3,1.0e6,Emax,Nsub,Nspl);
  Api =
MFQ.Hi_multifreq_photoionisation_rate(630.42,100.0,4*M_PI*3.086e18*3.086e18*3.086e16);
  //cout <<"api="<<Api<<endl;

  Tau0=6.3042e-3;
  double nH0 = 100.0;
  double Vshell = 4.0*M_PI*3.086e18*3.086e18*3.086e16;
  while (Tau0<6.3042e5) {
    Api = MFQ.Hi_multifreq_photoionisation_rate(Tau0,nH0,Vshell);
    Hpi = MFQ.Hi_multifreq_photoionisation_heating_rate(Tau0,nH0,Vshell);
    cout.setf( ios_base::scientific );
    cout.precision(6);
    //cout << NH0 <<"\t"<< Api <<"\t"<< Hpi <<"\t"<< Hpi/Api/1.602e-12 <<"\t";
    //cout << MFQ.photoion_rate_source_integral(Tstar, Rstar, NH0, Emax,
Nsub)*One_Over_NH_Vsh <<"\t";
    //cout << MFQ.photoheating_rate_source_integral(Tstar, Rstar, NH0, Emax,
Nsub)*One_Over_NH_Vsh <<"\n";
    //NH0 *=3.16227766; // half a decade
    //NH0 *=1.584893192; // one fifth of a decade
    Tau0 *=1.1; // finer
  }

  class hydrogen_photoion MFQ1, MFQ2, MFQ3, MFQ4;
  MFQ1.Setup_photoionisation_rate_table(Tstar,Rstar,1.0e39,0.001,1e6,Emax,Nsub,Nspl);
  MFQ2.Setup_photoionisation_rate_table(Tstar,Rstar,1.0e39,0.001,1e6,Emax,2*Nsub,Nspl);
  MFQ3.Setup_photoionisation_rate_table(Tstar,Rstar,1.0e39,0.001,1e6,Emax,4*Nsub,Nspl);
  MFQ4.Setup_photoionisation_rate_table(Tstar,Rstar,1.0e39,0.0001,1e6,2*Emax,8*Nsub,Nspl);
  Tau0=1.0e-4;
  while (Tau0<6.3042e5) {
    cout.setf( ios_base::scientific );
    cout.precision(6);
    cout << Tau0/6.3042e-18 <<"\t";
    cout << MFQ1.Hi_multifreq_photoionisation_rate(Tau0,nH0,Vshell) <<"\t";
    cout << MFQ2.Hi_multifreq_photoionisation_rate(Tau0,nH0,Vshell) <<"\t";
    cout << MFQ3.Hi_multifreq_photoionisation_rate(Tau0,nH0,Vshell) <<"\t";
    cout << MFQ4.Hi_multifreq_photoionisation_rate(Tau0,nH0,Vshell) <<"\t";
    //
    // "exact" value:
    //
    cout << MFQ.photoion_rate_source_integral(Tstar, Rstar, Tau0, Emax,
32*Nsub)*One_Over_NH_Vsh <<"\t";

    //
    // Heating rates
    //
    cout << MFQ1.Hi_multifreq_photoionisation_heating_rate(Tau0,nH0,Vshell)
<<"\t"; cout << MFQ2.Hi_multifreq_photoionisation_heating_rate(Tau0,nH0,Vshell)
<<"\t"; cout << MFQ3.Hi_multifreq_photoionisation_heating_rate(Tau0,nH0,Vshell)
<<"\t"; cout << MFQ4.Hi_multifreq_photoionisation_heating_rate(Tau0,nH0,Vshell)
<<"\t";
    //
    // "exact" value:
    //
    cout << MFQ.photoheating_rate_source_integral(Tstar, Rstar, Tau0, Emax,
32*Nsub)*One_Over_NH_Vsh <<"\n";
     //NH0 *=3.16227766; // half a decade
    //NH0 *=1.584893192; // one fifth of a decade
    Tau0 *=1.1; // finer
  }



  return 0;
}

*/
