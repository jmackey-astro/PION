///
/// \file HHe_photoion.cpp
/// \author Jonathan Mackey
/// \date 2013.08.19
///
/// This file has function definitions for the HHe_photoion class,
/// which calculates multifrequency photoionisation rates of H, He,
/// and He+ as a function of the threshold optical depths of each of
/// these species.
///
/// Modifications:
/// - 2013.08.19 JM: written and tested.
/// - 2013.08.23 JM: Debugging.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifndef EXCLUDE_MPV9

//#define TEST_HHe_PION  // for testing this class in isolation.

#include "HHe_photoion.h"
//#include "grid/cell_interface.h"
#include "global.h"

#include <iostream>
using namespace std;

#define REGION_A 10
#define REGION_B 11
#define REGION_C 12

//#define VERY_TINY_VALUE 1.0e-200
//#define ln10() 2.302585093


// ##################################################################
// ##################################################################



HHe_photoion::HHe_photoion()
{
  //
  // cross sections in cm^2
  //
  sigmaHN  = 6.304e-18;  // Standard textbook.
  sigmaHeN = 7.56e-18;    // Marr \& West (1976,ADNDT,18,497);
  sigmaHeP = 1.576e-18;   // Osterbrock (1989)

  //
  // Min and Max optical depths to fit with spline, and the number of
  // points to use in the spline fit.
  //
  MinTau = 1.0e-4;
  MaxTau = 1.0e4;
  PI_Nspl = 50;

  //
  // threshold frequencies in units of 13.6eV
  //
  x0=1.00;
  x1=1.81;
  x2=4.00;
  x3=7.353; // 100 eV cutoff at high freq.

  prTau=0;

  pirA=0;
  pirA_2=0;
  pirB=0;
  pirB_2=0;
  pirC=0;
  pirC_2=0;
  phrA=0;
  phrA_2=0;
  phrB=0;
  phrB_2=0;
  phrC=0;
  phrC_2=0;

  iltA=0;
  iltA_2=0;
  iltB=0;
  iltB_2=0;
  iltC=0;
  iltC_2=0;
  hltA=0;
  hltA_2=0;
  hltB=0;
  hltB_2=0;
  hltC=0;
  hltC_2=0;

  return;
}



// ##################################################################
// ##################################################################



double HHe_photoion::xsection_th(
        const int ion ///< integer ID of ion.
        )
{
  if      (ion==ION_H_N)  return sigmaHN;
  else if (ion==ION_HE_N) return sigmaHeN;
  else if (ion==ION_HE_P) return sigmaHeP;
  else if (ion==ION_DUST) return 5.0e-35;
  else                    return -1.0e99;
}



// ##################################################################
// ##################################################################



HHe_photoion::~HHe_photoion()
{
  //
  // Free arrays.
  //
  if (prTau) {
    prTau   = mem.myfree(prTau);

    pirA    = mem.myfree(pirA);
    pirB    = mem.myfree(pirB);
    pirC    = mem.myfree(pirC);
    pirA_2    = mem.myfree(pirA_2);
    pirB_2    = mem.myfree(pirB_2);
    pirC_2    = mem.myfree(pirC_2);

    phrA    = mem.myfree(phrA);
    phrB    = mem.myfree(phrB);
    phrC    = mem.myfree(phrC);
    phrA_2    = mem.myfree(phrA_2);
    phrB_2    = mem.myfree(phrB_2);
    phrC_2    = mem.myfree(phrC_2);

    iltA    = mem.myfree(iltA);
    iltB    = mem.myfree(iltB);
    iltC    = mem.myfree(iltC);
    iltA_2    = mem.myfree(iltA_2);
    iltB_2    = mem.myfree(iltB_2);
    iltC_2    = mem.myfree(iltC_2);

    hltA    = mem.myfree(hltA);
    hltB    = mem.myfree(hltB);
    hltC    = mem.myfree(hltC);
    hltA_2    = mem.myfree(hltA_2);
    hltB_2    = mem.myfree(hltB_2);
    hltC_2    = mem.myfree(hltC_2);
  }
}



// ##################################################################
// ##################################################################



double HHe_photoion::tau_total(
        const int r,    ///< which region we are in.
        const double tau0, ///< optical depth to H0 at x=1.00
        const double tau1, ///< optical depth to He0 at x=1.81
        const double tau2  ///< optical depth to He+ at x=4.00
        )
{
  //
  // Frank & Mellema (1994,A\&A,289,937) optical depth, constant term
  // in front of the power law index.  This function just adds the 
  // contributions to Tau in the way that is specified for each
  // spectral window.
  //
  if      (r==REGION_A) {
    return tau0; //*exp(-2.8*log(x));
  }
  else if (r==REGION_B) {
    return (0.4559*tau0+2.7419*tau1); //*exp(-1.7*log(x));
  }
  else if (r==REGION_C) {
    return (tau0+16.6439*tau1+48.5029*tau2); //*exp(-2.8*log(x));
  }
  else  {
    cout <<"requested optical depth for non-ionising photon!\n";
    return -1.0e99;
  }

}



// ##################################################################
// ##################################################################



void HHe_photoion::tau_frac(
        const int r,    ///< which region we are in.
        const double tau0, ///< optical depth to H0 at x=1.00
        const double tau1, ///< optical depth to He0 at x=1.81
        const double tau2, ///< optical depth to He+ at x=4.00
        double frac[]      ///< fractions for each species.
        )
{
  //
  // Frank & Mellema (1994,A\&A,289,937) optical depth, constant term
  // in front of the power law index.  This function returns the
  // fractional opacity in each region for each species.
  //
  if      (r==REGION_A) {
    frac[0] = 1.0;
    frac[1] = 0.0;
    frac[2] = 0.0;
  }
  else if (r==REGION_B) {
    frac[0] = 0.4559*tau0/(0.4559*tau0 +2.7419*tau1);
    frac[1] = 1.0-frac[0];
    frac[2] = 0.0;
  }
  else if (r==REGION_C) {
    frac[2] = tau0 +16.6439*tau1 +48.5029*tau2;
    frac[0] = tau0/frac[2];
    frac[1] = 16.6439*tau1/frac[2];
    frac[2] = 1.0-frac[0]-frac[1];
  }
  else  {
    cout <<"requested optical depth for non-ionising photon!\n";
  }
  return;
}



// ##################################################################
// ##################################################################


double HHe_photoion::Ix_BB_integral_ABC(
        const int region,  ///< Spectral region we are in.
        const double tau0, ///< Tau: constant prefactor in Tau.
        const double T, ///< T, the source temperature (K)
        const double R  ///< R, the source radius (cm)
        )
{
  //
  // integrate from x=x0 to x=x1 with Simpson's rule.
  //
  int Nsub=1000;
  double slope=0.0, xmin=0.0, xmax=0.0;

  if      (region==REGION_A) {
    xmin = x0;
    xmax = x1;
    slope = -2.8;
  }
  else if (region==REGION_B) {
    xmin = x1;
    xmax = x2;
    slope = -1.7;
  }
  else if (region==REGION_C) {
    xmin = x2;
    xmax = x3; // 100 eV
    slope = -2.8;
  }
  else {
    cout <<"(HHe_photoion::Ix_BB_integral_ABC) ";
    cout <<"Bad spectral region requested: "<<region<<"\n";
    return -1.0e99;
  }

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (xmax-xmin)/Nsub;
  //
  // f(0) lower limit, f(N) upper limit
  //
  double Api = Ix_integrand(xmin,
                            tau0*exp(slope*log(xmin)),
                            Lum_BB(xmin,T,R));
  Api += Ix_integrand(xmax,
                      tau0*exp(slope*log(xmax)),
                      Lum_BB(xmax,T,R));
  //
  // Intermediate points.
  //
  int wt = 4;
  double x=0.0;
  for (int i=1; i<Nsub; i++) {
    x = xmin +i*hh;
    Api += wt*Ix_integrand(x,
                           tau0*exp(slope*log(x)),
                           Lum_BB(x,T,R));
    wt = 6-wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Api *= hh/3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Api < VERY_TINY_VALUE) {
    Api=VERY_TINY_VALUE;
  }
  return Api;
}



// ##################################################################
// ##################################################################



double HHe_photoion::Hx_BB_integral_ABC(
        const int region,  ///< Spectral region we are in.
        const double tau0, ///< Tau: constant prefactor in Tau.
        const double T, ///< T, the source temperature (K)
        const double R  ///< R, the source radius (cm)
        )
{
  //
  // integrate from x=x0 to x=x1 with Simpson's rule.
  //
  int Nsub=1000;
  double slope=0.0, xmin=0.0, xmax=0.0;

  if      (region==REGION_A) {
    xmin = x0;
    xmax = x1;
    slope = -2.8;
  }
  else if (region==REGION_B) {
    xmin = x1;
    xmax = x2;
    slope = -1.7;
  }
  else if (region==REGION_C) {
    xmin = x2;
    xmax = x3; // 100 eV
    slope = -2.8;
  }
  else {
    cout <<"(HHe_photoion::Hx_BB_integral_ABC) ";
    cout <<"Bad spectral region requested: "<<region<<"\n";
    return -1.0e99;
  }

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (xmax-xmin)/Nsub;
  //
  // f(0) lower limit, f(N) upper limit
  //
  double Hx = Hx_integrand(xmin,
                          tau0*exp(slope*log(xmin)),
                          Lum_BB(xmin,T,R));
  Hx += Hx_integrand( xmax,
                      tau0*exp(slope*log(xmax)),
                      Lum_BB(xmax,T,R));
  //
  // Intermediate points.
  //
  int wt = 4;
  double x=0.0;
  for (int i=1; i<Nsub; i++) {
    x = xmin +i*hh;
    Hx += wt*Hx_integrand(x,
                          tau0*exp(slope*log(x)),
                          Lum_BB(x,T,R));
    wt = 6-wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Hx *= hh/3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Hx < VERY_TINY_VALUE) {
    Hx=VERY_TINY_VALUE;
  }
  return Hx;
}



// ##################################################################
// ##################################################################



double HHe_photoion::ILTx_BB_integral_ABC(
        const int region,  ///< Spectral region we are in.
        const double tau0, ///< Tau: constant prefactor in Tau.
        const double T, ///< T, the source temperature (K)
        const double R  ///< R, the source radius (cm)
        )
{
  //
  // integrate from x=x0 to x=x1 with Simpson's rule.
  //
  int Nsub=1000;
  double slope=0.0, xmin=0.0, xmax=0.0;

  if      (region==REGION_A) {
    xmin = x0;
    xmax = x1;
    slope = -2.8;
  }
  else if (region==REGION_B) {
    xmin = x1;
    xmax = x2;
    slope = -1.7;
  }
  else if (region==REGION_C) {
    xmin = x2;
    xmax = x3; // 100 eV
    slope = -2.8;
  }
  else {
    cout <<"(HHe_photoion::Ix_BB_integral_ABC) ";
    cout <<"Bad spectral region requested: "<<region<<"\n";
    return -1.0e99;
  }

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (xmax-xmin)/Nsub;
  //
  // f(0) lower limit, f(N) upper limit
  //
  double Api = exp(slope*log(xmin))*Ix_integrand(xmin,
                            tau0*exp(slope*log(xmin)),
                            Lum_BB(xmin,T,R));
  Api += exp(slope*log(xmax))*Ix_integrand(xmax,
                      tau0*exp(slope*log(xmax)),
                      Lum_BB(xmax,T,R));
  //
  // Intermediate points.
  //
  int wt = 4;
  double x=0.0;
  for (int i=1; i<Nsub; i++) {
    x = xmin +i*hh;
    Api += wt*exp(slope*log(x))*Ix_integrand(x,
                           tau0*exp(slope*log(x)),
                           Lum_BB(x,T,R));
    wt = 6-wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Api *= hh/3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Api < VERY_TINY_VALUE) {
    Api=VERY_TINY_VALUE;
  }
  return Api;
}



// ##################################################################
// ##################################################################



double HHe_photoion::HLTx_BB_integral_ABC(
        const int region,  ///< Spectral region we are in.
        const double tau0, ///< Tau: constant prefactor in Tau.
        const double T, ///< T, the source temperature (K)
        const double R  ///< R, the source radius (cm)
        )
{
  //
  // integrate from x=x0 to x=x1 with Simpson's rule.
  //
  int Nsub=1000;
  double slope=0.0, xmin=0.0, xmax=0.0;

  if      (region==REGION_A) {
    xmin = x0;
    xmax = x1;
    slope = -2.8;
  }
  else if (region==REGION_B) {
    xmin = x1;
    xmax = x2;
    slope = -1.7;
  }
  else if (region==REGION_C) {
    xmin = x2;
    xmax = x3; // 100 eV
    slope = -2.8;
  }
  else {
    cout <<"(HHe_photoion::Hx_BB_integral_ABC) ";
    cout <<"Bad spectral region requested: "<<region<<"\n";
    return -1.0e99;
  }

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (xmax-xmin)/Nsub;
  //
  // f(0) lower limit, f(N) upper limit
  //
  double Hx = exp(slope*log(xmin))*Hx_integrand(xmin,
                          tau0*exp(slope*log(xmin)),
                          Lum_BB(xmin,T,R));
  Hx += exp(slope*log(xmax))*Hx_integrand( xmax,
                      tau0*exp(slope*log(xmax)),
                      Lum_BB(xmax,T,R));
  //
  // Intermediate points.
  //
  int wt = 4;
  double x=0.0;
  for (int i=1; i<Nsub; i++) {
    x = xmin +i*hh;
    Hx += wt*exp(slope*log(x))*Hx_integrand(x,
                          tau0*exp(slope*log(x)),
                          Lum_BB(x,T,R));
    wt = 6-wt;
  }
  //
  // finally multiply by hh/3.0
  //
  Hx *= hh/3.0;

  //
  // Check for zero or very small values, because we have to take the log10
  // of this value in the spline fit.  So we set a tiny, but finite, value
  // as the floor value.
  //
  if (Hx < VERY_TINY_VALUE) {
    Hx=VERY_TINY_VALUE;
  }
  return Hx;
}


// ##################################################################
// ##################################################################


void HHe_photoion::Setup_photoionisation_rate_table(
        const double Teff,  ///< BB temperature (K)
        const double Rstar, ///< Radius of star (cm)
        const double, ///< UNUSED.
        const double, ///< Min optical depth at 13.6eV (UNUSED)
        const double, ///< Max optical depth at 13.6eV (UNUSED)
        const double, ///< UNUSED.
        const int,    ///< UNUSED.
        const int     ///< UNUSED.
        )
{
  // ----------------------------------------
  //
  // Setup arrays: pir for photoionisation, phr for photoheating,
  // ilt for low-dtau photoionisation, hlt for low-dtua photoheating.
  //
  if (prTau) {
    prTau   = mem.myfree(prTau);

    pirA    = mem.myfree(pirA);
    pirB    = mem.myfree(pirB);
    pirC    = mem.myfree(pirC);
    pirA_2    = mem.myfree(pirA_2);
    pirB_2    = mem.myfree(pirB_2);
    pirC_2    = mem.myfree(pirC_2);

    phrA    = mem.myfree(phrA);
    phrB    = mem.myfree(phrB);
    phrC    = mem.myfree(phrC);
    phrA_2    = mem.myfree(phrA_2);
    phrB_2    = mem.myfree(phrB_2);
    phrC_2    = mem.myfree(phrC_2);

    iltA    = mem.myfree(iltA);
    iltB    = mem.myfree(iltB);
    iltC    = mem.myfree(iltC);
    iltA_2    = mem.myfree(iltA_2);
    iltB_2    = mem.myfree(iltB_2);
    iltC_2    = mem.myfree(iltC_2);

    hltA    = mem.myfree(hltA);
    hltB    = mem.myfree(hltB);
    hltC    = mem.myfree(hltC);
    hltA_2    = mem.myfree(hltA_2);
    hltB_2    = mem.myfree(hltB_2);
    hltC_2    = mem.myfree(hltC_2);
  }
  prTau   = mem.myalloc(prTau,  PI_Nspl);

  pirA    = mem.myalloc(pirA,  PI_Nspl);
  pirB    = mem.myalloc(pirB,  PI_Nspl);
  pirC    = mem.myalloc(pirC,  PI_Nspl);
  pirA_2    = mem.myalloc(pirA_2, PI_Nspl);
  pirB_2    = mem.myalloc(pirB_2, PI_Nspl);
  pirC_2    = mem.myalloc(pirC_2, PI_Nspl);

  phrA    = mem.myalloc(phrA,  PI_Nspl);
  phrB    = mem.myalloc(phrB,  PI_Nspl);
  phrC    = mem.myalloc(phrC,  PI_Nspl);
  phrA_2    = mem.myalloc(phrA_2, PI_Nspl);
  phrB_2    = mem.myalloc(phrB_2, PI_Nspl);
  phrC_2    = mem.myalloc(phrC_2, PI_Nspl);

  iltA    = mem.myalloc(iltA,  PI_Nspl);
  iltB    = mem.myalloc(iltB,  PI_Nspl);
  iltC    = mem.myalloc(iltC,  PI_Nspl);
  iltA_2    = mem.myalloc(iltA_2, PI_Nspl);
  iltB_2    = mem.myalloc(iltB_2, PI_Nspl);
  iltC_2    = mem.myalloc(iltC_2, PI_Nspl);

  hltA    = mem.myalloc(hltA,  PI_Nspl);
  hltB    = mem.myalloc(hltB,  PI_Nspl);
  hltC    = mem.myalloc(hltC,  PI_Nspl);
  hltA_2    = mem.myalloc(hltA_2, PI_Nspl);
  hltB_2    = mem.myalloc(hltB_2, PI_Nspl);
  hltC_2    = mem.myalloc(hltC_2, PI_Nspl);

  // ----------------------------------------

  //
  // Calculate log10 of photoionisation and photoheating rates for PI_Nspl
  // logarithmically spaced values of Tau0.
  //
  double lTmax = log10(MaxTau);
  double lTmin = log10(MinTau);
  double hh = (lTmax-lTmin)/(PI_Nspl-1);
  double Tau0;
  for (int v=0; v<PI_Nspl; v++) {

    prTau[v] = lTmin +v*hh;
    Tau0 = exp(ln10()*prTau[v]);
    //cout <<"\t-- T*="<<Tstar<<", R*="<<Rstar<<", Tau0="<<Tau0<<", Emax="<<Emax<<", Nsub="<<Nsub<<"\n";
    //
    // First the normal integrals:
    //
    pirA[v] = log10(Ix_BB_integral_ABC(REGION_A,Tau0,Teff,Rstar));
    pirB[v] = log10(Ix_BB_integral_ABC(REGION_B,Tau0,Teff,Rstar));
    pirC[v] = log10(Ix_BB_integral_ABC(REGION_C,Tau0,Teff,Rstar));
    pirA_2[v] = 0.0;
    pirB_2[v] = 0.0;
    pirC_2[v] = 0.0;
    phrA[v] = log10(Hx_BB_integral_ABC(REGION_A,Tau0,Teff,Rstar));
    phrB[v] = log10(Hx_BB_integral_ABC(REGION_B,Tau0,Teff,Rstar));
    phrC[v] = log10(Hx_BB_integral_ABC(REGION_C,Tau0,Teff,Rstar));
    phrA_2[v] = 0.0;
    phrB_2[v] = 0.0;
    phrC_2[v] = 0.0;
    //
    // Now the low-dtau approximate integrals:
    //
    iltA[v] = log10(ILTx_BB_integral_ABC(REGION_A,Tau0,Teff,Rstar));
    iltB[v] = log10(ILTx_BB_integral_ABC(REGION_B,Tau0,Teff,Rstar));
    iltC[v] = log10(ILTx_BB_integral_ABC(REGION_C,Tau0,Teff,Rstar));
    iltA_2[v] = 0.0;
    iltB_2[v] = 0.0;
    iltC_2[v] = 0.0;
    hltA[v] = log10(HLTx_BB_integral_ABC(REGION_A,Tau0,Teff,Rstar));
    hltB[v] = log10(HLTx_BB_integral_ABC(REGION_B,Tau0,Teff,Rstar));
    hltC[v] = log10(HLTx_BB_integral_ABC(REGION_C,Tau0,Teff,Rstar));
    hltA_2[v] = 0.0;
    hltB_2[v] = 0.0;
    hltC_2[v] = 0.0;
  }

  // ----------------------------------------

  //
  // Fit a cubic spline to the Photoionisation and Photoheating rates.
  // Large values in the 4th,5th args tell it to use natural boundary conditions,
  // which means set the second derivative to zero at the endpoints.
  // A small value (<1.0e30) indicates that this is the actual value of the first
  // derivative at the boundary values (4th is lower limit, 5th is upper limit).
  //
  GS.spline(prTau, pirA, PI_Nspl, 1.e99, 1.e99, pirA_2);
  GS.spline(prTau, pirB, PI_Nspl, 1.e99, 1.e99, pirB_2);
  GS.spline(prTau, pirC, PI_Nspl, 1.e99, 1.e99, pirC_2);

  GS.spline(prTau, phrA, PI_Nspl, 1.e99, 1.e99, phrA_2);
  GS.spline(prTau, phrB, PI_Nspl, 1.e99, 1.e99, phrB_2);
  GS.spline(prTau, phrC, PI_Nspl, 1.e99, 1.e99, phrC_2);

  GS.spline(prTau, iltA, PI_Nspl, 1.e99, 1.e99, iltA_2);
  GS.spline(prTau, iltB, PI_Nspl, 1.e99, 1.e99, iltB_2);
  GS.spline(prTau, iltC, PI_Nspl, 1.e99, 1.e99, iltC_2);

  GS.spline(prTau, hltA, PI_Nspl, 1.e99, 1.e99, hltA_2);
  GS.spline(prTau, hltB, PI_Nspl, 1.e99, 1.e99, hltB_2);
  GS.spline(prTau, hltC, PI_Nspl, 1.e99, 1.e99, hltC_2);
  // ----------------------------------------

  return;
}



// ##################################################################
// ##################################################################


int HHe_photoion::HHe_photoion_rate(
        const double Tau0,  ///< Tau(H0) to front edge of cell
        const double dTau0, ///< Delta-Tau(H0) through cell
        const double Tau1,  ///< Tau(He0) to front edge of cell
        const double dTau1, ///< Delta-Tau(He0) through cell
        const double Tau2,  ///< Tau(He+) to front edge of cell
        const double dTau2, ///< Delta-Tau(He+) through cell
        const double nH,    ///< Local number density of H (per cm3) n(H).
        const double Vsh,   ///< Shell volume (cm3).
        double pir[],  ///< ionisation rates (per H atom)
        double phr[]   ///< heating rates (per H atom)
        )
{

  double t1 = tau_total(REGION_C,Tau0, Tau1, Tau2);
  double t2 = tau_total(REGION_C,Tau0+dTau0, Tau1+dTau1, Tau2+dTau2);
  get_region_HHe_integral_diff(REGION_C,t1,t2, &pir[2], &phr[2]);

  t1 = tau_total(REGION_B,Tau0, Tau1, Tau2);
  t2 = tau_total(REGION_B,Tau0+dTau0, Tau1+dTau1, Tau2+dTau2);
  get_region_HHe_integral_diff(REGION_B,t1,t2, &pir[1], &phr[1]);

  t1 = tau_total(REGION_A,Tau0, Tau1, Tau2);
  t2 = tau_total(REGION_A,Tau0+dTau0, Tau1+dTau1, Tau2+dTau2);
  get_region_HHe_integral_diff(REGION_A,t1,t2, &pir[0], &phr[0]);

  // ----------------------------------------------------------------
  //
  // Now we add the RegionC rates to H0 and He0 rates, and add the 
  // RegionB rates to H0 rates.  But we need to divide the photons
  // between the three species in each region, with the tau_frac()
  // function.
  //
  // ----------------------------------------------------------------

  double frac[3];
  //
  // region B
  //
  tau_frac(REGION_B,Tau0, Tau1, Tau2, frac);

  pir[0] += frac[0]*pir[1];  // give to H0
  pir[1] *= frac[1];         // take from He0

  phr[0] += frac[0]*phr[1];
  phr[1] *= frac[1];


  //
  // region C
  //
  tau_frac(REGION_C,Tau0, Tau1, Tau2, frac);

  pir[0] += frac[0]*pir[2];  // give to H0
  pir[1] += frac[1]*pir[2];  // give to He0
  pir[2] *= frac[2];         // take from He+

  phr[0] += frac[0]*phr[2];
  phr[1] *= frac[1]*phr[2];
  phr[2] *= frac[2];


  //
  // Now divide everything by n(H)*V(shell)
  //
  pir[0] /= nH*Vsh;
  pir[1] /= nH*Vsh;
  pir[2] /= nH*Vsh;

  phr[0] /= nH*Vsh;
  phr[1] /= nH*Vsh;
  phr[2] /= nH*Vsh;

  pir[0] = std::max(pir[0],VERY_TINY_VALUE);
  pir[1] = std::max(pir[1],VERY_TINY_VALUE);
  pir[2] = std::max(pir[2],VERY_TINY_VALUE);

  phr[0] = std::max(phr[0],VERY_TINY_VALUE);
  phr[1] = std::max(phr[1],VERY_TINY_VALUE);
  phr[2] = std::max(phr[2],VERY_TINY_VALUE);

#ifdef MP9_TESTING
//  cout <<"  pir0 = "<<*pir0;
//  cout <<"  pir1 = "<<*pir1;
//  cout <<"  pir2 = "<<*pir2;
//  cout <<"  phr0 = "<<*phr0;
//  cout <<"  phr1 = "<<*phr1;
//  cout <<"  phr2 = "<<*phr2;
//  cout <<"\n";
#endif // MP9_TESTING

  return 0;
}



// ##################################################################
// ##################################################################



void HHe_photoion::get_region_HHe_integral_diff(
        const int region, ///< Which spectral region
        const double t1,  ///< tau_0
        const double t2,  ///< tau_0+dtau_0
        double *pir,      ///< photoionisation rate 
        double *phr       ///< photoheating rate
        )
{
  double ans=0.0;

  //
  // set arrays to look-up based on which region was requested.
  //
  double *ia1=0, *ia2=0, *ha1=0, *ha2=0,
         *il1=0, *il2=0, *hl1=0, *hl2=0;
  switch (region) {

    case REGION_A:
    ia1 = pirA;
    ia2 = pirA_2;
    ha1 = phrA;
    ha2 = phrA_2;
    il1 = iltA;
    il2 = iltA_2;
    hl1 = hltA;
    hl2 = hltA_2;
    break;

    case REGION_B:
    ia1 = pirB;
    ia2 = pirB_2;
    ha1 = phrB;
    ha2 = phrB_2;
    il1 = iltB;
    il2 = iltB_2;
    hl1 = hltB;
    hl2 = hltB_2;
    break;

    case REGION_C:
    ia1 = pirC;
    ia2 = pirC_2;
    ha1 = phrC;
    ha2 = phrC_2;
    il1 = iltC;
    il2 = iltC_2;
    hl1 = hltC;
    hl2 = hltC_2;
    break;

    default:
    *pir=-1.0e99;
    *phr=-1.0e99;
    return;
    break;
  }

  //
  // If we have dTau<<1, then we use the low-dtau integrals because
  // the normal integrals will subtract to give identically zero.
  //
  if      (t2-t1<MinTau) {
    ans = max(MinTau, min(MaxTau, t1));
    GS.splint(prTau, il1, il2, PI_Nspl, log10(ans), &ans);
    *pir = (t2-t1)*exp(ln10()*ans);
    
    ans = max(MinTau, min(MaxTau, t1));
    GS.splint(prTau, hl1, hl2, PI_Nspl, log10(ans), &ans);
    *phr = (t2-t1)*exp(ln10()*ans);
  }
  //
  // Otherwise, if t1>MaxTau, then set the rates to zero.
  //
  else if (t1>MaxTau) {
    *pir = 0.0;
    *phr = 0.0;
  }
  //
  // Otherwise, if t2>MaxTau, but t1<MaxTau, then the second term is
  // zero (no photons leave cell, all are absorbed).  But we don't
  // need to do anything about this -- the normal case will deal with
  // this automatically.
  //
  //else if (t2>MaxTau) {
  //}
  //
  // Otherwise, we look up the normal interpolation tables for the
  // pre-calculated integrals.
  //
  else {
    // -------------- PHOTOIONISATION RATES -----------------
    //
    // Add the term for Tau.
    //
    ans = max(MinTau, min(MaxTau, t1));
    GS.splint(prTau, ia1, ia2, PI_Nspl, log10(ans), &ans);
    *pir = exp(ln10()*ans);
#ifdef MP9_TESTING
    cout.precision(10);
    double temp=exp(ln10()*ans);
    cout <<"a = "<<exp(ln10()*ans);
#endif  // MP9_TESTING
    //
    // Subtract the term for Tau+DTau.
    //
    ans = max(MinTau, min(MaxTau, t2));
    GS.splint(prTau, ia1, ia2, PI_Nspl, log10(ans), &ans);
    *pir -= exp(ln10()*ans);
#ifdef MP9_TESTING
    cout <<", b = "<<exp(ln10()*ans)<<" rel.diff="<< (temp-exp(ln10()*ans))/temp<<"\n";
#endif  // MP9_TESTING
    // -------------- PHOTOIONISATION RATES -----------------

    // -------------- PHOTOHEATING    RATES -----------------
    //
    // Add the term for Tau.
    //
    ans = max(MinTau, min(MaxTau, t1));
    GS.splint(prTau, ha1, ha2, PI_Nspl, log10(ans), &ans);
    *phr = exp(ln10()*ans);
#ifdef MP9_TESTING
    temp=exp(ln10()*ans);
    cout <<"a = "<<exp(ln10()*ans);
#endif  // MP9_TESTING
    //
    // Subtract the term for Tau+DTau.
    //
    ans = max(MinTau, min(MaxTau, t2));
    GS.splint(prTau, ha1, ha2, PI_Nspl, log10(ans), &ans);
    *phr -= exp(ln10()*ans);
#ifdef MP9_TESTING
    cout <<", b = "<<exp(ln10()*ans)<<" rel.diff="<< (temp-exp(ln10()*ans))/temp<<"\n";
#endif  // MP9_TESTING
    // -------------- PHOTOHEATING    RATES -----------------
  }

  return;
}


// ##################################################################
// ##################################################################

#endif //  if not EXCLUDE_MPV9



#ifdef TEST_HHe_PION
int main (int argc, char **argv)
{

  class HHe_photoion H;

  double Teff  = atof(argv[1]);    // typical of O star.
  double Rstar = atof(argv[2])*6.96e10;  // 10 solar radii

  //double Tau=1.0e-10;
  //double hh=1.0;
  //for (int i=0; i<5000; i++) {
  //  Tau *= 10;
  //  cout <<Tau;
  //  cout <<"  "<< H.Ix_BB_integral_ABC(REGION_A,Tau,Teff,Rstar);
  //  cout <<"  "<< H.Ix_BB_integral_ABC(REGION_B,Tau,Teff,Rstar);
  //  cout <<"  "<< H.Ix_BB_integral_ABC(REGION_C,Tau,Teff,Rstar);
  //  cout <<"  "<< H.Hx_BB_integral_ABC(REGION_A,Tau,Teff,Rstar);
  //  cout <<"  "<< H.Hx_BB_integral_ABC(REGION_B,Tau,Teff,Rstar);
  //  cout <<"  "<< H.Hx_BB_integral_ABC(REGION_C,Tau,Teff,Rstar);
  //  cout <<"\n";
  //  if (Tau>1.0e6) break;
  //}
  
  H.Setup_photoionisation_rate_table(Teff,  Rstar, -1.0, -1.0, -1.0, -1.0, -1.0, -1);

  double r  = 3.086e18;
  double dr = 1.543e17;
  double Vshell = 4.0*M_PI/3.0*(pow(r+dr,3.0)-pow(r,3.0));
  double nH = 1.0;
  double i0,i1,i2,h0,h1,h2;
  //double tau  = atof(argv[3]);
  //double dtau = atof(argv[4]);
  //H.HHe_photoion_rate(tau,dtau, tau,dtau, tau,dtau, nH, Vshell, &i0, &i1, &i2, &h0, &h1, &h2);
  
  //
  // Loop over all tau, dtau variables to get a big table so I can
  // plot it and look for any problem areas of parameter space.
  //
  size_t Ntau  = 64;
  size_t Ndtau = 64;
  double lTauMin=-6.0, lTauMax=4.0;
  double t0=1.0, t1=1.0, t2=1.0, dt0=0.01, dt1=0.01, dt2=0.01;
  cout.precision(8);
  for (size_t i=0; i<Ntau; i++) {
    t0 = exp(ln10()*(lTauMin+ i*(lTauMax-lTauMin)/Ntau));
    //for (size_t j=0; j<Ntau; j++) {
    //  t1 = exp(ln10()*(lTauMin+ j*(lTauMax-lTauMin)/Ntau));
      //for (size_t k=0; k<Ntau; k++) {
      //  t2 = exp(ln10()*(lTauMin+ k*(lTauMax-lTauMin)/Ntau));
        for (size_t di=0; di<Ndtau; di++) {
          dt0 = exp(ln10()*(lTauMin+ di*(lTauMax-lTauMin)/Ndtau));
          for (size_t dj=0; dj<Ndtau; dj++) {
            dt1 = exp(ln10()*(lTauMin+ dj*(lTauMax-lTauMin)/Ndtau));
            //for (size_t dk=0; dk<Ndtau; dk++) {
            //  dt2 = exp(ln10()*(lTauMin+ dk*(lTauMax-lTauMin)/Ndtau));

              H.HHe_photoion_rate(t0,dt0, t1,dt1, t2,dt2, nH, Vshell, &i0, &i1, &i2, &h0, &h1, &h2);
              /*
              cout << t0;
              cout <<" "<<t1;
              cout <<" "<<t2;
              cout <<"  "<<dt0;
              cout <<" "<<dt1;
              cout <<" "<<dt2;
              cout <<"  "<<i0;
              cout <<" "<<i1;
              cout <<" "<<i2;
              cout <<"  "<<h0;
              cout <<" "<<h1;
              cout <<" "<<h2;
              cout <<"\n";
              */
            //}
          }
        }
        //cout <<"\n";
      //}
    //}
  }
      

  return 0;
}
    
#endif // TEST_HHe_PION

