///
/// \file HHe_photoion.h
/// \author Jonathan Mackey
/// \date 2013.08.19
///
/// This file has the class declaration for the HHe_photoion class,
/// which calculates multifrequency photoionisation rates of H, He,
/// and He+ as a function of the threshold optical depths of each of
/// these species.
///
/// Modifications:
/// - 2013.08.19 JM: written and tested.
///

#ifndef HHE_PHOTOION_H
#define HHE_PHOTOION_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifndef EXCLUDE_MPV9

#ifndef INTEL
#include <cmath>     // Header file from gcc
#else
#include <mathimf.h> // Header file from Intel Compiler
#endif

#include "constants.h"

class HHe_photoion :
  virtual public constants
{
  public:
  HHe_photoion();
  ~HHe_photoion();

  ///
  /// Threshold cross-section for ionisation of an ion (cm^2).
  ///
  double xsection_th(
        const int  ///< integer ID of ion.
        );

  ///
  /// For given stellar parameters, set up the photoionisation and
  /// photoheating spline interpolation tables for the rates as a
  /// function of optical depth.
  ///
  void Setup_photoionisation_rate_table(
        const double, ///< BB temperature (K)
        const double, ///< Radius of star (cm)
        const double, ///< UNUSED.
        const double, ///< UNUSED.
        const double, ///< UNUSED.
        const double, ///< UNUSED.
        const int,    ///< UNUSED.
        const int     ///< UNUSED.
        );

  int HHe_photoion_rate(
        const double, ///< Tau(H0) to front edge of cell
        const double, ///< Delta-Tau(H0) through cell
        const double, ///< Tau(He0) to front edge of cell
        const double, ///< Delta-Tau(He0) through cell
        const double, ///< Tau(He+) to front edge of cell
        const double, ///< Delta-Tau(He+) through cell
        const double, ///< Local number density of H (per cm3) n(H).
        const double, ///< Shell volume (cm3).
        double *,     ///< H0, He0, He+ ionisation rates (per H atom)
        double *      ///< H0, He0, He+ heating rates (per H atom)
        );

  protected:
  ///
  /// For each region, what fraction of the photons are absorbed
  /// by each species, based on the fraction of the optical depth
  /// in the region that is from each species.  Friedrich+,2012.
  ///
  void tau_frac(
        const int,    ///< which region we are in.
        const double, ///< optical depth to H0 at x=1.00
        const double, ///< optical depth to He0 at x=1.81
        const double, ///< optical depth to He+ at x=4.00
        double []     ///< fractions for each species.
        );

  private:

  double sigmaHN;  ///< H cross section at threshold frequency
  double sigmaHeN; ///< He0 cross section at threshold frequency
  double sigmaHeP; ///< He+ cross section at threshold frequency

  double
    x0, ///< threshold frequency of H in units of nu_th
    x1, ///< threshold frequency of He in units of nu_th
    x2, ///< threshold frequency of He+ in units of nu_th
    x3; ///< Max. frequency of source spectrum (100eV)

  int PI_Nspl;   ///< number of elements in spline integral.
  double MinTau, ///< Smallest optical depth in splint (tau<<1)
         MaxTau; ///< Largest optical depth in splint (tau>>1)

  double 
    *prTau,   ///< array for log10(tau)
    *pirA,    ///< array for photoionisation rate in region A.
    *pirA_2,  ///< supplementary array for ionisation in region A.
    *pirB,    ///< array for photoionisation rate in region B.
    *pirB_2,  ///< supplementary array for ionisation in region B.
    *pirC,    ///< array for photoionisation rate in region C.
    *pirC_2,  ///< supplementary array for ionisation in region C.
    *phrA,    ///< array for photoheating rate in region A.
    *phrA_2,  ///< supplementary array for heating in region A.
    *phrB,    ///< array for photoheating rate in region B.
    *phrB_2,  ///< supplementary array for heating in region B.
    *phrC,    ///< array for photoheating rate in region C.
    *phrC_2;  ///< supplementary array for heating in region C.
  double
    *iltA,    ///< array for low-dtau PI rate in region A.
    *iltA_2,  ///< supplementary array for low-dtau PI in region A.
    *iltB,    ///< array for low-dtau PI rate in region B.
    *iltB_2,  ///< supplementary array for low-dtau PI in region B.
    *iltC,    ///< array for low-dtau PI rate in region C.
    *iltC_2,  ///< supplementary array for low-dtau PI in region C.
    *hltA,    ///< array for low-dtau PH rate in region A.
    *hltA_2,  ///< supplementary array for low-dtau PH in region A.
    *hltB,    ///< array for low-dtau PH rate in region B.
    *hltB_2,  ///< supplementary array for low-dtau PH in region B.
    *hltC,    ///< array for low-dtau PH rate in region C.
    *hltC_2;  ///< supplementary array for low-dtau PH in region C.


  ///
  /// For a specified region, return the values of the differences:
  /// I(tau)-I(tau+dtau)  and H(tau)-H(tau+dtau), or else use the
  /// low dtau approximations if dtau<TauMin.
  ///
  void get_region_HHe_integral_diff(
        const int,    ///< Which spectral region to calculate for.
        const double, ///< Tau to front edge of cell
        const double, ///< Delta-Tau through cell
        double *,     ///< photoionisation rate result
        double *     ///< photoheating rate result
        );

  ///
  /// Total optical depth as the sum of all species contributing,
  /// for a given frequency range.  The equations are from 
  /// Frank & Mellema (1994,A\&A,289,937).  This function returns
  /// the constant term normalising the power law, choosing the
  /// appropriate combination for each spectral region.
  ///
  double tau_total(
        const int, ///< which region we are in.
        const double, ///< optical depth to H0 at x=1.00
        const double, ///< optical depth to He0 at x=1.81
        const double  ///< optical depth to He+ at x=4.00
        );

  ///
  /// Integrand of the ionisation rate integral we pre-calculate.
  ///
  inline double Ix_integrand(
        const double x,   ///< x = E/13.6eV
        const double tau, ///< Tau_total at x.
        const double L    ///< L_x, source luminosity (erg/s/Hz) at x.
        ) { return L*exp(-tau)/x/h(); }

  ///
  /// Integrand of the photoheating rate integral we pre-calculate.
  ///
  inline double Hx_integrand(
        const double x,   ///< x = E/13.6eV
        const double tau, ///< Tau_total at x.
        const double L    ///< L_x, source luminosity (erg/s/Hz) at x.
        ) { return NuTh_H()*L*exp(-tau)*(x-1.0)/x; }

  ///
  /// Integral of I(x) from x=1 to x=xmax for a given value of Tau in
  /// the integral region A/B/C, assuming a blackbody spectrum.  Tau
  /// is the constant term multiplying the power-law scaling in the
  /// approximation of Frank & Mellema (1994).
  ///
  double Ix_BB_integral_ABC(
        const int,    ///< Spectral region we are in.
        const double, ///< Tau: constant prefactor in Tau(x).
        const double, ///< T, the source temperature (K)
        const double  ///< R, the source radius (cm)
        );

  ///
  /// Integral of I(x) from x=1 to x=xmax for a given value of Tau in
  /// the integral region A/B/C, assuming a blackbody spectrum.  Tau
  /// is the constant term multiplying the power-law scaling in the
  /// approximation of Frank & Mellema (1994).
  ///
  double Hx_BB_integral_ABC(
        const int,    ///< Spectral region we are in.
        const double, ///< Tau: constant prefactor in Tau(x).
        const double, ///< T, the source temperature (K)
        const double  ///< R, the source radius (cm)
        );

  ///
  /// Low dtau version of I(x) from x=1 to x=xmax for a given value
  /// of Tau in region A/B/C, assuming a blackbody spectrum. 
  ///
  double ILTx_BB_integral_ABC(
        const int,    ///< Spectral region we are in.
        const double, ///< Tau: constant prefactor in Tau(x).
        const double, ///< T, the source temperature (K)
        const double  ///< R, the source radius (cm)
        );

  ///
  /// Low dtau version of H(x) from x=1 to x=xmax for a given value
  /// of Tau in region A/B/C, assuming a blackbody spectrum. 
  ///
  double HLTx_BB_integral_ABC(
        const int,    ///< Spectral region we are in.
        const double, ///< Tau: constant prefactor in Tau(x).
        const double, ///< T, the source temperature (K)
        const double  ///< R, the source radius (cm)
        );

  ///
  /// Luminosity of a blackbody with input Temperature and radius.
  /// L(x,T) = 4*pi^2*R^2*B(x,T)
  /// units of erg/s/Hz
  ///
  inline double Lum_BB(
        const double x, ///< x, the photon frequency (nu_th)
        const double T, ///< T, the source temperature (K)
        const double R  ///< R, the source radius (cm)
        )
  {
    return 8.0*pi()*pi()*h()*pow(NuTh_H()*x,3.0)*R*R/c()/c()/
                      (exp(h()*NuTh_H()*x/kB()/T)-1.0);
  }

};

#endif //  if not EXCLUDE_MPV9
#endif //  HHE_PHOTOION_H

