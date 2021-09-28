///
/// \file hydrogen_photoion.h
/// \author Jonathan Mackey
/// \date 2011.03.14
///
/// \description
/// This contains the class declaration for hydrogen_photoion, which calculates
/// the multifrequency photoionisation rate as a function of optical depth,
/// including spectral hardening due to preferential absorption of lower energy
/// photons.  The algorithm is the discretised photon-conserving formulation of
/// Mellema et al (2006,NewA,11,374).
///
/// Modifications:
/// - 2011.03.29 JM: Added units to the function descriptions.
/// - 2011.04.17 JM: Added discretised photoionisation rate for monochromatic
///    radiation source. (fixed it on 18/4)
/// - 2011.05.02 JM: Added MinCol, MaxCol to check for spline overruns on the
///    multifreq photoionisation and photoheating rates.
/// - 2011.05.04 JM: Added a discretised multifreq photoionisation rate with
///    an approximation for dtau<<1.
/// - 2011.10.08 JM: Added switch to use interpolate.spline/splint instead of
/// the local
///    STL vector one, because the vector functions are slower by about 2.5%.
/// - 2014.03.27 JM: fixed bug in discrete monochromatic PI rate.
/// - 2019.09.04 JM: removed vector interpolation because not clear wehre
///    code comes from.
///
#ifndef HYDROGEN_PHOTOION
#define HYDROGEN_PHOTOION

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"

#define JUST_IONISED                                                           \
  2.178721e-11  ///< This is (1.0+5e-7) times the ionisation energy of H

///
/// This class calculates
/// the multifrequency photoionisation rate as a function of optical depth,
/// including spectral hardening due to preferential absorption of lower energy
/// photons.  The algorithm is the discretised photon-conserving formulation of
/// Mellema et al (2006,NewA,11,374).
///
/// Nspl should be set to 75, Nsub to 800, to give 1 per cent errors over most
/// of the relevant column density range for Hot stars.  Error is up to 2 per
/// cent for a star with Teff=10000K, assuming it is a blackbody and emits
/// ionising photons.  More realistic spectra could be added for such cool
/// stars, since a blackbody is not a good approximation for the UV emitted
/// flux.  For hot stars an integration up to 1000eV is recommended, although
/// the results are not very sensitive to this.  Cutting off at 500 eV gives
/// almost indistinguishable results, but at 250 eV the difference is
/// noticeable.  For cool stars the cutoff could probably be lower.
///
class hydrogen_photoion : virtual public interpolate_arrays {
public:
  hydrogen_photoion();
  ~hydrogen_photoion();

  ///
  /// For given stellar parameters, set up the photoionisation and
  /// photoheating spline interpolation tables for the rates as a function of
  /// column density. Uses discretised photon-conserving formulation of
  /// Mellema et al (2006,NewA,11,374).
  ///
  void Setup_photoionisation_rate_table(
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double,  ///< Stellar luminosity (ergs/s) (normalises BB curve).
      const double,  ///< Min optical depth at 13.6eV
      const double,  ///< Max optical depth at 13.6eV
      const double,  ///< Max energy to integrate to.
      const int,     ///< Number of sub-points in energy integration.
      const int      ///< Number of sub-points in spline interpolation.
  );

  ///
  /// Get DISCRETISED multifrequency photoionisation rate for a given column
  /// density of H0, local number density of H0, and Shell volume
  /// ~4.Pi.[(R+)^3-(R-)^3]/3 Uses discretised photon-conserving formulation
  /// of Mellema et al (2006,NewA,11,374).
  ///
  /// N.B. The result is returned as the energy added per second per H
  /// nucleon.
  ///
  double Hi_discrete_multifreq_photoion_rate(
      const double,  ///< Optical depth of H0 to front edge of cell
                     ///< (at 13.6eV).
      const double,  ///< Optical depth of H0 through cell (at 13.6eV)
      const double,  ///< Local number density of H (per cm3) n(H).
      const double,  ///< path length through cell, ds
      const double   ///< Shell volume (cm3).
  );
  ///
  /// Get DISCRETISED multifrequency photoheating rate for a given optical
  /// depth of H0, local number density of H0, and Shell volume
  /// ~4.Pi.[(R+)^3-(R-)^3]/3 Uses discretised photon-conserving formulation
  /// of Mellema et al (2006,NewA,11,374).
  ///
  /// N.B. The result is returned as the number of ionisations per second per
  /// H nucleon.
  ///
  double Hi_discrete_multifreq_photoheating_rate(
      const double,  ///< Optical depth of H0 to front edge of cell
                     ///< (at 13.6eV).
      const double,  ///< Optical depth of H0 through cell (at 13.6eV)
      const double,  ///< Local number density of H (per cm3) n(H).
      const double,  ///< path length through cell, ds
      const double   ///< Shell volume (cm3).
  );

  ///
  /// Get multifrequency photoionisation rate for a given optical depth of H0,
  /// local number density of H0, and Shell volume ~4.Pi.[(R+)^3-(R-)^3]/3
  /// Uses discretised photon-conserving formulation of Mellema et al
  /// (2006,NewA,11,374).
  ///
  /// N.B. The result is returned as the number of ionisations per second per
  /// neutral H atom.
  ///
  double Hi_multifreq_photoionisation_rate(
      const double,  ///< Optical depth of H0 (at 13.6eV).
      const double,  ///< Local number density of H0 (per cm3).
      const double   ///< Shell volume (cm3).
  );
  ///
  /// Get multifrequency photoionisation heating rate for a given optical
  /// depth of H0, local number density of H0, and Shell volume
  /// ~4.Pi.[(R+)^3-(R-)^3]/3 Uses discretised photon-conserving formulation
  /// of Mellema et al (2006,NewA,11,374).
  ///
  /// N.B. The result is returned as the energy added per second per neutral H
  /// atom.
  ///
  double Hi_multifreq_photoionisation_heating_rate(
      const double,  ///< Optical depth of H0 (at 13.6eV).
      const double,  ///< Local number density of H0 (per cm3).
      const double   ///< Shell volume (cm3).
  );

  ///
  /// Photoionisation cross section of atomic Hydrogen.  The formula is from
  /// Sutherland & Dopita (2003,Textbook,eq.5.32), which I think is the same
  /// as that given in Osterbrock (1989).
  ///
  double Hi_monochromatic_photo_ion_xsection(const double);

  ///
  /// Photoionisation cross section of atomic Hydrogen.  The formula is from
  /// Sutherland & Dopita (2003,Textbook,eq.5.32), which I think is the same
  /// as that given in Osterbrock (1989).  This function gives the ratio of
  /// the cross section at energy E to that at the ionisation edge.
  ///
  double Hi_monochromatic_photo_ion_xsection_fractional(const double);

  ///
  /// Discretised photoionisation rate according to Eq.6 of Mellema et al.
  /// (2006,NewA,11,374) in the limit of a monochromatic ionising source.
  ///
  /// N.B. The result is returned as the number of ionisations per
  /// second per H nucleon.
  ///
  double Hi_discrete_mono_photoion_rate(
      const double,  ///< Optical depth of H0 (at 13.6eV) to front edge of
                     ///< cell.
      const double,  ///< Optical depth of H0 (at 13.6eV) through cell
      const double,  ///< Local number density of H0 (per cm3) n(H0).
      const double,  ///< ionising photon luminosity of source
      const double,  ///< Energy of photons (erg)
      const double,  ///< path length through cell, ds
      const double   ///< Shell volume (cm3).
  );

private:
  ///
  /// This returns the value of the photoionisation rate integrand for the
  /// discretised rate, without the n(H0)*Vshell in the denominator.
  ///
  double photoion_rate_source_integral(
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double,  ///< Optical depth of H0 (at 13.6eV)
      const double,  ///< Max energy to integrate to.
      const int      ///< Number of sub-points in integration.
  );

  ///
  /// This returns the value of the photoheating rate integral for the
  /// discretised rate, without the n(H0)*Vshell in the denominator.
  ///
  double photoheating_rate_source_integral(
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double,  ///< Optical depth of H0 (at 13.6eV)
      const double,  ///< Max energy to integrate to.
      const int      ///< Number of sub-points in integration.
  );

  ///
  /// This returns the value of the photoionisation rate integrand for the
  /// discretised rate, without the n(H0)*Vshell in the denominator.
  ///
  double photoion_rate_source_integrand(
      const double,  ///< photon energy (ergs).
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double   ///< Optical depth of H0 (at 13.6eV)
  );

  ///
  /// This returns the value of the photoheating rate integrand for the
  /// discretised rate, without the n(H0)*Vshell in the denominator.
  ///
  double photoheating_rate_source_integrand(
      const double,  ///< photon energy (ergs).
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double   ///< Optical depth of H0 (at 13.6eV)
  );

  ///
  /// discretised version for dtau<<1, containing sigma(E) in the integrand.
  ///
  double PI_LowTau_rate_source_integral(
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double,  ///< Optical depth of H0 (at 13.6eV)
      const double,  ///< Max energy to integrate to.
      const int      ///< Number of sub-points in integration.
  );

  ///
  /// discretised version for dtau<<1, containing sigma(E) in the integrand.
  ///
  double PH_LowTau_rate_source_integral(
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double,  ///< Optical depth of H0 (at 13.6eV)
      const double,  ///< Max energy to integrate to.
      const int      ///< Number of sub-points in integration.
  );

  ///
  /// discretised version for dtau<<1, containing sigma(E) in the integrand.
  ///
  double PI_LowTau_rate_source_integrand(
      const double,  ///< photon energy (ergs).
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double   ///< Optical depth of H0 (at 13.6eV)
  );

  ///
  /// discretised version for dtau<<1, containing sigma(E) in the integrand.
  ///
  double PH_LowTau_rate_source_integrand(
      const double,  ///< photon energy (ergs).
      const double,  ///< BB temperature (K)
      const double,  ///< Radius of star (cm)
      const double   ///< Optical depth of H0 (at 13.6eV)
  );
  // private:
  int PI_Nspl;    ///< number of elements in spline integral.
  double MinTau,  ///< Smallest optical depth in spline int (must be tau<<1)
      MaxTau;     ///< Largest optical depth in spline int (must be tau>>1)
  double *PI_Tau, *PIrate, *PIheat;
  double *LTPIrate, *LTPIheat;
  int PIrt_id, PIht_id, LTPIrt_id, LTPIht_id;  ///< ID in interpolation class.
};

#endif  // HYDROGEN_PHOTOION
