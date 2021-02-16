/// \file hydrogen_recomb_Hummer94.h
/// \author Jonathan Mackey
///
/// Set up Spline interpolation table for Hummer (1994) recombination rate
/// and recombination cooling rate for Hydrogen.
///
/// Modifications:
/// - 2011.03.04 JM: Wrote class.
/// - 2011.03.14 JM: Renamed file.
/// - 2019.09.04 JM: changed spline interpolation class
///

#ifndef HUMMER94_HRECOMB_H
#define HUMMER94_HRECOMB_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

///
/// Set up Spline interpolation table for Hummer (1994) recombination rate
/// and recombination cooling rate for Hydrogen.
/// Functions return the rate between the limits of this table, or extrapolated
/// with the logarithmic slope at the end of the table if the temperature is too
/// large or too small.
///
/// Reference: Hummer, 1994, MNRAS, 268, 109.
///
class Hummer94_Hrecomb {
public:
  Hummer94_Hrecomb();
  ~Hummer94_Hrecomb();

  ///
  /// Recombination rate from H+ to H0, in cm^{3}/s, according to Hummer
  /// (1994)
  ///
  double Hii_rad_recomb_rate(double  ///< electron temperature.
  );

  ///
  /// Cooling rate due to recombination from H+ to H0, according to Hummer
  /// (1994), in units erg.cm^{3}/s
  ///
  double Hii_rad_recomb_cooling(double  ///< electron temperature.
  );

  ///
  /// Cooling rate due to ionised Hydrogen, being the sum of radiative
  /// recombination and free-free (Bremsstrahlung) cooling, beta_B^tot (last
  /// column of Table 1 in Hummer (1994).
  ///
  double Hii_total_cooling(double  ///< electron temperature.
  );

private:
  const double kB;  ///< Boltzmann constant.
  int hr_Nspl;      ///< number of interpolation points (31)
  double *hr_t,     ///< array for Temperature values.
      *hr_alpha,    ///< array for recomb rate.
                    //   *hr_alpha2, ///< array for recomb rate.
      *hr_beta,     ///< array for energy loss rate.
                    //   *hr_beta2,  ///< array for energy loss rate.
      *hr_btot;     ///< array for total energy loss rate.
                    //   *hr_btot2;  ///< array for total energy loss rate.
  double MinTemp;   ///< Minimum gas temperature which can be fit with spline.
  double MaxTemp;   ///< Maximum gas temperature which can be fit with spline.
  double MinSlope_alpha,  ///< Logarithmic slope for extrapolating spline fit.
      MaxSlope_alpha,     ///< Logarithmic slope for extrapolating spline fit.
      MinSlope_beta,      ///< Logarithmic slope for extrapolating spline fit.
      MaxSlope_beta,      ///< Logarithmic slope for extrapolating spline fit.
      MinSlope_btot,      ///< Logarithmic slope for extrapolating spline fit.
      MaxSlope_btot;      ///< Logarithmic slope for extrapolating spline fit.
  int hr_alpha_id, hr_beta_id, hr_btot_id;  ///< id for interpolation class.
};

#endif  // HUMMER94_HRECOMB_H
