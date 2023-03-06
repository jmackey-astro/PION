///
/// \file cooling_SD93_cie.h
/// \author Jonathan Mackey
///
/// Modifications:
/// - 2011.03.04 JM: Added metal-only cooling for CIE.
/// - 2011.04.12 JM: Added Wiersma et al. (2009) newer CIE cooling, both
///    the full curve and the metals-only curve.
///
#ifndef COOLING_SD93_CIE_H
#define COOLING_SD93_CIE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"

///
/// Class which calculates the Sutherland & Dopita (1993) cooling
/// rates for collisional ionisation equilibrium.  This is a minimal
/// class which only sets up a spline interpolation and calculates
/// interpolated values.
///
/// UPDATE 2011.04.12 JM: Added the Wiersma
/// et al. (2009,MN,393,99) CIE cooling function (metals-only and total).
///
/// Update 2022.12.22 JM: Added two cooling rates for colliding wind systems
/// with the winds having different composition
///
/// Two curves can be set up, either the full solar metallicity CIE curve,
/// or the same curve but with the contribution of Hydrogen and Helium
/// subtracted off.  This is done by subtracting the metal-free CIE curve
/// from the solar metallicity curve.
///
/// N.B. for the metals-only curves.  Helium makes a significant contribution
/// to the Bremsstrahlung curve, being 0.4 of the aplitude of the H-curve.
/// Be sure to add back both H and He cooling when using metals-only.
///
/// References:
/// - Sutherland & Dopita (1993, ApJS, 88, 253)
/// - Wiersma, Schaye, Smith (2009,MN,393,99).
///
/// Data files:
/// - http://www.mso.anu.edu.au/~ralph/data/cool/m-00.cie
/// - http://www.mso.anu.edu.au/~ralph/data/cool/mzero.cie
/// - http://www.strw.leidenuniv.nl/WSS08/
/// - http://www.strw.leidenuniv.nl/WSS08/z_collis.txt
///
class cooling_function_SD93CIE : virtual public interpolate_arrays {
public:
  cooling_function_SD93CIE();
  ~cooling_function_SD93CIE();

  ///
  /// Set up the spline interpolation arrays.
  ///
  void setup_SD93_cie();

  ///
  /// Calculate rate for a given temperature.  Result is returned in
  /// units of erg.cm^{3}.s^{-1}
  ///
  double cooling_rate_SD93CIE(const double  ///< Input Temperature.
  );

  ///
  /// Set up the spline interpolation arrays for metals-only cooling.
  ///
  void setup_SD93_cie_OnlyMetals();

  ///
  /// Set up the spline interpolation arrays for metal-free cooling.
  /// This is probably never going to be useful...
  ///
  void setup_SD93_cie_MetalFree();

  ///
  /// Sets up spline interpolation for metals-only cooling from
  /// Wiersma, R.P.C., Schaye, J., and Smith, B.D. (2009,MNRAS,393,99)
  /// (arxiv:0807.3748).  This CIE function is from
  /// http://www.strw.leidenuniv.nl/WSS08/z_collis.txt
  /// and is a factor of a few lower than SD93-CIE.
  ///
  void setup_WSS09_CIE_OnlyMetals();

  ///
  /// Sets up spline interpolation for total cooling from
  /// Wiersma, R.P.C., Schaye, J., and Smith, B.D. (2009,MNRAS,393,99)
  /// (arxiv:0807.3748).  This CIE function is from
  /// http://www.strw.leidenuniv.nl/WSS08/z_collis.txt
  /// and is a factor of a few lower than SD93-CIE.
  ///
  void setup_WSS09_CIE();

  ///
  /// Set up spline interpolation for cooling based on Eatson+(2022,MNRAS,516,
  /// 6132) model for cooling from colliding wind binary system
  ///
  void setup_Eatson_CIE();

  /// Calculate rate for a given temperature for a two-gas mixture.
  /// Eatson+(2022,MNRAS,516,6132)
  /// Result is returned in units of erg.cm^{3}.s^{-1}
  double cooling_rate_Eatson(
      double T,       ///< Input Temperature.
      double wc_frac  ///< Fraction of gas (by mass) from WC star
  );

private:
  const int Nspl;         ///< number of tabulated values in spline table.
  double *Tarray;         ///< array for tabulated temperatures.
  double *Larray;         ///< Array for tabulated cooling rates.
  int spline_id;          ///< id of spline in interpolation class.
  vector<double> Tdata;   ///< vector for temperature data
  vector<double> L1data;  ///< vector for cooling rate of gas 1 (WC star)
  vector<double> L2data;  ///< vector for cooling rate of gas 2 (O star)
  int spline_id1;         ///< id of spline for gas 1.
  int spline_id2;         ///< id of spline for gas 2.

  double MaxTemp;         ///< Max. tabulated temperature.
  double MinTemp;         ///< Min. tabulated temperature.
  double MaxSlope;        ///< power law index for high-T extrapolation.
  double MinSlope;        ///< power law index for low-T  extrapolation.
  bool have_set_cooling;  ///< set to true if we have set up a cooling
                          ///< function.
};

#endif  // COOLING_SD93_CIE_H
