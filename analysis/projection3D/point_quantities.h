///
/// \file point_quantities.h
/// \author Jonathan Mackey
/// \date 2015.08.06
///
/// This file defines a number or functions to calculate physical
/// properties of a point in space, by bilinear interpolation from
/// four other coplanar points.
///
/// Modifications:
/// - 2015.08.06 JM: Trying to avoid code duplication in
///  point_velocity and image classes, so both can derive from this
///  class.
/// - 2015.08.19 JM: Added get_point_RotationMeasure()
/// - 2015.10.13 JM: added 6GHz Bremsstrahlung and Emission measure
/// - 2018.01.25 JM: added functions to request n(H+),n(H0),n(e-)

#ifndef POINT_QUANTITIES_H
#define POINT_QUANTITIES_H

#include "../xray/xray_emission.h"

///
/// integration point along line of sight of a pixel.  This point 
/// is to be placed along a plane of cells, so the state vector at the 
/// point can be obtained by linear interpolation between two cells to the 
/// left and right of the point (in that plane).
///
struct point_4cellavg {
  cell *ngb[4];  ///< pointers to four surrounding cells.
  pion_flt wt[4];  ///< weight of right cell value (wt of left=1-wt).
  //pion_flt pos[3]; ///< position of point, in image coords.
};


///
/// Class for evaluating quantities at integration points.
///
class point_quantities : public Xray_emission {
  public:
  ///
  /// constructor does nothing
  ///
  point_quantities() {}
  ///
  /// destructor does nothing
  ///
  virtual ~point_quantities() {}

  protected:

  ///
  /// Get the density at the point, based on 4 cell bilinear interpolation.
  ///
  double get_point_density(
      const struct point_4cellavg * ///< point
      );

  ///
  /// Get the e^- number density, using Microphysics function to get
  /// the electron number density at each simulation point and then
  /// interpolating.
  ///
  double get_point_electron_numberdensity(
      const struct point_4cellavg *
      );

  ///
  /// Get the H^+ number density, using Microphysics function to get
  /// the ionized H number density at each simulation point and then
  /// interpolating.
  ///
  double get_point_ionizedH_numberdensity(
      const struct point_4cellavg *
      );

  ///
  /// Get the H^0 number density, using Microphysics function to get
  /// the neutral H number density at each simulation point and then
  /// interpolating.
  ///
  double get_point_neutralH_numberdensity(
      const struct point_4cellavg *
      );

  ///
  /// Get the temperature at a point, based on 4-cell bilinear
  /// interpolation.  This uses a microphysics class to get T.
  ///
  double get_point_temperature(
      const struct point_4cellavg *, ///< point
      const double gamma   ///< EOS gamma
      );

  ///
  /// Get the Stokes Q component for the perpendicular magnetic field
  ///
  double get_point_StokesQ(
      struct point_4cellavg *, ///< pt
      const int, ///< ifrac
      const int, ///< bx index (image coords)
      const int, ///< by index (image coords)
      const int, ///< bz index (image coords)
      const int, ///< sign(xx)
      const int, ///< sign(yy)
      const int, ///< sign(zz)
      const double, ///< sin(theta)
      const double  ///< cos(theta)
      );

  ///
  /// Get the Stokes U component for the perpendicular magnetic field
  ///
  double get_point_StokesU(
      struct point_4cellavg *, ///< pt
      const int, ///< ifrac
      const int, ///< bx index (image coords)
      const int, ///< by index (image coords)
      const int, ///< bz index (image coords)
      const int, ///< sign(xx)
      const int, ///< sign(yy)
      const int, ///< sign(zz)
      const double, ///< sin(theta)
      const double  ///< cos(theta)
      );

  ///
  /// Get the |BX| component for the perpendicular magnetic field.
  /// This returns sqrt(n_H)*Bx^2/|B|, so it's a density weighted
  /// value and also diluted by the proportion of B on the los
  /// direction
  ///
  double get_point_BXabs(
      struct point_4cellavg *, ///< pt
      const int, ///< ifrac
      const int, ///< bx index (image coords)
      const int, ///< by index (image coords)
      const int, ///< bz index (image coords)
      const int, ///< sign(xx)
      const int, ///< sign(yy)
      const int, ///< sign(zz)
      const double, ///< sin(theta)
      const double  ///< cos(theta)
      );

  ///
  /// Get the |BY| component for the perpendicular magnetic field.
  /// This returns sqrt(n_H)*By^2/|B|, so it's a density weighted
  /// value and also diluted by the proportion of B on the los
  /// direction
  ///
  double get_point_BYabs(
      struct point_4cellavg *, ///< pt
      const int, ///< ifrac
      const int, ///< bx index (image coords)
      const int, ///< by index (image coords)
      const int, ///< bz index (image coords)
      const int, ///< sign(xx)
      const int, ///< sign(yy)
      const int, ///< sign(zz)
      const double, ///< sin(theta)
      const double  ///< cos(theta)
      );

  ///
  /// Get the Rotation measure along the line of sight.  When
  /// multiplied by the path length along each line segment (in pc)
  /// and by sqrt(4pi) this gives the RM in units of rad/m^2.
  /// It assumes n_e = n_H*y(H+), i.e. all electrons are from H-->H^+
  ///
  double get_point_RotationMeasure(
      struct point_4cellavg *, ///< pt
      const int, ///< bx index (image coords)
      const int, ///< bz index (image coords)
      const int, ///< sign(xx)
      const int, ///< sign(zz)
      const double, ///< sin(theta)
      const double  ///< cos(theta)
      );

  ///
  /// Get the absorption and emission coefficients for H-alpha
  /// radiation, according to a fit to Osterbrock's data table for
  /// photoionised nebulae.
  ///
  void get_point_Halpha_params(
      const struct point_4cellavg *, ///< point in question.
      const int, ///< ifrac index in prim.vec.
      const double,   ///< EOS gamma
      double *,  ///< absorption coefficient (/cm)
      double *   ///< emission coeff (erg/cm^3/s/sq.arcsec)
      );

  ///
  /// Get the absorption and emission coefficients for [N II] 6584AA
  /// radiation, according to a fit from Dopita (1973,A&A,29,387).
  ///
  void get_point_NII6584_params(
      const struct point_4cellavg *, ///< point in question.
      const int, ///< ifrac index in prim.vec.
      const double,   ///< EOS gamma
      double *,  ///< absorption coefficient (/cm)
      double *   ///< emission coeff (erg/cm^3/s/sq.arcsec)
      );

  ///
  /// Get the Emission Measure at a point, based on 4-cell bilinear
  /// interpolation.  Returns value in cm^{-6}
  ///
  double get_point_EmissionMeasure(
      const struct point_4cellavg *
      );

  ///
  /// Get the 6GHz Bremsstrahlung at a point, based on 4-cell bilinear
  /// interpolation.  Returns value in MJy/sr/cm.
  ///
  double get_point_Bremsstrahlung6GHz(
      const struct point_4cellavg *, ///< point
      const double   ///< EOS gamma
      );

  ///
  /// function to get the X-ray emissivity at a point, by reading
  /// from an interpolation table.
  ///
  void get_point_Xray_params(
      const struct point_4cellavg *, ///< point in question.
      const int, ///< ifrac index in prim.vec.
      const int, ///< which X-ray emissivity to use (index in array).
      const double,   ///< EOS gamma
      double *,  ///< absorption coefficient (/cm)
      double *   ///< emission coeff (erg/cm^3/s/sq.arcsec)
      );

};


#endif // POINT_QUANTITIES_H

