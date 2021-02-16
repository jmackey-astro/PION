///
/// \file MPv7.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running calculations with a simple
/// two-temperature isothermal equation of state, where T is T_low
/// when gas is neutral, and T_high when ionised, and linearly
/// interpolated for partial ionisation.
///
/// The file inherits from mpv3, and instead of the Wolfire et al.
/// (2003) and Henney et al. (2009) heating/cooling rates, it sets
/// the temperature based on Hydrogen ion fraction.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2018.03.20 JM: Renamed file.

#ifndef MPV7_H
#define MPV7_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "microphysics/MPv3.h"

class MPv7 : public MPv3 {
public:
  ///
  /// Constructor
  ///
  MPv7(
      const int,              ///< grid dimensions
      const int,              ///< Coordinate System flag
      const int,              ///< Total number of variables in state vector
      const int,              ///< Number of tracer variables in state vector.
      const std::string*,     ///< List of what the tracer variables mean.
      struct which_physics*,  ///< extra physics stuff.
      struct rad_sources*,    ///< radiation sources.
      const double            ///< EOS Gamma
  );

  ///
  /// Destructor
  ///
  ~MPv7();

  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW
  //-------------------
  //---------------------------------------------------------------------------
public:
  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  virtual int ydot(
      double,  ///< current time (probably not needed for rate equations)
      const N_Vector,  ///< current Y-value
      N_Vector,        ///< vector for Y-dot values
      const double*    ///< extra user-data vector, P, for evaluating
                       ///< ydot(y,t,p)
  );

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
      const int,      ///< ion index in tracer array (optional).
      const double*,  ///< input state vector (primitive).
      const double    ///< EOS gamma (optional)
  );

protected:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  virtual int convert_prim2local(
      const double*,  ///< primitive vector from grid cell (length nv_prim)
      double*         ///< local vector [x(H0),E](n+1).
  );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
      const double*,  ///< local (updated) vector [x(H0),E](n+1).
      const double*,  ///< input primitive vector from grid cell (length
                      ///< nv_prim)
      double*  ///< updated primitive vector for grid cell (length nv_prim)
  );

  ///
  /// returns gas temperature according to ion fraction
  ///
  virtual double get_temperature(
      const double,  ///< nH
      const double,  ///< E_int
      const double   ///< x(H+) N.B. This is ion fraction, not neutral
                     ///< fraction.
  );

  ///
  /// Returns total number of particles.
  ///
  double get_ntot(
      const double,  ///< nH
      const double   ///< x(H+) N.B. This is ion fraction, not neutral
                     ///< fraction.
  );

  /// Number of neutral particles per neutral H nucleon.
  double TTI_Nnt;
  double TTI_Mol;  ///< 0.5 if neutral H is molecular, otherwise 1.0.
  /// Low temperature for fully neutral medium.
  double TTI_Tlo;
  /// High temperature for fully ionised medium.
  double TTI_Thi;

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS
  //-------------------
  //---------------------------------------------------------------------------
};

#endif  // MPV7_H
