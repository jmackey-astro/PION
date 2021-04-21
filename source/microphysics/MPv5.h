///
/// \file MPv5.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for modelling photoevaporation of dense molecular
/// clouds.
/// The file inherits from MPv3, and instead of the Wolfire et al.
/// (2003) neutral gas heating/cooling rates, it uses the Henney et
/// al. (2009) molecular cooling rate and heating rates.  So the only
/// significant difference is in the Ydot() function.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2018.03.20 JM: Renamed file.

#ifndef MPV5_H
#define MPV5_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "microphysics/MPv3.h"

class MPv5 : public MPv3 {
public:
  ///
  /// Constructor
  ///
  MPv5(
      const int,               ///< grid dimensions
      const int,               ///< Coordinate System flag
      const int,               ///< Total number of variables in state vector
      const int,               ///< Number of tracer variables in state vector.
      const std::string *,     ///< List of what the tracer variables mean.
      struct which_physics *,  ///< extra physics stuff.
      struct rad_sources *,    ///< radiation sources.
      const double             ///< EOS Gamma
  );

  ///
  /// Destructor
  ///
  ~MPv5();

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
      const double *   ///< extra user-data vector, P, for evaluating
                       ///< ydot(y,t,p)
  );

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS
  //-------------------
  //---------------------------------------------------------------------------
};

#endif  // MPV5_H
