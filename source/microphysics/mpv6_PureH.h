///
/// \file mpv6_PureH.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running the cosmological radiative transfer
/// test problems I (Iliev et al.,2006,MNRAS,371,1057) and II
/// (Iliev et al.,2009,MNRAS,400,1283).
///
/// The file inherits from mpv3, and instead of the Wolfire et al.
/// (2003) and Henney et al. (2009) heating/cooling rates, it uses
/// pure atomic hydrogen heating/cooling, and also sets the He and
/// metal abundances to zero (hardcoded!) regardless of the settings
/// in the parameterfile.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15


#ifndef MPV6_PUREH_H
#define MPV6_PUREH_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef  NEW_METALLICITY

#include "microphysics/mp_explicit_H.h"

class mpv6_PureH
  :
  public mp_explicit_H
{
  public:
  ///
  /// Constructor
  ///
  mpv6_PureH(
          const int,           ///< Total number of variables in state vector
	  const int,           ///< Number of tracer variables in state vector.
	  const std::string &, ///< List of what the tracer variables mean.
          struct which_physics * ///< extra physics stuff.
	  );

  ///
  /// Destructor
  ///
  ~mpv6_PureH();


  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW -------------------
  //---------------------------------------------------------------------------
  public:
  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  virtual int ydot(
          double,         ///< current time (probably not needed for rate equations)
          const N_Vector, ///< current Y-value
          N_Vector,       ///< vector for Y-dot values
          const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          );

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS -------------------
  //---------------------------------------------------------------------------
};





#endif // NEW_METALLICITY

#endif // MPV6_PUREH_H



