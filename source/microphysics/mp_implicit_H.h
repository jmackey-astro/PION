///
/// \file mp_implicit_H.h
/// \author Jonathan Mackey
/// \date 2011.10.12
///
/// This file updates the implicit (C2-ray type) integrator in
/// microphysics.cc with a new version based on the explicit
/// integrator and capable of treating multi-frequency
/// photoionisation.
///
/// modifications:
/// - getting it written: mods up until 2011.10.XX
/// - 2011.10.17 JM: Debugging.
/// - 2013.02.14 JM: Tidied up file.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).



#ifndef MP_IMPLICIT_H_H
#define MP_IMPLICIT_H_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifndef EXCLUDE_MPV4


#include <vector>
#include "microphysics/microphysics_base.h"
#include "microphysics/cooling_SD93_cie.h"
#include "microphysics/hydrogen_mp.h"
#include "microphysics/hydrogen_recomb_Hummer94.h"
#include "microphysics/hydrogen_photoion.h"
#include "microphysics/cvode_integrator.h"
#include "microphysics/mp_explicit_H.h"

class mp_implicit_H
  :
  public mp_explicit_H
{
  public:
  ///
  /// Constructor
  ///
  mp_implicit_H(
      const int,  ///< grid dimensions
      const int,  ///< Coordinate System flag
      const int,  ///< Total number of variables in state vector
      const int,  ///< Number of tracer variables in state vector.
      const std::string *, ///< List of what the tracer variables mean.
      struct which_physics *, ///< extra physics stuff.
      struct rad_sources *,    ///< radiation sources.
      const double  ///< EOS Gamma
      );

  ///
  /// Destructor
  ///
  ~mp_implicit_H();

  ///
  /// This takes a copy of the primitive vector and advances it in
  /// time over the step requested, and at the end copies the updated
  /// vector into the destination vector.  For fully local
  /// microphysics but WITH radiative transfer, where the column
  /// densities for diffuse and direct radiation are included as
  /// parameters.  The input list of column densities is ordered by
  /// the number of sources in each category in the vector of
  /// integers.
  ///
  /// Integers refer to:
  /// - Number of diffuse ionising sources (at infinity), 
  /// - Number of diffuse UV sources (at infinity),
  /// - Number of ionising point sources,
  /// - Number of UV point sources.
  ///
  /// Returned data is the time-averaged column density
  /// <<\delta\rho\times\delta x>>
  ///
  virtual int TimeUpdateMP_RTnew(
          const pion_flt *, ///< Primitive Vector to be updated.
          const int,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &,
          ///< list of UV-heating column densities and source properties.
          const int,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &,
          ///< list of ionising src column densities and source properties.
          pion_flt *,       ///< Destination Vector for updated values
                         ///< (can be same as first Vector.
          const double,   ///< Time Step to advance by.
          const double,   ///< EOS gamma.
          const int, ///< Switch for what type of integration to use.
                    ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
          double *    ///< any returned data (final temperature?).
          );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  (photoionisation time is not implemented).
  ///
  double timescales(
          const pion_flt *, ///< Current cell.
          const double,   ///< EOS gamma.
          const bool, ///< set to 'true' if including cooling time.
          const bool, ///< set to 'true' if including recombination time.
          const bool  ///< set to 'true' if including photo-ionsation time.
          );

  ///
  /// This returns the minimum timescale of all microphysical 
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.  It requires the
  /// radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.
  ///
  double timescales_RT(
          const pion_flt *, ///< Current cell.
          const int,      ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &,
          ///< list of UV-heating column densities and source properties.
          const int,      ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &,
          ///< list of ionising src column densities and source properties.
          const double   ///< EOS gamma.
          );

  protected:
  ///
  /// convert state vector from grid cell into local microphysics
  /// vector.
  ///
  virtual int convert_prim2local(
          const pion_flt *, ///< primitive vector from grid cell (length nv_prim)
          double *        ///< local vector [x,E](n+1).
          );

  ///
  /// Convert int(exp(-dtau)dt) into time-averaged value of
  /// rho*(1-x)*ds
  ///
  double get_timeaveraged_rhodx(
            const double, ///< this is int(exp(-dtau)dt)
            const double  ///< this is dt
            );
          
  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  virtual void setup_local_vectors();

  int lv_dtau; ///< index of dtau variable in local vector.

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

  ///
  /// Get the number of extra parameters and the number of equations.
  ///
  virtual void get_problem_size(int *, ///< number of equations
                        int *  ///< number of parameters in user_data vector.
                        );
  protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  virtual void get_error_tolerances(
          double *, ///< relative error tolerance (single value)
          double []  ///< absolute error tolerance (array)
          );

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS -------------------
  //---------------------------------------------------------------------------
};

#endif // if not excluding MPv4


#endif // MP_IMPLICIT_H_H


