/// \file flux_base.h
/// \author Jonathan Mackey
/// 
/// Contains class declarations for the basic flux-solver class.
/// This is an interface class, and also defines some basic
/// functionality, such as how to deal with tracer variables tacked
/// onto the end of the state vectors.
///
/// Modifications:\n
/// - 2010.12.23 JM: Moved from eqns_base.h into its own file.
///   Added tracer_flux() function.
///
/// - 2010.12.27 JM: Modified viscosity functions for new structure.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)


#ifndef FLUX_BASE_H
#define FLUX_BASE_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"



///
/// Abstract base class for a Flux Solver Class.  This class should be
/// grid-ignorant and geometry-ignorant, so just for calculating the 
/// flux across an interface given a left and right state.
///
class flux_solver_base : virtual public eqns_base 
{
 public:
  flux_solver_base(
      const int, ///< number of variables in state vector
      const double,  ///< coefficient of artificial viscosity.
      const int     ///< Number of tracer variables.
      );

  ~flux_solver_base();
  
protected:

};

#endif //FLUX_BASE_H

