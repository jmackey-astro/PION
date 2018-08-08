/// \file update_boundaries.h
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for routines to update grid boundaries with
/// different boundary conditions.
///
/// Modifications:\n
/// - 2018.05.10 JM: moved code from uniform_grid.h
///

#ifndef UPDATE_BOUNDARIES_H
#define UPDATE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/setup_fixed_grid.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"



// ##################################################################
// ##################################################################


class update_boundaries : virtual public setup_fixed_grid
{
  public:
  update_boundaries();
  ~update_boundaries();
 
  protected:

};

#endif // UPDATE_BOUNDARIES_H

