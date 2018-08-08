/// \file update_boundaries_nested.h
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for routines to update grid boundaries with
/// different boundary conditions.
///
/// Modifications:\n
/// - 2018.05.18 JM: started
///

#ifndef UPDATE_BOUNDARIES_NESTED_H
#define UPDATE_BOUNDARIES_NESTED_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/setup_nested_grid.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"
#include "sim_control/update_boundaries.h"



// ##################################################################
// ##################################################################


class update_boundaries_nested :
  virtual public setup_nested_grid,
  virtual public update_boundaries
{
  public:
  update_boundaries_nested();
  ~update_boundaries_nested();
   
 
  protected:
};

#endif // UPDATE_BOUNDARIES_NESTED_H

