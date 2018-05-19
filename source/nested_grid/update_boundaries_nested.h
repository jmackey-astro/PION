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

#ifndef UPDATE_BOUNDARIES_H
#define UPDATE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "nested_grid/setup_nested_grid.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"



// ##################################################################
// ##################################################################


class update_boundaries_nested : virtual public setup_nested_grid
{
  public:
  update_boundaries_nested();
  ~update_boundaries_nested();
  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.
  ///
  virtual int TimeUpdateExternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const double,   ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );

  ///
  /// Runs through boundary cells which are grid cells and does
  /// the appropriate time update on them.
  ///
  virtual int TimeUpdateInternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const double,   ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );
   
 
  protected:


  /// Updates data to a nested grid from finer grid.
  virtual int BC_update_NEST_FINE( boundary_data *);

  /// Updates data to an external boundary from coarser grid.
  virtual int BC_update_NEST_COARSE(
      boundary_data *, ///< pointer to boudary struct
      const int,  ///< current fractional step being taken.
      const int   ///< full step.
      );

};

#endif // UPDATE_BOUNDARIES_H

