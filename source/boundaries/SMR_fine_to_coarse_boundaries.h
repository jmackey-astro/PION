/// \file SMR_fine_to_coarse_boundaries.h
/// \brief Class declaration for SMR_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef SMR_FINE_TO_COARSE_BOUNDARIES_H
#define SMR_FINE_TO_COARSE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"


///
/// Implements SMR_fine_to_coarse boundaries for a uniform grid.
///
class SMR_fine_to_coarse_bc {
  protected:
  
  /// Assigns data to a nested grid from finer grid.
  virtual int BC_assign_FINE_TO_COARSE(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to child grid.
      );

  /// Updates data to a nested grid from finer grid.
  virtual int BC_update_FINE_TO_COARSE(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in nested grid struct.
      struct boundary_data *,
      const int,
      const int
      );
};


#endif // SMR_FINE_TO_COARSE_BOUNDARIES_H

