/// \file outflow_boundaries.h
/// \brief Class declaration for outflow boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef OUTFLOW_BOUNDARIES_H
#define OUTFLOW_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements outflow boundaries for a uniform grid.
///
class outflow_bc {
protected:
  /// Assigns data on an outflow (zero gradient) boundary.
  int BC_assign_OUTFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *);

  /// Updates data on an outflow (zero gradient) boundary.
  int BC_update_OUTFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,        ///< Boundary to update.
      const int,              ///< current fractional step being taken.
      const int               ///< final step.
  );
};

#endif  // OUTFLOW_BOUNDARIES_H
