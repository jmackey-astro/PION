/// \file inflow_boundaries.h
/// \brief Class declaration for inflow boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef INFLOW_BOUNDARIES_H
#define INFLOW_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements inflow boundaries for a uniform grid.
///
class inflow_bc {
protected:
  /// Assigns data to a inflow (fixed) boundary.
  virtual int BC_assign_INFLOW(
      class SimParams&,      ///< pointer to simulation parameters
      class GridBaseClass*,  ///< pointer to grid.
      boundary_data*);

  /// Updates data on a inflow boundary, data fixed to initial values
  virtual int BC_update_INFLOW(
      class SimParams&,      ///< pointer to simulation parameters
      class GridBaseClass*,  ///< pointer to grid.
      boundary_data*,        ///< Boundary to update.
      const int,             ///< current fractional step being taken.
      const int              ///< final step.
  );
};

#endif  // INFLOW_BOUNDARIES_H
