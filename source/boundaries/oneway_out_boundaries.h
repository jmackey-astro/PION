/// \file oneway_out_boundaries.h
/// \brief Class declaration for oneway_out boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef ONEWAY_OUT_BOUNDARIES_H
#define ONEWAY_OUT_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "boundaries/outflow_boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements oneway_out boundaries for a uniform grid.
/// This is like outflow/zero-gradient, except that inflow is not
/// permitted.
///
class oneway_out_bc : virtual public outflow_bc {
protected:
  ///
  /// Assigns data on a one-way outflow (zero gradient) boundary.
  /// If the flow is off-domain, then I use zero-gradient, but if flow
  /// is onto domain then I set the boundary cells to have zero normal
  /// velocity.
  ///
  int BC_assign_ONEWAY_OUT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *);

  ///
  /// Update the one-way outflow (zero gradient) boundary.
  /// If the flow is off-domain, then use zero-gradient, but if flow
  /// is onto domain then set the boundary cells to have zero normal
  /// velocity.
  ///
  int BC_update_ONEWAY_OUT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,        ///< Boundary to update.
      const int,              ///< current fractional step being taken.
      const int               ///< final step.
  );
};

#endif  // ONEWAY_OUT_BOUNDARIES_H
