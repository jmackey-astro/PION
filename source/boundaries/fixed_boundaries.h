/// \file fixed_boundaries.h
/// \brief Class declaration for fixed boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef FIXED_BOUNDARIES_H
#define FIXED_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements fixed boundaries for a uniform grid.
///
class fixed_bc {
protected:
  /// Assigns data to a fixed boundary.
  virtual int BC_assign_FIXED(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *);

  /// Updates data on a fixed boundary (data don't change).
  virtual int BC_update_FIXED(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,        ///< Boundary to update.
      const int,              ///< current fractional step being taken.
      const int               ///< final step.
  );
};

#endif  // FIXED_BOUNDARIES_H
