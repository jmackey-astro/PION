/// \file periodic_boundaries.h
/// \brief Class declaration for periodic boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef PERIODIC_BOUNDARIES_H
#define PERIODIC_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"


///
/// Implements periodic boundaries for a uniform grid.
///
class periodic_bc {
  protected:
  
  /// Assigns data to a periodic boundary.
  virtual int BC_assign_PERIODIC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Updates data on a periodic boundary.
  virtual int BC_update_PERIODIC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );
};

#endif // PERIODIC_BOUNDARIES_H

