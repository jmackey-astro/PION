/// \file reflecting_boundaries.h
/// \brief Class declaration for reflecting boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef REFLECTING_BOUNDARIES_H
#define REFLECTING_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"


///
/// Implements reflecting boundaries for a uniform grid.
///
class reflecting_bc {
  protected:
  
  /// Assigns data to a reflecting boundary.
  virtual int BC_assign_REFLECTING(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Updates data on a reflecting boundary.
  virtual int BC_update_REFLECTING(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );
};

#endif // REFLECTING_BOUNDARIES_H

