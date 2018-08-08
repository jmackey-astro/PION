/// \file double_Mach_ref_boundaries.h
/// \brief Class declaration for double_Mach_ref boundaries, to run
//         the Double Mach Reflection test problem.
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef DOUBLE_MACH_REF_BOUNDARIES_H
#define DOUBLE_MACH_REF_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"


///
/// Implements double_Mach_ref boundaries for a uniform grid.
///
class double_Mach_ref_bc {
  protected:
  
  /// Assigns data on a boundary for the Double Mach Reflection Problem.
  int BC_assign_DMACH(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on The other DMR test problem boundary
  virtual int BC_assign_DMACH2(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Updates data on the double mach reflection (DMR) boundary.
  int BC_update_DMACH(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const double,   ///< current simulation time (for DMACH)
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on the other DMR test problem boundary.
  virtual int BC_update_DMACH2(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );
};

#endif // DOUBLE_MACH_REF_BOUNDARIES_H

