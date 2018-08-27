/// \file NG_fine_to_coarse_boundaries.h
/// \brief Class declaration for NG_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef NG_FINE_TO_COARSE_BOUNDARIES_H
#define NG_FINE_TO_COARSE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"


///
/// Implements NG_fine_to_coarse boundaries for a uniform grid.
///
class NG_fine_to_coarse_bc {
  protected:
  
  /// Assigns data to a NG grid from finer grid.
  virtual int BC_assign_FINE_TO_COARSE(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to child grid.
      );

  /// Updates data to a NG grid from finer grid.
  virtual int BC_update_FINE_TO_COARSE(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int
      );

  /// averages conserved vars of a list of cells, weighted by
  /// volume, using coarse-cell volume to normalise
  int average_cells(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      class GridBaseClass *, ///< fine-level grid
      const int,      ///< number of fine-level cells
      list<cell *> &, ///< list of cells
      pion_flt *      ///< [OUTPUT] averaged data (conserved var).
      );

};


#endif // NG_FINE_TO_COARSE_BOUNDARIES_H

