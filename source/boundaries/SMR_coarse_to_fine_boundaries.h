/// \file SMR_coarse_to_fine_boundaries.h
/// \brief Class declaration for SMR_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef SMR_COARSE_TO_FINE_BOUNDARIES_H
#define SMR_COARSE_TO_FINE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"


///
/// Implements SMR_coarse_to_fine boundaries for a uniform grid.
///
class SMR_coarse_to_fine_bc {
  protected:
  
  /// Assigns data to an external boundary from coarser grid.
  virtual int BC_assign_COARSE_TO_FINE(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to parent grid.
      );

  /// Updates data to an external boundary from coarser grid.
  virtual int BC_update_COARSE_TO_FINE(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in nested grid struct.
      boundary_data *, ///< pointer to boudary struct
      const int  ///< current step for fine grid (odd or even).
      );

  /// bilinear interpolation from a cell on the coarse grid to four
  /// cells on the fine grid (to the cell centres).
  virtual void bilinear_interp(
      class SimParams &,      ///< pointer to simulation parameters
      cell *,  ///< coarse level cell
      cell *,  ///< fine level cell
      const double *,  ///< prim. vec. at corner of coarse cell
      const double *,  ///< prim. vec. at corner of coarse cell
      const double *,  ///< prim. vec. at corner of coarse cell
      const double *   ///< prim. vec. at corner of coarse cell
      );
};


#endif // SMR_COARSE_TO_FINE_BOUNDARIES_H

