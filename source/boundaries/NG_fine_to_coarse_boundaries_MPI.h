/// \file NG_fine_to_coarse_boundaries_MPI.h
/// \brief Class declaration for NG_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.24 JM: writing parallel algorithms

#ifndef NG_FINE_TO_COARSE_BOUNDARIES_MPI_H
#define NG_FINE_TO_COARSE_BOUNDARIES_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"


///
/// Implements NG_fine_to_coarse boundaries for a uniform grid,
/// parallelised with MPI.  The fine-level data are averaged and
/// sent to the coarse-level grid, where they replace whatever
/// was calculated on the coarse grid.
///
class NG_fine_to_coarse_MPI_bc :
  virtual public NG_fine_to_coarse_bc
{
  protected:
  
  /// Assigns data to a NG grid from finer grid.
  virtual int BC_assign_FINE_TO_COARSE_SEND(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to child grid.
      );

  /// Assigns data to a NG grid from finer grid.
  virtual int BC_assign_FINE_TO_COARSE_RECV(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to child grid.
      );

  /// Updates data to a NG grid from finer grid.
  virtual int BC_update_FINE_TO_COARSE_RECV(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int
      );

  /// Updates data to a NG grid from finer grid.
  virtual int BC_update_FINE_TO_COARSE_SEND(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int
      );
};


#endif // NG_FINE_TO_COARSE_BOUNDARIES_MPI_H

