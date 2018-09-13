/// \file NG_MPI_coarse_to_fine_boundaries.h
/// \brief Class declaration for NG_MPI_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef NG_MPI_COARSE_TO_FINE_BOUNDARIES_H
#define NG_MPI_COARSE_TO_FINE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "boundaries/NG_coarse_to_fine_boundaries.h"

///
/// Implements NG_MPI_coarse_to_fine boundaries for a uniform grid.
///
class NG_MPI_coarse_to_fine_bc :
  virtual public NG_coarse_to_fine_bc
{
  protected:
  
  /// Assigns cells to lists for sending to external boundaries of
  /// child cells.
  virtual int BC_assign_COARSE_TO_FINE_SEND(
      class SimParams &,     ///< pointer to simulation parameters
      const int,  ///< level of this grid.
      boundary_data *  ///< boundary data
      );

  /// Sends data from a coarser grid to set the external boundaries
  /// of any/all child grids.  If a child is on the same MPI process
  /// then do nothing.
  virtual int BC_update_COARSE_TO_FINE_SEND(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int
      );

  /// Assigns data to an external boundary from coarser grid.  Sets
  /// up the book-keeping to allow the update to happen every timestep.
  virtual int BC_assign_COARSE_TO_FINE_RECV(
      class SimParams &,     ///< pointer to simulation parameters
      const int,  ///< level of this grid.
      boundary_data *  ///< boundary data
      );

  /// Updates data to an external boundary from coarser grid.
  /// Receives data if parent is on a different MPI process.  If parent
  /// is on the same process, just grab the data from the parent grid.
  virtual int BC_update_COARSE_TO_FINE_RECV(
      class SimParams &,      ///< pointer to simulation parameters
      class FV_solver_base *, ///< pointer to equations
      const int, ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int
      );

  /// do what it says
  void add_cells_to_C2F_send_list(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      );

  /// do what it says for 1D grids
  void add_cells_to_C2F_send_list_1D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      );

  /// do what it says for 2D grids
  void add_cells_to_C2F_send_list_2D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      );

  /// do what it says for 3D grids
  void add_cells_to_C2F_send_list_3D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      );
};


#endif // NG_MPI_COARSE_TO_FINE_BOUNDARIES_H
