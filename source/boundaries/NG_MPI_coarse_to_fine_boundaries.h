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

#include "boundaries/NG_coarse_to_fine_boundaries.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Implements NG_MPI_coarse_to_fine boundaries for a uniform grid.
///
class NG_MPI_coarse_to_fine_bc : virtual public NG_coarse_to_fine_bc {
public:
  /// Sends data from a coarser grid to set the external boundaries
  /// of any/all child grids.  If a child is on the same MPI process
  /// then do nothing.
  virtual int BC_update_COARSE_TO_FINE_SEND(
      class SimParams &,       ///< pointer to simulation parameters
      class GridBaseClass *,   ///< pointer to coarse-level grid
      class FV_solver_base *,  ///< pointer to equations
      const int,               ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int);

  /// Updates data to an external boundary from coarser grid.
  /// Receives data if parent is on a different MPI process or, if
  /// on the same process, just grab the data from the parent grid.
  virtual int BC_update_COARSE_TO_FINE_RECV(
      class SimParams &,       ///< pointer to simulation parameters
      class FV_solver_base *,  ///< pointer to equations
      const int,               ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int  ///< timestep on this (fine) grid
  );

  ///
  /// Delete the temporary arrays used to send data to another
  /// MPI process
  void BC_COARSE_TO_FINE_SEND_clear_sends(class Sub_domain &);

protected:
  /// Assigns cells to lists for sending to external boundaries of
  /// child cells.
  virtual int BC_assign_COARSE_TO_FINE_SEND(
      class SimParams &,  ///< pointer to simulation parameters
      const int,          ///< level of this grid.
      boundary_data *     ///< boundary data
  );

  /// Assigns data to an external boundary from coarser grid.  Sets
  /// up the book-keeping to allow the update to happen every timestep.
  virtual int BC_assign_COARSE_TO_FINE_RECV(
      class SimParams &,  ///< pointer to simulation parameters
      const int,          ///< level of this grid.
      boundary_data *     ///< boundary data
  );

  /// do what it says
  void add_cells_to_C2F_send_list(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to coarse-level grid
      struct c2f *bdata,      ///< pointer to list of cells
      int *ixmin,             ///< child grid xmin (integer)
      int *ixmax,             ///< child grid xmax (integer)
      const int,              ///< level of fine grid.
      const int *,            ///< level xmin of fine grid.
      const int *             ///< level xmax of fine grid.
  );

  /// do what it says for 1D grids
  void add_cells_to_C2F_send_list_1D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to coarse-level grid
      struct c2f *bdata,      ///< pointer to list of cells
      int *ixmin,             ///< child grid xmin (integer)
      int *ixmax,             ///< child grid xmax (integer)
      const int,              ///< level of fine grid.
      const int *,            ///< level xmin of fine grid.
      const int *             ///< level xmax of fine grid.
  );

  /// do what it says for 2D grids
  void add_cells_to_C2F_send_list_2D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to coarse-level grid
      struct c2f *bdata,      ///< pointer to list of cells
      int *ixmin,             ///< child grid xmin (integer)
      int *ixmax,             ///< child grid xmax (integer)
      const int,              ///< level of fine grid.
      const int *,            ///< level xmin of fine grid.
      const int *             ///< level xmax of fine grid.
  );

  /// do what it says for 3D grids
  void add_cells_to_C2F_send_list_3D(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to coarse-level grid
      struct c2f *bdata,      ///< pointer to list of cells
      int *ixmin,             ///< child grid xmin (integer)
      int *ixmax,             ///< child grid xmax (integer)
      const int,              ///< level of fine grid.
      const int *,            ///< level xmin of fine grid.
      const int *             ///< level xmax of fine grid.
  );

  /// List of IDs for all sends, should be cleared at the beginning
  /// of each timestep.
  std::vector<string> NG_C2F_send_list;
};

#endif  // NG_MPI_COARSE_TO_FINE_BOUNDARIES_H
