/// \file NG_MPI_fine_to_coarse_boundaries.h
/// \brief Class declaration for NG_fine_to_coarse boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.24 JM: writing parallel algorithms

#ifndef NG_MPI_FINE_TO_COARSE_BOUNDARIES_H
#define NG_MPI_FINE_TO_COARSE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_fine_to_coarse_boundaries.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Implements NG_fine_to_coarse boundaries for a uniform grid,
/// parallelised with MPI.  The fine-level data are averaged and
/// sent to the coarse-level grid, where they replace whatever
/// was calculated on the coarse grid.
/// This class is set up on both fine and coarse grids, and one
/// calls the send function and the other the receive function.
///
class NG_MPI_fine_to_coarse_bc : virtual public NG_fine_to_coarse_bc {
protected:
  /// Assigns cells to a nested grid boundary whose data are to
  /// be averaged and sent to a coarser grid.
  int BC_assign_FINE_TO_COARSE_SEND(
      class SimParams &,  ///< pointer to simulation parameters
      const int,          ///< level of this grid.
      boundary_data *     ///< boundary data
  );

  /// Assigns cells to a nested grid boundary whose data are to
  /// be overwritten by data from a finer grid.
  int BC_assign_FINE_TO_COARSE_RECV(
      class SimParams &,  ///< pointer to simulation parameters
      const int,          ///< level of this grid.
      boundary_data *     ///< boundary data
  );

  /// List of IDs for all sends, should be cleared at the beginning
  /// of each timestep.
  std::vector<string> NG_F2C_send_list;

public:
  /// Send data from this grid to a coarser grid (maybe on another
  /// MPI process)
  int BC_update_FINE_TO_COARSE_SEND(
      class SimParams &,       ///< pointer to simulation parameters
      class FV_solver_base *,  ///< pointer to equations
      const int,               ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int);

  /// Receive data from a finer grid (maybe on another MPI
  /// process) and overwrite data on this grid.
  int BC_update_FINE_TO_COARSE_RECV(
      class SimParams &,       ///< pointer to simulation parameters
      class FV_solver_base *,  ///< pointer to equations
      const int,               ///< level of grid in NG grid struct.
      struct boundary_data *,
      const int,
      const int);

  ///
  /// Delete the temporary arrays used to send data to another
  /// MPI process
  void BC_FINE_TO_COARSE_SEND_clear_sends(class Sub_domain &);
};

#endif  // NG_MPI_FINE_TO_COARSE_BOUNDARIES_H
