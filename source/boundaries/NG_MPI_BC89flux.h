/// \file NG_MPI_BC89flux.h
///
/// \brief Declares class that calculates fluxes leaving/entering
/// levels on a nested grid, to make sure that fluxes between
/// adjacent levels are consistent with each other.  This file
/// includes extra code needed for parallel execution.
///
/// \author Jonathan Mackey
///
/// Modifications:n
/// - 2019.12.04 JM: moved functions from uniform_grid and
///   sim_control_NG classes to gather everything in one file.

#ifndef NG_MPI_BC89FLUX_H
#define NG_MPI_BC89FLUX_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "NG_BC89flux.h"
#include "grid/uniform_grid.h"
#include "sim_constants.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Class to set up flux arrays at level boundaries in a nested grid
/// structure, to save the fluxes after a calculation of dU(), and to
/// make the fluxes between levels consistent after each timestep.
/// The implementation follows the paper of Berger & Colella
/// (1989,JCP,82,64).  This version inherits from NG_BC89flux and
/// adds code needed for parallel execution with domain decomposition
///
class NG_MPI_BC89flux : virtual public NG_BC89flux {
public:
  NG_MPI_BC89flux() {}
  ~NG_MPI_BC89flux() {}

  ///
  /// Setup the flux struct flux_update_recv with list of interfaces
  /// that need to be updated with fluxes from a finer level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels.
  /// This is an MPI-parallelised version which can deal with grids
  /// not on this process.
  ///
  int setup_flux_recv(
      class SimParams &,      ///< simulation params (including BCs)
      class GridBaseClass *,  ///< pointer to coarse grid.
      const int               ///< level to receive from
  );

  ///
  /// Setup the flux struct flux_update_send with list of interfaces
  /// that need to be sent to a coarser level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels.
  /// This is an MPI-parallelised version which can deal with grids
  /// not on this process.
  ///
  int setup_flux_send(
      class SimParams &,      ///< simulation params (including BCs)
      class GridBaseClass *,  ///< pointer to finer grid.
      const int               ///< level to send to
  );

  ///
  /// Send fine-level fluxes at level boundary to coarser parent
  /// grid(s) for static mesh refinement.
  ///
  int send_BC89_fluxes_F2C(
      class SimParams &,  ///< simulation params (including BCs)
      const int,          ///< My level in grid hierarchy.
      const int,          ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int           ///< Full order of accuracy of simulation
  );

  ///
  /// Receive fine-level fluxes at level boundary onto coarser parent
  /// grid(s) for static mesh refinement.
  ///
  int recv_BC89_fluxes_F2C(
      class FV_solver_base *,  ///< spatial solver, for gradients
      class SimParams &,       ///< simulation params (including BCs)
      const int,               ///< My level in grid hierarchy.
      const int,               ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int                ///< Full order of accuracy of simulation
  );

  /// clear the non-blocking MPI sends when they have been
  /// received.
  void clear_sends_BC89_fluxes(class Sub_domain &);

protected:
  /// List of IDs for MPI sends related to the BC89 Flux
  /// correction algorithms.  Should be cleared at the beginning
  /// of each timestep.
  std::vector<string> BC89_flux_send_list;
};

#endif  // NG_BC89FLUX_H
