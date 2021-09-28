/// \file assign_update_bcs_NG_MPI.h
/// \brief Declares a class that inherits boundary types related to
///   static mesh-refinement (NG) and MCMD domain decomposition and
///   implements assignment and update functions.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.28 JM: adapted from assign_update_bcs_NG.cpp

#ifndef ASSIGN_UPDATE_BCS_NG_MPI_H
#define ASSIGN_UPDATE_BCS_NG_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_MPI_coarse_to_fine_boundaries.h"
#include "boundaries/NG_MPI_fine_to_coarse_boundaries.h"
#include "boundaries/assign_update_bcs_MPI.h"
#include "boundaries/assign_update_bcs_NG.h"
#include "boundaries/boundaries.h"
#include "decomposition/MCMD_control.h"
#include "spatial_solvers/solver_eqn_base.h"

class assign_update_bcs_NG_MPI :
    virtual public assign_update_bcs_MPI,
    virtual public assign_update_bcs_NG,
    virtual public NG_MPI_fine_to_coarse_bc,
    virtual public NG_MPI_coarse_to_fine_bc {
public:
  ///
  /// Assigns data to each boundary.
  ///
  int assign_boundary_data(
      class SimParams &,         ///< pointer to simulation parameters
      const int,                 ///< level in grid hierarchy
      class GridBaseClass *,     ///< pointer to grid.
      class microphysics_base *  ///< pointer to microphysics
  );

  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.
  ///
  int TimeUpdateExternalBCs(
      class SimParams &,       ///< pointer to simulation parameters
      const int,               ///< level in grid hierarchy
      class GridBaseClass *,   ///< pointer to grid.
      class FV_solver_base *,  ///< pointer to equations
      const double,            ///< current simulation time
      const int,               ///< Current step number in the timestep.
      const int                ///< Maximum step number in timestep.
  );
};

#endif  // ASSIGN_UPDATE_BCS_NG_MPI_H
