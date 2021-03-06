/// \file assign_update_bcs_MPI.h
/// \brief Declares a class that inherits all boundary types
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef ASSIGN_UPDATE_BCS_MPI_H
#define ASSIGN_UPDATE_BCS_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "boundaries/assign_update_bcs.h"
#include "decomposition/MCMD_control.h"
#include "periodic_boundaries_MPI.h"
#include "MCMD_boundaries.h"


class assign_update_bcs_MPI : 
  virtual public assign_update_bcs,
  virtual public periodic_pllel_bc,
  virtual public MCMD_bc
{
  public:  
  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      const int,              ///< level in grid hierarchy
      class GridBaseClass *   ///< pointer to grid.
      );

  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.
  ///
  virtual int TimeUpdateExternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      const int,              ///< level in grid hierarchy
      class GridBaseClass *,  ///< pointer to grid.
      class FV_solver_base *, ///< pointer to equations
      const double ,          ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );

};

#endif // ASSIGN_UPDATE_BCS_MPI_H

