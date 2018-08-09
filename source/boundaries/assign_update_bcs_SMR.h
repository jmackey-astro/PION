/// \file assign_update_bcs_SMR.h
/// \brief Declares a class that inherits all boundary types
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef ASSIGN_UPDATE_BCS_SMR_H
#define ASSIGN_UPDATE_BCS_SMR_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "boundaries/assign_update_bcs.h"
#include "boundaries/SMR_coarse_to_fine_boundaries.h"
#include "boundaries/SMR_fine_to_coarse_boundaries.h"
#include "spatial_solvers/solver_eqn_base.h"


class assign_update_bcs_SMR : 
  virtual public assign_update_bcs,
  virtual public SMR_fine_to_coarse_bc,
  virtual public SMR_coarse_to_fine_bc
{
  public:  
  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      class GridBaseClass *,  ///< pointer to parent.
      class GridBaseClass *  ///< pointer to child.
      );


  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.
  ///
  virtual int TimeUpdateExternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const int, ///< level of grid in nested grid struct.
      class FV_solver_base *, ///< pointer to equations
      const double,   ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );

  ///
  /// Runs through boundary cells which are grid cells and does
  /// the appropriate time update on them.
  ///
  virtual int TimeUpdateInternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const int, ///< level of grid in nested grid struct.
      class FV_solver_base *, ///< pointer to equations
      const double,   ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );


};

#endif // ASSIGN_UPDATE_BCS_SMR_H

