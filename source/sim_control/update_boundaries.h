/// \file update_boundaries.h
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for routines to update grid boundaries with
/// different boundary conditions.
///
/// Modifications:\n
/// - 2018.05.10 JM: moved code from uniform_grid.h
///

#ifndef UPDATE_BOUNDARIES_H
#define UPDATE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/setup_fixed_grid.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"



// ##################################################################
// ##################################################################


class update_boundaries : virtual public setup_fixed_grid
{
  public:
  update_boundaries();
  ~update_boundaries();
  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.
  ///
  virtual int TimeUpdateExternalBCs(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
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
      const double,   ///< current simulation time
      const int, ///< Current step number in the timestep.
      const int  ///< Maximum step number in timestep.
      );
   
 
  protected:
  ///
  /// Pointer to equations to solve, and routines for calculating
  /// fluxes on the grid.
  ///
  class FV_solver_base *spatial_solver;

  /// Updates data on a periodic boundary.
  virtual int BC_update_PERIODIC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on an outflow (zero gradient) boundary.
  int BC_update_OUTFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  ///
  /// Update the one-way outflow (zero gradient) boundary.
  /// If the flow is off-domain, then I use zero-gradient, but if flow
  /// is onto domain then I set the boundary cells to have zero normal 
  /// velocity.
  ///
  int BC_update_ONEWAY_OUT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on inflow boundary (data fixed to initial values).
  int BC_update_INFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on reflecting boundary.
  int BC_update_REFLECTING(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on fixed boundary (data doesn't change).
  int BC_update_FIXED(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on boundary where jet enters domain (keeps inflow fixed).
  virtual int BC_update_JETBC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on JetReflect boundary (see BC_assign_JETREFLECT description).
  int BC_update_JETREFLECT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on the double mach reflection (DMR) boundary.
  int BC_update_DMACH(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const double,   ///< current simulation time (for DMACH)
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  /// Updates data on the other DMR test problem boundary.
  virtual int BC_update_DMACH2(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );

  ///
  /// Update internal stellar wind boundaries -- these are (possibly time-varying)
  /// winds defined by a mass-loss-rate and a terminal velocity.  If fixed in time
  /// the wind is updated with b->refval, otherwise with a (slower) call to the 
  /// stellar wind class SW
  ///
  int BC_update_STWIND(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const double,   ///< current simulation time
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );
};

#endif // UPDATE_BOUNDARIES_H

