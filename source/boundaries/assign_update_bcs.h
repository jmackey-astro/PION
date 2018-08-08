/// \file assign_update_bcs.h
/// \brief Declares a class that inherits all boundary types
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef ASSIGN_UPDATE_BCS_H
#define ASSIGN_UPDATE_BCS_H

#include "boundaries/boundaries.h"
#include "boundaries/periodic_boundaries.h"
#include "boundaries/outflow_boundaries.h"
#include "boundaries/oneway_out_boundaries.h"
#include "boundaries/inflow_boundaries.h"
#include "boundaries/reflecting_boundaries.h"
#include "boundaries/fixed_boundaries.h"
#include "boundaries/jet_boundaries.h"
#include "boundaries/jetreflect_boundaries.h"
#include "boundaries/double_Mach_ref_boundaries.h"
#include "boundaries/stellar_wind_boundaries.h"


class assign_update_bcs : 
  virtual public periodic_bc,
  virtual public outflow_bc,
  virtual public oneway_out_bc,
  virtual public inflow_bc,
  virtual public reflecting_bc,
  virtual public fixed_bc,
  virtual public jet_bc,
  virtual public jetreflect_bc,
  virtual public double_Mach_ref_bc,
  virtual public stellar_wind_bc
{
  public:  
  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *  ///< pointer to grid.
      );

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
};

#endif // ASSIGN_UPDATE_BCS_H

