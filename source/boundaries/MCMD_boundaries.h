/// \file periodic_boundaries.h
/// \brief Class declaration for periodic boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef MCMD_BOUNDARIES_H
#define MCMD_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "decomposition/MCMD_control.h"

///
/// Implements periodic boundaries for a uniform grid.
///
class MCMD_bc {
  protected:
  
  ///
  /// Get boundary data from other processor.
  /// Int tag is to distinguish between periodic and internal boundaries, 
  /// as the domain can be split between two processors, so that proc 0 is
  /// getting a periodic and an internal boundary from proc 1, and vice
  /// versa.
  ///
  virtual int BC_assign_BCMPI(
      class SimParams &,      ///< simulation parameters
      const int,              ///< level in grid hierarchy
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< pointer to boundary we are assigning.
      int              ///< tag, either BC_MPItag or BC_PERtag.
      );

  ///
  /// Updates data on an inter-process communicating boundary.
  ///
  virtual int BC_update_BCMPI(
      class SimParams &,      ///< simulation parameters
      const int,              ///< level in grid hierarchy
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int,  ///< final step.
      int         ///< tag, either BC_MPItag or BC_PERtag.
      );

  /// Makes a list of cells from which we grab data to send to
  /// a neighbouring MPI process during boundary updates.  This is
  /// set up before the simulation starts in order to save time
  /// during the update.
  int BC_select_data2send(
      class SimParams &, ///< simulation parameters
      class GridBaseClass *, ///< pointer to grid.
      list<cell *> *,   ///< list of cells (returned by this func.)
      int *,            ///< number of cells in list.
      boundary_data *,  ///< pointer to boundary data.
      int         ///< tag, either BC_MPItag or BC_PERtag.
      );

};

#endif // MCMD_BOUNDARIES_H

