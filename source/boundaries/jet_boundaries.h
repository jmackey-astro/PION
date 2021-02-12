/// \file jet_boundaries.h
/// \brief Class declaration for jet boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef JET_BOUNDARIES_H
#define JET_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements jet boundaries for a uniform grid.
///
class jet_bc {
  protected:
    /// Sets some boundary cells to be fixed to the Jet inflow
    /// condition.
    virtual int BC_assign_JETBC(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*);

    /// Updates data on boundary where jet enters domain (keeps inflow
    /// fixed).
    virtual int BC_update_JETBC(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*,        ///< Boundary to update.
        const int,             ///< current fractional step being taken.
        const int              ///< final step.
    );
};

#endif  // JET_BOUNDARIES_H
