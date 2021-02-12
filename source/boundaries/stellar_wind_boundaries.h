/// \file stellar_wind_boundaries.h
/// \brief Class declaration for stellar_wind boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef STELLAR_WIND_BOUNDARIES_H
#define STELLAR_WIND_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements stellar_wind boundaries for a uniform grid.  This is
/// a helper class that uses the stellar_wind class in
/// grid/stellar_wind_BC.h to set cell properties on the grid.
///
class stellar_wind_bc {
  protected:
    ///
    /// Add internal stellar wind boundaries -- these are (possibly
    /// time-varying) winds defined by a mass-loss-rate and a terminal velocity.
    /// A region within the domain is given fixed values corresponding to a
    /// freely expanding wind from a cell-vertex-located source. Optionally, the
    /// wind can be latitude-dependent.
    ///
    int BC_assign_STWIND(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*);

    ///
    /// Update internal stellar wind boundaries -- these are (possibly
    /// time-varying) winds defined by a mass-loss-rate and a terminal velocity.
    /// If fixed in time the wind is updated with b->refval, otherwise with a
    /// (slower) call to the stellar wind class SW
    ///
    int BC_update_STWIND(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        const double,          ///< current simulation time
        boundary_data*,        ///< Boundary to update.
        const int,             ///< current fractional step being taken.
        const int              ///< final step.
    );

    ///
    /// Add cells to both the Wind class, and to the boundary data list
    /// of cells.  This is re-defined for cylindrical and spherical
    /// coords below.
    ///
    virtual int BC_assign_STWIND_add_cells2src(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        // boundary_data *, ///< boundary ptr.
        const int  ///< source id
    );
};

#endif  // STELLAR_WIND_BOUNDARIES_H
