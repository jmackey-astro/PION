/// \file periodic_boundaries_MPI.h
/// \brief Class declaration for periodic boundaries for parallel
///    code with domain decomposition.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.09 JM: moved code.

#ifndef PERIODIC_BOUNDARIES_MPI_H
#define PERIODIC_BOUNDARIES_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/MCMD_boundaries.h"
#include "boundaries/boundaries.h"
#include "boundaries/periodic_boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements periodic boundaries for a uniform grid.
///
class periodic_pllel_bc : virtual public periodic_bc, virtual public MCMD_bc {
  protected:
    ///
    /// Assigns data to a periodic boundary, getting data from another
    /// process if necessary.
    ///
    int BC_assign_PERIODIC(
        class SimParams&,      ///< pointer to simulation parameters
        const int,             ///< level in grid hierarchy
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*         ///< pointer to boundary data.
    );

    /// Updates data on a periodic boundary.
    int BC_update_PERIODIC(
        class SimParams&,      ///< pointer to simulation parameters
        const int,             ///< level in grid hierarchy
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*,        ///< Boundary to update.
        const int,             ///< current fractional step being taken.
        const int              ///< final step.
    );
};

#endif  // PERIODIC_BOUNDARIES_MPI_H
