/// \file axisymmetric_boundaries.h
/// \brief Class declaration for axisymmetric boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef AXISYMMETRIC_BOUNDARIES_H
#define AXISYMMETRIC_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements axisymmetric boundaries for a uniform grid.
///
class axisymmetric_bc {
  protected:
    /// Assigns data to a axisymmetric boundary.
    virtual int BC_assign_AXISYMMETRIC(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*);

    /// Updates data on a axisymmetric boundary.
    virtual int BC_update_AXISYMMETRIC(
        class SimParams&,      ///< pointer to simulation parameters
        class GridBaseClass*,  ///< pointer to grid.
        boundary_data*,        ///< Boundary to update.
        const int,             ///< current fractional step being taken.
        const int              ///< final step.
    );
};

#endif  // AXISYMMETRIC_BOUNDARIES_H
