/// \file reflecting_boundaries.h
/// \brief Class declaration for reflecting boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef RADIATIVE_SHOCK_TEST_H
#define RADIATIVE_SHOCK_TEST_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

///
/// Implements reflecting boundaries for a uniform grid.
///
class radiative_shock_test {
protected:
  /// Assigns data to a reflecting boundary of radiative shock test.
  virtual int BC_assign_RADSHOCK(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *);

  /// Updates data on a reflecting boundary of radiative shock test.
  /// It removes very dense material from near the wall to prevent
  /// pile-up of gas.
  virtual int BC_update_RADSHOCK(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,        ///< Boundary to update.
      const int,              ///< current fractional step being taken.
      const int               ///< final step.
  );
};

#endif  // RADIATIVE_SHOCK_TEST_H
