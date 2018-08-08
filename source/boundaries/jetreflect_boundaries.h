/// \file jetreflect_boundaries.h
/// \brief Class declaration for jetreflect boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef JETREFLECT_BOUNDARIES_H
#define JETREFLECT_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"


///
/// Implements jetreflect boundaries for a uniform grid.
/// This is the same as a reflecting boundary, except the B-field is
/// reflected differently. It is the base of a jet, so there is
/// assumed to be another jet in the opposite direction, so the
/// normal B-field is unchanged across the boundary and the
/// tangential field is reversed.
///
class jetreflect_bc {
  protected:
  
  ///
  /// Assigns data on a JetReflect boundary, which is the same
  /// as a reflecting boundary, except the B-field is reflected differently.
  /// It is the base of a jet, so there is assumed to be another jet
  /// in the opposite direction, so the normal B-field is unchanged across
  /// the boundary and the tangential field is reversed.
  ///
  virtual int BC_assign_JETREFLECT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Updates data on JetReflect boundary (see
  /// BC_assign_JETREFLECT description).
  virtual int BC_update_JETREFLECT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int   ///< final step.
      );
};

#endif // JETREFLECT_BOUNDARIES_H

