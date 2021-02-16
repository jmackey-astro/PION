/// \file NG_coarse_to_fine_boundaries.h
/// \brief Class declaration for NG_coarse_to_fine boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef NG_COARSE_TO_FINE_BOUNDARIES_H
#define NG_COARSE_TO_FINE_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Implements NG_coarse_to_fine boundaries for a uniform grid.
///
class NG_coarse_to_fine_bc {
protected:
  /// Assigns data to an external boundary from coarser grid.
  virtual int BC_assign_COARSE_TO_FINE(
      class SimParams&,      ///< pointer to simulation parameters
      class GridBaseClass*,  ///< pointer to grid.
      boundary_data*,        ///< boundary data
      class GridBaseClass*   ///< pointer to parent grid.
  );

  /// Updates data to an external boundary from coarser grid.
  virtual int BC_update_COARSE_TO_FINE(
      class SimParams&,       ///< pointer to simulation parameters
      class FV_solver_base*,  ///< pointer to equations
      const int,              ///< level of grid in NG grid struct.
      boundary_data*,         ///< pointer to boudary struct
      const int               ///< current step for fine grid (odd or even).
  );

  /// interpolate data from one coarse cell onto 2 fine cells in 1D
  virtual void interpolate_coarse2fine1D(
      class SimParams&,       ///< simulation parameters
      class GridBaseClass*,   ///< pointer to fine grid
      class FV_solver_base*,  ///< pointer to equations
      const pion_flt*,        ///< state vector of coarse cell.
      const pion_flt,         ///< volume of coarse cell.
      pion_flt*,              ///< dP/dx in coarse cell.
      cell*,                  ///< pointer to first fine cell  (XN)
      cell*                   ///< pointer to second fine cell (XP)
  );

  /// interpolate data from one coarse cell onto 4 fine cells in 2D
  virtual void interpolate_coarse2fine2D(
      class SimParams&,       ///< simulation parameters
      class GridBaseClass*,   ///< pointer to fine grid
      class FV_solver_base*,  ///< pointer to equations
      const pion_flt*,        ///< state vector of coarse cell.
      const int*,             ///< position of coarse cell.
      const pion_flt,         ///< volume of coarse cell.
      pion_flt*,              ///< dP/dx in coarse cell.
      pion_flt*,              ///< dP/dy in coarse cell.
      cell*,                  ///< pointer to first fine cell  (XN,YN)
      cell*,                  ///< pointer to second fine cell (XP,YN)
      cell*,                  ///< pointer to third fine cell  (XN,YP)
      cell*                   ///< pointer to fourth fine cell (XP,YP)
  );

  /// interpolate data from one coarse cell onto 8 fine cells in 3D
  virtual void interpolate_coarse2fine3D(
      class SimParams&,       ///< simulation parameters
      class GridBaseClass*,   ///< pointer to fine grid
      class FV_solver_base*,  ///< pointer to equations
      const pion_flt*,        ///< state vector of coarse cell.
      const int*,             ///< position of coarse cell.
      const pion_flt,         ///< volume of coarse cell.
      pion_flt*,              ///< dP/dx in coarse cell.
      pion_flt*,              ///< dP/dy in coarse cell.
      pion_flt*,              ///< dP/dz in coarse cell.
      cell**                  ///< pointers to 8 fine cells  (XN,YN)
  );

  /// For a coarse-grid cell with given position and optical depths,
  /// assign optical depths to a list of fine-grid child-cells, using
  /// the relative position of the source and coarse cell.
  void get_C2F_Tau(
      class SimParams&,     ///< simulation parameters
      std::vector<cell*>&,  ///< list of cells
      const pion_flt*,      ///< centre of coarse cell.
      pion_flt* T           ///< coarse-cell optical depths
  );

  int C2F_Nxd;             ///< number of extra data variables to send.
  vector<int> C2F_tauoff;  ///< offsets of optical depths from 0.
};

#endif  // NG_COARSE_TO_FINE_BOUNDARIES_H
