/// \file setup_NG_grid.h
///
/// \brief Declares a class for setting up NG grids.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.02.09 JM: Split sim_control class into a setup class and
///   a derived class for running simulations.
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef SETUP_NG_GRID_H
#define SETUP_NG_GRID_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_BC89flux.h"
#include "boundaries/assign_update_bcs_NG.h"
#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Set up a static NG grid structure.  Serial code, so each
/// level of the NG has a single grid.
///
class setup_NG_grid :
    virtual public setup_fixed_grid,
    virtual public assign_update_bcs_NG,
    virtual public NG_BC89flux {
public:
  setup_NG_grid();
  virtual ~setup_NG_grid();

  ///
  /// Populate the array SimPM.levels with Xmin,Xmax,Range,dx,etc.
  ///
  virtual void setup_NG_grid_levels(
      class SimParams &  ///< pointer to simulation parameters
  );

  ///
  /// Sets up a NG grid.
  ///
  virtual int setup_grid(
      vector<class GridBaseClass *> &,  ///< grid pointers.
      class SimParams &                 ///< pointer to simulation parameters
  );

  ///
  /// Decide if I need to setup RT class and, if so, set up a
  /// raytracer associated with each grid.
  ///
  using setup_fixed_grid::setup_raytracing; /* TODO: is this intended to be
                                               brought here? */
  virtual int setup_raytracing(
      class SimParams &,               ///< simulation parameters
      vector<class GridBaseClass *> &  ///< grid pointers.
  );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  virtual int boundary_conditions(
      class SimParams &,               ///< pointer to simulation parameters
      vector<class GridBaseClass *> &  ///< grid pointers.
  );

  //---------------------------------------
protected:
  //---------------------------------------

  /// function to setup data-I/O class.
  virtual void setup_dataio_class(
      class SimParams &,  ///< pointer to simulation parameters
      const int           ///< type of I/O: 1=text,2=fits,5=silo
  );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,      ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      const int               ///< level of grid in NG
  );

  /// set flag for cells if they are not leaf cells (i.e. if there is
  /// a finer-level grid that covers the same volume).
  virtual void set_leaf_cells(
      vector<class GridBaseClass *> &,  ///< grid pointers.
      class SimParams &                 ///< pointer to simulation parameters
  );

};  // setup_NG_grid

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif  // if not SETUP_NESTED_GRID_H
