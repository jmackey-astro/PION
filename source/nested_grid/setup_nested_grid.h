/// \file setup_nested_grid.h
/// 
/// \brief Declares a class for setting up nested grids.
/// 
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2015.02.09 JM: Split sim_control class into a setup class and
///   a derived class for running simulations.
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef SETUP_NESTED_GRID_H
#define SETUP_NESTED_GRID_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "decomposition/MCMD_control.h"

///
/// Set up a static nested grid structure.  Serial code, so each
/// level of the nest has a single grid.
///
class setup_nested_grid : virtual public setup_fixed_grid
{
  public:
  setup_nested_grid();  ///< Simple constructor, initialises value.
  virtual ~setup_nested_grid(); ///< Deletes any dynamic memory, if not already done.


  ///
  /// Populate the array SimPM.levels with Xmin,Xmax,Range,dx,etc.
  ///
  void setup_nested_grid_levels(
      class SimParams &  ///< pointer to simulation parameters
      );

  /// 
  /// Sets up a nested grid.
  ///
  int setup_grid(
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol *     ///< address of MCMD controller class.
      );

  ///
  /// Decide if I need to setup RT class and, if so, set up a
  /// raytracer associated with each grid.
  ///
  int setup_raytracing(
      class SimParams &,    ///< pointer to simulation parameters
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  int boundary_conditions(
      class SimParams &,  ///< pointer to simulation parameters
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );   

  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      class GridBaseClass *,  ///< pointer to parent.
      class GridBaseClass *  ///< pointer to child.
      );


  //---------------------------------------
  protected:
  //---------------------------------------

  /// equations class for assigning boundary data
  class eqns_base *eqn;

  /// Assigns data to a nested grid from finer grid.
  virtual int BC_assign_NEST_FINE(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to child grid.
      );

  /// Assigns data to an external boundary from coarser grid.
  virtual int BC_assign_NEST_COARSE(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *,  ///< boundary data
      class GridBaseClass *  ///< pointer to parent grid.
      );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      const int          ///< level of grid in nest
      );


}; // setup_nested_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_NESTED_GRID_H
