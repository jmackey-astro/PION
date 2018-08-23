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

#ifndef SETUP_NESTED_GRID_H
#define SETUP_NESTED_GRID_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "decomposition/MCMD_control.h"
#include "boundaries/assign_update_bcs_NG.h"

///
/// Set up a static NG grid structure.  Serial code, so each
/// level of the NG has a single grid.
///
class setup_NG_grid :
  virtual public setup_fixed_grid,
  virtual public assign_update_bcs_NG
{
  public:
  setup_NG_grid();  ///< Simple constructor, initialises value.
  virtual ~setup_NG_grid(); ///< Deletes any dynamic memory, if not already done.


  ///
  /// Populate the array SimPM.levels with Xmin,Xmax,Range,dx,etc.
  ///
  void setup_NG_grid_levels(
      class SimParams &  ///< pointer to simulation parameters
      );

  /// 
  /// Sets up a NG grid.
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
      class MCMDcontrol &,  ///< unused for serial code
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );   



  //---------------------------------------
  protected:
  //---------------------------------------

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      const int          ///< level of grid in NG
      );


}; // setup_NG_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_NESTED_GRID_H
