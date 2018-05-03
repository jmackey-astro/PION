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
/// The simplest finite volume grid - a uniform grid with cells that
/// are cubes in the chosen coordinates.  This class sets up the grid and
/// other things to get a simulation ready to run, so it is useful
/// for simulation analysis.  PION itself uses a derived class to
/// setup and run simulations.
///
class setup_nested_grid
{
  public:
  setup_nested_grid();  ///< Simple constructor, initialises value.
  virtual ~setup_nested_grid(); ///< Deletes any dynamic memory, if not already done.


  ///
  /// Populate the array SimPM.nest_levels with Xmin,Xmax,Range,dx,etc.
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
  /// Determines what kind of boundary conditions are needed.
  /// Sets gp.Nbc to the appropriate value for the order of accuracy used.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int boundary_conditions(
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &,  ///< pointer to simulation parameters
      );   


  //---------------------------------------
  protected:
  //---------------------------------------

}; // setup_nested_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_NESTED_GRID_H
