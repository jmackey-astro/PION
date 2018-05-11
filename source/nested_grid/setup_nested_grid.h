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
  /// Decide if I need to setup RT class and, if so, set up a
  /// raytracer associated with each grid.
  ///
  int setup_raytracing(
      class SimParams &,    ///< pointer to simulation parameters
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      vector<class RayTracingBase *> &  ///< address of vector of grid pointers.
      );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      const int,          ///< level of grid in nest
      std::vector<struct boundary_data *> &  ///< pointer to boundary data vector for this level
      );

  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      std::vector<struct boundary_data *> & ///< pointer to boundary structs
      );


  //---------------------------------------
  protected:
  //---------------------------------------
  /// vector of all boundaries, at each level.
  std::vector<std::vector<struct boundary_data *> > bdata_nest;  

  /// Assigns data to a nested grid from finer grid.
  virtual int BC_assign_NEST_FINE( boundary_data *);

  /// Assigns data to an external boundary from coarser grid.
  virtual int BC_assign_NEST_COARSE( boundary_data *);


}; // setup_nested_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_NESTED_GRID_H
