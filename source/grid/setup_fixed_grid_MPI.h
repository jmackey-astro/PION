/// \file setup_fixed_grid_MPI.h
/// 
/// \brief Declares setup_fixed_grid class for MPI-parallelised grids.
/// 
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.26 JM: split off from setup_fixed_grid.h

#ifndef SETUP_FIXED_GRID_MPI_H
#define SETUP_FIXED_GRID_MPI_H


#include "decomposition/MCMD_control.h"
#include "setup_fixed_grid.h"
#include "boundaries/assign_update_bcs_MPI.h"


#ifdef PARALLEL


/// 
/// This class reimplements some functions of setup_fixed_grid, so that they
/// work on multiple processors with the domain split between them.
/// 
class setup_fixed_grid_pllel :
  virtual public setup_fixed_grid,
  virtual public assign_update_bcs_MPI
{
  public:
   setup_fixed_grid_pllel();
   ~setup_fixed_grid_pllel();

  ///
  /// initialise the grid class with appropriate parameters.
  /// 
  /// This function sets up the appropriate grid; for parallel execution
  /// I need to define the domain of my grid, and then pass the appropriate
  /// parameters to the UniformGrid class.
  ///
  int setup_grid(
      vector<class GridBaseClass *> &,  ///< grid pointers.
      class SimParams &  ///< simulation parameters
      );

  ///
  /// Decide if I need to setup RT class, and do it if i need to.
  ///
  virtual int setup_raytracing(
      class SimParams &,  ///< pointer to simulation parameters
      class GridBaseClass * ///< pointer to computational grid 
      );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  virtual int boundary_conditions(
      class SimParams &,  ///< simulation parameters
      vector<class GridBaseClass *> &  ///< grid pointers.
      );

  protected:

  /// function to setup data-I/O class.
  virtual void setup_dataio_class(
      class SimParams &,  ///< pointer to simulation parameters
      const int  ///< type of I/O: 1=text,2=fits,5=silo
      );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *grid, ///< pointer to grid.
      const int          ///< unused
      );

}; // setup_fixed_grid_pllel

#endif // PARALLEL
#endif // if not SETUP_FIXED_GRID_MPI_H

