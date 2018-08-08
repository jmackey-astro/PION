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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "decomposition/MCMD_control.h"
#include "setup_fixed_grid.h"


#ifdef PARALLEL

//
// integer flags for MPI communication labels.
//
#define BC_ANYtag 0 ///< works for either sort of communication.
#define BC_MPItag 1 ///< This is an integer tag on send/receive operations, to label that this communicates MPI boundary data.
#define BC_PERtag 2 ///< Integer tag to say it is for periodic BC.
#define BC_RTtag  3 ///< Integer tag to say we are transferring a radiative transfer column density tag.


/// 
/// This class reimplements some functions of setup_fixed_grid, so that they
/// work on multiple processors with the domain split between them.
/// 
class setup_fixed_grid_pllel : virtual public setup_fixed_grid
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
      class GridBaseClass **, ///< address of pointer to computational grid.
      class SimParams &,  ///< pointer to simulation parameters
      class MCMDcontrol *     ///< address of MCMD controller class.
      );

  ///
  /// Decide if I need to setup RT class, and do it if i need to.
  ///
  virtual int setup_raytracing(
      class SimParams &,  ///< pointer to simulation parameters
      class GridBaseClass *, ///< pointer to computational grid 
      class RayTracingBase * ///< pointer to raytracing class
      );

  protected:
  
  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *grid ///< pointer to grid.
      );

  ///
  /// Assigns data to each boundary.  Called by SetupBCs().
  ///
  virtual int assign_boundary_data(
      const double,   ///< current simulation time (for DMACH)
      const double, ///< Simulation start time.
      const double,  ///< Simulation finish time.
      const double ///< minimum temperature allowed
      );
  
  ///
  /// Assigns data to a periodic boundary, getting data from another
  /// process if necessary.
  ///
  virtual int BC_assign_PERIODIC(
      class SimParams &,     ///< pointer to simulation parameters
      class MCMDcontrol &,   ///< domain decomposition info
      class GridBaseClass *, ///< pointer to grid.
      boundary_data *        ///< pointer to boundary data.
      );
  

  /// Makes a list of cells from which we grab data to send to
  /// a neighbouring MPI process during boundary updates.  This is
  /// set up before the simulation starts in order to save time
  /// during the update.
  int BC_select_data2send(
      list<cell *> *,  ///< list of cells (returned by this func.)
      int *,         ///< number of cells in list.
      boundary_data *  ///< pointer to boundary data.
      );

  ///
  /// Get boundary data from other processor.
  /// Int tag is to distinguish between periodic and internal boundaries, 
  /// as the domain can be split between two processors, so that proc 0 is
  /// getting a periodic and an internal boundary from proc 1, and vice
  /// versa.
  ///
  virtual int BC_assign_BCMPI(
      class SimParams &,     ///< pointer to simulation parameters
      class MCMDcontrol &,   ///< domain decomposition info
      class GridBaseClass *, ///< pointer to grid.
      boundary_data *, ///< pointer to boundary we are assigning.
      int              ///< tag, either BC_MPItag or BC_PERtag.
      );

  ///
  /// Updates data on an inter-process communicating boundary.
  ///
  virtual int BC_update_BCMPI(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *, ///< Boundary to update.
      const int,  ///< current fractional step being taken.
      const int,  ///< final step.
      int         ///< tag, either BC_MPItag or BC_PERtag.
      );

}; // setup_fixed_grid_pllel

#endif // PARALLEL
#endif // if not SETUP_FIXED_GRID_MPI_H

