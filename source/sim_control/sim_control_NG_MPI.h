/// \file sim_control_MPI.h
/// 
/// \brief Declares sim_control class for MPI-parallelised grids.
/// 
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.26 JM: split off from sim_control.h
/// - 2015.02.18 JM: moved setup functions to setup_fixed_grid_MPI.h.

#ifndef SIM_CONTROL_NG_MPI_H
#define SIM_CONTROL_NG_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_grid_NG_MPI.h"
#include "sim_control/sim_control.h"
#include "sim_control/sim_control_NG.h"
#include "sim_control/sim_control_MPI.h"


#ifdef PARALLEL


/// The Parallel implementation of the Uniform FV Integrator.
/// 
/// This class reimplements some functions of sim_control, so that they
/// work on multiple processors with the domain split between them.
/// 
class sim_control_NG_MPI : 
  virtual public sim_control_pllel,
  virtual public sim_control_NG,
  virtual public setup_grid_NG_MPI
{
  public:
   sim_control_NG_MPI();
   ~sim_control_NG_MPI();

  ///
  /// initialisation of a PION simulation.
  ///
  /// This function checks if the input file exists, reads the header
  /// and sets up the grids, then reads the input file and puts the
  /// data on the grid.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int Init(
      string,   ///< Name of input file.
      int,      ///< Type of File (1=ASCII, 2=FITS, 5=Silo, ...).
      int,      ///< Number of command-line arguments.
      string *, ///< Pointer to array of command-line arguments.
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );

  ///
  /// Time integration
  ///
  /// Steps forward in time until the stopping condition is reached
  /// and then returns. It calls a sequence of functions to advance
  /// the time by one timestep, all in a loop which runs until 
  /// end-of-sim is reached.
  /// 
  int Time_Int(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );

  protected:

  ///
  /// Calculate the appropriate timestep for all processors
  /// 
  /// For a uniform grid, all cells have the same timestep equal to the minimum
  /// of all calculated steps.  This function calls the calc_timestep() function
  /// for the local grid, and then gets the min of all processor's local
  /// timesteps, and uses that as the timestep.
  ///
  //int calculate_timestep(
  //    class SimParams &,      ///< pointer to simulation parameters
  //    class GridBaseClass *, ///< pointer to grid.
  //    class FV_solver_base *, ///< solver/equations class
  //    const int ///< level in NG grid (if applicable)
  //    );

}; // sim_control_NG_MPI

#endif // PARALLEL
#endif // if not SIM_CONTROL_NG_MPI_H
