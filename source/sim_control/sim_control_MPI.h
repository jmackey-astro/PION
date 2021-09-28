/// \file sim_control_MPI.h
///
/// \brief Declares sim_control class for MPI-parallelised grids.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.26 JM: split off from sim_control.h
/// - 2015.02.18 JM: moved setup functions to setup_fixed_grid_MPI.h.

#ifndef SIM_CONTROL_MPI_H
#define SIM_CONTROL_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_fixed_grid_MPI.h"
#include "sim_control/sim_control.h"

#ifdef PARALLEL

/// The Parallel implementation of the Uniform FV Integrator.
///
/// This class reimplements some functions of sim_control, so that they
/// work on multiple processors with the domain split between them.
///
class sim_control_pllel :
    virtual public sim_control,
    virtual public setup_fixed_grid_pllel {
public:
  sim_control_pllel();
  ~sim_control_pllel();

  ///
  /// initialisation.
  ///
  /// This function checks if we are reading from single or multiple files,
  /// modifies the input file string accordingly, checks the file exists,
  /// and then calls the sim_control::Init() function.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  virtual int Init(
      string,    ///< Name of input file.
      int,       ///< Type of File (1=ASCII, 2=FITS, 5=Silo, ...).
      int,       ///< Number of command-line arguments.
      string *,  ///< Pointer to array of command-line arguments.
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
  );

  ///
  /// Time integration
  ///
  /// This is the main driver of the code -- It does the time integration
  /// until the stopping condition is reached and then returns.
  /// It calls a sequence of functions to advance the time by one timestep,
  /// all in a loop which runs until end-of-sim is reached.
  ///
  /// Parallel version has an AllReduce operation, where I check if the
  /// runtime of any processor is more than a fixed walltime, and if so set
  /// eosim to true and finish.  Some supercomputers have a maximum runtime
  /// limit for their simulations, and this allows to have a snapshot near the
  /// end of the allowed runtime.
  ///
  virtual int Time_Int(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
  );

protected:
  ///
  /// Calculate the appropriate timestep for all processes
  ///
  /// For a uniform grid, all cells have the same timestep equal to the
  /// minimum of all calculated steps.  This function calls the
  /// calc_timestep() function for the local grid, and then gets the min of
  /// all processes' local timesteps, and uses that as the timestep.
  ///
  int calculate_timestep(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      const int               ///< level in NG grid (if applicable)
  );

};  // sim_control_pllel

#endif  // PARALLEL
#endif  // if not SIM_CONTROL_MPI_H
