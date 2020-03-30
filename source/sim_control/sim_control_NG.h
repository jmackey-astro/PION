/// \file sim_control_NG.h
/// 
/// \brief Declares grid parameter class, and grid methods classes.
/// 
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.05.04 JM: worked on class.

#ifndef SIM_CONTROL_NG_H
#define SIM_CONTROL_NG_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "grid/setup_NG_grid.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"
#include "sim_control/sim_control.h" 
#include "sim_control/calc_timestep.h"
#include "sim_control/time_integrator.h"

/// 
/// This can integrate any system of equations if given the right solver class.
/// It can solve the equations in 1st or 2nd order accuracy in space and time.
///
class sim_control_NG :
  virtual public sim_control,
  virtual public setup_NG_grid
{
  public:
  sim_control_NG();  ///< Simple constructor
  virtual ~sim_control_NG(); ///< Destructor

  ///
  /// initialisation.
  ///
  /// This function calls a sequence of other functions to set up the grid
  /// and populate it with the initial conditions, and give it the appropriate
  /// boundary conditions.  It gets the simulation ready to start, and checks 
  /// that everything is ready to start before returning.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  virtual int Init(
      string,   ///< Name of input file.
      int,      ///< Type of File (1=ASCII, 2=FITS, 5=Silo, ...)
      int,      ///< Number of command-line arguments.
      string *, ///< Pointer to array of command-line arguments.
      vector<class GridBaseClass *> &  ///< grid pointers.
      );

  ///
  /// Time integration
  ///
  /// This is the main part of the code -- It does all the time integration
  /// until the stopping condition is reached and then returns.
  /// It calls a sequence of functions to advance the time by one timestep,
  /// all in a loop which runs until end-of-sim is reached.
  ///
  virtual int Time_Int(
      vector<class GridBaseClass *> &  ///< grid pointers.
      );

  ///
  /// finalise the simulation, clean up, delete data.
  /// This function finished the simulation gracefully (hopefully!).
  ///
  int Finalise(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );
  
  //---------------------------------------
  protected:

  ///
  /// Calculates total values of conserved quantities.
  ///
  virtual int initial_conserved_quantities(
      vector<class GridBaseClass *> &  ///< grid pointers.
      );

#ifdef CHECK_MAGP
  ///
  /// This is only for a test problem -- it checks the magnetic
  /// pressure on the full domain and outputs it to screen
  ///
  void calculate_magnetic_pressure(
      vector<class GridBaseClass *> &  ///< grid pointers.
      );
#endif // CHECK_MAGP

#ifdef BLAST_WAVE_CHECK
  ///
  /// If running a spherical blast wave, calculate the shock
  /// position and output to screen.
  ///
  void calculate_blastwave_radius(
      vector<class GridBaseClass *> &  ///<  grid pointers.
      );
#endif // BLAST_WAVE_CHECK


  ///
  /// Advance a grid by one time step for a given level in a NG grid.
  /// This is a recursive function that calls itself on the next finer level
  /// if it exists.  This is a wrapper function that calls either the
  /// first-order or second-order update.
  ///
  virtual double advance_time(
      const int, ///< level in NG grid.
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// First-order-accurate time integration for a step on level l.
  /// Returns the sum of delta-t for the timestep just completed and
  /// the step to come.
  ///
  virtual double advance_step_OA1(
      const int ///< level in NG grid.
      );

  ///
  /// Second-order-accurate time integration for a step on level l.
  /// Returns the sum of delta-t for the timestep just completed, and
  /// the step to come.
  ///
  virtual double advance_step_OA2(
      const int ///< level in NG grid.
      );

  ///
  /// Checks Total energy relative to initial value, and prints a
  /// message if not.
  ///
  virtual int check_energy_cons(
      vector<class GridBaseClass *> &  ///< grid pointers.
      );

  

  ///
  /// Run through all radiation sources and calculate column
  /// densities through each grid for each one.  This exchanges data
  /// between grid levels so that the rays traverse all levels.
  ///
  virtual int RT_all_sources_levels(
      class SimParams &  ///< simulation parameters
      );

  ///
  /// Run through all radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell
  /// are stored in extra_data[i] for each cell.
  ///
  //virtual int RT_all_sources(
  //    class SimParams &,      ///< simulation parameters
  //    class GridBaseClass *,  ///< grid to trace rays on.
  //    const int               ///< level of NG grid.
  //    );

  ///
  /// Run through all on-grid radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell
  /// are stored in extra_data[i] for each cell.
  ///
  virtual int do_ongrid_raytracing(
      class SimParams &,      ///< simulation parameters
      class GridBaseClass *,  ///< grid to trace rays on.
      const int               ///< level of NG grid.
      );

  ///
  /// Run through all off-grid radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell
  /// are stored in extra_data[i] for each cell.
  ///
  virtual int do_offgrid_raytracing(
      class SimParams &,      ///< simulation parameters
      class GridBaseClass *,  ///< grid to trace rays on.
      const int               ///< level of NG grid.
      );

}; // sim_control_NG
   


#endif // if not SIM_CONTROL_NG_H
