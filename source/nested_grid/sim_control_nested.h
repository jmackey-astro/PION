/// \file sim_control_nested.h
/// 
/// \brief Declares grid parameter class, and grid methods classes.
/// 
/// \author Jonathan Mackey
/// 
/// IntegratorBaseFV is an abstract base class for finite volume grids.
/// 
/// IntUniformFV is a 1st/2nd order accurate finite volume solver, modelled on 
/// the solver presented in Falle, Komissarov, \& Joarder (1998) MNRAS, 297, 265.
/// 
/// Modifications :\n

#ifndef SIM_CONTROL_NESTED_H
#define SIM_CONTROL_NESTED_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "grid/uniform_grid.h"
#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"
#include "setup_fixed_grid.h"

/// The simplest finite volume grid -- a uniform grid with cells that
/// are cube-shaped in the chosen coordinates.
/// 
/// This can integrate any system of equations if given the right solver class.
/// It can solve the equations in 1st or 2nd order accuracy in space and time.
///
class sim_control : virtual public setup_fixed_grid
{
  public:
  sim_control();  ///< Simple constructor, initialises value.
  virtual ~sim_control(); ///< Deletes any dynamic memory, if not already done.

  ///
  /// Function to print command-line options for PION.
  ///
  void print_command_line_options(int, char **);

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
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
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
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
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

#ifdef BLAST_WAVE_CHECK
  ///
  /// If running a spherical blast wave, calculate the shock
  /// position and output to screen.
  ///
  void calculate_blastwave_radius(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );
#endif // BLAST_WAVE_CHECK


  ///
  /// For a nested grid we need to start at the finest level and work
  /// upwards through the levels.  This is a wrapper function that
  /// calls the uniform_grid version of advance_time on each level.
  ///
  int advance_time(
      vector<class GridBaseClass *> &, ///< grid pointer
      vector<class RayTracingBase *> & ///< raytracer for this grid.
      );

}; // sim_control_nestedgrid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SIM_CONTROL_NESTED_H
