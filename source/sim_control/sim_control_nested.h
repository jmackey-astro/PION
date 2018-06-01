/// \file sim_control_nested.h
/// 
/// \brief Declares grid parameter class, and grid methods classes.
/// 
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.05.04 JM: worked on class.

#ifndef SIM_CONTROL_NESTED_H
#define SIM_CONTROL_NESTED_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"
#include "sim_control/sim_control.h" 
#include "sim_control/sim_init_nested.h"
#include "sim_control/update_boundaries_nested.h"
#include "sim_control/calc_timestep.h"
#include "sim_control/time_integrator.h"

/// 
/// This can integrate any system of equations if given the right solver class.
/// It can solve the equations in 1st or 2nd order accuracy in space and time.
///
class sim_control_nestedgrid : virtual public sim_control, virtual public time_integrator
{
  public:
  sim_control_nestedgrid();  ///< Simple constructor
  virtual ~sim_control_nestedgrid(); ///< Destructor

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




}; // sim_control_nestedgrid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SIM_CONTROL_NESTED_H
