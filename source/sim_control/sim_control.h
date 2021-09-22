/// \file sim_control.h
///
/// \brief Declares grid parameter class, and grid methods classes.
///
/// \author Jonathan Mackey
///
/// IntegratorBaseFV is an abstract base class for finite volume grids.
///
/// IntUniformFV is a 1st/2nd order accurate finite volume solver, modelled on
/// the solver presented in Falle, Komissarov, \& Joarder (1998) MNRAS, 297,
/// 265.
///
/// Modifications :\n
/// - 2007-07-12 Started to move to new structure.
/// - 2007-07-13 New Class structure implemented.
/// - 2007-07-24 Added tracer variables.
/// - 2007-10-11 Updated Documentation.
/// - 2007-10-26 New data i/o class.
/// - 2010-04-21 JM: removed parallel output_data() function b/c not needed
/// anymore.
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.10.04 JM: Moved the field-loop magnetic pressure output to
///     an ifdeffed function.  Added a new function to calculate the
///     1D blast wave radius in spherical coordinates (also ifdeffed).
/// - 2010.11.15 JM: Modified update_dynamics() so it has pre- and
///   post-processing functions before and after calc_dU().
///   Added routine to calculate the H-correction eta values for
///   pre-processing.
/// - 2011.01.03 JM: Moved preprocess_data() and calc_Hcorrection from
///   gridMethods.cc to the base FV solver.
/// - 2011.03.21 JM: moved cell setup_extra_data() to its own function, to save
///    on copy-paste for the parallel version.
///    Rewrote setup_raytracing() and atomised update_microphysics() so there
///    are now a few new functions to deal with the different kinds of radiative
///    transfer we want to do.
/// - 2011.04.06 JM: Added thermal-conduction timestep limiting and flux
/// calculation.
/// - 2011.10.22 JM: Added lists of heating and ionising sources as class data.
/// - 2012.01.16 JM: Added setup_evolving_RT_sources() and
///    update_evolving_RT_sources() for stellar evolution models.
/// - 2012.08.16 JM: Added functions for new 2nd order time update.
/// - 2013.02.19 JM: Made some initialisation functions public.
/// - 2015.01.08 JM: Made some initialisation functions public. Moved
///    grid pointer from global to main, so it now needs to be passed
///    around in some functions.
/// - 2015.01.26 JM: changed name to sim_control.cpp, and class names
///    to sim_control_XX, to better reflect what the class does!
///    Moved MPI-parallelised class to its own new file.
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global
/// - 2018.05.** JM: moved most functions to new classes for calculating
/// timestep,
///    updating boundaries, and time integration.

#ifndef SIM_CONTROL_H
#define SIM_CONTROL_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "dataIO/dataio_base.h"
#include "grid/grid_base_class.h"
#include "grid/setup_fixed_grid.h"
#include "grid/uniform_grid.h"
#include "sim_control/calc_timestep.h"
#include "sim_control/sim_init.h"
#include "sim_control/time_integrator.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// Simulation
///
class sim_control : virtual public time_integrator {
public:
  sim_control();           ///< Simple constructor, initialises value.
  virtual ~sim_control();  ///< Deletes any dynamic memory, if not already
                           ///< done.

  ///
  /// Time integration
  ///
  /// This is the main part of the code -- It does the time integration
  /// until the stopping condition is reached and then returns.
  /// It calls a sequence of functions to advance the time by one timestep
  /// repeatedly until end-of-sim is reached.
  ///
  virtual int Time_Int(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
  );

  ///
  /// finalise the simulation, clean up, delete data.
  /// This function finished the simulation gracefully (hopefully!).
  ///
  virtual int Finalise(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
  );

  //---------------------------------------
protected:
#ifdef CHECK_MAGP
  ///
  /// This is only for a test problem -- it checks the magnetic
  /// pressure on the full domain and outputs it to screen
  ///
  virtual void calculate_magnetic_pressure(
      class GridBaseClass *  ///< address of grid pointer.
  );
#endif  // CHECK_MAGP

#ifdef BLAST_WAVE_CHECK
  ///
  /// If running a spherical blast wave, calculate the shock
  /// position and output to screen.
  ///
  virtual void calculate_blastwave_radius(
      class GridBaseClass *  ///< address of grid pointer.
  );
#endif  // BLAST_WAVE_CHECK

  /// Check if sim should stop.
  ///
  /// For shock tube problems, I stop the simulation whenever a disturbance
  /// reaches the edge of the domain.  Most problems set the finishtime to
  /// some value, and the end-of-sim criteria is if the current simtime has
  /// reached finishtime or not. \retval 0 success \retval 1 failure
  ///
  int check_eosim();  // Checks for end of simulation.

  ///
  /// Checks Total energy relative to initial value, and prints a
  /// message if not.
  ///
  int check_energy_cons(
      class GridBaseClass *  ///< address of vector of grid pointers.
  );

};  // sim_control

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif  // if not SIM_CONTROL_H
