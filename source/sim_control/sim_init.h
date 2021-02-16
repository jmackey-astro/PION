/// \file sim_init.h
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for sim_init, which sets up a PION simulation
/// and gets everything ready to run.
///
/// Modifications:\n
/// - 2018.05.11 JM: moved code from sim_control.h
///

#ifndef SIM_INIT_H
#define SIM_INIT_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_base.h"
#include "grid/grid_base_class.h"
#include "grid/setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"

// ##################################################################
// ##################################################################

///
/// Class to set up a simulation so that everything is ready to run.
/// Inherits from setup_fixed_grid.
///
class sim_init : virtual public setup_fixed_grid {
public:
  sim_init();
  ~sim_init();

  ///
  /// Function to print command-line options for PION.
  ///
  void print_command_line_options(int, char**);

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
      string*,  ///< Pointer to array of command-line arguments.
      vector<class GridBaseClass*>&  ///< grid pointers.
  );

  ///
  /// Set the maximum runtime to a new value. Should be set in main()
  ///
  void set_max_walltime(double  ///< New Max. runtime in seconde.
  );
  ///
  /// Get the maximum runtime in seconds.
  ///
  double get_max_walltime();

  ///
  /// information about the simulation
  ///
  class SimParams SimPM;

protected:
  //---------------------------------------
  //---------------------------------------
  // Class Member data:
  //
#ifdef TEST_CONSERVATION
  double initERG, initMMX, initMMY, initMMZ, initMASS;
  double nowERG, nowMMX, nowMMY, nowMMZ, nowMASS;
#endif

  ///
  /// Max. walltime to run for, in seconds, after which we save
  /// data and finish.
  ///
  double max_walltime;

  //---------------------------------------
  //---------------------------------------
  ///
  /// See if any command-line arguments should override those
  /// specified in the IC file, and if so, reset the parameters.
  ///
  int override_params(
      int,     ///< Number of command-line arguments.
      string*  ///< Pointer to array of command-line arguments.
  );

  ///
  /// Save simulation snapshot to file if required.
  ///
  /// This checks if I want to output data in this timestep, then
  /// checks what format to write in, and calls the appropriate
  /// function to write the data.
  ///
  virtual int output_data(vector<class GridBaseClass*>&  ///< grid pointers.
  );

  ///
  /// Calculates total values of conserved quantities.
  ///
  int initial_conserved_quantities(class GridBaseClass*  ///< grid pointers.
  );

  ///
  /// Run through all radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell
  /// are stored in extra_data[i] for each cell.
  ///
  virtual int RT_all_sources(
      class SimParams&,      ///< simulation parameters
      class GridBaseClass*,  ///< grid to trace rays on.
      const int              ///< level of NG grid.
  );
};

#endif  // SIM_INIT_H
