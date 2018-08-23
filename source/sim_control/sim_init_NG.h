/// \file sim_init_NG.h
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for sim_init_NG, which sets up a PION simulation
/// and gets everything ready to run.
///
/// Modifications:\n
/// - 2018.05.11 JM: moved code from sim_control.h
///

#ifndef SIM_INIT_NESTED_H
#define SIM_INIT_NESTED_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/assign_update_bcs_NG.h"
#include "grid/setup_NG_grid.h"
#include "grid/grid_base_class.h"
#include "sim_control/sim_init.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"



// ##################################################################
// ##################################################################

///
/// Class to set up a simulation so that everything is ready to run.
/// Inherits from setup_NG_grid and update_boundaries_NG.
///
class sim_init_NG
  :
  virtual public sim_init,
  virtual public setup_NG_grid
{
  public:
  sim_init_NG();
  ~sim_init_NG();
   
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


  protected:

  /// function to setup data-I/O class.
  virtual void setup_dataio_class(
      const int  ///< type of I/O: 1=text,2=fits,5=silo
      );

  ///
  /// Calculates total values of conserved quantities.
  ///
  int initial_conserved_quantities(
      vector<class GridBaseClass *> &  ///< address of vector of grid pointers.
      );


};




#endif // SIM_INIT_NESTED_H

