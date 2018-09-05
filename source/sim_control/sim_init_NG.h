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
   

};




#endif // SIM_INIT_NESTED_H

