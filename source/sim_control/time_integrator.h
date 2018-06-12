/// \file time_integrator.h
/// \brief time integration routines for sim_control class.
/// \author Jonathan Mackey
/// 
/// Time integration class definitions, for the 1st/2nd order Finite
/// Volume Solver following Falle, Komissarov, \& Joarder
/// (1998,MNRAS,297,265).
/// 
/// Modifications:
/// - 2018.05.11 JM: changed time-integration routines to be their
///    own class.

#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/setup_fixed_grid.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "equations/eqns_base.h"
#include "sim_control/sim_init.h"
#include "sim_control/update_boundaries.h"
#include "sim_control/calc_timestep.h"




// ##################################################################
// ##################################################################

///
/// Class to integrate a finite-volume grid of data forward in time,
/// for PION simulations of hydro/MHD.
///
class time_integrator :
  virtual public sim_init,
  virtual public calc_timestep
{
  public:
  time_integrator();
  ~time_integrator();
   

  ///
  /// Advance a grid by one time step.  This is the main time
  /// integration function, and calls calculate_timestep, followed by
  /// either a first-order or second-order time update.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  virtual int advance_time(
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This performs a first-order-accurate (in time) timestep for
  /// dynamics, microphysics, thermal conduction.
  ///
  /// If the order-of-accuracy parameter is OA1 then it assumes this
  /// is a full step and updates both P[] and Ph[].
  /// If it is OA2, then the functions assumes this is the half-step
  /// as part of a full second-order step, so it only updates Ph[].
  /// It advances by an interval of dt regardless of OA1 or OA2, so 
  /// if this is a half step 0.5*dt is passed to the function.
  ///
  int first_order_update(
      const double,  ///< dt, time interval to advance by.
      const int,     ///< time order of accuracy OA1/OA2.
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This performs a second-order-accurate (in time) timestep for
  /// dynamics, microphysics, thermal conduction.
  /// This performs the second part of a full second-order step, so
  /// the half-step must have been already called before this one.
  ///
  int second_order_update(
      const double, ///< dt, time interval to advance by.
      const int,    ///< time order of accuracy (must be OA2).
      class GridBaseClass * ///< grid pointer
      );
  
  ///
  /// This function does some checking on radiation sources to see
  /// what microphysics update to call, then calls one of 
  /// calc_RT_microphysics_dU() or
  /// calc_microphysics_dU().
  ///
  int calc_microphysics_dU(
      const double, ///< dt, timestep to integrate MP spatial_solvers.
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This calculates the change in internal energy and ion fractions
  /// for a timestep dt, by integrating the microphysics equations
  /// for the full timestep, storing the result in a temporary array,
  /// and differencing the initial and final states.
  /// This version is for microphysics integrations where there are
  /// radiation sources involved.
  ///
  int calc_RT_microphysics_dU(
      const double,   ///< dt, timestep to integrate
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This calculates the change in internal energy and ion fractions
  /// for a timestep dt, by integrating the microphysics equations
  /// for the full timestep, storing the result in a temporary array,
  /// and differencing the initial and final states.
  /// This version is for microphysics integrations where there are
  /// no radiation sources (e.g. pure heating+cooling, or collisional
  /// processes only).
  ///
  int calc_noRT_microphysics_dU(
      const double, ///< dt, timestep to integrate
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This calculates the change in the state vector for each point
  /// due to the dynamics, for a timestep dt, using either 1st or 
  /// 2nd order accuracy in space.
  /// It calls spatial_solver->preprocess_data(), then set_dynamics_dU(), and
  /// finally spatial_solver->PostProcess_dU().
  /// set_dynamics_dU() is the function that used to be called
  /// calc_dU().
  ///
  int calc_dynamics_dU(
      const double, ///< dt, timestep to integrate
      const int,    ///< spatial order of accuracy for update.
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This function used to be called calc_dU -- for every column of
  /// simulation data on the grid this calls dU_Column() to get the
  /// changes in the state vectors arising from the hydrodynamics
  /// over the timestep dt.  It uses the requested spatial order of
  /// accuracy.  This function also loops over all directions on the 
  /// grid that are active.
  ///
  int set_dynamics_dU(
      const double,    ///< dt, timestep for this calculation
      const int,       ///< space OOA for this calculation
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// Calculate dU, rate of change of conserved variables, in a 1D
  /// column of grid cells, according to the fluid dynamics equations.
  ///
  /// This runs through every cell in a 1D column in turn, and calculates the flux
  /// between the cell in question and its neighbour to the right, by obtaining
  /// an interface flux.
  /// 
  /// It then calculates dU for each cell according to the exact formula (if the 
  /// flux calculation were exact) given by Toro eq.5.76\n
  /// \f$ U_i^{n+1}-U_i^n =dU = \frac{\Delta t}{\Delta x}(F_{i-\frac{1}{2}} -F_{i+\frac{1}{2}}) \f$.
  ///
  int dynamics_dU_column(const class cell *, ///< starting point for column.
      const enum direction, ///< direction to traverse column in. 
      const enum direction, ///< opposite direction.
      const double,    ///< dt, timestep for this calculation
#ifdef TESTING
      const int,       ///< Time Order of accuracy to use.
#endif
      const int,        ///< Spatial Order of accuracy to use.
      class GridBaseClass * ///< grid pointer
      );

  ///
  /// This function takes the contents of each cell->dU[] vector and
  /// updates Ph[] the changes.  If we are on the full-step then it
  /// also updates P[] so that it and Ph[] are identical.
  ///
  int grid_update_state_vector(
      const double ,  ///< dt, timestep
      const int,      ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int,       ///< Full order of accuracy of simulation
      class GridBaseClass * ///< grid pointer
      );

};

#endif // TIME_INTEGRATOR_H
