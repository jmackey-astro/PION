/// \file calc_timestep.h
/// \brief routines for calculating the timestep on a grid.
/// \author Jonathan Mackey
/// \date 2018.05.10
/// 
/// Description:\n
/// Has a set of routines for calculating the timestep for fluid
/// dynamics simulations in PION.
/// 
/// Modifications:
/// - 2018.05.10 JM: moved from sim_control into its own class that
///    inherits from setup_fixed_grid.

#ifndef CALC_TIMESTEP_H
#define CALC_TIMESTEP_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "grid/setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"

///
/// This class operates on finite-volume computational grid for
/// fluid dynamics, calculating the allowed timesteps based on
/// various criteria such as local flow velocities, chemical
/// timescales, and heat conduction timescales.
///
/// It inherits from setup_fixed_grid, and is inherited by
/// sim_control.
///
class calc_timestep :
  virtual public setup_fixed_grid
{
  public:

  calc_timestep();
  
  ~calc_timestep();

  ///
  /// Calculate the appropriate timestep.
  /// 
  /// For a uniform grid, all cells have the same timestep equal to the minimum
  /// of all calculated steps.  This function calls two functions, one to get 
  /// the microphysics timestep (if needed), and another to get the dynamics
  /// timestep.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  virtual int calculate_timestep(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *, ///< pointer to grid.
      class FV_solver_base * ///< solver/equations class
      );

  protected:

  ///
  /// Calculate the microphysics timestep, based on heating/cooling and reaction
  /// rates.  Returns the minimum timestep of the local grid (negative if error).
  /// 
  double calc_microphysics_dt(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass * ///< pointer to grid.
      );

  ///
  /// Old microphysics timescales calculation with no radiation field.
  ///
  double get_mp_timescales_no_radiation(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass * 
      );

  ///
  /// New microphysics timescales calculation with pre-calculated radiation field.
  ///
  double get_mp_timescales_with_radiation(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass * 
      );

  ///
  /// Run through all diffuse and direct radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell are stored
  /// in extra_data[i] for each cell.
  ///
  int calculate_raytracing_column_densities(
      class SimParams &,      ///< pointer to simulation parameters
      class RayTracingBase * ///< raytracer for this grid.
      );

  ///
  /// Calculate the dynamics timestep, based on the Courant condition that
  /// the fastest signals cannot cross a full cell in a single step.  Returns
  /// the minimum timestep on the local grid, or negative if an error occurs.
  ///
  double calc_dynamics_dt(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *, 
      class FV_solver_base * ///< solver/equations class
      );


#ifdef THERMAL_CONDUCTION
  ///
  /// If doing thermal conduction, this calculates the max. timestep we can use
  /// without conduction changing the gas temperature by more than 30%.
  /// It also calls the solver function set_thermal_conduction_Edot() so it 
  /// knows what the flux in and out of cells is.  This Edot value is multiplied
  /// by the timestep dt in spatial_solver->preprocess_data().
  ///
  double calc_conduction_dt_and_Edot(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *, ///< pointer to grid.
      class FV_solver_base * ///< solver/equations class
      );
#endif // THERMAL CONDUCTION

  ///
  /// Limits the timestep based on a number of criteria.  Makes sure we don't
  /// overshoot the finish-time, or the next output-time, and that we don't 
  /// increase the timestep by a large factor from one step to the next (this
  /// can affect stability).
  ///
  void timestep_checking_and_limiting(
      class SimParams &      ///< pointer to simulation parameters
      );

  // ----------------------------------------------------------------
  // ---- Data Members ----
  // ----------------------------------------------------------------
  protected:
  // ----------------------------------------------------------------
  // ---- Data Members ----
  // ----------------------------------------------------------------

};

#endif // CALC_TIMESTEP_H
