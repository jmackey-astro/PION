/// \file sim_control_nested.cpp
/// 
/// \brief Simulation Control Class for Nested Grids.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for
/// the nested-grid simulation control class.  This is built on top
/// of the control class for uniform grids, and so doesn't add too
/// much, just the moving up and down between levels.
/// 
/// Modifications:
/// - 2018.05.03 JM: Started on nested grid simulation control.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "sim_control/sim_control_nested.h"

//#include "microphysics/microphysics_base.h"
//#include "raytracing/raytracer_SC.h"
//#include "dataIO/dataio_base.h"
//#ifdef SILO
//#include "dataIO/dataio_silo.h"
//#endif // if SILO
//#ifdef FITS
//#include "dataIO/dataio_fits.h"
//#endif // if FITS

//#include "spatial_solvers/solver_eqn_hydro_adi.h"
//#include "spatial_solvers/solver_eqn_mhd_adi.h"

//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <sys/time.h>
//#include <time.h>
//#include <climits>
using namespace std;

#define NEST_INT_TEST


// ##################################################################
// ##################################################################


sim_control_nestedgrid::sim_control_nestedgrid()
{
#ifdef TESTING
  cout << "(sim_control_nestedgrid::Constructor)\n";
#endif
}



// ##################################################################
// ##################################################################


sim_control_nestedgrid::~sim_control_nestedgrid()
{
#ifdef TESTING
  cout << "(sim_control_nestedgrid::Destructor)\n";
#endif
}




// ##################################################################
// ##################################################################




int sim_control_nestedgrid::Time_Int(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control_nestedgrid::Time_Int) STARTING TIME INTEGRATION\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  clk.start_timer("Time_Int"); double tsf=0;

  // make sure all levels start at the same time.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }


  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid[0]);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid);
#endif
    //
    // ----------------------------------------------------------------
    // Update RT sources and boundaries.
    //
    for (int l=0; l<SimPM.grid_nlevels; l++) {
      err = update_evolving_RT_sources(SimPM,grid[l]->RT);
      rep.errorTest("nested TIME_INT::update_evolving_RT_sources error",0,err);

      //cout <<"updating external boundaries for level "<<l<<"\n";
      err += TimeUpdateExternalBCs(SimPM, grid[l], l, SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
    }
    rep.errorTest("sim_control_nestedgrid: error from bounday update",0,err);
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
      //cout <<"updating internal boundaries for level "<<l<<"\n";
      err += TimeUpdateInternalBCs(SimPM, grid[l], l, SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
    }
    rep.errorTest("sim_control_nestedgrid: error from bounday update",0,err);
    // ----------------------------------------------------------------


    //
    // Get timestep on each level
    //
    int scale = 1;
    double mindt = 1.0e99;
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
      spatial_solver->set_dx(SimPM.levels[l].dx);
      //cout <<"dx="<<SimPM.levels[l].dx<<"\n";
      err += calculate_timestep(SimPM, grid[l],spatial_solver);
      rep.errorTest("TIME_INT::calc_timestep()",0,err);
      mindt = std::min(mindt, SimPM.dt/scale);
      //cout <<"level "<<l<<" got dt="<<SimPM.dt<<" and "<<SimPM.dt/scale <<"\n";
      SimPM.levels[l].dt = SimPM.dt;
      scale *= 2;
    }
    // make sure all levels use the same step (scaled by factors of 2).
    scale = 1;
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
      //cout <<"level "<<l<<", orig dt="<<SimPM.levels[l].dt;
      SimPM.levels[l].dt = mindt*scale;
      scale *= 2;
      //cout <<", new dt="<<SimPM.levels[l].dt<<"\n";
    }

    //clk.start_timer("advance_time");
    //
    // Use a recursive algorithm to update the coarsest level.  This
    // function also updates the next level twice, by calling itself
    // for the finer level, and so on.
    //
    advance_time(0);
    SimPM.simtime = SimPM.levels[0].simtime;

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.levels[0].dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
    tsf=clk.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err+= output_data(grid);
    rep.errorTest("TIME_INT::output_data()",0,err);
    err+= check_eosim();
    rep.errorTest("TIME_INT::check_eosim()",0,err);
  }

  cout <<"(sim_control_nestedgrid::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

  tsf=clk.time_so_far("Time_Int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep<<" wall-time=";
  cout <<tsf<<" time/step="<<tsf/static_cast<double>(SimPM.timestep)<<"\n";
  cout <<"STEPS "<<SimPM.timestep;
  cout.setf( ios_base::scientific );
  cout.precision(6);
  cout <<"\t"<<tsf<<"\t"<<tsf/static_cast<double>(SimPM.timestep);
  cout <<"\t"<<static_cast<double>(SimPM.timestep*SimPM.Ncell)/tsf<<"\n";
  cout <<"------------------------------------------------------------\n";

  return(0);
}



// ##################################################################
// ##################################################################





#ifdef CHECK_MAGP
///
/// This is only for a test problem -- it checks the magnetic
/// pressure on the full domain and outputs it to screen
///
void sim_control_nestedgrid::calculate_magnetic_pressure(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the total magnetic pressure on the domain, normalised to the
  // initial value.
  //
  double magp=0.0, cellvol=0.0;
  static double init_magp=-1.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    spatial_solver->set_dx(SimPM.levels[l].dx);
    
    cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd) 
        magp += (spatial_solver->Ptot(c->P,0.0) - c->P[PG]) * spatial_solver->CellVolume(c);
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }
  if (init_magp<0) init_magp = magp;
  cout <<SimPM.simtime<<"\t"<<magp/init_magp<<"\t"<<magp<<"\n";
  return;
}
#endif // CHECK_MAGP



// ##################################################################
// ##################################################################



#ifdef BLAST_WAVE_CHECK
///
/// If running a 1D spherical blast wave, calculate the shock position
/// and output to screen.
///
void sim_control_nestedgrid::calculate_blastwave_radius(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the blast wave outer shock position.
  // If a nested grid, start on the finest grid and work outwards
  //
  double shockpos=0.0;
  static double old_pos=0.0;
  bool shock_found = false;
  //  static double last_dt=0.0;
  for (int l=SimPM.grid_nlevels-1; l>=0; l++) {
    spatial_solver->set_dx(SimPM.levels[l].dx);

    if (shock_found) continue;
    cell *c=grid->LastPt();
    if (fabs(c->P[VX])>=1.0e4) {
      cout<<"level "<<l<<" does not contain shock.\n";
    }
    else {
      do {
        c = grid->NextPt(c,RNsph);
        //cout <<c->id<<", vx="<<c->P[VX]<<"\n";
      } while ( c!=0 && fabs(c->P[VX])<1.0e4);
      if (c && (c->P[VX] >= 1.0e4)) {
        shockpos = CI.get_dpos(c,Rsph);
        shock_found=true;
      }
    }
  }
  if (pconst.equalD(old_pos,0.0))
    old_pos = shockpos;
  cout <<SimPM.simtime<<"\t"<<shockpos;
  //cout <<"\t"<<(shockpos-old_pos)/(SimPM.dt+TINYVALUE);
  cout <<"\n";
  old_pos=shockpos;
  return;
}
#endif // BLAST_WAVE_CHECK



// ##################################################################
// ##################################################################



double sim_control_nestedgrid::advance_time(
      const int l       ///< level to advance.
      )
{
#ifdef TESTING
  cout <<"advance_time, level="<<l<<", starting.\n";
#endif
  double step=0.0;
  if (SimPM.tmOOA==1) {
    step = advance_step_OA1(l);
  }
  else if (SimPM.tmOOA==2) {
#ifdef NEST_INT_TEST
    //cout <<"Calling advance_step_OA2: level "<<l<<"\n";
#endif
    step = advance_step_OA2(l);
  }
  return step;
}




// ##################################################################
// ##################################################################



double sim_control_nestedgrid::advance_step_OA1(
      const int l       ///< level to advance.
      )
{
#ifdef TESTING
  cout <<"advance_step_OA1, level="<<l<<", starting.\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  class GridBaseClass *grid = SimPM.levels[l].grid;

  // take the first finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA1(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;

  // now calculate dU, the change in conserved variables on this grid,
  // for this step.
  spatial_solver->set_dx(SimPM.levels[l].dx);
  spatial_solver->Setdt(SimPM.levels[l].dt);
  // May need to do raytracing
  if (!FVI_need_column_densities_4dt && grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid->RT);
    rep.errorTest("scn::advance_time: calc_rt_cols()",0,err);
  }
  err += calc_microphysics_dU(SimPM.levels[l].dt, grid);
  err += calc_dynamics_dU(SimPM.levels[l].dt,OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(SimPM.levels[l].dt,OA1, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA1: calc_x_dU",0,err);

  // take the second finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA1(l+1);
  }

  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(SimPM.levels[l].dt,OA1,OA1, grid);
  rep.errorTest("scn::advance_step_OA1: state-vec update",0,err);  

  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }

  //
  // update internal and external boundaries.
  //
  err += TimeUpdateInternalBCs(SimPM, grid, l, SimPM.simtime, OA2, OA2);
  err += TimeUpdateExternalBCs(SimPM, grid, l, SimPM.simtime, OA2, OA2);

  // Now calculate next timestep: function stores dt in SimPM.dt
  err += calculate_timestep(SimPM, grid,spatial_solver);
  rep.errorTest("scn::advance_step_OA1: calc_timestep",0,err);

  // make sure step is not more than half of the coarser grid step.
  if (l>0) {
    SimPM.levels[l].dt = min(SimPM.dt, 0.5*SimPM.levels[l-1].dt);
  }
  else {
    SimPM.levels[l].dt = SimPM.dt;
  }

#ifdef TESTING
  cout <<"advance_step_OA1, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step<<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}




// ##################################################################
// ##################################################################




double sim_control_nestedgrid::advance_step_OA2(
      const int l       ///< level to advance.
      )
{
#ifdef NEST_INT_TEST
  cout <<"advance_step_OA2, level="<<l<<", starting. ";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step<<"\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  class GridBaseClass *grid = SimPM.levels[l].grid;

  // take the first finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA2(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;

  spatial_solver->set_dx(SimPM.levels[l].dx);
  spatial_solver->Setdt(SimPM.levels[l].dt);
  // May need to do raytracing
  if (!FVI_need_column_densities_4dt && grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid->RT);
    rep.errorTest("scn::advance_time: calc_rt_cols()",0,err);
  }

  //
  // now calculate dU, the change in conserved variables on this grid
  // for the half step of the 2nd order step.
  //
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" half step, start calc_microphysics_dU\n";
#endif
  err += calc_microphysics_dU(SimPM.levels[l].dt, grid);
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" half step, start calc_dynamics_dU\n";
#endif
  err += calc_dynamics_dU(SimPM.levels[l].dt, OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(SimPM.levels[l].dt, OA1, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA1",0,err);
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" half step, done with dU, grid_update_state_vector\n";
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt,OA1,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA1",0,err);  
  // Update boundary data.
  err += TimeUpdateInternalBCs(SimPM, grid, l, SimPM.simtime, OA1, OA2);
  err += TimeUpdateExternalBCs(SimPM, grid, l, SimPM.simtime, OA1, OA2);
  rep.errorTest("scn::advance_step_OA2: bounday update OA1",0,err);

  //
  // Now calculate dU for the full step (OA2)
  //
  if (grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid->RT);
    rep.errorTest("scn::advance_time: calc_rt_cols() OA2",0,err);
  }
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" full step, start calc_microphysics_dU\n";
#endif
  err += calc_microphysics_dU(SimPM.levels[l].dt, grid);
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" full step, start calc_dynamics_dU\n";
#endif
  err += calc_dynamics_dU(SimPM.levels[l].dt, OA2, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(SimPM.levels[l].dt, OA2, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA2",0,err);

  // take the second finer grid step, if there is a finer grid.
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" full step, call 2nd l+1 update\n";
#endif
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA2(l+1);
  }

  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
#ifdef NEST_INT_TEST
  cout <<"l="<<l<<" full step, done with dU, grid_update_state_vector\n";
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt,OA2,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA2",0,err);  

  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }

  //
  // update internal and external boundaries.
  //
  err += TimeUpdateInternalBCs(SimPM, grid, l, SimPM.simtime, OA2, OA2);
  err += TimeUpdateExternalBCs(SimPM, grid, l, SimPM.simtime, OA2, OA2);

  // Now calculate next timestep: function stores dt in SimPM.dt
  err += calculate_timestep(SimPM, grid,spatial_solver);
  rep.errorTest("scn::advance_step_OA2: calc_timestep",0,err);

  // make sure step is not more than half of the coarser grid step.
  if (l>0) {
    SimPM.levels[l].dt = min(SimPM.dt, 0.5*SimPM.levels[l-1].dt);
  }
  else {
    SimPM.levels[l].dt = SimPM.dt;
  }

#ifdef NEST_INT_TEST
  cout <<"advance_step_OA2, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step<<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}



// ##################################################################
// ##################################################################



