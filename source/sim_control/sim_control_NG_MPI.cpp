/// \file sim_control_NG_MPI.cpp
/// 
/// \brief Parallel Grid Methods Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for 
/// the "sim_control_NG_MPI" class, which is a 1st/2nd order Finite
/// Volume Solver according to the method outlined in Falle, 
/// Komissarov, \& Joarder (1998,MNRAS,297,265).
/// It includes extensions for a nested grid, with the grid on each
/// level decomposed into blocks with communication via MPI.
/// 
/// Modifications:
/// - 2018.09.04 JM: Modified from sim_control_MPI.cpp. 

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "tools/mem_manage.h"

#include "constants.h"

#include "decomposition/MCMD_control.h"
#include "sim_control/sim_control_NG_MPI.h"
#include "raytracing/raytracer_SC.h"


#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#ifdef PARALLEL


//#define TEST_INT

// ##################################################################
// ##################################################################


sim_control_NG_MPI::sim_control_NG_MPI()
{
#ifdef TESTING
  cout <<"sim_control_NG_MPI constructor.\n";
#endif
}



// ##################################################################
// ##################################################################


sim_control_NG_MPI::~sim_control_NG_MPI()
{
#ifdef TESTING    
  cout <<"sim_control_NG_MPI destructor.\n";
#endif
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
#ifdef TESTING
  cout <<"(sim_control_NG_MPI::init) Initialising grid: infile = "<<infile<<"\n";
#endif
  int err=0;

  //
  // Setup the MCMDcontrol class with rank and nproc.
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);

  // ----------------------------------------------------------------
  SimPM.typeofip=typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  if (!dataio->file_exists(infile))
    rep.error("infile doesn't exist!",infile);

  // ----------------------------------------------------------------
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("NG_MPI Init(): failed to read header",0,err);

  // Check if any commandline args override the file parameters.
  // ----------------------------------------------------------------
  err = override_params(narg, args);
  rep.errorTest("(NG_MPI INIT::override_params)",0,err);
  
  // setup the nested grid levels, and decompose the domain on each
  // level
  // ----------------------------------------------------------------
  setup_NG_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
  cout <<"NG_MPI Init: grid setup\n";
  err = setup_grid(grid,SimPM);
  cout <<"NG_MPI Init: grid setup finished\n";
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(NG_MPI INIT::setup_grid) error",0,err);

  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  // ----------------------------------------------------------------
  err = set_equations(SimPM);
  rep.errorTest("(NG_MPI INIT::set_equations)",0,err);
  spatial_solver->SetEOS(SimPM.gamma);

  // set up Microphysics, if needed.
  // ----------------------------------------------------------------
  err = setup_microphysics(SimPM);
  rep.errorTest("(NG_MPI INIT::setup_microphysics)",0,err);
  
  // assign data to the grid from snapshot file.
  // ----------------------------------------------------------------
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(NG_MPI INIT::assign_initial_data)",0,err);

  // ----------------------------------------------------------------
  // Set Ph[] = P[], and then implement the boundary conditions.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    cell *c = grid[l]->FirstPt();
    do {
      for(int v=0;v<SimPM.nvar;v++) c->Ph[v]=c->P[v];
    } while ((c=grid[l]->NextPt(c))!=0);

    if (SimPM.eqntype==EQGLM && SimPM.timestep==0) {
#ifdef TESTING
      cout <<"Initial state, zero-ing glm variable.\n";
#endif
      c = grid[l]->FirstPt(); do {
        c->P[SI] = c->Ph[SI] = 0.;
      } while ( (c=grid[l]->NextPt(c)) !=0);
    }
  } // loop over levels

  // ----------------------------------------------------------------
  err = boundary_conditions(SimPM, grid);
  rep.errorTest("(NG_MPI INIT::boundary_conditions) err!=0",0,err);

  // ----------------------------------------------------------------
  err += setup_raytracing(SimPM, grid);
  rep.errorTest("Failed to setup raytracer",0,err);

  // ----------------------------------------------------------------
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err = assign_boundary_data(SimPM, l, grid[l]);
    rep.errorTest("NG_MPI INIT::assign_boundary_data",0,err);
    COMM->barrier("level assign boundary data");
  }
  // ----------------------------------------------------------------


  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"NG_MPI updating external boundaries for level "<<l<<"\n";
    cout <<"@@@@@@@@@@@@  UPDATING EXTERNAL BOUNDARIES FOR LEVEL ";
    cout <<l<<"\n";
#endif
    err += TimeUpdateExternalBCs(SimPM,l,grid[l], spatial_solver,
                            SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"NG_MPI updating C2F boundaries for level "<<l<<"\n";
    cout <<"@@@@@@@@@@@@  UPDATING C2F BOUNDARIES FOR LEVEL ";
    cout <<l<<"\n";
#endif
    if (l<SimPM.grid_nlevels-1) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid[l],
                spatial_solver, l, grid[l]->BC_bd[i], 2,2);
        }
      }
    }
    if (l>0) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                      l,grid[l]->BC_bd[i],SimPM.levels[l].step);
        }
      }
    }
  }
  BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"NG_MPI updating external boundaries for level "<<l<<"\n";
    cout <<"@@@@@@@@@@@@  UPDATING EXTERNAL BOUNDARIES FOR LEVEL ";
    cout <<l<<"\n";
#endif
    err += TimeUpdateExternalBCs(SimPM,l,grid[l], spatial_solver,
                            SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------


  // ----------------------------------------------------------------
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TESTING
    cout <<"NG_MPI updating internal boundaries for level "<<l<<"\n";
    cout <<"@@@@@@@@@@@@  UPDATING INTERNAL BOUNDARIES FOR LEVEL ";
    cout <<l<<"\n";
#endif
    err += TimeUpdateInternalBCs(SimPM,l,grid[l], spatial_solver,
                            SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);


  //
  // If testing the code, this calculates the momentum and energy
  // on the domain.
  //
  // ----------------------------------------------------------------
  initial_conserved_quantities(grid);


  //
  // If using opfreq_time, set the next output time correctly.
  //
  // ----------------------------------------------------------------
  if (SimPM.op_criterion==1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!",
                SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
    double tmp = 
      ((SimPM.simtime/SimPM.opfreq_time)-
       floor(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
    SimPM.next_optime-= tmp;
  }

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  // ----------------------------------------------------------------
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    setup_dataio_class(SimPM, SimPM.typeofop);
    if (!dataio) rep.error("NG_MPI INIT:: dataio",SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

  if (SimPM.timestep==0) {
#ifdef TESTING
     cout << "(NG_MPI INIT) Writing initial data.\n";
#endif
     output_data(grid);
  }
  
  // ----------------------------------------------------------------
//#ifdef TESTING
  cell *c=0;
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
    //cout <<"LEVEL-ZERO-CHECK L="<<l<<"\n";
    c = (grid[l])->FirstPt_All();
    do {
      if (pconst.equalD(c->P[RO],0.0)) {
        cout <<"ZERO DATA IN CELL: ";
        CI.print_cell(c);
      }
    } while ( (c=(grid[l])->NextPt_All(c)) !=0 );
  }
//#endif // TESTING
  cout <<"-------------------------------------------------------\n";
  return(err);
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::Time_Int(
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      )
{
  cout <<"-------------------------------------------------------\n";
  cout <<"--- sim_control_NG_MPI:: STARTING TIME INTEGRATION. ---\n";
  cout <<"-------------------------------------------------------\n";
  int err=0;
  int log_freq=10;
  SimPM.maxtime=false;

  bool first_step=true;
  if (SimPM.timestep!=0) first_step=false;

  clk.start_timer("time_int"); double tsf=0.0;

  // make sure all levels start at the same time.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }

  // ----------------------------------------------------------------
  // update coarse-to-fine level boundaries
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"NG_MPI updating C2F boundaries for level "<<l<<"\n";
    cout <<l<<"\n";
#endif
    if (l<SimPM.grid_nlevels-1) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid[l],
                spatial_solver, l, grid[l]->BC_bd[i], 2,2);
        }
      }
    }
    if (l>0) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                      l,grid[l]->BC_bd[i],SimPM.levels[l].step);
        }
      }
    }
  }
  BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------
  // --------------------------------------------------------------
  // Update internal and external boundaries.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TEST_INT
    cout <<"updating external boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateExternalBCs(SimPM, l, grid[l], spatial_solver,
                  SimPM.levels[l].simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("sim_control_NG_MPI: external boundary",0,err);

  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TEST_INT
    cout <<"updating internal boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateInternalBCs(SimPM, l, grid[l], spatial_solver,
                  SimPM.levels[l].simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("sim_control_NG_MPI: internal boundary",0,err);
  // --------------------------------------------------------------
  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TESTING
    cout <<"NG_MPI updating F2C boundaries for level "<<l<<"\n";
    cout <<l<<"\n";
#endif
    if (l>0) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += BC_update_FINE_TO_COARSE_SEND(SimPM,
                spatial_solver, l, grid[l]->BC_bd[i], 2,2);
        }
      }
    }
    if (l<SimPM.grid_nlevels-1) {
      for (size_t i=0;i<grid[l]->BC_bd.size();i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += BC_update_FINE_TO_COARSE_RECV(SimPM,spatial_solver,
                      l,grid[l]->BC_bd[i],2,2);
        }
      }
    }
  }
  BC_FINE_TO_COARSE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------
  err+= output_data(grid);
  rep.errorTest("MPI_NG TIME_INT::output_data()",0,err);


  
  while (SimPM.maxtime==false) {
    // --------------------------------------------------------------
    // Get timestep on each level
    int scale = 1;
    double mindt = 1.0e99;
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TEST_INT
      cout <<"Calculate timestep, level "<<l<<", dx=";
      cout <<SimPM.levels[l].dx<<"\n";
#endif
      if (!first_step) SimPM.last_dt = SimPM.levels[l].dt;

      err += calculate_timestep(SimPM, grid[l],spatial_solver,l);
      rep.errorTest("TIME_INT::calc_timestep()",0,err);
      
      mindt = std::min(mindt, SimPM.dt/scale);
      mindt = COMM->global_operation_double("MIN",mindt);
#ifdef TEST_INT
      cout <<"level "<<l<<" got dt="<<SimPM.dt<<" and ";
      cout <<SimPM.dt/scale <<"\n";
#endif
      SimPM.levels[l].dt = SimPM.dt;
      scale *= 2;
    }
    // make sure all levels use same step (scaled by factors of 2).
    scale = 1;
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
      SimPM.levels[l].dt = mindt*scale;
      scale *= 2;
#ifdef TEST_INT
      cout <<"new dt="<<SimPM.levels[l].dt<<", t=";
      cout <<SimPM.levels[l].simtime<<"\n";
#endif
    }
    if (first_step) {
      // take a 10x smaller timestep for the first timestep.
      for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
        //cout <<"level "<<l<<", orig dt="<<SimPM.levels[l].dt;
        SimPM.levels[l].dt *=0.1;
      }
      first_step=false;
    }
    // --------------------------------------------------------------
    
    // --------------------------------------------------------------
#ifdef TESTING
    cout <<"NG_MPI time_int: stepping forward in time\n";
#endif
    // Use a recursive algorithm to update the coarsest level.  This
    // function also updates the next level twice, by calling itself
    // for the finer level, and so on.
    //
    advance_time(0,grid[0]);
    //advance_step_OA1(0);
    SimPM.simtime = SimPM.levels[0].simtime;
    COMM->barrier("step");
#ifdef TESTING
    cout <<"MPI time_int: finished timestep\n";
#endif

    //if ( (SimPM.levels[0].MCMD.get_myrank()==0) &&
    //     (SimPM.timestep%log_freq)==0) {
      cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime;
      cout <<"\t timestep: "<<SimPM.timestep;
      tsf=clk.time_so_far("time_int");
      cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
//#ifdef TESTING
      cout.flush();
//#endif // TESTING
    //}
    // --------------------------------------------------------------
	
    err += check_energy_cons(grid);

    //
    // check if we are at time limit yet.
    //
    tsf=clk.time_so_far("time_int");
    double maxt = COMM->global_operation_double("MAX", tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime=true;
      cout <<"RUNTIME>"<<get_max_walltime()<<" SECS.\n";
    }
	
    err+= output_data(grid);
    rep.errorTest("MPI_NG TIME_INT::output_data()",0,err);

    err+= check_eosim();
    rep.errorTest("MPI_NG TIME_INT::check_eosim()",0,err);
  }
  cout <<"sim_control_NG_MPI:: TIME_INT FINISHED.  MOVING ON TO ";
  cout <<"FINALISE SIM.\n";
  tsf=clk.time_so_far("time_int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep;
  cout <<", sim-time="<<SimPM.simtime;
  cout <<", wall-time=" <<tsf;
  cout <<", time/step="<<tsf/static_cast<double>(SimPM.timestep);
  cout <<"\n";
  if (grid[0]->RT!=0) {
    // print raytracing timing info.  Start and stop timers to get 
    // the correct runtime
    string t1="totalRT", t2="waitingRT", t3="doingRT";
    double total=0.0, wait=0.0, run=0.0;
    clk.start_timer(t1); total = clk.pause_timer(t1);
    clk.start_timer(t2); wait  = clk.pause_timer(t2);
    clk.start_timer(t3); run   = clk.pause_timer(t3);
    cout <<"TOTALS RT#: active="<<run<<" idle="<<wait;
    cout <<" total="<<total<<"\n";
  }
  cout <<"                *************************************\n\n";
  return(0);
}



// ##################################################################
// ##################################################################



double sim_control_NG_MPI::advance_step_OA1(
      const int l  ///< level to advance.
      )
{
#ifdef TEST_INT
  cout <<"NG-MPI advance_step_OA1, level="<<l<<", starting.\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  double ctime = SimPM.levels[l].simtime; // current time
  class GridBaseClass *grid = SimPM.levels[l].grid;
  class GridBaseClass *child = SimPM.levels[l].child;
  bool finest_level = (l<(SimPM.grid_nlevels-1)) ? false : true;

  err = update_evolving_RT_sources(
            SimPM,SimPM.levels[l].simtime,grid->RT);
  rep.errorTest("NG TIME_INT::update_RT_sources error",0,err);

#ifdef TEST_INT
  cout <<"advance_step_OA1: child="<<child<<"\n";
  cout <<"finest_level="<<finest_level<<", l="<<l<<", max="<<SimPM.grid_nlevels<<"\n";
#endif

  // 0. See if there are coarse-to-fine boundary data to send
  int c2f = -1, f2cs=-1, f2cr=-1;
  if (!finest_level) {
    // C2F data to send:
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_SEND) c2f=i;
    }
    // F2C data to receive
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_RECV) f2cr=i;
    }
  }
  if (l>0) {
    // F2C data to send
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_SEND) f2cs=i;
    }
    // C2F data to recv can be more than one external BC, so find
    // them later in a loop.
  }
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" c2f = "<<c2f;
  cout <<", f2c send="<<f2cs<<", f2c recv="<<f2cr<<"\n";
#endif

  // --------------------------------------------------------
  // 1. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" first fine step\n";
#endif
    dt2_fine = advance_step_OA1(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" dt="<<dt2_this<<"\n";
#endif

  // --------------------------------------------------------
  // 2. Calculate dU for this level (1st order)
  // --------------------------------------------------------
  spatial_solver->Setdt(dt2_this);
  if (grid->RT && (!FVI_need_column_densities_4dt ||
    (SimPM.levels[l].step%SimPM.levels[l].multiplier !=0) )) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" raytracing\n";
#endif
    err += calculate_raytracing_column_densities(SimPM,grid,l);
    rep.errorTest("NG-MPI::advance_step_OA1: calc_rt_cols()",0,err);
  }
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" microphysics\n";
#endif
  err += calc_microphysics_dU(dt2_this, grid);
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" dynamics\n";
#endif
  err += calc_dynamics_dU(dt2_this,TIMESTEP_FIRST_PART, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt2_this,TIMESTEP_FIRST_PART,grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("NG-MPI scn::advance_step_OA1: calc_x_dU",0,err);


  // --------------------------------------------------------
  // 3. Send external boundary data to finer grid (C2F)
  //    The (1,2) tells it that we are only half way
  //    through the coarse-level step.
  //    Then Receive the data on level l+1
  //    Finally clear the C2F sends.
  // --------------------------------------------------------
  if (c2f>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" C2F SEND\n";
#endif
    // send C2F data from l
    err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid,spatial_solver,
                                         l, grid->BC_bd[c2f], 1,2);
    // receive C2F data to l+1
    for (size_t i=0;i<child->BC_bd.size();i++) {
      if (child->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                    l+1,child->BC_bd[i],SimPM.levels[l+1].step);
      }
    }
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" C2F CLEAR SEND\n";
#endif
    // clear sends.
    BC_COARSE_TO_FINE_SEND_clear_sends();
  }

  // --------------------------------------------------------
  // 4. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" second fine step\n";
#endif
    dt2_fine = advance_step_OA1(l+1);
  }

  // --------------------------------------------------------
  // 5. Update grid and boundaries on level l:
  // --------------------------------------------------------
  spatial_solver->Setdt(dt2_this);
  //  - Receive level fluxes from finer grid (FLUX)
  //  - update grid on level l to new time
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" update grid\n";
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt,OA1,OA1, grid);
  rep.errorTest("scn::advance_step_OA1: state-vec update",0,err);

  //  - clear level flux sends
  clear_sends_BC89_fluxes();

  //  - update boundary conditions on level l
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" update BCs\n";
#endif
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                              ctime+dt2_this, OA1, OA1);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                              ctime+dt2_this, OA1, OA1);
  
  //  - Recv F2C data from l+1
  if (!finest_level && f2cr>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" receive F2C data\n";
    cout <<"F2C Recv, l="<<l<<", f2cr="<<f2cr<<", size=";
    cout <<grid->BC_bd.size()<<"\n";
#endif
    err += BC_update_FINE_TO_COARSE_RECV(
                SimPM,spatial_solver,l,grid->BC_bd[f2cr],OA1,OA1);
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" clear F2C sends\n";
#endif
    //  - Clear F2C sends
    BC_FINE_TO_COARSE_SEND_clear_sends();
  }


  // --------------------------------------------------------
  // 6. Send external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are all the way through
  //    the coarse-level step.
  // --------------------------------------------------------
  if (c2f>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" C2F SEND\n";
#endif
    err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid,spatial_solver,
                                         l, grid->BC_bd[c2f], 2,2);
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" update l+1 C2F BCs\n";
#endif
    for (size_t i=0;i<child->BC_bd.size();i++) {
      if (child->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                    l+1,child->BC_bd[i],SimPM.levels[l+1].step);
      }
    }
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" C2F CLEAR SEND\n";
#endif
    BC_COARSE_TO_FINE_SEND_clear_sends();
  }

  // --------------------------------------------------------
  // 7. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  if (l>0 && SimPM.levels[l].step%2!=0) {
    // - send level fluxes
    err += send_BC89_fluxes_F2C(l,OA1,OA1);
#ifdef TEST_INT
    cout <<"advance_step_OA1: l="<<l<<" send_BC89_fluxes_F2C()\n";
#endif

    // - send F2C data
#ifdef TEST_MPI_NG
    cout <<"LEVEL "<<l<<": update_bcs_NG_MPI: updating bc ";
    cout <<f2cs<<" with type "<<grid->BC_bd[f2cs]->type<<"\n";
    cout <<"found FINE_TO_COARSE_SEND boundary to update\n";
    cout <<"... step="<<SimPM.levels[l].step<<"\n";
#endif
    err += BC_update_FINE_TO_COARSE_SEND(
              SimPM,spatial_solver,l,grid->BC_bd[f2cs],OA1,OA1);
  }


  // --------------------------------------------------------
  // 8. increment time and timestep for this level
  // --------------------------------------------------------
#ifdef TEST_INT
  cout <<"advance_step_OA1: l="<<l<<" advancing time and return\n";
#endif
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }


#ifdef TESTING
  cout <<"NG-MPI advance_step_OA1, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step;
  cout <<", next dt="<<SimPM.levels[l].dt<<" next time=";
  cout << SimPM.levels[l].simtime + SimPM.levels[l].dt <<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}




// ##################################################################
// ##################################################################



double sim_control_NG_MPI::advance_step_OA2(
      const int l  ///< level to advance.
      )
{
#ifdef TEST_INT
  cout <<"NG-MPI advance_step_OA2, level="<<l<<", starting.\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  double ctime = SimPM.levels[l].simtime; // current time
  class GridBaseClass *grid = SimPM.levels[l].grid;
  class GridBaseClass *child = SimPM.levels[l].child;
  bool finest_level = (l<SimPM.grid_nlevels-1) ? false : true;

  err = update_evolving_RT_sources(
            SimPM,SimPM.levels[l].simtime,grid->RT);
  rep.errorTest("NG TIME_INT::update_RT_sources error",0,err);

#ifdef TEST_INT
  cout <<"advance_step_OA2: child="<<child<<"\n";
#endif

  // 0. See if there are coarse-to-fine boundary data to send
  int c2f = -1, f2cs=-1, f2cr=-1;
  if (!finest_level) {
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_SEND) c2f=i;
    }
    // F2C data to receive
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_RECV) f2cr=i;
    }
  }
  if (l>0) {
    // F2C data to send
    for (size_t i=0;i<grid->BC_bd.size();i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_SEND) f2cs=i;
    }
    // C2F data to recv can be more than one external BC, so find
    // them later in a loop.
  }
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" c2f = "<<c2f<<"\n";
#endif



  // --------------------------------------------------------
  // 1. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" first fine step\n";
#endif
    dt2_fine = advance_step_OA2(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" dt="<<dt2_this<<"\n";
#endif

  // --------------------------------------------------------
  // 2. Calculate dU for this level (1st order)
  // --------------------------------------------------------
  // Predictor step: use 0.5*dt to get to time-centred state
  double dt_now = dt2_this*0.5;             // half of the timestep
  spatial_solver->Setdt(dt_now);
  if (grid->RT && (!FVI_need_column_densities_4dt ||
    (SimPM.levels[l].step%SimPM.levels[l].multiplier !=0) )) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" raytracing\n";
#endif
    err += calculate_raytracing_column_densities(SimPM,grid,l);
    rep.errorTest("NG-MPI::advance_step_OA2: calc_rt_cols()",0,err);
  }
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" microphysics\n";
#endif
  err += calc_microphysics_dU(dt_now, grid);
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" dynamics\n";
#endif
  err += calc_dynamics_dU(dt_now,TIMESTEP_FIRST_PART, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now,TIMESTEP_FIRST_PART,grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("NG-MPI scn::advance_step_OA2: calc_x_dU",0,err);
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" 1/2 step, grid_update\n";
#endif

  // --------------------------------------------------------
  // 3. Update grid and boundary data on this level at 1/2 step
  // --------------------------------------------------------
  // update state vector Ph to half-step values
  err += grid_update_state_vector(dt_now,
                                    TIMESTEP_FIRST_PART,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA2",0,err);

  // - update internal and external boundaries locally
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                                    ctime+dt_now, OA1, OA2);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                                    ctime+dt_now, OA1, OA2);
  rep.errorTest("scn::advance_step_OA2: bounday update OA2",0,err);
  
  // - Receive F2C data from l+1
  //if (!finest_level && f2cr>=0) {
#ifdef TEST_INT
  //  cout <<"advance_step_OA2: l="<<l<<" receive F2C data\n";
  //  cout <<"F2C Recv, l="<<l<<", f2cr="<<f2cr<<", size=";
  //  cout <<grid->BC_bd.size()<<"\n";
#endif
  //  err += BC_update_FINE_TO_COARSE_RECV(
  //              SimPM,spatial_solver,l,grid->BC_bd[f2cr],OA1,OA2);
#ifdef TEST_INT
  //  cout <<"advance_step_OA2: l="<<l<<" clear F2C sends\n";
#endif
  //  //  - Clear F2C sends
  //  BC_FINE_TO_COARSE_SEND_clear_sends();
  //}

  // --------------------------------------------------------
  // 4. Calculate dU for the full step (OA2) on this level
  // --------------------------------------------------------
  //
  dt_now = dt2_this;  // full step
  spatial_solver->Setdt(dt_now);
  if (grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid,l);
    rep.errorTest("scn::advance_time: calc_rt_cols() OA2",0,err);
  }
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l;
  cout <<" full step, start calc_microphysics_dU\n";
#endif
  err += calc_microphysics_dU(dt_now, grid);
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l;
  cout <<" full step, start calc_dynamics_dU\n";
#endif
  err += calc_dynamics_dU(dt_now, TIMESTEP_FULL, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now,TIMESTEP_FULL, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA2",0,err);

  // --------------------------------------------------------
  // 5. Send/receive external boundary data to finer grid (C2F)
  //    The (1,2) tells it that we are half-way through
  //    the coarse-level step.
  // NB: This doesn't have corrected fluxes, so there is a
  //     possible inconsistency here.  Need to make sure that
  //     it doesn't cause conservation issues...
  // --------------------------------------------------------
  if (c2f>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" C2F SEND\n";
#endif
    err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid,spatial_solver,
                                         l, grid->BC_bd[c2f], 1,2);
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" update l+1 C2F BCs\n";
#endif
    for (size_t i=0;i<child->BC_bd.size();i++) {
      if (child->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                    l+1,child->BC_bd[i],SimPM.levels[l+1].step);
      }
    }
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" C2F CLEAR SEND\n";
#endif
    BC_COARSE_TO_FINE_SEND_clear_sends();
  }

  // --------------------------------------------------------
  // 6. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" second fine step\n";
#endif
    dt2_fine = advance_step_OA2(l+1);
  }

  // --------------------------------------------------------
  // 7. Update local grid to new time values
  // --------------------------------------------------------
  spatial_solver->Setdt(dt2_this);
  //  - Receive level fluxes from finer grid (FLUX)
  //  - update grid on level l to new time
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" update grid\n";
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt,OA2,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update",0,err);
  //  - clear level flux sends
  clear_sends_BC89_fluxes();

  //  - update boundary conditions on level l
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" update BCs\n";
#endif
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                              ctime+dt2_this, OA2, OA2);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                              ctime+dt2_this, OA2, OA2);
  
  //  - Recv F2C data from l+1
  if (!finest_level && f2cr>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" receive F2C data\n";
    cout <<"F2C Recv, l="<<l<<", f2cr="<<f2cr<<", size=";
    cout <<grid->BC_bd.size()<<"\n";
#endif
    err += BC_update_FINE_TO_COARSE_RECV(
                SimPM,spatial_solver,l,grid->BC_bd[f2cr],OA2,OA2);
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" clear F2C sends\n";
#endif
    //  - Clear F2C sends
    BC_FINE_TO_COARSE_SEND_clear_sends();
  }

  // --------------------------------------------------------
  // 8. Send/receive external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are all the way through
  //    the coarse-level step.
  // --------------------------------------------------------
  if (c2f>=0) {
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" C2F SEND\n";
#endif
    err += BC_update_COARSE_TO_FINE_SEND(SimPM,grid,spatial_solver,
                                         l, grid->BC_bd[c2f], 2,2);
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" update l+1 C2F BCs\n";
#endif
    for (size_t i=0;i<child->BC_bd.size();i++) {
      if (child->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(SimPM,spatial_solver,
                    l+1,child->BC_bd[i],SimPM.levels[l+1].step);
      }
    }
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" C2F CLEAR SEND\n";
#endif
    BC_COARSE_TO_FINE_SEND_clear_sends();
  }

  // --------------------------------------------------------
  // 9. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  // - send F2C data
  if (l>0 && SimPM.levels[l].step%2!=0) {
    // - send level fluxes
#ifdef TEST_INT
    cout <<"step="<<SimPM.levels[l].step<<"\n";
#endif
    err += send_BC89_fluxes_F2C(l,OA2,OA2);
#ifdef TEST_INT
    cout <<"advance_step_OA2: l="<<l<<" send_BC89_fluxes_F2C()\n";
#endif

#ifdef TEST_MPI_NG
    cout <<"LEVEL "<<l<<": update_bcs_NG_MPI: updating bc ";
    cout <<f2cs<<" with type "<<grid->BC_bd[f2cs]->type<<"\n";
    cout <<"found FINE_TO_COARSE_SEND boundary to update\n";
#endif
    err += BC_update_FINE_TO_COARSE_SEND(
              SimPM,spatial_solver,l,grid->BC_bd[f2cs],OA2,OA2);
  }


  // --------------------------------------------------------
  // 14. increment time and timestep for this level
  // --------------------------------------------------------
#ifdef TEST_INT
  cout <<"advance_step_OA2: l="<<l<<" advancing time and return\n";
#endif
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }


#ifdef TESTING
  cout <<"NG-MPI advance_step_OA2, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step;
  cout <<", next dt="<<SimPM.levels[l].dt<<" next time=";
  cout << SimPM.levels[l].simtime + SimPM.levels[l].dt <<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}



// ##################################################################
// ##################################################################



void sim_control_NG_MPI::clear_sends_BC89_fluxes()
{
  for (unsigned int ib=0; ib<BC89_flux_send_list.size(); ib++) {
    COMM->wait_for_send_to_finish(BC89_flux_send_list[ib]);
  }
  BC89_flux_send_list.clear();
  return;
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::send_BC89_fluxes_F2C(
      const int l,      ///< My level in grid hierarchy.
      const int step,   ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa     ///< Full order of accuracy of simulation
      )
{
  if (step != ooa) {
    cout <<"l="<<l<<": don't send fluxes on half step\n";
    return 1;
  }
  if (SimPM.levels[l].step%2 ==0 ) {
    rep.error("Don't call BC89 Flux-SEND on odd steps",1);
  }
  if (l==0) {
    rep.error("Coarsest level trying to send flux data",l);
  }
  int err=0;

  // else we have to send data to at least one other MPI process:
  class GridBaseClass *grid = SimPM.levels[l].grid;
  class MCMDcontrol *MCMD = &(SimPM.levels[l].MCMD);
  int n_send = grid->flux_update_send.size();
  for (int isend=0;isend<n_send;isend++) {
    // loop over boundaries (some send to more than one process)
    struct flux_update *fup = &(grid->flux_update_send[isend]);
    // some sends may be null, so we skip them:
    if (fup->fi[0]==0) {
#ifdef TEST_BC89FLUX
      cout <<"send "<<isend<<" is null, continuing...\n";
#endif
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": send "<<isend<<" is not null: sending data.\n";
#endif
    }
    size_t n_el = fup->fi.size();
    size_t n_data = n_el * SimPM.nvar;
    pion_flt *data = 0; // buffer for sending by MPI
    data = mem.myalloc(data,n_data);
    size_t iel=0;
    for (size_t ii=0;ii<n_el;ii++) {
      for (int v=0;v<SimPM.nvar;v++) {
        data[iel] = fup->fi[ii]->flux[v];
        iel++;
      }
    }
    //
    // Send data using a non-blocking MPI send
    //
    string id;
    int comm_tag = BC_MPI_FLUX_tag +100*isend +l;

    // loop over procs on level l-1 to send to.
    for (unsigned int ii=0; ii<fup->rank.size();ii++) {
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": isend="<<isend<<", send "<<ii<<" of ";
      cout <<fup->rank.size()<<" to rank "<<fup->rank[ii]<<".\n";
#endif
      if (fup->rank[ii] == MCMD->get_myrank()) {
#ifdef TEST_BC89FLUX
        cout <<"l="<<l<<": ignoring BC89 for isend="<<isend<<" (to myrank).\n";
#endif
        continue;
      }
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": BC89 FLUX: Sending "<<n_data;
      cout <<" doubles from proc "<<MCMD->get_myrank();
      cout <<" to parent proc "<<fup->rank[ii]<<"\n";
#endif
      err += COMM->send_double_data(
            fup->rank[ii],n_data,data,id,comm_tag);
      if (err) rep.error("FLUX_F2C send_data failed.",err);
#ifdef TEST_MPI_NG_F2C
      cout <<"l="<<l<<": BC89 FLUX: returned with id="<<id;
      cout <<"\n";
#endif
      // store ID to clear the send later
      BC89_flux_send_list.push_back(id);
    } // loop over ranks
    data = mem.myfree(data);
  } // loop over send boundaries.
  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::recv_BC89_fluxes_F2C(
      const int l,      ///< My level in grid hierarchy.
      const int step,   ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa     ///< Full order of accuracy of simulation
      )
{
  if (step != ooa) {
    cout <<"don't receive fluxes on half step\n";
    return 1;
  }
  if (l==SimPM.grid_nlevels-1) {
    rep.error("finest level trying to receive data from l+1",l);
  }
  int err=0;

  // loop over all boundaries that we need to receive
  double ftmp[SimPM.nvar], utmp[SimPM.nvar];
  class GridBaseClass *grid = SimPM.levels[l].grid;
  class MCMDcontrol *MCMD = &(SimPM.levels[l].MCMD);
  int n_bd = grid->flux_update_recv.size();
  
  for (int irecv=0;irecv<n_bd;irecv++) {
    struct flux_update *fup = &(grid->flux_update_recv[irecv]);
    // some recvs may be null, so we skip them:
    if (fup->fi[0]==0) {
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": recv "<<irecv<<" is null, continuing...\n";
#endif
      continue;
    }
    else if (fup->rank[0] == MCMD->get_myrank()) {
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": calling serial BC89 for irecv="<<irecv;
      cout <<" (same rank). dir="<<fup->dir<<endl;
#endif
      int d = fup->dir;
      struct flux_update *fsend = 
                  &(SimPM.levels[l].child->flux_update_send[d]);
      err = recv_BC89_flux_boundary(grid,*fsend,*fup,d,
                      static_cast<axes>(fup->ax));
      rep.errorTest("serial BC89 call",0,err);
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      cout <<"l="<<l<<": recv "<<irecv<<" is not null: recving data."<<endl;
#endif
    }

    // Normally we have nchild*2*ndim boundaries, so the direction
    // associated with a given "irecv" is irecv%(2*ndim).
    // But if only the grid boundary is in common with level l+1, 
    // then there are only 2^(ndim-1) elements, one for each l+1
    // grid that is adjacent to the level l grid.
    int dir = -1;
    if (n_bd >= 2*SimPM.ndim) dir = irecv % 2*SimPM.ndim;
    else                      dir = irecv;

    struct flux_interface *fi=0;
    size_t n_el = fup->fi.size();
    size_t n_data = n_el * SimPM.nvar;
    pion_flt *data = 0; // buffer for sending by MPI
    data = mem.myalloc(data,n_data);

    // receive data: comm_tag ensures that we select the data
    // relating to this value of "irecv".
    string recv_id;
    int recv_tag=-1;
    int from_rank=-1;
    int comm_tag = BC_MPI_FLUX_tag +100*dir +l+1;
    err = COMM->look_for_data_to_receive(
          &from_rank, // rank of sender (output)
          recv_id,    // identifier for receive (output).
          &recv_tag,  // comm_tag associated with data (output)
          comm_tag,
          COMM_DOUBLEDATA // type of data we want.
          );
    if (err) rep.error("FLUX look for double data failed",err);
#ifdef TEST_BC89FLUX
    cout <<"l="<<l<<": BC89 FLUX: found data from rank ";
    cout <<from_rank<<", and ID "<<recv_id<<endl;
#endif

    pion_flt *buf = 0;
    buf = mem.myalloc(buf,n_data);
    //
    // Receive data into buffer.  Data stored for each coarse cell:
    // Ph[nv],cellvol,cellpos[nd],slopeX[nv],slopeY[nv],slopeZ[nv]
    //
    err = COMM->receive_double_data(
                          from_rank, recv_tag, recv_id, n_data, buf);
    if (err) rep.error("(BC_update_C2F_RECV) getdata failed",err);

    size_t iel=0;
    for (size_t ii=0;ii<n_el;ii++) {
      fi = fup->fi[ii];
      // construct error in flux:
      for (int v=0;v<SimPM.nvar;v++) {
        // make flux[v] be the difference of coarse and fine.
        fi->flux[v] += buf[iel];
        // divide by face area so that it is a flux.
        fi->flux[v] /=  fi->area[0];
#ifdef TEST_BC89FLUX
        if (!isfinite(buf[iel])) {
          cout <<"l="<<l<<": element "<<ii<<" of FLUX C2F RECV: ";
          cout <<" var "<<iel<<" is "<<buf[iel]<<endl;
          rep.error("infinite",buf[ii]);
        }
#endif
        iel++;
      }
      // re-calculate dU based on error in flux.
      if (fup->dir%2 ==0) {
        spatial_solver->DivStateVectorComponent(fi->c[0], grid,
          static_cast<axes>(fup->ax), SimPM.nvar,ftmp,fi->flux,utmp);
      }
      else {
        spatial_solver->DivStateVectorComponent(fi->c[0], grid,
          static_cast<axes>(fup->ax), SimPM.nvar,fi->flux,ftmp,utmp);
      }
#ifdef TEST_BC89FLUX
      rep.printVec("**********  Error",utmp, SimPM.nvar);
#endif
      // correct dU so that coarse level is consistent with fine.
      for (int v=0;v<SimPM.nvar;v++) {
        fi->c[0]->dU[v] += utmp[v];
      }

    } // loop over elements
    
 } // loop over faces in this direction.

      
  return 0;
}



// ##################################################################
// ##################################################################




int sim_control_NG_MPI::grid_update_state_vector(
      const double dt,  ///< timestep
      const int step, ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa,   ///< Full order of accuracy of simulation
      class GridBaseClass *grid ///< Computational grid.
      )
{
  // get current level of grid in hierarchy.
  int level=0;
  int err=0;
  for (int v=0;v<SimPM.grid_nlevels;v++) {
    if (grid == SimPM.levels[v].grid) level = v;
  }

  //
  // loop through the faces, and compare the fine with coarse level
  // fluxes.  Subtract the flux difference from the cells adjacent to
  // the face.
  //
  if (step==ooa && level < SimPM.grid_nlevels-1) {
#ifdef DEBUG_NG
    cout <<"level "<<level<<": correcting fluxes from finer grid\n";
#endif
    err = recv_BC89_fluxes_F2C(level,step,ooa);
    rep.errorTest("return from recv_fluxes_F2C",0,err);
  }

  err = time_integrator::grid_update_state_vector(dt,step,
                                                   ooa,grid);
  rep.errorTest("return from grid_update_state_vec",0,err);
  return err;
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::initial_conserved_quantities(
      vector<class GridBaseClass *> &grid
      )
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION 
  pion_flt u[SimPM.nvar];
  //dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange = 0.0;
  //  cout <<"initERG: "<<dp.initERG<<"\n";
  initERG = 0.;  initMMX = initMMY = initMMZ = 0.; initMASS = 0.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    double dx = SimPM.levels[l].dx;
    double dv = 0.0;
    class cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd && c->isgd) {
        dv = spatial_solver->CellVolume(c,dx);
        spatial_solver->PtoU(c->P,u,SimPM.gamma);
        initERG += u[ERG]*dv;
        initMMX += u[MMX]*dv;
        initMMY += u[MMY]*dv;
        initMMZ += u[MMZ]*dv;
        initMASS += u[RHO]*dv;
      }
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }

  initERG = COMM->global_operation_double("SUM", initERG);
  initMMX = COMM->global_operation_double("SUM", initMMX);
  initMMY = COMM->global_operation_double("SUM", initMMY);
  initMMZ = COMM->global_operation_double("SUM", initMMZ);
  initMASS = COMM->global_operation_double("SUM", initMASS);

  cout <<"(conserved quantities) ["<< initERG <<", ";
  cout << initMMX <<", ";
  cout << initMMY <<", ";
  cout << initMMZ <<", ";
  cout << initMASS <<"]\n";

#endif // TEST_CONSERVATION
  return(0);
}






// ##################################################################
// ##################################################################



int sim_control_NG_MPI::check_energy_cons(
      vector<class GridBaseClass *> &grid
      )
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION 
  pion_flt u[SimPM.nvar];
  nowERG=0.;
  nowMMX = 0.;
  nowMMY = 0.;
  nowMMZ = 0.;
  nowMASS = 0.0;
  double totmom=0.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    double dx = SimPM.levels[l].dx;
    double dv = 0.0;
    class cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd && c->isgd) {
        dv = spatial_solver->CellVolume(c,dx);
        spatial_solver->PtoU(c->P,u,SimPM.gamma);
        nowERG += u[ERG]*dv;
        nowMMX += u[MMX]*dv;
        nowMMY += u[MMY]*dv;
        nowMMZ += u[MMZ]*dv;
        nowMASS += u[RHO]*dv;
        totmom += sqrt(u[MMX]*u[MMX] + u[MMY]*u[MMY] + u[MMZ]*u[MMZ])
                   *dv;
      }
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }

  //cout <<totmom;
  nowERG = COMM->global_operation_double("SUM", nowERG);
  nowMMX = COMM->global_operation_double("SUM", nowMMX);
  nowMMY = COMM->global_operation_double("SUM", nowMMY);
  nowMMZ = COMM->global_operation_double("SUM", nowMMZ);
  nowMASS = COMM->global_operation_double("SUM", nowMASS);
  totmom = COMM->global_operation_double("SUM", totmom);
  cout <<" totmom="<<totmom<<" initMMX="<<initMMX;
  cout <<", nowMMX="<<nowMMX<<"\n";

  cout <<"(conserved quantities) ["<< nowERG <<", ";
  cout << nowMMX <<", ";
  cout << nowMMY <<", ";
  cout << nowMMZ <<", ";
  cout << nowMASS <<"]\n";
  cout <<"(relative error      ) [";
  cout << (nowERG-initERG)/(initERG) <<", ";
  cout << (nowMMX-initMMX)/(totmom) <<", ";
  cout << (nowMMY-initMMY)/(totmom) <<", ";
  cout << (nowMMZ-initMMZ)/(totmom) <<", ";
  cout << (nowMASS-initMASS)/initMASS <<"]\n";

#endif // TEST_CONSERVATION
  return(0);
}



// ##################################################################
// ##################################################################






#endif // PARALLEL


