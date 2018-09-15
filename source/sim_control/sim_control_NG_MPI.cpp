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
#include "constants.h"

#include "decomposition/MCMD_control.h"
#include "sim_control/sim_control_NG_MPI.h"
#include "raytracing/raytracer_SC.h"

//#ifdef SILO
//#include "dataIO/dataio_silo.h"
//#include "dataIO/dataio_silo_utility.h"
//#endif // if SILO
//#ifdef FITS
//#include "dataIO/dataio_fits.h"
//#include "dataIO/dataio_fits_MPI.h"
//#endif // if FITS

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#ifdef PARALLEL




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
  }
  // ----------------------------------------------------------------



  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"NG_MPI updating external boundaries for level "<<l<<"\n";
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
#ifdef TESTING
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
    c = (grid[l])->FirstPt_All();
    do {
      if (pconst.equalD(c->P[RO],0.0)) {
        cout <<"zero data in cell: ";
        CI.print_cell(c);
      }
    } while ( (c=(grid[l])->NextPt_All(c)) !=0 );
  }
#endif // TESTING
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
  clk.start_timer("time_int"); double tsf=0.0;

  // make sure all levels start at the same time.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }
  
  while (SimPM.maxtime==false) {
    
    // --------------------------------------------------------------
    // Update RT sources and boundaries.
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
      mindt = COMM->global_op_double("MIN",mindt);
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
    advance_time(0);
    SimPM.simtime = SimPM.levels[0].simtime;
#ifdef TESTING
    cout <<"MPI time_int: finished timestep\n";
#endif

    if ( (SimPM.levels[0].MCMD.get_myrank()==0) &&
         (SimPM.timestep%log_freq)==0) {
      cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime;
      cout <<"\t timestep: "<<SimPM.timestep;
      tsf=clk.time_so_far("time_int");
      cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#ifdef TESTING
      cout.flush();
#endif // TESTING
    }
    // --------------------------------------------------------------
	
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


int sim_control_NG_MPI::calculate_timestep(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      class FV_solver_base *sp_solver, ///< solver/equations class
      const int l       ///< level to advance (for NG grid)
      )
{
  //
  // First get the local grid dynamics and microphysics timesteps.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt(par,grid,sp_solver);
  t_mp  = calc_microphysics_dt(par,grid,l);
  
  //
  // Now get global min over all grids for dynamics and microphysics timesteps.
  // We only need both if we are doing Dedner et al. 2002, mixed-GLM divergence
  // cleaning of the magnetic field, since there the dynamical dt is an important
  // quantity (see next block of code).
  //
  //par.dt = t_dyn;
  t_dyn = COMM->global_operation_double("MIN", t_dyn);
  //cout <<"proc "<<SimPM.levels[0].MCMD.get_myrank();
  //cout<<":\t my t_dyn="<<par.dt<<" and global t_dyn="<<t_dyn<<"\n";
  //par.dt = t_mp;
  t_mp = COMM->global_operation_double("MIN", t_mp);
  //cout <<"proc "<<SimPM.levels[0].MCMDM.get_myrank();
  //cout<<":\t my t_mp ="<<par.dt<<" and global t_mp ="<<t_mp<<"\n";
  
  // Write step-limiting info every tenth timestep.
  if (t_mp<t_dyn && (par.timestep%10)==0)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";

#ifdef THERMAL_CONDUCTION
  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to actually calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in eqn->preprocess_data()
  //
  double t_cond = calc_conduction_dt_and_Edot();
  t_cond = COMM->global_operation_double("MIN", t_cond);
  if (t_cond<par.dt) {
    cout <<"PARALLEL CONDUCTION IS LIMITING TIMESTEP: t_c="<<t_cond<<", t_m="<<t_mp;
    cout <<", t_dyn="<<t_dyn<<"\n";
  }
  par.dt = min(par.dt, t_cond);
#endif // THERMAL CONDUCTION

  //
  // if using MHD with GLM divB cleaning, the following sets the hyperbolic wavespeed.
  // If not, it does nothing.  By setting it here and using t_dyn, we ensure that the
  // hyperbolic wavespeed is equal to the maximum signal speed on the grid, and not
  // an artificially larger speed associated with a shortened timestep.
  //
  sp_solver->GotTimestep(t_dyn,grid->DX());

  //
  // Now the timestep is the min of the global microphysics and Dynamics timesteps.
  //
  SimPM.dt = min(t_dyn,t_mp);

  //
  // Check that the timestep doesn't increase too much between step, and that it 
  // won't bring us past the next output time or the end of the simulation.
  // This function operates on SimPM.dt, resetting it to a smaller value if needed.
  //
  timestep_checking_and_limiting(par);
  
  //
  // Tell the solver class what the resulting timestep is.
  //
  sp_solver->Setdt(par.dt);
  
#ifdef TESTING
  //
  // Check (sanity) that if my processor has modified dt to get to either
  // an output time or finishtime, then all processors have done this too!
  // This is really paranoid testing, and involves communication, so only 
  // do it when testing.
  //
  t_dyn=par.dt;
  t_mp = COMM->global_operation_double("MIN", t_dyn);
  if (!pconst.equalD(t_dyn,t_mp))
    rep.error("synchonisation trouble in timesteps!",t_dyn-t_mp);
#endif // TESTING

  return 0;
}


// ##################################################################
// ##################################################################



#endif // PARALLEL


