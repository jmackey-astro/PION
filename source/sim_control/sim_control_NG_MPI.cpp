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
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include "tools/timer.h"

#include "constants.h"

#include "decomposition/MCMD_control.h"
#include "raytracing/raytracer_SC.h"
#include "sim_control/sim_control_NG_MPI.h"

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef PARALLEL

//#define TEST_INT

// ##################################################################
// ##################################################################

sim_control_NG_MPI::sim_control_NG_MPI()
{
#ifdef TESTING
  cout << "sim_control_NG_MPI constructor.\n";
#endif
}

// ##################################################################
// ##################################################################

sim_control_NG_MPI::~sim_control_NG_MPI()
{
#ifdef TESTING
  cout << "sim_control_NG_MPI destructor.\n";
#endif
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::Init(
    string infile,
    int typeOfFile,
    int narg,
    string* args,
    vector<class GridBaseClass*>& grid  ///< address of vector of grid pointers.
)
{
  cout << "(pion) Init: infile = " << infile << "\n";
  int err = 0;

  //
  // Setup the MCMDcontrol class with rank and nproc.
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);

  // ----------------------------------------------------------------
  SimPM.typeofip = typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  if (!dataio->file_exists(infile)) rep.error("infile doesn't exist!", infile);

  // ----------------------------------------------------------------
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("NG_MPI Init(): failed to read header", 0, err);

  // Check if any commandline args override the file parameters.
  // ----------------------------------------------------------------
  err = override_params(narg, args);
  rep.errorTest("(NG_MPI INIT::override_params)", 0, err);

  // setup the nested grid levels, and decompose the domain on each
  // level
  // ----------------------------------------------------------------
  setup_NG_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
#ifdef TESTING
  cout << "NG_MPI Init: grid setup\n";
#endif
  err = setup_grid(grid, SimPM);
#ifdef TESTING
  cout << "NG_MPI Init: grid setup finished\n";
#endif
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(NG_MPI INIT::setup_grid) error", 0, err);

  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  // ----------------------------------------------------------------
  err = set_equations(SimPM);
  rep.errorTest("(NG_MPI INIT::set_equations)", 0, err);
  spatial_solver->SetEOS(SimPM.gamma);

  // set up Microphysics, if needed.
  // ----------------------------------------------------------------
  err = setup_microphysics(SimPM);
  rep.errorTest("(NG_MPI INIT::setup_microphysics)", 0, err);

  // assign data to the grid from snapshot file.
  // ----------------------------------------------------------------
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(NG_MPI INIT::assign_initial_data)", 0, err);

  // ----------------------------------------------------------------
  // Set Ph[] = P[], and then implement the boundary conditions.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    cell* c = grid[l]->FirstPt();
    double u[SimPM.nvar];
    // rep.printVec("First Point",c->P,SimPM.nvar);
    do {
      // make sure temperature of the gas is reasonable
      spatial_solver->PtoU(c->P, u, SimPM.gamma);
      spatial_solver->UtoP(u, c->P, SimPM.EP.MinTemperature, SimPM.gamma);
      // set Ph[]=P[]
      for (int v = 0; v < SimPM.nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = grid[l]->NextPt(c)) != 0);

    if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
#ifdef TESTING
      cout << "Initial state, zero-ing glm variable.\n";
#endif
      c = grid[l]->FirstPt();
      do {
        c->P[SI] = c->Ph[SI] = 0.;
      } while ((c = grid[l]->NextPt(c)) != 0);
    }
  }  // loop over levels

  // ----------------------------------------------------------------
  err = boundary_conditions(SimPM, grid);
  rep.errorTest("(NG_MPI INIT::boundary_conditions) err!=0", 0, err);

  // ----------------------------------------------------------------
  err += setup_raytracing(SimPM, grid);
  rep.errorTest("Failed to setup raytracer", 0, err);

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err = assign_boundary_data(SimPM, l, grid[l]);
    rep.errorTest("NG_MPI INIT::assign_boundary_data", 0, err);
    COMM->barrier("level assign boundary data");
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout << "NG_MPI updating external boundaries for level " << l << "\n";
    cout << "UPDATING EXTERNAL BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout << "NG_MPI updating C2F boundaries for level " << l << "\n";
    cout << "UPDATING C2F BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += BC_update_COARSE_TO_FINE_SEND(
              SimPM, grid[l], spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
#ifdef TESTING
        cout << "Init: l=" << l << ", C2F recv i=" << i
             << ", type=" << grid[l]->BC_bd[i]->type << endl;
#endif
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i],
              SimPM.levels[l].step);
        }
      }
    }
  }
  BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout << "NG_MPI updating external boundaries for level " << l << "\n";
    cout << "@@@@@@@@@@@@  UPDATING EXTERNAL BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TESTING
    cout << "NG_MPI updating internal boundaries for level " << l << "\n";
    cout << "@@@@@@@@@@@@  UPDATING INTERNAL BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    err += TimeUpdateInternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  rep.errorTest("NG_MPI INIT: error from bounday update", 0, err);

  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TESTING
    cout << "NG_MPI updating F2C boundaries for level " << l << "\n";
    cout << l << "\n";
#endif
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += BC_update_FINE_TO_COARSE_SEND(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += BC_update_FINE_TO_COARSE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
  BC_FINE_TO_COARSE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from bounday update", 0, err);
  // ----------------------------------------------------------------

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
  if (SimPM.op_criterion == 1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!", SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime + SimPM.opfreq_time;
    double tmp        = ((SimPM.simtime / SimPM.opfreq_time)
                  - floor(SimPM.simtime / SimPM.opfreq_time))
                 * SimPM.opfreq_time;
    SimPM.next_optime -= tmp;
  }

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  // ----------------------------------------------------------------
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {
      delete dataio;
      dataio = 0;
    }
    if (textio) {
      delete textio;
      textio = 0;
    }
    setup_dataio_class(SimPM, SimPM.typeofop);
    if (!dataio) rep.error("NG_MPI INIT:: dataio", SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

  if (SimPM.timestep == 0) {
#ifdef TESTING
    cout << "(NG_MPI INIT) Writing initial data.\n";
#endif
    output_data(grid);
  }

  // ----------------------------------------------------------------
  //#ifdef TESTING
  cell* c = 0;
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    // cout <<"LEVEL-ZERO-CHECK L="<<l<<"\n";
    c = (grid[l])->FirstPt_All();
    do {
      if (pconst.equalD(c->P[RO], 0.0)) {
        cout << "ZERO DATA IN CELL: ";
        CI.print_cell(c);
      }
    } while ((c = (grid[l])->NextPt_All(c)) != 0);
  }
  //#endif // TESTING
  cout << "-------------------------------------------------------\n";
  return (err);
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::Time_Int(
    vector<class GridBaseClass*>& grid  ///< grid pointers.
)
{
  cout << "-------------------------------------------------------\n";
  cout << "(pion-ng-mpi)  STARTING TIME INTEGRATION. \n";
  int err       = 0;
  SimPM.maxtime = false;

  bool first_step = true;
  bool restart    = true;
  if (SimPM.timestep != 0) first_step = false;
  if (SimPM.timestep == 0) restart = false;

  clk.start_timer("time_int");
  double tsf = 0.0;

  // make sure all levels start at the same time.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt      = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }

  // --------------------------------------------------------------
  // Update internal and external boundaries.

  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TEST_INT
    cout << "updating internal boundaries for level " << l << "... ";
#endif
    err += TimeUpdateInternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.levels[l].simtime, SimPM.tmOOA,
        SimPM.tmOOA);
#ifdef TEST_INT
    cout << "... done \n";
#endif
  }
  rep.errorTest("sim_control_NG_MPI: internal boundary", 0, err);
  // --------------------------------------------------------------
  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    if (l < SimPM.grid_nlevels - 1) {
#ifdef TEST_INT
      cout << "NG_MPI Receiving F2C boundaries for level " << l << "\n";
#endif
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += BC_update_FINE_TO_COARSE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }

#ifdef TEST_INT
    cout << "NG_MPI raytracing level " << l << "\n";
#endif
    err =
        update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid[l]->RT);
    rep.errorTest("NG TIME_INT::update_RT_sources error", 0, err);
    do_ongrid_raytracing(SimPM, grid[l], l);

    if (l > 0) {
#ifdef TEST_INT
      cout << "NG_MPI Sending F2C boundaries for level " << l << "\n";
#endif
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += BC_update_FINE_TO_COARSE_SEND(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
#ifdef TEST_INT
  cout << "NG_MPI updated F2C boundaries for all levels.\n";
#endif
  BC_FINE_TO_COARSE_SEND_clear_sends();
#ifdef TEST_INT
  cout << "NG_MPI F2C cleared all sends.\n";
#endif
  rep.errorTest("NG_MPI time-int: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // update coarse-to-fine level boundaries
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TEST_INT
    cout << "NG_MPI updating C2F boundaries for level " << l << "\n";
    cout << l << "\n";
#endif
#ifdef TEST_INT
    cout << "updating external boundaries for level " << l << "\n";
#endif
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i],
              SimPM.levels[l].step);
        }
      }
    }

    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.levels[l].simtime, SimPM.tmOOA,
        SimPM.tmOOA);
    rep.errorTest("sim_control_NG_MPI: external boundary", 0, err);

    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += BC_update_COARSE_TO_FINE_SEND(
              SimPM, grid[l], spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
  BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("NG_MPI time-int: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  while (SimPM.maxtime == false) {
    // --------------------------------------------------------------
    // Get timestep on each level
    int scale    = 1;
    double mindt = 1.0e99;

    for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TEST_INT
      cout << "Calculate timestep, level " << l << ", dx=";
      cout << SimPM.levels[l].dx << endl;
#endif

      if (!restart && !first_step) {
        SimPM.last_dt = SimPM.levels[l].last_dt;
      }
      else {
        SimPM.levels[l].last_dt = SimPM.last_dt / SimPM.levels[l].multiplier;
      }
      err += calculate_timestep(SimPM, grid[l], spatial_solver, l);
      rep.errorTest("TIME_INT::calc_timestep()", 0, err);

      mindt = std::min(mindt, SimPM.dt / scale);
      mindt = COMM->global_operation_double("MIN", mindt);
#ifdef TEST_INT
      cout << "level " << l << " got dt=" << SimPM.dt << " and ";
      cout << SimPM.dt / scale << "... mindt=" << mindt << "\n";
#endif
      SimPM.levels[l].dt = SimPM.dt;
      scale *= 2;
    }
    // make sure all levels use same step (scaled by factors of 2).
    scale = 1;
    for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
      SimPM.levels[l].dt = mindt * scale;
      scale *= 2;
#ifdef TEST_INT
      cout << "new dt=" << SimPM.levels[l].dt << ", t=";
      cout << SimPM.levels[l].simtime << "\n";
#endif
    }
    if (first_step) {
      // take a ~3x smaller timestep for the first timestep in case
      // a gentle start is needed.
      for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
        // cout <<"level "<<l<<", orig dt="<<SimPM.levels[l].dt;
        SimPM.levels[l].dt *= 0.3;
      }
      first_step = false;
    }
    if (restart) restart = false;
    SimPM.last_dt = SimPM.levels[0].last_dt;
    rep.errorTest("TIME_INT::calc_timestep()", 0, err);
    // --------------------------------------------------------------

    // --------------------------------------------------------------
#ifdef TESTING
    cout << "NG_MPI time_int: stepping forward in time\n";
#endif
    // Use a recursive algorithm to update the coarsest level.  This
    // function also updates the next level twice, by calling itself
    // for the finer level, and so on.
    //
    advance_time(0, grid[0]);
    SimPM.simtime = SimPM.levels[0].simtime;
    COMM->barrier("step");
#ifdef TESTING
    cout << "MPI time_int: finished timestep\n";
#endif

    if (SimPM.levels[0].MCMD.get_myrank() == 0) {
      cout << "New time: " << SimPM.simtime;
      cout << "\tdt: " << SimPM.levels[SimPM.grid_nlevels - 1].dt;
      cout << "\t steps: " << SimPM.timestep;
      cout << "\tl0 steps: "
           << SimPM.timestep / static_cast<int>(pow(2, SimPM.grid_nlevels - 1));
      tsf = clk.time_so_far("time_int");
      cout << "\t runtime: " << tsf << " s"
           << "\n";
#ifdef TESTING
      cout.flush();
#endif  // TESTING
    }
    // --------------------------------------------------------------

    err += output_data(grid);
    rep.errorTest("MPI_NG TIME_INT::output_data()", 0, err);

    err += check_energy_cons(grid);

    //
    // check if we are at time limit yet.
    //
    tsf         = clk.time_so_far("time_int");
    double maxt = COMM->global_operation_double("MAX", tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime = true;
      cout << "RUNTIME>" << get_max_walltime() << " SECS.\n";
    }

    err += check_eosim();
    rep.errorTest("MPI_NG TIME_INT::check_eosim()", 0, err);
  }
  cout << "sim_control_NG_MPI:: TIME_INT FINISHED.  MOVING ON TO ";
  cout << "FINALISE SIM.\n";
  tsf = clk.time_so_far("time_int");
  cout << "TOTALS ###: Nsteps=" << SimPM.timestep;
  cout << ", sim-time=" << SimPM.simtime;
  cout << ", wall-time=" << tsf;
  cout << ", time/step=" << tsf / static_cast<double>(SimPM.timestep);
  cout << "\n";
  if (grid[0]->RT != 0) {
    // print raytracing timing info.  Start and stop timers to get
    // the correct runtime
    string t1 = "totalRT", t2 = "waitingRT", t3 = "doingRT";
    double total = 0.0, wait = 0.0, run = 0.0;
    clk.start_timer(t1);
    total = clk.pause_timer(t1);
    clk.start_timer(t2);
    wait = clk.pause_timer(t2);
    clk.start_timer(t3);
    run = clk.pause_timer(t3);
    cout << "TOTALS RT#: active=" << run << " idle=" << wait;
    cout << " total=" << total << "\n";
  }
  cout << "                *************************************\n\n";
  return (0);
}

// ##################################################################
// ##################################################################

double sim_control_NG_MPI::advance_step_OA1(const int l  ///< level to advance.
)
{
#ifdef TEST_INT
  cout << "NG-MPI advance_step_OA1, level=" << l << ", starting.\n";
#endif
  int err                   = 0;
  double dt2_fine           = 0.0;  // timestep for two finer level steps.
  double dt2_this           = 0.0;  // two timesteps for this level.
  class GridBaseClass* grid = SimPM.levels[l].grid;
  bool finest_level         = (l < (SimPM.grid_nlevels - 1)) ? false : true;

#ifdef TEST_INT
  cout << "advance_step_OA1: child=" << SimPM.levels[l].child << "\n";
  cout << "finest_level=" << finest_level << ", l=" << l
       << ", max=" << SimPM.grid_nlevels << "\n";
#endif

  // --------------------------------------------------------
  // 0. See if there are coarse-to-fine boundary data to send
  int c2f = -1, f2cs = -1, f2cr = -1;
  if (!finest_level) {
    // C2F data to send:
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_SEND) c2f = i;
    }
    // F2C data to receive
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_RECV) f2cr = i;
    }
  }
  if (l > 0) {
    // F2C data to send
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_SEND) f2cs = i;
    }
    // C2F data to recv can be more than one external BC, so find
    // them later in a loop.
  }
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " c2f = " << c2f;
  cout << ", f2c send=" << f2cs << ", f2c recv=" << f2cr << "\n";
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 0. Coarse to fine recv:
#ifdef C2F_FULLSTEP
  // only receive every 2nd step (every full step on coarse grid)
  if (SimPM.levels[l].step % 2 == 0) {
#endif
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " recv C2F BCs\n";
#endif
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(
            SimPM, spatial_solver, l, grid->BC_bd[i], SimPM.levels[l].step);
      }
    }
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " C2F CLEAR SEND\n";
#endif
    BC_COARSE_TO_FINE_SEND_clear_sends();
#ifdef C2F_FULLSTEP
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 1. Update external boundary conditions on level l
  // We have received interpolated data from the coarser level grid
  // already from the advance_step_OA1() for level l-1, if it exists.
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " update external BCs\n";
#endif
  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime, OA1, OA1);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 2. Send/receive external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are on a full step at the
  //    coarse-level.
  // --------------------------------------------------------
  if (c2f >= 0) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " C2F send\n";
#endif
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 2, 2);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 3. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " advance l+1 step 1\n";
#endif
    dt2_fine = advance_step_OA1(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 4. Calculate dU for this level (1st order)
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " calc DU\n";
#endif
  spatial_solver->Setdt(dt2_this);
  err += calc_microphysics_dU(dt2_this, grid);
  err += calc_dynamics_dU(dt2_this, OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt2_this, OA1, grid);
#endif  // THERMAL_CONDUCTION
  rep.errorTest("NG-MPI scn::advance_step_OA1: calc_x_dU", 0, err);
  if (l > 0) save_fine_fluxes(SimPM, l);
  if (l < SimPM.grid_nlevels - 1) save_coarse_fluxes(SimPM, l);
    // --------------------------------------------------------

    // --------------------------------------------------------
    // 5. Send external boundary data to finer grid (C2F)
    //    The (1,2) tells it that we are only half way
    //    through the coarse-level step.
    //    Then Receive the data on level l+1 and clear the C2F sends.
    // --------------------------------------------------------
#ifndef C2F_FULLSTEP
  if (c2f >= 0) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " C2F Send\n";
#endif
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 1, 2);
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 6. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " advance l+1 step 2\n";
#endif
    dt2_fine = advance_step_OA1(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 7. Update grid and boundaries on level l:
  //  - Receive level fluxes from finer grid (BC89)
  //  - update grid on level l to new time
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " update state vec\n";
#endif
  spatial_solver->Setdt(dt2_this);
#ifndef SKIP_BC89_FLUX
  if (l < SimPM.grid_nlevels - 1) {
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA1, OA1);
    rep.errorTest("scn::advance_step_OA1: recv_BC89_flux", 0, err);
  }
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt, OA1, OA1, grid);
  rep.errorTest("scn::advance_step_OA1: state-vec update", 0, err);
#ifndef SKIP_BC89_FLUX
  clear_sends_BC89_fluxes();
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 8. increment time and timestep for this level
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " \n";
#endif
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step++;
  if (l == SimPM.grid_nlevels - 1) {
    SimPM.timestep++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l == 0) SimPM.last_dt = SimPM.levels[l].dt;
    // --------------------------------------------------------

    // --------------------------------------------------------
    // 9. update internal boundary conditions on level l
    // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA1: l=" << l << " update internal BCS\n";
#endif
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime, OA1, OA1);
  //  - Recv F2C data from l+1
  if (!finest_level && f2cr >= 0) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " F2C Receive\n";
#endif
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA1, OA1);

    BC_FINE_TO_COARSE_SEND_clear_sends();
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 10. Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  if (grid->RT) {
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " RT at end of step\n";
#endif
    update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
    rep.errorTest("NG TIME_INT::update_RT_sources error", 0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    rep.errorTest("NG-MPI::advance_step_OA1: calc_rt_cols()", 0, err);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 11. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  if (l > 0 && SimPM.levels[l].step % 2 == 0) {
    // - send level fluxes
#ifndef SKIP_BC89_FLUX
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " \n";
#endif
    err += send_BC89_fluxes_F2C(SimPM, l, OA1, OA1);
#endif
#ifdef TEST_INT
    cout << "advance_step_OA1: l=" << l << " F2C SEND at tend of step\n";
#endif
    err += BC_update_FINE_TO_COARSE_SEND(
        SimPM, spatial_solver, l, grid->BC_bd[f2cs], OA1, OA1);
  }
  // --------------------------------------------------------

#ifdef TESTING
  cout << "NG-MPI advance_step_OA1, level=" << l << ", returning. t=";
  cout << SimPM.levels[l].simtime << ", step=" << SimPM.levels[l].step;
  cout << ", next dt=" << SimPM.levels[l].dt << " next time=";
  cout << SimPM.levels[l].simtime + SimPM.levels[l].dt << "\n";
#endif

  return dt2_this + SimPM.levels[l].dt;
}

// ##################################################################
// ##################################################################

double sim_control_NG_MPI::advance_step_OA2(const int l  ///< level to advance.
)
{
#ifdef TEST_INT
  cout.flush();
  cout << "NG-MPI advance_step_OA2, level=" << l << ", starting." << endl;
#endif
  int err                   = 0;
  double dt2_fine           = 0.0;  // timestep for two finer level steps.
  double dt2_this           = 0.0;  // two timesteps for this level.
  double ctime              = SimPM.levels[l].simtime;  // current time
  class GridBaseClass* grid = SimPM.levels[l].grid;
  bool finest_level         = (l < SimPM.grid_nlevels - 1) ? false : true;

#ifdef TEST_INT
  cout << "advance_step_OA2: child=" << SimPM.levels[l].child << endl;
#endif

  // --------------------------------------------------------
  // 0. See if there are coarse-to-fine boundary data to send
  int c2f = -1, f2cs = -1, f2cr = -1;
  if (!finest_level) {
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      // cout <<"i="<<i<<", BD type = "<<grid->BC_bd[i]->type<<"\n";
      // there's only one C2F boundary
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_SEND) c2f = i;
    }
    // F2C data to receive
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_RECV) f2cr = i;
    }
  }
  if (l > 0) {
    // F2C data to send
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == FINE_TO_COARSE_SEND) f2cs = i;
    }
    // C2F data to recv can be more than one external BC, so find
    // them later in a loop.
  }
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " c2f = " << c2f << endl;
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 0. Coarse to fine recv:
#ifdef C2F_FULLSTEP
  // only receive every 2nd step (every full step on coarse grid)
  if (SimPM.levels[l].step % 2 == 0) {
#endif
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " recv C2F BCs" << endl;
#endif
    for (size_t i = 0; i < grid->BC_bd.size(); i++) {
      if (grid->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
        err += BC_update_COARSE_TO_FINE_RECV(
            SimPM, spatial_solver, l, grid->BC_bd[i], SimPM.levels[l].step);
      }
    }
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " C2F CLEAR SEND" << endl;
#endif
    BC_COARSE_TO_FINE_SEND_clear_sends();
#ifdef C2F_FULLSTEP
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 1. Update external boundary conditions on level l
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " update external BCs" << endl;
#endif
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver, ctime, OA2, OA2);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 2. Send/receive external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are on a full step at the
  //    coarse-level.
  // --------------------------------------------------------
  if (c2f >= 0) {
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " C2F SEND" << endl;
#endif
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 2, 2);
  }

  // --------------------------------------------------------
  // 3. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " first fine step" << endl;
#endif
    dt2_fine = advance_step_OA2(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " dt=" << dt2_this << endl;
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 4. Calculate dU for this level (1st order predictor step)
  //    and advance to time-centred state at 0.5*dt
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " calc DU half step" << endl;
#endif
  double dt_now = dt2_this * 0.5;  // half of the timestep
  spatial_solver->Setdt(dt_now);
  err += calc_microphysics_dU(dt_now, grid);
  err += calc_dynamics_dU(dt_now, OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now, OA1, grid);
#endif  // THERMAL_CONDUCTION
  rep.errorTest("NG-MPI scn::advance_step_OA2: calc_x_dU", 0, err);

  // update state vector Ph to half-step values
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " update cell half step" << endl;
#endif
  err += grid_update_state_vector(dt_now, OA1, OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA2", 0, err);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 5. Update grid and boundary data on this level at 1/2 step
  // - update internal and F2C boundaries
  // - update external boundaries
  // - don't update C2F (b/c it is 1/4 step on l-1 level)
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " update boundaries" << endl;
#endif
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, OA1, OA2);

  if (!finest_level && f2cr >= 0) {
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA1, OA2);

    BC_FINE_TO_COARSE_SEND_clear_sends();
  }

  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, OA1, OA2);
  rep.errorTest("scn::advance_step_OA2: bounday update OA2", 0, err);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 6. Calculate dU for the full step (OA2) on this level
  // --------------------------------------------------------
  dt_now = dt2_this;  // full step
  spatial_solver->Setdt(dt_now);
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " raytracing" << endl;
#endif
  if (grid->RT) {
    update_evolving_RT_sources(
        SimPM, SimPM.levels[l].simtime + 0.5 * dt_now, grid->RT);
    rep.errorTest("NG TIME_INT::update_RT_sources error", 0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    rep.errorTest("scn::advance_time: calc_rt_cols() OA2", 0, err);
  }
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " full step calc dU" << endl;
#endif
  err += calc_microphysics_dU(dt_now, grid);
  err += calc_dynamics_dU(dt_now, OA2, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now, OA2, grid);
#endif  // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA2", 0, err);
  if (l > 0) save_fine_fluxes(SimPM, l);
  if (l < SimPM.grid_nlevels - 1) save_coarse_fluxes(SimPM, l);
    // --------------------------------------------------------

    // --------------------------------------------------------
    // 7. Send/receive external boundary data to finer grid (C2F)
    //    The (1,2) tells it that we are half-way through
    //    the coarse-level step.
    // --------------------------------------------------------
#ifndef C2F_FULLSTEP
  if (c2f >= 0) {
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 1, 2);
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 8. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " second fine step" << endl;
#endif
    dt2_fine = advance_step_OA2(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 9. Update local grid to new time values
  //  - Receive level fluxes from finer grid (FLUX)
  //  - update grid on level l to new time
  // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " full step update" << endl;
#endif
  spatial_solver->Setdt(dt_now);
#ifndef SKIP_BC89_FLUX
  if (l < SimPM.grid_nlevels - 1) {
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA2, OA2);
    rep.errorTest("scn::advance_step_OA1: recv_BC89_flux", 0, err);
  }
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt, OA2, OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update", 0, err);
#ifndef SKIP_BC89_FLUX
  clear_sends_BC89_fluxes();
#endif
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " updated, sends cleared" << endl;
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 10. increment time and timestep for this level
  // --------------------------------------------------------
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step++;
  if (l == SimPM.grid_nlevels - 1) {
    SimPM.timestep++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l == 0) SimPM.last_dt = SimPM.levels[l].dt;
    // --------------------------------------------------------

    // --------------------------------------------------------
    // 11. update internal boundary conditions on level l
    // --------------------------------------------------------
#ifdef TEST_INT
  cout << "advance_step_OA2: l=" << l << " update internal BCs" << endl;
#endif
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, OA2, OA2);
  //  - Recv F2C data from l+1
  if (!finest_level && f2cr >= 0) {
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " F2C recv start" << endl;
#endif
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA2, OA2);
    //  - Clear F2C sends
    BC_FINE_TO_COARSE_SEND_clear_sends();
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " F2C recv done" << endl;
#endif
  }

  // --------------------------------------------------------
  // 12. Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  if (grid->RT) {
    update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
    rep.errorTest("NG TIME_INT::update_RT_sources error", 0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    rep.errorTest("NG-MPI::advance_step_OA2: raytracing()", 0, err);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 13. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  if (l > 0) {
#ifdef TEST_INT
    cout << "step=" << SimPM.levels[l].step << endl;
#endif

#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " send BC89 fluxes" << endl;
#endif
    if (SimPM.levels[l].step % 2 == 0) {
      // only send level fluxes every 2nd step (coarse grid is only
      // updated at the full-step, not at the half-step).
#ifndef SKIP_BC89_FLUX
      err += send_BC89_fluxes_F2C(SimPM, l, OA2, OA2);
#endif
    }
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " sent BC89 fluxes" << endl;
#endif

#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " update F2C" << endl;
#endif
    err += BC_update_FINE_TO_COARSE_SEND(
        SimPM, spatial_solver, l, grid->BC_bd[f2cs], OA2, OA2);
#ifdef TEST_INT
    cout << "advance_step_OA2: l=" << l << " updated F2C" << endl;
#endif
  }
  // --------------------------------------------------------

#ifdef TESTING
  cout << "NG-MPI advance_step_OA2, level=" << l << ", returning. t=";
  cout << SimPM.levels[l].simtime << ", step=" << SimPM.levels[l].step;
  cout << ", next dt=" << SimPM.levels[l].dt << " next time=";
  cout << SimPM.levels[l].simtime + SimPM.levels[l].dt << endl;
#endif
  return dt2_this + SimPM.levels[l].dt;
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::initial_conserved_quantities(
    vector<class GridBaseClass*>& grid)
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION
  pion_flt u[SimPM.nvar];
  // dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange =
  // 0.0;
  //  cout <<"initERG: "<<dp.initERG<<"\n";
  initERG = 0.;
  initMMX = initMMY = initMMZ = 0.;
  initMASS                    = 0.0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    double dx     = SimPM.levels[l].dx;
    double dv     = 0.0;
    class cell* c = grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        // cout <<"*** LEVEL "<<l<<", cell is a leaf: "<<c->isdomain<<"
        // "<<c->isleaf<<": "; rep.printVec("pos",c->pos,SimPM.ndim);
        dv = spatial_solver->CellVolume(c, dx);
        spatial_solver->PtoU(c->P, u, SimPM.gamma);
        initERG += u[ERG] * dv;
        initMMX += u[MMX] * dv;
        initMMY += u[MMY] * dv;
        initMMZ += u[MMZ] * dv;
        initMASS += u[RHO] * dv;
      }
      else {
        // cout <<"level "<<l<<", cell is not leaf: "<<c->isdomain<<"
        // "<<c->isleaf<<": "; rep.printVec("pos",c->pos,SimPM.ndim);
      }
    } while ((c = grid[l]->NextPt(c)) != 0);
  }

  // cout <<"(local quantities) ["<< initERG <<", ";
  // cout << initMMX <<", ";
  // cout << initMMY <<", ";
  // cout << initMMZ <<", ";
  // cout << initMASS <<"]\n";

  initERG  = COMM->global_operation_double("SUM", initERG);
  initMMX  = COMM->global_operation_double("SUM", initMMX);
  initMMY  = COMM->global_operation_double("SUM", initMMY);
  initMMZ  = COMM->global_operation_double("SUM", initMMZ);
  initMASS = COMM->global_operation_double("SUM", initMASS);

  cout << "(conserved quantities) [" << initERG << ", ";
  cout << initMMX << ", ";
  cout << initMMY << ", ";
  cout << initMMZ << ", ";
  cout << initMASS << "]\n";

#endif  // TEST_CONSERVATION
  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::check_energy_cons(vector<class GridBaseClass*>& grid)
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION
  pion_flt u[SimPM.nvar];
  nowERG        = 0.;
  nowMMX        = 0.;
  nowMMY        = 0.;
  nowMMZ        = 0.;
  nowMASS       = 0.0;
  double totmom = 0.0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    double dx     = SimPM.levels[l].dx;
    double dv     = 0.0;
    class cell* c = grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        dv = spatial_solver->CellVolume(c, dx);
        spatial_solver->PtoU(c->P, u, SimPM.gamma);
        nowERG += u[ERG] * dv;
        nowMMX += u[MMX] * dv;
        nowMMY += u[MMY] * dv;
        nowMMZ += u[MMZ] * dv;
        nowMASS += u[RHO] * dv;
        totmom +=
            sqrt(u[MMX] * u[MMX] + u[MMY] * u[MMY] + u[MMZ] * u[MMZ]) * dv;
      }
    } while ((c = grid[l]->NextPt(c)) != 0);
  }

  // cout <<"(local quantities) ["<< nowERG <<", ";
  // cout << nowMMX <<", ";
  // cout << nowMMY <<", ";
  // cout << nowMMZ <<", ";
  // cout << nowMASS <<"]\n";

  nowERG  = COMM->global_operation_double("SUM", nowERG);
  nowMMX  = COMM->global_operation_double("SUM", nowMMX);
  nowMMY  = COMM->global_operation_double("SUM", nowMMY);
  nowMMZ  = COMM->global_operation_double("SUM", nowMMZ);
  nowMASS = COMM->global_operation_double("SUM", nowMASS);
  totmom  = COMM->global_operation_double("SUM", totmom);
  // cout <<" totmom="<<totmom<<" initMMX="<<initMMX;
  // cout <<", nowMMX="<<nowMMX<<"\n";

  cout << "(conserved quantities) [" << nowERG << ", ";
  cout << nowMMX << ", ";
  cout << nowMMY << ", ";
  cout << nowMMZ << ", ";
  cout << nowMASS << "]\n";
  cout << "(relative error      ) [";
  cout << (nowERG - initERG) / (initERG) << ", ";
  cout << (nowMMX - initMMX) / (totmom) << ", ";
  cout << (nowMMY - initMMY) / (totmom) << ", ";
  cout << (nowMMZ - initMMZ) / (totmom) << ", ";
  cout << (nowMASS - initMASS) / initMASS << "]\n";

#endif  // TEST_CONSERVATION
  return (0);
}

// ##################################################################
// ##################################################################

#endif  // PARALLEL
