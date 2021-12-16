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

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "constants.h"

#include "raytracing/raytracer_SC.h"
#include "sim_control/sim_control_NG_MPI.h"
#include "sub_domain/sub_domain.h"

#include <fstream>
#include <iomanip>

#include <sstream>
using namespace std;

#ifdef PARALLEL

//#define TEST_INT

// ##################################################################
// ##################################################################

sim_control_NG_MPI::sim_control_NG_MPI()
{
  spdlog::debug("sim_control_NG_MPI constructor");
}

// ##################################################################
// ##################################################################

sim_control_NG_MPI::~sim_control_NG_MPI()
{
  spdlog::debug("sim_control_NG_MPI destructor");
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::Init(
    string infile,
    int typeOfFile,
    int narg,
    string *args,
    vector<class GridBaseClass *>
        &grid  ///< address of vector of grid pointers.
)
{
  spdlog::debug("(pion) Init: infile = {}", infile);
  int err = 0;

  // ----------------------------------------------------------------
  SimPM.typeofip = typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  if (!dataio->file_exists(infile))
    spdlog::error("{}: {}", "infile doesn't exist!", infile);

  // ----------------------------------------------------------------
  err = dataio->ReadHeader(infile, SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI Init(): failed to read header", 0,
        err);

  // Check if any commandline args override the file parameters.
  // ----------------------------------------------------------------
  err = override_params(narg, args);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_MPI INIT::override_params)", 0, err);

  // setup the nested grid levels, and decompose the domain on each
  // level
  // ----------------------------------------------------------------
  setup_NG_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
  spdlog::info("NG_MPI Init: grid setup");
  err = setup_grid(grid, SimPM);
  spdlog::info("NG_MPI Init: grid setup finished");
  SimPM.dx = grid[0]->DX();
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_MPI INIT::setup_grid) error", 0,
        err);

  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  // ----------------------------------------------------------------
  err = set_equations(SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_MPI INIT::set_equations)", 0, err);
  spatial_solver->SetEOS(SimPM.gamma);

  // set up Microphysics, if needed.
  // ----------------------------------------------------------------
  err = setup_microphysics(SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_MPI INIT::setup_microphysics)", 0,
        err);

  // assign data to the grid from snapshot file.
  // ----------------------------------------------------------------
  err = dataio->ReadData(infile, grid, SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_MPI INIT::assign_initial_data)", 0,
        err);

  // ----------------------------------------------------------------
  // Set Ph[] = P[], and then implement the boundary conditions.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    cell *c = grid[l]->FirstPt();
    std::vector<double> u(SimPM.nvar);
    do {
      // make sure temperature of the gas is reasonable
      if (SimPM.timestep == 0) {
        spatial_solver->PtoU(c->P, u.data(), SimPM.gamma);
        spatial_solver->UtoP(
            u.data(), c->P, SimPM.EP.MinTemperature, SimPM.gamma);
      }
      // set Ph[]=P[]
      for (int v = 0; v < SimPM.nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = grid[l]->NextPt(c)) != 0);

    if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
      spdlog::info("Initial state, zero-ing glm variable");
      c = grid[l]->FirstPt();
      do {
        c->P[SI] = c->Ph[SI] = 0.;
      } while ((c = grid[l]->NextPt(c)) != 0);
    }
  }  // loop over levels

  // ----------------------------------------------------------------
  err = boundary_conditions(SimPM, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "(NG_MPI INIT::boundary_conditions) err!=0", 0, err);

  // ----------------------------------------------------------------
  err += setup_raytracing(SimPM, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "Failed to setup raytracer", 0, err);

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err = assign_boundary_data(SimPM, l, grid[l], MP);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG_MPI INIT::assign_boundary_data", 0,
          err);
    SimPM.levels[0].sub_domain.barrier("level assign boundary data");
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "NG_MPI updating external boundaries for level {0}\nUPDATING EXTERNAL BOUNDARIES FOR LEVEL {0}",
        l);
    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI INIT: error from bounday update",
        0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "NG_MPI updating C2F boundaries for level {0}\nUPDATING C2F BOUNDARIES FOR LEVEL {0}",
        l);
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
        spdlog::debug(
            "Init: l={}, C2F recv i={}, type={}", l, i,
            grid[l]->BC_bd[i]->type);
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i],
              SimPM.levels[l].step);
        }
      }
#ifndef NDEBUG
      spdlog::debug("CLEAR C2F send from {} to {}...", l - 1, l);
#endif
      BC_COARSE_TO_FINE_SEND_clear_sends(SimPM.levels[l - 1].sub_domain);
#ifndef NDEBUG
      spdlog::debug("  ... done");
#endif
    }
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI INIT: error from bounday update",
        0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "NG_MPI updating external boundaries for level {0}@@@@@@@@@@@@  UPDATING EXTERNAL BOUNDARIES FOR LEVEL {0}",
        l);
    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI INIT: error from bounday update",
        0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    spdlog::debug(
        "NG_MPI updating internal boundaries for level {0}\n@@@@@@@@@@@@  UPDATING INTERNAL BOUNDARIES FOR LEVEL {0}",
        l);
    err += TimeUpdateInternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, 0.0, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI INIT: error from bounday update",
        0, err);

  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    spdlog::debug("NG_MPI updating F2C boundaries for level {}", l);
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
#ifndef NDEBUG
      spdlog::debug("CLEAR F2C send from {} to {}...", l + 1, l);
#endif
      BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[l + 1].sub_domain);
#ifndef NDEBUG
      spdlog::debug("  ... done");
#endif
    }
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG_MPI INIT: error from bounday update",
        0, err);
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
      spdlog::error(
          "{}: {}", "opfreq_time not set right and is needed!",
          SimPM.opfreq_time);
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
    if (!dataio)
      spdlog::error("{}: {}", "NG_MPI INIT:: dataio", SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  dataio->SetMicrophysics(MP);
  if (textio) {
    textio->SetSolver(spatial_solver);
    textio->SetMicrophysics(MP);
  }

  if (SimPM.timestep == 0) {
    spdlog::info("(NG_MPI INIT) Writing initial data");
    output_data(grid);
  }

  // ----------------------------------------------------------------
  //#ifndef NDEBUG
  cell *c = 0;
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
  //#endif // NDEBUG
  return (err);
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::Time_Int(
    vector<class GridBaseClass *> &grid  ///< grid pointers.
)
{
  spdlog::info("(pion-ng-mpi)  STARTING TIME INTEGRATION");
  int err       = 0;
  SimPM.maxtime = false;

  bool first_step = true;
  bool restart    = true;
  if (SimPM.timestep != 0) first_step = false;
  if (SimPM.timestep == 0) restart = false;

  // start all the timers for identifying performance bottlenecks.
  clk.start_timer("time_int");
  double tsf            = 0.0;
  vector<string> timing = {"dt",   "ibc", "ebc", "f2c", "c2f",
                           "bc89", "rt",  "dyn", "mp",  "upd"};
  for (auto i : timing) {
    clk.start_timer(i);
    clk.pause_timer(i);
  }

  // make sure all levels start at the same time.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt      = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }

  // --------------------------------------------------------------
  // Update internal and external boundaries.

  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    spdlog::debug("updating internal boundaries for level {}", l);
    err += TimeUpdateInternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.levels[l].simtime, 0.0,
        SimPM.tmOOA, SimPM.tmOOA);
    spdlog::debug("... done");
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "sim_control_NG_MPI: internal boundary",
        0, err);
  // --------------------------------------------------------------
  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    if (l < SimPM.grid_nlevels - 1) {
      spdlog::debug("NG_MPI Receiving F2C boundaries for level {}", l);
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += BC_update_FINE_TO_COARSE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
      BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[l + 1].sub_domain);
#ifdef TEST_INT
      spdlog::debug("NG_MPI F2C cleared sends from l={}", l + 1);
#endif
    }

    spdlog::debug("NG_MPI raytracing level {}", l);
    err =
        update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid[l]->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error",
          0, err);
    do_ongrid_raytracing(SimPM, grid[l], l);

    if (l > 0) {
      spdlog::debug("NG_MPI Sending F2C boundaries for level {}", l);
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += BC_update_FINE_TO_COARSE_SEND(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
#ifdef TEST_INT
  spdlog::info("NG_MPI updated F2C boundaries for all levels");
#endif
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "NG_MPI time-int: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // update coarse-to-fine level boundaries
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "NG_MPI updating C2F boundaries for level {0}\nupdating external boundaries for level {0}",
        l);
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid[l]->BC_bd[i],
              SimPM.levels[l].step);
        }
      }
      BC_COARSE_TO_FINE_SEND_clear_sends(SimPM.levels[l - 1].sub_domain);
    }

    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.levels[l].simtime, SimPM.tmOOA,
        SimPM.tmOOA);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "sim_control_NG_MPI: external boundary",
          0, err);

    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += BC_update_COARSE_TO_FINE_SEND(
              SimPM, grid[l], spatial_solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "NG_MPI time-int: error from bounday update", 0, err);
  // ----------------------------------------------------------------

  cout.setf(ios_base::scientific);
  cout.precision(7);
  while (SimPM.maxtime == false) {
    // --------------------------------------------------------------
    // Get timestep on each level
    int scale    = 1;
    double mindt = 1.0e99;

    clk.start_timer("dt");
    for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
      spdlog::debug(
          "Calculate timestep, level {}, dx={}", l, SimPM.levels[l].dx);

      if (!restart && !first_step) {
        SimPM.last_dt = SimPM.levels[l].last_dt;
      }
      else {
        SimPM.levels[l].last_dt = SimPM.last_dt / SimPM.levels[l].multiplier;
      }
      err += calculate_timestep(SimPM, grid[l], l);
      if (0 != err)
        spdlog::error(
            "{}: Expected {} but got {}", "TIME_INT::calc_timestep()", 0, err);

      mindt = std::min(mindt, SimPM.dt / scale);
      mindt = SimPM.levels[0].sub_domain.global_operation_double("MIN", mindt);
      spdlog::debug(
          "level {} got dt={} and {}... mindt={}", l, SimPM.dt,
          SimPM.dt / scale, mindt);
      SimPM.levels[l].dt = SimPM.dt;
      scale *= 2;
    }
    // make sure all levels use same step (scaled by factors of 2).
    scale = 1;
    for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
      SimPM.levels[l].dt = mindt * scale;
      scale *= 2;
      spdlog::debug(
          "new dt={}, t={}", SimPM.levels[l].dt, SimPM.levels[l].simtime);
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
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "TIME_INT::calc_timestep()", 0, err);
    clk.pause_timer("dt");
    // --------------------------------------------------------------

    // --------------------------------------------------------------
    spdlog::debug("NG_MPI time_int: stepping forward in time");
    // Use a recursive algorithm to update the coarsest level.  This
    // function also updates the next level twice, by calling itself
    // for the finer level, and so on.
    //
    advance_time(0, grid[0]);
    SimPM.simtime = SimPM.levels[0].simtime;
    SimPM.levels[0].sub_domain.barrier("step");
    spdlog::debug("MPI time_int: finished timestep");

    if (SimPM.levels[0].sub_domain.get_myrank() == 0) {
      tsf = clk.time_so_far("time_int");
      spdlog::info(
          "New time: {:12.6e}   dt: {:12.6e}   steps: {:8d}   l0 steps: {:6d}   runtime: {:12.4e} s",
          SimPM.simtime, SimPM.levels[SimPM.grid_nlevels - 1].dt,
          SimPM.timestep,
          SimPM.timestep / static_cast<int>(pow(2, SimPM.grid_nlevels - 1)),
          tsf);
      // cout <<"\tTimings: ";
      // for (auto i : timing) {
      //  //cout << i <<"  "<< clk.time_so_far_paused(i) <<"  ";
      //  cout << clk.time_so_far_paused(i) <<"  ";
      //}
      // cout <<"\n";
      // spdlog::debug("\t runtime: {}s", tsf);
    }
    // --------------------------------------------------------------

    err += output_data(grid);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "MPI_NG TIME_INT::output_data()", 0,
          err);

#ifdef TEST_CONSERVATION
    err += check_energy_cons(grid);
#endif

    //
    // check if we are at time limit yet.
    //
    tsf = clk.time_so_far("time_int");
    double maxt =
        SimPM.levels[0].sub_domain.global_operation_double("MAX", tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime = true;
      spdlog::debug("RUNTIME>{} SECS", get_max_walltime());
    }

    err += check_eosim();
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "MPI_NG TIME_INT::check_eosim()", 0,
          err);
  }
  spdlog::info(
      "sim_control_NG_MPI:: TIME_INT FINISHED.  MOVING ON TO FINALISE SIM");
  tsf = clk.time_so_far("time_int");
  spdlog::info(
      "TOTALS ###: Nsteps={}, sim-time={}, wall-time={}, time/step={}",
      SimPM.timestep, SimPM.simtime, tsf,
      tsf / static_cast<double>(SimPM.timestep));
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
    spdlog::info("TOTALS RT#: active={} idle={} total={}", run, wait, total);
  }
  ostringstream tm;
  tm << "#";
  for (auto i : timing) {
    tm << setw(11) << i;
  }
  tm << setw(11) << "sum";
  spdlog::info(tm.str());
  tm.str("");
  double t = 0.0;
  tm << "   ";
  tm.setf(ios_base::scientific);
  tm.precision(3);
  for (auto i : timing) {
    tsf = clk.time_so_far_paused(i);
    t += tsf;
    tm << tsf << "  ";
  }
  tm << t;
  spdlog::info(tm.str());
  spdlog::info("                *************************************\n");
  return (0);
}



// ##################################################################
// ##################################################################



double sim_control_NG_MPI::advance_step_OA1(const int l  ///< level to advance.
)
{
  spdlog::debug("NG-MPI advance_step_OA1, level={}, starting", l);
  int err                   = 0;
  double dt2_this           = 0.0;  // two timesteps for this level.
  class GridBaseClass *grid = SimPM.levels[l].grid;
  bool finest_level         = (l < (SimPM.grid_nlevels - 1)) ? false : true;

  spdlog::debug(
      "advance_step_OA1: child={}\nfinest_level={}, l={}, max={}",
      fmt::ptr(SimPM.levels[l].child), finest_level, l, SimPM.grid_nlevels);

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
  spdlog::debug(
      "advance_step_OA1: l={} c2f = {}, f2c send={}, f2c recv={}", l, c2f, f2cs,
      f2cr);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 0. Coarse to fine recv:
#ifdef C2F_FULLSTEP
  // only receive every 2nd step (every full step on coarse grid)
  if (SimPM.levels[l].step % 2 == 0) {
#endif
    if (l > 0) {
#ifdef TEST_INT
      spdlog::debug("advance_step_OA1: l={} recv C2F BCs", l);
#endif
      for (size_t i = 0; i < grid->BC_bd.size(); i++) {
        if (grid->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid->BC_bd[i], SimPM.levels[l].step);
        }
      }
#ifdef TEST_INT
      spdlog::debug("advance_step_OA1: l={} C2F CLEAR SEND", l);
#endif
      BC_COARSE_TO_FINE_SEND_clear_sends(SimPM.levels[l - 1].sub_domain);
    }
#ifdef C2F_FULLSTEP
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 1. Update external boundary conditions on level l
  // We have received interpolated data from the coarser level grid
  // already from the advance_step_OA1() for level l-1, if it exists.
  spdlog::debug("advance_step_OA1: l={} update external BCs", l);
  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime, OA1, OA1);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 2. Send/receive external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are on a full step at the
  //    coarse-level.
  // --------------------------------------------------------
  if (c2f >= 0) {
    spdlog::debug("advance_step_OA1: l={} C2F send", l);
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 2, 2);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 3. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
    spdlog::debug("advance_step_OA1: l={} advance l+1 step 1", l);
    advance_step_OA1(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 4. Calculate dU for this level (1st order)
  // --------------------------------------------------------
  spdlog::debug("advance_step_OA1: l={} calc DU", l);
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt2_this);
#ifdef PION_OMP
  }
#endif
  err += calc_microphysics_dU(dt2_this, grid);
  err += calc_dynamics_dU(dt2_this, OA1, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG-MPI scn::advance_step_OA1: calc_x_dU",
        0, err);
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
    spdlog::debug("advance_step_OA1: l={} C2F Send", l);
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 1, 2);
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 6. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
    spdlog::debug("advance_step_OA1: l={} advance l+1 step 2", l);
    advance_step_OA1(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 7. Update grid and boundaries on level l:
  //  - Receive level fluxes from finer grid (BC89)
  //  - update grid on level l to new time
  // --------------------------------------------------------
  spdlog::debug("advance_step_OA1: l={} update state vec", l);
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt2_this);
#ifdef PION_OMP
  }
#endif
#ifndef SKIP_BC89_FLUX
  if (l < SimPM.grid_nlevels - 1) {
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA1, OA1);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "scn::advance_step_OA1: recv_BC89_flux",
          0, err);
  }
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt, OA1, OA1, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA1: state-vec update",
        0, err);
#ifndef SKIP_BC89_FLUX
  if (l < SimPM.grid_nlevels - 1) {
    clear_sends_BC89_fluxes(SimPM.levels[l + 1].sub_domain);
  }
#endif
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 8. increment time and timestep for this level
  // --------------------------------------------------------
  spdlog::debug("advance_step_OA1: l={}", l);
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
  spdlog::debug("advance_step_OA1: l={} update internal BCS\n", l);
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime,
      SimPM.levels[l].dt, OA1, OA1);
  //  - Recv F2C data from l+1
  if (!finest_level && f2cr >= 0) {
    spdlog::debug("advance_step_OA1: l={} F2C Receive", l);
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA1, OA1);

    BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[l + 1].sub_domain);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 10. Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  if (grid->RT) {
    spdlog::debug("advance_step_OA1: l={} RT at end of step", l);
    update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error",
          0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG-MPI::advance_step_OA1: calc_rt_cols()", 0, err);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 11. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  if (l > 0 && SimPM.levels[l].step % 2 == 0) {
    // - send level fluxes
#ifndef SKIP_BC89_FLUX
    spdlog::debug("advance_step_OA1: l={}", l);
    err += send_BC89_fluxes_F2C(SimPM, l, OA1, OA1);
#endif
    spdlog::debug("advance_step_OA1: l={} F2C SEND at tend of step", l);
    err += BC_update_FINE_TO_COARSE_SEND(
        SimPM, spatial_solver, l, grid->BC_bd[f2cs], OA1, OA1);
  }
  // --------------------------------------------------------

  spdlog::debug(
      "NG-MPI advance_step_OA1, level={}, returning. t={}, step={}, next dt={} next time={}",
      l, SimPM.levels[l].simtime, SimPM.levels[l].step, SimPM.levels[l].dt,
      SimPM.levels[l].simtime + SimPM.levels[l].dt);

  return dt2_this + SimPM.levels[l].dt;
}



// ##################################################################
// ##################################################################



double sim_control_NG_MPI::advance_step_OA2(const int l  ///< level to advance.
)
{
  spdlog::debug("NG-MPI advance_step_OA2, level={}, starting", l);
  int err                   = 0;
  double dt2_this           = 0.0;  // two timesteps for this level.
  double ctime              = SimPM.levels[l].simtime;  // current time
  class GridBaseClass *grid = SimPM.levels[l].grid;
  bool finest_level         = (l < SimPM.grid_nlevels - 1) ? false : true;

  spdlog::debug("advance_step_OA2: child={}", fmt::ptr(SimPM.levels[l].child));

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
  spdlog::debug("advance_step_OA2: l={} c2f = {}", l, c2f);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 0. Coarse to fine recv:
  clk.start_timer("c2f");
#ifdef C2F_FULLSTEP
  // only receive every 2nd step (every full step on coarse grid)
  if (SimPM.levels[l].step % 2 == 0) {
#endif
    if (l > 0) {
#ifdef TEST_INT
      spdlog::debug("advance_step_OA2: l={} recv C2F BCs", l);
#endif
      for (size_t i = 0; i < grid->BC_bd.size(); i++) {
        if (grid->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += BC_update_COARSE_TO_FINE_RECV(
              SimPM, spatial_solver, l, grid->BC_bd[i], SimPM.levels[l].step);
        }
      }
#ifdef TEST_INT
      spdlog::debug("advance_step_OA2: l={} C2F CLEAR SEND", l);
#endif
      BC_COARSE_TO_FINE_SEND_clear_sends(SimPM.levels[l - 1].sub_domain);
    }
#ifdef C2F_FULLSTEP
  }
#endif
  clk.pause_timer("c2f");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 1. Update external boundary conditions on level l
  // --------------------------------------------------------
  clk.start_timer("ebc");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} update external BCs", l);
#endif
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver, ctime, OA2, OA2);
  clk.pause_timer("ebc");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 2. Send/receive external boundary data to finer grid (C2F)
  //    The (2,2) tells it that we are on a full step at the
  //    coarse-level.
  // --------------------------------------------------------
  clk.start_timer("c2f");
  if (c2f >= 0) {
    spdlog::debug("advance_step_OA2: l={} C2F SEND", l);
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 2, 2);
  }
  clk.pause_timer("c2f");

  // --------------------------------------------------------
  // 3. advance finer-level by one step, if it exists
  // --------------------------------------------------------
  if (!finest_level) {
    spdlog::debug("advance_step_OA2: l={} first fine step", l);
    advance_step_OA2(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
  spdlog::debug("advance_step_OA2: l={} dt={}", l, dt2_this);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 4. Calculate dU for this level (1st order predictor step)
  //    and advance to time-centred state at 0.5*dt
  // --------------------------------------------------------
  spdlog::debug("advance_step_OA2: l={} calc DU half step", l);
  double dt_now = dt2_this * 0.5;  // half of the timestep
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt_now);
#ifdef PION_OMP
  }
#endif
  clk.start_timer("mp");
  err += calc_microphysics_dU(dt_now, grid);
  clk.pause_timer("mp");
  clk.start_timer("dyn");
  err += calc_dynamics_dU(dt_now, OA1, grid);
  clk.pause_timer("dyn");
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "NG-MPI scn::advance_step_OA2: calc_x_dU",
        0, err);

  // update state vector Ph to half-step values
  clk.start_timer("upd");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} update cell half step", l);
#endif
  err += grid_update_state_vector(dt_now, OA1, OA2, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "scn::advance_step_OA2: state-vec update OA2", 0, err);
  clk.pause_timer("upd");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 5. Update grid and boundary data on this level at 1/2 step
  // - update internal and F2C boundaries
  // - update external boundaries
  // - don't update C2F (b/c it is 1/4 step on l-1 level)
  // --------------------------------------------------------
  clk.start_timer("ibc");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} update boundaries", l);
#endif
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, dt_now, OA1, OA2);
  clk.pause_timer("ibc");

  clk.start_timer("f2c");
  if (!finest_level && f2cr >= 0) {
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA1, OA2);

    BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[l + 1].sub_domain);
  }
  clk.pause_timer("f2c");

  clk.start_timer("ebc");
  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, OA1, OA2);
  clk.pause_timer("ebc");
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "scn::advance_step_OA2: bounday update OA2", 0, err);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 6. Calculate dU for the full step (OA2) on this level
  // --------------------------------------------------------
  dt_now = dt2_this;  // full step
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt_now);
#ifdef PION_OMP
  }
#endif
  clk.start_timer("rt");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} raytracing", l);
#endif
  if (grid->RT) {
    update_evolving_RT_sources(
        SimPM, SimPM.levels[l].simtime + 0.5 * dt_now, grid->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error",
          0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "scn::advance_time: calc_rt_cols() OA2",
          0, err);
  }
  clk.pause_timer("rt");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} full step calc dU", l);
#endif
  clk.start_timer("mp");
  err += calc_microphysics_dU(dt_now, grid);
  clk.pause_timer("mp");
  clk.start_timer("dyn");
  err += calc_dynamics_dU(dt_now, OA2, grid);
  clk.pause_timer("dyn");
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: calc_x_dU OA2", 0,
        err);
  clk.start_timer("bc89");
  if (l > 0) save_fine_fluxes(SimPM, l);
  if (l < SimPM.grid_nlevels - 1) save_coarse_fluxes(SimPM, l);
  clk.pause_timer("bc89");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 7. Send/receive external boundary data to finer grid (C2F)
  //    The (1,2) tells it that we are half-way through
  //    the coarse-level step.
  // --------------------------------------------------------
  clk.start_timer("c2f");
#ifndef C2F_FULLSTEP
  if (c2f >= 0) {
    err += BC_update_COARSE_TO_FINE_SEND(
        SimPM, grid, spatial_solver, l, grid->BC_bd[c2f], 1, 2);
  }
#endif
  clk.pause_timer("c2f");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 8. Take another step on finer grid
  // --------------------------------------------------------
  if (!finest_level) {
    spdlog::debug("advance_step_OA2: l={} second fine step", l);
    advance_step_OA2(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 9. Update local grid to new time values
  //  - Receive level fluxes from finer grid (FLUX)
  //  - update grid on level l to new time
  // --------------------------------------------------------
  spdlog::debug("advance_step_OA2: l={} full step update", l);
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt_now);
#ifdef PION_OMP
  }
#endif
#ifndef SKIP_BC89_FLUX
  clk.start_timer("bc89");
  if (l < SimPM.grid_nlevels - 1) {
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA2, OA2);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "scn::advance_step_OA1: recv_BC89_flux",
          0, err);
  }
  clk.pause_timer("bc89");
#endif
  clk.start_timer("upd");
  err += grid_update_state_vector(SimPM.levels[l].dt, OA2, OA2, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: state-vec update",
        0, err);
  clk.pause_timer("upd");
#ifndef SKIP_BC89_FLUX
  clk.start_timer("bc89");
  if (l < SimPM.grid_nlevels - 1) {
    clear_sends_BC89_fluxes(SimPM.levels[l + 1].sub_domain);
  }
  clk.pause_timer("bc89");
#endif
  spdlog::debug("advance_step_OA2: l={} updated, sends cleared", l);
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
  clk.start_timer("ibc");
#ifdef TEST_INT
  spdlog::debug("advance_step_OA2: l={} update internal BCs", l);
#endif
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, 0.5 * dt_now, OA2, OA2);
  clk.pause_timer("ibc");

  //  - Recv F2C data from l+1
  clk.start_timer("f2c");
  if (!finest_level && f2cr >= 0) {
    spdlog::debug("advance_step_OA2: l={} F2C recv start", l);
    err += BC_update_FINE_TO_COARSE_RECV(
        SimPM, spatial_solver, l, grid->BC_bd[f2cr], OA2, OA2);
    //  - Clear F2C sends
    BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[l + 1].sub_domain);
#ifdef TEST_INT
    spdlog::debug("advance_step_OA2: l={} F2C recv done", l);
#endif
  }
  clk.pause_timer("f2c");

  // --------------------------------------------------------
  // 12. Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  clk.start_timer("rt");
  if (grid->RT) {
    update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error",
          0, err);
    err += do_ongrid_raytracing(SimPM, grid, l);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG-MPI::advance_step_OA2: raytracing()", 0, err);
  }
  clk.pause_timer("rt");
  // --------------------------------------------------------

  // --------------------------------------------------------
  // 13. Send level fluxes and F2C data to coarser grid
  // --------------------------------------------------------
  if (l > 0) {
    spdlog::debug("step={}", SimPM.levels[l].step);

    clk.start_timer("bc89");
#ifdef TEST_INT
    spdlog::debug("advance_step_OA2: l={} send BC89 fluxes", l);
#endif
    if (SimPM.levels[l].step % 2 == 0) {
      // only send level fluxes every 2nd step (coarse grid is only
      // updated at the full-step, not at the half-step).
#ifndef SKIP_BC89_FLUX
      err += send_BC89_fluxes_F2C(SimPM, l, OA2, OA2);
#endif
    }
    clk.pause_timer("bc89");
    spdlog::debug("advance_step_OA2: l={} sent BC89 fluxes", l);

    spdlog::debug("advance_step_OA2: l={} update F2C", l);
    clk.start_timer("f2c");
    err += BC_update_FINE_TO_COARSE_SEND(
        SimPM, spatial_solver, l, grid->BC_bd[f2cs], OA2, OA2);
    clk.pause_timer("f2c");
    spdlog::debug("advance_step_OA2: l={} updated F2C", l);
  }
  // --------------------------------------------------------

  spdlog::debug(
      "NG-MPI advance_step_OA2, level={}, returning. t={}, step={}, next dt={} next time={}",
      l, SimPM.levels[l].simtime, SimPM.levels[l].step, SimPM.levels[l].dt,
      SimPM.levels[l].simtime + SimPM.levels[l].dt);
  return dt2_this + SimPM.levels[l].dt;
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::initial_conserved_quantities(
    vector<class GridBaseClass *> &grid)
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
    class cell *c = grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        // cout <<"*** LEVEL "<<l<<", cell is a leaf: "<<c->isdomain<<"
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

  initERG = SimPM.levels[0].sub_domain.global_operation_double("SUM", initERG);
  initMMX = SimPM.levels[0].sub_domain.global_operation_double("SUM", initMMX);
  initMMY = SimPM.levels[0].sub_domain.global_operation_double("SUM", initMMY);
  initMMZ = SimPM.levels[0].sub_domain.global_operation_double("SUM", initMMZ);
  initMASS =
      SimPM.levels[0].sub_domain.global_operation_double("SUM", initMASS);

  spdlog::debug(
      "(conserved quantities) [{}, {}, {}, {}, {}]", initERG, initMMX, initMMY,
      initMMZ, initMASS);

#endif  // TEST_CONSERVATION
  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_NG_MPI::check_energy_cons(vector<class GridBaseClass *> &grid)
{
  // Energy, and Linear Momentum in x-direction.
  std::vector<pion_flt> u(SimPM.nvar);
  double nowERG  = 0.;
  double nowMMX  = 0.;
  double nowMMY  = 0.;
  double nowMMZ  = 0.;
  double nowMASS = 0.0;
  double totmom  = 0.0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    double dx     = SimPM.levels[l].dx;
    double dv     = 0.0;
    class cell *c = grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        dv = spatial_solver->CellVolume(c, dx);
        spatial_solver->PtoU(c->P, u.data(), SimPM.gamma);
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

  nowERG  = SimPM.levels[0].sub_domain.global_operation_double("SUM", nowERG);
  nowMMX  = SimPM.levels[0].sub_domain.global_operation_double("SUM", nowMMX);
  nowMMY  = SimPM.levels[0].sub_domain.global_operation_double("SUM", nowMMY);
  nowMMZ  = SimPM.levels[0].sub_domain.global_operation_double("SUM", nowMMZ);
  nowMASS = SimPM.levels[0].sub_domain.global_operation_double("SUM", nowMASS);
  totmom  = SimPM.levels[0].sub_domain.global_operation_double("SUM", totmom);
  // cout <<" totmom="<<totmom<<" initMMX="<<initMMX;
  // cout <<", nowMMX="<<nowMMX<<"\n";

  // spdlog::debug(
  //    "(conserved quantities) [{}, {}, {}, {}, {}]\n"
  //    "(relative error      ) [{}, {}, {}, {}, {}]" nowERG,
  //    nowMMX, nowMMY, nowMMZ, nowMASS, (nowERG - initERG) / (initERG),
  //    (nowMMX - initMMX) / (totmom), (nowMMY - initMMY) / (totmom),
  //    (nowMMZ - initMMZ) / (totmom), (nowMASS - initMASS) / initMASS);

  return (0);
}



// ##################################################################
// ##################################################################



int sim_control_NG_MPI::RT_all_sources_levels(
    class SimParams &par  ///< simulation parameters
)
{
  /// Do this in 3 passes: 1st we go from coarse to fine, tracing the
  /// off-grid sources.  This gets 1/2 of those rays right.
  /// Then go from fine to coarse, tracing all sources and
  /// updating column densities as we go.
  /// Finally go from coarse to fine again, updating boundary data.
  int err                   = 0;
  class GridBaseClass *grid = 0;

  // --------------------------------------------------------------
  // Update off-grid sources and external boundaries.
  // for (int l=0; l<par.grid_nlevels; l++) {
#ifdef TEST_INT
  // cout <<"updating external boundaries for level "<<l<<"\n";
#endif
  // grid = par.levels[l].grid;
  // err = TimeUpdateExternalBCs(par, l, grid,
  //          spatial_solver, par.simtime,par.tmOOA,par.tmOOA);
  // rep.errorTest("NG RT_all_sources_levels: pass 1 BC-ext",0,err);
  // err = do_offgrid_raytracing(par,grid,l);
  // rep.errorTest("NG RT_all_sources_levels: pass 1 RT",0,err);
  //}
  // --------------------------------------------------------------

  // --------------------------------------------------------------
  // update internal boundaries and then all sources (fine first)
  for (int l = par.grid_nlevels - 1; l >= 0; l--) {
#ifdef TEST_INT
    cout << "RT: Receiving data for level " << l << "\n";
#endif
    grid = par.levels[l].grid;
    // receive column densities from finer grid
    if (l < par.grid_nlevels - 1) {
      for (size_t i = 0; i < grid->BC_bd.size(); i++) {
        if (grid->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += BC_update_FINE_TO_COARSE_RECV(
              par, spatial_solver, l, grid->BC_bd[i], 2, 2);
        }
      }
#ifndef NDEBUG
      spdlog::debug("CLEAR F2C send from {} to {}", l + 1, l);
#endif
      BC_FINE_TO_COARSE_SEND_clear_sends(par.levels[l + 1].sub_domain);
    }

#ifdef TEST_INT
    cout << "Doing raytracing for level " << l << "\n";
#endif
    err = do_ongrid_raytracing(par, grid, l);
    // err = do_offgrid_raytracing(par,grid,l);
    if (err)
      spdlog::error(
          "{}: Expected {} but got {}", "NG RT_all_sources_levels: pass 2 RT",
          0, err);

    // send column densities to coarser grid
    if (l > 0) {
      for (size_t i = 0; i < grid->BC_bd.size(); i++) {
        if (grid->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += BC_update_FINE_TO_COARSE_SEND(
              par, spatial_solver, l, grid->BC_bd[i], 2, 2);
        }
      }
    }

#ifdef TEST_INT
    cout << "moving on to next level.\n";
#endif
  }
  if (err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "sim_control_NG: internal boundary update", 0, err);
  // --------------------------------------------------------------

  // --------------------------------------------------------------
  // Update external boundaries.
  for (int l = 0; l < par.grid_nlevels; l++) {
#ifdef TEST_INT
    cout << "Pass 3: external boundaries for level " << l << "\n";
#endif
    grid = par.levels[l].grid;
    // C2F gets data from parent grid onto this grid.
    err = TimeUpdateExternalBCs(
        par, l, grid, spatial_solver, par.simtime, par.tmOOA, par.tmOOA);
    if (err)
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG RT_all_sources_levels: pass 3 BC-ext", 0, err);
  }
  // --------------------------------------------------------------

  return err;
}

// ##################################################################
// ##################################################################

#endif  // PARALLEL
