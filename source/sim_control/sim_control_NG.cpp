/// \file sim_control_NG.cpp
///
/// \brief Simulation Control Class for Nested Grids.
///
/// \author Jonathan Mackey
///
/// This file contains the definitions of the member functions for
/// the NG-grid simulation control class.  This is built on top
/// of the control class for uniform grids, and so doesn't add too
/// much, just the moving up and down between levels.
///
/// Modifications:
/// - 2018.05.03 JM: Started on NG grid simulation control.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "sim_control/sim_control_NG.h"

using namespace std;

//#define TEST_BC89FLUX

// ##################################################################
// ##################################################################

sim_control_NG::sim_control_NG()
{
#ifndef NDEBUG
  spdlog::debug("(sim_control_NG::Constructor)");
#endif
}

// ##################################################################
// ##################################################################

sim_control_NG::~sim_control_NG()
{
#ifndef NDEBUG
  spdlog::debug("(sim_control_NG::Destructor)");
#endif
}

// ##################################################################
// ##################################################################

int sim_control_NG::Init(
    string infile,
    int typeOfFile,
    int narg,
    string *args,
    vector<class GridBaseClass *>
        &grid  ///< address of vector of grid pointers.
)
{
#ifndef NDEBUG
  spdlog::info("(sim_control_NG::Init) Initialising grid");
#endif
  int err = 0;

  SimPM.typeofip = typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  err = dataio->ReadHeader(infile, SimPM);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::get_parameters) error", 0,
        err);
    exit(1);
  }

  // Check for commandline args that override the file parameters.
  err = override_params(narg, args);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::override_params) error", 0,
        err);
    exit(1);
  }

  //
  // Set up the Xmin/Xmax/Range/dx of each level in the NG grid
  //
  setup_NG_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
  err      = setup_grid(grid, SimPM);
  SimPM.dx = grid[0]->DX();
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::setup_grid) Something went wrong",
        0, err);
    exit(1);
  }

  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  // ----------------------------------------------------------------
  err = set_equations(SimPM);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::set_equations)", 0, err);
    exit(1);
  }
  spatial_solver->SetEOS(SimPM.gamma);

  // ----------------------------------------------------------------
  err = setup_microphysics(SimPM);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::setup_microphysics)", 0, err);
    exit(1);
  }

  // assign data to the grid from snapshot file.
  // ----------------------------------------------------------------
  err = dataio->ReadData(infile, grid, SimPM);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::assign_initial_data)", 0, err);
    exit(1);
  }

  // For each grid in the NG grid, set Ph[] = P[],
  // and then implement the boundary conditions on the grid and
  // ghost cells.
  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // Set Ph=P in every cell.
    cell *c = grid[l]->FirstPt();
    do {
      for (int v = 0; v < SimPM.nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = grid[l]->NextPt(*c)) != 0);

    if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
#ifndef NDEBUG
      spdlog::debug("Initial state, zero-ing glm variable");
#endif
      c = grid[l]->FirstPt();
      do {
        c->P[SI] = c->Ph[SI] = 0.;  // grid->divB(c);
      } while ((c = grid[l]->NextPt(*c)) != 0);
    }
  }  // loop over levels

  // Assign boundary conditions to boundary points.
  // ----------------------------------------------------------------
  err = boundary_conditions(SimPM, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(NG_INIT::boundary_conditions) err!=0",
        0, err);
    exit(1);
  }

  // Setup Raytracing on each grid, if needed.
  // ----------------------------------------------------------------
  err += setup_raytracing(SimPM, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "Failed to setup raytracer", 0, err);
    exit(1);
  }
  //  cout <<"Setting up RT sources\n";

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err = assign_boundary_data(SimPM, l, grid[l], MP);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "NG_INIT::assign_boundary_data", 0,
          err);
      exit(1);
    }
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifndef NDEBUG
    spdlog::debug("updating external boundaries for level {}", l);
#endif
    err += TimeUpdateExternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "NG_INIT: error from bounday update", 0,
        err);
    exit(1);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifndef NDEBUG
    spdlog::debug("updating internal boundaries for level {}", l);
#endif
    err += TimeUpdateInternalBCs(
        SimPM, l, grid[l], spatial_solver, SimPM.simtime, 0.0, SimPM.tmOOA,
        SimPM.tmOOA);
  }
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "NG_INIT: error from bounday update", 0,
        err);
    exit(1);
  }
  // ----------------------------------------------------------------

  //
  // If testing the code, this calculates the momentum and energy on
  // the domain.
  //
  initial_conserved_quantities(grid);

  //
  // If using opfreq_time, set the next output time correctly.
  //
  if (SimPM.op_criterion == 1) {
    if (SimPM.opfreq_time < TINYVALUE) {
      spdlog::error(
          "{}: {}", "opfreq_time not set right and is needed!",
          SimPM.opfreq_time);
      exit(1);
    }
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
    if (!dataio) {
      spdlog::error(
          "{}: {}", "NG_INIT:: dataio initialisation", SimPM.typeofop);
      exit(1);
    }
  }
  dataio->SetSolver(spatial_solver);
  dataio->SetMicrophysics(MP);
  if (textio) {
    textio->SetSolver(spatial_solver);
    textio->SetMicrophysics(MP);
  }

#ifndef NDEBUG
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    cell *c = (grid[l])->FirstPt_All();
    do {
      if (pconst.equalD(c->P[RO], 0.0)) {
        cout << "zero data in cell: ";
        CI.print_cell(*c);
      }
    } while ((c = (grid[l])->NextPt_All(*c)) != 0);
  }
#endif  // NDEBUG

  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_NG::initial_conserved_quantities(
    vector<class GridBaseClass *> &grid)
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION
  pion_flt u[SimPM.nvar];
  initERG = 0.;
  initMMX = initMMY = initMMZ = 0.;
  initMASS                    = 0.0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    double dx     = SimPM.levels[l].dx;
    double dv     = 0.0;
    class cell *c = grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        dv = spatial_solver->CellVolume(c, dx);
        spatial_solver->PtoU(c->P, u, SimPM.gamma);
        initERG += u[ERG] * dv;
        initMMX += u[MMX] * dv;
        initMMY += u[MMY] * dv;
        initMMZ += u[MMZ] * dv;
        initMASS += u[RHO] * dv;
      }
    } while ((c = grid[l]->NextPt(*c)) != 0);
  }

  spdlog::debug(
      "(conserved quantities) [{}, {}, {}, {}, {}]\n", initERG, initMMX,
      initMMY, initMMZ, initMASS);

#endif  // TEST_CONSERVATION
  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_NG::Time_Int(
    vector<class GridBaseClass *> &grid  ///< vector of grids
)
{
  spdlog::info(
      "-------------------------------------------------------\n(sim_control_NG::Time_Int) STARTING TIME INTEGRATION\n-------------------------------------------------------\n");
  int err         = 0;
  SimPM.maxtime   = false;
  bool first_step = true;
  bool restart    = true;
  if (SimPM.timestep != 0) first_step = false;
  if (SimPM.timestep == 0) restart = false;
  clk.start_timer("time_int");
  double tsf = 0;

  // make sure all levels start at the same time.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt      = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }

  // cout <<"raytracing all levels...\n";
  // Do raytracing on all levels, and update internal and external
  // boundaries to populate the column densities correctly.
  // Even if there is not RT, this updates the boundaries.
  err = RT_all_sources_levels(SimPM);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "sim_control_NG: RT_all_sources_levels",
        0, err);
    exit(1);
  }
  // cout <<"raytracing all levels... finished\n";
  if (SimPM.timestep == 0) {
    // cout << "(step=0) Writing initial data.\n";
    err = output_data(grid);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "Failed to write file... path?", 0,
          err);
      exit(1);
    }
  }

  while (SimPM.maxtime == false) {

#if defined(CHECK_MAGP)
    calculate_magnetic_pressure(grid);
#elif defined(BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid);
#endif

    // Get timestep on each level
    int scale    = 1;
    double mindt = 1.0e99;
    // err = RT_all_sources_levels(SimPM);
    // rep.errorTest("sim_control_NG: RT_all_sources_levels",0,err);

    for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TEST_INT
      spdlog::debug(
          "Calculate timestep, level {}, dx={}", l, SimPM.levels[l].dx);
#endif
      if (!restart && !first_step) {
        SimPM.last_dt = SimPM.levels[l].last_dt;
      }
      else {
        SimPM.levels[l].last_dt = SimPM.last_dt / SimPM.levels[l].multiplier;
      }

      err += calculate_timestep(SimPM, grid[l], l);
      if (0 != err) {
        spdlog::error(
            "{}: Expected {} but got {}", "TIME_INT::calc_timestep()", 0, err);
        exit(1);
      }

      mindt = std::min(mindt, SimPM.dt / scale);
#ifdef TEST_INT
      spdlog::debug("level {} got dt={} and {}", l, SimPM.dt, SimPM.dt / scale);
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
      spdlog::debug(
          "new dt={}, t={}", SimPM.levels[l].dt, SimPM.levels[l].simtime);
#endif
    }
    if (first_step) {
      // take a ~3x smaller timestep for the first timestep.
      for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
        // cout <<"level "<<l<<", orig dt="<<SimPM.levels[l].dt;
        SimPM.levels[l].dt *= 0.3;
      }
      first_step = false;
    }
    if (restart) restart = false;
    SimPM.last_dt = SimPM.levels[0].last_dt;

    //
    // Use a recursive algorithm to update the coarsest level.
    //
    advance_time(0, SimPM.levels[0].grid);
    SimPM.simtime = SimPM.levels[0].simtime;

#if !defined(CHECK_MAGP)
#if !defined(BLAST_WAVE_CHECK)
    tsf = clk.time_so_far("time_int");
    spdlog::info(
        "New time: {:12.6e}   dt: {:12.6e}   steps: {:8d}   l0 steps: {:6d}   runtime: {:12.4e} s",
        SimPM.simtime, SimPM.levels[SimPM.grid_nlevels - 1].dt, SimPM.timestep,
        SimPM.timestep / static_cast<int>(pow(2, SimPM.grid_nlevels - 1)), tsf);
#endif
#endif

#ifdef TEST_CONSERVATION
    err += check_energy_cons(grid);
#endif /* TEST_CONSERVATION */

    err += output_data(grid);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "TIME_INT::output_data()", 0, err);
      exit(1);
    }
    err += check_eosim();
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "TIME_INT::check_eosim()", 0, err);
      exit(1);
    }
  }

  spdlog::debug(
      "(sim_control_NG::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n");

  tsf = clk.time_so_far("time_int");
  spdlog::info(
      "TOTALS: Nsteps: {} wall-time: {} time/step: {}\nSTEPS: {}\t{}\t{}\t{}\n-------------------------------------------------------\n",
      SimPM.timestep, tsf, tsf / static_cast<double>(SimPM.timestep),
      SimPM.timestep, tsf, tsf / static_cast<double>(SimPM.timestep),
      static_cast<double>(SimPM.timestep * SimPM.Ncell) / tsf);

  // cout.setf(ios_base::scientific);
  // cout.precision(6);

  return (0);
}

// ##################################################################
// ##################################################################

#ifdef CHECK_MAGP
///
/// This is only for a test problem -- it checks the magnetic
/// pressure on the full domain and outputs it to screen
///
void sim_control_NG::calculate_magnetic_pressure(
    vector<class GridBaseClass *> &grid  ///< grid pointers.
)
{
  //
  // Calculate the total magnetic pressure on the domain, normalised to the
  // initial value.
  //
  double magp = 0.0, cellvol = 0.0;
  static double init_magp = -1.0;

  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {

    cell *c = grid[l]->FirstPt();
    do {
      if (!c->isbd && c->isleaf)
        magp += (spatial_solver->Ptot(c->P, 0.0) - c->P[PG])
                * spatial_solver->CellVolume(c, grid[l]->DX());
    } while ((c = grid[l]->NextPt(*c)) != 0);
  }

  if (init_magp < 0) init_magp = magp;
  spdlog::debug("{}\t{}\t{}", SimPM.simtime, magp / init_magp, magp);
}
#endif  // CHECK_MAGP

// ##################################################################
// ##################################################################

#ifdef BLAST_WAVE_CHECK
///
/// If running a 1D spherical blast wave, calculate the shock position
/// and output to screen.
///
void sim_control_NG::calculate_blastwave_radius(
    vector<class GridBaseClass *> &grid  ///< grid pointers.
)
{
  //
  // Calculate the blast wave outer shock position.
  // If a NG grid, start on the finest grid and work outwards
  //
  double shockpos       = 0.0;
  static double old_pos = 0.0;
  bool shock_found      = false;
  //  static double last_dt=0.0;
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l++) {

    if (shock_found) continue;
    cell *c = grid->LastPt();
    if (fabs(c->P[VX]) >= 1.0e4) {
      spdlog::debug("level {} does not contain shock.", l);
    }
    else {
      do {
        c = grid->NextPt(*c, RNsph);
        // cout <<c->id<<", vx="<<c->P[VX]<<"\n";
      } while (c != 0 && fabs(c->P[VX]) < 1.0e4);
      if (c && (c->P[VX] >= 1.0e4)) {
        shockpos    = CI.get_dpos(c, Rsph);
        shock_found = true;
      }
    }
  }
  if (pconst.equalD(old_pos, 0.0)) old_pos = shockpos;
  spdlog::debug("{}\t{}", SimPM.simtime, shockpos);
  // cout <<"\t"<<(shockpos-old_pos)/(SimPM.dt+TINYVALUE);
  old_pos = shockpos;
  return;
}
#endif  // BLAST_WAVE_CHECK

// ##################################################################
// ##################################################################

int sim_control_NG::Finalise(vector<class GridBaseClass *>
                                 &grid  ///< address of vector of grid pointers.
)
{
  int err = 0;
  spdlog::info(
      "------------------------------------------------------------\n(sim_control::Finalise) FINALISING SIMULATION.");
#ifdef TEST_CONSERVATION
  err += check_energy_cons(grid);
#endif /* TEST_CONSERVATION */
  err += output_data(grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "(FINALISE::output_data) Something went wrong", 0, err);
    exit(1);
  }
  spdlog::info(
      "\tSimTime = {}   #timesteps = {}", SimPM.simtime, SimPM.timestep);
#ifndef NDEBUG
  spdlog::info("(sim_control::Finalise) DONE.");
#endif
  spdlog::info("------------------------------------------------------------");
  return (0);
}

// ##################################################################
// ##################################################################

double sim_control_NG::advance_time(
    const int l,               ///< level to advance.
    class GridBaseClass *grid  ///< grid pointer
)
{
#ifndef NDEBUG
  spdlog::debug("advance_time, level={}, starting.\n", l);
#endif

  double step = 0.0;
  if (SimPM.tmOOA == 1) {
    step = advance_step_OA1(l);
  }
  else if (SimPM.tmOOA == 2) {
#ifdef TEST_INT
    // cout <<"Calling advance_step_OA2: level "<<l<<"\n";
#endif
    step = advance_step_OA2(l);
  }
  return step;
}



// ##################################################################
// ##################################################################



double sim_control_NG::advance_step_OA1(const int l  ///< level to advance.
)
{
#ifndef NDEBUG
  spdlog::debug("advance_step_OA1, level={}, starting.", l);
#endif
  int err = 0;
  // double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this           = 0.0;  // two timesteps for this level.
  class GridBaseClass *grid = SimPM.levels[l].grid;

  err = update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error", 0,
        err);
    exit(1);
  }

  // --------------------------------------------------------
  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime, OA1, OA1);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // take the first finer grid step, if there is a finer grid.
  if (l < SimPM.grid_nlevels - 1) {
    // dt2_fine = advance_step_OA1(l+1);
    advance_step_OA1(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // now calculate dU, change in conserved variables on this grid,
  // for this step.
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(SimPM.levels[l].dt);
#ifdef PION_OMP
  }
#endif
  err += calc_microphysics_dU(SimPM.levels[l].dt, grid);
  err += calc_dynamics_dU(SimPM.levels[l].dt, OA1, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA1: calc_x_dU", 0,
        err);
    exit(1);
  }
  if (l > 0) save_fine_fluxes(SimPM, l);
  if (l < SimPM.grid_nlevels - 1) save_coarse_fluxes(SimPM, l);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // take the second finer grid step, if there is a finer grid.
  if (l < SimPM.grid_nlevels - 1) {
    // dt2_fine = advance_step_OA1(l+1);
    advance_step_OA1(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  // First correct fluxes
  //
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(SimPM.levels[l].dt);
#ifdef PION_OMP
  }
#endif
#ifndef SKIP_BC89_FLUX
  if (l < SimPM.grid_nlevels - 1) {
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA1, OA1);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "scn::advance_step_OA1: recv_BC89_flux",
          0, err);
      exit(1);
    }
  }
#endif
  err += grid_update_state_vector(SimPM.levels[l].dt, OA1, OA1, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA1: state-vec update",
        0, err);
    exit(1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step++;
  if (l == SimPM.grid_nlevels - 1) {
    SimPM.timestep++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l == 0) SimPM.last_dt = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // update internal boundaries.
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime,
      SimPM.levels[l].dt, OA1, OA1);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  if (grid->RT) {
    err += do_ongrid_raytracing(SimPM, grid, l);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG-MPI::advance_step_OA1: calc_rt_cols()", 0, err);
      exit(1);
    }
  }
  // --------------------------------------------------------

#ifndef NDEBUG
  spdlog::debug(
      "advance_step_OA1, level={}, returning. t={}, step={}, next dt={} next time={}",
      l, SimPM.levels[l].simtime, SimPM.levels[l].step, SimPM.levels[l].dt,
      SimPM.levels[l].simtime + SimPM.levels[l].dt);
#endif
  return dt2_this + SimPM.levels[l].dt;
}



// ##################################################################
// ##################################################################



double sim_control_NG::advance_step_OA2(const int l  ///< level to advance.
)
{
#ifdef TEST_INT
  spdlog::debug(
      "advance_step_OA2, level={}, starting. {}, step={}", l,
      SimPM.levels[l].simtime, SimPM.levels[l].step);
#endif
  int err = 0;
  // double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this           = 0.0;  // two timesteps for this level.
  double ctime              = SimPM.levels[l].simtime;  // current time
  class GridBaseClass *grid = SimPM.levels[l].grid;

  err = update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid->RT);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "NG TIME_INT::update_RT_sources error", 0,
        err);
    exit(1);
  }

  // --------------------------------------------------------
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver, ctime, OA2, OA2);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // take the first finer grid step, if there is a finer grid.
  if (l < SimPM.grid_nlevels - 1) {
    // dt2_fine = advance_step_OA2(l+1);
    advance_step_OA2(l + 1);
  }
  dt2_this = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Predictor step: use 0.5*dt to get to time-centred state
  double dt_now = dt2_this * 0.5;  // half of the timestep
  spatial_solver->Setdt(dt_now);
  err += calc_microphysics_dU(dt_now, grid);
  err += calc_dynamics_dU(dt_now, OA1, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: calc_x_dU OA1", 0,
        err);
    exit(1);
  }

  err += grid_update_state_vector(dt_now, OA1, OA2, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: update OA1", 0,
        err);
    exit(1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Update boundary data.
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, dt_now, OA1, OA2);
  err += TimeUpdateExternalBCs(
      SimPM, l, grid, spatial_solver, ctime + dt_now, OA1, OA2);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "scn::advance_step_OA2: bounday update OA1", 0, err);
    exit(1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Now calculate dU for the full step (OA2)
  //
  dt_now = dt2_this;  // full step
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt_now);
#ifdef PION_OMP
  }
#endif
  err += do_ongrid_raytracing(SimPM, grid, l);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_time: calc_rt_cols() OA2",
        0, err);
    exit(1);
  }
  err += calc_microphysics_dU(dt_now, grid);
  err += calc_dynamics_dU(dt_now, OA2, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: calc_x_dU OA2", 0,
        err);
    exit(1);
  }
  // save fluxes at level boundaries
  if (l > 0) save_fine_fluxes(SimPM, l);
  if (l < SimPM.grid_nlevels - 1) save_coarse_fluxes(SimPM, l);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // take the second finer grid step, if there is a finer grid.
  if (l < SimPM.grid_nlevels - 1) {
    // dt2_fine = advance_step_OA2(l+1);
    advance_step_OA2(l + 1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Now update Ph[i] to new values (and P[i] also if full step).
  // First correct fluxes
  //
#ifdef TEST_INT
  spdlog::debug("l={} full step, grid_update_state_vector", l);
#endif
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
    err += recv_BC89_fluxes_F2C(spatial_solver, SimPM, l, OA2, OA2);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "scn::advance_step_OA1: recv_BC89_flux",
          0, err);
      exit(1);
    }
  }
#endif
  err += grid_update_state_vector(dt_now, OA2, OA2, grid);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "scn::advance_step_OA2: update OA2", 0,
        err);
    exit(1);
  }
  // --------------------------------------------------------

  // --------------------------------------------------------
  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step++;
  if (l == SimPM.grid_nlevels - 1) {
    SimPM.timestep++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l == 0) SimPM.last_dt = SimPM.levels[l].dt;
  // --------------------------------------------------------

  // --------------------------------------------------------
  // update internal boundaries.
  //
  err += TimeUpdateInternalBCs(
      SimPM, l, grid, spatial_solver, SimPM.levels[l].simtime, 0.5 * dt_now,
      OA2, OA2);
  // --------------------------------------------------------

  // --------------------------------------------------------
  // Do raytracing for next step, to send with F2C BCs.
  // --------------------------------------------------------
  err += do_ongrid_raytracing(SimPM, grid, l);
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}", "NG-MPI::advance_step_OA2: raytracing()",
        0, err);
    exit(1);
  }
  // --------------------------------------------------------

#ifdef TEST_INT
  spdlog::debug(
      "advance_step_OA2, level={}, returning. t={}, step={}", l,
      SimPM.levels[l].simtime, SimPM.levels[l].step);
#endif
  return dt2_this + SimPM.levels[l].dt;
}

// ##################################################################
// ##################################################################

int sim_control_NG::check_energy_cons(vector<class GridBaseClass *> &grid)
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
      if (!c->isbd && c->isgd) {
        dv = spatial_solver->CellVolume(*c, dx);
        spatial_solver->PtoU(c->P.data(), u.data(), SimPM.gamma);
        nowERG += u[ERG] * dv;
        nowMMX += u[MMX] * dv;
        nowMMY += u[MMY] * dv;
        nowMMZ += u[MMZ] * dv;
        nowMASS += u[RHO] * dv;
        totmom +=
            sqrt(u[MMX] * u[MMX] + u[MMY] * u[MMY] + u[MMZ] * u[MMZ]) * dv;
      }
    } while ((c = grid[l]->NextPt(*c)) != 0);
  }

  // spdlog::debug(
  //    "(conserved quantities) [{}, {}, {}, {}, {}]\n(relative error      )
  //    [{}, {}, {}, {}, {}]", nowERG, nowMMX, nowMMY, nowMMZ, nowMASS, (nowERG
  //    - initERG) / (initERG), (nowMMX - initMMX) / (totmom), (nowMMY -
  //    initMMY) / (totmom), (nowMMZ - initMMZ) / (totmom), (nowMASS - initMASS)
  //    / initMASS);

  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_NG::do_ongrid_raytracing(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< Computational grid.
    const int l                 ///< level in NG
)
{
  if (!grid->RT) return 0;
  int err = 0;
  //
  // If we have raytracing, we call the ray-tracing routines
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
    if (!SimPM.RS.sources[isrc].ongrid) continue;
#ifdef RT_TESTING
    spdlog::debug("calc_raytracing_col_dens: SRC-ID: {}", isrc);
#endif
    err += grid->RT->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      spdlog::debug("isrc={}\t", isrc);
      spdlog::error("{}: {}", "ongrid RT: RT return", err);
      exit(1);
    }  // if error
  }    // loop over sources
  return err;
}

// ##################################################################
// ##################################################################

int sim_control_NG::do_offgrid_raytracing(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< Computational grid.
    const int)
{
  if (!grid->RT) return 0;
  int err = 0;
  //
  // If we have raytracing, we call the ray-tracing routines
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].ongrid) {
      // cout <<"skipping source "<<isrc<<" b/c ongrid\n";
      continue;
    }
#ifdef RT_TESTING
    spdlog::debug("calc_raytracing_col_dens: SRC-ID: {}", isrc);
#endif
    err += grid->RT->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      spdlog::debug("isrc={}\t", isrc);
      spdlog::error("{}: {}", "offgrid RT: RT return", err);
      exit(1);
    }  // if error
  }    // loop over sources
  return err;
}

// ##################################################################
// ##################################################################

int sim_control_NG::RT_all_sources_levels(
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
    spdlog::debug("updating internal boundaries for level {}", l);
#endif
    grid = par.levels[l].grid;
    // F2C gets data from child grid onto this grid.
    err = TimeUpdateInternalBCs(
        par, l, grid, spatial_solver, par.simtime, 0.0, par.tmOOA, par.tmOOA);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG RT_all_sources_levels: pass 2 BC-int", 0, err);
      exit(1);
    }
#ifdef TEST_INT
    spdlog::debug("doing raytracing for level {}", l);
#endif
    err = do_ongrid_raytracing(par, grid, l);
    // err = do_offgrid_raytracing(par,grid,l);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "NG RT_all_sources_levels: pass 2 RT",
          0, err);
      exit(1);
    }
#ifdef TEST_INT
    spdlog::debug("moving on to next level.");
#endif
  }
  if (0 != err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "sim_control_NG: internal boundary update", 0, err);
    exit(1);
  }
  // --------------------------------------------------------------

  // --------------------------------------------------------------
  // Update external boundaries.
  for (int l = 0; l < par.grid_nlevels; l++) {
#ifdef TEST_INT
    spdlog::debug("Pass 3: external boundaries for level {}", l);
#endif
    grid = par.levels[l].grid;
    // C2F gets data from parent grid onto this grid.
    err = TimeUpdateExternalBCs(
        par, l, grid, spatial_solver, par.simtime, par.tmOOA, par.tmOOA);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}",
          "NG RT_all_sources_levels: pass 3 BC-ext", 0, err);
      exit(1);
    }
  }
  // --------------------------------------------------------------

  return err;
}

// ##################################################################
// ##################################################################
