/// \file sim_control_MPI.cpp
///
/// \brief Parallel Grid Methods Class Member Function definitions.
///
/// \author Jonathan Mackey
///
/// This file contains the definitions of the member functions for
/// the "sim_control_pllel" class, which is a modification
/// of the basic 1st/2nd order Finite Volume Solver according to the
/// method outlined in Falle, Komissarov, \& Joarder (1998),MNRAS,297,265.
///
/// Modifications:
/// - 2007-10-11 Started writing file.
/// - 2007-11-11 Basically working with uniform grid and all sorts of
/// boundaries.
/// - 2010-01-05 JM: Put RT step limiting in an #ifdef so that it is
///    not set for test problems.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits.  This also
///    meant I don't need a new parallel output_data() function so I
///    got rid of it.  All dataio classes now choose the filename
///   themselves.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///    extra_data.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.15 JM: Added CI.setup_extra_data() call in setup_grid()
///     function.
/// - 2011.01.12 JM: Set so only proc 0 displays RT timestep.
///     Commented out some stuff in calc_timestep which was used for
///     testing and made the code less efficient.
/// - 2011.02.17 JM: Raytracer header file moved to raytracing/ subdir.
///     More ray-tracing options for CI.setup_extra_data
/// - 2011.02.25 JM: removed HCORR ifdef around new code. Modified
/// setup_extra_data
///     for ray-tracing again.  Now we have N column-density vars for N sources.
/// - 2011.03.01 JM: fixed timings bug in time_int()
/// - 2011.03.02 JM: Added parallelised raytracer-shielding class setup.
/// - 2011.03.21 JM: moved cell setup_extra_data() to its own function, to save
///     on copy-paste for the parallel version.
///     Deleted commented out output_data() function, which is now the same as
///     the serial version.
/// - 2011.03.22 JM: Updated setup_raytracing() for improved raytracer
/// functionality.
/// - 2011.10.24 JM: Updated setup_raytracing() again.  It's better now.
/// - 2012.01.16 JM: Added update_evolving_RT_source() function to time_int().
/// - 2012.05.14 JM: Added check in setup_raytracing() for sources at infinity
///    in different directions (must set up the parallelised raytracer then).
/// - 2012.05.16 JM: fixed bug from last change.
/// - 2012.08.06 JM: Added separators between functions for clarity.
/// - 2013.04.16 JM: Fixed FITS read functions for new filename convention.
/// - 2013.09.05 JM: changed RS position[] to pos[].
/// - 2013.10.13 JM: Tidied up a bit.
/// - 2015.[01.26-02.03] JM: CHANGED FILENAME TO SIM_CONTROL_MPI.CPP,
///    added ParallelParams class, and fixing code for non-global mpiPM.
/// - 2015.02.18 JM: moved setup functions to setup_fixed_grid_MPI
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).  Changed default I/O to DOUBLE precision.
/// - 2017.08.03 JM: Changed silo dataio class to be the utility class.
/// - 2021.07.13 JM: updated Time_Int() to remove extra boundary updates.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "raytracing/raytracer_SC.h"
#include "sim_control/sim_control_MPI.h"
#include "sub_domain/sub_domain.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"
#endif  // if FITS

#include <fstream>

#include <sstream>
using namespace std;

#ifdef PARALLEL

// ##################################################################
// ##################################################################

sim_control_pllel::sim_control_pllel() : sim_control()
{
#ifndef NDEBUG
  spdlog::debug("sim_control_pllel constructor");
#endif
}

// ##################################################################
// ##################################################################

sim_control_pllel::~sim_control_pllel()
{
#ifndef NDEBUG
  spdlog::debug("sim_control_pllel destructor");
#endif
}

// ##################################################################
// ##################################################################

int sim_control_pllel::Init(
    string infile,
    int typeOfFile,
    int narg,
    string *args,
    vector<class GridBaseClass *> &grid  ///< vector of grids.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "(sim_control_pllel::init) Initialising grid: infile = {}", infile);
#endif
  int err = 0;

  //
  // Setup dataI/O class and check if we read from a single file or
  // multiple files.  Should be a single file, but for FITS it might
  // not be.
  //
  setup_dataio_class(SimPM, typeOfFile);
  if (!dataio->file_exists(infile))
    spdlog::error("{}: {}", "infile doesn't exist!", infile);

  if (typeOfFile == 2) {
    // FITS
    string::size_type pos = infile.find("_0000.");
    if (pos == string::npos) {
      SimPM.levels[0].sub_domain.set_ReadSingleFile(true);
    }
    else {
      SimPM.levels[0].sub_domain.set_ReadSingleFile(false);
      ostringstream t;
      t.str("");
      t << "_";
      t.width(4);
      t.fill('0');
      t << SimPM.levels[0].sub_domain.get_myrank() << ".";
      string t2 = t.str();
      infile.replace(pos, 6, t2);
    }
  }
  else if (typeOfFile == 5) {
    // SILO
    string::size_type pos = infile.find("_0000.");
    if (pos == string::npos) {
      SimPM.levels[0].sub_domain.set_ReadSingleFile(true);
    }
    else {
      SimPM.levels[0].sub_domain.set_ReadSingleFile(false);
    }
  }
  else
    spdlog::error("{}: {}", "bad input file type in Init", typeOfFile);

  //
  // We need to decompose the domain here, because setup_grid() needs
  // this, but this means we need to read the header to find out what
  // the grid dimensions are.  So the header is read twice, but this
  // should be ok because it only happens during initialisation.
  //
  err = dataio->ReadHeader(infile, SimPM);
  if (0 != err) {
    spdlog::error("PLLEL Init(): failed to read header: {}", err);
    exit(err);
  }

  // have to do something with SimPM.levels[0] because this
  // is used to set the local domain size in decomposeDomain
  SimPM.levels[0].parent = 0;
  SimPM.levels[0].child  = 0;
  SimPM.levels[0].Ncell  = SimPM.Ncell;
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].NG[v] = SimPM.NG[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Range[v] = SimPM.Range[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Xmin[v] = SimPM.Xmin[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Xmax[v] = SimPM.Xmax[v];
  SimPM.levels[0].dx         = SimPM.Range[XX] / SimPM.NG[XX];
  SimPM.levels[0].simtime    = SimPM.simtime;
  SimPM.levels[0].dt         = 0.0;
  SimPM.levels[0].multiplier = 1;

  std::vector<int> pbc = SimPM.get_pbc_bools();
  err                  = SimPM.levels[0].sub_domain.decomposeDomain(
      SimPM.ndim, SimPM.levels[0], std::move(pbc));
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "PLLEL Init():Couldn't Decompose Domain!",
        0, err);

  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  if (0 != err) {
    spdlog::error("(INIT::override_params) err = {}", err);
    exit(err);
  }

  // Now set up the grid structure.
  grid.resize(1);
  // cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  err = setup_grid(grid, SimPM);
  // cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  SimPM.dx             = grid[0]->DX();
  SimPM.levels[0].grid = grid[0];
  if (0 != err) {
    spdlog::error("(INIT::setup_grid) err = {}g", err);
    exit(err);
  }

  //
  // All grid parameters are now set, so I can set up the appropriate
  // equations/solver class.
  //
  err = set_equations(SimPM);
  if (0 != err) {
    spdlog::error("(INIT::set_equations) err = {} Fix me!", err);
    exit(err);
  }
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->SetEOS(SimPM.gamma);
#ifdef PION_OMP
  }
#endif

  //
  // Now setup Microphysics, if needed.
  //
  err = setup_microphysics(SimPM);
  if (0 != err) {
    spdlog::error("(INIT::setup_microphysics) err!=0 {} {}", 0, err);
    exit(1);
  }
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->SetMicrophysics(MP);
#ifdef PION_OMP
  }
#endif

  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, grid, SimPM);
  if (0 != err) {
    spdlog::error("(INIT::assign_initial_data) err = {}", err);
    exit(err);
  }
  //  cout <<"Read data finished\n";

  //
  // Set Ph[] = P[], and then implement the boundary conditions.
  //
  cell *c = grid[0]->FirstPt();
  do {
    for (int v = 0; v < SimPM.nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = grid[0]->NextPt(*c)) != 0);

  //
  // If I'm using the GLM method, make sure Psi variable is
  // initialised to zero.
  //
  if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
#ifndef NDEBUG
    spdlog::debug("Initial state, zero-ing glm variable");
#endif
    c = grid[0]->FirstPt();
    do {
      c->P[SI] = c->Ph[SI] = 0.;
    } while ((c = grid[0]->NextPt(*c)) != 0);
  }

  //
  // Assign boundary conditions to boundary points.
  //
  //  cout <<"starting BC setup\n";
  err = boundary_conditions(SimPM, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::boundary_conditions) err!=0", 0,
        err);

  //  cout <<"assigning BC data\n";
  err = assign_boundary_data(SimPM, 0, grid[0], MP);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::assign_boundary_data) err!=0", 0,
        err);

  //
  // Setup Raytracing on each grid, if needed.
  //
  //  cout <<"Setting up raytracing\n";
  err += setup_raytracing(SimPM, grid[0]);
  // cout <<"Setting up RT sources\n";
  err += setup_evolving_RT_sources(SimPM);
  // cout <<"Updating evolving RT sources\n";
  err += update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "Failed to setup raytracer and/or microphysics", 0, err);

  //
  // If testing the code, this calculates the momentum and energy on the
  // domain.
  //
  //  cout <<"initial conserved quantities\n";
  initial_conserved_quantities(grid);

  err += TimeUpdateInternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, 0.0, SimPM.tmOOA,
      SimPM.tmOOA);
  err += TimeUpdateExternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, SimPM.tmOOA,
      SimPM.tmOOA);
  if (err)
    spdlog::error(
        "{}: {}", "first_order_update: error from bounday update", err);

  //
  // If using opfreq_time, set the next output time correctly.
  //
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
      spdlog::error("{}: {}", "INIT:: dataio initialisation", SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  dataio->SetMicrophysics(MP);
  if (textio) {
    textio->SetSolver(spatial_solver);
    textio->SetMicrophysics(MP);
  }

  if (SimPM.timestep == 0) {
#ifndef NDEBUG
    spdlog::info("(PARALLEL INIT) Writing initial data");
#endif
    output_data(grid);
  }

#ifndef NDEBUG
  c = (grid[0])->FirstPt_All();
  do {
    if (pconst.equalD(c->P[RO], 0.0)) {
      cout << "zero data in cell: ";
      CI.print_cell(*c);
    }
  } while ((c = (grid[0])->NextPt_All(*c)) != 0);
#endif  // NDEBUG
  return (err);
}

// ##################################################################
// ##################################################################

/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int sim_control_pllel::Time_Int(
    vector<class GridBaseClass *> &grid  ///< vector of grids.
)
{
  spdlog::info("------------ (sim_control_pllel::time_int)"
               "STARTING TIME INTEGRATION ------------");
  int err = 0;
  err     = update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  if (err) {
    spdlog::error("TIME_INT:: initial RT src update() err = {}", err);
    exit(1);
  }
  err = RT_all_sources(SimPM, grid[0], 0);
  if (err) {
    spdlog::error("TIME_INT:: initial RT() err = {}", err);
    exit(1);
  }
  int log_freq  = 10;
  SimPM.maxtime = false;
  clk.start_timer("time_int");
  double tsf = 0.0;
  while (SimPM.maxtime == false) {
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
    if (err) {
      spdlog::error("TIME_INT::update_RT_sources() err = {}", err);
      exit(1);
    }
    err = RT_all_sources(SimPM, grid[0], 0);
    if (err) {
      spdlog::error("TIME_INT:: loop RT() err = {}", err);
      exit(1);
    }
    // clk.start_timer("advance_time");
#ifndef NDEBUG
    spdlog::debug("MPI time_int: calculating dt");
#endif
    SimPM.levels[0].last_dt = SimPM.last_dt;
    err += calculate_timestep(SimPM, grid[0], 0);
    if (err) {
      spdlog::error("TIME_INT::calc_timestep() err = {}", err);
      exit(1);
    }

#ifndef NDEBUG
    spdlog::debug("MPI time_int: stepping forward in time");
#endif
    advance_time(0, grid[0]);
    // cout <<"advance_time took "<<clk.stop_timer("advance_time")<<"
    // secs.\n";
#ifndef NDEBUG
    spdlog::debug("MPI time_int: finished timestep");
    log_freq = 1;
#endif

    if ((SimPM.levels[0].sub_domain.get_myrank() == 0)
        && (SimPM.timestep % log_freq) == 0) {
      tsf = clk.time_so_far("time_int");
      spdlog::info(
          "New time: {:12.6e}   dt: {:12.6e}   steps: {:8d}   runtime: {:12.4e} s",
          SimPM.simtime, SimPM.dt, SimPM.timestep, tsf);
#ifndef NDEBUG
#endif  // NDEBUG
    }

#ifdef TEST_CONSERVATION
    if ((SimPM.timestep % log_freq) == 0) {
      err += check_energy_cons(grid[0]);
    }
#endif

    //
    // check if we are at time limit yet.
    //
    tsf         = clk.time_so_far("time_int");
    double maxt = SimPM.levels[0].sub_domain.global_operation_double(MAX, tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime = true;
      spdlog::debug("RUNTIME>{} SECS.\n", get_max_walltime());
    }

    err += output_data(grid);
    if (err) {
      spdlog::error("(TIME_INT::output_data) err!=0 Something went bad");
      return (1);
    }

    err += check_eosim();
    if (err) {
      spdlog::error("(TIME_INT::) err!=0 Something went bad");
      return (1);
    }
  }
  spdlog::info("(sim_control_pllel::time_int) FINISHED - FINALISING SIM");
  tsf = clk.time_so_far("time_int");
  spdlog::info(
      "TOTALS: Nsteps: {}, sim-time: {}, wall-time: {}, time/step: {}",
      SimPM.timestep, SimPM.simtime, tsf,
      tsf / static_cast<double>(SimPM.timestep));
#ifdef RT_TESTING
  if (grid[0]->RT != 0) {
    //
    // output raytracing timing info.  Have to start and stop timers to get
    // the correct runtime (this is sort of a bug... there is no function to
    // get the elapsed time of a non-running timer.  I should add that.
    //
    string t1 = "totalRT", t2 = "waitingRT", t3 = "doingRT";
    double total = 0.0, wait = 0.0, run = 0.0;
    clk.start_timer(t1);
    total = clk.pause_timer(t1);
    clk.start_timer(t2);
    wait = clk.pause_timer(t2);
    clk.start_timer(t3);
    run = clk.pause_timer(t3);
    spdlog::debug("TOTALS RT: active: {} idle: {} total: {}", run, wait, total);
  }
#endif
  spdlog::info("**************************************");
  return (0);
}



// ##################################################################
// ##################################################################



int sim_control_pllel::calculate_timestep(
    class SimParams &par,       ///< simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level to advance (for NG grid)
)
{
  //
  // First get the local grid dynamics and microphysics timesteps.
  //
  double t_dyn = 0.0, t_mp = 0.0;
  t_dyn = calc_dynamics_dt(par, grid);
  t_mp  = calc_microphysics_dt(par, grid);
#ifndef NDEBUG
  spdlog::info("calc_time: local t_dyn= {}, t_mp= {}", t_dyn, t_mp);
#endif


  //
  // Get global min over all grids on this level.
  //
  t_dyn = SimPM.levels[0].sub_domain.global_operation_double(MIN, t_dyn);
  t_mp  = SimPM.levels[0].sub_domain.global_operation_double(MIN, t_mp);

#ifndef NDEBUG
  spdlog::info("calc_time: global t_dyn= {}, t_mp= {}", t_dyn, t_mp);
#endif

#ifndef NDEBUG
  // Write step-limiting info every tenth timestep.
  if (t_mp < t_dyn && (par.timestep % 10) == 0)
    spdlog::info("Limiting timestep by MP: mp_t={}\thydro_t={}", t_mp, t_dyn);
#endif

  par.dt = min(t_dyn, t_mp);

  //
  // Calculate the timestep limit imposed by thermal conduction,
  // and calcuate the multidimensional energy fluxes
  // associated with it.  Store Edot in c->dU[ERG], to be multiplied
  // by dt later (since at this stage we don't know dt).  This
  // later multiplication is done in eqn->preprocess_data()
  //
  double t_cond = set_conduction_dt_and_Edot(par, grid);
  t_cond = SimPM.levels[0].sub_domain.global_operation_double(MIN, t_cond);
#ifndef NDEBUG
  if (t_cond < par.dt) {
    spdlog::info(
        "PARALLEL CONDUCTION IS LIMITING TIMESTEP: t_c={12.6e}, t_m={12.6e}, t_dyn={12.6e}",
        t_cond, t_mp, t_dyn);
  }
#endif
  par.dt = min(par.dt, t_cond);

  //
  // If using MHD with GLM divB cleaning, the following sets the
  // hyperbolic wavespeed.  If not, it does nothing.  By setting it
  // here and using t_dyn, we ensure that the hyperbolic wavespeed is
  // equal to the maximum signal speed on the grid, and not an
  // artificially larger speed associated with a shortened timestep.
  //
  // we always calculate timestep on finest level first.
  static double td = 0.0;
  if (l == par.grid_nlevels - 1)
    td = t_dyn;
  else
    td = min(td, t_dyn / pow(2.0, par.grid_nlevels - 1 - l));

  double cr = 0.0;
#ifdef PION_OMP
  #pragma omp parallel private(cr)
  {
#endif
    if (par.grid_nlevels == 1) {
      cr = 0.25 / par.dx;
      spatial_solver->Set_GLM_Speeds(td, par.dx, cr);
    }
    else {
      cr = 0.25 / par.levels[par.grid_nlevels - 1].dx;
      if (l == 0)
        spatial_solver->Set_GLM_Speeds(
            td, par.levels[par.grid_nlevels - 1].dx, cr);
    }
#ifdef PION_OMP
  }
#endif

  //
  // Check that the timestep doesn't increase too much between step,
  // and that it  won't bring us past the next output time or the end
  // of the simulation. This function operates on SimPM.dt, resetting
  // it to a smaller value if needed.
  //
  timestep_checking_and_limiting(par, l);

  //
  // Tell the solver class what the resulting timestep is.
  //
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(par.dt);
#ifdef PION_OMP
  }
#endif

#ifndef NDEBUG
  // Check that if my process has modified dt to get to either
  // an output time or finishtime, then all procs have done this too!
  t_dyn = par.dt;
  t_mp  = SimPM.levels[0].sub_domain.global_operation_double(MIN, t_dyn);
  if (!pconst.equalD(t_dyn, t_mp)) {
    spdlog::error(
        "synchonisation trouble in timesteps! {:12.6e} {:12.6e}", t_dyn, t_mp);
    exit(2);
  }
#endif  // NDEBUG

  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_pllel::initial_conserved_quantities(
    vector<class GridBaseClass *> &g)
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION
  class GridBaseClass *grid = g[0];
  std::vector<pion_flt> u(SimPM.nvar, 0.0);
  initERG = 0.;
  initMMX = initMMY = initMMZ = 0.;
  initMASS                    = 0.0;
  double dx                   = grid->DX();
  double dv                   = 0.0;
  class cell *c               = grid->FirstPt();
  do {
    if (c->isgd) {
      // cout <<"*** LEVEL "<<l<<", cell is a leaf: "<<c->isdomain<<"
      dv = spatial_solver->CellVolume(*c, dx);
      spatial_solver->PtoU(c->P.data(), u.data(), SimPM.gamma);
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
  } while ((c = grid->NextPt(*c)) != 0);

  // cout <<"(local quantities) ["<< initERG <<", ";
  // cout << initMMX <<", ";
  // cout << initMMY <<", ";
  // cout << initMMZ <<", ";
  // cout << initMASS <<"]\n";

  initERG  = SimPM.levels[0].sub_domain.global_operation_double(SUM, initERG);
  initMMX  = SimPM.levels[0].sub_domain.global_operation_double(SUM, initMMX);
  initMMY  = SimPM.levels[0].sub_domain.global_operation_double(SUM, initMMY);
  initMMZ  = SimPM.levels[0].sub_domain.global_operation_double(SUM, initMMZ);
  initMASS = SimPM.levels[0].sub_domain.global_operation_double(SUM, initMASS);

  spdlog::info(
      "(conserved quantities) [{}, {}, {}, {}, {}]", initERG, initMMX, initMMY,
      initMMZ, initMASS);

#endif  // TEST_CONSERVATION
  return (0);
}



// ##################################################################
// ##################################################################



int sim_control_pllel::check_energy_cons(vector<class GridBaseClass *> &g)
{
#ifdef TEST_CONSERVATION
  // Energy, and Linear Momentum in x-direction.
  class GridBaseClass *grid = g[0];
  std::vector<pion_flt> u(SimPM.nvar, 0.0);
  double nowERG  = 0.;
  double nowMMX  = 0.;
  double nowMMY  = 0.;
  double nowMMZ  = 0.;
  double nowMASS = 0.0;
  double totmom  = 0.0;
  double dx      = grid->DX();
  double dv      = 0.0;
  class cell *c  = grid->FirstPt();
  do {
    if (c->isgd) {
      dv = spatial_solver->CellVolume(*c, dx);
      spatial_solver->PtoU(c->P.data(), u.data(), SimPM.gamma);
      nowERG += u[ERG] * dv;
      nowMMX += u[MMX] * dv;
      nowMMY += u[MMY] * dv;
      nowMMZ += u[MMZ] * dv;
      nowMASS += u[RHO] * dv;
      totmom += sqrt(u[MMX] * u[MMX] + u[MMY] * u[MMY] + u[MMZ] * u[MMZ]) * dv;
    }
  } while ((c = grid->NextPt(*c)) != 0);

  // cout <<"(local quantities) ["<< nowERG <<", ";
  // cout << nowMMX <<", ";
  // cout << nowMMY <<", ";
  // cout << nowMMZ <<", ";
  // cout << nowMASS <<"]\n";

  nowERG  = SimPM.levels[0].sub_domain.global_operation_double(SUM, nowERG);
  nowMMX  = SimPM.levels[0].sub_domain.global_operation_double(SUM, nowMMX);
  nowMMY  = SimPM.levels[0].sub_domain.global_operation_double(SUM, nowMMY);
  nowMMZ  = SimPM.levels[0].sub_domain.global_operation_double(SUM, nowMMZ);
  nowMASS = SimPM.levels[0].sub_domain.global_operation_double(SUM, nowMASS);
  totmom  = SimPM.levels[0].sub_domain.global_operation_double(SUM, totmom);
  // cout <<" totmom="<<totmom<<" initMMX="<<initMMX;
  // cout <<", nowMMX="<<nowMMX<<"\n";

  spdlog::info(
      "(conserved quantities) [{}, {}, {}, {}, {}]", nowERG, nowMMX, nowMMY,
      nowMMZ, nowMASS);
  spdlog::info(
      "(relative error      ) [{}, {}, {}, {}, {}]",
      (nowERG - initERG) / (initERG), (nowMMX - initMMX) / (totmom),
      (nowMMY - initMMY) / (totmom), (nowMMZ - initMMZ) / (totmom),
      (nowMASS - initMASS) / initMASS);
#endif  // TEST_CONSERVATION
  return (0);
}



#endif  // PARALLEL
