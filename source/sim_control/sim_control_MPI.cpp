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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"

#include "decomposition/MCMD_control.h"
#include "raytracing/raytracer_SC.h"
#include "sim_control/sim_control_MPI.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"
#endif  // if FITS

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef PARALLEL

// ##################################################################
// ##################################################################

sim_control_pllel::sim_control_pllel() : sim_control()
{
#ifdef TESTING
  cout << "sim_control_pllel constructor.\n";
#endif
}

// ##################################################################
// ##################################################################

sim_control_pllel::~sim_control_pllel()
{
#ifdef TESTING
  cout << "sim_control_pllel destructor.\n";
#endif
}

// ##################################################################
// ##################################################################

int sim_control_pllel::Init(
    string infile,
    int typeOfFile,
    int narg,
    string* args,
    vector<class GridBaseClass*>& grid  ///< address of vector of grid pointers.
)
{
#ifdef TESTING
  cout << "(sim_control_pllel::init) Initialising grid: infile = " << infile
       << "\n";
#endif
  int err = 0;

  //
  // Setup the MCMDcontrol class with rank and nproc.
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);

  //
  // Setup dataI/O class and check if we read from a single file or
  // multiple files.  Should be a single file, but for FITS it might
  // not be.
  //
  setup_dataio_class(SimPM, typeOfFile);
  if (!dataio->file_exists(infile)) rep.error("infile doesn't exist!", infile);

  if (typeOfFile == 2) {
    // FITS
    string::size_type pos = infile.find("_0000.");
    if (pos == string::npos) {
      SimPM.levels[0].MCMD.ReadSingleFile = true;
    }
    else {
      SimPM.levels[0].MCMD.ReadSingleFile = false;
      ostringstream t;
      t.str("");
      t << "_";
      t.width(4);
      t.fill('0');
      t << SimPM.levels[0].MCMD.get_myrank() << ".";
      string t2 = t.str();
      infile.replace(pos, 6, t2);
    }
  }
  else if (typeOfFile == 5) {
    // SILO
    string::size_type pos = infile.find("_0000.");
    if (pos == string::npos) {
      SimPM.levels[0].MCMD.ReadSingleFile = true;
    }
    else {
      SimPM.levels[0].MCMD.ReadSingleFile = false;
    }
  }
  else
    rep.error("bad input file type in Init", typeOfFile);

  //
  // We need to decompose the domain here, because setup_grid() needs
  // this, but this means we need to read the header to find out what
  // the grid dimensions are.  So the header is read twice, but this
  // should be ok because it only happens during initialisation.
  //
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("PLLEL Init(): failed to read header", 0, err);

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

  err = SimPM.levels[0].MCMD.decomposeDomain(SimPM, SimPM.levels[0]);
  rep.errorTest("PLLEL Init():Couldn't Decompose Domain!", 0, err);

  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  rep.errorTest("(INIT::override_params) err!=0 Something went wrong", 0, err);

  // Now set up the grid structure.
  grid.resize(1);
  // cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  err = setup_grid(grid, SimPM);
  // cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  SimPM.dx             = grid[0]->DX();
  SimPM.levels[0].grid = grid[0];
  rep.errorTest("(INIT::setup_grid) err!=0 Something went wrong", 0, err);

  //
  // All grid parameters are now set, so I can set up the appropriate
  // equations/solver class.
  //
  err = set_equations(SimPM);
  rep.errorTest("(INIT::set_equations) err!=0 Fix me!", 0, err);
  spatial_solver->SetEOS(SimPM.gamma);

  //
  // Now setup Microphysics, if needed.
  //
  err = setup_microphysics(SimPM);
  rep.errorTest("(INIT::setup_microphysics) err!=0", 0, err);

  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest(
      "(INIT::assign_initial_data) err!=0 Something went wrong", 0, err);
  //  cout <<"Read data finished\n";

  //
  // Set Ph[] = P[], and then implement the boundary conditions.
  //
  cell* c = grid[0]->FirstPt();
  do {
    for (int v = 0; v < SimPM.nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = grid[0]->NextPt(c)) != 0);

  //
  // If I'm using the GLM method, make sure Psi variable is
  // initialised to zero.
  //
  if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
#ifdef TESTING
    cout << "Initial state, zero-ing glm variable.\n";
#endif
    c = grid[0]->FirstPt();
    do {
      c->P[SI] = c->Ph[SI] = 0.;
    } while ((c = grid[0]->NextPt(c)) != 0);
  }

  //
  // Assign boundary conditions to boundary points.
  //
  //  cout <<"starting BC setup\n";
  err = boundary_conditions(SimPM, grid);
  rep.errorTest("(INIT::boundary_conditions) err!=0", 0, err);

  //  cout <<"assigning BC data\n";
  err = assign_boundary_data(SimPM, 0, grid[0]);
  rep.errorTest("(INIT::assign_boundary_data) err!=0", 0, err);

  //
  // Setup Raytracing on each grid, if needed.
  //
  //  cout <<"Setting up raytracing\n";
  err += setup_raytracing(SimPM, grid[0]);
  // cout <<"Setting up RT sources\n";
  err += setup_evolving_RT_sources(SimPM);
  // cout <<"Updating evolving RT sources\n";
  err += update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  rep.errorTest("Failed to setup raytracer and/or microphysics", 0, err);

  //
  // If testing the code, this calculates the momentum and energy on the
  // domain.
  //
  //  cout <<"initial conserved quantities\n";
  initial_conserved_quantities(grid[0]);

  err += TimeUpdateInternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, SimPM.tmOOA,
      SimPM.tmOOA);
  err += TimeUpdateExternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, SimPM.tmOOA,
      SimPM.tmOOA);
  if (err) rep.error("first_order_update: error from bounday update", err);

  //
  // If using opfreq_time, set the next output time correctly.
  //
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
    if (!dataio) rep.error("INIT:: dataio initialisation", SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

  if (SimPM.timestep == 0) {
#ifdef TESTING
    cout << "(PARALLEL INIT) Writing initial data.\n";
#endif
    output_data(grid);
  }

#ifdef TESTING
  c = (grid[0])->FirstPt_All();
  do {
    if (pconst.equalD(c->P[RO], 0.0)) {
      cout << "zero data in cell: ";
      CI.print_cell(c);
    }
  } while ((c = (grid[0])->NextPt_All(c)) != 0);
#endif  // TESTING
  cout << "------------------------------------------------------------\n";
  return (err);
}

// ##################################################################
// ##################################################################

/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int sim_control_pllel::Time_Int(
    vector<class GridBaseClass*>& grid  ///< vector of grids.
)
{
  cout << "-------------------------------------------------------\n";
  cout << "(sim_control_pllel::time_int) STARTING TIME INTEGRATION."
       << "\n";
  cout << "-------------------------------------------------------\n";
  int err = 0;
  err     = update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  rep.errorTest("TIME_INT:: initial RT src update()", 0, err);
  err = RT_all_sources(SimPM, grid[0], 0);
  rep.errorTest("TIME_INT:: initial RT()", 0, err);
  int log_freq  = 10;
  SimPM.maxtime = false;
  clk.start_timer("time_int");
  double tsf = 0.0;
  while (SimPM.maxtime == false) {
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
    rep.errorTest("TIME_INT::update_RT_sources()", 0, err);
    err = RT_all_sources(SimPM, grid[0], 0);
    rep.errorTest("TIME_INT:: loop RT()", 0, err);

    //
    // Update boundary data.
    //
#ifdef TESTING
    cout << "MPI time_int: updating internal boundaries\n";
#endif
    err += TimeUpdateInternalBCs(
        SimPM, 0, grid[0], spatial_solver, SimPM.simtime, OA2, OA2);
#ifdef TESTING
    cout << "MPI time_int: updating external boundaries\n";
#endif
    err += TimeUpdateExternalBCs(
        SimPM, 0, grid[0], spatial_solver, SimPM.simtime, OA2, OA2);
    if (err) rep.error("Boundary update at start of full step", err);

      // clk.start_timer("advance_time");
#ifdef TESTING
    cout << "MPI time_int: calculating dt\n";
#endif
    SimPM.levels[0].last_dt = SimPM.last_dt;
    err += calculate_timestep(SimPM, grid[0], spatial_solver, 0);
    rep.errorTest("TIME_INT::calc_timestep()", 0, err);

#ifdef TESTING
    cout << "MPI time_int: stepping forward in time\n";
#endif
    advance_time(0, grid[0]);
    // cout <<"advance_time took "<<clk.stop_timer("advance_time")<<"
    // secs.\n";
#ifdef TESTING
    cout << "MPI time_int: finished timestep\n";
    log_freq = 1;
#endif

    if ((SimPM.levels[0].MCMD.get_myrank() == 0)
        && (SimPM.timestep % log_freq) == 0) {
      cout << "New time: " << SimPM.simtime;
      cout << "\t dt=" << SimPM.dt;
      cout << "\t steps: " << SimPM.timestep;
      tsf = clk.time_so_far("time_int");
      cout << "\t runtime: " << tsf << " s"
           << "\n";
#ifdef TESTING
      cout.flush();
#endif  // TESTING
    }

    //
    // check if we are at time limit yet.
    //
    tsf         = clk.time_so_far("time_int");
    double maxt = COMM->global_operation_double("MAX", tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime = true;
      cout << "RUNTIME>" << get_max_walltime() << " SECS.\n";
    }

    err += output_data(grid);
    if (err != 0) {
      cerr << "(TIME_INT::output_data) err!=0 Something went bad"
           << "\n";
      return (1);
    }

    err += check_eosim();
    if (err != 0) {
      cerr << "(TIME_INT::) err!=0 Something went bad"
           << "\n";
      return (1);
    }
  }
  cout << "(sim_control_pllel::time_int) FINISHED - FINALISING SIM"
       << "\n";
  tsf = clk.time_so_far("time_int");
  cout << "TOTALS: Nsteps: " << SimPM.timestep;
  cout << ", sim-time: " << SimPM.simtime;
  cout << ", wall-time: " << tsf;
  cout << ", time/step: " << tsf / static_cast<double>(SimPM.timestep) << "\n";
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
    cout << "TOTALS RT: active: " << run << " idle: " << wait
         << " total: " << total << "\n";
  }
#endif
  cout << "                               "
          "**************************************\n\n";
  return (0);
}

// ##################################################################
// ##################################################################

int sim_control_pllel::calculate_timestep(
    class SimParams& par,             ///< simulation parameters
    class GridBaseClass* grid,        ///< pointer to grid.
    class FV_solver_base* sp_solver,  ///< solver/equations class
    const int l                       ///< level to advance (for NG grid)
)
{
  //
  // First get the local grid dynamics and microphysics timesteps.
  //
  double t_dyn = 0.0, t_mp = 0.0;
  t_dyn = calc_dynamics_dt(par, grid, sp_solver);
  t_mp  = calc_microphysics_dt(par, grid, l);

#ifdef TESTING
  cout << "calc_time: local t_dyn= " << t_dyn;
#endif

  //
  // Get global min over all grids on this level.
  //
  t_dyn = COMM->global_operation_double("MIN", t_dyn);
  t_mp  = COMM->global_operation_double("MIN", t_mp);

#ifdef TESTING
  cout << " , global t_dyn= " << t_dyn << endl;
  // Write step-limiting info every tenth timestep.
  if (t_mp < t_dyn && (par.timestep % 10) == 0)
    cout << "Limiting timestep by MP: mp_t=" << t_mp << "\thydro_t=" << t_dyn
         << "\n";
#endif

  par.dt = min(t_dyn, t_mp);

#ifdef THERMAL_CONDUCTION
  //
  // Calculate the timestep limit imposed by thermal conduction,
  // and calcuate the multidimensional energy fluxes
  // associated with it.  Store Edot in c->dU[ERG], to be multiplied
  // by dt later (since at this stage we don't know dt).  This
  // later multiplication is done in eqn->preprocess_data()
  //
  double t_cond = calc_conduction_dt_and_Edot();
  t_cond        = COMM->global_operation_double("MIN", t_cond);
  if (t_cond < par.dt) {
    cout << "PARALLEL CONDUCTION IS LIMITING TIMESTEP: t_c=";
    cout << t_cond << ", t_m=" << t_mp;
    cout << ", t_dyn=" << t_dyn << "\n";
  }
  par.dt = min(par.dt, t_cond);
#endif  // THERMAL CONDUCTION

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
  sp_solver->Setdt(par.dt);

#ifdef TESTING
  //
  // Check that if my processor has modified dt to get to either
  // an output time or finishtime, then all procs have done this too!
  //
  t_dyn = par.dt;
  t_mp  = COMM->global_operation_double("MIN", t_dyn);
  if (!pconst.equalD(t_dyn, t_mp))
    rep.error("synchonisation trouble in timesteps!", t_dyn - t_mp);
#endif  // TESTING

  return 0;
}

// ##################################################################
// ##################################################################

#endif  // PARALLEL
