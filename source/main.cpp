///
/// \file main.cpp
///
/// \brief Main program which sets up a uniform grid and runs PION.
///
/// \author Jonathan Mackey
///
/// This file just contains the main() function, which sets a number of
/// parameters based on command line arguments, then initialised the grid,
/// starts the time integration, and cleans up when the simulation is
/// finished.  This is for the serial version of PION without a nested grid
///
/// Arguments: \<main\> \<icfile\> [override options]\n
/// Parameters:
/// - \<icfile\> Can be an ASCII text parameter-file for 1D and 2D shocktube
/// test problems; otherwise should be an initial-condition file or a
/// restartfile in fits/silo format.
/// - [override options] are optional and of the format \<name\>=\<value\> with
/// no spaces.
///    - redirect=PATH: Path to redirect stdout/stderr to, (with trailing
///    forward slash).
///    - opfreq=N  : modify output frequency to every Nth timestep.
///    - optype=N  : modify type of output file,
///    1=TXT,2=FITS,3=FitsTable,4=TXT+FITS.
///    - outfile=NAME : Replacement output filename, with path.
///    - ooa=N     : modify order of accuracy.
///    - eqntype=N : modify type of equations, 1=Euler, 2=idealMHD, 3=glmMHD
///    [- artvisc=D : modify artificial viscosity, =0 none, Otherwise AVFalle
///    with eta=D,]
///    - AVtype=N  : modify type of AV: 0=none, 1=FKJ98, 2=CW84
///    - EtaVisc=D : modify viscosity parameter to double precision value.
///    - noise=D   : add noise to icfile if desired, at fractional level of D
///    - finishtime=D : set time to finish simulation, in code time units.
///    - coordsys=NAME : set coordinate system to [cartesian,cylindrical].
///    DANGEROUS TO OVERRIDE!
///    - cfl=D     : change the CFL no. for the simulation, in range (0,1).
///    - cooling=N : cooling=0 for no cooling, 1 for Sutherland&Dopita1993.
///
/// Written on 2006-12-22-Friday
///
/// Modifications:
/// - 2007-06-22 Updated command-line args.
/// - 2007-06-26 Updated documentation and some args.
/// - 2007-07-13 New Class structure implemented.
/// - 2010.09.30 JM: Added new flags for Artificial viscosity.
/// - 2010.10.13 JM: Moved print_commandline_options to global function.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.08 JM: Moved grid definition to this file from global.h
///    and added link to reporting class.
/// - 2015.01.(10-26) JM: New include statements for new file
///    structure, and non-global grid class.
/// - 2015.04.30 JM: tidying up.
/// - 2018.04.27 JM: removed some args (simpler command-line running).

#include <sstream>
using namespace std;

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid/grid_base_class.h"
#include "sim_constants.h"

#ifdef PION_NESTED
#include "sim_control/sim_control_NG.h"
#ifdef PARALLEL
#include "sim_control/sim_control_NG_MPI.h"
#endif /* PARALLEL */
#else
#include "sim_control/sim_control.h"
#ifdef PARALLEL
#include "sim_control/sim_control_MPI.h"
#endif /* PARALLEL */
#endif /* PION_NESTED */

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#ifdef PION_OMP
#include <omp.h>
#endif

int main(int argc, char **argv)
{
  int err = 0;

#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  //
  // Set up simulation controller class.
  //
#ifdef PION_NESTED
  string grid_type = "NG";
#ifndef PARALLEL
  class sim_control_NG *sim_control = new class sim_control_NG();
#else
  class sim_control_NG_MPI *sim_control = new class sim_control_NG_MPI();
#endif /* PARALLEL */
#else
  string grid_type               = "UG";
#ifndef PARALLEL
  class sim_control *sim_control = new class sim_control();
#else
  class sim_control *sim_control = new class sim_control_pllel();
#endif /* PARALLEL */
#endif /* PION_NESTED */

  if (!sim_control)
    spdlog::error(
        "{}: {}", "(pion) Couldn't initialise sim_control",
        fmt::ptr(sim_control));

#ifdef PARALLEL
  int myrank = sim_control->SimPM.levels[0].sub_domain.get_myrank();
  int nproc  = sim_control->SimPM.levels[0].sub_domain.get_nproc();
  /* turn off logging for not root processes for Release build */
#ifdef NDEBUG
  if (myrank > 0) {
    spdlog::set_level(spdlog::level::off);
    spdlog::flush_on(spdlog::level::off);
  }
#endif /* NDEBUG */
#endif /* PARALLEL */

  //
  // Check that command-line arguments are sufficient.
  //
  if (argc < 2) {
    sim_control->print_command_line_options(argc, argv);
    spdlog::error("{}: {}", "Bad arguments", argc);
  }

  string *args = new string[argc];
  for (int i = 0; i < argc; ++i)
    args[i] = argv[i];
#ifndef NDEBUG
  for (int i = 0; i < argc; ++i) {
    spdlog::info("args: i= {}, arg = {}", i, args[i]);
  }
#endif
  for (int i = 0; i < argc; ++i) {
    if (args[i].find("redirect=") != string::npos) {
      ostringstream path;
      path << args[i].substr(9);
#ifdef PARALLEL
      path << "_" << myrank;
#endif
      path << ".log";
      auto max_logfile_size = 1048576 * 5;
      auto max_logfiles     = 3;
      spdlog::set_default_logger(spdlog::rotating_logger_mt(
          "pion", path.str(), max_logfile_size, max_logfiles));
    }
  }

#ifdef PION_OMP
  // set number of OpenMP threads, if included
  int nth        = 100;  // set to large number initially
  bool found_omp = false;
  for (int i = 0; i < argc; i++) {
    if (args[i].find("omp-nthreads=") != string::npos) {
      nth = atoi((args[i].substr(13)).c_str());
      if (nth > omp_get_num_procs()) {
        spdlog::warn("override: requested too many threads");
        nth = min(nth, omp_get_num_procs());
      }
      spdlog::warn("override: setting OpenMP N-threads to {}", nth);
      found_omp = true;
    }
  }
  if (found_omp) {
    omp_set_num_threads(nth);
  }
  else {
    omp_set_num_threads(1);
    nth = 1;
  }
#endif

#ifdef PARALLEL
  string parallelism = "MPI";
#ifdef PION_OMP
  parallelism += "/OpenMP";
#endif /* PION_OMP */
#elif PION_OMP
  string parallelism             = "OpenMP";
#else
  string parallelism = "SERIAL";
#endif /* PARALLEL */
#ifndef PARALLEL
  int nproc = 1;
#endif
#ifndef PION_OMP
  int nth = 1;
#endif

  spdlog::info("-------------------------------------------------------");
  spdlog::info(
      "---------   pion {} {} v2.0  running   ---------------", grid_type,
      parallelism);
  spdlog::info(
      "-------------- n-proc = {}, n-thread = {} --------------\n", nproc, nth);

  // Set what type of file to open: 1=parameterfile, 2/5=restartfile.
  int ft;
  if (args[1].find(".silo") != string::npos) {
    spdlog::info("(pion) reading ICs from SILO IC file {}", args[1]);
    ft = 5;
  }
  else if (args[1].find(".fits") != string::npos) {
    spdlog::info("(pion) reading ICs from Fits ICfile {}", args[1]);
    ft = 2;
  }
  else {
    spdlog::info(
        "(pion) IC file not fits/silo: assuming text parameterfile {}",
        args[1]);
    ft = 1;
  }
  //
  // set up vector of grids.
  //
  vector<class GridBaseClass *> grid;

  //
  // Reset max. walltime to run the simulation for, if needed.
  // Input should be in hours.
  //
  for (int i = 0; i < argc; i++) {
    if (args[i].find("maxwalltime=") != string::npos) {
      double tmp = atof((args[i].substr(12)).c_str());
      if (!isfinite(tmp) || tmp < 0.0)
        spdlog::error(
            "{}: {}", "Don't recognise max walltime as a valid runtime!", tmp);

      sim_control->set_max_walltime(tmp * 3600.0);

#ifdef PARALLEL
      if (myrank == 0) {
#endif
        spdlog::info(
            "\tResetting MAXWALLTIME to {} seconds, or {} hours.\n",
            sim_control->get_max_walltime(),
            sim_control->get_max_walltime() / 3600.0);
#ifdef PARALLEL
      }
#endif
    }
  }

  //
  // Initialise the grid.
  // inputs are infile_name, infile_type, nargs, *args[]
  //
  err = sim_control->Init(args[1], ft, argc, args, grid);
  if (err != 0) {
    spdlog::error("(PION) err!=0 from Init");
    delete sim_control;
    return 1;
  }
  //
  // Integrate forward in time until the end of the calculation.
  //
  err += sim_control->Time_Int(grid);
  if (err != 0) {
    spdlog::error("(PION) err!=0 from Time_Int");
    delete sim_control;
    return 1;
  }
  //
  // Finalise and exit.
  //
  err += sim_control->Finalise(grid);
  if (err != 0) {
    spdlog::error("(PION) err!=0 from Finalise");
    delete sim_control;
    return 1;
  }

  delete sim_control;
  sim_control = 0;
  for (auto g : grid)
    delete g;
  delete[] args;
  args = 0;

  spdlog::info(
      "-------------------------------------------------------\n---------   pion {} {} v2.0  finished ---------------\n-------------------------------------------------------",
      grid_type, parallelism);

  return 0;
}
