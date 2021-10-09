/// \file main_NG.cpp
/// \brief Main program which sets up a NG grid and runs the simulation.
/// \author Jonathan Mackey
///
/// This file just contains the main() function, which sets a number of
/// parameters based on command line arguments, then initialised the grid,
/// starts the time integration, and cleans up when the simulation is
/// finished.
///
/// Arguments: \<pion_NG\> \<icfile\> [override options]\n
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
/// - 2018.06.12 JM: wrote code.

#include <iostream>
using namespace std;

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid/grid_base_class.h"
#include "sim_constants.h"
#include "sim_control/sim_control_NG.h"
#ifdef PARALLEL
#include "sim_control/sim_control_NG_MPI.h"
#endif
#include "tools/reporting.h"
#ifdef PION_OMP
#include <omp.h>
#endif

int main(int argc, char **argv)
{

  int err = 0;
#ifdef PARALLEL
  err        = COMM->init(&argc, &argv);
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
#endif

  //
  // Set up simulation controller class.
  //
#ifndef PARALLEL
  class sim_control_NG *sim_control = new class sim_control_NG();
#else
  class sim_control_NG_MPI *sim_control = new class sim_control_NG_MPI();
#endif
  if (!sim_control)
    rep.error("(pion) Couldn't initialise sim_control", sim_control);

  //
  // Check that command-line arguments are sufficient.
  //
  if (argc < 2) {
    sim_control->print_command_line_options(argc, argv);
    rep.error("Bad arguments", argc);
  }

  string *args = new string[argc];
  for (int i = 0; i < argc; i++)
    args[i] = argv[i];

  // Set up reporting class.
  // Redirects cout and cerr to text files in the directory specified.
  for (int i = 0; i < argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
#ifdef PARALLEL
      ostringstream path;
      path << outpath << "_" << myrank << "_";
      outpath = path.str();
      if (myrank == 0) {
        cout << "\tRedirecting stdout to " << outpath << "info.txt"
             << "\n";
      }
#else
      cout << "Redirecting stdout to " << outpath << "info.txt"
           << "\n";
#endif
      rep.redirect(outpath);
    }
  }
#ifdef PARALLEL
#ifdef NDEBUG
  rep.kill_stdout_from_other_procs(0);
#endif
#endif

  cout << "-------------------------------------------------------\n";
#ifdef PARALLEL
  cout << "---------   pion NG MPI v2.0 running   ----------------\n";
#else
  cout << "---------   pion NG SERIAL v.2.0  running   -----------\n";
#endif
  cout << "-------------------------------------------------------\n\n";

#ifdef PION_OMP
  // set number of OpenMP threads, if included
  int nth = 100;  // set to large number initially
  for (int i = 0; i < argc; i++) {
    if (args[i].find("omp-nthreads=") != string::npos) {
      nth = atoi((args[i].substr(13)).c_str());
      if (nth > omp_get_num_procs()) {
        cout << "\toverride: requested too many threads.\n";
        nth = min(nth, omp_get_num_procs());
      }
      cout << "\toverride: setting OpenMP N-threads to " << nth << "\n";
    }
    else
      nth = 1;
  }
  nth = min(nth, omp_get_num_procs());
  omp_set_num_threads(nth);
#endif

  // Set what type of file to open: 1=parameterfile, 2/5=restartfile.
  int ft;
  if (args[1].find(".silo") != string::npos) {
    cout << "(pion) reading ICs from SILO IC file " << args[1] << "\n";
    ft = 5;
  }
  else if (args[1].find(".fits") != string::npos) {
    cout << "(pion) reading ICs from Fits ICfile " << args[1] << "\n";
    ft = 2;
  }
  else {
    cout << "(pion) IC file not fits/silo: assuming text parameterfile "
         << args[1] << "\n";
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
      if (isnan(tmp) || isinf(tmp) || tmp < 0.0)
        rep.error("Don't recognise max walltime as a valid runtime!", tmp);

      sim_control->set_max_walltime(tmp * 3600.0);

#ifdef PARALLEL
      if (myrank == 0) {
#endif
        cout << "\tResetting MAXWALLTIME to ";
        cout << sim_control->get_max_walltime() << " seconds, or ";
        cout << sim_control->get_max_walltime() / 3600.0 << " hours.\n";
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
    cerr << "(*pion*) err!=0 from Init"
         << "\n";
    delete sim_control;
#ifdef PARALLEL
    cout << "rank: " << myrank << " nproc: " << nproc << "\n";
    delete COMM;
#endif
    return 1;
  }
  //
  // Integrate forward in time until the end of the calculation.
  //
  err += sim_control->Time_Int(grid);
  if (err != 0) {
    cerr << "(*pion*) err!=0 from Time_Int"
         << "\n";
    delete sim_control;
#ifdef PARALLEL
    cout << "rank: " << myrank << " nproc: " << nproc << "\n";
    delete COMM;
#endif
    return 1;
  }
  //
  // Finalise and exit.
  //
  err += sim_control->Finalise(grid);
  if (err != 0) {
    cerr << "(*pion*) err!=0 from Finalise"
         << "\n";
    delete sim_control;
#ifdef PARALLEL
    cout << "rank: " << myrank << " nproc: " << nproc << "\n";
    delete COMM;
#endif
    return 1;
  }

  delete sim_control;
  sim_control = 0;
  for (unsigned int v = 0; v < grid.size(); v++) {
    delete grid[v];
  }
  delete[] args;
  args = 0;
#ifdef PARALLEL
  cout << "rank: " << myrank << " nproc: " << nproc << "\n";
  delete COMM;
  COMM = 0;
#endif

  cout << "-------------------------------------------------------\n";
#ifdef PARALLEL
  cout << "---------  pion NG MPI v2.0 finished   ----------------\n";
#else
  cout << "---------  pion NG v.2.0  finished --------------------\n";
#endif
  cout << "-------------------------------------------------------\n";

  return 0;
}
