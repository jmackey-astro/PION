/// \file mainMPI.cpp
/// 
/// \brief Main program which sets up a Parallel uniform grid and runs the simulation.
/// 
/// \author Jonathan Mackey
/// 
/// This file just contains the main() function, which sets a number of 
/// parameters based on command line arguments, then initialised the grid, 
/// starts the time integration, and cleans up when the simulation is 
/// finished.
/// 
/// Modifications:
/// - 2007-11-08 made it better.
/// - 2010.09.30 JM: Added new flags for Artificial viscosity.
/// - 2010.10.13 JM: Moved print commandline options to global
///    function.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2013.01.17 JM: Made simulation initialisation less verbose.
/// - 2015.01.08 JM: Moved grid definition to this file from global.h
/// - 2015.01.26 JM: updates, moving mpiPM from global to sim_control
/// - 2015.04.30 JM: tidying up.
/// - 2018.04.27 JM: removed some args (simpler command-line running).

#include <iostream>
#include <sstream>
using namespace std;

//
// These tell code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

//
// grid base class
//
#include "grid/grid_base_class.h"
//
// simulation control toolkit class.
//
#include "sim_control/sim_control.h"
#include "sim_control/sim_control_MPI.h"


#ifndef PARALLEL
# error "PARALLEL not defined, but trying to compile the mpi version!"
#endif

int main(int argc, char **argv)
{

  int err = COMM->init(&argc,&argv);
  //cout <<"argc="<<argc<<"\n";
  
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);

  //
  // Set up simulation controller class.
  //
  class sim_control_pllel *sim_control = 0;
  sim_control = new class sim_control_pllel();
  if (!sim_control)
    rep.error("(PION) Couldn't initialise sim_control_pllel", sim_control);

  //
  // Check that command-line arguments are sufficient.
  //
  if (argc<2) {
    sim_control->print_command_line_options(argc,argv);
    rep.error("Bad arguments",argc);
  }

  //
  // copy cmd-line args into an array of strings (for ease of use.
  //
  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) {
    args[i] = argv[i];
  }
  
  //
  // Redirect stdout/stderr if required.
  //
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      ostringstream path; path << outpath <<"_"<<myrank<<"_";
      outpath = path.str();
      if (myrank==0) {
        cout <<"\tRedirecting stdout to "<<outpath<<"info.txt"<<"\n";
      }
      // Redirects cout and cerr to text files in the directory specified.
      rep.redirect(outpath);
    }
  }
#ifndef TESTING
  rep.kill_stdout_from_other_procs(0);
#endif
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion v.1.0  running   ---------------------\n";
  cout <<"-------------------------------------------------------\n\n";


  cout << "rank: " << myrank << " nproc: " << nproc << "\n";
  for (int i=0;i<argc;i++) {
    cout <<"arg "<<i<<" = "<<args[i]<<"\n";
  }
  
  // Set what type of file to open: 1=parameterfile, 2/5=restartfile.
  int ft;
  if      (args[1].find(".silo") != string::npos) {
    cout <<"(pion) reading ICs from SILO IC file "<<args[1]<<"\n";
    ft=5;
  }
  else if (args[1].find(".fits") != string::npos) {
    cout <<"(pion) reading ICs from Fits ICfile "<<args[1]<<"\n";
    ft=2;
  }
  else {
    cout <<"(pion) IC file not fits/silo: assuming text parameterfile "<<args[1]<<"\n";
    ft=1;
  }
  
  //
  // set up pointer to grid base class.
  //
  vector<class GridBaseClass *> grid;

  //
  // Reset max. walltime to run the simulation for, if needed.
  // Input should be in hours.
  //
  for (int i=0;i<argc; i++) {
    if (args[i].find("maxwalltime=") != string::npos) {
      double tmp = atof((args[i].substr(12)).c_str());
      if (isnan(tmp) || isinf(tmp) || tmp<0.0)
	rep.error("Don't recognise max walltime as a valid runtime!",tmp);

      sim_control->set_max_walltime(tmp*3600.0);

      if (myrank==0) {
        cout <<"\tResetting MAXWALLTIME to ";
        cout <<sim_control->get_max_walltime()<<" seconds, or ";
        cout <<sim_control->get_max_walltime()/3600.0<<" hours.\n";
      }
    }
  }
  
  //
  // Initialise code.
  //
  err = sim_control->Init(args[1], ft, argc, args, grid);
  if (err!=0) {
    cerr<<"(PION) err!=0 from Init"<<"\n";
    delete sim_control;
    return(1);
  }
  //
  // Step forward in time until the end of the simulation.
  //
  err+= sim_control->Time_Int(grid);
  if (err!=0) {
    cerr<<"(PION) err!=0 from Time_Int"<<"\n";
    delete sim_control;
    return(1);
  }
  //
  // Finalise the simulation.
  //
  err+= sim_control->Finalise(grid);
  if (err!=0) {
    cerr<<"(PION) err!=0 from Finalise"<<"\n";
    delete sim_control;
    return(1);
  }

  delete grid[0];
  delete sim_control; sim_control=0;
  delete [] args; args=0;

  COMM->finalise();
  cout << "rank: " << myrank << " nproc: " << nproc << "\n";
  delete COMM; COMM=0;
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion v.1.0  finsihed  ---------------------\n";
  cout <<"-------------------------------------------------------\n";

  return(0);
}

