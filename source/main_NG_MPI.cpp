/// \file main_NG_MPI.cpp
/// \brief Main program which sets up a NG grid and runs a MCMD sim.
/// \author Jonathan Mackey
/// 
/// 
/// Arguments: \<pion_NG_mpi\> \<icfile\> [override options]\n
/// Parameters:
/// - \<icfile\> Can be an ASCII text parameter-file for 1D and 2D shocktube
/// test problems; otherwise should be an initial-condition file or a
/// restartfile in fits/silo format.
/// - [override options] are optional and of the format \<name\>=\<value\> with no spaces.
///    - redirect=PATH: Path to redirect stdout/stderr to, (with trailing forward slash).
///    - opfreq=N  : modify output frequency to every Nth timestep.
///    - optype=N  : modify type of output file, 1=TXT,2=FITS,3=FitsTable,4=TXT+FITS.
///    - outfile=NAME : Replacement output filename, with path.
///    - ooa=N     : modify order of accuracy.
///    - eqntype=N : modify type of equations, 1=Euler, 2=idealMHD, 3=glmMHD
///    [- artvisc=D : modify artificial viscosity, =0 none, Otherwise AVFalle with eta=D,]
///    - AVtype=N  : modify type of AV: 0=none, 1=FKJ98, 2=CW84
///    - EtaVisc=D : modify viscosity parameter to double precision value.
///    - noise=D   : add noise to icfile if desired, at fractional level of D
///    - finishtime=D : set time to finish simulation, in code time units.
///    - coordsys=NAME : set coordinate system to [cartesian,cylindrical]. DANGEROUS TO OVERRIDE!
///    - cfl=D     : change the CFL no. for the simulation, in range (0,1).
///    - cooling=N : cooling=0 for no cooling, 1 for Sutherland&Dopita1993.
/// 
/// Modifications:
/// - 2018.08.30 JM: wrote code.

#include <iostream>
#include <sstream>
using namespace std;

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_constants.h"
#include "tools/reporting.h"
#include "grid/grid_base_class.h"
#include "sim_control/sim_control_NG_MPI.h"


int main(int argc, char **argv)
{
  

  int err = COMM->init(&argc,&argv);
  //cout <<"argc="<<argc<<"\n";
  
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);

  //
  // Set up simulation controller class.
  //
  class sim_control_NG_MPI *sim_control = 0;

  sim_control = new class sim_control_NG_MPI();
  if (!sim_control)
    rep.error("(pion) Couldn't initialise sim_control", sim_control);

  //
  // Check that command-line arguments are sufficient.
  //
  if (argc<2) {
    sim_control->print_command_line_options(argc,argv);
    rep.error("Bad arguments",argc);
  }
  
  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) args[i] = argv[i];

  // Set up reporting class.
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
  //rep.kill_stdout_from_other_procs(0);
#endif
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion NG MPI v1.0 running   ----------------\n";
  cout <<"-------------------------------------------------------\n";

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
    cout <<"(pion) IC file not fits/silo: assuming parameter file";
    cout <<args[1]<<"\n";
    ft=1;
  }

  //
  // set up vector of grids.
  //
  vector<class GridBaseClass *> grid;

  //
  // Initialise the grid.
  // inputs are infile_name, infile_type, nargs, *args[]
  //
  err = sim_control->Init(args[1], ft, argc, args, grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Init"<<"\n";
    delete sim_control;
    return 1;
  }
  //
  // Integrate forward in time until the end of the calculation.
  //
  err+= sim_control->Time_Int(grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Time_Int"<<"\n";
    delete sim_control;
    //delete grid;
    return 1;
  }
  //
  // Finalise and exit.
  //
  err+= sim_control->Finalise(grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Finalise"<<"\n";
    delete sim_control;
    return 1;
  }

  delete sim_control; sim_control=0;
  for (unsigned int v=0; v<grid.size(); v++) {
    delete grid[v];
  }
  delete [] args; args=0;
  
  COMM->finalise();
  cout << "rank: " << myrank << " nproc: " << nproc << "\n";
  delete COMM; COMM=0;
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion NG MPI v1.0 finished  ----------------\n";
  cout <<"-------------------------------------------------------\n";

  return 0;
}





