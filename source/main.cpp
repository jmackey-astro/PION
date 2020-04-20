///
/// \file main.cpp
/// 
/// \brief Main program which sets up a uniform grid and runs the simulation.
/// 
/// \author Jonathan Mackey
/// 
/// This file just contains the main() function, which sets a number of 
/// parameters based on command line arguments, then initialised the grid, 
/// starts the time integration, and cleans up when the simulation is 
/// finished.
/// 
/// Arguments: \<main\> \<icfile\> [override options]\n
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

#include <iostream>
using namespace std;

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_constants.h"
#include "tools/reporting.h"
#include "grid/grid_base_class.h"
#include "sim_control/sim_control.h"


int main(int argc, char **argv)
{
  
  //
  // Set up simulation controller class.
  //
  class sim_control *sim_control = 0;
  sim_control = new class sim_control();
  if (!sim_control)
    rep.error("(pion) Couldn't initialise sim_control", sim_control);

  //
  // Check that command-line arguments are sufficient.
  //
  if (argc<2) {
    sim_control->print_command_line_options(argc,argv);
    rep.error("Bad arguments",argc);
  }
  
  int err=0;
  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) args[i] = argv[i];

  // Set up reporting class.
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      cout <<"Redirecting stdout to "<<outpath<<"info.txt"<<"\n";
      // Redirects cout and cerr to text files in the directory specified.
      rep.redirect(outpath);
    }
  }
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion UG SERIAL v2.0  running   ------------\n";
  cout <<"-------------------------------------------------------\n\n";

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
    //delete grid;
    return 1;
  }

  delete sim_control; sim_control=0;
  delete grid[0];
  delete [] args; args=0;
  
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion v.1.0  finished  ---------------------\n";
  cout <<"-------------------------------------------------------\n";

  return 0;
}





