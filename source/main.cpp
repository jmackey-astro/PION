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
/// Arguments: \<main\> \<icfile\> \<typeoffile\> \<solver-type\> [override options]\n
/// Parameters:
/// - \<icfile\> Can be an ASCII text parameter-file for 1D and 2D shocktube
/// test problems; otherwise should be an initial-condition file or a
/// restartfile in fits format.
/// - \<typeoffile\> Integer flag to tell me what type of file I am starting from.
/// Can be one of [1=text paramfile, 2=fits restart/ICfile, 3=fitstable file].
/// - \<solvetype\> Integer =1 for uniform finite-volume, no other options.
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
/// - 2015.01.(10-16) JM: New include statements for new file
///    structure, and non-global grid class.

#include <iostream>
using namespace std;

//
// These tell code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_constants.h"

//
// Global variables
//
#include "global.h"

//
// reporting class, for dealing with stdio/stderr
//
#include "tools/reporting.h"

//
// grid base class
//
#include "grid/grid_base_class.h"
//
// simulation control toolkit class.
//
#include "grid.h"


int main(int argc, char **argv)
{
  
  if (argc<4) {
    print_command_line_options(argc,argv);
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
  cout <<"---------   pion v.0.1  running   ---------------------\n";
  cout <<"-------------------------------------------------------\n\n";

  // Set what type of file to open: 1=parameterfile, 2=restartfile.
  int ft;
  ft=atoi(argv[2]);
  if(ft <1 || ft>5) {cerr<<"(*pion*) Bad file type specifier.\n";return(1);}
  switch (ft) {
   case 1:
    cout <<"(pion) ft = "<<ft<<" so reading ICs from text parameterfile "<<argv[2]<<"\n";
    break;
   case 2:
   case 3:
    cout <<"(pion) ft = "<<ft<<" so reading ICs from Fits ICfile "<<argv[2]<<"\n";
    break;
  case 5:
    cout <<"(pion) ft = "<<ft<<" so reading ICs from SILO IC file "<<argv[2]<<"\n";
    break;
   default:
    rep.error("Bad filetype input to main",ft);
  }

  int type = atoi(argv[3]);
  if (type!=1) rep.error("Only know uniform FV solver (=1)",type);
 
  //
  // set up pointer to grid base class.
  //
  class GridBaseClass *grid = 0;

  //
  // Set up simulation controller class.
  //
  class IntegratorBaseFV *sim_control = 0;

  sim_control = new class IntUniformFV();
  if (!sim_control)
    rep.error("(pion) Couldn't initialise IntUniformFV sim_control", sim_control);

  // inputs are infile_name, infile_type, nargs, *args[]
  err = sim_control->Init(argv[1], ft, argc, args, &grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Init"<<"\n";
    delete sim_control;
    return 1;
  }
  err+= sim_control->Time_Int(grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Time_Int"<<"\n";
    delete sim_control;
    delete grid;
    return 1;
  }
  err+= sim_control->Finalise(grid);
  if (err!=0){
    cerr<<"(*pion*) err!=0 from Finalise"<<"\n";
    delete sim_control;
    delete grid;
    return 1;
  }

  delete sim_control; sim_control=0;
  if (grid) {delete grid; grid=0;}

  delete [] args; args=0;
  
  cout <<"-------------------------------------------------------\n";
  cout <<"---------   pion v.0.1  finsihed  ---------------------\n";
  cout <<"-------------------------------------------------------\n";

  return 0;
}
