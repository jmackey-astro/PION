/// \file mainMPI.cc
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
/// Arguments: \<main_parallel\> \<icfile\> \<typeoffile\> \<solver-type\> [override options]\n
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
///    [- artvisc=D : modify artificial viscosity, =0 none, Otherwise AVFalle with eta=D]
///    - AVtype=N  : modify type of AV: 0=none, 1=FKJ98, 2=CW84
///    - EtaVisc=D : modify viscosity parameter to double precision value.
///    - noise=D   : add noise to icfile if desired, at fractional level of D
///    - finishtime=D : set time to finish simulation, in code time units.
///    - coordsys=NAME : set coordinate system to [cartesian,cylindrical]. DANGEROUS TO OVERRIDE!
///    - cfl=D     : change the CFL no. for the simulation, in range (0,1).
///    - maxwalltime=D : change the max. runtime to D in seconds.
/// 
/// Written on 2007-10-20 (or so).
/// 
/// Modifications:
/// - 2007-11-08 made it better.
/// - 2010.09.30 JM: Added new flags for Artificial viscosity.
/// - 2010.10.13 JM: Moved print commandline options to global
///    function.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2013.01.17 JM: Made simulation initialisation less verbose.

#include <iostream>
#include <sstream>
using namespace std;
#include "global.h"
#include "grid.h"

#ifndef PARALLEL
# error "PARALLEL not defined, but trying to compile the mpi version!"
#endif

int main(int argc, char **argv)
{

  int err = COMM->init(&argc,&argv);
  //cout <<"argc="<<argc<<"\n";

  if (argc<4) {
    print_command_line_options(argc,argv);
    rep.error("Bad arguments",argc);
  }

  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) {
    args[i] = argv[i];
    //  cout <<"arg "<<i<<" = "<<args[i]<<"\n";
  }
  
  //
  // Redirect stdout/stderr if required.
  //
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      ostringstream path; path << outpath <<"_"<<mpiPM.myrank<<"_";
      outpath = path.str();
      if (mpiPM.myrank==0) {
        cout <<"\tRedirecting stdout to "<<outpath<<"info.txt"<<"\n";
      }
      rep.redirect(outpath); // Redirects cout and cerr to text files in the directory specified.
    }
  }

  //
  // Reset max. walltime to run the simulation for, if needed.
  //
  for (int i=0;i<argc; i++) {
    if (args[i].find("maxwalltime=") != string::npos) {
      double tmp = atof((args[i].substr(12)).c_str());
      if (isnan(tmp) || isinf(tmp) || tmp<0.0)
	rep.error("Don't recognise max walltime as a valid runtime!",tmp);
      mpiPM.set_max_walltime(tmp);
      if (mpiPM.myrank==0) {
        cout <<"\tResetting MAXWALLTIME to ";
        cout <<mpiPM.get_max_walltime()<<" seconds, or ";
        cout <<mpiPM.get_max_walltime()/3600.0<<" hours.\n";
      }
    }
  }

  cout << "rank: " << mpiPM.myrank << " nproc: " << mpiPM.nproc << "\n";
  for (int i=0;i<argc;i++) {
    cout <<"arg "<<i<<" = "<<args[i]<<"\n";
  }
  
  //char *tmp = new char [3];
  //char tmp[3];
  //strcpy(tmp,"hi jono, hows it going");
  //cout <<"tmp="<<tmp<<"\n";

  int ft;
  ft=atoi(argv[2]);
  if(ft <1 || ft>5) {cerr<<"(MAIN) Bad file type specifier.\n";return(1);}
  switch (ft) {
  case 1:
    cout <<"(MAIN) ft = "<<ft<<" so reading ICs from text parameterfile "<<argv[2]<<"\n";
    break;
#ifdef FITS
  case 2:
  case 3:
    cout <<"(MAIN) ft = "<<ft<<" so reading ICs from Fits ICfile "<<argv[2]<<"\n";
    break;
#endif // if FITS
#ifdef SILO
  case 5:
    cout <<"(MAIN) ft = "<<ft<<" so reading ICs from Silo ICfile "<<argv[2]<<"\n";
    break;
#endif // if SILO
  default:
    rep.error("Bad filetype input to main",ft);
  }
  
  if (integrator) rep.error("integrator already set up!",integrator);
  integrator = new class ParallelIntUniformFV();
  if (!integrator) rep.error("(MAIN) Couldn't initialise ParallelIntUniformFV integrator", integrator);
  //IntegratorBaseFV *solver=0;
  //solver = new class ParallelIntUniformFV();
  //if (!solver) rep.error("Couldn't initialise ParallelIntUniformFV integrator",solver);

  if(!integrator) {cerr<<"(MAIN) Error initialising grid...exiting."<<"\n";delete integrator; return(1);}
  err = integrator->Init(argv[1], ft, argc, args);
  if (err!=0){cerr<<"(MAIN) err!=0 Something went bad"<<"\n";delete integrator;return(1);}
  err+= integrator->Time_Int();
  if (err!=0){cerr<<"(MAIN) err!=0 Something went bad"<<"\n";delete integrator;return(1);}
  err+= integrator->Finalise();
  if (err!=0){cerr<<"(MAIN) err!=0 Something went bad"<<"\n";delete integrator;return(1);}

  delete integrator; integrator=0;
  
  COMM->finalise();
  cout << "rank: " << mpiPM.myrank << " nproc: " << mpiPM.nproc << "\n";
  delete COMM; COMM=0;

  delete [] args; args=0;
  return(0);
}

