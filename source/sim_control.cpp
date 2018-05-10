/// \file sim_control.cpp
/// 
/// \brief Simulation Control Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for sim_control 
/// class, which is the basic 1st/2nd order Finite Volume Solver according to
/// the method outlined in Falle, Komissarov, \& Joarder (1998), MNRAS, 297, 265.
/// 
/// 
/// Modifications:
/// - 2007-06-26 Hunted for bugs, but they were in another file.
/// - 2007-07-13 New Class structure implemented.
/// - 2007-07-16 New Class structure works (modified from Friday a little).  Same as old, but a bit faster.
/// - 2007-07-17 2D riemann problem initial setup improved (area averaged interface cells).
/// - 2007-07-24 Added passive tracer variable support.
/// - 2007-08-07 cylindrical coordinates working for 2d axi-symmetry.
/// - 2007-10-11 Cleaning up.
/// - 2007-11-01 Added dataio class last week, cleaning up today.
/// - 2008-09-20 Removed text I/O into its own class. ifdeffed silo/fits.
///
/// - JM 2009-12-16 Added ifdef in sim_control::Time_Int() so that I
///      can get the code to output magnetic pressure instead of
///      timing info every timestep.  This is purely to make a plot of
///      magnetic pressure for the Field Loop Advection Test and
///      should be switched off in general (although it makes no
///      difference to the running of the code!).
/// - 2009-12-18 JM: Added Axisymmetric Class
///    (cyl_FV_solver_Hydro_Euler) in set_equations().
///
/// - 2010-01-05 JM: Added check for microphysics timestep limiting in
///     calc_timestep().  This is controlled by the flag
///     SimPM.EP.MP_timestep_limit, which is set to true to turn it
///     on.
/// - 2010-01-21 JM: Added override option for SimPM.gamma, the
///    equation of state parameter.
/// - 2010-04-10 JM: fixed width timestep in text and fits filenames.
/// - 2010-04-21 JM: Changed filename setup so that i can write
/// checkpoint files with fname.999999.txt/silo/fits.  Added a check
/// for checkpointing in output_data() and removed all of the outfile
/// string generation (moved to dataio classes).
/// - 2010-04-25 JM: Added an extra command-line parameter to set the
///  minimum allowed timestep.  If the step is shorter than this then
///  in calc_timestep() we will bug out because something has gone
///  seriously wrong.
/// - 2010-07-20 JM: Changed spOOA,tmOOA to integers.
/// - 2010-07-21 JM: Set SimPM.next_optime in ready_to_start() instead
///    of ReadHeader().
/// - 2010-07-23 JM: Swtiched checkpointing to use two files and overwrite
///    alternate files, ensuring there will always be at least one valid
///    checkpoint file to restart from.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.09.27 JM: took out comments and old ifdefs from dU_column().
/// - 2010.09.30 JM: Added div(V) calculation to calc_dU() function for viscosity.
/// - 2010.10.01 JM: Spherical coordinates(1D only) for Euler equations.
///       Cut out testing myalloc/myfree
/// - 2010.10.04 JM: Moved the field-loop magnetic pressure output to
///       an ifdeffed function.  Added a new function to calculate the
///       1D blast wave radius in spherical coordinates (also ifdeffed).
/// - 2010.10.05 JM: Changed boundary point timestep calculation to
///       work for all sims, not just jet sims.  Note this only works
///       for the XN point from the grid FirstPt(), not all boundary
///       points (which it really should work for!)
/// - 2010.10.13 JM: Added option to setup new MPv9 class
///       in setup_microphysics() function.
///       Added MP_timestep_limit override option.
/// - 2010.11.03 JM: Changed "endl" to "\n" for JUROPA.  Added a
///       digit to output file counters.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
/// - 2010.11.15 JM: Modified update_dynamics() so it has pre- and
///   post-processing functions before and after calc_dU().
///   Added routine to calculate the H-correction eta values for
///   pre-processing.  Added extra_data setting in setup_grid().
///   2010.11.19 JM: Debugged the H-corr stuff.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  Will probably keep it since it is the way
///   things will go eventually.
/// - 2010.12.27 JM: Enclosed isothermal solver in an ifdef (it's
///   broken at the moment and I have no time to fix it).
/// - 2010.12.28 JM: Added Internal energy integrators to the
///   set_equations() function, enclosed in a hash-ifdef for testing.
///   (30.12 and 31.12) More on internal energy integration.
/// - 2011.01.03 JM: Moved preprocess_dU() to solver_eqn_base.h/.cc
///   so that it can be re-defined for the internal energy update.
/// - 2011.01.17 JM: Added override for checkpt_freq=N steps.
///                  Added new mp_only_cooling() class in MP setup.
/// - 2011.02.17 JM: Added new optype==6 for outputting text+silo files.
///                  Raytracer header file moved to raytracing/ subdir.
///                  More ray-tracing options for CI.setup_extra_data
/// - 2011.02.25 JM: New setup_raytracing() and update_microphysics() functions.
///    The interface is simpler, and the logic is much clearer, and it should
///    now work for multiple sources.
/// - 2011.03.21 JM: moved cell setup_extra_data() to its own function, to save 
///    on copy-paste for the parallel version.
///    Rewrote setup_raytracing() and atomised update_microphysics() so there are
///    now a few new functions to deal with the different kinds of radiative
///    transfer we want to do.
/// - 2011.04.06 JM: Added thermal-conduction timestep limiting and flux calculation.
/// - 2011.04.14 JM: Added MPv2 microphysics integration class.
/// - 2011.04.17 JM: minor debugging additions/subtractions for the new RT update.
/// - 2011.04.18 JM: more debugging.  nearly done now.
/// - 2011.04.22 JM: bugs for multiple point sources fixed.  Also removed isfinite
///    checks for ints.  the intel compiler bugs out with them!  I guess they are
///    redundant if an int doesn't have a bit combination for NAN or INF.
/// - 2011.04.29 JM: changed logic: replaced some if (c->isbd) to if (!c->isgd)
///  because the two are no longer mutually exclusive (grid data can be internal
///  boundary data also.  Added a check in update_microphysics so that if a cell
///  is boundary data it is not updated (and in microphysics calc_timestep).
/// - 2011.05.02 JM: Added set_multifreq_source_properties() call to 
///    setup_microphysics() function.
/// - 2011.06.21 JM: Added option of 2nd-order-in-time ray-tracing/microphysics with
///    two microphysics updates per step.
/// - 2011.10.13 JM: Added MPv4 class.
/// - 2011.10.14 JM: Removed raytracer_shielding class, updated RT interface a lot.
/// - 2011.10.22 JM: The new MP/RT interface is now default (old #deffed code
///    is now removed as of SVN369).  Improved RT interface, and simplified the 
///    logic again, so it should now be easier to add to.
/// - 2012.01.16 JM: Added setup_evolving_RT_sources() and
///    update_evolving_RT_sources() for stellar evolution models.
/// - 2012.01.20 JM: Updated setup_evolving_RT_sources() to scale luminosity to
///    match luminosities of Diaz-Miller+(1998,Table 1).
/// - 2012.01.23 JM: Added ifdefs for microphysics classes.
/// - 2012.07.24 JM: Changed time-update to improve stability (photoionisation
///    models with R-type I-fronts were developing ripples/waves in solution).
///
/// - 2012.08.05 JM: Cut out lots of code, moved to time_integrators/
///    but for now it is just ifdeffed out.  Will move more once the
///    new time-integration scheme is validated/tested.
/// - 2012.08.16 JM: Debugging.  It seems to be working well now, but
///    there is still more testing to do.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.02.15 JM: Added NEW_METALLICITY flag for testing the new
///    microphysics classes.
/// - 2013.04.15 JM: Moved microphysics setup to early in Init from
///    ready_to_start() function, so that the stellar wind boundary
///    setup functions can call MP->Set_Temp().
/// - 2013.04.18 JM: Removed NEW_METALLICITY flag.
/// - 2013.08.23 JM: Added new mpv9_HHe module code.
/// - 2015.01.(10-16) JM: New include statements for new file
///    structure, and non-global grid class.
/// - 2015.01.26 JM: CHANGED FILENAME TO SIM_CONTROL.CPP
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "sim_control.h"

#include "microphysics/microphysics_base.h"

#include "raytracing/raytracer_SC.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <climits>
using namespace std;


#define TIMESTEP_FULL 1
#define TIMESTEP_FIRST_PART 2


// ##################################################################
// ##################################################################


sim_control::sim_control()
{
  spatial_solver=0;
  dataio=0;
  textio=0;
  SimPM.checkpoint_freq=INT_MAX;
  max_walltime = 1.0e100;
}



// ##################################################################
// ##################################################################


sim_control::~sim_control()
{
#ifdef TESTING
  cout << "(sim_control::Destructor) Deleting Grid Class..." <<"\n";
#endif
  if (spatial_solver) {delete spatial_solver; spatial_solver=0;}
  if (dataio) {delete dataio; dataio=0;}
  if (textio) {delete textio; textio=0;}
  if (RT[0]) delete RT[0];
#ifdef TESTING
  cout << "(sim_control::Destructor) Done." <<"\n";
#endif
}


// ##################################################################
// ##################################################################


double sim_control::get_max_walltime()
{
  return max_walltime;
}


// ##################################################################
// ##################################################################


void sim_control::set_max_walltime(
        double t ///< New Max. runtime in seconds.
        )
{
  cout <<"\tResetting max. walltime from "<<max_walltime;
  max_walltime = t;
  cout <<" to new value: "<<max_walltime<<"\n";
}


/*****************************************************************/
/********************* SIM. INITIALISATION ***********************/
/*****************************************************************/



// ##################################################################
// ##################################################################


//---------------------------------------------------------
//
// Function to output commandline options for code.
//
void sim_control::print_command_line_options(
        int argc,
        char **argv
        )
{
  cout <<"PION: You ran:\n";
  for (int v=0;v<argc;v++)
    cout <<"  "<<argv[v];
  cout <<"\n      ************************         \n";
  cout << argv[0] <<": must call with at least 3 arguments...\n";
  cout <<" <main> <icfile> <typeoffile> <solver-type> [optional args]\n";
  cout <<"Parameters:\n";
  cout <<"<icfile> \n";
  cout <<"\tCan be an ASCII parameter-file for 1D and 2D shocktubes.\n";
  cout <<"\tOtherwise should be an initial-condition file or restart-file\n";
  cout <<"\tin FITS or Silo format.\n";
  cout <<"<typeoffile> \n";
  cout <<"\tInteger flag to tell me what type of file I am starting from.\n";
  cout <<"\tCan be one of [1=text paramfile, 2=FITS, 5=Silo file].\n";
  cout <<"<solvetype> \n";
  cout <<"\tInteger =1 for uniform finite-volume, no other options.\n";

  cout <<"\n";
  cout <<"[optional args] are in the format <name>=<value> with no spaces.\n\n";

  cout <<"\n*********** DATA I/O OPTIONS ************\n";
  cout <<"\t redirect=string : filename with path to redirect stdout/stderr to\n";

  cout <<"\t checkpt_freq=N  : checkpoint every N timesteps (default is never).\n";
  cout <<"\t op_criterion=N  : 0=output every I steps, 1=output every D time units.\n";
  cout <<"\t opfreq=N        : Output data every Nth timestep  (if op_criterion=0).\n";
  cout <<"\t opfreq_time=D   : Output data every Dth time unit (if op_criterion=1).\n";
  cout <<"\t finishtime=D    : set time to finish simulation, in code time units.\n";
  cout <<"\t optype=S        : Specify type of output file,";
  cout <<             " [1,text]=TEXT,[2,fits]=FITS,[4,both]=FITS+TEXT,[5,silo]=SILO,[6]=SILO+TEXT.\n";
  cout <<"\t outfile=NAME    : Replacement output filename, with path.\n";

  cout <<"\n*********** PHYSICS/Grid OPTIONS *************\n";
  cout <<"\t ooa=N         : modify order of accuracy (either 1 or 2).\n";
  cout <<"\t AVtype=N      : modify type of artificial viscosity:";
  cout <<" 0=none, 1=Falle,Komissarov,Joarder(1998), 2=Colella+Woodward(1984), 3=Sanders et al.(1998)[H-correction].\n";
  cout <<"\t\t\t WARNING -- AVtype=2 IS NOT WORKING WELL.  ONLY USE FKJ98/H-corr.";
  cout <<"\t EtaVisc=D     : modify viscosity parameter to the given double precision value.\n";

  //cout <<"\t coordsys=NAME : override coordinate system to [cartesian,cylindrical]. !DANGEROUS!\n";
  cout <<"\t cfl=D         : change the CFL no. for the simulation, in range (0,1).\n";
  cout <<"\t cooling=N     : cooling=0 for no cooling, 1 for Sutherland&Dopita1993.\n";
  cout <<"\t\t\t For other cooling functions see cooling.cc/cooling.h.\n";

  cout <<"\t solver=N      :\n";
  cout <<"\t\t 0 = Lax-Friedrichs Flux, TESTING ONLY!\n";
  cout <<"\t\t 1 = Linear Riemann Solver : EULER/MHD (mean value average)";
  cout  <<" (Falle, Komissarov, Joarder, 1998),\n";
  cout <<"\t\t 2 = Exact Riemann Solver  : EULER ONLY (Hirsch (199X), Toro, 1999)\n";
  cout <<"\t\t 3 = Hybrid Riemann Solver (1+2)         : EULER ONLY \n";
  cout <<"\t\t 4 = Roe Conserved Variables flux solver : EULER/MHD";
  cout  <<" (e.g. Toro, 1999, Stone, Gardiner et al. 2008)\n";
  cout <<"\t\t 5 = Roe Primitive Variables flux solver : EULER ONLY";
  cout  <<" (e.g. Stone, Gardiner et al. 2008)\n";
  cout <<"\t\t 6 = Flux vector splitting : EULER ONLY (van Leer, 1982) \n";
  cout <<"\t\t 7 = HLLD solver : MHD only (REF) \n";

  cout <<"\t timestep_limit=N:\n";
  cout <<"\t\t 0 = only dynamical Courant condition.\n";
  cout <<"\t\t 1 = dynamical + cooling time limits\n";
  cout <<"\t\t 2 = dyn +cool +recombination time limits\n";
  cout <<"\t\t 3 = dyn +cool +recomb +ionisation time limits\n";
  cout <<"\t\t 4 = dyn +recomb (NO cool, NO ion)\n";

  cout <<"\n*********** PARALLEL CODE ONLY *************\n";
  cout <<"\t maxwalltime=D : change the max. runtime to D in seconds.\n";
  cout <<"\n";
  cout <<"\n*********** NESTED GRID CODE ONLY *************\n";
  cout <<"\t nlevels=N     : modify number of levels in nested grid.";
  cout <<"\n";
  cout <<"********* DEPRECATED -- STILL HERE FOR LEGACY SCRIPTS... ******\n";
  cout <<"\t artvisc=D : modify artificial viscosity, 0=none, Otherwise FalleAV with eta=D,\n";
  cout <<"\t noise=D   : add noise to initial conditions if desired, at fractional level of D.\n";
  cout <<"     *********************************************\n\n";
  return;
}



// ##################################################################
// ##################################################################


int sim_control::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
#ifdef TESTING
  cout <<"(sim_control::Init) Initialising grid"<<"\n";
#endif
  int err=0;
  
#ifdef SERIAL
  SimPM.typeofip=typeOfFile;
  setup_dataio_class(typeOfFile);
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("(INIT::get_parameters) err!=0 Something went wrong",0,err);
#endif // SERIAL

  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  rep.errorTest("(INIT::override_params) err!=0 Something went wrong",0,err);
  
  //
  // Set up the Xmin/Xmax/Range/dx of each level in the nested grid
  //
  grid.resize(1);
  RT.resize(1);
  CI.set_dx(SimPM.dx);

  // Now set up the grid structure.
  cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  err = setup_grid(&(grid[0]),SimPM,&mpiPM);
  cout <<"Init:  &grid="<< &(grid[0])<<", and grid="<< grid[0] <<"\n";
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(INIT::setup_grid) err!=0 Something went wrong",0,err);

  //
  // All grid parameters are now set, so I can set up the appropriate
  // equations/solver class.
  //
  err = set_equations();
  rep.errorTest("(INIT::set_equations) err!=0 Fix me!",0,err);
  spatial_solver->set_dx(SimPM.dx);
  spatial_solver->SetEOS(SimPM.gamma);

  //
  // Now setup Microphysics, if needed.
  //
  err = setup_microphysics(SimPM);
  rep.errorTest("(INIT::setup_microphysics) err!=0",0,err);
  
  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(INIT::assign_initial_data) err!=0 Something went wrong",0,err);

  //
  // Set Ph[] = P[], and then implement the boundary conditions.
  //
  cell *c = grid[0]->FirstPt();
  do {
    for(int v=0;v<SimPM.nvar;v++) c->Ph[v]=c->P[v];
  } while ((c=grid[0]->NextPt(c))!=0);

  //
  // If I'm using the GLM method, make sure Psi variable is
  // initialised to zero.
  //
  if (SimPM.eqntype==EQGLM && SimPM.timestep==0) {
#ifdef TESTING
    cout <<"Initial state, zero-ing glm variable.\n";
#endif
    c = grid[0]->FirstPt(); do {
      c->P[SI] = c->Ph[SI] = 0.;
    } while ( (c=grid[0]->NextPt(c)) !=0);
  }

  //
  // Assign boundary conditions to boundary points.
  //
  err = boundary_conditions(SimPM, grid[0]);
  rep.errorTest("(INIT::boundary_conditions) err!=0",0,err);
  err = assign_boundary_data(SimPM, grid[0]);
  rep.errorTest("(INIT::assign_boundary_data) err!=0",0,err);

  //
  // Setup Raytracing on each grid, if needed.
  //
  err += setup_raytracing(SimPM, grid[0], &(RT[0]));
  err += setup_evolving_RT_sources(SimPM, RT[0]);
  rep.errorTest("Failed to setup raytracer and/or microphysics",0,err);

  //
  // If testing the code, this calculates the momentum and energy on the domain.
  //
  initial_conserved_quantities(grid[0]);

  err += TimeUpdateInternalBCs(SimPM, grid[0], SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  err += TimeUpdateExternalBCs(SimPM, grid[0], SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  if (err) 
    rep.error("first_order_update: error from bounday update",err);



  //
  // If using opfreq_time, set the next output time correctly.
  //
  if (SimPM.op_criterion==1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!",SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
    double tmp = 
      ((SimPM.simtime/SimPM.opfreq_time)-
       floor(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
    SimPM.next_optime-= tmp;
  }

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    setup_dataio_class(SimPM.typeofop);
    if (!dataio) rep.error("INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

#ifdef SERIAL
  if (SimPM.timestep==0) {
    cout << "(INIT) Writing initial data.\n";
    err=output_data(grid);
    if (err)
      rep.error("Failed to write file!","maybe dir does not exist?");
  }
  cout <<"------------------------------------------------------------\n";
#endif // SERIAL
  
#ifdef TESTING
  c = (grid[0])->FirstPt_All();
  do {
    if (pconst.equalD(c->P[RO],0.0)) {
      cout <<"zero data in cell: ";
      CI.print_cell(c);
    }
  } while ( (c=(grid[0])->NextPt_All(c)) !=0 );
#endif // TESTING
  
  return(0);
}



// ##################################################################
// ##################################################################



void sim_control::setup_dataio_class(
      const int typeOfFile ///< type of I/O: 1=text,2=fits,5=silo
      )
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

  case 1: // Start From ASCII Parameterfile.
    dataio = new dataio_text(SimPM);
    rep.error("dataio_text initialisation",dataio);
    break;

#ifdef FITS
  case 2: // Start from FITS restartfile.
  case 3: // Fits restartfile in table format (slower I/O than image...)
    dataio = new DataIOFits(SimPM);
    break;
  case 4: // fits +ascii
    dataio = new DataIOFits(SimPM);
    textio = new dataio_text(SimPM);
    break;
#endif // if FITS

#ifdef SILO
  case 5: // Start from Silo snapshot.
    dataio = new dataio_silo (SimPM, "DOUBLE");
    break; 
  case 6: // silo + text
    dataio = new dataio_silo (SimPM, "DOUBLE");
    textio = new dataio_text (SimPM);
    if (!textio) rep.error("INIT:: textio initialisation",SimPM.typeofop);
    break;
#endif // if SILO
  default:
    rep.error("sim_control::Init unhandled filetype",typeOfFile);
  }
  return;
}


// ##################################################################
// ##################################################################



int sim_control::override_params(int narg, string *args)
{
  cout <<"------------------------------------------------------\n";
  cout <<"--------  Overriding parameters if requested ---------\n";
  // Find command line params and assign them to SimPM.whatever.
  for (int i=2;i<narg;i++) {

    if      (args[i].find("ooa=") != string::npos) {
      // Assign order of accuracy;  string is 'ooa=N', where N=1 or 2.
      int tmp=SimPM.spOOA;
      SimPM.spOOA = atoi((args[i].substr(4)).c_str());
      SimPM.tmOOA = SimPM.spOOA;
      cout <<"\tOVERRIDE PARAMS: Resetting OOA from ooa="<<tmp;
      cout <<" to command-line value = "<<SimPM.spOOA<<"\n";
    }

    else if (args[i].find("nlevels=") != string::npos) {
      // Assign number of grid levels;  string is 'nlevels=N', where N>=1.
      int tmp=SimPM.grid_nlevels;
      SimPM.grid_nlevels = atoi((args[i].substr(8)).c_str());
      cout <<"\tOVERRIDE PARAMS: Resetting nlevels from nlevels="<<tmp;
      cout <<" to command-line value = "<<SimPM.grid_nlevels<<"\n";
    }

    else if (args[i].find("AVtype=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: old AV="<<SimPM.artviscosity<<" ... overriding!\n";
      // Assign art.viscosity parameter. String is 'artvisc=I' with I in [0,N].
      int v = atoi((args[i].substr(7)).c_str());
      if      (v == 0) {
        cout <<"\t\tNot using artificial viscosity.\n";
        SimPM.artviscosity=0; SimPM.etav=0.;
      }
      else if (v == 1) {
        cout <<"\t\tUsing Falle, Komissarov, Joarder (1998) AV prescription.\n";
        SimPM.artviscosity = AV_FKJ98_1D; // ==1
        SimPM.etav = 0.1;
      }
      else if (v == 2) {
        cout <<"\t\tUsing Colella and Woodward (1984) AV prescription (Lapidus).\n";
        cout <<"\t\t****** WARNING, THIS NEEDS TESTING, EXPERIMENTAL CODE!! ****\n";
        SimPM.artviscosity=AV_LAPIDUS; // ==2 (NEEDS TESTING!!!)
        SimPM.etav = 0.1;
      }
      else if (v == 3) {
        cout <<"\t\tUsing the H-correction of Sanders et al. (1998,JCP,145,511).\n";
        cout <<"\t\t****** WARNING, THIS NEEDS TESTING, EXPERIMENTAL CODE!! ****\n";
        SimPM.artviscosity=AV_HCORRECTION;
        SimPM.etav = 0.1; // This parameter is redundant for the H-correction.
      }
      else if (v == 4) {
        cout <<"\t\tUsing the H-correction of Sanders et al. (1998,JCP,145,511)\n";
        cout <<"\t\twith the 1D viscosity of Falle, Komissarov, Joarder (1998)\n";
        cout <<"\t\t****** WARNING, THIS NEEDS TESTING, EXPERIMENTAL CODE!! ****\n";
        SimPM.artviscosity=AV_HCORR_FKJ98; // ==4 (NEEDS TESTING!!!)
        SimPM.etav = 0.1;
      }
      else if (v == AV_VonNeuRicht) {
        // AV_VonNeuRicht==5
        cout <<"\t\tUsing Multi-D von Neumann & Richtmeyer (1950) viscosity.\n";
        cout <<"\t\tSee Tscharnuter & Winkler (1979), Stone & Norman (1992).\n";
        cout <<"\t\tWARNING -- THIS ONLY WORKS WITH EQNTYPE==9(EQEUL_EINT).\n";
        SimPM.artviscosity=AV_VonNeuRicht;
        SimPM.etav = 1.0;
      }
      else {
        cout <<"\t\t********************** FIX ME!!!! **************************\n";
        cout <<"\t\tDIDN'T UNDERSTAND AV="<<v<<", SETTING TO FALLE et al (1998).\n";
        cout <<"\t\t********************** FIX ME!!!! **************************\n";
        SimPM.artviscosity = 1;
        SimPM.etav = 0.1;
        rep.error("Bad viscosity flag from command-line",v);
      }  
      cout <<"\tOVERRIDE PARAMS: setting AV = "<<SimPM.artviscosity<<" and eta = "<<SimPM.etav<<"\n"; 
    }

    else if (args[i].find("EtaVisc=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: old and eta="<<SimPM.etav<<" ... overriding!\n";
      // Assign art.viscosity parameter. String is 'artvisc=D' with D in [0,N].
      double visc = atof((args[i].substr(8)).c_str());
      cout <<"\tOVERRIDE PARAMS: Resetting eta_visc from ";
      cout <<SimPM.etav<<" to "<<visc<<"\n";
      if (visc<0.0 || !isfinite(visc))
        rep.error("Error: eta viscosity parameter outside allowed range!",visc);
      SimPM.etav = visc;
    }

    else if (args[i].find("artvisc=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: old AV="<<SimPM.artviscosity<<" and eta="<<SimPM.etav<<" ... overriding!\n";
      // Assign art.viscosity parameter. String is 'artvisc=D' with D in [0,N].
      double visc = atof((args[i].substr(8)).c_str());
      if(fabs(visc) <= 1.e-6) {
        cout <<"\t\tNot using artificial viscosity.\n";
        SimPM.artviscosity=0; SimPM.etav=0.;
      }
      else if (visc <=0) {
        SimPM.artviscosity = 1;
        SimPM.etav = 0.15;
      }
      else {
        SimPM.artviscosity=1;
        SimPM.etav = visc;
      }
      cout <<"\tOVERRIDE PARAMS: setting AV = "<<SimPM.artviscosity<<" and eta = "<<SimPM.etav<<"\n"; 
    }

    else if (args[i].find("opfreq=") != string::npos) {
      // Assign output frequency to new value. String is 'opfreq=N' with N=[0..Nmax].
      int tmp = atoi((args[i].substr(7)).c_str());
      cout <<"\tOVERRIDE PARAMS: Resetting opfreq from "<<SimPM.opfreq<<" to "<<tmp<<"\n";
      SimPM.opfreq = tmp;
    }
    else if (args[i].find("outfile=") != string::npos) {
      // Assign a new output file.  string is outfile=char[128max]
      string tmp = args[i].substr(8);
      cout <<"\tOVERRIDE PARAMS: Resetting output file from "<<SimPM.outFileBase<<"\n";
      cout <<"\t\t\tto new name: "<<tmp<<".xxx.ftype\n";
      SimPM.outFileBase = tmp; tmp.clear();
    }
    
    else if (args[i].find("optype=") != string::npos) {
      // assign new op-type; 1=TXT,2=FITS,3=FitsTable,4=TXT+FITS, 5=SILO
      string now;
      if      (SimPM.typeofop==1) now="text";
      else if (SimPM.typeofop==2) now="fits";
      else if (SimPM.typeofop==3) now="ftab";
      else if (SimPM.typeofop==4) now="f+tt";
      else if (SimPM.typeofop==5) now="silo";
      else if (SimPM.typeofop==6) now="silo+text";
      else rep.error("What kind of output is this?",SimPM.typeofop);

      string chg = args[i].substr(7);
      int tmp=-1;
      if      (chg=="text" || chg=="txt" || chg=="TEXT" || chg=="TXT" || chg=="1")
        tmp=1;
      else if (chg=="fits" || chg=="FITS" || chg=="2")
        tmp=2;
      else if (chg=="ftab" || chg=="FTAB" || chg=="3")
        tmp=3;
      else if (chg=="txtfits" || chg=="both" || chg=="BOTH" || chg=="4")
        tmp=4;
      else if (chg=="silo" || chg=="SILO" || chg=="5")
        tmp=5;
      else if (chg=="txtsilo" || (chg=="textsilo") || (chg=="6"))
        tmp=6;
      else rep.error("What kind of output do you want?",chg);
      cout <<"\tOVERRIDE PARAMS: Resetting output file type from "<<now<<" to "<<chg<<"\n";
      SimPM.typeofop = tmp;
    }
    
    else if (args[i].find("noise=") != string::npos) {
      // assign new value to addnoise, double frac = fractional noise level
      double tmp = atof((args[i].substr(6)).c_str());
      cout <<"\tOVERRIDE PARAMS: Resetting addnoise value from "<<tmp<<" to "<<SimPM.addnoise<<"\n";
      SimPM.addnoise = tmp;
    }
    
    else if (args[i].find("finishtime=") != string::npos) {
      // assign new value to finishtime.
      double tmp = atof((args[i].substr(11)).c_str());
      cout <<"\tOVERRIDE PARAMS: Resetting finishtime value from "<<SimPM.finishtime<<" to "<<tmp<<"\n";
      SimPM.finishtime = tmp;
    }
    
    else if (args[i].find("redirect=") != string::npos) {
      // this is already handled by main
      cout <<"\tOVERRIDE PARAMS: already redirecting stdout, continuing...\n";
    }

    else if (args[i].find("maxwalltime=") != string::npos) {
      /** \section walltime Max. Walltime
       * Max. Walltime is only used by parallel code for now, so it is
       * read in mainMPI.cc, and resets the value of mpiPM.max_walltime,
       * but it is ignored by the serial code.  I might want to change
       * this at some time in the future, because there is no reason why
       * the walltime limit should be only a parallel code feature.
       */
      // this is already handled by mainMPI.cc, and ignored for serial code.
    }
    
    else if (args[i].find("coordsys=") != string::npos) {
      // Change the coordinate system!
      cout <<"\tOVERRIDE PARAMS: Resetting the coordinate system from value="<<SimPM.coord_sys;
      string t = args[i].substr(9);
      if (t=="cartesian" || t=="cart" || t=="crt") {
  cout <<" to cartesian coords.\n";
  SimPM.coord_sys = COORD_CRT;
      }
      else if (t=="cylindrical" || t=="cyl") {
  cout <<" to cylindrical coords.\n";
  SimPM.coord_sys = COORD_CYL;
      }
      else if (t=="spherical"   || t=="sph") {
  cout <<" to spherical coords.\n";
  SimPM.coord_sys = COORD_SPH;
      }
      else rep.error("don't know this coordinate system",t);
      cout <<"\t\t\tTHIS IS DANGEROUS!!! PROBABLY FATAL.\n";
    }
    
    else if (args[i].find("cfl=") != string::npos) {
      // Assign cfl no., where string is 'cfl=0.X'
      cout <<"\tOVERRIDE PARAMS: Resetting CFL from original value of "<<SimPM.CFL;
      double c = atof((args[i].substr(4)).c_str());
      if (c<0. || !isfinite(c)) rep.error("Bad CFL no.",c);
      else if (c>1.) cout <<"\tWARNING: CFL no. >1, so results will be unstable!!!\n";
      SimPM.CFL = c;
      cout <<" to command-line value = "<<SimPM.CFL<<"\n";
    }
    
    else if (args[i].find("gamma=") != string::npos) {
      // Assign new value to EOS gamma, where string is 'gamma=X.XXXXX'
      cout <<"\tOVERRIDE PARAMS: Resetting EOS gamma from original value of "<<SimPM.gamma;
      double c = atof((args[i].substr(6)).c_str());
      if (c<=1.0 || !isfinite(c)) rep.error("Bad EOS gamma no.",c);
      else if (c>2.0) cout <<"\tWARNING: gamma >2 ?\n";
      SimPM.gamma = c;
      cout <<" to command-line value = "<<SimPM.gamma<<"\n";
    }

    else if (args[i].find("cooling=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting cooling";
      int c = atoi((args[i].substr(8)).c_str());
      if (c<0 || c>100) rep.error("Bad cooling flag (only 0-11 allowed",c);
      cout <<" flag from "<<SimPM.EP.cooling;
      SimPM.EP.cooling = c;
      cout <<" to "<<SimPM.EP.cooling<<"\n";
    }
    
    else if (args[i].find("dynamics=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting dynamics";
      int c = atoi((args[i].substr(9)).c_str());
      if (c<0 || c>1 ) rep.error("Bad dynamics flag (only 0,1, allowed",c);
      cout <<" flag from "<<SimPM.EP.dynamics;
      SimPM.EP.dynamics = c;
      cout <<" to "<<SimPM.EP.dynamics<<"\n";
    }
    
    else if (args[i].find("solver=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting solver: ";
      cout <<"0=LF,1=RSlin,2=RSexact,3=RShybrid,4=RSroe,5=RSroePV,6=FVS:";
      int c = atoi((args[i].substr(7)).c_str());
      if (c<0 || c>7 )
  rep.error("Bad solver flag (only 0,1,2,3,4,5,6 allowed",c);
      cout <<" solver from "<<SimPM.solverType;
      SimPM.solverType = c;
      cout <<" to "<<SimPM.solverType<<"\n";
    }
    
    else if (args[i].find("op_criterion=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting op_criterion from "<<SimPM.op_criterion<<" to ";
      int c = atoi((args[i].substr(13)).c_str());
      if (c<0 || c>1) rep.error("Bad op_criterion flag:",c);
      SimPM.op_criterion = c;
      cout << SimPM.op_criterion <<"\n";
    }
    else if (args[i].find("opfreq_time=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting opfreq_time from ";
      cout <<SimPM.opfreq_time<<" units [NOT YEARS!!!] to ";
      double c = atof((args[i].substr(12)).c_str());
      if (c<0.0 || c>1.e50) rep.error("Bad opfreq_time flag:",c);
      SimPM.opfreq_time = c;
      cout << SimPM.opfreq_time<<" units." <<"\n";
      SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
      if (SimPM.timestep>0) {
        double tmp = ((SimPM.simtime/SimPM.opfreq_time)-static_cast<int>(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
        SimPM.next_optime-= tmp;
      }
    }

    else if (args[i].find("min_timestep=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: resetting min_timestep from ";
      cout <<SimPM.min_timestep<<" units [NOT YEARS!!!] to ";
      double c = atof((args[i].substr(13)).c_str());
      if (c<0.0 || c>1.e50) rep.error("Bad min_timestep flag:",c);
      SimPM.min_timestep = c;
      cout << SimPM.min_timestep<<" units." <<"\n";
    }

    else if (args[i].find("limit_timestep=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: limiting timestep: changing from ";
      cout <<SimPM.EP.MP_timestep_limit <<" to ";
      int c = atoi((args[i].substr(15)).c_str());
      if (c<0 || c>5) rep.error("Bad MP_timestep_limit flag:",c);
      SimPM.EP.MP_timestep_limit = c;
      cout << SimPM.EP.MP_timestep_limit <<"\n";
    }

    else if (args[i].find("checkpt_freq=") != string::npos) {
      cout <<"\tOVERRIDE PARAMS: checkpointing freq.: changing from ";
      cout <<SimPM.checkpoint_freq <<" to ";
      int c = atoi((args[i].substr(13)).c_str());
      if (c<0) rep.error("Bad checkpoint_freq flag:",c);
      SimPM.checkpoint_freq = c;
      cout << SimPM.checkpoint_freq <<"\n";
    }


    else rep.error("Don't recognise this optional argument, please fix.",args[i]);
  }
  cout <<"------------------------------------------------------\n\n";
  return(0);
} // sim_control::override_params



// ##################################################################
// ##################################################################



int sim_control::set_equations()
{
  cout <<"------------------------------------------------------\n";
  cout <<"--------  Setting up solver for eqauations------------\n";

  if (SimPM.solverType <0) rep.error("set_equations: solverType not set yet.",SimPM.solverType);
  if (SimPM.eqntype <=0) rep.error("set_equations: eqntype not set yet.",SimPM.eqntype);
  if (SimPM.eqnNDim !=3) rep.error("set_equations: eqnNDim not set right.",SimPM.eqnNDim);
  if (SimPM.nvar <=0) rep.error("set_equations: nvar not set yet.",SimPM.nvar);
  if (SimPM.ndim <=0) rep.error("set_equations: ndim not set yet.",SimPM.ndim);
  if (SimPM.artviscosity <0) rep.error("set_equations: artviscosity not set yet.",SimPM.artviscosity);
  if (SimPM.coord_sys <0) rep.error("set_equations: coordinate system not set.",SimPM.coord_sys);
  
  // Set the mean state vector to unity for now (passed into Riemann Solver, if needed.)
  //  double *meanvec = new double [SimPM.nvar];
  // for (int v=0;v<SimPM.nvar;v++) meanvec[v] = 1.;

  if (spatial_solver) spatial_solver=0;

  if (SimPM.coord_sys==COORD_CRT) {
    cout <<"set_equations() Using Cartesian coord. system.\n";
    switch (SimPM.eqntype) {
    case EQEUL:
      cout <<"set_equations() Using Euler Equations.\n";
      spatial_solver = new class FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
    case EQMHD:
      cout <<"set_equations() Using Ideal MHD Equations.\n";
      spatial_solver = new class FV_solver_mhd_ideal_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQMHD);
      break;
    case EQGLM:
      cout <<"set_equations() Using GLM MHD Equations.\n";
      spatial_solver = new class FV_solver_mhd_mixedGLM_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQGLM);
      break;
    case EQFCD:
      cout <<"set_equations() Using Field-CD MHD Equations.\n";
      rep.error("Field CD got lost in some code updates -- can re-do it from svn rev.100 or so",EQFCD);
      break;
    default:
      rep.error("Don't know the specified equations...",SimPM.eqntype);
      break;
    }
  } // cartesian

  else if (SimPM.coord_sys==COORD_CYL) {
    cout <<"set_equations() Using Cylindrical coord. system.\n";
    switch (SimPM.eqntype) {
    case EQEUL:
      cout <<"set_equations() Using Euler Equations.\n";
      spatial_solver = new class cyl_FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
#ifdef INCLUDE_EINT_ADI_HYDRO
    case EQEUL_EINT:
      cout <<"set_equations() Using Euler Equations, integrating ***INTERNAL ENERGY***.\n";
      spatial_solver = new class cyl_FV_solver_Hydro_Euler_Eint(SimPM.nvar,
            SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma,
            SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQEUL_EINT);
      break;      
#endif // if INCLUDE_EINT_ADI_HYDRO
    case EQMHD:
      cout <<"set_equations() Using Ideal MHD Equations.\n";
      spatial_solver = new class cyl_FV_solver_mhd_ideal_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQMHD);
      break;
    case EQGLM:
      cout <<"set_equations() Using GLM MHD Equations.\n";
      spatial_solver = new class cyl_FV_solver_mhd_mixedGLM_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQGLM);
      break;
    default:
      rep.error("not implemented yet for axisymmetry (only adiabatic hydro!)!!!!!!!!!!!",SimPM.eqntype);
    }
  } // axisymmetric
  
  else if (SimPM.coord_sys==COORD_SPH) {
    cout <<"set_equations() Using Spherical coordinate system.\n";
    switch (SimPM.eqntype) {
    case EQEUL:
      cout <<"set_equations() Using Euler Equations.\n";
      spatial_solver = new class sph_FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL,
            SimPM.dx, SimPM.gamma, SimPM.RefVec,
            SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
#ifdef INCLUDE_EINT_ADI_HYDRO
    case EQEUL_EINT:
      cout <<"set_equations() Using Euler Equations, integrating ***INTERNAL ENERGY***.\n";
      spatial_solver = new class sph_FV_solver_Hydro_Euler_Eint(SimPM.nvar,
            SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma,
            SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!spatial_solver) rep.error("Couldn't set up solver/equations class.",EQEUL_EINT);
      break;      
#endif // if INCLUDE_EINT_ADI_HYDRO
     default:
      rep.error("not implemented yet for spherical (only hydro in 1D!)",SimPM.eqntype);
    }
  } // axisymmetric


  //
  // Check that we set up an equations class!
  //
  if (!spatial_solver)
    rep.error("sim_control::set_equations() Failed to initialise an equations class.",spatial_solver);

  cout <<"------------------------------------------------------\n\n";
  return(0);
} // set_equations.



// ##################################################################
// ##################################################################




int sim_control::output_data(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  ///
  /// \section freq Output Frequency.
  /// I'm going to set it to output every nth timestep.  Can modify this if needed.
  /// This is a parameter in the parameterfile.  If opfreq is set to a negative
  /// number then only output the final state.
  /// We can also set it so that it outputs every D time units, and there is
  /// automatic checkpointing to filebase.9999999.[txt/fits/silo] every 100 steps.
  ///
  int err=0;
  //
  // First see if we want to write a checkpoint file, and if so then do it!
  //
  int checkpoint_freq;
  if (SimPM.checkpoint_freq>0) {
    checkpoint_freq = SimPM.checkpoint_freq;
  }
  else {
   checkpoint_freq = 250;
  }
  long int checkpoint_id=99999999;
  if ((SimPM.timestep !=0) && ((SimPM.timestep % checkpoint_freq)==0)) {
    cout <<"Checkpointing...\n";
    if ((SimPM.timestep % (2*checkpoint_freq))==0)
      checkpoint_id=99999998;
    else
      checkpoint_id=99999999;
    
    err = dataio->OutputData(SimPM.outFileBase, grid, SimPM, checkpoint_id);
    if (textio) err += textio->OutputData(SimPM.outFileBase, grid, SimPM, checkpoint_id);
    if (err) {cerr<<"\t Error writing data for checkpointing.\n"; return(1);}
  }
    
  //
  // Now we move on the various output criteria.  The if statements
  // call return(0) if we decide we don't need to ouput, so if we get
  // past them then just output the current timestep.
  //

  //
  // Always output at the start
  //
  if (SimPM.timestep==0) {}
  //
  // If outputting every nth step, see if we are on an output step.
  // If not, return.
  //
  else if (SimPM.op_criterion==0) {
    if( (SimPM.opfreq==0) && (SimPM.maxtime==false) &&(SimPM.timestep!=0) ) return 0;
    // Next check if we are in an outputting timestep.
    else if( (SimPM.maxtime==false) && (SimPM.timestep%SimPM.opfreq)!=0 &&(SimPM.timestep!=0) ) return(0);
  }
  //
  // If outputting every D time units, see if we are at an output time, and if not
  // then return.
  //
  else if (SimPM.op_criterion==1) {
    if (!pconst.equalD(SimPM.simtime,SimPM.next_optime) && (SimPM.maxtime==false))
      return 0;
    else 
      // we are to output data, so advance the 'next' counter and continue.
      SimPM.next_optime += SimPM.opfreq_time;
  }
  else rep.error("op_criterion must be 0 or 1",SimPM.op_criterion);

  //
  // Since we got past all that, we are in a timestep that should be outputted, so go and do it...
  //
  cout <<"\tSaving data, at simtime: "<<SimPM.simtime << " to file "<<SimPM.outFileBase<<"\n";
  err = dataio->OutputData(SimPM.outFileBase, grid, SimPM, SimPM.timestep);
  if (textio) err += textio->OutputData(SimPM.outFileBase, grid, SimPM, SimPM.timestep);
  if (err) {cerr<<"\t Error writing data.\n"; return(1);}
  return(0);
}


// ##################################################################
// ##################################################################



int sim_control::initial_conserved_quantities(
      class GridBaseClass *grid
      )
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TESTING 
  // Only track the totals if I am testing the code.
  pion_flt u[SimPM.nvar];
  dp.initERG = 0.;  dp.initMMX = dp.initMMY = dp.initMMZ = 0.;
  dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange = 0.0;
  //  cout <<"initERG: "<<dp.initERG<<"\n";
  class cell *cpt=grid->FirstPt();
  do {
     spatial_solver->PtoU(cpt->P,u,SimPM.gamma);
     dp.initERG += u[ERG]*spatial_solver->CellVolume(cpt);
     dp.initMMX += u[MMX]*spatial_solver->CellVolume(cpt);
     dp.initMMY += u[MMY]*spatial_solver->CellVolume(cpt);
     dp.initMMZ += u[MMZ]*spatial_solver->CellVolume(cpt);
  } while ( (cpt = grid->NextPt(cpt)) !=0);
  //cout <<"!!!!! cellvol="<<spatial_solver->CellVolume(cpt)<< "\n";
  cout <<"(LFMethod::InitialconservedQuantities) Total Energy = "<< dp.initERG <<"\n";
  cout <<"(LFMethod::InitialconservedQuantities) Total x-Momentum = "<< dp.initMMX <<"\n";
  cout <<"(LFMethod::InitialconservedQuantities) Total y-Momentum = "<< dp.initMMY <<"\n";
  cout <<"(LFMethod::InitialconservedQuantities) Total z-Momentum = "<< dp.initMMZ <<"\n";
#endif //TESTING
  return(0);
} //initial_conserved_quantities()



// ##################################################################
// ##############     SIMULATION CONTROL FUNCTIONS     ##############
// ##################################################################


/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int sim_control::Time_Int(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control::Time_Int) STARTING TIME INTEGRATION."<<"\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  clk.start_timer("Time_Int"); double tsf=0;
  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid[0]);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid[0]);
#endif
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources(SimPM,RT[0]);
    if (err) {
      cerr <<"(TIME_INT::update_evolving_RT_sources()) something went wrong!\n";
      return err;
    }

    //clk.start_timer("advance_time");
    err+= advance_time(grid[0],RT[0]);
    //cout <<"advance_time took "<<clk.stop_timer("advance_time")<<" secs.\n";
    if (err!=0){cerr<<"(TIME_INT::advance_time) err!=0 Something went wrong\n";return(1);}

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
    tsf=clk.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err+= output_data(grid);
    if (err!=0){cerr<<"(TIME_INT::output_data) err!=0 Something went wrong\n";return(1);}

    err+= check_eosim();
    if (err!=0){cerr<<"(TIME_INT::) err!=0 Something went wrong\n";return(1);}
  }

  cout <<"(sim_control::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

  tsf=clk.time_so_far("Time_Int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep<<" wall-time=";
  cout <<tsf<<" time/step="<<tsf/static_cast<double>(SimPM.timestep)<<"\n";
  cout <<"STEPS "<<SimPM.timestep;
  cout.setf( ios_base::scientific );
  cout.precision(6);
  cout <<"\t"<<tsf<<"\t"<<tsf/static_cast<double>(SimPM.timestep);
  cout <<"\t"<<static_cast<double>(SimPM.timestep*SimPM.Ncell)/tsf<<"\n";
  cout <<"------------------------------------------------------------\n";

  return(0);
}



// ##################################################################
// ##################################################################



#ifdef CHECK_MAGP
///
/// This is only for a test problem -- it checks the magnetic
/// pressure on the full domain and outputs it to screen
///
void sim_control::calculate_magnetic_pressure(
      class GridBaseClass *grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the total magnetic pressure on the domain, normalised to the
  // initial value.
  //
  double magp=0.0, cellvol=0.0;
  static double init_magp=-1.0;
    
  cell *c=grid->FirstPt();
  do {
    if (!c->isbd) 
      magp += (spatial_solver->Ptot(c->P,0.0) - c->P[PG]) *
                spatial_solver->CellVolume(c);
  } while ( (c =grid->NextPt(c)) !=0);
  if (init_magp<0) init_magp = magp;
  cout <<SimPM.simtime<<"\t"<<magp/init_magp<<"\t"<<magp<<"\n";
  return;
}
#endif // CHECK_MAGP



// ##################################################################
// ##################################################################



#ifdef BLAST_WAVE_CHECK
///
/// If running a 1D spherical blast wave, calculate the shock position
/// and output to screen.
///
void sim_control::calculate_blastwave_radius(
      class GridBaseClass *grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the blast wave outer shock position.
  //
  double shockpos=0.0;
  static double old_pos=0.0;
  //bool shock_found = false;
  //  static double last_dt=0.0;

  //if (shock_found) continue;
  cell *c=grid->LastPt();
  if (fabs(c->P[VX])>=1.0e4) {
    cout<<"grid "<<l<<" does not contain shock.\n";
    shockpos = CI.get_dpos(c,Rsph);
  }
  else {
    do {
      c = grid->NextPt(c,RNsph);
      //cout <<c->id<<", vx="<<c->P[VX]<<"\n";
    } while ( c!=0 && fabs(c->P[VX])<1.0e4);
    if (c && fabs(c->P[VX] >= 1.0e4)) {
      shockpos = CI.get_dpos(c,Rsph);
      //shock_found=true;
    }
  }
  
  if (pconst.equalD(old_pos,0.0))
    old_pos = shockpos;
  cout <<SimPM.simtime<<"\t"<<shockpos;
  //cout <<"\t"<<(shockpos-old_pos)/(SimPM.dt+TINYVALUE);
  cout <<"\n";
  old_pos=shockpos;
  return;
}
#endif // BLAST_WAVE_CHECK



// ##################################################################
// ##################################################################




int sim_control::calc_timestep(
      class GridBaseClass *grid, ///< pointer to grid.
      class RayTracingBase *raytracer ///< raytracer for this grid.
      )
{
#ifdef TESTING
  cout <<"sim_control::calc_timestep(): g="<<grid<<", rt="<<raytracer<<"\n";
#endif
  //
  // This is a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt(grid);
  t_mp  = calc_microphysics_dt(grid,raytracer);
#ifdef TESTING
  if (t_mp<t_dyn)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";
#endif
  SimPM.dt = min(t_dyn,t_mp);

#ifdef THERMAL_CONDUCTION
  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in spatial_solver->preprocess_data()
  //
  double t_cond = calc_conduction_dt_and_Edot();
#ifdef TESTING
  if (t_cond<t_dyn && t_cond<t_mp) {
    cout <<"CONDUCTION IS LIMITING TIMESTEP: t_c="<<t_cond<<", t_m="<<t_mp;
    cout <<", t_dyn="<<t_dyn<<"\n";
  }
#endif
  SimPM.dt = min(SimPM.dt, t_cond);
#endif // THERMAL CONDUCTION

  //
  // If using MHD with GLM divB cleaning, the following sets the hyperbolic wavespeed.
  // If not, it does nothing.  By setting it here and using t_dyn, we ensure that the
  // hyperbolic wavespeed is equal to the maximum signal speed on the grid, and not
  // an artificially larger speed associated with a shortened timestep.
  //
  spatial_solver->GotTimestep(t_dyn);

  //
  // Check that the timestep doesn't increase too much between steps, and that it 
  // won't bring us past the next output time or the end of the simulation.
  // This function operates on SimPM.dt, resetting it to a smaller value if needed.
  //
  timestep_checking_and_limiting();

  // sets the timestep info in the solver class.
  spatial_solver->Setdt(SimPM.dt);

#ifdef TESTING
  cout <<"sim_control::calc_timestep() finished.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int sim_control::check_eosim()
{
  //  cout <<"Checking eosim.";
  //  cout <<"finishtime="<<SimPM.finishtime<<"\n";

  if (SimPM.finishtime >0) {
    if (SimPM.simtime >= SimPM.finishtime ) {
      SimPM.maxtime=true;
      cout <<"finishtime="<<SimPM.finishtime<<"\n";
      return(0);
    }
  }
  else {
    cout <<"finishtime="<<SimPM.finishtime<<"\n";
    rep.error("Don't know how to check for end of simulation.",2);
  }  

  // Diagnose boundary cells... comment out if not needed.
  //  for (int bc=0;bc<SimPM.Nbc;bc++) {
  //    rep.printVec("boundary vec, P :",bd[bc].P );
  //    rep.printVec("boundary vec, Ph:",bd[bc].Ph);
  //  }
  // diagnose bcs
  return(0);
}

// ##################################################################
// ##################################################################



int sim_control::check_energy_cons(
      class GridBaseClass *grid
      )
{
#ifdef TESTING
  // Energy, and Linear Momentum in x-direction.
  pion_flt u[SimPM.nvar];
  double ergNow=0., mmxNow = 0., mmyNow = 0., mmzNow = 0.;
  class cell *cpt = grid->FirstPt();
  do {
     spatial_solver->PtoU(cpt->P,u,SimPM.gamma);
     ergNow += u[ERG]*spatial_solver->CellVolume(cpt);
     mmxNow += u[MMX]*spatial_solver->CellVolume(cpt);
     mmyNow += u[MMY]*spatial_solver->CellVolume(cpt);
     mmzNow += u[MMZ]*spatial_solver->CellVolume(cpt);
  } while ( (cpt =grid->NextPt(cpt)) !=0);
  //cout <<"!!!!! cellvol="<<spatial_solver->CellVolume(cpt)<< "\n";
  double relerror=0.0;
//  cout <<"(LFMethod::check_energy_cons) Initial Monentum: "<<dp.initMMX<<"\n";
//  cout <<"(LFMethod::check_energy_cons) Initial Energy:   "<<dp.initERG<<"\n";
//  cout <<"dp = "<<&dp<<"\n";
  
  relerror = fabs(ergNow-dp.initERG)/(dp.initERG+TINYVALUE);
  if (relerror>1.e5*MACHINEACCURACY) { // && ergNow>2*MACHINEACCURACY) {
    cout <<"(LFMethod::check_energy_cons) Total Energy = "<<ergNow;
    cout <<" and relative error is "<<relerror<<" at timestep "<<SimPM.timestep<<"\n";
    cout <<"(LFMethod::check_energy_cons) accounting says it is "<<dp.initERG<<"\n";
  }

  relerror = fabs(mmxNow-dp.initMMX)/(dp.initMMX+TINYVALUE);
  if (relerror>1.e5*MACHINEACCURACY && fabs(mmxNow)>1.e5*MACHINEACCURACY) {
    cout <<"(LFMethod::check_energy_cons) Total x-Momentum = "<<mmxNow;
    cout <<" and relative error is "<<relerror<<" at timestep "<<SimPM.timestep<<"\n";
    cout <<"(LFMethod::check_energy_cons) accounting says it is "<<dp.initMMX<<"\n";
  }

  relerror = fabs(mmyNow-dp.initMMY)/(dp.initMMY+TINYVALUE);
//  if (relerror>1.e5*MACHINEACCURACY && fabs(mmyNow)>1.e5*MACHINEACCURACY) {
//    cout <<"(LFMethod::check_energy_cons) Total y-Momentum = "<<mmyNow;
//    cout <<" and relative error is "<<relerror<<" at timestep "<<SimPM.timestep<<"\n";
//    cout <<"(LFMethod::check_energy_cons) accounting says it is "<<dp.initMMY<<"\n";
//  }

  relerror = fabs(mmzNow-dp.initMMZ)/(dp.initMMZ+TINYVALUE);
//  if (relerror>1.e5*MACHINEACCURACY && fabs(mmzNow)>1.e5*MACHINEACCURACY) {
//    cout <<"(LFMethod::check_energy_cons) Total z-Momentum = "<<mmzNow;
//    cout <<" and relative error is "<<relerror<<" at timestep "<<SimPM.timestep<<"\n";
//  }

#endif //TESTING
  return(0);
}


// ##################################################################
// ##################################################################


/*****************************************************************/
/*********************** FINISH SIMULATION ***********************/
/*****************************************************************/
int sim_control::Finalise(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  int err=0;
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control::Finalise) FINALISING SIMULATION."<<"\n";
  err += check_energy_cons(grid[0]);
  err+= output_data(grid);
  rep.errorTest("(FINALISE::output_data) Something went wrong",0,err);
  cout <<"\tSimTime = "<<SimPM.simtime<<"   #timesteps = "<<SimPM.timestep<<"\n";
#ifdef TESTING
  cout <<"(sim_control::Finalise) DONE.\n";
#endif
  cout <<"------------------------------------------------------------\n";
  return(0);
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
  


