/// \file sim_control.cpp
/// 
/// \brief Simulation Control Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for sim_control_fixedgrid 
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
/// - JM 2009-12-16 Added ifdef in sim_control_fixedgrid::Time_Int() so that I
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
/// - 2010.10.13 JM: Added option to setup new microphysics_lowZ class
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
/// - 2011.04.14 JM: Added mp_v2_aifa microphysics integration class.
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
/// - 2011.10.13 JM: Added mp_implicit_H class.
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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "sim_control.h"

#include "dataIO/dataio.h"
#include "microphysics/microphysics_base.h"



#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_hydro_adi_Eint.h"
#include "spatial_solvers/solver_eqn_hydro_iso.h"
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


sim_control_fixedgrid::sim_control_fixedgrid()
{
  eqn=0;
  dataio=0;
  textio=0;
  SimPM.checkpoint_freq=INT_MAX;
  max_walltime = 1.0e100;
}



// ##################################################################
// ##################################################################


sim_control_fixedgrid::~sim_control_fixedgrid()
{
#ifdef TESTING
  cout << "(sim_control_fixedgrid::Destructor) Deleting Grid Class..." <<"\n";
#endif
  if(eqn !=0) {delete eqn; eqn=0;}
  if (dataio) {delete dataio; dataio=0;}
  if (textio) {delete textio; textio=0;}
#ifdef TESTING
  cout << "(sim_control_fixedgrid::Destructor) Done." <<"\n";
#endif
}


// ##################################################################
// ##################################################################


double sim_control_fixedgrid::get_max_walltime()
{
  return max_walltime;
}


// ##################################################################
// ##################################################################


void sim_control_fixedgrid::set_max_walltime(
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
void sim_control_fixedgrid::print_command_line_options(
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

  cout <<"\n*********** PHYSICS OPTIONS *************\n";
  cout <<"\t ooa=N         : modify order of accuracy (either 1 or 2).\n";

  cout <<"\t eqntype=N     : modify type of equations,";
  cout <<" Euler=1, Ideal-MHD=2, GLM-MHD=3, FCD-MHD=4, IsoHydro=5.\n";
  cout <<"\t\t\t WARNING -- IT IS DANGEROUS TO CHANGE EQUATIONS!\n";

  cout <<"\t AVtype=N      : modify type of artificial viscosity:";
  cout <<" 0=none, 1=Falle,Komissarov,Joarder(1998), 2=Colella+Woodward(1984), 3=Sanders et al.(1998)[H-correction].\n";
  cout <<"\t\t\t WARNING -- AVtype=2 IS NOT WORKING WELL.  ONLY USE FKJ98.";
  cout <<" The H-correction works well in serial, if needed.\n";
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

  cout <<"\t timestep_limit=N:\n";
  cout <<"\t\t 0 = only dynamical Courant condition.\n";
  cout <<"\t\t 1 = dynamical + cooling time limits\n";
  cout <<"\t\t 2 = dyn +cool +recombination time limits\n";
  cout <<"\t\t 3 = dyn +cool +recomb +ionisation time limits\n";
  cout <<"\t\t 4 = dyn +recomb (NO cool, NO ion)\n";

  cout <<"\n*********** PARALLEL CODE ONLY *************\n";
  cout <<"\t maxwalltime=D : change the max. runtime to D in seconds.\n";
  

  cout <<"\n";
  cout <<"********* DEPRECATED -- STILL HERE FOR LEGACY SCRIPTS... ******\n";
  cout <<"\t artvisc=D : modify artificial viscosity, 0=none, Otherwise FalleAV with eta=D,\n";
  cout <<"\t noise=D   : add noise to initial conditions if desired, at fractional level of D.\n";
  cout <<"     *********************************************\n\n";
  return;
}



// ##################################################################
// ##################################################################


int sim_control_fixedgrid::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      class GridBaseClass **grid
      )
{
#ifdef TESTING
  cout <<"(UniformFV::Init) Initialising grid"<<"\n";
#endif
  int err=0;
  
  SimPM.typeofip=typeOfFile;
  // First get grid parameters from a file.
  switch (typeOfFile) {
  case 1: // Start From ASCII Parameterfile.
    if (!dataio) dataio = new dataio_text();
    if (!dataio) rep.error("dataio_text initialisation",dataio);
    break;
#ifdef FITS
  case 2: // Start from FITS restartfile, which may be initial conditions or a snapshot.
  case 3: // Fits restartfile in table format (slower I/O than image...)
    if (!dataio) dataio = new DataIOFits();
    if (!dataio) rep.error("DataIOFits initialisation",dataio);
    break;
#endif // if FITS
#ifdef SILO
  case 5: // Start from Silo ICfile or restart file.
    if (!dataio) dataio = new dataio_silo ("DOUBLE");
    if (!dataio) rep.error("dataio_silo initialisation",dataio);
    break; 
#endif // if SILO
    //   case : // Start from HDF5 restartfile, which may be initial conditions or a snapshot.
    //    err = getParametersHDF5(infile);
    //    break;
  default:
    cerr <<"(UniformFV::Init) Do not understand typeOfFile="<<typeOfFile<<", so exiting.\n";
    return(1);
  }
  err = dataio->ReadHeader(infile);
  rep.errorTest("(INIT::get_parameters) err!=0 Something went bad",0,err);
  
  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  rep.errorTest("(INIT::override_params) err!=0 Something went bad",0,err);
  
  // Now set up the grid structure.
  cout <<"Init: &grid="<< grid<<", and grid="<< *grid <<"\n";
  err = setup_grid((grid),&mpiPM);
  cout <<"Init: &grid="<< grid<<", and grid="<< *grid <<"\n";
  err += get_cell_size(*grid);
  if (err!=0) {
    cerr<<"(INIT::setup_grid) err!=0 Something went bad"<<"\n";
    return(1);
  }

  //
  // All grid parameters are now set, so I can set up the appropriate
  // equations/solver class.
  //
  err = set_equations();
  rep.errorTest("(INIT::set_equations) err!=0 Fix me!",0,err);

  //
  // Now setup Microphysics, if needed.
  //
  err = setup_microphysics();
  rep.errorTest("(INIT::setup_microphysics) err!=0",0,err);
  
  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, *grid);
  rep.errorTest("(INIT::assign_initial_data) err!=0 Something went bad",0,err);
  //
  // Set Ph=P in every cell.
  //
  cell *cpt = (*grid)->FirstPt();
  do {for(int v=0;v<SimPM.nvar;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=(*grid)->NextPt(cpt))!=0);

  // If we are to add noise to data, do it.
//  if (SimPM.addnoise >0) {
//    cout <<"\tAdding noise to data, at fractional level "<<SimPM.addnoise<<"\n";
//    addNoise2Data(2,SimPM.addnoise); // 1=to pressure, 2=adiabatic+random, 3=adiabatic+wave
//  }


  
  // Assign boundary conditions to boundary points.
  err = boundary_conditions(*grid);
  if (err!=0){cerr<<"(INIT::boundary_conditions) err!=0 Something went bad"<<"\n";return(1);}


  // Now do some checks that everything is ready to start.
#ifdef TESTING
  cout <<"(UniformFV::Init) Ready to start? \n";
#endif
  err = ready_to_start(*grid);
  rep.errorTest("(INIT::ready_to_start) err!=0 Something went bad",0,err);

#ifdef TESTING
  cout <<"(UniformFV::Init) Ready to start. Starting simulation."<<"\n";
#endif
#ifdef SERIAL
  // If outfile-type is different to infile-type, we need to delete dataio and set it up again.
  // This is ifdeffed because parallel version of init() will do it itself,
  // with parallel I/O classes.
  if (SimPM.typeofip != SimPM.typeofop) {
#ifdef TESTING
    cout <<"(UniformFV::INIT) infile-type="<<SimPM.typeofip;
    cout <<" and outfile-type="<<SimPM.typeofop;
    cout <<", so deleting and renewing dataio.\n";
#endif
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    switch (SimPM.typeofop) {
    case 1: // ascii, so do nothing.
      dataio = new dataio_text ();
      break;
#ifdef FITS
    case 2: // fits
    case 3: // fits
      dataio = new DataIOFits();
      break;
    case 4: // fits +ascii
      dataio = new DataIOFits();
      textio = new dataio_text();
      if (!textio) rep.error("INIT:: textio initialisation",SimPM.typeofop);
      break;
#endif // if FITS
#ifdef SILO
    case 5: // silo
      dataio = new dataio_silo ("DOUBLE");
      break;
    case 6: // silo + text
      dataio = new dataio_silo ("DOUBLE");
      textio = new dataio_text ();
      if (!textio) rep.error("INIT:: textio initialisation",SimPM.typeofop);
      break;
#endif // if SILO
    default:
      rep.error("Bad type of Output",SimPM.typeofop);
    }
    if (!dataio) rep.error("INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(eqn);
  if (textio) textio->SetSolver(eqn);

  if (SimPM.timestep==0) {
    cout << "(INIT) Outputting initial data.\n";
    err=output_data(*grid);
    if (err)
      rep.error("Failed to write file!","maybe dir does not exist?");
  }
  cout <<"------------------------------------------------------------\n";
#endif // SERIAL
  
#ifdef TESTING
  cell *c = (*grid)->FirstPt_All();
  do {
    if (pconst.equalD(c->P[RO],0.0))
      CI.print_cell(c);
  } while ( (c=(*grid)->NextPt_All(c)) !=0 );
#endif // TESTING
  
  return(0);
}



// ##################################################################
// ##################################################################




int sim_control_fixedgrid::override_params(int narg, string *args)
{
  cout <<"------------------------------------------------------\n";
  cout <<"--------  Overriding parameters if requested ---------\n";
  // Find command line params and assign them to SimPM.whatever.
  for (int i=4;i<narg;i++) {

    if ( args[i].find("eqntype=") != string::npos) {
      /** \section eqntype Changing Equations
       * Over-riding the equations is tricky, and only some changes are allowed,
       * so each possibility has to be gone through one by one and the number
       * of variables changed, etc.  In short, this is potentially buggy, and 
       * should be done with care, especially when the reference state vector
       * is not just a list of ones.  In fact it's better not to do it, unless
       * you are testing the MHD divergence cleaning algorithm.  One situation 
       * that should work well, is when you are just solving the Euler equations
       * on a dimensionless problem with no tracers or microphysics, and want to
       * solve the same problem with the MHD Riemann Solver.  In this case I have
       * tested it quite a lot and it works fine.
       */
      // assign eqntype 1=hd, 2=mhd, 3=glm-mhd. string is 'eqntype=N'
      int e =  atoi((args[i].substr(8)).c_str());
      if      ((SimPM.eqntype==EQMHD || SimPM.eqntype==EQGLM) && e==1)
  rep.error("Can't override eqntype from mhd to euler.",e);
      else if (SimPM.eqntype ==e) {
  // Don't do anything if not changing eqntype
      }

#ifdef INCLUDE_EINT_ADI_HYDRO
      else if (SimPM.eqntype==EQEUL && e==EQEUL_EINT) {
  //
  // This shouldn't cause any problems.
  //
  cout <<"\tOVERRIDE PARAMS: Resetting eqntype from Euler";
  cout <<" (Etot) to Euler-Eint -- non-conservative.\n";
  SimPM.eqntype=e;
#ifdef        EINT_ETOT_PARALLEL
  SimPM.nvar+=1; SimPM.ftr +=1;
#endif //     EINT_ETOT_PARALLEL
      }
      else if (SimPM.eqntype==EQEUL_EINT && e==EQEUL) {
  //
  // This shouldn't cause any problems,.
  //
  cout <<"\tOVERRIDE PARAMS: Resetting eqntype from Euler-Eint";
  cout <<" (non-conservative) to Euler (total energy).\n";
  SimPM.eqntype=e;
#ifdef        EINT_ETOT_PARALLEL
  SimPM.nvar-=1; SimPM.ftr -=1;
#endif //     EINT_ETOT_PARALLEL
      }
#endif // if INCLUDE_EINT_ADI_HYDRO

      else {
  cout <<"\tOVERRIDE PARAMS: Resetting eqntype from ";
  cout <<SimPM.eqntype<<" to "<<e<<", (1=HD,2=MHD,3=GLM-MHD,9=HD-EINT)\n";
  cout <<"\t WARNING: DON'T OVERRIDE EQNTYPE FOR JET SIMS, AS IT MESSES UP THE JETSTATE VEC FOR THE BCs.\n";
  cout <<"\t GENERIC WARNING: THIS COULD HAVE BIZARRE/UNPREDICTABLE CONSEQUENCES, SO BEST NOT TO DO IT!\n";

  if (SimPM.eqntype==EQEUL ) {
    if      (e==EQMHD || e==EQFCD) {
      SimPM.eqntype=e; SimPM.nvar+=3; SimPM.ftr +=3;
    }
    else if (e==EQGLM) {
      SimPM.eqntype=e; SimPM.nvar+=4; SimPM.ftr +=4;
    }
    else rep.error("What am i changing eqntype to?",e);
  } // if EQEUL

  else if (SimPM.eqntype==EQMHD) {
    if (e==EQGLM) {
      SimPM.eqntype=e; SimPM.nvar+=1; SimPM.ftr +=1;
    }
    else if (e==EQFCD) {
      cout <<"\t\t setting equations to FieldCD method.\n"; SimPM.eqntype=e;
    }
    else rep.error("Can only change eqntype from i-mhd to glm-mhd or fcd-mhd!",e);
  } // if EQMHD

  else if (SimPM.eqntype==EQFCD) {
    if (e==EQGLM) {
      SimPM.eqntype=e; SimPM.nvar+=1; SimPM.ftr +=1;
    }
    else if (e==EQMHD) {
      cout <<"\t\t setting equations to ideal-mhd method.\n"; SimPM.eqntype=e;
    }
    else rep.error("Can only change eqntype from fcd-mhd to glm-mhd or i-mhd!",e);
  } // if EQFCD

  else if (SimPM.eqntype==EQGLM) {
    if (e==EQMHD) {
      SimPM.eqntype=e; SimPM.nvar-=1; SimPM.ftr -=1;
    }
    else if (e==EQFCD) {
      cout <<"\t\t setting equations to FieldCD method.\n";
      SimPM.eqntype=e; SimPM.nvar-=1; SimPM.ftr -=1;
    }
    else rep.error("Can only change eqntype from glm-mhd to i-mhd or FieldCD-mhd!",e);
  } // if EQGLM
      } // if we need to reset eqntype and modify nvar, ftr.
    } // if "eqntype" found in args.
    
    else if (args[i].find("ooa=") != string::npos) {
      // Assign order of accuracy;  string is 'ooa=N', where N=1 or 2.
      int tmp=SimPM.spOOA;
      SimPM.spOOA = atoi((args[i].substr(4)).c_str());
      SimPM.tmOOA = SimPM.spOOA;
      cout <<"\tOVERRIDE PARAMS: Resetting OOA from ooa="<<tmp;
      cout <<" to command-line value = "<<SimPM.spOOA<<"\n";
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
  SimPM.artviscosity=AV_HCORRECTION; // ==3 (NEEDS TESTING!!!)
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
      if (c<0 || c>6 )
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
} // sim_control_fixedgrid::override_params



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::get_cell_size(
      class GridBaseClass *grid 
      )
{
//  cout <<"\tGetting Cell Dimensions from UniformGrid.\n";
  SimPM.dx = grid->DX();
//  SimPM.dA = grid->DA();
//  SimPM.dV = grid->DV();
  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::set_equations()
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

  if (eqn) eqn=0;

  if (SimPM.coord_sys==COORD_CRT) {
    cout <<"set_equations() Using Cartesian coord. system.\n";
    switch (SimPM.eqntype) {
    case EQEUL:
      cout <<"set_equations() Using Euler Equations.\n";
      eqn = new class FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
#ifdef ISOTHERMAL_SOLVERS_ENABLED
    case EQEUL_ISO:
      cout <<"set_equations() Using Isothermal Hydrodynamics Equations.\n";
      eqn = new class FV_solver_Hydro_iso(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL_ISO);
      break;
#endif // ISOTHERMAL_SOLVERS_ENABLED
#ifdef INCLUDE_EINT_ADI_HYDRO
    case EQEUL_EINT:
      cout <<"set_equations() Using Euler Equations, integrating ***INTERNAL ENERGY***.\n";
      eqn = new class FV_solver_Hydro_Euler_Eint(SimPM.nvar,
            SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma,
            SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL_EINT);
      break;      
#endif // if INCLUDE_EINT_ADI_HYDRO
    case EQMHD:
      cout <<"set_equations() Using Ideal MHD Equations.\n";
      eqn = new class FV_solver_mhd_ideal_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQMHD);
      break;
    case EQGLM:
      cout <<"set_equations() Using GLM MHD Equations.\n";
      eqn = new class FV_solver_mhd_mixedGLM_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQGLM);
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
      eqn = new class cyl_FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
#ifdef INCLUDE_EINT_ADI_HYDRO
    case EQEUL_EINT:
      cout <<"set_equations() Using Euler Equations, integrating ***INTERNAL ENERGY***.\n";
      eqn = new class cyl_FV_solver_Hydro_Euler_Eint(SimPM.nvar,
            SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma,
            SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL_EINT);
      break;      
#endif // if INCLUDE_EINT_ADI_HYDRO
    case EQMHD:
      cout <<"set_equations() Using Ideal MHD Equations.\n";
      eqn = new class cyl_FV_solver_mhd_ideal_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQMHD);
      break;
    case EQGLM:
      cout <<"set_equations() Using GLM MHD Equations.\n";
      eqn = new class cyl_FV_solver_mhd_mixedGLM_adi(SimPM.nvar, SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma, SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQGLM);
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
      eqn = new class sph_FV_solver_Hydro_Euler(SimPM.nvar, SimPM.ndim, SimPM.CFL,
            SimPM.dx, SimPM.gamma, SimPM.RefVec,
            SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL);
      break;
#ifdef INCLUDE_EINT_ADI_HYDRO
    case EQEUL_EINT:
      cout <<"set_equations() Using Euler Equations, integrating ***INTERNAL ENERGY***.\n";
      eqn = new class sph_FV_solver_Hydro_Euler_Eint(SimPM.nvar,
            SimPM.ndim, SimPM.CFL, SimPM.dx, SimPM.gamma,
            SimPM.RefVec, SimPM.etav, SimPM.ntracer);
      if (!eqn) rep.error("Couldn't set up solver/equations class.",EQEUL_EINT);
      break;      
#endif // if INCLUDE_EINT_ADI_HYDRO
     default:
      rep.error("not implemented yet for spherical (only hydro in 1D!)",SimPM.eqntype);
    }
  } // axisymmetric


  //
  // Check that we set up an equations class!
  //
  if (!eqn)
    rep.error("sim_control_fixedgrid::set_equations() Failed to initialise an equations class.",eqn);

  cout <<"------------------------------------------------------\n\n";
  return(0);
} // set_equations.



// ##################################################################
// ##################################################################




int sim_control_fixedgrid::output_data(
      class GridBaseClass *grid
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
    err = dataio->OutputData(SimPM.outFileBase, grid, checkpoint_id);
    if (textio) err += textio->OutputData(SimPM.outFileBase, grid, checkpoint_id);
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
  cout <<"\tOutputting data, at simtime: "<<SimPM.simtime << " to file "<<SimPM.outFileBase<<"\n";
  err = dataio->OutputData(SimPM.outFileBase, grid, SimPM.timestep);
  if (textio) err += textio->OutputData(SimPM.outFileBase, grid, SimPM.timestep);
  if (err) {cerr<<"\t Error writing data.\n"; return(1);}
  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::ready_to_start(
      class GridBaseClass *grid
      )
{
  //
  // If I'm using the GLM method, make sure Psi variable is initialised.
  //
  cell *c=0;
  if (SimPM.eqntype==EQGLM && SimPM.timestep==0) {
#ifdef TESTING
    cout <<"Initial state, zero-ing glm variable.\n";
#endif
    c = grid->FirstPt(); do {
      c->P[SI] = c->Ph[SI] = 0.;//grid->divB(c);
    } while ( (c=grid->NextPt(c)) !=0);
  }
  
  //
  // Set the value of gamma in the equations.
  //
  eqn->SetEOS(SimPM.gamma);
  
  //
  // Setup Raytracing, if they are needed.
  //
  int err=0;
  err += setup_raytracing(grid);
  err += setup_evolving_RT_sources();
  if (err) rep.error("Failed to setup raytracer and/or microphysics",err);
  //
  // If we are doing raytracing, we probably want to limit the timestep by
  // the ionisation timescale, so we call a very short time update here (to
  // advance time by 1 second) in order to set the column densities in each 
  // cell so that calc_timestep() can get the correct reaction rates.
  //

  //
  // If testing the code, this calculates the momentum and energy on the domain.
  //
  initial_conserved_quantities(grid);
  
  //
  // If using opfreq_time, set the next output time correctly.
  // 'floor' is a cmath function returning the largest integer NOT LARGER
  // THAN the double precision argument (should be identical to (int (x)) 
  // at least for x>0).
  //
  if (SimPM.op_criterion==1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!",SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
    double tmp = ((SimPM.simtime/SimPM.opfreq_time)-floor(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
    SimPM.next_optime-= tmp;
  }

  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::setup_evolving_RT_sources()
{
  //
  // Loop through list of sources, and see if any of them have an evolution
  // file (if none, then the string is set to NOFILE in the data I/O stage).
  //
  int err=0;
  int Nevo=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].EvoFile == "NOFILE") {
#ifdef TESTING
      cout <<"setup_evolving_RT_sources() Source "<<isrc<<" has no evolution file.\n";
#endif
    }
    else {
      if (SimPM.RS.sources[isrc].effect != RT_EFFECT_PION_MULTI) {
        rep.error("setup_evolving_RT_sources() Source is not multifreq but has EvoFile",isrc);
      }
      Nevo++;
      struct star istar;
#ifdef TESTING
      cout <<"setup_evolving_RT_sources() Source "<<isrc<<" has EvoFile "<<istar.file_name<<"\n";
#endif
      istar.file_name = SimPM.RS.sources[isrc].EvoFile;
      istar.src_id    = isrc;
      SimPM.STAR.push_back(istar);
    }
  }
  //
  // Now go through each one we found and read the evolution file into arrays
  // and set the rest of the data in the 'star' struct.
  //
  for (int isrc=0; isrc<Nevo; isrc++) {
    struct star *istar = &(SimPM.STAR[isrc]);
    //
    // Open file
    //
    ifstream infile(istar->file_name.c_str());
    if (!infile.is_open())
      rep.error("setup_evolving_RT_sources() Opening EvoFile",istar->file_name);
      
    //
    // Read Nlines, and put data into time, Log_L, Log_T, Log_R, Log_V.
    //
    string line;
    getline(infile, line);
    while (line.empty() == true || line.substr(0,1) == "#") {
      getline(infile, line);
    }
    //
    // Now first non-comment line should tell me Nlines, which we use to resize
    // the arrays.
    //
    istringstream iss2(line);
    string junk;
    iss2 >> junk >> istar->Nlines;
    //cout <<"\t\tgetting Nlines:: "<<junk<<": "<<istar->Nlines<<"\n";
    if (istar->Nlines>1000000 || istar->Nlines<2) {
      rep.error("setup_evolving_RT_sources() Bad Nlines in stellar radiation evolution",istar->Nlines);
    }
    istar->time.resize(istar->Nlines);
    istar->Log_L.resize(istar->Nlines);
    istar->Log_T.resize(istar->Nlines);
    istar->Log_R.resize(istar->Nlines);
    istar->Log_V.resize(istar->Nlines);
    size_t i=0;
    while (!infile.eof()) {
      getline(infile, line);
      istringstream iss(line);
      if (i>istar->Nlines)
        rep.error("setup_evolving_RT_sources() Too many lines in EvoFile",iss.str());
      if (!line.empty()) {
        iss >> istar->time[i] >> istar->Log_L[i] >> istar->Log_T[i] >> istar->Log_R[i] >> istar->Log_V[i];
        //
        // Rescale time to seconds
        //
        istar->time[i] *= pconst.year();
        
        //
        // For ionisation rate, we need to rescale the Blackbody luminosity so
        // that it is much smaller for T<30000K, b/c the actual ionising photon
        // luminosity of these stars is much less than indicated by BB curve.
        // I took data from Table 1 of Diaz-Miller, Franco, & Shore,
        // (1998,ApJ,501,192), compared them to the ionising photon luminosity
        // of a BB with the same radius and Teff, and got the following scaling
        // factor:
        //
        if (istar->Log_T[i]<4.55555) {
          //cout <<"L(BB) ="<<exp(pconst.ln10()*istar->Log_L[i])<<", T=";
          //cout <<exp(pconst.ln10()*istar->Log_T[i])<<", scale-factor=";
          //cout << exp(pconst.ln10()*(9.0*istar->Log_T[i] -41.0));
          istar->Log_L[i] += 9.0*istar->Log_T[i] -41.0;
          //cout <<", new L = "<<exp(pconst.ln10()*istar->Log_L[i])<<"\n";
        }
        //
        // calculate radius from L=4.pi.R^2.sigma.T^4 (I want to get rid of
        // this altogether soon).
        //
        istar->Log_R[i] =
          0.5*(istar->Log_L[i]-4.0*istar->Log_T[i]+log10(pconst.Lsun()/
              (4.0*M_PI*pconst.StefanBoltzmannConst()))) -log10(pconst.Rsun());

        //cout <<istar->time[i]<<"\t"<< istar->Log_L[i] <<"\t"<< istar->Log_T[i] <<"\t"<< istar->Log_R[i] <<"\t"<< istar->Log_V[i] <<"\n";
        i++;
      }
    }
    if (i!=istar->Nlines) {
      rep.error("setup_evolving_RT_sources() Too few lines in file!",istar->Nlines-i);
    }
    infile.close();

    //
    // Finally set the last_line counter to be the array index nearest to
    // (but less than) the current time.
    //
    i=0;
    while (istar->time[i] < SimPM.simtime) i++;
    istar->last_line = i;

    // initialise to zero.
    istar->Lnow = istar->Tnow = istar->Rnow = istar->Vnow = 0.0;
  }

  //
  // All done setting up the source.  Now we need to update the SimPM.RS.
  // source properties and send the changes to the raytracing class.  Need
  // time in secs, L,T,V in cgs and R in Rsun.
  //
  err = update_evolving_RT_sources();
  return err;
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::update_evolving_RT_sources()
{
  int err=0;
  bool updated=false;
  //
  // Loop over all sources with Evolution files.
  //
  for (unsigned int isrc=0; isrc<SimPM.STAR.size(); isrc++) {
    struct star *istar = &(SimPM.STAR[isrc]);
    istar->t_now = SimPM.simtime;
    size_t i = istar->last_line;

    //
    // Check if we have reached the last line of the file!
    //
    if (i==(istar->Nlines-1)) {
      cout <<"\n\n*#*#*#*#* WARNING #*#*#*#*#*#* update_evolving_RT_sources()";
      cout <<" Last line, assuming star is constant luminosity from now on!\n\n";
      return 0;
    }
    //
    // Check if we have moved forward one line in the code, in which
    // case we need to increment i.
    //
    while (istar->t_now > istar->time[i+1]) {
      //cout <<"update_evolving_RT_sources() Source has moved to next line. i="<<i<<" time="<<istar->time[i]<<"\n";
      i++;
      istar->last_line = i;
    }

    //
    // Now we know the star properties are bracketed by line i and line i+1,
    // so we can do a simple linear interpolation between them.
    //
    // First interpolate in log space.
    //
    double Lnow, Tnow, Rnow, Vnow;
    Lnow = istar->Log_L[i] +(istar->t_now-istar->time[i])*
        (istar->Log_L[i+1]-istar->Log_L[i])/(istar->time[i+1]-istar->time[i]);
    Tnow = istar->Log_T[i] +(istar->t_now-istar->time[i])*
        (istar->Log_T[i+1]-istar->Log_T[i])/(istar->time[i+1]-istar->time[i]);
    Rnow = istar->Log_R[i] +(istar->t_now-istar->time[i])*
        (istar->Log_R[i+1]-istar->Log_R[i])/(istar->time[i+1]-istar->time[i]);
    Vnow = istar->Log_V[i] +(istar->t_now-istar->time[i])*
        (istar->Log_V[i+1]-istar->Log_V[i])/(istar->time[i+1]-istar->time[i]);
    //
    // Now convert units (Radius is ok, but all others need conversion).
    //
    Lnow = exp(pconst.ln10()*(Lnow))*pconst.Lsun();
    Tnow = exp(pconst.ln10()*(Tnow));
    Rnow = exp(pconst.ln10()*(Rnow));
    Vnow = exp(pconst.ln10()*(Vnow));

    //
    // If L or T change by more than 1% then update them; otherwise leave as they are.
    //
    if ( fabs(Lnow-istar->Lnow)/istar->Lnow >0.01 || fabs(Tnow-istar->Tnow)/istar->Tnow >0.01 ) {
      cout <<"update_evolving_RT_sources() NOW: t="<<istar->t_now<<"\t"<< Lnow <<"\t"<< Tnow <<"\t"<< Rnow <<"\t"<< Vnow <<"\n";
      istar->Lnow = Lnow;
      istar->Tnow = Tnow;
      istar->Rnow = Rnow;
      istar->Vnow = Vnow;

      //
      // Copy new data into SimPM.RS, and send updates to RayTracing and
      // microphysics classes.
      //
      struct rad_src_info *rs = &(SimPM.RS.sources[istar->src_id]);
      
      rs->strength = istar->Lnow;
      rs->Rstar    = istar->Rnow;
      rs->Tstar    = istar->Tnow;
      
      RT->update_RT_source_properties(rs);

      err += MP->set_multifreq_source_properties(rs);
      if (err) rep.error("update_evolving_RT_sources() failed to update MP for source id",rs->id);

      updated=true;
    }

  }
  

  //
  // Finally get the data back from RT into the structs for MP updates.
  //
  if (updated) {
    RT->populate_UVheating_src_list(FVI_heating_srcs);
    RT->populate_ionising_src_list( FVI_ionising_srcs);
  }
  
  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::initial_conserved_quantities(
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
     eqn->PtoU(cpt->P,u,SimPM.gamma);
     dp.initERG += u[ERG]*eqn->CellVolume(cpt);
     dp.initMMX += u[MMX]*eqn->CellVolume(cpt);
     dp.initMMY += u[MMY]*eqn->CellVolume(cpt);
     dp.initMMZ += u[MMZ]*eqn->CellVolume(cpt);
  } while ( (cpt = grid->NextPt(cpt)) !=0);
  //cout <<"!!!!! cellvol="<<eqn->CellVolume(cpt)<< "\n";
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
int sim_control_fixedgrid::Time_Int(
      class GridBaseClass *grid
      )
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control_fixedgrid::Time_Int) STARTING TIME INTEGRATION."<<"\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  clk.start_timer("Time_Int"); double tsf=0;
  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid);
#endif
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources();
    if (err) {
      cerr <<"(TIME_INT::update_evolving_RT_sources()) something went wrong!\n";
      return err;
    }

    //clk.start_timer("advance_time");
    err+= advance_time(grid);
    //cout <<"advance_time took "<<clk.stop_timer("advance_time")<<" secs.\n";
    if (err!=0){cerr<<"(TIME_INT::advance_time) err!=0 Something went bad\n";return(1);}

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
    tsf=clk.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err+= output_data(grid);
    if (err!=0){cerr<<"(TIME_INT::output_data) err!=0 Something went bad\n";return(1);}

    err+= check_eosim();
    if (err!=0){cerr<<"(TIME_INT::) err!=0 Something went bad\n";return(1);}
  }

  cout <<"(sim_control_fixedgrid::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

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
void sim_control_fixedgrid::calculate_magnetic_pressure(
      class GridBaseClass *grid
      )
{
  //
  // Calculate the total magnetic pressure on the domain, normalised to the
  // initial value.
  //
  cell *c=grid->FirstPt();
  double magp=0.0;
  static double init_magp=-1.0;
  do {
    magp += eqn->Ptot(c->P,0.0) - c->P[PG];
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
/// If running a spherical blast wave, calculate the shock position
/// and output to screen.
///
void sim_control_fixedgrid::calculate_blastwave_radius(
      class GridBaseClass *grid
      )
{
  //
  // Calculate the blast wave outer shock position.
  //
  cell *c=grid->LastPt();
  double shockpos=0.0;
  static double old_pos=0.0;
  //  static double last_dt=0.0;
  do {
    c = grid->NextPt(c,RNsph);
    //cout <<c->id<<", vx="<<c->P[VX]<<"\n";
  } while ( c!=0 && fabs(c->P[VX])<1.0e4);
  if (!c)
    shockpos = 0.0;
  else
    shockpos = CI.get_dpos(c,Rsph);
  if (old_pos<0.4*grid->DX()) old_pos = shockpos;
  cout <<SimPM.simtime<<"\t"<<shockpos;
  //cout <<"\t"<<(shockpos-old_pos)/(SimPM.dt+TINYVALUE);
  cout <<"\n";
  old_pos=shockpos;
  return;
}
#endif // BLAST_WAVE_CHECK



// ##################################################################
// ##################################################################




int sim_control_fixedgrid::calc_timestep(
      class GridBaseClass *grid
      )
{
  //
  // This is now basically a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt(grid);
  t_mp  = calc_microphysics_dt(grid);
#ifdef TESTING
  if (t_mp<t_dyn)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";
#endif
  SimPM.dt = min(t_dyn,t_mp);

#ifdef THERMAL_CONDUCTION
  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to actually calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in eqn->preprocess_data()
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
  // if using MHD with GLM divB cleaning, the following sets the hyperbolic wavespeed.
  // If not, it does nothing.  By setting it here and using t_dyn, we ensure that the
  // hyperbolic wavespeed is equal to the maximum signal speed on the grid, and not
  // an artificially larger speed associated with a shortened timestep.
  //
  eqn->GotTimestep(t_dyn);

  //
  // Check that the timestep doesn't increase too much between steps, and that it 
  // won't bring us past the next output time or the end of the simulation.
  // This function operates on SimPM.dt, resetting it to a smaller value if needed.
  //
  timestep_checking_and_limiting();

  // sets the timestep info in the solver class.
  eqn->Setdt(SimPM.dt);

  return 0;
}



// ##################################################################
// ###########    END OF SIMULATION CONTROL FUNCTIONS     ###########
// ##################################################################




// ##################################################################
// ##################################################################



int sim_control_fixedgrid::check_eosim()
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



int sim_control_fixedgrid::check_energy_cons(
      class GridBaseClass *grid
      )
{
#ifdef TESTING
  // Energy, and Linear Momentum in x-direction.
  pion_flt u[SimPM.nvar];
  double ergNow=0., mmxNow = 0., mmyNow = 0., mmzNow = 0.;
  class cell *cpt = grid->FirstPt();
  do {
     eqn->PtoU(cpt->P,u,SimPM.gamma);
     ergNow += u[ERG]*eqn->CellVolume(cpt);
     mmxNow += u[MMX]*eqn->CellVolume(cpt);
     mmyNow += u[MMY]*eqn->CellVolume(cpt);
     mmzNow += u[MMZ]*eqn->CellVolume(cpt);
  } while ( (cpt =grid->NextPt(cpt)) !=0);
  //cout <<"!!!!! cellvol="<<eqn->CellVolume(cpt)<< "\n";
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
int sim_control_fixedgrid::Finalise(
      class GridBaseClass *grid
      )
{
  int err=0;
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control_fixedgrid::Finalise) FINALISING SIMULATION."<<"\n";
  err += check_energy_cons(grid);
  err+= output_data(grid);
  if (err!=0){
    cerr<<"(FINALISE::output_data) final state data output. err!=0 Something went bad"<<"\n";
    return(1);
  }
  cout <<"\tSimTime = "<<SimPM.simtime<<"   #timesteps = "<<SimPM.timestep<<"\n";
#ifdef TESTING
  cout <<"(sim_control_fixedgrid::Finalise) DONE.\n";
#endif
  cout <<"------------------------------------------------------------\n";
  return(0);
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
  


