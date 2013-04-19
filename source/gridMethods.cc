/// \file gridMethods.cc
/// 
/// \brief Grid Methods Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for IntUniformFV 
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
/// - JM 2009-12-16 Added ifdef in IntUniformFV::Time_Int() so that I
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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid.h"
#include "dataIO/dataio.h"
#include "microphysics/microphysics_base.h"

#ifndef EXCLUDE_MPV1
#include "microphysics/microphysics.h"
#endif 

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/microphysics_lowZ.h"
#endif

#include "microphysics/mp_only_cooling.h"

#ifndef EXCLUDE_MPV2
#ifdef MP_V2_AIFA
#include "microphysics/mp_v2_aifa.h"
#endif
#endif 

#ifndef EXCLUDE_MPV3
#include "microphysics/mp_explicit_H.h"
#endif

#ifndef EXCLUDE_MPV4
#include "microphysics/mp_implicit_H.h"
#endif 

#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"



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


IntUniformFV::IntUniformFV()
{
  eqn=0;
  dataio=0;
  textio=0;
  SimPM.checkpoint_freq=INT_MAX;
  FVI_nheat = FVI_nion = 0;
  FVI_heating_srcs.clear();
  FVI_ionising_srcs.clear();
  FVI_need_column_densities_4dt = false;
}



// ##################################################################
// ##################################################################


IntUniformFV::~IntUniformFV()
{
#ifdef TESTING
  cout << "(IntUniformFV::Destructor) Deleting Grid Class..." <<"\n";
#endif
  if(grid!=0) {
#ifdef TESTING
    cout << "\t Deleting Grid Data..." <<"\n";
#endif
    delete grid; grid=0;
  }
  if(eqn !=0) {
#ifdef TESTING
    cout <<"\t Deleting Solver/Equations class...\n";
#endif
    delete eqn; eqn=0;
  }
  if (dataio) {delete dataio; dataio=0;}
  if (textio) {delete textio; textio=0;}
  if (MP)     {delete MP; MP=0;}
  if (RT)     {delete RT; RT=0;}
#ifdef TESTING
  cout << "(IntUniformFV::Destructor) Done." <<"\n";
#endif
}

/*****************************************************************/
/********************* SIM. INITIALISATION ***********************/
/*****************************************************************/



// ##################################################################
// ##################################################################


int IntUniformFV::Init(string infile, int typeOfFile, int narg, string *args)
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
    if (!dataio) dataio = new dataio_silo ();
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
  
  // Now determine what to do at boundaries, and setup the grid structure with boundaries.
  //  cout <<"(UniformFV::Init) Setting up grid... \n";
  err = setup_grid();
  err += get_cell_size();
  if (err!=0){cerr<<"(INIT::setup_grid) err!=0 Something went bad"<<"\n";return(1);}

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
  err = dataio->ReadData(infile,grid);
  rep.errorTest("(INIT::assign_initial_data) err!=0 Something went bad",0,err);
  //
  // Set Ph=P in every cell.
  //
  cell *cpt = grid->FirstPt();
  do {for(int v=0;v<SimPM.nvar;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=grid->NextPt(cpt))!=0);

  // If we are to add noise to data, do it.
//  if (SimPM.addnoise >0) {
//    cout <<"\tAdding noise to data, at fractional level "<<SimPM.addnoise<<"\n";
//    addNoise2Data(2,SimPM.addnoise); // 1=to pressure, 2=adiabatic+random, 3=adiabatic+wave
//  }


  
  // Assign boundary conditions to boundary points.
  err = boundary_conditions();
  if (err!=0){cerr<<"(INIT::boundary_conditions) err!=0 Something went bad"<<"\n";return(1);}


  // Now do some checks that everything is ready to start.
#ifdef TESTING
  cout <<"(UniformFV::Init) Ready to start? \n";
#endif
  err = ready_to_start();
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
      dataio = new dataio_silo ();
      break;
    case 6: // silo + text
      dataio = new dataio_silo ();
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
    err=output_data();
    if (err)
      rep.error("Failed to write file!","maybe dir does not exist?");
  }
  cout <<"------------------------------------------------------------\n";
#endif // SERIAL
  
  return(0);
}



// ##################################################################
// ##################################################################




int IntUniformFV::override_params(int narg, string *args)
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
} // IntUniformFV::override_params



// ##################################################################
// ##################################################################



int IntUniformFV::get_cell_size()
{
//  cout <<"\tGetting Cell Dimensions from UniformGrid.\n";
  SimPM.dx = grid->DX();
//  SimPM.dA = grid->DA();
//  SimPM.dV = grid->DV();
  return(0);
}



// ##################################################################
// ##################################################################


int IntUniformFV::set_equations()
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
    rep.error("IntUniformFV::set_equations() Failed to initialise an equations class.",eqn);

  cout <<"------------------------------------------------------\n\n";
  return(0);
} // set_equations.



// ##################################################################
// ##################################################################



void IntUniformFV::setup_cell_extra_data()
{
  //
  // Cells can need extra data for ray-tracing optical depths, eta-values for the
  // H-correction or div(v) for some time-updates and/or viscosity corrections.
  //

  int hc_flag = 0, dv_flag=0;
  if (SimPM.artviscosity==AV_LAPIDUS ||
      SimPM.eqntype==EQEUL_EINT) {
    // Need one var. for Div(v)
    dv_flag = 1;
  }
  if (SimPM.artviscosity==AV_HCORRECTION ||
      SimPM.artviscosity==AV_HCORR_FKJ98 ||
      SimPM.eqntype==EQEUL_EINT) {
    //
    // need one var for each dimension here.  For H-correction they
    // are for the eta values.  For EQEUL_EINT we need von Neunamm-Richtmeyer
    // viscosity which needs storage for the diagonal Q-values along each axis.
    //
    hc_flag = SimPM.ndim;
  }

  CI.setup_extra_data(SimPM.RS, hc_flag, dv_flag);
  return;
}



// ##################################################################
// ##################################################################



int IntUniformFV::setup_grid()
{
  cout <<"------------------------------------------------------\n";
  cout <<"--------  Setting up computational grid --------------\n";
  if (SimPM.gridType!=1) {
    rep.warning("gridType not set correctly: Only know Uniform finite\
                 volume grid, so resetting to 1!",1,SimPM.gridType);
    SimPM.gridType=1;
  }
  if (SimPM.ndim <1 || SimPM.ndim>3)  rep.error("Only know 1D,2D,3D methods!",SimPM.ndim);

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  setup_cell_extra_data();

  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
#endif
  if (grid) rep.error("Grid already set up!",grid);

#ifdef GEOMETRIC_GRID
  if      (SimPM.coord_sys==COORD_CRT)
    grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_CYL)
    grid = new uniform_grid_cyl (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_SPH)
    grid = new uniform_grid_sph (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else 
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
#else  // GEOMETRIC_GRID
  grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
#endif // GEOMETRIC_GRID

  if (grid==0) rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
#ifdef TESTING
  cout <<"(UniformFV::setup_grid) Done. g="<<grid<<"\n";
  cout <<"DX = "<<grid->DX()<<"\n";
#endif
  cout <<"------------------------------------------------------\n\n";

  return(0);
} // setup_grid()




// ##################################################################
// ##################################################################


int IntUniformFV::boundary_conditions()
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"(UniformFV::boundary_conditions)";
  
  if(SimPM.typeofbc=="FIXED") {
    cout << "\t Using fixed boundary conditions; only useful for shock tube."<<"\n";
  }
  else if (SimPM.typeofbc=="PERIODIC") {
    cout <<"\t Using periodic BCs on all sides of grid.\n";
  }
  else if (SimPM.typeofbc=="ABSORBING") {
    cout <<"\t Using absorbing BCs on all sides of grid.\n";
  }
  else if (SimPM.typeofbc=="REFLECTING") {
    cout <<"\t Using reflecting BCs on all sides of grid.\n";
  }
  else {
    cout <<"\t Using the following BCs: "<<SimPM.typeofbc<<"\n";
  }
#endif

  // Nbc is the depth of the boundary layer.
  //  cout <<"Setting depth of boundary cells to be equal to spatial OOA.\n";
  if      (SimPM.spOOA==OA2) SimPM.Nbc = 2;
  else if (SimPM.spOOA==OA1) SimPM.Nbc = 1;
  else rep.error("Spatial order of accuracy unhandled by boundary conditions!",SimPM.spOOA);
  
  if (SimPM.solverType==FLUX_LF) {SimPM.spOOA = SimPM.tmOOA = OA1; SimPM.Nbc=1;} // force this if LF Method.
  
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<SimPM.Nbc<<"\n";
#endif
  int err = grid->SetupBCs(SimPM.Nbc,SimPM.typeofbc);
  if (err) rep.error("IntUniformFV::boundary_conditions() Couldn't \
                      set up boundary conditions class.",err);
#ifdef TESTING
  cout <<"(IntUniformFV::boundary_conditions) Done.\n";
#endif
  return 0;
}




// ##################################################################
// ##################################################################




int IntUniformFV::output_data()
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
    if (!GS.equalD(SimPM.simtime,SimPM.next_optime) && (SimPM.maxtime==false))
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



int IntUniformFV::ready_to_start()
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
  err += setup_raytracing();
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
  initial_conserved_quantities();
  
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



int IntUniformFV::setup_evolving_RT_sources()
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
        istar->time[i] *= GS.s_per_yr();
        
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
          //cout <<"L(BB) ="<<exp(GS.ln10()*istar->Log_L[i])<<", T=";
          //cout <<exp(GS.ln10()*istar->Log_T[i])<<", scale-factor=";
          //cout << exp(GS.ln10()*(9.0*istar->Log_T[i] -41.0));
          istar->Log_L[i] += 9.0*istar->Log_T[i] -41.0;
          //cout <<", new L = "<<exp(GS.ln10()*istar->Log_L[i])<<"\n";
        }
        //
        // calculate radius from L=4.pi.R^2.sigma.T^4 (I want to get rid of
        // this altogether soon).
        //
        istar->Log_R[i] =
          0.5*(istar->Log_L[i]-4.0*istar->Log_T[i]+log10(GS.Lsun()/
              (4.0*M_PI*GS.StefanBoltzmannConst()))) -log10(GS.Rsun());

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



int IntUniformFV::update_evolving_RT_sources()
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
    Lnow = exp(GS.ln10()*(Lnow))*GS.Lsun();
    Tnow = exp(GS.ln10()*(Tnow));
    Rnow = exp(GS.ln10()*(Rnow));
    Vnow = exp(GS.ln10()*(Vnow));

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



int IntUniformFV::setup_microphysics()
{
  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  //
  // Setup Microphysics class, if needed.
  // First see if we want the only_cooling class (much simpler), and if
  // not then check for the one of the bigger microphysics classes.
  //
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout <<"\t******* Requested cooling but no chemistry... setting";
    cout <<" up mp_only_cooling() class, with timestep-limiting.\n";
    SimPM.EP.MP_timestep_limit = 1;
    MP = new mp_only_cooling(SimPM.nvar, &(SimPM.EP));
    if (!MP) rep.error("mp_only_cooling() init",MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    cout <<"TRTYPE: "<<SimPM.trtype<<"\n";
    string mptype;
    if (SimPM.trtype.size() >=6)
      mptype = SimPM.trtype.substr(0,6); // Get first 6 chars for type of MP.
    else mptype = "None";
    bool have_set_MP=false;


#ifndef EXCLUDE_MPV1
    if      (mptype=="ChAH__" || mptype=="onlyH_") {
      cout <<"\t******* setting up MP_Hydrogen microphysics module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MP_Hydrogen(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      cout <<"\t**---** WARNING, THIS MODULE HAS BEEN SUPERSEDED BY mp_implicit_H. **--**\n";
      have_set_MP=true;
    }
#endif // exclude MPv1


#ifndef EXCLUDE_HD_MODULE
    if (mptype=="lowZ__") {
      cout <<"\t******* setting up microphysics_lowz module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new microphysics_lowz(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }
#endif // exclude Harpreet's module

#ifndef EXCLUDE_MPV2
    if (mptype=="MPv2__") {
#ifdef MP_V2_AIFA
      cout <<"\t******* setting up mp_v2_aifa module *********\n";
      cout <<"\t******* N.B. Timestep limiting is enforced. **\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_v2_aifa(SimPM.nvar, SimPM.ntracer, SimPM.trtype);
      SimPM.EP.MP_timestep_limit = 1;
#else
      rep.error("Enable mp_v2_aifa as an ifdef if you really want to use it",2);
#endif
      have_set_MP=true;
    }
#endif // exclude MPv2


#ifndef EXCLUDE_MPV3
    if (mptype=="MPv3__") {
      cout <<"\t******* setting up mp_explicit_H module *********\n";
#if MPV3_DTLIMIT>=0 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. Timestep limiting is enforced by #def";
      cout <<" MPV3_DTLIMIT="<<MPV3_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif

      MP = new mp_explicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP)
      );
      //if (SimPM.EP.MP_timestep_limit != 1)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv3


#ifndef EXCLUDE_MPV4
    if (mptype=="MPv4__") {
      cout <<"\t******* setting up mp_implicit_H module *********\n";
#if MPV4_DTLIMIT>=5 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. dt05-12 Timestep limiting is enforced by #def";
      cout <<" DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =5;
#elif MPV4_DTLIMIT>=0 && MPV4_DTLIMIT<=4
      cout <<"\t******* N.B. dt00-04 Timestep limiting is enforced by #def";
      cout <<" MPV4_DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =4;
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_implicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      //SimPM.EP.MP_timestep_limit = 4;  // limit by recombination time only
      //if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv4


    if (mptype=="MPv5__") {
      cout <<"\t******* setting up mpv5_molecular module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv5_molecular(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv6__") {
      cout <<"\t******* setting up mpv6_PureH module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv6_PureH(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv7__") {
      cout <<"\t******* setting up mpv7_TwoTempIso module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv7_TwoTempIso(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }




#ifndef EXCLUDE_MPV1
    //
    // Finally, if MP has not been set up yet, try to set up the v0
    // microphysics integrator, which is slow, but can model a number
    // of elements and ions.
    //
    if (!have_set_MP) {
      cout <<"\t******* setting up MicroPhysics (v0) module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MicroPhysics(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
        rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv1/0

    if (!MP) rep.error("microphysics init",MP);
    if (!have_set_MP) rep.error("HUH? have_set_MP",have_set_MP);
  }
  else {
    cout <<"\t******** not doing microphysics.\n";
    MP=0;
  }

  //
  // If we have a multifrequency ionising source, we can set its properties here.
  // We can only have one of these, so safe to just loop through...
  //
  int err=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE &&
        SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI &&
        MP!=0
        ) {
      err = MP->set_multifreq_source_properties(&SimPM.RS.sources[isrc]);
    }
  }
  if (err) rep.error("Setting multifreq source properties",err);
  

  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  return 0;
}



// ##################################################################
// ##################################################################



int IntUniformFV::setup_raytracing()
{
  //
  // If not doing raytracing, return immediately.
  //
  if (!SimPM.EP.raytracing) {
    return 0;
  }

  //
  // Now we are doing raytracing, so set up a raytracer and add sources to it.
  //
  if (!MP) rep.error("can't do raytracing without microphysics",MP);
  cout <<"\n----------------- RAYTRACER SETUP STARTING -----------------------\n";
  RT=0;
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays=true;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++)
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays=false;
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    RT = new raytracer_USC_infinity(grid,MP);
    if (!RT) rep.error("init pllel-rays raytracer error",RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    RT = new raytracer_USC(grid,MP);
    if (!RT) rep.error("init raytracer error 2",RT);
  }

  //
  // Now add the sources to the tracer.  Note that both the implicit and explicit
  // integrators can still only handle a single ionising source, so we do a check
  // for this and bug out if there is more than one.
  //
  int ion_count=0, uv_count=0, dif_count=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) {
      //
      // single sources have a flux (if at infinity) or a luminosity (if point
      // sources.
      //
      cout <<"Adding IONISING or UV single-source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MONO ||
          SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI)
        ion_count++;
      else 
        uv_count++;
    } // if ionising source
    else {
      // note that diffuse radiation must be at infinity, and the strength is assumed to
      // be an intensity not a flux, so it is multiplied by a solid angle appropriate
      // to its location in order to get a flux.
      cout <<"Adding DIFFUSE radiation source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      uv_count++;
      dif_count++;
    } // if diffuse source
  } // loop over sources
  if (ion_count>1) {
    rep.error("Can only have one ionising source for currently implemented method",ion_count);
  }
  cout <<"Added "<<ion_count<<" ionising and "<<uv_count<<" non-ionising";
  cout <<" radiation sources, of which "<<dif_count<<" are diffuse radiation.\n";
  RT->Print_SourceList();

  //
  // Now that we have added all of the sources, we query the raytracer to get
  // all of the source properties into structs for the microphysics calls.
  // NOTE THAT IF THE NUMBER OF SOURCES OR THEIR PROPERTIES CHANGE OVER TIME,
  // I WILL HAVE TO WRITE NEW CODE TO UPDATE THIS!
  //
  FVI_nheat = RT->N_heating_sources();
  FVI_nion  = RT->N_ionising_sources();
  FVI_heating_srcs.resize(FVI_nheat);
  FVI_ionising_srcs.resize(FVI_nion);
  RT->populate_UVheating_src_list(FVI_heating_srcs);
  RT->populate_ionising_src_list( FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (RT->type_of_RT_integration()==RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT
            && SimPM.EP.MP_timestep_limit==5) {
    // For implicit updates to limit by xdot and/or edot
    // Here the raytracing has not already been done, so we call it here.
    FVI_need_column_densities_4dt = true;
  }
  else {
    FVI_need_column_densities_4dt = false;
  }

  cout <<"----------------- RAYTRACER SETUP COMPLETE -----------------------\n";
  return 0;
}



// ##################################################################
// ##################################################################



int IntUniformFV::initial_conserved_quantities()
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TESTING 
  // Only track the totals if I am testing the code.
  double u[SimPM.nvar];
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
int IntUniformFV::Time_Int()
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(IntUniformFV::Time_Int) STARTING TIME INTEGRATION."<<"\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  GS.start_timer("Time_Int"); double tsf=0;
  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure();
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius();
#endif
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources();
    if (err) {
      cerr <<"(TIME_INT::update_evolving_RT_sources()) something went wrong!\n";
      return err;
    }

    //GS.start_timer("advance_time");
    err+= advance_time();
    //cout <<"advance_time took "<<GS.stop_timer("advance_time")<<" secs.\n";
    if (err!=0){cerr<<"(TIME_INT::advance_time) err!=0 Something went bad\n";return(1);}

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
    tsf=GS.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err+= output_data();
    if (err!=0){cerr<<"(TIME_INT::output_data) err!=0 Something went bad\n";return(1);}

    err+= check_eosim();
    if (err!=0){cerr<<"(TIME_INT::) err!=0 Something went bad\n";return(1);}
  }

  cout <<"(IntUniformFV::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

  tsf=GS.time_so_far("Time_Int");
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
void IntUniformFV::calculate_magnetic_pressure()
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
void IntUniformFV::calculate_blastwave_radius()
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




int IntUniformFV::calc_timestep()
{
  //
  // This is now basically a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt();
  t_mp  = calc_microphysics_dt();
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
#ifndef NEW_TIME_UPDATE


#ifdef THERMAL_CONDUCTION
double IntUniformFV::calc_conduction_dt_and_Edot()
{
  //
  // First we need to set Edot() in every cell.  This is stored in
  // c->dU[ERG] since this is where it will be updated.
  //
  //cout <<"\tCCdt: setting Edot.\n";
  eqn->set_thermal_conduction_Edot();
  //cout <<"\tCCdt: Done with div(Q).  Now getting timestep.\n";

  //
  // Now run through the grid and calculate E_int/Edot to get the smallest
  // timestep we can allow.  Also reset dU[RHO] to zero (was used as temporary
  // variable to hold gas temperature).
  // Note that we can allow a huge amount of energy into a cell, we just don't 
  // want more energy to leave than is already in it.  So we only count a cell
  // in the timestep limiting if Edot is negative.
  //
  double dt=1.0e200, gm1=SimPM.gamma-1.0, tempdt=0.0, minP=SimPM.RefVec[PG]*1.0e-3;
  cell *c = grid->FirstPt();
  do {
    // DIDN'T WORK -- CHECKERBOARDING!
    //if (c->dU[ERG]<0.0) tempdt = c->Ph[PG]/(gm1*(fabs(c->dU[ERG] +TINYVALUE)));
    //else                tempdt = 1.0e200;
    if (fabs(c->Ph[PG]>minP)) tempdt = c->Ph[PG]/(gm1*(fabs(c->dU[ERG] +TINYVALUE)));
    else                tempdt = 1.0e200;
    //if (c->id==39783)  {
    //  cout <<"\tCCdt: id="<<c->id<<", Edot="<<c->dU[ERG]<<", E="<<c->Ph[PG]/gm1;
    //  cout <<", this dt="<<tempdt<<",  min="<<dt<<"\n";
    //}
    dt = min(dt,tempdt);
  } while ((c=grid->NextPt(c)) !=0);
  //cout <<"Final conduction dt="<<dt<<", to be multiplied by 0.3.\n";

  if (SimPM.timestep==0) dt=1.0e3; // first timestep can be dodgy with conduction and stellar winds.

  return 0.1*dt;
} // calc_conduction_dt_and_Edot()
#endif // THERMAL CONDUCTION


// ##################################################################
// ##################################################################


void IntUniformFV::timestep_checking_and_limiting()
{
  //
  // If the timestep is less than the minimum allowed, report an error
  // and bug out!  This must be before we check for the next output
  // time because that can limit the step arbitrarily.
  //
  if (SimPM.dt < SimPM.min_timestep) {
    ostringstream temp; temp.str("");
    temp <<"Timestep too short! dt="<<SimPM.dt<<"  min-step="<<SimPM.min_timestep;
    rep.error(temp.str(),SimPM.dt);
  }

  //
  // So that we don't increase the timestep by more than 30% over last step:
  //
#ifdef TIMESTEP_LIMITING
  if (SimPM.dt > 1.3*SimPM.last_dt)
    cout <<"limiting step from "<<SimPM.dt<<" to "<<1.3*SimPM.last_dt<<"\n";
  SimPM.dt = min(SimPM.dt,1.3*SimPM.last_dt);
#endif // TIMESTEP_LIMITING

  //TESTING RT
  //SimPM.dt = 0.1*SimPM.finishtime*SimPM.CFL;
  // so we can do between 10 and 10^4 timesteps per sim.
  //TESTING RT

  //
  // If we are outputting every n-years, then check if we need to adjust dt.
  //
  if (SimPM.op_criterion==1) {
    SimPM.dt = min(SimPM.dt, SimPM.next_optime-SimPM.simtime);
    if (SimPM.dt <= 0.0)
      rep.error("Went past output time without outputting!",SimPM.dt);
  }

  //
  // Make sure we end up exactly at finishtime:
  //
  SimPM.dt = min(SimPM.dt, SimPM.finishtime-SimPM.simtime);
  if (SimPM.dt <= 0.0)
    rep.error("Went past Finish time without Stopping!",SimPM.dt);

  return;
}


// ##################################################################
// ##################################################################


double IntUniformFV::calc_dynamics_dt()
{
  double tempdt=0.0;
  double dt=1.e100; // Set it to very large no. initially.

  class cell *c=grid->FirstPt();
#ifdef TESTING
  dp.c = c;
#endif

  //
  // Some simulations have rapid inflow which is not on the domain for
  // the first timestep, so we calculate a timestep for the boundary
  // cell also, in the -x direction only! So make sure the inflow is
  // always in the x-direction from the left.  We do this for the
  // first few steps only, and only if doing a stellar-jet or
  // stellar-wind simulation.
  //
  // [NOTE: I WOULD LIKE TO MAKE THIS BETTER.  THERE SHOULD BE CORNER
  // BOUNDARY CELLS SO I SHOULD BE ABLE TO TRACE OUT THE WHOLE GRID
  // INCLUDING BOUNDARY DATA IN ONE GO.  IT'S ON THE TO-DO LIST!]

  //
  // First calculate for a boundary point (in case of incoming jet or
  // stellar wind!).
  //
#ifndef RT_TEST_PROBS
  if (SimPM.timestep<=10) {
     c = grid->NextPt(c,XN);
     tempdt = eqn->CellTimeStep(c,SimPM.gamma,SimPM.dx);
     c = grid->NextPt(c,XP);
#ifdef TESTING
     cout <<"\tBoundary point timestep! ";
#endif
     dt = min(dt,tempdt);
#ifdef TESTING
     cout <<"\tdt = "<<tempdt<<"\n";  
#endif
  }
#endif // not RT_TEST_PROBS

  //
  // Now go through all of the cells on the local grid.
  // The CellTimeStep() function returns a value which is
  // already multiplied by the CFL coefficient.
  //
  do {
#ifdef TESTING
    dp.c = c;
#endif
    tempdt = eqn->CellTimeStep(c, ///< pointer to cell
             SimPM.gamma, ///< gas EOS gamma.
             SimPM.dx  ///< Cell size dx.
             );
    if(tempdt<=0.0)
      rep.error("CellTimeStep function returned failing value",c->id);
    //    commandline.console("timestep -> ");
    dt = min(dt, tempdt);
    //cout <<"(get_min_timestep) i ="<<i<<"  min-dt="<<mindt<<"\n";    
  } while ( (c =grid->NextPt(c)) !=0);
  if (dt <= 0.0)
    rep.error("Got zero timestep!!!",dt);
  //  cout <<"(get_min_timestep)  min-dt="<<dt<<"\n";

  return dt;
}


// ##################################################################
// ##################################################################


double IntUniformFV::calc_microphysics_dt()
{
  //
  // If we have microphysics, we may want to limit the timestep by
  // cooling/heating/chemistry timescales, so we do that here.
  //
  // First check that we are doing microphysics, and if not return a very
  // large timestep so it cannot be the limiting factor.
  //
  if (!MP) {
    return 1.0e99; 
  }
  //
  // Now we know we have microphysics, so check that MP_timestep_limit is set
  // to a non-zero value.
  //
  if (SimPM.EP.MP_timestep_limit==0) {
    return 1.0e99;
  }

  double dt = -1.0e99; // initialise to negative value

  //
  // We see if we need to use the column densities and source properties in the 
  // timestep calculation, or if we just use the local microphysical cooling
  // and/or recombination times.
  //
  if (FVI_need_column_densities_4dt) {
    //
    // need column densities, so do raytracing, and then get dt.
    //
    //cout <<"calc_timestep, getting column densities.\n";
    int err = calculate_raytracing_column_densities();
    if (err) rep.error("calc_MP_dt: bad return value from calc_rt_cols()",err);
    dt = get_mp_timescales_with_radiation();
    if (dt<=0.0)
      rep.error("get_mp_timescales_with_radiation() returned error",dt);
  }
  else {
    //
    // don't need column densities, so call the no-RT version
    //
    //cout <<" getting timestep with no radiation\n";
    dt = get_mp_timescales_no_radiation();
    if (dt<=0.0)
      rep.error("get_mp_timescales_no_radiation() returned error",dt);
  }

#ifdef TESTING
  cout <<"(calc_microphysics_dt)  min-dt="<<dt<<"\n";
#endif

  return dt;
}


// ##################################################################
// ##################################################################


double IntUniformFV::get_mp_timescales_no_radiation()
{
#ifdef TESTING
  //
  // paranoid checking...
  //
  if (SimPM.EP.MP_timestep_limit==0) {
    cout <<"IntUniformFV::get_mp_timescales_no_radiation() called, but no MP-dt limiting!\n";
    return -1.0;
  }
#endif // TESTING
  //
  // So now we know we need to go through every cell and see what the limit is.
  //
  double tempdt=0.0, dt=1.0e99;
  class cell *c=grid->FirstPt();
  do {
#ifdef TESTING
    dp.c = c;
#endif
    //
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data)  If it is
    // boundary data then we skip it.

    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in get_mp_timescales_no_radiation() c->isbd.\n";
#endif
    }
    else {
      //
      // timescales(state-vec, gamma, t_cool?, t_rec?, t_photoion?)
      //
      switch (SimPM.EP.MP_timestep_limit) {
      case 1: // cooling time only.
        tempdt = MP->timescales(c->Ph, SimPM.gamma, true, false, false);
        break;
      case 4: // recomb only
        tempdt = MP->timescales(c->Ph, SimPM.gamma, false, true, false);
        break;
      case 2: // cooling+recomb
        tempdt = MP->timescales(c->Ph, SimPM.gamma, true,  true, false);
        break;
      case 3: // cooling+recomb+ionisation (not working!)
        tempdt = MP->timescales(c->Ph, SimPM.gamma, true,  true, true);
        break;
      default:
        rep.error("Bad MP_timestep_limit",SimPM.EP.MP_timestep_limit);
      }
      dt = min(dt, tempdt);
      //cout <<"(get_min_timestep) i ="<<i<<"  min-dt="<<dt<<"\n";
    } // if not boundary data.
  } while ( (c =grid->NextPt(c)) !=0);

#ifndef RT_TEST_PROBS
  //
  // If doing photo-ionisation, can underestimate 1st timestep because gas is
  // all neutral initially, but only for a few seconds!!!
  // (DON'T WANT TO SET THIS FOR NON-DYNAMICS TEST PROBLEMS)
  //
  if (SimPM.timestep<3 && (RT) && SimPM.EP.phot_ionisation) {
    //
    // adjust first timestep so that it corresponds to ionised cell.
    // 
    dt = min(dt, grid->DX()*SimPM.CFL/2.0e6); // 20 km/s wavespeed.
    //
    // Also make sure it is not larger than 1.0e6 seconds (0.03 years).
    // With photoionisation we must be using CGS units, so it is ok to 
    // hardcode the time.
    //
    dt = min(dt,1.0e7);
    cout <<"\tRT timestep: \t\t\tdt="<<dt<<"\n";
  }
  //
  // Make sure the first timestep is short if doing RT, here set to be 
  // 0.3333* the recombination time.
  // approx = 0.33333/(2.59e-13*rho/2.338e-24) = 3.009e-12/rho
  // Since we must have microphysics set up here, it is safe to assume
  // the code is using cgs units for density.
  //
  if ((SimPM.timestep==0) && (RT)) {
    c=grid->FirstPt();
    cout <<"rho="<<c->Ph[RO]<<", old dt="<<dt;
    dt = min(dt, 3.009e-12/c->Ph[RO]);
    cout <<", updated dt="<<dt<<"\n";
  }
#endif // RT_TEST_PROBS

  return dt;
}


// ##################################################################
// ##################################################################


double IntUniformFV::get_mp_timescales_with_radiation()
{
#ifdef TESTING
  //
  // paranoid checking...
  //
  if (SimPM.EP.MP_timestep_limit==0) {
    cout <<"IntUniformFV::get_mp_timescales_with_radiation() called, but no MP-dt limiting!\n";
    return -1.0;
  }
  if (!RT) rep.error("Called IntUniformFV::get_mp_timescales_with_radiation() but RT=0",1);
#endif // TESTING

  //
  // RT source properties are already in structs for the microphysics calls.
  // So now we need to go through every cell and see what the limit is.
  //
  double tempdt=0.0, dt=1.0e99;
  class cell *c=grid->FirstPt();
  do {
#ifdef TESTING
    dp.c = c;
#endif
    //
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data then we skip it.
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in get_mp_timescales_with_radiation() c->isbd.\n";
#endif
    }
    else {
      //
      // Get column densities and Vshell in struct for each source.
      //
      for (int v=0; v<FVI_nheat; v++) {
        FVI_heating_srcs[v].Vshell  = CI.get_cell_Vshell(c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].dS      = CI.get_cell_deltaS(c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].DelCol  = CI.get_cell_col(   c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].Column  = CI.get_col(c, FVI_heating_srcs[v].id) -FVI_heating_srcs[v].DelCol;
      }
      for (int v=0; v<FVI_nion; v++) {
        FVI_ionising_srcs[v].Vshell = CI.get_cell_Vshell(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].dS     = CI.get_cell_deltaS(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].DelCol = CI.get_cell_col(   c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].Column = CI.get_col(c, FVI_ionising_srcs[v].id)-FVI_ionising_srcs[v].DelCol;
      }
      //
      // For the new update we assume we want to limit by all relevant timescales.
      //
      tempdt = MP->timescales_RT(c->Ph, FVI_nheat, FVI_heating_srcs, FVI_nion, FVI_ionising_srcs, SimPM.gamma);
#ifdef TESTING
      //if (tempdt<dt) {
      //  cout <<"(get_min_timestep) id="<<c->id<<":  dt="<<tempdt<<", min-dt="<<dt;
      //  cout <<".\t 1-x="<<1.0-c->Ph[SimPM.ftr]<<", pg="<<c->Ph[PG]<<"\n";
      //}
#endif
      if (tempdt<=0.0) {
        cout <<"get_mp_timescales_with_radiation() negative timestep... ";
        cout <<"c->id="<<c->id<<"\tdt="<<tempdt<<"\n";
        rep.printVec("Ph",c->Ph,SimPM.nvar);
        rep.error("Negative timestep from microphysics with RT!",tempdt);
      }
      dt = min(dt, tempdt);
    } // if not boundary data.
  } while ( (c =grid->NextPt(c)) !=0);

  return dt;
}


// ##################################################################
// ##################################################################


int IntUniformFV::advance_time()
{
  int err=0;

  //
  // The type of Raytracing update determines what kind of time-update we need
  // to do.
  //
  if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT) {
#ifdef RT_TESTING
    cout <<"--- timestep_dynamics_then_microphysics() selected for implicit RT.\n";
#endif
    err += timestep_dynamics_then_microphysics();
    if (err)
      rep.error("advance_time_dynamics_then_microphysics() returned error",err);
  }
  else {
#ifdef TEST_SECOND_ORDER
    err += timestep_mp_dyn_second(); 
    //cout <<"Testing 2nd order MP update\n";
#else
    err += timestep_microphysics_then_dynamics();
#endif
    if (err)
      rep.error("advance_time_microphysics_then_dynamics() returned error",err);
  }

  return 0;
}


#ifdef TEST_SECOND_ORDER
int IntUniformFV::timestep_mp_dyn_second()
{
#ifdef TESTING
  cout <<"Using  IntUniformFV::timestep_mp_dyn_second() update.\n";
#endif // TESTING

  //
  // This algorithm is for small-timestep microphysics, with raytracing.
  //
  int err=0;
  //err += calculate_raytracing_column_densities();
  //if (err) 
  //  rep.error("Step_mp_dyn_second: bad return value from first calc_rt_cols()",err);

  //
  // First we calculate the timestep.  The microphysics timescales may depend
  // on the ray-tracing column densities, and if so all the column densities
  // will be calculated with raytracing calls in calc_mp_timestep()
  //
  err += calc_timestep();
  if (err) 
    rep.error("Step_mp_dyn_second: bad return value from calc_timestep()",err);

  if      (SimPM.tmOOA==OA1 && SimPM.spOOA==OA1 ) {
    eqn->Setdt(SimPM.dt);
    if (!FVI_need_column_densities_4dt) {
      err += calculate_raytracing_column_densities();
      if (err) 
        rep.error("Step_mp_dyn_second: bad return value from O1 calc_rt_cols()",err);
    }
    err += update_microphysics(SimPM.dt, OA1, OA1);
    err += grid->TimeUpdateInternalBCs(OA1,OA1);
    err += grid->TimeUpdateExternalBCs(OA1,OA1);
    if (err) 
      rep.error("Step_mp_dyn_second: bad return value from update_MP(half-step)",err);
    err  = update_dynamics(OA1,OA1); // space and time accuracy
    err += grid->TimeUpdateInternalBCs(OA2,OA2);
    err += grid->TimeUpdateExternalBCs(OA1,OA2);
    if (err) rep.error("O1 time update loop generated errors",err);
  }

  else {
    //
    // Now we do the 2nd order dynamics update.
    //cout <<"second order!!!\n";
    // Update microphysics to the half step.
    // Check that calc_mp_dt() did the raytracing already.
    //
    SimPM.dt /=2.;
    eqn->Setdt(SimPM.dt);
    if (!FVI_need_column_densities_4dt) {
      err += calculate_raytracing_column_densities();
      if (err) 
        rep.error("Step_mp_dyn_second: bad return value from first calc_rt_cols()",err);
    }
    err += update_microphysics(SimPM.dt, OA1, OA2);
    err += grid->TimeUpdateInternalBCs(OA2,OA2);
    err += grid->TimeUpdateExternalBCs(OA1,OA2);
    if (err) 
      rep.error("Step_mp_dyn_second: bad return value from update_MP(half-step)",err);
    //
    // Now do the full-step microphysics.
    // In Mackey (2012) this is in the middle of the dynamics update, but that
    // introduces some instabilities in the flow, so it is better to keep it
    // outside.
    //
    SimPM.dt *=2.;
    eqn->Setdt(SimPM.dt);
    err += calculate_raytracing_column_densities();
    if (err)
      rep.error("Step_mp_dyn_second: bad return value from second calc_rt_cols()",err);
    err += update_microphysics(SimPM.dt, OA2, OA2);
    err += grid->TimeUpdateInternalBCs(OA2,OA2);
    err += grid->TimeUpdateExternalBCs(OA2,OA2);
    if (err) 
      rep.error("Step_mp_dyn_second: bad return value from update_MP(half-step)",err);

    //
    // Update dynamics for the half step.
    //
    SimPM.dt /=2.;
    eqn->Setdt(SimPM.dt);
    err  = update_dynamics(OA1,OA1); // space and time accuracy
    err += grid->TimeUpdateInternalBCs(OA2,OA2);
    err += grid->TimeUpdateExternalBCs(OA1,OA2);
    if (err) rep.error("O2 half time update loop generated errors",err);
    //
    // Now do the update again, for the full timestep.
    //
    SimPM.dt *=2.;
    eqn->Setdt(SimPM.dt);
    err = update_dynamics(OA2,OA2); // space and time accuracy
    if (err) rep.error("O2 dynamics time update loop generated errors",err);
    err += grid->TimeUpdateInternalBCs(OA2,OA2);
    err += grid->TimeUpdateExternalBCs(OA2,OA2);
    if (err) rep.error("O2 full time update-boundaries generated errors",err);
  }
 
#ifdef TESTING
  if (SimPM.timestep%20 ==0) {
    //    cout <<"dp.initERG = "<<dp.initERG<<"\n";
    check_energy_cons();
  }
#endif // TESTING
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime +=SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;
  return(0);
}

#endif // TEST_SECOND_ORDER




int IntUniformFV::timestep_microphysics_then_dynamics()
{
#ifdef TESTING
  cout <<"Using  IntUniformFV::timestep_microphysics_then_dynamics() update.\n";
#endif // TESTING
  cout <<"***WARNING***, USING 2ND ORDER UPDATE WHERE MICROPHYSICS IS NOT UPDATED";
  cout <<" IN THE HALF STEP!  THIS IS NOT RECOMMENDED!!\n";
  
  //
  // This algorithm is for small-timestep microphysics, with raytracing.
  //
  int err=0;

  //
  // First we calculate the timestep.  If the microphysics timescales depend
  // on the ray-tracing column densities, then they will be calculated within
  // this call.
  //
  err += calc_timestep();
  if (err) 
    rep.error("Step_MP-then-DYN: bad return value from calc_timestep()",err);

  //
  // Now we update microphysics based on the timestep just calculated.
  // We also need to update internal and external boundary data to account
  // for changes on the grid-data.
  // Also check that calc_mp_dt() did the raytracing already.
  //
  if (!FVI_need_column_densities_4dt) {
    err += calculate_raytracing_column_densities();
    if (err) 
      rep.error("Step_mp_dyn_second: bad return value from first calc_rt_cols()",err);
  }
  err += update_microphysics(SimPM.dt, SimPM.tmOOA, SimPM.tmOOA);
  err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
  err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
  if (err) 
    rep.error("Step_MP-then-DYN: bad return value from update_MP()",err);

  //
  // Now we do the dynamics update, either 1st or 2nd order.  This code is 
  // mostly ripped from timestep_dynamics_then_microphysics().
  // 
  int sp = OA1;
  int tm = OA1;
  if (SimPM.tmOOA==OA1 && SimPM.spOOA==OA2) sp = OA2;

  if (SimPM.tmOOA ==OA1) {
    //
    // First order time accurate, so just a single step update.
    //
    eqn->Setdt(SimPM.dt);
    err  = update_dynamics(sp,tm);
    // Update Boundary Values
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    if (err) rep.error("O1 time update loop generated errors",err);
  } // if 1st order accurate
 
  else if (SimPM.tmOOA==OA2) {
    //cout <<"second order!!!\n";
    //
    // This is the Falle, Komissarov, Joarder (1998) 2nd order (time
    // and space) accurate time update.
    // First we halve the timestep to get the half-step values.
    //
    SimPM.dt /=2.;
    eqn->Setdt(SimPM.dt);
    err  = update_dynamics(sp,tm);
    // Update Boundary Values
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(tm,SimPM.tmOOA);
    if (err) rep.error("O2 half time update loop generated errors",err);
    
    //
    // Now do the update again, for the full timestep.
    //
    SimPM.dt *=2.;
    eqn->Setdt(SimPM.dt);
    //    cout <<"\tsecond pass.... \n";
    err = update_dynamics(SimPM.spOOA,SimPM.tmOOA);
#ifdef TESTING
    if (err) rep.error("O2 dynamics time update loop generated errors",err);
#endif // TESTING
    // Update Boundary Values
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    if (err) rep.error("O2 full time update loop generated errors",err);
  } // If 2nd order accurate.
  
  else rep.error("Only know first and second order accuracy",SimPM.tmOOA);
 
#ifdef TESTING
  if (SimPM.timestep%20 ==0) {
    //    cout <<"dp.initERG = "<<dp.initERG<<"\n";
    check_energy_cons();
  }
#endif // TESTING
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime +=SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;
  return(0);
}


int IntUniformFV::timestep_dynamics_then_microphysics()
{
#ifdef TESTING
  cout <<"Using  IntUniformFV::timestep_dynamics_then_microphysics() update.\n";
#endif // TESTING

  // This is a general conservative scheme,
  // with all the details in the calculation of the intercell fluxes.

  /** \section Boundaries
   * The boundary cells are updated by calling the function 
   * grid->TimeUpdate[Internal/External]BCs().  This function knows how each
   * boundary is to be updated.
   * */
  
  /** \section Description
   * This has to handle four cases so far, Space,Time accuracy of 
   * {[1,1], [1,2], [2,1], [2,2]}, leaving aside that some of these
   * may be unstable.
   * 
   * If time accuracy =1, then I only call dU once, and update the 
   * state vector directly.\n
   * If time accuracy =2, then I call dU twice, and update an intermediate
   * state vector the first time, use the intermediate state for the
   * second calculation of dU, and then update the main state vector.
   * If I ever put in time acc. =3, then I'll need to do two intermediate
   * steps and then a final update.  Given that it's unlikely I'll ever do 
   * more than this, I think it's safe to not do anything too fancy with this
   * algorithm.
   * 
   * */

  //
  // First get the timestep
  //
  int err=0;
  err += calc_timestep();
  if (err!=0) {
    cerr<<"(advance_time->calc_timestep) err!=0 Something went bad\n";
    return err;
  }

  int sp = OA1;
  int tm = OA1;
  if (SimPM.tmOOA==OA1 && SimPM.spOOA==OA2) sp = OA2;
  //  cout <<"dt = "<<SimPM.dt<<"\n";
  
  if (SimPM.tmOOA ==OA1) { // First order time time
    eqn->Setdt(SimPM.dt);
    err  = update_dynamics(sp,tm);
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    //     cout <<"updating microphysics.\n";
    err += update_microphysics(SimPM.dt, tm, SimPM.tmOOA);
    //    cout <<"done with mp.\n";
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    //     cout <<"updating external bcs.\n";
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    //    cout <<"done with external bcs.\n";
    if (err) rep.error("O1 time update loop generated errors",err);
  } // if 1st order accurate
  
  else if (SimPM.tmOOA==OA2) {
    //cout <<"second order!!!\n";
    //
    // This is the Falle, Komissarov, Joarder (1998) 2nd order (time
    // and space) accurate time update.
    //

    //
    // Halve the timestep for the first pass
    //
    SimPM.dt /=2.;
    eqn->Setdt(SimPM.dt);
    err  = update_dynamics(sp,tm);
    // Update Boundary Values
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    // Update MicroPhysics, if present (not on 1st half step though).
    //err += update_microphysics(SimPM.dt, SimPM.tmOOA, SimPM.tmOOA);
    // Update Boundary Values
    //err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(tm,SimPM.tmOOA);
    if (err) rep.error("O2 half time update loop generated errors",err);
    
    //
    // Now calculate dU again, for the full timestep.
    //
    SimPM.dt *=2.;
    eqn->Setdt(SimPM.dt);
    //    cout <<"\tsecond pass.... \n";
    err = update_dynamics(SimPM.spOOA,SimPM.tmOOA);
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    // Update MicroPhysics, if present
    err += update_microphysics(SimPM.dt, SimPM.tmOOA, SimPM.tmOOA);
    //  output_data();
    //  rep.error("finished.",0);
    // Update Boundary Values
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    if (err) rep.error("O2 full time update loop generated errors",err);
  } // If 2nd order accurate.
  
  else rep.error("Only know first and second order accuracy",SimPM.tmOOA);
  
#ifdef TESTING
  if (SimPM.timestep%20 ==0) {
    //    cout <<"dp.initERG = "<<dp.initERG<<"\n";
    check_energy_cons();
  }
#endif // TESTING
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime +=SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;
  return(0);
}



// ##################################################################
// ##################################################################



int IntUniformFV::calculate_raytracing_column_densities()
{
  int err=0;
  //
  // If we have raytracing, we call the new ray-tracing function 
  // to get Tau0, dTau, Vshell in cell->extra_data[n].
  //
  if (RT) {
    for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
#ifdef RT_TESTING
      cout <<"IntUniformFV::calc_RT_col_dens: SRC-ID: "<<isrc<<"\n";
#endif
      err += RT->RayTrace_Column_Density(isrc, 0.0, SimPM.gamma);
      if (err) {
        cout <<"isrc="<<isrc<<"\t"; 
        rep.error("calc_RT_col_dens step in returned error",err);
      } // if error
    } // loop over sources
  }
  return err;
}


// ##################################################################
// ##################################################################



int IntUniformFV::update_microphysics(
            const double delt, // timestep to integrate
            const int ctmoa,
            const int ttmoa
            )
{
  //cout <<"\tupdate_microphysics starting.\n";

  //
  // first if we aren't doing microphysics or raytracing, the modules will not be
  // initialised, so just return zero.  RT requires MP, so just check MP.
  //
  if (!MP) return 0;
  
  // Secondly, if we are on a fractional timestep, we don't want to 
  // do the update. (this function shouldn't be called in this case,
  // but just to be sure...).  So just return.
  //if (ctmoa != ttmoa) {
  //  cout <<"Don't call microphysics on half-step update!\n";
  //  return 0;
  //}
  int step = 0;
  if (ctmoa==ttmoa) {
#ifdef TESTING
    cout <<"Full time update for microphysics\n";
#endif
    step=TIMESTEP_FULL;
  }
  else {
#ifdef TESTING
    cout <<"Half time update for microphysics\n";
#endif
    step=TIMESTEP_FIRST_PART;
  }
  
#ifdef TESTING
  cout <<"update_microphysics() Updating MicroPhysics (OOA="<<ctmoa<<").  RT-Nsrc="<<SimPM.RS.Nsources<<"\n";
#endif // TESTING
  int err = 0;

  if (!RT) {
    //
    // If no radiative transfer, then just do a simple MP update.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling update_microphysics_no_RT()\n";
#endif // RT_TESTING
    err += update_microphysics_no_RT(delt,step);
  }

  else if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT) {
    //
    // This is the C2-ray single-source update which I did for my thesis.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling update_microphysics_JMs_C2ray_RT()\n";
#endif // RT_TESTING
    err += update_microphysics_JMs_C2ray_RT(delt,step);
  }

  else {
    //
    // must have at least one radiation source, and the RT and microphysics
    // are updated separately, and we may have diffuse radiation, so we do the
    // new RT update:
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling update_microphysics_general_RT()\n";
#endif // RT_TESTING
    err += update_microphysics_general_RT(delt,step);
  }
    
  //cout <<"\tupdate_microphysics finished.\n";
  return err;
}



// ##################################################################
// ##################################################################



int IntUniformFV::update_microphysics_general_RT(const double delt, // timestep to integrate
                            const int step   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                                )
{
#ifdef RT_TESTING
  if (!RT) rep.error("Logic error: must have RT unless i'm an idiot","GENERAL-RT");
  if (SimPM.RS.Nsources<1)
    rep.error("Need at least one source for diffuse-RT update",SimPM.RS.Nsources);
#endif // RT_TESTING
  int err=0;
  //
  // Now update microphyiscs with new interface
  // RT source properties are already in structs for the microphysics calls.
  //

  //
  // Now do MP update.  Only new MP classes have the  XX_RTnew() defined,
  // so if we try to call this with old code then it should return with 
  // a non-zero error code.
  //
  cell *c = grid->FirstPt();
  double p[SimPM.nvar];
  double tt=0.; // temperature returned at end of microphysics step.
  do {
    //
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in update_microphysics_general_RT() c->isbd.\n";
#endif
    }
    else {
      //
      // Get column densities and Vshell in struct for each source.
      //
      for (int v=0; v<FVI_nheat; v++) {
        FVI_heating_srcs[v].Vshell  = CI.get_cell_Vshell(c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].dS      = CI.get_cell_deltaS(c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].DelCol  = CI.get_cell_col(   c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].Column  = CI.get_col(c, FVI_heating_srcs[v].id) -FVI_heating_srcs[v].DelCol;
#ifdef TESTING
        //cout <<"HEAT: Vs="<<FVI_heating_srcs[v].Vshell<<", dS="<<FVI_heating_srcs[v].dS<<", dC="<<FVI_heating_srcs[v].DelCol<<", Col="<<FVI_heating_srcs[v].Column<<"\n";
#endif
      }
      for (int v=0; v<FVI_nion; v++) {
        FVI_ionising_srcs[v].Vshell = CI.get_cell_Vshell(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].dS     = CI.get_cell_deltaS(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].DelCol = CI.get_cell_col(   c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].Column = CI.get_col(c, FVI_ionising_srcs[v].id) -FVI_ionising_srcs[v].DelCol;
      }
      //
      // integer 9th argument determines type of integration substepping:
      // 0 = adaptive RK5 Cash-Karp method (use this!).
      // 4th and 5th args are for unused ionising sources.
      //
      err += MP->TimeUpdateMP_RTnew(c->P, FVI_nheat, FVI_heating_srcs, FVI_nion, FVI_ionising_srcs,
                                    p, delt, SimPM.gamma, 0, &tt);
      //
      // update Ph (and P if full step)
      //
      if (step==TIMESTEP_FULL) {
        for (int v=0;v<SimPM.nvar;v++) c->P[v] = c->Ph[v] = p[v];
      }
      else if (step==TIMESTEP_FIRST_PART) {
        for (int v=0;v<SimPM.nvar;v++) c->Ph[v] = p[v];
      }
      else {
        rep.error("Bad step variable in microphysics partial update",step);
      }
    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"update_microphysics() Updating MicroPhysics Done!\n";
  return err;
} // newest general-RT microphysics update.



// ##################################################################
// ##################################################################



int IntUniformFV::update_microphysics_JMs_C2ray_RT(const double delt, // timestep to integrate
                            const int step   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                                )
{
#ifdef RT_TESTING
  if (!RT) rep.error("Logic error: must have RT unless i'm an idiot","C2RAY");
  if (step != TIMESTEP_FULL) {
    cout <<"C2Ray update only done once per step.  returning 0.\n";
    return 0;
  }
#endif // RT_TESTING

  int err=0;
  //
  // This is the C2Ray-style implicit radiative transfer update, which does the
  // microphysics update as the rays are traced outwards. 
  // We have to get the UV heating rates first, so if there are any UV-heating
  // sources we get the column densities for these sources first, and then call
  // the ionising source update.
  // RT source properties are already in structs for the microphysics calls.
  //
  //
  // Set column densities for UV-heating sources
  //
  int isrc=-1;
  for (int s=0; s<FVI_nheat; s++) {
    isrc = FVI_heating_srcs[s].id;
#ifdef RT_TESTING
    cout <<" -- update_mp_JMs_C2ray_RT(): Tracing UV-heating src: id="<<isrc<<"\n";
#endif
    err += RT->RayTrace_Column_Density(isrc, 0.0, SimPM.gamma);
    if (err) {
      cout <<"isrc="<<isrc<<"\t"; 
      rep.error("RT_col_dens step in returned error",err);
    }
  }
#ifdef RT_TESTING
  if (FVI_nion != 1) 
    rep.error("Can't do implicit update unless exactly one ionising src",FVI_nion);
  cout <<" -- update_mp_JMs_C2ray_RT(): Tracing ionising src: id="<<isrc<<".\n";
#endif
  err += RT->RayTrace_SingleSource(FVI_ionising_srcs[0].id, delt, SimPM.gamma);
  if (err) {
    cout <<"isrc="<<SimPM.RS.sources[0].id<<"\t"; 
    rep.error("JMs C2ray RT, returned error",err);
  } // if error


  return err;
}



// ##################################################################
// ##################################################################



int IntUniformFV::update_microphysics_no_RT(const double delt, // timestep to integrate
                            const int step   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                            )
{
  //if (SimPM.timestep <50) return 0;
#ifdef TESTING
  cout <<"update_microphysics_no_RT starting.\n";
#endif
  //
  // No radiation sources at all, and no diffuse radiation optical depths,
  // so just do my original microphysics update.
  //
  cell *c = grid->FirstPt();
  double p[SimPM.nvar];
  double tt=0.; // temperature returned at end of microphysics step.
  int err=0;
  do {
    //
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in update_microphysics_no_RT() c->isbd.\n";
#endif
    }
    else {
      //
      // integer 5th argument determines type of integration substepping:
      // 0 = adaptive RK5 Cash-Karp method.
      // 1 = dumb adaptive euler integration.
      // 2 = single step RK4 method (at your own risk!)
      //
      err += MP->TimeUpdateMP(c->P, p, delt, SimPM.gamma, 0, &tt);
      //err += MP->TimeUpdateMP(c->P, p, delt, SimPM.gamma, 2);
      //err += ump->TimeUpdateMP(c->P, p, delt, gamma, MP_ADAPTIVE);
      //rep.printVec("Original vector P",c->P,nvar);
      //rep.printVec("Updated  vector p",p   ,nvar);
      if (err) rep.error("update_microphysics_no_RT returned error: cell id",c->id);
      //
      // update Ph (and P if full step)
      //
      if (step==TIMESTEP_FULL) {
        for (int v=0;v<SimPM.nvar;v++) c->P[v] = c->Ph[v] = p[v];
      }
      else if (step==TIMESTEP_FIRST_PART) {
        for (int v=0;v<SimPM.nvar;v++) c->Ph[v] = p[v];
      }
      else {
        rep.error("Bad step variable in microphysics partial update",step);
      }
    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"update_microphysics() Updating MicroPhysics Done!\n";
  return err;
} 



// ##################################################################
// ##################################################################


  
int IntUniformFV::update_dynamics(const int csp, ///< spatial order of accuracy to use for update.
          const int ctm  ///< time order of accuracy to use for update.
          )
{
  //cout <<"\tupdate_dynamics starting.\n";
  //
  // first check if we are doing dynamics, and return if not.
  //
  if (!SimPM.EP.dynamics) return 0;
  int err=0;

#ifdef TESTING
  if (csp==OA1)
    cout <<"*****Updating dynamics: OA1\n";
  else if (csp==OA2)
    cout <<"*****Updating dynamics: OA2\n";
  else rep.error("Bad ooa",csp);
#endif //TESTING

  //
  // First we pre-process the cells, if needed.  This is required for
  // genuinely multi-dimensional viscosity such as Lapidus-like AV or
  // the H-Correction.
  //
  err = eqn->preprocess_data(csp,ctm);

  //
  // Now calculate the directionally-unsplit time update for the
  // conserved variables:
  //
  err = calc_dU(csp,ctm);
  rep.errorTest("IntUniformFV::update_dynamics() eqn->calc_dU returned error.",0,err);

  //
  // Post-processing is for if we are doing something like Constrained
  // Transport, where we have to change the B-field update.  At the
  // moment there is *NO* solver which does anything here since I
  // found the Dedner et al. (2002) divergence cleaning to be more
  // robust than e.g. Toth (2000) Field-CT method.  (well the internal
  // energy solver uses it, but it's not really worth using).
  //
  err = eqn->PostProcess_dU(csp,ctm);
  rep.errorTest("IntUniformFV::update_dynamics() eqn->PostProcess_dU() returned error.",0,err);


  double temperg =0.; // temp variable to handle change of energy when correcting for negative pressure.
  class cell* c = grid->FirstPt();
  do {

#ifdef TESTING
     dp.ergTotChange = 0.;temperg =0.;
     dp.c = c;
#endif
     err += eqn->CellAdvanceTime(c,
         c->P,        // Initial Primitive State Vector.
         c->dU,       // Update vector dU
         c->Ph,       // Final Primitive state vector (can be same as initial vec.).
         &temperg,    // Tracks change of energy if I have to correct for negative pressure
         SimPM.gamma, // gas EOS gamma.
         SimPM.dt     // Cell timestep dt.
         );
#ifdef TESTING
    if (err) {
      cout <<"______ Error in Cell-advance-time: ";
      CI.print_cell(c);
      err=0;
    }
#endif // TESTING
     if (ctm==SimPM.tmOOA) {
       for (int v=0;v<SimPM.nvar;v++) c->P[v] = c->Ph[v]; // set new-time value to real value.
#ifdef TESTING
       // Update Total Energy from fixing negative pressures. Reset
       // update variables.
       dp.ergTotChange = temperg;
       dp.initERG += dp.ergTotChange*eqn->CellVolume(c);
#endif // TESTING
     }
  } while ( (c =grid->NextPt(c)) !=0);

#ifdef TESTING
  cout <<"\tupdate_dynamics done. error="<<err<<"\n";
#endif // TESTING
  return err; 
}



// ##################################################################
// ##################################################################


int IntUniformFV::calc_dU(const int csp, const int ctm)
{
  //  cout <<"\t\t\tCALC_DU: Starting calc_dU: ndim = "<<SimPM.ndim<<"\n";
  int retval=0;

  // 
  // Allocate arrays
  //
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP; posdirs[1] = YP; posdirs[2] = ZP;
  negdirs[0] = XN; negdirs[1] = YN; negdirs[2] = ZN;
  axis[0] = XX; axis[1] = YY; axis[2] = ZZ;


  //
  // Loop over all directions, and in each direction, calculate fluxes
  // in all columns of cells in that direction (it does work!).
  // This function depends on cells being labelled as on-grid or as
  // boundary cells, and also on dU_column returning 0 on successful
  // completion, and -1 if the column ends up at the last cell in the domain.
  // Any other return value will signal an error in this function, stopping the code.
  //
  // 2011.04.29 JM: changed logic here so we check for cells which are not grid
  // cells.  Checking for boundary data is not correct, since we can have 
  // internal boundaries which are also grid data.
  //
  for (int i=0;i<SimPM.ndim;i++) {
    //    cout <<"\t\t\tidim="<<i<<"\n";
    eqn->SetDirection(axis[i]);
    class cell *cpt    = grid->FirstPt();
    class cell *marker = grid->FirstPt();
    
    //    rep.printVec("cpt",cpt->x,SimPM.ndim);
    //    rep.printVec("+XX",(grid->NextPt(cpt,XP))->x,SimPM.ndim);
    //    rep.printVec("+YY",(grid->NextPt(cpt,YP))->x,SimPM.ndim);
    //    rep.printVec("+ZZ",(grid->NextPt(cpt,ZP))->x,SimPM.ndim);
    
    while ( (retval = dU_column(cpt,posdirs[i],negdirs[i], csp, ctm)) ==0) {
      //       cout <<"next dir= "<<(i+1)%SimPM.ndim<<"\n";
      if ( !(cpt=grid->NextPt(cpt,posdirs[(i+1)%SimPM.ndim]))->isgd ) {
        //   cout <<"next dir= "<<(i+2)%SimPM.ndim<<"\n";
        if ( !(cpt=grid->NextPt(marker,posdirs[(i+2)%SimPM.ndim]))->isgd ) {
          rep.error("Got to edge of box before last point!",cpt);
        }
        marker = cpt;
      } // if null pointer.
    } // Loop over columns.
    if (retval!=-1) rep.error("dUdtColumn returned abnormally.",retval);
  } // Loop over three directions.
  eqn->SetDirection(axis[0]); // Reset fluxes to x-dir, (just to be safe!).
  //  cout <<"\t\t\tCALC_DU done.\n";

  //
  // all done, so return
  //
  return 0;
}   // calc_dU()




// ##################################################################
// ##################################################################


int IntUniformFV::dU_column
        (
        const class cell *startingPt, 
        const enum direction posdir, 
        const enum direction negdir,
        const int csp,
        const int ctm
        )
{
  //  cout <<"Starting dU_column in direction "<<posdir<<" at ";
  //  rep.printVec("starting position",startingPt->x,SimPM.ndim);
  if ( (SimPM.spOOA>2) || (SimPM.tmOOA>2) || (csp>2) || (ctm>2) ) {
    cerr<<"(RSMethod::calc_dUdt) Error, only know 1st and 2nd order accurate methods.\n";
    return(1);
  }
  int err = 0;
#ifdef TESTING
  int ct=0;
#endif
  enum axes axis = eqn->GetDirection();

  //
  // Run through all cells in grid, to calculate (up to) second order time and 
  // space update.
  //

  // Calculate Flux at positive (right) boundary of cell and store in temporary arrays
  // for the current cell (Fr_this) and the negative neighbour (Fr_prev)
  // Have to do it this way b/c ISO C++ forbids re-assignment of arrays.
  double *Fr_this=0, *Fr_prev=0, *temp=0, *slope_cpt=0, *slope_npt=0, *edgeR=0, *edgeL=0, *pstar=0;
  Fr_prev   = mem.myalloc(Fr_prev,   SimPM.nvar);
  Fr_this   = mem.myalloc(Fr_this,   SimPM.nvar);
  slope_cpt = mem.myalloc(slope_cpt, SimPM.nvar);
  slope_npt = mem.myalloc(slope_npt, SimPM.nvar);
  edgeL     = mem.myalloc(edgeL,     SimPM.nvar);
  edgeR     = mem.myalloc(edgeR,     SimPM.nvar);
  pstar     = mem.myalloc(pstar,     SimPM.nvar);

  //
  // Set starting point, and next two points in the column.
  //
  cell *cpt = grid->NextPt(startingPt,posdir); //=grid->NextPt(startingPt,negdir);
  while (grid->NextPt(cpt,negdir)) {cpt = grid->NextPt(cpt,negdir);} // grid->PrintCell(cpt);}
  if(cpt==0) {cerr<<"(RSMethod::calc_dUdt) error finding left boundary cell.\n";return(1);}
  cell *npt  = grid->NextPt(cpt,posdir);
  cell *n2pt = grid->NextPt(npt,posdir);
  if (npt==0 || n2pt==0) rep.error("Couldn't find two real cells in column",0);
  //  cout<<"First Cell:"; grid->PrintCell(cpt);
  //  cout<<"Next Cell: "; grid->PrintCell(npt);
  
  //
  // Left Ghost Cell (doesn't get updated)
  //
  for (int v=0;v<SimPM.nvar;v++) { 
    slope_cpt[v] = 0.;
    slope_npt[v] = 0.;
    Fr_prev[v]   = 0.;
    Fr_this[v]   = 0.;
    edgeL[v]     = 0.;
    edgeR[v]     = 0.;
  }
  
  //
  // Now go through all cells in the column and calculate fluxes and add to dU vector.
  //
  do {
#ifdef TESTING
    dp.c = cpt;
    //    if (SimPM.timestep==2959 && dp.c->id==93) commandline.console("-ve density> ");
#endif
    // Get the flux from left and right states, adding artificial viscosity if needed.
    err += eqn->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp);
    err += eqn->SetSlope(npt, axis, SimPM.nvar, slope_npt, csp);
    err += eqn->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp);
    //    rep.printVec("El",edgeL,SimPM.nvar); rep.printVec("Er",edgeR,SimPM.nvar);
    //    rep.errorTest("Edge States not obtained!",0,err);
    err += eqn->InterCellFlux(cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType, SimPM.artviscosity, SimPM.gamma, SimPM.dx);
    //    rep.printVec("Fr",Fr_this,SimPM.nvar);    rep.errorTest("Intercell Flux not obtained!",0,err);
    err += eqn->dU_Cell(cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, SimPM.dx, SimPM.dt);
    //    rep.errorTest("dU not obtained!",0,err);
#ifdef TESTING
    for (int v=0;v<SimPM.nvar;v++) {
      if(!isfinite(cpt->dU[v])) {
  rep.printVec("Fl",Fr_prev,SimPM.nvar); rep.printVec("Fr",Fr_this,SimPM.nvar);
  cout <<"dt:"<<SimPM.dt<<"\tdx="<<SimPM.dx<<"\n";
//  rep.printVec("dU",&cpt->dU[v],1);
  grid->PrintCell(cpt); grid->PrintCell(npt);
  rep.error("nans!!!",2);
      }
      //if (dp.c->id==114337) grid->PrintCell(cpt);
    }
    // Track energy, momentum entering domain.
    if(ctm==SimPM.tmOOA && !(cpt->isgd) && npt->isgd) {
      ct++; if (ct>1) rep.error("Entering domain more than once! (dUcolumn)",ct);
      //      cout <<"Entering Domain Dir = "<<posdir<<" and interface area = "<<eqn->CellInterface(cpt,posdir)<<"\n";
      dp.initERG += Fr_this[ERG]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMX += Fr_this[MMX]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMY += Fr_this[MMY]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMZ += Fr_this[MMZ]*SimPM.dt*eqn->CellInterface(cpt,posdir);
//      if (posdir==YP && fabs(Fr_this[MMY])>2*MACHINEACCURACY) {
//  cout <<"R-momentum flux entering domain from R=0(!) = "<<Fr_this[MMY]<<"\n";
//  cout <<"v_R in first cell = "<<npt->Ph[VY]<<", "<<npt->P[VY]<<"\n";
//      }
    }
    else if (ctm==SimPM.tmOOA && !(npt->isgd) && cpt->isgd) {
      ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
      //      cout <<"Leaving Domain Dir = "<<posdir<<" and interface area = "<<eqn->CellInterface(cpt,posdir)<<"\n";
      dp.initERG -= Fr_this[ERG]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMX -= Fr_this[MMX]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMY -= Fr_this[MMY]*SimPM.dt*eqn->CellInterface(cpt,posdir);
      dp.initMMZ -= Fr_this[MMZ]*SimPM.dt*eqn->CellInterface(cpt,posdir);
    }
//    if(posdir==YP && ctm==SimPM.tmOOA && cpt->x[YY]<SimPM.dx && cpt->x[YY]>0 && fabs(Fr_prev[MMY])>2*MACHINEACCURACY) {
//      cout <<"R-Momentum Flux leaving first cell at R=dR = "<<Fr_this[MMY]<<"\n";
//    }
#endif //TESTING
    //
    // Now move temp arrays to prepare for moving on to the next cell.
    //
    temp=Fr_prev;
    Fr_prev = Fr_this;
    Fr_this = temp; // just point to the free memory to be overwritten next step.
    temp = slope_cpt;
    slope_cpt = slope_npt;
    slope_npt = temp;
    cpt = npt; npt = n2pt;
  }  while ( (n2pt=grid->NextPt(n2pt,posdir)) );
  
  
  //
  // Now n2pt=null. npt = bd-data, cpt= (gd/bd)-data. So have to do something different.
  //
#ifdef TESTING
  dp.c = cpt;
#endif
  err += eqn->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp);
  for (int v=0;v<SimPM.nvar;v++) slope_npt[v] = 0.; // last cell must be 1st order.
  err += eqn->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp);
  err += eqn->InterCellFlux(cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType, SimPM.artviscosity, SimPM.gamma, SimPM.dx);
  err += eqn->dU_Cell(cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, SimPM.dx, SimPM.dt);
#ifdef TESTING
  if (ctm==SimPM.tmOOA && cpt->isgd && !(npt->isgd)) {
    ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
    // cout <<"Leaving Domain Dir = "<<posdir<<" and interface area = "<<eqn->CellInterface(cpt,posdir)<<"\n";
    dp.initERG -= Fr_this[ERG]*SimPM.dt*eqn->CellInterface(cpt,posdir);
    dp.initMMX -= Fr_this[MMX]*SimPM.dt*eqn->CellInterface(cpt,posdir);
    dp.initMMY -= Fr_this[MMY]*SimPM.dt*eqn->CellInterface(cpt,posdir);
    dp.initMMZ -= Fr_this[MMZ]*SimPM.dt*eqn->CellInterface(cpt,posdir);
  }
#endif //TESTING
 
  //
  // Right Ghost Cell-- have to calculate it's left interface differently,
  //
  cpt=npt;
#ifdef TESTING
  dp.c = cpt;
#endif
  for (int v=0;v<SimPM.nvar;v++) cpt->dU[v] +=0.; // nothing to calculate for it.


  if (err!=0) {cerr<<"(RSMethod::calc_dUdt) Riemann solver returned abnormally... exiting.\n";return(err);}
  
  Fr_this   = mem.myfree(Fr_this);
  Fr_prev   = mem.myfree(Fr_prev);
  slope_cpt = mem.myfree(slope_cpt);
  slope_npt = mem.myfree(slope_npt);
  edgeL     = mem.myfree(edgeL);
  edgeR     = mem.myfree(edgeR);
  pstar     = mem.myfree(pstar);

  //
  // Check if this is the last column or not. (first track back from bd to grid).
  // If it is the last column, return -1 instead of 0 to indicate this. (A positive
  // return value from errors is picked up in the calling function).
  //
  do{} while( (cpt=grid->NextPt(cpt,negdir))->isbd );
  if (cpt->id == grid->LastPt()->id) return(-1);
  else return(0);
}


// ##################################################################
// ##################################################################

#endif // not NEW_TIME_UPDATE

// ##################################################################
// ##################################################################



int IntUniformFV::check_eosim()
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



int IntUniformFV::check_energy_cons()
{
#ifdef TESTING
  // Energy, and Linear Momentum in x-direction.
  double u[SimPM.nvar];
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
int IntUniformFV::Finalise()
{
  int err=0;
  cout <<"------------------------------------------------------------\n";
  cout <<"(IntUniformFV::Finalise) FINALISING SIMULATION."<<"\n";
  err += check_energy_cons();
  err+= output_data();
  if (err!=0){cerr<<"(FINALISE::output_data) final state data output. err!=0 Something went bad"<<"\n";return(1);}
  cout <<"\tSimTime = "<<SimPM.simtime<<"   #timesteps = "<<SimPM.timestep<<"\n";
#ifdef TESTING
  cout <<"(IntUniformFV::Finalise) DONE.\n";
#endif
  cout <<"------------------------------------------------------------\n";
  return(0);
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
  


