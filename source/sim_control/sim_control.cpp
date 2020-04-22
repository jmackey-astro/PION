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
/// - 2018.05.** JM: moved most functions to new classes for calculating timestep,
///    updating boundaries, and time integration.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"
#include "dataIO/dataio_base.h"
#include "sim_control.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
using namespace std;


// ##################################################################
// ##################################################################


sim_control::sim_control()
{
}



// ##################################################################
// ##################################################################


sim_control::~sim_control()
{
#ifdef TESTING
  cout << "(sim_control::Destructor)\n";
#endif
}



// ##################################################################
// ##################################################################




/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int sim_control::Time_Int(
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      )
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control::Time_Int) STARTING TIME INTEGRATION."<<"\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  clk.start_timer("time_int"); double tsf=0;
  class MCMDcontrol ppar; // unused for serial code.
  err = update_evolving_RT_sources(SimPM,SimPM.simtime,grid[0]->RT);
  rep.errorTest("TIME_INT:: initial RT src update()",0,err);
  err = RT_all_sources(SimPM,grid[0],0);
  rep.errorTest("TIME_INT:: initial RT()",0,err);
  err+= output_data(grid);
  rep.errorTest("TIME_INT:: initial save",0,err);

  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid[0]);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid[0]);
#endif
    //
    // Update RT sources and do raytracing.
    //
    err = update_evolving_RT_sources(SimPM,SimPM.simtime,grid[0]->RT);
    rep.errorTest("TIME_INT::update_RT_sources()",0,err);
    err = RT_all_sources(SimPM,grid[0],0);
    rep.errorTest("TIME_INT:: loop RT()",0,err);

    //clk.start_timer("advance_time");
    // step forward by dt.
    SimPM.levels[0].last_dt = SimPM.last_dt;
    err += calculate_timestep(SimPM, grid[0],spatial_solver,0);
    rep.errorTest("TIME_INT::calc_timestep()",0,err);

    advance_time(0, grid[0]);
    //cout <<"advance_time took "<<clk.stop_timer("advance_time")<<" secs.\n";

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"New time: "<<SimPM.simtime;
    cout <<"\t dt="<<SimPM.dt;
    cout <<"\t steps: "<<SimPM.timestep;
    tsf=clk.time_so_far("time_int");
    cout <<"\t runtime: "<<tsf<<" s"<<"\n";
#endif
#endif
    err += check_energy_cons(grid[0]);

    err+= output_data(grid);
    if (err!=0) {
      cerr<<"(TIME_INT::output_data) err!=0 Something went wrong\n";
      return(1);
    }
    
    err+= check_eosim();
    if (err!=0) {
      cerr<<"(TIME_INT::) err!=0 Something went wrong\n";
      return(1);
    }
  }

  cout <<"(sim_control::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

  tsf=clk.time_so_far("time_int");
  cout <<"TOTALS ###: Nsteps: "<<SimPM.timestep<<" wall-time: ";
  cout <<tsf<<" time/step: "<<tsf/static_cast<double>(SimPM.timestep)<<"\n";
  cout <<"STEPS: "<<SimPM.timestep;
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
                spatial_solver->CellVolume(c,grid->DX());
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
    cout<<"grid does not contain shock.\n";
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
#ifdef TEST_CONSERVATION 
  // Energy, and Linear Momentum in x-direction.
  pion_flt u[SimPM.nvar];
  nowERG=0.;
  nowMMX = 0.;
  nowMMY = 0.;
  nowMMZ = 0.;
  nowMASS = 0.0;
  double totmom=0.0;

  class cell *cpt = grid->FirstPt();
  double dR=grid->DX();
  double dv = 0.0;
  do {
    dv = spatial_solver->CellVolume(cpt,dR);
    spatial_solver->PtoU(cpt->P,u,SimPM.gamma);
    nowERG += u[ERG]*dv;
    nowMMX += u[MMX]*dv;
    nowMMY += u[MMY]*dv;
    nowMMZ += u[MMZ]*dv;
    nowMASS += u[RHO]*dv;
    totmom += sqrt( u[MMX]*u[MMX]  + u[MMY]*u[MMY] + u[MMZ]*u[MMZ] )
              *dv;
  } while ( (cpt =grid->NextPt(cpt)) !=0);
  cout <<"(conserved quantities) ["<< nowERG <<", ";
  cout << nowMMX <<", ";
  cout << nowMMY <<", ";
  cout << nowMMZ <<", ";
  cout << nowMASS <<"]\n";
  cout <<"(relative error      ) ["<< (nowERG-initERG)/(initERG) <<", ";
  cout << (nowMMX-initMMX)/(totmom) <<", ";
  cout << (nowMMY-initMMY)/(totmom) <<", ";
  cout << (nowMMZ-initMMZ)/(totmom) <<", ";
  cout << (nowMASS-initMASS)/initMASS <<"]\n";
  
#endif // TEST_CONSERVATION
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
  


