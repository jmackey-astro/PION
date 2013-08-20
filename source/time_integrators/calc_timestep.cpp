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
///  - 2007-06-26 Hunted for bugs, but they were in another file.
///  - 2007-07-13 New Class structure implemented.
///  - 2007-07-16 New Class structure works (modified from Friday a little).  Same as old, but a bit faster.
///  - 2007-07-17 2D riemann problem initial setup improved (area averaged interface cells).
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08-07 cylindrical coordinates working for 2d axi-symmetry.
///  - 2007-10-11 Cleaning up.
///  - 2007-11-01 Added dataio class last week, cleaning up today.
///  - 2008-09-20 Removed text I/O into its own class. ifdeffed silo/fits.
///
///  - JM 2009-12-16 Added ifdef in IntUniformFV::Time_Int() so that I
///      can get the code to output magnetic pressure instead of
///      timing info every timestep.  This is purely to make a plot of
///      magnetic pressure for the Field Loop Advection Test and
///      should be switched off in general (although it makes no
///      difference to the running of the code!).
///
/// - 2009-12-18 JM: Added Axisymmetric Class
///    (cyl_FV_solver_Hydro_Euler) in set_equations().
///
/// - 2010-01-05 JM: Added check for microphysics timestep limiting in
///     calc_timestep().  This is controlled by the flag
///     SimPM.EP.MP_timestep_limit, which is set to true to turn it
///     on.
///
/// - 2010-01-21 JM: Added override option for SimPM.gamma, the
///    equation of state parameter.
///
/// - 2010-04-10 JM: fixed width timestep in text and fits filenames.
///
/// - 2010-04-21 JM: Changed filename setup so that i can write
/// checkpoint files with fname.999999.txt/silo/fits.  Added a check
/// for checkpointing in output_data() and removed all of the outfile
/// string generation (moved to dataio classes).
///
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
///  - 2010.07.23 JM: New RSP source position class interface.
///  - 2010.09.27 JM: took out comments and old ifdefs from dU_column().
///  - 2010.09.30 JM: Added div(V) calculation to calc_dU() function for viscosity.
///  - 2010.10.01 JM: Spherical coordinates(1D only) for Euler equations.
///       Cut out testing myalloc/myfree
///  - 2010.10.04 JM: Moved the field-loop magnetic pressure output to
///       an ifdeffed function.  Added a new function to calculate the
///       1D blast wave radius in spherical coordinates (also ifdeffed).
///  - 2010.10.05 JM: Changed boundary point timestep calculation to
///       work for all sims, not just jet sims.  Note this only works
///       for the XN point from the grid FirstPt(), not all boundary
///       points (which it really should work for!)
///  - 2010.10.13 JM: Added option to setup new microphysics_lowZ class
///       in setup_microphysics() function.
///       Added MP_timestep_limit override option.
///  - 2010.11.03 JM: Changed "endl" to "\n" for JUROPA.  Added a
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
/// - 2012.08.05 JM: Moved timestep-calculation functions from gridMethods.cc
///    to time_integrators/calc_timestep.cpp.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.08.20 JM: Changed cell_interface for radiative transfer
///    variables, so heating/ionising source syntax has changed.

#include "../defines/functionality_flags.h"
#include "../defines/testing_flags.h"

#ifdef NEW_TIME_UPDATE

#include "grid.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"
#include "spatial_solvers/solver_eqn_base.h"



#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <climits>
using namespace std;



// ##################################################################
// ##################################################################

/*****************************************************************/
/********************* TIMESTEP CALCULATION **********************/
/*****************************************************************/

// ##################################################################
// ##################################################################


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
    // SEEMS TO WORK BETTER
    if (fabs(c->Ph[PG]>minP)) tempdt = c->Ph[PG]/(gm1*(fabs(c->dU[ERG] +TINYVALUE)));
    else                tempdt = 1.0e200;
    dt = min(dt,tempdt);
  } while ((c=grid->NextPt(c)) !=0);
  //cout <<"Final conduction dt="<<dt<<", to be multiplied by 0.3.\n";

  //
  // first timestep can be dodgy with conduction and stellar winds,
  // so set it to a very small value
  //
  if (SimPM.timestep==0) dt=min(dt,1.0e3);

  //
  // Set timestep to be 1/10th of the thermal conduction timescale.
  //
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
#ifdef TESTING
  cout <<"(calc_dynamics_dt)  min-dt="<<dt<<"\n";
#endif

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
        CI.get_cell_col(c, FVI_heating_srcs[v].id, FVI_heating_srcs[v].DelCol);
        CI.get_col(     c, FVI_heating_srcs[v].id, FVI_heating_srcs[v].Column);
        for (short unsigned int iC=0; iC<FVI_heating_srcs[v].NTau; iC++)
          FVI_heating_srcs[v].Column[iC] -= FVI_heating_srcs[v].DelCol[iC];
      }
      for (int v=0; v<FVI_nion; v++) {
        FVI_ionising_srcs[v].Vshell = CI.get_cell_Vshell(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].dS     = CI.get_cell_deltaS(c, FVI_ionising_srcs[v].id);
        CI.get_cell_col(c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].DelCol);
        CI.get_col(     c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].Column);
        for (short unsigned int iC=0; iC<FVI_ionising_srcs[v].NTau; iC++)
          FVI_ionising_srcs[v].Column[iC] -= FVI_ionising_srcs[v].DelCol[iC];
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
      //if (tempdt<dt)
      //  cout <<"c->id="<<c->id<<"\tdt="<<tempdt<<"\n";
      dt = min(dt, tempdt);
    } // if not boundary data.
  } while ( (c =grid->NextPt(c)) !=0);

  return dt;
}


// ##################################################################
// ############   END OF TIMESTEP CALCULATION FUNCTIONS  ############
// ##################################################################

#endif // NEW_TIME_UPDATE

