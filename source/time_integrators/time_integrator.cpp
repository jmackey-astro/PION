/// \file gridMethods.cc
/// 
/// \brief Grid Methods Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for sim_control_fixedgrid 
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
///  - JM 2009-12-16 Added ifdef in sim_control_fixedgrid::Time_Int() so that I
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
/// - 2012.08.05 JM: Moved time-integration functions from
///    gridMethods.cc to time_integrators/time_integrator.cpp (with
///    old code ifdeffed out).  Started working on new integration
///    scheme which should be more accurate (and truly 2nd order).
/// - 2012.08.16 JM: Debugging.  It seems to be working well now, but
///    there is still more testing to do.
/// - 2013.08.19 JM: Some cosmetic changes only.
/// - 2013.08.20 JM: Modified cell_interface for optical depth vars.
/// - 2013.10.13 JM: Fixed bug in dU_Column relating to internal
///    boundaries; seems it never arose before.
/// - 2013.12.03 JM: Modified NO_COOLING_ON_AXIS hack.
/// - 2015.01.12/13 JM: Modified for new code structure; began adding
///    the grid pointer everywhere.
/// - 2015.01.26 JM: Renamed class to sim_control_fixedgrid.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"

#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

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

#include "spatial_solvers/solver_eqn_base.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <climits>
using namespace std;


#define TIMESTEP_FULL OA2
#define TIMESTEP_FIRST_PART OA1



// ##################################################################
// #################    TIME UPDATE FUNCTIONS     ###################
// ##################################################################

// ##################################################################
// ##################################################################

int sim_control_fixedgrid::advance_time(
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;

  //
  // First we calculate the timestep.  The microphysics timescales may depend
  // on the ray-tracing column densities, and if so all the column densities
  // will be calculated with raytracing calls in calc_mp_timestep()
  //
  err += calc_timestep(grid);
  if (err) 
    rep.error("advance_time: bad return value from calc_timestep()",err);

  //
  // Check order-of-accuracy (OOA) requested, and perform the
  // appropriate update.
  //
  if      (SimPM.tmOOA==OA1 && SimPM.spOOA==OA1 ) {
    //
    // Send in full timestep, and the order-of-accuracy, so the 
    // function knows whether to update P[i] as well as Ph[i]
    //
    //cout <<"First order update\n";
    err += first_order_update(SimPM.dt, SimPM.tmOOA, grid);
    if (err)
      rep.error("first_order_update() returned error",err);
  }

  else if (SimPM.tmOOA==OA2 && SimPM.spOOA==OA2) {
    //
    // C2ray-type time-update requires full splitting of microphysics
    // and dynamics because it is implicit.  This is legacy code.
    //
    if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT) {
#ifdef RT_TESTING
      cout <<"--- tstep_dyn_then_mp() selected for implicit RT.\n";
#endif
      err += timestep_dynamics_then_microphysics(grid);
      if (err)
        rep.error("tstep_dyn_then_mp() returned error",err);
    }
    //
    // Normal update is not split; first a full 1st order step is
    // performed, then a full second-order step.
    //
    else {
      //cout <<"Second order update\n";
      err += first_order_update( 0.5*SimPM.dt, SimPM.tmOOA, grid);
      err += second_order_update(SimPM.dt,     SimPM.tmOOA, grid);
      if (err)
        rep.error("Second order time-update returned error",err);
    }
  }
  //
  // Add in 3rd order PPM at some stage???
  //
  else {
    rep.error("Bad OOA requests; choose (1,1) or (2,2)",SimPM.tmOOA);
  }

  //
  // Update timestepping variables to new state.
  //
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime +=SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;

  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::first_order_update(
      const double dt,
      const int   ooa,
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // Set dt for equations class
  // MULTITHREADING RISK!
  //
  eqn->Setdt(dt);

  //
  // May need to do raytracing, if it wasn't needed for calculating
  // the timestep.
  //
  if (!FVI_need_column_densities_4dt) {
    err += calculate_raytracing_column_densities();
    if (err) 
      rep.error("first_order_update: error from first calc_rt_cols()",err);
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(dt, grid);
  err += calc_dynamics_dU(dt,OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt,OA1, grid);
#endif // THERMAL_CONDUCTION
  if (err) 
    rep.error("first_order_update: error from calc_*_dU",err);
  
  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(dt,TIMESTEP_FIRST_PART,ooa, grid);
  if (err) 
    rep.error("first_order_update: error from state-vec update",err);

  //
  // Update boundary data.
  //
  err += grid->TimeUpdateInternalBCs(OA1,ooa);
  err += grid->TimeUpdateExternalBCs(OA1,ooa);
  if (err) 
    rep.error("first_order_update: error from bounday update",err);

  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::second_order_update(
      const double dt,
      const int   ooa,
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // Set dt for equations class
  // MULTITHREADING RISK!
  //
  eqn->Setdt(dt);

  //
  // Raytracing, to get column densities for microphysics update.
  //
  err += calculate_raytracing_column_densities();
  if (err) {
    rep.error("second_order_update: error from first calc_rt_cols()",err);
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(      dt,      grid);
  err += calc_dynamics_dU(          dt, OA2, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt, OA2, grid);
#endif // THERMAL_CONDUCTION
  if (err) 
    rep.error("second_order_update: error from calc_*_dU",err);
  
  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(  dt, TIMESTEP_FULL, ooa, grid);
  if (err) 
    rep.error("second_order_update: error from state-vec update",err);

  //
  // Update boundary data.
  //
  err += grid->TimeUpdateInternalBCs(   OA2, ooa);
  err += grid->TimeUpdateExternalBCs(   OA2, ooa);
  if (err) 
    rep.error("second_order_update: error from bounday update",err);

  return 0;
}




// ##################################################################
// ##################################################################





// ##################################################################
// ##################################################################



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::timestep_dynamics_then_microphysics(
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef TESTING
  cout <<"Using  sim_control_fixedgrid::timestep_dynamics_then_microphysics() update.\n";
#endif // TESTING
  
  // ----------------------------------------------------------------
  // THIS IS PROBABLY ONLY 1ST ORDER IN TIME.  USE A DIFFERENT   ----
  // INTEGRATOR!                                                 ----
  // ----------------------------------------------------------------

  // This is a general conservative scheme, with all the details in
  // the calculation of the intercell fluxes.

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

  int err=0;

  //int sp = OA1;
  //int tm = OA1;
  //if (SimPM.tmOOA==OA1 && SimPM.spOOA==OA2) sp = OA2;
  //  cout <<"dt = "<<SimPM.dt<<"\n";
  
  double dt = SimPM.dt;

  if (SimPM.tmOOA ==OA1) { // First order time time
    eqn->Setdt(dt);
    err  = calc_dynamics_dU(dt,OA1, grid);
    //     cout <<"updating microphysics.\n";
    err += calc_microphysics_dU(dt, grid);
    //    cout <<"done with mp.\n";

    err += grid_update_state_vector(dt,TIMESTEP_FIRST_PART,OA1, grid);
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    //     cout <<"updating external bcs.\n";
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    //    cout <<"done with external bcs.\n";
    if (err) rep.error("O1 time update loop generated errors",err);
  } // if 1st order accurate
  
  else if (SimPM.tmOOA==OA2) {
    //cout <<"second order!!!\n";
    //
    // The Falle, Komissarov, Joarder (1998) 2nd order-accurate
    // (time and space) time update.  BUT THERE ARE MISUNDERSTANDINGS
    // AND THIS IS PROBABLY ONLY 1ST ORDER IN TIME.  USE A DIFFERENT
    // INTEGRATOR!
    //

    //
    // Halve the timestep for the first pass
    //
    dt /= 2.0;
    eqn->Setdt(dt);
    err  = calc_dynamics_dU(dt,OA1, grid);

    err += grid_update_state_vector(dt,TIMESTEP_FIRST_PART,OA2, grid);
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(OA1,SimPM.tmOOA);
    if (err) rep.error("O2 half time update loop generated errors",err);
    
    //
    // Now calculate dU again, for the full timestep.
    //
    dt *= 2.0;
    eqn->Setdt(dt);
    //    cout <<"\tsecond pass.... \n";
    err = calc_dynamics_dU(dt,OA2, grid); //,SimPM.tmOOA);
    // Update MicroPhysics, if present
    err += calc_microphysics_dU(dt, grid);

    err += grid_update_state_vector(dt,TIMESTEP_FULL,OA2, grid);
    err += grid->TimeUpdateInternalBCs(SimPM.tmOOA,SimPM.tmOOA);
    err += grid->TimeUpdateExternalBCs(SimPM.tmOOA,SimPM.tmOOA);

    if (err) rep.error("O2 full time update loop generated errors",err);
  } // If 2nd order accurate.
  
  else rep.error("Only know first and second order accuracy",SimPM.tmOOA);
  
#ifdef TESTING
  if (SimPM.timestep%20 ==0) {
    //    cout <<"dp.initERG = "<<dp.initERG<<"\n";
    check_energy_cons(grid);
  }
#endif // TESTING
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime += SimPM.dt;
  SimPM.last_dt  = SimPM.dt;
  SimPM.timestep++;
  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::calculate_raytracing_column_densities(
      //class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // If we have raytracing, we call the new ray-tracing function 
  // to get Tau0, dTau, Vshell in cell->extra_data[n].
  //
  if (RT) {
    for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
#ifdef RT_TESTING
      cout <<"sim_control_fixedgrid::calc_RT_col_dens: SRC-ID: "<<isrc<<"\n";
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



int sim_control_fixedgrid::calc_microphysics_dU(
      const double delt, ///< timestep to integrate MP eqns.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //cout <<"\tcalc_microphysics_dU starting.\n";

  //
  // If we are not doing microphysics or raytracing, the modules will
  // not be initialised, so just return zero.  (RT requires MP)
  //
  if (!MP) return 0;
  
#ifdef TESTING
  cout <<"calc_microphysics_dU() Updating MicroPhysics. ";
  cout <<"  RT-Nsrc="<<SimPM.RS.Nsources<<"\n";
#endif // TESTING
  int err = 0;

  if (!RT) {
    //
    // If no radiative transfer, then just do a simple MP update.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling calc_microphysics_dU_no_RT()\n";
#endif // RT_TESTING
    err += calc_microphysics_dU_no_RT(delt, grid);
  }

  else if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT) {
    //
    // This is the C2-ray single-source update which I did for my thesis.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling calc_microphysics_dU_JMs_C2ray_RT()\n";
#endif // RT_TESTING
    err += calc_microphysics_dU_JMs_C2ray_RT(delt, grid);
  }

  else {
    //
    // must have at least one radiation source, and the RT and microphysics
    // are run separately, and we may have diffuse radiation, so we do the
    // new RT update:
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling calc_microphysics_dU_general_RT()\n";
#endif // RT_TESTING
    err += calc_microphysics_dU_general_RT(delt, grid);
  }
    
  //cout <<"\tcalc_microphysics_dU finished.\n";
  return err;
}


// ##################################################################
// ##################################################################



int sim_control_fixedgrid::calc_microphysics_dU_general_RT(
      const double delt, // timestep to integrate
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef RT_TESTING
  if (!RT)
    rep.error("Logic error: must have RT unless i'm an idiot",
              "GENERAL-RT");
  if (SimPM.RS.Nsources<1)
    rep.error("Need at least one source for diffuse-RT update",
              SimPM.RS.Nsources);
#endif // RT_TESTING

  int err=0;
  //
  // Update microphyiscs with new interface.
  // RT source properties are already in structs for the microphysics calls.
  //

  //
  // Now do MP update.  Only new MP classes have the  XX_RTnew() defined,
  // so if we try to call this with old code then it should return with 
  // a non-zero error code.
  //
  cell *c = grid->FirstPt();
  double p[SimPM.nvar]; // temporary state vector for output state.
  double ui[SimPM.nvar], uf[SimPM.nvar]; // conserved variable states.

  double tt=0.; // temperature returned at end of microphysics step.
  do {
    //
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in calc_microphysics_dU_general_RT() c->isbd.\n";
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
#ifdef TESTING
        //cout <<"HEAT: Vs="<<FVI_heating_srcs[v].Vshell<<", dS="<<FVI_heating_srcs[v].dS<<", dC="<<FVI_heating_srcs[v].DelCol<<", Col="<<FVI_heating_srcs[v].Column<<"\n";
#endif
      }
      for (int v=0; v<FVI_nion; v++) {
        FVI_ionising_srcs[v].Vshell = CI.get_cell_Vshell(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].dS     = CI.get_cell_deltaS(c, FVI_ionising_srcs[v].id);
        CI.get_cell_col(c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].DelCol);
        CI.get_col(     c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].Column);
        for (short unsigned int iC=0; iC<FVI_ionising_srcs[v].NTau; iC++)
          FVI_ionising_srcs[v].Column[iC] -= FVI_ionising_srcs[v].DelCol[iC];
      }
      //if (c->id<=3) {
      //  cout <<"*-*-*-*-*- i="<<c->id<<", --Tau="<<FVI_ionising_srcs[0].Column[0];
      //  cout <<", dTau="<<FVI_ionising_srcs[0].DelCol[0];
      //  cout <<", NTau="<<FVI_ionising_srcs[0].NTau<<"\n";
      //}
      //
      // integer 9th argument determines type of integration substepping:
      // 0 = adaptive RK5 Cash-Karp method (use this!).
      // 4th and 5th args are for ionising sources.
      //
      err += MP->TimeUpdateMP_RTnew(c->P, FVI_nheat, FVI_heating_srcs, FVI_nion, FVI_ionising_srcs,
                                    p, delt, SimPM.gamma, 0, &tt);

//#define NO_COOLING_ON_AXIS
#ifdef NO_COOLING_ON_AXIS
      //cout <<"hello\n";
//#error "Fix HACK in time_integrator.cpp"
      if (SimPM.coord_sys==COORD_CYL && 
          //!grid->NextPt(c,YN)->isgd &&
          p[RO] > 1.4e-20 &&    // density more than twice background density
          //p[RO] > 0.77e-20 &&    // density more than 1.1x background density
          c->pos[Rcyl] < 5 &&   // get the first three radial cells (R=0,2,4)
          c->pos[Zcyl] > 0 &&   // only consider cells with z>0 (upstream)
          p[PG] < c->P[PG] &&   // only consider cells that were cooled (not heated)
          p[SimPM.ftr] < 0.5    // only consider mostly neutral gas
          ) {
        //tt = MP->Temperature(p,SimPM.gamma);
        //if (tt < 1.0e3) {	
        //  MP->Set_Temp(p,1.0e3,SimPM.gamma);
        //}
        // Just set the pressure equal to what it was before cooling,
        // so that the gas is adiabatic.
        p[PG] = c->P[PG];
      }
#endif // NO_COOLING_ON_AXIS

      //
      // New state is p[], old state is c->P[].  Get dU from these.
      //
      eqn->PtoU(c->P,ui,SimPM.gamma);
      eqn->PtoU(p,   uf,SimPM.gamma);
      for (int v=0;v<SimPM.nvar;v++) c->dU[v] += uf[v]-ui[v];

    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"calc_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
} // newest general-RT microphysics update.



// ##################################################################
// ##################################################################



int sim_control_fixedgrid::calc_microphysics_dU_JMs_C2ray_RT(
      const double delt, ///< timestep to integrate
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef RT_TESTING
  if (!RT) rep.error("Logic error: must have RT unless i'm an idiot","C2RAY");
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



int sim_control_fixedgrid::calc_microphysics_dU_no_RT(
      const double delt, ///< timestep to integrate
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef TESTING
  cout <<"calc_microphysics_dU_no_RT starting.\n";
#endif
  //
  // No radiation sources at all, and no diffuse radiation optical depths,
  // so just do my original microphysics update.
  //
  cell *c = grid->FirstPt();
  double p[SimPM.nvar]; // temporary state vector for output state.
  double ui[SimPM.nvar], uf[SimPM.nvar]; // conserved variable states.
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
      cout <<"skipping cell "<<c->id<<" in calc_microphysics_dU_no_RT() c->isbd.\n";
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
      if (err) rep.error("calc_microphysics_dU_no_RT returned error: cell id",c->id);

      //
      // New state is p[], old state is c->P[].  Get dU from these.
      //
      eqn->PtoU(c->P,ui,SimPM.gamma);
      eqn->PtoU(p,   uf,SimPM.gamma);
      for (int v=0;v<SimPM.nvar;v++) c->dU[v] += uf[v]-ui[v];

    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"calc_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
} 



// ##################################################################
// ##################################################################


  
int sim_control_fixedgrid::calc_dynamics_dU(
      const double dt, ///< timestep to integrate
      const int space_ooa, ///< spatial order of accuracy for update.
      //const int time_ooa,   ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //cout <<"\tcalc_dynamics_dU starting.\n";
  //
  // first check if we are doing dynamics, and return if not.
  //
  if (!SimPM.EP.dynamics) return 0;
  int err=0;

#ifdef TESTING
  if (space_ooa==OA1)
    cout <<"*****Updating dynamics: OA1\n";
  else if (space_ooa==OA2)
    cout <<"*****Updating dynamics: OA2\n";
  else rep.error("Bad ooa",space_ooa);
#endif //TESTING

  //
  // First we pre-process the cells, if needed.  This is required for
  // genuinely multi-dimensional viscosity such as Lapidus-like AV or
  // the H-Correction.
  //
  err = eqn->preprocess_data(dt, space_ooa, grid); //,time_ooa);

  //
  // Now calculate the directionally-unsplit time update for the
  // conserved variables:
  //
  err = set_dynamics_dU(dt, space_ooa, grid); //,time_ooa);
  rep.errorTest("calc_dynamics_dU() eqn->set_dynamics_dU returned error.",
                0,err);

  //
  // Post-processing is for if we are doing something like Constrained
  // Transport, where we have to change the B-field update.  At the
  // moment there is *NO* solver which does anything here since I
  // found the Dedner et al. (2002) divergence cleaning to be more
  // robust than e.g. Toth (2000) Field-CT method.  (well the internal
  // energy solver uses it, but it's not really worth using).
  //
  err = eqn->PostProcess_dU(dt, space_ooa, grid); //,time_ooa);
  rep.errorTest("calc_dynamics_dU() eqn->PostProcess_dU()",0,err);

  return 0;
}


// ##################################################################
// ##################################################################



int sim_control_fixedgrid::set_dynamics_dU(
      const double dt,     ///< timestep for this calculation
      const int space_ooa, ///< space OOA for this calculation
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //  cout <<"\t\t\tStarting set_dynamics_dU: ndim = "<<SimPM.ndim<<"\n";
  int return_value=0;

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
  // boundary cells, and also on dynamics_dU_column returning 0 on successful
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
    
#ifdef TESTING
    rep.printVec("cpt",cpt->pos,SimPM.ndim);
    rep.printVec("+XX",(grid->NextPt(cpt,XP))->pos,SimPM.ndim);
    if (SimPM.ndim>1) rep.printVec("+YY",(grid->NextPt(cpt,YP))->pos,SimPM.ndim);
    if (SimPM.ndim>2) rep.printVec("+ZZ",(grid->NextPt(cpt,ZP))->pos,SimPM.ndim);
#endif
    
    while (
      (return_value = dynamics_dU_column(cpt,posdirs[i],negdirs[i], dt,
#ifdef TESTING
      // this is a hack, assuming spatial o-o-a is the same as the
      // time o-o-a.  But this is only needed for checking energy
      // and momentum conservation (to know if we are on the full
      // or half step) so it is not too important.
                                        space_ooa,
#endif
                                        space_ooa, grid)) ==0) {
      //cout <<"next dir= "<<(i+1)%SimPM.ndim<<"\n";
      //rep.printVec("cpt",cpt->pos,SimPM.ndim);
      if ( !(cpt=grid->NextPt(cpt,posdirs[(i+1)%SimPM.ndim]))->isgd ) {
        //cout <<"next dir= "<<(i+2)%SimPM.ndim<<"\n";
        if ( !(cpt=grid->NextPt(marker,posdirs[(i+2)%SimPM.ndim]))->isgd ) {
          CI.print_cell(cpt);
          rep.error("set_dynamics_dU: Got to edge of box before last point!",cpt);
        }
        marker = cpt;
      } // if null pointer.
    } // Loop over columns.
    if (return_value!=-1) rep.error("dUdtColumn returned abnormally.",return_value);
  } // Loop over three directions.
  eqn->SetDirection(axis[0]); // Reset fluxes to x-dir, (just to be safe!).
  //  cout <<"\t\t\tCALC_DU done.\n";

  //
  // all done, so return
  //
  return 0;
}   // set_dynamics_dU()



// ##################################################################
// ##################################################################

  

int sim_control_fixedgrid::dynamics_dU_column
      (
      const class cell *startingPt, ///< sterting point of column.
      const enum direction posdir, ///< direction to trace column.
      const enum direction negdir, ///< reverse direction
      const double dt, ///< timestep we are advancing by.
#ifdef TESTING
      const int ctm, ///< time order-of-accuracy (for conservation)
#endif
      const int csp,  ///< spatial order-of-accuracy for this step.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //  cout <<"Starting dU_column in direction "<<posdir<<" at ";
  //  rep.printVec("starting position",startingPt->x,SimPM.ndim);
  if ( (SimPM.spOOA>2) || (SimPM.tmOOA>2) || (csp>2)  ) {
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
    err += eqn->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
    err += eqn->SetSlope(npt, axis, SimPM.nvar, slope_npt, csp, grid);
    err += eqn->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
    //    rep.printVec("El",edgeL,SimPM.nvar); rep.printVec("Er",edgeR,SimPM.nvar);
    //    rep.errorTest("Edge States not obtained!",0,err);
    err += eqn->InterCellFlux(grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType, SimPM.artviscosity, SimPM.gamma, SimPM.dx);
    //    rep.printVec("Fr",Fr_this,SimPM.nvar);    rep.errorTest("Intercell Flux not obtained!",0,err);
    err += eqn->dU_Cell(grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, SimPM.dx, dt);
    //    rep.errorTest("dU not obtained!",0,err);
#ifdef TESTING
    for (int v=0;v<SimPM.nvar;v++) {
      if(!isfinite(cpt->dU[v])) {
  rep.printVec("Fl",Fr_prev,SimPM.nvar); rep.printVec("Fr",Fr_this,SimPM.nvar);
  cout <<"dt:"<<dt<<"\tdx="<<SimPM.dx<<"\n";
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
      dp.initERG += Fr_this[ERG]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMX += Fr_this[MMX]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMY += Fr_this[MMY]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMZ += Fr_this[MMZ]*dt*eqn->CellInterface(cpt,posdir);
//      if (posdir==YP && fabs(Fr_this[MMY])>2*MACHINEACCURACY) {
//  cout <<"R-momentum flux entering domain from R=0(!) = "<<Fr_this[MMY]<<"\n";
//  cout <<"v_R in first cell = "<<npt->Ph[VY]<<", "<<npt->P[VY]<<"\n";
//      }
    }
    else if (ctm==SimPM.tmOOA && !(npt->isgd) && cpt->isgd) {
      ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
      //      cout <<"Leaving Domain Dir = "<<posdir<<" and interface area = "<<eqn->CellInterface(cpt,posdir)<<"\n";
      dp.initERG -= Fr_this[ERG]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMX -= Fr_this[MMX]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMY -= Fr_this[MMY]*dt*eqn->CellInterface(cpt,posdir);
      dp.initMMZ -= Fr_this[MMZ]*dt*eqn->CellInterface(cpt,posdir);
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
  err += eqn->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
  for (int v=0;v<SimPM.nvar;v++) slope_npt[v] = 0.; // last cell must be 1st order.
  err += eqn->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
  err += eqn->InterCellFlux(grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType, SimPM.artviscosity, SimPM.gamma, SimPM.dx);
  err += eqn->dU_Cell(grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, SimPM.dx, dt);
#ifdef TESTING
  if (ctm==SimPM.tmOOA && cpt->isgd && !(npt->isgd)) {
    ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
    // cout <<"Leaving Domain Dir = "<<posdir<<" and interface area = "<<eqn->CellInterface(cpt,posdir)<<"\n";
    dp.initERG -= Fr_this[ERG]*dt*eqn->CellInterface(cpt,posdir);
    dp.initMMX -= Fr_this[MMX]*dt*eqn->CellInterface(cpt,posdir);
    dp.initMMY -= Fr_this[MMY]*dt*eqn->CellInterface(cpt,posdir);
    dp.initMMZ -= Fr_this[MMZ]*dt*eqn->CellInterface(cpt,posdir);
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
  do{} while( !(cpt=grid->NextPt(cpt,negdir))->isgd );
  if (cpt->id == grid->LastPt()->id) return(-1);
  else return(0);
}



// ##################################################################
// ##################################################################

  
int sim_control_fixedgrid::grid_update_state_vector(
      const double dt,  ///< timestep
      const int step, ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa,   ///< Full order of accuracy of simulation
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // temp variable to handle change of energy when correcting for negative pressure.
  //
  double temperg =0.0;

  //
  // Loop through grid, updating Ph[] with CellAdvanceTime function.
  //
  class cell* c = grid->FirstPt();
  do {

#ifdef TESTING
    dp.ergTotChange = 0.;temperg =0.;
    dp.c = c;
#endif
    err += eqn->CellAdvanceTime(c, c->P, c->dU, c->Ph, &temperg,
                                SimPM.gamma, dt);
#ifdef TESTING
    if (err) {
      cout <<"______ Error in Cell-advance-time: "; CI.print_cell(c);
      err=0;
    }
#endif // TESTING
    
    //
    // If the current step is the full update, then also set
    // the updated base-state-vector to updated value.
    //
    if (step==ooa) {
      for (int v=0;v<SimPM.nvar;v++) c->P[v] = c->Ph[v];
#ifdef TESTING
      //
      // Update Total Energy from fixing negative pressures. Reset
      // update variables.
      //
      dp.ergTotChange = temperg;
      dp.initERG += dp.ergTotChange*eqn->CellVolume(c);
#endif // TESTING
    }

  } while ( (c =grid->NextPt(c)) !=0);

#ifdef TESTING
  cout <<"\tcalc_dynamics_dU done. error="<<err<<"\n";
#endif // TESTING
  return err; 
}


// ##################################################################
// ##################################################################








