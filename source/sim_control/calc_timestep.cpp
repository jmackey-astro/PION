/// \file calc_timestep.cpp
/// \brief routines for calculating the timestep on a grid.
/// \author Jonathan Mackey
/// \date 2018.05.10
/// 
/// Description:\n
/// Has a set of routines for calculating the timestep for fluid
/// dynamics simulations in PION.
/// 
/// Modifications:
/// - 2018.05.10 JM: moved from sim_control into its own class that
///    inherits from setup_fixed_grid.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "sim_control/calc_timestep.h"
#include "setup_fixed_grid.h"
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



calc_timestep::calc_timestep()
{
  return;
}



// ##################################################################
// ##################################################################



calc_timestep::~calc_timestep()
{
  return;
}



// ##################################################################
// ##################################################################




int calc_timestep::calculate_timestep(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      class FV_solver_base *sp_solver ///< solver/equations class
      )
{
#ifdef TESTING
  cout <<"calc_timestep::calc_timestep(): g="<<grid<<", rt="<<grid->RT<<"\n";
#endif
  //
  // This is a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt(par,grid,sp_solver);
  t_mp  = calc_microphysics_dt(par,grid);
#ifdef TESTING
  if (t_mp<t_dyn)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";
#endif
  par.dt = min(t_dyn,t_mp);

#ifdef THERMAL_CONDUCTION
  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in sp_solver->preprocess_data()
  //
  double t_cond = calc_conduction_dt_and_Edot();
#ifdef TESTING
  if (t_cond<t_dyn && t_cond<t_mp) {
    cout <<"CONDUCTION IS LIMITING TIMESTEP: t_c="<<t_cond<<", t_m="<<t_mp;
    cout <<", t_dyn="<<t_dyn<<"\n";
  }
#endif
  par.dt = min(par.dt, t_cond);
#endif // THERMAL CONDUCTION

  //
  // If using MHD with GLM divB cleaning, the following sets the hyperbolic wavespeed.
  // If not, it does nothing.  By setting it here and using t_dyn, we ensure that the
  // hyperbolic wavespeed is equal to the maximum signal speed on the grid, and not
  // an artificially larger speed associated with a shortened timestep.
  //
  sp_solver->GotTimestep(t_dyn);

  //
  // Check that the timestep doesn't increase too much between steps, and that it 
  // won't bring us past the next output time or the end of the simulation.
  // This function operates on par.dt, resetting it to a smaller value if needed.
  //
  timestep_checking_and_limiting(par);

  // sets the timestep info in the solver class.
  sp_solver->Setdt(par.dt);

#ifdef TESTING
  cout <<"calc_timestep::calc_timestep() finished.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



#ifdef THERMAL_CONDUCTION
double calc_timestep::calc_conduction_dt_and_Edot(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      class FV_solver_base *sp_solver ///< solver/equations class
      )
{
  //
  // N.B. THIS IS VERY AD-HOC CODE THAT HAS NEVER BEEN USED FOR
  // PRODUCTION SIMULATIONS.  USE AT YOUR OWN RISK!
  //
  // First we need to set Edot() in every cell.  This is stored in
  // c->dU[ERG] since this is where it will be updated.
  //
  //cout <<"\tCCdt: setting Edot.\n";
  sp_solver->set_thermal_conduction_Edot();
  //cout <<"\tCCdt: Done with div(Q).  Now getting timestep.\n";

  //
  // Now run through the grid and calculate E_int/Edot to get the smallest
  // timestep we can allow.  Also reset dU[RHO] to zero (was used as temporary
  // variable to hold gas temperature).
  // Note that we can allow a huge amount of energy into a cell, we just don't 
  // want more energy to leave than is already in it.  So we only count a cell
  // in the timestep limiting if Edot is negative.
  //
  double dt=1.0e200, gm1=par.gamma-1.0, tempdt=0.0, minP=par.RefVec[PG]*1.0e-3;
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
  if (par.timestep==0) dt=min(dt,1.0e3);

  //
  // Set timestep to be 1/10th of the thermal conduction timescale.
  //
  return 0.1*dt;
} // calc_conduction_dt_and_Edot()
#endif // THERMAL CONDUCTION


// ##################################################################
// ##################################################################


void calc_timestep::timestep_checking_and_limiting(
      class SimParams &par      ///< pointer to simulation parameters
      )
{
  //
  // If the timestep is less than the minimum allowed, report an error
  // and bug out!  This must be before we check for the next output
  // time because that can limit the step arbitrarily.
  //
  if (par.dt < par.min_timestep) {
    ostringstream temp; temp.str("");
    temp <<"Timestep too short! dt="<<par.dt<<"  min-step="<<par.min_timestep;
    rep.error(temp.str(),par.dt);
  }

  //
  // So that we don't increase the timestep by more than 30% over last step:
  //
#ifdef TIMESTEP_LIMITING
  if (par.dt > 1.3*par.last_dt)
    cout <<"limiting step from "<<par.dt<<" to "<<1.3*par.last_dt<<"\n";
  par.dt = min(par.dt,1.3*par.last_dt);
#endif // TIMESTEP_LIMITING

  //
  // If we are outputting every n-years, then check if we need to adjust dt.
  //
  if (par.op_criterion==1) {
    par.dt = min(par.dt, par.next_optime-par.simtime);
    if (par.dt <= 0.0)
      rep.error("Went past output time without outputting!",par.dt);
  }

  //
  // Make sure we end up exactly at finishtime:
  //
  par.dt = min(par.dt, par.finishtime-par.simtime);
  if (par.dt <= 0.0)
    rep.error("Went past Finish time without Stopping!",par.dt);

  return;
}


// ##################################################################
// ##################################################################


double calc_timestep::calc_dynamics_dt(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,
      class FV_solver_base *sp_solver ///< solver/equations class
      )
{
  double tempdt=0.0;
  double dt=1.e100; // Set it to very large no. initially.
  double dx=grid->DX();

  class cell *c=grid->FirstPt();
#ifdef TESTING
  dp.c = c;
#endif

  //
  // For the first 10 timesteps, we calculate the timestep for
  // boundary data too, because we can have some setups with rapid
  // inflow (such as a jet), which doesn't appear on the grid
  // immediately.
  //
#ifndef RT_TEST_PROBS
  if (par.timestep<=10) {
     c = grid->FirstPt_All();
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
    tempdt = sp_solver->CellTimeStep(
              c, ///< pointer to cell
              par.gamma, ///< gas EOS gamma.
              dx  ///< Cell size dx.
              );
    if(tempdt<=0.0) {
      CI.print_cell(c);
      cout <<"celltimestep="<<tempdt<<", gamma="<<par.gamma<<", dx="<<dx<<"\n";
      rep.error("CellTimeStep function returned failing value",c->id);
    }
    //    commandline.console("timestep -> ");
    dt = min(dt, tempdt);
    //cout <<"(get_min_timestep) i ="<<i<<"  min-dt="<<mindt<<"\n";

    //
    // Move to next cell.  For the first few steps, we include
    // boundary data in this, using NextPt_All(); otherwise use the
    // NextPt() command to only take grid cells.
    //
#ifndef RT_TEST_PROBS
    if (par.timestep<=10) {
      c = grid->NextPt_All(c);
    }
    else {
#endif // not RT_TEST_PROBS
      c = grid->NextPt(c);
#ifndef RT_TEST_PROBS
    }
#endif // not RT_TEST_PROBS

  } while (c != 0);

  if (dt <= 0.0)
    rep.error("Got zero timestep!!!",dt);
#ifdef TESTING
  cout <<"(calc_dynamics_dt)  min-dt="<<dt<<"\n";
#endif

  return dt;
}


// ##################################################################
// ##################################################################


double calc_timestep::calc_microphysics_dt(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid ///< pointer to grid.
      )
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
  if (par.EP.MP_timestep_limit==0) {
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
    //cout <<"calc_timestep, getting column densities rt="<<raytracer<<".\n";
    int err = calculate_raytracing_column_densities(par,grid->RT);
    if (err) rep.error("calc_MP_dt: bad return value from calc_rt_cols()",err);
    dt = get_mp_timescales_with_radiation(par,grid);
    if (dt<=0.0)
      rep.error("get_mp_timescales_with_radiation() returned error",dt);
  }
  else {
    //
    // don't need column densities, so call the no-RT version
    //
    //cout <<" getting timestep with no radiation\n";
    dt = get_mp_timescales_no_radiation(par,grid);
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


double calc_timestep::get_mp_timescales_no_radiation(
      class SimParams &par,      ///< pointer to simulation parameters
        class GridBaseClass *grid
        )
{
#ifdef TESTING
  //
  // paranoid checking...
  //
  if (par.EP.MP_timestep_limit==0) {
    cout <<"calc_timestep::get_mp_timescales_no_radiation() called, but no MP-dt limiting!\n";
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
      switch (par.EP.MP_timestep_limit) {
      case 1: // cooling time only.
        tempdt = MP->timescales(c->Ph, par.gamma, true, false, false);
        break;
      case 4: // recomb only
        tempdt = MP->timescales(c->Ph, par.gamma, false, true, false);
        break;
      case 2: // cooling+recomb
        tempdt = MP->timescales(c->Ph, par.gamma, true,  true, false);
        break;
      case 3: // cooling+recomb+ionisation (not working!)
        tempdt = MP->timescales(c->Ph, par.gamma, true,  true, true);
        break;
      default:
        rep.error("Bad MP_timestep_limit",par.EP.MP_timestep_limit);
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
  if ((par.timestep<3) && (par.RS.Nsources>0) &&
      (par.EP.phot_ionisation)) {
    //
    // adjust first timestep so that it corresponds to ionised cell.
    // 
    dt = min(dt, grid->DX()*par.CFL/2.0e6); // 20 km/s wavespeed.
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
  if ((par.timestep==0) && (par.RS.Nsources>0)) {
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


double calc_timestep::get_mp_timescales_with_radiation(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid
      )
{
#ifdef TESTING
  //
  // paranoid checking...
  //
  if (par.EP.MP_timestep_limit==0) {
    cout <<"calc_timestep::get_mp_timescales_with_radiation() no MP-dt limiting\n";
    return -1.0;
  }
  if (par.RS.Nsources==0)
    rep.error("calc_timestep::get_mp_timescales_with_radiation() no sources",1);
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
      tempdt = MP->timescales_RT(c->Ph, FVI_nheat, FVI_heating_srcs, FVI_nion, FVI_ionising_srcs, par.gamma);
#ifdef TESTING
      if (tempdt<dt) {
        cout <<"(get_min_timestep) id="<<c->id<<":  dt="<<tempdt<<", min-dt="<<dt;
        cout <<".\t 1-x="<<1.0-c->Ph[par.ftr]<<", pg="<<c->Ph[PG]<<"\n";
      }
#endif
      if (tempdt<=0.0) {
        cout <<"get_mp_timescales_with_radiation() negative timestep... ";
        cout <<"c->id="<<c->id<<"\tdt="<<tempdt<<"\n";
        rep.printVec("Ph",c->Ph,par.nvar);
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
// ##################################################################




int calc_timestep::calculate_raytracing_column_densities(
      class SimParams &par,      ///< pointer to simulation parameters
      //class GridBaseClass *grid       ///< Computational grid.
      class RayTracingBase *raytracer ///< raytracer for this grid.
      )
{
  int err=0;
  if (!raytracer) rep.error("calculate_raytracing_column_densities() no RT",0);
  //
  // If we have raytracing, we call the ray-tracing routines 
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
#ifdef raytracer_TESTING
    cout <<"calc_raytracing_col_dens: SRC-ID: "<<isrc<<"\n";
#endif
    err += raytracer->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      cout <<"isrc="<<isrc<<"\t"; 
      rep.error("calc_raytracing_col_dens step in returned error",err);
    } // if error
  } // loop over sources
  return err;
}



// ##################################################################
// ############   END OF TIMESTEP CALCULATION FUNCTIONS  ############
// ##################################################################



