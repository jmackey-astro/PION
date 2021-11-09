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
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "grid/setup_fixed_grid.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"
#include "sim_control/calc_timestep.h"
#include "spatial_solvers/solver_eqn_base.h"

#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <time.h>

#ifdef PION_OMP
#include <omp.h>
#endif

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
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level to advance (for NG grid)
)
{
#ifndef NDEBUG
  cout << "calc_timestep::calc_timestep(): g=" << grid << ", rt=" << grid->RT
       << "\n";
#endif
  //
  // This is a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn = 0.0, t_mp = 0.0;
  t_dyn = calc_dynamics_dt(par, grid);
  t_mp  = calc_microphysics_dt(par, grid, l);
  // cout <<"l="<<l<<", \t t_dyn="<<t_dyn<<"and t_mp ="<<t_mp<<"\n";

#ifndef NDEBUG
  if (t_mp < t_dyn)
    cout << "Limiting timestep by MP: mp_t=" << t_mp << "\thydro_t=" << t_dyn
         << "\n";
  cout << "t_dyn = " << t_dyn << "  and t_mp = " << t_mp << "\n";
#endif

  par.dt = min(t_dyn, t_mp);

  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in spatial_solver->preprocess_data()
  //
  double t_cond = set_conduction_dt_and_Edot(par, grid);
#ifndef NDEBUG
  if (t_cond < t_dyn && t_cond < t_mp) {
    cout << "CONDUCTION IS LIMITING TIMESTEP: t_c=" << t_cond
         << ", t_m=" << t_mp;
    cout << ", t_dyn=" << t_dyn << "\n";
  }
#endif
  par.dt = min(par.dt, t_cond);

  //
  // If using MHD with GLM divB cleaning, the following sets the
  // hyperbolic wavespeed.  If not, it does nothing.  By setting it
  // here and using t_dyn, we ensure that the hyperbolic wavespeed is
  // equal to the maximum signal speed on the grid, and not an
  // artificially larger speed associated with a shortened timestep.
  //
  // we always calculate timestep on finest level first.
  static double td = 0.0;
  if (l == par.grid_nlevels - 1)
    td = t_dyn;
  else
    td = min(td, t_dyn / pow(2.0, par.grid_nlevels - 1 - l));

  double cr = 0.0;
#ifdef PION_OMP
  #pragma omp parallel private(cr)
  {
#endif
    if (par.grid_nlevels == 1) {
      cr = 0.25 / par.dx;
      spatial_solver->Set_GLM_Speeds(td, par.dx, cr);
    }
    else {
      cr = 0.25 / par.levels[par.grid_nlevels - 1].dx;
      if (l == 0) {
        spatial_solver->Set_GLM_Speeds(
            td, par.levels[par.grid_nlevels - 1].dx, cr);
      }
    }
#ifdef PION_OMP
  }
#endif

  //
  // Check that the timestep doesn't increase too much between steps, and that
  // it won't bring us past the next output time or the end of the simulation.
  // This function operates on par.dt, resetting it to a smaller value if
  // needed.
  //
  timestep_checking_and_limiting(par, l);

#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    // sets the timestep info in the solver class.
    spatial_solver->Setdt(par.dt);
#ifdef PION_OMP
  }
#endif

#ifndef NDEBUG
  cout << "calc_timestep::calc_timestep() finished.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



double calc_timestep::set_conduction_dt_and_Edot(
    class SimParams &par,      ///< pointer to simulation parameters
    class GridBaseClass *grid  ///< pointer to grid.
)
{
  if (par.EP.sat_thermal_cond == 0) return 1.0e200;
  //
  // N.B. THIS IS VERY AD-HOC CODE THAT HAS NEVER BEEN USED FOR
  // PRODUCTION SIMULATIONS.  USE AT YOUR OWN RISK!
  //
  // First we need to set Edot() in every cell.  This is stored in
  // c->dU[ERG] since this is where it will be updated.
  //
  // cout <<"\tCCdt: setting Edot.\n";
  set_thermal_conduction_Edot(par, OA1, grid);
  // cout <<"\tCCdt: Done with div(Q).  Now getting timestep.\n";

  //
  // Now run through the grid and calculate E_int/Edot to get the smallest
  // timestep we can allow.  Also reset dU[RHO] to zero (was used as temporary
  // variable to hold gas temperature).
  // Note that we can allow a huge amount of energy into a cell, we just don't
  // want more energy to leave than is already in it.  So we only count a cell
  // in the timestep limiting if Edot is negative.
  //
  double dt = 1.0e200, gm1 = par.gamma - 1.0, tempdt = 0.0,
         minP = par.RefVec[PG] * 1.0e-3;
  cell *c     = grid->FirstPt_All();
  do {
    // DIDN'T WORK -- CHECKERBOARDING!
    // if (c->dU[ERG]<0.0) tempdt = c->Ph[PG]/(gm1*(fabs(c->dU[ERG]
    // +TINYVALUE))); else                tempdt = 1.0e200;
    // SEEMS TO WORK BETTER
    // if (fabs(c->Ph[PG] > minP))
    //  tempdt = c->Ph[PG] / (gm1 * (fabs(c->dU[ERG] + TINYVALUE)));
    // else
    //  tempdt = 1.0e200;
    // This is appropriate for saturated TC, where dU[VX] is set to
    // the max speed * dx, to cancel one of the dx factors in the
    // numerator (saturated TC is hyperbolic).
    tempdt = par.CFL * par.dx * par.dx * 0.5 / c->dU[VX];
    // cout <<par.dx<<"  "<<par.CFL<<"  "<<c->dU[VX]<<"\n";
    c->dU[VX] = 0.0;
    dt        = min(dt, tempdt);
  } while ((c = grid->NextPt_All(c)) != 0);
  // cout <<"Final conduction dt="<<dt<<", to be multiplied by 0.3.\n";

  //
  // first timestep can be dodgy with conduction and stellar winds,
  // so set it to a very small value
  //
  if (par.timestep == 0) dt = min(dt, 1.0e0);

  //
  // Set timestep to be 1/10th of the thermal conduction timescale.
  //
  return dt;
}



// ##################################################################
// ##################################################################



int calc_timestep::set_thermal_conduction_Edot(
    class SimParams &par,      ///< pointer to simulation parameters
    const int oa,              ///< order of accuracy required.
    class GridBaseClass *grid  ///< pointer to grid.
)
{
  if (par.EP.sat_thermal_cond == 0) return 0;
  //
  // This function is quite similar to calc_dU() and dU_column()
  // because it steps through the grid in the same way.
  //
  // cout <<"calc_conduction_dt_and_Edot()\n\tcalculating
  // Temperature.\n";
  //
  // First we need to calculate the temperature in every cell.  This
  // is stored in dU[RHO] -- reset to zero at the end of function.
  //
  if (!MP) rep.error("No Microphysics == No Conduction", MP);
  cell *c = grid->FirstPt_All();
  do {
    // cout <<"dU[RHO]="<<c->dU[RHO];
    c->dU[RHO] = MP->Temperature(c->Ph, par.gamma);
    // cout <<", replaced with T="<<c->dU[RHO]<<"\n";
  } while ((c = grid->NextPt_All(c)) != 0);

  // cout <<"\tT calculated, now calculating divQ.\n";
  //
  // Allocate arrays
  //
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP;
  posdirs[1] = YP;
  posdirs[2] = ZP;
  negdirs[0] = XN;
  negdirs[1] = YN;
  negdirs[2] = ZN;
  axis[0]    = XX;
  axis[1]    = YY;
  axis[2]    = ZZ;

  // Loop through each dimension.
  for (int idim = 0; idim < par.ndim; idim++) {
#ifndef NDEBUG
    cout << "\t\t\tidim=" << idim << "\n";
#endif  // NDEBUG

#ifdef PION_OMP
    #pragma omp parallel
    {
      spatial_solver->SetDirection(axis[idim]);
    }
#else
    spatial_solver->SetDirection(axis[idim]);
#endif  // PION_OMP
    class cell *cpt    = grid->FirstPt_All();
    class cell *marker = cpt;
#ifdef TEST_INT
    cout << "Direction=" << axis[idim] << ", i=" << idim << "\n";
    rep.printVec("cpt", cpt->pos, par.ndim);
#endif
    // loop over the number of cells in the line/plane of starting
    // cells.
    enum axes x1 = axis[(idim + 1) % 3];
    enum axes x2 = axis[(idim + 2) % 3];
    enum axes x3 = axis[idim];
    // loop over the two perpendicular axes, to trace out a plane of
    // starting cells for calculating fluxes along columns along this
    // axis.  The plane can be a single cell
    // (in 1D) or a line (in 2D) or a plane (in 3D).
    int index[3];
#ifdef PION_OMP
    #pragma omp parallel
    {
      #pragma omp for collapse(2) private(index,cpt)
#endif
      for (int ax2 = 0; ax2 < grid->NG_All(x2); ax2++) {
        for (int ax1 = 0; ax1 < grid->NG_All(x1); ax1++) {
          index[x1]    = static_cast<int>(ax1);
          index[x2]    = static_cast<int>(ax2);
          index[x3]    = 0;
          cpt          = grid->get_cell_all(index[0], index[1], index[2]);
          cell *npt    = grid->NextPt(cpt, posdirs[idim]);
          cell *lpt    = 0;
          double q_neg = 0.0, q_pos = 0.0, gradT = 0.0, Qclassical = 0.0,
                 Qsaturated = 0.0, T = 0.0;
          double dx    = grid->DX();
          double kappa = 0.0, vc = 0.0, cn = 0.0;
          if (npt == 0) rep.error("Couldn't find two cells in column", 0);
          q_neg = 0.0;  // no flux coming in from non-existent boundary data.
          q_pos = 0.0;
          // Run through column, calculating conductive flux
          do {
            // Now use the Slavin & Cox (1992) formula for conduction to
            // get the conductive heat flux from cpt to npt in direction
            // posdir[idim].
            gradT = (npt->dU[RHO] - cpt->dU[RHO])
                    / (grid->idifference_cell2cell(cpt, npt, axis[idim])
                       * CI.phys_per_int());
            // If flow is from npt to cpt, we use npt's values for
            // calculating Q. Else we use cpt's values. (if
            // gradT>0, then T2>T1, flow from 2->1 in the *negative*
            // direction).
            if (gradT > 0.0) {
              c = npt;
              // vc = cpt->Ph[VX + idim];
              // cn = spatial_solver->maxspeed(npt->Ph, par.gamma);
            }
            else {
              c = cpt;
              // vc = -npt->Ph[VX + idim];
              // cn = spatial_solver->maxspeed(cpt->Ph, par.gamma);
            }
            // get ln(Lambda) and then the classical and
            // saturated fluxes. For ln(Lambda) the formula is only
            // valid for T>4.2e5.  The value of 4.2735e23 is
            // (1.4*m_p)^{-1}.
            T = c->dU[RHO];
            // if (T < 4.2e5)
            //  Qclassical = 29.7;
            // else
            //  Qclassical = 29.7 + log(T / (1.0e6 * sqrt(c->Ph[RO]
            //  * 4.2735e23)));
            // Qclassical = -1.84e-5 * pow(T, 2.5) * gradT / Qclassical;
            // kappa =
            //    max(0.0, max(6e-7, 6.0e-7 * pow(T, 2.5) * exp(-5e3 / T))
            //                 * (1.0 - vc / cn));
            // Qclassical = -kappa * gradT;

            // For saturated Q we follow S&C(1992) and use phi_s=0.3
            Qsaturated = -1.5 * c->Ph[RO] * pow(c->Ph[PG] / c->Ph[RO], 1.5);
            if (gradT < 0.0) Qsaturated *= -1.0;
            // Qsaturated *= max(0.0, min(1.0, 1.0 - 1.0 * vc / cn)); //
            // upwinding?
            Qsaturated *= par.EP.tc_strength;
            // now Q = Qs*(1-exp(-Qc/Qs)).
            //  - (Qs>>Qc)=>(Q->Qc).
            //  - (Qc>>Qs)=>(Q->Qs).
            q_pos = Qsaturated;  // * (1.0 - exp(-Qclassical / Qsaturated));
            if (!c->isdomain) q_pos = 0.0;
            // Finally cpt needs an updated -div(q) value from the
            // current direction. This is a hack, because there is no
            // VectorOps function to do this for me. I should write a
            // function in VectorOps which I can call to do this... I
            // should also make the base grid derive from
            // base-VectorOps, so that the functions are accessible!
            if (par.coord_sys == COORD_CYL && axis[idim] == Rcyl) {
              double rp = CI.get_dpos(c, Rcyl) + 0.5 * dx;
              double rn = rp - dx;
              cpt->dU[ERG] +=
                  2.0 * (rn * q_neg - rp * q_pos) / (rp * rp - rn * rn);
            }
            else if (par.coord_sys == COORD_SPH && axis[idim] == Rsph) {
              double rc = CI.get_dpos(c, Rsph);
              double rp = rc + 0.5 * dx;
              double rn = rp - dx;
              rc        = (pow(rp, 3.0) - pow(rn, 3.0)) / 3.0;
              cpt->dU[ERG] += (rn * rn * q_neg - rp * rp * q_pos) / rc;
            }
            else {
              cpt->dU[ERG] += (q_neg - q_pos) / dx;
            }
            // cout <<"\tQc="<<Qclassical<<", Qs="<<Qsaturated<<",
            // Q="<<q_pos<<", Edot="<<cpt->dU[ERG]<<"\n";

            // set diffusion coefficient D = kappa /(n kB)
            if (oa == OA1)
              // saturated flux is hyperbolic, so multiply by dx to remove a
              // factor of dx in the dt calculation
              if (1 == 1)  //(Qsaturated >= Qclassical)
                cpt->dU[VX] = max(
                    cpt->dU[VX],
                    2.0 * spatial_solver->maxspeed(c->Ph, par.gamma) * par.dx);
              else  // classical flux --> diffusion equation
                cpt->dU[VX] = max(
                    cpt->dU[VX], kappa * 2.3e-24 / (c->Ph[RO] * pconst.kB()));

            // Set npt to cpt, set current q_pos to q_neg for next cell.
            // Move to next cell.
            q_neg = q_pos;
            lpt   = cpt;
            cpt   = npt;
          } while ((npt = grid->NextPt(npt, posdirs[idim])) != 0);
          cpt->dU[VX] = max(cpt->dU[VX], lpt->dU[VX]);
        }  // ax1
      }    // ax2
#ifdef PION_OMP
    }  // OMP loop
#endif
  }  // Loop over three directions.

#ifdef PION_OMP
  #pragma omp parallel
  {
    spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
  }
#else
  spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
#endif  // PION_OMP

  //
  // Now reset the temporary storage of Temperature in dU[RHO] to zero.
  //
  c = grid->FirstPt_All();
  do {
    c->dU[RHO] = 0.0;
  } while ((c = grid->NextPt_All(c)) != 0);

  return 0;
}



// ##################################################################
// ##################################################################



void calc_timestep::timestep_checking_and_limiting(
    class SimParams &par,  ///< pointer to simulation parameters
    const int l            ///< level to advance (for NG grid)
)
{
  //
  // If the timestep is less than the minimum allowed, report an error
  // and bug out!  This must be before we check for the next output
  // time because that can limit the step arbitrarily.
  //
  if (par.dt < par.min_timestep) {
    ostringstream temp;
    temp.str("");
    temp << "Timestep too short! dt=" << par.dt
         << "  min-step=" << par.min_timestep;
    rep.error(temp.str(), par.dt);
  }

  //
  // So that we don't increase the timestep by more than 30% over last step:
  //
#ifdef TIMESTEP_LIMITING
  par.dt = min(par.dt, 1.3 * par.levels[l].last_dt);
#endif  // TIMESTEP_LIMITING

  //
  // If we are outputting every n-years, then check if we need to adjust dt.
  //
  if (par.op_criterion == 1) {
    par.dt = min(par.dt, par.next_optime - par.simtime);
    if (par.dt <= 0.0)
      rep.error("Went past output time without outputting!", par.dt);
  }

  //
  // Make sure we end up exactly at finishtime:
  //
  par.dt = min(par.dt, par.finishtime - par.simtime);
  if (par.dt <= 0.0) {
    cout << "dt=" << par.dt << ", finish=" << par.finishtime << ", now=";
    cout << par.simtime << "\n";
    rep.error("Negative timestep!", par.dt);
  }

  return;
}



// ##################################################################
// ##################################################################



double calc_timestep::calc_dynamics_dt(
    class SimParams &par,  ///< pointer to simulation parameters
    class GridBaseClass *grid)
{
  double tempdt = 0.0;
  double dt     = 1.e100;  // Set it to very large no. initially.
  double dx     = grid->DX();

  class cell *c = grid->FirstPt_All();
#ifndef NDEBUG
  dp.c = c;
#endif

  //
  // Now go through all of the cells on the local grid.
  // The CellTimeStep() function returns a value which is
  // already multiplied by the CFL coefficient.
  //
  // cout << "timestep!\n";
  int index[3];
  tempdt = 1.0e100;
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        tempdt   = 1.0e100;
        // get to first non-boundary cell:
        while (c && !c->isgd)
          c = grid->NextPt(c, XP);
        if (!c) {
          continue;  // end iteration if there are no non-boundary cells
        }
        // loop over all cells in this x-column.
        do {
          if (c->timestep && !c->isbd) {
            tempdt =
                min(tempdt, spatial_solver->CellTimeStep(c, par.gamma, dx));
          }
        } while ((c = grid->NextPt(c, XP)));
        dt = min(dt, tempdt);
      }
    }
#ifdef PION_OMP
  }
#endif

  // if on first step and a jet is inflowing, then set dt based on
  // the jet properties.
  if (par.timestep == 0 && JP.jetic != 0) {
    dt = std::min(dt, 0.1 * par.CFL * grid->DX() / JP.jetstate[VX]);
  }
  // if on first step and stellar winds are present, limit dt based
  // on the wind speed (convert from km/s to cm/s)
  if (par.timestep == 0 && SWP.Nsources != 0) {
    for (int v = 0; v < SWP.Nsources; v++)
      dt = std::min(
          dt, 0.1 * par.CFL * grid->DX() / (SWP.params[v]->Vinf * 1.0e5));
  }

  if (dt <= 0.0) rep.error("Got zero timestep!!!", dt);
#ifndef NDEBUG
  cout << "(calc_dynamics_dt)  min-dt=" << dt << "\n";
#endif

  return dt;
}



// ##################################################################
// ##################################################################



double calc_timestep::calc_microphysics_dt(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level to advance (for NG grid)
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
  if (par.EP.MP_timestep_limit == 0) {
    return 1.0e99;
  }

  double dt = -1.0e99;  // initialise to negative value

  //
  // We see if we need to use the column densities and source properties in
  // the timestep calculation, or if we just use the local microphysical
  // cooling and/or recombination times.
  //
  if (FVI_need_column_densities_4dt) {
    //
    // need column densities, so do raytracing, and then get dt.
    //
    // cout <<"calc_timestep, getting column densities
    // rt="<<par.RS.Nsources<<".\n"; int err = RT_all_sources(par,grid,l);
    dt = get_mp_timescales_with_radiation(par, grid);
    if (dt <= 0.0)
      rep.error("get_mp_timescales_with_radiation() returned error", dt);
  }
  else {
    //
    // don't need column densities, so call the no-RT version
    //
    // cout <<" getting timestep with no radiation\n";
    dt = get_mp_timescales_no_radiation(par, grid);
    if (dt <= 0.0)
      rep.error("get_mp_timescales_no_radiation() returned error", dt);
  }

#ifndef NDEBUG
  cout << "(calc_microphysics_dt)  min-dt=" << dt << "\n";
#endif

  return dt;
}



// ##################################################################
// ##################################################################



double calc_timestep::get_mp_timescales_no_radiation(
    class SimParams &par,  ///< pointer to simulation parameters
    class GridBaseClass *grid)
{
#ifndef NDEBUG
  // paranoid checking...
  if (par.EP.MP_timestep_limit == 0) {
    cout << "calc_timestep::get_mp_timescales_no_radiation() called, ";
    cout << "but no MP-dt limiting!\n";
    return -1.0;
  }
#endif  // NDEBUG
  //
  // Loop through every cell and see what the limit is.
  //
  double tempdt = 0.0, dt = 1.0e99, t = 0.0;
  class cell *c = grid->FirstPt_All();
  int index[3];
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index,t) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        tempdt   = 1.0e100;
        // get to first non-boundary cell:
        while (c && !c->isgd)
          c = grid->NextPt(c, XP);
        if (!c) {
          continue;  // end iteration if there are no non-boundary cells
        }
        do {
#ifndef NDEBUG
          dp.c = c;
#endif
          // Check if cell is boundary data or not (can only be an internal
          // boundary, such as a stellar wind, since we are looping over cells
          // which are grid data)  If it is boundary data then we skip it.
          if (c->isbd || !c->isleaf) {
#ifndef NDEBUG
            cout << "skipping cell " << c->id
                 << " in get_mp_timescales_no_radiation() c->isbd.\n";
#endif
          }
          else {
            //
            // timescales(state-vec, gamma, t_cool?, t_rec?, t_photoion?)
            //
            switch (par.EP.MP_timestep_limit) {
              case 1:  // cooling time only.
                t = MP->timescales(c->Ph, par.gamma, true, false, false);
                break;
              case 4:  // recomb only
                t = MP->timescales(c->Ph, par.gamma, false, true, false);
                break;
              case 2:  // cooling+recomb
                t = MP->timescales(c->Ph, par.gamma, true, true, false);
                break;
              case 3:  // cooling+recomb+ionisation (not working!)
                t = MP->timescales(c->Ph, par.gamma, true, true, true);
                break;
              default:
                rep.error("Bad MP_timestep_limit", par.EP.MP_timestep_limit);
            }
            tempdt = min(tempdt, t);
            // cout <<"(get_min_timestep) i ="<<i<<"  min-dt="<<dt<<"\n";
          }  // if not boundary data.
        } while ((c = grid->NextPt(c, XP)));
        dt = min(dt, tempdt);
      }
    }
#ifdef PION_OMP
  }
#endif

#ifndef RT_TEST_PROBS
  //
  // If doing photo-ionisation, can underestimate 1st timestep because gas is
  // all neutral initially, but only for a few seconds!
  // (DON'T WANT TO SET THIS FOR NON-DYNAMICS TEST PROBLEMS)
  //
  if ((par.timestep < 3) && (par.EP.raytracing > 0)
      && (par.EP.phot_ionisation)) {
    //
    // adjust first timestep so that it corresponds to ionised cell.
    //
    dt = min(dt, grid->DX() * par.CFL / 2.0e6);  // 20 km/s wavespeed.
    //
    // Also make sure it is not larger than 1.0e6 seconds (0.03 years).
    // With photoionisation we must be using CGS units, so it is ok to
    // hardcode the time.
    //
    dt = min(dt, 1.0e7);
#ifdef DEBUG_MP
    cout << "\tRT timestep: \t\t\tdt=" << dt << "\n";
#endif
  }
  //
  // Make sure the first timestep is short if doing RT, here set to be
  // 0.3333* the recombination time.
  // approx = 0.33333/(2.59e-13*rho/2.338e-24) = 3.009e-12/rho
  // Since we must have microphysics set up here, it is safe to assume
  // the code is using cgs units for density.
  //
  if ((par.timestep == 0) && (par.RS.Nsources > 0)) {
    c = grid->FirstPt();
#ifdef DEBUG_MP
    cout << "rho=" << c->Ph[RO] << ", old dt=" << dt;
#endif
    dt = min(dt, 3.009e-12 / c->Ph[RO]);
#ifdef DEBUG_MP
    cout << ", updated dt=" << dt << "\n";
#endif
  }
#endif  // RT_TEST_PROBS

  return dt;
}



// ##################################################################
// ##################################################################



double calc_timestep::get_mp_timescales_with_radiation(
    class SimParams &par,  ///< pointer to simulation parameters
    class GridBaseClass *grid)
{
#ifndef NDEBUG
  //
  // paranoid checking...
  //
  if (par.EP.MP_timestep_limit == 0) {
    cout << "calc_timestep::get_mp_timescales_with_radiation() no MP-dt "
            "limiting\n";
    return -1.0;
  }
  if (par.RS.Nsources == 0)
    rep.error(
        "calc_timestep::get_mp_timescales_with_radiation() no sources", 1);
#endif  // NDEBUG

  //
  // RT source properties are already in structs for the microphysics calls.
  // So now we need to go through every cell and see what the limit is.
  //
  double tempdt = 0.0, dt = 1.0e99, t = 0.0;
  class cell *c = grid->FirstPt_All();
  int index[3];
  // copies of the source data (one copy for each thread)
  // seems to be memory issues when using class member data.
  vector<struct rt_source_data> heating;
  vector<struct rt_source_data> ionize;

#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index,t,heating,ionize) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        heating  = FVI_heating_srcs;
        ionize   = FVI_ionising_srcs;
        tempdt   = 1.0e100;
        // get to first non-boundary cell:
        while (c && !c->isgd)
          c = grid->NextPt(c, XP);
        if (!c) {
          continue;  // end iteration if there are no non-boundary cells
        }
        do {
#ifndef NDEBUG
          dp.c = c;
#endif
          // Check if cell is boundary data (can only be an internal
          // boundary, such as a stellar wind, since we are looping over
          // cells which are grid data).  If it is boundary data then we
          // skip it.
          if (c->isbd || !c->isleaf) {
#ifndef NDEBUG
            cout << "skipping cell " << c->id
                 << " in get_mp_timescales_with_radiation() c->isbd.\n";
#endif
          }
          else {
            //
            // Get column densities and Vshell in struct for each source.
            //
            for (int v = 0; v < FVI_nheat; v++) {
              heating[v].Vshell = CI.get_cell_Vshell(c, heating[v].id);
              heating[v].dS     = CI.get_cell_deltaS(c, heating[v].id);
              CI.get_cell_col(c, heating[v].id, heating[v].DelCol);
              CI.get_col(c, heating[v].id, heating[v].Column);
              for (short unsigned int iC = 0; iC < heating[v].NTau; iC++)
                heating[v].Column[iC] -= heating[v].DelCol[iC];
            }
            for (int v = 0; v < FVI_nion; v++) {
              ionize[v].Vshell = CI.get_cell_Vshell(c, ionize[v].id);
              ionize[v].dS     = CI.get_cell_deltaS(c, ionize[v].id);
              CI.get_cell_col(c, ionize[v].id, ionize[v].DelCol);
              CI.get_col(c, ionize[v].id, ionize[v].Column);
              for (short unsigned int iC = 0; iC < ionize[v].NTau; iC++)
                ionize[v].Column[iC] -= ionize[v].DelCol[iC];
              if (ionize[v].Column[0] < 0.0) {
#ifdef RT_TESTING
                cout << "dx=" << grid->DX() << "  ";
                CI.print_cell(c);
                rep.error("time_int:calc_RT_microphysics_dU tau<0", 1);
#endif
                for (short unsigned int iC = 0; iC < ionize[v].NTau; iC++)
                  ionize[v].Column[iC] = max(0.0, ionize[v].Column[iC]);
              }
            }
            //
            // We assume we want to limit by all relevant timescales.
            //
            t = MP->timescales_RT(
                c->Ph, FVI_nheat, heating, FVI_nion, ionize, par.gamma);
#ifndef NDEBUG
            if (t < dt) {
              cout << "(get_min_timestep) id=" << c->id << ":  dt=" << tempdt
                   << ", min-dt=" << dt;
              cout << ".\t 1-x=" << 1.0 - c->Ph[par.ftr] << ", pg=" << c->Ph[PG]
                   << "\n";
            }
#endif
            if (t <= 0.0) {
              cout
                  << "get_mp_timescales_with_radiation() negative timestep... ";
              cout << "c->id=" << c->id << "\tdt=" << t << "\n";
              rep.printVec("Ph", c->Ph, par.nvar);
              rep.error("Negative timestep from microphysics with RT!", t);
            }
            tempdt = min(tempdt, t);
            // cout <<"(get_min_timestep) i ="<<i<<"  min-dt="<<dt<<"\n";
          }  // if not boundary data.
        } while ((c = grid->NextPt(c, XP)));
        dt = min(dt, tempdt);
      }
    }
#ifdef PION_OMP
  }
#endif

  return dt;
}
