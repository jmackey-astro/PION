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

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

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
  spdlog::debug(
      "calc_timestep::calc_timestep(): g={}, rt={}", fmt::ptr(grid),
      fmt::ptr(grid->RT));
#endif
  //
  // This is a wrapper function.  First we get the dynamics
  // timestep, and then the microphysics timestep.
  //
  double t_dyn = 0.0, t_mp = 0.0;
  t_dyn = calc_dynamics_dt(par, grid);
  t_mp  = calc_microphysics_dt(par, grid);
  // cout <<"l="<<l<<", \t t_dyn="<<t_dyn<<"and t_mp ="<<t_mp<<"\n";

#ifndef NDEBUG
  if (t_mp < t_dyn)
    spdlog::debug(
        "Limiting timestep by MP: mp_t={}\thydro_t={}t_dyn = {}  and t_mp = {}",
        t_mp, t_dyn, t_dyn, t_mp);
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
    spdlog::debug(
        "CONDUCTION IS LIMITING TIMESTEP: t_c={}, t_m={}, t_dyn={}", t_cond,
        t_mp, t_dyn);
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
  spdlog::info("calc_timestep::calc_timestep() finished");
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
  // First we need to set Edot() in every cell.  This is stored in
  // c->dU[ERG] since this is where it will be updated.
  //
  // cout <<"\tCCdt: setting Edot.\n";
  set_thermal_conduction_Edot(par, OA1, grid);
  // cout <<"\tCCdt: Done with div(Q).  Now getting timestep.\n";

  //
  // Now run through the grid and calculate E_int/Edot to get the smallest
  // timestep we can allow.  Also reset dU[VX] to zero (was used as temporary
  // variable).
  //
  double dt = 1.0e200, tempdt = 0.0;
  enum axes x1 = XX;
  enum axes x2 = YY;
  enum axes x3 = ZZ;
  int nx2      = grid->NG_All(x2);
  int nx3      = grid->NG_All(x3);
  // loop over the two perpendicular axes, to trace out a plane of
  // starting cells for calculating fluxes along columns along this
  // axis.  The plane can be a single cell
  // (in 1D) or a line (in 2D) or a plane (in 3D).
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(tempdt) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < nx3; ax3++) {
      for (int ax2 = 0; ax2 < nx2; ax2++) {
        int index[3];
        index[x1] = 0;
        index[x2] = static_cast<int>(ax2);
        index[x3] = static_cast<int>(ax3);
        cell *c   = grid->get_cell_all(index[0], index[1], index[2]);
        tempdt    = 1.0e100;
        do {
          // This is appropriate for saturated TC, where dU[VX] is set to
          // the max speed temporarily in set_thermal_conduction_Edot().
          tempdt    = min(tempdt, par.CFL * par.dx / c->dU[VX]);
          c->dU[VX] = 0.0;
        } while ((c = grid->NextPt(c, XP)) != 0);
        dt = min(dt, tempdt);
      }
    }
#ifdef PION_OMP
  }
#endif

  // first timestep can be dodgy with conduction and stellar winds,
  // so set it to a very small value
  if (par.timestep == 0) dt = min(dt, 1.0e0);
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
  // First we need to calculate the temperature in every cell.  This
  // is stored in dU[RHO] -- reset to zero at the end of function.
  //
  if (!MP)
    spdlog::error("{}: {}", "No Microphysics == No Conduction", fmt::ptr(MP));
  enum axes x1 = XX;
  enum axes x2 = YY;
  enum axes x3 = ZZ;
  int nx2      = grid->NG_All(x2);
  int nx3      = grid->NG_All(x3);
#ifdef PION_OMP
  #pragma omp parallel for collapse(2)
#endif
  for (int ax3 = 0; ax3 < nx3; ax3++) {
    for (int ax2 = 0; ax2 < nx2; ax2++) {
      int index[3];
      index[x1] = 0;
      index[x2] = static_cast<int>(ax2);
      index[x3] = static_cast<int>(ax3);
      cell *cpt = grid->get_cell_all(index[0], index[1], index[2]);
      do {
        cpt->dU[RHO] = MP->Temperature(cpt->Ph, par.gamma);
        cpt->dU[ERG] = 0.0;
      } while ((cpt = grid->NextPt(cpt, XP)) != 0);
    }
  }

  // cout <<"\tT calculated, now calculating divQ.\n";
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
    spdlog::debug("\t\t\tidim={}", idim);
#endif  // NDEBUG

#ifdef PION_OMP
    #pragma omp parallel
    {
      spatial_solver->SetDirection(axis[idim]);
    }
#else
    spatial_solver->SetDirection(axis[idim]);
#endif  // PION_OMP

#ifdef TEST_INT
    spdlog::debug("Direction={}, i={}", axis[idim], idim);
#endif
    // loop over the number of cells in the line/plane of starting
    // cells.
    x1      = axis[(idim + 1) % 3];
    x2      = axis[(idim + 2) % 3];
    x3      = axis[idim];
    nx2     = grid->NG_All(x2);
    int nx1 = grid->NG_All(x1);
    // loop over the two perpendicular axes, to trace out a plane of
    // starting cells for calculating fluxes along columns along this
    // axis.  The plane can be a single cell
    // (in 1D) or a line (in 2D) or a plane (in 3D).
#ifdef PION_OMP
    #pragma omp parallel
    {
      #pragma omp for collapse(2)
#endif
      for (int ax2 = 0; ax2 < nx2; ax2++) {
        for (int ax1 = 0; ax1 < nx1; ax1++) {
          int index[3];
          cell *cpt, *c2;
          index[x1]    = static_cast<int>(ax1);
          index[x2]    = static_cast<int>(ax2);
          index[x3]    = 0;
          cpt          = grid->get_cell_all(index[0], index[1], index[2]);
          cell *npt    = grid->NextPt(cpt, posdirs[idim]);
          cell *lpt    = 0;
          double q_neg = 0.0, q_pos = 0.0, gradT = 0.0, T = 0.0;
          double q1 = 0.0, q2 = 0.0;
          double dx = grid->DX();
          if (npt == 0)
            spdlog::error("{}: {}", "Couldn't find two cells in column", 0);
          q_neg = 0.0;  // no flux coming in from non-existent boundary data.
          q_pos = 0.0;
          // Run through column, calculating conductive flux following
          // Slavin & Cox (1992)
          do {
            // Calculate saturated heat flux from cpt to npt in direction
            // posdir[idim].  Only depends on sign of gradT, so no need
            // to get normalisation correct with denominator:
            gradT = (npt->dU[RHO] - cpt->dU[RHO]);
            // if gradT>0, then T2>T1, flow from 2->1 in the *negative*
            // direction.
            if (gradT > 0.0)
              c2 = npt;
            else
              c2 = cpt;
            T = c2->dU[RHO];
            // For saturated Q we follow S&C(1992) and use phi_s=0.3
            q_pos = 1.5 * c2->Ph[RO] * pow(c2->Ph[PG] / c2->Ph[RO], 1.5);
            if (gradT > 0.0) q_pos *= -1.0;
            q_pos *= par.EP.tc_strength;
            if (!c2->isdomain) q_pos = 0.0;  // don't update wind cells
            // Finally cpt needs an updated -div(q) value from the
            // current direction.
            if (par.coord_sys == COORD_CYL && axis[idim] == Rcyl) {
              double rp = CI.get_dpos(cpt, Rcyl) + 0.5 * dx;
              double rn = rp - dx;
              cpt->dU[ERG] +=
                  2.0 * (rn * q_neg - rp * q_pos) / (rp * rp - rn * rn);
            }
            else if (par.coord_sys == COORD_SPH && axis[idim] == Rsph) {
              double rc = CI.get_dpos(cpt, Rsph);
              double rp = rc + 0.5 * dx;
              double rn = rp - dx;
              rc        = (pow(rp, 3.0) - pow(rn, 3.0)) / 3.0;
              cpt->dU[ERG] += (rn * rn * q_neg - rp * rp * q_pos) / rc;
            }
            else {
              cpt->dU[ERG] += (q_neg - q_pos) / dx;
            }

            // set max speed in dU[VX] to send back to the timestep calculation,
            // if we are on the half-step update.
            if (oa == OA1) {
              // saturated flux is hyperbolic
              cpt->dU[VX] =
                  max(cpt->dU[VX],
                      4.0 * spatial_solver->maxspeed(c2->Ph, par.gamma));
            }
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

  // Reset the temporary storage of Temperature in dU[RHO] to zero.
  x1  = XX;
  x2  = YY;
  x3  = ZZ;
  nx2 = grid->NG_All(x2);
  nx3 = grid->NG_All(x3);
#ifdef PION_OMP
  #pragma omp parallel for collapse(2)
#endif
  for (int ax3 = 0; ax3 < nx3; ax3++) {
    for (int ax2 = 0; ax2 < nx2; ax2++) {
      int index[3];
      index[x1] = 0;
      index[x2] = static_cast<int>(ax2);
      index[x3] = static_cast<int>(ax3);
      cell *cpt = grid->get_cell_all(index[0], index[1], index[2]);
      do {
        cpt->dU[RHO] = 0.0;
      } while ((cpt = grid->NextPt(cpt, XP)) != 0);
    }
  }

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
    spdlog::error("{}: {}", temp.str(), par.dt);
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
      spdlog::error(
          "{}: {}", "Went past output time without outputting!", par.dt);
  }

  //
  // Make sure we end up exactly at finishtime:
  //
  par.dt = min(par.dt, par.finishtime - par.simtime);
  if (par.dt <= 0.0) {
    spdlog::debug(
        "dt={}, finish={}, now={}", par.dt, par.finishtime, par.simtime);
    spdlog::error("{}: {}", "Negative timestep!", par.dt);
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
  tempdt  = 1.0e100;
  int nx2 = grid->NG_All(YY);
  int nx3 = grid->NG_All(ZZ);
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < nx3; ax3++) {
      for (int ax2 = 0; ax2 < nx2; ax2++) {
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

  if (dt <= 0.0) spdlog::error("{}: {}", "Got zero timestep!!!", dt);
#ifndef NDEBUG
  spdlog::debug("(calc_dynamics_dt)  min-dt={}", dt);
#endif

  return dt;
}



// ##################################################################
// ##################################################################



double calc_timestep::calc_microphysics_dt(
    class SimParams &par,      ///< pointer to simulation parameters
    class GridBaseClass *grid  ///< pointer to grid.
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
      spdlog::error(
          "{}: {}", "get_mp_timescales_with_radiation() returned error", dt);
  }
  else {
    //
    // don't need column densities, so call the no-RT version
    //
    // cout <<" getting timestep with no radiation\n";
    dt = get_mp_timescales_no_radiation(par, grid);
    if (dt <= 0.0)
      spdlog::error(
          "{}: {}", "get_mp_timescales_no_radiation() returned error", dt);
  }

#ifndef NDEBUG
  spdlog::debug("(calc_microphysics_dt)  min-dt={}", dt);
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
    spdlog::error(
        "calc_timestep::get_mp_timescales_no_radiation() called, but no MP-dt limiting!\n");
    return -1.0;
  }
#endif  // NDEBUG
  //
  // Loop through every cell and see what the limit is.
  //
  double tempdt = 0.0, dt = 1.0e99, t = 0.0;
  class cell *c = grid->FirstPt_All();
  int index[3];
  int nx2 = grid->NG_All(YY);
  int nx3 = grid->NG_All(ZZ);
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index,t) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < nx3; ax3++) {
      for (int ax2 = 0; ax2 < nx2; ax2++) {
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
            spdlog::debug(
                "skipping cell {} in get_mp_timescales_no_radiation() c->isbd",
                c->id);
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
                spdlog::error(
                    "{}: {}", "Bad MP_timestep_limit",
                    par.EP.MP_timestep_limit);
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
    spdlog::debug("\tRT timestep: \t\t\tdt={}", dt);
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
    spdlog::debug("rho={}, old dt={}", c->Ph[RO], dt);
#endif
    dt = min(dt, 3.009e-12 / c->Ph[RO]);
#ifdef DEBUG_MP
    spdlog::debug(", updated dt={}", dt);
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
    spdlog::error(
        "calc_timestep::get_mp_timescales_with_radiation() no MP-dt limiting");
    return -1.0;
  }
  if (par.RS.Nsources == 0)
    spdlog::error(
        "{}: {}",
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
  int nx2 = grid->NG_All(YY);
  int nx3 = grid->NG_All(ZZ);

#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,tempdt,index,t,heating,ionize) reduction(min:dt)
#endif
    for (int ax3 = 0; ax3 < nx3; ax3++) {
      for (int ax2 = 0; ax2 < nx2; ax2++) {
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
            spdlog::debug(
                "skipping cell {} in get_mp_timescales_with_radiation() c->isbd.\n",
                c->id);
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
                spdlog::debug("dx={} ", grid->DX());
                CI.print_cell(c);
                spdlog::error(
                    "{}: {}", "time_int:calc_RT_microphysics_dU tau<0", 1);
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
              spdlog::debug(
                  "(get_min_timestep) id={}:  dt={}, min-dt={}.\t 1-x={}, pg={}",
                  c->id, tempdt, dt, 1.0 - c->Ph[par.ftr], c->Ph[PG]);
            }
#endif
            if (t <= 0.0) {
              spdlog::error(
                  "get_mp_timescales_with_radiation() negative timestep... c->id={}\tdt={}",
                  c->id, t);
              spdlog::debug(
                  "Ph : {}", std::vector<double>(c->Ph, c->Ph + par.nvar));
              spdlog::error(
                  "{}: {}", "Negative timestep from microphysics with RT!", t);
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
