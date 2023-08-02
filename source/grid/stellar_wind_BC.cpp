///
/// \file stellar_wind_BC.cpp
/// \author Jonathan Mackey
/// \date 2010.10.05
///
/// Stellar wind boundary condition class.
/// Moved from global.cc b/c it got too big.
///
/// Modifications:
///
/// - 2010.10.05 JM: Added spherical coordinate possibility.  Also
///    changed pressure calculation to calculate it at the stellar
///    radius and scale outwards assuming adiabatic expansion.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  So VX,VY,VZ are now set by using calls to
///   distance() functions in the grid class.  It doesn't seem to help
///   much.
/// - 2011.01.07 JM: I debugged the geometric grid functions, and now
///   it works very well!  I have a nice spherical expansion.
/// - 2011.01.18 JM: Added #def which sets pressure so that Tmin=10K
/// - 2011.02.14 JM: Added stellar_wind_evolution class for winds with
///    evolving properties, determined by a stellar evolution model, and
///    fitted with spline functions.
///    02.15 JM: Debugged. 02.16 JM: Debugged
/// - 2011.04.29 JM: Now in add_cell(), the c->isbd bool is set to true to
///    indicate that it has become boundary data (microphysics updates
///    will skip it in this case and, more importantly, microphysics
///    timescales calculations.
/// - 2011.06.20 JM: Got rid of non-ANSI-C exp10 functions
/// - 2011.11.22 JM: Added t_scalefactor parameter for stellar winds.
/// - 2011.12.01 JM: Switched from spline to linear interpolation for
///    winds.
/// - 2012.12.07/10 JM: Changed min. ion frac. in wind from 0 to 1e-7.
/// - 2013.04.15 JM: removed lots of comments (or put in TESTING def)
/// - 2013.04.16 JM: Fixed bug where Set_Temp() was called when
///    tracer variables were still (potentially) unset in wind cells.
/// - 2013.08.19 JM: got rid of cm_per_km() function.
/// - 2013.09.06 JM: removed the slowly-expanding switch-on wind
///    because it didn't help with anything, ever.
///    Removed the integer positions because they created potential
///    errors in the wind properties from rounding errors.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.10.19 JM: Fixed wind-tracer to always use pion_flt.
/// - 2017.07.26 JM: cleaned up code a bit.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG

#include "grid/grid_base_class.h"
#include "grid/stellar_wind_BC.h"
#include "microphysics/microphysics_base.h"
#include <array>
#include <sstream>
using namespace std;

//

// ##################################################################
// ##################################################################

stellar_wind::stellar_wind(
    const int nd,           ///< ndim
    const int nv,           ///< nvar
    const int nt,           ///< ntracer
    const int ft,           ///< ftr
    const std::string *tr,  ///< List of tracer variable names.
    const int cs,           ///< coord_sys
    const int eq,           ///< eqn_type
    const double mt         ///< Minimum temperature allowed on grid
    ) :
    ndim(nd),
    nvar(nv), ntracer(nt), ftr(ft), coordsys(cs), eqntype(eq), Tmin(mt)
{
  nsrc = 0;
  MP   = 0;

  for (int v = 0; v < nt; v++)
    tracers.push_back(tr[v]);
}



// ##################################################################
// ##################################################################



stellar_wind::~stellar_wind()
{
  //
  // Need to delete the wind_cell structs in each wind_source struct,
  // and then delete the wind_source structs.
  //
  struct wind_source *ws;

  for (int n = 0; n < nsrc; n++) {
    ws = wlist[n];
    for (int m = 0; m < ws->ncell; m++) {
      ws->wcells[m]->p = mem.myfree(ws->wcells[m]->p);
      ws->wcells[m]    = mem.myfree(ws->wcells[m]);
    }
    ws->wcells.clear();
    ws = mem.myfree(ws);
  }
  return;
}



// ##################################################################
// ##################################################################



int stellar_wind::add_source(
    struct stellarwind_params *wp  ///< pointer to wind parameters struct
)
{
  struct wind_source *ws = 0;
  ws                     = mem.myalloc(ws, 1);
  ws->pars               = wp;
  ws->ncell              = 0;
  switch (wp->type) {
    case WINDTYPE_CONSTANT:
    case WINDTYPE_EVOLVING:
      break;
    default:
      spdlog::error(
          "{}: {}", "What type of source is this?  add a new type?", wp->type);
      break;
  }

#ifdef ANALYTIC_ORBITS
  // calculate initial velocity:
  double r = sqrt(
      wp->PeriastronX * wp->PeriastronX + wp->PeriastronY * wp->PeriastronY);
  double v =
      2 * pconst.pi() * sqrt(1.0 - wp->eccentricity * wp->eccentricity)
      / (wp->OrbPeriod * (1.0 - wp->eccentricity) * (1.0 - wp->eccentricity));
  wp->velocity[XX] = v * wp->PeriastronY;
  wp->velocity[YY] = -v * wp->PeriastronX;
  wp->velocity[ZZ] = 0.0;
  for (int v = 0; v < ndim; v++)
    ws->dpos_init[v] = wp->dpos[v];
  for (int v = ndim; v < MAX_DIM; v++)
    ws->dpos_init[v] = 0.0;
#endif  // ANALYTIC_ORBITS

  // if using microphysics, find H+ tracer variable, if it exists.
  ws->Hplus  = -1;
  ws->iHplus = -1;
  int hplus  = -1;
  if (MP) {
    hplus = MP->Tr("H1+");
  }
  ws->Hplus = hplus;
  if (hplus >= 0) ws->iHplus = hplus - nvar + ntracer;

  ws->cells_added = false;
  if (!ws->wcells.empty())
    spdlog::error(
        "{}: {}", "wind_source: wcells not empty!", ws->wcells.size());

  // Make sure the source position is compatible with the geometry:
  if (coordsys == COORD_SPH) {
    if (!pconst.equalD(wp->dpos[Rsph], 0.0))
      spdlog::error(
          "{}: {}", "Spherical symmetry but source not at origin!",
          wp->dpos[Rsph]);
  }
  if (coordsys == COORD_CYL && ndim == 2) {
    if (!pconst.equalD(wp->dpos[Rcyl], 0.0))
      spdlog::error(
          "{}: {}", "Axisymmetry but source not at R=0!", wp->dpos[Rcyl]);
  }

  // initialise mypos for moving sources
  for (int v = 0; v < ndim; v++)
    ws->mypos[v] = wp->dpos[v];
  for (int v = ndim; v < MAX_DIM; v++)
    ws->mypos[v] = 0.0;

  wlist.push_back(ws);
  nsrc++;
  spdlog::debug(
      "\tAdded wind source id={} to list of {} elements", nsrc - 1, nsrc);
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind::Nsources()
{
  return nsrc;
}



// ##################################################################
// ##################################################################



int stellar_wind::add_cell(
    class GridBaseClass *grid,
    const int id,  ///< src id
    cell &c        ///< cell to add to list.
)
{
  if (id < 0 || id >= nsrc) spdlog::error("{}: {}", "bad src id", id);
  struct wind_source *WS        = wlist[id];
  struct stellarwind_params *WP = wlist[id]->pars;

  //
  // Setup a wind_cell struct
  //
  array<double, MAX_DIM> cpos, wpos;
  double dx = grid->DX();
  CI.get_dpos(c, cpos);
  for (int v = 0; v < MAX_DIM; v++)
    wpos[v] = WS->mypos[v];
  struct wind_cell *wc = 0;
  wc                   = mem.myalloc(wc, 1);

  // in 1D we try to account for geometry in setting density
  if (ndim == 1)
    wc->dist = grid->distance_vertex2cell(wpos, c);
  else
    wc->dist = grid->distance(wpos, cpos);

  //#ifndef NDEBUG
  if (wc->dist > WP->current_radius) {
    if (ndim > 1) {
      // only print warning messages if ndim>1. This happens with
      // cells near the origin in 1D where we use the distance to the
      // cell centre-of-volume, not the midpoing.
      spdlog::warn(
          "{}: Expected {} but got {}",
          "stellar_wind::add_cell() cell is outside radius", WP->current_radius,
          wc->dist);
      spdlog::info("wind pos {}", wpos);
      CI.print_cell(c);
    }
    return 1;
  }
  //#endif

  // Now set wc cell pointer to this one.  Also set c->isbd to indicate
  // that it is now boundary data (while also grid data).
  c.isbd     = true;
  c.isdomain = false;
  if (wc->dist < WP->Rstar)
    c.timestep = false;
  else
    c.timestep = true;
  wc->c = &c;

  wc->cfac = 1.0;  // correction factor for curvilinear coords

  // Calculate the polar angle theta
  // Set theta to 0 if 1D - no angle dependent wind in this case (should add
  // exit if angle + 1D)
  if (ndim == 1) {
    wc->theta = 0;
    // correction factor for distance should be set here
    // should be (3*(i+0.5)^2)/((i+1)^3-i^3)
    double i = (grid->distance(wpos, cpos) - 0.5 * dx) / dx;
    // wc->cfac = pow((pow(i+1.0,3) - pow(i,3)) / (3.0*(i+0.5)*(i+0.5)), 1.0);
    // This is a bit of a hack, correcting the wind density in cells near the
    // origin, obtained by trial and error.
    if (i <= 10.0) wc->cfac = 1.0 + 0.24 * exp(-pow(i, 0.7));
  }

  // Polar angle in 2D
  else if (ndim == 2) {
    // Opposite and adjacent of cell angle
    double opp = grid->difference_vertex2cell(wpos, c, Rcyl);
    double adj = grid->difference_vertex2cell(wpos, c, Zcyl);
    wc->theta  = atan(fabs(opp / adj));
  }

  // Polar angle in 3D
  else if (ndim == 3) {
    // Opposite and adjacent in X-Y plane
    double opp1 = grid->difference_vertex2cell(wpos, c, XX);
    double adj1 = grid->difference_vertex2cell(wpos, c, YY);
    // Opposite and adjacent in Y-Z plane
    double opp2 = grid->difference_vertex2cell(wpos, c, ZZ);
    double adj2 = sqrt(opp1 * opp1 + adj1 * adj1);
    wc->theta   = atan(fabs(adj2 / opp2));
    // DEBUG
    if (!isfinite(wc->theta)) {
      // spdlog::debug(
      //    "inf theta={}, opp1={}, adj1={}, opp2={}, adj2={}, arg={}",
      //    wc->theta, opp1, adj1, opp2, adj2, fabs(adj2 / opp2));
      spdlog::error("{}: {}", "theta is not finite.", wc->theta);
      exit(1);
    }
    // DEBUG
  }

  // Allocate memory for wind_cell reference state vector.
  wc->p = 0;
  wc->p = mem.myalloc(wc->p, nvar);

  // Now assign values to the state vector:
  // NB: This function has the EOS Gamma hardcoded to 5/3
  set_wind_cell_reference_state(grid, wc, WS, 5. / 3.);

  WS->wcells.push_back(wc);
  WS->ncell += 1;

#ifndef NDEBUG
  // spdlog::debug("*** dist={}", wc->dist);
  // spdlog::debug("Wind BC cell pos : {}", wc->c->pos);
  // spdlog::debug(
  //    "Wind BC cell values : {}", std::vector<double>(wc->p, wc->p + nvar));
  // CI.print_cell(c);
  // spdlog::debug(
  //    "Added cell: array size={}\tncell={}", WS->wcells.size(), WS->ncell);
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind::remove_cells(const int id  ///< src id
)
{
  if (id < 0 || id >= nsrc) spdlog::error("{}: {}", "bad src id", id);
  struct wind_source *WS = wlist[id];
  // Set former wind-cells to normal domain cells
  // Clear list of Windcells
  for (int i = 0; i < WS->ncell; i++) {
    WS->wcells[i]->p = mem.myfree(WS->wcells[i]->p);
    if (WS->wcells[i]->c->isedge >= 0) WS->wcells[i]->c->isbd = false;
    WS->wcells[i]->c->isdomain = true;
    WS->wcells[i]->c->timestep = true;
    WS->wcells[i]              = mem.myfree(WS->wcells[i]);
  }
  // WS->wcells.clear();
  WS->wcells.erase(WS->wcells.begin(), WS->wcells.begin() + WS->ncell);
  WS->ncell = WS->wcells.size();
  // Set counter for windcells to zero
  return 0;
}


// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_reference_state(
    class GridBaseClass *grid,
    struct wind_cell *wc,
    const struct wind_source *WS,
    const double gamma  ///< EOS gamma
)
{
  if (WS->pars->acc) {
    set_wind_cell_reference_state_acc(*grid, *wc, *WS, gamma);
  }
  else {
    set_wind_cell_reference_state_vinf(*grid, *wc, *WS, gamma);
  }
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_reference_state_vinf(
    class GridBaseClass &grid,
    struct wind_cell &wc,
    const struct wind_source &WS,
    const double gamma  ///< EOS gamma
)
{
  struct stellarwind_params *WP = WS.pars;

  // if inside the star, set values to very small values and return
  if (wc.dist < WP->Rstar) {
    set_stellar_interior_values(wc, *WP);
    return;
  }

  // In this function we set the density, pressure, velocity, and tracer
  // values for the reference state of the cell.  Every timestep the
  // cell-values will be reset to this reference state.
#ifndef NDEBUG
  cell *c = wc.c;
#endif
  double vcell = WP->Vinf;
  set_wind_cell_density_pressure(grid, wc, *WP, vcell, gamma);

  // set velocity and magnetic field components
  std::array<double, MAX_DIM> r;
  set_wind_cell_offset(grid, wc, WS, r);
  set_wind_cell_velocity_components_vinf(r, wc, *WP, vcell, WP->Vrot);
  if (eqntype == EQMHD || eqntype == EQGLM) {
    set_wind_cell_B_components_vinf(r, wc, *WP, vcell, WP->Vrot);
  }

#ifndef NDEBUG
  // spdlog::debug(
  //    "Set Wind Ref. State: id {}, dist {:9.3e}, rho {:9.3e}, vx {:9.3e} vcell
  //    {:9.3e}, pg {:9.3e}", c->id, wc.dist, wc.p[RO], wc.p[VX], vcell,
  //    wc.p[PG]);
#endif

  set_wind_cell_tracers(wc, WS);

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  // Set a minimum temperature in the wind
  if (MP) {
    if (MP->Temperature(wc.p, gamma) < Tmin) {
      MP->Set_Temp(wc.p, Tmin, gamma);
    }
  }
  else {
    // appropriate for a neutral gas, He+M mass fraction 0.285.
    wc.p[PG] =
        max(static_cast<double>(wc.p[PG]),
            Tmin * wc.p[RO] * pconst.kB() * 0.78625 / pconst.m_p());
  }
#endif  // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return;
}


// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_reference_state_acc(
    class GridBaseClass &grid,
    struct wind_cell &wc,
    const struct wind_source &WS,
    const double gamma  ///< EOS gamma
)
{
  struct stellarwind_params *WP = WS.pars;

  // if inside the star, set values to very small values and return
  if (wc.dist < WP->Rstar) {
    set_stellar_interior_values(wc, *WP);
    return;
  }

  // hardcode beta=1 for now.  All the following equations assume beta==1

  // hardcode injection velocity to 4x the surface isothermal sound speed
  double v0 = 4.0 * sqrt(pconst.kB() * WP->Tstar / pconst.m_p());

  // radial component of velocity follows beta law
  double vcell = v0 + (WP->Vinf - v0) * (1.0 - WP->Rstar / wc.dist);

  // initialise Alfven radius to stellar radius, and Alfven velocity to
  // small fraction of sound speed.
  double ra = WP->Rstar;
  double va = 0.01 * v0;
  // if MHD, then set ra, va appropriately
  if (eqntype == EQMHD || eqntype == EQGLM) {
    // Alfven radius and velocity can be obtained analytically for beta=1 and
    // assuming v0 << vinf (where v0 is the injection velocity = 4c(Rstar) )
    ra = 0.5 * WP->Rstar
         * (1.0
            + sqrt(
                1.0
                + 4.0 * pow(WP->Bstar * WP->Rstar, 2) / (WP->Mdot * WP->Vinf)));
    // va = Alfven speed at ra == wind speed at ra by definition.
    va = v0 + (WP->Vinf - v0) * (1.0 - WP->Rstar / ra);
  }

  // set gas density
  set_wind_cell_density_pressure(grid, wc, *WP, vcell, gamma);

  // set velocity components
  std::array<double, MAX_DIM> r;
  set_wind_cell_offset(grid, wc, WS, r);
  double vphi =
      set_wind_cell_velocity_components_acc(r, wc, *WP, vcell, va, ra);
  // set B-field components
  if (eqntype == EQMHD || eqntype == EQGLM) {
    set_wind_cell_B_components_acc(r, wc, *WP, vcell, vphi, va, ra);
  }

  // set tracer values
  set_wind_cell_tracers(wc, WS);

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  // Set a minimum temperature in the wind
  if (MP) {
    if (MP->Temperature(wc.p, gamma) < Tmin) {
      MP->Set_Temp(wc.p, Tmin, gamma);
    }
  }
  else {
    // appropriate for a neutral gas, He+M mass fraction 0.285.
    wc.p[PG] =
        max(static_cast<double>(wc.p[PG]),
            Tmin * wc.p[RO] * pconst.kB() * 0.78625 / pconst.m_p());
  }
#endif  // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_stellar_interior_values(
    struct wind_cell &wc,                ///< cell to calculate for
    const struct stellarwind_params &wp  ///< wind source struct
)
{
  // if inside the star, set values to large values and return
  wc.p[RO] = 1.0e-9;
  wc.p[PG] = 1.0e4;
  wc.p[VX] = 0.0;
  wc.p[VY] = 0.0;
  wc.p[VZ] = 0.0;
  if (eqntype == EQMHD || eqntype == EQGLM) {
    wc.p[BX] = 0.0;
    wc.p[BY] = 0.0;
    wc.p[BZ] = wp.Bstar;
    if (eqntype == EQGLM) {
      wc.p[SI] = 0.0;
    }
  }
  // update tracers: should be set already.
  for (int v = 0; v < ntracer; v++)
    wc.p[ftr + v] = wp.tr[v];

  // include stellar space velocity if appropriate
  if (wp.moving_star) {
    for (int v = 0; v < ndim; v++)
      wc.p[VX + v] = wp.velocity[v];
  }
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_density_pressure(
    class GridBaseClass &,
    struct wind_cell &wc,
    const struct stellarwind_params &wp,
    const double vcell,  ///< radial velocity of wind at cell
    const double gamma   ///< EOS gamma
)
{
  // Density at cell position: rho = Mdot/(4.pi.R^2.v_inf) (for 3D)
  // or in 2D (slab-symmetry) rho = Mdot/(2.pi.R.v_inf)

  if (ndim == 2 && coordsys == COORD_CRT) {
    // 2D slab symmetry --> 1/r force laws and density profile
    wc.p[RO] = wp.Mdot / (vcell * 2.0 * M_PI * wc.dist);
    // Set pressure based on wind density/temperature at the stellar radius:
    // rho_star = Mdot/(2.pi.R_star.v_inf), p_star =
    // rho_star.k.T_star/(mu.m_p) and then p(r) = p_star
    // (rho(r)/rho_star)^gamma
    wc.p[PG] = pconst.kB() * wp.Tstar / pconst.m_p();
    wc.p[PG] *=
        exp((gamma - 1.0) * log(2.0 * M_PI * wp.Rstar * vcell / wp.Mdot));
    wc.p[PG] *= exp((gamma)*log(wc.p[RO]));
  }

  else {
    // 3D geometry, so either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
    // rho = Mdot/(4.pi.R^2.v_inf)
    wc.p[RO] = 1.0 / (wc.dist * wc.cfac);
    wc.p[RO] *= wc.p[RO];
    wc.p[RO] *= wp.Mdot / (vcell * 4.0 * M_PI);
    //
    // Set pressure based on wind density/temperature at the stellar
    // radius, assuming adiabatic expansion outside Rstar
    wc.p[PG] = pconst.kB() * wp.Tstar / pconst.m_p();
    wc.p[PG] *=
        exp((gamma - 1.0)
            * log(4.0 * M_PI * wp.Rstar * wp.Rstar * vcell / wp.Mdot));
    wc.p[PG] *= exp((gamma)*log(wc.p[RO]));
  }
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_tracers(
    struct wind_cell &wc,         ///< cell to calculate for
    const struct wind_source &ws  ///< wind source struct
)
{
  // set H+ tracer value based on stellar temperature.
  if (ws.Hplus >= 0) {
    if (ws.pars->Tstar < 1.0e4)
      ws.pars->tr[ws.iHplus] = 1.0e-10;
    else if (ws.pars->Tstar > 1.5e4)
      ws.pars->tr[ws.iHplus] = 1.0;
    else
      ws.pars->tr[ws.iHplus] =
          1.0e-10 + (ws.pars->Tstar - 1.0e4) * (1.0 - 1.0e-10) / 0.5e4;
  }
  // update tracers: should be set already.
  for (int v = 0; v < ntracer; v++)
    wc.p[ftr + v] = ws.pars->tr[v];

  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_offset(
    class GridBaseClass &grid,
    struct wind_cell &wc,
    const struct wind_source &ws,
    std::array<double, MAX_DIM> &r)
{
  // Velocities and magnetic fields: get coordinates relative to star
  // for calculating sin/cos angles in theta and phi.
  cell *c = wc.c;

  switch (ndim) {
    case 1:
      // in 1D, v_r = v_infty, so need r[XX] = wc->dist.
      r[XX] = wc.dist;
      r[YY] = 0.0;
      r[ZZ] = 0.0;
      break;
    case 2:
      r[XX] = grid.difference_vertex2cell(ws.mypos, *c, XX);
      r[YY] = grid.difference_vertex2cell(ws.mypos, *c, YY);
      r[ZZ] = 0.0;
      break;
    case 3:
      r[XX] = grid.difference_vertex2cell(ws.mypos, *c, XX);
      r[YY] = grid.difference_vertex2cell(ws.mypos, *c, YY);
      r[ZZ] = grid.difference_vertex2cell(ws.mypos, *c, ZZ);
      break;
    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_velocity_components_vinf(
    const std::array<double, MAX_DIM> &r,
    struct wind_cell &wc,
    const struct stellarwind_params &wp,
    const double vcell,  ///< radial velocity of wind at cell
    const double vrot    ///< rotation velocity at equator
)
{
  // Velocities: cell-average values, i.e. values at the
  // centre-of-volume.
  // TODO: for general J vector in 3D, what is rotational component.
  switch (ndim) {
    case 1:
      wc.p[VX] = vcell * r[XX] / wc.dist;
      wc.p[VY] = 0.0;
      wc.p[VZ] = 0.0;
      break;

    case 2:
      if (coordsys == COORD_CRT) {
        wc.p[VX] = vcell * r[XX] / wc.dist;
        wc.p[VY] = vcell * r[YY] / wc.dist;
        wc.p[VX] += -vrot * r[YY] / wc.dist;
        wc.p[VY] += vrot * r[XX] / wc.dist;
        wc.p[VZ] = 0.0;
      }
      else {
        wc.p[VX] = vcell * r[XX] / wc.dist;
        wc.p[VY] = vcell * r[YY] / wc.dist;
        // J is hardcoded to be parallel to positive z-axis
        wc.p[VZ] = vrot * wp.Rstar * r[YY] / pconst.pow_fast(wc.dist, 2);
      }
      break;

    case 3:
      wc.p[VX] = vcell * r[XX] / wc.dist;
      wc.p[VY] = vcell * r[YY] / wc.dist;
      wc.p[VZ] = vcell * r[ZZ] / wc.dist;

      // add non-radial component to x/y-dir from rotation.
      // J is hardcoded to be parallel to positive z-axis
      wc.p[VX] += -vrot * wp.Rstar * r[YY] / pconst.pow_fast(wc.dist, 2);
      wc.p[VY] += vrot * wp.Rstar * r[XX] / pconst.pow_fast(wc.dist, 2);
      break;

    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }

  // include stellar space velocity if appropriate
  if (wp.moving_star) {
    for (int v = 0; v < ndim; v++)
      wc.p[VX + v] += wp.velocity[v];
  }

  return;
}



// ##################################################################
// ##################################################################



double stellar_wind::set_wind_cell_velocity_components_acc(
    const std::array<double, MAX_DIM> &r,
    struct wind_cell &wc,
    const struct stellarwind_params &wp,
    const double vcell,  ///< radial velocity of wind at cell
    const double va,     ///< Alfven velocity at Alfven radius
    const double ra      ///< Alfven radius
)
{
  // Velocities: cell-average values, i.e. values at the
  // centre-of-volume.
  // TODO: for general J vector in 3D, what is rotational component.

  // This is eq. 9.37 in Lamers & Cassinelli (1997), x and u defined in
  // eq. 9.35
  double u    = vcell / va;
  double x2   = pconst.pow_fast(wc.dist / ra, 2);
  double vphi = wp.Vrot * wc.dist / wp.Rstar * (1.0 - u) / (1.0 - x2 * u);

  switch (ndim) {
    case 1:
      // no rotation in 1D
      wc.p[VX] = vcell * r[XX] / wc.dist;
      wc.p[VY] = 0.0;
      wc.p[VZ] = 0.0;
      break;

    case 2:
      if (coordsys == COORD_CRT) {
        wc.p[VX] = vcell * r[XX] / wc.dist;
        wc.p[VY] = vcell * r[YY] / wc.dist;
        wc.p[VX] += -vphi * r[YY] / wc.dist;
        wc.p[VY] += vphi * r[XX] / wc.dist;
        wc.p[VZ] = 0.0;
      }
      else {
        wc.p[VX] = vcell * r[XX] / wc.dist;
        wc.p[VY] = vcell * r[YY] / wc.dist;
        // J is hardcoded to be parallel to positive z-axis
        wc.p[VZ] = vphi * r[YY] / wc.dist;
      }
      break;

    case 3:
      wc.p[VX] = vcell * r[XX] / wc.dist;
      wc.p[VY] = vcell * r[YY] / wc.dist;
      wc.p[VZ] = vcell * r[ZZ] / wc.dist;

      // add non-radial component to x/y-dir from rotation.
      // J is hardcoded to be parallel to positive z-axis
      wc.p[VX] += -vphi * r[YY] / wc.dist;
      wc.p[VY] += vphi * r[XX] / wc.dist;
      break;

    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }

  // include stellar space velocity if appropriate
  if (wp.moving_star) {
    for (int v = 0; v < ndim; v++)
      wc.p[VX + v] += wp.velocity[v];
  }

  return vphi;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_B_components_vinf(
    const std::array<double, MAX_DIM> &r,
    struct wind_cell &wc,
    const struct stellarwind_params &wp,
    const double vcell,  ///< radial velocity of wind at cell
    const double vrot    ///< rotation velocity at equator
)
{
  // Magnetic field: cell-average values, i.e. values at the
  // centre-of-volume.
  // Use a split monopole plus a rotational term adding toroidal
  // component.
  // TODO: for general J vector, what is rotational component.
  double B_s = wp.Bstar / sqrt(4.0 * M_PI);  // code units for B_surf
  double D_s = wp.Rstar / wc.dist;           // 1/d in stellar radii
  double D_2 = D_s * D_s;                    // 1/d^2 in stellar radii
  // this multiplies the toroidal component:
  double beta_B_sint = (vrot / vcell) * B_s * D_s;

  switch (ndim) {
    case 1:
      spdlog::error("{}: {}", "1D spherical but MHD?", ndim);
      break;
    case 2:
      if (coordsys == COORD_CYL) {
        // split monopole
        wc.p[BX] = B_s * D_2 * fabs(r[XX]) / wc.dist;
        wc.p[BY] = B_s * D_2 / wc.dist;
        wc.p[BY] = (r[XX] > 0.0) ? r[YY] * wc.p[BY] : -r[YY] * wc.p[BY];
        // toroidal component
        beta_B_sint = beta_B_sint * r[YY] / wc.dist;
        wc.p[BZ]    = (r[XX] > 0.0) ? -beta_B_sint : beta_B_sint;
      }
      else {
        // Cartesian: take vertical field in z direction
        beta_B_sint = wp.Bstar / sqrt(4.0 * M_PI) * D_s;
        // wc.p[BX]   = -beta_B_sint * r[YY] / wc.dist;
        // wc.p[BY]   = beta_B_sint * r[XX] / wc.dist;
        wc.p[BX] = 0.0;
        wc.p[BY] = 0.0;
        wc.p[BZ] = beta_B_sint;
      }

      break;

    case 3:
      // split monopole along z-axis, parallel to J
      wc.p[BX] = B_s * D_2 / wc.dist;
      wc.p[BX] = (r[ZZ] > 0.0) ? r[XX] * wc.p[BX] : -r[XX] * wc.p[BX];

      wc.p[BY] = B_s * D_2 / wc.dist;
      wc.p[BY] = (r[ZZ] > 0.0) ? r[YY] * wc.p[BY] : -r[YY] * wc.p[BY];

      wc.p[BZ] = B_s * D_2 * fabs(r[ZZ]) / wc.dist;

      // toroidal component in r[XX]-y plane from rotation, such that
      // we have a Parker spiral, inward winding for z<0 and
      // outwards for z>0.
      beta_B_sint *= sqrt(r[XX] * r[XX] + r[YY] * r[YY]) / wc.dist;
      beta_B_sint = (r[ZZ] > 0.0) ? -beta_B_sint : beta_B_sint;
      wc.p[BX] += -beta_B_sint * r[YY] / wc.dist;
      wc.p[BY] += beta_B_sint * r[XX] / wc.dist;
      break;

    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }

  if (eqntype == EQGLM) {
    wc.p[SI] = 0.0;
  }
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind::set_wind_cell_B_components_acc(
    const std::array<double, MAX_DIM> &r,
    struct wind_cell &wc,
    const struct stellarwind_params &wp,
    const double vcell,  ///< radial velocity of wind at cell
    const double vphi,   ///< phi component of wind velocity at cell
    const double va,     ///< Alfven velocity at Alfven radius
    const double ra      ///< Alfven radius
)
{
  // Magnetic field: cell-average values, i.e. values at the
  // centre-of-volume.
  // Use a split monopole plus a rotational term adding toroidal
  // component.
  // TODO: for general J vector, what is rotational component.
  double D_s = wp.Rstar / wc.dist;  // 1/d in stellar radii
  // convert to code units (divide by sqrt(4pi))  L&C eq. 9.11
  double Br = wp.Bstar / sqrt(4.0 * M_PI) * D_s * D_s;
  // From Lamers & Cassinelli (eq. 9.12)
  // Note Bphi here is always negative for r>Rstar, opposite sign to the term
  // "beta_B_sint" in the vinf case, so the sign is reversed in the below
  // expresssions compared with the _vinf() function.
  // double Bphi = Br * (vphi - wp.Vrot / D_s) / vcell;
  // L&C eq. 9.38
  double x    = wc.dist / ra;
  double u    = vcell / va;
  double Bphi = -Br * vphi * (1.0 - x * x) / (va * (1.0 - u));

  switch (ndim) {
    case 1:
      spdlog::error("{}: {}", "1D spherical but MHD?", ndim);
      break;
    case 2:
      if (coordsys == COORD_CYL) {
        // split monopole
        wc.p[BX] = Br * fabs(r[XX]) / wc.dist;
        wc.p[BY] = Br / wc.dist;
        wc.p[BY] = (r[XX] > 0.0) ? r[YY] * wc.p[BY] : -r[YY] * wc.p[BY];
        // toroidal component
        Bphi     = Bphi * r[YY] / wc.dist;
        wc.p[BZ] = (r[XX] > 0.0) ? Bphi : -Bphi;
      }
      else {
        // Cartesian: take vertical field in z direction.
        Bphi = wp.Bstar / sqrt(4.0 * M_PI) * D_s;
        // wc.p[BX]   = -Bphi * r[YY] / wc.dist;
        // wc.p[BY]   = Bphi * r[XX] / wc.dist;
        wc.p[BX] = 0.0;
        wc.p[BY] = 0.0;
        wc.p[BZ] = Bphi;
      }

      break;

    case 3:
      // split monopole along z-axis, parallel to J
      wc.p[BX] = Br / wc.dist;
      wc.p[BX] = (r[ZZ] > 0.0) ? r[XX] * wc.p[BX] : -r[XX] * wc.p[BX];

      wc.p[BY] = Br / wc.dist;
      wc.p[BY] = (r[ZZ] > 0.0) ? r[YY] * wc.p[BY] : -r[YY] * wc.p[BY];

      wc.p[BZ] = Br * fabs(r[ZZ]) / wc.dist;

      // toroidal component in r[XX]-y plane from rotation, such that
      // we have a Parker spiral, inward winding for z<0 and
      // outwards for z>0.
      Bphi *= sqrt(r[XX] * r[XX] + r[YY] * r[YY]) / wc.dist;
      Bphi = (r[ZZ] > 0.0) ? Bphi : -Bphi;
      wc.p[BX] += -Bphi * r[YY] / wc.dist;
      wc.p[BY] += Bphi * r[XX] / wc.dist;
      break;

    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }

  if (eqntype == EQGLM) {
    wc.p[SI] = 0.0;
  }
  return;
}


// ##################################################################
// ##################################################################



int stellar_wind::get_num_cells(const int id  ///< src id
)
{
  if (id < 0 || id >= nsrc) spdlog::error("{}: {}", "bad src id", id);
  return wlist[id]->ncell;
}



// ##################################################################
// ##################################################################



int stellar_wind::set_cell_values(
    class GridBaseClass *,
    const int id,  ///< src id
    const double   ///< simulation time
)
{
  if (id < 0 || id >= nsrc) spdlog::error("{}: {}", "bad src id", id);

  //
  // go through every cell in one go and update them all
  //
  struct wind_source *WS = wlist[id];
  spdlog::debug("updating source {} which has {} cells", id, WS->ncell);
  for (int i = 0; i < WS->ncell; i++) {
    for (int v = 0; v < nvar; v++)
      WS->wcells[i]->c->P[v] = WS->wcells[i]->p[v];
    for (int v = 0; v < nvar; v++)
      WS->wcells[i]->c->Ph[v] = WS->wcells[i]->p[v];
  }
  spdlog::debug("updated source {} which has {} cells", id, WS->ncell);
  return 0;
}



// ##################################################################
// ##################################################################



void stellar_wind::get_src_posn(
    const int id,                   ///< src id
    std::array<double, MAX_DIM> &x  ///< position vector (output)
)
{
  for (int v = 0; v < ndim; v++)
    x[v] = wlist[id]->pars->dpos[v];
}



// ##################################################################
// ##################################################################



#ifdef ANALYTIC_ORBITS
void stellar_wind::get_src_orbit(
    const int id,  ///< src id
    int *moving,   ///< is star moving? 1=yes, 0=no (output)
    double *x,     ///< eccentricity (output)
    double *y1,    ///< PeriastronX Vec (output)
    double *y2,    ///< PeriastronY Vec (output)
    double *z,     ///< Orbital period (output)
    double *w      ///< initial position
)
{
  *x  = wlist[id]->pars->eccentricity;
  *y1 = wlist[id]->pars->PeriastronX;
  *y2 = wlist[id]->pars->PeriastronY;
  *z  = wlist[id]->pars->OrbPeriod;
  for (int v = 0; v < ndim; v++)
    w[v] = wlist[id]->dpos_init[v];
  *moving = wlist[id]->pars->moving_star;
}
#endif  // ANALYTIC_ORBITS


// ##################################################################
// ##################################################################



void stellar_wind::set_src_posn(
    const int id,  ///< src id
    double *x      ///< position vector (output)
)
{
  for (int v = 0; v < ndim; v++)
    wlist[id]->pars->dpos[v] = x[v];
}

//------------------------------------------------

// ##################################################################
// ##################################################################

double stellar_wind::beta(const double Teff)
{
  //
  // Eldridge et al. (2006, MN, 367, 186).
  // Their eq. 2 has a typo, because the hot star values from Vink+
  // (2001) are for v_inf/v_esc, not for (v_inf/v_esc)^2.
  // So the beta value scales the escape velocity itself, not the
  // square.
  //
  double beta;
  double rsg = 0.125;  // Eldridge value
  // double rsg=0.04;  // Mackey+2012 Betelgeuse value

  if (Teff <= 3600.0)
    beta = rsg;
  else if (Teff >= 22000.0)
    beta = 2.6;
  else {
    //
    // Linear interpolation for beta from Eldridge et al. Table 1.
    //
    double b0, b1, T0, T1;
    if (Teff < 6000.0) {
      T0 = 3600.0;
      b0 = rsg;
      T1 = 6000.0;
      b1 = 0.5;
    }
    else if (Teff < 8000.0) {
      T0 = 6000.0;
      b0 = 0.5;
      T1 = 8000.0;
      b1 = 0.7;
    }
    else if (Teff < 10000.0) {
      T0 = 8000.0;
      b0 = 0.7;
      T1 = 10000.0;
      b1 = 1.3;
    }
    else if (Teff < 20000.0) {
      T0 = 10000.0;
      b0 = 1.3;
      T1 = 20000.0;
      b1 = 1.3;
    }
    else {
      T0 = 20000.0;
      b0 = 1.3;
      T1 = 22000.0;
      b1 = 2.6;
    }
    beta = b0 + (Teff - T0) * (b1 - b0) / (T1 - T0);
  }

  return beta;
}



// ##################################################################
// ##################################################################



double stellar_wind::get_min_wind_radius(
    const int id,    ///< src id
    const double dx  ///< cell size
)
{
  // make sure boundary radius is larger than Alfven radius to avoid numerical
  // problems.
  // For this we can use Lamers & Cassinelli, ch. 9 equations, assuming
  // beta-law wind acceleration with beta = 1
  struct stellarwind_params *p = wlist[id]->pars;
  double ra                    = p->Rstar;
  if ((eqntype == EQMHD || eqntype == EQGLM) && p->acc == 1) {
    ra = 0.5 * p->Rstar
         * (1.0
            + sqrt(
                1.0 + 4.0 * pow(p->Bstar * p->Rstar, 2) / (p->Mdot * p->Vinf)));
    // spdlog::info("pars: {:12.6e}, {:12.6e}, {:12.6e}, {:12.6e}, {:12.6e}",
    //            ra, p->Bstar, p->Rstar, p->Mdot, p->Vinf);
  }

  // wind boundary should be 4 cells larger than R_Alfven, and at least 2 cells
  // larger than R_star
  ra = max(ra + 2.0 * dx, wlist[id]->pars->Rstar + 2.0 * dx);
  // wind boundary should be at least MIN_WIND_RAD cells in radius.
  return max(ra, MIN_WIND_RAD * dx);
}



// ##################################################################
// ##################################################################

// ------------------------------------------------------------------
// ----------  STELLAR WIND WITH STELLAR EVOLUTION ------------------
// ------------------------------------------------------------------

stellar_wind_evolution::stellar_wind_evolution(
    const int nd,           ///< ndim
    const int nv,           ///< nvar
    const int nt,           ///< ntracer
    const int ft,           ///< ftr
    const std::string *tr,  ///< List of tracer variable names.
    const int cs,           ///< coord_sys
    const int eq,           ///< eqn_type
    const double mt,        ///< Minimum temperature allowed on grid
    const double ss,        ///< Simulation start time.
    const double sf         ///< Simulation finish time.
    ) :
    stellar_wind(nd, nv, nt, ft, tr, cs, eq, mt),
    sim_start(ss), sim_finish(sf)
{
  spdlog::info("Stellar wind with time evolution, constructor");
}

// ##################################################################
// ##################################################################

stellar_wind_evolution::~stellar_wind_evolution()
{
  spdlog::info("Stellar wind with time evolution, destructor");

  //
  // Delete arrays
  //
  // delete each wdata_evol[] element, making sure to delete the
  // arrays in each element first (arrays, pos, tracers).
  //
  struct evolving_wind_data *wd = 0;
  while (wdata_evol.size() > 0) {
    wd = wdata_evol.back();
    wdata_evol.pop_back();
    wd = mem.myfree(wd);
  }
}

// ##################################################################
// ##################################################################

int stellar_wind_evolution::add_source(
    struct stellarwind_params *pars  ///< pointer to wind parameters struct
)
{
  //
  // This function sets up an evolving_wind_data struct to hold the
  // constant wind.  It sets all the evolving wind stuff to zero, and
  // then calls the stellar_wind:: version to set up the wind source.
  //
  struct evolving_wind_data *temp = 0;
  temp                            = mem.myalloc(temp, 1);
  //
  // First set all the stuff we don't need for a constant wind.
  //
  temp->Npt           = 0;
  temp->offset        = 0.0;
  temp->tstart        = sim_start - 1.0e10;  // so it has already started.
  temp->tfinish       = 2.0 * sim_finish;    // so it keeps going to end of sim.
  temp->update_freq   = 1.0e99;              // so it never updates.
  temp->t_next_update = 1.0e99;              // so it never updates.

  //
  // Set source to be active
  //
  temp->is_active = true;

  // find indices of elements and dust in tracers list.
  set_element_indices(temp);
  // Now add source using constant wind version.
  stellar_wind::add_source(pars);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);

  return 0;
}



// ##################################################################
// ##################################################################



void stellar_wind_evolution::set_element_indices(struct evolving_wind_data *ewd)
{
  ewd->i_XH  = -1;
  ewd->i_XHe = -1;
  ewd->i_XC  = -1;
  ewd->i_XN  = -1;
  ewd->i_XO  = -1;
  ewd->i_XZ  = -1;
  ewd->i_XD  = -1;

  // Check for elements if using microphysics
  if (MP) {
    for (int v = 0; v < ntracer; v++) {
      if (tracers[v] == "X_H")
        ewd->i_XH = v;
      else if (tracers[v] == "X_He")
        ewd->i_XHe = v;
      else if (tracers[v] == "X_C")
        ewd->i_XC = v;
      else if (tracers[v] == "X_N")
        ewd->i_XN = v;
      else if (tracers[v] == "X_O")
        ewd->i_XO = v;
      else if (tracers[v] == "X_Z")
        ewd->i_XZ = v;
      else if (tracers[v] == "X_D")
        ewd->i_XD = v;
    }
  }
  return;
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::read_evolution_file(
    const string infile,             ///< file name to read data from.
    struct evolving_wind_data *data  ///< where to put the data
)
{

  //
  // Read in stellar evolution data
  // Format: time M L Teff Mdot vrot vcrit vinf
  //     X_H X_He X_C X_N X_O X_Z X_D
  //
  FILE *wf = 0;
  wf       = fopen(infile.c_str(), "r");
  if (!wf)
    spdlog::error("can't open wind file, stellar_wind_evo", fmt::ptr(wf));

  // Skip first two lines
  char line[512];
  char *rval = 0;
  rval       = fgets(line, 512, wf);
  if (!rval)
    spdlog::error("{}: {}", "stwind_angle: failed to get line 1", line);
  rval = fgets(line, 512, wf);
  if (!rval)
    spdlog::error("{}: {}", "stwind_angle: failed to get line 2", line);

  // read file line by line and add to struct vectors.
  // Everthing must be in CGS units already.
  // format:
  // time  mass  luminosity  T_eff  Mdot  v_rot  v_crit  v_inf
  //  X_H  X_He  X_C  X_N  X_O  X_Z  X_D
  double time = 0.0, mass = 0.0, lumi = 0.0, teff = 0.0, radi = 0.0, mdot = 0.0,
         vrot = 0.0, vcrt = 0.0, vinf = 0.0;
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;
  while ((rval = fgets(line, 512, wf)) != 0) {
    sscanf(
        line,
        "   %lE   %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
        &time, &mass, &lumi, &teff, &mdot, &vrot, &vcrt, &vinf, &xh, &xhe, &xc,
        &xn, &xo, &xz, &xd);

    // spdlog::debug(
    //    "{} {} {} {} {} {} {} {}", time, mass, lumi, teff, mdot, vrot, vcrt,
    //    vinf);

    // Set vector value
    data->time_evo.push_back(time);
    data->M_evo.push_back(mass);
    data->L_evo.push_back(lumi);
    data->Teff_evo.push_back(teff);

    // Stellar radius
    radi = sqrt(
        lumi
        / (4.0 * pconst.pi() * pconst.StefanBoltzmannConst()
           * pconst.pow_fast(teff, 4.0)));
    data->R_evo.push_back(radi);

    data->Mdot_evo.push_back(mdot);
    data->vrot_evo.push_back(vrot);
    data->vcrt_evo.push_back(vcrt);
    data->vinf_evo.push_back(vinf);
    data->X_H_evo.push_back(xh);
    data->X_He_evo.push_back(xhe);
    data->X_C_evo.push_back(xc);
    data->X_N_evo.push_back(xn);
    data->X_O_evo.push_back(xo);
    data->X_Z_evo.push_back(xz);
    data->X_D_evo.push_back(xd);
  }
  fclose(wf);

  // Column length
  size_t Npt = data->time_evo.size();
  data->Npt  = Npt;

  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::add_evolving_source(
    const double t_now,              ///< current time.
    struct stellarwind_params *pars  ///< pointer to wind parameters struct
)
{
  if (pars->type != WINDTYPE_EVOLVING) {
    spdlog::error(
        "{}: {}", "Bad wind type for evolving stellar wind!", pars->type);
  }
  //
  // First we will read the file, and see when the source should
  // switch on in the simulation (it may not be needed for a while).
  //
  spdlog::debug(
      "\t\tsw-evo: adding source from file {}", pars->evolving_wind_file);

  //
  // Wind source struct, to be added to class vector wdata_evol
  //
  struct evolving_wind_data *temp = 0;
  temp                            = mem.myalloc(temp, 1);
  int err = read_evolution_file(pars->evolving_wind_file, temp);
  if (err)
    spdlog::error(
        "{}: {}", "couldn't read wind evolution file",
        pars->evolving_wind_file);

  //
  // Optional time offset between simulation time and evolutionary
  // time.  Also optional scaling.
  //
  for (int i = 0; i < temp->Npt; i++) {
    temp->time_evo[i] += pars->time_offset;
    temp->time_evo[i] /= pars->t_scalefactor;
  }

  //
  // Offset is not used in the code past here.  All times are in
  // seconds.
  //
  temp->offset        = pars->time_offset / pars->t_scalefactor;
  temp->tstart        = temp->time_evo[0];
  temp->tfinish       = temp->time_evo[temp->Npt - 1];
  temp->update_freq   = pars->update_freq / pars->t_scalefactor;
  temp->t_next_update = max(temp->tstart, t_now);

  spdlog::debug(
      "\t\t tstart={}, next update={}, and tfinish={}", temp->tstart,
      temp->t_next_update, temp->tfinish);


  //
  // Decide if the wind src is active yet.  If it is, then
  // set up a constant wind source for updating its properties.
  //
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  if (((t_now + temp->update_freq) > temp->tstart
       || pconst.equalD(temp->tstart, t_now))
      && t_now < temp->tfinish) {
    temp->is_active = true;
    //
    // Get the current values for mdot, vinf, Teff, and setup a wind
    // source using the constant-wind function.
    //
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Teff_evo, t_now, pars->Tstar);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Mdot_evo, t_now, pars->Mdot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vinf_evo, t_now, pars->Vinf);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vrot_evo, t_now, pars->Vrot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->R_evo, t_now, pars->Rstar);

    // get tracer values for elements.
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_H_evo, t_now, xh);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->X_He_evo, t_now, xhe);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_C_evo, t_now, xc);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_N_evo, t_now, xn);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_O_evo, t_now, xo);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_Z_evo, t_now, xz);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_D_evo, t_now, xd);

    spdlog::info("Source is Active");
    spdlog::debug(
        "T = {}, mdot={}, vinf={}, rstar={}", pars->Tstar, pars->Mdot,
        pars->Vinf, pars->Rstar);
  }
  else {
    spdlog::warn(
        "Source is not yet active: tnow={}, tstart={}. Setting wind source to INACTIVE",
        t_now, temp->tstart);
    temp->is_active = false;
  }

  // set tracer values for elements
  set_element_indices(temp);
  if (temp->i_XH >= 0) pars->tr[temp->i_XH] = xh;
  if (temp->i_XHe >= 0) pars->tr[temp->i_XHe] = xhe;
  if (temp->i_XC >= 0) pars->tr[temp->i_XC] = xc;
  if (temp->i_XN >= 0) pars->tr[temp->i_XN] = xn;
  if (temp->i_XO >= 0) pars->tr[temp->i_XO] = xo;
  if (temp->i_XZ >= 0) pars->tr[temp->i_XZ] = xz;
  if (temp->i_XD >= 0) pars->tr[temp->i_XD] = xd;

  // Now add source using constant wind version.
  stellar_wind::add_source(pars);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_evolution::update_source(
    class GridBaseClass *grid,
    struct evolving_wind_data *wd,
    const double t_now,
    const double gamma)
{
  struct stellarwind_params *wp = wd->ws->pars;
  if (!wd->is_active) {
    array<double, MAX_DIM> wpos;
    for (int v = 0; v < MAX_DIM; v++)
      wpos[v] = wd->ws->mypos[v];
    spdlog::debug(
        "stellar_wind_evo::update_source activating source id={} at Simulation time t={}",
        wp->id, t_now);
    // spdlog::debug("Source position : {}", wd->ws->dpos);
    wd->is_active = true;
  }

  if (t_now < wd->tstart) {
    spdlog::warn("Updating source, not yet active! {} {}", wd->tstart, t_now);
    return 0;
  }
  else if (t_now >= wd->tfinish) {
    spdlog::warn(
        "Updating source: source no longer active! {} {}", wd->tfinish, t_now);
    return 1;
  }


  wd->t_next_update = t_now;  // (We update every timestep now)
  wd->t_next_update = min(wd->t_next_update, wd->tfinish);

  //
  // Now we update Mdot, Vinf, Teff by linear interpolation.
  //
  double mdot = 0.0, vinf = 0.0, vrot = 0.0, Twind = 0.0, rstar = 0.0;
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  interpolate.root_find_linear_vec(wd->time_evo, wd->Teff_evo, t_now, Twind);
  interpolate.root_find_linear_vec(wd->time_evo, wd->Mdot_evo, t_now, mdot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vrot_evo, t_now, vrot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vinf_evo, t_now, vinf);
  interpolate.root_find_linear_vec(wd->time_evo, wd->R_evo, t_now, rstar);

  wp->Mdot  = mdot;   // already cgs.
  wp->Vinf  = vinf;   // this is in cm/s already.
  wp->Vrot  = vrot;   // this is in cm/s already.
  wp->Tstar = Twind;  // This is in K.
  wp->Rstar = rstar;

  // get tracer values for elements.
  if (wd->i_XH >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_H_evo, t_now, xh);
    wp->tr[wd->i_XH] = xh;
  }
  if (wd->i_XHe >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_He_evo, t_now, xhe);
    wp->tr[wd->i_XHe] = xhe;
  }
  if (wd->i_XC >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_C_evo, t_now, xc);
    wp->tr[wd->i_XC] = xc;
  }
  if (wd->i_XN >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_N_evo, t_now, xn);
    wp->tr[wd->i_XN] = xn;
  }
  if (wd->i_XO >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_O_evo, t_now, xo);
    wp->tr[wd->i_XO] = xo;
  }
  if (wd->i_XZ >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_Z_evo, t_now, xz);
    wp->tr[wd->i_XZ] = xz;
  }
  if (wd->i_XD >= 0) {
    interpolate.root_find_linear_vec(wd->time_evo, wd->X_D_evo, t_now, xd);
    wp->tr[wd->i_XD] = xd;
  }

  //
  // Now re-assign state vector of each wind-boundary-cell with
  // updated values.
  //
  for (int i = 0; i < wd->ws->ncell; i++) {
    set_wind_cell_reference_state(grid, wd->ws->wcells[i], wd->ws, gamma);
  }

  return 0;
}

// ##################################################################
// ##################################################################

int stellar_wind_evolution::set_cell_values(
    class GridBaseClass *grid,
    const int id,       ///< src id
    const double t_now  ///< simulation time
)
{
  int err = 0;
  if (id < 0 || id >= nsrc) {
    spdlog::error("{}: {}", "bad src id", id);
  }

  struct evolving_wind_data *wd = wdata_evol[id];

  //
  // Check if the source needs to be updated, and if so, update it
  //
  if (t_now >= wd->t_next_update) {
    /// NB: This function has hardcoded EOS Gamma to 5/3
    err = update_source(grid, wd, t_now, 5. / 3.);
    if (err) {
      spdlog::info("update_source() star's life is over {}", err);
      return -1;
    }
  }

  //
  // Now we have an up-to-date source.  We check if it is active,
  // and if so we update the wind cell values using the constant wind
  // functions.  If not, then we ignore it and return.
  //
  if (wd->is_active) {
    err += stellar_wind::set_cell_values(grid, id, t_now);
  }

  return err;
}

// ##################################################################
// ##################################################################

double stellar_wind_evolution::Mdot_Vink_hot(
    const double L,  ///< luminosity (Lsun)
    const double M,  ///< mass (Msun)
    const double T,  ///< T_eff (K)
    const double Z,  ///< Metallicity wrt solar
    const double b   ///< beta of wind on hot side of BSJ
)
{
  double md = -6.697 + 2.194 * log10(L * 1e-5) - 1.313 * log10(M / 30.0)
              - 1.226 * log10(b / 2.0) + 0.933 * log10(T / 4.0e4)
              - 10.92 * pow(log10(T / 4.0e4), 2) + 0.85 * log10(Z);
  return exp(pconst.ln10() * md);
}

// ##################################################################
// ##################################################################

double stellar_wind_evolution::Mdot_Vink_cool(
    const double L,  ///< luminosity (Lsun)
    const double M,  ///< mass (Msun)
    const double T,  ///< T_eff (K)
    const double Z,  ///< Metallicity wrt solar
    const double b   ///< beta of wind on hot side of BSJ
)
{
  double md = -6.688 + 2.210 * log10(L * 1e-5) - 1.339 * log10(M / 30.0)
              - 1.601 * log10(b / 2.0) + 1.070 * log10(T / 2.0e4)
              + 0.850 * log10(Z);
  return exp(pconst.ln10() * md);
}

// ##################################################################
// ##################################################################

double stellar_wind_evolution::Mdot_Nieuwenhuijzen(
    const double L,  ///< luminosity (Lsun)
    const double M,  ///< mass (Msun)
    const double R,  ///< Radius (Rsun)
    const double Z   ///< Metallicity wrt solar
)
{
  double md = -14.02 + 1.24 * log10(L) + 0.16 * log10(M) + 0.81 * log10(R);
  return exp(pconst.ln10() * md) * Z;
}



// ##################################################################
// ##################################################################



double stellar_wind_evolution::Mdot_Brott(
    const double L,  ///< luminosity (Lsun)
    const double M,  ///< mass (Msun)
    const double T,  ///< T_eff (K)
    const double R,  ///< Radius (Rsun)
    const double Z   ///< Metallicity wrt solar
)
{
  double Tn    = 22.5e3;
  double Tp    = 27.1e3;
  double betaH = 2.6;  // wind velocity multiplier on hot side
  double md = 0.0, mdc = 0.0, mdh = 0.0;
  if (T > Tp) {
    md = Mdot_Vink_hot(L, M, T, Z, betaH);
  }
  else if (T > Tn) {
    mdh = Mdot_Vink_hot(L, M, T, Z, betaH);
    mdc = Mdot_Vink_cool(L, M, T, Z, betaH);
    md  = mdc + (T - Tn) * (mdh - mdc) / (Tp - Tn);
  }
  else {
    mdh = Mdot_Vink_cool(L, M, T, Z, betaH);
    mdc = Mdot_Nieuwenhuijzen(L, M, R, Z);
    md  = std::max(mdh, mdc);
  }
  return md * pconst.Msun() / pconst.year();
}



// ##################################################################
// ##################################################################
