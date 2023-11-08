///
/// \file solver_eqn_hydro_adi.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 written, based on old solver.cc classes.
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via either
/// Lax-Friedrichs or Riemann solver (linear and/or exact).  Adds viscosity if
/// asked for, and tracks flux of N passive tracers.
///
/// - 2009-12-18 JM: Added Axisymmetric Class (cyl_FV_solver_Hydro_Euler)
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux
/// functions).
/// - 2010.10.01 JM: Added spherical coordinate system.
/// - 2010.11.03 JM: Fixed source terms for spherical coordinates
///   (although results don't change much).
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.22 JM: Added new Riemann solver classes for Hydro.
/// - 2010.12.23 JM: Removed references to riemann_base class.
///    added extra variable to inviscid_flux() function.
///    Moved UtoP() etc. from solver to flux-solver.
/// - 2010.12.28 JM: removed some debugging comments.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.08.19 JM: tested a bunch of approximations, but nothing
///    was an improvement so I left it the way it was.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.04.14 JM: Moved flux solver to FV_solver

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "solver_eqn_hydro_adi.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

using namespace std;

// *********************************
// ***** FV SOLVER HYDRO EULER *****
// *********************************

// ##################################################################
// ##################################################################

FV_solver_Hydro_Euler::FV_solver_Hydro_Euler(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt *state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr), eqns_Euler(nv),
    riemann_Euler(nv, state, gam), Riemann_FVS_Euler(nv, gam),
    Riemann_Roe_Hydro_PV(nv, gam), Riemann_Roe_Hydro_CV(nv, gam),
    HLL_hydro(nv, gam), VectorOps_Cart(nd)
{
  spdlog::debug("FV_solver_Hydro_Euler::FV_solver_Hydro_Euler() constructor");
  spdlog::debug(
      "FV_solver_Hydro_Euler::FV_solver_Hydro_Euler() gamma = {}", eq_gamma);
  counter = 1;
  utemp.resize(nv);
  return;
}

// ##################################################################
// ##################################################################

FV_solver_Hydro_Euler::~FV_solver_Hydro_Euler()
{
  spdlog::debug("FV_solver_Hydro_Euler::~FV_solver_Hydro_Euler() destructor");
  return;
}

// ##################################################################
// ##################################################################

int FV_solver_Hydro_Euler::inviscid_flux(
    class SimParams &par,   ///< simulation parameters
    class GridBaseClass *,  ///< pointer to grid
    const double dx,        ///< cell-size dx (for LF method)
    class cell &Cl,         ///< Left state cell pointer
    class cell &Cr,         ///< Right state cell pointer
    const pion_flt *Pl,     ///< Left Primitive state vector.
    const pion_flt *Pr,     ///< Right Primitive state vector.
    pion_flt *flux,         ///< Resultant Flux state vector.
    pion_flt *pstar,        ///< State vector at interface.
    const int solve_flag,   ///< Solver to use
    const double g          ///< Gas constant gamma.
)
{
#ifdef TEST_INF
  //
  // Check input density and pressure are 'reasonably large'
  //
  if (Pl[eqRO] < TINYVALUE || Pl[eqPG] < TINYVALUE || Pr[eqRO] < TINYVALUE
      || Pr[eqPG] < TINYVALUE) {
    spdlog::error("left  : {}", std::vector<double>(Pl, Pl + eq_nvar));
    spdlog::error("left  : {}", std::vector<double>(Pr, Pr + eq_nvar));
    spdlog::error(
        "FV_solver_Hydro_Euler::calculate_flux() rho/P too small {}", Pr[eqPG]);
    exit_pion(1);
  }

  for (int v = 0; v < eq_nvar; v++)
    if (!isfinite(Pl[v])) spdlog::error("{}: {} {}", "flux hydro Pl", v, Pl[v]);
  for (int v = 0; v < eq_nvar; v++)
    if (!isfinite(Pr[v])) spdlog::error("{}: {} {}", "flux hydro Pr", v, Pr[v]);

  bool problem = false;
  for (int t = 0; t < FV_ntr; t++) {
    if (Pl[eqTR[t]] > 1.0 || Pl[eqTR[t]] < 0.0) problem = true;
    if (Pr[eqTR[t]] > 1.0 || Pr[eqTR[t]] < 0.0) problem = true;
  }
  if (problem) {
    spdlog::error("Tracers out of range in flux solver!");
    spdlog::error("Pl={}", std::vector<double>(Pl, Pl + eq_nvar));
    spdlog::error("Pr={}", std::vector<double>(Pr, Pr + eq_nvar));
    exit_pion(2);
  }
#endif

  int err = 0;
  double ustar[eq_nvar];
  for (int v = 0; v < eq_nvar; v++)
    ustar[v] = 0.0;
  for (int v = 0; v < eq_nvar; v++)
    flux[v] = 0.0;
  for (int v = 0; v < eq_nvar; v++)
    pstar[v] = 0.0;
  eq_gamma = g;

  //
  // Choose which Solver to use.  Each of these must set the values of
  // the flux[] and pstar[] arrays.
  //
  if (solve_flag == FLUX_LF) {
    // Lax-Friedrichs Method, so just get the flux.
    err += get_LaxFriedrichs_flux(Pl, Pr, flux, dx, eq_gamma);
    for (int v = 0; v < eq_nvar; v++)
      pstar[v] = 0.5 * (Pl[v] + Pr[v]);
  }

  else if (solve_flag == FLUX_FVS) {
    // Flux Vector Splitting (van Leer, 1982): This function takes the
    // left and right state and returns the FVS flux in 'flux' and the
    // Roe-average state in 'pstar'.
    err += FVS_flux(Pl, Pr, flux, pstar, eq_gamma);
  }

  else if (
      solve_flag == FLUX_RSlinear || solve_flag == FLUX_RSexact
      || solve_flag == FLUX_RShybrid) {
    // These are all Riemann Solver Methods, so call the solver:
    err += JMs_riemann_solve(Pl, Pr, pstar, solve_flag, eq_gamma);
    PtoFlux(pstar, flux, eq_gamma);
    if (err != 0) {
      spdlog::error("Error Flux Riemann Solver: {}", err);
      exit_pion(1);
    }
  }

  else if (solve_flag == FLUX_RSroe) {
    // Roe conserved variables flux solver (Toro 1999), using the
    // symmetric calculation.
    err += Roe_flux_solver_symmetric(Pl, Pr, eq_gamma, HC_etamax, pstar, flux);
  }

  else if (solve_flag == FLUX_RSroe_pv) {
    // Roe primitive variables linear solver:
    err += Roe_prim_var_solver(Pl, Pr, eq_gamma, pstar);
    // Convert pstar to a flux:
    PtoFlux(pstar, flux, eq_gamma);
  }

  // This is like the MHD HLLD/HLL solver, including compressive motion check
  // and strong-gradient zones check (Migone et al. 2011 ), except that it uses
  // the hydro Roe-CV solver instead of HLLD.  We only check for density ratio
  // here, to smooth out contact discontinuities that are too strong
  else if (solve_flag == FLUX_RCV_HLL) {
    // double dr = Pl[RO] / Pr[RO], drlim = 4.0;
    // if ((dr > drlim) || (dr < 1.0 / drlim)) {
    if (Cl.hll || Cr.hll) {
      // strong density jump across the cell
      // HLL solver -- Miyoshi and Kusano (2005) (m05)
      // spdlog::info("HLL! Cl {} {}, Cr {} {}",Cl.id, Cl.hll, Cr.id, Cr.hll);
      err += hydro_HLL_flux_solver(Pl, Pr, eq_gamma, flux, ustar);
      err += UtoP(ustar, pstar, par.EP.MinTemperature, eq_gamma);
    }
    else {
      // Roe conserved variables flux solver (Toro 1999), using the
      // symmetric calculation.
      err +=
          Roe_flux_solver_symmetric(Pl, Pr, eq_gamma, HC_etamax, pstar, flux);
    }
    if (0 != err) {
      spdlog::error("HLL/RoeCV Flux {}", err);
      exit_pion(1);
    }
  }

  // HLL solver, very diffusive 2 wave solver (Migone et al. 2011 )
  else if (solve_flag == FLUX_RS_HLL) {
    err += hydro_HLL_flux_solver(Pl, Pr, eq_gamma, flux, ustar);
    // HLL solver returns Ustar, so convert to Pstar
    err += UtoP(ustar, pstar, par.EP.MinTemperature, eq_gamma);
  }

  else {
    spdlog::error("what sort of flux solver do you mean? {}", solve_flag);
    exit_pion(1);
  }

  return err;
}



// ##################################################################
// ##################################################################



void FV_solver_Hydro_Euler::PtoU(const pion_flt *p, pion_flt *u, const double g)
{
  for (int t = 0; t < FV_ntr; t++)
    u[eqTR[t]] = p[eqTR[t]] * max(MIN_DENS, p[eqRO]);
  eqns_Euler::PtoU(p, u, g);
  return;
}



// ##################################################################
// ##################################################################



int FV_solver_Hydro_Euler::UtoP(
    const pion_flt *u,
    pion_flt *p,
    const double MinTemp,  ///< minimum temperature/pressure allowed
    const double g)
{
  for (int t = 0; t < FV_ntr; t++)
    p[eqTR[t]] = u[eqTR[t]] / max(MIN_DENS, u[eqRHO]);
  int err = eqns_Euler::UtoP(u, p, MinTemp, g);
  return err;
}



// ##################################################################
// ##################################################################



void FV_solver_Hydro_Euler::PUtoFlux(
    const pion_flt *p, const pion_flt *u, pion_flt *f)
{
  for (int t = 0; t < FV_ntr; t++)
    f[eqTR[t]] = p[eqTR[t]] * f[eqRHO];
  eqns_Euler::PUtoFlux(p, u, f);
  return;
}



// ##################################################################
// ##################################################################



void FV_solver_Hydro_Euler::UtoFlux(
    const pion_flt *u, pion_flt *f, const double g)
{
  eqns_Euler::UtoFlux(u, f, g);
  for (int t = 0; t < FV_ntr; t++)
    f[eqTR[t]] = u[eqTR[t]] * f[eqRHO] / u[eqRHO];
  return;
}



// ##################################################################
// ##################################################################



int FV_solver_Hydro_Euler::AVFalle(
    const pion_flt *Pl,
    const pion_flt *Pr,
    const pion_flt *pstar,
    pion_flt *flux,
    const double etav,
    const double gam)
{
  /// \section Equations
  /// Currently the equations are as follows (for flux along x-direction):
  /// \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*}_f (v_{xR}-v_{xL})  \,,\f]
  /// \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*}_f (v_{yR}-v_{yL})  \,,\f]
  /// \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*}_f (v_{zR}-v_{zL})  \,,\f]
  /// \f[ F[e]   = F[e]   - \eta \rho^{*} c^{*}_f \left( v^{*}_x
  /// (v_{xR}-v_{xL})
  /// + v^{*}_y (v_{yR}-v_{yL}) + v^{*}_z (v_{zR}-v_{zL}) \right)\,.\f] This
  /// follows from assuming a 1D problem, with non-zero bulk viscosity, and a
  /// non-zero shear viscosity.  The bulk viscosity gives the viscosity in the
  /// x-direction, opposing stretching and compression.  The shear viscosity
  /// given the y and z terms, opposing slipping across boundaries.
  ///
  /// The identification of the velocities and density with the
  /// resultant state from the Riemann Solver is a bit ad-hoc, and could
  /// easily be changed to using the mean of the left and right states, or
  /// something else. (See my notes on artificial viscosity).
  ///
  double prefactor = maxspeed(pstar, gam) * etav * pstar[eqRO];

  //
  // x-dir first.
  //
  double momvisc = prefactor * (Pr[eqVX] - Pl[eqVX]);
  double ergvisc = momvisc * pstar[eqVX];
  flux[eqMMX] -= momvisc;

  // y-dir
  momvisc = prefactor * (Pr[eqVY] - Pl[eqVY]);
  flux[eqMMY] -= momvisc;
  ergvisc += momvisc * pstar[eqVY];

  // z-dir
  momvisc = prefactor * (Pr[eqVZ] - Pl[eqVZ]);
  flux[eqMMZ] -= momvisc;
  ergvisc += momvisc * pstar[eqVZ];

  // update energy flux
  flux[eqERG] -= ergvisc;
  return 0;
}



// ##################################################################
// ##################################################################



///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_Hydro_Euler::dU_Cell(
    class GridBaseClass *grid,
    cell &c,                             // Current cell.
    const axes d,                        // Which axis we are looking along.
    const std::vector<pion_flt> &fn,     // Negative direction flux.
    const std::vector<pion_flt> &fp,     // Positive direction flux.
    const std::vector<pion_flt> &slope,  // slope vector for cell c.
    const int cstep,                     // spatial order of accuracy.
    const double dx,                     // cell length dx.
    const double dt                      // cell TimeStep, dt.
)
{
  if (!c.isdomain) return 0;
  // This calculates -dF/dx
  // utemp is overwritten, so don't need to initialize
  int err = DivStateVectorComponent(c, grid, d, eq_nvar, fn, fp, utemp);
  // add source terms
  geometric_source(c, d, slope, cstep, dx, utemp);
  wind_acceleration_source(grid, c, d, fn, fp, cstep, utemp);
  for (int v = 0; v < eq_nvar; v++)
    c.dU[v] += dt * utemp[v];
  return (err);
}



// ##################################################################
// ##################################################################



int FV_solver_Hydro_Euler::CellAdvanceTime(
    class cell &c,
    const pion_flt *Pin,  // Initial State Vector.
    pion_flt *dU,         // Update vector dU
    pion_flt *Pf,         // Final state vector (can be same as initial vec.).
    pion_flt *dE,  // NDEBUG Tracks change of energy for negative pressure
                   // correction.
    const double,  // gas EOS gamma.
    const double MinTemp,  ///< Min Temperature allowed on grid.
    const double           // Cell timestep dt.
)
{
  if (!c.isdomain) return 0;
  //
  // First convert from Primitive to Conserved Variables
  //
  pion_flt u1[eq_nvar];
  pion_flt Pintermediate[eq_nvar];
  pion_flt corrector[eq_nvar];
  // for (int v=0;v<eq_nvar;v++) corrector[v]=1.0;
  if (mp) {
    mp->sCMA(corrector, Pin);
    for (int t = 0; t < eq_nvar; t++)
      Pintermediate[t] = Pin[t] * corrector[t];
    PtoU(Pintermediate, u1, eq_gamma);
  }
  else {
    PtoU(Pin, u1, eq_gamma);
  }

  //
  // Now add dU[] to U[], and change back to primitive variables.
  // This can give negative pressures, so check for that and fix it if needed.
  //
  for (int v = 0; v < eq_nvar; v++) {
    u1[v] += dU[v];  // Update conserved variables
  }

  if (u1[RHO] < 0.0) {
#ifndef NDEBUG
    spdlog::info(
        "Pc {}", std::vector<double>(Pintermediate, Pintermediate + eq_nvar));
    spdlog::info("U  {}", std::vector<double>(u1, u1 + eq_nvar));
    spdlog::info("dU {}", std::vector<double>(dU, dU + eq_nvar));
    CI.print_cell(c);
#endif
    spdlog::warn(
        "CellAdvanceTime: Reset -ve density {:12.6e} to original {:12.6e}",
        u1[RHO], Pin[RO]);
    spdlog::warn("PION will exit if this is a leaf cell on the domain");
    if (c.isleaf && c.isdomain && !c.isbd && c.isgd) {
      spdlog::error("negative density in leaf grid cell, bugging out.");
      CI.print_cell(c);
      exit_pion(1);
    }
    for (int t = 0; t < FV_ntr; t++)
      u1[eqTR[t]] *= Pin[RO] / u1[RHO];
    u1[RHO] = Pin[RO];
  }

  int err;
  if ((err = UtoP(u1, Pf, MinTemp, eq_gamma)) != 0) {
#ifndef NDEBUG
    spdlog::error("(FV_solver_Hydro_Euler::CellAdvanceTime) UtoP complained"
                  " (maybe about negative pressure...) fixing");
    spdlog::debug("pin : {}", std::vector<double>(Pin, Pin + eq_nvar));
    spdlog::debug("dU  : {}", std::vector<double>(dU, dU + eq_nvar));
    spdlog::debug("u1  : {}", std::vector<double>(u1, u1 + eq_nvar));
    PtoU(Pin, u1, eq_gamma);
    spdlog::debug("Uin : {}", std::vector<double>(u1, u1 + eq_nvar));
    spdlog::debug("Pf  : {}", std::vector<double>(Pf, Pf + eq_nvar));
    exit_pion(3);
#endif
  }

  // Reset the dU array for the next timestep.
  for (int v = 0; v < eq_nvar; v++)
    dU[v] = 0.;

#ifdef TEST_INF
  for (int v = 0; v < eq_nvar; v++) {
    if (!isfinite(Pf[v])) {
      spdlog::debug(
          "NAN/INF in FV_solver_Hydro_Euler::CellAdvanceTime: var={}, val={}",
          v, Pf[v]);
      CI.print_cell(c);
      spdlog::error("NAN hydro cell update {} {}", v, Pf[v]);
      exit_pion(1);
    }
  }
#endif
  // int final_correct_flag = 0;
  // for (int v=0;v<eq_nvar;v++) corrector[v]=1.0;
  if (mp) {
    mp->sCMA(corrector, Pf);
    for (int t = 0; t < eq_nvar; t++)
      Pf[t] = Pf[t] * corrector[t];
  }

  return err;
}

// ##################################################################
// ##################################################################

///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_Hydro_Euler::CellTimeStep(
    const cell &c,   ///< pointer to cell
    const double,    ///< gas EOS gamma.
    const double dx  ///< Cell size dx.
)
{
  /// \section Algorithm
  /// First Get the maximum fluid velocity in each of the three directions.
  /// Add the sound speed to the abs. value of this to get the max signal
  /// propagation speed.  Then dt=dx/max.vel.
  /// Then multiply by the CFl no. and return.
  ///

  pion_flt temp = 0.0, l_dt = 0.0;
  for (int v = 0; v < FV_gndim; v++)
    temp += c.P[eqVX + v] * c.P[eqVX + v];
  temp = sqrt(temp);
  temp += chydro(c.P.data(), eq_gamma);

  //
  // Get Max velocity along a grid direction.
  //
  // pion_flt temp = fabs(c->P[eqVX]);
  // if (FV_gndim>1) temp = max(temp,static_cast<pion_flt>(fabs(c->P[eqVY])));
  // if (FV_gndim>2) temp = max(temp,static_cast<pion_flt>(fabs(c->P[eqVZ])));
  //
  // Add the sound speed to this, and it is the max wavespeed.
  //
  // temp += chydro(c->P,eq_gamma);
  l_dt = dx / temp;
  //
  // Now scale the max. allowed timestep by the CFL number we are using (<1).
  //
  l_dt *= FV_cfl;

#ifdef TEST_INF
  if (!isfinite(l_dt) || l_dt <= 0.0) {
    spdlog::error("cell has invalid timestep");
    CI.print_cell(c);
    exit_pion(5);
  }
#endif
  return l_dt;
}

// ##################################################################
// ##################################################################

////////////////////////////////////////////////////////////////////////
/// CYLINDRICAL (AXISYMMETRIC) COORDINATES
////////////////////////////////////////////////////////////////////////

// ##################################################################
// ##################################################################

cyl_FV_solver_Hydro_Euler::cyl_FV_solver_Hydro_Euler(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt *state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr), eqns_Euler(nv),
    riemann_Euler(nv, state, gam), Riemann_FVS_Euler(nv, gam),
    Riemann_Roe_Hydro_PV(nv, gam), Riemann_Roe_Hydro_CV(nv, gam),
    HLL_hydro(nv, gam), VectorOps_Cart(nd),
    FV_solver_Hydro_Euler(nv, nd, cflno, gam, state, avcoeff, ntr),
    VectorOps_Cyl(nd)
{
  if (nd != 2) {
    spdlog::error("Cylindrical coordinates only 2d", nd);
    exit_pion(nd);
  }
  return;
}

// ##################################################################
// ##################################################################

cyl_FV_solver_Hydro_Euler::~cyl_FV_solver_Hydro_Euler() {}

// ##################################################################
// ##################################################################

void cyl_FV_solver_Hydro_Euler::geometric_source(
    cell &c,                            ///< Current cell.
    const axes d,                       ///< Which axis we are looking along.
    const std::vector<pion_flt> &dpdx,  ///< slope vector for cell c.
    const int OA,                       ///< spatial order of accuracy.
    const double dR,                    ///< cell length dx.
    std::vector<pion_flt> &dU  ///< update vector to add source term to [OUTPUT]
)
{

  if (d == Rcyl) {
    switch (OA) {
      case OA1:
        dU[eqMMX] += c.Ph[eqPG] / CI.get_dpos(c, Rcyl);
        break;
      case OA2:
        dU[eqMMX] +=
            (c.Ph[eqPG] + (CI.get_dpos(c, Rcyl) - R_com(c, dR)) * dpdx[eqPG])
            / CI.get_dpos(c, Rcyl);
        break;
      default:
        spdlog::error("Bad OOA in cyl_IdealMHD_RS::dU: {}", OA);
        exit_pion(7);
        break;
    }
  }

  return;
}

// ##################################################################
// ##################################################################

////////////////////////////////////////////////////////////////////////
/// SPHERICAL (POLAR) COORDINATES
////////////////////////////////////////////////////////////////////////

// ##################################################################
// ##################################################################

sph_FV_solver_Hydro_Euler::sph_FV_solver_Hydro_Euler(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt *state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr), eqns_Euler(nv),
    riemann_Euler(nv, state, gam), Riemann_FVS_Euler(nv, gam),
    Riemann_Roe_Hydro_PV(nv, gam), Riemann_Roe_Hydro_CV(nv, gam),
    HLL_hydro(nv, gam), VectorOps_Cart(nd),
    FV_solver_Hydro_Euler(nv, nd, cflno, gam, state, avcoeff, ntr),
    VectorOps_Cyl(nd), VectorOps_Sph(nd)
{
  if (nd != 1) {
    spdlog::error("Spherical coordinates only 1D {}", nd);
    exit_pion(nd);
  }
  return;
}

// ##################################################################
// ##################################################################

sph_FV_solver_Hydro_Euler::~sph_FV_solver_Hydro_Euler() {}

// ##################################################################
// ##################################################################

void sph_FV_solver_Hydro_Euler::geometric_source(
    cell &c,                            ///< Current cell.
    const axes d,                       ///< Which axis we are looking along.
    const std::vector<pion_flt> &dpdx,  ///< slope vector for cell c.
    const int OA,                       ///< spatial order of accuracy.
    const double dR,                    ///< cell length dx.
    std::vector<pion_flt> &dU  ///< update vector to add source term to [OUTPUT]
)
{

  if (d == Rsph) {
    switch (OA) {
      //
      case OA1:
        dU[eqMMX] += 2.0 * c.Ph[eqPG] / R3(c, dR);
        break;
      //
      case OA2:
        dU[eqMMX] += 2.0
                     * ((c.Ph[eqPG] - dpdx[eqPG] * R_com(c, dR)) / R3(c, dR)
                        + dpdx[eqPG]);
        break;
      //
      default:
        spdlog::error("Bad OOA in sph_hydro_RS::dU: {}", OA);
        exit_pion(9);
        break;
    }
  }

  //#define GRAVITY
#ifdef GRAVITY
  //
  // Here we impose an external gravitational acceleration, with a
  // radial power law and a normalisation K=G.M, so that: a=-K/r^alpha
  //
#define GRAV_NORM 1.33e29  // 1.989e37*6.67e-8 = 10^4 Msun*G
#define GRAV_ALPHA 2.0     // point mass has 1/r potential, 1/r^2 acc.

  //
  // Mean value of k.r^{-a} in a cell: for a!=3, cell centred at r=ri,
  // and cell diameter delta, this is
  // k.[(ri+delta/2)^{3-alpha}-(ri-delta/2)^{3-alpha}]/[(3-alpha)*delta*ri*R3]
  //
  // Note this is the same for 1st and 2nd order because we are
  // calculating it from the exact potential integrated through the
  // cell.
  //
  if (d == Rsph) {
    double rc   = CI.get_dpos(c, Rsph);
    double temp = c->Ph[eqRO] * GRAV_NORM / (rc * R3(c));  // GM/r^2
    if (rc < 0.0) temp *= -1.0;  // so that we always point towards the origin!
    // cout <<"  r="<<rc<<"  GM/r^2="<<temp;
    // cout <<" d(MMX)/d(MMX,grav) = "<<c->dU[eqMMX]<<" / "<< -FV_dt*temp;
    c->dU[eqMMX] -= temp;
    // cout <<" d(ERG)/d(ERG,grav) = "<<c->dU[eqERG]<<" / "<<
    // -FV_dt*temp*c->Ph[eqVX]<<"\n";
    c->dU[eqERG] -= c->Ph[eqVX] * temp;
  }

#endif  // GRAVITY
  return;
}

// ##################################################################
// ##################################################################
