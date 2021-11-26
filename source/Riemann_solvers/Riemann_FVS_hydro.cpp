///
/// \file Riemann_FVS_hydro.cc
/// \author Jonathan Mackey
///
/// Flux vector splitting for the Euler Equations.  This is a fairly
/// self-contained class, depending only on the Euler Equations class.
/// It doesn't use the class variables from the other Riemann solvers.
///
/// References:
///  - Van Leer, B., 1982, LNP, 170, 507.
///    "Flux vector splitting for the Euler equations."
///  - Van Leer, B., 1991, "Flux Vector Splitting for the 1990s".
///  - http://www.chimeracfd.com/programming/gryphon/fluxvanleer.html
///    This ref. has a typo in eq.24 -- the f_b^{pm} term should only have one
///    +/-
///
/// History:
/// - 2010-09-21 JM: Started writing.
/// - 2010.12.22 JM: Moved to Riemann solvers sub-directory.
/// - 2010.12.23 JM: Removed tracer flux calculation -- now calculated
///   by calling function in flux_solver_...
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer
///    fluxes are dealt with by the flux-solver classes.
///    For this to work I had to explicitly call the Euler Eqns class for PtoU()
///    and PUtoFlux(); otherwise the flux-solver class (with tracers) would be
///    used.
/// - 2015.01.14 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include "Riemann_solvers/Riemann_FVS_hydro.h"

using namespace std;

// ##################################################################
// ##################################################################

Riemann_FVS_Euler::Riemann_FVS_Euler(
    const int nv,   ///< Length of State Vectors, nvar
    const double g  ///< Gamma for state vector.
    ) :
    eqns_base(nv),
    eqns_Euler(nv), rs_nvar(5)
{
  //
  // eq_gamma, eq_nvar are defined in eqns_base class
  //
  eq_gamma = g;

  //
  // Allocate memory for fpos,fneg
  //
  fpos = 0;
  fneg = 0;
  fpos = mem.myalloc(fpos, rs_nvar);
  fneg = mem.myalloc(fneg, rs_nvar);
}

// ##################################################################
// ##################################################################

Riemann_FVS_Euler::~Riemann_FVS_Euler()
{
  fpos = mem.myfree(fpos);
  fneg = mem.myfree(fneg);
}

// ##################################################################
// ##################################################################

int Riemann_FVS_Euler::FVS_flux(
    const pion_flt *pl,  ///< Left Primitive var. state vec.
    const pion_flt *pr,  ///< Right Primitive var. state vec.
    pion_flt *flux,      ///< Result Flux vector.
    pion_flt *pstar,     ///< Interface state (for viscosity)
    const double         ///< Gas constant gamma (unused)
)
{
#ifdef TEST_INF
  //
  // Check that inputs make sense!
  //
  for (int v = 0; v < rs_nvar; v++) {
    if (!isfinite(pl[v]) || !isfinite(pr[v])) {
      spdlog::info("NAN's detected FVS flux solver");
      rep.printVec("left ", pl, rs_nvar);
      rep.printVec("right", pr, rs_nvar);
      return 1;
    }
  }
#endif

  //
  // We need the left and right sound speeds and mach numbers,
  // and declare the variables f1,f2 for the split flux calculations.
  //
  double cl = chydro(pl, eq_gamma), cr = chydro(pr, eq_gamma),
         Ml = pl[eqVX] / cl, Mr = pr[eqVX] / cr, f1 = 0.0, f2 = 0.0;

  //
  // First the positive flux (Table 1, van Leer, 1982)
  // Supersonic flux in negative dir --> zero flux.
  // Supersonic flux in positive dir --> pure advection.
  // Else --> Use split flux formulae.
  //
  if (Ml < -1.0) {
    for (int v = 0; v < rs_nvar; v++)
      fpos[v] = 0.0;
  }
  else if (Ml > 1.0) {
    pion_flt utemp[rs_nvar];
    eqns_Euler::PtoU(pl, utemp, eq_gamma);
    eqns_Euler::PUtoFlux(pl, utemp, fpos);
  }
  else {
    f1          = 0.25 * pl[eqRO] * cl * (1.0 + Ml) * (1.0 + Ml);
    f2          = cl * ((eq_gamma - 1.0) * Ml + 2);
    fpos[eqRHO] = f1;
    fpos[eqMMX] = f1 * f2 / eq_gamma;
    fpos[eqMMY] = f1 * pl[eqVY];
    fpos[eqMMZ] = f1 * pl[eqVZ];
    fpos[eqERG] = f1
                  * (f2 * f2 * 0.5 / (eq_gamma * eq_gamma - 1.0)
                     + 0.5 * (pl[eqVY] * pl[eqVY] + pl[eqVZ] * pl[eqVZ]));
  }

  //
  // Now the negative flux:
  // Supersonic flux in positive dir --> zero flux.
  // Supersonic flux in negative dir --> pure advection.
  // Else --> Use split flux formulae.
  //
  if (Mr > 1.0) {
    for (int v = 0; v < rs_nvar; v++)
      fneg[v] = 0.0;
  }
  else if (Mr < -1.0) {
    pion_flt utemp[rs_nvar];
    eqns_Euler::PtoU(pr, utemp, eq_gamma);
    eqns_Euler::PUtoFlux(pr, utemp, fneg);
    // PtoFlux(pr,fneg,eq_gamma);
  }
  else {
    f1          = -0.25 * pr[eqRO] * cr * (1.0 - Mr) * (1.0 - Mr);
    f2          = cr * ((eq_gamma - 1.0) * Mr - 2);
    fneg[eqRHO] = f1;
    fneg[eqMMX] = f1 * f2 / eq_gamma;
    fneg[eqMMY] = f1 * pr[eqVY];
    fneg[eqMMZ] = f1 * pr[eqVZ];
    fneg[eqERG] = f1
                  * (f2 * f2 * 0.5 / (eq_gamma * eq_gamma - 1)
                     + 0.5 * (pr[eqVY] * pr[eqVY] + pr[eqVZ] * pr[eqVZ]));
  }

  //
  // Now add the split fluxes to get the total flux.
  //
  for (int v = 0; v < rs_nvar; v++)
    flux[v] = fpos[v] + fneg[v];

  //
  // Also calculate an average state for the viscosity
  //
  Roe_average_state(pl, pr, eq_gamma, pstar);

#ifdef TEST_INF
  //
  // Make sure the answer is finite
  //
  for (int v = 0; v < rs_nvar; v++) {
    if (!isfinite(fpos[v]) || !isfinite(fneg[v]) || !isfinite(flux[v])) {
      spdlog::info("NAN's detected!!!");
      rep.printVec("fpos ", fneg, rs_nvar);
      rep.printVec("fneg ", fpos, rs_nvar);
      rep.printVec("flux ", flux, rs_nvar);
      return 1;
    }
  }
#endif

  return 0;
}

// ##################################################################
// ##################################################################

void Riemann_FVS_Euler::Roe_average_state(
    const pion_flt *p1,  ///< state 1
    const pion_flt *p2,  ///< state 2
    const double,        ///< gamma (unused)
    pion_flt *ans        ///< Roe-averaged state
)
{
  //
  // Get Roe-average values for rho, v_x, v_y, v_z, H, a
  // Toro eq. 11.60
  //
  double RoeAvg_rl = sqrt(p1[eqRO]), RoeAvg_rr = sqrt(p2[eqRO]),
         RoeAvg_denom = 1.0 / (RoeAvg_rl + RoeAvg_rr);

  ans[eqRO] = RoeAvg_rl * RoeAvg_rr;
  ans[eqVX] = (RoeAvg_rl * p1[eqVX] + RoeAvg_rr * p2[eqVX]) * RoeAvg_denom;
  ans[eqVY] = (RoeAvg_rl * p1[eqVY] + RoeAvg_rr * p2[eqVY]) * RoeAvg_denom;
  ans[eqVZ] = (RoeAvg_rl * p1[eqVZ] + RoeAvg_rr * p2[eqVZ]) * RoeAvg_denom;
  //
  // Get the enthalpy of the mean state: we'll convert back to a
  // pressure below.
  //
  // Enthalpy per unit mass: H= v*v/2 +g*p/(g-1)/rho
  // Hence Adiabatic sound speed a^2 = (H-v*v/2)*(g-1)
  // Then pressure = rho*a^2/g
  //
  ans[eqPG] = RoeAvg_denom
              * (RoeAvg_rl * Enthalpy(p1, eq_gamma)
                 + RoeAvg_rr * Enthalpy(p2, eq_gamma));
  ans[eqPG] = (eq_gamma - 1.0)
              * (ans[eqPG]
                 - 0.5
                       * (ans[eqVX] * ans[eqVX] + ans[eqVY] * ans[eqVY]
                          + ans[eqVZ] * ans[eqVZ]));
  ans[eqPG] = ans[eqRO] * ans[eqPG] / eq_gamma;

  return;
}

// ##################################################################
// ##################################################################
