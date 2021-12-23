//      references:
//      * Stone et al. 2008     (s08)
//      * Miyoshi and Kusano 2005 (m05)
//      * Migone and Bodo 2008  (mb08)

#include "Riemann_solvers/HLL_hydro.h"
#include "Riemann_solvers/riemann.h"
#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */


#ifndef NDEBUG
#endif  // NDEBUG

#include "equations/eqns_hydro_adiabatic.h"
#include "microphysics/microphysics_base.h"
using namespace std;

// ##################################################################
// ##################################################################

HLL_hydro::HLL_hydro(
    const int nv,          ///< Length of State Vectors, nvar
    const double eq_gamma  ///< Gamma for state vector.
    ) :
    eqns_base(nv),
    eqns_Euler(nv)
{
  spdlog::debug("(HLL_hydro::HLL_hydro) Initialising HLL Solver Class");
  if (eq_nvar < 5) {
    spdlog::error("{}: {}", "#elements!=5, QUIT.", eq_nvar);
  }

  HD_lambda = mem.myalloc(HD_lambda, 2);    // wave speeds (only 2 for HLL)
  HD_UL     = mem.myalloc(HD_UL, eq_nvar);  // conserved
  HD_UR     = mem.myalloc(HD_UR, eq_nvar);
  HD_FL     = mem.myalloc(HD_FL, eq_nvar);  // flux
  HD_FR     = mem.myalloc(HD_FR, eq_nvar);

  spdlog::debug("(HLL_hydro::HLL_hydro) All set");
  return;
}

// ##################################################################
// ##################################################################

HLL_hydro::~HLL_hydro()
{

  spdlog::debug("(HLL_hydro::HLL_hydro) Commencing Destruction");

  HD_lambda = mem.myfree(HD_lambda);  // wave speeds (one entropy, two Alfvén)
  HD_UL     = mem.myfree(HD_UL);      // conserved
  HD_UR     = mem.myfree(HD_UR);
  HD_FL     = mem.myfree(HD_FL);  // flux
  HD_FR     = mem.myfree(HD_FR);

  spdlog::debug("(riemann_MHD::riemann_MHD) Mission Accomplished");
}

// ###################################################################
// ###################################################################

void HLL_hydro::HLL_signal_speeds(
    const pion_flt *Pl,  ///< inputs
    const pion_flt *Pr,
    const double eq_gamma,
    double &Sl,  ///< outputs
    double &Sr)
{
  //
  // compute wave speeds (m05 eq 3)
  //
  double cf_l   = chydro(Pl, eq_gamma);
  double cf_r   = chydro(Pr, eq_gamma);
  double cf_max = max(cf_l, cf_r);
  Sl            = min(Pl[eqVX], Pr[eqVX]) - cf_max;  // (m05 eq 67)
  Sr            = max(Pl[eqVX], Pr[eqVX]) + cf_max;  //
  return;
}

// ###################################################################
// ###################################################################

int HLL_hydro::hydro_HLL_flux_solver(
    const pion_flt *Pl,     ///< input left state
    const pion_flt *Pr,     ///< input right state
    const double eq_gamma,  ///< input gamma
    pion_flt *out_flux,     ///< output flux
    pion_flt *out_ustar     ///< output interface state (cons.var.)
)
{

  //
  // compute conserved U and Flux (m05 eq 3)
  //
  PtoU(Pl, HD_UL, eq_gamma);
  PtoU(Pr, HD_UR, eq_gamma);

  PUtoFlux(Pl, HD_UL, HD_FL);
  PUtoFlux(Pr, HD_UR, HD_FR);

  //
  // compute wave speeds (m05 eq 3)
  //
  HLL_signal_speeds(
      Pl, Pr, eq_gamma, HD_lambda[0],
      HD_lambda[1]);  // Pl,Pr,g,Sl,Sr
  if (HD_lambda[0] > 0) {
    for (int v = 0; v < eq_nvar; v++)
      out_flux[v] = HD_FL[v];
  }
  else if (HD_lambda[1] < 0) {
    for (int v = 0; v < eq_nvar; v++)
      out_flux[v] = HD_FR[v];
  }
  else {
    for (int v = 0; v < eq_nvar; v++) {
      out_flux[v] = (HD_lambda[1] * HD_FL[v] - HD_lambda[0] * HD_FR[v]
                     + HD_lambda[1] * HD_lambda[0] * (HD_UR[v] - HD_UL[v]))
                    / (HD_lambda[1] - HD_lambda[0]);
    }
  }

  // inteface state: first in conserved variables.
  for (int v = 0; v < eq_nvar; v++) {
    out_ustar[v] = (HD_lambda[1] * HD_UR[v] - HD_lambda[0] * HD_UL[v] + HD_FL[v]
                    - HD_FR[v])
                   / (HD_lambda[1] - HD_lambda[0]);
  }

  return 0;
}

// ###################################################################
// ###################################################################
