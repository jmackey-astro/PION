/// \file radiative_shock.cc
/// \author Jonathan Mackey
///
/// File for setting up radiative shock test problems.
///
/// - 2010-04-09 JM: added support for upstream and downstream tracer
/// variables in the RadiativeShockOutflow case.
/// - 2013.01.11 JM: Added tracer variables for RSH test.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "icgen.h"
#include "ics/icgen_base.h"
#include <sstream>

// ##################################################################
// ##################################################################

IC_radiative_shock::IC_radiative_shock()
{
  gg  = 0;
  rp  = 0;
  vsh = rho0 = T0 = gam = 0.0;
  return;
}

// ##################################################################
// ##################################################################

IC_radiative_shock::~IC_radiative_shock()
{
  return;
}

// ##################################################################
// ##################################################################

int IC_radiative_shock::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) spdlog::error("{}: {}", "null pointer to grid!", fmt::ptr(ggg));

  ICsetup_base::rp = rrp;
  if (!rp) spdlog::error("{}: {}", "null pointer to ReadParams", fmt::ptr(rp));

  string seek, str;

  IC_radiative_shock::eqns = SimPM->eqntype;
  if (eqns == EQEUL)
    eqns = 1;
  else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
    eqns = 2;
  else
    spdlog::error("{}: {}", "Bad equations", eqns);

  // Find shock velocity (cm/s)
  seek = "RADSH_vs";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_radiative_shock::vsh = atof(str.c_str());

  // Find initial gas density(g/cm^3)
  seek = "RADSH_r0";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_radiative_shock::rho0 = atof(str.c_str());

  // Find initial gas temperature (K)
  seek = "RADSH_T0";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_radiative_shock::T0 = atof(str.c_str());

  // Find initial gas B-field (Gauss)
  seek = "RADSH_B0";
  str  = rp->find_parameter(seek);
  if (str == "") IC_radiative_shock::B0 = 0.0;
  IC_radiative_shock::B0 = atof(str.c_str());
#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  B0 /= sqrt(4.0 * M_PI);
#endif

  IC_radiative_shock::gam = SimPM->gamma;

  // now make sure we are to do a radiative shock sim...
  string ics = rp->find_parameter("ics");

  if (ics == "")
    spdlog::error("{}: {}", "didn't get any ics to set up.", ics);
  else if (ics == "RadiativeShock") {
    spdlog::debug("\t\tsetting up test problem: {}", ics);
    err += setup_RadiativeShock();
  }
  else if (ics == "RadiativeShockOutflow") {
    spdlog::debug("\t\tsetting up test problem: {}", ics);
    err += setup_OutflowRadiativeShock();
  }
  else
    spdlog::error("{}: {}", "Don't know what Initial Condition is!", ics);

  // Add noise to data?  Smooth data?
  double noise = 0.0;
  int smooth   = 0;
  ics          = rp->find_parameter("noise");
  if (ics != "")
    noise = atof(ics.c_str());
  else
    noise = -1;
  if (isnan(noise))
    spdlog::error("{}: {}", "noise parameter is not a number", noise);
  if (noise > 0) err += AddNoise2Data(gg, *SimPM, 2, noise);

  ics = rp->find_parameter("smooth");
  if (ics != "")
    smooth = atoi(ics.c_str());
  else
    smooth = -1;
  if (isnan(smooth))
    spdlog::error("{}: {}", "Smooth parameter not a number", smooth);
  if (smooth > 0) err += SmoothData(smooth);

  return err;
}

// ##################################################################
// ##################################################################

int IC_radiative_shock::setup_RadiativeShock()
{

  spdlog::debug(
      "\t\tSetting up radiative shock problem with v={}, rho={}, T={}", vsh,
      rho0, T0);
  double mu = 1.27;   // 1.2; // mean mass per particle -- rough guess, good
                      // for neutral H, He.
  double x  = 0.101;  // initial ionisation fraction...
  double pg = rho0 * (1. + x) * pconst.kB() * T0 / mu / pconst.m_p();

  //
  // Tracer values: upstream and downstream (only upstream used for Radiative
  // Shock)
  //
  int ntr = SimPM->ntracer;
  double trup[ntr];
  string seek, str;
  for (int t = 0; t < ntr; t++) {
    ostringstream temp;
    temp.str("");
    temp << "RADSH_upTR" << t;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      trup[t] = atof(str.c_str());
    else
      trup[t] = 0.5;
  }

  class cell *c = gg->FirstPt();
  do {
    c->P[RO] = rho0;
    c->P[PG] = pg;
    c->P[VX] = -vsh;
    c->P[VY] = c->P[VZ] = 0.0;
    // magnetic field???
    if (eqns == 2) {
      c->P[BY] = B0;
      c->P[BX] = c->P[BZ] = 0.0;
    }
    // tracers (fractional abundances!)
    for (int i = 0; i < SimPM->ntracer; i++)
      c->P[SimPM->ftr + i] = trup[i];  // set later in MP.
                                       // done.
  } while ((c = gg->NextPt(c)) != 0);
  return 0;
}

// ##################################################################
// ##################################################################

int IC_radiative_shock::setup_OutflowRadiativeShock()
{
  spdlog::debug(
      "\t\tSetting up OUTFLOW radiative shock problem with v={}, rho={}, T={}",
      vsh, rho0, T0);

  double mu =
      1.22;  // mean mass per particle -- rough guess, good for neutral H, He.
  // mu /=2.0; // for ionised gas.
  double pg        = rho0 * pconst.kB() * T0 / mu / pconst.m_p();
  double xboundary = (SimPM->Xmax[XX] - SimPM->Xmin[XX]) / 5.;
  if (vsh <= 1.01e7)
    xboundary *= 2.5;  // stable shock should be near centre of grid.
  double range = (SimPM->Xmax[XX] - SimPM->Xmin[XX]) * 5.0 / SimPM->NG[XX];
  double mach0 = vsh / sqrt(gam * pg / rho0);
  spdlog::debug("shock mach no. = {}", mach0);
  double rho1 = rho0 * mach0 * mach0;  // isothermal shock jump condition.
  //  double vel1 = rho0*vsh/rho1;

  // This divisor factor is basically arbitrary -- I have used it at 3.0,
  // 10.0, 30.0 The idea is to get it so that it produces a reasonable shock
  // i.e. that the shock doesn't move off the grid one way or the other too
  // quickly, and that the waves advected downstream aren't too seriously
  // messing the solution.  3.0 seems to be a good value for Raymond's 1979
  // model E shock (100km/s).
  string seek, str;
  seek           = "RADSH_divisor";
  str            = rp->find_parameter(seek);
  double divisor = 0.0;
  if (str != "")
    divisor = atof(str.c_str());
  else
    divisor = 3.0;
  if (eqns == 2) rho1 /= divisor;
  double pg1 = rho1 * pconst.kB() * T0 / mu / pconst.m_p();

  //
  // Tracer values: upstream and downstream
  //
  int ntr = SimPM->ntracer;
  double trup[ntr], trdn[ntr];
  for (int t = 0; t < ntr; t++) {
    ostringstream temp;
    temp.str("");
    temp << "RADSH_upTR" << t;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      trup[t] = atof(str.c_str());
    else
      trup[t] = 0.0;

    temp.str("");
    temp << "RADSH_dnTR" << t;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      trdn[t] = atof(str.c_str());
    else
      trdn[t] = 0.0;
  }

  class cell *c = gg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  do {
    CI.get_dpos(c, dpos);
    if (dpos[XX] >= xboundary + range) {
      //
      // Upstream of shock -- low density fast gas flowing onto domain.
      //
      c->P[RO] = rho0;
      c->P[PG] = pg;
      c->P[VX] = -vsh;
      c->P[VY] = c->P[VZ] = 0.0;
      // magnetic field???
      if (eqns == 2) {
        c->P[BY] = c->P[BZ] = B0 / sqrt(2.);
        c->P[BX]            = 0.0;
      }
      // tracers (fractional abundances!)
      for (int i = 0; i < SimPM->ntracer; i++)
        c->P[SimPM->ftr + i] = trup[i];
    }
    else if (dpos[XX] <= xboundary - range) {
      //
      // We are in dense boundary layer.
      //
      c->P[RO] = rho1;
      c->P[PG] = pg1;
      c->P[VX] = 0.0;  //-vel1;
      c->P[VY] = c->P[VZ] = 0.0;
      // magnetic field???
      if (eqns == 2) {
        c->P[BY] = c->P[BZ] = B0 * (rho1 / rho0) / sqrt(2.);
        c->P[BX]            = 0.0;
      }
      // tracers (fractional abundances!)
      for (int i = 0; i < SimPM->ntracer; i++)
        c->P[SimPM->ftr + i] = trdn[i];
    }
    else {
      //
      // in boundary region, so we interpolate values.
      //
      double frac = (dpos[XX] - xboundary + range) / 2. / range;
      c->P[RO]    = frac * rho0 + (1. - frac) * rho1;
      c->P[PG]    = frac * pg + (1. - frac) * pg1;
      c->P[VX]    = -frac * vsh;
      // magnetic field???
      if (eqns == 2) {
        c->P[BY] = c->P[BZ] =
            frac * B0 / sqrt(2.) + (1. - frac) * B0 * (rho1 / rho0) / sqrt(2.);
        c->P[BX] = 0.0;
      }
      //
      // tracers (fractional abundances!) -- linearly interpolate
      // between upstream and downstream values.
      //
      for (int i = 0; i < SimPM->ntracer; i++)
        c->P[SimPM->ftr + i] = frac * trup[i] + (1. - frac) * trdn[i];
    }
    // done.
  } while ((c = gg->NextPt(c)) != 0);
  return 0;
}

// ##################################################################
// ##################################################################
