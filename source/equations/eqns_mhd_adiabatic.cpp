///
/// \file eqns_mhd_adiabatic.cc
///
/// \brief Class definition for eqns_mhd_ideal
/// \author Jonathan Mackey
///
/// This file contains the class definitions for the Ideal MHD Equations,
/// with the GLM divergence cleaning extensions.
///
/// Modifications:
/// - 2007-10-16 Moved the equations classes in here from global.cc
/// - 2009-10-20 renamed eqns_mhd_adiabatic.cc and moved all other equations
/// classes to their own files.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.23 JM: Added SetAvgState() function to eqns_mhd_ideal
///    class.
/// - 2011.01.18 JM: Added BASE_RHO parameter for correcting negative
///   densities.
/// - 2011.04.15 JM: UtoP() again -- added recalculation of all prim.
///    vars if a negative density is encountered.  Also added the
///    ifdef for setting negative pressure to a fixed temperature.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifndef NDEBUG
#endif  // NDEBUG

#include "eqns_mhd_adiabatic.h"
#include "microphysics/microphysics_base.h"

/* prevent clang-format reordering */
#include <fmt/ranges.h>

using namespace std;

/*******************************************************************/
// Member function definitions for eqns_mhd_ideal class

// ##################################################################
// ##################################################################

eqns_mhd_ideal::eqns_mhd_ideal(int nv) : eqns_base(nv)
{
#ifdef FUNCTION_ID
  spdlog::debug("eqns_mhd_ideal::eqns_mhd_ideal ...starting");
#endif  // FUNCTION_ID

#ifndef NDEBUG
  spdlog::debug(
      "(eqns_mhd_ideal::eqns_mhd_ideal) Setting up Ideal MHD Equations class");
#endif
  if (eq_nvar < 8) {
    spdlog::error("\tError!! Class expects (at least) 8 MHD variables");
    exit(1);
  }

#ifdef FUNCTION_ID
  spdlog::debug("eqns_mhd_ideal::eqns_mhd_ideal ...returning.");
#endif  // FUNCTION_ID
}

// ##################################################################
// ##################################################################

eqns_mhd_ideal::~eqns_mhd_ideal()
{
#ifndef NDEBUG
  spdlog::debug(
      "(eqns_mhd_ideal::~eqns_mhd_ideal) Deleting  Ideal MHD Equations class");
#endif
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::PtoU(const pion_flt *p, pion_flt *u, const double gamma)
{
  u[eqRHO] = p[eqRO];
  u[eqMMX] = p[eqRO] * p[eqVX];
  u[eqMMY] = p[eqRO] * p[eqVY];
  u[eqMMZ] = p[eqRO] * p[eqVZ];
  u[eqBBX] = p[eqBX];
  u[eqBBY] = p[eqBY];
  u[eqBBZ] = p[eqBZ];
  // E = p/(g-1) +rho*V^2/2 + B^2/2
  u[eqERG] =
      (p[eqRO] * (p[eqVX] * p[eqVX] + p[eqVY] * p[eqVY] + p[eqVZ] * p[eqVZ])
       * 0.5)
      + (p[eqPG] / (gamma - 1.))
      + ((u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
         * 0.5);

  // cout <<"gamma="<<gamma<<"\n";
  return;
}

// ##################################################################
// ##################################################################

int eqns_mhd_ideal::UtoP(
    const pion_flt *u,
    pion_flt *p,
    const double MinTemp,  ///< minimum temperature/pressure allowed
    const double gamma)
{
  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX] / u[eqRHO];
  p[eqVY] = u[eqMMY] / u[eqRHO];
  p[eqVZ] = u[eqMMZ] / u[eqRHO];
  p[eqPG] =
      (gamma - 1)
      * (u[eqERG]
         - p[eqRO] * (p[eqVX] * p[eqVX] + p[eqVY] * p[eqVY] + p[eqVZ] * p[eqVZ])
               / 2.
         - (u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
               / 2.);
  p[eqBX] = u[eqBBX];
  p[eqBY] = u[eqBBY];
  p[eqBZ] = u[eqBBZ];

  int err = check_pressure(u, p, MinTemp, gamma);
  return err;
}

// ##################################################################
// ##################################################################

int eqns_mhd_ideal::check_pressure(
    const pion_flt *u,
    pion_flt *p,           ///< Primitive State Vector.
    const double MinTemp,  ///< minimum temperature/pressure allowed
    const double gamma)
{
#ifndef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  static long int ct_pg = 0;
#endif
  static long int ct_rho = 0;

  //
  // First check for negative density.
  // This is usually fatal to a simulation, so we quit here.
  //
  int err = 0;
  if (p[eqRO] <= 0.0) {
    spdlog::debug("u : {}", std::vector<pion_flt>(u, u + eq_nvar));
    spdlog::debug("p : {}", std::vector<pion_flt>(p, p + eq_nvar));
    spdlog::error("{}: {}", "Negative Density! Bugging out", p[eqRO]);
    exit(1);
    if (ct_rho < 1000) {
      ct_rho++;
      spdlog::debug("(eqns_mhd_ideal::check_pressure) negative density!");
      spdlog::debug("u : {}", std::vector<pion_flt>(u, u + eq_nvar));
      spdlog::debug("p : {}", std::vector<pion_flt>(p, p + eq_nvar));
    }
    // reset all variables because a negative density will change the
    // sign of all of the velocities
    p[eqRO] = BASE_RHO * eq_refvec[eqRO];
    p[eqVX] *= u[eqRHO] / p[eqRO];
    p[eqVY] *= u[eqRHO] / p[eqRO];
    p[eqVZ] *= u[eqRHO] / p[eqRO];
    p[eqPG] =
        (gamma - 1)
        * (u[eqERG]
           - p[eqRO]
                 * (p[eqVX] * p[eqVX] + p[eqVY] * p[eqVY] + p[eqVZ] * p[eqVZ])
                 / 2.
           - (u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
                 / 2.);
    err += 1;
  }

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // Check for negative pressure, fixing it if needed.
  // Then if there was no negative pressure, check for pressure
  // being too low, and fix that.
  //
  if (p[eqPG] <= 0.0) {
    // Set minimum temperature to be 10K
    // cout <<"UtoP() mhd set-neg-press-to-fixed-T.  P<0\n";
    spdlog::debug(
        "UtoP() mhd set-neg-press-to-fixed-T.  P = {}, Tmin = {}", p[eqPG],
        MinTemp);
    if (mp) {
      mp->Set_Temp(p, MinTemp, gamma);
    }
    else {
      // or set p=0.01*rho
      spdlog::debug(
          "UtoP() mhd fixing negative pressure from {}  to {}", p[eqPG],
          0.01 * p[eqRO]);
      p[eqPG] = 0.01 * p[eqRO];
    }
    err += 1;
  }
  else if (mp && (mp->Temperature(p, gamma) < MinTemp)) {
    // If we have microphysics, just set T=MinTemp
    spdlog::debug(
        "UtoP() mhd set-small-press-to-fixed-T.  T = {}, Tmin = {}",
        mp->Temperature(p, gamma), MinTemp);
    // rep.printVec("U",u,eq_nvar);
    // rep.printVec("p0",p,eq_nvar);
    mp->Set_Temp(p, MinTemp, gamma);
    // rep.printVec("p1",p,eq_nvar);
  }

#else   // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  if (p[eqPG] <= 0.) {
    if (ct_pg < 1000) {
      ct_pg++;
      spdlog::debug(
          "(eqns_mhd_ideal::check_pressure) -ve p_g= {} , correcting, count={}",
          [eqPG], ct_pg);
    }
    p[eqPG] = eq_refvec[eqPG] * 1.0e-6;
    err += 1;
  }
#endif  // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return err;
}

// ##################################################################
// ##################################################################

double eqns_mhd_ideal::chydro(
    const pion_flt *p,  ///< Pointer to primitive variables.
    const double gamma  ///< Gas constant gamma.
)
{
  return sqrt(gamma * p[eqPG] / p[eqRO]);
}

// ##################################################################
// ##################################################################

/* Calculate fast magnetic wavespeed, in direction we are looking in. */
double eqns_mhd_ideal::cfast(const pion_flt *p, const double gamma)
{
  // cout <<"cfast!\n";
  double ch = chydro(p, gamma);
  double temp1 =
      ch * ch
      + (p[eqBX] * p[eqBX] + p[eqBY] * p[eqBY] + p[eqBZ] * p[eqBZ]) / p[eqRO];
  double temp2 = 4. * ch * ch * p[eqBX] * p[eqBX] / p[eqRO];
  temp2        = max(MACHINEACCURACY, temp1 * temp1 - temp2);
  return (sqrt((temp1 + sqrt(temp2)) / 2.));
}

// ##################################################################
// ##################################################################

double eqns_mhd_ideal::cfast_components(
    const double cfRO,
    const double cfPG,
    const double cfBX,
    const double cfBY,
    const double cfBZ,
    const double g)
{
  double ch    = sqrt(g * cfPG / cfRO);
  double temp1 = ch * ch + (cfBX * cfBX + cfBY * cfBY + cfBZ * cfBZ) / cfRO;
  double temp2 = 4. * ch * ch * cfBX * cfBX / cfRO;
  temp2        = max(MACHINEACCURACY, temp1 * temp1 - temp2);
  return (sqrt((temp1 + sqrt(temp2)) / 2.));
}

// ##################################################################
// ##################################################################

// Calculate slow magnetic wavespeed.
double eqns_mhd_ideal::cslow(const pion_flt *p, const double gamma)
{
  double temp1;
  double temp2;
  double ch = chydro(p, gamma);
  temp1 =
      ch * ch
      + (p[eqBX] * p[eqBX] + p[eqBY] * p[eqBY] + p[eqBZ] * p[eqBZ]) / p[eqRO];
  temp2 = 4. * ch * ch * p[eqBX] * p[eqBX] / p[eqRO];
  temp2 = max(MACHINEACCURACY, temp1 * temp1 - temp2);
  temp2 = max(MACHINEACCURACY, temp1 - sqrt(temp2));
  return (sqrt(temp2 / 2.));
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::PUtoFlux(const pion_flt *p, const pion_flt *u, pion_flt *f)
{
  /// \section Equations
  /// The equations for the flux are in Falle, Komissarov, Joarder,
  /// (1998,MNRAS,297,265) Equation (2).
  ///
  double pm = (u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
              / 2.;  // Magnetic pressure.
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX] * p[eqVX] + p[eqPG] + pm - u[eqBBX] * u[eqBBX];
  f[eqMMY] = u[eqMMX] * p[eqVY] - u[eqBBX] * u[eqBBY];
  f[eqMMZ] = u[eqMMX] * p[eqVZ] - u[eqBBX] * u[eqBBZ];
  f[eqERG] =
      p[eqVX] * (u[eqERG] + p[eqPG] + pm)
      - u[eqBBX]
            * (p[eqVX] * u[eqBBX] + p[eqVY] * u[eqBBY] + p[eqVZ] * u[eqBBZ]);
  f[eqBBX] = 0.;
  f[eqBBY] = p[eqVX] * p[eqBY] - p[eqVY] * p[eqBX];
  f[eqBBZ] = p[eqVX] * p[eqBZ] - p[eqVZ] * p[eqBX];
  return;
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::UtoFlux(const pion_flt *u, pion_flt *f, const double gamma)
{
  double pm = (u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
              / 2.;  // Magnetic pressure.
  double pg =
      (gamma - 1.)
      * (u[eqERG]
         - (u[eqMMX] * u[eqMMX] + u[eqMMY] * u[eqMMY] + u[eqMMZ] * u[eqMMZ])
               / (2. * u[eqRHO])
         - pm);

  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX] * u[eqMMX] / u[eqRHO] + pg + pm - u[eqBBX] * u[eqBBX];
  f[eqMMY] = u[eqMMX] * u[eqMMY] / u[eqRHO] - u[eqBBX] * u[eqBBY];
  f[eqMMZ] = u[eqMMX] * u[eqMMZ] / u[eqRHO] - u[eqBBX] * u[eqBBZ];
  f[eqERG] =
      u[eqMMX] * (u[eqERG] + pg + pm) / u[eqRHO]
      - u[eqBBX]
            * (u[eqMMX] * u[eqBBX] + u[eqMMY] * u[eqBBY] + u[eqMMZ] * u[eqBBZ])
            / u[eqRHO];
  f[eqBBX] = 0.;
  f[eqBBY] = (u[eqMMX] * u[eqBBY] - u[eqMMY] * u[eqBBX]) / u[eqRHO];
  f[eqBBZ] = (u[eqMMX] * u[eqBBZ] - u[eqMMZ] * u[eqBBX]) / u[eqRHO];
  return;
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::PtoFlux(const pion_flt *p, pion_flt *f, const double gamma)
{
  pion_flt u[eq_nvar];
  eqns_mhd_ideal::PtoU(p, u, gamma);
  eqns_mhd_ideal::PUtoFlux(p, u, f);
  return;
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::rotate(
    pion_flt *vec,      ///< State vector
    enum axes initdir,  ///< Initial orientation.
    enum axes finaldir  ///< Final Orientation.
)
{
  /** \section Directions
   * This only rotates between the three positive directions, XP,YP,ZP,
   * so it never introduces sign changes to the elements, just reordering
   * of the vector components.
   * */
  if (initdir == finaldir) return;
  pion_flt v[eq_nvar];
  for (int i = 0; i < eq_nvar; i++)
    v[i] = vec[i];
  int offset = (static_cast<int>(finaldir - initdir + 3)) % 3;
  if (offset == 1) {
    v[eqVX] = vec[eqVY];
    v[eqVY] = vec[eqVZ];
    v[eqVZ] = vec[eqVX];
    v[eqBX] = vec[eqBY];
    v[eqBY] = vec[eqBZ];
    v[eqBZ] = vec[eqBX];
  }
  else if (offset == 2) {
    v[eqVX] = vec[eqVZ];
    v[eqVY] = vec[eqVX];
    v[eqVZ] = vec[eqVY];
    v[eqBX] = vec[eqBZ];
    v[eqBY] = vec[eqBX];
    v[eqBZ] = vec[eqBY];
  }
  else
    spdlog::error("{}: {}", "rotate function broken.", offset);
  for (int i = 0; i < eq_nvar; i++)
    vec[i] = v[i];
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::rotateXY(pion_flt *v, double theta)
{
  double ct = cos(theta);
  double st = sin(theta);

  pion_flt vx = v[eqVX] * ct - v[eqVY] * st;
  pion_flt vy = v[eqVX] * st + v[eqVY] * ct;
  v[eqVX]     = vx;
  v[eqVY]     = vy;

  vx      = v[eqBX] * ct - v[eqBY] * st;
  vy      = v[eqBX] * st + v[eqBY] * ct;
  v[eqBX] = vx;
  v[eqBY] = vy;
}

// ##################################################################
// ##################################################################

///  Returns Internal Energy (per unit mass, so 'Temperature'), given primitive
///  variable vector.
double eqns_mhd_ideal::eint(
    const pion_flt *p,  ///< Primitive State Vector.
    const double g      ///< gas EOS gamma.
)
{
  return p[eqPG] / (g - 1.) / p[eqRO];
}

// ##################################################################
// ##################################################################

/// Returns Total Energy (per unit volume), given primitive variable vector.
double eqns_mhd_ideal::Etot(
    const pion_flt *p,  ///< State Vector.
    const double g      ///< gas EOS gamma.
)
{
  return (
      (p[eqRO] * (p[eqVX] * p[eqVX] + p[eqVY] * p[eqVY] + p[eqVZ] * p[eqVZ])
       / 2.)
      + (p[eqPG] / (g - 1.))
      + ((p[eqBX] * p[eqBX] + p[eqBY] * p[eqBY] + p[eqBZ] * p[eqBZ]) / 2.));
}

// ##################################################################
// ##################################################################

/// Returns Total Pressure (per unit Volume), given primitive variable vector.
double eqns_mhd_ideal::Ptot(
    const pion_flt *p,  ///< Primitive State Vector.
    const double        ///< gas EOS gamma.
)
{
  return (
      p[eqPG]
      + 0.5 * (p[eqBX] * p[eqBX] + p[eqBY] * p[eqBY] + p[eqBZ] * p[eqBZ]));
}

// ##################################################################
// ##################################################################

/// Given a pressure ratio and initial density, calculate adiabatic final
/// density.
double eqns_mhd_ideal::AdiabaticRho(
    const double pr,  ///< New to Old pressure ratio
    const double ri,  ///< Old Density
    const double g    ///< gas EOS gamma.
)
{
  return (ri * exp(log(pr) / g));
}

// ##################################################################
// ##################################################################

void eqns_mhd_ideal::SetAvgState(
    const pion_flt *state,  ///< Mean Primitive var. state vector
    const double g          ///< Gas constant gamma.
)
{
  //
  // Set typical values to be used for reference:
  // (only really care about refB, refvel, refvec[PG,RO])
  //
  eq_refvec[eqRO] = state[eqRO];
  eq_refvec[eqPG] = state[eqPG];
  eq_refvec[eqVX] = state[eqVX];
  eq_refvec[eqVY] = state[eqVY];
  eq_refvec[eqVZ] = state[eqVZ];
  eq_refvec[eqBX] = state[eqBX];
  eq_refvec[eqBY] = state[eqBY];
  eq_refvec[eqBZ] = state[eqBZ];

  double angle = 0.0, refvel = 0.0, refB = 0.0;
  angle = eq_refvec[eqBY] * eq_refvec[eqBY] + eq_refvec[eqBX] * eq_refvec[eqBX];
  if (angle > 10. * MACHINEACCURACY) {
    angle = M_PI / 2. - asin(eq_refvec[eqBY] / sqrt(angle));
    if (eq_refvec[eqBX] < 0) angle = -angle;
    rotateXY(&eq_refvec[0], angle);
    refvel = cfast(&eq_refvec[0], eq_gamma);  // Fast speed
    rotateXY(&eq_refvec[0], -angle);
  }
  else
    refvel = maxspeed(&eq_refvec[0], eq_gamma);

  refB = sqrt(
      eq_refvec[eqBX] * eq_refvec[eqBX] + eq_refvec[eqBY] * eq_refvec[eqBY]
      + eq_refvec[eqBZ] * eq_refvec[eqBZ]);
  //
  // Now reset reference vector velocities and B-fields to be refvel,
  // and refB
  //
  eq_refvec[eqVX] = eq_refvec[eqVY] = eq_refvec[eqVZ] = 0.1 * refvel;
  eq_refvec[eqBX] = eq_refvec[eqBY] = eq_refvec[eqBZ] = refB;

#ifndef NDEBUG
  spdlog::debug("eq_refvec : {}", eq_refvec);
#endif
}

// ##################################################################
// ##################################################################

// ******************************************************************
// eqns_mhd_mixedGLM class, for the Dedner-GLM divergence cleaning method.
// ******************************************************************

// ##################################################################
// ##################################################################

eqns_mhd_mixedGLM::eqns_mhd_mixedGLM(int nv) : eqns_base(nv), eqns_mhd_ideal(nv)
{
  //  cout <<"eqns_mhd_mixedGLM constructor!\n";
  eqSI  = SI;
  eqPSI = PSI;

  // default values
  m_chyp_scale = 1.0;
  m_cpar_limit = 0.3;
}

// ##################################################################
// ##################################################################

eqns_mhd_mixedGLM::~eqns_mhd_mixedGLM() {}

// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::GLMsetPsiSpeed(const double ch, const double crel)
{
  GLM_chyp =
      m_chyp_scale * ch;  // hyperbolic wavespeed is equal to max. fast speed
  GLM_cr = crel;          // crel = 1/(cp^2/ch) has units of 1/length.
  return;
}

// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::PtoU(
    const pion_flt *P,  ///< pointer to Primitive variables.
    pion_flt *U,        ///< pointer to conserved variables.
    const double gamma  ///< Gas constant gamma.
)
{
  //  cout <<"glm ptou\n";
  U[eqPSI] = P[eqSI];
  eqns_mhd_ideal::PtoU(P, U, gamma);
  U[eqERG] += 0.5 * U[eqPSI] * U[eqPSI];
  return;
}

// ##################################################################
// ##################################################################

int eqns_mhd_mixedGLM::UtoP(
    const pion_flt *u,     ///< pointer to conserved variables.
    pion_flt *p,           ///< pointer to Primitive variables.
    const double MinTemp,  ///< minimum temperature/pressure allowed
    const double g         ///< Gas constant gamma.
)
{
  //  cout <<"glm utop\n";
  p[eqSI] = u[eqPSI];
  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX] / u[eqRHO];
  p[eqVY] = u[eqMMY] / u[eqRHO];
  p[eqVZ] = u[eqMMZ] / u[eqRHO];
  p[eqPG] =
      (g - 1.0)
      * (u[eqERG]
         - p[eqRO] * (p[eqVX] * p[eqVX] + p[eqVY] * p[eqVY] + p[eqVZ] * p[eqVZ])
               * 0.5
         - 0.5 * u[eqPSI] * u[eqPSI]
         - (u[eqBBX] * u[eqBBX] + u[eqBBY] * u[eqBBY] + u[eqBBZ] * u[eqBBZ])
               * 0.5);
  p[eqBX] = u[eqBBX];
  p[eqBY] = u[eqBBY];
  p[eqBZ] = u[eqBBZ];

  int err = check_pressure(u, p, MinTemp, g);
  return err;
}

// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::GLMsource(
    pion_flt *psivar,  ///< Primitive Psi variable.
    const double delt  ///< timestep
)
{
#ifndef NDEBUG
  // TESTING ONLY
  static double temp = 0.0;
  if (!pconst.equalD(delt * GLM_chyp * GLM_cr, temp)) {
    spdlog::info(
        "glmsource: exp factor = {:12.4e}, lim = {}, scale = {}",
        delt * GLM_chyp * GLM_cr, m_cpar_limit, m_chyp_scale);
  }
  temp = delt * GLM_chyp * GLM_cr;
#endif

  // need to limit this because the exponent doubles for each coarser
  // level that we reach
  *psivar *= exp(-min(delt * GLM_chyp * GLM_cr, m_cpar_limit));
  return;
}

// ##################################################################
// ##################################################################
