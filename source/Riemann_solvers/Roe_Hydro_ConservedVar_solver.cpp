///
/// \file Roe_Hydro_ConservedVar_solver.cc
/// \author Jonathan Mackey
///
/// This is a linearised Riemann solver for the Euler Equations.
/// It is the Roe solver in conserved variables, and it can use
/// the H-correction to add multi-dimensional viscosity to the
/// fluxes.  The symmetric solver is recommended, the one-sided is
/// not.  Both will return a flux in the first 5 state variables (i.e.
/// the hydro vars) and also a state vector for the interface which
/// can be used for adding viscosity.
///
/// References:
///  - Toro (1999) Riemann Solvers Textbook, Chapter 11.2.2, pp.350-353.
///  - Sanders, Morano, Druguet, (1998, JCP, 145, 511).
///
/// History:
/// - 2010-12.22 JM: Moved functions from flux_hydro_adiabatic.h/.cc
///   and generated new class.  Split everything into functions.
///   Defined a bunch of private class variables.
///
/// - 2010.12.27 JM: Added check in Hcorr for non-zero etamax when
///   AVtype!=3 (only when code run with TESTING).
///
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer
///    fluxes are dealt with by the flux-solver classes.
///    For this to work I had to explicitly call the Euler Eqns class for PtoU()
///    and PUtoFlux(); otherwise the flux-solver class (with tracers) would be
///    used.
/// - 2015.01.14 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2016.05.21 JM: removed H-correction ifdefs (it should be always enabled).

#include "Riemann_solvers/Roe_Hydro_ConservedVar_solver.h"
#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"

using namespace std;

// ##################################################################
// ##################################################################

Riemann_Roe_Hydro_CV::Riemann_Roe_Hydro_CV(
    const int nv,   ///< Length of State Vectors, nvar
    const double g  ///< Gamma for state vector.
    ) :
    eqns_base(nv),
    eqns_Euler(nv), rs_nvar(5)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::Riemann_Roe_Hydro_CV ...starting.\n";
#endif  // FUNCTION_ID
  //
  // eq_gamma, eq_nvar are defined in eqns_base class
  //
  eq_gamma = g;

  //
  // Set the enthalpy index in primitive vector to be the same as the
  // pressure index.
  //
  eqHH = eqPG;

  //
  // Allocate memory for evalues, evectors, mean-state,
  // difference-state, left and right conserved var states.
  //
  RCV_meanp    = mem.myalloc(RCV_meanp, rs_nvar);
  RCV_ul       = mem.myalloc(RCV_ul, rs_nvar);
  RCV_ur       = mem.myalloc(RCV_ur, rs_nvar);
  RCV_eval     = mem.myalloc(RCV_eval, rs_nvar);
  RCV_strength = mem.myalloc(RCV_strength, rs_nvar);
  RCV_udiff    = mem.myalloc(RCV_udiff, rs_nvar);
  RCV_evec     = mem.myalloc(RCV_evec, rs_nvar);
  for (int v = 0; v < 5; v++)
    RCV_evec[v] = mem.myalloc(RCV_evec[v], rs_nvar);

  RCV_HC_etamax = 0.0;

  RCV_v2_mean = 0.0;
  RCV_a_mean  = 0.0;

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::Riemann_Roe_Hydro_CV ...returning.\n";
#endif  // FUNCTION_ID
}

// ##################################################################
// ##################################################################

Riemann_Roe_Hydro_CV::~Riemann_Roe_Hydro_CV()
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::~Riemann_Roe_Hydro_CV ...starting.\n";
#endif  // FUNCTION_ID

  RCV_meanp    = mem.myfree(RCV_meanp);
  RCV_ul       = mem.myfree(RCV_ul);
  RCV_ur       = mem.myfree(RCV_ur);
  RCV_eval     = mem.myfree(RCV_eval);
  RCV_strength = mem.myfree(RCV_strength);
  RCV_udiff    = mem.myfree(RCV_udiff);
  for (int v = 0; v < rs_nvar; v++)
    RCV_evec[v] = mem.myfree(RCV_evec[v]);
  RCV_evec = mem.myfree(RCV_evec);

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::~Riemann_Roe_Hydro_CV ...returning.\n";
#endif  // FUNCTION_ID
}

// ##################################################################
// ##################################################################

int Riemann_Roe_Hydro_CV::Roe_flux_solver_symmetric(
    const pion_flt *left,
    const pion_flt *right,
    const double g,
    const double hc_eta,
    pion_flt *out_pstar,
    pion_flt *out_flux)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::Roe_flux_solver_symmetric ...starting.\n";
#endif  // FUNCTION_ID
  //
  // Roe's linear flux solver, from Toro 11.2.2 pp.350-353.  This is
  // basically the same as the asymmetric one,
  // Roe_flux_solver_onesided(), except that we always take the mean
  // of all waves rather than just traversing from left to right.  In
  // theory this should be give a perfectly symmetric solution when
  // the initial conditions are symmetric.  It is much
  // better than the asymmetric version -- especially with the
  // H-correction the quality of the solution is *much* better.
  //
  eq_gamma = g;

  //
  // If using the H-correction, this will be non-zero, and the
  // eigenvalues will be modified according to the prescription in
  // Sanders, Morano, Druguet, (1998, JCP, 145, 511, eq.10).  If this
  // parameter is zero, then it has no effect on the solution.
  //
  RCV_HC_etamax = hc_eta;

  //
  // (1) Get Roe-average values for rho,vx,vy,vz,H,a
  // Toro eq. 11.60
  //
  // This assigns a mean state to RCV_meanp and also sets the vars
  // RCV_a_mean and RCV_v2_mean
  //
  set_Roe_mean_state(left, right);

  //
  // Eigenvalues: sets the values of eval[] based on meanp,a_mean.
  //
  set_eigenvalues();

  //
  // calculate the eigenvectors and wavestrengths.
  //
  //
  // Eigenvectors (eq.11.59 Toro)
  // Based on RCV_meanp, RCV_a_mean, RCV_v2_mean.
  //
  set_eigenvectors();

  //
  // Set left and right states in conservative variables, and then get
  // the difference vector.  Store results in ul[], ur[], udiff[].
  //
  set_ul_ur_udiff(left, right);

  //
  // this uses udiff, RCV_meanp, eq_gamma, RCV_a_mean to
  // calculate the wavestrengths according to Toro (1999).
  // Strengths are set in the class variable str[].
  //
  set_wave_strengths();

  //
  // Now get the flux by stepping across waves
  //
  // Flux = 0.5(F(left)+F(right)
  //            -sum_{waves}[strength_i*|evalue_i|*right_evec_i])
  //
  // This overwrites ul[], so don't use it past this point!
  //
  calculate_symmetric_flux(out_flux);

  set_pstar_from_meanp(out_pstar);

#ifdef DEBUG
  //  if (fabs(out_flux[0])>1.0e-50) {
  //    rep.printVec("RCV_meanp",RCV_meanp,5);
  //    cout <<"RCV_a_mean="<<RCV_a_mean<<", RCV_v2_mean="<<RCV_v2_mean<<"\n";
  //    rep.printVec("RCV_eval",RCV_eval,5);
  //    rep.printVec("RCV_udiff",RCV_udiff,5);
  //    rep.printVec("RCV_strength",RCV_strength,5);
  //    rep.printVec("out_flux",out_flux,5);
  //    rep.printVec("out_pstar",out_pstar,5);
  //  }

  if (fabs(out_flux[1]) > 1.0e50) {
    cout << "Very big Energy Flux!\n";
    rep.printVec("\tleft", left, 5);
    rep.printVec("\tright", right, 5);
    rep.printVec("\tRCV_meanp", RCV_meanp, 5);
    cout << "\tRCV_a_mean=" << RCV_a_mean << ",  RCV_v2_mean=" << RCV_v2_mean
         << "\n";
    rep.printVec("\tRCV_eval", RCV_eval, 5);
    rep.printVec("\tRCV_udiff", RCV_udiff, 5);
    rep.printVec("\tRCV_strength", RCV_strength, 5);
    rep.printVec("\tout_flux", out_flux, 5);
    rep.printVec("\tout_pstar", out_pstar, 5);
  }
#endif  // DEBUG

  //
  // Note that out_pstar[] and out_flux[] only have the first 5
  // physical elements set.  The tracer values need to be set
  // somewhere else (i.e. the calling function).
  //

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::Roe_flux_solver_symmetric ...returning.\n";
#endif  // FUNCTION_ID
  return 0;
}

// -------------------------------------------------------------------
// ------------------- END OF PUBLIC FUNCTIONS -----------------------
// -------------------------------------------------------------------

// ##################################################################
// ##################################################################

// -------------------------------------------------------------------
// --------------------- PRIVATE FUNCTIONS ---------------------------
// -------------------------------------------------------------------

int Riemann_Roe_Hydro_CV::test_left_right_equality(
    const pion_flt *left, const pion_flt *right)
{
  //
  // Return true if the left and right input states are almost
  // identical.
  //
  double diff = 0.0;
  diff += fabs(right[eqRO] - left[eqRO]) / (right[eqRO] + left[eqRO]);
  diff += fabs(right[eqPG] - left[eqPG]) / (right[eqPG] + left[eqPG]);
  diff += fabs(right[eqVX] - left[eqVX])
          / (fabs(right[eqVX]) + fabs(left[eqVX]) + SMALLVALUE);
  diff += fabs(right[eqVY] - left[eqVY])
          / (fabs(right[eqVY]) + fabs(left[eqVY]) + SMALLVALUE);
  diff += fabs(right[eqVZ] - left[eqVZ])
          / (fabs(right[eqVZ]) + fabs(left[eqVZ]) + SMALLVALUE);

  if (diff < 1.e-6) {
#ifdef RoeHYDRO_TESTING
    cout << "same states...\n";
#endif  // RoeHYDRO_TESTING
    return 1;
  }

  return 0;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_Roe_mean_state(
    const pion_flt *left, const pion_flt *right)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_Roe_mean_state ...starting.\n";
#endif  // FUNCTION_ID
  //
  // (1) Get Roe-average values for vx,vy,vz,H,a
  // Toro eq. 11.60
  //
  double rl = sqrt(left[eqRO]), rr = sqrt(right[eqRO]),
         lH = Enthalpy(left, eq_gamma), rH = Enthalpy(right, eq_gamma),
         denom = 1.0 / (rl + rr);

  RCV_a_mean  = 0.0;
  RCV_v2_mean = 0.0;

  RCV_meanp[eqRO] = rl * rr;
  RCV_meanp[eqVX] = (rl * left[eqVX] + rr * right[eqVX]) * denom;
  RCV_meanp[eqVY] = (rl * left[eqVY] + rr * right[eqVY]) * denom;
  RCV_meanp[eqVZ] = (rl * left[eqVZ] + rr * right[eqVZ]) * denom;
  RCV_meanp[eqHH] = (rl * lH + rr * rH) * denom;

  RCV_v2_mean = RCV_meanp[eqVX] * RCV_meanp[eqVX]
                + RCV_meanp[eqVY] * RCV_meanp[eqVY]
                + RCV_meanp[eqVZ] * RCV_meanp[eqVZ];
  // make sure sound speed is positive, and at least 10^{-6} of the
  // flow velocities
  RCV_a_mean = sqrt(
      (eq_gamma - 1.0)
      * max(RCV_meanp[eqHH] - 0.5 * RCV_v2_mean, 1.0e-12 * RCV_v2_mean));

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_Roe_mean_state ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_eigenvalues()
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_eigenvalues ...starting.\n";
#endif  // FUNCTION_ID

  RCV_eval[0] = RCV_meanp[eqVX] - RCV_a_mean;
  RCV_eval[1] = RCV_eval[2] = RCV_eval[3] = RCV_meanp[eqVX];
  RCV_eval[4]                             = RCV_meanp[eqVX] + RCV_a_mean;

  //
  // Modify the eigenvalues by the H-correction eta value.  Note that
  // HC_etamax is set to zero in the simulation initialisation, and it
  // is only changed if AVtype==3, so this code has no effect if we
  // are not using the H-correction.  Unless of course the eigenvalue
  // is _very_ close to zero and the eval changes within the machine
  // precision.  This shouldn't ever have a significant effect on
  // anything!
  //
  for (int v = 0; v < rs_nvar; v++) {
    //
    // if eval[v] <0  then set it to min(eval[v],-eta)
    // else                set it to max(eval[v], eta)
    //
    if (RCV_eval[v] < 0.0) {
      RCV_eval[v] = min(RCV_eval[v], -RCV_HC_etamax);
    }
    else {
      RCV_eval[v] = max(RCV_eval[v], RCV_HC_etamax);
    }
  }

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_eigenvalues ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_eigenvectors()
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_eigenvectors ...starting.\n";
#endif  // FUNCTION_ID
  //
  // Eigenvectors (eq.11.59 Toro)
  //
  RCV_evec[0][eqRHO] = 1.0;
  RCV_evec[0][eqMMX] = RCV_meanp[eqVX] - RCV_a_mean;
  RCV_evec[0][eqMMY] = RCV_meanp[eqVY];
  RCV_evec[0][eqMMZ] = RCV_meanp[eqVZ];
  RCV_evec[0][eqERG] = RCV_meanp[eqHH] - RCV_meanp[eqVX] * RCV_a_mean;

  RCV_evec[1][eqRHO] = 1.0;
  RCV_evec[1][eqMMX] = RCV_meanp[eqVX];
  RCV_evec[1][eqMMY] = RCV_meanp[eqVY];
  RCV_evec[1][eqMMZ] = RCV_meanp[eqVZ];
  RCV_evec[1][eqERG] = 0.5 * RCV_v2_mean;

  RCV_evec[2][eqRHO] = 0.0;
  RCV_evec[2][eqMMX] = 0.0;
  RCV_evec[2][eqMMY] = 1.0;
  RCV_evec[2][eqMMZ] = 0.0;
  RCV_evec[2][eqERG] = RCV_meanp[eqVY];

  RCV_evec[3][eqRHO] = 0.0;
  RCV_evec[3][eqMMX] = 0.0;
  RCV_evec[3][eqMMY] = 0.0;
  RCV_evec[3][eqMMZ] = 1.0;
  RCV_evec[3][eqERG] = RCV_meanp[eqVZ];

  RCV_evec[4][eqRHO] = 1.0;
  RCV_evec[4][eqMMX] = RCV_meanp[eqVX] + RCV_a_mean;
  ;
  RCV_evec[4][eqMMY] = RCV_meanp[eqVY];
  RCV_evec[4][eqMMZ] = RCV_meanp[eqVZ];
  RCV_evec[4][eqERG] = RCV_meanp[eqHH] + RCV_meanp[eqVX] * RCV_a_mean;
  ;

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_eigenvectors ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_ul_ur_udiff(
    const pion_flt *left, const pion_flt *right)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_ul_ur_udiff ...starting.\n";
#endif  // FUNCTION_ID

  eqns_Euler::PtoU(left, RCV_ul, eq_gamma);
  eqns_Euler::PtoU(right, RCV_ur, eq_gamma);
  //
  // Wave Strengths Toro eq.11.68/.69/.70
  // We need the difference states for this:
  //
  //  rep.printVec("ud",udiff,rs_nvar);
  for (int v = 0; v < rs_nvar; v++) {
    if (pconst.equalD(RCV_ur[v], RCV_ul[v]))
      RCV_udiff[v] = 0.0;
    else
      RCV_udiff[v] = RCV_ur[v] - RCV_ul[v];
  }
  // rep.printVec("ul",RCV_ul,rs_nvar);
  // rep.printVec("ur",RCV_ur,rs_nvar);
  // rep.printVec("ud",RCV_udiff,rs_nvar);

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_ul_ur_udiff ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_wave_strengths()
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_wave_strengths ...starting.\n";
#endif  // FUNCTION_ID
  //
  // calculate the wavestrengths according to Toro (1999,
  // eq.11.68/.69/.70).
  //
  // Uses the difference state udiff[] and the mean state RCV_meanp[]
  //

  RCV_strength[2] = RCV_udiff[eqMMY] - RCV_meanp[eqVY] * RCV_udiff[eqRO];
  RCV_strength[3] = RCV_udiff[eqMMZ] - RCV_meanp[eqVZ] * RCV_udiff[eqRO];

  double u5bar = RCV_udiff[eqERG] - RCV_strength[2] * RCV_meanp[eqVY]
                 - RCV_strength[3] * RCV_meanp[eqVZ];

  RCV_strength[1] =
      (RCV_udiff[eqRHO] * (RCV_meanp[eqHH] - RCV_meanp[eqVX] * RCV_meanp[eqVX])
       + RCV_meanp[eqVX] * RCV_udiff[eqMMX] - u5bar)
      * (eq_gamma - 1.0) / RCV_a_mean / RCV_a_mean;
  RCV_strength[0] = 0.5
                    * (RCV_udiff[eqRHO] * (RCV_meanp[eqVX] + RCV_a_mean)
                       - RCV_udiff[eqMMX] - RCV_a_mean * RCV_strength[1])
                    / RCV_a_mean;
  RCV_strength[4] = RCV_udiff[eqRHO] - RCV_strength[0] - RCV_strength[1];

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_wave_strengths ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::calculate_symmetric_flux(pion_flt *out_flux)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::calculate_symmetric_flux ...starting.\n";
#endif  // FUNCTION_ID
  //
  // Now get the flux by stepping across waves We also get the starred
  // state in pstar[] using the mean state vector for pstar, since
  // it has positive-definite pressure and density (this is actually a
  // better idea than using the actual pstar obtained from ustar).
  //
  // Flux = 0.5(F(left)+F(right)
  //            -sum_{waves}[strength_i*|evalue_i|*right_evec_i])
  //
  // We use the left state vector ul[] as a temp array once it has
  // made its contribution to the total flux.
  //
  eqns_Euler::UtoFlux(RCV_ul, out_flux, eq_gamma);
  eqns_Euler::UtoFlux(RCV_ur, RCV_ul, eq_gamma);
  for (int v = 0; v < rs_nvar; v++)
    out_flux[v] += RCV_ul[v];

  //
  // Now add the contributions of the states between all the waves:
  //
  for (int iwave = 0; iwave < 5; iwave++) {
    out_flux[eqRHO] -=
        RCV_strength[iwave] * fabs(RCV_eval[iwave]) * RCV_evec[iwave][eqRHO];
    out_flux[eqMMX] -=
        RCV_strength[iwave] * fabs(RCV_eval[iwave]) * RCV_evec[iwave][eqMMX];
    out_flux[eqMMY] -=
        RCV_strength[iwave] * fabs(RCV_eval[iwave]) * RCV_evec[iwave][eqMMY];
    out_flux[eqMMZ] -=
        RCV_strength[iwave] * fabs(RCV_eval[iwave]) * RCV_evec[iwave][eqMMZ];
    out_flux[eqERG] -=
        RCV_strength[iwave] * fabs(RCV_eval[iwave]) * RCV_evec[iwave][eqERG];
  }

  //
  // Divide by 2 to get the correct flux
  //
  for (int v = 0; v < rs_nvar; v++)
    out_flux[v] *= 0.5;

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::calculate_symmetric_flux ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::set_pstar_from_meanp(pion_flt *out_pstar)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_pstar_from_meanp ...starting.\n";
#endif  // FUNCTION_ID
  //
  // We also get the starred state in pstar[] using the mean state
  // vector for pstar, since it has positive-definite pressure and
  // density (this is actually a better idea than using the actual
  // pstar obtained from ustar).
  //
  // Convert the enthalpy variable to pressure for the data
  // assignment.
  //
  for (int v = 0; v < rs_nvar; v++)
    out_pstar[v] = RCV_meanp[v];
  out_pstar[eqPG] = RCV_meanp[eqRO] * RCV_a_mean * RCV_a_mean / eq_gamma;

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::set_pstar_from_meanp ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################

void Riemann_Roe_Hydro_CV::calculate_asymmetric_flux(
    const pion_flt *left, const pion_flt *right, pion_flt *out_flux)
{
#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::calculate_asymmetric_flux ...starting.\n";
#endif  // FUNCTION_ID
  //
  // We need to step across one or more waves:
  // First get the left state flux:
  //
  // rep.printVec("left  flux:",out_flux,5);
  pion_flt utemp[rs_nvar];
  eqns_Euler::PtoU(left, utemp, eq_gamma);
  eqns_Euler::PUtoFlux(left, utemp, out_flux);
  // eqns_Euler::PtoFlux(left,out_flux,eq_gamma);

  //
  // Next add the contribution from all the waves with negative
  // wavespeeds:
  //
  int i = 0;
  while ((i < 5) && (RCV_eval[i] <= 0.0)) {
    out_flux[eqRHO] += RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqRHO];
    out_flux[eqMMX] += RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMX];
    out_flux[eqMMY] += RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMY];
    out_flux[eqMMZ] += RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMZ];
    out_flux[eqERG] += RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqERG];
    i++;
  }

#ifdef RoeHYDRO_TESTING
  //
  // We should have calculated the flux across the boundary.  Check
  // this by going back from the right to the left, and making sure
  // we get the same answer!
  //
  pion_flt ftemp[rs_nvar];
  // double utemp[rs_nvar];
  eqns_Euler::PtoU(right, utemp, eq_gamma);
  eqns_Euler::PUtoFlux(right, utemp, ftemp);
  // eqns_Euler::PtoFlux(right,ftemp,eq_gamma);

  i = 4;
  while ((i >= 0) && (RCV_eval[i] > 0.0)) {
    ftemp[eqRHO] -= RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqRHO];
    ftemp[eqMMX] -= RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMX];
    ftemp[eqMMY] -= RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMY];
    ftemp[eqMMZ] -= RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqMMZ];
    ftemp[eqERG] -= RCV_strength[i] * RCV_eval[i] * RCV_evec[i][eqERG];
    i--;
  }

  //
  // Now compare the difference:
  //
  double diff = 0.0;
  for (int v = 0; v < rs_nvar; v++) {
    // cout <<"flux["<<v<<"]: left="<<out_flux[v]<<" and
    // right="<<ftemp[v]<<"\n";
    diff += (out_flux[v] - ftemp[v])
            / (fabs(out_flux[v]) + fabs(ftemp[v]) + TINYVALUE)
  }
  if (diff > 1e-3) {
    cout << "*** FLUX CALCULATION ERROR IN "
            "flux_solver_hydro_adi::Roe_flux_solver(): diff = "
         << diff << "\n";
    rep.printVec("left  flux:", out_flux, 5);
    rep.printVec("right flux:", ftemp, 5);
    rep.printVec("e-values", RCV_eval, 5);
    rep.printVec("strengths:", RCV_strength, 5);
    rep.printVec(" left", left, 5);
    rep.printVec("right", right, 5);
    rep.printVec("meanp", RCV_meanp, 5);
    rep.printVec("udiff", RCV_udiff, 5);
    for (int i = 0; i < 5; ++i) {
      cout << "rightevec[" << i << "] = [ ";
      for (int j = 0; j < 5; j++) {
        cout.width(9);
        cout << RCV_evec[i][j] << ", ";
      }
      cout << "]"
           << "\n";
    }
    cout << "*** FLUX CALCULATION ERROR IN "
            "flux_solver_hydro_adi::Roe_flux_solver()\n\n";
  }  // if large difference

#endif  // RoeHYDRO_TESTING

#ifdef FUNCTION_ID
  cout << "Riemann_Roe_Hydro_CV::acalculate_symmetric_flux ...returning.\n";
#endif  // FUNCTION_ID
  return;
}

// ##################################################################
// ##################################################################
