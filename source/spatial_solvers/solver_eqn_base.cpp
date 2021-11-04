/// \file solver_eqn_base.cc
///
/// \brief Class definition for various solvers implemented in my FV grid-code.
/// \author Jonathan Mackey
///
/// Modifications:
///  - 2007-07-10 Part Way through writing it.
///  - 2007-07-11 Still writing it and working out class heirarchy.
///  - 2007-07-12 Got the class heirarchy working (I think...).
///  - 2007-07-13 New Class structure implemented.
///  - 2007-07-16 Reverted to less complicated class hierarchy.  Runs fast.
///  - 2007-07-23 Started to add passive tracer variables
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08-01 cylindrical coordinates support.
///  - 2007-08-08 cylindrical coordinates hd/i-mhd/glm-mhd working.
///  - 2007-11-30 Worked on Toth's Field-CD method for keeping divB=0.
///  - 2009-10-20 New structure built on equations and flux classes...
///  - 2009-10-24 Cut out all the old solver classes.
///  - 2010.09.30 JM: Worked on Lapidus AV (added set_div_v() function for
///  cells)
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.  Added Pre- and
///   Post-flux viscosity functions for the H-correction and Lapidus
///   viscosity functions.
///   Added H-correction functions for calculating eta and getting
///   the maximum value of eta on a given H-stencil.
///   Made InterCellFlux general for all classes.
/// - 2010.11.19 JM: Debugged.
/// - 2010.12.08 JM: Got the H-correction to work for 1D grids.
/// - 2010.12.23 JM: Added new set_interface_tracer_flux(lp,rp,flux);
///   function.  Removed riemann_base calls.  Added pstar[] array
///   to intercell_flux() function since it is no longer inherited.
/// - 2010.12.27 JM: modified args to post_calc_viscous_terms to pass
///   in left, right, pstar (since these are no longer class vars).
/// - 2011.01.03 JM: Moved preprocess_data() and calc_Hcorrection from
///   gridMethods.cc
/// - 2011.02.25 JM: removed HCORR ifdef around new code
/// - 2011.04.06 JM: Added thermal-conduction to preprocess_data.
/// - 2012.08.05 JM: Added spacers between functions for readability.
///    Thermal Conduction module is unusable now, needs re-writing.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2016.05.21 JM: Tidied up H-correction terms.
/// - 2016.08.25 JM: Changed H-correction loop to start with
///    FirstPt_All() instead of FirstPt().  Now serial/parallel code
///    produces identical results for the 3D blastwave test.
/// - 2018.01.24 JM: worked on making SimPM non-global
/// - 2018.04.14 JM: Moved flux solver to FV_solver
/// - 2018.08.03 JM: Added Berger & Colella (1989) flux correction
///    algorithm to PostProcess_dU()

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "solver_eqn_base.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"

#ifdef PION_OMP
#include <omp.h>
#endif

using namespace std;

// ##################################################################
// ##################################################################

FV_solver_base::FV_solver_base(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_gndim(nd), FV_cfl(cflno), FV_etav(avcoeff), FV_etaB(avcoeff), FV_ntr(ntr)
{
  eq_gamma = gam;
  eqTR     = 0;
  if (FV_ntr > 0) {
    eqTR = mem.myalloc(eqTR, FV_ntr);
    for (int i = 0; i < FV_ntr; i++)
      eqTR[i] = eq_nvar - FV_ntr + i;
  }
  HC_etamax = 0.0;
  return;
}

FV_solver_base::~FV_solver_base()
{
  if (FV_ntr > 0) {
    eqTR = mem.myfree(eqTR);
  }
  return;
}

// ##################################################################
// ##################################################################

int FV_solver_base::get_LaxFriedrichs_flux(
    const pion_flt *l,
    const pion_flt *r,
    pion_flt *f,
    const double dx,  ///< cell size dx
    const double)
{
  //
  // This has to be in this solver because we need dx,dt,ndim
  //
  pion_flt u1[eq_nvar], u2[eq_nvar], f1[eq_nvar], f2[eq_nvar];
  PtoU(l, u1, eq_gamma);
  PtoU(r, u2, eq_gamma);
  UtoFlux(u1, f1, eq_gamma);
  UtoFlux(u2, f2, eq_gamma);
  // Then get inter-cell flux from this.
  for (int v = 0; v < eq_nvar; v++)
    f[v] = 0.5 * (f1[v] + f2[v] + dx / FV_dt * (u1[v] - u2[v]) / FV_gndim);

  //
  // Calculate tracer flux based on whether flow is to left or right.
  //
  if (FV_ntr > 0) {
    if (f[eqRHO] >= 0.) {
      for (int t = 0; t < FV_ntr; t++)
        f[eqTR[t]] = l[eqTR[t]] * f[eqRHO];
    }
    else {
      for (int t = 0; t < FV_ntr; t++)
        f[eqTR[t]] = r[eqTR[t]] * f[eqRHO];
    }
  }
  return 0;
}

// ##################################################################
// ##################################################################

///
/// Calculate Flux between a left and right state.
///
int FV_solver_base::InterCellFlux(
    class SimParams &par,  ///< simulation parameters
    class GridBaseClass *grid,
    class cell *Cl,  ///< Left state cell pointer
    class cell *Cr,  ///< Right state cell pointer
    pion_flt *lp,    ///< Left Primitive State Vector.
    pion_flt *rp,    ///< Right Primitive State Vector.
    pion_flt *f,     ///< Flux Vector. (written to).
    const double g,  ///< gas EOS gamma.
    const double dx  ///< Cell size dx.
)
{
  // cout <<"FV_solver_base::InterCellFlux() gamma="<<eq_gamma<<" and passed
  // in was g="<<g<<"\n";
  eq_gamma = g;
  pion_flt pstar[eq_nvar];

  //
  // Pre-calcualate anything needed for the viscosity (H-correction).
  // Data is stored in each cell, so no values returned.
  //
  pre_calc_viscous_terms(grid, Cl, Cr, par.artviscosity);

  //
  // Get the flux from the flux solver:
  //
  int err =
      inviscid_flux(par, grid, dx, Cl, Cr, lp, rp, f, pstar, par.solverType, g);

#ifdef DEBUG
  if (fabs(f[0]) > 1.0e-50) {
    rep.printVec("flux=", f, eq_nvar);
    rep.printVec("left=", lp, eq_nvar);
    rep.printVec("rght=", rp, eq_nvar);
    rep.printVec("pstr=", pstar, eq_nvar);
  }
#endif

  //
  // Post-calculate anthing needed for the viscosity: calls the FKJ98
  // viscosity function which acts after the flux has been calculated
  //
  post_calc_viscous_terms(Cl, Cr, lp, rp, pstar, f, par.artviscosity);

  //
  // Calculate tracer flux based on whether flow is to left or right.
  //
  set_interface_tracer_flux(lp, rp, f);
  // rep.printVec("left=",lp,eq_nvar);
  // rep.printVec("rght=",rp,eq_nvar);
  // rep.printVec("flux",f,eq_nvar);

  return err;
}

// ##################################################################
// ##################################################################

void FV_solver_base::pre_calc_viscous_terms(
    class GridBaseClass *grid,
    const cell *cl,    ///< left-of-interface cell
    const cell *cr,    ///< right-of-interface cell
    const int av_flag  ///< what kind of AV?
)
{
  //
  // Only the H-correction acts before the actual flux calculation:
  //
  switch (av_flag) {
    case AV_HCORRECTION:  // H-correction
    case AV_HCORR_FKJ98:  // H-correction +FKJ98 viscosity
      FV_solver_base::HC_etamax = select_Hcorr_eta(cl, cr, grid);
      break;
    default:
      // Just silently continue if not doing the H-correction.
      break;
  }

  return;
}

// ##################################################################
// ##################################################################

void FV_solver_base::post_calc_viscous_terms(
    const cell *cl,  ///< left-of-interface cell
    const cell *cr,  ///< right-of-interface cell
    const pion_flt *Pl,
    const pion_flt *Pr,
    const pion_flt *Pstar,
    pion_flt *flux,    ///< flux vector
    const int av_flag  ///< what kind of AV?
)
{
  //  cout <<"etav="<<FV_etav<<"\t";  rep.printVec("flux",flux,eq_nvar);
  //  rep.printVec("flux",flux,eq_nvar);

  int err = 0;
  switch (av_flag) {
    case AV_FKJ98_1D:     // FKJ98 Viscosity
    case AV_HCORR_FKJ98:  // FKJ98+HCORR
      // cout <<"FKJ98 AV being applied!\n";
      err += AVFalle(Pl, Pr, Pstar, flux, FV_etav, eq_gamma);
      break;

    default:
      // Silently continue if av_flag==0 or 3.
      break;
  }

  //  rep.printVec("flux",flux,eq_nvar);
  return;
}

// ##################################################################
// ##################################################################

void FV_solver_base::set_interface_tracer_flux(
    const pion_flt *left,   // prim.var.
    const pion_flt *right,  // prim.var.
    pion_flt *flux)
{
  pion_flt corrector[eq_nvar];
  for (int v = 0; v < eq_nvar; v++)
    corrector[v] = 1.0;

#ifdef FUNCTION_ID
  cout << "FV_solver_base::set_interface_tracer_flux ...starting.\n";
#endif  // FUNCTION_ID
        //
        // Calculate tracer flux here -- if mass flux is positive then
        // contact is at x>0 and we advect the left state tracer across
        // the boundary.  Otherwise we advect the right state to the left.
        //
        // N.B. Without introducing minimum density and velocity scales, it
        // is impossible to avoid introducing some asymmetry here.
#ifdef TEST_SYMMETRY
  if (FV_ntr > 0) {
    if (flux[eqRHO] > 1.0e-28)
      for (int t = 0; t < FV_ntr; t++)
        flux[eqTR[t]] = left[eqTR[t]] * flux[eqRHO];
    else if (flux[eqRHO] < -1.0e-28)
      for (int t = 0; t < FV_ntr; t++)
        flux[eqTR[t]] = right[eqTR[t]] * flux[eqRHO];
    else
      for (int t = 0; t < FV_ntr; t++)
        flux[eqTR[t]] = 0.0;
  }
#else
  if (FV_ntr > 0) {
#ifdef TEST_INF
    if (!isfinite(flux[eqRHO])) {
      cout << "FV_solver_base::set_interface_tracer_flux: ";
      rep.printVec("inf flux", flux, eq_nvar);
    }
#endif
    if (flux[eqRHO] > 0.0) {
      if (mp) mp->sCMA(corrector, left);
      for (int t = 0; t < FV_ntr; t++) {
        flux[eqTR[t]] = left[eqTR[t]] * flux[eqRHO] * corrector[eqTR[t]];
      }
    }
    else if (flux[eqRHO] < 0.0) {
      if (mp) mp->sCMA(corrector, right);
      for (int t = 0; t < FV_ntr; t++) {
        flux[eqTR[t]] = right[eqTR[t]] * flux[eqRHO] * corrector[eqTR[t]];
      }
    }
    else {
      for (int t = 0; t < FV_ntr; t++)
        flux[eqTR[t]] = 0.0;
    }
  }
#endif

#ifdef FUNCTION_ID
  cout << "FV_solver_base::set_interface_tracer_flux ...returning.\n";
#endif  // FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################



void FV_solver_base::set_Hcorrection(
    cell *c,                ///< cell to operate on
    const axes axis,        ///< axis normal to interface.
    const pion_flt *edgeL,  ///< Left state
    const pion_flt *edgeR,  ///< right state
    const double g          ///< gamma
)
{
  if (axis != GetDirection()) {
    cout << GetDirection() << "\t";
    rep.error("bad direction in FV_solver_base::set_Hcorrection()", axis);
  }
  //
  // Sanders, Morano, Druguet, (1998, JCP, 145, 511)  eq. 10
  //
  double eta = 0.5
               * (fabs(edgeR[eqVX] - edgeL[eqVX])
                  + fabs(maxspeed(edgeR, g) - maxspeed(edgeL, g)));

  CI.set_Hcorr(c, axis, eta);
  return;
}

// ##################################################################
// ##################################################################

double FV_solver_base::select_Hcorr_eta(
    const cell *cl,  ///< cell to left of interface
    const cell *cr,  ///< cell to right of interface
    class GridBaseClass *grid)
{
  //
  // Based on Sanders, Morano, Druguet, (1998, JCP, 145, 511)
  // Fig. 9 and Eq. 16 (with obvious extension to 3D).
  //

  //
  // First get the current direction.
  //
  double eta     = 0.0;
  enum axes axis = GetDirection();
  //
  // First the interface between left and right cells.
  //
  eta = CI.get_Hcorr(cl, axis);

  //
  // If we are in 1D this is the only interface, so we can return.
  // Otherwise we have another 4 or 8 more interfaces to calculate.
  //
  if (FV_gndim == 1) return eta;

  //
  // Now the two in the positive direction of the 1st perpendicular
  // axis.  (i,j+1/2) and (i+1,j+1/2)
  //
  enum axes perp = static_cast<axes>((static_cast<int>(axis) + 1) % FV_gndim);
  eta            = max(eta, CI.get_Hcorr(cl, perp));
  eta            = max(eta, CI.get_Hcorr(cr, perp));
  //
  // Positive direction of the 2nd perp. axis, if 3D grid:
  // (i,j,k+1/2) (i+1,j,k+1/2)
  //
  if (FV_gndim > 2) {
    perp = static_cast<axes>((static_cast<int>(axis) + 2) % FV_gndim);
    eta  = max(eta, CI.get_Hcorr(cl, perp));
    eta  = max(eta, CI.get_Hcorr(cr, perp));
  }
  //
  // Negative direction of the (up to) two perp. axes.  We need to
  // make sure the cells exist, and if they don't we just ignore the
  // non-existent interface.  It's up to the grid to have all the
  // cells it needs, with the correct connectivity.
  //
  enum direction negdir = NO;
  cell *cneg            = 0;
  for (int idim = 1; idim < FV_gndim; idim++) {
    perp   = static_cast<axes>((static_cast<int>(axis) + idim) % FV_gndim);
    negdir = static_cast<direction>(static_cast<int>(axis) * 2);
    cneg   = grid->NextPt(cl, negdir);
    if (cneg) eta = max(eta, CI.get_Hcorr(cneg, perp));
    cneg = grid->NextPt(cr, negdir);
    if (cneg) eta = max(eta, CI.get_Hcorr(cneg, perp));
  }

  //
  // Will want to comment this out later...
  //
#ifndef NDEBUG
  cout << "cell id=" << cl->id << " axis=" << axis << ", eta_max=" << eta
       << "\n";
#endif  // NDEBUG

  return eta;
}

// ##################################################################
// ##################################################################
