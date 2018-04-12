///
/// \file solver_eqn_mhd_adi.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 written, based on old solver.cc classes.
///
/// Solver for the adiabatic ideal MHD Equations.
/// Calculates flux via either Lax-Friedrichs or Riemann solver.
/// Adds viscosity if asked for, and tracks flux of N passive tracers.
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.27 JM: Removed riemann_base references.  Added
///   function_id identifiers.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2011.04.15 JM: Change in UtoP() for tracers (to try to correct for negative density!).
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.04.14 JM: Moved flux solver to FV_solver

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"

#include "solver_eqn_mhd_adi.h"
using namespace std;

// *****************************************
// ***** FV SOLVER MHD Ideal Adiabatic *****
// *****************************************


// ##################################################################
// ##################################################################


FV_solver_mhd_ideal_adi::FV_solver_mhd_ideal_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno, ///< CFL number
      const double cellsize, ///< dx, cell size.
      const double gam, ///< gas eos gamma.
      pion_flt *state, ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr ///< Number of tracer variables.
      )
  : eqns_base(nv),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    HLLD_MHD(nv,gam),
    VectorOps_Cart(nd,cellsize)
{
#ifdef TESTING
  cout <<"::FV_solver_mhd_ideal_adi() constructor.\n";
  //cout <<"::FV_solver_mhd_ideal_adi() gamma = "<<eq_gamma<<"\n";
#endif

  negPGct = negROct = 0;
  return;
}


// ##################################################################
// ##################################################################

FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  cout <<"FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi() destructor.\n";
#endif

#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################



int FV_solver_mhd_ideal_adi::inviscid_flux(
      const cell *Cl, ///< Left state cell pointer
      const cell *Cr, ///< Right state cell pointer
      const pion_flt *Pl, ///< Left Primitive vector.
      const pion_flt *Pr, ///< Right Primitive vector.
      pion_flt *flux,///< Resultant Flux vector.
      pion_flt *pstar, ///< State vector at interface.
      const int solve_flag, ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS4=RoeRS)
      const double g        ///< Gas constant gamma.
      )
{
#ifdef TESTING
  //
  // Check input density and pressure are 'reasonably large'
  //
  if (Pl[eqRO]<TINYVALUE || Pl[eqPG]<TINYVALUE ||
      Pr[eqRO]<TINYVALUE || Pr[eqPG]<TINYVALUE) {
    rep.printVec("left ",Pl,eq_nvar);
    rep.printVec("right",Pr,eq_nvar);
    rep.error("FV_solver_mhd_ideal_adi::calculate_flux() Density/Pressure too small",Pl[eqRO]);
  }
#endif //TESTING
  int err=0;


  //
  // Set flux and pstar vector to zero.
  //
  for (int v=0;v<eq_nvar;v++) flux[v]  = 0.0;
  for (int v=0;v<eq_nvar;v++) pstar[v]  = 0.0;
  eq_gamma = g;
  

  //
  // Choose which Solver to use:
  //
  if      (solve_flag==FLUX_LF) {
    //
    // Lax-Friedrichs Method, so just get the flux
    //
    err += get_LaxFriedrichs_flux(Pl,Pr,flux,eq_gamma);
    for (int v=0;v<eq_nvar;v++) pstar[v] = 0.5*(Pl[v]+Pr[v]);
  }

  else if (solve_flag==FLUX_RSroe) {
    //
    // Roe Flux solver in conserved variables (Cargo and Gallice, 1997).
    // This is the symmetric version which sums over all waves.
    //
    err += MHD_Roe_CV_flux_solver_symmetric(Pl,Pr,eq_gamma,
					    HC_etamax,
					    pstar,flux);

    //
    // If we get an error, try the linear solver:
    //
    if (err) {
      //err=0;
      cout <<"ROE SOLVER FAILED -- TRYING FKJ98 SOLVER. err="<<err<<"\n";
      err = JMs_riemann_solve(Pl,Pr,pstar,1,eq_gamma);
      PtoFlux(pstar, flux, eq_gamma);
    }
  }

  else if (solve_flag==FLUX_RSlinear ||
	   solve_flag==FLUX_RSexact  ||
	   solve_flag==FLUX_RShybrid) {
    //
    // JM's hybrid Riemann Solver -- Falle et al. (1998) with the Roe
    // and Balsara (1996) eigenvector normalisation.
    //
    err += JMs_riemann_solve(Pl,Pr,pstar,solve_flag,eq_gamma);
    PtoFlux(pstar, flux, eq_gamma);
  }

  else if (solve_flag==FLUX_RS_HLLD) {
    err += MHD_HLLD_flux_solver(Pl, Pr, eq_gamma, flux);
      //const pion_flt *left,  ///< input left state
      //const pion_flt *right, ///< input right state
      //const double gamma,    ///< input gamma
      ////const pion_flt etamax, ///< H-correction eta-max value.
      ////pion_flt *out_ps,       ///< output p*
      ////pion_flt *out_pss,       ///< output p**
      //pion_flt *out_flux         ///< output flux
      //)
  }


  else {
    rep.error("what sort of flux solver do you mean???",solve_flag);
  }

  //
  // That should do it for the flux
  //
  return err;
}


// ##################################################################
// ##################################################################



int FV_solver_mhd_ideal_adi::AVFalle(
      const pion_flt *Pleft,
      const pion_flt *Pright,
      const pion_flt *Pstar,
      pion_flt *flux, 
      const double eta, ///< already set as FS_etav
      const double gamma ///< already set as eq_gamma
      )
{
  /// \section Equations
  /// The equations are as follows (for flux along x-direction):
  /// \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*}_f (v_{xR}-v_{xL})  \,,\f]
  /// \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*}_f (v_{yR}-v_{yL})  \,,\f]
  /// \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*}_f (v_{zR}-v_{zL})  \,,\f]
  /// \f[ F[e]   = F[e]   - \eta \rho^{*} c^{*}_f \left[ v^{*}_x (v_{xR}-v_{xL}) + v^{*}_y (v_{yR}-v_{yL}) + v^{*}_z (v_{zR}-v_{zL}) \right]\,.\f]
  /// This follows from assuming a 1D problem, with non-zero bulk viscosity, and
  /// a non-zero shear viscosity.  The bulk viscosity gives the viscosity 
  /// in the x-direction, opposing stretching and compression.  The shear
  /// viscosity gives the y and z terms, opposing slipping across boundaries.
  /// 
  /// The identification of the velocities and density with the
  /// resultant state from the Riemann Solver is a bit ad-hoc, and could easily
  /// be changed to using the mean of the left and right states, or something
  /// else. (See my notes on artificial viscosity).
  /// 
  /// Andy says Sam also has similar terms for the Magnetic field, although 
  /// it isn't in his paper.  Andy uses this too.
  /// \f[ f[B_y] = f[B_y] - \eta c^{*}_f (B_{yR}-B_{yL}) \;, \f]
  /// \f[ f[B_z] = f[B_z] - \eta c^{*}_f (B_{zR}-B_{zL}) \;, \f]
  /// with associated terms added to energy flux:
  /// \f[ F[e]   = F[e]   - \eta c^{*}_f \left[ B^{*}_z (B_{zR}-B_{zL}) +B^{*}_y (B_{yR}-B_{yL})\right] \;. \f]
  /// 
  /// Andy uses the mean of the left and right state velocities/field to calculate
  /// the energy flux; I am using the starred state velocity/field.  Not sure
  /// if that is bad, or if it will make any discernible difference at all.
  /// After some testing, it makes very little difference.  It may help with
  /// some crazy cases, but I haven't encountered them yet.
  /// 
  /// I have switched from using \f$c^{*}_f\f$ to using the fast speed from 
  /// the mean vector of the left and right states.  This is because the 
  /// Riemann Solver can return negative density, which then gets set to a 
  /// density floor value, but the B-field is still normal, so the fast speed
  /// blows up.  If I use the mean state though, it represents a typical speed
  /// between the two states, so it is in principle no less realistic.
  ///
  double prefactor = cfast_components(0.5*(Pleft[eqRO]+Pright[eqRO]),
				      0.5*(Pleft[eqPG]+Pright[eqPG]),
				      0.5*(Pleft[eqBX]+Pright[eqBX]),
				      0.5*(Pleft[eqBY]+Pright[eqBY]),
				      0.5*(Pleft[eqBZ]+Pright[eqBZ]),
				      eq_gamma)*FS_etav*Pstar[eqRO];
				      
  //
  // Momentum flux
  //
  double momvisc = prefactor*(Pright[eqVX]-Pleft[eqVX]);//*4./3.;
  double ergvisc = momvisc*Pstar[eqVX];
  flux[eqMMX] -= momvisc;
  momvisc = prefactor*(Pright[eqVY]-Pleft[eqVY]);
  flux[eqMMY] -= momvisc;
  ergvisc += momvisc*Pstar[eqVY];
  momvisc = prefactor*(Pright[eqVZ]-Pleft[eqVZ]);
  flux[eqMMZ] -= momvisc;
  ergvisc += momvisc*Pstar[eqVZ];

  //
  // Magnetic field flux.
  //
  prefactor *= FS_etaB/(FS_etav*Pstar[eqRO]);
  momvisc = prefactor*(Pright[eqBY]-Pleft[eqBY]);
  flux[eqBBY] -= momvisc;
  ergvisc += momvisc*Pstar[eqBY];
  momvisc = prefactor*(Pright[eqBZ]-Pleft[eqBZ]);
  flux[eqBBZ] -= momvisc;
  ergvisc += momvisc*Pstar[eqBZ];
  flux[eqERG] -= ergvisc;
  return(0);  
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::PtoU(
      const pion_flt *p,
      pion_flt *u,
      const double g
      )
{
  eqns_mhd_ideal::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];
  return;
}


// ##################################################################
// ##################################################################


int FV_solver_mhd_ideal_adi::UtoP(
      const pion_flt *u,
      pion_flt *p,
      const double MinTemp, ///< Min Temperature allowed on grid.
      const double g
      )
{
  int err=eqns_mhd_ideal::UtoP(u,p,MinTemp,g);
  // we use u[eqRO] because if there was a negative density, then usually
  // the tracer sign will follow, and this way we get a positive primitive
  // variable tracer back.
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/u[eqRO];
  return err;
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::PUtoFlux(
      const pion_flt *p,
      const pion_flt *u,
      pion_flt *f
      )
{
  eqns_mhd_ideal::PUtoFlux(p,u,f);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = p[eqTR[t]]*f[eqRHO];
  return;
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::UtoFlux(
      const pion_flt *u,
      pion_flt *f,
      const double g
      )
{
  eqns_mhd_ideal::UtoFlux(u,f,g);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = u[eqTR[t]]*f[eqRHO]/u[eqRHO];
  return;
}


// ##################################################################
// ##################################################################



///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_mhd_ideal_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c,          // Current cell.
        const axes d,     // Which axis we are looking along.
        const pion_flt *fn, // Negative direction flux.
        const pion_flt *fp, // Positive direction flux.
        const pion_flt *,   // slope vector for cell c.
        const int,        // spatial order of accuracy.
        const double,     // cell length dx.
        const double      // cell TimeStep, dt.
        )
{
  pion_flt u1[eq_nvar];
  //
  // This calculates -dF/dx
  //
  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  return(err);
}


// ##################################################################
// ##################################################################

///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_mhd_ideal_adi::CellAdvanceTime(
      class cell *c, // cell to update.
      const pion_flt *Pin, // Initial State Vector.
      pion_flt *dU, // Update vector dU
      pion_flt *Pf, // Final state vector (can be same as initial vec.).
      pion_flt *dE, // Tracks change of energy if I have to correct for negative pressure
      const double, // gas EOS gamma.
      const double MinTemp, ///< Min Temperature allowed on grid.
      const double  // Cell timestep dt.
      )
{
  pion_flt u1[eq_nvar], u2[eq_nvar];
  //
  // First convert from Primitive to Conserved Variables
  //
  PtoU(Pin,u1, eq_gamma);

  //
  // Now add dU[] to U[], and change back to primitive variables.
  // This can give negative pressures, so check for that and fix it if needed.
  //
  for (int v=0;v<eq_nvar;v++) {
    u1[v] += dU[v];   // Update conserved variables
    dU[v] = 0.;       // Reset the dU array for the next timestep.
  }
  if(UtoP(u1,Pf,MinTemp, eq_gamma)!=0) {
    cout<<"(FV_solver_mhd_ideal_adi::CellAdvanceTime) UtoP complained (maybe about negative pressure...) fixing\n";
    PtoU(Pf, u2, eq_gamma);
    *dE += (u2[ERG]-u1[ERG]);
    UtoP(u2,Pf,MinTemp, eq_gamma);
  }

  return 0;
}



// ##################################################################
// ##################################################################


///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_mhd_ideal_adi::CellTimeStep(
        const cell *c, ///< pointer to cell
        const double, ///< gas EOS gamma.
        const double  ///< Cell size dx.
        )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellTimeStep ...starting.\n";
#endif //FUNCTION_ID

  //
  // Get Max velocity along a grid direction.
  //
  pion_flt u1[eq_nvar];
  pion_flt temp = fabs(c->P[eqVX]);
  if (FV_gndim>1) temp = max(temp,static_cast<pion_flt>(fabs(c->P[eqVY])));
  if (FV_gndim>2) temp = max(temp,static_cast<pion_flt>(fabs(c->P[eqVZ])));
  
  /*  if (fabs(c->P[VX])>fabs(c->P[VY])) {temp = fabs(c->P[VX]);}
   else {temp = fabs(c->P[VY]);}
   if (fabs(c->P[VZ])>temp2) temp2=c->P[VZ];
   */
  
  //
  // First rotate the state vector to the fastest directions,
  // then add the fast speed to the max-velocity, and this is
  // the max wavespeed.
  //
  enum axes newdir;
  if (FV_gndim==1) temp += cfast(c->P,eq_gamma);
  else { // We may have to rotate the state vector to find the fastest fast speed.
    newdir = XX;
    if (fabs(c->P[BY])<fabs(c->P[BX])) {
      newdir =YY;
      if (fabs(c->P[BZ])<fabs(c->P[BY])) newdir=ZZ;
    }
    else if (fabs(c->P[BZ])<fabs(c->P[BX])) newdir=ZZ;

    if (newdir !=XX) {
      for (int v=0;v<eq_nvar;v++) u1[v] = c->P[v];
      rotate(u1,XX,newdir);
      temp += cfast(u1, eq_gamma);
    }
    else temp += cfast(c->P, eq_gamma);
  }
  FV_dt = FV_dx/temp;

  //
  // Check the gradient of pressure with neighbouring cells, since this can 
  // dramatically shorten the timestep.  (CAN'T REMEMBER WHY!!!)
  //
  //if (c->isgd) {
  //  int pg = static_cast<int>(eqPG);
  //  double grad = max_grad_abs(c,0,pg, grid)/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) FV_dt /= temp;
  //}
  //else if (grid->NextPt(c,XP)) {
  //  double grad = fabs(c->P[PG]-grid->NextPt(c,XP)->P[PG])/FV_dx/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) {
  //    FV_dt /= temp;
  //  }
  //}
  //else rep.error("No neighbour to gradient test",c);    

  
  //
  // Now scale the max. allowed timestep by the CFL number we are using (<1).
  //
  FV_dt *= FV_cfl;

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellTimeStep ...returning.\n";
#endif //FUNCTION_ID
  return FV_dt;
}



// ##################################################################
// ##################################################################



// ********************************************************************************
// FV_solver_mhd_mixedGLM_adi class, for the Dedner-GLM divergence cleaning method.
// ********************************************************************************


// ##################################################################
// ##################################################################


FV_solver_mhd_mixedGLM_adi::FV_solver_mhd_mixedGLM_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno,   ///< CFL number
      const double cellsize,    ///< dx, cell size.
      const double gam,     ///< gas eos gamma.
      pion_flt *state,     ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr         ///< Number of tracer variables.
      )
  : eqns_base(nv),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    HLLD_MHD(nv,gam),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_mhd_mixedGLM(nv),
{
  return;
}


// ##################################################################
// ##################################################################

FV_solver_mhd_mixedGLM_adi::~FV_solver_mhd_mixedGLM_adi()
{
  return;
}


// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::inviscid_flux(
      const cell *Cl, ///< Left state cell pointer
      const cell *Cr, ///< Right state cell pointer
      const pion_flt *Pl, ///< Left Primitive state vector.
      const pion_flt *Pr, ///< Right Primitive state vector.
      pion_flt *flux, ///< Resultant Flux state vector.
      pion_flt *pstar, ///< State vector at interface.
      const int solve_flag, ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS4=RoeRS)
      const double g ///< Gas constant gamma.
      )
{
#ifdef TESTING
  //
  // Check input density and pressure are positive
  //
  if (Pl[eqRO]<TINYVALUE || Pl[eqPG]<TINYVALUE ||
      Pr[eqRO]<TINYVALUE || Pr[eqPG]<TINYVALUE) {
    rep.printVec("left ",Pl,eq_nvar);
    rep.printVec("right",Pr,eq_nvar);
    rep.error("FV_solver_mhd_mixedGLM_adi::calculate_flux() Density/Pressure too small",Pl[eqRO]);
  }
#endif //TESTING
  
  int err=0;

  //
  // Set flux and pstar vector to zero.
  //
  for (int v=0;v<eq_nvar;v++) flux[v]  = 0.0;
  for (int v=0;v<eq_nvar;v++) pstar[v]  = 0.0;
  eq_gamma = g;
  //
  // Need temporary left and right state vectors b/c we need to change
  // the left and right state values of eqBX to the resolved state of
  // the Dedner et al. (2002) Riemann problem for (Bx,Psi).
  //
  pion_flt left[eq_nvar], right[eq_nvar];
  for (int v=0;v<eq_nvar;v++) left[v]  = Pl[v];
  for (int v=0;v<eq_nvar;v++) right[v] = Pr[v];
  
  ///
  /// First get the resolved state in the 2x2 Dedner Riemann Problem,
  /// which is described well in their paper, and/or my thesis.
  ///
  /// Uses Dedner eq.41 for the flux in Bx and Psi:
  /// \f[ \partial_t B_x + \partial_x \psi = 0 \;, \qquad 
  ///     \partial_t \psi + \partial_x (c_h^2 B_x) = 0 \;. \f]
  /// where the source term has been omitted, as it is calculated separately.
  /// 
  /// The GLM method has Bx and Psi decoupled from all the other variables
  /// in the Riemann Problem, so they can be solved separately as a two
  /// variable system (Dedner eq.42)
  /// 
  /// \f[ F(\psi) = c_h^2 B_x^* = c_h^2 \left( \frac{1}{2}(B_x(L)+B_x(R))
  ///                             - \frac{1}{2c_h}(\psi_R-\psi_L) \right) \f]
  /// \f[ F(B_x) = \psi_* = \frac{1}{2}(\psi_L+\psi_R) - \frac{c_h}{2}(B_x(R)-B_X(L)) \f]
  /// 
  /// Bx(*) is then used for both the left and right state to calculate
  /// the flux.
  /// 
  double psistar = 0.5*(left[eqSI]+right[eqSI]
			-GLM_chyp*(right[eqBX]-left[eqBX]));
  double bxstar  = 0.5*(left[eqBX]+right[eqBX]
			-(right[eqSI]-left[eqSI])/GLM_chyp);
  left[eqBX] = right[eqBX] = bxstar;

  //
  // Now continue on in an identical manner as the ideal MHD solver, 
  // by calling its flux solver:
  //
  err=FV_solver_mhd_ideal_adi::inviscid_flux(Cl,Cr,left,right,
					       flux,pstar,
					       solve_flag,eq_gamma);


  //
  // Now we have to add in the flux in BX and PSI, based on Dedner et
  // al.'s method.
  //
  // NOTE: Dedner doesn't say what to do about the energy flux, so this is a
  // modification to ensure consistency (Mackey & Lim, 2011).
  // Derigs et al. (2017/2018) have a more consistent calculation...
  // 
  flux[eqPSI]  = GLM_chyp*GLM_chyp*bxstar;
  flux[eqBBX]  = psistar;
  flux[eqERG] += bxstar*psistar;
  return err;
}



// ##################################################################
// ##################################################################



///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_mhd_mixedGLM_adi::CellAdvanceTime(
      class cell *c,
      const pion_flt *Pin, // Initial State Vector.
      pion_flt *dU, // Update vector dU
      pion_flt *Pf, // Final state vector (can be same as initial vec.).
      pion_flt *dE, // Tracks change of energy if I have to correct for negative pressure
      const double, // gas EOS gamma.
      const double MinTemp, ///< Min Temperature allowed on grid.
      const double // Cell timestep dt.
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...starting.\n";
#endif //FUNCTION_ID

  int err=FV_solver_mhd_ideal_adi::CellAdvanceTime(c,Pin,dU,Pf,dE,0,MinTemp,0);
  GLMsource(&(Pf[eqSI]),FV_dt);

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...returning.\n";
#endif //FUNCTION_ID
  return(err);
}


// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::PtoU(
      const pion_flt *p,
      pion_flt *u,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::PtoU ...starting.\n";
#endif //FUNCTION_ID

  eqns_mhd_mixedGLM::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::UtoP(
      const pion_flt *u,
      pion_flt *p,
      const double MinTemp, ///< Min Temperature allowed on grid.
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::UtoP ...starting.\n";
#endif //FUNCTION_ID

  int err=eqns_mhd_mixedGLM::UtoP(u,p,MinTemp,g);
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/p[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::UtoP ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::GotTimestep(
        const double delt ///< timestep dt.
	)
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep ...starting.\n";
#endif //FUNCTION_ID

  //     cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep() setting wave speeds.\n";
  GLMsetPsiSpeed(FV_cfl,FV_dx,delt);

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep ...returning.\n";
#endif //FUNCTION_ID
}



// ##################################################################
// ##################################################################




/// ---------------------------------------------------------------------
/// -------------------  AXI-SYMMETRIC EQUATIONS ------------------------
/// ---------------------------------------------------------------------




// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_ideal_adi::cyl_FV_solver_mhd_ideal_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno, ///< CFL number
      const double cellsize, ///< dx, cell size.
      const double gam, ///< gas eos gamma.
      pion_flt *state,   ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr    ///< Number of tracer variables.
      )
  : eqns_base(nv),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    HLLD_MHD(nv,gam),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

//  cout <<"cyl_FV_solver_mhd_ideal_adi CONSTRUCTOR\n";
//  cout <<"IdealMHD Equations; Riemann Solver Method; Cylindrical Coordinates.\n";
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for 2d axial symmetry \
			 so far.  Sort it out!",nd);

#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


cyl_FV_solver_mhd_ideal_adi::~cyl_FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################

int cyl_FV_solver_mhd_ideal_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c, ///< Current cell.
        const axes d, ///< Which axis we are looking along.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        const pion_flt *dpdx, ///< slope vector for cell c.
        const int OA,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        )
{
#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_ideal_adi::dU_Cell ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar];

  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  if (d==Rcyl) {
    double pm = (c->Ph[eqBX]*c->Ph[eqBX] +c->Ph[eqBY]*c->Ph[eqBY] +c->Ph[eqBZ]*c->Ph[eqBZ])/2.;
    switch (OA) {
     case OA1:
       c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm)/CI.get_dpos(c,Rcyl);
      break;
     case OA2:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm + 
			  (CI.get_dpos(c,Rcyl)-R_com(c))*
		   (dpdx[eqPG] +c->Ph[eqBX]*dpdx[eqBX] +c->Ph[eqBY]*dpdx[eqBY] +c->Ph[eqBZ]*dpdx[eqBZ]))
                   /CI.get_dpos(c,Rcyl);
      break;
     default:
      rep.error("Bad OOA in cyl_IdealMHD_RS::dU, only know 1st,2nd",OA);
    }
  }

#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_ideal_adi::dU_Cell ...returning.\n";
#endif //FUNCTION_ID
  return(err);
} 


// ##################################################################
// ##################################################################

//------------------------------------------------------------------//
cyl_FV_solver_mhd_mixedGLM_adi::cyl_FV_solver_mhd_mixedGLM_adi(
        const int nv, ///< number of variables in state vector.
        const int nd, ///< number of space dimensions in grid.
        const double cflno, ///< CFL number
        const double cellsize, ///< dx, cell size.
        const double gam, ///< gas eos gamma.
        pion_flt *state,   ///< State vector of mean values for simulation.
        const double avcoeff, ///< Artificial Viscosity Parameter etav.
        const int ntr    ///< Number of tracer variables.
        )
  : eqns_base(nv),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam), 
    Riemann_Roe_MHD_CV(nv,gam),
    HLLD_MHD(nv,gam),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_mhd_mixedGLM(nv),
    FV_solver_mhd_mixedGLM_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

//  cout <<"cyl_FV_solver_mhd_mixedGLM_adi CONSTRUCTOR\n";
//  cout <<"glmMHD Equations; Riemann Solver Method; Cylindrical Coordinates.\n";
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for \
                        2d axial symmetry so far.  Sort it out!",nd);
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_mixedGLM_adi::~cyl_FV_solver_mhd_mixedGLM_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################

int cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c, ///< Current cell.
        const axes d, ///< Which axis we are looking along.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        const pion_flt *dpdx, ///< slope vector for cell c.
        const int OA,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        )
{
#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar];
  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  if (d==Rcyl) {
    double pm = (c->Ph[eqBX]*c->Ph[eqBX] +c->Ph[eqBY]*c->Ph[eqBY] +c->Ph[eqBZ]*c->Ph[eqBZ])/2.;
    switch (OA) {
     case OA1:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm)/CI.get_dpos(c,Rcyl);
      c->dU[eqBBX] += FV_dt*c->Ph[eqSI]/CI.get_dpos(c,Rcyl);
      break;
     case OA2:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm + 
			 (CI.get_dpos(c,Rcyl)-R_com(c))*
			 (dpdx[eqPG] +c->Ph[eqBX]*dpdx[eqBX] +c->Ph[eqBY]*dpdx[eqBY] +c->Ph[eqBZ]*dpdx[eqBZ]))
	                /CI.get_dpos(c,Rcyl);
      c->dU[eqBBX] += FV_dt*(c->Ph[eqSI] +(CI.get_dpos(c,Rcyl)-R_com(c))*dpdx[eqSI])/CI.get_dpos(c,Rcyl);
      break;
     default:
      rep.error("Bad OOA in cyl_glmMHD_RS::dU, only know 1st,2nd",OA);
    }
  }

#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell ...returning.\n";
#endif //FUNCTION_ID
  return(err);
} 


// ##################################################################
// ##################################################################

