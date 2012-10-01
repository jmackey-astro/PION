///
/// \file flux_hydro_isothermal.cc
/// \author Jonathan Mackey
/// 
/// This file has the class definitions for the isothermal hydrodynamics flux solver.
/// Solvers implemented: Linear Riemann Solver (plus Lax-Friedrichs in derived class)
///
/// - 2010-09-21 JM: THIS CODE SHOULD BE REGARDED AS UNTESTED, AND SO SHOULD BE 
///    TESTED AND TESTED AND TESTED SOME MORE BEFORE BEING USED.  IT HAS NEVER BEEN
///    USED FOR ANY SCIENCE.
///
/// - 2010-02-10 JM: put in a warning not to use the linear solver b/c
///    it's rubbish in its present form.
///
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///    Renamed calculate_flux() to inviscid_flux() and moved AV
///    calculation to FV_solver_base class.
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
///
#ifdef ISOTHERMAL_SOLVERS_ENABLED

#include "flux_hydro_isothermal.h"
using namespace std;


// *************************************************
// ***** ISOTHERMAL HYDRO RIEMANN SOLVER CLASS *****
// *************************************************
riemann_hydro_iso::riemann_hydro_iso(const int nv,       ///< number of variables
				     const double *state ///< mean state for reference
				     )
  :
  eqns_base(nv), riemann_base(nv), eqns_IsoEuler(nv)
{
  cout <<"riemann_hydro_iso::riemann_hydro_iso constructor.\n";
  SetAvgState(state,0);
}

riemann_hydro_iso::~riemann_hydro_iso()
{
  cout <<"riemann_hydro_iso::riemann_hydro_iso destructor.\n";
}

void riemann_hydro_iso::SetAvgState(const double *input_state, ///< input state vector == avg state.
				   const double               ///< unused (gamma)
				   )
{
  RS_refvec[eqRO] = input_state[eqRO];
  RS_refvec[eqAA] = input_state[eqAA];
  //
  // Set reference velocities to be the isothermal sound speed.
  //
  RS_refvec[eqVX] = input_state[eqAA];
  RS_refvec[eqVY] = input_state[eqAA];
  RS_refvec[eqVZ] = input_state[eqAA];
  
  return;

}


int riemann_hydro_iso::riemann_solve(const double *,  ///< Class already knows about left,right,pstar, and where to find them
				    const double *,
				    double *,
				    const int solver_flag, ///< type of solve.
				    const double           ///< isothermal equations have no gamma.
				    )
{
  int err=0;
  
  //
  // First test if the left and right states are the same, and if they are
  // then return the simple answer.
  //
  double diff=0.;
  for(int i=0;i<eq_nvar; i++)
    diff += fabs(RS_right[i]-RS_left[i])/(fabs(RS_refvec[i])+TINYVALUE);
  if (diff <1.e-6) {
    for(int i=0;i<eq_nvar; i++) RS_pstar[i] = (RS_left[i]+RS_right[i])/2.;
    return 0;
  }

  //
  // If we get this far we have a calculation to do, so call the appropriate solver:
  //
  if      (solver_flag==FLUX_RSlinear) { // Linear Solver:
    err += linear_solver_hydro_iso();
  }
  else if (solver_flag==FLUX_RSexact) {
    err += two_shock_solver_hydro_iso();
  }
  else rep.error("Only Linear solver written for Isothermal Hydro so far",solver_flag);
    
  //
  // Check for negative pressures and densities!
  //
  if (RS_pstar[eqAA]<=TINYVALUE) {
    rep.error("Got negative isothermal sound speed in riemann_hydro_iso!!!", RS_pstar[eqAA]);
    RS_pstar[eqAA] = BASEPG*RS_refvec[eqAA];
  }
  if (RS_pstar[eqRO]<=TINYVALUE) {
    rep.error("Got negative DENSITY in riemann_hydro_iso!!!", RS_pstar[eqRO]);
    RS_pstar[eqRO] = BASEPG*RS_refvec[eqRO];
  }

 return err;
}
 
  
int riemann_hydro_iso::linear_solver_hydro_iso()
{
  //
  // If Using straight averaging, can do it this way.  Roe Average might be better?
  //
  for (int v=0;v<eq_nvar;v++)
    RS_meanp[v] = 0.5*(RS_left[v] +RS_right[v]);
  cout <<"WARNING! LINEAR ISOTHERMAL SOLVER IS RUBBISH! USE A BETTER AVERAGING!\n";

  //
  // Set the eigenvalues:
  //
  double 
    lneg=RS_meanp[eqVX]-RS_meanp[eqAA],
    lpos=RS_meanp[eqVX]+RS_meanp[eqAA];
  //
  // If the negative e-value is positive, then the solution is the left state.
  // Conversely if the positive e-value is negative, the solution is the right state.
  //
  if      (lneg >= 0.0) {  // if(u_av-c_av>0) pstar = left state
    for (int i=0; i<eq_nvar; i++) RS_pstar[i] = RS_left[i];
    return 0;
  }
  else if (lpos <=0.) {  //else if(u_av+c_av<=0) pstar = right state
    for (int i=0; i<eq_nvar; i++) RS_pstar[i] = RS_right[i];
    return 0;
  }
  else {  // else we are in the starred region, so p_star, u_star are right
    RS_pstar[eqRO] = (lpos*RS_right[eqRO] -lneg*RS_left[eqRO]  -         (RS_right[eqVX]-RS_left[eqVX])) /(2.0*RS_meanp[eqAA]);
    RS_pstar[eqVX] = (lpos*lneg*(RS_right[eqRO]-RS_left[eqRO]) -lneg*RS_right[eqVX] +lpos*RS_left[eqVX]) /(2.0*RS_meanp[eqAA]);
  }
  
  //
  // Need to assign values to v_y,v_z if present.
  // They only change across the contact discontinuity.
  //
  if (RS_pstar[eqVX]>0) {
    RS_pstar[eqVY] = RS_left[eqVY];		    
    RS_pstar[eqVZ] = RS_left[eqVZ];
  }
  else {
    RS_pstar[eqVY] = RS_right[eqVY];		    
    RS_pstar[eqVZ] = RS_right[eqVZ];
  }
  //  cout << "(Isothermal Linear Riemann Solver) Success!" << "\n";
  return 0;
}




// This version is specific to solving the equation in the Exact Riemann Solver.
// Other versions could be specified, taking more or fewer parameters.
int riemann_hydro_iso::FR_find_root(double *ans, /**< pointer to result */
				    const double p1,  ///< rho_left
				    const double p2,  ///< rho_right
				    const double,  ///< parameter 3
				    const double,  ///< parameter 4
				    const double   ///< parameter 5
				    )
{
  int err=0;
 
  FR_param1 = p1;  //  rho_l
  FR_param2 = p2;  //  rho_r
  
  //
  // Set initial values for x1, x2
  // Modify this for your equations
  //  double x1 = 1.;
  //  double x2 = 2.;
  // My initial guesses are 1/3 of and 3 times the arithmetic mean 
  // of the left and right densities.
  // 
  double x1 = (FR_param1+ FR_param2)/6.0;
  double x2 = x1*9.0;

  //
  // This guess is just [0,1], which is not a bad starting point.
  //  double x1 = 0.;
  //  double x2 = 1.;
  // 
  
  //
  // Call the common positive-definite root-finder,
  // now that parameters are set properly.
  // 
  err = findroot::solve_pos(x1, x2, ans);
  
  return err;
}

double riemann_hydro_iso::FR_root_function(double pp)
{
  //
  // Wrapper function to solve across both waves given p*=pp, which returns
  // the velocity difference (u*(R) - u*(L)) = *uu assuming both waves are
  // shocks.
  // 
  // This takes in the test value of p*, and calculates the resulting
  // values of u* for the left and right waves, and returns the 
  // difference.
  // 
  
  //
  // First get the square root of the density ratios pstar/pleft, pstar/pright
  // 
  double 
    KL = sqrt(pp/RS_left[eqRO]),
    KR = sqrt(pp/RS_right[eqRO]);
  
  //
  // Next return the velocity difference:
  // v*R-v*L = VR-VL -a*[(KR*KR-1)/KR -(KL*KL-1)/KL]
  // I derived this from the jump conditions, not sure if there is a reference for 
  // it -- look in my thesis if interested!!!
  // 
  double ustarL = (KL-1.0/KL);
  if (KL>=1.0)
    ustarL = RS_left[eqVX] - RS_left[eqAA]*ustarL;
  else
    ustarL = RS_left[eqVX] + RS_left[eqAA]*ustarL;
  
  double ustarR = (KR-1.0/KR);
  if (KR>=1.0)
    ustarR = RS_right[eqVX] +RS_right[eqAA]*ustarR;
  else 
    ustarR = RS_right[eqVX] -RS_right[eqAA]*ustarR;
 
  return ustarR-ustarL;  
//  return ( (RS_right[eqVX]-RS_left[eqVX]) -
//	   RS_right[eqAA]*((KR-1.0/KR)-(KL-1.0/KL)) );
}

int riemann_hydro_iso::two_shock_solver_hydro_iso()
{
  int err=0;
  
  //
  // Rootfinder. This should have found the value of rho* which matches u*(L) and u*(R)
  // across the contact discontinuity.
  //
  err += FR_find_root(&(RS_pstar[eqRO]), RS_left[eqRO], RS_right[eqRO], 0, 0, 0);

  //
  // Given a converged value for rho*, we can get u*
  // For testing I'm going to get u*L and u*R and compare them.
  // 
  double K = sqrt(RS_pstar[eqRO]/RS_left[eqRO]);
  double ustarL = RS_left[eqVX]  -  RS_left[eqAA] *(K-1.0/K);
  K = sqrt(RS_pstar[eqRO]/RS_right[eqRO]);
  double ustarR = RS_right[eqVX] -  RS_right[eqAA]*(K-1.0/K);
  if (!GS.equalD(ustarL,ustarR)) {
    cout <<"Left and right velocities don't match!!! u*L="<<ustarL<<" and u*R="<<ustarR<<" rho*="<<RS_pstar[eqRO]<<"\n";
    err += 1;
  }
  
  //
  // So now we have rho*,u*, and we know a=const, so we are done.
  // 
  RS_pstar[eqVX] = ustarL;
  
  //
  // Need to assign values to v_y,v_z if present.
  // They only change across the contact discontinuity.
  //
  if (RS_pstar[eqVX]>0) {
    RS_pstar[eqVY] = RS_left[eqVY];		    
    RS_pstar[eqVZ] = RS_left[eqVZ];
  }
  else {
    RS_pstar[eqVY] = RS_right[eqVY];		    
    RS_pstar[eqVZ] = RS_right[eqVZ];
  }
 
  return err;
}

// **********************************************
// ***** ISOTHERMAL HYDRO FLUX SOLVER CLASS *****
// **********************************************



flux_solver_hydro_iso::flux_solver_hydro_iso(const int nv,        ///< Number of variables in state vector
					     const double *state, ///< state vector which is 'typical' in the problem being solved. 
					     const double eta,     ///< coefficient of artificial viscosity (may or may not be used!)
					     const int ntr         ///< Number of tracer variables.
					     )
  :
  eqns_base(nv), riemann_base(nv), flux_solver_base(nv, eta, ntr),
  eqns_IsoEuler(nv), riemann_hydro_iso(nv,state)
{
  cout <<"Setting up the Isothemal Hydro Flux Solver.\n";
  cout <<"Default solver set to "<<SimPM.solverType<<" where ";
  cout <<"0=LF,1=RSlin,2=RSexact,3=RShybrid,4=RSRoe,5=Roe-PV,6=FVS.\n";
}

flux_solver_hydro_iso::~flux_solver_hydro_iso()
{
}

int flux_solver_hydro_iso::inviscid_flux(const cell *Cl, ///< Left state cell pointer
					 const cell *Cr, ///< Right state cell pointer
					 const double *l,      ///< Left Primitive state vector.
					 const double *r,      ///< Right Primitive state vector.
					 double *out_flux,         ///< Resultant Flux state vector.
					 const int solve_flag, ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
					 const double          ///< Gas constant gamma.
					 )
{
  //
  // Check input density and pressure are 'reasonably large'
  //
  if (l[eqRO]<TINYVALUE || l[eqAA]<TINYVALUE ||
      r[eqRO]<TINYVALUE || r[eqAA]<TINYVALUE) {
    rep.printVec("left ",l,eq_nvar);
    rep.printVec("right",r,eq_nvar);
    rep.error("flux_solver_hydro_iso::calculate_flux() Density/Pressure too small",l[eqRO]);
  }
  int err=0;

  //
  // Copy over left and right states, and calculate the average state.
  //  
  for (int v=0;v<eq_nvar;v++) {
    RS_left[v] = l[v];
    RS_right[v] = r[v];
  }

  //
  // Choose which Solver to use:
  //
  if      (solve_flag==FLUX_LF) {
    //
    // Lax-Friedrichs Method, so just get the flux
    //
    err += get_LaxFriedrichs_flux(RS_left,RS_right,FS_flux,0);
  }

  else if (solve_flag==FLUX_RSlinear ||
	   solve_flag==FLUX_RSexact  ||
	   solve_flag==FLUX_RShybrid ||
	   solve_flag==FLUX_RSroe) {
    
    //
    // First see if we are getting the Roe Flux:
    // 
    if (solve_flag==FLUX_RSroe) {
      err += Roe_solver_hydro_iso();
    }
    
    else {
      //
      // All Riemann Solver Methods, so call the solver:
      //
      err += riemann_solve(RS_left,RS_right,RS_pstar,solve_flag,0);
      
      //
      // Convert pstar to a flux:
      //
      PtoFlux(RS_pstar, FS_flux, 0);
    }
    
    //
    // That should do it for the Riemann Solver flux
  }

  else {
    rep.error("what sort of flux solver do you mean???",solve_flag);
  }

  for (int v=0;v<eq_nvar;v++) out_flux[v] = FS_flux[v];

  return err;
}

#ifdef LAPIDUS_VISCOSITY_ENABLED
int flux_solver_hydro_iso::AVLapidus(const cell *Cl, ///< Left state cell pointer
				     const cell *Cr, ///< Right state cell pointer
				     double *flux, 
				     const double etav,
				     const double
				     )
{
  ///
  /// This is not working!!! Don't use it!
  /// THIS USES A VERY BAD DIV(V) APPROXIMATION -- REALLY NEED TO CALCULATE IT
  /// PROPERLY, SO MAYBE THE SOLVER CLASS CAN DO THIS???
  /// 2010.09.30: SEE ADIABATIC SOLVER!
  ///
  rep.error("Fix Lapidus for Isothermal!!!",1);

  double ul[5], ur[5]; // Euler equations only affect first 5 vars.
  PtoU(RS_left,ul,0); PtoU(RS_right,ur,0);
  double divu = ((RS_right[eqVX]-RS_left[eqVX])+(RS_right[eqVY]-RS_left[eqVY])+(RS_right[eqVZ]-RS_left[eqVZ]));
  double vt; int i=0;
  for (int v=0;v<5;v++) {
    vt = etav*divu*(ur[v]-ul[v]);
    if (fabs(vt)>= fabs(flux[v]) && fabs(flux[v])>BASEPG) {
      cout <<"\t\t\t\t\tlarge visc: var="<<v<<" vt= "<<vt<<" c.f. actual flux: "<<flux[v];
      cout <<" , ratio = "<<vt/flux[v]<<" i="<<i<<"\n";i++;
    }
    flux[v] -= vt;
  }
  return 0;
}
#endif // LAPIDUS_VISCOSITY_ENABLED

int flux_solver_hydro_iso::AVFalle(const double *Pl, ///< Left Primitive state vector.
				   const double *Pr, ///< Right Primitive state vector.
				   const double *Pstar, ///< Resolved (P*) state vector.
				   double *flux, 
				   const double etav,
				   const double
				   )
{
  /** \section Equations
   * Currently the equations are as follows (for flux along x-direction):
   * \f[ F[p_x] = F[p_x] - \eta \rho^{*} a^{*} (v_{xR}-v_{xL})  \,,\f]
   * \f[ F[p_y] = F[p_y] - \eta \rho^{*} a^{*} (v_{yR}-v_{yL})  \,,\f]
   * \f[ F[p_z] = F[p_z] - \eta \rho^{*} a^{*} (v_{zR}-v_{zL})  \,,\f]
   * This follows from assuming a 1D problem, with non-zero bulk viscosity, and
   * a non-zero shear viscosity.  The bulk viscosity gives the viscosity 
   * in the x-direction, opposing stretching and compression.  The shear
   * viscosity gives the y and z terms, opposing slipping across boundaries.
   * 
   * The identification of the velocities and density with the
   * resultant state from the Riemann Solver is a bit ad-hoc, and could easily
   * be changed to using the mean of the left and right states, or something
   * else. (See my notes on artificial viscosity).
   * */
  ///
  /// The isothermal sound speed is the max-speed, and this is RS_pstar[eqAA].
  /// 
  double prefactor = Pstar[eqAA]*etav*Pstar[eqRO];
  flux[eqMMX] -= prefactor*(Pr[eqVX]-Pl[eqVX]);
  flux[eqMMY] -= prefactor*(Pr[eqVY]-Pl[eqVY]);
  flux[eqMMZ] -= prefactor*(Pr[eqVZ]-Pl[eqVZ]);
  return 0;
}


int flux_solver_hydro_iso::Roe_solver_hydro_iso()
{
  //
  // Roe's linear solver, from Toro 11.2.1 pp.346-349.
  //
  
  
  //
  // Now get Roe-averaged Velocity:
  // 
  double rrl = sqrt(RS_left[eqRO]);
  double rrr = sqrt(RS_right[eqRO]);
  RS_meanp[eqVX] = (rrl*RS_left[eqVX] +rrr*RS_right[eqVX])/(rrl+rrr);
  RS_meanp[eqAA] = 0.5*(RS_left[eqAA] +RS_right[eqAA]);
  
  //
  // And eigenvalues for this velocity
  // 
  double 
    lneg=RS_meanp[eqVX]-RS_meanp[eqAA],
    lpos=RS_meanp[eqVX]+RS_meanp[eqAA];
  
  //
  // If the negative e-value is positive, then the flux is just the left state flux.
  // Conversely if the positive e-value is negative, the solution is the right state.
  //
  if      (lneg >= 0.0) {  // if(u_av-c_av>0) pstar = left state
    PtoFlux(RS_left,FS_flux,0);
    return 0;
  }
  else if (lpos <=0.) {  //else if(u_av+c_av<=0) pstar = right state
    PtoFlux(RS_right,FS_flux,0);
    return 0;
  }
  else {  // else we are in the starred region, so we need to get the flux the long way
    //
    // First get the left state flux:
    // 
    PtoFlux(RS_left,FS_flux,0);
    
    //
    // Next add the contribution from the first wave
    // 
    double alpha1 = 0.5*((RS_right[eqRO]-RS_left[eqRO])*lpos - 
			 (RS_right[eqRO]*RS_right[eqVX]-RS_left[eqRO]*RS_left[eqVX]))/RS_meanp[eqAA];
    FS_flux[eqRHO] += alpha1*lneg;
    FS_flux[eqMMX] += alpha1*lneg*lneg;
    
    //
    // Need to assign values to v_y,v_z if present.
    // They only change across the contact discontinuity.
    //
    if (FS_flux[eqRHO]>0) {
      FS_flux[eqMMY] = FS_flux[eqRHO]*RS_left[eqVY];		    
      FS_flux[eqMMZ] = FS_flux[eqRHO]*RS_left[eqVZ];
      FS_flux[eqAAA] = FS_flux[eqRHO]*RS_left[eqAA];
    }
    else {
      FS_flux[eqMMY] = FS_flux[eqRHO]*RS_right[eqVY];		    
      FS_flux[eqMMZ] = FS_flux[eqRHO]*RS_right[eqVZ];
      FS_flux[eqAAA] = FS_flux[eqRHO]*RS_right[eqAA];
    }
#ifdef RSTESTING
    //
    // Now go back the other way, for testing...
    // 
    double temp[eq_nvar];
    PtoFlux(RS_right,temp,0);
    double alpha2 = 0.5*(-(RS_right[eqRO]-RS_left[eqRO])*lneg + 
			 (RS_right[eqRO]*RS_right[eqVX]-RS_left[eqRO]*RS_left[eqVX]))/RS_meanp[eqAA];
    temp[eqRHO] -= alpha2*lpos;
    temp[eqMMX] -= alpha2*lpos*lpos;
    
    if (!GS.equalD(temp[eqRHO],FS_flux[eqRHO]) ||
	!GS.equalD(temp[eqMMX],FS_flux[eqMMX])) {
      cout <<"rho: l*="<<FS_flux[eqRHO]<<"\tr*="<<temp[eqRHO]<<"\n";
      cout <<"p_x: l*="<<FS_flux[eqMMX]<<"\tr*="<<temp[eqMMX]<<"\n";
    }
#endif //RSTESTING
  }

  //
  // Calculate tracer flux here -- if mass flux is positive then
  // contact is at x>0 and we advect the left state tracer across
  // the boundary.  Otherwise we advect the right state to the left.
  // 
  if (FS_ntr>0) {
    if (FS_flux[eqRHO]>=0.0)
      for (int t=0;t<FS_ntr;t++)
	FS_flux[eqTR[t]] =  RS_left[eqTR[t]]*FS_flux[eqRHO];
    else 
      for (int t=0;t<FS_ntr;t++)
	FS_flux[eqTR[t]] = RS_right[eqTR[t]]*FS_flux[eqRHO];
  }
  
  //  cout << "(Isothermal Linear Riemann Solver) Success!" << "\n";
  return 0;
}

 
#endif // ISOTHERMAL_SOLVERS_ENABLED
