///
/// \file flux_hydro_adiabatic.cc
/// \author Jonathan Mackey
/// Function Definitions of the adiabatic hydrodynamics flux solver class
///
/// ------------------------------------------------------------------
///
/// Wish-List:
///
/// - I think Lapidus viscosity could be in the base flux solver.  Or
///   maybe, since it is the only function here with grid cell pointers,
///   it should be in the solver class built on top of this flux class?
///   Or I can re-define it so that each direction has a div(v) for its
///   own interface, like the H-correction has an etamax?  That way I 
///   can avoid any grid-awareness.
///
/// ------------------------------------------------------------------
///
/// History: \n
///  - 2009-10-20 Started on the file
///  - 2010-01-13 JM: Added Resolved state calculation to Roe Solver (Toro,1999,eq.11.17).
///  - 2010-07-31 JM: Added Primitive variable linear solver with Roe average.
///  - 2010-09-21 JM: Shortened some long lines.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.11.12 JM: Changed ->col to use cell interface for 
///   extra_data.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///    Renamed calculate_flux() to inviscid_flux() and moved AV
///    calculation to FV_solver_base class.
///  - 2010.11.19 JM: Added H-correction to Roe conserved Var solver.
/// - 2010.11.21 JM: ifdef option to use the mean state for Pstar in
///   Roe conserved variable solver (4).  This is only relevant for
///   the FKJ98 viscosity function, and make sure that pstar has
///   positive definite pressure and density.  The ifdef statement is
///   at the top of this file.
/// - 2010.11.21 JM: Added symmetric version of the Roe conserved
///   variables flux solver, called Roe_flux_solver_symmetric().  It
///   is MUCH BETTER than the one-sided one, which I have renamed to
///   Roe_flux_solver_onesided().  For an axi-symmetric blast wave
///   with the H-correction, the one-sided solver developed spikes at
///   90 degree intervals, but the symmetric solver is perfectly
///   clean.  SO, DON'T USE ONE-SIDED ROE FLUX SOLVER ANYMORE!
/// - 2010.12.22 JM: Moved Roe PV and CV solvers to Riemann_solvers/
///   directory.  Got rid of the ROE_CV_MEANP ifdef (it is default).
/// - 2010.12.23 JM: Moved UtoP() etc. from solver to flux-solver.
///   Now all tracer stuff is dealt with here.  Added pstar[] array
///   to inviscid_flux() function since it is no longer inherited.
/// - 2010.12.27 JM: Enclosed Lapidus AV in an ifdef.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "flux_hydro_adiabatic.h"
using namespace std;

flux_solver_hydro_adi::flux_solver_hydro_adi(const int nv,        ///< Number of variables in state vector
					     const double *state, ///< state vector which is 'typical' in the problem being solved. 
					     const double eta,     ///< coefficient of artificial viscosity
					     const double g,       ///< gamma (EOS).
					     const int ntr         ///< Number of tracer variables.
					     )
  : eqns_base(nv),
    // riemann_base(nv),
    flux_solver_base(nv,eta,ntr),
    eqns_Euler(nv),
    riemann_Euler(nv,state,g),
    Riemann_FVS_Euler(nv,g),
    Riemann_Roe_Hydro_PV(nv,g),
    Riemann_Roe_Hydro_CV(nv,g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::flux_solver_hydro_adi ...starting.\n";
#endif //FUNCTION_ID

  eq_gamma = g;
#ifdef TESTING
  cout <<"flux_solver_hydro_adi::flux_solver_hydro_adi() constructor: gamma="<<eq_gamma<<"\n";
  cout <<"Default solver set to "<<SimPM.solverType<<" where ";
  cout <<"0=LF,1=RSlin,2=RSexact,3=RShybrid,4=RSRoe,5=Roe-PV,6=FVS.\n";
#endif

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::flux_solver_hydro_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}

flux_solver_hydro_adi::~flux_solver_hydro_adi()
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::~flux_solver_hydro_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::~flux_solver_hydro_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


int flux_solver_hydro_adi::inviscid_flux(const cell *Cl, ///< Left state cell pointer
					 const cell *Cr, ///< Right state cell pointer
					 const double *Pl,///< Left Primitive state vector.
					 const double *Pr,///< Right Primitive state vector.
					 double *flux,   ///< Resultant Flux state vector.
					 double *pstar,  ///< State vector at interface.
					 const int solve_flag,
	       ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
					 const double g  ///< Gas constant gamma.
					 )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::inviscid_flux ...starting.\n";
#endif //FUNCTION_ID

  //
  // Check input density and pressure are 'reasonably large'
  //
  if (Pl[eqRO]<TINYVALUE || Pl[eqPG]<TINYVALUE ||
      Pr[eqRO]<TINYVALUE || Pr[eqPG]<TINYVALUE) {
    rep.printVec("left ",Pl,eq_nvar);
    rep.printVec("right",Pr,eq_nvar);
    rep.error("flux_solver_hydro_adi::calculate_flux() Density/Pressure too small",
	      Pl[eqRO]);
  }
  int err=0;

  //
  // Set flux and pstar vector to zero.
  //
  for (int v=0;v<eq_nvar;v++) flux[v]  = 0.0;
  for (int v=0;v<eq_nvar;v++) pstar[v]  = 0.0;
  //
  // Set EOS gamma in riemann solver class:
  //
  eq_gamma = g;


  //
  // Choose which Solver to use.  Each of these must set the values of
  // the flux[] and pstar[] arrays.  I've commented out the calls to
  // set_interface_tracer_flux(Pl,Pr,flux) b/c this is called by the
  // solver which calls this function, after applying any requested
  // viscous corrections.  So it is redundant work to calculate it
  // now.
  //
  if      (solve_flag==FLUX_LF) {
    //
    // Lax-Friedrichs Method, so just get the flux.
    //
    //cout <<"using LF flux\n";
    err += get_LaxFriedrichs_flux(Pl,Pr,flux,eq_gamma);
    //set_interface_tracer_flux(Pl,Pr,flux);
    for (int v=0;v<eq_nvar;v++) pstar[v] = 0.5*(Pl[v]+Pr[v]);
    // pstar? yes.  flux? yes.
 }

  else if (solve_flag==FLUX_FVS) {
    //
    // Flux Vector Splitting (van Leer, 1982): This function takes the
    // left and right state and returns the FVS flux in 'flux' and the
    // Roe-average state in 'pstar'.
    //
    err += FVS_flux(Pl,Pr, flux, pstar, eq_gamma);
    //set_interface_tracer_flux(Pl,Pr,flux);
    // pstar? yes.  flux? yes.
  }

  else if (solve_flag==FLUX_RSlinear ||
	   solve_flag==FLUX_RSexact ||
	   solve_flag==FLUX_RShybrid) {
    //
    // These are all Riemann Solver Methods, so call the solver:
    //
    //cout <<"using Riemann Solver Flux: flag="<<solve_flag<<"\n";
    err += JMs_riemann_solve(Pl,Pr,pstar,solve_flag,eq_gamma);
    //
    // Convert pstar to a flux.
    PtoFlux(pstar, flux, eq_gamma);
    //set_interface_tracer_flux(Pl,Pr,flux);
    // pstar? yes.  flux? yes.
  }

  else if (solve_flag==FLUX_RSroe) {
    //
    // Roe conserved variables flux solver (Toro 1999), either the
    // symmetric or one-sided calculation (Symmetric is dramatically
    // better -- the asymmetric one is just there for debugging, or in
    // case i need it for something in the future!).
    //
    //    err += Roe_flux_solver_onesided(Pl,Pr,eq_gamma,
    //                                   HC_etamax,
    //				         pstar,FS_flu
    //				         ); // DON'T USE ASYMMETRIC VERSION!!!
    err += Roe_flux_solver_symmetric(Pl,Pr,eq_gamma,
				     HC_etamax,
				     pstar,flux);
    //set_interface_tracer_flux(Pl,Pr,flux);

    //#define RoeHYDRO_TESTING
#ifdef RoeHYDRO_TESTING
    err += riemann_solve(Pl,Pr,pstar,1,eq_gamma);
    double temp[eq_nvar], diff=0.0;
    PtoFlux(pstar, temp, eq_gamma);
    for (int v=0;v<5;v++)
      diff += fabs(temp[v]-flux[v])/
	(fabs(temp[v])+fabs(flux[v])+TINYVALUE);
    if (diff >1.e-5) {
      cout <<"diff="<<diff<<"\n";
      rep.printVec("Roe flux : ",flux ,eq_nvar);
      rep.printVec("Old flux : ",temp ,eq_nvar);
    }
#endif
    // pstar? yes.  flux? yes.
  }

  else if (solve_flag==FLUX_RSroe_pv) {
    //
    // Roe primitive variables linear solver:
    //
    //cout <<"using Riemann Solver Flux: flag="<<solve_flag<<"\n";
    err += Roe_prim_var_solver(Pl,Pr,eq_gamma,pstar);
    //
    // Convert pstar to a flux:
    //
    PtoFlux(pstar, flux, eq_gamma);
    //set_interface_tracer_flux(Pl,Pr,flux);
    // pstar? yes.  flux? yes.
  }

  else {
    rep.error("what sort of flux solver do you mean???",solve_flag);
  }

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::inviscid_flux ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


int flux_solver_hydro_adi::AVFalle(const double *Pl,
				   const double *Pr,
				   const double *pstar,
				   double *flux, 
				   const double etav,
				   const double gam
				   )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::AVFalle ...starting.\n";
#endif //FUNCTION_ID
  /** \section Directionality
   * This function assumes the boundary normal direction has been set.
   * */
  /** \section Equations
   * Currently the equations are as follows (for flux along x-direction):
   * \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*}_f (v_{xR}-v_{xL})  \,,\f]
   * \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*}_f (v_{yR}-v_{yL})  \,,\f]
   * \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*}_f (v_{zR}-v_{zL})  \,,\f]
   * \f[ F[e]   = F[e]   - \eta \rho^{*} c^{*}_f \left( v^{*}_x (v_{xR}-v_{xL}) + v^{*}_y (v_{yR}-v_{yL}) + v^{*}_z (v_{zR}-v_{zL}) \right)\,.\f]
   * This follows from assuming a 1D problem, with non-zero bulk viscosity, and
   * a non-zero shear viscosity.  The bulk viscosity gives the viscosity 
   * in the x-direction, opposing stretching and compression.  The shear
   * viscosity given the y and z terms, opposing slipping across boundaries.
   * 
   * The identification of the velocities and density with the
   * resultant state from the Riemann Solver is a bit ad-hoc, and could easily
   * be changed to using the mean of the left and right states, or something
   * else. (See my notes on artificial viscosity).
   * */
  double prefactor = maxspeed(pstar, gam)*etav*pstar[eqRO];

  //
  // x-dir first.
  //
  double momvisc = prefactor*(Pr[eqVX]-Pl[eqVX]);//*4./3.;
  double ergvisc = momvisc*pstar[eqVX];
  //cout <<"falle AV: eta="<<etav<<" momvisc="<<momvisc<<"\n";
  flux[eqMMX] -= momvisc;

  // y-dir
  momvisc = prefactor*(Pr[eqVY]-Pl[eqVY]);
  flux[eqMMY] -= momvisc;
  ergvisc += momvisc*pstar[eqVY];

  // z-dir
  momvisc = prefactor*(Pr[eqVZ]-Pl[eqVZ]);
  flux[eqMMZ] -= momvisc;
  ergvisc += momvisc*pstar[eqVZ];

  // update energy flux
  flux[eqERG] -= ergvisc;

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::AVFalle ...returning.\n";
#endif //FUNCTION_ID
  return 0;
}

#ifdef LAPIDUS_VISCOSITY_ENABLED
int flux_solver_hydro_adi::AVLapidus(const cell *Cl, ///< Left state cell pointer
				     const cell *Cr, ///< Right state cell pointer
				     const double *Pl, ///< Left state
				     const double *Pr, ///< Right state
				     double *flux, 
				     const double etav,
				     const double gam
				     )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::AVLapidus( ...starting.\n";
#endif //FUNCTION_ID
  //
#ifdef TEST_LAPIDUS
  //

  double ul[5], ur[5]; // Euler equations only affect first 5 vars.
  PtoU(Pl,ul,gam); PtoU(Pr,ur,gam);
  //
  // I need a way to get div(v) for the edge point based on the two
  // cell-centred divergences... this is a bit of a hack at the mo.
  //
  double divu = 0.5*(CI.get_DivV(Cl) + CI.get_DivV(Cr));
  //double divu = 0.5*(fabs(CI.get_DivV(Cl,XX)) + fabs(CI.get_DivV(Cr,XX)));
  //double divu = max(fabs(CI.get_DivV(Cl,XX)),fabs(CI.get_DivV(Cr,XX)));
  for (int v=0;v<5;v++) {
    //
    // When it has been tested, this is more efficient:
    // This is the CW84 formula, but it only adds viscosity in shocks.
    //
    flux[v] -= etav*max(-divu,0.0)*(ur[v]-ul[v]);
    //flux[v] -= etav*divu*(ur[v]-ul[v]);
    //    
    // vt = etav*divu*(ur[v]-ul[v]);
    // if (fabs(vt)>= fabs(flux[v]) && fabs(flux[v])>BASEPG) {
    //   cout <<"\t\t\t\t\tlarge visc: var="<<v<<" vt= "<<vt<<" c.f. actual flux: "<<flux[v];
    //   cout <<" , ratio = "<<vt/flux[v]<<" i="<<i<<"\n";i++;
    // }
    // flux[v] -= vt;
  }

  //
#else // instead of TEST_LAPIDUS
  //

  ///
  /// This is not working!!! Don't use it!
  /// THIS USES A VERY BAD DIV(V) APPROXIMATION -- REALLY NEED TO CALCULATE IT
  /// PROPERLY, SO MAYBE THE SOLVER CLASS CAN DO THIS???
  ///
  rep.error("Don't use Lapidus!!!",1);

  double ul[5], ur[5]; // Euler equations only affect first 5 vars.
  PtoU(Pl,ul,gam); PtoU(Pr,ur,gam);
  double divu = ((Pr[eqVX]-Pl[eqVX])+
		 (Pr[eqVY]-Pl[eqVY])+
		 (Pr[eqVZ]-Pl[eqVZ]));
  double vt; int i=0;
  for (int v=0;v<5;v++) {
    vt = etav*divu*(ur[v]-ul[v]);
    if (fabs(vt)>= fabs(flux[v]) && fabs(flux[v])>BASEPG) {
      cout <<"\t\t\t\t\tlarge visc: var="<<v<<" vt= "<<vt<<" c.f. actual flux: "<<flux[v];
      cout <<" , ratio = "<<vt/flux[v]<<" i="<<i<<"\n";i++;
    }
    flux[v] -= vt;
  }
  //
#endif // not TEST_LAPIDUS
  //

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::AVLapidus ...returning.\n";
#endif //FUNCTION_ID
  return(0);
}
#endif // LAPIDUS_VISCOSITY_ENABLED



void flux_solver_hydro_adi::PtoU(const double* p, double* u, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::PtoU ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}

int flux_solver_hydro_adi::UtoP(const double *u, double *p, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::UtoP ...starting.\n";
#endif //FUNCTION_ID

  int err=eqns_Euler::UtoP(u,p,g);
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/p[eqRO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return err;
}

void flux_solver_hydro_adi::PUtoFlux(const double *p, const double *u, double *f)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::PUtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler::PUtoFlux(p,u,f);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = p[eqTR[t]]*f[eqRHO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::PUtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}

void flux_solver_hydro_adi::UtoFlux(const double *u, double *f, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::UtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler::UtoFlux(u,f,g);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = u[eqTR[t]]*f[eqRHO]/u[eqRHO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi::UtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}


