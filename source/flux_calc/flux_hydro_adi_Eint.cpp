///
/// \file flux_hydro_adi_Eint.cc
/// \author Jonathan Mackey
///
/// Declaration of the adiabatic hydrodynamics flux solver class,
/// which integrates the INTERNAL ENERGY and NOT the total energy.
///
/// Modifications:\n
///
/// - 2010.12.28 JM: Started file.
///
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.

#include "../defines/functionality_flags.h"
#ifdef INCLUDE_EINT_ADI_HYDRO

../#include "../global.h"
//#include "../equations/eqns_base.h"
//#include "../equations/eqns_hydro_adiabatic.h"
#include "../equations/eqns_hydro_adi_Eint.h"
//#include "flux_hydro_adiabatic.h"
#include "flux_hydro_adi_Eint.h"
using namespace std;

flux_solver_hydro_adi_Eint::flux_solver_hydro_adi_Eint(const int nv,
     ///< Number of variables in state vector
						       const double *state,
     ///< state vector which is 'typical' in the problem being solved. 
						       const double eta,
     ///< coefficient of artificial viscosity
						       const double g,
     ///< gamma (EOS).
						       const int ntr 
     ///< Number of tracer variables.
						       )
  : eqns_base(nv),
    flux_solver_base(nv,eta,ntr),
    eqns_Euler(nv),
    riemann_Euler(nv,state,g),
    Riemann_FVS_Euler(nv,g),
    Riemann_Roe_Hydro_PV(nv,g),
    Riemann_Roe_Hydro_CV(nv,g),
    flux_solver_hydro_adi(nv,state,eta,g,ntr),
    eqns_Euler_Eint(nv)
{
#ifdef FUNCTION_ID
  cout <<"::flux_solver_hydro_adi_Eint ...starting.\n";
#endif //FUNCTION_ID

  eq_gamma = g;
  cout <<"flux_solver_hydro_adi_Eint::flux_solver_hydro_adi_Eint() constructor\n";
  cout <<"Default solver set to "<<SimPM.solverType<<" where ";
  cout <<"0=LF,1=RSlin,2=RSexact,3=RShybrid,..,5=Roe-PV.\n";

#ifdef FUNCTION_ID
  cout <<"::flux_solver_hydro_adi_Eint ...returning.\n";
#endif //FUNCTION_ID
  return;
}

flux_solver_hydro_adi_Eint::~flux_solver_hydro_adi_Eint()
{
#ifdef FUNCTION_ID
  cout <<"::~flux_solver_hydro_adi_Eint ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"::~flux_solver_hydro_adi_Eint ...returning.\n";
#endif //FUNCTION_ID
  return;
}

int flux_solver_hydro_adi_Eint::inviscid_flux(const cell *Cl, ///< Left state cell pointer
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
  cout <<"flux_solver_hydro_adi_Eint::inviscid_flux ...starting.\n";
#endif //FUNCTION_ID

  //
  // Check input density and pressure are 'reasonably large'
  //
  if (Pl[eqRO]<TINYVALUE || Pl[eqPG]<TINYVALUE ||
      Pr[eqRO]<TINYVALUE || Pr[eqPG]<TINYVALUE) {
    rep.printVec("left ",Pl,eq_nvar);
    rep.printVec("right",Pr,eq_nvar);
    rep.error("flux_solver_hydro_adi_Eint::calculate_flux() Density/Pressure too small",
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
  cout <<"flux_solver_hydro_adi_Eint::inviscid_flux ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


int flux_solver_hydro_adi_Eint::AVFalle(const double *Pl,
				   const double *Pr,
				   const double *pstar,
				   double *flux, 
				   const double etav,
				   const double gam
				   )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::AVFalle ...starting.\n";
#endif //FUNCTION_ID
  /** \section Equations
   * Currently the equations are as follows (for flux along x-direction):
   * \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*} (v_{xR}-v_{xL})  \,,\f]
   * \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*} (v_{yR}-v_{yL})  \,,\f]
   * \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*} (v_{zR}-v_{zL})  \,,\f]
   *
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
#ifdef        EINT_ETOT_PARALLEL
  double prefactor = maxspeed(pstar, gam)*etav*pstar[eqRO];
  // x-dir first.
  double momvisc = prefactor*(Pr[eqVX]-Pl[eqVX]);
  double ergvisc = momvisc*pstar[eqVX];
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

#else  // not EINT_ETOT_PARALLEL
  double prefactor = maxspeed(pstar, gam)*etav*pstar[eqRO];
  // x-dir first.
  flux[eqMMX] -= prefactor*(Pr[eqVX]-Pl[eqVX]);
  // y-dir
  flux[eqMMY] -= prefactor*(Pr[eqVY]-Pl[eqVY]);
  // z-dir
  flux[eqMMZ] -= prefactor*(Pr[eqVZ]-Pl[eqVZ]);
  // No need to update energy flux b/c Internal erg is just advected.
#endif //     EINT_ETOT_PARALLEL

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::AVFalle ...returning.\n";
#endif //FUNCTION_ID
  return 0;
}

#ifdef LAPIDUS_VISCOSITY_ENABLED
int flux_solver_hydro_adi_Eint::AVLapidus(const cell *Cl, ///< Left state cell pointer
				     const cell *Cr, ///< Right state cell pointer
				     const double *Pl, ///< Left state
				     const double *Pr, ///< Right state
				     double *flux, 
				     const double etav,
				     const double gam
				     )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::AVLapidus( ...starting.\n";
#endif //FUNCTION_ID
  //
#ifdef TEST_LAPIDUS
  //

  rep.error("Write this properly so that it works!","LAPIDUS");

  double ul[5], ur[5]; // Euler equations only affect first 5 vars.
  PtoU(Pl,ul,gam); PtoU(Pr,ur,gam);
  //
  // I need a way to get div(v) for the edge point based on the two
  // cell-centred divergences... this is a bit of a hack at the mo.
  //
  double divu = 0.5*(CI.get_DivV(Cl,XX) + CI.get_DivV(Cr,XX));
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
#error "USE NEW TEST-LAPIDUS FLAG"
#endif // not TEST_LAPIDUS
  //

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::AVLapidus ...returning.\n";
#endif //FUNCTION_ID
  return(0);
}
#endif // LAPIDUS_VISCOSITY_ENABLED



void flux_solver_hydro_adi_Eint::PtoU(const double* p, double* u, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::PtoU ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler_Eint::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}

int flux_solver_hydro_adi_Eint::UtoP(const double *u, double *p, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::UtoP ...starting.\n";
#endif //FUNCTION_ID

  int err=eqns_Euler_Eint::UtoP(u,p,g);
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/p[eqRO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return err;
}

void flux_solver_hydro_adi_Eint::PUtoFlux(const double *p, const double *u, double *f)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::PUtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler_Eint::PUtoFlux(p,u,f);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = p[eqTR[t]]*f[eqRHO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::PUtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}

void flux_solver_hydro_adi_Eint::UtoFlux(const double *u, double *f, const double g)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::UtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_Euler_Eint::UtoFlux(u,f,g);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = u[eqTR[t]]*f[eqRHO]/u[eqRHO];

#ifdef FUNCTION_ID
  cout <<"flux_solver_hydro_adi_Eint::UtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}


#endif // if INCLUDE_EINT_ADI_HYDRO
