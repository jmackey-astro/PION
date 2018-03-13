///
/// \file flux_mhd_adiabatic.cc
/// \author Jonathan Mackey
/// Function Definitions of the adiabatic hydrodynamics flux solver class
///
/// History: 
///  - 2009-10-20 Started on the file
///  - 2009-10-24 Got it working.
///  - 2010-01-15 JM: Added in calculation of Pstar for Roe Flux solver
///  - 2010-01-18 JM: Fixed error in Roe Solver (Pstar not calculated when star state
///      is the left or right state).
///  - 2010-02-19 JM: Today and yesterday I changed some of the error reporting in
///     the Roe solver.  Now it should only complain about the first 1000 negative 
///     pressures and 1000 negative densities in the starred state.  Also the Roe solver()
///     will return 1 on negative pressure and 2 on negative density, so calculate_flux()
///     can then decide what to do about the solution.  Currently if an error is returned
///     I try Sam Falle's solver and use that solution regardless.
///  - 2010-02-19 JM: Removed pointless PtoU()->UtoP() calculation to get Pstar when
///     the resolved state is either the left or right state.
///  - 2010-02-20 JM: Set whether to use Ustar or meanP for the Pstar value to return to 
///     be in an ifdef.  Ustar is not great for stability, and I think meanP will be much
///     better, as well as quicker to calculate.  So I have it unset now.
///  - 2010-09-02 JM: Commented out unnecessary initialisation of
///     variables in the Roe solver.  Identified ways to speed up the
///     eigenvector calculation by avoiding doubling the calculation
///     for the left and right moving waves.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///    Renamed calculate_flux() to inviscid_flux() and moved AV
///    calculation to FV_solver_base class.
/// - 2010.12.07 JM: AVFalle() now uses the passed--in pointer to the
///   flux vector.  I need to get rid of the class variables -- they 
///   are just a recipe for disaster (AV was having no effect at all!)
/// - 2010.12.27 JM: Moved Roe flux solver to own class in Riemann_solvers/
///   Got rid of inherited class flux/left/right/pstar variables.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include "flux_mhd_adiabatic.h"
using namespace std;


// ##################################################################
// ##################################################################



flux_solver_mhd_ideal_adi::flux_solver_mhd_ideal_adi(
       const int nv, ///< Number of variables in state vector
       const pion_flt *state, ///< state vector which is 'typical' in the problem being solved. 
       const double eta, ///< coefficient of artificial viscosity
       const double g, ///< gamma (EOS).
       const int ntr ///< Number of tracer variables.
       )
: eqns_base(nv),
    //    riemann_base(nv),
    flux_solver_base(nv,eta,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,g),
    Riemann_Roe_MHD_CV(nv,g)
{
  negPGct = negROct = 0;
  return;
}



// ##################################################################
// ##################################################################



flux_solver_mhd_ideal_adi::~flux_solver_mhd_ideal_adi()
{
#ifdef TESTING
  cout <<"flux_solver_mhd_ideal_adi::flux_solver_mhd_ideal_adi()  destructor.\n";
#endif
}


// ##################################################################
// ##################################################################



int flux_solver_mhd_ideal_adi::inviscid_flux(
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
    rep.error("flux_solver_mhd_ideal_adi::calculate_flux() Density/Pressure too small",Pl[eqRO]);
  }
#endif //TESTING
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
  // Choose which Solver to use:
  //
  if      (solve_flag==FLUX_LF) {
    //
    // Lax-Friedrichs Method, so just get the flux
    //
    err += get_LaxFriedrichs_flux(Pl,Pr,flux,eq_gamma);
    //set_interface_tracer_flux(Pl,Pr,flux);
    // This is not really needed for LF b/c shouldn't need to add
    // viscosity.
    for (int v=0;v<eq_nvar;v++) pstar[v] = 0.5*(Pl[v]+Pr[v]);
    // pstar? yes.  flux? yes.
  }

  else if (solve_flag==FLUX_RSroe) {
    //
    // Roe Flux solver in conserved variables (Cargo and Gallice, 1997).
    // This is the symmetric version which sums over all waves.
    //
    err += MHD_Roe_CV_flux_solver_symmetric(Pl,Pr,eq_gamma,
					    HC_etamax,
					    pstar,flux);

    // The one-sided solver (not recommended).
    //     err += MHD_Roe_CV_flux_solver_onesided(Pl,Pr,eq_gamma,
    // 				     HC_etamax,
    // 				     pstar,flux);

#ifdef RoeMHD_TESTING
    err += JMs_riemann_solve(Pl,Pr,pstar,1,eq_gamma);
    rep.printVec("Roe flux : ",flux ,eq_nvar);
    PtoFlux(pstar, flux, eq_gamma);
    rep.printVec("Old flux : ",flux ,eq_nvar);
#endif
    //
    // If we get an error, try the linear solver:
    //
    if (err) {
      //err=0;
      cout <<"ROE SOLVER FAILED -- TRYING FKJ98 SOLVER. err="<<err<<"\n";
      err = JMs_riemann_solve(Pl,Pr,pstar,1,eq_gamma);
      //
      // convert pstar to a flux:
      //
      PtoFlux(pstar, flux, eq_gamma);
    }
  }

  else if (solve_flag==FLUX_RSlinear ||
	   solve_flag==FLUX_RSexact  ||
	   solve_flag==FLUX_RShybrid) {
    
#ifdef TESTING
# ifdef RS_TESTING
    //rep.printVec("pstar: ",pstar,eq_nvar);
    //rep.printVec("left : ",Pl ,eq_nvar);
    //rep.printVec("right: ",Pr,eq_nvar);
    pion_flt l[eq_nvar], r[eq_nvar];
    for (int v=0;v<eq_nvar;v++) {
      l[v] = Pl[v];
      r[v] = Pr[v];
    }
    err += JMs_riemann_solve(l,r,pstar,solve_flag,eq_gamma);
    double diff=0;
    for (int v=0;v<7;v++){
      diff += fabs(Pl[v] -l[v]);
      diff += fabs(Pr[v]-r[v]);
    }
    if (!pconst.equalD(diff,0.0)) { 
      cout <<"************ dir="<<eq_dir<<" *************************************\n";
      rep.printVec("pstar: ",pstar,eq_nvar);
      rep.printVec("left : ",Pl ,eq_nvar);
      rep.printVec("right: ",Pr,eq_nvar);
      rep.printVec("templ: ",l ,eq_nvar);
      rep.printVec("tempr: ",r,eq_nvar);
      cout <<"*************************************************\n\n";
    }
# endif // RS_TESTING
#endif // TESTING
    
    //
    // JM's hybrid Riemann Solver -- Falle et al. (1998) with the Roe
    // and Balsara (1996) eigenvector normalisation.
    //
    err += JMs_riemann_solve(Pl,Pr,pstar,solve_flag,eq_gamma);
    //    rep.printVec("pstar: ",pstar,eq_nvar);
    //
    // Convert pstar to a flux:
    //
    PtoFlux(pstar, flux, eq_gamma);
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



int flux_solver_mhd_ideal_adi::AVFalle(
      const pion_flt *Pleft,
      const pion_flt *Pright,
      const pion_flt *Pstar,
      pion_flt *flux, 
      const double eta, ///< already set as FS_etav
      const double gamma ///< already set as eq_gamma
      )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_ideal_adi::AVFalle ...starting.\n";
#endif //FUNCTION_ID
  /** \section Directionality
   * This function assumes the boundary normal direction has been set.
   * */
  /** \section Equations
   * The equations are as follows (for flux along x-direction):
   * \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*}_f (v_{xR}-v_{xL})  \,,\f]
   * \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*}_f (v_{yR}-v_{yL})  \,,\f]
   * \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*}_f (v_{zR}-v_{zL})  \,,\f]
   * \f[ F[e]   = F[e]   - \eta \rho^{*} c^{*}_f \left[ v^{*}_x (v_{xR}-v_{xL}) + v^{*}_y (v_{yR}-v_{yL}) + v^{*}_z (v_{zR}-v_{zL}) \right]\,.\f]
   * This follows from assuming a 1D problem, with non-zero bulk viscosity, and
   * a non-zero shear viscosity.  The bulk viscosity gives the viscosity 
   * in the x-direction, opposing stretching and compression.  The shear
   * viscosity gives the y and z terms, opposing slipping across boundaries.
   * 
   * The identification of the velocities and density with the
   * resultant state from the Riemann Solver is a bit ad-hoc, and could easily
   * be changed to using the mean of the left and right states, or something
   * else. (See my notes on artificial viscosity).
   * 
   * Andy says Sam also has similar terms for the Magnetic field, although 
   * it isn't in his paper.  Andy uses this too.
   * \f[ f[B_y] = f[B_y] - \eta c^{*}_f (B_{yR}-B_{yL}) \;, \f]
   * \f[ f[B_z] = f[B_z] - \eta c^{*}_f (B_{zR}-B_{zL}) \;, \f]
   * with associated terms added to energy flux:
   * \f[ F[e]   = F[e]   - \eta c^{*}_f \left[ B^{*}_z (B_{zR}-B_{zL}) +B^{*}_y (B_{yR}-B_{yL})\right] \;. \f]
   * 
   * Andy uses the mean of the left and right state velocities/field to calculate
   * the energy flux; I am using the starred state velocity/field.  Not sure
   * if that is bad, or if it will make any discernible difference at all.
   * After some testing, it makes very little difference.  It may help with
   * some crazy cases, but I haven't encountered them yet.
   * 
   * I have switched from using \f$c^{*}_f\f$ to using the fast speed from 
   * the mean vector of the left and right states.  This is because the 
   * Riemann Solver can return negative density, which then gets set to a 
   * density floor value, but the B-field is still normal, so the fast speed
   * blows up.  If I use the mean state though, it represents a typical speed
   * between the two states, so it is in principle no less realistic! (says he
   * hopefully...)
   * */

//  rep.printVec(" left",Pleft,eq_nvar);
//  rep.printVec("right",Pright,eq_nvar);
//  rep.printVec("pstar",Pstar,eq_nvar);
//  cout <<"FKJ98 Artificial viscosity for MHD!\n";
  double prefactor = cfast_components(0.5*(Pleft[eqRO]+Pright[eqRO]),
				      0.5*(Pleft[eqPG]+Pright[eqPG]),
				      0.5*(Pleft[eqBX]+Pright[eqBX]),
				      0.5*(Pleft[eqBY]+Pright[eqBY]),
				      0.5*(Pleft[eqBZ]+Pright[eqBZ]),
				      eq_gamma)*FS_etav*Pstar[eqRO];
				      
//  double prefactor = cfast(pstar, eq_gamma)*etav*Pstar[eqRO]; // This causes trouble at low density in resolved state.
//  double prefactor = Pstar[eqVX]*etav*Pstar[eqRO]; // This is stupid... no shear viscosity!!!

  // Momentum flux
  double momvisc = prefactor*(Pright[eqVX]-Pleft[eqVX]);//*4./3.;
  double ergvisc = momvisc*Pstar[eqVX]; // This makes very little difference as far as I can tell (compared to line below).
  //  double ergvisc = momvisc*0.5*(Pleft[eqVX]+Pright[eqVX]);
#ifdef AVFALLE_TESTING
  int i=0;
  if (fabs(momvisc)>= fabs(flux[eqMMX]) && fabs(flux[eqMMX])>BASEPG) {
    cout <<"\tlarge MMX visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqMMX];
    cout <<" , ratio = "<<momvisc/flux[eqMMX]<<" i="<<i<<"\n";i++;
  }
#endif
  flux[eqMMX] -= momvisc;
  
  momvisc = prefactor*(Pright[eqVY]-Pleft[eqVY]);
#ifdef AVFALLE_TESTING
  if (fabs(momvisc)>= fabs(flux[eqMMY]) && fabs(flux[eqMMY])>BASEPG) {
    cout <<"\t\tlarge MMY visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqMMY];
    cout <<" , ratio = "<<momvisc/flux[eqMMY]<<" i="<<i<<"\n";i++;
 }
#endif
  flux[eqMMY] -= momvisc;
  ergvisc += momvisc*Pstar[eqVY];
//  ergvisc += momvisc*0.5*(Pleft[eqVY]+Pright[eqVY]);
  
  momvisc = prefactor*(Pright[eqVZ]-Pleft[eqVZ]);
#ifdef AVFALLE_TESTING
  if (fabs(momvisc)>= fabs(flux[eqMMZ]) && fabs(flux[eqMMZ])>BASEPG) {
    cout <<"\t\t\tlarge MMZ visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqMMZ];
    cout <<" , ratio = "<<momvisc/flux[eqMMZ]<<" i="<<i<<"\n";i++;
  }
#endif
  flux[eqMMZ] -= momvisc;
  ergvisc += momvisc*Pstar[eqVZ];
//  ergvisc += momvisc*0.5*(Pleft[eqVZ]+Pright[eqVZ]);

  //
  // Magnetic field flux.
  //
  prefactor *= FS_etaB/(FS_etav*Pstar[eqRO]);
  momvisc = prefactor*(Pright[eqBY]-Pleft[eqBY]);
#ifdef AVFALLE_TESTING
  if (fabs(momvisc)>= fabs(flux[eqBBY]) && fabs(flux[eqBBY])>1) {
    cout <<"\t\t\t\t\tlarge BBY visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqBBY];
    cout <<" , ratio = "<<momvisc/flux[eqBBY]<<" i="<<i<<"\n";i++;
    cout <<"\t\t\t\t\t\tPrefactor: "<<prefactor<<" and Br: "<<Pright[eqBY]<<" , Bl: "<<Pleft[eqBY]<<"\n";
    cout <<"\t\t\t\t\t\tcf: "<<cfast(Pstar, eq_gamma)<<" bx: "<<Pstar[eqBX]<<" by: "<<Pstar[eqBY]<<" bz: "<<Pstar[eqBZ]<<" pg: "<<Pstar[eqPG]<<" ro: "<<Pstar[eqRO]<<"\n";
  }
#endif

  flux[eqBBY] -= momvisc;
  ergvisc += momvisc*Pstar[eqBY];
  
  momvisc = prefactor*(Pright[eqBZ]-Pleft[eqBZ]);
#ifdef AVFALLE_TESTING
  if (fabs(momvisc)>= fabs(flux[eqBBZ]) && fabs(flux[eqBBZ])>1) {
    cout <<"\t\t\t\t\tlarge BBZ visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqBBZ];
    cout <<" , ratio = "<<momvisc/flux[eqBBZ]<<" i="<<i<<"\n";i++;
    cout <<"\t\t\t\t\t\tPrefactor: "<<prefactor<<" and Br: "<<Pright[eqBZ]<<" , Bl: "<<Pleft[eqBZ]<<"\n";
    cout <<"\t\t\t\t\t\tcf: "<<cfast(Pstar, eq_gamma)<<" bx: "<<Pstar[eqBX]<<" by: "<<Pstar[eqBY]<<" bz: "<<Pstar[eqBZ]<<" pg: "<<Pstar[eqPG]<<" ro: "<<Pstar[eqRO]<<"\n";
  }
#endif

  flux[eqBBZ] -= momvisc;
  ergvisc += momvisc*Pstar[eqBZ];

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!! EXPERIMENTAL CODE -- REMOVE THIS FLAG IF KEEPING !!!!!!!!!  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //#define AVFALLE_BXVISC
#ifdef AVFALLE_BXVISC
 cout <<"using BXVISC\n";
  momvisc = prefactor*(Pright[eqBX]-Pleft[eqBX]);
#ifdef AVFALLE_TESTING
  if (fabs(momvisc)>= fabs(flux[eqBBX]) && fabs(flux[eqBBX])>1) {
    cout <<"\t\t\t\t\tlarge BBX visc: "<<momvisc<<" c.f. actual flux: "<<flux[eqBBX];
    cout <<" , ratio = "<<momvisc/flux[eqBBX]<<" i="<<0<<"\n";
    cout <<"\t\t\t\t\t\tPrefactor: "<<prefactor<<" and Br: "<<Pright[eqBX]<<" , Bl: "<<Pleft[eqBX]<<"\n";
    cout <<"\t\t\t\t\t\tcf: "<<cfast(Pstar, eq_gamma)<<" bx: "<<Pstar[eqBX]<<" by: "<<Pstar[eqBY]<<" bz: "<<Pstar[eqBX]<<" pg: "<<Pstar[eqPG]<<" ro: "<<Pstar[eqRO]<<"\n";
  }
#endif
  flux[eqBBX] -= momvisc;
  ergvisc += momvisc*Pstar[eqBX];
#endif
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!! EXPERIMENTAL CODE -- REMOVE THIS FLAG IF KEEPING !!!!!!!!!  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Energy flux.
#ifdef AVFALLE_TESTING
  if (fabs(ergvisc)>= fabs(flux[eqERG]) && fabs(flux[eqERG])>1) {
    cout <<"large ERG visc: "<<ergvisc<<" c.f. actual flux: "<<flux[eqERG];
    cout <<" , ratio = "<<ergvisc/flux[eqERG]<<" i="<<i<<"\n";i++;
  }
#endif

  flux[eqERG] -= ergvisc;

#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_ideal_adi::AVFalle ...returning.\n";
#endif //FUNCTION_ID
  return(0);  
}


// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################





// **********************************************************************************
// flux_solver_mhd_mixedGLM_adi class, for the Dedner-GLM divergence cleaning method.
// **********************************************************************************

flux_solver_mhd_mixedGLM_adi::flux_solver_mhd_mixedGLM_adi(
      const int nv,
      ///< Number of variables in state vector
      const pion_flt *state,
      ///< state vector which is 'typical' in the problem being solved. 
      const double eta,
      ///< coefficient of artificial viscosity
      const double g,
      ///< gamma (EOS).
      const int ntr
      ///< number of tracer variables.
      )
  : eqns_base(nv),
    //    riemann_base(nv),
    flux_solver_base(nv,eta, ntr), 
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,g),
    Riemann_Roe_MHD_CV(nv,g),
    flux_solver_mhd_ideal_adi(nv,state,eta,g,ntr),
    eqns_mhd_mixedGLM(nv)
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::flux_solver_mhd_mixedGLM_adi() ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::flux_solver_mhd_mixedGLM_adi() ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################




flux_solver_mhd_mixedGLM_adi::~flux_solver_mhd_mixedGLM_adi()
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::~flux_solver_mhd_mixedGLM_adi() ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::~flux_solver_mhd_mixedGLM_adi() ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################



int flux_solver_mhd_mixedGLM_adi::inviscid_flux(
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
#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::inviscid_flux ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  //
  // Check input density and pressure are 'reasonably large'
  //
  if (Pl[eqRO]<TINYVALUE || Pl[eqPG]<TINYVALUE ||
      Pr[eqRO]<TINYVALUE || Pr[eqPG]<TINYVALUE) {
    rep.printVec("left ",Pl,eq_nvar);
    rep.printVec("right",Pr,eq_nvar);
    rep.error("flux_solver_mhd_mixedGLM_adi::calculate_flux() Density/Pressure too small",Pl[eqRO]);
  }
#endif //TESTING
  
  int err=0;

  //
  // Set flux and pstar vector to zero.
  //
  for (int v=0;v<eq_nvar;v++) flux[v]  = 0.0;
  for (int v=0;v<eq_nvar;v++) pstar[v]  = 0.0;
  //
  // Need temporary left and right state vectors b/c we need to change
  // the left and right state values of eqBX to the resolved state of
  // the Dedner et al. (2002) Riemann problem for (Bx,Psi).
  //
  pion_flt left[eq_nvar], right[eq_nvar];
  for (int v=0;v<eq_nvar;v++) left[v]  = Pl[v];
  for (int v=0;v<eq_nvar;v++) right[v] = Pr[v];

  //
  // Set EOS gamma in riemann solver class:
  //
  eq_gamma = g;
  
  //
  // First get the resolved state in the 2x2 Dedner Riemann Problem,
  // which is described well in their paper, and/or my thesis.
  //
  // Uses Dedner eq.41 for the flux in Bx and Psi:
  // \f[ \partial_t B_x + \partial_x \psi = 0 \;, \qquad 
  //     \partial_t \psi + \partial_x (c_h^2 B_x) = 0 \;. \f]
  // where the source term has been omitted, as it is calculated separately.
  // 
  // The GLM method has Bx and Psi decoupled from all the other variables
  // in the Riemann Problem, so they can be solved separately as a two
  // variable system (Dedner eq.42)
  // 
  // \f[ F(\psi) = c_h^2 B_x^* = c_h^2 \left( \frac{1}{2}(B_x(L)+B_x(R))
  //                             - \frac{1}{2c_h}(\psi_R-\psi_L) \right) \f]
  // \f[ F(B_x) = \psi_* = \frac{1}{2}(\psi_L+\psi_R) - \frac{c_h}{2}(B_x(R)-B_X(L)) \f]
  // 
  // Bx(*) is then used for both the left and right state to calculate
  // the flux.
  // 
  double psistar = 0.5*(left[eqSI]+right[eqSI]
			-GLM_chyp*(right[eqBX]-left[eqBX]));
  double bxstar  = 0.5*(left[eqBX]+right[eqBX]
			-(right[eqSI]-left[eqSI])/GLM_chyp);
  left[eqBX] = right[eqBX] = bxstar;

  //
  // Now continue on in an identical manner as the ideal MHD solver, 
  // by calling its flux solver:
  //
  err=flux_solver_mhd_ideal_adi::inviscid_flux(Cl,Cr,left,right,
					       flux,pstar,
					       solve_flag,eq_gamma);


  //
  // Now we have to add in the flux in BX and PSI, based on Dedner et
  // al.'s method.
  //
  // NOTE: Dedner doesn't say what to do about the energy flux, so this is my
  // best guess to ensure consistency.  Andy does something very similar, if
  // not identical.
  // 
  flux[eqPSI]  = GLM_chyp*GLM_chyp*bxstar;
  flux[eqBBX]  = psistar;
  flux[eqERG] += bxstar*psistar;

  //
  // That should do it for the flux
  //

#ifdef FUNCTION_ID
  cout <<"flux_solver_mhd_mixedGLM_adi::inviscid_flux ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


// ##################################################################
// ##################################################################



