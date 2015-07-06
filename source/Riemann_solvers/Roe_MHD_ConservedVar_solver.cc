///
/// \file Roe_MHD_ConservedVar_solver.cc
/// \author Jonathan Mackey
///
/// Linearised Riemann solver for the ideal MHD equations.
/// This is a Roe solver in conserved variables.  It has a one-sided
/// solver and a symmetric solver which maintains symmetry
/// in a problem much more successfully.
/// The symmetric solver is recommended over the one-sided solver.
/// Both will return a flux in the first 8 state variables (i.e.
/// the physical vars) and also a state vector for the interface which
/// can be used for adding viscosity.
/// 
/// References:\n
/// - Cargo & Gallice (1997) JCP, 136, 446
/// - Stone et al. (2009), ApJS, 178, 137
///
/// History:
/// - 2010.12.27 JM: Moved from flux_mhd_adiabatic.h
/// - 2015.01.14 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include "Riemann_solvers/Roe_MHD_ConservedVar_solver.h"
#include "Riemann_solvers/riemannMHD.h"

using namespace std;

// **********************************************************************************
//  ROE SOLVER FOR ADIABATIC MHD, FROM CARGO & GALLICE, (1997) JCP, 136, 446
// **********************************************************************************


// ##################################################################
// ##################################################################



Riemann_Roe_MHD_CV::Riemann_Roe_MHD_CV(const int nv,
     ///< Length of State Vectors, nvar
				       const double g
     ///< Gamma for state vector.
				       )
  : eqns_base(nv), eqns_mhd_ideal(nv)
{
#ifdef FUNCTION_ID
  cout <<"::Riemann_Roe_MHD_CV ...starting.\n";
#endif //FUNCTION_ID
  //
  // eq_gamma, eq_nvar are defined in eqns_base class
  //
  eq_gamma = g;

  //
  // Allocate memory for evalues, wave-strengths, evectors,
  // mean-state, difference-state, left and right conserved var states.
  //
  Roe_evalues = mem.myalloc(Roe_evalues,7);
  Roe_strengths = mem.myalloc(Roe_strengths,7);
  Roe_right_evecs = mem.myalloc(Roe_right_evecs,7);
  for (int v=0;v<7;v++) 
    Roe_right_evecs[v] = mem.myalloc(Roe_right_evecs[v],7);
  Roe_meanp = mem.myalloc(Roe_meanp,eq_nvar);
  Roe_udiff = mem.myalloc(Roe_udiff,eq_nvar);
  Roe_pdiff = mem.myalloc(Roe_pdiff,eq_nvar);
  Roe_UL = mem.myalloc(Roe_UL,eq_nvar);
  Roe_UR = mem.myalloc(Roe_UR,eq_nvar);


  //   for (int v=0;v<7;v++)
  //     for (int w=0;w<7;w++) {
  //       Roe_right_evecs[v][w] = v+w+0.5;
  //       cout <<"v+w="<<Roe_right_evecs[v][w]<<"\n";
  //     }
  
  //
  // Set the enthalpy index to replace pressure.
  //
  eqHH = eqPG;

#ifdef FUNCTION_ID
  cout <<"::Riemann_Roe_MHD_CV ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################



Riemann_Roe_MHD_CV::~Riemann_Roe_MHD_CV()
{
#ifdef FUNCTION_ID
  cout <<"::~Riemann_Roe_MHD_CV ...starting.\n";
#endif //FUNCTION_ID

  //
  // Free Memory:
  //
  Roe_evalues    = mem.myfree(Roe_evalues);
  Roe_strengths  = mem.myfree(Roe_strengths);
  for (int v=0;v<7;v++) 
    Roe_right_evecs[v] = mem.myfree(Roe_right_evecs[v]);
  Roe_right_evecs = mem.myfree(Roe_right_evecs);
  Roe_meanp = mem.myfree(Roe_meanp);
  Roe_UL = mem.myfree(Roe_UL);
  Roe_UR = mem.myfree(Roe_UR);
  Roe_udiff = mem.myfree(Roe_udiff);
  Roe_pdiff = mem.myfree(Roe_pdiff);

#ifdef FUNCTION_ID
  cout <<"::~Riemann_Roe_MHD_CV ...returning.\n";
#endif //FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################



///
/// Roe's approximate flux solver (returning the flux calculated from
/// a one-sided calculation (just jump across waves with v<0).
///
int Riemann_Roe_MHD_CV::MHD_Roe_CV_flux_solver_onesided(const double *left,
							const double *right,
							const double g,
#ifdef HCORR
							const double hc_etamax,
#endif // HCORR
							double *out_pstar,
							double *out_flux
							)
{

  //
  // First test if the left and right states are the same, and if they are
  // then return the flux of the left state:
  //
  double diff=0.0;
  for(int i=0;i<eq_nvar; i++)
    diff += fabs(right[i]-left[i])/(fabs(eq_refvec[i])+TINYVALUE);

  if (diff <1.e-6) {
#ifdef RoeMHD_TESTING
    cout <<"same states...\n";
#endif
    for (int v=0;v<eq_nvar;v++)
      out_pstar[v] = 0.5*(left[v]+right[v]);
    PtoFlux(left,out_flux,eq_gamma);
    return 0;
  }


  int err=0; 
#ifdef RoeMHD_TESTING
  cout <<"\t************\n";
  rep.printVec(" left",left,8);
  rep.printVec("right",right,8);
#endif

  set_UL_and_UR(left,right);

  err += Roe_get_average_state(left,right);

  err += Roe_get_difference_states(left,right);

  err += Roe_get_wavespeeds();

  err += Roe_get_eigenvalues(
#ifdef HCORR
			     hc_etamax
#endif // HCORR
			     );
 
  err += Roe_get_wavestrengths();
 
  err += Roe_get_right_evectors();

  //
  // The flux function calculates both out_pstar[] and out_flux[], the
  // resolved state and its flux.
  //
  err += Roe_get_flux_onesided(left,right,
#ifdef MHD_ROE_USE_USTAR
			       out_pstar,
#endif // MHD_ROE_USE_USTAR
			       out_flux);

#ifndef MHD_ROE_USE_USTAR
  //
  // We need to put something into pstar for the AV calculation, so 
  // if it's not the starred state we use the Roe-average state.
  //
  set_pstar_from_meanp(out_pstar);
#endif // not MHD_ROE_USE_USTAR

  
  //  rep.printVec(" left",left,8);
  //  rep.printVec("right",right,8);
  //rep.printVec("meanp",Roe_meanp,8);
  //rep.printVec("pdiff",Roe_pdiff,8);
  //rep.printVec("udiff",Roe_udiff,8);
  //  for (int i=0; i<7; ++i) {
  //   cout << "rightevec["<<i<<"] = [ ";
  //    for (int j=0; j<7; j++) {
  //      cout.width(9);
  //      cout << Roe_right_evecs[i][j] <<", ";
  //    }
  //    cout << "]" << "\n";
  //   }
  //cout <<"\t*************************************\n";
  //  rep.error("bugging out deliberately!!!",24);
  
  //if (err)
  //  cout <<"Riemann_Roe_MHD_CV::Roe_Conserved_flux_solver() ";
  //  cout <<"Picked up an error calculating fluxes...\n";

  return err;
}


// ##################################################################
// ##################################################################



///
/// Roe's approximate flux solver (returning the flux calculated from
/// a symmetric calculation (all waves contribute).
///
int Riemann_Roe_MHD_CV::MHD_Roe_CV_flux_solver_symmetric(const double *left,
							 const double *right,
							 const double g,
#ifdef HCORR
							 const double hc_etamax,
#endif // HCORR
							 double *out_pstar,
							 double *out_flux
							 )
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::MHD_Roe_CV_flux_solver_symmetric ...starting.\n";
#endif //FUNCTION_ID

  int err=0; 
#ifdef RoeMHD_TESTING
  cout <<"\t************\n";
  rep.printVec(" left",left,8);
  rep.printVec("right",right,8);
#endif

  eq_gamma=g;

  set_UL_and_UR(left,right);
  err += Roe_get_average_state(left,right);
  err += Roe_get_difference_states(left,right);
  err += Roe_get_wavespeeds();
  err += Roe_get_eigenvalues(
#ifdef HCORR
			     hc_etamax
#endif // HCORR
			     );
  err += Roe_get_wavestrengths();
  err += Roe_get_right_evectors();

  //
  // This calculates the flux by summing over all waves, and we put
  // meanP[] as the resolved state in a function call below.
  //
  err += calculate_symmetric_flux(left, right, out_flux);
  
  //
  // We need something for out_pstar, so we use the mean state.
  //
  set_pstar_from_meanp(out_pstar);

#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::MHD_Roe_CV_flux_solver_symmetric ...returning.\n";
#endif //FUNCTION_ID
  return err;
}



// ##################################################################
// ##################################################################




double Riemann_Roe_MHD_CV::Enthalpy(const double *p, ///< primitive State Vector.
					   const double g ///< gas EOS gamma.
					   )
{
  //
  // This returns the enthalpy per unit mass
  // H = (e+p_g+p_b)/rho = (0.5rhoV^2 + p/(g-1)+B^2/2 + p+B^2/2)/rho
  //
  return ( (p[eqRO]*(p[eqVX]*p[eqVX]+p[eqVY]*p[eqVY]+p[eqVZ]*p[eqVZ])/2.0
	    +(g*p[eqPG]/(g-1.0))
	    +(p[eqBX]*p[eqBX] +p[eqBY]*p[eqBY] +p[eqBZ]*p[eqBZ]) 
	    )/p[eqRO] );
}


// ##################################################################
// ##################################################################



///
/// Set UL[] and UR[] from PL[] and PR[].
///
void Riemann_Roe_MHD_CV::set_UL_and_UR(const double *left, ///< left primitive vec.
				       const double *right  ///< right primitive vec.
				       )
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::set_UL_and_UR ...starting.\n";
#endif //FUNCTION_ID

  PtoU(left ,Roe_UL,eq_gamma);
  PtoU(right,Roe_UR,eq_gamma);

#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::set_UL_and_UR ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################



///
/// Set Pstar[] from Roe_meanp[] (need to replace enthalpy with
/// pressure).
///
void Riemann_Roe_MHD_CV::set_pstar_from_meanp(double *out_pstar)
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::set_pstar_from_meanp ...starting.\n";
#endif //FUNCTION_ID

  for (int v=0;v<eq_nvar;v++)
    out_pstar[v] = Roe_meanp[v];
  //
  // Remember the pressure is actually enthalpy in the mean state, so
  // reset it to pressure, where we do LL+RR averaging, as oppose to
  // RL+LR.
  //
  //double 
  //  rl = sqrt(left[eqRO]),
  //  rr = sqrt(right[eqRO]);
  //out_pstar[eqPG] = (rl*left[eqPG]+rr*right[eqPG])/(rl+rr);
  //
  // Instead just convert from enthalpy to pressure:
  // p_g = [(g-1)/g]*[H*rho -rho*v^2/2 -B^2]
  //
  out_pstar[eqPG] = ((eq_gamma-1.0)/eq_gamma)*
    (out_pstar[eqRO]*(out_pstar[eqHH]-0.5*Roe_V*Roe_V)-Roe_B*Roe_B);
  //if (out_pstar[eqPG] <0.0) {
  // cout <<"pressure: pl="<<left[eqPG]<<" pr="<<right[eqPG]<<" pm="<<out_pstar[eqPG];
  // cout <<" roe-avg="<< (sqrt(left[eqRO])*left[eqPG]
  // +sqrt(right[eqRO])*right[eqPG])/(sqrt(left[eqRO])+sqrt(right[eqRO]))<<"\n";
  //}

#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::set_pstar_from_meanp ...returning.\n";
#endif //FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################




///
/// Get the Roe averages for the conserved variables:
/// From Stone et al. (2009), ApJS, 178, 137, eq.65.
///
int Riemann_Roe_MHD_CV::Roe_get_average_state(const double *left,
						     const double *right
						     )
{
  // 
  // This gets the average of the primitive variables
  // rho, vx, vy, vz, Bx, By, Bz, and H (not P!!!)
  double 
    rl = sqrt(left[eqRO]),
    rr = sqrt(right[eqRO]),
    lH = Enthalpy(left ,eq_gamma),
    rH = Enthalpy(right,eq_gamma);

  Roe_denom = 1.0/(rl+rr),

  Roe_meanp[eqRO] = rl*rr;
  Roe_meanp[eqVX] = (rl*left[eqVX]+rr*right[eqVX])*Roe_denom;
  Roe_meanp[eqVY] = (rl*left[eqVY]+rr*right[eqVY])*Roe_denom;
  Roe_meanp[eqVZ] = (rl*left[eqVZ]+rr*right[eqVZ])*Roe_denom;

  Roe_meanp[eqBY] = (rr*left[eqBY]+rl*right[eqBY])*Roe_denom;
  Roe_meanp[eqBZ] = (rr*left[eqBZ]+rl*right[eqBZ])*Roe_denom;
  //
  // Bx is just a parameter, and GLM method sets both left and right
  // states to be Bx(star)
  //
  Roe_meanp[eqBX] = 0.5*(left[eqBX] + right[eqBX]);
  if (Roe_meanp[eqBX]>=0.0)
    Roe_signBX = 1;
  else
    Roe_signBX = -1;

  Roe_meanp[eqHH] = (rl*lH + rr*rH)*Roe_denom;

  //
  // Set velocity and field magnitudes for mean state:
  //
  Roe_V = sqrt(Roe_meanp[eqVX]*Roe_meanp[eqVX]+
	       Roe_meanp[eqVY]*Roe_meanp[eqVY]+
	       Roe_meanp[eqVZ]*Roe_meanp[eqVZ]);

  Roe_B = sqrt(Roe_meanp[eqBX]*Roe_meanp[eqBX]+
	       Roe_meanp[eqBY]*Roe_meanp[eqBY]+
	       Roe_meanp[eqBZ]*Roe_meanp[eqBZ]);

  //
  // Now the parameters Bt, betay, betaz:
  //
  Roe_Bt = sqrt(Roe_meanp[eqBY]*Roe_meanp[eqBY]+Roe_meanp[eqBZ]*Roe_meanp[eqBZ]);
  if (Roe_Bt >= TINYVALUE) {
    Roe_betay = Roe_meanp[eqBY]/Roe_Bt;
    Roe_betaz = Roe_meanp[eqBZ]/Roe_Bt;
  }
  else {
    Roe_betay = 1.0/sqrt(2.0);
    Roe_betaz = 1.0/sqrt(2.0);
  }
//  cout <<"meanp: v="<<Roe_V<<" B="<<Roe_B<<" Bt="<<Roe_Bt<<" betay="<<Roe_betay<<" betaz="<<Roe_betaz<<"\n";
  return 0;
}


// ##################################################################
// ##################################################################



///
/// Get the Roe averages differences for the primitive and conserved variables:
/// From Stone et al. (2009), ApJS, 178, 137, eq.65.
///
int Riemann_Roe_MHD_CV::Roe_get_difference_states(const double *left,
						  const double *right
						  )
{
  // double Ul[eq_nvar], Ur[eq_nvar];
  // PtoU(left,  Ul, eq_gamma);
  // PtoU(right, Ur, eq_gamma);
  //
  // The conserved vectors UL,UR should already be set.
  //
  for (int v=0; v<eq_nvar; v++) {
    Roe_udiff[v] = Roe_UR[v]-Roe_UL[v];
    Roe_pdiff[v] = right[v]-left[v];
  }
  //  rep.printVec("pdiff",Roe_pdiff,8);
  //  rep.printVec("udiff",Roe_udiff,8);

  //
  // Enforce no jump in BX, just in case I try to calculate it anywhere:
  //
  Roe_udiff[eqBBX] = Roe_pdiff[eqBX] = 0.0;

  //
  // Set GC97's parameter X:
  //
  Roe_CGparamX = (Roe_pdiff[eqBY]*Roe_pdiff[eqBY] 
		  +Roe_pdiff[eqBZ]*Roe_pdiff[eqBZ])*0.5*Roe_denom*Roe_denom;

  //
  // delta(p) is different though...
  // See CG97 equation 4.15+1
  //
  Roe_pdiff[eqPG] = ( (0.5*Roe_V*Roe_V -Roe_CGparamX)*Roe_pdiff[eqRO]
		      -(Roe_meanp[eqVX]*Roe_udiff[eqMMX]+
			Roe_meanp[eqVY]*Roe_udiff[eqMMY]+
			Roe_meanp[eqVZ]*Roe_udiff[eqMMZ])
		      +Roe_udiff[eqERG]
		      -(Roe_meanp[eqBY]*Roe_pdiff[eqBY]+
			Roe_meanp[eqBZ]*Roe_pdiff[eqBZ])
		      )*(eq_gamma-1.0);
  
  //  rep.printVec("pdiff",Roe_pdiff,8);
  //  rep.printVec("udiff",Roe_udiff,8);
  return 0;
}


// ##################################################################
// ##################################################################



///
/// Get the Roe-averaged wavespeeds
///
int Riemann_Roe_MHD_CV::Roe_get_wavespeeds()
{
  //
  // Mean wave speeds from CG97 eq.4.17
  // 
  double b2 = Roe_B*Roe_B/Roe_meanp[eqRO];

  Roe_a = sqrt((2.0-eq_gamma)*Roe_CGparamX +
	       (eq_gamma-1.0)*(Roe_meanp[eqHH] -0.5*Roe_V*Roe_V -b2));
  double astar2 = Roe_a*Roe_a +b2;

  Roe_ca = sqrt(Roe_meanp[eqBX]*Roe_meanp[eqBX]/Roe_meanp[eqRO]);
//  cout <<"a="<<Roe_a<<"  c_a="<<Roe_ca;

  //
  // Need to do some checks for fast and slow speeds to make sure
  // we're not getting square roots of negative numbers.
  //
  Roe_cs = astar2*astar2 -4.0*Roe_a*Roe_a*Roe_ca*Roe_ca;
  if (Roe_cs <= 0.0) {
    //    cout <<"\t\tfirst sqrt <0 in fast/slow speeds! value="<<Roe_cs<<"\n";
    Roe_cs = 0.0;
  }
  else 
    Roe_cs = sqrt(Roe_cs);
  
  Roe_cf = sqrt(0.5*(astar2 + Roe_cs));

  Roe_cs = astar2 - Roe_cs;
  if (Roe_cs <= 0.0) {
    //cout <<"\t\tsecond sqrt <0 for slow speed! value="<<Roe_cs<<"\n";
    Roe_cs = 0.0;
  }
  else 
    Roe_cs = sqrt(0.5*Roe_cs);

  //cout <<"  c_s="<<Roe_cs<<"  c_f="<<Roe_cf;//<<"\n";

  //
  // So we have a, c_a, c_s, c_f from the mean state vector.
  // Check that they are in the right order:
  //
  if (Roe_ca > Roe_cf) {
    //cout <<"Alfven speed faster than fast speed! (1-ca/cf)="<<1.0-Roe_ca/Roe_cf;
    //cout <<", Resetting so they are equal.\n";
    Roe_ca = Roe_cf;
  }
  if (Roe_cs > Roe_ca) {
    //cout <<"slow speed faster than Alfven speed! (1-cs/ca)="<<1.0-Roe_cs/Roe_ca;
    //cout <<", Resetting so they are equal.\n";
    Roe_cs = Roe_ca;
  }

  //
  // Now get the Roe-Balsara normalisation parameters, taking care that they
  // are finite and sensible...
  //
  double cf2diff;
  if ((cf2diff= Roe_cf*Roe_cf-Roe_cs*Roe_cs) > MACHINEACCURACY) { 
    if ((Roe_alphaf = Roe_a*Roe_a-Roe_cs*Roe_cs) <0.0) Roe_alphaf=0.;
    if ((Roe_alphas = Roe_cf*Roe_cf-Roe_a*Roe_a) <0.0) Roe_alphas=0.;
    if ((Roe_alphaf = sqrt(Roe_alphaf/cf2diff)) >1.0) {
      // cout<<"Roe_alpha_f = "<<Roe_alphaf<<" !!!\n";
      Roe_alphaf =1.0;
    }
    if ((Roe_alphas = sqrt(Roe_alphas/cf2diff)) >1.0) {
      // cout<<"Roe_alpha_s>1!!!\n";
      Roe_alphas =1.0;
    }
  }
  else {
    Roe_alphaf = Roe_alphas = 1.0/sqrt(2.0);
  }
//  cout <<"\talpha_s="<<Roe_alphas<<" alpha_f="<<Roe_alphaf<<"\n";
  return 0;
}



// ##################################################################
// ##################################################################



///
/// Get the Roe-averaged eigenvalues
///
int Riemann_Roe_MHD_CV::Roe_get_eigenvalues(
#ifdef HCORR
					    const double Hcorr_etamax
#endif // HCORR
					    )
{
  //
  // We can just write these down now since all the hard work is done:
  //
  Roe_evalues[FN] = Roe_meanp[eqVX] - Roe_cf;
  Roe_evalues[AN] = Roe_meanp[eqVX] - Roe_ca;
  Roe_evalues[SN] = Roe_meanp[eqVX] - Roe_cs;
  Roe_evalues[CT] = Roe_meanp[eqVX];
  Roe_evalues[SP] = Roe_meanp[eqVX] + Roe_cs;
  Roe_evalues[AP] = Roe_meanp[eqVX] + Roe_ca;
  Roe_evalues[FP] = Roe_meanp[eqVX] + Roe_cf;
#ifdef RoeMHD_TESTING
  rep.printVec("e-values",Roe_evalues,7);
#endif

#ifdef HCORR
#ifdef TESTING
  //
  // Paranoid test!  Make sure eta=0 if not using H-correction.
  //
  if (SimPM.artviscosity!=3 && !pconst.equalD(Hcorr_etamax,0.0))
    rep.error("H-correction is non-zero but we're not using it!",Hcorr_etamax);
#endif // TESTING

  //
  // Modify the eigenvalues by the H-correction eta value.  Note that
  // HC_etamax is set to zero in the simulation initialisation, and it
  // is only changed if AVtype==3, so this code has no effect if we
  // are not using the H-correction.  Unless of course the eigenvalue
  // is _very_ close to zero and the eval changes within the machine
  // precision.  This shouldn't ever have a significant effect on
  // anything!
  //
  for (int v=0; v<7; v++) {
    //
    // if eval[v] <0  then set it to min(eval[v],-eta)
    // else                set it to max(eval[v], eta)
    //
    if (Roe_evalues[v]<0.0) {
      Roe_evalues[v] = min(Roe_evalues[v],-Hcorr_etamax);
    }
    else {
      Roe_evalues[v] = max(Roe_evalues[v], Hcorr_etamax);
    }
  }
#ifdef RoeMHD_TESTING
  rep.printVec("e-values after H-correction",Roe_evalues,7);
#endif
#endif // HCORR

  return 0;
}


// ##################################################################
// ##################################################################



int Riemann_Roe_MHD_CV::Roe_get_wavestrengths()
{
  ///
  /// Get the Roe-averaged wave strengths,
  /// from Cargo & Gallice (1997) JCP, 136, 446, eq.4.20
  ///
 //cout <<"alphaf,alphas="<<Roe_alphaf<<", "<<Roe_alphas<<"\n";
 //cout <<"X="<<Roe_CGparamX<<" by,bz="<<Roe_betay<<", "<<Roe_betaz<<"\n";
 //cout <<"a="<<Roe_a<<" cs="<<Roe_cs<<" ca="<<Roe_ca<<" cf="<<Roe_cf<<"\n";
 //rep.printVec("meanp",Roe_meanp,8);
 //rep.printVec("pdiff",Roe_pdiff,8);

  Roe_strengths[FN] = 
    0.5*( Roe_alphaf*(Roe_CGparamX*Roe_pdiff[eqRO] +Roe_pdiff[eqPG])
	  +Roe_meanp[eqRO]*Roe_alphas*Roe_cs*Roe_signBX*(Roe_betay*Roe_pdiff[eqVY] +
							Roe_betaz*Roe_pdiff[eqVZ])
	  -Roe_meanp[eqRO]*Roe_alphaf*Roe_cf*Roe_pdiff[eqVX]
	  +sqrt(Roe_meanp[eqRO])*Roe_alphas*Roe_a*(Roe_betay*Roe_pdiff[eqBY] +
						  Roe_betaz*Roe_pdiff[eqBZ])
	  );

  Roe_strengths[FP] = 
    0.5*( Roe_alphaf*(Roe_CGparamX*Roe_pdiff[eqRO] +Roe_pdiff[eqPG])
	  -Roe_meanp[eqRO]*Roe_alphas*Roe_cs*Roe_signBX*(Roe_betay*Roe_pdiff[eqVY] +
							Roe_betaz*Roe_pdiff[eqVZ])
	  +Roe_meanp[eqRO]*Roe_alphaf*Roe_cf*Roe_pdiff[eqVX]
	  +sqrt(Roe_meanp[eqRO])*Roe_alphas*Roe_a*(Roe_betay*Roe_pdiff[eqBY] +
						  Roe_betaz*Roe_pdiff[eqBZ])
	  );

  Roe_strengths[SN] = 
    0.5*( Roe_alphas*(Roe_CGparamX*Roe_pdiff[eqRO] +Roe_pdiff[eqPG])
	  -Roe_meanp[eqRO]*Roe_alphaf*Roe_cf*Roe_signBX*(Roe_betay*Roe_pdiff[eqVY] +
							Roe_betaz*Roe_pdiff[eqVZ])
	  -Roe_meanp[eqRO]*Roe_alphas*Roe_cs*Roe_pdiff[eqVX]
	  -sqrt(Roe_meanp[eqRO])*Roe_alphaf*Roe_a*(Roe_betay*Roe_pdiff[eqBY] +
						  Roe_betaz*Roe_pdiff[eqBZ])
	  );
  
  Roe_strengths[SP] = 
    0.5*( Roe_alphas*(Roe_CGparamX*Roe_pdiff[eqRO] +Roe_pdiff[eqPG])
	  +Roe_meanp[eqRO]*Roe_alphaf*Roe_cf*Roe_signBX*(Roe_betay*Roe_pdiff[eqVY] +
							Roe_betaz*Roe_pdiff[eqVZ])
	  +Roe_meanp[eqRO]*Roe_alphas*Roe_cs*Roe_pdiff[eqVX]
	  -sqrt(Roe_meanp[eqRO])*Roe_alphaf*Roe_a*(Roe_betay*Roe_pdiff[eqBY] +
						  Roe_betaz*Roe_pdiff[eqBZ])
	  );

  Roe_strengths[AN] =
    0.5*( +Roe_betay*Roe_pdiff[eqVZ] -Roe_betaz*Roe_pdiff[eqVY]
	  +Roe_signBX*(Roe_betay*Roe_pdiff[eqBZ] -Roe_betaz*Roe_pdiff[eqBY])/sqrt(Roe_meanp[eqRO])
	  );

  Roe_strengths[AP] =
    0.5*( -Roe_betay*Roe_pdiff[eqVZ] +Roe_betaz*Roe_pdiff[eqVY]
	  +Roe_signBX*(Roe_betay*Roe_pdiff[eqBZ] -Roe_betaz*Roe_pdiff[eqBY])/sqrt(Roe_meanp[eqRO])
	  );

  Roe_strengths[CT] = (Roe_a*Roe_a - Roe_CGparamX)*Roe_pdiff[eqRO] - Roe_pdiff[eqPG];

#ifdef RoeMHD_TESTING
  rep.printVec("strengths:",Roe_strengths,7);
  //cout <<"str[FN] = "<<Roe_strengths[FN]<<"\n";
  //cout <<"str[AN] = "<<Roe_strengths[AN]<<"\n";
  //cout <<"str[SN] = "<<Roe_strengths[SN]<<"\n";
  //cout <<"str[CT] = "<<Roe_strengths[CT]<<"\n";
  //cout <<"str[SP] = "<<Roe_strengths[SP]<<"\n";
  //cout <<"str[AP] = "<<Roe_strengths[AP]<<"\n";
  //cout <<"str[FP] = "<<Roe_strengths[FP]<<"\n";
#endif
  return 0;
}


// ##################################################################
// ##################################################################




///
/// Calculate the Right Eigenvectors of the average state,
/// from Cargo & Gallice (1997) JCP, 136, 446, eq.4.18,4.19
///
int Riemann_Roe_MHD_CV::Roe_get_right_evectors()
{
  //
  // Elements in the eigenvectors go as:
  // [rho,vx,vy,vz,By,Bz,Energy] = [0,1,2,3,4,5,6]
  //
  double rootrho = sqrt(Roe_meanp[eqRO]);
  //
  // contact
  //
  Roe_right_evecs[CT][0] = 1;
  Roe_right_evecs[CT][1] = Roe_meanp[eqVX];
  Roe_right_evecs[CT][2] = Roe_meanp[eqVY];
  Roe_right_evecs[CT][3] = Roe_meanp[eqVZ];
  Roe_right_evecs[CT][4] = 0.0;
  Roe_right_evecs[CT][5] = 0.0;
  Roe_right_evecs[CT][6] = 0.5*Roe_V*Roe_V +Roe_CGparamX*(eq_gamma-2)/(eq_gamma-1);
  for (int v=0;v<7;v++)
    Roe_right_evecs[CT][v] /= Roe_a*Roe_a;

  //
  // alfven negative
  //
  Roe_right_evecs[AN][0] = 0.0;
  Roe_right_evecs[AN][1] = 0.0;
  Roe_right_evecs[AN][2] = -Roe_meanp[eqRO]*Roe_betaz;
  Roe_right_evecs[AN][3] = +Roe_meanp[eqRO]*Roe_betay;
  Roe_right_evecs[AN][4] = -Roe_signBX*rootrho*Roe_betaz;
  Roe_right_evecs[AN][5] = +Roe_signBX*rootrho*Roe_betay;
  Roe_right_evecs[AN][6] = -Roe_meanp[eqRO]*(Roe_meanp[eqVY]*Roe_betaz -Roe_meanp[eqVZ]*Roe_betay);

  //
  // alfven positive
  //
  Roe_right_evecs[AP][0] = 0.0;
  Roe_right_evecs[AP][1] = 0.0;
  Roe_right_evecs[AP][2] = -Roe_right_evecs[AN][2]; //+Roe_meanp[eqRO]*Roe_betaz;
  Roe_right_evecs[AP][3] = -Roe_right_evecs[AN][3]; //-Roe_meanp[eqRO]*Roe_betay;
  Roe_right_evecs[AP][4] =  Roe_right_evecs[AN][4]; //-Roe_signBX*rootrho*Roe_betaz;
  Roe_right_evecs[AP][5] =  Roe_right_evecs[AN][5]; //+Roe_signBX*rootrho*Roe_betay;
  Roe_right_evecs[AP][6] = -Roe_right_evecs[AN][6]; //+Roe_meanp[eqRO]*(Roe_meanp[eqVY]*Roe_betaz -Roe_meanp[eqVZ]*Roe_betay);


  double dens_alphas = Roe_meanp[eqRO]*Roe_alphas;
  double dens_alphaf = Roe_meanp[eqRO]*Roe_alphaf;
  //
  // slow negative
  //
  Roe_right_evecs[SN][0] = dens_alphas;
  Roe_right_evecs[SN][1] = dens_alphas*(Roe_meanp[eqVX]-Roe_cs);
  Roe_right_evecs[SN][2] = dens_alphas*Roe_meanp[eqVY] -dens_alphaf*Roe_cf*Roe_betay*Roe_signBX;
  Roe_right_evecs[SN][3] = dens_alphas*Roe_meanp[eqVZ] -dens_alphaf*Roe_cf*Roe_betaz*Roe_signBX;
  Roe_right_evecs[SN][4] = -rootrho*Roe_alphaf*Roe_a*Roe_betay;
  Roe_right_evecs[SN][5] = -rootrho*Roe_alphaf*Roe_a*Roe_betaz;
  Roe_right_evecs[SN][6] = dens_alphas*(Roe_meanp[eqHH]-Roe_B*Roe_B/Roe_meanp[eqRO]-Roe_meanp[eqVX]*Roe_cs)
    -dens_alphaf*Roe_cf*Roe_signBX*(Roe_meanp[eqVY]*Roe_betay+Roe_meanp[eqVZ]*Roe_betaz)
    -rootrho*Roe_alphaf*Roe_a*Roe_Bt;

  //
  // slow positive
  //
  Roe_right_evecs[SP][0] = dens_alphas;            //Roe_meanp[eqRO]*Roe_alphas;
  Roe_right_evecs[SP][1] = dens_alphas*(Roe_meanp[eqVX]+Roe_cs);
  Roe_right_evecs[SP][2] = dens_alphas*Roe_meanp[eqVY] +dens_alphaf*Roe_cf*Roe_betay*Roe_signBX;
  Roe_right_evecs[SP][3] = dens_alphas*Roe_meanp[eqVZ] +dens_alphaf*Roe_cf*Roe_betaz*Roe_signBX;
  Roe_right_evecs[SP][4] = Roe_right_evecs[SN][4]; //-rootrho*Roe_alphaf*Roe_a*Roe_betay;
  Roe_right_evecs[SP][5] = Roe_right_evecs[SN][5]; //-rootrho*Roe_alphaf*Roe_a*Roe_betaz;
  Roe_right_evecs[SP][6] = dens_alphas*(Roe_meanp[eqHH]-Roe_B*Roe_B/Roe_meanp[eqRO]+Roe_meanp[eqVX]*Roe_cs)
    +dens_alphaf*Roe_cf*Roe_signBX*(Roe_meanp[eqVY]*Roe_betay+Roe_meanp[eqVZ]*Roe_betaz)
    -rootrho*Roe_alphaf*Roe_a*Roe_Bt;

  //
  // fast negative
  //
  Roe_right_evecs[FN][0] = dens_alphaf;
  Roe_right_evecs[FN][1] = dens_alphaf*(Roe_meanp[eqVX]-Roe_cf);
  Roe_right_evecs[FN][2] = dens_alphaf*Roe_meanp[eqVY] +dens_alphas*Roe_cs*Roe_betay*Roe_signBX;
  Roe_right_evecs[FN][3] = dens_alphaf*Roe_meanp[eqVZ] +dens_alphas*Roe_cs*Roe_betaz*Roe_signBX;
  Roe_right_evecs[FN][4] = rootrho*Roe_alphas*Roe_a*Roe_betay;
  Roe_right_evecs[FN][5] = rootrho*Roe_alphas*Roe_a*Roe_betaz;
  Roe_right_evecs[FN][6] = dens_alphaf*(Roe_meanp[eqHH]-Roe_B*Roe_B/Roe_meanp[eqRO]-Roe_meanp[eqVX]*Roe_cf)
    +dens_alphas*Roe_cs*Roe_signBX*(Roe_meanp[eqVY]*Roe_betay+Roe_meanp[eqVZ]*Roe_betaz)
    +rootrho*Roe_alphas*Roe_a*Roe_Bt;

  //
  // fast positive
  //
  Roe_right_evecs[FP][0] = dens_alphaf;            //Roe_meanp[eqRO]*Roe_alphaf;
  Roe_right_evecs[FP][1] = dens_alphaf*(Roe_meanp[eqVX]+Roe_cf);
  Roe_right_evecs[FP][2] = dens_alphaf*Roe_meanp[eqVY] -dens_alphas*Roe_cs*Roe_betay*Roe_signBX;
  Roe_right_evecs[FP][3] = dens_alphaf*Roe_meanp[eqVZ] -dens_alphas*Roe_cs*Roe_betaz*Roe_signBX;
  Roe_right_evecs[FP][4] = Roe_right_evecs[FN][4]; //rootrho*Roe_alphas*Roe_a*Roe_betay;
  Roe_right_evecs[FP][5] = Roe_right_evecs[FN][5]; //rootrho*Roe_alphas*Roe_a*Roe_betaz;
  Roe_right_evecs[FP][6] = dens_alphaf*(Roe_meanp[eqHH]-Roe_B*Roe_B/Roe_meanp[eqRO]+Roe_meanp[eqVX]*Roe_cf)
    -dens_alphas*Roe_cs*Roe_signBX*(Roe_meanp[eqVY]*Roe_betay+Roe_meanp[eqVZ]*Roe_betaz)
    +rootrho*Roe_alphas*Roe_a*Roe_Bt;

  //
  // Normalise Fast and Slow Wave Vectors:
  //
  double norm = Roe_meanp[eqRO]*Roe_a*Roe_a;
  for (int v=0;v<7;v++)
    Roe_right_evecs[SN][v] /= norm;
  for (int v=0;v<7;v++)
    Roe_right_evecs[SP][v] /= norm;
  for (int v=0;v<7;v++)
    Roe_right_evecs[FN][v] /= norm;
  for (int v=0;v<7;v++)
    Roe_right_evecs[FP][v] /= norm;

#ifdef RoeMHD_TESTING
//  for (int vec1=0; vec1<7; vec1++) {
//    for (int vec2=0; vec2<7; vec2++) {
//      double dp=0.0;
//      for (int v=0;v<7;v++)
//	dp += Roe_right_evecs[vec1][v]*Roe_right_evecs[vec2][v];
//      cout <<"vec["<<vec1<<"].vec["<<vec2<<"] = "<<dp<<"\n";
//    }
//  }
#endif
  
  return 0;
}



// ##################################################################
// ##################################################################



///
/// Using the evalues,wave-strengths,evectors, calculate the
/// Roe-average Flux from the left state across to zero.
///
int Riemann_Roe_MHD_CV::Roe_get_flux_onesided(const double *left,
					      const double *right,
#ifdef MHD_ROE_USE_USTAR
					      double *out_pstar,
#endif // MHD_ROE_USE_USTAR
					      double *out_flux)
{
  //
  // We also need to get the starred state by jumping across waves, so 
  // we use this variable to track that.
  //
#ifdef MHD_ROE_USE_USTAR
  double Ustar[eq_nvar];
#ifdef RoeMHD_TESTING
  double UstarR[eq_nvar];  
#endif
#endif // MHD_ROE_USE_USTAR
  //
  // If the negative e-value is positive, then the flux is just the left state flux.
  // Conversely if the positive e-value is negative, the solution is the right state.
  //
  if      (Roe_evalues[FN] >= 0.0) {  // if(u_av-c_av>0) pstar = left state
    PUtoFlux(left,Roe_UL,out_flux);
#ifdef MHD_ROE_USE_USTAR
    for (int v=0;v<eq_nvar;v++)
      out_pstar[v] = left[v];
#endif // MHD_ROE_USE_USTAR
    return 0;
  }

  else if (Roe_evalues[FP] <=0.) {  //else if(u_av+c_av<=0) pstar = right state
    PUtoFlux(right,Roe_UR,out_flux);
#ifdef MHD_ROE_USE_USTAR
    for (int v=0;v<eq_nvar;v++)
      out_pstar[v] = right[v];
#endif // MHD_ROE_USE_USTAR
    return 0;
  }

  else {  // else we are in the starred region, so we need to get the flux the long way
    //
    // First get the left state flux:
    // 
    //rep.printVec("left  flux:",out_flux,8);
    PUtoFlux(left, Roe_UL, out_flux);
#ifdef MHD_ROE_USE_USTAR
    //PtoU(left,Ustar,eq_gamma);
    for (int v=0;v<eq_nvar;v++)
      Ustar[v] = Roe_UL[v];
#endif // MHD_ROE_USE_USTAR
#ifdef RoeMHD_TESTING
    //rep.printVec("left  flux:",out_flux,8);
#endif
    
    //
    // Next add the contribution from all the waves with negative 
    // wavespeeds:
    //    
    int i=0;
    while ((i<7) && (Roe_evalues[i]<=0.0)) {
      //
      // Flux:
      //
      out_flux[eqRHO] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][0];
      out_flux[eqMMX] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][1];
      out_flux[eqMMY] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][2];
      out_flux[eqMMZ] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][3];
      out_flux[eqBBY] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][4];
      out_flux[eqBBZ] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][5];
      out_flux[eqERG] += Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][6];
#ifdef MHD_ROE_USE_USTAR
      //
      // Conserved vector state:
      //
      Ustar[eqRHO] += Roe_strengths[i]*Roe_right_evecs[i][0];
      Ustar[eqMMX] += Roe_strengths[i]*Roe_right_evecs[i][1];
      Ustar[eqMMY] += Roe_strengths[i]*Roe_right_evecs[i][2];
      Ustar[eqMMZ] += Roe_strengths[i]*Roe_right_evecs[i][3];
      Ustar[eqBBY] += Roe_strengths[i]*Roe_right_evecs[i][4];
      Ustar[eqBBZ] += Roe_strengths[i]*Roe_right_evecs[i][5];
      Ustar[eqERG] += Roe_strengths[i]*Roe_right_evecs[i][6];
#endif // MHD_ROE_USE_USTAR
      i++;
    }
    out_flux[eqBBX] = 0.0; // by definition!

#ifdef RoeMHD_TESTING
    rep.printVec("left  flux:",out_flux,8);
#endif
    //cout <<"0-7 ="<<eqRHO<<","<<eqERG<<","<<eqMMX<<","<<eqMMY<<","<<eqMMZ<<","<<eqBBX<<","<<eqBBY<<","<<eqBBZ<<"\n";
    //rep.printVec("  evalues:",Roe_evalues,7);
    //rep.printVec("strengths:",Roe_strengths,7);
    //rep.printVec("left  flux:",out_flux,8);

#ifdef RoeMHD_TESTING
    //
    // Now we should be at the flux across the boundary.  Check this
    // by going back from the right to the left:
    //
    double ftemp[eq_nvar];
    PUtoFlux(right, Roe_UR, ftemp);
#ifdef MHD_ROE_USE_USTAR
    //PtoU(right,UstarR,eq_gamma);
    for (int v=0;v<eq_nvar;v++)
      UstarR[v] = Roe_UR[v];
#endif // MHD_ROE_USE_USTAR
    i=6;
    while ((i>=0) && (Roe_evalues[i]>0.0)) {
      //
      // Flux:
      //
      ftemp[eqRHO] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][0];
      ftemp[eqMMX] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][1];
      ftemp[eqMMY] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][2];
      ftemp[eqMMZ] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][3];
      ftemp[eqBBY] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][4];
      ftemp[eqBBZ] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][5];
      ftemp[eqERG] -= Roe_strengths[i]*Roe_evalues[i]*Roe_right_evecs[i][6];
#ifdef MHD_ROE_USE_USTAR
      //
      // Conserved vector state:
      //
      UstarR[eqRHO] -= Roe_strengths[i]*Roe_right_evecs[i][0];
      UstarR[eqMMX] -= Roe_strengths[i]*Roe_right_evecs[i][1];
      UstarR[eqMMY] -= Roe_strengths[i]*Roe_right_evecs[i][2];
      UstarR[eqMMZ] -= Roe_strengths[i]*Roe_right_evecs[i][3];
      UstarR[eqBBY] -= Roe_strengths[i]*Roe_right_evecs[i][4];
      UstarR[eqBBZ] -= Roe_strengths[i]*Roe_right_evecs[i][5];
      UstarR[eqERG] -= Roe_strengths[i]*Roe_right_evecs[i][6];
#endif // MHD_ROE_USE_USTAR
      i--;
    }
    ftemp[eqBBX] = 0.0; // by definition!
//#ifdef RoeMHD_TESTING
    rep.printVec("right flux:",ftemp,8);
//#endif
    
    double diff=0.0;
    for (int v=0;v<8;v++) {
      //cout <<"flux["<<v<<"]: left="<<out_flux[v]<<" and right="<<ftemp[v]<<"\n";
      diff += (out_flux[v]-ftemp[v])/(fabs(out_flux[v])+fabs(ftemp[v])+TINYVALUE);
    }
    if (diff>1e-3) {
      cout <<"*** FLUX CALCULATION ERROR IN Riemann_Roe_MHD_CV::Roe_get_flux(): diff = "<<diff<<"\n";
      rep.printVec("left  flux:",out_flux,8);
      rep.printVec("right flux:",ftemp,8);
      rep.printVec("e-values",Roe_evalues,7);
      rep.printVec("strengths:",Roe_strengths,7);
      rep.printVec(" left",left,8);
      rep.printVec("right",right,8);
      rep.printVec("meanp",Roe_meanp,8);
      rep.printVec("pdiff",Roe_pdiff,8);
      rep.printVec("udiff",Roe_udiff,8);
      for (int i=0; i<7; ++i) {
	cout << "rightevec["<<i<<"] = [ ";
	for (int j=0; j<7; j++) {
	  cout.width(9);
	  cout << Roe_right_evecs[i][j] <<", ";
	}
	cout << "]" << "\n";
      }
      //cout <<"\t*************************************\n";
      return 1;
    }
#ifdef MHD_ROE_USE_USTAR
    //
    // Now check that Ustar[] and UstarR[] are the same:
    //
    diff=0.0;
    for (int v=0;v<8;v++) {
      diff += (Ustar[v]-UstarR[v])/eq_refvec[v];
    }
    if (diff>1e-3) {
      cout <<"*** FLUX CALCULATION ERROR IN Riemann_Roe_MHD_CV::Roe_get_flux(): diff = "<<diff<<"\n";
      rep.printVec(" left",Ustar, 8);
      rep.printVec("right",UstarR,8);
    }
#endif // MHD_ROE_USE_USTAR
#endif // RoeMHD_TESTING
  }


#ifdef MHD_ROE_USE_USTAR
  //
  // Finally assign Pstar with the values from Ustar[]:
  //
  int err = UtoP(Ustar,out_pstar,eq_gamma);
  //
  // Check for Errors!
  //
  switch (err) {
  case 0:
    break;
  case 1:
    negPGct ++;
    if (negPGct<1000) {
      cout <<"Got negative pressure in Roe starred state! err="<<err<<"\n";
      rep.printVec("U*L",Ustar ,8);
#ifdef RoeMHD_TESTING
      rep.printVec("U*R",UstarR,8);
#endif // RoeMHD_TESTING
      rep.printVec("P* ",out_pstar,8);
    }
    break;
  case 2:
    negROct ++;
    if (negROct<1000) {
      cout <<"Got negative pressure in Roe starred state! err="<<err<<"\n";
      rep.printVec("U*L",Ustar ,8);
#ifdef RoeMHD_TESTING
      rep.printVec("U*R",UstarR,8);
#endif // RoeMHD_TESTING
      rep.printVec("P* ",out_pstar,8);
    }
    //
    // Use Roe-average density instead of starred state negative density
    //
    Ustar[eqRHO] = Roe_meanp[eqRO];
    break;
  default:
    rep.error("unhandled return error code from UtoP()",err);
    break;
  }

  return err;
#else  // not MHD_ROE_USE_USTAR

  return 0;
#endif // MHD_ROE_USE_USTAR

}



// ##################################################################
// ##################################################################



int Riemann_Roe_MHD_CV::calculate_symmetric_flux(const double *left,
						 const double *right,
						 double *out_flux
						 )
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::calculate_symmetric_flux ...starting.\n";
#endif //FUNCTION_ID
  //
  // Get the flux by stepping across waves:
  //
  // Flux = 0.5(F(left)+F(right)
  //            -sum_{waves}[strength_i*|evalue_i|*right_evec_i])
  //
  // We use the left state vector Roe_UL[] as a temp array once it has
  // made its contribution to the total flux.
  //

  //
  // First add the left and right state fluxes to out_flux[]:
  //
  PUtoFlux(left, Roe_UL, out_flux);
  PUtoFlux(right,Roe_UR, Roe_UL);
  for (int v=0;v<8;v++) 
    out_flux[v] += Roe_UL[v];

  //
  // Now add the contributions of the states between all the waves:
  //
  for (int iwave=0;iwave<7;iwave++) {
    out_flux[eqRHO] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][0];
    out_flux[eqMMX] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][1];
    out_flux[eqMMY] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][2];
    out_flux[eqMMZ] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][3];
    out_flux[eqBBY] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][4];
    out_flux[eqBBZ] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][5];
    out_flux[eqERG] -= Roe_strengths[iwave]*fabs(Roe_evalues[iwave])
      *Roe_right_evecs[iwave][6];
  }

  //
  // Divide by 2 to get the correct flux
  //
  for (int v=0;v<8;v++)
    out_flux[v] *= 0.5;


#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_MHD_CV::calculate_symmetric_flux ...returning.\n";
#endif //FUNCTION_ID

  return 0;
}
  

// ##################################################################
// ##################################################################



// **********************************************************************************
//  ROE SOLVER FOR ADIABATIC MHD, FROM CARGO & GALLICE, (1997) JCP, 136, 446
// **********************************************************************************

