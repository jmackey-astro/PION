///
/// \file Roe_Hydro_PrimitiveVar_solver.cc
/// \author Jonathan Mackey
///
/// Takes in a left and right state in primitive variables, and
/// returns the intermediate state Pstar from which the flux across
/// the interface can be calculated.  This uses the Roe average for
/// the mean state, and is a linear Riemann solver.
///
/// References:
///  - Toro (1999) Riemann Solvers Textbook, Chapter 11.
///
/// History:
/// - 2010-12.22 JM: Moved functions from flux_hydro_adiabatic.h/.cc
///   and generated new class.
/// - 2011.01.19 JM: Simplified calculation of rho*L/R.
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer 
///    fluxes are dealt with by the flux-solver classes.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "Roe_Hydro_PrimitiveVar_solver.h"

using namespace std;

Riemann_Roe_Hydro_PV::Riemann_Roe_Hydro_PV(
      const int nv,///< Length of State Vectors, nvar
      const double g ///< Gamma for state vector.
      )
  : eqns_base(nv), eqns_Euler(nv), rs_nvar(5)
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::Riemann_Roe_Hydro_PV ...starting.\n";
#endif //FUNCTION_ID
  //
  // eq_gamma, eq_nvar are defined in eqns_base class
  //
  eq_gamma = g;
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::Riemann_Roe_Hydro_PV ...returning.\n";
#endif //FUNCTION_ID
  
}

Riemann_Roe_Hydro_PV::~Riemann_Roe_Hydro_PV()
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::~Riemann_Roe_Hydro_PV ...starting.\n";
#endif //FUNCTION_ID
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::~Riemann_Roe_Hydro_PV ...returning.\n";
#endif //FUNCTION_ID
}


int Riemann_Roe_Hydro_PV::Roe_prim_var_solver(
      const pion_flt *rpv_left,
      const pion_flt *rpv_right,
      const double rpv_g,
      pion_flt *rpv_pstar
      )
{
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::Roe_prim_var_solver ...starting.\n";
#endif //FUNCTION_ID
  int err=0;

  //
  // Linearised Riemann solver using the Roe average state.  This
  // doesn't refer to a reference vector, in the hope that it will be
  // better for stellar wind shocks.
  //

#ifdef RSTESTING
  //
  // Make sure we have finite inputs.
  //
  for (int v=0;v<rs_nvar;v++) {
    if (!isfinite(rpv_left[v]) || !isfinite(rpv_right[v]) || !isfinite(rpv_g)) {
      rep.printVec(" left",rpv_left, rs_nvar);
      rep.printVec("right",rpv_right,rs_nvar);
      rep.printVec("pstar",rpv_pstar,rs_nvar);
      rep.error("NAN input states to linear solver!","NANANANANA");
    }
  }
#endif // RSTESTING

  //
  // First construct average state:
  //
  enum primitive eqHH = eqPG;
  double 
    rl = sqrt(rpv_left[eqRO]),
    rr = sqrt(rpv_right[eqRO]),
    lH = Enthalpy(rpv_left ,rpv_g),
    rH = Enthalpy(rpv_right,rpv_g),
    denom = 1.0/(rl+rr),
    a_mean =0.0, v2_mean = 0.0;
  pion_flt rpv_meanp[rs_nvar];
  rpv_meanp[eqRO] = rl*rr;
  rpv_meanp[eqVX] = (rl*rpv_left[eqVX]+rr*rpv_right[eqVX])*denom;
  rpv_meanp[eqVY] = (rl*rpv_left[eqVY]+rr*rpv_right[eqVY])*denom;
  rpv_meanp[eqVZ] = (rl*rpv_left[eqVZ]+rr*rpv_right[eqVZ])*denom;
  //
  // Get the enthalpy of the mean state: we'll convert back to a
  // pressure below.
  //
  // Enthalpy per unit mass: H= v*v/2 +g*p/(g-1)/rho
  // Hence Adiabatic sound speed a = sqrt[(H-v*v/2)*(g-1)]
  // Then pressure = rho*a*a/g
  //
  rpv_meanp[eqHH] = (rl*lH + rr*rH)*denom;
  
  //
  // v2_mean = V^2 for the mean state.
  // a_mean  = adiabatic sound speed for mean state.
  //
  v2_mean = rpv_meanp[eqVX]*rpv_meanp[eqVX] 
    +rpv_meanp[eqVY]*rpv_meanp[eqVY] 
    +rpv_meanp[eqVZ]*rpv_meanp[eqVZ];
  a_mean = sqrt((rpv_g-1.0)*(rpv_meanp[eqHH]-0.5*v2_mean));
  rpv_meanp[eqPG] = rpv_meanp[eqRO]*a_mean*a_mean/rpv_g;

  //
  // Old code to see if the Roe-average pressure is different from the
  // pressure obtained from the Roe-average enthalpy (IT IS!).
  //
  // double temp=(rr*rpv_left[eqPG]+rl*rpv_right[eqPG])*denom;
  // if (fabs(temp/rpv_meanp[eqPG]-1.0)>1.e-8) {
  //   cout <<"rpv_meanp(H)="<<rpv_meanp[eqPG];
  //   cout <<"  roe-avg="<<temp<<"\n";
  // }

  //
  // Now we copy the code from riemann_Euler::linear_solver(), with a
  // few modifications to get rid of the average state.
  //
  // TESTER FOR THE EIGENVALUE SIGNS... (u_av-c_av, u_av, u_av+c_av)
  //
  if (rpv_meanp[eqVX]-a_mean >=0.) {
    //
    // if(u_av-c_av>0) pstar = left state
    //
    for (int i=0; i<rs_nvar; i++) rpv_pstar[i] = rpv_left[i];
  }
  else if (rpv_meanp[eqVX]+a_mean <=0.) {
    //
    // else if(u_av+c_av<=0) pstar = right state
    //
    for (int i=0; i<rs_nvar; i++) rpv_pstar[i] = rpv_right[i];
  }
  else {
    //
    // else we are in the starred region, so p_star, u_star give the
    // solution.
    //
    rpv_pstar[eqPG] =
      0.5*(rpv_left[eqPG]+rpv_right[eqPG]
	   -rpv_meanp[eqRO]*a_mean*(rpv_right[eqVX]-rpv_left[eqVX]));
    rpv_pstar[eqVX] = 
      0.5*(rpv_left[eqVX]+rpv_right[eqVX]
	   -(rpv_right[eqPG]-rpv_left[eqPG])/rpv_meanp[eqRO]/a_mean);

    //
    // Now just need to determine rho_star, based on the sign of
    // vx* (which determines where the contact discontinuity is.
    //
    if (rpv_pstar[eqVX] >0.0) {
      //
      // if(u*>0) pstar = left starred state
      //
      rpv_pstar[eqRO] =	rpv_left[eqRO]+
     	rpv_meanp[eqRO]*(rpv_left[eqVX]-rpv_pstar[eqVX])/a_mean; // rho_(*L)
    }
    else {
      //
      // else if(u*<0) pstar = right starred state
      //
      rpv_pstar[eqRO] = rpv_right[eqRO]+
     	rpv_meanp[eqRO]*(rpv_pstar[eqVX]-rpv_right[eqVX])/a_mean; // rho_(*R)
    }

    //
    // Need to assign values to v_y,v_z if present.
    // They only change across the contact discontinuity.
    //
    if (rpv_pstar[eqVX]>0.0) {
      rpv_pstar[eqVY] = rpv_left[eqVY];		    
      rpv_pstar[eqVZ] = rpv_left[eqVZ];
    }
    else {
      rpv_pstar[eqVY] = rpv_right[eqVY];		    
      rpv_pstar[eqVZ] = rpv_right[eqVZ];
    }
    //  cout << "(riemann_Euler::solve) Success!" << "\n";
  }

  //
  // Now we are finished, so return 0 for success.
  //
#ifdef FUNCTION_ID
  cout <<"Riemann_Roe_Hydro_PV::Roe_prim_var_solver ...returning.\n";
#endif //FUNCTION_ID
  return err;

}


