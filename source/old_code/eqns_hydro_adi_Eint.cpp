///
/// \file eqns_hydro_adi_Eint.cc
///
/// \author Jonathan Mackey
///
/// Declaration of the adiabatic hydrodynamics equations class, which 
/// uses internal energy as a 'conserved' variable instead of the 
/// total energy.
///
/// Modifications:\n
/// - 2010.12.28 JM: Created (moved from eqns_hydro_adiabatic.cc)
/// - 2010.12.30 JM: Can now integrate both internal and total energy.
/// - 2010.12.31 JM: updated logic in correcting negative pressures in
///   UtoP().  I think it's better now.
///
#include "../defines/functionality_flags.h"
#ifdef INCLUDE_EINT_ADI_HYDRO

#include "eqns_hydro_adi_Eint.h"
#include "eqns_hydro_adiabatic.h"
#include "eqns_base.h"
#include "../global.h"
using namespace std;


// --------------------------------------------------------------------
// ------------- EULER EQUATIONS WITH **INTERNAL ENERGY** -------------
// --------------------------------------------------------------------
eqns_Euler_Eint::eqns_Euler_Eint(int nv)
  : eqns_base(nv),
    eqns_Euler(nv)
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint:eqns_Euler_Eint ...starting.\n";
#endif //FUNCTION_ID
  //
  // Don't need to do anything much b/c everything is set in the Euler
  // equations class.
  //
#ifdef EINT_ETOT_PARALLEL
  eqEINT = 5;
  if (nv<6) rep.error("Set extra variable for Eint!",nv);
#endif // EINT_ETOT_PARALLEL

  cout <<"Internal Energy version of Euler Equations class!\n";
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint:eqns_Euler_Eint ...returning.\n";
#endif //FUNCTION_ID
}


eqns_Euler_Eint::~eqns_Euler_Eint()
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::~eqns_Euler_Eint ...starting.\n";
#endif //FUNCTION_ID
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::~eqns_Euler_Eint ...returning.\n";
#endif //FUNCTION_ID
}

void eqns_Euler_Eint::PtoU(const double* p,
			   double* u,
			   const double gamma)
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::PtoU ...starting.\n";
#endif //FUNCTION_ID
  // First rho
  u[eqRHO] = p[eqRO];
  // Second rho*u
  u[eqMMX] = p[eqRO]*p[eqVX];
  u[eqMMY] = p[eqRO]*p[eqVY];
  u[eqMMZ] = p[eqRO]*p[eqVZ];

#ifdef EINT_ETOT_PARALLEL
  //
  // Third E_int = p/(g-1), E_tot = E_int +0.5*rho*v^2
  //
  u[eqEINT] = p[eqPG]/(gamma-1.0);
  u[eqERG]  = u[eqEINT] +
    0.5*p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]);

#else // not EINT_ETOT_PARALLEL
  //
  // Third E_int = p/(g-1)
  //
  u[eqERG] = p[eqPG]/(gamma-1.0);
#endif // EINT_ETOT_PARALLEL

#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}

int eqns_Euler_Eint::UtoP(const double *u,
			  double *p,
			  const double gamma)
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::UtoP ...starting.\n";
#endif //FUNCTION_ID
  int err=0;
  static long int ct_pg=0, ct_rho=0;
  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX]/u[eqRHO];
  p[eqVY] = u[eqMMY]/u[eqRHO];
  p[eqVZ] = u[eqMMZ]/u[eqRHO];

#ifdef        EINT_ETOT_PARALLEL
  p[eqPG] = (gamma-1.0)*
    (u[eqERG] -0.5*p[eqRO]*
     (p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]));

  p[eqEINT] = (gamma-1.0)*u[eqEINT];

  bool test=false;
  if (p[eqEINT] <=0.0) {
    test=true;
    if (ct_pg<1000) {
      cout <<"(eqns_Euler_Eint::UtoP) negative pressure from internal";
      cout <<" energy!...p="<<p[eqEINT];
    }
    p[eqEINT]= max(p[eqPG],BASEPG*SimPM.RefVec[PG]);
    if (ct_pg<1000) {
      cout <<" ... will replace with max(p[eqPG],basevalue)=";
      cout <<p[eqEINT]<<"\n";
    }
    err += 1;
  }
  if (p[eqPG] <=0.0) {
    test=true;
    if (ct_pg<1000) {
      cout <<"(eqns_Euler_Eint::UtoP) negative pressure...p="<<p[eqPG];
      cout <<" ... will replace with Eint*(gamma-1) or basevalue="<<p[eqEINT]<<"\n";
    }
    p[eqPG] = p[eqEINT];
  }

  if (test) {
    ct_pg ++;
    if (ct_pg<1000) {
      cout <<"(eqns_Euler_Eint::UtoP) negative pressure...fixed: p="<<p[eqPG];
      cout <<", u[eqERG]="<< u[eqERG]<<", rho v^2=";
      cout << p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]);
#ifdef TESTING
      cout <<"NEG.PRES.CELL:";CI.print_cell(dp.c);
#endif
    }
  }
  
#else  // not EINT_ETOT_PARALLEL
  //
  // Here we are calculating pressure from internal energy, so it
  // should always be positive.
  //
  p[eqPG] = (gamma-1.0)*u[eqERG];
  if (p[eqPG] <=0.0) {
    if (ct_pg<1000) {
      ct_pg ++;
      cout <<"(eqns_Euler_Eint::UtoP) negative pressure...p="<<p[eqPG];
      cout <<", u[eqERG]="<< u[eqERG]<<", rho v^2=";
      cout << p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]);
      cout << "... correcting\n"; 
#ifdef TESTING
      cout <<"NEG.PRES.CELL:";CI.print_cell(dp.c);
#endif
    }
    p[eqPG] = BASEPG*SimPM.RefVec[PG]; // BASEPG is defined in global.h
    err +=1;
  }
#endif //     EINT_ETOT_PARALLEL

  //
  // Check for negative density and correct if it is found.
  // Note this is usually fatal to the simulation...
  //
  if (p[eqRO] <=0.0) {
    if (ct_rho<1000) {
      ct_rho ++;
      cout <<"(eqns_Euler_Eint::UtoP) negative density!  ";
      rep.printVec("u",u,eq_nvar); rep.printVec("p",p,eq_nvar);
#ifdef TESTING
      cout <<"NEG.DENS.CELL:";CI.print_cell(dp.c);
#endif
    }
    p[eqRO] = BASEPG*SimPM.RefVec[RO]; //100*MACHINEACCURACY;
    err += 1;
  }


  //  if (err) cout <<"Error in UtoP()!!\n";
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::UtoP ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


void eqns_Euler_Eint::PUtoFlux(const double *p, const double *u, double *f)
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::PUtoFlux ...starting.\n";
#endif //FUNCTION_ID
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*p[eqVX] +p[eqPG];
  f[eqMMY] = u[eqMMX]*p[eqVY];
  f[eqMMZ] = u[eqMMX]*p[eqVZ];

#ifdef        EINT_ETOT_PARALLEL
  //
  // Total Energy Flux u(E+p) = u (rho*u^2/2 +g*p/(g-1)) for ideal gas
  //
  f[eqERG] = p[eqVX]*(u[eqERG]+p[eqPG]);
  //
  // Internal energy flux is just Eint*vx
  //
  f[eqEINT] = p[eqVX]*u[eqEINT];

#else  // not EINT_ETOT_PARALLEL
  //
  // Energy Flux for internal energy is vx*Eint
  // The PdV term is calculated as a source term.
  //
  f[eqERG] = p[eqVX]*u[eqERG];
 
#endif //     EINT_ETOT_PARALLEL

 //      cout <<"Fluxes: "<<c->F[0]<<"  "<<c->F[1]<<"  "<<c->F[2]<<"\n";

#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::PUtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}

void eqns_Euler_Eint::UtoFlux(const double *u, double *f, const double gamma)
{
#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::UtoFlux ...starting.\n";
#endif //FUNCTION_ID
  //
  // This is easier than the total energy Euler equations.
  //
  f[eqRHO] = u[eqMMX];
  f[eqMMY] = u[eqMMX]*u[eqMMY]/u[eqRHO];
  f[eqMMZ] = u[eqMMX]*u[eqMMZ]/u[eqRHO];

#ifdef        EINT_ETOT_PARALLEL
  double pg = (gamma-1.) *
    (u[eqERG] -
     (u[eqMMX]*u[eqMMX] +u[eqMMY]*u[eqMMY] +u[eqMMZ]*u[eqMMZ])
     /(2.0*u[eqRHO]));
  f[eqMMX] = u[eqMMX]*u[eqMMX]/u[eqRHO] + pg;
  f[eqERG] = u[eqMMX]*(u[eqERG]+pg)/u[eqRHO];
  f[eqEINT] = u[eqMMX]*u[eqEINT]/u[eqRHO];

#else  // not EINT_ETOT_PARALLEL
  f[eqMMX] = u[eqMMX]*u[eqMMX]/u[eqRHO] +(gamma-1.0)*u[eqERG];
  f[eqERG] = u[eqMMX]*u[eqERG]/u[eqRHO];
#endif //     EINT_ETOT_PARALLEL

#ifdef FUNCTION_ID
  cout <<"eqns_Euler_Eint::UtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// --------------------------------------------------------------------
// ------------- EULER EQUATIONS WITH **INTERNAL ENERGY** -------------
// --------------------------------------------------------------------


#endif // if INCLUDE_EINT_ADI_HYDRO
