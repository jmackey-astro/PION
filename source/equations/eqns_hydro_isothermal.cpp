///
/// \file equations_isothermal.cc
/// \author Jonathan Mackey
/// 
/// This file has the class definitions for the isothermal hydrodynamics equations.
///
/// History: 2009-10-23 written and tested.
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
/// - 2018.01.24 JM: worked on making SimPM non-global
///

#ifdef ISOTHERMAL_SOLVERS_ENABLED


#include "eqns_hydro_isothermal.h"
//#include <cmath>
using namespace std;


eqns_IsoEuler::eqns_IsoEuler(const int nv)
  : eqns_base(nv)
{
  cout <<"eqns_IsoEuler constructor.\n";
  eqAA  = PG;  // This is where we'll put the sound speed.
  eqAAA = ERG;
}

eqns_IsoEuler::~eqns_IsoEuler()
{}

void eqns_IsoEuler::PtoU(const double *p, ///< pointer to Primitive variables.
			 double *u,       ///< pointer to conserved variables.
			 const double     ///< Gas constant gamma.
			 )
{
  u[eqRHO] = p[eqRO];
  u[eqMMX] = p[eqRO]*p[eqVX];
  u[eqMMY] = p[eqRO]*p[eqVY];
  u[eqMMZ] = p[eqRO]*p[eqVZ];
  u[eqAAA] = p[eqRO]*p[eqAA]; /// sound speed advected as passive tracer
  return;
}

int eqns_IsoEuler::UtoP(
      const double *u, ///< pointer to conserved variables.
      double *p,       ///< pointer to Primitive variables.
      const double, ///< minimum temperature/pressure allowed
      const double   ///< unused (for gamma)
      )
{
  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX]/u[eqRHO];
  p[eqVY] = u[eqMMY]/u[eqRHO];
  p[eqVZ] = u[eqMMZ]/u[eqRHO];
  p[eqAA] = u[eqAAA]/u[eqRHO];
  if (p[eqRO]<TINYVALUE) {
    cout <<"eqns_IsoEuler::UtoP() samll or negative density!!! rho="<<p[eqRO];
    p[eqRO] = BASE_RHO:;
    cout <<"\t setting to base value. rho="<<p[eqRO]<<"\n";
    return 1;
  }
  return 0;
}

void eqns_IsoEuler::PtoFlux(const double *p,  ///< primitive vector
			    double *f,        ///< flux vector (output)
			    const double
			    )
{
   double u[eq_nvar];
   PtoU(p, u, 0);
   PUtoFlux(p,u,f);
}

void eqns_IsoEuler::PUtoFlux(const double *p, ///< pointer to Primitive variables.
			     const double *u, ///< pointer to conserved variables.
			     double *f  ///< Pointer to flux variable.
			     )
{
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*p[eqVX] +u[eqAAA]*p[eqAA];
  f[eqMMY] = u[eqMMX]*p[eqVY];
  f[eqMMZ] = u[eqMMX]*p[eqVZ];
  f[eqAAA] = u[eqMMX]*p[eqAA];
  return;
}

void eqns_IsoEuler::UtoFlux(const double*u, ///< Pointer to conserved variables state vector.
			    double*f,       ///< Pointer to flux variable state vector.
			    const double    ///< Unused Gas constant gamma.
			    )
{
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*u[eqMMX]/u[eqRHO] + u[eqAAA]*u[eqAAA]/u[eqRHO];
  f[eqMMY] = u[eqMMX]*u[eqMMY]/u[eqRHO];
  f[eqMMZ] = u[eqMMX]*u[eqMMZ]/u[eqRHO];
  f[eqAAA] = u[eqMMX]*u[eqAAA]/u[eqRHO];
  return;
}

/// Returns Internal Energy, which for isothermal equations is just p_g/rho=a^2
double eqns_IsoEuler::eint(const double *p, ///< Primitive State Vector.
				const double ///< gas EOS gamma.
				)
{
  return p[eqAA]*p[eqAA];
}

/// Returns Total Energy, for isothermal is 0.5*rho*v^2 +rho*a^2
double eqns_IsoEuler::Etot(const double *p, ///< State Vector.
				const double ///< gas EOS gamma.
				)
{
  return p[eqRO]*( 0.5*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]) +
		   p[eqAA]*p[eqAA]);
}

/// Returns Total Pressure (per unit Volume), which for hydro is just P_g
double eqns_IsoEuler::Ptot(const double *p, ///< Primitive State Vector.
				const double ///< gas EOS gamma.
				)
{
  return p[eqRO]*p[eqAA]*p[eqAA];
}


double eqns_IsoEuler::AdiabaticRho(const double pr, ///< New to Old pressure ratio
				   const double ri, ///< Old Density
				   const double     ///< gas EOS gamma.
				   )
{
  return ri*pr;
}

#endif // ISOTHERMAL_SOLVERS_ENABLED
