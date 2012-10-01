/** \file eqns_hydro_adiabatic.cc
 * 
 * \brief Class definition for eqns_Euler
 * \author Jonathan Mackey
 * 
 * This file contains the class definitions for the Euler Equations.
 * 
 * Modifications:
 *  - 2007-10-16 Moved the equations classes in here from global.cc
 *  - 2009-10-20 renamed eqns_hydro_adiabatic.cc and moved all other equations classes to their own files.
 * */
///
/// - 2010.07.26 JM: added static int counters to suppress negative
///    pressure/density messages after a certain number have been reported.
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///
/// - 2010.11.21 JM: Changed UtoP() so that it corrects both pressure
///   and density before returning!
///
/// - 2010.12.21 JM: Added new equations class which uses internal
///   energy and not total energy as a conserved variable.
///
/// - 2010.12.23 JM: Added SetAvgState() function to eqns_Euler class.
///                  (fixed bug in this 27/12)
///
/// - 2010.12.28 JM: moved internal energy Euler eqns to own file.
///
/// - 2011.01.16 JM: Added ifdef to set pressure so that temperature 
///   is 10K when UtoP gives a negative pressure.  This can be used 
///   for stellar wind. (modified 2011.01.18 so Tmin=10K). modified
///   2011.02.17 JM so negative density is fixed before pressure, so 
///   that a "meaningful" temperature can be calculated.
/// - 2011.04.06 JM: minor mods to UtoP() function again in correcting
///    for negative and "small" temperatures.
/// - 2011.04.15 JM: UtoP() again -- added recalculation of all prim.
///    vars if a negative density is encountered.
///
#include "../global.h"
#include "eqns_hydro_adiabatic.h"
using namespace std;



// Member function definitions for eqns_Euler class
eqns_Euler::eqns_Euler(int nv)
  : eqns_base(nv)
{
  cout <<"(eqns_Euler::eqns_Euler) Setting up Euler Equations Class.\n";
//  cout <<"\tVector lengths: "<<eqnvar<<"\n";
  if(eq_nvar<5) rep.error("eqns_Euler initialised with eq_nvar<5.",eq_nvar);
  //cout <<"Setting Flux functions to X-dir: pu2f(), u2f()\n";
//  pu2flux = &eqns_Euler::pu2f; 
//  u2flux = &eqns_Euler::u2f;

  SetDirection(XX);
  eqRO = RO; eqPG = PG;
  eqRHO=RHO; eqERG=ERG;
  //  cout <<"(eqns_Euler::eqns_Euler) Done.\n";
}

eqns_Euler::~eqns_Euler()
{
//  cout <<"(eqns_Euler::~eqns_Euler) Deleting Euler Equations Class.\n";
}

void eqns_Euler::PtoU(const double* p, double* u, const double gamma)
{
  // First rho
  u[eqRHO] = p[eqRO];
  // Second rho*u
  u[eqMMX] = p[eqRO]*p[eqVX];
  u[eqMMY] = p[eqRO]*p[eqVY];
  u[eqMMZ] = p[eqRO]*p[eqVZ];
  // Third E = p/(g-1) +rho*V^2/2
  u[eqERG] = p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ])/2.0 +p[eqPG]/(gamma-1.);
  return;
}

int eqns_Euler::UtoP(const double *u, double *p, const double gamma)
{
  int err=0;
#ifndef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  static long int ct_pg=0;
#endif
  static long int ct_rho=0;

  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX]/u[eqRHO];
  p[eqVY] = u[eqMMY]/u[eqRHO];
  p[eqVZ] = u[eqMMZ]/u[eqRHO];
  p[eqPG] = (gamma-1.0) *
    (u[eqERG] -p[eqRO]*
     (p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ])/2.0);


  //
  // First check for negative density, and fix it if present.
  // Note this is usually fatal to a simulation, so we print out messages so
  // that we know this is what has happened (for the first 1000 instances).
  //
  if (p[eqRO] <=0.0) {
    if (ct_rho<1000) {
      ct_rho ++;
      cout <<"(eqns_Euler::UtoP) negative density!  ";
      rep.printVec("u",u,eq_nvar);
      rep.printVec("p",p,eq_nvar);
#ifdef TESTING
      cout <<"NEG.DENS.CELL:";CI.print_cell(dp.c);
#endif
    }
    // reset all variables because a negative density will change the sign of 
    // all of the velocities!
    p[eqRO] = BASE_RHO*SimPM.RefVec[RO];
    p[eqVX] = u[eqMMX]/p[eqRO];
    p[eqVY] = u[eqMMY]/p[eqRO];
    p[eqVZ] = u[eqMMZ]/p[eqRO];
    p[eqPG] = (gamma-1.0)*(u[eqERG]-p[eqRO]*(p[eqVX]*p[eqVX]+p[eqVY]*p[eqVY]+p[eqVZ]*p[eqVZ])/2.0);
    err += 1;
  }

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // First check for negative pressure, fixing it if needed.
  // Then if there was no negative pressure, check for pressure
  // being too low, and fix that.
  //
  if (p[eqPG] <=0.0) {
    //
    // Set minimum temperature to be 10K
    //
    if (MP) {
      //if (MP->Temperature(p,gamma) <SimPM.EP.MinTemperature) {
#ifdef TESTING
      cout <<"fixing negative pressure from p="<<p[eqPG]<<" to p=";
#endif
      MP->Set_Temp(p,SimPM.EP.MinTemperature,gamma);
#ifdef TESTING
      cout <<p[eqPG]<<"\n";
#endif
    }
    else {
      //
      // If not, assume mu*m_p=2.34e-24g, and set p=k*rho
      //
#ifdef TESTING
      cout <<"fixing negative pressure from p="<<p[eqPG]<<" to p=";
#endif
      p[eqPG] = 5.8974e8*p[eqRO];
#ifdef TESTING
      cout <<p[eqPG]<<"\n";
#endif
    }
  }
  else if (MP && (MP->Temperature(p,gamma) <SimPM.EP.MinTemperature)) {
    //
    // If we have microphysics, just set T=10K.
    //
#ifdef TESTING
    cout <<"fixing low pressure from p="<<p[eqPG]<<" to p=";
#endif
    MP->Set_Temp(p,SimPM.EP.MinTemperature,gamma);
#ifdef TESTING
    cout <<p[eqPG]<<"\n";
#endif
  }

#else // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  if (p[eqPG] <=0.0) {
    if (ct_pg<1000) {
      ct_pg ++;
      cout <<"(eqns_Euler::UtoP) negative pressure...p="<<p[eqPG];
      cout <<", u[eqERG]="<< u[eqERG]<<", 0.5*rho*v^2=";
      cout << 0.5*p[eqRO]*
	(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ]);
      cout << "... correcting\n"; 
#ifdef TESTING
      cout <<"NEG.PRES.CELL:";CI.print_cell(dp.c);
#endif
    }
    p[eqPG] = BASEPG*SimPM.RefVec[PG]; //100*MACHINEACCURACY;
    //err +=1;
  }
#endif // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  //  if (err) cout <<"Error in UtoP()!!\n";
  return err;
}

double eqns_Euler::chydro(const double *p, ///< Pointer to primitive variables.
			  const double g   ///< Gas constant gamma.
			  )
{
  return(sqrt(g*p[eqPG]/p[eqRO]));
}

int eqns_Euler::HydroWave(int lr, const double pp, const double *prewave, double *u, const double gamma) 
{
  double pratio = pp/prewave[eqPG];
  // First set the appropriate prewave sound speed.
  double c0;
  c0 = sqrt(gamma*prewave[eqPG]/prewave[eqRO]);
  // Variable from the solver are gamma, prewave (points to *left/ *right), pstar
  if (pratio<1) { // Rarefaction
    // We have |u*-u_0| = 2c_0/(g-1) *(1- (p* /p_0)^((g-1)/2g) )
    // In order to disentangle the sign ambiguity we need to know if it is a left of right wave
    *u = 2.*c0/(gamma-1.)*(1 - exp((gamma-1.)/2./gamma*log(pratio)));
    if (lr == XN) // Left moving rarefaction has u*-u0 >0
      *u = prewave[eqVX] + (*u);
    else
      *u = prewave[eqVX] - (*u);
  }
  else if (pratio>1) { // Shock
    // Here we have |u*-u0| = c0/sqrt(g(g-1)/2) *(pr-1)/sqrt(1+(g+1)/(g-1)pr)
    // In order to disentangle the sign ambiguity we need to know if it is a left of right wave
    *u = c0 *(pratio-1.) /sqrt(gamma*(gamma-1.)/2. *(1.+pratio*(gamma+1.)/(gamma-1.)));
    if (lr == XN) // Left moving shock has u*-u0 <0
      *u = prewave[eqVX] - (*u);
    else               // Right moving has u*-u0 >0
      *u = prewave[eqVX] + (*u);
  }
  else {
    // In this case there is no wave, because there is no p*/p0=1.
    *u = prewave[eqVX];
  }
  return(0);
}

						  
int eqns_Euler::HydroWaveFull(int lr, const double pp, const double *prewave, double *u, double *rho, const double gamma)
{
  double pratio = pp/prewave[eqPG];
  // First get the appropriate velocity u*
  int err = HydroWave(lr, pp, prewave, u, gamma);
  if (err !=0) {
    cerr << "(riemann::hydro_wave_full) hydro wave call didn't work, exiting." << "\n";
    return(1);
  }
  // All that is left to do is solve for the density.
  // We don't care if it's left or right moving, b/c rho only depends on P, and both are scalar.
  if (pratio<1) // Here we can use the fact that rarefactions are adiabatic (p/rho^gamma conserved).
    *rho = prewave[eqRO] *exp(log(pratio)/gamma);
  else if (pratio>1) // Here we use the shock jump conditions from Hirsch p.206 (vol.2).
    *rho = prewave[eqRO] *(1+pratio*(gamma+1)/(gamma-1.)) /((gamma+1.)/(gamma-1.) +pratio);
  else *rho = prewave[eqRO]; // There is no wave.
  return(0);
}


void eqns_Euler::PUtoFlux(const double *p, const double *u, double *f)
{
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*p[eqVX] +p[eqPG];
  f[eqMMY] = u[eqMMX]*p[eqVY];
  f[eqMMZ] = u[eqMMX]*p[eqVZ];
  // Energy Flux u(E+p) = u (rho*u^2/2 +g*p/(g-1)) for ideal gas
  f[eqERG] = p[eqVX]*(u[eqERG]+p[eqPG]);
  //      cout <<"Fluxes: "<<c->F[0]<<"  "<<c->F[1]<<"  "<<c->F[2]<<"\n";
  return;
}

void eqns_Euler::UtoFlux(const double *u, double *f, const double gamma)
{
  double pg = (gamma-1.) *(u[eqERG] -(u[eqMMX]*u[eqMMX] +u[eqMMY]*u[eqMMY] +u[eqMMZ]*u[eqMMZ])/(2.*u[eqRHO]));
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*u[eqMMX]/u[eqRHO] + pg;
  f[eqMMY] = u[eqMMX]*u[eqMMY]/u[eqRHO];
  f[eqMMZ] = u[eqMMX]*u[eqMMZ]/u[eqRHO];
  f[eqERG] = u[eqMMX]*(u[eqERG]+pg)/u[eqRHO];
  return;
}


///  Returns Enthalpy (per unit mass), given primitive variable vector. 
double eqns_Euler::Enthalpy(const double *p, ///< Primitive State Vector.
			    const double g   ///< gas EOS gamma.
			    )
{
  //cout <<"Enthalpy!\n";
  return (0.5*(p[eqVX]*p[eqVX]+p[eqVY]*p[eqVY]+p[eqVZ]*p[eqVZ]) +g*p[eqPG]/(g-1.0)/p[eqRO]);
}

///  Returns Internal Energy (per unit mass, so 'Temperature'), given primitive variable vector. 
double eqns_Euler::eint(const double *p, ///< Primitive State Vector.
			const double g   ///< gas EOS gamma.
			)
{
  return p[eqPG]/(g-1.)/p[eqRO];
}


///  Returns Total Energy (per unit volume), given primitive variable vector. 
double eqns_Euler::Etot(const double *p, ///< State Vector.
			const double g   ///< gas EOS gamma.
			) 
{
  return p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ])/2. +p[eqPG]/(g-1.);
}

///  Returns Total Pressure (per unit Volume), given primitive variable vector. 
double eqns_Euler::Ptot(const double *p, ///< Primitive State Vector.
			const double    ///< gas EOS gamma.
			) 
{
  return p[eqPG];
}

///  Given a pressure ratio and initial density, calculate adiabatic final density.
double eqns_Euler::AdiabaticRho(const double pr, ///< New to Old pressure ratio
				const double ri, ///< Old Density
				const double g ///< gas EOS gamma.
				)
{
  return ri*exp(log(pr)/g);
}

void eqns_Euler::SetAvgState(const double *ms,  ///< Mean Prim. var. state vector
			     const double g ///< Gas constant gamma.
			     )
{
  for (int v=0; v<eq_nvar; v++)
    eq_refvec[v] = ms[v];
  //
  // Now reset reference velocities to be 1/10th of the sound speed.
  //
  double refvel = chydro(eq_refvec, g);
  eq_refvec[eqVX] = eq_refvec[eqVY] = eq_refvec[eqVZ] = 0.1*refvel;

  return;
}

