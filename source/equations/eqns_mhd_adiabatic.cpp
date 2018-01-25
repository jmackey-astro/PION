///
/// \file eqns_mhd_adiabatic.cc
/// 
/// \brief Class definition for eqns_mhd_ideal
/// \author Jonathan Mackey
/// 
/// This file contains the class definitions for the Ideal MHD Equations,
/// with the GLM divergence cleaning extensions.
/// 
/// Modifications:
/// - 2007-10-16 Moved the equations classes in here from global.cc
/// - 2009-10-20 renamed eqns_mhd_adiabatic.cc and moved all other equations classes to their own files.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.23 JM: Added SetAvgState() function to eqns_mhd_ideal
///    class.
/// - 2011.01.18 JM: Added BASE_RHO parameter for correcting negative
///   densities.
/// - 2011.04.15 JM: UtoP() again -- added recalculation of all prim.
///    vars if a negative density is encountered.  Also added the
///    ifdef for setting negative pressure to a fixed temperature.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "constants.h"

#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/microphysics_base.h"
#include "eqns_mhd_adiabatic.h"
#include <iostream>
using namespace std;


/*******************************************************************/
// Member function definitions for eqns_mhd_ideal class

// ##################################################################
// ##################################################################

eqns_mhd_ideal::eqns_mhd_ideal(int nv)
  : eqns_base(nv)
{
#ifdef FUNCTION_ID
  cout <<"eqns_mhd_ideal::eqns_mhd_ideal ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  cout <<"(eqns_mhd_ideal::eqns_mhd_ideal) Setting up Ideal MHD Equations class.\n";
#endif
  if(eq_nvar<8) {cerr<<"\tError!! Class expects (at least) 8 MHD variables.\n"; exit(1);}

#ifdef FUNCTION_ID
  cout <<"eqns_mhd_ideal::eqns_mhd_ideal ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################

eqns_mhd_ideal::~eqns_mhd_ideal()
{
#ifdef TESTING
  cout <<"(eqns_mhd_ideal::~eqns_mhd_ideal) Deleting  Ideal MHD Equations class.\n";
#endif
}


// ##################################################################
// ##################################################################

void eqns_mhd_ideal::PtoU(
      const pion_flt *p,
      pion_flt *u,
      const double gamma
      )
{
  u[eqRHO] = p[eqRO];
  u[eqMMX] = p[eqRO]*p[eqVX];
  u[eqMMY] = p[eqRO]*p[eqVY];
  u[eqMMZ] = p[eqRO]*p[eqVZ];
  u[eqBBX] = p[eqBX];
  u[eqBBY] = p[eqBY];
  u[eqBBZ] = p[eqBZ];
  // E = p/(g-1) +rho*V^2/2 + B^2/2
  u[eqERG] = (p[eqRO]*(p[eqVX]*p[eqVX]+p[eqVY]*p[eqVY]+p[eqVZ]*p[eqVZ])/2.) +(p[eqPG]/(gamma-1.)) +((u[eqBBX]*u[eqBBX] +u[eqBBY]*u[eqBBY] +u[eqBBZ]*u[eqBBZ])/2.);
  //cout <<"gamma="<<gamma<<"\n";
  return;
}


// ##################################################################
// ##################################################################

int eqns_mhd_ideal::UtoP(
      class SimParams &SimPM, ///< pointer to simulation parameters
      const pion_flt *u,
      pion_flt *p,
      const double gamma
      )
{
#ifndef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  static long int ct_pg=0;
#endif
  static long int ct_rho=0;

  p[eqRO] = u[eqRHO];
  p[eqVX] = u[eqMMX]/u[eqRHO];
  p[eqVY] = u[eqMMY]/u[eqRHO];
  p[eqVZ] = u[eqMMZ]/u[eqRHO];
  p[eqPG] = (gamma-1) *(u[eqERG] 
		      -p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ])/2. 
		      -(u[eqBBX]*u[eqBBX] +u[eqBBY]*u[eqBBY] +u[eqBBZ]*u[eqBBZ])/2.);
  p[eqBX] = u[eqBBX];
  p[eqBY] = u[eqBBY];
  p[eqBZ] = u[eqBBZ];

  
  //
  // First check for negative density, and fix it if present.
  // Note this is usually fatal to a simulation, so we print out messages so
  // that we know this is what has happened (for the first 1000 instances).
  //
  int err=0;
  if (p[eqRO] <=0.0) {
    if (ct_rho<1000) {
      ct_rho ++;
      cout <<"(eqns_mhd_ideal::UtoP) negative density!  ";
      rep.printVec("u",u,eq_nvar);
      rep.printVec("p",p,eq_nvar);
#ifdef TESTING
      cout <<"NEG.DENS.CELL:";CI.print_cell(dp.c);
#endif
    }
    // reset all variables because a negative density will change the sign of 
    // all of the velocities!
    p[eqRO] = BASE_RHO*SimPM.RefVec[RO];
    p[eqVX] *= u[eqRHO]/p[eqRO];
    p[eqVY] *= u[eqRHO]/p[eqRO];
    p[eqVZ] *= u[eqRHO]/p[eqRO];
    p[eqPG] = (gamma-1)
        *(u[eqERG] -p[eqRO]*(p[eqVX]*p[eqVX] +p[eqVY]*p[eqVY] +p[eqVZ]*p[eqVZ])/2. 
                     -(u[eqBBX]*u[eqBBX] +u[eqBBY]*u[eqBBY] +u[eqBBZ]*u[eqBBZ])/2.);
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
    //cout <<"UtoP() mhd set-neg-press-to-fixed-T.  P<0\n";
    if (MP) {
      MP->Set_Temp(p,SimPM.EP.MinTemperature,gamma);
    }
    else {
      //
      // If not, assume mu*m_p=2.34e-24g, and set p=k*rho
      //
      p[eqPG] = 5.8974e8*p[eqRO];
    }
  }
  else if (MP && (MP->Temperature(p,gamma) <SimPM.EP.MinTemperature)) {
    //
    // If we have microphysics, just set T=10K.
    //
    //cout <<"UtoP() mhd set-neg-press-to-fixed-T.  T<Tmin\n";
    MP->Set_Temp(p,SimPM.EP.MinTemperature,gamma);
  }

#else // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  if (p[eqPG] <=0.) {
    if (ct_pg<1000) {
      ct_pg ++;
      cout <<"(eqns_mhd_ideal::UtoP) negative pressure...p="<<p[eqPG];
      cout <<", correcting, count="<<ct_pg<<"\n";
    }
    p[eqPG] =  BASEPG*SimPM.RefVec[PG];
    err += 1;
  }
#endif // don't SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  return err;
}


// ##################################################################
// ##################################################################


double eqns_mhd_ideal::chydro(
      const pion_flt *p,   ///< Pointer to primitive variables.
      const double gamma ///< Gas constant gamma.
      )
{
  return sqrt(gamma*p[eqPG]/p[eqRO]);
}


// ##################################################################
// ##################################################################


/* Calculate fast magnetic wavespeed, in direction we are looking in. */
double eqns_mhd_ideal::cfast(
      const pion_flt *p,
      const double gamma
      )
{
  //cout <<"cfast!\n";
  double ch = chydro(p, gamma);
  double temp1 = ch*ch + (p[eqBX]*p[eqBX] +p[eqBY]*p[eqBY] +p[eqBZ]*p[eqBZ])/p[eqRO];
  double temp2 = 4.*ch*ch*p[eqBX]*p[eqBX]/p[eqRO];
  if ((temp2=temp1*temp1-temp2) <MACHINEACCURACY) temp2=MACHINEACCURACY; // This is as good as the computer can get.
  return( sqrt( (temp1 + sqrt(temp2))/2.) );
}


// ##################################################################
// ##################################################################


double eqns_mhd_ideal::cfast_components(
      const double cfRO,
      const double cfPG,
      const double cfBX,
      const double cfBY,
      const double cfBZ,
      const double g
      )
{
  double ch = sqrt(g*cfPG/cfRO);
  double temp1 = ch*ch + (cfBX*cfBX +cfBY*cfBY +cfBZ*cfBZ)/cfRO;
  double temp2 = 4.*ch*ch*cfBX*cfBX/cfRO;
  //
  // This subtraction has to be done carefully to avoid sqrt(-ve number).
  // This is as good as the computer can get.
  //
  if ((temp2=temp1*temp1-temp2) <MACHINEACCURACY) temp2=MACHINEACCURACY;
  return( sqrt( (temp1 + sqrt(temp2))/2.) );
}


// ##################################################################
// ##################################################################


// Calculate slow magnetic wavespeed.
double eqns_mhd_ideal::cslow(
      const pion_flt *p,
      const double gamma
      )
{
  double temp1;
  double temp2;
  double ch = chydro(p, gamma);
  temp1 = ch*ch + (p[eqBX]*p[eqBX] +p[eqBY]*p[eqBY] +p[eqBZ]*p[eqBZ])/p[eqRO];
  temp2 = 4.*ch*ch*p[eqBX]*p[eqBX]/p[eqRO];
  if ((temp2=temp1*temp1-temp2) <MACHINEACCURACY) temp2=MACHINEACCURACY; // This is as good as the computer can get.
  if ((temp2=temp1-sqrt(temp2)) <MACHINEACCURACY) temp2=MACHINEACCURACY; // Again, as good as the machine can get.
  /** \section Accuracy
   * Because I'm taking the square root of a difference of two numbers, 
   * the slow magnetic speed can only be determined accurately down to
   * the square root of the machine accuracy (\f$\simeq 3 \times 10^{-8}\f$ 
   * for double precision), so I have it hardcoded to never be smaller than this.
   * */
  return( sqrt(temp2/2.) );
}


// ##################################################################
// ##################################################################


void eqns_mhd_ideal::PUtoFlux(
      const pion_flt *p,
      const pion_flt *u,
      pion_flt *f
      )
{
  /** \section Equations
   * The equations for the flux are in Falle, Komissarov, Joarder, 1998, MNRAS,297,265.
   * Equation (2).
   * */
  double pm = (u[eqBBX]*u[eqBBX] +u[eqBBY]*u[eqBBY] +u[eqBBZ]*u[eqBBZ])/2.; // Magnetic pressure.
  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*p[eqVX] +p[eqPG] +pm -u[eqBBX]*u[eqBBX];
  f[eqMMY] = u[eqMMX]*p[eqVY] -u[eqBBX]*u[eqBBY];
  f[eqMMZ] = u[eqMMX]*p[eqVZ] -u[eqBBX]*u[eqBBZ];
  f[eqERG] = p[eqVX]*(u[eqERG]+p[eqPG]+pm) -u[eqBBX]*(p[eqVX]*u[eqBBX] +p[eqVY]*u[eqBBY] +p[eqVZ]*u[eqBBZ]);
  f[eqBBX] = 0.;
  f[eqBBY] = p[eqVX]*p[eqBY] - p[eqVY]*p[eqBX];
  f[eqBBZ] = p[eqVX]*p[eqBZ] - p[eqVZ]*p[eqBX];
  return;
}


// ##################################################################
// ##################################################################


void eqns_mhd_ideal::UtoFlux(
      const pion_flt *u,
      pion_flt *f,
      const double gamma
      )
{
  double pm = (u[eqBBX]*u[eqBBX] +u[eqBBY]*u[eqBBY] +u[eqBBZ]*u[eqBBZ])/2.; // Magnetic pressure.
  double pg = (gamma-1.) *(u[eqERG] -(u[eqMMX]*u[eqMMX] +u[eqMMY]*u[eqMMY] +u[eqMMZ]*u[eqMMZ])/(2.*u[eqRHO]) -pm);

  f[eqRHO] = u[eqMMX];
  f[eqMMX] = u[eqMMX]*u[eqMMX]/u[eqRHO] + pg +pm -u[eqBBX]*u[eqBBX];
  f[eqMMY] = u[eqMMX]*u[eqMMY]/u[eqRHO] - u[eqBBX]*u[eqBBY];
  f[eqMMZ] = u[eqMMX]*u[eqMMZ]/u[eqRHO] - u[eqBBX]*u[eqBBZ];
  f[eqERG] = u[eqMMX]*(u[eqERG]+pg+pm)/u[eqRHO] - u[eqBBX]*(u[eqMMX]*u[eqBBX]+u[eqMMY]*u[eqBBY]+u[eqMMZ]*u[eqBBZ])/u[eqRHO];
  f[eqBBX] = 0.;
  f[eqBBY] = (u[eqMMX]*u[eqBBY]-u[eqMMY]*u[eqBBX])/u[eqRHO];
  f[eqBBZ] = (u[eqMMX]*u[eqBBZ]-u[eqMMZ]*u[eqBBX])/u[eqRHO];
  return;
}


// ##################################################################
// ##################################################################


void eqns_mhd_ideal::rotate(
      pion_flt *vec, ///< State vector
      enum axes initdir, ///< Initial orientation.
      enum axes finaldir ///< Final Orientation.
      )
{
  /** \section Directions
   * This only rotates between the three positive directions, XP,YP,ZP, 
   * so it never introduces sign changes to the elements, just reordering
   * of the vector components.
   * */
  if (initdir==finaldir) return;
  pion_flt v[eq_nvar];
  for (int i=0;i<eq_nvar;i++) v[i] = vec[i];
  int offset = (static_cast<int>(finaldir-initdir+3))%3;
  if (offset ==1) {
    v[eqVX] = vec[eqVY];
    v[eqVY] = vec[eqVZ];
    v[eqVZ] = vec[eqVX];
    v[eqBX] = vec[eqBY];
    v[eqBY] = vec[eqBZ];
    v[eqBZ] = vec[eqBX];    
  }
  else if (offset==2) {
    v[eqVX] = vec[eqVZ];
    v[eqVY] = vec[eqVX];
    v[eqVZ] = vec[eqVY];
    v[eqBX] = vec[eqBZ];
    v[eqBY] = vec[eqBX];
    v[eqBZ] = vec[eqBY];
  }
  else rep.error("rotate function broken.",offset);
  for (int i=0;i<eq_nvar;i++) vec[i] = v[i];
}


// ##################################################################
// ##################################################################


void eqns_mhd_ideal::rotateXY(
      pion_flt *v,
      double theta
      )
{
  double ct=cos(theta); double st=sin(theta);

  pion_flt vx = v[eqVX]*ct - v[eqVY]*st;
  pion_flt vy = v[eqVX]*st + v[eqVY]*ct;
  v[eqVX] = vx; v[eqVY]=vy;

  vx = v[eqBX]*ct - v[eqBY]*st;
  vy = v[eqBX]*st + v[eqBY]*ct;
  v[eqBX] = vx; v[eqBY]=vy;
}


// ##################################################################
// ##################################################################


///  Returns Internal Energy (per unit mass, so 'Temperature'), given primitive variable vector. 
double eqns_mhd_ideal::eint(
      const pion_flt *p, ///< Primitive State Vector.
      const double g   ///< gas EOS gamma.
      )
{
  return p[eqPG]/(g-1.)/p[eqRO];
}


// ##################################################################
// ##################################################################

/// Returns Total Energy (per unit volume), given primitive variable vector. 
double eqns_mhd_ideal::Etot(
      const pion_flt *p, ///< State Vector.
      const double g   ///< gas EOS gamma.
      ) 
{
  return( (p[eqRO]*(p[eqVX]*p[eqVX]+p[eqVY]*p[eqVY]+p[eqVZ]*p[eqVZ])/2.)
	  +(p[eqPG]/(g-1.))
	  +((p[eqBX]*p[eqBX] +p[eqBY]*p[eqBY] +p[eqBZ]*p[eqBZ])/2.));
}


// ##################################################################
// ##################################################################


/// Returns Total Pressure (per unit Volume), given primitive variable vector. 
double eqns_mhd_ideal::Ptot(
      const pion_flt *p, ///< Primitive State Vector.
      const double    ///< gas EOS gamma.
      )
{
return( p[eqPG]+ 0.5*(p[eqBX]*p[eqBX] +p[eqBY]*p[eqBY] +p[eqBZ]*p[eqBZ]) );
}


// ##################################################################
// ##################################################################


/// Given a pressure ratio and initial density, calculate adiabatic final density.
double eqns_mhd_ideal::AdiabaticRho(
      const double pr, ///< New to Old pressure ratio
      const double ri, ///< Old Density
      const double g ///< gas EOS gamma.
      )
{
  return(ri*exp(log(pr)/g));
}


// ##################################################################
// ##################################################################

void eqns_mhd_ideal::SetAvgState(
      const pion_flt *state, ///< Mean Primitive var. state vector
      const double g ///< Gas constant gamma.
      )
{
  //
  // Set typical values to be used for reference:
  // (only really care about refB, refvel, refvec[PG,RO])
  //
  eq_refvec[eqRO] = state[eqRO];
  eq_refvec[eqPG] = state[eqPG];
  eq_refvec[eqVX] = state[eqVX];
  eq_refvec[eqVY] = state[eqVY];
  eq_refvec[eqVZ] = state[eqVZ];
  eq_refvec[eqBX] = state[eqBX];
  eq_refvec[eqBY] = state[eqBY];
  eq_refvec[eqBZ] = state[eqBZ];
  
  double angle = 0.0, refvel=0.0, refB=0.0;
  angle = eq_refvec[eqBY]*eq_refvec[eqBY]+eq_refvec[eqBX]*eq_refvec[eqBX];
  if (angle > 10.*MACHINEACCURACY) {
    angle = M_PI/2. -asin(eq_refvec[eqBY]/sqrt(angle));
    if(eq_refvec[eqBX]<0) angle = -angle;
    rotateXY(eq_refvec,angle);
    refvel = cfast(eq_refvec, eq_gamma); // Fast speed
    rotateXY(eq_refvec,-angle);
  }
  else refvel = maxspeed(eq_refvec, eq_gamma);
  
  refB = sqrt(eq_refvec[eqBX]*eq_refvec[eqBX]+
	      eq_refvec[eqBY]*eq_refvec[eqBY]+
	      eq_refvec[eqBZ]*eq_refvec[eqBZ]);
  //
  // Now reset reference vector velocities and B-fields to be refvel,
  // and refB
  //
  eq_refvec[eqVX] = eq_refvec[eqVY] = eq_refvec[eqVZ] = 0.1*refvel;
  eq_refvec[eqBX] = eq_refvec[eqBY] = eq_refvec[eqBZ] = refB;
  
  rep.printVec("eq_refvec",eq_refvec,eq_nvar);
  return;
}



// ##################################################################
// ##################################################################



// ******************************************************************
// eqns_mhd_mixedGLM class, for the Dedner-GLM divergence cleaning method.
// ******************************************************************


// ##################################################################
// ##################################################################

eqns_mhd_mixedGLM::eqns_mhd_mixedGLM(int nv)
  : eqns_base(nv), eqns_mhd_ideal(nv)
{
  //  cout <<"eqns_mhd_mixedGLM constructor!\n";
  eqSI = SI;
  eqPSI= PSI;
}


// ##################################################################
// ##################################################################

eqns_mhd_mixedGLM::~eqns_mhd_mixedGLM()
{
}


// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::GLMsetPsiSpeed(
      const double cfl,
      const double delx,
      const double delt
      )
{
  ///
  /// \section chyp Hyperbolic Wave Speed
  /// See the Code Algorithms page \ref algorithms for details of my
  /// implementation.
  ///

  //  cout <<"calculating GLM_chyp!\n";
  GLM_chyp = cfl*delx/delt; /// hyperbolic wavespeed is equal to max. allowed value for given CFL no.
  GLM_cr = 4.0*delx; // This works well for general use.
  //GLM_cr *=20.5; // This is for when getting negative pressure near outflow boundaries.
  //GLM_cr = 0.3*delx;
  
  //
  // The following method is from Dedner's thesis (eq.8.22, p.121), and is 
  // larger than using 4dx for the advection problem (i get 14.051 for divBpeak
  // problem, as oppose to 4 for what I have above).
  //
  //  GLM_cr=0.;
  //  for (int i=0;i<SimPM.ndim;i++)
  //    GLM_cr += 4./(SimPM.Xmax[i]-SimPM.Xmin[i])/(SimPM.Xmax[i]-SimPM.Xmin[i]);
  //  GLM_cr = 1./M_PI/sqrt(GLM_cr);


  //  cout <<"GLM_cr="<<GLM_cr<<"\n";
  //  cout <<"GLM_chyp = "<<GLM_chyp<<" and cfl*dx/dt="<<cfl*delx/delt<<"\n";
  return;
}
			    

// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::PtoU(
      const pion_flt *P,    ///< pointer to Primitive variables.
      pion_flt *U,          ///< pointer to conserved variables.
      const double gamma  ///< Gas constant gamma.
      )
{
//  cout <<"glm ptou\n";
  U[eqPSI] = P[eqSI];
  eqns_mhd_ideal::PtoU(P,U,gamma);
  return;
}


// ##################################################################
// ##################################################################

int eqns_mhd_mixedGLM::UtoP(
      const pion_flt *U, ///< pointer to conserved variables.
      pion_flt *P,       ///< pointer to Primitive variables.
      const double g   ///< Gas constant gamma.
      )
{
  //  cout <<"glm utop\n";
  P[eqSI]=U[eqPSI];
  return(eqns_mhd_ideal::UtoP(U,P,g));
}


// ##################################################################
// ##################################################################

void eqns_mhd_mixedGLM::GLMsource(
      pion_flt *psivar, ///< Primitive Psi variable.
      const double delt ///< timestep
      )
{
  //cout <<"\tglmsource: exp factor="<<-delt*GLM_chyp/GLM_cr<<"\n";
  *psivar *= exp(-delt*GLM_chyp/GLM_cr);
  return;
}


// ##################################################################
// ##################################################################

