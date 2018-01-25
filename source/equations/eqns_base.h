/// \file eqns_base.h
/// \brief Contains Class declarations for equations classes.
/// 
/// \author Jonathan Mackey
/// 
/// Modifications:
///  - 2007-10-16 Moved equations classes from global.h into this file.
///  - 2009-10-20 Changed name to eqns_base.h and moved euler and mhd equations to other files.
///  - 2009-12-19 JM: Added MaxSpeed() function.
///  - 2010.09.30 JM: Worked on Lapidus AV: added Cl,Cr pointers to flux 
///     functions for base solver class.
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.
///  - 2010.11.19 JM: Added H-correction eta variable to flux solver.
/// - 2010.12.23 JM: Removed riemann_base.  Moved SetAvgState to
///   eqns_base class.  Renamed RS_meanvec to eq_meanvec.
///   Moved flux_solver_base class to flux_base.h
/// - 2011.01.03 JM: Added eq_posdir and eq_negdir direction variables.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef EQNS_BASE_H
#define EQNS_BASE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "constants.h"  // for axes enum in function arguments.

/// Abstract Base Class for equations, from which others are derived. 
///
/// The public functions declared here are basically the only ones accessible
/// from the main Grid Code, because it declares a pointer to the simplest
/// solver, BasicFVSolver, which inherits from : virtual public eqns_base.
///
/// This class doesn't actually solve any equations itself, but any (at least
/// hyperbolic) set of equations should be deriveable from this class.
///
class eqns_base {
 public:
  eqns_base(const int ///< Number of Variables in State Vector
	    );
  virtual ~eqns_base();

  ///  Converts from primitive to conserved variables. 
  virtual void PtoU(
      const pion_flt *, ///< pointer to Primitive variables.
      pion_flt *,       ///< pointer to conserved variables.
      const double    ///< Gas constant gamma.
      ) =0;

  ///
  ///  convert from conserved to primitive variables.  This needs a return
  /// value in case of negative pressure, so we can tell dU_Cell() that the 
  /// update went badly.
  ///
  virtual int UtoP(
      class SimParams &, ///< pointer to simulation parameters
      const pion_flt *, ///< pointer to conserved variables.
      pion_flt *, ///< pointer to Primitive variables.
      const double    ///< Gas constant gamma.
      ) =0;

  ///  Converts from primitive and conserved variables to corresponding flux.
  ///This assumes that the direction has been set correctly.
  ///
  virtual void PUtoFlux(
      const pion_flt *, ///< pointer to Primitive variables.
      const pion_flt *, ///< pointer to conserved variables.
      pion_flt *  ///< Pointer to flux variable.
      ) =0;

  ///  convert direct from primitive variables to flux.
  ///Creates conserved variable array as an intermediate step, 
  ///and then calls PUtoFlux(). 
  virtual void PtoFlux(
      const pion_flt *, ///< pointer to Primitive variables.
      pion_flt *,       ///< Pointer to flux variable.
      const double    ///< Gas constant gamma.
      );

  ///  Converts from conserved variables to flux. 
  virtual void UtoFlux(
      const pion_flt *, ///< Pointer to conserved variables state vector.
      pion_flt *,       ///< Pointer to flux variable state vector.
      const double   ///< Gas constant gamma.
      ) =0;

  ///
  /// Returns Internal Energy (per unit mass, so 'Temperature'), given
  /// a primitive variable vector.
  ///
  virtual double eint(
      const pion_flt *, ///< Primitive State Vector.
      const double ///< gas EOS gamma.
      ) =0;

  ///
  /// Returns Total Energy (per unit volume), given primitive variable
  /// vector.
  ///
  virtual double Etot(
      const pion_flt *, ///< Primitive State Vector.
      const double ///< gas EOS gamma.
      ) =0;

  ///
  /// Returns Enthalpy (per unit mass), given primitive variable
  /// vector.
  ///
  virtual double Enthalpy(
      const pion_flt *, ///< State Vector.
      const double ///< gas EOS gamma.
      ) =0;   

  ///
  /// Returns Total Pressure (per unit Volume), given primitive
  /// variable vector.
  ///
  virtual double Ptot(
      const pion_flt *, ///< Primitive State Vector.
      const double ///< gas EOS gamma.
      ) =0;

  ///
  /// Given a pressure ratio and initial density, calculate adiabatic
  /// final density.
  ///
  virtual double AdiabaticRho(
      const double, ///< New to Old pressure ratio
      const double, ///< Old Density
      const double ///< gas EOS gamma.
      ) =0;

  ///
  ///  Returns the fastest wavespeed for the relevant equations.
  ///  For eqns_Euler it returns the hydro speed, for MHD the fast speed.
  ///
  virtual double maxspeed(
        const pion_flt *p, ///< Pointer to primitive variables.
        const double g   ///< Gas constant gamma.
        ) =0;

  ///  Sets the direction I'm looking, to decide which flux function to call.
  virtual void SetDirection(
    const enum axes ///< Direction, XX,YY or ZZ.
    );

  ///
  ///  Returns the currently set direction of equations. 
  ///
  virtual enum axes GetDirection();

  ///
  ///  Rearranges a vector for a rotation by pi/2 
  ///
  virtual void rotate(
      pion_flt *, ///< State vector
      enum axes, ///< Initial orientation.
      enum axes  ///< Final Orientation.
      );

  ///
  /// Rotates a state vector through theta radians.
  /// Note this is equivalent to rotating the axes through -theta radians.
  /// The rotation is happening in the x-y plane, so if it is called in a 
  /// 3d setting it is equivalent to a rotation about the z-axis.  This
  /// version just rotates the fluid velocity components.
  ///
  virtual void rotateXY(
      pion_flt *, ///< State vector
      double    ///< angle to rotate through in 2d plane.
      );

  ///
  /// Set Values for mean velocity, pressure, density, B-field, etc.
  /// (velocity is set by sound speed)
  ///
  virtual void SetAvgState(
      const pion_flt *,    ///< Mean Primitive var. state vector
      const double       ///< Gas constant gamma.
      ) =0;

  protected:

  int eq_nvar;      ///< number of elements in state vector (total, including tracers, etc.)
  enum axes eq_dir; ///< Which axis we are looking along XX,YY,ZZ.
  enum direction eq_posdir; ///< positive direction along current axis.
  enum direction eq_negdir; ///< negative direction along current axis.
  enum primitive eqRO,   eqPG; ///< Here so tracer variables can be advected by generic functions.
  enum conserved eqRHO, eqERG; ///< Here so tracer variables can be advected by generic functions.
  enum primitive eqVX, eqVY, eqVZ, eqBX, eqBY, eqBZ;
  enum conserved eqMMX, eqMMY, eqMMZ, eqBBX, eqBBY, eqBBZ;
  double eq_gamma;  ///< Gas equation of state gamma.

  ///
  /// state vector of typical values in the simulation (for testing
  /// if a variable value is small or large).
  ///
  pion_flt *eq_refvec; 
   
};

#endif //EQNS_BASE_H
