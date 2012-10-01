///
/// \file eqns_hydro_isothermal.h
/// \author Jonathan Mackey
/// Declaration of the isothermal hydrodynamics equations class.
///
/// History: 2009-10-20 Started on the file.
///          2009-10-21 Moved flux solver to new file.
///          2009-12-21 JM: Added maxspeed() function b/c it's in the base class.
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
///

#ifdef ISOTHERMAL_SOLVERS_ENABLED

#ifndef EQNS_HYDRO_ISOTHERMAL_H
#define EQNS_HYDRO_ISOTHERMAL_H

#include "eqns_base.h"

/// Equations Class for Isothermal Hydrodynamics
///
/// \author Jonathan Mackey
/// Written 2009-10-17
class eqns_IsoEuler : virtual public eqns_base {
public:
  eqns_IsoEuler(const int ///< Number of variables in state vector.
		     );
  ~eqns_IsoEuler();

   /// Converts from primitive to conserved variables.
   virtual void PtoU(const double *, ///< pointer to Primitive variables.
		     double *,       ///< pointer to conserved variables.
		     const double    ///< Gas constant gamma.
		     );
   /// convert from conserved to primitive variables.
   virtual int UtoP(const double *, ///< pointer to conserved variables.
		    double *, ///< pointer to Primitive variables.
		    const double    ///< Gas constant gamma.
		    );
   /// Converts from primitive and conserved variables to corresponding flux.
   /// This assumes that the direction has been set correctly.
   
   virtual void PUtoFlux(const double *, ///< pointer to Primitive variables.
		const double *, ///< pointer to conserved variables.
		double *  ///< Pointer to flux variable.
		);
   /// convert direct from primitive variables to flux.
   /// Creates conserved variable array as an intermediate step, 
   /// and then calls PUtoFlux().
   virtual void PtoFlux(const double *, ///< pointer to Primitive variables.
		       double *,       ///< Pointer to flux variable.
		       const double    ///< Gas constant gamma.
		       );   
   /// Converts from conserved variables to flux.
   virtual void UtoFlux(const double*, ///< Pointer to conserved variables state vector.
	       double*,       ///< Pointer to flux variable state vector.
	       const double   ///< Gas constant gamma.
	       );
   ///  Returns the fastest wavespeed for the relevant equations.
   ///For eqns_Euler it returns the hydro speed.
   ///
   virtual inline double maxspeed(const double *p, ///< Pointer to primitive variables.
				  const double g   ///< Gas constant gamma.
				  ) {return(p[eqAA]);}
   /// Returns Internal Energy (per unit mass, so 'Temperature'), given primitive variable vector.
   virtual double eint(const double *, ///< Primitive State Vector.
		       const double ///< gas EOS gamma.
		       );
   ///  Returns Enthalpy (per unit mass), given primitive variable vector. 
   virtual double Enthalpy(const double *, ///< State Vector.
			   const double ///< gas EOS gamma.
			   ) {rep.error("ENTHALPY FOR ISOTHERMAL EQNS",456); return 456.789;}
   /// Returns Total Energy (per unit volume), given primitive variable vector.
   virtual double Etot(const double *, ///< State Vector.
		       const double ///< gas EOS gamma.
		       );
   /// Returns Total Pressure (per unit Volume), given primitive variable vector.
   virtual double Ptot(const double *, ///< Primitive State Vector.
		       const double ///< gas EOS gamma.
		       );
   /// Given a pressure ratio and initial density, calculate adiabatic final density.*/
   virtual double AdiabaticRho(const double, ///< New to Old pressure ratio
			       const double, ///< Old Density
			       const double ///< gas EOS gamma.
			       );
protected:
  // INHERITED: enum axes dir; ///< Direction we are looking in.
  // INHERITED: const int nvar; ///< Vector length.
  enum primitive eqAA;
  enum conserved eqAAA;
   
};


#endif // EQNS_HYDRO_ISOTHERMAL_H

#endif // ISOTHERMAL_SOLVERS_ENABLED
