///
/// \file eqns_hydro_adi_Eint.h
///
/// \author Jonathan Mackey
///
/// Declaration of the adiabatic hydrodynamics equations class, which 
/// uses internal energy as a 'conserved' variable instead of the 
/// total energy.
///
/// Modifications:\n
/// - 2010.12.28 JM: Created (moved from eqns_hydro_adiabatic.h)
/// - 2010.12.30 JM: Added eqEINT local variable.

#ifndef EQNS_HYDRO_ADI_EINT_H
#define EQNS_HYDRO_ADI_EINT_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef INCLUDE_EINT_ADI_HYDRO

#include "eqns_hydro_adiabatic.h"
#include "eqns_base.h"

///
/// Class describing the Euler Equations of Inviscid Hydrodynamics.
/// 
/// This class uses [density, pressure, velocity] as primitive variables
/// and [density, *INTERNAL ENERGY*, momentum] as conserved variables.
///
class eqns_Euler_Eint :  virtual public eqns_Euler {
 public:
  eqns_Euler_Eint(int); ///< Constructor, sets state vector length, and initial flux direction to XX.
  ~eqns_Euler_Eint(); ///< Destructor.
  
  ///
  /// Converts from primitive to conserved variables. 
  ///  
  virtual void PtoU(const double *, ///< pointer to Primitive variables.
		    double *,       ///< pointer to conserved variables.
		    const double    ///< Gas constant gamma.
		    );
  
  ///
  ///  convert from conserved to primitive variables. 
  ///
  virtual int UtoP(const double *, ///< pointer to conserved variables.
		    double *, ///< pointer to Primitive variables.
		    const double    ///< Gas constant gamma.
		    );
  ///
  /// Converts from primitive and conserved variables to corresponding flux.
  ///
  /// This assumes that the direction has been set correctly!
  /// Note that the energy flux is just vx*Eint, and P.div(V) is 
  /// added to dEint/dt in post-processing.
  ///
  virtual void PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);

  ///
  ///  Converts from conserved variables to flux. 
  ///
  virtual void UtoFlux(const double *, ///< Pointer to conserved variables state vector.
		       double *,       ///< Pointer to flux variable state vector.
		       const double    ///< Gas constant gamma.
		       );
  
  protected:
#ifdef EINT_ETOT_PARALLEL
  int eqEINT;
#endif // EINT_ETOT_PARALLEL
};

#endif // if INCLUDE_EINT_ADI_HYDRO


#endif // EQNS_HYDRO_ADI_EINT_H
