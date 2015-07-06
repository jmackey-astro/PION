///
/// \file eqns_hydro_adiabatic.h
/// \author Jonathan Mackey
/// Declaration of the adiabatic hydrodynamics equations class
///
/// History: 2009-10-20 Started on the file (moved from old equations.h)
///
/// - 2010.12.21 JM: Added Euler equations class which uses internal
///   energy as a conserved variable, and not the total energy.
///   Tidied up comments on Euler Equations class.
///
/// - 2010.12.23 JM: Added SetAvgState() function to eqns_Euler class.
///
/// - 2010.12.28 JM: moved internal energy Euler eqns to own file.
///

#ifndef EQNS_HYDRO_ADIABATIC_H
#define EQNS_HYDRO_ADIABATIC_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "eqns_base.h"
///
/// Class describing the Euler Equations of Inviscid Hydrodynamics.
/// 
/// This class uses [density, pressure, velocity] as primitive variables
/// and [density, *TOTAL ENERGY*, momentum] as conserved variables.
///
class eqns_Euler :  virtual public eqns_base {
 public:
  /// Constructor, sets state vector length, and initial flux
  /// direction to XX.
  eqns_Euler(int);
  ~eqns_Euler(); ///< Destructor.
  
  ///
  ///  Converts from primitive to conserved variables. 
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
  /// Converts from primitive and conserved variables to corresponding
  /// flux.
  ///
  /// This assumes that the direction has been set correctly.
  ///
  virtual void PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);

  ///
  /// Converts from conserved variables to flux. 
  ///
  virtual void UtoFlux(const double *, ///< Pointer to conserved variables state vector.
		       double *,       ///< Pointer to flux variable state vector.
		       const double    ///< Gas constant gamma.
		       );

  ///
  /// Calculate hydrodynamic sound speed.
  /// \retval ch(>0) Success
  /// \retval -1     Failure
  ///
  virtual double chydro(const double *, ///< Pointer to primitive variables.
			const double    ///< Gas constant gamma.
			);
  
  ///
  /// Returns the fastest wavespeed for the relevant equations.
  /// For eqns_Euler it returns the hydro speed.
  ///
  virtual inline double maxspeed(const double *p, ///< Pointer to primitive variables.
				 const double g   ///< Gas constant gamma.
				 ) {return(chydro(p,g));}

  ///
  ///  Calculates u* given p* and a pre-wave state.
  ///
  /// Hydro wave takes in a pre-wave state, a postwave pressure (pp), 
  ///and returns a postwave velocity, u*, for a pure hydrodynamic wave,
  ///rarefaction or shock.  The int L/R tells it if it's a left or right
  ///moving wave, so it knows how to calculate the velocity correctly.
  ///\retval 0 success
  ///\retval 1 failure
  ///
  int HydroWave(int,            ///< whether it is a left (XN) or right (XP) moving wave.
		const double,   ///< The post-wave pressure, p*
		const double *, ///< pointer to the pre-wave Primitive state vector.
		double *,       ///< pointer to post-wave velocity variable.
		const double    ///< Gas constant gamma.
		);
  
  ///
  ///  Calculates u* and rho* given p* and a pre-wave state.
  ///
  ///HydroWaveFull takes in a prewave state and a postwave pressure (pp), and 
  ///then calls HydroWave() to get the postwave velocity, and then uses this to
  ///calculate the postwave density.  The density calculation is separate 
  ///because in the code we only want that at the end, whereas we use the 
  ///velocity calculation to match waves and find the correct starred state.
  ///This makes things more efficient.
  ///\retval 0 success
  ///\retval 1 failure
  ///
  int HydroWaveFull(int,            ///< whether it is a left (XN) or right (XP) moving wave.
		    const double,   ///< The post-wave pressure, p*
		    const double *, ///< pointer to the pre-wave Primitive state vector.
		    double *,       ///< pointer to post-wave velocity variable.
		    double *,       ///< pointer to post-wave density variable.
		    const double    ///< Gas constant gamma.
		    );
  
  ///
  /// Returns Internal Energy (per unit mass, so 'Temperature'), given
  /// primitive variable vector.
  ///
  virtual double eint(const double *, ///< Primitive State Vector.
		      const double    ///< gas EOS gamma.
		      );

  ///
  /// Returns Enthalpy (per unit mass), given primitive variable
  /// vector.
  ///
  virtual double Enthalpy(const double *, ///< Primitive State Vector.
			  const double ///< gas EOS gamma.
			  );

  ///
  /// Returns Total Energy (per unit volume), given primitive variable
  /// vector.
  ///
  virtual double Etot(const double *, ///< Primitive State Vector.
		      const double    ///< gas EOS gamma.
		      );

  ///
  /// Returns Total Pressure, given primitive variable vector.
  ///
  virtual double Ptot(const double *, ///< Primitive State Vector.
		      const double    ///< gas EOS gamma.
		      );

  ///
  /// Given a pressure ratio and initial density, calculate adiabatic
  /// final density.
  ///
  virtual double AdiabaticRho(const double, ///< New to Old pressure ratio
			      const double, ///< Old Density
			      const double  ///< gas EOS gamma.
			      );

  ///
  /// Set Values for mean velocity, pressure, density.
  ///
  virtual void SetAvgState(const double *, ///< Mean Primitive var. state vector
			   const double ///< Gas constant gamma.
			   );


 protected:
};


#endif // EQNS_HYDRO_ADIABATIC_H
