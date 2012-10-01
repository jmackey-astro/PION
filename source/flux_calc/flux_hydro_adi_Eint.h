///
/// \file flux_hydro_adi_Eint.h
/// \author Jonathan Mackey
///
/// Declaration of the adiabatic hydrodynamics flux solver class,
/// which integrates the INTERNAL ENERGY and NOT the total energy.
///
/// Modifications:\n
///
/// - 2010.12.28 JM: Started file.
/// - 2010.12.30 JM: Can now integrate both internal and total energy.
///

#ifndef FLUX_HYDRO_ADI_EINT_H
#define FLUX_HYDRO_ADI_EINT_H

#ifdef INCLUDE_EINT_ADI_HYDRO


#include "flux_base.h"
#include "flux_hydro_adiabatic.h"
#include "../Riemann_solvers/riemann.h"
#include "../Riemann_solvers/Roe_Hydro_PrimitiveVar_solver.h"

///
/// Flux Solver Class for Adiabatic Hydrodynamics (Euler Equations), 
/// modified to calculate the internal energy flux and NOT the total
/// energy flux.
///
/// This uses the adiabatic equations, and calculates the flux across
/// an interface.  It can't use Flux-vector-splitting or the Roe
/// conserved variables solvers since they return a flux and not an 
/// interface state.  The Hybrid and Roe-PV solvers return a primitive
/// variabel interface state which I can then use to get the flux in
/// [rho,rho*v,Eint].
///
/// \author Jonathan Mackey
/// Written 2010.12.28
///
class flux_solver_hydro_adi_Eint
: virtual public flux_solver_hydro_adi,
  virtual public riemann_Euler,
  virtual public Riemann_Roe_Hydro_PV,
  virtual public eqns_Euler_Eint
{
public:
  /// Constructor: doesn't do much except call the correct
  /// constructors for the classes it is derivedf from.
  flux_solver_hydro_adi_Eint(const int,
      ///< Number of variables in state vector
			     const double *,
      ///< state vector which is 'typical' in the problem being solved. 
			     const double,
      ///< coefficient of artificial viscosity
			     const double,
      ///< gamma (EOS).
			     const int
      ///< Number of tracer variables.
			);

  /// Doesn't do anything much.
  ~flux_solver_hydro_adi_Eint();

  /// Calculates Flux based on a left and right state vector (primitive).
  int inviscid_flux(const cell *, ///< Left state cell pointer
		    const cell *, ///< Right state cell pointer
		    const double *, ///< Left Primitive state vector.
		    const double *, ///< Right Primitive state vector.
		    double *,       ///< Resultant Flux state vector.
		    double *,      ///< Resultant Pstar state vector.
		    const int,
      ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
		    const double    ///< Gas constant gamma.
		    );

  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The 
  /// conserved variable is the mass density of this.
  ///
  virtual void PtoU(const double *, ///< pointer to Primitive variables.
		    double *,       ///< pointer to conserved variables.
		    const double    ///< Gas constant gamma.
		    );

  ///
  /// This calls the original version and then adds conversion of
  /// tracer variables.
  /// 
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The 
  /// conserved variable is the mass density of this.
  ///
  virtual int UtoP(const double *, ///< pointer to conserved variables.
		   double *, ///< pointer to Primitive variables.
		   const double    ///< Gas constant gamma.
		   );

  ///
  /// This calls the original version and then adds conversion of
  /// tracer variables.
  /// 
  /// The flux of a passive tracer is equal to the mass flux times 
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);

  ///
  /// This calls the original version and then adds conversion of
  /// tracer variables.
  /// 
  /// The flux of a passive tracer is equal to the mass flux times 
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void UtoFlux(const double*, ///< Pointer to conserved variables state vector.
		       double*,       ///< Pointer to flux variable state vector.
		       const double   ///< Gas constant gamma.
		       );

protected:
  ///
  /// Sam Falle's Artificial Viscosity Calculation.
  ///
  int AVFalle(const double *, ///< Left Primitive state vector.
	      const double *, ///< Right Primitive state vector.
	      const double *, ///< Resolved (P*) state vector.
	      double *, ///< Pointer to associated Flux Vector.
	      const double, ///< Artificial Viscosity parameter, etav.
	      const double  ///< gamma
	      );

#ifdef LAPIDUS_VISCOSITY_ENABLED
  ///
  ///Lapidus Artificial Viscosity Calculation. 
  /// This calculates div(v) at the interface and subtracts a fraction of it from
  /// the flux \f$ F = F + \eta \mbox{div}(v)(U_R-U_L) \f$.
  /// THIS DOESN'T WORK! I STOPPED WRITING IT BEFORE I GOT IT WORKING.
  ///
  int AVLapidus(const cell *, ///< Left state cell pointer
		const cell *, ///< Right state cell pointer
		const double *, ///< Left Primitive state vector.
		const double *, ///< Right Primitive state vector.
		double *,       ///< Pointer to associated Flux Vector.
		const double, ///< Artificial Viscosity parameter, etav.
		const double  ///< gamma
		);
#endif // LAPIDUS_VISCOSITY_ENABLED

 private:
};





#endif // if INCLUDE_EINT_ADI_HYDRO


#endif // FLUX_HYDRO_ADI_EINT_H
