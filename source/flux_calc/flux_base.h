/// \file flux_base.h
/// \author Jonathan Mackey
/// 
/// Contains class declarations for the basic flux-solver class.
/// This is an interface class, and also defines some basic
/// functionality, such as how to deal with tracer variables tacked
/// onto the end of the state vectors.
///
/// Modifications:\n
/// - 2010.12.23 JM: Moved from eqns_base.h into its own file.
///   Added tracer_flux() function.
///
/// - 2010.12.27 JM: Modified viscosity functions for new structure.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.


#ifndef FLUX_BASE_H
#define FLUX_BASE_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "equations/eqns_base.h"
#include "grid/cell_interface.h" // to get the 'cell' class.

///
/// Abstract base class for a Flux Solver Class.  This class should be
/// grid-ignorant and geometry-ignorant, so just for calculating the 
/// flux across an interface given a left and right state.
///
class flux_solver_base : virtual public eqns_base 
{
 public:
  flux_solver_base(const int, ///< number of variables in state vector
		   const double,  ///< coefficient of artificial viscosity.
		   const int     ///< Number of tracer variables.
		   );
  ~flux_solver_base();
  ///
  /// Calculates Flux based on a left and right state vector (primitive).
  ///
  virtual int inviscid_flux(const cell *,  ///< Left state cell pointer
			    const cell *,  ///< Right state cell pointer
			    const double *,///< Left Primitive state vector.
			    const double *,///< Right Primitive state vector.
			    double *,      ///< Resultant Flux state vector.
			    double *,      ///< State vector at interface.
			    const int, 
	   ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
			    const double    ///< Gas constant gamma.
			    ) =0;
protected:
  ///
  /// Falle, Komissarov & Joarder (1998,MNRAS,297,265) Artificial
  /// Viscosity Calculation (one-dimensional).
  ///
  virtual int AVFalle(const double *, ///< Left Primitive state vector.
		      const double *, ///< Right Primitive state vector.
		      const double *, ///< Resolved (P*) state vector.
		      double *,     ///< Pointer to associated Flux Vector.
		      const double, ///< Artificial Viscosity parameter, etav.
		      const double  ///< gamma
		      )=0;

#ifdef LAPIDUS_VISCOSITY_ENABLED
  ///
  ///Lapidus Artificial Viscosity Calculation. 
  /// This calculates div(v) at the interface and subtracts a fraction of it from
  /// the flux \f$ F = F + \eta \mbox{div}(v)(U_R-U_L) \f$.
  ///
  int AVLapidus(const cell *, ///< Left state cell pointer
		const cell *, ///< Right state cell pointer
		double *,       ///< Pointer to associated Flux Vector.
		const double, ///< Artificial Viscosity parameter, etav.
		const double  ///< gamma
		);
#endif // LAPIDUS_VISCOSITY_ENABLED

  ///
  /// This calculates the Lax-Friedrichs flux across an interface, but the 
  /// implementation will have to be in the general solver classes because LF flux
  /// requires dx, dt, and ndim, which the equations/riemann solvers don't know 
  /// about and don't care about.  I'm putting in a dummy implementation here so that
  /// I can set up a flux_solver class if I want to go without it.
  ///
  virtual int get_LaxFriedrichs_flux(const double *, ///< Left  Prim. var. state vector.
				     const double *, ///< Right Prim. var. state vector.
				     double *,       ///< Resulting Flux vector.
				     const double    ///< gamma
				     );

  const double FS_etav; ///< coefficient of (artificial) viscosity for velocity field.
  const double FS_etaB; ///< coefficient of (artificial) viscosity for magnetic field.
  double HC_etamax; ///< max. value of eta used for H-correction.

  /// Number of passive tracer variables.
  const int FS_ntr;
  /// Pointer to array of indices of tracer variables in the state vector.
  int *eqTR;

  ///
  /// Calculate tracer flux based on a flux vector and left and right
  /// states.  If mass flux is positive, then tracer flux is advected
  /// from the left state; if negative then from the right state.
  /// Flux must be already calculated for the physical variables!
  ///
  void set_interface_tracer_flux(const double *, ///< input left state.
				 const double *, ///< input right state
				 double *        ///< Flux vector.
				 );


};

#endif //FLUX_BASE_H

