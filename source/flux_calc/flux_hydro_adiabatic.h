///
/// \file flux_hydro_adiabatic.h
/// \author Jonathan Mackey
///
/// Declaration of the adiabatic hydrodynamics flux solver class.
///
/// In the heirarchy of classes, the equations are at the bottom,
/// the Riemann solvers are next up, then the flux solver defined
/// here, and above this is the full solver for the Euler equations,
/// which knows about the grid, the timestep etc.
/// 
/// The flux solver is where tracer variables are introduced (Riemann
/// solvers don't know about them, nor do equations classes).  So we
/// redefine PtoU(), PUtoFlux(), etc. to also convert tracer values
/// after converting the hydrodynamic variables.
///
/// History:
///  - 2009-10-20 Started on the file
///
///  - 2010-07-31 JM: Added Primitive variable linear solver with Roe average.
///
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.
///
/// - 2010.11.21 JM: Added symmetric version of the Roe conserved
///   variables flux solver, called Roe_flux_solver_symmetric().  It
///   is MUCH BETTER than the one-sided one, which I have renamed to
///   Roe_flux_solver_onesided().  For an axi-symmetric blast wave
///   with the H-correction, the one-sided solver developed spikes at
///   90 degree intervals, but the symmetric solver is perfectly
///   clean.  SO, DON'T USE ONE-SIDED ROE FLUX SOLVER ANYMORE!
///
/// - 2010.12.22 JM: Moved Roe PV and CV solvers to Riemann_solvers/
///   directory.  Got rid of the ROE_CV_MEANP ifdef (it is default).
///
/// - 2010.12.23 JM: Moved UtoP() etc. from solver to flux-solver.
///   Now all tracer stuff is dealt with here.
///
/// - 2010.12.27 JM: Enclosed Lapidus AV in ifdef.

#ifndef FLUX_HYDRO_ADIABATIC_H
#define FLUX_HYDRO_ADIABATIC_H

#include "flux_base.h"
#include "../Riemann_solvers/riemann.h"
#include "../Riemann_solvers/Riemann_FVS_hydro.h"
#include "../Riemann_solvers/Roe_Hydro_PrimitiveVar_solver.h"
#include "../Riemann_solvers/Roe_Hydro_ConservedVar_solver.h"

///
/// Flux Solver Class for Adiabatic Hydrodynamics (Euler Equations)
///
/// This uses the adiabatic equations, and calculates the flux across an interface.
///
/// \author Jonathan Mackey
/// Written 2009-10-20
///
class flux_solver_hydro_adi
: virtual public flux_solver_base,
  virtual public riemann_Euler,
  virtual public Riemann_FVS_Euler,
  virtual public Riemann_Roe_Hydro_PV,
  virtual public Riemann_Roe_Hydro_CV
{
public:

  /// Constructor: Allocates memory for arrays, and sets a state vector of typical values for
  /// primitive variables (can be overwritten later).
  flux_solver_hydro_adi(const int,      ///< Number of variables in state vector
			const double *, ///< state vector which is 'typical' in the problem being solved. 
			const double,   ///< coefficient of artificial viscosity
			const double,    ///< gamma (EOS).
			const int     ///< Number of tracer variables.
			);

  /// Deletes any arrays.
  ~flux_solver_hydro_adi();

  /// Calculates Flux based on a left and right state vector (primitive).
  int inviscid_flux(const cell *, ///< Left state cell pointer
		    const cell *, ///< Right state cell pointer
		    const double *, ///< Left Primitive state vector.
		    const double *, ///< Right Primitive state vector.
		    double *,       ///< Resultant Flux state vector.
		    double *,      ///< Resultant Pstar state vector.
		    const int,      ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
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
  /// This calls the original version and then adds conversion of tracer variables.
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
  /// This calls the original version and then adds conversion of tracer variables.
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
  /// This calls the original version and then adds conversion of tracer variables.
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


#endif //FLUX_HYDRO_ADIABATIC_H
