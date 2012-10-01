///
/// \file flux_mhd_adiabatic.h
/// \author Jonathan Mackey
/// Declaration of the adiabatic MHD flux solver class
///
/// History: 2009-10-20 Started on the file
///
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.
///
/// - 2010.12.27 JM: Moved Roe flux solver to own class in Riemann_solvers/
///   Got rid of inherited class flux/left/right/pstar variables.
///


#ifndef FLUX_MHD_ADIABATIC_H
#define FLUX_MHD_ADIABATIC_H

#include "../equations/eqns_mhd_adiabatic.h"
#include "../Riemann_solvers/riemannMHD.h"
#include "../Riemann_solvers/Roe_MHD_ConservedVar_solver.h"
#include "flux_base.h"

///
/// Flux Solver Class for Adiabatic MHD
///
/// This uses the adiabatic equations, and calculates the flux across an interface.
///
/// \author Jonathan Mackey
/// Written 2009-10-20
///

class flux_solver_mhd_ideal_adi :
  virtual public flux_solver_base,
  virtual public riemann_MHD,
  virtual public Riemann_Roe_MHD_CV
{
public:
  ///
  /// Constructor: Allocates memory for arrays, and sets a state
  /// vector of typical values for primitive variables (can be
  /// overwritten later).
  ///
  flux_solver_mhd_ideal_adi(const int,      ///< Number of variables in state vector
			    const double *, ///< state vector which is 'typical' in the problem being solved. 
			    const double,   ///< coefficient of artificial viscosity
			    const double,    ///< gamma (EOS).
			    const int     ///< Number of tracer variables.
			    );

  /// Deletes any arrays.
  ~flux_solver_mhd_ideal_adi();
  
  /// Calculates Flux based on a left and right state vector (primitive).
  virtual int inviscid_flux(const cell *, ///< Left state cell pointer
			    const cell *, ///< Right state cell pointer
			    const double *, ///< Left Primitive state vector.
			    const double *, ///< Right Primitive state vector.
			    double *,       ///< Resultant Flux state vector.
			    double *,      ///< State vector at interface.
			    const int, 
     ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS,4=RoeRS)
			    const double    ///< Gas constant gamma.
			    );

protected:
  ///
  /// Falle, Komissarov & Joarder (1998,MNRAS,297,265) Artificial
  /// Viscosity Calculation (one-dimensional).
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
		double *,       ///< Pointer to associated Flux Vector.
		const double, ///< Artificial Viscosity parameter, etav.
		const double  ///< gamma
		);
#endif // LAPIDUS_VISCOSITY_ENABLED

  /// shut off reporting if we get more than 1000 negative pressures.
  long int negPGct;
  /// shut off reporting if we get more than 1000 negative densities.
  long int negROct;
};



// **********************************************************************************
// flux_solver_mhd_mixedGLM_adi class, for the Dedner-GLM divergence cleaning method.
// **********************************************************************************


class flux_solver_mhd_mixedGLM_adi
: virtual public flux_solver_mhd_ideal_adi,
  virtual public eqns_mhd_mixedGLM
{
  public:
   flux_solver_mhd_mixedGLM_adi(const int,      ///< Number of variables in state vector
				const double *, ///< state vector which is 'typical' in the problem being solved. 
				const double,   ///< coefficient of artificial viscosity
				const double,    ///< gamma (EOS).
				const int     ///< Number of tracer variables.
				);
   ~flux_solver_mhd_mixedGLM_adi();
   
   ///
   /// Calculates Flux based on a left and right state vector (primitive).
   /// This is the same as for ideal MHD except that we use Dedner et al. (2002)'s 
   /// method to calculate the flux in BX and in the extra scalar field Psi.
   ///
   /// Uses Dedner eq.41 for the flux in Bx and Psi:
   /// \f[ \partial_t B_x + \partial_x \psi = 0 \;, \qquad 
   ///     \partial_t \psi + \partial_x (c_h^2 B_x) = 0 \;. \f]
   /// where the source term has been omitted, as it is calculated separately.
   /// 
   /// The GLM method has Bx and Psi decoupled from all the other variables
   /// in the Riemann Problem, so they can be solved separately as a two
   /// variable system (Dedner eq.42)
   /// 
   /// \f[ F(\psi) = c_h^2 B_x^* = c_h^2 \left( \frac{1}{2}(B_x(L)+B_x(R)) - \frac{1}{2c_h}(\psi_R-\psi_L) \right) \f]
   /// \f[ F(B_x) = \psi_* = \frac{1}{2}(\psi_L+\psi_R) - \frac{c_h}{2}(B_x(R)-B_X(L)) \f]
   /// 
   virtual int inviscid_flux(const cell *, ///< Left state cell pointer
			     const cell *, ///< Right state cell pointer
			     const double *, ///< Left Primitive state vector.
			     const double *, ///< Right Primitive state vector.
			     double *,       ///< Resultant Flux state vector.
			     double *,      ///< State vector at interface.
			     const int,
      ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS,4=RoeRS)
			     const double    ///< Gas constant gamma.
			     );
};

#endif //FLUX_MHD_ADIABATIC_H
