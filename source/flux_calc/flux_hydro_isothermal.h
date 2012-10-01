///
/// \file flux_hydro_isothermal.h
/// \author Jonathan Mackey
/// Declaration of the isothermal hydrodynamics flux solver class
///
/// History: 
/// - 2009-10-21 Moved flux solver from eqns_hydro_isothermal.h
///
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
///

#ifdef ISOTHERMAL_SOLVERS_ENABLED


#ifndef FLUX_HYDRO_ISOTHERMAL_H
#define FLUX_HYDRO_ISOTHERMAL_H

#include "../equations/eqns_hydro_isothermal.h"
#include "../Riemann_solvers/findroot.h"

///
/// Isothermal Hydrodynamics Riemann Solver Class.
/// This can do a linear solver, and hopefully an exact solver before the day is out...
///
class riemann_hydro_iso
: //virtual public riemann_base,
  virtual public eqns_IsoEuler,
  virtual public findroot
{
  public:
   /// 
   /// Assumes an isothermal equation of state.
   ///
   riemann_hydro_iso(const int,     ///< Length of State Vectors, nvar
		     const double * ///< Mean values of primitive variables on grid [vector, length nvar]
		     );
   ///
   ///   Destructor: deletes dynamically allocated data.
   ///
   ~riemann_hydro_iso();

   ///
   /// Gets the solution to the Riemann Problem.
   /// This will probably become private eventually, since calculate_flux should become
   /// the main interface function...
   ///  
   int riemann_solve(const double *, ///< Left Primitive var. state vector.
		     const double *, ///< Right Primitive var. state vector.
		     double *,       ///< Result Primitive var. state vector.
		     const int,      ///< Solve Type (1=LinearRS,2=ExactRS,3=HybridRS)
		     const double    ///< Gas constant gamma.
		     );
   
   ///
   /// Prints out info on what the solver has done, and how well it's done it.
   ///
   ///void testing();
   
   ///
   /// Set Values for mean velocity, pressure, density.
   ///
   virtual void SetAvgState(const double *, ///< Mean Primitive var. state vector
			    const double    ///< Gas constant gamma.
			    );
  protected:

  private:
   ///
   /// Isothermal Hydro linear Riemann Solver
   ///
   int linear_solver_hydro_iso();

   ///
   /// Isothermal Hydro exact Riemann Solver
   ///
   int two_shock_solver_hydro_iso();

   ///
   /// Re-definition of the public root-finding function.
   ///
   virtual int FR_find_root(double *, /**< pointer to result */
			    const double,  ///< parameter 1
			    const double,  ///< parameter 2
			    const double,  ///< parameter 3
			    const double,  ///< parameter 4
			    const double   ///< parameter 5
			    );
   ///
   /// Re-definition of the funtion to find the root of.
   /// Input an x value and this returns f(x).  The aim of this class is to find the x
   /// value for which f(x)=0, so we are trying to minimize the return value of this function.
   ///
   virtual double FR_root_function(const double  ///< x-value
				   );
   

};



//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//



///
/// Flux Solver Class for Isothermal Hydrodynamics
///
/// This uses the isothermal equations, and calculates the flux across an interface.
///
/// \author Jonathan Mackey
/// Written 2009-10-17
///
class flux_solver_hydro_iso
: virtual public flux_solver_base,
  virtual public riemann_hydro_iso
{
public:

  /// Constructor: Allocates memory for arrays, and sets a state vector of typical values for
  /// primitive variables (can be overwritten later).
  flux_solver_hydro_iso(const int,      ///< Number of variables in state vector
			const double *, ///< state vector which is 'typical' in the problem being solved. 
			const double,   ///< coefficient of artificial viscosity (may or may not be used!)
			const int     ///< Number of tracer variables.
			);

  /// Deletes any arrays.
  ~flux_solver_hydro_iso();

  /// Calculates Flux based on a left and right state vector (primitive).
  int inviscid_flux(const cell *, ///< Left state cell pointer
		    const cell *, ///< Right state cell pointer
		    const double *, ///< Left Primitive state vector.
		    const double *, ///< Right Primitive state vector.
		    double *,       ///< Resultant Flux state vector.
		    double *,       ///< Resultant Pstar state vector.
		    const int,
      ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
		    const double    ///< Gas constant gamma.
		    );



protected:
  ///
  /// Falle, Komissarov & Joarder (1998) Artificial Viscosity
  /// Calculation.  This does not update the energy, of course,
  /// just the momentum.
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
  ///
  int AVLapidus(const cell *, ///< Left state cell pointer
		const cell *, ///< Right state cell pointer
		double *,       ///< Pointer to associated Flux Vector.
		const double,    ///< Artificial Viscosity parameter, etav.
		const double  ///< gamma
		);
#endif // LAPIDUS_VISCOSITY_ENABLED
   
   ///
   /// Roe's approximate flux solver:
   ///
   int Roe_solver_hydro_iso();
};

#endif // FLUX_HYDRO_ISOTHERMAL_H


#endif // ISOTHERMAL_SOLVERS_ENABLED
