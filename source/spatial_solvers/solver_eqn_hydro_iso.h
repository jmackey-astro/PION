///
/// \file solver_eqn_hydro_iso.h
/// \author Jonathan Mackey
///
/// Solver for the isothermal Euler Equations.  Calculates flux via either Lax-Friedrichs
/// or Riemann solver (linear and/or exact).  Adds viscosity if asked for, and tracks flux
/// of N passive tracers.
///
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
///  - 2010.11.15 JM:
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
///

#ifdef ISOTHERMAL_SOLVERS_ENABLED


#ifndef SOLVER_EQN_HYDRO_ISO_H
#define SOLVER_EQN_HYDRO_ISO_H

#include "solver_eqn_base.h"
#include "../flux_calc/flux_hydro_isothermal.h"

///
/// The main solver the code uses for integrating the Euler Equations (isothermal).
/// This solver can get 2 flux solutions so far:
/// 0=Lax-Friedrichs, 1=Linear Riemann
///
class FV_solver_Hydro_iso : virtual public FV_solver_base,
  virtual public flux_solver_hydro_iso,
  virtual public VectorOps_Cart
{
 public:
  FV_solver_Hydro_iso(const int, ///< number of variables in state vector.
		      const int, ///< number of space dimensions in grid.
		      const double, ///< CFL number
		      const double, ///< dx, cell size.
		      const double, ///< gas eos gamma.
		      double *,     ///< State vector of mean values for simulation.
		      const double, ///< Artificial Viscosity Parameter etav.
		      const int     ///< Number of tracer variables.
		      );
  ~FV_solver_Hydro_iso();
  
   ///
   /// Adds the contribution from flux in the current direction to dU.
   ///
   virtual int dU_Cell(cell *, ///< Current cell.
		       const axes, ///< Which axis we are looking along.
		       const double *, ///< Negative direction flux.
		       const double *, ///< Positive direction flux.
		       const double *, ///< slope vector for cell c.
		       const int,      ///< spatial order of accuracy.
		       const double, ///< cell length dx.
		       const double  ///< cell TimeStep, dt.
		       );
   ///
   /// General Finite volume scheme for updating a cell's
   /// primitive state vector, for homogeneous equations.
   ///
   virtual int CellAdvanceTime(const double *, ///< Initial Primitive State Vector.
			       double *, ///< Update vector dU
			       double *, ///< Final Primitive state vector (can be same as initial vec.).
			       double *,  ///< Tracks change of energy if I have to correct for negative pressure
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell timestep dt.
			       );
   ///
   /// Given a cell, calculate the hydrodynamic timestep.
   ///
   virtual double CellTimeStep(const cell *, ///< pointer to cell
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell size dx.
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
};

#endif // SOLVER_EQN_HYDRO_ISO_H

#endif // ISOTHERMAL_SOLVERS_ENABLED
