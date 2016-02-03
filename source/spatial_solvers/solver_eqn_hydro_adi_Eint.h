///
/// \file solver_eqn_hydro_adi_Eint.h
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations which integrates the
/// INTERNAL ENERGY and NOT the total energy.  Calculates flux via
/// either Lax-Friedrichs or Riemann solver (linear and/or exact).
/// Also tracks flux of N passive tracers.
///
/// Modifications:\n
/// - 2010.12.28 JM: Started.
/// - 2010.12.30 JM: Can now integrate both internal and total energy.
/// - 2011.01.03 JM: Added PostProcess_dU() for artificial viscosity.
///   Added functions to calculate the Q-viscosity, including a 
///   redefinition of preprocess_data()
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
///

#ifndef SOLVER_EQN_HYDRO_ADI_EINT_H
#define SOLVER_EQN_HYDRO_ADI_EINT_H

#ifdef INCLUDE_EINT_ADI_HYDRO


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "equations/eqns_hydro_adi_Eint.h"
#include "flux_calc/flux_hydro_adi_Eint.h"
#include "solver_eqn_hydro_adi.h"


///
/// The main solver the code uses for integrating the Euler Equations (adiabatic),
/// but integrating the INTERNAL ENERGY and NOT the total energy.
///
/// This can make a few types of flux calculation:
/// 0= Lax-Friedrichs flux: Very diffusive, 1st order.
/// 1= Linear Riemann:      Uses arithmetic mean of left and right states.  Not great.
/// 2= Exact Riemann:       Very slow, but the best solver.
/// 3= Hybrid Riemann (Linear/Exact): uses linear solver when possible.
/// 5= Roe primitive variable solver: linear solver using the Roe-average state.
///
class FV_solver_Hydro_Euler_Eint :
 virtual public FV_solver_Hydro_Euler,
 virtual public flux_solver_hydro_adi_Eint,
 virtual public VectorOps_Cart
{
 public:
  FV_solver_Hydro_Euler_Eint(const int,
  ///< number of variables in state vector.
			     const int,
  ///< number of space dimensions in grid.
			     const double,
  ///< CFL number
			     const double,
  ///< dx, cell size.
			     const double,
  ///< gas eos gamma.
			     double *,    
  ///< State vector of mean values for simulation.
			     const double,
  ///< Artificial Viscosity Parameter etav.
			     const int    
  ///< Number of tracer variables.
			);
  ~FV_solver_Hydro_Euler_Eint();

  ///
  /// Adds the contribution from flux in the current direction to dU.
  ///
  virtual int dU_Cell(
        class GridBaseClass *,
        cell *, ///< Current cell.
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
#ifdef        EINT_ETOT_PARALLEL
  ///
  /// This function calculates div(v) and if it is positive then uses the
  /// internal energy to update the pressure to a new value.  In this way
  /// rarefactions will use the internal energy update while shocks will
  /// use the total energy update.
  ///
  virtual int CellAdvanceTime(class cell *, ///< Cell to update.
			      const double *, ///< Initial Primitive State Vector.
			      double *, ///< Update vector dU
			      double *, ///< Final Primitive state vector (can be same as initial vec.).
			      double *,  ///< Tracks change of energy if I have to correct for negative pressure
			      const double, ///< gas EOS gamma.
			      const double  ///< Cell timestep dt.
			      );
#endif //     EINT_ETOT_PARALLEL

  ///
  /// Given a cell, calculate the hydrodynamic timestep.  This augments
  /// the Courant condition with a viscous timestep requirement based
  /// on div(v).
  ///
  virtual double CellTimeStep(const cell *, ///< pointer to cell
			      const double, ///< gas EOS gamma.
			      const double  ///< Cell size dx.
			      );

  ///
  /// Pre-process the grid before calc_dU() in this function.
  /// Multi-dimensional calculations to be performed on every cell
  /// before the flux calculations.  This function calls the base
  /// implementation and then calculates div(v) and the Qii values
  /// for the von Neumann/Richtmeyer (1950) viscosity.
  ///
  /// Refs:
  /// Tscharnuter & Winkler (1979) Computer Physics Communications, 18,171.
  /// Stone & Norman (1992) ApJS, 80, 753.
  ///
  virtual int preprocess_data(const int, ///< Space order of acc for this call.
			      const int  ///< Time order of acc for this call.
			      );
  
  ///
  /// This function is a place-holder for integrations which need to 
  /// add a multi-dimensional correction/addition to the update
  /// vector. Here we add the P.div(V) term in the internal--energy
  /// solver for the Euler equations, and also the von Neumann type
  /// artificial viscosity.
  ///
  virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
			     const int  ///< Time order of acc for this call.
			     );
 protected:
  ///
  /// Set the Q-viscosity parameter for a given direction on the grid.
  /// This is genuinely multi-dimensional bulk viscosity, but only the
  /// diagonal elements are calculated here so there is no shear
  /// viscosity.
  /// Refs:
  /// - Tscharnuter & Winkler (1979) Comp. Phys. Comms., 18, 171.
  /// - Stone & Norman (1992) ApJS, 80, 753.
  ///
  virtual void set_Qvisc(cell *,  ///< cell to operate on.
		 const axes, ///< axis we are looking along (Q1/2/3).
		 const double ///< value of divv (should be pre-calculated).
		 );

  ///
  /// This gets the Q-viscosity parameter from the cell's extra-data
  /// list, for the requested direction axis.
  ///
  double get_Qvisc(const cell *, ///< cell to get Q from.
		   const axes    ///< get Q for this axis.
		   );

  double Qvisc_eta; ///< viscosity parameter = C*deltaX, with C~1.
};

///
/// Solver for Euler equations in axial symmetry with AV and tracers, 
/// integrating the INTERNAL ENERGY and NOT the total energy.
///
class cyl_FV_solver_Hydro_Euler_Eint
: virtual public FV_solver_Hydro_Euler_Eint,
  virtual public VectorOps_Cyl
{
 public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  cyl_FV_solver_Hydro_Euler_Eint(const int, ///< number of variables in state vector.
			    const int, ///< number of space dimensions in grid.
			    const double, ///< CFL number
			    const double, ///< dx, cell size.
			    const double, ///< gas eos gamma.
			    double *,     ///< State vector of mean values for simulation.
			    const double, ///< Artificial Viscosity Parameter etav.
			    const int     ///< Number of tracer variables.
			    );

  ~cyl_FV_solver_Hydro_Euler_Eint();

  ///
  /// Adds the contribution from flux in the current direction to dU.
  /// Includes source terms for cylindrical coordinates.
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
  /// Set the Q-viscosity parameter for a given direction on the grid.
  /// This is genuinely multi-dimensional bulk viscosity, but only the
  /// diagonal elements are calculated here so there is no shear
  /// viscosity.  This is based on the Cartesian version, but radial
  /// gradients must be calculated differently.
  /// Refs:
  /// - Tscharnuter & Winkler (1979) Comp. Phys. Comms., 18, 171.
  /// - Stone & Norman (1992) ApJS, 80, 753.
  ///
  void set_Qvisc(cell *,  ///< cell to operate on.
		 const axes, ///< axis we are looking along (Q1/2/3).
		 const double ///< value of divv (should be pre-calculated).
		 );

  ///
  /// Multi-dimensional correction/addition to the update
  /// vector. Here we add the P.div(V) term in the internal--energy
  /// solver for the Euler equations, and also the von Neumann type
  /// artificial viscosity.
  ///
  virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
  			     const int  ///< Time order of acc for this call.
  			     );
};

///
/// Solver for Euler equations in spherical coordinates with viscosity
/// and tracers, integrating the INTERNAL ENERGY and NOT the total energy.
///
class sph_FV_solver_Hydro_Euler_Eint
: virtual public FV_solver_Hydro_Euler_Eint,
  virtual public VectorOps_Sph
{
 public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  sph_FV_solver_Hydro_Euler_Eint(const int, ///< number of variables in state vector.
			    const int, ///< number of space dimensions in grid.
			    const double, ///< CFL number
			    const double, ///< dx, cell size.
			    const double, ///< gas eos gamma.
			    double *,     ///< State vector of mean values for simulation.
			    const double, ///< Artificial Viscosity Parameter etav.
			    const int     ///< Number of tracer variables.
			    );
  ~sph_FV_solver_Hydro_Euler_Eint();

  ///
  /// Adds the contribution from flux in the current direction to dU.
  /// Includes source terms for spherical polar coordinates.
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
  /// Multi-dimensional correction/addition to the update
  /// vector. Here we add the P.div(V) term in the internal--energy
  /// solver for the Euler equations, and also the von Neumann type
  /// artificial viscosity.
  ///
  virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
  			     const int  ///< Time order of acc for this call.
  			     );

  ///
  /// Set the Q-viscosity parameter for a given direction on the grid.
  /// This is genuinely multi-dimensional bulk viscosity, but only the
  /// diagonal elements are calculated here so there is no shear
  /// viscosity.  This is based on the Cartesian version, but radial
  /// gradients must be calculated differently.  Of course in 1D it is
  /// just 1D diffusion.
  /// Refs:
  /// - Tscharnuter & Winkler (1979) Comp. Phys. Comms., 18, 171.
  /// - Stone & Norman (1992) ApJS, 80, 753.
  /// - Boss (2006) ApJ, 641,1148.
  ///
  void set_Qvisc(cell *,  ///< cell to operate on.
		 const axes, ///< axis we are looking along (Q1/2/3).
		 const double ///< value of divv (should be pre-calculated).
		 );
};






#endif // if INCLUDE_EINT_ADI_HYDRO

#endif // SOLVER_EQN_HYDRO_ADI_EINT_H
