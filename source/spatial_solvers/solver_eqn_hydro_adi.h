///
/// \file solver_eqn_hydro_adi.h
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via either
/// Lax-Friedrichs or Riemann solver (linear and/or exact).  Also tracks flux of
/// N passive tracers.
///
/// - 2009-12-18 JM: Added Axisymmetric Class (cyl_FV_solver_Hydro_Euler)
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux
///  functions).
/// - 2010.10.01 JM: Added spherical coordinate system.
/// - 2010.11.15 JM:
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.21 JM: updated documentation.
/// - 2010.12.23 JM: Removed references to riemann_base class.
///    added extra variable to inviscid_flux() function.
///    Moved UtoP() etc. from solver to flux-solver.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.04.14 JM: Moved flux solver to FV_solver

#ifndef SOLVER_EQN_HYDRO_ADI_H
#define SOLVER_EQN_HYDRO_ADI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "Riemann_solvers/HLL_hydro.h"
#include "Riemann_solvers/Riemann_FVS_hydro.h"
#include "Riemann_solvers/Roe_Hydro_ConservedVar_solver.h"
#include "Riemann_solvers/Roe_Hydro_PrimitiveVar_solver.h"
#include "Riemann_solvers/riemann.h"
#include "coord_sys/VectorOps_spherical.h"
#include "equations/eqns_base.h"
#include "grid/cell_interface.h"  // to get the 'cell' class.
#include "spatial_solvers/solver_eqn_base.h"

///
/// The main solver the code uses for integrating the Euler Equations
/// (adiabatic). This is the most functional solver, and can do many types of
/// flux calculation: 0= Lax-Friedrichs flux: Very diffusive, 1st order. 1=
/// Linear Riemann:      Uses arithmetic mean of left and right states.  Not
/// great. 2= Exact Riemann:       Very slow, but the best solver. 3= Hybrid
/// Riemann (Linear/Exact): uses linear solver when possible. 4= Roe conserved
/// variable solver: calculates flux from e-values/e-vectors. 5= Roe primitive
/// variable solver: linear solver using the Roe-average state. 6= van Leer's
/// Flux Vector Splitting: very robust, diffusive contact discontinuity.
///
class FV_solver_Hydro_Euler :
    virtual public FV_solver_base,
    virtual public riemann_Euler,
    virtual public Riemann_FVS_Euler,
    virtual public Riemann_Roe_Hydro_PV,
    virtual public Riemann_Roe_Hydro_CV,
    virtual public HLL_hydro,
    virtual public VectorOps_Cart {
public:
  FV_solver_Hydro_Euler(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt*,     ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );

  ~FV_solver_Hydro_Euler();
  long int counter;

  ///
  /// This calls the equations version and then adds conversion of tracer
  /// variables.
  ///
  /// For passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The
  /// conserved variable is the mass density of this.
  ///
  virtual void PtoU(
      const pion_flt*,  ///< pointer to Primitive variables.
      pion_flt*,        ///< pointer to conserved variables.
      const double      ///< Gas constant gamma.
  );

  ///
  /// This calls the equations version and then adds conversion of tracer
  /// variables.
  ///
  /// For passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The
  /// conserved variable is the mass density of this.
  ///
  virtual int UtoP(
      const pion_flt*,  ///< pointer to conserved variables.
      pion_flt*,        ///< pointer to Primitive variables.
      const double,     ///< minimum temperature/pressure allowed
      const double      ///< Gas constant gamma.
  );

  ///
  /// This calls the equations version and then adds conversion of tracer
  /// variables.
  ///
  /// The flux of a passive tracer is equal to the mass flux times
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void PUtoFlux(
      const pion_flt*,  ///< pointer to Primitive variables.
      const pion_flt*,  ///< pointer to conserved variables.
      pion_flt*         ///< Pointer to flux variable.
  );

  ///
  /// This calls the equations version and then adds conversion of tracer
  /// variables.
  ///
  /// The flux of a passive tracer is equal to the mass flux times
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void UtoFlux(
      const pion_flt*,  ///< Pointer to conserved variables state vector.
      pion_flt*,        ///< Pointer to flux variable state vector.
      const double      ///< Gas constant gamma.
  );

  /// Calculates Flux based on a left and right state vector (primitive).
  int inviscid_flux(
      class SimParams&,      ///< simulation parameters
      class GridBaseClass*,  ///< pointer to grid
      const double,          ///< cell-size dx (for LF method)
      class cell*,           ///< Left state cell pointer
      class cell*,           ///< Right state cell pointer
      const pion_flt*,       ///< Left Primitive state vector.
      const pion_flt*,       ///< Right Primitive state vector.
      pion_flt*,             ///< Resultant Flux state vector.
      pion_flt*,             ///< Resultant Pstar state vector.
      const int,             ///< Which Riemann solver
      const double           ///< Gas constant gamma.
  );

  ///
  /// Adds the contribution from flux in the current direction to dU.
  ///
  virtual int dU_Cell(
      class GridBaseClass*,
      cell*,            ///< Current cell.
      const axes,       ///< Which axis we are looking along.
      const pion_flt*,  ///< Negative direction flux.
      const pion_flt*,  ///< Positive direction flux.
      const pion_flt*,  ///< slope vector for cell c.
      const int,        ///< spatial order of accuracy.
      const double,     ///< cell length dx.
      const double      ///< cell TimeStep, dt.
  );

  ///
  /// Geometric source terms (does nothing for Cartesian geometry).
  ///
  virtual void geometric_source(
      cell*,            ///< Current cell.
      const axes,       ///< Which axis we are looking along.
      const pion_flt*,  ///< slope vector for cell c.
      const int,        ///< spatial order of accuracy.
      const double,     ///< cell length dx.
      pion_flt*         ///< update vector to add source term to [OUTPUT]
  )
  {
    return;
  }

  ///
  /// General Finite volume scheme for updating a cell's
  /// primitive state vector, for homogeneous equations.
  ///
  virtual int CellAdvanceTime(
      class cell*,      ///< current cell.
      const pion_flt*,  ///< Initial Primitive State Vector.
      pion_flt*,        ///< Update vector dU
      pion_flt*,     ///< Final Primitive state vector (can be same as initial
                     ///< vec.).
      pion_flt*,     ///< Tracks change of energy if I have to correct for
                     ///< negative pressure
      const double,  ///< gas EOS gamma.
      const double,  ///< Min Temperature allowed on grid.
      const double   ///< Cell timestep dt.
  );

  ///
  /// Given a cell, calculate the hydrodynamic timestep.
  ///
  virtual double CellTimeStep(
      const cell*,   ///< pointer to cell
      const double,  ///< gas EOS gamma.
      const double   ///< Cell size dx.
  );

protected:
  ///
  /// Falle et al. (1998) Artificial Viscosity Calculation.
  ///
  int AVFalle(
      const pion_flt*,  ///< Left Primitive state vector.
      const pion_flt*,  ///< Right Primitive state vector.
      const pion_flt*,  ///< Resolved (P*) state vector.
      pion_flt*,        ///< Pointer to associated Flux Vector.
      const double,     ///< Artificial Viscosity parameter, etav.
      const double      ///< gamma
  );
};

/// Solver for Euler equations in cylindrical coordinates with AV
/// and tracers.
class cyl_FV_solver_Hydro_Euler :
    virtual public FV_solver_Hydro_Euler,
    virtual public VectorOps_Cyl {
public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  cyl_FV_solver_Hydro_Euler(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt*,     ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );

  ~cyl_FV_solver_Hydro_Euler();

  ///
  /// Geometric source terms.
  /// Includes geometric source term p/R for 1st and 2nd order
  /// spatial accuracy.
  ///
  virtual void geometric_source(
      cell*,            ///< Current cell.
      const axes,       ///< Which axis we are looking along.
      const pion_flt*,  ///< slope vector for cell c.
      const int,        ///< spatial order of accuracy.
      const double,     ///< cell length dx.
      pion_flt*         ///< update vector to add source term to [OUTPUT]
  );
};

///
/// Solver for Euler equations in spherical coordinates with viscosity and
/// tracers.
///
class sph_FV_solver_Hydro_Euler :
    virtual public FV_solver_Hydro_Euler,
    virtual public VectorOps_Sph {
public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  sph_FV_solver_Hydro_Euler(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt*,     ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );

  ~sph_FV_solver_Hydro_Euler();

  ///
  /// Geometric source terms.
  /// Includes geometric source term p/R for 1st and 2nd order
  /// spatial accuracy.
  ///
  virtual void geometric_source(
      cell*,            ///< Current cell.
      const axes,       ///< Which axis we are looking along.
      const pion_flt*,  ///< slope vector for cell c.
      const int,        ///< spatial order of accuracy.
      const double,     ///< cell length dx.
      pion_flt*         ///< update vector to add source term to [OUTPUT]
  );

  ///	Vector for corrector values (used to modify flux according to sCMA)
  std::vector<double> corrector;
};

#endif  // SOLVER_EQN_HYDRO_ADI_H
