///
/// \file solver_eqn_hydro_adi.h
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via either Lax-Friedrichs
/// or Riemann solver (linear and/or exact).  Also tracks flux
/// of N passive tracers.
///
/// - 2009-12-18 JM: Added Axisymmetric Class (cyl_FV_solver_Hydro_Euler)
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
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

#ifndef SOLVER_EQN_HYDRO_ADI_H
#define SOLVER_EQN_HYDRO_ADI_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "spatial_solvers/solver_eqn_base.h"
#include "flux_calc/flux_hydro_adiabatic.h"
#include "coord_sys/VectorOps_spherical.h"

///
/// The main solver the code uses for integrating the Euler Equations (adiabatic).
/// This is the most functional solver, and can do many types of flux calculation:
/// 0= Lax-Friedrichs flux: Very diffusive, 1st order.
/// 1= Linear Riemann:      Uses arithmetic mean of left and right states.  Not great.
/// 2= Exact Riemann:       Very slow, but the best solver.
/// 3= Hybrid Riemann (Linear/Exact): uses linear solver when possible.
/// 4= Roe conserved variable solver: calculates flux from e-values/e-vectors.
/// 5= Roe primitive variable solver: linear solver using the Roe-average state.
/// 6= van Leer's Flux Vector Splitting: very robust, diffusive contact discontinuity.
///
class FV_solver_Hydro_Euler : virtual public FV_solver_base, virtual public flux_solver_hydro_adi, virtual public VectorOps_Cart
{
 public:
  FV_solver_Hydro_Euler(
      const int, ///< number of variables in state vector.
      const int, ///< number of space dimensions in grid.
      const double, ///< CFL number
      const double, ///< dx, cell size.
      const double, ///< gas eos gamma.
      pion_flt *,     ///< State vector of mean values for simulation.
      const double, ///< Artificial Viscosity Parameter etav.
      const int     ///< Number of tracer variables.
      );

  ~FV_solver_Hydro_Euler();

  ///
  /// Adds the contribution from flux in the current direction to dU.
  ///
  virtual int dU_Cell(
        class GridBaseClass *,
        cell *, ///< Current cell.
        const axes, ///< Which axis we are looking along.
        const pion_flt *, ///< Negative direction flux.
        const pion_flt *, ///< Positive direction flux.
        const pion_flt *, ///< slope vector for cell c.
        const int,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        );

  ///
  /// General Finite volume scheme for updating a cell's
  /// primitive state vector, for homogeneous equations.
  ///
  virtual int CellAdvanceTime(
      class cell *,   ///< current cell.
      const pion_flt *, ///< Initial Primitive State Vector.
      pion_flt *, ///< Update vector dU
      pion_flt *, ///< Final Primitive state vector (can be same as initial vec.).
      pion_flt *,  ///< Tracks change of energy if I have to correct for negative pressure
      const double, ///< gas EOS gamma.
      const double  ///< Cell timestep dt.
      );

  ///
  /// Given a cell, calculate the hydrodynamic timestep.
  ///
  virtual double CellTimeStep(
      const cell *, ///< pointer to cell
      const double, ///< gas EOS gamma.
      const double  ///< Cell size dx.
      );
};


/** \brief Solver for Euler equations in axial symmetry with AV and tracers. */
class cyl_FV_solver_Hydro_Euler
: virtual public FV_solver_Hydro_Euler,
  virtual public VectorOps_Cyl
{
 public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  cyl_FV_solver_Hydro_Euler(
      const int, ///< number of variables in state vector.
      const int, ///< number of space dimensions in grid.
      const double, ///< CFL number
      const double, ///< dx, cell size.
      const double, ///< gas eos gamma.
      pion_flt *,     ///< State vector of mean values for simulation.
      const double, ///< Artificial Viscosity Parameter etav.
      const int     ///< Number of tracer variables.
      );

  ~cyl_FV_solver_Hydro_Euler();
  ///
  /// Adds the contribution from flux in the current direction to dU.
  ///
  virtual int dU_Cell(
        class GridBaseClass *,
        cell *, ///< Current cell.
        const axes, ///< Which axis we are looking along.
        const pion_flt *, ///< Negative direction flux.
        const pion_flt *, ///< Positive direction flux.
        const pion_flt *, ///< slope vector for cell c.
        const int,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        );
};

///
/// Solver for Euler equations in spherical coordinates with viscosity and tracers.
///
class sph_FV_solver_Hydro_Euler
: virtual public FV_solver_Hydro_Euler,
  virtual public VectorOps_Sph
{
 public:
  ///
  /// sets indices for tracer variables in state vector.
  ///
  sph_FV_solver_Hydro_Euler(
      const int, ///< number of variables in state vector.
      const int, ///< number of space dimensions in grid.
      const double, ///< CFL number
      const double, ///< dx, cell size.
      const double, ///< gas eos gamma.
      pion_flt *,     ///< State vector of mean values for simulation.
      const double, ///< Artificial Viscosity Parameter etav.
      const int     ///< Number of tracer variables.
      );

  ~sph_FV_solver_Hydro_Euler();
  ///
  /// Adds the contribution from flux in the current direction to dU.
  /// Includes source terms for spherical polar coordinates.
  ///
  virtual int dU_Cell(
        class GridBaseClass *,
        cell *, ///< Current cell.
        const axes, ///< Which axis we are looking along.
        const pion_flt *, ///< Negative direction flux.
        const pion_flt *, ///< Positive direction flux.
        const pion_flt *, ///< slope vector for cell c.
        const int,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        );
};

#endif // SOLVER_EQN_HYDRO_ADI_H

