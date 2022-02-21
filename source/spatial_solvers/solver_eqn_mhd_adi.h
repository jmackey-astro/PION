///
/// \file solver_eqn_mhd_adi.h
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via
/// either Lax-Friedrichs or Riemann solver (linear and/or exact).
/// Adds viscosity if asked to, and tracks flux of N passive tracers.
///
/// Modifications:
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux
///  functions).
///  - 2010.11.15 JM:
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2015.01.15 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.04.14 JM: Moved flux solver to FV_solver

#ifndef SOLVER_EQN_MHD_ADI_H
#define SOLVER_EQN_MHD_ADI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "Riemann_solvers/HLLD_MHD.h"
#include "Riemann_solvers/Roe_MHD_ConservedVar_solver.h"
#include "Riemann_solvers/riemannMHD.h"
#include "equations/eqns_mhd_adiabatic.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"

///
/// The main solver the code uses for integrating the ideal MHD Equations
/// (adiabatic). This will eventually do a few types of flux solution:
///  - 0=Lax-Friedrichs,
///  - 1=Linear Riemann,
///  - 2=Exact Riemann,
///  - 3=Hybrid Riemann (Linear/Exact)
///  - 4=Roe-Averaged approximate solver.
///
///
class FV_solver_mhd_ideal_adi :
    virtual public FV_solver_base,
    virtual public riemann_MHD,
    virtual public Riemann_Roe_MHD_CV,
    virtual public HLLD_MHD,
    virtual public VectorOps_Cart {
public:
  FV_solver_mhd_ideal_adi(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt *,    ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );
  ~FV_solver_mhd_ideal_adi();

  /// Calculates Flux based on a left and right state vector (primitive).
  virtual int inviscid_flux(
      class SimParams &,      ///< simulation parameters
      class GridBaseClass *,  ///< pointer to grid
      const double,           ///< cell-size dx (for LF method)
      class cell &,           ///< Left state cell pointer
      class cell &,           ///< Right state cell pointer
      const pion_flt *,       ///< Left Primitive state vector.
      const pion_flt *,       ///< Right Primitive state vector.
      pion_flt *,             ///< Resultant Flux state vector.
      pion_flt *,             ///< State vector at interface.
      const int,              ///< Which Riemann solver
      const double            ///< Gas constant gamma.
  );

  ///
  /// Adds the contribution from flux in the current direction to dU.
  ///
  virtual int dU_Cell(
      class GridBaseClass *,
      cell &,            ///< Current cell.
      const axes,        ///< Which axis we are looking along.
      const pion_flt *,  ///< Negative direction flux.
      const pion_flt *,  ///< Positive direction flux.
      const pion_flt *,  ///< slope vector for cell c.
      const int,         ///< spatial order of accuracy.
      const double,      ///< cell length dx.
      const double       ///< cell TimeStep, dt.
  );

  ///
  /// calculate Powell and GLM source terms for multi-D MHD
  ///
  virtual int MHDsource(
      class GridBaseClass *,  ///< pointer to grid.
      class cell &,           ///< pointer to cell of left state
      class cell &,           ///< pointer to cell of right state
      pion_flt *,             ///< left edge state
      pion_flt *,             ///< right edge state
      const axes,             ///< Which axis we are looking along.
      enum direction,         ///< positive direction normal to interface
      enum direction,         ///< negative direction normal to interface
      const double            ///< timestep dt
  );

  ///
  /// Geometric source terms (does nothing for Cartesian geometry).
  ///
  virtual void geometric_source(
      cell &,            ///< Current cell.
      const axes,        ///< Which axis we are looking along.
      const pion_flt *,  ///< slope vector for cell c.
      const int,         ///< spatial order of accuracy.
      const double,      ///< cell length dx.
      pion_flt *         ///< update vector to add source term to [OUTPUT]
  )
  {
    return;
  }

  ///
  /// General Finite volume scheme for updating a cell's
  /// primitive state vector, for homogeneous equations.
  ///
  virtual int CellAdvanceTime(
      class cell &c,     ///< cell to update.
      const pion_flt *,  ///< Initial Primitive State Vector.
      pion_flt *,        ///< Update vector dU
      pion_flt *,    ///< Final Primitive state vector (can be same as initial
                     ///< vec.).
      pion_flt *,    ///< Tracks change of energy if I have to correct for
                     ///< negative pressure
      const double,  ///< gas EOS gamma.
      const double,  ///< Min Temperature allowed on grid.
      const double   ///< Cell timestep dt.
  );

  ///
  /// Given a cell, calculate the MHD/hydrodynamic timestep.
  ///
  virtual double CellTimeStep(
      const cell &,  ///< pointer to cell
      const double,  ///< gas EOS gamma.
      const double   ///< Cell size dx.
  );

  ///
  /// This calls the original version and then adds conversion of tracer
  /// variables.
  ///
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The
  /// conserved variable is the mass density of this.
  ///
  virtual void PtoU(
      const pion_flt *,  ///< pointer to Primitive variables.
      pion_flt *,        ///< pointer to conserved variables.
      const double       ///< Gas constant gamma.
  );

  ///
  /// This calls the original version and then adds conversion of tracer
  /// variables.
  ///
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The
  /// conserved variable is the mass density of this.
  ///
  virtual int UtoP(
      const pion_flt *,  ///< pointer to conserved variables.
      pion_flt *,        ///< pointer to Primitive variables.
      const double,      ///< Min Temperature allowed on grid.
      const double       ///< Gas constant gamma.
  );

  ///
  /// This calls the original version and then adds conversion of tracer
  /// variables.
  ///
  /// The flux of a passive tracer is equal to the mass flux times
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void PUtoFlux(
      const pion_flt *,  ///< pointer to Primitive variables.
      const pion_flt *,  ///< pointer to conserved variables.
      pion_flt *         ///< Pointer to flux variable.
  );
  ///
  /// This calls the original version and then adds conversion of tracer
  /// variables.
  ///
  /// The flux of a passive tracer is equal to the mass flux times
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void UtoFlux(
      const pion_flt *,  ///< Pointer to conserved variables state vector.
      pion_flt *,        ///< Pointer to flux variable state vector.
      const double       ///< Gas constant gamma.
  );

protected:
  double max_speed;  ///< max. fast magnetosonic speed on the domain.

  ///
  /// Falle, Komissarov & Joarder (1998,MNRAS,297,265) Artificial
  /// Viscosity Calculation (one-dimensional).
  ///
  int AVFalle(
      const pion_flt *,  ///< Left Primitive state vector.
      const pion_flt *,  ///< Right Primitive state vector.
      const pion_flt *,  ///< Resolved (P*) state vector.
      pion_flt *,        ///< Pointer to associated Flux Vector.
      const double,      ///< Artificial Viscosity parameter, etav.
      const double       ///< gamma
  );

  /// shut off reporting if we get more than 1000 negative pressures.
  long int negPGct;
  /// shut off reporting if we get more than 1000 negative densities.
  long int negROct;
};

// **********************************************************************************
// For the Dedner-GLM divergence cleaning method.
// **********************************************************************************

///
/// The main solver the code uses for integrating the ideal MHD Equations
/// (adiabatic), with the mixed GLM method of divergence cleaning (Dedner et
/// al., 2002).
///
/// This will eventually do a few types of flux solution:
///  - 0=Lax-Friedrichs,
///  - 1=Linear Riemann,
///  - 2=Exact Riemann,
///  - 3=Hybrid Riemann (Linear/Exact)
///  - 4=Roe-Averaged approximate solver.
///
///
class FV_solver_mhd_mixedGLM_adi :
    virtual public FV_solver_mhd_ideal_adi,
    virtual public eqns_mhd_mixedGLM,
    virtual public VectorOps_Cart {
public:
  FV_solver_mhd_mixedGLM_adi(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt *,    ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );
  ~FV_solver_mhd_mixedGLM_adi();

  ///
  /// Given a cell, calculate the MHD/hydrodynamic timestep.
  /// This calculates the ideal MHD step, and also tracks the maximum
  /// velocity on the domain.
  ///
  virtual double CellTimeStep(
      const cell &,  ///< pointer to cell
      const double,  ///< gas EOS gamma.
      const double   ///< Cell size dx.
  );

  ///
  /// General Finite volume scheme for updating a cell's
  /// primitive state vector, for homogeneous equations.
  /// For the mixedGLM equations, we need to add a source term to the
  /// update of Psi and the total energy.
  ///
  virtual int CellAdvanceTime(
      class cell &c,     ///< cell to update.
      const pion_flt *,  ///< Initial Primitive State Vector.
      pion_flt *,        ///< Update vector dU
      pion_flt *,    ///< Final Primitive state vector (can be same as initial
                     ///< vec.).
      pion_flt *,    ///< Tracks change of energy if I have to correct for
                     ///< negative pressure
      const double,  ///< gas EOS gamma.
      const double,  ///< Min Temperature allowed on grid.
      const double   ///< Cell timestep dt.
  );

  ///
  /// This sets the values of GLM_chyp and GLM_cr based on the timestep.
  ///
  virtual void Set_GLM_Speeds(
      const double,  ///< timestep dt.
      const double,  ///< grid cell size dx.
      const double   ///< GLM damping coefficient c_r
  );

  ///
  /// Calculates Flux based on a left and right state vector (primitive).
  /// This is the same as for ideal MHD except that we use Dedner et al.
  /// (2002)'s method to calculate the flux in BX and in the extra scalar
  /// field Psi.
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
  /// \f[ F(\psi) = c_h^2 B_x^* = c_h^2 \left( \frac{1}{2}(B_x(L)+B_x(R)) -
  /// \frac{1}{2c_h}(\psi_R-\psi_L) \right) \f] \f[ F(B_x) = \psi_* =
  /// \frac{1}{2}(\psi_L+\psi_R) - \frac{c_h}{2}(B_x(R)-B_X(L)) \f]
  ///
  virtual int inviscid_flux(
      class SimParams &,      ///< simulation parameters
      class GridBaseClass *,  ///< pointer to grid
      const double,           ///< cell-size dx (for LF method)
      class cell &,           ///< Left state cell pointer
      class cell &,           ///< Right state cell pointer
      const pion_flt *,       ///< Left Primitive state vector.
      const pion_flt *,       ///< Right Primitive state vector.
      pion_flt *,             ///< Resultant Flux state vector.
      pion_flt *,             ///< State vector at interface.
      const int,              ///< Which Riemann Solver
      const double            ///< Gas constant gamma.
  );

  virtual int MHDsource(
      class GridBaseClass *,  ///< pointer to grid.
      class cell &,           ///< pointer to cell of left state
      class cell &,           ///< pointer to cell of right state
      pion_flt *,             ///< left edge state
      pion_flt *,             ///< right edge state
      const axes,             ///< Which axis we are looking along.
      enum direction,         ///< positive direction normal to interface
      enum direction,         ///< negative direction normal to interface
      const double            ///< timestep dt
  );

  ///
  /// Same as ideal MHD version except that total energy contains
  /// contribution from psi, and includes psi conversion.
  ///
  virtual void PtoU(
      const pion_flt *,  ///< pointer to Primitive variables.
      pion_flt *,        ///< pointer to conserved variables.
      const double       ///< Gas constant gamma.
  );

  ///
  /// Same as ideal MHD version except that total energy contains
  /// contribution from psi.
  /// includes psi conversion.
  ///
  virtual int UtoP(
      const pion_flt *,  ///< pointer to conserved variables.
      pion_flt *,        ///< pointer to Primitive variables.
      const double,      ///< Min Temperature allowed on grid.
      const double       ///< Gas constant gamma.
  );
};

/// ---------------------------------------------------------------------
/// -------------------  AXI-SYMMETRIC EQUATIONS ------------------------
/// ---------------------------------------------------------------------

// ##################################################################
// ##################################################################

/// Solver for Ideal MHD equations in axial symmetry with AV and
/// tracers.
class cyl_FV_solver_mhd_ideal_adi :
    virtual public FV_solver_mhd_ideal_adi,
    virtual public VectorOps_Cyl {
public:
  /// sets indices for tracer variables in state vector.
  cyl_FV_solver_mhd_ideal_adi(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt *,    ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );

  ~cyl_FV_solver_mhd_ideal_adi();

  ///
  /// Geometric source terms.
  /// Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
  /// spatial accuracy.
  ///
  virtual void geometric_source(
      cell &,            ///< Current cell.
      const axes,        ///< Which axis we are looking along.
      const pion_flt *,  ///< slope vector for cell c.
      const int,         ///< spatial order of accuracy.
      const double,      ///< cell length dx.
      pion_flt *         ///< update vector to add source term to [OUTPUT]
  );

  ///
  /// calculate Powell and GLM source terms for multi-D MHD,
  /// using the cylindrical coordinate divergence operator
  ///
  virtual int MHDsource(
      class GridBaseClass *,  ///< pointer to grid.
      class cell &,           ///< pointer to cell of left state
      class cell &,           ///< pointer to cell of right state
      pion_flt *,             ///< left edge state
      pion_flt *,             ///< right edge state
      const axes,             ///< Which axis we are looking along.
      enum direction,         ///< positive direction normal to interface
      enum direction,         ///< negative direction normal to interface
      const double            ///< timestep dt
  );
};

// ##################################################################
// ##################################################################

/// Solver for mixed-GLM MHD equations with AV and tracers, in
/// axial symmetry.
class cyl_FV_solver_mhd_mixedGLM_adi :
    virtual public FV_solver_mhd_mixedGLM_adi,
    virtual public cyl_FV_solver_mhd_ideal_adi,
    virtual public VectorOps_Cyl {
public:
  /// sets indices for tracer variables in state vector.
  cyl_FV_solver_mhd_mixedGLM_adi(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      pion_flt *,    ///< State vector of mean values for simulation.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );

  ~cyl_FV_solver_mhd_mixedGLM_adi();

  ///
  /// Geometric source terms.
  /// Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
  /// spatial accuracy, and the GLM Psi source terms.
  ///
  virtual void geometric_source(
      cell &,            ///< Current cell.
      const axes,        ///< Which axis we are looking along.
      const pion_flt *,  ///< slope vector for cell c.
      const int,         ///< spatial order of accuracy.
      const double,      ///< cell length dx.
      pion_flt *         ///< update vector to add source term to [OUTPUT]
  );

  ///
  /// calculate Powell and GLM source terms for multi-D MHD,
  /// using the cylindrical coordinate divergence operator in the
  /// Powell terms
  ///
  virtual int MHDsource(
      class GridBaseClass *,  ///< pointer to grid.
      class cell &,           ///< pointer to cell of left state
      class cell &,           ///< pointer to cell of right state
      pion_flt *,             ///< left edge state
      pion_flt *,             ///< right edge state
      const axes,             ///< Which axis we are looking along.
      enum direction,         ///< positive direction normal to interface
      enum direction,         ///< negative direction normal to interface
      const double            ///< timestep dt
  );
};

#endif  // SOLVER_EQN_MHD_ADI_H
