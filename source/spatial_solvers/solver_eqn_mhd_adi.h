///
/// \file solver_eqn_hydro_adi.h
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via
/// either Lax-Friedrichs or Riemann solver (linear and/or exact).
/// Adds viscosity if asked to, and tracks flux of N passive tracers.
///
/// Modifications:
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///  - 2010.11.15 JM:
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2015.01.15 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#ifndef SOLVER_EQN_MHD_ADI_H
#define SOLVER_EQN_MHD_ADI_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "spatial_solvers/solver_eqn_base.h"
#include "flux_calc/flux_mhd_adiabatic.h"

///
/// The main solver the code uses for integrating the ideal MHD Equations (adiabatic).
/// This will eventually do a few types of flux solution:
///  - 0=Lax-Friedrichs,
///  - 1=Linear Riemann,
///  - 2=Exact Riemann,
///  - 3=Hybrid Riemann (Linear/Exact)
///  - 4=Roe-Averaged approximate solver.
///
///
class FV_solver_mhd_ideal_adi
: virtual public FV_solver_base,
  virtual public flux_solver_mhd_ideal_adi,
  virtual public VectorOps_Cart
{
 public:
  FV_solver_mhd_ideal_adi(
        const int, ///< number of variables in state vector.
        const int, ///< number of space dimensions in grid.
        const double, ///< CFL number
        const double, ///< dx, cell size.
        const double, ///< gas eos gamma.
        pion_flt *,     ///< State vector of mean values for simulation.
        const double, ///< Artificial Viscosity Parameter etav.
        const int     ///< Number of tracer variables.
        );
  ~FV_solver_mhd_ideal_adi();

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
        class cell *c, ///< cell to update.
        const pion_flt *, ///< Initial Primitive State Vector.
        pion_flt *, ///< Update vector dU
        pion_flt *, ///< Final Primitive state vector (can be same as initial vec.).
        pion_flt *,  ///< Tracks change of energy if I have to correct for negative pressure
        const double, ///< gas EOS gamma.
        const double  ///< Cell timestep dt.
        );
  ///
  /// Given a cell, calculate the MHD/hydrodynamic timestep.
  ///
  virtual double CellTimeStep(
        const cell *, ///< pointer to cell
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
  virtual void PtoU(
        const pion_flt *, ///< pointer to Primitive variables.
        pion_flt *,       ///< pointer to conserved variables.
        const double    ///< Gas constant gamma.
        );
  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The 
  /// conserved variable is the mass density of this.
  ///
  virtual int UtoP(
        const pion_flt *, ///< pointer to conserved variables.
        pion_flt *, ///< pointer to Primitive variables.
        const double    ///< Gas constant gamma.
        );
  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// The flux of a passive tracer is equal to the mass flux times 
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void PUtoFlux(
        const pion_flt *, ///< pointer to Primitive variables.
        const pion_flt *, ///< pointer to conserved variables.
        pion_flt *  ///< Pointer to flux variable.
        );
  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// The flux of a passive tracer is equal to the mass flux times 
  /// the value of the primitive tracer variable.  I take the left
  /// state tracer var. if the mass flux is to the right, and vice versa.
  ///
  virtual void UtoFlux(
        const pion_flt *, ///< Pointer to conserved variables state vector.
        pion_flt *,       ///< Pointer to flux variable state vector.
        const double   ///< Gas constant gamma.
        );
};


///
/// The main solver the code uses for integrating the ideal MHD Equations (adiabatic),
/// with the mixed GLM method of divergence cleaning (Dedner et al., 2002).
///
/// This will eventually do a few types of flux solution:
///  - 0=Lax-Friedrichs,
///  - 1=Linear Riemann,
///  - 2=Exact Riemann,
///  - 3=Hybrid Riemann (Linear/Exact)
///  - 4=Roe-Averaged approximate solver.
///
///
class FV_solver_mhd_mixedGLM_adi
  : virtual public FV_solver_mhd_ideal_adi,
  virtual public flux_solver_mhd_mixedGLM_adi,
  virtual public VectorOps_Cart
{
  public:
   FV_solver_mhd_mixedGLM_adi(
        const int, ///< number of variables in state vector.
        const int, ///< number of space dimensions in grid.
        const double, ///< CFL number
        const double, ///< dx, cell size.
        const double, ///< gas eos gamma.
        pion_flt *,     ///< State vector of mean values for simulation.
        const double, ///< Artificial Viscosity Parameter etav.
        const int     ///< Number of tracer variables.
        );
   ~FV_solver_mhd_mixedGLM_adi();
   
   ///
   /// General Finite volume scheme for updating a cell's
   /// primitive state vector, for homogeneous equations.
   /// For the mixedGLM equations, we need to add a source term to the update of Psi.
   ///
   virtual int CellAdvanceTime(
        class cell *c, ///< cell to update.
        const pion_flt *, ///< Initial Primitive State Vector.
        pion_flt *, ///< Update vector dU
        pion_flt *, ///< Final Primitive state vector (can be same as initial vec.).
        pion_flt *,  ///< Tracks change of energy if I have to correct for negative pressure
        const double, ///< gas EOS gamma.
        const double  ///< Cell timestep dt.
        );
  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The 
  /// conserved variable is the mass density of this.
  ///
  virtual void PtoU(
        const pion_flt *, ///< pointer to Primitive variables.
        pion_flt *,       ///< pointer to conserved variables.
        const double    ///< Gas constant gamma.
        );
  ///
  /// This calls the original version and then adds conversion of tracer variables.
  /// 
  /// For purely passive tracers, the primitive variable is just a number,
  /// such as the 'colour' of the gas, or where it started out.  The 
  /// conserved variable is the mass density of this.
  ///
  virtual int UtoP(
        const pion_flt *, ///< pointer to conserved variables.
        pion_flt *, ///< pointer to Primitive variables.
        const double    ///< Gas constant gamma.
        );

   ///
   /// This sets the values of GLM_chyp and GLM_cr based on the timestep.
   ///
   virtual void GotTimestep(
        const double  ///< timestep dt.
        );
};


/// ---------------------------------------------------------------------
/// -------------------  AXI-SYMMETRIC EQUATIONS ------------------------
/// ---------------------------------------------------------------------


/** \brief Solver for Ideal MHD equations in axial symmetry with AV and tracers. */
class cyl_FV_solver_mhd_ideal_adi
  : virtual public FV_solver_mhd_ideal_adi, virtual public VectorOps_Cyl
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   cyl_FV_solver_mhd_ideal_adi(
        const int, ///< number of variables in state vector.
        const int, ///< number of space dimensions in grid.
        const double, ///< CFL number
        const double, ///< dx, cell size.
        const double, ///< gas eos gamma.
        pion_flt *,     ///< State vector of mean values for simulation.
        const double, ///< Artificial Viscosity Parameter etav.
        const int     ///< Number of tracer variables.
        );
   ~cyl_FV_solver_mhd_ideal_adi();
   /** \brief Adds the contribution from flux in the current direction to dU.
    * 
    * Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
    * spatial accuracy.
    */
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

/** \brief Solver for mixed-GLM MHD equations with AV and tracers, in 
 * axial symmetry. */
class cyl_FV_solver_mhd_mixedGLM_adi
  : virtual public FV_solver_mhd_mixedGLM_adi, virtual public VectorOps_Cyl
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   cyl_FV_solver_mhd_mixedGLM_adi(
        const int, ///< number of variables in state vector.
        const int, ///< number of space dimensions in grid.
        const double, ///< CFL number
        const double, ///< dx, cell size.
        const double, ///< gas eos gamma.
        pion_flt *,     ///< State vector of mean values for simulation.
        const double, ///< Artificial Viscosity Parameter etav.
        const int     ///< Number of tracer variables.
        );
   ~cyl_FV_solver_mhd_mixedGLM_adi();
   /** \brief Adds the contribution from flux in the current direction to dU.
    * 
    * Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
    * spatial accuracy.
    */
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


#endif // SOLVER_EQN_MHD_ADI_H
