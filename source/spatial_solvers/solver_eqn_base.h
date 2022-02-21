/// \file solver_eqn_base.h
///
/// \brief Class declarations for various solvers implemented in my FV
/// grid-code. \author Jonathan Mackey
///
/// Modifications:
///  - 2007-07-11 Still writing it and working out class hierarchy.
///  - 2007-07-12 Got the class heirarchy working (I think...).
///  - 2007-07-13 New Class structure implemented and working
///  - 2007-07-16 Reverted to less complicated class hierarchy. works well.
///  - 2007-07-23 Started to add passive tracer variables
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08-07 cylindrical coordinates support working.
///  - 2007-08-08 cyl.coords for iMHD and glmMHD.
///  - 2007-12-11 Field-CD method 'working' in 2D cartesian coords.
///  - 2009-10-20 New structure built on equations and flux classes...
///  - 2010.09.30 JM: Worked on Lapidus AV (added set_div_v() function for
///  cells)
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2010.11.15 JM: Made InterCellFlux general for all classes.
///   Added H-correction functions for calculating eta and getting
///   the maximum value of eta on a given H-stencil.
///   Added Pre- and Post-flux viscosity functions for the
///   H-correction and Lapidus viscosity functions.
/// - 2010.12.27 JM: changed interface to post_calc_viscous_terms so
///   that it takes left[]/right[]/pstar[] as parameters.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2011.01.03 JM: Documented postprocess_du() function.
///   Moved preprocess_data() and calc_Hcorrection from gridMethods.cc
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2015.01.14 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.07.06 JM: tidied up a bit (but much more to do!)
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.01.24 JM: worked on making SimPM non-global
/// - 2018.04.14 JM: Moved flux solver to FV_solver

#ifndef SOLVER_EQN_BASE_H
#define SOLVER_EQN_BASE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "coord_sys/VectorOps.h"
#include "equations/eqns_base.h"

///
/// Base Class for Flux-based Finite Volume Solvers.
/// This defines the interface to the main equations solvers.
///
class FV_solver_base : virtual public eqns_base, virtual public BaseVectorOps {
public:
  FV_solver_base(
      const int,     ///< number of variables in state vector.
      const int,     ///< number of space dimensions in grid.
      const double,  ///< CFL number
      const double,  ///< gas eos gamma.
      const double,  ///< Artificial Viscosity Parameter etav.
      const int      ///< Number of tracer variables.
  );                 ///< Constructor.

  ~FV_solver_base();  ///< Destructor.

  /// \brief Calculates GLM wave-speed -- only used by GLM-MHD equations.
  virtual void Set_GLM_Speeds(
      const double,  ///< Current timestep value.
      const double,  ///< grid cell size dx.
      const double   ///< GLM damping coefficient c_r
  )
  {
    return;
  }

  /// calculate Powell and GLM source terms for multi-D MHD
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
  )
  {
    return 0;
  }

  /// \brief sets current timestep value in the solver class.
  virtual void Setdt(const double delt  ///< Current timestep value.
  )
  {
    FV_dt = delt;
    //  cout <<"dt="<<FV_dt<<endl;
    return;
  }

  /// \brief sets value of gamma, the ideal gas EOS parameter.
  virtual void SetEOS(const double g  ///< gas EOS gamma.
  )
  {
    eq_gamma = g;
    return;
  }

  ///
  /// Set the H-correction eta value in a cell at the positive
  /// boundary along a given axis, based on the previously
  /// calculated left and right edge states at the boundary.
  ///
  void set_Hcorrection(
      cell &,            ///< cell to assign eta value to.
      const axes,        ///< axis we are looking along
      const pion_flt *,  ///< left state (from current cell).
      const pion_flt *,  ///< right state (from next cell).
      const double       ///< EOS gamma.
  );

  ///
  /// Calculate Flux between a left and right state.
  ///
  int InterCellFlux(
      class SimParams &,          ///< simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid
      class cell &,               ///< Left state cell pointer
      class cell &,               ///< Right state cell pointer
      pion_flt *,                 ///< Left Primitive State Vector.
      pion_flt *,                 ///< Right Primitive State Vector.
      pion_flt *,                 ///< Flux Vector. (written to).
      const double,               ///< gas EOS gamma.
      const double                ///< Cell size dx.
  );

  ///
  /// Calculates Flux based on a left and right state vector (prim).
  ///
  virtual int inviscid_flux(
      class SimParams &,      ///< simulation parameters
      class GridBaseClass *,  ///<  pointer to grid
      const double,           ///< cell-size dx (for LF method)
      class cell &,           ///< Left state cell pointer
      class cell &,           ///< Right state cell pointer
      const pion_flt *,       ///< Left Primitive state vector.
      const pion_flt *,       ///< Right Primitive state vector.
      pion_flt *,             ///< Resultant Flux state vector.
      pion_flt *,             ///< State vector at interface.
      const int,              ///< Which Riemann solver
      const double            ///< Gas constant gamma.
      ) = 0;

  /// Adds the contribution from flux in the current direction to dU.
  virtual int dU_Cell(
      class GridBaseClass *grid,
      cell &,            ///< Current cell.
      const axes,        ///< Which axis we are looking along.
      const pion_flt *,  ///< Negative direction flux.
      const pion_flt *,  ///< Positive direction flux.
      const pion_flt *,  ///< slope vector for cell c.
      const int,         ///< spatial order of accuracy.
      const double,      ///< cell length dx.
      const double       ///< cell TimeStep, dt.
      ) = 0;

  /// \brief General Finite volume scheme for updating a cell's
  /// primitive state vector, for homogeneous equations.
  ///
  virtual int CellAdvanceTime(
      class cell &,      ///< cell to update.
      const pion_flt *,  ///< Initial Primitive State Vector.
      pion_flt *,        ///< Update vector dU
      pion_flt *,    ///< Final Primitive state vector (can be same as initial
                     ///< vec.).
      pion_flt *,    ///< Tracks change of energy if I have to correct for
                     ///< negative pressure
      const double,  ///< gas EOS gamma.
      const double,  ///< Min Temperature allowed on grid.
      const double   ///< Cell timestep dt.
      ) = 0;

  /// \brief Given a cell, calculate the hydrodynamic timestep.
  virtual double CellTimeStep(
      const cell &,  ///< pointer to cell
      const double,  ///< gas EOS gamma.
      const double   ///< Cell size dx.
      ) = 0;


  ///
  /// This function is for integrations which need to
  /// add a multi-dimensional correction/addition to the update
  /// vector.  Examples are the Field-CD method which processes the
  /// update vector to ensure that the B-field remains divergence free
  /// (not really working), and the Pdiv(V) term in the internal--
  /// energy solver for the Euler equations.
  /// Base implementation does nothing.
  ///
  virtual int PostProcess_dU(
      const double,          ///< current timestep, dt.
      const int,             ///< TIMESTEP_FIRST_PART orTIMESTEP_FULL
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *  ///< pointer to computational grid.
  )
  {
    return 0;
  }

protected:
  const int FV_gndim;  ///< number of spatial directions in grid.
  /// Courant-Friedrichs-Levy parameter (<1 for stability).
  const double FV_cfl;
  double FV_dt;  ///< Timestep
  /// coefficient of (artificial) viscosity for velocity field.
  const double FV_etav;
  /// coefficient of (artificial) viscosity for magnetic field.
  const double FV_etaB;
  /// max. value of eta used for H-correction, if required.
  double HC_etamax;
  /// Number of passive tracer variables.
  const int FV_ntr;
  /// Pointer to array of indices of tracer variables in the state vector.
  int *eqTR;

  ///
  /// This calculates the Lax-Friedrichs flux across an interface, but the
  /// implementation has to be in the general solver classes because LF flux
  /// requires dx, dt, and ndim, which the equations/riemann solvers don't
  /// know about and don't care about. This function should work for most
  /// equations, but glm-mhd will need some extra work.
  ///
  virtual int get_LaxFriedrichs_flux(
      const pion_flt *,  ///< Left  Primitive var. state vector.
      const pion_flt *,  ///< Right Primitive var. state vector.
      pion_flt *,        ///< Resulting Flux vector.
      const double,      ///< cell size dx
      const double       ///< gamma
  );

  ///
  /// Falle, Komissarov & Joarder (1998,MNRAS,297,265) Artificial
  /// Viscosity Calculation (one-dimensional).
  ///
  virtual int AVFalle(
      const pion_flt *,  ///< Left Primitive state vector.
      const pion_flt *,  ///< Right Primitive state vector.
      const pion_flt *,  ///< Resolved (P*) state vector.
      pion_flt *,        ///< Pointer to associated Flux Vector.
      const double,      ///< Artificial Viscosity parameter, etav.
      const double       ///< gamma
      ) = 0;

  ///
  /// Function to calculate viscosity-related quantities before the
  /// flux calculation.  So far this is only for the H-correction.
  ///
  void pre_calc_viscous_terms(
      class GridBaseClass *,  ///< pointer to computational grid.
      const cell &,           ///< left-of-interface cell
      const cell &,           ///< right-of-interface cell
      const int               ///< What kind of AV?
  );

  ///
  /// Function to calculate viscous modifications to the flux (which
  /// has already been calculated).  This does nothing for the
  /// H-correction, but everything for FKJ98 and LAPIDUS viscosities.
  ///
  void post_calc_viscous_terms(
      const cell &,      ///< left-of-interface cell
      const cell &,      ///< right-of-interface cell
      const pion_flt *,  ///< left state Prim.vec. (Pl)
      const pion_flt *,  ///< right state Prim.vec. (Pr)
      const pion_flt *,  ///< interface state Prim.vec. (P*)
      pion_flt *,        ///< flux vector.
      const int          ///< what kind of AV?
  );

  ///
  /// This calculates the maximum of the H-shaped list of interface
  /// eta values for the H-correction of Sanders et al
  /// (1998,JCP,14,511).  Uses Eq. 16 and Fig. 9 from paper.
  ///
  double select_Hcorr_eta(
      const cell &,          ///< cell to left of interface
      const cell &,          ///< cell to right of interface
      class GridBaseClass *  ///< pointer to computational grid.
  );

  ///
  /// Calculate tracer flux based on a flux vector and left and right
  /// states.  If mass flux is positive, then tracer flux is advected
  /// from the left state; if negative then from the right state.
  /// Flux must be already calculated for the physical variables!
  ///
  void set_interface_tracer_flux(
      const pion_flt *,  ///< input left state.
      const pion_flt *,  ///< input right state
      pion_flt *         ///< Flux vector.
  );

  ///	Vector for corrector values (used to modify flux according to sCMA)
  std::vector<double> corrector;
};

#endif  // SOLVER_EQN_BASE_H
