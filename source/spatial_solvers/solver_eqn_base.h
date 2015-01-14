/** \file solver_eqn_base.h
 * 
 * \brief Class declarations for various solvers implemented in my FV grid-code.
 * \author Jonathan Mackey
 * 
 * Modifications:
 *  - 2007-07-11 Still writing it and working out class hierarchy.
 *  - 2007-07-12 Got the class heirarchy working (I think...).
 *  - 2007-07-13 New Class structure implemented and working
 *  - 2007-07-16 Reverted to less complicated class hierarchy. works well.
 *  - 2007-07-23 Started to add passive tracer variables
 *  - 2007-07-24 Added passive tracer variable support.
 *  - 2007-08-07 cylindrical coordinates support working.
 *  - 2007-08-08 cyl.coords for iMHD and glmMHD.
 *  - 2007-12-11 Field-CD method 'working' in 2D cartesian coords.
 *  - 2009-10-20 New structure built on equations and flux classes...
 * */
///
///  - 2010.09.30 JM: Worked on Lapidus AV (added set_div_v() function for cells)
///
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
///
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

#ifndef SOLVER_EQN_BASE_H
#define SOLVER_EQN_BASE_H



#include "global.h"
#include "coord_sys/VectorOps.h"
#include "equations/eqns_base.h"
#include "flux_calc/flux_base.h"

///
/// Base Class for Flux-based Finite Volume Solvers.
/// This defines the interface to the main equations solvers.
///
class FV_solver_base : virtual public flux_solver_base, virtual public BaseVectorOps
{
 public :
  FV_solver_base(const int, ///< number of variables in state vector.
		 const int, ///< number of space dimensions in grid.
		 const double, ///< CFL number
		 const double, ///< dx, cell size.
		 const double, ///< gas eos gamma.
		 const double, ///< Artificial Viscosity Parameter etav.
		 const int     ///< Number of tracer variables.
		 ); ///< Constructor.
  ~FV_solver_base(); ///< Destructor.
  /** \brief Calculates GLM wave-speed -- only used by GLM-MHD equations. */
  virtual void GotTimestep(const double ///< Current timestep value.
			   ) {return;}
  /** \brief sets current timestep value in the solver class. */
  virtual void Setdt(const double delt ///< Current timestep value.
		     ) 
   {
     FV_dt =delt;
   //  cout <<"dt="<<FV_dt<<endl;
     return;
   }
  /** \brief sets value of gamma, the ideal gas EOS parameter. */
  virtual void SetEOS(const double g ///< gas EOS gamma.
		      ) {eq_gamma = g; return;}

  ///
  /// Set the velocity divergence in a cell for viscosity calculation.
  ///
  void set_div_v(
        class cell *,
        class GridBaseClass *  ///< pointer to computational grid.
        );

  ///
  /// Set the H-correction eta value in a cell at the positive
  /// boundary along a given axis, based on the previously
  /// calculated left and right edge states at the boundary.
  ///
  void set_Hcorrection(cell *,         ///< cell to assign eta value to.
		       const axes,     ///< axis we are looking along
		       const double *, ///< left state (from current cell).
		       const double *, ///< right state (from next cell).
		       const double    ///< EOS gamma.
		       );

#ifdef THERMAL_CONDUCTION
  ///
  /// Calculate the rate of heating and cooling in every cell due to thermal
  /// conduction.  This sweeps through the data Ndim+2 times (once to get
  /// the temperature (stored in dU[RHO] temporarily, so this function must
  /// be called at the beginning of a dynamics update where I am not 
  /// overwriting any information!!.  This is probably not safe for long-term).
  /// It uses c->Ph[] as the source data for calculating T, n(H), etc.
  /// It is virtual because with MHD the flux becomes dependent on the magnetic
  /// field direction so I may want to add that later.
  ///
  virtual int set_thermal_conduction_Edot();
#endif // THERMAL CONDUCTION


  ///
  /// Calculate Flux between a left and right state.
  ///
  int InterCellFlux(
        class GridBaseClass *grid,
        const cell *, ///< Left state cell pointer
        const cell *, ///< Right state cell pointer
        double *, ///< Left Primitive State Vector.
        double *, ///< Right Primitive State Vector.
        double *, ///< Flux Vector. (written to).
        const int,      ///< Solve Type (0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS)
        const int,      ///< Viscosity Flag (0=none,1=Falle's,2=Lapidus(broken),etc.)
        const double, ///< gas EOS gamma.
        const double  ///< Cell size dx.
        );

  /// Adds the contribution from flux in the current direction to dU.
  virtual int dU_Cell(
        class GridBaseClass *grid,
        cell *, ///< Current cell.
        const axes, ///< Which axis we are looking along.
        const double *, ///< Negative direction flux.
        const double *, ///< Positive direction flux.
        const double *, ///< slope vector for cell c.
        const int,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        )=0;

  /** \brief General Finite volume scheme for updating a cell's
   * primitive state vector, for homogeneous equations.
   * */
  virtual int CellAdvanceTime(class cell *, ///< cell to update.
			      const double *, ///< Initial Primitive State Vector.
			      double *, ///< Update vector dU
			      double *, ///< Final Primitive state vector (can be same as initial vec.).
			      double *,  ///< Tracks change of energy if I have to correct for negative pressure
			      const double, ///< gas EOS gamma.
			      const double  ///< Cell timestep dt.
			      )=0;
  /** \brief Given a cell, calculate the hydrodynamic timestep. */
  virtual double CellTimeStep(const cell *, ///< pointer to cell
			      const double, ///< gas EOS gamma.
			      const double  ///< Cell size dx.
			      ) =0;
  
  ///
  /// Pre-process the grid before calc_dU() in this function.
  /// Multi-dimensional calculations to be performed on every cell
  /// before the flux calculations. e.g. it calculates div(v) for
  /// Lapidus viscosity and the Ndim 'eta' parameters for the
  /// H-correction viscosity.
  ///
  virtual int preprocess_data(
        const int, ///< Space order of acc for this call.
        const int, ///< Time order of acc for this call.
        class GridBaseClass *  ///< pointer to computational grid.
        );
  

   ///
   /// Runs through the grid in all directions and calculates 'eta'
   /// for each grid point in each direction.  These values get stored
   /// in the cell data.
   /// 
   int calc_Hcorrection(
        const int, ///< Space order of acc for this call.
        const int, ///< Time order of acc for this call.
        class GridBaseClass *  ///< pointer to computational grid.
        );


  ///
  /// This function is a place-holder for integrations which need to 
  /// add a multi-dimensional correction/addition to the update
  /// vector.  Examples are the Field-CD method which processes the
  /// update vector to ensure that the B-field remains divergence free
  /// (not really working), and the Pdiv(V) term in the internal--
  /// energy solver for the Euler equations.
  /// For most solvers it just returns immediately.
  ///
  virtual int PostProcess_dU(
        const int, ///< Space order of acc for this call.
        const int, ///< Time order of acc for this call.
        class GridBaseClass *  ///< pointer to computational grid.
        ) {return 0;}

 protected:
  const int FV_gndim;   ///< number of spatial directions in grid.
  //  const int FV_ntr;     ///< Number of tracer Variables.
  const double FV_cfl;  ///< Courant-Friedrichs-Levy parameter (<1 for stability).
  // INHERITED   const double FV_etav; ///< Artificial Viscosity parameter.
  // INHERITED  double FV_gamma;      ///< EOS gamma
  double FV_dx;         ///< Grid spacing delta-x
  double FV_dt;         ///< Timestep
  ///
  /// This calculates the Lax-Friedrichs flux across an interface, but the 
  /// implementation has to be in the general solver classes because LF flux
  /// requires dx, dt, and ndim, which the equations/riemann solvers don't know 
  /// about and don't care about.
  /// This function should work for most equations, but glm-mhd will need some extra work.
  ///
  virtual int get_LaxFriedrichs_flux(const double *, ///< Left  Primitive var. state vector.
				     const double *, ///< Right Primitive var. state vector.
				     double *,       ///< Resulting Flux vector.
				     const double    ///< gamma
				     );

  ///
  /// Function to calculate viscosity-related quantities before the
  /// flux calculation.  So far this is only for the H-correction.
  ///
  void pre_calc_viscous_terms(
        class GridBaseClass *,  ///< pointer to computational grid.
        const cell *, ///< left-of-interface cell
        const cell *, ///< right-of-interface cell
        const int     ///< What kind of AV?
        );

  ///
  /// Function to calculate viscous modifications to the flux (which
  /// has already been calculated).  This does nothing for the 
  /// H-correction, but everything for FKJ98 and LAPIDUS viscosities.
  ///
  void post_calc_viscous_terms(const cell *, ///< left-of-interface cell
			       const cell *, ///< right-of-interface cell
			       const double *, ///< left state Prim.vec. (Pl)
			       const double *, ///< right state Prim.vec. (Pr)
			       const double *, ///< interface state Prim.vec. (P*)
			       double *,     ///< flux vector.
			       const int     ///< what kind of AV?
			       );


  ///
  /// This calculates the maximum of the H-shaped list of interface
  /// eta values for the H-correction of Sanders et al
  /// (1998,JCP,14,511).  Uses Eq. 16 and Fig. 9 from paper.
  ///
  double select_Hcorr_eta(
        const cell *, ///< cell to left of interface
        const cell *, ///< cell to right of interface
        class GridBaseClass *  ///< pointer to computational grid.
        );
};



#ifdef NOT_NEW_SOLVER_STRUCTURE
//-----------------------------------------------------------
// This is old code, ifdeffed out so it should never compile.
// May be useful for updating all the solvers later...
//-----------------------------------------------------------

#include "../global.h"
#include "equations.h"
#include "../coord_sys/VectorOps.h"
#include "riemann.h"
#include "riemannMHD.h"
#include "microphysics.h"


/** \brief Base Class for Flux-based Finite Volume Solvers.
 * This is a pure virtual base class. */
class BaseFVSolver : virtual public BaseEqn, virtual public BaseVectorOps
{
  public :
   ~BaseFVSolver() {} ///< Destructor.
   /** \brief Calculates GLM wave-speed -- only used by GLM-MHD equations. */
   virtual void GotTimestep(const double delt ///< Current timestep value.
		      ) =0;
   /** \brief sets current timestep value in the solver class. */
   virtual void Setdt(const double delt ///< Current timestep value.
		      ) =0;
   /** \brief sets value of gamma, the ideal gas EOS parameter. */
   virtual void SetEOS(const double g ///< gas EOS gamma.
		       ) =0;
   /** \brief Sets which axes we are looking along. */
   virtual void SetDir(const enum axes d) =0;
   /** \brief Sets Average State Vector for anything that needs it (e.g. Riemann Solver)*/
   virtual void SetAvgState(const double *mv, ///< Average primitive state vector.
			    const double g   ///< EOS gamma.
			    ) =0;
   /** \brief Returns which axes we are looking along. */
   virtual enum axes GetDir() =0;
   /** \brief Calculate Flux between a left and right state.
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar Vector. (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     ) =0;
   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		       const axes, ///< Which axis we are looking along.
		       const double *, ///< Negative direction flux.
		       const double *, ///< Positive direction flux.
		       const double *, ///< slope vector for cell c.
		       const int,      ///< spatial order of accuracy.
		       const double, ///< cell length dx.
		       const double  ///< cell TimeStep, dt.
		       )=0;
   /** \brief General Finite volume scheme for updating a cell's
    * primitive state vector, for homogeneous equations.
    * */
   virtual int CellAdvanceTime(const double *, ///< Initial Primitive State Vector.
			       double *, ///< Update vector dU
			       double *, ///< Final Primitive state vector (can be same as initial vec.).
			       double *,  ///< Tracks change of energy if I have to correct for negative pressure
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell timestep dt.
			       )=0;
   /** \brief Given a cell, calculate the hydrodynamic timestep. */
   virtual double CellTimeStep(const cell *, ///< pointer to cell
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell size dx.
			       ) =0;
   /** \brief This function only does anything for the Field-CD method which 
    * processes the update vector to ensure that the B-field remains divergence free.
    * For all other solvers it just returns immediately.
    * */
   virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
			      const int  ///< Time order of acc for this call.
			      )=0;
};

/** \brief Basic Lax-Friedrichs Solver.
 * 
 * The InterCellFlux() function calculates the LF Flux according to the
 * scheme outlined in Toro (1999) p.184-185. (eqn.5.77)
 * 
 * The CellAdvanceTime() function is a general FV time update for a 
 * system of homogeneous conservation equations (i.e. without source
 * terms), so it should be general for many Solvers.
 *
 * \section References
 * Toro, E.F., 1999, Riemann Solvers and Numerical Methods for Fluid Dynamics,
 * (Text Book), Springer, ISBN: 3-540-65966-8
 * */
class LF_FVSolver : public BaseFVSolver
{
  public:
   LF_FVSolver(const int, ///< number of variables in state vector.
	       const int, ///< number of space dimensions in grid.
	       const double, ///< dx, cell size.
	       const double ///< gas eos gamma.
	       ); ///< Constructor.
   ~LF_FVSolver(); ///< Destructor.
   /** \brief This function only does anything for the Field-CD method which 
    * processes the update vector to ensure that the B-field remains divergence free.
    * For all other solvers it just returns immediately.
    * */
   virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
			      const int  ///< Time order of acc for this call.
			      );
   /** \brief Lax Friedrichs Flux between a left and right state.
    * Lax-Friedrichs Scheme, Toro (1999) p.184-185. (eqn.5.77)
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar Vector. (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     );
   /** \brief General Finite volume scheme for updating a cell's
    * primitive state vector, for homogeneous equations.
    * 
    * Notes:\n
    *  - Currently gamma and dt are not passed in, but are assumed to 
    * be set correctly as class variables.
    *  - Source Terms must be added into a derived class.
    * */
   virtual int CellAdvanceTime(const double *, ///< Initial State Vector.
			       double *,       ///< Update vector dU
			       double *,       ///< Final state vector (can be same as initial vec.).
			       double *,       ///< Tracks change of energy if I have to correct for negative pressure
			       const double,   ///< gas EOS gamma.
			       const double    ///< Cell timestep dt.
			       );
   /** \brief sets value of gamma, the ideal gas EOS parameter. */
   virtual void SetEOS(const double g ///< gas EOS gamma.
		       ) {gamma = g; return;}
   /** \brief Does nothing until re-implemented. */
   virtual void GotTimestep(const double delt ///< Current timestep value.
		      ) {return;}
   /** \brief sets current timestep. */
   virtual void Setdt(const double delt ///< Current timestep value.
		      ) {dt =delt; return;}
   /** \brief Sets which axes we are looking along. */
   virtual void SetDir(const enum axes d) {SetDirection(d);}
   /** \brief Returns which axes we are looking along. */
   virtual enum axes GetDir() {return(GetDirection());}
   /** \brief Sets mean values for pressure,density,velocity,B-field, etc.,
    * for anything that needs it (e.g. Riemann Solver)*/
   virtual void SetAvgState(const double *, ///< Average primitive state vector.
			    const double   ///< EOS gamma.
			    );
   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		       const axes, ///< Which axis we are looking along.
		       const double *, ///< Negative direction flux.
		       const double *, ///< Positive direction flux.
		       const double *, ///< slope vector for cell c.
		       const int,      ///< spatial order of accuracy.
		       const double, ///< cell length dx.
		       const double  ///< cell TimeStep, dt.
		       )=0;
  protected:
   const int nvar;
   const int gndim;
   double gamma;
   double dx;
   double dt;
#ifdef USE_MM
#else
   double *u1, *u2, *f1, *f2;
#endif
};

/** \brief Abstract class containing some generic functions used by
 * Riemann-solver based solvers.
 * */
class RS_FVSolver : virtual public LF_FVSolver
{
  public:
   /** \brief constructor. */
   RS_FVSolver(const int, ///< number of variables in state vector.
	       const int, ///< number of space dimensions in grid.
	       const double, ///< CFL number
	       const double, ///< dx, cell size.
	       const double, ///< gas eos gamma.
	       double *      ///< State vector of mean values for simulation.
	       );
   ~RS_FVSolver(); ///< destructor.
   /** \brief Simple Riemann Solve Calculation of Flux between left and right state.
    * 
    * The axis we are looking along (normal direction) should have already been set.
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar Vector. (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     );
   /** \brief Sets which axes we are looking along. */
   virtual void SetDir(const enum axes);
   /** \brief Sets Average State Vector for Riemann Solver)*/
   virtual void SetAvgState(const double *, ///< Average primitive state vector.
			    const double   ///< EOS gamma.
			    );
  protected:
   class BaseRiemann *rs; ///< pointer to the riemann solver class.
   double *meanvec;  ///< vector of typical values for primitive variables in the simulation.
};

/** \brief Riemann solver class with Artificial Viscosity, and the ability to handle
 * tracer variables.  Abstract class with generic functions.
 * */
class RS_AV_TR_FVSolver : virtual public LF_FVSolver, virtual public RS_FVSolver
{
  public:
   /** \brief constructor. */
   RS_AV_TR_FVSolver(const int, ///< number of variables in state vector.
		     const int, ///< number of space dimensions in grid.
		     const double, ///< CFL number
		     const double, ///< dx, cell size.
		     const double, ///< gas eos gamma.
		     double *,     ///< State vector of mean values for simulation.
		     const double, ///< Artificial Viscosity Parameter etav.
		     const int     ///< Number of tracer variables.
		     );
   ~RS_AV_TR_FVSolver(); ///< destructor.
   /** \brief Riemann Solve Calculation of Flux between left and right state, with Artificial Viscosity.
    * 
    * The axis we are looking along (normal direction) should have already been set.
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar Vector (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     );
  protected :
   const double etav; ///< Artificial Viscosity parameter.
   int *eqTR; ///< Pointer to array of indices of tracer variables in the state vector.
   const int ntr; ///< Number of tracer Variables.
};
   

/** \brief The LF Method for the Euler Equations. */
class Euler_LF : virtual public LF_FVSolver, public EulerEqn, virtual public VectorOps_Cart
{
  public:
   Euler_LF(const int, ///< number of variables in state vector.
	   const int, ///< number of space dimensions in grid.
	   const double, ///< CFL number
	   const double, ///< dx, cell size.
	   const double ///< gas eos gamma.
	   );
   ~Euler_LF();
   /** \brief Calculate the shortest allowable timestep for this cell, 
    * based on the Courant condition. */
   virtual double CellTimeStep(const cell *, ///< pointer to cell
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell size dx.
			       );
   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
  protected:
   const int eqndim; ///< number of equations to solve.
   const double cfl; ///< Courant-Friedrichs-Levy parameter (<1 for stability).
};

/** \brief The Riemann Solver (Godunov) Method for the Euler Equations. */
class Euler_RS : virtual public Euler_LF, virtual public RS_FVSolver
{
  public:
   Euler_RS(const int, ///< number of variables in state vector.
	   const int, ///< number of space dimensions in grid.
	   const double, ///< CFL number
	   const double, ///< dx, cell size.
	   const double, ///< gas eos gamma.
	   double *      ///< State vector of mean values for simulation.
	   );
   ~Euler_RS();
};

/** \brief The LF Method for the Ideal MHD Equations (no divB cleaning). */
class IdealMHD_LF : virtual public LF_FVSolver, public IdealMHDEqn, virtual public VectorOps_Cart
{
  public:
   IdealMHD_LF(const int, ///< number of variables in state vector.
	       const int, ///< number of space dimensions in grid.
	       const double, ///< CFL number
	       const double, ///< dx, cell size.
	       const double ///< gas eos gamma.
	       );
   ~IdealMHD_LF();
   virtual double CellTimeStep(const cell *, ///< pointer to cell
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell size dx.
			       );
   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
  protected:
   const int eqndim;
   const double cfl;
};

/** \brief  The Riemann Solver (Godunov) Method for the Ideal MHD Equations.
 * 
 * This improves on the LF method by using a linear Roe-type riemann solver
 * to solve the ideal MHD equations (with no divB cleaning).  It's not 
 * recommended for full use as it has no artificial diffusion which can lead
 * to streaking and other numerical artefacts.  It is good for testing though,
 * as the fluxes are just from the riemann solver, so if there is a problem 
 * it is easier to find.
 */
class IdealMHD_RS : virtual public IdealMHD_LF, virtual public RS_FVSolver
{
  public:
   IdealMHD_RS(const int, ///< number of variables in state vector.
	       const int, ///< number of space dimensions in grid.
	       const double, ///< CFL number
	       const double, ///< dx, cell size.
	       const double, ///< gas eos gamma.
	       double *      ///< State vector of mean values for simulation.
	       );
   ~IdealMHD_RS();
};


/** \brief Solver Based on Gabor Toth's (2000) Field-CD method.
 * 
 * This method is the same as the usual finite volume method, but then 
 * discards the update vector of the magnetic field for a new calculation.
 * The Electric field is calculated from the cross product of the velocity
 * and the Magnetic field (negative of this).  The the curl of the Electric
 * field gives the negative of the Magnetic field update, so we multiply the 
 * curl by minus dt and this is the new update.  This guarantees that the 
 * magnetic field remains divergence free, but we no longer conserve magnetic
 * energy, and we have changed the magnetic energy in the cell without a 
 * corresponding change in the total energy.  In theory they are independent
 * variables so this is not a problem, but it may well be a problem.
 */
class FieldCD_MHD_RS : virtual public IdealMHD_RS
{
  public:
   FieldCD_MHD_RS(const int, ///< number of variables in state vector.
		  const int, ///< number of space dimensions in grid.
		  const double, ///< CFL number
		  const double, ///< dx, cell size.
		  const double, ///< gas eos gamma.
		  double *      ///< State vector of mean values for simulation.
		  );
   ~FieldCD_MHD_RS();
   /** \brief This function only does anything for the Field-CD method which 
    * processes the update vector to ensure that the B-field remains divergence free.
    * For all other solvers it just returns immediately.
    * */
   virtual int PostProcess_dU(const int, ///< Space order of acc for this call.
			      const int  ///< Time order of acc for this call.
			      );
};


/** \brief MHD solver using Dedner's GLM method for divergence cleaning.
 * 
 * This is the one to use for multi-dimensional simulations with magnetic
 * fields, as it advects the monopoles to the boundary as fast as possible,
 * as well as smoothing them out as much as possible.
 * */
class glmMHD_RS : virtual public RS_FVSolver,
  public glmMHDEqn,
  virtual public VectorOps_Cart
{
  public:
   glmMHD_RS(const int, ///< number of variables in state vector.
	      const int, ///< number of space dimensions in grid.
	      const double, ///< CFL number
	      const double, ///< dx, cell size.
	      const double, ///< gas eos gamma.
	      double *      ///< State vector of mean values for simulation.
	      );
   ~glmMHD_RS();
   virtual double CellTimeStep(const cell *, ///< pointer to cell
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell size dx.
			       );
   /** \brief Calculate TimeStep for a cell based on the max. wave speed, and grad(Pg)
    * The pressure gradient is passed in as a parameter.  Sets GLM speed.  Note that the
    * solver variable delt is not set by this function; call Setdt(delt) instead.
    * */
   virtual void GotTimestep(const double ///< timestep dt.
			    );
   /** \brief Updates previous FV cell-update by adding in source term from GLM method.
    * 
    * Notes:\n
    *  - Currently gamma and dt are not passed in, but are assumed to 
    * be set correctly as class variables.
    *  - The GLM source terms is calculated from glmMHDEqn::GLMSource()
    * */
   virtual int CellAdvanceTime(const double *, ///< Initial State Vector.
			       double *, ///< Update vector dU
			       double *, ///< Final state vector (can be same as initial vec.).
			       double *, ///< Tracks change of energy if I have to correct for negative pressure
			       const double, ///< gas EOS gamma.
			       const double  ///< Cell timestep dt.
			       );
   /** \brief Riemann Solve Calculation of Flux between left and right state, including GLM solve.
    * The axis we are looking along (normal direction) should have already been set.
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar vector (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     );

   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
  protected:
   const int eqndim; ///< Number of elements in vector quantities: velocity, B-field.
   const double cfl; ///< Courant-Friedrich-Levy factor (<1 for stable scheme).
};


/** \brief Solver for Euler equations with AV and tracers. */
class Euler_RS_AV_TR : virtual public Euler_RS, virtual public RS_AV_TR_FVSolver
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   Euler_RS_AV_TR(const int, ///< number of variables in state vector.
		 const int, ///< number of space dimensions in grid.
		 const double, ///< CFL number
		 const double, ///< dx, cell size.
		 const double, ///< gas eos gamma.
		 double *,     ///< State vector of mean values for simulation.
		 const double, ///< Artificial Viscosity Parameter etav.
		 const int     ///< Number of tracer variables.
		 );
   ~Euler_RS_AV_TR();
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int PtoU(const double *, ///< pointer to Primitive variables.
		    double *,       ///< pointer to conserved variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int UtoP(const double *, ///< pointer to conserved variables.
		    double *, ///< pointer to Primitive variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief convert direct from primitive variables to flux.
    * Creates conserved variable array as an intermediate step, 
    * and then calls PUtoFlux(). */
   //   virtual int PtoFlux(const double *, ///< pointer to Primitive variables.
   //		       double *,       ///< Pointer to flux variable.
   //	       const double    ///< Gas constant gamma.
   //	       );   
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int UtoFlux(const double*, ///< Pointer to conserved variables state vector.
		       double*,       ///< Pointer to flux variable state vector.
		       const double   ///< Gas constant gamma.
		       );
};

/** \brief Solver for Euler equations in axial symmetry with AV and tracers. */
class cyl_Euler_RS_AV_TR 
  : virtual public Euler_RS_AV_TR,
  virtual public VectorOps_Cyl
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   cyl_Euler_RS_AV_TR(const int, ///< number of variables in state vector.
		      const int, ///< number of space dimensions in grid.
		      const double, ///< CFL number
		      const double, ///< dx, cell size.
		      const double, ///< gas eos gamma.
		      double *,     ///< State vector of mean values for simulation.
		      const double, ///< Artificial Viscosity Parameter etav.
		      const int     ///< Number of tracer variables.
		      );
   ~cyl_Euler_RS_AV_TR();
   /** \brief Adds the contribution from flux in the current direction to dU. */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
};

/** \brief Solver for Ideal MHD equations with AV and tracers. */
class IdealMHD_RS_AV_TR : virtual public IdealMHD_RS, virtual public RS_AV_TR_FVSolver
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   IdealMHD_RS_AV_TR(const int, ///< number of variables in state vector.
		 const int, ///< number of space dimensions in grid.
		 const double, ///< CFL number
		 const double, ///< dx, cell size.
		 const double, ///< gas eos gamma.
		 double *,     ///< State vector of mean values for simulation.
		 const double, ///< Artificial Viscosity Parameter etav.
		 const int     ///< Number of tracer variables.
		 );
   ~IdealMHD_RS_AV_TR();
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int PtoU(const double *, ///< pointer to Primitive variables.
		    double *,       ///< pointer to conserved variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int UtoP(const double *, ///< pointer to conserved variables.
		    double *, ///< pointer to Primitive variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief convert direct from primitive variables to flux.
    * Creates conserved variable array as an intermediate step, 
    * and then calls PUtoFlux(). */
   //   virtual int PtoFlux(const double *, ///< pointer to Primitive variables.
   //	       double *,       ///< Pointer to flux variable.
   //	       const double    ///< Gas constant gamma.
   //	       );   
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int UtoFlux(const double*, ///< Pointer to conserved variables state vector.
		       double*,       ///< Pointer to flux variable state vector.
		       const double   ///< Gas constant gamma.
		       );
};

/** \brief Solver for Ideal MHD equations in axial symmetry with AV and tracers. */
class cyl_IdealMHD_RS_AV_TR
  : virtual public IdealMHD_RS_AV_TR, virtual public VectorOps_Cyl
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   cyl_IdealMHD_RS_AV_TR(const int, ///< number of variables in state vector.
		      const int, ///< number of space dimensions in grid.
		      const double, ///< CFL number
		      const double, ///< dx, cell size.
		      const double, ///< gas eos gamma.
		      double *,     ///< State vector of mean values for simulation.
		      const double, ///< Artificial Viscosity Parameter etav.
		      const int     ///< Number of tracer variables.
		      );
   ~cyl_IdealMHD_RS_AV_TR();
   /** \brief Adds the contribution from flux in the current direction to dU.
    * 
    * Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
    * spatial accuracy.
    */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
};

/** \brief Solver for Ideal MHD equations with AV and tracers, and the 
 * Field-CD divergence cleaning method.
 * 
 * Note this doesn't really work well -- use GLM method instead.
 */
class FieldCD_MHD_RS_AV_TR 
  : virtual public FieldCD_MHD_RS,
  virtual public IdealMHD_RS_AV_TR
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   FieldCD_MHD_RS_AV_TR(const int, ///< number of variables in state vector.
			const int, ///< number of space dimensions in grid.
			const double, ///< CFL number
			const double, ///< dx, cell size.
			const double, ///< gas eos gamma.
			double *,     ///< State vector of mean values for simulation.
			const double, ///< Artificial Viscosity Parameter etav.
			const int     ///< Number of tracer variables.
			);
   ~FieldCD_MHD_RS_AV_TR(); ///< destructor.
};

/** \brief Solver for Ideal MHD equations with AV and tracers, using Dedner
 * et al(2003) method for divB cleaning (their mixed-GLM method).  See the 
 * extra page on code algorithms.
 */
class glmMHD_RS_AV_TR : virtual public glmMHD_RS, virtual public RS_AV_TR_FVSolver
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   glmMHD_RS_AV_TR(const int, ///< number of variables in state vector.
		 const int, ///< number of space dimensions in grid.
		 const double, ///< CFL number
		 const double, ///< dx, cell size.
		 const double, ///< gas eos gamma.
		 double *,     ///< State vector of mean values for simulation.
		 const double, ///< Artificial Viscosity Parameter etav.
		 const int     ///< Number of tracer variables.
		 );
   ~glmMHD_RS_AV_TR();
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int PtoU(const double *, ///< pointer to Primitive variables.
		    double *,       ///< pointer to conserved variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief adds conversion of tracer variables.
    * 
    * For purely passive tracers, the primitive variable is just a number,
    * such as the 'colour' of the gas, or where it started out.  The 
    * conserved variable is the mass density of this.
    * */
   virtual int UtoP(const double *, ///< pointer to conserved variables.
		    double *, ///< pointer to Primitive variables.
		    const double    ///< Gas constant gamma.
		    );
   /** \brief convert direct from primitive variables to flux.
    * Creates conserved variable array as an intermediate step, 
    * and then calls PUtoFlux(). */
   //   virtual int PtoFlux(const double *, ///< pointer to Primitive variables.
   //	       double *,       ///< Pointer to flux variable.
   //	       const double    ///< Gas constant gamma.
   //	       );   
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int PUtoFlux(const double *, ///< pointer to Primitive variables.
			const double *, ///< pointer to conserved variables.
			double *  ///< Pointer to flux variable.
			);
   /** \brief adds conversion of tracer variables.
    * 
    * The flux of a passive tracer is equal to the mass flux times 
    * the value of the primitive tracer variable.  I take the left
    * state tracer var. if the mass flux is to the right, and vice versa.
    * */
   virtual int UtoFlux(const double*, ///< Pointer to conserved variables state vector.
		       double*,       ///< Pointer to flux variable state vector.
		       const double   ///< Gas constant gamma.
		       );
   /** \brief Simple Riemann Solve Calculation of Flux between left and right state.
    * 
    * The axis we are looking along (normal direction) should have already been set.
    * 
    * This gets the flux from a riemann solver, modifies it with artificial
    * viscosity, and then calculates the tracer flux according to the sign
    * of the mass flux.
    * */
   virtual int InterCellFlux(double *, ///< Left Primitive State Vector.
			     double *, ///< Right Primitive State Vector.
			     double *, ///< Pstar Vector. (written to).
			     double *, ///< Flux Vector. (written to).
			     const double, ///< gas EOS gamma.
			     const double  ///< Cell size dx.
			     );
};

/** \brief Solver for mixed-GLM MHD equations with AV and tracers, in 
 * axial symmetry. */
class cyl_glmMHD_RS_AV_TR
  : virtual public glmMHD_RS_AV_TR, virtual public VectorOps_Cyl
{
  public:
   /** \brief sets indices for tracer variables in state vector.*/
   cyl_glmMHD_RS_AV_TR(const int, ///< number of variables in state vector.
		      const int, ///< number of space dimensions in grid.
		      const double, ///< CFL number
		      const double, ///< dx, cell size.
		      const double, ///< gas eos gamma.
		      double *,     ///< State vector of mean values for simulation.
		      const double, ///< Artificial Viscosity Parameter etav.
		      const int     ///< Number of tracer variables.
		      );
   ~cyl_glmMHD_RS_AV_TR();
   /** \brief Adds the contribution from flux in the current direction to dU.
    * 
    * Includes geometric source term (p^2+B^2/2)/R for 1st and 2nd order
    * spatial accuracy.
    */
   virtual int dU_Cell(cell *, ///< Current cell.
		  const axes, ///< Which axis we are looking along.
		  const double *, ///< Negative direction flux.
		  const double *, ///< Positive direction flux.
		  const double *, ///< slope vector for cell c.
		  const int,      ///< spatial order of accuracy.
		  const double, ///< cell length dx.
		  const double  ///< cell TimeStep, dt.
		  );
};

#endif // NOT_NEW_SOLVER_STRUCTURE

#endif //SOLVER_EQN_BASE_H

