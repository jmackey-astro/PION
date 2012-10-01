/** \file grid.h
 * 
 * \brief Declares grid parameter class, and grid methods classes.
 * 
 * \author Jonathan Mackey
 * 
 * IntegratorBaseFV is an abstract base class for finite volume grids.
 * 
 * IntUniformFV is a 1st/2nd order accurate finite volume solver, modelled on 
 * the solver presented in Falle, Komissarov, \& Joarder (1998) MNRAS, 297, 265.
 * 
 * Modifications :\n
 *  - 2007-07-12 Started to move to new structure.
 *  - 2007-07-13 New Class structure implemented.
 *  - 2007-07-24 Added tracer variables.
 *  - 2007-10-11 Updated Documentation.
 *  - 2007-10-26 New data i/o class.
 *
 *  - 2010-04-21 JM: removed parallel output_data() function b/c not needed anymore.
 *
 * */
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///
/// - 2010.10.04 JM: Moved the field-loop magnetic pressure output to
///     an ifdeffed function.  Added a new function to calculate the
///     1D blast wave radius in spherical coordinates (also ifdeffed).
///
/// - 2010.11.15 JM: Modified update_dynamics() so it has pre- and
///   post-processing functions before and after calc_dU().
///   Added routine to calculate the H-correction eta values for
///   pre-processing.
///
/// - 2011.01.03 JM: Moved preprocess_data() and calc_Hcorrection from
///   gridMethods.cc to the base FV solver.
/// - 2011.03.21 JM: moved cell setup_extra_data() to its own function, to save 
///    on copy-paste for the parallel version.
///    Rewrote setup_raytracing() and atomised update_microphysics() so there are
///    now a few new functions to deal with the different kinds of radiative
///    transfer we want to do.
///
/// - 2011.04.06 JM: Added thermal-conduction timestep limiting and flux calculation.
/// - 2011.10.22 JM: Added lists of heating and ionising sources as class data.
/// - 2012.01.16 JM: Added setup_evolving_RT_sources() and
///    update_evolving_RT_sources() for stellar evolution models.
/// - 2012.08.16 JM: Added functions for new 2nd order time update.

#ifndef GRID_H
#define GRID_H

#include "global.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "grid/uniform_grid.h"
#include "dataIO/dataio.h"


/** \brief The simplest finite volume grid -- a uniform grid with cubic cells
 * in the chosen coordinates.
 * 
 * This can integrate any system of equations if given the right solver class.
 * It can solve the equations in 1st or 2nd order accuracy in space and time.
 * */
class IntUniformFV : public IntegratorBaseFV
{
  public:
   IntUniformFV();  ///< Simple constructor, initialises value.
   ~IntUniformFV(); ///< Deletes any dynamic memory, if not already done.
   /** \brief initialisation.
    * This function calls a sequence of other functions to set up the grid
    * and populate it with the initial conditions, and give it the appropriate
    * boundary conditions.  It gets the simulation ready to start, and checks 
    * that everything is ready to start before returning.
    * \retval 0 success
    * \retval 1 failure
    * */
   virtual int Init(string, ///< Name of input file.
		    int,    ///< Type of File (1=ASCII, 2=FITS, 3=fitstable, 4=fits and ascii, ...).
		    int,    ///< Number of command-line arguments.
		    string * ///< Pointer to array of command-line arguments.
       );
   /** \brief Time integration
    * This is the main part of the code -- It does all the time integration
    * until the stopping condition is reached and then returns.
    * It calls a sequence of functions to advance the time by one timestep,
    * all in a loop which runs until end-of-sim is reached.
    * */
   virtual int Time_Int();
   /** \brief finalise the simulation, clean up, delete data.
    * This function finished the simulation gracefully (hopefully!).
    * */
   int Finalise();
   //---------------------------------------
  protected:
   //---------------------------------------
   //---------------------------------------
   // Data Variables common to all implementations.
   class DataIOBase *dataio; ///< pointer to class for reading/writing data.
   class DataIOBase *textio; ///< pointer to class for reading/writing textdata.
   class FV_solver_base *eqn;  ///< Pointer to equations to solve, initialised to
                               ///< some derived class at runtime when the equations are set.

#ifdef CHECK_MAGP
   ///
   /// This is only for a test problem -- it checks the magnetic
   /// pressure on the full domain and outputs it to screen
   ///
   void calculate_magnetic_pressure();
#endif // CHECK_MAGP

#ifdef BLAST_WAVE_CHECK
   ///
   /// If running a spherical blast wave, calculate the shock position and output to screen.
   ///
   void calculate_blastwave_radius();
#endif // BLAST_WAVE_CHECK

   /** \brief Get cell dimensions from UniformGrid class.
    * 
    * This requests the cell dimensions from UniformGrid and assigns values to 
    * the appropriate variables in GridParams (dx,dA,dV).
    */
   int get_cell_size();
   /** \brief See if any command-line arguments should override those
    * specified in the IC file, and if so, reset the parameters.
    * */
   int override_params(int,      ///< Number of command-line arguments.
		      string *  ///< Pointer to array of command-line arguments.
		      );
   /** \brief Initialise the correct Equations to solve, based on paramters.
    * */
   int set_equations();

  ///
  /// Setup cell extra data through the cell_interface class CI.
  ///
  void setup_cell_extra_data();
    
   /** \brief initialise the grid class with appropriate parameters.
    * 
    * This function sets up the appropriate grid; so far I only have a 
    * UniformGrid class -- uniform, cartesian, finite volume grid.
    * */
   virtual int setup_grid();
   /** \brief Determines what kind of boundary conditions are needed.
    * Sets gp.Nbc to the appropriate value for the order of accuracy used.
    * \retval 0 success
    * \retval 1 failure
    * */
   int boundary_conditions();
   /** \brief Delete any init data and make sure things are ready to go.
    *  */
   virtual int ready_to_start();
   
  ///
  /// Decide if I need to setup MP class, and do it if i need to.
  ///
  virtual int setup_microphysics();

  ///
  /// Decide if I need to setup RT class, and do it if i need to.
  ///
  virtual int setup_raytracing();

  ///
  /// Check for any time-evolving radiation sources, and read the evolution
  /// file if there are any.  Data is stored in global struct SimPM.STAR[v]
  ///
  virtual int setup_evolving_RT_sources();



   /********************************************************************/
   // Functions called by time_int()

  // ***********************************************************************
  // *********************** TIMESTEP CALCULATION **************************
  // ***********************************************************************

  /// flag: true if timestep limit needs raytracing column densities
  bool FVI_need_column_densities_4dt;

  int FVI_nheat; ///< number of RT heating sources
  int FVI_nion;  ///< number of ionising sources
  /// vector of RT heating sources, of size FVI_nheat
  std::vector<struct rt_source_data> FVI_heating_srcs;
  /// vector of RT ionising sources, of size FVI_nion
  std::vector<struct rt_source_data> FVI_ionising_srcs;

  ///
  /// Calculate the appropriate timestep.
  /// 
  /// For a uniform grid, all cells have the same timestep equal to the minimum
  /// of all calculated steps.  This function calls two functions, one to get 
  /// the microphysics timestep (if needed), and another to get the dynamics
  /// timestep.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  virtual int calc_timestep();

  ///
  /// Check for any time-evolving radiation sources, and update source
  /// properties from global struct SimPM.STAR[v] if needed
  ///
  virtual int update_evolving_RT_sources();

  ///
  /// Calculate the microphysics timestep, based on heating/cooling and reaction
  /// rates.  Returns the minimum timestep of the local grid (negative if error).
  /// 
  double calc_microphysics_dt();

  ///
  /// Old microphysics timescales calculation with no radiation field.
  ///
  double get_mp_timescales_no_radiation();

  ///
  /// New microphysics timescales calculation with pre-calculated radiation field.
  ///
  double get_mp_timescales_with_radiation();

  ///
  /// Calculate the dynamics timestep, based on the Courant condition that
  /// the fastest signals cannot cross a full cell in a single step.  Returns
  /// the minimum timestep on the local grid, or negative if an error occurs.
  ///
  double calc_dynamics_dt();
 

#ifdef THERMAL_CONDUCTION
  ///
  /// If doing thermal conduction, this calculates the max. timestep we can use
  /// without conduction changing the gas temperature by more than 30%.
  /// It also calls the solver function set_thermal_conduction_Edot() so it 
  /// knows what the flux in and out of cells is.  This Edot value is multiplied
  /// by the timestep dt in eqn->preprocess_data().
  ///
  double calc_conduction_dt_and_Edot();
#endif // THERMAL CONDUCTION

  ///
  /// Limits the timestep based on a number of criteria.  Makes sure we don't
  /// overshoot the finish-time, or the next output-time, and that we don't 
  /// increase the timestep by a large factor from one step to the next (this
  /// can affect stability).
  ///
  void timestep_checking_and_limiting();

  // ***********************************************************************
  // *********************** TIMESTEP CALCULATION **************************
  // ***********************************************************************
  ///
  /// Advance the simulation by one time step.  This is a wrapper function,
  /// which calls one of two specific algorithms.  If we have no raytracing,
  /// or if we are doing the C2-ray type ray-tracing which does not limit 
  /// the timestep by the ionisation time, then the original algorithm is
  /// used.  If we are doing the newer update, where only column densities
  /// are calculated in the ray-tracing update, then the timestep is limited
  /// so that no microphysics quantities change by more than a certain 
  /// percentage, then the microphysics is updated and finally the dynamics 
  /// is updated.
  ///
  /// This is the main time integration function.  If we're first order,
  /// it calls calc_dU once, which calculates the update for each cell in
  /// all directions (dimensionally unsplit method), and then adds dU to U
  /// to get the new state vectors.  It then increments the time and timestep
  /// counter.
  ///
  /// If we're second order, it halves the timestep, calls calc_dU once, and 
  /// adds dU to U to get the half-time state vector cpt->Ph, and then 
  /// doubles the timestep (to the full value), calls calc_dU again, and
  /// updates the full timestep to get the new state vector, then increments
  /// the time and timestep.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int advance_time();

  ///
  /// This is the original time integration algorithm used for JM's PhD thesis.
  /// The timestep is calculated first, based on the Courant condition and
  /// possibly microphysics timescales (excluding photoionisation time).
  /// The dynamics is updated first, and then the microphysics update is called.
  /// This *must* use the C2-Ray type algorithm to integrate the microphysics as
  /// rays are being traced (if, of course, we are doing raytracing).
  ///
  /// This is really a legacy function now, to enable me to run the 
  /// implicit integrator if I want to.
  ///
  int timestep_dynamics_then_microphysics();


#ifdef NEW_TIME_UPDATE
  ///
  /// This performs a first-order-accurate (in time) timestep for
  /// dynamics, microphysics, thermal conduction, everything.
  /// If the order-of-accuracy parameter is OA1 then it assumes this
  /// is a full step and updates both P[] and Ph[].
  /// If it is OA2, then the functions assumes this is the half-step
  /// as part of a full second-order step, so it only updates Ph[].
  /// It advances by an interval of dt regardless of OA1 or OA2, so 
  /// if this is a half step, you should pass 0.5*dt to the function.
  ///
  int first_order_update(
          const double, ///< dt, time interval to advance by.
          const int     ///< time order of accuracy OA1/OA2.
          );

  ///
  /// This performs a second-order-accurate (in time) timestep for
  /// dynamics, microphysics, thermal conduction, everything.
  /// This performs the second part of a full second-order step, so
  /// the half-step must have been already called before this one.
  ///
  int second_order_update(
          const double, ///< dt, time interval to advance by.
          const int     ///< time order of accuracy (must be OA2).
          );
  
  ///
  /// This function does some checking on radiation sources to see
  /// what microphysics update to call, then calls one of 
  /// calc_microphysics_dU_general_RT(),
  /// calc_microphysics_dU_JMs_C2ray_RT(), or
  /// calc_microphysics_dU_no_RT().
  ///
  int calc_microphysics_dU(
            const double ///< dt, timestep to integrate MP eqns.
            );

  ///
  /// This calculates the change in internal energy and ion fractions
  /// for a timestep dt, by integrating the microphysics equations
  /// for the full timestep, storing the result in a temporary array,
  /// and differencing the initial and final states.
  /// This version is for microphysics integrations where there are
  /// radiation sources involved.
  ///
  int calc_microphysics_dU_general_RT(
        const double   ///< dt, timestep to integrate
        );


  ///
  /// This calculates the change in internal energy and ion fractions
  /// for a timestep dt, by integrating the microphysics equations
  /// for the full timestep, storing the result in a temporary array,
  /// and differencing the initial and final states.
  /// This version is for microphysics integrations where I am using
  /// my old implicit photoionisation integrator based on C2-ray.
  ///
  int calc_microphysics_dU_JMs_C2ray_RT(
          const double ///< dt, timestep to integrate
          );


  ///
  /// This calculates the change in internal energy and ion fractions
  /// for a timestep dt, by integrating the microphysics equations
  /// for the full timestep, storing the result in a temporary array,
  /// and differencing the initial and final states.
  /// This version is for microphysics integrations where there are
  /// no radiation sources (e.g. pure heating+cooling, or collisional
  /// processes only).
  ///
  int calc_microphysics_dU_no_RT(
          const double ///< dt, timestep to integrate
          );

  ///
  /// This calculates the change in the state vector for each point
  /// due to the dynamics, for a timestep dt, using either 1st or 
  /// 2nd order accuracy in space.
  /// It calls eqn->preprocess_data(), then set_dynamics_dU(), and
  /// finally eqn->PostProcess_dU().
  /// set_dynamics_dU() is the function that used to be called
  /// calc_dU().
  ///
  int calc_dynamics_dU(
          const double, ///< dt, timestep to integrate
          const int     ///< spatial order of accuracy for update.
          );

  ///
  /// This function used to be called calc_dU -- for every column of
  /// simulation data on the grid this calls dU_Column() to get the
  /// changes in the state vectors arising from the hydrodynamics
  /// over the timestep dt.  It uses the requested spatial order of
  /// accuracy.  This function also loops over all directions on the 
  /// grid that are active.
  ///
  int set_dynamics_dU(
            const double,    ///< dt, timestep for this calculation
            const int        ///< space OOA for this calculation
            );

  ///
  /// Calculate dU, rate of change of conserved variables, in a 1D
  /// column of grid cells, according to the fluid dynamics equations.
  ///
  /// This runs through every cell in a 1D column in turn, and calculates the flux
  /// between the cell in question and its neighbour to the right, by obtaining
  /// an interface flux.
  /// 
  /// It then calculates dU for each cell according to the exact formula (if the 
  /// flux calculation were exact) given by Toro eq.5.76\n
  /// \f$ U_i^{n+1}-U_i^n =dU = \frac{\Delta t}{\Delta x}(F_{i-\frac{1}{2}} -F_{i+\frac{1}{2}}) \f$.
  ///
  int dynamics_dU_column(const class cell *, ///< starting point for column.
              const enum direction, ///< direction to traverse column in. 
              const enum direction, ///< opposite direction.
              const double,    ///< dt, timestep for this calculation
#ifdef TESTING
              const int,       ///< Time Order of accuracy to use.
#endif
              const int        ///< Spatial Order of accuracy to use.
              );

  ///
  /// This function takes the contents of each cell->dU[] vector and
  /// updates Ph[] the changes.  If we are on the full-step then it
  /// also updates P[] so that it and Ph[] are identical.
  ///
  int grid_update_state_vector(
            const double ,  ///< dt, timestep
            const int,      ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
            const int       ///< Full order of accuracy of simulation
            );
#endif // NEW_TIME_UPDATE


#ifndef NEW_TIME_UPDATE
  ///
  /// This is the newer time integration scheme, where column densities are
  /// calculated for each radiation source first and stored as cell data.  This
  /// allows a microphysics timestep to be calculated which includes ionisation
  /// and radiative heating, so we can limit the timestep correctly.
  /// Then the microphysics is updated, and then the dynamics is updated in a 1st
  /// or 2nd order algorithm, as described in Falle, Kommisarov, & Joarder (1998).
  ///
  int timestep_microphysics_then_dynamics();
#ifdef TEST_SECOND_ORDER
  /// Testing second order microphysics for the raytracing.
  int timestep_mp_dyn_second();
#endif



   /** \brief Updates the microphysics of each cell in turn, if required.
    * 
    * If we aren't doing microphysics, then this returns 0;
    * If we are not on the full timestep, this returns 0;
    * If we are on the full timestep, and are doing microphysics, then this
    * runs through all the cells, updating the tracers and internal energy
    * appropriately, and then returns 0 on successful updating.
    * 
    * \retval 0 success.
    * \retval 1 failure.
    */
   int update_microphysics(const double,   ///< timestep to take
			   const int, ///< fraction of timestep we are on
			   const int  ///< full no. of fractional timesteps.
			   );
  

  ///
  /// microphysics update where there are no heating or ionising radiation fields
  /// explicitly calculated.
  ///
  int update_microphysics_no_RT(const double, // timestep to integrate
                                const int   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                );

  ///
  /// microphysics update used in JM's PhD thesis, based on the C2-ray method, for a
  /// single monochromatic source of ionising radiation, where the microphysics is
  /// updated as the rays are traced outwards.
  ///
  int update_microphysics_JMs_C2ray_RT(const double, // timestep to integrate
                                       const int   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                       );


  ///
  /// General microphysics update including (for sure) at least one ionising source, and 
  /// optionally some UV heating sources, either diffuse or associated with an ionising
  /// source.  Diffuse and direct column densities were pre-calculated before calling
  /// this function, so these are then organised and passed into the MPnew() interface
  /// function.
  ///
  int update_microphysics_general_RT(const double, // timestep to integrate
                                     const int   ///< Whether it is TIMESTEP_FULL, TIMESTEP_FIRST_PART
                                     );

   ///
   /// This function encapsulates the calc_dU function and the time update.
   ///
   /// Called by advance_time(), calls calc_dU(),
   /// eqn->CellAdvanceTime().  Also tracks energy and momentum
   /// conservation properties of simulation if testing the code.  If
   /// there is any pre- or post-processing to be done (e.g. viscosity)
   /// then that is done in appropriately named functions called
   /// before and after the calc_dU() function.
   ///
   int update_dynamics(const int, ///< Space order of acc for this call.
		       const int  ///< Time order of acc for this call.
		       );
   


   /** \brief Calculate dU, rate of change of conserved variables.
    * 
    * This is a multidimensional routine, which identifies 1D columns of
    * points and calls dU_column() on them all, in each grid direction.
    * 
    * It is a dimensionally unsplit routine, calculating updates in each 
    * dimension and then letting advance_time() do the actual time update 
    * later.
    * 
    * It uses the exact formula (if the flux calculation were exact):
    * \f$ U_i^{n+1}-U_i^n =dU = \frac{\Delta t}{\Delta x}\left[
    *     (F_{i-\frac{1}{2}} -F_{i+\frac{1}{2}}) +
    *     (G_{j-\frac{1}{2}} -G_{j+\frac{1}{2}}) +
    *     (H_{k-\frac{1}{2}} -H_{k+\frac{1}{2}}) \right] \f$, where it
    * is assumed that we are in three (Cartesian) spatial dimensions.  Terms are 
    * removed from this as appropriate for the number of dimensions specified.
    * This formula is appropriate for cartesian geometry, but the code calculates
    * the divergence differently if we are using cylindrical or other coordinates.
    * I don't think it would work for non-orthonormal coords.
    * */
   int calc_dU(const int, ///< Space order of acc for this call.
	       const int  ///< Time order of acc for this call.
	       );

   /** \brief Calculate dU, rate of change of conserved variables, in a 1D column
    * 
    * This runs through every cell in a 1D column in turn, and calculates the flux
    * between the cell in question and its neighbour to the right, by obtaining
    * an interface flux.
    * 
    * It then calculates dU for each cell according to the exact formula (if the 
    * flux calculation were exact) given by Toro eq.5.76\n
    * \f$ U_i^{n+1}-U_i^n =dU = \frac{\Delta t}{\Delta x}(F_{i-\frac{1}{2}} -F_{i+\frac{1}{2}}) \f$.
    * This is for cartesian geometry, and the form is different for cylindrical.  The 
    * code calls a different function for different geometries.
    * */
   int dU_column(const class cell *, ///< starting point for column.
		 const enum direction, ///< direction to traverse column in. 
		 const enum direction, ///< opposite direction.
		 const int,       ///< Spatial Order of accuracy to use.
		 const int        ///< Time Order of accuracy to use.
		 );
#endif // not NEW_TIME_UPDATE


  ///
  /// Run through all diffuse and direct radiation sources and calculate column
  /// densities through the grid for each one.  Tau, DTau, and Vshell are stored
  /// in extra_data[i] for each cell.
  ///
  int calculate_raytracing_column_densities();

   /**  \brief Output the data to file if required.
    * This checks if I want to output data in this timestep, then
    * checks what format to write in, and calls the appropriate 
    * function to write the data.
    * \retval 0 success
    * \retval 1 failure
    * */
   virtual int output_data();



   /**  \brief Check if sim should stop.
    * 
    * For shock tube problems, I stop the simulation whenever a disturbance
    * reaches the edge of the domain.  Most problems set the finishtime to 
    * some value, and the end-of-sim criteria is if the current simtime has reached
    * finishtime or not.
    * \retval 0 success
    * \retval 1 failure
    * */
   int check_eosim();    // Checks for end of simulation.

   /** \brief Checks Total energy relative to initial value, and prints a message if not.*/
   int check_energy_cons();

   /** \brief Calculates total values of conserved quantities. */
   int initial_conserved_quantities();
}; // IntUniformFV
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#ifdef PARALLEL
/** \brief The Parallel implementation of the Uniform FV Integrator.
 * 
 * This class reimplements some functions of IntUniformFV, so that they
 * work on multiple processors with the domain split between them.
 * */
class ParallelIntUniformFV : public IntUniformFV
{
  public:
   ParallelIntUniformFV();
   ~ParallelIntUniformFV();
   /** \brief initialisation.
    * This function checks if we are reading from single or multiple files,
    * modifies the input file string accordingly, checks the file exists, 
    * and then calls the IntUniformFV::Init() function.
    *
    * \retval 0 success
    * \retval 1 failure
    * */
   int Init(string, ///< Name of input file.
	    int,    ///< Type of File (1=ASCII, 2=FITS, 3=fitstable, 4=fits and ascii, ...).
	    int,    ///< Number of command-line arguments.
	    string * ///< Pointer to array of command-line arguments.
       );
  /** \brief Time integration
    * This is the main part of the code -- It does all the time integration
    * until the stopping condition is reached and then returns.
    * It calls a sequence of functions to advance the time by one timestep,
    * all in a loop which runs until end-of-sim is reached.
    * 
    * Parallel version has an AllReduce operation, where I check if the runtime of 
    * any processor is more than a fixed walltime, and if so set eosim to true and
    * finish.  This is because ICHEC machines have a maximum runtime limit for their
    * simulations on some of the queues, and I want to make sure I have an output 
    * near the end of the allowed runtime.
    * */
   int Time_Int();
  protected:
   /** \brief initialise the grid class with appropriate parameters.
    * 
    * This function sets up the appropriate grid; for parallel execution
    * I need to define the domain of my grid, and then pass the appropriate
    * parameters to the UniformGrid class.
    * */
   int setup_grid();

   /**  \brief   Calculate the appropriate timestep for all processors
    * 
    * For a uniform grid, all cells have the same timestep equal to the minimum
    * of all calculated steps.  This function calls the calc_timestep() function
    * for the local grid, and then gets the min of all processor's local
    * timesteps, and uses that as the timestep.
    * \retval 0 success
    * \retval 1 failure
    * */
   int calc_timestep();
   /** \brief Decide if I need to setup RT class, and do it if i need to. */
   virtual int setup_raytracing();
}; // ParallelIntUniformFV
#endif // PARALLEL

#endif // if not GRID_H
