/// \file raytracer_SC.h
/// \brief Contains Declaration of Short Characteristic Tracer Class.
/// 
/// Author: Jonathan Mackey
/// 
/// Modifications:
///  - 2008-03-19 Wrote file.
///  - 2010-01-15 JM: Turned off NO_SOURCE_CELL_GEOMETRY (seems to make very little difference)
///  - 2010-01-19 JM: put weighting schemes in new function to avoid code duplication in parallel version.
///  - 2010-01-20 JM: comments.
///  - 2010-01-22 JM: Moved #defs for (NON_)CELL_CENTRED_SRC from raytracer_SC.h to global.h
///  - 2010-01-24 JM: comments only, also forced either cell or corner centring to be defined.
///
///  - 2010.07.23 JM: New RSP source position class interface.
///
/// - 2011.02.25 JM: Changed Add_Source() interface, so that all relevant parameters
///    are passed in to function as pointer to element in SimPM.RS
///    Removed HCORR ifdef around new code.
///    Removed NEW_RT_MP_INTERFACE flags.
///
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function. 
/// - 2011.04.15 JM: Moved old ProcessCell() to ProcessCell_TimeUpdate(), and the 
///    new ProcessCell() function is a chooser based on source properties.  It can
///    just calculate column-densities, or do the TimeUpdate.
/// - 2011.04.18 JM: Added interpolate_2D_RHO() function
/// - 2011.04.23 JM: Added TauMin[] array of values of TauMin for each source, ordered
///    by source id.  This removes the need for interpolate_2D_RHO() so I got rid of it.
///    Now the interpolate functions read TauMin[src_id] and apply that to the interpolation.
///
/// - 2011.10.14 JM: Added new interface for rt_source_data structs.
/// - 2011.10.17 JM: Debugging.  Getting new implicit integrator working.
/// - 2011.12.15 JM: Removed the NEW_RT_STRUCTURE flag since it is not default.
/// - 2012.01.16 JM: Gave update_RT_source_properties() a real function.
/// - 2012.03.31 JM: Made the Add_Source() function more modular, so I don't need the
///    ifdefs in the serial version for setting Vshell.
/// - 2013.08.20 JM: Modifications so a given source can have an
///    array of optical depths rather than a single value per source.
///    This is tough going, and is only half-way done so far.

#ifndef RAYTRACING_H
#define RAYTRACING_H


#include "global.h"

#ifndef CELL_CENTRED_SRC
#ifndef NON_CELL_CENTRED_SRC
#error "MUST define one of CELL_CENTRED_SRC and NON_CELL_CENTRED_SRC in global.h"
#endif
#endif

#ifdef CELL_CENTRED_SRC
#error "CELL_CENTRED_SRC has not been used for a LONG time and may be buggy."
#error "USE AT YOUR OWN RISK.  DELETE THIS ERROR IN source/raytracing/raytracer_SC.h"
#endif

#define NO_SOURCE_CELL_GEOMETRY ///< this stops the code changing tau in the
                                ///< source cell for different angles (makes no difference)
#define INTERPOLATE_METHOD 0 /// 0=C2Ray mintau=0.7, 1=Geometric,
                             /// 2=Geometric in exp(-tau), 3=quadratic/exp(-tau) !!!ONLY C2-RAY IS ANY GOOD!!!

/** \brief Related quadrants in x-y plane to indices for the order they are to be traversed.
 * The source is in Q1, so it comes first, then Q2, then Q4, and finally Q3.  Q2 and Q4 could
 * be interchanged, but this is the way I'm doing it.
 * */
enum Quadrants {Q1=0, Q2=1, Q4=2, Q3=3};
#define NQUADS 4
/** \brief This arranges the octants in integers according to how I am tracing them.
 * The first four octants are in ZP dir, last four in ZN, and each group of four 
 * traces an anticlockwise circle in the x-y plane, starting with (XP,YP).
 * */
enum Octants {OCT1=0, OCT2=1, OCT4=2, OCT3=3, OCT5=4, OCT6=5, OCT8=6, OCT7=7};


///
/// Struct to hold info on radiation sources.  An extension of rad_src_info in global.h
///
struct rad_source {
  ///
  /// pointer to source (set to SimPM.RS.source[id]) with the basic
  /// info about the source.
  ///
  struct rad_src_info *s;

  class cell *sc; ///< nearest cell to source.
  bool src_on_grid; ///< true if source is at a grid cell.
  int ipos[MAX_DIM];    ///< source position in integer form (grid units, dx=2).

  //int id;   ///< source id.
  //double *pos;  ///< source position.
  //double strength;  ///< source strength.

  //bool at_infinity; ///< True if source is at infinity.
  //int type;   ///< Type of source.  RT_SRC_DIFFUSE or RT_SRC_SINGLE.
  //int effect;      ///< Either UV heating or photoionisation+heating.
  //int opacity_src; ///< What provides the opacity: RT_OPACITY_TOTAL, RT_OPACITY_MINUS, RT_OPACITY_TRACER.
  //int opacity_var; ///< optional tracer variable index in state vector, for opacity calculation.
  //int update; ///< how the source is updated: RT_UPDATE_IMPLICIT=1, RT_UPDATE_EXPLICIT=2

  ///
  /// This struct is used by the code to pass cell and source data
  /// to the microphysics integrator.  It contains the relevant source
  /// information, and also some cell-source geometry information which 
  /// must be set on a cell-by-cell basis as rays are traced.
  /// The struct is declared in source/microphysics/microphysics_base.h
  ///
  struct rt_source_data data;
};




/** \brief Special simpler raytracer for if we have parallel rays coming in 
 * along one axis.  This is appropriate for a single source at infinity.
 *
 * For parallel rays, the source strength is the flux entering the grid. from the
 * direction of the source.
 *
 * c->col tracks the optical depth from the source, along a ray through the cell
 * centre, and up to the point the ray exits the cell.
 */
class raytracer_USC_infinity : public RayTracingBase {
 public:
  /** \brief Constructor. */
  raytracer_USC_infinity(class GridBaseClass *, ///< Pointer to grid
			 class MicroPhysicsBase * ///< Pointer to MicroPhysics Class.
			 );
  ~raytracer_USC_infinity(); ///< destructor.


  /** \brief Adds a source to the list of sources to trace.
   * 
   * Returns the sources id, which starts at zero and increases 
   * linearly.  So if we have 10 sources, and add another source,
   * it's id with be 10.  
   *
   * For a source at infinity, the source strength refers to the flux
   * arriving at the grid, as we can't use luminosity.
   */
  virtual int Add_Source(
        struct rad_src_info * ///< ptr to source info.
        );


  /** \brief Processes a source's effect on the grid over a timestep.
   * This is the 1D algorithm for parallel rays along one coordinate axis.
   */
   virtual int RayTrace_SingleSource(const int,   ///< Source id
				     const double, ///< Timestep
				     const double  ///< EOS gamma
				     );
  ///
  /// Just calculate the column densities required for RT.
  ///
  virtual int RayTrace_Column_Density(
                const int,    ///< Source id
                const double, ///< Timestep
                const double  ///< EOS Gamma.
                );

   /** \brief Prints list of sources with id, location, strength. */
   void Print_SourceList();
   /** \brief Returns the number of sources to track. */
   int NSources();

  ///
  /// Returns whether we are doing an implicit (==0) or an explicit (==1)
  /// integration of the raytracing/microphysics.
  ///
  int type_of_RT_integration() {return type_of_RT_int;}

  ///
  /// This sets the number of ionising and UV heating sources of radiation,
  /// and makes sure the rt_source_data structs are populated correctly.
  /// It can be used to change the radiation sources if e.g. the luminosity
  /// changes over time.
  ///
  virtual void update_RT_source_properties(
                  const struct rad_src_info * ///< ptr to source info.
                  );

  /// Returns the number of ionising sources
  int N_ionising_sources() {return N_ion_srcs;}
  /// Returns the number of UV heating sources
  int N_heating_sources() {return N_uvh_srcs;}

  ///
  /// This function copies the ionising source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_ionising_src_list(
                std::vector<struct rt_source_data> & ///< list of data for ionising sources
                );

  ///
  /// This function copies the UV-heating source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_UVheating_src_list(
                std::vector<struct rt_source_data> & ///< list of data for UV-heating sources
                );

 protected:
   int ndim;  ///< dimensionality of grid to traverse.
   double delt; ///< Timestep to step by.
   double gamma; ///< EOS gamma.

  /// Flag for the type of RT integration we are performing 1=IMP,2=EXP.
  int type_of_RT_int;
  /// Number of ionising radiation sources.
  int N_ion_srcs;
  /// Number of UV-heating sources.
  int N_uvh_srcs;

  /// cell/source data on IONISING srcs.
  std::vector<struct rt_source_data> ION_data;
  /// cell/source data on UV-HEATING srcs.
  std::vector<struct rt_source_data> UVH_data;

  std::vector<rad_source> SourceList; ///< List of sources, with id, position, strength.

  ///
  /// This function increments N_ion_srcs or N_uvh_srcs, and checks whether
  /// to change type_of_RT_int to IMPLICIT based on the properties of the 
  /// newly added source.
  ///
  void update_local_variables_for_new_source(
                const struct rad_source ///< newly added source
                );


  class GridBaseClass *gridptr; ///< pointer to grid, which should be setup
                                 ///< independently of the raytracer.

  class MicroPhysicsBase *mpptr; ///< pointer to microphysics class, which 
                                  ///< should be set up independently.

  struct rad_source *current_src; ///< pointer to source we are currently working with.
   
   /** \brief If a source is at infinity, we can trace each column one by one with 
    * a much simpler algorithm, so that's what this function does.
    */
  int trace_parallel_rays(const rad_source *,   ///< source we are dealing with.
			   const enum direction ///< direction to source at infinity.
			   );
   /** \brief Traces a 1D column from a starting cell, in a direction, to the edge of the grid. 
    * This function explicitly assumes rays parallel to the given direction.
    **/
  int trace_column_parallel(const rad_source *,  ///< source we are tracing from.
			     cell *,              ///< cell to start from.
			     const enum direction ///< direction we are looking.
			     );
   /** \brief Get column density to current cell from source, and path length through cell.
    * This assumes parallel rays along the column, appropriate for 1D or source at infinity.
    * */
  int cell_cols_1d(const rad_source *, ///< pointer to source struct.
		    class cell *,       ///< Current Cell
		    const enum direction, ///< direction to source. 
		    double *,           ///< column to cell.
		    double *            ///< path length of ray in cell.
		    );

  ///
  /// Do whatever is to be done to the cell due to the ray passing through.
  /// If doing the C2RAY update, this calls ProcessCell_TimeUpdate().  The
  /// behaviour of this function is set by rs->update and rs->opacity_src.
  ///
  virtual int ProcessCell(cell *,            ///< Current cell.
		   double[],            ///< Column to cell.
		   double,            ///< Path Length through cell.
		   const rad_source *, ///< pointer to source struct.
		   const double       ///< Timestep
		   );

  ///
  /// This does most of the work of Add_Source() for sources at infinity.
  /// It finds the direction of the source, adds its properties to the 
  /// local source-list maintained by the raytracing class, and makes sure
  /// that its properties are appropriate.
  ///
  virtual void add_source_to_list(
              struct rad_src_info * ///< source info.
              );

  ///
  /// Set Vshell in all cells for the current source (See Mellema et al. 2006,
  /// NewAst.).  This function works for serial, parallel, and sources at
  /// infinity.  Vshell is used for the photon-conserving photoionisation rate
  /// calculation.
  ///
  void set_Vshell_for_source(
              struct rad_source *
              );

  ///
  /// Set Vshell in the current cell (See Mellema et al. 2006, NewAst.).
  /// For parallel rays this is just ds.  It is used for the photon-conserving
  /// photoionisation rate calculation.
  ///
  virtual void set_Vshell_in_cell(
            cell *, ///< current cell.
            double,            ///< Path Length through cell.
            const rad_source * ///< pointer to source struct.
            );

  ///
  /// For the C2RAY update, where we integrate the microphysics as we trace
  /// rays, we need to do quite a bit of work, so this is in its own function.
  ///
  virtual int ProcessCell_TimeUpdate(
            cell *,            ///< Current cell.
            double[],            ///< Column to cell.
            double,            ///< Path Length through cell.
            const rad_source *, ///< pointer to source struct.
            const double       ///< Timestep
            );


};





/** \brief More compicated RayTracer, for a source which can be on grid, off grid,
 * or at infinity in one coordinate direction.
 */
class raytracer_USC : public raytracer_USC_infinity {
 public:
  /** \brief Constructor. */
  raytracer_USC(class GridBaseClass *, ///< Pointer to grid
		class MicroPhysicsBase * ///< Pointer to MicroPhysics Class.
		);
  ~raytracer_USC(); ///< destructor.

  /// \brief Adds a source to the list of sources to trace.
  ///
  /// Returns the sources id, which starts at zero and increases 
  /// linearly.  So if we have 10 sources, and add another source,
  /// it's id with be 10.
  /// 
  /// ifdef CELL_CENTRED_SRC
  /// Currently if the source is on the grid, it will move it to the centre of
  /// the nearest cell to it.  This makes the algorithm simpler, but it would
  /// be advantageous to relax this requirement, for e.g. axisymmetric runs
  /// with an on-axis source.
  ///
  /// ifdef NON_CELL_CENTRED_SRC
  /// This ifdef moves a source to the nearest cell VERTEX if it is not there
  /// already.  This works for axisymmetric grids, improves parallel scaling
  /// for a small number of cores, and seems to be more accurate for Cartesian
  /// grids, so it is recommended to always set this!
  ///
  virtual int Add_Source(struct rad_src_info * ///< ptr to source info.
			  );

   /** \brief Processes a source's effect on the grid over a timestep.
    * This is a generic algorithm for 1D,2D,3D.  It calls different 
    * functions depending on the dimensionality of the simulation.
    */
   virtual int RayTrace_SingleSource(const int,    ///< Source id
				     const double, ///< Timestep
				     const double  ///< eos gamma.
				     );


  /** \brief Prints list of sources with id, location, strength. */
  void Print_SourceList();

  /** \brief Returns the number of sources to track. */
  int NSources();
   
  protected:
  enum direction 
    dir1[8], ///< list of x-directions for tracing octants in order.
    dir2[8], ///< list of y-directions for tracing octants in order.
    dir3[8]; ///< list of z-directions for tracing octants in order.

  ///
  /// Direction from cells in question to source.  To be set
  /// for each line/quadrant/octant in RayTraceSource() function.
  ///
  std::vector<enum direction> SrcDir;
  ///
  /// min-Tau values for source column density interpolation, ordered by src id.
  ///
  std::vector<double> TauMin;
  
	
  /** \brief Find the source cell, or if the source is off the grid, find the
   * nearest cell to the source. */
  class cell * find_source_cell(double * ///< position of source.
															 	);

   /** \brief Set the source position to be the centre of a cell, even if off grid.*/
  void centre_source_on_cell(
          double *, ///< position of source.
          enum axes              ///< axis to find source along.
          );
   
   /** \brief This goes along an axis until it gets to the source location. 
    * 
    * This is called by findSourceCell() if the source is definitely on the
    * grid along a certain coordinate axis, and it just goes along the axis until
    * it finds the cell the source is in, and moves the source to the cell centre, if
    * it is not there already.
    * */
  virtual class cell * find_src_on_grid(
        double *, ///< position of source.
        cell *,                ///< starting cell.
        enum axes              ///< axis to find source along.
        );

  /** \brief This will return a pointer to the source cell, or the on-grid cell nearest the off-grid source.
   * Works for 2D and 3D.
   */
  void find_closest_cell(
        const rad_source *, ///< pointer to source
        cell *,             ///< cell to move to cell nearest to source.
        enum direction *    ///< array to put in directions to source from cell nearest.
        );
   /** \brief Given a source cell, assign a list of start cells for each line/quad/octant
    * we need to trace along (some may be null).
    */
  void set_startcells(
        cell *,          ///< source cell, or cell nearest to source if off grid (const).
        cell **,         ///< list of startcells for each line/quadrant/octant (to be assigned).
        enum direction * ///< list of directions to source if it is off grid (const).
        );


   /** \brief Traces a 1D column from a starting cell, in a direction, to the edge of the grid. */
   int trace_column(const rad_source *,  ///< source we are tracing from.
		    cell *,              ///< cell to start from.
		    const enum direction ///< direction we are looking.
		    );


   /** \brief This ray-traces a quadrant of the grid plane starting at the input cell.
    * 
    * It calls the trace_column() function to trace x-columns and moves in the y-direction
    * after each column, so that the whole quadrant gets covered.
    * */
   int trace_plane(const rad_source *,   ///< source we are dealing with.
		   cell *,               ///< starting cell in quadrant
		   const enum direction, ///< x-direction from starting cell to go in.
		   const enum direction  ///< y-direction from starting cell to go in.
		   );


   /** \brief This ray-traces an octant of the grid volume starting at the input cell.
    * 
    * It calls the trace_plane() function to trace x-y planes at each z-cell, and
    * moves through the column of z-cells in the octant.
    * */
   int trace_octant(const rad_source *,   ///< source we are dealing with.
		    cell *,               ///< starting cell in octant
		    const enum direction, ///< x-direction from starting cell to go in.
		    const enum direction, ///< y-direction from starting cell to go in.
		    const enum direction  ///< z-direction from starting cell to go in.
		    );
   /** \brief Get column density to current cell from source, and path length through cell.
    * This uses interpolation for 2D and 3D grids.  For parallel rays use 
    * get_cell_columns_parallel()
    * */


   int get_cell_columns(const rad_source *, ///< pointer to source struct.
			class cell *,       ///< Current Cell
			double *,           ///< column to cell.
			double *            ///< path length of ray in cell.
			);


   /** \brief Get column density to current cell from source, and path length through cell.
    * This uses interpolation for 2D grids.  Various weighting schemes have been tried, 
    * and the best one is the one that isn't commented out!
    * */
   int cell_cols_2d(const rad_source *, ///< pointer to source struct.
		    class cell *,       ///< Current Cell
		    double *,           ///< column to cell.
		    double *            ///< path length of ray in cell.
		    );


   /** \brief Get column density to current cell from source, and path length through cell.
    * This uses interpolation for 3D grids.
    * */
   int cell_cols_3d(const rad_source *, ///< pointer to source struct.
		    class cell *,       ///< Current Cell
		    double *,           ///< column to cell.
		    double *            ///< path length of ray in cell.
		    );


   ///
   /// Short Characteristic Method of getting column density to cell.
   ///
   virtual void col2cell_2d(
          const rad_source *,     ///< source we are working on.
          const cell *,           ///< cell to get column to.
          const enum direction,   ///< face ray enters cell through.
          const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
          const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
          double [] ///< Column densities.
          );


   ///
   /// Short Characteristic Method of getting column density to cell.
   ///
   virtual void col2cell_3d(
          const rad_source *,     ///< source we are working on.
          const cell *,           ///< cell to get column to.
          const enum direction,   ///< face ray enters cell through.
          const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
          const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
          double [] ///< Column densities.
          );


   ///
   /// Apply the appropriate weighting scheme to get an interpolated optical depth for 2D
   /// This version assumes the inputs are actually optical depths Tau1, Tau2.
   ///
   double interpolate_2D(
            const int, ///< source id (to get TauMin value).
            const double, ///< delta = min(abs(dy/dx),abs(dx/dy));
            const double, ///< first optical depth tau_1
            const double  ///< second optical depth tau_2
            );
   
   ///
   /// Apply the appropriate weighting scheme to get an interpolated optical depth for 3D
   ///
   double interpolate_3D(
            const int, ///< source id (to get TauMin value).
            const double, ///< delta0 = abs(dy/dx)
            const double, ///< delta1 = abs(dz/dx)
            const double, ///< first optical depth tau_1
            const double, ///< second optical depth tau_2
            const double, ///< third optical depth tau_3
            const double  ///< fourth optical depth tau_4
            );
   

  ///
  /// This sets the TauMin value for each radiation source, depending on what
  /// kind of source it is.  Now the raytracer just calculates projected mass
  /// along rays, so it is not exactly Tau, so TauMin is a projected density
  /// which corresponds to Tau=0.6 or 0.7, for the appropriate source opacity.
  ///
  /// This is a bit of a hack, and should eventually be replaced by something
  /// more flexible and automatic, by e.g. getting info from microphysics
  /// about the opacity per unit mass.
  ///
  void set_TauMin_for_source(
                const struct rad_source ///< new source
                );


  ///
  /// This does most of the work of Add_Source() for sources at infinity.
  /// It finds the direction of the source, adds its properties to the 
  /// local source-list maintained by the raytracing class, and makes sure
  /// that its properties are appropriate.
  ///
  virtual void add_source_to_list(
              struct rad_src_info * ///< source info.
              );

  ///
  /// Set Vshell in the current cell (See Mellema et al. 2006, NewAst.).
  /// For parallel rays this is just ds.  It is used for the photon-conserving
  /// photoionisation rate calculation.
  ///
  virtual void set_Vshell_in_cell(
            cell *, ///< current cell.
            double,            ///< Path Length through cell.
            const rad_source * ///< pointer to source struct.
            );

#ifdef CELL_CENTRED_SRC
   /** \brief Source Cell is processed differently, as there is no column density
    * to the cell, and the column through the cell is different.
    * */
   int ProcessSourceCell(cell *,             ///< Current cell.
			 const rad_source *, ///< pointer to source struct.
			 const double       ///< Timestep
			 );
#endif // CELL_CENTRED_SRC
};





#ifdef PARALLEL

/** \brief Distributed Memory Raytracer, for a source which can be on grid, off grid,
 * or at infinity in one coordinate direction.
 */
class raytracer_USC_pllel : public raytracer_USC {
 public:
  /** \brief Constructor. */
  raytracer_USC_pllel(class GridBaseClass *, ///< Pointer to grid
		class MicroPhysicsBase * ///< Pointer to MicroPhysics Class.
		);
  ~raytracer_USC_pllel(); ///< destructor.

  /// \brief Adds a source to the list of sources to trace.
  ///
  /// Differs from serial version in that it also calls 
  /// gridptr->setup_RT_boundaries() so that the grid can decide which (if any)
  /// extra corner boundaries it needs to set up for sending and receiving data.
  ///
  virtual int Add_Source(struct rad_src_info * ///< ptr to source info.
                        );

   /** \brief Processes a source's effect on the grid over a timestep.
    * This is a generic algorithm for 1D,2D,3D.  It calls different 
    * functions depending on the dimensionality of the simulation.
    *
    * This parallel version receives boundary column densities, calls
    * the serial raytracer, and then sends outgoing boundary column densities.
    */
   virtual int RayTrace_SingleSource(const int,    ///< Source id
				     const double, ///< Timestep
				     const double  ///< eos gamma.
				     );


  protected:
  /** \brief Short Characteristic Method of getting column density to cell. 
   * The parallel version assumes boundary cells exist, so doesn't check that we
   * are on the grid.
   * */
  virtual void col2cell_2d(
        const rad_source *,     ///< source we are working on.
        const cell *,           ///< cell to get column to.
        const enum direction,   ///< face ray enters cell through.
        const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
        const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
        double []               ///< Column densities.
        );

   /** \brief Short Characteristic Method of getting column density to cell.
    * The parallel version assumes boundary cells exist, so doesn't check that we
    * are on the grid.
    * */
  virtual void col2cell_3d(
        const rad_source *,     ///< source we are working on.
        const cell *,           ///< cell to get column to.
        const enum direction,   ///< face ray enters cell through.
        const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
        const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
        double []               ///< Column densities.
        );

};
#endif //PARALLEL


#endif // RAYTRACING_H
