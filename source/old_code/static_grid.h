/// \file static_grid.h
///
/// \brief Declares static uniform grid class
///
/// \author Jonathan Mackey
///
///
/// UniformGrid is designed to contain what every uniform finite
/// volume grid should need.  It is not intended to integrate
/// or solve PDEs by itself, but it should be set up by a class
/// which can operate on the grid to solve the PDEs.  In this way
/// UniformGrid should be suitable for solving any set of equations
/// which can be put on an n-dimensional finite volume grid (where
/// n<=3).
///
/// Modifications:\n
/// - 2007-11-06 Got going on adding in parallel grid and getting it
///    working.
/// - 2010-01-26 JM: Added get_iXmin/iXmax/iRange() functions to grid
///    class to give integer positions.
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; added
///    oneway-outflow boundary.
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010-07-24 JM: Added stellar wind class, and boundary
///    setup/update functions.
/// - 2010.10.04 JM: Added last-point-set flag.
/// - 2010.12.04 JM: Added distance functions which change depending
///   the grid geometry.  Added cylindrical and spherical grids.
///   This is all in an ifdef for now (GEOMETRIC_GRID)
/// - 2011.01.06 JM: New stellar wind interface.
/// - 2011.01.07 JM: I debugged the geometric grid functions, and now
///    it works very well!  Added iR_com(c) function to get the cell
///    centre-of-volume radius in integer code units.  This was a
///    problem before b/c the code and physical systems have 
///    different origins -- integer r=0 is at ix=-1.
/// - 2011.02.15 JM: Made 'Wind' a pointer so that it can be either
///    constant or evolving wind source class.
/// - 2011.02.24 JM: Worked on generalising multi-source radiative
///    transfer. Still some way to go to getting it functional.
/// - 2011.02.25 JM: Got rid of some vestigial data.
/// - 2011.04.06 JM: Added idifference_cell2cell() function for
///    thermal conduction calculation.
/// - 2012.05.15 JM: Added global iXmin,iXmax,iRange functions for
///    parallel grids.
/// - 2013.01.14 JM: NEW FILE FOR GRID V2, INHERITED FROM GRID V1.
///    Cut out parallel grid declaration (can go in a new file).
///    Tidied up comments, and added some new functions for the new
///    linked list covering all cells, including boundary data.

#ifndef STATIC_GRID_H
#define STATIC_GRID_H
#ifdef GRIDV2


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


/// I should work on this more.

#include <list>
using namespace std;

#include "stellar_wind_BC.h"
#include "coord_sys/VectorOps.h"
#include "coord_sys/VectorOps_spherical.h"



// ##################################################################
// ##################################################################


enum BoundaryTypes {
    PERIODIC   = 1,
      ///< periodic bcs.
    OUTFLOW    = 2,
      ///< outflow or absorbing bcs.
    INFLOW     = 3,
      ///< inflow bcs., value on boundary doesn't change.
    REFLECTING = 4,
      ///< reflecting bcs., hard wall.
    FIXED      = 5,
      ///< fixed bcs, means every point on boundary has same
      ///< unchanging value.
    JETBC      = 6,
      ///< A jet inflow boundary, internal boundary.
    JETREFLECT = 7,
      ///< Sort-of reflection for bi-directional jet, 
      ///< where normal B field passes through, but tangential is
      ///< reversed.
    DMACH      = 8,
      ///< Outflow boundary for double mach reflection test problem
      ///< only.
    DMACH2     = 9,
      ///< Fixed boundary on y=0, x in [0,1/6] fixed to postshock
      ///< state.
    BCMPI      =10,
      ///< boundary between processor domains in parallel grid.
    RADSHOCK   =11,
      ///< Boundary condition adjacent to cold wall for radiative
      ///< shock test problem.
    RADSH2     =12,
      ///< Outflow augmented with fixed outflow speed.
    ONEWAY_OUT =13,
      ///< One-way valve -- allows outflow but not inflow (zero 
      ///< gradient OR reflecting).
    STWIND     =14
      ///< Stellar wind sources exist, so apply internal boundaries.
};




// ##################################################################
// ##################################################################



///
/// Struct to contain all the information for a grid boundary.
///
struct boundary_data {
   enum direction dir; ///< Outward Normal direction of boundary (NO dir if internal).
   string type; ///< What type of boundary it is (Periodic, Absorbing, Fixed, Reflective, etc.).
   int itype;  ///< Integer flag for boundary type (per=1,out=2,inf=3,ref=4,fix=5).
   list<cell *> data; ///< STL linked list for boundary data cells.
   double *refval; ///< Optional reference state vector (e.g. for fixed BCs.).
};



// ##################################################################
// ##################################################################


/// The uniform finite-volume grid.
/// 
/// Can be 1,2, or 3 dimensional.  The grid cells must be cubic, so the number
/// of cells in each dimension should be in the right proportion to the 
/// length of the grid in each dimension.
/// 
/// A number of boundary conditions can be implemented; functions relating to
/// the boundaries are prefixed with \'BC\_\' as there are quite a few of them.
/// In many ways the boundaries are the most complicated part of the grid, once
/// it is set up.
///
class UniformGrid
: virtual public GridBaseClass,
  virtual public VectorOps_Cart
{
 protected:
  class stellar_wind *Wind; ///< stellar wind boundary condition class.
   int G_firstPtID;  ///< ID of the first point in the grid.
   int G_lastPtID;   ///< ID of the last point in the grid.
   cell *G_lpt;      ///< Pointer to last point on grid (opposite corner of grid to FirstPt())
   cell *G_fpt;      ///< Grid Data, pointer to first grid cell.
   const int G_ndim; ///< Dimensionality of Grid; public b/c it's constant so can't be messed with.
   const int G_nvar; ///< Number of variables in state vector;  public b/c it's constant so can't be messed with.
   const int G_eqntype; ///< Which equations we are using (used in reflective BCs to get it right).
   int *G_ng;       ///< number of grid points in each direction.

  ///< Total number of grid points.
  size_t G_ncell;
  /// Total number of grid points, including boundary points.
  size_t G_ncell_all;

  double *G_range; ///< Size of domain in x,y,z-direction.
  double *G_xmin;  ///< Min value of x,y,z in domain.
  double *G_xmax;  ///< Max value of x,y,z in domain.

  int *G_irange; ///< Size of domain in x,y,z-direction (integer coordinates).
  int *G_ixmin;  ///< Min value of x,y,z in domain (integer coordinates).
  int *G_ixmax;  ///< Max value of x,y,z in domain (integer coordinates).

  double G_dx;     ///< Linear side length of (uniform, cubic, cartesian) grid cells.
  double G_dA;     ///< Area of one surface of the (uniform, cubic, cartesian) grid cells.
  double G_dV;     ///< Volume of one cubic, cartesian grid cell (same for all cells).


  ///
  /// Set cell dimensions based on grid properties.
  ///
  int setCellSize();

  ///
  /// Allocate memory for grid data.
  /// This is for allocating the initial array and putting it into
  /// whatever structure is required.  Assigns IDs to the grid
  /// points. Sets the npt and npt_all pointers to the next cell.
  ///
  int allocateGridData();

   /** \brief Assign positions to each grid point, pointers to 
    *  neighbours. */
   int assignGridStructure();
   
   // Functions Relating to Boundary Conditions.
   /** \brief Set the boundary conditions string and initialise BC_bd */
   virtual int BC_setBCtypes(string ///< typeofbc string from initial conditions file.
			     );
   /** \brief Create new boundary cells for an edge cell. */
   cell * BC_setupBCcells(cell *, ///< pointer to edge cell.
			 const enum direction,    ///< direction from edge cell to boundary cell.
			 double, ///< size of cell, dx.
			 int ///< number of cells off the edge we are assigning,
			     ///< =1 for edge cell, =2 for 1st boundary cell.
			 );
   /** \brief Assigns data to a periodic boundary. */
   virtual int BC_assign_PERIODIC(  boundary_data *);
   /** \brief Assigns data on an outflow (zero gradient) boundary. */
   int         BC_assign_OUTFLOW(   boundary_data *);
   ///
   /// Assigns data on a one-way outflow (zero gradient) boundary.
   /// If the flow is off-domain, then I use zero-gradient, but if flow
   /// is onto domain then I set the boundary cells to have zero normal 
   /// velocity.
   ///
   int         BC_assign_ONEWAY_OUT(boundary_data *);
   /** \brief Assigns data on an inflow (fixed) boundary. */
   int         BC_assign_INFLOW(    boundary_data *);
   /** \brief Assigns data on a reflecting boundary. */
   int         BC_assign_REFLECTING(boundary_data *);
   /** \brief Assigns data on a fixed boundary. */
   int         BC_assign_FIXED(     boundary_data *);
   /** \brief Sets some boundary cells to be fixed to the Jet inflow
    * condition.
    * */
   virtual int BC_assign_JETBC(     boundary_data *);
   /** \brief Assigns data on a JetReflect boundary, which is the same
    * as a reflecting boundary, except the B-field is reflected differently.
    * It is the base of a jet, so there is assumed to be another jet
    * in the opposite direction, so the normal B-field is unchanged across
    * the boundary and the tangential field is reversed.
    */
   int         BC_assign_JETREFLECT(boundary_data *);
   /** \brief Assigns data on a boundary for the Double Mach Reflection Problem. */
   int         BC_assign_DMACH(     boundary_data *);
   /** \brief Assigns data on The other DMR test problem boundary */
   virtual int BC_assign_DMACH2(    boundary_data *);
   /** \brief Assigns data for the Radiative Shock Test XN boundary. */
   virtual int BC_assign_RADSHOCK(  boundary_data *);
   /** \brief Assigns data for the alternative Radiative Shock XN boundary. */
   virtual int BC_assign_RADSH2(    boundary_data *);
   ///
   /// Add internal stellar wind boundaries -- these are (possibly time-varying)
   /// winds defined by a mass-loss-rate and a terminal velocity.  A region within
   /// the domain is given fixed values corresponding to a freely expanding wind
   /// from a cell-vertex-located source.
   ///
   int         BC_assign_STWIND(    boundary_data *);
   ///
   /// Add cells to both the Wind class, and to the boundary data list
   /// of cells.  This is re-defined for cylindrical and spherical
   /// coords below.
   ///
   virtual int BC_assign_STWIND_add_cells2src(const int,    ///< source id
					      boundary_data * ///< boundary ptr.
					      );

   /** \brief Updates data on a periodic boundary. */
   virtual int BC_update_PERIODIC(  boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on an outflow (zero gradient) boundary. */
   int         BC_update_OUTFLOW(   boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   ///
   /// Update the one-way outflow (zero gradient) boundary.
   /// If the flow is off-domain, then I use zero-gradient, but if flow
   /// is onto domain then I set the boundary cells to have zero normal 
   /// velocity.
   ///
   int         BC_update_ONEWAY_OUT(boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on inflow boundary (data fixed to initial values). */
   int         BC_update_INFLOW(    boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on reflecting boundary. */
   int         BC_update_REFLECTING(boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on fixed boundary (data doesn't change). */
   int         BC_update_FIXED(     boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on boundary where jet enters domain (keeps inflow fixed). */
   virtual int BC_update_JETBC(     boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on JetReflect boundary (see BC_assign_JETREFLECT description). */
   int         BC_update_JETREFLECT(boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on the double mach reflection (DMR) boundary. */
   int         BC_update_DMACH(     boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on the other DMR test problem boundary. */
   virtual int BC_update_DMACH2(    boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on the Radiative Shock Test Problem XN boundary. */
   virtual int BC_update_RADSHOCK(  boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   /** \brief Updates data on the alternative RadShock boundary. */
   virtual int BC_update_RADSH2(    boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );
   ///
   /// Update internal stellar wind boundaries -- these are (possibly time-varying)
   /// winds defined by a mass-loss-rate and a terminal velocity.  If fixed in time
   /// the wind is updated with b->refval, otherwise with a (slower) call to the 
   /// stellar wind class SW
   ///
   int         BC_update_STWIND(    boundary_data *, ///< Boundary to update.
				    const int,  ///< current fractional step being taken.
				    const int   ///< final step.
				    );

   /** \brief prints all the cells in the given boundary data pointer. */
   int BC_printBCdata(boundary_data *);
   /** \brief Destructor for all the boundary data, BC_bd.  Needed because 
    * we need to delete some cells in the list, and others we just need to
    * delete the pointers to them.
    * */
   void BC_deleteBoundaryData();
   struct boundary_data *BC_bd; ///< pointer to array of all boundaries.
   int BC_nbd; ///< Number of Boundaries.
   const int BC_nbc; ///< Depth of boundary cells (1 or 2 cells so far).

  public:
  /// Constructor.
  ///
  /// - Checks input parameters are reasonable.
  /// - Allocates data for the grid.
  /// - Sets up pointers to first element, and all neighbours.
  ///
  UniformGrid(
    int, ///< ndim, length of position vector.
    int, ///< nvar, length of state vectors.
    int, ///< eqntype, which equations we are using (needed by BCs).
    int, ///< number of boundary cells to use.
    double *, ///< array of minimum values of x,y,z for this grid.
    double *, ///< array of maximum values of x,y,z for this grid.
    int *,    ///< array of number of cells in x,y,z directions.
    double *, ///< array of min. x/y/z for full simulation.
    double *  ///< array of max. x/y/z for full simulation.
    );

  ///
  /// Destructor, deletes boundaries and grid data.
  ///
  ~UniformGrid();


  // ---------- QUERY BASIC GRID PROPERTIES -------------------------

  /// Returns dx (assuming cubic cells).
  virtual double DX() const {return(G_dx);}
  /// Returns dA (assuming cubic cells).
  virtual double DA() const {return(G_dA);}
  /// Returns dV (assuming cubic cells).
  virtual double DV() const {return(G_dV);}

  ///
  /// Returns cell diameter for a given cell along a given axis.
  ///
  virtual double DX(
    const cell *,   ///< cell to get dx for
    const enum axes ///< axis along which to get dx.
    ) const {return G_dx;}

  ///
  /// Returns cell boundary area for a given cell in a given
  /// direction.
  ///
  virtual double DA(
    const cell *,   ///< cell to get dA for
    const enum direction ///< direction in which to get dA.
    ) const {return G_dA;}

  ///
  /// Returns cell volume for a given cell.
  ///
  virtual double DV(const cell *) const {return G_dV;}

  /// Returns dimensionality of grid.
  virtual int Ndim() const {return(G_ndim);}
  /// Returns length of state vectors.
  virtual int Nvar() const {return(G_nvar);}


  /// Returns x/y/z lower boundary of grid in code units.
  virtual double  Xmin(enum axes a) const {return(G_xmin[a] );}
  /// Returns x/y/z upper boundary of grid in code units.
  virtual double  Xmax(enum axes a) const {return(G_xmax[a] );}
  /// Returns x/y/z range of grid in code units.
  virtual double Range(enum axes a) const {return(G_range[a]);}

  /// Returns x/y/z lower boundary of grid in integer coords (dx=2).
  virtual int  iXmin(enum axes a) const {return(G_ixmin[a] );}
  /// Returns x/y/z upper boundary of grid in integer coords (dx=2).
  virtual int  iXmax(enum axes a) const {return(G_ixmax[a] );}
  /// Returns x/y/z range of grid in integer coords (dx=2).
  virtual int iRange(enum axes a) const {return(G_irange[a]);}

  /// Returns Simulation x,y,z lower bounds in cell integer coords (dx=2)
  virtual int  SIM_iXmin(enum axes a) const {return(G_ixmin[a] );}
  /// Returns Simulation x,y,z upper bounds in cell integer coords (dx=2)
  virtual int  SIM_iXmax(enum axes a) const {return(G_ixmax[a] );}
  /// Returns Simulation x,y,z range in cell integer coords (dx=2)
  virtual int SIM_iRange(enum axes a) const {return(G_irange[a]);}

  /// Returns Simulation xyz lower bounds (code units)
  virtual int  SIM_Xmin(enum axes a) const {return(G_xmin[a] );}
  /// Returns Simulation xyz upper bounds (code units)
  virtual int  SIM_Xmax(enum axes a) const {return(G_xmax[a] );}
  /// Returns Simulation range (code units)
  virtual int SIM_Range(enum axes a) const {return(G_range[a]);}

  // ---------- QUERY BASIC GRID PROPERTIES -------------------------

  // ---------- ACCESSING AND MOVING FROM CELL TO CELL --------------
  ///
  /// Returns pointer to the first grid point (ex. boundary data).
  ///
  virtual cell * FirstPt();

  ///
  /// Returns pointer to the first grid point (inc. boundary data).
  ///
  virtual cell * FirstPt_All();

  ///
  /// Returns pointer to the last grid point (ex. boundary data).
  ///
  virtual cell * LastPt();

  ///
  /// Returns pointer to the last grid point (inc. boundary data).
  ///
  virtual cell * LastPt_All();

  ///
  /// returns a pointer to a neighbouring cell in the given
  /// direction.
  ///
  /// This function is designed to be fast, so it does no range
  /// checking on the direction, to make sure the pointer exists.  It
  /// will do something unexpected if 2*Ndim<dir, as it will access a
  /// random piece of memory.  It returns whatever is in 
  /// cell->ngb[direction] i.e. the neighbour if it exists or 0 if it
  /// doesn't.
  ///
  virtual cell * NextPt(
    const cell *c,           ///< Current cell.
    const enum direction dir ///< direction of neighbour.
    )   {return(c->ngb[dir]);}

  ///
  /// returns a pointer to the next cell (ex. boundary data).
  ///
  virtual cell * NextPt(const cell *c) {return(c->npt);}

  ///
  /// returns a pointer to the next cell (inc. boundary data).
  ///
  virtual cell * NextPt_All(const cell *c) {return(c->npt_all);}

  ///
  /// Like nextPt(cell,dir), but in reverse direction.
  ///
  virtual class cell* PrevPt(
    const class cell*c, ///< Current Point.
    enum direction dir  ///< Direction to go in.
    );
  // ---------- ACCESSING AND MOVING FROM CELL TO CELL --------------

   /** \brief Returns the opposite direction to the one passed in to the function. */
   virtual enum direction OppDir(enum direction );

   /** \brief Runs through ghost boundary cells and does the appropriate time update on them.*/
   virtual int TimeUpdateExternalBCs(const int, ///< Current step number in the timestep.
				     const int  ///< Maximum step number in timestep.
				     );
   /** \brief Runs through boundary cells which are grid cells and does the appropriate time update on them.*/
   virtual int TimeUpdateInternalBCs(const int, ///< Current step number in the timestep.
				     const int  ///< Maximum step number in timestep.
				     );
   
   /** \brief Sets up boundary data on each boundary, and any extra boundaries
    * specified in the input string.  Also assigns data to each boundary.
    * */
   virtual int SetupBCs(int,   ///< Depth of Boundary cells, 1,2,etc.
			string ///< string containing info on types of BCs on all sides.
			);
#ifdef PLLEL_RT
   /** \brief Setup lists of processors to receive data from and send data to, 
    * and setup extra boundaries at corners. */
   virtual int Setup_RT_Boundaries(const int  ///< source id
				   ) {cerr<<"DONT CALL ME!\n"; return 0;}
   /** \brief Receive all optical depths for boundaries closer to source. */
   virtual int Receive_RT_Boundaries(const int ///< source id
				     ) {cerr<<"DONT CALL ME!\n"; return 0;}
   /** \bried Send all optical depths for boundaries to domains further from source. */
   virtual int Send_RT_Boundaries(const int ///< source id
				  ) {cerr<<"DONT CALL ME!\n"; return 0;}
#endif // PLLEL_RT

   ///
   /// Calculate distance between two points, where the two position
   /// are interpreted in the appropriate geometry.
   /// This function takes input in physical units, and outputs in 
   /// physical units.
   ///
   virtual double distance(const double *, ///< position 1 (physical)
			   const double *  ///< position 2 (physical)
			   );
   
   ///
   /// Calculate distance between two points, where the two position
   /// are interpreted in the appropriate geometry.
   /// This function takes input in code integer units, and outputs in
   /// integer units (but obviously the answer is not an integer).
   ///
   virtual double idistance(const int *, ///< position 1 (integer)
			    const int *  ///< position 2 (integer)
			    );
   
   ///
   /// Calculate distance between two cell--centres (will be between
   /// centre-of-volume of cells if non-cartesian geometry).
   /// Result returned in physical units (e.g. centimetres).
   ///
   virtual double distance_cell2cell(const cell *, ///< cell 1
				     const cell *  ///< cell 2
				     );
   ///
   /// Calculate distance between two cell--centres (will be between
   /// centre-of-volume of cells if non-cartesian geometry).
   /// Result returned in grid--integer units (one cell has a diameter
   /// two units).
   ///
   virtual double idistance_cell2cell(const cell *, ///< cell 1
				      const cell *  ///< cell 2
				      );
   
   ///
   /// Calculate distance between a cell-vertex and a cell--centres
   /// (will be between centre-of-volume of cells if non-cartesian
   /// geometry).  Here both input and output are physical units.
   ///
   virtual double distance_vertex2cell(const double *, ///< vertex (physical)
				       const cell *    ///< cell
				       );

   ///
   /// Calculate distance between a cell-vertex and a cell--centres
   /// (will be between centre-of-volume of cells if non-cartesian
   /// geometry).  Here both input and output are code-integer units.
   ///
   virtual double idistance_vertex2cell(const int *, ///< vertex (integer)
					const cell * ///< cell
					);

   ///
   /// As idistance_vertex2cell(int,cell) but for a single component
   /// of the position vector, and not the absolute value.  It returns
   /// the *cell* coordinate minus the *vertex* coordinate.
   ///
   virtual double idifference_vertex2cell(const int *,  ///< vertex (integer)
					  const cell *, ///< cell
					  const axes    ///< Axis to calculate.
					  );
  ///
  /// As idifference_vertex2cell(int,cell,axis) but for the coordinate
  /// difference between two cell positions along a given axis.
  /// It returns *cell2* coordinate minus *cell1* coordinate.
  ///
  virtual double idifference_cell2cell(
              const cell *, ///< cell 1
              const cell *, ///< cell 2
              const axes    ///< Axis.
              );
   
};
  


// ##################################################################
// ##################################################################


///
/// Now we define a derive class cylindrical grid.  For now only the
/// distance--between--two--points functions will be re-defined.
///
class uniform_grid_cyl
: virtual public UniformGrid,
  virtual public VectorOps_Cyl {
 public:
  ///
  /// Constructor:
  ///  - Checks input parameters are reasonable.
  ///  - All grid setup is done in the cartesian version.
  ///
  uniform_grid_cyl(int, ///< ndim, length of position vector.
		   int, ///< nvar, length of state vectors.
		   int, ///< eqntype, which equations we are using (needed by BCs).
		   double *, ///< array of minimum values of x,y,z.
		   double *, ///< array of maximum values of x,y,z.
		   int *  ///< array of number of cells in x,y,z directions.
		   );
  ///
  /// Destructor: does nothing so far.
  ~uniform_grid_cyl();
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in physical units, and outputs in 
  /// physical units.
  ///
  virtual double distance(const double *, ///< position 1 (physical)
			  const double *  ///< position 2 (physical)
			  );
   
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(const int *, ///< position 1 (integer)
			   const int *  ///< position 2 (integer)
			   );
   
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(const cell *, ///< cell 1
				    const cell *  ///< cell 2
				    );
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has a diameter
  /// two units).
  ///
  virtual double idistance_cell2cell(const cell *, ///< cell 1
				     const cell *  ///< cell 2
				     );
   
  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(const double *, ///< vertex (physical)
				      const cell *    ///< cell
				      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(const int *, ///< vertex (integer)
				       const cell * ///< cell
				       );

   ///
   /// As idistance_vertex2cell(int,cell) but for a single component
   /// of the position vector, and not the absolute value.  It returns
   /// the *cell* coordinate minus the *vertex* coordinate.
   ///
   virtual double idifference_vertex2cell(const int *,  ///< vertex (integer)
					  const cell *, ///< cell
					  const axes    ///< Axis to calculate.
					  );
  ///
  /// As idifference_vertex2cell(int,cell,axis) but for the coordinate
  /// difference between two cell positions along a given axis.
  /// It returns *cell2* coordinate minus *cell1* coordinate.
  ///
  virtual double idifference_cell2cell(
              const cell *, ///< cell 1
              const cell *, ///< cell 2
              const axes    ///< Axis.
              );
 protected:
   ///
   /// Add cells to both the Wind class, and to the boundary data list
   /// of cells.  This is re-defined for cylindrical and spherical
   /// coords below.
   ///
   virtual int BC_assign_STWIND_add_cells2src(const int,    ///< source id
					      boundary_data * ///< boundary ptr.
					      );

   ///
   /// Returns the centre of volume of a cell (in the radial
   /// direction) in the dimensionless integer coordinate system.
   /// WILL NEED RE-DEFINING IN PARALLEL GRIDS B/C IXMIN REFERS TO THE
   /// LOCAL GRID.
   ///
   virtual double iR_cov(const cell *);
};



// ##################################################################
// ##################################################################


///
/// Now we define a derive class spherical grid.  For now only the
/// distance--between--two--points functions will be re-defined.
///
class uniform_grid_sph 
: virtual public UniformGrid,
  virtual public VectorOps_Sph {
 public:
  ///
  /// Constructor:
  ///  - Checks input parameters are reasonable.
  ///  - All grid setup is done in the cartesian version.
  ///
  uniform_grid_sph(int, ///< ndim, length of position vector.
		   int, ///< nvar, length of state vectors.
		   int, ///< eqntype, which equations we are using (needed by BCs).
		   double *, ///< array of minimum values of x,y,z.
		   double *, ///< array of maximum values of x,y,z.
		   int *  ///< array of number of cells in x,y,z directions.
		   );
  ///
  /// Destructor: does nothing so far.
  ~uniform_grid_sph();
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in physical units, and outputs in 
  /// physical units.
  ///
  virtual double distance(const double *, ///< position 1 (physical)
			  const double *  ///< position 2 (physical)
			  );
   
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(const int *, ///< position 1 (integer)
			   const int *  ///< position 2 (integer)
			   );
   
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(const cell *, ///< cell 1
				    const cell *  ///< cell 2
				    );
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has a diameter
  /// two units).
  ///
  virtual double idistance_cell2cell(const cell *, ///< cell 1
				     const cell *  ///< cell 2
				     );
   
  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(const double *, ///< vertex (physical)
				      const cell *    ///< cell
				      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(const int *, ///< vertex (integer)
				       const cell * ///< cell
				       );

   ///
   /// As idistance_vertex2cell(int,cell) but for a single component
   /// of the position vector, and not the absolute value.  It returns
   /// the *cell* coordinate minus the *vertex* coordinate.
   ///
   virtual double idifference_vertex2cell(const int *,  ///< vertex (integer)
					  const cell *, ///< cell
					  const axes    ///< Axis to calculate.
					  );
  ///
  /// As idifference_vertex2cell(int,cell,axis) but for the coordinate
  /// difference between two cell positions along a given axis.
  /// It returns *cell2* coordinate minus *cell1* coordinate.
  ///
  virtual double idifference_cell2cell(
              const cell *, ///< cell 1
              const cell *, ///< cell 2
              const axes    ///< Axis.
              );
 protected:
   ///
   /// Add cells to both the Wind class, and to the boundary data list
   /// of cells.  This is re-defined for cylindrical and spherical
   /// coords below.
   ///
   virtual int BC_assign_STWIND_add_cells2src(const int,    ///< source id
					      boundary_data * ///< boundary ptr.
					      );

   ///
   /// Returns the centre of volume of a cell (in the radial
   /// direction) in the dimensionless integer coordinate system.
   /// WILL NEED RE-DEFINING IN PARALLEL GRIDS B/C IXMIN REFERS TO THE
   /// LOCAL GRID.
   ///
   virtual double iR_cov(const cell *);
};



// ##################################################################
// ##################################################################






#endif // GRIDV2
#endif // STATIC_GRID_H
