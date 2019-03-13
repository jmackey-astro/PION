/// \file uniform_grid.h
///
/// \brief Declares Uniform Grid Class, and parallel extension of it.
///
/// \author Jonathan Mackey
///
///
/// UniformGrid is designed to contain what every uniform finite
/// volume grid should need.  It is not intended to do any integration
/// or solving of PDEs by itself, but it should be set up by a class
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
///   setup/update functions.
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
/// - 2013.06.13 JM: Added STARBENCH1 internal boundary and functions.
/// - 2013.09.06 JM: Added difference_vertex2cell() functions.
/// - 2015.01.10 JM: New include statements for new file structure.
/// - 2015.01.28 JM: updated parallel grid declaration; tidied up the
///    tabbing and comments.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2016.02.12 JM: included static_grid.h modifications, worked on
///    new grid structure.
/// - 2016.02.19 JM: new grid structure finished, compiles and runs
///    the DMR test.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries). 03.24:fixed some bugs, redefined dir_XP
/// - 2016.05.02 JM: fixed bug: parallel grid had no SIM_Xmin()/max/
///    range() functions, so it was returning G_xmin!
/// - 2017.11.07-22 JM: updating boundary setup.
/// - 2017.12.09 JM: updated function args to get rid of SimPM references.
/// - 2018.05.10 JM: Moved boundary assignment to setup_fixed_grid,
///    and boundary updates to update_boundaries.cpp

#ifndef UNIFORM_GRID_H
#define UNIFORM_GRID_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"



#include <list>
using namespace std;

#include "sim_constants.h"
#include "sim_params.h"


#include "grid/grid_base_class.h"
#include "grid/stellar_wind_BC.h"
#include "grid/stellar_wind_angle.h"
#include "coord_sys/VectorOps.h"
#include "coord_sys/VectorOps_spherical.h"





// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################

///
/// The uniform finite-volume grid.
/// 
/// Can be 1,2, or 3 dimensional.  The grid cells must be cubes, so the number
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

  ///
  /// Pointer to last point on grid (opposite corner of grid to
  /// FirstPt()); the point with the most positive coordinates in 
  /// each direction.
  ///
  cell *G_lpt;

  ///
  /// Pointer to last point on grid (including ghost cells) (opposite
  /// corner of grid to FirstPt())
  ///
  cell *G_lpt_all;

  ///
  /// Grid Data, pointer to first grid cell; the point with the most
  /// negative coordinates in each direction.
  ///
  cell *G_fpt;

  ///
  /// Grid Data, pointer to first grid cell (including ghost cells).
  ///
  cell *G_fpt_all;

  const int G_ndim;   ///< Dimensionality of Grid
  const int G_nvar;   ///< Number of variables in state vector
  const int G_eqntype; ///< Which equations we are using.
  int G_coordsys;     ///< coordinate system (cart/cyl/sph)
  int G_ntracer;      ///< number of tracer fields.
  int G_ftr;          ///< index of first tracer field.
  int *G_ng;          ///< number of grid points in each direction.
  size_t G_ncell;   ///< Total number of grid points.
  /// number of grid points in each direction, including boundary cells
  int *G_ng_all;
  /// Total number of grid points, including boundary points.
  size_t G_ncell_all;

  double *G_range; ///< Size of domain in x,y,z-direction.
  double *G_xmin;  ///< Min value of x,y,z in domain.
  double *G_xmax;  ///< Max value of x,y,z in domain.
  double *G_range_all; ///< as G_range but including boundary data.
  double *G_xmin_all;  ///< as G_xmin but including boundary data.
  double *G_xmax_all;  ///< as G_xmax but including boundary data.

  int *G_irange; ///< Size of domain in x,y,z-direction (integer coordinates).
  int *G_ixmin;  ///< Min value of x,y,z in domain (integer coordinates).
  int *G_ixmax;  ///< Max value of x,y,z in domain (integer coordinates).
  int *G_irange_all; ///< as G_irange but including boundary data.
  int *G_ixmin_all;  ///< as G_ixmin but including boundary data.
  int *G_ixmax_all;  ///< as G_ixmax but including boundary data.

  double G_dx;  ///< Linear side length of (uniform, cube-shaped, cartesian) grid cells.
  int G_idx;  ///< diameter of grid cells in integer units.
  double G_dA;  ///< Area of one surface of the (uniform, cube-shaped, cartesian) grid cells.
  double G_dV;  ///< Volume of one cube-shaped, cartesian grid cell (same for all cells).

  double *Sim_range; ///< Size of full domain in x,y,z-direction.
  double *Sim_xmin;  ///< Min value of x,y,z in full domain.
  double *Sim_xmax;  ///< Max value of x,y,z in full domain.
  int *Sim_irange; ///< Size of full domain: int coords
  int *Sim_ixmin;  ///< Min value of x,y,z in full domain: int coords
  int *Sim_ixmax;  ///< Max value of x,y,z in full domain: int coords

  double *L_range; ///< Size of level domain in x,y,z-direction.
  double *L_xmin;  ///< Min value of x,y,z in level domain.
  double *L_xmax;  ///< Max value of x,y,z in level domain.
  int *L_irange; ///< Size of level domain: int coords
  int *L_ixmin;  ///< Min value of x,y,z in level domain: int coords
  int *L_ixmax;  ///< Max value of x,y,z in level domain: int coords

  ///
  /// Set cell dimensions based on grid properties.
  ///
  int set_cell_size();

  ///
  /// Allocate memory for grid data.
  /// This is for allocating the initial array and putting it into
  /// whatever structure is required.  Assigns IDs to the grid
  /// points. Sets the npt and npt_all pointers to the next cell.
  ///
  int allocate_grid_data();

  ///
  /// Assign positions to each grid point, pointers to neighbours.
  ///
  int assign_grid_structure();




  int BC_nbd; ///< Number of Boundaries.
  ///
  /// Depth of boundary cells (1 or 2 cells so far).
  ///
  const int BC_nbc;

  public:
  /// 
  /// Constructor.
  ///  - Checks input parameters are reasonable.
  ///  - Allocates data for the grid.
  ///  - Sets up pointers to first element, and all neighbours.
  ///
  UniformGrid(
    int, ///< ndim, length of position vector.
    int, ///< nvar, length of state vectors.
    int, ///< eqntype, which equations we are using (needed by BCs).
    int, ///< number of boundary cells to use.
    double *, ///< array of minimum values of x,y,z for this grid.
    double *, ///< array of maximum values of x,y,z for this grid.
    int *,    ///< array of number of cells in x,y,z directions.
    double *, ///< array of min. x/y/z for level.
    double *, ///< array of max. x/y/z for level.
    double *, ///< array of min. x/y/z for full simulation.
    double *  ///< array of max. x/y/z for full simulation.
    );

  ///
  /// Destructor, deletes boundaries and grid data.
  ///
  virtual ~UniformGrid();
   
  // ---------- QUERY BASIC GRID PROPERTIES -------------------------

  /// Returns the diameter of a grid cell.
  virtual double DX() const {return(G_dx);}
  size_t Ncell() const {return G_ncell;} ///< number of cells on grid
  size_t Ncell_all() const {return G_ncell_all;} ///< number of grid+ghost cells

  ///
  /// Returns cell diameter for a given cell along a given axis.
  ///
  virtual double DX(
    const cell *,   ///< cell to get dx for
    const enum axes ///< axis along which to get dx.
    ) const {return G_dx;}

  /// Returns dimensionality of grid.
  virtual int Ndim() const {return(G_ndim);}

  /// Returns length of state vectors.
  virtual int Nvar() const {return(G_nvar);}

  /// Returns number of grid cells along given axis, excluding
  /// boundary cells.
  virtual int NG(
    const enum axes ax ///< axis of grid.
    ) const {return G_ng[ax];} 
  
  /// Returns number of grid cells along given axis, including
  /// boundary cells.
  virtual int NG_All(
    const enum axes ax ///< axis of grid.
    ) const {return G_ng_all[ax];}

  /// Returns x/y/z lower boundary of grid in code units.
  virtual double  Xmin(enum axes a) const
  {return(G_xmin[a] );}

  /// Returns x/y/z upper boundary of grid in code units.
  virtual double  Xmax(enum axes a) const
  {return(G_xmax[a] );}

  /// Returns x/y/z range of grid in code units.
  virtual double Range(enum axes a) const
  {return(G_range[a]);}

  /// Returns Simulation xyz lower bounds (code units)
  virtual double SIM_Xmin(enum axes a) const
  {return(Sim_xmin[a] );}

  /// Returns Simulation xyz upper bounds (code units)
  virtual double SIM_Xmax(enum axes a) const
  {return(Sim_xmax[a] );}

  /// Returns Simulation range (code units)
  virtual double SIM_Range(enum axes a) const
  {return(Sim_range[a]);}

  /// Returns Level xyz lower bounds (code units)
  virtual double level_Xmin(enum axes a) const
  {return(L_xmin[a] );}

  /// Returns Level xyz upper bounds (code units)
  virtual double level_Xmax(enum axes a) const
  {return(L_xmax[a] );}

  /// Returns Level range (code units)
  virtual double level_Range(enum axes a) const
  {return(L_range[a]);}


  /// Returns x/y/z min of grid in integer coords.
  virtual double  Xmin_all(enum axes a) const
  {return(G_xmin_all[a] );}

  /// Returns x/y/z max of grid in integer coords.
  virtual double  Xmax_all(enum axes a) const
  {return(G_xmax_all[a] );}

  /// Returns x/y/z range of grid in integer coords.
  virtual double Range_all(enum axes a) const
  {return(G_range_all[a]);}


  /// Returns x/y/z lower boundary of grid in integer coords.
  virtual int  iXmin(enum axes a) const
  {return(G_ixmin[a] );}

  /// Returns x/y/z upper boundary of grid in integer coords.
  virtual int  iXmax(enum axes a) const
  {return(G_ixmax[a] );}

  /// Returns x/y/z range of grid in integer coords.
  virtual int iRange(enum axes a) const
  {return(G_irange[a]);}

  /// Returns x/y/z min of grid in integer coords.
  virtual int  iXmin_all(enum axes a) const
  {return(G_ixmin_all[a] );}

  /// Returns x/y/z max of grid in integer coords.
  virtual int  iXmax_all(enum axes a) const
  {return(G_ixmax_all[a] );}

  /// Returns x/y/z range of grid in integer coords.
  virtual int iRange_all(enum axes a) const
  {return(G_irange_all[a]);}

  /// Returns Simulation x,y,z lower bounds in cell integer coords
  virtual int  SIM_iXmin(enum axes a) const
  {return(Sim_ixmin[a] );}

  /// Returns Simulation x,y,z upper bounds in cell integer coords
  virtual int  SIM_iXmax(enum axes a) const
  {return(Sim_ixmax[a] );}

  /// Returns Simulation x,y,z range in cell integer coords
  virtual int SIM_iRange(enum axes a) const
  {return(Sim_irange[a]);}

  /// Returns Level x,y,z lower bounds in cell integer coords
  virtual int  level_iXmin(enum axes a) const
  {return(L_ixmin[a] );}

  /// Returns Level x,y,z upper bounds in cell integer coords
  virtual int  level_iXmax(enum axes a) const
  {return(L_ixmax[a] );}

  /// Returns Level x,y,z range in cell integer coords
  virtual int level_iRange(enum axes a) const
  {return(L_irange[a]);}

  /// Return cell size in integer units
  int idx() const
  {return G_idx;}

  virtual double CellVolume(
      const cell *c, const double
      ) {return VectorOps_Cart::CellVolume(c,G_dx);}

  virtual double CellInterface(
      const cell *c, ///< Cell
      const direction dir, ///< outward normal to interface.
      const double
      ) {return VectorOps_Cart::CellInterface(c,dir,G_dx);}

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
  /// Like NextPt(cell,dir), but in reverse direction.
  ///
  virtual class cell* PrevPt(
        const class cell *, ///< Current Point.
        enum direction      ///< Direction to go in.
        );

  // ---------- ACCESSING AND MOVING FROM CELL TO CELL --------------

  ///
  /// Returns the opposite direction to the one passed in to the function.
  ///
  virtual enum direction OppDir(enum direction );
  
  ///
  /// Sets up grid cells on each boundary, and any extra boundaries
  /// specified (internal, e.g. stellar wind).
  ///
  virtual int SetupBCs(
      class SimParams &  ///< List of simulation params (including BCs)
      );

  ///
  /// Setup the flux struct flux_update_recv with list of interfaces
  /// that need to be updated with fluxes from a finer level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual int setup_flux_recv(
      class SimParams &,  ///< simulation params (including BCs)
      const int           ///< level to receive from
      );

  ///
  /// Setup the flux struct flux_update_send with list of interfaces
  /// that need to be sent to a coarser level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual int setup_flux_send(
      class SimParams &,  ///< simulation params (including BCs)
      const int           ///< level to receive from
      );

  ///
  /// Add cells to the flux_update_recv and flux_update_send lists,
  /// for keeping track of flux entering/leaving a grid on one level.
  /// This is book-keeping for NG/AMR to ensure fluxes are
  /// consistent across all levels, following Berger & Colella
  /// (1989,JCP,82,64).
  /// 
  int add_cells_to_face(
      enum direction,  ///< which direction we're facing
      int *,   ///< xmin of interface region (integer units)
      int *,   ///< xmax of interface region (integer units)
      int *,   ///< number of elements in interface region
      const int,     ///< number of cells per face, per dim.
      struct flux_update & ///< struct with list to populate
      );

  ///
  /// Add fluxes from boundary cells to grid structures, for cells
  /// at the grid boundary, to be sent to the coarser level grid.
  /// Assumes that fluxes have already been calculated and saved in
  /// array "F" of all relevant cells.
  ///  
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual void save_fine_fluxes(
      const int,   ///< step number for this grid level
      const double ///< dt for this grid level
      );

  ///
  /// Add fluxes from internal cells to grid structures, for cells
  /// that sit above a grid boundary at a finer level.
  /// Assumes that fluxes have already been calculated and saved in
  /// array "F" of all relevant cells.
  ///
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual void save_coarse_fluxes(
      const double ///< dt for this grid level
      );

  ///
  /// prints all the cells in the given boundary data pointer.
  ///
  int BC_printBCdata(boundary_data *);

  ///
  /// Destructor for all the boundary data, BC_bd. 
  ///
  void BC_deleteBoundaryData();

  ///
  /// Destructor for a boundary data struct.
  ///
  void BC_deleteBoundaryData(boundary_data *);

  ///
  /// Setup lists of processors to receive data from and send data to, 
  /// and setup extra boundaries at corners.
  ///
  virtual int Setup_RT_Boundaries(
      const int,  ///< source id
      struct rad_src_info &
      ) {cerr<<"DONT CALL ME 1!\n"; return 0;}

  ///
  /// Receive all optical depths for boundaries closer to source.
  ///
  virtual int Receive_RT_Boundaries(
      const int, ///< source id
      struct rad_src_info &
      ) {cerr<<"DONT CALL ME 2!\n"; return 0;}

  ///
  /// Send all optical depths for boundaries to domains further from source.
  ///
  virtual int Send_RT_Boundaries(
      const int, ///< source id
      struct rad_src_info &
      ) {cerr<<"DONT CALL ME 3!\n"; return 0;}

  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in physical units, and outputs in 
  /// physical units.
  ///
  virtual double distance(
      const double *, ///< position 1 (physical)
      const double *  ///< position 2 (physical)
      );
   
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(
      const int *, ///< position 1 (integer)
      const int *  ///< position 2 (integer)
      );
   
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );

  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has a diameter
  /// two units).
  ///
  virtual double idistance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );
   
  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(
      const double *, ///< vertex (physical)
      const cell *    ///< cell
      );

  ///
  /// As distance_vertex2cell(double[],cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double difference_vertex2cell(
      const double *,  ///< vertex (double)
      const cell *, ///< cell
      const axes    ///< Axis to calculate.
      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(
      const int *, ///< vertex (integer)
      const cell * ///< cell
      );

  ///
  /// As idistance_vertex2cell(int,cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double idifference_vertex2cell(
      const int *,  ///< vertex (integer)
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

  ///
  /// Returns true if the position is on the grid, false otherwise
  ///
  bool point_on_grid(
    const double * ///< position
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
  uniform_grid_cyl(
    int, ///< ndim, length of position vector.
    int, ///< nvar, length of state vectors.
    int, ///< eqntype, which equations we are using (needed by BCs).
    int, ///< number of boundary cells to use.
    double *, ///< array of minimum values of x,y,z for this grid.
    double *, ///< array of maximum values of x,y,z for this grid.
    int *,    ///< array of number of cells in x,y,z directions.
    double *, ///< array of min. x/y/z for level.
    double *, ///< array of max. x/y/z for level.
    double *, ///< array of min. x/y/z for full simulation.
    double *  ///< array of max. x/y/z for full simulation.
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
  virtual double distance(
      const double *, ///< position 1 (physical)
      const double *  ///< position 2 (physical)
      );
   
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(
      const int *, ///< position 1 (integer)
      const int *  ///< position 2 (integer)
      );
   
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );
  
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has a diameter
  /// two units).
  ///
  virtual double idistance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );
   
  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(
      const double *, ///< vertex (physical)
      const cell *    ///< cell
      );

   ///
   /// As distance_vertex2cell(double[],cell) but for a single component
   /// of the position vector, and not the absolute value.  It returns
   /// the *cell* coordinate minus the *vertex* coordinate.
   ///
   virtual double difference_vertex2cell(
      const double *,  ///< vertex (double)
      const cell *, ///< cell
      const axes    ///< Axis to calculate.
      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(
      const int *, ///< vertex (integer)
      const cell * ///< cell
      );

  ///
  /// As idistance_vertex2cell(int,cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double idifference_vertex2cell(
      const int *,  ///< vertex (integer)
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

  virtual double CellVolume(
      const cell *c, const double
      ) {return VectorOps_Cyl::CellVolume(c,G_dx);}

  virtual double CellInterface(
      const cell *c, ///< Cell
      const direction dir, ///< outward normal to interface.
      const double
      ) {return VectorOps_Cyl::CellInterface(c,dir,G_dx);}

 protected:

  ///
  /// Returns the centre of volume of a cell (in the radial
  /// direction) in the dimensionless integer coordinate system.
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
  virtual public VectorOps_Sph
{
  public:

  ///
  /// Constructor:
  ///  - Checks input parameters are reasonable.
  ///  - All grid setup is done in the cartesian version.
  ///
  uniform_grid_sph(
    int, ///< ndim, length of position vector.
    int, ///< nvar, length of state vectors.
    int, ///< eqntype, which equations we are using (needed by BCs).
    int, ///< number of boundary cells to use.
    double *, ///< array of minimum values of x,y,z for this grid.
    double *, ///< array of maximum values of x,y,z for this grid.
    int *,    ///< array of number of cells in x,y,z directions.
    double *, ///< array of min. x/y/z for level.
    double *, ///< array of max. x/y/z for level.
    double *, ///< array of min. x/y/z for full simulation.
    double *  ///< array of max. x/y/z for full simulation.
    );

  ///
  /// Destructor:
  ///
  ~uniform_grid_sph();
  
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in physical units, and outputs in 
  /// physical units.
  ///
  virtual double distance(
      const double *, ///< position 1 (physical)
      const double *  ///< position 2 (physical)
      );
   
  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(
      const int *, ///< position 1 (integer)
      const int *  ///< position 2 (integer)
      );
   
  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );

  ///
  /// Calculate distance between two cell--centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has a diameter
  /// two units).
  ///
  virtual double idistance_cell2cell(
      const cell *, ///< cell 1
      const cell *  ///< cell 2
      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(
      const double *, ///< vertex (physical)
      const cell *    ///< cell
      );

  ///
  /// As distance_vertex2cell(double[],cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double difference_vertex2cell(
      const double *,  ///< vertex (double)
      const cell *, ///< cell
      const axes    ///< Axis to calculate.
      );

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(
      const int *, ///< vertex (integer)
      const cell * ///< cell
      );

  ///
  /// As idistance_vertex2cell(int,cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double idifference_vertex2cell(
      const int *,  ///< vertex (integer)
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

  virtual double CellVolume(
      const cell *c, const double
      ) {return VectorOps_Sph::CellVolume(c,G_dx);}

  virtual double CellInterface(
      const cell *c, ///< Cell
      const direction dir, ///< outward normal to interface.
      const double
      ) {return VectorOps_Sph::CellInterface(c,dir,G_dx);}

 protected:

  ///
  /// Returns the centre of volume of a cell (in the radial
  /// direction) in the dimensionless integer coordinate system.
  ///
  virtual double iR_cov(const cell *);
};



#endif // UNIFORM_GRID_H
