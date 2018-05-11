/// \file nested_grid.h
///
/// \brief Declares Nested Grid Class (serial code).
///
/// \author Jonathan Mackey
///
/// Description:\n
/// This sets up a grid that is part of a hierarchy of nested grids.
/// It inherits from UniformGrid, adding on functions to point to the
/// parent and child grids in the hierarchy, and functions for
/// communicating boundary data between levels.  Cells that are not
/// leaf cells (i.e. there is a cell at the same position on a finer
/// grid) are labelled boundary data and are then only updated by
/// communication of data from the more refined level.
///
/// Modifications:\n
/// 2018.05.03 JM: Started on the class

#ifndef NESTED_GRID_H
#define NESTED_GRID_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"



#include <list>
using namespace std;

#include "sim_constants.h"
#include "sim_params.h"


#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "grid/stellar_wind_BC.h"
#include "grid/stellar_wind_angle.h"
#include "coord_sys/VectorOps.h"
#include "coord_sys/VectorOps_spherical.h"





// ##################################################################
// ##################################################################

///
/// Nested finite-volume grid class.
/// 
/// Can be 1,2, or 3 dimensional.  The grid cells must be cubes, so the number
/// of cells in each dimension should be in the right proportion to the 
/// length of the grid in each dimension.
///
class nested_grid 
: virtual public UniformGrid,
  virtual public GridBaseClass,
  virtual public VectorOps_Cart
{



  public:
  /// 
  /// Constructor.
  ///  - Calls the UniformGrid constructor.
  ///
  nested_grid(
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
  virtual ~nested_grid();
   
  ///
  /// Runs through ghost boundary cells and does the appropriate
  /// time update on them.  This is different from the UniformGrid
  /// because finer levels get data from coarser levels, rather than
  /// implementing the outer boundary conditions.
  ///
  //virtual int TimeUpdateExternalBCs(
  //      const double,   ///< current simulation time
  //      const int, ///< Current step number in the timestep.
  //      const int  ///< Maximum step number in timestep.
  //      );

  ///
  /// Runs through boundary cells which are grid cells and does
  /// the appropriate time update on them.  This is different from
  /// the UniformGrid because some on-grid cells are not leaf cells
  /// and so are labelled "boundary data" and are updated by getting
  /// data from the finer level grid.
  ///
  //virtual int TimeUpdateInternalBCs(
  //      const double,   ///< current simulation time
  //      const int, ///< Current step number in the timestep.
  //      const int  ///< Maximum step number in timestep.
  //      );
   
  /// Set pointer to parent grid.
  //void set_parent_grid(
  //    class GridBaseClass *  ///< pointer to parent grid.
  //    );

  /// Set pointer to child grid.
  //void set_child_grid(
  //    class GridBaseClass *  ///< pointer to child grid.
  //    );

  ///
  /// Calls UniformGrid SetupBCs function and then adds nested-grid
  /// boundaries to internal data for which finer-level data exists,
  /// and external data for which coarser-level data exists.
  ///
  //virtual int BC_setBCtypes(
  //      class SimParams &  ///< reference to SimParams list.
  //      );



  protected:

  ///
  /// Set grid cells that are not leaves (i.e. there is a child grid
  /// underneath them with finer-scale data) to be boundary data so
  /// that they are not updated.  Called by BC_setBCtypes().
  ///
  void set_nonleaf_cells_as_BD();

  ///
  /// External boundaries on refined levels may get data from the 
  /// coarser level above them, not from the pre-defined external
  /// boundary conditions.  This function checks this and resets
  /// the boundaries appropriately.  Called by BC_setBCtypes().
  ///
  void set_external_bcs_from_parent();

  ///
  /// Assigns data to each boundary.  Called by SetupBCs().  This
  /// version calls the UniformGrid version and then checks for 
  /// nested grid boundaries and assigns data for them if present.
  ///
  virtual int assign_boundary_data(
      const double,   ///< current simulation time (for DMACH)
      const double, ///< Simulation start time.
      const double,  ///< Simulation finish time.
      const double ///< minimum temperature allowed
      );

  /// Assigns data to a nested grid from finer grid.
  virtual int BC_assign_NEST_FINE( boundary_data *);

  /// Assigns data to an external boundary from coarser grid.
  virtual int BC_assign_NEST_COARSE( boundary_data *);


  /// Updates data to a nested grid from finer grid.
  virtual int BC_update_NEST_FINE( boundary_data *);

  /// Updates data to an external boundary from coarser grid.
  virtual int BC_update_NEST_COARSE(
      boundary_data *, ///< pointer to boudary struct
      const int,  ///< current fractional step being taken.
      const int   ///< full step.
      );

};
  



// ##################################################################
// ##################################################################



///
/// This is a nested cylindrical grid.  It derives from nested_grid
/// and from uniform_grid_cyl, and is a hybrid of the two.  It
/// doesn't need any new functions.
///
class nested_grid_cyl
: virtual public nested_grid,
  virtual public uniform_grid_cyl,
  virtual public UniformGrid,
  virtual public VectorOps_Cyl {
 public:
  nested_grid_cyl(
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
  /// Destructor: does nothing so far.
  ~nested_grid_cyl();
};


// ##################################################################
// ##################################################################


///
/// This is a nested cylindrical grid.  It derives from nested_grid
/// and from uniform_grid_cyl, and is a hybrid of the two.  It
/// doesn't need any new functions.
///
class nested_grid_sph
: virtual public nested_grid,
  virtual public uniform_grid_sph,
  virtual public UniformGrid,
  virtual public VectorOps_Sph 
{
 public:
  nested_grid_sph(
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
  /// Destructor: does nothing so far.
  ~nested_grid_sph();
};



#endif // NESTED_GRID_H
