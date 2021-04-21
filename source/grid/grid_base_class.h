///
/// \file grid_base_class.h
/// \author Jonathan Mackey
///
/// This file declares the abstract base class for the computational
/// grid (moved out from global.h).
///
/// Modifications:
/// - 2013.01.14 JM: Moved from global.h and tidied comments.
///    Added new functions for NextPt_All() and other functions
///    related to the npt_all linked list for grid+boundary data.
/// - 2015.01.08 JM: Moved old grid (v1) definition from global.h to
///    this file.
/// - 2015.07.06 JM: got rid of GEOMETRIC_GRID ifdef
/// - 2017.12.09 JM: updated function args for boundary data.
/// - 2018.05.10 JM: moved boundary functions from grid to setup-grid
///    and sim-control.
/// - 2018.07.30 JM: added interface flux calculation data structures

#ifndef GRID_BASE_CLASS_H
#define GRID_BASE_CLASS_H

//
// These tells code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "constants.h"
#include "grid/cell_interface.h"
#include "grid/stellar_wind_BC.h"
#include "raytracing/raytracer_base.h"

#include <list>
#include <string>
#include <vector>

///
/// Abstract Base Class to define interface for the Grid data and
/// member functions.
///
class GridBaseClass {
public:
  ///
  /// virtual destructor (to stop compiler complaining).
  ///
  virtual ~GridBaseClass() {}

  // ---------- ACCESSING AND MOVING FROM CELL TO CELL --------------
  ///
  /// Returns pointer to the first grid point (ex. boundary data).
  ///
  virtual cell *FirstPt() = 0;

  ///
  /// Returns pointer to the first grid point (inc. boundary data).
  ///
  virtual cell *FirstPt_All() = 0;

  ///
  /// Returns pointer to the last grid point (ex. boundary data).
  ///
  virtual cell *LastPt() = 0;

  ///
  /// Returns pointer to the last grid point (inc. boundary data).
  ///
  virtual cell *LastPt_All() = 0;

  ///
  /// returns a pointer to a neighbouring cell in the given direction.
  ///
  virtual cell *NextPt(
      const cell *,         ///< Current cell.
      const enum direction  ///< direction of neighbour.
      ) = 0;

  ///
  /// returns a pointer to the next cell (ex. boundary data).
  ///
  virtual cell *NextPt(const cell *) = 0;

  ///
  /// returns a pointer to the next cell (inc. boundary data).
  ///
  virtual cell *NextPt_All(const cell *) = 0;

  ///
  /// Like nextPt(cell,dir), but in reverse direction.
  ///
  virtual class cell *PrevPt(
      const class cell *,  ///< Current Point.
      enum direction       ///< Direction to go in.
      ) = 0;
  // ---------- ACCESSING AND MOVING FROM CELL TO CELL --------------

  ///
  /// Returns the opposite direction to the one passed in to the
  /// function (why is this here?)
  ///
  virtual enum direction OppDir(enum direction) = 0;

  // ---------- QUERY BASIC GRID PROPERTIES -------------------------
  virtual double DX() const    = 0;  ///< Returns dx (assuming cells are cubes).
  virtual int idx() const      = 0;  ///< Returns idx.
  virtual size_t Ncell() const = 0;  ///< number of cells on grid
  virtual size_t Ncell_all() const = 0;  ///< number of grid+ghost cells

  virtual double CellVolume(
      const cell *,  ///< Cell
      const double   /// unused
      ) = 0;         ///< Returns Volume of cell.

  virtual double CellInterface(
      const cell *,     ///< Cell
      const direction,  ///< outward normal to interface.
      const double      ///< unused
      ) = 0;            ///< Returns Surface area of interface.

  /// Return depth of boundary data in given direction
  virtual int boundary_depth(enum direction) const = 0;

  ///
  /// Returns cell diameter for a given cell along a given axis.
  ///
  virtual double DX(
      const cell *,    ///< cell to get dx for
      const enum axes  ///< axis along which to get dx.
      ) const = 0;

  virtual int Ndim() const = 0;  ///< Returns dimensionality of grid.
  virtual int Nvar() const = 0;  ///< Returns length of state vectors.

  /// Returns number of grid cells along given axis, excluding
  /// boundary cells.
  virtual int NG(const enum axes  ///< axis of grid.
                 ) const = 0;

  /// Returns number of grid cells along given axis, including
  /// boundary cells.
  virtual int NG_All(const enum axes  ///< axis of grid.
                     ) const = 0;

  /// Returns x/y/z lower boundary of grid in code units.
  virtual double Xmin(enum axes) const = 0;
  /// Returns x/y/z upper boundary of grid in code units.
  virtual double Xmax(enum axes) const = 0;
  /// Returns x/y/z range of grid in code units.
  virtual double Range(enum axes) const = 0;

  /// Returns x/y/z lower boundary of grid in code units
  /// including boundary data.
  virtual double Xmin_all(enum axes) const = 0;
  /// Returns x/y/z upper boundary of grid in code units
  /// including boundary data.
  virtual double Xmax_all(enum axes) const = 0;
  /// Returns x/y/z range of grid in code units including
  /// boundary data.
  virtual double Range_all(enum axes) const = 0;

  /// Returns Simulation x,y,z lower bounds in code units.
  virtual double SIM_Xmin(enum axes) const = 0;
  /// Returns Simulation x,y,z upper bounds in code units.
  virtual double SIM_Xmax(enum axes) const = 0;
  /// Returns Simulation x,y,z range in code units.
  virtual double SIM_Range(enum axes) const = 0;

  /// Returns level x,y,z lower bounds in code units.
  virtual double level_Xmin(enum axes) const = 0;
  /// Returns level x,y,z upper bounds in code units.
  virtual double level_Xmax(enum axes) const = 0;
  /// Returns lelel x,y,z range in code units.
  virtual double level_Range(enum axes) const = 0;

  /// Returns x/y/z lower boundary of grid in integer coords.
  virtual int iXmin(enum axes) const = 0;
  /// Returns x/y/z upper boundary of grid in integer coords.
  virtual int iXmax(enum axes) const = 0;
  /// Returns x/y/z range of grid in integer coords.
  virtual int iRange(enum axes) const = 0;

  /// Returns x/y/z lower boundary of grid in integer coords
  /// including boundary data.
  virtual int iXmin_all(enum axes) const = 0;
  /// Returns x/y/z upper boundary of grid in integer coords
  /// including boundary data.
  virtual int iXmax_all(enum axes) const = 0;
  /// Returns x/y/z range of grid in integer coords including
  /// boundary data.
  virtual int iRange_all(enum axes) const = 0;

  /// Returns Simulation x,y,z lower bounds in cell integer coords
  virtual int SIM_iXmin(enum axes) const = 0;
  /// Returns Simulation x,y,z upper bounds in cell integer coords
  virtual int SIM_iXmax(enum axes) const = 0;
  /// Returns Simulation x,y,z range in cell integer coords
  virtual int SIM_iRange(enum axes) const = 0;

  /// Returns level x,y,z lower bounds in integer coords.
  virtual int level_iXmin(enum axes) const = 0;
  /// Returns level x,y,z upper bounds in integer coords.
  virtual int level_iXmax(enum axes) const = 0;
  /// Returns lelel x,y,z range in integer coords.
  virtual int level_iRange(enum axes) const = 0;

  // ---------- QUERY BASIC GRID PROPERTIES -------------------------

  // ----------- SETUP AND UPDATE BOUNDARY DATA ---------------------
  ///
  /// Add boundary cells to a grid.
  ///
  virtual int SetupBCs(
      class SimParams &  ///< List of simulation params (including BCs)
      ) = 0;

  /// array of all boundaries.
  std::vector<struct boundary_data *> BC_bd;

  class stellar_wind *Wind;  ///< stellar wind boundary condition.
  class RayTracingBase *RT;  ///< pointer to raytracing class

  ///
  /// prints all the cells in the given boundary data pointer.
  ///
  virtual int BC_printBCdata(boundary_data *) = 0;

  ///
  /// Destructor for all the boundary data, BC_bd.  Needed because
  /// we need to delete some cells in the list, and others we just need to
  /// delete the pointers to them.
  ///
  virtual void BC_deleteBoundaryData() = 0;

  ///
  /// Destructor for a boundary data struct.  Needed because
  /// we need to delete some cells in the list, and others we just need to
  /// delete the pointers to them.
  ///
  virtual void BC_deleteBoundaryData(boundary_data *) = 0;

  // ----------- SETUP AND UPDATE BOUNDARY DATA ---------------------

  ///
  /// Calculate distance between two points, where the two position
  /// are interpreted in the appropriate geometry.
  /// This function takes input in physical units, and outputs in
  /// physical units.
  ///
  virtual double distance(
      const double *,  ///< position 1 (physical)
      const double *   ///< position 2 (physical)
      ) = 0;

  ///
  /// Calculate distance between a cell-vertex and a cell-centre
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are physical units.
  ///
  virtual double distance_vertex2cell(
      const double *,  ///< vertex (physical)
      const cell *     ///< cell
      ) = 0;

  ///
  /// Calculate distance between two cell-centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in physical units (e.g. centimetres).
  ///
  virtual double distance_cell2cell(
      const cell *,  ///< cell 1
      const cell *   ///< cell 2
      ) = 0;

  ///
  /// As distance_vertex2cell(double[],cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double difference_vertex2cell(
      const double *,  ///< vertex (double)
      const cell *,    ///< cell
      const axes       ///< Axis to calculate.
      ) = 0;

  ///
  /// Calculate distance between two points, where the two position
  /// are NOT interpreted in the appropriate geometry.
  /// This function takes input in code integer units, and outputs in
  /// integer units (but obviously the answer is not an integer).
  ///
  virtual double idistance(
      const int *,  ///< position 1 (integer)
      const int *   ///< position 2 (integer)
      ) = 0;

  ///
  /// Calculate distance between two cell-centres (will be between
  /// centre-of-volume of cells if non-cartesian geometry).
  /// Result returned in grid--integer units (one cell has dx=2).
  ///
  virtual double idistance_cell2cell(
      const cell *,  ///< cell 1
      const cell *   ///< cell 2
      ) = 0;

  ///
  /// Calculate distance between a cell-vertex and a cell--centres
  /// (will be between centre-of-volume of cells if non-cartesian
  /// geometry).  Here both input and output are code-integer units.
  ///
  virtual double idistance_vertex2cell(
      const int *,  ///< vertex (integer)
      const cell *  ///< cell
      ) = 0;

  ///
  /// As idistance_vertex2cell(int,cell) but for a single component
  /// of the position vector, and not the absolute value.  It returns
  /// the *cell* coordinate minus the *vertex* coordinate.
  ///
  virtual double idifference_vertex2cell(
      const int *,   ///< vertex (integer)
      const cell *,  ///< cell
      const axes     ///< Axis to calculate.
      ) = 0;

  ///
  /// As idifference_vertex2cell(int,cell,axis) but for the coordinate
  /// difference between two cell positions along a given axis.
  /// It returns *cell2* coordinate minus *cell1* coordinate.
  ///
  virtual double idifference_cell2cell(
      const cell *,  ///< cell 1
      const cell *,  ///< cell 2
      const axes     ///< Axis.
      ) = 0;

  ///
  /// Returns true if the position is on the grid, false otherwise
  ///
  virtual bool point_on_grid(const double *  ///< position
                             ) = 0;

  /// Set pointer to child grid (for NG grids)
  virtual void set_child_grid(class GridBaseClass *  ///< pointer to child grid.
  )
  {
    return;
  }

  /// Set pointer to parent grid (for NG grids)
  virtual void set_parent_grid(
      class GridBaseClass *  ///< pointer to parent grid.
  )
  {
    return;
  }

protected:
};

#endif  // GRID_BASE_CLASS_H
