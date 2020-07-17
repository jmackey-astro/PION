
///
/// file:    perp_projection.h
/// author:  Jonathan Mackey
/// date:    2019-11-26
///
/// Description: Routines for calculating quantities along rays that
/// are perpendicular to the grid in 2D simulations.

#ifndef PERP_PROJECTION_H
#define PERP_PROJECTION_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/cell_interface.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"
#include "../xray/xray_emission.h"
#include <vector>


///
/// Integrate all the lines of sight for data on the grid, for
/// projection perpendicular to the line-of-sight.
///
int generate_perpendicular_image(
    class SimParams &,      ///< simulation parameters
    class GridBaseClass *,  ///< computational grid
    class Xray_emission &,  ///< pointer to class.
    int [],   ///< Number of pixels in each direction
    size_t,   ///< total number of pixels
    size_t,      ///< number of images to make
    double ** ///< pointer to the image arrays.
    );


///
/// calculate projection for a column of cells in the R-direction
/// This is a driver function that calls calc_projection_column()
/// and calc_projectionRT_column() after setting up data arrays and
/// calling get_emission_absorption_data().
///
int calculate_column(
    class SimParams &,    ///< simulation parameters
    class cell *,           ///< cell at start of column.
    class Xray_emission &,  ///< Xray emission class.
    class GridBaseClass *,  ///< pointer to grid.
    int,      ///< counter of number of columns completed.
    int,      ///< Number of cells in radial direction
    int,      ///< Number of images.
    int *,    ///< Number of pixels 
    double,   ///< Cell diameter
    double ** ///< image array.
    );

///
/// Project scalar quantities onto plane of the sky.
///
double calc_projection_column(
      const double *, ///< radius array
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius
      );

///
/// Project quantities with emission and absorption onto sky.
///
double calc_projectionRT_column(
      const double *, ///< radius array.
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius.
      );

///
/// From the simulation data arrays, calculate various emission
/// quantities.  Operates on a column of data in the R-direction.
///
int get_emission_absorption_data(
    class SimParams &,    ///< simulation parameters
    class GridBaseClass *,  ///< pointer to grid.
    class cell *,           ///< cell at start of column.
    double const* const*, ///< raw data to get variable from
    const int,    ///< number of images to write
    const size_t, ///< Number of radial data elements
    class Xray_emission &,  ///< pointer to X-ray emission class.
    double **,    ///< array for emission[img][rad] data.
    double **     ///< array for absorption[img][rad] data.
    );

///
/// Add emission (or whatever projected quantity is) from current cell
/// assuming data within the cell is constant, and that the geometric
/// part of the integral is pre-calculated in "integral" variable.
/// So the integral reduces to the sum:
///   SUM_{cells} integral * emissivity(P)
///
void add_cell_emission_to_ray(
    class SimParams &,      ///< simulation parameters
    double,           ///< path-length/geometry part of integral
    double *,         ///< State vector of cell (primitive vars)
    class Xray_emission &, ///< class for getting emission
    vector<double> &  ///< resulting array of emissivities.
    );


#endif // PERP_PROJECTION_H
