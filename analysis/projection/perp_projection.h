
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
#include "xray_emission.h"
#include <vector>

///
/// calculate projection for a column of cells in the R-direction
///
void calculate_column(
    class cell *,           ///< cell at start of column.
    class Xray_emission &,  ///< Xray emission class.
    class GridBaseClass *,  ///< pointer to grid.
    int,      ///< counter of number of columns completed.
    int,      ///< Number of cells in radial direction
    int,      ///< Number of images.
    int *,    ///< Number of pixels 
    double,   ///< Cell diameter
    size_t,   ///< H+ fraction is 1st tracer.
    size_t,   ///< Wind fraction is 2nd tracer.
    double ** ///< image array.
    );

///
/// From the simulation data arrays, calculate various emission
/// quantities.  Operates on a column of data in the R-direction.
///
int get_emission_absorption_data(
    class cell *,           ///< cell at start of column.
    double const* const*, ///< raw data to get variable from
    const int,    ///< number of images to write
    const size_t, ///< Number of radial data elements
    class Xray_emission &,  ///< pointer to X-ray emission class.
    double **,    ///< array for emission[img][rad] data.
    double **     ///< array for absorption[img][rad] data.
    );



//
// Project scalar quantities onto plane of the sky.
//
double calc_projection(
      const double *, ///< radius array
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius
      );

//
// Project quantities with emission and absorption onto sky.
//
double calc_projectionRT(
      const double *, ///< radius array.
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius.
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

//
// Integrate all the lines of sight for data on the grid.
//
int generate_perpendicular_image(
    class SimParams &,      ///< simulation parameters
    class GridBaseClass *,  ///< computational grid
    class Xray_emission &,  ///< pointer to class.
    double,     ///< angle between LOS and symmetry axis [1,89]
    int [],   ///< Number of pixels in each direction
    size_t,   ///< total number of pixels
    double ** ///< pointer to the image arrays.
    );


#endif // PERP_PROJECTION_H
