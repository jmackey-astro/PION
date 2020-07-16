
///
/// file:    angle_projection.h
/// author:  Jonathan Mackey
/// date:    2019-11-26
///
/// Description: Routines for calculating quantities along rays that
/// are not perpendicular to the grid in 2D simulations.

#ifndef ANGLE_PROJECTION_H
#define ANGLE_PROJECTION_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/cell_interface.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"
#include "../xray/xray_emission.h"
#include <vector>


///
/// pixel coords: i=x-axis (Z), j=y-axis (R)
/// Xmin = SimPM.Xmin + n_extra*DX.
/// Xmax = SimPM.Xmax + n_extra*DX.
/// 
void pixel_centre(
    double *, // Global Xmin of grid
    double ,  // grid size, dx
    int,      // n_extra: number of extra cells in Z-dir off end of grid
    int,      // i, pixel index
    int,      // j, pixel index
    double *  // pixel centre [OUTPUT]
    );


///
/// Get Z-coordinate of ray, given pixel coordinates and an R-coord.
///
double get_Z_from_R(
    double *, ///< position vector of pixel/cell
    double,   ///< R coordinate of point sought.
    double,   ///< angle of LOS with respect to +ve z-axis (radians)
    int       ///< =1 for incoming ray, -1 for outgoing.
    );

///
/// Get R-coordinate of ray, given pixel coordinates and a Z-coord.
///
double get_R_from_Z(
    double *, ///< position vector of pixel/cell
    double,   ///< Z coordinate of point sought.
    double    ///< angle of LOS with respect to +ve z-axis (radians)
    );


///
/// Get cell containing start of ray (most distant from observer),
/// given pixel coordinates and an angle.
///
cell * get_start_of_ray(
    class GridBaseClass *,  ///< pointer to grid
    double *,  ///< position vector of pixel/cell
    double     ///< angle of ray with respect to radial direction
    );


///
/// Calculate the positions where a (curved) ray enters and exits a
/// cell in the (z,R) plane for the given input parameters.
///
void get_entry_exit_points(
    class GridBaseClass *,  ///< pointer to grid
    double *, ///< position vector of pixel.
    double,   ///< angle of ray with respect to radial direction
    int,      ///< =1 for incoming ray, -1 for outgoing.
    cell *,   ///< cell we are at.
    std::vector<double> &, ///< entry point of ray into cell  (output)
    std::vector<double> &, ///< exit point of ray from cell  (output)
    cell **           ///< cell the ray exits into.              (output)
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
int generate_angle_image(
    class SimParams &,      ///< simulation parameters
    class GridBaseClass *,  ///< computational grid
    class Xray_emission &,  ///< pointer to class.
    double,     ///< angle between LOS and symmetry axis [1,89]
    int [],   ///< Number of pixels in each direction
    size_t,   ///< total number of pixels
    int,       ///< number of extra pixels w.r.t. cells.
    size_t,      ///< number of images to make
    double ** ///< pointer to the image arrays.
    );


#endif // ANGLE_PROJECTION_H
