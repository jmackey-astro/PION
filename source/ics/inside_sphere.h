/// \file inside_sphere.h
/// \author Jonathan Mackey
///
/// Tests what fraction of a cell is within a spherical region.
///
/// Modifications:
/// - 2015.03.26 JM: added include statement.

#ifndef INSIDE_SPHERE_H
#define INSIDE_SPHERE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <iostream>
//#include <cmath>
using namespace std;

#include "grid/cell_interface.h"

class inside_sphere {
  public:
    inside_sphere(
        double*,  ///< centre coordinates of circle/sphere
        double,   ///< radius of circle/sphere.
        double,   ///< size of cell (length of one side)
        int,      ///< number of subcells to use per cell dimension.
        int       ///< number of spatial dimensions on grid.
    );
    ~inside_sphere() {}            ///< trivial destructor
    double volumeFraction(cell*);  ///< calculates fraction of cell that is
                                   ///< within radius r of a point.
  private:
    const double sr;    ///< radius of sphere.
    double spos[3];     ///< centre coords of sphere.
    const double clen;  ///< Length of side of square/cubic cell.
    const double del;   ///< half the side length.
    double diag;        ///< length of cell diagonal (from centre to corner).
    double cpos[3];     ///< centre coords of cubic cell.
    const int nint;  ///< nint=number of subcells to split the cells into (per
                     ///< dimension)
    const int ndim;  ///< dimensionality of space.
    bool equalD(const double, const double);  ///< test for equality.
    double distance(
        const double*,
        const double*,
        int);  ///< distance between two points.
};

#endif  // INSIDE_SPHERE_H
