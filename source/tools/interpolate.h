///
/// \file interpolate.h
///
/// \brief General purpose class for interpolating tables.
/// \author Jonathan Mackey
///
/// Modifications:
/// - 2015.03.04 JM: moved from global.h GeneralStuff class.
/// - 2015.03.23 JM: added bilinear interpolation function.

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include <vector>
using namespace std;

///
/// Interpolation class, including cubic spline (adapted from Numerical
/// Recipes) and bisection with linear interpolation.
///
class interpolate_arrays {
  public:
    interpolate_arrays();
    ~interpolate_arrays();

    ///
    /// Sets up Cubic Spline interpolation
    /// Use GSL interpolation routines
    ///
    void spline(
        const double*,  ///< Array of x values.
        const double*,  ///< Array of y values.
        const int,      ///< Length of arrays.
        double,         ///< First Derivative of interpolating function at x[1]
                        ///< (unit offset array) (>1.e30 for natural spline)
        double,         ///< First Derivative of interpolating function at x[n]
                        ///< (unit offset array) (>1.e30 for natural spline)
        int&            ///< reference id for this spline interpolation
    );

    ///
    /// Performs cubic spline interpolation to get y(x)
    /// Use GSL interpolation routines
    ///
    void splint(
        const double*,  ///< Array of x values.
        const double*,  ///< Array of y values.
        const int,      ///< reference id for this spline interpolation
        const int,      ///< nspl
        const double,   ///< x we are searching for.
        double*         ///< pointer to result.
    );

    ///
    /// Given a vector of x-values, and corresponding y-values, and an input
    /// x value, find the corresponding y-value by bisection and then linear
    /// interopolation.
    ///
    void root_find_linear(
        const double*,  ///< Array of x values.
        const double*,  ///< Array of y values.
        const size_t,   ///< Array sizes
        const double,   ///< x we are searching for.
        double*         ///< pointer to result.
    );

    ///
    /// This brackets an (x,y) value with 2 function values in each
    /// direction, and does bilinear interpolation on them.
    ///
    void root_find_bilinear(
        const double*,  ///< Array of x values.
        const double*,  ///< Array of y values.
        double**,       ///< Array of function values
        const size_t*,  ///< Array sizes
        const double*,  ///< (x,y) we are searching for.
        double*         ///< pointer to result.
    );

    ///
    /// Given a vector of x-values, and corresponding y-values, and an input
    /// x value, find the corresponding y-value by bisection and then linear
    /// interopolation.
    ///
    void root_find_linear_vec(
        const vector<double>&,  ///< Array of x values.
        const vector<double>&,  ///< Array of y values.
        const double,           ///< x we are searching for.
        double&                 ///< pointer to result.
    );

    ///
    /// This brackets an (x,y) value with 2 function values in each
    /// direction, and does bilinear interpolation on them.  Same as
    /// root_find_bilinear(), but uses STL vector objects.
    ///
    double root_find_bilinear_vec(
        const vector<double>&,          ///< Array of x values.
        const vector<double>&,          ///< Array of y values.
        const vector<vector<double>>&,  ///< Array of function values
        const vector<size_t>&,          ///< Array sizes
        const vector<double>&           ///< (x,y) we are searching for.
    );

    //
    // Function to find the value of f(x,y,z) for an inputted (x,y,z),
    // by interpolating vectors of x, y, and z.
    // See: spie.org/samples/PM159.pdf - 9.2.2 Trilinear interpolation
    //
    double root_find_trilinear_vec(
        const vector<double>&,                  ///< Array of x values
        const vector<double>&,                  ///< Array of y values
        const vector<double>&,                  ///< Array of z values
        const vector<vector<vector<double>>>&,  ///< Array of function values
                                                ///< f(x,y,z)
        const vector<size_t>&,                  ///< Array sizes
        const vector<double>&                   ///< (x,y,z) we want f for
    );

  protected:
    std::vector<gsl_spline*> slist;  ///< list of spline interpolations.
};

///
/// global instance of class, defined in interpolate.cpp.
///
extern class interpolate_arrays interpolate;

#endif  // INTERPOLATE_H
