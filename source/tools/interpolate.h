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

#include <vector>
using namespace std;

///
/// Interpolation class, including cubic spline (from Numerical
/// Recipes) and bisection with linear interpolation.
///
class interpolate_arrays {
  public:
  interpolate_arrays();
  ~interpolate_arrays();

  ///
  /// Sets up Cubic Spline interpolation
  /// (from Martin White's Code, from NR)
  ///
  void spline(
        const double *, ///< Array of x values.
        const double *, ///< Array of y values.
        const int ,     ///< Length of arrays.
        double ,  ///< First Derivative of interpolating function at x[1] (unit offset array) (>1.e30 for natural spline)
        double ,  ///< First Derivative of interpolating function at x[n] (unit offset array) (>1.e30 for natural spline)
        double *  ///< Empty array to store d2y/dx2 in.
        );

  ///
  /// Performs cubic spline interpolation to get y(x)
  /// (from Martin White's Code, from NR) 
  ///
  void splint(
        const double *, ///< Array of x values.
        const double *, ///< Array of y values.
        const double *, ///< Array of d2y/dx2 values.
        const int ,     ///< nspl
        const double ,  ///< x we are searching for.
        double *  ///< pointer to result.
        );

  ///
  /// NR92 spline function for C++ STL vectors.
  ///
  void spline_vec(
        const std::vector<double> &,
        const std::vector<double> &,
        const int ,
        double ,
        double ,
        std::vector<double> &
        );

  ///
  /// NR92 splint function for C++ STL vectors.
  ///
  void splint_vec(
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const int,
        const double,
        double *
        );

  ///
  /// Given a vector of x-values, and corresponding y-values, and an input
  /// x value, find the corresponding y-value by bisection and then linear
  /// interopolation.
  ///
  void root_find_linear(
        const double *, ///< Array of x values.
        const double *, ///< Array of y values.
        const size_t,   ///< Array sizes
        const double ,  ///< x we are searching for.
        double *  ///< pointer to result.
        );

  ///
  /// This brackets an (x,y) value with 2 function values in each
  /// direction, and does bilinear interpolation on them.
  ///
  void root_find_bilinear(
        const double *,    ///< Array of x values.
        const double *,    ///< Array of y values.
        double **, ///< Array of function values
        const size_t *,   ///< Array sizes
        const double *,   ///< (x,y) we are searching for.
        double *           ///< pointer to result.
        );
};

///
/// global instance of class, defined in interpolate.cpp.
///
extern class interpolate_arrays interpolate;

#endif // INTERPOLATE_H
