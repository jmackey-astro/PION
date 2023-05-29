///
/// \file interpolate.cpp
///
/// \brief General purpose class for interpolating tables.
/// \author Jonathan Mackey
///
/// Modifications:
/// - 2015.03.04 JM: moved from global.h GeneralStuff class.
/// - 2015.03.23 JM: added bilinear interolation, added bounds
///   checking for linear/bilinear int., fixed bug in linear int.
/// - 2017.07.26 JM: added root_find_bilinear_vec() function.
/// - 2019.09.03 JM: changed spline interpolation to use GSL
/// - 2021.04.09 JM: switched GSL for Boost library.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "tools/interpolate.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

///
/// global instance of class, defined in interpolate.cpp.
///
class interpolate_arrays interpolate;

// ##################################################################
// ##################################################################

interpolate_arrays::interpolate_arrays() {}

// ##################################################################
// ##################################################################

interpolate_arrays::~interpolate_arrays()
{
  // cout <<"timers.size: "<<timers.size()<<"\n";
}



// ##################################################################
// ##################################################################



void interpolate_arrays::spline(
    const double *x,
    const double *y,
    const int n,
    double yp1,
    double ypn,
    int &id  ///< reference id for this spline interpolation.
)
{
  vector<double> vx(x, x + n);
  vector<double> vy(y, y + n);
  using boost::math::interpolators::makima;
  class makima<vector<double> > *data =
      new makima<std::vector<double> >(move(vx), move(vy), yp1, ypn);
  id = slist.size();
  slist.push_back(data);

  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::splint(
    const double[],
    const double[],
    const int id,  ///< reference id for this spline interpolation.
    const int,
    const double x,
    double *y)
{
  if (id >= static_cast<int>(slist.size()))
    spdlog::error("{}: {}", "bad splint request", id);
  using boost::math::interpolators::makima;
  class makima<vector<double> > data = *(slist[id]);
  *y                                 = data(x);
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::spline_vec(
    vector<double> x,
    vector<double> y,
    double,
    double,
    int &id  ///< reference id for this spline interpolation.
)
{
  using boost::math::interpolators::makima;
  class makima<vector<double> > *data =
      new makima<std::vector<double> >(move(x), move(y));

  id = slist.size();
  slist.push_back(data);
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::splint_vec(
    const int id,  ///< reference id for this spline interpolation.
    const double x,
    double *y)
{
  if (id >= static_cast<int>(slist.size()))
    spdlog::error("{}: {}", "bad splint request", id);
  using boost::math::interpolators::makima;
  class makima<vector<double> > data = *(slist[id]);
  *y                                 = data(x);
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::root_find_linear_vec(
    const vector<double> &xarr,  ///< Array of x values.
    const vector<double> &yarr,  ///< Array of y values.
    const double xreq,           ///< x we are searching for.
    double &yreq                 ///< pointer to result.
)
{
  //
  // Given a vector of x-values, and corresponding y-values, and an input
  // x value, find the corresponding y-value by bisection and then linear
  // interopolation.
  //
  // First we find the two x-points in the array which bracket the requested
  // x-value, with bisection.
  //
  size_t len = xarr.size();
  size_t ihi = len - 1,  // upper bracketing value
      ilo    = 0,        // lower bracketing value
      imid   = 0;        // midpoint
  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (xarr[imid] < xreq)
      ilo = imid;
    else
      ihi = imid;
  } while (ihi - ilo > 1);

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  double xval = 0.0;
  if (xreq > xarr[ihi])
    xval = xarr[ihi];
  else if (xreq < xarr[ilo])
    xval = xarr[ilo];
  else
    xval = xreq;

  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  yreq =
      yarr[ilo]
      + (yarr[ihi] - yarr[ilo]) * (xval - xarr[ilo]) / (xarr[ihi] - xarr[ilo]);
  return;
}

// ##################################################################
// ##################################################################

void interpolate_arrays::root_find_linear(
    const double *xarr,  ///< Array of x values.
    const double *yarr,  ///< Array of y values.
    const size_t len,    ///< Array sizes
    const double xreq,   ///< x we are searching for.
    double *yreq         ///< pointer to result.
)
{
  //
  // Given a vector of x-values, and corresponding y-values, and an input
  // x value, find the corresponding y-value by bisection and then linear
  // interopolation.
  //
  // First we find the two x-points in the array which bracket the requested
  // x-value, with bisection.
  //
  size_t ihi = len - 1,  // upper bracketing value
      ilo    = 0,        // lower bracketing value
      imid   = 0;        // midpoint
  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (xarr[imid] < xreq)
      ilo = imid;
    else
      ihi = imid;
  } while (ihi - ilo > 1);

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  double xval = 0.0;
  if (xreq > xarr[ihi])
    xval = xarr[ihi];
  else if (xreq < xarr[ilo])
    xval = xarr[ilo];
  else
    xval = xreq;

  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  *yreq =
      yarr[ilo]
      + (yarr[ihi] - yarr[ilo]) * (xval - xarr[ilo]) / (xarr[ihi] - xarr[ilo]);
  return;
}

// ##################################################################
// ##################################################################

void interpolate_arrays::root_find_bilinear(
    const double *x,     ///< Array of x values.
    const double *y,     ///< Array of y values.
    double **f,          ///< Array of function values
    const size_t *len,   ///< Array sizes
    const double *xreq,  ///< (x,y) we are searching for.
    double *res          ///< pointer to result.
)
{
  //
  // Given a vector of (x,y)-values, and corresponding f-values,
  // and an input (x,y) value, find the corresponding f-value by
  // bisection and then bi-linear interpolation.
  //
  // First we find the two x/y-points in the array which bracket the requested
  // x/y-value, with bisection.
  //
  size_t ihi = len[0] - 1, jhi = len[1] - 1,  // upper bracketing value
      ilo = 0, jlo = 0,                       // lower bracketing value
      imid = 0, jmid = 0;                     // midpoint
  int count   = 0;
  double xval = 0.0, yval = 0.0;

  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (x[imid] < xreq[0])
      ilo = imid;
    else
      ihi = imid;
    count++;
  } while (ihi - ilo > 1);

  count = 0;
  do {
    jmid = jlo + floor((jhi - jlo) / 2.0);
    if (y[jmid] < xreq[1])
      jlo = jmid;
    else
      jhi = jmid;
    count++;
  } while (jhi - jlo > 1);

  if (ihi - ilo != 1) {
    spdlog::error(
        "root_find_bilinear: Couldn't bracket root i: {}, {}", ihi, ilo);
  }

  if (jhi - jlo != 1) {
    spdlog::error(
        "root_find_bilinear: Couldn't bracket root j: {}, {}", jhi, jlo);
  }

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  if (xreq[0] > x[ihi])
    xval = x[ihi];
  else if (xreq[0] < x[ilo])
    xval = x[ilo];
  else
    xval = xreq[0];
  if (xreq[1] > y[jhi])
    yval = y[jhi];
  else if (xreq[1] < y[jlo])
    yval = y[jlo];
  else
    yval = xreq[1];

  //
  // Now we use bilinear interpolation to get the result
  // f(x,y)=( f(lo,lo)(xhi-x)(yhi-y)+
  //          f(hi,lo)(x-xlo)(yhi-y)+
  //          f(lo,hi)(xhi-x)(y-ylo)+
  //          f(hi,hi)(x-xlo)(y-ylo) )/(dx*dy)
  //
  *res =
      (f[ilo][jlo] * (x[ihi] - xval) * (y[jhi] - yval)
       + f[ihi][jlo] * (xval - x[ilo]) * (y[jhi] - yval)
       + f[ilo][jhi] * (x[ihi] - xval) * (yval - y[jlo])
       + f[ihi][jhi] * (xval - x[ilo]) * (yval - y[jlo]));
  *res /= ((x[ihi] - x[ilo]) * (y[jhi] - y[jlo]));
  return;
}

// ##################################################################
// ##################################################################

double interpolate_arrays::root_find_bilinear_vec(
    const vector<double> &x,           ///< Array of x values.
    const vector<double> &y,           ///< Array of y values.
    const vector<vector<double> > &f,  ///< Array of function values
    const vector<size_t> &len,         ///< Array sizes
    const vector<double> &xreq         ///< (x,y) we are searching for.
)
{
  //
  // Given a vector of (x,y)-values, and corresponding f-values,
  // and an input (x,y) value, find the corresponding f-value by
  // bisection and then bi-linear interopolation.
  //
  // First we find the two x/y-points in the array which bracket the requested
  // x/y-value, with bisection.
  //
  size_t ihi = len[0] - 1, jhi = len[1] - 1,  // upper bracketing value
      ilo = 0, jlo = 0,                       // lower bracketing value
      imid = 0, jmid = 0;                     // midpoint
  int count   = 0;
  double xval = 0.0, yval = 0.0, result = 0.0;

  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (x[imid] < xreq[0])
      ilo = imid;
    else
      ihi = imid;
    count++;
  } while (ihi - ilo > 1);
  // cout <<"count="<<count<<", ihi="<<ihi<<" and ilo="<<ilo<<", x[ihi]=";
  // cout <<x[ihi]<<" and x[ilo]="<<x[ilo]<<" and xreq="<<xreq[0];

  count = 0;
  do {
    jmid = jlo + floor((jhi - jlo) / 2.0);
    if (y[jmid] < xreq[1])
      jlo = jmid;
    else
      jhi = jmid;
    count++;
  } while (jhi - jlo > 1);
  // cout <<"count="<<count<<", jhi="<<jhi<<" and jlo="<<jlo<<", y[jhi]=";
  // cout <<y[jhi]<<" and y[jlo]="<<y[jlo]<<" and yreq="<<xreq[1]<<"\n";

  if (ihi - ilo != 1) {
    spdlog::error(
        "root_find_bilinear: Couldn't bracket root i: {}, {}", ihi, ilo);
  }

  if (jhi - jlo != 1) {
    spdlog::error(
        "root_find_bilinear: Couldn't bracket root j: {}, {}", jhi, jlo);
  }

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  if (xreq[0] > x[ihi])
    xval = x[ihi];
  else if (xreq[0] < x[ilo])
    xval = x[ilo];
  else
    xval = xreq[0];
  // cout <<", xval="<<xval<<"\n";
  if (xreq[1] > y[jhi])
    yval = y[jhi];
  else if (xreq[1] < y[jlo])
    yval = y[jlo];
  else
    yval = xreq[1];

  //
  // Now we use bilinear interpolation to get the result
  // f(x,y)=( f(lo,lo)(xhi-x)(yhi-y)+
  //          f(hi,lo)(x-xlo)(yhi-y)+
  //          f(lo,hi)(xhi-x)(y-ylo)+
  //          f(hi,hi)(x-xlo)(y-ylo) )/(dx*dy)
  //
  result =
      (f[ilo][jlo] * (x[ihi] - xval) * (y[jhi] - yval)
       + f[ihi][jlo] * (xval - x[ilo]) * (y[jhi] - yval)
       + f[ilo][jhi] * (x[ihi] - xval) * (yval - y[jlo])
       + f[ihi][jhi] * (xval - x[ilo]) * (yval - y[jlo]));
  result /= ((x[ihi] - x[ilo]) * (y[jhi] - y[jlo]));
  return result;
}

// ##################################################################
// ##################################################################

double interpolate_arrays::root_find_trilinear_vec(
    const vector<double> &x_vec,  ///< Array of x values
    const vector<double> &y_vec,  ///< Array of y values
    const vector<double> &z_vec,  ///< Array of z values
    const vector<vector<vector<double> > >
        &f,                          ///< Array of function values f(x,y,z)
    const vector<size_t> &vec_size,  ///< Array sizes
    const vector<double> &input      ///< (x,y,z) we want f for.
)
{
  //
  // Determine (x0,y0,z0) and (x1,y1,z1) - nearest neighbours to input (x,y,z)
  //

  // Initialise nearest neighbour variables and x, y, z
  double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
  double x = input[0], y = input[1], z = input[2];

  // Loops to determine the value of the nearest neighbours of x, y and z
  size_t x_index = 0, y_index = 0, z_index = 0;

  while (x > x_vec[x_index])
    x_index++;
  while (y > y_vec[y_index])
    y_index++;
  while (z > z_vec[z_index])
    z_index++;

  // Set nearest neighbours
  x0 = x_vec[x_index - 1], x1 = x_vec[x_index];
  y0 = y_vec[y_index - 1], y1 = y_vec[y_index];
  z0 = z_vec[z_index - 1], z1 = z_vec[z_index];

  if (x_index <= 0 || x_index >= vec_size[0]) {
    spdlog::error("x out of range: x_index={}, {}", x_index, vec_size[0]);
    spdlog::error("xvec : {}", x_vec);
    spdlog::error("x={}", x);
    spdlog::error("Bug {}", 1);
    exit(1);
  }
  if (y_index <= 0 || y_index >= vec_size[1]) {
    spdlog::error(
        "y out of range: y_index={}, y={}, {}", y_index, y, vec_size[1]);
    spdlog::error("yvec {}", y_vec);
    spdlog::error("Bug {}", 2);
    exit(2);
  }
  if (z_index <= 0 || z_index >= vec_size[2]) {
    spdlog::error("z out of range: z_index={}, {}", z_index, vec_size[2]);
    spdlog::error("Bug {}", 3);
    exit(3);
  }

  // Calculate delta x, delta y and delta z terms for trilinear
  // interpolation
  double dx = (x - x0) / (x1 - x0);
  double dy = (y - y0) / (y1 - y0);
  double dz = (z - z0) / (z1 - z0);

  // Calculate f000, f001 etc. terms for coefficients
  double f000 = f[x_index - 1][y_index - 1][z_index - 1];  // f(x0, y0, z0)
  double f001 = f[x_index - 1][y_index - 1][z_index];      // f(x0, y0, z1)
  double f010 = f[x_index - 1][y_index][z_index - 1];      // f(x0, y1, z0)
  double f100 = f[x_index][y_index - 1][z_index - 1];      // f(x1, y0, z0)
  double f110 = f[x_index][y_index][z_index - 1];          // f(x1, y1, z0)
  double f011 = f[x_index - 1][y_index][z_index];          // f(x0, y1, z1)
  double f101 = f[x_index][y_index - 1][z_index];          // f(x1, y0, z1)
  double f111 = f[x_index][y_index][z_index];              // f(x1, y1, z1)

  // Calculate c coefficients for trilinear interpolation
  double c0 = f000;
  double c1 = f100 - f000;
  double c2 = f010 - f000;
  double c3 = f001 - f000;
  double c4 = f110 - f010 - f100 + f000;
  double c5 = f011 - f001 - f010 + f000;
  double c6 = f101 - f001 - f100 + f000;
  double c7 = f111 - f011 - f101 - f110 + f100 + f001 + f010 - f000;

  // Return f(x,y,z)
  return c0 + c1 * dx + c2 * dy + c3 * dz + c4 * dx * dy + c5 * dy * dz
         + c6 * dz * dx + c7 * dx * dy * dz;
}

// ##################################################################
// ##################################################################
