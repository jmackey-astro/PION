/// \file inside_sphere.cc
/// \author Jonathan Mackey
///
/// Test if a point is inside a sphere.
///
/// Modifications:
/// - 2015.03.11 JM: tidied up code. (also 03.26)

#include "constants.h"
#include "inside_sphere.h"

#include <spdlog/spdlog.h>

// ##################################################################
// ##################################################################

inside_sphere::inside_sphere(
    double *cen, double r, double size, int N, int nd) :
    sr(r),
    clen(size), del(size / 2.), nint(N), ndim(nd)
{
  diag = sqrt(static_cast<double>(ndim)) * clen;
  for (int i = 0; i < ndim; i++)
    spos[i] = cen[i];
}

bool inside_sphere::equalD(const double a, const double b)
{
  if (a == b) return (true);
  if (fabs(a + b) < 1.e-100) {
    spdlog::debug(
        "tiny numbers in equalD(a,b); a,b <1.e-100... a={}, b={}; returning true",
        a, b);
    return (true);
  }
  if ((fabs(a - b) / fabs(a + b + 1.e-100)) < 10. * MACHINEACCURACY)
    return (true);  // true is 1
  else
    return (false);  // false is zero.
}

// ##################################################################
// ##################################################################

double inside_sphere::distance(const double *d1, const double *d2, int nd)
{
  double temp = 0.;
  for (int i = 0; i < nd; i++)
    temp += (d1[i] - d2[i]) * (d1[i] - d2[i]);
  return (sqrt(temp));
}

// ##################################################################
// ##################################################################

double inside_sphere::volumeFraction(cell *cpt)
{
  // Assign centre of square positions.
  for (int v = 0; v < ndim; v++)
    cpos[v] = CI.get_dpos(cpt, v);
  //  printVec("Cell Centre",cpos,ndim);

  double pos[ndim];
  // First calculate nearest corner
  for (int v = 0; v < ndim; v++) {
    if ((pos[v] = cpos[v] - del - spos[v]) > 0) {
      pos[v] = pos[v] * pos[v];
    }
    else if ((pos[v] = cpos[v] + del - spos[v]) < 0) {
      pos[v] = pos[v] * pos[v];
    }
    else {
      pos[v] = 0.;
    }
  }
  double dist = 0.;
  for (int v = 0; v < ndim; v++)
    dist += pos[v];
  dist = sqrt(dist);  // distance from sphere centre to nearest edge of cell.

  // Test what this is
  if (dist > sr) {
    return (0.);
  }
  else if (dist + diag < sr) {
    return (1.);
  }

  double dv  = 1.;
  double vol = 1.;
  for (int v = 0; v < ndim; v++) {
    dv *= clen / nint;
    vol *= clen;
  }
  double dx = clen / nint;
  double startpt[ndim];
  double frac = 0.;
  // Set first position.
  int ntot = 1;
  for (int v = 0; v < ndim; v++) {
    startpt[v] = cpos[v] - del + dx / 2.;
    ntot *= nint;
  }

  for (int i = 0; i < ntot; i++) {
    pos[XX] = startpt[XX] + (i % nint) * dx;
    if (ndim > 1) pos[YY] = startpt[YY] + ((i / nint) % nint) * dx;
    if (ndim > 2) pos[ZZ] = startpt[ZZ] + ((i / nint / nint) % nint) * dx;
    if (distance(pos, spos, ndim) <= sr) {
      frac += dv;
    }
  }
  frac /= vol;
  return (frac);
}

// ##################################################################
// ##################################################################
