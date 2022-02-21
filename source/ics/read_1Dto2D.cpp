/// \file read_1Dto2D.cpp
/// \author Jonathan Mackey
///
/// File for mapping a 1D simulation of a spherical blast wave onto a 2D
/// axisymmetric grid, with a spherical clump about to be hit by a blast wave.
///
/// - 2012.02.07 JM: Started on file.
///

#include "coord_sys/VectorOps.h"
#include "ics/icgen.h"
#include <array>
#include <fstream>
#include <sstream>

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

using namespace std;

// ##################################################################
// ##################################################################

IC_read_1Dto2D::IC_read_1Dto2D() {}

// ##################################################################
// ##################################################################

IC_read_1Dto2D::~IC_read_1Dto2D() {}

// ##################################################################
// ##################################################################

int IC_read_1Dto2D::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) spdlog::error("{}: {}", "null pointer to grid!", fmt::ptr(ggg));
  ICsetup_base::rp = rrp;
  if (!rp) spdlog::error("{}: {}", "null pointer to ReadParams", fmt::ptr(rp));
  string seek, str;

  seek = "1D_InputFile";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  string inputfile = str;

  //
  // set up a vector of vectors for the simulation data.
  // The outer index is the element of the state vector.
  //
  std::vector<double> radius;
  std::vector<std::vector<double> > data;
  data.resize(SimPM->nvar);

  ifstream infile;
  infile.open(inputfile.c_str());
  if (!infile.is_open()) {
    spdlog::error("Error opening file");
    return 1;
  }
  double r, x[SimPM->nvar];
  string x2, line;
  int i = 0;
  while (!infile.eof()) {
    getline(infile, line);
    // If it's not a line containing a parameter, continue and read another
    // line.
    if ((line.empty() == true) || (line.substr(0, 1) == "#")) {
      continue;
    }
    else {
      // We have found a line of data, so read it into variables.
      istringstream ff(line);
      ff >> r;
      radius.push_back(r);
      for (int v = 0; v < SimPM->nvar; v++) {
        ff >> x[v];
        data[v].push_back(x[v]);
      }
      i++;
    }
  }
  spdlog::debug("read {} lines", i);
  infile.close();

  if (SimPM->ndim == 2 && SimPM->coord_sys != COORD_CYL)
    spdlog::error("{}: {}", "Wrong coords, use cylindrical!", SimPM->coord_sys);

  //
  // Now write the data to the grid... calculate r=sqrt(R^2+Z^2)
  // and assume spherical symmetry with properties at that value of r.
  //
  cell *c = ggg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  std::array<double, MAX_NVAR> data_vals;
  for (int v = 0; v < SimPM->nvar; v++)
    data_vals[v] = 0.0;
  // 1D sims have no B field...
  int nvar = 0, mhd = 0;
  if (SimPM->eqntype == EQEUL) {
    nvar = SimPM->nvar;
    mhd  = 0;
  }
  else {
    nvar = 5 + SimPM->ntracer;
    mhd  = 1;
  }

  do {
    CI.get_dpos(*c, dpos);
    switch (SimPM->ndim) {
      case 2:
        get_data_vals(dpos, radius, data, nvar, &data_vals[0]);
        break;
      case 3:
        get_3D_data_vals(dpos, radius, data, nvar, &data_vals[0]);
        break;
      default:
        spdlog::error("{}: {}", "read1d2d: dims", SimPM->ndim);
        break;
    }

    for (int v = 0; v < 5; v++)
      c->P[v] = data_vals[v];
    for (int v = 0; v < SimPM->ntracer; v++)
      c->P[SimPM->ftr + v] = data_vals[5 + v];
    if (mhd == 1) {
      c->P[BX] = 1.0e-8;  // low field.
      c->P[BY] = 0.0;
      c->P[BZ] = 0.0;
    }
    for (int v = 0; v < SimPM->nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = ggg->NextPt(*c)) != 0);

  radius.clear();
  for (int v = 0; v < SimPM->nvar; v++)
    data[v].clear();
  data.clear();
  return 0;
}

// ##################################################################
// ##################################################################

void IC_read_1Dto2D::get_data_vals(
    std::array<double, MAX_DIM> &dpos,  ///< Cell centre
    vector<double> &radius,             ///< radius vector
    vector<vector<double> > &data,      ///< arrays of variable data.
    const int nvar,                     ///< number of variables.
    double *out  ///< array for output data values at pos.
)
{
  int imin = 0, imax = 0, len = radius.size();
  std::array<double, MAX_DIM> origin;
  double seek = 0.0, sx = 0.0;
  for (int v = 0; v < MAX_DIM; v++)
    origin[v] = 0;
  seek = gg->distance(dpos, origin);
  // bracket this value in the radius[] array.
  while (radius[imax] < seek && imax < len - 1) {
    imax++;
  }
  imin = max(0, imax - 1);
  if (imin < 0 || imax > len - 1)
    spdlog::error("{}: {}", "position out of range.", imax);
  //
  // Now linearly interpolate to get the correct value.
  //
  if (imax == len - 1)
    sx = 1.0;  // zero slope for extrapolation
  else if (imax == 0)
    sx = 1.0;  // hope this is inside the boundary
  else
    sx = (seek - radius[imin]) / (radius[imax] - radius[imin]);
  for (int v = 0; v < nvar; v++) {
    out[v] = data[v][imin] + (data[v][imax] - data[v][imin]) * sx;
  }
  // rotate velocities of clump to point in radial direction.
  // assumes VX=radial velocity, VY=0, VZ=0
  out[VZ] = 0.0;
  out[VY] = out[VX] * dpos[Rcyl] / seek;
  out[VX] = out[VX] * dpos[Zcyl] / seek;
  return;
}

// ##################################################################
// ##################################################################

void IC_read_1Dto2D::get_3D_data_vals(
    std::array<double, MAX_DIM> &dpos,  ///< Cell centre
    vector<double> &radius,             ///< radius vector
    vector<vector<double> > &data,      ///< arrays of variable data.
    const int nvar,                     ///< number of variables.
    double *out  ///< array for output data values at pos.
)
{
  int imin = 0, imax = 0, len = radius.size();
  std::array<double, MAX_DIM> origin;
  double seek = 0.0, sx = 0.0;
  for (int v = 0; v < MAX_DIM; v++)
    origin[v] = 0;
  seek = gg->distance(dpos, origin);
  // bracket this value in the radius[] array.
  while (radius[imax] < seek && imax < len - 1) {
    imax++;
  }
  imin = max(0, imax - 1);
  if (imin < 0 || imax > len - 1)
    spdlog::error("{}: {}", "position out of range.", imax);
  //
  // Now linearly interpolate to get the correct value.
  //
  if (imax == len - 1)
    sx = 1.0;  // zero slope for extrapolation
  else if (imax == 0)
    sx = 1.0;  // hope this is inside the boundary
  else
    sx = (seek - radius[imin]) / (radius[imax] - radius[imin]);
  for (int v = 0; v < nvar; v++) {
    out[v] = data[v][imin] + (data[v][imax] - data[v][imin]) * sx;
  }
  // rotate velocities to point in radial direction.
  // assumes VX=radial velocity, VY=0, VZ=0
  out[VZ] = out[VX] * dpos[ZZ] / seek;
  out[VY] = out[VX] * dpos[YY] / seek;
  out[VX] = out[VX] * dpos[XX] / seek;
  return;
}

// ##################################################################
// ##################################################################
