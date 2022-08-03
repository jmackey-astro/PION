/// \file IC_3D_ShockCloud.cc
/// \author Jonathan Mackey
///
/// File for mapping a 1D simulation of a spherical blast wave onto a
/// 3D Cartesian grid, with a ~spherical clump about to be hit by a
/// blast wave.
///
/// - 2012.10.22 JM: Started on file.
/// - 2012.11.20 JM: Got it working.

#include "coord_sys/VectorOps.h"
#include "ics/icgen.h"
#include <sstream>
using namespace std;
#ifndef EXCLUDE_HD_MODULE

IC_HD_3D_ShockCloud::IC_HD_3D_ShockCloud() {}
IC_HD_3D_ShockCloud::~IC_HD_3D_ShockCloud() {}

int IC_HD_3D_ShockCloud::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) spdlog::error("{}: {}", "null pointer to grid!", fmt::ptr(ggg));

  ICsetup_base::rp = rrp;
  if (!rp) spdlog::error("{}: {}", "null pointer to ReadParams", fmt::ptr(rp));
  ;

  string seek, str;

  //
  // First get the radii of interest in the simulation:
  //  - Rmin is the radius in the 1D model where we start the 3D grid.
  //  - Rmax is the ending radius in the 1D model.
  //    Note that Rmax-Rmin must correspond to Xmax-Xmin in the 3D pfile.
  //  - Rclump is the position of the clump centre.
  //  - DRclump is the (initial) radius of the clump.
  //  - Rshock is the max. radius of the shock's impact.
  //
  double Rmin, Rmax, Rclump, Rshock, DRclump;

  seek = "HD_SC3D_Rmin";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  Rmin = atof(str.c_str()) * GS.parsec();

  seek = "HD_SC3D_Rmax";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  Rmax = atof(str.c_str()) * GS.parsec();

  seek = "HD_SC3D_Rclump";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  Rclump = atof(str.c_str()) * GS.parsec();

  seek = "HD_SC3D_Rshock";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  Rshock = atof(str.c_str()) * GS.parsec();

  seek = "HD_SC3D_DRclump";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  DRclump = atof(str.c_str()) * GS.parsec();

  //
  // Make sure shock is at least 2 clump radii away from clump!
  // also do some other paranoid checks.
  //
  if (Rclump - Rshock < 2.0 * DRclump)
    spdlog::error(
        "{}: {}", "Shock is too close to clump!",
        (Rclump - Rshock) / (DRclump));
  if (Rclump < Rshock)
    spdlog::error("{}: {}", "Clump is behind shock", Rclump - Rshock);
  if (Rshock < Rmin)
    spdlog::error("{}: {}", "Shock is at smaller radius than Rmin", Rmin);
  if (Rclump > Rmax)
    spdlog::error("{}: {}", "Clump is at larger radius than Rmax", Rmax);

  //
  // Now read data from the input file.
  //
  seek = "HD_SC3D_InputFile";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  string inputfile = str;

  //
  // set up a vector of vectors for the simulation data.
  // The outer index is the element of the state vector.
  //
  std::vector<double> radius;
  std::vector<std::vector<double> > data;
  data.resize(SimPM.nvar);

  ifstream infile;
  infile.open(inputfile.c_str());
  if (!infile.is_open()) {
    spdlog::error("Error opening file");
    return 1;
  }
  double r, x[SimPM.nvar];
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
      // cout <<r<< "\t";
      radius.push_back(r);
      for (int v = 0; v < SimPM.nvar; v++) {
        ff >> x[v];
        // cout <<x[v]<<"  ";
        data[v].push_back(x[v]);
      }
      // cout <<"\n";
      i++;
    }
  }
  cout << "read " << i << " lines.\n";
  infile.close();

  cout << "Rmin/Rmax=" << Rmin << "/" << Rmax << ", Rclump=" << Rclump
       << ", Rshock=" << Rshock << "\n";

  //
  // Now delete all the data outside of Rmin-->Rmax, but keep radius[0] and
  // radius[N-1] outside the simulation domain, so that we can still
  // interpolate.
  //
  while (radius[1] < Rmin && radius.size() > 1) {
    radius.erase(radius.begin());
    for (int v = 0; v < SimPM.nvar; v++)
      data[v].erase(data[v].begin());
  }
  if (radius.size() == 2)
    spdlog::error("{}: {}", "Need radii greater than Rmin", Rmin);

  size_t sz = radius.size();
  while (radius[sz - 2] > Rmax && sz > 2) {
    radius.pop_back();
    for (int v = 0; v < SimPM.nvar; v++)
      data[v].pop_back();
    sz = radius.size();
  }
  if (radius.size() == 2)
    spdlog::error(
        "{}: {}", "Need interior data, sampling not fine enough.", Rmax);

  //
  // Reset radii so that origin is at the clump centre.
  //
  for (size_t i = 0; i < radius.size(); i++) {
    radius[i] -= Rclump;
    cout << "i=" << i << ", rad=" << radius[i] << "\n";
  }

  //
  // Make sure Simulation range is compatible with Rmin/Rmax
  //
  if (!GS.equalD(SimPM.Range[XX], Rmax - Rmin))
    spdlog::error(
        "{}: {}", "Sim range doesn't match Rmax-Rmin", SimPM.Range[XX]);
  if (!GS.equalD(SimPM.Xmin[XX], Rmin - Rclump))
    spdlog::error("{}: {}", "Sim Rmin doesn't match Rmin", SimPM.Xmin[XX]);
  if (!GS.equalD(SimPM.Xmax[XX], Rmax - Rclump))
    spdlog::error("{}: {}", "Sim Rmax doesn't match Rmax", SimPM.Xmax[XX]);
  if (SimPM.coord_sys != COORD_CRT)
    spdlog::error(
        "{}: {}", "Wrong coordinate system, use Cartesian!", SimPM.coord_sys);

  //
  // Now write the data to the grid...  this is a bit complicated!
  // - For XX positions more than 2*DR from the clump centre we take the
  //  the spherical symmetry value and make it a planar slab at that value of
  //  z.
  // - For all other positions we calculate r=sqrt(R^2+Z^2) and assume a
  //  spherical clump with properties at that value of r.
  //
  cell *c = ggg->FirstPt();
  std::array<double, SimPM.ndim> dpos, data_vals;
  do {
    CI.get_dpos(c, dpos);
    get_data_vals(dpos, Rshock - Rclump, radius, data, SimPM.nvar, data_vals);
    for (int v = 0; v < SimPM.nvar; v++)
      c->P[v] = data_vals[v];
    for (int v = 0; v < SimPM.nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = ggg->NextPt(c)) != 0);

  radius.clear();
  for (int v = 0; v < SimPM.nvar; v++)
    data[v].clear();
  data.clear();
  return 0;
}

void IC_HD_3D_ShockCloud::get_data_vals(
    std::array<double, MAX_DIM> &dpos,  ///< Cell centre
    const double Rshock,                ///< Shock position (negative number).
    vector<double> &radius,             ///< radius vector
    vector<vector<double> > &data,      ///< arrays of variable data.
    const int nvar,                     ///< number of variables.
    double *out  ///< array for output data values at pos.
)
{
  //
  // It is assumed that Rcloud=0, and Rshock is a negative number.
  // i.e. the radial coordinate has been offset so that zero is
  // now centred on the cloud.
  //
  // cout <<"Rshock="<<Rshock<<"
  // Rmin/Rmax="<<radius.front()<<"/"<<radius.back()<<"\n";
  //
  // First see if z<Rshock, in which case just use slab symmetry.
  //
  double seek = dpos[XX];
  int imin = 0, imax = 0;
  if (seek < Rshock) {
  }
  //
  // Else we want to use circular symmetry about the origin.
  //
  else {
    double origin[MAX_DIM];
    for (int v = 0; v < MAX_DIM; v++)
      origin[v] = 0;
    seek = GS.distance(dpos, origin, SimPM.ndim);
    //
    // if distance is less than |Rshock| then can just get value
    // from radius array, but if larger, then we need to fill in the
    // blank region with the value at |Rshock|-Epsilon.
    // N.B. seek must be a negative number between Rshock and zero!
    //
    if (seek >= fabs(Rshock)) {
      seek = Rshock * ONE_MINUS_EPS;
    }
    else {
      seek = -seek;
    }
  }

  // cout <<"z="<<seek<<", min,max="<<radius.front()<<",
  // "<<radius.back()<<"\n";
  //
  // Now we have the "radius" (negative number) we are loodking for, so we
  // bracket this value in the radius[] array.
  //
  while (radius[imax] < seek) {
    imax++;
  }
  imin = imax - 1;
  if (imin < 0 || imax > (static_cast<int>(radius.size()) - 1))
    spdlog::error("{}: {}", "shock position out of range.", imax);
  //
  // Now linearly interpolate to get the correct value.
  //
  for (int v = 0; v < nvar; v++) {
    out[v] = data[v][imin]
             + (data[v][imax] - data[v][imin]) * (seek - radius[imin])
                   / (radius[imax] - radius[imin]);
  }

  //
  // Now get the bulk velocity of the clump, assuming this is the
  // same as the pre-shock velocity.
  //
  double seek2 = Rshock * ONE_MINUS_EPS;
  int ilo = 0, ihi = 0;
  while (radius[ihi] < seek2)
    ihi++;
  ilo           = ihi - 1;
  double bulk_V = data[VX][ilo]
                  + (data[VX][ihi] - data[VX][ilo]) * (seek2 - radius[ilo])
                        / (radius[ihi] - radius[ilo]);

  //
  // Now linearly interpolate Primitive variables to get the correct
  // values.
  //
  for (int v = 0; v < nvar; v++) {
    out[v] = data[v][imin]
             + (data[v][imax] - data[v][imin]) * (seek - radius[imin])
                   / (radius[imax] - radius[imin]);
  }

  if (dpos[XX] > Rshock) {
    //
    // rotate velocities of clump to point in radial direction.
    // But also maintain the bulk motion of the cloud, so that we
    // don't make it artificially collapsing.
    //
    out[VZ] =
        -(out[VX] - bulk_V) * dpos[ZZ]
        / sqrt(dpos[ZZ] * dpos[ZZ] + dpos[YY] * dpos[YY] + dpos[XX] * dpos[XX]);
    out[VY] =
        -(out[VX] - bulk_V) * dpos[YY]
        / sqrt(dpos[ZZ] * dpos[ZZ] + dpos[YY] * dpos[YY] + dpos[XX] * dpos[XX]);
    out[VX] = bulk_V
              - (out[VX] - bulk_V) * dpos[XX]
                    / sqrt(
                        dpos[ZZ] * dpos[ZZ] + dpos[YY] * dpos[YY]
                        + dpos[XX] * dpos[XX]);
  }

  return;
}

#endif  // don't EXCLUDE_HD_MODULE
