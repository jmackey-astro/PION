/// \file blast_wave.cc
/// \author Jonathan Mackey
///
/// File for setting up blast wave test problems for all possible grids
/// and geometries -- 1d,2d,3d, cartesian,cylindrical,spherical, hydro,mhd.
///
///
/// - 2010-10-01 JM: Added spherically symmetric BW (1D) in polar coordinates.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include <spdlog/spdlog.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "coord_sys/VectorOps.h"
#include "ics/icgen.h"
#include "ics/icgen_base.h"
#include <fstream>
#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

IC_blastwave::IC_blastwave()
{
  ambient   = 0;
  interface = 0;
  BW_tr     = 0;
}
IC_blastwave::~IC_blastwave()
{
  if (ambient) ambient = mem.myfree(ambient);
  if (BW_tr) BW_tr = mem.myfree(BW_tr);
}

// ##################################################################
// ##################################################################

int IC_blastwave::setup_data(
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

  seek = "BW_energy";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_blastwave::bw_energy = 1.0e51;
  else
    IC_blastwave::bw_energy = atof(str.c_str());

  seek = "BW_nzones";
  str  = rp->find_parameter(seek);
  if (str == "") IC_blastwave::bw_nzones = 4;
  IC_blastwave::bw_nzones = atof(str.c_str());

  //  seek = "BWradius";
  //  str = rp->find_parameter(seek);
  //  if (str=="") spdlog::error("{}: {}", "didn't find parameter",seek);
  //  IC_blastwave::bw_rad = atof(str.c_str());

  seek = "BWpressure";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_blastwave::bw_PG = atof(str.c_str());

  seek = "BWdensity";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_blastwave::bw_RO = atof(str.c_str());

  //  seek = "BWpgRatio";
  //  str = rp->find_parameter(seek);
  //  if (str=="") spdlog::error("{}: {}", "didn't find parameter",seek);
  //  IC_blastwave::bw_PGratio = atof(str.c_str());

  //  seek = "BWroRatio";
  //  str = rp->find_parameter(seek);
  //  if (str=="") spdlog::error("{}: {}", "didn't find parameter",seek);
  //  IC_blastwave::bw_ROratio = atof(str.c_str());

  seek = "BW_blast_dens";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  IC_blastwave::bw_blastRO = atof(str.c_str());

  seek = "BWmagfieldX";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_blastwave::bw_BX = 0.0;  // could have no field.
  else
    IC_blastwave::bw_BX = atof(str.c_str());
  seek = "BWmagfieldY";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_blastwave::bw_BY = 0.0;  // could have no field.
  else
    IC_blastwave::bw_BY = atof(str.c_str());
  seek = "BWmagfieldZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    IC_blastwave::bw_BZ = 0.0;  // could have no field.
  else
    IC_blastwave::bw_BZ = atof(str.c_str());

#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  bw_BX /= sqrt(4.0 * M_PI);
  bw_BY /= sqrt(4.0 * M_PI);
  bw_BZ /= sqrt(4.0 * M_PI);
#endif

  // tracer variables
  BW_tr = mem.myalloc(BW_tr, SimPM->ntracer);
  ostringstream temp;
  string v;
  for (int t = 0; t < SimPM->ntracer; t++) {
    temp.str("");
    temp << "BW_amb_TR" << t;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      BW_tr[t] = atof(str.c_str());
    else
      BW_tr[t] = -1.0e99;
    spdlog::debug("seek={}, str={}, val={}", seek, str, BW_tr[t]);
  }

  IC_blastwave::gam = SimPM->gamma;

  IC_blastwave::eqns = SimPM->eqntype;
  if (eqns == EQEUL)
    eqns = 1;
  else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
    eqns = 2;
  else
    spdlog::error("{}: {}", "Bad equations", eqns);

  //
  // Now get optional second ambient medium params, which define an interface
  // location, and an outer ambient medium.
  //
  get_amb2_params();

  // now make sure we are to do a blast wave sim.
  string ics = rp->find_parameter("ics");

  if (ics == "")
    spdlog::error("{}: {}", "didn't get any ics to set up.", ics);
  else if (ics == "BlastWave" && SimPM->coord_sys == COORD_CRT) {
    spdlog::info("Setting up Cartesian BlastWave problem");
    err += setup_cart_bw();
  }
  else if (ics == "BlastWave" && SimPM->coord_sys == COORD_CYL) {
    spdlog::info("Setting up Cylindrical BlastWave problem");
    err += setup_cyl_bw();
  }
  else if (ics == "BlastWave" && SimPM->coord_sys == COORD_SPH) {
    spdlog::info("Setting up Spherical BlastWave problem");
    err += setup_sph_bw();
  }
  else if (ics == "BlastWave_File" && SimPM->coord_sys == COORD_SPH) {
    spdlog::info("Setting up Spherical Blastwave with input density field");
    //
    // get text file to read density, pressure etc. from.
    //
    string fname;
    seek = "BWfilename";
    str  = rp->find_parameter(seek);
    if (str == "")
      spdlog::error("{}: {}", "Need input filename to read data from", str);
    else
      fname = str;
    //
    // index of tracer in the list of tracer variables.
    //
    int tr_id = 0;
    seek      = "BW_Hplus_tracer";
    str       = rp->find_parameter(seek);
    if (str == "")
      spdlog::error("{}: {}", "Need ID of tracer variable", str);
    else
      tr_id = atoi(str.c_str());
    //
    // now set up the data...
    //
    err += setup_sph_bw_File(fname, tr_id);
  }
  /*
  else if (ics=="BlastWaveReflect" && SimPM->coord_sys = COORD_CRT) {
    err += setup_cart_symmetric_bw();
  }
  else if (ics=="BlastWaveReflect" && SimPM->coord_sys = COORD_CYl) {
    err += setup_cyl_symmetric_bw();
  }
  */
  else
    spdlog::error("{}: {}", "Don't know what Initial Condition is!", ics);

  // Add noise to data?  Smooth data?
  double noise = 0.0;
  int smooth   = 0;
  ics          = rp->find_parameter("noise");
  if (ics != "")
    noise = atof(ics.c_str());
  else
    noise = -1;
  if (isnan(noise))
    spdlog::error("{}: {}", "noise parameter is not a number", noise);
  if (noise > 0) err += AddNoise2Data(gg, *SimPM, 2, noise);

  ics = rp->find_parameter("smooth");
  if (ics != "")
    smooth = atoi(ics.c_str());
  else
    smooth = -1;
  if (isnan(smooth))
    spdlog::error("{}: {}", "Smooth parameter not a number", smooth);
  if (smooth > 0) err += SmoothData(smooth);

  return err;
}

// ##################################################################
// ##################################################################

void IC_blastwave::get_amb2_params()
{
  interface = 0;
  ambient   = 0;
  ambient   = mem.myalloc(ambient, SimPM->nvar);

  //
  // 2011.11.11 JM: Added Second ambient medium params list, for splitting
  // domain.
  //
  string seek, str;
  ostringstream temp;
  string v;
  v = "RO";
  temp.str("");
  temp << "BW_amb2_" << v;
  seek = temp.str();
  str  = rp->find_parameter(seek);
  if (str != "")
    ambient[RO] = atof(str.c_str());
  else
    ambient[RO] = -1.0e99;

  v = "PG";
  temp.str("");
  temp << "BW_amb2_" << v;
  seek = temp.str();
  str  = rp->find_parameter(seek);
  if (str != "")
    ambient[PG] = atof(str.c_str());
  else
    ambient[PG] = -1.0e99;

  v = "VX";
  temp.str("");
  temp << "BW_amb2_" << v;
  seek = temp.str();
  str  = rp->find_parameter(seek);
  if (str != "")
    ambient[VX] = atof(str.c_str());
  else
    ambient[VX] = -1.0e99;

  v = "VY";
  temp.str("");
  temp << "BW_amb2_" << v;
  seek = temp.str();
  str  = rp->find_parameter(seek);
  if (str != "")
    ambient[VY] = atof(str.c_str());
  else
    ambient[VY] = -1.0e99;

  v = "VZ";
  temp.str("");
  temp << "BW_amb2_" << v;
  seek = temp.str();
  str  = rp->find_parameter(seek);
  if (str != "")
    ambient[VZ] = atof(str.c_str());
  else
    ambient[VZ] = -1.0e99;

  if (eqns == 2) {  // mhd sim
    v = "BX";
    temp.str("");
    temp << "BW_amb2_" << v;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      ambient[BX] = atof(str.c_str());
    else
      ambient[BX] = -1.0e99;

    v = "BY";
    temp.str("");
    temp << "BW_amb2_" << v;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      ambient[BY] = atof(str.c_str());
    else
      ambient[BY] = -1.0e99;

    v = "BZ";
    temp.str("");
    temp << "BW_amb2_" << v;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      ambient[BZ] = atof(str.c_str());
    else
      ambient[BZ] = -1.0e99;

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    ambient[BX] /= sqrt(4.0 * M_PI);
    ambient[BY] /= sqrt(4.0 * M_PI);
    ambient[BZ] /= sqrt(4.0 * M_PI);
#endif

    if (SimPM->eqntype == EQGLM) ambient[SI] = 0.;
  }  // if mhd vars

  // tracer variables
  for (int t = 0; t < SimPM->ntracer; t++) {
    temp.str("");
    temp << "BW_amb2_TR" << t;
    seek = temp.str();
    str  = rp->find_parameter(seek);
    if (str != "")
      ambient[t + SimPM->ftr] = atof(str.c_str());
    else
      ambient[t + SimPM->ftr] = -1.0e99;
    spdlog::debug(
        "AMB2: seek={}, str={}, val={}", seek, str, ambient[t + SimPM->ftr]);
  }

  seek = "BW_interface";
  str  = rp->find_parameter(seek);
  if (str != "")
    interface = atof(str.c_str());
  else
    interface = 1.0e99;

  return;
}

// ##################################################################
// ##################################################################

//
// Function for reading a text file with the simulation data from the last
// output from a previous simulation (of e.g. the HII region).
//
int IC_blastwave::setup_sph_bw_File(const std::string fname, const int tr_id)
{
  int ndim = gg->Ndim();
  if (ndim != 1)
    spdlog::error("{}: {}", "Bad ndim in setup spherical BW", ndim);
  spdlog::debug(
      "Setting up a 1D spherically symmetric simulation of a blast wave with"
      "energy = {}, density ={} in {} zones, and no magnetic field, and nvar={}\n"
      "Reading input density, pressure, and y(H+) from file {}",
      bw_energy, bw_blastRO, bw_nzones, SimPM->nvar, fname);

  //
  // open input file.
  //
  ifstream infile;
  infile.open(fname.c_str());
  if (!infile.is_open()) {
    spdlog::error("Error opening file");
    return 1;
  }
  //
  // Set up data vectors
  //
  std::vector<double> radius;
  std::vector<std::vector<double> > data;
  data.resize(SimPM->nvar);
  //
  // now read in the data
  //
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
      //
      // Read in the fluid variables plus the first tracer (this is H+).
      //
      for (int v = 0; v < SimPM->nvar - SimPM->ntracer + 1; v++) {
        ff >> x[v];
        data[v].push_back(x[v]);
      }
      i++;
    }
  }
  spdlog::debug("read {} lines", i);
  infile.close();

  //
  // Now see if the number of lines of data equals the number of grid cells
  //
  if (i != SimPM->Ncell) {
    spdlog::error("nlines={}; Ncell={}", i, SimPM->Ncell);
    return 1;
  }

  //
  // Centred at [0]
  //
  std::array<double, MAX_DIM> centre;
  centre[Rsph] = 0.0;
  //
  // Set values for the region with the SN energy:
  //
  double bw_rad = gg->DX() * bw_nzones;  // radius of blast wave region
  double Pin = 3.0 * bw_energy * (gam - 1.0) / (4.0 * M_PI * pow(bw_rad, 3.0));

  //
  // Now add the primitive variables to the grid.
  //
  spdlog::info("Assigning primitive vectors");
  class cell *c = gg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  i = 0;
  do {
    // gg->Xmin(XX)="<<gg->Xmin(XX)<<"\n";
    if (radius[i] > gg->Xmin(XX) && radius[i] < gg->Xmax(XX)) {
      // Set values of primitive variables.
      c->P[RO] = data[RO][i];
      c->P[PG] = data[PG][i];
      c->P[VX] = data[VX][i];
      c->P[VY] = data[VY][i];
      c->P[VZ] = data[VZ][i];
      if (eqns == 2) {
        c->P[BX] = data[BX][i];
        c->P[BY] = data[BY][i];
        c->P[BZ] = data[BZ][i];
      }
      //
      // Now the tracer: The ion fraction of H+ is one of the tracer
      // values in Harpreet's module, but I'm not sure which one...
      //
      for (int itr = 0; itr < SimPM->ntracer; itr++) {
        if (itr == tr_id)
          c->P[SimPM->ftr + itr] = data[SimPM->ftr][i];
        else if (itr == 4)
          c->P[SimPM->ftr + itr] = data[SimPM->ftr][i] * 0.08;
        else
          c->P[SimPM->ftr + itr] = 0.0;
      }
      // This is where I set the state inside the blast radius.
      // This is exact for spherical coords.
      // The ion fraction should be 1 already, so ignore it.
      CI.get_dpos(c, dpos);
      // dist="<<gg->distance(centre,dpos)<<"  bw_rad="<<bw_rad<<"\n";
      if (gg->distance(centre, dpos) <= bw_rad) {
        c->P[PG] = Pin;
        c->P[RO] = bw_blastRO;
        spdlog::debug("Setting cell {} to internal value", c->id);
      }
      c = gg->NextPt(c);
    }
    i++;
  } while (c != 0);
  spdlog::info("Got through data successfully");
  // Data done.

  return 0;
}

// ##################################################################
// ##################################################################

int IC_blastwave::setup_sph_bw()
{
  int ndim = gg->Ndim();
  if (ndim != 1)
    spdlog::error("{}: {}", "Bad ndim in setup spherical BW", ndim);
  spdlog::debug(
      "Setting up a 1D spherically symmetric simulation of a blast wave with energy = {}, "
      "density = {} in {} zones, and no magnetic field, and nvar=",
      bw_energy, bw_blastRO, bw_nzones, SimPM->nvar);

  //
  // Centred at [0]
  //
  std::array<double, MAX_DIM> centre;
  centre[Rsph] = 0.0;
  //
  // Set values for the region with the SN energy:
  //
  double bw_rad = gg->DX() * bw_nzones;  // radius of blast wave region
  double Pin = 3.0 * bw_energy * (gam - 1.0) / (4.0 * M_PI * pow(bw_rad, 3.0));

  //
  // Data.
  //
  spdlog::info("Assigning primitive vectors");
  class cell *cpt = gg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  do {
    // Set values of primitive variables.
    cpt->P[RO] = bw_RO;
    cpt->P[PG] = bw_PG;
    cpt->P[VX] = cpt->P[VY] = cpt->P[VZ] =
        0.0;  // Stationary gas to start with.
    if (eqns == 2) {
      cpt->P[BX] = bw_BX;
      cpt->P[BY] = bw_BY;
      cpt->P[BZ] = bw_BZ;
    }
    for (int i = 0; i < SimPM->ntracer; i++) {
      cpt->P[SimPM->ftr + i] = BW_tr[i];
    }
    // This is where I set the state inside the blast radius.
    // This is exact for spherical coords.
    CI.get_dpos(cpt, dpos);
    if (gg->distance(centre, dpos) <= bw_rad) {
      cpt->P[PG] = Pin;
      cpt->P[RO] = bw_blastRO;
      // for (int i=0;i<SimPM->ntracer;i++) {
      // cpt->P[SimPM->ftr+i] = 1.0;
      //}
    }
    //
    // Now we could have an outer ambient medium, in which case we check
    // here and overwrite the data if we are at a radius greater than the
    // interface value.
    //
    // dpos="<<dpos[Rsph]<<"\n";
    if (ambient != 0 && dpos[Rsph] >= interface) {
      for (int v = 0; v < SimPM->nvar; v++)
        cpt->P[v] = ambient[v];
    }

  } while ((cpt = gg->NextPt(cpt)) != NULL);
  spdlog::info("Got through data successfully");
  // Data done.
  return (0);
}

// ##################################################################
// ##################################################################

int IC_blastwave::setup_cyl_bw()
{
  int ndim = gg->Ndim();
  if (ndim != 2)
    spdlog::error("{}: {}", "Bad ndim in setupNDSimpleBlastWave", ndim);
  spdlog::debug(
      "Setting up a 2D cylindrically symmetric simulation of a blast wave with energy = {}, "
      "density = {} in {} zones, and magnetic field BZ = {} and nvar = {}",
      bw_energy, bw_blastRO, bw_nzones, bw_BZ, SimPM->nvar);
  // Centred at (0,0)
  std::array<double, MAX_DIM> centre;
  centre[Zcyl] = 0.0;
  centre[Rcyl] = 0.0;
  //
  // Set values for the region with the SN energy:
  //
  double bw_rad = gg->DX() * bw_nzones;  // radius of blast wave region
  double Pin = 3.0 * bw_energy * (gam - 1.0) / (4.0 * M_PI * pow(bw_rad, 3.0));

  // Set up the inside_sphere class, with 100 subpoints per cell.
  class inside_sphere stest(centre.data(), bw_rad, gg->DX(), 100, ndim);
  double vfrac;

  // Data.
  spdlog::info("Assigning primitive vectors");
  class cell *cpt = gg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  do {
    // Set values of primitive variables.
    cpt->P[RO] = bw_RO;
    cpt->P[PG] = bw_PG;
    cpt->P[VX] = cpt->P[VY] = cpt->P[VZ] =
        0.0;  // Stationary gas to start with.
    if (eqns == 2) {
      cpt->P[BX] = bw_BX;
      cpt->P[BY] = bw_BY;
      cpt->P[BZ] = bw_BZ;
    }
    for (int i = 0; i < SimPM->ntracer; i++) {
      cpt->P[SimPM->ftr + i] = 1.0;
    }
    // This is where I set the state inside the blast radius.
    // This is not so exact for cyl.coords, as I assume space inside the
    // cell is cartesian!  But it's better than nothing.
    // if( (vfrac=stest.volumeFraction(cpt)) >0) {
    CI.get_dpos(cpt, dpos);
    if (gg->distance(centre, dpos) <= bw_rad) {
      vfrac      = 1.0;
      cpt->P[PG] = vfrac * (Pin) + (1. - vfrac) * cpt->P[PG];
      cpt->P[RO] = vfrac * (bw_blastRO) + (1. - vfrac) * cpt->P[RO];
      for (int i = 0; i < SimPM->ntracer; i++) {
        cpt->P[SimPM->ftr + i] = -vfrac + (1. - vfrac);
      }
    }
  } while ((cpt = gg->NextPt(cpt)) != NULL);
  spdlog::info("Got through data successfully");
  // Data done.
  return (0);
}

// ##################################################################
// ##################################################################

int IC_blastwave::setup_cart_bw()
{
  int ndim = gg->Ndim();
  if (ndim < 2 || ndim > 3)
    spdlog::error("{}: {}", "Bad ndim in setupNDSimpleBlastWave", ndim);

  double centre[ndim];
  for (int i = 0; i < ndim; i++) {
    centre[i] = 0.0;
    // centre[i] = (SimPM->Xmax[i]-SimPM->Xmin[i])/2.;
  }
  //
  // Set values for the region with the SN energy:
  //
  double bw_rad = gg->DX() * bw_nzones;  // radius of blast wave region
  double Pin    = 0.0;
  if (ndim == 1)
    Pin = bw_energy;
  else if (ndim == 2)
    Pin = bw_energy * (gam - 1.0) / (M_PI * pow(bw_rad, 2.0));
  else if (ndim == 3)
    Pin = 3.0 * bw_energy * (gam - 1.0) / (4.0 * M_PI * pow(bw_rad, 3.0));
  else
    spdlog::error("{}: {}", "bad ndim", ndim);

  // Set up the inside_sphere class, with 100 subpoints per cell.
  int nsub = 0;
  if (ndim == 2)
    nsub = 100;
  else
    nsub = 32;
  class inside_sphere stest(centre, bw_rad, gg->DX(), nsub, ndim);
  double vfrac;

  // Data.
  spdlog::info("Assigning primitive vectors");
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    cpt->P[RO] = bw_RO;
    cpt->P[PG] = bw_PG;
    cpt->P[VX] = cpt->P[VY] = cpt->P[VZ] = 0.;  // Stationary gas to start with.
    if (eqns == 2) {
      cpt->P[BX] = bw_BX;
      cpt->P[BY] = bw_BY;
      cpt->P[BZ] = bw_BZ;
    }
    for (int i = 0; i < SimPM->ntracer; i++) {
      cpt->P[SimPM->ftr + i] = 1.0;
    }
    // This is where I set the state inside the blast radius.
    if ((vfrac = stest.volumeFraction(cpt)) > 0) {
      cpt->P[PG] = vfrac * (Pin) + (1.0 - vfrac) * cpt->P[PG];
      cpt->P[RO] = vfrac * (bw_blastRO) + (1.0 - vfrac) * cpt->P[RO];
      for (int i = 0; i < SimPM->ntracer; i++) {
        cpt->P[SimPM->ftr + i] = -vfrac + (1.0 - vfrac);
      }
    }
  } while ((cpt = gg->NextPt(cpt)) != NULL);
  spdlog::info("Got through data successfully");
  // Data done.
  return (0);
}

// ##################################################################
// ##################################################################
