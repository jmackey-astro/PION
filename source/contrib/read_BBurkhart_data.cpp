/// \file read_BBurkhart_data.cpp
/// \author Jonathan Mackey
///
/// File for setting up initial conditions where data from a
/// turbulence simulation from Blakesley Burkhart are read in.
///
/// - 2012.10.15 JM: Started on file.
/// - 2012.10.16 JM: Debugged and got a stable version working.
/// - 2013.02.27 JM: minor tidying up.
/// - 2015.10.27 JM: updated for new pion
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/mem_manage.h"


#include "constants.h"
#include "coord_sys/VectorOps.h"
#include "dataIO/dataio.h"
#include "dataIO/dataio_fits.h"
#include "ics/icgen.h"
#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

IC_read_BBurkhart_data::IC_read_BBurkhart_data() {}

// ##################################################################
// ##################################################################

IC_read_BBurkhart_data::~IC_read_BBurkhart_data() {}

// ##################################################################
// ##################################################################

int IC_read_BBurkhart_data::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  ICsetup_base::gg = ggg;
  if (!gg) spdlog::error("{}: {}", "null pointer to grid!", fmt::ptr(ggg));

  ICsetup_base::rp = rrp;
  if (!rp) spdlog::error("{}: {}", "null pointer to ReadParams", fmt::ptr(rp));
  ;

  string seek, str;

  //
  // First get data file names and check that they exist.
  //
  seek = "BB_filename_RO";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  File_RO = str;

  //
  // There is no pressure file (because isothermal), so this search
  // must be allowed to fail.
  //
  seek = "BB_filename_PG";
  str  = rp->find_parameter(seek);
  if (str == "")
    File_PG = "NONE";
  else
    File_PG = str;

  seek = "BB_filename_VX";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  File_VX = str;

  seek = "BB_filename_VY";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  File_VY = str;

  seek = "BB_filename_VZ";
  str  = rp->find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
  File_VZ = str;

  if (SimPM.eqntype == EQMHD || SimPM.eqntype == EQGLM
      || SimPM.eqntype == EQFCD) {
    seek = "BB_filename_BX";
    str  = rp->find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
    File_BX = str;

    seek = "BB_filename_BY";
    str  = rp->find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
    File_BY = str;

    seek = "BB_filename_BZ";
    str  = rp->find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "didn't find parameter", seek);
    File_BZ = str;
  }
  else {
    File_BX = "NOFILE";
    File_BY = "NOFILE";
    File_BZ = "NOFILE";
  }

  //
  // See if there is a velocity offset to the data:
  //
  seek = "BB_offset_VX";
  str  = rp->find_parameter(seek);
  if (str == "")
    VXoffset = 0.0;
  else
    VXoffset = atof(str.c_str());
  seek = "BB_offset_VY";
  str  = rp->find_parameter(seek);
  if (str == "")
    VYoffset = 0.0;
  else
    VYoffset = atof(str.c_str());
  seek = "BB_offset_VZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    VZoffset = 0.0;
  else
    VZoffset = atof(str.c_str());
  cout << "Velocity offset: [" << VXoffset << ", " << VYoffset;
  cout << ", " << VZoffset << "].\n";

  //
  // Now get scaling parameters, for converting code units to
  // CGS units.
  //
  double CodePDratio = 0.0, ISMdensity = 0.0, Tism = 0.0,
         mass_per_particle = 0.0;

  seek = "BB_code_PDratio";
  str  = rp->find_parameter(seek);
  if (str == "")
    spdlog::error("{}: {}", "No pressure/density ratio", "BB_code_PDratio");
  else
    CodePDratio = atof(str.c_str());
  cout << "Pressure/density ratio = " << CodePDratio << "\n";

  seek = "BB_cgs_density";
  str  = rp->find_parameter(seek);
  if (str == "")
    spdlog::error("{}: {}", "No ism density", "BB_cgs_density");
  else
    ISMdensity = atof(str.c_str());

  seek = "BB_cgs_temperature";
  str  = rp->find_parameter(seek);
  if (str == "")
    spdlog::error("{}: {}", "No ism temperature", "BB_cgs_temperature");
  else
    Tism = atof(str.c_str());

  seek = "BB_cgs_mass_per_particle";
  str  = rp->find_parameter(seek);
  if (str == "")
    spdlog::error("{}: {}", "No mass_per_particle", "BB_cgs_mass_per_particle");
  else
    mass_per_particle = atof(str.c_str());

  //
  // Loop over variables and read in data
  //
  string infile, file_seek, off_seek, scale_seek, temp;
  double offset;
  double scaling;
  for (int v = 0; v < (SimPM.nvar - SimPM.ntracer); v++) {
    switch (v) {
      case RO:
        file_seek = "BB_filename_RO";
        // scale_seek = "BB_scaling_RO";
        off_seek = "BB_offset_RO";
        break;
      case PG:
        file_seek = "BB_filename_PG";
        // scale_seek = "BB_scaling_PG";
        off_seek = "BB_offset_PG";
        break;
      case VX:
        file_seek = "BB_filename_VX";
        // scale_seek = "BB_scaling_VX";
        off_seek = "BB_offset_VX";
        break;
      case VY:
        file_seek = "BB_filename_VY";
        // scale_seek = "BB_scaling_VY";
        off_seek = "BB_offset_VY";
        break;
      case VZ:
        file_seek = "BB_filename_VZ";
        // scale_seek = "BB_scaling_VZ";
        off_seek = "BB_offset_VZ";
        break;
      case BX:
        file_seek = "BB_filename_BX";
        // scale_seek = "BB_scaling_BX";
        off_seek = "BB_offset_BX";
        break;
      case BY:
        file_seek = "BB_filename_BY";
        // scale_seek = "BB_scaling_BY";
        off_seek = "BB_offset_BY";
        break;
      case BZ:
        file_seek = "BB_filename_BZ";
        // scale_seek = "BB_scaling_BZ";
        off_seek = "BB_offset_BZ";
        break;
      default:
        file_seek  = "";
        scale_seek = "";
        off_seek   = "";
        break;
    }

    infile = rp->find_parameter(file_seek);

    seek = rp->find_parameter(off_seek);
    if (seek == "")
      offset = 0.0;
    else
      offset = atof(seek.c_str());

    cout << "var=" << v << ", file=" << infile << " and offset=" << offset
         << "\n";

    switch (v) {
      case RO:
        //
        // code units has mean density =1, so to get to cgs units we
        // multiply by the required ISM mass density.
        //
        scaling = ISMdensity;
        break;

      case PG:
        //
        // code units has p=rho*c^2, so to get cgs units we multiply
        // by ISMdensity * v0^2 where v0 is the ratio of cgs sound
        // speed to code sound speed: v0^2 = kT/(mu*P_init).
        // It is obtained by multiplying by density, which is already
        // in cgs units, so we don't need ISMdensity in the scaling.
        //
        scaling = pconst.kB() * Tism / (mass_per_particle * CodePDratio);
        break;

      case VX:
      case VY:
      case VZ:
        //
        // Velocity scaling is the ratio of cgs-to-code sound speeds.
        //
        scaling = sqrt(pconst.kB() * Tism / (mass_per_particle * CodePDratio));
        break;

      case BX:
      case BY:
      case BZ:
        //
        // B-field scaling is with sqrt(Ro)*v0, because the Alfven
        // velocity in code units is defined to be unity.  My code uses
        // cgs units divided by sqrt(4pi), so take that into account...
        //
        scaling =
            sqrt(ISMdensity / (4.0 * M_PI))
            * sqrt(pconst.kB() * Tism / (mass_per_particle * CodePDratio));
        break;

      case SI:
        // The glm-mhd psi variable has no scaling, and should be set to
        // zero.
        scaling = 0.0;
        break;

      default:
        spdlog::error("{}: {}", "Bad variable", v);
        break;
    }

    //
    // Now either read data from source file, or else set it from
    // stuff we already know about the simulation.
    //
    if (infile == "" || infile == "NONE") {
      cout << "infile for var " << v << " doesn't exist.\n";
      if (v == PG) {
        //
        // Because simulations are isothermal, pressure is a fixed
        // multiple of density, which we get from the parameterfile.
        //
        cell *c = ggg->FirstPt();
        do {
          c->P[v] = (c->P[RO] * CodePDratio) * scaling;
        } while ((c = ggg->NextPt(c)) != 0);
      }
      else if (v == SI) {
        cell *c = ggg->FirstPt();
        do {
          c->P[v] = 0.0;
        } while ((c = ggg->NextPt(c)) != 0);
      }
      else {
        spdlog::error("{}: {}", "Unhandled variable with no source file", v);
      }
    }
    else {
      read_file(infile, v, scaling, offset, fmt::ptr(ggg));
    }
  }

  //
  // Now set values for tracers:
  //
  double data_vals[SimPM.ntracer];
  ostringstream ttt;
  for (int t = 0; t < SimPM.ntracer; t++) {
    ttt.str("");
    ttt << "BB_TR";
    ttt.width(2);
    ttt.fill('0');
    ttt << t;
    temp = ttt.str();
    cout << " looking for tracer " << temp << ";";
    seek = rp->find_parameter(temp);
    if (seek != "")
      data_vals[t] = atof(seek.c_str());
    else
      data_vals[t] = 0.0;
    cout << " value = " << data_vals[t] << "\n";
  }
  //
  // Assign tracer values everywhere on grid.
  //
  cell *c = ggg->FirstPt();
  do {
    for (int v = 0; v < SimPM.ntracer; v++)
      c->P[v + SimPM.ftr] = data_vals[v];
  } while ((c = ggg->NextPt(c)) != 0);

  return 0;
}

// ##################################################################
// ##################################################################

void IC_read_BBurkhart_data::read_file(
    const string infile,    ///< input fits file to read
    const int var,          ///< variable in state vector
    const double scaling,   ///< optional scaling of value.
    const double offset,    ///< optional offset of value, after scaling.
    class GridBaseClass *g  ///< pointer to grid
)
{
  //
  // We assume that NGrid[X/Y/Z] are set appropriately.
  // Set up the utility fitsio class to read data.
  //
  // class utility_fitsio fio ();
  class utility_fitsio *fio = 0;
  fio                       = new utility_fitsio();
  fitsfile *ff              = 0;
  int status = 0, err = 0;
  err = fits_open_file(&ff, infile.c_str(), READONLY, &status);
  if (status) {
    fits_report_error(stderr, status);
    spdlog::error("{}: {}", "open fits file", err);
  }

  int ndim = g->Ndim();
  int NG[ndim];
  size_t ncell = 1;
  for (int v = 0; v < ndim; v++) {
    NG[v] = g->iRange(static_cast<axes>(v)) / CI.get_integer_cell_size();
    ncell *= NG[v];
  }

  string HDU_NAME = "";  //"HDU1";

  // err = fio->check_fits_image_dimensions(ff,HDU_NAME,ndim,NG);
  // if (err) {
  //  //cout <<"image doesn't match grid size... parallelise me!\n";
  //  err=0;
  //}

  double *data = 0;
  // size_t Npix = SimPM.Ncell;
  // int npix[ndim];
  // for (int v=0; v<ndim;v++) npix[v] = SimPM.NG[v];

  double l_xmin[ndim], g_xmin[ndim];
  for (int v = 0; v < ndim; v++) {
    l_xmin[v] = g->Xmin(static_cast<axes>(v));
    g_xmin[v] = SimPM.Xmin[v];
  }
  double pix_size = g->DX();
  // data = mem.myalloc(data,Npix);

  err = fio->read_fits_image_to_data(
      ff, HDU_NAME, ndim, l_xmin, g_xmin, pix_size, NG, ncell, TDOUBLE, &data);
  if (err) spdlog::error("{}: {}", "Failed to read fits subset", err);

  cell *c      = g->FirstPt();
  size_t count = 0;
  do {
    c->P[var] = data[count] * scaling - offset;
    count++;
  } while ((c = g->NextPt(c)) != 0 && count < ncell);

  fits_close_file(ff, &status);
  data = mem.myfree(data);

  return;
}

// ##################################################################
// ##################################################################
