/// \file overwrite_snapshot.cpp
/// \author Jonathan Mackey
///
/// File for overwriting part of a snapshot with a supernova
///
/// - 2022.09.23 JM: Started on file.
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

IC_overwrite_snapshot::IC_overwrite_snapshot() {}

// ##################################################################
// ##################################################################

IC_overwrite_snapshot::~IC_overwrite_snapshot() {}

// ##################################################################
// ##################################################################

int IC_overwrite_snapshot::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  spdlog::info("Setting up supernova ICs on finest grid");
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) {
    spdlog::error("{}: {}", "null pointer to grid!", fmt::ptr(ggg));
    exit(1);
  }
  ICsetup_base::rp = rrp;
  if (!rp) {
    spdlog::error("{}: {}", "null pointer to ReadParams", fmt::ptr(rp));
    exit(1);
  }
  string seek, str;

  // Get SN radius, energy, mass, X_H, X_He from parameterfile
  seek = "SN-radius-Ncells";
  str  = rp->find_parameter(seek);
  if (str == "") {
    spdlog::error("{}: {}", "didn't find parameter", seek);
    exit(1);
  }
  // spdlog::info("Radius: {} {}",str, atoi(str.c_str()));
  int radius = atoi(str.c_str());

  seek = "SN-energy";
  str  = rp->find_parameter(seek);
  if (str == "") {
    spdlog::error("{}: {}", "didn't find parameter", seek);
    exit(1);
  }
  double energy = atof(str.c_str());

  seek = "SN-ejecta-mass";
  str  = rp->find_parameter(seek);
  if (str == "") {
    spdlog::error("{}: {}", "didn't find parameter", seek);
    exit(1);
  }
  double ejecta = atof(str.c_str());

  seek = "SN-ejecta-H-frac";
  str  = rp->find_parameter(seek);
  if (str == "") {
    spdlog::error("{}: {}", "didn't find parameter", seek);
    exit(1);
  }
  double ejecta_X = atof(str.c_str());

  seek = "SN-ejecta-He-frac";
  str  = rp->find_parameter(seek);
  if (str == "") {
    spdlog::error("{}: {}", "didn't find parameter", seek);
    exit(1);
  }
  double ejecta_Y = atof(str.c_str());

  double ejecta_Z = 1.0 - ejecta_X - ejecta_Y;

  spdlog::info("Supernova: got required parameters");

  // overwrite data if within radius but only on finest level.
  if (ggg->level() != SimPM->grid_nlevels - 1) {
    spdlog::info("Supernova: not on finest level, returning");
    return 0;
  }


  spdlog::info("Supernova: on finest level");
  // simple supernova prescription: uniform density and energy within
  // a spherical radius.  Pressure = internal energy * (gamma-1)
  double mean_dens, mean_pres;
  double dx  = ggg->DX();
  double vol = 4.0 * pconst.pi() * pow(radius * dx, 3.0) / 3.0;
  mean_dens  = ejecta / vol;
  mean_pres  = energy * (SimPM->gamma - 1.0) / vol;
  spdlog::info(
      "Supernova: dx {} vol {} dens {} pres {}", dx, vol, mean_dens, mean_pres);

  // find tracer variables for H, He
  int tr_xh = 0, tr_xhe = 0, tr_xz = 0;
  for (int v = 0; v < SimPM->ntracer; v++) {
    if (SimPM->tracers[v] == "X_H")
      tr_xh = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_He")
      tr_xhe = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_Z")
      tr_xz = SimPM->ftr + v;
  }
  if (tr_xh == 0 || tr_xhe == 0 || tr_xz == 0) {
    spdlog::error(
        "overwrite-snapshot: H {}, He {} or Z {} tracers", tr_xh, tr_xhe,
        tr_xz);
    exit(1);
  }
  spdlog::info("Supernova: found required tracers");

  // hardcode SN position at the origin
  std::array<double, MAX_DIM> snpos;
  for (int v = 0; v < MAX_DIM; v++)
    snpos[v] = 0.0;

  // loop over cells.  All within r should have values overwritten
  cell *c = ggg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  double distance = 0.0;
  do {
    CI.get_dpos(*c, dpos);
    distance = ggg->distance(snpos, dpos);
    // spdlog::info("Supernova: cell id {}, distance =
    // {:12.6e}",c->id,distance);

    if (distance < radius * dx) {
      // spdlog::info("\twithin radius {:12.6e}, settting vals {} {} {} {} {}",
      //  radius*dx, mean_dens, mean_pres, ejecta_X, ejecta_Y, ejecta_Z);
      c->P[RO]     = mean_dens;
      c->P[PG]     = mean_pres;
      c->P[tr_xh]  = ejecta_X;
      c->P[tr_xhe] = ejecta_Y;
      c->P[tr_xz]  = ejecta_Z;
    }
    for (int v = 0; v < SimPM->nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = ggg->NextPt(*c)) != 0);

  return 0;
}



// ##################################################################
// ##################################################################
