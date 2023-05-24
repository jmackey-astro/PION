/// \file overwrite_snapshot.cpp
/// \author Jonathan Mackey, Yvonne Fichtner
///
/// File for overwriting part of a snapshot with a supernova
///
/// - 2022.09.23 JM: Started on file.
/// - 2022.11.03 YF: Editing profile.
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

  // Edited by Yvonne Fichtner

  if (radius < 1e-5) {
    spdlog::info("Supernova: SN radius set to 0, no SN created, returning");
    return 0;
  }

  // Check if it is a normal SN
  // or if we have additionally a stellar wind (from a secondary)
  bool wind_source;
  double min_radius, max_radius;
  if (SWP.Nsources > 0) {
    wind_source = true;
    min_radius  = 2;  // SWP.params[isw]->radius !?! CHANGE+copy to correct file
    max_radius  = 2 + radius;  // SWP.params[isw]->radius + radius
    spdlog::info("Supernova: Creating SN with stellar wind source");
  }
  else {
    wind_source = false;
    min_radius  = 0;
    max_radius  = radius;
  }


  spdlog::info("Supernova: on finest level");
  // supernova prescription: 2 component from Dan Whalen et al. 2008
  //
  double test_ejecta_before;  //, test_ejecta, test_energy;
  double iter_ejecta, iter_energy, F_par, v_core;
  bool found_v_core;
  double dx    = ggg->DX();
  int n_exp    = 9;
  double v_max = 3e4 * 1e5;
  // time the ejected needed to reach the outer boundary
  double T_start     = max_radius * dx / v_max;
  int steps_per_cell = 1000000;                  // 10000;
  int N_steps        = steps_per_cell * radius;  // cells SN is build in
  test_ejecta_before = 0;
  iter_ejecta        = 0;
  iter_energy        = 0;
  found_v_core       = false;


  spdlog::info(
      "Supernova: dx {} v_max {}, T_start {}, StartTime {}, SimTime {}", dx,
      v_max, T_start, SimPM->starttime, SimPM->starttime);

  SimPM->starttime = SimPM->starttime + T_start;
  SimPM->simtime   = SimPM->simtime + T_start;

  spdlog::info(
      "Supernova: dx {} v_max {}, T_start {}, StartTime {}, SimTime {}", dx,
      v_max, T_start, SimPM->starttime, SimPM->starttime);


  // find tracer variables for H, He, Z and wind
  int tr_xh = 0, tr_xhe = 0, tr_xz = 0, tr_xc = 0, tr_xn = 0, tr_xo = 0,
      tr_hp = 0, tr_w = 0, tr_xd = 0;
  for (int v = 0; v < SimPM->ntracer; v++) {
    if (SimPM->tracers[v] == "X_H")
      tr_xh = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_He")
      tr_xhe = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_Z")
      tr_xz = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_C")
      tr_xc = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_N")
      tr_xn = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_O")
      tr_xo = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "H1+")
      tr_hp = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "WIND")
      tr_w = SimPM->ftr + v;
    else if (SimPM->tracers[v] == "X_D") {
      tr_xd             = SimPM->ftr + v;
      SimPM->tracers[v] = "X_ZS";
    }
  }
  if (tr_xh == 0 || tr_xhe == 0 || tr_xz == 0 || tr_xc == 0 || tr_xn == 0
      || tr_xo == 0 || tr_hp == 0 || tr_w == 0 || tr_xd == 0) {
    spdlog::error(
        "overwrite-snapshot: H {}, He {}, Z {}, C {}, N {}, O {}, H1p {}, wind {} or ZS {} tracers",
        tr_xh, tr_xhe, tr_xz, tr_xc, tr_xn, tr_xo, tr_hp, tr_w, tr_xd);
    exit(1);
  }
  spdlog::info(
      "Supernova: H {}, He {}, Z {}, C {}, N {}, O {}, H1p {}, wind {} or ZS {} tracers",
      tr_xh, tr_xhe, tr_xz, tr_xc, tr_xn, tr_xo, tr_hp, tr_w, tr_xd);
  spdlog::info("Supernova: found required tracers");

  // hardcode SN position at the origin
  std::array<double, MAX_DIM> snpos;
  for (int v = 0; v < MAX_DIM; v++)
    snpos[v] = 0.0;

  // Find ejecta before to compare
  // loop over cells.  All within r should have values overwritten
  cell *c = ggg->FirstPt();
  std::array<double, MAX_DIM> dpos;
  double distance = 0.0;
  do {
    CI.get_dpos(*c, dpos);
    distance = ggg->distance(snpos, dpos) - min_radius * dx;

    if ((distance > 0) && (distance < radius * dx)) {
      // distance: center of cell from real OR shifted 0 point
      test_ejecta_before =
          test_ejecta_before + (c->P[RO]) * ggg->CellVolume(*c, dx);
    }

    for (int v = 0; v < SimPM->nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = ggg->NextPt(*c)) != 0);

  if (test_ejecta_before > 0.05 * ejecta) {
    spdlog::error(
        "Supernova: test_ejecta_before {}, SN ejecta {}", test_ejecta_before,
        ejecta);
    exit(1);
  }


  for (int v = 1; v < N_steps; v++) {
    v_core =
        dx / T_start * (double)v
        / (double)steps_per_cell;  // from center to outer edge of SN region
    F_par = 1.0;

    iter_ejecta = 0.0;
    iter_energy = 0.0;
    // loop over cells.  All within r should have values overwritten
    c        = ggg->FirstPt();
    distance = 0.0;
    do {
      CI.get_dpos(*c, dpos);
      distance = ggg->distance(snpos, dpos) - min_radius * dx;


      if ((distance > 0) && (distance < radius * dx)) {
        c->P[VX] = (distance - 0.45 * dx) / T_start;
        if (c->P[VX] <= v_core) {
          c->P[RO] = F_par * pow(T_start, -3.0);
        }
        else {
          c->P[RO] =
              F_par * pow(T_start, -3.0) * pow(c->P[VX] / v_core, -n_exp);
        }
        mp->Set_Temp(c->P.data(), SimPM->EP.MinTemperature, SimPM->gamma);
        c->P[tr_xh]  = ejecta_X;
        c->P[tr_xhe] = ejecta_Y;
        // c->P[tr_xz]  = ejecta_Z;
        c->P[tr_xz] = 0.0;
        c->P[tr_xc] = 0.0;
        c->P[tr_xn] = 0.0;
        c->P[tr_xo] = 0.0;
        c->P[tr_hp] = 0.0;
        c->P[tr_w]  = 0.0;
        c->P[tr_xd] = ejecta_Z;
      }


      if ((distance > 0) && (distance < radius * dx)) {
        iter_ejecta += c->P[RO] * ggg->CellVolume(*c, dx);
        iter_energy +=
            0.5 * c->P[RO] * c->P[VX] * c->P[VX] * ggg->CellVolume(*c, dx);
      }

      for (int v = 0; v < SimPM->nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = ggg->NextPt(*c)) != 0);


    F_par = ejecta / iter_ejecta;

    iter_ejecta = 0.0;
    iter_energy = 0.0;
    // loop over cells.  All within r should have values overwritten
    c        = ggg->FirstPt();
    distance = 0.0;
    do {
      CI.get_dpos(*c, dpos);
      distance = ggg->distance(snpos, dpos) - min_radius * dx;


      if ((distance > 0) && (distance < radius * dx)) {
        c->P[VX] = (distance - 0.45 * dx) / T_start;
        if (c->P[VX] <= v_core) {
          c->P[RO] = F_par * pow(T_start, -3.0);
        }
        else {
          c->P[RO] =
              F_par * pow(T_start, -3.0) * pow(c->P[VX] / v_core, -n_exp);
        }
        mp->Set_Temp(c->P.data(), SimPM->EP.MinTemperature, SimPM->gamma);
        c->P[tr_xh]  = ejecta_X;
        c->P[tr_xhe] = ejecta_Y;
        // c->P[tr_xz]  = ejecta_Z;
        c->P[tr_xz] = 0.0;
        c->P[tr_xc] = 0.0;
        c->P[tr_xn] = 0.0;
        c->P[tr_xo] = 0.0;
        c->P[tr_hp] = 0.0;
        c->P[tr_w]  = 0.0;
        c->P[tr_xd] = ejecta_Z;
      }


      if ((distance > 0) && (distance < radius * dx)) {
        iter_ejecta += c->P[RO] * ggg->CellVolume(*c, dx);
        iter_energy +=
            0.5 * c->P[RO] * c->P[VX] * c->P[VX] * ggg->CellVolume(*c, dx);
      }

      for (int v = 0; v < SimPM->nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = ggg->NextPt(*c)) != 0);


    // spdlog::info(
    //   "Supernova: v {}, Nsteps {} vcore {} F_par {}",   v, N_steps, v_core,
    //   F_par);



    if (fabs(iter_energy - energy) < 0.0001 * energy) {
      found_v_core = true;
      break;
    }
  }

  spdlog::info("Supernova: found_v_core {}", found_v_core);
  if (found_v_core == false) {
    spdlog::error("No matching vcore and F could be found");
    exit(1);
  }

  spdlog::info("Supernova: vcore {} F_par {}", v_core, F_par);

  spdlog::info(
      "Supernova: test_ejecta_before {} iter_ejecta {} ejecta {} iter_energy {} energy {}",
      test_ejecta_before, iter_ejecta, ejecta, iter_energy, energy);


  return 0;
}



// ##################################################################
// ##################################################################
