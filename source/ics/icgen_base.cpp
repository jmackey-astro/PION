/// \file icgen_base.cpp
/// \author Jonathan Mackey
/// \date 2018-05-04
///
/// Description:\n
/// This implements a set of routines that are common to serial,
/// parallel, and NG grid code for generating initial conditions.
///
/// Modifications:\n
/// 2018.05.04 JM: moved code from icgen.cpp/.h
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG

#include "icgen.h"
#include "icgen_base.h"
#include "microphysics/microphysics_base.h"

#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

void setup_ics_type(
    string ics,              ///< string giving type of ICs
    class ICsetup_base **ic  ///< pointer to address of IC class
)
{
  // invoke an appropriate class for whatever 'ics' is.
  if (ics == "ShockTube") *ic = new IC_shocktube();
  // some basic tests...
  else if (
      ics == "OrszagTang" || ics == "Uniform" || ics == "Advection"
      || ics == "AdvectSineWave" || ics == "KelvinHelmholtz"
      || ics == "KelvinHelmholtzStone" || ics == "FieldLoop"
      || ics == "FieldLoopVz" || ics == "FieldLoopStatic"
      || ics == "LiskaWendroffImplosion")
    *ic = new IC_basic_tests();

  else if (ics == "Jet" || ics == "JET" || ics == "jet") {
    *ic = new IC_jet();
  }

  else if (ics == "DoubleMachRef")
    *ic = new IC_basic_tests();

  else if (ics == "RadiativeShock" || ics == "RadiativeShockOutflow")
    *ic = new IC_radiative_shock();

  else if (ics == "LaserAblationAxi" || ics == "LaserAblation3D")
    *ic = new IC_laser_ablation();

  else if (ics == "ShockCloud")
    *ic = new IC_shock_cloud();

  else if (ics == "BlastWave" || ics == "BlastWave_File")
    *ic = new IC_blastwave();

  else if (
      ics == "PhotoEvaporatingClump" || ics == "PhotoEvaporatingClump2"
      || ics == "PhotoEvap_radial" || ics == "PhotoEvap_powerlaw"
      || ics == "PhotoEvap_paralleltest" || ics == "PhotoEvap_CloudClump")
    *ic = new IC_photoevaporatingclump();

  else if (
      ics == "PhotEvap_RandomClumps" || ics == "PERC" || ics == "PERC2"
      || ics == "PhotEvap_RandomClumps2")
    *ic = new IC_photevap_random_clumps();

  else if (
      ics == "PhotEvap_MultiClumps_FixNum" || ics == "PE_MC_FN"
      || ics == "PE_MC_FM" || ics == "PhotEvap_MultiClumps_FixMass")
    *ic = new IC_photevap_multi_clumps();

  else if (ics == "Clump_Spherical" || ics == "Clump_Axisymmetric")
    *ic = new IC_spherical_clump();

#ifdef CODE_EXT_SBII
  else if (
      ics == "StarBench_ContactDiscontinuity1"
      || ics == "StarBench_ContactDiscontinuity2"
      || ics == "StarBench_ContactDiscontinuity3"
      || ics == "StarBench_ContactDiscontinuity4"
      || ics == "StarBench_IFI_testA" || ics == "StarBench_IFI_testB"
      || ics == "StarBench_IFI_testC" || ics == "StarBench_IFI_V2"
      || ics == "StarBench_IrrCloud_Uniform"
      || ics == "StarBench_IrrCloud_IsoSph"
      || ics == "StarBench_TremblinCooling" || ics == "StarBench_Cone") {
    *ic = new IC_StarBench_Tests();
  }
#endif  // CODE_EXT_SBII

#ifdef HARPREETS_CODE_EXT
#ifndef EXCLUDE_HD_MODULE
  else if (ics == "HD_2D_ShockCloud")
    *ic = new IC_HD_2D_ShockCloud();
//  else if (ics=="HD_3D_ShockCloud")
//    *ic = new IC_HD_3D_ShockCloud();
#endif
#endif  // HARPREETS_CODE_EXT

#ifdef BBTURBULENCE_CODE_EXT
  else if (ics == "ReadBBTurbulence") {
    *ic = new IC_read_BBurkhart_data();
  }
#endif  // BBTURBULENCE_CODE_EXT

  else if (ics == "1Dto2D") {
    *ic = new IC_read_1Dto2D();
  }

  else
    spdlog::error("{}: {}", "BAD IC identifier", ics);
  if (!*ic) spdlog::error("{}: {}", "failed to init", ics);

  return;
}

// ##################################################################
// ##################################################################

int ICsetup_base::equilibrate_MP(
    class GridBaseClass *gg,
    class microphysics_base *mp,
    class ReadParams *rp,
    class SimParams &SimPM)
{

  if (!mp || !gg || !rp)
    spdlog::error(
        "{}: {}", "microphysics or grid not initialised.", fmt::ptr(mp));
  spdlog::debug("Init left  vec : {}", gg->FirstPt()->P);
  spdlog::debug("Init right vec : {}", gg->LastPt()->P);

  string seek = "InitIons";
  string s    = rp->find_parameter(seek);
  if (s == "" || s == "YES" || s == "Y" || s == "y") {
    // integrate ion fractions to equilibrium
    cell *c = gg->FirstPt();
    do {
      //
      // Check if cell is boundary data or not (can only be an internal
      // boundary, such as a stellar wind, since we are looping over cells
      // which are grid data).  If it is boundary data, then we don't want
      // to update anything, so we skip it
      //
      if (c->isbd) {
        spdlog::debug("skipping cell {} in equilibrate_MP() c->isbd", c->id);
      }
      else {
        mp->Init_ionfractions(c->P.data(), SimPM.gamma, -1);
      }
    } while ((c = gg->NextPt(*c)) != 0);
    spdlog::debug("init left  vec : {}", gg->FirstPt()->P);
    spdlog::debug("init right vec : {}", gg->LastPt()->P);

    // now do a long time integration to get to equilibrium.
    c           = gg->FirstPt();
    double tint = sqrt(SimPM.gamma * c->P[PG] / c->P[RO]);
    tint        = 50. * gg->DX()
           / tint;  // gives us 50 times the dynamical time for a cell.
    spdlog::debug("time to step by={}", tint);
    double tt = 0.;
    c         = gg->FirstPt();
    do {
      if (!c->isbd) {
        for (int i = 0; i < 50; i++) {
          mp->TimeUpdateMP(c->P.data(), c->P.data(), tint, SimPM.gamma, 0, &tt);
        }
      }
    } while ((c = gg->NextPt(*c)) != 0);
    spdlog::debug("Final left  vec : {}", gg->FirstPt()->P);
    spdlog::debug("Final right vec : {}", gg->LastPt()->P);
    c = gg->FirstPt();
    do {
      if (!c->isbd) {
        for (int i = 0; i < 50; i++) {
          mp->TimeUpdateMP(c->P.data(), c->P.data(), tint, SimPM.gamma, 0, &tt);
        }
      }
    } while ((c = gg->NextPt(*c)) != 0);
    spdlog::debug("Final left  vec : {}", gg->FirstPt()->P);
    spdlog::debug("Final right vec : {}", gg->LastPt()->P);
  }
  else if (s == "NO" || s == "N" || s == "n" || s == "no") {
    // initial values should be read from paramfile.
    string vb = "Tracer";
    cell *c   = 0;
    for (int i = 0; i < SimPM.ntracer; i++) {
      ostringstream t;
      t << vb << i;
      string var = t.str();
      spdlog::debug("var: {}", var);
      s = rp->find_parameter(var);
      if (s == "")
        spdlog::error("{}: {}", "Don't know what to set tracer to.", s);
      else {
        double tr = atof(s.c_str());
        c         = gg->FirstPt();
        do {
          c->P[SimPM.ftr + i] = tr;
        } while ((c = gg->NextPt(*c)) != 0);
      }
    }
    spdlog::debug("Final left  vec : {}", gg->FirstPt()->P);
    spdlog::debug("Final right vec : {}", gg->LastPt()->P);
  }
  else if (s == "LEAVE") {
    // do nothing! hopefully a subroutine has set them already.
  }
  else
    spdlog::error("{}: {}", "Bad InitIons specifier:", s);

  return 0;
}

// ##################################################################
// ##################################################################

int ICsetup_base::AddNoise2Data(
    class GridBaseClass *grid, class SimParams &SimPM, int n, double frac)
{
  int seed = 975;
#ifdef PARALLEL
  seed += sub_domain->get_myrank();
  bool true_edge = false;
#endif
  srand(seed);
  class cell *cpt;
  double avg  = 0.;
  long int ct = 0;
  switch (n) {
    case 1:
      spdlog::debug(
          "\tAdding random noise to pressure at fractional level of {}", frac);
      cpt = grid->FirstPt();
      do {
        if (!cpt->isbd) {
          avg += cpt->P[PG];
          ct++;
        }
      } while ((cpt = grid->NextPt(*cpt)) != 0);
#ifdef PARALLEL
      avg = sub_domain->global_operation_double(SUM, avg);
#endif
      spdlog::debug("avg = {}\t ct = {}", avg, ct);
      avg /= static_cast<double>(SimPM.Ncell);
      avg *= frac;  // avg is now a fraction 'frac' of the mean pressure
                    // on the grid.
      spdlog::debug("avg = {}", avg);
      cpt = grid->FirstPt();
      do {
#ifdef SERIAL
        if (!cpt->isedge && !cpt->isbd)
#endif
#ifdef PARALLEL
          //
          // We want to exclude edge cells, but only those at the edge
          // of the full domain, not the local domain.
          //
          true_edge = false;
        if (cpt->isedge) {
          //
          // find out which direction the edge is (may be more than
          // one!). If any edge has no neighbour process, then it must
          // be a full-domain boundary.  If the neigbour doesn't
          // exist, ngbproc(dir)<0.
          //
          // x-dir
          if (cpt->isedge % 3 == 1) {  // XN boundary
            if (sub_domain->get_neighbour_rank(XN) < 0) true_edge = true;
          }
          else if (cpt->isedge % 3 == 2) {  // XP boundary
            if (sub_domain->get_neighbour_rank(XP) < 0) true_edge = true;
          }
          // y-dir
          if (SimPM.ndim > 1) {
            if ((cpt->isedge % 9) / 3 == 1) {  // YN boundary
              if (sub_domain->get_neighbour_rank(YN) < 0) true_edge = true;
            }
            else if ((cpt->isedge % 9) / 3 == 2) {  // YP boundary
              if (sub_domain->get_neighbour_rank(YP) < 0) true_edge = true;
            }
          }
          // z-dir
          if (SimPM.ndim > 2) {
            if (cpt->isedge / 9 == 1) {  // ZN boundary
              if (sub_domain->get_neighbour_rank(ZN) < 0) true_edge = true;
            }
            else if (cpt->isedge / 9 == 2) {  // ZP boundary
              if (sub_domain->get_neighbour_rank(ZP) < 0) true_edge = true;
            }
          }
        }
        if (true_edge == false && !cpt->isbd)
#endif
        {  // Don't want to alter any edge cells.
          cpt->P[PG] += avg * (static_cast<double>(rand()) / RAND_MAX - 0.5);
        }
      } while ((cpt = grid->NextPt(*cpt)) != 0);
      break;

    case 2:
      spdlog::debug(
          "Adding adiabatic random noise to cells at fractional level of {}",
          frac);
      cpt = grid->FirstPt();
      do {
#ifdef SERIAL
        if (!cpt->isedge && !cpt->isbd)
#endif
#ifdef PARALLEL
          //
          // We want to exclude edge cells, but only those at the edge
          // of the full domain, not the local domain.
          //
          true_edge = false;
        if (cpt->isedge) {
          //
          // find out which direction the edge is (may be more than
          // one!). If any edge has no neighbour process, then it must
          // be a full-domain boundary.  If the neigbour doesn't
          // exist, ngbproc(dir)<0.
          //
          // x-dir
          if (cpt->isedge % 3 == 1) {  // XN boundary
            if (sub_domain->get_neighbour_rank(XN) < 0) true_edge = true;
          }
          else if (cpt->isedge % 3 == 2) {  // XP boundary
            if (sub_domain->get_neighbour_rank(XP) < 0) true_edge = true;
          }
          // y-dir
          if (SimPM.ndim > 1) {
            if ((cpt->isedge % 9) / 3 == 1) {  // YN boundary
              if (sub_domain->get_neighbour_rank(YN) < 0) true_edge = true;
            }
            else if ((cpt->isedge % 9) / 3 == 2) {  // YP boundary
              if (sub_domain->get_neighbour_rank(YP) < 0) true_edge = true;
            }
          }
          // z-dir
          if (SimPM.ndim > 2) {
            if (cpt->isedge / 9 == 1) {  // ZN boundary
              if (sub_domain->get_neighbour_rank(ZN) < 0) true_edge = true;
            }
            else if (cpt->isedge / 9 == 2) {  // ZP boundary
              if (sub_domain->get_neighbour_rank(ZP) < 0) true_edge = true;
            }
          }
        }
        if (true_edge == false && !cpt->isbd)
#endif
        {  // Don't want to alter any edge cells.

          double temp;
          // if(cpt->isedge==0) {    // Don't want to alter any edge
          //
          // final pressure = initial pressure * (1 +/- frac)
          // rho = pressure^(1/gamma)
          //
          temp = 2. * frac * (static_cast<double>(rand()) / RAND_MAX - 0.5);
          cpt->P[PG] *= 1 + temp;
          cpt->P[RO] *= exp(log(1 + temp) / SimPM.gamma);
          //}
        }
      } while ((cpt = grid->NextPt(*cpt)) != 0);
      break;

    case 3:
      spdlog::debug(
          "Adding adiabatic wave to cells at fractional level of {}", frac);
#ifdef PARALLEL
      spdlog::error(
          "{}: {}", "adiabatic wave noise not implemented in parallel", 3);
#endif
      // First get average pressure value.
      cpt = grid->FirstPt();
      do {
        if (!cpt->isedge && !cpt->isbd) avg += cpt->P[PG];
        ct++;
      } while ((cpt = grid->NextPt(*cpt)) != 0);
      spdlog::debug("avg = {}\t ct = {}", avg, ct);
      avg /= static_cast<double>(ct);
      // Now add wave to preshock state.
      cpt = grid->FirstPt();
      do {
        double temp;
        if (cpt->isedge == 0 && !cpt->isbd
            && cpt->P[PG] < avg) {  // Don't want to alter any edge cells.
          temp = frac
                 * sin(2. * M_PI * (CI.get_dpos(*cpt, YY) / SimPM.Range[YY])
                       * (SimPM.NG[YY] / 50));
          cpt->P[PG] *= 1 + temp;
          cpt->P[RO] *= exp(log(1 + temp) / SimPM.gamma);
        }
      } while ((cpt = grid->NextPt(*cpt)) != 0);
      break;

      //
      // Isothermal perturbation.
      //
    case 4:
      spdlog::debug(
          "Adding isothermal random noise to cells at fractional level of {}",
          frac);
      cpt = grid->FirstPt();
      do {
#ifdef SERIAL
        if (!cpt->isedge && !cpt->isbd)
#endif
#ifdef PARALLEL
          //
          // We want to exclude edge cells, but only those at the edge
          // of the full domain, not the local domain.
          //
          true_edge = false;
        if (cpt->isedge) {
          //
          // find out which direction the edge is (may be more than
          // one!). If any edge has no neighbour process, then it must
          // be a full-domain boundary.  If the neigbour doesn't
          // exist, ngbproc(dir)<0.
          //
          // x-dir
          if (cpt->isedge % 3 == 1) {  // XN boundary
            if (sub_domain->get_neighbour_rank(XN) < 0) true_edge = true;
          }
          else if (cpt->isedge % 3 == 2) {  // XP boundary
            if (sub_domain->get_neighbour_rank(XP) < 0) true_edge = true;
          }
          // y-dir
          if (SimPM.ndim > 1) {
            if ((cpt->isedge % 9) / 3 == 1) {  // YN boundary
              if (sub_domain->get_neighbour_rank(YN) < 0) true_edge = true;
            }
            else if ((cpt->isedge % 9) / 3 == 2) {  // YP boundary
              if (sub_domain->get_neighbour_rank(YP) < 0) true_edge = true;
            }
          }
          // z-dir
          if (SimPM.ndim > 2) {
            if (cpt->isedge / 9 == 1) {  // ZN boundary
              if (sub_domain->get_neighbour_rank(ZN) < 0) true_edge = true;
            }
            else if (cpt->isedge / 9 == 2) {  // ZP boundary
              if (sub_domain->get_neighbour_rank(ZP) < 0) true_edge = true;
            }
          }
        }
        if (true_edge == false && !cpt->isbd)
#endif
        {
          // Don't want to alter any edge cells.
          double temp;
          //
          // final pressure = initial pressure * (1 +/- frac)
          // final density  = initial density  * (1 +/- frac)
          //
          temp = 2.0 * frac * (static_cast<double>(rand()) / RAND_MAX - 0.5);
          cpt->P[PG] *= 1.0 + temp;
          cpt->P[RO] *= 1.0 + temp;
        }
      } while ((cpt = grid->NextPt(*cpt)) != 0);
      break;

    default:
      spdlog::error(
          "\t Error, don't know what type of noise corresponds to {}...", n);
      return (1);
  }
  return (0);
}  // AddNoise2Data()

// ##################################################################
// ##################################################################

int ICsetup_base::SmoothData(int smooth)
{
  return (0);
}
