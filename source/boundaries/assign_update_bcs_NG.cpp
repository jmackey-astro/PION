/// \file assign_update_bcs_NG.cpp
/// \brief Defines a class that inherits boundary types related to
///   static mesh-refinement (NG) and implements assignment and
///   update functions.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "assign_update_bcs_NG.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

using namespace std;



// ##################################################################
// ##################################################################



int assign_update_bcs_NG::assign_boundary_data(
    class SimParams &par,        ///< pointer to sim parameters
    const int level,             ///< level in grid hierarchy
    class GridBaseClass *grid,   ///< pointer to grid.
    class microphysics_base *mp  ///< pointer to microphysics
)
{
  // first call the Uniform Grid version.
  int err = assign_update_bcs::assign_boundary_data(par, level, grid, mp);
  if (err != 0) {
    spdlog::error(
        "{}: Expected {} but got {}", "assign_update_bcs::assign_boundary_data",
        err, 0);
    exit_pion(1);
  }

  //
  // Then check for NG-grid boundaries and assign data for them.
  //
  for (size_t i = 0; i < grid->BC_bd.size(); i++) {
#ifndef NDEBUG
    spdlog::debug("NG grid assign BCs: BC[{}] starting", i);
#endif
    switch (grid->BC_bd[i]->itype) {
      case FINE_TO_COARSE:
#ifndef NDEBUG
        spdlog::debug("NG grid setup: Assigning FINE_TO_COARSE BC");
#endif
        err += BC_assign_FINE_TO_COARSE(
            par, grid, grid->BC_bd[i], par.levels[level].child, 0);
        break;

      case COARSE_TO_FINE:
#ifndef NDEBUG
        spdlog::debug("assign_update_bcs_NG:: Assigning COARSE_TO_FINE BC");
#endif
        err += BC_assign_COARSE_TO_FINE(
            par, grid, grid->BC_bd[i], par.levels[level].parent);
        break;

      default:
#ifndef NDEBUG
        spdlog::debug("leaving BC {} alone in NG grid assign fn", i);
#endif
        break;
    }
  }
  return err;
}



// ##################################################################
// ##################################################################



int assign_update_bcs_NG::TimeUpdateInternalBCs(
    class SimParams &par,          ///< pointer to simulation parameters
    const int level,               ///< level in grid hierarchy
    class GridBaseClass *grid,     ///< pointer to grid.
    class FV_solver_base *solver,  ///< pointer to equations
    const double simtime,          ///< current simulation time
    const double dt,               ///< timestep
    const int cstep,
    const int maxstep)
{
#ifdef TEST_MPI_NG
  spdlog::info("assign_update_bcs_NG::TimeUpdateInternalBCs() running.");
#endif
  struct boundary_data *b;
  int err  = 0;
  size_t i = 0;
  // cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i = 0; i < grid->BC_bd.size(); i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
      case STWIND:
        err +=
            BC_update_STWIND(par, level, grid, simtime, dt, b, cstep, maxstep);
        break;
      case RADSHOCK:
        err += BC_update_RADSHOCK(par, grid, b, cstep, maxstep);
        break;

      case PERIODIC:
      case OUTFLOW:
      case ONEWAY_OUT:
      case INFLOW:
      case REFLECTING:
      case AXISYMMETRIC:
      case FIXED:
      case JETBC:
      case JETREFLECT:
      case DMACH:
      case DMACH2:
      case BCMPI:
      case COARSE_TO_FINE:
      case COARSE_TO_FINE_SEND:
      case COARSE_TO_FINE_RECV:
      case FINE_TO_COARSE_SEND:
      case FINE_TO_COARSE_RECV:
        //
        // boundaries not affected by NG grid are updated elsewhere
        //
        break;

      case FINE_TO_COARSE:
#ifdef TEST_MPI_NG
        spdlog::debug("found FINE_TO_COARSE boundary to update");
#endif
        err +=
            BC_update_FINE_TO_COARSE(par, solver, level, b, 0, cstep, maxstep);
        break;

      default:
#ifdef TEST_MPI_NG
        spdlog::error("no internal boundaries to update");
#endif
        spdlog::error(
            "{}: {}", "Unhandled BC: serial NG update internal", b->itype);
        break;
    }
  }
#ifdef TEST_MPI_NG
  spdlog::info("assign_update_bcs_NG::TimeUpdateInternalBCs() returns.");
#endif
  return err;
}

// ##################################################################
// ##################################################################

int assign_update_bcs_NG::TimeUpdateExternalBCs(
    class SimParams &par,          ///< pointer to sim parameters
    const int level,               ///< level in grid hierarchy
    class GridBaseClass *grid,     ///< pointer to grid.
    class FV_solver_base *solver,  ///< pointer to equations
    const double simtime,          ///< current simulation time
    const int cstep,
    const int maxstep)
{
#ifdef TEST_MPI_NG
  // cout <<"update_bcs_NG: external boundary update"<<endl;
#endif

  struct boundary_data *b;
  int err  = 0;
  size_t i = 0;
  // cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i = 0; i < grid->BC_bd.size(); i++) {
    b = grid->BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
      case PERIODIC:
        err += BC_update_PERIODIC(par, level, grid, b, cstep, maxstep);
        break;
      case OUTFLOW:
        err += BC_update_OUTFLOW(par, grid, b, cstep, maxstep);
        break;
      case ONEWAY_OUT:
        err += BC_update_ONEWAY_OUT(par, grid, b, cstep, maxstep);
        break;
      case INFLOW:
        err += BC_update_INFLOW(par, grid, b, cstep, maxstep);
        break;
      case REFLECTING:
        err += BC_update_REFLECTING(par, grid, b, cstep, maxstep);
        break;
      case AXISYMMETRIC:
        err += BC_update_AXISYMMETRIC(par, grid, b, cstep, maxstep);
        break;
      case FIXED:
        err += BC_update_FIXED(par, grid, b, cstep, maxstep);
        break;
      case JETBC:
        err += BC_update_JETBC(par, grid, b, cstep, maxstep);
        break;
      case JETREFLECT:
        err += BC_update_JETREFLECT(par, grid, b, cstep, maxstep);
        break;
      case DMACH:
        err += BC_update_DMACH(par, grid, simtime, b, cstep, maxstep);
        break;
      case DMACH2:
        err += BC_update_DMACH2(par, grid, b, cstep, maxstep);
        break;

      // skip these:
      case STWIND:
      case RADSHOCK:
      case BCMPI:
      case FINE_TO_COARSE:
        break;

        // get outer boundary of this grid from coarser grid.
      case COARSE_TO_FINE:
        // only update if at the start of a full step.
        if (cstep == maxstep
            // && (par.levels[level].step+2)%2==0
        ) {
          // cout <<"found COARSE_TO_FINE boundary to update\n";
          err += BC_update_COARSE_TO_FINE(
              par, solver, level, b, par.levels[level].step);
        }
        break;

      default:
        spdlog::error(
            "{}: {}", "Unhandled BC: serial NG update external", b->itype);
        break;
    }
  }
#ifdef TEST_NEST
  spdlog::info("updated NG-grid serial external BCs");
#endif
  return err;
}

// ##################################################################
// ##################################################################
