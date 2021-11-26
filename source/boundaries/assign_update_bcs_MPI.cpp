/// \file assign_update_bcs_MPI.cpp
/// \brief Defines a class that inherits all boundary types and
///   assigns and updated the boundaries.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "assign_update_bcs_MPI.h"


#include <spdlog/spdlog.h>

#include <sstream>
using namespace std;



// ##################################################################
// ##################################################################



int assign_update_bcs_MPI::assign_boundary_data(
    class SimParams &par,        ///< pointer to simulation parameters
    const int level,             ///< level in grid hierarchy
    class GridBaseClass *grid,   ///< pointer to grid.
    class microphysics_base *mp  ///< pointer to microphysics
)
{
#ifdef TEST_MPI_BC
  class Sub_domain *ppar = &(par.levels[level].sub_domain);
  spdlog::debug("{} Setting up MPI boundaries...", ppar->get_myrank());
#endif
  int err = assign_update_bcs::assign_boundary_data(par, level, grid, mp);

  if (err != 0)
    spdlog::error(
        "{}: Expected {} but got {}", "assign_update_bcs::assign_boundary_data",
        err, 0);
  //
  // Loop through all boundaries, and assign data to MPI one.
  //
  for (size_t i = 0; i < grid->BC_bd.size(); i++) {
    switch (grid->BC_bd[i]->itype) {
      case BCMPI:
#ifdef TEST_MPI_BC
        spdlog::debug(
            "{} assigning MPI boundary in dir {}", ppar->get_myrank(), i);
#endif
        err += BC_assign_BCMPI(par, level, grid, grid->BC_bd[i], BC_MPItag);
        break;
      case PERIODIC:
        break;
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
      case STWIND:
      case FINE_TO_COARSE:
      case COARSE_TO_FINE:
      case FINE_TO_COARSE_SEND:
      case FINE_TO_COARSE_RECV:
      case COARSE_TO_FINE_SEND:
      case COARSE_TO_FINE_RECV:
        break;  // assigned in NG grid class
      default:
        spdlog::error(
            "{}: {}", "assign_update_bcs_MPI::Unhandled BC: assign",
            grid->BC_bd[i]->itype);
        break;
    }
  }
#ifdef TEST_MPI_BC
  spdlog::debug("{}: finished assigning boundaries\n", ppar->get_myrank());
#endif
  return (err);
}



// ##################################################################
// ##################################################################



int assign_update_bcs_MPI::TimeUpdateExternalBCs(
    class SimParams &par,          ///< pointer to sim parameters
    const int level,               ///< level in grid hierarchy
    class GridBaseClass *grid,     ///< pointer to grid.
    class FV_solver_base *solver,  ///< pointer to equations
    const double simtime,          ///< current simulation time
    const int cstep,
    const int maxstep)
{
#ifdef TEST_MPI_NG
  spdlog::info("update_bcs_MPI: external boundary update");
#endif
  struct boundary_data *b;
  int err  = 0;
  size_t i = 0;

  // make a map so that boundaries are updated pairwise.
  // ix[] is the position of this process in each direction of the
  // block of MPI processes.
  int nb = grid->BC_bd.size();
  int map[nb];
  int rank = par.levels[level].sub_domain.get_myrank();
  int ix[MAX_DIM];
  par.levels[level].sub_domain.get_domain_coordinates(rank, ix);
  // x-direction
  for (int j = 0; j < par.ndim; j++) {
    if (ix[j] % 2 == 0) {
      map[2 * j]     = 2 * j;
      map[2 * j + 1] = 2 * j + 1;
    }
    else {
      map[2 * j]     = 2 * j + 1;
      map[2 * j + 1] = 2 * j;
    }
  }
  for (int j = 2 * par.ndim; j < nb; j++) {
    map[j] = j;
  }

  for (i = 0; i < grid->BC_bd.size(); i++) {
    b = grid->BC_bd[map[i]];
#ifdef TEST_MPI_BC
    spdlog::debug("updating bc {} with type {}", map[i], b->type);
#endif
    switch (b->itype) {

      case PERIODIC:
        err += BC_update_PERIODIC(par, level, grid, b, cstep, maxstep);
        break;
      case BCMPI:
        err += BC_update_BCMPI(par, level, grid, b, cstep, maxstep, BC_MPItag);
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

      case RADSHOCK:
      case RADSH2:
      case STWIND:
        //
        // internal bcs updated separately
        //
        break;

      default:
        spdlog::error(
            "{}: {}", "assign_update_bcs_MPI::Unhandled BC: external",
            b->itype);
        break;
    }

    par.levels[0].sub_domain.barrier("External BC update");

    // if (i==XP || (i==YP && par.ndim>=2) || (i==ZP && par.ndim==3)) {
    //  // Need to make sure all processors are updating boundaries along
    //  // each axis together, so that there is no deadlock from one
    //  // process sending data in y/z-dir into a processor that is still
    //  // working on x-dir.
    //  //      cout <<"Barrier synching for dir="<<i<<" with
    //  1=XP,2=YP,3=ZP.\n"; ostringstream tmp; tmp
    //  <<"assign_update_bcs_MPI_TimeUpdateExternalBCs_"<<i;
    //  par.levels[0].sub_domain.barrier(tmp.str());
    //  tmp.str("");
    //}
  }

  return 0;
}

// ##################################################################
// ##################################################################
