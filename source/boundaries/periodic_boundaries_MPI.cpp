/// \file periodic_boundaries_MPI.cpp
/// \brief Class definitions for periodic boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "boundaries/periodic_boundaries_MPI.h"
#include "sub_domain/sub_domain.h"
using namespace std;

// ##################################################################
// ##################################################################

int periodic_pllel_bc::BC_assign_PERIODIC(
    class SimParams &par,       ///< pointer to simulation parameters
    const int level,            ///< level in grid hierarchy
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
  //
  // For parallel grid, periodic data may be on a different proc.,
  // which is already pointed to by ppar->get_neighbour_ranks(b->dir)
  // So I just have to call BC_assign_BCMPI and it will do the job.
  int err                = 0;
  class Sub_domain *ppar = &(par.levels[level].sub_domain);
  if (ppar->get_neighbour_rank(b->dir) < 0
      || ppar->get_neighbour_rank(b->dir) == ppar->get_myrank()) {
    // cout <<"BC_assign_PERIODIC: non comm periodic in direction ";
    // cout <<b->dir<<"\n";
    err = periodic_bc::BC_assign_PERIODIC(par, level, grid, b);
  }
  else {
    // cout<<"BC_assign_PERIODIC: communicating periodic bc in ";
    // cout<<"direction "<<b->dir<<"\n";
    // cout<<"BC_assign_PERIODIC: calling mpi assign BC function\n";
    err = BC_assign_BCMPI(par, level, grid, b, BC_PERtag);
  }
  return err;
}

// ##################################################################
// ##################################################################

int periodic_pllel_bc::BC_update_PERIODIC(
    class SimParams &par,       ///< pointer to simulation parameters
    const int level,            ///< level in grid hierarchy
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  //
  // For parallel grid, periodic data can be on a different proc.,
  // which is already pointed to by ppar->get_neighbour_rank(b->dir)
  // So I just have to call BC_update_BCMPI and it will do the job.
  //
  int err                = 0;
  class Sub_domain *ppar = &(par.levels[level].sub_domain);
  if (ppar->get_neighbour_rank(b->dir) < 0
      || ppar->get_neighbour_rank(b->dir) == ppar->get_myrank()) {
#ifndef NDEBUG
    spdlog::debug(
        "BC_update_PERIODIC: non-communicating periodic BC in direction {}",
        static_cast<int>(b->dir));
#endif
    err = periodic_bc::BC_update_PERIODIC(par, level, grid, b, cstep, maxstep);
  }
  else {
#ifndef NDEBUG
    spdlog::debug(
        "BC_update_PERIODIC: communicating periodic BC in direction {}",
        static_cast<int>(b->dir));
    spdlog::debug("BC_update_PERIODIC: calling mpi update BC function");
#endif
    err = BC_update_BCMPI(par, level, grid, b, cstep, maxstep, BC_PERtag);
  }
  return err;
}

// ##################################################################
// ##################################################################
