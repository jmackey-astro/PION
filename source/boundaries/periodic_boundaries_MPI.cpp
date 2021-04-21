/// \file periodic_boundaries_MPI.cpp
/// \brief Class definitions for periodic boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/periodic_boundaries_MPI.h"
#include "decomposition/MCMD_control.h"
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
  // which is already pointed to by ppar->ngbprocs[b->dir]
  // So I just have to call BC_assign_BCMPI and it will do the job.
  int err                 = 0;
  class MCMDcontrol *ppar = &(par.levels[level].MCMD);
  if (ppar->ngbprocs[b->dir] < 0) {
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
  // which is already pointed to by ppar->ngbprocs[b->dir]
  // So I just have to call BC_update_BCMPI and it will do the job.
  //
  int err                 = 0;
  class MCMDcontrol *ppar = &(par.levels[level].MCMD);
  if (ppar->ngbprocs[b->dir] < 0) {
#ifdef TESTING
    cout << "BC_update_PERIODIC: non-communicating periodic BC in ";
    cout << "direction " << b->dir << "\n";
#endif
    err = periodic_bc::BC_update_PERIODIC(par, level, grid, b, cstep, maxstep);
  }
  else {
#ifdef TESTING
    cout << "BC_update_PERIODIC: communicating periodic BC in ";
    cout << "direction " << b->dir << "\n";
    cout << "BC_update_PERIODIC: calling mpi update BC function\n";
#endif
    err = BC_update_BCMPI(par, level, grid, b, cstep, maxstep, BC_PERtag);
  }
  return err;
}

// ##################################################################
// ##################################################################
