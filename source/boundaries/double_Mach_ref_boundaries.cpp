/// \file double_Mach_ref_boundaries.cpp
/// \brief Class definitions for double_Mach_ref boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include <spdlog/spdlog.h>

#include "boundaries/double_Mach_ref_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

// ##################################################################
// ##################################################################

int double_Mach_ref_bc::BC_assign_DMACH(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
#ifndef NDEBUG
  spdlog::info("Setting up DMACH boundary... starting.");
#endif  // NDEBUG

  if (b->data.empty()) {
    spdlog::error("{}: {}", "BC_assign_DMACH: empty boundary data", b->itype);
  }
  //
  // Set reference value to be downstream values.
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  b->refval[RO] = 1.4;
  b->refval[PG] = 1.0;
  b->refval[VX] = 0.0;
  b->refval[VY] = 0.0;
  b->refval[VZ] = 0.0;
  for (int v = par.ftr; v < par.nvar; v++)
    b->refval[v] = -1.0;

  //
  // Run through all boundary cells, and give them either upstream or
  // downstream value, depending on their position.
  //
  list<cell *>::iterator bpt = b->data.begin();
  unsigned int ct            = 0;
  double bpos                = 0.0;
  do {
    (*bpt)->isdomain = false;
    //
    // This is the boundary position:
    //
    bpos = 10.0 * par.simtime / sin(M_PI / 3.0) + 1.0 / 6.0
           + CI.get_dpos(*bpt, YY) / tan(M_PI / 3.0);

    if (CI.get_dpos(*bpt, XX) <= bpos) {
      (*bpt)->P[RO] = 8.0;
      (*bpt)->P[PG] = 116.5;
      (*bpt)->P[VX] = 7.14470958;
      (*bpt)->P[VY] = -4.125;
      (*bpt)->P[VZ] = 0.0;
      for (int v = par.ftr; v < par.nvar; v++)
        (*bpt)->P[v] = 1.0;
      (*bpt)->Ph[RO] = 8.0;
      (*bpt)->Ph[PG] = 116.5;
      (*bpt)->Ph[VX] = 7.14470958;
      (*bpt)->Ph[VY] = -4.125;
      (*bpt)->Ph[VZ] = 0.0;
      for (int v = par.ftr; v < par.nvar; v++)
        (*bpt)->Ph[v] = 1.0;
    }
    else {
      for (int v = 0; v < par.nvar; v++)
        (*bpt)->P[v] = b->refval[v];
      for (int v = 0; v < par.nvar; v++)
        (*bpt)->Ph[v] = b->refval[v];
    }
    ++bpt;
    ct++;
  } while (bpt != b->data.end());

  if (ct != b->data.size()) {
    spdlog::error(
        "{}: {}", "BC_assign_: missed some cells!", ct - b->data.size());
  }
#ifndef NDEBUG
  spdlog::info("Setting up DMACH boundary... finished.");
#endif  // NDEBUG
  return 0;
}

// ##################################################################
// ##################################################################

int double_Mach_ref_bc::BC_assign_DMACH2(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
#ifndef NDEBUG
  spdlog::info("Setting up DMACH2 boundary... starting.");
#endif  // NDEBUG

  if (b->dir != NO) {
    spdlog::error("{}: {}", "DMACH2 not internal boundary!", b->dir);
  }
  spdlog::info("DMACH2 boundary, from x=0 to x=1/6 at y=0, fixed bd.");
  if (b->refval) {
    spdlog::error(
        "{}: {}", "Already initialised memory in DMACH2 boundary refval",
        fmt::ptr(b->refval));
  }
  b->refval = mem.myalloc(b->refval, par.nvar);

  b->refval[RO] = 8.0;
  b->refval[PG] = 116.5;
  b->refval[VX] = 7.14470958;
  b->refval[VY] = -4.125;
  b->refval[VZ] = 0.0;
  for (int v = par.ftr; v < par.nvar; v++)
    b->refval[v] = 1.0;

  //
  // Now have to go from first point onto boundary and across to
  // x<=1/6
  //
  if (!b->data.empty()) {
    spdlog::error(
        "{}: {}", "BC_assign_DMACH2: Not empty boundary data", b->itype);
  }
  cell *c = grid->FirstPt();
  //
  // if running in parallel, need to check that YN boundary of grid
  // is also YN boundary of full domain:
  //
  if (c->pos[YY] > grid->idx()) {
    spdlog::info("domain is not at YN boundary, returning.");
    return 0;
  }
  cell *temp = 0;
  do {
    //
    // check if we are <1/6 in case we are in a parallel grid where
    // the domain is outside the DMR region.  This saves having to
    // rewrite the function in the parallel uniform grid.
    //
    if (CI.get_dpos(c, XX) <= 1. / 6.) {
      temp = c;
      while ((temp = grid->NextPt(temp, YN)) != 0) {
        for (int v = 0; v < par.nvar; v++)
          temp->P[v] = b->refval[v];
        for (int v = 0; v < par.nvar; v++)
          temp->Ph[v] = b->refval[v];
        b->data.push_back(temp);
        if (temp->isgd) {
          spdlog::error(
              "{}: {}", "BC_assign_DMACH2: Looking for Boundary cells!",
              fmt::ptr(temp));
        }
      }
    }
  } while ((c = grid->NextPt(c, XP)) && (CI.get_dpos(c, XX) <= 1. / 6.));

#ifndef NDEBUG
  spdlog::info("Setting up DMACH2 boundary... finished.");
#endif  // NDEBUG
  return 0;
}

// ##################################################################
// ##################################################################

int double_Mach_ref_bc::BC_update_DMACH(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const double simtime,       ///< current simulation time
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  double bpos              = 0.0;
  list<cell *>::iterator c = b->data.begin();
  for (c = b->data.begin(); c != b->data.end(); ++c) {
    //
    // This is the boundary position:
    //
    bpos = 10.0 * simtime / sin(M_PI / 3.0) + 1.0 / 6.0
           + CI.get_dpos(*c, YY) / tan(M_PI / 3.0);

    if (CI.get_dpos(*c, XX) <= bpos) {
      (*c)->Ph[RO] = 8.0;
      (*c)->Ph[PG] = 116.5;
      (*c)->Ph[VX] = 7.14470958;
      (*c)->Ph[VY] = -4.125;
      (*c)->Ph[VZ] = 0.0;
      for (int v = par.ftr; v < par.nvar; v++)
        (*c)->Ph[v] = 1.0;
    }
    else {
      for (int v = 0; v < par.nvar; v++)
        (*c)->Ph[v] = b->refval[v];
    }
    for (int v = 0; v < par.nvar; v++)
      (*c)->dU[v] = 0.0;
    if (cstep == maxstep) {
      for (int v = 0; v < par.nvar; v++)
        (*c)->P[v] = (*c)->Ph[v];
    }
  }
  return 0;
}

// ##################################################################
// ##################################################################

int double_Mach_ref_bc::BC_update_DMACH2(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int,
    const int)
{
  //
  // Fixed at all times, so no difference between full and half step.
  //
  list<cell *>::iterator c = b->data.begin();
  for (c = b->data.begin(); c != b->data.end(); ++c) {
    for (int v = 0; v < par.nvar; v++)
      (*c)->dU[v] = 0.;
    for (int v = 0; v < par.nvar; v++)
      (*c)->P[v] = b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*c)->Ph[v] = b->refval[v];
  }  // all cells.
  return 0;
}

// ##################################################################
// ##################################################################
