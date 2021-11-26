/// \file fixed_boundaries.cpp
/// \brief Class definitions for fixed boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include <spdlog/spdlog.h>

#include "boundaries/fixed_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

// ##################################################################
// ##################################################################

int fixed_bc::BC_assign_FIXED(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
#ifndef NDEBUG
  spdlog::info(" setup_fixed_grid::BC_assign_FIXED starting");
#endif
  enum direction ondir = b->ondir;
  if (b->data.empty()) {
    spdlog::error("{}: {}", "BC_assign_FIXED: empty boundary data", b->itype);
  }
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  list<cell *>::iterator bpt = b->data.begin();
  cell *temp                 = 0;
  unsigned int ct            = 0;
  //
  // First find an on-grid cell near a boundary point.  Because of
  // corner cells, we can't guarantee that every boundary cell will
  // reach an on-grid cell by moving in the on-grid direction.
  //
#ifndef NDEBUG
  spdlog::debug("Finding first on-grid cell, size={}", b->data.size());
#endif
  do {
    (*bpt)->isdomain = false;
    temp             = (*bpt);
    CI.print_cell(temp);
    for (int v = 0; v > (*bpt)->isedge; v--) {
      temp = grid->NextPt(temp, ondir);
    }
  } while (!temp->isgd);
  if (!temp) spdlog::error("{}: {}", "Got lost assigning FIXED bcs.", temp->id);

    //
    // Now set reference value to be the on-grid value.
    //
#ifndef NDEBUG
  spdlog::debug("Setting reference value");
#endif
  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = temp->P[v];
  //
  // Initialise all the values to be the fixed value.
  //
  bpt = b->data.begin();
  do {
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->P[v] = b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->Ph[v] = b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->dU[v] = 0.;
    ++bpt;
    ct++;
  } while (bpt != b->data.end());

  if (ct != b->data.size()) {
    spdlog::error(
        "{}: {}", "BC_assign_FIXED: missed some cells!", ct - b->data.size());
  }

  return 0;
}

// ##################################################################
// ##################################################################

int fixed_bc::BC_update_FIXED(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  //
  // Fixed means all boundary points have same fixed value, stored
  // in refval.
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
