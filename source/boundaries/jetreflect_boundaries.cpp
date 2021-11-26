/// \file jetreflect_boundaries.cpp
/// \brief Class definitions for jetreflect boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/jetreflect_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

// ##################################################################
// ##################################################################

int jetreflect_bc::BC_assign_JETREFLECT(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = b->ondir;
  if (b->data.empty()) {
    spdlog::error("{}: {}", "BC_assign_: empty boundary data", b->itype);
  }

  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = 1.0;
  //
  // Set Normal velocity multiplier to -1 for reflection.
  //
  switch (offdir) {
    case XN:
    case XP:
      b->refval[VX] = -1.0;
      break;
    case YN:
    case YP:
      b->refval[VY] = -1.0;
      break;
    case ZN:
    case ZP:
      b->refval[VZ] = -1.0;
      break;
    default:
      spdlog::error("{}: {}", "BAD DIRECTION", offdir);
  }  // Set Normal velocity direction.
  //
  // Set tangential B-field multiplier to -1 for this boundary, to
  // allow a dipole to exist
  //
  if (par.eqntype == EQMHD || par.eqntype == EQGLM || par.eqntype == EQFCD) {
    switch (offdir) {
      case XN:
      case XP:
        b->refval[BY] = b->refval[BZ] = -1.0;
        break;
      case YN:
      case YP:
        b->refval[BX] = b->refval[BZ] = -1.0;
        break;
      case ZN:
      case ZP:
        b->refval[BX] = b->refval[BY] = -1.0;
        break;
      default:
        spdlog::error("{}: {}", "BAD DIRECTION", offdir);
    }  // Set normal b-field direction.
  }    // if B-field exists

  //
  // Now go through each column of boundary points and assign values
  // to them.
  //
  list<cell *>::iterator bpt = b->data.begin();
  cell *temp                 = 0;
  unsigned int ct            = 0;

  do {
    temp = (*bpt);
    for (int v = 0; v > (*bpt)->isedge; v--) {
      temp = grid->NextPt(temp, ondir);
    }
    if (!temp) {
      spdlog::error(
          "{}: {}", "Got lost assigning jet-reflecting bcs.", temp->id);
    }
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->P[v] = temp->P[v] * b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->Ph[v] = temp->P[v] * b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->dU[v] = 0.0;
    (*bpt)->npt = temp;
    ++bpt;
    ct++;
  } while (bpt != b->data.end());

  if (ct != b->data.size()) {
    spdlog::error(
        "{}: {}", "BC_assign_: missed some cells!", ct - b->data.size());
  }
  return 0;
}

// ##################################################################
// ##################################################################

int jetreflect_bc::BC_update_JETREFLECT(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  //
  // same routine as for reflecting, except the normal B is
  // unchanged, but the tangential is reversed.
  //
  list<cell *>::iterator c = b->data.begin();
  for (c = b->data.begin(); c != b->data.end(); ++c) {
    for (int v = 0; v < par.nvar; v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v] * b->refval[v];
    }
    for (int v = 0; v < par.nvar; v++)
      (*c)->dU[v] = 0.;
    if (cstep == maxstep) {
      for (int v = 0; v < par.nvar; v++) {
        (*c)->P[v] = (*c)->npt->P[v] * b->refval[v];
      }
    }
  }  // all cells.
  return 0;
}

// ##################################################################
// ##################################################################
