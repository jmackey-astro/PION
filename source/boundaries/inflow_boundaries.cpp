/// \file inflow_boundaries.cpp
/// \brief Class definitions for inflow boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/inflow_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

// ##################################################################
// ##################################################################

int inflow_bc::BC_assign_INFLOW(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
  enum direction ondir = b->ondir;

  if (b->data.empty()) {
    rep.error("BC_assign_INFLOW: empty boundary data", b->itype);
  }
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  list<cell *>::iterator bpt = b->data.begin();
  cell *temp;
  unsigned int ct = 0;
  do {
    //
    // Find the cell to point to.  This should be the first cell on
    // the other side of the boundary, located at x[baxis]=bloc.
    // Can't just use the ->isgd property because in y and z dirs
    // the "ondir" direction will sometimes never hit the grid.
    //
    temp = (*bpt);
    //
    // If nbc==2, then we have two pointers to two sheets of cells,
    // and we want both of them to point to the first on-grid cell,
    // so we use "isedge" to move 1 or 2 cells on-grid to get the
    // on-grid value.
    //
    for (int v = 0; v > (*bpt)->isedge; v--) {
      temp = grid->NextPt(temp, ondir);
    }

    //
    // Now set inflow data to be the first on-grid cell's values.
    //
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->P[v] = temp->P[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->Ph[v] = temp->P[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->dU[v] = 0.0;

    (*bpt)->isdomain = false;
    // rep.printVec("Setting inflow boundary values:",temp->P,par.nvar);
    ct++;
    ++bpt;
  } while (bpt != b->data.end());

  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = temp->P[v];

  if (ct != b->data.size())
    rep.error("BC_assign_INFLOW: missed some cells!", ct - b->data.size());

  return 0;
}

// ##################################################################
// ##################################################################

int inflow_bc::BC_update_INFLOW(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int,
    const int)
{
  //
  // Inflow means BC is constant at the initial value of the
  // neighbouring edge cell, so we set dU=0.
  //
  list<cell *>::iterator c = b->data.begin();
  for (c = b->data.begin(); c != b->data.end(); ++c) {
    for (int v = 0; v < par.nvar; v++)
      (*c)->dU[v] = 0.0;
    for (int v = 0; v < par.nvar; v++)
      (*c)->P[v] = b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*c)->Ph[v] = b->refval[v];
    // rep.printVec("inflow boundary values:",(*c)->P,par.nvar);
  }  // all cells.
  return 0;
}

// ##################################################################
// ##################################################################
