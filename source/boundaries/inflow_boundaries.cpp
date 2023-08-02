/// \file inflow_boundaries.cpp
/// \brief Class definitions for inflow boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

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
    spdlog::error("{}: {}", "BC_assign_INFLOW: empty boundary data", b->itype);
  }
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  double edge = 0.0;
  enum axes a = XX;
  switch (b->dir) {
    case XN:
      edge = grid->Xmin(XX);
      a    = XX;
      break;
    case XP:
      edge = grid->Xmax(XX);
      a    = XX;
      break;
    case YN:
      edge = grid->Xmin(YY);
      a    = YY;
      break;
    case YP:
      edge = grid->Xmax(YY);
      a    = YY;
      break;
    case ZN:
      edge = grid->Xmin(ZZ);
      a    = ZZ;
      break;
    case ZP:
      edge = grid->Xmax(ZZ);
      a    = ZZ;
      break;
    case NO:
    default:
      spdlog::error("bad direction in inflow boundary");
      exit_pion(1);
      break;
  }

  list<cell *>::iterator bpt = b->data.begin();
  cell *temp;
  std::array<double, MAX_DIM> dpos;

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
    // for (int v = 0; v > (*bpt)->isedge; v--) {
    //  temp = grid->NextPt(*temp, ondir);
    //}
    bool brk = false;
    do {
      temp = grid->NextPt(*temp, ondir);
      CI.get_dpos(*temp, dpos);
      switch (b->dir) {
        case XN:
        case YN:
        case ZN:
          if (dpos[a] > edge) brk = true;
          break;
        case XP:
        case YP:
        case ZP:
          if (dpos[a] < edge) brk = true;
          break;
        case NO:
        default:
          spdlog::error("bad direction in inflow boundary 2");
          exit_pion(1);
          break;
      }
    } while (!brk);

    //
    // Now set inflow data to be the first on-grid cell's values.
    //
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->P[v] = temp->P[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->Ph[v] = temp->P[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->dU[v] = 0.0;

    (*bpt)->npt      = temp;
    (*bpt)->isdomain = false;
    // rep.printVec("Setting inflow boundary values:",temp->P,par.nvar);
    ct++;
    ++bpt;
  } while (bpt != b->data.end());

  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = temp->P[v];

#ifndef NDEBUG
  std::vector<pion_flt> tt(par.nvar);
  for (int v = 0; v < par.nvar; v++)
    tt[v] = temp->P[v];
  spdlog::info("INFLOW REF: {}", tt);
#endif  // NDEBUG

  if (ct != b->data.size()) {
    spdlog::error(
        "{}: {}", "BC_assign_INFLOW: missed some cells!", ct - b->data.size());
    exit_pion(1);
  }

  return 0;
}

// ##################################################################
// ##################################################################

int inflow_bc::BC_update_INFLOW(
    class SimParams &par,   ///< pointer to simulation parameters
    class GridBaseClass *,  ///< pointer to grid.
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
    for (int v = 0; v < par.nvar; v++) {
      (*c)->P[v] = (*c)->npt->P[v];
      //(*c)->P[v] = b->refval[v];
    }
    for (int v = 0; v < par.nvar; v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v];
      //(*c)->Ph[v] = b->refval[v];
    }
    // rep.printVec("inflow boundary values:",(*c)->P,par.nvar);
  }  // all cells.
  return 0;
}

// ##################################################################
// ##################################################################
