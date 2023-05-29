/// \file radiative_shock_test.cpp
/// \brief Class definitions for boundary conditions for a radiative
///  shock test problem.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/radiative_shock_test.h"
#include "tools/mem_manage.h"
using namespace std;



// ##################################################################
// ##################################################################



int radiative_shock_test::BC_assign_RADSHOCK(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
  enum direction offdir = b->dir;
  enum direction ondir  = b->ondir;

  if (!b->data.empty()) {
    spdlog::error(
        "{} : {}", "BC_assign_RADSHOCK: not empty boundary data", b->itype);
    exit(1);
  }
  // set reference state
  b->refval = mem.myalloc(b->refval, par.nvar);
  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = 1.0;
  cell *c = grid->FirstPt();
  // This boundary requires density to be <10x initial density.
  b->refval[RO] = c->P[RO] * 10.0;
  int ct        = 0;
  if (par.ndim == 1) {
    do {
      c->P[RO] = c->Ph[RO] = min(c->P[RO], b->refval[RO]);
      b->data.push_back(c);
      ct++;
    } while (CI.get_dpos(*(c = grid->NextPt(*c)), XX) <= par.Xmax[XX] / 256.0);
  }
  spdlog::info(
      "{} {} {}", "******** Added ", ct,
      " cells to RADSHOCK internal boundary.");
  return 0;
}



// ##################################################################
// ##################################################################



int radiative_shock_test::BC_update_RADSHOCK(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  //
  // same routine as for outflow, except multiply v_n,B_n by -1.
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
