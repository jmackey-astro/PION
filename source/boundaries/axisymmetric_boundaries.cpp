/// \file axisymmetric_boundaries.cpp
/// \brief Class definitions for axisymmetric boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/axisymmetric_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

// ##################################################################
// ##################################################################

int axisymmetric_bc::BC_assign_AXISYMMETRIC(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b)
{
  if (b->data.empty()) {
    spdlog::error("BC_assign_AXISYMMETRIC empty boundary data {}", b->itype);
    exit_pion(1);
  }

  if (b->dir != RNcyl) {
    spdlog::error(
        "Axisymmetric boundary not at R=0 : {}", static_cast<int>(b->dir));
    exit_pion(2);
  }

  //
  // set reference state so that it is mostly zeros but has some +/-1
  // entries to flip the signs of the velocity and B-field (as
  // appropriate).
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
    for (int v = 0; v < par.nvar; v++)
      b->refval[v] = 1.0;
    //
    // velocity:
    //
    b->refval[VY] = -1.0;  // reflect the radial component
    b->refval[VZ] = -1.0;  // rotate the theta component

    //
    // B-field:
    //
    if (par.eqntype == EQMHD || par.eqntype == EQGLM || par.eqntype == EQFCD) {
      b->refval[BY] = -1.0;  // reflect the radial component
      b->refval[BZ] = -1.0;  // rotate the theta component
    }                        // Setting up reference value.
  }                          // if we needed to set up refval.

  //
  // Now go through each of the boundary points and assign values
  // to them, multiplying the relevant entries by -1.
  //
  list<cell *>::iterator bpt = b->data.begin();
  cell *temp                 = 0;
  unsigned int ct            = 0;
  // cout <<"\t\t**** PRINTING CELLS FOR BOUNDARY DIR = "<<b->dir<<"
  // ****\n\n";
  do {
    temp = (*bpt);
    for (int v = -(*bpt)->isedge; v > (*bpt)->isedge + 1; v--) {
      temp = grid->NextPt(*temp, b->ondir);
    }
    if (!temp) {
      spdlog::error("Got lost assigning axisymmetric bcs. {}", temp->id);
      exit_pion(3);
    }
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->P[v] = temp->P[v] * b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->Ph[v] = temp->Ph[v] * b->refval[v];
    for (int v = 0; v < par.nvar; v++)
      (*bpt)->dU[v] = 0.0;
    (*bpt)->npt = temp;
    // CI.print_cell((*bpt));
    // axisymmetric boundary is a conformal mapping of the domain, so
    // isdomain should be true??
    (*bpt)->isdomain = false;
    ++bpt;
    ct++;
  } while (bpt != b->data.end());

  if (ct != b->data.size()) {
    spdlog::error(
        "{}: {}", "BC_assign_AXISYMMETRIC: missed some cells! {} {}", ct,
        b->data.size());
    exit_pion(4);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int axisymmetric_bc::BC_update_AXISYMMETRIC(
    class SimParams &par,   ///< pointer to simulation parameters
    class GridBaseClass *,  ///< pointer to grid.
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
