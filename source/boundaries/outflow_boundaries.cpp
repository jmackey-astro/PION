/// \file outflow_boundaries.cpp
/// \brief Class definitions for outflow boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/outflow_boundaries.h"
using namespace std;

// ##################################################################
// ##################################################################

int outflow_bc::BC_assign_OUTFLOW(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    boundary_data* b)
{
    //
    // Zero order extrapolation, if edge cell is index k, and boundary cells
    // k+1 and k+2, then P_{k+1} = P_{k+2} = P_{k}
    // First order extrapolation would be P_{k+1} = 2P_{k} - P_{k-1} but is
    // said to have stability problems sometimes (LeVeque, S.7.2.1, p131-132)
    //
    enum direction ondir = b->ondir;

    if (b->data.empty()) {
        rep.error("BC_assign_OUTFLOW: empty boundary data", b->itype);
    }
    list<cell*>::iterator bpt = b->data.begin();
    cell* temp                = 0;
    unsigned int ct           = 0;  // counter (for accounting).

    //
    // loop through all boundary cells, set npt to point to the grid
    // cell where we get the data from.
    //
    do {
        //
        // Find the cell to point to.  This should be the first cell on
        // the other side of the boundary, located at x[baxis]=bloc.
        // Can't just use the ->isgd property because in y and z dirs
        // the "ondir" direction will sometimes never hit the grid.
        //
        temp = (*bpt);
        //
        // If nbc>1, then we have pointers to two or more sheets of cells
        // and we want both of them to point to the first on-grid cell,
        // so we use "isedge" to move 1 or 2 cells on-grid to get the
        // on-grid value.
        //
        for (int v = 0; v > (*bpt)->isedge; v--) {
            // CI.print_cell(temp);
            // cout <<"ondir="<<ondir<<"\n";
            temp = grid->NextPt(temp, ondir);
        }

        for (int v = 0; v < par.nvar; v++)
            (*bpt)->P[v] = temp->P[v];
        for (int v = 0; v < par.nvar; v++)
            (*bpt)->Ph[v] = temp->P[v];
        (*bpt)->npt = temp;
        // CI.print_cell((*bpt));
        (*bpt)->isdomain = false;

        //
        // The GLM boundary is somewhat different, because
        // zero-gradient didn't work well (Mackey & Lim, 2011,
        // MNRAS,412,2079).  So we switch the sign instead.
        //
#ifdef GLM_ZERO_BOUNDARY
        if (par.eqntype == EQGLM) {
            (*bpt)->P[SI] = (*bpt)->Ph[SI] = 0.0;
        }
#endif  // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
        if (par.eqntype == EQGLM) {
            if (par.ndim == 1) {
                rep.error(
                    "Psi outflow boundary condition doesn't work for 1D!", 99);
            }
            // get to nth cell on grid.
            for (int v = (*bpt)->isedge + 1; v < 0; v++)
                temp = grid->NextPt(temp, ondir);
            (*bpt)->P[SI]  = -temp->P[SI];
            (*bpt)->Ph[SI] = -temp->Ph[SI];
        }
#endif  // GLM_NEGATIVE_BOUNDARY
        ct++;
        ++bpt;
    } while (bpt != b->data.end());

    if (ct != b->data.size()) {
        rep.error("BC_assign_OUTFLOW: missed some cells!", ct - b->data.size());
    }
    return 0;
}

// ##################################################################
// ##################################################################

int outflow_bc::BC_update_OUTFLOW(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    struct boundary_data* b,
    const int cstep,
    const int maxstep)
{
    //
    // Outflow or Absorbing BCs; boundary cells are same as edge cells.
    // This is zeroth order outflow bcs.
    //
    list<cell*>::iterator c = b->data.begin();
    cell* gc;

    for (c = b->data.begin(); c != b->data.end(); ++c) {
        //
        // gc is the on-grid cell.
        //
        gc = (*c)->npt;
        for (int v = 0; v < par.nvar; v++)
            (*c)->Ph[v] = gc->Ph[v];
        for (int v = 0; v < par.nvar; v++)
            (*c)->dU[v] = 0.;
        if (cstep == maxstep) {
            for (int v = 0; v < par.nvar; v++)
                (*c)->P[v] = gc->P[v];
        }

        //
        // The GLM boundary is somewhat different, because I found
        // that zero-gradient didn't work well (Mackey & Lim, 2011,
        // MNRAS,412,2079).  So we switch the sign instead.
        //
#ifdef GLM_ZERO_BOUNDARY
        if (par.eqntype == EQGLM) {
            (*c)->P[SI]  = 0.0;
            (*c)->Ph[SI] = 0.0;
        }
#endif  // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
        if (par.eqntype == EQGLM) {
            for (int v = (*c)->isedge + 1; v < 0; v++)
                gc = grid->NextPt(gc, b->ondir);
            (*c)->P[SI]  = -gc->P[SI];
            (*c)->Ph[SI] = -gc->Ph[SI];
        }
#endif  // GLM_NEGATIVE_BOUNDARY

    }  // all cells.
    return 0;
}

// ##################################################################
// ##################################################################
