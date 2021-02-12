/// \file periodic_boundaries.cpp
/// \brief Class definitions for periodic boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/periodic_boundaries.h"
using namespace std;

// ##################################################################
// ##################################################################

int periodic_bc::BC_assign_PERIODIC(
    class SimParams& par,       ///< pointer to simulation parameters
    const int,                  ///< level in grid hierarchy
    class GridBaseClass* grid,  ///< pointer to grid.
    boundary_data* b)
{
    enum direction ondir = b->ondir;
    if (b->data.empty())
        rep.error("BC_assign_PERIODIC: empty boundary data", b->itype);

    // loop through all cells in boundary, and make the cell's "npt"
    // pointer point to the corresponding cell wrapped around on the
    // other side of the grid.  Assign data from that cell too.
    list<cell*>::iterator bpt = b->data.begin();
    cell* temp;
    unsigned int ct = 0;
    do {
        // boundary cells are part of the domain, so set isdomain
        (*bpt)->isdomain = true;
        //
        // Go across the grid NG[baxis] cells, so that we map onto the
        // cell on the other side of the grid.
        //
        temp = (*bpt);
        for (int v = 0; v < grid->NG(b->baxis); v++) {
            temp = grid->NextPt(temp, ondir);
        }
        for (int v = 0; v < par.nvar; v++)
            (*bpt)->P[v] = temp->P[v];
        for (int v = 0; v < par.nvar; v++)
            (*bpt)->Ph[v] = temp->P[v];
        for (int v = 0; v < par.nvar; v++)
            (*bpt)->dU[v] = 0.0;
        (*bpt)->npt = temp;
        //
        // increment counters.
        //
        ct++;
        ++bpt;
    } while (bpt != b->data.end());

    if (ct != b->data.size())
        rep.error(
            "BC_assign_PERIODIC: missed some cells!", ct - b->data.size());
    return 0;
}

// ##################################################################
// ##################################################################

int periodic_bc::BC_update_PERIODIC(
    class SimParams& par,       ///< pointer to simulation parameters
    const int l,                ///< level in grid hierarchy
    class GridBaseClass* grid,  ///< pointer to grid.
    struct boundary_data* b,
    const int cstep,
    const int maxstep)
{
    // cout <<"updating "<<b->data.size()<<" cells for periodic on level
    // "<<l<<"\n";
    list<cell*>::iterator c = b->data.begin();
    for (c = b->data.begin(); c != b->data.end(); ++c) {
        // cout <<"cell and npt: ";
        // rep.printVec("cell",(*c)->pos,par.ndim);
        // rep.printVec("npt",(*c)->npt->pos,par.ndim);
        for (int v = 0; v < par.nvar; v++)
            (*c)->Ph[v] = (*c)->npt->Ph[v];
        for (int v = 0; v < par.nvar; v++)
            (*c)->dU[v] = 0.;
        if (cstep == maxstep) {
            for (int v = 0; v < par.nvar; v++)
                (*c)->P[v] = (*c)->npt->P[v];
        }
    }  // all cells.
    return 0;
}

// ##################################################################
// ##################################################################
