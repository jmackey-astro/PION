/// \file oneway_out_boundaries.cpp
/// \brief Class definitions for oneway_out boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/oneway_out_boundaries.h"
using namespace std;

// ##################################################################
// ##################################################################

int oneway_out_bc::BC_assign_ONEWAY_OUT(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    boundary_data* b)
{
    //
    // The setup for this is identical to outflow.
    //
    int err = BC_assign_OUTFLOW(par, grid, b);
    return err;
}

// ##################################################################
// ##################################################################

int oneway_out_bc::BC_update_ONEWAY_OUT(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    struct boundary_data* b,
    const int cstep,
    const int maxstep)
{
    //
    // Outflow or Absorbing BCs; boundary cells are same as edge cells.
    // This is zeroth order outflow bcs.
    // For the one-way boundary we set the normal velocity to be zero if
    // the flow wants to be onto the domain.
    //

    //
    // First get the normal velocity component, and whether the offdir is
    // positive or negative.
    //
    enum direction offdir = b->dir;
    enum primitive Vnorm  = RO;
    int norm_sign         = 0;
    switch (offdir) {
        case XN:
        case XP:
            Vnorm = VX;
            break;
        case YN:
        case YP:
            Vnorm = VY;
            break;
        case ZN:
        case ZP:
            Vnorm = VZ;
            break;
        default:
            rep.error("bad dir", offdir);
    }
    if (Vnorm == RO)
        rep.error("Failed to set normal velocity component", Vnorm);

    switch (offdir) {
        case XN:
        case YN:
        case ZN:
            norm_sign = -1;
            break;
        case XP:
        case YP:
        case ZP:
            norm_sign = 1;
            break;
        default:
            rep.error("bad dir", offdir);
    }

    //
    // Now run through all cells in the boundary
    //
    list<cell*>::iterator c = b->data.begin();
    cell* gc;
    for (c = b->data.begin(); c != b->data.end(); ++c) {
        (*c)->isdomain = false;
        gc             = (*c)->npt;
        //
        // exactly same routine as for periodic.
        // ONEWAY_OUT: overwrite the normal velocity if it is inflow:
        //
        for (int v = 0; v < par.nvar; v++)
            (*c)->Ph[v] = (*c)->npt->Ph[v];
        for (int v = 0; v < par.nvar; v++)
            (*c)->dU[v] = 0.;
        (*c)->Ph[Vnorm] =
            norm_sign
            * max(static_cast<pion_flt>(0.0), (*c)->Ph[Vnorm] * norm_sign);
        if (cstep == maxstep) {
            for (int v = 0; v < par.nvar; v++)
                (*c)->P[v] = (*c)->npt->P[v];
            (*c)->P[Vnorm] =
                norm_sign
                * max(static_cast<pion_flt>(0.0), (*c)->P[Vnorm] * norm_sign);
        }

        //
        // The GLM boundary is somewhat different, because I found
        // that zero-gradient didn't work well (Mackey & Lim, 2011,
        // MNRAS,412,2079).  So we switch the sign instead.  Same as
        // periodic BCs.
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
