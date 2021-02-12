/// \file jet_boundaries.cpp
/// \brief Class definitions for jet boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/jet_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

//#define SOFTJET
#define JETPROFILE

// ##################################################################
// ##################################################################

int jet_bc::BC_assign_JETBC(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    boundary_data* b)
{
    if (!JP.jetic) {
        rep.error("BC_assign_JETBC: not a jet simulation!", JP.jetic);
    }
    if (b->dir != NO) {
        rep.error("BC_assign_JETBC: boundary is not an internal one!", b->dir);
    }
    cell* c    = grid->FirstPt();
    cell *temp = 0, *cy = 0;
    int ct    = 0;
    int ctot  = 0;
    int maxnv = 0;

    //
    // Set the physical radius of jet.
    //
    double jr = JP.jetradius * grid->DX();
#ifdef TESTING
    cout << "jetrad=" << JP.jetradius << " dx=" << grid->DX() << "\n";
#endif  // TESTING

    //
    // Assign reference values, containing Jet parameters:
    //
    if (!b->refval) {
        b->refval = mem.myalloc(b->refval, par.nvar);
    }

    if (par.eqntype == EQEUL || par.eqntype == EQEUL_ISO
        || par.eqntype == EQEUL_EINT || par.eqntype == EQMHD
        || par.eqntype == EQGLM || par.eqntype == EQFCD) {
        rep.printVec("JetState", JP.jetstate, MAX_NVAR);
        b->refval[RO] = JP.jetstate[RO];
        b->refval[PG] = JP.jetstate[PG];
        b->refval[VX] = JP.jetstate[VX];
        b->refval[VY] = JP.jetstate[VY];
        b->refval[VZ] = JP.jetstate[VZ];
        maxnv         = 5;
    }
    else {
        rep.error("BC_assign_JETBC: bad equation type", par.eqntype);
    }

    if (par.eqntype == EQMHD || par.eqntype == EQGLM || par.eqntype == EQFCD) {
        // jetstate[BX] is the axial field, and jetstate[BY] is the toroidal
        // field
        if (par.ndim == 2) {
            b->refval[BX] = JP.jetstate[BX];
            b->refval[BY] = 0.0;  // this is the divergence component, must be 0
            b->refval[BZ] = JP.jetstate[BY];
        }
        else {
            rep.error("Need to code B-field within jet for 3D", par.ndim);
        }
        maxnv = 8;
    }

    if (par.eqntype == EQGLM) {
        b->refval[SI] = JP.jetstate[SI];
        maxnv         = 9;
    }

    for (int v = maxnv; v < par.nvar; v++) {
        b->refval[v] = JP.jetstate[v];
    }

    if (!b->data.empty()) {
        rep.error("BC_assign_JETBC: boundary data exists!", b->itype);
    }

    //
    // Simplest jet -- zero opening angle.  Check that this domain
    // is at the edge of the full simulation domain, if not then we
    // don't need to add any cells so skip this whole loop.
    //
    if (JP.jetic == 1 && pconst.equalD(grid->Xmin(XX), par.Xmin[XX])) {
        //
        // Axi-symmetry first -- this is relatively easy to set up.
        //
        if (par.ndim == 2 && par.coord_sys == COORD_CYL) {
            c = grid->FirstPt();
            do {
                temp = c;
                while ((temp = grid->NextPt(temp, XN)) != 0) {
                    for (int v = 0; v < par.nvar; v++)
                        temp->P[v] = b->refval[v];
                    for (int v = 0; v < par.nvar; v++)
                        temp->Ph[v] = b->refval[v];
#ifdef SOFTJET
                    //
                    // jetradius is in number of cells, jr is in physical units.
                    // Jet centre is along Z_cyl axis, centred on origin.
                    // Set the outer 25% of the jet radius to have a linearly
                    // decreasing velocity from Vjet to zero in that range.
                    //
                    temp->P[VX] =
                        b->refval[VX]
                        * min(1., 4.0 - 4.0 * CI.get_dpos(temp, YY) / jr);
                    temp->Ph[VX] = temp->P[VX];
#endif  // SOFTJET
#ifdef JETPROFILE
                    double r     = fabs(CI.get_dpos(temp, YY));
                    double rm    = 0.9 * jr;
                    double B_phi = b->refval[BZ];
                    double p0    = b->refval[PG];
                    double rm_jr = rm / jr;
                    double r_jr  = r / jr;
                    double r_rm  = r / rm;
                    double pm    = p0
                                - pow(B_phi, 2) / pow(1 - rm_jr, 2)
                                      * (3 * (1 - rm_jr) - (1 - pow(rm_jr, 2))
                                         + log(rm_jr));
                    if (r == rm) {
                        temp->P[PG] = temp->Ph[PG] =
                            p0;  // surface gas pressure
                    }
                    if (r < rm) {
                        temp->P[BZ] = temp->Ph[BZ] = B_phi * r_rm;
                        temp->P[PG]                = temp->Ph[PG] =
                            pow(B_phi, 2) * (1 - pow(r_rm, 2)) + pm;
                    }
                    else if (r < jr) {
                        temp->P[BZ] = temp->Ph[BZ] =
                            B_phi * (jr - r) / (jr - rm);
                        temp->P[PG] = temp->Ph[PG] =
                            p0
                            + pow(B_phi, 2) / pow(1 - rm_jr, 2)
                                  * (3 * (1 - r_jr) - (1 - pow(r_jr, 2))
                                     + log(r_jr));
                    }
                    else {
                        temp->P[BZ] = temp->Ph[BZ] = 0;
                    }
#endif  // JETPROFILE
                    b->data.push_back(temp);
                    ctot++;
                    if (temp->isgd) {
                        rep.error("Looking for Boundary cells! setupjet", temp);
                    }
                }
                ct++;
            } while ((c = grid->NextPt(c, YP)) && ct < JP.jetradius);
            if (ct != JP.jetradius) {
                rep.error("Not enough cells for jet", ct);
            }
            cout << "Got " << ctot << " Cells in total for jet boundary.\n";
        }  // 2D Axial Symmetry

        //
        // 3D now, more difficult.
        //
        else if (par.ndim == 3 && par.coord_sys == COORD_CRT) {
            double dist = 0.0;
            c           = grid->FirstPt();
            //
            // Also, the jet will come in at the centre of the XN boundary,
            // which must be the origin.
            //
            do {  // loop over ZZ axis
                cy = c;
                do {  // loop over YY axis
                    dist = sqrt(
                        CI.get_dpos(cy, YY) * CI.get_dpos(cy, YY)
                        + CI.get_dpos(cy, ZZ) * CI.get_dpos(cy, ZZ));
                    //
                    // if dist <= jr, then we are within the jet inflow, and we
                    // add the cells to the boundary.
                    //
                    if (dist <= jr) {
                        temp = cy;
                        while ((temp = grid->NextPt(temp, XN)) != 0) {
                            for (int v = 0; v < par.nvar; v++)
                                temp->P[v] = b->refval[v];
                            for (int v = 0; v < par.nvar; v++)
                                temp->Ph[v] = b->refval[v];
#ifdef SOFTJET
                            // Set the outer 25% of the jet radius to have a
                            // linearly decreasing velocity from Vjet to zero in
                            // that range.
                            temp->P[VX] =
                                b->refval[VX] * min(1., 4.0 - 4.0 * dist / jr);
                            temp->Ph[VX] = temp->P[VX];
#endif  // SOFTJET
                            b->data.push_back(temp);
                            ctot++;
                            if (temp->isgd) {
                                rep.error(
                                    "Looking for Boundary cells! setupjet",
                                    temp);
                            }
                        }
                    }  // if within jet radius
                } while (
                    (cy = grid->NextPt(cy, YP)));  // through cells on YY axis.
            } while ((c = grid->NextPt(c, ZP)));   // through cells on ZZ axis.
                                                   //      BC_printBCdata(b);
        }                                          // 3D Cartesian

        else {
            rep.error(
                "Only know how to set up jet in 2Dcyl or 3Dcart", par.ndim);
        }

    }  // jetic==1

    return 0;
}

// ##################################################################
// ##################################################################

int jet_bc::BC_update_JETBC(
    class SimParams& par,       ///< pointer to simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    struct boundary_data* b,
    const int,
    const int)
{
#ifdef SOFTJET
    double dist = 0.0;
    double jr   = JP.jetradius * grid->DX();
#endif  // SOFTJET

    list<cell*>::iterator c = b->data.begin();
    for (c = b->data.begin(); c != b->data.end(); ++c) {
        for (int v = 0; v < par.nvar; v++)
            (*c)->dU[v] = 0.0;
        for (int v = 0; v < par.nvar; v++)
            (*c)->P[v] = b->refval[v];
        for (int v = 0; v < par.nvar; v++)
            (*c)->Ph[v] = b->refval[v];
#ifdef SOFTJET
        dist = 0.0;
        if (par.ndim == 2)
            dist = CI.get_dpos(*c, YY);
        else if (par.ndim == 3) {
            dist = sqrt(
                CI.get_dpos(*c, YY) * CI.get_dpos(*c, YY)
                + CI.get_dpos(*c, ZZ) * CI.get_dpos(*c, ZZ));
        }
        else
            rep.error("Jet BC, but not 2d or 3d!!!", par.ndim);
        (*c)->P[VX]  = b->refval[VX] * min(1., 4 - 4.0 * dist / jr);
        (*c)->Ph[VX] = (*c)->P[VX];
#endif  // SOFTJET
    }   // all cells.
    return 0;
}

// ##################################################################
// ##################################################################
