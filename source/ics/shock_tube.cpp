/// \file shock_tube.cc
///
/// File for setting up shock-tube problems.
/// Modified:
///  - 2009-12-22 JM:
///    Changed interface position for 2D rotated shock tubes.
///    Started working on 2D AW test.
///  - 2009-12-23 JM:
///    Changed interface position again, so that it's always at the
///    same position at y=0.5.
///    Added check for negative angle, which is converted to arctan(0.5).
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.05 JM: Added pion_flt datatype.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "equations/eqns_hydro_adiabatic.h"
#include "equations/eqns_mhd_adiabatic.h"
#include "ics/icgen.h"
#include "ics/icgen_base.h"
#include <sstream>

// ##################################################################
// ##################################################################

IC_shocktube::IC_shocktube()
{
    preshock = postshock = 0;
    gg                   = 0;
    rp                   = 0;
    ndim = coords = eqns = -1;
    angleXY = angleXZ = 0;
    number            = -10;
    return;
}

// ##################################################################
// ##################################################################

IC_shocktube::~IC_shocktube()
{
    if (preshock) {
        delete[] preshock;
        preshock = 0;
    }
    if (postshock) {
        delete[] postshock;
        postshock = 0;
    }
    return;
}

// ##################################################################
// ##################################################################

int IC_shocktube::setup_data(
    class ReadParams* rrp,    ///< pointer to parameter list.
    class GridBaseClass* ggg  ///< pointer to grid
)
{
    int err = 0;

    ICsetup_base::gg = ggg;
    if (!gg) rep.error("null pointer to grid!", ggg);

    ICsetup_base::rp = rrp;
    if (!rp) rep.error("null pointer to ReadParams", rp);

    IC_shocktube::ndim = SimPM->ndim;
    if (ndim != 1 && ndim != 2 && ndim != 3)
        rep.error("Shock-Tube problem must be 1d or 2d or 3d", ndim);
    IC_shocktube::coords = SimPM->coord_sys;
    if (coords != COORD_CRT && coords != COORD_CYL)
        rep.error("Bad coord sys", coords);
    IC_shocktube::eqns = SimPM->eqntype;
    if (eqns == EQEUL)
        eqns = 1;
    else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
        eqns = 2;
    else
        rep.error("Bad equations", eqns);

    // initialise pre and post shock vectors to zero.
    IC_shocktube::preshock  = 0;
    IC_shocktube::postshock = 0;
    preshock                = new double[SimPM->nvar];
    postshock               = new double[SimPM->nvar];
    if (!preshock || !postshock)
        rep.error("malloc pre/post shock vecs", preshock);
    for (int v = 0; v < SimPM->nvar; v++)
        preshock[v] = postshock[v] = 0.0;

    cout << "Note that shocktube B-field values are in units where ";
    cout << "magnetic pressure is 0.5*B^2, i.e. no sqrt(4pi).\n";

    IC_shocktube::shockpos = 0.0;  // initial value

    // set pre-/post- shock and cloud states.
    string seek, str;
    // Shock Tube number:
    seek = "STnumber";
    str  = rp->find_parameter(seek);
    if (str == "") {
        seek = "RIEMANN";  // try this legacy shock tube flag...
        str  = rp->find_parameter(seek);
        if (str == "") rep.error("didn't find parameter", seek);
    }
    IC_shocktube::number = atoi(str.c_str());

    seek = "STshockpos";
    str  = rp->find_parameter(seek);
    if (str == "")
        IC_shocktube::shockpos = 0.0;
    else
        IC_shocktube::shockpos = atof(str.c_str());
    cout << "shock pos:" << shockpos << endl;

    //
    // angle in XY plane (degrees) (angle normal makes with x-axis, projected
    // onto xy plane) Read in as an integer angle in degrees and converted to
    // radians. A negative value is interpreted as a request for arctan(0.5)
    //
    seek = "STangleXY";
    str  = rp->find_parameter(seek);
    if (str == "")
        IC_shocktube::angleXY = 0.0;
    else
        IC_shocktube::angleXY = static_cast<double>(atoi(str.c_str()));
    angleXY *= M_PI / 180.0;
    if (angleXY < 0.0) angleXY = atan(0.5);

    // angle in XZ plane (degrees) (angle normal makes with x-axis, projected
    // onto xz plane)
    seek = "STangleXZ";
    str  = rp->find_parameter(seek);
    if (str == "")
        IC_shocktube::angleXZ = 0.0;
    else
        IC_shocktube::angleXZ = static_cast<double>(atoi(str.c_str()));
    angleXZ *= M_PI / 180.0;

    if (number <= 0) {
        // pre- and post-shock state vectors
        ostringstream temp;
        string v;

        v = "RO";
        temp.str("");
        temp << "STprevec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::preshock[RO] = atof(str.c_str());
        else
            IC_shocktube::preshock[RO] = -1.0e99;
        temp.str("");
        temp << "STpostvec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::postshock[RO] = atof(str.c_str());
        else
            IC_shocktube::postshock[RO] = -1.0e99;

        v = "PG";
        temp.str("");
        temp << "STprevec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::preshock[PG] = atof(str.c_str());
        else
            IC_shocktube::preshock[PG] = -1.0e99;
        temp.str("");
        temp << "STpostvec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::postshock[PG] = atof(str.c_str());
        else
            IC_shocktube::postshock[PG] = -1.0e99;

        v = "VX";
        temp.str("");
        temp << "STprevec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::preshock[VX] = atof(str.c_str());
        else
            IC_shocktube::preshock[VX] = -1.0e99;
        temp.str("");
        temp << "STpostvec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::postshock[VX] = atof(str.c_str());
        else
            IC_shocktube::postshock[VX] = -1.0e99;

        v = "VY";
        temp.str("");
        temp << "STprevec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::preshock[VY] = atof(str.c_str());
        else
            IC_shocktube::preshock[VY] = -1.0e99;
        temp.str("");
        temp << "STpostvec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::postshock[VY] = atof(str.c_str());
        else
            IC_shocktube::postshock[VY] = -1.0e99;

        v = "VZ";
        temp.str("");
        temp << "STprevec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::preshock[VZ] = atof(str.c_str());
        else
            IC_shocktube::preshock[VZ] = -1.0e99;
        temp.str("");
        temp << "STpostvec" << v;
        seek = temp.str();
        str  = rp->find_parameter(seek);
        if (str != "")
            IC_shocktube::postshock[VZ] = atof(str.c_str());
        else
            IC_shocktube::postshock[VZ] = -1.0e99;

        if (eqns == 2) {  // mhd sim
            v = "BX";
            temp.str("");
            temp << "STprevec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::preshock[BX] = atof(str.c_str());
            else
                IC_shocktube::preshock[BX] = -1.0e99;
            temp.str("");
            temp << "STpostvec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::postshock[BX] = atof(str.c_str());
            else
                IC_shocktube::postshock[BX] = -1.0e99;

            v = "BY";
            temp.str("");
            temp << "STprevec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::preshock[BY] = atof(str.c_str());
            else
                IC_shocktube::preshock[BY] = -1.0e99;
            temp.str("");
            temp << "STpostvec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::postshock[BY] = atof(str.c_str());
            else
                IC_shocktube::postshock[BY] = -1.0e99;

            v = "BZ";
            temp.str("");
            temp << "STprevec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::preshock[BZ] = atof(str.c_str());
            else
                IC_shocktube::preshock[BZ] = -1.0e99;
            temp.str("");
            temp << "STpostvec" << v;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::postshock[BZ] = atof(str.c_str());
            else
                IC_shocktube::postshock[BZ] = -1.0e99;

            if (SimPM->eqntype == EQGLM) preshock[SI] = postshock[SI] = 0.;
        }  // if mhd vars

        // tracer variables
        for (int t = 0; t < SimPM->ntracer; t++) {
            temp.str("");
            temp << "STprevecTR" << t;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::preshock[t + SimPM->ftr] = atof(str.c_str());
            else
                IC_shocktube::preshock[t + SimPM->ftr] = -1.0e99;
            temp.str("");
            temp << "STpostvecTR" << t;
            seek = temp.str();
            str  = rp->find_parameter(seek);
            if (str != "")
                IC_shocktube::preshock[t + SimPM->ftr] = atof(str.c_str());
            else
                IC_shocktube::preshock[t + SimPM->ftr] = -1.0e99;
        }
    }  // if number<=0, read in vectors
    else
        err += get_riemann_ics(number, postshock, preshock, &shockpos);
    if (err) rep.error("get_riemann_ics error", err);
    cout << "shockpos: " << shockpos << endl;

    rep.printVec("preshock", preshock, SimPM->nvar);
    rep.printVec("postshock", postshock, SimPM->nvar);
    // now we have left and right state vectors, and the shock position and
    // orientation.
    err += assign_data(postshock, preshock, shockpos);
    return err;
}

// ##################################################################
// ##################################################################

int IC_shocktube::assign_data(
    double* l_in,     ///< input left state
    double* r_in,     ///< input right state
    double interface  ///< location of interface
)
{
    int nvar = gg->Nvar();
    if (ndim < 1 || ndim > 3) rep.error("Bad ndim in setupNDWave", ndim);

    //
    // Cell data could be float or double, depending on pion_flt, so we
    // copy the input vectors to an array of the right type.
    //
    pion_flt left[nvar], right[nvar];
    for (int v = 0; v < nvar; v++) {
        left[v]  = l_in[v];
        right[v] = r_in[v];
    }

    class eqns_base* eqn = 0;
    if (eqns == 1)
        eqn = new class eqns_Euler(nvar);
    else if (eqns == 2)
        eqn = new class eqns_mhd_ideal(nvar);
    else
        rep.error("Bad equations in assign_data()", eqns);
    if (!eqn) rep.error("Bad equations in assign_data()", eqns);

    cout << "Setting up a " << ndim
         << "-D simulation with a ShockTube problem.\n";
    cout << "Initial shock position is at x=" << interface << " and angle is t="
         << angleXY * 180.0 / M_PI << " degrees.\n";
    rep.printVec("Left : ", left, nvar);
    rep.printVec("Right: ", right, nvar);

    // preshock state vector.
    class cell* cpt = gg->FirstPt();
    double dx       = 2.0 * gg->DX();
    if (ndim == 1) {
        do {
            // Set values of primitive variables.
            if (CI.get_dpos(cpt, XX) < (interface - dx))
                for (int v = 0; v < nvar; v++)
                    cpt->P[v] = left[v];
            else if (CI.get_dpos(cpt, XX) > (interface + dx))
                for (int v = 0; v < nvar; v++)
                    cpt->P[v] = right[v];
            else {
                for (int v = 0; v < nvar; v++)
                    cpt->P[v] = 0.5 * (left[v] + right[v])
                                + 0.5 * (CI.get_dpos(cpt, XX) - interface)
                                      * (right[v] - left[v]) / dx;
            }
        } while ((cpt = gg->NextPt(cpt)) != NULL);
    }
    else if (ndim == 2 || ndim == 3) {
        //
        // Rotate state vectors about the z-axis; setting position of interface.
        //
        eqn->rotateXY(left, angleXY);
        eqn->rotateXY(right, angleXY);
        if (ndim == 3) rep.error("Please re-code me for 3D shock tubes!", ndim);
        double tt   = tan(angleXY);  // st=sin(angleXY), ct=cos(angleXY);
        double xmax = interface + (0.5 - SimPM->Xmin[YY]) * tt;
        double dpos[ndim], x0;
        do {
            //
            // Get cell position, and position of discontinuity at this y-value.
            //
            CI.get_dpos(cpt, dpos);
            x0 = xmax - (dpos[YY] - SimPM->Xmin[YY]) * tt;
            if (dpos[XX] <= x0) {
                for (int v = 0; v < nvar; v++)
                    cpt->P[v] = cpt->Ph[v] = left[v];
                // This case puts in a smoothing of the shock over one cell
                // instead of the sharp jump. Not sure if it's a good idea --
                // the total pressure dips and comes back up again, which seems
                // wrong.  P \propto B^2 which is why it happens, but it may be
                // a bad idea. else if ( (dpos[YY]-SimPM->dx/2.) <
                // ymax*(1.-dpos[XX]/xmax) ) { cout <<"boundary point!\n";
                // printCell(cpt);
                //  for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = cpt->Ph[v] =
                //  (left[v]+right[v])/2.;
            }
            else
                for (int v = 0; v < nvar; v++)
                    cpt->P[v] = cpt->Ph[v] = right[v];
        } while ((cpt = gg->NextPt(cpt)) != 0);

        // Now enforce divB=0 if needed.
        /*    if (nvar>=8) {
              cell *c=gg->FirstPt(); cout <<"Bx="<<c->P[BX]<<endl; do {
              CI.get_dpos(cpt,dpos);
              c->Ph[BX] =0.;
              c->Ph[BY] = c->P[BZ]*dpos[XX];
              c->Ph[BZ] = c->P[BX]*dpos[YY] - c->P[BY]*dpos[XX];
              }  while( (c=gg->NextPt(c))!=0);

              double temp[3];
              c=gg->FirstPt(); cout <<"Bx="<<c->P[BX]<<endl; do {
              if (!c->isedge) {
              gg->VecCurl(c,2,1,temp);
              cout <<"curl A=["<<temp[0]<<", "<<temp[1]<<", "<<temp[2]<<"
           ]"<<endl; c->P[BX] = temp[0]; c->P[BY] = temp[1]; c->P[BZ] = temp[2];
              }
              }  while( (c=gg->NextPt(c))!=0);
              c=gg->FirstPt(); cout <<"Bx="<<c->P[BX]<<endl;
              }
              * */
    }  // ndim==2,3
    else {
        rep.error("bad ndim", ndim);
    }

    //
    // If we have the Alfven wave test, we need to do extra work
    // since it isn't a discontinuity.
    //
    if (IC_shocktube::number == 8) {
        cout << "Alfven wave: switching to periodic boundaries.\n";
        if (ndim == 1) {
            SimPM->BC_XN   = "periodic";
            SimPM->BC_XP   = "periodic";
            SimPM->BC_Nint = 0;
            double len = 0.3, dpos[ndim], amp = 1.0;

            cpt = gg->FirstPt();
            do {
                CI.get_dpos(cpt, dpos);
                //
                // If we are in [0.5,0.8] then add in a rotation by 2pi
                //
                if (dpos[XX] >= interface && dpos[XX] < (interface + len)) {
                    cpt->P[VY] = cpt->Ph[VY] =
                        amp * cos(2.0 * M_PI * (dpos[XX] - interface) / len);
                    cpt->P[BY] = cpt->Ph[BY] =
                        amp * cos(2.0 * M_PI * (dpos[XX] - interface) / len);
                    cpt->P[VZ] = cpt->Ph[VZ] =
                        amp
                        * (1.0
                           + sin(2.0 * M_PI * (dpos[XX] - interface) / len));
                    cpt->P[BZ] = cpt->Ph[BZ] =
                        amp * sin(2.0 * M_PI * (dpos[XX] - interface) / len);
                }
            } while ((cpt = gg->NextPt(cpt)) != 0);
        }  // 1D
        else if (ndim == 2) {
            // rep.error("AW test not set up in 2D yet.",ndim);
            cout << "Alfven Wave test in 2D -- note not the same as 1D!!!\n";
            SimPM->BC_XN   = "periodic";
            SimPM->BC_XP   = "periodic";
            SimPM->BC_YN   = "periodic";
            SimPM->BC_YP   = "periodic";
            SimPM->BC_Nint = 0;
            double theta = atan(2.0), dpos[ndim], amp = 0.1;

            cpt = gg->FirstPt();
            do {
                CI.get_dpos(cpt, dpos);
                //
                // This uses the test from Stone's code test page.
                //
                cpt->P[RO] = 1.0;
                cpt->P[PG] = 0.1;
                cpt->P[VX] = 0.0;
                cpt->P[VY] =
                    amp
                    * sin(2.0 * M_PI
                          * (dpos[XX] * cos(theta) + dpos[YY] * sin(theta)));
                cpt->P[VZ] =
                    amp
                    * cos(2.0 * M_PI
                          * (dpos[XX] * cos(theta) + dpos[YY] * sin(theta)));
                cpt->P[BX] = 1.0;
                cpt->P[BY] = cpt->P[VY];
                cpt->P[BZ] = cpt->P[VZ];
                eqn->rotateXY(cpt->P, theta);
            } while ((cpt = gg->NextPt(cpt)) != 0);
        }  // 2D
        else
            rep.error("AW test not set up in 3D yet.", ndim);
    }

    delete eqn;
    eqn = 0;
    return 0;
}

// ##################################################################
// ##################################################################

int IC_shocktube::get_riemann_ics(int sw, double* l, double* r, double* xm)
{
    // for these, we are dimensionless, so no microphysics allowed!
    SimPM->EP.raytracing        = 0;
    SimPM->EP.cooling           = 0;
    SimPM->EP.chemistry         = 0;
    SimPM->EP.coll_ionisation   = 0;
    SimPM->EP.rad_recombination = 0;
    SimPM->EP.phot_ionisation   = 0;
    SimPM->EP.update_erg        = false;

    // These are Toro's five tests on p.225 of his book.
    switch (sw) {
        case 1:
            /** case 1: Toro's test no.1 on p.225 of his book.\n*/
            l[RO] = 1.;
            l[PG] = 1.;
            l[VX] = 0.75;
            l[VY] = l[VZ] = 0.;
            r[RO]         = 0.125;
            r[PG]         = 0.1;
            r[VX]         = 0.0;
            r[VY] = r[VZ] = 0.;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.3;
            SimPM->gamma      = 1.4;
            SimPM->finishtime = 0.2;
            break;
        case 2:
            /** case 2: Toro's test no.2 on p.225 of his book.\n*/
            l[RO] = 1.;
            l[PG] = 0.4;
            l[VX] = -2.0;
            l[VY] = l[VZ] = 0.;
            r[RO]         = 1.;
            r[PG]         = 0.4;
            r[VX]         = 2.0;
            r[VY] = r[VZ] = 0.;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.5;
            SimPM->gamma      = 1.4;
            SimPM->finishtime = 0.15;
            break;
        case 3:
            /** case 3: Toro's test no.3 on p.225 of his book.\n*/
            l[RO] = 1.;
            l[PG] = 1000.;
            l[VX] = 0.0;
            l[VY] = l[VZ] = 0.;
            r[RO]         = 1.;
            r[PG]         = 0.01;
            r[VX]         = 0.0;
            r[VY] = r[VZ] = 0.;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.5;
            SimPM->gamma      = 1.4;
            SimPM->finishtime = 0.012;
            break;
        case 4:
            /** case 4: Toro's test no.4 on p.225 of his book.\n*/
            l[RO] = 5.99924;
            l[PG] = 460.894;
            l[VX] = 19.5975;
            l[VY] = l[VZ] = 0.;
            r[RO]         = 5.99242;
            r[PG]         = 46.0950;
            r[VX]         = -6.19633;
            r[VY] = r[VZ] = 0.;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.4;
            SimPM->gamma      = 1.4;
            SimPM->finishtime = 0.035;
            break;
        case 5:
            /** case 5: Toro's test no.5 on p.225 of his book.\n*/
            cout << "Setting up Toro5 shock tube problem\n";
            l[RO] = 1.;
            l[PG] = 1000.;
            l[VX] = -19.59745;
            l[VY] = l[VZ] = 0.;
            r[RO]         = 1.;
            r[PG]         = 0.01;
            r[VX]         = -19.59745;
            r[VY] = r[VZ] = 0.;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.8;
            SimPM->gamma      = 1.4;
            SimPM->finishtime = 0.012;
            break;
        // other wierd cases
        // Left                     [rho,v,p] = [0.604543, 1.876, 1.69426 ]
        // Right                    [rho,v,p] = [1, 2, 1 ]
        // This fools the linear solver...
        case 6:
            /** case 6: Slightly difficult case with an almost stationary
             * rarefaction.*/
            l[RO] = 0.604543;
            l[PG] = 1.69426;
            l[VX] = 1.876;
            l[VY] = l[VZ] = 0.4;
            r[RO]         = 1;
            r[PG]         = 1;
            r[VX]         = 2;
            r[VY] = r[VZ] = 0.5;
            if (SimPM->eqntype == 2 || SimPM->eqntype == EQGLM) {
                l[BX] = l[BY] = l[BZ] = r[BX] = r[BY] = r[BZ] = 0.;
            }
            *xm               = 0.3;
            SimPM->gamma      = 5. / 3.;
            SimPM->finishtime = 0.15;
            break;
            // From here on I am using MHD test cases, so it won't work for
            // hydro.
        case 7:
            /** case 7: Sam Falle's test 'BW', the Brio and Wu problem.\n */
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = 1.;
            l[PG] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 0.75;
            l[BY]                 = 1.;
            l[BZ]                 = 0.;
            r[RO]                 = 0.125;
            r[PG]                 = 0.1;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 0.75;
            r[BY]                 = -1.;
            r[BZ]                 = 0.;
            *xm                   = 0.5;
            SimPM->gamma          = 2.0;
            SimPM->finishtime     = 0.12;
            break;
        case 8:
            /// case 8: Sam Falle's test 'AW', an Alfven wave.  We put in the
            /// rotation between left and right states later -- this is not a
            /// discontinuity.
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = l[PG] = 1.;
            l[VX]         = 0.;
            l[VY] = l[VZ] = 1.;
            l[BX] = l[BY] = 1.;
            l[BZ]         = 0.;
            r[RO] = r[PG] = 1.;
            r[VX]         = 0.;
            r[VY] = r[VZ] = 1.;
            r[BX] = r[BY]     = 1.;
            r[BZ]             = 0.;
            *xm               = 0.5;
            SimPM->gamma      = 5. / 3.;
            SimPM->finishtime = 5.0;
            break;
        case 9:
            /** case 9: Sam Falle's test 'FS', a fast shock.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            // l[RO]=3.; l[PG]=49./3.; l[VX]=1.0-sqrt(3.0); l[VY]=-4./3.;
            // l[VZ]=0.; l[BX]=3.; l[BY]=sqrt(16.0/3.0); l[BZ]=0.0; //1.;
            // r[RO]=1.; r[PG]=1.; r[VX]=-sqrt(12.0)-sqrt(3)+1.0; r[VY]=0.;
            // r[VZ]=0.; r[BX]=3.; r[BY]=0.; r[BZ]=0.;
            l[RO]             = 3.;
            l[PG]             = 16.33;
            l[VX]             = -0.732;
            l[VY]             = -1.333;
            l[VZ]             = 0.;
            l[BX]             = 3.;
            l[BY]             = 2.309;
            l[BZ]             = 0.0;  // 1.;
            r[RO]             = 1.;
            r[PG]             = 1.;
            r[VX]             = -4.196;
            r[VY]             = 0.;
            r[VZ]             = 0.;
            r[BX]             = 3.;
            r[BY]             = 0.;
            r[BZ]             = 0.;
            *xm               = 0.3;
            SimPM->gamma      = 5. / 3.;
            SimPM->finishtime = 0.4;
            break;
        case 10:
            /** case 10: Sam Falle's test 'SS', a slow shock.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = 1.368;
            l[PG] = 1.769;
            l[VX] = 0.269;
            l[VY] = 1.;
            l[VZ] = 0.;
            l[BX] = 1.;
            l[BY] = l[BZ] = 0.;
            r[RO]         = 1.;
            r[PG]         = 1.;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 1.;
            r[BY]                 = 1.;
            r[BZ]                 = 0.;
            *xm                   = 0.3;
            SimPM->gamma          = 5. / 3.;
            SimPM->finishtime     = 0.5;
            break;
        case 11:
            /** case 11: Sam Falle's test 'FR', a fast rarefaction.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = 1.;
            l[PG] = 2.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 1.;
            l[BY]                 = 3.;
            l[BZ]                 = 0.;
            r[RO]                 = 0.2641;
            r[PG]                 = 0.2175;
            r[VX]                 = 3.6;
            r[VY]                 = -2.551;
            r[VZ]                 = 0.;
            r[BX]                 = 1.;
            r[BY] = r[BZ]     = 0.;
            *xm               = 0.5;
            SimPM->gamma      = 5. / 3.;
            SimPM->finishtime = 0.1;
            break;
        case 12:
            /** case 12: Sam Falle's test 'SR', a slow rarefaction.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = 1.;
            l[PG] = 2.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 1.;
            l[BY] = l[BZ]     = 0.;
            r[RO]             = 0.2;
            r[PG]             = 0.1368;
            r[VX]             = 1.186;
            r[VY]             = 2.967;
            r[VZ]             = 0.;
            r[BX]             = 1.;
            r[BY]             = 1.6405;
            r[BZ]             = 0.;
            *xm               = 0.5;
            SimPM->gamma      = 5. / 3.;
            SimPM->finishtime = 0.3;
            break;
        case 13:
            /** case 13: Sam Falle's test 'OFS', an oblique fast shock.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            l[RO] = 1.;
            l[PG] = 1.;
            l[VX] = 6.505;
            l[VY] = 1.;
            l[VZ] = 0.;
            l[BX] = l[BY] = l[BZ] = 1.;
            r[RO]                 = 3.;
            r[PG]                 = 20.268;
            r[VX]                 = 2.169;
            r[VY]                 = 1.331;
            r[VZ]                 = 0.331;
            r[BX]                 = 1.;
            r[BY]                 = 3.153;
            r[BZ]                 = 3.153;
            *xm                   = 0.5;
            SimPM->gamma          = 5. / 3.;
            SimPM->finishtime     = 0.15;
            break;
        case 14:
            /** case 14: Trivial case, only call it if you intend to add noise
             * later.\n*/
            for (int v = 0; v < SimPM->nvar; v++) {
                l[v] = r[v] = 1.;
                *xm         = 0.5;
            }
            break;
            /** Ryu and Jones (1995) Shock Tube tests, 1a-5b follow. [Ryu \&
             * Jones, 1995, ApJ,442,228].\n */
        case 15:
            /** case 15: Ryu and Jones test 1a.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = 10.;
            l[VY] = l[VZ] = 0.;
            l[BX] = l[BY] = 5. / sqrt(4 * M_PI);
            l[BZ]         = 0.;
            l[PG]         = 20.;
            r[RO]         = 1.;
            r[VX]         = -10.;
            r[VY] = r[VZ] = 0.;
            r[BX] = r[BY] = 5. / sqrt(4 * M_PI);
            r[BZ]         = 0.;
            r[PG]         = 1.;
            *xm           = 0.5;
            break;
        case 16:
            /** case 16: Ryu and Jones test 1b.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 3. / sqrt(4 * M_PI);
            l[BY]                 = 5. / sqrt(4 * M_PI);
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.1;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 3. / sqrt(4 * M_PI);
            r[BY]                 = 2. / sqrt(4 * M_PI);
            r[BZ]                 = 0.;
            r[PG]                 = 10.;
            *xm                   = 0.5;
            break;
        case 17:
            /** case 17: Ryu and Jones test 2a.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.08;
            l[VX] = 1.2;
            l[VY] = 0.01;
            l[VZ] = 0.5;
            l[BX] = 2. / sqrt(4. * M_PI);
            l[BY] = 3.6 / sqrt(4. * M_PI);
            l[BZ] = 2. / sqrt(4. * M_PI);
            l[PG] = 0.95;
            r[RO] = 1.;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 2. / sqrt(4. * M_PI);
            r[BY]                 = 4. / sqrt(4. * M_PI);
            r[BZ]                 = 2. / sqrt(4. * M_PI);
            r[PG]                 = 1.;
            *xm                   = 0.5;
            break;
        case 18:
            /** case 18: Ryu and Jones test 2b.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 3. / sqrt(4. * M_PI);
            l[BY]                 = 6. / sqrt(4. * M_PI);
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.1;
            r[VX]                 = 0.;
            r[VY]                 = 2.;
            r[VZ]                 = 1.;
            r[BX]                 = 3. / sqrt(4. * M_PI);
            r[BY]                 = 1. / sqrt(4. * M_PI);
            r[BZ]                 = 0.;
            r[PG]                 = 10.;
            *xm                   = 0.5;
            break;
        case 19:
            /** case 19: Ryu and Jones test 3a.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 0.1;
            l[VX] = 50.;
            l[VY] = l[VZ] = 0.;
            l[BX]         = 0.;
            l[BY]         = -1. / sqrt(4. * M_PI);
            l[BZ]         = -2. / sqrt(4. * M_PI);
            l[PG]         = 0.4;
            r[RO]         = 0.1;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 0.;
            r[BY]                 = 1. / sqrt(4. * M_PI);
            r[BZ]                 = 2. / sqrt(4. * M_PI);
            r[PG]                 = 0.2;
            *xm                   = 0.5;
            break;
        case 20:
            /** case 20: Ryu and Jones test 3b.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = -1.;
            l[VY] = l[VZ] = 0.;
            l[BX]         = 0.;
            l[BY]         = 1.;
            l[BZ]         = 0.;
            l[PG]         = 1.;
            r[RO]         = 1.;
            r[VX]         = 1.;
            r[VY] = r[VZ] = 0.;
            r[BX]         = 0.;
            r[BY]         = 1.;
            r[BZ]         = 0.;
            r[PG]         = 1.;
            *xm           = 0.5;
            break;
        case 21:
            /** case 21: Ryu and Jones test 4a.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 1.;
            l[BY]                 = 1.;
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.2;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 1.;
            r[BY]                 = 0.;
            r[BZ]                 = 0.;
            r[PG]                 = 0.1;
            *xm                   = 0.5;
            break;
        case 22:
            /** case 22: Ryu and Jones test 4b.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 0.4;
            l[VX] = -0.66991;
            l[VY] = 0.98263;
            l[VZ] = 0.;
            l[BX] = 1.3;
            l[BY] = 0.0025293;
            l[BZ] = 0.;
            l[PG] = 0.52467;
            r[RO] = 1.;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 1.3;
            r[BY]                 = 1.;
            r[BZ]                 = 0.;
            r[PG]                 = 1.;
            *xm                   = 0.5;
            break;
        case 23:
            /** case 23: Ryu and Jones test 4c.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 0.65;
            l[VX] = 0.667;
            l[VY] = -0.257;
            l[VZ] = 0.;
            l[BX] = 0.75;
            l[BY] = 0.55;
            l[BZ] = 0.;
            l[PG] = 0.5;
            r[RO] = 1.;
            r[VX] = 0.4;
            r[VY] = -0.94;
            r[VZ] = 0.;
            r[BX] = 0.75;
            r[BY] = 0.;
            r[BZ] = 0.;
            r[PG] = 0.75;
            *xm   = 0.5;
            break;
        case 24:
            /** case 24: Ryu and Jones test 4d.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 0.7;
            l[BY]                 = 0.;
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.3;
            r[VX] = r[VY] = 0.;
            r[VZ]         = 1.;
            r[BX]         = 0.7;
            r[BY]         = 1.;
            r[BZ]         = 0.;
            r[PG]         = 0.2;
            *xm           = 0.5;
            break;
        case 25:
            /** case 25: Ryu and Jones test 5a.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 0.75;
            l[BY]                 = 1.;
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.125;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 0.75;
            r[BY]                 = -1.;
            r[BZ]                 = 0.;
            r[PG]                 = 0.1;
            *xm                   = 0.5;
            break;
        case 26:
            /** case 26: Ryu and Jones test 5b.\n*/
            if (eqns != 2) {
                cerr
                    << "(IC_shocktube::get_riemann_ics) Not using MHD but asking for "
                       "MHD test problem. Exiting.\n";
                return (1);
            }
            SimPM->gamma = 5. / 3.;
            cout << "\t Forcing gamma=5/3 for test problem.\n";
            l[RO] = 1.;
            l[VX] = l[VY] = l[VZ] = 0.;
            l[BX]                 = 1.3;
            l[BY]                 = 1.;
            l[BZ]                 = 0.;
            l[PG]                 = 1.;
            r[RO]                 = 0.4;
            r[VX] = r[VY] = r[VZ] = 0.;
            r[BX]                 = 1.3;
            r[BY]                 = -1.;
            r[BZ]                 = 0.;
            r[PG]                 = 0.4;
            *xm                   = 0.5;
            break;
        default:
            cout << "Error: only know 26 tests, but sw!={1,..,26}" << endl;
            return (1);
    }
    // tracers: not much use for these dimensionless problems.
    for (int t = 0; t < SimPM->ntracer; t++) {
        l[SimPM->ftr + t] = 1.0;
        r[SimPM->ftr + t] = -1.0;
    }

    cout << "(IC_shocktube::get_riemann_ics) Got test number: " << sw << endl;
    return (0);
}

// ##################################################################
// ##################################################################
