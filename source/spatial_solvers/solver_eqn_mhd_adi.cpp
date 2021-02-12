///
/// \file solver_eqn_mhd_adi.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 written, based on old solver.cc classes.
///
/// Solver for the adiabatic ideal MHD Equations.
/// Calculates flux via either Lax-Friedrichs or Riemann solver.
/// Adds viscosity if asked for, and tracks flux of N passive tracers.
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux
/// functions).
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.27 JM: Removed riemann_base references.  Added
///   function_id identifiers.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2011.04.15 JM: Change in UtoP() for tracers (to try to correct for
/// negative density!).
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2018.04.14 JM: Moved flux solver to FV_solver

// change
#include "coord_sys/VectorOps.h"

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"

#include "solver_eqn_mhd_adi.h"
using namespace std;

// *****************************************
// ***** FV SOLVER MHD Ideal Adiabatic *****
// *****************************************

// ##################################################################
// ##################################################################

FV_solver_mhd_ideal_adi::FV_solver_mhd_ideal_adi(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt* state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv, state, gam),
    Riemann_Roe_MHD_CV(nv, gam),
    HLLD_MHD(nv, gam),
    VectorOps_Cart(nd)
{
#ifdef TESTING
    cout << "::FV_solver_mhd_ideal_adi() constructor.\n";
    // cout <<"::FV_solver_mhd_ideal_adi() gamma = "<<eq_gamma<<"\n";
#endif
    max_speed = 0.0;
    negPGct = negROct = 0;
    return;
}

// ##################################################################
// ##################################################################

FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
    cout << "::~FV_solver_mhd_ideal_adi ...starting.\n";
#endif  // FUNCTION_ID

#ifdef TESTING
    cout << "FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi() destructor.\n";
#endif

#ifdef FUNCTION_ID
    cout << "::~FV_solver_mhd_ideal_adi ...returning.\n";
#endif  // FUNCTION_ID
    return;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::inviscid_flux(
    class SimParams& par,       ///< simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid
    const double dx,            ///< cell-size dx (for LF method)
    class cell* Cl,             ///< Left state cell pointer
    class cell* Cr,             ///< Right state cell pointer
    const pion_flt* Pl,         ///< Left Primitive vector.
    const pion_flt* Pr,         ///< Right Primitive vector.
    pion_flt* flux,             ///< Resultant Flux vector.
    pion_flt* pstar,            ///< State vector at interface.
    const int solve_flag,       ///< Solver to use
    const double eq_gamma       ///< Gas constant gamma.
)
{
#ifdef TESTING
    // Check input density and pressure are 'reasonably large'
    if (Pl[eqRO] < TINYVALUE || Pl[eqPG] < TINYVALUE || Pr[eqRO] < TINYVALUE
        || Pr[eqPG] < TINYVALUE) {
        rep.printVec("left ", Pl, eq_nvar);
        rep.printVec("right", Pr, eq_nvar);
        rep.error(
            "FV_solver_mhd_ideal_adi::calculate_flux() Density/Pressure too small",
            Pl[eqRO]);
    }
#endif  // TESTING
    int err = 0;
    double ustar[eq_nvar];
    for (int v = 0; v < eq_nvar; v++)
        ustar[v] = 0.0;
    for (int v = 0; v < eq_nvar; v++)
        flux[v] = 0.0;
    for (int v = 0; v < eq_nvar; v++)
        pstar[v] = 0.0;

    // Choose which Solver to use:
    if (solve_flag == FLUX_LF) {
        // Lax-Friedrichs Method, so just get the flux
        err += get_LaxFriedrichs_flux(Pl, Pr, flux, dx, eq_gamma);
        for (int v = 0; v < eq_nvar; v++)
            pstar[v] = 0.5 * (Pl[v] + Pr[v]);
    }

    else if (solve_flag == FLUX_RSroe) {
        // Roe Flux solver in conserved variables (Cargo and Gallice, 1997).
        // This is the symmetric version which sums over all waves.
        err += MHD_Roe_CV_flux_solver_symmetric(
            Pl, Pr, eq_gamma, HC_etamax, pstar, flux);

        //
        // If we get an error, try the linear solver:
        //
        if (err) {
            // err=0;
            cout << "ROE SOLVER FAILED -- TRYING FKJ98 SOLVER. err=" << err
                 << "\n";
            err = JMs_riemann_solve(Pl, Pr, pstar, 1, eq_gamma);
            PtoFlux(pstar, flux, eq_gamma);
        }
    }

    else if (
        solve_flag == FLUX_RSlinear || solve_flag == FLUX_RSexact
        || solve_flag == FLUX_RShybrid) {
        // JM's linear Riemann Solver -- Falle et al. (1998) with the Roe
        // and Balsara (1996) eigenvector normalisation.
        err += JMs_riemann_solve(Pl, Pr, pstar, solve_flag, eq_gamma);
        PtoFlux(pstar, flux, eq_gamma);
    }

    // HLLD/HLL solver, including compressive motion check and
    // strong-gradient zones check (Migone et al. 2011 )
    else if (solve_flag == FLUX_RS_HLLD) {
        double DivVl = CI.get_DivV(Cl);
        double DivVr = CI.get_DivV(Cr);
        double Gradl = CI.get_MagGradP(Cl);
        double Gradr = CI.get_MagGradP(Cr);
        if ((DivVl < 0. && Gradl > 5.) || (DivVr < 0. && Gradr > 5.)) {
            // compressive motion & strong-gradient zones check
            // Migone et al 2012
            // HLL solver -- Miyoshi and Kusano (2005) (m05)
            err += MHD_HLL_flux_solver(Pl, Pr, eq_gamma, flux, ustar);
        }
        else {
            // HLLD solver -- Miyoshi and Kusano (2005) (m05)
            err += MHD_HLLD_flux_solver(Pl, Pr, eq_gamma, flux, ustar);
        }
        rep.errorTest("HLL/HLLD Flux", 0, err);
        err = UtoP(ustar, pstar, par.EP.MinTemperature, eq_gamma);
        rep.errorTest("HLL/HLLD UtoP", 0, err);
    }

    // HLL solver, diffusive 2 wave solver (Migone et al. 2011 )
    else if (solve_flag == FLUX_RS_HLL) {
        err += MHD_HLL_flux_solver(Pl, Pr, eq_gamma, flux, ustar);
        rep.errorTest("HLL Flux", 0, err);
        err = UtoP(ustar, pstar, par.EP.MinTemperature, eq_gamma);
        rep.errorTest("HLL UtoP", 0, err);
    }

    else {
        rep.error("what sort of flux solver do you mean???", solve_flag);
    }

    return err;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::AVFalle(
    const pion_flt* Pleft,
    const pion_flt* Pright,
    const pion_flt* Pstar,
    pion_flt* flux,
    const double eta,   ///< already set as FV_etav
    const double gamma  ///< already set as eq_gamma
)
{
    /// \section Equations
    /// The equations are as follows (for flux along x-direction):
    /// \f[ F[p_x] = F[p_x] - \eta \rho^{*} c^{*}_f (v_{xR}-v_{xL})  \,,\f]
    /// \f[ F[p_y] = F[p_y] - \eta \rho^{*} c^{*}_f (v_{yR}-v_{yL})  \,,\f]
    /// \f[ F[p_z] = F[p_z] - \eta \rho^{*} c^{*}_f (v_{zR}-v_{zL})  \,,\f]
    /// \f[ F[e]   = F[e]   - \eta \rho^{*} c^{*}_f \left[ v^{*}_x
    /// (v_{xR}-v_{xL})
    /// + v^{*}_y (v_{yR}-v_{yL}) + v^{*}_z (v_{zR}-v_{zL}) \right]\,.\f] This
    /// follows from assuming a 1D problem, with non-zero bulk viscosity, and a
    /// non-zero shear viscosity.  The bulk viscosity gives the viscosity in the
    /// x-direction, opposing stretching and compression.  The shear viscosity
    /// gives the y and z terms, opposing slipping across boundaries.
    ///
    /// The identification of the velocities and density with the
    /// resultant state from the Riemann Solver is a bit ad-hoc, and could
    /// easily be changed to using the mean of the left and right states, or
    /// something else. (See my notes on artificial viscosity).
    ///
    /// Andy says Sam also has similar terms for the Magnetic field, although
    /// it isn't in his paper.  Andy uses this too.
    /// \f[ f[B_y] = f[B_y] - \eta c^{*}_f (B_{yR}-B_{yL}) \;, \f]
    /// \f[ f[B_z] = f[B_z] - \eta c^{*}_f (B_{zR}-B_{zL}) \;, \f]
    /// with associated terms added to energy flux:
    /// \f[ F[e]   = F[e]   - \eta c^{*}_f \left[ B^{*}_z (B_{zR}-B_{zL})
    /// +B^{*}_y (B_{yR}-B_{yL})\right] \;. \f]
    ///
    /// Andy uses the mean of the left and right state velocities/field to
    /// calculate the energy flux; I am using the starred state velocity/field.
    /// Not sure if that is bad, or if it will make any discernible difference
    /// at all. After some testing, it makes very little difference.  It may
    /// help with some crazy cases, but I haven't encountered them yet.
    ///
    /// I have switched from using \f$c^{*}_f\f$ to using the fast speed from
    /// the mean vector of the left and right states.  This is because the
    /// Riemann Solver can return negative density, which then gets set to a
    /// density floor value, but the B-field is still normal, so the fast speed
    /// blows up.  If I use the mean state though, it represents a typical speed
    /// between the two states, so it is in principle no less realistic.
    ///
    double prefactor = cfast_components(
                           0.5 * (Pleft[eqRO] + Pright[eqRO]),
                           0.5 * (Pleft[eqPG] + Pright[eqPG]),
                           0.5 * (Pleft[eqBX] + Pright[eqBX]),
                           0.5 * (Pleft[eqBY] + Pright[eqBY]),
                           0.5 * (Pleft[eqBZ] + Pright[eqBZ]), eq_gamma)
                       * FV_etav * Pstar[eqRO];

    //
    // Momentum flux
    //
    double momvisc = prefactor * (Pright[eqVX] - Pleft[eqVX]);  //*4./3.;
    double ergvisc = momvisc * Pstar[eqVX];
    flux[eqMMX] -= momvisc;
    momvisc = prefactor * (Pright[eqVY] - Pleft[eqVY]);
    flux[eqMMY] -= momvisc;
    ergvisc += momvisc * Pstar[eqVY];
    momvisc = prefactor * (Pright[eqVZ] - Pleft[eqVZ]);
    flux[eqMMZ] -= momvisc;
    ergvisc += momvisc * Pstar[eqVZ];

    //
    // Magnetic field flux.
    //
    prefactor *= FV_etaB / (FV_etav * Pstar[eqRO]);
    momvisc = prefactor * (Pright[eqBY] - Pleft[eqBY]);
    flux[eqBBY] -= momvisc;
    ergvisc += momvisc * Pstar[eqBY];
    momvisc = prefactor * (Pright[eqBZ] - Pleft[eqBZ]);
    flux[eqBBZ] -= momvisc;
    ergvisc += momvisc * Pstar[eqBZ];
    flux[eqERG] -= ergvisc;
    return (0);
}

// ##################################################################
// ##################################################################

void FV_solver_mhd_ideal_adi::PtoU(
    const pion_flt* p, pion_flt* u, const double g)
{
    eqns_mhd_ideal::PtoU(p, u, g);
    for (int t = 0; t < FV_ntr; t++)
        u[eqTR[t]] = p[eqTR[t]] * p[eqRO];
    return;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::UtoP(
    const pion_flt* u,
    pion_flt* p,
    const double MinTemp,  ///< Min Temperature allowed on grid.
    const double g)
{
    for (int t = 0; t < FV_ntr; t++)
        p[eqTR[t]] = u[eqTR[t]] / u[eqRHO];
    int err = eqns_mhd_ideal::UtoP(u, p, MinTemp, g);
    return err;
}

// ##################################################################
// ##################################################################

void FV_solver_mhd_ideal_adi::PUtoFlux(
    const pion_flt* p, const pion_flt* u, pion_flt* f)
{
    eqns_mhd_ideal::PUtoFlux(p, u, f);
    for (int t = 0; t < FV_ntr; t++)
        f[eqTR[t]] = p[eqTR[t]] * f[eqRHO];
    return;
}

// ##################################################################
// ##################################################################

void FV_solver_mhd_ideal_adi::UtoFlux(
    const pion_flt* u, pion_flt* f, const double g)
{
    eqns_mhd_ideal::UtoFlux(u, f, g);
    for (int t = 0; t < FV_ntr; t++)
        f[eqTR[t]] = u[eqTR[t]] * f[eqRHO] / u[eqRHO];
    return;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::dU_Cell(
    class GridBaseClass* grid,
    cell* c,                // Current cell.
    const axes d,           // Which axis we are looking along.
    const pion_flt* fn,     // Negative direction flux.
    const pion_flt* fp,     // Positive direction flux.
    const pion_flt* slope,  // slope vector for cell c.
    const int ooa,          // spatial order of accuracy.
    const double dx,        // cell length dx.
    const double dt         // cell TimeStep, dt.
)
{
    pion_flt u1[eq_nvar];
    // This calculates -dF/dx
    int err = DivStateVectorComponent(c, grid, d, eq_nvar, fn, fp, u1);
    geometric_source(c, d, slope, ooa, dx, u1);

    for (int v = 0; v < eq_nvar; v++)
        c->dU[v] += FV_dt * u1[v];
    return (err);
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::MHDsource(
    class GridBaseClass* grid,  ///< pointer to grid.
    class cell* Cl,             ///< pointer to cell of left state
    class cell* Cr,             ///< pointer to cell of right state
    pion_flt* Pl,               ///< left edge state
    pion_flt* Pr,               ///< right edge state
    const axes d,               ///< Which axis we are looking along.
    enum direction pos,         ///< positive normal direction
    enum direction neg,         ///< negative normal direction
    const double dt             ///< timestep dt
)
{
    // The Powell source terms from Powell's paper (1999)
    // called by time_integrator::dynamics_dU_column()
    double dx   = grid->DX();
    double bm   = 0.5 * (Cl->Ph[eqBX] + Cr->Ph[eqBX]);
    double uB_l = Cl->Ph[eqBX] * Cl->Ph[eqVX] + Cl->Ph[eqBY] * Cl->Ph[eqVY]
                  + Cl->Ph[eqBZ] * Cl->Ph[eqVZ];
    double uB_r = Cr->Ph[eqBX] * Cr->Ph[eqVX] + Cr->Ph[eqBY] * Cr->Ph[eqVY]
                  + Cr->Ph[eqBZ] * Cr->Ph[eqVZ];
    pion_flt Powell_l[eq_nvar], Powell_r[eq_nvar];
    for (int v = 0; v < eq_nvar; v++) {
        Powell_l[v] = Powell_r[v] = 0.0;
    }
    Powell_l[eqMMX] = Cl->Ph[eqBX];
    Powell_l[eqMMY] = Cl->Ph[eqBY];
    Powell_l[eqMMZ] = Cl->Ph[eqBZ];
    Powell_l[eqERG] = uB_l;
    Powell_l[eqBBX] = Cl->Ph[eqVX];
    Powell_l[eqBBY] = Cl->Ph[eqVY];
    Powell_l[eqBBZ] = Cl->Ph[eqVZ];

    Powell_r[eqMMX] = Cr->Ph[eqBX];
    Powell_r[eqMMY] = Cr->Ph[eqBY];
    Powell_r[eqMMZ] = Cr->Ph[eqBZ];
    Powell_r[eqERG] = uB_r;
    Powell_r[eqBBX] = Cr->Ph[eqVX];
    Powell_r[eqBBY] = Cr->Ph[eqVY];
    Powell_r[eqBBZ] = Cr->Ph[eqVZ];

    for (int v = 0; v < eq_nvar; v++) {
        Cl->dU[v] -= dt * bm * (Powell_l[v]) / dx;
        Cr->dU[v] += dt * bm * (Powell_r[v]) / dx;
    }
    return 0;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_ideal_adi::CellAdvanceTime(
    class cell* c,         // cell to update.
    const pion_flt* Pin,   // Initial State Vector.
    pion_flt* dU,          // Update vector dU
    pion_flt* Pf,          // Final state vector (can be same as initial vec.).
    pion_flt* dE,          // Tracks dE (to correct for negative pressure)
    const double,          // gas EOS gamma.
    const double MinTemp,  ///< Min Temperature allowed on grid.
    const double           // Cell timestep dt.
)
{
    pion_flt u1[eq_nvar], u2[eq_nvar];
    pion_flt Pintermediate[eq_nvar];
    pion_flt corrector[eq_nvar];
    //
    // First convert from Primitive to Conserved Variables
    //
    if (MP) {
        MP->sCMA(corrector, Pin);
        for (int t = 0; t < eq_nvar; t++)
            Pintermediate[t] = Pin[t] * corrector[t];
        PtoU(Pintermediate, u1, eq_gamma);
    }
    else {
        PtoU(Pin, u1, eq_gamma);
    }

    // Now add dU[] to U[], and change back to primitive variables.
    // This can give negative pressures, so check for that and fix it
    // if needed.
    for (int v = 0; v < eq_nvar; v++)
        u1[v] += dU[v];

    if (u1[RHO] < 0.0) {
        cout << "celladvancetime, negative density. rho=" << u1[RHO] << "\n";
        CI.print_cell(c);
    }

    if (UtoP(u1, Pf, MinTemp, eq_gamma) != 0) {
        cout << "(FV_solver_mhd_ideal_adi::CellAdvanceTime) UtoP ";
        cout << "complained (maybe about negative pressure...) fixing\n";
        PtoU(Pf, u2, eq_gamma);
        *dE += (u2[ERG] - u1[ERG]);
        UtoP(u2, Pf, MinTemp, eq_gamma);
    }

    // Reset the dU array for the next timestep.
    for (int v = 0; v < eq_nvar; v++)
        dU[v] = 0.;
    if (MP) {
        MP->sCMA(corrector, Pf);
        for (int t = 0; t < eq_nvar; t++)
            Pf[t] = Pf[t] * corrector[t];
    }
    return 0;
}

// ##################################################################
// ##################################################################

///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_mhd_ideal_adi::CellTimeStep(
    const cell* c,   ///< pointer to cell
    const double,    ///< gas EOS gamma.
    const double dx  ///< Cell size dx.
)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_ideal_adi::CellTimeStep ...starting.\n";
#endif  // FUNCTION_ID

    //
    // Get Max velocity along a grid direction.
    //
    pion_flt u1[eq_nvar];
    pion_flt temp = fabs(c->P[eqVX]);
    if (FV_gndim > 1) temp = max(temp, fabs(c->P[eqVY]));
    if (FV_gndim > 2) temp = max(temp, fabs(c->P[eqVZ]));

    //
    // First rotate the state vector to the fastest directions,
    // then add the fast speed to the max-velocity, and this is
    // the max wavespeed.
    //
    enum axes newdir;
    double cf = 0.0;
    if (FV_gndim == 1)
        temp += cfast(c->P, eq_gamma);
    else {
        // We may have to rotate the state vector to find the fastest
        // fast speed.
        newdir = XX;
        if (fabs(c->P[BY]) < fabs(c->P[BX])) {
            newdir = YY;
            if (fabs(c->P[BZ]) < fabs(c->P[BY])) newdir = ZZ;
        }
        else if (fabs(c->P[BZ]) < fabs(c->P[BX]))
            newdir = ZZ;

        if (newdir != XX) {
            for (int v = 0; v < eq_nvar; v++)
                u1[v] = c->P[v];
            rotate(u1, XX, newdir);
            cf = cfast(u1, eq_gamma);
            temp += cf;
        }
        else {
            cf += cfast(c->P, eq_gamma);
            temp += cf;
        }
    }

    max_speed = max(max_speed, cf);

    FV_dt = dx / temp;
    FV_dt *= FV_cfl;

#ifdef TEST_INF
    if (!isfinite(FV_dt) || FV_dt <= 0.0) {
        cout << "cell has invalid timestep\n";
        CI.print_cell(c);
        cout.flush();
    }
#endif
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_ideal_adi::CellTimeStep ...returning.\n";
#endif  // FUNCTION_ID
    return FV_dt;
}

// ##################################################################
// ##################################################################

// ********************************************************************************
// FV_solver_mhd_mixedGLM_adi class, for the Dedner-GLM divergence cleaning
// method.
// ********************************************************************************

// ##################################################################
// ##################################################################

FV_solver_mhd_mixedGLM_adi::FV_solver_mhd_mixedGLM_adi(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt* state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv, state, gam),
    Riemann_Roe_MHD_CV(nv, gam),
    HLLD_MHD(nv, gam),
    VectorOps_Cart(nd),
    FV_solver_mhd_ideal_adi(nv, nd, cflno, gam, state, avcoeff, ntr),
    eqns_mhd_mixedGLM(nv)
{
    return;
}

// ##################################################################
// ##################################################################

FV_solver_mhd_mixedGLM_adi::~FV_solver_mhd_mixedGLM_adi()
{
    return;
}

// ##################################################################
// ##################################################################

double FV_solver_mhd_mixedGLM_adi::CellTimeStep(
    const cell* c,   ///< pointer to cell
    const double,    ///< gas EOS gamma.
    const double dx  ///< Cell size dx.
)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::CellTimeStep ...starting.\n";
#endif  // FUNCTION_ID
    double dt = FV_solver_mhd_ideal_adi::CellTimeStep(c, eq_gamma, dx);
    return dt;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::inviscid_flux(
    class SimParams& par,       ///< simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid
    const double dx,            ///< cell-size dx (for LF method)
    class cell* Cl,             ///< Left state cell pointer
    class cell* Cr,             ///< Right state cell pointer
    const pion_flt* Pl,         ///< Left Primitive state vector.
    const pion_flt* Pr,         ///< Right Primitive state vector.
    pion_flt* flux,             ///< Resultant Flux state vector.
    pion_flt* pstar,            ///< State vector at interface.
    const int solve_flag,       ///< Solver to use
    const double eq_gamma       ///< Gas constant gamma.
)
{
#ifdef TESTING
    //
    // Check input density and pressure are positive
    //
    if (Pl[eqRO] < TINYVALUE || Pl[eqPG] < TINYVALUE || Pr[eqRO] < TINYVALUE
        || Pr[eqPG] < TINYVALUE) {
        rep.printVec("left ", Pl, eq_nvar);
        rep.printVec("right", Pr, eq_nvar);
        rep.error(
            "FV_solver_mhd_mixedGLM_adi::calculate_flux() Density/Pressure "
            "too small",
            Pl[eqRO]);
    }
#endif  // TESTING

    int err = 0;

    //
    // Set flux and pstar vector to zero.
    //
    for (int v = 0; v < eq_nvar; v++)
        flux[v] = 0.0;
    for (int v = 0; v < eq_nvar; v++)
        pstar[v] = 0.0;
    //
    // Need temporary left and right state vectors b/c we need to change
    // the left and right state values of eqBX to the resolved state of
    // the Dedner et al. (2002) Riemann problem for (Bx,Psi).
    //
    pion_flt left[eq_nvar], right[eq_nvar];
    for (int v = 0; v < eq_nvar; v++)
        left[v] = Pl[v];
    for (int v = 0; v < eq_nvar; v++)
        right[v] = Pr[v];

    ///
    /// First get the resolved state in the 2x2 Dedner Riemann Problem,
    /// which is described well in their paper, and/or my thesis.
    ///
    /// Uses Dedner eq.41 for the flux in Bx and Psi:
    /// \f[ \partial_t B_x + \partial_x \psi = 0 \;, \qquad
    ///     \partial_t \psi + \partial_x (c_h^2 B_x) = 0 \;. \f]
    /// where the source term has been omitted, as it is calculated separately.
    ///
    /// The GLM method has Bx and Psi decoupled from all the other variables
    /// in the Riemann Problem, so they can be solved separately as a two
    /// variable system (Dedner eq.42)
    ///
    /// \f[ F(\psi) = c_h^2 B_x^* = c_h^2 \left( \frac{1}{2}(B_x(L)+B_x(R))
    ///                             - \frac{1}{2c_h}(\psi_R-\psi_L) \right) \f]
    /// \f[ F(B_x) = \psi_* = \frac{1}{2}(\psi_L+\psi_R) -
    /// \frac{c_h}{2}(B_x(R)-B_X(L)) \f]
    ///
    /// Bx(*) is then used for both the left and right state to calculate
    /// the flux.
    ///
    double psistar = 0.0, bxstar = 0.0;
    // set Psi to zero in left and right states, so that Riemann solvers
    // don't get confused (because otherwise it will contribute to the
    // total energy.
    // psistar = 0.5*(left[eqSI]+right[eqSI]);  // Derigs
    // bxstar  = 0.5*(left[eqBX]+right[eqBX]);  // Derigs
    psistar =
        0.5
        * (left[eqSI] + right[eqSI] - (right[eqBX] - left[eqBX]));  // Dedner
    bxstar =
        0.5
        * (left[eqBX] + right[eqBX] - (right[eqSI] - left[eqSI]));  // Dedner

    double psi_L = left[eqSI], psi_R = right[eqSI];
    double bxl = left[eqBX], bxr = right[eqBX];
    left[eqSI] = right[eqSI] = 0.0;
    left[eqBX] = right[eqBX] = bxstar;

    // Now call the ideal-MHD flux solver:
    err = FV_solver_mhd_ideal_adi::inviscid_flux(
        par, grid, dx, Cl, Cr, left, right, flux, pstar, solve_flag, eq_gamma);

    //
    // Now we have to add in the flux in BX and PSI, based on Dedner et
    // al.'s method.
    //
    // NOTE: Dedner doesn't say what to do about the energy flux, so
    // it is modified to ensure consistency (Mackey & Lim, 2011).
    //
    // See Derigs et al. (2018) eq. 4.45: 3 terms f6*, f9*, last term.
    // energy (ERG) is f5, Bx (BBX) is f6, PSI is f9.
    // Derigs et al. (2018) eq. 4.43 f6* and f9*
    // These should be a better solution, but JM could not get them to
    // produce good results.
    //
    // *** N.B. using the Dedner solution here, with Mackey & Lim
    // (2011) correction to the flux, not Dominik's equations, but
    // with Psi defined as in Dominik's paper. ***
    flux[eqERG] += GLM_chyp * bxstar * psistar;
    flux[eqBBX] = GLM_chyp * psistar;
    flux[eqPSI] = GLM_chyp * bxstar;
    left[eqSI]  = psi_L;
    right[eqSI] = psi_R;
    left[eqBX]  = bxl;
    right[eqBX] = bxr;

    return err;
}

// ##################################################################
// ##################################################################

//
// calculate GLM source terms for multi-D MHD and add to Powell source
// Not exactly as indicated in Dominik's paper, but it works.
//
int FV_solver_mhd_mixedGLM_adi::MHDsource(
    class GridBaseClass* grid,  ///< pointer to grid.
    class cell* Cl,             ///< pointer to cell of left state
    class cell* Cr,             ///< pointer to cell of right state
    pion_flt* Pl,               ///< left edge state
    pion_flt* Pr,               ///< right edge state
    const axes d,               ///< Which axis we are looking along.
    enum direction pos,         ///< positive direction normal to interface
    enum direction neg,         ///< negative direction normal to interface
    const double dt             ///< timestep dt
)
{
    double dx = grid->DX();
    double sm = 0.5 * (Cl->Ph[eqSI] + Cr->Ph[eqSI]);
    FV_solver_mhd_ideal_adi::MHDsource(grid, Cl, Cr, Pl, Pr, d, pos, neg, dt);

    pion_flt psi_l[eq_nvar], psi_r[eq_nvar];
    for (int v = 0; v < eq_nvar; v++) {
        psi_l[v] = psi_r[v] = 0.0;
    }
    psi_l[eqERG] = Cl->Ph[eqVX] * Cl->Ph[eqSI];
    psi_l[eqPSI] = Cl->Ph[eqVX];

    psi_r[eqERG] = Cr->Ph[eqVX] * Cr->Ph[eqSI];
    psi_r[eqPSI] = Cr->Ph[eqVX];

    for (int v = 0; v < eq_nvar; v++) {
        Cl->dU[v] -= dt * sm * psi_l[v] / dx;
        Cr->dU[v] += dt * sm * psi_r[v] / dx;
    }
    return 0;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::CellAdvanceTime(
    class cell* c,
    const pion_flt* Pin,  // Initial State Vector.
    pion_flt* dU,         // Update vector dU
    pion_flt* Pf,         // Final state vector (can be same as initial vec.).
    pion_flt* dE,  // Tracks change of energy if I have to correct for negative
                   // pressure
    const double,  // gas EOS gamma.
    const double MinTemp,  ///< Min Temperature allowed on grid.
    const double           // Cell timestep dt.
)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...starting.\n";
#endif  // FUNCTION_ID

    int err = FV_solver_mhd_ideal_adi::CellAdvanceTime(
        c, Pin, dU, Pf, dE, 0, MinTemp, 0);
    GLMsource(&(Pf[eqSI]), FV_dt);

#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...returning.\n";
#endif  // FUNCTION_ID
    return (err);
}

// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::PtoU(
    const pion_flt* p, pion_flt* u, const double g)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::PtoU ...starting.\n";
#endif  // FUNCTION_ID

    eqns_mhd_mixedGLM::PtoU(p, u, g);
    for (int t = 0; t < FV_ntr; t++)
        u[eqTR[t]] = p[eqTR[t]] * p[eqRO];

#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::PtoU ...returning.\n";
#endif  // FUNCTION_ID
    return;
}

// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::UtoP(
    const pion_flt* u,
    pion_flt* p,
    const double MinTemp,  ///< Min Temperature allowed on grid.
    const double g)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::UtoP ...starting.\n";
#endif  // FUNCTION_ID

    for (int t = 0; t < FV_ntr; t++)
        p[eqTR[t]] = u[eqTR[t]] / u[eqRO];
    int err = eqns_mhd_mixedGLM::UtoP(u, p, MinTemp, g);

#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::UtoP ...returning.\n";
#endif  // FUNCTION_ID
    return err;
}

// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::Set_GLM_Speeds(
    const double delt,  ///< timestep dt.
    const double delx,  ///< cell size dx.
    const double cr     ///< GLM damping coefficient c_r
)
{
#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::Set_GLM_Speeds ...starting.\n";
#endif  // FUNCTION_ID

    GLMsetPsiSpeed(FV_cfl * delx / delt, cr);

#ifdef FUNCTION_ID
    cout << "FV_solver_mhd_mixedGLM_adi::Set_GLM_Speeds ...returning.\n";
#endif  // FUNCTION_ID
    return;
}

// ##################################################################
// ##################################################################

/// ---------------------------------------------------------------------
/// -------------------  AXI-SYMMETRIC EQUATIONS ------------------------
/// ---------------------------------------------------------------------

// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_ideal_adi::cyl_FV_solver_mhd_ideal_adi(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt* state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv, state, gam),
    Riemann_Roe_MHD_CV(nv, gam),
    HLLD_MHD(nv, gam),
    VectorOps_Cart(nd),
    FV_solver_mhd_ideal_adi(nv, nd, cflno, gam, state, avcoeff, ntr),
    VectorOps_Cyl(nd)
{
#ifdef FUNCTION_ID
    cout << "::cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif  // FUNCTION_ID

    if (nd != 2) rep.error("Cylindrical coordinates only 2D", nd);

#ifdef FUNCTION_ID
    cout << "::cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif  // FUNCTION_ID
    return;
}

// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_ideal_adi::~cyl_FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
    cout << "::~cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif  // FUNCTION_ID

#ifdef FUNCTION_ID
    cout << "::~cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif  // FUNCTION_ID
}

// ##################################################################
// ##################################################################

void cyl_FV_solver_mhd_ideal_adi::geometric_source(
    cell* c,               ///< Current cell.
    const axes d,          ///< Which axis we are looking along.
    const pion_flt* dpdx,  ///< slope vector for cell c.
    const int OA,          ///< spatial order of accuracy.
    const double dR,       ///< cell length dx.
    pion_flt* dU           ///< add to update vector [OUTPUT]
)
{

    if (d == Rcyl) {
        double pm = (c->Ph[eqBX] * c->Ph[eqBX] + c->Ph[eqBY] * c->Ph[eqBY]
                     + c->Ph[eqBZ] * c->Ph[eqBZ])
                    / 2.;
        switch (OA) {
            case OA1:
                dU[eqMMX] += (c->Ph[eqPG] + pm) / CI.get_dpos(c, Rcyl);
                break;
            case OA2:
                dU[eqMMX] += (c->Ph[eqPG] + pm
                              + (CI.get_dpos(c, Rcyl) - R_com(c, dR))
                                    * (dpdx[eqPG] + c->Ph[eqBX] * dpdx[eqBX]
                                       + c->Ph[eqBY] * dpdx[eqBY]
                                       + c->Ph[eqBZ] * dpdx[eqBZ]))
                             / CI.get_dpos(c, Rcyl);
                break;
            default:
                rep.error(
                    "Bad OOA in cyl_IdealMHD_RS::dU, only know 1st,2nd", OA);
        }
    }

    return;
}

// ##################################################################
// ##################################################################

int cyl_FV_solver_mhd_ideal_adi::MHDsource(
    class GridBaseClass* grid,  ///< pointer to grid.
    class cell* Cl,             ///< pointer to cell of left state
    class cell* Cr,             ///< pointer to cell of right state
    pion_flt* Pl,               ///< left edge state
    pion_flt* Pr,               ///< right edge state
    const axes d,               ///< Which axis we are looking along.
    enum direction pos,         ///< positive normal direction
    enum direction neg,         ///< negative normal direction
    const double dt             ///< timestep dt
)
{
    // The Powell source terms from Powell's paper (1999)
    // called by time_integrator::dynamics_dU_column()
    double dx = grid->DX(), rp = 0.0, rn = 0.0;
    double bm   = 0.5 * (Cl->Ph[eqBX] + Cr->Ph[eqBX]);
    double uB_l = Cl->Ph[eqBX] * Cl->Ph[eqVX] + Cl->Ph[eqBY] * Cl->Ph[eqVY]
                  + Cl->Ph[eqBZ] * Cl->Ph[eqVZ];
    double uB_r = Cr->Ph[eqBX] * Cr->Ph[eqVX] + Cr->Ph[eqBY] * Cr->Ph[eqVY]
                  + Cr->Ph[eqBZ] * Cr->Ph[eqVZ];
    pion_flt Powell_l[eq_nvar], Powell_r[eq_nvar];
    for (int v = 0; v < eq_nvar; v++) {
        Powell_l[v] = Powell_r[v] = 0.0;
    }
    Powell_l[eqMMX] = Cl->Ph[eqBX];
    Powell_l[eqMMY] = Cl->Ph[eqBY];
    Powell_l[eqMMZ] = Cl->Ph[eqBZ];
    Powell_l[eqERG] = uB_l;
    Powell_l[eqBBX] = Cl->Ph[eqVX];
    Powell_l[eqBBY] = Cl->Ph[eqVY];
    Powell_l[eqBBZ] = Cl->Ph[eqVZ];

    Powell_r[eqMMX] = Cr->Ph[eqBX];
    Powell_r[eqMMY] = Cr->Ph[eqBY];
    Powell_r[eqMMZ] = Cr->Ph[eqBZ];
    Powell_r[eqERG] = uB_r;
    Powell_r[eqBBX] = Cr->Ph[eqVX];
    Powell_r[eqBBY] = Cr->Ph[eqVY];
    Powell_r[eqBBZ] = Cr->Ph[eqVZ];

    switch (d) {
        case Zcyl:
            for (int v = 0; v < eq_nvar; v++) {
                Cl->dU[v] -= dt * bm * (Powell_l[v]) / dx;
                Cr->dU[v] += dt * bm * (Powell_r[v]) / dx;
            }
            break;
        case Rcyl:
            rp = CI.get_dpos(Cl, Rcyl) + dx * 0.5;
            rn = rp - dx;
            for (int v = 0; v < eq_nvar; v++) {
                Cl->dU[v] -=
                    dt * bm * (Powell_l[v]) * 2.0 * rp / (rp * rp - rn * rn);
            }
            rn = rp;
            rp += dx;
            for (int v = 0; v < eq_nvar; v++) {
                Cr->dU[v] +=
                    dt * bm * (Powell_r[v]) * 2.0 * rn / (rp * rp - rn * rn);
            }
            break;
        case Tcyl:
            rep.error("3D cylindrical GLM-MHD Source", d);
            break;
        default:
            rep.error("GLM-MHD Source bad direction", d);
            break;
    }

    return 0;
}

// ##################################################################
// ##################################################################

//------------------------------------------------------------------//
cyl_FV_solver_mhd_mixedGLM_adi::cyl_FV_solver_mhd_mixedGLM_adi(
    const int nv,          ///< number of variables in state vector.
    const int nd,          ///< number of space dimensions in grid.
    const double cflno,    ///< CFL number
    const double gam,      ///< gas eos gamma.
    pion_flt* state,       ///< State vector of mean values for simulation.
    const double avcoeff,  ///< Artificial Viscosity Parameter etav.
    const int ntr          ///< Number of tracer variables.
    ) :
    eqns_base(nv),
    FV_solver_base(nv, nd, cflno, gam, avcoeff, ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv, state, gam),
    Riemann_Roe_MHD_CV(nv, gam),
    HLLD_MHD(nv, gam),
    VectorOps_Cart(nd),
    FV_solver_mhd_ideal_adi(nv, nd, cflno, gam, state, avcoeff, ntr),
    eqns_mhd_mixedGLM(nv),
    FV_solver_mhd_mixedGLM_adi(nv, nd, cflno, gam, state, avcoeff, ntr),
    VectorOps_Cyl(nd),
    cyl_FV_solver_mhd_ideal_adi(nv, nd, cflno, gam, state, avcoeff, ntr)
{
#ifdef FUNCTION_ID
    cout << "::cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif  // FUNCTION_ID

    //  cout <<"cyl_FV_solver_mhd_mixedGLM_adi CONSTRUCTOR\n";
    //  cout <<"glmMHD Equations; Riemann Solver Method; Cylindrical
    //  Coordinates.\n";
    if (nd != 2)
        rep.error(
            "Cylindrical coordinates only implemented for \
                        2d axial symmetry so far.  Sort it out!",
            nd);
#ifdef FUNCTION_ID
    cout << "::cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif  // FUNCTION_ID
    return;
}

// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_mixedGLM_adi::~cyl_FV_solver_mhd_mixedGLM_adi()
{
#ifdef FUNCTION_ID
    cout << "::~cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif  // FUNCTION_ID

#ifdef FUNCTION_ID
    cout << "::~cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif  // FUNCTION_ID
}

// ##################################################################
// ##################################################################

void cyl_FV_solver_mhd_mixedGLM_adi::geometric_source(
    cell* c,               ///< Current cell.
    const axes d,          ///< Which axis we are looking along.
    const pion_flt* dpdx,  ///< slope vector for cell c.
    const int OA,          ///< spatial order of accuracy.
    const double dR,       ///< cell length dx.
    pion_flt* dU           ///< update vector to add source term to [OUTPUT]
)
{

    if (d == Rcyl) {
        double pm = (c->Ph[eqBX] * c->Ph[eqBX] + c->Ph[eqBY] * c->Ph[eqBY]
                     + c->Ph[eqBZ] * c->Ph[eqBZ])
                    / 2.;
        switch (OA) {
            case OA1:
                // if (c->pos[Rcyl]<4 && c->pos[Rcyl]>0) {
                //  cout <<"pos="<<c->pos[Rcyl]<<", dU = "<<c->dU[eqMMX] <<"  ";
                //  cout << (c->Ph[eqPG]+pm)/CI.get_dpos(c,Rcyl);
                //}
                dU[eqMMX] += (c->Ph[eqPG] + pm) / CI.get_dpos(c, Rcyl);
                // if (c->pos[Rcyl]<4 && c->pos[Rcyl]>0) {
                //  cout <<"  "<<c->dU[eqMMX]<<"\n";
                //}
                dU[eqBBX] += GLM_chyp * c->Ph[eqSI] / CI.get_dpos(c, Rcyl);
                break;
            case OA2:
                dU[eqMMX] += (c->Ph[eqPG] + pm
                              + (CI.get_dpos(c, Rcyl) - R_com(c, dR))
                                    * (dpdx[eqPG] + c->Ph[eqBX] * dpdx[eqBX]
                                       + c->Ph[eqBY] * dpdx[eqBY]
                                       + c->Ph[eqBZ] * dpdx[eqBZ]))
                             / CI.get_dpos(c, Rcyl);
                dU[eqBBX] +=
                    GLM_chyp
                    * (c->Ph[eqSI]
                       + (CI.get_dpos(c, Rcyl) - R_com(c, dR)) * dpdx[eqSI])
                    / CI.get_dpos(c, Rcyl);
                break;
            default:
                rep.error("Bad OOA in cyl_glmMHD::dU, only know 1st,2nd", OA);
        }
    }

    return;
}

// ##################################################################
// ##################################################################

///
/// calculate GLM source terms for multi-D MHD and add to Powell source
/// Not exactly as indicated in Dominik's paper, but it works.
///
int cyl_FV_solver_mhd_mixedGLM_adi::MHDsource(
    class GridBaseClass* grid,  ///< pointer to grid.
    class cell* Cl,             ///< pointer to cell of left state
    class cell* Cr,             ///< pointer to cell of right state
    pion_flt* Pl,               ///< left edge state
    pion_flt* Pr,               ///< right edge state
    const axes d,               ///< Which axis we are looking along.
    enum direction pos,         ///< positive direction normal to interface
    enum direction neg,         ///< negative direction normal to interface
    const double dt             ///< timestep dt
)
{
    double dx = grid->DX();
    double sm = 0.5 * (Cl->Ph[eqSI] + Cr->Ph[eqSI]);
    cyl_FV_solver_mhd_ideal_adi::MHDsource(
        grid, Cl, Cr, Pl, Pr, d, pos, neg, dt);

    pion_flt psi_l[eq_nvar], psi_r[eq_nvar];
    for (int v = 0; v < eq_nvar; v++) {
        psi_l[v] = psi_r[v] = 0.0;
    }
    psi_l[eqERG] = Cl->Ph[eqVX] * Cl->Ph[eqSI];
    psi_l[eqPSI] = Cl->Ph[eqVX];

    psi_r[eqERG] = Cr->Ph[eqVX] * Cr->Ph[eqSI];
    psi_r[eqPSI] = Cr->Ph[eqVX];

    for (int v = 0; v < eq_nvar; v++) {
        Cl->dU[v] -= dt * sm * psi_l[v] / dx;
        Cr->dU[v] += dt * sm * psi_r[v] / dx;
    }
    return 0;
}

// ##################################################################
// ##################################################################
