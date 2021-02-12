///
/// \file Riemann_FVS_hydro.h
/// \author Jonathan Mackey
///
/// Flux vector splitting for the Euler Equations.  This is a fairly
/// self-contained class, depending only on the Euler Equations class.
/// It doesn't use the class variables from the other Riemann solvers.
///
/// References:
///  - Van Leer, B., 1982, LNP, 170, 507.
///    "Flux vector splitting for the Euler equations."
///  - Van Leer, B., 1991, "Flux Vector Splitting for the 1990s".
///  - http://www.chimeracfd.com/programming/gryphon/fluxvanleer.html
///    This ref. has a typo in eq.24 -- the f_b^{pm} term should only have one
///    +/-
///
/// History:
///  - 2010-09-21 JM: Started writing.
/// - 2010.12.22 JM: Moved to sub-directory containing Riemann solvers.
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer
///    fluxes are dealt with by the flux-solver classes.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#ifndef RIEMANN_FVS_HYDRO_H
#define RIEMANN_FVS_HYDRO_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_hydro_adiabatic.h"

class Riemann_FVS_Euler : virtual public eqns_Euler {
  public:
    Riemann_FVS_Euler(
        const int,    ///< Length of State Vectors, nvar
        const double  ///< Gamma for state vector.
    );

    ///
    /// Destructor: deletes dynamically allocated data.
    ///
    ~Riemann_FVS_Euler();

    ///
    /// Gets the solution to the Riemann Problem.
    ///
    int FVS_flux(
        const pion_flt*,  ///< Left Primitive var. state vector.
        const pion_flt*,  ///< Right Primitive var. state vector.
        pion_flt*,        ///< Result Flux vector.
        pion_flt*,        ///< Interface state (for viscosity)
        const double      ///< Gas constant gamma (unused b/c g=constant)
    );

  private:
    const int
        rs_nvar;  ///< length of state vectors in solver (ignore tracers!).

    pion_flt *fpos, *fneg;  ///< split flux vectors.

    void Roe_average_state(
        const pion_flt*,  ///< state 1
        const pion_flt*,  ///< state 2
        const double,     ///< gamma
        pion_flt*         ///< Roe-averaged state
    );
};

#endif  // RIEMANN_FVS_HYDRO_H
