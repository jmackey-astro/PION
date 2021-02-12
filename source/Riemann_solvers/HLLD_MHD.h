
#ifndef HLLD_MHD_H
#define HLLD_MHD_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_mhd_adiabatic.h"

class HLLD_MHD : virtual public eqns_mhd_ideal {
  public:
    ///
    /// Constructor: passes (eq_nvar, gamma) to equations setup.
    ///
    HLLD_MHD(
        const int,    ///< Length of State Vectors, nvar
        const double  ///< Gamma for state vector.
    );

    ///
    /// Destructor: Deletes dynamically allocated private data.
    ///
    ~HLLD_MHD();

  protected:
    ///
    /// The HLLD MHD solver (reference!)
    ///
    int MHD_HLLD_flux_solver(
        const pion_flt*,  ///< input left state
        const pion_flt*,  ///< input right state
        const double,     ///< input eq_gamma
        pion_flt*,        ///< output flux
        pion_flt*         ///< output Ustar vec
    );

    double *HD_lambda, *HD_UL, *HD_UR, *HD_FL, *HD_FR, *HD_ULs, *HD_URs,
        *HD_FLs, *HD_FRs, *HD_ULss, *HD_URss, *HD_FLss, *HD_FRss;
    double HD_nvar;

    void HLLD_signal_speeds(
        const pion_flt*,  ///< inputs
        const pion_flt*,
        const double,
        double&,  ///< outputs
        double&);

    int MHD_HLL_flux_solver(
        const pion_flt*,  ///< input left state
        const pion_flt*,  ///< input right state
        const double,     ///< input gamma
        pion_flt*,        ///< output flux
        pion_flt*         ///< output Ustar vec
    );
};

#endif  // HLLD_MHD_H
