
#ifndef HLLD_HYDRO_H
#define HLLD_HYDRO_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "equations/eqns_hydro_adiabatic.h"

class HLL_hydro : virtual public eqns_Euler {
public:
  ///
  /// Constructor: passes (eq_nvar, gamma) to equations setup.
  ///
  HLL_hydro(
      const int,    ///< Length of State Vectors, nvar
      const double  ///< Gamma for state vector.
  );

  ///
  /// Destructor: Deletes dynamically allocated private data.
  ///
  ~HLL_hydro();

protected:
  double *HD_lambda, *HD_UL, *HD_UR, *HD_FL, *HD_FR;

  /// Calculate the positive and negative max. signal speeds based
  /// on the left and right states.
  void HLL_signal_speeds(
      const pion_flt*,  ///< inputs
      const pion_flt*,
      const double,
      double&,  ///< outputs
      double&);

  /// Calculate the HLL flux and interface state.
  int hydro_HLL_flux_solver(
      const pion_flt*,  ///< input left state
      const pion_flt*,  ///< input right state
      const double,     ///< input gamma
      pion_flt*,        ///< output flux
      pion_flt*         ///< output interface state (cons.var.)
  );
};

#endif  // HLL_HYDRO_H
