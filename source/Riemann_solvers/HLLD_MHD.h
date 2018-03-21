
#ifndef HLLD_MHD_H
#define HLLD_MHD_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_mhd_adiabatic.h"

class HLLD_MHD :  virtual public eqns_mhd_ideal
{
 public:
  ///
  /// Constructor: passes (eq_nvar, gamma) to equations setup.
  ///
  HLLD_MHD(
      const int,      ///< Length of State Vectors, nvar
      const double    ///< Gamma for state vector.
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
      const pion_flt *, ///< input left state
      const pion_flt *, ///< input right state
      const double,     ///< input gamma
      const pion_flt,   ///< H-correction eta-max value.
      pion_flt *,       ///< output pstar
      pion_flt *        ///< output flux
      );

};




#endif // HLLD_MHD_H
