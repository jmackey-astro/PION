



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include "Riemann_solvers/HLLD_MHD.h"
#include "Riemann_solvers/riemannMHD.h"

using namespace std;



// ##################################################################
// ##################################################################



HLLD_MHD::HLLD_MHD(
      const int nv,      ///< Length of State Vectors, nvar
      const double g    ///< Gamma for state vector.
      )
 :  eqns_base(nv), eqns_mhd_ideal(nv)
{
  return;
}



// ##################################################################
// ##################################################################


HLLD_MHD::~HLLD_MHD()
{
  return;
}



// ##################################################################
// ##################################################################


int HLLD_MHD::MHD_HLLD_flux_solver(
      const pion_flt *left,  ///< input left state
      const pion_flt *right, ///< input right state
      const double gamma,    ///< input gamma
      const pion_flt etamax, ///< H-correction eta-max value.
      pion_flt *pstar,       ///< output pstar
      pion_flt *flux         ///< output flux
      )
{
  // calculate pstar and flux...

  return 0;
}



