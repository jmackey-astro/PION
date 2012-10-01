/// \file flux_base.cc
/// \author Jonathan Mackey
/// 
/// Contains class declarations for the basic flux-solver class.
/// This is an interface class, and also defines some basic
/// functionality, such as how to deal with tracer variables tacked
/// onto the end of the state vectors.
///
/// Modifications:
///
/// - 2010.12.23 JM: Moved from eqns_base.h into its own file.
///   Added tracer_flux() function.
///
/// - 2010.12.27 JM: Enclosed Lapidus AV in an ifdef.
///
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.


#include "../global.h"
#include "flux_base.h"
using namespace std;

// **********************************
// ***** BASE FLUX SOLVER CLASS *****
// **********************************

flux_solver_base::flux_solver_base(const int nv,    ///< length of state vector.
				   const double eta, ///< coefficient of artificial viscosity.
				   const int ntr     ///< Number of tracer variables.
				   )
  : eqns_base(nv),
    //    riemann_base(nv),
    FS_etav(eta), FS_etaB(eta), FS_ntr(ntr)
{

  //
  // Allocate memory for tracer indices, and set their values:
  //
  eqTR = 0;
  if (FS_ntr>0) {
    eqTR = mem.myalloc(eqTR, FS_ntr);
  //  cout <<"\tSetting tracer variables to be last "<<ntr<<" elements of state vector.\n";
    for (int i=0;i<FS_ntr;i++) eqTR[i] = eq_nvar-FS_ntr+i;
  }

  HC_etamax=0.0;
  return;
}

flux_solver_base::~flux_solver_base()
{
  //
  // Delete the array of tracer indices
  //
  if (FS_ntr>0) {
    eqTR = mem.myfree(eqTR);
  }
  return;
}

int flux_solver_base::get_LaxFriedrichs_flux(const double *l,
					     const double *r,
					     double *f,
					     const double gamma
					     )
{
  //
  // Need some temporary arrays:
  //
  rep.error("Don't call get_LaxFriedrichs_flux() without an implementation!!!",99);
  return 0;
}

#ifdef LAPIDUS_VISCOSITY_ENABLED
int flux_solver_base::AVLapidus(const cell *Cl, ///< Left state cell pointer
				const cell *Cr, ///< Right state cell pointer
				double *flux, 
				const double etav,
				const double gamma
				)
{
  /// This is not working! Don't use it!
  rep.error("Fix Lapidus Viscosity!",99);
  return 0;
}
#endif // LAPIDUS_VISCOSITY_ENABLED

void flux_solver_base::set_interface_tracer_flux(const double *left, //prim.var.
						 const double *right,//prim.var.
						 double *flux
						 )
{
#ifdef FUNCTION_ID
  cout <<"flux_solver_base::set_interface_tracer_flux ...starting.\n";
#endif //FUNCTION_ID
  //
  // Calculate tracer flux here -- if mass flux is positive then
  // contact is at x>0 and we advect the left state tracer across
  // the boundary.  Otherwise we advect the right state to the left.
  // 
  if (FS_ntr>0) {
    if (flux[eqRHO]>=0.0)
      for (int t=0;t<FS_ntr;t++)
	flux[eqTR[t]] =  left[eqTR[t]]*flux[eqRHO];
    else 
      for (int t=0;t<FS_ntr;t++)
	flux[eqTR[t]] = right[eqTR[t]]*flux[eqRHO];
  }

#ifdef FUNCTION_ID
  cout <<"flux_solver_base::set_interface_tracer_flux ...returning.\n";
#endif //FUNCTION_ID
  return;
}
