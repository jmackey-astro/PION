///
/// \file Roe_Hydro_PrimitiveVar_solver.h
/// \author Jonathan Mackey
///
/// Takes in a left and right state in primitive variables, and
/// returns the intermediate state Pstar from which the flux across
/// the interface can be calculated.  This uses the Roe average for
/// the mean state, and is a linear Riemann solver.
///
/// References:
///  - Toro (1999) Riemann Solvers Textbook, Chapter 11.
///
/// History:
/// - 2010-12.22 JM: Moved functions from flux_hydro_adiabatic.h/.cc
///   and generated new class.
///
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer 
///    fluxes are dealt with by the flux-solver classes.
///

#ifndef ROE_HYDRO_PRIMITIVEVAR_SOLVER_H
#define ROE_HYDRO_PRIMITIVEVAR_SOLVER_H

#include "../equations/eqns_hydro_adiabatic.h"

class Riemann_Roe_Hydro_PV : virtual public eqns_Euler 
{
public:
  ///
  /// Constructor: doesn't do much.
  ///
  Riemann_Roe_Hydro_PV(const int,      ///< Length of State Vectors, nvar
		       const double    ///< Gamma for state vector.
		       );
  ///
  /// Destructor: deletes dynamically allocated data.
  ///
  ~Riemann_Roe_Hydro_PV();

  ///
  /// Primitive variable linearised solver for the Euler Equations,
  /// using the Roe avergae, see Toro (1999).
  ///
  int Roe_prim_var_solver(const double *, ///< input left state
			  const double *, ///< input right state
			  const double,   ///< input EOS gamma
			  double *        ///< output pstar.
			  );

private:
  const int rs_nvar; ///< length of state vectors in solver (ignore tracers!).
};



#endif // ROE_HYDRO_PRIMITIVEVAR_SOLVER_H
