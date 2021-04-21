///
/// \file Roe_Hydro_ConservedVar_solver.h
/// \author Jonathan Mackey
///
/// This is a linearised Riemann solver for the Euler Equations.
/// It is the Roe solver in conserved variables, and it can use
/// the H-correction to add multi-dimensional viscosity to the
/// fluxes.  The symmetric solver is recommended, the one-sided is
/// not.  Both will return a flux in the first 5 state variables (i.e.
/// the hydro vars) and also a state vector for the interface which
/// can be used for adding viscosity.
///
/// References:
///  - Toro (1999) Riemann Solvers Textbook, Chapter 11.2.2, pp.350-353.
///  - Sanders, Morano, Druguet, (1998, JCP, 145, 511).
///
/// History:
/// - 2010-12.22 JM: Moved functions from flux_hydro_adiabatic.h/.cc
///   and generated new class.  Split everything into functions.
///   Defined a bunch of private class variables.
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer
///    fluxes are dealt with by the flux-solver classes.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2016.05.21 JM: removed H-correction ifdefs (it should be always enabled).

#ifndef ROE_HYDRO_CONSERVEDVAR_SOLVER_H
#define ROE_HYDRO_CONSERVEDVAR_SOLVER_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_hydro_adiabatic.h"

class Riemann_Roe_Hydro_CV : virtual public eqns_Euler {
public:
  ///
  /// Constructor: doesn't do much.
  ///
  Riemann_Roe_Hydro_CV(
      const int,    ///< Length of State Vectors, nvar
      const double  ///< Gamma for state vector.
  );

  ///
  /// Destructor: deletes dynamically allocated data.
  ///
  ~Riemann_Roe_Hydro_CV();

  ///
  /// Symmetric version of Roe's Flux solver for the Euler Equations,
  /// from Toro (1999), chapter 11.2.2, pp.350-353.
  ///
  int Roe_flux_solver_symmetric(
      const pion_flt *,  ///< input left state
      const pion_flt *,  ///< input right state
      const double,      ///< input gamma
      const double,      ///< H-correction eta-max value.
      pion_flt *,        ///< output pstar
      pion_flt *         ///< output flux
  );

private:
  /// H-correction eta-max value (pre-calculated for this interface)
  pion_flt RCV_HC_etamax;

  const int rs_nvar;  ///< length of state vectors in solver (ignore tracers).

  /// Square of the gas velocity in the mean state.
  double RCV_v2_mean;

  /// Adiabatic sound speed in mean state vector.
  double RCV_a_mean;

  /// Index for the enthalpy variable in state vectors.
  enum primitive eqHH;

  /// Local copy of mean state vector (5 elements only!)
  pion_flt *RCV_meanp;

  /// Eigenvalues
  pion_flt *RCV_eval;

  /// Wave strengths.
  pion_flt *RCV_strength;

  /// Difference vector
  pion_flt *RCV_udiff;

  /// Left and right state conserved variable state.
  pion_flt *RCV_ul, *RCV_ur;

  /// matrix of eigenvectors
  pion_flt **RCV_evec;

  ///
  /// Returns true if the left and right states are very similar.
  ///
  int test_left_right_equality(
      const pion_flt *,  ///< input left state
      const pion_flt *   ///< input right state
  );

  ///
  /// Set the Roe-average mean state vector RCV_meanp[].
  /// Also set the variables RCV_v2_mean and RCV_a_mean, which are
  /// the square of the velocity and the sound speed (adiabatic) in
  /// the mean state.
  ///
  void set_Roe_mean_state(
      const pion_flt *,  ///< input left state
      const pion_flt *   ///< input right state
  );

  ///
  /// Set left and right states in conservative variables, and then get
  /// the difference vector.  Store results in ul[], ur[], udiff[].
  ///
  void set_ul_ur_udiff(
      const pion_flt *,  ///< input left state
      const pion_flt *   ///< input right state
  );

  ///
  /// Assigns eigenvalues to eval[], based on meanp,a_mean.
  /// If using the H-correction, this alters the evalues here.
  ///
  void set_eigenvalues();

  ///
  /// Assigns eigenvectors to evec[][], based on meanp,a_mean,v2_mean.
  ///
  void set_eigenvectors();

  ///
  /// this uses udiff, RCV_meanp, eq_gamma, RCV_a_mean to
  /// calculate the wavestrengths according to Toro (1999, eq.11.68/.69/.70).
  /// Strengths are set in the class variable str[].
  ///
  void set_wave_strengths();

  ///
  /// Calculate the symmetric Roe flux (for the 5 physical variables
  /// only!) by stepping across all waves.
  ///
  void calculate_symmetric_flux(pion_flt *  ///< output flux vector
  );

  ///
  /// Calculate asymmetric (one-sided) Roe flux by stepping across
  /// waves from the left state until we get to the starred state.
  /// Again this is for the 5 hydrodynamical variables only.  This
  /// function is really superceded by the symmetric version since
  /// this version does not maintain symmetry in symmetric problems!
  ///
  void calculate_asymmetric_flux(
      const pion_flt *, const pion_flt *, pion_flt *);

  ///
  /// Set Pstar[] from RCV_meanp[] (need to replace enthalpy with
  /// pressure).
  ///
  void set_pstar_from_meanp(pion_flt *  ///< output pstar vector
  );
};

#endif  // ROE_HYDRO_CONSERVEDVAR_SOLVER_H
