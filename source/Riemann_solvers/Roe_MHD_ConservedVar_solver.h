///
/// \file Roe_MHD_ConservedVar_solver.h
/// \author Jonathan Mackey
///
/// Linearised Riemann solver for the ideal MHD equations.
/// This is a Roe solver in conserved variables.  It has a one-sided
/// solver and a symmetric solver which maintains symmetry
/// in a problem much more successfully.
/// The symmetric solver is recommended, the one-sided solver is
/// not.  Both will return a flux in the first 8 state variables (i.e.
/// the physical vars) and also a state vector for the interface which
/// can be used for adding viscosity.
/// 
/// References:\n
/// - Cargo & Gallice (1997) JCP, 136, 446
/// - Stone et al. (2009), ApJS, 178, 137
///
/// History:
/// - 2010.12.27 JM: Moved from flux_mhd_adiabatic.h
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
///

#ifndef ROE_MHD_CONSERVEDVAR_SOLVER_H
#define ROE_MHD_CONSERVEDVAR_SOLVER_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "equations/eqns_mhd_adiabatic.h"

//
// Flag for testing the Roe solver:
//
//#define RoeMHD_TESTING

//
// Set this if we want to use the interface state rather than the
// mean state for post-processing (adding diffusion etc.).  This is 
// probably not a good idea since Ustar does not have
// positive--definite pressure.
//
//#define MHD_ROE_USE_USTAR


///
/// Linearised Riemann solver for the ideal MHD equations.
/// This is a Roe solver in conserved variables.  It has a one-sided
/// solver and a symmetric solver which maintains symmetry
/// in a problem much more successfully.
/// The symmetric solver is recommended, the one-sided solver is
/// not.  Both will return a flux in the first 8 state variables (i.e.
/// the physical vars) and also a state vector for the interface which
/// can be used for adding viscosity.
/// 
/// References:\n
/// - Cargo & Gallice (1997) JCP, 136, 446
/// - Stone et al. (2009), ApJS, 178, 137
///
class Riemann_Roe_MHD_CV :  virtual public eqns_mhd_ideal
{
 public:
  ///
  /// Constructor: passes (eq_nvar, gamma) to equations setup.
  /// Allocated private data arrays.
  ///
  Riemann_Roe_MHD_CV(
      const int,      ///< Length of State Vectors, nvar
      const double    ///< Gamma for state vector.
      );

  ///
  /// Destructor: Deletes dynamically allocated private data.
  ///
  ~Riemann_Roe_MHD_CV();


 protected:
  ///
  /// Roe's approximate flux solver: input left and right states with
  /// gamma (and the H-correction eta-max value if needed), and the
  /// function returns the flux across the interface, and the
  /// Roe-average state (or the interface state if MHD_ROE_USE_USTAR
  /// is set -- it isn't usually because Ustar can have negative 
  /// pressure whereas meanP is positive--definite).
  /// This is the one-sided flux calculation (only waves with v<0).
  ///
  int MHD_Roe_CV_flux_solver_onesided(
      const pion_flt *, ///< input left state
      const pion_flt *, ///< input right state
      const double,   ///< input gamma
#ifdef HCORR
      const pion_flt, ///< H-correction eta-max value.
#endif // HCORR
      pion_flt *, ///< output pstar
      pion_flt *  ///< output flux
      );
  
  ///
  /// Roe's approximate flux solver: input left and right states with
  /// gamma (and the H-correction eta-max value if needed), and the
  /// function returns the flux across the interface, and the
  /// Roe-average state as Pstar[].
  /// This is the symmetric flux calculation (all waves summed over).
  ///
  int MHD_Roe_CV_flux_solver_symmetric(
      const pion_flt *, ///< input left state
      const pion_flt *, ///< input right state
      const double,   ///< input gamma
#ifdef HCORR
      const pion_flt, ///< H-correction eta-max value.
#endif // HCORR
      pion_flt *, ///< output pstar
      pion_flt *  ///< output flux
      );
  
 private:
  ///
  /// Set UL[] and UR[] from PL[] and PR[].
  ///
  void set_UL_and_UR(
      const pion_flt *, ///< left primitive vec.
      const pion_flt *  ///< right primitive vec.
      );

  ///
  /// Set Pstar[] from Roe_meanp[] (need to replace enthalpy with
  /// pressure).
  ///
  void set_pstar_from_meanp(
      pion_flt * ///< output pstar vector
      );

  ///
  /// Get the Roe averages for the primitive variables:
  /// From Stone et al. (2009), ApJS, 178, 137, eq.65.
  ///
  int Roe_get_average_state(
      const pion_flt *, ///< left primitive vec.
      const pion_flt *  ///< right primitive vec.
      );

  ///
  /// Get the Roe averages differences for the primitive and conserved variables:
  /// From Stone et al. (2009), ApJS, 178, 137, eq.65.
  ///
  int Roe_get_difference_states(
      const pion_flt *, ///< left primitive vec.
      const pion_flt *  ///< right primitive vec.
      );
  ///
  /// Get the Roe-averaged wavespeeds
  ///
  int Roe_get_wavespeeds();

  ///
  /// Get the Roe-averaged eigenvalues, modified by the H-correction
  /// if required.  Note that eta-max is set to zero if not using the
  /// H-correction.
  ///
  int Roe_get_eigenvalues(
#ifdef HCORR
	const pion_flt ///< H-correction eta-max.
#endif // HCORR
);

  ///
  /// Get the Roe-averaged wave strengths,
  /// from Cargo & Gallice (1997) JCP, 136, 446, eq.4.20
  ///
  int Roe_get_wavestrengths();

  ///
  /// Calculate the Right Eigenvectors of the average state,
  /// from Cargo & Gallice (1997) JCP, 136, 446, eq.4.18,4.19
  ///
  int Roe_get_right_evectors();

  ///
  /// Using the evalues,wave-strengths,evectors, calculate the
  /// Roe-average Flux from the left state across to zero.
  ///
  int Roe_get_flux_onesided(
      const pion_flt *, ///< left primitive vec.
      const pion_flt *,  ///< right primitive vec.
#ifdef MHD_ROE_USE_USTAR
      pion_flt *, ///< output Pstar vector
#endif // MHD_ROE_USE_USTAR
      pion_flt * ///< output flux vector.
      );

  ///
  /// Using the evalues,wave-strengths,evectors, calculate the
  /// Roe-average Flux by summing across all waves (this is the 
  /// symmetric summation.
  ///
  int calculate_symmetric_flux(
      const pion_flt *, ///< left primitive vec.
      const pion_flt *, ///< right primitive vec.
      pion_flt *        ///< output flux vector.
      );


  pion_flt *Roe_evalues;   ///< eigenvalues vector.
  pion_flt *Roe_strengths; ///< wave strengths.
  pion_flt *Roe_meanp;   ///< Mean state vector (prim.var.).
  pion_flt *Roe_UL;  ///< left state (conserved var.).
  pion_flt *Roe_UR;  ///< right state (conserved var.).
  pion_flt *Roe_udiff; ///< conserved variable differences.
  pion_flt *Roe_pdiff; ///< primitive variable differences.
  pion_flt **Roe_right_evecs; ///< the seven right eigenvectors.
  double 
     Roe_V,  ///< mean state velocity magnitude
     Roe_B,  ///< mean state B-field magnitude
     Roe_a,  ///< sound speed
     Roe_cf, ///< fast speed
     Roe_cs, ///< slow speed
     Roe_ca; ///< Alfven speed.
  double
     Roe_alphas, ///< Roe-Balsara slow norm. 
     Roe_alphaf, ///< Roe-Balsara fast norm.
     Roe_betay,  ///< By/Bt
     Roe_betaz,  ///< Bz/Bt
     Roe_Bt,     ///< sqrt(By^2+Bz^2)
     Roe_signBX, ///< sign of Bx.
     Roe_denom,  ///<  = 1/(sqrt(rho_l)+sqrt(rho_r))
     Roe_CGparamX; ///< CG97's parameter X=deltaB^2/(2(rrl+rrr)^2)
  enum primitive eqHH; ///< location of enthalpy in primitive var vector.
  double Enthalpy(
      const pion_flt *, ///< State Vector.
      const double ///< gas EOS gamma.
      );

};

#endif // ROE_MHD_CONSERVEDVAR_SOLVER_H
