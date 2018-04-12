///
/// \file flux_mhd_adiabatic.h
/// \author Jonathan Mackey
/// Declaration of the adiabatic MHD flux solver class
///
/// History:
/// - 2009-10-20 Started on the file
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.11.15 JM: Renamed calculate_flux() to inviscid_flux() and
///   moved AV calculation to FV_solver_base class.
/// - 2010.12.27 JM: Moved Roe flux solver to own class in Riemann_solvers/
///   Got rid of inherited class flux/left/right/pstar variables.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#ifndef FLUX_MHD_ADIABATIC_H
#define FLUX_MHD_ADIABATIC_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "flux_base.h"

///
/// Flux Solver Class for Adiabatic MHD
///
/// This uses the adiabatic equations, and calculates the flux across an interface.
///
/// \author Jonathan Mackey
/// Written 2009-10-20
///

class flux_solver_mhd_ideal_adi :
  virtual public flux_solver_base,
{
public:
  ///
  /// Constructor: Allocates memory for arrays, and sets a state
  /// vector of typical values for primitive variables (can be
  /// overwritten later).
  ///
  flux_solver_mhd_ideal_adi(
        const int,      ///< Number of variables in state vector
        const pion_flt *, ///< state vector which is 'typical' in the problem being solved. 
        const double,   ///< coefficient of artificial viscosity
        const double,    ///< gamma (EOS).
        const int     ///< Number of tracer variables.
        );

  /// Deletes any arrays.
  ~flux_solver_mhd_ideal_adi();
  

};



class flux_solver_mhd_mixedGLM_adi
: virtual public flux_solver_mhd_ideal_adi,
{
  public:
  flux_solver_mhd_mixedGLM_adi(
      const int,      ///< Number of variables in state vector
      const pion_flt *, ///< state vector which is 'typical' in the problem being solved. 
      const double,   ///< coefficient of artificial viscosity
      const double,    ///< gamma (EOS).
      const int     ///< Number of tracer variables.
      );

  ~flux_solver_mhd_mixedGLM_adi();
   

#endif //FLUX_MHD_ADIABATIC_H

