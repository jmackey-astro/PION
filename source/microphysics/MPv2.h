///
/// \file MPv2.h
/// \author Jonathan Mackey
/// \date 2011.03.15
///
/// modifications:
/// - getting it written: mods up until 2011.03.XX
/// - 2011.03.21 JM: Updated  RTnew() interface for more sources.  It is now simpler.
/// - 2011.03.29 JM: Added intermediate class mp_rates_ExpH_ImpMetals (Explicit-Hydrogen,
///    Implicit-Metal treatments).  This is so it will work with CVodes.
/// - 2011.03.31 JM: Finished coding, fixed a lot of bugs, need to test it now.
/// - 2011.04.14 JM: In process of testing.  fixed a few things.
///   2011.04.17 JM: Debugging.
/// - 2011.05.02 JM: Added set_multifreq_source_properties() function
/// - 2011.05.04 JM: Added bounding values for ion fraction [eps,1-eps].
/// - 2018.03.20 JM: renamed class

#ifndef MPV2_H
#define MPV2_H

/// Description:
/// This class is an update on the microphysics class used for JM's thesis.  
/// THIS IS LEGACY CODE THAT IS NO LONGER USED.
///
/// - It uses multi-frequency photoionisation including spectral hardening with
///   optical depth, using the method outlined in Frank & Mellema
///   (1994,A&A,289,937), and updated by Mellema et al. (2006,NewA,11,374).
///
/// - It uses the Hummer (1994,MN,268,109) rates for radiative recombination and
///   its associated cooling, and Bremsstrahlung cooling.
///
/// - For collisional ionisation the function of Voronov (1997,ADNDT,65,1) is used.
///
/// - Collisional excitation of neutral H, table from Raga, Mellema, & Lundqvist
///   (1997,ApJS,109,517), using data from Aggarwal (1993).
///
/// - For cooling due to heavy elements, which are not explicitly included, we use
/// the CIE cooling curve of Sutherland & Dopita (1993,ApJS,88,253), but this time I
/// subtract the zero-metals curve from the solar-metallicity curve so that I only
/// have cooling due to metals.  This eliminates the potential for double-counting
/// which was in my previous cooling function.
///
/// - Then I use various formulae from Henney et al. (2009,MN,398,157) for cooling due
/// to collisional excitation of photoionised O,C,N (eq.A9), collisional excitation of
/// neutral metals (eq.A10), and the SD93 CIE curve.  I think I will take the max of
/// the SD93 curve and the Henney et al. functions, to avoid double counting.  For 
/// neutral gas I can use the cooling curve of Henney et al 2009 eq.A14.
///
/// - Photoheating from ionisation is discussed above.  Cosmic ray heating will use a
/// standard value, X-ray heating I'm not yet sure about.  UV heating due to the 
/// interstellar radiation field (ISRF) will be according to the optical depth from the
/// edge of the domain to the cell in question, using e.g. HEA09 eq.A3.  UV heating
/// from the star can use the same equation, but with the optical depth from the 
/// source (using a total nucleon column density).
///


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifdef LEGACY_CODE

#include <vector>
#include "microphysics_base.h"
#include "cooling_SD93_cie.h"
#include "hydrogen_mp.h"
#include "hydrogen_recomb_Hummer94.h"
#include "hydrogen_photoion.h"


// ##################################################################
// ##################################################################



///
/// This class cannot be part of MPv2 because it cannot be the MP global pointer.
/// It needs a local initialisation (with different interface functions) which CVodes
/// can use.
///
class mp_rates_ExpH_ImpMetals
  :
  public Hydrogen_chem,
  public cooling_function_SD93CIE
{
  public:
  mp_rates_ExpH_ImpMetals();
  ~mp_rates_ExpH_ImpMetals();

  ///
  /// This should be called by MPv2 constructor.  We can't set it up
  /// in this class's constructor b/c the global gamma won't be set at that
  /// stage.  If we have an ionising source then this functions gets the 
  /// source luminosity from SimPM.RS.sources (and Rstar,Teff if multifreq).
  ///
  void set_gamma_and_srcs(
      const double, ///< gamma value
      const int,    ///< No diffuse  sources --> 0, otherwise --> 1
      const int     ///< No ionising sources --> 0, otherwise --> 1
      );

  ///
  /// This returns the rate of change of the Hydrogen ion fraction x, and the
  /// internal energy density Eint.
  ///
  int dYdt(
      const double, ///< (1-x_in): input neutral fraction
      const double, ///< E_in: input internal energy density
      double *, ///< xdot returned rate
      double *  ///< Edot returned rate
      );

  ///
  /// This should be called before the microphysics integration, to set the density,
  /// column densities, Vshell, etc.
  ///
  void set_parameters_for_current_step(
      const std::vector<double> & ///< list of parameters in an array.
      );

  ///
  /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
  /// appropriate for a gas with 10% Helium by number, and if He is singly ionised
  /// whenever H is.
  ///
  double get_temperature(
      const double, ///< nH
      const double, ///< E_int
      const double  ///< x(H+)
      );
        
  ///
  /// Get the number of extra parameters and the number of equations.
  /// Only needed by the CVodes Ydot function.
  ///
  void get_problem_size(
      int *, ///< number of equations
      int *  ///< number of parameters in user_data vector.
      );

  ///
  /// returns the relative error tolerance required for this system of equations,
  /// and absolute error tolerances for both x(H+) and E_int
  ///
  void get_error_tolerances(
      double &, ///< relative error tolerance.
      std::vector<double> & ///< absolute error tolerances
      );
  
  //
  // Public functions inherited from Hummer94_Hrecomb:
  //  Hii_rad_recomb_rate(T)
  //  Hii_rad_recomb_cooling(T)
  //  Hii_total_cooling(T) (case B recomb +Bremsstrahlung)
  //
  // Public functions inherited from hydrogen_photoion:
  //  void Setup_photoionisation_rate_table(Tstar,Rstar,etc) (cgs)
  //  double Hi_multifreq_photoionisation_rate(NH,nH,Vshell) (cgs)
  //  double Hi_multifreq_photoionisation_heating_rate(NH,nH,Vshell) (cgs)
  //  double Hi_monochromatic_photo_ion_xsection(E) (cgs)
  //
  // Public functions inherited from Hydrogen_chem:
  //  double Hi_coll_ion_rate(T);
  //  double Hi_coll_ion_cooling_rate(T);
  //  void Hi_coll_ion_rates(T,&CIR,&CICR);
  //  double Hi_coll_excitation_cooling_rate(T);
  //  double Hi_monochromatic_photo_ion_heating(E_phot);
  //
  // Public functions from SD93 cooling class:
  //  void setup_SD93_cie();
  //  double cooling_rate_SD93CIE(T);
  //  void setup_SD93_cie_OnlyMetals();
  //  void setup_SD93_cie_MetalFree();
  //
  // Public functions from Hummer1994 recombination/cooling class:
  //  double Hii_rad_recomb_rate(T);
  //  double Hii_rad_recomb_cooling(T);
  //  double Hii_total_cooling(T);

  
  private:
  double
    nH, ///< total hydrogen number density at current cell.
    Vshell, ///< geometric factor in point source flux calculation (\sim 4\pi R^2 dR).
    NH0,    ///< column density of neutral hydrogen to front edge of cell.
    dNH0,   ///< column density of neutral hydrogen through cell (=nH*dS*(1-x)).
    G0_UV,  ///< UV heating flux, including attenuation, F*exp(-1.9Av).
    G0_IR,  ///< Heating due to UV flux re-radiated in IR and re-absorbed, F*exp(-0.05Av).
    NIdot,  ///< photon luminosity of monochromatic ionising source (ionising photons/s).
    delta_S,///< path length through current cell.
    gamma,  ///< EOS gamma (usually 5/3).
    gamma_minus_one, k_B, m_p;
  int
    n_eq, ///< number of equations (and hence variables) to solve.
    n_xd, ///< number of variables in extra-data void pointer for cvodes (zero here).
    diff, ///< No diffuse sources --> 0, otherwise --> 1
    ion,  ///< No ionising sources --> 0, otherwise --> 1
    ion_src_type; ///< Either RT_EFFECT_MFION or RT_EFFECT_PION_MONO.

};


// ##################################################################
// ##################################################################



///
/// We need an instance of this class.
///
extern class mp_rates_ExpH_ImpMetals MPR;


// ##################################################################
// ##################################################################





// ##################################################################
// ##################################################################



//
// Header files for CVodes solver.
//
#include <cvode/cvode.h>           // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h>  // serial N_Vector types, fcts., and macros
#include <cvode/cvode_dense.h>     // prototype for CVDense
#include <sundials/sundials_dense.h> // definitions CVMatrix DENSE_ELEM
#include <sundials/sundials_types.h> // definition of type realtype
#include <cstdio>
#include <cmath>

#define MAX_TRYS  50
#define MAX_STEPS 100


// ##################################################################
// ##################################################################



//
// ydot needs to be a stand-alone function for cvodes, so this function will
// just provide the interface to call solver::ydot();
//
int Ydot_for_cvodes(
      double, ///< current time
      N_Vector, ///< current Y-value
      N_Vector,  ///< vector for Y-dot values
      void *    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
      );



// ##################################################################
// ##################################################################


///
/// This class will now just be a driver for CVodes, and using the dYdt() function from 
/// mp_rates_ExpH_ImpMetals MPR.  It is essentially an interface class like MPv9
/// which interfaces between my code and Harpreet's.
///
class MPv2
:
  public microphysics_base
{
  public:
  ///
  /// Constructor.  Sets index for H ion fraction in primitive variable state vector.
  ///
  MPv2(
      const int,  ///< grid dimensions
      const int,  ///< Coordinate System flag
      const int,          ///< Total number of variables in state vector
      const int,          ///< Number of tracer variables in state vector.
      const std::string *, ///< List of what the tracer variables mean.
      struct which_physics *, ///< pointer to extra physics flags.
      struct rad_sources *,    ///< radiation sources.
      const double  ///< EOS Gamma
      );

  ~MPv2();

  ///
  /// The main update function.  A primitive vector is input, and lots of optional 
  /// extra data in the last argument, and the ion fraction(s) and internal energy
  /// are updated by the requested timestep.  Results are output to the destination
  /// vector, which can be the same pointer as the initial vector.
  ///
  /// The integration switch uses the following values:
  /// 0=adaptive-RK5, 1=adaptive-Euler, 2=single-step-RK4
  /// For strongly ionised cells an analytic solution is used, where the ion fraction
  /// relaxes exponentially to its equilibrium value.
  ///
  int TimeUpdateMP(
      const double *, ///< Primitive Vector to be updated.
      double *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use.
      double *        ///< Vector of extra data (column densities, etc.).
      );

  ///
  /// If doing ray-tracing, the tracer can call this function to
  /// integrate the microphysics variables forward one timestep
  /// given an external radiation flux.
  ///
  ///
  int TimeUpdate_RTsinglesrc(
      const double *, ///< Primitive Vector to be updated.
      double *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use.
      const double,   ///< flux in per unit length along ray (F/ds or L/dV)
      const double,   ///< path length ds through cell.
      const double,   ///< Optical depth to entry point of ray into cell.
      double *        ///< return optical depth through cell in this variable.
      );

  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative transfer,
  /// where the column densities for diffuse and direct radiation are included as 
  /// parameters.  The input list of column densities is ordered by the number of 
  /// sources in each category in the vector of integers.
  ///
  /// Integers refer to:
  /// - Number of diffuse ionising sources (at infinity), 
  /// - Number of diffuse UV sources (at infinity),
  /// - Number of ionising point sources,
  /// - Number of UV point sources.
  ///
  int TimeUpdateMP_RTnew(
      const double *, ///< Primitive Vector to be updated.
      const int,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      double *,       ///< Destination Vector for updated values
                      ///< (can be same as first Vector.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use.
                  ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature (not strictly needed).
      );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  /// - Is threadsafe.
  ///
  double Temperature(
      const double *, ///< primitive vector
      const double    ///< eos gamma
      );

  ///
  /// Set the gas temperature to a specified value.
  /// Again only needed if you want this feature in the initial condition generator.
  ///
  int Set_Temp(
      double *,     ///< primitive vector.
      const double, ///< temperature
      const double  ///< eos gamma.
      );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Not implemented for this class, so this will just print
  /// out a warning and return an unconstrainingly large timescale.  Use
  /// the newer timescales interface.
  ///
  double timescales(
      const double *, ///< Current cell.
      const double,   ///< EOS gamma.
      const bool, ///< set to 'true' if including cooling time.
      const bool, ///< set to 'true' if including recombination time.
      const bool  ///< set to 'true' if including photo-ionsation time.
      );

  ///
  /// This returns the minimum timescale of all microphysical processes, including
  /// reaction times for each species and the total heating/cooling time for the gas.
  /// It requires the radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.
  ///
  double timescales_RT(
      const double *, ///< Current cell.
      const int,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      );

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium value.
  /// This is optionally used in the initial condition generator.  Not implemented here.
  ///
  int Init_ionfractions(
      double *, ///< Primitive vector to be updated.
      const double, ///< eos gamma.
      const double  ///< optional gas temperature to end up at. (negative means use pressure)
      );

  ///
  /// Return index of tracer for a given string. (only hydrogen for this class!)
  ///
  int Tr(
      const string ///< name of tracer we are looking for.
      );

  ///
  /// Set the properties of a multifrequency ionising radiation source.
  ///
  int set_multifreq_source_properties(
      const struct rad_src_info *
      );

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
      const int,      ///< ion index in tracer array (optional).
      const pion_flt *, ///< input state vector (primitive).
      const double    ///< EOS gamma (optional)
      );
  
  private:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  int convert_prim2local(
      const double *, ///< primitive vector from grid cell (length nv_prim)
      double *        ///< local vector [x,E](n+1).
      );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  int convert_local2prim(
      const double *, ///< local (updated) vector [x,E](n+1).
      const double *, ///< input primitive vector from grid cell (length nv_prim)
      double *       ///< updated primitive vector for grid cell (length nv_prim)
      );

  ///
  /// This initialises the CVodes solver without using a Jacobian function.
  ///
  int setup_cvodes_solver_without_Jacobian();

  ///
  /// This takes a step dt, returning a non-zero error code if the error fails.
  ///
  int integrate_cvodes_step(
      const std::vector<double> &,  ///< input vector
      double *, ///< parameters for user_data
      double,  ///< start time.
      double,  ///< time-step.
      std::vector<double> &  ///< output vector.
      );

  ///
  /// For each cell, dYdt() needs to know the local radiation field.  This function
  /// takes all of the input radiation sources and calculates what dYdt() needs to
  /// know for both heating and ionisation sources.
  ///
  void setup_radiation_source_parameters(
      const double *, ///< primitive input state vector.
      double *,  ///< local input state vector (x_in,E_int)
      const int , ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &
      ///< list of ionising src column densities and source properties.
      );

  N_Vector 
    y_in,  ///< current y-vector
    y_out, ///< output y-vector
    abstol; ///< vector of absolute error tolerances for y-elements.
  void *cvode_mem; ///< pointer to memory allocation for the solver.
  int n_eq; ///< number of equations to solve.
  int n_xd; ///< number of elements in user-data array.


  std::vector<double> diff_angle; ///< solid angle of each diffuse radiation source.
  //
  // Any data which can be calculated/set at the start of the simulation
  // can be defined here.
  //
  double k_B; ///<  Boltzmanns constant.
  double m_p; ///< Mass of proton.
  double gamma;     ///< EOS gamma for ideal gas.
  double gamma_minus_one; ///< as named.
  double lv_nH;  ///< current number density of H nucleons.
  struct which_physics ep; ///< struct with flags for which extra physics we are or aren't doing.
  double Min_IonFrac; ///< minimum H+ fraction allowed (eps=1.0e-12)
  double Max_IonFrac; ///< Maximum H+ fraction allowed (1-eps)
  double mean_mass_per_H; ///< mean mass per hydrogen nucleon, should be about 2.34e-24;

  const int ndim; ///< Number of dimensions in grid.
  const double eos_gamma; ///< EOS gamma for ideal gas.
  const int coord_sys; ///< Coordinate System flag
  int nvl;     ///< number of variables in local state vector.
  int lv_eint; ///< internal energy local variable index. 
  int lv_Hp;   ///< ionised hydrogeen local variable index.
  int pv_Hp;   ///< index for element of Primitive vector that holds ionisation pot.
  bool have_setup_cvodes; ///< flag to make sure we only set up CVODES once.
};


// ##################################################################
// ##################################################################

#endif // LEGACY_CODE
#endif // MPV2_H

