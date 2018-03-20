/// \file MPv1.h
/// 
/// \brief Class Definitions for MicroPhysics Routines.
/// \author Jonathan Mackey
///
/// \description
/// This microphysics class was used in Mackey & Lim (2010,2011) for
/// studying the expansion of HII regions and the formation of
/// pillars, globules, elephant trunks at the borders of HII regions.
/// It uses a modified version of the C2-ray algorithm developed by
/// G.Mellema and collaborators (Mellema et al., 2006).
/// It has since been superseded by updated algorithms that scale
/// better on supercomputers (see Mackey, 2012, A&A).
/// 
/// Modifications:
///  - 2007-12-17 Started File with CoolingFn, UpdateMicroPhysics classes.
///  - 2008-01-21 Moved Cooling to its own file.  updated interface.
///  - 2010-01-05 JM: Added public function for getting microphysics timescales
///      for heating/cooling/recombination etc.
/// - 2010-04-10 JM: Added flag for Raga's recombination rates (not
///     used at the moment)
///
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2011.03.23 JM: Added new interface functions (they are just placeholders for now).
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.08.05 JM: tidied up code; added pion_flt datatype.
/// - 2018.03.20 JM: Renamed file.

#ifndef MPV1_H
#define MPV1_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifndef EXCLUDE_MPV1

#include <map>
#include <vector>
#include "microphysics_base.h"
#include "integrator.h"
#include "cooling.h"


#define MP_EULERINT 1
#define MP_RK4ORDER 2
#define MP_ADAPTIVE 3

#define NEW_RATES ///< Rates for microphysics...
//#define RAGA_RATES ///< Rates for microphysics...



//#define ISOTHERMAL_MP
#define HUMMER_RECOMB
//#define NO_DOUBLECOUNTING
//#define MP_DEBUG
/// Class to update the chemistry and cooling for Atomic Hydrogen.
///
/// This is a simplified version of the general solver, streamlined for 
/// updating microphysics of a pure atomic hydrogen gas.  This means we
/// only need to deal with one tracer variable, as the ionisation fraction
/// of hydrogen is the same as the electron fraction.
/// 
class MPv1 : public MicroPhysicsBase, public Integrator_Base {
 public:
  /// Constructor.  Takes in info on tracers to determine what sort of
  /// chemistry we are using.  The string has a specific format described in 
  /// page \ref userguide.
  /// 
  MPv1(
      const int,          ///< Total number of variables in state vector
      const int,          ///< Number of tracer variables in state vector.
      const std::string *, ///< List of what the tracer variables mean.
      struct which_physics *, ///< pointer to extra physics flags.
      struct rad_sources *    ///< radiation sources.
      );

  /// Destructor deletes dynamically allocated member data. 
  ~MPv1();

  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics (no R-T!)
  /// 
  int TimeUpdateMP(
      const pion_flt *, ///< Primitive Vector to be updated.
      pion_flt *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,       ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature.
      );

  /// Same as TimeUpdateMP, except that a incident photon flux
  /// is given, and photoionisation is calculated. optical depth
  /// through the cell being processed is returned.
  /// 
  int TimeUpdate_RTsinglesrc(
      const pion_flt *, ///< Primitive Vector to be updated.
      pion_flt *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
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
  /// parameters.  The input list of column densities is ordered in a list of UV
  /// heating sources and a list of ionising sources.
  /// THIS IS A DUMMY FUNCTION; IT JUST PRINTS A MESSAGE AND RETURNS 0
  ///
  virtual int TimeUpdateMP_RTnew(
      const pion_flt *, ///< Primitive Vector to be updated.
      const int,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      pion_flt *,       ///< Destination Vector for updated values
      ///< (can be same as first Vector.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int, ///< Switch for what type of integration to use.
      ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature (not strictly needed).
        ) {cout <<"MP_Hydrogem::TimeUpdateMP_RTnew: I don't do anything!\n";return 0;}

  /// Returns element number of named tracer variable in state vector. 
  int Tr(
      string ///< tracer we want to get index for;
      );

  /// Initialise microphysics ionisation fractions to an equilibrium value. 
  int Init_ionfractions(
      pion_flt *, ///< Primitive vector to be updated.
      const double, ///< eos gamma.
      const double  ///< optional gas temperature to end up at. (negative means use pressure)
      );

  /// Set the gas temperature to a specified value. 
  int Set_Temp(
      pion_flt *,     ///< primitive vector.
      const double, ///< temperature
      const double  ///< eos gamma.
      );

  /// \brief Returns the gas temperature (not very optimized though) 
  ///
  /// - Assumes primitive vector is in cgs units.
  /// - Is threadsafe.
  ///
  double Temperature(
      const pion_flt *, ///< primitive vector
      const double    ///< eos gamma
      );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Time is returned in seconds.
  ///
  double timescales(
      const pion_flt *, ///< Current cell.
      const double,   ///< EOS gamma.
      const bool, ///< set to true if including cooling time.
      const bool, ///< set to true if including recombination time.
      const bool  ///< set to true if including photo-ionsation time.
      );

  ///
  /// This returns the minimum timescale of all microphysical processes, including
  /// reaction times for each species and the total heating/cooling time for the gas.
  /// It requires the radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.
  ///
  virtual double timescales_RT(
      const pion_flt *, ///< Current cell.
      const int,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      ) {cout <<"Don't call timescales for old MP class!\n";return 1.0e99;}

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
      const int,      ///< ion index in tracer array (optional).
      const pion_flt *, ///< input state vector (primitive).
      const double    ///< EOS gamma (optional)
      )
  {
    cout <<"MPv1::get_recombination_rate: not written!\n";
    return 1.0e99;
  }

 private:
  /// convert state vector from grid cell into local microphysics vector. 
  int convert_prim2local(
      const pion_flt *, ///< primitive vector from grid cell (length nv_prim)
      double *,       ///< local vector (length nvl)
      const double    ///< eos gamma.
      );

  /// convert local microphysics vector into state vector for grid cell. 
  int convert_local2prim(
      const double *, ///< local vector (length nvl)
      const pion_flt *, ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *,       ///< updated primitive vector for grid cell (length nv_prim)
      const double    ///< eos gamma.
      );

  /// Calculate rate of change of local state vector. Note this is 
  //// certainly a different vector to the main code state vector, and 
  //// consists of n_h, E_int, and if needed, x_e and all the ions x_i.
  /////
  int dPdt(
      const int ,     ///< length of state vector (for checking).
      const double *, ///< current state vector P.
      double *        ///< Rate vector to write to, R=dPdt(P)
      );

  /// Calculate rate of change of local state vector. Note this is 
  //// certainly a different vector to the main code state vector, and 
  //// consists only of n_h, E_int.  
  ////
  //// Special function for when there are no ions, so we assume fully ionised
  //// hydrogen, and just do the cooling function.
  /////
  int dPdt_OnlyCooling(
      const int ,     ///< length of state vector (for checking).
      const double *, ///< current state vector P.
      double *        ///< Rate vector to write to, R=dPdt(P)
      );

  ///  This takes a copy of the primitive vector and advances it in time over
  //// the step requested, and at the end copies the updated vector into the
  //// destination vector.  For fully local microphysics (no R-T!).
  ////
  //// This function is for when we have no ions, but just do the cooling assuming the 
  //// gas is fully ionised hydrogen.
  //// 
  int TimeUpdate_OnlyCooling(
      const pion_flt *, ///< Primitive Vector to be updated.
      double *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature.
      );

#ifdef COUNT_ENERGETICS
  bool have_counted_ergs;
#endif
  const double kB;  ///<  Boltzmanns constant cgs units
  const double m_p; ///< Mass of proton, cgs units.
  double gamma;     ///< EOS gamma for ideal gas.
  double min_elecf; ///< minimum electron fraction n(e-)/n(H)
  class CoolingFn      *cool; ///< Pointer to generic cooling function.
  struct which_physics ep; ///< struct with flags for which extra physics we are or aren't doing.

  const int nv_prim; ///< Number of variables in state vector.
  int       nvl;     ///< number of variables in local state vector.
  int lv_nh;   ///< neutral hydrogen local variable index.
  int lv_eint; ///< internal energy local variable index. 
  int lv_Hp;   ///< ionised hydrogeen local variable index.
  int lv_dtau; ///< tracker variable for integral of optical depth.
  int pv_Hp;   ///< tracker for element of Primitive vector that holds ionisation pot.
  double ion_pot; ///< ionisation potential of Hydrogen (13.6eV).
  double path_length; ///< for ionisation, the path length to use to calculate tau.
  double tau_cell; ///< tracker variable for the optical depth.
  double photons_in; ///< number of incident photons (F/ds or L/dV)

  //
  // These variables are for the FUV heating rate, calculated according to equation A3 in 
  // Henney, Arthur, et al. 2009.
  //
  double FUV_unattenuated_flux; ///< Incident flux of FUV photons (6-13eV) in units of 1.2e7 per cm2 per sec.
  double FUV_extinction;        ///< Extinction A_v due to dust.
  double FUV_attenuated_flux;   ///< FUV_unattenuated_flux*exp(-1.9*FUV_extinction)

  double i_crit; ///< Max. ionisation fraction to get to in explicit integration.
  double errtol; ///< error tolerance for integrations.
#ifdef HUMMER_RECOMB
  int hr_nspl; ///< number of interpolation points (31)
  double 
   *hr_t,
   *hr_alpha,
   *hr_alpha2,
   *hr_beta,
   *hr_beta2; ///< arrays for recomb rate and energy loss rate.
#endif // HUMMER_RECOMB
  /// \brief Returns Collisional Ionisation rate from H to H+, in [cm^3/s] units. 
  double coll_ion_rate(
      double ///< Precalculated Temperature.
      );

  ///\brief Returns energy loss (as positive number) in ergs per collisional ionisation. 
  double coll_ion_energy(
      double ///< Precalculated Temperature.
      );

  /// \brief Returns Radiative Recombination rate from H+ to H, in [cm^3/s] units. 
  double rad_recomb_rate(
      double ///< Precalculated Temperature.
      );

  /// \brief Returns energy loss per recombination in ergs (positive number) for Case B. 
  double rad_recomb_energy(
      double ///< Precalculated Temperature.
      );

  /// \brief returns photoionisation cross section of ion. 
  double phot_xsection(
      double ///< Precalculated Temperature.
      );

  /// \brief Returns energy gain per photo-ionisation in ergs. 
  double phot_ion_energy(
      double ///< Precalculated Temperature.
      );

  /// \brief for strongly driven photo-ionisation, take an implicit
  /// step by dt, and get to a specified error tolerance.
  ///
  int implicit_step(
      const int,      ///< number of variables we are expecting.
      const double,   ///< timestep to use.
      const double *, ///< Current state vector.
      double *,       ///< Final State Vector.
      const double    ///< error tolerance
      );

  /// \brief Do an adaptive 5th order Cash-Karp integration to a given relative accuracy. 
  /// 
  /// Pointer to final state can be same as pointer to input state.
  /// This is based on the algorithm description in Numerical Recipes in C (1992), ch16.2.
  /// It is assumed that the Chemistry class contains a function to calculate the 
  /// rate dPdt().
  ///
  int Int_Adaptive_RKCK(
      const int,      ///< number of elements in P array.
      const double *, ///< value of P at initial value of t.
      const double,   ///< initial value of t.
      const double,   ///< Total step dt to take.
      const double, ///< Required fractional accuracy.
      double *,  ///< pointer to final P value.
      double *  ///< pointer to final t value.
      );

  /// This is for if we are solving the rate equation, and returns the
  /// creation rate of some quantity. xdot=A*(1-x)-B*x, so this returns A(x).
  ///
  int C_rate(
      const int,      ///< length of state vector.
      const double *, ///< current state vector P.
      double *        ///< Creation rate vector to write to.
      ) {return -99;}

  /// This is for if we are solving the rate equation, and returns the
  /// destruction rate of some quantity. xdot=A*(1-x)-B*x, so this returns A(x)+B(x).
  ///
  int D_rate(
      const int,      ///< length of state vector.
      const double *, ///< current state vector P.
      double *        ///< Destruction rate vector to write to.
      ) {return -99;}
};


#endif // if not excluding MPv1

#endif // MPV1_H
