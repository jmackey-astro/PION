///
/// \file mpv8_TwoTempIso.h
/// \author Jonathan Mackey
/// \date 2013.08.13
///
/// This class is for running calculations with a simple 
/// two-temperature isothermal equation of state, where T is T_low
/// when gas is neutral, and T_high when ionised, and linearly
/// interpolated for partial ionisation.  It also assumes ionisation
/// equilibrium using the Stroemgren integral approach to calculate
/// the ionisation fraction for a monochromatic point source of
/// radiation.
///
/// - THE CODE IS NOT TESTED AGAINST STANDARD SOLUTIONS.
/// - THE CODE IS NOT TESTED IN MULTI-D, SO THE INTERPOLATION MAY BE
///    SKETCHY.
/// - THE CODE IS VERY UNSTABLE AT THE I-FRONT, DRIVING A FORCED
///    OSCILLATION INTO THE HII REGION INTERIOR!
///
/// Modifications:
/// - getting it written: mods up until 2013.08.13


#ifndef MPV8_PUREH_H
#define MPV8_PUREH_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "microphysics_base.h"

//#include "microphysics/mp_explicit_H.h"
#define JM_MINNEU 1.0e-12  ///< minimum neutral fraction i care about.
#define JM_MINERG 1.0e-17  ///< Minimum internal energy density I care about.

class mpv8_TwoTempIso
: public microphysics_base
//  :
//  public mp_explicit_H
{
  public:
  ///
  /// Constructor
  ///
  mpv8_TwoTempIso(
          const int,           ///< Total number of variables in state vector
	  const int,           ///< Number of tracer variables in state vector.
	  const std::string &, ///< List of what the tracer variables mean.
          struct which_physics * ///< extra physics stuff.
	  );

  ///
  /// Destructor
  ///
  ~mpv8_TwoTempIso();


  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW -------------------
  //---------------------------------------------------------------------------
  public:



  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative transfer,
  /// where the column densities for diffuse and direct radiation are included as 
  /// parameters.  The input list of column densities is ordered in a list of UV
  /// heating sources and a list of ionising sources.
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
          const int, ///< Switch for what type of integration to use.
          ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
          double *    ///< final temperature (not strictly needed).
          );

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
          const int,      ///< ion index in tracer array (optional).
          const double *, ///< input state vector (primitive).
          const double    ///< EOS gamma (optional)
          );


  ///
  /// The NON-RT MICROPHYSICS update function.
  /// Should never be called with this class.
  ///
  int TimeUpdateMP(
        const double *, ///< Primitive Vector to be updated.
        double *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use.
        double *        ///< Vector of extra data (column densities, etc.).
        ) {return DONT_CALL_ME;}
  
  ///
  /// This is only for the implicit time-integrator.
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
        ) {return DONT_CALL_ME;}

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
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
  /// arguments.  Not implemented here because it is only for sims
  /// without radiative transfer.
  ///
  virtual double timescales(
          const double *, ///< Current cell.
          const double,   ///< EOS gamma.
          const bool, ///< set to 'true' if including cooling time.
          const bool, ///< set to 'true' if including recombination time.
          const bool  ///< set to 'true' if including photo-ionsation time.
          ) {return DONT_CALL_ME;}

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, but since we assume ionisation and thermal equil. we
  /// can just return a very large number.
  ///
  virtual double timescales_RT(
        const double *, ///< Current cell.
        const int,      ///< Number of UV heating sources.
        const std::vector<struct rt_source_data> &,
        ///< list of UV-heating column densities and source properties.
        const int,      ///< number of ionising radiation sources.
        const std::vector<struct rt_source_data> &,
        ///< list of ionising src column densities and source properties.
        const double   ///< EOS gamma.
        ) {return 1.0e99;}

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value.  Not implemented here.
  ///
  int Init_ionfractions(
        double *, ///< Primitive vector to be updated.
        const double, ///< eos gamma.
        const double  ///< optional gas temperature to end up at.
        ) {return DONT_CALL_ME;}

  ///
  /// Return index of tracer for a given string. (only hydrogen!)
  ///
  int Tr(const string ///< name of tracer we are looking for.
	) {return lv_H0;}

  protected:

  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  virtual int convert_prim2local(
            const double *, ///< primitive vector from grid cell (length nv_prim)
            double *        ///< local vector [x(H0),E](n+1).
            );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
            const double *, ///< local (updated) vector [x(H0),E](n+1).
            const double *, ///< input primitive vector from grid cell (length nv_prim)
            double *       ///< updated primitive vector for grid cell (length nv_prim)
            );

  ///
  /// returns gas temperature according to ion fraction
  ///
  virtual double get_temperature(
    const double, ///< nH
    const double, ///< E_int
    const double  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
    );

  ///
  /// Returns total number of particles.
  ///
  double get_ntot(
    const double, ///< nH
    const double  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
    );

  /// Number of neutral particles per neutral H nucleon.
  double TTI_Nnt;
  double TTI_Mol; ///< 0.5 if neutral H is molecular, otherwise 1.0.
  /// Low temperature for fully neutral medium.
  double TTI_Tlo;
  /// High temperature for fully ionised medium.
  double TTI_Thi;

  double k_B; ///<  Boltzmanns constant.
  double m_p; ///< Mass of proton.
  double gamma;     ///< EOS gamma for ideal gas.
  double gamma_minus_one; ///< as named.
  struct which_physics *EP; ///< struct with flags for which extra physics we are or aren't doing.
  double Min_NeutralFrac; ///< minimum H0 fraction allowed (eps=1.0e-12)
  double Max_NeutralFrac; ///< Maximum H0 fraction allowed (1-eps)
  double mean_mass_per_H; ///< mean mass per hydrogen nucleon, should be about 2.34e-24;
  double JM_NELEC; ///< Number of electrons per ionised H atom.
  double JM_NION;  ///< Number of ions per ionised H atom.

  const int nv_prim; ///< Number of variables in state vector.
  int       nvl;     ///< number of variables in local state vector.
  int lv_nH;   ///< number density of H nucleons, local var. index.
  int lv_eint; ///< internal energy local variable index. 
  int lv_H0;   ///< neutral hydrogen local variable index.
  int pv_Hp;   ///< primitive vector index for H+ fraction.


  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS -------------------
  //---------------------------------------------------------------------------
};

#endif // MPV8_PUREH_H



