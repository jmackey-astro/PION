///
/// \file mp_only_cooling.h
/// \author Jonathan Mackey
/// \date 14.01.2011
///
/// This is a simple microphysics class which only has heating and
/// cooling processes (i.e. no species are tracked explicitly).
///
/// Modifications:\n
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE flag. 
/// - 2011.04.12 JM: Added WSS09 cooling function.  Cooling_flag==8
///    is the recommended cooling function, since it 
///    self-consistently models the hydrogen ioniation-heating and
///    recombination-cooling.  The normalisation is lower than the
///    SD93 curves, possibly because the Oxygen abundance is now
///    lower than it was 15 years ago (Lodders, 2003, ApJ).  Removed
///    the DO_HEATING ifdef -- now we have 5 different cooling
///    functions.
/// - 2013.02.14 JM: Tidied up file.
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2018.01.25 JM: added functions to request n(H+),n(H0),n(e-)

#ifndef MP_ONLY_COOLING_H
#define MP_ONLY_COOLING_H



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "microphysics/cooling.h"
#include "microphysics/cooling_SD93_cie.h"
#include "microphysics/microphysics_base.h"
#include "microphysics/hydrogen_recomb_Hummer94.h"


///
/// This is a simple microphysics class which only has heating and
/// cooling processes (i.e. no species are tracked explicitly).
///
class mp_only_cooling 
 :
  virtual public microphysics_base,
  virtual public cooling_function_SD93CIE,
  virtual public Hummer94_Hrecomb,
  virtual public CoolingFn
{
  public :
  ///
  /// constructor.
  ///
  mp_only_cooling(
      const int, ///< length of prim. state vec.
      struct which_physics *, ///< pointer to extra physics flags.
      struct rad_sources *    ///< radiation sources.
      );

  ///
  /// destructor.
  ///
  ~mp_only_cooling();

  ///
  /// Non-RT microphysics update, so cooling and heating.
  /// This uses various integration methods to update the elements
  /// and the internal energy.
  ///
  virtual int TimeUpdateMP(
        const pion_flt *, ///< Primitive Vector to be updated.
        pion_flt *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use (not used here!)
        double *        ///< final temperature.
        );

  ///
  /// This is just a dummy function which bugs out if called.
  ///
  virtual int TimeUpdate_RTsinglesrc(
        const pion_flt *,  ///< Primitive Vector to be updated.
        pion_flt *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
        const double,   ///< flux in per unit length along ray (F/ds or L/dV)
        const double,   ///< path length ds through cell.
        const double,   ///< Optical depth to entry point of ray into cell.
        double *        ///< return optical depth through cell in this variable.
        ) {return DONT_CALL_ME;}

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value.  This function does nothing here.
  ///
  virtual int Init_ionfractions(
        pion_flt *, ///< Primitive vector to be updated.
        const double, ///< eos gamma.
        const double  ///< optional gas temperature to end up at. (negative means use pressure)
        ) {return DONT_CALL_ME;}

  ///
  /// \brief Given a cell, calculated converged time-averaged rates
  /// of change for species abundances (x(E)) and heating/cooling,
  /// including radiative effects.
  /// 
  /// The rates are to be stored in the Ph vector of the cell.  This
  /// means that this function must not be called on a half-timestep,
  /// and is supposed to be completely operator split from the
  /// dynamics.  It should only be called by the raytracing class.
  /// 
  /// The idea is that for a infinitesimally weak radiation source,
  /// this will produce an identical update to the non-radiative
  /// TimeUpdateMP function, once the time-averaged rates are added
  /// to the primitive vector to get the updated value.
  ///
  virtual int Get_MP_Rates(
        cell *,       ///< current cell.
        const double, ///< timestep
        const double, ///< column to cell
        const double, ///< Flux entering cell?
        double *      ///< initial column through cell, to be changed to time-averaged column.
        ) {return DONT_CALL_ME;}

  ///
  /// Given a string, return the index of the tracer that it refers
  /// to. Not used here.
  ///
  virtual int Tr(
        const string ///< name of tracer we are looking for.
        ) {return -1;}

  ///
  /// Sets the pressure in p-vec to be such that temperature is what
  /// you want it to be.  Here we assume ionised gas with mu=0.7m_p.
  ///
  virtual int Set_Temp(
        pion_flt *, ///< primitive vector.
        const double, ///< temperature requested.
        const double ///< eos gamma.
        );

  ///
  /// Returns the gas temperature.  Assumes primitive vector is in
  /// cgs units and ionised gas with mu=0.7m_p.
  /// - Is threadsafe.
  ///
  virtual double Temperature(
        const pion_flt *, ///< primitive vector
        const double    ///< eos gamma
        );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Time is returned in seconds.  Only cooling flag
  /// has any effect here.
  ///
  virtual double timescales(
        const pion_flt *, ///< Current cell.
        const double,   ///< EOS gamma.
        const bool, ///< set to true if including cooling time.
        const bool, ///< set to true if including recombination time.
        const bool  ///< set to true if including photo-ionsation time.
        );

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.  It requires the
  /// radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.  But for this
  /// class, which by definition has no radiative transfer
  /// capability, it doesn't do anything!
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
        )
  {
    cout <<"Don't call timescales for only-cooling class!\n";
    return 1.0e99;
  }

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
    cout <<"get_recombination_rate: mp_only_cooling has no ions!\n";
    return 1.0e99;
  }

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_elec(
        const pion_flt * ///< primitive state vector.
        );

  ///
  /// Get H+ number density (cm^{-3})
  ///
  double get_n_Hplus(
        const pion_flt * ///< primitive state vector.
        );

  ///
  /// Get neutral H number density (cm^{-3})
  ///
 double get_n_Hneutral(
        const pion_flt * ///< primitive state vector.
        );

 

 private:
  int cooling_flag;
  double Mu; ///< mean mass per Hydrogen atom.
  double inv_Mu2; ///< One over Mu squared.
  double Mu_elec; ///< mean mass per electron.
  double inv_Mu2_elec_H; ///< One over (Mu*Mu_elec)
  double Mu_ion;  ///< mean mass per ion.
  double Mu_tot;  ///< mean mass per particle: Mu_tot = rho/(n[e-]+n[ions])
  double Mu_tot_over_kB; ///< Mu_tot/kB.
  const int nv_prim; ///< number of variables in primitive state vec
  double MinT_allowed; ///< minimum sensible temperature (from pfile).
  double MaxT_allowed; ///< maximum sensible temperature (from pfile).

  ///
  /// Calculate the rate of change of internal energy density (result
  /// is returned in units of erg/cm3/s).  This function just goes
  /// through a switch/case list of functions for the different
  /// cooling functions used.
  ///
  double Edot(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );

  ///
  /// only cooling, uses SD93-CIE curve only.
  ///
  double Edot_SD93CIE_cool(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );

  ///
  /// cooling from SD93-CIE plus heating assuming full ionisation of
  /// H, where the heating rate equals recombination rate times
  /// 5eV/ionisation.
  ///
  double Edot_SD93CIE_heat_cool(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );

  ///
  /// only cooling, uses WSS09-CIE cooling function.
  ///
  double Edot_WSS09CIE_cool(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );

  /// 
  /// cooling from WSS09-CIE plus heating assuming full ionisation of
  /// H, where the heating rate equals recombination rate times
  /// 5eV/ionisation.
  ///
  double Edot_WSS09CIE_heat_cool(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );

  ///
  /// A bit more self-consistent: metal cooling due to WSS09-CIE or
  /// photoionised metal lines (from Henney et al 2009) -- we take
  /// the max of these two.  Hydrogen cooling from Hummer 1994 
  /// recombination+Bremsstrahlung rate.  Heating from 5eV times the
  /// Hummer 1994 recombination rate.
  ///
  double Edot_WSS09CIE_heat_cool_metallines(
        const double, ///< mass density (g/cm3)
        const double  ///< Temperature (K)
        );
  
};

/************************ MICROPHYSICS ***************************/
#endif // MP_ONLY_COOLING_H
