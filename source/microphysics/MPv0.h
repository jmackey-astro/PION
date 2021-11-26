/// \file MPv0.h
///
/// \brief Class Definitions for MicroPhysics Routines.
/// \author Jonathan Mackey
///
/// This is the original microphysics class written for PION by
/// J.Mackey, but it had some problems (speed, reliability) and was
/// never used for research resulting in publications.  Will probably
/// be removed in future versions of PION.
///
/// Modifications:
///  - 2007-12-17 Started File with CoolingFn, UpdateMicroPhysics classes.
///  - 2008-01-21 Moved Cooling to its own file.  updated interface.
///
///  - 2010-01-05 JM: Added public function for getting microphysics timescales
///      for heating/cooling/recombination etc.
///
///
/// - 2010-04-10 JM: Added flag for Raga's recombination rates (not
///     used at the moment)
///
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2011.03.23 JM: Added new interface functions (they are just placeholders
/// for now).
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.08.05 JM: tidied up code; added pion_flt datatype.
/// - 2018.03.20 JM: changed filename

#ifndef MPV0_H
#define MPV0_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include <spdlog/spdlog.h>

#ifdef LEGACY_CODE

#include "cooling.h"
#include "integrator.h"
#include "microphysics_base.h"
#include <map>
#include <vector>

#define MP_EULERINT 1
#define MP_RK4ORDER 2
#define MP_ADAPTIVE 3

#define NEW_RATES  ///< Rates for microphysics...
//#define RAGA_RATES ///< Rates for microphysics...

/// Assign an int to each species.
enum species {
  NONE = -1,  ///< For if we want to specify no ion.
  H_0  = 0,   ///< neutral hydrogen.
  H_1p = 1,   ///< ionised hydrogen.
  He0  = 2,   ///< neutral helium.
  He1p = 3,   ///< singly ionised helium.
  He2p = 4,   ///< doubly ionised helium.
  C0   = 5,   ///< Carbon
  C1p  = 6,
  C2p  = 7,
  C3p  = 8,
  C4p  = 9,
  C5p  = 10,
  C6p  = 11,
  N0   = 12,  ///< Nitrogen
  N1p  = 13,
  N2p  = 14,
  N3p  = 15,
  N4p  = 16,
  N5p  = 17,
  N6p  = 18,
  N7p  = 19,
  O0   = 20,  ///< Oxygen
  O1p  = 21,
  O2p  = 22,
  O3p  = 23,
  O4p  = 24,
  O5p  = 25,
  O6p  = 26,
  O7p  = 27,
  O8p  = 28
};

/// \brief Contains data for each ion I am using in the code.
///
/// Includes ion's name, and its element's name, name of the lower and higher
/// ionisation stages (if they exist), its ionisation potential (in ergs), its
/// charge (in units of electron charge), an enum identifier (for faster
/// identification than strcmp()), a pointer to its element's data, and it's
/// index in the local state vector (if relevant... not for neutral species).
///
struct ion_struct {
  std::string ion,          ///< Name of ion.
      ip1,                  ///< name of ion with higher ionisation stage.
      im1,                  ///< name of ion with lower ionisation stage.
      el;                   ///< Name of ion's element.
  struct ion_struct *iip1,  ///< pointer to higher stage ion, if present.
      *iim1;                ///< pointer to lower stage ion, if present.
  double ion_pot;  ///< Ionisation energy of ion to next ionisation stage, in
                   ///< ergs.
  int charge;      ///< Charge of ion in units of the electron charge.
  enum species i;  ///< ion id number.
  struct element_struct *e;  ///< pointer to associated element info.
  int index;                 ///< where in the local state vector the ion is.
  int pv_index;  ///< index in primitive variable state vector from/to main
                 ///< code.
  int g_stat;    ///< statistical weight of this ion (H0=2,H1+=1,etc.)
  bool ch_ex;    ///< true if there is a charge exchange reaction with H0/H1+,
                 ///< false if not.
};

/// \brief Contains data for each element I am using in the code.
///
/// Contains the name of the element, a vector of the names of all its
/// species including neutral (e.g. He0,He1+,He2+), the number of
/// species included in state vector, its mass and its number fraction
/// relative to hydrogen (n(E)/n(H)), and a vector listing the index of
/// each ion's number fraction in the local state vector.
///
struct element_struct {
  std::string el;            ///< Name of element.
  std::vector<string> ions;  ///< List of ions (including neutral).
  int nspecies;      ///< Number of ionisation states included in state vec
                     ///< (may not have all ionisation states).
  int *ion_indices;  ///< list of indices of ions in the local state vector.
  double mass;       ///< Mass in units of the proton mass.
  double numfrac;    ///< Number fraction of element in relation to hydrogen.
};

/// \brief Class to hold the info for updating the chemistry and cooling.
///
/// This should have only one publicly callable function, to update a
/// state vector from time t to t+dt.  All the details of how this is
/// done should be hidden.
///
class MPv0 : public microphysics_base, public Integrator_Base {
public:
  /// Constructor.  Takes in info on tracers to determine what sort of
  /// chemistry we are using.
  MPv0(
      const int,               ///< Total number of variables in state vector
      const int,               ///< Number of tracer variables in state vector.
      const std::string,       ///< type of chemistry we are running.
      const std::string *,     ///< List of what the tracer variables mean.
      struct which_physics *,  ///< pointer to extra physics flags.
      struct rad_sources *     ///< radiation sources.
  );

  /// \brief Destructor deletes dynamically allocated member data.
  ~MPv0();

  /// \brief  This takes a copy of the primitive vector and advances it in
  /// time over the step requested, and at the end copies the updated vector
  /// into the destination vector.  For fully local microphysics (no R-T!)
  ///
  int TimeUpdateMP(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use. (0=adaptive
                  ///< RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature.
  );

  /// \brief Same as TimeUpdateMP, except that a incident photon flux
  /// is given, and photoionisation is calculated.
  /// Optical depth through the cell being processed is returned.
  ///
  int TimeUpdate_RTsinglesrc(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use. (0=adaptive
                  ///< RK5, 1=adaptive Euler,2=onestep o4-RK)
      const double,  ///< flux in per unit length along ray (F/ds or L/dV)
      const double,  ///< path length ds through cell.~
      const double,  ///< Optical depth to entry point of ray into cell.
      double *       ///< return optical depth through cell in this variable.
  );

  /// \brief Returns element number of named tracer variable in state vector.
  int Tr(string  ///< tracer we want to get index for;
  );

  /// \brief Initialise microphysics ionisation fractions to an equilibrium
  /// value.
  int Init_ionfractions(
      pion_flt *,    ///< Primitive vector to be updated.
      const double,  ///< eos gamma.
      const double   ///< optional gas temperature to end up at.
                     ///< (negative means use pressure)
  );

  int Set_Temp(
      pion_flt *,    ///< primitive vector.
      const double,  ///< temperature
      const double   ///< eos gamma.
  );

  /// \brief Returns the gas temperature (not very optimized though)
  ///
  /// Assumes primitive vector is in cgs units.
  ///
  double Temperature(
      const pion_flt *,  ///< primitive vector
      const double       ///< eos gamma
  );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Time is returned in seconds.
  ///
  double timescales(
      const pion_flt *,  ///< Current cell.
      const double,      ///< EOS gamma.
      const bool,        ///< set to true if including cooling time.
      const bool,        ///< set to true if including recombination time.
      const bool         ///< set to true if including photo-ionsation time.
  );

  ///
  /// This returns the minimum timescale of all microphysical processes,
  /// including reaction times for each species and the total heating/cooling
  /// time for the gas. It requires the radiation field as an input, so it has
  /// substantially greater capability than the other timescales function.
  ///
  virtual double timescales_RT(
      const pion_flt *,  ///< Current cell.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
  )
  {
    spdlog::info("Don't call timescales for old MP class!");
    return 1.0e99;
  }

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  double get_recombination_rate(
      const int,         ///< ion index in tracer array (optional).
      const pion_flt *,  ///< input state vector (primitive).
      const double       ///< EOS gamma (optional)
  )
  {
    spdlog::info("MPv0::get_recombination_rate: not written!");
    return 1.0e99;
  }

  ///
  /// Return the H mass fraction
  ///
  virtual inline double get_X_H() { return ep.H_MassFrac; }

private:
  /// \brief convert state vector from grid cell into local microphysics
  /// vector.
  int convert_prim2local(
      const pion_flt *,  ///< primitive vector (nv_prim)
      double *,          ///< local vector (length nvl)
      const double       ///< eos gamma.
  );

  /// \brief convert local microphysics vector into state vector for grid
  /// cell.
  int convert_local2prim(
      const double *,    ///< local vector (length nvl)
      const pion_flt *,  ///< input primitive vector(length nv_prim)
      pion_flt *,        ///< updated primitive vector(length nv_prim)
      const double       ///< eos gamma.
  );

  /// \brief Calculate rate of change of local state vector. Note this is
  /// a different vector to the main code state vector, and
  /// consists of n_h, E_int, and if needed, x_e and all the ions x_i.
  ///
  int dPdt(
      const int,       ///< length of state vector (for checking).
      const double *,  ///< current state vector P.
      double *         ///< Rate vector to write to, R=dPdt(P)#
  );

  /// \brief Calculate rate of change of local state vector. Note this is
  /// a different vector to the main code state vector, and
  /// consists only of n_h, E_int.
  ///
  /// Special function for when there are no ions, so we assume fully ionised
  /// hydrogen, and just do the cooling function.
  ///
  int dPdt_OnlyCooling(
      const int,       ///< length of state vector (for checking).
      const double *,  ///< current state vector P.
      double *         ///< Rate vector to write to, R=dPdt(P)
  );

  /// \brief  This takes a copy of the primitive vector and advances it in
  /// time over the step requested, and at the end copies the updated vector
  /// into the destination vector.  For fully local microphysics (no R-T!).
  ///
  /// This function is for when we have no ions, but just do the cooling
  /// assuming the gas is fully ionised hydrogen.
  ///
  int TimeUpdate_OnlyCooling(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,  ///< Switch for what type of integration to use. (0=adaptive
                  ///< RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *    ///< final temperature.
  );

  /// \brief This is for if we are solving the rate equation, and returns the
  /// creation rate of some quantity. xdot=A*(1-x)-B*x, so this returns A(x).
  ///
  int C_rate(
      const int,       ///< length of state vector.
      const double *,  ///< current state vector P.
      double *         ///< Creation rate vector to write to.
  );

  /// \brief This is for if we are solving the rate equation, and returns the
  /// destruction rate of some quantity. xdot=A*(1-x)-B*x, so this returns
  /// A(x)+B(x).
  ///
  int D_rate(
      const int,       ///< length of state vector.
      const double *,  ///< current state vector P.
      double *         ///< Destruction rate vector to write to.
  );

  void set_atomic_data();  ///< sets atomic data arrays to right values.
  const double kB;         ///<  Boltzmanns constant cgs units
  const double m_p;        ///< Mass of proton, cgs units.
  double gamma;            ///< EOS gamma for ideal gas.
  double min_elecf;        ///< minimum electron fraction n(e-)/n(H)
  std::string chemtype;    ///< type of chemistry network we are using.
  //   class Chemistry_Base *chem; ///< Chemistry class, which knows how to
  //   calculate things.
  class CoolingFn *cool;    ///< Pointer to generic cooling function.
  struct which_physics ep;  ///< struct with flags for which extra physics we
                            ///< are or aren't doing.

  int nvl;                     ///< number of variables in local state vector.
  std::map<string, int> pvar;  ///< List of tracer variables in primitive
                               ///< vector matched to their index.
  std::map<string, int> lvar;  ///< List of variables in local state vector
                               ///< (only microphysics vars).
  int lv_nh;                   ///< neutral hydrogen local variable index.
  int lv_eint;                 ///< internal energy local variable index.
  int lv_hp;                   ///< ionised hydrogeen local variable index.
  int lv_elec;                 ///< electron fraction local variable index.
  int pv_elec;                 ///< electron fraction state vector index.
  int lv_dtau;         ///< tracker variable for integral of optical depth.
  double path_length;  ///< for ionisation, the path length to use to
                       ///< calculate tau.
  double tau_cell;     ///< tracker variable for the optical depth.
  double photons_in;   ///< number of incident photons (F/ds or L/dV)

  //
  // These variables are for the FUV heating rate, calculated according to
  // equation A3 in Henney, Arthur, et al. 2009.
  //
  double FUV_unattenuated_flux;  ///< Incident flux of FUV photons (6-13eV) in
                                 ///< units of 1.2e7 per cm2 per sec.
  double FUV_extinction;         ///< Extinction A_v due to dust.
  double
      FUV_attenuated_flux;  ///< FUV_unattenuated_flux*exp(-1.9*FUV_extinction)

  std::vector<string> ions,  ///< list of ions used.
      ion_list;              ///< list of all ions.
  std::vector<string> els,   ///< list of elements used.
      el_list;               ///< list of all elements.
  int nions;                 ///< number of ions used.
  int nels;                  ///< number of elements used.
  std::map<std::string, struct ion_struct>
      ion_props;  ///< properties of each ion.
  std::map<std::string, struct element_struct>
      el_props;                ///< properties of each element.
  struct ion_struct **ii;      ///< pointer to array of ions used.
  struct element_struct **ee;  ///< pointer to array of elements used.

  void copy_element_struct(
      const element_struct,  ///< source.
      element_struct *       /// destination.
  );

  void copy_ion_struct(
      const ion_struct,  ///< source.
      ion_struct *       /// destination.
  );

  void show_ion_struct(const ion_struct *);

  /// Returns number density of hydrogen in all its forms.
  double Get_nH(const double  ///< gas density.
  );

  /// Get number density of all ions from state vector.
  double Get_nIons(const double *  ///< State vector to get number density from
  );

  /// Get number density of all particles from state vector.
  double Get_nTot(const double *  ///< State vector to get number density from
  );

  /// For a give ion, calculate neutral fraction of its element.
  double neutral_fraction(
      const double *,     ///< State vector to get ionised fractions from.
      const ion_struct *  ///< pointer to ion we want to get neutral fraction
                          ///< of.
  );

  /// Returns temperature for a given state vector.
  double Get_Temp(const double *  ///< State vector to get Temperature from.
  );

  /// Given a temperature, set internal energy to consistent value.
  int Set_Eint(
      double *,     ///< local state vector
      const double  ///< Temperature we want to set to.
  );

  /// Returns Collisional Ionisation rate from ion passed in to next stage up,
  /// in [cm^3/s] units.
  double Coll_Ion_rate(
      double,             ///< Precalculated Temperature.
      const ion_struct *  ///< pointer to current ion data.
  );

  /// Returns Radiative Recombination rate from ion passed in to next stage
  /// down, in [cm^3/s] units.
  double Rad_Recomb_rate(
      double,             ///< Precalculated Temperature.
      const ion_struct *  ///< pointer to current ion data.
  );

  /// Some ions have dielectronic recombination; this returns the DR rate,
  /// calculated using rates from Mazzotta et al. (1998).
  double dielec_recomb(
      double,       ///< Temperature
      enum species  ///< which ion we are working on.
  );

  /// Radiative recombination rate, calculated using rates from
  /// Verner & Ferland (1996) and their fits to Pequignot et al. (1991).
  double rad_recomb(
      double,       ///< Temperature
      enum species  ///< which ion we are working on.
  );

  /// returns photoionisation cross section of ion.
  double phot_xsection(const ion_struct *  ///< pointer to current ion data.
  );

  /// Returns charge transfer rate between current ion and hydrogen.
  double charge_exchange_rate(
      double,                     ///< Precalculated Temperature.
      const struct ion_struct *,  ///< current ion.
      const int  ///< 0 for recombination from ion, 1 for ionisation of ion.
  );
};

#endif  // LEGACY_CODE

#endif  // MPV0_H
