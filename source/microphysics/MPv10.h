///
/// \file MPv10.h
/// \author Jonathan Mackey
/// \date 2011.10.06
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
///
/// Next generation modifications:
/// - 2011.10.06 JM: NEW FILE!! Getting it written.  Based on mp_v2_aifa.h
/// - 2011.10.12 JM: Changed local variable from ion fraction to neutral fraction.
/// - 2011.10.17 JM: Debugging.
/// - 2012.04.19 JM: Added "PUREHYDROGEN" ifdef for the Iliev et al. 2009 tests.
/// - 2013.02.14 JM: Started modifying this to include He/Metal mass
///    fractions as EP parameters, to make metallicity and mu into
///    variables that can be set from the parameterfile.
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2014.09.22 JM: Added  total_cooling_rate() function to get the
///    cooling rates per cell for postprocessing.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2018.03.20 JM: Renamed file.

#ifndef MPv10_H
#define MPv10_H

/// Description:
/// This class is an update on the microphysics class used for JM's
/// thesis.  It uses an explicit raytracing scheme, so the photon
/// transmission through the cell does not have to be integrated with
/// the rate and energy equations.
///
/// - It uses multi-frequency photoionisation including spectral
///   hardening with optical depth, using the method outlined in
///   Frank & Mellema (1994,A&A,289,937), and updated by Mellema et
///   al. (2006,NewA,11,374).
/// - It uses the Hummer (1994,MN,268,109) rates for radiative
///   recombination and its associated cooling, and Bremsstrahlung
///   cooling.
/// - For collisional ionisation the function of Voronov
///   (1997,ADNDT,65,1) is used.
/// - Collisional excitation of neutral H, table from Raga, Mellema,
///   Lundqvist (1997,ApJS,109,517), using data from Aggarwal (1993).
/// - Then I use various formulae from Henney et al. (2009,MN,398,
///  157) for cooling due to collisional excitation of photoionised
///  O,C,N (eq.A9), collisional excitation of neutral metals
///  (eq.A10), and the Wiersma+2009 CIE metals-only curve.  I take
///  the max of the WSS09 curve and the Henney et al. functions, to
///  avoid double counting.  For neutral gas I use heating and
///  cooling rates from Wolfire et al. (2003).
/// - Photoheating from ionisation is discussed above.  Cosmic ray
///  heating will use a standard value, X-ray heating is ignored.
///  UV heating due to the interstellar radiation field (ISRF) is
///  according to the optical depth from the edge of the domain to
///  the cell in question, using e.g. HEA09 eq.A3, if requested.  UV
///  heating from the star uses the same equation, but with the
///  optical depth from the source (using a total H-nucleon column
///  density).
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from 
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <vector>
#include <set>
#include <sstream>
#include <cstring>
#include <iostream>
#include "microphysics_base.h"
#include "cooling_SD93_cie.h"
#include "hydrogen_mp.h"
#include "hydrogen_recomb_Hummer94.h"
#include "hydrogen_photoion.h"
#include "cvode_integrator.h"

/// this is the max change in x or E which is integrated with Euler integration.
/// A larger change will be integrated with backward-differencing and Newton
/// iteration, which is more accurate, stable, but expensive.
#define EULER_CUTOFF 0.05


#define JM_RELTOL 1.0e-4   ///< relative-error tolerance (actual error can be larger).
#define JM_MINNEU 1.0e-20  ///< minimum neutral fraction i care about.
#define JM_MINERG 1.0e-17  ///< Minimum internal energy density I care about.


///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of hydrogen, and for the internal energy density.
/// It is explicit in the sense of the raytracer -- the column density is 
/// not time-averaged by the microphysics integrator but rather is an 
/// instantaneous value.
///
class MPv10
  :
  public Hydrogen_chem,
  public cooling_function_SD93CIE,
  public microphysics_base,
  public cvode_solver
{
  public:
  ///
  /// Constructor
  ///
  MPv10(
      const int,  ///< grid dimensions
      const int,  ///< Coordinate System flag
      const int,  ///< Total number of variables in state vector
      const int,  ///< Number of tracer variables in state vector.
      const std::string *, ///< List of what the tracer variables mean.
      struct which_physics *, ///< pointer to extra physics flags.
      struct rad_sources *,    ///< radiation sources.
      const double  ///< EOS Gamma
      );

  ///
  /// Destructor
  ///
  ~MPv10();


  ///
  /// The NON-RT MICROPHYSICS update function.
  /// THIS FUNCTION JUST CALLS TimeUpdateMP_RTnew() WITH NO RADIATION SOURCES.
  ///
  /// A primitive vector is input, and lots of optional extra data in the
  /// last argument, and the ion fraction(s) and internal energy are updated
  /// by the requested timestep.  Results are output to the destination vector,
  /// which can be the same pointer as the initial vector.
  ///
  int TimeUpdateMP(
      const pion_flt *, ///< Primitive Vector to be updated.
      pion_flt *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use.
      double *        ///< Vector of extra data (column densities, etc.).
      );

  ///
  /// UNUSED FUNCTION!!
  ///
  int TimeUpdate_RTsinglesrc(
      const pion_flt *, ///< Primitive Vector to be updated.
      pion_flt *,       ///< Destination Vector for updated values.
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int,      ///< Switch for what type of integration to use.
      const double,   ///< flux in per unit length along ray (F/ds or L/dV)
      const double,   ///< path length ds through cell.
      const double,   ///< Optical depth to entry point of ray into cell.
      double *        ///< return optical depth through cell in this variable.
      )
  {cout <<"MPv10::TimeUpdate_RTsinglesrc() is not implemented!\n";return 1;}

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
      double *    ///< any returned data (final temperature?).
      );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  /// - Is threadsafe.
  ///
  double Temperature(
      const pion_flt *, ///< primitive vector
      const double    ///< eos gamma
      );

  ///
  /// Set the gas temperature to a specified value.
  /// Again only needed if you want this feature in the initial condition generator.
  ///
  int Set_Temp(
      pion_flt *,     ///< primitive vector.
      const double, ///< temperature
      const double  ///< eos gamma.
      );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Not implemented for this class, so this will just print
  /// out a warning and return an unconstrainingly large timescale.  Use
  /// the newer timescales interface.
  ///
  virtual double timescales(
      const pion_flt *, ///< Current cell.
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
  virtual double timescales_RT(
      const pion_flt *, ///< Current cell.
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
      pion_flt *, ///< Primitive vector to be updated.
      const double, ///< eos gamma.
      const double  ///< optional gas temperature to end up at. (negative means use pressure)
      )
  {cout <<"MPv10::Init_ionfractions() is not implemented! Write me!\n";return 1;}

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
  /// Get the total cooling rate.  This is for postprocessing the
  /// simulation data only -- IT IS NOT OPTIMISED FOR SPEED.
  ///
  virtual double total_cooling_rate(
      const pion_flt *, ///< Current cell values.
      const int,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      );

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  virtual double get_recombination_rate(
      const int,      ///< ion index in tracer array (optional).
      const pion_flt *, ///< input state vector (primitive).
      const double    ///< EOS gamma (optional)
      );
  
  ///
  /// Return the H mass fraction
  ///
  virtual inline double get_X_H()
    {return EP->H_MassFrac;}

  protected:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  virtual int convert_prim2local(
      const pion_flt *, ///< primitive vector from grid cell (length nv_prim)
      double *        ///< local vector [x(H0),E](n+1).
      );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
      const double *, ///< local (updated) vector [x(H0),E](n+1).
      const pion_flt *, ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *       ///< updated primitive vector for grid cell (length nv_prim)
      );

  ///
  /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
  /// appropriate for a gas with 10% Helium by number, and if He is singly ionised
  /// whenever H is.
  ///
  virtual double get_temperature(
      const double, ///< nH
      const double, ///< E_int
      const double  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
      );

  ///
  /// Returns total number of particles.
  ///
  virtual double get_ntot(
    const double, ///< nH
    const double  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
    );

  ///
  /// For each cell, ydot() needs to know the local radiation field.  This function
  /// takes all of the input radiation sources and calculates what ydot() needs to
  /// know for both heating and ionisation sources.
  ///
  void setup_radiation_source_parameters(
      const pion_flt *, ///< primitive input state vector.
      double *,  ///< local input state vector (x_in,E_int)
      const int , ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &
      ///< list of ionising src column densities and source properties.
      );

  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  virtual void setup_local_vectors();

  ///
  /// Sets the solid angles associated with diffuse sources at infinity (they are not 
  /// identical for non-cartesian geometries).
  ///
  void setup_diffuse_RT_angle();

  //
  // ********* STUFF FROM THE mp_v2_aifa CLASS **********
  //
  N_Vector 
    y_in,  ///< current y-vector
    y_out; ///< output y-vector
  int N_equations; ///< number of equations i.e. number of species in Y.
  int N_extradata; ///< number of elements in user-data array.


  std::vector<double> diff_angle; ///< solid angle of each diffuse radiation source.
  //
  // Any data which can be calculated/set at the start of the simulation
  // can be defined here.
  //
  double k_B; ///<  Boltzmanns constant.
  double m_p; ///< Mass of proton.
  double m_H; ///< Mass of hydrogen (grams).
  double m_He; ///< Mass of helium.
  const int ndim; ///< Number of dimensions in grid.
  const int nv_prim; ///< Number of variables in state vector.
  const double eos_gamma; ///< EOS gamma for ideal gas.
  const int coord_sys; ///< Coordinate System flag
  double gamma_minus_one; ///< as named.
  double Min_NeutralFrac; ///< minimum H0 fraction allowed (eps=1.0e-12)
  double Max_NeutralFrac; ///< Maximum H0 fraction allowed (1-eps)
  double mean_mass_per_H; ///< mean mass per hydrogen nucleon, should be about 2.34e-24;
  double JM_NELEC; ///< Number of electrons per ionised H atom.
  double JM_NION;  ///< Number of ions per ionised H atom.
  double METALLICITY; ///< Metallicity of gas, in units of solar.

  int       nvl;     ///< number of variables in local state vector.
  int lv_eint; ///< internal energy local variable index.
  int lv_H0;   ///< neutral hydrogeen local variable index. 
  
  // NOTE \Maggie{ can obtain analogous to lv_H0 by using y_ion_index[i] - N_species.
  int N_elem;
  int N_species;
  int N_eqns; ///< := n species
  int N_cons_eqns; /// < := n_species + 1, due to E_int.
  
  
  // NOTE \Maggie{ Added vector variables to keep track of tracers.}
  
  //**********************************************************************
  //******************* PRIMITIVE VECTOR INDICES *************************
  //**********************************************************************
  std::vector<int> y_ion_index; ///<index matching y_ion fration, analogous to pv_Hp before.
  std::vector<int> y_ion_num_elec; ///<index matching number of electrons corresponding to each y_ion (e.g. C6+ = 6)
  std::vector<int> y_ion_mass_frac_index; ///<same length as y_ion_index, records X_elem index as corresponding to each y_ion
  std::vector<pion_flt> y_ion_number_density; ///<X number density at current cell, organised to match y_ion.
  std::vector<pion_flt> y_elem_atom_mass; ///<same length as y_ion_index, records the proper mass of each element, analogous to m_p

  std::vector<float> X_atom_mass; /// < vector of atomic masses corresponding to X_elem_index, e.g. for ["H", "He"], X_atom_mass=[1.6738e-24, 6.6464764e-24 
  std::vector<int> X_elem_index; /// < primitive vector ineices, used to trace X_H etc, like pv_Hp. 
  std::vector<int> N_species_by_elem; ///< records # species in each element, to iterate over later.
  std::set<std::string> set_elem; ///<set of element characters e.g. {"C", "He"}. Used to get Nelem from user's tracer input.
  
  //Need following arrays for MPv10::Tr(), i.e. to get the index of a species based on user input string tracer.
  std::vector<int> H_ion_index; //Locates position of ion with N+1 electrons missing, e.g. H_ion_index[0] -> H+ position.
  std::vector<int> He_ion_index;
  
  int pv_H1p;    ///< legacy, should ideally remove

  int
    N_diff_srcs, ///< No diffuse sources --> 0, otherwise --> 1
    N_ion_srcs,  ///< No ionising sources --> 0, otherwise --> 1
    ion_src_type; ///< Either RT_EFFECT_PION_MULTI or RT_EFFECT_PION_MONO.

  //---------------------------------------------------------------------------
  //-------------- FUNCTIONS DERIVED FROM BASE CLASS FOLLOW -------------------
  //---------------------------------------------------------------------------
  public:
  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  virtual int ydot(
      double,         ///< current time (probably not needed for rate equations)
      const N_Vector, ///< current Y-value
      N_Vector,       ///< vector for Y-dot values
      const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
      );

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  ///
  virtual int Jacobian(
      int,            ///< N (not sure what this is for! Must be internal)
      double,         ///< time, t
      const N_Vector, ///< current Y-value
      const N_Vector, ///< vector for Y-dot values
      const double *, ///< extra user-data vector, P, for evaluating ydot(y,t,p)
      DlsMat          ///< Jacobian matrix
      ) {cout <<"Jacobian not implemented in MPv10!\n"; return 1;}

  ///
  /// Get the number of extra parameters and the number of equations.
  ///
  virtual void get_problem_size(
      int *, ///< number of equations
      int *  ///< number of parameters in user_data vector.
      );

  protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  virtual void get_error_tolerances(
      double *, ///< relative error tolerance (single value)
      double []  ///< absolute error tolerance (array)
      );

  //---------------------------------------------------------------------------
  //-------------- END OF FUNCTIONS DERIVED FROM BASE CLASS -------------------
  //---------------------------------------------------------------------------

  //
  // Constant data in the cell, received from cell data.
  //
  double
    mpv_nH,     ///< total hydrogen number density at current cell.
    mpv_Vshell, ///< geometric factor in point source flux calculation (\sim 4\pi R^2 dR).
    mpv_Tau0,   ///< Optical depth of neutral hydrogen to front edge of cell at 13.6eV
    mpv_dTau0,  ///< Optical depth of neutral hydrogen through cell (=nH*dS*(1-x)*sigma) at 13.6eV.
    mpv_G0_UV,  ///< UV heating flux, including attenuation, F*exp(-1.9Av).
    mpv_G0_IR,  ///< Heating due to UV flux re-radiated in IR and re-absorbed, F*exp(-0.05Av).
    mpv_NIdot,  ///< photon luminosity of monochromatic ionising source (ionising photons/s).
    mpv_delta_S;///< path length through current cell.

};

#endif // MPv10_H




