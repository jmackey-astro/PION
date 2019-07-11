///
/// \file MPv10.h
/// \author Maggie Goulden
/// \date 2018.10
///
/// modifications:
/// - 2018.10 .

#ifndef MPv10_H
#define MPv10_H

/// Description:
/// A multi-ion chemistry solver with non-equilibrium ionization and
/// recombination of different species, and with elemental abundances
/// that are input for each solve (i.e. can vary within simulation).
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
#include <map>
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


#define MPv10_RELTOL 1.0e-4   ///< relative-error tolerance (actual error can be larger).
#define MPv10_ABSTOL 1.0e-12  ///< minimum neutral fraction i care about.
#define MPv10_MINERG 1.0e-17  ///< Minimum internal energy density I care about.


///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of various ions, and for the internal energy density.
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

  ///MPv10::MPv10
  /// Destructor
  ///
  ~MPv10();
  
  ///
  /// Function for updating vectors according to species found in tracer list.
  ///
  void species_tracer_initialise(
      const std::string *,  ///< List of what the tracer variables mean.
      int, ///< index of current tracer in for loop (i in s=tracers[i])
      std::string , /// < current tracer in for loop
      std::string, ///< element symbol, e.g. "He", "H"
      int, ///< element length, e.g. "H" is of length 1, "He" of length 2.
      int, ///< element index in N_species_by_elem, used in for loops for densities etc
      int ///< length of tracers vector
      );

  ///
  /// The NON-RT MICROPHYSICS update function.
  ///
  /// A primitive vector is input, and lots of optional extra data in
  /// the last argument, and the ion fraction(s) and internal energy
  /// are updated by the requested timestep.  Results are saved to
  /// the destination vector, which can be the same pointer as the
  /// input vector.
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
  {cout <<"MPv10::TimeUpdate_RTsinglesrc() not implemented\n";return 1;}

  ///
  /// Not used for this class, so far
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
      const double,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int, ///< Switch for what type of integration to use.
      double *    ///< any returned data (final temperature?).
      )
  {cout <<"MPv10::TimeUpdateMP_RTnew() not implemented\n";return 1;}

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
  /// Needed if you want this feature to set initial conditions.
  ///
  int Set_Temp(
      pion_flt *,     ///< primitive vector.
      const double, ///< temperature
      const double  ///< eos gamma.
      );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments. 
  ///
  virtual double timescales(
      const pion_flt *, ///< Current cell.
      const double,   ///< EOS gamma.
      const bool, ///< set to 'true' if including cooling time.
      const bool, ///< set to 'true' if including recombination time.
      const bool  ///< set to 'true' if including photo-ionsation time.
      );

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.
  /// It requires the radiation field as an input, but this is not
  /// used for now so the vectors can be null.
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
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value. This is optionally used in the initial condition
  /// generator.  Not implemented here.
  ///
  int Init_ionfractions(
      pion_flt *, ///< Primitive vector to be updated.
      const double, ///< eos gamma.
      const double  ///< optional gas temperature to end up at.
      )
  {cout <<"MPv10::Init_ionfractions() not implemented\n";return 1;}

  ///
  /// Return index of tracer for a given string.
  ///
  int Tr(
      const string ///< name of tracer we are looking for.
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
      )
  {cout <<"MPv10::total_cooling_rate() not implemented\n";return 1;}

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  virtual double get_recombination_rate(
      const int,      ///< ion index in tracer array (optional).
      const pion_flt *, ///< input state vector (primitive).
      const double    ///< EOS gamma (optional)
      )
  {cout <<"MPv10::get_recombination_rate() not implemented\n";return 1;}
  
  ///
  /// Return the H mass fraction
  ///
  virtual inline double get_X_H()
    {return EP->H_MassFrac;}

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_elec(
      const pion_flt * ///< primitive state vector.
      );

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_ion(
      std::string, ///< ion name
      const pion_flt * ///< primitive state vector.
      );

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
      double *,//const double, ///< nH //< y_ion_fraction (by y_ion_index)
      vector<pion_flt>&, ///< y_ion_number_density
      const double  ///< E_int (per unit volume)
      );

  ///
  /// Returns total number of particles.
  ///
  virtual double get_ntot(
      double *,//const double, ///< nH //< y_ion_fraction (by y_ion_index)
      vector<pion_flt>& ///< y_ion_number_density
    );

  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  virtual void setup_local_vectors();


	/// 
	/// Updates the corrector vector
	/// according to sCMA (simple Consistent Multi-fluid Advection, Plewa + Muller, 1999).
	/// Used for modifying tracer fluxes.
	///
	virtual void sCMA(
			std::vector<double>, ///< input corrector vector
			const pion_flt *); ///< input primitive vector from grid cell (length nv_prim)

  //
  // ********* STUFF FROM THE mp_v2_aifa CLASS **********
  //
  N_Vector 
    y_in,  ///< current y-vector
    y_out; ///< output y-vector
  int N_equations; ///< number of equations i.e. number of species in Y.
  int N_extradata; ///< number of elements in user-data array.


  //
  // Any data which can be calculated/set at the start of the simulation
  // can be defined here.
  //
  double k_B; ///<  Boltzmanns constant.
  double m_p; ///< Mass of proton.
  double m_H; ///< Mass of hydrogen (grams).
  double m_He; ///< Mass of helium.
  double m_C; ///< Mass of Carbon
  double m_N; ///< Mass of Nitrogen
  double m_O; ///< Mass of Oxygen
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

  int nvl;     ///< number of variables in local state vector.
  int lv_eint; ///< internal energy local variable index.
  int ftr; //lv_y_ion_index_offset; ///< gives the index at which ions first occur in primitive vector, maps to first index of local vector
  
  int N_elem;
  int N_species;
  int N_eqns; ///< := n species
  int N_cons_eqns; /// < := n_elem + 1, due to E_int.
  
	///	Vector to store initial corrector values (used to modify flux according to sCMA)
	std::vector<double> corrector;
  std::vector<double> init_corrector;
  /// ===========================================================================
  ///               Vectors to Access Primitive / Local vectors
  /// ===========================================================================
  std::vector<int> X_mass_frac_index; /// < primitive vector indices, used to trace X_H etc, like pv_Hp. 
  std::vector<int> y_ion_index_prim; ///<index matching y_ion mass fraction in prim vector, analogous to pv_Hp before.
  std::vector<int> y_ion_index_local; ///<index matching y_ion fraction in local vector.
  
  int H_ion_index; ///index of X_elem in primitive
  int He_ion_index;
  int C_index;
  int N_index;
  int O_index;
  std::map<string,int> tracer_list;
  std::map<string,int> element_list;

	int print_flag;

  /// ===========================================================================
  ///               Vectors to Access Adjacent Ions
  /// ===========================================================================
  std::vector<int> y_ip1_index_local; ///< local index of ion p1 (plus 1); ion one stage higher. If ion doesn't exist, sets to -1.
  std::vector<int> y_im1_index_local; ///< local index of ion m1 (minus 1); ion one stage lower. if lower stage ion is neutral species, sets to -2.
                                      ///   also, if we're just not tracking lower stage (i.e. "doesn't exist"), also sets to -1.
  std::vector<int> y_ip1_index_tables; ///< index of the ip1 species in tables. This follows the order "H0", "H1+", "He0", "He1+", "He2+" etc.
  std::vector<int> y_ion_index_tables;
  std::vector<int> y_im1_index_tables; ///< index of the im1 species in tables. Follows same order as above.
  
  /// ===========================================================================
  ///           Vectors to Store Other Ion / Element info
  /// ===========================================================================
  std::vector<int> y_ion_num_elec; ///<index matching number of electrons corresponding to each y_ion (e.g. C6+ = 6)
  std::vector<double> X_elem_atomic_mass; /// < vector of atomic masses corresponding to X_elem_index, e.g. for ["H", "He"], X_atom_mass=[1.6738e-24, 6.6464764e-24 
  std::vector<double> X_elem_number_density;
  std::vector<int> N_species_by_elem; ///< records # species in each element, to iterate over later.

  /// ===========================================================================
  ///  Variables to store table (temp, ionisation / recomb rates, etc)information
  /// ===========================================================================
  std::vector<double> Temp_Table;
  std::vector< std::vector<double> > recomb_rate_table;
  std::vector< std::vector<double> > ionise_rate_table;
  std::vector< std::vector<double> > recomb_slope_table;
  std::vector< std::vector<double> > ionise_slope_table;
  
  std::vector<double> ionisation_potentials; //array of ionisation potentials, in the same order as ionise_slope
  
  /// Assign an int to each species.
  enum species {
    NONE =-1, ///< For if we want to specify no ion.
    H_0  = 0, ///< neutral hydrogen.
    H_1p = 1, ///< ionised hydrogen.
    He0  = 2, ///< neutral helium.
    He1p = 3, ///< singly ionised helium.
    He2p = 4, ///< doubly ionised helium.
    C0  = 5, ///< Carbon
    C1p = 6,
    C2p = 7,
    C3p = 8,
    C4p = 9,
    C5p = 10,
    C6p = 11,
    N0  = 12, ///< Nitrogen
    N1p = 13,
    N2p = 14,
    N3p = 15,
    N4p = 16,
    N5p = 17,
    N6p = 18,
    N7p = 19,
    O0  = 20, ///< Oxygen
    O1p = 21,
    O2p = 22,
    O3p = 23,
    O4p = 24,
    O5p = 25,
    O6p = 26,
    O7p = 27,
    O8p = 28
  };

  
  //establishes the ion_struct structure, which contains things like ionisation potential
  struct ion_struct {
  std::string ion, ///< Name of ion.
    ip1, ///< name of ion with higher ionisation stage.
    im1, ///< name of ion with lower ionisation stage.
    el;  ///< Name of ion's element.
    struct ion_struct 
      *iip1, ///< pointer to higher stage ion, if present.
      *iim1; ///< pointer to lower stage ion, if present.
  double ion_pot; ///< Ionisation energy of ion to next ionisation stage, in ergs.
  int charge; ///< Charge of ion in units of the electron charge.
  enum species i; ///< ion id number.
  int g_stat; ///< statistical weight of this ion (H0=2,H1+=1,etc.)
  };


  std::map<std::string,struct ion_struct> ion_props; ///< properties of each ion.
  struct ion_struct **ii; ///< pointer to array of ions used.
  
  void set_atomic_data();
  double Coll_Ion_rate(
       double,   ///< Precalculated Temperature.
       const ion_struct *///< current ion.
       );
  double Rad_Recomb_rate(
      double,   ///< Precalculated Temperature.
      const struct ion_struct * ///< current ion.
      );
  
  double dielec_recomb(
      double,
      enum species
      );
  double rad_recomb(
      double,
      enum species
      );
  void generate_lookup_tables();
  
  const double T_min;
  const double T_max;
  const int Num_temps; //NB this needs to be const so can initialise arrays with it. If this is >=1e4, get stack overflow errors.
  double delta_log_temp;
  
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




