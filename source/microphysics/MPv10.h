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
/// * This code is not yet working and should not be used.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */


#include "constants.h"
#include <algorithm>
//#include "hydrogen_mp.h"
//#include "hydrogen_photoion.h"


#include "cvode_integrator.h"
#include "microphysics_base.h"
#include "sim_params.h"
//#include "mellema_cooling.h"
#include "atomic_physics_database.h"
#include "chianti_cooling.h"
#include "coll_ionise_recomb.h"
//#include "photo_xsections.h"
#include <fstream>
#include <sys/stat.h>

/// this is the max change in x or E which is integrated with Euler integration.
/// A larger change will be integrated with backward-differencing and Newton
/// iteration, which is more accurate, stable, but expensive.
#define EULER_CUTOFF 0.05


#define MPv10_RELTOL                                                           \
  1.0e-04  ///< relative-error tolerance (actual error can be larger).
#define MPv10_ABSTOL 1.0e-12  ///< minimum neutral fraction I care.
#define MPv10_MINERG 1.0e-17  ///< Minimum internal energy density I care.

/// Assign an int to each species.
enum species_1 {
  NONE = -1,  ///< For if we want to specify no ion.
  H_0  = 0,   ///< neutral hydrogen.
  H_1p = 1,   ///< ionised hydrogen.
  He0  = 2,   ///< neutral helium.
  He1p = 3,   ///< singly ionised helium.
  He2p = 4,   ///< doubly ionised helium.
  C_0  = 5,   ///< Carbon
  C_1p = 6,
  C_2p = 7,
  C_3p = 8,
  C_4p = 9,
  C_5p = 10,
  C_6p = 11,
  N_0  = 12,  ///< Nitrogen
  N_1p = 13,
  N_2p = 14,
  N_3p = 15,
  N_4p = 16,
  N_5p = 17,
  N_6p = 18,
  N_7p = 19,
  O_0  = 20,  ///< Oxygen
  O_1p = 21,
  O_2p = 22,
  O_3p = 23,
  O_4p = 24,
  O_5p = 25,
  O_6p = 26,
  O_7p = 27,
  O_8p = 28
};


// establishes the ion_struct structure, which contains things like ionisation
// potential
struct ion_struct {
  std::string ion,          ///< Name of ion.
      ip1,                  ///< name of ion with higher ionisation stage.
      im1,                  ///< name of ion with lower ionisation stage.
      el;                   ///< Name of ion's element.
  struct ion_struct *iip1,  ///< pointer to higher stage ion, if present.
      *iim1;                ///< pointer to lower stage ion, if present.
  double ion_pot;    ///< Ionisation energy of ion to next ionisation stage, in
                     ///< ergs.
  int charge;        ///< Charge of ion in units of the electron charge.
  enum species_1 i;  ///< ion id number.
  int g_stat;        ///< statistical weight of this ion (H0=2,H1+=1,etc.)
};



//############################################################################
// MPv10 CLASS VARIABLES AND FUNCTIONS
///
/// Integrator for the microphysics equations for the non-equilibrium ion
/// fraction of various ions, and for the internal energy density.
///
class MPv10 :
    public microphysics_base,
    public cvode_solver,
    public coll_ionise_recomb,
    // public mellema_cooling
    public chianti_cooling
//    public photo_xsections
{
public:
  // ****** CONSTRUCTOR ******
  MPv10(
      const int,               ///< grid dimensions
      const int,               ///< Coordinate System flag
      const int,               ///< Total number of variables in state vector
      const int,               ///< Number of tracer variables in state vector.
      const std::string *,     ///< List of what the tracer variables mean.
      struct which_physics *,  ///< pointer to extra physics flags.
      struct rad_sources *,    ///< radiation sources.
      const double             ///< EOS Gamma
  );

  // ****** DESTRUCTOR ******
  ~MPv10();

  // ****** CATALOGER ******
  ///< MPv10 species tracer list.
  std::vector<vector<string> > MPv10_tracer_list;
  ///< MPv10 ion list.
  std::vector<string> MPv10_ion_list;

  ///< Record which element the ion belong to.
  std::vector<int> ions_tracer_elem;

  // Map of ion string and their primitive index
  std::map<string, int> ions_primindex_map;

  ///< Vector to record position of elements in the primitive vector.
  std::vector<int> elem_prim_index;
  ///< Vector for atomic masses of elements in the tracer list.
  std::vector<double> elem_atomic_mass;
  ///< Vector to record position of ions in the primitive vector.
  std::vector<int> ions_prim_index;
  ///< Vector to record position of ions in the local vector.
  std::vector<int> ions_local_index;
  ///< Vector to record number of ionic states in tracer elements.
  std::vector<int> N_sp_by_elem;
  ///< Vector to record number of species (neutral and ion )in tracer elements.
  std::vector<int> N_ions_by_elem;
  ///< Vector to record number of electron contributed by each ion.
  std::vector<int> ions_electron_num;
  ///< No of species in local y vector
  int N_sp;
  /// Record the lesser ion's local index if exist in MPv10_ion_list.
  std::vector<int> minus_ions_local_index;

  int N_elements;

  int N_prim;

  int N_local;

  int E_index;

  // Index of the first tracer variable in the primitive vector.
  int first_tracer_index;

  // ****** LOOKUP TABLE ******
  /// list of species collisional ionisation LUT
  std::vector<string> ci_tracer_list;
  /// Record which element the collisional ionisation tracer species belong to.
  std::vector<int> ci_tracer_elem;
  /// list of species recombination LUT
  std::vector<string> recomb_tracer_list;
  /// Record which element the recombination tracer species belong to.
  std::vector<int> recomb_tracer_elem;
  /// Uses this steps to generate lookup tables
  double delta_logT;
  double delta_logne;
  std::vector<double> T_table;   ///< Temperature table
  std::vector<double> ne_table;  ///< electron number density table


  const double T_min;
  const double T_max;
  const int Num_temps;  // NB this needs to be const so can initialise arrays
                        // with it. If this is >=1e4, get stack overflow errors.

  const int Num_ne;
  const double ne_min;
  const double ne_max;

  /// CIE cooling function (for CIE test)
  std::ofstream mc_outfile;

  // ****** GET RATES ******
  /// Function to return collisional ionisation rate of the given species for
  /// relevant temperature.
  double collisional_ionisation_rate(
      const int,    // Species identifier
      const double  // Temperature
  );
  /// Function to return recombination rate.
  double recombination_rate(
      const int,    // Species identifier
      const double  // Temperature
  );
  /// Fucntion to return mellema colling rate.
  virtual double mellema_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );
  /// Fucntion to return chianti colling rate.
  virtual double chianti_cooling_rate(
      const int,     //  species identifier
      const double,  // Temperature
      const double   // ne
  );


protected:
  /// Function will generate all the indexing necessary for MPv10 module.
  void tracer_species_cataloger(
      const int,  ///< Number of variables in the primitive vector.
      const int,  ///< Number of tracer variables in the primitive vector.
      const std::string *  ///< List of tracer variables.
  );

  // END OF MPv10 CLASS VARIABLES AND FUNCTIONS
  // ***************************************************************************



  //############################################################################
  // FUNCTIONS DERIVED FROM BASE CLASS FOLLOW
public:
  N_Vector y_in,    ///< current y-vector
      y_out;        ///< output y-vector
  int N_equations;  ///< number of equations i.e. number of species in Y.
  int N_extradata;  ///< number of elements in user-data array.

  ///
  /// calculate dy/dt for the vector of y-values.
  ///
  virtual int ydot(
      double,  ///< current time (probably not needed for rate equations)
      const N_Vector,  ///< current Y-value
      N_Vector,        ///< vector for Y-dot values
      const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
  );

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  ///
  virtual int Jacobian(
#if defined CVODE2
      int,  ///< N (not sure what this is for!)
#endif
      double,          ///< time, t
      const N_Vector,  ///< current Y-value
      const N_Vector,  ///< vector for Y-dot values
      const double
          *,    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
      CVMatrix  ///< Jacobian matrix
  )
  {
    spdlog::info("Jacobian not implemented in MPv10!");
    return 1;
  }

  ///
  /// Get the number of extra parameters and the number of equations.
  ///
  virtual void get_problem_size(
      int *,  ///< number of equations
      int *   ///< number of parameters in user_data vector.
  );


protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  virtual void get_error_tolerances(
      double *,  ///< relative error tolerance (single value)
      double[]   ///< absolute error tolerance (array)
  );
  // END OF FUNCTIONS DERIVED FROM BASE CLASS
  // ***************************************************************************



  //############################################################################
  //############################################################################


  ///
  /// Function for updating vectors according to species found in tracer list.
  ///
  // s pecies_tracer_initialise

  ///
  /// A primitive vector is input, and lots of optional extra data in
  /// the last argument, and the ion fraction(s) and internal energy
  /// are updated by the requested timestep.  Results are saved to
  /// the destination vector, which can be the same pointer as the
  /// input vector.
  ///
  int TimeUpdateMP(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,         ///< Switch for what type of integration to use.
      double *           ///< Vector of extra data (column densities, etc.).
  );

  ///
  ///
  ///
  int TimeUpdate_RTsinglesrc(
      const pion_flt *,  ///< Primitive Vector to be updated.
      pion_flt *,        ///< Destination Vector for updated values.
      const double,      ///< Time Step to advance by.
      const double,      ///< EOS gamma.
      const int,         ///< Switch for what type of integration to use.
      const double,      ///< flux in per unit length along ray (F/ds or L/dV)
      const double,      ///< path length ds through cell.
      const double,      ///< Optical depth to entry point of ray into cell.
      double *  ///< return optical depth through cell in this variable.
  )
  {
    spdlog::info("MPv10::TimeUpdate_RTsinglesrc() not implemented");
    return 1;
  }

  ///
  /// Update microphysics including photoionization and photoheating
  ///
  virtual int TimeUpdateMP_RTnew(
      const pion_flt *,  ///< Primitive Vector to be updated.
      const int,         ///< unused
      const std::vector<struct rt_source_data> &,  ///< unused
      const int,                                   ///< unused
      std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      pion_flt *,    ///< Destination Vector for updated values
      const double,  ///< Time Step to advance by.
      const double,  ///< EOS gamma.
      const int,     ///< unused
      double *       ///< any returned data (final temperature?).
  );

  ///
  /// Set the properties of a multifrequency ionising radiation source.
  ///
  int set_multifreq_source_properties(
      // int tempin,
      // double *
      const struct rad_src_info *,  ///< source data
      double *  ///< O/P source strength in different energy bins.
  );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  /// - Is threadsafe.
  ///
  double Temperature(
      const pion_flt *,  ///< primitive vector
      const double       ///< eos gamma
  );

  ///
  /// Set the gas temperature to a specified value.
  /// Needed if you want this feature to set initial conditions.
  ///
  int Set_Temp(
      pion_flt *,    ///< primitive vector.
      const double,  ///< temperature
      const double   ///< eos gamma.
  );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.
  ///
  virtual double timescales(
      const pion_flt *,  ///< Current cell.
      const double,      ///< EOS gamma.
      const bool,        ///< set to 'true' if including cooling time.
      const bool,        ///< set to 'true' if including recombination time.
      const bool         ///< set to 'true' if including photo-ionsation time.
  );

  ///
  /// This returns the minimum timescale of all microphysical
  /// processes, including reaction times for each species and the
  /// total heating/cooling time for the gas.
  /// It requires the radiation field as an input, but this is not
  /// used for now so the vectors can be null.
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
  );

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium
  /// value. This is optionally used in the initial condition
  /// generator.  Not implemented here.
  ///
  int Init_ionfractions(
      pion_flt *,    ///< Primitive vector to be updated.
      const double,  ///< eos gamma.
      const double   ///< optional gas temperature to end up at.
  )
  {
    spdlog::info("MPv10::Init_ionfractions() not implemented");
    return 1;
  }

  ///
  /// Return index of tracer for a given string.
  ///
  int Tr(const string  ///< name of tracer we are looking for.
  );

  ///
  /// Get the total cooling rate.  This is for postprocessing the
  /// simulation data only -- IT IS NOT OPTIMISED FOR SPEED.
  ///
  virtual double total_cooling_rate(
      const pion_flt *,  ///< Current cell values.
      const int,         ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &,
      ///< list of UV-heating column densities and source properties.
      const int,  ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &,
      ///< list of ionising src column densities and source properties.
      const double  ///< EOS gamma.
  )
  {
    spdlog::info("MPv10::total_cooling_rate() not implemented");
    return 1;
  }

  ///
  /// Get the total recombination rate for an ion, given the input
  /// state vector.
  ///
  virtual double get_recombination_rate(
      const int,         ///< ion index in tracer array (optional).
      const pion_flt *,  ///< input state vector (primitive).
      const double       ///< EOS gamma (optional)
  )
  {
    spdlog::info("MPv10::get_recombination_rate() not implemented");
    return 1;
  }

  ///
  /// Return the H mass fraction
  ///
  virtual double get_X_H();

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_elec(const pion_flt *  ///< primitive state vector.
  );

  ///
  /// Get electron number density (cm^{-3})
  ///
  double get_n_ion(
      std::string,      ///< ion name
      const pion_flt *  ///< primitive state vector.
  );

  // ###########################################################################
  // ###########################################################################

protected:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  ///
  virtual int convert_prim2local(
      const pion_flt *,       ///< primitive vector (nv_prim)
      std::vector<double> &,  ///< local vector [x(H0),E](n+1).
      int                     /// flags the function
  );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  virtual int convert_local2prim(
      const std::vector<double> &,  ///< local (updated) vector [x(H0),E](n+1).
      const pion_flt *,             ///< input primitive vector (nv_prim)
      pion_flt *,                   ///< updated primitive vector (nv_prim)
      int                           /// flags the function
  );

  ///
  /// Set local variables (y_ion_frac, y_neutral_frac, ne) by calling
  /// set_y_ion_frac(), set_y_neutral_frac), set_ne().
  ///
  virtual void set_localvariables(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &,        ///< y_ion_frac
      std::vector<double> &,        ///< y_neutral_frac
      double &                      ///< ne
  );

  /// Set y_ion_frac
  virtual void set_y_ion_frac(
      const std::vector<double> &,  ///< local vector
      std::vector<double> &         ///< y_ion_frac
  );

  /// Set y_neutral_frac
  virtual void set_y_neutral_frac(
      const std::vector<double> &,  ///< y_ion_frac
      std::vector<double> &         ///< y_neutral_frac
  );

  /// Calculate electron number density
  virtual void set_ne(
      const std::vector<double> &,  ///< y_ion_frac
      double &                      ///< electron number density
  );

  /// Correct local vector
  /// + check that ion fractions are >= ABSTOL.
  /// + check that energy is within temperature limits
  virtual int correct_localvector(
      std::vector<double> &,  /// local vector
      int                     /// flags the function
  );


  // ###########################################################################
  // ###########################################################################



  ///
  /// returns gas temperature according to E=nkT/(g-1) with n=1.1*nH*(1+x_in),
  /// appropriate for a gas with 10% Helium by number, and if He is singly
  /// ionised whenever H is.
  ///
  virtual double get_temperature(
      const std::vector<double> &,    ///< ion fractions
      const std::vector<pion_flt> &,  ///< element number densities
      const double                    ///< E_int (per unit volume)
  );

  ///
  /// Returns total number of particles.
  ///
  virtual double get_ntot(
      const std::vector<double> &,  ///< ion fractions (input)
      const std::vector<double> &   ///< element number density (input)
  );

  ///
  /// Set the size of local vectors, and index them with integers.
  ///
  virtual void setup_local_vectors();


  ///
  /// Updates the corrector vector
  /// according to sCMA (simple Consistent Multi-fluid Advection, Plewa +
  /// Muller, 1999). Used for modifying tracer fluxes.
  ///
  virtual void sCMA(
      std::vector<double>
          &,            ///< input corrector vector (changed by this function)
      const pion_flt *  ///< input primitive vector (nv_prim)
  );



  ///
  /// Calculates the change in optical depth over the cell
  /// for use in the Raytracing module
  ///
  virtual void get_dtau(
      const pion_flt,  ///< ds, thickness of the cell
      const pion_flt
          *,      ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *  ///< output dtau vector
  );



  //
  // Any data which can be calculated/set at the start of the simulation
  // can be defined here.
  //


  // double erg_per_eV;       /// convert eV to erg.

  const int ndim;          ///< Number of dimensions in grid.
  const double eos_gamma;  ///< EOS gamma for ideal gas.
  const int coord_sys;     ///< Coordinate System flag
  double gamma_minus_one;  ///< as named.
  double Min_NeutralFrac;  ///< minimum H0 fraction allowed (eps=1.0e-12)
  double Max_NeutralFrac;  ///< Maximum H0 fraction allowed (1-eps)

  // int nvl;      ///< number of variables in local state vector.
  // int lv_eint;  ///< internal energy local variable index.
  // int ftr;  // lv_y_ion_index_offset; ///< gives the index at which ions
  // first occur in primitive vector, maps to first index of local vector



  ///***********************************************

  ///	Vector to store initial corrector values (used to modify flux according to
  /// sCMA)
  std::vector<double> corrector;


  int print_flag;


  std::vector<double> elem_number_density;



  /// ===========================================================================
  ///  unknown variables
  /// ===========================================================================


  int N_diff_srcs,   ///< No diffuse sources --> 0, otherwise --> 1
      N_ion_srcs,    ///< No ionising sources --> 0, otherwise --> 1
      ion_src_type;  ///< Either RT_EFFECT_MFION or RT_EFFECT_PION_MONO.

  /// For each cell, ydot() needs to know the local radiation field.  This
  /// function takes all of the input radiation sources and calculates what
  /// ydot() needs to know for both heating and ionisation sources.
  ///
  void setup_radiation_source_parameters(
      const pion_flt *,  ///< primitive input state vector.
      vector<double> &,  ///< local input state vector (x_in,E_int)
      std::vector<struct rt_source_data> &
      ///< list of ionising src column densities and source properties.
  );

  // Constant data in the cell, received from cell data.
  double mpv_nH;  ///< total hydrogen number density at current cell.
  int N_rad_src;  ///< number of radiation sources to include.
  std::vector<struct rt_source_data> rt_data;  ///< data on radiation sources.

  /*double Emax[15] = {13.6*1e-3, 14.5*1e-3, 24.4*1e-3,
  24.58741*1e-3, 29.6*1e-3, 47.5*1e-3, 47.9*1e-3, 54.41778*1e-3,
  64.5*1e-3, 77.5*1e-3, 97.9*1e-3, 392.1*1e-3, 490.0*1e-3, 552.1*1e-3,
  667.0*1e-3};//bin edges correspond to ionisation energies double Emin[15] =
  {11.3*1e-3, 13.6*1e-3, 14.5*1e-3, 24.4*1e-3, 24.58741*1e-3,
  29.6*1e-3, 47.5*1e-3, 47.9*1e-3, 54.41778*1e-3, 64.5*1e-3,
  77.5*1e-3, 97.9*1e-3, 392.1*1e-3, 490.0*1e-3, 552.1*1e-3}; int Nbins = 15;*/


  /// index matching photo_ion xsection to ion in local vector
  std::vector<std::vector<double> > y_ion_xsections;



  // #########################################################################
  // *********************** DEBUG  FUNCTIONS ********************************
  // #########################################################################

  /// This will print collisional ionisation or recombination rate for a
  /// different set of T points in the range T_min<T<T_max using already
  /// setup rate and slope lookup table.
  void print_CIR_LUT(
      const std::vector<double> &,  // T table
      // rate table
      const std::vector<std::vector<double> > &,
      // rate slope table
      const std::vector<std::vector<double> > &,
      // which rate (collisional or recombination )
      const std::string filename);

  /// Print 2D vector to a file
  void print_2D_vector(
      std::vector<std::vector<double> > &,  // Vector to be written
      const std::string,                    // Name of the ion
      std::string                           // Base name/file path
  );

  /// Write 1D vector to a file
  void print_1D_vector(
      std::vector<double> &,  // Vector to be written
      const std::string,      // Name of the variable
      std::string             // Base name/file path
  );
};
#endif  // MPv10_H



/** \section numfrac Number Fractions
 * Number fractions are from Lodders,K., 2003, ApJ, 591, 1220-1247
 *  \section Masses
 * Masses are from Wikipedia, in Atomic Mass Units (amu), converted to units
 * of proton mass.
 * */
