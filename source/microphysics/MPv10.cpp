/// \file MPv10.cpp
/// \author Maggie Celeste Goulden
/// \date 2018.10
///
/// Description:
/// - multi-species non-equilibrium chemistry solver.

/// ionization/recombination
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method use backwards differencing (i.e. implicit) with Newton
/// iteration.
///
/// Modifications:
/// - 2018.10.09 JM: edited header
///



// ##################################################################
// List of Abbreviations used.
// LUT - Look Up Table
// CLV  - Correct Local Vector
// CIR - collisional ionisation and recombination
// ******************************************************************



// ##################################################################
// Note:
// When referring to element, we use variable 'elem' in loops
// When referring to ions, we use 'i' in loops
// When referring to neutral atoms + ions, we use 'species' and
// variable 's' in loops
// ******************************************************************



#include "microphysics/MPv10.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG



using namespace std;

//#define MPv10_CIE_TEST
//#define MPv10_CIE_COOLING
//#define MPv10_DEBUG

// Timestep-limiting is important for making chemistry consistent
// with hydrodynamics.
// This is a good value for MPv3 (see Mackey,2012,A&A,539,A147)
// Will have to do some tests for MPv10...

#define DTFRAC 0.25


// GET ERROR TOLERANCES #################################################
void MPv10::get_error_tolerances(
    double *reltol,  ///< relative error tolerance.
    double atol[]    ///< absolute error tolerances
)
{
  *reltol = MPv10_RELTOL;
  for (int i = 0; i < N_equations - 1; i++)
    atol[i] = MPv10_ABSTOL;  ///< minimum neutral fraction I care about.
  atol[N_equations - 1] = MPv10_MINERG;  ///< E_int: for n=1.0, T=1.0e4, ==>
                                         ///< E=2.07e-12, so say 1e-17?
  return;
}
// END OF GET ERROR TOLERANCES



// GET PROBLEM SIZE  ########################################
void MPv10::get_problem_size(
    int *n_eqn,  ///< number of equations
    int *n_para  ///< number of parameters in user_data vector.
)
{
  *n_eqn  = N_equations;
  *n_para = N_extradata;
  return;
}
// END OF GET PROBLEM SIZE  ########################################



// CONSTRUCTOR ############################################################
MPv10::MPv10(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr,   ///< Number of tracer variables in state vector.
    const std::string *tracers,   ///< List of what the tracer variables mean.
    struct which_physics *ephys,  ///< pointer to extra physics flags.
    struct rad_sources *rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    microphysics_base(nv, ntr, tracers, ephys, rsrcs),
    coll_ionise_recomb(),
    // mellema_cooling(),
    chianti_cooling(), ndim(nd), eos_gamma(g), coord_sys(csys), T_min(1e1),
    T_max(1e9), Num_temps(300), Num_ne(100), ne_min(1.0), ne_max(1.0e6)

//  , photo_xsections()  // photo_xsections(&Emin[0],&Emax[0],Nbins)
{

  spdlog::info("Setting up microphysics MPv10 module");
  spdlog::debug(
      "MPv10: Setting up tracer variables assuming tracers are last {} "
      "variables in primitive vector",
      ntracer);
  spdlog::debug(
      "MPv10: Tracers variables {} ",
      std::vector<string>(tracers, tracers + ntr));

  // IONISATION AND RECOMBINATION **************************************
  // This function will generate all necessary indices for MPv10 module.
  tracer_species_cataloger(nv, ntr, tracers);

  spdlog::debug("MPv10: Total {} variables in the primitive vector", N_prim);
  spdlog::debug("MPv10: MPv10 tracer list {}", MPv10_tracer_list);
  spdlog::debug("MPv10: Element's primitive index {}", elem_prim_index);
  spdlog::debug("MPv10: Element's atomic mass {}", elem_atomic_mass);
  spdlog::debug("MPv10: MPv10 ion list {}", MPv10_ion_list);
  spdlog::debug("MPv10: ion's element index {}", ions_tracer_elem);
  spdlog::debug(
      "MPv10: ion's string - primitive index map {}", ions_primindex_map);
  spdlog::debug("MPv10: ions primitive index {}", ions_prim_index);
  spdlog::debug("MPv10: ions local index {}", ions_local_index);
  spdlog::debug("MPv10: Number of species in local vector {}", N_sp);
  spdlog::debug("MPv10: Number of species by elements {}", N_sp_by_elem);
  spdlog::debug("MPv10: Number of ions by elements {}", N_ions_by_elem);
  spdlog::debug("MPv10: ion's electron numbers {}", ions_electron_num);
  spdlog::debug("MPv10: minus ions local index {}", minus_ions_local_index);
  spdlog::debug("MPv10: ci tracer list = {}", ci_tracer_list);
  spdlog::debug("MPv10: ci tracer elements = {}", ci_tracer_elem);
  spdlog::debug("MPv10: recomb tracer list = {}", recomb_tracer_list);
  spdlog::debug("MPv10: recomb tracer elements = {}", recomb_tracer_elem);


  // lOOKUP TABLE SETTINGS *********************************************
  // Setting up step size for log temperature (T) and log electron number
  // density (ne). This step size is used to generate look up tables. The same
  // is used to calculate temperature index and electron number density index to
  // get the rate for a specific value of T or ne. Setting up uniform step size
  // for log T.
  delta_logT = (log10(T_max) - log10(T_min)) / (Num_temps - 1);
  // Make temperature table for all LUTs.
  T_table.resize(Num_temps);
  for (int i = 0; i < Num_temps; i++) {
    // Temperature table for coll_ionise and recomb LUTs:
    T_table[i] = pow(10.0, (log10(T_min) + i * delta_logT));
  }
  // Setting up uniform step size for log ne.
  delta_logne = (log10(ne_max) - log10(ne_min)) / (Num_ne - 1);
  // Make electron number density table for all LUTs.
  ne_table.resize(Num_ne);
  for (size_t i = 0; i < Num_ne; i++) {
    ne_table[i] = pow(10, log10(ne_min) + i * delta_logne);
  }

  // IONISATION AND RECOMBINATION *********************************************
  // setting up collisional ionisation
  setup_collisional_ionisation(ci_tracer_list);
  // Generating collisional ionisation rate lookup tables
  generate_collisional_ionisation_LUT(T_table);
  // setting up recombination
  setup_recombination(recomb_tracer_list);
  // Generating recombination rate lookup tables
  generate_recombination_LUT(T_table);

#ifdef MPv10_DEBUG
  mkdir("MPV10_coll-recom_tables", 0777);
  // PRINT TO FILE COLLISIONAL AND RECOMBINATION TABLES ***********************
  print_CIR_LUT(
      T_table, collisional_rate_table, collisional_slope_table,
      "coll-recom_tables/coll-ionise-rate");
  print_CIR_LUT(
      T_table, recombination_rate_table, recombination_slope_table,
      "coll-recom_tables/recomb-rate");
#endif

  /*
  // MELLEMA COOLING FUNCTION *************************************************
  setup_mellema_cooling(MPv10_tracer_list);
  spdlog::debug("mc species locator = {}", mc_tracer_locator);
  generate_mellema_table(T_table, ne_table);
  */

  // CHIANTI COOLING FUNCTION *************************************************
  setup_chianti_cooling(MPv10_tracer_list);
  spdlog::debug("chianti species locator = {}", chianti_tracer_locator);
  generate_chianti_table(T_table, ne_table);

#ifdef MPv10_DEBUG
  mkdir("MPV10_cooling_tables", 0777);
  // PRINT MELLEMA COOLING TABLES ********************************************
  print_1D_vector(T_table, "logT", "cooling_tables/chianti_rate");
  print_1D_vector(ne_table, "logne", "cooling_tables/chianti_rate");

  // print resized chianti cooling rate table to files
  // print the log10 of actual rate
  spdlog::debug("Writting mellema cooling look up tables into files");
  int chianti_species_index = 0;
  int location;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      location = chianti_tracer_locator[i][j];  // database location
      spdlog::debug(
          "Writting LUT({}) as Mellema-{}", MPv10_tracer_list[i][j],
          chianti_table[chianti_species_index].Name);
      print_2D_vector(
          chianti_table[chianti_species_index].rate, MPv10_tracer_list[i][j],
          "cooling_tables/chianti_rates");
      chianti_species_index++;
    }
  }

  /*
  // PRINT MELLEMA COOLING TABLES ********************************************
  //print log10 values.
  print_1D_vector(T_table, "logT", "cooling_tables/mellema_rate");
  print_1D_vector(ne_table, "logne", "cooling_tables/mellema_rate");

  // print resized mellema cooling rate table to files
  // print the log10 of actual rate
  spdlog::debug("Writting mellema cooling look up tables into files");
  int mc_species_index = 0;
  int location;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      location = mc_tracer_locator[i][j];  // database location
      spdlog::debug("Writting LUT({}) as Mellema-{}",
                    MPv10_tracer_list[i][j],
  mellema_table[mc_species_index].Name); print_2D_vector(
  mellema_table[mc_species_index].rate, MPv10_tracer_list[i][j],
  "cooling_tables/mellema_rates"); mc_species_index++;
    }
  }
  */
#endif


#ifdef MPv10_CIE_COOLING
  // create a text file to write the cooling as the function of temperature.
  /// (remove this part when it is done)
  std::ofstream mc_outfile("cooling_function.txt");
  mc_outfile
      << "# Collisional ionisation equilibrium (CIE) cooling function for \n";
  mc_outfile << "# all species in the MPv10 tracer list.\n";
  mc_outfile
      << "# Column-1: Temperature (K), Other columns: Cooling function (erg cm^3/s).\n";
#endif

  // INITIALIZING IMPORTANT VARIABLES AND ARRAY ************************
  // resizing element number density vector
  elem_number_density.resize(N_elements);
  // resizing the corrector vector
  corrector.resize(N_prim, 1.0);
  // --- Set up local variables
  setup_local_vectors();
  gamma_minus_one = eos_gamma - 1.0;
  Min_NeutralFrac = MPv10_ABSTOL;
  Max_NeutralFrac = 1.0 - MPv10_ABSTOL;


  // CVODE SOLVER *****************************************************
  // Initialise the CVODES solver memory etc.
  setup_cvode_solver_without_Jacobian();


  // Set flags for whether we have radiation sources *******************
  N_rad_src = 0;
  for (int isrc = 0; isrc < RS->Nsources; isrc++) {
    if (RS->sources[isrc].type == RT_SRC_SINGLE
        && RS->sources[isrc].effect == RT_EFFECT_MFION) {
      N_rad_src++;
      rt_data.resize(N_rad_src);
      spdlog::debug("RT_data size {}", rt_data.size());
      int err = set_multifreq_source_properties(
          &RS->sources[isrc], rt_data[N_rad_src - 1].strength);
      if (err) {
        spdlog::error("{}: {}", "multifreq photoionisation setup MPv10", err);
        exit(1);
      }
    }
  }
  spdlog::debug("MPv10: got {} radiation sources", N_rad_src);

  // CONSTRUCTOR END MESSAGE  *****************************************
  spdlog::debug("MPv10: Constructor finished and returning");
}
// END OF CONSTRUCTOR ##########################################################



// ##################################################################

void MPv10::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  N_local     = N_sp + 1;  // number of species plus energy
  N_extradata = 0;
  N_equations = N_local;
#if defined(CVODE6)
  y_in  = N_VNew_Serial(N_equations, sunctx);
  y_out = N_VNew_Serial(N_equations, sunctx);
#else
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
#endif
  E_index = N_sp;
  return;
}
// ##################################################################



// MPV10 DESTRUCTOR #########################################################
MPv10::~MPv10()
{
  // Free vector memory
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}
// END OF MPV10 DESTRUCTOR ##################################################



// TRACER SPECIES CATALOGER ################################################
void MPv10::tracer_species_cataloger(
    const int N_primitive, const int N_tracers, const std::string *tracers)
{

  // Beginning of cataloging tracers
  spdlog::info("MPv10: Cataloging microphysics tracers");

  // Make a copy of tracer list
  string tracer_list_replica[N_tracers];
  for (int s = 0; s < N_tracers; s++)
    tracer_list_replica[s] = tracers[s];

  // Number of elements in the primitive vector.
  N_prim = N_primitive;

  // Index of the first tracer variable in the primitive vector.
  first_tracer_index = N_prim - N_tracers;


  // Creating related indexing of neutral species in the tracer ***************
  N_elements = 0;
  // Tracer loop for identifying neutral species.
  for (int s = 0; s < N_tracers; s++) {
    if (tracer_list_replica[s].substr(0, 2) == "X_") {
      // Index vector to record position of elements in the primitive vector.
      elem_prim_index.push_back(first_tracer_index + N_elements);
      // Enumerate number of elements in the tracer.
      N_elements++;
      // Rename neutral ions by removing "X_".
      tracer_list_replica[s].erase(0, 2);
      // Make MPv10 tracer species list
      MPv10_tracer_list.resize(N_elements);
      MPv10_tracer_list[s].push_back(tracer_list_replica[s]);
      // Atomic masses of elements in the tracer list.
      elem_atomic_mass.push_back(get_atomic_mass(tracer_list_replica[s]));
      // No of ions in each element
      N_sp_by_elem.resize(N_elements, 0);  // rename this ***
      // No of species in each element
      N_ions_by_elem.resize(N_elements, 1);  // rename this ***
    }
  }
  // End of tracer loop for identifying neutral specie*************************


  // Creating related indexing of charged species in the tracer list **********
  int ion_count = 0;
  int delete_length;
  int string_length;
  // Note: The tracer loop is restarted at s = N_elements.
  // Tracer loop for identifying ion species.
  for (int s = N_elements; s < N_tracers; s++) {

    // Extract the element symbol of the ion
    delete_length = 1;
    string_length = tracer_list_replica[s].size();
    for (int i = 0; i < string_length; i++) {
      if (isdigit(tracer_list_replica[s][i])) delete_length++;
    }
    std::string element_symbol = tracer_list_replica[s].substr(
        0, tracer_list_replica[s].size() - delete_length);

    // Find which element the ion belong to.
    for (int e = 0; e < N_elements; e++) {
      if (element_symbol == MPv10_tracer_list[e][0]) {
        // MPv10 species tracer list (2D vector).
        MPv10_tracer_list[e].push_back(tracer_list_replica[s]);
        // MPv10 ion list
        MPv10_ion_list.push_back(tracer_list_replica[s]);
        // Record which element the ion belong to.
        ions_tracer_elem.push_back(e);
        // Map ion string in MPv10_ion_list to their primitive index.
        ions_primindex_map[tracer_list_replica[s]] = first_tracer_index + s;
        // Number of electron contributed by each ion.
        ions_electron_num.push_back(get_electron_num(tracer_list_replica[s]));
        // Position of ions in the primitive vector.
        ions_prim_index.push_back(first_tracer_index + s);
        // Position of ions in the local vector.
        ions_local_index.push_back(ion_count);
        // Increment the counter when ion of the same element is found.
        ion_count++;
        // No of ions each element
        N_sp_by_elem[e]++;
        // No of species in each element
        N_ions_by_elem[e]++;
      }
    }
  }
  // End of tracer loop for identifying ion species ************************

  // Number of species in local (y) vector
  N_sp = MPv10_ion_list.size();  // rename ****

  // record the lesser ions local index for MPv10_ion_list ********************
  for (int i = 0; i < N_sp; i++) {
    if (i != 0
        && get_element_name(MPv10_ion_list[i])
               == get_element_name(MPv10_ion_list[i - 1])) {
      if (get_charge(MPv10_ion_list[i]) > get_charge(MPv10_ion_list[i - 1])) {
        // spdlog::debug("{}: lesser ion exist {}", MPv10_ion_list[i],
        // MPv10_ion_list[i-1]);
        minus_ions_local_index.push_back(i - 1);
      }
    }
    else {
      // spdlog::debug("{} lesser ion do not exist ...", MPv10_ion_list[i]);
      minus_ions_local_index.push_back(-1);
    }
  }
  // End of recording lesser ions local index ********************************

  spdlog::debug("MPv10: Read {} elements and {} Species", N_elements, N_sp);

  // Make species list to generate LUTs *************************************
  for (int e = 0; e < MPv10_tracer_list.size(); e++) {
    for (int i = 0; i < MPv10_tracer_list[e].size(); i++) {

      // Make collisional ionisation tracer list
      if (i < MPv10_tracer_list[e].size() - 1) {
        ci_tracer_list.push_back(MPv10_tracer_list[e][i]);
        ci_tracer_elem.push_back(e);
      }

      // Make recombination tracer list
      if (i != 0) {
        recomb_tracer_list.push_back(MPv10_tracer_list[e][i]);
        recomb_tracer_elem.push_back(e);
      }
    }
  }
  // End of making species list to generate LUTs ****************************
}
// END OF TRACER SPECIES CATALOGER ###########################################



// GET ELECTRON NUMBER DENSITY ################################################
// Calculates electron number density from primitive vector
double MPv10::get_n_elec(const pion_flt *P  ///< primitive state vector array.
)
{
  // get electron number density from P
  int species_counter = 0;
  double ne;
  for (int elem = 0; elem < N_elements; elem++) {  // loop over every element
    int N_elem_species    = N_sp_by_elem[elem];
    double X_elem_density = P[first_tracer_index + elem];

    // loop over every species in THIS element
    for (int s = 0; s < N_elem_species; s++) {
      double y_ion_frac = P[first_tracer_index + N_elements + species_counter];

      // add to ne based on the number of electrons liberated to obtain this ion
      pion_flt number_density = X_elem_density * y_ion_frac;

      int num_elec = ions_electron_num[species_counter];

      ne += num_elec * number_density;
      species_counter++;
    }
  }
  return 10.0;
}
// END OF GET ELECTRON NUMBER DENSITY #########################################



// GET ION NUMBER DENSITY ##############################################
// Calculate given ion's number density from the primitive vector.
double MPv10::get_n_ion(
    string s,                 ///< ion name
    const pion_flt *prim_vec  ///< primitive state vector.
)
{
  if (Tr(s) == -1)
    return -1;

  else {
    // index of the ion in primitive vector
    int ion_primitive_index = Tr(s);
    // index of ion's element
    int element_index =
        ions_tracer_elem[ion_primitive_index - N_elements - first_tracer_index];
    // element number density
    double elem_number_density = prim_vec[RO]
                                 * prim_vec[first_tracer_index + element_index]
                                 / elem_atomic_mass[element_index];
    // ion fraction y_ion = X_ion/X_element;
    double ion_fraction = prim_vec[ion_primitive_index];
    // number denisty of the ion sepcies s.
    return elem_number_density * ion_fraction;
  }
}
// END OF GET ION NUMBER DENSITY ##############################################


// ##################################################################
double MPv10::get_X_H()
{
  spdlog::error("don't call get_X_H() in MPv10");
  exit(1);
  return 0.0;
}
// ##################################################################



// ##################################################################
// Returns the index of the ion in the primitive vector
int MPv10::Tr(const string s)
{
  if (ions_primindex_map.find(s) == ions_primindex_map.end())
    return -1;
  else
    return ions_primindex_map[s];
}
// ##################################################################



// GET TEMPERATURE ##################################################
// Calculate temperature from ion fraction, element number density
// and internal energy using T = E*(gamma-1)/(k_B*n_total)
double MPv10::get_temperature(
    const std::vector<double> &y_ion_frac,
    const std::vector<pion_flt> &elem_number_density,
    const double E)
{
  return gamma_minus_one * E
         / (pconst.kB() * get_ntot(y_ion_frac, elem_number_density));
}
// END OF GET TEMPERATURE ####################################################



// GET TOTAL NUMBER DENSITY ##################################################
// Calculate total number density from y_ion fraction and element
// number density.
double MPv10::get_ntot(
    const std::vector<double> &y_ion_frac,  ///< y_ion_fraction
    const std::vector<pion_flt>
        &element_number_density  ///< element number density
)
{
  int species_counter = 0;
  pion_flt n_tot      = 0;

  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species = N_sp_by_elem[e];
    // neutral frac, got by subtracting off the ion fractions in the next loop.
    pion_flt neutral_frac = 1;

    // loop over every species in THIS element
    for (int s = 0; s < N_elem_species; s++) {
      pion_flt number_density =
          element_number_density[e] * y_ion_frac[species_counter];
      int N_electrons =
          ions_electron_num[species_counter];  // corr: y_ion_num_elem ->
                                               // ions_electron_num

      // add on number of particles got by electrons + ionised atom
      n_tot += (1 + N_electrons) * number_density;

      neutral_frac -= y_ion_frac[species_counter];
      species_counter++;
    }
    n_tot += max(0.0, neutral_frac) * element_number_density[e];
  }
  return n_tot;
}
// END OF GET TOTAL NUMBER DENSITY
// ##################################################



// #############################################################################
// CONVERT PRIMITIVE VECTOR TO LOCAL VECTOR
int MPv10::convert_prim2local(
    const pion_flt *p_in,  ///< primitive vector from grid cell (length N_prim)
    std::vector<double> &p_local,
    int function_flag  /// < flag to say which function called this function
)
{

  // function flag 1 -> Set_Temp()
  // function flag 2 -> TimeUpdateMP_RTnew()
  // function flag 3 -> timescales_RT()
  // function flag 4 -> Temperature()

  if (p_in[PG] > 1.0) exit(0);

  // make sure that p_in[] is a self-consistent state:
  sCMA(corrector, p_in);


  // ==============================================================
  //  Set elemental number density from the current primitive vector
  // ==============================================================
  for (int i = 0; i < N_elements; i++) {  // loop over every element
    elem_number_density[i] =
        p_in[RO]
        * (p_in[elem_prim_index[i]] * corrector[elem_prim_index[i]]
           / elem_atomic_mass[i]);  // CORR: X_elem_atomic_mass ->
                                    // elem_atomic_mass.
  }
  // **************************************************************


  // ==============================================================
  // Set Ion fractions (local vector variables) of each element in
  // the local vector.
  // ==============================================================
  // Note: Ion fractions are obtained by loop over every species with
  // local_vec[i] = prim_vec[i]/X_e, where i runs over mass fraction
  // of ions in the primitive vectors.
  int species_counter = 0;
  // loop over every element
  for (int e = 0; e < N_elements; e++) {
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    // loop over every species in this element
    for (int s = 0; s < N_elem_species; s++) {

      p_local[ions_local_index[species_counter]] =  // CORR: y_ion_index_local
                                                    // -> ions_local_index
          p_in[ions_prim_index[species_counter]]    // CORR: y_ion_index_prim ->
                                                    // ions_prim_index
          * corrector[ions_prim_index[species_counter]]  // CORR:
                                                         // y_ion_index_prim ->
                                                         // ions_prim_index
          / p_in[elem_prim_index[e]];  // CORR: X_mass_frac_index ->
                                       // elem_prim_index

      p_local[ions_local_index[species_counter]] =  // CORR: y_ion_index_local
                                                    // -> ions_local_index
          min(1.0,
              p_local[ions_local_index[species_counter]]);  // CORR:
                                                            // y_ion_index_local
                                                            // ->
                                                            // ions_local_index

      p_local[ions_local_index[species_counter]] =  // CORR: y_ion_index_local
                                                    // -> ions_local_index
          max(0.0,
              p_local[ions_local_index[species_counter]]);  // CORR:
                                                            // y_ion_index_local
                                                            // ->
                                                            // ions_local_index

      species_counter++;
    }
  }
  // **************************************************************


  // ==============================================================
  // Set internal energy density in local vector.
  // ==============================================================
  p_local[E_index] = p_in[PG] / (gamma_minus_one);
  // **************************************************************


#ifdef MPv10_DEBUG
  /*
  //
  // Check for NAN/INF
  //
  for (int v = 0; v < 2; v++) {
    if (!isfinite(p_local[v])) {
      spdlog::error("{}: {}", "INF/NAN input to microphysics", p_local[v]);
      exit(1);
    }
    if (mpv_nH < 0.0 || !isfinite(mpv_nH)) {
      spdlog::error(
          "{}: {}", "Bad density input to MPv10::convert_prim2local", mpv_nH);
      exit(1);
    }
  }
   */
#endif  // MPv10_DEBUG



  return 0;
}
// END OF CONVERT PRIMITIVE VECTOR TO LOCAL VECTOR #############################



// #############################################################################
// CONVERT LOCAL VECTOR TO PRIMITIVE VECTOR
int MPv10::convert_local2prim(
    const std::vector<double> &p_local,
    const pion_flt *p_in,  ///< input primitive vector
    pion_flt *p_out,       ///< updated primitive vector
    int function_flag      /// < flag to say which function called me
)
{
  for (int v = 0; v < N_prim; v++)
    p_out[v] = p_in[v];
  //
  //=================================================================
  //     Update output primitive vector according to local vector
  //     Also, define y_ion_frac while you're at it.
  //=================================================================
  //
  p_out[PG] = p_local[E_index] * (gamma_minus_one);
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp);  // CORR N_species -> N_sp


  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    for (int s = 0; s < N_elem_species;
         s++) {  // loop over every species in THIS element
      p_out[ions_prim_index[species_counter]] =  // CORR: y_ion_index_prim ->
                                                 // ions_prim_index
          p_local[ions_local_index[species_counter]]  // CORR: y_ion_index_local
                                                      // -> ions_local_index
          * p_in[elem_prim_index[e]];  // CORR: X_mass_frac_index ->
                                       // elem_prim_index
      y_ion_frac[species_counter] =
          p_local[ions_local_index[species_counter]];  // CORR:
                                                       // y_ion_index_local ->
                                                       // ions_local_index
      species_counter++;
    }
  }

  // Set mass fraction tracers to be within the required range (not too close
  // to zero or 1)

  species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    p_out[elem_prim_index[e]] =
        max(  // CORR: X_mass_frac_index -> elem_prim_index
            Min_NeutralFrac,
            min(Max_NeutralFrac,
                static_cast<double>(
                    p_out[elem_prim_index[e]])));  // CORR: X_mass_frac_index ->
                                                   // elem_prim_index

    for (int s = 0; s < N_elem_species;
         s++) {  // loop over every species in THIS element

#ifdef MPv10_DEBUG
      // Introducing sense checks -- make sure value is positive, less than 1
      if (static_cast<double>(
              p_out[ions_prim_index[species_counter]])  // CORR:
                                                        // y_ion_index_prim ->
                                                        // ions_prim_index
          < (-2 * MPv10_ABSTOL)) {
        spdlog::debug(
            "convert_local2prim: {} mass fraction goes negative here. \n [",
            function_flag);
        for (int v = 0; v < N_prim; v++) {
          spdlog::debug("{}, ", p_out[v]);
        }
        spdlog::debug("]");
        print_flag = 1;
      }

      else if (
          static_cast<double>(
              p_out[ions_prim_index[species_counter]])  // CORR:
                                                        // y_ion_index_prim ->
                                                        // ions_prim_index
          > (1 + MPv10_ABSTOL)
                * static_cast<double>(
                    p_out[elem_prim_index[e]])) {  // CORR: X_mass_frac_index ->
                                                   // elem_prim_index
        spdlog::debug(
            "convert_local2prim: {} mass frac too large for species {}: X = {}\nPrim vector: \n [",
            function_flag, s,
            p_out[ions_prim_index[species_counter]]);  // CORR: y_ion_index_prim
                                                       // -> ions_prim_index
        for (int v = 0; v < N_prim; v++) {
          spdlog::debug("{}, ", p_out[v]);
        }
        spdlog::debug("]");
        print_flag = 1;
      }
#endif

      p_out[ions_prim_index[species_counter]] = max(  // CORR: y_ion_index_prim
                                                      // -> ions_prim_index
          Min_NeutralFrac,
          min(static_cast<double>(
                  p_out[elem_prim_index[e]])  // CORR: X_mass_frac_index ->
                                              // elem_prim_index
                  * Max_NeutralFrac,
              static_cast<double>(
                  p_out[ions_prim_index
                            [species_counter]])));  // CORR:
                                                    // y_ion_index_prim
                                                    // ->
                                                    // ions_prim_index
      species_counter++;
    }
  }

  // Set output pressure to be within required temperature range (use the
  // possibly corrected x(H+) from p_out[]).
  //
  double T = get_temperature(y_ion_frac, elem_number_density, p_local[E_index]);
  // cout <<"current temp = " << T<<"\n";
  if (T > 1.01 * EP->MaxTemperature) {
    Set_Temp(p_out, EP->MaxTemperature, 0);
#ifdef MPv10_DEBUG
    spdlog::debug(
        "MPv10::convert_local2prim() HIGH temperature encountered. T={}, obtained from nH={}, eint={} x={}...  limiting to T={}",
        T, mpv_nH, p_local[E_index], 10  // p_out[pv_H1p]
        ,
        EP->MaxTemperature);
#endif  // MPv10_DEBUG
  }
  if (T < 0.99 * EP->MinTemperature) {
    spdlog::warn("Low temperature!!!");
    spdlog::warn("T={}", T);
    Set_Temp(p_out, EP->MinTemperature, 0);
    spdlog::warn("MPv10::convert_local2prim() LOW temperature encountered.");
    spdlog::debug(
        "T={}, obtained from nH={}, eint={}, ...  limiting to T={}", T, mpv_nH,
        p_local[E_index], EP->MinTemperature);
  }
  p_out[PG] = p_local[E_index] * (gamma_minus_one);



  return 0;
}
// END OF CONVERT LOCAL VECTOR TO PRIMITIVE VECTOR #############################



// #############################################################################
// SET LOCAL ION FRACTIONS
void MPv10::set_y_ion_frac(
    const std::vector<double> &p_local, std::vector<double> &y_ion_frac)
{
  for (int s = 0; s < N_sp; s++) {  // CORR N_species -> N_sp
    y_ion_frac[s] = p_local[s];
    // Set y_ion_frac within physically acceptable range.
    y_ion_frac[s] = max(1e-20, y_ion_frac[s]);
    y_ion_frac[s] = min(1.0, y_ion_frac[s]);
  }
  return;
}
// END OF SET LOCAL ION FRACTIONS ##############################################



// #############################################################################
// SET LOCAL NEUTRAL FRACTIONS
void MPv10::set_y_neutral_frac(
    const std::vector<double> &y_ion_frac, std::vector<double> &y_neutral_frac)
{
  // Calculate neutral fraction from the current local vector.
  int sct = 0;
  // loop over every element
  for (int elem = 0; elem < N_elements; elem++) {
    y_neutral_frac[elem] = 1.0;
    // loop over every species in this element
    for (int s = 0; s < N_sp_by_elem[elem];
         s++) {  // CORR N_species_by_elem -> N_sp_by_elem
      y_neutral_frac[elem] -= y_ion_frac[sct];
      sct++;
    }
    // Ensure the neutral frac is within the acceptable range
    y_neutral_frac[elem] = max(1.0e-20, y_neutral_frac[elem]);
    y_neutral_frac[elem] = min(1.0, y_neutral_frac[elem]);
  }
  return;
}
// END OF LOCAL NEUTRAL FRACTIONS ##############################################



// #############################################################################
// SET ELECTRON NUMBER DENSITY
void MPv10::set_ne(const std::vector<double> &y_ion_frac, double &ne)
{


  // Calculate electron number density
  int sct = 0;
  // loop over every element
  for (int elem = 0; elem < N_elements; elem++) {
    // loop over every species in THIS element
    for (int s = 0; s < N_sp_by_elem[elem];
         s++) {  // CORR N_species_by_elem -> N_sp_by_elem
      ne += ions_electron_num[sct] * elem_number_density[elem]
            * y_ion_frac[sct];  // CORR y_ion_num_elec -> ions_electron_num
      sct++;
    }
  }
  // Ensure electron number density is not negative.
  if (ne < 0.0) {
    spdlog::info("Warning: Negative ne = {}, resetting to {}", ne, 0.0);
    ne = max(0.0, ne);
  }

  return;
}
// END OF SET ELECTRON NUMBER DENSITY ##########################################



// #############################################################################
// SET LOCAL VARIABLES

// Set local variables local ion fraction (y_ion_frac), neutral fraction
// (y_neutral_frac), electron number density.
void MPv10::set_localvariables(
    const std::vector<double> &p_local,
    std::vector<double> &y_ion_frac,
    std::vector<double> &y_neutral_frac,
    double &ne)
{
  // Set y ion fraction from the local vector.
  set_y_ion_frac(p_local, y_ion_frac);

  // Set local neutral fraction from local ion fraction.
  set_y_neutral_frac(y_ion_frac, y_neutral_frac);

  // Calculate electron number density from local ion fraction.
  set_ne(y_ion_frac, ne);

  return;
}
// END OF CONVERT LOCAL VECTOR TO PRIMITIVE VECTOR #############################



// ##################################################################
// ##################################################################


double MPv10::Temperature(
    const pion_flt *pv,  ///< primitive vector
    const double         ///< eos gamma
)
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  if (pv[RO] <= 0.0 || pv[PG] <= 0.0) {
    // cout <<"MPv10::Temperature() negative rho="<<pv[RO]<<" or
    // p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (y,nE,Eint), and get Temperature from it.
  //
  std::vector<double> p_local;
  p_local.resize(N_local);
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp);  // CORR N_species -> N_sp

  convert_prim2local(pv, p_local, 4);

  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    // loop over every species in THIS element
    for (int s = 0; s < N_elem_species; s++) {
      // CORR: y_ion_index_prim -> ions_prim_index
      // CORR: X_mass_frac_index -> elem_prim_index
      y_ion_frac[species_counter] =
          pv[ions_prim_index[species_counter]] / pv[elem_prim_index[e]];
      species_counter++;
    }
  }
  return (get_temperature(y_ion_frac, elem_number_density, p_local[E_index]));
}



// ##################################################################
// ##################################################################
int MPv10::Set_Temp(
    pion_flt *p_pv,  ///< primitive vector.
    const double T,  ///< temperature
    const double     ///< eos gamma.
)
{

  //
  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.
  //
  // cout << "set temperature\n";
  if (p_pv[PG] <= 0.0) {
    spdlog::debug("MP_Hydrogen::Set_Temp() correcting negative pressure");
    p_pv[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }

  std::vector<double> p_local;
  p_local.resize(N_local);
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp);  // CORR N_species -> N_sp

  int err = convert_prim2local(p_pv, p_local, 1);

  // Determine y_ion_frac from the primitive vector
  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {   // loop over every element
    int N_elem_species = N_sp_by_elem[e];  // CORR N_species -> N_sp
    // loop over every species in THIS element
    // CORR: y_ion_index_prim -> ions_prim_index
    // CORR: X_mass_frac_index -> elem_prim_index
    for (int s = 0; s < N_elem_species; s++) {
      y_ion_frac[species_counter] =
          p_pv[ions_prim_index[species_counter]] / p_pv[elem_prim_index[e]];
      species_counter++;
    }
  }


  // Determine internal energy using get_ntot
  p_local[E_index] = get_ntot(y_ion_frac, elem_number_density) * pconst.kB() * T
                     / (gamma_minus_one);


  // Call convert_local2prim with the new local vector; this will
  // generate a new temperature value;
  err += convert_local2prim(p_local, p_pv, p_pv, 1);
  return err;
}
// ##################################################################



// ##################################################################
// TIME UPDATE MP
int MPv10::TimeUpdateMP(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    pion_flt *p_out,       ///< Destination Vector for updated values.
    const double dt,       ///< Time Step to advance by.
    const double,          ///< EOS gamma.
    const int,             ///< Switch for what type of integration to use.
    double *random_stuff   ///< Vector of extra data (column densities, etc.).
)
{
  //
  // Call the new update function, but with zero radiation sources.
  //
  std::vector<struct rt_source_data> temp;
  int err =
      TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}
// END OF TIME UPDATE MP ############################################



// ##################################################################
// TIME UPDATE MP RT NEW
int MPv10::TimeUpdateMP_RTnew(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    const int,             ///< unused.
    const std::vector<struct rt_source_data> &,  ///< unused.
    const int,                                   ///< unused.
    std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    pion_flt *p_out,      ///< Destination Vector for updated values
                          ///< (can be same as first Vector.
    const double dt,      ///< Time Step to advance by.
    const double,         ///< EOS gamma.
    const int,            ///< Switch (unused)
    double *random_stuff  ///< final temperature (debugging).
)
{

  int err = 0;
  std::vector<double> P;
  P.resize(N_local);

  err = convert_prim2local(p_in, P, 2);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to MPv10::TimeUpdateMP()", err);
    exit(1);
  }

  correct_localvector(P, 1);

  setup_radiation_source_parameters(p_in, P, ion_src);

  // Populates CVODE vector with initial conditions (input)
  for (int v = 0; v < N_local; v++)
    NV_Ith_S(y_in, v) = P[v];


  // Calculate y-dot[] to see if anything is changing significantly over dt
  // Here y_out is y_dot calculated from the ydot function.
  double maxdelta = 0.0;
  err             = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in MPv10::TimeUpdateMP_RTnew()",
        err);
    exit(1);
  }

  // Here y_out is y_dot calculating maxdelta.
  for (int v = 0; v < N_local; v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out, v) * dt / NV_Ith_S(y_in, v)));
  }

  // Now if nothing is changing much, just to a forward Euler integration.
  // Here y_out on the RHS is y_dot, which is over-written by y_out
  // calculated using euler integration.
  if (maxdelta < EULER_CUTOFF) {
    // cout <<"maxdelta="<<maxdelta<<", Doing euler integration...\n";
    for (int v = 0; v < N_local; v++) {
      NV_Ith_S(y_out, v) = NV_Ith_S(y_in, v) + dt * NV_Ith_S(y_out, v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration.
  // Here y_out (which is y_dot) is input, and the output is over-written onto
  // y_out by integrate_cvode_step.
  else {
    // cout <<"maxdelta="<<maxdelta<<", Doing CVODE integration...\n";
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      spdlog::error(
          "{}: {}", "integration failed: MPv10::TimeUpdateMP_RTnew()", err);
      exit(1);
    }
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v = 0; v < N_local; v++)
    P[v] = NV_Ith_S(y_out, v);

  correct_localvector(P, 2);

  err = convert_local2prim(P, p_in, p_out, 2);

#ifdef TEST_INF
  for (int v = 0; v < N_prim; v++) {
    if (!isfinite(p_in[v]) || !isfinite(p_out[v])) {
      spdlog::debug("NAN in MPv3 update: {}", v);
      spdlog::debug("Pin  : {}", p_in);
      spdlog::debug("Pout : {}", p_out);
      spdlog::debug("Ploc  : {}", P);
      // spdlog::error("{}: {}", "NAN in MPv3",P[2]);
      return 1;
    }
  }
#endif

  return err;
}
// END OF TIME UPDATE MP RT NEW #####################################



// CORRECT LOCAL VECTOR (CLV) #####################################
int MPv10::correct_localvector(
    std::vector<double> &p_local,  ///< check the local vector
    int function_flag              /// flags the function
)
{
  // function_flag 1 -> towards begining of TimeUpdateMP_RTnew
  // function_flag 2 -> towards the end of TimeUpdateMP_RTnew
  // function_flag 3 -> timescales_RT


  // ==============================================================
  // Check if y_ion_frac is within the acceptable range
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp);

  for (int i = 0; i < N_sp; i++) {  // loop over ions in p_local
    y_ion_frac[i] = p_local[i];
    // Set y_ion_frac within physically acceptable range.
    if (y_ion_frac[i] < 0.0) {
      if (y_ion_frac[i] < -MPv10_ABSTOL) {
        spdlog::warn(
            "MPv10::CLV flag:{} - Negative ion fraction for vector index ={}, {}"
            ", resetting to {}",
            function_flag, i, y_ion_frac[i], MPv10_ABSTOL);
      }
      // resetting
      y_ion_frac[i] = MPv10_ABSTOL;
    }

    if (y_ion_frac[i] > 1.0) {
      if (y_ion_frac[i] > 1.0 + MPv10_RELTOL) {
        spdlog::warn(
            "MPv10::CLV flag: {} - Ion fraction >1 for vector index = {}, {}"
            ", resetting to {}",
            function_flag, i, y_ion_frac[i], 1.0);
      }
      // resetting
      y_ion_frac[i] = 1.0;
      // exit(0);
    }
  }
  // **************************************************************

  // Feeding the corrected y_ion_frac back in to p_local
  for (int i = 0; i < N_sp; i++)
    p_local[i] = y_ion_frac[i];

  // ==============================================================
  // Check for negative pressure
  // Note: This shouldn't happen, so we output a warning) and set to
  // 10K if we find it.
  if (p_local[E_index] <= 0.0) {
    spdlog::warn(
        "MPv10::CLV flag: {} - Negative pressure input: e = {}, "
        "setting to {} K",
        function_flag, p_local[E_index], EP->MinTemperature);

    // reset the internal energy (requires using y_ion_frac in get_ntot)
    p_local[E_index] = get_ntot(y_ion_frac, elem_number_density) * pconst.kB()
                       * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************


  // ==============================================================
  // if the temperature is below 10 K, set it back to 10 K by set the internal
  // energy accordingly
  double T = get_temperature(y_ion_frac, elem_number_density, p_local[E_index]);

  if (T < EP->MinTemperature) {
    spdlog::warn(
        "MPv10::CLV flag: {} - Temperature = {} K, below T_MIN = {} K",
        function_flag, T, EP->MinTemperature);

    // Resetting the internal energy with Minimum Temperature
    p_local[E_index] = get_ntot(y_ion_frac, elem_number_density) * pconst.kB()
                       * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************


  return 0;
}
// END OF CORRECT LOCAL VECTOR ######################################


// ##################################################################
double MPv10::timescales(
    const pion_flt *p_in,  ///< Current cell state vector.
    const double,          ///< EOS gamma.
    const bool,            ///< set to 'true' if including cooling time.
    const bool,            ///< set to 'true' if including recombination time.
    const bool             ///< set to 'true' if including photo-ionsation time.
)
{
#ifdef MPv10_DEBUG
  if (RS->Nsources != 0) {
    spdlog::debug("WARNING: MPv10::timescales() using non-RT version!");
  }
#endif  // MPv10_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin = timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}
// ##################################################################


// TIME SCALE RT ###########################################################
/// This returns the minimum timescale of all microphysical processes,
/// including reaction times for each species and the total heating/cooling
/// time for the gas. It requires the radiation field as an input, so it has
/// substantially greater capability than the other timescales function.
/// Default setting is DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double MPv10::timescales_RT(
    const pion_flt *p_in,  ///< Current cell state vector.
    const int N_heat,      ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    const double  ///< EOS gamma.
)
{

  // Note: Minimum value of micro-physics trace variables in the primitive
  // vector is set to 1e-12.
  int err = 0;
  std::vector<double> p_local;
  p_local.resize(N_local);

  err = convert_prim2local(p_in, p_local, 3);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to MPv10::timescales_RT()", err);
    exit(1);
  }

  // correct local vector
  correct_localvector(p_local, 3);

  // v runs over all ions and internal energy
  for (int v = 0; v < N_local; v++)
    NV_Ith_S(y_in, v) = p_local[v];


  // Now calculate y-dot[]...
  // ydot function returns RHS of ydot equation. Here y_out is RHS of ydot.
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in MPv10::timescales_RT()", err);
    exit(1);
  }

  double tt = Temperature(p_in, 0);
  // spdlog::info("{:12.6e}  {:12.6e}", tt, NV_Ith_S(y_out,E_index));

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t = HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the
  // energy change (heating timescale) is not required for accurately tracking
  // ionisation fronts (although it may be needed for cooling!). For testing
  // purposes there are ifdefs to allow the code to use a relative change in
  // neutral fraction and/or the relative change in energy as the timestep
  // criterion, rather than the default of absolute change in neutral
  // fraction.
  //

  std::vector<double> Y_dot;
  Y_dot.resize(N_local);
  // Since y_out is ydot.
  for (int v = 0; v < N_local; v++)
    Y_dot[v] = NV_Ith_S(y_out, v);

  for (int v = 0; v < N_equations; v++) {
    // t = min(t, DTFRAC / (fabs(NV_Ith_S(y_out, v)) + TINYVALUE));
    t = min(t, DTFRAC / (fabs(Y_dot[v]) + TINYVALUE));
    if (t < 1e5) {
      cout << "limit by dx: dt=" << t << "\n";
      for (int v = 0; v < N_local; v++)
        cout << "p_local[" << v << "]= " << p_local[v] << ",  ";
      cout << "\n";
      for (int v = 0; v < N_local; v++)
        cout << "y_out[" << v << "]= " << Y_dot[v] << ",  ";
      cout << "\n";
      //   spdlog::debug("MPv3:: EP and RS: {}\t{}", fmt::ptr(EP),
      //   fmt::ptr(RS));
    }
  }

  return t;
}



// ##################################################################
// ##################################################################
void MPv10::sCMA(
    std::vector<double> &corrector,  ///< input corrector vector
    const pion_flt *p_in             ///< input primitive vector from grid cell
)
{
  // "Consistent multi-species advection" correction step from
  // Plewa and Muller (1999).

  //  Re-initialise corrector every step
  for (int i = 0; i < N_prim; i++)
    corrector[i] = 1.0;
  int print_flagg        = 0;
  double total_mass_frac = 0;

  // loop over every species and get the sum
  int species_counter = 0;

  // Calculate all-element correction
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    total_mass_frac +=
        p_in[elem_prim_index[e]];  // CORR: X_mass_frac_index -> elem_prim_index
  }
  double e_correction = 1.0 / total_mass_frac;
  species_counter     = 0;
  // apply all-element correction, calculate species correction, apply species
  // correction
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species =
        N_sp_by_elem[e];  // CORR N_species_by_elem -> N_sp_by_elem
    corrector[elem_prim_index[e]] =
        e_correction;  // correct THIS element //CORR: X_mass_frac_index ->
                       // elem_prim_index
    // Calculate all-species-pr-element correction, if needed, i.e.
    double s_frac = 0;

    for (int s = 0; s < N_elem_species; s++) {
      s_frac +=
          p_in[ions_prim_index[species_counter]];  // CORR: y_ion_index_prim ->
                                                   // ions_prim_index
      species_counter++;
    }

    if (s_frac
        > ((p_in[elem_prim_index[e]] * e_correction)
           - Min_NeutralFrac)) {  // CORR: X_mass_frac_index -> elem_prim_index
      print_flagg = 1;
      double s_correction =
          ((p_in[elem_prim_index[e]] * e_correction)
           - Min_NeutralFrac)  // CORR: X_mass_frac_index -> elem_prim_index
          / s_frac;
      int inner_species_counter = (species_counter - N_elem_species);
      for (int s = 0; s < N_elem_species; s++) {
        corrector[ions_prim_index[inner_species_counter]] =
            s_correction;  // CORR: y_ion_index_prim -> ions_prim_index
      }
    }
  }
  return;
}



// ##################################################################
// YDOT FUNCTION
// This function calculate the RHS of ydot equations.
int MPv10::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double *         ///< extra user-data vector (UNUSED)
)
{

  //========================================================
  // Determine ne, y_ion_frac and neutral fraction from the current local vector
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp);  // CORR N_species -> N_sp
  std::vector<double> y_neutral_frac;
  y_neutral_frac.resize(N_elements);
  double ne = 0.0;

  std::vector<double> p_local;
  p_local.resize(N_local);
  // make a copy of y_now into p_local
  for (int v = 0; v < N_local; v++)
    p_local[v] = NV_Ith_S(y_now, v);

  // Set local quantities: y_ion_frac, y_neutral_frac and ne
  set_localvariables(p_local, y_ion_frac, y_neutral_frac, ne);
  // TODO: set species number density
  // set_species_number_density(y_ion_frac, y_neutral_frac,
  // species_number_density);

#ifdef MPv10_CIE_TEST
  // Restricting ne value not less than 1e-4 to perform CIE test.
  // CORR: Remove this line, this line is unphysical.
  ne = max(1e-4, ne);
#endif

  // Get internal energy from the current local vector (y_now)
  double E_in = NV_Ith_S(y_now, E_index);

  // Initialise all elements of ydot vector to zero.
  for (int v = 0; v < N_equations; v++)
    NV_Ith_S(y_dot, v) = 0.0;

  // Calculate temperature from E_in, y_ion_frac and number density.
  double T = get_temperature(y_ion_frac, elem_number_density, E_in);
  // cout<< "ydot: " << "ne = " << ne << " | "<< "T = " << T << " | " <<endl;

  // Initializing Edot to zero
  double Edot = 0.0;

  // Restricting temperature within T_min and T_max
  if (T > T_max)
    T = T_max;
  else if (T < T_min)
    T = T_min;


  // SECTION: Collisional Ionisation *******************************************

  // Calculating RHS terms corresponding to collisional ionisation for
  // species in ci_tracer_list. These terms are calculated only once and
  // attached at relevant ydot equations.
  //
  // For example, for species s, if ion, the term ci_rate(s)y(s)ne is subtracted
  // for ydot(s) while the same terms is added for the ydot(s+1) equation,
  // provided the specie s+1 exist. Whereas if the species s is neutral, then
  // only addition is performed.

  double ci_rate    = 0.0;
  double this_y_dot = 0.0;

  for (int s = 0; s < ci_tracer_list.size(); s++)  // loop over ci species
  {
    // corr y_im1_index_tables with ci tracer list index
    ci_rate = collisional_ionisation_rate(s, T);

    if (minus_ions_local_index[s] != -1) {
      // if the less ionised species is not neutral
      this_y_dot = ci_rate * NV_Ith_S(y_now, minus_ions_local_index[s])
                   * ne;  // corr: y_im1_index_local -> minus_ions_local_index
      NV_Ith_S(y_dot, minus_ions_local_index[s]) -= this_y_dot;
    }
    else {
      // if the less ionised species IS neutral
      this_y_dot = ci_rate * y_neutral_frac[ci_tracer_elem[s]] * ne;
    }

    NV_Ith_S(y_dot, ions_local_index[s]) +=
        this_y_dot;  // corr: y_ion_index_local -> ions_local_index

    // Cooling due to collisional ionisation of this species.
    Edot -= ci_ion_pot[s] * this_y_dot * elem_number_density[ci_tracer_elem[s]];
  }
  // END OF SECTION: Collisional Ionisation ************************************


  // SECTION: Recombination (Radiative + Dielectronic) *************************

  // RHS terms corresponding to Recombination (R) are calculated for
  // all species. In this case, for specie s, the term rec_rate(s)y(s)ne is
  // substrated for ydot(s-1) and the same terms is added for the ydot(s),
  // provided species s is a ion of same element. Whereas
  // if the species s is neutral, then only addition performed.

  double recomb_rate = 0.0;
  this_y_dot         = 0.0;

  for (int s = 0; s < recomb_tracer_list.size();
       s++)  // loop over recomb species
  {
    recomb_rate = recombination_rate(s, T);

    this_y_dot = recomb_rate * NV_Ith_S(y_now, ions_local_index[s])
                 * ne;  // CORR: y_ion_index_local -> ions_local_index

    // Subtract this term to the current species equation
    NV_Ith_S(y_dot, ions_local_index[s]) -=
        this_y_dot;  // CORR: y_ion_index_local -> ions_local_index
    // add this term to less ionised species equation provided it
    // is not neutral species
    if (minus_ions_local_index[s] != -1)
      NV_Ith_S(y_dot, minus_ions_local_index[s]) +=
          this_y_dot;  // corr: y_im1_index_local -> minus_ions_local_index

    // Heating due to recombination of this species.
    Edot -= (3. / 2.) * T * pconst.kB() * this_y_dot
            * elem_number_density[recomb_tracer_elem[s]];
  }
  // END OF SECTION: Recombination *********************************************


  // SECTION: Mellema cooling function *****************************************

  // Edot source term associated with mellema Cooling is calculated in this
  // Section.
  double Lambda         = 0.0;
  double L              = 0.0;
  double number_density = 0.0;

#ifdef MPv10_CIE_COOLING
  mc_outfile.open("cooling_function.txt", ios::app);
  mc_outfile << T << "  ";
#endif

  int sct           = 0;
  int species_index = 0;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      // Calculating the current number density of neutral atoms.
      if (j == 0) {
        number_density = y_neutral_frac[i] * elem_number_density[i];
      }
      // Calculating the current number density of ions.
      else {
        number_density = y_ion_frac[sct] * elem_number_density[i];
        sct++;
      }

      // mellema cooling rate
      // L = mellema_cooling_rate(species_index, T, ne);
      // chianti cooling rate
      L = chianti_cooling_rate(species_index, T, ne);
      Lambda += ne * number_density * L;
      species_index++;

#ifdef MPv10_CIE_COOLING
      // cooling function of individual species ne*number_density*L
      mc_outfile << ne * number_density * L << "  ";
#endif
    }
  }
#ifdef MPv10_CIE_COOLING
  mc_outfile << "\n";
  mc_outfile.close();
#endif

  Edot -= Lambda;
  // END OF SECTION: Mellema Cooling ******************************************


  //
  // We want to limit cooling as we approach the minimum temperature, so we
  // scale the rate to linearly approach zero as we reach T_min.
  //
  if (Edot < 0.0 && T < 2.0 * EP->MinTemperature) {
    Edot = min(0.0, (Edot) * (T - EP->MinTemperature) / EP->MinTemperature);
  }


  if (!EP->update_erg) Edot = 0.0;
  NV_Ith_S(y_dot, E_index) = Edot;

  return 0;
}
// END OF YDOT FUNCTION ########################################################



// GET COLLISIONAL IONISATION RATE FROM LUT ################################
double MPv10::collisional_ionisation_rate(
    const int species_index,  // species index (first index)
    const double T            // temperature
)
{

  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];

  double rate = collisional_rate_table[species_index][T_index]
                + dT * collisional_slope_table[species_index][T_index];

  return rate;
}
// END OF GET COLLISIONAL IONISATION RATE FROM LUT #############################



// GET RECOMBINATION RATE ####################################################
double MPv10::recombination_rate(
    const int species_index,  // species identifier
    const double T            // temperature
)
{

  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  double dT   = T - T_table[T_index];

  double rate = recombination_rate_table[species_index][T_index]
                + dT * recombination_slope_table[species_index][T_index];

  return rate;
}
// END OF GET RECOMBINATION RATE ############################################



// GET MELLEMA COOLING RATE #################################################
double MPv10::mellema_cooling_rate(
    const int location,  // species identifier (database location)
    double T,            // Temperature
    double ne            // electron number density
)
{

  // 1. Set up log10(T) and log10(ne) indices to access LUT.
  // A conditional statement is added to access the Mellema Cooling table only
  // if ne and T range lie within the table in the database of
  // Mellema_/// rate.cpp

  if (T > T_max) T = T_max;
  if (T < T_min) T = T_min;
  if (ne < ne_min) ne = ne_min;

  // Note: Since the Mellema cooling rates are given in log scales,
  // i.e, L(Log T, Log ne), we calculate the differential change dLog_T
  // (step) and dLog_ne along with their indices. Delta_Log_T and
  // Delta_Log_ne is a constant step size calculated in the
  // Mellema_cooling_rate.cpp

  // Get temperature vector index
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  // Get the difference in temperature.
  double dT = T - T_table[T_index];

  // Get electron number density vector index.
  int ne_index = int((log10(ne) - log10(ne_min)) / delta_logne);
  // Get the difference in electron number density.
  double dne = ne - ne_table[ne_index];

  // cout<<"T= "<<T<< ", "<<"T_index= "<<T_Index<<", "<<"T= "<<Vec_T[T_Index]
  // << endl; cout<<"ne= "<<ne<< ", "<<"ne_index= "<<ne_Index<<", "<<"ne=
  // "<<Vec_ne[ne_Index] << endl;

  // Calculating Mellema Cooling Rate.
  double L = 0.0;
  // Calculating cooling from the LookUp Table (LUT) with the location of
  // each Ion and atom in the database.
  /*
  L = mellema_table[location].rate[T_index][ne_index]
          + dT * mellema_table[location].T_rateslope[T_index][ne_index]
          + dne * mellema_table[location].ne_rateslope[T_index][ne_index];
  */
  return L;
}
// END OF GET MELLEMA COOLING RATE ############################################

// GET CHIANTI COOLING RATE #################################################
double MPv10::chianti_cooling_rate(
    const int location,  // species identifier (database location)
    double T,            // Temperature
    double ne            // electron number density
)
{

  // 1. Set up log10(T) and log10(ne) indices to access LUT.
  // A conditional statement is added to access the Mellema Cooling table only
  // if ne and T range lie within the table in the database of
  // Mellema_/// rate.cpp

  if (T > T_max) T = T_max;
  if (T < T_min) T = T_min;
  if (ne < ne_min) ne = ne_min;

  // Note: Since the Mellema cooling rates are given in log scales,
  // i.e, L(Log T, Log ne), we calculate the differential change dLog_T
  // (step) and dLog_ne along with their indices. Delta_Log_T and
  // Delta_Log_ne is a constant step size calculated in the
  // Mellema_cooling_rate.cpp

  // Get temperature vector index
  int T_index = int((log10(T) - log10(T_min)) / delta_logT);
  // Get the difference in temperature.
  double dT = T - T_table[T_index];

  // Get electron number density vector index.
  int ne_index = int((log10(ne) - log10(ne_min)) / delta_logne);
  // Get the difference in electron number density.
  double dne = ne - ne_table[ne_index];

  // cout<<"T= "<<T<< ", "<<"T_index= "<<T_Index<<", "<<"T= "<<Vec_T[T_Index]
  // << endl; cout<<"ne= "<<ne<< ", "<<"ne_index= "<<ne_Index<<", "<<"ne=
  // "<<Vec_ne[ne_Index] << endl;

  // Calculating Mellema Cooling Rate.
  double L = 0.0;
  // Calculating cooling from the LookUp Table (LUT) with the location of
  // each Ion and atom in the database.
  L = chianti_table[location].rate[T_index][ne_index]
      + dT * chianti_table[location].T_rateslope[T_index][ne_index]
      + dne * chianti_table[location].ne_rateslope[T_index][ne_index];

  return L;
}
// END OF GET MELLEMA COOLING RATE ############################################

// #############################################################################
// ***************************** END OF MPv10 **********************************
// #############################################################################



// #############################################################################
// UNUSED FUNCTIONS


// ##################################################################
void MPv10::setup_radiation_source_parameters(
    const pion_flt *p_in,  ///< primitive input state vector.
    vector<double> &P,     ///< local input state vector (x_in,E_int)
    std::vector<struct rt_source_data> &ion_src
    ///< list of ionising src column densities and source properties.
)
{
  for (unsigned int v = 0; v < ion_src.size(); v++)
    rt_data[v] = ion_src[v];

  // struct rt_source_data {
  //  double Vshell;   ///< Shell volume for discrete
  //  photo-ionisation/-heating rates. double dS;       ///< Path length
  //  through cell. double strength[MAX_TAU]; ///< Luminosity (or flux if
  //  source at infinity). double Column[MAX_TAU];  ///< integral of
  //  quantities along LOS to near edge of cell. double DelCol[MAX_TAU];  ///<
  //  integral of quantities along LOS through cell. int id;   ///< source id.
  //  int type; ///< diffuse-radiation or a real source. short unsigned int
  //  NTau; ///< Number of LOS quantities traced for the source.
  //};

  // we want to use these data in ydot to loop over each radiation
  // source, and calculate the ionization rate (using the photon-
  // conserving formula) for each species.
  // - Strength = erg/s/bin luminosity of source.
  // - Column = Tau to edge of cell
  // - DelCol = dTau through cell (we'll re-calculate this in ydot)
  // - Vshell and dS are obvious from Mellema et al. paper.

  return;
}



// ##################################################################
void MPv10::get_dtau(
    const pion_flt ds,     ///< ds, thickness of the cell
    const pion_flt *p_in,  ///< input primitive vector
    pion_flt *dtau_vec     ///< output dtau vector
)
{
  // get optical depth through cell for each photon energy range.
  return;
  /*  for (int bin = 0; bin < get_nbins(); bin++) {
      double dtau         = 0;  // sum dtau across all species within this bin
      int species_counter = 0;

      for (int e = 0; e < N_elements; e++) {  // loop over every element
        int N_elem_species = N_species_by_elem[e];

        for (int s = 0; s < N_elem_species; s++) {  // loop over every species
          double n_s =
              p_in[RO]
              * (p_in[y_ion_index_prim[species_counter]] /
    X_elem_atomic_mass[e]); double xsec =
    y_ion_xsections[species_counter][bin]; dtau += n_s * xsec * ds;
          species_counter++;
        }
      }
      dtau_vec[bin] = dtau;
    }
    return;
    */
}



// ##################################################################
int MPv10::set_multifreq_source_properties(
    // int tempin,
    // double *output
    const struct rad_src_info *rsi,  ///< source data
    double *output  ///< O/P source luminosity per energy bin (erg/s/bin).
)
{
  if (rsi->effect != RT_EFFECT_MFION) {
    spdlog::error("{}: {}", "Wrong source type for id", rsi->id);
    exit(1);
  }
  // Create empty table to store temperature / flux values -- nb There are 65
  // temperatures recorded for logg=4
  double flux_table[65][12];
  double temps[65];

  // Grid of fluxes from Castelli, compiled into one text file of fractional
  // luminosity / bin in a separate python script Original data from:
  // http://www.oact.inaf.it/castelli/castelli/grids/gridp00k2odfnew/fp00k2tab.html
  FILE *pf = fopen("./source/microphysics/castelli.txt", "r");

  /*
  fscanf(pf, "%*[^\n]");  // Read and discard first line, as it's just a header
  for (int i = 0; i < 65; i++) {  // Loop over remaining lines to populate table
    fscanf(pf, "%lf", &temps[i]);  // record first entry, i.e. temperature

    for (int j = 0; j < 12;
         j++) {  // record every other column as an entry in the flux table
      fscanf(pf, "\t%lf", &flux_table[i][j]);
    }
    fscanf(pf, "\n");  // progress to next line...
  }
  */
  // To interpolate, first find which temperatures our star's temperature lies
  // between Start by finding the first T in temps greater than our star's temp,
  // i.e. upper value to interpolate from
  int i = 0;
  while (rsi->Tstar > temps[i] && i < 65)
    i++;

  double Tlower = temps[i - 1];
  double Tupper = temps[i];

  // Then linearly interpolate between flux bins above / below to determine flux
  // bins for the given temperature
  for (int j = 0; j < 12; j++) {
    output[j] =
        flux_table[i - 1][j] * (1 - (rsi->Tstar - Tlower) / (Tupper - Tlower))
        + flux_table[i][j] * (1 - (Tupper - rsi->Tstar) / (Tupper - Tlower));
    // multiply by the total luminosity to have output in erg/s/bin:
    output[j] = output[j] * rsi->strength;
    spdlog::debug("Output[{}] = {}", j, output[j]);
  }

  // TODO: \Maggie{I think this is now done -- need to verify, though!}
  // Now we need to figure out how to get the luminosity of the star
  // in each frequency bin, in erg/s/bin.
  // * rsi->strength gives the luminosity of the star in erg/s
  // * rsi->Tstar    gives the effective temperature of the star.
  // * rsi->Rstar    gives the Radius of the star.
  // If the star were a blackbody, then this would be enough to
  // calculate the luminosity in each bin, if we have the bin ranges
  // set (which we do).  Unfortunately a BB is a bad approximation.
  // Maybe it is the best we can do for now.
  //
  // We want to add the luminosity in each bin to the array "str".
  return 0;
}



// #############################################################################
// ************************** DEBUG  FUNCTIONS *********************************
// #############################################################################


// #########################################################################
// Print collisional ionisation and recombination LUT to file
void MPv10::print_CIR_LUT(
    const std::vector<double> &X,
    const std::vector<std::vector<double> > &Fx,
    const std::vector<std::vector<double> > &Fx_slope,
    const std::string filename)
{

  std::ofstream outfile(filename + ".txt");
  outfile
      << "# MPv10-" + filename
             + ": column-1 = Temperature(K), other columns = Rates(cm^3/s) \n";
  if (!outfile.is_open()) {
    spdlog::error("MPv10 {} : {}", "couldn't open outfile", 1);
    exit(1);
  }

  outfile.setf(ios_base::scientific);
  outfile.precision(4);

  double x_min      = X.front();
  double x_max      = X.back();
  double delta_logx = log10(X[1]) - log10(X[0]);


  // Make a new set of X array, call it X_new
  int N_points = 500;  // No of new x-points, set it to 500
  double delta = (log10(x_max) - log10(x_min)) / (N_points - 1);
  std::vector<double> X_new;
  for (int i = 0; i < N_points; i++) {
    double logx_new = log10(x_min) + i * delta;
    X_new.push_back(pow(10.0, logx_new));
  }


  for (int k = 0; k < X_new.size(); k++) {
    int x_index = int((log10(X_new[k]) - log10(x_min)) / delta_logx);
    double dx   = X_new[k] - X[x_index];
    outfile << X_new[k];

    for (int i = 0; i < Fx.size(); i++) {
      outfile << "  " << Fx[i][x_index] + dx * Fx_slope[i][x_index];
    }

    outfile << "\n";
  }
}
// End of print collisional ionisation and recombination LUT to file
//*************************************************************************


// #########################################################################
// Print 2D vector to file
void MPv10::print_2D_vector(
    std::vector<std::vector<double> > &vector,
    const std::string ion_name,
    std::string filename)
{
  filename = filename + "_" + ion_name + ".txt";
  ofstream file;
  file.open(filename);

  for (size_t j = 0; j < vector.size(); j++) {
    for (size_t k = 0; k < vector[j].size(); k++) {
      // output is log10 of actual value.
      file << log10(vector[j][k]);
      if (k == vector[j].size() - 1) {
        file << endl;
      }
      else
        file << ",";
    }
  }

  file.close();
}
// End of print 2D vector to file
//*************************************************************************


// #########################################################################
// Print 1D vector to file
void MPv10::print_1D_vector(
    std::vector<double> &vector,
    const std::string variable_name,
    std::string filename)
{
  filename = filename + "_" + variable_name + ".txt";
  ofstream file;
  file.open(filename);

  for (size_t i = 0; i < vector.size(); i++) {
    // output is log10 of actual value.
    file << log10(vector[i]);
    if (i == vector.size() - 1) {
      break;
    }
    else
      file << endl;
  }

  file.close();
}
// End of print 1D vector to file
//*************************************************************************



/*
//###########################################################################
// Print mellema cooling rates to files.
void MPv10::print_to_file_MC_rate(
        int DBLocation  // DataBase Location of the Ion.
        ) {

std::vector<double> logT_points;
logT_points.resize(81);
for (int i = 0; i < 81; i++)
  logT_points[i] =  1.0 + 0.1 * i;

std::vector<double> logne_points;
logne_points.resize(13);
for (int i = 0; i < 13; i++)
  logne_points[i] =  0.5 * i;



std::vector <std::vector<double>> Mellema_Rate;
// Our Original Mellema cooling tabale has the follwing size.
Mellema_Rate.resize(81);
for (int i = 0; i < 81; i++)
  Mellema_Rate[i].resize(13);


double T;
double ne;

cout << "Species Name = " << LUT[DBLocation].Name << endl;

for (int i = 0; i < 81; i++) {
  for (int j = 0; j < 13; j++) {
    T = pow(10, 1.0 + 0.1 * i);
    ne = pow(10, 0.5 * j);
    //cout << "T = " << T << " | " << "ne = " << ne << endl;
    Mellema_Rate[i][j] = mellema_cooling_rate(DBLocation, T, ne);
  }
}

print_2D_vector(Mellema_Rate, LUT[DBLocation].Name,"Mellema_Rate");
print_1D_vector(logT_points, "logT", "Mellema_Rate");
print_1D_vector(logne_points, "logne", "Mellema_Rate");

}
// End of Print mellema cooling rates to files *******************************

*/