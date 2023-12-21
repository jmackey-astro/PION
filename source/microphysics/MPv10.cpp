/// \file MPv10.cpp
/// \authors Arun Mathew, Jonathan Mackey, Maggie Celeste Goulden
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
// Note:
// When iterating over elements, we use the integer variable 'elem'
// as the iterator in loops.
// When iterating over ions, we use 'i' as the iterator in loops.
// When referring to neutral atoms and ions collectively, we use the
// term 'species', and variable 's' as the iterator in loops.
// When working with energy bins, we use the variable 'b' as the
// iterator in the loop.
//
// List of Abbreviations used.
// LUT - Look Up Table
// CLV  - Correct Local Vector
// CIR - collisional ionisation and recombination
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

//#define DTFRAC 0.25
#define DTFRAC 0.01

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



// GET PROBLEM SIZE  ######################################################
void MPv10::get_problem_size(
    int *n_eqn,  ///< number of equations
    int *n_para  ///< number of parameters in user_data vector.
)
{
  *n_eqn  = N_equations;
  *n_para = N_extradata;
  return;
}



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
#ifdef MELLEMA
    mellema_cooling(),
#elif defined CHIANTI
    chianti_cooling(),
#else
#error "must define MELLEMA or CHIANTI"
#endif
    ndim(nd), eos_gamma(g), coord_sys(csys), T_min(1e1), T_max(1e9),
    Num_temps(300), Num_ne(100), ne_min(1.0), ne_max(1.0e6),
    stellar_atmosphere_models(), multi_ion_photo_xsections(), photoionisation()
{

  spdlog::debug(
      "MPv10: Setting up tracer variables, assuming tracers"
      " are last {} variables in primitive vector",
      ntracer);

  // INVOKING MPV10 TRACER SPECIES CATALOGER *****************************
  tracer_species_cataloger(nv, ntr, tracers);
  // This function will generate following indices for MPv10 module.
  spdlog::debug("MPv10: Sorted tracer list {}", MPv10_tracer_list);
  spdlog::debug("MPv10: Element's primitive index {}", elem_prim_index);
  spdlog::debug("MPv10: Element's atomic mass {}", elem_atomic_mass);
  spdlog::debug("MPv10: MPv10 ion list {}", MPv10_ion_list);
  spdlog::debug("MPv10: Ion's element index {}", ions_tracer_elem);
  spdlog::debug("MPv10: Ion-string prim-index map {}", ions_primindex_map);
  spdlog::debug("MPv10: Ion's primitive index {}", ions_prim_index);
  spdlog::debug("MPv10: Ion's local index {}", ions_local_index);
  spdlog::debug("MPv10: Number of ions by elements {}", N_ions_by_elem);
  spdlog::debug("MPv10: Number of species by elements {}", N_species_by_elem);
  spdlog::debug("MPv10: Ion's electron numbers {}", ions_electron_num);
  spdlog::debug("MPv10: Minus ions local index {}", minus_ions_local_index);
  spdlog::debug("MPv10: Colli tracer list = {}", ci_tracer_list);
  spdlog::debug("MPv10: Colli tracer elements = {}", ci_tracer_elem);
  spdlog::debug("MPv10: Recomb tracer list = {}", recomb_tracer_list);
  spdlog::debug("MPv10: Recomb tracer elements = {}", recomb_tracer_elem);


  // SETTING UP LOOKUP TABLE VARIABLES ***********************************
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

  // COLLISIONAL IONISATION AND RECOMBINATION ****************************
  // setting up collisional ionisation
  setup_collisional_ionisation(ci_tracer_list);
  // Generating collisional ionisation rate lookup tables
  generate_collisional_ionisation_LUT(T_table);
  // setting up recombination
  setup_recombination(recomb_tracer_list);
  // Generating recombination rate lookup tables
  generate_recombination_LUT(T_table);

#ifdef MPv10_DEBUG
  mkdir("MPV10-Tables", 0777);
  mkdir("MPV10-Tables/coll-recom-tables", 0777);
  // PRINT TO FILE COLLISIONAL AND RECOMBINATION TABLES
  spdlog::debug("MPv10: writting collisional ionisation rate table into file");
  std::string ci_info =
      "#FILE: Generated by PION - MPv10 Module\n"
      "#DATA DESCRIPTION: Collisional ionisation rate table for the collisional-tracer species.";
  print_CIR_LUT(
      ci_tracer_list, ci_info, T_table, collisional_rate_table,
      collisional_slope_table,
      "MPV10-Tables/coll-recom-tables/coll-ionise-rate");
  spdlog::debug("MPv10: writting recombination rate table into txt file");
  std::string recomb_info =
      "#FILE: Generated by PION - MPv10 Module\n"
      "#DATA DESCRIPTION: Recombination rate table for the recombination-tracer species.";
  print_CIR_LUT(
      recomb_tracer_list, recomb_info, T_table, recombination_rate_table,
      recombination_slope_table, "MPV10-Tables/coll-recom-tables/recomb-rate");
#endif

#ifdef MELLEMA
  // MELLEMA COOLING FUNCTION *********************************************
  setup_mellema_cooling(MPv10_tracer_list);
  spdlog::debug(
      "MPv10: Mellema cooling species locator = {}", mc_tracer_locator);
  generate_mellema_table(T_table, ne_table);
#elif defined CHIANTI
  // CHIANTI COOLING FUNCTION ********************************************
  setup_chianti_cooling(MPv10_tracer_list);
  spdlog::debug("MPv10: chianti species locator = {}", chianti_tracer_locator);
  generate_chianti_table(T_table, ne_table);
#endif

#ifdef MPv10_DEBUG
#ifdef MELLEMA
  mkdir("MPV10-LUTs/cooling-tables", 0777);
  mkdir("MPV10-LUTs/cooling-tables/mellema", 0777);
  spdlog::debug(
      "MPv10: writting Mellema cooling look up tables into txt files");
  std::string cooling_table_info;

  int species_index = 0;
  int location;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      location = mc_tracer_locator[i][j];  // database location
      print_cooling_table(
          ne_table, T_table, mellema_table[species_index].rate,
          MPv10_tracer_list[i][j],
          "MPV10-LUTs/cooling-tables/mellema/mellema_rates");
      species_index++;
    }
  }
#elif defined CHIANTI
  mkdir("MPV10-Tables/cooling-tables", 0777);
  mkdir("MPV10-Tables/cooling-tables/chianti", 0777);
  spdlog::debug("MPv10: writting Chianti cooling lookup tables into txt files");
  std::string cooling_table_info;

  int species_index = 0;
  int location;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      location = chianti_tracer_locator[i][j];  // database location
      print_cooling_table(
          ne_table, T_table, chianti_table[species_index].rate,
          MPv10_tracer_list[i][j],
          "MPV10-Tables/cooling-tables/chianti/chianti_rates");
      species_index++;
    }
  }
#endif
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


  // INITIALIZING IMPORTANT VARIABLES AND VECTORS ************************
  // resizing element number density vector
  elem_number_density.resize(N_elements);
  // resizing the corrector vector
  corrector.resize(N_prim, 1.0);
  // Set up local variables
  setup_local_vectors();
  // Set gamma - 1 from EOS gamma
  gamma_minus_one = eos_gamma - 1.0;
  // Set tolerance
  Min_NeutralFrac = MPv10_ABSTOL;
  Max_NeutralFrac = 1.0 - MPv10_ABSTOL;



  // INITIALIZING CVODE SOLVER *******************************************
  // Initialise the CVODES solver memory etc.
  setup_cvode_solver_without_Jacobian();



  // PHOTO-IONIZATION AND PHOTO-HEATING **********************************
  if (EP->phot_ionisation) {

    // Setup stellar atmosphere models for stellar luminosity spectrum
    setup_stellar_atmosphere_models(ebins, N_ebins, photo_emax);
    // This will set up the following
    spdlog::debug("MPv10: Energy bins {}", ebins);
    spdlog::debug("MPv10: Number of energy bins {}", N_ebins);
    spdlog::debug("MPv10: Maximum photo-ionisation energy {} eV", photo_emax);

    // Invoking photo-ionisation species cataloger
    photo_species_cataloger();
    // This will generate the following indexing
    spdlog::debug("MPv10: Photo tracer list {}", photo_tracer_list);
    spdlog::debug("MPv10: Number of photo species {}", N_photo_species);
    spdlog::debug("MPv10: Photo prim index {}", photo_prim_index);
    spdlog::debug("MPv10: Photo tracer elements {}", photo_tracer_elem);
    spdlog::debug("MPv10: Photo tracer charge {}", photo_tracer_charge);
    spdlog::debug("MPv10: Photo ions local index {}", photo_ions_local_index);
    spdlog::debug(
        "MPv10: Photo species higher ion local index {}",
        photo_species_higher_ion_local_index);

    // Setup photo-ionisation x-section module: multi_ion_photo_xsection class.
    // This will determine the appropriate size for the x-section and bin
    // fraction tables. Additionally, it will return the photo tracer species'
    // ionisation threshold vector.
    photo_tracer_Eth = setup_photo_xsections(photo_tracer_list);
    spdlog::debug("MPv10: ionisation threshold {}", photo_tracer_Eth);

    // Generating mean cross-section and bin fraction table for the energy bins
    generate_ebin_mean_xsection(ebins);

#ifdef MPv10_DEBUG
    mkdir("MPV10-Tables/photo-tables", 0777);
    // print to file photo-ionisation cross-section tables
    spdlog::debug("MPv10: writting photo cross-section table into file");
    const std::string xsection_info = "Mean photo-ionization cross-section"
                                      " for the photo-tracer species.";
    print_to_file_photoionisation(
        xsection_info, photo_tracer_list, ebins, photo_xsection_table,
        "MPV10-Tables/photo-tables/mean-photo-xsection");

    // print to file bin weight tables
    spdlog::debug("MPv10: writting bin weight table into file");
    const std::string weight_info = "Bin weight for the photo-tracer species.";
    print_to_file_photoionisation(
        weight_info, photo_tracer_list, ebins, bin_weight_table,
        "MPV10-Tables/photo-tables/bin-weights");
#endif

    // Identify the radiation sources
    // generate corresponding stellar luminosity per bin, and
    // create corresponding tables for photoionization rates
    // and photoheating rates.
    N_rsrc = 0;
    for (int q = 0; q < RS->Nsources; q++) {

      if (RS->sources[q].type == RT_SRC_SINGLE) {
        N_rsrc++;
        // IMPO: proceed if they are ionising sources.

        // IMPO: why the following line ?
        rt_data.resize(N_rsrc);  // setting the size of the vector rt_data
        // spdlog::debug("RT_data size {}", rt_data.size());

        // Generate the energy flux in different energy bin for the given
        // radiation source.
        int err = generate_bin_luminosity(
            RS->sources[q].effect, RS->sources[q].Tstar, RS->sources[q].Rstar,
            RS->sources[q].strength, bin_luminosity);
        if (err) {
          spdlog::error(
              "Attempt to generate stellar spectrum for source-id: {} "
              "failed {}",
              RS->sources[q].id, bin_luminosity);
          exit(1);
        }
      }
    }
    spdlog::debug("MPv10: {} radiation source(s) detected", N_rsrc);


    // Setup photo-ionisation module: photoionisation
    setup_photoionisation(
        ebins, N_photo_species, photo_xsection_table, bin_weight_table,
        bin_luminosity, photo_tracer_Eth);


  }  // end of photo-ionisation section


  // CONSTRUCTOR END MESSAGE **********************************************
  spdlog::debug("MPv10: Constructor finished and returning");
}



// #######################################################################
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



// MPV10 DESTRUCTOR #########################################################
MPv10::~MPv10()
{
  // Free vector memory
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}



// TRACER SPECIES CATALOGER ################################################
void MPv10::tracer_species_cataloger(
    const int N_primitive, const int N_tracers, const std::string *tracers)
{

  // Beginning of cataloging tracers
  spdlog::info("Cataloging microphysics tracers");

  spdlog::debug(
      "MPv10: Tracers variables {} ",
      std::vector<string>(tracers, tracers + N_tracers));


  // Make a copy of tracer list
  tracer_list_replica.resize(N_tracers);
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
      elem_atomic_mass.push_back(
          atomic_physics_data.get_atomic_mass(tracer_list_replica[s]));
      // No of ions in each element
      N_ions_by_elem.resize(N_elements, 0);  // initializing
      // No of species in each element
      N_species_by_elem.resize(N_elements, 1);  // initializing
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
        ions_electron_num.push_back(
            atomic_physics_data.get_electron_num(tracer_list_replica[s]));
        // Position of ions in the primitive vector.
        ions_prim_index.push_back(first_tracer_index + s);
        // Position of ions in the local vector.
        ions_local_index.push_back(ion_count);
        // Increment the counter when ion of the same element is found.
        ion_count++;
        // No of ions each element
        N_ions_by_elem[e]++;
        // No of species in each element
        N_species_by_elem[e]++;
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
      if (atomic_physics_data.get_charge(MPv10_ion_list[i])
          > atomic_physics_data.get_charge(MPv10_ion_list[i - 1])) {
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

  spdlog::debug(
      "MPv10: Found {} elements and {} ions in primitive vector", N_elements,
      N_sp);

  // Make species list and resize species number density vector **************
  // Collisional ionisation tracer list to generate LUTs
  // Photoionisation tracer list to generate LUTs
  // Recombination to generate LUTs
  // resizing species_number_density 2D-vector

  species_number_density.resize(MPv10_tracer_list.size());  // resizing
  for (int e = 0; e < MPv10_tracer_list.size(); e++) {
    species_number_density[e].resize(MPv10_tracer_list[e].size());  // resizing
    for (int i = 0; i < MPv10_tracer_list[e].size(); i++) {

      // Make collisional/photo ionisation tracer list
      if (i < MPv10_tracer_list[e].size() - 1) {
        ci_tracer_list.push_back(MPv10_tracer_list[e][i]);
        // photo_tracer_list.push_back(MPv10_tracer_list[e][i]);
        ci_tracer_elem.push_back(e);
        // photo_tracer_elem.push_back(e);
      }

      // Make recombination tracer list
      if (i != 0) {
        recomb_tracer_list.push_back(MPv10_tracer_list[e][i]);
        recomb_tracer_elem.push_back(e);
      }
    }
  }
  // End of making species list and resizing species number density **********
}



// PHOTO-IONISATION SPECIES CATALOGER #####################################
void MPv10::photo_species_cataloger()
{
  spdlog::info("Cataloging photo-ionisation tracers");

  N_photo_species = 0;
  for (int e = 0; e < MPv10_tracer_list.size(); e++) {
    for (int i = 0; i < MPv10_tracer_list[e].size(); i++) {

      // Make photo ionisation tracer list
      if (i < MPv10_tracer_list[e].size() - 1) {
        // if the ionisation energy is less than the maximum photon energy,
        // then the tracer is qualified photo species
        if (atomic_physics_data.get_ionisation_energy(MPv10_tracer_list[e][i])
            < photo_emax) {
          // record the photo tracer species in a vector
          photo_tracer_list.push_back(MPv10_tracer_list[e][i]);
          // record photo species's element
          photo_tracer_elem.push_back(e);
          // count number of photo species
          N_photo_species++;

          // find the photo tracer position in primitive vector
          auto it = std::find(
              tracer_list_replica.begin(), tracer_list_replica.end(),
              MPv10_tracer_list[e][i]);
          if (it != tracer_list_replica.end()) {
            int position = static_cast<int>(
                std::distance(tracer_list_replica.begin(), it));
            // record photo species primitive index
            // this could be done by adding the position to first tracer index
            photo_prim_index.push_back(first_tracer_index + position);
            // record the charge of photo species
            photo_tracer_charge.push_back(i);

            // find the photo ion's local index and record -1 otherwise
            if (i == 0) {
              // As neutral species are not present in the local vector,
              // record -1 to indicate 'not found'
              photo_ions_local_index.push_back(-1);
            }
            // else, the local index is obtained by subtracting N_elements
            // from the position of the photo species
            else
              photo_ions_local_index.push_back(position - N_elements);
          }
        }  // End of qualified photo  species
      }
    }  // loop over ions of the specific element
  }    // loop over elements in the tracer list

  // Find the index of next ionised state of photo species in the
  // local vector
  for (int s = 0; s < N_photo_species; s++) {
    std::string element_name =
        atomic_physics_data.get_element_name(photo_tracer_list[s]);
    int charge = atomic_physics_data.get_charge(photo_tracer_list[s]);

    // search for the next higher ionised state of the ion in the
    // MPv10_ion_list.
    // Instead of 's', we use 'j' for the species iteration in this loop.
    for (int j = 0; j < N_sp; j++) {
      if (atomic_physics_data.get_element_name(MPv10_ion_list[j])
              == element_name
          && atomic_physics_data.get_charge(MPv10_ion_list[j]) == charge + 1) {
        photo_species_higher_ion_local_index.push_back(j);
        break;
      }
      // IMPO: what if the next higher ion do not exist in the tracer list.
    }
  }
}



// GET ELECTRON NUMBER DENSITY ################################################
// Calculates electron number density from primitive vector
double MPv10::get_n_elec(const pion_flt *P  ///< primitive state vector array.
)
{
  // get electron number density from P
  int species_counter = 0;
  double ne;
  for (int elem = 0; elem < N_elements; elem++) {  // loop over every element
    int N_elem_species    = N_ions_by_elem[elem];
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


// ##################################################################
double MPv10::get_X_H()
{
  spdlog::error("don't call get_X_H() in MPv10");
  exit_pion(1);
  return 0.0;
}



// ##################################################################
// Returns the index of the ion in the primitive vector
int MPv10::Tr(const string s)
{
  if (ions_primindex_map.find(s) == ions_primindex_map.end())
    return -1;
  else
    return ions_primindex_map[s];
}



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
    int N_elem_species = N_ions_by_elem[e];
    // neutral frac, got by subtracting off the ion fractions in the next loop.
    pion_flt neutral_frac = 1;

    // loop over every species in THIS element
    for (int s = 0; s < N_elem_species; s++) {
      pion_flt number_density =
          element_number_density[e] * y_ion_frac[species_counter];
      int N_electrons = ions_electron_num[species_counter];
      // add on number of particles got by electrons + ionised atom
      n_tot += (1 + N_electrons) * number_density;

      neutral_frac -= y_ion_frac[species_counter];
      species_counter++;
    }
    n_tot += max(0.0, neutral_frac) * element_number_density[e];
  }
  return n_tot;
}


// ##################################################################
// GET OPTICAL DEPTH
// Calculate the optical depth of the current cell for individual energy bins
void MPv10::get_dtau(
    const rad_source *,    ///< pointer to radiation source struct
    const pion_flt ds,     // ds, thickness of the cell
    const pion_flt *p_in,  // input primitive vector
    pion_flt *dtau_vec     // output dtau vector
)
{
  // spdlog::info("get_dtau: starting");
  // 1. Calculating the neutral mass fraction for the current cell
  std::vector<double> neutral_mass_fraction(0);
  double neutral_fraction;
  int sct = 0;
  for (int s = 0; s < photo_tracer_list.size(); ++s) {
    neutral_fraction = 0.0;
    if (photo_tracer_charge[s] == 0) {
      // initial value
      neutral_fraction = p_in[photo_prim_index[s]];
      // substracting the corresponding ion fraction of the same element
      for (int i = 0; i < N_ions_by_elem[photo_tracer_elem[s]]; i++) {
        neutral_fraction -= p_in[ions_prim_index[sct]];
        sct++;
      }
      // if (neutral_fraction < 0.0) {
      //  spdlog::error("***nfrac! s {}, frac {}",s, neutral_fraction);
      //}
      neutral_mass_fraction.push_back(max(0.0, neutral_fraction));
    }
  }

  // 2. Calculate optical depth through cell for each photon energy bin.
  double species_mass_fraction;  // mass fraction of a specific species
  double bin_optical_depth;      // cell optical depth in a specific energy bin
  int neu_species_counter;       // counter for neutral species

  for (int bin = 0; bin < N_ebins; bin++) {
    // loop over energy bin
    // calculate optical depth in bin-th energy bin
    bin_optical_depth   = 0.0;
    neu_species_counter = 0;

    // loop over photo species
    for (int s = 0; s < photo_tracer_list.size(); s++) {
      // mass fraction of each photo species
      species_mass_fraction = 0.0;
      // if the photo species is neutral
      if (photo_tracer_charge[s] == 0) {
        species_mass_fraction = neutral_mass_fraction[neu_species_counter];
        neu_species_counter++;
      }
      // if the photo species is ionised
      else {
        species_mass_fraction = p_in[photo_prim_index[s]];
      }

      // equation: X_{species} sigma_{species, bin} / atomic_mass_{species}
      bin_optical_depth += species_mass_fraction * photo_xsection_table[s][bin]
                           * bin_weight_table[s][bin]
                           / elem_atomic_mass[photo_tracer_elem[s]];
      // spdlog::info("get_dtau: smf {}, psx {}, bwt {}, am"
      //            "{}",species_mass_fraction, photo_xsection_table[s][bin],
      // bin_weight_table[s][bin],  elem_atomic_mass[photo_tracer_elem[s]]);

    }                                    // end of species loop
    bin_optical_depth *= p_in[RO] * ds;  // equation: * rho * ds
    // spdlog::info("get_dtau: bin {}, bin_optical_depth {}", bin,
    // bin_optical_depth);
    dtau_vec[bin] = bin_optical_depth;
  }  // end of loop over energy bin
  // spdlog::info("get_dtau: ending");
}



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

  if (p_in[PG] > 1.0) exit_pion(0);

  // make sure that p_in[] is a self-consistent state:
  // sCMA(corrector, p_in);

  // ==============================================================
  //  Set elemental number density from the current primitive vector
  // ==============================================================
  for (int i = 0; i < N_elements; i++) {  // loop over every element
    elem_number_density[i] =
        p_in[RO]
        * (p_in[elem_prim_index[i]] * corrector[elem_prim_index[i]]
           / elem_atomic_mass[i]);
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
    int N_elem_species = N_ions_by_elem[e];
    // loop over every species in this element
    for (int s = 0; s < N_elem_species; s++) {

      p_local[ions_local_index[species_counter]] =
          p_in[ions_prim_index[species_counter]]
          * corrector[ions_prim_index[species_counter]]
          / p_in[elem_prim_index[e]];

      p_local[ions_local_index[species_counter]] =
          min(1.0, p_local[ions_local_index[species_counter]]);

      p_local[ions_local_index[species_counter]] =
          max(0.0, p_local[ions_local_index[species_counter]]);
      species_counter++;
    }
  }
  // **************************************************************

  // ==============================================================
  // Set internal energy density in local vector.
  // ==============================================================
  p_local[E_index] = p_in[PG] / (gamma_minus_one);
  if (!isfinite(p_local[E_index])) {
    spdlog::error("mpv10: input pressure is not finite: {}", p_in[PG]);
    exit_pion(4);
  }
  // **************************************************************

#ifdef MPv10_DEBUG
  /*
  //
  // Check for NAN/INF
  //
  for (int v = 0; v < 2; v++) {
    if (!isfinite(p_local[v])) {
      spdlog::error("{}: {}", "INF/NAN input to microphysics", p_local[v]);
      exit_pion(1);
    }
    if (mpv_nH < 0.0 || !isfinite(mpv_nH)) {
      spdlog::error(
          "{}: {}", "Bad density input to MPv10::convert_prim2local", mpv_nH);
      exit_pion(1);
    }
  }
   */
#endif  // MPv10_DEBUG

  return 0;
}



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
  y_ion_frac.resize(N_sp);


  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species = N_ions_by_elem[e];
    for (int s = 0; s < N_elem_species;
         s++) {  // loop over every species in THIS element

      p_out[ions_prim_index[species_counter]] =
          p_local[ions_local_index[species_counter]] * p_in[elem_prim_index[e]];

      y_ion_frac[species_counter] = p_local[ions_local_index[species_counter]];
      species_counter++;
    }
  }

  // Set mass fraction tracers to be within the required range (not too close
  // to zero or 1)
  /*
    species_counter = 0;
    for (int e = 0; e < N_elements; e++) {  // loop over every element
      int N_elem_species        = N_ions_by_elem[e];
      p_out[elem_prim_index[e]] = max(
          0.0,
          min(1.0, static_cast<double>(p_out[elem_prim_index[e]])));

      for (int s = 0; s < N_elem_species;
           s++) {  // loop over every species in THIS element


        p_out[ions_prim_index[species_counter]] = max(
            Min_NeutralFrac,
            min(static_cast<double>(p_out[elem_prim_index[e]]) *
    Max_NeutralFrac,
                static_cast<double>(p_out[ions_prim_index[species_counter]])));
        species_counter++;
      }
    }
  */

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

  return 0;
}



// #############################################################################
// SET LOCAL ION FRACTIONS
void MPv10::set_y_ion_frac(
    const std::vector<double> &p_local, std::vector<double> &y_ion_frac)
{
  for (int s = 0; s < N_sp; s++) {
    y_ion_frac[s] = p_local[s];
    // Set y_ion_frac within physically acceptable range.
    y_ion_frac[s] = max(MPv10_ABSTOL, y_ion_frac[s]);
    y_ion_frac[s] = min(1.0 - MPv10_ABSTOL, y_ion_frac[s]);
  }
  return;
}



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
    for (int s = 0; s < N_ions_by_elem[elem]; s++) {
      y_neutral_frac[elem] -= y_ion_frac[sct];
      // spdlog::info("y-ion-frac: {:12.9e}, y-n-frac
      // {:12.3e}",y_ion_frac[sct],y_neutral_frac[elem]);
      sct++;
    }
    // Ensure the neutral frac is within the acceptable range
    y_neutral_frac[elem] = max(MPv10_ABSTOL, y_neutral_frac[elem]);
    y_neutral_frac[elem] = min(1.0 - MPv10_ABSTOL, y_neutral_frac[elem]);
  }
  return;
}



// #############################################################################
// SET ELECTRON NUMBER DENSITY
void MPv10::set_ne(const std::vector<double> &y_ion_frac, double &ne)
{
  // Calculate electron number density
  int sct = 0;
  // loop over every element
  for (int elem = 0; elem < N_elements; elem++) {
    // loop over every species in THIS element
    for (int s = 0; s < N_ions_by_elem[elem]; s++) {
      ne +=
          ions_electron_num[sct] * elem_number_density[elem] * y_ion_frac[sct];
      sct++;
    }
  }
  // Ensure electron number density is not negative.
  if (ne < 0.0) {
    spdlog::info(
        "Warning: Negative ne = {:12.3e}, resetting to {:12.3e}", ne, 1.0e-10);
    ne = max(1.0e-10, ne);
  }

  return;
}



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
  y_ion_frac.resize(N_sp);

  convert_prim2local(pv, p_local, 4);

  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species = N_ions_by_elem[e];
    // loop over every species in THIS element
    for (int s = 0; s < N_elem_species; s++) {
      y_ion_frac[species_counter] =
          pv[ions_prim_index[species_counter]] / pv[elem_prim_index[e]];
      species_counter++;
    }
  }
  return (get_temperature(y_ion_frac, elem_number_density, p_local[E_index]));
}



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
  y_ion_frac.resize(N_sp);

  int err = convert_prim2local(p_pv, p_local, 1);

  // Determine y_ion_frac from the primitive vector
  int species_counter = 0;
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species = N_ions_by_elem[e];
    // loop over every species in THIS element
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
    const int N_ion,  ///< number of ionising radiation sources.
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
  std::vector<double> p_local(N_local, 0.0);
  err = convert_prim2local(p_in, p_local, 2);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to MPv10::TimeUpdateMP()", err);
    exit_pion(1);
  }

  err = correct_localvector(p_local, 1);
  if (err) {
    spdlog::error("Bad input vector to TimeUpdateMP_RTnew");
    exit(1);
  }

  // issue is here, local vector P and NV_Ith_S(y_in, v) are different

  // setup_radiation_source_parameters(p_in, P, ion_src);
  setup_radiation_source_parameters(N_ion, ion_src);

  // Populates CVODE vector with initial conditions (input)
  for (int v = 0; v < N_local; v++) {
    NV_Ith_S(y_in, v) = p_local[v];
    // spdlog::info("v {} y_in {}", v, NV_Ith_S(y_in, v));
  }

  //============
  // Calculate y-dot[] to see if anything is changing significantly over dt
  // Here y_out is y_dot calculated from the ydot function.
  // err = ydot(0, y_in, y_out, 0);
  // if (err) {
  //  spdlog::error(
  //      "{}: {}", "dYdt() returned an error in MPv10::TimeUpdateMP_RTnew()",
  //      err);
  //  exit(1);
  //}
  //============

  err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
  if (err) {
    spdlog::error(
        "integration failed: MPv10::TimeUpdateMP_RTnew() error {}", err);
    for (int v = 0; v < N_local; v++)
      p_local[v] = NV_Ith_S(y_in, v);
    spdlog::error("y_in: {}", p_local);
    for (int v = 0; v < N_local; v++)
      p_local[v] = NV_Ith_S(y_out, v);
    spdlog::error("y_out: {}", p_local);
    exit(1);
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v = 0; v < N_local; v++)
    p_local[v] = NV_Ith_S(y_out, v);

  err = correct_localvector(p_local, 2);
  if (err) {
    spdlog::error("bad output vector from TimeUpdateMP_RTnew");
    exit(2);
  }
  err = convert_local2prim(p_local, p_in, p_out, 2);


#ifdef TEST_INF
  for (int v = 0; v < N_prim; v++) {
    if (!isfinite(p_in[v]) || !isfinite(p_out[v])) {
      spdlog::debug("NAN in MPv3 update: {}", v);
      spdlog::debug("Pin  : {}", std::vector<double>(p_in, p_in + N_prim));
      spdlog::debug("Pout : {}", std::vector<double>(p_out, p_out + N_prim));
      spdlog::debug("Ploc  : {}", p_local);
      // spdlog::error("{}: {}", "NAN in MPv3",P[2]);
      exit_pion(v);
    }
  }
#endif

  return err;
}



// CORRECT LOCAL VECTOR (CLV) #####################################
int MPv10::correct_localvector(
    std::vector<double> &p_local,  ///< check the local vector
    int function_flag              /// flags the function
)
{
  // function_flag 1 -> towards begining of TimeUpdateMP_RTnew
  // function_flag 2 -> towards the end of TimeUpdateMP_RTnew
  // function_flag 3 -> timescales_RT
  // function_flag 4 -> ydot

  int flag = 0;

  // *********************************************************************
  // Check if y_ion_fracs are within the acceptable range
  std::vector<double> y_ion_frac(N_sp);

  for (int i = 0; i < N_sp; i++) {  // loop over ions in p_local
    y_ion_frac[i] = p_local[i];
    // Set y_ion_frac within physically acceptable range.
    if (y_ion_frac[i] < MPv10_ABSTOL) {
      // spdlog::warn("correct localvec {}: negative ion frac [{}] = {:12.3e}",
      //              function_flag,i,y_ion_frac[i]);
      y_ion_frac[i] = MPv10_ABSTOL;  // resetting
    }
    if (y_ion_frac[i] > 1.0 - MPv10_ABSTOL) {
      // spdlog::warn("correct localvec {}: too large ion frac [{}] by
      // {:12.3e}",
      //              function_flag,i,y_ion_frac[i]-1);
      y_ion_frac[i] = 1.0 - MPv10_ABSTOL;  // resetting
    }
  }
  // *********************************************************************


  // *********************************************************************
  // Scaling the y_ion_frac if the total elemental fraction exceed 1.

  // 1. Calculating total ion fraction for each element.
  std::vector<double> y_total_ion_frac(N_elements, 0.0);
  int local_index = 0;
  for (int elem = 0; elem < N_elements; elem++) {
    for (int i = 0; i < N_ions_by_elem[elem]; i++) {
      y_total_ion_frac[elem] += y_ion_frac[local_index];
      local_index++;
    }
  }

  // 2. Rescale the y_ion_frac if total ion fraction is > 1
  local_index = 0;
  for (int elem = 0; elem < N_elements; elem++) {
    // spdlog::info("correct localvec {}: checking el {}",function_flag, elem);
    if (1.0 - y_total_ion_frac[elem] < MPv10_ABSTOL) {
      // if total ion fraction is nearly one or more than one,
      // the neutral fraction must be minumum (say, MPv10_ABSTOL)
      if (y_total_ion_frac[elem] - 1.0 > 1.0e4 * MPv10_RELTOL
          && function_flag != 4) {
        spdlog::warn(
            "correct localvec {}: too many ions! el {}, y-1 = {:12.3e}",
            function_flag, elem, y_total_ion_frac[elem] - 1);
      }
      // double total_fraction = y_total_ion_frac[elem] + 1.0e-30;
      for (int i = 0; i < N_ions_by_elem[elem]; i++) {
        y_ion_frac[local_index] /=
            y_total_ion_frac[elem] + MPv10_ABSTOL;  // total_fraction;
        local_index++;
      }
    }
    else {
      local_index += N_ions_by_elem[elem];
    }
  }
  // **************************************************************

  // **************************************************************
  // Feeding the corrected y_ion_frac back in to p_local
  for (int i = 0; i < N_sp; i++)
    p_local[i] = y_ion_frac[i];
  // **************************************************************

  // **************************************************************
  // Check for negative pressure
  // Note: This shouldn't happen, so we output a warning) and set to
  // 10K if we find it.
  if (p_local[E_index] <= 0.0) {
    if (function_flag != 4)
      spdlog::warn(
          "MPv10::CLV flag: {} - Negative pressure input: e = {}, "
          "setting to {} K",
          function_flag, p_local[E_index], EP->MinTemperature);

    // reset the internal energy (requires using y_ion_frac in get_ntot)
    p_local[E_index] = get_ntot(y_ion_frac, elem_number_density) * pconst.kB()
                       * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************


  // **************************************************************
  // if the temperature is below Minimum Temperature, set it back to
  // Minimum Temperature by set the internal energy accordingly
  double T = get_temperature(y_ion_frac, elem_number_density, p_local[E_index]);

  if (T < 0.95 * EP->MinTemperature) {
    if (function_flag != 4)
      spdlog::warn(
          "MPv10::CLV flag: {} - Temperature = {} K, below T_MIN = {} K",
          function_flag, T, EP->MinTemperature);

    // Resetting the internal energy with Minimum Temperature
    p_local[E_index] = get_ntot(y_ion_frac, elem_number_density) * pconst.kB()
                       * EP->MinTemperature / (gamma_minus_one);
  }
  // **************************************************************

  return flag;
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

  // spdlog::info("got past blockage 0");
  // Note: Minimum value of micro-physics trace variables in the primitive
  // vector is set to 1e-12.
  int err = 0;
  std::vector<double> p_local;
  p_local.resize(N_local, 0.0);

  // First convert to local variables.
  err = convert_prim2local(p_in, p_local, 3);
  if (err) {
    spdlog::error("{}: {}", "Bad input state to MPv10::timescales_RT()", err);
    exit_pion(1);
  }

  // correct local vector
  err = correct_localvector(p_local, 3);
  // spdlog::debug("after correcting - p_local = {}", p_local);
  if (err) {
    spdlog::error("bad input vector to timescales_RT");
    exit_pion(1);
  }

  // v runs over all ions and internal energy
  for (int v = 0; v < N_local; v++) {
    NV_Ith_S(y_in, v) = p_local[v];
  }

  //
  // Next set the radiation properties of the current cell.
  //
  setup_radiation_source_parameters(N_ion, ion_src);
  // spdlog::info("got past blockage 1");
  // spdlog::info("p_in {}",p_local);
  // spdlog::info("E = {}",p_local[E_index]);

  // Now calculate y-dot[]...
  // ydot function returns RHS of ydot equation. Here y_out is RHS of ydot.
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in MPv10::timescales_RT()", err);
    exit_pion(1);
  }
  // spdlog::info("got past blockage 2");

  double tt = Temperature(p_in, 0);
  // spdlog::info("{:12.6e}  {:12.6e}", tt, NV_Ith_S(y_out,E_index));
  // spdlog::info("got past blockage 3");

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
  Y_dot.resize(N_local, 0.0);
  // Since y_out is ydot.
  for (int v = 0; v < N_local; v++)
    Y_dot[v] = NV_Ith_S(y_out, v);
  // spdlog::info("got past blockage 4");

  for (int v = 0; v < N_equations; v++) {
    // t = min(t, DTFRAC / (fabs(NV_Ith_S(y_out, v)) + TINYVALUE));
    t = min(t, DTFRAC / (fabs(Y_dot[v]) + TINYVALUE));
    /*
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
    */
  }
  // spdlog::info("got past blockage N");

  return t;
}



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
    int N_elem_species = N_ions_by_elem[e];
    total_mass_frac += p_in[elem_prim_index[e]];
  }
  double e_correction = 1.0 / total_mass_frac;
  species_counter     = 0;
  // apply all-element correction, calculate species correction, apply species
  // correction
  for (int e = 0; e < N_elements; e++) {  // loop over every element
    int N_elem_species            = N_ions_by_elem[e];
    corrector[elem_prim_index[e]] = e_correction;
    // Calculate all-species-pr-element correction, if needed, i.e.
    double s_frac = 0;

    for (int s = 0; s < N_elem_species; s++) {
      s_frac += p_in[ions_prim_index[species_counter]];
      species_counter++;
    }

    if (s_frac
        > ((p_in[elem_prim_index[e]] * e_correction) - Min_NeutralFrac)) {
      print_flagg = 1;
      double s_correction =
          ((p_in[elem_prim_index[e]] * e_correction) - Min_NeutralFrac)
          / s_frac;
      int inner_species_counter = (species_counter - N_elem_species);
      for (int s = 0; s < N_elem_species; s++) {
        corrector[ions_prim_index[inner_species_counter]] = s_correction;
      }
    }
  }
  return;
}



// #######################################################################
// YDOT FUNCTION
// This function calculate the RHS of ydot equations.
int MPv10::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double *         ///< extra user-data vector (UNUSED)
)
{

  // *********************************************************************
  // SECTION: DETERMINIE THE FOLLOWING FROM THE CURRENT LOCAL VECTOR
  // 1. ne (electron number density)
  // 2. y_ion_frac (ion fractions in the local vector)
  // 3. neutral fraction (neutral fractions in the local vector)
  std::vector<double> y_ion_frac;
  y_ion_frac.resize(N_sp, 0.0);
  std::vector<double> y_neutral_frac;
  y_neutral_frac.resize(N_elements, 0.0);
  double ne = 0.0;

  std::vector<double> p_local;
  p_local.resize(N_local, 0.0);
  // make a copy of y_now into p_local
  for (int v = 0; v < N_local; v++) {
    p_local[v] = NV_Ith_S(y_now, v);
    // spdlog::info("v {} ynow {}", v, NV_Ith_S(y_now, v));
  }

  correct_localvector(p_local, 4);
  // Set local quantities: y_ion_frac, y_neutral_frac and ne
  set_localvariables(p_local, y_ion_frac, y_neutral_frac, ne);
  // *********************************************************************


#ifdef MPv10_CIE_TEST
  // Restricting ne value not less than 1e-4 to perform CIE test.
  ne = max(1e-2, ne);
#endif

  // Initialise all elements of ydot vector to zero.
  for (int v = 0; v < N_equations; v++) {
    NV_Ith_S(y_dot, v) = 0.0;
  }


  // Initializing Edot to zero
  double Edot = 0.0;

  // Get internal energy from the current local vector (y_now)
  // double E_in = NV_Ith_S(y_now, E_index);
  double E_in = p_local[E_index];

  // Calculate temperature from E_in, y_ion_frac and number density.
  double T = get_temperature(y_ion_frac, elem_number_density, E_in);
  // cout<< "ydot: " << "ne = " << ne << " | "<< "T = " << T << " | " <<endl;

  // Restricting temperature within T_min and T_max
  if (T > T_max)
    T = T_max;
  else if (T < T_min)
    T = T_min;

  // spdlog::info("y: {} : {}",y_ion_frac, y_neutral_frac);

  // *********************************************************************
  // SECTION: SPECIES NUMBER DENSITY
  int ion_index = 0;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    double excess_elem_frac = 0.0;
    std::vector<double> species_frac;
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      // Calculating the current number density of neutral atoms.
      if (j == 0) {
        species_number_density[i][j] =
            y_neutral_frac[i] * elem_number_density[i];
        excess_elem_frac += y_neutral_frac[i];
        species_frac.push_back(y_neutral_frac[i]);
      }
      // Calculating the current number density of ions.
      else {
        species_number_density[i][j] =
            y_ion_frac[ion_index] * elem_number_density[i];
        excess_elem_frac += y_ion_frac[ion_index];
        species_frac.push_back(y_ion_frac[ion_index]);
        ion_index++;
      }
    }
    if (fabs(excess_elem_frac - 1.0) > 1e-7) {
      spdlog::error(
          "species counting is wrong {}: species {}, \n\t\t values {}",
          excess_elem_frac, MPv10_tracer_list[i], species_frac);
      if (fabs(excess_elem_frac - 1.0) > 1e-5) {
        exit(1);
      }
    }
  }
  // END OF SECTION: SPECIES NUMBER DENSITY ******************************



  // SECTION: Collisional Ionisation *************************************
  if (EP->coll_ionisation) {
    // Calculating RHS terms corresponding to collisional ionisation for
    // species in ci_tracer_list. These terms are calculated only once and
    // attached at relevant ydot equations.
    //
    // For example, for species s, if ion, the term ci_rate(s)y(s)ne is
    // subtracted for ydot(s) while the same terms is added for the ydot(s+1)
    // equation, provided the specie s+1 exist. Whereas if the species s is
    // neutral, then only addition is performed.

    double ci_rate    = 0.0;
    double this_y_dot = 0.0;

    for (int s = 0; s < ci_tracer_list.size(); s++)  // loop over ci species
    {
      ci_rate = collisional_ionisation_rate(s, T);

      if (minus_ions_local_index[s] != -1) {
        // if the less ionised species is not neutral
        // this_y_dot = ci_rate * NV_Ith_S(y_now, minus_ions_local_index[s]) *
        // ne;
        this_y_dot = ci_rate * p_local[minus_ions_local_index[s]] * ne;
        NV_Ith_S(y_dot, minus_ions_local_index[s]) -= this_y_dot;
      }
      else {
        // if the less ionised species is neutral
        this_y_dot = ci_rate * y_neutral_frac[ci_tracer_elem[s]] * ne;
      }

      NV_Ith_S(y_dot, ions_local_index[s]) += this_y_dot;

      // Cooling due to collisional ionisation of this species.
      Edot -=
          ci_ion_pot[s] * this_y_dot * elem_number_density[ci_tracer_elem[s]];
    }
  }
  // END OF SECTION: Collisional Ionisation ******************************


  // SECTION: Recombination (Radiative + Dielectronic) *******************
  if (EP->recombination) {
    // RHS terms corresponding to Recombination (R) are calculated for
    // all species. In this case, for specie s, the term rec_rate(s)y(s)ne is
    // substrated for ydot(s-1) and the same terms is added for the ydot(s),
    // provided species s is a ion of same element. Whereas
    // if the species s is neutral, then only addition performed.

    double recomb_rate = 0.0;
    double this_y_dot  = 0.0;
    // std::vector<double> temp_rr(recomb_tracer_list.size());  //DEBUG

    for (int s = 0; s < recomb_tracer_list.size(); s++) {
      recomb_rate = recombination_rate(s, T);
      // spdlog::info("rec: {}, i={}, i-1={}", s, ions_local_index[s],
      //              minus_ions_local_index[s]);

      // this_y_dot = recomb_rate * NV_Ith_S(y_now, ions_local_index[s]) * ne;
      this_y_dot = recomb_rate * p_local[ions_local_index[s]] * ne;

      // Subtract this term to the current species equation
      NV_Ith_S(y_dot, ions_local_index[s]) -= this_y_dot;
      // add this term to less ionised species equation provided it
      // is not neutral species
      if (minus_ions_local_index[s] != -1)
        NV_Ith_S(y_dot, minus_ions_local_index[s]) += this_y_dot;

      // Cooling due to recombination of this species.
      Edot -= (3. / 2.) * T * pconst.kB() * this_y_dot
              * elem_number_density[recomb_tracer_elem[s]];
      // if (T>1e6)
      //  spdlog::info("rec-rate {:12.3e} , y {:12.9e}, ne {:12.3e}, cooling
      //  {:12.3e}", this_y_dot, NV_Ith_S(y_now, ions_local_index[s]), ne,(3.
      //  / 2.) * T * pconst.kB() * this_y_dot *
      //  elem_number_density[recomb_tracer_elem[s]]);
      // temp_rr[s] = this_y_dot;
    }
    // spdlog::info("rr: {}", temp_rr);
  }

  // END OF SECTION: Recombination **************************************


  // SECTION: Photoionisation *******************************************
  if (EP->phot_ionisation) {
    // Calculating RHS terms corresponding to photoionisation for species in
    // photo_tracer_list. These terms are calculated only once and attached
    // at relevant ydot equations.
    //
    // For example, for species s, if ion, the term pi_rate(s)y(s) is
    // subtracted for ydot(s) while the same terms is added for the
    // ydot(s+1) equation.
    // Whereas if the species s is neutral, then pi_rate(s)y(s) is added to
    // ydot(s+1) equation.
    int elem_index             = 0;
    int charge_index           = 0;
    double this_y_dot          = 0.0;
    double pi_rate             = 0.0;
    double species_num_density = 0.0;
    double element_num_density = 0.0;

    // std::vector<double> temp_pir(N_photo_species);  //DEBUG

    for (int s = 0; s < N_photo_species; s++)  // loop over pi species
    {
      // get element and charge index
      elem_index          = photo_tracer_elem[s];    // element index
      charge_index        = photo_tracer_charge[s];  // charge index
      species_num_density = species_number_density[elem_index][charge_index];
      element_num_density = elem_number_density[elem_index];
      // if (T>1.0e6)
      //  spdlog::info("{} {:8.3e} {:8.3e}",
      //  MPv10_tracer_list[elem_index][charge_index], species_num_density,
      //  element_num_density);

      // get photo-ionisation rate for the species s
      pi_rate =
          photoionisation_rate(s, species_num_density, element_num_density);

      // spdlog::info("index ={}, n_s ={}, pi_rate = {}", s,
      // species_num_density, pi_rate);

      // if the species is neutral
      if (photo_ions_local_index[s] == -1) {
        // add pi_s * y^{neu}_s to the first ionised state
        this_y_dot = pi_rate;
        // if (T>1e6)
        //  spdlog::info("index ={}, photo_species_higher_ion_local_index[s] =
        //  {}, {:12.3e}", s,
        //  photo_species_higher_ion_local_index[s],this_y_dot);
        NV_Ith_S(y_dot, photo_species_higher_ion_local_index[s]) += this_y_dot;
      }
      // otherwise
      else {
        this_y_dot = pi_rate;
        // if (T>1e6)
        //  spdlog::info("index ={}, photo_ions_local_index[s] = {}, {:12.3e}",
        //  s, photo_ions_local_index[s],this_y_dot);

        NV_Ith_S(y_dot, photo_ions_local_index[s]) -= this_y_dot;
        // IMPO: what if the next higher ion do not exist in the tracer list,
        // then we should skip next line.
        // spdlog::info("index ={}, photo_species_higher_ion_local_index[s] =
        // {}", s, photo_species_higher_ion_local_index[s]);
        NV_Ith_S(y_dot, photo_species_higher_ion_local_index[s]) += this_y_dot;
      }
      // Heating due to photo-ionisation of this species.

      Edot += photoheating_rate(s, species_num_density);

      // temp_pir[s] = pi_rate;
      // spdlog::info("{}",MPv10_tracer_list[elem_index][charge_index]);
      // if (MPv10_tracer_list[elem_index][charge_index] == "H")
      //  spdlog::info("{} {:8.3e} {:8.3e} {:8.3e} {:8.3e}",
      //    MPv10_tracer_list[elem_index][charge_index], pi_rate,
      //    number_density, photoheating_rate(s, number_density) / (pi_rate *
      //    number_density * pconst.eV()), T);
    }
    // spdlog::info("pir: {}",temp_pir);
  }
  // END OF SECTION: Photoionisation *************************************



  // SECTION: Cooling Function ********************************************
  // Edot source term associated with Cooling is calculated in this Section.
  double Lambda         = 0.0;
  double L              = 0.0;
  double number_density = 0.0;

#ifdef MPv10_CIE_COOLING
  mc_outfile.open("cooling_function.txt", ios::app);
  mc_outfile << T << "  ";
#endif

  int species_index = 0;
  for (int i = 0; i < MPv10_tracer_list.size(); i++) {
    for (int j = 0; j < MPv10_tracer_list[i].size(); j++) {
      number_density = species_number_density[i][j];
#ifdef MELLEMA
      // mellema cooling rate
      L = mellema_cooling_rate(species_index, T, ne);
#elif defined CHIANTI
      // chianti cooling rate
      L = chianti_cooling_rate(species_index, T, ne);
#endif
      Lambda += ne * number_density * L;
      // if (L*number_density > 1.0e-22)
      //  spdlog::info("{} {}  {:8.3e} {:8.3e} {:8.3e}",
      //    MPv10_tracer_list[i][j], species_index, L*number_density, T, ne);
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
  // END OF SECTION: Cooling Function **************************************

  // spdlog::debug("Edot = {}", Edot);

  //
  // We want to limit cooling as we approach the minimum temperature, so we
  // scale the rate to linearly approach zero as we reach T_min.
  //
  if (Edot < 0.0 && T < 2.0 * EP->MinTemperature) {
    Edot = min(0.0, (Edot) * (T - EP->MinTemperature) / EP->MinTemperature);
  }

  if (!EP->update_erg) Edot = 0.0;
  NV_Ith_S(y_dot, E_index) = Edot;

  for (int v = 0; v < N_sp; v++) {
    if (y_ion_frac[v] < MPv10_ABSTOL && NV_Ith_S(y_dot, v) < 0.0)
      NV_Ith_S(y_dot, v) = 0.0;
    if (y_ion_frac[v] > 1.0 - MPv10_ABSTOL && NV_Ith_S(y_dot, v) > 0.0)
      NV_Ith_S(y_dot, v) = 0.0;
  }

  // for (int v = 0; v < N_local; v++) p_local[v] = NV_Ith_S(y_dot, v);
  // spdlog::info("ydot: {}",p_local);

  return 0;
}
// END OF YDOT FUNCTION ####################################################



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



// GET RECOMBINATION RATE ###################################################
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

// GET PHOTO-IONISATION RATE  ###############################################
double MPv10::photoionisation_rate(
    const int species_index,           // photo species index (first index)
    const double species_num_density,  // species number density
    const double element_num_density   // species element number density
)
{
  // spdlog::info("index = {}, n_H = {}, n_s = {}", species_index,
  // element_num_density, species_num_density);
  return calculate_photoionisation_rate(
      species_index, species_num_density, element_num_density, vshell, ds,
      tau_bins);
}

// GET PHOTO-HEATING RATE  ###############################################
double MPv10::photoheating_rate(
    const int species_index,          // photo species index (first index)
    const double species_num_density  // species number density
)
{
  return calculate_photoheating_rate(
      species_index, species_num_density, vshell, ds, tau_bins);
}



#ifdef MELLEMA
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
  L = mellema_table[location].rate[T_index][ne_index]
      + dT * mellema_table[location].T_rateslope[T_index][ne_index]
      + dne * mellema_table[location].ne_rateslope[T_index][ne_index];
  return L;
}
#endif


#ifdef CHIANTI
// GET CHIANTI COOLING RATE ###################################################
double MPv10::chianti_cooling_rate(
    const int location,  // species identifier (database location)
    double T,            // Temperature
    double ne            // electron number density
)
{

  // 1. Set up log10(T) and log10(ne) indices to access LUT.
  // A conditional statement is added to access the chinati cooling table only
  // if ne and T range lie within the table in the database of
  // Mellema_/// rate.cpp

  if (T > T_max) T = T_max;
  if (T < T_min) T = T_min;
  if (ne < ne_min) ne = ne_min;

  // Note: Since the chianti cooling rates are given in log scales,
  // i.e, L(Log T, Log ne), we calculate the differential change dLog_T
  // (step) and dLog_ne along with their indices. Delta_Log_T and
  // Delta_Log_ne is a constant step size calculated in the
  // chianti_cooling_rate.cpp

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

  // Calculating chianti cooling rate.
  double L = 0.0;
  // Calculating cooling from the LookUp Table (LUT) with the location of
  // each Ion and atom in the database.
  L = chianti_table[location].rate[T_index][ne_index]
      + dT * chianti_table[location].T_rateslope[T_index][ne_index]
      + dne * chianti_table[location].ne_rateslope[T_index][ne_index];

  return L;
}
#endif


// SET RADIATION SOURCE PARAMETERS ############################################
void MPv10::setup_radiation_source_parameters(
    // const pion_flt *p_in,  ///< primitive input state vector.
    // vector<double> &P,     ///< local input state vector (x_in,E_int)
    const int N_ion,
    const std::vector<struct rt_source_data> &ion_src
    ///< list of ionising src column densities and source properties.
)
{

  if (ion_src.size() != static_cast<unsigned int>(N_ion)) {
    spdlog::error(
        "{}: {}",
        "Timescales: N_ionising_srcs doesn't match vector size in MP3",
        ion_src.size());
    exit(1);
  }

  // necessary? size of energy bin in parameter file should match the
  // energy bin hardcoded

  // proceed if there are any ionizing sources present.
  if (N_ion > 0) {
    // IMPO: this is hard coded for a single source.
    // path length through the current cell
    ds = ion_src[0].dS;
    // shell volume in the current cell
    vshell = ion_src[0].Vshell;
    // optical depth in each energy bin at the edge of the cell
    tau_bins.resize(N_ebins, 0.0);
    for (int b = 0; b < N_ebins; b++)
      tau_bins[b] = ion_src[0].Column[b];
  }
  // spdlog::debug("ds {:12.6e}, Vsh {:12.6e}, col {}", ds, vshell, tau_bins);
}



// ############################################################################
// ***************************** END OF MPv10 *********************************
// ############################################################################



// ############################################################################
//  DEBUG  FUNCTIONS

// ############################################################################
// Print collisional ionisation and recombination LUT to file
void MPv10::print_CIR_LUT(
    const std::vector<string> &tracer_list,
    const std::string info,
    const std::vector<double> &X,
    const std::vector<std::vector<double> > &Fx,
    const std::vector<std::vector<double> > &Fx_slope,
    const std::string filename)
{

  std::ofstream outfile(filename + ".txt");
  outfile << info + "\n";
  // Column names
  outfile << "#ATTRIBUTE LABEL: Temperature";
  for (int j = 0; j < tracer_list.size(); j++) {
    outfile << "\t" << tracer_list[j];
  }
  outfile << "\n";
  // column units
  outfile << "#UNIT: K";
  for (int j = 0; j < tracer_list.size(); j++) {
    outfile << "\t"
            << "cm^3/s";
  }
  outfile << "\n";
  outfile << "#DATA: \n";

  if (!outfile.is_open()) {
    spdlog::error("MPv10 {} : {}", "couldn't open outfile", 1);
    exit_pion(1);
  }

  outfile.setf(ios_base::scientific);
  outfile.precision(4);

  double x_min      = X.front();
  double x_max      = X.back();
  double delta_logx = log10(X[1]) - log10(X[0]);

  // Make a new set of X array, call it X_new
  int N_points = 200;  // No of new x-points, set it to 200
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

// ############################################################################
// Use the following function to print photo-ionisation mean cross-section
// and bin_fraction to a txt file
void MPv10::print_to_file_photoionisation(
    const std::string info,
    const std::vector<string> &species_list,
    const std::vector<std::vector<double> > &X,
    const std::vector<std::vector<double> > &Fx,
    const std::string filename)
{
  std::ofstream outfile(filename + ".txt");
  if (!outfile.is_open()) {
    spdlog::error("MPv10 {} : {}", "couldn't open outfile", 1);
    exit(1);
  }

  outfile << "#FILE: Generated by PION - MPv10 Module \n";
  outfile << "#DATA DESCRIPTION: " + info + "\n";
  // Column names
  outfile << "#ATTRIBUTE LABEL: Bin_Min\tBin_Max";
  for (int j = 0; j < species_list.size(); j++) {
    outfile << "\t" << species_list[j];
  }
  outfile << "\n";
  // column units
  outfile << "#UNIT: eV\teV";
  size_t found = info.find("cross-section");
  if (found != std::string::npos) {
    for (int j = 0; j < species_list.size(); j++)
      outfile << "\t"
              << "cm^2";
  }
  else {
    for (int j = 0; j < species_list.size(); j++)
      outfile << "\t"
              << "1";
  }
  outfile << "\n";


  outfile.setf(ios_base::scientific);
  outfile.precision(4);
  outfile << "#DATA: \n";
  for (int k = 0; k < X.size(); k++) {
    outfile << X[k][0] << "\t" << X[k][1];
    for (int j = 0; j < Fx.size(); j++) {
      outfile << "\t" << Fx[j][k];
    }
    outfile << "\n";
  }
}

// ###########################################################################
// Print cooling tables to file
void MPv10::print_cooling_table(
    std::vector<double> &logne,  // log10(ne)
    std::vector<double> &logT,   // log10(T)
    std::vector<std::vector<double> > &cooling_rate,
    const std::string ion_name,
    std::string filename)
{

  std::ofstream outfile(filename + "_" + ion_name + ".txt");
  if (!outfile.is_open()) {
    spdlog::error("MPv10 {} : {}", "couldn't open outfile", 1);
    exit(1);
  }

  outfile << "#FILE: Generated by PION - MPv10 Module \n";
  outfile << "#DATA DESCRIPTION: " + ion_name + " cooling rate table\n";
  outfile << "#ATTRIBUTE LABEL: Temperature";  // Column names

  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t"
            << "Rate";
  }
  outfile << endl;
  outfile << "#UNIT: K";
  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t"
            << "ergs cm^3/s";
  }
  outfile << endl;

  outfile << "#DATA: \n";
  outfile << "#Log(ne): ";
  for (int j = 0; j < logne.size(); j++) {
    outfile << "\t" << log10(logne[j]);
  }
  outfile << endl;

  for (size_t j = 0; j < cooling_rate.size(); j++) {
    outfile << log10(logT[j]) << ",";
    for (size_t k = 0; k < cooling_rate[j].size(); k++) {
      outfile << log10(cooling_rate[j][k]);
      if (k == cooling_rate[j].size() - 1) {
        outfile << endl;
      }
      else
        outfile << ",";
    }
  }

  outfile.close();
}
