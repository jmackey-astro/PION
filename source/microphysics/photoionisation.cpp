/**
 *  \file     photoionisation.cpp
 *  \author   Arun Mathew
 *  \date     03-07-2023
 **/
#include "photoionisation.h"


photoionisation::photoionisation() {}
photoionisation::~photoionisation() {}


// SETUP PHOTO-IONISATION ################################################
void photoionisation::setup_photoionisation(
    const std::vector<std::vector<double> > &ebins,  // energy bins
    const int N_photo_species,                       // number of photo species
    const std::vector<std::vector<double> > &photo_xsection_table,  //
    const std::vector<std::vector<double> >
        &bin_weight_table,                       // bin weight table
    const std::vector<double> &bin_luminosity,   // stellar spectrum
    const std::vector<double> &photo_tracer_Eth  // photoionisation threshold
)
{

  spdlog::info("Setting up photoionisation module");

  // replicating energy bin
  E_bins = ebins;
  // get the bin size
  bin_size = ebins.size();
  // make a replica of the following
  xsection_table   = photo_xsection_table;
  bin_frac_table   = bin_weight_table;
  bin_spectrum     = bin_luminosity;
  ionisation_thres = photo_tracer_Eth;

  // 1. Making the species threshold bin index table
  // This table contains the index of the bin where the threshold energy
  // lies in the energy bin vector for each photo-species.
  species_threshold_bin_index.resize(N_photo_species, -1);

  int b;
  for (int s = 0; s < N_photo_species; s++) {
    for (b = 0; b < bin_size; b++) {
      if (bin_frac_table[s][b] > 0.0) break;
    }
    species_threshold_bin_index[s] = b;
  }
  spdlog::debug(
      "photoionisation: species_threshold_bin_pos {}",
      species_threshold_bin_index);

  // Verifying the correctness of the bin indexing.
  auto it = std::find(
      species_threshold_bin_index.begin(), species_threshold_bin_index.end(),
      -1);
  if (it != species_threshold_bin_index.end()) {
    // If any entry in the species bin index table is zero, the
    // indexing has failed. It is necessary to reset the energy
    // bin properly.
    spdlog::error("photoionisation: species bin indexing failed.");
    exit(1);
  }
  else {
    spdlog::debug("photoionisation: species bin indexing successful");
  }

  // 2. Making the mean energy table
  mean_energy_table.resize(N_photo_species);
  double mean_bin_energy;
  for (int s = 0; s < N_photo_species; s++) {
    for (int b = 0; b < bin_size; b++) {

      // if the bin weight is zero then skip the iteration with zero mean value
      if (bin_frac_table[s][b] == 0.0) {
        mean_energy_table[s].push_back(0.0);
        continue;
      }
      // if the bin weight is <1, then set the e_bin_min as ionisation threshold
      // of the specific species while calculating the mean bin energy
      if (bin_frac_table[s][b] < 1.0) {
        E_bins[b][0]    = ionisation_thres[s];
        mean_bin_energy = 0.5 * (E_bins[b][0] + E_bins[b][1]);
      }
      // otherwise, that is, if bin weight is 1, then use the following
      // calculation to estimate mean bin energy
      else
        mean_bin_energy = 0.5 * (E_bins[b][0] + E_bins[b][1]);
      // append the mean bin energy to the table for the species
      mean_energy_table[s].push_back(mean_bin_energy * pconst.eV());
    }  // bin loop
  }    // species loop
}


// CALCULATING PHOTO-IONISATION RATE #####################################
double photoionisation::calculate_photoionisation_rate(
    const int species,                   // Species_identifier
    const double species_num_density,    // Species number density,
    const double element_num_density,    // Element number density,
    const double vshell,                 // vshell
    const double ds,                     // ds
    const std::vector<double> &tau_bins  // Optical depth in each energy bin
)
{
  double rate              = 0.0;
  double bin_rate          = 0.0;
  double delta_species_tau = 0.0;

  // spdlog::debug("species = {}", species);

  // spdlog::debug("tau_bins ={}", tau_bins);

  // We calculate the photo-ionization rate for a specific species by summing
  // the following expression over all energy bins:
  // Gamma = (L_b * exp(-tau_b) * (1 - exp(-tau_b^s))) / (E_b * vshell * n_s)
  // where the subscript "b" represents the value of the corresponding
  // quantities in the bin.
  // double localrate = 0.0;
  for (int b = species_threshold_bin_index[species]; b < bin_size; b++) {
    // spdlog::debug("bin {}, species_threshold_bin_index {}", b,
    // species_threshold_bin_index[species]);
    bin_rate = 0.0;
    bin_rate =
        bin_spectrum[b] * exp(-tau_bins[b]) / mean_energy_table[species][b];
    // calculating dtau for the specific species
    delta_species_tau = species_num_density * bin_frac_table[species][b]
                        * xsection_table[species][b] * ds;
    // spdlog::debug("delta tau_bins ={}", delta_species_tau);
    // if dtau << 1 then take the following value
    if (delta_species_tau < 1.0E-02)
      bin_rate *= delta_species_tau / (vshell * element_num_density);
    // if dtau approx 1 then take the following value
    else
      bin_rate *=
          (1.0 - exp(-delta_species_tau)) / (vshell * element_num_density);
    // adding bin rate to the total rate
    rate += bin_rate;
    // localrate += bin_spectrum[b] * exp(-tau_bins[b]) *
    // bin_frac_table[species][b] * xsection_table[species][b] * ds *
    // species_num_density / (vshell * element_num_density *
    // mean_energy_table[species][b]);
  }
  // spdlog::info("discrete rate: {:12.3e}, local rate:
  // {:12.3e}",rate,localrate);
  return rate;
}


// CALCULATING PHOTO-HEATING RATE #####################################
double photoionisation::calculate_photoheating_rate(
    const int species,                   ///< Species_identifier
    const double species_num_density,    ///< Species number density,
    const double vshell,                 ///< vshell
    const double ds,                     ///< ds
    const std::vector<double> &tau_bins  ///< Optical depth in each energy bin
)
{
  double heat_rate;
  double bin_heat_rate;
  double delta_species_tau;
  double frac_excess_energy;
  double exponential;
  // We calculate the photo-heating rate for a specific species by summing
  // the following expression over all energy bins:
  // Gamma = (L_b * exp(-tau_b) * (1 - exp(-tau_b^s))) (E_b - E_th) / (E_b *
  // vshell * n_s) where the subscript "b" represents the value of the
  // corresponding quantities in the bin.
  heat_rate = 0.0;
  for (int b = species_threshold_bin_index[species]; b < bin_size; b++) {

    exponential = exp(-tau_bins[b]);

    // if (exponent  < 1.0e-50) {
    //  heat_rate += 0.0;
    //  continue; // Skip to the next iteration
    //}

    bin_heat_rate = bin_spectrum[b] * exponential;

    // the fractional excess energy
    frac_excess_energy = 1.0
                         - ionisation_thres[species] * pconst.eV()
                               / mean_energy_table[species][b];

    bin_heat_rate *= frac_excess_energy;

    // calculating dtau for the specific species
    delta_species_tau = species_num_density * bin_frac_table[species][b]
                        * xsection_table[species][b] * ds;

    // if dtau << 1 then take the following value
    if (delta_species_tau < 1.0E-02) {
      bin_heat_rate *= delta_species_tau / vshell;
    }
    // if dtau approx > 1.0E-03 then take the following value
    else {
      bin_heat_rate *= (1.0 - exp(-delta_species_tau)) / vshell;
    }
    // spdlog::debug("bin ={}, bin_heat_rate = {} ", b, bin_heat_rate);
    // summing over energy bins
    heat_rate += bin_heat_rate;
  }

  // spdlog::debug("{}");
  return heat_rate;
}
