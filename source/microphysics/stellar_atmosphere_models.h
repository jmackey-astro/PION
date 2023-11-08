/**
 *  \file     stellar_atmosphere_models.h
 *  \author   Arun Mathew
 *  \date     05-07-2023
 **/

#ifndef PION_STELLAR_ATMOSPHERE_MODELS_H
#define PION_STELLAR_ATMOSPHERE_MODELS_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include "microphysics_base.h"
#include <fmt/ranges.h>
#include <fstream>
#include <sys/stat.h>
#include <variant>


struct stellar_model_datatype {
  string info;
  std::vector<std::vector<double> > data;
};

struct star_datatype {
  int Type;                 // identifier
  double Teff;              // Effective temperature of the star
  double Radius;            // Stellar radius
  double Luminosity;        // Total luminosity
  double Stefan_Boltzmann;  // Luminosity
};


class stellar_atmosphere_models {

public:
  stellar_atmosphere_models();
  ~stellar_atmosphere_models();

  void setup_stellar_atmosphere_models(
      std::vector<std::vector<double> >
          &,    ///< Energy bins (in the units of eV)
      int &,    ///< Number of energy bins
      double &  ///< Maximum photo-ionisation energy (in the units of eV)
  );

  /** \brief This function calculate the luminosity number in different
   * energy bin for the given radiation source.
   **/
  int generate_bin_luminosity(
      const int,             ///< Radiation source id
      const double,          ///< Effective star temperature
      const double,          ///< Stellar radius
      const double,          /// Total luminosity of the star
      std::vector<double> &  ///< O/P luminosity per energy bin (number/s/bin).
  );

private:
  std::vector<std::vector<double> > energy_bins;
  int N_energy_bins;

  void stellar_atmosphere_model_database();

  // get stellar atmosphere model
  stellar_model_datatype get_stellar_atmosphere_model(int);

  star_datatype star_info;

  // stellar atmosphere models;
  stellar_model_datatype atlas_ckp00_g45;  ///< O3V STAR
  stellar_model_datatype powr_mw_wc;       ///< WC
  stellar_model_datatype powr_mw_wne;      ///< WNE
  // add models here


  /// Debug tools: write bin luminosity (monochromatic luminosity) to a file
  void print_bin_luminosity(
      const std::string,                          // model info
      const star_datatype,                        // star info
      const std::vector<double> &,                // bin luminosity
      const std::vector<std::vector<double> > &,  // energy bin
      std::string                                 // Base name/file path
  );
};



#endif  // PION_STELLAR_ATMOSPHERE_MODELS_H
