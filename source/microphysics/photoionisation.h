/**
 *  \file     photoionisation.h
 *  \author   Arun Mathew
 *  \date     03-07-2023
 **/

#ifndef PION_PHOTOIONISATION_H
#define PION_PHOTOIONISATION_H
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "constants.h"
#include "tools/interpolate.h"


class photoionisation : virtual public interpolate_arrays {
public:
  photoionisation();
  ~photoionisation();

  /** \brief This function will adjust the number of photo-ionization
   * rate species tables and the photo-ionization heating rate species
   * tables based on the photo tracer list.
   **/
  void setup_photoionisation(
      const std::vector<std::vector<double> > &,  ///< energy bins
      const int,                                  ///< Number of photo species
      const std::vector<std::vector<double> > &,  ///< mean xsection table
      const std::vector<std::vector<double> > &,  ///< bin fraction table
      const std::vector<double> &,                ///< stellar spectrum
      const std::vector<double> &                 ///< photoionisation threshold
  );

  double calculate_photoionisation_rate(
      const int,                   ///< Species_identifier
      const double,                ///< Species number density,
      const double,                ///< Element number density,
      const double,                ///< vshell
      const double,                ///< ds
      const std::vector<double> &  ///< Optical depth in each energy bin
  );

  double calculate_photoheating_rate(
      const int,                   ///< Species_identifier
      const double,                ///< Species number density,
      const double,                ///< vshell
      const double,                ///< ds
      const std::vector<double> &  ///< Optical depth in each energy bin
  );


private:
  int bin_size;
  std::vector<std::vector<double> > E_bins;
  std::vector<std::vector<double> > mean_energy_table;

  /// This table record the index of the bin where the threshold energy
  /// lies for each photo-species.
  std::vector<int> species_threshold_bin_index;
  std::vector<std::vector<double> > xsection_table;
  std::vector<std::vector<double> > bin_frac_table;
  std::vector<double> bin_spectrum;
  std::vector<double> ionisation_thres;
};
#endif  // PION_PHOTOIONISATION_H
