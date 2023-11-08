//
// Created by tony on 19/06/2023.
//

#ifndef PION_MULTI_ION_PHOTO_XSECTIONS_H
#define PION_MULTI_ION_PHOTO_XSECTIONS_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "atomic_physics_database.h"
#include <map>


/**
 * \brief Struct to hold fit parameters of photo-ionisation cross-sections.
 **/ ///< photo-ionisation fit parameter X (dimensionless)
struct photo_xsection_datatype {
  double fit_emax;  ///< Fit parameters should not be extrapolated above this
                    ///< value (in the units of eV)
  double E0;        ///< photo-ionisation fit parameter E0 (in the units of eV)
  double sigma0;    ///< photo-ionisation fit parameter sigma0 (in the units of
                    ///< 10^{-8} cm^2)
  double ya;        ///< photo-ionisation fit parameter ya (dimensionless)
  double P;         ///< photo-ionisation fit parameter P (dimensionless)
  double yw;        ///< photo-ionisation fit parameter yw (dimensionless)
  double y0;        ///< photo-ionisation fit parameter y0 (dimensionless)
  double y1;        ///< photo-ionisation fit parameter y1 (dimensionless)
};

// MULTI ION PHOTOIONISATION CROSS SECTION CLASS ############################
/**
 * \brief Class to generate
 * \details
 **/
class multi_ion_photo_xsections : virtual public atomic_physics_database {
  // PUBLIC METHODS AND VARIABLE **** *******************
public:
  /** \brief Constructor: Setup photo-ionisation cross-section database. **/
  multi_ion_photo_xsections();

  /** \brief Destructor **/
  ~multi_ion_photo_xsections();

  /** \brief ????? **/
  std::vector<double> setup_photo_xsections(
      const std::vector<std::string> &  ///< photo tracer species
  );

  /** \brief ????? **/
  void generate_ebin_mean_xsection(
      const std::vector<std::vector<double> > &  ///< array of energy bins
  );

  /** \brief Record of mean photo-ionisation cross-section of
   * photo tracer list.
   * **/
  std::vector<std::vector<double> > photo_xsection_table;

  /** \brief Record proportion of the energy bin that results in a non-zero
   * cross-section for the photo tracer list.
   * **/
  std::vector<std::vector<double> > bin_weight_table;


private:
  ///
  double species_xsection(
      const std::string,  ///< Name of the species
      const double &      ///< Energy of the photon in the units of eV
  );

  ///
  std::vector<double> mean_bin_xsection(
      const std::string,  ///< Name of the species
      const double &,
      const std::vector<double> &  ///< Range of particular energy bin
  );


  /// Map coll_ionise_datatype to the species string (species symbol)
  std::map<std::string, photo_xsection_datatype> photo_xsection_parameters;

  ///
  void multi_ion_photo_xsections_database();


  std::vector<std::string> photo_species_list;  ///<

  std::vector<double> photo_threshold;  ///< Ionisation threshold energy vector
};



#endif  // PION_MULTI_ION_PHOTO_XSECTIONS_H
