/**
 *  \file     mellema_cooling.h
 *  \author   Arun Mathew
 *  \author   Jonathan Mackey,
 *  \date     18-01-2022
 **/

#ifndef MELLEMA_COOLING_H
#define MELLEMA_COOLING_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "tools/interpolate.h"


// USER DEFINED DATATYPE #################################################

/**
 * \brief Struct to hold original Mellema table attributes.
 **/
struct database_datatype {
  std::string Name;           ///< Species name
  std::vector<double> logT;   ///< Original temperature table (log10)
  std::vector<double> logne;  ///< Electron number density table (log10)
  std::vector<std::vector<double> >
      rate;  ///< Original Mellema cooling rate table
};

/**
 * \brief Struct to hold Mellema cooling LUT attributes.
 **/
struct LUT_datatype {
  /// Name of the species
  std::string Name;
  /// New Mellema cooling rate lookup table
  std::vector<std::vector<double> > rate;
  /// Mellema rate slope along T-axis
  std::vector<std::vector<double> > T_rateslope;
  /// Mellema rate slope along ne-axis
  std::vector<std::vector<double> > ne_rateslope;
};
// END OF USER DEFINED DATATYPE #########################################



// MELLEMA COOLING CLASS ####################################################
/** \brief Class to generate Mellema cooling LUTs
 * \details This class stack the original tabulated mellema cooling table in
 * the database and use bilinear interpolation method to generate the lookup
 * tables of arbitrary size for the mc_tracer_list. To perform bilinear
 * interpolation this class inherent root_find_bilinear_vec() function from
 * the class interpolate_array.
 * */
class mellema_cooling : virtual public interpolate_arrays {

  // PUBLIC METHODS AND VARIABLE ***********************
public:
  /**
   * \brief Constructor: Initializing Mellema colling database
   **/
  mellema_cooling();

  /**
   * \brief Destructor: Erase database, once the LUTs are generated.
   **/
  ~mellema_cooling();

  /**
   * \brief This function will set up the species locator for the given tracer
   * list. By looping through its database against the tracer list, this
   * function will record the location (integer value) of the species in the
   * tracer list in the database into 2D vector mc_tracer_locator.
   **/
  void setup_mellema_cooling(const std::vector<std::vector<std::string> >
                                 &  ///< MP tracer list for cooling.
  );

  /**
   * \brief This function generates Mellema cooling LUTs for a given temperature
   * and electron number density tables using bilinear interpolation for the
   * species in the mc_tracer_locator. Since Mellema cooling rate coefficients
   * are both temperature and electron number density functions, this function
   * will also generate the first partial derivative of the rate coefficients
   * along temperature and electron number density.
   **/
  void generate_mellema_table(
      const std::vector<double> &,  ///< Temperature table
      const std::vector<double> &   ///< Electron number density table
  );

  std::vector<std::vector<int> >
      mc_tracer_locator;  ///< Location of tracer species in the database.
  std::vector<LUT_datatype> mellema_table;  ///< LUT datatype variable

  // PRIVATE METHODS AND VARIABLES *************************
private:
  /** \brief  This function stash the original cooling rates for the species H,
   * H1+, He, He1+, He2+, C, C1+, C2+, C3+, C4+, C5+, C6+, N, N1+, N2+, N3+,
   * N4+, N5+, N6+, N7+, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+, Ne, Ne1+,
   * Ne2+, Ne3+, Ne4+, Ne5+, Ne6+, Ne7+, Ne8+, Ne9+, Ne10+, S, S1+, S2+, S3+,
   * S4+.
   *
   * These rate tables are taken from Raga (1997)
   * <a href="https://iopscience.iop.org/article/10.1086/312987/meta">
   *[ref]</a>.
   *
   * Cooling table of Hydrogen include collisional excitation of Lyman alpha,
   * radiative recombination and free-free emission. For other species
   * collisional excitation of optical and ultraviolet lines are included.
   **/
  void cooling_table_database();

  std::vector<database_datatype>
      species;  ///< Cooling table database data type.
};
// END OF MELLEMA COOLING CLASS #############################################
#endif  // MELLEMA_COOLING_H
