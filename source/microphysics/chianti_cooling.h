/**
 *  \file     chianti_cooling.cpp
 *  \author   Arun Mathew
 *  \author   Jonathan Mackey,
 *  \date     10-02-2023
 **/

#ifndef CHIANTI_COOLING_H
#define CHIANTI_COOLING_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "tools/interpolate.h"


// USER DEFINED DATATYPE #################################################

/**
 * \brief Struct to hold chianti cooling table attributes.
 **/
struct database_datatype {
  std::string Name;                        ///< Species name
  std::vector<std::vector<double> > rate;  ///< chianti cooling rate table
};

/**
 * \brief Struct to hold chianti cooling LUT attributes.
 **/
struct LUT_datatype {
  /// Name of the species
  std::string Name;
  /// New chianti cooling rate lookup table
  std::vector<std::vector<double> > rate;
  /// Chianti rate slope along T-axis
  std::vector<std::vector<double> > T_rateslope;
  /// Chianti rate slope along ne-axis
  std::vector<std::vector<double> > ne_rateslope;
};
// END OF USER DEFINED DATATYPE #########################################



// CHIANTI COOLING CLASS ####################################################
/** \brief Class to generate Chianti cooling LUTs
 * \details This class stack the tabulated chianti cooling table in
 * the database and use bilinear interpolation method to generate the lookup
 * tables of arbitrary size for the mc_tracer_list. To perform bilinear
 * interpolation this class inherent root_find_bilinear_vec() function from
 * the class interpolate_array.
 * */
class chianti_cooling : virtual public interpolate_arrays {

  // PUBLIC METHODS AND VARIABLE ***********************
public:
  /**
   * \brief Constructor: Initializing Chianti colling database
   **/
  chianti_cooling();

  /**
   * \brief Destructor: Erase database, once the LUTs are generated.
   **/
  ~chianti_cooling();

  /**
   * \brief This function will set up the species locator for the given tracer
   * list. By looping through its database against the tracer list, this
   * function will record the location (integer value) of the species in the
   * tracer list in the database into 2D vector chianti_tracer_locator.
   **/
  void setup_chianti_cooling(const std::vector<std::vector<std::string> >
                                 &  ///< MP tracer list for cooling.
  );

  /**
   * \brief This function generates Chianti cooling LUTs for a given temperature
   * and electron number density tables using bilinear interpolation for the
   * species in the chianti_tracer_locator. Since Chianti cooling rate
   *coefficients are both functions of temperature and electron number density,
   *this function will also return the first partial derivative of the rate
   *coefficients with respect to temperature and electron number density for the
   *same tables of temperature and electron density.
   **/
  void generate_chianti_table(
      const std::vector<double> &,  ///< Temperature table
      const std::vector<double> &   ///< Electron number density table
  );

  std::vector<double> chianti_logT_table;
  std::vector<double> chianti_logne_table;
  std::vector<std::vector<int> >
      chianti_tracer_locator;  ///< Location of tracer species in the database.
  std::vector<LUT_datatype> chianti_table;  ///< LUT datatype variable

  // PRIVATE METHODS AND VARIABLES *************************
private:
  /** \brief  This function stash the original cooling rates for the species H,
   * H1+, He, He1+, He2+, C, C1+, C2+, C3+, C4+, C5+, C6+, N, N1+, N2+, N3+,
   * N4+, N5+, N6+, N7+, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+, Ne, Ne1+,
   * Ne2+, Ne3+, Ne4+, Ne5+, Ne6+, Ne7+, Ne8+, Ne9+, Ne10+, Si, Si1+, Si2+,
   * Si3+, Si4+, Si5+, Si6+, Si7+, Si8+, Si9+, Si10+, Si11+, Si12+, Si13+,
   * Si14+, S, S1+, S2+, S3+, S4+, S5+, S6+, S7+, S8+, S9+, S10+, S11+, S12+,
   * S13+, S14+, S15+, S16+, Fe, Fe1+, Fe2+, Fe3+, Fe4+, Fe5+, Fe6+, Fe7+,
   * Fe8+, Fe9+, Fe10+, Fe11+, Fe12+, Fe13+, Fe14+, Fe15+, Fe16+, Fe17+, Fe18+,
   * Fe19+, Fe20+, Fe21+, Fe22+, Fe23+, Fe24+, Fe25+, Fe26+.
   *
   * The rate tables were calculated using CHIANTIPy version 10.0 [1, 2],
   * without make assumptions about ionization equilibrium or abundances.
   * These values are valid  for 10K to 10^9K, and electron densities ranging
   * from 1cm^-3 to 10^6cm^-3.
   *
   * Cooling rates include bremsstrahlung, radiative recombination, line
   * radiation and two-photon emission. The radiative loss include
   * line-emissions of all spectral lines of the specific ion [2, 3]. The
   * free-bound emissions following the prescription of [4]. Bremsstrahlung
   * emission using the expression given in [5] with the Gaunt factor calculated
   * using the method of [6]. The two-photon emission of hydrogen [7] and
   * helium [8] isoelectronic sequence ions with the spectral distribution
   * function taken from [9] for hydrogen isoelectronic ions and [10] for
   * helium isoelectronic ions.
   *
   * Referance
   * [1] Zanna G. D., Dere K. P., Young P. R., Landi E., 2021, ApJ, 909, 38
   * [2] Young P. R. et al., 2003, ApJS, 144, 135
   * [3] Zanna G. D., Dere K. P., Young P. R., Landi E., Mason H. E., 2015, A&A,
   *582, A56 [4] Mewe R., Lemen J. R., van den Oord G. H. J., 1986, A&AS, 65,
   *511 [5] Rybicki G. B., Lightman A. P., 1986, Radiative Processes in
   *Astrophysics [6] Karzas W. J., Latter R., 1961, ApJS, 6, 167 [7] Parpia F.
   *A., Johnson W. R., 1982, Phys. Rev. A, 26, 1142 [8] Drake G. W. F., 1986,
   *Phys. Rev. A, 34, 2871 [9] Goldman S. P., Drake G. W. F., 1981, Phys. Rev.
   *A, 24, 183 [10] Drake G. W. F., Victor G. A., Dalgarno A., 1969, Phys. Rev.,
   *180, 25
   *
   **/
  void cooling_table_database();

  std::vector<database_datatype>
      species;  ///< Cooling table database data type.
};
// END OF CHIANTI COOLING CLASS #############################################
#endif  // CHIANTI_COOLING_H
