/**
 *  \file     coll_ionise_recomb.h
 *  \author   Arun Mathew
 *  \date     19-12-2022
 **/

#ifndef PION_COLL_IONISE_RECOMB_H
#define PION_COLL_IONISE_RECOMB_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "atomic_physics_database.h"
#include "constants.h"
#include <map>

#define HYDROGEN_CASE_B

/// Switch to use case B radiative recombination rate for H1+
//#define HYDROGEN_CASE_B

/**
 * \brief Struct to hold fit parameters of collisional ionisation rate
 * coefficient.
 **/
struct coll_ionise_datatype {
  int P;     ///< Collisional ionisation fit parameter P (dimensionless)
  double A;  ///< Collisional ionisation fit parameter A (in units of cm^3/s)
  double X;  ///< Collisional ionisation fit parameter X (dimensionless)
  double K;  ///< Collisional ionisation fit parameter K (dimensionless)
};

/**
 * \brief Struct to hold fit parameters of Radiative and Di-electronic
 * recombination rate coefficient.
 **/
struct recomb_datatype {
  double A;   ///< Radiative recombination fit parameter A  (in units of cm^3/s)
  double B;   ///< Radiative recombination fit parameter B  (dimensionless)
  double T0;  ///< Radiative recombination fit parameter T0 (in units of K)
  double T1;  ///< Radiative recombination fit parameter T1 (in units of K)
  double C;   ///< Radiative recombination fit parameter C  (dimensionless)
  double T2;  ///< Radiative recombination fit parameter T2 (in units of K)
  std::vector<double> C_coefs{std::vector<double>(9)};
  ///< Dielectronic recombination C_i coefficients (in units of cm^3K^1.5/s)
  std::vector<double> E_coefs{std::vector<double>(9)};
  ///< Dielectronic recombination E_i coefficients (in units of K)
};


// COLLISIONAL IONISATION AND RECOMBINATION CLASS ############################
/**
 * \brief Class to generate Ionisation and Recombination LUTs
 * \details This class serves as the repository for the fit parameters
 * of collisional ionisation, radiative recombination, and dielectronic
 * recombination fitting formulae. This class can be inherited
 * and set up using setup functions with the preferred species list to
 * generate the corresponding look-up tables.
 **/
class coll_ionise_recomb : virtual public atomic_physics_database {
  // PUBLIC METHODS AND VARIABLE ***********************
public:
  /** \brief Constructor: Setup both collisional and recombination database. **/
  coll_ionise_recomb();

  /** \brief Destructor **/
  ~coll_ionise_recomb();

  /** \brief This function will set the species size of the collisional
   *ionisation rate lookup tables and slope lookup tables to the number of
   *species in collisional ionisation tracer list. The function also set up the
   * collisional ionisation tracer's ionisation potential vector for the
   * species in collisional ionisation tracer list.
   **/
  void setup_collisional_ionisation(
      const std::vector<std::string> &  ///< collisional ionisation tracer list
  );

  /** \brief This function will set the species size of the recombination rate
   * lookup tables and slope lookup tables.
   **/
  void setup_recombination(
      const std::vector<std::string> &  ///< recombination tracer list
  );

  /** \brief Function generating collisional ionisation rate and slope lookup
   * tables for a given array of temperature values with ci_tracer species list.
   **/
  void generate_collisional_ionisation_LUT(
      const std::vector<double> &  ///< Temperature table
  );

  /** \brief Function generating recombination rate and slope lookup tables
   * for a given array of temperature values with recomb_tracer species list.
   **/
  void generate_recombination_LUT(
      const std::vector<double> &  ///< Temperature table
  );

  std::vector<std::vector<double> >
      collisional_rate_table;  ///< Collisional ionisation rate LUT
  std::vector<std::vector<double> >
      collisional_slope_table;  ///< Collisional ionisation rate slope LUT
  std::vector<std::vector<double> >
      recombination_rate_table;  ///< Recombination rate LUT
  std::vector<std::vector<double> >
      recombination_slope_table;  ///< Recombination rate slope LUT
  std::vector<double>
      ci_ion_pot;  ///< Ionisation potential vector for ci_tracer_list


  // PRIVATE METHODS AND VARIABLES *************************
private:
  /** \brief Collisional ionisation fit parameter database for H, He, He1+, C,
   * C1+, C2+, C3+, C4+, C5+, N, N1+, N2+, N3+, N4+, N5+, N6+, O, O1+, O2+,
   * O3+, O4+, O5+, O6+, O7+, Ne, Ne1+, Ne2+, Ne3+, Ne4+, Ne5+, Ne6+, Ne7+,
   * Ne8+, Ne9+, Si, Si1+, Si2+, Si3+, Si4+, Si5+, Si6+, Si7+, Si8+, Si9+,
   * Si10+, Si11+, Si12+, Si13+, S, S1+, S2+, S3+, S4+, S5+, S6+, S7+, S8+, S9+,
   * S10+, S11+, S12+, S13+, S14+, S15+, Fe, Fe1+, Fe2+, Fe3+, Fe4+, Fe5+, Fe6+,
   * Fe7+, Fe8+, Fe9+, Fe10+, Fe11+, Fe12+, Fe13+, Fe14+, Fe15+, Fe16+, Fe17+,
   * Fe18+, Fe19+, Fe20+, Fe21+, Fe22+, Fe23+, Fe24+, Fe25+.
   *
   * Collisional ionisation rate coefficient is described by the fit formula
   * (with parameters P, A, X, and K) given by Voronov (1997)
   * <a href="https://doi.org/10.1006/adnd.1997.0732"> [ref]</a>.
   */
  void collisional_ionisation_database();

  std::vector<std::string>
      ci_species_list;  ///< Collisional ionisation tracer list
  std::vector<std::string> recomb_species_list;  ///< Recombination tracer list

  /** \brief Radiative and Dielectronic recombination fit parameter database
   * for H1+, He1+, He2+, C1+, C2+, C3+, C4+, C5+, C6+, N1+, N2+, N3+, N4+,
   * N5+, N6+, N7+, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+, Ne1+, Ne2+, Ne3+,
   * Ne4+, Ne5+, Ne6+, Ne7+, Ne8+, Ne9+, Ne10+, Si1+, Si2+, Si3+, Si4+, Si5+,
   * Si6+, Si7+, Si8+, Si9+, Si10+, Si11+, Si12+, Si13+, Si14+, S2+, S3+, S4+,
   * S5+, S6+, S7+, S8+, S9+, S10+, S11+, S12+, S13+, S14+, S15+, S16+, Fe8+,
   * Fe12+, Fe13+, Fe14+, Fe15+, Fe16+, Fe17+, Fe18+, Fe19+, Fe20+, Fe21+,
   * Fe22+, Fe23+, Fe24+, Fe25+, Fe26+.
   *
   * The parameters A, B, T0, T1, C, T2 for radiative recombination rate
   * coefficient are taken from Badnell (2006)
   * <a href="https://iopscience.iop.org/article/10.1086/508465"> [ref]</a>.
   * Dielectronic recombination rate coefficients with parameters C_coefs
   * and E_coefs are taken from Zatsarinny (2003)
   * <a href="https://www.aanda.org/articles/aa/pdf/2003/48/aa4001.pdf">
   * [ref]</a>. All species dielectronic recombination rate coefficients rates
   * can be found in <a href="http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/">
   * [ref]</a>.
   *
   *
   * We use the old values for Fe11+, Fe10+, Fe9+, Fe7+, Fe6+, Fe5+, Fe4+,
   * Fe3+, Fe2+, Fe1+. This can be found at
   * <a http://www.pa.uky.edu/~verner/fortran.html> [ref]</a>
   *
   * Verner-Ferland fit parameters for the radiative recombination of S1+
   * species is missing in the literature. However, the Aldrovandi-Pequignot
   * fit parameter for the same is available at
   * <a https://www.pa.uky.edu/~verner/dima/rec//hdr.txt.>[ref]</a>.
   * We use a Python code to obtain appropriate Verner-Ferland fit parameters
   * for S1+ from the already known Aldrovandi-Pequignot fit parameters.
   * Input Aldrovandi-Pequignot fit parameter: A_rad = 4.1E-13, eta = 0.630
   * Obtained Verner-Ferland fit parameters:
   * A = 7.50E-12, B = 6.45E-01, T_0 = 1.75E+02, T_0 = 3.70e+08
   *
   * The dielectronic recombination is also taken an old data set at
   * <a https://www.pa.uky.edu/~verner/dima/rec//hdr.dat > [ref]</a>
   * C_1 = 1.62E-03, E 1.25E+05
   *
   * Relevent References:
   * S. M. V. Aldrovandi & D. Pequignot (1973, A&A, 25, 137),
   * J. M. Shull & M. Van Steenberg (1982, ApJS, 48, 95),
   * M. Arnaud & R. Rothenflug (1985, A&AS, 60, 425).
   *
   */
  void recombination_database();

  /// Map coll_ionise_datatype to the species string (species symbol)
  std::map<std::string, coll_ionise_datatype> coll_ionise_parameters;
  /// Map recomb_datatype to the species string (species symbol)
  std::map<std::string, recomb_datatype> recomb_parameters;

  /** \brief Calculate collisional ionisation rate of a species for a given
   * temperature T.
   **/
  double get_collisional_ionisation_rate(
      const std::string,  ///< species name
      const double        ///< Temperature
  );

  /** \brief Calculate recombination rate of a species for a given temperature
   * T.
   **/
  double get_recombination_rate(
      const std::string,  ///< species name
      const double        ///< Temperature
  );
};
// END OF COLLISIONAL IONISATION AND RECOMBINATION CLASS #####################


#endif  // PION_COLL_IONISE_RECOMB_H
