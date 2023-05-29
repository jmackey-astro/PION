/**
 *  \file     atomic_physics_database.h
 *  \author   Arun Mathew
 *  \date     12/12/2022
 **/

#ifndef PION_ATOMIC_PHYSICS_DATABASE_H
#define PION_ATOMIC_PHYSICS_DATABASE_H

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include "constants.h"
#include <map>
#include <spdlog/spdlog.h>


/** \brief
 * Struct to hold atomic physics species attributes.
 **/
struct species_datatype {
  std::string element;  ///< Element to which the ion belong to.
  int charge;           ///< Charge of the species.
  double mass;          ///< Mass of the element.
  double ionise_E;  ///< Ionisation energy of species (in the units of ergs).
};


// ATOMIC PHYSICS DATABASE CLASS #########################################
/**
 * \brief Class for atomic physics data repository
 * \details This class serves as the repository for atomic physics attributes
 * of neutral atoms and their ions. The attributes include element name, charge,
 * mass, ionisation potential energy.
 **/
class atomic_physics_database {

  // PUBLIC METHODS AND VARIABLE ***********************
public:
  /** \brief Constructor **/
  atomic_physics_database();

  /** \brief Destructor **/
  ~atomic_physics_database();

  /** \brief Returns the element name string of the species. **/
  std::string get_element_name(const std::string species);

  /** \brief Returns the atomic mass of the element. **/
  double get_atomic_mass(const std::string  ///< element name
  );

  /** \brief Returns the charge of the species. **/
  int get_charge(const std::string  ///< species name
  );

  /** \brief Returns the number of electron released by the species. **/
  int get_electron_num(const std::string  ///< species name
  );

  /** \brief Returns the ionisation potential of the species in the units of
   * ergs.
   **/
  double get_ionisation_potential(const std::string  ///< species name
  );

  // PRIVATE METHODS AND VARIABLES *************************
private:
  /** \brief  Atomic physics data of the following species: H, H1+, He, He1+,
   * He2+, C, C1+, C2+, C3+, C4+, C5+, C6+, N, N1+, N2+, N3+, N4+, N5+, N6+,
   * N7+, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+, Ne, Ne1+, Ne2+, Ne3+, Ne4+,
   * Ne5+, Ne6+, Ne7+, Ne8+, Ne9+, Ne10+, Si, Si1+, Si2+, Si3+, Si4+, Si5+,
   * Si6+, Si7+, Si8+, Si9+, Si10+, Si11+, Si12+, Si13+, Si14+, S, S1+, S2+,
   * S3+, S4+, S5+, S6+, S7+, S8+, S9+, S10+, S11+, S12+, S13+, S14+, S15+,
   * S16+, Fe, Fe1+, Fe2+, Fe3+, Fe4+, Fe5+, Fe6+, Fe7+, Fe8+, Fe9+, Fe10+,
   * Fe11+, Fe12+, Fe13+, Fe14+, Fe15+, Fe16+, Fe17+, Fe18+, Fe19+, Fe20+,
   * Fe21+, Fe22+, Fe23+, Fe24+, Fe25+, Fe26+.
   *
   * Atomic masses are taken from Juris Meija et al. (2015), "Atomic weights
   * of the elements 2013 (IUPAC Technical Report)" Pure and Applied
   * Chemistry, vol. 88, no. 3, 2016, pp. 265-291.
   * <a href="https://doi.org/10.1515/pac-2015-0305"> [ref]</a>
   *
   * Values of ionisation potential energies are taken from Voronov (1997)
   * <a href="https://doi.org/10.1006/adnd.1997.0732"> [ref]</a>
   *
   **/
  void database();

  /// Map species_datatype to the species (species symbol)
  std::map<std::string, species_datatype> species_attributes;
};
// END OF ATOMIC PHYSICS DATABASE CLASS ####################################

#endif  // PION_ATOMIC_PHYSICS_DATABASE_H