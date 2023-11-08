/**
 *  \file     atomic_physics_database.h
 *  \author   Arun Mathew
 *  \date     12/12/2022
 **/

#include "atomic_physics_database.h"

// Global instance of class
class atomic_physics_database atomic_physics_data;


// ##################################################################
// Constructor
atomic_physics_database::atomic_physics_database()
{
  database();
}


// ##################################################################
// Destructor
atomic_physics_database::~atomic_physics_database() {}



// ##################################################################
// Get atomic element from the database
std::string atomic_physics_database::get_element_name(const std::string species)
{
  return species_attributes[species].element;
}



// ##################################################################
// Get atomic mass from the database
double atomic_physics_database::get_atomic_mass(const std::string species)
{
  return species_attributes[species].mass;
}



// ##################################################################
// Get the charge of the species
int atomic_physics_database::get_charge(const std::string species)
{
  return species_attributes[species].charge;
}



// ##################################################################
// Get number of electron released by the species
int atomic_physics_database::get_electron_num(const std::string species)
{
  return species_attributes[species].charge;
}



// ##################################################################
// Get ionisation potential in units of erg
double atomic_physics_database::get_ionisation_potential(
    const std::string species)
{
  return species_attributes[species].ionise_E * pconst.eV();
}



// ##################################################################
// Get ionisation potential in units of eV
double atomic_physics_database::get_ionisation_energy(const std::string species)
{
  return species_attributes[species].ionise_E;
}



// ##################################################################
// Atomic physics database
void atomic_physics_database::database()
{
  spdlog::info("Setting up atomic physics database");

  // Ionisation energies are taken from NIST Atomic Spectra
  // Database (ver.5.10), Kramida, A., Ralchenko, Yu., Reader,
  // J., and NIST ASD Team (2022).
  // Available at DOI: https://doi.org/10.18434/T4W30F
  // Ionisation energies are given in the units of eV.

  // Entry: H
  species_attributes["H"] = species_datatype{"H", 0, 1.6738e-24, 13.59844};

  // Entry: H1+
  species_attributes["H1+"] = species_datatype{
      "H", 1,
      0,       // mass of H1+ is set equal to zero.
      1.0e-99  // small dummy value.
  };

  // Entry: He
  species_attributes["He"] = species_datatype{"He", 0, 6.6464768e-24, 24.58741};

  // Entry: He1+
  species_attributes["He1+"] = species_datatype{"He", 1, 0, 54.41778};

  // Entry: He2+
  species_attributes["He2+"] = species_datatype{"He", 2, 0, 1.0e-99};

  // Entry: C
  species_attributes["C"] = species_datatype{"C", 0, 1.994374e-23, 11.2602880};

  // Entry: C1+
  species_attributes["C1+"] = species_datatype{"C", 1, 0, 24.383143};

  // Entry: C2+
  species_attributes["C2+"] = species_datatype{"C", 2, 0, 47.88778};

  // Entry: C3+
  species_attributes["C3+"] = species_datatype{"C", 3, 0, 64.49352};

  // Entry: C4+
  species_attributes["C4+"] = species_datatype{"C", 4, 0, 392.090518};

  // Entry: C5+
  species_attributes["C5+"] = species_datatype{"C", 5, 0, 489.99320779};

  // Entry: C6+
  species_attributes["C6+"] = species_datatype{"C", 6, 0, 1.0e-99};

  // Entry: N
  species_attributes["N"] = species_datatype{"N", 0, 2.325892e-23, 14.53413};

  // Entry: N1+
  species_attributes["N1+"] = species_datatype{"N", 1, 0, 29.60125};

  // Entry: N2+
  species_attributes["N2+"] = species_datatype{"N", 2, 0, 47.4453};

  // Entry: N3+
  species_attributes["N3+"] = species_datatype{"N", 3, 0, 77.4735};

  // Entry: N4+
  species_attributes["N4+"] = species_datatype{"N", 4, 0, 97.8901};

  // Entry: N5+
  species_attributes["N5+"] = species_datatype{"N", 5, 0, 552.06733};

  // Entry: N6+
  species_attributes["N6+"] = species_datatype{"N", 6, 0, 667.0461377};

  // Entry: N7+
  species_attributes["N7+"] = species_datatype{"N", 7, 0, 1.0e-99};

  // Entry: O
  species_attributes["O"] = species_datatype{"O", 0, 2.6567628e-23, 13.618055};

  // Entry: O1+
  species_attributes["O1+"] = species_datatype{"O", 1, 0, 35.12112};

  // Entry: O2+
  species_attributes["O2+"] = species_datatype{"O", 2, 0, 54.93554};

  // Entry: O3+
  species_attributes["O3+"] = species_datatype{"O", 3, 0, 77.41350};

  // Entry: O4+
  species_attributes["O4+"] = species_datatype{"O", 4, 0, 113.8990};

  // Entry: O5+
  species_attributes["O5+"] = species_datatype{"O", 5, 0, 138.1189};

  // Entry: O6+
  species_attributes["O6+"] = species_datatype{"O", 6, 0, 739.32683};

  // Entry: O7+
  species_attributes["O7+"] = species_datatype{"O", 7, 0, 871.4099138};

  // Entry: O8+
  species_attributes["O8+"] = species_datatype{"O", 8, 0, 1.e-99};

  // Entry: Ne
  species_attributes["Ne"] =
      species_datatype{"Ne", 0, 3.3509177e-23, 21.564541};

  // Entry: Ne1+
  species_attributes["Ne1+"] = species_datatype{"Ne", 1, 0, 40.96297};

  // Entry: Ne2+
  species_attributes["Ne2+"] = species_datatype{"Ne", 2, 0, 63.4233};

  // Entry: Ne3+
  species_attributes["Ne3+"] = species_datatype{"Ne", 3, 0, 97.1900};

  // Entry: Ne4+
  species_attributes["Ne4+"] = species_datatype{"Ne", 4, 0, 126.247};

  // Entry: Ne5+
  species_attributes["Ne5+"] = species_datatype{"Ne", 5, 0, 157.934};

  // Entry: Ne6+
  species_attributes["Ne6+"] = species_datatype{"Ne", 6, 0, 207.271};

  // Entry: Ne7+
  species_attributes["Ne7+"] = species_datatype{"Ne", 7, 0, 239.0970};

  // Entry: Ne8+
  species_attributes["Ne8+"] = species_datatype{"Ne", 8, 0, 1195.80784};

  // Entry: Ne9+
  species_attributes["Ne9+"] = species_datatype{"Ne", 9, 0, 1362.199256};

  // Entry: Ne10+
  species_attributes["Ne10+"] = species_datatype{"Ne", 10, 0, 1.0e-99};

  // Entry: Si
  species_attributes["Si"] = species_datatype{"Si", 0, 4.6637066e-23, 8.15168};

  // Entry: Si1+
  species_attributes["Si1+"] = species_datatype{"Si", 1, 0, 16.34585};

  // Entry: Si2+
  species_attributes["Si2+"] = species_datatype{"Si", 2, 0, 33.49300};

  // Entry: Si3+
  species_attributes["Si3+"] = species_datatype{"Si", 3, 0, 45.14179};

  // Entry: Si4+
  species_attributes["Si4+"] = species_datatype{"Si", 4, 0, 166.767};

  // Entry: Si5+
  species_attributes["Si5+"] = species_datatype{"Si", 5, 0, 205.279};

  // Entry: Si6+
  species_attributes["Si6+"] = species_datatype{"Si", 6, 0, 246.57};

  // Entry: Si7+
  species_attributes["Si7+"] = species_datatype{"Si", 7, 0, 303.59};

  // Entry: Si8+
  species_attributes["Si8+"] = species_datatype{"Si", 8, 0, 351.28};

  // Entry: Si9+
  species_attributes["Si9+"] = species_datatype{"Si", 9, 0, 401.38};

  // Entry: Si10+
  species_attributes["Si10+"] = species_datatype{"Si", 10, 0, 476.273};

  // Entry: Si11+
  species_attributes["Si11+"] = species_datatype{"Si", 11, 0, 523.415};

  // Entry: Si12+
  species_attributes["Si12+"] = species_datatype{"Si", 12, 0, 2437.65815};

  // Entry: Si13+
  species_attributes["Si13+"] = species_datatype{"Si", 13, 0, 2673.177958};

  // Entry: Si14+
  species_attributes["Si14+"] = species_datatype{"Si", 14, 0, 1.0e-99};

  // Entry: S
  species_attributes["S"] = species_datatype{"S", 0, 5.3245181e-23, 10.36001};

  // Entry: S1+
  species_attributes["S1+"] = species_datatype{"S", 1, 0, 23.33788};

  // Entry: S2+
  species_attributes["S2+"] = species_datatype{"S", 2, 0, 34.86};

  // Entry: S3+
  species_attributes["S3+"] = species_datatype{"S", 3, 0, 47.222};

  // Entry: S4+
  species_attributes["S4+"] = species_datatype{"S", 4, 0, 72.5945};

  // Entry: S5+
  species_attributes["S5+"] = species_datatype{"S", 5, 0, 88.0529};

  // Entry: S6+
  species_attributes["S6+"] = species_datatype{"S", 6, 0, 280.954};

  // Entry: S7+
  species_attributes["S7+"] = species_datatype{"S", 7, 0, 328.794};

  // Entry: S8+
  species_attributes["S8+"] = species_datatype{"S", 8, 0, 379.84};

  // Entry: S9+
  species_attributes["S9+"] = species_datatype{"S", 9, 0, 447.7};

  // Entry: S10+
  species_attributes["S10+"] = species_datatype{"S", 10, 0, 504.55};

  // Entry: S11+
  species_attributes["S11+"] = species_datatype{"S", 11, 0, 564.41};

  // Entry: S12+
  species_attributes["S12+"] = species_datatype{"S", 12, 0, 651.96};

  // Entry: S13+
  species_attributes["S13+"] = species_datatype{"S", 13, 0, 706.994};

  // Entry: S14+
  species_attributes["S14+"] = species_datatype{"S", 14, 0, 3223.7807};

  // Entry: S15+
  species_attributes["S15+"] = species_datatype{"S", 15, 0, 3494.188518};

  // Entry: S16+
  species_attributes["S16+"] = species_datatype{"S", 16, 0, 1.0e-99};

  // Entry: Fe
  species_attributes["Fe"] =
      species_datatype{"Fe", 0, 9.2732796e-23, 7.9024681};

  // Entry: Fe1+
  species_attributes["Fe1+"] = species_datatype{"Fe", 1, 0, 16.19921};

  // Entry: Fe2+
  species_attributes["Fe2+"] = species_datatype{"Fe", 2, 0, 30.651};

  // Entry: Fe3+
  species_attributes["Fe3+"] = species_datatype{"Fe", 3, 0, 54.91};

  // Entry: Fe4+
  species_attributes["Fe4+"] = species_datatype{"Fe", 4, 0, 75.00};

  // Entry: Fe5+
  species_attributes["Fe5+"] = species_datatype{"Fe", 5, 0, 98.985};

  // Entry: Fe6+
  species_attributes["Fe6+"] = species_datatype{"Fe", 6, 0, 124.9671};

  // Entry: Fe7+
  species_attributes["Fe7+"] = species_datatype{"Fe", 7, 0, 151.060};

  // Entry: Fe8+
  species_attributes["Fe8+"] = species_datatype{"Fe", 8, 0, 233.6};

  // Entry: Fe9+
  species_attributes["Fe9+"] = species_datatype{"Fe", 9, 0, 262.10};

  // Entry: Fe10+
  species_attributes["Fe10+"] = species_datatype{"Fe", 10, 0, 290.9};

  // Entry: Fe11+
  species_attributes["Fe11+"] = species_datatype{"Fe", 11, 0, 330.8};

  // Entry: Fe12+
  species_attributes["Fe12+"] = species_datatype{"Fe", 12, 0, 361.0};

  // Entry: Fe13+
  species_attributes["Fe13+"] = species_datatype{"Fe", 13, 0, 392.2};

  // Entry: Fe14+
  species_attributes["Fe14+"] = species_datatype{"Fe", 14, 0, 456.2};

  // Entry: Fe15+
  species_attributes["Fe15+"] = species_datatype{"Fe", 15, 0, 489.312};

  // Entry: Fe16+
  species_attributes["Fe16+"] = species_datatype{"Fe", 16, 0, 1262.7};

  // Entry: Fe17+
  species_attributes["Fe17+"] = species_datatype{"Fe", 17, 0, 1357.8};

  // Entry: Fe18+
  species_attributes["Fe18+"] = species_datatype{"Fe", 18, 0, 1460};

  // Entry: Fe19+
  species_attributes["Fe19+"] = species_datatype{"Fe", 19, 0, 1575.6};

  // Entry: Fe20+
  species_attributes["Fe20+"] = species_datatype{"Fe", 20, 0, 1687.0};

  // Entry: Fe21+
  species_attributes["Fe21+"] = species_datatype{"Fe", 21, 0, 1798.4};

  // Entry: Fe22+
  species_attributes["Fe22+"] = species_datatype{"Fe", 22, 0, 1950.4};

  // Entry: Fe23+
  species_attributes["Fe23+"] = species_datatype{"Fe", 23, 0, 2045.759};

  // Entry: Fe24+
  species_attributes["Fe24+"] = species_datatype{"Fe", 24, 0, 8828.1879};

  // Entry: Fe25+
  species_attributes["Fe25+"] = species_datatype{"Fe", 25, 0, 9277.6886};

  // Entry: Fe26+
  species_attributes["Fe26+"] = species_datatype{"Fe", 26, 0, 1.0e-99};
}


// **************************************************************
