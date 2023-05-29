/**
 *  \file      coll_ionise_recomb.cpp
 *  \brief     Collisional ionisation and recombination class definition.
 *  \author    Arun Mathew
 *  \date      19-12-2022
 **/

#include "coll_ionise_recomb.h"
#include <iostream>


// ##################################################################
// Constructor
coll_ionise_recomb::coll_ionise_recomb() : atomic_physics_database()

{
  collisional_ionisation_database();
  recombination_database();
}
// End of constructor ###############################################



// ##################################################################
// Destructor
coll_ionise_recomb::~coll_ionise_recomb() {}
// End of destructor ################################################



// ##################################################################
// SETUP COLLISIONAL IONISATION AND RECOMBINATION CLASS
void coll_ionise_recomb::setup_collisional_ionisation(
    const std::vector<std::string> &ci_tracer_list)
{

  // Resize CI rate table based on ci tracer list.
  collisional_rate_table.resize(ci_tracer_list.size());

  // Resize CI rate slope table based on ci tracer list.
  collisional_slope_table.resize(ci_tracer_list.size());

  // ionisation potential of species in ci_tracer_list
  for (int s = 0; s < ci_tracer_list.size(); s++)
    ci_ion_pot.push_back(get_ionisation_potential(ci_tracer_list[s]));

  // make a copy of ci_tracer_list for this class.
  ci_species_list = ci_tracer_list;
}
// END OF SETUP COLLISIONAL IONISATION  ##############################



// ##################################################################
// SETUP RECOMBINATION
void coll_ionise_recomb::setup_recombination(
    const std::vector<std::string> &recomb_tracer_list)
{

  // Resize recombination rate table based on recomb tracer list.
  recombination_rate_table.resize(recomb_tracer_list.size());

  // Resize recombination rate slope table based on MPv10 tracer list.
  recombination_slope_table.resize(recomb_tracer_list.size());

  recomb_species_list = recomb_tracer_list;
}
// END OF SETUP RECOMBINATION #######################################


// ##################################################################
// Generate LUT for collisional ionisation rate and rate slope
void coll_ionise_recomb::generate_collisional_ionisation_LUT(
    const std::vector<double> &T_table)
{
  spdlog::info("Setting up collisional ionisation LUTs");
  // Generate collisional ionisation rate table for all species
  // in the species list.
  // loop over ci species list
  for (int s = 0; s < collisional_rate_table.size(); s++) {
    for (int k = 0; k < T_table.size(); k++) {  // loop over T_table
      collisional_rate_table[s].push_back(
          get_collisional_ionisation_rate(ci_species_list[s], T_table[k]));
    }
  }


  // Generate collisional ionisation  rate slope table for all
  // species in the species list.
  for (int s = 0; s < collisional_slope_table.size(); s++) {
    for (int k = 0; k < T_table.size() - 1; k++) {
      collisional_slope_table[s].push_back(
          (collisional_rate_table[s][k + 1] - collisional_rate_table[s][k])
          / (T_table[k + 1] - T_table[k]));
    }
    collisional_slope_table[s].push_back(0.0);
  }
}
// End of generate LUT for CI rate and rate slope ##################


// ##################################################################
// Generate LUT for recombination rate and rate slope
void coll_ionise_recomb::generate_recombination_LUT(
    const std::vector<double> &T_table)
{
  spdlog::info("Setting up recombination LUTs");
  // Generate recombination rate table for all species in the species
  // list.
  for (int s = 0; s < recombination_rate_table.size(); s++) {
    for (int k = 0; k < T_table.size(); k++) {
      recombination_rate_table[s].push_back(
          get_recombination_rate(recomb_species_list[s], T_table[k]));
    }
  }

  // Generate recombination rate slope table for all species in
  // the species list.
  for (int s = 0; s < recombination_slope_table.size(); s++) {
    for (int k = 0; k < T_table.size() - 1; k++) {
      recombination_slope_table[s].push_back(
          (recombination_rate_table[s][k + 1] - recombination_rate_table[s][k])
          / (T_table[k + 1] - T_table[k]));
    }
    recombination_slope_table[s].push_back(0.0);
  }
}
// End of generate LUT for recombination rate and rate slope #######



// ##################################################################
// Calculate collisional ionisation rate of a species for a given
// temperature T.
double coll_ionise_recomb::get_collisional_ionisation_rate(
    const std::string species, const double T)
{
  double k_B  = pconst.kB();  // Boltzmann constant.
  double U    = 0.0;
  double rate = 0.0;
  U           = get_ionisation_potential(species) / k_B / T;
  rate        = coll_ionise_parameters[species].A
         * (1. + coll_ionise_parameters[species].P * sqrt(U))
         * exp(coll_ionise_parameters[species].K * log(U) - U)
         / (coll_ionise_parameters[species].X + U);

  return rate;
}
// End of get collisional ionisation rate ***************************


// ##################################################################
// Calculate recombination rate of a species for a given temperature.
double coll_ionise_recomb::get_recombination_rate(
    const std::string species, const double T)
{
  double rate = 0.0;

  // Calculating di-electronic recombination rate:
  for (int i = 0; i < 9; i++) {
    rate += recomb_parameters[species].C_coefs[i]
            * exp(-recomb_parameters[species].E_coefs[i] / T);
  }

  rate *= pow(T, -1.5);

  // Calculating radiative recombination rate
  double B_new = recomb_parameters[species].B;
  if (recomb_parameters[species].C > 0.0) {
    B_new = recomb_parameters[species].B
            + recomb_parameters[species].C
                  * exp(-recomb_parameters[species].T2 / T);
  }

  rate += recomb_parameters[species].A
          * pow(
              pow(T / recomb_parameters[species].T0, 0.5)
                  * pow(
                      1.0 + pow(T / recomb_parameters[species].T0, 0.5),
                      1.0 - B_new)
                  * pow(
                      1.0 + pow(T / recomb_parameters[species].T1, 0.5),
                      1.0 + B_new),
              -1.0);

  return rate;
}
// End of get recombination rate ###################################



// ##################################################################
// Database for collisional ionisation fitting parameters
void coll_ionise_recomb::collisional_ionisation_database()
{

  spdlog::info("Setting up collisional ionisation database");

  // This uses parametric fitting formulae from Voronov (1997)
  // Atomic data and nuclear data tables, 65, 1.
  // rate returned in cm^3/s
  // Parameters for this fit formulae are given below for each species.

  // HYDROGEN ###################################################
  coll_ionise_parameters["H"] = coll_ionise_datatype{0, 2.91e-8, 0.232, 0.39};

  coll_ionise_parameters["H1+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // HELIUM ###################################################
  coll_ionise_parameters["He"] = coll_ionise_datatype{0, 1.75e-8, 0.180, 0.35};

  coll_ionise_parameters["He1+"] =
      coll_ionise_datatype{1, 2.05e-9, 0.265, 0.25};

  coll_ionise_parameters["He2+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // CARBON ###################################################
  coll_ionise_parameters["C"] = coll_ionise_datatype{0, 0.685e-7, 0.193, 0.25};

  coll_ionise_parameters["C1+"] =
      coll_ionise_datatype{1, 0.186e-7, 0.286, 0.24};

  coll_ionise_parameters["C2+"] =
      coll_ionise_datatype{1, 0.635e-8, 0.427, 0.21};

  coll_ionise_parameters["C3+"] =
      coll_ionise_datatype{1, 0.150e-8, 0.416, 0.13};

  coll_ionise_parameters["C4+"] =
      coll_ionise_datatype{1, 0.299e-9, 0.666, 0.02};

  coll_ionise_parameters["C5+"] =
      coll_ionise_datatype{1, 0.123e-9, 0.620, 0.16};

  coll_ionise_parameters["C6+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // NITROGEN ###################################################
  coll_ionise_parameters["N"] = coll_ionise_datatype{0, 0.482e-7, 0.0652, 0.42};

  coll_ionise_parameters["N1+"] =
      coll_ionise_datatype{0, 0.298e-7, 0.310, 0.30};

  coll_ionise_parameters["N2+"] =
      coll_ionise_datatype{1, 0.810e-8, 0.350, 0.24};

  coll_ionise_parameters["N3+"] =
      coll_ionise_datatype{1, 0.371e-8, 0.549, 0.18};

  coll_ionise_parameters["N4+"] =
      coll_ionise_datatype{0, 0.151e-8, 0.0167, 0.74};

  coll_ionise_parameters["N5+"] =
      coll_ionise_datatype{0, 0.371e-9, 0.546, 0.29};

  coll_ionise_parameters["N6+"] =
      coll_ionise_datatype{1, 0.777e-10, 0.624, 0.16};

  coll_ionise_parameters["N7+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // OXYGEN ###################################################
  coll_ionise_parameters["O"] = coll_ionise_datatype{0, 0.359e-7, 0.073, 0.34};

  coll_ionise_parameters["O1+"] =
      coll_ionise_datatype{1, 0.139e-7, 0.212, 0.22};

  coll_ionise_parameters["O2+"] =
      coll_ionise_datatype{1, 0.931e-8, 0.270, 0.27};

  coll_ionise_parameters["O3+"] =
      coll_ionise_datatype{0, 0.102e-7, 0.614, 0.27};

  coll_ionise_parameters["O4+"] =
      coll_ionise_datatype{1, 0.219e-8, 0.630, 0.17};

  coll_ionise_parameters["O5+"] =
      coll_ionise_datatype{0, 0.195e-8, 0.360, 0.54};

  coll_ionise_parameters["O6+"] =
      coll_ionise_datatype{0, 0.212e-9, 0.396, 0.35};

  coll_ionise_parameters["O7+"] =
      coll_ionise_datatype{1, 0.521e-10, 0.629, 0.16};

  coll_ionise_parameters["O8+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // NEON ##################################################
  coll_ionise_parameters["Ne"] =
      coll_ionise_datatype{1, 0.150E-07, 0.0329, 0.43};

  coll_ionise_parameters["Ne1+"] =
      coll_ionise_datatype{0, 0.198E-07, 0.295, 0.20};

  coll_ionise_parameters["Ne2+"] =
      coll_ionise_datatype{1, 0.703E-08, 0.0677, 0.39};

  coll_ionise_parameters["Ne3+"] =
      coll_ionise_datatype{1, 0.424E-08, 0.0482, 0.58};

  coll_ionise_parameters["Ne4+"] =
      coll_ionise_datatype{1, 0.279E-08, 0.305, 0.25};

  coll_ionise_parameters["Ne5+"] =
      coll_ionise_datatype{0, 0.345E-08, 0.581, 0.28};

  coll_ionise_parameters["Ne6+"] =
      coll_ionise_datatype{1, 0.956E-09, 0.749, 0.14};

  coll_ionise_parameters["Ne7+"] =
      coll_ionise_datatype{1, 0.473E-09, 0.992, 0.04};

  coll_ionise_parameters["Ne8+"] =
      coll_ionise_datatype{1, 0.392E-10, 0.262, 0.20};

  coll_ionise_parameters["Ne9+"] =
      coll_ionise_datatype{1, 0.277E-10, 0.661, 0.13};

  coll_ionise_parameters["Ne10+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // SILICON ###################################################
  coll_ionise_parameters["Si"] =
      coll_ionise_datatype{1, 0.188E-06, 0.376, 0.25};

  coll_ionise_parameters["Si1+"] =
      coll_ionise_datatype{1, 0.643E-07, 0.632, 0.20};

  coll_ionise_parameters["Si2+"] =
      coll_ionise_datatype{1, 0.201E-07, 0.473, 0.22};

  coll_ionise_parameters["Si3+"] =
      coll_ionise_datatype{1, 0.494E-08, 0.172, 0.23};

  coll_ionise_parameters["Si4+"] =
      coll_ionise_datatype{1, 0.176E-08, 0.102, 0.31};

  coll_ionise_parameters["Si5+"] =
      coll_ionise_datatype{1, 0.174E-08, 0.180, 0.29};

  coll_ionise_parameters["Si6+"] =
      coll_ionise_datatype{1, 0.123E-08, 0.518, 0.07};

  coll_ionise_parameters["Si7+"] =
      coll_ionise_datatype{1, 0.827E-09, 0.239, 0.28};

  coll_ionise_parameters["Si8+"] =
      coll_ionise_datatype{1, 0.601E-09, 0.305, 0.25};

  coll_ionise_parameters["Si9+"] =
      coll_ionise_datatype{1, 0.465E-09, 0.666, 0.04};

  coll_ionise_parameters["Si10+"] =
      coll_ionise_datatype{1, 0.263E-09, 0.666, 0.16};

  coll_ionise_parameters["Si11+"] =
      coll_ionise_datatype{1, 0.118E-09, 0.734, 0.16};

  coll_ionise_parameters["Si12+"] =
      coll_ionise_datatype{0, 0.336E-10, 0.336, 0.37};

  coll_ionise_parameters["Si13+"] =
      coll_ionise_datatype{1, 0.119E-10, 0.989, 0.08};

  coll_ionise_parameters["Si14+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};

  // SULFUR ###################################################
  coll_ionise_parameters["S"] = coll_ionise_datatype{1, 0.549E-07, 0.100, 0.25};

  coll_ionise_parameters["S1+"] =
      coll_ionise_datatype{1, 0.681E-07, 0.693, 0.21};

  coll_ionise_parameters["S2+"] =
      coll_ionise_datatype{1, 0.214E-07, 0.353, 0.24};

  coll_ionise_parameters["S3+"] =
      coll_ionise_datatype{1, 0.166E-07, 1.030, 0.14};

  coll_ionise_parameters["S4+"] =
      coll_ionise_datatype{1, 0.612E-08, 0.580, 0.19};

  coll_ionise_parameters["S5+"] =
      coll_ionise_datatype{1, 0.133E-08, 0.0688, 0.35};

  coll_ionise_parameters["S6+"] =
      coll_ionise_datatype{0, 0.493E-08, 1.130, 0.16};

  coll_ionise_parameters["S7+"] =
      coll_ionise_datatype{1, 0.873E-09, 0.193, 0.28};

  coll_ionise_parameters["S8+"] =
      coll_ionise_datatype{0, 0.135E-08, 0.431, 0.32};

  coll_ionise_parameters["S9+"] =
      coll_ionise_datatype{1, 0.459E-09, 0.242, 0.28};

  coll_ionise_parameters["S10+"] =
      coll_ionise_datatype{1, 0.349E-09, 0.305, 0.25};

  coll_ionise_parameters["S11+"] =
      coll_ionise_datatype{0, 0.523E-09, 0.428, 0.35};

  coll_ionise_parameters["S12+"] =
      coll_ionise_datatype{0, 0.259E-09, 0.854, 0.12};

  coll_ionise_parameters["S13+"] =
      coll_ionise_datatype{1, 0.750E-10, 0.734, 0.16};

  coll_ionise_parameters["S14+"] =
      coll_ionise_datatype{0, 0.267E-10, 0.572, 0.28};

  coll_ionise_parameters["S15+"] =
      coll_ionise_datatype{1, 0.632E-11, 0.585, 0.17};

  coll_ionise_parameters["S16+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};


  // IRON ###################################################
  coll_ionise_parameters["Fe"] =
      coll_ionise_datatype{0, 0.252E-06, 0.701, 0.25};

  coll_ionise_parameters["Fe1+"] =
      coll_ionise_datatype{1, 0.221E-07, 0.033, 0.45};

  coll_ionise_parameters["Fe2+"] =
      coll_ionise_datatype{0, 0.410E-07, 0.366, 0.17};

  coll_ionise_parameters["Fe3+"] =
      coll_ionise_datatype{0, 0.353E-07, 0.243, 0.39};

  coll_ionise_parameters["Fe4+"] =
      coll_ionise_datatype{1, 0.104E-07, 0.285, 0.17};

  coll_ionise_parameters["Fe5+"] =
      coll_ionise_datatype{1, 0.123E-07, 0.411, 0.21};

  coll_ionise_parameters["Fe6+"] =
      coll_ionise_datatype{1, 0.947E-08, 0.458, 0.21};

  coll_ionise_parameters["Fe7+"] =
      coll_ionise_datatype{1, 0.471E-08, 0.280, 0.28};

  coll_ionise_parameters["Fe8+"] =
      coll_ionise_datatype{1, 0.302E-08, 0.697, 0.15};

  coll_ionise_parameters["Fe9+"] =
      coll_ionise_datatype{1, 0.234E-08, 0.764, 0.14};

  coll_ionise_parameters["Fe10+"] =
      coll_ionise_datatype{1, 0.176E-08, 0.805, 0.14};

  coll_ionise_parameters["Fe11+"] =
      coll_ionise_datatype{1, 0.114E-08, 0.773, 0.15};

  coll_ionise_parameters["Fe12+"] =
      coll_ionise_datatype{1, 0.866E-09, 0.805, 0.14};

  coll_ionise_parameters["Fe13+"] =
      coll_ionise_datatype{1, 0.661E-09, 0.762, 0.14};

  coll_ionise_parameters["Fe14+"] =
      coll_ionise_datatype{1, 0.441E-09, 0.698, 0.16};

  coll_ionise_parameters["Fe15+"] =
      coll_ionise_datatype{1, 0.118E-09, 0.211, 0.15};

  coll_ionise_parameters["Fe16+"] =
      coll_ionise_datatype{1, 0.361E-09, 1.160, 0.09};

  coll_ionise_parameters["Fe17+"] =
      coll_ionise_datatype{1, 0.245E-09, 0.978, 0.13};

  coll_ionise_parameters["Fe18+"] =
      coll_ionise_datatype{1, 0.187E-09, 0.988, 0.14};

  coll_ionise_parameters["Fe19+"] =
      coll_ionise_datatype{1, 0.133E-09, 1.030, 0.12};

  coll_ionise_parameters["Fe20+"] =
      coll_ionise_datatype{1, 0.784E-10, 0.848, 0.14};

  coll_ionise_parameters["Fe21+"] =
      coll_ionise_datatype{0, 0.890E-10, 1.200, 0.35};

  coll_ionise_parameters["Fe22+"] =
      coll_ionise_datatype{1, 0.229E-10, 0.936, 0.12};

  coll_ionise_parameters["Fe23+"] =
      coll_ionise_datatype{0, 0.112E-10, 0.034, 0.81};

  coll_ionise_parameters["Fe24+"] =
      coll_ionise_datatype{1, 0.246E-11, 1.02, 0.02};

  coll_ionise_parameters["Fe25+"] =
      coll_ionise_datatype{1, 0.979E-12, 0.664, 0.14};

  coll_ionise_parameters["Fe26+"] = coll_ionise_datatype{0, 0.0, 0.0, 0.0};
}
// END OF COLLISIONAL IONISATION DATABASE ################################



// ##################################################################
// Database for recombination fitting parameters
void coll_ionise_recomb::recombination_database()
{

  // Radiative Recombination Data from
  // http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/
  // badnell@phys.strath.ac.uk from August 2021.
  // Dielectronic Recombination Data from
  // http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/
  // badnell@phys.strath.ac.uk from August 2021.

  spdlog::info("Setting up recombination database");

  // Hydrogen and its ion.
  recomb_parameters["H"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["H1+"] = recomb_datatype{
      8.318E-11,
      0.7472,
      2.965E+00,
      7.001E+05,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // Helium and its ion.
  recomb_parameters["He"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["He1+"] = recomb_datatype{
      5.235E-11,
      0.6988,
      7.301E+00,
      4.475E+06,
      0.0829,
      1.682E+05,
      {1.417E-03, 2.235E-04, -2.185E-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {4.633E+05, 5.532E+05, 8.887E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["He2+"] = recomb_datatype{
      1.818E-10,
      0.7492,
      1.017E+01,
      2.786E+06,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // Carbon and its ions.
  recomb_parameters["C"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["C1+"] = recomb_datatype{
      2.995E-09,
      0.7849,
      6.670E-03,
      1.943E+06,
      0.1597,
      4.955E+04,
      {6.346E-09, 9.793E-09, 1.634E-06, 8.369E-04, 3.355E-04, 0.0, 0.0, 0.0,
       0.0},
      {1.217E+01, 7.380E+01, 1.523E+04, 1.207E+05, 2.144E+05, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["C2+"] = recomb_datatype{
      2.067E-09,
      0.8012,
      1.643E-01,
      2.172E+06,
      0.0427,
      6.341E+04,
      {3.489E-06, 2.222E-07, 1.954E-05, 4.212E-03, 2.037E-04, 2.936E-04, 0.0,
       0.0, 0.0},
      {2.660E+03, 3.756E+03, 2.566E+04, 1.400E+05, 1.801E+06, 4.307E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["C3+"] = recomb_datatype{
      1.120E-10,
      0.6737,
      1.115E+02,
      5.938E+06,
      0.0,
      0.0,
      {4.673E-07, 1.887E-05, 1.305E-05, 3.099E-03, 3.001E-04, 2.553E-03, 0.0,
       0.0, 0.0},
      {7.233E+02, 2.847E+03, 1.054E+04, 8.915E+04, 2.812E+05, 3.254E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["C4+"] = recomb_datatype{
      4.798E-11,
      0.4834,
      1.355E+03,
      1.872E+07,
      0.0,
      0.0,
      {2.646E-03, 1.762E-02, -7.843E-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {2.804E+06, 3.485E+06, 4.324E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["C5+"] = recomb_datatype{
      2.044E-10,
      0.6742,
      2.647E+02,
      2.773E+07,
      0.0,
      0.0,
      {1.426E-03, 3.046E-02, 8.373E-04, 0.0, 0.0, 0.0, 0.0, 0.0},
      {3.116E+06, 4.075E+06, 5.749E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["C6+"] = recomb_datatype{
      5.337E-10,
      0.7485,
      9.502E+01,
      2.517E+07,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // Nitrogen and its ions.
  recomb_parameters["N"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["N1+"] = recomb_datatype{
      6.387E-10,
      0.7308,
      9.467E-02,
      2.954E+06,
      0.2440,
      6.739E+04,
      {1.658E-08, 2.760E-08, 2.391E-09, 7.585E-07, 3.012E-04, 7.132E-04, 0.0,
       0.0, 0.0},
      {1.265E+01, 8.425E+01, 2.964E+02, 5.923E+03, 1.278E+05, 2.184E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["N2+"] = recomb_datatype{
      2.410E-09,
      0.7948,
      1.231E-01,
      3.016E+06,
      0.0774,
      1.016E+05,
      {7.712E-08, 4.839E-08, 2.218E-06, 1.536E-03, 3.647E-03, 4.234E-05, 0.0,
       0.0, 0.0},
      {7.113E+01, 2.765E+02, 1.439E+04, 1.347E+05, 2.496E+05, 2.204E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["N3+"] = recomb_datatype{
      7.923E-10,
      0.7768,
      3.750E+00,
      3.468E+06,
      0.0223,
      7.206E+0,
      {3.386E-06, 3.036E-05, 5.945E-05, 1.195E-03, 6.462E-03, 1.358E-03, 0.0,
       0.0, 0.0},
      {1.406E+03, 6.965E+03, 2.604E+04, 1.304E+05, 1.965E+05, 4.466E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["N4+"] = recomb_datatype{
      1.533E-10,
      0.6682,
      1.823E+02,
      7.751E+06,
      0.0,
      0.0,
      {2.040E-06, 6.986E-05, 3.168E-04, 4.353E-03, 7.765E-04, 5.101E-03, 0.0,
       0.0, 0.0},
      {3.084E+03, 1.332E+04, 6.475E+04, 1.181E+05, 6.687E+05, 4.778E+06, 0.0,
       0.0, 0.0}};


  recomb_parameters["N5+"] = recomb_datatype{
      6.245E-11,
      0.4985,
      1.957E+03,
      2.177E+07,
      0.0,
      0.0,
      {5.761E-03, 3.434E-02, -1.660E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {3.860E+06, 4.883E+06, 6.259E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};


  recomb_parameters["N6+"] = recomb_datatype{
      2.388E-10,
      0.6732,
      3.960E+02,
      3.583E+07,
      0.0,
      0.0,
      {2.801E-03, 4.362E-02, 1.117E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {4.198E+06, 5.516E+06, 8.050E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["N7+"] = recomb_datatype{
      6.170E-10,
      0.7481,
      1.316E+02,
      3.427E+07,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // Oxygen and its ions.
  recomb_parameters["O"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["O1+"] = recomb_datatype{
      6.622E-11,
      0.6109,
      4.136E+00,
      4.214E+06,
      0.4093,
      8.770E+04,
      {5.629E-08, 2.550E-07, 6.173E-04, 1.627E-04, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.395E+03, 1.770E+04, 1.671E+05, 2.687E+05, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["O2+"] = recomb_datatype{
      2.096E-09,
      0.7668,
      1.602E-01,
      4.377E+06,
      0.1070,
      1.392E+05,
      {1.627E-07, 1.262E-07, 6.663E-07, 3.925E-06, 2.406E-03, 1.146E-03, 0.0,
       0.0, 0.0},
      {4.535E+01, 2.847E+02, 4.166E+03, 2.877E+04, 1.953E+05, 3.646E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["O3+"] = recomb_datatype{
      2.501E-09,
      0.7844,
      5.235E-01,
      4.470E+06,
      0.0447,
      1.642E+05,
      {3.932E-07, 2.523E-07, 3.447E-05, 5.776E-03, 5.101E-03, 0.0, 0.0, 0.0,
       0.0},
      {1.509E+02, 6.211E+02, 1.562E+04, 1.936E+05, 4.700E+05, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["O4+"] = recomb_datatype{
      3.955E-09,
      0.7813,
      6.821E-01,
      5.076E+06,
      0.0,
      0.0,
      {1.615E-05, 9.299E-06, 1.530E-04, 6.616E-04, 1.080E-02, 7.503E-04,
       2.892E-03, 0.0, 0.0},
      {7.569E+02, 3.659E+03, 1.984E+04, 8.429E+04, 2.294E+05, 1.161E+06,
       6.137E+06, 0.0, 0.0, 0.0}};

  recomb_parameters["O5+"] = recomb_datatype{
      1.724E-10,
      0.6556,
      3.372E+02,
      1.030E+07,
      0.0,
      0.0,
      {2.389E-05, 1.355E-04, 5.885E-03, 2.163E-03, 6.341E-04, 1.348E-02, 0.0,
       0.0, 0.0},
      {2.326E+04, 3.209E+04, 1.316E+05, 6.731E+05, 1.892E+06, 6.150E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["O6+"] = recomb_datatype{
      8.193E-11,
      0.5165,
      2.392E+03,
      2.487E+07,
      0.0,
      0.0,
      {6.135E-02, 1.968E-04, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.113E+06, 3.656E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["O7+"] = recomb_datatype{
      2.652E-10,
      0.6705,
      5.842E+02,
      4.559E+07,
      0.0,
      0.0,
      {4.925E-03, 5.837E-02, 1.359E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.440E+06, 7.170E+06, 1.152E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["O8+"] = recomb_datatype{
      6.552E-10,
      0.7470,
      1.951E+02,
      4.483E+07,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  //##################################################################

  // NEON
  recomb_parameters["Ne"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Ne1+"] = recomb_datatype{
      1.295E-11,
      0.3556,
      6.739E+01,
      7.563E+06,
      0.6472,
      1.598E+05,
      {4.152E-09, 4.656E-09, 1.310E-08, 1.417E-09, 7.968E-04, 1.271E-05, 0.0,
       0.0, 0.0},
      {2.689E+01, 2.021E+02, 7.200E+02, 4.892E+04, 3.144E+05, 6.738E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["Ne2+"] = recomb_datatype{
      1.773E-10,
      0.6434,
      9.924E+00,
      8.878E+06,
      0.2205,
      2.292E+05,
      {2.980E-08, 1.257E-07, 1.122E-06, 2.626E-03, 8.802E-04, 1.231E+00, 0.0,
       0.0, 0.0},
      {4.579E+01, 4.753E+02, 1.481E+04, 2.810E+05, 4.763E+05, 4.677E+08, 0.0,
       0.0, 0.0}};

  recomb_parameters["Ne3+"] = recomb_datatype{
      8.321E-10,
      0.7254,
      3.332E+00,
      8.696E+06,
      0.0921,
      3.044E+05,
      {2.763E-06, 1.053E-05, 4.453E-05, 6.244E-03, 3.146E-04, 4.465E-03, 0.0,
       0.0, 0.0},
      {6.393E+02, 1.499E+03, 3.227E+04, 2.561E+05, 4.505E+05, 7.934E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["Ne4+"] = recomb_datatype{
      1.861E-09,
      0.7593,
      2.504E+00,
      8.037E+06,
      0.0406,
      3.255E+05,
      {2.922E-06, 7.144E-06, 2.836E-05, 9.820E-05, 8.379E-03, 1.009E-02, 0.0,
       0.0, 0.0},
      {2.050E+02, 2.205E+03, 9.271E+03, 4.988E+04, 2.904E+05, 8.782E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["Ne5+"] = recomb_datatype{
      1.127E-09,
      0.7556,
      1.311E+01,
      8.047E+06,
      0.0250,
      2.771E+05,
      {5.653E-06, 4.344E-05, 1.086E-04, 5.980E-04, 1.457E-02, 1.601E-02,
       5.365E-04, 0.0, 0.0},
      {6.280E+02, 2.812E+03, 1.324E+04, 8.064E+04, 3.052E+05, 1.032E+06,
       2.388E+06, 0.0, 0.0}};

  recomb_parameters["Ne6+"] = recomb_datatype{
      2.557E-09,
      0.7601,
      6.293E+00,
      8.091E+06,
      0.0,
      0.0,
      {8.460E-05, 1.817E-04, 1.176E-03, 1.397E-02, 5.566E-03, 6.229E-03,
       1.031E-02, 0.0, 0.0},
      {1.049E+03, 3.829E+03, 6.133E+04, 2.568E+05, 4.600E+05, 1.324E+06,
       9.353E+06, 0.0, 0.0}};

  recomb_parameters["Ne7+"] = recomb_datatype{
      2.755E-10,
      0.6586,
      5.102E+02,
      1.535E+07,
      0.0,
      0.0,
      {2.399E-04, 3.532E-04, 8.928E-03, 5.427E-03, 5.342E-03, 3.981E-02, 0.0,
       0.0, 0.0},
      {2.536E+04, 4.800E+04, 1.770E+05, 1.030E+06, 1.859E+06, 9.743E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Ne8+"] = recomb_datatype{
      1.186E-10,
      0.5354,
      3.647E+03,
      3.365E+07,
      0.0,
      0.0,
      {1.552E-02, 9.008E-02, 1.182E-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.845E+06, 9.803E+06, 1.209E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Ne9+"] = recomb_datatype{
      3.415E-10,
      0.6706,
      9.552E+02,
      6.778E+07,
      0.0,
      0.0,
      {1.183E-02, 9.011E-02, 1.828E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {8.405E+06, 1.111E+07, 1.812E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Ne10+"] = recomb_datatype{
      8.278E-10,
      0.7470,
      2.991E+02,
      7.006E+07,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // SILICON
  recomb_parameters["Si"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Si1+"] = recomb_datatype{
      3.262E-11,
      0.6270,
      1.590E+01,
      4.237E+07,
      0.2333,
      5.828E+04,
      {3.408E-08, 1.913E-07, 1.679E-07, 7.523E-07, 8.386E-05, 4.083E-03, 0.0,
       0.0, 0.0},
      {2.431E+01, 1.293E+02, 4.272E+02, 3.729E+03, 5.514E+04, 1.295E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["Si2+"] = recomb_datatype{
      1.964E-10,
      0.6287,
      7.712E+00,
      2.951E+07,
      0.1523,
      4.804E+05,
      {2.930E-06, 2.803E-06, 9.023E-05, 6.909E-03, 2.582E-05, 0.0, 0.0, 0.0,
       0.0},
      {1.162E+02, 5.721E+03, 3.477E+04, 1.176E+05, 3.505E+06, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["Si3+"] = recomb_datatype{
      6.739E-11,
      0.4931,
      2.166E+02,
      4.491E+07,
      0.1667,
      9.046E+05,
      {3.819E-06, 2.421E-05, 2.283E-04, 8.604E-03, 2.617E-03, 0.0, 0.0, 0.0,
       0.0},
      {3.802E+03, 1.280E+04, 5.953E+04, 1.026E+05, 1.154E+06, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["Si4+"] = recomb_datatype{
      5.134E-11,
      0.3678,
      1.009E+03,
      8.514E+07,
      0.1646,
      1.084E+06,
      {1.422E-04, 9.474E-03, 1.650E-030, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.685E+05, 1.208E+06, 1.839E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Si5+"] = recomb_datatype{
      2.468E-10,
      0.6113,
      1.649E+02,
      3.231E+07,
      0.0636,
      9.837E+05,
      {7.163E-07, 2.656E-06, 1.119E-06, 4.796E-05, 4.052E-03, 6.101E-03,
       2.366E-02, 0.0, 0.0},
      {5.625E+02, 2.952E+03, 9.682E+03, 1.473E+05, 5.064E+05, 8.047E+05,
       1.623E+06, 0.0, 0.0}};

  recomb_parameters["Si6+"] = recomb_datatype{
      4.615E-10,
      0.6753,
      1.143E+02,
      2.377E+07,
      0.0356,
      8.595E+05,
      {2.086E-06, 9.423E-06, 3.423E-05, 3.950E-04, 1.535E-02, 4.986E-02,
       4.067E-04, 0.0, 0.0},
      {3.708E+02, 3.870E+03, 2.226E+04, 1.318E+05, 5.285E+05, 1.742E+06,
       8.392E+06, 0.0, 0.0}};

  recomb_parameters["Si7+"] = recomb_datatype{
      7.532E-10,
      0.7072,
      8.860E+01,
      1.997E+07,
      0.0185,
      6.949E+05,
      {5.272E-05, 2.282E-04, 1.345E-03, 2.246E-02, 9.606E-02, 8.366E-04, 0.0,
       0.0, 0.0},
      {2.829E+03, 2.617E+04, 1.374E+05, 4.520E+05, 2.072E+06, 7.808E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Si8+"] = recomb_datatype{
      2.100E-09,
      0.7401,
      2.523E+01,
      1.842E+07,
      0.0,
      0.0,
      {5.845E-04, 6.600E-04, 7.180E-04, 1.714E-02, 2.958E-02, 1.075E-01, 0.0,
       0.0, 0.0},
      {9.856E+02, 6.577E+03, 4.281E+04, 3.944E+05, 1.296E+06, 2.441E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Si9+"] = recomb_datatype{
      1.688E-09,
      0.7390,
      5.549E+01,
      1.716E+07,
      0.0,
      0.0,
      {1.246E-04, 6.649E-04, 2.912E-03, 2.912E-02, 3.049E-02, 1.056E-01, 0.0,
       0.0, 0.0},
      {1.167E+03, 9.088E+03, 5.332E+04, 4.219E+05, 1.571E+06, 2.721E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Si10+"] = recomb_datatype{
      1.851E-09,
      0.7384,
      6.906E+01,
      1.644E+07,
      0.0,
      0.0,
      {2.214E-04, 1.260E-03, 3.832E-03, 1.158E-02, 2.871E-02, 5.456E-02,
       5.035E-02, 0.0, 0.0},
      {3.586E+03, 1.324E+04, 5.636E+04, 2.498E+05, 5.363E+05, 2.809E+06,
       1.854E+07, 0.0, 0.0}};

  recomb_parameters["Si11+"] = recomb_datatype{
      4.633E-10,
      0.6602,
      1.088E+03,
      2.896E+07,
      0.0,
      0.0,
      {1.205E-03, 1.309E-02, 5.333E-03, 2.858E-02, 3.195E-02, 1.433E-01, 0.0,
       0.0, 0.0},
      {4.824E+04, 2.137E+05, 5.492E+05, 2.460E+06, 3.647E+06, 1.889E+07, 0.0,
       0.0, 0.0}};

  recomb_parameters["Si12+"] = recomb_datatype{
      2.017E-10,
      0.5588,
      6.494E+03,
      5.693E+07,
      0.0,
      0.0,
      {5.318E-02, 1.874E-01, 1.227E-02, 7.173E-04, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.552E+07, 1.969E+07, 2.532E+07, 2.696E+08, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Si13+"] = recomb_datatype{
      4.870E-10,
      0.6697,
      2.026E+03,
      1.265E+08,
      0.0,
      0.0,
      {3.846E-02, 1.491E-01, 2.779E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.627E+07, 2.154E+07, 3.827E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Si14+"] = recomb_datatype{
      1.261E-09,
      0.7488,
      5.068E+02,
      1.365E+08,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // SULFUR
  recomb_parameters["S"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["S1+"] = recomb_datatype  // DATA NOT AVAILABLE
      {7.50E-12,
       0.645,
       1.75E+02,
       3.70E+08,
       0.0,
       0.0,
       {1.62E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.25E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["S2+"] = recomb_datatype{
      2.478E-11,
      0.4642,
      3.294E+02,
      2.166E+07,
      0.3351,
      7.630E+05,
      {3.040E-07, 4.393E-07, 1.609E-06, 4.980E-06, 3.457E-05, 8.617E-03,
       9.284E-04, 0.0, 0.0},
      {5.016E+01, 3.266E+02, 3.102E+03, 1.210E+04, 4.969E+04, 2.010E+05,
       2.575E+05, 0.0, 0.0}};

  recomb_parameters["S3+"] = recomb_datatype{
      3.043E-10,
      0.6947,
      1.678E+01,
      2.050E+07,
      0.0795,
      6.868E+04,
      {5.817E-07, 1.391E-06, 1.123E-05, 1.521E-04, 1.875E-03, 2.097E-02, 0.0,
       0.0, 0.0},
      {3.628E+02, 1.058E+03, 7.160E+03, 3.260E+04, 1.235E+05, 2.070E+05, 0.0,
       0.0, 0.0}};

  recomb_parameters["S4+"] = recomb_datatype{
      2.615E-10,
      0.6343,
      6.238E+01,
      2.803E+07,
      0.0773,
      1.059E+06,
      {9.571E-06, 6.268E-05, 3.807E-04, 1.874E-02, 5.526E-03, 0.0, 0.0, 0.0,
       0.0},
      {1.180E+03, 6.443E+03, 2.264E+04, 1.530E+05, 3.564E+05, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["S5+"] = recomb_datatype{
      1.588E-10,
      0.5584,
      3.350E+02,
      5.188E+07,
      0.0591,
      1.656E+06,
      {2.816E-06, 3.172E-05, 1.832E-04, 4.360E-03, 1.618E-02, 7.707E-03, 0.0,
       0.0, 0.0},
      {7.590E+03, 1.558E+04, 4.013E+04, 1.156E+05, 1.601E+05, 1.839E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["S6+"] = recomb_datatype{
      9.565E-11,
      0.4517,
      1.599E+03,
      9.252E+07,
      0.0612,
      1.986E+06,
      {3.931E-05, 4.431E-03, 5.156E-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {9.455E+05, 1.365E+06, 2.169E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["S7+"] = recomb_datatype{
      3.384E-10,
      0.6175,
      3.380E+02,
      4.347E+07,
      0.0225,
      1.672E+06,
      {2.585E-06, 9.517E-06, 5.194E-06, 3.715E-04, 1.553E-02, 1.013E-01, 0.0,
       0.0, 0.0},
      {1.277E+03, 5.978E+03, 2.097E+04, 2.292E+05, 7.974E+05, 2.376E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["S8+"] = recomb_datatype{
      8.773E-10,
      0.6853,
      1.115E+02,
      3.386E+07,
      0.0,
      0.0,
      {2.924E-05, 3.366E-04, 2.104E-04, 1.910E-02, -4.017E-04, 6.541E-02,
       9.546E-02, 0.0, 0.0},
      {7.192E+02, 7.490E+03, 4.510E+04, 5.586E+05, 3.830E+05, 1.976E+06,
       2.894E+06, 0.0, 0.0}};

  recomb_parameters["S9+"] = recomb_datatype{
      1.137E-09,
      0.7080,
      1.099E+02,
      2.745E+07,
      0.0,
      0.0,
      {1.830E-04, 1.551E-03, 1.717E-03, 2.798E-02, 6.933E-02, 1.727E-01, 0.0,
       0.0, 0.0},
      {1.122E+04, 3.121E+04, 9.348E+04, 4.792E+05, 2.052E+06, 3.291E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["S10+"] = recomb_datatype{
      1.518E-09,
      0.7246,
      9.845E+01,
      2.390E+07,
      0.0,
      0.0,
      {1.371E-04, 3.096E-04, 1.782E-03, 1.278E-02, 3.323E-02, 2.593E-01, 0.0,
       0.0, 0.0},
      {1.505E+03, 1.100E+04, 6.720E+04, 3.364E+05, 9.686E+05, 3.195E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["S11+"] = recomb_datatype{
      1.702E-09,
      0.7301,
      1.138E+02,
      2.253E+07,
      0.0,
      0.0,
      {9.126E-05, 1.691E-04, 3.050E-03, 2.604E-02, 3.245E-02, 2.511E-01,
       4.459E-04, 0.0, 0.0},
      {3.575E+03, 1.491E+04, 9.388E+04, 3.856E+05, 9.630E+05, 3.490E+06,
       1.743E+08, 0.0, 0.0}};

  recomb_parameters["S12+"] = recomb_datatype{
      1.740E-09,
      0.7303,
      1.494E+02,
      2.193E+07,
      0.0,
      0.0,
      {2.610E-04, 9.442E-04, 8.190E-03, 4.271E-02, 1.997E-02, 9.590E-02,
       8.444E-02, 0.0, 0.0},
      {5.114E+03, 1.750E+04, 1.061E+05, 4.691E+05, 1.821E+06, 4.033E+06,
       2.476E+07, 0.0, 0.0}};

  recomb_parameters["S13+"] = recomb_datatype{
      5.511E-10,
      0.6598,
      1.492E+03,
      3.755E+07,
      0.0,
      0.0,
      {1.173E-03, 1.718E-03, 1.657E-02, 7.474E-03, 9.698E-02, 1.399E-01,
       7.029E-02, 0.0, 0.0},
      {1.366E+04, 6.371E+04, 2.540E+05, 6.493E+05, 3.868E+06, 2.186E+07,
       2.983E+07, 0.0, 0.0}};

  recomb_parameters["S14+"] = recomb_datatype{
      2.362E-10,
      0.5615,
      8.776E+03,
      7.208E+07,
      0.0,
      0.0,
      {8.410E-02, 2.381E-01, 1.065E-02, 1.049E-03, 0.0, 0.0, 0.0, 0.0, 0.0},
      {2.032E+07, 2.592E+07, 3.206E+07, 2.016E+08, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["S15+"] = recomb_datatype{
      5.546E-10,
      0.6692,
      2.754E+03,
      1.633E+08,
      0.0,
      0.0,
      {6.659E-02, 1.762E-01, -6.522E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {2.122E+07, 2.897E+07, 5.786E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["S16+"] = recomb_datatype{
      1.432E-09,
      0.7485,
      6.688E+02,
      1.793E+08,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  // IRON
  recomb_parameters["Fe"] = recomb_datatype{
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe1+"] = recomb_datatype{
      8.945E-09,
      0.2156,
      4.184E-02,
      5.353E+13,
      0.0,
      0.0,
      {0.17609E-09, 0.80041E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.942E+04, 1.497E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe2+"] = recomb_datatype{
      6.076E-09,
      0.3112,
      3.401E-01,
      1.960E+12,
      0.0,
      0.0,
      {0.18409E-08, 0.21611E-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.938E+05, 3.644E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe3+"] = recomb_datatype{
      5.422E-09,
      0.5067,
      4.998E-01,
      8.079E+11,
      0.0,
      0.0,
      {0.12006E-07, 0.37619E-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {3.319E+05, 6.046E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe4+"] = recomb_datatype{
      6.330E-09,
      0.6355,
      5.457E-01,
      8.545E+11,
      0.0,
      0.0,
      {0.30416E-07, 0.12807E-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {4.328E+05, 7.821E+05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe5+"] = recomb_datatype{
      8.288E-09,
      0.6012,
      1.216E+00,
      1.182E+12,
      0.0,
      0.0,
      {0.64033E-07, 0.19210E-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.290E+05, 1.160E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe6+"] = recomb_datatype{
      1.036E-08,
      0.5428,
      2.743E+00,
      1.014E+12,
      0.0,
      0.0,
      {0.73638E-07, 0.32817E-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.280E+05, 4.178E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe7+"] = recomb_datatype{
      1.184E-08,
      0.4798,
      5.883E+00,
      2.582E+12,
      0.0,
      0.0,
      {0.12807E-06, 0.28815E-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.740E+05, 1.427E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe8+"] = recomb_datatype{
      3.341E-10,
      0.4865,
      1.891E+03,
      5.181E+07,
      0.0575,
      2.734E+04,
      {3.786E-04, 1.494E-04, 9.447E-04, 2.100E-03, 3.520E-03, 2.785E-02,
       2.734E-01, 0.0, 0.0},
      {2.021E+03, 2.187E+03, 7.063E+03, 2.074E+04, 4.974E+04, 2.743E+05,
       8.348E+05, 0.0, 0.0}};

  recomb_parameters["Fe9+"] = recomb_datatype{
      1.907E-08,
      0.3768,
      1.216E+01,
      5.431E+12,
      0.0,
      0.0,
      {0.11206E-06, 0.20811E-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {2.507E+05, 1.578E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe10+"] = recomb_datatype{
      2.460E-08,
      0.3387,
      1.487E+01,
      5.228E+12,
      0.0,
      0.0,
      {0.80041E-07, 0.22412E-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {2.576E+05, 1.671E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe11+"] = recomb_datatype{
      2.165E-08,
      0.3403,
      2.195E+01,
      6.383E+12,
      0.0,
      0.0,
      {0.18009E-06, 0.18489E-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.916E+05, 4.201E+06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe12+"] = recomb_datatype{
      1.984E-09,
      0.7101,
      1.158E+02,
      4.400E+07,
      0.0,
      0.0,
      {4.469E-03, 8.538E-03, 1.741E-02, 1.630E-01, 8.680E-02, 0.0, 0.0, 0.0,
       0.0},
      {2.462E+03, 1.261E+04, 9.330E+04, 4.887E+05, 1.312E+06, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["Fe13+"] = recomb_datatype{
      1.121E-09,
      0.6984,
      4.071E+02,
      4.264E+07,
      0.0,
      0.0,
      {1.753E-03, 1.038E-02, 2.573E-02, 1.189E-01, 1.070E-01, 3.080E-02,
       6.324E-04, 0.0, 0.0},
      {4.287E+03, 1.679E+04, 9.992E+04, 3.841E+05, 6.758E+05, 1.476E+06,
       3.334E+08, 0.0, 0.0}};

  recomb_parameters["Fe14+"] = recomb_datatype{
      5.398E-10,
      0.6295,
      1.881E+03,
      5.429E+07,
      0.0,
      0.0,
      {5.636E-04, 7.860E-03, 5.063E-02, 1.753E-01, 1.209E-01, 1.934E-01, 0.0,
       0.0, 0.0},
      {3.628E+03, 2.489E+04, 1.405E+05, 5.133E+05, 5.018E+06, 8.689E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Fe15+"] = recomb_datatype{
      3.133E-10,
      0.5507,
      6.295E+03,
      9.035E+07,
      0.0,
      0.0,
      {7.676E-04, 5.587E-03, 1.152E-01, 4.929E-02, 7.274E-01, 7.347E-03, 0.0,
       0.0, 0.0},
      {2.935E+04, 8.158E+04, 3.591E+05, 1.735E+06, 7.545E+06, 4.634E+07, 0.0,
       0.0, 0.0}};

  recomb_parameters["Fe16+"] = recomb_datatype{
      2.034E-10,
      0.4548,
      1.751E+04,
      1.579E+08,
      0.0,
      0.0,
      {6.342E-04, 8.350E-02, 1.045E+00, 3.663E-01, -3.955E-02, 0.0, 0.0, 0.0,
       0.0},
      {2.770E+06, 3.978E+06, 7.052E+06, 1.300E+07, 3.579E+07, 0.0, 0.0, 0.0,
       0.0}};

  recomb_parameters["Fe17+"] = recomb_datatype{
      4.791E-10,
      0.5823,
      4.967E+03,
      9.535E+07,
      0.0,
      0.0,
      {8.207E-05, 2.766E-04, 1.897E-03, 2.842E-02, 4.022E-01, 1.434E+00, 0.0,
       0.0, 0.0},
      {6.365E+03, 3.842E+04, 2.002E+05, 1.150E+06, 4.736E+06, 8.891E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Fe18+"] = recomb_datatype{
      7.556E-10,
      0.6351,
      2.800E+03,
      7.742E+07,
      0.0,
      0.0,
      {2.033E-04, 1.159E-03, 5.567E-03, 5.482E-02, 3.370E-01, 1.518E+00, 0.0,
       0.0, 0.0},
      {3.502E+03, 3.557E+04, 2.177E+05, 1.078E+06, 4.515E+06, 9.050E+06, 0.0,
       0.0, 0.0}};

  recomb_parameters["Fe19+"] = recomb_datatype{
      1.135E-09,
      0.6705,
      1.691E+03,
      6.809E+07,
      0.0,
      0.0,
      {2.106E-03, 6.569E-03, 1.532E-02, 3.799E-02, 7.669E-02, 6.701E-01,
       1.298E+00, 0.0, 0.0},
      {4.463E+03, 3.545E+04, 1.944E+05, 6.148E+05, 1.635E+06, 6.100E+06,
       1.026E+07, 0.0, 0.0}};

  recomb_parameters["Fe20+"] = recomb_datatype{
      1.659E-09,
      0.6958,
      1.061E+03,
      6.253E+07,
      0.0,
      0.0,
      {2.565E-03, 1.685E-02, 1.827E-02, 6.957E-02, 3.254E-01, 5.101E-01,
       7.325E-01, 0.0, 0.0},
      {6.553E+03, 4.351E+04, 2.059E+05, 9.250E+05, 4.610E+06, 1.019E+07,
       1.019E+07, 0.0, 0.0}};


  recomb_parameters["Fe21+"] = recomb_datatype{
      2.199E-09,
      0.7118,
      7.810E+02,
      5.946E+070,
      0.0,
      0.0,
      {8.382E-03, 7.897E-03, 3.157E-02, 1.159E-01, 3.919E-01, 1.017E+00, 0.0,
       0.0, 0.0},
      {5.297E+03, 5.829E+04, 2.454E+05, 9.663E+05, 5.580E+06, 1.110E+07, 0.0,
       0.0, 0.0}};

  recomb_parameters["Fe22+"] = recomb_datatype{
      3.322E-09,
      0.7264,
      4.563E+02,
      5.746E+07,
      0.0,
      0.0,
      {4.954E-03, 8.183E-03, 9.028E-03, 2.566E-02, 3.470E-01, 7.455E-01,
       3.559E-01, 0.0, 0.0},
      {1.376E+04, 8.251E+04, 2.794E+05, 9.378E+05, 4.688E+06, 1.106E+07,
       6.543E+07, 0.0, 0.0}};

  recomb_parameters["Fe23+"] = recomb_datatype{
      1.186E-09,
      0.6713,
      3.253E+03,
      9.392E+07,
      0.0,
      0.0,
      {7.882E-03, 1.636E-02, 4.868E-02, 4.230E-02, 4.151E-01, 5.339E-01,
       4.544E-03, 0.0, 0.0},
      {5.977E+04, 2.074E+05, 6.170E+05, 3.789E+06, 1.086E+07, 6.363E+07,
       1.599E+08, 0.0, 0.0}};

  recomb_parameters["Fe24+"] = recomb_datatype{
      4.458E-10,
      0.5802,
      2.155E+04,
      1.701E+08,
      0.0,
      0.0,
      {2.676E-01, 4.097E-01, 2.990E-02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.394E+07, 6.854E+07, 9.651E+07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe25+"] = recomb_datatype{
      9.983E-10,
      0.6754,
      6.651E+03,
      4.017E+08,
      0.0,
      0.0,
      {1.984E-01, 2.676E-01, -2.293E-03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {5.552E+07, 7.475E+07, 1.236E+08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

  recomb_parameters["Fe26+"] = recomb_datatype{
      2.275E-09,
      0.7481,
      1.836E+03,
      4.736E+08,
      0.0,
      0.0,
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
}
// End of database for recombination fitting parameters ################