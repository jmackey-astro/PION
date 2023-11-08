/// \file multi_ion_photo_xsections.cpp
/// \author Arun Mathew
/// \date 2023.06.17
///
/// Description:
/// This class calculates the photoionisation cross-sections of multi-ion
/// module species using Verner-Ferland fit formula
///
/// References:
/// - VF96: Verner and Ferland (1996,ApJ,465,187) Table 1 & Eqn 1

#include "multi_ion_photo_xsections.h"



// ##################################################################
// Constructor
multi_ion_photo_xsections::multi_ion_photo_xsections()
{
  multi_ion_photo_xsections_database();
}


// ##################################################################
// Destructor
multi_ion_photo_xsections::~multi_ion_photo_xsections() {}



// ##################################################################
//
std::vector<double> multi_ion_photo_xsections::setup_photo_xsections(
    const std::vector<std::string> &photo_tracer_list  // photo tracer species
)
{
  spdlog::debug("Setting up multi-ion photoionisation cross-section module");

  // Resize photo cross-section table based on photo tracer list.
  photo_xsection_table.resize(photo_tracer_list.size());
  // Resizing bin fraction table
  bin_weight_table.resize(photo_tracer_list.size());
  // making replica
  photo_species_list = photo_tracer_list;

  // Get the threshold energy of every photo tracer species and the
  // energy bin size.
  for (int s = 0; s < photo_tracer_list.size(); s++) {
    photo_threshold.push_back(
        atomic_physics_data.get_ionisation_energy(photo_tracer_list[s]));
  }

  // returning ionisation threshold of photo tracer species
  return photo_threshold;
}


// ##################################################################
// Calculate cross-section of a species for a specific energy.
double multi_ion_photo_xsections::species_xsection(
    const std::string species,  // Name of the species
    const double &E             // Energy of photons
)
{
  if (E < atomic_physics_data.get_ionisation_energy(species)) return 0.0;

  // show error if the energy is larger than emax
  if (E > 1.2 * photo_xsection_parameters[species].fit_emax) {
    spdlog::warn(
        "multi-ion photo-xsection: energy exceeds fit range"
        " for {}, E: {} eV",
        species, E);
  }

  double x = E / photo_xsection_parameters[species].E0
             - photo_xsection_parameters[species].y0;
  double y = sqrt(x * x + pow(photo_xsection_parameters[species].y1, 2.0));
  double Fy =
      (pow(x - 1.0, 2.0) + pow(photo_xsection_parameters[species].yw, 2.0))
      * pow(y, 0.5 * photo_xsection_parameters[species].P - 5.5)
      * pow(
          1.0 + sqrt(y / photo_xsection_parameters[species].ya),
          -photo_xsection_parameters[species].P);
  double sigmaE = photo_xsection_parameters[species].sigma0 * 1e-18 * Fy;
  return sigmaE;
}



// ##################################################################
void multi_ion_photo_xsections::generate_ebin_mean_xsection(
    const std::vector<std::vector<double> > &ebins  // energy bins
)
{
  spdlog::info("Generating mean photo-xsection tables");

  std::vector<double> get_bin_mean_xsection_fraction = {0.0, 0.0};
  for (int s = 0; s < photo_species_list.size(); s++) {
    spdlog::debug(
        "multi-ion photo-xsection: generating"
        " mean-photo-xsection {:s}",
        photo_species_list[s]);
    for (int b = 0; b < ebins.size(); b++) {
      get_bin_mean_xsection_fraction = mean_bin_xsection(
          photo_species_list[s], photo_threshold[s], ebins[b]);
      photo_xsection_table[s].push_back(get_bin_mean_xsection_fraction[0]);
      bin_weight_table[s].push_back(get_bin_mean_xsection_fraction[1]);
    }
  }
}


// MEAN BIN CROSS-SECTION ###############################################
std::vector<double> multi_ion_photo_xsections::mean_bin_xsection(
    const std::string species,       // Name of the species
    const double &Eth,               //
    const std::vector<double> &Ebin  // Range of particular energy bin
)
{
  // This function returns a vector with two elements. The first element
  // represents the mean cross-section for the given energy bin, while
  // the last element indicates the proportion of the energy bin that
  // results in a non-zero cross-section.

  // If the ionization threshold exceeds the upper energy limit of the bin,
  // the function return a zero mean cross-section.
  // double set_precision  = 0.3; || Ebin[1] - Eth < set_precision
  if (Ebin[1] - Eth < 0.0) return {0.0, 0.0};

  // Calculate the bin_fraction for the given bin
  double bin_fraction = 0.0;
  double bin_width    = Ebin[1] - Ebin[0];
  // if the lower bin edge is larger than ionisation threshold of the species,
  // then return bin fraction as 1.0
  if ((Ebin[1] - Eth) / bin_width > 1) bin_fraction = 1.0;
  // if ionisation threshold of the species lies within the bin width, then
  // calculate the bin fraction
  else
    bin_fraction = (Ebin[1] - Eth) / bin_width;

  // if the bin fraction is less than 1.0, resetting the lower bin edge
  // to ionisation threshold of the species
  double lower_edge  = Ebin[0];
  double higher_edge = Ebin[1];
  if (bin_fraction < 1.0) lower_edge = Eth;

  // Calculating the mean cross-section for the given bin
  double mean_xsection = 0.0;
  double mean_energy   = 0.5 * (higher_edge + lower_edge);
  // for Composite Simpson's 3/8 rule, the number of sub-intervals must be
  // a multiple of 3. Here we choose this value as 99.
  int n = 99;
  // step size
  double h = bin_width / n;

  // Inverse of the cross-section for mean energy of the bin (inverse area)
  double IA = 1.0 / species_xsection(species, mean_energy);

  double sum = exp(-IA * species_xsection(species, lower_edge))
               + exp(-IA * species_xsection(species, higher_edge));

  for (int i = 1; i < n; i++) {
    double E = lower_edge + i * h;

    if (i % 3 == 0) {
      sum += 2 * exp(-IA * species_xsection(species, E));
    }
    else {
      sum += 3 * exp(-IA * species_xsection(species, E));
    }
  }
  mean_xsection = -1.0 * log(((3 * h / 8) * sum) / bin_width) / IA;

  // If the bin_fraction is less than 1.0, the calculation of mean
  // cross-section yields the value for the fraction of the bin where
  // the integration is carried out.
  return {mean_xsection, bin_fraction};
}



// ##################################################################
// Database for photo-ionisation cross-section
void multi_ion_photo_xsections::multi_ion_photo_xsections_database()
{

  spdlog::info("Setting up photo-ionisation cross-section database");

  // This database contain the analytic fit parameters for photoionisation
  // cross-sections of atoms and ions.
  // Ref: Verner D. A., Ferland G. J., Korista K. T., Yakovlev D.~G.,
  // 1996, ApJ, 465, 487. doi:10.1086/177435
  // Parameters for fit formulae are given below for each species.

  // HYDROGEN *******************************************
  photo_xsection_parameters["H"] =
      photo_xsection_datatype{5.000E+04, 4.298E-01, 5.475E+04, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["H1+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // HELIUM *******************************************
  photo_xsection_parameters["He"] =
      photo_xsection_datatype{5.000E+04, 1.361E+01, 9.492E+02, 1.469E+00,
                              3.188E+00, 2.039E+00, 4.434E-01, 2.136E+00};

  photo_xsection_parameters["He1+"] =
      photo_xsection_datatype{5.000E+04, 1.720E+00, 1.369E+04, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["He2+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // CARBON *******************************************
  photo_xsection_parameters["C"] =
      photo_xsection_datatype{2.910E+02, 2.144E+00, 5.027E+02, 6.216E+01,
                              5.101E+00, 9.157E-02, 1.133E+00, 1.607E+00};

  photo_xsection_parameters["C1+"] =
      photo_xsection_datatype{3.076E+02, 4.058E-01, 8.709E+00, 1.261E+02,
                              8.578E+00, 2.093E+00, 4.929E+01, 3.234E+00};

  photo_xsection_parameters["C2+"] =
      photo_xsection_datatype{3.289E+02, 4.614E+00, 1.539E+04, 1.737E+00,
                              1.593E+01, 5.922E+00, 4.378E-03, 2.528E-02};

  photo_xsection_parameters["C3+"] =
      photo_xsection_datatype{3.522E+02, 3.506E+00, 1.068E+02, 1.436E+01,
                              7.457E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["C4+"] =
      photo_xsection_datatype{5.000E+04, 4.624E+01, 2.344E+02, 2.183E+01,
                              2.581E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["C5+"] =
      photo_xsection_datatype{5.000E+04, 1.548E+01, 1.521E+03, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["C6+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // NITROGEN *******************************************
  photo_xsection_parameters["N"] =
      photo_xsection_datatype{4.048E+02, 4.034E+00, 8.235E+02, 8.033E+01,
                              3.928E+00, 9.097E-02, 8.598E-01, 2.325E+00};

  photo_xsection_parameters["N1+"] =
      photo_xsection_datatype{4.236E+02, 6.128E-02, 1.944E+00, 8.163E+02,
                              8.773E+00, 1.043E+01, 4.280E+02, 2.030E+01};

  photo_xsection_parameters["N2+"] =
      photo_xsection_datatype{4.473E+02, 2.420E-01, 9.375E-01, 2.788E+02,
                              9.156E+00, 1.850E+00, 1.877E+02, 3.999E+00};

  photo_xsection_parameters["N3+"] =
      photo_xsection_datatype{4.753E+02, 5.494E+00, 1.690E+04, 1.714E+00,
                              1.706E+01, 7.904E+00, 6.415E-03, 1.937E-02};

  photo_xsection_parameters["N4+"] =
      photo_xsection_datatype{5.043E+02, 4.471E+00, 8.376E+01, 3.297E+01,
                              6.003E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["N5+"] =
      photo_xsection_datatype{5.000E+04, 6.943E+01, 1.519E+02, 2.627E+01,
                              2.315E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["N6+"] =
      photo_xsection_datatype{5.000E+04, 2.108E+01, 1.117E+03, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["N7+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // OXYGEN *******************************************
  photo_xsection_parameters["O"] =
      photo_xsection_datatype{5.380E+02, 1.240E+00, 1.745E+03, 3.784E+00,
                              1.764E+01, 7.589E-02, 8.698E+00, 1.271E-01};

  photo_xsection_parameters["O1+"] =
      photo_xsection_datatype{5.581E+02, 1.386E+00, 5.967E+01, 3.175E+01,
                              8.943E+00, 1.934E-02, 2.131E+01, 1.503E-02};

  photo_xsection_parameters["O2+"] =
      photo_xsection_datatype{5.840E+02, 1.723E-01, 6.753E+02, 3.852E+02,
                              6.822E+00, 1.191E-01, 3.839E-03, 4.569E-01};

  photo_xsection_parameters["O3+"] =
      photo_xsection_datatype{6.144E+02, 2.044E-01, 8.659E-01, 4.931E+02,
                              8.785E+00, 3.143E+00, 3.328E+02, 4.285E+01};

  photo_xsection_parameters["O4+"] =
      photo_xsection_datatype{6.491E+02, 2.854E+00, 1.642E+04, 1.792E+00,
                              2.647E+01, 2.836E+01, 3.036E-02, 5.554E-02};

  photo_xsection_parameters["O5+"] =
      photo_xsection_datatype{6.837E+02, 7.824E+00, 6.864E+01, 3.210E+01,
                              5.495E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["O6+"] =
      photo_xsection_datatype{5.000E+04, 8.709E+01, 1.329E+02, 2.535E+01,
                              2.336E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["O7+"] =
      photo_xsection_datatype{5.000E+04, 2.754E+01, 8.554E+02, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["O8+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // NEON *******************************************
  photo_xsection_parameters["Ne"] =
      photo_xsection_datatype{8.701E+02, 4.870E+00, 4.287E+03, 5.798E+00,
                              8.355E+00, 2.434E-01, 4.236E-02, 5.873E+00};

  photo_xsection_parameters["Ne1+"] =
      photo_xsection_datatype{8.831E+02, 1.247E+01, 1.583E+03, 3.935E+00,
                              7.810E+00, 6.558E-02, 1.520E+00, 1.084E-01};

  photo_xsection_parameters["Ne2+"] =
      photo_xsection_datatype{9.131E+02, 7.753E-01, 5.708E+00, 6.725E+01,
                              1.005E+01, 4.633E-01, 7.654E+01, 2.023E+00};

  photo_xsection_parameters["Ne3+"] =
      photo_xsection_datatype{9.480E+02, 5.566E+00, 1.685E+03, 6.409E+02,
                              3.056E+00, 8.290E-03, 5.149E+00, 6.687E+00};

  photo_xsection_parameters["Ne4+"] =
      photo_xsection_datatype{9.873E+02, 1.248E+00, 2.430E+00, 1.066E+02,
                              8.999E+00, 6.855E-01, 9.169E+01, 3.702E-01};

  photo_xsection_parameters["Ne5+"] =
      photo_xsection_datatype{1.031E+03, 1.499E+00, 9.854E-01, 1.350E+02,
                              8.836E+00, 1.656E+00, 1.042E+02, 1.435E+00};

  photo_xsection_parameters["Ne6+"] =
      photo_xsection_datatype{1.078E+03, 4.888E+00, 1.198E+04, 1.788E+00,
                              2.550E+01, 2.811E+01, 2.536E-02, 4.417E-02};

  photo_xsection_parameters["Ne7+"] =
      photo_xsection_datatype{1.125E+03, 1.003E+01, 5.631E+01, 3.628E+01,
                              5.585E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Ne8+"] =
      photo_xsection_datatype{5.000E+04, 1.586E+02, 6.695E+01, 3.352E+01,
                              2.002E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Ne9+"] =
      photo_xsection_datatype{5.000E+04, 4.304E+01, 5.475E+02, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Ne10+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // SILICON *******************************************
  photo_xsection_parameters["Si"] =
      photo_xsection_datatype{1.060E+02, 2.317E+01, 2.506E+01, 2.057E+01,
                              3.546E+00, 2.837E-01, 1.672E-05, 4.207E-01};

  photo_xsection_parameters["Si1+"] =
      photo_xsection_datatype{1.186E+02, 2.556E+00, 4.140E+00, 1.337E+01,
                              1.191E+01, 1.570E+00, 6.634E+00, 1.272E-01};

  photo_xsection_parameters["Si2+"] =
      photo_xsection_datatype{1.311E+02, 1.659E-01, 5.790E-04, 1.474E+02,
                              1.336E+01, 8.626E-01, 9.613E+01, 6.442E-01};

  photo_xsection_parameters["Si3+"] =
      photo_xsection_datatype{1.466E+02, 1.288E+01, 6.083E+00, 1.356E+06,
                              3.353E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Si4+"] =
      photo_xsection_datatype{1.887E+03, 7.761E-01, 8.863E-01, 1.541E+02,
                              9.980E+00, 1.303E+00, 2.009E+02, 4.537E+00};

  photo_xsection_parameters["Si5+"] =
      photo_xsection_datatype{1.946E+03, 6.305E+01, 7.293E+01, 1.558E+02,
                              2.400E+00, 2.989E-03, 1.115E+00, 8.051E-02};

  photo_xsection_parameters["Si6+"] =
      photo_xsection_datatype{2.001E+03, 3.277E-01, 6.680E-02, 4.132E+01,
                              1.606E+01, 3.280E+00, 1.149E-02, 6.396E-01};

  photo_xsection_parameters["Si7+"] =
      photo_xsection_datatype{2.058E+03, 7.655E-01, 3.477E-01, 3.733E+02,
                              8.986E+00, 1.476E-03, 3.850E+02, 8.999E-02};

  photo_xsection_parameters["Si8+"] =
      photo_xsection_datatype{2.125E+03, 3.343E-01, 1.465E-01, 1.404E+03,
                              8.503E+00, 1.646E+00, 1.036E+03, 2.936E-01};

  photo_xsection_parameters["Si9+"] =
      photo_xsection_datatype{2.194E+03, 8.787E-01, 1.950E-01, 7.461E+02,
                              8.302E+00, 4.489E-01, 4.528E+02, 1.015E+00};

  photo_xsection_parameters["Si10+"] =
      photo_xsection_datatype{2.268E+03, 1.205E+01, 1.992E+04, 1.582E+00,
                              2.425E+01, 2.392E+01, 1.990E-02, 1.007E-02};

  photo_xsection_parameters["Si11+"] =
      photo_xsection_datatype{2.336E+03, 3.560E+01, 2.539E+01, 3.307E+01,
                              4.728E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Si12+"] =
      photo_xsection_datatype{5.000E+04, 2.752E+02, 4.754E+01, 2.848E+01,
                              2.135E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Si13+"] =
      photo_xsection_datatype{5.000E+04, 8.447E+01, 2.793E+02, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Si14+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  // SULFUR *******************************************
  photo_xsection_parameters["S"] =
      photo_xsection_datatype{1.700E+02, 1.808E+01, 4.564E+04, 1.000E+00,
                              1.361E+01, 6.385E-01, 9.935E-01, 2.486E-01};

  photo_xsection_parameters["S1+"] =
      photo_xsection_datatype{1.846E+02, 8.787E+00, 3.136E+02, 3.442E+00,
                              1.281E+01, 7.354E-01, 2.782E+00, 1.788E-01};

  photo_xsection_parameters["S2+"] =
      photo_xsection_datatype{1.995E+02, 2.027E+00, 6.666E+00, 5.454E+01,
                              8.611E+00, 4.109E+00, 1.568E+01, 9.421E+00};

  photo_xsection_parameters["S3+"] =
      photo_xsection_datatype{2.164E+02, 2.173E+00, 2.606E+00, 6.641E+01,
                              8.655E+00, 1.863E+00, 1.975E+01, 3.361E+00};

  photo_xsection_parameters["S4+"] =
      photo_xsection_datatype{2.350E+02, 1.713E-01, 5.072E-04, 1.986E+02,
                              1.307E+01, 7.880E-01, 9.424E+01, 6.265E-01};

  photo_xsection_parameters["S5+"] =
      photo_xsection_datatype{2.557E+02, 1.413E+01, 9.139E+00, 1.656E+03,
                              3.626E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["S6+"] =
      photo_xsection_datatype{2.569E+03, 3.757E-01, 5.703E-01, 1.460E+02,
                              1.135E+01, 1.503E+00, 2.222E+02, 4.606E+00};

  photo_xsection_parameters["S7+"] =
      photo_xsection_datatype{2.641E+03, 1.462E+01, 3.161E+01, 1.611E+01,
                              8.642E+00, 1.153E-03, 1.869E+01, 3.037E-01};

  photo_xsection_parameters["S8+"] =
      photo_xsection_datatype{2.705E+03, 1.526E-01, 9.646E+03, 1.438E+03,
                              5.977E+00, 1.492E+00, 1.615E-03, 4.049E-01};

  photo_xsection_parameters["S9+"] =
      photo_xsection_datatype{2.782E+03, 1.040E+01, 5.364E+01, 3.641E+01,
                              7.090E+00, 2.310E+00, 1.775E+01, 1.663E+00};

  photo_xsection_parameters["S10+"] =
      photo_xsection_datatype{2.859E+03, 6.485E+00, 1.275E+01, 6.583E+01,
                              7.692E+00, 1.678E+00, 3.426E+01, 1.370E-01};

  photo_xsection_parameters["S11+"] =
      photo_xsection_datatype{2.941E+03, 2.443E+00, 3.490E-01, 5.411E+02,
                              7.769E+00, 7.033E-01, 2.279E+02, 1.172E+00};

  photo_xsection_parameters["S12+"] =
      photo_xsection_datatype{3.029E+03, 1.474E+01, 2.294E+04, 1.529E+00,
                              2.568E+01, 2.738E+01, 2.203E-02, 1.073E-02};

  photo_xsection_parameters["S13+"] =
      photo_xsection_datatype{3.107E+03, 3.310E+01, 2.555E+01, 3.821E+01,
                              5.037E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["S14+"] =
      photo_xsection_datatype{5.000E+04, 4.390E+02, 2.453E+01, 4.405E+01,
                              1.765E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["S15+"] =
      photo_xsection_datatype{5.000E+04, 1.104E+02, 2.139E+02, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["S16+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};


  // IRON *******************************************
  photo_xsection_parameters["Fe"] =
      photo_xsection_datatype{6.600E+01, 5.461E-02, 3.062E-01, 2.671E+07,
                              7.923E+00, 2.069E+01, 1.382E+02, 2.481E-01};

  photo_xsection_parameters["Fe1+"] =
      photo_xsection_datatype{7.617E+01, 1.761E-01, 4.365E+03, 6.298E+03,
                              5.204E+00, 1.141E+01, 9.272E+01, 1.075E+02};

  photo_xsection_parameters["Fe2+"] =
      photo_xsection_datatype{8.705E+01, 1.698E-01, 6.107E+00, 1.555E+03,
                              8.055E+00, 8.698E+00, 1.760E+02, 1.847E+01};

  photo_xsection_parameters["Fe3+"] =
      photo_xsection_datatype{1.067E+02, 2.544E+01, 3.653E+02, 8.913E+00,
                              6.538E+00, 5.602E-01, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Fe4+"] =
      photo_xsection_datatype{1.288E+02, 7.256E-01, 1.523E-03, 3.736E+01,
                              1.767E+01, 5.064E+01, 8.871E+01, 5.280E-02};

  photo_xsection_parameters["Fe5+"] =
      photo_xsection_datatype{1.527E+02, 2.656E+00, 5.259E-01, 1.450E+01,
                              1.632E+01, 1.558E+01, 3.361E+01, 3.743E-03};

  photo_xsection_parameters["Fe6+"] =
      photo_xsection_datatype{1.783E+02, 5.059E+00, 2.420E+04, 4.850E+04,
                              2.374E+00, 2.516E-03, 4.546E-01, 2.683E+01};

  photo_xsection_parameters["Fe7+"] =
      photo_xsection_datatype{2.055E+02, 7.098E-02, 1.979E+01, 1.745E+04,
                              6.750E+00, 2.158E+02, 2.542E+03, 4.672E+02};

  photo_xsection_parameters["Fe8+"] =
      photo_xsection_datatype{9.211E+02, 6.741E+00, 2.687E+01, 1.807E+02,
                              6.290E+00, 2.387E-04, 2.494E+01, 8.251E+00};

  photo_xsection_parameters["Fe9+"] =
      photo_xsection_datatype{9.590E+02, 6.886E+01, 6.470E+01, 2.062E+01,
                              4.111E+00, 2.778E-04, 1.190E-05, 6.570E-03};

  photo_xsection_parameters["Fe10+"] =
      photo_xsection_datatype{9.983E+02, 8.284E+00, 3.281E+00, 5.360E+01,
                              8.571E+00, 3.279E-01, 2.971E+01, 5.220E-01};

  photo_xsection_parameters["Fe11+"] =
      photo_xsection_datatype{1.039E+03, 6.295E+00, 1.738E+00, 1.130E+02,
                              8.037E+00, 3.096E-01, 4.671E+01, 1.425E-01};

  photo_xsection_parameters["Fe12+"] =
      photo_xsection_datatype{1.081E+03, 1.317E-01, 2.791E-03, 2.487E+03,
                              9.791E+00, 6.938E-01, 2.170E+03, 6.852E-03};

  photo_xsection_parameters["Fe13+"] =
      photo_xsection_datatype{1.125E+03, 8.509E-01, 1.454E-01, 1.239E+03,
                              8.066E+00, 4.937E-01, 4.505E+02, 2.504E+00};

  photo_xsection_parameters["Fe14+"] =
      photo_xsection_datatype{1.181E+03, 5.555E-02, 2.108E+02, 2.045E+04,
                              6.033E+00, 1.885E-03, 2.706E-04, 1.628E+00};

  photo_xsection_parameters["Fe15+"] =
      photo_xsection_datatype{1.216E+03, 2.873E+01, 1.207E+01, 5.150E+02,
                              3.846E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Fe16+"] =
      photo_xsection_datatype{7.651E+03, 3.444E-01, 1.452E+00, 3.960E+02,
                              1.013E+01, 1.264E+00, 2.891E+01, 3.404E+00};

  photo_xsection_parameters["Fe17+"] =
      photo_xsection_datatype{7.769E+03, 3.190E+01, 2.388E+00, 2.186E+01,
                              9.589E+00, 2.902E-02, 3.805E+01, 4.805E-01};

  photo_xsection_parameters["Fe18+"] =
      photo_xsection_datatype{7.918E+03, 7.519E-04, 6.066E-05, 1.606E+06,
                              8.813E+00, 4.398E+00, 1.915E+06, 3.140E+01};

  photo_xsection_parameters["Fe19+"] =
      photo_xsection_datatype{8.041E+03, 2.011E+01, 4.455E-01, 4.236E+01,
                              9.724E+00, 2.757E+00, 6.847E+01, 3.989E+00};

  photo_xsection_parameters["Fe20+"] =
      photo_xsection_datatype{8.184E+03, 9.243E+00, 1.098E+01, 7.637E+01,
                              7.962E+00, 1.748E+00, 4.446E+01, 3.512E+00};

  photo_xsection_parameters["Fe21+"] =
      photo_xsection_datatype{8.350E+03, 9.713E+00, 7.204E-02, 1.853E+02,
                              8.843E+00, 9.551E-03, 1.702E+02, 4.263E+00};

  photo_xsection_parameters["Fe22+"] =
      photo_xsection_datatype{8.484E+03, 4.575E+01, 2.580E+04, 1.358E+00,
                              2.604E+01, 2.723E+01, 3.582E-02, 8.712E-03};

  photo_xsection_parameters["Fe23+"] =
      photo_xsection_datatype{8.638E+03, 7.326E+01, 1.276E+01, 4.914E+01,
                              4.941E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Fe24+"] =
      photo_xsection_datatype{5.000E+04, 1.057E+03, 1.195E+01, 5.769E+01,
                              1.718E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Fe25+"] =
      photo_xsection_datatype{5.000E+04, 2.932E+02, 8.099E+01, 3.288E+01,
                              2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00};

  photo_xsection_parameters["Fe26+"] =
      photo_xsection_datatype{0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
                              0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00};
}
// END OF COLLISIONAL IONISATION DATABASE ***********************************
