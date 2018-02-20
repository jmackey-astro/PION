/// \file xray-emission.h
/// \author Jonathan Mackey
/// \date 2016.07.04
///
/// Purpose:
/// - setup table for calculating X-ray emission.
///
/// Modifications:
///
/// - 2018-01-11 SG: Added in more energy bands (0.2, 2, 10 keV)


#ifndef XRAY_EMISSION_H
#define XRAY_EMISSION_H

// ##################################################################
// ##################################################################


class Xray_emission {
  public:
  Xray_emission();

  ~Xray_emission()
  {
    free_xray_tables_priv(&LT, &L1, &L2, &L3, &L4, &L5, &L6, &L7);
  }

  ///
  /// Function to calculate the X-ray emissivity given an input gas
  /// temperature.  Returns the emissivity in units of ergs.cm^3/sec,
  /// and should be multiplied by n_e*n_H.
  ///
  void get_xray_emissivity(
        const double, ///< Temperature (K)
        double *  ///< Results.
        );


  private:
  size_t XNel;  ///< Number of rows in table.

  double 
    *LT,  ///< log temperature
    *L1,  ///< log emissivity in band 1
    *L2,  ///< log emissivity in band 2
    *L3,  ///< log emissivity in band 3
    *L4,  ///< log emissivity in band 4
    *L5,  ///< log emissivity in band 5
    *L6,  ///< log emissivity in band 6
    *L7;  ///< log emissivity in band 7

  ///
  /// Function to set up the X-ray emissivity tables, from XSPEC
  ///
  int setup_xray_tables_priv(
      size_t *, ///< Number of elements read
      double **, ///< Log(T) array
      double **, ///< Log(L>0.1keV) array
      double **, ///< Log(L>0.2keV) array
      double **, ///< Log(L>0.5keV) array
      double **, ///< Log(L>1.0keV) array
      double **, ///< Log(L>2.0keV) array
      double **,  ///< Log(L>5.0keV) array
      double ** ///< Log(L>10.0keV) array
      );

  ///
  /// Function to free memory allocated for the X-ray emissivity tables
  ///
  void free_xray_tables_priv(
      double **, ///< Log(T) array
      double **, ///< Log(L>0.1keV) array
      double **, ///< Log(L>0.2keV) array
      double **, ///< Log(L>0.5keV) array
      double **, ///< Log(L>1.0keV) array
      double **, ///< Log(L>2.0keV) array
      double **,  ///< Log(L>5.0keV) array
      double ** ///< Log(L>10.0keV) array
      );

};

// ##################################################################
// ##################################################################

#endif // XRAY_EMISSION_H

