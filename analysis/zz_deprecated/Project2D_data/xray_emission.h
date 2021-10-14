/// \file xray-emission.h
/// \author Jonathan Mackey
/// \date 2016.07.04
///
/// Purpose:
/// - setup table for calculating X-ray emission.
///


#ifndef XRAY_EMISSION_H
#define XRAY_EMISSION_H

// ##################################################################
// ##################################################################


class Xray_emission {
public:
  Xray_emission();

  ~Xray_emission() { free_xray_tables_priv(&LT, &L1, &L2, &L3, &L4); }

  ///
  /// Function to calculate the X-ray emissivity given an input gas
  /// temperature.  Returns the emissivity in units of ergs.cm^3/sec,
  /// and should be multiplied by n_e*n_H.
  ///
  void get_xray_emissivity(
      const double,  ///< Temperature (K)
      double *       ///< Results.
  );


private:
  size_t XNel;  ///< Number of rows in table.

  double *LT,  ///< log temperature
      *L1,     ///< log emissivity in band 1
      *L2,     ///< log emissivity in band 2
      *L3,     ///< log emissivity in band 3
      *L4;     ///< log emissivity in band 4

  ///
  /// Function to set up the X-ray emissivity tables, from XSPEC
  ///
  int setup_xray_tables_priv(
      size_t *,   ///< Number of elements read
      double **,  ///< Log(T) array
      double **,  ///< Log(L>0.1keV) array
      double **,  ///< Log(L>0.5keV) array
      double **,  ///< Log(L>1.0keV) array
      double **   ///< Log(L>5.0keV) array
  );

  ///
  /// Function to free memory allocated for the X-ray emissivity tables
  ///
  void free_xray_tables_priv(
      double **,  ///< Log(T) array
      double **,  ///< Log(L>0.1keV) array
      double **,  ///< Log(L>0.5keV) array
      double **,  ///< Log(L>1.0keV) array
      double **   ///< Log(L>5.0keV) array
  );
};

// ##################################################################
// ##################################################################

#endif  // XRAY_EMISSION_H
