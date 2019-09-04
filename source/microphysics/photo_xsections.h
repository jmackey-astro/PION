///
/// \file photo_xsections.h
/// \author Jonathan Mackey
/// \date 2015.05.19
///
/// This class calculates ion heating and ionization rates for a gas
/// with a mixture of H, H2, and He, where n(He)/(n(H)+2n(H2))=0.1 is
/// hardcoded.
///
/// References:
/// - DYL99: Dalgarno, Yan, Liu (1999,ApJS,125,237) Table 7
/// - MS05:  Meijerink & Spaans (2005,A&A,436,397)  App. B
/// - PCP12: Panoglou,Cabrit,+  (2012,A&A,538,A2)
/// - Yan (1997) PhD thesis:
///     http://adsabs.harvard.edu/abs/1997PhDT.........4Y
/// - Meijerink (2006) PhD thesis:
///     
///
/// Modifications:
/// -2015.08.18 JM: removed NH from parameter list (not needed).
/// -2015.10.08 JM: Added absorbed power per H-atom as a parameter in
///  ionization_heating_rates().  Renamed some variables.
/// - 2015.12.01 JM: started to implement user-defined energy bins.
/// - 2015.12.02 JM: updated to use user-specified energy bins.

//#define TEST_XDR
//#define NEON
#include <iostream>
#include <string>
#include <stdio.h>
#include <map>
#include <stdlib.h>
using namespace std;



/// Struct for storing parameter data:
  struct params_struct {
  int Z, N; // atomic number, number of electrons
  double E_th; 
  double E_max;
  double E_0;
  double sigma_0;
  double y_a;
  double P;
  double y_w;
  double y_0;
  double y_1;
  double *xsections;
  };



///
/// x-ray microphysics class.
/// This class uses keV as energy/frequency units, and CGS units for
/// everything else.  All function arguments are assumed to be in
/// these units.
///
class photo_xsections
{
  public:
  ///
  /// constructor: inputs are Emin [], Emax [], nbins.
  /// Emin and Emax are arrays of length nbins, with the lower
  /// and upper energy limits of each bin.
  /// Emin must be >= 0.03 keV
  /// Emax must be <= 50   keV  ???
  /// 
  photo_xsections(
    /*const double *, ///< Emin
    const double *, ///< Emax
    const int       ///< nbins*/
    );

  ///
  /// Destructor: cleans up memory.
  ///
  ~photo_xsections();
  
  ///
  /// returns the number of energy bins.
  ///
  unsigned int get_nbins()
  {return Nbin;}

  ///
  /// returns an array of energy-bin central energies (keV)
  ///
  void get_energies(double *E)
  {for (unsigned int v=0;v<Nbin;v++) E[v] = E0[v]; return;}

  ///
  /// returns an array of cross sections for each energy (cm^2).
  ///
  void get_xsections(double *x)
  {for (unsigned int v=0;v<Nbin;v++) x[v] = xs[v]; return;}

  ///
  /// Calculate cross-section of the gas, per H nucleus, as a 
  /// function of energy (cm^2)
  ///
  double ion_xsection(
    const double,  ///< energy of photons (keV).
    string ///< name of ion
    );

  ///
  /// Calculates the mean cross-section of an energy bin between
  /// two limits, using a mean attenuation.
  ///
  double set_bin_mean_xsection(
    const double, ///< Minimum energy in bin (keV)
    const double, ///< Maximum energy in bin (keV)
    string       ///< ion name
    );

  protected:
  unsigned int Nbin;  ///< number of energy bins.
  double *EBmin;      ///< list of min-energies of energy bins (keV).
  double *EBmax;      ///< list of max-energies of energy bins (keV).
  double *dE;         ///< bin size (keV)
  double *E0;         ///< array of central energies of bins (keV).
  double *xs;         ///< array of cross-sections for each bin (cm^2).
  double Emin;        ///< min energy of the spectrum (keV)
  double Emax;        ///< max energy of the spectrum (keV)
  
  /// Calculate photoionization cross section of Neon
  /// Return value in cm^2
  ///
  double get_Ne0_xsection(
    const double   ///< Energy of photon (keV)
    );
  
  ///
  /// Calculates the mean cross-section of an energy bin between
  /// two limits, using a mean attenuation, for Neon.
  ///
  double set_Ne_bin_mean_xsection(
    const double, ///< Minimum energy in bin (keV)
    const double  ///< Maximum energy in bin (keV)
    );

  public:
  double get_Ne_primary_photoion_rate(
    double *   ///< input fluxes for each bin (already attenuated)
    );

  //78
  /// list of ions in the order they appear in VF96.
  /// Note this is photoionisation rate INTOOO these ions.
  string *ions;
  int n_ions; ///< number of ions in list.
  /// Map for parameters by ion name
  std::map<string, params_struct> parameters;

};



