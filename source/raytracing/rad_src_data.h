/// \file rad_src_data.h
/// \author: Jonathan Mackey
///
/// Contains structures used for keeping track of radiation sources.
///
/// Modifications:
/// - 2018.03.19 JM: Moved from various files to get everything in
///    one place.
///

#ifndef RAD_SRC_DATA_H
#define RAD_SRC_DATA_H

#include <array>
#include <defines/functionality_flags.h>
#include <defines/testing_flags.h>
#include <sim_constants.h>
#include <string>
#include <vector>

// ##################################################################
// ##################################################################

///
/// Radiation source struct.
///
struct rad_src_info {
  std::array<double, MAX_DIM> pos;  ///< src position (physical units).
  double strength;  ///< src strength (photons/sec, or ergs/sec for multifreq.)
  double Rstar;     ///< stellar radius in solar radii (for multifreq.
                    ///< photoionisation).
  double Tstar;     ///< stellar effective temperature (for multifreq.
                    ///< photoionisation).
  int id;           ///< src identifier
  int type;         ///< src type: either RT_SRC_DIFFUSE or RT_SRC_SINGLE.
  int update;       ///< how the source is updated: RT_UPDATE_IMPLICIT=1,
                    ///< RT_UPDATE_EXPLICIT=2
  int at_infinity;  ///< set to true if source is at infinity.
  int ongrid;       ///< true if source is within the simulation domain.
  ///
  /// "effect" is what the source does, and this defines many of its
  /// properties implicitly.  Options are:
  /// - RT_EFFECT_UV_HEATING,
  /// - RT_EFFECT_PION_MONO,
  /// - RT_EFFECT_MFION,
  /// - RT_EFFECT_PHOTODISS, (UNUSED)
  /// - RT_EFFECT_PION_EQM, (UNUSED--for photoionisation equilibrium)
  /// - RT_EFFECT_HHE_MFQ  (NOT YET IMPLEMENTED)
  ///
  int effect;
  ///
  /// "NTau" sets the number of quantities traced from the source.
  ///
  int NTau;
  ///
  /// What provides the opacity: RT_OPACITY_TOTAL, RT_OPACITY_MINUS,
  /// RT_OPACITY_TRACER.
  int opacity_src;
  ///
  /// optional tracer variable index in state vector, for opacity
  /// calculation.
  int opacity_var;
  ///
  /// Optional text file with output from stellar evolution code for
  /// time-varying source.
  std::string EvoFile;
};

// ##################################################################
// ##################################################################

///
/// List of all radiation sources.  Used in SimParams class.
///
struct rad_sources {
  int Nsources;
  std::vector<struct rad_src_info> sources;
};

// ##################################################################
// ##################################################################

#define MAX_TAU 20
///
/// Radiation Source data struct, used for passing info to microphysics classes.
///
struct rt_source_data {
  double
      Vshell;  ///< Shell volume for discrete photo-ionisation/-heating rates.
  double dS;   ///< Path length through cell.
  double strength[MAX_TAU];  ///< Luminosity (or flux if source at infinity).
  double Column[MAX_TAU];    ///< integral of quantities along LOS to near edge
                             ///< of cell.
  double DelCol[MAX_TAU];    ///< integral of quantities along LOS through cell.
  int id;                    ///< source id.
  int type;                  ///< diffuse-radiation or a real source.
  short unsigned int NTau;  ///< Number of LOS quantities traced for the source.
};

// ##################################################################
// ##################################################################

///
/// Struct to hold info on radiation sources.  An extension of
/// rad_src_info
///
struct rad_source {
  ///
  /// pointer to source (set to SimPM.RS.source[id]) with the basic
  /// info about the source.
  ///
  struct rad_src_info *s;

  class cell *sc;    ///< nearest cell to source.
  bool src_on_grid;  ///< true if source is at a grid cell.
  std::array<int, MAX_DIM>
      ipos;  ///< source position in integer form (grid units, dx=2).

  ///
  /// This struct is used by the code to pass cell and source data
  /// to the microphysics integrator.  It contains the relevant source
  /// information, and also some cell-source geometry information which
  /// must be set on a cell-by-cell basis as rays are traced.
  /// The struct is declared in source/microphysics/microphysics_base.h
  ///
  struct rt_source_data data;
};

// ##################################################################
// ##################################################################

#endif  // RAD_SRC_DATA_H
