

#ifndef GET_SIM_INFO_H
#define GET_SIM_INFO_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "dataIO/readparams.h"
#include <iostream>

class get_sim_info {
public:
  get_sim_info();
  ~get_sim_info();
  int read_gridparams(
      std::string,       ///< parameter file.
      class SimParams &  ///< pointer to simulation paramters.
  );

private:
  /** \brief Reads in extra physics params from the text file, if they are
   * there. */
  int read_extra_physics(
      class SimParams &  ///< pointer to simulation paramters.
  );
  /** \brief if we are doing raytracing, read in the source list. */
  int read_radsources(class SimParams &  ///< pointer to simulation paramters.
  );
  /** \brief Reads in units params from the text file, if they are there. */
  int read_units();
  ///
  /// If we got one or more stellar wind sources, this function reads
  /// their properties
  ///
  int read_wind_sources(class SimParams &  ///< pointer to simulation paramters.
  );

  ///
  /// Read in parameters for a stellar jet, if requested by N_jet !=0
  ///
  int read_jet_params(
      class SimParams &spar,  ///< pointer to simulation paramters.
      class JetParams &jpar   ///< pointer to jet parameters class.
  );

  class ReadParams *rp;
};

#endif  // GET_SIM_INFO_H
