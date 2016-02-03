

#ifndef GET_SIM_INFO_H
#define GET_SIM_INFO_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include <iostream>
#include "dataIO/readparams.h"


class get_sim_info {
public:
  get_sim_info();
  ~get_sim_info();
  int read_gridparams(std::string ///< paramfile.
		      );
private:
  /** \brief Reads in extra physics params from the text file, if they are there. */
  int read_extra_physics();
  /** \brief if we are doing raytracing, read in the source list. */
  int read_radsources();
  /** \brief Reads in units params from the text file, if they are there. */
  int read_units();
  ///
  /// If we got one or more stellar wind sources, this function reads
  /// their properties
  ///
  int read_wind_sources();

  class ReadParams *rp;
};

#endif // GET_SIM_INFO_H

