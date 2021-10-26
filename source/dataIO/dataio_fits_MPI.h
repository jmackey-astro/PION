/// \file dataio_fits.h
/// \author Jonathan Mackey
///
/// This file contains class declarations for the serial and parallel
/// versions of the Fits data I/O class.  Info on FITS can be found at
/// http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html at
/// NASA.  It's advantage is that it is simple, but I am using very
/// basic parallel I/O, where each processor writes its own file.
///
/// modified:\n
/// - 2015.01.28 JM: wrote file, for new code structure.
///
#ifndef DATAIO_FITS_MPI_H
#define DATAIO_FITS_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef FITS

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_fits.h"
#include "fitsio.h"
#include "sub_domain/sub_domain.h"

///
/// Class for reading and writing uniform grid data to fits images.
/// This deals with parallel I/O, and is not very well developed or
/// efficient.  For production simulations it is much better to use
/// the SILO format.
///
class DataIOFits_pllel : public DataIOFits {
public:
  ///
  /// Constructor.
  ///
  DataIOFits_pllel(
      class SimParams &,  ///< pointer to simulation parameters
      class Sub_domain *  ///< address of sub_domain controller class.
  );

  /// Destructor (doensn't have much to do).
  ~DataIOFits_pllel();

  ///
  /// Choose filename based on counter and base-name.
  ///
  std::string choose_filename(
      const std::string,  ///< filebase passed in from main code.
      const int           ///< file counter to use (e.g. timestep).
  );

  ///
  /// This writes the fits header for the simulation parameters,
  /// and then the data, with a separate image for each variable.
  ///
  /// If the solver pointer is not null, it also writes an image for the
  /// internal energy (and the Magnetic Field divergence and total Pressure
  /// if the solver is an MHD solver).
  ///
  int OutputData(
      const string,                     ///< File to write to
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &,                ///< pointer to simulation parameters
      const long int  ///< number to stamp file with (e.g. timestep)
  );

  ///
  /// This writes the header and data for the simulation parameters, for
  /// a given level in the nested-grid structure.
  ///
  /// If the solver pointer is not null, it also writes some derived
  /// variables such as Temperature, Div(B), etc.
  ///
  int SaveLevelData(
      const string,           ///< File-base to write to
      const int,              ///< level in nested grid to write.
      class GridBaseClass *,  ///< grid pointer.
      class SimParams &,      ///< simulation parameters
      const long int          ///< timestep
  );

  ///
  /// This reads the fits images in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the fits header, which should be read first.
  ///
  int ReadData(
      string,                           ///< file to read from
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &                 ///< pointer to simulation parameters
  );

protected:
  class Sub_domain *mpiPM;
};

#endif  // if FITS

#endif  // DATAIO_FITS_MPI_H
