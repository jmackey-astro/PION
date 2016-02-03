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

#include "dataIO/dataio.h"
#include "dataIO/dataio_fits.h"
#include "fitsio.h"
#include "MCMD_control.h"


///
/// Class for reading and writing uniform grid data to fits images.
/// This deals with parallel I/O exclusively.
///
class DataIOFits_pllel : public DataIOFits {
  public:
  ///
  /// Constructor.
  ///
  DataIOFits_pllel(
      class MCMDcontrol *  ///< address of MCMD controller class.
  );

  /// Destructor (doensn't have much to do).
  ~DataIOFits_pllel();

  ///
  /// Choose filename based on counter and base-name.
  ///
  std::string choose_filename(
        const std::string, ///< filebase passed in from main code.
        const int          ///< file counter to use (e.g. timestep).
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
        const string, ///< File to write to
        class GridBaseClass *, ///< pointer to data.
        const long int ///< number to stamp file with (e.g. timestep)
        );

  ///
  /// This reads the fits images in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the fits header, which should be read first.
  ///
  int ReadData(
        string, ///< file to read from
        class GridBaseClass * ///< pointer to data.
        );

  protected:
  class MCMDcontrol *mpiPM;

};

#endif // if FITS

#endif // DATAIO_FITS_MPI_H
