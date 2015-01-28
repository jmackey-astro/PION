/// \file dataio_fits.h
/// \author Jonathan Mackey
///
/// This file contains class declarations for the serial and parallel
/// versions of the Fits data I/O class.  Info on FITS can be found at
/// http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html at 
/// NASA.
/// This file is forked off from dataio.h
///
/// modified:\n 
/// - 2010-02-03 JM: made ~utility_fitsio() virtual function.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010-07-20/21 JM: Work on new dataio structure with a list of
///    parameters to read and write.  Made a start.
/// - 2010.07.23 JM: removed obselete read_fits_header(),
///    write_fits_header() functions.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2011.06.02 JM: Added WriteHeader() function so I can over-write header
///    parameters and restart a simulation with e.g. different microphysics.
/// - 2013.04.16 JM: choose_filename() made public.
/// - 2015.01.28 JM: Tidied up, added "virtual" to some functions.

#ifndef DATAIO_FITS_H
#define DATAIO_FITS_H

#ifdef FITS

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "dataio.h"
#include "fitsio.h"

///
/// Utility class for various fits operations.
///
class utility_fitsio {
 public:
  utility_fitsio();
  virtual ~utility_fitsio() {}
  
  ///
  /// Called by OutputData() to create an empty fits image in the 
  /// current output file, to hold the simulation data for some variable.
  ///
  int create_fits_image(
        fitsfile *,   ///< fitsfile pointer.
        const string, ///< extname of hdu to create.
        const int,    ///< Dimensionality of image
        const int *   ///< vector npix in each dimension
        );
  
  ///
  /// Called by OutputData() once an image has been created; this
  /// file writes data to a subset of an image (can be the full image if local
  /// and global xmin are the same, and ntot is the total number of pixels).
  /// Note the subset has to be a contiguous region, brick-shaped.
  /// 
  /// For parallel I/O, if we are writing multiple files, the image is just
  /// the size of the local grid.  There is an option to write a single file,
  /// but this has synchronization problems, so should not be used without
  /// modification.  In this case the image size is the whole simulation, but
  /// this function only writes the part of the image relating to the local
  /// domain.
  ///
  int write_fits_image(
        fitsfile *,     ///< fitsfile pointer.
        const string,   ///< extname of hdu to write to.
        const double *, ///< vector local xmin (subdomain).
        const double *, ///< vector global xmin (full domain).
        const double,   ///< pixel size.   
        const int,      ///< dimensionality of image.
        const int *,    ///< vector npix to write in each direction.
        long int,       ///< total npix to write (size of data array)
        double *  ///< array of data elements to write to subset of image
        );

  ///
  /// This checks that the fits image to be read or written has the
  /// dimensions and number of pixels we expect.  If not, it says why and the
  /// code exits, as this would be a serious problem.
  ///
  int check_fits_image_dimensions(
        fitsfile *,   ///< fitsfile pointer.
        const string, ///< extname of hdu to check.
        const int,    ///< dimensionality we are expecting.
        const int *   ///< npix in each direction to compare image to.
        );

  ///
  /// Given that we have opened a fits file, read in the header, 
  /// checked the image is the right size, this reads the data into an array, 
  /// which is initialised to be the right size for the data.
  ///
  /// Note that this function can read in the full image, if global and local xmin
  /// are the same, and if npix is the full number of image pixels.  Alternately it
  /// can read in a brick-shaped subset based on the min starting point and npix.
  ///
  int read_fits_image_to_data(
        fitsfile *,     ///< fitsfile pointer.
        const string,   ///< name of hdu to read from.
        const int,      ///< dimensionality of image. 
        const double *, ///< local xmin (subdomain).
        const double *, ///< global xmin (full domain).
        const double,   ///< pixel size
        const int *,     ///< number of pixels to read in each direction
        const long int ,      ///< total number of pixels to be read.
        double **       ///< data array to write to.
        );
};



// ##################################################################
// ##################################################################



///
/// Class for reading and writing uniform grid data to fits images.
///
class DataIOFits : public utility_fitsio, public DataIOBase {
  public:
  ///
  /// Constructor.
  ///
  DataIOFits();

  /// Destructor (doensn't have much to do).
  virtual ~DataIOFits();

  ///
  /// Class can run with or without a solver, but this function
  /// allows you to set a pointer to the solver.
  ///
  void SetSolver(
        FV_solver_base * ///< Pointer to the solver (to get Eint,divB,Ptot)
        );

  ///
  /// Choose filename based on counter and base-name.
  ///
  virtual std::string choose_filename(
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
  virtual int OutputData(
        const string, ///< File to write to
        class GridBaseClass *, ///< pointer to data.
        const long int ///< number to stamp file with (e.g. timestep)
        );

  ///
  /// Reads the fits header from the file specified, and puts the
  /// simulation parameters in the global data class declared in global.h.
  /// This should be called first if we want to set up a uniform grid with
  /// the right dimensions to take the data from the fits file.
  ///
  virtual int ReadHeader(
        std::string ///< file to read from
        );

  ///
  /// This reads the fits images in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the fits header, which should be read first.
  ///
  virtual int ReadData(
        std::string, ///< file to read from
        class GridBaseClass * ///< pointer to data.
        );

  ///
  /// Write simulation header info to file (file must already exist!)
  ///
  int WriteHeader(
          const std::string ///< file to write to (full, exact filename).
          );

  protected:
  class GridBaseClass *gp; ///< pointer to computational grid.
  class FV_solver_base *eqn; ///< pointer to the solver, which knows the equations we are solving.
  fitsfile *file_ptr;

  ///
  /// Function which defines how to get the data from a silo file.
  ///
  int read_header_param(class pm_base *);

  ///
  /// Function which defines how to write a parameter to a silo file.
  ///
  int write_header_param(class pm_base *);

  ///
  /// Puts a grid variable into a 1D data array.
  ///
  int put_variable_into_data_array(
        const string,   ///< variable name to put in array.
        const long int, ///< size of data array to be initialised.
        double **       ///< pointer to uninitialised data.
        );

  ///
  /// Given that we have opened a fits file, read in the header, 
  /// checked the image is the right size, this reads the data into the
  /// local grid.
  /// 
  /// For parallel I/O, the local xmin is the local processor domain, while
  /// global is the xmin for the full image, and npix is not necessarily the
  /// size of the full image if I am only reading part of it for the local
  /// processor.
  ///
  /// This function calls the utility_fitsio() function of the same name, which 
  /// reads the image into an array.  This function then puts that array into
  /// a grid data variable.
  ///
  int read_fits_image(
        fitsfile *, ///< fitsfile pointer.
        string,     ///< extname of hdu to read from.
        double *,   ///<  vector local xmin (subdomain).
        double *,   ///< vector global xmin (full domain).
        int *,      ///< vector npix to read in each direction.
        long int    ///< total npix to read.
        );
};

#endif // if FITS

#endif // DATAIO_FITS_H
