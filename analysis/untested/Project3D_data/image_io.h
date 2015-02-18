///
/// \file image_io.h
/// \author Jonathan Mackey
/// \date 20XX
///
/// Description: Writes images in ASCII, FITS, VTK format.
///
/// Modifications:
/// - 2013.10.15 JM: Added VTK output support.

#ifndef IMAGE_IO
#define IMAGE_IO


//
// Includes for inheritance, and fits definitions:
//
#include "dataIO/dataio_fits.h"
#include <fitsio.h>

//
// Image I/O class, for writing text and fits images.
//
class image_io: public file_status {
public:
  image_io() {
    ff=0;
    fvtk=0;
    vtk_header_written=false;
  }
  ~image_io() {}

  string get_output_filename(
        const string , ///< outfile base
        const int,     ///< multi_opfiles: 0=single file, 1=multiple files.
        const int,     ///< op_filetype: 0=textfile, 1=fitsfile
        const int      ///< integer count of file number to put in name.
        );

  int open_image_file(
        const string, ///< this_outfile
        const int,    ///< op_filetype: 0=textfile, 1=fitsfile
        string *      ///< pointer to string handle for file.
        );

  int close_image_file(
        const string ///< string handle for file to close.
        );

  int write_image_to_file(
        const string,   ///< string handle for file.
        const int,      ///< type of file (0,1,3)
        double *,       ///< image array
        const int,      ///< number of pixels: image.num_pixels
        const int,      ///< dimension of image [0,1,2]: image.dim
        const int *,    ///< array of npix in each dimension: image.npix
        const string,   ///< Name of image, e.g. TotMass (8 char max. for fits).
        const double *, ///< array with minimum dimensions of image.
        const double *, ///< array with pixel size in each dimension.
        const double,   ///< simulation time (optional)
        const long int  ///< simulation timestep (optional)
        );

protected:
  string filehandle; ///< just one for now, but should be able to handle many.
  fitsfile *ff;  ///< Fits file pointer.
  FILE *fvtk;    ///< VTK file pointer (standard c pointer).
  bool vtk_header_written; ///< whether the VTK header has been written or not.
  ofstream outf; ///< Ascii text file I/O stream.
  class utility_fitsio utfits; ///< utility fits i/o class for writing images.

  ///
  /// Swap the endian-ness of data (from PLUTO, Mignone).
  ///
  void SWAP_ENDIAN (
        void *,   ///< data.
        const int ///< number of bytes.
        );
};


#endif // IMAGE_IO
