/// \file utility_fits_class.cc
///
/// This file contains the class definitions for the utility_fitsio class, which
/// uses FITS, from
/// http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html at NASA.  I
/// have tried to make the functions atomic as much as possible, so bits of code
/// can be reused.  This is mostly used by DataIOFits class, in dataio_fits.cc,
/// with whom it shares a header file (dataio_fits.h).  This class is intended
/// to be grid and solver independent -- i.e. it just does things with fits
/// files and arrays of data.
///
/// FITS functions have a short name, typically ffxxxx(), and a long name,
/// typically fits_do_some_task().  I always use the long name, except in cases
/// where it isn't defined and I had to use ffmahd(ff,1,0,&status), which is
/// fits_move_absolute_hdu() and moves to the numbered HDU in the argument. Also
/// ffmrhd(ff,N,0,&status) moves forward by N HDU's if possible.
///
/// modified:\n
///  - 2009-06-06 Split into two classes: a base utility class that knows
///  nothing about the grid,
///     and a DataIOFits class which interfaces with the grid.
///
///  - 2010-02-03 JM: made ~utility_fitsio() virtual and put in header file.
///    Re-named repeated variable definitions and removed unused variables.
///
/// - 2010.10.01 JM: Got rid of testing myalloc/myfree commands.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2012.10.15 JM: minor mods to reading fits data, so that I don't
///    need the HDU's name; without a name it will default to hdu1.
/// - 2013.04.16 JM: some debugging messages and new comments added.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.02.13 JM: improved read_fits_image_to_data() quite a lot.

#ifdef FITS

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "dataIO/dataio_fits.h"
#include "fitsio.h"
#include <cstring>
#include <sstream>
#include <vector>
using namespace std;

// ##################################################################
// ##################################################################

utility_fitsio::utility_fitsio() {}

// ##################################################################
// ##################################################################

int utility_fitsio::create_fits_image(
    fitsfile* ff, const string name, const int ndim, const int* npt)
{
    string fname = "utility_fitsio::create_fits_image";
    int status   = 0;
    long int pix[ndim];
    for (int j = 0; j < ndim; j++)
        pix[j] = npt[j];

    //  for (int j=0;j<ndim;j++) cout <<"  pix[j]="<<pix[j];
    // cout <<"\tff: "<<ff<<" imtype: "<<DOUBLE_IMG<<" ndim:"<<ndim<<"
    // status:"<<status<<" name:"<<name<<"\n";
    fits_create_img(ff, DOUBLE_IMG, ndim, pix, &status);
    // cout <<"\tff: "<<ff<<" imtype: "<<DOUBLE_IMG<<" ndim:"<<ndim<<"
    // status:"<<&status<<"\n";

    char temp[256];
    strcpy(temp, name.c_str());
    fits_write_key(ff, TSTRING, "extname", temp, "Image Name", &status);

    return status;
}

// ##################################################################
// ##################################################################

int utility_fitsio::write_fits_image(
    fitsfile* ff,
    const string,  /// This is the name of the ext, but it's not used at the mo.
    const double* localxmin,
    const double* globalxmin,
    const double pix_size,
    const int ndim,
    const int* npt,
    long int ntot,
    double* data)
{
    int status     = 0;
    long int *fpix = 0, *lpix = 0;
    fpix = mem.myalloc(fpix, ndim);
    lpix = mem.myalloc(lpix, ndim);

    //
    // First convert xmin,xmax into pixels, and create bottom-left and top-right
    // pixels.
    //
    long int npix = 1;
    for (int i = 0; i < ndim; i++) {
        // FITS uses unit offset values to address pixels. First pixel is
        // [1,1,1]
        fpix[i] = static_cast<long int>(
                      (localxmin[i] - globalxmin[i]) * ONE_PLUS_EPS / pix_size)
                  + 1;
        lpix[i] = fpix[i] + npt[i]
                  - 1;  // -1 because it's inclusive: fpix,pfix+1,...,lpix
        npix *= (lpix[i] - fpix[i] + 1);  // +1 because of previous line.
        //    cout <<"fpix[i],lpix[i] = "<<fpix[i]<<", "<<lpix[i]<<"\n";
    }
    if (npix != ntot)
        rep.error("Pixel counting failed in write_fits_image()", npix - ntot);

    // finally write the data.
    //  int fits_write_subset/ffpss(ff,TDOUBLE,long *fpix,
    //  long*lpix,*data,&status)
    // fpix is a vector of the bottom left pixel (unit offset).  lpix is the top
    // right pixel.  The function then writes data into that subset of the
    // image.
    fits_write_subset(ff, TDOUBLE, fpix, lpix, data, &status);
    if (status) fits_report_error(stderr, status);

    fpix = mem.myfree(fpix);
    lpix = mem.myfree(lpix);
    return status;
}

// ##################################################################
// ##################################################################

int utility_fitsio::check_fits_image_dimensions(
    fitsfile* ff,       ///< file pointer.
    const string name,  ///< Name of hdu image to check
    const int ndim,     ///< dimensionality we are expecting.
    const int* npix     ///< number of pixels in each direction.
)
{
    int status = 0, num1 = 0, num = 0;
    fits_get_hdu_num(ff, &num);
    cout << "Current hdu: " << num << "\t and extname = " << name << "\n";

    char* keyval = 0;
    keyval       = mem.myalloc(keyval, 256);
    strcpy(keyval, name.c_str());
    fits_movnam_hdu(ff, ANY_HDU, keyval, 0, &status);
    if (status) fits_report_error(stderr, status);
    fits_get_hdu_num(ff, &num1);
    cout << "Current hdu: " << num1 << "\t and extname = " << name << "\n";

    if (num1 != num) rep.error("Not in correct hdu for given extname", name);
    fits_read_keyword(ff, "extname", keyval, 0, &status);
    if (status) fits_report_error(stderr, status);
    printf("keyval= %s\n", keyval);
    //  string temp=keyval; cout <<"temp keyval = "<<temp<<"\n";
    //  if (extname != temp) rep.error("not in correct hdu!",(extname+=temp));

    int bitpix      = -1;
    int naxis       = -1;
    long int* naxes = 0;
    naxes           = mem.myalloc(naxes, ndim);
    fits_get_img_param(ff, 0, &bitpix, &naxis, naxes, &status);
    // cout <<"naxis="<<naxis<<", axes=["<<naxes[0]<<", "<<naxes[1]<<"],
    // status="<<status<<"\n";
    fits_read_keys_lng(
        ff, "naxis", 1, naxis, naxes, &num,
        &status);  // reads ndim keys matching naxis, returns number found.
    cout << "naxis=" << naxis << ", axes=[" << naxes[0] << ", " << naxes[1]
         << "], status=" << status << "\n";
    if (status) {
        fits_report_error(stderr, status);
        return (status);
    }
    //
    // Check that the image HDU has the right size:
    //
    if (bitpix != DOUBLE_IMG) rep.error("Bad image type", bitpix);
    if (naxis != ndim) rep.error("Bad image dimensionality!", naxis);
    for (int j = 0; j < naxis; j++)
        if (naxes[j] != npix[j]) {
            cout << "j=" << j << "  axes=" << naxes[j] << " npix=" << npix[j]
                 << "\n";
            rep.error(
                "Bad image length in at least one direction delta(N) follows:",
                naxes[j] - npix[j]);
        }

    naxes  = mem.myfree(naxes);
    keyval = mem.myfree(keyval);
    return status;
}

// ##################################################################
// ##################################################################

int utility_fitsio::read_fits_image_to_data(
    fitsfile* ff,           ///< fitsfile pointer.
    const string name,      ///< name of hdu to read from.
    const int ndim,         ///< dimensionality of image.
    const double* l_xmin,   ///< local xmin (subdomain).
    const double* g_xmin,   ///< global xmin (full domain).
    const double pix_size,  ///< pixel size
    const int* npt,         ///< number of pixels to read in each direction
    const long int ntot,    ///< total number of pixels to be read.
    const int datatype,     ///< FITS datatype (e.g. TDOUBLE or TFLOAT)
    void* input_data        ///< data array to write to.
)
{
    if (datatype != TDOUBLE && datatype != TFLOAT) {
        cout << "read_fits_image_to_data() only double or float.";
        return 1;
    }

    int hdu_num = 1;  // default to the first hdu.

    //
    // Set first and last pixel to be read from image.
    //
    int status = 0;
    // First convert xmin,xmax into pixels, and create bottom-left and top-right
    // pixels.
    long int fpix[ndim];
    long int lpix[ndim];
    long int inc[ndim];
    long int npix = 1;

    for (int i = 0; i < ndim; i++) {
        inc[i]  = 1;
        fpix[i] = static_cast<long int>(
                      (l_xmin[i] - g_xmin[i]) * ONE_PLUS_EPS / pix_size)
                  + 1;
        lpix[i] = fpix[i] + npt[i]
                  - 1;  // -1 because it's inclusive: fpix,fpix+1,...,lpix
        npix *= (lpix[i] - fpix[i] + 1);  // +1 because of previous line.
        cout << "fpix[i],lpix[i] = " << fpix[i] << ", " << lpix[i] << "\n";
    }
    if (npix != ntot)
        rep.error("Pixel counting failed in read_fits_image()", npix - ntot);

    //
    // Read subset of hdu to array.
    //
    // char keyval[256];
    // fits_read_keyword(ff,"extname",keyval,0,&status);
    // if (status) fits_report_error(stderr,status);
    //  cout <<"hdu extname: "<<keyval<<", and we are expecting: "<<name<<"\n";

    //
    // We shouldn't need to move to the hdu, b/c we should already be there, but
    // this is just an extra layer of checking.
    //
    int err = 0;
    if (name != "") {
        char keyval[256];
        strcpy(keyval, name.c_str());
        err = fits_movnam_hdu(ff, ANY_HDU, keyval, hdu_num, &status);
        if (status) {
            fits_report_error(stderr, status);
        }
        if (err) {
            // If can't find variable, set them all to zero.
            if (status) {
                fits_report_error(stderr, status);
            }
            fits_clear_errmsg();
            status = 0;
            cout << "utility_fitsio::read_fits_image() couldn't get data ";
            cout << "for variable " << name;
            cout << "; will return error code.\n";
            return err;
        }
    }
    else {
        err = ffmahd(ff, hdu_num, 0, &status);
        if (status) {
            fits_report_error(stderr, status);
            return status;
        }
        if (err) {
            cout << "problem with ffmahd, err=" << err << "\n";
            return err;
        }
    }

    void* nulval = 0;
    double dnull = -1.0e20;
    float fnull  = -1.0e20;
    if (datatype == TDOUBLE)
        nulval = static_cast<void*>(&dnull);
    else
        nulval = static_cast<void*>(&fnull);

    int anynul = 0;
    fits_read_subset(
        ff, datatype, fpix, lpix, inc, nulval, input_data, &anynul, &status);
    if (status) {
        fits_report_error(stderr, status);
        return status;
    }
    cout << "anynul = " << anynul << "\n";

    return 0;
}

// ##################################################################
// ##################################################################

#endif  // ifdef FITS
