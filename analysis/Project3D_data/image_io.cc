


#include "image_io.h"
#include "global.h"
#include <sstream>
using namespace std;

// -------------------------------------------------------------
// *************************************************************
// ********************** IMAGE I/O CLASS **********************
// *************************************************************
// -------------------------------------------------------------


string image_io::get_output_filename(const string outfile,    ///< outfile base
				     const int multi_opfiles, ///< multi_opfiles: 0=single file, 1=multiple files.
				     const int op_filetype,   ///< op_filetype: 0=textfile, 1=fitsfile
				     const int ifile          ///< integer count of file number to put in name.
				     )
{
  ostringstream temp;
  temp.str(""); temp << outfile;
  if (multi_opfiles) {
    temp.fill('0'); temp <<"_step"; temp.width(3); temp<<ifile;
  }
  if      (op_filetype==0) temp << ".txt";
  else if (op_filetype==1) temp << ".fits";
  else rep.error("bad op_filetype",op_filetype);
  string this_outfile = temp.str();
  filehandle = this_outfile;
  return this_outfile;
}

int image_io::open_image_file(const string this_outfile, ///< this_outfile
			      const int op_filetype,     ///< op_filetype: 0=textfile, 1=fitsfile
			      string *f                  ///< file handle
			      )
{
  (*f).clear();
  *f = filehandle;
  int status=0;
  string fitsfile;

  switch (op_filetype) {
  case 0:
    //
    // Text File
    //
    if (outf.is_open())
      rep.error("Output text file is open! Close file before opening a new one!",1);
    outf.open(this_outfile.c_str());
    if (!outf.is_open()) rep.error("Failed to open text file for writing",this_outfile);
    outf <<"# writing: simtime timestep M_tot M_neutral M_ionised M_ambient M_clump M_n1e3 M_n3e3 M_n1e4 M_n3e4 Nc1e3 Nc1e4 Nc3e4\n#\n";
    break;
  case 1:
    //
    // fits file
    //
    if (ff) rep.error("fits file is open! close it before opening a new one!",1);
    fitsfile = this_outfile;
    if (file_exists(fitsfile)) {
      cout << "\tfits file exists: " <<fitsfile<<" ...overwriting!\n";
      ostringstream temp; temp.str(""); temp <<"!"<<fitsfile; fitsfile=temp.str();
    }
    fits_create_file(&ff, fitsfile.c_str(), &status);
    if (status) rep.error("Creating new file went bad.\n",status);
    cout <<"***** Created fits file: "<< fitsfile<<", file pointer="<<ff<<endl;
    break;
  default:
    rep.error("Bad op_filetype in open_image_file",op_filetype);
  }

  return 0;
}

int image_io::close_image_file(const string f)
{
  //
  // Check that f is equal to filehandle
  //
  if (f!=filehandle) rep.error("bad file handle in close_image_file()",f);
  
  int err=0, status=0;
  //
  // Close file
  //
  if (outf.is_open())
    outf.close();
  
  if (ff) {
    err += fits_close_file(ff,&status);
    ff=0;
  }

  return 0;
}

int image_io::write_image_to_file(const string f, ///< file handle
				  const int op_filetype, ///< fits or text?
				  double *im,  ///< image array
				  const int num_pix, ///< number of pixels: image.num_pixels
				  const int im_dim,  ///< dimension of image [0,1,2]: image.dim
				  const int *npix,   ///< array of npix in each dimension: image.npix
				  const string name  ///< Name of image, e.g. TotMass (8 char max. for fits).
				  )
{
  //
  // Check that f is equal to filehandle
  //
  if (f!=filehandle) rep.error("bad file handle",f);
  double min[im_dim];
  for (int v=0;v<im_dim;v++) min[v]=0.0;

  switch (op_filetype) {
  case 0:
    if (!outf.is_open()) rep.error("can't write image to unopened file.",f);
    outf <<"****************** image name: "<<name<< "***********************\n";
    outf <<"numpix="<<num_pix<<"\tndim="<<im_dim<<endl;
    for (int i=0; i<num_pix; i++) {
      if ( ((i+npix[0])%npix[0])==0) outf<<"\n";
      outf << (i+npix[0])%npix[0] <<"\t"<< i/npix[0] <<"\t"<< im[i] <<"\n";
    }
    outf <<"\n\n*************************************************\n\n";
    break;
  case 1:
    if (!ff) rep.error("Don't call write_image_to_file() without having an open file!!!",ff);
    cout <<"WRITING IMAGE TO FITS FILE: "<<f<<endl;
    utfits.create_fits_image(ff, name, im_dim, npix);
    utfits.write_fits_image(ff, name, min, min, 1.0, im_dim, npix, num_pix, im);
    break;
  default:
    rep.error("Bad op_filetype to write_image_to_file()",op_filetype);
  }
  return 0;
}

// -------------------------------------------------------------
// *************************************************************
// ********************** IMAGE I/O CLASS **********************
// *************************************************************
// -------------------------------------------------------------
