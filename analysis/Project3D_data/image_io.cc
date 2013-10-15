///
/// \file image_io.cc
/// \author Jonathan Mackey
/// \date 20XX
///
/// Description: Writes images in ASCII, FITS, VTK format.
///
/// Modifications:
/// - 2013.10.15 JM: Added VTK output support.




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
    temp.fill('0'); temp <<"_img."; temp.width(8); temp<<ifile;
  }
  if      (op_filetype==0) temp << ".txt";
  else if (op_filetype==1) temp << ".fits";
  else if (op_filetype==3) temp << ".vtk";
  else rep.error("bad op_filetype",op_filetype);
  string this_outfile = temp.str();
  filehandle = this_outfile;
  return this_outfile;
}



// ##################################################################
// ##################################################################



int image_io::open_image_file(
        const string this_outfile, ///< this_outfile
        const int op_filetype,     ///< op_filetype: 0=textfile, 1=fits,3=vtk
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
    outf <<"# projection code.\n#\n";
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

  case 3:
    //
    // VTK File, open using standard files.
    //
    if (fvtk)
      rep.error("Output VTK file pointer is not null!",fvtk);
    fvtk = fopen(this_outfile.c_str(),"wb");
    if (!fvtk) rep.error("error opening VTK file",this_outfile);
    break;

  default:
    rep.error("Bad op_filetype in open_image_file",op_filetype);
  }

  return 0;
}



// ##################################################################
// ##################################################################



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

  if (fvtk) {
    //
    // Set bool to false, so we write the header for the next file.
    //
    vtk_header_written=false;
    //
    // close file.
    //
    fclose(fvtk);
    fvtk=0;
  }


  return 0;
}



// ##################################################################
// ##################################################################



int image_io::write_image_to_file(
        const string f,      ///< file handle
        const int op_filetype, ///< text,fits,vtk
        double *im,          ///< image array
        const int num_pix,   ///< number of pixels: image.num_pixels
        const int im_dim,    ///< dimension of image [0,1,2]: image.dim
        const int *npix,     ///< array of npix in each dimension: image.npix
        const string name,   ///< Name of image, e.g. TotMass (8 char max. for fits).
        const double *xmin,  ///< array with minimum dimensions of image.
        const double *dx,    ///< array with pixel size in each dimension.
        const double time,   ///< simulation time (optional)
        const long int cycle ///< simulation timestep (optional)
        )
{
  //
  // Check that f is equal to filehandle
  //
  if (f!=filehandle) rep.error("bad file handle",f);
  double min[im_dim];
  for (int v=0;v<im_dim;v++) min[v]=0.0;
  ostringstream ofn;

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

  case 3:
    if (!fvtk) {
      rep.error("image_io::write_image_to_file() vtk file not open",
                fvtk);
    }
    //
    // write header if needed.
    //
    if (!vtk_header_written) {
      ofn.str("");
      ofn << "# vtk DataFile Version 2.0\n";
      ofn << "JM's VTK Data\n";
      ofn << "BINARY\n";
      ofn << "DATASET STRUCTURED_POINTS\n";

      ofn << "DIMENSIONS ";
      ofn <<               static_cast<int>(npix[0]+1) <<" ";
      if (im_dim>1) ofn << static_cast<int>(npix[1]+1) <<" ";
      else          ofn << "2 ";
      if (im_dim>2) ofn << static_cast<int>(npix[2]+1) <<"\n";
      else          ofn << "2 \n";
      
      ofn << "ORIGIN "  << xmin[0] <<" ";
      if (im_dim>1) ofn << xmin[1] <<" ";
      else          ofn << "0.0 ";
      if (im_dim>2) ofn << xmin[2] <<"\n";
      else          ofn << 0.0 <<"\n";

      ofn << "SPACING "      << dx[0];
      if (im_dim>1) ofn <<" "<< dx[1];
      else          ofn <<" "<< dx[0];
      if (im_dim>2) ofn <<" "<< dx[2];
      else          ofn <<" "<< dx[0];

      ofn << "FIELD FieldData 2\n";
      ofn << "TIME 1 1 double\n";
      fprintf (fvtk, "%s", ofn.str().c_str());

      double tempd=time;
      SWAP_ENDIAN(&tempd, sizeof(double));
      fwrite(&tempd, sizeof(double),1, fvtk);

      ofn.str("");
      ofn << "\nCYCLE 1 1 int\n";
      fprintf (fvtk, "%s", ofn.str().c_str());

      int tempi=cycle;
      SWAP_ENDIAN(&tempi, sizeof(int));
      fwrite(&tempi, sizeof(int),1, fvtk);

      ofn.str("");
      ofn << "\nCELL_DATA "<< static_cast<int>(num_pix) <<"\n";
      fprintf (fvtk, "%s", ofn.str().c_str());

      //
      // Set bool to true, so we only write it once.
      //
      vtk_header_written=true;
    }

    //
    // Now write image
    //
    ofn.str("");
    ofn << "SCALARS "<< name.c_str() <<" double\n";
    ofn << "LOOKUP_TABLE default\n";
    fprintf (fvtk, "%s", ofn.str().c_str());

    for (long int v=0;v<num_pix;v++)
      SWAP_ENDIAN(&(im[v]), sizeof(double));
    fwrite(im,sizeof(double),num_pix,fvtk);
    
    break;


  default:
    rep.error("Bad op_filetype to write_image_to_file()",op_filetype);
  }
  return 0;
}



// ##################################################################
// ##################################################################



void image_io::SWAP_ENDIAN(
        void *x,
        const int nbytes
        ) 
///
/// PURPOSE: Swap the byte order of x (legacy binary vtk must be
/// big-endian!).
/// \author Andrea Mignone (PLUTO code).
/// 
{
  if (nbytes>16) {cerr <<"swap-endian buffer is too large! Max=16bytes.\n";return;}
  int k;
  static char Swapped[16];
  char *c;
  c = (char *) x;
  for (k = nbytes; k--; ) Swapped[k] = *(c + nbytes - 1 - k);
  for (k = nbytes; k--; ) c[k] = Swapped[k];
  return;
}



// ##################################################################
// ##################################################################


// -------------------------------------------------------------
// *************************************************************
// ********************** IMAGE I/O CLASS **********************
// *************************************************************
// -------------------------------------------------------------
