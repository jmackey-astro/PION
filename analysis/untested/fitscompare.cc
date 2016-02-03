/** \file fitscompare.cc
 * \author Jonathan Mackey
 * 
 * This file reads in two fits files and compares them,
 * outputting the difference to a new fits file.
 * */

#include "fitsio.h"
using namespace std;
#include "../source/global.h"
//#include "../testing/grid.h"

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc<4) {
    cerr << "Error: must call with at least 3 arguments...\n";
    cerr << "fitscompare: <fitscompare> <infile1> <infile2> <outfile>\n";
    exit(1);
  }
  string infile1 = argv[1];
  string infile2 = argv[2];
  string outfile = "!"; outfile += argv[3];
  cout <<"reading from files "<<infile1<<" and "<<infile2<<endl;
  cout <<"Writing to file "<<outfile<<endl;
  cout <<"Note the output file contains the absolute error, as opposed to relative error.\n";
  cout <<"You should use \'fdiff\' from the command line to see if the headers differ.\n";
  // Open two input files and output file.
  fitsfile *ff1,*ff2,*ffout;
  int status=0,err=0;
  cout <<"opening fits files...";
  err = fits_open_file(&ff1, infile1.c_str(), READONLY, &status);
  if(status) {fits_report_error(stderr,status); return(err);}
  err = fits_open_file(&ff2, infile2.c_str(), READONLY, &status);
  if(status) {fits_report_error(stderr,status); return(err);}
  fits_create_file(&ffout, outfile.c_str(), &status);
  if(status) {cerr<<"outfile open went bad.\n";exit(1);}
  cout <<"done.\n";
  
  // Get both input files to hdu2:
  int num;
  fits_get_hdu_num(ff1, &num);
  if (num !=1) err += ffmahd(ff1,1,0,&status);
  err += ffmrhd(ff1,1,0,&status);   fits_get_hdu_num(ff1, &num);
  if (status) {fits_report_error(stderr,status); return(err);} 
  //  cout <<"ff1: Current hdu: "<<num<<"\t and err="<<err<<endl;

  fits_get_hdu_num(ff2, &num);
  if (num !=1) err += ffmahd(ff2,1,0,&status);
  err += ffmrhd(ff2,1,0,&status);   fits_get_hdu_num(ff2, &num);
  if (status) {fits_report_error(stderr,status); return(err);} 
  //  cout <<"ff2: Current hdu: "<<num<<"\t and err="<<err<<endl;

  
  int ndim;
  int bitpix1, bitpix2;
  int naxis1; long int *naxes1 = new long int [3];
  int naxis2; long int *naxes2 = new long int [3];
  double nulval = -1.e99; int anynul=0;
  
  bool arraysetup=false;
  double *array1=0, *array2=0, *array3=0;
  
  // loop over image hdus.
  do {
    // open hdu i.  get dimensions. check that they match.
    err += fits_get_img_param(ff1,0, &bitpix1,&naxis1,naxes1,&status);
    err += fits_read_keys_lng(ff1,"naxis",1,naxis1,naxes1,&num,&status); // reads ndim keys matching naxis, returns number found.
    if (status) {fits_report_error(stderr,status); return(err);}
    if (bitpix1 != DOUBLE_IMG) rep.error("Bad image type",bitpix1);
    ndim = naxis1;
    err += fits_get_img_param(ff2,0, &bitpix2,&naxis2,naxes2,&status);
    err += fits_read_keys_lng(ff2,"naxis",1,naxis2,naxes2,&num,&status); // reads ndim keys matching naxis, returns number found.
    if (status) {fits_report_error(stderr,status); return(err);}
    if (bitpix2 != DOUBLE_IMG) rep.error("Bad image type",bitpix1);
    if (naxis1 != naxis2)  rep.errorTest("Images not same dimensionality!",naxis1,naxis2);
    long int npix=1;
    for (int j=0;j<ndim;j++){ 
      if (naxes1[j]!=naxes2[j])
	rep.error("Bad image length in at least one direction",naxes1[j]);
      npix *=naxes1[j]; // cout <<"npix="<<npix<<endl;
    }
    
    if (!arraysetup) {
      array1 = new double[npix];
      array2 = new double[npix];
      array3 = new double[npix];
      if (!array1 || !array2 || !array3) rep.error("memory allocation.",array3);
      arraysetup = true;
    }
    
    // images match, so set output hdu to be same size.
    // This should start a new hdu and give it the header of the CHDU of ff1.
    err += fits_copy_header(ff1,ffout,&status);
    if (status) {fits_report_error(stderr,status); return(err);}
    
    // read data from file1, hdu i, into array1.
    err += fits_read_img(ff1,TDOUBLE,1,npix,&nulval,array1,&anynul,&status);
    if(status) {fits_report_error(stderr,status); return(err);}
  // read data from file2, hdu i, into array2.
    err += fits_read_img(ff2,TDOUBLE,1,npix,&nulval,array2,&anynul,&status);
    if(status) {fits_report_error(stderr,status); return(err);}
  // write (d1-d2)/(d1+d2) (or 0 if d1,d2=0) into array3.
    double maxval=-1.e99, minval=1.e99; double meanval=0;
    for (long int j=0; j<npix; j++) {
//      array3[j] = (fabs(array1[j])-fabs(array2[j]))/(fabs(array1[j])+fabs(array2[j]));
//      if ( (fabs(array1[j])+fabs(array2[j])) < 1.e-50) array3[j] = 0.;
      meanval += (array1[j]+array2[j])/2.;
      array3[j] = array1[j]-array2[j];
      if (array3[j]>maxval) maxval=array3[j];
      if (array3[j]<minval) minval=array3[j];
    }
    cout.setf( ios_base::scientific,ios_base::floatfield ); cout.precision(6);
    cout <<"max diff = "<<maxval<<"\t and min diff = "<<minval<<"\t and meanval="<<meanval/npix<<endl;
  // write array3 to hdu i in output file3.
    for (int i=0;i<ndim;i++) naxes1[i]=1;
    err += fits_write_pix(ffout, TDOUBLE, naxes1, npix, array3, &status);
    if(status) fits_report_error(stderr,status);
    
    // Move to next hdu
    err += ffmrhd(ff1,1,0,&status);
    err += ffmrhd(ff2,1,0,&status);
    fits_get_hdu_num(ff1, &num);
    // cout <<"ff1: Current hdu: "<<num<<"\t and err="<<err<<endl;
    fits_get_hdu_num(ff2, &num);
    // cout <<"ff2: Current hdu: "<<num<<"\t and err="<<err<<endl;
  } while (!status); 
  // Close files.
  
  err += fits_close_file(ff1,&status);
  err += fits_close_file(ff2,&status);
  err += fits_close_file(ffout,&status);
  
  if(arraysetup) {
    delete [] array1;
    delete [] array2;
    delete [] array3;
  }
  delete [] naxes1;
  delete [] naxes2;
  
  return(0);
}




