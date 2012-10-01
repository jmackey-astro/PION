/** \file stitchfits.cc
 * \author Jonathan Mackey
 * 
 * This file reads in fits files generated by parallel code and stitches the
 * partial files together into a single file.
 * input files are of name filename_procno.timestep.fits
 * */

#ifndef PARALLEL
# error "define PARALLEL so this will work!"
#endif

#include <sstream>
#include "fitsio.h"
using namespace std;
#include "../testing/global.h"
#include "../testing/dataio.h"
#include "../testing/dataio_fits.h"

int main(int argc, char **argv)
{
  // Get an input file, an output file, and a step.
  if (argc!=6) {
    cerr << "Error: must call as follows...\n";
    cerr << "stitchfits: <stitchfits> <infilebase> <outfile> <nproc> <startstep> <step>\n";
    exit(1);
  }
  string infilebase  = argv[1];
  string outfilebase = argv[2];
  mpiPM.nproc = atoi(argv[3]);
  int start= atoi(argv[4]);
  int step = atoi(argv[5]);
  cout <<"reading from file base "<<infilebase<<endl;
  cout <<"Writing to file "<<outfilebase<<endl;
  cout <<"Step between outputs is "<<step<<" timesteps.\n";
  cout <<"number of processors: "<<mpiPM.nproc<<endl;
  cout <<"**********************************************\n";
  // Open two input files and output file.
  class DataIOFits dataio;
  class file_status fs;
  
  fitsfile *ffin,*ffout;
  int status=0,err=0;
  ostringstream temp; string infile,outfile;
  
  // Need to loop this over all timesteps, incrementing 'start' by 'step' each time
  // until there are no more files to stitch together.
  do {
     temp<<infilebase<<"_0."<<start<<".fits";
     infile = temp.str();
     temp.str(""); temp <<outfilebase<<"."<<start<<".fits";
     outfile = temp.str();
     cout <<"\n*****************************************************\n";
     cout <<"Opening fits files for timestep index "<<start<<"\n";
     err = fits_open_file(&ffin, infile.c_str(), READONLY, &status);
     if(status) {fits_report_error(stderr,status); return(err);}
     if (fs.file_exists(outfile)) {
       temp.str(""); temp << "!" << outfile; outfile=temp.str();
       cout <<"Output file exists!  hopefully this is ok.\n";
     }
     fits_create_file(&ffout, outfile.c_str(), &status);
     if(status) {cerr<<"outfile open went bad.\n";exit(1);}
     cout <<"infile[0]: "<<infile<<"\nand outfile: "<<outfile<<endl;
     
     // This should set mpiPM.LocalXmin and SimPM.Xmin/max/range/NG
     mpiPM.myrank = 0;
     err = dataio.ReadHeader(infile);
     if (err) rep.error("Didn't read header",err);
     
     // Outfile:
     //  - copy header from first infile.
     err = fits_copy_header(ffin,ffout,&status);
     if (status) {fits_report_error(stderr,status); return(err);}
     
     //  - for each hdu in proc 0 infile, create full size image hdu in outfile
     int num; err = fits_get_num_hdus(ffin,&num,&status);
     if (status) {fits_report_error(stderr,status); return(err);}
     long int *pix = new long int [SimPM.ndim];
     for (int j=0;j<SimPM.ndim;j++) pix[j]=SimPM.NG[j];
     cout <<"\tCreating "<<num<<" hdus in outfile (including header)...";
     for (int i=2;i<=num;i++) { // hdu1 is header.
       fits_create_img(ffout,DOUBLE_IMG,SimPM.ndim,pix,&status);
       // get hdu name from ffin
       ffmahd(ffin,i,0,&status);
       if (status) {
	 cout <<"move to infile hdu "<<i<<endl;
	 fits_report_error(stderr,status); return(err);
       }
       char keyval[32];
       fits_read_keyword(ffin,"extname",keyval,0,&status);
       string sss=keyval;
       //cout <<"sss="<<sss;
       int length = sss.length(); sss=sss.substr(1,length-2);
       //cout <<" sss="<<sss<<endl;
       strcpy(keyval,sss.c_str());
       fits_write_key(ffout, TSTRING, "extname",  keyval, "Image Name", &status);
       if (status) {
	 cout <<"reading/writing extname for hdu "<<i<<endl;
	 fits_report_error(stderr,status); return(err);
       }
     }
     cout <<" done.\n";
     delete [] pix;
     err += fits_close_file(ffin,&status);
     if(status) {fits_report_error(stderr,status); return(err);}
     
     cout <<"\tNow copying infile hdus into outfile.\n";
     //  - for each hdu in proc i infile, copy image into appropriate place in oufile image.
     for (int proc=0; proc<mpiPM.nproc; proc++) {
       cout <<"\t\tproc "<<proc<<": ";
       mpiPM.myrank = proc;
       // read infile header to get local ng/xmin/xmax
       temp.str("");temp<<infilebase<<"_"<<proc<<"."<<start<<".fits";
       infile = temp.str();
       temp.str(""); temp <<outfilebase<<"."<<start<<".fits";
       if (!fs.file_exists(infile)) rep.error("infile doesn't exist",infile);
       err = dataio.ReadHeader(infile);
       if (err) rep.error("Didn't read header",err);
       err = fits_open_file(&ffin, infile.c_str(), READONLY, &status);
       if(status) {fits_report_error(stderr,status); return(err);}
       
       double *array=0;
       array = new double[mpiPM.LocalNcell];
       if (!array) rep.error("mem alloc",array);
       
       for (int im=2;im<=num;im++) { // hdu1 is header.
	 ffmahd(ffin, im,0,&status);
	 ffmahd(ffout,im,0,&status);
	 if (status) {
	   cout <<"move to infile/outfile hdu "<<im<<endl;
	   fits_report_error(stderr,status); return(err);
	 }
	 cout <<" hdu"<<im;
	 // read image from infile into array.
	 long int *fpix = new long int [SimPM.ndim];
	 long int *lpix = new long int [SimPM.ndim];
	 long int *inc  = new long int [SimPM.ndim]; // I think this is the increment in num. pix. per read.
	 long int npix=1; double dx = mpiPM.LocalRange[XX]/mpiPM.LocalNG[XX];
	 for (int i=0;i<SimPM.ndim;i++) {
	   inc[i] = 1;
	   fpix[i] = 1;
	   lpix[i] = mpiPM.LocalNG[i]; // It's inclusive: fpix,fpix+1,...,lpix
	   npix *= (lpix[i]-fpix[i]+1);  // +1 because of previous line.
	   //    cout <<"fpix[i],lpix[i] = "<<fpix[i]<<", "<<lpix[i]<<endl;
	 }
	 if (npix != mpiPM.LocalNcell) {
	   cout <<"ncell = "<<mpiPM.LocalNcell<<" and counted "<<npix<<" cells.\n";
	   rep.error("Pixel counting failed in Image Read",npix-mpiPM.LocalNcell);
	 }
	 double nulval = -1.e99; int anynul=0;
	 // Read data from image.
	 fits_read_subset(ffin, TDOUBLE, fpix, lpix, inc, &nulval, array, &anynul, &status);
	 if (status) {fits_report_error(stderr,status); return status;}
	 
	 // write image into subset of output file.
	 npix = 1;
	 for (int i=0;i<SimPM.ndim;i++) {
	   fpix[i] = static_cast<long int>((mpiPM.LocalXmin[i]-SimPM.Xmin[i])/dx*1.00000001) +1;
	   lpix[i] = fpix[i] + mpiPM.LocalNG[i] -1; // -1 because it's inclusive: fpix,pfix+1,...,lpix
	   npix *= (lpix[i]-fpix[i]+1);  // +1 because of previous line.
	   //	   cout <<"fpix[i],lpix[i] = "<<fpix[i]<<", "<<lpix[i]<<endl;
	 }
	 if (npix != mpiPM.LocalNcell) {
	   cout <<"ncell = "<<mpiPM.LocalNcell<<" and counted "<<npix<<" cells.\n";
	   rep.error("Pixel counting failed in Image Write",npix-mpiPM.LocalNcell);
	 }
	 fits_write_subset(ffout, TDOUBLE, fpix, lpix, array, &status);
	 if(status) fits_report_error(stderr,status);
	 delete [] fpix; fpix=0;
	 delete [] lpix; lpix=0;
	 delete [] inc;  inc=0;
       } // copy all images. (int im)
       cout <<"\n";

       err += fits_close_file(ffin,&status);
       if(status) {fits_report_error(stderr,status); return(err);}
       delete [] array; array=0;
     } // loop over all infiles (int proc)
     cout <<"\tOutfile "<<outfile<<" written.  Moving to next step.\n";
     //  - close outfile.
     err += fits_close_file(ffout,&status);
     if(status) {fits_report_error(stderr,status); return(err);}
     
     start += step;
     temp.str("");
     temp<<infilebase<<"_0."<<start<<".fits";
     infile = temp.str(); temp.str("");
  } while (fs.file_exists(infile));
  // loop over all timesteps.
  cout <<"\n***************************************************\n";
  cout <<"couldn't find file "<<infile<<" for step "<<start<<"... assuming i'm finished!\n";

  if (COMM) {delete COMM; COMM=0;}

  return 0;
}
