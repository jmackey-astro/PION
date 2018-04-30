/// \file dataio_fits.cc
/// \author Jonathan Mackey
///
/// This file contains the class definitions for the DataIOFits class, which 
/// uses FITS, from http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
/// at NASA.  I have tried to make the functions atomic as much as possible,
/// so bits of code can be reused.
/// 
/// FITS functions have a short name, typically ffxxxx(), and a long
/// name, typically fits_do_some_task().  I always use the long name,
/// except in cases where it isn't defined and I had to use
/// ffmahd(ff,1,0,&status), which is fits_move_absolute_hdu() and moves
/// to the numbered HDU in the argument.  Also ffmrhd(ff,N,0,&status)
/// moves forward by N HDU's if possible.
///
/// modified:\n 
/// - 2007-10-25 got the parallel fits-io class to compile and write
///     data.
/// - 2007-10-26 parallel fits-io class reads data from single and
///     multiple files. Put in serial ifdefs, so that it should work
///     for serial code too.
/// - 2008-09-19 Renamed to dataio_fits.cc, and added in ifdefs so
///     code can be compiled on machines without cfitsio installed.
/// - 2009-06-06 Split into two classes: a base utility class that
///     knows nothing about the grid, and a DataIOFits class which
///     interfaces with the grid.
/// - 2010-02-03 JM: removed unused variables; renamed some variables
///    where i used 'i' twice in a function.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010.07.21 JM: order of accuracy flags are now integers.
/// - 2010.07.22 JM: new parameter read/write routines.
/// - 2010.07.23 JM: removed obselete read_fits_header(),
///    write_fits_header() functions.
/// - 2010.10.01 JM: Got rid of testing myalloc/myfree commands.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2010.11.03 JM: Added Ndigits variable for the step-counter in filename.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.02.25 JM: removed HCORR ifdef flags for column density.
/// - 2011.06.02 JM: Added WriteHeader() function so I can over-write header
///    parameters and restart a simulation with e.g. different microphysics.
///    Replaced fits_write_key() calls with fits_update_key() calls.  The new
///    function updates a key if it exists, and creates it otherwise.  The
///    old function would duplicate keys, which is not good for replacing the
///    header.
/// - 2013.04.16 JM: Added spacer lines between functions.
/// - 2013.04.18 JM: wrapped column-density output in ifdef so it is
///    not written to file by default.
/// - 2013.08.20 JM: Modified cell_interface for optical depth vars.
/// - 2014.10.21 JM: Fixed memory leak where array was freed in the
///    wrong place.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.01.28 JM: Removed parallel code, put into new class.
/// - 2015.02.13 JM: read_fits_image_to_data() no longer allocates
///    memory for the image, so I do it in read_fits_image().
/// - 2018.01.24 JM: worked on making SimPM non-global


#ifdef FITS

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"

#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "dataIO/dataio_fits.h"

#include "microphysics/microphysics_base.h"

#ifdef RT_TESTING_OUTPUTCOL
#include "raytracing/raytracer_base.h"
#endif // RT_TESTING_OUTPUTCOL

#include "fitsio.h"
#include <cstring>
#include <sstream>
#include <vector>
using namespace std;


// ##################################################################
// ##################################################################



//-------------------------------------------------
//-------------  FITS DATA I/O CLASS  -------------
//-------------------------------------------------


// ##################################################################
// ##################################################################


DataIOFits::DataIOFits(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
: DataIOBase(SimPM)
{
  cout <<"Setting up DataIOFits class.\n";
  DataIOFits::eqn =0;
  DataIOFits::gp=0;
  DataIOFits::file_ptr=0;

}


// ##################################################################
// ##################################################################



DataIOFits::~DataIOFits()
{
  cout <<"Deleting DataIOFits class.\n";
  DataIOFits::eqn =0;
  DataIOFits::gp=0;
}


// ##################################################################
// ##################################################################



void DataIOFits::SetSolver(
      FV_solver_base *solver
      )
{
  cout <<"DataIOFits::SetSolver() Setting solver pointer.\n";
  DataIOFits::eqn = solver;
}



// ##################################################################
// ##################################################################


int DataIOFits::OutputData(
      string outfilebase,          ///< base filename
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM, ///< pointer to simulation parameters
//      class MCMDcontrol  *mpiPM,   ///< pointer to multi-core params
      const long int file_counter  ///< number to stamp file with (e.g. timestep)
      )
{
  string fname="DataIOFits::OutputData";
  int nvar = SimPM.nvar;
  string *extname=0;
  if (SimPM.ntracer>5) 
    rep.error("OutputFitsData:: only accepts <=5 tracers!",SimPM.ntracer);
#ifdef RT_TESTING_OUTPUTCOL
  // save column densities
  if (SimPM.RS.Nsources>0) {
    for (int si=0; si<SimPM.RS.Nsources; si++) {
      nvar += SimPM.RS.sources[si].NTau;
    }
  }
#endif // RT_TESTING_OUTPUTCOL
  if (SimPM.eqntype==EQEUL || SimPM.eqntype==EQEUL_ISO || SimPM.eqntype==EQEUL_EINT) {
    extname=mem.myalloc(extname,nvar+1);
    string pvar[10] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","TR0","TR1","TR2","TR3","TR4"};
    for (int i=0;i<SimPM.nvar;i++) extname[i] = pvar[i];
    extname[nvar  ] = "Eint";
    if (DataIOFits::eqn!=0 || MP!=0) { // only write eint/divb,ptot if eqn is there to calculate it!
      nvar +=1;
    }
  }
  else if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQFCD) {
    extname=mem.myalloc(extname,nvar+3);
    string pvar[13] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","Bx","By","Bz","TR0","TR1","TR2","TR3","TR4"};
    for (int i=0;i<SimPM.nvar;i++) extname[i] = pvar[i];
    extname[nvar  ] = "Eint";
    extname[nvar+1] = "divB";
    extname[nvar+2] = "Ptot";
    if (DataIOFits::eqn!=0) {
      nvar +=3;
    }
  }
  else if (SimPM.eqntype==EQGLM) {
    extname=mem.myalloc(extname,nvar+3);
    string pvar[14] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","Bx","By","Bz","psi","TR0","TR1","TR2","TR3","TR4"};
    for (int i=0;i<SimPM.nvar;i++) extname[i] = pvar[i];
    extname[nvar  ] = "Eint";
    extname[nvar+1] = "divB";
    extname[nvar+2] = "Ptot";
    if (DataIOFits::eqn!=0) {
      nvar +=3;
    }
  } 
  else {
    extname=mem.myalloc(extname,10);
    rep.error("What equations?!",SimPM.eqntype);
  }
#ifdef RT_TESTING_OUTPUTCOL
  // save column densities
  if (SimPM.RS.Nsources>0) {
    if (extname[SimPM.nvar]!="")
      rep.error("Tau not writeable!",extname[SimPM.nvar]);
    //
    // Loop over all sources, and all variables for each source:
    //
    ostringstream var;
    unsigned int ivar=SimPM.nvar;
    for (int v=0;v<SimPM.RS.Nsources;v++) {
      for (int iT=0; iT<SimPM.RS.sources[v].NTau; iT++) {
        var.str("");
        var << "Col_Src_" << v << "_T"<<iT;
        extname[ivar] = var.str();
        ivar++;
      } // loop over Tau variables for source v.
    } // loop over Nsources
  } // if RT
#endif // RT_TESTING_OUTPUTCOL





  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the nested grid.
    CI.set_dx(SimPM.nest_levels[l].dx);
    ostringstream temp; temp << outfile << "_level";
    temp.width(2); temp.fill('0');
    temp << l;

  if (!cg[l])
    rep.error("DataIOFits::OutputData() null pointer to grid!",cg[l]);
  DataIOFits::gp = cg[l];

  fitsfile *ff=0;
  int status=0,err=0;



  //
  // Choose filename based on the basename and the counter passed to
  // this function.
  //
  string outfile = choose_filename(outfilebase, file_counter);
  ostringstream temp; temp.str("");

  // -------------------------------------------------------
  // -------------------------------------------------------
  //cout <<"DataIOFits::OutputData() writing file.";
  //cout <<":\t writing to file "<<outfile;

  if (file_exists(outfile)) {
    //cout <<":\t file exists... overwriting!... ";
    temp.str("");
    temp <<"!"<<outfile; outfile=temp.str();
  }

  // Create fits file.
  fits_create_file(&ff, outfile.c_str(), &status);
  if(status) rep.error("Creating new file went bad.",status);
  status=0;

  // write fits header
  //  err += write_fits_header(ff);
  //if(err) rep.error("DataIOFits::OutputData() couldn't write fits header",err);

  // --------------------------------------------------------
  //
  // create HDU for header
  //
  fits_create_img(ff,DOUBLE_IMG,0,0,&status);
  if(status) {fits_report_error(stderr,status);}
  //
  // set file pointer for the header writing function.
  //
  file_ptr=ff;
  err = write_simulation_parameters(SimPM);
  if (err) rep.error("DataIOFits::OutputData() couldn't write fits header",err);
  ff=file_ptr;
  // --------------------------------------------------------

  //
  // for each image, create image and write my portion of it.
  //
  for (int i=0;i<nvar;i++) {
    //cout <<extname[i]<<"\t"<<SimPM.Xmin[0]<<"\t"<<SimPM.Xmax[0]<<"\t"<<SimPM.NG[0]<<"\t"<<SimPM.Ncell<<"\n";
    //cout <<extname[i]<<"\t"<<SimPM.Xmin[1]<<"\t"<<SimPM.Xmax[1]<<"\t"<<SimPM.NG[1]<<"\t"<<SimPM.Ncell<<"\n";
    //cout <<extname[i]<<"\t"<<SimPM.Xmin[2]<<"\t"<<SimPM.Xmax[2]<<"\t"<<SimPM.NG[2]<<"\t"<<SimPM.Ncell<<"\n";
    err += create_fits_image(ff,extname[i], SimPM.ndim, SimPM.NG);
    double *data=0;
    err += put_variable_into_data_array(SimPM, extname[i], SimPM.Ncell, &data);
    err += write_fits_image(ff, extname[i], SimPM.Xmin, SimPM.Xmin, gp->DX(), SimPM.ndim, SimPM.NG, SimPM.Ncell, data);
    data = mem.myfree(data);
  }
  if (err) rep.error("DataIOFits::OutputData() Serial Image Writing went bad",err);

  // Close file
  err += fits_close_file(ff,&status);
  //  release_lock(outfile);
  //  cout <<": file written.\n";
  // -------------------------------------------------------
  // -------------------------------------------------------

  extname = mem.myfree(extname);
  
  if (status) {
    fits_report_error(stderr,status);
    fits_clear_errmsg();
    return(status);
  }

  return err;
}



// ##################################################################
// ##################################################################


int DataIOFits::ReadHeader(
      string infile, ///< file to read from
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  int err=0; int status=0; fitsfile *ff;
//  cout <<"DataIOFits::ReadHeader() opening fits file to read header...";
  err = fits_open_file(&ff, infile.c_str(), READONLY, &status);
  if(status) {fits_report_error(stderr,status); return(err);}

  //  cout <<"done. Now reading header.\n";
  //  err += read_fits_header(ff);
  //  if (err) 
  //   rep.error("DataIOFits::ReadHeader() read fits header failed.",err);

  //
  // Get to first hdu, which is always the header
  //
  int num=-1; fits_get_hdu_num(ff, &num);
  /* int fits_moveabs_hdu(ff, hdu-wanted, int *hdutype [IMAGE_HDU], &status); */
  if (num !=1) ffmahd(ff,1,0,&status);
  if (status) {fits_report_error(stderr,status); return(err);} 
  //
  // set file pointer for the header reading function.
  //
  file_ptr=ff;
  err += read_simulation_parameters(SimPM);
  if (err) 
    rep.error("DataIOFits::ReadHeader() read fits header failed.",err);
  ff=file_ptr;
  //
  // Now should be done, so close file and move on.
  //
  err += fits_close_file(ff,&status);
  if(status) {fits_report_error(stderr,status); return(err);}

  return (err);  
}



// ##################################################################
// ##################################################################




int DataIOFits::WriteHeader(
      const string fname, ///< file to write to (full, exact filename).
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  int err=0; int status=0; fitsfile *ff;
  err = fits_open_file(&ff, fname.c_str(), READWRITE, &status);
  if(status) {fits_report_error(stderr,status); return(err);}
  //
  // Get to first hdu, which is always the header
  //
  int num=-1; fits_get_hdu_num(ff, &num);
  /* int fits_moveabs_hdu(ff, hdu-wanted, int *hdutype [IMAGE_HDU], &status); */
  if (num !=1) ffmahd(ff,1,0,&status);
  if (status) {
    fits_report_error(stderr,status);
    //
    // Maybe header doesn't exist?  In this case create HDU for header
    //
    fits_create_img(ff,DOUBLE_IMG,0,0,&status);
    if(status) {fits_report_error(stderr,status); return status;}
  }

  //
  // set file pointer for the header writing function; then write the header.
  //
  file_ptr=ff;
  err = write_simulation_parameters(SimPM);
  if (err) rep.error("DataIOFits::OutputData() couldn't write fits header",err);
  ff=file_ptr;

  //
  // close file and return.
  //
  err = fits_close_file(ff,&status);
  if(status) {fits_report_error(stderr,status); return(err);}
  return (err);
}



// ##################################################################
// ##################################################################



int DataIOFits::ReadData(
      string infile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  string fname="DataIOFits::ReadData";

  if (!cg[0])
    rep.error("DataIOFits::ReadData() null pointer to grid!",cg[0]);
  DataIOFits::gp = cg[0];

  int err=0; int status=0; fitsfile *ff;
  //cout <<"DataIOFits::ReadData() opening fits file to read data...";
  err = fits_open_file(&ff, infile.c_str(), READONLY, &status);
  if(status) {fits_report_error(stderr,status); return(err);}
  //cout <<"done.\n";
  // Move to first data hdu; should be 2nd hdu;
  int num;
  fits_get_hdu_num(ff, &num);
  if (num !=1) err += ffmahd(ff,1,0,&status);
  if (status) {fits_report_error(stderr,status); return(err);} 
  err += ffmrhd(ff,1,0,&status);   fits_get_hdu_num(ff, &num);
  //  cout <<"Current hdu: "<<num<<"\t and err="<<err<<"\n";
  // -------------------------------------------------------

  int nvar = SimPM.nvar;
  string *var=0;
  if (SimPM.ntracer>5) rep.error("DataIOFits::ReadData() only handles up to 5 tracer variables! Add more if needed.",SimPM.ntracer);
  if (SimPM.eqntype==EQEUL || SimPM.eqntype==EQEUL_ISO || SimPM.eqntype==EQEUL_EINT) {
    var = mem.myalloc(var,10);
    if (SimPM.nvar>10) rep.error("DataIOFits::ReadData() need more tracers.",SimPM.nvar);
    string t[10] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","TR0","TR1","TR2","TR3","TR4"};
    for (int i=0;i<10;i++) var[i] = t[i];
  }
  else if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQGLM || SimPM.eqntype==EQFCD) {
    var = mem.myalloc(var,14);
    if (SimPM.nvar>14) rep.error("DataIOFits::ReadData() need more tracers.",SimPM.nvar);
    string t[14] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","Bx","By","Bz","psi","TR0","TR1","TR2","TR3","TR4"};
    for (int i=0;i<14;i++) var[i] = t[i];
  }
  else {
    var = mem.myalloc(var,10);
    rep.error("What equations?!",SimPM.eqntype);
  }
  // -------------------------------------------
  // Loop over all Variables and read from file.
  // -------------------------------------------
  char temp[32];
  for (int i=0; i<nvar; i++) {
    strcpy(temp,var[i].c_str());
    int v=0;
    if      (var[i]=="GasDens") v=static_cast<int>(RO);
    else if (var[i]=="GasPres") v=static_cast<int>(PG);
    else if (var[i]=="GasVX")   v=static_cast<int>(VX);
    else if (var[i]=="GasVY")   v=static_cast<int>(VY);
    else if (var[i]=="GasVZ")   v=static_cast<int>(VZ);
    else if (var[i]=="Bx")      v=static_cast<int>(BX);
    else if (var[i]=="By")      v=static_cast<int>(BY);
    else if (var[i]=="Bz")      v=static_cast<int>(BZ);
    else if (var[i]=="psi")     v=static_cast<int>(SI);
    else if (var[i]=="TR0")     {v=SimPM.ftr; cout <<"reading from first tracer var: "<<v<<"\n";}
    else if (var[i]=="TR1")     v=SimPM.ftr+1;
    else if (var[i]=="TR2")     v=SimPM.ftr+2;
    else if (var[i]=="TR3")     v=SimPM.ftr+3;
    else if (var[i]=="TR4")     v=SimPM.ftr+4;
    else rep.error("Bad variable index in fits read routine",var[i]);
    err += fits_movnam_hdu(ff,ANY_HDU,temp,0,&status);
    
    if (err!=0) {
      // If can't find variable, set them all to zero.
      cell *c = gp->FirstPt();
      do {c->P[v] = 0.;} while ( (c=gp->NextPt(c))!=0 );
      if(status) {fits_report_error(stderr,status);}
      err =0; fits_clear_errmsg(); status=0;
      cout <<"couldn't get data for variable "<<temp<<"; will set data to zero and hope for the best.\n";
    }
    
    else {
      // Variable found, check we're at the right hdu and read data.
      fits_get_hdu_num(ff, &num);
      //cout <<"Current hdu: "<<num<<"\t i="<<i<<" and var[i] = "<<var[i]<<"\n";
      //  cout <<"\t\tDataIOFits::ReadData() Reading fits image.\n";
      //  cout <<"\t\t reading from file "<<infile<<"\n";
      err += check_fits_image_dimensions(ff, var[i],  SimPM.ndim, SimPM.NG);
      if (err) rep.error("image wrong size.",err);
      //      cout <<"***************ncell = "<<SimPM.Ncell<<"\n";
      err += read_fits_image(SimPM, ff, var[i], SimPM.Xmin, SimPM.Xmin, SimPM.NG, SimPM.Ncell);
      if (err) rep.error("error reading image.",err);
      //  cout <<"\t\tDataIOFits::ReadData() Got fits image.\n";
    } // got real hdu and read data.

  } // Loop over all Primitive Variables


  var = mem.myfree(var);
  //  cout <<"Closing fits file. err="<<err<<"\n";
  fits_close_file(ff,&status);
  if (status) {
    fits_report_error(stderr,status);
    fits_clear_errmsg();
  }
  //  cout <<"Closed fits file. err="<<err<<"\n";

  // Now assign Ph to be equal to P for each cell.
  cell *cpt = gp->FirstPt();
  do {for(int v=0;v<nvar;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=gp->NextPt(cpt))!=0);
    
  return err;
}



// ##################################################################
// ##################################################################




///
/// Choose filename based on counter and base-name.
///
std::string DataIOFits::choose_filename(
          const std::string fbase,   ///< filebase passed in from main code.
          const int         file_counter ///< file counter to use (e.g. timestep).
          )
{
  //
  // Choose filename based on the basename and the counter passed to
  // this function.
  //
  string outfile;
  ostringstream temp; temp.str("");
  temp <<fbase;
  temp <<".";
  if (file_counter >=0) {
    temp.width(Ndigits); temp.fill('0');
    temp <<file_counter<<".";
  }
  temp <<"fits";
  outfile = temp.str();
  temp.str("");
  return outfile;
}



// ##################################################################
// ##################################################################







int DataIOFits::read_header_param(class pm_base *p)
{
  int err=0, status=0;
  char key[128];
  //
  // We read data to a temp var and copy it into the address held
  // by the parameter class.
  //
  int i=p->type;
  strcpy(key,p->name.c_str());
  if      (i==MY_INT) {
    int x;
    err += fits_read_key(file_ptr, TINT,    key, &x,   0, &status);
    p->assign_val(&x);
  }
  else if (i==MY_DOUBLE) {
    double x;
    err += fits_read_key(file_ptr, TDOUBLE, key, &x, 0, &status);
    p->assign_val(&x);
  }
  else if (i==MY_FLOAT) {
    float x;
    err += fits_read_key(file_ptr, TFLOAT,  key, &x, 0, &status);
    p->assign_val(&x);
  }
  else if (i==MY_LONG) {
    long int x;
    err += fits_read_key(file_ptr, TLONG,   key, &x, 0, &status);
    p->assign_val(&x);
  }
  else if (i==MY_STRING) {
    char x[128];
    err += fits_read_key(file_ptr, TSTRING, key, x,  0, &status);
    string temp(x);
    p->assign_val(&temp);
  }
  else if (i==MY_DDIMARR) {
    //
    // Easier to give each element of an array its own numbered name.
    //
    double x[MAX_DIM];
    for (int v=0;v<MAX_DIM; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_read_key(file_ptr, TDOUBLE, key, &(x[v]),0, &status);
    }
    p->assign_val(x);
  }
  else if (i==MY_IDIMARR) {
    int x[MAX_DIM];
    for (int v=0;v<MAX_DIM; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_read_key(file_ptr, TINT,    key, &(x[v]),0, &status);
    }
    p->assign_val(x);
  }
  else if (i==MY_DVARARR) {
    double x[MAX_NVAR];
    for (int v=0;v<MAX_NVAR; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_read_key(file_ptr, TDOUBLE, key, &(x[v]),0, &status);
    }
    p->assign_val(x);
  }

  if (status) {
    cout <<"\t"<<p->name<<":  ERROR READING VAR!";
    cout <<" err="<<err<<" status="<<status<<"\n";
  }
  // else {
  //   cout <<"\t"<<p->name<<":  "; p->show_val(); cout <<"\n";
  // }
  return status;
}



// ##################################################################
// ##################################################################


int DataIOFits::write_header_param(class pm_base *p)
{
  int err=0, status=0;
  char key[128];
  int i=p->type;
  strcpy(key,p->name.c_str());

  if      (i==MY_INT) {
    int *x = static_cast<int *>(p->get_ptr());
    err += fits_update_key(file_ptr, TINT,    key, x, 0, &status);
  }
  else if (i==MY_DOUBLE) {
    double *x = static_cast<double *>(p->get_ptr());
    err += fits_update_key(file_ptr, TDOUBLE, key, x, 0, &status);
  }
  else if (i==MY_FLOAT) {
    float *x = static_cast<float *>(p->get_ptr());
    err += fits_update_key(file_ptr, TFLOAT,  key, x, 0, &status);
  }
  else if (i==MY_LONG) {
    long int *x = static_cast<long int *>(p->get_ptr());
    err += fits_update_key(file_ptr, TLONG,   key, x, 0, &status);
  }
  else if (i==MY_STRING) {
    //
    // strings are harder -- need to get pointer and copy to char[]
    //
    string x(*(static_cast<string *>(p->get_ptr())));
    char temp[128]; strcpy(temp,x.c_str());
    err += fits_update_key(file_ptr, TSTRING, key, temp, 0, &status);
  }
  else if (i==MY_DDIMARR) {
    double *x = static_cast<double *>(p->get_ptr());
    //
    // Easier to give each element of an array its own numbered name.
    //
    for (int v=0;v<MAX_DIM; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_update_key(file_ptr, TDOUBLE, key, &(x[v]), 0, &status);
    }
   }
  else if (i==MY_IDIMARR) {
    int *x = static_cast<int *>(p->get_ptr());
    for (int v=0;v<MAX_DIM; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_update_key(file_ptr, TINT, key, &(x[v]), 0, &status);
    }
  }
  else if (i==MY_DVARARR) {
    double *x = static_cast<double *>(p->get_ptr());
    for (int v=0;v<MAX_NVAR; v++) {
      ostringstream temp2; temp2 << p->name << v;
      strcpy(key,(temp2.str()).c_str());
      err += fits_update_key(file_ptr, TDOUBLE, key, &(x[v]), 0, &status);
    }
  }

  if (status) {
    cout <<"\t"<<p->name<<":  ERROR WRITING VAR!";
    cout <<" err="<<err<<" status="<<status<<"\n";
  }
  // else {
  //   cout <<"\t"<<p->name<<":  "; p->show_val(); cout <<"\n";
  // }
  return status;
}


// ##################################################################
// ##################################################################




int DataIOFits::put_variable_into_data_array(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      const string name,   ///< variable name to put in array.
      const long int ntot, ///< size of data array to be initialised.
      double **data        ///< pointer to uninitialised data.
      )
{
  (*data) = mem.myalloc((*data),ntot);
  
  // Choose variable to write to, based on name string.
  int v = 0;
  if      (name=="GasDens") v=static_cast<int>(RO);
  else if (name=="GasPres") v=static_cast<int>(PG);
  else if (name=="GasVX")   v=static_cast<int>(VX);
  else if (name=="GasVY")   v=static_cast<int>(VY);
  else if (name=="GasVZ")   v=static_cast<int>(VZ);
  else if (name=="Bx")      v=static_cast<int>(BX);
  else if (name=="By")      v=static_cast<int>(BY);
  else if (name=="Bz")      v=static_cast<int>(BZ);
  else if (name=="psi")     v=static_cast<int>(SI);
  else if (name=="TR0")     {v=SimPM.ftr; /*cout <<"first tracer, v="<<v<<"\n";*/ }
  else if (name=="TR1")     v=SimPM.ftr+1;
  else if (name=="TR2")     v=SimPM.ftr+2;
  else if (name=="TR3")     v=SimPM.ftr+3;
  else if (name=="TR4")     v=SimPM.ftr+4;
  else if (name=="Eint")     v=-1;
  else if (name=="divB")     v=-2;
  else if (name=="Ptot")     v=-3;
  else if (name.find("Col_Src") !=string::npos) v=-4;
  else rep.error("Bad variable index in fits write routine",name);

  long int ct=0;
  cell *c = gp->FirstPt();
  if (v>=0) {
    do {(*data)[ct] = c->P[v]; ct++;} while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-1) { // internal energy (or temperature if we have microphysics)
    do {
      if (MP) {
	(*data)[ct] = MP->Temperature(c->P,SimPM.gamma);
	//cout <<"temp="<<(*data)[ct]<<"\n";
      }
      else    (*data)[ct] = eqn->eint(c->P,SimPM.gamma);
//      cout <<"(*data) ["<<ct<<"] = "<<(*data)[ct] <<"\n";
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
  }

  else if (v==-2) { // divB
    int *vars=0;
    vars = mem.myalloc(vars,3);
    vars[0] = static_cast<int>(BX);
    vars[1] = static_cast<int>(BY);
    vars[2] = static_cast<int>(BZ);
    do {(*data)[ct] = eqn->Divergence(c,0,vars, gp); ct++;}
    while ( (c=gp->NextPt(c))!=0 );
    vars = mem.myfree(vars);
  }

  else if (v==-3) { // total pressure.
    do {
      (*data)[ct] = eqn->Ptot(c->P,SimPM.gamma);
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
  }

  else if (v==-4) {
    //
    // optical depth variable: parse source id, and which Tau var we
    // want to output from the variable name.
    //
    int Si = atoi(name.substr(8).c_str());
    int Ti = atoi(name.substr(11).c_str());
    double Tau[MAX_TAU];    
    do {
      CI.get_col(c, Si, Tau);      
      (*data)[ct] = Tau[Ti];
      ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }

  else rep.error("Don't understand what variable to write.",v);

  if (ct!=ntot) rep.error("Counting cells in put_variable_into_data()",ct-ntot);
  return 0;
}



// ##################################################################
// ##################################################################



int DataIOFits::read_fits_image(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      fitsfile *ff,
      string name,
      double *localxmin,
      double *globalxmin,
      int *npt,
      long int ntot
      )
{
  double *data=0;
  data = mem.myalloc(data,ntot);

  int err = utility_fitsio::read_fits_image_to_data(ff, name,
      SimPM.ndim, localxmin, globalxmin, gp->DX(), npt, ntot,
      TDOUBLE, static_cast<void *>(data));
  if (err) rep.error(" DataIOFits::read_fits_image() Failed to read image from file",err);

  // Choose variable to read to, based on name string, whose hdu is the currently open hdu.
  int v = 0;
  if      (name=="GasDens") v=static_cast<int>(RO);
  else if (name=="GasPres") v=static_cast<int>(PG);
  else if (name=="GasVX")   v=static_cast<int>(VX);
  else if (name=="GasVY")   v=static_cast<int>(VY);
  else if (name=="GasVZ")   v=static_cast<int>(VZ);
  else if (name=="Bx")      v=static_cast<int>(BX);
  else if (name=="By")      v=static_cast<int>(BY);
  else if (name=="Bz")      v=static_cast<int>(BZ);
  else if (name=="psi")     v=static_cast<int>(SI);
  else if (name=="TR0")     {v=SimPM.ftr; /*cout <<"first tracer, v="<<v<<"\n";*/ }
  else if (name=="TR1")     v=SimPM.ftr+1;
  else if (name=="TR2")     v=SimPM.ftr+2;
  else if (name=="TR3")     v=SimPM.ftr+3;
  else if (name=="TR4")     v=SimPM.ftr+4;
//  else if (name=="Eint")     v=-1;
//  else if (name=="divB")     v=-2;
//  else if (name=="Ptot")     v=-3;
  else rep.error("Bad variable index in fits write routine",name);
  if (v>=SimPM.nvar) rep.error("reading variable, but no element in vector free for it.",v-SimPM.nvar);
     
  // assign data to grid points, to the variable determined above.
  double nulval = -1.e99;
  long int ct=0; cell *c = gp->FirstPt(); 
  do {
    c->P[v] = data[ct];
    if (pconst.equalD(data[ct],nulval)) {
      cout <<"(dataio::read data: ERROR: var="<<v<<" and data="<<data[ct]<<"\n";
      rep.error("Read null value from file! pixel count follows",ct);
    }
    ct++;
  } while ( (c=gp->NextPt(c))!=0 );
  if (ct!=ntot) rep.error("Counting cells in read_fits_image()",ct-ntot);
  
  data = mem.myfree(data);

  return 0;
}



// ##################################################################
// ##################################################################





/*
  // From icgen.cc
  // contains code for writing a fits binary table... may be useful for AMR.
  int outputFitsData(string opfile, int flag)
  {
  cout <<"Writing data to file "<<opfile<<"\n";
  int ndim=gp->Ndim(); int nvar=gp->Nvar();
  fitsfile *ff;
  int status=0,err=0;
  //  cout <<"opening fits file...";
  // This creates a new file with name given.  status is non-zero if it fails,
  // which happens if the file exists.
  fits_create_file(&ff, opfile.c_str(), &status);
  if(status) {cerr<<"file open went bad.\n";exit(1);}
  //  cout <<"done.\n";
  
  char temp[256];
  cout <<"Writing Fits header.\n";
  //  int count=0;
  // Write header with all my parameters to the first HDU of the file.
  fits_create_img(ff,DOUBLE_IMG,0,0,&status);
  err += fits_write_key(ff, TINT, "gridtype",  &SimPM.gridType, "UniformGrid=1", &status);
  err += fits_write_key(ff, TINT, "eqntype",   &SimPM.eqntype,  "HD=1,MHD=2", &status);
  err += fits_write_key(ff, TINT, "solver",    &SimPM.solverType, "LF=1,GD=2", &status);
  err += fits_write_key(ff, TINT, "gridndim",  &SimPM.ndim, "Num. spatial dimensions of grid", &status);
  err += fits_write_key(ff, TINT, "eqnndim",   &SimPM.eqnNDim,  "Num. spatial dimensions in vectors", &status);
  err += fits_write_key(ff, TINT, "eqnnvar",   &SimPM.nvar,  "Num. variables in equations", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TINT, "ntracer",   &SimPM.ntracer, "Num. tracer variables in equations", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TINT, "coordsys",   &SimPM.coord_sys, "Coordinate system, cart=1,cyl=2", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "tsim",   &SimPM.simtime, "current time in simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "tstart", &SimPM.starttime, "initial time to start simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "tfinish",&SimPM.finishtime, "Time to finish the simulation", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TINT, "tstep",&SimPM.timestep, "TIMESTEP SINCE START OF SIM", &status);
  err += fits_write_key(ff, TINT, "ngx",&SimPM.NG[XX], "Number of grid-points in X-direction", &status);
  err += fits_write_key(ff, TINT, "ngy",&SimPM.NG[YY], "Number of grid-points in Y-direction", &status);
  err += fits_write_key(ff, TINT, "ngz",&SimPM.NG[ZZ], "Number of grid-points in Z-direction", &status);
  err += fits_write_key(ff, TINT, "Ncell",&SimPM.Ncell, "Total Number of Grid Points", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "xmin",   &SimPM.Xmin[XX], "min. value of X", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"xmin count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "ymin",   &SimPM.Xmin[YY], "min. value of Y", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"ymin count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "zmin",   &SimPM.Xmin[ZZ], "min. value of Z", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"zmin count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "xmax",   &SimPM.Xmax[XX], "max. value of X", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "ymax",   &SimPM.Xmax[YY], "max. value of Y", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  err += fits_write_key(ff, TDOUBLE, "zmax",   &SimPM.Xmax[ZZ], "max. value of Z", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  strcpy(temp,SimPM.typeofbc.c_str());
  err += fits_write_key(ff, TSTRING, "typeofbc", temp, "Type Of Boundary Conditions", &status);
  err += fits_write_key(ff, TINT, "nbc", &SimPM.Nbc, "Depth of boundary cells (redundant with spOOA)", &status);
  int spooa = static_cast<int>(SimPM.spOOA);
  int tmooa = static_cast<int>(SimPM.tmOOA);
  err += fits_write_key(ff, TINT, "spooa", &spooa, "Spatial Order of Accuracy", &status);
  err += fits_write_key(ff, TINT, "tmooa", &tmooa, "Time Order of Accuracy", &status);
  err += fits_write_key(ff, TDOUBLE, "gamma", &SimPM.gamma, "Gas Constant Gamma", &status);
  err += fits_write_key(ff, TDOUBLE, "cfl", &SimPM.CFL, "CFL number", &status);
  err += fits_write_key(ff, TINT, "artvisc", &SimPM.artviscosity, "Type of Art.Visc. 0=none", &status);
  err += fits_write_key(ff, TDOUBLE, "etav", &SimPM.etav, "Artificial Viscosity Parameter, eta", &status);
  err += fits_write_key(ff, TINT, "typeofop", &SimPM.typeofop, "Type of Output File to write", &status);
  strcpy(temp,SimPM.outFileBase.c_str());
  err += fits_write_key(ff, TSTRING, "outfile", temp, "Base filename, with path, for data outputs", &status);
  err += fits_write_key(ff, TINT, "opfreq", &SimPM.opfreq, "File Output Frequency", &status);
  err += fits_write_key(ff, TINT, "jetic",  &JP.jetic    , "Is a jet sim? 0=no,1=yes", &status);
  err += fits_write_key(ff, TINT, "jetrad", &JP.jetradius, "Radius of jet (no. of cells)", &status);
  for (int v=0;v<MAX_NVAR; v++) {
    ostringstream temp;
    temp << "jstate" << v;
    char t2[32]; strcpy(t2,(temp.str()).c_str());
    err += fits_write_key(ff, TDOUBLE, t2, &JP.jetstate[v],"Jet State Vector, element v", &status);
  }
  // Write Units.
  strcpy(temp,uc.unitsys.c_str());  err += fits_write_key(ff, TSTRING, "unitsys",  temp, "System of Units File is using", &status);
  strcpy(temp,uc.density.c_str());  err += fits_write_key(ff, TSTRING, "unitdens", temp, "Density   Units", &status);
  strcpy(temp,uc.length.c_str());   err += fits_write_key(ff, TSTRING, "unitlen",  temp, "Length    Units", &status);
  strcpy(temp,uc.velocity.c_str()); err += fits_write_key(ff, TSTRING, "unitvel",  temp, "Velocity  Units", &status);
  strcpy(temp,uc.bfield.c_str());   err += fits_write_key(ff, TSTRING, "unitmagf", temp, "Mag.Field Units", &status);
  err += fits_write_key(ff, TDOUBLE, "rhoval", &uc.rhoVal, "Number of Ref. units per code unit", &status);
  err += fits_write_key(ff, TDOUBLE, "lenval", &uc.lenVal, "Number of Ref. units per code unit", &status);
  err += fits_write_key(ff, TDOUBLE, "velval", &uc.velVal, "Number of Ref. units per code unit", &status);
  err += fits_write_key(ff, TDOUBLE, "magval", &uc.magVal, "Number of Ref. units per code unit", &status);
  cout <<"Header written moving on to data.\n";
  //  err += fits_write_key(ff, TDOUBLE, "", &SimPM., "", &status);
  //  if(status) {fits_report_error(stderr,status);cout <<"count="<<count<<"\n";exit(1);} else count++;
  
  if (flag==2) { // write data to fits images.
    // Now write all the data in sequence, with each variable written to its own HDU image.
    double *vv = new double [SimPM.Ncell]; int ct;
    long int pix[ndim];
    string *var=0;
    if (SimPM.ntracer>5) rep.error("OutputFitsData:: only handles up to 5 tracer variables! Add more if needed.",SimPM.ntracer);
    if (SimPM.eqntype==EQEUL) {
      var  = new string [10];
      string pvar[10] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","TR0","TR1","TR2","TR3","TR4"};
      for (int i=0;i<10;i++) var[i] = pvar[i];
    }
    else if (SimPM.eqntype==EQMHD) {
      var  = new string [13];
      string pvar[13] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","Bx","By","Bz","TR0","TR1","TR2","TR3","TR4"};
      for (int i=0;i<13;i++) var[i] = pvar[i];
    }
    else if (SimPM.eqntype==EQGLM) {
      var  = new string [14];
      string pvar[14] = {"GasDens","GasPres","GasVX","GasVY","GasVZ","Bx","By","Bz","psi","TR0","TR1","TR2","TR3","TR4"};
      for (int i=0;i<14;i++) var[i] = pvar[i];
    } 
    else {
      var = new string [10]; 
      rep.error("What equations?!",SimPM.eqntype);
    }
    for (int i=0;i<nvar;i++) {
      for (int j=0;j<ndim;j++) pix[j]=SimPM.NG[j];
      fits_create_img(ff,DOUBLE_IMG,ndim,pix,&status);
      char temp[32]; strcpy(temp,var[i].c_str());
      err += fits_write_key(ff, TSTRING, "extname",  temp, "Image Name", &status);
      //  err += fits_write_key(ff, TINT, "bzero",  &SimPM.Xmin[XX], "Xmin", &status);
      //  err += fits_write_key(ff, TINT, "bscale",  &SimPM.dx, "Size of pixel", &status);
      cell *cpt = gp->FirstPt();
      int v;
      if      (var[i]=="GasDens") v=static_cast<int>(RO);
      else if (var[i]=="GasPres") v=static_cast<int>(PG);
      else if (var[i]=="GasVX")   v=static_cast<int>(VX);
      else if (var[i]=="GasVY")   v=static_cast<int>(VY);
      else if (var[i]=="GasVZ")   v=static_cast<int>(VZ);
      else if (var[i]=="Bx")      v=static_cast<int>(BX);
      else if (var[i]=="By")      v=static_cast<int>(BY);
      else if (var[i]=="Bz")      v=static_cast<int>(BZ);
      else if (var[i]=="psi")     v=static_cast<int>(SI);
      else if (var[i]=="TR0")     v=SimPM.ftr;
      else if (var[i]=="TR1")     v=SimPM.ftr+1;
      else if (var[i]=="TR2")     v=SimPM.ftr+2;
      else if (var[i]=="TR3")     v=SimPM.ftr+3;
      else if (var[i]=="TR4")     v=SimPM.ftr+4;
      else rep.error("Bad variable index in fits write routine",var[i]);
      ct=0;
      do {
      //	if (isnan(cpt->P[v])) {
      //	  cout <<"NAN value in Primitive vector: ct="<<ct<<" and var="<<v<<" and P[v]="<<cpt->P[v]<<"\n";
      //	}
	vv[ct] = cpt->P[v];
	ct++;
      } while ( (cpt=gp->NextPt(cpt))!=0 );
      for (int i=0;i<ndim;i++) pix[i]=1;
      cout <<"Writing data for variable: "<<temp<<"\n";
      err += fits_write_pix(ff, TDOUBLE, pix, SimPM.Ncell, vv, &status);
      if(status) fits_report_error(stderr,status);
    } // Primitive Variables
    
  }// End of writing image.
  
  else if (flag==3){ // Write Binary Table
  //  int  fits_create_tbl(fitsfile ptr.to.file, ///< fits file pointer.
  //                int tbltype=bin/acs, ///< 
  //		  int num.reserved.rows,///< should be set to zero unless i have a good reason not to.
  //		  int num.cols, ///< number of columns of table.
  //		  char** col.name.arr, ///< array of c_strings with names of columns
  //		  char** col.datattype.arr, ///< array of c_strings with datatypes J=int,E=float,D=double
  //		  char** col.phys.units.arr, ///< array of phsyical units, can be null
  //		  char* table.name, ///< name of table, can be null.
  //		  int &status ///< error code for returning.
  //		  )
  int nels = 32;
  int ncols = 3;
  //    if (SimPM.eqntype ==2) ncols +=2; // For divB and Ptot
    string colnames[] = {"Posn","PrimVec"};
    string coltypes[] = {"1D","1D","1D","1D","1D"};
    char *ttype[ncols]; char *tform[ncols];
    for (int i=0; i<ncols; i++) {
      ttype[i] = new char [nels]; tform[i] = new char [nels];
    }
    for (int i=0; i<ncols; i++) {
      strcpy(ttype[i],colnames[i].c_str());
      strcpy(tform[i],coltypes[i].c_str());
    }
    char tblname[nels]; strcpy(tblname,"datatable");
    err += fits_create_tbl(ff, BINARY_TBL, 0, ncols,ttype,tform,0,tblname,&status);
    err += fits_modify_vector_len(ff,1,gp->Ndim(),&status); // change position column to ndim vector.
    err += fits_modify_vector_len(ff,2,gp->Nvar(),&status); // change primitive variable column to nvar vector.
    // Insert rows in the binary table, and write data to the rows.
    fits_insert_rows(ff,0,SimPM.Ncell,&status);
    cell *cpt = gp->FirstPt();  int row=0;
    do {
      row++; // fits uses unit offset.
      //fits_write_col(*fptr, datatype, colnum, firstrow, firstelem, nelements, data, &status)
      err += fits_write_col(ff, TDOUBLE, 1, row, 1, gp->Ndim(), cpt->x, &status);
      err += fits_write_col(ff, TDOUBLE, 2, row, 1, gp->Nvar(), cpt->P, &status);
     } while ( (cpt=gp->NextPt(cpt))!=0 );
  } // end write data to table.

  else rep.error("Bad output type passed to outputFitsData()",flag);
  
  err += fits_close_file(ff,&status);
  return(err);
  }
*/


#endif // if FITS
