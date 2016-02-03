/// \file twoD2oneD_shocktubes.cc
/// \author Jonathan Mackey
/// Written 2009-12-22
///
/// Updated code to convert 2D shocktube data to a 1D string of points which
/// can be plotted with gnuplot.
///
/// - 2010.12.07 JM: Added setup_extra_data() call to cell_interface
///   before setting up grid.  Need to add Geomtric grid options!
///
/// - 2010.12.31 JM: Updated setup_extra_data() function.  Updated equations
///   to allow EQEUL_EINT

#ifdef FITS
#include <fitsio.h>
#include "../../../source/dataIO/dataio_fits.h"
#endif
#ifdef SILO
#include <silo.h>
#include "../../../source/dataIO/dataio_silo.h"
#endif
#include "../../../source/global.h"
#include "../../../source/grid/uniform_grid.h"
#include "../../../source/dataIO/dataio.h"
#include <cmath>
#include <sstream>
using namespace std;



void rotateXY(const double *v_in, ///< x,y pair of data
	      const double theta, ///< rotation angle (radians)
	      double *v_out       ///< x,y pair of output data.
	      )
{
  //
  // This rotates an [x,y] pair via an angle +theta.  This rotates the
  // coordinate system, so it is equivalent to moving a vector through 
  // -theta in a given system.  It is the opposite of eqns_base::rotateXY()
  //
  double ct=cos(theta); double st=sin(theta);
  v_out[XX] =  v_in[XX]*ct + v_in[YY]*st;
  v_out[YY] = -v_in[XX]*st + v_in[YY]*ct;
  return;
}

bool cell_in_region(const double *cpos,
		    const double *xmin,
		    const double *xmax
		    )
{
  for (int v=0;v<SimPM.ndim;v++) {
    if (cpos[v] < xmin[v]) return false;
    if (cpos[v] > xmax[v]) return false;
  }
  return true;
}



int main(int argc, char **argv)
{
  //
  // Get input files and output files from cmd-line args.
  //
  if (argc!=9) {
    cerr << "Error: must call with 8 arguments...\n";
    cerr << "twoD2oneD_shocktubes: <executable> <Infile> <infile-type[fits/silo]> <rotation angle[int,degrees]> <outfile.txt> xmin xmax ymin ymax\n";
    rep.error("Bad number of Args",argc+1);
  }

  string infile(argv[1]);
  string ftype(argv[2]);
  string outfile(argv[4]);
  //
  // Get angle and convert to radians.
  //
  int ang = atoi(argv[3]);
  double theta = 0.0;
  if (ang ==-1) {
    //
    // Assume we are running a test with theta=atan(0.5) 
    //
    theta = atan(0.5);
  }
  else if (ang ==-2) {
    //
    // now it's atan(2)
    //
    theta = atan(2.0);
  }
  else {
    theta = static_cast<double>(ang)*M_PI/180.0;
  }

  //
  // Set up dataio class, either fits or silo.
  //
  class DataIOBase *dataio=0;
#ifdef FITS
  if (ftype=="fits" || ftype=="FITS") {
    dataio = new DataIOFits ();
    ftype  = "fits";
  }
#endif
#ifdef SILO
  if (ftype=="silo" || ftype=="SILO") {
    dataio = new dataio_silo ();
    ftype  = "silo";
  }
#endif
  if (!dataio) rep.error("Failed to init dataio class... check ifdefs.",ftype);

  
  //
  // First we need to open the first file in the list, and get the grid dimensions,
  // so we can set it up once and use it for all the infiles.
  //
  cout <<"**********************************************\n";
  cout <<"\tReading from input file "<<infile;
  cout <<" and writing new positions to file "<<outfile <<endl;
  cout <<"**********************************************\n";

  //ostringstream temp; temp.str("");
  //temp << infilebase <<"."<<startct <<".fits";
  //infile=temp.str();
  cout <<"Initially reading from file "<<infile<<endl;

  int err=0;
  if (!dataio->file_exists(infile)) rep.error("######################## File not found!",infile);
  err = dataio->ReadHeader(infile);
  if (err) rep.error("read header went bad",err);
  
  int eqns=SimPM.eqntype;
  if      (eqns==EQEUL ||
	   eqns==EQEUL_EINT)
    eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD)
    eqns=2;
  else rep.error("Bad equations",eqns);

  //
  // check dimensionality is ok.
  //
  if (SimPM.ndim!=2)
    rep.error("need 2D sim for shock tube test",SimPM.ndim);

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  int hc_flag = 0, dv_flag=0;
  if (SimPM.artviscosity==AV_LAPIDUS ||
      SimPM.eqntype==EQEUL_EINT) {
    // Need one var. for Div(v)
    dv_flag = 1;
  }
  if (SimPM.artviscosity==AV_HCORRECTION ||
      SimPM.artviscosity==AV_HCORR_FKJ98) {
    // need one var for each dimension here.
    hc_flag = SimPM.ndim;
  }
  CI.setup_extra_data(SimPM.RS, hc_flag, dv_flag);

  //
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  //
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (grid==0) rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout <<"(setup_grid) Done. g="<<grid<<" and dx="<<grid->DX()<<"\n";

  
  //
  // Set up and open text outfile
  //
  if (dataio->file_exists(outfile)) cout <<"WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",outfile);
  cout <<"writing to file "<<outfile<<endl;
  outf.setf( ios_base::scientific );
  outf.precision(6);
  outf << "# twoD2oneD_shocktubes.cc: input file: "<<infile<<" and input angle="<<theta*180.0/M_PI<<endl;
  outf << "# Columns are rotated x-pos, rho, p_g, v_x, v_y, v_z, [B_x, B_y, B_z, [Psi]].\n#\n";

  //
  // read data onto grid.
  //
  err += dataio->ReadData(infile,grid);
  if (err) rep.error("read data went bad for file",err);

  //
  // Choose the subdomain we want to take data from:
  // Set from the commandline.
  //
  double xmin[SimPM.ndim], xmax[SimPM.ndim];
  xmin[XX] = atof(argv[5]);
  xmax[XX] = atof(argv[6]);
  xmin[YY] = atof(argv[7]);
  xmax[YY] = atof(argv[8]);

  cell *c = grid->FirstPt();
  double cpos[SimPM.ndim];   // cell position
  double xprime[SimPM.ndim]; // rotated cell position.
  double vec[3];      // vector field
  double vecprime[3]; // rotated vector field.

  do {
    CI.get_dpos(c,cpos);
    //
    // If cell is in the central region, add it to the list, otherwise skip it.
    //
    if (cell_in_region(cpos,xmin,xmax)) {
      rotateXY(cpos,theta,xprime);
      //
      // Rotate the velocity field.
      //
      vec[0]=c->P[VX]; vec[1]=c->P[VY]; vec[2]=c->P[VZ];
      rotateXY(vec,theta,vecprime);
      c->P[VX] = vecprime[0]; c->P[VY] = vecprime[1]; //c->P[VZ] = vecprime[2];
      //
      // Rotate the Magnetic Field (if present)
      //
      if (eqns==2) {
	vec[0]=c->P[BX]; vec[1]=c->P[BY]; vec[2]=c->P[BZ];
	rotateXY(vec,theta,vecprime);
	c->P[BX] = vecprime[0]; c->P[BY] = vecprime[1]; //c->P[BZ] = vecprime[2];
      }
      //
      // Now output data
      //
      outf << xprime[XX] <<"\t"; //<< xprime[YY] <<"\t";
      for (int v=0;v<SimPM.nvar;v++) outf << c->P[v]<<"\t";
      outf <<endl;
    }
  } while ((c=grid->NextPt(c)) !=0);
    
  outf.close();

  delete grid; grid=0;
  delete dataio; dataio=0;
  return 0;
}
