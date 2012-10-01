/// \file RadShockPosn.cc
/// \author Jonathan Mackey
/// 
/// This file reads in a series of files, which are assumed to be outputs
/// of radiative shock simulations in 1D or 2D.  For each file it determines
/// the position of the shock in the x-direction, in the 1st y-column, and 
/// writes the (time,position) pair to a tab delimited text file.  The first line
/// of the file should contain info about the simulation that the data is from.
/// 
/// This is explicitly serial code, so if there are parallel outputs to multiple
/// files this code won't work.
/// 
/// Compile with:\n
/// make -f Makefile.RadShockPosn
/// 
/// Run with (example):\n
/// ./RadShockPosn v140.txt ../results/RadShock2D_n128x32_v140_n10_T1e4 0 50
/// 
///
/// 2009-12-17 JM:
///  Modified it to work with new code structure, and with silo files.  Seems
///  to have been last changed in March 2009, but that might have been the svn 
///  import, so maybe more than a year since I used it.
/// - 2011.04.14 JM: Modified to work with new code. 

#ifdef FITS
#include <fitsio.h>
#include "../../source/dataIO/dataio_fits.h"
#endif
#ifdef SILO
#include <silo.h>
#include "../../source/dataIO/dataio_silo.h"
#endif
#include "../../source/global.h"
#include "../../source/grid/uniform_grid.h"
#include "../../source/dataIO/dataio.h"
#include <cmath>
#include <sstream>
using namespace std;

void get_file_name(const string , ///< file base
		   const string , ///< file type (fits,silo)
		   const int ,    ///< count.
		   string &       ///< write to filename.
		   );

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc!=6) {
    cerr << "Error: must call with 5 arguments...\n";
    cerr << "RadShockPosn: <executable> <OutfileName> <infile-base> <infile-type[fits,silo]> <FirstOutput> <OutputFreq>\n";
    rep.error("Bad number of Args",argc);
  }
  string outfile = argv[1];
  string infilebase  = argv[2];
  int startct = atoi(argv[4]);
  int opfreq  = atoi(argv[5]);
  if (isnan(startct) || isnan(opfreq) || opfreq==0) rep.error("Bad ints in args",opfreq);
  
  cout <<"reading from first file "<<infilebase<<"."<<startct<<"."<<argv[3]<<"\n";
  cout <<"Writing shock position to file "<<outfile<<"\n";
  cout <<"**********************************************\n";

  string ftype(argv[3]);
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

  
  // First we need to open the first file in the list, and get the grid dimensions,
  // so we can set it up once and use it for all the infiles.
  string infile;
  get_file_name(infilebase, ftype, startct, infile);
  
  int err=0;
  if (!dataio->file_exists(infile)) rep.error("First file not found!",infile);
  err = dataio->ReadHeader(infile);
  if (err) rep.error("read header went bad",err);
  
  // check dimensionality is ok.
  if (SimPM.ndim!=1 && SimPM.ndim!=2)
    rep.error("need 1D or 2D sim for rad.shock test",SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
  CI.setup_extra_data(SimPM.RS,  0, 0);
  grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (!grid) rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout <<"(setup_grid) Done. g="<<grid<<"\n";

  // Set up and open outfile
  if (dataio->file_exists(outfile)) cout <<"WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",outfile);
  outf.setf( ios_base::scientific );
  outf.precision(6);
  outf << "# Radiative Shock Test Problem outputs.  First file: "<<infile<<endl;
  outf << "# Columns are time and shock position, and should be in cgs units (s,cm).\n\n";
  
  int count = startct; double refvel=0.0;
  // Need to loop this over all timesteps, incrementing 'start' by 'step' each time
  // until there are no more files to analyse (note the last one might get left out).
  do {
     cout <<"Reading from file "<<infile<<endl;
     // read data onto grid.
     err += dataio->ReadHeader(infile);
     err += dataio->ReadData(infile,grid);
     if (err) rep.error("read data went bad for file",err);

     //
     // get first point, and move to XP end of grid, and get the 
     // inflow velocity.
     //
     cell *c = grid->FirstPt();
     do {c=grid->NextPt(c,XP);} while (grid->NextPt(c,XP) !=0);
     cell *c2 = grid->NextPt(c,XN); if (!c2) {rep.error("Lost on grid",c2); grid->PrintCell(c);}
     refvel = c->P[VX];

     //
     // find the shock position by locating where VX first changes by >30%
     //
     while ( fabs(fabs(c2->P[VX]/refvel)-1.0) <= 0.3) {
       c = c2;
       c2 = grid->NextPt(c2,XN);
       if (!c2) { cout <<"no shock found!\n"; c2=c; break; }
     }

     //
     // Write (x_sh,t_sim) to file.
     //
     //outf <<SimPM.simtime<<"\t"<<c2->x[XX]<<"\n";
     outf <<SimPM.simtime<<"\t"<<CI.get_dpos(c2,static_cast<int>(XX))<<"\n";

     //
     // increment filename
     //
     count += opfreq;
     get_file_name(infilebase, ftype, count, infile);

  } while (dataio->file_exists(infile));
  // loop over all timesteps.

  cout <<"\n***************************************************\n";
  cout <<"couldn't find file "<<infile<<" for step "<<count<<"... assuming i'm finished!\n";
  outf.close();

  delete dataio; dataio=0;
  delete grid; grid=0;

  return 0;
}

  
void get_file_name(const string base, ///< file base
		   const string type, ///< file type (fits,silo)
		   const int ct,      ///< count.
		   string &fname      ///< write to filename.
		   )
{
    ostringstream temp; temp.str("");
    temp << base;
#ifdef SERIAL
    temp <<".";
#endif
#ifdef PARALLEL
    temp <<"_0000.";
#endif

    if (ct >= 0) {
      //
      // Only write the count if it's positive.
      // Silo files also have a six digit count.
      //
      if (type=="silo") {
	temp.width(8); temp.fill('0');
      }
      temp << ct << ".";
    }

    temp << type;
    fname = temp.str();
    return;
}
