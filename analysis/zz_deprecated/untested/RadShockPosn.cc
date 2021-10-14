/** \file RadShockPosn.cc
 * \author Jonathan Mackey
 *
 * This file reads in a series of fits files, which are assumed to be outputs
 * of radiative shock simulations in 1D or 2D.  For each file it determines
 * the position of the shock in the x-direction, in the 1st y-column, and
 * writes the (time,position) pair to a tab delimited text file.  The first line
 * of the file should contain info about the simulation that the data is from.
 *
 * This is explicitly serial code, so if there are parallel outputs to multiple
 * files this code won't work.
 *
 * Compile with:\n
 * g++ -Wall -DSERIAL RadShockPosn.cc ../testing/global.cc
 * ../testing/uniformGrid.cc ../testing/dataio.cc -lreadline -lcfitsio
 *
 * Run with (example):\n
 * ./a.out v140.txt ../results/RadShock2D_n128x32_v140_n10_T1e4 0 50
 *
 * */

#include "fitsio.h"
using namespace std;
#include "../testing/dataio.h"
#include "../testing/global.h"
#include "../testing/uniformGrid.h"

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc != 5) {
    cerr << "Error: must call with 4 arguments...\n";
    cerr
        << "RadShockPosn: <executable> <OutfileName> <file-base> <FirstOutput> <OutputFreq>\n";
    rep.error("Bad number of Args", argc);
  }
  string outfile    = argv[1];
  string infilebase = argv[2];
  int startct       = atoi(argv[3]);
  int opfreq        = atoi(argv[4]);
  if (isnan(startct) || isnan(opfreq) || opfreq == 0)
    rep.error("Bad ints in args", opfreq);

  cout << "reading from first file " << infilebase << "." << startct
       << ".fits\n";
  cout << "Writing shock position to file " << outfile << "\n";
  cout << "**********************************************\n";

  class DataIOFits dataio;
  class file_status fs;

  // First we need to open the first file in the list, and get the grid
  // dimensions, so we can set it up once and use it for all the infiles.
  string infile;
  ostringstream temp;
  temp.str("");
  temp << infilebase << "." << startct << ".fits";
  infile = temp.str();
  cout << "Initially reading from file " << infile << endl;

  int err = 0;
  if (!fs.file_exists(infile)) rep.error("First file not found!", infile);
  err = dataio.ReadHeader(infile);
  if (err) rep.error("read header went bad", err);

  // check dimensionality is ok.
  if (SimPM.ndim != 1 && SimPM.ndim != 2)
    rep.error("need 1D or 2D sim for rad.shock test", SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout << "(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid(
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (grid == 0)
    rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout << "(setup_grid) Done. g=" << grid << "\n";

  // Set up and open outfile
  if (fs.file_exists(outfile))
    cout << "WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if (!outf.is_open()) rep.error("couldn't open outfile", outfile);
  outf.setf(ios_base::scientific);
  outf.precision(6);
  outf << "# Radiative Shock Test Problem outputs.  First file: " << infile
       << endl;
  outf
      << "# Columns are time and shock position, and should be in cgs units (s,cm).\n\n";

  int count     = startct;
  double refvel = 0.0;
  // Need to loop this over all timesteps, incrementing 'start' by 'step' each
  // time until there are no more files to analyse (note the last one might get
  // left out).
  do {
    cout << "Reading from file " << infile << endl;
    // read data onto grid.
    err += dataio.ReadHeader(infile);
    err += dataio.ReadData(infile);
    if (err) rep.error("read data went bad for file", err);
    // get first point, and move to XP end of grid.
    cell *c = grid->FirstPt();
    do {
      c = grid->NextPt(c, XP);
    } while (grid->NextPt(c, XP) != 0);
    cell *c2 = grid->NextPt(c, XN);
    if (!c2) {
      rep.error("Lost on grid", c2);
      grid->PrintCell(c);
    }
    refvel = c->P[VX];
    // find the shock position by locating where VX first changes by >10%
    while (fabs(fabs(c2->P[VX] / refvel) - 1.0) <= 0.3) {
      c  = c2;
      c2 = grid->NextPt(c2, XN);
      if (!c2) {
        cout << "no shock found!\n";
        c2 = c;
        break;
      }
    }
    // Write (x_sh,t_sim) to file.
    outf << SimPM.simtime << "\t" << c2->x[XX] << "\n";
    // increment filename
    count += opfreq;
    temp.str("");
    temp << infilebase << "." << count << ".fits";
    infile = temp.str();
  } while (fs.file_exists(infile));
  // loop over all timesteps.
  cout << "\n***************************************************\n";
  cout << "couldn't find file " << infile << " for step " << count
       << "... assuming i'm finished!\n";
  outf.close();
  delete grid;
  grid = 0;
  return 0;
}
