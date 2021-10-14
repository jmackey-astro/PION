/** \file parallelionisation.cc
 * \author Jonathan Mackey
 *
 * This file reads in a series of fits files, which are assumed to be outputs
 * of a photoevaporation problem with a single ionising source at infinity in
 * the XN direction.
 *
 * This is explicitly serial code, so if there are parallel outputs to multiple
 * files this code won't work.
 *
 * Compile with:\n
 * g++ -Wall -g -DSERIAL parallelionisation.cc ../../testing/global.cc
 * ../../testing/uniformGrid.cc ../../testing/dataio.cc -lreadline -lcfitsio
 *
 * Run with (example):\n
 * ./a.out positions ../results/pec3D_n51_recomb 0 10 6
 * plot_radius: <executable> <Outfile-base> <file-base> <FirstO/P-Num>
 * <O/P-Freq> <x(HII)index>
 * */

#include "fitsio.h"
using namespace std;
#include "../../testing/dataio.h"
#include "../../testing/global.h"
#include "../../testing/uniformGrid.h"

// Given a source and a sim-time, we can calculate the position of the I-front.
int get_theoretical_position(
    const double,  ///< density in this column
    double *       ///< variable to put radius in.
);


// Find the best-fitting position for the I-front without any prior assumptions.
int fit_data(
    long int,  // npt
    double *,  // array of distances from edge of grid
    double *,  // array of ion-fractions.
    double,    // a_min (min value of fit parameter)
    double,    // a_max (max value of fit parameter)
    double *,  // fit value of a, returned.
    double *,  // goodness of fit, returned (reduced chi-squared)
    double *   // optional second goodness of fit parameter
);

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc != 6) {
    cerr << "Error: must call with 5 arguments...\n";
    cerr
        << "parallelionisation: <executable> <Outfile-base> <Infile-base> <First-File-Num> <File-Num-Freq> <x(HII) var>\n";
    rep.error("Bad number of Args", argc);
  }
  string outfilebase = argv[1];
  string infilebase  = argv[2];
  int startct        = atoi(argv[3]);
  int opfreq         = atoi(argv[4]);
  int var            = atoi(argv[5]);
  if (isnan(startct) || isnan(opfreq) || opfreq == 0)
    rep.error("Bad ints in args", opfreq);

  cout << "\n\nreading from first file " << infilebase << "." << startct
       << ".fits\n";
  cout << "**********************************************\n";

  class DataIOFits dataio;
  class file_status fs;

  // First we need to open the first file in the list, and get the grid
  // dimensions, so we can set it up once and use it for all the infiles.
  string infile, outfile, outfitsfile;
  ostringstream temp;
  temp.str("");
  temp << infilebase << "." << startct << ".fits";
  infile = temp.str();
  cout << "Initially reading from file " << infile << endl;

  int err = 0;
  if (!fs.file_exists(infile))
    rep.error("########################First file not found!", infile);
  err = dataio.ReadHeader(infile);
  if (err) rep.error("read header went bad", err);

  // check dimensionality is ok.
  if (SimPM.ndim != 1 && SimPM.ndim != 2 && SimPM.ndim != 3)
    rep.error("need 1D/2D/3D sim for parallelrays test", SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout << "(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid(
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (grid == 0)
    rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout << "(setup_grid) Done. g=" << grid << " and dx=" << grid->DX() << "\n";
  err += dataio.ReadData(infile);
  if (err) rep.error("Read initial data", err);

  int count = startct;
  // Need to loop this over all timesteps, incrementing 'start' by 'step' each
  // time until there are no more files to analyse (note the last one might get
  // left out).

  if (RSP.nsources != 1) rep.error("need exactly one source...", RSP.nsources);

  double srcpos[SimPM.ndim];
  // ofstream of2("time.txt");
  for (int i = 0; i < SimPM.ndim; i++)
    srcpos[i] = RSP.position[0][i];

  // Set up and open text outfile
  temp.str("");
  temp << outfilebase << ".txt";
  outfile = temp.str();
  if (fs.file_exists(outfile))
    cout << "WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if (!outf.is_open()) rep.error("couldn't open outfile", outfile);
  cout << "writing to file " << outfile << endl;
  outf.setf(ios_base::scientific);
  outf.precision(2);
  outf << "# Parallelrays:  First file: " << infile << endl;
  outf
      << "# Columns are [time, I-front positions/calc.posn.], and should be in cgs units (s,cm).\n\n";

  cell *c = grid->FirstPt();
  outf << "# column number densities follow.\n# ";

  if (SimPM.ndim == 2) {
    do {
      outf << c->P[RO] / GS.m_p() << " ";
    } while ((c = grid->NextPt(c, YP)) != 0);
    outf << "\n\n";
  }
  else if (SimPM.ndim == 1)
    outf << c->P[RO] / GS.m_p() << "\n\n";

  outf.precision(6);

  double *dist = 0, *ifrac = 0;
  dist  = new double[SimPM.NG[XX]];
  ifrac = new double[SimPM.NG[XX]];
  if (!dist || !ifrac) rep.error("mem. alloc. dist/ifrac", dist);

  // increment filename
  count += opfreq;
  temp.str("");
  temp << infilebase << "." << count << ".fits";
  infile = temp.str();
  if (!fs.file_exists(infile))
    rep.error("Can't find second file to open", infile);

  // loop over all files.
  do {
    cout << "Reading from file " << infile << endl;
    // read data onto grid.
    err += dataio.ReadHeader(infile);
    err += dataio.ReadData(infile);
    if (err) rep.error("read data went bad for file", err);

    if (SimPM.ndim == 2) {
      // get first point, and for each point calculate distance to source (at
      // origin), and also the number of ionised atoms.
      outf << SimPM.simtime << "  ";
      cell *startpt = grid->FirstPt();
      do {
        int i = 0;
        c     = startpt;
        do {
          dist[i]  = c->x[XX] - SimPM.Xmin[XX];
          ifrac[i] = c->P[var];
          i++;
        } while ((c = grid->NextPt(c, XP)) != 0);
        // if (SimPM.simtime==0.0 && grid->NextPt(startpt,YN)==0)
        // for (int v=0;v<SimPM.NG[XX];v++) cout <<"  "<<dist[v]<<endl;
        // cout <<"got "<<i<<"cells\n";
        if (i > SimPM.NG[XX]) rep.error("too many cells in x-dir!", i);
        // now fit a position to the data
        // pass in (npt, dist, ifrac, a_min, a_max, &a_fit, &rchisq, &err2)
        double afit = 0.0, acalc = 0.0, rchisq = 0.0, rchisq2 = 0.0;
        err += fit_data(
            SimPM.NG[XX], dist, ifrac, grid->DX() / 2., grid->Range(XX), &afit,
            &rchisq, &rchisq2);
        // cout <<"afit="<<afit<<endl;
        err += get_theoretical_position(startpt->P[RO] / GS.m_p(), &acalc);
        outf << afit << "  " << acalc << "\t";
      } while ((startpt = grid->NextPt(startpt, YP)) != 0);
      outf << "\n";
    }  // if 2D
    else if (SimPM.ndim == 1) {
      int i = 0;
      outf << SimPM.simtime << "  ";
      c = grid->FirstPt();
      do {
        dist[i]  = c->x[XX] - SimPM.Xmin[XX];
        ifrac[i] = c->P[var];
        i++;
      } while ((c = grid->NextPt(c, XP)) != 0);
      if (i > SimPM.NG[XX]) rep.error("too many cells in x-dir!", i);
      // now fit a position to the data
      // pass in (npt, dist, ifrac, a_min, a_max, &a_fit, &rchisq, &err2)
      double afit = 0.0, acalc = 0.0, rchisq = 0.0, rchisq2 = 0.0;
      err += fit_data(
          SimPM.NG[XX], dist, ifrac, grid->DX() / 2., grid->Range(XX), &afit,
          &rchisq, &rchisq2);
      // cout <<"afit="<<afit<<endl;
      err +=
          get_theoretical_position(grid->FirstPt()->P[RO] / GS.m_p(), &acalc);
      outf << afit << "  " << acalc << "  ";
      afit = 0.0;
      for (int v = 0; v < i; v++)
        afit += ifrac[v];
      outf << afit * grid->DX() << "\n";
    }  // if 1D

    // increment filename
    count += opfreq;
    temp.str("");
    temp << infilebase << "." << count << ".fits";
    infile = temp.str();
  } while (fs.file_exists(infile));
  // loop over all timesteps.


  delete[] dist;
  dist = 0;
  delete[] ifrac;
  ifrac = 0;
  cout << "\n***************************************************\n";
  cout << "couldn't find file " << infile << " for step " << count
       << "... assuming i'm finished!\n";
  outf.close();
  delete grid;
  grid = 0;
  return 0;
}


// Given a source, a geometry, and a sim-time, we can calculate the position of
// the I-front.
int get_theoretical_position(
    const double nh,  ///< number density of gas in this column
    double *a         ///< variable to put position in.
)
{
  double rad = 0.0;

  // with recombs
  if (SimPM.EP.rad_recombination) {
    double alpha = 2.59e-13;  // recombination coeff.
    rad          = RSP.strength[0] / nh / nh / alpha;
    // cout <<"rad="<<rad<<" \texp="<<nh*alpha*SimPM.simtime<<endl;
    rad *= (1.0 - exp(-nh * alpha * SimPM.simtime));
    *a = rad;
  }
  else {
    // assume constant density!
    //  cout <<"calculating position!\n";
    rad = RSP.strength[0] * SimPM.simtime / nh;
    *a  = rad;
    if (*a > SimPM.Range[XX]) *a = SimPM.Range[XX];
  }
  // cout <<"got position = "<<*a<<endl;
  return 0;
}

int fit_data(
    long int n,     // npt
    double *x,      // array of radii from centre
    double *f,      // array of ion-fractions.
    double a_min,   // a_min (min value of fit parameter)
    double a_max,   // a_max (max value of fit parameter)
    double *a_fit,  // fit value of a, returned.
    double *rcsq,   // goodness of fit, returned
    double *rcsq2   // optional second goodness of fit parameter
)
{
  // take average position of cells with 0.1<x<0.9
  long int ct = 0;
  double rad = 0.0, max_rad = -1.e199, min_rad = 1.e199;
  for (int i = 0; i < n; i++) {
    if (f[i] >= 0.1 && f[i] <= 0.9) {
      rad += x[i];
      ct++;
      max_rad = max(max_rad, x[i]);
      min_rad = min(min_rad, x[i]);
    }
  }
  if (ct > 0) {
    *a_fit = rad / static_cast<double>(ct);
    *rcsq  = max_rad;
    *rcsq2 = min_rad;
  }
  else {
    if (f[0] < 0.1) {
      *a_fit = 0.0;
      *rcsq  = 0.0;
      *rcsq2 = 0.0;
    }
    else if (f[n - 1] > 0.9) {
      *a_fit = x[n - 1] + 0.5 * grid->DX();
      *rcsq  = x[n - 1];
      *rcsq2 = x[n - 1];
    }
    else {
      // there are no transition cells...
      bool found = false;
      for (int i = 0; i < n; i++) {
        if (f[i] < 0.5 && found == false) {
          found  = true;
          *a_fit = 0.5 * (x[i] + x[i - 1]);
          *rcsq  = 0.5 * (x[i] - x[i - 1]);
          *rcsq2 = 0.5 * (x[i + 1] - x[i]);
        }
      }
    }
  }
  //  cout <<"\tFound "<<ct<<" cells near I-front, with mean position "<<*a_fit;
  //  cout <<" and width "<<*rcsq-*rcsq2<<"\n";
  return 0;

  cout << "********* trying new method for weighted average.************\n";
  double cct = 0.0, wt = 0.0;
  rad     = 0.0;
  max_rad = -1.e199;
  min_rad = 1.e199;
  for (int i = 0; i < n; i++) {
    if (f[i] >= 0.1 && f[i] <= 0.9) {
      wt = pow(1.0 - 2.0 * fabs(0.5 - f[i]), 2.0);
      rad += x[i] * wt;
      cct += wt;
      max_rad = max(max_rad, x[i]);
      min_rad = min(min_rad, x[i]);
      // cout <<"f="<<f[i]<<" and wt="<<wt<<endl;
    }
  }
  if (cct > 0.01) {
    *a_fit = rad / cct;
    *rcsq  = max_rad;
    *rcsq2 = min_rad;
  }
  else {
    *a_fit = 0.0;
    *rcsq  = 0.0;
    *rcsq2 = 0.0;
  }
  cout << "\tFound " << ct << " cells near I-front, with mean position "
       << *a_fit;
  cout << " and width " << *rcsq - *rcsq2 << "\n";
}
