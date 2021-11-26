/** \file stromgensphere.cc
 * \author Jonathan Mackey
 *
 * This file reads in a series of fits files, which are assumed to be outputs
 * of a photoevaporation problem with a source at the origin.
 *
 * This is explicitly serial code, so if there are parallel outputs to multiple
 * files this code won't work.
 *
 * Compile with:\n
 * g++ -Wall -DSERIAL stromgensphere.cc ../testing/global.cc
 * ../testing/uniformGrid.cc ../testing/dataio.cc -lreadline -lcfitsio
 *
 * Run with (example):\n
 * ./a.out positions ../results/pec3D_n51_recomb 0 10 6
 * stromgensphere: <executable> <Outfile-base> <file-base> <FirstO/P-Num>
 * <O/P-Freq> <x(HII)index>
 * */

#include "fitsio.h"
using namespace std;
#include "../testing/dataio.h"
#include "../testing/global.h"
#include "../testing/uniformGrid.h"


int fit_data(
    long int,  // npt
    double *,  // array of radii from centre
    double *,  // array of ion-fractions.
    double,    // a_min (min value of fit parameter)
    double,    // a_max (max value of fit parameter)
    double *,  // fit value of a, returned.
    double *   // goodness of fit, returned (reduced chi-squared)
);

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc != 6) {
    cerr << "Error: must call with 5 arguments...\n";
    cerr
        << "stromgensphere: <executable> <OutfileName> <file-base> <FirstOutput> <OutputFreq> <x(HII) var>\n";
    spdlog::error("{}: {}", "Bad number of Args", argc);
  }
  string outfilebase = argv[1];
  string infilebase  = argv[2];
  int startct        = atoi(argv[3]);
  int opfreq         = atoi(argv[4]);
  int var            = atoi(argv[5]);
  if (isnan(startct) || isnan(opfreq) || opfreq == 0)
    spdlog::error("{}: {}", "Bad ints in args", opfreq);

  cout << "reading from first file " << infilebase << "." << startct
       << ".fits\n";
  cout << "Writing I-front position to file " << outfilebase << "\n";
  cout << "**********************************************\n";

  class DataIOFits dataio;
  class file_status fs;

  // First we need to open the first file in the list, and get the grid
  // dimensions, so we can set it up once and use it for all the infiles.
  string infile, outfile;
  ostringstream temp;
  temp.str("");
  temp << infilebase << "." << startct << ".fits";
  infile = temp.str();
  cout << "Initially reading from file " << infile << endl;

  int err = 0;
  if (!fs.file_exists(infile))
    spdlog::error("{}: {}", "First file not found!", infile);
  err = dataio.ReadHeader(infile);
  if (err) spdlog::error("{}: {}", "read header went bad", err);

  // check dimensionality is ok.
  if (SimPM.ndim != 1 && SimPM.ndim != 2 & SimPM.ndim != 3)
    spdlog::error("{}: {}", "need 1D/2D/3D sim for stromgen test", SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout << "(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid(
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin.data(), SimPM.Xmax,
      SimPM.NG);
  if (grid == 0)
    spdlog::error(
        "{}: {}", "(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout << "(setup_grid) Done. g=" << grid << "\n";


  int count = startct;
  // Need to loop this over all timesteps, incrementing 'start' by 'step' each
  // time until there are no more files to analyse (note the last one might get
  // left out).
  double srcpos[SimPM.ndim];
  // ofstream of2("time.txt");
  for (int i = 0; i < SimPM.ndim; i++)
    srcpos[i] = 0.0;
  // Set up and open outfile
  temp.str("");
  temp << outfilebase << ".txt";
  outfile = temp.str();
  if (fs.file_exists(outfile))
    cout << "WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if (!outf.is_open())
    spdlog::error("{}: {}", "couldn't open outfile", outfile);
  cout << "writing to file " << outfile << endl;
  outf.setf(ios_base::scientific);
  outf.precision(6);
  outf << "# Radiative Shock Test Problem outputs.  First file: " << infile
       << endl;
  outf
      << "# Columns are time and I-front position (and from nions), and should be in cgs units (s,cm).\n\n";

  // loop over all files.
  do {
    cout << "Reading from file " << infile << endl;
    // read data onto grid.
    err += dataio.ReadHeader(infile);
    err += dataio.ReadData(infile);
    if (err) spdlog::error("{}: {}", "read data went bad for file", err);

    // allocate arrays for distance and ion-fraction.
    double dist[SimPM.Ncell], ifrac[SimPM.Ncell], nions = 0.0, cvol = 0.0;
    double dx = grid->DX();
    // get first point, and for each point calculate distance to source (at
    // origin), and also the number of ionised atoms.
    cell *c = grid->FirstPt();
    int i   = 0;
    do {
      dist[i]  = GS.distance(c->x, srcpos, SimPM.ndim);
      ifrac[i] = c->P[var];
      if (SimPM.coord_sys == COORD_CYL)
        cvol = 2. * M_PI * c->x[Rcyl] * dx * dx;
      else {
        cvol = 1.0;
        for (int d = 0; d < SimPM.ndim; d++)
          cvol *= dx;
      }
      nions += (c->P[var]) * cvol;  // this is volume occupied by ions in cell.
      i++;
      if (i > SimPM.Ncell)
        spdlog::error("{}: {}", "bad ncell", SimPM.Ncell - i);
    } while ((c = grid->NextPt(c)) != 0);

    // now do a least squares fit to the data (assume uniform absolute errors)
    // pass in (npt, dist, ifrac, a_min, a_max, &a_fit, &rchisq)
    double a = 0.0, rchisq = 0.0;
    err += fit_data(
        SimPM.Ncell, dist, ifrac, grid->DX() / 2., grid->Range(XX), &a,
        &rchisq);

    // calculate a radius from total number of ions:
    if (SimPM.coord_sys == COORD_CYL)
      nions =
          pow(nions
                  / (2. * M_PI / 3.0
                     + M_PI * dx / 2. / pow(3. * nions / 2. / M_PI, 1. / 3.)),
              1. / 3.);
    else {
      if (SimPM.ndim == 2) nions = sqrt(4. * nions / M_PI);
      if (SimPM.ndim == 3) nions = pow(6. * nions / M_PI, 1. / 3.);
    }
    cout << "******* " << SimPM.simtime << "\t" << a << "\t" << rchisq << "\t"
         << nions << endl;
    outf << SimPM.simtime << "\t" << a << "\t" << rchisq << "\t" << nions
         << endl;

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


double fit_func(double r, double a)
{
  double temp = exp(15.0 * log(a / r));
  // if (a>1.e16 && r<2.e16) cout <<"r="<<r<<" and
  // f(r)="<< 1.0-exp(-temp)<<endl;
  return 1.0 - exp(-temp);
}

int fit_data(
    long int n,     // npt
    double *x,      // array of radii from centre
    double *f,      // array of ion-fractions.
    double a_min,   // a_min (min value of fit parameter)
    double a_max,   // a_max (max value of fit parameter)
    double *a_fit,  // fit value of a, returned.
    double *rcsq    // goodness of fit, returned (reduced chi-squared)
)
{
  int Na = 1000;
  double a_try, temp, chisq, min_chisq, min_a;

  min_chisq = 1.e300;
  min_a     = -1.0;
  int i     = 0;
  for (i = 0; i < Na; i++) {
    a_try = exp(log(a_min) + i * log(a_max / a_min) / (Na - 1.));
    chisq = 0.0;
    for (int j = 0; j < n; j++) {
      temp = f[j] - fit_func(x[j], a_try);
      // if (f[j]>0.0001 && f[j]<0.9999) cout
      // <<"r="<<x[j]<<"\td="<<f[j]<<"\tfit="<<f[j]-temp<<endl;
      chisq += temp * temp;
    }
    // chisq /= n-1;
    if (chisq < min_chisq) {
      min_chisq = chisq;
      min_a     = a_try;
    }
    // cout <<"a_try="<<a_try<<"\t chisq="<<chisq<<endl;
    if (chisq > 100 * min_chisq) break;
  }
  cout << "after " << i << " tries, best fit is " << min_a
       << " with chisq=" << min_chisq << endl;
  *a_fit = min_a;
  *rcsq  = min_chisq;
  return 0;
}
