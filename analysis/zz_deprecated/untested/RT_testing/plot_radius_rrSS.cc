/** \file plot_radius_rrSS.cc
 * \author Jonathan Mackey
 *
 * This file reads in a series of fits files, which are assumed to be outputs
 * of a photoevaporation problem with a single ionising source.  The data is
 * assumed to be from a simulation with radiative recombinations on, but
 * dynamics off, so that the result should follow the classical stromgen sphere
 * analysis.
 *
 * This is explicitly serial code, so if there are parallel outputs to multiple
 * files this code won't work.
 *
 * Compile with:\n
 * g++ -Wall -g -DSERIAL plot_radius_rrSS.cc ../../testing/global.cc
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

// Given a source, a geometry, and a sim-time, we can calculate the radius of
// the I-front.
int get_theoretical_radius(double *  ///< variable to put radius in.
);

// Given the radius, we replace cell values at this radius with some artificial
// value.
int overplot_radius(
    double,    ///< radius of I-front
    double *,  ///< list of radii from source.
    double *   ///< list of maxvals.
);

// Find the best-fitting radius for the I-front without any prior assumptions.
int fit_data(
    long int,  // npt
    double *,  // array of radii from centre
    double *,  // array of ion-fractions.
    double,    // a_min (min value of fit parameter)
    double,    // a_max (max value of fit parameter)
    double *,  // fit value of a, returned.
    double *,  // goodness of fit, returned (reduced chi-squared)
    double *   // optional second goodness of fit parameter
);

// output IF radius as a function of angle.
int angle_data(
    int,       ///< ncells
    double,    ///< fit radius
    double,    ///< analytic radius
    double *,  ///< array of distances from source.
    double *,  ///< array of ionisation fracitons.
    string     ///< base filename of output.
);

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc != 6) {
    cerr << "Error: must call with 5 arguments...\n";
    cerr
        << "plot_radius: <executable> <Outfile-base> <Infile-base> <First-File-Num> <File-Num-Freq> <x(HII) var>\n";
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
  cout << "Writing I-front to files " << outfilebase << "." << startct
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
    rep.error("need 1D/2D/3D sim for stromgen test", SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout << "(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid(
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (grid == 0)
    rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout << "(setup_grid) Done. g=" << grid << " and dx=" << grid->DX() << "\n";


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
  outf.precision(6);
  outf << "# Stromgen Sphere angle-averaged outputs.  First file: " << infile
       << endl;
  outf
      << "# Columns are [time, I-front position, max.dev., min.dev., posn from Nions, Theoretical Posn.], and should be in cgs units (s,cm).\n\n";

  double *dist = 0, *ifrac = 0;
  dist  = new double[SimPM.Ncell];
  ifrac = new double[SimPM.Ncell];
  if (!dist || !ifrac) rep.error("mem. alloc. dist/ifrac", dist);
  double fit_radius, analytic_radius;

  // loop over all files.
  do {
    cout << "Reading from file " << infile << endl;
    // read data onto grid.
    err += dataio.ReadHeader(infile);
    err += dataio.ReadData(infile);
    if (err) rep.error("read data went bad for file", err);

    cout << "read data, now putting distance/ifrac into arrays.\n";
    // allocate arrays for distance and ion-fraction.
    double nions = 0.0, cvol = 0.0, maxvals[SimPM.nvar];
    double dx = grid->DX();

    for (int v = 0; v < SimPM.nvar; v++)
      maxvals[v] = -1.e199;
    // get first point, and for each point calculate distance to source (at
    // origin), and also the number of ionised atoms.
    cell *c = grid->FirstPt();
    int i   = 0;
    do {
      for (int v = 0; v < SimPM.nvar; v++)
        maxvals[v] = max(maxvals[v], c->P[v]);
      // for (int v=0;v<SimPM.nvar;v++) c->Ph[v] = c->P[v];
      dist[i]  = GS.distance(c->x, srcpos, SimPM.ndim);
      ifrac[i] = c->P[var];

      if (SimPM.coord_sys == COORD_CYL)
        cvol = 2. * M_PI * c->x[Rcyl] * dx * dx;
      else {
        cvol = 1.0;
        for (int d = 0; d < SimPM.ndim; d++)
          cvol *= dx;
      }

      nions += c->P[RO] / GS.m_p() * c->P[var]
               * cvol;  // num. ions in cell (for pure Hydrogen).
      i++;
      if (i > SimPM.Ncell) rep.error("bad ncell", SimPM.Ncell - i);
    } while ((c = grid->NextPt(c)) != 0);

    // now fit a radius to the data
    // pass in (npt, dist, ifrac, a_min, a_max, &a_fit, &rchisq)
    double a = 0.0, rchisq = 0.0, rchisq2;
    err += fit_data(
        SimPM.Ncell, dist, ifrac, grid->DX() / 2., grid->Range(XX), &a, &rchisq,
        &rchisq2);
    // calculate a radius from total number of ions:
    c = grid->FirstPt();
    if (SimPM.ndim == 2) cvol = sqrt(nions * GS.m_p() / M_PI / c->P[RO]);
    if (SimPM.ndim == 3) cvol = pow(0.75 * nions * GS.m_p() / M_PI, 1. / 3.);
    cout << "******* " << SimPM.simtime << "\t" << a << "\t" << rchisq << "\t"
         << cvol << endl;
    fit_radius = a;
    outf << SimPM.simtime << "\t" << a << "\t" << rchisq << "\t" << rchisq2
         << "\t" << cvol;

    // With the theoretical radius, we can 'overplot' a radius on the data.
    err += get_theoretical_radius(&a);
    err += overplot_radius(a, dist, maxvals);
    analytic_radius = a;
    outf << "\t" << analytic_radius << endl;

    // now output the data to the output fits file.
    temp.str("");
    temp << outfilebase << "." << count << ".fits";
    outfitsfile = temp.str();
    // comment this out if we don't want to output new fits files!!!
    // err += dataio.OutputData(outfitsfile);
    if (err) rep.error("failed to overplot radius", err);

    // increment filename
    count += opfreq;
    temp.str("");
    temp << infilebase << "." << count << ".fits";
    infile = temp.str();
  } while (fs.file_exists(infile));
  // loop over all timesteps.

  // now we are at the last timestep, so write out a list of IF-radius vs. angle
  err += angle_data(
      SimPM.Ncell, fit_radius, analytic_radius, dist, ifrac, outfilebase);

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

// Given a grid and distances, calculate the I-front radius as a function of
// angle. THIS ASSUMES SOURCE IS AT THE ORIGIN!!!
int angle_data(
    int Ncell,
    double r_fit,
    double r_calc,
    double *dist,
    double *ifrac,
    string outfilebase)
{
  // open outfile
  ostringstream temp;
  temp.str("");
  temp << outfilebase << "_angle.txt";
  string outfile = temp.str();
  ofstream outf(outfile.c_str());
  if (!outf.is_open()) rep.error("couldn't open outfile", outfile);
  cout << "writing to file " << outfile << endl;
  outf.setf(ios_base::scientific);
  outf.precision(6);
  outf << "# Stromgen Sphere outputs.  radius as function of angle. r_fit="
       << r_fit << ", r_calc=" << r_calc << endl;
  outf
      << "# Columns are [cell id, ifrac, radius, phi, phi%90deg, theta, theta%90deg, r_calc [cgs units (s,cm), angles in degrees].\n\n";

  double *x;
  double phi = 0.0, theta = 0.0;

  cell *c = grid->FirstPt();
  int i   = 0;
  do {
    if (ifrac[i] >= 0.1 && ifrac[i] <= 0.9) {
      outf << c->id << "\t" << ifrac[i] << "\t" << dist[i] << "\t";
      x   = c->x;
      phi = 180. * acos(x[XX] / sqrt(x[XX] * x[XX] + x[YY] * x[YY])) / M_PI;
      if (x[YY] < RSP.position[0][YY]) phi += 180.;
      if (isnan(phi)) phi = 0.0;
      outf << phi << "\t" << (static_cast<int>(phi) + 90) % 90 << "\t";
      if (SimPM.ndim > 2) {
        theta =
            180.
            * acos(x[ZZ] / sqrt(x[XX] * x[XX] + x[YY] * x[YY] + x[ZZ] * x[ZZ]))
            / M_PI;
        if (isnan(theta)) theta = 0.0;
        outf << theta << "\t" << (static_cast<int>(theta) + 90) % 90 << "\t";
      }  // if 3d
      outf << r_calc << endl;
    }  // if at the IF
    i++;
  } while ((c = grid->NextPt(c)) != 0);
  outf.close();
  cout << "written angle data. moving on to phi,theta histograms.\n";

  double *rphi = 0, *rtheta = 0, *pmax = 0, *pmin = 0, *tmax = 0, *tmin = 0;
  int *Np = 0, *Nt = 0;
  rphi   = new double[360];
  Np     = new int[360];
  pmax   = new double[360];
  pmin   = new double[360];
  rtheta = new double[180];
  Nt     = new int[180];
  tmax   = new double[180];
  tmin   = new double[180];
  if (!rphi || !rtheta || !Np || !Nt) rep.error("mem. alloc. angle", phi);
  if (!tmax || !tmin || !pmax || !pmin) rep.error("mem. alloc. 2 angle", tmax);

  for (i = 0; i < 360; i++) {
    rphi[i] = pmax[i] = 0.0;
    Np[i]             = 0;
    pmin[i]           = 1.e200;
  }
  for (i = 0; i < 180; i++) {
    rtheta[i] = tmax[i] = 0.0;
    Nt[i]               = 0;
    tmin[i]             = 1.e200;
  }
  c = grid->FirstPt();
  int angle;
  i = 0;
  cout << "going through data.\n";
  do {
    if (ifrac[i] >= 0.1 && ifrac[i] <= 0.9) {
      x   = c->x;
      phi = 180. * acos(x[XX] / sqrt(x[XX] * x[XX] + x[YY] * x[YY])) / M_PI;
      if (x[YY] < RSP.position[0][YY]) phi += 180.;
      if (isnan(phi)) phi = 0.0;
      angle = static_cast<int>(phi);  // integer part of phi, so [0,359]
      if (angle < 0 || angle > 359) rep.error("Casting", phi);
      rphi[angle] += dist[i];
      Np[angle] += 1;
      pmax[angle] = max(pmax[angle], dist[i]);
      pmin[angle] = min(pmin[angle], dist[i]);

      if (SimPM.ndim > 2) {
        theta =
            180.
            * acos(x[ZZ] / sqrt(x[XX] * x[XX] + x[YY] * x[YY] + x[ZZ] * x[ZZ]))
            / M_PI;
        if (theta == 180.0) theta -= 1.e-10;
        if (isnan(theta)) theta = 0.0;
        angle = static_cast<int>(theta);  // integer part of phi, so [0,179]
        if (angle < 0 || angle > 179) rep.error("Casting", theta);
        rtheta[angle] += dist[i];
        Nt[angle] += 1;
        tmax[angle] = max(tmax[angle], dist[i]);
        tmin[angle] = min(tmin[angle], dist[i]);
      }  // if 3d
    }    // if at the IF
    i++;
  } while ((c = grid->NextPt(c)) != 0);

  for (i = 0; i < 360; i++) {
    rphi[i] /= max(1.0, static_cast<double>(Np[i]));
    if (Np[i] == 0) pmin[i] = 0.0;
  }
  for (i = 0; i < 180; i++) {
    rtheta[i] /= max(1.0, static_cast<double>(Nt[i]));
    if (Nt[i] == 0) tmin[i] = 0.0;
  }

  cout << "writing phi data.\n";
  temp.str("");
  temp << outfilebase << "_phi.txt";
  outfile = temp.str();
  ofstream of2(outfile.c_str());
  if (!of2.is_open()) rep.error("couldn't open outfile", outfile);
  cout << "writing to file " << outfile << endl;
  of2.setf(ios_base::scientific);
  of2.precision(6);
  of2 << "# Stromgen Sphere outputs.  radius as function of angle. r_fit="
      << r_fit << ", r_calc=" << r_calc << endl;
  of2 << "# Columns are [phi, phi%90deg, radius, N(r), r_calc [cgs units (s,cm), angles in degrees].\n\n";
  for (i = 0; i < 360; i++)
    of2 << i << "\t" << (i + 90) % 90 << "\t" << rphi[i] << "\t" << Np[i]
        << "\t" << pmax[i] << "\t" << pmin[i] << "\t" << r_calc << endl;
  of2.close();

  if (SimPM.ndim > 2) {
    cout << "writing theta data.\n";
    temp.str("");
    temp << outfilebase << "_theta.txt";
    outfile = temp.str();
    ofstream of3(outfile.c_str());
    if (!of3.is_open()) rep.error("couldn't open outfile", outfile);
    cout << "writing to file " << outfile << endl;
    of3.setf(ios_base::scientific);
    of3.precision(6);
    of3 << "# Stromgen Sphere outputs.  radius as function of angle. r_fit="
        << r_fit << ", r_calc=" << r_calc << endl;
    of3 << "# Columns are [theta, theta%90deg, radius, N(theta), r_max, r_min, r_calc [cgs units (s,cm), angles in degrees].\n\n";
    for (i = 0; i < 180; i++)
      of3 << i << "\t" << (i + 90) % 90 << "\t" << rtheta[i] << "\t" << Nt[i]
          << "\t" << tmax[i] << "\t" << tmin[i] << "\t" << r_calc << endl;
    of3.close();
  }
  cout << "all done... quitting.\n";
  delete[] rphi;
  delete[] rtheta;
  delete[] Np;
  delete[] Nt;
  delete[] tmax;
  delete[] tmin;
  delete[] pmax;
  delete[] pmin;
  return 0;
}



// Given a source, a geometry, and a sim-time, we can calculate the radius of
// the I-front.
int get_theoretical_radius(double *a  ///< variable to put radius in.
)
{
  cell *c = grid->FirstPt();
  for (int i = 1; i < SimPM.NG[XX]; i += 2)
    c = grid->NextPt(c, XP);
  for (int i = 1; i < SimPM.NG[YY]; i += 2)
    c = grid->NextPt(c, YP);
  if (SimPM.ndim > 2)
    for (int i = 1; i < SimPM.NG[ZZ]; i += 2)
      c = grid->NextPt(c, ZP);
  // rep.printVec("x",c->x,SimPM.ndim);
  // should be at centre of grid now.

  // assume constant density, with recombinations
  double nh    = c->P[RO] / GS.m_p();  // number density of hydrogen.
  double alpha = 2.59e-13;             // recombination coeff.
  double rad   = RSP.strength[0] / nh / nh / alpha / M_PI;
  rad *= (1.0 - exp(-nh * alpha * SimPM.simtime));
  if (SimPM.coord_sys == COORD_CRT && SimPM.ndim == 2) {
    rad = sqrt(rad);
  }
  else if (
      (SimPM.coord_sys == COORD_CYL && SimPM.ndim == 2)
      || (SimPM.coord_sys == COORD_CRT && SimPM.ndim == 3)) {
    rad = pow(0.75 * rad, 1. / 3.);
  }
  else
    rep.error("Bad ndim/coords combo in get_theoretical_radius()", SimPM.ndim);
  /*
  // assume constant density!
  double rad = RSP.strength[0]*SimPM.simtime*GS.m_p()/M_PI/c->P[RO];
  if       (SimPM.coord_sys==COORD_CRT && SimPM.ndim==2) {
    rad = sqrt(rad);
  }
  else if ((SimPM.coord_sys==COORD_CYL && SimPM.ndim==2) ||
     (SimPM.coord_sys==COORD_CRT && SimPM.ndim==3)) {
    rad = pow(0.75*rad, 1./3.);
  }
  else rep.error("Bad ndim/coords combo in
  get_theoretical_radius()",SimPM.ndim);
  */
  /*
  // 1/r^2 density.
  double n0 = c->P[RO]/GS.m_p();
  double r0 = (SimPM.Xmax[XX]-SimPM.Xmin[XX])/10.0;
  if (SimPM.coord_sys==COORD_CRT && SimPM.ndim==2) {
    rad = r0*sqrt(exp(RSP.strength[0]*SimPM.simtime/M_PI/n0/r0/r0)-1.0);
  }
  else rep.error("only 2d radial profiles calculated!",SimPM.ndim);
  */

  *a = rad;
  return 0;
}


// Given the radius, we replace cell values at this radius with some artificial
// value.
int overplot_radius(
    double a,        ///< radius of I-front
    double *dist,    ///< list of radii of cells from source.
    double *maxvals  ///< list of max vals on grid.
)
{
  double del = 0.0;
  if (SimPM.ndim == 2)
    del = 0.7 * sqrt(2.0) * grid->DX() / 2.0;
  else if (SimPM.ndim == 3)
    del = 0.7 * sqrt(3.0) * grid->DX() / 2.0;
  for (int v = 0; v < SimPM.nvar; v++) {
    maxvals[v] *= 10;
    if (maxvals[v] < 0.0) maxvals[v] = 1.0;
  }

  cell *c     = grid->FirstPt();
  long int ct = 0, num = 0;
  do {
    if (fabs(a - dist[ct]) <= del) {
      for (int v = 0; v < SimPM.nvar; v++)
        c->P[v] = maxvals[v];
      num++;
    }
    ct++;
  } while ((c = grid->NextPt(c)) != 0);
  cout << "got " << num << " cells in radius.  surface is roughly ";
  if (SimPM.ndim == 2)
    del = 2.0 * M_PI * a / grid->DX();
  else if (SimPM.ndim == 3)
    del = 4.0 * M_PI * a * a / grid->DX() / grid->DX();
  cout << del << " cells in extent.\n";
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
    double *rcsq,   // goodness of fit, returned
    double *rcsq2   // optional second goodness of fit parameter
)
{
  // take average radius of cells with 0.1<x<0.9
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
    *a_fit = 0.0;
    *rcsq  = 0.0;
    *rcsq2 = 0.0;
  }
  cout << "\tFound " << ct << " cells near I-front, with mean radius "
       << *a_fit;
  cout << " and width " << *rcsq - *rcsq2 << "\n";
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
  cout << "\tFound " << ct << " cells near I-front, with mean radius "
       << *a_fit;
  cout << " and width " << *rcsq - *rcsq2 << "\n";


  /*
  // chi-squared minimization of what the radius is.
  int Na=1000;
  double a_try, temp, chisq, min_chisq, min_a;
  min_chisq=1.e300; min_a=-1.0;
  int i=0;
  for (i=0;i<Na;i++) {
    a_try=exp( log(a_min) + i*log(a_max/a_min)/(Na-1.) );
    chisq=0.0;
    for (int j=0;j<n;j++) {
      temp = f[j]-fit_func(x[j],a_try);
      //if (f[j]>0.0001 && f[j]<0.9999) cout
  <<"r="<<x[j]<<"\td="<<f[j]<<"\tfit="<<f[j]-temp<<endl; chisq += temp*temp;
    }
    //chisq /= n-1;
    if (chisq < min_chisq) {
      min_chisq=chisq; min_a=a_try;
    }
    //cout <<"a_try="<<a_try<<"\t chisq="<<chisq<<endl;
    if (chisq>100*min_chisq) break;
  }
  cout <<"after "<<i<<" tries, best fit is "<<min_a<<" with
  chisq="<<min_chisq<<endl; *a_fit = min_a; *rcsq = min_chisq; return 0;
  */
}
