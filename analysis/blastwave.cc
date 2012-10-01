/** \file blastwave.cc
 * \author Jonathan Mackey
 * 
 * This file reads in a fits file, which is assumed to be output
 * of a blastwave problem with a source at the origin.
 * 
 * This is explicitly serial code, so if there are parallel outputs to multiple
 * files this code won't work.
 * 
 * Compile with:\n
 * g++ -Wall -DSERIAL blastwave.cc ../testing/global.cc ../testing/uniformGrid.cc ../testing/dataio.cc -lreadline -lcfitsio
 * 
 * Run with (example):\n
 * ./a.out profile.txt ../results/BWfilename
 * blastwave: <executable> <Outfile> <infile>
 * */

#include "fitsio.h"
using namespace std;
#include "../testing/global.h"
#include "../testing/uniformGrid.h"
#include "../testing/dataio.h"


int fit_data(long int, // npt
	     double *, // array of radii from centre
	     double *, // array of ion-fractions.
	     double,   // a_min (min value of fit parameter)
	     double,   // a_max (max value of fit parameter)
	     double *, // fit value of a, returned.
	     double *  // goodness of fit, returned (reduced chi-squared)
	     );

int main(int argc, char **argv)
{
  // Get two input files and one output file from cmd-line args.
  if (argc!=3) {
    cerr << "Error: must call with 3 arguments...\n";
    cerr << "blastwave: <executable> <OutfileName> <infile> \n";
    rep.error("Bad number of Args",argc);
  }
  string outfile = argv[1];
  string infile  = argv[2];
  cout <<"reading from file "<<infile<<"\n";
  cout <<"Writing cell data to file "<<outfile<<"\n";
  cout <<"**********************************************\n";

  class DataIOFits dataio;
  class file_status fs;
  
  // First we need to open the first file in the list, and get the grid dimensions,
  // so we can set it up once and use it for all the infiles.
  int err=0;
  if (!fs.file_exists(infile)) rep.error("First file not found!",infile);
  err = dataio.ReadHeader(infile);
  if (err) rep.error("read header went bad",err);
  
  // check dimensionality is ok.
  if (SimPM.ndim!=1 && SimPM.ndim!=2 &SimPM.ndim!=3)
    rep.error("need 1D/2D/3D sim for BW test",SimPM.ndim);
  // Now the header should contain the sim dimensionality, number of vars,
  // size of box, so we can use these to set up the grid.
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
  grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  if (grid==0) rep.error("(IntUniformFV::setup_grid) Couldn't assign data!", grid);
  cout <<"(setup_grid) Done. g="<<grid<<"\n";

  
  double srcpos [SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++) srcpos[i]=0.0;
  // Set up and open outfile
  if (fs.file_exists(outfile)) cout <<"WARNING:: file exists, I am overwriting a text file.\n";
  ofstream outf(outfile.c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",outfile);
  cout <<"writing to file "<<outfile<<endl;
  outf.setf( ios_base::scientific );
  outf.precision(6);
  outf << "# BW output.  file: "<<infile<<endl;
  outf << "# simtime="<<SimPM.simtime<<endl;
  outf << "# Columns are distance density pressure velocity \n\n";

  cout <<"Reading from file "<<infile<<endl;
  // read data onto grid.
  err += dataio.ReadHeader(infile);
  err += dataio.ReadData(infile);
  if (err) rep.error("read data went bad for file",err);

  // allocate arrays for distance and ion-fraction.
  double dx=grid->DX(), d, vel;
  vector<double> dist, pg;
  // get first point, and for each point calculate distance to source (at origin).
  cell *c = grid->FirstPt();
  double p_in = c->P[PG];
  c= grid->LastPt();
  double p_out = c->P[PG];
  double pratio = p_in/p_out;
  if (pratio<50) rep.error("low pressure ratio",pratio);

  long int ct=0;
  c = grid->FirstPt();
  do {
    //    d = GS.distance(c->x,srcpos,SimPM.ndim);
    d=0.0;
    for (int i=0;i<SimPM.ndim;i++) d += c->x[i]*c->x[i];
    d = sqrt(d);
    if (d<SimPM.Range[XX]) {
      if (ct%11==0) {
	outf << d/3.086e18 <<"\t"<< c->P[RO]/2.0e-24 <<"\t"<< c->P[PG];
	vel = c->P[VX]*c->P[VX];
	if (SimPM.ndim>1) vel+= c->P[VY]*c->P[VY];
	if (SimPM.ndim>2) vel+= c->P[VZ]*c->P[VZ];
	vel=sqrt(vel);
	outf <<"\t"<< vel<<endl;
      }
    }
    
    if (c->P[PG]>3.0*p_out && c->P[PG]<0.5*p_in) {
      dist.push_back(d);
      pg.push_back(c->P[PG]);
      //      cout <<"d="<<d<<" and p="<<c->P[PG]<<endl;
    }
    ct++;
  } while ((c=grid->NextPt(c)) !=0);
  
  double sum=0.0, lv=1.e99, uv=0.0;
  for (unsigned int i=0; i<dist.size(); i++) {
    sum += dist[i];
    lv = min(lv,dist[i]);
    uv = max(uv,dist[i]);
  }
  sum /= dist.size();
  cout <<"****** "<<SimPM.simtime<<"\t"<<sum<<"\t"<<lv<<"\t"<<uv<<"\t"<<dx<<endl;
  cout <<"\n***************************************************\n";
  outf.close();
  delete grid; grid=0;
  return 0;
}


double fit_func(double r, double a)
{
  double temp=exp(15.0*log(a/r));
  //if (a>1.e16 && r<2.e16) cout <<"r="<<r<<" and f(r)="<< 1.0-exp(-temp)<<endl;
  return 1.0-exp(-temp);
}
 
int fit_data(long int n, // npt
	     double *x, // array of radii from centre
	     double *f, // array of ion-fractions.
	     double a_min,   // a_min (min value of fit parameter)
	     double a_max,   // a_max (max value of fit parameter)
	     double *a_fit, // fit value of a, returned.
	     double *rcsq  // goodness of fit, returned (reduced chi-squared)
	     )
{
  int Na=1000;
  double a_try, temp, chisq, min_chisq, min_a;

  min_chisq=1.e300; min_a=-1.0;
  int i=0;
  for (i=0;i<Na;i++) {
    a_try=exp( log(a_min) + i*log(a_max/a_min)/(Na-1.) );  
    chisq=0.0;
    for (int j=0;j<n;j++) {
      temp = f[j]-fit_func(x[j],a_try);
      //if (f[j]>0.0001 && f[j]<0.9999) cout <<"r="<<x[j]<<"\td="<<f[j]<<"\tfit="<<f[j]-temp<<endl;
      chisq += temp*temp;
    }
    //chisq /= n-1;
    if (chisq < min_chisq) {
      min_chisq=chisq; min_a=a_try;
    }
    //cout <<"a_try="<<a_try<<"\t chisq="<<chisq<<endl;
    if (chisq>100*min_chisq) break;
  }
  cout <<"after "<<i<<" tries, best fit is "<<min_a<<" with chisq="<<min_chisq<<endl;
  *a_fit = min_a;
  *rcsq = min_chisq;
  return 0;
}
