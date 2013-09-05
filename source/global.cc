///
/// \file global.cc
///  
/// \brief Class definitions for global classes and structures.
/// \author Jonathan Mackey
/// 
/// This file contains the class definitions for all global classes.
/// These include SimParams, ParallelParams, JetParams, reporting,
/// GeneralStuff, etc.
///  
/// Modifications:
/// - 2007-06-27 Worked on direction-dependent functions in the equations classes.
///  Think I have got it much better -- no function pointers, but use local vector indices instead.
/// - 2007-07-05 Added in glmMHDEqn class, for the divergence cleaning method.
/// - 2007-07-13 fixed up things in glmclass.
/// - 2007-07-16 simplified inheritance structure so noly baseeqn is virtual.
/// - 2007-10-11 Added interactive command line interface class.
/// - 2007-10-16 Moved equations classes into equations.cc
/// - 2008-08-27 moved domain decomposition to ParallelParams.
/// - 2010-02-05 JM: Added offsets to parallelparams.
/// - 2010-02-06 JM: Added function returning all abutting domains.
/// - 2010-04-25 JM: Added min_timestep parameter to SimPM to bug out
///   if dt gets too small.
/// - 2010-07-21 JM: EP_update_erg is now an int rather than a bool (KISS).
/// - 2010-07-24 JM: Added stellar wind class.
/// - 2010.07.26 JM: Fixed stellar wind divergence for 2D slab-symmetry
/// - 2010.10.01 JM: Spherical coordinates added (distance formulae).
///    Cut out testing myalloc/myfree
/// - 2010.10.05 JM: Moved stellar winds to their own file.
/// - 2010.10.13 JM: Added a function to display command-line options.
/// - 2010.11.03 JM: switched endl for c-style end of line.
/// - 2010.11.12 JM: Moved cell interface away from global.cc to
///   cell_interface.cc
/// - 2010.11.19 JM: Got rid of testing myfree(,str) in parallel params.   
/// - 2010.12.04 JM: Added text warning against using the distance
///   functions in GeneralStuff!
/// - 2011.01.12 JM: Commented out redirection of stderr -- too many
///   files for parallel simulations (undone 14/1/11).
/// - 2011.01.17 JM: Added checkpt_freq=N steps to commandline options.
/// - 2011.02.24 JM: Added Struct SimPM.RS to handle radiation sources more simply.
///     Also added a new function to RSP so that I can query it more easily.
/// - 2011.02.25 JM: Updated decomposedomain() so that it only inhibits decomposition
///     along a given axis if we only have one ionising source and no diffuse sources.
/// - 2011.02.28 JM: Got rid of RSP class.  It was too complicated.
/// - 2011.06.02 JM: RefVec is now a stack array rather than dynamically allocated.
/// - 2011.10.14 JM: Commented out RT_DIFF class
/// - 2011.12.01 JM: Added GeneralStuff::root_find_linear()
/// - 2012.05.14 JM: Updated decomposedomain() so it now also checks if there
///    are multiple sources at infinity in different directions, in which case
///    we can't make all of them fully parallel, so at_infinity --> false.
/// - 2012.05.15 JM: fixed the logic of decomposedomain()
/// - 2013.01.17 JM: Changed redirect so that there are many files
///    only when TESTING is defined; otherwise just rank 0 writes to
///    file, and all others have stdout supressed.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.02.14 JM: Added He/Metal mass fractions as EP parameters,
///    to make metallicity and mu into parameterfile settings.
/// - 2013.04.18 JM: Removed NEW_METALLICITY flag.
/// - 2013.08.19 JM: Changed constants treatment in GeneralStuff.
///    Will remove them from here eventually.
/// - 2013.09.05 JM: changed RS position[] to pos[].

#include "global.h"
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <stdexcept>
using namespace std;

class SimParams SimPM;
class JetParams JP;
class reporting rep;
class Units uc;
class GeneralStuff GS;
//class stellar_wind SW;
struct stellarwind_list SWP;

//---------------------------------------------------------
//---------------- GRID BASE CLASS POINTER ----------------
//---------------------------------------------------------
class GridBaseClass *grid=0;
class IntegratorBaseFV *integrator=0;
//---------------------------------------------------------

//---------------------------------------------------------
//
// Function to output commandline options for code.
//
void print_command_line_options(int argc, char **argv)
{
  cout <<"You ran:\n";
  for (int v=0;v<argc;v++)
    cout <<"  "<<argv[v];
  cout <<"\n      ************************         \n";
  cout << argv[0] <<": must call with at least 3 arguments...\n";
  cout <<" <main> <icfile> <typeoffile> <solver-type> [optional args]\n";
  cout <<"Parameters:\n";
  cout <<"<icfile> \n";
  cout <<"\tCan be an ASCII parameter-file for 1D and 2D shocktubes.\n";
  cout <<"\tOtherwise should be an initial-condition file or restart-file\n";
  cout <<"\tin FITS or Silo format.\n";
  cout <<"<typeoffile> \n";
  cout <<"\tInteger flag to tell me what type of file I am starting from.\n";
  cout <<"\tCan be one of [1=text paramfile, 2=FITS, 5=Silo file].\n";
  cout <<"<solvetype> \n";
  cout <<"\tInteger =1 for uniform finite-volume, no other options.\n";

  cout <<"\n";
  cout <<"[optional args] are in the format <name>=<value> with no spaces.\n\n";

  cout <<"\n*********** DATA I/O OPTIONS ************\n";
  cout <<"\t redirect=string : filename with path to redirect stdout/stderr to\n";

  cout <<"\t checkpt_freq=N  : checkpoint every N timesteps (default is 250).\n";
  cout <<"\t op_criterion=N  : 0=output every I steps, 1=output every D time units.\n";
  cout <<"\t opfreq=N        : Output data every Nth timestep  (if op_criterion=0).\n";
  cout <<"\t opfreq_time=D   : Output data every Dth time unit (if op_criterion=1).\n";
  cout <<"\t finishtime=D    : set time to finish simulation, in code time units.\n";
  cout <<"\t optype=S        : Specify type of output file,";
  cout <<             " [1,text]=TEXT,[2,fits]=FITS,[4,both]=FITS+TEXT,[5,silo]=SILO,[6]=SILO+TEXT.\n";
  cout <<"\t outfile=NAME    : Replacement output filename, with path.\n";

  cout <<"\n*********** PHYSICS OPTIONS *************\n";
  cout <<"\t ooa=N         : modify order of accuracy (either 1 or 2).\n";

  cout <<"\t eqntype=N     : modify type of equations,";
  cout <<" Euler=1, Ideal-MHD=2, GLM-MHD=3, FCD-MHD=4, IsoHydro=5.\n";
  cout <<"\t\t\t WARNING -- IT IS DANGEROUS TO CHANGE EQUATIONS!\n";

  cout <<"\t AVtype=N      : modify type of artificial viscosity:";
  cout <<" 0=none, 1=Falle,Komissarov,Joarder(1998), 2=Colella+Woodward(1984), 3=Sanders et al.(1998)[H-correction].\n";
  cout <<"\t\t\t WARNING -- AVtype=2 IS NOT WORKING WELL.  ONLY USE FKJ98.";
  cout <<" The H-correction works well in serial, if needed.\n";
  cout <<"\t EtaVisc=D     : modify viscosity parameter to the given double precision value.\n";

  cout <<"\t coordsys=NAME : override coordinate system to [cartesian,cylindrical]. !DANGEROUS!\n";
  cout <<"\t cfl=D         : change the CFL no. for the simulation, in range (0,1).\n";
  cout <<"\t cooling=N     : cooling=0 for no cooling, 1 for Sutherland&Dopita1993.\n";
  cout <<"\t\t\t For other cooling functions see cooling.cc/cooling.h.\n";

  cout <<"\t solver=N      :\n";
  cout <<"\t\t 0 = Lax-Friedrichs Flux, TESTING ONLY!\n";
  cout <<"\t\t 1 = Linear Riemann Solver : EULER/MHD (mean value average)";
  cout  <<" (Falle, Komissarov, Joarder, 1998),\n";
  cout <<"\t\t 2 = Exact Riemann Solver  : EULER ONLY (Hirsch (199X), Toro, 1999)\n";
  cout <<"\t\t 3 = Hybrid Riemann Solver (1+2)         : EULER ONLY \n";
  cout <<"\t\t 4 = Roe Conserved Variables flux solver : EULER/MHD";
  cout  <<" (e.g. Toro, 1999, Stone, Gardiner et al. 2008)\n";
  cout <<"\t\t 5 = Roe Primitive Variables flux solver : EULER ONLY";
  cout  <<" (e.g. Stone, Gardiner et al. 2008)\n";
  cout <<"\t\t 6 = Flux vector splitting : EULER ONLY (van Leer, 1982) \n";

  cout <<"\t timestep_limit=N:\n";
  cout <<"\t\t 0 = only dynamical Courant condition.\n";
  cout <<"\t\t 1 = dynamical + cooling time limits\n";
  cout <<"\t\t 2 = dyn +cool +recombination time limits\n";
  cout <<"\t\t 3 = dyn +cool +recomb +ionisation time limits\n";
  cout <<"\t\t 4 = dyn +recomb (NO cool, NO ion)\n";

  cout <<"\n*********** PARALLEL CODE ONLY *************\n";
  cout <<"\t maxwalltime=D : change the max. runtime to D in seconds.\n";
  

  cout <<"\n";
  cout <<"********* DEPRECATED -- STILL HERE FOR LEGACY SCRIPTS... ******\n";
  cout <<"\t artvisc=D : modify artificial viscosity, 0=none, Otherwise FalleAV with eta=D,\n";
  cout <<"\t noise=D   : add noise to initial conditions if desired, at fractional level of D.\n";
  cout <<"     *********************************************\n\n";
  return;
}
//---------------------------------------------------------


/************************ MULTI-PROCESS COMMS *************************/
#ifdef PARALLEL
#if   defined USE_MPI
class comms_base *COMM = new comm_mpi ();
#elif defined USE_FILE_COMMS
class comms_base *COMM = new comm_files ();
#else
#error "MUST DEFINE EITHER USE_MPI or USE_FILE_COMMS"
#endif

#endif // PARALLEL

/************************ MULTI-PROCESS COMMS *************************/

/************************ MICROPHYSICS ***************************/
class MicroPhysicsBase *MP=0;
/************************ MICROPHYSICS ***************************/

/************************  RAYTRACING  ***************************/
class RayTracingBase *RT=0;
/************************  RAYTRACING  ***************************/

/************************* MEMORY MANAGEMENT ***********************/
class memory_management mem;
/************************* MEMORY MANAGEMENT ***********************/


/************************* CELL INTERFACE ***********************/
#ifdef COUNT_ENERGETICS
struct energetics *GLOBAL_CE=0;
///< for tracking rates in microphysics/raytracing.
#endif
class cell_interface CI;
///< global class for accessing cell data, positions, neigbours.
/************************* CELL INTERFACE ***********************/




//------------------------------------------------
//--- General Stuff Class ------------------------
GeneralStuff::GeneralStuff()
{
  timers.clear();
}

GeneralStuff::~GeneralStuff() 
{
  //cout <<"timers.size: "<<timers.size()<<"\n";
  timers.clear();
}

bool GeneralStuff::equalD(const double a, const double b)
{
  if(a==b) return(true);
  if(fabs(a)+fabs(b)<TINYVALUE) {
    cout <<"tiny numbers in equalD(a,b); a,b <1.e-100... a="<<a<<", b="<<b<<"; returning true.\n";
    return(true);
  }
  if( (fabs(a-b)/(fabs(a)+fabs(b)+TINYVALUE)) <SMALLVALUE) return(true); // true is 1
  else return(false); // false is zero.
}

void GeneralStuff::spline(const double *x,
			  const double *y,
			  const int n,
			  double yp1,
			  double ypn,
			  double *y2
			  )
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double [n];
  
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  //cout <<"x[0]="<<x[0]<<" y[0]="<<y[0]<<" y2[0]="<<y2[0]<<"\n";
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    //cout <<"x[i]="<<x[i]<<" y[i]="<<y[i]<<" y2[i]="<<y2[i]<<"\n";
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  //rep.printVec("y2",y2,50);
  delete [] u;
  return;
}

void GeneralStuff::splint(const double xa[],
			  const double ya[],
			  const double y2a[],
			  const int n,
			  const double x,
			  double *y
			  )
{
  int klo,khi,k;
  double h,b,a;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
  }
  //cout <<"khi="<<khi<<" klo="<<klo;
  //cout <<"\t\txhi="<<xa[khi]<<" xlo="<<xa[klo];
  //cout <<"\t\tyhi="<<ya[khi]<<" ylo="<<ya[klo];
  //cout <<"\t\ty2hi="<<y2a[khi]<<" y2lo="<<y2a[klo]<<"\n";
  h=xa[khi]-xa[klo];
  //if (h == 0.0) { rep.error("Bad xa input to routine splint",h); }
  if (h < VERY_TINY_VALUE) { rep.error("Bad xa input to routine splint",h); }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}


// ##################################################################
// ##################################################################



void GeneralStuff::root_find_linear(
        const double *xarr, ///< Array of x values.
        const double *yarr, ///< Array of y values.
        const size_t len,     ///< Array sizes
        const double xreq,  ///< x we are searching for.
        double *yreq  ///< pointer to result.
        )
{
  //
  // Given a vector of x-values, and corresponding y-values, and an input
  // x value, find the corresponding y-value by bisection and then linear
  // interopolation.
  //
  // First we find the two x-points in the array which bracket the requested
  // x-value, with bisection.
  //
  size_t
    ihi = len,  // upper bracketing value
    ilo = 0,    // lower bracketing value
    imid= 0;    // midpoint
  //int count=0;
  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (xarr[imid] < xreq) ilo = imid;
    else                   ihi = imid;
    //count ++;
  } while (ihi-ilo >1);
  //cout.precision(12);
  //cout <<"count="<<count<<", ihi="<<ihi<<" and ilo="<<ilo<<", x[ihi]=";
  //cout <<xarr[ihi]<<" and x[ilo]="<<xarr[ilo]<<" and xreq="<<xreq<<"\n";
  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  *yreq = yarr[ilo] + (yarr[ihi]-yarr[ilo])*(xreq-xarr[ilo])/(xarr[ihi]-xarr[ilo]);
  return;
}


// ##################################################################
// ##################################################################



double GeneralStuff::idistance(const int *p1, ///< position 1
			       const int *p2, ///< position 2
			       const int nd     ///< number of spatial dimensions.
			       )
{
  //
  // THIS SHOULD BE A GRID MEMBER FUNCTION!  AND I SHOULD HAVE
  // SEPARATE GRIDS FOR DIFFERENT COORDINATE SYSTEMS!!!
  // (2011.11.24 JM: Grid-classes now have functions to call instead of this one).
  //
  double temp=0.0;
  if      (SimPM.coord_sys==COORD_CRT) {
    for (int i=0;i<nd;i++)
      temp += static_cast<double>((p1[i]-p2[i])*(p1[i]-p2[i]));
  }
  else if (SimPM.coord_sys==COORD_CYL) {
    if (nd<2) rep.error("cylindrical coords but ndim<2",nd);
    temp += static_cast<double>((p1[Zcyl]-p2[Zcyl])*(p1[Zcyl]-p2[Zcyl]));
    //
    // THIS SHOULD BE TO THE CENTRE-OF-VOLUME OF THE CELL, SINCE THAT
    // IS WHERE EVERYTHING IS CALCULATED!! THIS MAY BE WHY THE STELLAR
    // WINDS ARE NOT WORKING IN CYLINDRICAL COORDINATES!! I NEED TO
    // DISTINGUISH BETWEEN GETTING DISTANCES BETWEEN CELL--CENTRES AND
    // DISTANCES BETWEEN ARBITRARY POINTS, AND REALLY THE INTEGER
    // DISTANCE FUNCTION IS ONLY FOR EITHER CELL VERTICES OR CENTRES
    // WHICH ARE TWO SPECIAL CASES AND CAN BE TREATED SPECIALLY
    // (2011.11.24 JM: Grid-classes now have functions to call instead of this one).
    //
    temp += static_cast<double>((p1[Rcyl]-p2[Rcyl])*(p1[Rcyl]-p2[Rcyl]));
    if (nd>2) temp += 2.0*static_cast<double>(p1[Rcyl]*p2[Rcyl])*(1.0-cos(static_cast<double>(p1[Tcyl]-p2[Tcyl])));
  }
  else if (SimPM.coord_sys==COORD_SPH) {
    temp += static_cast<double>((p1[Rsph]-p2[Rsph])*(p1[Rsph]-p2[Rsph]));
    if (nd>1) rep.error("Define integer system for spherical coords!",nd);
  }
  else rep.error("Bad coord sys in SimPM",SimPM.coord_sys);
  return sqrt(temp);
}
  

// ##################################################################
// ##################################################################



double GeneralStuff::distance(const double *p1, ///< position 1
			      const double *p2, ///< position 2
			      const int nd     ///< number of spatial dimensions.
			      )
{
  //
  // THIS SHOULD BE A GRID MEMBER FUNCTION!  AND I SHOULD HAVE
  // SEPARATE GRIDS FOR DIFFERENT COORDINATE SYSTEMS!!! AND I NEED TO
  // DISTINGUISH BETWEEN CALLS TO CELL--CENTRES, WHICH SHOULD PASS
  // CELL--POINTERS, AND CALLS FOR DISTANCES BETWEEN RANDOM POINTS.
  // (2011.11.24 JM: Grid-classes now have functions to call instead of this one).
  //
  double temp=0.0;
  if      (SimPM.coord_sys==COORD_CRT) {
    for (int i=0;i<nd;i++)
      temp += (p1[i]-p2[i])*(p1[i]-p2[i]);
  }
  else if (SimPM.coord_sys==COORD_CYL) {
    if (nd<2) rep.error("cylindrical coords but ndim<2",nd);
    temp += (p1[Zcyl]-p2[Zcyl])*(p1[Zcyl]-p2[Zcyl]);
    temp += (p1[Rcyl]-p2[Rcyl])*(p1[Rcyl]-p2[Rcyl]);
    if (nd>2) temp += 2.0*p1[Rcyl]*p2[Rcyl]*(1.0-cos(p1[Tcyl]-p2[Tcyl]));
  }
  else if (SimPM.coord_sys==COORD_SPH) {
    if (nd==1)
      temp += (p1[Rsph]-p2[Rsph])*(p1[Rsph]-p2[Rsph]);
    if (nd>1) rep.error("Code distance formula for spherical coords!",nd);
    // 2D: sqrt( (R1*cosT1-R2*cosT2)^2 +(R1*sinT1-R2*sinT2)^2 )
    // 3D: Same, but with [RsinTcosP, RsinTsinP, RcosT]
  }
  else rep.error("Bad coord sys in SimPM",SimPM.coord_sys);
  return sqrt(temp);
}



// ##################################################################
// ##################################################################



void GeneralStuff::start_timer(string id)
{
  struct timeval s;
  gettimeofday(&s,0);
  double t = s.tv_sec +1.e-6*s.tv_usec;

  //
  // If timer exists, find it, and restart the timer.  It should have been
  // paused previously, in which case it's current value is the time it has
  // been running so far, so we set the value to the current time minus its
  // value, which is an 'effective' start time.
  // If the timer doesn't exist, initialise it to the current time in seconds.
  //
  if (timers.find(id)!=timers.end()) {
    timers[id] = t-timers[id];
  }
  else {
    timers[id] = t;
  }
  //cout << "id="<<id<<" start="<<timers[id]<<"\n";
  return;
}

// ##################################################################
// ##################################################################

double GeneralStuff::pause_timer(string id)
{
  struct timeval s;
  gettimeofday(&s,0);
  double t = s.tv_sec +1.e-6*s.tv_usec;
  //
  // Set timer to be the number of seconds it has been running.
  //
  timers[id] = t-timers[id];
  return timers[id];
}

// ##################################################################
// ##################################################################

double GeneralStuff::stop_timer(string id)
{
  struct timeval s;
  gettimeofday(&s,0);
  double t = s.tv_sec +1.e-6*s.tv_usec;
  //cout <<" start="<<timers[id];
  //cout <<" end  ="<<t;
  t -= timers[id];

  //
  // Delete the timer.
  //
  timers.erase(id);
  //cout <<" time="<<t<<"\n";
  return t;
}

// ##################################################################
// ##################################################################

double GeneralStuff::time_so_far(string id)
{
  struct timeval s;
  gettimeofday(&s,0);
  double t = s.tv_sec +1.e-6*s.tv_usec;
  //cout <<" start="<<timers[id];
  //cout <<" now ="<<t;
  //cout <<" timesofar="<<t-timers[id];
  return t-timers[id];
}

//------------------------------------------------


//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------
#ifdef PARALLEL
class ParallelParams mpiPM;

ParallelParams::ParallelParams()
{
   nproc = -1;
   myrank = -1;
   LocalNcell = -1;
   for (int i=0; i<MAX_DIM; i++) {
     LocalNG[i] = offsets[i] = ix[i] = nx[i] = -1;
     LocalXmin[i] = LocalXmax[i] = LocalRange[i] = -1.e99;
   }
   ngbprocs=0;
   ReadSingleFile  =true; ///< If the ICs are in a single file, set this to true.
   WriteSingleFile =false; ///< If you want all the processors to write to one file, set this (BUGGY!)
   WriteFullImage  =false; ///< If you want multiple fits files, but each one is the full domain size (bad!), set this.
   max_walltime = 8.64e6; ///< default runtime is 10 days. Can be reset at commandline.
}

ParallelParams::~ParallelParams() {
  if (ngbprocs!=0) ngbprocs = mem.myfree(ngbprocs);
  return;
}

double ParallelParams::get_max_walltime()
{return max_walltime;}

void ParallelParams::set_max_walltime(double t ///< New Max. runtime in seconds.
				      )
{
  cout <<"\tResetting max. walltime from "<<max_walltime;
  max_walltime = t;
  cout <<" to new value: "<<max_walltime<<"\n";
}

int ParallelParams::decomposeDomain()
{
  //  cout << "---ParallelParams::decomposeDomain() decomposing domain.\n";

  //  
  // First check if the Local domain has been read in from a restart file.
  // If it has, we don't want to override it, so just check it is a valid
  // subdomain, set up pointers to neighbouring processors, and return.
  //
  // This is defunct with Silo I/O, and there is no harm in doing the decomposition
  // again even if I've read in everything from the fits header, so I'm commenting
  // this out (and I sometimes want to override the decomposition in analysis code).
  //
//   if (mpiPM.LocalXmax[XX] > -1.e90) {
//     cout <<"\tdecomposeDomain() Domain already set. Must be read from paramfile.\n";
//     cout <<"\t Will check it is a sensible domain and continue.\n";
//     for (int i=0;i<SimPM.ndim; i++) {
//       if (mpiPM.LocalXmax[i] <= mpiPM.LocalXmin[i]) rep.error("Bad domain",mpiPM.LocalRange[i]);
//       if (!GS.equalD(mpiPM.LocalXmax[i]-mpiPM.LocalXmin[i],mpiPM.LocalRange[i]))
// 	rep.error("Bad domain", (mpiPM.LocalXmax[i]-mpiPM.LocalXmin[i])/mpiPM.LocalRange[i]);
//       if (mpiPM.LocalRange[i] > SimPM.Range[i]) rep.error("range",mpiPM.LocalRange[i]);
//       if (mpiPM.LocalXmin[i]  < SimPM.Xmin[i] ) rep.error("range",mpiPM.LocalXmin[i]);
//       if (mpiPM.LocalXmax[i]  > SimPM.Xmax[i] ) rep.error("range",mpiPM.LocalXmax[i]);
//       if (mpiPM.LocalNG[i]    > SimPM.NG[i]   ) rep.error("range",mpiPM.LocalNG[i]);
//       if (mpiPM.LocalNcell    > SimPM.Ncell   ) rep.error("range",mpiPM.LocalNcell);
//     }
//     pointToNeighbours();
//     return 0;
//   }
  
  //
  // If the domain isn't already determined, we set up the simplest one 
  // possible, which is to subdivide the domain in half recursively, where
  // the axis we cut along is always the one in which the domain is longest.
  // In this way we minimize domain interfaces, which should minimize 
  // communication.
  //
  if (SimPM.ndim==1) {
    // 1D is so simple it's worth putting in it's own section.
    LocalRange[XX] = SimPM.Range[XX]/nproc;
    LocalXmin[XX] = SimPM.Xmin[XX] + myrank*LocalRange[XX];
    LocalXmax[XX] = LocalXmin[XX] + LocalRange[XX];
    LocalNG[XX] = SimPM.NG[XX]/nproc;
    offsets[XX] = 0;
    LocalNcell = LocalNG[XX];
    // Set up ngbprocs array to point to neighbouring processors
    nx[XX] = nproc;
    ix[XX] = myrank;
    pointToNeighbours();
  } // If 1D
  
  else if (SimPM.ndim==2 || SimPM.ndim==3) {
    for (int i=0;i<SimPM.ndim;i++) {
      LocalRange[i] = SimPM.Range[i];
      LocalXmin[i]  = SimPM.Xmin[i];
      LocalXmax[i]  = SimPM.Xmax[i];
      LocalNG[i]    = SimPM.NG[i];
      LocalNcell    = SimPM.Ncell;
    }
    double maxrange=0.;
    enum axes maxdir=XX;
    int npcounter=1;
    
    // Loop, dividing range in half n times, until npcounter=nproc.
    if (nproc==1) { 
      //cout <<"Only one processor!!! Use serial version...\n";
      //rep.error("Only one processor!!! Use serial version",nproc);
    }
    while (npcounter < nproc) {
      // find longest axis of subrange.
      int i=0, dsrc=-1;
      
      // --- Check if we are doing raytracing with a source at infinity. ---
      if (SimPM.EP.raytracing && SimPM.RS.Nsources>0) {
        //
	// check if we have only one source at infinity, b/c then we decompose to keep
	// rays on one processor all the time.
        //
        bool at_infinity=true;
        for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
          if (SimPM.RS.sources[isrc].type==RT_SRC_DIFFUSE) at_infinity=false;
          if ((SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) &&
              (!SimPM.RS.sources[isrc].at_infinity))       at_infinity=false;
        }
        //
        // Now if at_infinity is still set, we want to get the direction along
        // which we don't want to decompose the domain.
        // Now also checks if there are multiple sources at infinity in
        // different directions, in which case we can't make all of them fully
        // parallel.
        //
        if (at_infinity) {
          for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
            if ((SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) && 
                (SimPM.RS.sources[isrc].at_infinity) ) {
              for (int ii=0;ii<SimPM.ndim;ii++) {
                if (fabs(SimPM.RS.sources[isrc].pos[ii])>1.e99) {
                  if (dsrc==-1) dsrc=ii;
                  // srcs at infinity in multiple directions
                  else if (dsrc!=ii) at_infinity=false;
                }
              } // loop over directions
            }   // if we are at the ionising source.
          }     // loop over sources.
        }       // if at infinity
        
        if (at_infinity) {
          //cout <<"\t\tFound a source at infinity in direction "<<dsrc;
          //cout <<", so not decomposing domain in this direction.\n";
          if (dsrc<0) rep.error("no direction to source at infinity",dsrc);
        }
        else {
          //cout <<"\t\tEither multiple or no sources at infinity.\n";
          dsrc=-1;
        }
      }
      // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
      // DOMAIN DECOMPOSITION.
      dsrc = -1;
      // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
      // DOMAIN DECOMPOSITION.

      // --- end of RT source at infinity bit ---
      
      maxrange=0.; i=0;
      while (i<SimPM.ndim) {
        if (LocalRange[i] > maxrange && dsrc!=i) {
          maxrange=LocalRange[i];
          maxdir = static_cast<axes>(i);
        }
        i++;
      }
      // Half that range and multiply nproc by 2.
      LocalRange[maxdir] /= 2.0;
      npcounter *=2;
    }
    if (npcounter != nproc) rep.error("nproc not a power of 2!",nproc);
    
    //
    // Now we know the range of each subcell, so calculate where I fit into
    // the hierarchy, defined by
    // \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
    // This requires myrank to count from zero!
    //
    for (int i=0;i<SimPM.ndim;i++)
      nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/LocalRange[i]);
    int temp=myrank;
    if (SimPM.ndim==3) {
      ix[ZZ] = temp/nx[XX]/nx[YY];
      temp -= ix[ZZ]*nx[XX]*nx[YY];
    }
    ix[YY] = temp/nx[XX];
    temp -= ix[YY]*nx[XX];
    ix[XX] = temp;
    
    LocalNcell = SimPM.Ncell;
    for (int i=0;i<SimPM.ndim;i++) {
      LocalXmin[i]  = SimPM.Xmin[i] +ix[i]*LocalRange[i];
      LocalXmax[i]  = SimPM.Xmin[i] +(ix[i]+1)*LocalRange[i];
      LocalNG[i]    = SimPM.NG[i]/nx[i];
      LocalNcell   /= nx[i];
      offsets[i]    = ix[i]*LocalNG[i];
    }
    // Set up ngbprocs array to point to neighbouring processors
    pointToNeighbours();    
  } // if 2D or 3D
  
  else rep.error("Bad NDIM in DecomposeDomain",SimPM.ndim);
  
  // Display some debugging info.
  //  if (myrank==0) {
  //    for (int i=0;i<SimPM.ndim;i++) {
  //     cout <<"Sim: idim="<<i<<"  \tRange="<<SimPM.Range[i];
  //     cout <<",\t  xmin,xmax = "<<SimPM.Xmin[i]<<", "<<SimPM.Xmax[i];
  //     cout <<"    \t Ncell = "<<SimPM.Ncell<<"\n";
  //   }
  // }
  // for (int i=0;i<SimPM.ndim;i++) {
  //   cout <<"Proc "<<myrank<<": idim="<<i<<"  \tRange="<<LocalRange[i];
  //   cout <<",\t  xmin,xmax = "<<LocalXmin[i]<<", "<<LocalXmax[i];
  //   cout <<"    \t Ncell = "<<LocalNcell;
  //   cout <<"\t neighbours : "<<ngbprocs[2*i]<<", "<<ngbprocs[2*i+1]<<"\n";
  // }
  // cout << "---ParallelParams::decomposeDomain() Domain decomposition done.\n\n";
  return(0);
}

int ParallelParams::pointToNeighbours()
{
  /** \section PBC Periodic Boundaries
   * Note that if there are periodic boundaries, the pointers to the 
   * wrapped around neighbouring processors are not set here.  They are
   * set when the grid boundaries are set up, in the function
   * UniformGridParallel::BC_setBCtypes
   * I'm not completely happy with this as it makes the 'modular-ness' of 
   * the code weaker, but it is something that will work.  The neighbouring
   * processor list is a global variable, and I don't think it's the worst
   * thing for it to be set in two places.
   * Here any simulation boundaries have neighbouring processor id set to 
   * -999.
   * */
  int nx[SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++) nx[i] =static_cast<int>(ONE_PLUS_EPS*SimPM.Range[i]/LocalRange[i]);
  // Point to neigbours
  if (ngbprocs) {
    //rep.error("neigbours already set up!",ngbprocs);
    delete [] ngbprocs; ngbprocs=0;
  }
  ngbprocs = new int [2*SimPM.ndim];
  if (!ngbprocs) rep.error("Memory Allocation of neighbours.",ngbprocs);
  ngbprocs[XN] = myrank -1;
  ngbprocs[XP] = myrank +1;
  if (GS.equalD(LocalXmin[XX],SimPM.Xmin[XX])) ngbprocs[XN] = -999;
  if (GS.equalD(LocalXmax[XX],SimPM.Xmax[XX])) ngbprocs[XP] = -999;
  if (SimPM.ndim >1) {
    ngbprocs[YN] = myrank -nx[XX];
    ngbprocs[YP] = myrank +nx[XX];
    if (GS.equalD(LocalXmin[YY],SimPM.Xmin[YY])) ngbprocs[YN] = -999;
    if (GS.equalD(LocalXmax[YY],SimPM.Xmax[YY])) ngbprocs[YP] = -999;
  }
  if (SimPM.ndim >2) {
    ngbprocs[ZN] = myrank -nx[XX]*nx[YY];
    ngbprocs[ZP] = myrank +nx[XX]*nx[YY];
    if (GS.equalD(LocalXmin[ZZ],SimPM.Xmin[ZZ])) ngbprocs[ZN] = -999;
    if (GS.equalD(LocalXmax[ZZ],SimPM.Xmax[ZZ])) ngbprocs[ZP] = -999;
  }
  return(0);
}

///
/// Get a list of all abutting domains, including corner/edge intersections.
///
void ParallelParams::get_abutting_domains(std::vector<int> &dl ///< write list to this vector.
					  )
{
  //
  // If we previously generated the list, then just copy elements.
  //
  //if (!full_ngb_list.empty()) {
  //  for (unsigned int v=0; v<full_ngb_list.size(); v++)
  //    dl.push_back(full_ngb_list[v]);
  //  return;
  //}
  full_ngb_list.clear();

  //
  // Otherwise we generate the list and then copy.
  //
  //
  // List ordering is as follows:
  //  XN,XP,YN,YP,ZN,ZP
  //  YNXN, YNXP, YPXN, YPXP,
  //  ZNXN, ZNXP, ZNYN, ZNYP,
  //  ZNYNXN, ZNYNXP, ZNYPXN, ZNYPXP
  //  ZPXN, ZPXP, ZPYN, ZPYP,
  //  ZPYNXN, ZPYNXP, ZPYPXN, ZPYPXP
  //

  //
  // Find which coordinate directions have neighbours and store in a
  // bool array.
  //
  bool d[2*MAX_DIM];
  for (int v=0;v<2*MAX_DIM;v++) d[v] = false;
  enum direction posdir, negdir;
  for (int v=0;v<SimPM.ndim;v++) {
    negdir = static_cast<direction>(2*v);
    posdir = static_cast<direction>(2*v+1);
    if (ngbprocs[negdir]>=0) d[negdir] = true;
    if (ngbprocs[posdir]>=0) d[posdir] = true;
  }

  //
  // Add coordinate directions to list.
  //
  for (int v=0; v<2*SimPM.ndim; v++) {
    if (d[v]) full_ngb_list.push_back(ngbprocs[v]);
  }

  //
  // X-Y plane -- can read off neighbour and add/subtract 1 from it.
  //
  if (d[YN] && d[XN]) full_ngb_list.push_back(ngbprocs[YN]-1);
  if (d[YN] && d[XP]) full_ngb_list.push_back(ngbprocs[YN]+1);
  if (d[YP] && d[XN]) full_ngb_list.push_back(ngbprocs[YP]-1);
  if (d[YP] && d[XP]) full_ngb_list.push_back(ngbprocs[YP]+1);

  //
  // X-Y-ZN plane -- need to calculate rank since offsets are not trivial.
  //
  if (d[ZN]) {
    if (d[XN]) 
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]) +(ix[XX]-1) );
    if (d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]) +(ix[XX]+1) );
    if (d[YN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]) );
    if (d[YP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]) );
    if (d[YN] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]-1) );
    if (d[YN] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]-1) +(ix[XX]+1) );
    if (d[YP] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]-1) );
    if (d[YP] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]-1) +nx[XX]*(ix[YY]+1) +(ix[XX]+1) );
  }
  //
  // X-Y-ZP plane
  //
  if (d[ZP]) {
    if (d[XN]) 
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]  ) +(ix[XX]-1) );
    if (d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]  ) +(ix[XX]+1) );
    if (d[YN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]  ) );
    if (d[YP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]  ) );
    if (d[YN] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]-1) );
    if (d[YN] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]-1) +(ix[XX]+1) );
    if (d[YP] && d[XN])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]-1) );
    if (d[YP] && d[XP])
      full_ngb_list.push_back( nx[XX]*nx[YY]*(ix[ZZ]+1) +nx[XX]*(ix[YY]+1) +(ix[XX]+1) );
  }

  //
  // Copy to return vector and return.
  //
  for (unsigned int v=0; v<full_ngb_list.size(); v++)
    dl.push_back(full_ngb_list[v]);

  return;
}

//
// Returns the ix array for any requested rank.
//
void ParallelParams::get_domain_ix(const int r, ///< rank ix requested for.
				   int *arr     ///< array to put ix into.
				   )
{
  int temp=r;
  if (SimPM.ndim==3) {
    arr[ZZ] = temp/nx[XX]/nx[YY];
    temp -= arr[ZZ]*nx[XX]*nx[YY];
  }
  if (SimPM.ndim>1) {
    arr[YY] = temp/nx[XX];
    temp -= arr[YY]*nx[XX];
  }
  arr[XX] = temp;
  return;
}

#endif //PARALLEL
//------------------------------------------------



//------------------------------------------------
//-------------- JET PARAMETERS ------------------
//------------------------------------------------
JetParams::JetParams()
{
   jetic =0; jetradius = -1;
   jetstate=0; jetstate = new double [MAX_NVAR];
   if (!jetstate) rep.error("Couldn't allocate memory for JP.jetstate[]",jetstate);
   for (int v=0; v<MAX_NVAR; v++) jetstate[v] = -1.e99;
}

JetParams::~JetParams()
{
  if (jetstate) {delete [] jetstate; jetstate=0;}
}
//------------------------------------------------


SimParams::SimParams()
{
//  cout <<"Setting up SimParams class... ";
  gridType = eqntype = solverType = -1;
  ndim = eqnNDim = nvar = -1;
  ntracer = ftr = -1; trtype="BAD-TRACER-TYPE";
  starttime = simtime = finishtime = dt = -1.e99; last_dt = 1.e100;
  timestep = -1; maxtime = false;
  for (int i=0;i<MAX_DIM;i++) {
    NG[i] = -1;
    Range[i] = Xmin[i] = Xmax[i] = -1.e99;
  }
  Ncell = Nbc = -1;
  spOOA = tmOOA = OA1;
  dx = gamma = CFL = etav = -1.e99;
  artviscosity = opfreq = -1;
  typeofip = typeofop = -1;
  typeofbc = "BAD-BC";
  outFileBase = "BAD-FILE";
  op_criterion=0; // default to per n-steps
  next_optime = opfreq_time = 0.0;
  addnoise=0;
  //RefVec=0; RefVec = new double [MAX_NVAR];
  //if (!RefVec) rep.error("Couldn't allocate memory for SimPM.RefVec[]",RefVec);
  for (int v=0; v<MAX_NVAR; v++) RefVec[v] = -1.e99;
  EP.dynamics          = 1;
  EP.raytracing        = 0;
  EP.cooling           = 0;
  EP.chemistry         = 0;
  EP.coll_ionisation   = 0;
  EP.phot_ionisation   = 0;
  EP.rad_recombination = 0;
  EP.update_erg        = 1;  ///< this is effectively a boolean value.
  EP.MP_timestep_limit = false; ///< by default only use hydro limit.

  EP.MinTemperature = 0.0;
  EP.MaxTemperature = 1.0e100;
  
  EP.Helium_MassFrac = 0.2703;
  EP.Metal_MassFrac  = 0.0142;

  RS.Nsources = -1;
  RS.sources.clear();

  STAR.clear();

  min_timestep = 0.0;
  //  cout <<"done!\n";
}

SimParams::~SimParams()
{
  //cout <<"Deleting SimParams class.\n";
  //if (RefVec) {delete [] RefVec; RefVec=0;}
  RS.sources.clear();

  for (size_t v=0; v<STAR.size(); v++) {
    STAR[v].time.clear();
    STAR[v].Log_L.clear();
    STAR[v].Log_T.clear();
    STAR[v].Log_R.clear();
    STAR[v].Log_V.clear();
  }
  STAR.clear();
}

// ----------------------------------------------------
// ----------------- REPORTING CLASS ------------------
// ----------------------------------------------------
reporting::reporting()
{
//  cout <<"Default reporting Constructor. O/P goes to cout/cerr.\n";
}

reporting::~reporting()
{
#if defined (SERIAL)
  cout.rdbuf(saved_buffer_cout);   
  infomsg.close();
  //cout<<"Deleting reporting class, This should be stdout.\n";
  cerr.rdbuf(saved_buffer_cerr);   
  errmsg.close();
  //cout<<"Deleting reporting class, This should be stderr.\n";
#elif defined (PARALLEL)
  if (mpiPM.myrank==0) {
    cout.rdbuf(saved_buffer_cout);   
    infomsg.close();
  }
  else {
    std::cout.clear();
  }
#else
#error "Must define either SERIAL or PARALLEL (reporting::~reporting)"
#endif
}

int reporting::redirect(const string &path)
{
  string temp;

#if defined (SERIAL)
#ifdef TESTING
  cout <<"(reporting::redirect): O/P goes to text files in "<<path<<"\n";
#endif
  temp = path+"errors.txt";
  errmsg.open(temp.c_str(), ios::trunc);
  if(!errmsg.is_open()) {
    cerr<<"Reporting: can't open errors.txt for writing.\n";
    exit(1);
  }
  errmsg.copyfmt( cerr );
  saved_buffer_cerr = cerr.rdbuf();
  cerr.rdbuf( errmsg.rdbuf() );

  temp = path+"info.txt";
  infomsg.open(temp.c_str(), ios::trunc);
  if(!infomsg.is_open()) {
    cerr<<"Reporting: can't open info.txt for writing.\n";
    exit(1);
  }
  infomsg.copyfmt( cout );
  saved_buffer_cout = cout.rdbuf();
  cout.rdbuf( infomsg.rdbuf() );
  
  cout.setf(ios_base::scientific); cout.precision(7);

#elif defined (PARALLEL)
#ifndef TESTING
  //
  // For parallel execution (production runs) we only want a single
  // log file, and errors should be printed to stderr.
  //
  //cout <<"myrank="<<mpiPM.myrank<<"\n";
  if (mpiPM.myrank==0) {
    cout <<"(reporting::redirect): O/P goes to text files in ";
    cout <<path<<"\n";
    cout <<"Note: not redirecting error messages, and suppressing ";
    cout <<"stdout from all processes except myrank=0.\n";
    
    temp = path+"info.txt";
    infomsg.open(temp.c_str(), ios::trunc);
    if(!infomsg.is_open()) {
      cerr<<"Reporting: can't open info.txt for writing.\n"; exit(1);
    }
    infomsg.copyfmt( cout );
    saved_buffer_cout = cout.rdbuf();
    cout.rdbuf( infomsg.rdbuf() );
    cout.setf(ios_base::scientific); cout.precision(7);
  }
  else {
    //saved_buffer_cout = cout.rdbuf(); // <-- save
    //cout.rdbuf (nullstream.rdbuf());  // <-- redirect
    std::cout.setstate(std::ios::failbit) ;
  }
#else
  //
  // for testing we want all processors to have their own log file.
  //
  temp = path+"info.txt";
  infomsg.open(temp.c_str(), ios::trunc);
  if(!infomsg.is_open()) {
    cerr<<"Reporting: can't open info.txt for writing.\n";
    exit(1);
  }
  infomsg.copyfmt( cout );
  saved_buffer_cout = cout.rdbuf();
  cout.rdbuf( infomsg.rdbuf() );
  cout.setf(ios_base::scientific); cout.precision(7);

  temp = path+"errors.txt";
  errmsg.open(temp.c_str(), ios::trunc);
  if(!errmsg.is_open()) {
    cerr<<"Reporting: can't open errors.txt for writing.\n";
    exit(1);
  }
  errmsg.copyfmt( cerr );
  saved_buffer_cerr = cerr.rdbuf();
  cerr.rdbuf( errmsg.rdbuf() );
#endif // TESTING
#else
#error "Must define either SERIAL or PARALLEL (reporting::redirect)"
#endif
  return(0);
}


#ifdef TESTING
class DebugParams dp;
class CommandLineInterface commandline;

DebugParams::DebugParams()
{
//   cout <<"initialising DebugParams...\n";
   initERG=initMMX=initMMY=initMMZ=0.0;
   ergTotChange=mmxTotChange=mmyTotChange=mmzTotChange=0.0;
   c=0;
}

DebugParams::~DebugParams()
{
}

/************************ COMMAND LINE INTERFACE ************************************/
CommandLineInterface::CommandLineInterface() 
{
//  cout <<"Setting up interactive debugging.\n";
}
CommandLineInterface::~CommandLineInterface() 
{
//  cout <<"Destructing interactive debugging.\n";
}

void CommandLineInterface::auto_console(char prompt[])
{
  cout<<"\n Welcome to the AUTOPILOT command line\n";
  cout<<"****************************************\n";
  cout <<prompt<<"\n";
  cout <<"going to first point\n";
  fpt();
  print_cell();
  end_of_col("YP");
  print_cell();
  next_point("YP");
  print_cell();
  next_point("XN");
  print_cell();
  next_point("YN");
  print_cell();
  next_point("XP");
  print_cell();

  cout<<"\n LEAVING the AUTOPILOT command line\n";
  cout<<"****************************************\n";
  return;
}


void CommandLineInterface::console(char prompt[])
{
  //
  // A command line interface which allows investigation of
  //  the grid during runtime.
  //
  fprintf(stderr,"\n Welcome to the command line - type \"help\"\n\n");
  read_history(".cmddemo_CLI.hist");
  //
  // Get a command. Exercise: figure out what the vicious bit of
  //  ptr shit (A) does! :) It's useful!
  //
  char *comm=NULL;
  for (;;) {
    comm=rl_gets(comm,prompt);
    char *c=comm;
    do {
      char comm1[512],*d=comm1;
      *d=0;
      while (*c&&((*(c))!=';')) {*d++=*c++;*(d+1)=0;}   // <--- (A) 
      if (((*c)==';'))  c++;
      if (execute(comm1)) return;
    } while (*c);
  }
}

char * CommandLineInterface::rl_gets(char *line, char *prompt)
{
  //
  // Reads a string, and return a pointer to it.  Returns NULL on EOF.
  //
  if (line) {
    //
    // If the buffer has already been allocated, return the memory
    //      to the free pool.
    //
    free (line);
    line = (char *)NULL;
  }
  //
  // Get a line from the user.
  //
  line=readline(prompt);
  //
  // If the line has any text in it, save it on the history.
  //
  if (line && *line &&line[0]!='q') add_history (line);
  return (line);
}



int CommandLineInterface::execute(char *com)
{
  //
  // Executes the command in com[], Could be smarter than to hardcode the 
  //  command lengths.
  //
  while (((*(com))==' ')) {*com++;} // strip off any leading white space
  if (!(*com)) return 0;            // don't execute NULL commands!!!
  else if  (!strncmp(com,"!"     ,1))  system(com+1);
  else if  (!strncmp(com,"cmd1"  ,4))  cmd1();
  else if  (!strncmp(com,"cmd2"  ,4))  cmd2(com+4);
  else if  (!strncmp(com,"bigcmd",6))  bigcmd();
  else if  (!strncmp(com,"print_cell",10))  print_cell();
  else if  (!strncmp(com,"print_flux",10))  print_flux();
  else if  (!strncmp(com,"next_point",10))  {
    string s=com;
    s=s.substr(11,2); cout <<s<<"\n";
    next_point(s);
  }
  else if  (!strncmp(com,"fpt",3)) fpt();
  else if  (!strncmp(com,"lpt",3)) lpt();
  else if  (!strncmp(com,"end_of_col",10)) {
    string s=com;
    s=s.substr(11,2); cout <<s<<"\n";
    end_of_col(s);
  }
  else if  (!strncmp(com,"help"  ,4)) {
    fprintf(stderr,"\n commands:\n");
    fprintf(stderr,"           cmd1 - demo cmd\n");
    fprintf(stderr,"           cmd2 - demo cmd with args\n");
    fprintf(stderr,"         bigcmd - no particular reason for this\n");
    cerr<<"print_cell - prints info about the current cell.\n";
    cerr<<"print_flux - If in dynamics solver, prints intercell fluxes.\n";
    cerr<<"next_point - move from current cell to next one, in a direction (XY,XN,etc.)\n";
    cerr<<"fpt() - moves dp.c to the first point.\n";
    cerr<<"lpt() - moves dp.c to the last point.\n";
    cerr<<"end_of_col(string) - moves dp.c to the edge of the grid in direction (XN,XP,ZN,ZP,YN,YP).\n";
    fprintf(stderr,"           help - prints this list\n");
    fprintf(stderr,"  q, quit, exit - end the program\n");
    fprintf(stderr,"\n leading \"!\" sends cmd to shell\n\n");
  }
  else if  (!strncmp(com,"cont"  ,4))  {
    fprintf(stderr," Continuing with code execution.\n");
    write_history(".cmddemo_CLI.hist");
    history_truncate_file (".cmddemo_CLI.hist",50);
    return 1;
  }
  else if ((!strncmp(com,"q"       ,1))||
           (!strncmp(com,"quit"    ,4))||
           (!strncmp(com,"exit"    ,4))) {
    write_history(".cmddemo_CLI.hist");
    history_truncate_file (".cmddemo_CLI.hist",50);
    rep.error("Quitting at user request",com);
  }
//   else if  (!strncmp(com,"shutdown",8)) system("halt -p"); Steph wouldn't like this...
  else {
    fprintf(stderr," Unknown command: %s\n\n",com);
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                        Commands are below.                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////
void CommandLineInterface::cmd1() 
{
  fprintf(stderr," This is command 1\n");
  double x[MAX_DIM];
  CI.get_dpos(dp.c,x);
  rep.printVec("position",x,2);
}

void CommandLineInterface::cmd2(char *args) 
{
  fprintf(stderr," This is command 2, args: \"%s\"\n",args);
  rep.printVec("P ",dp.c->P, 5);
  rep.printVec("Ph",dp.c->P, 5);
  
}

void CommandLineInterface::bigcmd() 
{
  fprintf(stderr," This is a big command (apparently).\n");
}

void CommandLineInterface::print_cell() 
{
  grid->PrintCell(dp.c);
  cout <<"\n";
  return;
}

void CommandLineInterface::next_point(string s) 
{
  enum direction dir = parse_dir(s);
  if (dir==NO) return;
  
  if (!(dp.c->ngb[dir])) {
    cout <<"no neighbour in direction "<<s<<"; try again.\n"; return;
  }
  else {
    cout <<"moving from cell "<<dp.c->id;
    dp.c = dp.c->ngb[dir];
    cout <<" to cell "<<dp.c->id<<"\n";
  }  
  return;
}

void CommandLineInterface::print_flux() 
{
  rep.printVec("left ",dp.vec2,SimPM.nvar);
  rep.printVec("right",dp.vec3,SimPM.nvar);
  rep.printVec("pstar",dp.vec1,SimPM.nvar);
  rep.printVec("flux ",dp.vec4,SimPM.nvar);
  return;
}

void CommandLineInterface::fpt()
{
  dp.c = grid->FirstPt();
  return;
}

void CommandLineInterface::lpt()
{
  dp.c = grid->LastPt();
  return;
}

void CommandLineInterface::end_of_col(const string s)
{
  if (!(dp.c)) {
    cout <<"dp.c is null pointer! can't get to end of column.\n";
    return;
  }
  enum direction dir = parse_dir(s);
  if (dir==NO) return;

  int ct=0;
  while (grid->NextPt(dp.c,dir)->isgd) {
    dp.c=grid->NextPt(dp.c,dir);
    ct++;
  }
  cout <<"moved "<<ct<<" cells in direction "<<s<<"\n";
  return;
}

enum direction CommandLineInterface::parse_dir(const string s)
{
  enum direction dir=NO;
  if      (s=="XP") dir=XP;
  else if (s=="XN") dir=XN;
  else if (s=="YP") dir=YP;
  else if (s=="YN") dir=YN;
  else if (s=="ZP") dir=ZP;
  else if (s=="ZN") dir=ZN;
  else {cout <<"bad direction string: "<<s<<", try again.\n";}
  return dir;
}

#endif //TESTING


