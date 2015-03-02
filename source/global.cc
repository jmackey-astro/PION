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
/// - 2015.01.08 JM: moved grid base class to grid/grid_base_class.h
/// - 2015.01.12 JM: moved reporting to tools/reporting.h.
///    Moved commandline interface to tools/command_line_interface.h.
///    Moved memory management to tools/mem_manage.h.
/// - 2015.01.26 JM: Removed ParallelParams class
///    (to sim_control_MPI.cpp), and COMMs class stuff.

#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <string>
#include <stdexcept>
using namespace std;

#include "global.h"
#include "tools/reporting.h"
#include "sim_constants.h"
#include "sim_params.h"


class SimParams SimPM;
class JetParams JP;
class Units uc;
class GeneralStuff GS;
//class stellar_wind SW;
struct stellarwind_list SWP;


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
//#ifdef PARALLEL
//#if   defined USE_MPI
//class comms_base *COMM = new comm_mpi ();
//#elif defined USE_FILE_COMMS
//class comms_base *COMM = new comm_files ();
//#else
//#error "MUST DEFINE EITHER USE_MPI or USE_FILE_COMMS"
//#endif

//#endif // PARALLEL

/************************ MULTI-PROCESS COMMS *************************/

/************************ MICROPHYSICS ***************************/
class MicroPhysicsBase *MP=0;
/************************ MICROPHYSICS ***************************/

/************************  RAYTRACING  ***************************/
class RayTracingBase *RT=0;
/************************  RAYTRACING  ***************************/



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
}

GeneralStuff::~GeneralStuff() 
{
  //cout <<"timers.size: "<<timers.size()<<"\n";
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



