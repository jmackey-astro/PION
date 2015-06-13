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
/// - 2015.03.10 JM: Removed GeneralStuff class, to constants.h and
///    other files.

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
#include "constants.h"
#include "sim_params.h"


class SimParams SimPM;
class JetParams JP;
class Units uc;
struct stellarwind_list SWP;



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



