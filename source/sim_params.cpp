///
/// \file sim_params.h
/// \author Jonathan Mackey
///
/// Structure with simulation parameters.
///
/// Modifications:
/// - 2015.07.06 JM: moved from global.cc
/// - 2017.11.22 JM: worked on boundary conditions string


#include "sim_params.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
using namespace std;

//class SimParams SimPM;
class JetParams JP;
class Units uc;
struct stellarwind_list SWP;



//------------------------------------------------
//-------------- JET PARAMETERS ------------------
//------------------------------------------------
JetParams::JetParams()
{
   jetic =0; jetradius = -1;
   jetstate=0; jetstate = new pion_flt [MAX_NVAR];
   if (!jetstate) rep.error("Couldn't allocate memory for JP.jetstate[]",jetstate);
   for (int v=0; v<MAX_NVAR; v++) jetstate[v] = -1.0e30;
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

  ntracer = ftr = -1;
  chem_code = "BAD_CHEM_CODE";
  TRTYPE = "";
  tracers=0;

  starttime = simtime = finishtime = dt = -1.e99; last_dt = 1.e100;
  timestep = -1; maxtime = false;
  for (int i=0;i<MAX_DIM;i++) {
    NG[i] = -1;
    Range[i] = Xmin[i] = Xmax[i] = -1.e99;
  }
  grid_nlevels = 1;
  levels.resize(1);
  Ncell = Nbc = -1;
  spOOA = tmOOA = OA1;
  dx = gamma = CFL = etav = -1.e99;
  artviscosity = opfreq = -1;
  typeofip = typeofop = -1;

  BC_XN = "";
  BC_XP = "";
  BC_YN = "";
  BC_YP = "";
  BC_ZN = "";
  BC_ZP = "";
  BC_INT = 0;
  BC_STRING = "";

  outFileBase = "BAD-FILE";
  op_criterion=0; // default to per n-steps
  next_optime = opfreq_time = 0.0;
  addnoise=0;
  for (int v=0; v<MAX_NVAR; v++) RefVec[v] = -1.e30;
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
  
  EP.H_MassFrac = 0.73;
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

  BC_INT = mem.myfree(BC_INT); BC_INT=0;

  tracers=mem.myfree(tracers);
}



