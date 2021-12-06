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
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include <sstream>

using namespace std;

// class SimParams SimPM;
class JetParams JP;
class Units uc;
struct stellarwind_list SWP;

//------------------------------------------------
//-------------- JET PARAMETERS ------------------
//------------------------------------------------
JetParams::JetParams()
{
  jetic     = 0;
  jetradius = -1;
  jetstate  = 0;
  jetstate  = new pion_flt[MAX_NVAR];
  if (!jetstate)
    spdlog::error(
        "{}: {}", "Couldn't allocate memory for JP.jetstate[]",
        fmt::ptr(jetstate));
  for (int v = 0; v < MAX_NVAR; v++)
    jetstate[v] = -1.0e30;
}

JetParams::~JetParams()
{
  if (jetstate) {
    delete[] jetstate;
    jetstate = 0;
  }
}
//------------------------------------------------

SimParams::SimParams()
{
  //  cout <<"Setting up SimParams class... ";
  gridType = eqntype = solverType = -1;
  ndim = eqnNDim = nvar = -1;

  ntracer = ftr = -1;
  chem_code     = "BAD_CHEM_CODE";
  TRTYPE        = ""; /* TODO */
  tracers       = 0;

  starttime = simtime = finishtime = dt = -1.e99;
  last_dt                               = 1.e100;
  timestep                              = -1;
  maxtime                               = false;
  for (int i = 0; i < MAX_DIM; i++) {
    NG[i]    = -1;
    Range[i] = Xmin[i] = Xmax[i] = -1.e99;
    NG_centre[i]                 = 0.0;
    NG_refine[i]                 = 0;
  }
  grid_nlevels = 1;
  levels.resize(grid_nlevels);
  Ncell = Nbc = -1;
  spOOA = tmOOA = OA1;
  dx = gamma = CFL = etav = -1.e99;
  artviscosity = opfreq = -1;
  typeofip = typeofop = -1;

  BC_XN     = "";
  BC_XP     = "";
  BC_YN     = "";
  BC_YP     = "";
  BC_ZN     = "";
  BC_ZP     = "";
  BC_INT    = 0;
  BC_STRING = "";

  outFileBase  = "BAD-FILE";
  op_criterion = 0;  // default to per n-steps
  next_optime = opfreq_time = 0.0;
  addnoise                  = 0;
  for (int v = 0; v < MAX_NVAR; v++)
    RefVec[v] = -1.e30;
  EP.dynamics          = 1;
  EP.raytracing        = 0;
  EP.cooling           = 0;
  EP.chemistry         = 0;
  EP.coll_ionisation   = 0;
  EP.phot_ionisation   = 0;
  EP.rad_recombination = 0;
  EP.update_erg        = 1;      ///< this is effectively a boolean value.
  EP.MP_timestep_limit = false;  ///< by default only use hydro limit.

  EP.MinTemperature = 0.0;
  EP.MaxTemperature = 1.0e100;

  EP.H_MassFrac      = 0.73;
  EP.Helium_MassFrac = 0.2703;
  EP.Metal_MassFrac  = 0.0142;

  EP.sat_thermal_cond = 0;
  EP.tc_strength      = 1.0;

  RS.Nsources = -1;
  RS.sources.clear();

  STAR.clear();

  min_timestep = 0.0;
  //  cout <<"done!\n";
}

// ##################################################################
// ##################################################################


SimParams::SimParams(const std::string pfile)
{
  int err = read_gridparams(pfile);
  if (err) spdlog::error("{}: {}", "Read Grid Params Error", err);
  last_dt = 1.e100;
  maxtime = false;
  levels.resize(grid_nlevels);
}

// ##################################################################
// ##################################################################


SimParams::~SimParams()
{
  // cout <<"Deleting SimParams class.\n";
  // if (RefVec) {delete [] RefVec; RefVec=0;}
  RS.sources.clear();

  for (size_t v = 0; v < STAR.size(); v++) {
    STAR[v].time.clear();
    STAR[v].Log_L.clear();
    STAR[v].Log_T.clear();
    STAR[v].Log_R.clear();
    STAR[v].Log_V.clear();
  }
  STAR.clear();

  if (BC_Nint) {
    BC_INT = mem.myfree(BC_INT);
    BC_INT = 0;
  }

  if (ntracer) tracers = mem.myfree(tracers);

#ifndef NDEBUG
  spdlog::info("SimParams Destructor: done");
#endif
}

// ##################################################################
// ##################################################################

int SimParams::read_gridparams(const string pfile  ///< paramfile.
)
{
  int err = 0;

  err += rp.read_paramfile(pfile);
  if (err) spdlog::error("{}: {}", "Error reading parameterfile", pfile);
  //   rp.write_out_parameters();

  gridType = 1;  // uniform grid. only option for now.

  // Sytem of Equations to solve:
  string str = rp.find_parameter("eqn");
  if (str == "hd" || str == "euler") {
    // inviscid euler equations.
    eqntype = EQEUL;
    nvar    = 5;
  }
  else if (str == "hd-Eint" || str == "EulerEint") {
    // inviscid Euler equations, integrating internal energy.
    eqntype = EQEUL_EINT;
    nvar    = 5;
  }
  else if (str == "iso-hd" || str == "iso-eul") {
    //
    // inviscid isothermal hydro equations.
    //
    eqntype = EQEUL_ISO;
    nvar    = 5;
  }
  else if (str == "i-mhd") {
    // ideal mhd with no divB cleaning (not great!)
    eqntype = EQMHD;
    nvar    = 8;
  }
  else if (str == "mhd" || str == "glm-mhd") {
    // dedner's glm method (GOOD!)
    eqntype = EQGLM;
    nvar    = 9;
  }
  else if (str == "fieldCD" || str == "fcd-mhd") {
    // toth's field CD method (BAD! and not working!)
    eqntype = EQFCD;
    nvar    = 8;
  }
  else
    spdlog::error("{}: {}", "Don't know what these equations are:", str);

  str = rp.find_parameter("solver");
  if (str == "")
    solverType = 4;  // Roe-type solver (as opposed to lax-friedrichs).
  else if (str == "LF" || str == "LAX-FRIEDRICHS" || str == "0")
    solverType = 0;
  else if (str == "LINEAR" || str == "linear" || str == "1")
    solverType = 1;
  else if (str == "EXACT" || str == "exact" || str == "2")
    solverType = 2;
  else if (str == "HYBRID" || str == "hybrid" || str == "3")
    solverType = 3;
  else if (str == "ROE" || str == "roe" || str == "4")
    solverType = 4;
  else if (str == "ROEPV" || str == "roepv" || str == "5")
    solverType = 5;
  else if (str == "FVS" || str == "fvs" || str == "6")
    solverType = 6;
  else if (str == "HLLD" || str == "hlld" || str == "7")
    solverType = 7;
  else if (str == "HLL" || str == "hll" || str == "8")
    solverType = 8;
  else
    spdlog::error("{}: {}", "what kind of solver is this???", str);

  // dimensionality of grid.
  str = rp.find_parameter("ndim");
  if (str == "")
    spdlog::error("{}: {}", "Bad ndim in pfile", str);
  else
    ndim = atoi(str.c_str());

  if (isnan(ndim)) spdlog::error("{}: {}", "ndim is not a number", ndim);
  eqnNDim = 3;

  //
  // tracer variables: number of tracers, and then each one gets a
  // string name, stored in tracers
  //
  str = rp.find_parameter("ntracer");
  if (str == "") {
#ifndef NDEBUG
    spdlog::debug("Not using tracer variablesn");
#endif
    ntracer = 0;
    ftr     = nvar;
  }
  else if (isfinite(ntracer = atoi(str.c_str()))) {
#ifndef NDEBUG
    spdlog::debug("using {} passive tracer variables", ntracer);
#endif
    ftr = nvar;
    nvar += ntracer;
    //
    // Get what type of chemistry we are doing:
    //
    chem_code = rp.find_parameter("chem_code");

    //
    // Each tracer has its own parameter, called Tracer000, Tracer001,
    // etc., so setup a string for each in the tracers array.
    //
    tracers = mem.myalloc(tracers, ntracer);
    for (int i = 0; i < ntracer; i++) {
      ostringstream temp;
      temp << "Tracer";
      temp.width(3);
      temp.fill('0');
      temp << i;
      tracers[i] = rp.find_parameter(temp.str());
#ifndef NDEBUG
      spdlog::debug("using tracer(s) described as {}", tracers[i]);
#endif
      if (tracers[i] == "") {
        spdlog::error("{}: {}", "Can't find tracer name for number", i);
      }
    }
  }
  else
    spdlog::error("{}: {}", "number of tracers is not finite!", ntracer);

  // coordinate system.
  str = rp.find_parameter("coordinates");
  if (str == "") {
    spdlog::debug(
        "no coordinate system specified.  using cartesian as default");
    coord_sys = COORD_CRT;
  }
  else if (str == "cartesian" || str == "cart")
    coord_sys = COORD_CRT;
  else if (str == "cylindrical" || str == "cyl")
    coord_sys = COORD_CYL;
  else if (str == "spherical" || str == "sph")
    coord_sys = COORD_SPH;
  else
    spdlog::error(
        "{}: {}", "Don't recognise coordinate system in GetParameters", str);

  // Domain of grid, and number of points.  These must all be present in the
  // pfile.
  string seek;
  seek = "NGridX";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  NG[XX] = atoi(str.c_str());
  Ncell  = NG[XX];

  seek = "Xmin";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  Xmin[XX] = atof(str.c_str());

  seek = "Xmax";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  Xmax[XX]  = atof(str.c_str());
  Range[XX] = Xmax[XX] - Xmin[XX];

  if (ndim > 1) {
    seek = "NGridY";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    NG[YY] = atoi(str.c_str());
    Ncell *= NG[YY];

    seek = "Ymin";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    Xmin[YY] = atof(str.c_str());

    seek = "Ymax";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    Xmax[YY]  = atof(str.c_str());
    Range[YY] = Xmax[YY] - Xmin[YY];
  }
  else {
    NG[YY]   = 1;
    Xmin[YY] = Xmax[YY] = Range[YY] = 0.;
  }

  if (ndim > 2) {
    seek = "NGridZ";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    NG[ZZ] = atoi(str.c_str());
    Ncell *= NG[ZZ];

    seek = "Zmin";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    Xmin[ZZ] = atof(str.c_str());

    seek = "Zmax";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    Xmax[ZZ]  = atof(str.c_str());
    Range[ZZ] = Xmax[ZZ] - Xmin[ZZ];
  }
  else {
    NG[ZZ]   = 1;
    Xmin[ZZ] = Xmax[ZZ] = Range[ZZ] = 0.;
  }

  dx = Range[XX] / NG[XX];

  if (ndim > 1) {
    if (fabs(Range[1] / (NG[1]) / dx - 1.) > 100. * MACHINEACCURACY) {
      spdlog::error(
          "{}: {}",
          "SimParams: Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          fabs(Range[1] / (NG[1]) / dx - 1.));
    }
  }
  if (ndim > 2) {
    if (fabs(Range[2] / (NG[2]) / dx - 1.) > 100. * MACHINEACCURACY) {
      spdlog::error(
          "{}: {}",
          "SimParams: Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          fabs(Range[2] / (NG[2]) / dx - 1.));
    }
  }

  // Nested grid info (all optional for now)
  seek = "grid_nlevels";
  str  = rp.find_parameter(seek);
  if (str == "")
    grid_nlevels = 1;
  else
    grid_nlevels = atoi(str.c_str());

  seek = "grid_aspect_ratio_XX";
  str  = rp.find_parameter(seek);
  if (str == "")
    grid_aspect_ratio[XX] = 1;
  else
    grid_aspect_ratio[XX] = atoi(str.c_str());

  seek = "grid_aspect_ratio_YY";
  str  = rp.find_parameter(seek);
  if (str == "")
    grid_aspect_ratio[YY] = 1;
  else
    grid_aspect_ratio[YY] = atoi(str.c_str());

  seek = "grid_aspect_ratio_ZZ";
  str  = rp.find_parameter(seek);
  if (str == "")
    grid_aspect_ratio[ZZ] = 1;
  else
    grid_aspect_ratio[ZZ] = atoi(str.c_str());

  seek = "NG_centre_XX";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_centre[XX] = 0.0;
  else
    NG_centre[XX] = atof(str.c_str());

  seek = "NG_centre_YY";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_centre[YY] = 0.0;
  else
    NG_centre[YY] = atof(str.c_str());

  seek = "NG_centre_ZZ";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_centre[ZZ] = 0.0;
  else
    NG_centre[ZZ] = atof(str.c_str());

  seek = "NG_refine_XX";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_refine[XX] = 1;
  else
    NG_refine[XX] = atoi(str.c_str());

  seek = "NG_refine_YY";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_refine[YY] = 1;
  else
    NG_refine[YY] = atoi(str.c_str());

  seek = "NG_refine_ZZ";
  str  = rp.find_parameter(seek);
  if (str == "")
    NG_refine[ZZ] = 1;
  else
    NG_refine[ZZ] = atoi(str.c_str());

  // output info
  str = rp.find_parameter("OutputPath");
  if (str == "") spdlog::error("{}: {}", "outputpath", str);
  ostringstream temp;
  temp << str;
  str = rp.find_parameter("OutputFile");
  if (str == "") spdlog::error("{}: {}", "outputfile", str);
  temp << str;
  outFileBase = temp.str();
  str         = rp.find_parameter("OutputFileType");

  if (str == "text")
    typeofop = 1;
  else if (str == "fits")
    typeofop = 2;
  else if (str == "ftab") {  // fits table (not implemented yet.
    typeofop = 3;
    spdlog::error(
        "{}: {}", "Fits Table not allowed as o/p format.. not implemented",
        str);
  }
  else if (str == "both")
    typeofop = 4;  // Both means fits and text.
  else if (str == "silo")
    typeofop = 5;  // Silo output.
  else
    spdlog::error("{}: {}", "Error, bad type of outputfile in pfile.", str);

  seek = "OutputFrequency";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  opfreq = atoi(str.c_str());
  // cout <<"\tOutFile: "<<outFileBase<< ".xxx\t Type="<<typeofop;
  // cout <<" every "<<opfreq<<" timesteps."<<"\n";

  // Optional output-by-years parameters:
  seek = "OutputCriterion";
  str  = rp.find_parameter(seek);
  if (str == "") {
    op_criterion = 0;
  }
  else
    op_criterion = atoi(str.c_str());
  // cout <<"\top_criterion = "<<op_criterion<<"\n";

  seek = "OPfreqTime";
  str  = rp.find_parameter(seek);
  if (str == "") {
    opfreq_time = 0.0;
    next_optime = 0.0;
  }
  else {
    opfreq_time = atof(str.c_str());
    next_optime = simtime + opfreq_time;
    double tmp =
        ((simtime / opfreq_time) - static_cast<int>(simtime / opfreq_time))
        * opfreq_time;
    next_optime -= tmp;
  }
  // cout <<" and opfreq="<<opfreq_time<<" code time units.\n";

  //
  // Boundary conditions.
  // First do the edges of the domain:
  //
  BC_XN = rp.find_parameter("BC_XN");
  BC_XP = rp.find_parameter("BC_XP");
  if (ndim > 1) {
    BC_YN = rp.find_parameter("BC_YN");
    BC_YP = rp.find_parameter("BC_YP");
  }
  else {
    BC_YN = "NONE";
    BC_YP = "NONE";
  }
  if (ndim > 2) {
    BC_ZN = rp.find_parameter("BC_ZN");
    BC_ZP = rp.find_parameter("BC_ZP");
  }
  else {
    BC_ZN = "NONE";
    BC_ZP = "NONE";
  }
  //
  // Now do the internal boundaries (if any).  Seek the string
  // for a given internal boundary, and add to a vector until
  // no more are found.
  //
  str = rp.find_parameter("BC_Ninternal");
  if (str == "")
    BC_Nint = 0;
  else
    BC_Nint = atoi(str.c_str());
  if (BC_Nint > 0) {
    BC_INT = mem.myalloc(BC_INT, BC_Nint);
  }
  for (int v = 0; v < BC_Nint; v++) {
    ostringstream intbc;
    intbc.str("");
    intbc << "BC_INTERNAL_";
    intbc.width(3);
    intbc.fill('0');
    intbc << v;
    // cout <<"Looking for internal boundary: "<<intbc.str();
    str = rp.find_parameter(intbc.str());
    // cout <<"   Found: "<<str<<"\n";
    BC_INT[v] = str;
  }
  // rep.printVec("BC_INT",BC_INT,BC_Nint);

  Nbc = -1;  // Set it to negative so I know it's not set.

  // Timing
  seek = "StartTime";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  starttime = atof(str.c_str());
  simtime   = starttime;
  seek      = "FinishTime";
  str       = rp.find_parameter(seek);
  if (str == "")
    finishtime = -1.;
  else
    finishtime = atof(str.c_str());
  timestep = 0;  // Counter for what timestep we are on.

  seek = "OrderOfAccSpace";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  spOOA = atoi(str.c_str());
  seek  = "OrderOfAccTime";
  str   = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  tmOOA = atoi(str.c_str());

  // Physics
  seek = "GAMMA";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  gamma = atof(str.c_str());
  seek  = "CFL";
  str   = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  CFL = atof(str.c_str());

  seek = "ArtificialViscosity";
  str  = rp.find_parameter(seek);
  if (str == "") spdlog::error("{}: {}", "param not found", seek);
  if ((artviscosity = atoi(str.c_str())) == 0) {
    etav = 0.;
  }
  else if (artviscosity == 1 || artviscosity == 4) {
    seek = "EtaViscosity";
    str  = rp.find_parameter(seek);
    if (str == "") spdlog::error("{}: {}", "param not found", seek);
    etav = atof(str.c_str());
  }
  else if (artviscosity == 3) {
    // using H-correction.
    etav = 0.1;
  }
  else
    spdlog::error("{}: {}", "\tUnknown viscosity requested... fix me.", str);
  // cout <<"\tArtificial Viscosity: eta="<<etav<<"\n";
  // Which Physics
  err += read_extra_physics();
  if (err) spdlog::error("{}: {}", "read_extra_physics", err);

  // Raytracing
  err += read_radsources();
  if (err) spdlog::error("{}: {}", "read_radsources", err);

  //
  // STELLAR WINDS???
  //
  seek = "WIND_NSRC";
  str  = rp.find_parameter(seek);
  if (str == "") {
    spdlog::warn("\tWIND_NSRC: param not found, setting to 0");
    SWP.Nsources = 0;
  }
  else {
    SWP.Nsources = atoi(str.c_str());
    // cout <<"\tWIND_NSRC: got "<<SWP.Nsources<<" sources.\n";
    if (SWP.Nsources > 0) {
      err += read_wind_sources();
      if (err) spdlog::error("{}: {}", "read_wind_sources", err);
    }
  }

  //
  // Jets?
  //
  seek = "N_JET";
  str  = rp.find_parameter(seek);
  if (str == "") {
    spdlog::warn("\tN_JET: param not found, setting to 0");
    JP.jetic = 0;
  }
  else {
    JP.jetic = atoi(str.c_str());
    if (JP.jetic) {
      err += read_jet_params(JP);
      if (0 != err)
        spdlog::error("{}: Expected {} but got {}", "read_jet_params", 0, err);
    }
  }

  //
  // Reference State Vector.  Only look for the first nvar elements and set
  // the rest to zero.
  //
  if ((str = rp.find_parameter("refvec0")) != "") {
    for (int v = 0; v < nvar; v++) {
      ostringstream temp2;
      temp2 << "refvec" << v;
      str = rp.find_parameter(temp2.str());
      if (str != "") {
        // cout <<"v="<<v<<" MAX_NVAR="<<MAX_NVAR<<"\n";
        RefVec[v] = atof(str.c_str());
      }
      else {
        spdlog::debug("no refvec[{}], setting to 1", v);
        RefVec[v] = 1.;
      }
    }
    for (int v = nvar; v < MAX_NVAR; v++)
      RefVec[v] = 1.0;
  }
  else {
    for (int v = 0; v < MAX_NVAR; v++)
      RefVec[v] = 1.;
  }

  // Initialise gridparams so that we assume we are not doing a jet sim.
  JP.jetic = 0;

  // Units Used
  err += read_units();
  if (err) spdlog::error("{}: {}", "read_units", err);

  // seek=""; str=rp.find_parameter(seek); if (str=="") spdlog::error("{}: {}",
  // "param not found",seek);
  return err;
}

int SimParams::read_radsources()
{
  string a;
  if ((a = rp.find_parameter("RT_Nsources")) != "")
    RS.Nsources = atoi(a.c_str());
  else {
    // cout <<"no nsources in pfile... assuming there are no RT sources.\n";
    RS.Nsources = 0;
  }

  for (int i = 0; i < RS.Nsources; i++) {
    struct rad_src_info rs_temp;
    ostringstream temp2;

    for (int v = 0; v < ndim; v++) {
      temp2.str("");
      temp2 << "RT_position_" << i << "_" << v;
      if ((a = rp.find_parameter(temp2.str())) != "")
        rs_temp.pos[v] = atof(a.c_str());
      else
        spdlog::error("{}: {}", "no src position in pfile", temp2.str());
    }
    for (int v = ndim; v < MAX_DIM; v++)
      rs_temp.pos[v] = 0.0;

    temp2.str("");
    temp2 << "RT_strength_" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.strength = atof(a.c_str());
    else
      spdlog::error("{}: {}", "no src strength in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_Rstar____" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.Rstar = atof(a.c_str());
    else {
      spdlog::debug(
          "{}: parameter not found in pfile.  Setting to -1.0e200",
          temp2.str());
      rs_temp.Rstar = -1.0e200;
    }

    temp2.str("");
    temp2 << "RT_Tstar____" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.Tstar = atof(a.c_str());
    else {
      spdlog::debug(
          "{}: parameter not found in pfile.  Setting to -1.0e200",
          temp2.str());
      rs_temp.Tstar = -1.0e200;
    }

    rs_temp.id = i;

    temp2.str("");
    temp2 << "RT_src_type_" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.type = atoi(a.c_str());
    else
      spdlog::error("{}: {}", "no src type in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_update___" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.update = atoi(a.c_str());
    else {
      // spdlog::error("{}: {}", "no src update-method in pfile",temp2.str());
      spdlog::error(
          "*** no src update-method in pfile: {}... using C2Ray method",
          temp2.str());
      rs_temp.update = RT_UPDATE_IMPLICIT;
    }

    temp2.str("");
    temp2 << "RT_at_infty_" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.at_infinity = atoi(a.c_str());
    else
      spdlog::error("{}: {}", "no src at-infinity in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_effect___" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.effect = atoi(a.c_str());
    else {
      // spdlog::error("{}: {}", "no src update-method in pfile",temp2.str());
      spdlog::error(
          "*** no src effect param in pfile: {}... using monochromatic photoionisation",
          temp2.str());
      rs_temp.effect = RT_EFFECT_PION_MONO;
    }

    temp2.str("");
    temp2 << "RT_Tau_src__" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.opacity_src = atoi(a.c_str());
    else {
      // spdlog::error("{}: {}", "no src update-method in pfile",temp2.str());
      spdlog::error(
          "*** no opacity src param in pfile: {}... using neutral hydrogen",
          temp2.str());
      rs_temp.opacity_src = RT_OPACITY_MINUS;
    }

    temp2.str("");
    temp2 << "RT_Tau_var__" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.opacity_var = atoi(a.c_str());
    else {
      // spdlog::error("{}: {}", "no src update-method in pfile",temp2.str());
      spdlog::error(
          "*** no opacity variable index in pfile: {}... using First tracer",
          temp2.str());
      rs_temp.opacity_var = ftr;
    }

    temp2.str("");
    temp2 << "RT_EVO_FILE_" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.EvoFile = a;
    else {
      spdlog::error(
          "*** no RS Evolution File in pfile: {}... using NOFILE, i.e. constant source",
          temp2.str());
      rs_temp.EvoFile = "NOFILE";
    }

    temp2.str("");
    temp2 << "RT_Nbins____" << i;
    if ((a = rp.find_parameter(temp2.str())) != "")
      rs_temp.NTau = atoi(a.c_str());
    else {
      spdlog::error(
          "*** no NTau in parameter file: {}... using 1", temp2.str());
      rs_temp.NTau = 1;
    }

    //
    // Add source i to RS list.
    //
    RS.sources.push_back(rs_temp);

    spdlog::debug(
        "\tRT-src: i={}: strength={} type={} at_inf={}, update={} and pos[]=",
        i, rs_temp.strength, rs_temp.type, rs_temp.at_infinity, rs_temp.update);
    spdlog::debug(" : {}", rs_temp.pos);

    // cout <<"\tRT-src: i="<<i<<": strength=";
    // cout <<RS.sources[i].strength;
    // cout <<" type="<<RS.sources[i].type<<" at_inf=";
    // cout <<RS.sources[i].at_infinity;
    // cout <<", update="<<RS.sources[i].update;
    // cout <<" and pos[]=";
    // rep.printVec("",RS.sources[i].pos,ndim);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int SimParams::read_wind_sources()
{
  if (SWP.Nsources == 0) return 2;
  int err = 0;

  for (int i = 0; i < SWP.Nsources; i++) {
    // cout <<"\tREADING WIND SOURCE "<<i<<"\n";
    string a;
    ostringstream temp;
    double Mdot = 0.0, posn[MAX_DIM], vel[MAX_DIM], Vinf = 0.0, Vrot = 0.0,
           Tw = 0.0, Rstar = 0.0, Bstar = 0.0, trcr[MAX_NVAR], xi = 0.0,
           rad = 0.0;
    int type   = 0;
    string evofile;
    int enhance_mdot = 0;
    double time_offset, update_freq, time_scalefac, eccentricity, OrbPeriod,
        PeriastronX, PeriastronY;

    for (int v = 0; v < ndim; v++) {
      temp.str("");
      temp << "WIND_" << i << "_pos" << v;
      if ((a = rp.find_parameter(temp.str())) != "")
        posn[v] = atof(a.c_str());
      else
        spdlog::error("{}: {}", "no src position in pfile", temp.str());
    }
    for (int v = ndim; v < MAX_DIM; v++)
      posn[v] = 0.0;

    for (int v = 0; v < ndim; v++) {
      temp.str("");
      temp << "WIND_" << i << "_velocity" << v;
      if ((a = rp.find_parameter(temp.str())) != "")
        vel[v] = atof(a.c_str());
      else
        spdlog::error("{}: {}", "no src velocity in pfile", temp.str());
    }
    for (int v = ndim; v < MAX_DIM; v++)
      posn[v] = 0.0;

    temp.str("");
    temp << "WIND_" << i << "_radius";
    if ((a = rp.find_parameter(temp.str())) != "")
      rad = atof(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_type";
    if ((a = rp.find_parameter(temp.str())) != "")
      type = atoi(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());
    // cout <<"\t*******wind type="<<type<<"\n";

    temp.str("");
    temp << "WIND_" << i << "_mdot";
    if ((a = rp.find_parameter(temp.str())) != "")
      // convert to g/s from Msun/yr
      Mdot = atof(a.c_str()) * pconst.Msun() / pconst.year();
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_vinf";
    if ((a = rp.find_parameter(temp.str())) != "")
      Vinf = atof(a.c_str()) * 1.0e5;  // convert from km/s to cm/s
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_vrot";
    if ((a = rp.find_parameter(temp.str())) != "")
      Vrot = atof(a.c_str()) * 1.0e5;  // convert from km/s to cm/s
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_temp";
    if ((a = rp.find_parameter(temp.str())) != "")
      Tw = atof(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_Rstr";
    if ((a = rp.find_parameter(temp.str())) != "")
      Rstar = atof(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_Bsrf";
    if ((a = rp.find_parameter(temp.str())) != "")
      Bstar = atof(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());

    for (int v = 0; v < ntracer; v++) {
      temp.str("");
      temp << "WIND_" << i << "_TR" << v;
      if ((a = rp.find_parameter(temp.str())) != "")
        trcr[v] = atof(a.c_str());
      else
        spdlog::error("{}: {}", "param not found in pfile", temp.str());
    }
    for (int v = ntracer; v < MAX_NVAR; v++)
      trcr[v] = 0.0;

    //
    // new stuff for evolving winds:
    //
    temp.str("");
    temp << "WIND_" << i << "_evofile";
    if ((a = rp.find_parameter(temp.str())) != "") {
      evofile = a;
    }
    else {
      evofile = "NOFILE";
    }

    temp.str("");
    temp << "WIND_" << i << "_enhance_mdot";
    if ((a = rp.find_parameter(temp.str())) != "") {
      enhance_mdot = atof(a.c_str());
    }
    else {
      enhance_mdot = 0;
    }

    temp.str("");
    temp << "WIND_" << i << "_t_offset";
    if ((a = rp.find_parameter(temp.str())) != "") {
      time_offset = atof(a.c_str()) * pconst.year();
    }
    else {
      time_offset = -1.0e99;
    }

    temp.str("");
    temp << "WIND_" << i << "_updatefreq";
    if ((a = rp.find_parameter(temp.str())) != "") {
      update_freq = atof(a.c_str()) * pconst.year();
    }
    else {
      update_freq = -1.0e99;
    }

    temp.str("");
    temp << "WIND_" << i << "_t_scalefac";
    if ((a = rp.find_parameter(temp.str())) != "") {
      time_scalefac = atof(a.c_str());
    }
    else {
      time_scalefac = 1.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_eccentricity";
    if ((a = rp.find_parameter(temp.str())) != "") {
      eccentricity = atof(a.c_str());
    }
    else {
      eccentricity = 0.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_orbital_period";
    if ((a = rp.find_parameter(temp.str())) != "") {
      OrbPeriod = atof(a.c_str()) * pconst.year();  // convert to sec
    }
    else {
      OrbPeriod = 0.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_periastron_vec_x";
    if ((a = rp.find_parameter(temp.str())) != "")
      PeriastronX = atof(a.c_str());
    else {
      spdlog::debug(
          "No periastron_vec_x in pfile for wind {}, setting to 0.0", i);
      PeriastronX = 0.0;
    }

    temp.str("");
    temp << "WIND_" << i << "_periastron_vec_y";
    if ((a = rp.find_parameter(temp.str())) != "")
      PeriastronY = atof(a.c_str());
    else {
      spdlog::debug(
          "No periastron_vec_y in pfile for wind {}, setting to 0.0", i);
      PeriastronY = 0.0;
    }

    temp.str("");
    temp << "WIND_" << i << "_xi";
    if ((a = rp.find_parameter(temp.str())) != "") {
      xi = atof(a.c_str());
    }
    else {
      xi = -0.43;
    }  // default value from Bjorkman & Cassinelli (1993)

    int moving_star = 0;
    temp.str("");
    temp << "WIND_" << i << "_moving_star";
    if ((a = rp.find_parameter(temp.str())) != "") {
      moving_star = atoi(a.c_str());
    }

    // stellar mass in solar masses
    double mass = 0;
    temp.str("");
    temp << "WIND_" << i << "_mass";
    if ((a = rp.find_parameter(temp.str())) != "") {
      mass = atof(a.c_str()) * pconst.Msun();
    }

    //
    // Now we should have got all the sources, so add the source to
    // the global list.
    //
    struct stellarwind_params *wind = 0;
    wind                            = mem.myalloc(wind, 1);
    wind->id                        = i;
    for (int v = 0; v < MAX_DIM; v++)
      wind->dpos[v] = posn[v];
    for (int v = 0; v < MAX_DIM; v++)
      wind->velocity[v] = vel[v];
    wind->radius      = rad;
    wind->Mass        = mass;
    wind->Mdot        = Mdot;
    wind->Vinf        = Vinf;
    wind->Vrot        = Vrot;
    wind->Tstar       = Tw;
    wind->Rstar       = Rstar;
    wind->Bstar       = Bstar;
    wind->type        = type;
    wind->moving_star = moving_star;
    for (int v = 0; v < MAX_NVAR; v++) {
      wind->tr[v] = trcr[v];
    }
    //
    // new stuff for evolving winds:
    //
    wind->evolving_wind_file = evofile;
    wind->enhance_mdot       = enhance_mdot;
    wind->xi                 = xi;
    wind->time_offset        = time_offset;
    wind->update_freq        = update_freq;
    wind->t_scalefactor      = time_scalefac;
    wind->eccentricity       = eccentricity;
    wind->OrbPeriod          = OrbPeriod;
    wind->PeriastronX        = PeriastronX;
    wind->PeriastronY        = PeriastronY;
    // cout <<"\tgot parameters, adding source! rad="<<rad<<"\n";
    SWP.params.push_back(wind);
    // cout <<"\tadded WIND source. returning.\n";
  }

  return err;
}

// ##################################################################
// ##################################################################

int SimParams::read_extra_physics()
{
  // read extra physics parameters if present, and set them to
  // zero if not.
  string a;
  if ((a = rp.find_parameter("EP_dynamics")) != "")
    EP.dynamics = atoi(a.c_str());
  else
    EP.dynamics = 1;

  if ((a = rp.find_parameter("EP_raytracing")) != "")
    EP.raytracing = atoi(a.c_str());
  else
    EP.raytracing = 0;

  if ((a = rp.find_parameter("EP_cooling")) != "")
    EP.cooling = atoi(a.c_str());
  else
    EP.cooling = 0;

  if ((a = rp.find_parameter("EP_chemistry")) != "")
    EP.chemistry = atoi(a.c_str());
  else
    EP.chemistry = 0;

  if ((a = rp.find_parameter("EP_coll_ionisation")) != "")
    EP.coll_ionisation = atoi(a.c_str());
  else
    EP.coll_ionisation = 0;

  if ((a = rp.find_parameter("EP_phot_ionisation")) != "")
    EP.phot_ionisation = atoi(a.c_str());
  else
    EP.phot_ionisation = 0;

  if ((a = rp.find_parameter("EP_rad_recombination")) != "")
    EP.rad_recombination = atoi(a.c_str());
  else
    EP.rad_recombination = 0;

  if ((a = rp.find_parameter("EP_update_erg")) != "")
    EP.update_erg = atoi(a.c_str());
  else
    EP.update_erg = 1;

  if ((a = rp.find_parameter("EP_MP_timestep_limit")) != "")
    EP.MP_timestep_limit = atoi(a.c_str());
  else
    EP.MP_timestep_limit = 1;

  if ((a = rp.find_parameter("EP_Min_Temperature")) != "")
    EP.MinTemperature = atof(a.c_str());
  else
    EP.MinTemperature = 0.0;

  if ((a = rp.find_parameter("EP_Max_Temperature")) != "")
    EP.MaxTemperature = atof(a.c_str());
  else
    EP.MaxTemperature = 1.0e100;

  //
  // Hydrogen abundance (by mass) X.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp.find_parameter("EP_Hydrogen_MassFrac")) != "")
    EP.H_MassFrac = atof(a.c_str());
  else
    EP.H_MassFrac = 0.7154;

  //
  // Helium abundance (by mass) Y.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp.find_parameter("EP_Helium_MassFrac")) != "")
    EP.Helium_MassFrac = atof(a.c_str());
  else
    EP.Helium_MassFrac = 0.2703;

  //
  // Metal abundance (by mass) Z.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp.find_parameter("EP_Metal_MassFrac")) != "")
    EP.Metal_MassFrac = atof(a.c_str());
  else
    EP.Metal_MassFrac = 0.0142;

  // default to not use thermal conduction:
  if ((a = rp.find_parameter("EP_sat_thermal_cond")) != "")
    EP.sat_thermal_cond = atoi(a.c_str());
  else
    EP.sat_thermal_cond = 0;

  // default to use the Cowie & McKee (1977) saturated flux with phi=0.3
  if ((a = rp.find_parameter("EP_tc_strength")) != "")
    EP.tc_strength = atof(a.c_str());
  else
    EP.tc_strength = 1.0;

  return 0;
}

// ##################################################################
// ##################################################################

int SimParams::read_units()
{
  string unit;
  if ((unit = rp.find_parameter("units")) == "") {
    spdlog::debug("No units specified in parameterfile.  Assuming cgs");
    uc.unitsys  = "cgs";
    uc.density  = "g/cm^3";
    uc.length   = "cm";
    uc.velocity = "cm/s";
    uc.bfield   = "gauss";
    uc.rhoVal = uc.lenVal = uc.velVal = uc.magVal = 1.;
  }
  else if (unit == "MKS" || unit == "SI") {
    spdlog::debug("\t Using MKS (SI) as reference units");
    uc.unitsys  = "mks";
    uc.density  = "kg/m^3";
    uc.length   = "m";
    uc.velocity = "m/s";
    uc.bfield   = "Tesla/4pi";
    uc.rhoVal   = atof((rp.find_parameter("rhoval")).c_str());
    uc.lenVal   = atof((rp.find_parameter("lenval")).c_str());
    uc.velVal   = atof((rp.find_parameter("velval")).c_str());
    uc.magVal   = atof((rp.find_parameter("magval")).c_str());
  }
  else if (unit == "cgs" || unit == "CGS") {
    spdlog::debug("\t Using cgs as reference units");
    uc.unitsys  = "cgs";
    uc.density  = "g/cm^3";
    uc.length   = "cm";
    uc.velocity = "cm/s";
    uc.bfield   = "gauss";
    uc.rhoVal   = atof((rp.find_parameter("rhoval")).c_str());
    uc.lenVal   = atof((rp.find_parameter("lenval")).c_str());
    uc.velVal   = atof((rp.find_parameter("velval")).c_str());
    uc.magVal   = atof((rp.find_parameter("magval")).c_str());
  }
  else {
    spdlog::error("{}: {}", "Don't recognise units system", unit);
  }
  return 0;
}

// ##################################################################
// ##################################################################



int SimParams::read_jet_params(class JetParams &jpar  ///< jet parameters class.
)
{
  string a;
  if ((a = rp.find_parameter("JETradius")) != "")
    jpar.jetradius = atoi(a.c_str());
  else
    spdlog::error("{}: {}", "failed to find par: JETradius", 0);

  if ((a = rp.find_parameter("JETdensity")) != "")
    jpar.jetstate[RO] = atof(a.c_str());
  else
    spdlog::error("{}: {}", "failed to find par: JETdensity", 0);

  if ((a = rp.find_parameter("JETpressure")) != "")
    jpar.jetstate[PG] = atof(a.c_str());
  else
    spdlog::error("{}: {}", "failed to find par: JETpressure", 0);

  if ((a = rp.find_parameter("JETvelocity")) != "")
    jpar.jetstate[VX] = atof(a.c_str());
  else
    spdlog::error("{}: {}", "failed to find par: JETvelocity", 0);
  jpar.jetstate[VY] = 0.0;
  jpar.jetstate[VZ] = 0.0;

  if ((a = rp.find_parameter("JET_Bax")) != "")
    jpar.jetstate[BX] = atof(a.c_str());
  if ((a = rp.find_parameter("JET_Btor")) != "")
    jpar.jetstate[BY] = atof(a.c_str());
  jpar.jetstate[BZ] = 0.0;
#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  jpar.jetstate[BX] /= sqrt(4.0 * M_PI);
  jpar.jetstate[BY] /= sqrt(4.0 * M_PI);
  jpar.jetstate[BZ] /= sqrt(4.0 * M_PI);
#endif

  if (eqntype == EQGLM) {
    jpar.jetstate[SI] = 0.0;
  }

  ostringstream temp;
  for (int v = 0; v < ntracer; v++) {
    temp.str("");
    temp << "JETjetTR" << v;
    if ((a = rp.find_parameter(temp.str())) != "")
      jpar.jetstate[ftr + v] = atof(a.c_str());
    else
      spdlog::error("{}: {}", "param not found in pfile", temp.str());
  }
  for (int v = ftr + ntracer; v < MAX_NVAR; v++)
    jpar.jetstate[v] = 0.0;

  return 0;
}



// ##################################################################
// ##################################################################


std::vector<int> SimParams::get_pbc_bools() const
{
  std::vector<int> pbc(2 * ndim, 0);
  if (BC_XN == "periodic") {
    pbc[0] = 1;
  }
  if (BC_XP == "periodic") {
    pbc[1] = 1;
  }
  if (ndim > 1) {
    if (BC_YN == "periodic") {
      pbc[2] = 1;
    }
    if (BC_YP == "periodic") {
      pbc[3] = 1;
    }
  }
  if (ndim > 2) {
    if (BC_ZN == "periodic") {
      pbc[4] = 1;
    }
    if (BC_ZP == "periodic") {
      pbc[5] = 1;
    }
  }
  return pbc;
}

// ##################################################################
// ##################################################################
