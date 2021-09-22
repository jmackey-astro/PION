/// \file get_sim_info.cc
///
/// \author Jonathan Mackey
///
/// Gets simulation info from an ASCII text parameter file.
///
/// modifications:
///
/// - 2010.10.01 JM: Added spherical coordinate option.  Fixed bug in
///   wind_src routines where the read function was called when
///   Nsrc==0.  Added comments at start of file!
/// - 2010.12.29 JM: Added Internal--energy--solving equations.
/// - 2010.01.06 JM: New stellar wind interface. (bug fix 01.18 JM).
/// - 2011.02.15 JM: Added new stellar wind parameters for evolving
///    wind sources.
/// - 2011.02.28 JM: Got rid of RSP radiation-source-parameters class
///    references.
/// - 2011.03.02 JM: Got rid of references to MAX_NTR (no longer used).
/// - 2011.03.21 JM: Added RT_s_update_N as another radiation source
///    parameter. (and 22.03, changed notation of RT flags 15.04)
/// - 2011.05.02 JM: new RT params.
/// - 2011.06.02 JM: some debugging text additions. read_radsources()
///    always called now.
/// - 2012.01.14 JM: Added RT_EVO_FILE_[i] (optional) for time-
///    varying radiation source.
/// - 2013.02.14 JM: Added He/Metal mass fractions as EP parameters,
///    to make metallicity and mu into parameterfile settings.
/// - 2013.04.15 JM: Removed lots of cout/cerr statements to clean up
///    std i/o messages.
/// - 2013.08.19 JM: Added Hydrogen MassFrac to EP parameter list
/// - 2013.08.20 JM: Modified cell_interface for optical depth vars.
/// - 2013.08.23 JM: Added new mpv9_HHe module code.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.06/07 JM: Change tracer setup in files, so that each
///    tracer has its own variable.
/// - 2017.03.07 JM: changed logic so that ArtificalViscosity=4 is
///    allowed.
/// - 2017.11.07-22 JM: updating boundary setup.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "ics/get_sim_info.h"

#include <sstream>
using namespace std;

get_sim_info::get_sim_info()
{
  rp = 0;
  return;
}

get_sim_info::~get_sim_info()
{
  if (rp) {
    delete rp;
    rp = 0;
  }
}

int get_sim_info::read_gridparams(
    string pfile,           ///< paramfile.
    class SimParams &SimPM  ///< pointer to simulation paramters.
)
{
  int err          = 0;
  get_sim_info::rp = 0;
  rp               = new ReadParams;
  if (!rp) rep.error("get_sim_info::read_gridparams() initialising RP", rp);

  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);
  //   rp->write_out_parameters();

  SimPM.gridType = 1;  // uniform grid. only option for now.

  // Sytem of Equations to solve:
  string str = rp->find_parameter("eqn");
  if (str == "hd" || str == "euler") {
    // inviscid euler equations.
    SimPM.eqntype = EQEUL;
    SimPM.nvar    = 5;
  }
  else if (str == "hd-Eint" || str == "EulerEint") {
    // inviscid Euler equations, integrating internal energy.
    SimPM.eqntype = EQEUL_EINT;
    SimPM.nvar    = 5;
  }
  else if (str == "iso-hd" || str == "iso-eul") {
    //
    // inviscid isothermal hydro equations.
    //
    SimPM.eqntype = EQEUL_ISO;
    SimPM.nvar    = 5;
  }
  else if (str == "i-mhd") {
    // ideal mhd with no divB cleaning (not great!)
    SimPM.eqntype = EQMHD;
    SimPM.nvar    = 8;
  }
  else if (str == "mhd" || str == "glm-mhd") {
    // dedner's glm method (GOOD!)
    SimPM.eqntype = EQGLM;
    SimPM.nvar    = 9;
  }
  else if (str == "fieldCD" || str == "fcd-mhd") {
    // toth's field CD method (BAD! and not working!)
    SimPM.eqntype = EQFCD;
    SimPM.nvar    = 8;
  }
  else
    rep.error("Don't know what these equations are:", str);

  str = rp->find_parameter("solver");
  if (str == "")
    SimPM.solverType = 4;  // Roe-type solver (as opposed to lax-friedrichs).
  else if (str == "LF" || str == "LAX-FRIEDRICHS" || str == "0")
    SimPM.solverType = 0;
  else if (str == "LINEAR" || str == "linear" || str == "1")
    SimPM.solverType = 1;
  else if (str == "EXACT" || str == "exact" || str == "2")
    SimPM.solverType = 2;
  else if (str == "HYBRID" || str == "hybrid" || str == "3")
    SimPM.solverType = 3;
  else if (str == "ROE" || str == "roe" || str == "4")
    SimPM.solverType = 4;
  else if (str == "ROEPV" || str == "roepv" || str == "5")
    SimPM.solverType = 5;
  else if (str == "FVS" || str == "fvs" || str == "6")
    SimPM.solverType = 6;
  else if (str == "HLLD" || str == "hlld" || str == "7")
    SimPM.solverType = 7;
  else if (str == "HLL" || str == "hll" || str == "8")
    SimPM.solverType = 8;
  else
    rep.error("what kind of solver is this???", str);

  // dimensionality of grid.
  str = rp->find_parameter("ndim");
  if (str == "")
    rep.error("Bad ndim in pfile", str);
  else
    SimPM.ndim = atoi(str.c_str());
  if (isnan(SimPM.ndim)) rep.error("ndim is not a number", SimPM.ndim);
  const int ndim = SimPM.ndim;
  SimPM.eqnNDim  = 3;

  //
  // tracer variables: number of tracers, and then each one gets a
  // string name, stored in SimPM.tracers
  //
  str = rp->find_parameter("ntracer");
  if (str == "") {
#ifndef NDEBUG
    cout << "Not using tracer variables.\n";
#endif
    SimPM.ntracer = 0;
    SimPM.ftr     = SimPM.nvar;
  }
  else if (isfinite(SimPM.ntracer = atoi(str.c_str()))) {
#ifndef NDEBUG
    cout << "using " << SimPM.ntracer << " passive tracer variables\n";
#endif
    SimPM.ftr = SimPM.nvar;
    SimPM.nvar += SimPM.ntracer;
    //
    // Get what type of chemistry we are doing:
    //
    SimPM.chem_code = rp->find_parameter("chem_code");

    //
    // Each tracer has its own parameter, called Tracer000, Tracer001,
    // etc., so setup a string for each in the tracers array.
    //
    SimPM.tracers = mem.myalloc(SimPM.tracers, SimPM.ntracer);
    for (int i = 0; i < SimPM.ntracer; i++) {
      ostringstream temp;
      temp << "Tracer";
      temp.width(3);
      temp.fill('0');
      temp << i;
      SimPM.tracers[i] = rp->find_parameter(temp.str());
#ifndef NDEBUG
      cout << "using tracer(s) described as " << SimPM.tracers[i] << "\n";
#endif
      if (SimPM.tracers[i] == "") {
        rep.error("Can't find tracer name for number", i);
      }
    }
  }
  else
    rep.error("number of tracers is not finite!", SimPM.ntracer);

  // coordinate system.
  str = rp->find_parameter("coordinates");
  if (str == "") {
    cout << "no coordinate system specified.  using cartesian as default.\n";
    SimPM.coord_sys = COORD_CRT;
  }
  else if (str == "cartesian" || str == "cart")
    SimPM.coord_sys = COORD_CRT;
  else if (str == "cylindrical" || str == "cyl")
    SimPM.coord_sys = COORD_CYL;
  else if (str == "spherical" || str == "sph")
    SimPM.coord_sys = COORD_SPH;
  else
    rep.error("Don't recognise coordinate system in GetParameters", str);

  // Domain of grid, and number of points.  These must all be present in the
  // pfile.
  string seek;
  seek = "NGridX";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.NG[XX] = atoi(str.c_str());
  SimPM.Ncell  = SimPM.NG[XX];

  seek = "Xmin";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.Xmin[XX] = atof(str.c_str());

  seek = "Xmax";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.Xmax[XX]  = atof(str.c_str());
  SimPM.Range[XX] = SimPM.Xmax[XX] - SimPM.Xmin[XX];

  if (ndim > 1) {
    seek = "NGridY";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.NG[YY] = atoi(str.c_str());
    SimPM.Ncell *= SimPM.NG[YY];

    seek = "Ymin";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.Xmin[YY] = atof(str.c_str());

    seek = "Ymax";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.Xmax[YY]  = atof(str.c_str());
    SimPM.Range[YY] = SimPM.Xmax[YY] - SimPM.Xmin[YY];
  }
  else {
    SimPM.NG[YY]   = 1;
    SimPM.Xmin[YY] = SimPM.Xmax[YY] = SimPM.Range[YY] = 0.;
  }

  if (ndim > 2) {
    seek = "NGridZ";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.NG[ZZ] = atoi(str.c_str());
    SimPM.Ncell *= SimPM.NG[ZZ];

    seek = "Zmin";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.Xmin[ZZ] = atof(str.c_str());

    seek = "Zmax";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.Xmax[ZZ]  = atof(str.c_str());
    SimPM.Range[ZZ] = SimPM.Xmax[ZZ] - SimPM.Xmin[ZZ];
  }
  else {
    SimPM.NG[ZZ]   = 1;
    SimPM.Xmin[ZZ] = SimPM.Xmax[ZZ] = SimPM.Range[ZZ] = 0.;
  }

  SimPM.dx = SimPM.Range[XX] / SimPM.NG[XX];

  if (ndim > 1) {
    if (fabs(SimPM.Range[1] / (SimPM.NG[1]) / SimPM.dx - 1.)
        > 100. * MACHINEACCURACY) {
      rep.error(
          "Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          fabs(SimPM.Range[1] / (SimPM.NG[1]) / SimPM.dx - 1.));
    }
  }
  if (ndim > 2) {
    if (fabs(SimPM.Range[2] / (SimPM.NG[2]) / SimPM.dx - 1.)
        > 100. * MACHINEACCURACY) {
      rep.error(
          "Cells must be same length in each direction! Set the range "
          "and number of points appropriately.",
          fabs(SimPM.Range[2] / (SimPM.NG[2]) / SimPM.dx - 1.));
    }
  }

  // Nested grid info (all optional for now)
  seek = "grid_nlevels";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.grid_nlevels = 1;
  else
    SimPM.grid_nlevels = atoi(str.c_str());

  seek = "grid_aspect_ratio_XX";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.grid_aspect_ratio[XX] = 1;
  else
    SimPM.grid_aspect_ratio[XX] = atoi(str.c_str());

  seek = "grid_aspect_ratio_YY";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.grid_aspect_ratio[YY] = 1;
  else
    SimPM.grid_aspect_ratio[YY] = atoi(str.c_str());

  seek = "grid_aspect_ratio_ZZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.grid_aspect_ratio[ZZ] = 1;
  else
    SimPM.grid_aspect_ratio[ZZ] = atoi(str.c_str());

  seek = "NG_centre_XX";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_centre[XX] = 0.0;
  else
    SimPM.NG_centre[XX] = atof(str.c_str());

  seek = "NG_centre_YY";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_centre[YY] = 0.0;
  else
    SimPM.NG_centre[YY] = atof(str.c_str());

  seek = "NG_centre_ZZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_centre[ZZ] = 0.0;
  else
    SimPM.NG_centre[ZZ] = atof(str.c_str());

  seek = "NG_refine_XX";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_refine[XX] = 1;
  else
    SimPM.NG_refine[XX] = atoi(str.c_str());

  seek = "NG_refine_YY";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_refine[YY] = 1;
  else
    SimPM.NG_refine[YY] = atoi(str.c_str());

  seek = "NG_refine_ZZ";
  str  = rp->find_parameter(seek);
  if (str == "")
    SimPM.NG_refine[ZZ] = 1;
  else
    SimPM.NG_refine[ZZ] = atoi(str.c_str());

  // output info
  str = rp->find_parameter("OutputPath");
  if (str == "") rep.error("outputpath", str);
  ostringstream temp;
  temp << str;
  str = rp->find_parameter("OutputFile");
  if (str == "") rep.error("outputfile", str);
  temp << str;
  SimPM.outFileBase = temp.str();
  str               = rp->find_parameter("OutputFileType");

  if (str == "text")
    SimPM.typeofop = 1;
  else if (str == "fits")
    SimPM.typeofop = 2;
  else if (str == "ftab") {  // fits table (not implemented yet.
    SimPM.typeofop = 3;
    rep.error("Fits Table not allowed as o/p format.. not implemented", str);
  }
  else if (str == "both")
    SimPM.typeofop = 4;  // Both means fits and text.
  else if (str == "silo")
    SimPM.typeofop = 5;  // Silo output.
  else
    rep.error("Error, bad type of outputfile in pfile.", str);

  seek = "OutputFrequency";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.opfreq = atoi(str.c_str());
  // cout <<"\tOutFile: "<<SimPM.outFileBase<< ".xxx\t Type="<<SimPM.typeofop;
  // cout <<" every "<<SimPM.opfreq<<" timesteps."<<"\n";

  // Optional output-by-years parameters:
  seek = "OutputCriterion";
  str  = rp->find_parameter(seek);
  if (str == "") {
    SimPM.op_criterion = 0;
  }
  else
    SimPM.op_criterion = atoi(str.c_str());
  // cout <<"\top_criterion = "<<SimPM.op_criterion<<"\n";

  seek = "OPfreqTime";
  str  = rp->find_parameter(seek);
  if (str == "") {
    SimPM.opfreq_time = 0.0;
    SimPM.next_optime = 0.0;
  }
  else {
    SimPM.opfreq_time = atof(str.c_str());
    SimPM.next_optime = SimPM.simtime + SimPM.opfreq_time;
    double tmp        = ((SimPM.simtime / SimPM.opfreq_time)
                  - static_cast<int>(SimPM.simtime / SimPM.opfreq_time))
                 * SimPM.opfreq_time;
    SimPM.next_optime -= tmp;
  }
  // cout <<" and opfreq="<<SimPM.opfreq_time<<" code time units.\n";

  //
  // Boundary conditions.
  // First do the edges of the domain:
  //
  SimPM.BC_XN = rp->find_parameter("BC_XN");
  SimPM.BC_XP = rp->find_parameter("BC_XP");
  if (SimPM.ndim > 1) {
    SimPM.BC_YN = rp->find_parameter("BC_YN");
    SimPM.BC_YP = rp->find_parameter("BC_YP");
  }
  else {
    SimPM.BC_YN = "NONE";
    SimPM.BC_YP = "NONE";
  }
  if (SimPM.ndim > 2) {
    SimPM.BC_ZN = rp->find_parameter("BC_ZN");
    SimPM.BC_ZP = rp->find_parameter("BC_ZP");
  }
  else {
    SimPM.BC_ZN = "NONE";
    SimPM.BC_ZP = "NONE";
  }
  //
  // Now do the internal boundaries (if any).  Seek the string
  // for a given internal boundary, and add to a vector until
  // no more are found.
  //
  str = rp->find_parameter("BC_Ninternal");
  if (str == "")
    SimPM.BC_Nint = 0;
  else
    SimPM.BC_Nint = atoi(str.c_str());
  if (SimPM.BC_Nint > 0) {
    SimPM.BC_INT = mem.myalloc(SimPM.BC_INT, SimPM.BC_Nint);
  }
  for (int v = 0; v < SimPM.BC_Nint; v++) {
    ostringstream intbc;
    intbc.str("");
    intbc << "BC_INTERNAL_";
    intbc.width(3);
    intbc.fill('0');
    intbc << v;
    // cout <<"Looking for internal boundary: "<<intbc.str();
    str = rp->find_parameter(intbc.str());
    // cout <<"   Found: "<<str<<"\n";
    SimPM.BC_INT[v] = str;
  }
  // rep.printVec("BC_INT",SimPM.BC_INT,SimPM.BC_Nint);

  SimPM.Nbc = -1;  // Set it to negative so I know it's not set.

  // Timing
  seek = "StartTime";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.starttime = atof(str.c_str());
  SimPM.simtime   = SimPM.starttime;
  seek            = "FinishTime";
  str             = rp->find_parameter(seek);
  if (str == "")
    SimPM.finishtime = -1.;
  else
    SimPM.finishtime = atof(str.c_str());
  SimPM.timestep = 0;  // Counter for what timestep we are on.

  seek = "OrderOfAccSpace";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.spOOA = atoi(str.c_str());
  seek        = "OrderOfAccTime";
  str         = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.tmOOA = atoi(str.c_str());

  // Physics
  seek = "GAMMA";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.gamma = atof(str.c_str());
  seek        = "CFL";
  str         = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  SimPM.CFL = atof(str.c_str());

  seek = "ArtificialViscosity";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("param not found", seek);
  if ((SimPM.artviscosity = atoi(str.c_str())) == 0) {
    SimPM.etav = 0.;
  }
  else if (SimPM.artviscosity == 1 || SimPM.artviscosity == 4) {
    seek = "EtaViscosity";
    str  = rp->find_parameter(seek);
    if (str == "") rep.error("param not found", seek);
    SimPM.etav = atof(str.c_str());
  }
  else if (SimPM.artviscosity == 3) {
    // using H-correction.
    SimPM.etav = 0.1;
  }
  else
    rep.error("\tUnknown viscosity requested... fix me.", str);
  // cout <<"\tArtificial Viscosity: eta="<<SimPM.etav<<"\n";
  // Which Physics
  err += read_extra_physics(SimPM);
  if (err) rep.error("read_extra_physics", err);

  // Raytracing
  err += read_radsources(SimPM);
  if (err) rep.error("read_radsources", err);

  //
  // STELLAR WINDS???
  //
  seek = "WIND_NSRC";
  str  = rp->find_parameter(seek);
  if (str == "") {
    cout << "\tWIND_NSRC: param not found, setting to 0\n";
    SWP.Nsources = 0;
  }
  else {
    SWP.Nsources = atoi(str.c_str());
    // cout <<"\tWIND_NSRC: got "<<SWP.Nsources<<" sources.\n";
    if (SWP.Nsources > 0) {
      err += read_wind_sources(SimPM);
      if (err) rep.error("read_wind_sources", err);
    }
  }

  //
  // Jets?
  //
  seek = "N_JET";
  str  = rp->find_parameter(seek);
  if (str == "") {
    cout << "\tN_JET: param not found, setting to 0\n";
    JP.jetic = 0;
  }
  else {
    JP.jetic = atoi(str.c_str());
    if (JP.jetic) {
      err += read_jet_params(SimPM, JP);
      rep.errorTest("read_jet_params", 0, err);
    }
  }

  //
  // Reference State Vector.  Only look for the first nvar elements and set
  // the rest to zero.
  //
  if ((str = rp->find_parameter("refvec0")) != "") {
    for (int v = 0; v < SimPM.nvar; v++) {
      ostringstream temp2;
      temp2 << "refvec" << v;
      str = rp->find_parameter(temp2.str());
      if (str != "") {
        // cout <<"v="<<v<<" MAX_NVAR="<<MAX_NVAR<<"\n";
        SimPM.RefVec[v] = atof(str.c_str());
      }
      else {
        cout << "no refvec[" << v << "], setting to 1.\n";
        SimPM.RefVec[v] = 1.;
      }
    }
    for (int v = SimPM.nvar; v < MAX_NVAR; v++)
      SimPM.RefVec[v] = 1.0;
  }
  else {
    for (int v = 0; v < MAX_NVAR; v++)
      SimPM.RefVec[v] = 1.;
  }

  // Initialise gridparams so that we assume we are not doing a jet sim.
  JP.jetic = 0;

  // Units Used
  err += read_units();
  if (err) rep.error("read_units", err);

  // seek=""; str=rp->find_parameter(seek); if (str=="") rep.error("param not
  // found",seek);
  return err;
}

// ##################################################################
// ##################################################################

int get_sim_info::read_radsources(
    class SimParams &SimPM  ///< pointer to simulation paramters.
)
{
  if (!rp) return 1;
  string a;
  if ((a = rp->find_parameter("RT_Nsources")) != "")
    SimPM.RS.Nsources = atoi(a.c_str());
  else {
    // cout <<"no nsources in pfile... assuming there are no RT sources.\n";
    SimPM.RS.Nsources = 0;
  }

  for (int i = 0; i < SimPM.RS.Nsources; i++) {
    struct rad_src_info rs_temp;
    ostringstream temp2;

    for (int v = 0; v < SimPM.ndim; v++) {
      temp2.str("");
      temp2 << "RT_position_" << i << "_" << v;
      if ((a = rp->find_parameter(temp2.str())) != "")
        rs_temp.pos[v] = atof(a.c_str());
      else
        rep.error("no src position in pfile", temp2.str());
    }
    for (int v = SimPM.ndim; v < MAX_DIM; v++)
      rs_temp.pos[v] = 0.0;

    temp2.str("");
    temp2 << "RT_strength_" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.strength = atof(a.c_str());
    else
      rep.error("no src strength in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_Rstar____" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.Rstar = atof(a.c_str());
    else {
      cout << temp2.str()
           << ": parameter not found in pfile.  Setting to -1.0e200.\n";
      rs_temp.Rstar = -1.0e200;
    }

    temp2.str("");
    temp2 << "RT_Tstar____" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.Tstar = atof(a.c_str());
    else {
      cout << temp2.str()
           << ": parameter not found in pfile.  Setting to -1.0e200.\n";
      rs_temp.Tstar = -1.0e200;
    }

    rs_temp.id = i;

    temp2.str("");
    temp2 << "RT_src_type_" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.type = atoi(a.c_str());
    else
      rep.error("no src type in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_update___" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.update = atoi(a.c_str());
    else {
      // rep.error("no src update-method in pfile",temp2.str());
      cerr << "*** no src update-method in pfile: " << temp2.str();
      cerr << "... using C2Ray method. \n";
      rs_temp.update = RT_UPDATE_IMPLICIT;
    }

    temp2.str("");
    temp2 << "RT_at_infty_" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.at_infinity = atoi(a.c_str());
    else
      rep.error("no src at-infinity in pfile", temp2.str());

    temp2.str("");
    temp2 << "RT_effect___" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.effect = atoi(a.c_str());
    else {
      // rep.error("no src update-method in pfile",temp2.str());
      cerr << "*** no src effect param in pfile: " << temp2.str();
      cerr << "... using monochromatic photoionisation.\n";
      rs_temp.effect = RT_EFFECT_PION_MONO;
    }

    temp2.str("");
    temp2 << "RT_Tau_src__" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.opacity_src = atoi(a.c_str());
    else {
      // rep.error("no src update-method in pfile",temp2.str());
      cerr << "*** no opacity src param in pfile: " << temp2.str();
      cout << "... using neutral hydrogen.\n";
      rs_temp.opacity_src = RT_OPACITY_MINUS;
    }

    temp2.str("");
    temp2 << "RT_Tau_var__" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.opacity_var = atoi(a.c_str());
    else {
      // rep.error("no src update-method in pfile",temp2.str());
      cerr << "*** no opacity variable index in pfile: " << temp2.str();
      cerr << "... using First tracer.\n";
      rs_temp.opacity_var = SimPM.ftr;
    }

    temp2.str("");
    temp2 << "RT_EVO_FILE_" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.EvoFile = a;
    else {
      cerr << "*** no RS Evolution File in pfile: " << temp2.str();
      cerr << "... using NOFILE, i.e. constant source.\n";
      rs_temp.EvoFile = "NOFILE";
    }

    temp2.str("");
    temp2 << "RT_Nbins____" << i;
    if ((a = rp->find_parameter(temp2.str())) != "")
      rs_temp.NTau = atoi(a.c_str());
    else {
      cerr << "*** no NTau in parameter file: " << temp2.str();
      cerr << "... using 1.\n";
      rs_temp.NTau = 1;
    }

    //
    // Add source i to SimPM.RS list.
    //
    SimPM.RS.sources.push_back(rs_temp);

    cout << "\tRT-src: i=" << i << ": strength=" << rs_temp.strength;
    cout << " type=" << rs_temp.type << " at_inf=" << rs_temp.at_infinity;
    cout << ", update=" << rs_temp.update;
    cout << " and pos[]=";
    rep.printVec("", rs_temp.pos, SimPM.ndim);

    // cout <<"\tRT-src: i="<<i<<": strength=";
    // cout <<SimPM.RS.sources[i].strength;
    // cout <<" type="<<SimPM.RS.sources[i].type<<" at_inf=";
    // cout <<SimPM.RS.sources[i].at_infinity;
    // cout <<", update="<<SimPM.RS.sources[i].update;
    // cout <<" and pos[]=";
    // rep.printVec("",SimPM.RS.sources[i].pos,SimPM.ndim);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int get_sim_info::read_wind_sources(
    class SimParams &SimPM  ///< pointer to simulation paramters.
)
{
  if (!rp) return 1;
  if (SWP.Nsources == 0) return 2;
  int err = 0;

  for (int i = 0; i < SWP.Nsources; i++) {
    // cout <<"\tREADING WIND SOURCE "<<i<<"\n";
    string a;
    ostringstream temp;
    double Mdot = 0.0, posn[MAX_DIM], Vinf = 0.0, Vrot = 0.0, Tw = 0.0,
           Rstar = 0.0, Bstar = 0.0, trcr[MAX_NVAR], xi = 0.0, rad = 0.0;
    int type = 0;
    //
    // new stuff for evolving winds:
    //
    string evofile;
    int enhance_mdot = 0;
    double time_offset, update_freq, time_scalefac, ecentricity, OrbPeriod,
        PeriastronX, PeriastronY;

    for (int v = 0; v < SimPM.ndim; v++) {
      temp.str("");
      temp << "WIND_" << i << "_pos" << v;
      if ((a = rp->find_parameter(temp.str())) != "")
        posn[v] = atof(a.c_str());
      else
        rep.error("no src position in pfile", temp.str());
    }
    for (int v = SimPM.ndim; v < MAX_DIM; v++)
      posn[v] = 0.0;

    temp.str("");
    temp << "WIND_" << i << "_radius";
    if ((a = rp->find_parameter(temp.str())) != "")
      rad = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_type";
    if ((a = rp->find_parameter(temp.str())) != "")
      type = atoi(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());
    // cout <<"\t*******wind type="<<type<<"\n";

    temp.str("");
    temp << "WIND_" << i << "_mdot";
    if ((a = rp->find_parameter(temp.str())) != "")
      Mdot = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_vinf";
    if ((a = rp->find_parameter(temp.str())) != "")
      Vinf = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_vrot";
    if ((a = rp->find_parameter(temp.str())) != "")
      Vrot = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_temp";
    if ((a = rp->find_parameter(temp.str())) != "")
      Tw = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_Rstr";
    if ((a = rp->find_parameter(temp.str())) != "")
      Rstar = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    temp.str("");
    temp << "WIND_" << i << "_Bsrf";
    if ((a = rp->find_parameter(temp.str())) != "")
      Bstar = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());

    for (int v = 0; v < SimPM.ntracer; v++) {
      temp.str("");
      temp << "WIND_" << i << "_TR" << v;
      if ((a = rp->find_parameter(temp.str())) != "")
        trcr[v] = atof(a.c_str());
      else
        rep.error("param not found in pfile", temp.str());
    }
    for (int v = SimPM.ntracer; v < MAX_NVAR; v++)
      trcr[v] = 0.0;

    //
    // new stuff for evolving winds:
    //
    temp.str("");
    temp << "WIND_" << i << "_evofile";
    if ((a = rp->find_parameter(temp.str())) != "") {
      evofile = a;
    }
    else {
      evofile = "NOFILE";
    }

    temp.str("");
    temp << "WIND_" << i << "_enhance_mdot";
    if ((a = rp->find_parameter(temp.str())) != "") {
      enhance_mdot = atof(a.c_str());
    }
    else {
      enhance_mdot = 0;
    }

    temp.str("");
    temp << "WIND_" << i << "_t_offset";
    if ((a = rp->find_parameter(temp.str())) != "") {
      time_offset = atof(a.c_str());
    }
    else {
      time_offset = -1.0e99;
    }

    temp.str("");
    temp << "WIND_" << i << "_updatefreq";
    if ((a = rp->find_parameter(temp.str())) != "") {
      update_freq = atof(a.c_str());
    }
    else {
      update_freq = -1.0e99;
    }

    temp.str("");
    temp << "WIND_" << i << "_t_scalefac";
    if ((a = rp->find_parameter(temp.str())) != "") {
      time_scalefac = atof(a.c_str());
    }
    else {
      time_scalefac = 1.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_ecentricity_fac";
    if ((a = rp->find_parameter(temp.str())) != "") {
      ecentricity = atof(a.c_str());
    }
    else {
      ecentricity = 0.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_orbital_period";
    if ((a = rp->find_parameter(temp.str())) != "") {
      OrbPeriod = atof(a.c_str());
    }
    else {
      OrbPeriod = 0.0;  // default value
    }

    temp.str("");
    temp << "WIND_" << i << "_periastron_vec_x";
    if ((a = rp->find_parameter(temp.str())) != "")
      PeriastronX = atof(a.c_str());
    else {
      cout << "No periastron_vec_x in pfile for wind " << i;
      cout << ", setting to 0.0\n";
      PeriastronX = 0.0;
    }

    temp.str("");
    temp << "WIND_" << i << "_periastron_vec_y";
    if ((a = rp->find_parameter(temp.str())) != "")
      PeriastronY = atof(a.c_str());
    else {
      cout << "No periastron_vec_y in pfile for wind " << i;
      cout << ", setting to 0.0\n";
      PeriastronY = 0.0;
    }


    temp.str("");
    temp << "WIND_" << i << "_xi";
    if ((a = rp->find_parameter(temp.str())) != "") {
      xi = atof(a.c_str());
    }
    else {
      xi = -0.43;
    }  // default value from Bjorkman & Cassinelli (1993)

    //
    // Now we should have got all the sources, so add the source to
    // the global list.
    //
    struct stellarwind_params *wind = 0;
    wind                            = mem.myalloc(wind, 1);
    wind->id                        = i;
    for (int v = 0; v < MAX_DIM; v++)
      wind->dpos[v] = posn[v];
    wind->radius = rad;
    wind->Mdot   = Mdot;
    wind->Vinf   = Vinf;
    wind->Vrot   = Vrot;
    wind->Tstar  = Tw;
    wind->Rstar  = Rstar;
    wind->Bstar  = Bstar;
    wind->type   = type;
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
    wind->ecentricity        = ecentricity;
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

int get_sim_info::read_extra_physics(
    class SimParams &SimPM  ///< pointer to simulation paramters.
)
{
  if (!rp) return 1;
  // read extra physics parameters if present, and set them to
  // zero if not.
  string a;
  if ((a = rp->find_parameter("EP_dynamics")) != "")
    SimPM.EP.dynamics = atoi(a.c_str());
  else
    SimPM.EP.dynamics = 1;

  if ((a = rp->find_parameter("EP_raytracing")) != "")
    SimPM.EP.raytracing = atoi(a.c_str());
  else
    SimPM.EP.raytracing = 0;

  if ((a = rp->find_parameter("EP_cooling")) != "")
    SimPM.EP.cooling = atoi(a.c_str());
  else
    SimPM.EP.cooling = 0;

  if ((a = rp->find_parameter("EP_chemistry")) != "")
    SimPM.EP.chemistry = atoi(a.c_str());
  else
    SimPM.EP.chemistry = 0;

  if ((a = rp->find_parameter("EP_coll_ionisation")) != "")
    SimPM.EP.coll_ionisation = atoi(a.c_str());
  else
    SimPM.EP.coll_ionisation = 0;

  if ((a = rp->find_parameter("EP_phot_ionisation")) != "")
    SimPM.EP.phot_ionisation = atoi(a.c_str());
  else
    SimPM.EP.phot_ionisation = 0;

  if ((a = rp->find_parameter("EP_rad_recombination")) != "")
    SimPM.EP.rad_recombination = atoi(a.c_str());
  else
    SimPM.EP.rad_recombination = 0;

  if ((a = rp->find_parameter("EP_update_erg")) != "")
    SimPM.EP.update_erg = atoi(a.c_str());
  else
    SimPM.EP.update_erg = 1;

  if ((a = rp->find_parameter("EP_MP_timestep_limit")) != "")
    SimPM.EP.MP_timestep_limit = atoi(a.c_str());
  else
    SimPM.EP.MP_timestep_limit = 1;

  if ((a = rp->find_parameter("EP_Min_Temperature")) != "")
    SimPM.EP.MinTemperature = atof(a.c_str());
  else
    SimPM.EP.MinTemperature = 0.0;

  if ((a = rp->find_parameter("EP_Max_Temperature")) != "")
    SimPM.EP.MaxTemperature = atof(a.c_str());
  else
    SimPM.EP.MaxTemperature = 1.0e100;

  //
  // Hydrogen abundance (by mass) X.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp->find_parameter("EP_Hydrogen_MassFrac")) != "")
    SimPM.EP.H_MassFrac = atof(a.c_str());
  else
    SimPM.EP.H_MassFrac = 0.7154;

  //
  // Helium abundance (by mass) Y.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp->find_parameter("EP_Helium_MassFrac")) != "")
    SimPM.EP.Helium_MassFrac = atof(a.c_str());
  else
    SimPM.EP.Helium_MassFrac = 0.2703;

  //
  // Metal abundance (by mass) Z.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  if ((a = rp->find_parameter("EP_Metal_MassFrac")) != "")
    SimPM.EP.Metal_MassFrac = atof(a.c_str());
  else
    SimPM.EP.Metal_MassFrac = 0.0142;

  return 0;
}

// ##################################################################
// ##################################################################

int get_sim_info::read_units()
{
  if (!rp) return 1;
  string unit;
  if ((unit = rp->find_parameter("units")) == "") {
    cout << "No units specified in parameterfile.  Assuming cgs.\n";
    uc.unitsys  = "cgs";
    uc.density  = "g/cm^3";
    uc.length   = "cm";
    uc.velocity = "cm/s";
    uc.bfield   = "gauss";
    uc.rhoVal = uc.lenVal = uc.velVal = uc.magVal = 1.;
  }
  else if (unit == "MKS" || unit == "SI") {
    cout << "\t Using MKS (SI) as reference units.\n";
    uc.unitsys  = "mks";
    uc.density  = "kg/m^3";
    uc.length   = "m";
    uc.velocity = "m/s";
    uc.bfield   = "Tesla/4pi";
    uc.rhoVal   = atof((rp->find_parameter("rhoval")).c_str());
    uc.lenVal   = atof((rp->find_parameter("lenval")).c_str());
    uc.velVal   = atof((rp->find_parameter("velval")).c_str());
    uc.magVal   = atof((rp->find_parameter("magval")).c_str());
  }
  else if (unit == "cgs" || unit == "CGS") {
    cout << "\t Using cgs as reference units.\n";
    uc.unitsys  = "cgs";
    uc.density  = "g/cm^3";
    uc.length   = "cm";
    uc.velocity = "cm/s";
    uc.bfield   = "gauss";
    uc.rhoVal   = atof((rp->find_parameter("rhoval")).c_str());
    uc.lenVal   = atof((rp->find_parameter("lenval")).c_str());
    uc.velVal   = atof((rp->find_parameter("velval")).c_str());
    uc.magVal   = atof((rp->find_parameter("magval")).c_str());
  }
  else {
    rep.error("Don't recognise units system", unit);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int get_sim_info::read_jet_params(
    class SimParams &s_par,  ///< simulation paramters.
    class JetParams &jpar    ///< jet parameters class.
)
{
  string a;
  if ((a = rp->find_parameter("JETradius")) != "")
    jpar.jetradius = atoi(a.c_str());
  else
    rep.error("failed to find par: JETradius", 0);

  if ((a = rp->find_parameter("JETdensity")) != "")
    jpar.jetstate[RO] = atof(a.c_str());
  else
    rep.error("failed to find par: JETdensity", 0);

  if ((a = rp->find_parameter("JETpressure")) != "")
    jpar.jetstate[PG] = atof(a.c_str());
  else
    rep.error("failed to find par: JETpressure", 0);

  if ((a = rp->find_parameter("JETvelocity")) != "")
    jpar.jetstate[VX] = atof(a.c_str());
  else
    rep.error("failed to find par: JETvelocity", 0);
  jpar.jetstate[VY] = 0.0;
  jpar.jetstate[VZ] = 0.0;

  if ((a = rp->find_parameter("JET_Bax")) != "")
    jpar.jetstate[BX] = atof(a.c_str());
  if ((a = rp->find_parameter("JET_Btor")) != "")
    jpar.jetstate[BY] = atof(a.c_str());
  jpar.jetstate[BZ] = 0.0;
#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  jpar.jetstate[BX] /= sqrt(4.0 * M_PI);
  jpar.jetstate[BY] /= sqrt(4.0 * M_PI);
  jpar.jetstate[BZ] /= sqrt(4.0 * M_PI);
#endif

  if (s_par.eqntype == EQGLM) {
    jpar.jetstate[SI] = 0.0;
  }

  ostringstream temp;
  for (int v = 0; v < s_par.ntracer; v++) {
    temp.str("");
    temp << "JETjetTR" << v;
    if ((a = rp->find_parameter(temp.str())) != "")
      jpar.jetstate[s_par.ftr + v] = atof(a.c_str());
    else
      rep.error("param not found in pfile", temp.str());
  }
  for (int v = s_par.ftr + s_par.ntracer; v < MAX_NVAR; v++)
    jpar.jetstate[v] = 0.0;

  return 0;
}

// ##################################################################
// ##################################################################
