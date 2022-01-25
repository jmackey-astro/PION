/// \file dataio_base.h
/// \author Jonathan Mackey
///
/// modified:\n
/// - 2018.05.01 JM: moved DataIOBase from dataio.h

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#ifndef NDEBUG
#endif  // NDEBUG

#include "dataio_base.h"
#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

DataIOBase::DataIOBase(
    class SimParams &SimPM  ///< pointer to simulation parameters
    ) :
    Ndigits(8)
{
  params.clear();
  jet_pm.clear();
  rt_src.clear();
  windsrc.clear();
  set_params(SimPM);
  return;
}

// ##################################################################
// ##################################################################

DataIOBase::~DataIOBase()
{
  clear_param_list(params);
  clear_param_list(jet_pm);
  clear_param_list(rt_src);
  clear_param_list(windsrc);

  return;
}

// ##################################################################
// ##################################################################

void DataIOBase::set_params(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  //
  // make sure list is empty:
  //
  if (!params.empty()) {
    spdlog::error(
        "WARNING! params list not empty. CLEARING IT, BUT THIS SHOULDN'T HAPPEN!");
    params.clear();
  }

  //
  // now create parameter structs and add them to the list.
  // CRITICAL parameters don't need a default value since the code will bug
  // out if they are missing.  Other params should get one.
  //
  // START WITH GRID PROPERTIES
  //
  pm_base *p;

  // pm_int     *p001 = new pm_int
  //  ("gridtype",    &SimPM.gridType);
  // p = p001; p->critical=true;
  // params.push_back(p);
  pm_int *p002 = new pm_int("gridndim", &SimPM.ndim);
  p            = p002;
  p->critical  = true;
  params.push_back(p);
  pm_idimarr *p003 = new pm_idimarr("NGrid", SimPM.NG.data());
  p                = p003;
  p->critical      = true;
  params.push_back(p);
  pm_long *p004 = new pm_long("Ncell", &SimPM.Ncell);
  p             = p004;
  p->critical   = true;
  params.push_back(p);
  pm_ddimarr *p005 = new pm_ddimarr("Xmin", SimPM.Xmin.data());
  p                = p005;
  p->critical      = true;
  params.push_back(p);
  pm_ddimarr *p006 = new pm_ddimarr("Xmax", SimPM.Xmax.data());
  p                = p006;
  p->critical      = true;
  params.push_back(p);

  //
  // Nested grid parameters
  //
  pm_int *n001 = new pm_int("grid_nlevels", &SimPM.grid_nlevels, 1);
  p            = n001;
  p->critical  = false;
  params.push_back(p);

  int *defaspect = new int[MAX_DIM];
  for (int v = 0; v < MAX_DIM; v++)
    defaspect[v] = 1;
  pm_idimarr *n002 =
      new pm_idimarr("grid_aspect_ratio", SimPM.grid_aspect_ratio, defaspect);
  p           = n002;
  p->critical = false;
  params.push_back(p);
  delete[] defaspect;

  double *defngcentre = new double[MAX_DIM];
  for (int v = 0; v < MAX_DIM; v++)
    defngcentre[v] = 0.0;
  pm_ddimarr *n003 = new pm_ddimarr("NG_centre", SimPM.NG_centre, defngcentre);
  p                = n003;
  p->critical      = false;
  params.push_back(p);
  delete[] defngcentre;

  int *defngrefine = new int[MAX_DIM];
  for (int v = 0; v < MAX_DIM; v++)
    defngrefine[v] = 0;
  pm_idimarr *n004 = new pm_idimarr("NG_refine", SimPM.NG_refine, defngrefine);
  p                = n004;
  p->critical      = false;
  params.push_back(p);
  delete[] defngrefine;

  //
  // Boundary conditions
  //
  // Number of internal boundaries
  pm_int *p124 = new pm_int("BC_Ninternal", &SimPM.BC_Nint);
  p            = p124;
  p->critical  = false;
  params.push_back(p);
  have_setup_bc_pm = false;  // so we know to populate list later.

  // LEGACY: look for a string parameter called BC, in case the
  // data file was written with a pre-2018 version of PION.
  pm_string *p125 = new pm_string("typeofbc_str", &SimPM.BC_STRING);
  p               = p125;
  p->critical     = false;
  params.push_back(p);

  //
  // EQUATIONS
  //
  pm_int *p008 = new pm_int("eqn_type", &SimPM.eqntype);
  p            = p008;
  p->critical  = true;
  params.push_back(p);
  pm_int *p009 = new pm_int("eqn_ndim", &SimPM.eqnNDim, 3);
  p            = p009;
  p->critical  = false;
  params.push_back(p);
  pm_int *p010 = new pm_int("eqn_nvar", &SimPM.nvar);
  p            = p010;
  p->critical  = true;
  params.push_back(p);

  //
  // Tracers
  //
  pm_int *p011 = new pm_int("num_tracer", &SimPM.ntracer);
  p            = p011;
  p->critical  = true;
  params.push_back(p);
  pm_string *p012a = new pm_string("chem_code", &SimPM.chem_code, "CHEM_CODE");
  p                = p012a;
  p->critical      = false;
  params.push_back(p);
  have_setup_tracers = false;  // so we know to populate list later.
  // LEGACY: look for a string parameter called trtype, in case the
  // data file was written with a pre-2018 version of PION.
  pm_string *p312 = new pm_string("tracer_str", &SimPM.TRTYPE);
  p               = p312;
  p->critical     = false;
  params.push_back(p);

  //
  // hydro solver parameters.
  //
  pm_int *p013 = new pm_int("solver", &SimPM.solverType);
  p            = p013;
  p->critical  = true;
  params.push_back(p);
  pm_int *p014 = new pm_int("coord_sys", &SimPM.coord_sys);
  p            = p014;
  p->critical  = true;
  params.push_back(p);
  pm_int *p015 = new pm_int("Space_OOA", &SimPM.spOOA, 2);
  p            = p015;
  p->critical  = false;
  params.push_back(p);
  pm_int *p016 = new pm_int("Time_OOA", &SimPM.tmOOA, 2);
  p            = p016;
  p->critical  = false;
  params.push_back(p);
  pm_double *p017 = new pm_double("Gamma", &SimPM.gamma);
  p               = p017;
  p->critical     = true;
  params.push_back(p);
  pm_double *p018 = new pm_double("CFL", &SimPM.CFL);
  p               = p018;
  p->critical     = true;
  params.push_back(p);
  pm_int *p019 = new pm_int("art_visc", &SimPM.artviscosity);
  p            = p019;
  p->critical  = true;
  params.push_back(p);
  pm_double *p020 = new pm_double("eta_visc", &SimPM.etav);
  p               = p020;
  p->critical     = true;
  params.push_back(p);

  pm_dvararr *p120 = new pm_dvararr("Ref_Vector", SimPM.RefVec);
  p                = p120;
  p->critical      = true;
  params.push_back(p);

  //
  // PHYSICS FLAGS
  //
  pm_int *p021 = new pm_int("EP_dynamics", &SimPM.EP.dynamics);
  p            = p021;
  p->critical  = true;
  params.push_back(p);
  pm_int *p022 = new pm_int("EP_raytracing", &SimPM.EP.raytracing);
  p            = p022;
  p->critical  = true;
  params.push_back(p);
  pm_int *p023 = new pm_int("EP_cooling", &SimPM.EP.cooling);
  p            = p023;
  p->critical  = true;
  params.push_back(p);
  pm_int *p024 = new pm_int("EP_chemistry", &SimPM.EP.chemistry);
  p            = p024;
  p->critical  = true;
  params.push_back(p);
  pm_int *p025 = new pm_int("EP_coll_ionisation", &SimPM.EP.coll_ionisation);
  p            = p025;
  p->critical  = true;
  params.push_back(p);
  pm_int *p026 = new pm_int("EP_phot_ionisation", &SimPM.EP.phot_ionisation);
  p            = p026;
  p->critical  = true;
  params.push_back(p);
  pm_int *p027 =
      new pm_int("EP_rad_recombination", &SimPM.EP.rad_recombination);
  p           = p027;
  p->critical = true;
  params.push_back(p);
  pm_int *p028 = new pm_int("EP_update_erg", &SimPM.EP.update_erg);
  p            = p028;
  p->critical  = true;
  params.push_back(p);
  pm_int *p029 =
      new pm_int("EP_MP_timestep_limit", &SimPM.EP.MP_timestep_limit, 0);
  p           = p029;
  p->critical = false;
  params.push_back(p);

  pm_double *pMNT =
      new pm_double("EP_Min_Temperature", &SimPM.EP.MinTemperature, 0.0);
  p           = pMNT;
  p->critical = false;
  params.push_back(p);

  pm_double *pMXT =
      new pm_double("EP_Max_Temperature", &SimPM.EP.MaxTemperature, 1.0e100);
  p           = pMXT;
  p->critical = false;
  params.push_back(p);

  //
  // Hydrogen abundance (by mass) X.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double *pXXX =
      new pm_double("EP_Hydrogen_MassFrac", &SimPM.EP.H_MassFrac, 0.7154);
  p           = pXXX;
  p->critical = false;
  params.push_back(p);

  //
  // Helium abundance (by mass) Y.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double *pYYY =
      new pm_double("EP_Helium_MassFrac", &SimPM.EP.Helium_MassFrac, 0.2703);
  p           = pYYY;
  p->critical = false;
  params.push_back(p);
  //
  // Metal abundance (by mass) Z.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double *pZZZ =
      new pm_double("EP_Metal_MassFrac", &SimPM.EP.Metal_MassFrac, 0.0142);
  p           = pZZZ;
  p->critical = false;
  params.push_back(p);

  pm_int *pSTC =
      new pm_int("EP_sat_thermal_cond", &SimPM.EP.sat_thermal_cond, 0);
  p           = pSTC;
  p->critical = false;
  params.push_back(p);

  pm_double *pTCs = new pm_double("EP_tc_strength", &SimPM.EP.tc_strength, 1.0);
  p               = pTCs;
  p->critical     = false;
  params.push_back(p);

  //
  // TIMESTEPS
  //
  pm_double *p030 = new pm_double("t_start", &SimPM.starttime);
  p               = p030;
  p->critical     = true;
  params.push_back(p);
  pm_double *p031 = new pm_double("t_finish", &SimPM.finishtime);
  p               = p031;
  p->critical     = true;
  params.push_back(p);
  pm_int *p032 = new pm_int("t_step", &SimPM.timestep, 0);
  p            = p032;
  p->critical  = false;
  params.push_back(p);
  pm_double *p033 = new pm_double("t_sim", &SimPM.simtime, 0.0);
  p               = p033;
  p->critical     = false;
  params.push_back(p);
  pm_double *p034 =
      new pm_double("min_timestep", &SimPM.min_timestep, TINYVALUE);
  p           = p034;
  p->critical = false;
  params.push_back(p);
  pm_double *p035 = new pm_double("last_dt", &SimPM.last_dt, 0.0);
  p               = p035;
  p->critical     = false;
  params.push_back(p);

  //
  // OUTPUT PARAMETERS
  //
  pm_int *p040 = new pm_int("typeofop", &SimPM.typeofop);
  p            = p040;
  p->critical  = true;
  params.push_back(p);
  pm_string *p041 = new pm_string("outfile", &SimPM.outFileBase);
  p               = p041;
  p->critical     = true;
  params.push_back(p);
  pm_int *p042 = new pm_int("op_freq", &SimPM.opfreq, 0);
  p            = p042;
  p->critical  = false;
  params.push_back(p);
  pm_double *p043 = new pm_double("opfreq_time", &SimPM.opfreq_time, -1.0);
  p               = p043;
  p->critical     = false;
  params.push_back(p);
  pm_int *p044 = new pm_int("op_criterion", &SimPM.op_criterion, 0);
  p            = p044;
  p->critical  = false;
  params.push_back(p);

  //
  // JET SIMULATION PARAMS (GET WIDTH/STAT LATER)
  //
  pm_int *p050 = new pm_int("JetSim", &JP.jetic, 0);
  p            = p050;
  p->critical  = false;
  params.push_back(p);
  have_setup_jet_pm = false;  // so we know to populate list later.

  //
  // RAY-TRACING PARAMS (GET STRENGTH/POSN LATER)
  //
  pm_int *p060 = new pm_int("RT_Nsources", &(SimPM.RS.Nsources), 0);
  p            = p060;
  p->critical  = false;
  params.push_back(p);
  have_setup_rt_src = false;  // so we know to populate list later.

  //
  // STELLAR WIND PARAMS (GET STRENGTH/POSN LATER)
  //
  pm_int *w200 = new pm_int("WIND_Nsources", &SWP.Nsources, 0);
  p            = w200;
  p->critical  = false;
  params.push_back(p);
  have_setup_windsrc = false;  // so we know to populate list later.

  //
  // UNITS (NOTE THESE ARE JUST FOR INFORMATION PURPOSES NOW!)
  //
  pm_string *u001 = new pm_string("unitsys", &uc.unitsys, "CGS");
  p               = u001;
  p->critical     = false;
  params.push_back(p);
  pm_string *u002 = new pm_string("unitdens", &uc.density, "g.cm-3");
  p               = u002;
  p->critical     = false;
  params.push_back(p);
  pm_string *u003 = new pm_string("unitlen", &uc.length, "cm");
  p               = u003;
  p->critical     = false;
  params.push_back(p);
  pm_string *u004 = new pm_string("unitvel", &uc.velocity, "cm.s-1");
  p               = u004;
  p->critical     = false;
  params.push_back(p);
  pm_string *u005 = new pm_string("unitmagf", &uc.bfield, "Gauss/sqrt(4pi)");
  p               = u005;
  p->critical     = false;
  params.push_back(p);
  //
  // For the values, I think the idea is that one code unit equals
  // this many units in the current unit system.  But I can't remember.
  //
  pm_double *u006 = new pm_double("rhoval", &uc.rhoVal, 1.0);
  p               = u006;
  p->critical     = false;
  params.push_back(p);
  pm_double *u007 = new pm_double("lenval", &uc.lenVal, 1.0);
  p               = u007;
  p->critical     = false;
  params.push_back(p);
  pm_double *u008 = new pm_double("velval", &uc.velVal, 1.0);
  p               = u008;
  p->critical     = false;
  params.push_back(p);
  pm_double *u009 = new pm_double("magval", &uc.magVal, 1.0);
  p               = u009;
  p->critical     = false;
  params.push_back(p);

  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::clear_param_list(std::list<class pm_base *> &listptr)
{
  if (listptr.empty())
    return;
  else {
    do {
      list<pm_base *>::iterator i = listptr.begin();
      pm_base *p                  = (*i);
      p->set_ptr(0);
      delete p;
      listptr.erase(i);
    } while (!listptr.empty());
  }
  return;
}



// ##################################################################
// ##################################################################



int DataIOBase::read_simulation_parameters(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  if (params.empty())
    spdlog::error(
        "{}: {}", "Parameter list is empty -- make sure it populates itself!!",
        0);

  //
  // loop over all parameters.
  //
  pm_base *p = 0;
  int err    = 0;
  for (list<pm_base *>::iterator iter = params.begin(); iter != params.end();
       ++iter) {
    p   = (*iter);
    err = read_header_param(p);
    if (err) {
      if (p->critical) {
        spdlog::error("{}: {}", "Error reading parameter", p->name);
      }
      else {
        // cout <<"parameter "<<p->name<<" not found. setting to default
        // val.\n";
        p->set_to_default();
        err = 0;
      }
    }
  }

  //
  // Read boundary conditions for each edge and internal BCs
  // Two options: lecacy code uses BC_STRING, new code uses individual
  // parameters.  Test for legacy code first.
  //
  if (SimPM.BC_STRING != "") {
    int i = 0;
    string::size_type pos;
    string temp;

    // First get the external boundaries
    string d[6] = {"XN", "XP", "YN", "YP", "ZN", "ZP"};
    string *par[6];
    par[0] = &(SimPM.BC_XN);
    par[1] = &(SimPM.BC_XP);
    par[2] = &(SimPM.BC_YN);
    par[3] = &(SimPM.BC_YP);
    par[4] = &(SimPM.BC_ZN);
    par[5] = &(SimPM.BC_ZP);
    for (i = 0; i < 2 * SimPM.ndim; i++) {
      if ((pos = SimPM.BC_STRING.find(d[i])) == string::npos)
        spdlog::error("{}: {}", "Couldn't find boundary condition for ", d[i]);
      else {
        temp = SimPM.BC_STRING.substr(pos + 2, 3);
        if (temp == "per")
          *par[i] = "periodic";
        else if (temp == "out")
          *par[i] = "outflow";
        else if (temp == "owo")
          *par[i] = "one-way-outflow";
        else if (temp == "inf")
          *par[i] = "inflow";
        else if (temp == "ref")
          *par[i] = "reflecting";
        else if (temp == "jrf")
          *par[i] = "equator-reflect";
        else if (temp == "fix")
          *par[i] = "fixed";
        else if (temp == "dmr")
          *par[i] = "DMR";
        else if (temp == "sb1")
          *par[i] = "SB1";
        else
          spdlog::error("{}: {}", "Unrecognised BC type", SimPM.BC_STRING);
      }
    }  // loop over external boundaries
    // cout <<SimPM.BC_XN <<"  "<<SimPM.BC_XP <<"  ";
    // cout <<SimPM.BC_YN <<"  "<<SimPM.BC_YP <<"  ";
    // cout <<SimPM.BC_ZN <<"  "<<SimPM.BC_ZP <<"\n";

    // Now look for any internal boundaries
    int len = SimPM.BC_STRING.length();
    len     = (len + 5) / 6;
    if (len > 2 * SimPM.ndim) {
      SimPM.BC_Nint = len - 2 * SimPM.ndim;
      if (!SimPM.BC_INT)
        SimPM.BC_INT = mem.myalloc(SimPM.BC_INT, SimPM.BC_Nint);
      // cout <<" got "<<SimPM.BC_Nint<<" internal boundaries\n";
    }
    for (i = 2 * SimPM.ndim; i < len; i++) {
      // cout <<"i="<<i<<", len="<<len<<"\n";
      if ((pos = SimPM.BC_STRING.find("IN", i * 6)) == string::npos) {
        spdlog::error(
            "{}: {}", "internal boundary condition not found", SimPM.BC_STRING);
      }
      else {
        temp = SimPM.BC_STRING.substr(pos + 2, 3);
        if (temp == "jet")
          SimPM.BC_INT[i - 2 * SimPM.ndim] = "jet";
        else if (temp == "dm2")
          SimPM.BC_INT[i - 2 * SimPM.ndim] = "DMR2";
        else if (temp == "rsh")
          SimPM.BC_INT[i - 2 * SimPM.ndim] = "RadShock";
        else if (temp == "rs2")
          SimPM.BC_INT[i - 2 * SimPM.ndim] = "RadShock2";
        else if (temp == "wnd")
          SimPM.BC_INT[i - 2 * SimPM.ndim] = "stellar-wind";
        else
          spdlog::error("{}: {}", "Unrecognised INT BC type", SimPM.BC_STRING);
      }
    }  // loop over internal boundaries

  }  // if OLD LEGACY BC STRING

  else {
    // New boundary conditions format
    if (!have_setup_bc_pm) set_bc_pm_params(SimPM);
    if (bc_pm.empty())
      spdlog::error("{}: {}", "Boundary parameter list is empty!!", 0);
    //
    // now read them:
    //
    int ct = 0;
    for (list<pm_base *>::iterator iter = bc_pm.begin(); iter != bc_pm.end();
         ++iter) {
      p   = (*iter);
      err = read_header_param(p);
      // cout <<"boundary list "<<ct<<", parameter "<<p->name<<"\n";
      if (err) spdlog::error("{}: {}", "Error reading parameter", p->name);
      ct++;
    }
  }

  //
  // We now use num_tracer to set the position of the first tracer:
  //
  if (SimPM.ntracer > 0)
    SimPM.ftr = SimPM.nvar - SimPM.ntracer;
  else
    SimPM.ftr = SimPM.nvar;
  //
  // Set up tracer parameters, based on ntracer and read them in
  //
  if (!have_setup_tracers) set_tracer_params(SimPM);

  //
  // Two options: first legacy code, which uses a string
  // called trtype to set all tracers, or new code, which
  // uses chem_code and then a list called TracerIJK.
  //
  if (SimPM.TRTYPE != "") {
    SimPM.chem_code = SimPM.TRTYPE.substr(0, 6);
    int len         = (SimPM.TRTYPE.length() + 5) / 6 - 1;
    if (len != SimPM.ntracer)
      spdlog::error("{}: {}", "bad tracer string (LEGACY)", SimPM.TRTYPE);
    for (int i = 0; i < len; i++) {
      SimPM.tracers[i] = SimPM.TRTYPE.substr(6 * (i + 1), 6);
      // cout <<"tracer["<<i<<"] = "<<SimPM.tracers[i] <<"\n";
    }

  }  // if LEGACY tracer variables used.

  else {
    for (list<pm_base *>::iterator iter = tr_pm.begin(); iter != tr_pm.end();
         ++iter) {
      p   = (*iter);
      err = read_header_param(p);
      if (err) spdlog::error("{}: {}", "Error reading parameter", p->name);
    }
  }

  //
  // Set SimPM.Range[] based on Xmin,Xmax
  //
  for (int i = 0; i < SimPM.ndim; i++)
    SimPM.Range[i] = SimPM.Xmax[i] - SimPM.Xmin[i];

  //
  //
  // Read Jet parameters if doing a JET SIM
  //
  if (JP.jetic) {
    if (!have_setup_jet_pm) set_jet_pm_params();
    if (jet_pm.empty())
      spdlog::error(
          "{}: {}",
          "Jet parameter list is empty -- make sure it populates itself!!", 0);
    //
    // now read them:
    //
    for (list<pm_base *>::iterator iter = jet_pm.begin(); iter != jet_pm.end();
         ++iter) {
      p   = (*iter);
      err = read_header_param(p);
      if (err) spdlog::error("{}: {}", "Error reading parameter", p->name);
    }
  }

  //
  // Now the radiation sources -- only read in the sources if we
  // haven't already done so!  This can happen when analysing multiple
  // timesteps.  First check if there are radiation sources --
  // SimPM.RS.Nsources was already assigned by reading the header data from
  // the input file. THIS CODE ASSUMES THE NUMBER OF SOURCES DOES NOT CHANGE
  // OVER TIME!!
  //
  // cout <<"RT-Nsrcs = "<<SimPM.RS.Nsources<<"\n";
  if (SimPM.RS.Nsources > 0) {
    //
    // First set up the list of parameters we need to read:
    //
    if (!have_setup_rt_src || SimPM.RS.sources.empty()) {
      have_setup_rt_src = false;
      clear_param_list(rt_src);
      set_rt_src_params(SimPM);
    }
    if (rt_src.empty()) {
      spdlog::error(
          "{}: {}",
          "RT-src parameter list is empty -- make sure it populates itself!!",
          0);
    }

    //
    // Now read properties for each source and populate SimPM.RS.sources
    // vector. Note the rt_src list has all parameters for Nsources
    // radiation sources, so Nsources*Nparams elements.
    //
    for (list<pm_base *>::iterator iter = rt_src.begin(); iter != rt_src.end();
         ++iter) {
      p   = (*iter);
      err = read_header_param(p);
      if (err && p->critical) {
        spdlog::error("{}: {}", "Error reading parameter", p->name);
      }
    }

    //
    // Now we have N sources in SimPM.RS.sources.  First assign an id
    // to each of them, and make sure their data got assigned.
    //
    for (int i = 0; i < SimPM.RS.Nsources; i++) {
      SimPM.RS.sources.at(i).id = i;
      if (SimPM.RS.sources.at(i).type < 0) {
        spdlog::error("{}: {}", "Failed to get source type for source id", i);
      }
      if (SimPM.RS.sources.at(i).at_infinity < 0) {
        spdlog::error(
            "{}: {}", "Failed to get source_at_infty for source id", i);
      }
      if (SimPM.RS.sources.at(i).pos[XX] < -9.0e200) {
        spdlog::error(
            "{}: {}", "Failed to get source position for source id", i);
      }
      if (SimPM.RS.sources.at(i).EvoFile == "NOFILE") {
        // cout <<"\t|+|+|+|+|+| Non-evolving radiation source "<<i<<"
        // found.\n";
      }
      else {
        // cout <<"\t|+|+|+|+|+| Evolving radiation source "<<i<<"
        // detected = "; cout <<SimPM.RS.sources.at(i).EvoFile<<"\n";
      }
      if (SimPM.RS.sources.at(i).opacity_var + SimPM.ftr >= SimPM.nvar) {
        spdlog::error(
            "{}: {}",
            "Opacity var for source is off end of array (ftr offset!)", i);
      }
      if (SimPM.RS.sources.at(i).NTau < 1) {
        spdlog::error("{}: {}", "Failed to set source NTau for source id", i);
      }
    }

  }  // read RT data.

  //
  // Now read in wind parameters, if doing a wind simulation
  // THIS CODE ASSUMES THE NUMBER OF SOURCES DOES NOT CHANGE OVER TIME!!
  //
  if (SWP.Nsources > 0) {
    if (!have_setup_windsrc) set_windsrc_params();
    if (windsrc.empty())
      spdlog::error(
          "{}: {}",
          "wind-src parameter list is empty -- make sure it populates itself!!",
          0);
    //
    // Now read each property for each wind source and add the source
    // Need to read from disk to temp variables first.
    //
    list<pm_base *>::iterator iter = windsrc.begin();
    //
    // Problem: for data analysis, we read the wind source from each output
    // in turn, so we can't simply add it every time without deleting all
    // the elements beforehand.
    //
    struct stellarwind_params *temp_wind = 0;
    while (SWP.params.size() > 0) {
      temp_wind = SWP.params.back();
      SWP.params.pop_back();
      mem.myfree(temp_wind);
    }

    for (int isw = 0; isw < SWP.Nsources; isw++) {
      // double Mdot, rad, posn[MAX_DIM], Vinf, Tw, trcr[MAX_NVAR], Rstar;
      // int type;
      ostringstream nm;
      struct stellarwind_params *wind = 0;
      wind                            = mem.myalloc(wind, 1);
      wind->id                        = isw;

      nm.str("");
      nm << "WIND_" << isw << "_posn";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(posn));
      (*iter)->set_ptr(static_cast<void *>(wind->dpos));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_velocity";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(wind->velocity));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_mass";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&rad));
      (*iter)->set_ptr(static_cast<void *>(&wind->Mass));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_rad_";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&rad));
      (*iter)->set_ptr(static_cast<void *>(&wind->radius));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Mdot";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Mdot));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Vinf";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Vinf));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Vrot";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Vrot));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->Vrot<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Tw__";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Tstar));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Rstr";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Rstar));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_Bstr";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->Bstar));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_type";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->type));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->type<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_trcr";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(wind->tr));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      //
      // zero extra variables to avoid uninitialised data...
      // tracers are stored in the first ntracer values of array.
      //
      for (int v = SimPM.ntracer; v < MAX_NVAR; v++)
        wind->tr[v] = 0.0;
      ++iter;
      // cout<<nm<<"\n";

      //
      // New stuff for evolving winds: data-file, time-offset,
      // update-frequency.
      //
      nm.str("");
      nm << "WIND_" << isw << "_evofile";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->evolving_wind_file));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->evolving_wind_file<<"\n";

      //
      // Whether to enhance Mdot based on Omega (over and above what
      // the evolutionary code does).
      //
      nm.str("");
      nm << "WIND_" << isw << "_enhance_mdot";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->enhance_mdot));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "failed to find WIND_{}_enhance_mdot parameter, setting to 0", isw);
      }
      ++iter;
      // cout<<nm.str()<<" = "<<wind->enhance_mdot<<"\n";

      //
      // Value of xi, power-law index of term that focusses stellar
      // wind onto the equator for rotating stars.
      /// (1-Omega*sin(theta))^xi
      /// See e.g. Langer, Garcia-Segura & Mac Low (1999,ApJ,520,L49).
      //
      nm.str("");
      nm << "WIND_" << isw << "_xi";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->xi));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "failed to find WIND_{}_xi parameter, setting to -0.43", isw);
        wind->xi = -0.43;
        // spdlog::error("{}: {}", "Error reading parameter",(*iter)->name);
      }
      ++iter;
      // cout<<nm.str()<<" = "<<wind->xi<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_t_offset";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->time_offset));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->time_offset<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_updatefreq";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->update_freq));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->update_freq<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_t_scalefac";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->t_scalefactor));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of 1.0",
            (*iter)->name);
        wind->t_scalefactor = 1.0;
        err                 = 0;
      }
      ++iter;
      // cout<<nm.str()<<" = "<<wind->t_scalefactor<<"\n";
      // Test for moving source
      nm.str("");
      nm << "WIND_" << isw << "_moving_star";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->moving_star));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of -1.",
            (*iter)->name);
        wind->moving_star = -1;
        err               = 0;
      }
      // cout<<nm.str()<<" = "<<wind->moving_star<<"\n";
      iter++;

      nm.str("");
      nm << "WIND_" << isw << "_eccentricity";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->eccentricity));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of 1.0",
            (*iter)->name);
        wind->eccentricity = 1.0;
        err                = 0;
      }
      // cout<<nm.str()<<" = "<<wind->eccentricity<<"\n";
      iter++;

      nm.str("");
      nm << "WIND_" << isw << "_orbital_period";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->OrbPeriod));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of 0.0",
            (*iter)->name);
        wind->OrbPeriod = 0.0;
        err             = 0;
      }
      // cout<<nm.str()<<" = "<<wind->OrbPeriod<<"\n";
      iter++;

      nm.str("");
      nm << "WIND_" << isw << "_periastron_vec_x";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->PeriastronX));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of 0.0",
            (*iter)->name);
        wind->PeriastronX = 0.0;
        err               = 0;
      }
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_periastron_vec_y";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(&wind->PeriastronY));
      err = read_header_param(*iter);
      if (err) {
        spdlog::debug(
            "Error reading parameter {} setting to default value of 0.0",
            (*iter)->name);
        wind->PeriastronY = 0.0;
        err               = 0;
      }
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_velocity";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(wind->velocity));
      err = read_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      // Now we should have got all the sources, so add the source to
      // the global list.
      SWP.params.push_back(wind);

    }  // loop over sources

    //
    // Check we got all the sources:
    //
    // if (SWP.Nsources != SW.Nsources())
    //   spdlog::error("{}: {}", "Got the wrong number of stellar wind
    //   sources!",
    //     SWP.Nsources-SW.Nsources());
    if (SWP.Nsources != static_cast<int>(SWP.params.size())) {
      spdlog::debug(
          "Num wind srcs={}, SWP.params.size()={}", SWP.Nsources,
          SWP.params.size());
      spdlog::error(
          "{}: {}", "Got the wrong number of stellar wind sources!",
          SWP.Nsources - static_cast<int>(SWP.params.size()));
    }
  }  // read WIND data.

  //
  // Finally run some checks to make sure parameters are sane
  //
  err += check_header_parameters(SimPM);

  return err;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_tracer_params(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  if (have_setup_tracers) {
    spdlog::error(
        "{}: {}", "Trying to setup tracer parameters twice!",
        have_setup_tracers);
  }
  if (SimPM.ntracer == 0) {
    // cout <<"\t *** No tracers to set up\n";
    have_setup_tracers = true;
    return;
  }
  //
  // So now we have tracers and we need to add them to the list.
  //
  if (!SimPM.tracers) {
    SimPM.tracers = mem.myalloc(SimPM.tracers, SimPM.ntracer);
  }

  class pm_base *p;

  for (int i = 0; i < SimPM.ntracer; i++) {
    ostringstream temp;
    temp << "Tracer";
    temp.width(3);
    temp.fill('0');
    temp << i;
    // cout <<"i="<<i<<", setting up tracer : "<<temp.str()<<"\n";

    pm_string *ptemp =
        new pm_string(temp.str(), &SimPM.tracers[i], "NEED_TRACER_VALUES!");
    p           = ptemp;
    p->critical = true;
    tr_pm.push_back(p);
  }
  have_setup_tracers = true;
  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_windsrc_params()
{
  //
  // Sanity checks!
  //
  if (have_setup_windsrc || !windsrc.empty() || !SWP.Nsources) {
    spdlog::debug(
        "FLAG: {} EMPTY?: {} Wind-SIM?: {}", have_setup_windsrc,
        windsrc.empty(), SWP.Nsources);
    spdlog::error(
        "{}: {}", "WARNING! why set wind parameters?!", have_setup_windsrc);
  }

  //
  // Loop over sources and add strength+position for each
  //
  for (int n = 0; n < SWP.Nsources; n++) {
    ostringstream temp01;
    temp01.str("");
    temp01 << "WIND_" << n << "_posn";
    ostringstream temp24;
    temp24.str("");
    temp24 << "WIND_" << n << "_velocity";
    ostringstream temp23;
    temp23.str("");
    temp23 << "WIND_" << n << "_mass";
    ostringstream temp02;
    temp02.str("");
    temp02 << "WIND_" << n << "_rad_";
    ostringstream temp03;
    temp03.str("");
    temp03 << "WIND_" << n << "_Mdot";
    ostringstream temp04;
    temp04.str("");
    temp04 << "WIND_" << n << "_Vinf";
    ostringstream temp09;
    temp09.str("");
    temp09 << "WIND_" << n << "_Vrot";
    ostringstream temp05;
    temp05.str("");
    temp05 << "WIND_" << n << "_Tw__";
    ostringstream temp08;
    temp08.str("");
    temp08 << "WIND_" << n << "_Rstr";
    ostringstream temp16;
    temp16.str("");
    temp16 << "WIND_" << n << "_Bstr";
    ostringstream temp06;
    temp06.str("");
    temp06 << "WIND_" << n << "_type";
    ostringstream temp07;
    temp07.str("");
    temp07 << "WIND_" << n << "_trcr";
    //
    // New stuff for evolving winds:
    //
    ostringstream temp11;
    temp11.str("");
    temp11 << "WIND_" << n << "_evofile";
    ostringstream temp13;
    temp13.str("");
    temp13 << "WIND_" << n << "_enhance_mdot";
    ostringstream temp14;
    temp14.str("");
    temp14 << "WIND_" << n << "_xi";
    ostringstream temp15;
    temp15.str("");
    temp15 << "WIND_" << n << "_t_offset";
    ostringstream temp10;
    temp10.str("");
    temp10 << "WIND_" << n << "_updatefreq";
    ostringstream temp12;
    temp12.str("");
    temp12 << "WIND_" << n << "_t_scalefac";
    ostringstream temp21;
    temp21.str("");
    temp21 << "WIND_" << n << "_moving_star";
    ostringstream temp17;
    temp17.str("");
    temp17 << "WIND_" << n << "_eccentricity";
    ostringstream temp18;
    temp18.str("");
    temp18 << "WIND_" << n << "_orbital_period";
    ostringstream temp19;
    temp19.str("");
    temp19 << "WIND_" << n << "_periastron_vec_x";
    ostringstream temp20;
    temp20.str("");
    temp20 << "WIND_" << n << "_periastron_vec_y";
    ostringstream temp22;
    temp22.str("");
    temp22 << "WIND_" << n << "_velocity";


    pm_ddimarr *w001 = new pm_ddimarr(temp01.str());  // position of source (cm)
    windsrc.push_back(w001);

    pm_ddimarr *w024 =
        new pm_ddimarr(temp24.str());  // velocity of source (cm/s)
    windsrc.push_back(w024);

    pm_double *w023 = new pm_double(temp23.str());  // mass of star (Msun)
    w023->critical  = false;
    windsrc.push_back(w023);

    pm_double *w002 = new pm_double(temp02.str());  // radius of wind BC (cm)
    windsrc.push_back(w002);

    pm_double *w003 = new pm_double(temp03.str());  // Mdot (Msun/yr)
    windsrc.push_back(w003);

    pm_double *w004 = new pm_double(temp04.str());  // v_inf (km/s)
    windsrc.push_back(w004);

    pm_double *w009 = new pm_double(temp09.str());  // v_rot (km/s)
    w009->critical  = false;
    windsrc.push_back(w009);

    pm_double *w005 = new pm_double(temp05.str());  // Twind (K)
    windsrc.push_back(w005);

    pm_double *w008 = new pm_double(temp08.str());  // radius of star (cm)
    windsrc.push_back(w008);

    pm_double *w016 = new pm_double(temp16.str());  // B-field of star (G)
    w016->critical  = false;
    windsrc.push_back(w016);

    pm_int *w006 = new pm_int(temp06.str());  // wind type flag
    windsrc.push_back(w006);

    pm_dvararr *w007 = new pm_dvararr(temp07.str());  // tracers
    windsrc.push_back(w007);

    // wind-evolution file.
    pm_string *w011 = new pm_string(temp11.str());
    w011->critical  = false;
    windsrc.push_back(w011);

    // enhance mdot based on rotation?
    pm_int *w013   = new pm_int(temp13.str());
    w013->critical = false;
    windsrc.push_back(w013);

    // power-law index xi for rotating stars
    pm_double *w014 = new pm_double(temp14.str());
    w014->critical  = false;
    windsrc.push_back(w014);

    // time offset
    pm_double *w015 = new pm_double(temp15.str());
    double dv       = 0.0;
    w015->critical  = false;
    w015->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w015);

    // update frequency (in years)
    pm_double *w010 = new pm_double(temp10.str());
    w010->critical  = false;
    dv              = 1000.0;
    w010->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w010);

    // scale factor (default must be 1, parameter must not be critical).
    pm_double *w012 = new pm_double(temp12.str());
    w012->critical  = false;
    dv              = 1.0;
    w012->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w012);

    // is star moving?  1=yes, 0=no
    pm_int *w021   = new pm_int(temp21.str());
    w021->critical = false;
    int *mv        = new int;
    *mv            = 0;
    w021->set_default_val(static_cast<void *>(mv));
    windsrc.push_back(w021);

    // eccentricity (default must be 1, parameter must not be critical).
    pm_double *w017 = new pm_double(temp17.str());
    w017->critical  = false;  // dv=1.0;
    // w017->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w017);

    // Orbital Period (default must be 0, parameter must not be critical).
    pm_double *w018 = new pm_double(temp18.str());
    w018->critical  = false;  // dv=1.0;
    // w018->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w018);

    // Periastron vectror (default must be 1, parameter must not be critical).
    pm_double *w019 = new pm_double(temp19.str());
    w019->critical  = false;  // dv=1.0;
    // w019->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w019);

    // Periastron vectror (default must be 1, parameter must not be critical).
    pm_double *w020 = new pm_double(temp20.str());
    w020->critical  = false;  // dv=1.0;
    // w020->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w020);

    double vel[MAX_DIM] = {0.0, 0.0, 0.0};
    pm_ddimarr *w022 = new pm_ddimarr(temp22.str());  // velocity of source (cm)
    w022->critical   = false;
    w022->set_default_val(static_cast<void *>(vel));
    windsrc.push_back(w022);
  }
  have_setup_windsrc = true;
  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_rt_src_params(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  //
  // Sanity checks!
  //
  if (have_setup_rt_src || !rt_src.empty() || !SimPM.RS.Nsources) {
    spdlog::debug(
        "FLAG: {} EMPTY?: {} RT-SIM?: {}", have_setup_rt_src, rt_src.empty(),
        SimPM.RS.Nsources);
    spdlog::error(
        "{}: {}", "WARNING! why set RT parameters?!", have_setup_rt_src);
  }

  //
  // First make sure that SimPM.RS.sources has the required number of
  // elements. If we are reading data, this should be empty, and we add
  // Nsources elements. If we are writing, it should have been set up before,
  // so we shouldn't need to add any elements.
  //
  if (SimPM.RS.sources.empty()) {
#ifdef RT_TESTING
    spdlog::debug(
        " DataIOBase::set_rt_src_params() setting up SimPM.RS.sources: {}  {}",
        SimPM.RS.sources.size(), SimPM.RS.Nsources);
#endif
    for (int n = 0; n < SimPM.RS.Nsources; n++) {
      struct rad_src_info temp;

      for (int v = 0; v < MAX_DIM; v++)
        temp.pos[v] = -1.0e200;
      temp.strength    = -1.e200;
      temp.Rstar       = -1.0e200;
      temp.Tstar       = 1.0e200;
      temp.id          = n;
      temp.type        = -1;
      temp.update      = -1;
      temp.at_infinity = -1;
      temp.effect      = -1;
      temp.NTau        = 1;
      temp.opacity_src = -1;
      temp.opacity_var = -1;

      SimPM.RS.sources.push_back(temp);
      if (SimPM.RS.sources.size() != static_cast<unsigned int>(n + 1)) {
        spdlog::error(
            "{}: {}",
            "Radiation source list size error, DataIOBase::set_rt_src_params",
            n);
      }
    }
  }
  if (SimPM.RS.sources.size() != static_cast<unsigned int>(SimPM.RS.Nsources)) {
    spdlog::error("{}: {}", "wrong no. of srcs.", SimPM.RS.sources.size());
  }

  //
  // Loop over sources and add strength, position, type, and (bool) location
  // for each
  //
#ifdef RT_TESTING
  spdlog::debug("DataIOBase::set_rt_src_params() Nsrc={}", SimPM.RS.Nsources);
#endif  // RT_TESTING
  for (int n = 0; n < SimPM.RS.Nsources; n++) {
    ostringstream temp2;
    temp2.str("");
    temp2 << "RT_position_" << n << "_";
    ostringstream temp3;
    temp3.str("");
    temp3 << "RT_strength_" << n;
    ostringstream temp4;
    temp4.str("");
    temp4 << "RT_src_type_" << n;
    ostringstream temp5;
    temp5.str("");
    temp5 << "RT_at_infty_" << n;
    ostringstream temp6;
    temp6.str("");
    temp6 << "RT_update___" << n;
    ostringstream temp7;
    temp7.str("");
    temp7 << "RT_Opacity__" << n;
    // -------------------------------------------------------------
    // include this because of typo in old version: Remove later.
    ostringstream tmp17;
    tmp17.str("");
    tmp17 << "RT_Opactiy__" << n;
    // -------------------------------------------------------------
    ostringstream temp8;
    temp8.str("");
    temp8 << "RT_Tau_var__" << n;
    ostringstream temp9;
    temp9.str("");
    temp9 << "RT_effect___" << n;
    ostringstream tmp10;
    tmp10.str("");
    tmp10 << "RT_Rstar____" << n;
    ostringstream tmp11;
    tmp11.str("");
    tmp11 << "RT_Tstar____" << n;
    ostringstream tmp12;
    tmp12.str("");
    tmp12 << "RT_EVO_FILE_" << n;
    ostringstream tmp13;
    tmp13.str("");
    tmp13 << "RT_Nbins____" << n;

// ADD SOURCE_EFFECT VARIABLE
#ifdef RT_TESTING
    spdlog::info("Adding Source Effect Variables");
#endif
    pm_ddimarr *rtpos =
        new pm_ddimarr(temp2.str(), SimPM.RS.sources[n].pos.data());
    rt_src.push_back(rtpos);
    pm_double *rtstr =
        new pm_double(temp3.str(), &(SimPM.RS.sources[n].strength));
    rt_src.push_back(rtstr);
    pm_int *rttyp = new pm_int(temp4.str(), &(SimPM.RS.sources[n].type));
    rt_src.push_back(rttyp);
    pm_int *rtinf = new pm_int(temp5.str(), &(SimPM.RS.sources[n].at_infinity));
    rt_src.push_back(rtinf);
    pm_int *rtupd = new pm_int(temp6.str(), &(SimPM.RS.sources[n].update));
    rt_src.push_back(rtupd);
    pm_int *rttsc = new pm_int(temp7.str(), &(SimPM.RS.sources[n].opacity_src));
    rttsc->critical = false;
    rt_src.push_back(rttsc);
    // -------------------------------------------------------------
    // include this because of typo in old version: Remove later.
    pm_int *rtts2 = new pm_int(tmp17.str(), &(SimPM.RS.sources[n].opacity_src));
    rtts2->critical = false;
    rt_src.push_back(rtts2);
    // -------------------------------------------------------------
    pm_int *rttvr = new pm_int(temp8.str(), &(SimPM.RS.sources[n].opacity_var));
    rt_src.push_back(rttvr);
    pm_int *rteff = new pm_int(temp9.str(), &(SimPM.RS.sources[n].effect));
    rt_src.push_back(rteff);
    pm_double *rtRst =
        new pm_double(tmp10.str(), &(SimPM.RS.sources[n].Rstar), 0.0);
    rtRst->critical = false;
    rt_src.push_back(rtRst);
    pm_double *rtTst =
        new pm_double(tmp11.str(), &(SimPM.RS.sources[n].Tstar), 0.0);
    rtTst->critical = false;
    rt_src.push_back(rtTst);
    pm_string *rtEvo =
        new pm_string(tmp12.str(), &(SimPM.RS.sources[n].EvoFile), "NOFILE");
    rtEvo->critical = false;
    rt_src.push_back(rtEvo);
    pm_int *rtNBn   = new pm_int(tmp13.str(), &(SimPM.RS.sources[n].NTau), 1);
    rtNBn->critical = false;
    rt_src.push_back(rtNBn);
  }
  have_setup_rt_src = true;
  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_jet_pm_params()
{
  if (have_setup_jet_pm || !jet_pm.empty() || !JP.jetic) {
    spdlog::debug(
        "FLAG: {} EMPTY?: {} JET-SIM?: {}", have_setup_jet_pm, jet_pm.empty(),
        JP.jetic);
    spdlog::error(
        "{}: {}", "WARNING! why set jet parameters?!", have_setup_jet_pm);
  }

  pm_int *p051 = new pm_int("JetRadius", &JP.jetradius);
  jet_pm.push_back(p051);
  pm_dvararr *p052 = new pm_dvararr("JetState", JP.jetstate);
  jet_pm.push_back(p052);

  have_setup_jet_pm = true;
  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_bc_pm_params(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  // cout <<"setting up BC parameters\n";
  if (have_setup_bc_pm || !bc_pm.empty()) {
    spdlog::debug("BCs FLAG: {} EMPTY?: {}\n", have_setup_bc_pm, bc_pm.empty());
    spdlog::error("{}: {}", "set bc parameters twice?!", have_setup_bc_pm);
  }

  pm_string *p007 = new pm_string("BC_XN", &SimPM.BC_XN);
  p007->critical  = true;
  bc_pm.push_back(p007);
  pm_string *p117 = new pm_string("BC_XP", &SimPM.BC_XP);
  p117->critical  = true;
  bc_pm.push_back(p117);
  pm_string *p118 = new pm_string("BC_YN", &SimPM.BC_YN);
  bc_pm.push_back(p118);
  pm_string *p119 = new pm_string("BC_YP", &SimPM.BC_YP);
  bc_pm.push_back(p119);
  pm_string *p122 = new pm_string("BC_ZN", &SimPM.BC_ZN);
  bc_pm.push_back(p122);
  pm_string *p123 = new pm_string("BC_ZP", &SimPM.BC_ZP);
  bc_pm.push_back(p123);

  //
  // Read internal boundary data
  //
  if (SimPM.BC_Nint > 0) {
    if (!SimPM.BC_INT) SimPM.BC_INT = mem.myalloc(SimPM.BC_INT, SimPM.BC_Nint);
  }
  for (int v = 0; v < SimPM.BC_Nint; v++) {
    ostringstream intbc;
    intbc.str("");
    intbc << "BC_INTERNAL_";
    intbc.width(3);
    intbc.fill('0');
    intbc << v;
    pm_string *p124 = new pm_string(intbc.str(), &(SimPM.BC_INT[v]));
    bc_pm.push_back(p124);
  }
  have_setup_bc_pm = true;
  return;
}



// ##################################################################
// ##################################################################



int DataIOBase::write_simulation_parameters(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  if (params.empty())
    spdlog::error(
        "{}: {}", "Parameter list is empty -- make sure it populates itself!!",
        0);

  //
  // loop over all *normal* parameters, writing one-by-one.
  //
  pm_base *p = 0;
  int err    = 0;
  for (list<pm_base *>::iterator iter = params.begin(); iter != params.end();
       ++iter) {
    p   = (*iter);
    err = write_header_param(p);
    if (err) spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
  }

  //
  // Write Boundary string parameters
  //
  if (!have_setup_bc_pm) set_bc_pm_params(SimPM);
  if (bc_pm.empty())
    spdlog::error("{}: {}", "bc_pm empty, Need boundaries!", "Huh?");
  //
  // now write them:
  //
  for (list<pm_base *>::iterator iter = bc_pm.begin(); iter != bc_pm.end();
       ++iter) {
    p = (*iter);
    // cout <<p->name<<"    "; p->show_val(); cout <<"\n";
    err = write_header_param(p);
    if (err) spdlog::error("{}: {}", "Error writing BC parameter", p->name);
  }

  //
  // Write tracer parameters
  //
  if (!have_setup_tracers) set_tracer_params(SimPM);
  // cout <<"Writing tracer names.\n";
  for (list<pm_base *>::iterator iter = tr_pm.begin(); iter != tr_pm.end();
       ++iter) {
    p = (*iter);
    // cout <<"tracer val: "; p->show_val(); cout <<"\n";
    err = write_header_param(p);
    if (err) spdlog::error("{}: {}", "Error writing tracer parameter", p->name);
  }

  //
  // Write Jet parameters if doing a JET SIM
  //
  if (JP.jetic) {
    if (!have_setup_jet_pm) set_jet_pm_params();
    for (list<pm_base *>::iterator iter = jet_pm.begin(); iter != jet_pm.end();
         ++iter) {
      p   = (*iter);
      err = write_header_param(p);
      if (err) spdlog::error("{}: {}", "Error writing Jet parameter", p->name);
    }
  }

  //
  // Now the radiation sources
  //
  if (SimPM.RS.Nsources > 0) {
#ifndef NDEBUG
    spdlog::debug(
        "WRITING Radiation Source PARAMETERS, Nsrc={}", SimPM.RS.Nsources);
#endif
    //
    // Check we have rt_src parameters
    //
    if (!have_setup_rt_src) set_rt_src_params(SimPM);
    if (rt_src.empty()) {
      spdlog::error(
          "{}: {}", "rt_src list empty, but running RT sim!?", "HMMM");
    }
    //
    // Data should be already set up and ready to go.
    //
    // cout <<"rt_src.size = "<<rt_src.size()<<"\n";
    for (list<pm_base *>::iterator iter = rt_src.begin(); iter != rt_src.end();
         ++iter) {
      p   = (*iter);
      err = write_header_param(p);
      if (err)
        spdlog::error(
            "{}: {}", "Error writing Radiation source parameter", p->name);
    }
  }  // write RT data.

  // --------------------------------------
  // ---- Now the stellar wind sources ----
  // --------------------------------------
  if (SWP.Nsources > 0) {
    //
    // Check we have windsrc parameters
    //
    if (!have_setup_windsrc) set_windsrc_params();
    if (windsrc.empty())
      spdlog::error(
          "{}: {}", "windsrc list empty, but we are apparently ouputting \
                 a wind sim!?",
          "HMMM");

    //
    // Now write each property from WS class.
    // Need to copy to a temp variable for the silo interface.
    // List elements are ordered src 0,1,2.
    //
    list<pm_base *>::iterator iter = windsrc.begin();

    for (int isw = 0; isw < SWP.Nsources; isw++) {
      // double xd[MAX_NVAR];
      // int xi[MAX_NVAR];
      ostringstream nm;
      //
      // Make sure the pm name is right, then set the pointer and write
      // the data for each property of source number isw
      //
      nm.str("");
      nm << "WIND_" << isw << "_posn";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->dpos));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_velocity";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->velocity));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_mass";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Mass));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_rad_";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->radius));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Mdot";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Mdot));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Vinf";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Vinf));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Vrot";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Vrot));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Tw__";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Tstar));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Rstr";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Rstar));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_Bstr";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Bstar));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_type";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->type));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_trcr";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->tr));
      for (int v = SimPM.ntracer; v < MAX_NVAR; v++) {
        SWP.params[isw]->tr[v] = 0.0;
      }
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing wIND parameter", (*iter)->name);
      ++iter;
      //
      // New stuff for evolving winds: data-file, time-offset,
      // update-frequency.
      //
      nm.str("");
      nm << "WIND_" << isw << "_evofile";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(
          static_cast<void *>(&SWP.params[isw]->evolving_wind_file));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->evolving_wind_file<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_enhance_mdot";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->enhance_mdot));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->enhance_mdot<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_xi";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->xi));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->xi<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_t_offset";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->time_offset));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->time_offset<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_updatefreq";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->update_freq));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error reading parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->update_freq<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_t_scalefac";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->t_scalefactor));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm<<" = "<<wind->t_scalefactor<<"\n";

      // Test for moving source
      nm.str("");
      nm << "WIND_" << isw << "_moving_star";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->moving_star));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->moving_star<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_eccentricity";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->eccentricity));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;
      // cout<<nm.str()<<" = "<<wind->eccentricity<<"\n";

      nm.str("");
      nm << "WIND_" << isw << "_orbital_period";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->OrbPeriod));
      err = write_header_param(*iter);
      if (err) {
        spdlog::error("Error writing parameter {}", (*iter)->name);
      }
      // cout<<nm.str()<<" = "<<wind->OrbPeriod<<"\n";
      iter++;

      nm.str("");
      nm << "WIND_" << isw << "_periastron_vec_x";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(posn));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->PeriastronX));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_periastron_vec_y";
      if ((*iter)->name.compare(nm.str()) != 0)
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(posn));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->PeriastronY));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing parameter", (*iter)->name);
      ++iter;

      nm.str("");
      nm << "WIND_" << isw << "_velocity";
      if ((*iter)->name.compare(nm.str()) != 0) {
        spdlog::error(
            "{}: {}", "Stellar wind parameters not ordered as expected!",
            (*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->velocity));
      err = write_header_param(*iter);
      if (err)
        spdlog::error("{}: {}", "Error writing WIND parameter", (*iter)->name);
      ++iter;
    }  // loop over sources
  }    // write Stellar wind data.

  return err;
}



// ##################################################################
// ##################################################################



int DataIOBase::check_header_parameters(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  // This is where we check that nvar,ndim,ntracer,outfile, etc. are all set
  // to something sensible.
  if (SimPM.simtime < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) simtime not set, setting to "
        "zero. THIS MAY BE A BIG ERROR",
        1., SimPM.simtime);
    SimPM.simtime = 0.;
  }
  if (SimPM.starttime < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) startime <0, PROBABLY AN ERROR!", 0.,
        SimPM.starttime);
    SimPM.starttime = 0.;
  }
  if (SimPM.finishtime <= 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) finishtime not set.", 1.,
        SimPM.finishtime);
    SimPM.finishtime = -1.;
  }
  if (SimPM.timestep < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) timestep <0, PROBABLY AN ERROR!", 1,
        SimPM.timestep);
    SimPM.timestep = 0;
  }
  if (SimPM.NG[0] < 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) NG[0]<0 -- must not have read "
        "from file",
        SimPM.NG[0]);
  if (SimPM.ndim > 1)
    if (SimPM.NG[1] < 0)
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) NG[1]<0 -- must not have "
          "read from file",
          SimPM.NG[1]);
  if (SimPM.ndim > 2)
    if (SimPM.NG[2] < 0)
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) NG[2]<0 -- must not have "
          "read from file",
          SimPM.NG[2]);

  if (SimPM.Ncell < 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) Ncell<0 -- Error reading from file",
        SimPM.Ncell);

  if (SimPM.Xmin[0] < -1.e50)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) Xmin[0]<0 -- Error reading "
        "from file",
        SimPM.Xmin[0]);
  if (SimPM.ndim > 1)
    if (SimPM.Xmin[1] < -1e50)
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) Xmin[1]<0 -- must not have "
          "read from file",
          SimPM.Xmin[1]);
  if (SimPM.ndim > 2)
    if (SimPM.Xmin[2] < -1e50)
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) Xmin[2]<0 -- must not have "
          "read from file",
          SimPM.Xmin[2]);

  if (SimPM.Xmax[0] <= SimPM.Xmin[0])
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) Xmax[0]<0 -- Error reading "
        "from file",
        SimPM.Xmax[0]);
  if (SimPM.ndim > 1)
    if (SimPM.Xmax[1] <= SimPM.Xmin[1])
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) Xmax[1]<0 -- must not have "
          "read from file",
          SimPM.Xmax[1]);
  if (SimPM.ndim > 2)
    if (SimPM.Xmax[2] <= SimPM.Xmin[2])
      spdlog::error(
          "{}: {}",
          "(Dataio::check_header_parameters) Xmax[2]<0 -- must not have "
          "read from file",
          SimPM.Xmax[2]);

  if (SimPM.spOOA != OA1 && SimPM.spOOA != OA2) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) Order of Accuracy not set "
        "from file. Default to second order",
        2, SimPM.spOOA);
    SimPM.spOOA = OA2;
  }
  if (SimPM.tmOOA != OA1 && SimPM.tmOOA != OA2) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) Order of Accuracy not set "
        "from file. Default to second order",
        2, SimPM.tmOOA);
    SimPM.tmOOA = OA2;
  }

  //  cout <<"(Dataio::check_header_parameters) spOOA: "<<SimPM.spOOA<<"
  //  tmOOA:
  //  "<<SimPM.tmOOA<<"\n";
  if (SimPM.gamma < 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) gamma<0 -- must not have read "
        "it from file",
        SimPM.gamma);
  if (SimPM.CFL < 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) CFL<0   -- must not have read "
        "it from file",
        SimPM.CFL);
  if (SimPM.artviscosity < 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) artviscosity<0 -- must not "
        "have read from file",
        SimPM.artviscosity);
  if (SimPM.etav < 0 && SimPM.artviscosity > 0)
    spdlog::error(
        "{}: {}",
        "(Dataio::check_header_parameters) etav<0 -- must not have read "
        "from file",
        SimPM.etav);
  else if (SimPM.etav < 0 && SimPM.artviscosity < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) etav<0 but art.visc. not "
        "enabled so it's ok.",
        static_cast<double>(SimPM.artviscosity), SimPM.etav);
    SimPM.etav = 0.15;
  }
  if (SimPM.typeofop < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) type of output file not "
        "set, setting to Fits",
        2, SimPM.typeofop);
    SimPM.typeofop = 2;
  }

  if (SimPM.opfreq < 0) {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) Output frequency not set.  "
        "Setting to 50",
        100, SimPM.opfreq);
    SimPM.opfreq = 50;
  }
  if (SimPM.outFileBase == "") {
    spdlog::warn(
        "{}: Expected {} but got {}",
        "(Dataio::check_header_parameters) output file base not set, "
        "setting to ./noname",
        1, 1);
    SimPM.outFileBase = "noname";
  }
  SimPM.maxtime = false;
  SimPM.dt      = 0.;

  if (SimPM.EP.dynamics != 1) {
    spdlog::warn(
        "WARNING: DYNAMICS SWITCHED OFF!!! *****************************");
  }

  //
  // Check that Helium_MassFrac and H_MassFrac are set.
  //
  if (SimPM.EP.H_MassFrac < 0.001 || SimPM.EP.H_MassFrac > 1.0) {
    spdlog::warn("H_MassFrac = {}, resetting to 0.7154", SimPM.EP.H_MassFrac);
    SimPM.EP.H_MassFrac = 0.7154;
  }
  if ((SimPM.EP.H_MassFrac + SimPM.EP.Helium_MassFrac) > 1.1
      || (SimPM.EP.H_MassFrac + SimPM.EP.Helium_MassFrac) < 0.9) {
    spdlog::warn(
        "H_MassFrac  = {}\nWarning: He_MassFrac = {}resetting to 0.7154 and 0.2846, respectively",
        SimPM.EP.H_MassFrac, SimPM.EP.H_MassFrac);
    SimPM.EP.H_MassFrac      = 0.7154;
    SimPM.EP.Helium_MassFrac = 0.2846;
  }
  if (SimPM.EP.Metal_MassFrac < 0.0 || SimPM.EP.Metal_MassFrac > 1.0) {
    spdlog::warn(
        "SimPM.EP.Metal_MassFrac = {} !! resetting to 0.0142",
        SimPM.EP.Metal_MassFrac);
    SimPM.EP.Metal_MassFrac = 0.0142;
  }

  return 0;
}



// -----------------------------------------------------
// --------- BASE DATA I/O CLASS DEFINITIONS -----------
// -----------------------------------------------------

// ##################################################################
// ##################################################################
