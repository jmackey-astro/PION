/// \file setup_fixed_grid.cpp
///
/// \brief Class for setting up grids.
///
/// \author Jonathan Mackey
///
///
/// Modifications:
/// - 2015.02.09 JM: Split sim_control class into a setup class and
///   a derived class for running simulations.
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"

#include "microphysics/microphysics_base.h"

#ifdef LEGACY_CODE
#include "microphysics/MPv0.h"
#include "microphysics/MPv1.h"
#include "microphysics/MPv2.h"
#include "microphysics/MPv4.h"
#endif

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/MPv9.h"
#endif

#include "microphysics/mp_only_cooling.h"

#ifndef EXCLUDE_MPV3
#include "microphysics/MPv3.h"
#endif

#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"
#include "microphysics/MPv8.h"

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif

#include "raytracing/raytracer_SC.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif  // if FITS

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <time.h>
using namespace std;

// ##################################################################
// ##################################################################

setup_fixed_grid::setup_fixed_grid()
{
  FVI_nheat = FVI_nion = 0;
  FVI_heating_srcs.clear();
  FVI_ionising_srcs.clear();
  FVI_need_column_densities_4dt = false;
  spatial_solver                = 0;
  dataio                        = 0;
  textio                        = 0;
}

// ##################################################################
// ##################################################################

setup_fixed_grid::~setup_fixed_grid()
{
#ifdef TESTING
  cout << "(setup_fixed_grid::Destructor) ..."
       << "\n";
#endif
  if (MP) {
    delete MP;
    MP = 0;
  }
  if (spatial_solver) {
    delete spatial_solver;
    spatial_solver = 0;
  }
#ifdef TESTING
  cout << "(setup_fixed_grid::Destructor) Done."
       << "\n";
#endif
}

// ##################################################################
// ##################################################################

void setup_fixed_grid::setup_cell_extra_data(
    class SimParams& SimPM  ///< pointer to simulation parameters
)
{
  //
  // Cells can need extra data for ray-tracing optical depths, eta-values for
  // the H-correction or div(v) for some time-updates and/or viscosity
  // corrections.
  //

  int hc_flag = 0, dv_flag = 0, gp_flag = 0;
  if (SimPM.artviscosity == AV_LAPIDUS || SimPM.eqntype == EQEUL_EINT) {
    // Need one var. for Div(v)
    dv_flag = 1;
  }
  if (SimPM.solverType == FLUX_RS_HLLD) {
    dv_flag = 1;  // Need DivV for shock switch in HLLD.
    gp_flag = 1;
  }
  if (SimPM.artviscosity == AV_HCORRECTION
      || SimPM.artviscosity == AV_HCORR_FKJ98 || SimPM.eqntype == EQEUL_EINT) {
    //
    // need one var for each dimension here.  For H-correction they
    // are for the eta values.  For EQEUL_EINT we need von
    // Neunamm-Richtmeyer viscosity which needs storage for the diagonal
    // Q-values along each axis.
    //
    hc_flag = SimPM.ndim;
  }

  CI.setup_extra_data(SimPM.RS, hc_flag, dv_flag, gp_flag);
  return;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::setup_grid(
    vector<class GridBaseClass*>& g,  ///< grid pointers.
    class SimParams& SimPM            ///< pointer to simulation parameters
)
{
  cout << "(pion ug)  Setting up computational grid\n";

#ifdef TESTING
  cout << "Init::setup_grid: \n";
#endif  // TESTING
  class GridBaseClass** grid = &(g[0]);

  if (SimPM.ndim < 1 || SimPM.ndim > 3)
    rep.error("Only know 1D,2D,3D methods!", SimPM.ndim);

    //
    // Nbc is the depth of the boundary layer.
    //
#ifdef TESTING
  cout << "Setting number of boundary cells == spatial OOA: ";
  cout << SimPM.spOOA << "\n";
#endif  // TESTING
  if (SimPM.spOOA == OA2) {
    SimPM.Nbc    = 2;
    SimPM.Nbc_DD = 2;
  }
  else if (SimPM.spOOA == OA1) {
    SimPM.Nbc    = 1;
    SimPM.Nbc_DD = 1;
  }
  else
    rep.error(
        "Spatial order of accuracy unhandled by boundary conditions!",
        SimPM.spOOA);

  // Force Nbc=1 if using Lax-Friedrichs flux.
  if (SimPM.solverType == FLUX_LF) {
    SimPM.spOOA = SimPM.tmOOA = OA1;
    SimPM.Nbc                 = 1;
    SimPM.Nbc_DD              = 1;
  }

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set Cell dx in cell interface class, and also xmin.
  //
  double dx((SimPM.Xmax[XX] - SimPM.Xmin[XX]) / SimPM.NG[XX]);
  CI.set_nlevels(dx, SimPM.grid_nlevels);  // sets dx on all levels.
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);

  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout << "(setup_fixed_grid::setup_grid) Setting up grid...\n";
#endif
  if (*grid) rep.error("Grid already set up!", *grid);

  if (SimPM.coord_sys == COORD_CRT)
    *grid = new UniformGrid(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin,
        SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax, SimPM.Xmin, SimPM.Xmax);
  else if (SimPM.coord_sys == COORD_CYL)
    *grid = new uniform_grid_cyl(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin,
        SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax, SimPM.Xmin, SimPM.Xmax);
  else if (SimPM.coord_sys == COORD_SPH)
    *grid = new uniform_grid_sph(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin,
        SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax, SimPM.Xmin, SimPM.Xmax);
  else
    rep.error("Bad Geometry in setup_grid()", SimPM.coord_sys);

  if (*grid == 0)
    rep.error("(setup_fixed_grid::setup_grid) Couldn't assign data!", *grid);
#ifdef TESTING
  cout << "(setup_fixed_grid::setup_grid) Done. &grid=" << grid;
  cout << ", and grid=" << *grid << "\n";
  cout << "DX = " << (*grid)->DX() << "\n";
  dp.grid = (*grid);
  cout << "------------------------------------------------------\n\n";
#endif

  return (0);
}  // setup_grid()

// ##################################################################
// ##################################################################

int setup_fixed_grid::setup_microphysics(
    class SimParams& SimPM  ///< pointer to simulation parameters
)
{
  // cout <<"------------------------------------------------------------\n";
  // cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  // cout <<"------------------------------------------------------------\n";
  cout << "(pion)  Setting up microphysics\n";
  //
  // Setup Microphysics class, if needed.
  // First see if we want the only_cooling class (much simpler), and if
  // not then check for the one of the bigger microphysics classes.
  //
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout << "\tRequested cooling but no chemistry... setting";
    cout << " up mp_only_cooling() class. \n";
    cout << "\tTimestep limit = " << SimPM.EP.MP_timestep_limit << "\n";
    MP = new mp_only_cooling(
        SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
    if (!MP) rep.error("mp_only_cooling() init", MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    string mptype;
    mptype           = SimPM.chem_code;
    bool have_set_MP = false;
    cout << "setting up MP type: " << mptype << "\n";

#ifdef LEGACY_CODE
    if (mptype == "MPv0") {
      cout << "\tsetting up MPv0 module\n";
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv0(
          SimPM.nvar, SimPM.ntracer, SimPM.chem_code, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS));
      if (SimPM.EP.MP_timestep_limit < 0 || SimPM.EP.MP_timestep_limit > 5)
        rep.error("BAD dt LIMIT", SimPM.EP.MP_timestep_limit);
      have_set_MP = true;
    }

    if (mptype == "MPv1") {
      cout << "\tsetting up MPv1 microphysics module\n";
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv1(
          SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      cout << "\t**---** WARNING, THIS MODULE HAS BEEN SUPERSEDED BY MPv4. "
              "**--**\n";
      have_set_MP = true;
    }

    if (mptype == "MPv2") {
      cout << "\tsetting up MPv2 module\n";
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv2(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS));
      SimPM.EP.MP_timestep_limit = 1;
      have_set_MP                = true;
    }

    if (mptype == "MPv4") {
      cout << "\tsetting up MPv4 module\n";
#if MPV4_DTLIMIT >= 5 && MPV4_DTLIMIT <= 12
      // cout <<"\t******* N.B. dt05-12 Timestep limiting is enforced by
      // #def"; cout <<" DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 5;
#elif MPV4_DTLIMIT >= 0 && MPV4_DTLIMIT <= 4
      // cout <<"\t******* N.B. dt00-04 Timestep limiting is enforced by
      // #def"; cout <<" MPV4_DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 4;
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv4(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP = true;
    }

    if (mptype == "MPv8") {
      cout << "\tsetting up MPv8 module\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv8(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP = true;
    }

#endif  // LEGACY_CODE

#ifndef EXCLUDE_HD_MODULE
    if (mptype == "MPv9") {
      cout << "\tsetting up microphysics_lowz module\n";
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new microphysics_lowz(
          SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      have_set_MP = true;
    }
#endif  // exclude Harpreet's module

    if (mptype == "MPv3" || mptype == "MPv3__") {
      cout << "\tsetting up MPv3 module\n";
#if MPV3_DTLIMIT >= 0 && MPV4_DTLIMIT <= 12
      // cout <<"\t******* N.B. Timestep limiting is enforced by #def";
      // cout <<" MPV3_DTLIMIT="<<MPV3_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif

      MP = new MPv3(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      // if (SimPM.EP.MP_timestep_limit != 1)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP = true;
    }

    if (mptype == "MPv5" || mptype == "MPv5__") {
      cout << "\tsetting up MPv5 module\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv5(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP = true;
    }

    if (mptype == "MPv6" || mptype == "MPv6__") {
      cout << "\tsetting up MPv6 module\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv6(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP = true;
    }

    if (mptype == "MPv7" || mptype == "MPv7__") {
      cout << "\tsetting up MPv7 module\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new MPv7(
          SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers,
          &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP = true;
    }

#ifdef CODE_EXT_HHE
    if (mptype == "MPv10") {
      cout << "\tsetting up MPv10 module\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised", mptype);
      MP = new mpv9_HHe(
          SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), SimPM.gamma);
      have_set_MP = true;
    }
#endif

    if (!MP) rep.error("microphysics init", MP);
    if (!have_set_MP) rep.error("have_set_MP", have_set_MP);
  }
  else {
    cout << "\tno microphysics.\n";
    MP = 0;
  }

  //
  // If we have a multifrequency ionising source, we can set its properties
  // here. We can only have one of these, so safe to just loop through...
  //
  int err = 0;
  double data[MAX_TAU];  // temp var not used
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type == RT_SRC_SINGLE
        && SimPM.RS.sources[isrc].effect == RT_EFFECT_MFION && MP != 0
        && SimPM.RS.sources[isrc].EvoFile == "NONE") {
      err = MP->set_multifreq_source_properties(&SimPM.RS.sources[isrc], data);
    }
  }
  if (err) rep.error("Setting multifreq source properties", err);

  // cout <<"------------------------------------------------------------\n";
  // cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  // cout <<"------------------------------------------------------------\n";
  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::setup_raytracing(
    class SimParams& SimPM,    ///< pointer to simulation parameters
    class GridBaseClass* grid  ///< pointer to grid
)
{
  if (!SimPM.EP.raytracing) {
    return 0;
  }
  cout << "(pion)  Setting up raytracing on level\n";

  if (!MP) rep.error("can't do raytracing without microphysics", MP);
#ifdef RT_TESTING
  cout << "\n------ RAYTRACER SETUP STARTING --------------\n";
#endif
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays = true;
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++)
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays = false;
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    grid->RT = new raytracer_USC_infinity(
        grid, MP, SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!(grid->RT)) rep.error("init pllel-rays raytracer error", grid->RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    grid->RT = new raytracer_USC(
        grid, MP, SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ftr,
        SimPM.RS.Nsources);
    if (!(grid->RT)) rep.error("init raytracer error 2", grid->RT);
  }

  //
  // Now add the sources to the tracer.  Note that both the implicit
  // and explicit integrators can still only handle a single ionising
  // source, so we do a check for this and bug out if there is more
  // than one.
  //
  int ion_count = 0, uv_count = 0, dif_count = 0;
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++) {

    // see if source is on the domain or not.  Set flag accordingly.
    SimPM.RS.sources[isrc].ongrid = true;
    for (int d = 0; d < SimPM.ndim; d++) {
      if (SimPM.RS.sources[isrc].pos[d] < SimPM.levels[0].Xmin[d]
          || SimPM.RS.sources[isrc].pos[d] > SimPM.levels[0].Xmax[d])
        SimPM.RS.sources[isrc].ongrid = false;
    }

    // source types
    if (SimPM.RS.sources[isrc].type == RT_SRC_SINGLE) {
      //
      // single sources have a flux (if at infinity) or a luminosity (if
      // point sources.
      //
      int s = grid->RT->Add_Source(&(SimPM.RS.sources[isrc]));
#ifdef RT_TESTING
      cout << "Adding IONISING or UV single-source with id: ";
      cout << s << "\n";
#endif
      if (SimPM.RS.sources[isrc].effect == RT_EFFECT_PION_MONO
          || SimPM.RS.sources[isrc].effect == RT_EFFECT_MFION)
        ion_count++;
      else
        uv_count++;
    }  // if ionising source
    else {
      // note that diffuse radiation must be at infinity, and the strength
      // is assumed to be an intensity not a flux, so it is multiplied by
      // a solid angle appropriate to its location in order to get a flux.
      int s = grid->RT->Add_Source(&(SimPM.RS.sources[isrc]));
#ifdef RT_TESTING
      cout << "Adding DIFFUSE radiation source with id: ";
      cout << s << "\n";
#endif
      uv_count++;
      dif_count++;
    }  // if diffuse source
  }    // loop over sources
  if (ion_count > 1) {
    rep.error(
        "Can only have one ionising source for currently implemented method",
        ion_count);
  }
#ifdef RT_TESTING
  cout << "Added " << ion_count << " ionising and " << uv_count
       << " non-ionising";
  cout << " radiation sources, of which " << dif_count
       << " are diffuse radiation.\n";
#endif
  grid->RT->Print_SourceList();

  //
  // Now that we have added all of the sources, we query the raytracer to get
  // all of the source properties into structs for the microphysics calls.
  // NOTE THAT IF THE NUMBER OF SOURCES OR THEIR PROPERTIES CHANGE OVER TIME,
  // I WILL HAVE TO WRITE NEW CODE TO UPDATE THIS!
  //
  FVI_nheat = grid->RT->N_heating_sources();
  FVI_nion  = grid->RT->N_ionising_sources();
  FVI_heating_srcs.resize(FVI_nheat);
  FVI_ionising_srcs.resize(FVI_nion);
  grid->RT->populate_UVheating_src_list(FVI_heating_srcs);
  grid->RT->populate_ionising_src_list(FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (grid->RT->type_of_RT_integration() == RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (
      grid->RT && grid->RT->type_of_RT_integration() == RT_UPDATE_IMPLICIT
      && SimPM.EP.MP_timestep_limit == 5) {
    // For implicit updates to limit by xdot and/or edot
    // The raytracing has not already been done, so we call it here.
    FVI_need_column_densities_4dt = true;
  }
  else {
    FVI_need_column_densities_4dt = false;
  }

#ifdef RT_TESTING
  cout << "------------- RAYTRACER SETUP COMPLETE ----------------\n";
#endif
  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::setup_evolving_RT_sources(
    class SimParams& SimPM  ///< pointer to simulation parameters
)
{
  //
  // Loop through list of sources, and see if any of them have an
  // evolution file (if none, then the string is set to NOFILE in
  // the data I/O stage).
  //
  int Nevo = 0;
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].EvoFile == "NOFILE") {
#ifdef RT_TESTING
      cout << "setup_evolving_RT_sources() Source " << isrc;
      cout << " has no evolution file.\n";
#endif
    }
    else {
      Nevo++;
      struct star istar;
#ifdef RT_TESTING
      cout << "setup_evolving_RT_sources() Source " << isrc;
      cout << " has EvoFile " << istar.file_name << "\n";
#endif
      istar.file_name = SimPM.RS.sources[isrc].EvoFile;
      istar.src_id    = isrc;
      SimPM.STAR.push_back(istar);
    }
  }
  //
  // Now go through each one we found and read the evolution file into arrays
  // and set the rest of the data in the 'star' struct.
  //
  for (int isrc = 0; isrc < Nevo; isrc++) {
    struct star* istar = &(SimPM.STAR[isrc]);
    istar->Nlines      = 0;
    istar->time.resize(0);
    istar->Log_L.resize(0);
    istar->Log_T.resize(0);
    istar->Log_R.resize(0);
    istar->Log_V.resize(0);
    //
    // Open file
    //
    FILE* infile = 0;
    infile       = fopen(istar->file_name.c_str(), "r");
    if (!infile)
      rep.error("can't open wind evolving radiation source file", infile);
    // Skip first two lines
    char line[512];
    char* rval = 0;
    rval       = fgets(line, 512, infile);
    if (!rval) rep.error("setup_fixed_grid, RS, file read 1", line);
    rval = fgets(line, 512, infile);
    if (!rval) rep.error("setup_fixed_grid, RS, file read 2", line);
    // Temporary variables for column values
    // Columns are time, M, L, Teff, Mdot, vrot, vcrit, vinf
    double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0, t6 = 0.0, t7 = 0.0;
    size_t iline = 0;
    while ((rval = fgets(line, 512, infile)) != 0) {
      // printf("rval= %s\n",rval);
      sscanf(
          line, "   %lE   %lE %lE %lE %lE %lE %lE", &t1, &t2, &t3, &t4, &t5,
          &t6, &t7);
      // cout.precision(16);
      // cout <<t1 <<"  "<<t2  <<"  "<< t3  <<"  "<< t4 <<"  "<< t5 <<"
      // "<< t6
      // <<"\n"; cout <<t1 <<"  "<< t3  <<"  "<< t4 <<"  "<< t5 <<"  "<<
      // t6
      // <<"\n";
      istar->Nlines++;
      istar->time.push_back(t1);
      istar->Log_L.push_back(log10(t3 / pconst.Lsun()));
      istar->Log_T.push_back(log10(t4));
      istar->Log_V.push_back(log10(t6 / 1.0e5));

      // Stellar radius, from Stefan Boltzmann Law.
      t6 = sqrt(
          pow(10.0, istar->Log_L[iline]) * pconst.Lsun()
          / (4.0 * pconst.pi() * pconst.StefanBoltzmannConst() * pow(t4, 4.0)));
      istar->Log_R.push_back(log10(t6 / pconst.Rsun()));
      iline++;
      // cout <<istar->Log_L.back()<<"\n";
    }
    fclose(infile);

    //
    // set the last_line counter to be the array index
    // nearest to (but less than) the current time.
    //
    iline = 0;
    while (istar->time[iline] < SimPM.simtime)
      iline++;
    istar->last_line = iline - 1;
    // cout <<"last line="<<iline-1<<"\n";

    // initialise to zero.
    istar->Lnow = istar->Tnow = istar->Rnow = istar->Vnow = 0.0;
  }

  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::update_evolving_RT_sources(
    class SimParams& SimPM,   ///< simulation parameters
    const double time,        ///< current simulation time
    class RayTracingBase* RT  ///< pointer to raytracing class
)
{
  bool updated = false;
  //
  // Loop over all sources with Evolution files.
  //
  for (unsigned int isrc = 0; isrc < SimPM.STAR.size(); isrc++) {
    // cout <<"isrc="<<isrc<<" of "<<SimPM.STAR.size()<<"\n";
    struct star* istar = &(SimPM.STAR[isrc]);
    istar->t_now       = time;
    size_t i           = istar->last_line;

    // Check if we have reached the last line of the file!
    if (i == (istar->Nlines - 1)) {
      cout << "\n\n*#*#* WARNING #*#*#* update_evolving_RT_sources()";
      cout << " Last line, assuming constant Lum from now on!\n\n";
      return 0;
    }
    // Check if we have moved forward one line in table, in which
    // case we need to increment i.
    while (istar->t_now > istar->time[i + 1]) {
      i++;
      istar->last_line = i;
    }

    // The star properties are bracketed by line i and line i+1,
    // so we can do a simple interpolation between them.
    double Lnow, Tnow, Rnow, Vnow;
    Lnow = istar->Log_L[i]
           + (istar->t_now - istar->time[i])
                 * (istar->Log_L[i + 1] - istar->Log_L[i])
                 / (istar->time[i + 1] - istar->time[i]);
    Tnow = istar->Log_T[i]
           + (istar->t_now - istar->time[i])
                 * (istar->Log_T[i + 1] - istar->Log_T[i])
                 / (istar->time[i + 1] - istar->time[i]);
    Rnow = istar->Log_R[i]
           + (istar->t_now - istar->time[i])
                 * (istar->Log_R[i + 1] - istar->Log_R[i])
                 / (istar->time[i + 1] - istar->time[i]);
    Vnow = istar->Log_V[i]
           + (istar->t_now - istar->time[i])
                 * (istar->Log_V[i + 1] - istar->Log_V[i])
                 / (istar->time[i + 1] - istar->time[i]);

    // convert units
    Lnow = exp(pconst.ln10() * (Lnow)) * pconst.Lsun();  // erg/s
    Tnow = exp(pconst.ln10() * (Tnow));                  // K
    Rnow = exp(pconst.ln10() * (Rnow));                  // Rsun
    Vnow = exp(pconst.ln10() * (Vnow));                  // km/s

    // If L or T change by more than 1% then update them; otherwise
    // leave as they are.
    if (fabs(Lnow - istar->Lnow) / istar->Lnow > 0.01
        || fabs(Tnow - istar->Tnow) / istar->Tnow > 0.01) {
      cout << "update_evolving_RT_sources() NOW: t=" << istar->t_now;
      cout << "\tL=" << Lnow;
      cout << "\tT=" << Tnow;
      cout << "\tR=" << Rnow;
      // cout <<"\tV="<< Vnow;
      cout << "\n";
      istar->Lnow = Lnow;
      istar->Tnow = Tnow;
      istar->Rnow = Rnow;
      istar->Vnow = Vnow;

      // Copy new data into SimPM.RS, send updates to RayTracing and
      // microphysics classes.
      struct rad_src_info* rs = &(SimPM.RS.sources[istar->src_id]);
      rs->strength            = istar->Lnow;
      rs->Rstar               = istar->Rnow;
      rs->Tstar               = istar->Tnow;

      // This is a hack, fix it to something sensible
      if (rs->effect == RT_EFFECT_UV_HEATING) {
        rs->strength =
            1.0e48 * (rs->strength / 1.989e38) * exp(-1e4 / rs->Tstar);
        // cout <<"FUV source: strength="<<rs->strength<<"\n";
      }

      RT->update_RT_source_properties(rs);
      updated = true;
    }
    else {
      // cout <<"not updating source: L="<<istar->Lnow<<",
      // T="<<istar->Tnow<<"\n";
    }
  }

  //
  // Get the data back from RT into the structs for MP updates.
  //
  if (updated) {
    RT->populate_UVheating_src_list(FVI_heating_srcs);
    RT->populate_ionising_src_list(FVI_ionising_srcs);
  }

  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::boundary_conditions(
    class SimParams& par,               ///< simulation parameters
    vector<class GridBaseClass*>& grid  ///< grid pointers.
    // class GridBaseClass *grid ///< pointer to grid.
)
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout << "Setting up BCs in Grid with Nbc=" << par.Nbc << "\n";
#endif
  //
  // Choose what BCs to set up based on BC strings.
  //
  int err = setup_boundary_structs(par, grid[0], 0);
  rep.errorTest("sfg::boundary_conditions::sb_structs", 0, err);

  //
  // Ask grid to set up data for external boundaries.
  //
  err = grid[0]->SetupBCs(par);
  rep.errorTest("sfg::boundary_conditions::SetupBCs", 0, err);

#ifdef TESTING
  cout << "(setup_fixed_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::setup_boundary_structs(
    class SimParams& par,       ///< simulation parameters
    class GridBaseClass* grid,  ///< pointer to grid.
    const int)
{
#ifdef TESTING
  cout << "Set BC types...\n";
#endif

  // Set number of boundaries: 2 for each dimension, plus internal.
  int len = 2 * par.ndim + par.BC_Nint;
#ifdef TESTING
  cout << "Got " << len << " boundaries to set up.\n";
#endif
  if (grid->BC_bd.size() == 0) {
#ifdef TESTING
    cout << "size=" << grid->BC_bd.size() << " allocating boundaries\n";
#endif
    for (int b = 0; b < len; b++) {
      struct boundary_data* bd = new boundary_data;
      if (b < 2 * par.ndim)
        bd->depth = grid->boundary_depth(static_cast<enum direction>(b));
      grid->BC_bd.push_back(bd);
    }
  }
  else {
    // assume this has already been called so quit (happens for
    // MPI Nested Grid algorithm)
#ifdef TESTING
    cout << "already set up boundaries, so returning here.\n";
#endif
    return 0;
  }

  //
  // First the 2N external boundaries.  Put the strings into an array.
  //
  std::vector<std::string> bc_strings(2 * par.ndim);
  bc_strings[0] = par.BC_XN;
  bc_strings[1] = par.BC_XP;
  if (par.ndim > 1) {
    bc_strings[2] = par.BC_YN;
    bc_strings[3] = par.BC_YP;
  }
  if (par.ndim > 2) {
    bc_strings[4] = par.BC_ZN;
    bc_strings[5] = par.BC_ZP;
  }

  //
  // Now go through each boundary and assign everything needed.
  //
  int i = 0;
  for (i = 0; i < 2 * par.ndim; i++) {
    grid->BC_bd[i]->dir =
        static_cast<direction>(i);  // XN=0,XP=1,YN=2,YP=3,ZN=4,ZP=5
    grid->BC_bd[i]->ondir = grid->OppDir(grid->BC_bd[i]->dir);
#ifdef TESTING
    cout << "i=" << i << ", dir = " << grid->BC_bd[i]->dir
         << ", ondir=" << grid->BC_bd[i]->ondir << "\n";
#endif
    grid->BC_bd[i]->baxis = static_cast<axes>(i / 2);
    //
    // odd values of i are positive boundaries, others are negative.
    //
    if ((i + 2) % 2 != 0) {
      grid->BC_bd[i]->bloc = grid->Xmax(grid->BC_bd[i]->baxis);
      grid->BC_bd[i]->bpos = true;
    }
    else {
      grid->BC_bd[i]->bloc = grid->Xmin(grid->BC_bd[i]->baxis);
      grid->BC_bd[i]->bpos = false;
    }
    //
    // find boundary condition specified:
    //
    grid->BC_bd[i]->type = bc_strings[i];

    if (grid->BC_bd[i]->type == "periodic") {
      grid->BC_bd[i]->itype = PERIODIC;
      grid->BC_bd[i]->type  = "PERIODIC";
    }
    else if (
        grid->BC_bd[i]->type == "outflow"
        || grid->BC_bd[i]->type == "zero-gradient") {
      grid->BC_bd[i]->itype = OUTFLOW;
      grid->BC_bd[i]->type  = "OUTFLOW";
    }
    else if (grid->BC_bd[i]->type == "one-way-outflow") {
      grid->BC_bd[i]->itype = ONEWAY_OUT;
      grid->BC_bd[i]->type  = "ONEWAY_OUT";
    }
    else if (grid->BC_bd[i]->type == "inflow") {
      grid->BC_bd[i]->itype = INFLOW;
      grid->BC_bd[i]->type  = "INFLOW";
    }
    else if (grid->BC_bd[i]->type == "reflecting") {
      grid->BC_bd[i]->itype = REFLECTING;
      grid->BC_bd[i]->type  = "REFLECTING";
    }
    else if (grid->BC_bd[i]->type == "axisymmetric") {
      grid->BC_bd[i]->itype = AXISYMMETRIC;
      grid->BC_bd[i]->type  = "AXISYMMETRIC";
    }
    else if (grid->BC_bd[i]->type == "equator-reflect") {
      grid->BC_bd[i]->itype = JETREFLECT;
      grid->BC_bd[i]->type  = "JETREFLECT";
    }
    else if (grid->BC_bd[i]->type == "fixed") {
      grid->BC_bd[i]->itype = FIXED;
      grid->BC_bd[i]->type  = "FIXED";
    }
    else if (grid->BC_bd[i]->type == "DMR") {
      grid->BC_bd[i]->itype = DMACH;
      grid->BC_bd[i]->type  = "DMACH";
    }
    else {
      rep.error("Don't know this BC type", grid->BC_bd[i]->type);
    }

    if (!grid->BC_bd[i]->data.empty())
      rep.error(
          "Boundary data not empty in constructor!",
          grid->BC_bd[i]->data.size());
    grid->BC_bd[i]->refval = 0;
#ifdef TESTING
    cout << "\tBoundary type " << i << " is " << grid->BC_bd[i]->type << "\n";
#endif
  }

  if (i < len) {
#ifdef TESTING
    cout << "Got " << i << " boundaries, but have " << len << " boundaries.\n";
    cout << "Must have extra BCs... checking for internal BCs\n";
#endif
    do {
      grid->BC_bd[i]->dir = NO;
      if (par.BC_Nint < i - 2 * par.ndim) {
        rep.error("Bad Number of boundaries", par.BC_Nint);
      }
      else {
        grid->BC_bd[i]->type = par.BC_INT[i - 2 * par.ndim];
      }

      if (grid->BC_bd[i]->type == "jet") {
        grid->BC_bd[i]->itype = JETBC;
        grid->BC_bd[i]->type  = "JETBC";
      }
      else if (grid->BC_bd[i]->type == "DMR2") {
        grid->BC_bd[i]->itype = DMACH2;
        grid->BC_bd[i]->type  = "DMACH2";
      }
      else if (grid->BC_bd[i]->type == "stellar-wind") {
        grid->BC_bd[i]->itype = STWIND;
        grid->BC_bd[i]->type  = "STWIND";
      }
      else {
        rep.error("Don't know this BC type", grid->BC_bd[i]->type);
      }

      if (!grid->BC_bd[i]->data.empty()) {
        rep.error(
            "Boundary data not empty in constructor!",
            grid->BC_bd[i]->data.size());
      }
      grid->BC_bd[i]->refval = 0;
#ifdef TESTING
      cout << "\tBoundary type " << i << " is " << grid->BC_bd[i]->type << "\n";
#endif
      i++;
    } while (i < len);
  }

#ifdef TESTING
  cout << len << " BC structs set up.\n";
#endif
  return 0;
}

// ##################################################################
// ##################################################################

void setup_fixed_grid::setup_dataio_class(
    class SimParams& par,  ///< simulation parameters
    const int typeOfFile   ///< type of I/O: 1=text,2=fits,5=silo
)
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

    case 1:  // Start From ASCII Parameterfile.
      dataio = new dataio_text(par);
      if (!dataio) rep.error("dataio_text initialisation", dataio);
      break;

#ifdef FITS
    case 2:  // Start from FITS restartfile.
      dataio = new DataIOFits(par);
      break;
    case 4:  // fits +ascii
      dataio = new DataIOFits(par);
      textio = new dataio_text(par);
      break;
#endif  // if FITS

#ifdef SILO
    case 5:  // Start from Silo snapshot.
      dataio = new dataio_silo(par, "DOUBLE");
      break;
    case 6:  // silo + text
      dataio = new dataio_silo(par, "DOUBLE");
      textio = new dataio_text(par);
      if (!textio) rep.error("INIT:: textio initialisation", par.typeofop);
      break;
#endif  // if SILO
    default:
      rep.error("setup_fixed_grid::Init unhandled filetype", typeOfFile);
  }
  return;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid::set_equations(
    class SimParams& par  ///< simulation parameters
)
{
  cout << "(pion)  Setting up solver for equations\n";

  if (par.solverType < 0)
    rep.error("set_equations: solverType not set yet.", par.solverType);
  if (par.eqntype <= 0)
    rep.error("set_equations: eqntype not set yet.", par.eqntype);
  if (par.eqnNDim != 3)
    rep.error("set_equations: eqnNDim not set right.", par.eqnNDim);
  if (par.nvar <= 0) rep.error("set_equations: nvar not set yet.", par.nvar);
  if (par.ndim <= 0) rep.error("set_equations: ndim not set yet.", par.ndim);
  if (par.artviscosity < 0)
    rep.error("set_equations: artviscosity not set yet.", par.artviscosity);
  if (par.coord_sys < 0)
    rep.error("set_equations: coordinate system not set.", par.coord_sys);

  if (spatial_solver) {
    delete spatial_solver;
    spatial_solver = 0;
  }

  if (par.coord_sys == COORD_CRT) {
    cout << "\tset_equations() Using Cartesian coord. system.\n";
    switch (par.eqntype) {
      case EQEUL:
        cout << "\tset_equations() Using Euler Equations.\n";
        spatial_solver = new class FV_solver_Hydro_Euler(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQEUL);
        break;
      case EQMHD:
        cout << "\tset_equations() Using Ideal MHD Equations.\n";
        spatial_solver = new class FV_solver_mhd_ideal_adi(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQMHD);
        break;
      case EQGLM:
        cout << "\tset_equations() Using GLM MHD Equations.\n";
        spatial_solver = new class FV_solver_mhd_mixedGLM_adi(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQGLM);
        break;
      case EQFCD:
        cout << "\tset_equations() Using Field-CD MHD Equations.\n";
        rep.error("Field CD got lost in some code updates", EQFCD);
        break;
      default:
        rep.error("Don't know the specified equations...", par.eqntype);
        break;
    }
  }  // cartesian

  else if (par.coord_sys == COORD_CYL) {
    cout << "\tset_equations() Using Cylindrical coord. system.\n";
    switch (par.eqntype) {
      case EQEUL:
        cout << "\tset_equations() Using Euler Equations.\n";
        spatial_solver = new class cyl_FV_solver_Hydro_Euler(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQEUL);
        break;
      case EQMHD:
        cout << "\tset_equations() Using Ideal MHD Equations.\n";
        spatial_solver = new class cyl_FV_solver_mhd_ideal_adi(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQMHD);
        break;
      case EQGLM:
        cout << "\tset_equations() Using GLM MHD Equations.\n";
        spatial_solver = new class cyl_FV_solver_mhd_mixedGLM_adi(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQGLM);
        break;
      default:
        rep.error("not implemented yet for axisymmetry", par.eqntype);
    }
  }  // axisymmetric

  else if (par.coord_sys == COORD_SPH) {
    cout << "\tset_equations() Using Spherical coordinate system.\n";
    switch (par.eqntype) {
      case EQEUL:
        cout << "\tset_equations() Using Euler Equations.\n";
        spatial_solver = new class sph_FV_solver_Hydro_Euler(
            par.nvar, par.ndim, par.CFL, par.gamma, par.RefVec, par.etav,
            par.ntracer);
        if (!spatial_solver)
          rep.error("Couldn't set up solver/equations class.", EQEUL);
        break;
      default:
        rep.error("not implemented yet for spherical", par.eqntype);
    }
  }  // spherically symmetric

  //
  // Check that we set up an equations class!
  //
  if (!spatial_solver)
    rep.error("setup_fixed_grid::set_equations() Failed", spatial_solver);

  // cout <<"------------------------------------------------------\n\n";
  return (0);
}  // set_equations.

// ##################################################################
// ##################################################################
