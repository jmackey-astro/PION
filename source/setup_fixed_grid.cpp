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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"

#include "tools/command_line_interface.h"

#include "sim_control.h"

#include "dataIO/dataio.h"
#include "microphysics/microphysics_base.h"

#ifndef EXCLUDE_MPV1
#include "microphysics/microphysics.h"
#endif 

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/microphysics_lowZ.h"
#endif

#include "microphysics/mp_only_cooling.h"

#ifndef EXCLUDE_MPV2
#ifdef MP_V2_AIFA
#include "microphysics/mp_v2_aifa.h"
#endif
#endif 

#ifndef EXCLUDE_MPV3
#include "microphysics/mp_explicit_H.h"
#endif

#ifndef EXCLUDE_MPV4
#include "microphysics/mp_implicit_H.h"
#endif 

#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"
#include "microphysics/mpv8_StarBench_heatcool.h"

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif


#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_hydro_adi_Eint.h"
#include "spatial_solvers/solver_eqn_hydro_iso.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <climits>
using namespace std;


#define TIMESTEP_FULL 1
#define TIMESTEP_FIRST_PART 2



// ##################################################################
// ##################################################################


setup_fixed_grid::setup_fixed_grid()
{
  FVI_nheat = FVI_nion = 0;
  FVI_heating_srcs.clear();
  FVI_ionising_srcs.clear();
  FVI_need_column_densities_4dt = false;
}



// ##################################################################
// ##################################################################


setup_fixed_grid::~setup_fixed_grid()
{
#ifdef TESTING
  cout << "(setup_fixed_grid::Destructor) Deleting Grid Class..." <<"\n";
#endif
  if (MP)     {delete MP; MP=0;}
  if (RT)     {delete RT; RT=0;}
#ifdef TESTING
  cout << "(setup_fixed_grid::Destructor) Done." <<"\n";
#endif
}




// ##################################################################
// ##################################################################



void setup_fixed_grid::setup_cell_extra_data()
{
  //
  // Cells can need extra data for ray-tracing optical depths, eta-values for the
  // H-correction or div(v) for some time-updates and/or viscosity corrections.
  //

  int hc_flag = 0, dv_flag=0;
  if (SimPM.artviscosity==AV_LAPIDUS ||
      SimPM.eqntype==EQEUL_EINT) {
    // Need one var. for Div(v)
    dv_flag = 1;
  }
  if (SimPM.artviscosity==AV_HCORRECTION ||
      SimPM.artviscosity==AV_HCORR_FKJ98 ||
      SimPM.eqntype==EQEUL_EINT) {
    //
    // need one var for each dimension here.  For H-correction they
    // are for the eta values.  For EQEUL_EINT we need von Neunamm-Richtmeyer
    // viscosity which needs storage for the diagonal Q-values along each axis.
    //
    hc_flag = SimPM.ndim;
  }

  CI.setup_extra_data(SimPM.RS, hc_flag, dv_flag);
  return;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_grid(
      class GridBaseClass **grid,
      class MCMDcontrol * ///< unused for serial code.
      )
{
  cout <<"------------------------------------------------------\n";
  cout <<"--------  Setting up computational grid --------------\n";

#ifdef TESTING
  cout <<"Init::setup_grid: &grid="<< grid<<", and grid="<<*grid<<"\n";
#endif // TESTING

  if (SimPM.gridType!=1) {
    rep.warning("gridType not set correctly: Only know Uniform finite\
                 volume grid, so resetting to 1!",1,SimPM.gridType);
    SimPM.gridType=1;
  }

  if (SimPM.ndim <1 || SimPM.ndim>3)
    rep.error("Only know 1D,2D,3D methods!",SimPM.ndim);

  //
  // Nbc is the depth of the boundary layer.
  //
#ifdef TESTING
  cout <<"Setting number of boundary cells == spatial OOA: ";
  cout <<SimPM.spOOA<<"\n";
#endif // TESTING
  if      (SimPM.spOOA==OA2) SimPM.Nbc = 2;
  else if (SimPM.spOOA==OA1) SimPM.Nbc = 1;
  else rep.error("Spatial order of accuracy unhandled by boundary \
                  conditions!",SimPM.spOOA);
  
  // Force Nbc=1 if using Lax-Friedrichs flux.
  if (SimPM.solverType==FLUX_LF)
  {SimPM.spOOA = SimPM.tmOOA = OA1; SimPM.Nbc=1;}

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  setup_cell_extra_data();

  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
#endif
  if (*grid) rep.error("Grid already set up!",*grid);

  if      (SimPM.coord_sys==COORD_CRT)
    *grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin, SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax);
  else if (SimPM.coord_sys==COORD_CYL)
    *grid = new uniform_grid_cyl (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin, SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax);
  else if (SimPM.coord_sys==COORD_SPH)
    *grid = new uniform_grid_sph (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc, SimPM.Xmin, SimPM.Xmax, SimPM.NG, SimPM.Xmin, SimPM.Xmax);
  else 
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);

  if (*grid==0) rep.error("(setup_fixed_grid::setup_grid) Couldn't assign data!", *grid);
#ifdef TESTING
  cout <<"(UniformFV::setup_grid) Done. &grid="<< grid<<", and grid="<<*grid<<"\n";
  cout <<"DX = "<<(*grid)->DX()<<"\n";
  dp.grid = (*grid);
#endif
  cout <<"------------------------------------------------------\n\n";

  return(0);
} // setup_grid()




// ##################################################################
// ##################################################################


int setup_fixed_grid::boundary_conditions(
      class GridBaseClass *grid 
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"(UniformFV::boundary_conditions)";
  
  if(SimPM.typeofbc=="FIXED") {
    cout << "\t Using fixed boundary conditions; only useful for shock tube."<<"\n";
  }
  else if (SimPM.typeofbc=="PERIODIC") {
    cout <<"\t Using periodic BCs on all sides of grid.\n";
  }
  else if (SimPM.typeofbc=="ABSORBING") {
    cout <<"\t Using absorbing BCs on all sides of grid.\n";
  }
  else if (SimPM.typeofbc=="REFLECTING") {
    cout <<"\t Using reflecting BCs on all sides of grid.\n";
  }
  else {
    cout <<"\t Using the following BCs: "<<SimPM.typeofbc<<"\n";
  }
#endif

#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<SimPM.Nbc<<"\n";
#endif
  int err = grid->SetupBCs(SimPM.Nbc,SimPM.typeofbc);
  if (err) rep.error("setup_fixed_grid::boundary_conditions() Couldn't \
                      set up boundary conditions class.",err);
#ifdef TESTING
  cout <<"(setup_fixed_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_microphysics()
{
  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  //
  // Setup Microphysics class, if needed.
  // First see if we want the only_cooling class (much simpler), and if
  // not then check for the one of the bigger microphysics classes.
  //
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout <<"\t******* Requested cooling but no chemistry... setting";
    cout <<" up mp_only_cooling() class, with timestep-limiting.\n";
    SimPM.EP.MP_timestep_limit = 1;
    MP = new mp_only_cooling(SimPM.nvar, &(SimPM.EP));
    if (!MP) rep.error("mp_only_cooling() init",MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    string mptype;

#ifdef OLD_TRACER

    cout <<"TRTYPE: "<<SimPM.trtype<<"\n";
    if (SimPM.trtype.size() >=6)
      mptype = SimPM.trtype.substr(0,6); // Get first 6 chars for type of MP.
    else mptype = "None";

# else
    
    mptype = SimPM.chem_code;

#endif // OLD_TRACER

    bool have_set_MP=false;


#ifndef EXCLUDE_MPV1
    if      (mptype=="ChAH__" || mptype=="onlyH_") {
      cout <<"\t******* setting up MP_Hydrogen microphysics module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MP_Hydrogen(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      cout <<"\t**---** WARNING, THIS MODULE HAS BEEN SUPERSEDED BY mp_implicit_H. **--**\n";
      have_set_MP=true;
    }
#endif // exclude MPv1


#ifndef EXCLUDE_HD_MODULE
    if (mptype=="lowZ__") {
      cout <<"\t******* setting up microphysics_lowz module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new microphysics_lowz(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }
#endif // exclude Harpreet's module

#ifndef EXCLUDE_MPV2
    if (mptype=="MPv2__") {
#ifdef MP_V2_AIFA
      cout <<"\t******* setting up mp_v2_aifa module *********\n";
      cout <<"\t******* N.B. Timestep limiting is enforced. **\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_v2_aifa(SimPM.nvar, SimPM.ntracer, SimPM.trtype);
      SimPM.EP.MP_timestep_limit = 1;
#else
      rep.error("Enable mp_v2_aifa as an ifdef if you really want to use it",2);
#endif
      have_set_MP=true;
    }
#endif // exclude MPv2


#ifndef EXCLUDE_MPV3
    if (mptype=="MPv3__") {
      cout <<"\t******* setting up mp_explicit_H module *********\n";
#if MPV3_DTLIMIT>=0 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. Timestep limiting is enforced by #def";
      cout <<" MPV3_DTLIMIT="<<MPV3_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif

      MP = new mp_explicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP)
      );
      //if (SimPM.EP.MP_timestep_limit != 1)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv3


#ifndef EXCLUDE_MPV4
    if (mptype=="MPv4__") {
      cout <<"\t******* setting up mp_implicit_H module *********\n";
#if MPV4_DTLIMIT>=5 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. dt05-12 Timestep limiting is enforced by #def";
      cout <<" DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =5;
#elif MPV4_DTLIMIT>=0 && MPV4_DTLIMIT<=4
      cout <<"\t******* N.B. dt00-04 Timestep limiting is enforced by #def";
      cout <<" MPV4_DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =4;
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_implicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      //SimPM.EP.MP_timestep_limit = 4;  // limit by recombination time only
      //if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv4


    if (mptype=="MPv5__") {
      cout <<"\t******* setting up mpv5_molecular module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv5_molecular(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv6__") {
      cout <<"\t******* setting up mpv6_PureH module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv6_PureH(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv7__") {
      cout <<"\t******* setting up mpv7_TwoTempIso module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv7_TwoTempIso(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

#ifdef CODE_EXT_SBII
    if (mptype=="MPSBHC") {
      cout <<"\t******* setting up mpv8_StarBench_heatcool module *********\n";
      cout <<"\t******* This is for StarBench test propblems with heating and cooling.\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv8_SBheatcool(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }
#endif // CODE_EXT_SBII

#ifdef CODE_EXT_HHE
    if (mptype=="MPv9__") {
      cout <<"\t******* setting up mpv9_HHe module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv9_HHe(SimPM.nvar, SimPM.ntracer, SimPM.trtype, 
                        &(SimPM.EP), SimPM.gamma);
      have_set_MP=true;
    }
#endif

#ifndef EXCLUDE_MPV1
    //
    // Finally, if MP has not been set up yet, try to set up the v0
    // microphysics integrator, which is slow, but can model a number
    // of elements and ions.
    //
    if (!have_set_MP) {
      cout <<"\t******* setting up MicroPhysics (v0) module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);

#ifdef OLD_TRACER

      MP = new MicroPhysics(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));

# else

      MP = new MicroPhysics(SimPM.nvar, SimPM.ntracer, SimPM.chem_code, SimPM.trtype, &(SimPM.EP));

#endif // OLD_TRACER

      if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
        rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv1/0

    if (!MP) rep.error("microphysics init",MP);
    if (!have_set_MP) rep.error("HUH? have_set_MP",have_set_MP);
  }
  else {
    cout <<"\t******** not doing microphysics.\n";
    MP=0;
  }

  //
  // If we have a multifrequency ionising source, we can set its properties here.
  // We can only have one of these, so safe to just loop through...
  //
  int err=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE &&
        SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI &&
        MP!=0
        ) {
      err = MP->set_multifreq_source_properties(&SimPM.RS.sources[isrc]);
    }
  }
  if (err) rep.error("Setting multifreq source properties",err);
  

  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_raytracing(
      class GridBaseClass *grid
      )
{
  //
  // If not doing raytracing, return immediately.
  //
  if (!SimPM.EP.raytracing) {
    return 0;
  }

  //
  // Now we are doing raytracing, so set up a raytracer and add sources to it.
  //
  if (!MP) rep.error("can't do raytracing without microphysics",MP);
  cout <<"\n----------------- RAYTRACER SETUP STARTING -----------------------\n";
  RT=0;
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays=true;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++)
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays=false;
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    RT = new raytracer_USC_infinity(grid,MP);
    if (!RT) rep.error("init pllel-rays raytracer error",RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    RT = new raytracer_USC(grid,MP);
    if (!RT) rep.error("init raytracer error 2",RT);
  }

  //
  // Now add the sources to the tracer.  Note that both the implicit and explicit
  // integrators can still only handle a single ionising source, so we do a check
  // for this and bug out if there is more than one.
  //
  int ion_count=0, uv_count=0, dif_count=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) {
      //
      // single sources have a flux (if at infinity) or a luminosity (if point
      // sources.
      //
      cout <<"Adding IONISING or UV single-source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MONO ||
          SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI)
        ion_count++;
      else 
        uv_count++;
    } // if ionising source
    else {
      // note that diffuse radiation must be at infinity, and the strength is assumed to
      // be an intensity not a flux, so it is multiplied by a solid angle appropriate
      // to its location in order to get a flux.
      cout <<"Adding DIFFUSE radiation source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      uv_count++;
      dif_count++;
    } // if diffuse source
  } // loop over sources
  if (ion_count>1) {
    rep.error("Can only have one ionising source for currently implemented method",ion_count);
  }
  cout <<"Added "<<ion_count<<" ionising and "<<uv_count<<" non-ionising";
  cout <<" radiation sources, of which "<<dif_count<<" are diffuse radiation.\n";
  RT->Print_SourceList();

  //
  // Now that we have added all of the sources, we query the raytracer to get
  // all of the source properties into structs for the microphysics calls.
  // NOTE THAT IF THE NUMBER OF SOURCES OR THEIR PROPERTIES CHANGE OVER TIME,
  // I WILL HAVE TO WRITE NEW CODE TO UPDATE THIS!
  //
  FVI_nheat = RT->N_heating_sources();
  FVI_nion  = RT->N_ionising_sources();
  FVI_heating_srcs.resize(FVI_nheat);
  FVI_ionising_srcs.resize(FVI_nion);
  RT->populate_UVheating_src_list(FVI_heating_srcs);
  RT->populate_ionising_src_list( FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (RT->type_of_RT_integration()==RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT
            && SimPM.EP.MP_timestep_limit==5) {
    // For implicit updates to limit by xdot and/or edot
    // Here the raytracing has not already been done, so we call it here.
    FVI_need_column_densities_4dt = true;
  }
  else {
    FVI_need_column_densities_4dt = false;
  }

  cout <<"----------------- RAYTRACER SETUP COMPLETE -----------------------\n";
  return 0;
}

