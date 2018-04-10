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

#include "tools/reporting.h"

#include "tools/command_line_interface.h"

#include "sim_control.h"

#include "dataIO/dataio.h"
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

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_hydro_adi.h"
//#include "spatial_solvers/solver_eqn_hydro_adi_Eint.h"
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



void setup_fixed_grid::setup_cell_extra_data(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
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
      class SimParams &SimPM,  ///< pointer to simulation parameters
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
  else
    rep.error("Spatial order of accuracy unhandled by boundary conditions!",SimPM.spOOA);
  
  // Force Nbc=1 if using Lax-Friedrichs flux.
  if (SimPM.solverType==FLUX_LF)
  {SimPM.spOOA = SimPM.tmOOA = OA1; SimPM.Nbc=1;}

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  setup_cell_extra_data(SimPM);

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
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class GridBaseClass *grid 
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<SimPM.Nbc<<"\n";
#endif
  int err = grid->SetupBCs(SimPM);
  if (err) rep.error("setup_fixed_grid::boundary_conditions() Couldn't \
                      set up boundary conditions class.",err);
#ifdef TESTING
  cout <<"(setup_fixed_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_microphysics(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
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
    MP = new mp_only_cooling(SimPM.nvar, &(SimPM.EP), &(SimPM.RS));
    if (!MP) rep.error("mp_only_cooling() init",MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    string mptype;
    mptype = SimPM.chem_code;
    bool have_set_MP=false;
    cout <<"setting up MP type: "<<mptype<<"\n";

#ifdef LEGACY_CODE
    if      (mptype=="MPv0") {
      cout <<"\t******* setting up MPv0 module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv0(SimPM.nvar, SimPM.ntracer, SimPM.chem_code, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
        rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }

    if      (mptype=="MPv1") {
      cout <<"\t******* setting up MPv1 microphysics module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv1(SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      cout <<"\t**---** WARNING, THIS MODULE HAS BEEN SUPERSEDED BY MPv4. **--**\n";
      have_set_MP=true;
    }

    if      (mptype=="MPv2") {
      cout <<"\t******* setting up MPv2 module *********\n";
      cout <<"\t******* N.B. Timestep limiting is enforced. **\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv2(SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      SimPM.EP.MP_timestep_limit = 1;
      have_set_MP=true;
    }

    if (mptype=="MPv4") {
      cout <<"\t******* setting up MPv4 module *********\n";
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
      MP = new MPv4(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

    if (mptype=="MPv8") {
      cout <<"\t******* setting up MPv8 module *********\n";
      cout <<"\t******* This is for StarBench test propblems with heating and cooling.\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv8(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

#endif // LEGACY_CODE


#ifndef EXCLUDE_HD_MODULE
    if (mptype=="MPv9") {
      cout <<"\t******* setting up microphysics_lowz module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new microphysics_lowz(SimPM.nvar, SimPM.ntracer, SimPM.tracers, &(SimPM.EP), &(SimPM.RS));
      have_set_MP=true;
    }
#endif // exclude Harpreet's module



    if (mptype=="MPv3") {
      cout <<"\t******* setting up MPv3 module *********\n";
#if MPV3_DTLIMIT>=0 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. Timestep limiting is enforced by #def";
      cout <<" MPV3_DTLIMIT="<<MPV3_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif

      MP = new MPv3(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      //if (SimPM.EP.MP_timestep_limit != 1)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }

    if (mptype=="MPv5") {
      cout <<"\t******* setting up MPv5 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv5(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

    if (mptype=="MPv6") {
      cout <<"\t******* setting up MPv6 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv6(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

    if (mptype=="MPv7") {
      cout <<"\t******* setting up MPv7 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv7(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }


#ifdef CODE_EXT_HHE
    if (mptype=="MPv10") {
      cout <<"\t******* setting up MPv10 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv9_HHe(SimPM.nvar, SimPM.ntracer, SimPM.tracers, 
                        &(SimPM.EP), SimPM.gamma);
      have_set_MP=true;
    }
#endif


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
      class SimParams &SimPM,  ///< pointer to simulation parameters
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
    RT = new raytracer_USC_infinity(grid, MP, SimPM.ndim,
                            SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!RT) rep.error("init pllel-rays raytracer error",RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    RT = new raytracer_USC(grid, MP, SimPM.ndim, SimPM.coord_sys,
                          SimPM.nvar, SimPM.ftr, SimPM.RS.Nsources);
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


// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_evolving_RT_sources(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // Loop through list of sources, and see if any of them have an evolution
  // file (if none, then the string is set to NOFILE in the data I/O stage).
  //
  int err=0;
  int Nevo=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].EvoFile == "NOFILE") {
#ifdef TESTING
      cout <<"setup_evolving_RT_sources() Source "<<isrc<<" has no evolution file.\n";
#endif
    }
    else {
      //if (SimPM.RS.sources[isrc].effect != RT_EFFECT_PION_MULTI) {
      //  rep.error("setup_evolving_RT_sources() Source is not multifreq but has EvoFile",isrc);
      //}
      Nevo++;
      struct star istar;
#ifdef TESTING
      cout <<"setup_evolving_RT_sources() Source "<<isrc<<" has EvoFile "<<istar.file_name<<"\n";
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
  for (int isrc=0; isrc<Nevo; isrc++) {
    struct star *istar = &(SimPM.STAR[isrc]);
    istar->Nlines = 0;
    istar->time.resize(0);
    istar->Log_L.resize(0);
    istar->Log_T.resize(0);
    istar->Log_R.resize(0);
    istar->Log_V.resize(0);
    //
    // Open file
    //
    FILE *infile=0;
    infile = fopen(istar->file_name.c_str(), "r");
    if (!infile) rep.error("can't open wind evolving radiation source file",infile);
    // Skip first two lines
    char line[512];
    char *rval = 0;
    rval = fgets(line,512,infile);
    if (!rval) rep.error("setup_fixed_grid, RS, file read 1",line);
    //printf("%s",line);
    rval = fgets(line,512,infile);
    if (!rval) rep.error("setup_fixed_grid, RS, file read 2",line);
    //printf("%s",line);
    // Temporary variables for column values
    // Columns are time, M, L, Teff, Mdot, vrot
    double t1=0.0, t2=0.0, t3=0.0, t4=0.0, t5=0.0, t6=0.0;
    size_t iline=0;
    while (fscanf(infile, "   %lE   %lE %lE %lE %lE %lE", &t1, &t2, &t3, &t4, &t5, &t6) != EOF){
      //cout.precision(16);
      //cout <<t1 <<"  "<<t2  <<"  "<< t3  <<"  "<< t4 <<"  "<< t5 <<"  "<< t6 <<"\n";
      //cout <<t1 <<"  "<< t3  <<"  "<< t4 <<"  "<< t5 <<"  "<< t6 <<"\n";
      istar->Nlines ++;
      istar->time.push_back(t1);
      istar->Log_L.push_back( log10(t3/pconst.Lsun()) );
      istar->Log_T.push_back( log10(t4) );
      //
      // For ionisation rate, we need to rescale the Blackbody luminosity so
      // that it is much smaller for T<30000K, b/c the actual ionising photon
      // luminosity of these stars is much less than indicated by BB curve.
      // I took data from Table 1 of Diaz-Miller, Franco, & Shore,
      // (1998,ApJ,501,192), compared them to the ionising photon luminosity
      // of a BB with the same radius and Teff, and got the following scaling
      // factor, using file conversion.py in code_misc/testing/planck_fn/
      //
      if (istar->Log_T[iline]<4.53121387658 &&
          SimPM.RS.sources[isrc].effect == RT_EFFECT_PION_MULTI) {
        //cout <<"L(BB) ="<<exp(pconst.ln10()*istar->Log_L[i])<<", T=";
        //cout <<exp(pconst.ln10()*istar->Log_T[i])<<", scale-factor=";
        double beta = -4.65513741*istar->Log_T[iline] + 21.09342323;
        //istar->Log_L[iline] -= 2.0*beta;
        // HACK!!! THIS REDUCES FLUX BY LESS THAN IT SHOULD...
        istar->Log_L[iline] -= 0.4*beta;
        //cout <<", new L = "<<exp(pconst.ln10()*istar->Log_L[i])<<"\n";
      }

      istar->Log_V.push_back( log10(t6/1.0e5) );

      // Stellar radius, from Stefan Boltzmann Law.
      t6 = sqrt( pow(10.0,istar->Log_L[iline])*pconst.Lsun()/ 
                (4.0*pconst.pi()*pconst.StefanBoltzmannConst()*pow(t4, 4.0)));
      istar->Log_R.push_back( log10(t6/pconst.Rsun() ));

      //if (SimPM.RS.sources[isrc].effect == RT_EFFECT_PION_MULTI) {
      //  cout <<t1 <<"!!"<< pow(10.0,istar->Log_L[iline])*pconst.Lsun()  <<"  "<< pow(10.0,istar->Log_T[iline]) <<"  "<< t5 <<"  "<< pow(10.0,istar->Log_V[iline]) <<"\n";
      //}
      iline ++;
    }
    fclose(infile);

    //
    // Finally set the last_line counter to be the array index nearest to
    // (but less than) the current time.
    //
    iline=0;
    while (istar->time[iline] < SimPM.simtime) iline++;
    istar->last_line = iline;

    // initialise to zero.
    istar->Lnow = istar->Tnow = istar->Rnow = istar->Vnow = 0.0;
  }

  //
  // All done setting up the source.  Now we need to update the SimPM.RS.
  // source properties and send the changes to the raytracing class.  Need
  // time in secs, L,T,V in cgs and R in Rsun.
  //
  err = update_evolving_RT_sources(SimPM);
  return err;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::update_evolving_RT_sources(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  int err=0;
  bool updated=false;
  //
  // Loop over all sources with Evolution files.
  //
  for (unsigned int isrc=0; isrc<SimPM.STAR.size(); isrc++) {
    struct star *istar = &(SimPM.STAR[isrc]);
    istar->t_now = SimPM.simtime;
    size_t i = istar->last_line;

    //
    // Check if we have reached the last line of the file!
    //
    if (i==(istar->Nlines-1)) {
      cout <<"\n\n*#*#*#*#* WARNING #*#*#*#*#*#* update_evolving_RT_sources()";
      cout <<" Last line, assuming star is constant luminosity from now on!\n\n";
      return 0;
    }
    //
    // Check if we have moved forward one line in the code, in which
    // case we need to increment i.
    //
    while (istar->t_now > istar->time[i+1]) {
      //cout <<"update_evolving_RT_sources() Source has moved to next line. i="<<i<<" time="<<istar->time[i]<<"\n";
      i++;
      istar->last_line = i;
    }

    //
    // Now we know the star properties are bracketed by line i and line i+1,
    // so we can do a simple linear interpolation between them.
    //
    // First interpolate in log space.
    //
    double Lnow, Tnow, Rnow, Vnow;
    Lnow = istar->Log_L[i] +(istar->t_now-istar->time[i])*
        (istar->Log_L[i+1]-istar->Log_L[i])/(istar->time[i+1]-istar->time[i]);
    Tnow = istar->Log_T[i] +(istar->t_now-istar->time[i])*
        (istar->Log_T[i+1]-istar->Log_T[i])/(istar->time[i+1]-istar->time[i]);
    Rnow = istar->Log_R[i] +(istar->t_now-istar->time[i])*
        (istar->Log_R[i+1]-istar->Log_R[i])/(istar->time[i+1]-istar->time[i]);
    Vnow = istar->Log_V[i] +(istar->t_now-istar->time[i])*
        (istar->Log_V[i+1]-istar->Log_V[i])/(istar->time[i+1]-istar->time[i]);
    //
    // Now convert units (Radius is ok, but all others need conversion).
    //
    Lnow = exp(pconst.ln10()*(Lnow))*pconst.Lsun();   // erg/s
    Tnow = exp(pconst.ln10()*(Tnow));                 // K
    Rnow = exp(pconst.ln10()*(Rnow));                 // Rsun  (!!!)
    Vnow = exp(pconst.ln10()*(Vnow));                 // km/s  (!!!)

    //
    // If L or T change by more than 1% then update them; otherwise leave as they are.
    //
    if ( fabs(Lnow-istar->Lnow)/istar->Lnow >0.01 || fabs(Tnow-istar->Tnow)/istar->Tnow >0.01 ) {
      cout <<"update_evolving_RT_sources() NOW: t="<<istar->t_now;
      cout <<"\tL="<< Lnow;
      cout <<"\tT="<< Tnow;
      cout <<"\tR="<< Rnow;
      cout <<"\tV="<< Vnow <<"\n";
      istar->Lnow = Lnow;
      istar->Tnow = Tnow;
      istar->Rnow = Rnow;
      istar->Vnow = Vnow;

      //
      // Copy new data into SimPM.RS, and send updates to RayTracing and
      // microphysics classes.
      //
      struct rad_src_info *rs = &(SimPM.RS.sources[istar->src_id]);
      
      rs->strength = istar->Lnow;
      rs->Rstar    = istar->Rnow;
      rs->Tstar    = istar->Tnow;
      
      //
      // This is a horrible hack, fix it to something sensible ASAP!!!
      //
      if (rs->effect == RT_EFFECT_UV_HEATING) {
        rs->strength = 1.0e48*(rs->strength/1.989e38)*exp(-1e4/rs->Tstar);
        cout <<"FUV source: strength="<<rs->strength<<"\n";
      }

      RT->update_RT_source_properties(rs);

      if (rs->effect==RT_EFFECT_PION_MULTI) {
        err += MP->set_multifreq_source_properties(rs);
        if (err) rep.error("update_evolving_RT_sources() failed to update MP for source id",rs->id);
      }

      updated=true;
    }

  }
  

  //
  // Finally get the data back from RT into the structs for MP updates.
  //
  if (updated) {
    RT->populate_UVheating_src_list(FVI_heating_srcs);
    RT->populate_ionising_src_list( FVI_ionising_srcs);
  }
  
  return 0;
}





