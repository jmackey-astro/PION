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
#include "tools/mem_manage.h"

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
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"
#include "setup_fixed_grid.h"
#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <climits>
using namespace std;



// ##################################################################
// ##################################################################


setup_fixed_grid::setup_fixed_grid()
{
  FVI_nheat = FVI_nion = 0;
  FVI_heating_srcs.clear();
  FVI_ionising_srcs.clear();
  FVI_need_column_densities_4dt = false;
  eqn=0;
}



// ##################################################################
// ##################################################################


setup_fixed_grid::~setup_fixed_grid()
{
#ifdef TESTING
  cout << "(setup_fixed_grid::Destructor) Deleting Grid Class..." <<"\n";
#endif
  if (MP)     {delete MP; MP=0;}
  if (eqn)    {delete eqn; eqn=0;}
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
  // Set Cell dx in cell interface class, and also xmin.
  //
  CI.set_dx((SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX]);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);

  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(setup_fixed_grid::setup_grid) Setting up grid...\n";
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
  cout <<"(setup_fixed_grid::setup_grid) Done. &grid="<< grid<<", and grid="<<*grid<<"\n";
  cout <<"DX = "<<(*grid)->DX()<<"\n";
  dp.grid = (*grid);
#endif
  cout <<"------------------------------------------------------\n\n";

  return(0);
} // setup_grid()



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



    if (mptype=="MPv3" || mptype=="MPv3__") {
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

    if (mptype=="MPv5" || mptype=="MPv5__") {
      cout <<"\t******* setting up MPv5 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv5(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

    if (mptype=="MPv6" || mptype=="MPv6__") {
      cout <<"\t******* setting up MPv6 module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MPv6(SimPM.ndim, SimPM.coord_sys, SimPM.nvar,
                            SimPM.ntracer, SimPM.tracers,
                            &(SimPM.EP), &(SimPM.RS), SimPM.gamma);
      have_set_MP=true;
    }

    if (mptype=="MPv7" || mptype=="MPv7__") {
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
      class SimParams &SimPM,    ///< pointer to simulation parameters
      class GridBaseClass *grid ///< pointer to grid
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
    grid->RT = new raytracer_USC_infinity(grid, MP, SimPM.ndim,
                            SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!(grid->RT)) rep.error("init pllel-rays raytracer error",grid->RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    grid->RT = new raytracer_USC(grid, MP, SimPM.ndim, SimPM.coord_sys,
                          SimPM.nvar, SimPM.ftr, SimPM.RS.Nsources);
    if (!(grid->RT)) rep.error("init raytracer error 2",grid->RT);
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
      cout << grid->RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
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
      cout << grid->RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      uv_count++;
      dif_count++;
    } // if diffuse source
  } // loop over sources
  if (ion_count>1) {
    rep.error("Can only have one ionising source for currently implemented method",ion_count);
  }
  cout <<"Added "<<ion_count<<" ionising and "<<uv_count<<" non-ionising";
  cout <<" radiation sources, of which "<<dif_count<<" are diffuse radiation.\n";
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
  grid->RT->populate_ionising_src_list( FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (grid->RT->type_of_RT_integration()==RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (grid->RT && grid->RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT
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
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class RayTracingBase *RT   ///< pointer to raytracing class
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
  err = update_evolving_RT_sources(SimPM,RT);
  return err;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::update_evolving_RT_sources(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class RayTracingBase *RT   ///< pointer to raytracing class
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




// ##################################################################
// ##################################################################


int setup_fixed_grid::boundary_conditions(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid ///< pointer to grid.
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<par.Nbc<<"\n";
#endif
  //
  // Choose what BCs to set up based on BC strings.
  //
  int err = setup_boundary_structs(par,grid);
  rep.errorTest("sfg::boundary_conditions::sb_structs",0,err);

  //
  // Ask grid to set up data for external boundaries.
  //
  err = grid->SetupBCs(par);
  rep.errorTest("sfg::boundary_conditions::SetupBCs",0,err);

#ifdef TESTING
  cout <<"(setup_fixed_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::setup_boundary_structs(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid ///< pointer to grid.
      )
{
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif

  // Set number of boundaries: 2 for each dimension, plus internal.
  int len = 2*par.ndim + par.BC_Nint;
#ifdef TESTING
  cout <<"Got "<<len<<" boundaries to set up.\n";
#endif
  grid->BC_bd.resize(len);  // class member data
  BC_nbd = len;       // class member data
  for (int b=0;b<len;b++) grid->BC_bd[b] = mem.myalloc(grid->BC_bd[b],1);

  //
  // First the 2N external boundaries.  Put the strings into an array.
  //
  std::vector<std::string> bc_strings(2*par.ndim);
  bc_strings[0] = par.BC_XN; bc_strings[1] = par.BC_XP;
  if (par.ndim>1) {
    bc_strings[2] = par.BC_YN; bc_strings[3] = par.BC_YP;
  }
  if (par.ndim>2) {
    bc_strings[4] = par.BC_ZN; bc_strings[5] = par.BC_ZP;
  }

  //
  // Now go through each boundary and assign everything needed.
  //
  int i=0;
  for (i=0; i<2*par.ndim; i++) {
    grid->BC_bd[i]->dir = static_cast<direction>(i); //XN=0,XP=1,YN=2,YP=3,ZN=4,ZP=5
    grid->BC_bd[i]->ondir = grid->OppDir(grid->BC_bd[i]->dir);
#ifdef TESTING
    cout <<"i="<<i<<", dir = "<<grid->BC_bd[i]->dir<<", ondir="<< grid->BC_bd[i]->ondir<<"\n";
#endif
    grid->BC_bd[i]->baxis = static_cast<axes>(i/2);
    //
    // odd values of i are positive boundaries, others are negative.
    //
    if ((i+2)%2 !=0) {
      grid->BC_bd[i]->bloc  = grid->Xmax(grid->BC_bd[i]->baxis);
      grid->BC_bd[i]->bpos  = true;
    }
    else {
      grid->BC_bd[i]->bloc  = grid->Xmin(grid->BC_bd[i]->baxis);
      grid->BC_bd[i]->bpos  = false;
    }
    //
    // find boundary condition specified:
    //
    grid->BC_bd[i]->type = bc_strings[i];

    if      (grid->BC_bd[i]->type=="periodic") {
      grid->BC_bd[i]->itype=PERIODIC;
      grid->BC_bd[i]->type="PERIODIC";
    }
    else if (grid->BC_bd[i]->type=="outflow" || grid->BC_bd[i]->type=="zero-gradient") {
      grid->BC_bd[i]->itype=OUTFLOW;
      grid->BC_bd[i]->type="OUTFLOW";
    }
    else if (grid->BC_bd[i]->type=="one-way-outflow") {
      grid->BC_bd[i]->itype=ONEWAY_OUT;
      grid->BC_bd[i]->type="ONEWAY_OUT";
    }
    else if (grid->BC_bd[i]->type=="inflow") {
      grid->BC_bd[i]->itype=INFLOW ;
      grid->BC_bd[i]->type="INFLOW";
    }
    else if (grid->BC_bd[i]->type=="reflecting") {
      grid->BC_bd[i]->itype=REFLECTING;
      grid->BC_bd[i]->type="REFLECTING";
    }
    else if (grid->BC_bd[i]->type=="equator-reflect") {
      grid->BC_bd[i]->itype=JETREFLECT;
      grid->BC_bd[i]->type="JETREFLECT";
    }
    else if (grid->BC_bd[i]->type=="fixed") {
      grid->BC_bd[i]->itype=FIXED;
      grid->BC_bd[i]->type="FIXED";
    }
    else if (grid->BC_bd[i]->type=="DMR") {
      grid->BC_bd[i]->itype=DMACH;
      grid->BC_bd[i]->type="DMACH";
    }
    else {
      rep.error("Don't know this BC type",grid->BC_bd[i]->type);
    }

    if(!grid->BC_bd[i]->data.empty())
      rep.error("Boundary data not empty in constructor!",grid->BC_bd[i]->data.size());
    grid->BC_bd[i]->refval=0;
#ifdef TESTING
    cout <<"\tBoundary type "<<i<<" is "<<grid->BC_bd[i]->type<<"\n";
#endif
  }

  if (i<BC_nbd) {
#ifdef TESTING
    cout <<"Got "<<i<<" boundaries, but have "<<BC_nbd<<" boundaries.\n";
    cout <<"Must have extra BCs... checking for internal BCs\n";
#endif
    do {
      grid->BC_bd[i]->dir = NO;
      if (par.BC_Nint < i-2*par.ndim) {
        rep.error("Bad Number of boundaries",par.BC_Nint);
      }
      else {
        grid->BC_bd[i]->type = par.BC_INT[i-2*par.ndim];
      }

      if      (grid->BC_bd[i]->type=="jet") {
        grid->BC_bd[i]->itype=JETBC;
        grid->BC_bd[i]->type="JETBC";
      }
      else if (grid->BC_bd[i]->type=="DMR2") {
        grid->BC_bd[i]->itype=DMACH2;
        grid->BC_bd[i]->type="DMACH2";
      }
      else if (grid->BC_bd[i]->type=="stellar-wind") {
        grid->BC_bd[i]->itype=STWIND;
        grid->BC_bd[i]->type="STWIND";
      }
      else {
        rep.error("Don't know this BC type",grid->BC_bd[i]->type);
      }

      if(!grid->BC_bd[i]->data.empty()) {
        rep.error("Boundary data not empty in constructor!",
                  grid->BC_bd[i]->data.size());
      }
      grid->BC_bd[i]->refval=0;
#ifdef TESTING
      cout <<"\tBoundary type "<<i<<" is "<<grid->BC_bd[i]->type<<"\n";
#endif
      i++;
    } while (i<BC_nbd);
  }

  BC_nbd = grid->BC_bd.size();
#ifdef TESTING
  cout <<BC_nbd<<" BC structs set up.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::assign_boundary_data(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid  ///< pointer to grid.
      )
{
  // ----------------------------------------------------------------
  // only needed for nested grid I think... maybe also for stellar
  // wind.
  switch (par.eqntype) {
  case EQEUL:
    eqn = new eqns_Euler(par.nvar);
    break;
  case EQMHD:
    eqn = new eqns_mhd_ideal(par.nvar);
    break;
  case EQGLM:
    eqn = new eqns_mhd_mixedGLM(par.nvar);
    break;
  default:
    rep.error("Don't know the specified equations...",par.eqntype);
    break;
  }
  // ----------------------------------------------------------------
  int err=0;
  //
  // Loop through all boundaries, and assign data to them.
  //
  for (int i=0; i<BC_nbd; i++) {
    switch (grid->BC_bd[i]->itype) {
     case PERIODIC:   err += BC_assign_PERIODIC(  par,grid,grid->BC_bd[i]); break;
     case OUTFLOW:    err += BC_assign_OUTFLOW(   par,grid,grid->BC_bd[i]); break;
     case ONEWAY_OUT: err += BC_assign_ONEWAY_OUT(par,grid,grid->BC_bd[i]); break;
     case INFLOW:     err += BC_assign_INFLOW(    par,grid,grid->BC_bd[i]); break;
     case REFLECTING: err += BC_assign_REFLECTING(par,grid,grid->BC_bd[i]); break;
     case FIXED:      err += BC_assign_FIXED(     par,grid,grid->BC_bd[i]); break;
     case JETBC:      err += BC_assign_JETBC(     par,grid,grid->BC_bd[i]); break;
     case JETREFLECT: err += BC_assign_JETREFLECT(par,grid,grid->BC_bd[i]); break;
     case DMACH:      err += BC_assign_DMACH(     par,grid,grid->BC_bd[i]); break;
     case DMACH2:     err += BC_assign_DMACH2(    par,grid,grid->BC_bd[i]); break;
     case STWIND:     err += BC_assign_STWIND(    par,grid,grid->BC_bd[i]); break;
     case NEST_FINE: break; // assigned in nested grid class
     case NEST_COARSE: break; // assigned in nested grid class
     default:
      rep.warning("Unhandled BC",grid->BC_bd[i]->itype,-1); err+=1; break;
    }

  }
  return(err);
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_PERIODIC(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  enum direction ondir  = b->ondir;
  if (b->data.empty())
    rep.error("BC_assign_PERIODIC: empty boundary data",b->itype);

  // loop through all cells in boundary, and make the cell's "npt"
  // pointer point to the corresponding cell wrapped around on the
  // other side of the grid.  Assign data from that cell too.
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp; unsigned int ct=0;
  do{
    //
    // Go across the grid NG[baxis] cells, so that we map onto the
    // cell on the other side of the grid.
    //
    temp=(*bpt);
    for (int v=0; v<par.NG[b->baxis]; v++) {
      temp=grid->NextPt(temp,ondir);
    }
    for (int v=0;v<par.nvar;v++) (*bpt)->P[v]  = temp->P[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->Ph[v] = temp->P[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->dU[v] = 0.0;
    (*bpt)->npt = temp;
    //
    // increment counters.
    //
    ct++;
    ++bpt;
  } while (bpt !=b->data.end());

  if (ct != b->data.size())
    rep.error("BC_assign_PERIODIC: missed some cells!",ct-b->data.size());
  return 0;
}





// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_OUTFLOW(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // Zero order extrapolation, if edge cell is index k, and boundary cells
  // k+1 and k+2, then P_{k+1} = P_{k+2} = P_{k}
  // First order extrapolation would be P_{k+1} = 2P_{k} - P_{k-1} but is
  // said to have stability problems sometimes (LeVeque, S.7.2.1, p131-132)
  //
  enum direction ondir  = b->ondir;

  if (b->data.empty()) {
    rep.error("BC_assign_OUTFLOW: empty boundary data",b->itype);
  }
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp = 0;
  unsigned int ct = 0;  // counter (for accounting).

  //
  // loop through all boundary cells, set npt to point to the grid
  // cell where we get the data.
  //
  //cout <<"\t\t**** PRINTING CELLS FOR BOUNDARY DIR = "<<b->dir<<" ****\n\n";
  do{
    //
    // Find the cell to point to.  This should be the first cell on
    // the other side of the boundary, located at x[baxis]=bloc.
    // Can't just use the ->isgd property because in y and z dirs
    // the "ondir" direction will sometimes never hit the grid.
    //
    temp = (*bpt);
    //
    // If nbc==2, then we have two pointers to two sheets of cells,
    // and we want both of them to point to the first on-grid cell,
    // so we use "isedge" to move 1 or 2 cells on-grid to get the 
    // on-grid value.
    //
    for (int v=0; v>(*bpt)->isedge; v--) {
      //CI.print_cell(temp);
      //cout <<"ondir="<<ondir<<"\n";
      temp = grid->NextPt(temp,ondir);
    }

    for (int v=0;v<par.nvar;v++) (*bpt)->P[v]  = temp->P[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->Ph[v] = temp->P[v];
    (*bpt)->npt = temp;
    //CI.print_cell((*bpt));

    //
    // The GLM boundary is somewhat different, because I found
    // that zero-gradient didn't work well (Mackey & Lim, 2011,
    // MNRAS,412,2079).  So we switch the sign instead.
    //

#ifdef GLM_ZERO_BOUNDARY
    if (par.eqntype==EQGLM) {
      (*bpt)->P[SI]=(*bpt)->Ph[SI]=0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (par.eqntype==EQGLM) {
      if (par.ndim==1) {
        rep.error("Psi outflow boundary condition doesn't work for 1D!",99);
      }
      if ((*bpt)->isedge == -1) {
        (*bpt)->P[SI]  = -temp->P[SI];
        (*bpt)->Ph[SI] = -temp->Ph[SI];
      }
      else if ((*bpt)->isedge == -2) {
        //
        // Get data from 2nd on-grid cell.
        //
        (*bpt)->P[SI]  = -grid->NextPt(temp,ondir)->P[SI];
        (*bpt)->Ph[SI] = -grid->NextPt(temp,ondir)->Ph[SI];
      }
      else rep.error("only know 1st/2nd order bcs",(*bpt)->id);
    }
#endif // GLM_NEGATIVE_BOUNDARY
    ct++;
    ++bpt;
  } while (bpt !=b->data.end());

  if (ct != b->data.size()) {
    rep.error("BC_assign_OUTFLOW: missed some cells!",
              ct-b->data.size());
  }
  return 0;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_ONEWAY_OUT(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // The setup for this is identical to outflow, so just call
  // outflow() setup function.
  //
  int err=BC_assign_OUTFLOW(par,grid,b);
  return err;
}





// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_INFLOW(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  enum direction ondir  = b->ondir;

  if (b->data.empty()) {
    rep.error("BC_assign_INFLOW: empty boundary data",b->itype);
  }

  list<cell*>::iterator bpt=b->data.begin();
  cell *temp;
  unsigned int ct=0;
  do{
    //
    // Find the cell to point to.  This should be the first cell on
    // the other side of the boundary, located at x[baxis]=bloc.
    // Can't just use the ->isgd property because in y and z dirs
    // the "ondir" direction will sometimes never hit the grid.
    //
    temp = (*bpt);

    //
    // If nbc==2, then we have two pointers to two sheets of cells,
    // and we want both of them to point to the first on-grid cell,
    // so we use "isedge" to move 1 or 2 cells on-grid to get the 
    // on-grid value.
    //
    for (int v=0; v>(*bpt)->isedge; v--) {
      temp = grid->NextPt(temp,ondir);
    }

    //
    // Now set inflow data to be the first on-grid cell's values.
    //
    for (int v=0;v<par.nvar;v++) (*bpt)->P[v]  = temp->P[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->Ph[v] = temp->P[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->dU[v] = 0.0;
    ct++;
    ++bpt;
  } while (bpt !=b->data.end());

  if (ct != b->data.size())
    rep.error("BC_assign_INFLOW: missed some cells!",ct-b->data.size());
  
  return 0;
}





// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_REFLECTING(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  enum direction offdir = b->dir;
  enum direction ondir  = b->ondir;

  if (b->data.empty()) {
    rep.error("BC_assign_REFLECTING: empty boundary data",b->itype);
  }
  //
  // set reference state so that it is mostly zeros but has some +/-1
  // entries to flip the signs of the velocity and B-field (as 
  // appropriate).
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
    for (int v=0;v<par.nvar;v++)
      b->refval[v] = 1.0;
    //
    // velocity:
    //
    switch (offdir) {
     case XN: case XP:
      b->refval[VX] = -1.0;
      break;
     case  YN: case YP:
      b->refval[VY] = -1.0;
      break;
     case  ZN: case ZP:
      b->refval[VZ] = -1.0;
      break;
     default:
      rep.error("BAD DIRECTION REFLECTING",offdir);
      break;
    } // Set Normal velocity direction.
    
    //
    // B-field:
    //
    if (par.eqntype==EQMHD || par.eqntype==EQGLM || par.eqntype==EQFCD) {
      switch (offdir) {
       case XN: case XP:
        b->refval[BX] = -1.0;
        break;
       case  YN: case YP:
        b->refval[BY] = -1.0;
        break;
       case  ZN: case ZP:
        b->refval[BZ] = -1.0;
        break;
       default:
        rep.error("BAD DIRECTION REFLECTING",offdir);
        break;
      } // Set normal b-field direction.
    } // Setting up reference value.
  } // if we needed to set up refval.

  //
  // Now go through each of the boundary points and assign values
  // to them, multiplying the relevant entries by -1.
  //
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp=0;
  unsigned int ct=0;
  //cout <<"\t\t**** PRINTING CELLS FOR BOUNDARY DIR = "<<b->dir<<" ****\n\n";
  do{
    temp = (*bpt);
    for (int v=0; v>(*bpt)->isedge; v--) {
      temp = grid->NextPt(temp,ondir);
    }
    if(!temp) {
      rep.error("Got lost assigning reflecting bcs.",temp->id);
    }
    for (int v=0;v<par.nvar;v++)
      (*bpt)->P[v]  = temp->P[v]*b->refval[v];
    for (int v=0;v<par.nvar;v++)
      (*bpt)->Ph[v] = temp->Ph[v]*b->refval[v];
    for (int v=0;v<par.nvar;v++)
      (*bpt)->dU[v] = 0.0;
    (*bpt)->npt = temp;
    //CI.print_cell((*bpt));
    ++bpt;
    ct++;
  } while (bpt !=b->data.end());

  if (ct != b->data.size()) {
    rep.error("BC_assign_REFLECTING: missed some cells!",
              ct-b->data.size());
  }
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_FIXED(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  cout <<" setup_fixed_grid::BC_assign_FIXED starting\n";
  enum direction ondir  = b->ondir;
  if (b->data.empty()) {
    rep.error("BC_assign_FIXED: empty boundary data",b->itype);
  }
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  list<cell*>::iterator bpt=b->data.begin();
  cell *temp=0;
  unsigned int ct=0;
  //
  // First find an on-grid cell near a boundary point.  Because of
  // corner cells, we can't guarantee that every boundary cell will
  // reach an on-grid cell by moving in the on-grid direction.
  //
  cout <<"Finding first on-grid cell, size="<<b->data.size()<<".\n";
  do {
    //++bpt;
    temp = (*bpt);
    CI.print_cell(temp);
    for (int v=0; v>(*bpt)->isedge; v--) {
      temp = grid->NextPt(temp,ondir);
    }
  } while (!temp->isgd);
  if(!temp) rep.error("Got lost assigning FIXED bcs.",temp->id);

  //
  // Now set reference value to be the on-grid value.
  //
  cout <<"Setting reference value.\n";
  for (int v=0;v<par.nvar;v++) b->refval[v] = temp->P[v];
  //
  // Initialise all the values to be the fixed value.
  //
  bpt=b->data.begin();
  do{
    for (int v=0;v<par.nvar;v++) (*bpt)->P[v]  = b->refval[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->Ph[v] = b->refval[v];
    for (int v=0;v<par.nvar;v++) (*bpt)->dU[v] = 0.;
    ++bpt;
    ct++;
  } while (bpt !=b->data.end());

  if (ct != b->data.size()) {
    rep.error("BC_assign_FIXED: missed some cells!",
              ct-b->data.size());
  }
  
  return 0;
}




// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_JETBC(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  if (!JP.jetic) {
    rep.error("BC_assign_JETBC: not a jet simulation!",JP.jetic);
  }
  if (b->dir != NO) {
    rep.error("BC_assign_JETBC: boundary is not an internal one!",
              b->dir);
  }
  cell *c = grid->FirstPt();
  cell *temp=0, *cy=0;
  int ct=0;
  int ctot=0;
  int maxnv=0;

  //
  // Set the physical radius of jet.
  //
  double jr = JP.jetradius*grid->DX();
#ifdef TESTING
  cout <<"jetrad="<<JP.jetradius<<" dx="<<grid->DX()<<"\n";
#endif // TESTING

  //
  // Assign reference values, containing Jet parameters:
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  if (par.eqntype==EQEUL || par.eqntype==EQEUL_ISO ||
      par.eqntype==EQEUL_EINT || par.eqntype==EQMHD ||
      par.eqntype==EQGLM || par.eqntype==EQFCD) {
    rep.printVec("JetState",JP.jetstate,MAX_NVAR);
    b->refval[RO] = JP.jetstate[RO];
    b->refval[PG] = JP.jetstate[PG];
    b->refval[VX] = JP.jetstate[VX];
    b->refval[VY] = JP.jetstate[VY];
    b->refval[VZ] = JP.jetstate[VZ];
    maxnv=5;
  }
  else {
    rep.error("BC_assign_JETBC: bad equation type",par.eqntype);
  }

  if (par.eqntype==EQMHD || par.eqntype==EQGLM || par.eqntype==EQFCD) {
    b->refval[BX] = JP.jetstate[BX];
    b->refval[BY] = JP.jetstate[BY];
    b->refval[BZ] = JP.jetstate[BZ];
    maxnv=8;
  }

  if (par.eqntype==EQGLM) {
    b->refval[SI] = JP.jetstate[SI];
    maxnv=9;
  }

  for (int v=maxnv; v<par.nvar; v++) {
    b->refval[v] = JP.jetstate[v];
  }
  
  if (!b->data.empty()) {
    rep.error("BC_assign_JETBC: boundary data exists!",b->itype);
  }

  //
  // Simplest jet -- zero opening angle.  Check that this domain
  // is at the edge of the full simulation domain, if not then we
  // don't need to add any cells so skip this whole loop.
  //
  if (JP.jetic==1  && pconst.equalD(grid->Xmin(XX),par.Xmin[XX])) {
    //
    // Axi-symmetry first -- this is relatively easy to set up.
    //
    if (par.ndim==2 && par.coord_sys==COORD_CYL) {
      c = grid->FirstPt();
      do {
        temp=c;
        while ( (temp=grid->NextPt(temp,XN))!=0 ) {
          for (int v=0;v<par.nvar;v++) temp->P[v]  = b->refval[v];
          for (int v=0;v<par.nvar;v++) temp->Ph[v] = b->refval[v];
#ifdef SOFTJET
          //
          // jetradius is in number of cells, jr is in physical units.
          // Jet centre is along Z_cyl axis, centred on origin.
          //
          temp->P[VX]  = b->refval[VX]
                          *min(1., 4.0-4.0*CI.get_dpos(temp,YY)/jr);
          temp->Ph[VX] = temp->P[VX];
#endif //SOFTJET
          b->data.push_back(temp);
          ctot++;
          if (temp->isgd){
            rep.error("Looking for Boundary cells! setupjet",temp);
          }
        }
        ct++;
      } while ( (c=grid->NextPt(c,YP)) && ct<JP.jetradius);
      if (ct!=JP.jetradius) {
        rep.error("Not enough cells for jet",ct);
      }
      cout <<"Got "<<ctot<<" Cells in total for jet boundary.\n";
    } // 2D Axial Symmetry

    //
    // 3D now, more difficult.
    //
    else if (par.ndim==3 && par.coord_sys==COORD_CRT) {
      double dist=0.0;
      c = grid->FirstPt();
      //
      // 3D, so we need to convert the jet radius to a real length.
      // Also, the jet will come in at the centre of the XN boundary,
      // which must be the origin.
      //
      do { // loop over ZZ axis
        cy=c;
        do { // loop over YY axis
          dist = sqrt(CI.get_dpos(cy,YY)*CI.get_dpos(cy,YY) +
                      CI.get_dpos(cy,ZZ)*CI.get_dpos(cy,ZZ)  );
          //
          // if dist <= jr, then we are within the jet inflow, and we
          // add the cells to the boundary.
          //
          if (dist <= jr) {
            temp = cy;
            while ( (temp=grid->NextPt(temp,XN))!=0 ) {
              for (int v=0;v<par.nvar;v++) temp->P[v]  = b->refval[v];
              for (int v=0;v<par.nvar;v++) temp->Ph[v] = b->refval[v];
# ifdef SOFTJET
              temp->P[VX]  = b->refval[VX] *min(1., 4.0-4.0*dist/jr);
              temp->Ph[VX] = temp->P[VX];
# endif //SOFTJET
              b->data.push_back(temp);
              ctot++;
              if (temp->isgd) {
                rep.error("Looking for Boundary cells! setupjet",
                          temp);
              }
            }
          } // if within jet radius
        } while ( (cy=grid->NextPt(cy,YP)) ); // through cells on YY axis.
      } while ( (c=grid->NextPt(c,ZP)) );     // through cells on ZZ axis.
      //      BC_printBCdata(b);
    } // 3D Cartesian

    else {
      rep.error("Only know how to set up jet in 2Dcyl or 3Dcart",
                par.ndim);
    }

  } // jetic==1
  else {
    rep.error("Only know simple jet with jetic=1",JP.jetic);
  }

  return 0;
}


// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_JETREFLECT(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  enum direction offdir = b->dir;
  enum direction ondir  = b->ondir;
  if (b->data.empty()) {
    rep.error("BC_assign_: empty boundary data",b->itype);
  }

  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  for (int v=0;v<par.nvar;v++) b->refval[v] = 1.0;
  //
  // Set Normal velocity multiplier to -1 for reflection.
  //
  switch (offdir) {
   case XN: case XP:
    b->refval[VX] = -1.0;
    break;
   case  YN: case YP:
    b->refval[VY] = -1.0;
    break;
   case  ZN: case ZP:
    b->refval[VZ] = -1.0;
    break;
   default:
    rep.error("BAD DIRECTION",offdir);
  } // Set Normal velocity direction.
  //
  // Set normal B-field multiplier to -1 for reflection.
  //
  if (par.eqntype==EQMHD || par.eqntype==EQGLM || par.eqntype==EQFCD) {
    switch (offdir) {
     case XN: case XP:
      b->refval[BY] = b->refval[BZ] = -1.0;
      break;
     case  YN: case YP:
      b->refval[BX] = b->refval[BZ] = -1.0;
      break;
     case  ZN: case ZP:
      b->refval[BX] = b->refval[BY] = -1.0;
      break;
     default:
      rep.error("BAD DIRECTION",offdir);
    } // Set normal b-field direction.
  } // if B-field exists

  //
  // Now go through each column of boundary points and assign values
  // to them.
  //
  list<cell*>::iterator bpt=b->data.begin();
  cell *temp=0;
  unsigned int ct=0;

  do{
    temp = (*bpt);
    for (int v=0; v>(*bpt)->isedge; v--) {
      temp = grid->NextPt(temp,ondir);
    }
    if(!temp) {
      rep.error("Got lost assigning jet-reflecting bcs.",temp->id);
    }
    for (int v=0;v<par.nvar;v++)
      (*bpt)->P[v]  = temp->P[v]*b->refval[v];
    for (int v=0;v<par.nvar;v++)
      (*bpt)->Ph[v] = temp->P[v]*b->refval[v];
    for (int v=0;v<par.nvar;v++)
      (*bpt)->dU[v] = 0.0;
    (*bpt)->npt = temp;
    ++bpt;
    ct++;
  } while (bpt !=b->data.end());
  
  if (ct != b->data.size()) {
    rep.error("BC_assign_: missed some cells!",ct-b->data.size());
  }
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_DMACH(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
#ifdef TESTING
  cout <<"Setting up DMACH boundary... starting.\n";
#endif // TESTING

  if (b->data.empty()) {
    rep.error("BC_assign_DMACH: empty boundary data",b->itype);
  }
  //
  // Set reference value to be downstream values.
  //
  if (!b->refval) {
    b->refval = mem.myalloc(b->refval, par.nvar);
  }

  b->refval[RO] = 1.4;
  b->refval[PG] = 1.0;
  b->refval[VX] = 0.0;
  b->refval[VY] = 0.0;
  b->refval[VZ] = 0.0;
  for (int v=par.ftr; v<par.nvar; v++) b->refval[v] = -1.0;
  cout <<"par.ftr="<<par.ftr<<", par.nvar="<<par.nvar<<"\n";

  //
  // Run through all boundary cells, and give them either upstream or
  // downstream value, depending on their position.
  //
  list<cell*>::iterator bpt=b->data.begin();
  unsigned int ct=0;
  double bpos = 0.0;
  do {
    //
    // This is the boundary position:
    //
    bpos =  10.0*par.simtime/sin(M_PI/3.0)
          + 1.0/6.0
          + CI.get_dpos(*bpt,YY)/tan(M_PI/3.0);

    if (CI.get_dpos(*bpt,XX) <= bpos) {
      (*bpt)->P[RO] = 8.0;
      (*bpt)->P[PG] = 116.5;
      (*bpt)->P[VX] = 7.14470958;
      (*bpt)->P[VY] = -4.125;
      (*bpt)->P[VZ] = 0.0;
      for (int v=par.ftr; v<par.nvar; v++) (*bpt)->P[v]  = 1.0;
      (*bpt)->Ph[RO] = 8.0;
      (*bpt)->Ph[PG] = 116.5;
      (*bpt)->Ph[VX] = 7.14470958;
      (*bpt)->Ph[VY] = -4.125;
      (*bpt)->Ph[VZ] = 0.0;
      for (int v=par.ftr; v<par.nvar; v++) (*bpt)->Ph[v] = 1.0;
    }
    else {
      for (int v=0;v<par.nvar;v++) (*bpt)->P[v]  = b->refval[v];
      for (int v=0;v<par.nvar;v++) (*bpt)->Ph[v] = b->refval[v];
    }
    ++bpt;
    ct++;
  } while (bpt !=b->data.end());

  if (ct != b->data.size()) {
    rep.error("BC_assign_: missed some cells!",ct-b->data.size());
  }
#ifdef TESTING
  cout <<"Setting up DMACH boundary... finished.\n";
#endif // TESTING
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_DMACH2( 
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
#ifdef TESTING
  cout <<"Setting up DMACH2 boundary... starting.\n";
#endif // TESTING

  if (b->dir != NO) {
    rep.error("DMACH2 not internal boundary!",b->dir);
  }
  cout <<"DMACH2 boundary, from x=0 to x=1/6 at y=0, fixed bd.\n";
  if (b->refval) {
    rep.error("Already initialised memory in DMACH2 boundary refval",
              b->refval);
  }
  b->refval = mem.myalloc(b->refval, par.nvar);

  b->refval[RO] = 8.0;
  b->refval[PG] = 116.5;
  b->refval[VX] = 7.14470958;
  b->refval[VY] = -4.125;
  b->refval[VZ] = 0.0;
  for (int v=par.ftr; v<par.nvar; v++) b->refval[v] = 1.0;

  //
  // Now have to go from first point onto boundary and across to
  // x<=1/6
  //
  if (!b->data.empty()) {
    rep.error("BC_assign_DMACH2: Not empty boundary data",b->itype);
  }
  cell *c = grid->FirstPt();
  cell *temp=0;
  do {
    //
    // check if we are <1/6 in case we are in a parallel grid where
    // the domain is outside the DMR region.  This saves having to
    // rewrite the function in the parallel uniform grid.
    //
    if (CI.get_dpos(c,XX)<=1./6.) {
      temp=c;
      while ( (temp=grid->NextPt(temp,YN)) !=0 ) {
        for (int v=0;v<par.nvar;v++) temp->P[v]  = b->refval[v];
        for (int v=0;v<par.nvar;v++) temp->Ph[v] = b->refval[v];
        b->data.push_back(temp);
        if (temp->isgd) {
          rep.error("BC_assign_DMACH2: Looking for Boundary cells!",
                    temp);
        }
      }
    }
  } while ( (c=grid->NextPt(c,XP)) && (CI.get_dpos(c,XX)<=1./6.));

#ifdef TESTING
  cout <<"Setting up DMACH2 boundary... finished.\n";
#endif // TESTING
  return 0;
}



//
// Add internal stellar wind boundaries -- these are (possibly
// time-varying) winds defined by a mass-loss-rate and a terminal
// velocity.  A region within the domain is given fixed values
// corresponding to a freely expanding wind from a
// cell-vertex-located source.
//
int setup_fixed_grid::BC_assign_STWIND(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // Check that we have an internal boundary struct, and that we have
  // a stellar wind source to set up.
  //
  if (b->dir != NO)
    rep.error("STWIND not external boundary!",b->dir);
#ifdef TESTING
  cout <<"Assigning data to STWIND boundary. Nsrc=";
  cout <<SWP.Nsources<<"\n";
#endif
  if (SWP.Nsources<1) {
    rep.error("BC_assign_STWIND() No Sources!",SWP.Nsources);
  }

  //
  // Setup reference state vector and initialise to zero.
  //
  if (b->refval) {
    rep.error("Initialised STWIND boundary refval",b->refval);
  }
  b->refval = mem.myalloc(b->refval, par.nvar);
  if (!b->data.empty()) {
    rep.error("BC_assign_STWIND: Not empty boundary data",b->itype);
  }
  for (int v=0;v<par.nvar;v++) b->refval[v]=0.0;

  //
  // New structure: we need to initialise the stellar wind class with
  // all of the wind sources in the global parameter list (this was
  // formerly done in DataIOBase::read_simulation_parameters()).
  //
  // The type of class we set up is determined first.
  // stellar_wind_evolution is derived from stellar_wind, and 
  // stellar_wind_angle is derived from stellar_wind_evolution.
  //
  int err=0;
  int Ns = SWP.Nsources;
  for (int isw=0; isw<Ns; isw++) {
    if (SWP.params[isw]->type ==WINDTYPE_EVOLVING) err=1;
  }
  for (int isw=0; isw<Ns; isw++) {
    if (SWP.params[isw]->type ==WINDTYPE_ANGLE) err=2;
  }
  if (Ns>0) {
    cout <<"\n----------- SETTING UP STELLAR WIND CLASS ----------\n";
    if      (err==0) {
      grid->Wind = new stellar_wind(par.ndim, par.nvar, par.ntracer, par.ftr,
                              par.coord_sys, par.eqntype,
                              par.EP.MinTemperature);
    }
    else if (err==1) {
      grid->Wind = new stellar_wind_evolution(par.ndim, par.nvar, par.ntracer,
            par.ftr, par.coord_sys, par.eqntype, par.EP.MinTemperature,
            par.starttime, par.finishtime);
      err=0;
    }
    else if (err==2) {
      cout <<"Setting up stellar_wind_angle class\n";
      grid->Wind = new stellar_wind_angle(par.ndim, par.nvar, par.ntracer, par.ftr,
                  par.coord_sys, par.eqntype, par.EP.MinTemperature,
                  par.starttime, par.finishtime);
    }
  }

  //
  // Run through sources and add sources.
  //
  for (int isw=0; isw<Ns; isw++) {
    cout <<"\tUniGrid::BC_assign_STWIND: Adding source "<<isw<<"\n";
    if (SWP.params[isw]->type==WINDTYPE_CONSTANT) {
      //
      // This is for spherically symmetric winds that are constant
      // in time.
      //
      err = grid->Wind->add_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->Mdot,
        SWP.params[isw]->Vinf,
        SWP.params[isw]->Tstar,
        SWP.params[isw]->Rstar,
        SWP.params[isw]->tr
        );
    }
    else  {
      //
      // This works for spherically symmetric winds and for
      // latitude-dependent winds that evolve over time.
      //
      cout <<"Adding source "<<isw<<" with filename "<<SWP.params[isw]->evolving_wind_file<<"\n";
      err = grid->Wind->add_evolving_source(
        SWP.params[isw]->dpos,
        SWP.params[isw]->radius,
        SWP.params[isw]->type,
        SWP.params[isw]->Rstar,
        SWP.params[isw]->tr,
        SWP.params[isw]->evolving_wind_file,
        SWP.params[isw]->enhance_mdot,
        SWP.params[isw]->time_offset,
        par.simtime,
        SWP.params[isw]->update_freq,
        SWP.params[isw]->t_scalefactor
        );
    }
    if (err) rep.error("Error adding wind source",isw);
  }

  //
  // loop over sources, adding cells to boundary data list in order.
  //
  for (int id=0;id<Ns;id++) {
    cout <<"\tUniGrid::BC_assign_STWIND: Adding cells to source ";
    cout <<id<<"\n";
    BC_assign_STWIND_add_cells2src(par,grid, id);
  }
  //
  // Now we should have set everything up, so we assign the boundary
  // cells with their boundary values.
  //
  //err += BC_update_STWIND(par.simtime, b,0,0);
  cout <<"\tFinished setting up wind parameters\n";
  cout <<"------ DONE SETTING UP STELLAR WIND CLASS ----------\n\n";
  return err;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid::BC_assign_STWIND_add_cells2src(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      //struct boundary_data *b,
      const int id ///< source id
      )
{
  //
  // We run through each cell, and if it is within the 
  // source's radius of influence, then we add it to the lists.
  //
  int err=0;
  int ncell=0;
  double srcpos[MAX_DIM];
  double srcrad;
  grid->Wind->get_src_posn(id,srcpos);
  grid->Wind->get_src_drad(id,&srcrad);

#ifdef TESTING
  cout <<"*** srcrad="<<srcrad<<"\n";
  rep.printVec("src", srcpos, par.ndim);
#endif

  cell *c = grid->FirstPt_All();
  do {
    //
    // GEOMETRY: This is to the centre-of-volume of cell.
    // It makes no difference for Cartesian grids b/c the centre-of-
    // volume coincides the midpoint, but for non-cartesian geometry
    // it is different.
    //
#ifdef TESTING
    cout <<"cell: "<<grid->distance_vertex2cell(srcpos,c)<<"\n";
#endif
    if (grid->distance_vertex2cell(srcpos,c) <= srcrad) {
      ncell++;
      err += grid->Wind->add_cell(grid, id,c);
    }
  } while ((c=grid->NextPt_All(c))!=0);
  
  err += grid->Wind->set_num_cells(id,ncell);

#ifdef TESTING
  cout <<"setup_fixed_grid: Added "<<ncell;
  cout <<" cells to wind boundary for WS "<<id<<"\n";
#endif
  return err;
}



// ##################################################################
// ##################################################################



