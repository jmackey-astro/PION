/// \file setup_fixed_grid_MPI.cpp
/// 
/// \brief MPI-Parallel Class for setting up fixed grids.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for 
/// the "setup_fixed_grid_pllel" class, which is for setting up grids
/// 
/// Modifications:
/// - 2015.02.18 JM: new file for setting up parallel grids.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"

#include "decomposition/MCMD_control.h"
#include "setup_fixed_grid.h"
#include "setup_fixed_grid_MPI.h"
#include "raytracing/raytracer_SC.h"


#include <iostream>
using namespace std;

#ifdef PARALLEL




// ##################################################################
// ##################################################################


setup_fixed_grid_pllel::setup_fixed_grid_pllel()
  : setup_fixed_grid()
{
#ifdef TESTING
  cout <<"setup_fixed_grid_pllel constructor.\n";
#endif
}



// ##################################################################
// ##################################################################


setup_fixed_grid_pllel::~setup_fixed_grid_pllel()
{
#ifdef TESTING    
  cout <<"setup_fixed_grid_pllel destructor.\n";
#endif
}



// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::setup_grid(
      class GridBaseClass **grid, ///< address of pointer to computational grid.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class MCMDcontrol *MCMD     ///< address of MCMD controller class.
      )
{
#ifdef TESTING
  cout <<"setup_fixed_grid_pllel: setting up parallel grid.\n";
#endif

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
  // May need to setup extra data in each cell for ray-tracing
  // optical depths and/or viscosity variables.  Cells cannot be
  // created unless this the number of such extra variables has been
  // set.
  //
  setup_cell_extra_data(SimPM);

  //
  // Now set up the parallel uniform grid.
  //
#ifdef TESTING
  cout <<"(setup_fixed_grid_pllel::setup_grid) Setting up grid...\n";
#endif

  if      (SimPM.coord_sys==COORD_CRT) {
    *grid = new UniformGridParallel (
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
      MCMD->LocalXmin, MCMD->LocalXmax, MCMD->LocalNG,
      SimPM.Xmin, SimPM.Xmax, MCMD);
  }
  else if (SimPM.coord_sys==COORD_CYL) {
    *grid = new uniform_grid_cyl_parallel (
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
      MCMD->LocalXmin, MCMD->LocalXmax, MCMD->LocalNG,
      SimPM.Xmin, SimPM.Xmax, MCMD);
  }
  else if (SimPM.coord_sys==COORD_SPH) {
    *grid = new uniform_grid_sph_parallel (
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
      MCMD->LocalXmin, MCMD->LocalXmax, MCMD->LocalNG,
      SimPM.Xmin, SimPM.Xmax, MCMD);
  }
  else {
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
  }


  if (*grid==0)
    rep.error("(setup_fixed_grid_pllel::setup_grid) Couldn't assign data!", *grid);

#ifdef TESTING
  cout <<"(setup_fixed_grid_pllel::setup_grid) Done. ";
  cout <<"&grid="<<grid<<", and grid="<<*grid<<", and";//<<"\n";
  cout <<"\t DX = "<<(*grid)->DX()<<"\n";
  dp.grid = (*grid);
#endif
  cout <<"DX = "<<(*grid)->DX()<<"\n";

  return(0);
}



// ##################################################################
// ##################################################################


int setup_fixed_grid_pllel::setup_raytracing(
      class SimParams &SimPM,  ///< pointer to simulation parameters
        class GridBaseClass *grid
        )
{
  //
  // This function is identical to the serial setup function, except
  // that it sets up parallelised versions of the raytracers.
  //
  // If not doing raytracing, return immediately.
  //
  if (!SimPM.EP.raytracing) {
    return 0;
  }

  //
  // Now we are doing raytracing, so set up a raytracer and add
  // sources to it.
  //
  if (!MP) rep.error("can't do raytracing without microphysics",MP);
  cout <<"\n***************** RAYTRACER SETUP STARTING ***********************\n";
  RT=0;
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays=true;
  int dir=-1;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays=false;
    //
    // source is at infinity, so make sure all sources at infinity have rays
    // travelling in the same direction (by checking direction to source).
    //
    else {
      for (int i=0;i<SimPM.ndim;i++) {
        if (fabs(SimPM.RS.sources[isrc].pos[i])>1.e99) {
          if (dir==-1) dir=i;
          else if (dir!=i) parallel_rays=false;
        }
      }
    }
  }     // loop over sources.
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  parallel_rays=false;
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    RT = new raytracer_USC_infinity(grid,MP, SimPM.ndim,
                            SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!RT) rep.error("init pllel-rays raytracer error",RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    RT = new raytracer_USC_pllel(grid,MP, SimPM.ndim, SimPM.coord_sys,
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
  
  cout <<"***************** RAYTRACER SETUP ***********************\n";
  return 0;
}



// ##################################################################
// ##################################################################




#endif // PARALLEL


