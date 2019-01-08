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
#include "raytracing/raytracer_SC_pllel.h"
#include "grid/uniform_grid_pllel.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo_utility.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits_MPI.h"
#endif // if FITS

#include <iostream>
using namespace std;

#ifdef PARALLEL




// ##################################################################
// ##################################################################


setup_fixed_grid_pllel::setup_fixed_grid_pllel()
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
      vector<class GridBaseClass *> &g,  ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
#ifdef TESTING
  cout <<"setup_fixed_grid_pllel: setting up parallel grid.\n";
#endif
  class GridBaseClass **grid = &g[0];
  class MCMDcontrol *MCMD = &(SimPM.levels[0].MCMD);

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
  // Set Cell dx in cell interface class, and also xmin.
  //
  double dx = (SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX];
  CI.set_nlevels(dx,1);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
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
      SimPM.Xmin, SimPM.Xmax,
      SimPM.Xmin, SimPM.Xmax);
  }
  else if (SimPM.coord_sys==COORD_CYL) {
    *grid = new uniform_grid_cyl_parallel (
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
      MCMD->LocalXmin, MCMD->LocalXmax, MCMD->LocalNG,
      SimPM.Xmin, SimPM.Xmax,
      SimPM.Xmin, SimPM.Xmax);
  }
  else if (SimPM.coord_sys==COORD_SPH) {
    *grid = new uniform_grid_sph_parallel (
      SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
      MCMD->LocalXmin, MCMD->LocalXmax, MCMD->LocalNG,
      SimPM.Xmin, SimPM.Xmax,
      SimPM.Xmin, SimPM.Xmax);
  }
  else {
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
  }


  if (*grid==0)
    rep.error("(setup_fixed_grid_pllel::setup_grid) Couldn't assign data!", *grid);

#ifdef TESTING
  cout <<"(setup_fixed_grid_pllel::setup_grid) Done. ";
  cout <<"grid="<<*grid<<", and";
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
      class GridBaseClass *grid ///< pointer to grid
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
  grid->RT=0;
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
    grid->RT = new raytracer_USC_infinity(grid,MP, SimPM.ndim,
                            SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!grid->RT) rep.error("init pllel-rays raytracer error",grid->RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    grid->RT = new raytracer_USC_pllel(grid,MP,&SimPM,&(SimPM.levels[0].MCMD), SimPM.ndim, SimPM.coord_sys,
                          SimPM.nvar, SimPM.ftr, SimPM.RS.Nsources);
    if (!grid->RT) rep.error("init raytracer error 2",grid->RT);
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
  
  cout <<"***************** RAYTRACER SETUP ***********************\n";
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::boundary_conditions(
      class SimParams &par,     ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      //class GridBaseClass *grid ///< pointer to grid.
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<par.Nbc<<"\n";
#endif
  //
  // Choose what BCs to set up based on BC strings.
  //
  int err = setup_boundary_structs(par,grid[0],0);
  rep.errorTest("sfg::boundary_conditions::sb_structs",0,err);

  //
  // Ask grid to set up data for external boundaries.
  //
  err = grid[0]->SetupBCs(par);
  rep.errorTest("sfg::boundary_conditions::SetupBCs",0,err);

#ifdef TESTING
  cout <<"(setup_fixed_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::setup_boundary_structs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      const int l
      )
{
  string fname="setup_fixed_grid_pllel::setup_boundary_structs";
#ifdef TESTING
  cout <<"PLLEL: Set BC types...\n";
#endif 
  //
  // call serial version of setBCtypes, to set up the boundaries
  //
  int err = setup_fixed_grid::setup_boundary_structs(par,grid,l);
  if (err) {
    rep.error("sfg_pllel::setup_boundary_structs:: serial",err);
  }
  
  //
  // Now go through the 6 directions, and see if we need to replace
  // any edge boundaries with MPI communication boundaries, if the
  // local grid does not reach the global grid boundary.
  //
  int i=0;
  string temp;
  // Loop through boundaries, and if local boundary is not sim boundary,
  // set it to be a parallel boundary.
  for (i=0; i<par.ndim; i++) {
    if (!pconst.equalD(grid->Xmin(static_cast<axes>(i)), par.Xmin[i])) {
      // local xmin is not Sim xmin, so it's an mpi boundary
      grid->BC_bd[2*i]->itype=BCMPI;
      grid->BC_bd[2*i]->type="BCMPI";
    }
    if (!pconst.equalD(grid->Xmax(static_cast<axes>(i)), par.Xmax[i])) {
      // local xmax is not Sim xmin, so it's an mpi boundary
      grid->BC_bd[2*i+1]->itype=BCMPI;
      grid->BC_bd[2*i+1]->type="BCMPI";
    }
  }
  
#ifdef TESTING
  cout <<"PLLEL: BC types and data set up.\n";
#endif 

  //
  // If we have periodic boundaries, need to set neighbouring processors to
  // wrap around.  So set the number of procs in each direction.
  //
  int nx[par.ndim];
  class MCMDcontrol *ppar = &(par.levels[0].MCMD);
  for (i=0;i<par.ndim;i++) {
    nx[i] =static_cast<int>(ONE_PLUS_EPS*par.Range[i]/
                            grid->Range(static_cast<axes>(i)));
  }
  for (i=0; i<2*par.ndim; i++) {
    if (grid->BC_bd[i]->itype == PERIODIC) {
      switch (i) {
       case XN:
        ppar->ngbprocs[XN] = ppar->get_myrank() +nx[XX] -1;
        break;
       case XP:
        ppar->ngbprocs[XP] = ppar->get_myrank() -nx[XX] +1;
        break;
       case YN:
        ppar->ngbprocs[YN] = ppar->get_myrank() +(nx[YY]-1)*nx[XX];
        break;
       case YP:
        ppar->ngbprocs[YP] = ppar->get_myrank() -(nx[YY]-1)*nx[XX];
        break;
       case ZN:
        ppar->ngbprocs[ZN] = ppar->get_myrank() +
                              (nx[ZZ]-1)*nx[YY]*nx[XX];
        break;
       case ZP:
        ppar->ngbprocs[ZP] = ppar->get_myrank() -
                              (nx[ZZ]-1)*nx[YY]*nx[XX];
        break;
       default:
        rep.error("sfg_pllel::setup_boundary_structs: Bad direction",i);
        break;
      } // set neighbour according to direction.

      if ( (ppar->ngbprocs[i]<0) ||
           (ppar->ngbprocs[i]>=ppar->get_nproc()) )
        rep.error("sfg_pllel::setup_boundary_structs: Bad periodic \
                   neighbour",ppar->ngbprocs[i]);
      if (ppar->ngbprocs[i] == ppar->get_myrank()) {
        //  cout <<"setup_fixed_grid_pllel::setup_boundary_structs: only one proc in dir [i]: "<<i<<"\n";
        //  cout <<"setup_fixed_grid_pllel::setup_boundary_structs: periodic on single proc, so setting ngb to -999.\n";
        ppar->ngbprocs[i] = -999;
      }
    } // if periodic  
#ifdef TESTING
    cout<<"Neighbouring processor in dir "<<i<<" = "<<ppar->ngbprocs[i]<<"\n";
#endif // TESTING
  } // loop over directions.
  return(0);
}




// ##################################################################
// ##################################################################



void setup_fixed_grid_pllel::setup_dataio_class(
      class SimParams &par,     ///< simulation parameters
      const int typeOfFile ///< type of I/O: 1=text,2=fits,5=silo
      )
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

  case 1: // Start From ASCII Parameterfile.
    rep.error("No text file for parallel I/O!",typeOfFile);
    break;

#ifdef FITS
  case 2: // Start from FITS restartfile
    dataio = new DataIOFits_pllel
              (par, &(par.levels[0].MCMD));
    break;
#endif // if FITS

#ifdef SILO
  case 5: // Start from Silo ICfile or restart file.
    dataio = new dataio_silo_utility 
              (par, "DOUBLE", &(par.levels[0].MCMD));
    break; 
#endif // if SILO
  default:
    rep.error("sim_control::Init unhandled filetype",typeOfFile);
  }
  return;
}



// ##################################################################
// ##################################################################






#endif // PARALLEL


