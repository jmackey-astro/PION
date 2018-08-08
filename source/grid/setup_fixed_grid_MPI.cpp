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
      class GridBaseClass *grid,
      class RayTracingBase *RT ///< pointer to raytracing class
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



int setup_fixed_grid_pllel::setup_boundary_structs(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid ///< pointer to grid.
      class MCMDcontrol &ppar   ///< domain decomposition info
      )
{
  string fname="setup_fixed_grid_pllel::setup_boundary_structs";
#ifdef TESTING
  cout <<"PLLEL: Set BC types...\n";
#endif 
  //
  // call serial version of setBCtypes, to set up the boundaries
  //
  int err = setup_fixed_grid::setup_boundary_structs(par,grid);
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
  for (i=0; i<SimPM.ndim; i++) {
    if (!pconst.equalD(G_xmin[i], Sim_xmin[i])) {
      // local xmin is not Sim xmin, so it's an mpi boundary
      BC_bd[2*i].itype=BCMPI;
      BC_bd[2*i].type="BCMPI";
    }
    if (!pconst.equalD(G_xmax[i], Sim_xmax[i])) {
      // local xmax is not Sim xmin, so it's an mpi boundary
      BC_bd[2*i+1].itype=BCMPI;
      BC_bd[2*i+1].type="BCMPI";
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
  for (i=0;i<SimPM.ndim;i++) {
    nx[i] =static_cast<int>(ONE_PLUS_EPS*Sim_range[i]/G_range[i]);
  }
  for (i=0; i<2*SimPM.ndim; i++) {
    if (BC_bd[i].itype == PERIODIC) {
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
        rep.error("setup_fixed_grid_pllel::setup_boundary_structs: Bad direction",i);
        break;
      } // set neighbour according to direction.

      if ( (ppar->ngbprocs[i]<0) ||
           (ppar->ngbprocs[i]>=ppar->get_nproc()) )
        rep.error("setup_fixed_grid_pllel::setup_boundary_structs: Bad periodic \
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



int setup_fixed_grid_pllel::assign_boundary_data(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      class MCMDcontrol &ppar   ///< domain decomposition info
      )
{
#ifdef TESTING
  cout<<"setup_fixed_grid_pllel::assign_boundary_data()\n";
#endif // TESTING
  //
  // Loop through all boundaries, and assign data to them.
  // This is only different to the serial version in that it includes
  // BCMPI boundaries (MPI communication of internal boundaries).
  //
  int err=0;
  for (int i=0; i<BC_nbd; i++) {
    switch (BC_bd[i].itype) {
    case PERIODIC:
      err += BC_assign_PERIODIC(  par,ppar,grid,grid->BC_bd[i]);
      break;
    case OUTFLOW:
      err += BC_assign_OUTFLOW(   par,grid,grid->BC_bd[i]);
      break;
    case ONEWAY_OUT:
      err += BC_assign_ONEWAY_OUT(par,grid,grid->BC_bd[i]);
      break;
    case INFLOW:
      err += BC_assign_INFLOW(    par,grid,grid->BC_bd[i]);
      break;
    case REFLECTING:
      err += BC_assign_REFLECTING(par,grid,grid->BC_bd[i]);
      break;
    case FIXED:
      err += BC_assign_FIXED(     par,grid,grid->BC_bd[i]);
      break;
    case JETBC:
      err += BC_assign_JETBC(     par,grid,grid->BC_bd[i]);
      break;
    case JETREFLECT:
      err += BC_assign_JETREFLECT(par,grid,grid->BC_bd[i]);
      break;
    case DMACH:
      err += BC_assign_DMACH(     par,grid,grid->BC_bd[i]);
      break;
    case DMACH2:
      err += BC_assign_DMACH2(    par,grid,grid->BC_bd[i]);
      break;
    case BCMPI:
      err += BC_assign_BCMPI(par,ppar,grid,grid->BC_bd[i],BC_MPItag);
      break;
    case STWIND:
      err += BC_assign_STWIND(    par,grid,grid->BC_bd[i]);
      break;
     case FINE_TO_COARSE:
     case COARSE_TO_FINE:
      break; // assigned in nested grid class
     default:
      rep.error("Unhandled BC",grid->BC_bd[i]->itype);
      break;
    }
    if (i==XP || (i==YP && SimPM.ndim>1) || (i==ZP && SimPM.ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      COMM->barrier("setup_fixed_grid_pllel__SetupBCs");
    }
  }
#ifdef TESTING
  cout<<"setup_fixed_grid_pllel::assign_boundary_data() finished\n";
#endif // TESTING
  return err;
}



// ##################################################################
// ##################################################################




// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::BC_assign_PERIODIC(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // For parallel grid, periodic data may be on a different proc.,
  // which is already pointed to by ppar->ngbprocs[b->dir]
  // So I just have to call BC_assign_BCMPI and it will do the job.
  int err=0;
  if (ppar->ngbprocs[b->dir] <0) {
    // cout <<"BC_assign_PERIODIC: non comm periodic in direction "<<b->dir<<"\n";
    err = setup_fixed_grid::BC_assign_PERIODIC(par,grid,b);
  }
  else {
    // cout<<"BC_assign_PERIODIC: communicating periodic bc in direction "<<b->dir<<"\n";
    // cout<<"BC_assign_PERIODIC: calling mpi assign BC function\n";
    err = BC_assign_BCMPI(par,ppar,grid,b,BC_PERtag);
  }
  return err;
}



// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::BC_assign_BCMPI(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,
      int comm_tag
      )
{
  //
  // first choose cells to send to the appropriate processor.  This
  // is stored in a separate list "send_data" that is part of the
  // boundary_data struct.
  //
  int err = 0;
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: sending data in dir: "<<b->dir<<"\n";
#endif
  int ncell =0;
  if (b->send_data.size() != 0) {
    rep.error("send_data is not empty!",b->send_data.size());
  }
  err += BC_select_data2send(&(b->send_data), &ncell, b);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: got "<<ncell<<" cells in send_data\n";
#endif

  //
  // This is the same as the update function, except that we want
  // to set P[] and Ph[] vectors to the same values, so we set cstep
  // equal to maxstep.
  //
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: starting\n";
#endif 
  err = BC_update_BCMPI(par,ppar,grid,b,2,2,comm_tag);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: finished\n";
  cout <<"*******************************************\n";
#endif 
  return err;
}






// ##################################################################
// ##################################################################



int setup_fixed_grid_pllel::BC_select_data2send(
      list<cell *> *l,
      int *nc,
      boundary_data *b
      )
{
  // Check inputs.
  if ( !(*l).empty() ) {
#ifdef TESTING
    rep.warning("BC_select_data2send: List not empty! Emptying it now.",0,(*l).size());
#endif
    (*l).clear(); if(!(*l).empty()) rep.error("emptying list.",(*l).empty());
  }
  if ( *nc !=0 ) {
#ifdef TESTING
    rep.warning("BC_select_data2send: uninitialized counter in input. setting it to zero.",0,*nc);
#endif
    *nc=0;
  }
  
  //
  // We want to get all cells on the grid that are adjacent to the
  // boundary data, so we just go BC_nbc cells in the "ondir"
  // direction from every boundary cell, and this is a "send" cell.
  //
  int count=0;
  list<cell*>::iterator c=b->data.begin();
  cell *temp =0;
  do {
    temp = *c;
    for (int v=0;v<BC_nbc;v++) temp = NextPt(temp,b->ondir);
    (*l).push_back(temp);
    count++;
    ++c;
  } while (c!=b->data.end());

#ifdef TESTING
  cout <<"Got "<<count<<" cells, expected "<<b->data.size();
  cout <<"  list size = "<< (*l).size() <<"\n";
#endif
  *nc = count;


  return 0;
}





#endif // PARALLEL


