/// \file setup_NG_grid.cpp
/// 
/// \brief Class for setting up NG grids.
/// 
/// \author Jonathan Mackey
/// 
/// 
/// Modifications:
/// - 2018.05.01 JM derived class from setup_fixed_grid.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "raytracing/raytracer_SC.h"

#include "grid/setup_NG_grid.h"
#include "grid/uniform_grid.h"

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"
#include "microphysics/MPv3.h"
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"
#include "microphysics/MPv8.h"

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/MPv9.h"
#endif

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif



#include "dataIO/dataio_base.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS


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



setup_NG_grid::setup_NG_grid()
{
}



// ##################################################################
// ##################################################################



setup_NG_grid::~setup_NG_grid()
{
}



// ##################################################################
// ##################################################################



void setup_NG_grid::setup_NG_grid_levels(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // populate "levels" struct in SimPM based on NG grid parameters.
  //
  SimPM.levels.clear();
  SimPM.levels.resize(SimPM.grid_nlevels);

  for (int i=0;i<SimPM.grid_nlevels;i++) {
    SimPM.levels[i].parent=0;
    SimPM.levels[i].child=0;
    for (int v=0;v<MAX_DIM;v++)
      SimPM.levels[i].NG[v] = SimPM.NG[v];
    SimPM.levels[i].Ncell = SimPM.Ncell;
    if (i==0) {
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Range[v] = SimPM.Range[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Xmin[v] = SimPM.Xmin[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Xmax[v] = SimPM.Xmax[v];
      SimPM.levels[i].dx = SimPM.Range[XX]/SimPM.NG[XX];
    }
    else {
      for (int v=0;v<MAX_DIM;v++) 
        SimPM.levels[i].Range[v] = 0.5*SimPM.levels[i-1].Range[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Xmin[v] = 0.5*(SimPM.levels[i-1].Xmin[v]+SimPM.NG_centre[v]);
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Xmax[v] = 0.5*(SimPM.levels[i-1].Xmax[v]+SimPM.NG_centre[v]);
      SimPM.levels[i].dx = 0.5*SimPM.levels[i-1].dx;
    }
    SimPM.levels[i].simtime = SimPM.simtime;
    SimPM.levels[i].dt = 0.0;

    if (i==0) SimPM.levels[i].multiplier = 1;
    else      SimPM.levels[i].multiplier = 2*SimPM.levels[i-1].multiplier;
  }

  for (int i=SimPM.grid_nlevels-1;i>=0; i--) {
    if (i==SimPM.grid_nlevels-1)
      SimPM.levels[i].step = SimPM.timestep;
    else
      SimPM.levels[i].step = SimPM.levels[i+1].step/2;
    
    ostringstream temp; temp<<i;
    string lv = "level "+temp.str();
    string t2=lv+"_Range";
    rep.printVec(t2,SimPM.levels[i].Range,SimPM.ndim);
    t2 = lv+"_Xmin";
    rep.printVec(t2,SimPM.levels[i].Xmin,SimPM.ndim);
    t2 = lv+"_Xmax";
    rep.printVec(t2,SimPM.levels[i].Xmax,SimPM.ndim);
    cout <<"dx="<<SimPM.levels[i].dx;
    cout <<", step="<<SimPM.levels[i].step<<"\n";
  }

  return;
}



// ##################################################################
// ##################################################################



int setup_NG_grid::setup_grid(
      vector<class GridBaseClass *> &grid,  ///< grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  cout <<"------------------------------------------------------\n";
  cout <<"----------------  Setting up NG grid -----------------\n";

  if (SimPM.ndim <1 || SimPM.ndim>3)
    rep.error("Only know 1D,2D,3D methods!",SimPM.ndim);

#ifdef TESTING
  cout <<"Setting number of boundary cells == spatial OOA: ";
  cout <<SimPM.spOOA<<"\n";
#endif // TESTING

  //
  // Nbc is the depth of the boundary layer around each grid.
  //
  if      (SimPM.spOOA==OA2) SimPM.Nbc = 4;
  else if (SimPM.spOOA==OA1) SimPM.Nbc = 2;
  else rep.error("unhandles spatial order of accuracy",SimPM.spOOA);
  
  //
  // May need to setup extra data in each cell.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set values cell interface class.
  //
  double dx=(SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX];
  CI.set_nlevels(dx,SimPM.grid_nlevels); // sets dx on all levels.
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
  
  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(setup_NG_grid::setup_grid) Setting up grid...\n";
#endif
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    cout <<"Init: level="<< l <<",  &grid="<< &(grid[l])<<", and grid="<< grid[l] <<"\n";
    
    if (grid[l]) rep.error("Grid already set up!",grid[l]);

    if      (SimPM.coord_sys==COORD_CRT)
      grid[l] = new UniformGrid (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].Xmin, SimPM.levels[l].Xmax,
              SimPM.levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_CYL)
      grid[l] = new uniform_grid_cyl (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].Xmin, SimPM.levels[l].Xmax,
              SimPM.levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_SPH)
      grid[l] = new uniform_grid_sph (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].Xmin, SimPM.levels[l].Xmax,
              SimPM.levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else 
      rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);

    if (grid[l]==0)
      rep.error("(setup_NG_grid::setup_grid) Couldn't assign data!", grid[l]);

#ifdef TESTING
    cout <<"(setup_NG_grid::setup_grid) Done. &grid=";
    cout << &(grid[l])<<", and grid="<<grid[l]<<"\n";
    cout <<"DX = "<<(grid[l])->DX()<<"\n";
    dp.grid = (grid[l]);
#endif
  }

  for (int l=0;l<SimPM.grid_nlevels;l++) {
    SimPM.levels[l].grid = grid[l];
    if (l==0) {
      SimPM.levels[l].parent = 0;
      SimPM.levels[l].child  = grid[l+1];
    }
    else if (l==SimPM.grid_nlevels-1) {
      SimPM.levels[l].parent = grid[l-1];
      SimPM.levels[l].child  = 0;
    }
    else {
      SimPM.levels[l].parent = grid[l-1];
      SimPM.levels[l].child  = grid[l+1];
    }
  }

  // setup arrays for fluxes into and out of fine grid, and
  // equivalent cells on coarse grid, for making the fluxes
  // consistent across levels.
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    if (l!=0) grid[l]->setup_flux_send(SimPM,l-1);
    if (l!=SimPM.grid_nlevels-1) grid[l]->setup_flux_recv(SimPM,l+1);
  }

cout <<"------------------------------------------------------\n\n";

  return(0);
} // setup_grid()



// ##################################################################
// ##################################################################



int setup_NG_grid::setup_raytracing(
      class SimParams &SimPM,    ///< simulation parameters
      vector<class GridBaseClass *> &grid  ///< vector of grids.
      )
{
  int err = 0;
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    cout <<"setting up raytracing for grid level "<<l<<"\n";
    err += setup_fixed_grid::setup_raytracing(SimPM,grid[l]);
    rep.errorTest("setup_NG_grid::setup_raytracing()",0,err);
  }
  
  err += setup_evolving_RT_sources(SimPM);
  rep.errorTest("setup_NG_grid::setup_evolving_RT_sources()",0,err);
  
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err += update_evolving_RT_sources(SimPM,SimPM.levels[l].simtime, 
                                                        grid[l]->RT);
    rep.errorTest("setup_NG_grid::update_RT_sources()",0,err);
  }
  return 0;
}



// ##################################################################
// ##################################################################



int setup_NG_grid::boundary_conditions(
      class SimParams &par,  ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<par.Nbc<<"\n";
#endif
  int err = 0;
  for (int l=0;l<par.grid_nlevels;l++) {
    cout <<"level "<<l<<", setting up boundaries\n";

    //
    // Choose what BCs to set up based on BC strings.
    //
    err = setup_boundary_structs(par,grid[l],l);
    rep.errorTest("sng::boundary_conditions sb_structs",0,err);

    err = grid[l]->SetupBCs(par);
    rep.errorTest("sng::boundary_conditions SetupBCs",0,err);
  }
#ifdef TESTING
  cout <<"(setup_NG_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int setup_NG_grid::setup_boundary_structs(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      const int l   ///< level in NG grid
      )
{
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif

  // first call fixed grid version
  int err = setup_fixed_grid::setup_boundary_structs(par,grid);
  rep.errorTest("sng::setup_boundary_structs fixed grid",0,err);

  //
  // Now check for NG grid boundaries if this grid has a parent
  // grid (i.e. if l > 0).
  //
  if (l>0) {
    for (int i=0; i<par.ndim; i++) {
      if (!pconst.equalD(par.levels[l-1].Xmin[i],par.levels[l].Xmin[i])) {
#ifdef TESTING
        cout <<"reassigning neg. bc for axis "<<i<<" to COARSE_TO_FINE\n";
#endif
        grid->BC_bd[2*i]->itype = COARSE_TO_FINE;
        grid->BC_bd[2*i]->type  = "COARSE_TO_FINE";
      }
      if (!pconst.equalD(par.levels[l-1].Xmax[i],par.levels[l].Xmax[i])) {
#ifdef TESTING
        cout <<"reassigning pos. bc for axis "<<i<<" to COARSE_TO_FINE\n";
#endif
        grid->BC_bd[2*i+1]->itype = COARSE_TO_FINE;
        grid->BC_bd[2*i+1]->type  = "COARSE_TO_FINE";
      }
    }
  }

  //
  // Now check for NG grid boundaries if this grid has a child
  // grid (i.e. if l < nlevels), and set non-leaf cells to boundary
  // data.
  //
  // for serial code, a grid can only have one child grid, because
  // a single grid covers the full domain at each level.  So we go
  // through all cells on the grid, and if they are within the domain
  // of the child grid (if it exists) then we label the cells as
  // boundary data and add them to a new FINE_TO_COARSE boundary.
  //
  if (l < par.grid_nlevels-1) {
    double cxmin[MAX_DIM], cxmax[MAX_DIM], cpos[MAX_DIM];
    bool within_child=true;
    size_t ct = 0;
    for (int v=0;v<par.ndim;v++) {
      cxmin[v] = par.levels[l+1].Xmin[v];
      cxmax[v] = par.levels[l+1].Xmax[v];
    }
    struct boundary_data *bd = new boundary_data;
    bd->itype = FINE_TO_COARSE;
    bd->type  = "FINE_TO_COARSE";
    bd->dir = NO;
    bd->ondir = NO;
    bd->refval=0;
    bd->NGrecvF2C.resize(1);
    cell *c = grid->FirstPt();
    do {
      within_child=true;
      CI.get_dpos(c,cpos);
      for (int v=0;v<par.ndim;v++) {
        if (cpos[v]<cxmin[v] || cpos[v]>cxmax[v]) within_child=false;
      }
      if (within_child) {
        c->isbd = true;
        //cout <<"Nested-grid: FINE_TO_COARSE setup, cell added.\n";
        bd->NGrecvF2C[0].push_back(c);
        ct++;
      }
    } while ((c=grid->NextPt(c)) !=0);
#ifdef TESTING
    cout <<"Got "<<ct<<" cells for FINE_TO_COARSE boundary, "<<bd->data.size() <<"\n";
#endif
    grid->BC_bd.push_back(bd);
#ifdef TESTING
    cout <<"BC_data: "<<;
    cout << grid->BC_bd[grid->BC_bd.size()-1]->NGrecvF2C[0].size();
    cout << <<"\n";
#endif
  }
  
  
#ifdef TESTING
  cout <<"BC structs set up.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



void setup_NG_grid::setup_dataio_class(
      class SimParams &par,     ///< simulation parameters
      const int typeOfFile ///< type of I/O: 1=text,2=fits,5=silo
      )
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

#ifdef SILO
  case 5: // Start from Silo snapshot.
    dataio = new dataio_silo (par, "DOUBLE");
    break; 
#endif // if SILO

  default:
    rep.error("sim_control_NG::Init unhandled filetype",typeOfFile);
  }
  return;
}




// ##################################################################
// ##################################################################



