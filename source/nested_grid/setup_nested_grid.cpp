/// \file setup_nested_grid.cpp
/// 
/// \brief Class for setting up nested grids.
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

#include "setup_nested_grid.h"

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


setup_nested_grid::setup_nested_grid()
{
}



// ##################################################################
// ##################################################################


setup_nested_grid::~setup_nested_grid()
{
}




// ##################################################################
// ##################################################################



void setup_nested_grid::setup_nested_grid_levels(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // populate "levels" struct in SimPM based on nested grid parameters.
  //
  SimPM.nest_levels.clear();
  SimPM.nest_levels.resize(SimPM.grid_nlevels);

  for (int i=0;i<SimPM.grid_nlevels;i++) {
    for (int v=0;v<MAX_DIM;v++)
      SimPM.nest_levels[i].NG[v] = SimPM.NG[v];
    SimPM.nest_levels[i].Ncell = SimPM.Ncell;
    if (i==0) {
      for (int v=0;v<MAX_DIM;v++)
        SimPM.nest_levels[i].Range[v] = SimPM.Range[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.nest_levels[i].Xmin[v] = SimPM.Xmin[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.nest_levels[i].Xmax[v] = SimPM.Xmax[v];
      SimPM.nest_levels[i].dx = SimPM.Range[XX]/SimPM.NG[XX];
    }
    else {
      for (int v=0;v<MAX_DIM;v++) 
        SimPM.nest_levels[i].Range[v] = 0.5*SimPM.nest_levels[i-1].Range[v];
      for (int v=0;v<MAX_DIM;v++)
        SimPM.nest_levels[i].Xmin[v] = 0.5*(SimPM.nest_levels[i-1].Xmin[v]+SimPM.grid_nest_centre[v]);
      for (int v=0;v<MAX_DIM;v++)
        SimPM.nest_levels[i].Xmax[v] = 0.5*(SimPM.nest_levels[i-1].Xmax[v]+SimPM.grid_nest_centre[v]);
      SimPM.nest_levels[i].dx = 0.5*SimPM.nest_levels[i-1].dx;
    }
    SimPM.nest_levels[i].simtime = SimPM.simtime;
    SimPM.nest_levels[i].dt = 0.0;

    
    ostringstream temp; temp<<i;
    string lv = "level data"+temp.str();
    rep.printVec(lv,SimPM.nest_levels[i].Range,SimPM.ndim);
    rep.printVec(lv,SimPM.nest_levels[i].Xmin,SimPM.ndim);
    rep.printVec(lv,SimPM.nest_levels[i].Xmax,SimPM.ndim);
    cout <<"dx="<<SimPM.nest_levels[i].dx<<"\n";
  }
  return;
}





// ##################################################################
// ##################################################################



int setup_nested_grid::setup_grid(
      vector<class GridBaseClass *> &grid,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class MCMDcontrol * ///< unused for serial code.
      )
{
  cout <<"------------------------------------------------------\n";
  cout <<"------------  Setting up nested grid -----------------\n";

  if (SimPM.ndim <1 || SimPM.ndim>3)
    rep.error("Only know 1D,2D,3D methods!",SimPM.ndim);

#ifdef TESTING
  cout <<"Setting number of boundary cells == spatial OOA: ";
  cout <<SimPM.spOOA<<"\n";
#endif // TESTING

  //
  // Nbc is the depth of the boundary layer around each grid.
  //
  if      (SimPM.spOOA==OA2) SimPM.Nbc = 2;
  else if (SimPM.spOOA==OA1) SimPM.Nbc = 1;
  else rep.error("unhandles spatial order of accuracy",SimPM.spOOA);
  
  // Force Nbc=1 if using Lax-Friedrichs flux.
  if (SimPM.solverType==FLUX_LF)
  {SimPM.spOOA = SimPM.tmOOA = OA1; SimPM.Nbc=1;}

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set values cell interface class; note dx changes with level.
  //
  CI.set_dx((SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX]);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
  
  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(setup_nested_grid::setup_grid) Setting up grid...\n";
#endif
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    cout <<"Init: level="<< l <<",  &grid="<< &(grid[l])<<", and grid="<< grid[l] <<"\n";
    CI.set_dx(SimPM.nest_levels[l].dx);
    
    if (grid[l]) rep.error("Grid already set up!",grid[l]);

    if      (SimPM.coord_sys==COORD_CRT)
      grid[l] = new nested_grid (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.nest_levels[l].Xmin, SimPM.nest_levels[l].Xmax,
              SimPM.nest_levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_CYL)
      grid[l] = new nested_grid_cyl (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.nest_levels[l].Xmin, SimPM.nest_levels[l].Xmax,
              SimPM.nest_levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_SPH)
      grid[l] = new nested_grid_sph (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.nest_levels[l].Xmin, SimPM.nest_levels[l].Xmax,
              SimPM.nest_levels[l].NG, SimPM.Xmin, SimPM.Xmax);
    else 
      rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);

    if (grid[l]==0)
      rep.error("(setup_nested_grid::setup_grid) Couldn't assign data!", grid[l]);

#ifdef TESTING
    cout <<"(setup_nested_grid::setup_grid) Done. &grid=";
    cout << &(grid[l])<<", and grid="<<grid[l]<<"\n";
    cout <<"DX = "<<(grid[l])->DX()<<"\n";
    dp.grid = (grid[l]);
#endif
  }
  cout <<"------------------------------------------------------\n\n";

  return(0);
} // setup_grid()




// ##################################################################
// ##################################################################


int setup_nested_grid::boundary_conditions(
      vector<class GridBaseClass *> &grid,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      )
{
  // For uniform fixed cartesian grid.
#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<SimPM.Nbc<<"\n";
#endif
  int err = 0;
  for (int i=0;i<SimPM.grid_nlevels;i++) {

    if (l!=0) grid[l]->set_parent_grid(grid[l-1]);
    else      grid[l]->set_parent_grid(0);

    if (l!=SimPM.grid_nlevels-1) grid[l]->set_child_grid(grid[l+1]);
    else                         grid[l]->set_child_grid(0);

    err = grid[l]->SetupBCs(SimPM);
    rep.errorTest("setup_nested_grid::boundary_conditions()",0,err);
  }
#ifdef TESTING
  cout <<"(setup_nested_grid::boundary_conditions) Done.\n";
#endif
  return 0;
}



