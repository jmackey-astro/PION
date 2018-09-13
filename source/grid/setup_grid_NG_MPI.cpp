/// \file setup_grid_NG_MPI.cpp
/// 
/// \brief Class for setting up parallel NG grids.
/// 
/// \author Jonathan Mackey
/// 
/// 
/// Modifications:
/// - 2018.08.24 JM: derived from setup_NG_grid and
///    setup_fixed_grid_MPI.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "raytracing/raytracer_SC.h"

#include "grid/setup_grid_NG_MPI.h"
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



setup_grid_NG_MPI::setup_grid_NG_MPI()
{
#ifdef TESTING
  cout <<"setup_grid_NG_MPI constructor called.\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



setup_grid_NG_MPI::~setup_grid_NG_MPI()
{
#ifdef TESTING
  cout <<"setup_grid_NG_MPI destructor called.\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



void setup_grid_NG_MPI::setup_NG_grid_levels(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  // call serial version to set global properties of each level.
  setup_NG_grid::setup_NG_grid_levels(SimPM);

  // For each level, do a domain decomposition
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err = SimPM.levels[l].MCMD.decomposeDomain(SimPM, SimPM.levels[l]);
    rep.errorTest("PLLEL Init():Decompose Domain!",0,err);
    SimPM.levels[l].MCMD.ReadSingleFile = true; // legacy option.
  }

  // For each level, get rank of process for parent grid and each
  // child grid.
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    SimPM.levels[l].MCMD.set_NG_hierarchy(SimPM,l);
  }
  
  return;
}



// ##################################################################
// ##################################################################



int setup_grid_NG_MPI::setup_grid(
      vector<class GridBaseClass *> &grid,  ///< grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  cout <<"------------------------------------------------------\n";
  cout <<"------------  Setting up MPI NG grid -----------------\n";

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
  else
    rep.error("unhandles spatial order of accuracy",SimPM.spOOA);
  
  //
  // May need to setup extra data in each cell.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set values cell interface class; note dx changes with level.
  //
  double dx=(SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX];
  CI.set_nlevels(dx,SimPM.grid_nlevels);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
  
  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(setup_grid_NG_MPI::setup_grid) Setting up grid...\n";
#endif
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    cout <<"Init: level="<< l <<",  &grid="<< &(grid[l]);
    cout <<", and grid="<< grid[l] <<"\n";
    
    if (grid[l]) rep.error("Grid already set up!",grid[l]);

    if      (SimPM.coord_sys==COORD_CRT)
      grid[l] = new UniformGridParallel (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].MCMD.LocalXmin,
              SimPM.levels[l].MCMD.LocalXmax,
              SimPM.levels[l].MCMD.LocalNG,
              SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_CYL)
      grid[l] = new uniform_grid_cyl_parallel (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].MCMD.LocalXmin,
              SimPM.levels[l].MCMD.LocalXmax,
              SimPM.levels[l].MCMD.LocalNG,
              SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys==COORD_SPH)
      grid[l] = new uniform_grid_sph_parallel (
              SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
              SimPM.levels[l].MCMD.LocalXmin,
              SimPM.levels[l].MCMD.LocalXmax,
              SimPM.levels[l].MCMD.LocalNG,
              SimPM.Xmin, SimPM.Xmax);
    else 
      rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);

    if (grid[l]==0)
      rep.error("(setup_grid_NG_MPI::setup_grid)", grid[l]);

#ifdef TESTING
    cout <<"(setup_grid_NG_MPI::setup_grid) Done. &grid=";
    cout << &(grid[l])<<", and grid="<<grid[l]<<"\n";
    cout <<"DX = "<<(grid[l])->DX()<<"\n";
    dp.grid = (grid[l]);
#endif
  }

  for (int l=0;l<SimPM.grid_nlevels;l++) {
    SimPM.levels[l].grid = grid[l];
    SimPM.levels[l].parent = 0; // this is used for serial code
    SimPM.levels[l].child = 0;  // this is used for serial code
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



int setup_grid_NG_MPI::setup_raytracing(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid ///< vec of grid pointers.
      )
{
  int err = 0;
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    cout <<"setting up raytracing for grid level "<<l<<"\n";
    err += setup_fixed_grid_pllel::setup_raytracing(SimPM,grid[l]);
    rep.errorTest("setup_grid_NG_MPI::setup_raytracing()",0,err);
  }
  
  err += setup_evolving_RT_sources(SimPM);
  rep.errorTest("setup_grid_NG_MPI::setup_evolving_RT_sources()",0,err);
  
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err += update_evolving_RT_sources(SimPM,grid[l]->RT);
    rep.errorTest("setup_grid_NG_MPI::update_evolving_RT_sources()",0,err);
  }
  return 0;
}




// ##################################################################
// ##################################################################



int setup_grid_NG_MPI::boundary_conditions(
      class SimParams &par,  ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid  ///< vec of grid pointers.
      )
{
#ifdef TESTING
  cout <<"Setting up BCs in MPI-NG Grid with Nbc="<<par.Nbc<<"\n";
#endif
  int err = setup_NG_grid::boundary_conditions(par,grid);
  rep.errorTest("setup_grid_NG_MPI::boundary_conditions",0,err);
#ifdef TESTING
  cout <<"(setup_grid_NG_MPI::boundary_conditions) Done.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int setup_grid_NG_MPI::setup_boundary_structs(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      const int l   ///< level in NG grid
      )
{
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif

  // first call fixed grid version
  int err = setup_fixed_grid_MPI::setup_boundary_structs(par,grid);
  rep.errorTest("sng::setup_boundary_structs fixed grid",0,err);



  //
  // Now check for NG grid boundaries if this grid has a parent
  // grid (i.e. if l > 0).
  //
  if (l>0) {
    // replace external boundary conditions with one that
    // receives data from a coarser level grid.
    for (int i=0; i<par.ndim; i++) {
      if (!pconst.equalD(par.levels[l-1].Xmin[i],
                         par.levels[l].Xmin[i])) {
#ifdef TESTING
        cout <<"reassigning neg. bc for axis "<<i<<" to COARSE_TO_FINE\n";
#endif
        grid->BC_bd[2*i]->itype = COARSE_TO_FINE_RECV;
        grid->BC_bd[2*i]->type  = "COARSE_TO_FINE_RECV";
      }
      if (!pconst.equalD(par.levels[l-1].Xmax[i],
                         par.levels[l].Xmax[i])) {
#ifdef TESTING
        cout <<"reassigning pos. bc for axis "<<i<<" to COARSE_TO_FINE\n";
#endif
        grid->BC_bd[2*i+1]->itype = COARSE_TO_FINE_RECV;
        grid->BC_bd[2*i+1]->type  = "COARSE_TO_FINE_RECV";
      }
    }

    // Add internal boundary to send averaged local data
    // to the parent grid.  Whether any data need to 
    // be sent/received is decided later in the function
    // assign_boundary_data().

#ifdef TESTING
    cout <<"Adding FINE_TO_COARSE_SEND boundary for level ";
    cout <<l<<", current # boundaries: "<<bd->data.size() <<"\n";
#endif
    struct boundary_data *bd = new boundary_data;
    bd->itype = FINE_TO_COARSE_SEND;
    bd->type  = "FINE_TO_COARSE_SEND";
    bd->dir = NO;
    bd->ondir = NO;
    bd->refval=0;
    grid->BC_bd.push_back(bd);
  }

  //
  // Now check for NG grid boundaries if this grid has a child
  // grid (i.e. if l < nlevels-1) and, if so, add boundary structs
  // to receive averaged data from a finer grid, and to send 
  // external boundary data to the finer grid.  Whether any data
  // need to be sent/received is decided later in the function
  // assign_boundary_data().
  //
  if (l < par.grid_nlevels-1) {
#ifdef TESTING
    cout <<"Adding FINE_TO_COARSE_RECV boundary for level ";
    cout <<l<<", current # boundaries: "<<bd->data.size() <<"\n";
#endif
    struct boundary_data *bd = new boundary_data;
    bd->itype = FINE_TO_COARSE_RECV;
    bd->type  = "FINE_TO_COARSE_RECV";
    bd->dir = NO;
    bd->ondir = NO;
    bd->refval=0;
    grid->BC_bd.push_back(bd);

#ifdef TESTING
    cout <<"Adding COARSE_TO_FINE_SEND boundary for level ";
    cout <<l<<", current # boundaries: "<<bd->data.size() <<"\n";
#endif
    struct boundary_data *bd2 = new boundary_data;
    bd2->itype = COARSE_TO_FINE_SEND;
    bd2->type  = "COARSE_TO_FINE_SEND";
    bd2->dir = NO;
    bd2->ondir = NO;
    bd2->refval=0;
    grid->BC_bd.push_back(bd2);
  }
  
  
#ifdef TESTING
  cout <<"BC structs set up.\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



void setup_grid_NG_MPI::setup_dataio_class(
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
    dataio = new dataio_silo_utility (par, "DOUBLE");
    break; 
#endif // if SILO

  default:
    rep.error("sim_control_NG::Init unhandled filetype",typeOfFile);
  }
  return;
}



// ##################################################################
// ##################################################################


