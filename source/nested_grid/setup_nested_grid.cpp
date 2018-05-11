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

#include "nested_grid/setup_nested_grid.h"
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
  SimPM.levels.clear();
  SimPM.levels.resize(SimPM.grid_nlevels);

  for (int i=0;i<SimPM.grid_nlevels;i++) {
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
        SimPM.levels[i].Xmin[v] = 0.5*(SimPM.levels[i-1].Xmin[v]+SimPM.grid_nest_centre[v]);
      for (int v=0;v<MAX_DIM;v++)
        SimPM.levels[i].Xmax[v] = 0.5*(SimPM.levels[i-1].Xmax[v]+SimPM.grid_nest_centre[v]);
      SimPM.levels[i].dx = 0.5*SimPM.levels[i-1].dx;
    }
    SimPM.levels[i].simtime = SimPM.simtime;
    SimPM.levels[i].dt = 0.0;

    
    ostringstream temp; temp<<i;
    string lv = "level data"+temp.str();
    rep.printVec(lv,SimPM.levels[i].Range,SimPM.ndim);
    rep.printVec(lv,SimPM.levels[i].Xmin,SimPM.ndim);
    rep.printVec(lv,SimPM.levels[i].Xmax,SimPM.ndim);
    cout <<"dx="<<SimPM.levels[i].dx<<"\n";
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
  double dx=(SimPM.Xmax[XX]-SimPM.Xmin[XX])/SimPM.NG[XX];
  CI.set_dx(dx);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
  CI.set_nlevels(dx,SimPM.grid_nlevels);
  
  //
  // Now we can setup the grid:
  //
#ifdef TESTING
  cout <<"(setup_nested_grid::setup_grid) Setting up grid...\n";
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
      class SimParams &par,  ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  // For uniform fixed cartesian grid.
//#ifdef TESTING
  cout <<"Setting up BCs in Grid with Nbc="<<par.Nbc<<"\n";
//#endif
  int err = 0;
  bdata_nest.resize(par.grid_nlevels);

  for (int l=0;l<par.grid_nlevels;l++) {
    cout <<"level "<<l<<", setting up boundaries\n";

    //
    // Choose what BCs to set up based on BC strings.
    //
    int err = setup_boundary_structs(par,grid[l],l,bdata_nest[l]);
    rep.errorTest("sng::boundary_conditions sb_structs",0,err);


    //if (l!=0) grid[l]->set_parent_grid(grid[l-1]);
    //else      grid[l]->set_parent_grid(0);

    //if (l!=par.grid_nlevels-1) grid[l]->set_child_grid(grid[l+1]);
    //else                         grid[l]->set_child_grid(0);

    err = grid[l]->SetupBCs(par,bdata_nest[l]);
    rep.errorTest("sng::boundary_conditions SetupBCs",0,err);
  }
//#ifdef TESTING
  cout <<"(setup_nested_grid::boundary_conditions) Done.\n";
//#endif
  return 0;
}



// ##################################################################
// ##################################################################




int setup_nested_grid::setup_boundary_structs(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< pointer to grid.
      const int l,
      std::vector<struct boundary_data *> &BC_bd  ///< pointer to boundary data vector for this level
      )
{
#ifdef TESTING
  cout <<"Set BC types...\n";
#endif

  // first call fixed grid version
  int err = setup_fixed_grid::setup_boundary_structs(par,grid,BC_bd);

  //
  // Now check for nested grid boundaries if this grid has a parent
  // grid (i.e. if l > 0).
  //
  if (l>0) {
    for (int i=0; i<par.ndim; i++) {
      if (!pconst.equalD(par.levels[l-1].Xmin[i],par.levels[l].Xmin[i])) {
        cout <<"reassigning neg. bc for axis "<<i<<" to NEST_COARSE\n";
        BC_bd[2*i]->itype = NEST_COARSE;
        BC_bd[2*i]->type  = "NEST_COARSE";
      }
      if (!pconst.equalD(par.levels[l-1].Xmax[i],par.levels[l].Xmax[i])) {
        cout <<"reassigning pos. bc for axis "<<i<<" to NEST_COARSE\n";
        BC_bd[2*i+1]->itype = NEST_COARSE;
        BC_bd[2*i+1]->type  = "NEST_COARSE";
      }
    }
  }

  //
  // Now check for nested grid boundaries if this grid has a child
  // grid (i.e. if l < nlevels), and set non-leaf cells to boundary
  // data.
  //
  // for serial code, a grid can only have one child grid, because
  // a single grid covers the full domain at each level.  So we go
  // through all cells on the grid, and if they are within the domain
  // of the child grid (if it exists) then we label the cells as
  // boundary data and add them to a new NEST_FINE boundary.
  //
  if (l < par.grid_nlevels) {
    double cxmin[MAX_DIM], cxmax[MAX_DIM], cpos[MAX_DIM];
    bool within_child=true;
    for (int v=0;v<par.ndim;v++) {
      cxmin[v] = par.levels[l+1].Xmin[v];
      cxmax[v] = par.levels[l+1].Xmax[v];
    }
    struct boundary_data *bd = new boundary_data;
    bd->itype = NEST_FINE;
    bd->type  = "NEST_FINE";
    bd->dir = NO;
    bd->ondir = NO;
    bd->refval=0;
    cell *c = grid->FirstPt();
    do {
      within_child=true;
      CI.get_dpos(c,cpos);
      for (int v=0;v<par.ndim;v++) {
        if (cpos[v]<cxmin[v] || cpos[v]>cxmax[v]) within_child=false;
      }
      if (within_child) {
        c->isbd = true;
        bd->data.push_back(c);
      }
    } while ((c=grid->NextPt(c)) !=0);

    BC_bd.push_back(bd);
  }
  
#ifdef TESTING
  cout <<"BC structs set up.\n";
#endif
  return 0;
}





// ##################################################################
// ##################################################################



int setup_nested_grid::assign_boundary_data(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      std::vector<struct boundary_data *> &BC_bd ///< pointer to boundary structs
      )
{
  // THIS NEEDS THE LEVELS AND POINTERS TO PARENT AND CHILD GRIDS!
  // first call the UniformGrid version.
  int err = UniformGrid::assign_boundary_data(par,grid,BC_bd);

  //
  // Then check for nested-grid boundaries and assign data for them.
  //
  for (int i=0; i<BC_nbd; i++) {
    switch (BC_bd[i]->itype) {
      case NEST_FINE:   err += BC_assign_NEST_FINE(  BC_bd[i]);break;
      case NEST_COARSE: err += BC_assign_NEST_COARSE(BC_bd[i]);break;
      default:
      break;
    }
  }
  return err;
}



// ##################################################################
// ##################################################################



int setup_nested_grid::BC_assign_NEST_FINE(
      boundary_data *b
      )
{
  //
  // Need to give each of these cells a pointer to a cell on the
  // finer level grid, if possible.
  //
  if (b->data.empty())
    rep.error("BC_assign_NEST_COARSE: empty boundary data",b->itype);
  b->nest.clear();

  list<cell*>::iterator bpt=b->data.begin();
  class GridBaseClass *cg = child; // child (finer) grid.
  cell *cc = cg->FirstPt_All(); // child cell.
  int cdx = 0.5*cg->idx();
  double distance =  0.0;

  // Map each bpt cell to a cell in b->nest list, which is the first
  // cell in the finer grid that is part of the coarse cell (i.e. the
  // one with the most negative coordinates).
  do{
    cc = cg->FirstPt_All();
    for (int v=0;v<G_ndim;v++) {
      while (cc && cc->pos[v] < (*bpt)->pos[v]-cdx)
        cc = cg->NextPt(cc,static_cast<direction>(2*v+1));
    }
    if (!cc) rep.error("BC_assign_NEST_FINE: lost on fine grid",0);
    b->nest.push_back(cc);
    
    ++bpt;
  }  while (bpt !=b->data.end());

  return 0;
}



// ##################################################################
// ##################################################################



int setup_nested_grid::BC_assign_NEST_COARSE(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *b
      )
{
  //
  // Make a list of pointers to cells in the coarser grid that map
  // onto this (finer) grid external boundary, and then write an
  // alogrithm to interpolate the coarse data onto the finer grid.
  //
  if (b->data.empty())
    rep.error("BC_assign_NEST_COARSE: empty boundary data",b->itype);
  b->nest.clear();

  list<cell*>::iterator bpt=b->data.begin();
  class GridBaseClass *pg = parent; // parent (coarser) grid.
  int pidx = pg->idx();
  //cout <<"BC_assign_NEST_COARSE: dx="<<G_idx<<", parent dx="<<pidx<<"\n";

  cell *pc = pg->FirstPt_All(); // parent cell.

  double distance =  0.0;
  bool loop;
  
  do{
    loop = false;
    distance = idistance(pc->pos, (*bpt)->pos);
    // Find parent cell that covers this boundary cell.  It should be
    // G_idx/2 away from the boundary cell in each direction.
    //rep.printVec("bpt pos",(*bpt)->pos,G_ndim);
    while (distance > G_idx && pc!=0) {
      //cout <<"distance="<<distance<<"; "; rep.printVec("pc pos",pc->pos,G_ndim);
      pc = pg->NextPt_All(pc);
      if (!pc && !loop) { // hack: if get to the end, then go back...
        pc = b->nest.front();
        loop = true;
      }
      distance = idistance(pc->pos, (*bpt)->pos);
    }
    if (!pc) rep.error("BC_assign_NEST_COARSE() left parent grid",0);
    
    // add this parent cell to the "parent" list of this boundary.
    b->nest.push_back(pc);
    (*bpt)->npt = pc;
    ++bpt;
  }  while (bpt !=b->data.end());

  // add data to boundary cells.
  //BC_update_NEST_COARSE(b,OA2,OA2);

  return 0;
}





// ##################################################################
// ##################################################################



int setup_nested_grid::setup_raytracing(
      class SimParams &SimPM,    ///< pointer to simulation parameters
      vector<class GridBaseClass *> &grid,  ///< address of vector of grid pointers.
      vector<class RayTracingBase *> &RT  ///< address of vector of grid pointers.
      )
{
  int err = 0;
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err += setup_fixed_grid::setup_raytracing(SimPM,grid[l],&(RT[l]));
    err += setup_evolving_RT_sources(SimPM,RT[l]);
    rep.errorTest("setup_nested_grid::setup_raytracing()",0,err);
  }
  return 0;
}


