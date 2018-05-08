///
/// \file nested_grid.cpp
/// 
/// \author Jonathan Mackey
/// 
/// Function definitions for nested_grid class:
/// - setup boundaries.
/// - update boundaries.
/// 
/// Modifications:\n
///  - 2018.05.03 JM: started on code.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include "nested_grid/nested_grid.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#include <fstream>
#include <iostream>
using namespace std;


//#define GLM_ZERO_BOUNDARY ///< Set this flag to make Psi=0 on boundary cells.
#define GLM_NEGATIVE_BOUNDARY ///< Set this flag for Psi[boundary cell]=-Psi[edge cell]



// ##################################################################
// ##################################################################



nested_grid::nested_grid(
    int nd,
    int nv,
    int eqt,
    int Nbc,       ///< Number of boundary cells to use.
    double *g_xn,
    double *g_xp,
    int *g_nc,
    double *sim_xn,
    double *sim_xp
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp)
{

#ifdef TESTING
  cout <<"Setting up nested_grid with G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif
#ifdef TESTING
  rep.printVec("nested_grid ixmin ", G_ixmin, G_ndim);
  rep.printVec("nested_grid ixmax ", G_ixmax, G_ndim);
  rep.printVec("nested_grid irange", G_irange,G_ndim);
  rep.printVec("nested_grid xmin ", G_xmin, G_ndim);
  rep.printVec("nested_grid xmax ", G_xmax, G_ndim);
  rep.printVec("nested_grid range", G_range,G_ndim);
#endif

  parent = 0;
  child  = 0;

#ifdef TESTING
  cout <<"nested_grid Constructor done.\n";
#endif
} //nested_grid Constructor



// ##################################################################
// ##################################################################


nested_grid::~nested_grid()
{
#ifdef TESTING
  cout <<"nested_grid Destructor.\n";
#endif
} // Destructor



// ##################################################################
// ##################################################################



/// Set pointer to parent grid.
void nested_grid::set_parent_grid(
      class GridBaseClass *gp  ///< pointer to parent grid.
      )
{
  parent = gp;
  return;
}



// ##################################################################
// ##################################################################



void nested_grid::set_child_grid(
      class GridBaseClass *gp  ///< pointer to child grid.
      )
{
  child = gp;
  return;
}



// ##################################################################
// ##################################################################



int nested_grid::BC_setBCtypes(
      class SimParams &par  ///< List of simulation params (including BCs)
      )
{
  int err = UniformGrid::BC_setBCtypes(par);
  rep.errorTest("nested_grid::BC_setBCtypes() unigrid call",0,err);

  if (parent) set_external_bcs_from_parent();

  if (child) set_nonleaf_cells_as_BD();
  return 0;
}



// ##################################################################
// ##################################################################



void nested_grid::set_external_bcs_from_parent()
{
  //
  // For serial code, we can grab the data from the parent level and
  // then refine it and assign it to the current level, if the level
  // boundary is not located at the global boundary.  We re-label
  // external boundaries as NEST_COARSE if we will get the boundary
  // data from a coarser grid.
  //
  double pxmin[MAX_DIM], pxmax[MAX_DIM];

  for (int v=0;v<G_ndim;v++) {
    pxmin[v] = parent->Xmin(static_cast<axes>(v));
    pxmax[v] = parent->Xmax(static_cast<axes>(v));
  }
  int i=0;
  for (i=0; i<G_ndim; i++) {
    if (!pconst.equalD(pxmin[i],G_xmin[i])) {
      cout <<"reassigning neg. bc for axis "<<i<<" to NEST_COARSE\n";
      BC_bd[2*i]->itype = NEST_COARSE;
      BC_bd[2*i]->type  = "NEST_COARSE";
    }
    if (!pconst.equalD(pxmax[i],G_xmax[i])) {
      cout <<"reassigning pos. bc for axis "<<i<<" to NEST_COARSE\n";
      BC_bd[2*i+1]->itype = NEST_COARSE;
      BC_bd[2*i+1]->type  = "NEST_COARSE";
    }
  }
  return;
}


// ##################################################################
// ##################################################################



void nested_grid::set_nonleaf_cells_as_BD()
{
  //
  // for serial code, a grid can only have one child grid, because
  // a single grid covers the full domain at each level.  So we go
  // through all cells on the grid, and if they are within the domain
  // of the child grid (if it exists) then we label the cells as
  // boundary data and add them to a new NEST_FINE boundary.
  //
  double cxmin[MAX_DIM], cxmax[MAX_DIM], cpos[MAX_DIM];
  bool within_child=true;
  struct boundary_data *bd = new boundary_data;

  for (int v=0;v<G_ndim;v++) {
    cxmin[v] = child->Xmin(static_cast<axes>(v));
    cxmax[v] = child->Xmax(static_cast<axes>(v));
  }

  bd->itype = NEST_FINE;
  bd->type  = "NEST_FINE";
  bd->dir = NO;
  bd->ondir = NO;
  bd->refval=0;
  
  cell *c = FirstPt();
  do {
    within_child=true;
    CI.get_dpos(c,cpos);
    for (int v=0;v<G_ndim;v++) {
      if (cpos[v]<cxmin[v] || cpos[v]>cxmax[v]) within_child=false;
    }
    if (within_child) {
      c->isbd = true;
      bd->data.push_back(c);
    }
  } while ((c=NextPt(c)) !=0);

  BC_bd.push_back(bd);

  return;
}




// ##################################################################
// ##################################################################



int nested_grid::assign_boundary_data(
        const double simtime,     ///< current simulation time
        const double sim_start,   ///< start time of simulation
        const double sim_finish,  ///< finish time of simulation
        const double Tmin         ///< minimum temperature allowed
        )
{
  // first call the UniformGrid version.
  int err = UniformGrid::assign_boundary_data(simtime, sim_start,
                                              sim_finish, Tmin);

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



int nested_grid::BC_assign_NEST_FINE(
      boundary_data *b
      )
{
  return 0;
}



// ##################################################################
// ##################################################################



int nested_grid::BC_assign_NEST_COARSE(
      boundary_data *b
      )
{
  return 0;
}



// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################




nested_grid_cyl::nested_grid_cyl(
    int nd,         ///< ndim, length of position vector.
    int nv,         ///< nvar, length of state vectors.
    int eqt,        ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,        ///< Number of boundary cells to use.
    double *g_xn,   ///< array of minimum values of x,y,z for this grid.
    double *g_xp,   ///< array of maximum values of x,y,z for this grid.
    int *g_nc,      ///< array of number of cells in x,y,z directions.
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  nested_grid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  uniform_grid_cyl(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp)
{
}



// ##################################################################
// ##################################################################



nested_grid_cyl::~nested_grid_cyl()
{
#ifdef TESTING
  cout <<"nested_grid_cyl destructor.\n";
#endif
}



// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################



nested_grid_sph::nested_grid_sph(
    int nd,         ///< ndim, length of position vector.
    int nv,         ///< nvar, length of state vectors.
    int eqt,        ///< eqntype, which equations we are using (needed by BCs).
    int Nbc,        ///< Number of boundary cells to use.
    double *g_xn,   ///< array of minimum values of x,y,z for this grid.
    double *g_xp,   ///< array of maximum values of x,y,z for this grid.
    int *g_nc,      ///< array of number of cells in x,y,z directions.
    double *sim_xn, ///< array of min. x/y/z for full simulation.
    double *sim_xp  ///< array of max. x/y/z for full simulation.
    )
  : 
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  nested_grid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd),
  uniform_grid_sph(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp)
{
}



// ##################################################################
// ##################################################################



nested_grid_sph::~nested_grid_sph()
{
#ifdef TESTING
  cout <<"nested_grid_sph destructor.\n";
#endif
}



// ##################################################################
// ##################################################################



