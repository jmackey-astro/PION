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
  // This code assumes that the UniformGrid class has setup the
  // boundary as a regular external boundary and added cells to the
  // list of data.
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



int nested_grid::BC_update_NEST_FINE(
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // This is relatively straighforward, in that we just weight each
  // fine cell by its volume.  Assume there are two cells in each
  // dimension in the fine grid, and so we can loop over this.
  //
  list<cell*>::iterator coarse=b->data.begin();
  list<cell*>::iterator fine=b->nest.begin();
  cell *c, *f;
  double c[G_nvar];

  for (coarse=b->data.begin(); coarse!=b->data.end(); ++coarse) {
    c = (*coarse);
    f = (*fine);

    // 1D
    for (int i=0;i<2;i++) {
      // Need to convert to conserved variables!!!!  Can only do this
      // in sim_control_nested.


    // constant data:
    (*c)->Ph[v] = (*c)->npt->Ph[v]
    for (int v=0;v<G_nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    }

    ++fine;
  }
  return 0;
}



// ##################################################################
// ##################################################################



int nested_grid::BC_assign_NEST_COARSE(
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
  BC_update_NEST_COARSE(b,OA2,OA2);

  return 0;
}



// ##################################################################
// ##################################################################



/// Updates data to an external boundary from coarser grid.
int nested_grid::BC_update_NEST_COARSE(
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // This is a complicated problem to use linear interpolation (or
  // higher order), because you have to conserve mass, momentum and
  // energy between the two levels.  It should also be monotonic and
  // produce cell-averaged values. Flash, amrvac and pluto all seem
  // to use very old Fortran code to accomplish this, and the code is
  // almost impenetrable.
  //
  // For now we just do contant data.  Will fix this eventually once
  // the rest of the code is working.  It actually doesn't really
  // matter for expanding nebulae because there is supersonic outflow
  // but this really needs to be improved for solving more general
  // problems.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    // set slope along x-direction in parent cell.
    // parent->SetSlope((*c)->npt,XX,G_nvar,sx,2,parent);
    // get physical offset distance between cell and parent.
    //dist = idifference_cell2cell(*c,(*c)->npt,XX)*CI.phys_per_int();
    // interpolate linearly to cell position.
    //for (int v=0;v<G_nvar;v++)
    //  (*c)->Ph[v] = (*c)->npt->Ph[v] + sx[v]*dist;

    // constant data:
    (*c)->Ph[v] = (*c)->npt->Ph[v]
    for (int v=0;v<G_nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    }
  }
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



