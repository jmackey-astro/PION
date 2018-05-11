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



// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################




// ##################################################################
// ##################################################################




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



