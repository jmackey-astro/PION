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
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cart(nd)
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
  cout <<"nested_grid Destructor\n";
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



/// Set pointer to child grid.
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

///
/// Constructor
///
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
  nested_grid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  uniform_grid_cyl(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd)
{
}


// ##################################################################
// ##################################################################


nested_grid_cyl::~nested_grid_cyl()
{
#ifdef TESTING
  cout <<"nested_grid_cyl destructor. Present and correct!\n";
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
  nested_grid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  UniformGrid(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  uniform_grid_sph(nd,nv,eqt,Nbc,g_xn,g_xp,g_nc,sim_xn,sim_xp),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd)
{
}



// ##################################################################
// ##################################################################



nested_grid_sph::~nested_grid_sph()
{
#ifdef TESTING
  cout <<"nested_grid_sph destructor. Present and correct!\n";
#endif
}



// ##################################################################
// ##################################################################



