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


#include "grid/nested_grid.h"
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
  chile  = 0;

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


