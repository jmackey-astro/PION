/// \file assign_update_bcs.cpp
/// \brief Defines a class that inherits all boundary types and 
///   assigns and updated the boundaries.
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "assign_update_bcs.h"
using namespace std;




// ##################################################################
// ##################################################################



int assign_update_bcs::assign_boundary_data(
      class SimParams &par,     ///< simulation parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid  ///< pointer to grid.
      )
{
  int err=0;
  struct boundary_data *b;
  //
  // Loop through all boundaries, and assign data to them.
  //
  for (size_t i=0; i<grid->BC_bd.size(); i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case PERIODIC:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
     err += BC_assign_PERIODIC(  par,level,grid,b);
     break;
    case OUTFLOW:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_OUTFLOW(   par,grid,b);
      break;
    case ONEWAY_OUT:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_ONEWAY_OUT(par,grid,b);
      break;
    case INFLOW:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_INFLOW(    par,grid,b);
      break;
    case REFLECTING:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_REFLECTING(par,grid,b);
      break;
    case AXISYMMETRIC:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_AXISYMMETRIC(par,grid,b);
      break;
    case FIXED:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_FIXED(     par,grid,b);
      break;
    case JETBC:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_JETBC(     par,grid,b);
      break;
    case JETREFLECT:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_JETREFLECT(par,grid,b);
      break;
    case DMACH:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_DMACH(     par,grid,b);
      break;
    case DMACH2:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_DMACH2(    par,grid,b);
      break;
    case STWIND:
#ifdef TEST_MPI_NG
      cout <<"assign_bcs: assigning bc "<<i<<" with type "<<b->type<<"\n";
#endif
      err += BC_assign_STWIND(    par,grid,grid->BC_bd[i]);
      break;
    case BCMPI:
    case FINE_TO_COARSE: case COARSE_TO_FINE:
    case FINE_TO_COARSE_SEND: case FINE_TO_COARSE_RECV:
    case COARSE_TO_FINE_SEND: case COARSE_TO_FINE_RECV:
      break; // assigned in NG grid class
    default:
      rep.error("assign_update_bcs::Unhandled BC: assign",
                                            grid->BC_bd[i]->itype);
      break;
    }

  }
  return(err);
}



// ##################################################################
// ##################################################################





int assign_update_bcs::TimeUpdateInternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid,  ///< pointer to grid.
      class FV_solver_base *solver, ///< pointer to equations
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  struct boundary_data *b;
  int err=0;
  for (size_t i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case STWIND:
      err += BC_update_STWIND(par,grid, simtime, b, cstep, maxstep);
      break;
    case PERIODIC:
    case OUTFLOW:
    case ONEWAY_OUT:
    case INFLOW:
    case REFLECTING:
    case AXISYMMETRIC:
    case FIXED:
    case JETBC:
    case JETREFLECT:
    case DMACH:
    case DMACH2:
    case BCMPI:
      //
      // External BCs updated elsewhere
      //     
      break;
      
    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.error("assign_update_bcs::Unhandled BC:  internal",b->itype);
      break;
    }
  }
  return(0);
}
  


// ##################################################################
// ##################################################################



int assign_update_bcs::TimeUpdateExternalBCs(
      class SimParams &par,     ///< simulation parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid,    ///< pointer to grid.
      class FV_solver_base *solver, ///< pointer to equations
      const double simtime,     ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
#ifdef TEST_MPI_NG
  cout <<"update_bcs: external boundary update"<<endl;
#endif
  struct boundary_data *b;
  int err=0;
  for (size_t i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case PERIODIC:
      err += BC_update_PERIODIC(  par,level,grid, b, cstep, maxstep);
      break;
    case OUTFLOW:
      err += BC_update_OUTFLOW(    par,grid, b, cstep, maxstep);
      break;
    case ONEWAY_OUT:
      err += BC_update_ONEWAY_OUT( par,grid, b, cstep, maxstep);
      break;
    case INFLOW:
      err += BC_update_INFLOW(     par,grid, b, cstep, maxstep);
      break;
    case REFLECTING:
      err += BC_update_REFLECTING( par,grid, b, cstep, maxstep);
      break;
    case AXISYMMETRIC:
      err += BC_update_AXISYMMETRIC( par,grid, b, cstep, maxstep);
      break;
    case FIXED:
      err += BC_update_FIXED(      par,grid, b, cstep, maxstep);
      break;
    case JETBC:
      err += BC_update_JETBC(      par,grid, b, cstep, maxstep);
      break;
    case JETREFLECT:
      err += BC_update_JETREFLECT( par,grid, b, cstep, maxstep);
      break;
    case DMACH:
      err += BC_update_DMACH(      par,grid, simtime, b, cstep, maxstep);
      break;
    case DMACH2:
      err += BC_update_DMACH2(     par,grid, b, cstep, maxstep);
      break;
    case RADSHOCK:
    case RADSH2:
    case STWIND:
    case BCMPI:
    case FINE_TO_COARSE:
    case COARSE_TO_FINE:
      //
      // internal bcs updated separately
      //
      break;
    default:
      rep.error("assign_update_bcs::Unhandled BC: external",b->itype);
      break;
    }
  }
  return(0);
}



// ##################################################################
// ##################################################################



