/// \file assign_update_bcs_MPI.cpp
/// \brief Defines a class that inherits all boundary types and 
///   assigns and updated the boundaries.
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "assign_update_bcs_MPI.h"

#include <sstream>
using namespace std;


// ##################################################################
// ##################################################################



int assign_update_bcs_MPI::assign_boundary_data(
      class SimParams &par,      ///< pointer to simulation parameters
      class MCMDcontrol &ppar,    ///< domain decomposition info
      class GridBaseClass *grid ///< pointer to grid.
      )
{
  cout <<ppar.get_myrank()<<" Setting up MPI boundaries...\n";

  int err= assign_update_bcs::assign_boundary_data(par,grid);
  rep.errorTest("assign_update_bcs::assign_boundary_data",err,0);
  //
  // Loop through all boundaries, and assign data to MPI one.
  //
  for (size_t i=0; i<grid->BC_bd.size(); i++) {
    switch (grid->BC_bd[i]->itype) {
    case BCMPI:
      cout <<ppar.get_myrank()<<" assigning MPI boundary in dir "<<i<<"\n";
      err += BC_assign_BCMPI(par,ppar,grid,grid->BC_bd[i],BC_MPItag);
      break;      
    case PERIODIC:
      // periodic was set up wrongly in serial call (1st line of this
      // function, so we delete the boundary that was set up and 
      // reset it here.
      grid->BC_deleteBoundaryData(grid->BC_bd[i]);
      err += BC_assign_PERIODIC(  par,ppar,grid,grid->BC_bd[i]);
      break;
    case OUTFLOW:
    case ONEWAY_OUT:
    case INFLOW:
    case REFLECTING:
    case FIXED:
    case JETBC:
    case JETREFLECT:
    case DMACH:
    case DMACH2:
    case STWIND:
    case FINE_TO_COARSE:
    case COARSE_TO_FINE:
      break; // assigned in nested grid class
    default:
      rep.error("assign_update_bcs_MPI::Unhandled BC: assign",grid->BC_bd[i]->itype);
      break;
    }

  }
  return(err);
}



// ##################################################################
// ##################################################################



int assign_update_bcs_MPI::TimeUpdateExternalBCs(
      class SimParams &par,   ///< pointer to simulation parameters
      class MCMDcontrol &ppar,    ///< domain decomposition info
      class GridBaseClass *grid,  ///< pointer to grid.
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
      
    case PERIODIC:
      err += BC_update_PERIODIC(par,ppar,grid, b, cstep, maxstep);
      break;
    case BCMPI:
      err += BC_update_BCMPI(par, ppar, grid, b, cstep, maxstep,
                                                        BC_MPItag);
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
    case FINE_TO_COARSE:
    case COARSE_TO_FINE:
      //
      // internal bcs updated separately
      //
      break;

    default:
      rep.error("assign_update_bcs_MPI::Unhandled BC: external",b->itype);
      break;
    }

    if (i==XP || (i==YP && par.ndim>=2) || (i==ZP && par.ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      //      cout <<"Barrier synching for dir="<<i<<" with 1=XP,2=YP,3=ZP.\n";
      ostringstream tmp;
      tmp <<"assign_update_bcs_MPI_TimeUpdateExternalBCs_"<<i;
      COMM->barrier(tmp.str());
      tmp.str("");
    }
  }
  return 0;
}



// ##################################################################
// ##################################################################



