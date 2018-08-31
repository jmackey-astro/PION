/// \file assign_update_bcs_NG_MPI.cpp
/// \brief Defines a class that inherits boundary types related to
///   static mesh-refinement (NG) and MCMD domain decomposition and
///   implements assignment and update functions.
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.28 JM: adapted from assign_update_bcs_NG.cpp

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "assign_update_bcs_NG_MPI.h"
using namespace std;


// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::assign_boundary_data(
      class SimParams &par,     ///< pointer to simulation parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid ///< pointer to grid.
      )
{
  // first call the Uniform Grid version.
  int err = assign_update_bcs::assign_boundary_data(par,level,grid);
  rep.errorTest("assign_update_bcs::assign_boundary_data",err,0);

  //
  // Then check for NG-grid boundaries and assign data for them.
  //
  for (size_t i=0; i<grid->BC_bd.size(); i++) {
#ifdef TESTING
    cout <<"NG grid assign BCs: BC["<<i<<"] starting.\n";
#endif
    switch (grid->BC_bd[i]->itype) {
      case FINE_TO_COARSE:
#ifdef TESTING
      cout <<"NG grid setup: Assigning FINE_TO_COARSE BC\n";
#endif
      err += BC_assign_FINE_TO_COARSE(par,grid,grid->BC_bd[i],child);
      break;

      case COARSE_TO_FINE:
#ifdef TESTING
      cout <<"assign_update_bcs_NG_MPI:: Assigning COARSE_TO_FINE BC\n";
#endif
      err += BC_assign_COARSE_TO_FINE(par,grid,grid->BC_bd[i],parent);
      break;

      default:
#ifdef TESTING
      cout <<"leaving BC "<<i<<" alone in NG grid assign fn.\n";
#endif
      break;
    }
  }
  return err;
}



// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::TimeUpdateInternalBCs(
      class SimParams &par,     ///< pointer to sim parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid,    ///< pointer to grid.
      class FV_solver_base *solver, ///< pointer to equations
      const double simtime,     ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  int err = assign_update_bcs::TimeUpdateInternalBCs(par,level,grid,
                                      solver,simtime,cstep,maxstep);
  rep.errorTest("assign_update_bcs_NG_MPI: uni-grid int. BC update",
                                                            0,err);
#ifdef TEST_NEST
  cout <<"updated unigrid serial internal BCs\n";
#endif

  struct boundary_data *b;
  size_t i=0;
  //cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case STWIND:
    case PERIODIC:
    case OUTFLOW:
    case ONEWAY_OUT:
    case INFLOW:
    case REFLECTING:
    case FIXED:
    case JETBC:
    case JETREFLECT:
    case DMACH:
    case DMACH2:
    case BCMPI:
    case COARSE_TO_FINE:
      //
      // boundaries not affected by NG grid are updated elsewhere
      //     
      break;

    case FINE_TO_COARSE:
      //cout <<"found FINE_TO_COARSE boundary to update\n";
      err += BC_update_FINE_TO_COARSE(par,solver,level,b,cstep,maxstep);
      break;

    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.error("Unhandled BC: serial NG update internal",b->itype);
      break;
    }
  }
#ifdef TEST_NEST
  cout <<"updated NG-grid serial internal BCs\n";
#endif
  return 0;
}
  


// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::TimeUpdateExternalBCs(
      class SimParams &par,     ///< pointer to sim parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid,    ///< pointer to grid.
      class FV_solver_base *solver, ///< pointer to equations
      const double simtime,     ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  int err = assign_update_bcs::TimeUpdateExternalBCs(par,level,grid,
                                      solver,simtime,cstep,maxstep);
  rep.errorTest("assign_update_bcs_NG_MPI: uni-grid ext. BC update",
                                                            0,err);
#ifdef TEST_NEST
  cout <<"updated unigrid serial external BCs\n";
#endif

  struct boundary_data *b;
  size_t i=0;
  //cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
      // skip all these:
      case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW:
      case REFLECTING: case FIXED: case JETBC: case JETREFLECT:
      case DMACH: case DMACH2:
      case STWIND: case BCMPI: case FINE_TO_COARSE:
      break;

      // get outer boundary of this grid from coarser grid.
      case COARSE_TO_FINE:
      // only update if at the start of a full step.
      if (cstep==maxstep) {
        //cout <<"found COARSE_TO_FINE boundary to update\n";
        err += BC_update_COARSE_TO_FINE(par,solver,level,b,
                                        par.levels[level].step);
      }
      break;

      default:
      rep.error("Unhandled BC: NG-MPI update external",b->itype);
      break;
    }
  }
#ifdef TEST_NEST
  cout <<"updated NG-MPI-grid serial external BCs\n";
#endif
  return(0);
}



// ##################################################################
// ##################################################################



