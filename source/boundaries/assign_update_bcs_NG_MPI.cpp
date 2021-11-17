/// \file assign_update_bcs_NG_MPI.cpp
/// \brief Defines a class that inherits boundary types related to
///   static mesh-refinement (NG) and sub_domain domain decomposition and
///   implements assignment and update functions.
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.28 JM: adapted from assign_update_bcs_NG.cpp

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "assign_update_bcs_NG_MPI.h"
#include "tools/reporting.h"
using namespace std;

//#define TEST_MPI_NG
// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::assign_boundary_data(
    class SimParams &par,        ///< simulation parameters
    const int level,             ///< level in grid hierarchy
    class GridBaseClass *grid,   ///< pointer to grid.
    class microphysics_base *mp  ///< pointer to microphysics
)
{
  int err = 0;
  struct boundary_data *b;

  //
  // assign data for each boundary.
  //
  for (size_t i = 0; i < grid->BC_bd.size(); i++) {
#ifdef TEST_MPI_NG
    cout << "LEVEL " << level << ": NG grid assign BCs: BC[" << i
         << "] starting.\n";
    cout.flush();
#endif
    b = grid->BC_bd[i];
    switch (b->itype) {
      case PERIODIC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_PERIODIC(par, level, grid, b);
        break;
      case BCMPI:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_BCMPI(par, level, grid, grid->BC_bd[i], BC_MPItag);
        break;
      case OUTFLOW:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_OUTFLOW(par, grid, b);
        break;
      case ONEWAY_OUT:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_ONEWAY_OUT(par, grid, b);
        break;
      case INFLOW:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_INFLOW(par, grid, b);
        break;
      case REFLECTING:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_REFLECTING(par, grid, b);
        break;
      case AXISYMMETRIC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_AXISYMMETRIC(par, grid, b);
        break;
      case FIXED:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_FIXED(par, grid, b);
        break;
      case JETBC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_JETBC(par, grid, b);
        break;
      case JETREFLECT:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_JETREFLECT(par, grid, b);
        break;
      case DMACH:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_DMACH(par, grid, b);
        break;
      case DMACH2:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_DMACH2(par, grid, b);
        break;
      case STWIND:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
#endif
        err += BC_assign_STWIND(par, grid, b, mp);
        break;
      case FINE_TO_COARSE_SEND:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_FINE_TO_COARSE_SEND(par, level, b);
        break;

      case FINE_TO_COARSE_RECV:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_FINE_TO_COARSE_RECV(par, level, b);
        break;

      case COARSE_TO_FINE_SEND:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_COARSE_TO_FINE_SEND(par, level, b);
        break;

      case COARSE_TO_FINE_RECV:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": NG_MPI_Assign: assigning bc " << i;
        cout << " with type " << b->type << "\n";
        cout.flush();
#endif
        err += BC_assign_COARSE_TO_FINE_RECV(par, level, b);
        break;

      default:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": leaving BC " << i
             << " alone in NG grid assign fn.\n";
#endif
        rep.error("Unhandled boundary", 2);
        break;
    }
  }
  return err;
}



// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::TimeUpdateExternalBCs(
    class SimParams &par,          ///< pointer to sim parameters
    const int level,               ///< level in grid hierarchy
    class GridBaseClass *grid,     ///< pointer to grid.
    class FV_solver_base *solver,  ///< pointer to equations
    const double simtime,          ///< current simulation time
    const int cstep,
    const int maxstep)
{
#ifdef TEST_MPI_NG
  cout << "LEVEL " << level << ": update_bcs_NG_MPI: ext BC update" << endl;
#endif

  struct boundary_data *b;
  int err  = 0;
  size_t i = 0;

  // make a map so that boundaries are updated pairwise.
  // ix[] is the position of this process in each direction of the
  // block of MPI processes.
  int nb = grid->BC_bd.size();
  int map[nb];
  int rank = par.levels[level].sub_domain.get_myrank();
  int ix[MAX_DIM];
  par.levels[level].sub_domain.get_domain_coordinates(rank, ix);
  for (int j = 0; j < par.ndim; j++) {
    if (ix[j] % 2 == 0) {
      map[2 * j]     = 2 * j;
      map[2 * j + 1] = 2 * j + 1;
    }
    else {
      map[2 * j]     = 2 * j + 1;
      map[2 * j + 1] = 2 * j;
    }
  }
  for (int j = 2 * par.ndim; j < nb; j++) {
    map[j] = j;
  }

#ifdef TEST_MPI_NG
  rep.printVec("TimeUpdateExternalBCs Map", map, nb);
  cout << "LEVEL " << level << ": BC_nbd = " << grid->BC_bd.size() << "\n";
#endif
  for (i = 0; i < grid->BC_bd.size(); i++) {
    b = grid->BC_bd[map[i]];
#ifdef TEST_MPI_NG
    cout << "updating bc " << map[i] << " with type " << b->type << "\n";
    cout.flush();
#endif
    switch (b->itype) {
      case PERIODIC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_PERIODIC(par, level, grid, b, cstep, maxstep);
        break;
      case BCMPI:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_BCMPI(par, level, grid, b, cstep, maxstep, BC_MPItag);
        break;
      case OUTFLOW:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_OUTFLOW(par, grid, b, cstep, maxstep);
        break;
      case ONEWAY_OUT:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_ONEWAY_OUT(par, grid, b, cstep, maxstep);
        break;
      case INFLOW:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_INFLOW(par, grid, b, cstep, maxstep);
        break;
      case REFLECTING:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_REFLECTING(par, grid, b, cstep, maxstep);
        break;
      case AXISYMMETRIC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_AXISYMMETRIC(par, grid, b, cstep, maxstep);
        break;
      case FIXED:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_FIXED(par, grid, b, cstep, maxstep);
        break;
      case JETBC:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_JETBC(par, grid, b, cstep, maxstep);
        break;
      case JETREFLECT:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_JETREFLECT(par, grid, b, cstep, maxstep);
        break;
      case DMACH:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_DMACH(par, grid, simtime, b, cstep, maxstep);
        break;
      case DMACH2:
#ifdef TEST_MPI_NG
        cout << "LEVEL " << level << ": update_bcs_NG_MPI: updating bc ";
        cout << map[i] << " with type " << b->type << "\n";
#endif
        err += BC_update_DMACH2(par, grid, b, cstep, maxstep);
        break;
      case STWIND:
      case FINE_TO_COARSE_SEND:
      case FINE_TO_COARSE_RECV:
#ifdef TEST_MPI_NG
        // cout <<"skipping this boundary in MPI_NG \n";
#endif
        break;

      case COARSE_TO_FINE_RECV:
      case COARSE_TO_FINE_SEND:
        break;

      default:
        rep.error("Unhandled BC: NG-MPI update external", b->itype);
        break;
    }
    par.levels[0].sub_domain.barrier("external bcs");
  }
#ifdef TEST_MPI_NG
  cout << "LEVEL " << level << ": updated NG-MPI-grid external BCs\n";
#endif

  return (0);
}

// ##################################################################
// ##################################################################



int assign_update_bcs_NG_MPI::TimeUpdateInternalBCs(
    class SimParams &par,          ///< pointer to simulation parameters
    const int level,               ///< level in grid hierarchy
    class GridBaseClass *grid,     ///< pointer to grid.
    class FV_solver_base *solver,  ///< pointer to equations
    const double simtime,          ///< current simulation time
    const int cstep,
    const int maxstep)
{
#ifdef TEST_MPI_NG
  cout << "assign_update_bcs_NG::TimeUpdateInternalBCs() running.\n";
#endif
  struct boundary_data *b;
  int err  = 0;
  size_t i = 0;
  // cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i = 0; i < grid->BC_bd.size(); i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
      case STWIND:
        err += BC_update_STWIND(par, grid, simtime, b, cstep, maxstep);
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
      case COARSE_TO_FINE:
      case FINE_TO_COARSE:
      case COARSE_TO_FINE_SEND:
      case COARSE_TO_FINE_RECV:
      case FINE_TO_COARSE_SEND:
      case FINE_TO_COARSE_RECV:
        //
        // boundaries not affected by NG grid are updated elsewhere
        //
        break;
    }
  }
#ifdef TEST_MPI_NG
  cout << "updated NG-grid serial internal BCs\n";
  cout << "assign_update_bcs_NG::TimeUpdateInternalBCs() returns.\n";
#endif
  return 0;
}

// ##################################################################
// ##################################################################
