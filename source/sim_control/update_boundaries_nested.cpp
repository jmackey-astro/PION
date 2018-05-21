/// \file update_boundaries_nested.cpp
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class definitions for routines to update grid boundaries with
/// different boundary conditions.
///
/// Modifications:\n
/// - 2018.05.10 JM: moved code from uniform_grid.cpp
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "sim_control/update_boundaries_nested.h"


// ##################################################################
// ##################################################################

update_boundaries_nested::update_boundaries_nested()
{
  return;
}



// ##################################################################
// ##################################################################

update_boundaries_nested::~update_boundaries_nested()
{
  return;
}



// ##################################################################
// ##################################################################



int update_boundaries_nested::TimeUpdateInternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{

  int err = update_boundaries::TimeUpdateInternalBCs(par,grid,simtime,cstep,maxstep);
  rep.errorTest("update_boundaries_nested: uni-grid int. BC update",0,err);

  struct boundary_data *b;
  int i=0;
  for (i=0;i<BC_nbd;i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case STWIND:
    case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
    case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2: case BCMPI:
      //
      // boundaries not affected by nested grid are updated elsewhere
      //     
      break;

    case NEST_FINE:
      err += BC_update_NEST_FINE( b, cstep, maxstep); break;
    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.warning("Unhandled BC: serial update internal",b->itype,-1); err+=1; break;
      break;
    }
  }
  return(0);
}
  


// ##################################################################
// ##################################################################



int update_boundaries_nested::TimeUpdateExternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  int err = update_boundaries::TimeUpdateExternalBCs(par,grid,simtime,cstep,maxstep);
  rep.errorTest("update_boundaries_nested: uni-grid ext. BC update",0,err);

  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<BC_nbd;i++) {
    b = grid->BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
      // skip all these:
      case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
      case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2:
      case STWIND: case BCMPI:
      break;

      // get outer boundary of this grid from coarser grid.
      case NEST_COARSE:
      err += BC_update_NEST_FINE( b, cstep, maxstep); break;

      default:
      rep.warning("Unhandled BC: serial update external",b->itype,-1); err+=1; break;
      break;
    }
  }
  return(0);
}



// ##################################################################
// ##################################################################



int update_boundaries_nested::BC_update_NEST_FINE(
      class SimParams &par,      ///< pointer to simulation parameters
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
  double cd[par.nvar], u[par.nvar], vol=0.0;

  for (coarse=b->data.begin(); coarse!=b->data.end(); ++coarse) {
    c = (*coarse);
    f = (*fine);

    // 1D
    // get conserved vars for cell 1 in fine grid, *cellvol.
    eqn->PtoU(
    for (int v=0;v<par.nvar;v++) {
    // get conserved vars for cell 2 in fine grid, *cellvol.

    // if 2D
    // get conserved vars for cell 3 in fine grid, *cellvol.
    // get conserved vars for cell 4 in fine grid, *cellvol.
    
    // if 3D
    // get conserved vars for cell 5 in fine grid, *cellvol.
    // get conserved vars for cell 6 in fine grid, *cellvol.
    // get conserved vars for cell 7 in fine grid, *cellvol.
    // get conserved vars for cell 8 in fine grid, *cellvol.


      


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
int update_boundaries_nested::BC_update_NEST_COARSE(
      class SimParams &par,      ///< pointer to simulation parameters
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




