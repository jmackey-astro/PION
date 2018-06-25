/// \file update_boundaries.cpp
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
#include "sim_control/update_boundaries.h"


// ##################################################################
// ##################################################################

update_boundaries::update_boundaries()
{
  spatial_solver=0;
  return;
}

update_boundaries::~update_boundaries()
{
  if (spatial_solver) {delete spatial_solver; spatial_solver=0;}
  return;
}



int update_boundaries::TimeUpdateInternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case STWIND:     err += BC_update_STWIND(  par,grid, simtime, b, cstep, maxstep); break;
    case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
    case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2: case BCMPI:
    case FINE_TO_COARSE: case COARSE_TO_FINE:
      //
      // External BCs updated elsewhere
      //     
      break;
      
    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.error("Unhandled BC: serial update internal",b->itype);
      break;
    }
  }
  return(0);
}
  


// ##################################################################
// ##################################################################



int update_boundaries::TimeUpdateExternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  // TEMP_FIX
  //BC_printBCdata(&BC_bd[0]);
  // TEMP_FIX
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
    case PERIODIC:   err += BC_update_PERIODIC(   par,grid, b, cstep, maxstep); break;
    case OUTFLOW:    err += BC_update_OUTFLOW(    par,grid, b, cstep, maxstep); break;
    case ONEWAY_OUT: err += BC_update_ONEWAY_OUT( par,grid, b, cstep, maxstep); break;
    case INFLOW:     err += BC_update_INFLOW(     par,grid, b, cstep, maxstep); break;
    case REFLECTING: err += BC_update_REFLECTING( par,grid, b, cstep, maxstep); break;
    case FIXED:      err += BC_update_FIXED(      par,grid, b, cstep, maxstep); break;
    case JETBC:      err += BC_update_JETBC(      par,grid, b, cstep, maxstep); break;
    case JETREFLECT: err += BC_update_JETREFLECT( par,grid, b, cstep, maxstep); break;
    case DMACH:      err += BC_update_DMACH(      par,grid, simtime, b, cstep, maxstep); break;
    case DMACH2:     err += BC_update_DMACH2(     par,grid, b, cstep, maxstep); break;
    case RADSHOCK: case RADSH2: case STWIND: case BCMPI: case FINE_TO_COARSE: case COARSE_TO_FINE:
      //
      // internal bcs updated separately
      //
      break;
    default:
      //      cout <<"do i have a bc to update?? no.\n";
      rep.error("Unhandled BC: serial update external",b->itype);
      break;
    }
  }
  return(0);
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_PERIODIC(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = (*c)->npt->Ph[v];
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) (*c)->P[v] = (*c)->npt->P[v];
    }
  } // all cells.
  return 0;
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_OUTFLOW(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // Outflow or Absorbing BCs; boundary cells are same as edge cells.
  // This is zeroth order outflow bcs.
  //
  list<cell*>::iterator c=b->data.begin();
  cell *gc;

  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //
    // gc is the on-grid cell.
    //
    gc = (*c)->npt;
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = gc->Ph[v];
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) (*c)->P[v] = gc->P[v];
    }
    
    //
    // The GLM boundary is somewhat different, because I found
    // that zero-gradient didn't work well (Mackey & Lim, 2011,
    // MNRAS,412,2079).  So we switch the sign instead.
    //
#ifdef GLM_ZERO_BOUNDARY
    if (par.eqntype==EQGLM) {
      (*c)->P[SI]  = 0.0;
      (*c)->Ph[SI] = 0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (par.eqntype==EQGLM) {

      if      ((*c)->isedge == -1) {
        (*c)->P[SI]  = -gc->P[SI];
        (*c)->Ph[SI] = -gc->Ph[SI];
      }
      else if ((*c)->isedge == -2) {
        (*c)->P[SI]  = -NextPt(gc,b->ondir)->P[SI];
        (*c)->Ph[SI] = -NextPt(gc,b->ondir)->Ph[SI];
      }
      else {
        rep.error("only know 1st/2nd order bcs",(*c)->id);
      }
    }
#endif // GLM_NEGATIVE_BOUNDARY

  } // all cells.
  return 0;
}
// ##################################################################
// ##################################################################



int update_boundaries::BC_update_ONEWAY_OUT(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // Outflow or Absorbing BCs; boundary cells are same as edge cells.
  // This is zeroth order outflow bcs.
  //
  // For the one-way boundary we set the normal velocity to be zero if
  // the flow wants to be onto the domain.
  //

  //
  // First get the normal velocity component, and whether the offdir is
  // positive or negative.
  //
  enum direction offdir = b->dir;
  enum primitive Vnorm = RO;
  int norm_sign = 0;
  switch (offdir) {
  case XN: case XP:
    Vnorm = VX; break;
  case YN: case YP:
    Vnorm = VY; break;
  case ZN: case ZP:
    Vnorm = VZ; break;
  default:
    rep.error("bad dir",offdir);
  }
  if (Vnorm == RO)
    rep.error("Failed to set normal velocity component",Vnorm);

  switch (offdir) {
  case XN: case YN: case ZN:
    norm_sign = -1; break;
  case XP: case YP: case ZP:
    norm_sign =  1; break;
  default:
    rep.error("bad dir",offdir);
  }

  //
  // Now run through all cells in the boundary
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //
    // exactly same routine as for periodic.
    // ONEWAY_OUT: overwrite the normal velocity if it is inflow:
    //
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = (*c)->npt->Ph[v];
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.;
    (*c)->Ph[Vnorm] = norm_sign*max(static_cast<pion_flt>(0.0),
                                    (*c)->Ph[Vnorm]*norm_sign);
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) (*c)->P[v] = (*c)->npt->P[v];
      (*c)->P[Vnorm] = norm_sign*max(static_cast<pion_flt>(0.0),
                                      (*c)->P[Vnorm]*norm_sign);
    }
    
    //
    // The GLM boundary is somewhat different, because I found
    // that zero-gradient didn't work well (Mackey & Lim, 2011,
    // MNRAS,412,2079).  So we switch the sign instead.  Same as
    // periodic BCs.
    //
#ifdef GLM_ZERO_BOUNDARY
    if (par.eqntype==EQGLM) {
      (*c)->P[SI]  = 0.0;
      (*c)->Ph[SI] = 0.0;
    }
#endif // GLM_ZERO_BOUNDARY
#ifdef GLM_NEGATIVE_BOUNDARY
    if (par.eqntype==EQGLM) {

      if      ((*c)->isedge == -1) {
        (*c)->P[SI]  = -(*c)->npt->P[SI];
        (*c)->Ph[SI] = -(*c)->npt->Ph[SI];
      }
      else if ((*c)->isedge == -2) {
        (*c)->P[SI]  = -NextPt(((*c)->npt),b->ondir)->P[SI];
        (*c)->Ph[SI] = -NextPt(((*c)->npt),b->ondir)->Ph[SI];
      }
      else {
        rep.error("only know 1st/2nd order bcs",(*c)->id);
      }
    }
#endif // GLM_NEGATIVE_BOUNDARY


  } // all cells.
  return 0;
  
}
// ##################################################################
// ##################################################################



int update_boundaries::BC_update_INFLOW(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int , /// current ooa step e.g. 1 (not used here)
      const int   /// overall ooa      e.g. 2 (not used here)
      )
{
  //
  // Inflow means BC is constant at the initial value of the
  // neighbouring edge cell, so we set dU=0.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.0;
  } // all cells.
  return 0;
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_REFLECTING(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // same routine as for outflow, except multiply v_n,B_n by -1.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v] *b->refval[v];
    }
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) {
        (*c)->P[v] = (*c)->npt->P[v] *b->refval[v];
      }
    }
  } // all cells.
  return 0;
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_FIXED(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int , // current ooa step e.g. 1 (unused here)
      const int   // overall ooa      e.g. 2 (unused here)
      )
{
  //
  // Fixed means all boundary points have same fixed value, stored
  // in refval.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) (*c)->dU[v]=0.;
    for (int v=0;v<par.nvar;v++) (*c)->P[v]  = b->refval[v];
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = b->refval[v];
  }    // all cells.   
  return 0;
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_JETBC(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int,
      const int
      )
{
# ifdef SOFTJET
  double dist=0.0;
  double jr = JP.jetradius*grid->DX();
# endif //SOFTJET

  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) (*c)->dU[v]=0.0;
    for (int v=0;v<par.nvar;v++) (*c)->P[v]  = b->refval[v];
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = b->refval[v];
# ifdef SOFTJET
    dist =0.0;
    if      (par.ndim==2) dist = CI.get_dpos(*c,YY);
    else if (par.ndim==3) {
      dist = sqrt(CI.get_dpos(*c,YY)*CI.get_dpos(*c,YY)+
                  CI.get_dpos(*c,ZZ)*CI.get_dpos(*c,ZZ) );
    }
    else rep.error("Jet BC, but not 2d or 3d!!!",par.ndim);
    (*c)->P[VX]  = b->refval[VX] *min(1., 4-4.0*dist/jr);
    (*c)->Ph[VX] = (*c)->P[VX];
# endif //SOFTJET
  }    // all cells.   
  return 0;
}



// ##################################################################
// ##################################################################



int update_boundaries::BC_update_JETREFLECT(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // same routine as for reflecting, except the normal B is
  // unchanged, but the tangential is reversed.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) {
      (*c)->Ph[v] = (*c)->npt->Ph[v] *b->refval[v];
    }
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.;
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) {
        (*c)->P[v] = (*c)->npt->P[v] *b->refval[v];
      }
    }
  } // all cells.
  return 0;
}

// ##################################################################
// ##################################################################



int update_boundaries::BC_update_DMACH(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  double bpos = 0.0;
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    //
    // This is the boundary position:
    //
    bpos =  10.0*simtime/sin(M_PI/3.0)
          + 1.0/6.0
          + CI.get_dpos(*c,YY)/tan(M_PI/3.0);

    if (CI.get_dpos(*c,XX) <= bpos) {
      (*c)->Ph[RO] = 8.0;
      (*c)->Ph[PG] = 116.5;
      (*c)->Ph[VX] = 7.14470958;
      (*c)->Ph[VY] = -4.125;
      (*c)->Ph[VZ] = 0.0;
      for (int v=par.ftr; v<par.nvar; v++) (*c)->Ph[v] = 1.0;
    }
    else {
      for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = b->refval[v];
    }
    for (int v=0;v<par.nvar;v++) (*c)->dU[v] = 0.0;
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) (*c)->P[v] = (*c)->Ph[v];
    }
  }
  return 0;
}



// ##################################################################
// ##################################################################



int update_boundaries::BC_update_DMACH2(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      struct boundary_data *b,
      const int,
      const int
      )
{
  //
  // Fixed at all times, so no difference between full and half step.
  //
  list<cell*>::iterator c=b->data.begin();
  for (c=b->data.begin(); c!=b->data.end(); ++c) {
    for (int v=0;v<par.nvar;v++) (*c)->dU[v]=0.;
    for (int v=0;v<par.nvar;v++) (*c)->P[v]  = b->refval[v];
    for (int v=0;v<par.nvar;v++) (*c)->Ph[v] = b->refval[v];
  }    // all cells.   
  return 0;
}



// ##################################################################
// ##################################################################



//
// Update internal stellar wind boundaries -- these are (possibly
// time-varying) winds defined by a mass-loss-rate and a terminal
// velocity.  If fixed in time the wind is updated with b->refval,
// otherwise with a (slower) call to the  stellar wind class SW
//
int update_boundaries::BC_update_STWIND(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const double simtime,   ///< current simulation time
      boundary_data *b, ///< Boundary to update.
      const int ,  ///< current fractional step being taken.
      const int    ///< final step (not needed b/c fixed BC).
      )
{
  //
  // The stellar_wind class already has a list of cells to update
  // for each source, together with pre-calculated state vectors,
  // so we just call the set_cell_values() function.
  //
#ifdef TESTING
  cout <<"update_boundaries: updating wind boundary\n";
#endif
  int err=0;
  for (int id=0;id<grid->Wind->Nsources();id++) {
#ifdef TESTING
    cout <<"update_boundaries: updating wind boundary for id="<<id<<"\n";
#endif
    err += grid->Wind->set_cell_values(grid, id,simtime);
  }

  return err;
}




