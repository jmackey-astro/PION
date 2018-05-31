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

#define TEST_NEST

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
      const int level, ///< level in the nested grid structure
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  int err = update_boundaries::TimeUpdateInternalBCs(par,grid,simtime,cstep,maxstep);
  rep.errorTest("update_boundaries_nested: uni-grid int. BC update",0,err);
#ifdef TEST_NEST
  cout <<"updated unigrid serial internal BCs\n";
#endif

  struct boundary_data *b;
  int i=0;
  cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    switch (b->itype) {
    case STWIND:
    case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
    case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2: case BCMPI:
    case COARSE_TO_FINE:
      //
      // boundaries not affected by nested grid are updated elsewhere
      //     
      break;

    case FINE_TO_COARSE:
      cout <<"found FINE_TO_COARSE boundary to update\n";
      err += BC_update_FINE_TO_COARSE(par, level, b, cstep, maxstep);
      break;

    default:
      //      cout <<"no internal boundaries to update.\n";
      rep.error("Unhandled BC: serial nested update internal",b->itype);
      break;
    }
  }
#ifdef TEST_NEST
  cout <<"updated nested-grid serial internal BCs\n";
#endif
  return(0);
}
  


// ##################################################################
// ##################################################################



int update_boundaries_nested::TimeUpdateExternalBCs(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      const int level, ///< level in the nested grid structure
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  int err = update_boundaries::TimeUpdateExternalBCs(par,grid,simtime,cstep,maxstep);
  rep.errorTest("update_boundaries_nested: uni-grid ext. BC update",0,err);
#ifdef TEST_NEST
  cout <<"updated unigrid serial external BCs\n";
#endif

  struct boundary_data *b;
  int i=0;
  cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
  for (i=0;i<grid->BC_bd.size();i++) {
    b = grid->BC_bd[i];
    //    cout <<"updating bc "<<i<<" with type "<<b->type<<"\n";
    switch (b->itype) {
      // skip all these:
      case PERIODIC: case OUTFLOW: case ONEWAY_OUT: case INFLOW: case REFLECTING:
      case FIXED: case JETBC: case JETREFLECT: case DMACH: case DMACH2:
      case STWIND: case BCMPI: case FINE_TO_COARSE:
      break;

      // get outer boundary of this grid from coarser grid.
      case COARSE_TO_FINE:
      cout <<"found COARSE_TO_FINE boundary to update\n";
      err += BC_update_COARSE_TO_FINE(par, level, b, cstep, maxstep);
      break;

      default:
      rep.error("Unhandled BC: serial nested update external",b->itype);
      break;
    }
  }
#ifdef TEST_NEST
  cout <<"updated nested-grid serial external BCs\n";
#endif
  return(0);
}



// ##################################################################
// ##################################################################



int update_boundaries_nested::BC_update_FINE_TO_COARSE(
      class SimParams &par,      ///< pointer to simulation parameters
      const int level, ///< level in the nested grid structure
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
  list<cell*>::iterator c_iter=b->data.begin();
  list<cell*>::iterator f_iter=b->nest.begin();
  cell *c, *f;
  cout <<"Fine to Coarse boundary update (internal).  list sizes: ";
  cout << b->data.size() <<",  :"<<b->nest.size()<<"\n";

  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].grid;
  class GridBaseClass *fine   = par.levels[level].child;

  // vol_sum is for testing only (make sure that fine grid cells
  // have the same cumulative volume as the coarse cell).
  double cd[par.nvar], u[par.nvar], vol=0.0, vol_sum=0.0;
  int cpos[MAX_DIM],fpos[MAX_DIM];

  for (c_iter=b->data.begin(); c_iter!=b->data.end(); ++c_iter) {
    c = (*c_iter);
    f = (*f_iter);
    vol=0.0, vol_sum=0.0;

#ifdef TEST_NEST
    //CI.get_ipos(c,cpos);
    //CI.get_ipos(f,fpos);
    //rep.printVec("coarse pos:",cpos,par.ndim);
    //rep.printVec("fine 1 pos:",fpos,par.ndim);
#endif
    
    // 1D
    // get conserved vars for cell 1 in fine grid, multiply by cell
    // volume.
    //
    eqn->PtoU(f->Ph, u, par.gamma);
    vol = fine->CellVolume(f);
    for (int v=0;v<par.nvar;v++) cd[v] = u[v]*vol;
#ifdef TEST_NEST
    //vol_sum += vol;
#endif

    // get conserved vars for cell 2 in fine grid, *cellvol.
    f = fine->NextPt(f,XP);
    eqn->PtoU(f->Ph, u, par.gamma);
    vol = fine->CellVolume(f);
    for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
    //CI.get_ipos(f,fpos);
    //rep.printVec("fine 2 pos:",fpos,par.ndim);
    //vol_sum += vol;
#endif

    // if 2D
    if (par.ndim>1) {
      // get conserved vars for cell 3 in fine grid, *cellvol.
      f =  fine->NextPt((*f_iter),YP);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 3 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 4 in fine grid, *cellvol.
      f = fine->NextPt(f,XP);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 4 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif
    }
    
    
    // if 3D
    if (par.ndim>2) {
      // get conserved vars for cell 5 in fine grid, *cellvol.
      f =  fine->NextPt((*f_iter),ZP);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 5 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 6 in fine grid, *cellvol.
      f = fine->NextPt(f,XP);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 6 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 7 in fine grid, *cellvol.
      f = fine->NextPt(f,YP);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 7 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 8 in fine grid, *cellvol.
      f = fine->NextPt(f,XN);
      eqn->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 8 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif
    }

    //
    // divide by coarse cell volume.
    // convert conserved averaged data to primitive
    //
    vol = coarse->CellVolume(c);
#ifdef TEST_NEST
    //cout <<"coarse vol="<<vol<<", fine vol="<<vol_sum<<"\n";
#endif
    for (int v=0;v<par.nvar;v++) cd[v] /= vol;
    eqn->UtoP(cd,c->Ph,par.EP.MinTemperature,par.gamma);
    //
    // if full step then assign to c->P as well as c->Ph.
    //
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) c->P[v] = c->Ph[v];
    }

#ifdef TEST_NEST
//    rep.printVec("coarse data",c->Ph,par.nvar);
//    rep.printVec("fine data",f->Ph,par.nvar);
#endif

    ++f_iter;
  }
  return 0;
}



// ##################################################################
// ##################################################################



/// Updates data to an external boundary from coarser grid.
int update_boundaries_nested::BC_update_COARSE_TO_FINE(
      class SimParams &par,      ///< pointer to simulation parameters
      const int level, ///< level in the nested grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // This is a complicated problem to use linear interpolation (or
  // higher order), because you have to conserve mass, momentum and
  // energy between the two levels.  It should also be monotonic and
  // produce cell-averaged values. Flash and pluto seem
  // to use very old Fortran code to accomplish this, and the code is
  // almost impenetrable.  AMRVAC uses linear/bilinear/trilinear
  // interpolation with renormalisation to conserve mass/P/E.
  //
  // For now we just do contant data.  Will fix this eventually once
  // the rest of the code is working.  It actually doesn't really
  // matter for expanding nebulae because there is supersonic outflow
  // but this needs to be improved for solving more general problems.
  //
  // piecewise-linear data can be done using slopes of primitive 
  // variables in the cell, using bi/tri/linear interpolation

  // number of fine-grid cells for each coarse grid cell.
  int nfine = 1;
  for (int v=0;v<par.ndim;v++) nfine *=2;


  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].parent;
  class GridBaseClass *fine   = par.levels[level].grid;
  // pointers to lists of cells in coarse and fine grids that are
  // part of this boundary.
  //list<cell*>::iterator c_iter=b->nest.begin();
  list<cell*>::iterator f_iter=b->data.begin();
  cell *c, *f;

  for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
    f = (*f_iter);
    c = f->npt;

    // set slope along x-direction in parent cell.
    // parent->SetSlope((*c)->npt,XX,G_nvar,sx,2,parent);
    // get physical offset distance between cell and parent.
    //dist = idifference_cell2cell(*c,(*c)->npt,XX)*CI.phys_per_int();
    // interpolate linearly to cell position.
    //for (int v=0;v<G_nvar;v++)
    //  (*c)->Ph[v] = (*c)->npt->Ph[v] + sx[v]*dist;

    // constant data:
    for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v];
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
    }
    
    // iterate to next coarse cell.
    //++c;
  }
  return 0;
}



// ##################################################################
// ##################################################################




