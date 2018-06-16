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

//#define TEST_NEST

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
  //cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
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
      //cout <<"found FINE_TO_COARSE boundary to update\n";
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
  //cout <<"BC_nbd = "<<grid->BC_bd.size()<<"\n";
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
      // only update if at the start of a full step.
      if (cstep==maxstep) {
        //cout <<"found COARSE_TO_FINE boundary to update\n";
        err += BC_update_COARSE_TO_FINE(par, level, b, par.levels[level].step);
      }
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
  //cout <<"Fine to Coarse boundary update (internal).  list sizes: ";
  //cout << b->data.size() <<",  :"<<b->nest.size()<<"\n";

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
      const int step
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
  double U[par.nvar], P[par.nvar];

  //
  // if on an odd-numbered step, then need to update the data on the
  // coarse grid to half way through a coarse step.  Assume dU has
  // been calculated for the coarse grid, but the state vector not
  // updated, so we can convert U+0.5*dU into a new primitive
  // state for the half step.
  // NB: We overwrite Ph[] in the coarse cell, assuming it is not
  // needed anymore on the coarse grid because dU[] is already
  // calculated.
  //
  if (step%2 != 0) {
    for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
      f = (*f_iter);
      c = f->npt;
      eqn->PtoU(c->P, U, par.gamma);
      for (int v=0;v<par.nvar;v++) U[v] += 0.5*c->dU[v];
      eqn->UtoP(U,c->Ph, par.EP.MinTemperature, par.gamma);
    }
  }
  
  if (par.spOOA == OA1) {
    for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
      f = (*f_iter);
      c = f->npt;
      // ----- constant data -----
      for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v];
      for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
      for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;
      // ----- constant data -----
    }
  }

  else if (par.spOOA == OA2) {
    //
    // Each dimension is sufficiently different that we have an if/else
    // loop for each dimension, and then do linear/bilinear/trilinear
    // interpolation as needed.
    //
    if (par.ndim ==1) {
      double sx[par.nvar]; // slope in x-dir
      double dx = fine->DX(); // dx
      //
      // Do two fine cells at a time: they have the same parent.
      // In 1D the geometry is very easy.
      //
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        f = (*f_iter);
        c = f->npt;
        coarse->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
        for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] * (1.0-0.25*dx*sx[v]);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;
        
        f_iter++;
        f = (*f_iter);
        for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] * (1.0+0.25*dx*sx[v]);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;

        // Now need to check mass/momentum/energy conservation between
        // coarse and fine levels!
        rep.error("Check conservation 1D",9);
      } // loop over fine cells
    } // 1D

    else if (par.ndim == 2) {
      //
      // Need to do bilinear interpolation.  Would ideally do 4
      // cells at a time, but for now just do 2 (ordering in the
      // list is not ideal).
      //
      double sx[par.nvar], sy[par.nvar]; // slope in x-dir
      double dxo2 = 0.5*fine->DX(); // dx
      double P00[par.nvar], P10[par.nvar], P01[par.nvar], P11[par.nvar];
      int signx=0, signy=0;
      //
      // Do two fine cells at a time: they have the same parent.
      // First get the values at the 4 corners of coarse-cell c.
      // Then interpolate with bilinear interpolation to the fine-cell
      // positions as 0.25*dx in each direction from the corners.
      //
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        f = (*f_iter);
        c = f->npt;
        coarse->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
        coarse->SetSlope(c,YY,par.nvar,sy,OA2,coarse);
        for (int v=0;v<par.nvar;v++) sx[v] *= dx02;
        for (int v=0;v<par.nvar;v++) sy[v] *= dx02;
        for (int v=0;v<par.nvar;v++) P00[v] = c->Ph[v] * (1.0-sx) * (1.0-sy);
        for (int v=0;v<par.nvar;v++) P10[v] = c->Ph[v] * (1.0+sx) * (1.0-sy);
        for (int v=0;v<par.nvar;v++) P01[v] = c->Ph[v] * (1.0-sx) * (1.0+sy);
        for (int v=0;v<par.nvar;v++) P11[v] = c->Ph[v] * (1.0+sx) * (1.0+sy);

        if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
          // 1/4,1/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (9.0*P00[v] + 3.0*P10[v] + 3.0*P01[v] +     P11[v]);
        }
        else if ( (f->pos[XX] > c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
          // 3/4,1/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (3.0*P00[v] + 9.0*P10[v] +     P01[v] + 3.0*P11[v]);
        }
        else if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] > c->pos[YY]) ) {
          // 1/4,3/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (3.0*P00[v] +     P10[v] + 9.0*P01[v] + 3.0*P11[v]);
        }
        else {
          // 3/4,3/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (    P00[v] + 3.0*P10[v] + 3.0*P01[v] + 9.0*P11[v]);
        }
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;
        
        // second cell in row, should also have c as parent if we always
        // update the fine cells two by two.
        f_iter++;
        f = (*f_iter);
        if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
          // 1/4,1/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (9.0*P00[v] + 3.0*P10[v] + 3.0*P01[v] +     P11[v]);
        }
        else if ( (f->pos[XX] > c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
          // 3/4,1/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (3.0*P00[v] + 9.0*P10[v] +     P01[v] + 3.0*P11[v]);
        }
        else if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] > c->pos[YY]) ) {
          // 1/4,3/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (3.0*P00[v] +     P10[v] + 9.0*P01[v] + 3.0*P11[v]);
        }
        else {
          // 3/4,3/4
          for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] /16.0 *
                        (    P00[v] + 3.0*P10[v] + 3.0*P01[v] + 9.0*P11[v]);
        }
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;
      } // loop over fine cells
    } // 2D
    else if (par.ndim == 3) {
      //
      // Need to do trilinear interpolation.  Would ideally do 8
      // cells at a time, but for now just do 2 (ordering in the
      // list is not ideal).
      //
      rep.error("3D coarse-to-fine interpolation at 2nd order!",3);
    } // 3D
  } // 2nd-order accuracy

  return 0;
}



// ##################################################################
// ##################################################################




