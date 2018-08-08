/// \file SMR_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for SMR_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/SMR_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################



int SMR_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,  ///< boundary data
      class GridBaseClass *parent  ///< pointer to parent grid.
      )
{
  //
  // Make a list of pointers to cells in the coarser grid that map
  // onto this (finer) grid external boundary, and then write an
  // alogrithm to interpolate the coarse data onto the finer grid.
  //
  if (b->data.empty())
    rep.error("BC_assign_COARSE_TO_FINE: empty boundary data",b->itype);
  b->nest.clear();

  list<cell*>::iterator bpt=b->data.begin();
  //int pidx = parent->idx();
  int gidx = grid->idx();
  //cout <<"BC_assign_COARSE_TO_FINE: dx="<<G_idx<<", parent dx="<<pidx<<"\n";

  cell *pc = parent->FirstPt_All(); // parent cell.

  double distance =  0.0;
  bool loop;
  
  do{
    loop = false;
    distance = grid->idistance(pc->pos, (*bpt)->pos);
    // Find parent cell that covers this boundary cell.  It should be
    // G_idx/2 away from the boundary cell in each direction.
    //rep.printVec("bpt pos",(*bpt)->pos,G_ndim);
    while (distance > gidx && pc!=0) {
      //cout <<"distance="<<distance<<"; "; rep.printVec("pc pos",pc->pos,G_ndim);
      pc = parent->NextPt_All(pc);
      if (!pc && !loop) { // hack: if get to the end, then go back...
        pc = b->nest.front();
        loop = true;
      }
      distance = grid->idistance(pc->pos, (*bpt)->pos);
    }
    if (!pc) rep.error("BC_assign_COARSE_TO_FINE() left parent grid",0);
    
    // add this parent cell to the "parent" list of this boundary.
    b->nest.push_back(pc);
    (*bpt)->npt = pc;
    ++bpt;
  }  while (bpt !=b->data.end());

  // add data to boundary cells.
  //BC_update_COARSE_TO_FINE(b,OA2,OA2);

  return 0;
}



// ##################################################################
// ##################################################################



int SMR_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
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
  double *U   = new double [par.nvar];
  double *P   = new double [par.nvar];
  double *U1  = new double [par.nvar];
  double *U2  = new double [par.nvar];
  double *P00 = new double [par.nvar];
  double *P10 = new double [par.nvar];
  double *P01 = new double [par.nvar];
  double *P11 = new double [par.nvar];
  double *sx  = new double [par.nvar];
  double *sy  = new double [par.nvar];

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
#ifdef TEST_C2F
    cout <<"C2F: odd step, interpolating coarse data in time.\n";
#endif
    for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
      f = (*f_iter);
      c = f->npt;
      solver->PtoU(c->P, U, par.gamma);
      for (int v=0;v<par.nvar;v++) U[v] += 0.5*c->dU[v];
      solver->UtoP(U,c->Ph, par.EP.MinTemperature, par.gamma);
#ifdef TEST_INF
      for (int v=0;v<par.nvar;v++) {
        if (!isfinite(c->Ph[v])) {
          rep.printVec("NAN c->P ",c->P,par.nvar);
          rep.printVec("NAN c->Ph",c->Ph,par.nvar);
        }
      }
#endif
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
      double f_vol[2], c_vol;
      double dx = fine->DX(); // dx
      //
      // Do two fine cells at a time: they have the same parent.
      // In 1D the geometry is very easy.
      //
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        f = (*f_iter);
        c = f->npt;
        solver->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
        for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] * (1.0-0.5*dx*sx[v]);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;
        
        f_iter++;
        f = (*f_iter);
        for (int v=0;v<par.nvar;v++) f->Ph[v] = c->Ph[v] * (1.0+0.5*dx*sx[v]);
        for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;

        // Now need to check mass/momentum/energy conservation between
        // coarse and fine levels!
        // sum energy of fine cells.
        solver->PtoU(f->Ph, U1, par.gamma);
        f_vol[0] = fine->CellVolume(f);
        f_iter--; f = (*f_iter);
        solver->PtoU(f->Ph, U2, par.gamma);
        f_vol[1] = fine->CellVolume(f);
        for (int v=0;v<par.nvar;v++) U[v] = U1[v]*f_vol[0] + U2[v]*f_vol[1];
        // compare with coarse cell.
        solver->PtoU(c->Ph, P, par.gamma);
        c_vol = coarse->CellVolume(c);
        for (int v=0;v<par.nvar;v++) P[v] *= c_vol;
#ifdef TEST_C2F
        rep.printVec("1D coarse", P,par.nvar);
        rep.printVec("1D fine  ", U,par.nvar);
#endif
        // scale U1, U2 by ratio of coarse to fine energy.
        for (int v=0;v<par.nvar;v++) {
          if (!pconst.equalD(P[v],U[v])) {
            U1[v] *= P[v] / U[v];
            U2[v] *= P[v] / U[v];
          }
        }
#ifdef TEST_C2F
        for (int v=0;v<par.nvar;v++) U[v] = U1[v]+U2[v];
        rep.printVec("1D fine 2", U, par.nvar); 
#endif
        solver->UtoP(U2,f->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
        f_iter++; f = (*f_iter);
        solver->UtoP(U1,f->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
      } // loop over fine cells
    } // 1D

    else if (par.ndim == 2) {
      //
      // Need to do bilinear interpolation, 4 cells at a time.
      //
      double dxo2 = 0.5*fine->DX(); // dx
      double c_vol=0.0, f_vol[4];
      //
      // Do 4 fine cells at a time: they have the same parent.
      // First get the values at the 4 corners of coarse-cell c.
      // Then interpolate with bilinear interpolation to the fine-cell
      // positions as 0.25*dx in each direction from the corners.
      //
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        f = (*f_iter);
        c = f->npt;

        // only do this on every second row (because we update 4
        // cells at a time.
        if (!fine->NextPt(f,YP) || fine->NextPt(f,YP)->npt != c) {
          continue;
        }

        // use slopes in each direction to get corner values for the
        // coarse cell.
        solver->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
        solver->SetSlope(c,YY,par.nvar,sy,OA2,coarse);
        for (int v=0;v<par.nvar;v++) sx[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
        for (int v=0;v<par.nvar;v++) sy[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
        for (int v=0;v<par.nvar;v++) P00[v] = c->Ph[v] -sx[v] -sy[v];
        for (int v=0;v<par.nvar;v++) P10[v] = c->Ph[v] +sx[v] -sy[v];
        for (int v=0;v<par.nvar;v++) P01[v] = c->Ph[v] -sx[v] +sy[v];
        for (int v=0;v<par.nvar;v++) P11[v] = c->Ph[v] +sx[v] +sy[v];

        // now interpolate all four cells using the 4 corner states.
        bilinear_interp(par, c, f, P00, P01, P10, P11);
        bilinear_interp(par, c, fine->NextPt(f,YP), P00, P01, P10, P11);
        f_iter++; f = (*f_iter);
        bilinear_interp(par, c, f, P00, P01, P10, P11);
        bilinear_interp(par, c, fine->NextPt(f,YP), P00, P01, P10, P11);

        f_iter--; f = (*f_iter);
        bilinear_interp(par, c, fine->NextPt(f,YP), P00, P01, P10, P11);
        f_iter++; f = (*f_iter);
        bilinear_interp(par, c, fine->NextPt(f,YP), P00, P01, P10, P11);

        // Now need to check mass/momentum/energy conservation between
        // coarse and fine levels!
        //
        c_vol = coarse->CellVolume(c);
        f_iter--; f = (*f_iter);

        // four have same parent, so calculate conservation.
        solver->PtoU(f->P, P00, par.gamma);
        f_vol[0] = fine->CellVolume(f);
        solver->PtoU(fine->NextPt(f,YP)->P, P10, par.gamma);
        f_vol[1] = fine->CellVolume(fine->NextPt(f,YP));
        f_iter++; f = (*f_iter);
        solver->PtoU(f->P, P01, par.gamma);
        f_vol[2] = fine->CellVolume(f);
        solver->PtoU(fine->NextPt(f,YP)->P, P11, par.gamma);
        f_vol[3] = fine->CellVolume(fine->NextPt(f,YP));

        for (int v=0;v<par.nvar;v++)
          U[v] = P00[v]*f_vol[0] + P10[v]*f_vol[1] + P01[v]*f_vol[2] + P11[v]*f_vol[3];
        // compare with coarse cell.
        solver->PtoU(c->Ph, P, par.gamma);
        for (int v=0;v<par.nvar;v++) P[v] *= c_vol;

#ifdef DEBUG_SMR
        for (int v=0;v<par.nvar;v++) {
          if (!isfinite(P00[v]) || !isfinite(P01[v]) || !isfinite(P10[v]) || !isfinite(P11[v])) {
            rep.printVec("Unscaled fine00",P00,par.nvar);
            rep.printVec("Unscaled fine01",P01,par.nvar);
            rep.printVec("Unscaled fine10",P10,par.nvar);
            rep.printVec("Unscaled fine11",P11,par.nvar);
          }
        }
#endif

        // scale fine conserved vec by adding the difference between
        // conserved quantities on the fine and coarse grids.
        for (int v=0;v<par.nvar;v++) {
          P[v] = 0.25*(P[v] - U[v])/c_vol;
          P00[v] += P[v];
          P01[v] += P[v];
          P10[v] += P[v];
          P11[v] += P[v];
        }
        
        // put scaled conserved variable vectors back into fine cells
        solver->UtoP(P01,f->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];

        solver->UtoP(P11,fine->NextPt(f,YP)->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++)
          fine->NextPt(f,YP)->P[v] = fine->NextPt(f,YP)->Ph[v];

        f_iter--; f = (*f_iter);
        solver->UtoP(P00,f->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];

        solver->UtoP(P10,fine->NextPt(f,YP)->Ph, par.EP.MinTemperature, par.gamma);
        for (int v=0;v<par.nvar;v++)
          fine->NextPt(f,YP)->P[v] = fine->NextPt(f,YP)->Ph[v];
        
#ifdef DEBUG_SMR
        for (int v=0;v<par.nvar;v++) {
          if (!isfinite(f->P[v]) || !isfinite(fine->NextPt(f,YP)->P[v]))
            rep.error("fine 3,4 not finite",f->P[v]);
        }
#endif

        f_iter++; f = (*f_iter);
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

  delete [] U;
  delete [] P;
  delete [] U1;
  delete [] U2;
  delete [] P00;
  delete [] P10;
  delete [] P01;
  delete [] P11;
  delete [] sx;
  delete [] sy;
  return 0;
}



// ##################################################################
// ##################################################################



void SMR_coarse_to_fine_bc::bilinear_interp(
      class SimParams &par,      ///< pointer to simulation parameters
      cell *c,  ///< coarse level cell
      cell *f,  ///< fine level cell
      const double *P00,  ///< prim. vec. at corner of coarse cell
      const double *P01,  ///< prim. vec. at corner of coarse cell
      const double *P10,  ///< prim. vec. at corner of coarse cell
      const double *P11   ///< prim. vec. at corner of coarse cell
      )
{
  if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
    // 1/4,1/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (9.0*P00[v] + 3.0*P10[v] + 3.0*P01[v] +     P11[v]);
  }
  else if ( (f->pos[XX] > c->pos[XX]) && (f->pos[YY] < c->pos[YY]) ) {
    // 3/4,1/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (3.0*P00[v] + 9.0*P10[v] +     P01[v] + 3.0*P11[v]);
  }
  else if ( (f->pos[XX] < c->pos[XX]) && (f->pos[YY] > c->pos[YY]) ) {
    // 1/4,3/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (3.0*P00[v] +     P10[v] + 9.0*P01[v] + 3.0*P11[v]);
  }
  else {
    // 3/4,3/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (    P00[v] + 3.0*P10[v] + 3.0*P01[v] + 9.0*P11[v]);
  }
  for (int v=0;v<par.nvar;v++) f->P[v] = f->Ph[v];
  for (int v=0;v<par.nvar;v++) f->dU[v] = 0.0;

  return;
}



// ##################################################################
// ##################################################################



