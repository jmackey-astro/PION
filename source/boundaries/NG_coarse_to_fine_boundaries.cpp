/// \file NG_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/NG_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

//#define TEST_C2F

// ##################################################################
// ##################################################################



int NG_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
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
  b->NG.clear();

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
        pc = b->NG.front();
        loop = true;
      }
      distance = grid->idistance(pc->pos, (*bpt)->pos);
    }
    if (!pc) rep.error("BC_assign_COARSE_TO_FINE() left parent grid",0);
    
    // add this parent cell to the "parent" list of this boundary.
    b->NG.push_back(pc);
    (*bpt)->npt = pc;
    ++bpt;
  }  while (bpt !=b->data.end());

  // add data to boundary cells.
  //BC_update_COARSE_TO_FINE(b,OA2,OA2);

  return 0;
}



// ##################################################################
// ##################################################################



int NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the NG grid structure
      struct boundary_data *b,
      const int step
      )
{
#ifdef TEST_C2F
    cout <<"C2F: updating boundary data from coarse grid to level ";
    cout <<level<<"\n";
#endif
  //
  // This is a complicated problem to use linear interpolation (or
  // higher order), because you have to conserve mass, momentum and
  // energy between the two levels.  It should also be monotonic and
  // produce cell-averaged values. Flash and pluto seem
  // to use very old Fortran code to accomplish this, and the code is
  // almost impenetrable.  AMRVAC uses linear/bilinear/trilinear
  // interpolation with renormalisation to conserve mass/P/E.
  // I'm doing what AMRVAC does, I think.

  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].parent;
  class GridBaseClass *fine   = par.levels[level].grid;
  list<cell*>::iterator f_iter=b->data.begin();
  cell *f, *c;

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
    double U[par.nvar];
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
  
  // if spatial order of accuracy is 1, then we have piecewise
  // constant data, so there is no interpolation to be done.
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
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        cell *f1, *f2, *c;
        f1 = (*f_iter);
        c = f1->npt;
        f_iter++;
        f2 = (*f_iter);
        interpolate_coarse2fine1D(par,coarse,fine,solver,c,f1,f2);
      }
    } // 1D

    else if (par.ndim == 2) {
#ifdef TEST_C2F
        cout <<"interpolating coarse to fine 2d: ncells="<< b->data.size()<<"\n";
#endif

      //
      // Need to do bilinear interpolation, 4 cells at a time.
      //
      // First get the values at the 4 corners of coarse-cell c.
      // Then interpolate with bilinear interp. to the fine-cell
      // positions as 0.25*dx in each direction from the corners.
      //
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        cell *f1, *f2, *f3, *f4, *c;
        c = (*f_iter)->npt;
        // only do this on every second row because we update 4
        // cells at a time.
        if (!fine->NextPt((*f_iter),YP) ||
             fine->NextPt((*f_iter),YP)->npt != c) {
          continue;
        }
        // get list of four fine cells.
        f1 = (*f_iter);
        f_iter++;
        f2 = (*f_iter);
        f3 = fine->NextPt(f1,YP);
        f4 = fine->NextPt(f2,YP);
        
#ifdef TEST_C2F
        cout <<"interpolating coarse to fine 2d: coarse="<<c->id;
        cout <<", f1="<<f1->id<<", f2="<<f2->id;
        cout <<", f3="<<f3->id<<", f4="<<f4->id<<"\n";
#endif
        interpolate_coarse2fine2D(par,coarse,fine,solver,c,f1,f2,f3,f4);

      } // loop over fine cells
    } // 2D

    else if (par.ndim == 3) {
      //
      // Need to do trilinear interpolation on 8 fine cells.
      //
      rep.error("3D coarse-to-fine interpolation at 2nd order!",3);
    } // 3D
  } // 2nd-order accuracy

  return 0;
}



// ##################################################################
// ##################################################################



void NG_coarse_to_fine_bc::bilinear_interp(
      class SimParams &par,      ///< pointer to simulation parameters
      cell *c,  ///< coarse level cell
      cell *f,  ///< fine level cell
      const double *P00,  ///< prim. vec. at XN,YN corner of cell
      const double *P01,  ///< prim. vec. at XP,YN corner of cell
      const double *P10,  ///< prim. vec. at YP,XN corner of cell
      const double *P11   ///< prim. vec. at XP,YP corner of cell
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



void NG_coarse_to_fine_bc::interpolate_coarse2fine1D(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *coarse,  ///< pointer to coarse grid
      class GridBaseClass *fine,    ///< pointer to fine grid
      class FV_solver_base *solver, ///< pointer to equations
      cell *c,  ///< pointer to cell on coarse grid
      cell *f1, ///< pointer to first fine cell  (XN)
      cell *f2  ///< pointer to second fine cell (XP)
      )
{
  double sx[par.nvar]; // slope in x-dir
  double fU[par.nvar], f1U[par.nvar], f2U[par.nvar], cU[par.nvar];
  double f_vol[2], c_vol;
  double dx = fine->DX(); // dx
  //
  // In 1D the geometry is very easy.
  //
  solver->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
  for (int v=0;v<par.nvar;v++) f1->Ph[v] = c->Ph[v] * (1.0-0.5*dx*sx[v]);
  for (int v=0;v<par.nvar;v++) f1->P[v] = f1->Ph[v];
  for (int v=0;v<par.nvar;v++) f1->dU[v] = 0.0;
    
  for (int v=0;v<par.nvar;v++) f2->Ph[v] = c->Ph[v] * (1.0+0.5*dx*sx[v]);
  for (int v=0;v<par.nvar;v++) f2->P[v] = f2->Ph[v];
  for (int v=0;v<par.nvar;v++) f2->dU[v] = 0.0;

  // Now need to check mass/momentum/energy conservation between
  // coarse and fine levels (Berger & Colella, 1989)
  // sum energy of fine cells.
  solver->PtoU(f1->Ph, f1U, par.gamma);
  f_vol[0] = fine->CellVolume(f1);
  solver->PtoU(f2->Ph, f2U, par.gamma);
  f_vol[1] = fine->CellVolume(f2);
  for (int v=0;v<par.nvar;v++) fU[v] = f1U[v]*f_vol[0] + f2U[v]*f_vol[1];
  // compare with coarse cell.
  solver->PtoU(c->Ph, cU, par.gamma);
  c_vol = coarse->CellVolume(c);
  for (int v=0;v<par.nvar;v++) cU[v] *= c_vol;
#ifdef TEST_C2F
  rep.printVec("1D coarse", cU,par.nvar);
  rep.printVec("1D fine  ", fU,par.nvar);
#endif
  // scale f1U, f2U by ratio of coarse to fine energy.
  // scale fine conserved vec by adding the difference between
  // conserved quantities on the fine and coarse grids.
  for (int v=0;v<par.nvar;v++) cU[v] = 0.5*(cU[v] - fU[v])/c_vol;
  for (int v=0;v<par.nvar;v++) f1U[v] += cU[v];
  for (int v=0;v<par.nvar;v++) f2U[v] += cU[v];
#ifdef TEST_C2F
  for (int v=0;v<par.nvar;v++) fU[v] = f1U[v]+f2U[v];
  rep.printVec("1D fine 2", fU, par.nvar); 
#endif
  solver->UtoP(f2U,f2->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f2->P[v] = f2->Ph[v];
  solver->UtoP(f1U,f1->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f1->P[v] = f1->Ph[v];

  return;
}



// ##################################################################
// ##################################################################



void NG_coarse_to_fine_bc::interpolate_coarse2fine2D(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *coarse,  ///< pointer to coarse grid
      class GridBaseClass *fine,    ///< pointer to fine grid
      class FV_solver_base *solver, ///< pointer to equations
      cell *c,  ///< pointer to cell on coarse grid
      cell *f1, ///< pointer to first fine cell  (XN,YN)
      cell *f2, ///< pointer to second fine cell (XP,YN)
      cell *f3, ///< pointer to third fine cell  (XN,YP)
      cell *f4  ///< pointer to fourth fine cell (XP,YP)
      )
{
  double sx[par.nvar], sy[par.nvar]; // slopes in coarse cell
  double fU[par.nvar], f1U[par.nvar], f2U[par.nvar];
  double f3U[par.nvar], f4U[par.nvar], cU[par.nvar];
  double dxo2 = 0.5*fine->DX(); // dx
  double c_vol=0.0, f_vol[4];
  //
  // Need to do bilinear interpolation, 4 cells at a time.
  // use slopes in each direction to get corner values for the
  // coarse cell.
  //
  solver->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
  solver->SetSlope(c,YY,par.nvar,sy,OA2,coarse);
  for (int v=0;v<par.nvar;v++) sx[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
  for (int v=0;v<par.nvar;v++) sy[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
  for (int v=0;v<par.nvar;v++) f1U[v] = c->Ph[v] -sx[v] -sy[v];
  for (int v=0;v<par.nvar;v++) f2U[v] = c->Ph[v] +sx[v] -sy[v];
  for (int v=0;v<par.nvar;v++) f3U[v] = c->Ph[v] -sx[v] +sy[v];
  for (int v=0;v<par.nvar;v++) f4U[v] = c->Ph[v] +sx[v] +sy[v];

  // interpolate all four cells using the 4 corner states.
  bilinear_interp(par, c, f1, f1U, f2U, f3U, f4U);
  bilinear_interp(par, c, f2, f1U, f2U, f3U, f4U);
  bilinear_interp(par, c, f3, f1U, f2U, f3U, f4U);
  bilinear_interp(par, c, f4, f1U, f2U, f3U, f4U);

  // Need to check mass/momentum/energy conservation between
  // coarse and fine levels
  //
  c_vol = coarse->CellVolume(c);
  solver->PtoU(f1->P, f1U, par.gamma);
  solver->PtoU(f2->P, f2U, par.gamma);
  solver->PtoU(f3->P, f3U, par.gamma);
  solver->PtoU(f4->P, f4U, par.gamma);
  f_vol[0] = fine->CellVolume(f1);
  f_vol[1] = fine->CellVolume(f2);
  f_vol[2] = fine->CellVolume(f3);
  f_vol[3] = fine->CellVolume(f4);

  for (int v=0;v<par.nvar;v++)
    fU[v] = f1U[v]*f_vol[0] + f3U[v]*f_vol[1] +
            f2U[v]*f_vol[2] + f4U[v]*f_vol[3];
  // compare with coarse cell.
  solver->PtoU(c->Ph, cU, par.gamma);
  for (int v=0;v<par.nvar;v++) cU[v] *= c_vol;

#ifdef DEBUG_NG
  for (int v=0;v<par.nvar;v++) {
    if (!isfinite(f1U[v]) || !isfinite(f1U[v]) ||
        !isfinite(f3U[v]) || !isfinite(f4U[v])) {
      rep.printVec("Unscaled fine00",f1U,par.nvar);
      rep.printVec("Unscaled fine10",f2U,par.nvar);
      rep.printVec("Unscaled fine01",f3U,par.nvar);
      rep.printVec("Unscaled fine11",f4U,par.nvar);
    }
  }
#endif

  // scale fine conserved vec by adding the difference between
  // conserved quantities on the fine and coarse grids.
  for (int v=0;v<par.nvar;v++) cU[v] = 0.25*(cU[v] - fU[v])/c_vol;
  for (int v=0;v<par.nvar;v++) f1U[v] += cU[v];
  for (int v=0;v<par.nvar;v++) f2U[v] += cU[v];
  for (int v=0;v<par.nvar;v++) f3U[v] += cU[v];
  for (int v=0;v<par.nvar;v++) f4U[v] += cU[v];
  
  // put scaled conserved variable vectors back into fine cells
  solver->UtoP(f1U,f1->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f1->P[v] = f1->Ph[v];
  solver->UtoP(f2U,f2->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f2->P[v] = f2->Ph[v];
  solver->UtoP(f3U,f3->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f3->P[v] = f3->Ph[v];
  solver->UtoP(f4U,f4->Ph, par.EP.MinTemperature, par.gamma);
  for (int v=0;v<par.nvar;v++) f4->P[v] = f4->Ph[v];

#ifdef DEBUG_NG
  for (int v=0;v<par.nvar;v++) {
    if (!isfinite(f3->P[v]) || !isfinite(f4->P[v]))
      rep.error("fine 3,4 not finite",f3->P[v]);
  }
#endif
  return;
}



// ##################################################################
// ##################################################################


