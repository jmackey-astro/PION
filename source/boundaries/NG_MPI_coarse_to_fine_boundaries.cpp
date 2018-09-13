/// \file NG_MPI_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_MPI_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.09.06 JM: started writing code.


#include "boundaries/NG_MPI_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

//#define TEST_C2F

// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_SEND(
      class SimParams &par,  ///< simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b,  ///< boundary data
      )
{

  class GridBaseClass *grid = par.levels[l].grid;
  int gidx = grid->idx();

  // see how many child grids I have
  class MCMDcontrol *MCMD = &(par.level[l].MCMD);
  int nchild = MCMD->child_procs.size();
  b->NGsend.clear();

  // loop over child grids
  for (int i=0;i<nchild;i++) {

    if (MCMD->get_myrank() == MCMD->child_procs[i].rank) {
      // if child is on my process, do nothing because child grid
      // can grab the data directly.
#ifdef TEST_MPI_NG
      cout <<"my rank ("<<MCMD->get_myrank()<<") == child rank (";
      cout <<MCMD->child_procs[i].rank<<"), no need to set up ";
      cout <<"COARSE_TO_FINE_SEND\n";
#endif
    }
    else {
      // if child is on another process, make a list of cells that
      // we need to send.  One list for each external boundary that
      // we need to update.  Only external boundaries of the whole
      // domain are included.
#ifdef TEST_MPI_NG
      cout <<"my rank != child rank, running parallel ";
      cout <<"COARSE_TO_FINE_SEND\n";
#endif
      // get dimensions of child grid from struct
      int ixmin[MAX_DIM], ixmax[MAX_DIM];
      CI.get_ipos_vec(MCMD->child_procs[i].xmin, ixmin);
      CI.get_ipos_vec(MCMD->child_procs[i].xmax, ixmax);

      // loop over dimensions
      for (int d=0;d<par.ndim;d++) {
        // if child xmin == its level xmin, but > my level xmin,
        // then we need to send data, so set up a list.
        if ( (pconst.equalD(MCMD->child_procs[i].xmin[d],
                            par.level[l+1].xmin[d]))        &&
             (MCMD->child_procs[i].xmin[d] > 
              par.level[l].xmin[d]*ONE_PLUS_EPS) ) {
          struct c2f *bdata = new struct c2f;
          bdata->rank = MCMD->child_procs[i].rank;
          bdata->dir  = 2*d;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(par,grid,bdata,ixmin,ixmax);
        }
        // if child xmax == its level xmax, but < my level xmax,
        // then we need to send data, so set up a list.
        else if ((pconst.equalD(MCMD->child_procs[i].xmax[d],
                                par.level[l+1].xmax[d]))    &&
                 (MCMD->child_procs[i].xmax[d] < 
                  par.level[l].xmax[d]*ONE_MINUS_EPS) ) {
          struct c2f *bd = new struct c2f;
          bdata->rank = MCMD->child_procs[i].rank;
          bdata->dir  = 2*d+1;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(par,grid,bdata,ixmin,ixmax);
          b->NGsend.push_back(bdata);
        }
      } // loop over dimensions
    } // if child is not on my process
  } // loop over child grids
  return 0;
}



// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_SEND(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int l, ///< level in the NG grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  return 0;
}



// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_RECV(
      class SimParams &par,     ///< pointer to simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b,  ///< boundary data
      )
{
  return 0;
}




// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_RECV(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int l, ///< level in the NG grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
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



void NG_MPI_coarse_to_fine_bc::bilinear_interp(
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



void NG_MPI_coarse_to_fine_bc::interpolate_coarse2fine1D(
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



void NG_MPI_coarse_to_fine_bc::interpolate_coarse2fine2D(
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



void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      )
{
  // In XN,XP direction we add cells with faces that touch the fine-
  // level grid from the outside.
  // In YN,YP direction we also add corner data
  // In ZN,ZP direction we also add edge data
  //
  // easier to have a function for different grid dimensions.
  if (par.ndim==1)
    add_cells_to_C2F_send_list_1D(par,grid,bdata,ixmin,ixmax);
  else if (par.ndim==2)
    add_cells_to_C2F_send_list_2D(par,grid,bdata,ixmin,ixmax);
  else
    add_cells_to_C2F_send_list_3D(par,grid,bdata,ixmin,ixmax);

  return;
}



// ##################################################################
// ##################################################################



void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_1D(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      )
{
  cell *c = grid->FirstPt();
  int bsize = grid->idx()*par.Nbc/2; // idx is >=2, Nbc is >=1.
  
  // define domain of boundary region
  int xn,xp;
  switch (bdata->dir) {
    case XN:
    xn = ixmin[XX]-bsize;
    xp = ixmin[XX];
    break;

    case XP:
    xn = ixmax[XX];
    xp = ixmax[XX]+bsize;
    break;

    default:
    rep.error("bad direction in 1D C2F",bdata->dir);
  }

  do {
    if (c->pos[XX]>xn && c->pos[XX]<xp)
      bdata->c.push_back(c);
  } while ((c=grid->NextPt(c)) !=0);

  return;
}



// ##################################################################
// ##################################################################



void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_2D(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      )
{
  cell *c = grid->FirstPt();
  int bsize = grid->idx()*par.Nbc/2; // idx is >=2, Nbc is >=1.
  
  // define domain of boundary region
  int xn,xp,yn,yp;
  switch (bdata->dir) {
    case XN:
    xn = ixmin[XX]-bsize;
    xp = ixmin[XX];
    yn = ixmin[YY];
    yp = ixmax[YY];
    break;

    case XP:
    xn = ixmax[XX];
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY];
    yp = ixmax[YY];
    break;

    case YN:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY]-bsize;
    yp = ixmin[YY];
    break;

    case YP:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmax[YY];
    yp = ixmax[YY]+bsize;
    break;


    default:
    rep.error("bad direction in 2D C2F",bdata->dir);
  }

  do {
    if (c->pos[XX]>xn && c->pos[XX]<xp &&
        c->pos[YY]>yn && c->pos[YY]<yp)
      bdata->c.push_back(c);
  } while ((c=grid->NextPt(c)) !=0);

  return;
}



// ##################################################################
// ##################################################################



void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_3D(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to coarse-level grid
      struct c2f *bdata,          ///< pointer to list of cells
      int *ixmin,                 ///< child grid xmin (integer)
      int *ixmax                  ///< child grid xmax (integer)
      )
{
  cell *c = grid->FirstPt();
  int bsize = grid->idx()*par.Nbc/2; // idx is >=2, Nbc is >=1.
  
  // define domain of boundary region
  int xn,xp,yn,yp,zn,zp;
  switch (bdata->dir) {
    case XN:
    xn = ixmin[XX]-bsize;
    xp = ixmin[XX];
    yn = ixmin[YY];
    yp = ixmax[YY];
    zn = ixmin[ZZ];
    zp = ixmax[ZZ];
    break;

    case XP:
    xn = ixmax[XX];
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY];
    yp = ixmax[YY];
    zn = ixmin[ZZ];
    zp = ixmax[ZZ];
    break;

    case YN:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY]-bsize;
    yp = ixmin[YY];
    zn = ixmin[ZZ];
    zp = ixmax[ZZ];
    break;

    case YP:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmax[YY];
    yp = ixmax[YY]+bsize;
    zn = ixmin[ZZ];
    zp = ixmax[ZZ];
    break;
    
    case ZN:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY]-bsize;
    yp = ixmax[YY]+bsize;
    zn = ixmin[ZZ]-bsize;
    zp = ixmin[ZZ];
    break;

    case ZP:
    xn = ixmin[XX]-bsize;
    xp = ixmax[XX]+bsize;
    yn = ixmin[YY]-bsize;
    yp = ixmax[YY]+bsize;
    zn = ixmax[ZZ];
    zp = ixmax[ZZ]+bsize;
    break;

    default:
    rep.error("bad direction in 3D C2F",bdata->dir);
  }

  do {
    if (c->pos[XX]>xn && c->pos[XX]<xp &&
        c->pos[YY]>yn && c->pos[YY]<yp &&
        c->pos[ZZ]>zn && c->pos[ZZ]<zp)
      bdata->c.push_back(c);
  } while ((c=grid->NextPt(c)) !=0);

  return;
}



// ##################################################################
// ##################################################################





