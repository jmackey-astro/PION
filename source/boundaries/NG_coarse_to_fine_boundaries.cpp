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
      class SimParams &par,     ///< simulation parameters
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
    rep.error("BC_assign_COARSE_TO_FINE: empty boundary data",
                                                        b->itype);

  //
  // If we are doing raytracing, then also send the column densities
  // from coarse to fine grids.  Here we set up the number of them.
  //
  struct rad_src_info *s;
  C2F_Nxd = 0;
  C2F_tauoff.resize(par.RS.Nsources);
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    s = &(par.RS.sources[isrc]);
    C2F_tauoff[isrc] = C2F_Nxd;
    C2F_Nxd += 2*s->NTau; // Need col2cell and cell_col for each Tau.
#ifdef TEST_MPI_NG
    cout <<"C2F_BC: RT Source "<<isrc<<": adding "<<2*s->NTau;
    cout <<" cols, running total = "<<C2F_Nxd<<"\n";
#endif
  }

  list<cell*>::iterator bpt=b->data.begin();
  int gidx = grid->idx();
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
      pc = parent->NextPt_All(pc);
      if (!pc && !loop) { // hack: if get to the end, then go back...
        pc =  parent->FirstPt_All();
        loop = true;
      }
      distance = grid->idistance(pc->pos, (*bpt)->pos);
    }
    if (!pc)
      rep.error("BC_assign_COARSE_TO_FINE() left parent grid",0);
    
    // set boundary cell's 'npt' pointer to point to the parent cell.
    (*bpt)->npt = pc;
    ++bpt;
  }  while (bpt !=b->data.end());

  return 0;
}



// ##################################################################
// ##################################################################



int NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(
      class SimParams &par,      ///< simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the NG grid structure
      struct boundary_data *b,
      const int step
      )
{
  if (level==0) return 0;
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
  // PION is doing what AMRVAC does, I think.
  //

  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].parent;
  class GridBaseClass *fine   = par.levels[level].grid;
#ifdef TEST_C2F
  cout <<"C2F: fine="<<fine<<", parent="<<coarse<<", level=";
  cout <<level<<"\n";
#endif

  list<cell*>::iterator f_iter=b->data.begin();
  cell *f, *c;

  // ----------------------------------------------------------------
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
#ifdef TEST_C2F
      for (int v=0;v<par.nvar;v++) {
        if (!isfinite(c->Ph[v])) {
          rep.printVec("NAN c->P ",c->P,par.nvar);
          rep.printVec("NAN c->Ph",c->Ph,par.nvar);
        }
      }
#endif
    }
  }
  // ----------------------------------------------------------------
  
  // ----------------------------------------------------------------
  // save coarse cell fluxes into array, and then interpolate onto
  // fine cells.
  for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
    cell *f1=0, *f2=0;
    struct rad_src_info *s;
    int off = 0, ncells=0;
    pion_flt Tau[MAX_TAU], T[C2F_Nxd], cpos[MAX_DIM];
    std::vector<cell *> fcl;
    c = (*f_iter)->npt;
    CI.get_dpos(c,cpos);

    // save Tau and dTau for each source in array T.
    for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      off = C2F_tauoff[isrc];
      CI.get_col(c, s->id, Tau);
      for (int iT=0; iT<s->NTau; iT++) T[off+iT] = Tau[iT];
      off += s->NTau;
      CI.get_cell_col(c, s->id, Tau);
      for (int iT=0; iT<s->NTau; iT++) T[off+iT] = Tau[iT];
    }
    //rep.printVec("T coarse",T,C2F_Nxd);
    
    // add fine cells to list.
    if      (par.ndim == 1) {
      fcl.push_back(*f_iter);
      f_iter++;
      fcl.push_back(*f_iter);
      ncells=2;
    }
    else if (par.ndim == 2) {
      if (!fine->NextPt((*f_iter),YP) ||
           fine->NextPt((*f_iter),YP)->npt != c) {
        //cout <<"skipping column\n";
        ncells=0;
        break;
      }
      // get list of four fine cells.
      fcl.push_back(*f_iter); f1 = (*f_iter);
      f_iter++;
      fcl.push_back(*f_iter); f2 = (*f_iter);
      fcl.push_back(fine->NextPt(f1,YP));
      fcl.push_back(fine->NextPt(f2,YP));
      ncells=4;
    }
    else {
      rep.error("C2F RT 3D not implemented",par.ndim);
    }
    if (ncells==0) {
      //cout <<"ncells=0: skipping column\n";
      continue; // must be on odd-numbered column
    }

    //rep.printVec("f0 pos",fcl[0]->pos,2);
    //rep.printVec("f1 pos",fcl[1]->pos,2);
    //rep.printVec("f2 pos",fcl[2]->pos,2);
    //rep.printVec("f3 pos",fcl[3]->pos,2);
    //rep.printVec("c  pos",cpos,2);

    // This function assigns interpolated Tau,dTau to fine cells.
    get_C2F_Tau(par,fcl,cpos,T);
  } // loop over cells and assign optical depths
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
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
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  else if (par.spOOA == OA2) {
    //
    // Dimensions is sufficiently different that we have an if/else
    // loop for each dimension, and then do linear/bilinear/trilinear
    // interpolation as needed.
    //
    if (par.ndim ==1) {
      pion_flt slope[par.nvar], c_vol=0;
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        cell *f1, *f2, *c;
        f1 = (*f_iter);
        f_iter++;
        f2 = (*f_iter);
        // coarse cell properties:
        c = f1->npt;
        c_vol = coarse->CellVolume(c,0);
        solver->SetSlope(c,XX,par.nvar,slope,OA2,coarse);
        interpolate_coarse2fine1D(
                            par,fine,solver,c->Ph,c_vol,slope,f1,f2);
      }
    } // 1D

    else if (par.ndim == 2) {
#ifdef TEST_C2F
      cout <<"interpolating coarse to fine 2d: ncells=";
      cout << b->data.size()<<"\n";
#endif
      pion_flt sx[par.nvar], sy[par.nvar], c_vol=0;

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
        c_vol = coarse->CellVolume(c,0);
        solver->SetSlope(c,XX,par.nvar,sx,OA2,coarse);
        solver->SetSlope(c,XX,par.nvar,sy,OA2,coarse);
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
        cout <<"BEFORE interpolating coarse to fine 2d: coarse="<<c->id;
        cout <<", f1="<<f1->id<<", f2="<<f2->id;
        cout <<", f3="<<f3->id<<", f4="<<f4->id<<"\n";
        CI.print_cell(c);
        //CI.print_cell(f1);
        //CI.print_cell(f2);
        //CI.print_cell(f3);
        //CI.print_cell(f4);
#endif
        interpolate_coarse2fine2D(
              par,fine,solver,c->Ph,c->pos,c_vol,sx,sy,f1,f2,f3,f4);
        
#ifdef TEST_C2F
        cout <<"AFTER interpolating coarse to fine 2d: coarse="<<c->id;
        cout <<", f1="<<f1->id<<", f2="<<f2->id;
        cout <<", f3="<<f3->id<<", f4="<<f4->id<<"\n";
        CI.print_cell(c);
        CI.print_cell(f1);
        //CI.print_cell(f2);
        //CI.print_cell(f3);
        //CI.print_cell(f4);
#endif

      } // loop over fine cells
    } // 2D

    else if (par.ndim == 3) {
      //
      // Need to do trilinear interpolation on 8 fine cells.
      //
      rep.error("3D coarse-to-fine interpolation at 2nd order!",3);
    } // 3D
  } // 2nd-order accuracy
  // ----------------------------------------------------------------

  return 0;
}



// ##################################################################
// ##################################################################



void NG_coarse_to_fine_bc::bilinear_interp(
      class SimParams &par,      ///< pointer to simulation parameters
      const int *cpos,  ///< coarse level cell integer position
      cell *f,  ///< fine level cell
      const double *P00,  ///< prim. vec. at XN,YN corner of cell
      const double *P01,  ///< prim. vec. at XP,YN corner of cell
      const double *P10,  ///< prim. vec. at YP,XN corner of cell
      const double *P11   ///< prim. vec. at XP,YP corner of cell
      )
{
  if ( (f->pos[XX] < cpos[XX]) && (f->pos[YY] < cpos[YY]) ) {
    // 1/4,1/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (9.0*P00[v] + 3.0*P10[v] + 3.0*P01[v] +     P11[v]);
  }
  else if ( (f->pos[XX] > cpos[XX]) && (f->pos[YY] < cpos[YY]) ) {
    // 3/4,1/4
    for (int v=0;v<par.nvar;v++) f->Ph[v] = 1.0/16.0 *
                  (3.0*P00[v] + 9.0*P10[v] +     P01[v] + 3.0*P11[v]);
  }
  else if ( (f->pos[XX] < cpos[XX]) && (f->pos[YY] > cpos[YY]) ) {
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
      class SimParams &par,         ///< simulation parameters
      class GridBaseClass *fine,    ///< pointer to fine grid
      class FV_solver_base *solver, ///< pointer to equations
      const pion_flt *P,    ///< state vector of coarse cell.
      const pion_flt c_vol, ///< volume of coarse cell.
      pion_flt *sx,   ///< dP/dx in coarse cell.
      cell *f1, ///< pointer to first fine cell  (XN)
      cell *f2  ///< pointer to second fine cell (XP)
      )
{
  double fU[par.nvar], f1U[par.nvar], f2U[par.nvar], cU[par.nvar];
  double f_vol[2];
  double dx = fine->DX(); // dx
  //
  // In 1D the geometry is very easy.
  //
  for (int v=0;v<par.nvar;v++) f1->Ph[v] = P[v] * (1.0-0.5*dx*sx[v]);
  for (int v=0;v<par.nvar;v++) f1->P[v] = f1->Ph[v];
  for (int v=0;v<par.nvar;v++) f1->dU[v] = 0.0;
    
  for (int v=0;v<par.nvar;v++) f2->Ph[v] = P[v] * (1.0+0.5*dx*sx[v]);
  for (int v=0;v<par.nvar;v++) f2->P[v] = f2->Ph[v];
  for (int v=0;v<par.nvar;v++) f2->dU[v] = 0.0;

  // Now need to check mass/momentum/energy conservation between
  // coarse and fine levels (Berger & Colella, 1989)
  // sum energy of fine cells.
  solver->PtoU(f1->Ph, f1U, par.gamma);
  f_vol[0] = fine->CellVolume(f1,0);
  solver->PtoU(f2->Ph, f2U, par.gamma);
  f_vol[1] = fine->CellVolume(f2,0);
  for (int v=0;v<par.nvar;v++) fU[v] = f1U[v]*f_vol[0] + f2U[v]*f_vol[1];
  // compare with coarse cell.
  solver->PtoU(P, cU, par.gamma);
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
      class GridBaseClass *fine,    ///< pointer to fine grid
      class FV_solver_base *solver, ///< pointer to equations
      const pion_flt *P,    ///< state vector of coarse cell.
      const int *cpos,      ///< position of coarse cell.
      const pion_flt c_vol, ///< volume of coarse cell.
      pion_flt *sx,   ///< dP/dx in coarse cell.
      pion_flt *sy,   ///< dP/dy in coarse cell.
      cell *f1, ///< pointer to first fine cell  (XN,YN)
      cell *f2, ///< pointer to second fine cell (XP,YN)
      cell *f3, ///< pointer to third fine cell  (XN,YP)
      cell *f4  ///< pointer to fourth fine cell (XP,YP)
      )
{
  double fU[par.nvar], f1U[par.nvar], f2U[par.nvar];
  double f3U[par.nvar], f4U[par.nvar], cU[par.nvar];
  double dxo2 = 0.5*fine->DX(); // dx
  double f_vol[4];
  //
  // Need to do bilinear interpolation, 4 cells at a time.
  // use slopes in each direction to get corner values for the
  // coarse cell.
  //
  for (int v=0;v<par.nvar;v++) sx[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
  for (int v=0;v<par.nvar;v++) sy[v] *= 2.0*dxo2; // coarse dx/2 = fine 2*(dx/2)
  for (int v=0;v<par.nvar;v++) f1U[v] = P[v] -sx[v] -sy[v];
  for (int v=0;v<par.nvar;v++) f2U[v] = P[v] +sx[v] -sy[v];
  for (int v=0;v<par.nvar;v++) f3U[v] = P[v] -sx[v] +sy[v];
  for (int v=0;v<par.nvar;v++) f4U[v] = P[v] +sx[v] +sy[v];

  // interpolate all four cells using the 4 corner states.
  bilinear_interp(par, cpos, f1, f1U, f2U, f3U, f4U);
  bilinear_interp(par, cpos, f2, f1U, f2U, f3U, f4U);
  bilinear_interp(par, cpos, f3, f1U, f2U, f3U, f4U);
  bilinear_interp(par, cpos, f4, f1U, f2U, f3U, f4U);

  // Need to check mass/momentum/energy conservation between
  // coarse and fine levels
  //
  solver->PtoU(f1->P, f1U, par.gamma);
  solver->PtoU(f2->P, f2U, par.gamma);
  solver->PtoU(f3->P, f3U, par.gamma);
  solver->PtoU(f4->P, f4U, par.gamma);
  f_vol[0] = fine->CellVolume(f1,0);
  f_vol[1] = fine->CellVolume(f2,0);
  f_vol[2] = fine->CellVolume(f3,0);
  f_vol[3] = fine->CellVolume(f4,0);

  for (int v=0;v<par.nvar;v++)
    fU[v] = f1U[v]*f_vol[0] + f3U[v]*f_vol[1] +
            f2U[v]*f_vol[2] + f4U[v]*f_vol[3];
  // compare with coarse cell.
  solver->PtoU(P, cU, par.gamma);
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

  //int d=par.nvar-1;
  //cout <<"interpolate_coarse2fine2D: "<<P[d]<<": "<<f1->P[d]<<", ";
  //cout <<f2->P[d]<<", "<<f3->P[d]<<", "<<f4->P[d]<<"\n";

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



void NG_coarse_to_fine_bc::get_C2F_Tau(
      class SimParams &par,   ///< pointer to simulation parameters
      std::vector<cell *> &c, ///< list of cells
      const pion_flt *cpos,   ///< centre of coarse cell.
      double *T               ///< coarse-cell optical depths
      )
{
  // data needed for getting coarse cell Taus onto coarse grid.
  double Tau[MAX_TAU], dTau[MAX_TAU];
  class cell *f0, *f1, *f2, *f3;
  struct rad_src_info *s;
  double diffx,diffy;

  if (par.ndim == 1) {
    f0 = c[0];
    f1 = c[1];

    for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      for (int iT=0; iT<s->NTau; iT++)
        Tau[iT] = T[C2F_tauoff[isrc]+iT];
      for (int iT=0; iT<s->NTau; iT++)
        dTau[iT] = T[C2F_tauoff[isrc]+s->NTau+iT];
      // column from source through cell depends on
      // which direction is to the source
      if (s->pos[XX] > cpos[XX]) {
        // f0 has same column as coarse, f1 has 0.5*dTau subtracted.
        // Both cells have dTau divided by 2.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
      }
      else {
        // Other way around.  Now f0 is closer to source, so we
        // subtract 0.5*dTau from it.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
      }
    } // loop over sources
  } // if 1D

  else if (par.ndim == 2) {
#ifdef RT_TESTING
    cout <<"2D C2F RT Routine\n";
#endif
    f0 = c[0];
    f1 = c[1];
    f2 = c[2];
    f3 = c[3];
#ifdef RT_TESTING
    cout <<"f: ["<<f0->pos[XX]<<","<<f0->pos[YY]<<"], c: [";
    cout <<cpos[XX]<<","<<cpos[YY]<<"] : ";
#endif
    for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      for (int iT=0; iT<s->NTau; iT++)
        Tau[iT] = T[C2F_tauoff[isrc]+iT];
      for (int iT=0; iT<s->NTau; iT++)
        dTau[iT] = T[C2F_tauoff[isrc]+s->NTau+iT];
      diffx = s->pos[XX] - cpos[XX];
      diffy = s->pos[YY] - cpos[YY];
#ifdef RT_TESTING
      cout <<"diffxy = "<<diffx<<"  "<<diffy<<"\n";
      cout <<"t1="<<*Tau1<<", t2="<<*Tau2<<", t3=";
      cout <<*Tau3<<", t4="<<*Tau4;
#endif
      // column from source through cell depends on
      // which direction is to the source.
      //
      /// We assume we are in the far-field limit, so that the
      /// rays through the coarse and 4 fine cells are essentially
      /// parallel to each other, and the ray segments throug the 
      /// fine cells have the same length.
      //
      if      (diffx>0 && fabs(diffx)>=fabs(diffy)) {
        // Source in Q1 coming from dir XP (-45 < theta < 45 deg)
        // f2,f0 have Tau same as coarse cell, and f1,f3 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
        CI.set_col(f2, s->id, Tau);
        CI.set_cell_col(f2, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
        CI.set_col(f3, s->id, Tau);
        CI.set_cell_col(f3, s->id, dTau);
      }
      else if (diffy>0 && fabs(diffx)<fabs(diffy)) {
        // source in Q2, coming from dir YP (45 < theta < 135 deg)
        // f0,f1 have Tau same as coarse cell, and f2,f3 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f2, s->id, Tau);
        CI.set_cell_col(f2, s->id, dTau);
        CI.set_col(f3, s->id, Tau);
        CI.set_cell_col(f3, s->id, dTau);
      }
      else if (diffx<0 && fabs(diffx)>=fabs(diffy)) {
        // source in Q3, coming from XN (135 < theta < 225 deg)
        // f1,f3 have Tau same as coarse cell, and f0,f2 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
        CI.set_col(f3, s->id, Tau);
        CI.set_cell_col(f3, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
        CI.set_col(f2, s->id, Tau);
        CI.set_cell_col(f2, s->id, dTau);
      }
      else {
        // source in Q4, coming from YN (225 < theta < 315 deg)
        // f2,f3 have Tau same as coarse cell, and f0,f1 have Tau
        // reduced by dTau/2
        // All 4 cells have dTau reduced by a factor of 2.
        for (int iT=0; iT<s->NTau; iT++) dTau[iT] *= 0.5;
        CI.set_col(f2, s->id, Tau);
        CI.set_cell_col(f2, s->id, dTau);
        CI.set_col(f3, s->id, Tau);
        CI.set_cell_col(f3, s->id, dTau);
        for (int iT=0; iT<s->NTau; iT++) Tau[iT] -= dTau[iT];
        CI.set_col(f0, s->id, Tau);
        CI.set_cell_col(f0, s->id, dTau);
        CI.set_col(f1, s->id, Tau);
        CI.set_cell_col(f1, s->id, dTau);
      }
#ifdef RT_TESTING
      cout <<"  tc="<<*tmp;
#endif
    } // loop over sources.
  } // if 2D

  else {
    rep.error("3D C2F RT not implemented in NG grid",par.ndim);
  } // if 3D
} // if there is a finer level



// ##################################################################
// ##################################################################







