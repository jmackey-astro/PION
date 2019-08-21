/// \file NG_fine_to_coarse_boundaries.cpp
/// \brief Class definitions for NG_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/NG_fine_to_coarse_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################



int NG_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE(
      class SimParams &par,     ///< simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,  ///< boundary data
      class GridBaseClass *child,  ///< pointer to child grid.
      const int i        ///< which element of NGrecvF2C to use
      )
{
  //
  // Make a list of child-grid cells to map onto the coarse grid
  //
  if (b->NGrecvF2C.size() <= static_cast<unsigned int>(i))
    rep.error("BC_assign_FINE_TO_COARSE: recv vec too small",i);
  if (b->NGrecvF2C[i].empty())
    rep.error("BC_assign_FINE_TO_COARSE: empty boundary data",
                                                        b->itype);

  //
  // If we are doing raytracing, then also send the column densities
  // from fine to coarse grids.  Here we set up the number of them.
  //
  struct rad_src_info *s;
  F2C_Nxd = 0;
  F2C_tauoff.resize(par.RS.Nsources);
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    s = &(par.RS.sources[isrc]);
    F2C_tauoff[isrc] = F2C_Nxd;
    F2C_Nxd += 2*s->NTau; // Need col2cell and cell_col for each Tau.
#ifdef TEST_MPI_NG
    cout <<"F2C_BC: RT Source "<<isrc<<": adding "<<2*s->NTau;
    cout <<" cols, running total = "<<F2C_Nxd<<"\n";
#endif
  }

  int nc=1;  // number of fine cells per coarse cell
  for (int id=0;id<par.ndim;id++) nc*=2;
  int nel = b->NGrecvF2C[i].size();
  b->avg.resize(nel);
  for (int v=0;v<nel;v++) {
    b->avg[v].c.resize(nc);
    b->avg[v].avg_state = mem.myalloc(b->avg[v].avg_state,
                                            par.nvar+F2C_Nxd);
  }

  // add each cell in the child grid to the "avg" vector:
#ifdef TEST_MPI_NG
  cout <<"F2C_SERIAL: adding cells to avg struct. nel="<<nel<<"\n";
#endif
  add_cells_to_avg(par.ndim,child,nel,b->avg);


  return 0;
}



// ##################################################################
// ##################################################################



void NG_fine_to_coarse_bc::add_cells_to_avg(
      int ndim,                   ///< grid dimension
      class GridBaseClass *grid, ///< pointer to child grid.
      int nel,                    ///< number of coarse cells
      std::vector<struct averaging> &avg ///< avg list
      )
{
#ifdef TEST_MPI_NG
  cout <<"NG_fine_to_coarse_bc::add_cells_to_avg()";
  cout <<" nd="<<ndim<<", grid="<<grid<<", nel="<<nel;
  cout <<", avg.size="<<avg.size()<<"\n";
  cout <<"NG="<<grid->NG(XX)<<"\n";
#endif
  // loop through avg vector and add cells and positions
  cell *f = grid->FirstPt(), *m=f;
  int v=0, ix=0, iy=0, iz=0;
  int ipos[MAX_DIM];
  int dxo2 = grid->idx()/2;
  
  for (v=0;v<nel;v++) {
#ifdef TEST_MPI_NG
    //cout <<"v="<<v<<" working! ix="<<ix<<"\n";
#endif
    // add fine cells to each avg struct
    avg[v].c[0] = f;
    avg[v].c[1] = grid->NextPt(f,XP);
    if (ndim>1) {
      avg[v].c[2] = grid->NextPt(avg[v].c[0],YP);
      avg[v].c[3] = grid->NextPt(avg[v].c[1],YP);
    }
    if (ndim>2) {
      avg[v].c[4] = grid->NextPt(avg[v].c[0],ZP);
      avg[v].c[5] = grid->NextPt(avg[v].c[1],ZP);
      avg[v].c[6] = grid->NextPt(avg[v].c[2],ZP);
      avg[v].c[7] = grid->NextPt(avg[v].c[3],ZP);
    }
    // add position of coarse cell centre (for testing)
    for (int i=0;i<ndim;i++)
      ipos[i] = avg[v].c[0]->pos[i]+dxo2;
    for (int i=ndim;i<MAX_DIM;i++)
      ipos[i] = 0.0;
    CI.get_dpos_vec(ipos,avg[v].cpos);
#ifdef TEST_MPI_NG
    rep.printVec("~AVG cellpos",ipos,ndim);
    for (unsigned int i=0;i<avg[0].c.size();i++) {
      rep.printVec("cellpos",avg[v].c[i]->pos,ndim);
    }
    //rep.printVec("fine cell pos",f->pos,ndim);
#endif
    // get to next cell.
    f = grid->NextPt(f);
    ix++;
    f = grid->NextPt(f);
    ix++;
    if (ix>=grid->NG(XX)) {
      // end of column, loop to next y-column
#ifdef TEST_MPI_NG
      cout <<"eoc: "<<ix<<","<<iy<<","<<iz<<"\n";
#endif
      ix = 0;
      if (ndim>1) {
        iy++;
        if (iy<grid->NG(YY)-1) {
          f = grid->NextPt(f,YP);
          iy++;
#ifdef TEST_MPI_NG
          cout <<"moving to next plane, iy="<<iy<<"\n";
#endif
        }
        else {
          // end of plane, loop to next z-column
#ifdef TEST_MPI_NG
          cout <<"eop: "<<ix<<","<<iy<<","<<iz<<"\n";
#endif
          iy = 0;
          if (ndim>2) {
            iz++;
            m = grid->NextPt(m,ZP);
            f=m;
            if (iz<grid->NG(ZZ)-1) {
              m = grid->NextPt(m,ZP);
              f=m;
              iz++;
            }
            else {
              // must be at end of grid.
#ifdef TEST_MPI_NG
              cout <<"eog: "<<ix<<","<<iy<<","<<iz<<"\n";
#endif
            }
          } // ndim==3
        } // end of plane
      } // ndim>1
    } // end of column
  } // go to next avg element.
  return;
}



// ##################################################################
// ##################################################################



int NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the NG grid structure
      struct boundary_data *b,
      const int i,        ///< which element of NGrecvF2C to use
      const int cstep,
      const int maxstep
      )
{
  //
  // This is relatively straighforward, in that we just weight each
  // fine cell by its volume.  Assume there are two cells in each
  // dimension in the fine grid, and so we can loop over this.
  //
  list<cell*>::iterator c_iter=b->NGrecvF2C[i].begin();
  cell *c;
  class GridBaseClass *fine   = par.levels[level].child;
  if (!fine) return 0; // if there is no child

  // averaged data contains primitive vector plus column densities.
  int nv = par.nvar + F2C_Nxd;
  double cd[nv];
  size_t i_el=0;
  int nc = b->avg[0].c.size();

  for (c_iter =b->NGrecvF2C[i].begin();
       c_iter!=b->NGrecvF2C[i].end();
                              ++c_iter) {
    c = (*c_iter);
    for (int v=0;v<nv;v++) cd[v]=0.0;
    
    average_cells(par, solver, fine, nc,
                  b->avg[i_el].c, b->avg[i_el].cpos, cd);
    
    // set coarse cell optical depths for any radiation sources by
    // taking values from array "cd" (see get_F2C_Tau() for ordering)
    struct rad_src_info *s;
    int off;
    for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
      s = &(par.RS.sources[isrc]);
      off = par.nvar + F2C_tauoff[isrc];
      CI.set_col(c, s->id, &(cd[off]));
      CI.set_cell_col(c, s->id, &(cd[off+s->NTau]));
    }
    
    // set primitive variables in coarse cell based on interpolated
    // values stored in first nvar elements of "cd".
    solver->UtoP(cd,c->Ph,par.EP.MinTemperature,par.gamma);
    //
    // if full step then assign to c->P as well as c->Ph.
    //
    if (cstep==maxstep) {
      for (int v=0;v<par.nvar;v++) c->P[v] = c->Ph[v];
    }
    i_el++;
  }
  return 0;
}





// ##################################################################
// ##################################################################



int NG_fine_to_coarse_bc::average_cells(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      class GridBaseClass *grid, ///< fine-level grid
      const int ncells,  ///< number of fine-level cells
      std::vector<cell *> &c,   ///< list of cells
      const pion_flt *cpos,  ///< centre of coarse cell.
      double *cd       ///< [OUTPUT] averaged data (conserved var).
      )
{
  pion_flt u[par.nvar];
  //
  // loop through list, adding conserved var * cell-vol,
  // then divide by coarse cell vol.
  //
  double sum_vol=0.0, vol=0.0;
  vector<cell*>::iterator c_iter;
  for (c_iter=c.begin(); c_iter!=c.end(); ++c_iter) {
    cell *f = (*c_iter);
#ifdef TEST_MPI_NG
    if (!f) rep.error("cell doesn't exist average_cells",f);
#endif
    //cout <<"cell "<<f->id<<" averaging. ";
    // get conserved vars for cell in fine grid, *cellvol.
    solver->PtoU(f->Ph, u, par.gamma);
    vol = grid->CellVolume(f,0);
    sum_vol += vol;
    //cout <<"vol="<<vol<<", sum="<<sum_vol<<"; ";
    for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
  }
  //cout <<"\n";
  for (int v=0;v<par.nvar;v++) cd[v] /= sum_vol;

  //
  // If doing RT, we also want to update column densities here.
  //
  int pos[MAX_DIM];
  CI.get_ipos_vec(cpos,pos);
  get_F2C_TauAvg(par, ncells, c, pos, &(cd[par.nvar]) );

  return 0;
}



// ##################################################################
// ##################################################################



void NG_fine_to_coarse_bc::get_F2C_TauAvg(
      class SimParams &par,   // pointer to simulation parameters
      const int ncells,       // number of fine-level cells
      std::vector<cell *> &c, // list of cells
      const int *cpos,  // centre of coarse cell (integer coords)
      double *T              // [OUTPUT] pointer to optical depths
      )
{

  if      (par.ndim==1) get_F2C_TauAvg_1D(par,ncells,c,cpos,T);
  else if (par.ndim==2) get_F2C_TauAvg_2D(par,ncells,c,cpos,T);
  else                  get_F2C_TauAvg_3D(par,ncells,c,cpos,T);
  return;
}



// ##################################################################
// ##################################################################


void NG_fine_to_coarse_bc::get_F2C_TauAvg_1D(
      class SimParams &par,   // pointer to simulation parameters
      const int ncells,       // number of fine-level cells
      std::vector<cell *> &c, // list of cells
      const int *cpos,  // centre of coarse cell (integer coords)
      double *T              // [OUTPUT] pointer to optical depths
      )
{
  // data needed for getting fine cell Taus onto coarse grid.
  double Tau1[MAX_TAU], Tau2[MAX_TAU], Tavg[MAX_TAU];
  class cell *f1, *f2;
  struct rad_src_info *s;
  //int diffx,diffy;
  int spos[MAX_DIM];
  for (int iT=0; iT<MAX_TAU; iT++) Tavg[iT] = 0.0;

  f1 = c[0];
  f2 = c[1];

  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    s = &(par.RS.sources[isrc]);
    CI.get_ipos_vec(s->pos,spos);
    CI.get_col(f1, s->id, Tau1);
    CI.get_col(f2, s->id, Tau2);
    // column from source through cell depends on
    // which direction is to the source
    if (spos[XX] > cpos[XX]) {
      for (int iT=0; iT<s->NTau; iT++) Tavg[iT] = Tau1[iT];
    }
    else {
      for (int iT=0; iT<s->NTau; iT++) Tavg[iT] = Tau2[iT];
    }
    // column through cell is sum of two fine cells.
    CI.get_cell_col(f1, s->id, Tau1);
    CI.get_cell_col(f2, s->id, Tau2);
    for (int v=0; v<s->NTau; v++) {
      Tau1[v] += Tau2[v];
    }
    // Tavg[] is col2cell for coarse cell for this source.
    // Tau1[] is cell_col for this source.
    // Copy them into array T.
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+iT] = Tavg[iT];
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+s->NTau+iT] = Tau1[iT];
  }

  return;
}



// ##################################################################
// ##################################################################



void NG_fine_to_coarse_bc::get_F2C_TauAvg_2D(
      class SimParams &par,   // pointer to simulation parameters
      const int ncells,       // number of fine-level cells
      std::vector<cell *> &c, // list of cells
      const int *cpos,  // centre of coarse cell (integer coords)
      double *T              // [OUTPUT] pointer to optical depths
      )
{
  // data needed for getting fine cell Taus onto coarse grid.
  double Tau1[MAX_TAU], Tau2[MAX_TAU], Tau3[MAX_TAU],
         Tau4[MAX_TAU], Tavg[MAX_TAU], dTavg[MAX_TAU];
  class cell *f1, *f2, *f3, *f4;
  struct rad_src_info *s;
  int diffx,diffy;
  int spos[MAX_DIM];
  for (int iT=0; iT<MAX_TAU; iT++) Tavg[iT] = 0.0;
  for (int iT=0; iT<MAX_TAU; iT++) dTavg[iT] = 0.0;

#ifdef RT_TESTING
  cout <<"2D F2C RT Routine\n";
#endif
  f1 = c[0];
  f2 = c[1];
  f3 = c[2];
  f4 = c[3];
#ifdef RT_TESTING
  cout <<"f: ["<<f1->pos[XX]<<","<<f1->pos[YY]<<"], c: [";
  cout <<cpos[XX]<<","<<cpos[YY]<<"] : ";
#endif
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    s = &(par.RS.sources[isrc]);
    CI.get_ipos_vec(s->pos,spos);
    //rep.printVec("spos",spos,2);
    //rep.printVec("cpos",cpos,2);
    CI.get_col(f1, s->id, Tau1);
    CI.get_col(f2, s->id, Tau2);
    CI.get_col(f3, s->id, Tau3);
    CI.get_col(f4, s->id, Tau4);
    diffx = spos[XX] - cpos[XX];
    diffy = spos[YY] - cpos[YY];
#ifdef RT_TESTING
    cout <<"diffxy = "<<diffx<<"  "<<diffy<<": ";
    cout <<"t1="<<*Tau1<<", t2="<<*Tau2<<", t3=";
    cout <<*Tau3<<", t4="<<*Tau4;
#endif
    // column from source through cell depends on
    // which direction is to the source
    if      (diffx>0 && diffx==diffy) {
      // src at 45 deg, take Tau 1.
      for (int v=0; v<s->NTau; v++) Tavg[v] = Tau1[v];
      // dtau = dtau1 + dtau4
      CI.get_cell_col(f1, s->id, Tau1);
      CI.get_cell_col(f4, s->id, Tau4);
      for (int v=0; v<s->NTau; v++) dTavg[v] = Tau1[v]+Tau4[v];
    }
    else if (diffx>0 && diffx==-diffy) {
      // src at -45 deg, take Tau 3.
      for (int v=0; v<s->NTau; v++) Tavg[v] = Tau3[v];
      // dtau = dtau2 + dtau3
      CI.get_cell_col(f2, s->id, Tau2);
      CI.get_cell_col(f3, s->id, Tau3);
      for (int v=0; v<s->NTau; v++) dTavg[v] = Tau2[v]+Tau3[v];
    }
    else if (diffx<0 && diffx==diffy) {
      // src at 225 deg, take Tau 4.
      for (int v=0; v<s->NTau; v++) Tavg[v] = Tau4[v];
      // dtau = dtau1 + dtau4
      CI.get_cell_col(f1, s->id, Tau1);
      CI.get_cell_col(f4, s->id, Tau4);
      for (int v=0; v<s->NTau; v++) dTavg[v] = Tau1[v]+Tau4[v];
    }
    else if (diffx<0 && diffx==-diffy) {
      // src at 135 deg, take Tau 2.
      for (int v=0; v<s->NTau; v++) Tavg[v] = Tau2[v];
      // dtau = dtau2 + dtau3
      CI.get_cell_col(f2, s->id, Tau2);
      CI.get_cell_col(f3, s->id, Tau3);
      for (int v=0; v<s->NTau; v++) dTavg[v] = Tau2[v]+Tau3[v];
    }
    else {
      // not at 45deg to grid, so do simpler averaging.
      if (diffx>0 && fabs(diffx)>=fabs(diffy)) {
        // Source in Q1 coming from dir XP (-45 < theta < 45 deg)
        for (int v=0; v<s->NTau; v++) {
          //Tavg[v] = 0.5*(Tau1[v]+Tau3[v]);
          Tavg[v] = std::max(Tau1[v],Tau3[v]);
        }
      }
      else if (diffy>0 && fabs(diffx)<fabs(diffy)) {
        // source in Q2, coming from dir YP (45 < theta < 135 deg)
        for (int v=0; v<s->NTau; v++) {
          //Tavg[v] = 0.5*(Tau1[v]+Tau2[v]);
          Tavg[v] = std::max(Tau1[v],Tau2[v]);
        }
      }
      else if (diffx<0 && fabs(diffx)>=fabs(diffy)) {
        // source in Q3, coming from XN (135 < theta < 225 deg)
        for (int v=0; v<s->NTau; v++) {
          //Tavg[v] = 0.5*(Tau2[v]+Tau4[v]);
          Tavg[v] = std::max(Tau2[v],Tau4[v]);
        }
      }
      else {
        // source in Q4, coming from YN (225 < theta < 315 deg)
        for (int v=0; v<s->NTau; v++) {
          //Tavg[v] = 0.5*(Tau3[v]+Tau4[v]);
          Tavg[v] = std::max(Tau3[v],Tau4[v]);
        }
      }
#ifdef RT_TESTING
      cout <<"  tc="<<*Tavg;
#endif
      // column through cell is sum of 4 fine cells divided by 2.
      // this is a fairly crude approx, and so it is checked in the
      // next code block.
      CI.get_cell_col(f1, s->id, Tau1);
      CI.get_cell_col(f2, s->id, Tau2);
      CI.get_cell_col(f3, s->id, Tau3);
      CI.get_cell_col(f4, s->id, Tau4);
      for (int v=0; v<s->NTau; v++) {
        dTavg[v] = 0.5*(Tau1[v] + Tau2[v] + Tau3[v] + Tau4[v]);
      }
    }

    // fix dTau if dTau>Tau
    if (dTavg[0]>Tavg[0]*ONE_PLUS_EPS) {
      //cout <<"dTau is bigger than Tau: "<<dTavg[0]<<"  "<<Tavg[0];
      //cout <<"  "<<diffx<<"  "<<diffy;
      //cout <<"  "<<*Tau1<<"  "<<*Tau2<<"  "<<*Tau3<<"  "<<*Tau4<<"\n";
      //CI.print_cell(f1);
      //CI.print_cell(f2);
      //CI.print_cell(f3);
      //CI.print_cell(f4);
      //rep.error("interpolation",99);
      for (int iT=0; iT<s->NTau; iT++) dTavg[iT] = 0.99*Tavg[iT];
    }
      
    // Tavg[] is col2cell for coarse cell for this source.
    // dTavg[] is cell_col for this source.
    // Copy them into array T.
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+iT] = Tavg[iT];
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+s->NTau+iT] = dTavg[iT];

#ifdef RT_TESTING
    cout <<"  dtc="<<*Tau1 <<"\n";
#endif
  }
  //rep.printVec("*** T ***",T,F2C_Nxd);

  return;
}


// ##################################################################
// ##################################################################


void NG_fine_to_coarse_bc::get_F2C_TauAvg_3D(
      class SimParams &par,   // pointer to simulation parameters
      const int ncells,       // number of fine-level cells
      std::vector<cell *> &c, // list of cells
      const int *cpos,    // centre of coarse cell (integer coords)
      double *T              // [OUTPUT] pointer to optical depths
      )
{
  // data needed for getting fine cell Taus onto coarse grid.
  double Tau[8][MAX_TAU], Tavg[MAX_TAU], dTavg[MAX_TAU];
  class cell *f[8]; // list of fine cells
  struct rad_src_info *s;
  //int diffx,diffy,diffz;
  int spos[MAX_DIM];
  for (int iT=0; iT<MAX_TAU; iT++) Tavg[iT] = 0.0;
  for (int iT=0; iT<MAX_TAU; iT++) dTavg[iT] = 0.0;

#ifdef RT_TESTING
  cout <<"3D F2C RT Routine\n";
#endif
  for (int i=0;i<8;i++) f[i] = c[i];

  // loop over sources
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    s = &(par.RS.sources[isrc]);
    CI.get_ipos_vec(s->pos,spos);
    //rep.printVec("spos",spos,2);
    //rep.printVec("cpos",cpos,2);
    // Get "Tau to cell" for each fine cell.
    for (int i=0;i<8;i++) CI.get_col(f[i], s->id, Tau[i]);
    // Get offsets of source from coarse cell centre.
    //diffx = spos[XX] - cpos[XX];
    //diffy = spos[YY] - cpos[YY];
    //diffz = spos[ZZ] - cpos[ZZ];
    
    // get Tau and dTau for coarse cell.
    for (int v=0; v<s->NTau; v++)  Tavg[v]=0.0;
    for (int v=0; v<s->NTau; v++) dTavg[v]=0.0;

    // Set Tau to be the max of all the fine cells.
    // This is a small overestimate for some cells, but is only a
    // significant error near a source (where there is a finer grid
    // anyway).
    for (int i=0;i<8;i++) {
      for (int v=0; v<s->NTau; v++) {
        Tavg[v] = std::max(Tavg[v],Tau[i][v]);
      }
    }

    // column through cell is sum of 8 fine cells divided by 4.
    // this is a fairly crude approx, and so it is checked and
    // corrected if needed.
    for (int i=0;i<8;i++) CI.get_cell_col(f[i], s->id, Tau[i]);
    for (int v=0; v<s->NTau; v++) {
      for (int i=0;i<8;i++) {
        dTavg[v] += Tau[i][v];
      }
      dTavg[v] *= 0.25;
      if (dTavg[v]>Tavg[v]) dTavg[v] = 0.9999*Tavg[v];
    }

    // Tavg[] is col2cell for coarse cell for this source.
    // dTavg[] is cell_col for this source.
    // Copy them into array T.
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+iT] = Tavg[iT];
    for (int iT=0; iT<s->NTau; iT++)
      T[F2C_tauoff[isrc]+s->NTau+iT] = dTavg[iT];

#ifdef RT_TESTING
    cout <<"  3D  dtc="<<*Tau1 <<"\n";
#endif
  }
  //rep.printVec("*** T ***",T,F2C_Nxd);

  return;
}




// ##################################################################
// ##################################################################



