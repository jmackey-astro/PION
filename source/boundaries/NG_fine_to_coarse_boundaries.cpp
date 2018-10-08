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
  //cell *cc = child->FirstPt_All(); // child cell.
  //int cdx = 0.5*child->idx();
  int nc=1;  // number of fine cells per coarse cell
  for (int i=0;i<par.ndim;i++) nc*=2;
  int nel = b->NGrecvF2C[i].size();
  b->avg.resize(nel);
  for (int v=0;v<nel;v++) {
    b->avg[v].c.resize(nc);
    b->avg[v].avg_state = mem.myalloc(b->avg[v].avg_state,
                                                par.nvar);
  }

  // add each cell in the child grid to the "avg" vector:
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
  cell *f = grid->FirstPt();
  int v=0, ix=0, iy=0, iz=0;
  int ipos[MAX_DIM];
  double dxo2 = 0.5*grid->DX();
  
  for (v=0;v<nel;v++) {
#ifdef TEST_MPI_NG
    cout <<"v="<<v<<" working! ix="<<ix<<"\n";
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
    for (unsigned int i=0;i<avg[0].c.size();i++) {
      rep.printVec("cellpos",avg[v].c[i]->pos,ndim);
    }
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
      //f = grid->NextPt(f);
      ix = 0;
      if (ndim>1) {
        iy++;
        if (iy<grid->NG(YY)) {
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
            if (iz<grid->NG(ZZ)) {
              f = grid->NextPt(f,ZP);
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
  double cd[par.nvar];
  size_t i_el=0;
  int nc = b->avg[0].c.size();

  for (c_iter =b->NGrecvF2C[i].begin();
       c_iter!=b->NGrecvF2C[i].end();
                              ++c_iter) {
    c = (*c_iter);
    for (int v=0;v<par.nvar;v++) cd[v]=0.0;
    
    average_cells(par,solver,fine,nc,b->avg[i_el].c,cd);
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
      pion_flt *cd       ///< [OUTPUT] averaged data (conserved var).
      )
{
  pion_flt u[par.nvar];
  //
  // simple: loop through list, adding conserved var * cell-vol,
  // then divide by coarse cell vol.
  //
  double sum_vol=0.0, vol=0.0;
  vector<cell*>::iterator c_iter;
  for (c_iter=c.begin(); c_iter!=c.end(); ++c_iter) {
    cell *f = (*c_iter);
#ifdef TEST_MPI_NG
    if (!f) rep.error("cell doesn't exist average_cells",f);
#endif
    // get conserved vars for cell in fine grid, *cellvol.
    solver->PtoU(f->Ph, u, par.gamma);
    vol = grid->CellVolume(f);
    sum_vol += vol;
    for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
  }
  for (int v=0;v<par.nvar;v++) cd[v] /= sum_vol;
  return 0;
}



// ##################################################################
// ##################################################################



