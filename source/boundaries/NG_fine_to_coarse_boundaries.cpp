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
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,  ///< boundary data
      class GridBaseClass *child  ///< pointer to child grid.
      )
{
  //
  // Make a list of child-grid cells to map onto the coarse grid
  //
  if (b->data.empty())
    rep.error("BC_assign_FINE_TO_COARSE: empty boundary data",b->itype);
  b->NG.clear();

  list<cell*>::iterator bpt=b->data.begin();
  cell *cc = child->FirstPt_All(); // child cell.
  int cdx = 0.5*child->idx();

  // Map each bpt cell to a cell in b->NG list, which is the first
  // cell in the finer grid that is part of the coarse cell (i.e. the
  // one with the most negative coordinates).
  do{
    cc = child->FirstPt_All();
    for (int v=0;v<par.ndim;v++) {
      while (cc && cc->pos[v] < (*bpt)->pos[v]-cdx)
        cc = child->NextPt(cc,static_cast<direction>(2*v+1));
    }
    if (!cc) rep.error("BC_assign_FINE_TO_COARSE: lost on fine grid",0);
    b->NG.push_back(cc);
    
    ++bpt;
  }  while (bpt !=b->data.end());

  return 0;
}



// ##################################################################
// ##################################################################



int NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the NG grid structure
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
  list<cell*>::iterator f_iter=b->NG.begin();
  cell *c, *f;
  //cout <<"Fine to Coarse boundary update (internal).  list sizes: ";
  //cout << b->data.size() <<",  :"<<b->NG.size()<<"\n";

  // pointers to coarse and fine grids:
  class GridBaseClass *coarse = par.levels[level].grid;
  class GridBaseClass *fine   = par.levels[level].child;

  // vol_sum is for testing only (make sure that fine grid cells
  // have the same cumulative volume as the coarse cell).
  double cd[par.nvar], u[par.nvar], vol=0.0;
#ifdef TEST_NEST
  double vol_sum=0.0;
  int cpos[MAX_DIM];
  int fpos[MAX_DIM];
#endif

  for (c_iter=b->data.begin(); c_iter!=b->data.end(); ++c_iter) {
    c = (*c_iter);
    f = (*f_iter);
    vol=0.0;

#ifdef TEST_NEST
    vol_sum=0.0;
    //CI.get_ipos(c,cpos);
    //CI.get_ipos(f,fpos);
    //rep.printVec("coarse pos:",cpos,par.ndim);
    //rep.printVec("fine 1 pos:",fpos,par.ndim);
#endif
    
    int nc = 1;
    for (int i=0;i<par.ndim;i++) nc*=2;
    list<cell *> c;

    c.push_back(f);
    c.push_back(fine->NextPt(f,XP));
    if (ndim>1) {
      c.push_back(fine->NextPt(f,YP));
      c.push_back(fine->NextPt(fine->NextPt(f,XP),YP));
    }
    if (ndim>2) {
      f=fine->NextPt(f,ZP);
      c.push_back(f);
      c.push_back(fine->NextPt(f,XP));
      c.push_back(fine->NextPt(f,YP));
      c.push_back(fine->NextPt(fine->NextPt(f,XP),YP));
    }

    for (int v=0;v<par.nvar;v++) cd[v]=0.0;
    
    average_cells(par,solver,fine,nc,c,cd);

    //vol = coarse->CellVolume(c);
    //for (int v=0;v<par.nvar;v++) cd[v] /= vol;
    solver->UtoP(cd,c->Ph,par.EP.MinTemperature,par.gamma);
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



int NG_fine_to_coarse_bc::average_cells(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      class GridBaseClass *grid, ///< fine-level grid
      const int ncells,  ///< number of fine-level cells
      list<cell *> &c,   ///< list of cells
      pion_flt *cd       ///< [OUTPUT] averaged data (conserved var).
      )
{
  pion_flt u[par.nvar];
  //
  // simple: loop through list, adding conserved var * cell-vol,
  // then divide by coarse cell vol.
  //
  double sum_vol=0.0;
  list<cell*>::iterator c_iter;
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

