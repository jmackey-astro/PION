/// \file SMR_fine_to_coarse_boundaries.cpp
/// \brief Class definitions for SMR_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/SMR_fine_to_coarse_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################



int SMR_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE(
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
  b->SMR.clear();

  list<cell*>::iterator bpt=b->data.begin();
  cell *cc = child->FirstPt_All(); // child cell.
  int cdx = 0.5*child->idx();

  // Map each bpt cell to a cell in b->SMR list, which is the first
  // cell in the finer grid that is part of the coarse cell (i.e. the
  // one with the most negative coordinates).
  do{
    cc = child->FirstPt_All();
    for (int v=0;v<par.ndim;v++) {
      while (cc && cc->pos[v] < (*bpt)->pos[v]-cdx)
        cc = child->NextPt(cc,static_cast<direction>(2*v+1));
    }
    if (!cc) rep.error("BC_assign_FINE_TO_COARSE: lost on fine grid",0);
    b->SMR.push_back(cc);
    
    ++bpt;
  }  while (bpt !=b->data.end());

  return 0;
}



// ##################################################################
// ##################################################################



int SMR_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the SMR grid structure
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
  list<cell*>::iterator f_iter=b->SMR.begin();
  cell *c, *f;
  //cout <<"Fine to Coarse boundary update (internal).  list sizes: ";
  //cout << b->data.size() <<",  :"<<b->SMR.size()<<"\n";

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
    
    // 1D
    // get conserved vars for cell 1 in fine grid, multiply by cell
    // volume.
    //
    solver->PtoU(f->Ph, u, par.gamma);
    vol = fine->CellVolume(f);
    for (int v=0;v<par.nvar;v++) cd[v] = u[v]*vol;
#ifdef TEST_NEST
    //vol_sum += vol;
#endif

    // get conserved vars for cell 2 in fine grid, *cellvol.
    f = fine->NextPt(f,XP);
    solver->PtoU(f->Ph, u, par.gamma);
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
      solver->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 3 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 4 in fine grid, *cellvol.
      f = fine->NextPt(f,XP);
      solver->PtoU(f->Ph, u, par.gamma);
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
      solver->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 5 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 6 in fine grid, *cellvol.
      f = fine->NextPt(f,XP);
      solver->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 6 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 7 in fine grid, *cellvol.
      f = fine->NextPt(f,YP);
      solver->PtoU(f->Ph, u, par.gamma);
      vol = fine->CellVolume(f);
      for (int v=0;v<par.nvar;v++) cd[v] += u[v]*vol;
#ifdef TEST_NEST
      //CI.get_ipos(f,fpos);
      //rep.printVec("fine 7 pos:",fpos,par.ndim);
      //vol_sum += vol;
#endif

      // get conserved vars for cell 8 in fine grid, *cellvol.
      f = fine->NextPt(f,XN);
      solver->PtoU(f->Ph, u, par.gamma);
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




