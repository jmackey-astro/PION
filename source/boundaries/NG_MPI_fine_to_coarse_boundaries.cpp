/// \file NG_fine_to_coarse_boundaries_MPI.cpp
/// \brief Class definitions for NG_fine_to_coarse_MPI boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/NG_fine_to_coarse_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################



int NG_fine_to_coarse_MPI_bc::BC_assign_FINE_TO_COARSE_SEND(
      class SimParams &par,     ///< pointer to simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b  ///< boundary data
      //class GridBaseClass *child  ///< pointer to child grid.
      )
{
  // Check if parent grid is on my MPI process
  int pproc = par.level[l].MCMD.parent_proc;
  int err=0;
  class MCMDcontrol *MCMD = &(par.level[l].MCMD);

  // If so, do nothing because the parent grid will call the 
  // serial NG setup function, which just grabs the data from 
  // this grid.
  if (MCMD->get_myrank() == pproc) {
#ifdef TEST_MPI_NG
    cout <<"my rank == parent rank, not setting up ";
    cout <<"FINE_TO_COARSE_SEND\n";
#endif
  }
  else {
    // Else, make a vector of structs with a list of cells contained
    // in each coarse cell (2,4,or 8) of the parent grid.  This grid
    // is (by design) entirely contained within the parent.
    class GridBaseClass *grid = par.level[l].grid;
    double dxo2 = 0.5*grid->DX();
    int ipos[MAX_DIM];
    int nc=1;  // number of fine cells per coarse cell
    for (int i=0;i<par.ndim;i++) nc*=2;
    int nel = grid->Ncell()/nc; // number of coarse cells
    // initialise the "avg" vector to hold the average state.
    b->avg.resize(nel);
    for (int v=0;v<nel;v++) {
      b->avg[v]->c.resize(nc);
      b->avg[v]->state = mem.myalloc(b->avg[v]->state,par.nvar);
    }
    // loop through avg vector and add cells and positions
    cell *f = grid->FirstPt();
    int v=0, ix=0, iy=0, iz=0;
    for (v=0;v<nel;v++) {
      // add fine cells to each avg struct
      b->avg[v]->c[0] = f;
      b->avg[v]->c[1] = grid->NextPt(f,XP);
      if (ndim>1) {
        b->avg[v]->c[2] = grid->NextPt(b->avg[v]->c[0],YP);
        b->avg[v]->c[3] = grid->NextPt(b->avg[v]->c[1],YP);
      }
      if (ndim>2) {
        b->avg[v]->c[4] = grid->NextPt(b->avg[v]->c[0],ZP);
        b->avg[v]->c[5] = grid->NextPt(b->avg[v]->c[1],ZP);
        b->avg[v]->c[6] = grid->NextPt(b->avg[v]->c[2],ZP);
        b->avg[v]->c[7] = grid->NextPt(b->avg[v]->c[3],ZP);
      }
      // add position of coarse cell centre (for testing)
      for (int i=0;i<par.ndim;i++)
        ipos[v] = b->avg[v]->c[0]->pos[v]+dxo2;
      for (int i=par.ndim;i<MAX_DIM;i++)
        ipos[v] = 0.0;
      CI.get_dpos(ipos,b->avg[v]->cpos);

      // get to next cell.
      f = grid->NextPt(f);
      ix++;
      if (ix<grid->NG(XX)) {
        f = grid->NextPt(f);
        ix++;
      }
      else {
        // end of column, loop to next y-column
#ifdef TEST_MPI_NG
        cout <<"eoc: "<<ix<<","<<iy<<","<<iz<<"\n";
#endif
        f = grid->NextPt(f);
        ix = 0;
        if (ndim>1) {
          iy++;
          if (iy<grid->NG(YY)) {
            f = grid->NextPt(f,YP);
            iy++;
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
  } // if parent is not on a different MPI process

  return 0;
}



// ##################################################################
// ##################################################################



int NG_fine_to_coarse_MPI_bc::BC_assign_FINE_TO_COARSE_RECV(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b  ///< boundary data
      //class GridBaseClass *parent  ///< pointer to parent grid.
      )
{
  return 0;
}



// ##################################################################
// ##################################################################



int NG_fine_to_coarse_MPI_bc::BC_update_FINE_TO_COARSE_SEND(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int level, ///< level in the NG grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{

  // Check if parent grid is on my MPI process
  int pproc = par.level[l].MCMD.parent_proc;
  int err=0;
  class MCMDcontrol *MCMD = &(par.level[l].MCMD);

  // If so, do nothing because the parent grid will call the 
  // serial NG setup function, which just grabs the data from 
  // this grid.
  if (MCMD->get_myrank() == pproc) {
#ifdef TEST_MPI_NG
    cout <<"my rank == parent rank, no need to send anything ";
    cout <<"FINE_TO_COARSE_SEND\n";
#endif
    //NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(par,solver,level,b,cstep,maxstep);
  }
  else {
    class GridBaseClass *grid = par.level[l].grid;
    double dxo2 = 0.5*grid->DX();
    int nc = b->avg[0]->c.size();
    int nel = b->avg.size();

    // loop through avg vector and add cells and positions
    int v=0, ix=0, iy=0, iz=0;
    double cd[par.nvar];
    for (v=0;v<nel;v++) {

      for (int v=0;v<par.nvar;v++) cd[v]=0.0;
      average_cells(par,solver,grid,nc,b->avg[v]->c,cd);
    } // go to next avg element.

    // pack data to send
    
  } // if parent is not on a different MPI process


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


void 

