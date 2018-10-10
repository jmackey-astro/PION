/// \file NG_MPI_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_MPI_coarse_to_fine boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.09.06 JM: started writing code.


#include "boundaries/NG_MPI_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
#include <sstream>
using namespace std;

//#define TEST_C2F

// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_SEND(
      class SimParams &par,  ///< simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b  ///< boundary data
      )
{

  class GridBaseClass *grid = par.levels[l].grid;
  //int gidx = grid->idx();

  // see how many child grids I have
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int nchild = MCMD->child_procs.size();
  b->NGsendC2F.clear();

#ifdef TEST_MPI_NG
  if (nchild==0) {
    cout <<"COARSE_TO_FINE_SEND: no children.\n";
  }
#endif

  // loop over child grids
  for (int i=0;i<nchild;i++) {

    if (MCMD->get_myrank() == MCMD->child_procs[i].rank) {
      // if child is on my process, do nothing because child grid
      // can grab the data directly.
#ifdef TEST_MPI_NG
      cout <<"C2F_SEND: child "<<i<<", ";
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
      cout <<"C2F_SEND: child "<<i<<", ";
      cout <<"my rank != child rank, running parallel ";
      cout <<"COARSE_TO_FINE_SEND\n";
#endif
      // get dimensions of child grid from struct
      int ixmin[MAX_DIM], ixmax[MAX_DIM];
      CI.get_ipos_vec(MCMD->child_procs[i].Xmin, ixmin);
      CI.get_ipos_vec(MCMD->child_procs[i].Xmax, ixmax);

      // loop over dimensions
      for (int d=0;d<par.ndim;d++) {
        // if child xmin == its level xmin, but > my level xmin,
        // then we need to send data, so set up a list.
        if ( (pconst.equalD(MCMD->child_procs[i].Xmin[d],
                            par.levels[l+1].Xmin[d]))        &&
             (MCMD->child_procs[i].Xmin[d] > 
              par.levels[l].Xmin[d]*ONE_PLUS_EPS) ) {
          cout <<"C2F_SEND: child "<<i<<", dim "<<d<<" NEG DIR\n";
          struct c2f *bdata = new struct c2f;
          bdata->rank = MCMD->child_procs[i].rank;
          bdata->dir  = 2*d;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(par,grid,bdata,ixmin,ixmax);
          b->NGsendC2F.push_back(bdata);
        }
        // if child xmax == its level xmax, but < my level xmax,
        // then we need to send data, so set up a list.
        else if ((pconst.equalD(MCMD->child_procs[i].Xmax[d],
                                par.levels[l+1].Xmax[d]))    &&
                 (MCMD->child_procs[i].Xmax[d] < 
                  par.levels[l].Xmax[d]*ONE_MINUS_EPS) ) {
          cout <<"C2F_SEND: child "<<i<<", dim "<<d<<" POS DIR\n";
          struct c2f *bdata = new struct c2f;
          bdata->rank = MCMD->child_procs[i].rank;
          bdata->dir  = 2*d+1;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(par,grid,bdata,ixmin,ixmax);
          b->NGsendC2F.push_back(bdata);
        }
      } // loop over dimensions
    } // if child is not on my process
  } // loop over child grids
  return 0;
}



// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_SEND(
      class SimParams &par,         ///< simulation parameters
      class GridBaseClass *grid,    ///< pointer to coarse-level grid
      class FV_solver_base *solver, ///< pointer to equations
      const int l,                  ///< level in the NG hierarchy
      struct boundary_data *b,      ///< pointer to boundary struct
      const int cstep,              ///< fractional step
      const int maxstep             ///< number of fractional steps
      )
{
#ifdef TEST_C2F
  cout <<"MPI C2F SEND: starting... ";
  if (b->NGsendC2F.size()==0) {
    cout <<"empty send list, so just returning now.\n";
    return 0;
  }
  cout <<endl;
#endif
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int err=0;
  //
  // if on an odd-numbered step, then need to update the data on the
  // coarse grid to half way through a coarse step.  Assume dU has
  // been calculated for the coarse grid, but not updated to Ph[]
  //
  if (cstep != maxstep) {
#ifdef TEST_C2F
    cout <<"MPI C2F SEND: odd step, interpolating data in time.\n";
#endif
    double U[par.nvar];
    for (unsigned int ib=0; ib<b->NGsendC2F.size(); ib++) {
      for (unsigned int c_iter=0; c_iter<b->NGsendC2F[ib]->c.size();
                                                          c_iter++) {
        cell *c = b->NGsendC2F[ib]->c[c_iter];
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
      } // loop over cells in boundary
    } // loop over send boundaries
  } // if not at a full step update

  // loop over send boundaries, pack and send the data.
  for (unsigned int ib=0; ib<b->NGsendC2F.size(); ib++) {
    //
    // if 1st order accuracy, then just need Ph[]+cell-vol.
    // if 2nd order accuracy, also a slope vector for each dimension
    //
    size_t n_cell = b->NGsendC2F[ib]->c.size();
    size_t n_el = 0;
    if      (par.spOOA == OA1) n_el = n_cell*(par.nvar+1+par.ndim);
    else if (par.spOOA == OA2) n_el = n_cell*(3*par.nvar+1+par.ndim);
    else rep.error("bad spOOA in MPI C2F",par.spOOA);
    pion_flt *buf = new pion_flt [n_el];
    double slope[par.nvar];

    // loop over cells, add Ph[], cell-vol, slopes to send buffer
    size_t ibuf=0;
    for (unsigned int c_iter=0; c_iter<b->NGsendC2F[ib]->c.size();
                                                        c_iter++) {
      cell *c = b->NGsendC2F[ib]->c[c_iter];
      for (int v=0;v<par.nvar;v++) buf[ibuf+v]= c->Ph[v];
      ibuf += par.nvar;
      buf[ibuf] = grid->CellVolume(c);
      ibuf++;
      for (int v=0;v<par.ndim;v++)
        buf[ibuf+v]= static_cast<double>(c->pos[v]);
      ibuf += par.ndim;

      if (par.spOOA == OA2) {
        for (int idim=0;idim<par.ndim; idim++) {
          enum axes a = static_cast<axes>(idim);
          solver->SetSlope(c,a,par.nvar,slope,OA2,grid);
          for (int v=0;v<par.nvar;v++) buf[ibuf+v]= slope[v];
          ibuf += par.nvar;
        } // idim
      } // if 2nd order
    } // loop over cells for this send boundary

    if (ibuf != n_el) rep.error("C2F MPI SEND counting",ibuf);

    //
    // Send data using a non-blocking MPI send
    //
    ostringstream tmp;
    tmp <<"C2F_"<<MCMD->get_myrank()<<"_to_"<<b->NGsendC2F[ib]->rank;
    string id = tmp.str();
#ifdef TEST_MPI_NG
    cout <<"BC_update_COARSE_TO_FINE_SEND: Sending "<<n_el;
    cout <<" doubles from proc "<<MCMD->get_myrank();
    cout <<" to child proc "<<b->NGsendC2F[ib]->rank<<"\n";
#endif
    err += COMM->send_double_data(b->NGsendC2F[ib]->rank,n_el,buf,
                                  id,BC_MPI_NGC2F_tag);
    if (err) rep.error("Send_C2F send_data failed.",err);
#ifdef TEST_MPI_NG
    cout <<"BC_update_COARSE_TO_FINE_SEND: returned with id="<<id;
    cout <<"\n";
#endif
    // store ID to clear the send later (and delete the MPI temp data)
    NG_C2F_send_list.push_back(id);
    buf = mem.myfree(buf);
  } // loop over send boundaries
  return 0;
}



// ##################################################################
// ##################################################################



void NG_MPI_coarse_to_fine_bc::BC_COARSE_TO_FINE_SEND_clear_sends()
{
  for (unsigned int ib=0; ib<NG_C2F_send_list.size(); ib++) {
#ifdef TEST_MPI_NG
    cout <<"C2F_send: clearing send # "<<ib<<" of ";
    cout <<NG_C2F_send_list.size()<<"\n";
    cout.flush();
#endif
    COMM->wait_for_send_to_finish(NG_C2F_send_list[ib]);
  }
  NG_C2F_send_list.clear();
  return;
}



// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_RECV(
      class SimParams &par,   ///< simulation parameters
      const int l,            ///< level of this grid.
      boundary_data *b        ///< boundary data
      )
{
  // The boundary is already a regular external boundary, so just
  // need to decide if parent is on the same MPI process (and call
  // the serial version) or on a different MPI process (and set up
  // some data structures to receive the coarse-grid data for
  // interpolation).
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  class GridBaseClass *grid = par.levels[l].grid;
  int pproc = MCMD->parent_proc;

  if (MCMD->get_myrank() == pproc) {
#ifdef TEST_MPI_NG
    cout <<"my rank == parent rank, setting up serial ";
    cout <<"COARSE_TO_FINE_RECV\n";
#endif
    int err = NG_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
                par,grid,b,par.levels[l].parent);
    rep.errorTest("serial C2F BC setup",0,err);
  }
  else {
#ifdef TEST_MPI_NG
    cout <<"my rank != parent rank, so will need MPI for ";
    cout <<"COARSE_TO_FINE_RECV";
    cout <<"me="<<MCMD->get_myrank();
    cout <<", parent="<<pproc<<"\n";
#endif

    //
    // group cells into lists that are child cells of a coarse cell
    //
    size_t n_cell = b->data.size();
    for (int idim=0;idim<par.ndim;idim++) n_cell /=2;
    b->NGrecvC2F.resize(n_cell);
    size_t ic = 0;
    list<cell*>::iterator f_iter=b->data.begin();

    if      (par.ndim==1) {
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        cell *c = (*f_iter);
        b->NGrecvC2F[ic].push_back(c);
        f_iter++;
        c = (*f_iter);
        b->NGrecvC2F[ic].push_back(c);
        ic++;
      } // loop over cells
    } // if 1D
    
    else if (par.ndim==2) {
      // 4 fine cells per coarse cell
      f_iter=b->data.begin();
      cell *c = (*f_iter), *temp=0;
      int row_y = c->pos[YY];
      int idx = grid->idx();
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        c = (*f_iter);
        if (c->pos[YY]-row_y == 2*idx) {
          // move to next row of cells
          row_y = c->pos[YY];
        }
        if (c->pos[YY] == row_y) {
          // on same row of cells so continue adding cells
          b->NGrecvC2F[ic].push_back(c);
          f_iter++;
          temp = (*f_iter);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(c,YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp,YP));
          ic++;
        }
        else if (c->pos[YY]-row_y == idx) {
          // odd-numbered row, already dealt with.
          continue;
        }
        else
          rep.error("error in 3d logic C2FR_setup",c->pos[YY]-row_y);
      } // loop over cells
    } // if 2D

    else if (par.ndim==3) {
      // 8 fine cells per coarse cell
      f_iter=b->data.begin();
      cell *c = (*f_iter), *temp=0;
      int row_y = c->pos[YY];
      int row_z = c->pos[ZZ];
      int idx = grid->idx();
      for (f_iter=b->data.begin(); f_iter!=b->data.end(); ++f_iter) {
        c = (*f_iter);
        if (c->pos[ZZ]-row_z == 2*idx) {
          // move to next plane of cells
          row_z = c->pos[ZZ];
        }
        if (c->pos[YY]-row_y == 2*idx) {
          // move to next row of cells in z-plane
          row_y = c->pos[YY];
        }
        if (c->pos[YY] == row_y && c->pos[ZZ] == row_z) {
          // on same row of cells so continue adding cells
          f_iter++;
          temp = (*f_iter);
          b->NGrecvC2F[ic].push_back(c);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(c,YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp,YP));
          c = grid->NextPt(c,ZP);
          temp = grid->NextPt(temp,YP);
          b->NGrecvC2F[ic].push_back(c);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(c,YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp,YP));
          ic++;
        }
        else if (c->pos[YY]-row_y == idx ||
                 c->pos[ZZ]-row_z == idx) {
          // odd-numbered row/plane, already dealt with.
          continue;
        }
        else
          rep.error("error in 3d logic C2FR_setup",c->pos[YY]-row_y);
      } // loop over cells
    } // if 3D
  } // if parent is on other MPI process

  return 0;
}




// ##################################################################
// ##################################################################



int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_RECV(
      class SimParams &par,      ///< simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int l,    ///< level in the NG grid structure
      struct boundary_data *b,  ///< boundary to update
      const int step  ///< timestep on this (fine) grid
      )
{
#ifdef TEST_C2F
  cout <<"C2F_MPI: receiving boundary data to level ";
  cout <<l<<"\n";
#endif
  int err=0;
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int pproc = MCMD->parent_proc;
  class GridBaseClass *grid   = par.levels[l].grid;

  if (MCMD->get_myrank() == pproc) {
#ifdef TEST_MPI_NG
    cout <<"my rank == parent rank, calling serial ";
    cout <<"COARSE_TO_FINE\n";
#endif
    NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(
                                              par,solver,l,b,step);
  }
  else {
#ifdef TEST_MPI_NG
    cout <<"my rank != parent rank, so updating ";
    cout <<"COARSE_TO_FINE_RECV";
    cout <<"me="<<MCMD->get_myrank();
    cout <<", parent="<<pproc<<"\n";
#endif

    //
    // receive data.
    //
    string recv_id; int recv_tag=-1; int from_rank=-1;
    err = COMM->look_for_data_to_receive(&from_rank, recv_id,
                        &recv_tag,BC_MPI_NGC2F_tag, COMM_DOUBLEDATA);
    if (err) rep.error("look for double data failed",err);
#ifdef TEST_MPI_NG
    cout <<"BC_update_COARSE_TO_FINE_RECV: found data from rank ";
    cout <<from_rank<<", with tag "<< recv_tag<<"\n";
#endif 

    // receive the data: nel is the number of coarse grid cells
    size_t n_cell = b->data.size();
    for (int idim=0;idim<par.ndim;idim++) n_cell /=2;
    size_t n_el = 0;
    if      (par.spOOA == OA1) n_el = n_cell*(par.nvar +1+par.ndim);
    else if (par.spOOA == OA2) n_el = n_cell*(3*par.nvar+1+par.ndim);
    else rep.error("bad spOOA in MPI C2F",par.spOOA);
    pion_flt *buf = 0;
    buf = mem.myalloc(buf,n_el);
#ifdef TEST_MPI_NG
    cout <<"BC_update_COARSE_TO_FINE_RECV: get "<<n_cell;
    cout <<" cells, and "<<n_el<<" doubles.\n";
#endif 
    //
    // Receive data into buffer.  Data stored for each coarse cell:
    // Ph[nv],cellvol,cellpos[nd],slopeX[nv],slopeY[nv],slopeZ[nv]
    //
    err = COMM->receive_double_data(
                              from_rank, recv_tag, recv_id, n_el, buf);
    if (err) rep.error("(BC_update_C2F_RECV) getdata failed",err);

    //
    // Get Ph and slopes from coarse cells, and interpolate into
    // groups of fine cells that are children of each coarse cell.
    //
    size_t ibuf=0;
    if (par.spOOA==OA1) {
      // 1st order, so no averaging.  Fine cells have exactly the
      // same state as the coarse one.
      double Ph[par.nvar];
#ifdef TEST_MPI_NG
      double cpos[par.ndim];
#endif
      //double c_vol=0.0;
      cell *c=0;
      for (unsigned int ic=0; ic<b->NGrecvC2F.size(); ic++) {
        // read data for this coarse cell into arrays
        for (int v=0;v<par.nvar;v++) Ph[v] = buf[ibuf+v];
        ibuf += par.nvar;
        //c_vol = buf[ibuf];
        ibuf++;
        for (int v=0;v<par.ndim;v++) cpos[v] = buf[ibuf+v];
        ibuf += par.ndim;

        list<cell*>::iterator f_iter=b->NGrecvC2F[ic].begin();
        for (f_iter=b->NGrecvC2F[ic].begin();
              f_iter!=b->NGrecvC2F[ic].end(); ++f_iter) {
          c = (*f_iter);
#ifdef TEST_MPI_NG
          for (int v=0;v<par.ndim;v++) cpos[v] -= c->pos[v];
          cout <<"ic="<<ic<<", cell is "<<c<<"  ";
          rep.printVec("offset is:",cpos,par.ndim);
#endif
          for (int v=0;v<par.nvar;v++) c->Ph[v] = Ph[v];
        } // loop over fine cells
      } // loop over coarse cells
    } // if 1st order accurate

    else if (par.spOOA==OA2) {
      // 2nd order needs a slope for each dimension, so split the 
      // put different ndim in if/else statements.
      if (par.ndim==1) {
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol=0.0;
        double sx[par.nvar];
        cell *f[2];
        for (unsigned int ic=0; ic<b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v=0;v<par.nvar;v++) Ph[v] = buf[ibuf+v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v=0;v<par.ndim;v++) cpos[v] = buf[ibuf+v];
          ibuf += par.ndim;
          for (int v=0;v<par.nvar;v++) sx[v] = buf[ibuf+v];
          ibuf += par.nvar;
          list<cell*>::iterator f_iter=b->NGrecvC2F[ic].begin();
          for (int v=0;v<2;v++) {
            f[v] = *f_iter;
            f_iter++;
          }
          interpolate_coarse2fine1D(
                        par,grid,solver,Ph,c_vol,sx,f[0],f[1]);
        } // loop over coarse cells
      }   // if 1D
      else if (par.ndim==2) {
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol=0.0;
        int ipos[par.ndim]; ipos[0]=0; ipos[1]=0;
        double sx[par.nvar], sy[par.nvar];
        cell *f[4];
        for (unsigned int ic=0; ic<b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v=0;v<par.nvar;v++) Ph[v] = buf[ibuf+v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v=0;v<par.ndim;v++) cpos[v] = buf[ibuf+v];
          ibuf += par.ndim;
          for (int v=0;v<par.nvar;v++) sx[v] = buf[ibuf+v];
          ibuf += par.nvar;
          for (int v=0;v<par.nvar;v++) sy[v] = buf[ibuf+v];
          ibuf += par.nvar;
          list<cell*>::iterator f_iter=b->NGrecvC2F[ic].begin();
          for (int v=0;v<4;v++) {
            f[v] = *f_iter;
            f_iter++;
          }
          rep.printVec("cpos",cpos,par.ndim);
          rep.printVec("Ph",Ph,par.nvar);
          CI.get_ipos_vec(cpos,ipos);
          interpolate_coarse2fine2D(
                  par,grid,solver,Ph,ipos,c_vol,sx,sy,f[0],f[1],f[2],f[3]);
        } // loop over coarse cells
      }   // if 2D
      else {
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol=0.0;
        int ipos[par.ndim];
        double sx[par.nvar], sy[par.nvar], sz[par.nvar];
        cell *f[8];
        for (unsigned int ic=0; ic<b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v=0;v<par.nvar;v++) Ph[v] = buf[ibuf+v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v=0;v<par.ndim;v++) cpos[v] = buf[ibuf+v];
          ibuf += par.ndim;
          for (int v=0;v<par.nvar;v++) sx[v] = buf[ibuf+v];
          ibuf += par.nvar;
          for (int v=0;v<par.nvar;v++) sy[v] = buf[ibuf+v];
          ibuf += par.nvar;
          for (int v=0;v<par.nvar;v++) sz[v] = buf[ibuf+v];
          ibuf += par.nvar;
          list<cell*>::iterator f_iter=b->NGrecvC2F[ic].begin();
          for (int v=0;v<8;v++) {
            f[v] = *f_iter;
            f_iter++;
          }
          CI.get_ipos_vec(cpos,ipos);
          interpolate_coarse2fine2D(
                  par,grid,solver,Ph,ipos,c_vol,sx,sy,f[0],f[1],f[2],f[3]);
          rep.error("write 3D interpolation routine C2F",par.ndim);
        } // loop over coarse cells
      }   // if 3D  
    } // if 2nd order accurate 
      

    buf = mem.myfree(buf);
    buf=0;
  } // if parent proc is different to my proc.

#ifdef TEST_MPI_NG
  cout <<"NG_MPI_C2F_bc::BC_update_COARSE_TO_FINE_RECV() done\n";
  cout.flush();
#endif
  return 0;
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
  int bsize = grid->idx()*par.Nbc/2; // idx is >=2, Nbc is >=1.
  cout <<"bsize="<<bsize<<", and idx="<<grid->idx()<<", Nbc=";
  cout <<par.Nbc<<"\n";
  rep.printVec("ixmin",ixmin,par.ndim);
  rep.printVec("ixmax",ixmax,par.ndim);
  
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

  cout <<"boundary: x in ["<<xn<<","<<xp<<"], y in["<<yn<<","<<yp<<"]\n";
  size_t ct=0;
  cell *c = grid->FirstPt_All();
  do {
    rep.printVec("cpos",c->pos,par.ndim);
    if (c->pos[XX]>xn && c->pos[XX]<xp &&
        c->pos[YY]>yn && c->pos[YY]<yp) {
      bdata->c.push_back(c);
      ct++;
    }
  } while ((c=grid->NextPt_All(c)) !=0);
#ifdef TEST_MPI_NG
  cout <<"add_cells_to_C2F_send_list_2D: added "<<ct<<" cells.\n";
#endif

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





