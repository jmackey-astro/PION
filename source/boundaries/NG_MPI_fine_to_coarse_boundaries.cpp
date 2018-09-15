/// \file NG_MPI_fine_to_coarse_boundaries.cpp
/// \brief Class definitions for NG_MPI_fine_to_coarse boundaries
/// \author Jonathan Mackey
/// 
/// \section Description
/// Coarse-to-Fine boundary updates send grid data from a coarse grid
/// to update the external boundaries of a fine grid contained within
/// the coarse grid.  This version is parallelised with MPI.
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/NG_MPI_fine_to_coarse_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;


// ##################################################################
// ##################################################################



int NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_SEND(
      class SimParams &par,     ///< pointer to simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b  ///< boundary data
      )
{
  // Check if parent grid is on my MPI process
  int pproc = par.levels[l].MCMD.parent_proc;
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);

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
    class GridBaseClass *grid = par.levels[l].grid;
    double dxo2 = 0.5*grid->DX();
    int ipos[MAX_DIM];
    int nc=1;  // number of fine cells per coarse cell
    for (int i=0;i<par.ndim;i++) nc*=2;
    int nel = grid->Ncell()/nc; // number of coarse cells
    // initialise the "avg" vector to hold the average state.
    b->avg.resize(nel);
    for (int v=0;v<nel;v++) {
      b->avg[v].c.resize(nc);
      b->avg[v].avg_state = mem.myalloc(b->avg[v].avg_state,
                                                  par.nvar);
    }
    // loop through avg vector and add cells and positions
    cell *f = grid->FirstPt();
    int v=0, ix=0, iy=0, iz=0;
    for (v=0;v<nel;v++) {
      // add fine cells to each avg struct
      b->avg[v].c[0] = f;
      b->avg[v].c[1] = grid->NextPt(f,XP);
      if (par.ndim>1) {
        b->avg[v].c[2] = grid->NextPt(b->avg[v].c[0],YP);
        b->avg[v].c[3] = grid->NextPt(b->avg[v].c[1],YP);
      }
      if (par.ndim>2) {
        b->avg[v].c[4] = grid->NextPt(b->avg[v].c[0],ZP);
        b->avg[v].c[5] = grid->NextPt(b->avg[v].c[1],ZP);
        b->avg[v].c[6] = grid->NextPt(b->avg[v].c[2],ZP);
        b->avg[v].c[7] = grid->NextPt(b->avg[v].c[3],ZP);
      }
      // add position of coarse cell centre (for testing)
      for (int i=0;i<par.ndim;i++)
        ipos[v] = b->avg[v].c[0]->pos[v]+dxo2;
      for (int i=par.ndim;i<MAX_DIM;i++)
        ipos[v] = 0.0;
      CI.get_dpos_vec(ipos,b->avg[v].cpos);

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
        if (par.ndim>1) {
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
            if (par.ndim>2) {
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



int NG_MPI_fine_to_coarse_bc::BC_update_FINE_TO_COARSE_SEND(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int l, ///< level in the NG grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{

  // Check if parent grid is on my MPI process
  int pproc = par.levels[l].MCMD.parent_proc;
  int err=0;
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);

  // If so, do nothing because the parent grid will call the 
  // serial NG setup function, which just grabs the data from 
  // this grid.
  if (MCMD->get_myrank() == pproc) {
#ifdef TEST_MPI_NG
    cout <<"my rank == parent rank, no need to send anything ";
    cout <<"FINE_TO_COARSE_SEND\n";
#endif
    return 0;
  }

  // else we have to send data to another MPI process:
  class GridBaseClass *grid = par.levels[l].grid;
  int nc = b->avg[0].c.size();
  int nel = b->avg.size();

  // data to send will be ordered as position,conserved-var
  // for each averaged cell.
  double *data = 0;
  data = mem.myalloc(data,nel*(par.nvar+par.ndim));

  // loop through avg vector and add cells and positions to
  // data array.
  int v=0;
  double cd[par.nvar];
  size_t ct=0;
  for (v=0;v<nel;v++) {
    for (int j=0;j<par.nvar;j++) cd[j]=0.0;
    average_cells(par,solver,grid,nc,b->avg[v].c,cd);
    for (int i=0;i<par.ndim;i++) data[ct+i] = b->avg[v].cpos[i];
    ct += par.ndim;
    for (int i=0;i<par.nvar;i++) data[ct+i] = cd[i];
    ct += par.nvar;
  } // go to next avg element.

  //
  // Send data using a non-blocking MPI send
  //
  string id;
  //id <<"F2C_"<<MCMD->get_myrank()<<"_to_"<<pproc;
#ifdef TEST_MPI_NG
  cout <<"BC_update_FINE_TO_COARSE_SEND: Sending "<<ct;
  cout <<" doubles from proc "<<MCMD->get_myrank();
  cout <<" to parent proc "<<pproc<<"\n";
#endif
  err += COMM->send_double_data(
        pproc,ct,data,id,BC_MPI_NGF2C_tag);
  if (err) rep.error("Send_F2C send_data failed.",err);
#ifdef TEST_MPI_NG
  cout <<"BC_update_FINE_TO_COARSE_SEND: returned with id="<<id[i];
  cout <<"\n";
#endif

  // store ID to clear the send later (and delete the MPI temp data)
  NG_F2C_send_list.push_back(id);
  data = mem.myfree(data);
  
  return 0;
}





// ##################################################################
// ##################################################################



void NG_MPI_fine_to_coarse_bc::BC_FINE_TO_COARSE_SEND_clear_sends()
{
  for (unsigned int i=0;i<NG_F2C_send_list.size();i++) {
    COMM->wait_for_send_to_finish(NG_F2C_send_list[i]);
  }
  NG_F2C_send_list.clear();
  return;
}





// ##################################################################
// ##################################################################



int NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_RECV(
      class SimParams &par,     ///< pointer to simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b  ///< boundary data
      )
{
  // Check if child grids exist or are on my MPI process
  int err=0;
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int nchild = MCMD->child_procs.size();
  b->NGrecvF2C.resize(nchild);

  // loop over children:
  for (int i=0;i<nchild;i++) {

    if (MCMD->get_myrank() == MCMD->child_procs[i].rank) {
      // If child is on my grid call serial version that just grabs data
      // directly from the child grid.
#ifdef TEST_MPI_NG
      cout <<"my rank == child rank, calling serial ";
      cout <<"FINE_TO_COARSE\n";
#endif
      err = NG_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE(
          par, par.levels[l].grid, b, par.levels[l].child);
      rep.errorTest("BC_assign_FINE_TO_COARSE_RECV serial",0,err);
    }
    else {
      //
      // else we have to create a list of cells for each
      // child, and add them to b->NGrec[i]
      //
#ifdef TEST_MPI_NG
      cout <<"my rank != child rank, running parallel ";
      cout <<"FINE_TO_COARSE_RECV\n";
#endif
      // get dimensions of child grid from struct
      int ixmin[MAX_DIM], ixmax[MAX_DIM];
      CI.get_ipos_vec(MCMD->child_procs[i].Xmin, ixmin);
      CI.get_ipos_vec(MCMD->child_procs[i].Xmax, ixmax);

      class GridBaseClass *grid = par.levels[l].grid;
      cell *c = grid->FirstPt();
      size_t ct=0;
      do {
        bool ongrid=true;
        for (int j=0;j<par.ndim;j++) {
          if (c->pos[j]<ixmin[j] || c->pos[j]>ixmax[j]) ongrid=false;
        }
        if (ongrid) {
          b->NGrecvF2C[i].push_back(c);
          ct++;
        }
      } while ((c=grid->NextPt(c)) !=0);

#ifdef TEST_MPI_NG
    cout <<"i="<<i<<" added "<<ct<<" cells to boundary, ";
    cout <<"Ncell/2^ndim="<<grid->Ncell()/pow(2,par.ndim)<<"\n";
#endif
    } // if child is on a different MPI process

  } // loop over child grids

  return 0;
}



// ##################################################################
// ##################################################################



int NG_MPI_fine_to_coarse_bc::BC_update_FINE_TO_COARSE_RECV(
      class SimParams &par,      ///< pointer to simulation parameters
      class FV_solver_base *solver, ///< pointer to equations
      const int l, ///< level in the NG grid structure
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  // Check if child grids exist or are on my MPI process
  int err=0;
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int nchild =  b->NGrecvF2C.size();
  int ichild = 0;

  // loop over children twice, once for child grids that are on my
  // MPI process, and once for grids that on other processes
  for (int i=0;i<nchild;i++) {

    if (MCMD->get_myrank() == MCMD->child_procs[i].rank) {
      // If child is on my grid call serial version that just grabs data
      // directly from the child grid.
#ifdef TEST_MPI_NG
      cout <<"my rank == child rank, calling serial ";
      cout <<"FINE_TO_COARSE update\n";
#endif
      err = NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
          par, solver, l, b, cstep, maxstep);
      rep.errorTest("BC_update_FINE_TO_COARSE_RECV serial",0,err);
      ichild++;
    }
  }

  while (ichild < nchild) {

    //
    // else we go through the list of cells for each on-grid
    // child (b->NGrec[i]), get the data from another MPI process,
    // and then update the cells
    //
#ifdef TEST_MPI_NG
    cout <<"my rank != child rank, running parallel ";
    cout <<"FINE_TO_COARSE_RECV\n";
#endif
    
    //
    // receive data.
    //
    string recv_id; int recv_tag=-1; int from_rank=-1;
    err = COMM->look_for_data_to_receive(
          &from_rank, // rank of sender (output)
          recv_id,    // identifier for receive (output).
          &recv_tag,  // comm_tag associated with data (output)
          COMM_DOUBLEDATA // type of data we want.
          );
    if (err) rep.error("look for double data failed",err);
#ifdef TEST_MPI_NG
    cout <<"BC_update_FINE_TO_COARSE_RECV: found data from rank ";
    cout <<from_rank<<"\n";
#endif 
    // associate data with one of the child grids:
    int i=-1;
    for (int v=0;v<nchild;v++) {
      if (MCMD->get_myrank() == MCMD->child_procs[v].rank) i = v;
    }
    if (i<0) rep.error("receiving data, but don't know why",i);

    // receive the data
    size_t nel = b->NGrecvF2C[i].size();
    size_t ct = nel*(par.nvar+par.ndim);
    pion_flt *buf = mem.myalloc(buf,ct);
#ifdef TEST_MPI_NG
    cout <<"BC_update_FINE_TO_COARSE_RECV: get "<<nel<<" cells.\n";
#endif
    //
    // Receive data into buffer.
    //
    err = COMM->receive_double_data(from_rank, recv_tag, recv_id,
                                    ct, buf);
    if (err) rep.error("(BC_update_F2C_RECV) getdata failed",err);

    // Go through list of cells in child i, and put the received
    // data onto these cells.
    list<cell*>::iterator c_iter=b->NGrecvF2C[i].begin();
    cell *c=0;
    pion_flt pos[MAX_DIM], prim[par.nvar];
    size_t i_el=0;
    for (c_iter=b->NGrecvF2C[i].begin();
         c_iter!=b->NGrecvF2C[i].end(); ++c_iter) {
      c = (*c_iter);
      CI.get_dpos(c,pos);
#ifdef TEST_MPI_NG
      rep.printVec("cell pos", pos, par.ndim);
      rep.printVec("recv pos", &(buf[i_el]), par.ndim);
#endif
      i_el += par.ndim;
      solver->UtoP(&(buf[i_el]),prim,
                   par.EP.MinTemperature,par.gamma);
#ifdef TEST_MPI_NG
      rep.printVec("cell Prim", c->Ph, par.nvar);
      rep.printVec("recv Prim", prim, par.nvar);
#endif
      for (int v=0;v<par.ndim;v++) c->Ph[v] = prim[v];
      if (cstep==maxstep) {
        for (int v=0;v<par.ndim;v++) c->P[v] = c->Ph[v];
      }
      i_el += par.nvar;
    } // loop over cells

#ifdef TEST_MPI_NG
    cout <<"(BC_update_F2C_RECV) i_el="<<i_el<<" of "<<ct;
    cout <<" total elements.\n";
#endif
    buf = mem.myfree(buf);
    ichild++;
  } // loop until we got through all child grids.

  return 0;
}




// ##################################################################
// ##################################################################




