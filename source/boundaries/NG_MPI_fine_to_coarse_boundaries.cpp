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
    return;
  }

  // else we have to send data to another MPI process:
  class GridBaseClass *grid = par.level[l].grid;
  double dxo2 = 0.5*grid->DX();
  int nc = b->avg[0]->c.size();
  int nel = b->avg.size();

  // data to send will be ordered as position,conserved-var
  // for each averaged cell.
  pion_flt *data = mem.myalloc(data,nel*(par.nvar+par.ndim));

  // loop through avg vector and add cells and positions to
  // data array.
  int v=0, ix=0, iy=0, iz=0;
  double cd[par.nvar];
  size_t ct=0;
  for (v=0;v<nel;v++) {
    for (int j=0;j<par.nvar;j++) cd[j]=0.0;
    average_cells(par,solver,grid,nc,b->avg[v]->c,cd);
    for (int i=0;i<par.ndim;i++) data[ct+i] = b->avg[v]->cpos[i];
    ct += par.ndim;
    for (int i=0;i<par.nvar;i++) data[ct+i] = cd[i];
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
        pproc, ///< rank to send to.
        ct,      ///< size of buffer, in number of doubles.
        data,    ///< pointer to double array.
        id.str(), ///< identifier for send, for tracking delivery.
        BC_MPI_NG_tag ///< comm_tag, to say what kind of send this is
        );
  if (err) rep.error("Send_F2C send_data failed.",err);
#ifdef TEST_MPI_NG
  cout <<"BC_update_FINE_TO_COARSE_SEND: returned with id="<<id[i];
  cout <<"\n";
#endif

  // store ID to clear the send later (and delete the MPI temp data)
  NG_send_list.push_back(id);
  data = mem.myfree(data);
  
  return 0;
}





// ##################################################################
// ##################################################################



void NG_fine_to_coarse_MPI_bc::BC_FINE_TO_COARSE_SEND_clear_sends()
{
  for (int i=0;i<NG_send_list.size();i++) {
    COMM->wait_for_send_to_finish(NG_send_list[i];
  }
  NG_send_list.clear();
  return;
}





// ##################################################################
// ##################################################################



int NG_fine_to_coarse_MPI_bc::BC_assign_FINE_TO_COARSE_RECV(
      class SimParams &par,     ///< pointer to simulation parameters
      const int l,  ///< level of this grid.
      boundary_data *b  ///< boundary data
      )
{
  // Check if child grids exist or are on my MPI process
  int err=0;
  class MCMDcontrol *MCMD = &(par.level[l].MCMD);
  int nchild = MCMD->child_procs.size();
  b->NGrecv.resize(nchild);

  // loop over children:
  for (int i=0;i<nchild;i++) {

    if (MCMD->get_myrank() == child_procs[i].rank) {
      // If child is on my grid call serial version that just grabs data
      // directly from the child grid.
#ifdef TEST_MPI_NG
      cout <<"my rank == child rank, calling serial ";
      cout <<"FINE_TO_COARSE\n";
#endif
      err = NG_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE(
          par, par.level[l].grid, b, par.level[l].child);
      rep.errorTest("BC_assign_FINE_TO_COARSE_RECV serial",0,err);
    }
    else {
      //
      // else we have to create a list of cells for each off-grid
      // child, and add them to b->NGrec[i]
      //
#ifdef TEST_MPI_NG
      cout <<"my rank != child rank, running parallel ";
      cout <<"FINE_TO_COARSE_RECV\n";
#endif
      // get dimensions of child grid from struct
      

  } // if child is on a different MPI process

  return 0;
}


  return 0;
}




// ##################################################################
// ##################################################################



