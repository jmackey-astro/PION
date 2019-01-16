/// \file MCMD_boundaries.cpp
/// \brief Class definitions for periodic boundaries
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.


#include "boundaries/MCMD_boundaries.h"
using namespace std;

//#define TEST_COMMS

// ##################################################################
// ##################################################################



int MCMD_bc::BC_assign_BCMPI(
      class SimParams &par,     ///< pointer to simulation parameters
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,
      int comm_tag
      )
{
  //
  // first choose cells to send to the appropriate processor.  This
  // is stored in a separate list "send_data" that is part of the
  // boundary_data struct.
  //
  int err = 0;

  // set isdomain to true because these data are part of the domain.
  list<cell*>::iterator c=b->data.begin();
  do {
    c++;
  } while (c!=b->data.end());

#ifdef TEST_COMMS
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: sending data in dir: "<<b->dir<<"\n";
#endif
  int ncell =0;
  if (b->send_data.size() != 0) {
    rep.error("send_data is not empty!",b->send_data.size());
  }
  err += BC_select_data2send(par,grid,&(b->send_data), &ncell, b);
#ifdef TEST_COMMS
  cout <<"BC_assign_BCMPI: got "<<ncell<<" cells in send_data\n";
#endif
  return err;
}




// ##################################################################
// ##################################################################



int MCMD_bc::BC_select_data2send(
      class SimParams &par,     ///< pointer to simulation parameters
      class GridBaseClass *grid,  ///< pointer to grid.
      list<cell *> *l,
      int *nc,
      boundary_data *b
      )
{
  // Check inputs.
  if ( !(*l).empty() ) {
#ifdef TEST_COMMS
    rep.warning("BC_select_data2send: List not empty! Emptying it now.",0,(*l).size());
#endif
    (*l).clear(); if(!(*l).empty()) rep.error("emptying list.",(*l).empty());
  }
  if ( *nc !=0 ) {
#ifdef TEST_COMMS
    rep.warning("BC_select_data2send: uninitialized counter in input. setting it to zero.",0,*nc);
#endif
    *nc=0;
  }
  
  //
  // We want to get all cells on the grid that are adjacent to the
  // boundary data, so we just go BC_nbc cells in the "ondir"
  // direction from every boundary cell, and this is a "send" cell.
  //
  int count=0;
  list<cell*>::iterator c=b->data.begin();
  cell *temp =0;
  do {
    temp = *c;
    //CI.print_cell(temp);
    for (int v=0;v<par.Nbc;v++) temp = grid->NextPt(temp,b->ondir);
    (*l).push_back(temp);
    count++;
    ++c;
  } while (c!=b->data.end());

#ifdef TEST_COMMS
  cout <<"Got "<<count<<" cells, expected "<<b->data.size();
  cout <<"  list size = "<< (*l).size() <<"\n";
#endif
  *nc = count;

  return 0;
}




// ##################################################################
// ##################################################################



int MCMD_bc::BC_update_BCMPI(
      class SimParams &par,       ///< simulation parameters
      const int level,            ///< level in grid hierarchy
      class GridBaseClass *grid,  ///< pointer to grid.
      boundary_data *b,
      const int cstep,
      const int maxstep,
      int comm_tag
      )
{
#ifdef TEST_COMMS
  cout <<"*******************************************\n";
  cout <<"BC_update_BCMPI: sending data in dir: "<<b->dir;
  cout <<": sending "<<b->send_data.size();
  cout <<" cells.  Boundary data contains "<<b->data.size();
  cout <<" cells.\n";
  cout.flush();
#endif

  //
  // send the data.
  //
  int err = 0;
  class MCMDcontrol *ppar = &(par.levels[level].MCMD);
  string send_id;
#ifdef TEST_COMMS
  cout <<"BC_update_BCMPI: sending data...\n";
#endif 
  err += COMM->send_cell_data(
        ppar->ngbprocs[b->dir], // to_rank
        &(b->send_data),        // cells list.
        b->send_data.size(),    // number of cells.
        par.ndim, par.nvar,
        send_id,     // identifier for send.
        comm_tag
        );
  if (err) rep.error("sending data failed",err);

  //
  // receive data.  Expect to receive data from the same
  // direction as where we sent data.
  //
  //struct boundary_data *recv_b = 0;
  //int seek_tag=0;
  //recv_b = grid->BC_bd[b->dir];
  //if      (recv_b->itype == PERIODIC) seek_tag = BC_PERtag;
  //else if (recv_b->itype == BCMPI)    seek_tag = BC_MPItag;
  //else {
  //  // opposite direction is not an MPI boundary, so just wait for
  //  // return message from where we sent our data to.
  //  recv_dir = b->dir;
  //  seek_tag = comm_tag;
  //  recv_b = grid->BC_bd[b->dir];
  //}

  string recv_id; int recv_tag=-1; int from_rank=-1;
  err = COMM->look_for_data_to_receive(
        &from_rank, ///< rank of sender
        recv_id,    ///< identifier for receive.
        &recv_tag,  ///< comm_tag associated with data.
        comm_tag,   ///< requested comm_tag
        COMM_CELLDATA ///< type of data we want.
        );
  if (err) rep.error("look for cell data failed",err);
#ifdef TEST_COMMS
  cout <<"BC_update_BCMPI: got data to receive from rank: "<<from_rank<<"\n";
#endif 


#ifdef TEST_COMMS
  cout <<ppar->get_myrank()<<"\tBC_update_BCMPI: Receiving Data type ";
  cout <<recv_tag<<" from rank: "<<from_rank<<" from direction ";
  cout <<b->dir<<", expecting rank "<<ppar->ngbprocs[b->dir]<<"\n";
#endif 


  //
  // Receive the data:
  //
  err = COMM->receive_cell_data(
        from_rank,  ///< rank of process we are receiving from.
        &(b->data),     ///< list of cells to get data for. 
        b->data.size(), ///< number of cells in list (extra checking!)
        par.ndim, par.nvar,
        recv_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
        recv_id     ///< identifier for receive, for any book-keeping.
        );
  if (err) rep.error("COMM->receive_cell_data() returned error",err);

  //
  // If on a full timestep update, set P=Ph for received boundary.
  //
  if (cstep==maxstep) {
    //struct boundary_data *b2 = grid->BC_bd[static_cast<int>(dir)];
    list<cell*>::iterator c=b->data.begin();
    for (c=b->data.begin(); c!=b->data.end(); ++c) {
      for (int v=0;v<par.nvar;v++) (*c)->P[v] = (*c)->Ph[v];
#ifdef TEST_COMMS
      rep.printVec("P",(*c)->P,par.nvar);
#endif
    } // all cells.
  } // if full timestep.
  
  //
  // Wait until send is finished, and then return.
  //
  //  cout <<"BC_update_BCMPI: waiting for finish "<<b->dir<<"\n";
  err = COMM->wait_for_send_to_finish(send_id);
  if (err) rep.error("Waiting for send to complete!", err);
  //cout <<"BC_update_BCMPI: send and receive finished.\n";
  //cout <<"BC_update_BCMPI: Sent BC "<<b->dir<<", received BC "<<dir<<"\n";
  //cout <<"*******************************************\n\n";

  return 0;
}





// ##################################################################
// ##################################################################




