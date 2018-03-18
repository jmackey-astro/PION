/// \file uniform_grid_pllel.cpp
/// \author Jonathan Mackey
/// 
/// Definitions for parallel uniform grid.
/// 
/// Modifications:\n
/// - 2007-11-06 started to get it going.
/// - 2010-01-22 JM: Changed RT source finding if we have non-cell-centred sources.
/// - 2010-01-26 JM: Debugging 3D RT for vertex-centred sources.
/// - 2010-02-03 JM: changed variable names in destructor ('i' was defined twice!)
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; Added oneway-outflow BC
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.  Switched 'endl' statements for '\n'.
/// - 2010.11.19 JM: Got rid of testing myalloc() myfree() functions.
/// - 2010.12.07 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  Will probably keep it since it is the way
///   things will go eventually.  The new grid classes have extra
///   function for the distance between two points, two cells, and
///   between a vertex and a cell.
/// - 2011.02.24 JM: Worked on generalising multi-source radiative transfer.
///    Still some way to go to getting it functional.  Wraps old setup() 
///    function for point sources into a new function.  The setup function now
///    calls either this new function, or another new one for sources at infinity.
/// - 2011.02.25 JM: removed HCORR ifdef around new code. 
///     Added possibility of multiple sources in send/recv.
/// - 2011.02.28 JM: Got rid of RSP radiation-sources-parameters class.
/// - 2011.03.02 JM: Fixed bugs in diffuse radiation boundary comms.
/// - 2011.03.21 JM: Got rid of zero-ing of column densities.  It is done in
///    the cell constructor, so we don't need to worry about it.
/// - 2011.03.22 JM: small bugfixes.
/// - 2011.04.22 JM: Fixed setup_recv_boundaries() functions so that new boundary cells
///    are only created if they don't already exist.  It now works for multiple point
///    sources, and I don't think they need to be at the same place.
/// - 2013.09.05 JM: Debugged for new get/set col functions.
/// - 2015.01.28 JM: Removed GEOMETRIC_GRID (it is default now), and
///    updated for the new code structure.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
/// - 2016.03.21 JM: Worked on simplifying RT boundaries (first hack!
///    compiles but buggy...) 03.24:fixed some bugs, redefined dir_XP
/// - 2017.11.07-22 JM: updating boundary setup.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"


#include "grid/uniform_grid.h"
#include "sim_control_MPI.h"
#include "microphysics/microphysics_base.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#ifdef PARALLEL



// ##################################################################
// ##################################################################



UniformGridParallel::UniformGridParallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn, ///< array of minimum values of x,y,z for this grid.
      double *xp, ///< array of maximum values of x,y,z for this grid.
      int *nc,    ///< array of number of cells in x,y,z directions.
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp, ///< array of max. x/y/z for full simulation.
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{
#ifdef TESTING
  cout <<"UniformGridParallel constructor.\n";
  rep.printVec("Local Xmin",xn,nd);
  rep.printVec("Local Xmax",xp,nd);
  rep.printVec("Local Npt ",nc,nd);
#endif

#ifdef PLLEL_RT
  //RT_nbd=0;
  //RT_nbc=0;
  //  RT_bd=0;
#endif // PLLEL_RT

#ifdef TESTING
  rep.printVec("SIM Xmin ", Sim_xmin, G_ndim);
  rep.printVec("SIM Xmax ", Sim_xmax, G_ndim);
  rep.printVec("SIM Range", Sim_range,G_ndim);
  rep.printVec("SIM iXmin ", Sim_ixmin, G_ndim);
  rep.printVec("SIM iXmax ", Sim_ixmax, G_ndim);
  rep.printVec("SIM iRange", Sim_irange,G_ndim);
#endif
  
  //
  // Set pointer to muli-core class.
  //
  mpiPM = p;

#ifdef TESTING
  cout <<"UniformGridParallel constructor done.\n";
#endif
  return;
}



// ##################################################################
// ##################################################################



UniformGridParallel::~UniformGridParallel()
{
#ifdef TESTING
  cout <<"UniformGridParallel destructor.\n";
#endif
#ifdef PLLEL_RT
  //
  // Delete cells created for RT receive boundaries.
  // Have to delete lists of cells, and boundary_data, and maybe cells.
  //
  for (unsigned int isrc=0;isrc<RT_source_list.size(); isrc++) {
    struct RT_source_comms_info *rc = &(RT_source_list[isrc]);
    if (!rc->RT_recv_list.empty()) {
      struct boundary_data *b;
      for (unsigned int ii=0; ii<rc->RT_recv_list.size(); ii++) {
        int n=rc->RT_recv_list[ii].dir;
        b = rc->RT_recv_list[ii].RT_bd;
        //
        // If a pre-existing boundary, just delete the list.
        //
        if (n==dir_XN || n==dir_XP || n==dir_YN || n==dir_YP || n==dir_ZN || n==dir_ZP) {
          //cout <<"\tNot deleting standard boundary data for dir: "<<n<<"\n";
          if      (!b) {
            //cout <<"\t    Zero boundary data pointer for dir: "<<n<<"\n";
          }
          else if (!b->data.empty()) {
            //cout <<"\t    Emptying list.\n";
            b->data.clear(); // No dynamically allocated memory in these boundaries.
          }
          else {
            //cout <<"\t    Empty boundary list (!)(?) Probably shouldn't happen.\n";
          }
        } // if pre-existing.
        

        //cout <<"\t\tNow deleting boundary data...\n";
        if (b) {
          b = mem.myfree(b);
          //RT_recv_list[ii].RT_bd = mem.myfree(RT_recv_list[ii].RT_bd);
        }
        //cout <<"\t\tFinished with this Receive boundary\n";
      }
      rc->RT_recv_list.clear();
    }
    
    //
    // Delete send lists.  No cells to delete, just the lists, and boundary_data
    //
    if (!rc->RT_send_list.empty()) {
      //cout <<"\tDeleting RT Send List.\n";
      struct boundary_data *b;
      for (unsigned int ii=0; ii<rc->RT_send_list.size(); ii++) {
        b = rc->RT_send_list[ii].RT_bd;
        if      (!b)  {
          //cout <<"\t\tZero send boundary data pointer for dir: "<<RT_send_list[i].dir<<"\n";
        }
        else if (!b->data.empty()) {
          //cout <<"\t\tEmptying send list.\n";
          b->data.clear(); // No dynamically allocated memory in these boundaries.
        }
        else {
          //cout <<"\t\tEmpty send boundary list (!)(?) Probably shouldn't happen.\n";
        }
        b = mem.myfree(b);
      }
      rc->RT_send_list.clear();
      //cout <<"\tDeleting RT Send List... done!\n";
    }
    
  }
#endif // PLLEL_RT

  return;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_setBCtypes(
  class SimParams &par
  )
{
  string fname="UniformGridParallel::BC_setBCtypes";
#ifdef TESTING
  cout <<"PLLEL: Set BC types...\n";
#endif 
  //
  // call serial version of setBCtypes, to set up the boundaries
  //
  int err = UniformGrid::BC_setBCtypes(par);
  if (err) {
    rep.error("UniGridPllel::BC_setBCtypes:: serial call error",err);
  }
  
  //
  // Now go through the 6 directions, and see if we need to replace
  // any edge boundaries with MPI communication boundaries, if the
  // local grid does not reach the global grid boundary.
  //
  int i=0;
  enum axes dir;
  string temp;
  // Loop through boundaries, and if local boundary is not sim boundary,
  // set it to be a parallel boundary.
  for (i=0; i<G_ndim; i++) {
    dir = static_cast<axes>(i);
    if (!pconst.equalD(G_xmin[i], Sim_xmin[i])) {
      // local xmin is not Sim xmin, so it's an mpi boundary
      BC_bd[2*i].itype=BCMPI;
      BC_bd[2*i].type="BCMPI";
    }
    if (!pconst.equalD(G_xmax[i], Sim_xmax[i])) {
      // local xmax is not Sim xmin, so it's an mpi boundary
      BC_bd[2*i+1].itype=BCMPI;
      BC_bd[2*i+1].type="BCMPI";
    }
  }
  
#ifdef TESTING
  cout <<"PLLEL: BC types and data set up.\n";
#endif 

  //
  // If we have periodic boundaries, need to set neighbouring processors to
  // wrap around.  So set the number of procs in each direction.
  //
  int nx[par.ndim];
  for (i=0;i<G_ndim;i++) {
    nx[i] =static_cast<int>(ONE_PLUS_EPS*Sim_range[i]/G_range[i]);
  }
  for (i=0; i<2*G_ndim; i++) {
    if (BC_bd[i].itype == PERIODIC) {
      switch (i) {
       case XN:
        mpiPM->ngbprocs[XN] = mpiPM->get_myrank() +nx[XX] -1;
        break;
       case XP:
        mpiPM->ngbprocs[XP] = mpiPM->get_myrank() -nx[XX] +1;
        break;
       case YN:
        mpiPM->ngbprocs[YN] = mpiPM->get_myrank() +(nx[YY]-1)*nx[XX];
        break;
       case YP:
        mpiPM->ngbprocs[YP] = mpiPM->get_myrank() -(nx[YY]-1)*nx[XX];
        break;
       case ZN:
        mpiPM->ngbprocs[ZN] = mpiPM->get_myrank() +
                              (nx[ZZ]-1)*nx[YY]*nx[XX];
        break;
       case ZP:
        mpiPM->ngbprocs[ZP] = mpiPM->get_myrank() -
                              (nx[ZZ]-1)*nx[YY]*nx[XX];
        break;
       default:
        rep.error("UniformGridParallel::BC_setBCtypes: Bad direction",i);
        break;
      } // set neighbour according to direction.

      if ( (mpiPM->ngbprocs[i]<0) ||
           (mpiPM->ngbprocs[i]>=mpiPM->get_nproc()) )
        rep.error("UniformGridParallel::BC_setBCtypes: Bad periodic \
                   neighbour",mpiPM->ngbprocs[i]);
      if (mpiPM->ngbprocs[i] == mpiPM->get_myrank()) {
        //  cout <<"UniformGridParallel::BC_setBCtypes: only one proc in dir [i]: "<<i<<"\n";
        //  cout <<"UniformGridParallel::BC_setBCtypes: periodic on single proc, so setting ngb to -999.\n";
        mpiPM->ngbprocs[i] = -999;
      }
    } // if periodic  
#ifdef TESTING
    cout<<"Neighbouring processor in dir "<<i<<" = "<<mpiPM->ngbprocs[i]<<"\n";
#endif // TESTING
  } // loop over directions.
  return(0);
}



// ##################################################################
// ##################################################################



int UniformGridParallel::assign_boundary_data(
        const double simtime,     ///< current simulation time
        const double sim_start,   ///< start time of simulation
        const double sim_finish,  ///< finish time of simulation
        const double Tmin         ///< minimum temperature allowed
        )
{
#ifdef TESTING
  cout<<"UniformGridParallel::assign_boundary_data()\n";
#endif // TESTING
  //
  // Loop through all boundaries, and assign data to them.
  // This is only different to the serial version in that it includes
  // BCMPI boundaries (MPI communication of internal boundaries).
  //
  int err=0;
  for (int i=0; i<BC_nbd; i++) {
    switch (BC_bd[i].itype) {
     case PERIODIC:   err += BC_assign_PERIODIC(  &BC_bd[i]); break;
     case OUTFLOW:    err += BC_assign_OUTFLOW(   &BC_bd[i]); break;
     case ONEWAY_OUT: err += BC_assign_ONEWAY_OUT(&BC_bd[i]); break;
     case INFLOW:     err += BC_assign_INFLOW(    &BC_bd[i]); break;
     case REFLECTING: err += BC_assign_REFLECTING(&BC_bd[i]); break;
     case FIXED:      err += BC_assign_FIXED(     &BC_bd[i]); break;
     case JETBC:      err += BC_assign_JETBC(     &BC_bd[i]); break;
     case JETREFLECT: err += BC_assign_JETREFLECT(&BC_bd[i]); break;
     case DMACH: err += BC_assign_DMACH(simtime,  &BC_bd[i]); break;
     case DMACH2:     err += BC_assign_DMACH2(    &BC_bd[i]); break;
     case BCMPI: err += BC_assign_BCMPI(&BC_bd[i],BC_MPItag); break;
     case STWIND:     err += BC_assign_STWIND(simtime, sim_start,
                                sim_finish, Tmin, &BC_bd[i]); break;
     case STARBENCH1: err += BC_assign_STARBENCH1(&BC_bd[i]); break;
     default:
      rep.warning("Unhandled BC",BC_bd[i].itype,-1); err+=1; break;
    }
    if (i==XP || (i==YP && G_ndim>1) || (i==ZP && G_ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      COMM->barrier("UniformGridParallel__SetupBCs");
    }
  }
#ifdef TESTING
  cout<<"UniformGridParallel::assign_boundary_data() finished\n";
#endif // TESTING
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::TimeUpdateExternalBCs(
      const double simtime,   ///< current simulation time
      const int cstep,
      const int maxstep
      )
{
  struct boundary_data *b;
  int i=0; int err=0;
  for (i=0;i<BC_nbd;i++) {
    b = &BC_bd[i];
    switch (b->itype) {
    case PERIODIC:   err += BC_update_PERIODIC(   b, cstep, maxstep); break;
    case OUTFLOW:    err += BC_update_OUTFLOW(    b, cstep, maxstep); break;
    case ONEWAY_OUT: err += BC_update_ONEWAY_OUT( b, cstep, maxstep); break;
    case INFLOW:     err += BC_update_INFLOW(     b, cstep, maxstep); break;
    case REFLECTING: err += BC_update_REFLECTING( b, cstep, maxstep); break;
    case FIXED:      err += BC_update_FIXED(      b, cstep, maxstep); break;
    case JETBC:      err += BC_update_JETBC(      b, cstep, maxstep); break;
    case JETREFLECT: err += BC_update_JETREFLECT( b, cstep, maxstep); break;
    case DMACH:      err += BC_update_DMACH(simtime, b, cstep, maxstep); break;
    case DMACH2:     err += BC_update_DMACH2(     b, cstep, maxstep); break;
    case BCMPI: err += BC_update_BCMPI( b, cstep, maxstep,BC_MPItag); break;
    case STARBENCH1: err += BC_update_STARBENCH1( b, cstep, maxstep); break;
    case STWIND: case RADSHOCK: case RADSH2:
      //
      // These are updated in TimeUpdateInternalBCs
      //
      break;
    default:
       rep.warning("Unhandled BC: parallel update external",b->itype,-1); err+=1; break;
    }
    if (i==XP || (i==YP && G_ndim>=2) || (i==ZP && G_ndim==3)) {
      // Need to make sure all processors are updating boundaries along
      // each axis together, so that there is no deadlock from one 
      // process sending data in y/z-dir into a processor that is still
      // working on x-dir.
      //      cout <<"Barrier synching for dir="<<i<<" with 1=XP,2=YP,3=ZP.\n";
      ostringstream tmp; tmp <<"UniformGridParallel_TimeUpdateExternalBCs_"<<i;
      COMM->barrier(tmp.str());
      tmp.str("");
    }
  }
  return 0;
} // time update BCs.



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_assign_PERIODIC(
      boundary_data *b
      )
{
  //
  // For parallel grid, periodic data may be on a different proc.,
  // which is already pointed to by mpiPM->ngbprocs[b->dir]
  // So I just have to call BC_assign_BCMPI and it will do the job.
  int err=0;
  if (mpiPM->ngbprocs[b->dir] <0) {
    // cout <<"BC_assign_PERIODIC: non comm periodic in direction "<<b->dir<<"\n";
    err = UniformGrid::BC_assign_PERIODIC(b);
  }
  else {
    // cout<<"BC_assign_PERIODIC: communicating periodic bc in direction "<<b->dir<<"\n";
    // cout<<"BC_assign_PERIODIC: calling mpi assign BC function\n";
    err = UniformGridParallel::BC_assign_BCMPI(b,BC_PERtag);
  }
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_update_PERIODIC(
      struct boundary_data *b,
      const int cstep,
      const int maxstep
      )
{
  //
  // For parallel grid, periodic data can be on a different proc.,
  // which is already pointed to by mpiPM->ngbprocs[b->dir]
  // So I just have to call BC_update_BCMPI and it will do the job.
  //
  int err=0;
  if (mpiPM->ngbprocs[b->dir] <0) {
#ifdef TESTING
    cout <<"BC_update_PERIODIC: non-communicating periodic BC in ";
    cout <<"direction "<<b->dir<<"\n";
#endif
    err = UniformGrid::BC_update_PERIODIC(b,cstep,maxstep);
  }
  else {
#ifdef TESTING
    cout<<"BC_update_PERIODIC: communicating periodic BC in ";
    cout<<"direction "<<b->dir<<"\n";
    cout<<"BC_update_PERIODIC: calling mpi update BC function\n";
#endif
    err = UniformGridParallel::BC_update_BCMPI(b,cstep,maxstep,BC_PERtag);
  }
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_assign_BCMPI(
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
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: sending data in dir: "<<b->dir<<"\n";
#endif
  int ncell =0;
  if (b->send_data.size() != 0) {
    rep.error("send_data is not empty!",b->send_data.size());
  }
  err += BC_select_data2send(&(b->send_data), &ncell, b);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: got "<<ncell<<" cells in send_data\n";
#endif

  //
  // This is the same as the update function, except that we want
  // to set P[] and Ph[] vectors to the same values, so we set cstep
  // equal to maxstep.
  //
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_assign_BCMPI: starting\n";
#endif 
  err = BC_update_BCMPI(b,2,2,comm_tag);
#ifdef TESTING
  cout <<"BC_assign_BCMPI: finished\n";
  cout <<"*******************************************\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::BC_update_BCMPI(
      boundary_data *b,
      const int cstep,
      const int maxstep,
      int comm_tag
      )
{
#ifdef TESTING
  cout <<"*******************************************\n";
  cout <<"BC_update_BCMPI: sending data in dir: "<<b->dir<<"\n";
  cout <<"BC_update_BCMPI: sending "<<b->send_data.size();
  cout <<" cells.  Boundary data contains "<<b->data.size();
  cout <<" cells.\n";
  cout.flush();
#endif // TESTING

  //
  // send the data.
  //
  int err = 0;
  string send_id;
#ifdef TESTING
  cout <<"BC_update_BCMPI: sending data...\n";
#endif 
  err += COMM->send_cell_data(
        mpiPM->ngbprocs[b->dir], // to_rank
        &(b->send_data),        // cells list.
        b->send_data.size(),    // number of cells.
        G_ndim, G_nvar,
        send_id,     // identifier for send.
        comm_tag
        );
  if (err) rep.error("sending data failed",err);

  //
  // receive data.
  //
  string recv_id; int recv_tag=-1; int from_rank=-1;
  err = COMM->look_for_data_to_receive(
        &from_rank, ///< rank of sender
        recv_id,    ///< identifier for receive.
        &recv_tag,  ///< comm_tag associated with data.
        COMM_CELLDATA ///< type of data we want.
        );
  if (err) rep.error("look for cell data failed",err);
#ifdef TESTING
  cout <<"BC_update_BCMPI: got data to receive from rank: "<<from_rank<<"\n";
#endif 

  //
  // See what direction we are getting data from:
  // This is complicated for the case where we get an internal and a
  // periodic boundary from the same process -- we need to decide
  // which one we are getting.
  //
#ifdef TESTING
  cout <<"BC_update_BCMPI: associating direction with rank.\n";
#endif 
  enum direction dir=NO;
  for (int i=0; i<2*G_ndim; i++) {
    if (from_rank == mpiPM->ngbprocs[i]) {
      if (dir==NO) dir = static_cast<direction>(i);
      else {
        //
        // must be receiving periodic and internal boundary from same proc.
        // In this case we need to decide which it is, and whether the source
        // direction is the positive or negative direction.
        // N.B. if none of the conditions below are satisfied, then the earlier 
        // assignment of dir was the correct one.  This is just for if the first
        // assignment was the wrong one.
        //
        if (recv_tag == BC_PERtag) {
          // Periodic boundary, so in neg.dir, from_rank should be > myrank, and vice versa.
          if      ( from_rank>mpiPM->get_myrank() && 
                    (i==XN || i==YN || i==ZN)) {
            dir = static_cast<direction>(i);
          }
          else if ( from_rank<mpiPM->get_myrank() &&
                    (i==XP || i==YP || i==ZP)) {
            dir = static_cast<direction>(i);
          }
        }
        else if (recv_tag == BC_MPItag) {
          // Internal boundary, so in neg.dir. from_rank should be < myrank.
          if      ( from_rank<mpiPM->get_myrank() &&
                    (i==XN || i==YN || i==ZN)) {
            dir = static_cast<direction>(i);
          }
          else if ( from_rank>mpiPM->get_myrank() &&
                    (i==XP || i==YP || i==ZP)) {
            dir = static_cast<direction>(i);
          }
        }
        else rep.error("bad tag",recv_tag);
      }
    }
  }
  if (dir==NO) rep.error("Message is not from a neighbour!",from_rank);
#ifdef TESTING
  cout <<mpiPM->get_myrank()<<"\tBC_update_BCMPI: Receiving Data type ";
  cout <<recv_tag<<" from rank: "<<from_rank<<" from direction "<<dir<<"\n";
#endif 

  //
  // Now choose the boundary data associated with boundary we are
  // receiving:
  //
  struct boundary_data *recv_b = &BC_bd[static_cast<int>(dir)];

  //
  // Receive the data:
  //
  err = COMM->receive_cell_data(
        from_rank,  ///< rank of process we are receiving from.
        &(recv_b->data),     ///< list of cells to get data for. 
        recv_b->data.size(), ///< number of cells in list (extra checking!)
        G_ndim, G_nvar,
        recv_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
        recv_id     ///< identifier for receive, for any book-keeping.
        );
  if (err) rep.error("COMM->receive_cell_data() returned error",err);

  //
  // If on a full timestep update, set P=Ph for received boundary.
  //
  if (cstep==maxstep) {
    struct boundary_data *b2 = &BC_bd[static_cast<int>(dir)];
    list<cell*>::iterator c=b2->data.begin();
    for (c=b2->data.begin(); c!=b2->data.end(); ++c) {
      for (int v=0;v<G_nvar;v++) (*c)->P[v] = (*c)->Ph[v];
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



int UniformGridParallel::BC_select_data2send(
      list<cell *> *l,
      int *nc,
      boundary_data *b
      )
{
  // Check inputs.
  if ( !(*l).empty() ) {
#ifdef TESTING
    rep.warning("BC_select_data2send: List not empty! Emptying it now.",0,(*l).size());
#endif
    (*l).clear(); if(!(*l).empty()) rep.error("emptying list.",(*l).empty());
  }
  if ( *nc !=0 ) {
#ifdef TESTING
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
    for (int v=0;v<BC_nbc;v++) temp = NextPt(temp,b->ondir);
    (*l).push_back(temp);
    count++;
    ++c;
  } while (c!=b->data.end());

#ifdef TESTING
  cout <<"Got "<<count<<" cells, expected "<<b->data.size();
  cout <<"  list size = "<< (*l).size() <<"\n";
#endif
  *nc = count;


  return 0;
}




// ##################################################################
// ##################################################################




/***************************************************************************************/
/***************************** RAYTRACING BOUNDARY STUFF *******************************/
/***************************************************************************************/
#ifdef PLLEL_RT


// ##################################################################
// ##################################################################



int UniformGridParallel::Setup_RT_Boundaries(
      const int src_id,
      struct rad_sources &RS
      )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::Setup_RT_Boundaries() starting!\n";
#endif 

  // ----------------------- SANITY CHECKS --------------------------
  // processor ids are numbered from zero; if they're not (by MPI implementation)
  // then my code will bug out before it gets to here...
  //if (!RT_recv_list.empty() || !RT_send_list.empty())
  //  rep.error("send/recv lists not empty!!! more than one source???",RT_recv_list.size());

  //
  // Check that send/recv lists for this source have not already been set up!
  // Element i of the list should NOT correspond to source i, unless by chance.
  // Look through the ids in RT_source_list, just to make sure it is not already there.
  //
  int sle=-1; // source list element.
  for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
    if (RT_source_list[v].source_id == src_id) sle = v;
  }
  if (sle>=0)
    rep.error("source is already present in RT_source_list!",sle);

  //
  // Check that we are not past the correct number of sources!
  //
  if (src_id >= RS.Nsources) {
      rep.error("Requested RT source which does not exist",
                src_id-RS.Nsources);
  }
  // ----------------------- SANITY CHECKS --------------------------

  // ------------------------ SETUP BOUNDARIES ----------------------
  //
  // Now we need to add a new send list and a new recv list for this source.
  //
  struct RT_source_comms_info  this_src_comms;
  this_src_comms.source_id = src_id;
  this_src_comms.RT_recv_list.clear();
  this_src_comms.RT_send_list.clear();
  
  //
  // Now we call the appropriate function, depending on what type of source it is
  // The only two relevant types of source, as far as boundary communication is
  // concerned, are sources at infinity and those at a finite position, because the
  // boundary communication differs between them.
  //
  int err=0;
  if (!RS.sources[src_id].at_infinity) {
    err += setup_RT_finite_ptsrc_BD(
          this_src_comms.source_id, RS,
          this_src_comms.RT_recv_list,
          this_src_comms.RT_send_list
          );
  }
  else {
    err += setup_RT_infinite_src_BD(
          this_src_comms.source_id, RS,
          this_src_comms.RT_recv_list,
          this_src_comms.RT_send_list
          );
  }

  //
  // Now we add this_src_comms to the list.
  //
  RT_source_list.push_back(this_src_comms);

  // ------------------------ SETUP BOUNDARIES ----------------------

#ifdef RT_TESTING
  cout <<"UniformGridParallel::Setup_RT_Boundaries() finished!\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



int UniformGridParallel::setup_RT_infinite_src_BD(
      const int src_id,
      struct rad_sources &RS,
      std::vector<struct RT_boundary_list_element>  &RT_recv_list,
      std::vector<struct RT_boundary_list_element>  &RT_send_list
      )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_diffuse_radn_src_BD() starting.\n";
#endif 
  int err=0;
  //
  // Make sure src is at infinity, and has the correct type.
  //
  if (!RS.sources[src_id].at_infinity) {
    rep.error("Source for diffuse radiation is not at infinity",src_id);
  }
  //
  // If the previous function returned at all, then the source
  // exists, so no need to check that.  For a source at infinity
  // there can be only one send and one receive domain (or none, if
  // current domain an edge domain).
  //
  
  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
    rep.error("diffuse radn src recv-list not empty!",
              RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
    rep.error("diffuse radn src send-list not empty!",
              RT_send_list.size());
  }

  //
  // Get source direction:
  //
  enum direction srcdir = RT_src_at_infty_direction(src_id, RS);
  if (srcdir==NO) {
    rep.error("src position is not at infinity",srcdir);
  }
  enum direction recv_dir = srcdir;
  enum direction send_dir = OppDir(srcdir);
  //
  // Recv boundary in the source direction, send boundary in opp.dir.
  // Can have only one send and one recv source, so the lists should
  // get exactly one element.
  //
  int recv_proc = mpiPM->ngbprocs[recv_dir];
  int send_proc = mpiPM->ngbprocs[send_dir];

  //
  // Set up the lone recv boundary list element.  Note we leave an
  // empty vector if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempR;
  tempR.RT_bd=0;
  switch (recv_dir) {
    case XN: tempR.dir = static_cast<int>(dir_XN); break;
    case XP: tempR.dir = static_cast<int>(dir_XP); break;
    case YN: tempR.dir = static_cast<int>(dir_YN); break;
    case YP: tempR.dir = static_cast<int>(dir_YP); break;
    case ZN: tempR.dir = static_cast<int>(dir_ZN); break;
    case ZP: tempR.dir = static_cast<int>(dir_ZP); break;
    default:
    cout <<"\t No processor in receive direction d="<<recv_dir;
    cout <<": proc="<<recv_proc<<".\n";
    tempR.dir=-1;
    break;
  }
  tempR.rank = recv_proc;
  if (tempR.rank >=0) {
#ifdef RT_TESTING
    cout <<"\t\tFound diffuse-recv-proc in dir="<<recv_dir;
    cout <<" (R.dir="<<tempR.dir<<")";
    cout <<", rank="<<tempR.rank<<", setting up boundary data.\n";
#endif 
    RT_recv_list.push_back(tempR);
    err = setup_RT_recv_boundary(RT_recv_list.front());
    if (err) 
      rep.error("failed to set up diffuse radn recv boundary",err);
  }
  else {
#ifdef RT_TESTING
    cout <<"\t\tDiffuse-recv-proc doesn't exist in dir="<<recv_dir;
    cout <<"; I must be at edge of domain.\n";
#endif 
  }

  //
  // Set up the lone send boundary list element.  Note we leave an
  // empty vector if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempS;
  tempS.RT_bd=0;
  switch (send_dir) {
    case XN: tempS.dir = static_cast<int>(dir_XN); break;
    case XP: tempS.dir = static_cast<int>(dir_XP); break;
    case YN: tempS.dir = static_cast<int>(dir_YN); break;
    case YP: tempS.dir = static_cast<int>(dir_YP); break;
    case ZN: tempS.dir = static_cast<int>(dir_ZN); break;
    case ZP: tempS.dir = static_cast<int>(dir_ZP); break;
    default:
#ifdef RT_TESTING
    cout <<"\t No processor in send direction d="<<send_dir;
    cout <<": proc="<<send_proc<<".\n";
#endif 
    tempS.dir=-1;
    break;
  }
  tempS.rank = send_proc;
  if (tempS.rank >=0) {
#ifdef RT_TESTING
    cout <<"\t\tFound diffuse-send-proc in dir="<<send_dir;
    cout <<", rank="<<tempS.rank<<", setting up boundary data.\n";
#endif 
    RT_send_list.push_back(tempS);
    err = setup_RT_send_boundary(RT_send_list.front());
    if (err) rep.error("failed to set up diffuse radn send boundary",err);
  }
  else {
#ifdef RT_TESTING
    cout <<"\t\tDiffuse-send-proc doesn't exist in dir="<<send_dir;
    cout <<"; I must be at edge of domain.\n";
#endif 
  }

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_diffuse_radn_src_BD() returning.\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



enum direction UniformGridParallel::RT_src_at_infty_direction(
      const int src_id,
      struct rad_sources &RS
      )
{
  //
  // Get source position vector; compare values to find the
  // unreasonable value which identifies the direction of the
  // infinite source.
  //
  double srcpos[MAX_DIM];
  for (int v=0;v<G_ndim;v++) {
    srcpos[v] = RS.sources[src_id].pos[v];
  }
  enum direction srcdir=NO;
  for (int v=0;v<G_ndim;v++) {
    if (srcpos[v] < -1.e99) srcdir = static_cast<direction>(2*v);
    if (srcpos[v] >  1.e99) srcdir = static_cast<direction>(2*v+1);
  }
  if (srcdir==NO) {
    cout <<"WARNING: RT source is not at infinity! RT_src_at_infty_direction()\n";
    rep.printVec("Srcpos",srcpos,G_ndim);
  }
  return srcdir;
}



// ##################################################################
// ##################################################################




int UniformGridParallel::setup_RT_finite_ptsrc_BD(
      const int src_id,
      struct rad_sources &RS,
      std::vector<struct RT_boundary_list_element>  &RT_recv_list,
      std::vector<struct RT_boundary_list_element>  &RT_send_list
      )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_finite_ptsrc_BD() starting.\n";
#endif
  //
  // NOTE THIS FUNCTION IS OLD CODE, AND PROBABLY NOT THE MOST EFFICIENT WAY OF 
  // DOING THIS!  IN PARTICULAR, IT IS DEFINITELY NOT MODULAR IN ANY WAY, AND 
  // IS VERY COMPLICATED, AND HAD LOTS OF BUGS -- TOOK ME A MONTH AT LEAST TO 
  // FIX ALL OF THEM.  SO I DIDN'T WANT TO TOUCH IT WHEN UPGRADING -- I JUST 
  // WRAPPED IT IN THIS FUNCTION.
  //

  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
          rep.error("Monochromatic point src recv-list not empty!",
                    RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
          rep.error("Monochromatic point src send-list not empty!",
                    RT_send_list.size());
  }

  //
  // First we find the source:
  //
  // Source should be at a cell corner, so we need to be careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processors!), so
  // this should work fine.
  //
  enum direction srcdir[G_ndim];
  double srcpos[MAX_DIM];
  for (int v=0;v<G_ndim;v++) {
    srcpos[v] = RS.sources[src_id].pos[v];
  }


  int i_srcpos[G_ndim];
  CI.get_ipos_vec(srcpos,i_srcpos);

  for (int i=0;i<G_ndim;i++) {
    if      (i_srcpos[i]<G_ixmin[i]) {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin=";
      cout <<G_xmin[i]<<" src is off grid in neg.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i);
    }
    else if (i_srcpos[i]>G_ixmax[i])  {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmax=";
      cout <<G_xmax[i]<<" src is off grid in pos.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i+1);
    }
    else {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin=";
      cout <<G_xmin[i]<<" xmax="<<G_xmax[i];
      cout <<" src is on local domain.\n";
#endif
      srcdir[i] = static_cast<direction>(-1);
    }
  }

  //
  // Choose processors I need to receive data from and send data to.
  //
  int nx[G_ndim];
  for (int i=0;i<G_ndim;i++)
    nx[i] =static_cast<int>(ONE_PLUS_EPS*Sim_range[i]/G_range[i]);

  //
  // First see if direction of source off grid has a neigbour processor and
  // if the opposite direction has a neighbour processor.
  //
  // recv procs:
  bool recv_proc_exists[G_ndim];
  for (int i=0;i<G_ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if      ((srcdir[i]==negdir) && (!pconst.equalD(G_xmin[i],Sim_xmin[i])))
      recv_proc_exists[i]=true;
    else if ((srcdir[i]==posdir) && (!pconst.equalD(G_xmax[i],Sim_xmax[i])))
      recv_proc_exists[i]=true;
    else
      recv_proc_exists[i]=false;
  }

  //
  // send procs:
  //
  // Source should be at a cell corner, so we need to be more careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processors!), so
  // this should work fine.
  //
  // So we only send data to a neighbour if no part of it is in a
  // source plane.
  //
  bool send_proc_exists[2*G_ndim];
  for (int i=0;i<G_ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if ( (i_srcpos[i]>G_ixmin[i]) && (!pconst.equalD(G_xmin[i],Sim_xmin[i])) )
      send_proc_exists[negdir] = true;
    else 
      send_proc_exists[negdir] = false; // either doesn't exist, or we don't need it.

    if ( (i_srcpos[i]<G_ixmax[i]) && (!pconst.equalD(G_xmax[i],Sim_xmax[i])) )
      send_proc_exists[posdir] = true;
    else 
      send_proc_exists[posdir] = false; // either doesn't exist, or we don't need it.
  }
#ifdef RT_TESTING
  rep.printVec("send_proc_exists",send_proc_exists,2*G_ndim);
#endif

  //
  // Now get id's and directions of receive processors.
  //
  struct RT_boundary_list_element temp;  temp.RT_bd=0;
  // x-dir
  if      (srcdir[XX]==XN && recv_proc_exists[XX]==true) {
    temp.rank=mpiPM->get_myrank()-1; temp.dir = dir_XN;
    RT_recv_list.push_back(temp);
  }
  else if (srcdir[XX]==XP && recv_proc_exists[XX]==true) {
    temp.rank=mpiPM->get_myrank()+1; temp.dir = dir_XP;
    RT_recv_list.push_back(temp);
  }

  // y-dir
  if (G_ndim>1) {
    if (srcdir[YY]==YN && recv_proc_exists[YY]==true) {
      temp.rank=mpiPM->get_myrank()-nx[XX]; temp.dir = dir_YN;
      RT_recv_list.push_back(temp);
    }
    else if (srcdir[YY]==YP && recv_proc_exists[YY]==true) {
      temp.rank=mpiPM->get_myrank()+nx[XX]; temp.dir = dir_YP;
      RT_recv_list.push_back(temp);
    }
  } // y-dir

  // z-dir
  if (G_ndim>2 && recv_proc_exists[ZZ]==true && srcdir[ZZ]!=NO) {
    int rank; int dir;
    if       (srcdir[ZZ]==ZN) {
      rank = mpiPM->get_myrank()-nx[XX]*nx[YY]; dir = dir_ZN;
    }
    else if (srcdir[ZZ]==ZP) {
      rank = mpiPM->get_myrank()+nx[XX]*nx[YY]; dir = dir_ZP;
    }
    else {
      rank=999; dir=999;
      rep.error("Bad Logic -- z-dir RT setup",srcdir[ZZ]);
    }
    temp.rank=rank; temp.dir=dir; temp.RT_bd=0;
    RT_recv_list.push_back(temp);
  } //z-dir
  // Up to Seven processors to receive from, or as few as None.

  //
  // Choose processors I need to send data to.
  //
  // x-dir
  int rank = mpiPM->get_myrank();
  int dir  = 0;
  if (send_proc_exists[XN]) {
    temp.rank = rank-1;
    temp.dir = dir_XN;
    RT_send_list.push_back(temp);
  }
  if (send_proc_exists[XP]) {
    temp.rank = rank+1;
    temp.dir = dir_XP;
    RT_send_list.push_back(temp);
  }

  // y-dir
  if (G_ndim>1) {
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[YN]) {
      rank -= nx[XX]; dir=dir_YN;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
    }
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[YP]) {
      rank += nx[XX]; dir=dir_YP;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
    }
  }
  // z-dir
  if (G_ndim>2) {
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[ZN]) {
      rank -= nx[XX]*nx[YY];
      dir = dir_ZN;
      temp.rank = rank;
      temp.dir=dir;
      RT_send_list.push_back(temp);
    }
    rank=mpiPM->get_myrank(); dir=0;
    if (send_proc_exists[ZP]) {
      rank += nx[XX]*nx[YY];
      dir = dir_ZP;
      temp.rank = rank;
      temp.dir=dir;
      RT_send_list.push_back(temp);
    }
  } // z-dir

  //
  // Figure out how many RT boundaries I need and set them up.
  //
#ifdef RT_TESTING
  cout <<"send list: "<<RT_send_list.size()<<" boundaries.\n";
  cout <<"recv list: "<<RT_recv_list.size()<<" boundaries.\n";
  cout <<"SEND LIST DIRECTIONS: ";
  for (unsigned int i=0;i<RT_send_list.size();i++) cout <<RT_send_list[i].dir<<" ";
  cout <<"\n";
  cout <<"RECV LIST DIRECTIONS: ";
  for (unsigned int i=0;i<RT_recv_list.size();i++) cout <<RT_recv_list[i].dir<<" ";
  cout <<"\n";
#endif

  //
  // Set up receive boundaries.  First initialise all boundary_data pointers
  // to zero, then setup faces, edges, and corners, in that order (so ngb 
  // pointers are set up to get to corners.
  //
  int err=0;
  for (unsigned int i=0;i<RT_recv_list.size();i++)
    RT_recv_list[i].RT_bd=0;
  //
  // First set up face boundaries.
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    int d=RT_recv_list[i].dir;
    if (d==dir_XN || d==dir_XP ||
        d==dir_YN || d==dir_YP ||
        d==dir_ZN || d==dir_ZP) {
      err += setup_RT_recv_boundary(RT_recv_list[i]);
    }
  }
  //
  // Make sure we initialised all the boundaries correctly!
  //
  for (unsigned int i=0;i<RT_recv_list.size();i++) {
    if (!(RT_recv_list[i].RT_bd)) {
      rep.error("Missed out on receive boundary!",RT_recv_list[i].dir);
    }
  }
  if (err) rep.error("failed to setup RT recv boundaries",err);

  for (unsigned int i=0;i<RT_send_list.size();i++) {
    RT_send_list[i].RT_bd=0;
    err += setup_RT_send_boundary(RT_send_list[i]);
  }
  for (unsigned int i=0;i<RT_send_list.size();i++) {
    if (!(RT_send_list[i].RT_bd)) {
      rep.error("Missed out on send boundary!",i);
    }
  }
  if (err) rep.error("failed to setup RT send boundaries",err);

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_monochromatic_ptsrc_BD() finished!\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################




int UniformGridParallel::setup_RT_recv_boundary(
      struct RT_boundary_list_element &b ///< boundary info.
      )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_recv_boundary() starting (dir="<<b.dir<<").\n";
#endif 
  int err=0;

  //
  // First set up boundary data, since it starts as an uninitialised pointer.
  //
  if (b.RT_bd) rep.error("Already set up RT boudary data!",b.RT_bd);
  b.RT_bd = mem.myalloc(b.RT_bd,1);
  enum direction d1;

  switch (b.dir) {
    //
    // First do the faces, which already have boundaries set up, so we only need
    // to connect the boundary cells to each other laterally.
    //
  case dir_XN:
    d1=XN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_XP:
    d1=XP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_YN:
    d1=YN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_YP:
    d1=YP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_ZN:
    d1=ZN;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;
  case dir_ZP:
    d1=ZP;
    if (!BC_bd) rep.error("setup regular boundaries before RT ones!!!",b.dir);
    err += RT_populate_recv_boundary(b.RT_bd, &(BC_bd[d1]), d1);
    break;

  default:
    rep.error("bad dir!",b.dir);
  }

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_recv_boundary() returning (dir="<<b.dir<<").\n";
#endif 
  return err;
}


// ##################################################################
// ##################################################################



int UniformGridParallel::RT_populate_recv_boundary(
        struct boundary_data *b,        ///< pointer to RT boundary data.
        const struct boundary_data *b2, ///< pointer to BC boundary data.
        const enum direction offdir     ///< face dir
        )
{
  //if (G_ndim==1) return 0; // no connections to make!
  if (!b || !b2)
    rep.error("UniformGridParallel::RT_populate_recv_boundary() Null b/b2 pointer",b);

#ifdef RT_TESTING
    cout <<"RT_populate_recv_boundary() offdir="<< offdir;
    cout <<",  ondir="<<b2->ondir<<"\n";
#endif // RT_TESTING

  //
  // Run through all the BC boundary points (may be 2 deep, so skip
  // 2nd ones).  Add 1st boundary cells to the RT boundary list.
  //
  list<cell*>::const_iterator bpt=b2->data.begin();
  do{
    //
    // if isedge==-1 then cell is one cell off-grid in offdir.
    //
    if ((*bpt)->isedge == -1) {
      //
      // Add to RT list.
      //
      b->data.push_back(*bpt);
#ifdef RT_TESTING
      cout <<"RT_populate_recv_boundary() cpos="<< (*bpt)->pos[0]<<"\n";
#endif // RT_TESTING
    }
    //
    // Move to next cell.
    //
    ++bpt;
  } while (bpt !=b2->data.end());
  return 0;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::setup_RT_send_boundary(
      struct RT_boundary_list_element &send_b ///< boundary info.
      )
{
#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_send_boundary() ";
  cout <<"starting (dir="<<send_b.dir<<").\n";
#endif 
  int err=0;
  //
  // get a pointer to the existing grid boundary.
  //
  boundary_data *grid_b = &(BC_bd[send_b.dir]);
  if (!grid_b) {
    rep.error("RT boundary in dir with no real boundary!",send_b.dir);
  }

  //
  // Set up boundary data.
  //
  if (send_b.RT_bd) {
    rep.error("Already set up RT boudary data!",send_b.RT_bd);
  }
  send_b.RT_bd = mem.myalloc(send_b.RT_bd,1);
  
  //
  // Now get every cell in the grid boundary for which isedge==-1,
  // indicating that it is one cell off-grid, and add the neighbour
  // cell in the on-grid-direction to the send-boundary list.
  //
  list<cell*>::const_iterator bpt=grid_b->data.begin();
  do{
    //
    // if isedge==-1 then cell is one cell off-grid in offdir.
    //
    if ((*bpt)->isedge == -1) {
      //
      // Add to RT list.
      //
      send_b.RT_bd->data.push_back(NextPt(*bpt,grid_b->ondir));
#ifdef RT_TESTING
      cout <<"setup_RT_send_boundary() cpos="<< (*bpt)->pos[0];
      if (G_ndim>1) cout <<", "<<(*bpt)->pos[1]<<"\n";
#endif // RT_TESTING
    }
    //
    // Move to next cell.
    //
    ++bpt;
  } while (bpt !=grid_b->data.end());

#ifdef RT_TESTING
  cout <<"UniformGridParallel::setup_RT_send_boundary() returning ";
  cout <<"(dir="<<send_b.dir<<").\n";
#endif 
  return err;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::Receive_RT_Boundaries(
      const int src_id, ///< source id
      struct rad_sources &RS
      )
{
#ifdef RT_TESTING
  cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": Starting\n";
#endif 
  int err=0;
  
  //
  // We need the recv-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle=-1; // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
   cout <<"\tRecv_RT_Boundaries() src="<<src_id<<": no procs to recv from.\n";
#endif 
   return 0;
  }
  else {
    for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle<0) rep.error("source not found in RT_source_list",sle);
#ifdef RT_TESTING
  else
   cout <<"\tRevcv_RT_Boundaries() src="<<src_id<<": receiving data.\n";
#endif 

//  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
//    rep.error("Can't receive boundary for src_id b/c boundary was never set up",src_id);
//  }
//  if (RT_source_list[src_id].source_id != src_id) {
//    rep.error("source id for boundary indexed by src_id is wrong",src_id);
//  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);

  //
  // See how many boundaries to receive, and create a list of them
  //
  int n_recv = i_src->RT_recv_list.size();
  if (n_recv==0) {
#ifdef RT_TESTING
    cout <<"UniformGridParallel::Receive_RT_Boundaries() nothing to receive... returning.\n";
#endif 
  }
  else {
    int  r_dirs[n_recv]; // List of directions.
    bool r_done[n_recv]; // Will set these to true as I receive them.
#ifdef RT_TESTING
    cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": to recv from: ";
#endif 
    for (int i=0;i<n_recv;i++) {
      r_dirs[i] = i_src->RT_recv_list[i].dir;
      r_done[i] = false;
#ifdef RT_TESTING
      cout <<" "<<r_dirs[i]<<"[rank="<<i_src->RT_recv_list[i].rank<<"]";
#endif 
    }
#ifdef RT_TESTING
    cout <<"\n";
#endif 
    
    //
    // Now receive the n_recv boundaries, in any order:
    //
    double *buf=0;
    struct boundary_data *b;
    int from_rank, comm_tag, ct, i_recv;
    string id;
    for (int i=0;i<n_recv;i++) {
      id.erase();
      //
      // Find a message that has been sent to us.
      //
      err = COMM->look_for_data_to_receive(&from_rank, id, &comm_tag, COMM_DOUBLEDATA);
      if (err) rep.error("Receive_RT_Boundaries() look4data returned abnormally",err);
      if (comm_tag!=BC_RTtag) {
        rep.error("Getting data, but not RT and waiting for RT data!",comm_tag);
      }

      // 
      // Set pointer to boundary we are receiving, and note which one it is in i_recv.
      //
      b=0; bool test=false; i_recv=-1;
      for (int v=0; v<n_recv; v++) {
        if (i_src->RT_recv_list[v].rank==from_rank) {
          b=i_src->RT_recv_list[v].RT_bd;
          i_recv = v;
          if (!test) test=true;
          else rep.error("Got two boundaries with same fromrank!!!",from_rank);
        }
      }
      if (!b) rep.error("No Receive Boundary corresponding to data from this proc!",from_rank);
      if (r_dirs[i_recv] != i_src->RT_recv_list[i_recv].dir) {
        rep.error("direction mismatch",i_recv); // this is paranoid checking!
      }
      if (r_done[i_recv]) rep.error("This boundary already received!!!",i_recv);
      
      //
      // See how many doubles we are expecting, and allocate memory.
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      ct = b->data.size()*RS.sources[src_id].NTau;
      if (ct<1) {
        cerr <<"data size = "<< b->data.size();
        cerr <<", NTau = "<<RS.sources[src_id].NTau<<"\n";
        rep.error("Empty boundary!",ct);
      }
      buf=0;
      buf = mem.myalloc(buf,ct);
#ifdef RT_TESTING
      cout <<"\tReceive_RT_Boundaries() rad_src="<<src_id<<": getting "<<ct<<" cells from proc. "<<from_rank<<"\n";
#endif
      
      //
      // Receive data into buffer.
      //
      err = COMM->receive_double_data(from_rank, ///< rank of process we are receiving from.
              comm_tag,  ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
              id, ///< identifier for receive, for any book-keeping that might be needed.
              ct, ///< number of doubles to receive
              buf ///< Pointer to array to write to (must be already initialised).
              );
      if (err) rep.error("Receive_RT_Boundaries() getdata failed",err);

      //
      // Loop over cells in boundary and copy received col values into c->col.
      // Hopefully they are ordered the same way in packing as unpacking!!!
      //
      list<cell*>::iterator c=b->data.begin();
      int count=0;
      double tau[MAX_TAU];
      for (short unsigned int v=0; v<MAX_TAU; v++) tau[v]=0.0;

      for (c=b->data.begin(); c!=b->data.end(); ++c) {
        if (count>=ct) rep.error("too many cells!!!",count-ct);
#ifdef RT_TESTING
        if (count<32) {
          cout <<"col = "<<buf[count]<<" for cell "<<count<<": cell-id="<<(*c)->id;
                cout <<", pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (G_ndim>2) cout<<","<<(*c)->pos[ZZ]<<"]"<<"\n";
          else cout <<"]\n";
        }
#endif  
        for (short unsigned int v=0; v<RS.sources[src_id].NTau; v++) {
          tau[v] = buf[count];
          count++;
        }
        CI.set_col(*c, src_id, tau);
      }
      if (count != ct) rep.error("BIG ERROR!",count-ct);
      
      //
      // Note that this boundary has been received, so we know if we get the same data twice.
      //
      r_done[i_recv] = true;
      buf = mem.myfree(buf);
      if (err) rep.error("Error receiving boundary data from rank: ",from_rank);
    } // Loop over boundaries.
  } // else have boundaries.

#ifdef RT_TESTING
  cout <<"\tReceive_RT_Boundaries() src="<<src_id<<": Finished\n";
#endif 
  return err;
}




// ##################################################################
// ##################################################################




int UniformGridParallel::Send_RT_Boundaries(
      const int src_id, ///< source id
      struct rad_sources &RS
      )
{
  int err=0;
#ifdef RT_TESTING
  cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Starting\n";
#endif 
  //
  // We need the send-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle=-1; // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
   cout <<"\tSend_RT_Boundaries() src="<<src_id<<": no procs to send to.\n";
#endif 
   return 0;
  }
  else {
    for (int v=0;v<static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle<0) rep.error("source not found in RT_source_list",sle);
#ifdef RT_TESTING
  else
   cout <<"\tSend_RT_Boundaries() src="<<src_id<<": sending data.\n";
#endif 

  //  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
  //    rep.error("Can't send boundary for src_id b/c boundary was never set up",src_id);
  //  }
  //  if (RT_source_list[src_id].source_id != src_id) {
  //    rep.error("source id for boundary indexed by src_id is wrong",src_id);
  //  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);


  //
  // See how many boundaries to send:
  // If none, just skip on to the COMM->barrier() call at the end.
  // If we have some, do the sends.
  //
  int n_send = i_src->RT_send_list.size();
  if (n_send==0) {
#ifdef RT_TESTING
    cout <<"UniformGridParallel::Send_RT_Boundaries() nothing to send... returning.\n";
#endif 
  }
  else {
    //
    // Allocate memory for send buffers.
    //
#ifdef RT_TESTING
    cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Allocating "<<n_send<<" data buffers.\n";
#endif 
    double *data=0;
    std::string *id=0;
    id = mem.myalloc(id, n_send);
    //
    // Loop over Boundaries, put c->col into data, send it.
    // I'm going to use non-blocking sends, so they all send at once.
    //
    for (int i=0; i<n_send; i++) {
      struct boundary_data *b = i_src->RT_send_list[i].RT_bd;
      data=0;

#ifdef RT_TESTING
      cout <<"\tSend_RT_Boundaries() src="<<src_id;
      cout <<": Sending boundary in dir "<<i_src->RT_send_list[i].dir<<"\n";
#endif 
      //
      // How many cell column densities to send:
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      int nc = b->data.size()*RS.sources[src_id].NTau;
      data = mem.myalloc(data, nc);

      //
      // Put columns into data
      //
      list<cell*>::iterator c=b->data.begin();
      int count=0;
      double tau[MAX_TAU];
      for (short unsigned int v=0; v<MAX_TAU; v++) tau[v]=0.0;
      
      for (c=b->data.begin(); c!=b->data.end(); ++c) {
        if (count>=nc) rep.error("too many cells!!!",count-nc);

        CI.get_col(*c, src_id, tau);
        for (short unsigned int v=0;
             v<RS.sources[src_id].NTau; v++) {
          data[count] = tau[v];
          count++;
        }
#ifdef RT_TESTING
        if (count<32) {
          cout <<"send data ["<<i<<"]: col[0] = ";
          CI.get_col(*c,src_id,tau);
          cout << tau[0] <<" for cell ";
          cout <<count<<": pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (G_ndim>2)cout<<","<<(*c)->pos[ZZ]<<"]\n";
          else cout <<"]\n";
        }
        //if (count<32) {
        //  cout <<"send data ["<<i<<"]: col = "<<(*c)->col<<" for cell ";
        //  cout <<count<<": "; rep.printVec("pos",(*c)->pos, G_ndim);
        //}
#endif 
      }

      //
      // Send data
      //
#ifdef RT_TESTING
      cout <<"Send_BC["<<i<<"]: Sending "<<nc<<" doubles to proc: ";
      cout <<i_src->RT_send_list[i].rank;
      cout <<" in direction "<<i_src->RT_send_list[i].dir<<"\n";
#endif
      err += COMM->send_double_data(
            i_src->RT_send_list[i].rank, ///< rank to send to.
            nc,      ///< size of buffer, in number of doubles.
            data,    ///< pointer to double array.
            id[i],   ///< identifier for send, for tracking delivery later.
            BC_RTtag ///< comm_tag, to say what kind of send this is.
            );
      if (err) rep.error("Send_BC[] send_data failed.",err);
#ifdef RT_TESTING
      cout <<"Send_BC["<<i<<"]: send returned with id="<<id[i]<<"\n";
#endif
    //
    // Delete data array for reuse in next send.
    //
    data = mem.myfree(data);
    } // loop over all sends to send data.

    //
    // Loop over all sends, and wait until they complete!
    //
#ifdef RT_TESTING
    cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Waiting for sends to complete.\n";
#endif 
    for (int i=0; i<n_send; i++) {
      err += COMM->wait_for_send_to_finish(id[i]);
    }

    //
    // delete array of ids
    //
    id = mem.myfree(id);
  } // else (have boundaries to send).

  //
  // Wait for all processors to finish RT for this source, so no messages get mixed up!
  // Not sure if this is strictly necessary... so should try to do without it, once I 
  // know the algorithm is sound.
  //
#ifdef RT_TESTING
  cout <<"\tSend_RT_Boundaries() src="<<src_id<<": Waiting for all procs to finish.\n";
#endif
  COMM->barrier("Send_RT_Boundaries_finished");
  return err;
}


// ##################################################################
// ##################################################################



#endif // PLLEL_RT
//-------------------------------------------------------------
//------------------- CARTESIAN GRID END ----------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID START ------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_cyl_parallel::uniform_grid_cyl_parallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn,
      double *xp,
      int *nc,
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp, ///< array of max. x/y/z for full simulation.
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp,p),
  VectorOps_Cyl(nd),
  uniform_grid_cyl(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel:: cyl. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim>2)
    rep.error("Need to write code for 3 dimensions",G_ndim);

#ifdef TESTING
  cout <<"cylindrical grid: dR="<<G_dx<<"\n";
#endif 
  set_dx(G_dx);
}



// ##################################################################
// ##################################################################


uniform_grid_cyl_parallel::~uniform_grid_cyl_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_cyl_parallel destructor. Present and correct!\n";
#endif 
}



// ##################################################################
// ##################################################################


double uniform_grid_cyl_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  //cout <<" Cell radius: "<< R_com(c)/CI.phys_per_int() +G_ixmin[Rcyl];
  //rep.printVec("  cell centre",c->pos,G_ndim);
  return (R_com(c)-Sim_xmin[Rcyl])/CI.phys_per_int() +Sim_ixmin[Rcyl];
}



// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID END --------------------
//-------------------------------------------------------------



//-------------------------------------------------------------
//------------------- SPHERICAL GRID START --------------------
//-------------------------------------------------------------


// ##################################################################
// ##################################################################


///
/// Constructor
///
uniform_grid_sph_parallel::uniform_grid_sph_parallel(
      int nd,
      int nv,
      int eqt,
      int nbc, ///< number of boundary cells to use.
      double *xn,
      double *xp,
      int *nc,
      double *sim_xn, ///< array of min. x/y/z for full simulation.
      double *sim_xp, ///< array of max. x/y/z for full simulation.
      class MCMDcontrol *p
      )
  :
  VectorOps_Cart(nd),
  UniformGrid(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp),
  UniformGridParallel(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp,p),
  VectorOps_Cyl(nd),
  VectorOps_Sph(nd),
  uniform_grid_sph(nd,nv,eqt,nbc,xn,xp,nc,sim_xn,sim_xp)
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel:: sph. uniform PARALLEL grid";
  cout <<" G_ndim="<<G_ndim<<" and G_nvar="<<G_nvar<<"\n";
#endif 
  if (G_ndim!=1)
    rep.error("Need to write code for >1 dimension",G_ndim);

#ifdef TESTING
  cout <<"spherical grid: dr="<<G_dx<<"\n";
#endif 
  set_dx(G_dx);
}


// ##################################################################
// ##################################################################



uniform_grid_sph_parallel::~uniform_grid_sph_parallel()
{
#ifdef TESTING
  cout <<"uniform_grid_sph_parallel destructor. Present and correct!\n";
#endif 
}


// ##################################################################
// ##################################################################



double uniform_grid_sph_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  return (R_com(c)-Sim_xmin[Rsph])/CI.phys_per_int() +Sim_ixmin[Rsph];
}


// ##################################################################
// ##################################################################


//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------
#endif //PARALLEL
