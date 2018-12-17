/// \file RT_MPI_boundaries.cpp
/// \brief Class definitions for raytracing across a parallel
///   decomposed domain, by communicating with neighbouring MPI
///   processes
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.09 JM: moved code from parallel uniform grid.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "boundaries/RT_MPI_boundaries.h"
#include <iostream>
#include <vector>
using namespace std;

//#define RT_TESTING

// ##################################################################
// ##################################################################


RT_MPI_bc::~RT_MPI_bc()
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc destructor.\n";
#endif
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

  return;
}




// ##################################################################
// ##################################################################




int RT_MPI_bc::Setup_RT_Boundaries(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid, ///< pointer to grid.
      const int src_id,
      struct rad_src_info &RS
      )
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc::Setup_RT_Boundaries() starting!\n";
#endif 

  // ----------------------- SANITY CHECKS --------------------------
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
  if (!RS.at_infinity) {
    err += setup_RT_finite_ptsrc_BD(par,ppar,grid,
          this_src_comms.source_id, RS,
          this_src_comms.RT_recv_list,
          this_src_comms.RT_send_list
          );
  }
  else {
    err += setup_RT_infinite_src_BD(par,ppar,grid,
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
  cout <<"RT_MPI_bc::Setup_RT_Boundaries() finished!\n";
#endif 
  return err;
}






// ##################################################################
// ##################################################################




int RT_MPI_bc::Receive_RT_Boundaries(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid, ///< pointer to grid.
      const int src_id, ///< source id
      struct rad_src_info &RS
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
   cout <<"\tRecv_RT_Boundaries() src="<<src_id<<": receiving data.\n";
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
    cout <<"RT_MPI_bc::Receive_RT_Boundaries() nothing to receive... returning.\n";
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
      err = COMM->look_for_data_to_receive(&from_rank, id, &comm_tag,
                                          BC_RTtag,COMM_DOUBLEDATA);
      if (err) rep.error("Receive_RT_Boundaries() look4data",err);
      if (comm_tag!=BC_RTtag) {
        rep.error("Getting data, but not RT",comm_tag);
      }

      // 
      // Set pointer to boundary we are receiving, and note which one
      // it is in i_recv.
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
      ct = b->data.size()*RS.NTau;
      if (ct<1) {
        cerr <<"data size = "<< b->data.size();
        cerr <<", NTau = "<<RS.NTau<<"\n";
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
      err = COMM->receive_double_data(
          from_rank, ///< rank of process we are receiving from.
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
        if (count<100) {
          cout <<"col = "<<buf[count]<<" for cell "<<count<<": cell-id="<<(*c)->id;
          cout <<", pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (par.ndim>2) cout<<","<<(*c)->pos[ZZ]<<"]"<<"\n";
          else cout <<"]\n";
        }
#endif  
        for (short unsigned int v=0; v<RS.NTau; v++) {
          tau[v] = buf[count];
          count++;
        }
        CI.set_col(*c, src_id, tau);
      }
      if (count != ct) rep.error("BIG ERROR!",count-ct);
      
      //
      // Record boundary as received, so we know if we get the
      // same data twice.
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




int RT_MPI_bc::Send_RT_Boundaries(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid, ///< pointer to grid.
      const int src_id, ///< source id
      struct rad_src_info &RS
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
    cout <<"RT_MPI_bc::Send_RT_Boundaries() nothing to send... returning.\n";
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
      int nc = b->data.size()*RS.NTau;
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
#ifdef RT_TESTING
        if (count<100) {
          cout <<"send data ["<<i<<"]: col[0] = ";
          CI.get_col(*c,src_id,tau);
          cout << tau[0] <<" for cell ";
          cout <<count<<": pos=["<<(*c)->pos[XX]<<","<<(*c)->pos[YY];
          if (par.ndim>2)cout<<","<<(*c)->pos[ZZ]<<"]\n";
          else cout <<"]\n";
        }
        //if (count<32) {
        //  cout <<"send data ["<<i<<"]: col = "<<(*c)->col<<" for cell ";
        //  cout <<count<<": "; rep.printVec("pos",(*c)->pos, par.ndim);
        //}
#endif 
        for (short unsigned int v=0;
             v<RS.NTau; v++) {
          data[count] = tau[v];
          count++;
        }
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



int RT_MPI_bc::setup_RT_infinite_src_BD(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid, ///< pointer to grid.
      const int src_id,
      struct rad_src_info &RS,
      std::vector<struct RT_boundary_list_element>  &RT_recv_list,
      std::vector<struct RT_boundary_list_element>  &RT_send_list
      )
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_diffuse_radn_src_BD() starting.\n";
#endif 
  int err=0;
  //
  // Make sure src is at infinity, and has the correct type.
  //
  if (!RS.at_infinity) {
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
  enum direction srcdir = RT_src_at_infty_direction(par,src_id, RS);
  if (srcdir==NO) {
    rep.error("src position is not at infinity",srcdir);
  }
  enum direction recv_dir = srcdir;
  enum direction send_dir = grid->OppDir(srcdir);
  //
  // Recv boundary in the source direction, send boundary in opp.dir.
  // Can have only one send and one recv source, so the lists should
  // get exactly one element.
  //
  int recv_proc = ppar.ngbprocs[recv_dir];
  int send_proc = ppar.ngbprocs[send_dir];

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
    err = setup_RT_recv_boundary(grid,RT_recv_list.front());
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
    err = setup_RT_send_boundary(grid,RT_send_list.front());
    if (err) rep.error("failed to set up diffuse radn send boundary",err);
  }
  else {
#ifdef RT_TESTING
    cout <<"\t\tDiffuse-send-proc doesn't exist in dir="<<send_dir;
    cout <<"; I must be at edge of domain.\n";
#endif 
  }

#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_diffuse_radn_src_BD() returning.\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



enum direction RT_MPI_bc::RT_src_at_infty_direction(
      class SimParams &par,     ///< pointer to simulation parameters
      const int src_id,
      struct rad_src_info &RS
      )
{
  //
  // Get source position vector; compare values to find the
  // very large value which identifies the direction of the
  // infinite source.
  //
  enum direction srcdir=NO;
  for (int v=0;v<par.ndim;v++) {
    if (RS.pos[v] < -1.e99) srcdir = static_cast<direction>(2*v);
    if (RS.pos[v] >  1.e99) srcdir = static_cast<direction>(2*v+1);
  }
  if (srcdir==NO) {
    cout <<"WARNING: RT source is not at infinity! RT_src_at_infty_direction()\n";
    rep.printVec("Srcpos",RS.pos,par.ndim);
  }
  return srcdir;
}



// ##################################################################
// ##################################################################




int RT_MPI_bc::setup_RT_finite_ptsrc_BD(
      class SimParams &par,     ///< pointer to simulation parameters
      class MCMDcontrol &ppar,   ///< domain decomposition info
      class GridBaseClass *grid, ///< pointer to grid.
      const int src_id,
      struct rad_src_info &RS,
      std::vector<struct RT_boundary_list_element>  &RT_recv_list,
      std::vector<struct RT_boundary_list_element>  &RT_send_list
      )
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_finite_ptsrc_BD() starting.\n";
#endif
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
  enum direction srcdir[par.ndim];
  double srcpos[MAX_DIM];
  for (int v=0;v<par.ndim;v++) {
    srcpos[v] = RS.pos[v];
  }


  int i_srcpos[par.ndim];
  CI.get_ipos_vec(srcpos,i_srcpos);

  for (int i=0;i<par.ndim;i++) {
    if      (i_srcpos[i]<grid->iXmin(static_cast<axes>(i))) {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin=";
      cout <<grid->Xmin(static_cast<axes>(i))<<" src is off grid in neg.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i);
    }
    else if (i_srcpos[i]>grid->iXmax(static_cast<axes>(i)))  {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmax=";
      cout <<grid->Xmax(static_cast<axes>(i))<<" src is off grid in pos.dir.\n";
#endif
      srcdir[i] = static_cast<direction>(2*i+1);
    }
    else {
#ifdef RT_TESTING
      cout <<"*** axis="<<i<<" srcpos="<<srcpos[i]<<" xmin=";
      cout <<grid->Xmin(static_cast<axes>(i))<<" xmax=";
      cout <<grid->Xmax(static_cast<axes>(i));
      cout <<" src is on local domain.\n";
#endif
      srcdir[i] = static_cast<direction>(-1);
    }
  }

  //
  // Choose processors I need to receive data from and send data to.
  //
  int nx[par.ndim];
  for (int i=0;i<par.ndim;i++)
    nx[i] =static_cast<int>(ONE_PLUS_EPS*grid->SIM_Range(static_cast<axes>(i))/grid->Range(static_cast<axes>(i)));

  //
  // First see if direction of source off grid has a neigbour processor and
  // if the opposite direction has a neighbour processor.
  //
  // recv procs:
  bool recv_proc_exists[par.ndim];
  for (int i=0;i<par.ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if      ((srcdir[i]==negdir) && (!pconst.equalD(grid->Xmin(static_cast<axes>(i)),grid->SIM_Xmin(static_cast<axes>(i)))))
      recv_proc_exists[i]=true;
    else if ((srcdir[i]==posdir) && (!pconst.equalD(grid->Xmax(static_cast<axes>(i)),grid->SIM_Xmax(static_cast<axes>(i)))))
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
  bool send_proc_exists[2*par.ndim];
  for (int i=0;i<par.ndim;i++) {
    enum direction posdir=static_cast<direction>(2*i+1);
    enum direction negdir=static_cast<direction>(2*i);
    if ( (i_srcpos[i]>grid->iXmin(static_cast<axes>(i))) && 
         (!pconst.equalD(grid->Xmin(static_cast<axes>(i)),
                        grid->SIM_Xmin(static_cast<axes>(i))) ) )
      send_proc_exists[negdir] = true;
    else 
      send_proc_exists[negdir] = false; // either doesn't exist, or we don't need it.

    if ( (i_srcpos[i]<grid->iXmax(static_cast<axes>(i))) &&
         (!pconst.equalD(grid->Xmax(static_cast<axes>(i)),
                         grid->SIM_Xmax(static_cast<axes>(i))) ) )
      send_proc_exists[posdir] = true;
    else 
      send_proc_exists[posdir] = false; // either doesn't exist, or we don't need it.
  }
#ifdef RT_TESTING
  rep.printVec("send_proc_exists",send_proc_exists,2*par.ndim);
#endif

  //
  // Now get id's and directions of receive processors.
  //
  struct RT_boundary_list_element temp;  temp.RT_bd=0;
  // x-dir
  if      (srcdir[XX]==XN && recv_proc_exists[XX]==true) {
    temp.rank=ppar.get_myrank()-1; temp.dir = dir_XN;
    RT_recv_list.push_back(temp);
  }
  else if (srcdir[XX]==XP && recv_proc_exists[XX]==true) {
    temp.rank=ppar.get_myrank()+1; temp.dir = dir_XP;
    RT_recv_list.push_back(temp);
  }

  // y-dir
  if (par.ndim>1) {
    if (srcdir[YY]==YN && recv_proc_exists[YY]==true) {
      temp.rank=ppar.get_myrank()-nx[XX]; temp.dir = dir_YN;
      RT_recv_list.push_back(temp);
    }
    else if (srcdir[YY]==YP && recv_proc_exists[YY]==true) {
      temp.rank=ppar.get_myrank()+nx[XX]; temp.dir = dir_YP;
      RT_recv_list.push_back(temp);
    }
  } // y-dir

  // z-dir
  if (par.ndim>2 && recv_proc_exists[ZZ]==true && srcdir[ZZ]!=NO) {
    int rank; int dir;
    if       (srcdir[ZZ]==ZN) {
      rank = ppar.get_myrank()-nx[XX]*nx[YY]; dir = dir_ZN;
    }
    else if (srcdir[ZZ]==ZP) {
      rank = ppar.get_myrank()+nx[XX]*nx[YY]; dir = dir_ZP;
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
  int rank = ppar.get_myrank();
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
  if (par.ndim>1) {
    rank=ppar.get_myrank(); dir=0;
    if (send_proc_exists[YN]) {
      rank -= nx[XX]; dir=dir_YN;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
    }
    rank=ppar.get_myrank(); dir=0;
    if (send_proc_exists[YP]) {
      rank += nx[XX]; dir=dir_YP;
      temp.rank = rank; temp.dir=dir; RT_send_list.push_back(temp);
    }
  }
  // z-dir
  if (par.ndim>2) {
    rank=ppar.get_myrank(); dir=0;
    if (send_proc_exists[ZN]) {
      rank -= nx[XX]*nx[YY];
      dir = dir_ZN;
      temp.rank = rank;
      temp.dir=dir;
      RT_send_list.push_back(temp);
    }
    rank=ppar.get_myrank(); dir=0;
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
  // to zero, then setup faces.
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
      err += setup_RT_recv_boundary(grid,RT_recv_list[i]);
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
    err += setup_RT_send_boundary(grid,RT_send_list[i]);
  }
  for (unsigned int i=0;i<RT_send_list.size();i++) {
    if (!(RT_send_list[i].RT_bd)) {
      rep.error("Missed out on send boundary!",i);
    }
  }
  if (err) rep.error("failed to setup RT send boundaries",err);

#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_monochromatic_ptsrc_BD() finished!\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################




int RT_MPI_bc::setup_RT_recv_boundary(
      class GridBaseClass *grid, ///< pointer to grid.
      struct RT_boundary_list_element &b ///< boundary info.
      )
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_recv_boundary() starting (dir="<<b.dir<<").\n";
#endif 
  int err=0;

  //
  // First set up boundary data, since it starts as an uninitialised pointer.
  //
  if (b.RT_bd) rep.error("Already set up RT boudary data!",b.RT_bd);
  b.RT_bd = mem.myalloc(b.RT_bd,1);
  enum direction d1;

  if (grid->BC_bd.size()==0)
    rep.error("setup regular boundaries before RT ones!!!",b.dir);

  switch (b.dir) {
    //
    // First do the faces, which already have boundaries set up, so we only need
    // to connect the boundary cells to each other laterally.
    //
  case dir_XN:
    d1=XN;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;
  case dir_XP:
    d1=XP;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;
  case dir_YN:
    d1=YN;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;
  case dir_YP:
    d1=YP;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;
  case dir_ZN:
    d1=ZN;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;
  case dir_ZP:
    d1=ZP;
    err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
    break;

  default:
    rep.error("bad dir!",b.dir);
  }

#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_recv_boundary() returning (dir="<<b.dir<<").\n";
#endif 
  return err;
}


// ##################################################################
// ##################################################################



int RT_MPI_bc::RT_populate_recv_boundary(
        struct boundary_data *b,        ///< pointer to RT boundary data.
        const struct boundary_data *b2, ///< pointer to BC boundary data.
        const enum direction offdir     ///< face dir
        )
{
  //if (par.ndim==1) return 0; // no connections to make!
  if (!b || !b2)
    rep.error("RT_MPI_bc::RT_populate_recv_boundary() Null b/b2 pointer",b);

#ifdef RT_TESTING
    cout <<"RT_populate_recv_boundary() offdir="<< offdir;
    cout <<",  ondir="<<b2->ondir<<"\n";
#endif // RT_TESTING

  //
  // Add all boundary cells to the RT boundary list.
  //
  list<cell*>::const_iterator bpt=b2->data.begin();
  do{
    b->data.push_back(*bpt);
#ifdef RT_TESTING
    cout <<"RT_populate_recv_boundary() cpos= ";
    rep.printVec("",(*bpt)->pos,MAX_DIM);
#endif // RT_TESTING
    ++bpt;
  } while (bpt !=b2->data.end());
  return 0;
}




// ##################################################################
// ##################################################################




int RT_MPI_bc::setup_RT_send_boundary(
      class GridBaseClass *grid, ///< pointer to grid.
      struct RT_boundary_list_element &send_b ///< boundary info.
      )
{
#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_send_boundary() ";
  cout <<"starting (dir="<<send_b.dir<<").\n";
#endif 
  int err=0;
  //
  // get a pointer to the existing grid boundary.
  //
  boundary_data *grid_b = grid->BC_bd[send_b.dir];
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
  // Now get every cell in the grid boundary, and add the 'Nbc'th
  // cell in the on-grid-direction to the send-boundary list.
  //
  list<cell*>::const_iterator bpt=grid_b->data.begin();
  cell *target=0;
  do{
    target = (*bpt);
    for (int v=0;v<grid_b->depth;v++)
      target = grid->NextPt(target,grid_b->ondir);

    send_b.RT_bd->data.push_back(target);
#ifdef RT_TESTING
    cout <<"setup_RT_send_boundary() cpos= ";
    rep.printVec(" ",target->pos,MAX_DIM);
#endif // RT_TESTING
    ++bpt;
  } while (bpt !=grid_b->data.end());

#ifdef RT_TESTING
  cout <<"RT_MPI_bc::setup_RT_send_boundary() returning ";
  cout <<"(dir="<<send_b.dir<<").\n";
#endif 
  return err;
}



// ##################################################################
// ##################################################################



// ##################################################################
// ##################################################################



