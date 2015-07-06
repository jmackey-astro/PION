/// \file comm_mpi.cc
///
/// \brief Contains comms class for multi-process communication using 
/// the Message Passing Interface (MPI).
///
/// \author Jonathan Mackey
/// \date 2009-01-23.
//
///  - 2010-02-02 JM: comments.
///
/// Modifications:\n
/// - 2010-02-03 JM: changed variable name (used 'i' twice in
///    function). comm_mpi::receive_double_data() wasn't checking
///    that I got the number of elements expected either, so I fixed
///    that. 
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.03.22 JM: removed old myalloc/myfree functions with string
///    args.
/// - 2011.06.02 JM: initialised two vars to values to get rid of
///    compiler warnings.
/// - 2012.05.15 JM: Added function for global-reduce (max/min/sum)
///    of arrays.
/// - 2013.01.17 JM: Made class less verbose; wrapped cout in TESTING
///    flags.
/// - 2015.01.26 JM: removed references to mpiPM (no longer global),
///    added COMM pointer setup, added get_rank_nproc()

#ifdef PARALLEL
#ifdef USE_MPI


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "constants.h"

#include "tools/command_line_interface.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"

#include "comms/comm_mpi.h"

#include <sstream>
using namespace std;

class comms_base *COMM = new comm_mpi ();


// ##################################################################
// ##################################################################


comm_mpi::comm_mpi()
{
#ifdef TESTING
  cout <<"*** comm_mpi constructor. ***\n";
#endif
  myrank = -1;
  nproc  = -1;
}


// ##################################################################
// ##################################################################


comm_mpi::~comm_mpi()
{
#ifdef TESTING
  cout <<"*** comm_mpi  destructor. ***\n";
#endif
}



// ##################################################################
// ##################################################################



int comm_mpi::init(int *argc,   ///< number of program arguments.
		   char ***argv ///< character list of arguments.
		   )
{
  int err = MPI_Init(argc, argv);
  if (err) {
    cerr << "error!!!" << "\n";
    MPI_Abort(MPI_COMM_WORLD, err);
  }

  err += MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  err += MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (err) {
    cerr << "Error!!! couldn't get rank/nproc.\n";
    MPI_Abort(MPI_COMM_WORLD, err);
  }

  return 0;
}



// ##################################################################
// ##################################################################



int comm_mpi::get_rank_nproc(
        int *r, ///< rank.
	int *n  ///< nproc
	)
{
  *r = myrank;
  *n = nproc;
#ifdef TESTING
  cout << "comm_mpi::get_rank_nproc():  rank: " << myrank;
  cout << " nproc: " << nproc << "\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int comm_mpi::finalise()
{
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef TESTING
  int i=-1;
  cout << "comm_mpi::finalise():  rank: ";
  MPI_Comm_rank(MPI_COMM_WORLD, &i);
  cout << i << " nproc: ";
  MPI_Comm_size(MPI_COMM_WORLD, &i);
  cout << i << "\n";
#endif
  MPI_Finalize();
  return 0;
}



// ##################################################################
// ##################################################################



int comm_mpi::abort()
{
  MPI_Abort(MPI_COMM_WORLD,999);
  return 0;
}



// ##################################################################
// ##################################################################



int comm_mpi::barrier(const std::string)
{
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}



// ##################################################################
// ##################################################################



double comm_mpi::global_operation_double(
        const string s, ///< Either Max, Min, or Sum
	const double d  ///< this process's max/min value.
	)
{
  double local=d, global=0.0;
  int err=0;

  if      (s=="MAX") {
    global = -HUGEVALUE;
    err += MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else if (s=="MIN") {
    global = HUGEVALUE;
    err += MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else if (s=="SUM") {
    global=0.0;
    err += MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else 
    rep.error("comm_mpi:global_operation_double: Bad identifier",s);
  
  if (err)
    rep.error("comm_mpi:global_operation_double: MPI call returned abnormally",err);

  return global;
}



// ##################################################################
// ##################################################################



void comm_mpi::global_op_double_array(
        const std::string s, ///< MAX,MIN,SUM
        const size_t  Nel,   ///< Number of elements in array.
	double *data         ///< pointer to this process's data array.
        )
{
  int err=0;
  double op_data[Nel];

  if      (s=="MAX") {
    //err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    err += MPI_Allreduce(static_cast<void *>(data), static_cast<void *>(op_data), Nel, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else if (s=="MIN") {
    //err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    err += MPI_Allreduce(static_cast<void *>(data), static_cast<void *>(op_data), Nel, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else if (s=="SUM") {
    //err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    err += MPI_Allreduce(static_cast<void *>(data), static_cast<void *>(op_data), Nel, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else 
    rep.error("comm_mpi:global_op_double_array: Bad identifier",s);
  
  if (err)
    rep.error("comm_mpi:global_op_double_array: MPI call returned abnormally",err);
  
  for (size_t v=0;v<Nel;v++) data[v] = op_data[v];
  return;
}




// ##################################################################
// ##################################################################



int comm_mpi::broadcast_data(
        const int sender,       ///< rank of sender.
	const std::string type, ///< Type of data INT,DOUBLE,etc.
	const int n_el,         ///< number of elements
	void *data              ///< pointer to data.
	)
{
  int err=0;
  if      (type=="DOUBLE") {
    err += MPI_Bcast(static_cast<double *>(data), n_el, MPI_DOUBLE, sender, MPI_COMM_WORLD);
  }
  else if (type=="INT") {
    err += MPI_Bcast(static_cast<int *>(data),    n_el, MPI_INT,    sender, MPI_COMM_WORLD);
  }
  else if (type=="FLOAT") {
    err += MPI_Bcast(static_cast<float *>(data),  n_el, MPI_FLOAT,  sender, MPI_COMM_WORLD);
  }
  else rep.error("Bad type of data to send",type);
  return err;
}





// ##################################################################
// ##################################################################



int comm_mpi::send_cell_data(
        const int to_rank,    ///< rank to send to.
	std::list<cell *> *l, ///< list of cells to get data from.
	long int nc,          ///< number of cells in list (extra checking!)
	string &id,           ///< identifier for send, for tracking delivery later.
	const int comm_tag    ///< comm_tag, to say what kind of send this is.
	)
{

  //
  // First initialise everything and check we have cells to get data from.
  //
  if(!id.empty()) id.erase();
  if (nc==0 || (l->empty()) ) {
    cerr <<myrank<<"\t Nothing to send to rank: "<<to_rank<<" !!!\n";
    return 1;
  }
  if (to_rank<0 || to_rank>nproc)
    rep.error("to_rank is out of bounds",to_rank);

  list<cell *>::iterator c=l->begin();
  int err=0;

  //
  // Determine size of send buffer needed
  //
  int unitsize = 
    SimPM.ndim*sizeof(CI.get_ipos(*c,0))
      + SimPM.nvar*sizeof((*c)->P[0])
	+ sizeof((*c)->id);
  long int totalsize = 0;
  totalsize = sizeof(long int) + nc*unitsize;
  //cout <<" Send buffer size= ndim*"<<sizeof(CI.get_ipos(*c,0))<<"\n";


  //
  // Allocate memory for send buffer, and for the record of the send.
  //
  char             *send_buff =0;
  struct sent_info *si        =0;
  send_buff = mem.myalloc(send_buff, totalsize);
  si = mem.myalloc(si,1);

  //
  // Pack data, using MPI.
  //
  // int MPI_Pack(void* inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm)
  //
  int position=0, ct=0, ipos[MAX_DIM];
  err += MPI_Pack(reinterpret_cast<void *>(&nc), 1, MPI_LONG, reinterpret_cast<void *>(send_buff), totalsize, &position, MPI_COMM_WORLD);
  do {
    CI.get_ipos(*c,ipos);
    err += MPI_Pack(reinterpret_cast<void *>(&((*c)->id)), 1,          MPI_INT,    reinterpret_cast<void *>(send_buff), totalsize, &position, MPI_COMM_WORLD);
    err += MPI_Pack(reinterpret_cast<void *>(ipos),        SimPM.ndim, MPI_INT,    reinterpret_cast<void *>(send_buff), totalsize, &position, MPI_COMM_WORLD);
    err += MPI_Pack(reinterpret_cast<void *>((*c)->Ph),    SimPM.nvar, MPI_DOUBLE, reinterpret_cast<void *>(send_buff), totalsize, &position, MPI_COMM_WORLD);
    ct++;
    ++c;  // next cell in list.
  } while ( c != l->end() );

  //
  // Check that we packed the right amount of data:
  //
  //  cout <<myrank<<"\tcomm_pack_send_data: bufsiz: ";
  //  cout <<totalsize<<" position in buf: "<<position<<" nc="<<nc<<" ct="<<ct<<"\n";
  if (ct != nc) rep.error("Length of list doesn't match nc",ct-nc);
  if (position > totalsize) rep.error("Overwrote buffer length",position-totalsize);
  if (err) rep.error("MPI_Pack returned abnormally",err);

  //
  // Create record of the send.
  //
  si->comm_tag = comm_tag;
  si->from_rank = myrank;
  si->to_rank = to_rank;
  si->data = reinterpret_cast<void *>(send_buff);
  si->type = COMM_CELLDATA;
  comm_mpi::sent_list.push_back(si);

  //
  // Send data. (non-blocking, so i can return and receive a boundary!!!)
  //
  //  BLOCKING:      int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  //  NON-BLOCKING: int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
  //err += MPI_Isend(reinterpret_cast<void *>(send_buff), position, MPI_PACKED, si->to_rank, si->comm_tag, MPI_COMM_WORLD, si->request);
  err += MPI_Isend(si->data, position, MPI_PACKED, si->to_rank, si->comm_tag, MPI_COMM_WORLD, &(si->request));
  if (err) rep.error("MPI_Send failed",err);
#ifdef TESTING
  cout <<"comm_mpi::send_cell_data: sending "<<position<<" bytes in buffer size: "<<totalsize<<" to rank "<<to_rank<<"\n";
#endif //TESTING

  ostringstream temp; temp.str("");
  temp <<"from "<<myrank<<" to "<<to_rank<<" tag="<<comm_tag;
  id = temp.str(); temp.str("");
  si->id = id;

  return 0;
}




// ##################################################################
// ##################################################################




int comm_mpi::wait_for_send_to_finish(
        string &id ///< identifier for the send we are waiting on.
        )
{
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() starting\n";
#endif //TESTING
  //
  // Find the send in the list of active sends.
  //
  int el=0;
  list<sent_info *>::iterator i;
  struct sent_info *si=0;

#ifdef TESTING
  cout <<"rank: "<<myrank;
  cout <<"  comm_mpi::wait_for_send_to_finish() more than one send,";
  cout <<" so finding in list.\n";
  cout <<"\t\tsend id="<<id<<"\n";
#endif //TESTING

  for (i=sent_list.begin(); i!=sent_list.end(); ++i) {
    si = (*i); el++;
    if (si->id==id) break;
  }
  if (i==sent_list.end())
    rep.error("Failed to find send with id:",id);

#ifdef TESTING
  cout <<"found send id="<<si->id<<" and looking for id="<<id<<"\n";
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() found this send.\n";
#endif //TESTING

  //
  // Use MPI to wait for the send to complete:
  //
  MPI_Status status;

#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() waiting for send.\n";
  //cout <<"request: "<<si->request<<" and address: "<<&(si->request)<<"\n";
#endif //TESTING

  int err = MPI_Wait(&(si->request), &status);
  if (err) rep.error("Waiting for send to complete!", err);

#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() send has finished\n";
#endif //TESTING

  //
  // Now free the data, and remove the send from the list.
  //
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() deleting send record\n";
#endif //TESTING

  if      (si->type==COMM_CELLDATA) {
    char *t = reinterpret_cast<char *>(si->data);
    t = mem.myfree(t);
  }
  else if (si->type==COMM_DOUBLEDATA) {
    double *t = reinterpret_cast<double *>(si->data);
    t = mem.myfree(t);
  }
  else rep.error("BAD type of data!",si->type);
  si = mem.myfree(si);
  sent_list.erase(i);


#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_mpi::wait_for_send_to_finish() returning\n";
#endif //TESTING
  return 0;
}



// ##################################################################
// ##################################################################




int comm_mpi::look_for_data_to_receive(
        int *from_rank, ///< rank of sender
        string &id,     ///< identifier for receive.
	int *comm_tag,  ///< comm_tag associated with data.
	const int type  ///< type of data we are looking for.
        )
{

#ifdef TESTING
  cout <<"rank: "<<myrank;
  cout <<"  comm_mpi::look_for_data_to_receive() starting\n";
#endif //TESTING
  //
  // Create a new received info record.
  //
  struct recv_info *ri=0;
  ri = mem.myalloc(ri,1);
  ri->to_rank   = myrank;
  ri->data = 0;
  if (type!=COMM_CELLDATA && type!=COMM_DOUBLEDATA)
    rep.error("only know two types of data to look for!",type);
  ri->type = type;

#ifdef TESTING
  cout <<"rank: "<<myrank;
  cout <<"  comm_mpi::look_for_data_to_receive() looking for source\n";
#endif //TESTING
  //
  // First Look for Data, and wait till we find some.
  // Find a message that has been sent to us.
  // int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status)
  //
  int err = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &(ri->status));
  if (err) rep.error("mpi probe failed",err);
  
#ifdef TESTING
  cout <<"rank: "<<myrank;
  cout <<"  comm_mpi::look_for_data_to_receive() found data, parsing...\n";
#endif //TESTING

  //
  // Now get the size and source of the data.
  //
  *from_rank = ri->status.MPI_SOURCE;
  *comm_tag  = ri->status.MPI_TAG;
  ostringstream temp; temp.str("");
  temp <<"from "<<*from_rank<<" to "<<myrank<<" tag="<<*comm_tag;
  id = temp.str(); temp.str("");

#ifdef TESTING
  cout <<"comm_mpi::look_for_data_to_receive: found data from ";
  cout <<*from_rank<<" with tag:"<<*comm_tag<<"\n";
#endif //TESTING
  
  //
  // Set record of where data is coming from.
  ri->id     = id;
  ri->comm_tag  = *comm_tag;
  ri->from_rank = *from_rank;
  comm_mpi::recv_list.push_back(ri);

#ifdef TESTING
  cout <<"comm_mpi::look_for_data_to_receive: returning.\n";
#endif //TESTING

  return err;
}



// ##################################################################
// ##################################################################



int comm_mpi::receive_cell_data(
        const int from_rank,  ///< rank of process we are receiving from.
        std::list<cell *> *l, ///< list of cells to get data for. 
	const long int nc,    ///< number of cells in list (extra checking!)
        const int comm_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
	const string &id      ///< identifier for receive, for book-keeping.
	)
{
  int err=0;

#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: starting.\n";
#endif //TESTING

  //
  // Find the recv in the list of active receives.
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: recv_list size="<<recv_list.size()<<"\n";
#endif //TESTING

  if (recv_list.empty()) rep.error("Call look4data before receive_data",recv_list.size());

  struct recv_info *info=0;
  list<recv_info *>::iterator iget;
  for (iget=recv_list.begin(); iget!=recv_list.end(); ++iget) {
    info = (*iget);
    if (info->id==id) break;
  }
  if (iget==recv_list.end())
    rep.error("Failed to find recv with id:",id);

#ifdef TESTING
  cout <<"found recv id="<<info->id<<" and looking for id="<<id<<"\n";
  cout <<"comm_mpi::receive_cell_data: found recv id\n";
#endif //TESTING

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank) rep.error("from_ranks don't match",from_rank-info->from_rank);
  if (comm_tag  != info->comm_tag ) rep.error("Bad comm_tag",comm_tag);
  if (id        != info->id       ) rep.error("Bad id",id);

  //
  // See how many bytes we are getting!
  // int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count) 
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: getting size of buffer.\n";
#endif //TESTING

  int ct = 0;
  err += MPI_Get_count(&(info->status), MPI_PACKED, &ct);
  if (err) rep.error("getting size of message to receive",err);

#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: buffer is "<<ct<<" bytes.\n.";
#endif //TESTING

  if (ct<0) rep.error("bad buffer size count",ct);

  //
  // Allocate memory for data.
  //
  char *buf = 0;
  buf = mem.myalloc(buf, ct);

  //
  // Now Receive the data in recv_buff
  // int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: receiving buffer of data from rank: "<<from_rank<<"\n";
#endif //TESTING

  err += MPI_Recv(buf, ct, MPI_PACKED, from_rank, comm_tag, MPI_COMM_WORLD, &(info->status));
  if (err) rep.error("MPI_Recv failed", err);
  //  cout <<"\tcomm_receive_any_data: status: "<<status.MPI_SOURCE<<", "<<status.MPI_TAG<<", ct="<<ct<<"\n";

  //  
  // First unpack the cells received count.
  //
  //   int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm)

#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: getting number of cells received\n";
#endif //TESTING

  int position=0;
  long int ncell=0;

  err = MPI_Unpack(buf, ct, &position, &ncell, 1, MPI_LONG, MPI_COMM_WORLD);

  if (err) rep.error("Unpack",err);

  if (ncell != static_cast<int>(l->size()) || (ncell !=nc)) {
    cerr <<myrank<<"\tcomm_mpi:recv: length of data = "<<l->size();
    cerr <<" or we're told it's "<<nc<<"\n";
    cerr <<myrank<<"\tcomm_mpi:recv: cells received = "<<ncell<<"\n";
    cerr <<myrank<<"\tcomm_mpi:recv: Bugging out!\n"; return(-99);
  }

#ifdef TESTING
    cout <<myrank<<"\tcomm_mpi:recv: got "<<ncell<<" cells, as expected.\n";
#endif //TESTING
  
  // 
  // Allocate arrays to unpack data into.
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: allocating arrays to unpack data into.\n";
#endif //TESTING

  int *ipos = 0;
  double *p = 0;
  ipos = mem.myalloc(ipos,SimPM.ndim);
  p = mem.myalloc(p,SimPM.nvar);

  //
  // Data should come in in the order the boundary cells are listed in.
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: unpacking data into list of cells.\n";
#endif //TESTING

  int cpos[MAX_DIM];
  int c_id=0;
  list<cell*>::iterator c=l->begin();

  for (int i=0; i<ncell; i++) {
    CI.get_ipos(*c,cpos);
    err += MPI_Unpack(buf, ct, &position, &c_id, 1,          MPI_INT, MPI_COMM_WORLD);
    err += MPI_Unpack(buf, ct, &position, ipos,  SimPM.ndim, MPI_INT, MPI_COMM_WORLD);
    err += MPI_Unpack(buf, ct, &position, p,     SimPM.nvar, MPI_DOUBLE, MPI_COMM_WORLD);
    if(err) rep.error("Unpack",err);
    // For a given boundary data iterator, put data into cells
    if (c==l->end()) rep.error("Got too many cells!",i);
    for (int v=0; v<SimPM.nvar; v++) (*c)->Ph[v] = p[v];
    //
    // This position checking doesn't work for periodic boundaries!
    //
    //for (int v=0; v<SimPM.ndim; v++)       
    //  if (ipos[v]!=cpos[v])
    //	cout <<"*** Position x["<<v<<"] for received cell "<<i<<": got x="<<ipos[v]<<" expected x="<<cpos[v]<<"; distance="<<ipos[v]-cpos[v]<<"\n";
    //rep.printVec("\trecvd",ipos,SimPM.ndim);
    //rep.printVec("\tlocal",cpos,SimPM.ndim);
    ++c;
  }

  // 
  // Check that we got to the end of the data:
  //
  //  cout <<myrank<<"\tcomm_unpack_data: length of data in boundary = "<<b->data.size()<<"\n";
  //  cout <<myrank<<"\tcomm_unpack_data: cells received for boundary= "<<ncell<<"\n";
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: checking we got the right amount of data.\n";
#endif //TESTING
  if (c!=l->end()) {
    cerr <<"length of data = "<<l->size()<<"\n";
    rep.error("Didn't get enough cells!",CI.get_ipos((*c),XX));
  }
  if (position != ct) rep.error("length of data doesn't match",position-ct);

  //
  // Free memory and delete entry from recv_list
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: freeing memory\n";
#endif //TESTING
  if (recv_list.size()==1) {
    recv_list.pop_front();
  }
  else
    rep.error("recv list is big!",recv_list.size());
  info= mem.myfree(info);
  buf = mem.myfree(buf);
  ipos = mem.myfree(ipos);
  p = mem.myfree(p);

#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: returning\n";
#endif //TESTING
  return 0;
}




// ##################################################################
// ##################################################################




int comm_mpi::send_double_data(const int to_rank,      ///< rank to send to.
			       const long int n_el, ///< size of buffer, in number of doubles.
			       const double *data, ///< pointer to double array.
			       string &id,       ///< identifier for send, for tracking delivery later.
			       const int comm_tag      ///< comm_tag, to say what kind of send this is.
			       )
{
  if (!data) rep.error("comm_mpi::send_double_data() null pointer!",data);
  int err=0;

  //
  // Allocate memory for a record of the send, and for send_buffer
  //
  double *send_buf=0;
  struct sent_info *si=0;
  send_buf = mem.myalloc(send_buf, n_el);
  si       = mem.myalloc(si,       1   );
  //
  // It's up to the caller to ensure there are no array over-runs.
  //
  for (int i=0;i<n_el;i++) send_buf[i] = data[i];

  //
  // Add info to send record.
  //
  si->comm_tag = comm_tag;
  si->from_rank = myrank;
  si->to_rank = to_rank;
  si->data = reinterpret_cast<void *>(send_buf);
  si->type = COMM_DOUBLEDATA;
  ostringstream temp; temp.str("");
  temp <<"from "<<myrank<<" to "<<to_rank<<" tag="<<comm_tag;
  id = temp.str(); temp.str("");
  si->id = id;
  comm_mpi::sent_list.push_back(si);

  //
  // Non-blocking send of data:
  // NON-BLOCKING: int MPI_Isend(void* buf, int count, MPI_Datatype datatype, 
  //                             int dest, int tag, MPI_Comm comm, MPI_Request *request)
  //
  err += MPI_Isend(si->data, n_el, MPI_DOUBLE, si->to_rank,
                   si->comm_tag, MPI_COMM_WORLD, &(si->request));

  //
  // All done, so return.
  //
  return err;
}




// ##################################################################
// ##################################################################



int comm_mpi::receive_double_data(
        const int from_rank, ///< rank of process we are receiving from.
	const int comm_tag,  ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
        const string &id,    ///< identifier for receive, for any book-keeping.
	const long int nel,  ///< number of doubles to receive
	double *data         ///< Pointer to array to write to (must be already initialised).
	)
{
  int err=0;
#ifdef TESTING
  cout <<"comm_mpi::receive_double_data: starting.\n";
#endif //TESTING

  //
  // Find the recv in the list of active receives.
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_double_data: recv_list size="<<recv_list.size()<<"\n";
#endif //TESTING
  if (recv_list.empty()) rep.error("Call look4data before receive_data",recv_list.size());

  struct recv_info *info=0;
  list<recv_info *>::iterator iget;
  for (iget=recv_list.begin(); iget!=recv_list.end(); ++iget) {
    info = (*iget);
    if (info->id==id) break;
  }
  if (iget==recv_list.end())
    rep.error("Failed to find recv with id:",id);
#ifdef TESTING
  cout <<"found recv id="<<info->id<<" and looking for id="<<id<<"\n";
  cout <<"comm_mpi::receive_double_data: found recv id\n";
#endif //TESTING

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank) rep.error("from_ranks don't match",from_rank-info->from_rank);
  if (comm_tag  != info->comm_tag ) rep.error("Bad comm_tag",comm_tag);
  if (id        != info->id       ) rep.error("Bad id",id);
  if (info->type!= COMM_DOUBLEDATA) rep.error("data is not double array!",info->type);

  //
  // See how many doubles we are getting.
  // int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count) 
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_double_data: getting size of buffer.\n";
#endif //TESTING
  int ct = 0;
  err += MPI_Get_count(&(info->status), MPI_DOUBLE, &ct);
  if (err) rep.error("getting size of message to receive",err);
#ifdef TESTING
  cout <<"comm_mpi::receive_double_data: buffer is "<<ct<<" elements.\n.";
#endif //TESTING
  if (ct<0) rep.error("bad buffer size count",ct);
  if (ct != nel) {
    cout <<"comm_mpi::receive_double_data: getting "<<ct<<" doubles, but expecting "<<nel<<"\n";
    if (ct > nel) rep.error("Getting more data than expected, may overflow buffer!!!",ct-nel);
  }

  //
  // Now Receive the data
  // int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_double_data: receiving buffer of data from rank: "<<from_rank<<"\n";
#endif //TESTING
  void *buf = data;
  err += MPI_Recv(buf, ct, MPI_DOUBLE, info->from_rank, info->comm_tag, MPI_COMM_WORLD, &(info->status));
  if (err) rep.error("MPI_Recv failed", err);
  //  cout <<"\tcomm_receive_any_data: status: "<<status.MPI_SOURCE<<", "<<status.MPI_TAG<<", ct="<<ct<<"\n";

  //
  // Free memory and delete entry from recv_list
  //
#ifdef TESTING
  cout <<"comm_mpi::receive_cell_data: freeing memory\n";
#endif //TESTING
  if (recv_list.size()==1) {
    recv_list.pop_front();
  }
  else
    rep.error("recv list is big!",recv_list.size());
  info= mem.myfree(info);

  return 0;
}




// ##################################################################
// ##################################################################



#ifdef SILO
int comm_mpi::silo_pllel_init(const int num_files,    ///< number of files to write
			      const std::string iotype, ///< READ or WRITE
			      const std::string session_id, ///< identifier for this read/write.
			      int *group_rank,    ///< rank of group (nth file).
			      int *rank_in_group  ///< rank in group (nth proc in file i).
			      )
{
  comm_mpi::silo_id = session_id;
  comm_mpi::bat=0;
  int tag=0;
  PMPIO_iomode_t iomode=PMPIO_READ;
  if      (iotype=="READ" ) {tag = PMPIO_READ;  iomode=PMPIO_READ; }
  else if (iotype=="WRITE") {tag = PMPIO_WRITE; iomode=PMPIO_WRITE;}
  else rep.error("Bad iotype",iotype);
  
  bat = PMPIO_Init(num_files, iomode, MPI_COMM_WORLD, tag, 
                   PMPIO_DefaultCreate, PMPIO_DefaultOpen,
                   PMPIO_DefaultClose, 0);
  if (!bat) rep.error("failed to init files",myrank);

  // Get the index of the group I am in.
  *group_rank = PMPIO_GroupRank(bat,myrank);

  // Get the rank of my processor within the group:
  *rank_in_group = PMPIO_RankInGroup(bat,myrank);

  return 0;
}



// ##################################################################
// ##################################################################


int comm_mpi::silo_pllel_wait_for_file(
        const std::string id, ///< identifier for this read/write.
        const std::string filename, ///< File Name
        const std::string dir, ///< Directory to open in file.
        DBfile **dbfile ///< pointer that file gets returned in.
        )
{
#ifdef TESTING
  cout <<"comm_mpi::silo_pllel_wait_for_file() opening file: ";
  cout <<filename<<" into directory: "<<dir<<"\n";
#endif

  if (id!=silo_id) rep.error(id, silo_id);
  if (*dbfile) rep.error("please pass in null file pointer!",*dbfile);
  if (!bat) rep.error("call init before wait_for_file!",bat);


  *dbfile = static_cast<DBfile *>(PMPIO_WaitForBaton(bat, filename.c_str(), dir.c_str()));
  if (!(*dbfile)) {
    rep.error("wait for baton failed to return file pointer.",*dbfile);
  }

  return 0;
}

int comm_mpi::silo_pllel_finish_with_file(
        const std::string id, ///< identifier for this read/write.
	DBfile **dbfile ///< pointer to file we have been working on.
	)
{
#ifdef TESTING
  cout <<"comm_mpi::silo_pllel_finish_with_file() passing file on to next proc.\n";
#endif
  if (id!=silo_id) rep.error(id, silo_id);
  if (!(*dbfile)) rep.error("file pointer is null!",*dbfile);

  PMPIO_HandOffBaton(bat,static_cast<void *>(*dbfile));
  PMPIO_Finish(bat);

  bat =0;
  silo_id.erase();
  return 0;
}



// ##################################################################
// ##################################################################


void comm_mpi::silo_pllel_get_ranks(const std::string id, ///< identifier for this read/write.
				    const int proc, ///< processor rank
				    int *group, ///< rank of group processor is in.
				    int *rank   ///< rank of processor within group.
				    )
{
  if (id!=silo_id) rep.error(id, silo_id);

  *group = PMPIO_GroupRank(bat,proc);
  *rank  = PMPIO_RankInGroup(bat,proc);
  return;
}

#endif // SILO


// ##################################################################
// ##################################################################



#endif // USE_MPI
#endif // PARALLEL
