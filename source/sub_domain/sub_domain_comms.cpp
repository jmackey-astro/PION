/// \file Sub_domain.cpp
///
/// \brief Contains Sub_domain class functions for multi-process communication
/// using the Message Passing Interface (MPI).
///
/// \author Jonathan Mackey
/// \date 2009-01-23.
//
///  - 2010-02-02 JM: comments.
///
/// Modifications:\n
/// - 2010-02-03 JM: changed variable name (used 'i' twice in
///    function). Sub_domain::receive_double_data() wasn't checking
///    that I got the number of elements expected either, so I fixed
///    that.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.03.22 JM: removed old myalloc/myfree functions with string
///    args.
/// - 2011.06.02 JM: initialised two vars to values to get rid of
///    compiler warnings.
/// - 2012.05.15 JM: Added function for global-reduce (max/min/sum)
///    of arrays.
/// - 2013.01.17 JM: Made class less verbose; wrapped cout in TEST_COMMS
///    flags.
/// - 2015.01.26 JM: removed references to mpiPM (no longer global),
///    added COMM pointer setup, added get_rank_nproc()
/// - 2015.08.05 JM: Added pion_flt to send/receive for cell data.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).

#ifdef PARALLEL

#include "tools/command_line_interface.h"
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>

#include <sstream>

#include "sub_domain.h"

using namespace std;

// ##################################################################
// ##################################################################

int Sub_domain::abort()
{
  MPI_Abort(MPI_COMM_WORLD, 999);
  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::barrier(const std::string)
{
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

// ##################################################################
// ##################################################################

double Sub_domain::global_operation_double(
    const string s,  ///< Either Max, Min, or Sum
    const double d   ///< this process's max/min value.
    ) const
{
  double local = d, global = 0.0;
  int err = 0;

  if (s == "MAX") {
    global = -HUGEVALUE;
    err +=
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else if (s == "MIN") {
    global = HUGEVALUE;
    err +=
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else if (s == "SUM") {
    global = 0.0;
    err +=
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else
    spdlog::error(
        "{}: {}", "Sub_domain:global_operation_double: Bad identifier", s);

  if (err)
    spdlog::error(
        "{}: {}",
        "Sub_domain:global_operation_double: MPI call returned abnormally",
        err);

  return global;
}

// ##################################################################
// ##################################################################

void Sub_domain::global_op_double_array(
    const std::string s,  ///< MAX,MIN,SUM
    const size_t Nel,     ///< Number of elements in array.
    double *data          ///< pointer to this process's data array.
    ) const
{
  int err = 0;
  double op_data[Nel];

  if (s == "MAX") {
    // err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel,
    // MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    err += MPI_Allreduce(
        static_cast<void *>(data), static_cast<void *>(op_data), Nel,
        MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else if (s == "MIN") {
    // err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel,
    // MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    err += MPI_Allreduce(
        static_cast<void *>(data), static_cast<void *>(op_data), Nel,
        MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else if (s == "SUM") {
    // err += MPI_Allreduce(MPI_IN_PLACE, static_cast<void *>(data), Nel,
    // MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    err += MPI_Allreduce(
        static_cast<void *>(data), static_cast<void *>(op_data), Nel,
        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else
    spdlog::error(
        "{}: {}", "Sub_domain:global_op_double_array: Bad identifier", s);

  if (err)
    spdlog::error(
        "{}: {}",
        "Sub_domain:global_op_double_array: MPI call returned abnormally", err);

  for (size_t v = 0; v < Nel; v++)
    data[v] = op_data[v];
  return;
}

// ##################################################################
// ##################################################################

int Sub_domain::broadcast_data(
    const int sender,        ///< rank of sender.
    const std::string type,  ///< Type of data INT,DOUBLE,etc.
    const int n_el,          ///< number of elements
    void *data               ///< pointer to data.
    ) const
{
  int err = 0;
  if (type == "DOUBLE") {
    err += MPI_Bcast(
        static_cast<double *>(data), n_el, MPI_DOUBLE, sender, MPI_COMM_WORLD);
  }
  else if (type == "INT") {
    err += MPI_Bcast(
        static_cast<int *>(data), n_el, MPI_INT, sender, MPI_COMM_WORLD);
  }
  else if (type == "FLOAT") {
    err += MPI_Bcast(
        static_cast<float *>(data), n_el, MPI_FLOAT, sender, MPI_COMM_WORLD);
  }
  else
    spdlog::error("{}: {}", "Bad type of data to send", type);
  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::send_cell_data(
    const int to_rank,     ///< rank to send to.
    std::list<cell *> *l,  ///< list of cells to get data from.
    long int nc,           ///< number of cells in list (extra checking!)
    const int ndim,        ///< ndim
    const int nvar,        ///< nvar
    string &id,         ///< identifier for send, for tracking delivery later.
    const int comm_tag  ///< comm_tag, to say what kind of send this is.
)
{

  //
  // First initialise everything and check we have cells to get data from.
  //
  if (!id.empty()) id.erase();
  if (nc == 0 || (l->empty())) {
    spdlog::error("{}\t Nothing to send to rank: {} !!!", myrank, to_rank);
    return 1;
  }
  if (to_rank < 0 || to_rank > nproc) {
    spdlog::error(
        "nc={}  ndim={}  nvar={}  id={}  tag={}", nc, ndim, nvar, id, comm_tag);
    spdlog::error("{}: {}", "to_rank is out of bounds", to_rank);
  }

  list<cell *>::iterator c = l->begin();
  int err                  = 0;

  //
  // Determine size of send buffer needed
  //
  int unitsize = ndim * sizeof(CI.get_ipos(*c, 0)) + nvar * sizeof((*c)->P[0])
                 + sizeof((*c)->id);
  long int totalsize = 0;
  totalsize          = sizeof(long int) + nc * unitsize;
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::send_cell_data() Send buffer unitsize= {} total size = {}",
      unitsize, totalsize);
#endif

  //
  // Allocate memory for send buffer, and for the record of the send.
  //
  char *send_buff      = 0;
  struct sent_info *si = 0;
  send_buff            = mem.myalloc(send_buff, totalsize);
  si                   = mem.myalloc(si, 1);

  //
  // Pack data, using MPI.
  //
  // int MPI_Pack(void* inbuf, int incount, MPI_Datatype datatype, void
  // *outbuf, int outsize, int *position, MPI_Comm comm)
  //
  int position = 0, ct = 0, ipos[MAX_DIM];
  err += MPI_Pack(
      reinterpret_cast<void *>(&nc), 1, MPI_LONG,
      reinterpret_cast<void *>(send_buff), totalsize, &position,
      MPI_COMM_WORLD);
  do {
    CI.get_ipos(*c, ipos);
    err += MPI_Pack(
        reinterpret_cast<void *>(&((*c)->id)), 1, MPI_LONG,
        reinterpret_cast<void *>(send_buff), totalsize, &position,
        MPI_COMM_WORLD);
    err += MPI_Pack(
        reinterpret_cast<void *>(ipos), ndim, MPI_INT,
        reinterpret_cast<void *>(send_buff), totalsize, &position,
        MPI_COMM_WORLD);

#ifndef NDEBUG
    // rep.printVec("cpos",(*c)->pos,ndim);
    // rep.printVec("P ",(*c)->P,nvar);
    // rep.printVec("Ph",(*c)->Ph,nvar);
#endif

#if defined PION_DATATYPE_DOUBLE
    err += MPI_Pack(
        reinterpret_cast<void *>((*c)->Ph), nvar, MPI_DOUBLE,
        reinterpret_cast<void *>(send_buff), totalsize, &position,
        MPI_COMM_WORLD);
#elif defined PION_DATATYPE_FLOAT
    err += MPI_Pack(
        reinterpret_cast<void *>((*c)->Ph), nvar, MPI_FLOAT,
        reinterpret_cast<void *>(send_buff), totalsize, &position,
        MPI_COMM_WORLD);
#else
#error "MUST define either PION_DATATYPE_FLOAT or PION_DATATYPE_DOUBLE"
#endif

    ct++;
    ++c;  // next cell in list.
  } while (c != l->end());

  //
  // Check that we packed the right amount of data:
  //
  //  cout <<myrank<<"\tcomm_pack_send_data: bufsiz: ";
  //  cout <<totalsize<<" position in buf: "<<position<<" nc="<<nc<<"
  //  ct="<<ct<<"\n";
  if (ct != nc) {
    spdlog::error("{}: {}", "Length of list doesn't match nc", ct - nc);
  }
  if (position > totalsize) {
    spdlog::error("{}: {}", "Overwrote buffer length", position - totalsize);
  }
  if (err) {
    spdlog::error("{}: {}", "MPI_Pack returned abnormally", err);
  }

  //
  // Create record of the send.
  //
  si->comm_tag  = comm_tag;
  si->from_rank = myrank;
  si->to_rank   = to_rank;
  si->data      = reinterpret_cast<void *>(send_buff);
  si->type      = COMM_CELLDATA;
  sent_list.push_back(si);

#ifndef NDEBUG
  // rep.printVec("send_buff",reinterpret_cast<pion_flt *>(send_buff),113);
#endif

  //
  // Send data. (non-blocking, so i can return and receive a boundary)
  //
  err += MPI_Isend(
      si->data, position, MPI_PACKED, si->to_rank, si->comm_tag, MPI_COMM_WORLD,
      &(si->request));
  if (err) spdlog::error("{}: {}", "MPI_Send failed", err);

#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::send_cell_data: sending {} bytes in buffer size: {} to rank {}",
      position, totalsize, to_rank);
#endif

  ostringstream temp;
  temp.str("");
  temp << "from " << myrank << " to " << to_rank << " tag=" << comm_tag;
  id = temp.str();
  temp.str("");
  si->id = id;

  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::wait_for_send_to_finish(
    string &id  ///< identifier for the send we are waiting on.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::wait_for_send_to_finish() starting\n", myrank);
  // "<<id<<"\n";
#endif  // TEST_COMMS
  //
  // Find the send in the list of active sends.
  //
  int el = 0;
  list<sent_info *>::iterator i;
  struct sent_info *si = 0;

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::wait_for_send_to_finish()  finding in list, send id={}",
      myrank, id);
#endif  // TEST_COMMS

  for (i = sent_list.begin(); i != sent_list.end(); ++i) {
    si = (*i);
    el++;
#ifndef NDEBUG
    spdlog::debug("\t\t el={}, id = {}", el, si->id);
#endif  // TEST_COMMS
    if (si->id == id) break;
  }
  if (i == sent_list.end())
    spdlog::error("{}: {}", "Failed to find send with id:", id);

#ifndef NDEBUG
  spdlog::debug("found send id={} and looking for id={}", si->id, id);
  // cout <<"rank: "<<myrank<<"  Sub_domain::wait_for_send_to_finish() found
  // this send.\n";
#endif  // TEST_COMMS

  //
  // Use MPI to wait for the send to complete:
  //
  MPI_Status status;

#ifndef NDEBUG
  // cout <<"rank: "<<myrank<<"  Sub_domain::wait_for_send_to_finish() waiting
  // for send.\n"; cout <<"request: "<<si->request<<" and address:
  // "<<&(si->request)<<"\n";
#endif  // TEST_COMMS

  int err = MPI_Wait(&(si->request), &status);
  if (err) spdlog::error("{}: {}", "Waiting for send to complete!", err);

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::wait_for_send_to_finish() send has finished",
      myrank);
#endif  // TEST_COMMS

  //
  // Now free the data, and remove the send from the list.
  //
#ifndef NDEBUG
  // cout <<"rank: "<<myrank<<"  Sub_domain::wait_for_send_to_finish() deleting
  // send record\n";
#endif  // TEST_COMMS

  if (si->type == COMM_CELLDATA) {
    char *t = reinterpret_cast<char *>(si->data);
    t       = mem.myfree(t);
  }
  else if (si->type == COMM_DOUBLEDATA) {
    double *t = reinterpret_cast<double *>(si->data);
    t         = mem.myfree(t);
  }
  else
    spdlog::error("{}: {}", "BAD type of data!", si->type);
  si = mem.myfree(si);
  sent_list.erase(i);

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::wait_for_send_to_finish() returning", myrank);
#endif  // TEST_COMMS
  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::look_for_data_to_receive(
    int *from_rank,      ///< rank of sender
    string &id,          ///< identifier for receive.
    int *recv_tag,       ///< comm_tag associated with data.
    const int comm_tag,  ///< comm_tag: (PER,MPI,F2C,C2F)
    const int type       ///< type of data we are looking for.
)
{

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::look_for_data_to_receive() starting\n", myrank);
#endif  // TEST_COMMS
  //
  // Create a new received info record.
  //
  struct recv_info *ri = 0;
  ri                   = mem.myalloc(ri, 1);
  ri->to_rank          = myrank;
  ri->data             = 0;
  if (type != COMM_CELLDATA && type != COMM_DOUBLEDATA)
    spdlog::error("{}: {}", "only know two types of data to look for!", type);
  ri->type = type;

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::look_for_data_to_receive() looking for source for comm_tag: {}",
      myrank, comm_tag);

  // list<sent_info *>::iterator i;
  // struct sent_info *si=0;
  // int el=0;
  // for (i=sent_list.begin(); i!=sent_list.end(); ++i) {
  //  si = (*i); el++;
  //  cout <<"SENT LIST: el="<<el<<", id = "<<si->id<<"\n";
  //  cout.flush();
  //}
#endif  // TEST_COMMS
  //
  // First Look for Data, and wait till we find some.
  // Find a message that has been sent to us.
  // int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status)
  //
  // int err = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
  // &(ri->status)); cout <<" max_tag = "<<MPI_TAG_UB<<"<<
  // max_int="<<INT_MAX<<endl;

  // from_rank is set to <0 if we want any source, otherwise we want
  // data from a specific rank.
  int err = 0;
  if (*from_rank >= 0) {
    err = MPI_Probe(*from_rank, comm_tag, MPI_COMM_WORLD, &(ri->status));
  }
  else {
    err = MPI_Probe(MPI_ANY_SOURCE, comm_tag, MPI_COMM_WORLD, &(ri->status));
  }
  if (err) spdlog::error("{}: {}", "mpi probe failed", err);
  if (*from_rank >= 0 && *from_rank != ri->status.MPI_SOURCE) {
    spdlog::error(
        "{}: {}", "looking for specific rank but got different",
        ri->status.MPI_SOURCE);
  }

#ifndef NDEBUG
  spdlog::debug(
      "rank: {}  Sub_domain::look_for_data_to_receive() found data from rank {}",
      myrank, *from_rank);
#endif  // TEST_COMMS

  //
  // Now get the size and source of the data.
  //
  *from_rank = ri->status.MPI_SOURCE;
  *recv_tag  = ri->status.MPI_TAG;
  if (*recv_tag != comm_tag)
    spdlog::error("{}: {}", "bad comm_tag received", *recv_tag);
  ostringstream temp;
  temp.str("");
  temp << "from " << *from_rank << " to " << myrank << " tag=" << *recv_tag;
  id = temp.str();
  temp.str("");

#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::look_for_data_to_receive: found data from {} with tag:{}",
      *from_rank, *recv_tag);
#endif  // TEST_COMMS

  //
  // Set record of where data is coming from.
  ri->id        = id;
  ri->comm_tag  = *recv_tag;
  ri->from_rank = *from_rank;
  recv_list.push_back(ri);

#ifndef NDEBUG
  spdlog::debug("Sub_domain::look_for_data_to_receive: returning");
#endif  // TEST_COMMS

  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::receive_cell_data(
    const int from_rank,   ///< rank of process we are receiving from.
    std::list<cell *> *l,  ///< list of cells to get data for.
    const long int nc,     ///< number of cells in list (extra checking!)
    const int ndim,        ///< ndim
    const int nvar,        ///< nvar
    const int comm_tag,    ///< comm_tag: (PER,MPI,F2C,C2F)
    const string &id       ///< identifier for receive, for book-keeping.
)
{
  int err = 0;

#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: starting");
#endif  // TEST_COMMS

  //
  // Find the recv in the list of active receives.
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: recv_list size={}", recv_list.size());
#endif  // TEST_COMMS

  if (recv_list.empty())
    spdlog::error(
        "{}: {}", "Call look4data before receive_data", recv_list.size());

  struct recv_info *info = 0;
  list<recv_info *>::iterator iget;
  for (iget = recv_list.begin(); iget != recv_list.end(); ++iget) {
    info = (*iget);
    if (info->id == id) break;
  }
  if (iget == recv_list.end())
    spdlog::error("{}: {}", "Failed to find recv with id:", id);

#ifndef NDEBUG
  spdlog::debug(
      "found recv id={} and looking for id={}\nSub_domain::receive_cell_data: found recv id\n",
      info->id, id);
#endif  // TEST_COMMS

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank)
    spdlog::error(
        "{}: {}", "from_ranks don't match", from_rank - info->from_rank);
  if (comm_tag != info->comm_tag)
    spdlog::error("{}: {}", "Bad comm_tag", comm_tag);
  if (id != info->id) spdlog::error("{}: {}", "Bad id", id);

    //
    // See how many bytes we are getting!
    // int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int
    // *count)
    //
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: getting size of buffer");
#endif  // TEST_COMMS

  int ct = 0;
  err += MPI_Get_count(&(info->status), MPI_PACKED, &ct);
  if (err) spdlog::error("{}: {}", "getting size of message to receive", err);

#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: buffer is {} bytes", ct);
#endif  // TEST_COMMS

  if (ct < 0) spdlog::error("{}: {}", "bad buffer size count", ct);

  //
  // Allocate memory for data.
  //
  char *buf = 0;
  buf       = mem.myalloc(buf, ct);

  //
  // Now Receive the data in recv_buff
  // int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int
  // tag, MPI_Comm comm, MPI_Status *status)
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: receiving buffer of data from rank: {}",
      from_rank);
#endif  // TEST_COMMS

  err += MPI_Recv(
      buf, ct, MPI_PACKED, from_rank, comm_tag, MPI_COMM_WORLD,
      &(info->status));
  if (err) spdlog::error("{}: {}", "MPI_Recv failed", err);
#ifndef NDEBUG
  spdlog::debug(
      "\tcomm_receive_any_data: status: {}, {}, ct={}", info->status.MPI_SOURCE,
      info->status.MPI_TAG, ct);
#endif  // TEST_COMMS

  //
  // First unpack the cells received count.
  //
  //   int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf,
  //   int outcount, MPI_Datatype datatype, MPI_Comm comm)

#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: getting number of cells received");
#endif  // TEST_COMMS

  int position   = 0;
  long int ncell = 0;

  err = MPI_Unpack(buf, ct, &position, &ncell, 1, MPI_LONG, MPI_COMM_WORLD);
  if (err) spdlog::error("{}: {}", "Unpack", err);
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: got {} cells", ncell);
#endif  // TEST_COMMS

  if (ncell != static_cast<int>(l->size()) || (ncell != nc)) {
    spdlog::error(
        "{}\tSub_domain:recv: length of data = {} or we're told it's {}\n{}\tSub_domain:recv: cells received = {}\n{}\tSub_domain:recv: Bugging out!\n",
        myrank, l->size(), nc, myrank, ncell, myrank);
    return (-99);
  }

#ifndef NDEBUG
  spdlog::debug(
      "{}\tSub_domain:recv: got {} cells, as expected", myrank, ncell);
#endif  // TEST_COMMS

  //
  // Allocate arrays to unpack data into.
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: allocating arrays to unpack data into");
#endif  // TEST_COMMS

  int *ipos   = 0;
  pion_flt *p = 0;
  ipos        = mem.myalloc(ipos, ndim);
  p           = mem.myalloc(p, nvar);

  //
  // Data should come in in the order the boundary cells are listed in.
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: unpacking data into list of cells");
#endif  // TEST_COMMS

  int cpos[MAX_DIM];
  long int c_id            = 0;
  list<cell *>::iterator c = l->begin();

  for (int i = 0; i < ncell; i++) {
#ifndef NDEBUG
    // cout <<"Sub_domain::receive_cell_data: unpacking cell "<<i<<".\n";
    // cout.flush();
#endif  // TEST_COMMS
    CI.get_ipos(*c, cpos);
#ifndef NDEBUG
    // cout <<"position="<<position<<" : cell id = "<<(*c)->id<<" : ";
    // rep.printVec("cpos",cpos,ndim);
    // cout.flush();
#endif  // TEST_COMMS
    err += MPI_Unpack(buf, ct, &position, &c_id, 1, MPI_LONG, MPI_COMM_WORLD);
#ifndef NDEBUG
    // cout <<"position="<<position<<" : ";
    // cout <<"cell id = "<<c_id<<"\n";
    // cout.flush();
#endif  // TEST_COMMS
    err += MPI_Unpack(buf, ct, &position, ipos, ndim, MPI_INT, MPI_COMM_WORLD);
#ifndef NDEBUG
    // cout <<"position="<<position<<" : ";
    // rep.printVec("recvd pos",ipos,ndim);
    // cout.flush();
#endif  // TEST_COMMS
#if defined PION_DATATYPE_DOUBLE
    err += MPI_Unpack(buf, ct, &position, p, nvar, MPI_DOUBLE, MPI_COMM_WORLD);
#elif defined PION_DATATYPE_FLOAT
    err += MPI_Unpack(buf, ct, &position, p, nvar, MPI_FLOAT, MPI_COMM_WORLD);
#else
#error "MUST define either PION_DATATYPE_FLOAT or PION_DATATYPE_DOUBLE"
#endif
#ifndef NDEBUG
    // rep.printVec("\tlocal",cpos,ndim);
    // cout <<"position="<<position<<" : ";
    // rep.printVec("data var",p,nvar);
    // cout.flush();
#endif  // TEST_COMMS
    if (err) spdlog::error("{}: {}", "Unpack", err);

    // For a given boundary data iterator, put data into cells
    if (c == l->end()) spdlog::error("{}: {}", "Got too many cells!", i);
    for (int v = 0; v < nvar; v++)
      (*c)->Ph[v] = p[v];
      //
      // This position checking doesn't work for periodic boundaries
      //
      // rep.printVec("ipos",ipos,ndim);
      // rep.printVec("cpos",cpos,ndim);
      // for (int v=0; v<ndim; v++)
      //  if (ipos[v]!=cpos[v]) {
      //    cout <<"*** Position x["<<v<<"] for received cell "<<i;
      //    cout <<": got x="<<ipos[v]<<" expected x="<<cpos[v];
      //    cout <<"; distance="<<ipos[v]-cpos[v]<<"\n";
      //    spdlog::error("{}: {}", "positions do not match",ipos[v]);
      //  }
#ifndef NDEBUG
      // rep.printVec("\trecvd",ipos,ndim);
      // rep.printVec("\tlocal",cpos,ndim);
      // cout.flush();
#endif
    ++c;
  }

  //
  // Check that we got to the end of the data:
  //
  //  cout <<myrank<<"\tcomm_unpack_data: length of data in boundary =
  //  "<<b->data.size()<<"\n"; cout <<myrank<<"\tcomm_unpack_data: cells
  //  received for boundary= "<<ncell<<"\n";
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: checking we got the right amount of data");
#endif  // TEST_COMMS
  if (c != l->end()) {
    spdlog::error("length of data = {}", l->size());
    spdlog::error("{}: {}", "Didn't get enough cells!", CI.get_ipos((*c), XX));
  }
  if (position != ct)
    spdlog::error("{}: {}", "length of data doesn't match", position - ct);

    //
    // Free memory and delete entry from recv_list
    //
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: freeing memory");
#endif  // TEST_COMMS
  if (recv_list.size() == 1) {
    recv_list.pop_front();
  }
  else
    spdlog::error("{}: {}", "recv list is big!", recv_list.size());
  info = mem.myfree(info);
  buf  = mem.myfree(buf);
  ipos = mem.myfree(ipos);
  p    = mem.myfree(p);

#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: returning");
#endif  // TEST_COMMS
  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::send_double_data(
    const int to_rank,    ///< rank to send to.
    const long int n_el,  ///< size of buffer, in number of doubles.
    const double *data,   ///< pointer to double array.
    string &id,           ///< identifier for send, for tracking delivery.
    const int comm_tag    ///< comm_tag, to say what kind of send.
)
{
  if (!data)
    spdlog::error(
        "{}: {}", "Sub_domain::send_double_data() null pointer!",
        fmt::ptr(data));
  int err = 0;

  //
  // Allocate memory for a record of the send, and for send_buffer
  //
  double *send_buf     = 0;
  struct sent_info *si = 0;
  send_buf             = mem.myalloc(send_buf, n_el);
  si                   = mem.myalloc(si, 1);
  //
  // It's up to the caller to ensure there are no array over-runs.
  //
  for (int i = 0; i < n_el; i++)
    send_buf[i] = data[i];

  //
  // Add info to send record.
  //
  si->comm_tag  = comm_tag;
  si->from_rank = myrank;
  si->to_rank   = to_rank;
  si->data      = reinterpret_cast<void *>(send_buf);
  si->type      = COMM_DOUBLEDATA;
  ostringstream temp;
  temp.str("");
  temp << "from " << myrank << " to " << to_rank << " tag=" << comm_tag;
  id = temp.str();
  temp.str("");
  si->id = id;
  sent_list.push_back(si);

  //
  // Non-blocking send of data:
  // NON-BLOCKING: int MPI_Isend(void* buf, int count, MPI_Datatype datatype,
  //                             int dest, int tag, MPI_Comm comm, MPI_Request
  //                             *request)
  //
  err += MPI_Isend(
      si->data, n_el, MPI_DOUBLE, si->to_rank, si->comm_tag, MPI_COMM_WORLD,
      &(si->request));

  //
  // All done, so return.
  //
  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::receive_double_data(
    const int from_rank,  ///< rank of process we are receiving from.
    const int comm_tag,   ///< comm_tag: what sort of comm we are looking for
                          ///< (PER,MPI,etc.)
    const string &id,     ///< identifier for receive, for any book-keeping.
    const long int nel,   ///< number of doubles to receive
    double
        *data  ///< Pointer to array to write to (must be already initialised).
)
{
  int err = 0;
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_double_data: starting");
#endif  // TEST_COMMS

  //
  // Find the recv in the list of active receives.
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_double_data: recv_list size={}", recv_list.size());
#endif  // TEST_COMMS
  if (recv_list.empty())
    spdlog::error(
        "{}: {}", "Call look4data before receive_data", recv_list.size());

  struct recv_info *info = 0;
  list<recv_info *>::iterator iget;
  for (iget = recv_list.begin(); iget != recv_list.end(); ++iget) {
    info = (*iget);
    if (info->id == id) break;
  }
  if (iget == recv_list.end())
    spdlog::error("{}: {}", "Failed to find recv with id:", id);
#ifndef NDEBUG
  spdlog::debug(
      "found recv id={} and looking for id={}\nSub_domain::receive_double_data: found recv id",
      info->id, id);
#endif  // TEST_COMMS

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank)
    spdlog::error(
        "{}: {}", "from_ranks don't match", from_rank - info->from_rank);
  if (comm_tag != info->comm_tag)
    spdlog::error("{}: {}", "Bad comm_tag", comm_tag);
  if (id != info->id) spdlog::error("{}: {}", "Bad id", id);
  if (info->type != COMM_DOUBLEDATA)
    spdlog::error("{}: {}", "data is not double array!", info->type);

    //
    // See how many doubles we are getting.
    // int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int
    // *count)
    //
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_double_data: getting size of buffer");
#endif  // TEST_COMMS
  int ct = 0;
  err += MPI_Get_count(&(info->status), MPI_DOUBLE, &ct);
  if (err) spdlog::error("{}: {}", "getting size of message to receive", err);
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_double_data: buffer is {} elements", ct);
#endif  // TEST_COMMS
  if (ct < 0) spdlog::error("{}: {}", "bad buffer size count", ct);
  if (ct != nel) {
    spdlog::debug(
        "Sub_domain::receive_double_data: getting {} doubles, but expecting {}",
        ct, nel);
    spdlog::error("{}: {}", "Getting wrong amount of data", ct - nel);
  }

  //
  // Now Receive the data
  // int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int
  // tag, MPI_Comm comm, MPI_Status *status)
  //
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_double_data: receiving buffer of data from rank: {}",
      from_rank);
#endif  // TEST_COMMS
  void *buf = data;
  err += MPI_Recv(
      buf, ct, MPI_DOUBLE, info->from_rank, info->comm_tag, MPI_COMM_WORLD,
      &(info->status));
  if (err) spdlog::error("{}: {}", "MPI_Recv failed", err);
    //  cout <<"\tcomm_receive_any_data: status: "<<status.MPI_SOURCE<<",
    //  "<<status.MPI_TAG<<", ct="<<ct<<"\n";

    //
    // Free memory and delete entry from recv_list
    //
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_double_data: freeing memory");
#endif  // TEST_COMMS
  if (recv_list.size() == 1) {
    recv_list.pop_front();
  }
  else
    spdlog::error("{}: {}", "recv list is big!", recv_list.size());
  info = mem.myfree(info);

  return 0;
}

// ##################################################################
// ##################################################################

#ifdef SILO
int Sub_domain::silo_pllel_init(
    const int num_files,           ///< number of files to write
    const std::string iotype,      ///< READ or WRITE
    const std::string session_id,  ///< identifier for this read/write.
    int *group_rank,               ///< rank of group (nth file).
    int *rank_in_group             ///< rank in group (nth proc in file i).
)
{
  silo_id               = session_id;
  bat                   = 0;
  int tag               = 0;
  PMPIO_iomode_t iomode = PMPIO_READ;
  if (iotype == "READ") {
    tag    = PMPIO_READ;
    iomode = PMPIO_READ;
  }
  else if (iotype == "WRITE") {
    tag    = PMPIO_WRITE;
    iomode = PMPIO_WRITE;
  }
  else
    spdlog::error("{}: {}", "Bad iotype", iotype);

  bat = PMPIO_Init(
      num_files, iomode, MPI_COMM_WORLD, tag, PMPIO_DefaultCreate,
      PMPIO_DefaultOpen, PMPIO_DefaultClose, 0);
  if (!bat) spdlog::error("{}: {}", "failed to init files", myrank);

  // Get the index of the group I am in.
  *group_rank = PMPIO_GroupRank(bat, myrank);

  // Get the rank of my processor within the group:
  *rank_in_group = PMPIO_RankInGroup(bat, myrank);

  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::silo_pllel_wait_for_file(
    const std::string id,        ///< identifier for this read/write.
    const std::string filename,  ///< File Name
    const std::string dir,       ///< Directory to open in file.
    DBfile **dbfile              ///< pointer that file gets returned in.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::silo_pllel_wait_for_file() opening file: {} into directory: {}",
      filename, dir);
#endif

  if (id != silo_id) spdlog::error("{}: {}", id, silo_id);
  if (*dbfile)
    spdlog::error(
        "{}: {}", "please pass in null file pointer!", fmt::ptr(*dbfile));
  if (!bat)
    spdlog::error("{}: {}", "call init before wait_for_file!", fmt::ptr(bat));

  *dbfile = static_cast<DBfile *>(
      PMPIO_WaitForBaton(bat, filename.c_str(), dir.c_str()));
  if (!(*dbfile)) {
    spdlog::error(
        "{}: {}", "wait for baton failed to return file pointer.",
        fmt::ptr(*dbfile));
  }

  return 0;
}

int Sub_domain::silo_pllel_finish_with_file(
    const std::string id,  ///< identifier for this read/write.
    DBfile **dbfile        ///< pointer to file we have been working on.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::silo_pllel_finish_with_file() passing file on to next proc");
#endif

  if (id != silo_id) spdlog::error("{}: {}", id, silo_id);
  if (!(*dbfile))
    spdlog::error("{}: {}", "file pointer is null!", fmt::ptr(*dbfile));

  PMPIO_HandOffBaton(bat, static_cast<void *>(*dbfile));
  PMPIO_Finish(bat);

  bat = 0;
  silo_id.erase();
  return 0;
}

// ##################################################################
// ##################################################################

void Sub_domain::silo_pllel_get_ranks(
    const std::string id,  ///< identifier for this read/write.
    const int proc,        ///< processor rank
    int *group,            ///< rank of group processor is in.
    int *rank              ///< rank of processor within group.
    ) const
{
  if (id != silo_id) spdlog::error("{}: {}", id, silo_id);

  *group = PMPIO_GroupRank(bat, proc);
  *rank  = PMPIO_RankInGroup(bat, proc);
}

#endif  // SILO

// ##################################################################
// ##################################################################

#endif  // PARALLEL
