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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_params.h"
#include "sub_domain.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include <algorithm>
#include <sstream>

using namespace std;

// ##################################################################
// ##################################################################

int Sub_domain::abort() const
{
  return MPI_Abort(cart_comm, 999);
}

// ##################################################################
// ##################################################################

int Sub_domain::barrier() const
{
  return MPI_Barrier(cart_comm);
}

// ##################################################################
// ##################################################################

double Sub_domain::global_operation_double(
    const mpi_op op,  ///< Either Max, Min, or Sum
    const double d    ///< this process's max/min value.
    ) const
{
  double global;
  int err = 0;

  switch (op) {
    case MAX:
      err += MPI_Allreduce(&d, &global, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
      break;
    case MIN:
      err += MPI_Allreduce(&d, &global, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
      break;
    case SUM:
      err += MPI_Allreduce(&d, &global, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
      break;
    default:
      spdlog::error(
          "Sub_domain:global_operation_double: Bad identifier: {} (0 = MAX, 1 = MIN, 2 = SUM)",
          op);
      exit(1);
  }

  if (err) {
    spdlog::error(
        "{}: {}",
        "Sub_domain:global_operation_double: MPI call returned abnormally",
        err);
    exit(1);
  }

  return global;
}

// ##################################################################
// ##################################################################

/* TODO: unused */
void Sub_domain::global_op_double_array(
    const mpi_op op,            ///< MAX,MIN,SUM
    const size_t num_elements,  ///< Number of elements in array.
    double *data                ///< pointer to this process's data array.
    ) const
{
  int err = 0;
  double op_data[num_elements];

  switch (op) {
    case MAX:
      err += MPI_Allreduce(
          static_cast<void *>(data), static_cast<void *>(op_data), num_elements,
          MPI_DOUBLE, MPI_MAX, cart_comm);
      break;
    case MIN:
      err += MPI_Allreduce(
          static_cast<void *>(data), static_cast<void *>(op_data), num_elements,
          MPI_DOUBLE, MPI_MIN, cart_comm);
      break;
    case SUM:
      err += MPI_Allreduce(
          static_cast<void *>(data), static_cast<void *>(op_data), num_elements,
          MPI_DOUBLE, MPI_SUM, cart_comm);
      break;
    default:
      spdlog::error(
          "Sub_domain:global_op_double_array: Bad identifier: {} (0 = MAX, 1 = MIN, 2 = SUM)",
          op);
      exit(1);
  }

  if (err) {
    spdlog::error(
        "{}: {}",
        "Sub_domain:global_op_double_array: MPI call returned abnormally", err);
    exit(1);
  }

  for (size_t v = 0; v < num_elements; v++)
    data[v] = op_data[v];
}

// ##################################################################
// ##################################################################

int Sub_domain::broadcast_data(
    const int sender,        ///< rank of sender.
    const mpi_type type,     ///< Type of data INT,DOUBLE,etc.
    const int num_elements,  ///< number of elements
    void *data               ///< pointer to data.
    ) const
{
  int err = 0;
  switch (type) {
    case INT:
      err += MPI_Bcast(
          static_cast<double *>(data), num_elements, MPI_INT, sender,
          cart_comm);
      break;
    case FLOAT:
      err += MPI_Bcast(
          static_cast<double *>(data), num_elements, MPI_FLOAT, sender,
          cart_comm);
      break;
    case DOUBLE:
      err += MPI_Bcast(
          static_cast<double *>(data), num_elements, MPI_DOUBLE, sender,
          cart_comm);
      break;
    default:
      spdlog::error(
          "Bad type of data to send: {} (0 = INT, 1 = FLOAT, 2 = DOUBLE)",
          type);
      exit(1);
  }
  if (err) {
    spdlog::error(
        "{}: {}", "Sub_domain:broadcast_data: MPI call returned abnormally",
        err);
    exit(1);
  }
  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::send_cell_data(
    const int to_rank,             ///< rank to send to.
    std::list<cell *> &cell_list,  ///< list of cells to get data from.
    const int nvar,                ///< nvar
    string &send_id,    ///< identifier for send, for tracking delivery later.
    const int comm_tag  ///< comm_tag, to say what kind of send this is.
)
{
  /* First initialise everything and check we have cells to get data from */
  if (cell_list.empty()) {
    spdlog::error("{}\t Nothing to send to rank: {} !!!", myrank, to_rank);
    exit(1);
    return 1;
  }

  m_num_send_cells = cell_list.size();
  if (to_rank < 0 || to_rank > nproc) {
    spdlog::error(
        "num_cells={}  nvar={}  id={}  tag={}", m_num_send_cells, nvar, send_id,
        comm_tag);
    spdlog::error("{}: {}", "to_rank is out of bounds", to_rank);
    exit(1);
  }

  int err = 0;

  /* Determine size of send buffer needed */
  const int unitsize = nvar * sizeof((*cell_list.begin())->Ph[0])
                       + sizeof((*cell_list.begin())->id);
  const long int totalsize = sizeof(long int) + m_num_send_cells * unitsize;
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::send_cell_data() Send buffer unitsize= {} total size = {}",
      unitsize, totalsize);
#endif

  /* Create record of the send */
  send_list.emplace_back();
  send_list.back().comm_tag  = comm_tag;
  send_list.back().from_rank = myrank;
  send_list.back().to_rank   = to_rank;
  send_list.back().type      = COMM_CELLDATA;
  /* Allocate memory for send buffer, and for the record of the send */
  send_list.back().send_buff.resize(totalsize);

  /*
   *  Pack data, using MPI.
   *
   * int MPI_Pack(void* inbuf, int incount, MPI_Datatype datatype, void
   * *outbuf, int outsize, int *position, MPI_Comm comm)
   *
   */
  int position = 0;
  err += MPI_Pack(
      &m_num_send_cells, 1, MPI_LONG, send_list.back().send_buff.data(),
      totalsize, &position, cart_comm);
  for (auto const &c : cell_list) {
#if defined PION_DATATYPE_DOUBLE
    err += MPI_Pack(
        c->Ph, nvar, MPI_DOUBLE, send_list.back().send_buff.data(), totalsize,
        &position, cart_comm);
#elif defined PION_DATATYPE_FLOAT
    err += MPI_Pack(
        c->Ph, nvar, MPI_FLOAT, send_list.back().send_buff.data(), totalsize,
        &position, cart_comm);
#else
#error "MUST define either PION_DATATYPE_FLOAT or PION_DATATYPE_DOUBLE"
#endif
  }

  if (position > totalsize) {
    spdlog::error("{}: {}", "Overwrote buffer length", position - totalsize);
    exit(1);
  }
  if (err) {
    spdlog::error("{}: {}", "MPI_Pack returned abnormally", err);
    exit(1);
  }

#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::send_cell_data: sending {} bytes in buffer size: {} to rank {}",
      position, totalsize, to_rank);
#endif

  /* create send tag */
  ostringstream temp;
  temp.str("");
  temp << "from " << myrank << " to " << to_rank << " tag=" << comm_tag;
  send_id = temp.str();
  temp.str("");
  send_list.back().id = send_id;

  /* Send data. (non-blocking, so i can return and receive a boundary) */
  err += MPI_Isend(
      send_list.back().send_buff.data(), position, MPI_PACKED,
      send_list.back().to_rank, send_list.back().comm_tag, cart_comm,
      &(send_list.back().request));
  if (err) {
    spdlog::error("{}: {}", "MPI_Send failed", err);
    exit(1);
  }

  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::wait_for_send_to_finish(
    string &send_id  ///< identifier for the send we are waiting on.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::wait_for_send_to_finish() finding in list, send id={}",
      send_id);
#endif  // TEST_COMMS

  /* Find the send in the list of active sends */
  auto si =
      find_if(send_list.begin(), send_list.end(), [send_id](const auto &si) {
        return si.id == send_id;
      });

  if (si == send_list.end()) {
    spdlog::error("{}: {}", "Failed to find send with id:", send_id);
    exit(1);
  }

#ifndef NDEBUG
  spdlog::debug("found send id={} while looking for id={}", si->id, send_id);
#endif  // TEST_COMMS

  spdlog::debug(
      "waiting for send to {}, with tag {} to finish", si->to_rank,
      si->comm_tag);
  MPI_Status status;
  int err = MPI_Wait(&(si->request), &status);
  if (err) {
    spdlog::error("{}: {}", "Waiting for send to complete!", err);
    exit(1);
  }

  send_list.erase(si);

#ifndef NDEBUG
  spdlog::debug("Sub_domain::wait_for_send_to_finish() returning");
#endif  // TEST_COMMS
  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::look_for_data_to_receive(
    int *from_rank,      ///< rank of sender
    string &recv_id,     ///< identifier for receive.
    int *recv_tag,       ///< comm_tag associated with data.
    const int comm_tag,  ///< comm_tag: (PER,MPI,F2C,C2F)
    const int type       ///< type of data we are looking for.
)
{
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::look_for_data_to_receive() looking for source for comm_tag: {}",
      comm_tag);
#endif  // TEST_COMMS

  /* Create a new received info record */
  recv_list.emplace_back();

  recv_list.back().to_rank = myrank;
  recv_list.back().data    = 0;
  if (type != COMM_CELLDATA && type != COMM_DOUBLEDATA) {
    spdlog::error("{}: {}", "only know two types of data to look for!", type);
    exit(1);
  }
  recv_list.back().type = type;

  /*
   * Find a message that has been sent to us.
   * int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status)
   *
   * from_rank is set to <0 if we want any source, otherwise we want
   * data from a specific rank.
   */
  int err = 0;
  if (*from_rank >= 0) {
    spdlog::debug(
        "looking to receive from rank {} with tag {}", *from_rank, comm_tag);
    err =
        MPI_Probe(*from_rank, comm_tag, cart_comm, &(recv_list.back().status));
  }
  else {
#ifndef NDEBUG
    spdlog::debug(
        "Sub_domain::look_for_data_to_receive() calling MPI_Probe, any src: {}",
        comm_tag);
#endif  // TEST_COMMS

    err = MPI_Probe(
        MPI_ANY_SOURCE, comm_tag, cart_comm, &(recv_list.back().status));
  }
  if (err) {
    spdlog::error("{}: {}", "mpi probe failed", err);
    exit(1);
  }

  /* Now get the size and source of the data */
  *from_rank = recv_list.back().status.MPI_SOURCE;
  *recv_tag  = recv_list.back().status.MPI_TAG;

  /* create communication id */
  ostringstream temp;
  temp.str("");
  temp << "from " << *from_rank << " to " << myrank << " tag=" << *recv_tag;
  recv_id = temp.str();
  temp.str("");

  /* Set record of where data is coming from */
  recv_list.back().id        = recv_id;
  recv_list.back().comm_tag  = *recv_tag;
  recv_list.back().from_rank = *from_rank;

#ifndef NDEBUG
  spdlog::debug("Sub_domain::look_for_data_to_receive: returning");
#endif  // TEST_COMMS

  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::receive_cell_data(
    const int from_rank,           ///< rank of process we are receiving from.
    std::list<cell *> &cell_list,  ///< list of cells to get data for.
    const int nvar,                ///< nvar
    const int comm_tag,            ///< comm_tag: (PER,MPI,F2C,C2F)
    const string &recv_id  ///< identifier for receive, for book-keeping.
)
{
  int err = 0;

#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: starting");
#endif  // TEST_COMMS

  auto info =
      find_if(recv_list.begin(), recv_list.end(), [recv_id](const auto &ri) {
        return ri.id == recv_id;
      });

  if (info == recv_list.end()) {
    spdlog::error("{}: {}", "Failed to find recv with id:", recv_id);
    exit(1);
  }

  /* Check that everthing matches */
  if (from_rank != info->from_rank) {
    spdlog::error(
        "{}: {}", "from_ranks don't match", from_rank - info->from_rank);
    exit(1);
  }
  if (comm_tag != info->comm_tag) {
    spdlog::error("{}: {}", "Bad comm_tag", comm_tag);
    exit(1);
  }
  if (recv_id != info->id) {
    spdlog::error("{}: {}", "Bad id", recv_id);
    exit(1);
  }

  /*
   * See how many bytes we are getting!
   * int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int
   * *count)
   */
  int buffer_size;
  err += MPI_Get_count(&(info->status), MPI_PACKED, &buffer_size);
  if (err) {
    spdlog::error("getting size of message to receive: {}", err);
    exit(1);
  }

  if (buffer_size < 0) {
    spdlog::error("bad buffer size count: {}", buffer_size);
    exit(1);
  }

#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: buffer is {} bytes", buffer_size);
#endif  // TEST_COMMS

  vector<char> buf(buffer_size);

  /*
   * Now Receive the data in recv_buff
   * int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int
   * tag, MPI_Comm comm, MPI_Status *status)
   */
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: receiving buffer of data from rank: {}",
      from_rank);
#endif  // TEST_COMMS

  err += MPI_Recv(
      buf.data(), buffer_size, MPI_PACKED, from_rank, comm_tag, cart_comm,
      &(info->status));
  if (err) {
    spdlog::error("{}: {}", "MPI_Recv failed", err);
    exit(1);
  }
#ifndef NDEBUG
  spdlog::debug(
      "\tcomm_receive_any_data: status: {}, {}, buffer_size={}",
      info->status.MPI_SOURCE, info->status.MPI_TAG, buffer_size);
#endif  // TEST_COMMS

  /*
   * First unpack the number of sent cells
   *
   *   int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf,
   *        int outcount, MPI_Datatype datatype, MPI_Comm comm)
   */

  int position = 0;
  long int ncell;

  err = MPI_Unpack(
      buf.data(), buffer_size, &position, &ncell, 1, MPI_LONG, cart_comm);
  if (err) {
    spdlog::error("{}: {}", "Unpack", err);
    exit(1);
  }

  const int num_cells = cell_list.size();
  if (ncell != num_cells) {
    spdlog::error(
        "{}\tSub_domain:recv: length of data = {} or we're told it's {}\n{}\tSub_domain:recv: cells received = {}\n{}\tSub_domain:recv: Bugging out!\n",
        myrank, cell_list.size(), num_cells, myrank, ncell, myrank);
    exit(1);
  }

  /* Data should come in in the order the boundary cells are listed in */
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_cell_data: unpacking data into list of cells");
#endif  // TEST_COMMS

  for (auto const &c : cell_list) {
#if defined PION_DATATYPE_DOUBLE
    err += MPI_Unpack(
        buf.data(), buffer_size, &position, c->Ph, nvar, MPI_DOUBLE, cart_comm);
#elif defined PION_DATATYPE_FLOAT
    err += MPI_Unpack(
        buf.data(), buffer_size, &position, c->Ph, nvar, MPI_FLOAT, cart_comm);
#else
#error "MUST define either PION_DATATYPE_FLOAT or PION_DATATYPE_DOUBLE"
#endif

    if (err) {
      spdlog::error("{}: {}", "Unpack", err);
      exit(1);
    }
  }

  if (position != buffer_size) {
    spdlog::error(
        "{}: {}", "length of data doesn't match", position - buffer_size);
    exit(1);
  }

  /* delete entry from recv_list */
  if (recv_list.size() == 1) {
    recv_list.pop_front();
  }
  else {
    spdlog::error("{}: {}", "recv list is big!", recv_list.size());
    exit(1);
  }

#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_cell_data: returning");
#endif  // TEST_COMMS
  return 0;
}

// ##################################################################
// ##################################################################

int Sub_domain::send_double_data(
    const int to_rank,            ///< rank to send to.
    const long int num_elements,  ///< size of buffer, in number of doubles.
    vector<double> &data,         ///< pointer to double array.
    string &send_id,    ///< identifier for send, for tracking delivery.
    const int comm_tag  ///< comm_tag, to say what kind of send.
)
{
  if (data.empty()) {
    spdlog::error("Sub_domain::send_double_data() empty send data!");
    exit(1);
  }
  int err = 0;

  send_list.emplace_back();
  send_list.back().comm_tag  = comm_tag;
  send_list.back().from_rank = myrank;
  send_list.back().to_rank   = to_rank;
  send_list.back().type      = COMM_DOUBLEDATA;
  /* store vector in member variable so it is not freed out of scope */
  send_list.back().send_buff_double = std::move(data);

  /* create send tag */
  ostringstream temp;
  temp.str("");
  temp << "from " << myrank << " to " << to_rank << " tag=" << comm_tag;
  send_id = temp.str();
  temp.str("");
  send_list.back().id = send_id;

  spdlog::debug(
      "sending to {} with tag {}", send_list.back().to_rank,
      send_list.back().comm_tag);
  /* Non-blocking send of data */
  err += MPI_Isend(
      send_list.back().send_buff_double.data(), num_elements, MPI_DOUBLE,
      send_list.back().to_rank, send_list.back().comm_tag, cart_comm,
      &(send_list.back().request));

  return err;
}

// ##################################################################
// ##################################################################

int Sub_domain::receive_double_data(
    const int from_rank,    ///< rank of process we are receiving from.
    const int comm_tag,     ///< comm_tag: what sort of comm we are looking for
                            ///< (PER,MPI,etc.)
    const string &recv_id,  ///< identifier for receive, for any book-keeping.
    const long int num_elements,  ///< number of doubles to receive
    vector<double>
        &data  ///< Pointer to array to write to (must be already initialised).
)
{
  int err = 0;
#ifndef NDEBUG
  spdlog::debug("Sub_domain::receive_double_data: starting");
#endif  // TEST_COMMS

  auto info =
      find_if(recv_list.begin(), recv_list.end(), [recv_id](const auto &ri) {
        return ri.id == recv_id;
      });

  if (info == recv_list.end()) {
    spdlog::error("{}: {}", "Failed to find recv with id:", recv_id);
    exit(1);
  }
#ifndef NDEBUG
  spdlog::debug(
      "found recv id={} and looking for id={}\nSub_domain::receive_double_data: found recv id",
      info->id, recv_id);
#endif  // TEST_COMMS

  /* Check that everthing matches */
  if (from_rank != info->from_rank) {
    spdlog::error(
        "{}: {}", "from_ranks don't match", from_rank - info->from_rank);
    exit(1);
  }
  if (comm_tag != info->comm_tag) {
    spdlog::error("{}: {}", "Bad comm_tag", comm_tag);
    exit(1);
  }
  if (recv_id != info->id) {
    spdlog::error("{}: {}", "Bad id", recv_id);
    exit(1);
  }
  if (info->type != COMM_DOUBLEDATA) {
    spdlog::error("{}: {}", "data is not double array!", info->type);
    exit(1);
  }

  /*
   * See how many doubles we are getting.
   * int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int
   * *count)
   */
  int buffer_size;
  err += MPI_Get_count(&(info->status), MPI_DOUBLE, &buffer_size);
  if (err) {
    spdlog::error("{}: {}", "getting size of message to receive", err);
    exit(1);
  }
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_double_data: buffer is {} elements", buffer_size);
#endif  // TEST_COMMS
  if (buffer_size != num_elements) {
    spdlog::debug(
        "Sub_domain::receive_double_data: getting {} doubles, but expecting {}",
        buffer_size, num_elements);
    spdlog::error(
        "{}: {}", "Getting wrong amount of data", buffer_size - num_elements);
    exit(1);
  }

  /*
   * Now Receive the data
   * int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int
   * tag, MPI_Comm comm, MPI_Status *status)
   */
#ifndef NDEBUG
  spdlog::debug(
      "Sub_domain::receive_double_data: receiving buffer of data from rank: {}",
      from_rank);
#endif  // TEST_COMMS
  err += MPI_Recv(
      data.data(), buffer_size, MPI_DOUBLE, info->from_rank, info->comm_tag,
      cart_comm, &(info->status));
  if (err) {
    spdlog::error("{}: {}", "MPI_Recv failed", err);
    exit(1);
  }

  if (recv_list.size() == 1) {
    recv_list.pop_front();
  }
  else {
    spdlog::error("{}: {}", "recv list is big!", recv_list.size());
    exit(1);
  }

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
  m_silo_id = session_id;
  int tag;
  PMPIO_iomode_t iomode = PMPIO_READ;
  if (iotype == "READ") {
    tag    = PMPIO_READ;
    iomode = PMPIO_READ;
  }
  else if (iotype == "WRITE") {
    tag    = PMPIO_WRITE;
    iomode = PMPIO_WRITE;
  }
  else {
    spdlog::error("{}: {}", "Bad iotype", iotype);
    exit(1);
  }

  m_baton = PMPIO_Init(
      num_files, iomode, cart_comm, tag, PMPIO_DefaultCreate, PMPIO_DefaultOpen,
      PMPIO_DefaultClose, 0);
  if (!m_baton) {
    spdlog::error("{}: {}", "failed to init files", myrank);
    exit(1);
  }

  /* Get the index of the group I am in */
  *group_rank = PMPIO_GroupRank(m_baton, myrank);

  /* Get the rank of my processor within the group */
  *rank_in_group = PMPIO_RankInGroup(m_baton, myrank);

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

  if (id != m_silo_id) {
    spdlog::error("{}: {}", id, m_silo_id);
    exit(1);
  }
  if (*dbfile) {
    spdlog::error("please pass in null file pointer!: {}", fmt::ptr(*dbfile));
    exit(1);
  }
  if (!m_baton) {
    spdlog::error("call init before wait_for_file!: {}", fmt::ptr(m_baton));
    exit(1);
  }

  *dbfile = static_cast<DBfile *>(
      PMPIO_WaitForBaton(m_baton, filename.c_str(), dir.c_str()));
  if (!(*dbfile)) {
    spdlog::error(
        "wait for baton failed to return file pointer: {}", fmt::ptr(*dbfile));
    exit(1);
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

  if (id != m_silo_id) {
    spdlog::error("{} != {}", id, m_silo_id);
    exit(1);
  }
  if (!(*dbfile)) {
    spdlog::error("file pointer is null!: {}", fmt::ptr(*dbfile));
    exit(1);
  }

  PMPIO_HandOffBaton(m_baton, static_cast<void *>(*dbfile));
  PMPIO_Finish(m_baton);

  m_baton = 0;
  m_silo_id.erase();
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
  if (id != m_silo_id) {
    spdlog::error("{}: {}", id, m_silo_id);
    exit(1);
  }

  *group = PMPIO_GroupRank(m_baton, proc);
  *rank  = PMPIO_RankInGroup(m_baton, proc);
}

#endif  // SILO

// ##################################################################
// ##################################################################

#endif  // PARALLEL
