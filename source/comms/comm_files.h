/* \file comm_files.h
 *
 * \brief Contains comms class for multi-process communication using
 * files.  This is only intended for debugging code -- all processes
 * must be started up by hand on the same machine in the same directory,
 * and this class creates and destroys files to pass information between
 * processes.
 *
 * \author Jonathan Mackey
 * \date 2009-01-27.
 */ 
/// Modifications:
/// - 2012.05.15 JM: Added function for global-reduce (max/min/sum) of arrays.
///

#ifdef PARALLEL
#ifdef USE_FILE_COMMS

#ifndef COMM_FILES_H
#define COMM_FILES_H

//
// These tells code what to compile and what to leave out.
//
#include "../defines/functionality_flags.h"
#include "../defines/testing_flags.h"

#include "comms.h"
#include "../dataIO/file_status.h"

#define FDELAY_USECS 500000

struct sent_info {
  std::string id;      ///< string that code uses as handle for send.
  int comm_tag;        ///< the tag used to describe the send.
  int from_rank;       ///< sender.
  int to_rank;         ///< recipient
  int type;            ///< type of data: COMM_CELLDATA=char array; COMM_DOUBLEDATA=double array.
};

struct recv_info {
  std::string id;    ///< string that code uses as handle for send.
  int comm_tag;      ///< the tag used to describe the send.
  int from_rank;     ///< sender.
  int to_rank;       ///< recipient
  int type;          ///< type of data: COMM_CELLDATA=char array; COMM_DOUBLEDATA=double array.
};

class comm_files : public comms_base {
 public:
  comm_files();
  ~comm_files();

  /** \brief Tell other processes that I exist, and wait until all of them exist too. 
   * Also set myrank and nproc in the global struct mpiPM.  The last two command line
   * args should be "myrank=i nproc=n", and this function strips off those args.
   **/
  int init(int *, ///< number of program arguments.
	   char *** ///< character list of arguments.
	   );

  /** \brief Tell other processes that I am exiting, and either exit, or wait for the
   * others and then exit. */
  int finalise();
  /** \brief Tell other processes to abort! */
  int abort();

  /** \brief Set up a barrier, and don't return until all other processes have also 
   * set up their own barrier. */
  int barrier(const std::string);

  /** \brief Do a global operation on local data, options are MAX,MIN,SUM.
   * Returns the global value when local value is passed in.
   */
  double global_operation_double(const std::string, ///< MAX,MIN,SUM
				 const double  ///< this process's local value.
				 );
 
  /// 
  /// Do a global operation on local data, options are MAX,MIN,SUM,but easy to
  /// add in others if needed.  Receives in local value, performs global 
  /// operation on it, and returns global value.  This function works with
  /// N-element arrays, and should work "in-place" i.e. send and recv buffers
  /// are the same.
  ///
  void global_op_double_array(
                      const std::string, ///< MAX,MIN,SUM
                      const size_t,      ///< Number of elements in array.
		      double * ///< pointer to this process's data array.
		      );

  /** \brief Broadcast data from one process to all others. */
  int broadcast_data(const int,         ///< rank of sender.
		     const std::string, ///< Type of data INT,DOUBLE,etc.
		     const int,         ///< number of elements
		     void *             ///< pointer to data.
		     );

  /** \brief Send cell data to another processor, and return immediately, but
   * keep a record of the send so that I can tell later when the send has been 
   * received.
   *
   * This receives a list of cells, extracts specific data from each cell into a 
   * buffer to send, and sends it to another processor, which is expected to 
   * know how to unpack the data from it's matching receive call.
   */
  int send_cell_data(const int,           ///< rank to send to.
		     std::list<cell *> *, ///< list of cells to get data from.
		     long int,            ///< number of cells in list (extra checking!)
		     std::string &,            ///< identifier for send, for tracking delivery later.
		     const int            ///< comm_tag, to say what kind of send this is.
		     );

  /** \brief Send an array of n doubles to another processor, and return immediately, but
   * keep a record of the send so that I can tell later when the send has been received.
   */
  int send_double_data(const int,      ///< rank to send to.
		       const long int, ///< size of buffer, in number of doubles.
		       const double *, ///< pointer to double array.
		       std::string &,       ///< identifier for send, for tracking delivery later.
		       const int       ///< comm_tag, to say what kind of send this is.
		       );
  /** \brief This function is called when we need to make sure a send has been received.
   * It only returns once the receiver has confirmed that it has got the data.
   */
  int wait_for_send_to_finish(std::string & ///< identifier for the send we are waiting on.
			      );

  /** \brief Look for some data that is being sent to us.  Does not return
   * until it finds some, so it is up to the programmer to prevent deadlock!
   **/
  int look_for_data_to_receive(int *,    ///< rank of sender
			       std::string &, ///< identifier for receive.
			       int *,     ///< comm_tag associated with data.
			       const int  ///< type of data to look for (COMM_CELLDATA,COMM_DOUBLEDATA)
			       );
  /** \brief Receive Cell data from a specific process rank. 
   *
   * It is up to the caller to make sure that data was packed in the right way, and that the cells
   * were sent in the same order that they are received.  No checking of this is performed.
   **/
  int receive_cell_data(const int,           ///< rank of process we are receiving from.
			std::list<cell *> *, ///< list of cells to get data for. 
			const long int,      ///< number of cells in list (extra checking!)
			const int,           ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
			const std::string &  ///< identifier for receive, for any book-keeping that might be needed.
			);

  /** \brief Receive array of doubles from a specific process rank. 
   *
   * It is up to the caller to make sure that it knows what to do with the list!
   **/
  int receive_double_data(const int,      ///< rank of process we are receiving from.
			  const int,      ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
			  const std::string &, ///< identifier for receive, for any book-keeping that might be needed.
			  const long int, ///< number of doubles to receive
			  double *        ///< Pointer to array to write to (must be already initialised).
			  );
#ifdef SILO
  int silo_pllel_init(const int,    ///< number of files to write
		      const std::string, ///< READ or WRITE
		      const std::string,  ///< identifier for this read/write.
		      int *, ///< group rank.
		      int * ///< rank in group
		      );
  int silo_pllel_wait_for_file(const std::string, ///< identifier for this read/write.
			       const std::string, ///< File Name
			       const std::string, ///< Directory to open in file.
			       DBfile ** ///< pointer that file gets returned in.
			       );
  int silo_pllel_finish_with_file(const std::string, ///< identifier for this read/write.
				  DBfile ** ///< pointer to file we have been working on.
				  );
  void silo_pllel_get_ranks(std::string id, ///< identifier for this read/write.
			    const int, ///< processor rank
			    int *, ///< rank of group processor is in.
			    int *  ///< rank of processor within group.
			    );
#endif
 private:
  std::list<struct sent_info *> sent_list;
  std::list<struct recv_info *> recv_list;
#ifdef SILO
  //PMPIO_baton_t *bat;
  std::string silo_id; ///< identifier for this silo session.
  int 
    num_files,     ///< how many files to write in total
    grp_rank,      ///< my rank in the list of files to write
    rank_in_grp,   ///< my rank in the file i am to write
    silo_filetype; ///< either DB_PDB or DB_HDF5
  std::string silo_iotype; ///< whether we are to read or write...
#endif
  int myrank;
  int nproc;
  int barrier_count; ///< counter to uniquely identify each barrier we come up against.
  std::string dir; ///< directory name for filecomms.
  std::string ini; ///< name for file that root generates on init().
  std::string bar; ///< name for barrier file.
  std::string go_order;    ///< name for file telling slaves to go.
  std::string go_response; ///< name for file slaves create to respond to go order.
  std::string celldata;    ///< identifier for file containing cell data.
  std::string doubledata;  ///< identifier for file containing double data.
  std::string comm_global; ///< identifier for file to be read globally (for max,min,sum,etc.)
  std::string finished_with_file; ///< identifier to be appended to file when a proc has finished reading a file.
  //std::string ; ///<

  class file_status fs; ///< for checking if file exists.

  void wait_for_file(const std::string);
  void wait_for_file_to_disappear(const std::string);
  void master_send_go_signal(const std::string ///< identifier for the current operation we are organising
			     );
  void slave_recv_go_signal(const std::string ///< identifier for the current operation we are organising
			    );
  void master_wait_on_slave_files(const std::string, ///< message.
				  const std::string  ///< what slave files are called.
				  );
  void wait_for_peer_to_read_file(const std::string ///< name of file
				  );
};



#endif // COMM_FILES_H

#endif // USE_FILE_COMMS
#endif //PARALLEL
