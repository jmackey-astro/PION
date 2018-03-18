/// \file comms.h
/// \brief Contains abstract base class for multi-process communication.
///
/// \author Jonathan Mackey
///
/// Modifications:
/// - 2012.05.15 JM: Added function for global-reduce (max/min/sum) of arrays.
/// - 2015.01.26 JM: added get_rank_nproc() function, and pointer.

#ifdef PARALLEL

#ifndef COMMS_H
#define COMMS_H

//
// These tells code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <list>

#ifdef SILO
#include <silo.h>
#endif

#include "grid/cell_interface.h"

#define COMM_CELLDATA   1
#define COMM_DOUBLEDATA 2

class comms_base {
 public:
  /** \brief Virtual destructor. */
  virtual ~comms_base() {}

  /** \brief Tell other processes that I exist, and wait until all of them exist too.  
   * Also set myrank and nproc in the global struct mpiPM.
   **/
  virtual int init(int *, ///< number of program arguments.
		   char *** ///< character list of arguments.
		   )=0;
  
  ///
  /// Get this process's rank, and total number of processes.
  ///
  virtual int get_rank_nproc(
        int *, ///< rank.
	int *  ///< nproc
	)=0;

  /** \brief Tell other processes that I am exiting, and either exit, or wait for the
   * others and then exit. */
  virtual int finalise()=0;

  /** \brief Tell other processes to abort! */
  virtual int abort()=0;

  /** \brief Set up a barrier, and don't return until all other processes have also 
   * set up their own barrier. */
  virtual int barrier(const std::string)=0;

  /** \brief Do a global operation on local data, options are MAX,MIN,SUM,but easy to
   * add in others if needed.  Receives in local value, performs global operation on it,
   * and returns global value.
   */
  virtual double global_operation_double(const std::string, ///< MAX,MIN,SUM
					 const double  ///< this process's max/min value.
					 )=0;

  /// 
  /// Do a global operation on local data, options are MAX,MIN,SUM,but easy to
  /// add in others if needed.  Receives in local value, performs global 
  /// operation on it, and returns global value.  This function works with
  /// N-element arrays, and should work "in-place" i.e. send and recv buffers
  /// are the same.
  ///
  virtual void global_op_double_array(
                      const std::string, ///< MAX,MIN,SUM
                      const size_t,      ///< Number of elements in array.
		      double * ///< pointer to this process's data array.
		      )=0;

  /** \brief Broadcast array of data from one process to all others. 
   *
   * Data can't be const because n-1 of the processes have to write to it!
   **/
  virtual int broadcast_data(const int,         ///< rank of sender.
			     const std::string, ///< Type of data INT,DOUBLE,etc.
			     const int,         ///< number of elements
			     void *             ///< pointer to data.
			     )=0;
			     
  /// 
  /// \brief Send cell data to another processor, and return immediately, but
  /// keep a record of the send so that I can tell later when the send has been 
  /// received.
  ///
  /// This receives a list of cells, extracts specific data from each cell into a 
  /// buffer to send, and sends it to another processor, which is expected to 
  /// know how to unpack the data from it's matching receive call.
  ///
  virtual int send_cell_data(
      const int,           ///< rank to send to.
      std::list<cell *> *, ///< list of cells to get data from.
      long int,            ///< number of cells in list
      const int,           ///< ndim
      const int,           ///< nvar
      std::string &,       ///< identifier for send, for tracking delivery later.
      const int            ///< comm_tag, to say what kind of send this is.
      )=0;

  /** \brief Send an array of n doubles to another processor, and return immediately, but
   * keep a record of the send so that I can tell later when the send has been received.
   *
   * The buffer can be deleted on return, since the comms class makes a copy of it if it
   * needs to keep it for any reason (e.g. MPI needs to store buffers on its own timescale).
   */
  virtual int send_double_data(const int,      ///< rank to send to.
			       const long int, ///< size of buffer, in number of doubles.
			       const double *, ///< pointer to double array.
			       std::string &,  ///< identifier for send, for tracking delivery later.
			       const int       ///< comm_tag, to say what kind of send this is.
			       )=0;

  /** \brief This function is called when we need to make sure a send has been received.
   * It only returns once the receiver has confirmed that it has got the data.
   */
  virtual int wait_for_send_to_finish(std::string & ///< identifier for the send we are waiting on.
				      )=0;

  /** \brief Look for some data that is being sent to us.  Does not return
   * until it finds some, so it is up to the programmer to prevent deadlock!
   *
   * When it finds some data of the type requested, it reads the comm_tag of the data,
   * the process who is sending it, and associates an identifier string with it.
   * This string can be used to receive the data by a subsequent function call.
   **/
  virtual int look_for_data_to_receive(int *,         ///< rank of sender
				       std::string &, ///< identifier for receive.
				       int *,         ///< comm_tag associated with data.
				       const int      ///< type of data to look for (COMM_CELLDATA,COMM_DOUBLEDATA)
				       )=0;

  ///
  /// \brief Receive Cell data from a specific process rank. 
  ///
  /// It is up to the caller to make sure that data was packed in the right way, and that the cells
  /// were sent in the same order that they are received.  No checking of this is performed.
  ///
  virtual int receive_cell_data(
      const int,           ///< rank of process we are receiving from.
      std::list<cell *> *, ///< list of cells to get data for. 
      const long int,      ///< number of cells in list (extra checking!)
      const int,           ///< ndim
      const int,           ///< nvar
      const int,           ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
      const std::string &  ///< identifier for receive, for any book-keeping that might be needed.
      )=0;

  /** \brief Receive array of doubles from a specific process rank. 
   *
   * It is up to the caller to make sure that it knows what to do with the data received!
   * No ordering is enforced, so the calling code needs to ensure that sender and receiver
   * order the data in the same way.
   **/
  virtual int receive_double_data(const int,           ///< rank of process we are receiving from.
				  const int,           ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
				  const std::string &, ///< identifier for receive, for any book-keeping that might be needed.
				  const long int,      ///< number of doubles to receive
				  double *             ///< Pointer to array to write to (must be already initialised).
				  )=0;

#ifdef SILO
  /** \brief Silo Parallel I/O interface.  init() assigns me a group and a rank within
   * the group.  Each group writes one file, one process at a time.  So the number of 
   * groups/files is the number of processes writing/reading concurrently.  Optimal write 
   * speed is hard to predict, but is probably achieved with 2-8 files.
   *
   * In theory this interface could work for any filetype with the pointer types changed and
   * new file create/open/close functions defined.
   */ 
  virtual int silo_pllel_init(const int,         ///< number of files to write
			      const std::string, ///< READ or WRITE
			      const std::string, ///< identifier for this read/write.
			      int *,             ///< group rank.
			      int *              ///< rank in group
			      )=0;

  /** \brief Each process in a group has to wait it's turn to read/write a file. This
   * function determines how this happens.  e.g. a 'baton' can be passed around from the
   * first to last process in the group.
   **/
  virtual int silo_pllel_wait_for_file(const std::string, ///< identifier for this read/write.
				       const std::string, ///< File Name
				       const std::string, ///< Directory to open in file.
				       DBfile ** ///< pointer that file gets returned in.
				       )=0;

  /** \brief Once I have finished operating on the file, hand it off to the next process
   * to operate on, or close the file if I am the last process.
   **/
  virtual int silo_pllel_finish_with_file(const std::string, ///< identifier for this read/write.
					  DBfile ** ///< pointer to file we have been working on.
					  )=0;

  /** \brief If I need a group+rank for an arbitrary process rank, this function
   * calculates it.
   **/
  virtual void silo_pllel_get_ranks(const std::string id, ///< identifier for this read/write.
				    const int, ///< processor rank to get group and rank_in_group for.
				    int *,     ///< rank of group processor is in.
				    int *      ///< rank of processor within group.
				    )=0;
#endif

};


extern class comms_base *COMM;

#endif // COMMS_H

#endif //PARALLEL
