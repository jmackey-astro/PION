///
/// \file reporting.cpp
/// \author Jonathan Mackey
///
/// This file defines the reporting class, for dealing with standard
/// input/output to screen.
///
/// Modifications:
/// - 2015.01.08 JM: created file, moved class from global.h
/// - 2017.12.09 JM: got rid of SimPM refs.

#include "comms/comms.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"

using namespace std;

///
/// global instance of the reporting class, for anyone that wants
/// to use it.
///
class reporting rep;

// ##################################################################
// ##################################################################

reporting::reporting()
{
  //  cout <<"Default reporting Constructor. O/P goes to cout/cerr.\n";
}

// ##################################################################
// ##################################################################

reporting::~reporting()
{
#if defined(SERIAL)
  cout.rdbuf(saved_buffer_cout);
  infomsg.close();
  //
  // cout<<"Deleting reporting class, This should be stdout.\n";
  //
  cerr.rdbuf(saved_buffer_cerr);
  errmsg.close();
  //
  // cout<<"Deleting reporting class, This should be stderr.\n";
  //
#elif defined(PARALLEL)
  if (myrank == 0) {
    cout.rdbuf(saved_buffer_cout);
    infomsg.close();
  }
  else {
    std::cout.clear();
  }
#else
#error "Must define either SERIAL or PARALLEL (reporting::~reporting)"
#endif
}

// ##################################################################
// ##################################################################

int reporting::redirect(const string& path)
{
  string temp;

#if defined(SERIAL)
#ifdef TESTING
  cout << "(reporting::redirect): O/P goes to text files in " << path << "\n";
#endif
  temp = path + "errors.txt";
  errmsg.open(temp.c_str(), ios::trunc);
  if (!errmsg.is_open()) {
    cerr << "Reporting: can't open errors.txt for writing.\n";
    exit(1);
  }
  errmsg.copyfmt(cerr);
  saved_buffer_cerr = cerr.rdbuf();
  cerr.rdbuf(errmsg.rdbuf());

  temp = path + "info.txt";
  infomsg.open(temp.c_str(), ios::trunc);
  if (!infomsg.is_open()) {
    cerr << "Reporting: can't open info.txt for writing.\n";
    exit(1);
  }
  infomsg.copyfmt(cout);
  saved_buffer_cout = cout.rdbuf();
  cout.rdbuf(infomsg.rdbuf());

  cout.setf(ios_base::scientific);
  cout.precision(7);

#elif defined(PARALLEL)
#ifndef TESTING
  //
  // For parallel execution (production runs) we only want a single
  // log file, and errors should be printed to stderr.
  //
  myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  // cout <<"myrank="<<myrank<<"\n";
#ifdef REPORT_RANK0
  if (myrank == 0) {
#endif
    // cout <<"(reporting::redirect): O/P goes to text files in ";
    // cout <<path<<"\n";
    // cout <<"Note: not redirecting error messages, and suppressing ";
    // cout <<"stdout from all processes except myrank=0.\n";

    temp = path + "info.txt";
    infomsg.open(temp.c_str(), ios::trunc);
    if (!infomsg.is_open()) {
      cerr << "Reporting: can't open info.txt for writing.\n";
      exit(1);
    }
    infomsg.copyfmt(cout);
    saved_buffer_cout = cout.rdbuf();
    cout.rdbuf(infomsg.rdbuf());
    cout.setf(ios_base::scientific);
    cout.precision(7);
#ifdef REPORT_RANK0
  }
#endif
#else
  //
  // for testing we want all processors to have their own log file.
  //
  temp = path + "info.txt";
  infomsg.open(temp.c_str(), ios::trunc);
  if (!infomsg.is_open()) {
    cerr << "Reporting: can't open info.txt for writing.\n";
    exit(1);
  }
  infomsg.copyfmt(cout);
  saved_buffer_cout = cout.rdbuf();
  cout.rdbuf(infomsg.rdbuf());
  cout.setf(ios_base::scientific);
  cout.precision(7);

  temp = path + "errors.txt";
  errmsg.open(temp.c_str(), ios::trunc);
  if (!errmsg.is_open()) {
    cerr << "Reporting: can't open errors.txt for writing.\n";
    exit(1);
  }
  errmsg.copyfmt(cerr);
  saved_buffer_cerr = cerr.rdbuf();
  cerr.rdbuf(errmsg.rdbuf());
#endif  // TESTING
#else
#error "Must define either SERIAL or PARALLEL (reporting::redirect)"
#endif
  return (0);
}

// ##################################################################
// ##################################################################

void reporting::kill_stdout_from_other_procs(
    const int core  ///< don't kill it from this core
)
{
#if defined(PARALLEL)
  //
  // For parallel execution (production runs) we only want ouput
  // from proc. 0, so we redirect all stdout/stderr from other
  // processes to some dead end like /dev/null.
  //
  myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  if (myrank != core) {
    // saved_buffer_cout = cout.rdbuf(); // <-- save
    // cout.rdbuf (nullstream.rdbuf());  // <-- redirect
    std::cout.setstate(std::ios::failbit);
  }
#endif
  return;
}

// ##################################################################
// ##################################################################
