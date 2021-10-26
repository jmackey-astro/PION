/// /// \file reporting.h /// \author Jonathan Mackey ///
/// This file declares the reporting class, for dealing with standard
/// input/output to screen.
///
/// Modifications:
/// - 2015.01.08 JM: created file, moved class from global.h
/// - 2015.01.16 JM: extra flush/endl statements in error().
/// - 2015.01.26 JM: added include for comms and cell classes.
/// - 2017.12.09 JM: got rid of SimPM refs.

#ifndef REPORTING_H
#define REPORTING_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_params.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

///
/// Global Class For Writing Error messages and Reports.
///
/// This determines where to write different messages to.
///
class reporting {
public:
  reporting();   ///< Default Constructor, write to std[out/err].
  ~reporting();  ///< Destructor.

  ///
  /// Redirects stdout/stderr to files in the path specified.
  ///
  int redirect(const string &  ///< Location of files to write reporting to.
  );

  ///
  /// This kills all output to screen from all cores except
  /// the one listed.
  ///
  void kill_stdout_from_other_procs(const int  ///< don't kill it from this core
  );

  ///
  /// This exits the code cleanly if I know I have an error.
  ///
  ///
  /// This exits from the code, printing an error message and the
  /// offending value.
  /// I would like to have this in reporting.cpp, but it's a template
  /// function, and the export keyword is not implemented in any
  /// compilers, so I have to have the definitions in the header.
  ///
  template<class T>
  inline void error(
      string msg,  ///< Error Message.
      T err        ///< Value which is wrong.
  )
  {
    cerr << msg << "\t error code: " << err << " ...exiting.\n";
    cout.flush();
    cout << endl;  // try to write all to stdout/file.
    cerr.flush();
    cerr << endl;  // try to write all to stderr.

    exit(1);
  }

  ///
  /// Prints a warning message but don't stop execution.
  ///
  template<class T1, class T2>
  inline void warning(
      string msg,   ///< Message.
      T1 expected,  ///< Expected Value
      T2 found      ///< Received Value
  )
  {
    cout << "WARNING: " << msg << "\t Expected: " << expected;
    cout << " but got " << found << endl;
    return;
  }

  ///
  /// Error if actual value is less than Expected Value.
  ///
  /// T1 and T2 must be compatible types for testing 'less-than' operation.
  ///
  template<class T1, class T2>
  inline void errorLT(
      string msg,  ///< Message.
      T1 exptd,    ///< Expected Value
      T2 recvd     ///< Actual Value
  )
  {
    if (recvd < exptd) {
      cout << "ERROR DETECTED: " << msg << "\t expected less than " << exptd;
      cout << " but got " << recvd << endl;
      error(msg, recvd);
    }
    return;
  }

  ///
  /// Tests if two values are the same, and error if not.
  ///
  /// If the expected value is not equal to the actual value, then
  /// print an error message and exit.
  ///
  inline void errorTest(
      string msg,  ///< Error Message.
      int exptd,   ///< Expected Value.
      int recvd    ///< Actual Value.
  )
  {
    if (exptd != recvd) {
      cerr << "ERROR DETECTED: " << msg << "\t expected " << exptd;
      cerr << " but got " << recvd << endl;
      error(msg, recvd);
    }
    return;
  }

  ///
  /// Print out a vector.
  ///
  template<class T>
  inline void printVec(
      string msg,  ///< Name of Vector
      T *vec,      ///< pointer to vector
      int nd       ///< length of vector.
  )
  {
    cout << "Vector " << msg << " : [";
    for (int i = 0; i < nd - 1; i++)
      cout << vec[i] << ", ";
    cout << vec[nd - 1] << " ]\n";
    return;
  }

  ///
  /// Print out a vector.
  ///
  template<class T>
  inline void printVec(
      string msg,                ///< Name of Vector
      const std::vector<T> &vec  ///< pointer to vector
  )
  {
    int nd = vec.size();
    cout << "Vector " << msg << " : [";
    for (int i = 0; i < nd - 1; i++)
      cout << vec[i] << ", ";
    cout << vec[nd - 1] << " ]\n";
    return;
  }

#ifdef PARALLEL
  void set_rank(const int rank) { myrank = rank; }
#endif

private:
  ofstream errmsg, infomsg;  //, iomsg;
  streambuf *saved_buffer_cout, *saved_buffer_cerr;
#ifdef PARALLEL
  int myrank;  ///< rank in list of MPI processes
#endif         // PARALLEL
};

extern class reporting rep;

#endif  // REPORTING_H
