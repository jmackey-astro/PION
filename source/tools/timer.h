/// \file timer.h
/// \author Jonathan Mackey
///
/// Class to managing timers, to measure how long the different parts
/// of the code are taking.
///
/// Modifications:
/// - 2015.03.02 JM: Moved from global GeneralStuff class.

#ifndef TIMER_H
#define TIMER_H

//
// These tell code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <iostream>
#include <map>
using namespace std;

///
/// Class for timing how long code spends doing certain tasks.
///
class timers {
public:
  timers();
  ~timers();

  ///
  /// start a timer, identified by a string. If timer exists
  /// already, this function assumes it has been paused and sets it
  /// counting again.
  ///
  void start_timer(string);

  ///
  /// Pause a timer, identified by a string, and return time in
  /// seconds.  Sets time value to number of seconds it has been
  /// running.
  ///
  double pause_timer(string);

  ///
  /// Return total time since timer started, and keep running.
  ///
  double time_so_far(string);

  ///
  /// stop a timer, identified by a string, and return time in
  /// seconds.  Deletes the timer.
  ///
  double stop_timer(string);

private:
  map<string, double> clocks;  ///< maps a string identifier for a timer to an
                               ///< index in the start vector.
};

extern class timers clk;  ///< the timer class.

#endif  // TIMER_H
