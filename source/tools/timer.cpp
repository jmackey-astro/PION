///
/// \file timer.cc
///
/// \author Jonathan Mackey
///
/// Class for measuring time spent in different parts of the code.
///
/// Modifications:
/// - 2015.03.02 JM: Copied from global.cc GeneralStuff class.

#include <sys/time.h>
//#include <time.h>
#include "timer.h"

using namespace std;

class timers clk;

//------------------------------------------------
//--- General Stuff Class ------------------------
timers::timers()
{
    clocks.clear();
}

timers::~timers()
{
    // cout <<"timers.size: "<<timers.size()<<"\n";
    clocks.clear();
}

// ##################################################################
// ##################################################################

void timers::start_timer(string id)
{
    struct timeval s;
    gettimeofday(&s, 0);
    double t = s.tv_sec + 1.e-6 * s.tv_usec;

    //
    // If timer exists, find it, and restart the timer.  It should have been
    // paused previously, in which case it's current value is the time it has
    // been running so far, so we set the value to the current time minus its
    // value, which is an 'effective' start time.
    // If the timer doesn't exist, initialise it to the current time in seconds.
    //
    if (clocks.find(id) != clocks.end()) {
        clocks[id] = t - clocks[id];
    }
    else {
        clocks[id] = t;
    }
    // cout << "id="<<id<<" start="<<clocks[id]<<"\n";
    return;
}

// ##################################################################
// ##################################################################

double timers::pause_timer(string id)
{
    struct timeval s;
    gettimeofday(&s, 0);
    double t = s.tv_sec + 1.e-6 * s.tv_usec;
    //
    // Set timer to be the number of seconds it has been running.
    //
    clocks[id] = t - clocks[id];
    return clocks[id];
}

// ##################################################################
// ##################################################################

double timers::stop_timer(string id)
{
    struct timeval s;
    gettimeofday(&s, 0);
    double t = s.tv_sec + 1.e-6 * s.tv_usec;
    // cout <<" start="<<clocks[id];
    // cout <<" end  ="<<t;
    t -= clocks[id];

    //
    // Delete the timer.
    //
    clocks.erase(id);
    // cout <<" time="<<t<<"\n";
    return t;
}

// ##################################################################
// ##################################################################

double timers::time_so_far(string id)
{
    struct timeval s;
    gettimeofday(&s, 0);
    double t = s.tv_sec + 1.e-6 * s.tv_usec;
    // cout <<" start="<<clocks[id];
    // cout <<" now ="<<t;
    // cout <<" timesofar="<<t-clocks[id];
    return t - clocks[id];
}

//------------------------------------------------
