//
// Threadpool demo
//
// Leak trace in MSCRT.
//
#include "msvc_constants.h"
#if defined(_DEBUG) &&  defined(_MSC_VER) &&  defined(MSVC_DEBUG_NEW_TRACE_ON)
  #define CRTDBG_MAP_ALLOC
  #include <stdlib.h> 
  #include <crtdbg.h> 
  #define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif
//
// Includes
//
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "constants.h"
#include "logmessages.h"
#include "threadpool/threadpool.h"
//
#include <new>
#include <iomanip>
#include <sstream>
//
using namespace std;
//
// Defines.
//
//
// Definitions.
//
struct testfunc_args {
  //
  // Arguments for Dispatch().
  //
  int id;
  int count;
};
//
// Global variables.
//
threadpool_t     tp;                // main threadpool
int monsecs_gl=0; // seconds since the start of the month
//
// Function prototypes
//
void testfuncW(void *arg);
int testfunc(const int count,const int id);
long get_random_int(long range);
void seed_rand_secs_month();
//
// External variables.
//
//
// External functions.
//
//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//        Main program flow begins here.                                   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//





int main(int argc,char **argv)
{
  //
  // Seed rand().
  //
  seed_rand_secs_month();
  //
  // Start threadpool.
  //
  tp_init(&tp,NUM_THREADS_MAIN,"Main Threadpool");
//   tp_close(&tp);    // <----- CAUSES SEGFAULT INVESTIGATE SOMETIME!!!!

//   int nprocs=NUM_THREADS_MAIN;
  int nprocs=10;

  DbgMsg(" main(): starting %i threads...",nprocs);
  for (int i=0;i<nprocs;i++) {
    //
    // Set a random count between 10 and 30 secs for each thread.
    //
    int count=10+get_random_int(21);
    //
    testfunc_args *ta=NULL;
    try   {ta=new testfunc_args;}
    catch (bad_alloc xa) {
      Abort(" main(): allocation failure.");
      return false;
    }
    ta->id=i;
    ta->count=count;
    //
    tp_addWork(&tp,testfuncW,(void *) ta,"main()");
  }
  DbgMsg(" main(): waiting for %i threads...",nprocs);
  tp_waitOnFinished(&tp,nprocs);
  DbgMsg(" main(): all threads finished.");
  return 0;
}


long get_random_int(long range) 
{
  //
  // The highest value we should get is range-1. 
  //
  return rand()%range;
}


void testfuncW(void *arg)
{
  testfunc_args *ta=(testfunc_args *) arg;
  testfunc(ta->count,ta->id);
  delete ta;
}
int testfunc(const int count,const int id)
{

  for (int i=0;i<count;i++) {
    DbgMsg(" I am thread %i at %i/%i secs",id,i,count);
    sleep(SLEEP_SEC);
  }
  DbgMsg(" -------------- I am thread %i finishing after %i secs",id,count);
  //
  return 1;
}

void seed_rand_secs_month()
{
  //
  // Calculates the number of seconds since the start of the month
  //  (assuming every month has 30 days in it!) and uses this value to
  //  seed rand().
  //
  time_t starttime=time(&starttime);
  time_t secsmon=30*((int) DAYS); // no of seconds in 30 days.
  time_t nmons  =starttime/secsmon;
  time_t starttime_mons=nmons*secsmon;
  if (starttime<starttime_mons) starttime_mons-=secsmon;
  int monsecs=starttime-starttime_mons;
  monsecs_gl=monsecs;
  srand(monsecs);
  DbgMsg(" Seeding rand() with: %i",monsecs);
}
