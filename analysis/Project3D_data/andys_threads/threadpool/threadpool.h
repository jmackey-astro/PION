//
// Threadpool for handling multiple execution threads.
//
#ifndef THREADPOOL_H
#define THREADPOOL_H
//
// Includes.
//
#include "../reefa_constants.h"
#include "../threads.h"
//
#include <string>
//
using namespace std;
//
// Defines.
//
#define TRUE_TP            1
#define FALSE_TP           !TRUE_TP
//
#define PRIORITY_LOWER  0
#define PRIORITY_NORMAL 1
//
#define THREADNAME_MAP_OFF // allows use of getThreadName() under MSVC dbg
//
// Definitions.
//
extern "C" {

struct threadWork_t {
#ifndef RELEASE_ON
  int dbg_id;
#endif
  string caller;
  string thread_tag;
  void (*func)(void *);
  void *arg;
  struct threadWork_t *next;	
};

struct threadpool_t {
  string dbg_name;     // label identifying this pool
  int priority;
  int numThread;
  pthread_t *threads;
  thread_mutex_t Qlock;
  thread_cond_t QNotEmpty;
  thread_cond_t TFinished;
  
  struct threadWork_t *Qhead;
  struct threadWork_t *Qtail;
  int shutdown;
  int numfinished;
};

//
// Wrapper for a conditional - we can use this to wait and
//  broadcast a wake up signal.
//
struct off_thread_cond_t {
  string dbg_name;     // label identifying this conditional
  thread_mutex_t lock;
  thread_cond_t cond;
};


int  tp_init(threadpool_t *tp,int nt,const char *name,
	     int priority=PRIORITY_LOWER,
	     bool name_thread=true);
int  tp_addWork(threadpool_t *tp,void (*func)(void *),void *args,
		const char *caller,const char *thrd_tag=NULL);
int  tp_waitOnFinished(threadpool_t *tp,int n);
void tp_reset(threadpool_t *tp);
void tp_close(threadpool_t *tp);

} // extern "C"
//
// Conditional functions.
//
off_thread_cond_t *ThreadCond_Create(const char *name,const long maxcount=1);
void ThreadCond_Destroy(off_thread_cond_t *tc);
void ThreadCond_Wait(off_thread_cond_t *tc);
void ThreadCond_TimedWait(off_thread_cond_t *tc,int secs);
void ThreadCond_Broadcast(off_thread_cond_t *tc);
//
// cygthreadpool - cygwin has a nasty memory leak that makes using
//  the threadpool a pain, memory gushes when a new task is added
//  to the threadpool. This pool sits on top of the normal threadpool
//  with each thread being a continuous task that monitors the a list
//  for new jobs. This will only be built under cygwin, which
//  is only used for dev - builds from MSVC, mac and linux should use
//  the a normal threadpool where this class is used under cygwin.
//
#if !defined(_CYGWIN)
//
// Non-cygwin build - alias the cygthreadpool to the normal threadpool.
//
#define cygthreadpool_t threadpool_t
#else
//
// cygthreadpool struct
//
struct cygthreadpool_t {
  string dbg_name;          // label identifying this pool
  long int wait_time;     // delay of pool loop in u-secs.
  int numThreads;
  thread_mutex_t Qlock;
  
  struct threadWork_t *Qhead;
  struct threadWork_t *Qtail;	
  int shutdown;
  thread_mutex_t shutdownLock;
  int numfinished;		  	  
  thread_mutex_t numfinishedLock;
};
//
// cygthreadpool functions
//
void CygTP_SetShutdown(cygthreadpool_t *tp);
int  CygTP_Shutdown(cygthreadpool_t *tp);
void CygTP_IncrNumFinished(cygthreadpool_t *tp,int incr);
int  CygTP_NumFinished(cygthreadpool_t *tp);
void CygTP_ThreadW(void *arg);
void CygTP_Thread(cygthreadpool_t *tp);
threadWork_t *CygTP_CheckForTask(cygthreadpool_t *tp);
threadWork_t *CygTP_GetTask(cygthreadpool_t *tp);
void CygTP_ExecuteTask(threadWork_t *tw,const char *dbg_name0);

#endif
//
// These funcs are used whether cygthreadpool_t is aliased or not.
//
int CygTP_Init(cygthreadpool_t *tp, int nt,long int wait_time,
	       const char *name,
	       int priority=PRIORITY_LOWER,bool name_thread=true);
int CygTP_Close(cygthreadpool_t *tp);
int CygTP_addWork(cygthreadpool_t *tp,void (*func)(void *),void *args,
		  const char *caller);
//
// MSVC specific debug funcs.
//
#if defined(_DEBUG) && defined(_MSC_VER)
#ifdef THREADNAME_MAP_ON
string getThreadName();
#endif
#endif


#endif



