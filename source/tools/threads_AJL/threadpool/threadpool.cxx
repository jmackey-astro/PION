//
// Threadpool for handling multiple execution threads.
//
// This source code is released under the GNU General Public Licence:
//    http://www.gnu.org/copyleft/gpl.html
//
//
// Leak trace in MSCRT.
//
#include "../msvc_constants.h"
#if defined(_DEBUG) &&  defined(_MSC_VER) &&  defined(MSVC_DEBUG_NEW_TRACE_ON)
  #define CRTDBG_MAP_ALLOC
  #include <stdlib.h> 
  #include <crtdbg.h> 
  #define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "threadpool.h"
#include "../logmessages.h"
#if defined(_WIN32)
  #if defined(_MSC_VER)
    #if !defined(_MT)
      #error You must compile for multithreaded environment!
    #endif
  #endif
  #include <process.h>
#endif
//
#include <new>
#include <sstream>
#include <string>
//#include <iostream>
//
using namespace std;
//
// Defines.
//
//
// External variables.
//
extern int monsecs_gl;
//
// External functions.
//
#ifndef RELEASE_ON
extern void dbg_incr_nthr(int i);
extern void dbg_print_nthr(char *str,char *caller);
#endif
//
// Function protoypes.
//
void* threadFunc(void *ptr);
//
// We can wrap the thread function to get some extra Win32 features.
//
#if defined(_WIN32)
//
// MS debugger support.
// It's a crying shame named threads isn't part of the Win32 API.
//
#if defined(_DEBUG) && defined(_MSC_VER)

#ifdef THREADNAME_MAP_ON
#include <map>
//
using namespace std;
//
map<DWORD,string> thread_names;
#endif
//
//
//
typedef struct tagTHREADNAME_INFO
{
  DWORD  dwType;     // must be 0x1000
  LPCSTR szName;     // pointer to name (in user addr space)
  DWORD  dwThreadID; // thread ID (-1=caller thread)
  DWORD  dwFlags;    // reserved for future use, must be zero
} THREADNAME_INFO;
//
// Function to set the name of a thread.
//
void __stdcall setThreadName(DWORD dwThreadID,const char* szThreadName)
{


  THREADNAME_INFO info;
  info.dwType     = 0x1000;
  info.szName     = szThreadName;
  info.dwThreadID = dwThreadID;
  info.dwFlags    = 0;

  __try
    {
      RaiseException(0x406d1388,0,sizeof(info)/sizeof(DWORD),(DWORD*)&info);
    }
  __except (EXCEPTION_CONTINUE_EXECUTION)
    {
    }
}
//
// Flag is a debugger is present.
//
bool bDebuggerPresent=IsDebuggerPresent()!=0;
#endif // _DEBUG && _MSC_VER
//
// Win32 thread startup function wrapper. Basically just a signature adapter
// with some whistles and bells added for debugging support.
//
void set_msvc_debug_info(const string &tag,const bool use_num)
{
#if defined(_DEBUG) && defined(_MSC_VER)
  if (!bDebuggerPresent) return;
  static int n = 0;
  ++n;
  // NOTE: MAX 8 chars for MSVC6 debugger!
  stringstream ss;
  ss << tag[0] << tag[1] << tag[2];
  if (use_num) ss << n;
  DWORD dwThreadID=GetCurrentThreadId();
  setThreadName(dwThreadID,ss.str().c_str());
#ifdef THREADNAME_MAP_ON
  map<DWORD,string>::iterator t=thread_names.find(dwThreadID);
  if (t!=thread_names.end()) thread_names.erase(t); // rename thrd
  thread_names.insert(pair<DWORD,string>(dwThreadID,ss.str()));
#endif
#endif
}

#ifdef THREADNAME_MAP_ON
string getThreadName()
{
  if (!bDebuggerPresent) return "no dbgr";
  DWORD dwThreadID=GetCurrentThreadId();
  map<DWORD,string>::iterator t=thread_names.find(dwThreadID);
  if (t!=thread_names.end()) return t->second;
  return "no name!";
}
#endif



unsigned __stdcall threadFuncWin32(void* p)
{
  //
  // Set priority - usually lower for worker threads
  //
  threadpool_t *tp=(threadpool_t *) p;
  if (tp->priority==PRIORITY_LOWER) {
    SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_BELOW_NORMAL);
  }
  //
  // Some debugger stuff.
  //
  string tag="wrk";
  if      (tp->dbg_name=="Retr") tag="ret"; 
  else if (tp->dbg_name=="Stor") tag="sto"; 
  else if (tp->dbg_name=="Dist") tag="dst"; 
  else if (tp->dbg_name=="User") tag="usr"; 
  else if (tp->dbg_name=="DLMa") tag="dlm"; 
  set_msvc_debug_info(tag,true);
  //
  // Execute the task.
  //
  threadFunc(p);
  //
  // Win32: _endthread() or _endthreadex() is (apparently) called
  //  implicitly when the threaded func ends.
  //
  _endthreadex(0); // let's not trust them.
  return 0;
}

#endif // _WIN32


//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//   Main threadpool functions.                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//


int tp_init(threadpool_t *tp,int nt,const char *name,int priority,
	    bool name_thread) 
{
  //
  // Initialization function- creates the threads and sets up the 
  //  tpool struct.
  //
  tp->dbg_name   =name;
  tp->priority   =priority;
  tp->numThread  =nt;
  tp->numfinished=0;
  tp->shutdown   =0;
  tp->Qhead=tp->Qtail=NULL;
  //
  // Make array of threads.
  //
  try   {tp->threads=new pthread_t[nt];}
  catch (bad_alloc xa) {
    Abort(" tp_init(): can't create thread array - this is terminal :(");
    return false;
  }
  //
  // Tnitialise thread locks.
  //
  tp->QNotEmpty.init();
  tp->TFinished.init();
  //
  // Name the main thread too, else it's hell to find among all
  // threads in the threads list when debugging.
  // Note: As mentioned elsewhere, limit to 8 chars for MSVC6 debugger.
  //
#if defined(_WIN32) && defined(_DEBUG) && defined(_MSC_VER)
  if (name_thread&&bDebuggerPresent) {
    setThreadName(GetCurrentThreadId(),"MAINthrd");
  }
#endif
  //
  // Create the individual threads in the pool.
  //
  for (int i=0;i<nt;i++) {
#if defined(_WIN32) && defined(_MSC_VER)
    //
    // assume Microsoft's runtime library too. NOTE: It is VITAL
    // to use the runtime library, and not Win32, to create the thread
    // as it can (and do) use e.g. stdio that needs to be initialized
    // for every thread using it.
    //
    tp->threads[i]=(HANDLE)_beginthreadex(
			      NULL,  // default security attributes 
			      0,     // 0 => same stack size as main thread
			      &threadFuncWin32, // threaded fn
			      tp,    // args for threaded fn
			      0,     // 0 => thread created running
			      NULL   // thrdaddr, NULL => unused
			      );
    //
#else
    pthread_create(&(tp->threads[i]),NULL,(threadFunc),(void*)tp);
#endif
  }
  return TRUE_TP;
}

void tp_close(threadpool_t *tp) 
{
  //
  // SOMETHING IS VERY FUCKED UP WITH THIS FUNCTION AND I KNOW NOT
  //  WHAT IT IS!!!!
  //
  // Set the shutdown flag.
  //
  tp->Qlock.lock();
  tp->shutdown=1;
  tp->Qlock.unlock();
  //
  // Wake up any sleeping threads
  //
#if !defined(_WIN32) || !defined(_MSC_VER)
  tp->QNotEmpty.broadcast();
#else
  for (int i=0;i<tp->numThread;i++) {
    //
    // It seems that, under win32 we need to broadcast to wake up each
    //  thread. MAYBE WE COULD FIX THAT MY MODIFYING THE PARAMETERS OF
    //  THE threadpool_t::QNotEmpty CONDITIONAL!!!!
    //
    tp->QNotEmpty.broadcast();
  }
#endif
  //
  // Destroy the threads.
  //
  for (int i=0;i<tp->numThread;i++) {
#if defined(_WIN32) && defined(_MSC_VER)
    //
    // NOTE: We cant use WaitForMultipleObjects(.. ALL) due to fact there are
    // more than MAXIMUM_WAIT_OBJECTS (64) worker threads.
    // TODO: Consider refactor this to instead use I/O completion ports
    // or asynch I/O.
    //
//     WaitForSingleObject(HANDLE(tp->threads[i]),INFINITE);
    WaitForSingleObject(HANDLE(tp->threads[i]),30000); // 30 second timeout
    //
    // Threads started with _beginthreadex() require CloseHandle() to be
    //   called.
    //
    CloseHandle(HANDLE(tp->threads[i]));
    tp->threads[i] = 0;
    tp->QNotEmpty.broadcast();
#else
    pthread_join(tp->threads[i],NULL);
#endif
  }
  //
  // Free the threads. SEEMS TO CAUSE THE SEGFAULT!!!! 
  //
#if defined(__WXGTK__)
  DbgMsg(" tp_close(): (gtk) pool: %s - could be leaking some memory here!",
	 tp->dbg_name.c_str());
#elif defined(_WIN32) && defined(_MSC_VER)
  DbgMsg(" tp_close(): (win32) pool: %s - closed %i threads",
	 tp->dbg_name.c_str(),tp->numThread);
#endif
  //
  // Free the thread array.
  //
  delete [] tp->threads;
  return;	
}

int tp_addWork(threadpool_t *tp,void (*func)(void *),void *args,
	       const char *caller,const char *thrd_tag) 
{
  //
  // Adds work to the queue 
  //
  //
  // Create a work struct.
  //
  threadWork_t *work=new threadWork_t;
  //
  // Copy caller name for debugging.
  //
  work->caller=caller;
  if (thrd_tag) work->thread_tag=thrd_tag;
  //
  // Put the function and arguments into work.
  //	
  work->func=func;
  work->arg =args;
  work->next=NULL;
  //
  // Lock the queue.
  //
  tp->Qlock.lock();
  //
  // Add work to the queue
  //
  if (tp->Qhead==NULL) {
    //
    // The queue is empty.
    //
    tp->Qhead=work;
    tp->Qtail=work;
    //
    // the queue was empty before, so there should be at least 1 worker
    // thread blocking on the QNotEmpty condition. Let's wake it up!
    //
    tp->QNotEmpty.broadcast();
  }
  else {
    //
    // The queue is not empty add this job to the tail.
    //
    tp->Qtail->next=work;	  
    tp->Qtail=work;
    //
    // Inform anything that is waiting that the queue isn't empty
    // Note that this might not wake anything up!
    //
    tp->QNotEmpty.broadcast();
  }
  //
  // Unlock the queue and return success
  //
  tp->Qlock.unlock();
  return TRUE_TP;
}

/* This is the body of the worker thread */

void* threadFunc(void *ptr) 
{
  //
  // Set a ptr to the parent threadpool object.
  //
  threadpool_t *tp=(threadpool_t *) ptr;
  //
  // Seed rand() for this thread - it seems rand() needs to be seeded
  //  for each thread (at least inder windows) so we do it here. We don't
  //  want to use the same seed for all threads - just incrementing the
  //  seed should ensure a different sequence for each thread.
  // We cheat here and use Qlock to control access to monsecs_gl
  //  while we increment it. Ugly, but hopefully harmless.
  //
  tp->Qlock.lock();
  srand(monsecs_gl++);
  tp->Qlock.unlock();
  //
  // Loop waiting for work items to be added to the queue.
  //
  for (;;) {	
    //
    // Queue empty? Wait until something is added to it
    //
    tp->Qlock.lock();
    if (tp->Qhead==NULL) {
#if defined(_WIN32)
      tp->Qlock.unlock();
      tp->QNotEmpty.wait();
      tp->Qlock.lock();
#else
      tp->QNotEmpty.wait(tp->Qlock);
#endif
    }
    //
    // Is the TP being closed? If so end this thread 
    //
    if (tp->shutdown) {
      tp->Qlock.unlock();
#if defined(_WIN32) && defined(_MSC_VER)
      return 0; // The CRT handles the rest.
#else
      pthread_exit(NULL);
#endif
    }
    //	 
    // Get work from head of queue and execute it's function 
    //
    if (tp->Qhead) {
      //
      // Save the caller name and increment counter - for debugging.
      //
      threadWork_t tw;
      tw.caller=tp->Qhead->caller.c_str();
// #ifndef RELEASE_ON
//       if (tp->dbg_name=="Main Threadpool") {
// 	dbg_incr_nthr(1);
// 	dbg_print_nthr(" threadFunc(): starting job",tw.caller.c_str());
//       }
// #endif
      //
      // Copy work item info from the queue head into tw and
      //  delete the queue head item. 
      //
      tw.func  =tp->Qhead->func;
      tw.arg   =tp->Qhead->arg;
#if defined(_WIN32) && defined(_MSC_VER)
      if (!tp->Qhead->thread_tag.empty()) {
	set_msvc_debug_info(tp->Qhead->thread_tag,false);
      }
#endif
      threadWork_t *tmp=tp->Qhead; // used to hold the head ptr for delete.
      tp->Qhead=tp->Qhead->next;
      if (tp->Qhead==NULL) tp->Qtail=NULL;
      delete tmp;
      //
      // Execute the work item. We release the queue lock while
      //  we do this and lock it again when we are ready to inform
      //  the queue that we are done
      //
      tp->Qlock.unlock();
      (*tw.func)((void *)tw.arg);			
      tp->Qlock.lock();
      //
      // If the threadpool was closed while we were executing the work item,
      //  we will call QNotEmpty.wait() above and may deadlock. This
      //  seems to be a particular problem under linux where pthread
      //  broadcast has no effect if no threads are waiting on the 
      //  conditional. Hence we check the shutdown flag again here.
      //
      if (tp->shutdown) {
        tp->Qlock.unlock();
#if defined(_WIN32) && defined(_MSC_VER)
        return 0; // The CRT handles the rest.
#else
        pthread_exit(NULL);
#endif
      }
      //
      // Decrement debug counter.
      //
// #ifndef RELEASE_ON
//       if (tp->dbg_name=="Main Threadpool") {
// 	dbg_incr_nthr(-1);
// 	dbg_print_nthr(" threadFunc(): finishing job",tw.caller.c_str());
//       }
// #endif
      //
      // Inform tp_waitOnFinished that another thread has finished 
      //
      tp->numfinished++;
      tp->TFinished.broadcast();
    }
    tp->Qlock.unlock();
  }
}


void tp_reset(threadpool_t *tp) 
{
  //
  // Not currently used.
  //
  tp->Qlock.lock();
  tp->numfinished=0;
  tp->Qlock.unlock();
}


int tp_waitOnFinished(threadpool_t *tp,int n) 
{
  //
  // Not currently used.
  //
  // Wait on all current tasks to complete.
  //
  //std::cout <<"waiting on threads... "<<tp->numfinished<<", out of "<<n<<"\n";
  //std::cout.flush();
  tp->Qlock.lock();
  //std::cout <<"waiting on threads part 2 ... "<<tp->numfinished<<", out of "<<n<<"\n";
  //std::cout.flush();
  while(tp->numfinished<n) {
#if defined(_WIN32)
    tp->Qlock.unlock();
    tp->TFinished.wait();
    tp->Qlock.lock();
#else
    //std::cout <<"waiting on threads part 3 ... "<<tp->numfinished<<", out of "<<n<<"\n";
    //std::cout.flush();
    tp->TFinished.wait(tp->Qlock);
    //std::cout <<"waiting on threads part 4 ... "<<tp->numfinished<<", out of "<<n<<"\n";
    //std::cout.flush();
#endif
  }
  tp->numfinished=0;
  tp->Qlock.unlock();
  return n;
}

//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//     Cygwin threadpool - a bodge to get around the nasty mem leak        //
//      If the leak is fixed for we stop using cygwin for dev, then        //
//      these funcs can be dispensed with.                                 //
//      For non cygwin builds, this will just wrap the normal              //
//      threadpool.                                                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//
#include "../reefa_constants.h"
#include "../logmessages.h"


#if !defined(_CYGWIN)

int CygTP_Init(cygthreadpool_t *ctp,int nt,long int wait_time,
	       const char *name,int priority,bool name_thread) 
{
  //
  // Non-cygwin build - wrap the normal threapool init function.
  // 
  tp_init(ctp,nt,name,priority,name_thread);
  return TRUE_TP;
}

int CygTP_Close(cygthreadpool_t *ctp) 
{
  //
  // Non-cygwin build - DO NOTHING UNTIL WE FIX THE TP_CLOSE CRASH!!!! 
  //
#if defined(__WXGTK__)
  tp_close(ctp); // LEAKS SOME MEMORY!!!!
#elif defined(_WIN32) && defined(_MSC_VER)
  tp_close(ctp); // MAYBE LEAKS SOME MEMORY!!!! :)
#else
  DbgMsg(" Cygpool (%s) CygTP_Close(): (NON-CYGWIN BUILD - CLOSE NOT IMPLEMENTED!!!!)",
	 ctp->dbg_name.c_str());
#endif
  return TRUE_TP;
}

int CygTP_addWork(cygthreadpool_t *ctp,void (*func)(void *),void *args,
		  const char *caller)
{
  //
  // Non-cygwin build - wrap the normal threapool addWork function 
  // 
  tp_addWork(ctp,func,args,caller);
  return TRUE_TP;
}
#else
//
// Cygwin build - use the cygwin threadpool bodge.
//

extern threadpool_t tp;


int CygTP_Init(cygthreadpool_t *ctp,int nt,long int wait_time,
	       const char *name,
	       int priority,bool name_thread) 
{
  //
  // Initialization function - creates the threads and sets up the 
  //  tpool struct. Priority is not supported in this cygpool hack!
  //
  ctp->numThreads=nt;
  ctp->wait_time=wait_time;
  ctp->Qhead=ctp->Qtail=NULL;
  ctp->numfinished=0;
  ctp->shutdown=0;
  ctp->dbg_name=name;
  //
  // 
  //
  for (int i=0;i<nt;i++) {
    tp_addWork(&tp,CygTP_ThreadW,(void *) ctp,"CygTP_Thread");
  }
  return TRUE_TP;
}

int CygTP_Close(cygthreadpool_t *tp) 
{
  //
  // Wait for the threads to finish.
  //
  CygTP_SetShutdown(tp);
  int t=0;
  int n=CygTP_NumFinished(tp);
  while (n<tp->numThreads) {
    DbgMsg("\twaiting for %i cygpool threads (%s) to terminate:  %i/30 secs",
	   tp->numThreads-n,tp->dbg_name.c_str(),t);
    sleep(SLEEP_SEC);
    if (t++>30) {
      DbgMsg("\t\tlooks like we couldn't terminate some threads cleanly.");
      break;
    }
    n=CygTP_NumFinished(tp);
  }
  if (t<=30) DbgMsg(" Cygpool (%s): %i threads  terminated.",
		    tp->dbg_name.c_str(),CygTP_NumFinished(tp));
  return TRUE_TP;
}

void CygTP_SetShutdown(cygthreadpool_t *tp)
{
  tp->shutdownLock.lock();
  tp->shutdown=1;
  tp->shutdownLock.unlock();
}
int CygTP_Shutdown(cygthreadpool_t *tp)
{
  tp->shutdownLock.lock();
  int val=tp->shutdown;
  tp->shutdownLock.unlock();
  return val;
}
void CygTP_IncrNumFinished(cygthreadpool_t *tp,int incr)
{
  tp->numfinishedLock.lock();
  tp->numfinished+=incr;
  tp->numfinishedLock.unlock();
}
int CygTP_NumFinished(cygthreadpool_t *tp)
{
  tp->numfinishedLock.lock();
  int val=tp->numfinished;
  tp->numfinishedLock.unlock();
  return val;
}


void CygTP_ThreadW(void *arg)
{
  CygTP_Thread((cygthreadpool_t *) arg);
}
void CygTP_Thread(cygthreadpool_t *tp)
{
  //
  // Wrapper for server() used in threading. 
  //
  threadWork_t *tw;
  while ((tw=CygTP_GetTask(tp))) {
    CygTP_ExecuteTask(tw,tp->dbg_name.c_str());
    delete tw;
    if (CygTP_Shutdown(tp)==1) break;
  }
  CygTP_IncrNumFinished(tp,1);
}

threadWork_t *CygTP_CheckForTask(cygthreadpool_t *tp)
{
  //
  // If the queue has members return, otherwise return 0.
  //
  threadWork_t *tw=NULL;  
  tp->Qlock.lock();
  if (tp->Qhead!=NULL) {
    tw=tp->Qhead;
    tp->Qhead=tw->next;
    if (tp->Qhead==NULL) tp->Qtail=NULL; // list is now empty
  }
  tp->Qlock.unlock();
  return tw;
}

threadWork_t *CygTP_GetTask(cygthreadpool_t *tp)
{
  //
  // When a task appears in tp remove it from the queue
  //  and return a ptr to it.
  //
  threadWork_t *tw=NULL;
  while (!(tw=CygTP_CheckForTask(tp))) {
    if (CygTP_Shutdown(tp)==1) return NULL;
    usleep(tp->wait_time);
  }
  return tw;
}


void CygTP_ExecuteTask(threadWork_t *tw,const char *dbg_name0)
{
  //
  // Execute the function in tw.
  //
  (*tw->func)((void *)tw->arg);
}


int CygTP_addWork(cygthreadpool_t *tp,void (*func)(void *),void *args,
		  const char *caller) 
{
  //
  // Adds work to the queue 
  //
  //
  // Create a work struct.
  //
  threadWork_t *work;
  work=new threadWork_t;
  work->caller=caller;            // Copy caller name for debugging.
  work->func=func;                // Put the function and arguments into work.
  work->arg =args;
  work->next=NULL;
  //
  // Add work to the queue
  //
  tp->Qlock.lock();
  if (tp->Qhead==NULL) {
    //
    // The queue is empty.
    //
    tp->Qhead=work;
    tp->Qtail=work;
  }
  else {
    //
    // The queue is not empty add this job to the tail.
    //
    tp->Qtail->next=work;	  
    tp->Qtail=work;
  }
  tp->Qlock.unlock();
  return TRUE_TP;
}

#endif


//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//      Wrapper for a conditional - we can use this to wait and            //
//  broadcast a wake up signal.                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//


off_thread_cond_t *ThreadCond_Create(const char *name,const long maxcount)
{
  //
  // Create and initialise a conditional - returns a ptr to it.
  //
  off_thread_cond_t *tc=new off_thread_cond_t;
  tc->dbg_name=name;
  tc->cond.init(maxcount);
  return tc;
}

void ThreadCond_Destroy(off_thread_cond_t *tc)
{
  //
  // "Destructor" function - just wraps delete().
  //
  delete tc;
}

void ThreadCond_Wait(off_thread_cond_t *tc)
{
  //
  // Wait on a wakeup signal.
  //
#if defined(_WIN32)
  tc->cond.wait();
#else
  tc->lock.lock();
  tc->cond.wait(tc->lock);
  tc->lock.unlock();
#endif
}

void ThreadCond_TimedWait(off_thread_cond_t *tc,int secs)
{
  //
  // Wait on a wakeup signal.
  //
#if defined(_WIN32)
  tc->cond.timed_wait(secs);
#else
  tc->lock.lock();
  tc->cond.timed_wait(tc->lock,secs);
  tc->lock.unlock();
#endif
}

void ThreadCond_Broadcast(off_thread_cond_t *tc)
{
  //
  // Send a wakeup signal to any threads that are waiting on cond.
  //
  tc->lock.lock();
  tc->cond.broadcast();
  tc->lock.unlock();
}
