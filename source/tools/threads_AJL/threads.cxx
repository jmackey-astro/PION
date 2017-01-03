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
// Thread and locking support - platform abstraction
#include <assert.h>
#include "threads.h"
//#include <iostream>

thread_mutex_t::thread_mutex_t() throw()
#if defined(_WIN32) &&!defined(R_WIN32_USE_CRITICAL_SECTION)
  : mutex(CreateMutex(0, FALSE, 0))
#endif
{
#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  InitializeCriticalSection(&critsect);
# else
  assert(mutex);
# endif
#else // !_WIN32
  pthread_mutex_init(&mutex, 0);
#endif
}


thread_mutex_t::~thread_mutex_t() throw()
{
#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  DeleteCriticalSection(&critsect);
# else
  CloseHandle(mutex);
# endif
#else // !_WIN32
  pthread_mutex_destroy(&mutex);
#endif
}


void thread_mutex_t::lock() throw()
{
#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  EnterCriticalSection(&critsect);
# else
  WaitForSingleObject(mutex, INFINITE);
# endif
#else
  pthread_mutex_lock(&mutex);
#endif
}


// Attempts to get the lock immediately.
// Returns true if it got the lock.
bool thread_mutex_t::try_lock() throw()
{
#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  return TryEnterCriticalSection(&critsect) != FALSE;
# else
  const DWORD ret = WaitForSingleObject(mutex, 1);
  assert(ret == WAIT_OBJECT_0 || ret == WAIT_ABANDONED || ret == STATUS_TIMEOUT);
  // returns true if we got the lock
  return (ret == WAIT_OBJECT_0) || (ret == WAIT_ABANDONED);
# endif
#else
  return pthread_mutex_trylock(&mutex) == 0; // return true if we got the lock
#endif
}


void thread_mutex_t::unlock() throw()
{
#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  LeaveCriticalSection(&critsect);
# else
  ReleaseMutex(mutex);
# endif
#else
  pthread_mutex_unlock(&mutex);
#endif
}


void thread_cond_t::init(long maxcount) throw()
{
#if defined(_WIN32)
  cond = CreateSemaphore(
			 0,        // No security attributes
			 0,        // Initial count
			 maxcount, // Maximum count.
			 0         // unnamed
			 );
#else
  pthread_cond_init(&cond, 0);
#endif
}


void thread_cond_t::broadcast() throw()
{
#if defined(_WIN32)
  ReleaseSemaphore(      // increases the count by of the semaphore
		   cond, //  increase count of cond
		   1,    //  increase the count by 1
		   NULL  //  we're not interested in the previous value
		   );
#else
  pthread_cond_broadcast(&cond);
#endif
}


#if defined(_WIN32)

void thread_cond_t::wait() throw()
{
  WaitForSingleObject(cond, INFINITE); // waits while the count is zero
                                       // decreases the count on wait end
}

#else

void thread_cond_t::wait(thread_mutex_t& mutex) throw()
{
  //std::cout <<"thread_cond_t::wait() starting \n";
  //std::cout.flush();
  pthread_cond_wait(&cond, &mutex.mutex);
  //std::cout <<"thread_cond_t::wait() finishing \n";
  //std::cout.flush();
}

#endif

#if defined(_WIN32)

void thread_cond_t::timed_wait(int secs) throw()
{
  WaitForSingleObject(cond,secs*1000); // wants wait time in millisecs.
}

#else


void thread_cond_t::timed_wait(thread_mutex_t& mutex, int secs) throw()
{
  time_t now=time(&now);
  struct timespec tv;
  tv.tv_sec =  now+secs;      // this many seconds
  tv.tv_nsec = 0;             // and this many nano-seconds
  pthread_cond_timedwait(&cond, &mutex.mutex,&tv);
}

#endif

