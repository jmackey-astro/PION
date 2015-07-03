// Thread and locking support types and functions.
#if !defined(R_THREADS_H)
#define R_THREADS_H

#if defined(_CYGWIN)
#  undef _WIN32
#endif

#if defined(_WIN32)
# if defined _MSC_VER
#  pragma once
# endif
# if defined(_WIN32_WINNT) && (_WIN32_WINNT < 0x0400)
#  error Dude... _WIN32_WINNT < 0x0400? You have got to be kidding!
# endif
# if !defined(_WIN32_WINNT)
#  define _WIN32_WINNT 0x0400 // intended for Win9x+ and NT4+
# endif
# include <windows.h>
#else
# include <pthread.h>
#endif

#include "constants.h"

#if defined(_WIN32) && !defined(_CYGWIN)
// Please note that the reason to conditionally use mutex vs.
// critical section is that we need to attempt to get a lock, but
// not stall in case it's already taken. As TryEnterCriticalSection
// is only for later NT-based operating system it would prevent
// running the resulting binary on Win9x, where we by conditional
// compilation support them using a (much slower) mutex.
#define R_WIN32_USE_CRITICAL_SECTION
#endif


struct thread_mutex_t
{
  thread_mutex_t() throw();
  ~thread_mutex_t() throw();

  void lock() throw();
  void unlock() throw();

  // Attempts to get the lock. Always returns immediately.
  // Returns true if it got the lock, else false.
  bool try_lock() throw();

#if defined(_WIN32)
# if defined(R_WIN32_USE_CRITICAL_SECTION)
  CRITICAL_SECTION critsect;
# else
  HANDLE mutex;
# endif
#else
  pthread_mutex_t mutex;
#endif

private:
  thread_mutex_t(const thread_mutex_t&); // no impl
  void operator=(const thread_mutex_t&); // no impl

};


struct thread_cond_t
{
  void init(long maxcount=NUM_THREADS_MAIN) throw();
  void broadcast() throw();
#if defined(_WIN32)
  void wait() throw();
  void timed_wait(int secs) throw();
#else
  void wait(thread_mutex_t& mutex) throw();
  void timed_wait(thread_mutex_t& mutex, int secs) throw();
#endif

#if defined(_WIN32)
  HANDLE cond;
#else
  pthread_cond_t cond;
#endif
};


#if defined(_WIN32)
typedef HANDLE pthread_t;

static inline void usleep(unsigned int us)
{
  Sleep(us / 1000);
}
#endif


#endif // R_THREADS_H
