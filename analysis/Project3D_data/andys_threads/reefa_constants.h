//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//            Defined constants used by various modules.                   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//
#ifndef REEFA_CONSTS_H
#define REEFA_CONSTS_H
//
// Check operating system flags.
//
#if !defined(__WXMSW__) && !defined(__WXGTK__) && !defined(__WXMAC__)
  #error bad operating system definition flag (must be __WXMSW__, __WXGTK__ or __WXMAC__.
#endif
//
//  Cygwin specific global defines.
//
#ifdef _CYGWIN
  #include <windows.h> // needed for INT64
  #undef _WIN32 // Cygwin build needs _WIN32 undefined
  #define END_LINE "\n"
#endif
//
//  MSVC specific global defines (msvc).
//
#if defined(_MSC_VER)
  #pragma warning( disable : 4786)   // disable MSVC long name warnings
  #if _MSC_VER <= 1200
    #define for if (0) {} else for   // MSVC has crappy scope rules
  #else
    #define strtoll _strtoi64
  #endif
  #define finite    _finite
  #define vsnprintf _vsnprintf
  #define snprintf  _snprintf
  #define INT64 __int64
  #define END_LINE "\r\n"
  #include <float.h>
#endif
//
//  Windows specific global defines (msvc and cygwin).
//
#ifdef __WXMSW__
  #ifndef END_LINE
    #define END_LINE "\r\n"
  #endif
  #define PATH_DELIM_S "\\"
  #define PATH_DELIM_C '\\'
#endif
//
//  MAC specific global defines.
//
#ifdef __WXMAC__
  #define INT64 long long
  #define UINT32 uint32_t
  #define END_LINE "\n"
  #define PATH_DELIM_S "/"
  #define PATH_DELIM_C '/'
#endif
//
//  Linux specific global defines.
//
#ifdef __WXGTK__
  #define INT64 long long
  #define UINT32 uint32_t
  #define END_LINE "\n"
  #define PATH_DELIM_S "/"
  #define PATH_DELIM_C '/'
#endif
//
// This gets around some VS2005 warnings about vsnprintf, rmdir, etc.
//
#define VSN_STR_SIZE STR_MAX
#if defined(_MSC_VER) && _MSC_VER > 1200
  #define VSNPRINTF vsnprintf_s
  #define VSN_STR_ARG STR_MAX,STR_MAX-10
  #define RMDIR     _rmdir
#else
  #define VSNPRINTF vsnprintf
  #define VSN_STR_ARG STR_MAX-10
  #define RMDIR     rmdir
#endif
//
// Number of threads.
//
#define NUM_THREADS_MAIN  7
#define THREADNAME_MAP_OFF // allows use of getThreadName() under MSVC dbg
//
// Release or debug build.
//
#define RELEASE_ON
//
// Debug versions output timestamps on GMT.
//
#ifndef RELEASE_ON
  #define LOCALTIME gmtime
#else
  #define LOCALTIME localtime
#endif
//
// Maximum string length constants.
//
#define STR_MAX       1024   // max large string length
//
// Time units for convenience - in seconds.
//
#define MINS  (60.0)
#define HOURS (3600.0)
#define DAYS  (86400.0)
#if defined(_WIN32)
#define SLEEP_SEC 1000        // win32 takes milli-seconds
#else
#define SLEEP_SEC 1           // linux takes seconds
#endif
//
// End of line character constants.
//
#define CR 13
#define LF 10

#endif // #ifndef REEFA_CONSTS_H
