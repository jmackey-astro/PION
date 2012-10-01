//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//            Functions for writing msgs to logs and stderr.               //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include "constants.h"
#include "logmessages.h"
//
#include <string>
#include <sstream>
//
using namespace std;
//
// File scope variables.
//
//
// Function prototypes
//
static void ShowMessage(const char *str,int iretbef,int iretaft,
			   bool no_post=false);
void PrintStrToStderr(const char *str);
//
// External variables.
//
//
// External functions.
//
//
// Functions
//

#if defined(_WIN32)
static void flockfile(FILE* f) {
/* We _could_ do something like
  LockFile(_get_osfhandle(...), ...) */
}
static void funlockfile(FILE* f) {
/* We _could_ do something like
  UnlockFile(_get_osfhandle(_get_osfhandle(...), ...) */
}
#endif // _WIN32


void strip_white_space(string &str)
{
  //
  // Remove any white space from the start and the end of the string.   
  //
  //   i) trim from right.
  //
  string::size_type i=str.find_last_not_of(" \t\r\n");
  if (i!=string::npos) str.erase(i+1);
  //
  //  ii) trim from left
  //
  str.erase(0,str.find_first_not_of(" \t\r\n"));
}

string dbg_ascii_time(time_t at)
{
  //
  // Wraps asctime(). Uses LOCALTIME(), which is gmtime() for debug builds.
  //
#if defined(_MSC_VER) && _MSC_VER > 1200
  char buf[128]={"<error getting time>"};
  struct tm _tm;
 #ifndef RELEASE_ON
  gmtime_s(&_tm,&at);
 #else
  localtime_s(&_tm,&at);
 #endif
  asctime_s(buf,128,&_tm);
  string label=buf;
#else
  string label=asctime(LOCALTIME(&at));
#endif
  strip_white_space(label);
  return label;
}

string time_prefix()
{
  //
  // Returns a string containing the time stamp at the start of a 
  //  logfile line.
  //
  time_t now=time(&now);
  string tmp=dbg_ascii_time(now);
  tmp+=": ";
  return tmp;
}


void add_line_break(char *str)
{
  //
  // Append a '\n' to str if necessary.
  //
  int l=strlen(str)-1;
  if (l<0) {
    str[0]='\n';
    str[1]='\0';
    return;
  }
  if (str[l]!='\n') {
    str[l+1]='\n';
    str[l+2]='\0';
  }
}


void DbgMsg(const char *fmt, ...)
{
  //
  // Writes to stderr or a log file.
  //
  va_list args;
  va_start(args,fmt);
  char str[VSN_STR_SIZE];
  VSNPRINTF(str,VSN_STR_ARG,fmt,args);
  add_line_break(str);
  PrintStrToStderr(str);
  va_end(args);
}


static void ShowMessage(const char *str,int iretbef,int iretaft,
			   bool no_post)
{
  //
  // Sends the message, str[], to the GUI-log or the screen.
  //  Prints iretbef blank lines before the message and iretaft blank
  //  lines after the message. 
  //
#ifdef GUI_ON
 #ifndef _WIN32
  DbgMsg("%s",str);
 #endif
#else
  for (int i=0;i<iretaft;i++) PrintStrToStderr("\n");
  DbgMsg("%s",str);
  for (int i=0;i<iretaft;i++) PrintStrToStderr("\n");
#endif
}


void PrintStrToStderr(const char *str)
{
#ifndef _WIN32
  flockfile(stderr);
  fprintf(stderr,"%s%s",time_prefix().c_str(),str);
  funlockfile(stderr);
#endif
}


void Abort(const char *fmt, ...)
{
  //
  // Handles fatal errors - writes any message in the call to
  //  a file and terminates the program.
  //
  va_list args;
  va_start(args,fmt);
  char err_msg[VSN_STR_SIZE];
  VSNPRINTF(err_msg,VSN_STR_ARG,fmt,args);
  time_t now=time(&now);
  string time_str=dbg_ascii_time(now);
  DbgMsg("\n\n  Terminal error at: %s\n\n",time_str.c_str());
  PrintStrToStderr(err_msg);
  DbgMsg("\n\n\n\n");
  va_end(args);
  exit(1);
}


