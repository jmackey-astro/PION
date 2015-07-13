//
//    Functions which handle messaging between the core and the gui
//     or logging to files.
//
#ifndef R_LOGMESSAGES_H
#define R_LOGMESSAGES_H

#include "reefa_constants.h"

void Abort(const char *fmt, ...);
void DbgMsg(const char *fmt, ...);

#endif // #ifndef R_LOGMESSAGES_H
