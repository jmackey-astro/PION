/// \file flux_base.cc
/// \author Jonathan Mackey
/// 
/// Contains class declarations for the basic flux-solver class.
/// This is an interface class, and also defines some basic
/// functionality, such as how to deal with tracer variables tacked
/// onto the end of the state vectors.
///
/// Modifications:
/// - 2010.12.23 JM: Moved from eqns_base.h into its own file.
///   Added tracer_flux() function.
/// - 2010.12.27 JM: Enclosed Lapidus AV in an ifdef.
/// - 2011.02.25 JM: removed HCORR ifdef around new code; it is solid now.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING



#include "flux_calc/flux_base.h"
using namespace std;

// **********************************
// ***** BASE FLUX SOLVER CLASS *****
// **********************************


// ##################################################################
// ##################################################################




