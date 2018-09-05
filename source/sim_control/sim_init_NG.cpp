/// \file sim_init_NG.cpp
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for sim_init_NG, which sets up a PION simulation
/// and gets everything ready to run.
///
/// Modifications:\n
/// - 2018.05.11 JM: moved code from sim_control.cpp
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "sim_control/sim_init_NG.h"

#include "raytracing/raytracer_SC.h"
#include "microphysics/microphysics_base.h"
#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"


#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo_NG.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include <climits>
using namespace std;



// ##################################################################
// ##################################################################


sim_init_NG::sim_init_NG()
{
#ifdef TESTING
  cout << "(sim_init_NG::Constructor)\n";
#endif
  return;
}



// ##################################################################
// ##################################################################


sim_init_NG::~sim_init_NG()
{
#ifdef TESTING
  cout << "(sim_init_NG::Destructor)\n";
#endif
  return;
}

