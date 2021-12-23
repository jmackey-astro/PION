/// \file raytracer_SC_pllel.cc
///
/// \author Jonathan Mackey
///
/// Definitions for parallel uniform grid.
///
/// Modifications:\n
///  - 2010-01-19 JM: call new weighting scheme functions instead of duplicating
///  serial code.
///  - 2010-01-22 JM: worked on sources at cell-corners to try and get it
///  working.
///  - 2010-01-23 JM: made sure source cell geometry is only for cell-centred
///  sources.
///  - 2010-01-26 JM: Debugging sources at cell-vertices in 3D.
///  - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
///    Changed source indexing so IDs are set by SimPM.RS and not assigned
///    locally. Changed Add_Source() interface so that all info gets passed to
///    class. Should handle multiple sources now, and sources at infinity.  Of
///    course I need a "rates-only" update function to actually use multiple
///    emitting sources.
/// - 2011.03.21 JM: Added RayTrace_Column_Density() interface function.  Just a
///    dummy function for now.  (2011.04.15 JM: it points to the old update fn
///    now! The deciding of whether to update or calculate column-density is in
///    ProcessCell()).
/// - 2011.04.22 JM: Updated Add_Source() to set Vshell in cells if
/// Column-density
///    update is required.
/// - 2011.04.23 JM: Added TauMin[] array of values of TauMin for each source,
/// ordered
///    by source id.  This removes the need for interpolate_2D_RHO() so I got
///    rid of it. Now the interpolate functions read TauMin[src_id] and apply
///    that to the interpolation.
/// - 2012.02.27 JM: Only output timings for source with id=0 (reduces
/// filesize).
/// - 2012.03.31 JM: Added separators between functions.  Updated Add_Source()
///    function to call add_source_to_list() and set_Vshell_for_source() at the
///    appropriate places.
/// - 2013.09.05 JM: Debugged for new get/set col functions.
/// - 2015.01.28 JM: New include statements for new file structure.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/mem_manage.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#ifndef NDEBUG
#endif  // NDEBUG

#include "raytracing/raytracer_SC_pllel.h"
#include <fstream>
using namespace std;

#ifdef PARALLEL

// ##################################################################
// ##################################################################

raytracer_USC_pllel::raytracer_USC_pllel(
    class GridBaseClass *ggg,      ///< Pointer to grid
    class microphysics_base *mmm,  ///< Pointer to MicroPhysics Class.
    class SimParams *sp,           ///< simulation parameters
    class Sub_domain *mcmd,        ///< domain decomposition info
    int nd,                        ///< number of dimensions of grid
    int csys,                      ///< coordinate system
    int nv,                        ///< number of variables in state vector
    int ftr,      ///< index of first tracer variable in state vector
    int Nsources  ///< Number of radiation sources
    ) :
    raytracer_USC(ggg, mmm, nd, csys, nv, ftr, Nsources)
{
#ifdef RT_TESTING
  spdlog::info("SC PARALLEL raytracer class constructor!");
#endif

  par        = sp;
  sub_domain = mcmd;
  return;
}

raytracer_USC_pllel::~raytracer_USC_pllel()
{
#ifdef RT_TESTING
  spdlog::info("SC PARALLEL raytracer class destructor!");
#endif
}

// ##################################################################
// ##################################################################

int raytracer_USC_pllel::Add_Source(struct rad_src_info *src  ///< source info.
)
{
  //
  // First call serial version.  This finds the source, and centres
  // it on a cell vertex if needed.
  //
#ifdef RT_TESTING
  spdlog::info("\n--BEGIN-----raytracer_USC_pllel::AddSource()------------\n"
               "\t**** PARALLEL Add_Source: calling serial version.\n");
#endif
  if (src->at_infinity) {
    raytracer_USC_infinity::add_source_to_list(src);
  }
  else {
    add_source_to_list(src);
  }
  int id = SourceList.back().s->id;
  spdlog::debug("\t**** PARALLEL Add_Source: serial fn ret id={}", id);

  //
  // Now tell the parallel grid to decide which boundaries it needs
  // to receive data from before processing this source, and which it
  // needs to send data to after processing.
  //
#ifdef RT_TESTING
  spdlog::info("\t**** PARALLEL Add_Source: Setup extra RT boundaries");
#endif
  int err = Setup_RT_Boundaries(*par, *sub_domain, gridptr, id, *src);
  if (err) spdlog::error("{}: {}", "Failed to setup RT Boundaries", err);

    //
    // Set Vshell for every cell on the grid.
    //
#ifdef RT_TESTING
  spdlog::info("\t**** PARALLEL: setting Vshell for source");
#endif
  set_Vshell_for_source(&SourceList.back());

#ifdef RT_TESTING
  spdlog::info("\t**** PARALLEL Add_Source: all done..\n"
               "--END-----raytracer_USC_pllel::AddSource()------------\n");
#endif
  return id;
}

// ##################################################################
// ##################################################################

int raytracer_USC_pllel::RayTrace_SingleSource(
    const int s_id,   ///< Source id
    const double dt,  ///< Timestep
    const double g    ///< eos gamma.
)
{
  int err = 0;
#ifdef RT_TESTING
  spdlog::debug("RT_MPI: Starting Raytracing for source: {}", s_id);
#endif

  // Find source in list.
  struct rad_src_info *RS = 0;
  for (vector<rad_source>::iterator i = SourceList.begin();
       i != SourceList.end(); ++i) {
    if ((*i).s->id == s_id) RS = (*i).s;
  }
  if (!RS) {
    spdlog::error(
        "{}: {}", "RayTrace_SingleSource() source not in source list.", s_id);
  }

  string t1 = "totalRT", t2 = "waitingRT", t3 = "doingRT";  //, t4="tempRT";
  // double total = 0.0, wait = 0.0, run = 0.0;

  clk.start_timer(t1);
  //
  // First Receive RT boundaries from processors nearer source.
  //
  clk.start_timer(t2);
  // clk.start_timer(t4);
#ifdef RT_TESTING
  spdlog::info("RT_MPI: receiving RT boundaries");
#endif
  err += Receive_RT_Boundaries(*par, *sub_domain, gridptr, s_id, *RS);

#ifdef RT_TESTING
  spdlog::info("RT_MPI: received RT boundaries");
#endif

  clk.pause_timer(t2);

  //
  // Now we have the boundary conditions, so call the serial Raytracer.
  //
  clk.start_timer(t3);
  // clk.start_timer(t4);
#ifdef RT_TESTING
  spdlog::info("RT_MPI: calling serial raytrace function");
#endif

  err += raytracer_USC::RayTrace_SingleSource(s_id, dt, g);

#ifdef RT_TESTING
  spdlog::info("RT_MPI: serial raytrace done");
#endif

  // run = clk.pause_timer(t3);

  //
  // Finally, send the new column densities to processors further from source.
  //
  clk.start_timer(t2);
  // clk.start_timer(t4);
#ifdef RT_TESTING
  spdlog::info("RT_MPI: sending RT boundaries");
#endif

  err += Send_RT_Boundaries(*par, *sub_domain, gridptr, s_id, *RS);
#ifdef RT_TESTING
  spdlog::info("RT_MPI: sent RT boundaries");
#endif

  // wait  = clk.pause_timer(t2);
  // total = clk.pause_timer(t1);

  return err;
}

// ##################################################################
// ##################################################################

#endif  // PARALLEL
