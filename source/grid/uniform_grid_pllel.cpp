/// \file uniform_grid_pllel.cpp
/// \author Jonathan Mackey
///
/// Definitions for parallel uniform grid.
///
/// Modifications:\n
/// - 2007-11-06 started to get it going.
/// - 2010-01-22 JM: Changed RT source finding if we have non-cell-centred
/// sources.
/// - 2010-01-26 JM: Debugging 3D RT for vertex-centred sources.
/// - 2010-02-03 JM: changed variable names in destructor ('i' was defined
/// twice!)
/// - 2010-03-13 JM: moved BoundaryTypes enum to uniformgrid.h; Added
/// oneway-outflow BC
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.  Switched 'endl' statements for '\n'.
/// - 2010.11.19 JM: Got rid of testing myalloc() myfree() functions.
/// - 2010.12.07 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef.  Will probably keep it since it is the way
///   things will go eventually.  The new grid classes have extra
///   function for the distance between two points, two cells, and
///   between a vertex and a cell.
/// - 2011.02.24 JM: Worked on generalising multi-source radiative transfer.
///    Still some way to go to getting it functional.  Wraps old setup()
///    function for point sources into a new function.  The setup function now
///    calls either this new function, or another new one for sources at
///    infinity.
/// - 2011.02.25 JM: removed HCORR ifdef around new code.
///     Added possibility of multiple sources in send/recv.
/// - 2011.02.28 JM: Got rid of RSP radiation-sources-parameters class.
/// - 2011.03.02 JM: Fixed bugs in diffuse radiation boundary comms.
/// - 2011.03.21 JM: Got rid of zero-ing of column densities.  It is done in
///    the cell constructor, so we don't need to worry about it.
/// - 2011.03.22 JM: small bugfixes.
/// - 2011.04.22 JM: Fixed setup_recv_boundaries() functions so that new
/// boundary cells
///    are only created if they don't already exist.  It now works for multiple
///    point sources, and I don't think they need to be at the same place.
/// - 2013.09.05 JM: Debugged for new get/set col functions.
/// - 2015.01.28 JM: Removed GEOMETRIC_GRID (it is default now), and
///    updated for the new code structure.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
/// - 2016.03.21 JM: Worked on simplifying RT boundaries (first hack!
///    compiles but buggy...) 03.24:fixed some bugs, redefined dir_XP
/// - 2017.11.07-22 JM: updating boundary setup.
/// - 2018.08.09 JM: moved all boundary stuff to new classes.
/// - 2018.10.12 JM: added flux_update functions for BC89 refinement.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include "grid/uniform_grid_pllel.h"
#include "microphysics/microphysics_base.h"
#include <fstream>
#include <sstream>
using namespace std;

#ifdef PARALLEL

//#define TEST_BC89FLUX

//-------------------------------------------------------------
//------------------- CARTESIAN GRID START --------------------
//-------------------------------------------------------------

// ##################################################################
// ##################################################################

UniformGridParallel::UniformGridParallel(
    int nd,
    int nv,
    int eqt,
    int nbc,  ///< number of boundary cells to use.
    const std::array<double, MAX_DIM>
        &xn,  ///< array of minimum values of x,y,z for this grid.
    const std::array<double, MAX_DIM>
        &xp,  ///< array of maximum values of x,y,z for this grid.
    const std::array<int, MAX_DIM>
        &nc,  ///< array of number of cells in x,y,z directions.
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM>
        &sim_xn,  ///< array of min. x/y/z for full simulation.
    std::array<double, MAX_DIM>
        &sim_xp  ///< array of max. x/y/z for full simulation.
    ) :
    VectorOps_Cart(nd),
    UniformGrid(nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp)
{

#ifndef NDEBUG
  spdlog::info("UniformGridParallel constructor");
  spdlog::debug("Local Xmin : {}", xn);
  spdlog::debug("Local Xmax : {}", xp);
  spdlog::debug("Local Npt  : {}", nc);

  spdlog::debug("SIM Xmin  : {}", Sim_xmin);
  spdlog::debug("SIM Xmax  : {}", Sim_xmax);
  spdlog::debug("SIM Range : {}", Sim_range);

  spdlog::debug("SIM iXmin  : {}", Sim_ixmin);
  spdlog::debug("SIM iXmax  : {}", Sim_ixmax);
  spdlog::debug("SIM iRange : {}", Sim_irange);

  spdlog::info("UniformGridParallel constructor done");
#endif

  return;
}

//-------------------------------------------------------------
//------------------- CARTESIAN GRID END ----------------------
//-------------------------------------------------------------

//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID START ------------------
//-------------------------------------------------------------

// ##################################################################
// ##################################################################

///
/// Constructor
///
uniform_grid_cyl_parallel::uniform_grid_cyl_parallel(
    int nd,
    int nv,
    int eqt,
    int nbc,  ///< number of boundary cells to use.
    const std::array<double, MAX_DIM> &xn,
    const std::array<double, MAX_DIM> &xp,
    const std::array<int, MAX_DIM> &nc,
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM>
        &sim_xn,  ///< array of min. x/y/z for full simulation.
    std::array<double, MAX_DIM>
        &sim_xp  ///< array of max. x/y/z for full simulation.
    ) :
    VectorOps_Cart(nd),
    UniformGrid(nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp),
    UniformGridParallel(
        nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp),
    VectorOps_Cyl(nd),
    uniform_grid_cyl(
        nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp)
{
  spdlog::debug(
      "uniform_grid_cyl_parallel:: cyl. uniform PARALLEL grid\n G_ndim={} and G_nvar={}",
      G_ndim, G_nvar);
  if (G_ndim > 2)
    spdlog::error("{}: {}", "Need to write code for 3 dimensions", G_ndim);

  spdlog::debug("cylindrical grid: dR={}", G_dx);
  return;
}

// ##################################################################
// ##################################################################

uniform_grid_cyl_parallel::~uniform_grid_cyl_parallel()
{
  spdlog::info("uniform_grid_cyl_parallel destructor. Present and correct!");
}

// ##################################################################
// ##################################################################

double uniform_grid_cyl_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  return (R_com(c, G_dx) - Sim_xmin[Rcyl]) / CI.phys_per_int()
         + Sim_ixmin[Rcyl];
}

// ##################################################################
// ##################################################################

//-------------------------------------------------------------
//------------------- CYLINDRICAL GRID END --------------------
//-------------------------------------------------------------

//-------------------------------------------------------------
//------------------- SPHERICAL GRID START --------------------
//-------------------------------------------------------------

// ##################################################################
// ##################################################################

///
/// Constructor
///
uniform_grid_sph_parallel::uniform_grid_sph_parallel(
    int nd,
    int nv,
    int eqt,
    int nbc,  ///< number of boundary cells to use.
    const std::array<double, MAX_DIM> &xn,
    const std::array<double, MAX_DIM> &xp,
    const std::array<int, MAX_DIM> &nc,
    std::array<double, MAX_DIM> &lev_xn,  // level xmin
    std::array<double, MAX_DIM> &lev_xp,  // level xmax
    std::array<double, MAX_DIM>
        &sim_xn,  ///< array of min. x/y/z for full simulation.
    std::array<double, MAX_DIM>
        &sim_xp  ///< array of max. x/y/z for full simulation.
    ) :
    VectorOps_Cart(nd),
    UniformGrid(nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp),
    UniformGridParallel(
        nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp),
    VectorOps_Cyl(nd), VectorOps_Sph(nd),
    uniform_grid_sph(
        nd, nv, eqt, nbc, xn, xp, nc, lev_xn, lev_xp, sim_xn, sim_xp)
{
  spdlog::debug(
      "uniform_grid_sph_parallel:: sph. uniform PARALLEL grid G_ndim={} and G_nvar={}",
      G_ndim, G_nvar);
  if (G_ndim != 1)
    spdlog::error("{}: {}", "Need to write code for >1 dimension", G_ndim);

  spdlog::debug("spherical grid: dr={}", G_dx);
  return;
}

// ##################################################################
// ##################################################################

uniform_grid_sph_parallel::~uniform_grid_sph_parallel()
{
  spdlog::info("uniform_grid_sph_parallel destructor. Present and correct!");
}

// ##################################################################
// ##################################################################

double uniform_grid_sph_parallel::iR_cov(const cell *c)
{
  //
  // integer and physical units have different origins, so I need a
  // function to get R_com() in integer units, measured from the
  // integer coord-sys. origin.
  //
  return (R_com(c, G_dx) - Sim_xmin[Rsph]) / CI.phys_per_int()
         + Sim_ixmin[Rsph];
}

// ##################################################################
// ##################################################################

//-------------------------------------------------------------
//-------------------- SPHERICAL GRID END ---------------------
//-------------------------------------------------------------
#endif  // PARALLEL
