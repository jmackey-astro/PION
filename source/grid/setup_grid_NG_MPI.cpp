/// \file setup_grid_NG_MPI.cpp
///
/// \brief Class for setting up parallel NG grids.
///
/// \author Jonathan Mackey
///
///
/// Modifications:
/// - 2018.08.24 JM: derived from setup_NG_grid and
///    setup_fixed_grid_MPI.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "raytracing/raytracer_SC.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "grid/setup_grid_NG_MPI.h"
#include "grid/uniform_grid.h"
#include "grid/uniform_grid_pllel.h"

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include "microphysics/MPv3.h"
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"
#include "microphysics/MPv8.h"
#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"

#ifdef FITS
#include "dataIO/dataio_fits_MPI.h"
#endif  // if FITS

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/MPv9.h"
#endif

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif

#include <climits>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <time.h>
using namespace std;

// ##################################################################
// ##################################################################

setup_grid_NG_MPI::setup_grid_NG_MPI()
{
  spdlog::debug("setup_grid_NG_MPI constructor called");
  return;
}

// ##################################################################
// ##################################################################

setup_grid_NG_MPI::~setup_grid_NG_MPI()
{
  spdlog::debug("setup_grid_NG_MPI destructor called");
  return;
}

// ##################################################################
// ##################################################################

void setup_grid_NG_MPI::setup_NG_grid_levels(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  // call serial version to set global properties of each level.
  int err = 0;
  setup_NG_grid::setup_NG_grid_levels(SimPM);

  // For each level, do a domain decomposition
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (l == 0) {
      err = SimPM.levels[l].sub_domain.decomposeDomain(
          SimPM.ndim, SimPM.levels[l], SimPM.get_pbc_bools());
    }
    else {
      // err = SimPM.levels[l].sub_domain.decomposeDomain(
      //    SimPM.ndim, SimPM.levels[l]);
      SimPM.levels[l].sub_domain = SimPM.levels[0].sub_domain;
      err                        = SimPM.levels[l].sub_domain.decomposeDomain(
          SimPM.ndim, SimPM.levels[l]);
    }
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "PLLEL Init():Decompose Domain!", 0,
          err);
    SimPM.levels[l].sub_domain.set_ReadSingleFile(true);  // legacy option.
  }

  // For each level, get rank of process for parent grid and each
  // child grid.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].sub_domain.set_NG_hierarchy(SimPM, l);
  }

  return;
}

// ##################################################################
// ##################################################################

int setup_grid_NG_MPI::setup_grid(
    vector<class GridBaseClass *> &grid,  ///< grid pointers.
    class SimParams &SimPM                ///< pointer to simulation parameters
)
{
  spdlog::debug("(pion)  Setting up MPI NG grid");

  if (SimPM.ndim < 1 || SimPM.ndim > 3)
    spdlog::error("{}: {}", "Only know 1D,2D,3D methods!", SimPM.ndim);

  spdlog::debug(
      "Setting number of boundary cells == spatial OOA: {}", SimPM.spOOA);

  //
  // Nbc is the depth of the boundary layer around each grid.
  //
  if (SimPM.spOOA == OA2) {
    SimPM.Nbc    = 6;
    SimPM.Nbc_DD = 4;
  }
  else if (SimPM.spOOA == OA1) {
    SimPM.Nbc    = 4;
    SimPM.Nbc_DD = 2;
  }
  else
    spdlog::error("{}: {}", "unhandles spatial order of accuracy", SimPM.spOOA);

  //
  // May need to setup extra data in each cell.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set values cell interface class; note dx changes with level.
  //
  double dx = (SimPM.Xmax[XX] - SimPM.Xmin[XX]) / SimPM.NG[XX];
  CI.set_nlevels(dx, SimPM.grid_nlevels);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);

  //
  // Now we can setup the grid:
  //
  spdlog::debug("(setup_grid_NG_MPI::setup_grid) Setting up grid...");

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "Init: level={},  &grid={}, and grid={}", l, fmt::ptr(&(grid[l])),
        fmt::ptr(grid[l]));

    if (grid[l])
      spdlog::error("{}: {}", "Grid already set up!", fmt::ptr(grid[l]));

    if (SimPM.coord_sys == COORD_CRT)
      grid[l] = new UniformGridParallel(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].sub_domain.get_Xmin(),
          SimPM.levels[l].sub_domain.get_Xmax(),
          SimPM.levels[l].sub_domain.get_directional_Ncells(),
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys == COORD_CYL)
      grid[l] = new uniform_grid_cyl_parallel(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].sub_domain.get_Xmin(),
          SimPM.levels[l].sub_domain.get_Xmax(),
          SimPM.levels[l].sub_domain.get_directional_Ncells(),
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys == COORD_SPH)
      grid[l] = new uniform_grid_sph_parallel(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].sub_domain.get_Xmin(),
          SimPM.levels[l].sub_domain.get_Xmax(),
          SimPM.levels[l].sub_domain.get_directional_Ncells(),
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else
      spdlog::error("{}: {}", "Bad Geometry in setup_grid()", SimPM.coord_sys);

    if (grid[l] == 0)
      spdlog::error(
          "{}: {}", "(setup_grid_NG_MPI::setup_grid)", fmt::ptr(grid[l]));

    spdlog::debug(
        "(setup_grid_NG_MPI::setup_grid) Done. &grid={}, and grid={}, DX={}",
        fmt::ptr(&(grid[l])), fmt::ptr(grid[l]), (grid[l])->DX());
  }

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].grid = grid[l];
    if (l == 0) {
      SimPM.levels[l].parent = 0;
      if (SimPM.grid_nlevels > 1) SimPM.levels[l].child = grid[l + 1];
    }
    else if (l == SimPM.grid_nlevels - 1) {
      if (SimPM.grid_nlevels > 1) SimPM.levels[l].parent = grid[l - 1];
      SimPM.levels[l].child = 0;
    }
    else {
      SimPM.levels[l].parent = grid[l - 1];
      SimPM.levels[l].child  = grid[l + 1];
    }
  }

  set_leaf_cells(grid, SimPM);

  // setup arrays for fluxes into and out of fine grid, and
  // equivalent cells on coarse grid, for making the fluxes
  // consistent across levels.
  setup_flux_vectors(SimPM.grid_nlevels);
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (l != 0) setup_flux_send(SimPM, grid[l], l - 1);
    if (l != SimPM.grid_nlevels - 1) setup_flux_recv(SimPM, grid[l], l + 1);
  }

  return (0);
}  // setup_grid()

// ##################################################################
// ##################################################################

int setup_grid_NG_MPI::setup_raytracing(
    class SimParams &SimPM,              ///< pointer to simulation parameters
    vector<class GridBaseClass *> &grid  ///< vec of grid pointers.
)
{
  if (!SimPM.EP.raytracing) {
    return 0;
  }
  int err = 0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug("setting up raytracing for grid level {}", l);

    err += setup_fixed_grid_pllel::setup_raytracing(SimPM, grid[l]);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "setup_grid_NG_MPI::setup_raytracing()",
          0, err);
  }

  spdlog::debug("NG-MPI setting up evolving RT sources from setup_raytracing");
  err += setup_evolving_RT_sources(SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "setup_grid_NG_MPI::setup_evolving_RT_sources()", 0, err);

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug(
        "NG-MPI l={}: updating evolving RT sources from setup_raytracing", l);

    err +=
        update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid[l]->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}",
          "setup_grid_NG_MPI::update_evolving_RT_sources()", 0, err);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int setup_grid_NG_MPI::boundary_conditions(
    class SimParams &par,                ///< pointer to simulation parameters
    vector<class GridBaseClass *> &grid  ///< vec of grid pointers.
)
{
  spdlog::debug("Setting up BCs in MPI-NG Grid with Nbc={}", par.Nbc);

  int err = setup_NG_grid::boundary_conditions(par, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "setup_grid_NG_MPI::boundary_conditions",
        0, err);
  spdlog::debug("(setup_grid_NG_MPI::boundary_conditions) Done");
  return 0;
}

// ##################################################################
// ##################################################################

int setup_grid_NG_MPI::setup_boundary_structs(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level in NG grid
)
{
  spdlog::debug("Set BC types...");

  // first call fixed grid version
  int err = 0;
  spdlog::debug("setting up serial boundary structs");
  err = setup_fixed_grid::setup_boundary_structs(par, grid, l);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "png::setup_boundary_structs fixed grid",
        0, err);

    //
    // Now check for NG grid boundaries if this grid has a parent
    // grid (i.e. if l > 0).
    //
#ifdef SKIP_C2F_BC
  if (1 == 0) {
#endif
    if (l > 0) {
      spdlog::debug("replacing external BCs with C2F as needed");
      // replace external boundary conditions with one that
      // receives data from a coarser level grid.
      for (int i = 0; i < par.ndim; i++) {
        if (!pconst.equalD(par.levels[l - 1].Xmin[i], par.levels[l].Xmin[i])
            && pconst.equalD(
                   par.levels[l].Xmin[i], grid->Xmin(static_cast<axes>(i)))) {
          spdlog::debug("reassigning neg. bc for axis {} to COARSE_TO_FINE", i);
          grid->BC_bd[2 * i]->itype = COARSE_TO_FINE_RECV;
          grid->BC_bd[2 * i]->type  = "COARSE_TO_FINE_RECV";
        }
        if (!pconst.equalD(par.levels[l - 1].Xmax[i], par.levels[l].Xmax[i])
            && pconst.equalD(
                   par.levels[l].Xmax[i], grid->Xmax(static_cast<axes>(i)))) {
          spdlog::debug("reassigning pos. bc for axis {} to COARSE_TO_FINE", i);
          grid->BC_bd[2 * i + 1]->itype = COARSE_TO_FINE_RECV;
          grid->BC_bd[2 * i + 1]->type  = "COARSE_TO_FINE_RECV";
        }
      }
    }  // if l>0 add C2F recv
#ifdef SKIP_C2F_BC
  }
#endif

  //
  // Now check for NG grid boundaries if this grid has a child
  // grid (i.e. if l < nlevels-1) and, if so, add boundary structs
  // to receive averaged data from a finer grid, and to send
  // external boundary data to the finer grid.  Whether any data
  // need to be sent/received is decided later in the function
  // assign_boundary_data().
  //
  if (l < par.grid_nlevels - 1) {
#ifdef SKIP_F2C_BC
    if (1 == 0) {
#endif
      spdlog::debug(
          "Adding FINE_TO_COARSE_RECV boundary for level {}, current # boundaries: {}",
          l, grid->BC_bd.size());
      struct boundary_data *bd = new boundary_data;
      bd->itype                = FINE_TO_COARSE_RECV;
      bd->type                 = "FINE_TO_COARSE_RECV";
      bd->dir                  = NO;
      bd->ondir                = NO;
      bd->refval               = 0;
      grid->BC_bd.push_back(bd);
#ifdef SKIP_F2C_BC
    }
#endif

#ifdef SKIP_C2F_BC
    if (1 == 0) {
#endif
      spdlog::debug(
          "Adding COARSE_TO_FINE_SEND boundary for level {}, current # boundaries: {}",
          l, grid->BC_bd.size());
      struct boundary_data *bd2 = new boundary_data;
      bd2->itype                = COARSE_TO_FINE_SEND;
      bd2->type                 = "COARSE_TO_FINE_SEND";
      bd2->dir                  = NO;
      bd2->ondir                = NO;
      bd2->refval               = 0;
      grid->BC_bd.push_back(bd2);
#ifdef SKIP_C2F_BC
    }
#endif
  }

#ifdef SKIP_F2C_BC
  if (1 == 0) {
#endif
    if (l > 0) {
      // Add internal boundary to send averaged local data
      // to the parent grid.  Whether any data need to
      // be sent/received is decided later in the function
      // assign_boundary_data().
      // Add it here, because we need to receive F2C data first, then
      // delete any previous F2C send buffers, and then send any new
      // F2C data after the old buffers are gone.
      spdlog::debug(
          "Adding FINE_TO_COARSE_SEND boundary for level {}, current # boundaries: {}",
          l, grid->BC_bd.size());

      struct boundary_data *bd = new boundary_data;
      bd->itype                = FINE_TO_COARSE_SEND;
      bd->type                 = "FINE_TO_COARSE_SEND";
      bd->dir                  = NO;
      bd->ondir                = NO;
      bd->refval               = 0;
      grid->BC_bd.push_back(bd);
    }
#ifdef SKIP_F2C_BC
  }
#endif

  spdlog::debug("BC structs set up");
  for (unsigned int v = 0; v < grid->BC_bd.size(); v++) {
    spdlog::debug(
        "i={}, BC type={}, BC itype={}", v, grid->BC_bd[v]->type,
        grid->BC_bd[v]->itype);
  }

  spdlog::debug("calling pll fixed grid setup function");

  err = setup_fixed_grid_pllel::setup_boundary_structs(par, grid, l);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "png::setup_boundary_structs pll fixed grid", 0, err);

  spdlog::debug("BC structs set up");
  for (unsigned int v = 0; v < grid->BC_bd.size(); v++) {
    spdlog::debug("i={}, BC type={}", v, grid->BC_bd[v]->type);
  }
  return 0;
}

// ##################################################################
// ##################################################################

void setup_grid_NG_MPI::setup_dataio_class(
    class SimParams &par,  ///< simulation parameters
    const int typeOfFile   ///< type of I/O: 1=text,2=fits,5=silo
)
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

#ifdef SILO
    case 5:  // Start from Silo snapshot.
      dataio =
          new dataio_silo_utility(par, "DOUBLE", &(par.levels[0].sub_domain));
      break;
#endif  // if SILO

    default:
      spdlog::error(
          "{}: {}", "sim_control_NG_MPI::Init unhandled filetype", typeOfFile);
  }
  return;
}

// ##################################################################
// ##################################################################
