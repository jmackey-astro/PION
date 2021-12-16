/// \file setup_NG_grid.cpp
///
/// \brief Class for setting up NG grids.
///
/// \author Jonathan Mackey
///
///
/// Modifications:
/// - 2018.05.01 JM derived class from setup_fixed_grid.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "raytracing/raytracer_SC.h"
#include "tools/command_line_interface.h"


#include "grid/setup_NG_grid.h"
#include "grid/uniform_grid.h"
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include "microphysics/MPv3.h"
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"
#include "microphysics/MPv8.h"
#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/MPv9.h"
#endif

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif

#include "dataIO/dataio_base.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif  // if FITS

#include <climits>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <time.h>
using namespace std;

// ##################################################################
// ##################################################################

setup_NG_grid::setup_NG_grid() {}

// ##################################################################
// ##################################################################

setup_NG_grid::~setup_NG_grid() {}

// ##################################################################
// ##################################################################

void setup_NG_grid::setup_NG_grid_levels(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  spdlog::debug("(pion-ng)  Setting up nested grid parameters");
  // first make sure that NG_centre[] is oriented such that an
  // oct-tree structure works.  Centre should be xmin + i/4 of the
  // full domain, where i is in [0,4]
  for (int d = 0; d < SimPM.ndim; d++) {
    double f = 4.0 * (SimPM.NG_centre[d] - SimPM.Xmin[d]) / SimPM.Range[d];
    f        = fmod(f, 1.0);
    if (!pconst.equalD(f, 0.0)) {
      spdlog::debug(
          "setup_NG_grid_levels:  axis={}, resetting NG_centre to i/4 of the domain... current NG_centre={}",
          d, SimPM.NG_centre[d]);
      if (f > 0.5) {
        SimPM.NG_centre[d] += (1.0 - f) * SimPM.Range[d] / 4.0;
      }
      else {
        SimPM.NG_centre[d] -= f * SimPM.Range[d] / 4.0;
      }
      spdlog::debug(" reset to {}", SimPM.NG_centre[d]);
    }
  }

  //
  // populate "levels" struct in SimPM based on NG grid parameters.
  //

  for (int i = 0; i < SimPM.grid_nlevels; i++) {
    SimPM.levels[i].parent = 0;
    SimPM.levels[i].child  = 0;
    // Refine only along directions specified (NG_refine[dir]=1)
    // Otherwise need to double number of cells in each refined level
    SimPM.levels[i].Ncell = 1;
    for (int v = 0; v < MAX_DIM; v++) {
      if (SimPM.NG_refine[v] == 1)
        SimPM.levels[i].NG[v] = SimPM.NG[v];
      else
        SimPM.levels[i].NG[v] = SimPM.NG[v] * pow(2, i);
      SimPM.levels[i].Ncell *= SimPM.levels[i].NG[v];
    }
    if (i == 0) {
      for (int v = 0; v < MAX_DIM; v++)
        SimPM.levels[i].Range[v] = SimPM.Range[v];
      for (int v = 0; v < MAX_DIM; v++)
        SimPM.levels[i].Xmin[v] = SimPM.Xmin[v];
      for (int v = 0; v < MAX_DIM; v++)
        SimPM.levels[i].Xmax[v] = SimPM.Xmax[v];
      SimPM.levels[i].dx = SimPM.Range[XX] / SimPM.NG[XX];
    }
    else {
      for (int v = 0; v < MAX_DIM; v++) {
        if (SimPM.NG_refine[v] == 1) {
          SimPM.levels[i].Range[v] = 0.5 * SimPM.levels[i - 1].Range[v];
          SimPM.levels[i].Xmin[v] =
              0.5 * (SimPM.levels[i - 1].Xmin[v] + SimPM.NG_centre[v]);
          SimPM.levels[i].Xmax[v] =
              0.5 * (SimPM.levels[i - 1].Xmax[v] + SimPM.NG_centre[v]);
        }
        else {
          SimPM.levels[i].Range[v] = SimPM.levels[i - 1].Range[v];
          SimPM.levels[i].Xmin[v]  = SimPM.levels[i - 1].Xmin[v];
          SimPM.levels[i].Xmax[v]  = SimPM.levels[i - 1].Xmax[v];
        }
      }
      SimPM.levels[i].dx = 0.5 * SimPM.levels[i - 1].dx;
    }
    SimPM.levels[i].simtime = SimPM.simtime;
    SimPM.levels[i].dt      = 0.0;

    if (i == 0)
      SimPM.levels[i].multiplier = 1;
    else
      SimPM.levels[i].multiplier = 2 * SimPM.levels[i - 1].multiplier;
  }

  for (int i = SimPM.grid_nlevels - 1; i >= 0; i--) {
    if (i == SimPM.grid_nlevels - 1)
      SimPM.levels[i].step = SimPM.timestep;
    else
      SimPM.levels[i].step = SimPM.levels[i + 1].step / 2;

#ifndef NDEBUG
    ostringstream temp;
    temp << i;
    string lv = "level " + temp.str();
    string t2 = lv + "_Range";
    spdlog::debug("{} : {}", t2, SimPM.levels[i].Range);
    t2 = lv + "_Xmin";
    spdlog::debug("{} : {}", t2, SimPM.levels[i].Xmin);
    t2 = lv + "_Xmax";
    spdlog::debug("{} : {}", t2, SimPM.levels[i].Xmax);
    spdlog::debug(
        "\t\tdx={}, step={}", SimPM.levels[i].dx, SimPM.levels[i].step);
#endif
  }

  return;
}

// ##################################################################
// ##################################################################

int setup_NG_grid::setup_grid(
    vector<class GridBaseClass *> &grid,  ///< grid pointers.
    class SimParams &SimPM                ///< pointer to simulation parameters
)
{
  spdlog::info("(pion ng)  Setting up computational grid");

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
  // Set values cell interface class.
  //
  double dx = (SimPM.Xmax[XX] - SimPM.Xmin[XX]) / SimPM.NG[XX];
  CI.set_nlevels(dx, SimPM.grid_nlevels);  // sets dx on all levels.
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);

  //
  // Now we can setup the grid:
  //
  spdlog::info("(setup_NG_grid::setup_grid) Setting up grid...");
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (grid[l])
      spdlog::error("{}: {}", "Grid already set up!", fmt::ptr(grid[l]));

    if (SimPM.coord_sys == COORD_CRT)
      grid[l] = new UniformGrid(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.levels[l].NG,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys == COORD_CYL)
      grid[l] = new uniform_grid_cyl(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.levels[l].NG,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else if (SimPM.coord_sys == COORD_SPH)
      grid[l] = new uniform_grid_sph(
          SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.levels[l].NG,
          SimPM.levels[l].Xmin, SimPM.levels[l].Xmax, SimPM.Xmin, SimPM.Xmax);
    else
      spdlog::error("{}: {}", "Bad Geometry in setup_grid()", SimPM.coord_sys);

    if (grid[l] == 0)
      spdlog::error(
          "{}: {}", "(setup_NG_grid::setup_grid) Couldn't assign data!",
          fmt::ptr(grid[l]));

#ifndef NDEBUG
    spdlog::debug(
        "(setup_NG_grid::setup_grid) Done. &grid={}, and grid={}, DX={}",
        fmt::ptr(&(grid[l])), fmt::ptr(grid[l]), (grid[l])->DX());
    dp.grid = (grid[l]);
#endif
  }

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    SimPM.levels[l].grid = grid[l];
    if (l == 0) {
      SimPM.levels[l].parent = 0;
      SimPM.levels[l].child  = grid[l + 1];
    }
    else if (l == SimPM.grid_nlevels - 1) {
      SimPM.levels[l].parent = grid[l - 1];
      SimPM.levels[l].child  = 0;
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

void setup_NG_grid::set_leaf_cells(
    vector<class GridBaseClass *> &grid,  ///< grid pointers.
    class SimParams &par                  ///< pointer to simulation parameters
)
{
  //
  // if there is an interface region, set a flag on the grid cells
  // that are not leaf cells.
  //
  std::array<int, MAX_DIM> sxmin, sxmax, lxmin, lxmax;
  bool notleaf;
  CI.get_ipos_vec(par.levels[0].Xmin, sxmin);
  CI.get_ipos_vec(par.levels[0].Xmax, sxmax);

  for (int l = 0; l < par.grid_nlevels; l++) {
    class cell *c = grid[l]->FirstPt_All();
    do {
      // if outside domain, then cell is not a leaf.
      for (int v = 0; v < par.ndim; v++) {
        if (c->pos[v] < sxmin[v] || c->pos[v] > sxmax[v]) c->isleaf = false;
      }
      if (l < par.grid_nlevels - 1) {
        notleaf = true;
        CI.get_ipos_vec(par.levels[l + 1].Xmin, lxmin);
        CI.get_ipos_vec(par.levels[l + 1].Xmax, lxmax);
        for (int v = 0; v < par.ndim; v++) {
          if (c->pos[v] > lxmax[v] || c->pos[v] < lxmin[v]) notleaf = false;
        }
        if (notleaf) c->isleaf = false;
      }
      if (!c->isleaf) c->timestep = false;

    } while ((c = grid[l]->NextPt_All(c)) != 0);
  }
  return;
}

// ##################################################################
// ##################################################################

int setup_NG_grid::setup_raytracing(
    class SimParams &SimPM,              ///< simulation parameters
    vector<class GridBaseClass *> &grid  ///< vector of grids.
)
{
  if (!SimPM.EP.raytracing) {
    return 0;
  }
  int err = 0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err += setup_fixed_grid::setup_raytracing(SimPM, grid[l]);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "setup_NG_grid::setup_raytracing()", 0,
          err);
  }

  err += setup_evolving_RT_sources(SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "setup_NG_grid::setup_evolving_RT_sources()", 0, err);

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err +=
        update_evolving_RT_sources(SimPM, SimPM.levels[l].simtime, grid[l]->RT);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "setup_NG_grid::update_RT_sources()", 0,
          err);
  }
  return 0;
}

// ##################################################################
// ##################################################################

int setup_NG_grid::boundary_conditions(
    class SimParams &par,                ///< pointer to simulation parameters
    vector<class GridBaseClass *> &grid  ///< grid pointers.
)
{
  // For uniform fixed cartesian grid.
  spdlog::debug("Setting up BCs in Grid with Nbc={}", par.Nbc);
  int err = 0;
  for (int l = 0; l < par.grid_nlevels; l++) {
    //
    // Choose what BCs to set up based on BC strings.
    //
    err = setup_boundary_structs(par, grid[l], l);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "sng::boundary_conditions sb_structs",
          0, err);

    err = grid[l]->SetupBCs(par);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "sng::boundary_conditions SetupBCs", 0,
          err);
  }
  spdlog::info("(setup_NG_grid::boundary_conditions) Done");
  return 0;
}

// ##################################################################
// ##################################################################

int setup_NG_grid::setup_boundary_structs(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level in NG grid
)
{
  spdlog::debug("Set BC types...");

  // first call fixed grid version
  int err = setup_fixed_grid::setup_boundary_structs(par, grid, l);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "sng::setup_boundary_structs fixed grid",
        0, err);

  //
  // Now check for NG grid boundaries if this grid has a parent
  // grid (i.e. if l > 0).
  //
  if (l > 0) {
    for (int i = 0; i < par.ndim; i++) {
      if (!pconst.equalD(par.levels[l - 1].Xmin[i], par.levels[l].Xmin[i])) {
        spdlog::debug("reassigning neg. bc for axis {} to COARSE_TO_FINE", i);
        grid->BC_bd[2 * i]->itype = COARSE_TO_FINE;
        grid->BC_bd[2 * i]->type  = "COARSE_TO_FINE";
      }
      if (!pconst.equalD(par.levels[l - 1].Xmax[i], par.levels[l].Xmax[i])) {
        spdlog::debug("reassigning pos. bc for axis {} to COARSE_TO_FINE", i);
        grid->BC_bd[2 * i + 1]->itype = COARSE_TO_FINE;
        grid->BC_bd[2 * i + 1]->type  = "COARSE_TO_FINE";
      }
    }
  }

  //
  // Now check for NG grid boundaries if this grid has a child
  // grid (i.e. if l < nlevels), and set non-leaf cells to boundary
  // data.
  //
  // for serial code, a grid can only have one child grid, because
  // a single grid covers the full domain at each level.  So we go
  // through all cells on the grid, and if they are within the domain
  // of the child grid (if it exists) then we label the cells as
  // boundary data and add them to a new FINE_TO_COARSE boundary.
  //
  if (l < par.grid_nlevels - 1) {
    std::array<double, MAX_DIM> cxmin, cxmax, cpos;
    bool within_child = true;
    size_t ct         = 0;
    for (int v = 0; v < par.ndim; v++) {
      cxmin[v] = par.levels[l + 1].Xmin[v];
      cxmax[v] = par.levels[l + 1].Xmax[v];
    }
    struct boundary_data *bd = new boundary_data;
    bd->itype                = FINE_TO_COARSE;
    bd->type                 = "FINE_TO_COARSE";
    bd->dir                  = NO;
    bd->ondir                = NO;
    bd->refval               = 0;
    bd->NGrecvF2C.resize(1);
    cell *c = grid->FirstPt();
    do {
      within_child = true;
      CI.get_dpos(c, cpos);
      for (int v = 0; v < par.ndim; v++) {
        if (cpos[v] < cxmin[v] || cpos[v] > cxmax[v]) within_child = false;
      }
      if (within_child) {
        c->isbd = true;
        bd->NGrecvF2C[0].push_back(c);
        ct++;
      }
    } while ((c = grid->NextPt(c)) != 0);
    spdlog::debug(
        "Got {} cells for FINE_TO_COARSE boundary, ", ct, bd->data.size());

    grid->BC_bd.push_back(bd);

    spdlog::debug(
        "BC_data: {}",
        grid->BC_bd[grid->BC_bd.size() - 1]->NGrecvF2C[0].size());
  }

  spdlog::debug("BC structs set up");
  return 0;
}

// ##################################################################
// ##################################################################

void setup_NG_grid::setup_dataio_class(
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
      dataio = new dataio_silo(par, "DOUBLE");
      break;
#endif  // if SILO

    default:
      spdlog::error(
          "{}: {}", "sim_control_NG::Init unhandled filetype", typeOfFile);
  }
  return;
}

// ##################################################################
// ##################################################################
