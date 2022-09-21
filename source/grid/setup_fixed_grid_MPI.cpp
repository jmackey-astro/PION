/// \file setup_fixed_grid_MPI.cpp
///
/// \brief MPI-Parallel Class for setting up fixed grids.
///
/// \author Jonathan Mackey
///
/// This file contains the definitions of the member functions for
/// the "setup_fixed_grid_pllel" class, which is for setting up grids
///
/// Modifications:
/// - 2015.02.18 JM: new file for setting up parallel grids.
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "grid/uniform_grid_pllel.h"
#include "raytracing/raytracer_SC.h"
#include "raytracing/raytracer_SC_pllel.h"
#include "setup_fixed_grid.h"
#include "setup_fixed_grid_MPI.h"
#include "sub_domain/sub_domain.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo_utility.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits_MPI.h"
#endif  // if FITS

using namespace std;

#ifdef PARALLEL

// ##################################################################
// ##################################################################

setup_fixed_grid_pllel::setup_fixed_grid_pllel()
{
  spdlog::debug("setup_fixed_grid_pllel constructor");
}

// ##################################################################
// ##################################################################

setup_fixed_grid_pllel::~setup_fixed_grid_pllel()
{
  spdlog::debug("setup_fixed_grid_pllel destructor");
}

// ##################################################################
// ##################################################################

int setup_fixed_grid_pllel::setup_grid(
    vector<class GridBaseClass *> &g,  ///< grid pointers.
    class SimParams &SimPM             ///< simulation parameters
)
{
  spdlog::debug("setup_fixed_grid_pllel: setting up parallel grid");

  class GridBaseClass **grid   = &g[0];
  class Sub_domain *sub_domain = &(SimPM.levels[0].sub_domain);

  if (SimPM.gridType != 1) {
    SimPM.gridType = 1;
  }
  if (SimPM.ndim < 1 || SimPM.ndim > 3)
    spdlog::error("{}: {}", "Only know 1D,2D,3D methods!", SimPM.ndim);

  //
  // Nbc is the depth of the boundary layer.
  //
  spdlog::debug(
      "Setting number of boundary cells == spatial OO: {}", SimPM.spOOA);

  if (SimPM.spOOA == OA2) {
    SimPM.Nbc    = 2;
    SimPM.Nbc_DD = 2;
  }
  else if (SimPM.spOOA == OA1) {
    SimPM.Nbc    = 1;
    SimPM.Nbc_DD = 1;
  }
  else
    spdlog::error(
        "{}: {}", "Spatial order of accuracy unhandled by boundary conditions!",
        SimPM.spOOA);

  // Force Nbc=1 if using Lax-Friedrichs flux.
  if (SimPM.solverType == FLUX_LF) {
    SimPM.spOOA = SimPM.tmOOA = OA1;
    SimPM.Nbc                 = 1;
  }

  //
  // May need to setup extra data in each cell for ray-tracing
  // optical depths and/or viscosity variables.  Cells cannot be
  // created unless this the number of such extra variables has been
  // set.
  //
  setup_cell_extra_data(SimPM);

  //
  // Set Cell dx in cell interface class, and also xmin.
  //
  double dx = (SimPM.Xmax[XX] - SimPM.Xmin[XX]) / SimPM.NG[XX];
  CI.set_nlevels(dx, 1);
  CI.set_ndim(SimPM.ndim);
  CI.set_nvar(SimPM.nvar);
  CI.set_xmin(SimPM.Xmin);
  //
  // Now set up the parallel uniform grid.
  //
  spdlog::debug("(setup_fixed_grid_pllel::setup_grid) Setting up grid..");

  if (SimPM.coord_sys == COORD_CRT) {
    *grid = new UniformGridParallel(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
        sub_domain->get_Xmin(), sub_domain->get_Xmax(),
        sub_domain->get_directional_Ncells(), SimPM.Xmin, SimPM.Xmax,
        SimPM.Xmin, SimPM.Xmax);
  }
  else if (SimPM.coord_sys == COORD_CYL) {
    *grid = new uniform_grid_cyl_parallel(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
        sub_domain->get_Xmin(), sub_domain->get_Xmax(),
        sub_domain->get_directional_Ncells(), SimPM.Xmin, SimPM.Xmax,
        SimPM.Xmin, SimPM.Xmax);
  }
  else if (SimPM.coord_sys == COORD_SPH) {
    *grid = new uniform_grid_sph_parallel(
        SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Nbc,
        sub_domain->get_Xmin(), sub_domain->get_Xmax(),
        sub_domain->get_directional_Ncells(), SimPM.Xmin, SimPM.Xmax,
        SimPM.Xmin, SimPM.Xmax);
  }
  else {
    spdlog::error("{}: {}", "Bad Geometry in setup_grid()", SimPM.coord_sys);
  }

  (*grid)->set_level(0);

  if (*grid == 0)
    spdlog::error(
        "{}: {}", "(setup_fixed_grid_pllel::setup_grid) Couldn't assign data!",
        fmt::ptr(*grid));

  spdlog::debug("(setup_fixed_grid_pllel::setup_grid) Done");
  spdlog::debug("grid DX = {0}", (*grid)->DX());
  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid_pllel::setup_raytracing(
    class SimParams &SimPM,    ///< pointer to simulation parameters
    class GridBaseClass *grid  ///< pointer to grid
)
{
  //
  // This function is identical to the serial setup function, except
  // that it sets up MPI-aware versions of the raytracers.
  //
  if (!SimPM.EP.raytracing || SimPM.RS.Nsources == 0) {
    return 0;
  }
  spdlog::debug("(pion-mpi)  Setting up raytracing on leve");

  if (!MP) {
    spdlog::error(
        "{}: {}", "can't do raytracing without microphysics", fmt::ptr(MP));
    exit(1);
  }
  grid->RT = 0;
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays = true;
  int dir            = -1;
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++) {
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays = false;
    //
    // source is at infinity, so make sure all sources at infinity have rays
    // travelling in the same direction (by checking direction to source).
    //
    else {
      for (int i = 0; i < SimPM.ndim; i++) {
        if (fabs(SimPM.RS.sources[isrc].pos[i]) > 1.e99) {
          if (dir == -1)
            dir = i;
          else if (dir != i)
            parallel_rays = false;
        }
      }
    }
  }  // loop over sources.
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  parallel_rays = false;
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    grid->RT = new raytracer_USC_infinity(
        grid, MP, SimPM.ndim, SimPM.coord_sys, SimPM.nvar, SimPM.ftr);
    if (!grid->RT) {
      spdlog::error(
          "{}: {}", "init pllel-rays raytracer error", fmt::ptr(grid->RT));
      exit(1);
    }
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    grid->RT = new raytracer_USC_pllel(
        grid, MP, &SimPM, &(SimPM.levels[0].sub_domain), SimPM.ndim,
        SimPM.coord_sys, SimPM.nvar, SimPM.ftr, SimPM.RS.Nsources);
    if (!grid->RT)
      spdlog::error("{}: {}", "init raytracer error 2", fmt::ptr(grid->RT));
  }

  //
  // Now add the sources to the raytracer.
  //
  int ion_count = 0, uv_count = 0, dif_count = 0;
  for (int isrc = 0; isrc < SimPM.RS.Nsources; isrc++) {

    // see if source is on the domain or not.  Set flag accordingly.
    SimPM.RS.sources[isrc].ongrid = true;
    for (int d = 0; d < SimPM.ndim; d++) {
      if (SimPM.RS.sources[isrc].pos[d] < SimPM.levels[0].Xmin[d]
          || SimPM.RS.sources[isrc].pos[d] > SimPM.levels[0].Xmax[d])
        SimPM.RS.sources[isrc].ongrid = false;
    }

    if (SimPM.RS.sources[isrc].type == RT_SRC_SINGLE) {
      //
      // single sources have a flux (if at infinity) or a luminosity (if
      // point sources.
      //
      int s = grid->RT->Add_Source(&(SimPM.RS.sources[isrc]));
      spdlog::debug("Adding IONISING or UV single-source with id: {}", s);
      if (SimPM.RS.sources[isrc].effect == RT_EFFECT_PION_MONO
          || SimPM.RS.sources[isrc].effect == RT_EFFECT_MFION)
        ion_count++;
      else
        uv_count++;
    }  // if ionising source
    else {
      // note that diffuse radiation must be at infinity, and the strength
      // is assumed to be an intensity not a flux, so it is multiplied by
      // a solid angle appropriate to its location in order to get a flux.
      int s = grid->RT->Add_Source(&(SimPM.RS.sources[isrc]));

      spdlog::debug("Adding DIFFUSE radiation source with id: {}", s);

      uv_count++;
      dif_count++;
    }  // if diffuse source
  }    // loop over sources
  if (ion_count > 1) {
    spdlog::error(
        "{}: {}",
        "Can only have one ionising source for currently implemented method",
        ion_count);
  }

  spdlog::debug(
      "Added {} ionising and {} non-ionising radiation sources, of which {} are diffuse radiation",
      ion_count, uv_count, dif_count);

  grid->RT->Print_SourceList();

  //
  // Now that we have added all of the sources, we query the raytracer to get
  // all of the source properties into structs for the microphysics calls.
  // NOTE THAT IF THE NUMBER OF SOURCES OR THEIR PROPERTIES CHANGE OVER TIME,
  // I WILL HAVE TO WRITE NEW CODE TO UPDATE THIS!
  //
  FVI_nheat = grid->RT->N_heating_sources();
  FVI_nion  = grid->RT->N_ionising_sources();
  FVI_heating_srcs.resize(FVI_nheat);
  FVI_ionising_srcs.resize(FVI_nion);
  grid->RT->populate_UVheating_src_list(FVI_heating_srcs);
  grid->RT->populate_ionising_src_list(FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (grid->RT->type_of_RT_integration() == RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (
      grid->RT && grid->RT->type_of_RT_integration() == RT_UPDATE_IMPLICIT
      && SimPM.EP.MP_timestep_limit == 5) {
    // For implicit updates to limit by xdot and/or edot
    // Here the raytracing has not already been done, so we call it here.
    FVI_need_column_densities_4dt = true;
  }
  else {
    FVI_need_column_densities_4dt = false;
  }

  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid_pllel::boundary_conditions(
    class SimParams &par,                ///< pointer to simulation parameters
    vector<class GridBaseClass *> &grid  ///< grid pointers.
    // class GridBaseClass *grid ///< pointer to grid.
)
{
  // For uniform fixed cartesian grid.
  spdlog::debug("Setting up BCs in Grid with Nbc={}", par.Nbc);

  //
  // Choose what BCs to set up based on BC strings.
  //
  int err = setup_boundary_structs(par, grid[0], 0);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "sfg::boundary_conditions::sb_structs", 0,
        err);

  //
  // Ask grid to set up data for external boundaries.
  //
  err = grid[0]->SetupBCs(par);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "sfg::boundary_conditions::SetupBCs", 0,
        err);

  spdlog::debug("(setup_fixed_grid::boundary_conditions) Done");
  return 0;
}

// ##################################################################
// ##################################################################

int setup_fixed_grid_pllel::setup_boundary_structs(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int l                 ///< level
)
{
  string fname = "setup_fixed_grid_pllel::setup_boundary_structs";
  spdlog::debug("PLLEL: Set BC types..");
  //
  // call serial version of setBCtypes, to set up the boundaries
  //
  int err = setup_fixed_grid::setup_boundary_structs(par, grid, l);
  if (err) {
    spdlog::error("{}: {}", "sfg_pllel::setup_boundary_structs:: serial", err);
  }

  //
  // Now go through the 6 directions, and see if we need to replace
  // any edge boundaries with MPI communication boundaries, if the
  // local grid does not reach the global grid boundary.
  //
  int i = 0;
  string temp;
  // Loop through boundaries, and if local boundary is not sim boundary,
  // set it to be a parallel boundary.
  for (i = 0; i < par.ndim; i++) {
    if (!pconst.equalD(
            grid->Xmin(static_cast<axes>(i)), par.levels[l].Xmin[i])) {
      // local xmin is not Sim xmin, so it's an mpi boundary
      grid->BC_bd[2 * i]->itype = BCMPI;
      grid->BC_bd[2 * i]->type  = "BCMPI";
    }
    if (!pconst.equalD(
            grid->Xmax(static_cast<axes>(i)), par.levels[l].Xmax[i])) {
      // local xmax is not Sim xmin, so it's an mpi boundary
      grid->BC_bd[2 * i + 1]->itype = BCMPI;
      grid->BC_bd[2 * i + 1]->type  = "BCMPI";
    }
  }

  spdlog::debug("PLLEL: BC types and data set up");
  spdlog::debug(
      "Neighbouring processors: {}",
      par.levels[0].sub_domain.get_neighbour_ranks());

  return (0);
}

// ##################################################################
// ##################################################################

void setup_fixed_grid_pllel::setup_dataio_class(
    class SimParams &par,  ///< simulation parameters
    const int typeOfFile   ///< type of I/O: 1=text,2=fits,5=silo
)
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

    case 1:  // Start From ASCII Parameterfile.
      spdlog::error("{}: {}", "No text file for parallel I/O!", typeOfFile);
      break;

#ifdef FITS
    case 2:  // Start from FITS restartfile
      dataio = new DataIOFits_pllel(par, &(par.levels[0].sub_domain));
      break;
#endif  // if FITS

#ifdef SILO
    case 5:  // Start from Silo ICfile or restart file.
      dataio =
          new dataio_silo_utility(par, "DOUBLE", &(par.levels[0].sub_domain));
      break;
#endif  // if SILO
    default:
      spdlog::error(
          "{}: {}", "sim_control::Init unhandled filetype", typeOfFile);
  }
  return;
}

// ##################################################################
// ##################################################################

#endif  // PARALLEL
