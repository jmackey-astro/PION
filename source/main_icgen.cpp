/// \author Jonathan Mackey
///
/// Modifications:
///  - 2007-06-28 Started to write an ND shocktube function.  Need to test it.
///  - 2007-07-17 ND MHD shock-cloud interaction problem added to list.  Works
///  well.
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08 to 2008-02 Added more functions.
///  - 2009-12-18 Added in Laser Ablation check.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef. (only serial so far).
/// - 2010.01.06 JM: New stellar wind interface.
/// - 2011.04.14 JM: Added mp_v2_aifa microphysics integration class.
/// - 2011.04.26 JM: removed endl statements.
/// - 2011.04.29 JM: in AddNoise2data() and equilibrate_MP() functions, check
///    for edge cells and internal boundary data, and don't add noise to these
///    cells.  If we have e.g. inflow boundaries, which are set by the edge
///    cell value, we don't that to be a random number!
///
/// - 2011.05.02 JM: AddNoise2data() -- don't add noise to MPI boundaries.
///
/// - 2011.10.13 JM: Added mp_[ex,im]plicit_H classes. (updated 2011.10.22)
/// - 2012.01.23 JM: Added ifdefs for microphysics classes.
/// - 2012.02.07 JM: Added class for Harpreet's 1D to 2D mapping.
/// - 2012.07.25 JM: Fixed bug where noise was added to edge cells in the
///    YZ-face for parallel grids.
/// - 2013.01.10 JM: Added setup lines for StarBench Tests.
/// - 2013.02.15 JM: Added NEW_METALLICITY flag for testing the new
///    microphysics classes.
/// - 2013.02.27 JM: Added IC_read_BBurkhart_data() class for
///    turbulent simulations.
/// - 2013.02.27 JM: Added class for Harpreet's 1D to 3D mapping.
/// - 2013.03.23 JM: Added setup lines for StarBench Tests.
/// - 2013.03.24 JM: Added another StarBench test.
/// - 2013.06.13 JM: Added StarBench_TremblinCooling test.
/// - 2013.08.23 JM: Added new mpv9_HHe module code.
/// - 2014.07.11 JM: Added isothermal noise perturbation option.
/// - 2015.01.(15-26) JM: Added new include statements for new PION version.
/// - 2015.02.03 JM: changed to use IC_base class sub_domain pointer.
/// - 2016.05.02 JM: Changed order of code so that MP is initialised
///    before the "setup" function is called.
/// - 2017.08.03 JM: Don't write IC_filename.silo, just name it like
///    a normal snapshot.
/// - 2018.05.04 JM: moved parallel code to icgen_parallel.cpp

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif /* NDEBUG */

#include "ics/icgen.h"
#include "ics/icgen_base.h"

#include "dataIO/dataio_base.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#ifdef PARALLEL
#ifdef PION_NESTED
#include "dataIO/dataio_silo_utility.h"
#else
#include "dataIO/dataio_silo_MPI.h"
#endif /* PION_NESTED */
#endif /* PARALLEL */
#elif FITS
#include "dataIO/dataio_fits.h"
#ifdef PARALLEL
#include "dataIO/dataio_fits_MPI.h"
#endif /* PARALLEL */
#endif /* SILO */

#ifdef PION_NESTED
#ifdef PARALLEL
#include "grid/setup_grid_NG_MPI.h"
#else
#include "grid/setup_NG_grid.h"
#endif /* PARALLEL */
#else
#ifdef PARALLEL
#include "grid/setup_fixed_grid_MPI.h"
#else
#include "grid/setup_fixed_grid.h"
#endif /* PARALLEL */
#endif /* PION_NESTED */


#include "grid/uniform_grid.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

int main(int argc, char **argv)
{
  int err = 0;

  if (argc < 2) {
    spdlog::error(
        "Error, please give a filename to read IC parameters from.\nUsage <icgen> <paramfile> [ic-filetype]");
    return (1);
  }

  auto max_logfile_size = 1048576 * 5;
  auto max_logfiles     = 3;
#ifdef PARALLEL
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "icgen_pre_mpi", "icgen.log", max_logfile_size, max_logfiles));
#else
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "icgen", "icgen.log", max_logfile_size, max_logfiles));
#endif /* PARALLEL */

#ifdef NDEBUG
  spdlog::set_level(spdlog::level::err);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::info);
#endif

  string pfile = argv[1];
  string icftype;
  if (argc > 2)
    icftype = argv[2];
  else
    icftype = "silo";

  class SimParams SimPM(pfile);
#ifdef PARALLEL
  int r  = SimPM.levels[0].sub_domain.get_myrank();
  int np = SimPM.levels[0].sub_domain.get_nproc();

  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "icgen", "icgen_process_" + to_string(r) + ".log", max_logfile_size,
      max_logfiles));
#endif /* PARALLEL */

  string *args = 0;
  args         = new string[argc];
  for (int i = 0; i < argc; i++)
    args[i] = argv[i];

#ifdef PION_OMP
  // set number of OpenMP threads, if included
  int nth        = 100;  // set to large number initially
  bool found_omp = false;
  for (int i = 0; i < argc; i++) {
    if (args[i].find("omp-nthreads=") != string::npos) {
      nth = atoi((args[i].substr(13)).c_str());
      if (nth > omp_get_num_procs()) {
          spdlog::warn("override: requested too many threads";
        nth = min(nth, omp_get_num_procs());
      }
      spdlog::warn("override: setting OpenMP N-threads to {}", nth);
      found_omp = true;
    }
  }
  if (found_omp)
    omp_set_num_threads(nth);
  else
    omp_set_num_threads(1);
#endif /* PION_OMP */

  class DataIOBase *dataio    = 0;
  class ICsetup_base *ic      = 0;
  class ReadParams *rp        = 0;
  class microphysics_base *MP = 0;

#ifdef PION_NESTED
#ifdef PARALLEL
  class setup_grid_NG_MPI *SimSetup = new setup_grid_NG_MPI();
#else
  class setup_NG_grid *SimSetup = new setup_NG_grid();
#endif /* PARALLEL */
  SimSetup->setup_NG_grid_levels(SimPM);
#else
  /*
   * have to do something with SimPM.levels[0] because this
   * is used to set the local domain size in decomposeDomain
   */
  SimPM.grid_nlevels     = 1;
  SimPM.levels[0].parent = 0;
  SimPM.levels[0].child  = 0;
  SimPM.levels[0].Ncell  = SimPM.Ncell;
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].NG[v] = SimPM.NG[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Range[v] = SimPM.Range[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Xmin[v] = SimPM.Xmin[v];
  for (int v = 0; v < MAX_DIM; v++)
    SimPM.levels[0].Xmax[v] = SimPM.Xmax[v];
  SimPM.levels[0].dx                     = SimPM.Range[XX] / SimPM.NG[XX];
  SimPM.levels[0].simtime                = SimPM.simtime;
  SimPM.levels[0].dt                     = 0.0;
  SimPM.levels[0].multiplier             = 1;
#ifdef PARALLEL
  class setup_fixed_grid_pllel *SimSetup = new setup_fixed_grid_pllel();
  err = SimPM.levels[0].sub_domain.decomposeDomain(SimPM.ndim, SimPM.levels[0]);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "Couldn't Decompose Domain!", 0, err);
#else
  class setup_fixed_grid *SimSetup = new setup_fixed_grid();
#endif /* PARALLEL */
#endif /* PION_NESTED */

  //
  // Set up a grid(s).
  //
  vector<class GridBaseClass *> grid(SimPM.grid_nlevels);
  err      = SimSetup->setup_grid(grid, SimPM);
  SimPM.dx = grid[0]->DX();
  if (!grid[0]) spdlog::error("{}: {}", "Grid setup failed", fmt::ptr(grid[0]));

  //
  // read in what kind of ICs we are setting up.
  //
  rp = new ReadParams;
  if (!rp) spdlog::error("{}: {}", "icgen:: initialising RP", fmt::ptr(rp));
  ;
  err += rp->read_paramfile(pfile);
  if (err) spdlog::error("{}: {}", "Error reading parameterfile", pfile);
  string seek = "ics";
  string ics  = rp->find_parameter(seek);
  setup_ics_type(ics, &ic);
  ic->set_SimPM(&SimPM);
#ifdef PARALLEL
  ic->set_sub_domain_pointer(&SimPM.levels[0].sub_domain);
#endif

#ifdef PION_NESTED
  err = SimSetup->set_equations(SimPM);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(icgen::set_equations) err!=0 Fix me!",
        0, err);
  class FV_solver_base *solver = SimSetup->get_solver_ptr();
#endif

  spdlog::info("MAIN: setting up microphysics module");
  SimSetup->setup_microphysics(SimPM);
  MP = SimSetup->get_mp_ptr();
  ic->set_mp_pointer(MP);
  // ----------------------------------------------------------------

  // have to setup jet simulation before setting up boundary
  // conditions because jet boundary needs some grid data values.
#ifndef PION_NESTED
  if (ics == "Jet" || ics == "JET" || ics == "jet") {
#endif /* PION_NESTED */
    for (int l = 0; l < SimPM.grid_nlevels; l++) {
      err += ic->setup_data(rp, grid[l]);
      if (err) spdlog::error("{}: {}", "Initial conditions setup failed.", err);
    }
#ifndef PION_NESTED
  }
#else
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // Set Ph=P in every cell.
    cell *c = grid[l]->FirstPt();
    do {
      for (int v = 0; v < SimPM.nvar; v++)
        c->Ph[v] = c->P[v];
    } while ((c = grid[l]->NextPt(c)) != 0);
    //
    // If I'm using the GLM method, make sure Psi variable is initialised.
    //
    if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
      for (int l = 0; l < SimPM.grid_nlevels; l++) {
        c = grid[l]->FirstPt();
        do {
          c->P[SI] = c->Ph[SI] = 0.;  // grid->divB(c);
        } while ((c = grid[l]->NextPt(c)) != 0);
      }
    }
  }
#endif /* PION_NESTED */

  //
  // Set up the boundary conditions and data
  //
  SimSetup->boundary_conditions(SimPM, grid);
  if (err) spdlog::error("{}: {}", "icgen Couldn't set up boundaries.", err);

#ifndef PION_NESTED
  err += SimSetup->setup_raytracing(SimPM, grid[0]);
  err += SimSetup->setup_evolving_RT_sources(SimPM);
  err +=
      SimSetup->update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  if (err)
    spdlog::error(
        "{}: {}", "icgen: Failed to setup raytracer and/or microphysics", err);
  // ----------------------------------------------------------------
  // call "setup" to set up the data on the computational grid.
  err += ic->setup_data(rp, grid[0]);
  if (err) spdlog::error("{}: {}", "Initial conditions setup failed.", err);
    // ----------------------------------------------------------------

#ifndef PARALLEL
  err = SimSetup->assign_boundary_data(SimPM, 0, grid[0], MP);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen::assign_boundary_data", 0, err);
#endif /* PARALLEL */
#else
  err += SimSetup->setup_raytracing(SimPM, grid);
  if (err) spdlog::error("{}: {}", "icgen-ng: Failed to setup raytracer", err);

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // cout <<"icgen_NG: assigning boundary data for level "<<l<<"\n";
    err = SimSetup->assign_boundary_data(SimPM, l, grid[l], MP);
#ifdef PARALLEL
    SimPM.levels[0].sub_domain.barrier("level assign boundary data");
#endif /* PARALLEL */
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "icgen-ng::assign_boundary_data", 0,
          err);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // cout <<"updating external boundaries for level "<<l<<"\n";
    err += SimSetup->TimeUpdateExternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen-ng: error from bounday update", 0,
        err);
    // ----------------------------------------------------------------

#ifdef PARALLEL
    // ----------------------------------------------------------------
#ifndef NDEBUG
  spdlog::info("icgen-ng: updating C2F boundaries");
#endif /* NDEBUG */
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifndef NDEBUG
    spdlog::debug("icgen-ng updating C2F boundaries for level {}", l);
#endif /* NDEBUG */
    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += SimSetup->BC_update_COARSE_TO_FINE_SEND(
              SimPM, grid[l], solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
#ifndef NDEBUG
    spdlog::info("icgen-ng: updating C2F boundaries, sent, now recv");
#endif /* NDEBUG */
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += SimSetup->BC_update_COARSE_TO_FINE_RECV(
              SimPM, solver, l, grid[l]->BC_bd[i], SimPM.levels[l].step);
        }
      }
    }
#ifndef NDEBUG
    spdlog::info("icgen-ng: updating C2F boundaries: done with level");
#endif /* NDEBUG */
  }
  spdlog::info("icgen-ng: updating C2F boundaries done");
  SimSetup->BC_COARSE_TO_FINE_SEND_clear_sends(SimPM.levels[0].sub_domain);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen-ng: error from boundary update", 0,
        err);
    // ----------------------------------------------------------------

    // ----------------------------------------------------------------
#ifndef NDEBUG
  spdlog::info("icgen-ng: updating external boundaries");
#endif /* NDEBUG */
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err += SimSetup->TimeUpdateExternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen-ng: error from boundary update", 0,
        err);
    // ----------------------------------------------------------------
#endif /* PARALLEL */

  // ----------------------------------------------------------------
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
    // cout <<"updating internal boundaries for level "<<l<<"\n";
    err += SimSetup->TimeUpdateInternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen-ng: error from bounday update", 0,
        err);
    // ----------------------------------------------------------------

#ifdef PARALLEL
    // ----------------------------------------------------------------
    // update fine-to-coarse level boundaries
#ifndef NDEBUG
  spdlog::info("icgen-ng: updating F2C boundaries");
#endif /* NDEBUG */
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifndef NDEBUG
    spdlog::debug("icgen-ng updating F2C boundaries for level {}", l);
#endif /* NDEBUG */
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_SEND) {
          err += SimSetup->BC_update_FINE_TO_COARSE_SEND(
              SimPM, solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == FINE_TO_COARSE_RECV) {
          err += SimSetup->BC_update_FINE_TO_COARSE_RECV(
              SimPM, solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
  }
  SimSetup->BC_FINE_TO_COARSE_SEND_clear_sends(SimPM.levels[0].sub_domain);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "icgen-ng: error from boundary update", 0,
        err);
    // ----------------------------------------------------------------
#endif /* PARALLEL */
#endif /* PION_NESTED */

  // ----------------------------------------------------------------
  // if data initialised ok, maybe we need to equilibrate the
  // chemistry...
  //
  if (SimPM.ntracer > 0 && (SimPM.EP.chemistry)) {
    spdlog::info("MAIN: equilibrating the chemical species.");
    if (!MP) spdlog::error("{}: {}", "microphysics init", fmt::ptr(MP));

    // first avoid cooling the gas in getting to equilbrium, by
    // setting update_erg to false.
    bool uerg           = SimPM.EP.update_erg;
    SimPM.EP.update_erg = false;
    err                 = ic->equilibrate_MP(grid[0], MP, rp, SimPM);
    if (err)
      spdlog::error(
          "{}: {}", "setting chemical states to equilibrium failed", err);

    SimPM.EP.update_erg = uerg;
    spdlog::info("MAIN: finished equilibrating the chemical species.");
  }
  // ----------------------------------------------------------------

  spdlog::debug("IC file-type is ", icftype);
  seek          = "OutputFile";
  string icfile = rp->find_parameter(seek);
  if (icfile == "") {
    spdlog::warn("WARNING: no filename for ic file.  writing to IC_temp.****");
    icfile = "IC_temp";
  }

  dataio = 0;  // zero the class pointer.
#ifdef SILO
  if (icftype == "silo") {
    spdlog::info("WRITING SILO FILE: {}", icfile);
#ifdef PARALLEL
#ifdef PION_NESTED
    dataio =
        new dataio_silo_utility(SimPM, "DOUBLE", &SimPM.levels[0].sub_domain);
#else
    dataio =
        new dataio_silo_pllel(SimPM, "DOUBLE", &SimPM.levels[0].sub_domain);
#endif /* PION_NESTED */
#else
    dataio = new dataio_silo(SimPM, "DOUBLE");
#endif /* PARALLEL */
  }
#elif !defined(PION_NESTED) /* elifndef would be nice */
#ifndef PARALLEL
  if (icftype == "text") {
    spdlog::info("WRITING ASCII TEXT FILE: {}", icfile);
  }
#endif

#ifdef FITS
  if (icftype == "fits") {
    spdlog::info("WRITING FITS FILE: {}", icfile);
#ifdef PARALLEL
    dataio = new DataIOFits_pllel(SimPM, &SimPM.levels[0].sub_domain);
#else
    dataio = new DataIOFits(SimPM);
#endif /* PARALLEL */
  }
#endif /* FITS */
#endif /* PION_NESTED */

  if (!dataio) spdlog::error("{}: {}", "IO class initialisation: ", icftype);
  err = dataio->OutputData(icfile, grid, SimPM, 0);
  if (err) spdlog::error("{}: {}", "File write error", err);
  delete dataio;
  spdlog::info("{} FILE WRITTEN... exiting", icftype);

  // delete everything and return
  if (MP) {
    delete MP;
    MP = 0;
  }
  if (rp) {
    delete rp;
    rp = 0;
  }  // Delete the read_parameters class.
  if (ic) {
    delete ic;
    ic = 0;
  }
  if (SimSetup) {
    delete SimSetup;
    SimSetup = 0;
  }
  for (auto g : grid)
    delete g;

  //
  // Also delete any dynamic memory in the stellarwind_list in
  // global.h, since we may have added wind sources in the
  // read-params functions.
  //
  while (SWP.params.size() > 0) {
    int i                           = static_cast<int>(SWP.params.size()) - 1;
    struct stellarwind_params *temp = SWP.params[i];
    SWP.params.pop_back();    // remove struct from list.
    temp = mem.myfree(temp);  // delete struct.
  }

  delete[] args;
  args = 0;
  return err;
}

// ##################################################################
// ##################################################################
