/// \file icgen_NG.cpp
/// \brief Program to generate Initial Conditions for a NG grid.
/// \author Jonathan Mackey
///
/// Modifications:
/// - 2018.05.04 JM: adapted from icgen.cpp

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "ics/get_sim_info.h"
#include "ics/icgen.h"
#include "ics/icgen_base.h"

#include "dataIO/dataio_base.h"
#ifdef SILO
#include "dataIO/dataio_silo_utility.h"
#endif  // if SILO

#include "grid/setup_grid_NG_MPI.h"
#include "grid/uniform_grid.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;

//#define TESTING

// ##################################################################
// ##################################################################

int main(int argc, char **argv)
{
  int err = COMM->init(&argc, &argv);
  if (err) rep.error("comms init error", err);

  int r = -1, np = -1;
  COMM->get_rank_nproc(&r, &np);
  class SimParams SimPM;

  if (argc < 2) {
    cerr << "Error, please give a filename to read IC parameters from.\n";
    cerr << "Usage <icgen_NG_serial> <paramfile> [ic-filetype]\n";
    return (1);
  }

  string *args = 0;
  args         = new string[argc];
  for (int i = 0; i < argc; i++)
    args[i] = argv[i];
  // Redirect stdout/stderr if required.
  for (int i = 0; i < argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      ostringstream path;
      path << outpath << "_" << r << "_";
      outpath = path.str();
      if (r == 0) {
        cout << "Redirecting stdout to " << outpath << "info.txt"
             << "\n";
      }
      rep.redirect(outpath);  // Redirects cout and cerr to text file.
    }
  }
#ifndef TESTING
  rep.kill_stdout_from_other_procs(0);
#endif

  class DataIOBase *dataio          = 0;
  class get_sim_info *siminfo       = 0;
  class ICsetup_base *ic            = 0;
  class ReadParams *rp              = 0;
  class setup_grid_NG_MPI *SimSetup = 0;
  MP                                = 0;  // global microphysics class pointer.

  SimSetup = new setup_grid_NG_MPI();

  string pfile = argv[1];
  string icftype;
  if (argc > 2)
    icftype = argv[2];
  else
    icftype = "silo";  // This is the default for now.

  siminfo = 0;
  siminfo = new class get_sim_info();
  if (!siminfo) rep.error("Sim Info class init error", siminfo);
  err = 0;
  err += siminfo->read_gridparams(pfile, SimPM);
  if (err) rep.error("Read Grid Params Error", err);
  delete siminfo;
  siminfo = 0;

  // cout <<"ICGEN_NG_MPI setting up grid levels\n";
  SimSetup->setup_NG_grid_levels(SimPM);
  cout << "ICGEN_NG_MPI: grid levels set up.\n";
  vector<class GridBaseClass *> grid;
  grid.resize(SimPM.grid_nlevels);

  //
  // Set up the grids.
  //
  err      = SimSetup->setup_grid(grid, SimPM);
  SimPM.dx = grid[0]->DX();
  if (!grid[0]) rep.error("Grid setup failed", grid[0]);

  //
  // read in what kind of ICs we are setting up.
  //
  rp = new ReadParams;
  if (!rp) rep.error("icgen:: initialising RP", rp);
  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);
  string seek = "ics";
  string ics  = rp->find_parameter(seek);
  setup_ics_type(ics, &ic);
  ic->set_SimPM(&SimPM);
  ic->set_MCMD_pointer(&(SimPM.levels[0].MCMD));

  err = SimSetup->set_equations(SimPM);
  rep.errorTest("(icgen::set_equations) err!=0 Fix me!", 0, err);
  class FV_solver_base *solver = SimSetup->get_solver_ptr();

  cout << "ICGEN_NG_MPI: setting up microphysics module\n";
  SimSetup->setup_microphysics(SimPM);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // call "setup" to set up the data on the computational grid.
  cout << "ICGEN_NG_MPI: putting data on grids\n";
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // cout <<"ICGEN_NG_MPI: putting data on grids\n";
    err += ic->setup_data(rp, grid[l]);
    if (err) rep.error("Initial conditions setup failed.", err);
  }

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    // Set Ph=P in every cell.
    cell *c = grid[l]->FirstPt();
    // rep.printVec("First Point",c->P,SimPM.nvar);
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

  //
  // Set up the boundary conditions, since internal boundary data
  // should be already set to its correct value in the initial
  // conditions file.
  //
  //#ifdef TESTING
  cout << "ICGEN_NG_MPI: Setting up boundaries\n";
  //#endif
  SimSetup->boundary_conditions(SimPM, grid);
  if (err) rep.error("icgen: Couldn't set up boundaries.", err);

  //#ifdef TESTING
  cout << "ICGEN_NG_MPI: Setting up raytracing\n";
  //#endif
  err += SimSetup->setup_raytracing(SimPM, grid);
  if (err) rep.error("icgen: Failed to setup raytracer", err);

  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout << "icgen_NG_MPI: assigning boundary data for level " << l << "\n";
#endif
    err = SimSetup->assign_boundary_data(SimPM, l, grid[l]);
    COMM->barrier("level assign boundary data");
    rep.errorTest("icgen_NG_MPI::assign_boundary_data", 0, err);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  //#ifdef TESTING
  cout << "icgen_NG_MPI: updating boundary data\n";
//#endif
#ifdef TESTING
  cout << "icgen_NG_MPI: updating external boundaries\n";
#endif
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err += SimSetup->TimeUpdateExternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  rep.errorTest("sim_init_NG: error from boundary update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
#ifdef TESTING
  cout << "icgen_NG_MPI: updating C2F boundaries\n";
#endif
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout << "NG_MPI updating C2F boundaries for level " << l << "\n";
    cout << "@@@@@@@@@@@@  UPDATING C2F BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    if (l < SimPM.grid_nlevels - 1) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_SEND) {
          err += SimSetup->BC_update_COARSE_TO_FINE_SEND(
              SimPM, grid[l], solver, l, grid[l]->BC_bd[i], 2, 2);
        }
      }
    }
#ifdef TESTING
    cout << "icgen_NG_MPI: updating C2F boundaries, sent, now recv\n";
#endif
    if (l > 0) {
      for (size_t i = 0; i < grid[l]->BC_bd.size(); i++) {
        if (grid[l]->BC_bd[i]->itype == COARSE_TO_FINE_RECV) {
          err += SimSetup->BC_update_COARSE_TO_FINE_RECV(
              SimPM, solver, l, grid[l]->BC_bd[i], SimPM.levels[l].step);
        }
      }
    }
#ifdef TESTING
    cout << "icgen_NG_MPI: updating C2F boundaries: done with level\n";
#endif
  }
  cout << "icgen_NG_MPI: updating C2F boundaries done\n";
  SimSetup->BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from boundary update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
#ifdef TESTING
  cout << "icgen_NG_MPI: updating external boundaries\n";
#endif
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    err += SimSetup->TimeUpdateExternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  // SimSetup->BC_COARSE_TO_FINE_SEND_clear_sends();
  rep.errorTest("sim_init_NG: error from boundary update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
#ifdef TESTING
  cout << "icgen_NG_MPI: updating internal boundaries\n";
#endif
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TESTING
    cout << "@@@@@@@@@@@@  UPDATING INTERNAL BOUNDARIES FOR LEVEL ";
    cout << l << "\n";
#endif
    err += SimSetup->TimeUpdateInternalBCs(
        SimPM, l, grid[l], solver, SimPM.simtime, SimPM.tmOOA, SimPM.tmOOA);
  }
  rep.errorTest("sim_init_NG: error from boundary update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // update fine-to-coarse level boundaries
#ifdef TESTING
  cout << "icgen_NG_MPI: updating F2C boundaries\n";
#endif
  for (int l = SimPM.grid_nlevels - 1; l >= 0; l--) {
#ifdef TESTING
    cout << "NG_MPI updating F2C boundaries for level " << l << "\n";
    cout << l << "\n";
#endif
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
  SimSetup->BC_FINE_TO_COARSE_SEND_clear_sends();
  rep.errorTest("NG_MPI INIT: error from boundary update", 0, err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // if data initialised ok, maybe we need to equilibrate the
  // chemistry...
  //
  if (SimPM.ntracer > 0 && (SimPM.EP.chemistry)) {
    cout << "MAIN: equilibrating the chemical species.\n";
    if (!MP) rep.error("microphysics init", MP);

    // first avoid cooling the gas in getting to equilbrium, by
    // setting update_erg to false.
    bool uerg           = SimPM.EP.update_erg;
    SimPM.EP.update_erg = false;
    err                 = ic->equilibrate_MP(grid[0], MP, rp, SimPM);
    if (err) rep.error("setting chemical states to equilibrium failed", err);

    SimPM.EP.update_erg = uerg;
#ifdef TESTING
    cout << "MAIN: finished equilibrating the chemical species.\n";
#endif
  }
  // ----------------------------------------------------------------

  cout << "IC file-type is " << icftype << "\n";
  seek          = "OutputFile";
  string icfile = rp->find_parameter(seek);
  if (icfile == "") {
    cout << "WARNING: no filename for ic file.  writing to IC_temp.****\n";
    icfile = "IC_temp";
  }

  dataio = 0;  // zero the class pointer.
  if (icftype == "text") {
    cout << "WRITING ASCII TEXT FILE: ";
    cout << icfile << "\n";
  }

#ifdef SILO
  if (icftype == "silo") {
    cout << "WRITING SILO FILE: ";
    //    icfile = icfile+".silo";
    cout << icfile << "\n";
    dataio = 0;
    dataio = new dataio_silo_utility(SimPM, "DOUBLE", &(SimPM.levels[0].MCMD));
  }
#endif  // if SILO defined.
  if (!dataio) rep.error("IO class initialisation: ", icftype);
  err = dataio->OutputData(icfile, grid, SimPM, 0);
  if (err) rep.error("File write error", err);
  delete dataio;
  dataio = 0;
  cout << icftype << " FILE WRITTEN\n";

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

  for (unsigned int v = 0; v < grid.size(); v++) {
    class GridBaseClass *g = grid[v];
    delete g;
  }

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
  cout << "rank: " << SimPM.levels[0].MCMD.get_myrank();
  cout << " nproc: " << SimPM.levels[0].MCMD.get_nproc() << "\n";
  COMM->finalise();
  delete COMM;
  COMM = 0;
  return err;
}

// ##################################################################
// ##################################################################
