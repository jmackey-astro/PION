/// \file sim_init.cpp
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for sim_init, which sets up a PION simulation
/// and gets everything ready to run.
///
/// Modifications:\n
/// - 2018.05.11 JM: moved code from sim_control.cpp
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "sim_control/sim_init.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"
#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif  // if FITS

#include <climits>
#include <sstream>
using namespace std;

// ##################################################################
// ##################################################################

sim_init::sim_init()
{
#ifndef NDEBUG
  spdlog::debug("(sim_init::Constructor)");
#endif
  SimPM.checkpoint_freq = INT_MAX;
  max_walltime          = 1.0e100;
  return;
}

// ##################################################################
// ##################################################################

sim_init::~sim_init()
{
#ifndef NDEBUG
  spdlog::debug("(sim_init::Destructor)");
#endif
  if (dataio) {
    delete dataio;
    dataio = 0;
  }
  if (textio) {
    delete textio;
    textio = 0;
  }
  return;
}

// ##################################################################
// ##################################################################

double sim_init::get_max_walltime()
{
  return max_walltime;
}

// ##################################################################
// ##################################################################

void sim_init::set_max_walltime(double t  ///< New Max. runtime in seconds.
)
{
  spdlog::debug("\tResetting max. walltime from {}", max_walltime);
  max_walltime = t;
  spdlog::debug(" to new value: {}", max_walltime);
}

// ##################################################################
// ##################################################################

//---------------------------------------------------------
//
// Function to output commandline options for code.
//
void sim_init::print_command_line_options(int argc, char **argv)
{
  spdlog::info("PION: You ran:");
  std::stringstream invocation;
  for (int v = 0; v < argc; v++)
    invocation << argv[v] << " ";
  spdlog::info(invocation.get());
  spdlog::info("      ************************         ");
  spdlog::info("{}: must call with at least 1 argument...", argv[0]);
  spdlog::info(" <main> <icfile> [optional args]");
  spdlog::info("Parameters:");
  spdlog::info("<icfile> ");
  spdlog::info("\tCan be an ASCII parameter-file for 1D and 2D shocktubes.");
  spdlog::info("\tOtherwise should be a restart-file in FITS or Silo format.");
  spdlog::info(
      "[optional args] are in the format <name>=<value> with no spaces.");
  spdlog::info("*********** DATA I/O OPTIONS ************");
  spdlog::info(
      "\t redirect=string : filename with path to redirect stdout/stderr to");
  spdlog::info(
      "\t op_criterion=N  : 0=output every I steps, 1=output every D time units.");
  spdlog::info(
      "\t opfreq=N        : Save snapshot every Nth timestep  (if op_criterion=0).");
  spdlog::info(
      "\t opfreq_time=D   : Save snapshot every Dth time unit (if op_criterion=1).");
  spdlog::info(
      "\t finishtime=D    : set time to finish simulation, in code time units.");
  spdlog::info(
      "\t optype=S        : Specify type of output file, [1,text]=TEXT,[2,fits]=FITS,[4,both]=FITS+TEXT,[5,silo]=SILO,[6]=SILO+TEXT.");
  spdlog::info(
      "\t outfile=NAME    : Replacement snapshot filename, with path.");
  spdlog::info("*********** PHYSICS/Grid OPTIONS *************");
  spdlog::info("\t ooa=N         : modify order of accuracy (either 1 or 2).");
  spdlog::info(
      "\t AVtype=N      : modify type of artificial viscosity: 0=none, 1=Falle,Komissarov,Joarder(1998), 3=Sanders et al.(1998)[H-correction], 4=both 1+3.");
  spdlog::info(
      "\t EtaVisc=D     : modify viscosity parameter to the given double precision value.");
  spdlog::info(
      "\t cfl=D         : change the CFL no. for the simulation, in range (0,1).");
  spdlog::info(
      "\t cooling=N     : cooling=0 for no cooling, >0 for different prescriptions.\t solver=N      :");
  spdlog::info("\t\t 0 = Lax-Friedrichs Flux");
  spdlog::info(
      "\t\t 1 = Linear Riemann Solver : HD/MHD (Falle, Komissarov, Joarder, 1998),");
  spdlog::info(
      "\t\t 2 = Exact Riemann Solver  : HD (Hirsch (199X), Toro, 1999)");
  spdlog::info("\t\t 3 = Hybrid Riemann Solver (1+2)         : HD ");
  spdlog::info(
      "\t\t 4 = Roe Conserved Variables flux solver : HD/MHD (e.g. Toro, 1999, Stone, Gardiner et al. 2008)");
  spdlog::info(
      "\t\t 5 = Roe Primitive Variables flux solver : HD (e.g. Stone, Gardiner et al. 2008)");
  spdlog::info("\t\t 6 = Flux vector splitting : HD only (van Leer, 1982) ");
  spdlog::info("\t\t 7 = HLLD solver : MHD only");
  spdlog::info("\t\t 8 = HLL  solver : HD/MHD ");
  spdlog::info("*********** PARALLEL CODE ONLY *************");
  spdlog::info("\t maxwalltime=D : change the max. runtime to D in hours.");
  spdlog::info("*********** NESTED GRID CODE ONLY *************");
  spdlog::info("\t nlevels=N     : modify number of levels in NG grid.");
  spdlog::info(
      "\t wind_radius_N=D : change radius of boundary for wind source N to value D (in cm)");
  spdlog::info("     *********************************************");
}



// ##################################################################
// ##################################################################



int sim_init::Init(
    string infile,
    int typeOfFile,
    int narg,
    string *args,
    vector<class GridBaseClass *>
        &grid  ///< address of vector of grid pointers.
)
{
  spdlog::info("(pion)  Initialising");
  int err = 0;

#ifdef SERIAL
  SimPM.typeofip = typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  err = dataio->ReadHeader(infile, SimPM);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "(INIT::get_parameters) err!=0 Something went wrong", 0, err);
    exit_pion(1);
  }
#endif  // SERIAL

  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "(INIT::override_params) err!=0 Something went wrong", 0, err);
    exit_pion(1);
  }

  // Now set up the grid structure.
  grid.resize(1);
  err      = setup_grid(grid, SimPM);
  SimPM.dx = grid[0]->DX();
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "(INIT::setup_grid) err!=0 Something went wrong", 0, err);
    exit_pion(1);
  }

  //
  // All grid parameters are now set, so I can set up the appropriate
  // equations/solver class.
  //
  err = set_equations(SimPM);
  if (err) {
    spdlog::error("(INIT::set_equations) err!=0 Fix me! {} {}", 0, err);
    exit_pion(1);
  }
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->SetEOS(SimPM.gamma);
#ifdef PION_OMP
  }
#endif

  //
  // Now setup Microphysics, if needed.
  //
  err = setup_microphysics(SimPM);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::setup_microphysics) err!=0", 0,
        err);
    exit_pion(1);
  }
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->SetMicrophysics(MP);
#ifdef PION_OMP
  }
#endif

  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, grid, SimPM);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "(INIT::assign_initial_data) err!=0 Something went wrong", 0, err);
    exit_pion(1);
  }

  //
  // Set Ph[] = P[], and then implement the boundary conditions.
  //
  cell *c = grid[0]->FirstPt();
  do {
    for (int v = 0; v < SimPM.nvar; v++)
      c->Ph[v] = c->P[v];
  } while ((c = grid[0]->NextPt(*c)) != 0);

  //
  // If I'm using the GLM method, make sure Psi variable is
  // initialised to zero.
  //
  if (SimPM.eqntype == EQGLM && SimPM.timestep == 0) {
#ifndef NDEBUG
    spdlog::debug("Initial state, zero-ing glm variable.");
#endif
    c = grid[0]->FirstPt();
    do {
      c->P[SI] = c->Ph[SI] = 0.;
    } while ((c = grid[0]->NextPt(*c)) != 0);
  }

  //
  // Assign boundary conditions to boundary points.
  //
  err = boundary_conditions(SimPM, grid);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::boundary_conditions) err!=0", 0,
        err);
    exit_pion(1);
  }
  err = assign_boundary_data(SimPM, 0, grid[0], MP);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}", "(INIT::assign_boundary_data) err!=0", 0,
        err);
    exit_pion(1);
  }

  //
  // Setup Raytracing on each grid, if needed.
  //
  err += setup_raytracing(SimPM, grid[0]);
  err += setup_evolving_RT_sources(SimPM);
  err += update_evolving_RT_sources(SimPM, SimPM.simtime, grid[0]->RT);
  if (err) {
    spdlog::error(
        "{}: Expected {} but got {}",
        "Failed to setup raytracer and/or microphysics", 0, err);
    exit_pion(1);
  }

  //
  // If testing the code, this calculates the momentum and energy on the
  // domain.
  //
  initial_conserved_quantities(grid[0]);

  err += TimeUpdateInternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, 0.0, SimPM.tmOOA,
      SimPM.tmOOA);
  err += TimeUpdateExternalBCs(
      SimPM, 0, grid[0], spatial_solver, SimPM.simtime, SimPM.tmOOA,
      SimPM.tmOOA);
  if (err) {
    spdlog::error(
        "{}: {}", "first_order_update: error from bounday update", err);
    exit_pion(1);
  }

  //
  // If using opfreq_time, set the next output time correctly.
  //
  if (SimPM.op_criterion == 1) {
    if (SimPM.opfreq_time < TINYVALUE) {
      spdlog::error("{}: {}", "opfreq_time not set right!", SimPM.opfreq_time);
      exit_pion(1);
    }
    SimPM.next_optime = SimPM.simtime + SimPM.opfreq_time;
    double tmp        = ((SimPM.simtime / SimPM.opfreq_time)
                  - floor(SimPM.simtime / SimPM.opfreq_time))
                 * SimPM.opfreq_time;
    SimPM.next_optime -= tmp;
  }

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {
      delete dataio;
      dataio = 0;
    }
    if (textio) {
      delete textio;
      textio = 0;
    }
    setup_dataio_class(SimPM, SimPM.typeofop);
    if (!dataio) {
      spdlog::error("{}: {}", "INIT:: dataio initialisation", SimPM.typeofop);
      exit_pion(1);
    }
  }
  dataio->SetSolver(spatial_solver);
  dataio->SetMicrophysics(MP);
  if (textio) {
    textio->SetSolver(spatial_solver);
    textio->SetMicrophysics(MP);
  }

#ifdef SERIAL
  if (SimPM.timestep == 0) {
    spdlog::info("(INIT) Writing initial data.");
    err = output_data(grid);
    if (err) {
      spdlog::error("Failed to write file! maybe dir does not exist?");
      exit_pion(1);
    }
  }
  spdlog::info("-------------------------------------------------------");
#endif  // SERIAL

#ifndef NDEBUG
  c = (grid[0])->FirstPt_All();
  do {
    if (pconst.equalD(c->P[RO], 0.0)) {
      cout << "zero data in cell: ";
      CI.print_cell(*c);
    }
  } while ((c = (grid[0])->NextPt_All(*c)) != 0);
#endif  // NDEBUG

  return (0);
}

// ##################################################################
// ##################################################################

int sim_init::override_params(int narg, string *args)
{

  spdlog::info("(pion)  Overriding parameters if requested...");

  // Find command line params and change them
  for (int i = 2; i < narg; i++) {

    if (args[i].find("ooa=") != string::npos) {
      // Assign order of accuracy;  string is 'ooa=N', where N=1 or 2.
      int tmp     = SimPM.spOOA;
      SimPM.spOOA = atoi((args[i].substr(4)).c_str());
      SimPM.tmOOA = SimPM.spOOA;
      spdlog::info(
          "OVERRIDE PARAMS: Resetting OOA from ooa={} to command-line value = {}",
          tmp, SimPM.spOOA);
    }

    else if (args[i].find("nlevels=") != string::npos) {
      // Assign number of grid levels;  string is 'nlevels=N', where N>=1.
      int tmp            = SimPM.grid_nlevels;
      SimPM.grid_nlevels = atoi((args[i].substr(8)).c_str());
      spdlog::info(
          "OVERRIDE PARAMS: Resetting nlevels from nlevels={} to command-line value = {}",
          tmp, SimPM.grid_nlevels);
    }

    else if (args[i].find("AVtype=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: old AV={} ... overriding!\n", SimPM.artviscosity);
      // Assign art.viscosity parameter. String is 'artvisc=I' with I in
      // [0,N].
      int v = atoi((args[i].substr(7)).c_str());
      if (v == 0) {
        spdlog::info("Not using artificial viscosity.");
        SimPM.artviscosity = 0;
        SimPM.etav         = 0.;
      }
      else if (v == 1) {
        spdlog::info(
            "Using Falle, Komissarov, Joarder (1998) AV prescription.");
        SimPM.artviscosity = AV_FKJ98_1D;  // ==1
        SimPM.etav         = 0.1;
      }
      else if (v == 2) {
        spdlog::error(
            "{}: {}", "divv viscosity not working, use AVtype=1 or 3", v);
        exit_pion(1);
        spdlog::info(
            "Using Colella and Woodward (1984) AV prescription (Lapidus).****** WARNING, THIS NEEDS TESTING, EXPERIMENTAL CODE )****");
        SimPM.artviscosity = AV_LAPIDUS;  // ==2 (NEEDS TESTING!!!)
        SimPM.etav         = 0.1;
      }
      else if (v == 3) {
        spdlog::info(
            "Using the H-correction of Sanders et al. (1998,JCP,145,511).");
        SimPM.artviscosity = AV_HCORRECTION;
        SimPM.etav = 0.1;  // This parameter is redundant for the H-correction.
      }
      else if (v == 4) {
        spdlog::info(
            "Using the H-correction of Sanders et al. (1998,JCP,145,511) with the 1D viscosity of Falle, Komissarov, Joarder (1998)");
        SimPM.artviscosity = AV_HCORR_FKJ98;  // ==4 (NEEDS TESTING!!!)
        SimPM.etav         = 0.1;
      }
      else if (v == AV_VonNeuRicht) {
        spdlog::error(
            "{}: {}", "von Neumann & Richtmeyer viscosity not working", v);
        exit_pion(1);
        // AV_VonNeuRicht==5
        spdlog::info(
            "Using Multi-D von Neumann & Richtmeyer (1950) viscosity. See Tscharnuter & Winkler (1979), Stone & Norman (1992).");
        spdlog::info("WARNING -- THIS ONLY WORKS WITH EQNTYPE==9(EQEUL_EINT).");
        SimPM.artviscosity = AV_VonNeuRicht;
        SimPM.etav         = 1.0;
      }
      else {
        spdlog::info(
            "DIDN'T UNDERSTAND AV={}, SETTING TO FALLE et al (1998).\n", v);
        SimPM.artviscosity = 1;
        SimPM.etav         = 0.1;
        spdlog::error("{}: {}", "Bad viscosity flag from command-line", v);
        exit_pion(1);
      }
      spdlog::info(
          "OVERRIDE PARAMS: setting AV = {} and eta = {}", SimPM.artviscosity,
          SimPM.etav);
    }

    else if (args[i].find("EtaVisc=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: old and eta={} ... overriding!\n", SimPM.etav);
      // Assign art.viscosity parameter. String is 'artvisc=D' with D in
      // [0,N].
      double visc = atof((args[i].substr(8)).c_str());
      spdlog::info(
          "OVERRIDE PARAMS: Resetting eta_visc from {} to {}", SimPM.etav,
          visc);
      if (visc < 0.0 || !isfinite(visc)) {
        spdlog::error(
            "{}: {}", "Error: eta viscosity parameter outside allowed range!",
            visc);
        exit_pion(1);
      }
      SimPM.etav = visc;
    }

    else if (args[i].find("opfreq=") != string::npos) {
      // Assign output frequency to new value. String is 'opfreq=N' with
      // N=[0..Nmax].
      int tmp = atoi((args[i].substr(7)).c_str());
      spdlog::info(
          "OVERRIDE PARAMS: Resetting opfreq from {} to {}", SimPM.opfreq, tmp);
      SimPM.op_criterion = 0;
      SimPM.opfreq       = tmp;
    }
    else if (args[i].find("outfile=") != string::npos) {
      // Assign a new output file.  string is outfile=char[128max]
      string tmp = args[i].substr(8);
      spdlog::info(
          "OVERRIDE PARAMS: Resetting output file from {} to new name: {}.xxx.ftype\n",
          SimPM.outFileBase, tmp);
      SimPM.outFileBase = tmp;
      tmp.clear();
    }

    else if (args[i].find("optype=") != string::npos) {
      // assign new op-type; 1=TXT,2=FITS,3=FitsTable,4=TXT+FITS, 5=SILO
      string now;
      if (SimPM.typeofop == 1)
        now = "text";
      else if (SimPM.typeofop == 2)
        now = "fits";
      else if (SimPM.typeofop == 3)
        now = "ftab";
      else if (SimPM.typeofop == 4)
        now = "f+tt";
      else if (SimPM.typeofop == 5)
        now = "silo";
      else if (SimPM.typeofop == 6)
        now = "silo+text";
      else {
        spdlog::error("{}: {}", "What kind of output is this?", SimPM.typeofop);
        exit_pion(1);
      }

      string chg = args[i].substr(7);
      int tmp    = -1;
      if (chg == "text" || chg == "txt" || chg == "TEXT" || chg == "TXT"
          || chg == "1")
        tmp = 1;
      else if (chg == "fits" || chg == "FITS" || chg == "2")
        tmp = 2;
      else if (chg == "ftab" || chg == "FTAB" || chg == "3")
        tmp = 3;
      else if (chg == "txtfits" || chg == "both" || chg == "BOTH" || chg == "4")
        tmp = 4;
      else if (chg == "silo" || chg == "SILO" || chg == "5")
        tmp = 5;
      else if (chg == "txtsilo" || (chg == "textsilo") || (chg == "6"))
        tmp = 6;
      else {
        spdlog::error("{}: {}", "What kind of output do you want?", chg);
        exit_pion(1);
      }
      spdlog::info(
          "OVERRIDE PARAMS: Resetting output file type from {} to {}", now,
          chg);
      SimPM.typeofop = tmp;
    }

    else if (args[i].find("noise=") != string::npos) {
      // assign new value to addnoise, double frac = fractional noise
      // level
      double tmp = atof((args[i].substr(6)).c_str());
      spdlog::info(
          "OVERRIDE PARAMS: Resetting addnoise value from {} to {}", tmp,
          SimPM.addnoise);
      SimPM.addnoise = tmp;
    }

    else if (args[i].find("finishtime=") != string::npos) {
      // assign new value to finishtime.
      double tmp = atof((args[i].substr(11)).c_str());
      spdlog::info(
          "OVERRIDE PARAMS: Resetting finishtime value from {} to {}",
          SimPM.finishtime, tmp);
      SimPM.finishtime = tmp;
    }

    else if (args[i].find("redirect=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: already redirecting stdout, continueing...");
    }

    else if (args[i].find("maxwalltime=") != string::npos) {
      // this is already handled by mainMPI.cc, and ignored for serial
      // code.
    }

    else if (args[i].find("omp-nthreads=") != string::npos) {
      // this is already handled by mainMPI.cc, and ignored for serial
      // code.
    }

    else if (args[i].find("coordsys=") != string::npos) {
      // Change the coordinate system!
      spdlog::info(
          "OVERRIDE PARAMS: Resetting the coordinate system from value={}",
          SimPM.coord_sys);
      string t = args[i].substr(9);
      if (t == "cartesian" || t == "cart" || t == "crt") {
        spdlog::info(" to cartesian coords.");
        SimPM.coord_sys = COORD_CRT;
      }
      else if (t == "cylindrical" || t == "cyl") {
        spdlog::info(" to cylindrical coords.");
        SimPM.coord_sys = COORD_CYL;
      }
      else if (t == "spherical" || t == "sph") {
        spdlog::info(" to spherical coords.");
        SimPM.coord_sys = COORD_SPH;
      }
      else {
        spdlog::error("{}: {}", "don't know this coordinate system", t);
        exit_pion(1);
      }
      spdlog::warn("THIS IS DANGEROUS!!!");
    }

    else if (args[i].find("cfl=") != string::npos) {
      // Assign cfl no., where string is 'cfl=0.X'
      spdlog::info(
          "OVERRIDE PARAMS: Resetting CFL from original value of {}",
          SimPM.CFL);
      double c = atof((args[i].substr(4)).c_str());
      if (c < 0. || !isfinite(c)) {
        spdlog::error("Bad CFL no. {}", c);
        exit_pion(1);
      }
      else if (c > 1.)
        spdlog::warn("WARNING: CFL no. >1, so results will be unstable!!!");
      SimPM.CFL = c;
      spdlog::info(" to command-line value = {}", SimPM.CFL);
    }

    else if (args[i].find("gamma=") != string::npos) {
      // Assign new value to EOS gamma, where string is 'gamma=X.XXXXX'
      spdlog::info(
          "OVERRIDE PARAMS: Resetting EOS gamma from original value of {}",
          SimPM.gamma);
      double c = atof((args[i].substr(6)).c_str());
      if (c <= 1.0 || !isfinite(c)) {
        spdlog::error("{}: {}", "Bad EOS gamma no.", c);
        exit_pion(1);
      }
      else if (c > 2.0)
        spdlog::warn("WARNING: gamma >2 ?");
      SimPM.gamma = c;
      spdlog::info(" to command-line value = {}", SimPM.gamma);
    }

    else if (args[i].find("cooling=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting cooling");
      int c = atoi((args[i].substr(8)).c_str());
      if (c < 0 || c > 100) {
        spdlog::error("{}: {}", "Bad cooling flag (only 0-11 allowed", c);
        exit_pion(1);
      }
      spdlog::info(" flag from {}", SimPM.EP.cooling);
      SimPM.EP.cooling = c;
      spdlog::info(" to {}", SimPM.EP.cooling);
    }

    else if (args[i].find("dynamics=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting dynamics");
      int c = atoi((args[i].substr(9)).c_str());
      if (c < 0 || c > 1) {
        spdlog::error("{}: {}", "Bad dynamics flag (only 0,1, allowed", c);
        exit_pion(1);
      }
      spdlog::info(" flag from {}", SimPM.EP.dynamics);
      SimPM.EP.dynamics = c;
      spdlog::info(" to {}", SimPM.EP.dynamics);
    }

    else if (args[i].find("simulation_time=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting simulation time");
      double c = atof((args[i].substr(16)).c_str());
      if (c < 0 || !isfinite(c)) {
        spdlog::error("Bad time value: {}", c);
        exit_pion(1);
      }
      spdlog::info("from {}", SimPM.simtime);
      SimPM.simtime = c;
      spdlog::info(" to {}", SimPM.simtime);
    }

    else if (args[i].find("raytracing=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting raytracing");
      int c = atoi((args[i].substr(11)).c_str());
      if (c < 0 || c > 1) {
        spdlog::error("{}: {}", "Bad raytracing flag (only 0,1, allowed", c);
        exit_pion(1);
      }
      spdlog::info(" flag from {}", SimPM.EP.raytracing);
      SimPM.EP.raytracing = c;
      spdlog::info(" to {}", SimPM.EP.raytracing);
    }

    else if (args[i].find("chemistry=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting chemistry");
      int c = atoi((args[i].substr(10)).c_str());
      if (c < 0 || c > 1) {
        spdlog::error("{}: {}", "Bad chemistry flag (only 0,1, allowed", c);
        exit_pion(1);
      }
      spdlog::info(" flag from {}", SimPM.EP.chemistry);
      SimPM.EP.chemistry = c;
      spdlog::info(" to {}", SimPM.EP.chemistry);
    }

    else if (args[i].find("microphysics=") != string::npos) {
      spdlog::info("OVERRIDE PARAMS: resetting microphysics");
      string t = (args[i].substr(13));
      spdlog::info(" flag from {}", SimPM.chem_code);
      SimPM.chem_code = t;
      spdlog::info(" to {}", SimPM.chem_code);
    }

    else if (args[i].find("solver=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting solver: 0=LF,1=RSlin,2=RSexact,3=RShybrid,4=RSroe,5=RSroePV,6=FVS,7=HLLD,8=HLL,9=Hydro-hybrid:");
      int c = atoi((args[i].substr(7)).c_str());
      if (c < 0 || c > 9) {
        spdlog::error("Bad solver flag (only 0,1,2,3,4,5,6,7,8,9 allowed", c);
        exit_pion(1);
      }
      spdlog::info(" solver from {}", SimPM.solverType);
      SimPM.solverType = c;
      spdlog::info(" to {}", SimPM.solverType);
    }

    else if (args[i].find("op_criterion=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting op_criterion from {} to ",
          SimPM.op_criterion);
      int c = atoi((args[i].substr(13)).c_str());
      if (c < 0 || c > 1) {
        spdlog::error("{}: {}", "Bad op_criterion flag:", c);
        exit_pion(1);
      }
      SimPM.op_criterion = c;
      spdlog::info("{}", SimPM.op_criterion);
    }
    else if (args[i].find("opfreq_time=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting opfreq_time from {} units to ",
          SimPM.opfreq_time);
      double c = atof((args[i].substr(12)).c_str());
      if (c < 0.0 || c > 1.e50) {
        spdlog::error("{}: {}", "Bad opfreq_time flag:", c);
        exit_pion(1);
      }
      SimPM.op_criterion = 1;
      SimPM.opfreq_time  = c;
      spdlog::info("{} units.", SimPM.opfreq_time);
      SimPM.next_optime = SimPM.simtime + SimPM.opfreq_time;
      if (SimPM.timestep > 0) {
        double tmp = ((SimPM.simtime / SimPM.opfreq_time)
                      - static_cast<int>(SimPM.simtime / SimPM.opfreq_time))
                     * SimPM.opfreq_time;
        SimPM.next_optime -= tmp;
      }
    }

    else if (args[i].find("min_timestep=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting min_timestep from {} units [NOT YEARS!!!] to ",
          SimPM.min_timestep);
      double c = atof((args[i].substr(13)).c_str());
      if (c < 0.0 || c > 1.e50) {
        spdlog::error("{}: {}", "Bad min_timestep flag:", c);
        exit_pion(1);
      }
      SimPM.min_timestep = c;
      spdlog::info("{} units.", SimPM.min_timestep);
    }

    else if (args[i].find("limit_timestep=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: limiting timestep: changing from {} to ",
          SimPM.EP.MP_timestep_limit);
      int c = atoi((args[i].substr(15)).c_str());
      if (c < 0 || c > 5) {
        spdlog::error("{}: {}", "Bad MP_timestep_limit flag:", c);
        exit_pion(1);
      }
      SimPM.EP.MP_timestep_limit = c;
      spdlog::info("{}", SimPM.EP.MP_timestep_limit);
    }

    else if (args[i].find("checkpt_freq=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: checkpointing freq.: changing from {} to ",
          SimPM.checkpoint_freq);
      int c = atoi((args[i].substr(13)).c_str());
      if (c < 0) {
        spdlog::error("{}: {}", "Bad checkpoint_freq flag:", c);
        exit_pion(1);
      }
      SimPM.checkpoint_freq = c;
      spdlog::info("{}", SimPM.checkpoint_freq);
    }

    else if (args[i].find("max_T=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting MaxTemperature from {} K to ",
          SimPM.EP.MaxTemperature);
      double c = atof((args[i].substr(6)).c_str());
      if (c < 0.0 || c > 1.e50) {
        spdlog::error("{}: {}", "Bad Max_T flag:", c);
        exit_pion(1);
      }
      SimPM.EP.MaxTemperature = c;
      spdlog::info("{} K.", SimPM.EP.MaxTemperature);
    }

    else if (args[i].find("min_T=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting MinTemperature from {} K to ",
          SimPM.EP.MinTemperature);
      double c = atof((args[i].substr(6)).c_str());
      if (c < 0.0 || c > 1.e50) {
        spdlog::error("{}: {}", "Bad Min_T flag:", c);
        exit_pion(1);
      }
      SimPM.EP.MinTemperature = c;
      spdlog::info("{} K.", SimPM.EP.MinTemperature);
    }

    else if (args[i].find("wind_radius=") != string::npos) {
      if (SWP.Nsources < 1) {
        spdlog::error("reset wind radius without a source {}", SWP.Nsources);
        exit_pion(1);
      }
      string q = (args[i].substr(11, 1));
      if (q == "=") {
        spdlog::error("must ID wind source, e.g. wind_radius_0 {}", args[i]);
        exit_pion(1);
      }
      int src = atoi((args[i].substr(12)).c_str());
      if (src < 0 || src > 9 || !isfinite(src)) {
        spdlog::error("expect format wind_radius_0=1.2e17 {}", src);
        exit_pion(1);
      }
      else if (static_cast<size_t>(src) >= SWP.params.size()) {
        spdlog::error(
            "change wind radius for source that doesn't exist {}", src);
        exit_pion(1);
      }
      spdlog::info(
          "OVERRIDE PARAMS: resetting radius of wind src {} from {} cm to ",
          src, SWP.params[src]->radius);
      double c = atof((args[i].substr(14)).c_str());
      if (c < 0.0 || c > 1.e50) spdlog::error("{}: {}", "Bad radius flag:", c);
      SWP.params[src]->radius = c;
      spdlog::info("{} cm.", SWP.params[src]->radius);
    }

    else if (args[i].find("tc_strength=") != string::npos) {
      spdlog::info(
          "OVERRIDE PARAMS: resetting EP.tc_strength from {} to ",
          SimPM.EP.tc_strength);
      double c = atof((args[i].substr(12)).c_str());
      if (c < 0.0 || c > 1.01) {
        spdlog::error("{}: {}", "Bad tc_strength flag:", c);
        exit_pion(1);
      }
      SimPM.EP.tc_strength = c;
      spdlog::info("{}", SimPM.EP.tc_strength);
    }


    else {
      spdlog::error(
          "Don't recognise this optional argument, please fix. {}", args[i]);
      exit_pion(i);
    }
  }
  return (0);
}  // sim_init::override_params



// ##################################################################
// ##################################################################



int sim_init::output_data(vector<class GridBaseClass *>
                              &grid  ///< address of vector of grid pointers.
)
{
  int err = 0;
  //
  // Now we move on the various output criteria.  The if statements
  // call return(0) if we decide we don't need to ouput, so if we get
  // past them then just output the current timestep.
  //
  // Always output at the start or the end.
  //
  if (SimPM.timestep == 0 || SimPM.maxtime == true) {
  }
  //
  // If outputting every nth step, see if we are on an output step.
  // If not, return.
  //
  else if (SimPM.op_criterion == 0) {
    if ((SimPM.opfreq == 0) && (SimPM.maxtime == false)
        && (SimPM.timestep != 0))
      return 0;
    // Next check if we are in an outputting timestep.
    else if (
        (SimPM.maxtime == false) && (SimPM.timestep % SimPM.opfreq) != 0
        && (SimPM.timestep != 0))
      return (0);
  }
  //
  // If outputting every D time units, see if we are at an output
  // time, and if not then return.
  //
  else if (SimPM.op_criterion == 1) {
    if (!pconst.equalD(SimPM.simtime, SimPM.next_optime)
        && (SimPM.maxtime == false))
      return 0;
    else
      // we are to output data, so advance the 'next' counter and
      //  continue.
      SimPM.next_optime += SimPM.opfreq_time;
  }
  else {
    spdlog::error("op_criterion must be 0 or 1: {}", SimPM.op_criterion);
    exit_pion(1);
  }

  //
  // Since we got past all that, we are in a timestep that should be
  // saved, so go and do it...
  //
  spdlog::info(
      "Saving data at step {} and sim-time {:12.6e} to file {}", SimPM.timestep,
      SimPM.simtime, SimPM.outFileBase);
  err = dataio->OutputData(SimPM.outFileBase, grid, SimPM, SimPM.timestep);
  if (textio)
    err += textio->OutputData(SimPM.outFileBase, grid, SimPM, SimPM.timestep);
  if (err) {
    spdlog::error("Error writing data");
    exit_pion(1);
  }
  return (0);
}

// ##################################################################
// ##################################################################

int sim_init::initial_conserved_quantities(class GridBaseClass *grid)
{
#ifdef TEST_CONSERVATION
  // Energy, and Linear Momentum in x-direction.
  std::vector<pion_flt> u(SimPM.nvar, 0.0);
  double dx = grid->DX();
  double dv = 0.0;
  initERG   = 0.;
  initMMX = initMMY = initMMZ = 0.;
  initMASS                    = 0.0;
  class cell *cpt             = grid->FirstPt();
  do {
    if (cpt->isdomain) {
      spatial_solver->PtoU(cpt->P.data(), u.data(), SimPM.gamma);
      dv = spatial_solver->CellVolume(*cpt, dx);
      initERG += u[ERG] * dv;
      initMMX += u[MMX] * dv;
      initMMY += u[MMY] * dv;
      initMMZ += u[MMZ] * dv;
      initMASS += u[RHO] * dv;
    }
  } while ((cpt = grid->NextPt(*cpt)) != 0);
  spdlog::info(
      "(sim_init::InitialconservedQuantities) [{}, {}, {}, {}, {}]\n", initERG,
      initMMX, initMMY, initMMZ, initMASS);
#endif  // TEST_CONSERVATION
  return (0);
}  // initial_conserved_quantities()

// ##################################################################
// ##################################################################

int sim_init::RT_all_sources(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< Computational grid.
    const int)
{
  int err = 0;
  if (!grid->RT) return 0;
  //
  // If we have raytracing, we call the ray-tracing routines
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
#ifdef RT_TESTING
    spdlog::debug("calc_raytracing_col_dens: SRC-ID: {}", isrc);
#endif
    err += grid->RT->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      spdlog::debug("isrc={}", isrc);
      spdlog::error("calc_raytracing_col_dens step in returned error {}", err);
      exit_pion(1);
    }  // if error
  }    // loop over sources
  return err;
}

// ##################################################################
// ##################################################################
