///
/// file:    silo2ascii.cpp
/// author:  Jonathan Mackey
/// date:    2013.04.18
///
/// Adapted from templates/HIIregion_KineticE.cpp
///
/// Description: This is a template file for analysing simulation
/// data with N cores, but reading data which may have been written
/// from M!=N processors.  It is modified to read in silo data and
/// write out fits files.  Run with a single core to get a single big
/// fits file.
///
/// Modifications:
///
///  - 2010-01-28 JM: Calculate a bunch of things for the Lim &
///   Mellema (2003) simulations.
/// - 2010-02-03 JM: Made it work for multi-core use, so can run with
///  8, 16, 32 cores etc. (and the number of cores doesn't have to
///  match the number of cores used to write the data i.e. independent
///  of the number of quadmeshes in the silo file).
/// - 2012.05.18 JM: Adapted template for calculating kinetic energy
///    of gas in simulations of HII regions around runaway stars.
/// - 2012.10.17 JM: Added more statistics.
/// - 2012.11.06 JM: Added yet more statistics, and output
///    mass-function files to their own directory.
/// - 2013.03.10 JM: Modified for the StarBench Spitzer test problem.
/// - 2013.04.18 JM: Modified for pure format conversion.
/// - 2013.04.19 JM: Added option to resize the computational domain.
/// - 2014.05.05 JM: Fixed accounting bug in 2D data.
/// - 2014.10.14 JM: Added option to enlarge the domain with LastPt()
/// - 2014.11.28 JM: Added option to only convert any Nth file.
/// - 2016.06.16 JM: Updated for new PION.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "sim_params.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#include "grid/setup_fixed_grid_MPI.h"
#include "grid/uniform_grid.h"
#include "sub_domain/sub_domain.h"
#ifdef PION_NESTED
#include "grid/setup_grid_NG_MPI.h"
#endif /* PION_NESTED */

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"

#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"
#ifndef EXCLUDE_MPV3
#include "microphysics/MPv3.h"
#endif
#ifndef EXCLUDE_MPV4
#include "microphysics/MPv4.h"
#endif
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"


#include "raytracing/raytracer_SC.h"

#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/spdlog.h>

#include <fitsio.h>

#include <fstream>
#include <silo.h>
#include <sstream>
using namespace std;


// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  int err = 0;

  auto max_logfile_size = 1048576 * 5;
  auto max_logfiles     = 3;
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "silo2fits", "silo2fits.log", max_logfile_size, max_logfiles));

#ifdef NDEBUG
  spdlog::set_level(spdlog::level::err);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  //
  // Also initialise the sub_domain class with myrank and nproc.
  // Get nproc from command-line (number of fits files for each
  // snapshot)
  //
  class SimParams SimPM;

  int myrank = SimPM.levels[0].sub_domain.get_myrank();
  int nproc  = SimPM.levels[0].sub_domain.get_nproc();

  if (nproc > 1) spdlog::error("{}: {}", "This is serial code", nproc);

  //
  // Get input files and an output file.
  //
  if (argc < 5) {
    spdlog::error("{}: {}", "Bad number of args", argc);
    spdlog::error(
        "Use as follows:\n<executable-filename> <input-path> <input-silo-file-base> <output-path> <output-fits-file-base> [Nskip] [<xmin> <xmax> [<ymin> <ymax> [<zmin> <zmax>]]]\n******************************************\ninput path:  path to input files.\ninput file:  base filename of sequence of files including _0000 if parallel.\noutput path: directory to write output files to.\noutput file: filename for output FITS file(s).\n[Nskip]:     skip files each iteration.\nxmin/xmax:   optional parameters to resize the domain");
  }

  //
  // The code will get a list of files matching 'input_file' in the
  // directory 'input_path' and do the same operation on all files in
  // thie list.
  //
  string input_path = argv[1];
  string input_file = argv[2];

  //
  // outfile should contain the path as well (relative or absolute)
  //
  string op_path = argv[3];
  string outfile = argv[4];

  //
  // Only convert every Nskip file
  //
  size_t Nskip = 0;
  if (argc > 5) {
    Nskip = static_cast<size_t>(atoi(argv[5]));
  }

  //
  // Redirect output to a text file if you want to:
  //
  ostringstream redir;
  redir.str("");
  redir << op_path << "/msg_" << outfile << "_rank" << myrank << "_";

  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  spdlog::info(
      "-------------------------------------------------------\n--------------- Getting List of Files to read ---------");

  //
  // set up dataio_utility class and fits-writer class.
  //
  class dataio_silo_utility dataio(
      SimPM, "DOUBLE", &(SimPM.levels[0].sub_domain));

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file, &files);
  if (err) spdlog::error("{}: {}", "failed to get list of files", err);
  //
  // remove non-silo files
  //
  for (list<string>::iterator s = files.begin(); s != files.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo") == string::npos) {
      spdlog::debug("removing file {} from list", *s);
      files.erase(s);
      s = files.begin();
    }
    else {
      spdlog::debug("files: {}", *s);
    }
  }
  size_t nfiles = files.size();
  if (nfiles < 1)
    spdlog::error("{}: {}", "Need at least one file, but got none", nfiles);

  spdlog::info(
      "--------------- Got list of Files ---------------------\n-------------------------------------------------------\n--------------- Setting up Grid -----------------------");

  //
  // Set low-memory cells
  //
  CI.set_minimal_cell_data();

  //
  // Set up an iterator to run through all the files.
  //
  list<string>::iterator ff = files.begin();
  //
  // Open first file, read header, and setup grid
  //
  ostringstream temp;
  temp << input_path << "/" << *ff;
  string first_file = temp.str();
  temp.str("");
  err = dataio.ReadHeader(first_file, SimPM);
  if (err) spdlog::error("{}: {}", "Didn't read header", err);


  //
  // write simulation xmin/xmax and radiation source position to a
  // text file.
  //
  redir.str("");
  redir << op_path << "/data_" << outfile << ".txt";
  ofstream outf;
  if (outf.is_open())
    spdlog::error("{}: {}", "Output text file is already open!", 1);
  outf.open(redir.str().c_str());
  if (!outf.is_open())
    spdlog::error(
        "{}: {}", "Failed to open text file for writing", redir.str());
  outf.setf(ios_base::scientific);
  outf.precision(6);
  outf << "Text data for simulation " << input_file << "\n\n";
  outf << "## GRID PROPERTIES ##\n";
  outf << "Xmin " << SimPM.Xmin[XX] << "  cm\n";
  outf << "Xmax " << SimPM.Xmax[XX] << "  cm\n";
  outf << "Ymin " << SimPM.Xmin[YY] << "  cm\n";
  outf << "Ymax " << SimPM.Xmax[YY] << "  cm\n";
  outf << "Zmin " << SimPM.Xmin[ZZ] << "  cm\n";
  outf << "Zmax " << SimPM.Xmax[ZZ] << "  cm\n";
  outf << "N_X " << SimPM.NG[XX] << "\n";
  outf << "N_Y " << SimPM.NG[YY] << "\n";
  outf << "N_Z " << SimPM.NG[ZZ] << "\n";
  outf << "#\n";
  if (SimPM.RS.Nsources > 0) {
    outf << "## RADIATION SOURCE ##\n";
    outf << "POS_X " << SimPM.RS.sources[0].pos[XX] << "  cm\n";
    outf << "POS_Y " << SimPM.RS.sources[0].pos[YY] << "  cm\n";
    outf << "POS_Z " << SimPM.RS.sources[0].pos[ZZ] << "  cm\n";
    outf << "Strength " << SimPM.RS.sources[0].strength << " erg/s";
    outf << "  blackbody source\n";
    outf << "T_star   " << SimPM.RS.sources[0].Tstar << " K. ";
    outf << "  Star's effective T.\n";
    outf << "#\n";
  }
  outf.close();

  //
  // get a setup_grid class, and use it to set up the grid!
  //
#ifdef PION_NESTED
  SimPM.levels.resize(SimPM.grid_nlevels);
  class setup_grid_NG_MPI *SimSetup = new setup_grid_NG_MPI();
  SimSetup->setup_NG_grid_levels(SimPM);
#else
  class setup_fixed_grid_pllel *SimSetup = new setup_fixed_grid_pllel();
#endif /* PION_NESTED */
  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  vector<class GridBaseClass *> G;
  G.resize(SimPM.grid_nlevels);
  SimSetup->setup_grid(G, SimPM);
  class GridBaseClass *grid = G[0];
  if (!grid) spdlog::error("{}: {}", "Grid setup failed", fmt::ptr(grid));
  SimPM.dx = grid->DX();
  spdlog::debug("\t\tg={}\tDX = {}", fmt::ptr(grid), grid->DX());

  //
  // setup microphysics class
  //
  err += SimSetup->setup_microphysics(SimPM);
  // err += setup_raytracing();
  if (err) spdlog::error("{}: {}", "Setup of microphysics", err);
  class microphysics_base *MP = SimSetup->get_mp_ptr();

  spdlog::info(
      "--------------- Finished Setting up Grid --------------\n-------------------------------------------------------");

  size_t ifile = 0;
  class DataIOFits_pllel writer(SimPM, &(SimPM.levels[0].sub_domain));

  spdlog::info(
      "-------------------------------------------------------\n--------------- Starting Loop over all input files ----\n-------------------------------------------------------");
  clk.start_timer("analyse_data");

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile = 0; ifile < nfiles; ifile++) {
    ifile += Nskip;
    spdlog::debug(
        "------ Starting Next Loop: ifile={}, time so far={}", ifile,
        clk.time_so_far("analyse_data ----"));
    // cout <<"-------------------------------------------------------\n";
    // cout <<"--------------- Reading Simulation data to grid -------\n";
    spdlog::info(" reading file {}", *ff);

    // *******************
    // * Read Input Data *
    // *******************
    //
    // Get filename with path (ff is a list of string filenames)
    //
    temp.str("");
    temp << input_path << "/" << *ff;
    string infile = temp.str();
    temp.str("");
    ff++;
    for (size_t skip = 0; skip < Nskip; skip++)
      ff++;

    // Read header to get timestep info.
    err = dataio.ReadHeader(infile, SimPM);
    if (err) spdlog::error("{}: {}", "Didn't read header", err);
    // cout.flush();

    // read data onto grid.
    err = dataio.ReadData(infile, G, SimPM);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "(main) Failed to read data", 0, err);

    spdlog::info(
        "--------------- Finished Reading Data  ----------------\n{}",
        "-------------------------------------------------------\n{}",
        "--------------- Starting Writing Data  ----------------");

    temp.str("");
    temp << op_path << "/" << outfile;
    writer.OutputData(temp.str(), G, SimPM, SimPM.timestep);
    spdlog::info("--------------- Finished Writing Data  ----------------");

  }  // Loop over all files.


  spdlog::info(
      "-------------------------------------------------------\n---- Finised with all Files: time={}--------\n-------------------------------------------------------\n--------------- Clearing up and Exiting ---------------",
      clk.stop_timer("analyse_data"));

  if (grid != 0) {
    spdlog::info("\t Deleting Grid Data...");
    delete grid;
    grid = 0;
  }
  if (MP) {
    delete MP;
    MP = 0;
  }

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------



// ##################################################################
// ##################################################################
