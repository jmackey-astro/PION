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
#include "tools/reporting.h"
#include "tools/timer.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_fixed_grid_MPI.h"
#include "grid/uniform_grid.h"
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

#include <fitsio.h>
#include <iostream>
#include <silo.h>
#include <sstream>
using namespace std;


// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  //
  // First initialise the comms class
  //
  int err = COMM->init(&argc, &argv);
  //
  // Also initialise the MCMD class with myrank and nproc.
  // Get nproc from command-line (number of fits files for each
  // snapshot)
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);
  if (nproc > 1) rep.error("This is serial code", nproc);

  //
  // Get input files and an output file.
  //
  if (argc < 5) {
    cout << "Use as follows:\n";
    cout << "<executable-filename> <input-path>";
    cout << " <input-silo-file-base>";
    cout << " <output-path> <output-fits-file-base> [Nskip]";
    cout << " [<xmin> <xmax> [<ymin> <ymax> [<zmin> <zmax>]]]";
    cout << "\n";
    cout << "******************************************\n";
    cout << "input path:  path to input files.\n";
    cout
        << "input file:  base filename of sequence of files including _0000 if parallel.\n";
    cout << "output path: directory to write output files to.\n";
    cout << "output file: filename for output FITS file(s).\n";
    cout << "[Nskip]:     skip files each iteration.\n";
    cout << "xmin/xmax:   optional parameters to resize the domain.\n";
    rep.error("Bad number of args", argc);
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
  // rep.redirect(redir.str());

  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Getting List of Files to read ---------\n";

  //
  // set up dataio_utility class and fits-writer class.
  //
  class dataio_silo_utility dataio(SimPM, "DOUBLE", &(SimPM.levels[0].MCMD));

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file, &files);
  if (err) rep.error("failed to get list of files", err);
  //
  // remove non-silo files
  //
  for (list<string>::iterator s = files.begin(); s != files.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo") == string::npos) {
      cout << "removing file " << *s << " from list.\n";
      files.erase(s);
      s = files.begin();
    }
    else {
      cout << "files: " << *s << endl;
    }
  }
  size_t nfiles = files.size();
  if (nfiles < 1) rep.error("Need at least one file, but got none", nfiles);

  cout << "--------------- Got list of Files ---------------------\n";
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Setting up Grid -----------------------\n";

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
  if (err) rep.error("Didn't read header", err);


  //
  // write simulation xmin/xmax and radiation source position to a
  // text file.
  //
  redir.str("");
  redir << op_path << "/data_" << outfile << ".txt";
  ofstream outf;
  if (outf.is_open()) rep.error("Output text file is already open!", 1);
  outf.open(redir.str().c_str());
  if (!outf.is_open())
    rep.error("Failed to open text file for writing", redir.str());
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

  cout.flush();

  //
  // get a setup_grid class, and use it to set up the grid!
  //
#ifdef PION_NESTED
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
  if (!grid) rep.error("Grid setup failed", grid);
  SimPM.dx = grid->DX();
  cout << "\t\tg=" << grid << "\tDX = " << grid->DX() << endl;

  //
  // setup microphysics class
  //
  err += SimSetup->setup_microphysics(SimPM);
  // err += setup_raytracing();
  if (err) rep.error("Setup of microphysics", err);
  class microphysics_base *MP = SimSetup->get_mp_ptr();

  cout << "--------------- Finished Setting up Grid --------------\n";
  cout << "-------------------------------------------------------\n";

  size_t ifile = 0;
  class DataIOFits_pllel writer(SimPM, &(SimPM.levels[0].MCMD));

  cout << "-------------------------------------------------------\n";
  cout << "--------------- Starting Loop over all input files ----\n";
  cout << "-------------------------------------------------------\n";
  clk.start_timer("analyse_data");

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile = 0; ifile < nfiles; ifile++) {
    ifile += Nskip;
    cout << "------ Starting Next Loop: ifile=" << ifile << ", time so far=";
    cout << clk.time_so_far("analyse_data") << " ----\n";
    // cout <<"-------------------------------------------------------\n";
    // cout <<"--------------- Reading Simulation data to grid -------\n";
    cout << " reading file " << *ff << "\n";
    cout.flush();

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
    if (err) rep.error("Didn't read header", err);
    // cout.flush();

    // read data onto grid.
    err = dataio.ReadData(infile, G, SimPM);
    rep.errorTest("(main) Failed to read data", 0, err);

    cout << "--------------- Finished Reading Data  ----------------\n";
    cout << "-------------------------------------------------------\n";
    cout << "--------------- Starting Writing Data  ----------------\n";

    temp.str("");
    temp << op_path << "/" << outfile;
    writer.OutputData(temp.str(), G, SimPM, SimPM.timestep);
    cout << "--------------- Finished Writing Data  ----------------\n";

  }  // Loop over all files.


  cout << "-------------------------------------------------------\n";
  cout << "---- Finised with all Files: time=";
  cout << clk.stop_timer("analyse_data") << "--------\n";
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Clearing up and Exiting ---------------\n";

  if (grid != 0) {
    cout << "\t Deleting Grid Data..." << endl;
    delete grid;
    grid = 0;
  }
  if (MP) {
    delete MP;
    MP = 0;
  }

  delete COMM;
  COMM = 0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------



// ##################################################################
// ##################################################################
