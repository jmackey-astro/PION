///
/// file:    project2D.cpp
/// author:  Jonathan Mackey
/// date:    2019-11-26
///
/// Description: Read in Silo data outputs for 2D axisymmetric (z,R)
/// simulations and generate projected 3D maps of the simulation in
/// various quantities.  This code merges the two previous codes for
/// projection perpendicular to the axes and at a non-trivial angle
/// to the axes.
///

#include <cmath>
#include <fitsio.h>

#include <silo.h>
#include <sstream>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "sim_params.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"
#include "image_io.h"

#include "grid/setup_fixed_grid_MPI.h"
#include "sub_domain/sub_domain.h"

#include "microphysics/MPv3.h"
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"
#include "microphysics/MPv8.h"
#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"

#include "../xray/xray_emission.h"
#include "angle_projection.h"
#include "perp_projection.h"
#include "projection_constants.h"

#ifdef PROJ_OMP
#include <omp.h>
#endif

// ##################################################################
// ##################################################################


void setup_image_array(
    double ***img_array,          ///< pointer to be initialised.
    size_t NIMG,                  ///< number of images to make
    size_t,                       ///< total number of grid cells (unused)
    std::array<int, MAX_DIM> NG,  ///< number of Grid Cells in each direction
    int n_extra,                  ///< number of extra pixels in z-direction
    int npix[],     ///< OUTPUT: Number of pixels in each direction
    size_t &numpix  ///< OUTPUT: total number of pixels
)
{
  npix[0]    = NG[Zcyl] + 2 * n_extra;
  npix[1]    = NG[Rcyl];
  numpix     = npix[0] * npix[1];
  *img_array = mem.myalloc(*img_array, NIMG);
  for (size_t v = 0; v < NIMG; v++) {
    (*img_array)[v] =
        mem.myalloc((*img_array)[v], static_cast<long int>(numpix));
  }
}



// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  int err = 0;

  auto max_logfile_size = 1048576 * 5;
  auto max_logfiles     = 20;
#ifdef PARALLEL
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "project2D_pre_mpi", "project2D.log", max_logfile_size, max_logfiles));
#else
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "project2D", "project2D.log", max_logfile_size, max_logfiles));
#endif

#ifdef NDEBUG
  spdlog::set_level(spdlog::level::info);
  spdlog::flush_on(spdlog::level::err);
#else
  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_on(spdlog::level::trace);
#endif

  //
  // Also initialise the sub_domain class with myrank and nproc.
  //
  class SimParams SimPM;
  class Sub_domain *sub_domain = &(SimPM.levels[0].sub_domain);
  int myrank                   = sub_domain->get_myrank();
  // int nproc                    = SimPM.levels[0].sub_domain.get_nproc();

#ifdef PARALLEL
  spdlog::set_default_logger(spdlog::rotating_logger_mt(
      "project2D", "project2D_process_" + to_string(myrank) + ".log",
      max_logfile_size, max_logfiles));
#endif /* PARALLEL */

#ifdef PROJ_OMP
  omp_set_dynamic(0);
  omp_set_num_threads(12);
#endif

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc < 6) {
    spdlog::info("{}: {}", "Bad number of args", argc);
    spdlog::info(
        "Use as follows:\nexecutable-filename: <executable-filename> <input-path> <input-silo-file-base> <angle>");
    spdlog::info("\t\t <output-file> <op-file-type> <skip>");
    spdlog::info("******************************************");
    spdlog::info("input path:   path to input files.");
    spdlog::info("input file:   base filename of sequence of files.");
    spdlog::info(
        "angle:        Angle with respect to symmetry-axis for projection (degrees, float)");
    spdlog::info("output file:  filename for output file(s).");
    spdlog::info("op-file-type: integer/string [1,fits,FITS], [3,vtk,VTK]");
    spdlog::info(
        "skip:         will skip this number of input files each loop. (0 means it will calculate every file)");
    exit(1);
  }

  //
  // The code will get a list of files matching 'input_file' in the
  // directory 'input_path' and do the same operation on all files in
  // this list.
  //
  string input_path = argv[1];
  string input_file = argv[2];
  double angle      = atof(argv[3]);
  if (!isfinite(angle) || angle >= 180.0 || angle <= 0.0
      //|| pconst.equalD(angle,90.0)
  ) {
    spdlog::error("bad angle: {}, must be in range (0,180)", angle);
    exit(1);
  }
  int angle_int = static_cast<int>(angle * ONE_PLUS_EPS);
  angle *= M_PI / 180.0;

  //
  // outfile should contain the path as well (relative or absolute)
  //
  string outfile = argv[4];

  //
  // Save data as VTK or FITS images
  //
  string optype = argv[5];
  //
  // Set output file; if multiple files, append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to label output
  // file.
  //
  int op_filetype;
  if (optype == "1" || optype == "fits" || optype == "FITS") {
    op_filetype = OP_FITS;
    spdlog::info("\t\tsaving data to fits files");
  }
  else if (optype == "2" || optype == "silo" || optype == "SILO") {
    op_filetype = OP_SILO;
    spdlog::error("silo output files not implemented yet... fix me please!");
    exit(1);
  }
  else if (optype == "3" || optype == "vtk" || optype == "VTK") {
    op_filetype = OP_VTK;
    spdlog::info("\t\tsaving data to vtk files");
  }
  else {
    spdlog::error("{}: {}", "What sort of output is this?", optype);
    exit(1);
  }

  // write one image file for each snapshot
  int multi_opfiles = 1.0;

  // how sparsely to sample the data files.
  size_t skip = static_cast<size_t>(atoi(argv[6]));


  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  spdlog::info("-------------------------------------------------------");
  spdlog::info("--------------- Getting List of Files to read ---------");

  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio(SimPM, "DOUBLE", sub_domain);

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file, &files);
  if (err) spdlog::error("{}: {}", "failed to get list of files", err);
  for (list<string>::iterator s = files.begin(); s != files.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo") == string::npos) {
      // cout <<"removing file "<<*s<<" from list.\n";
      files.erase(s);
      s = files.begin();
    }
    else {
      // cout <<"files: "<<*s<<endl;
    }
  }
  int nfiles = static_cast<int>(files.size());
  if (nfiles < 1)
    spdlog::error("{}: {}", "Need at least one file, but got none", nfiles);

  spdlog::info("--------------- Got list of Files ---------------------");
  spdlog::info("-------------------------------------------------------");
  spdlog::info("--------------- Setting up Grid -----------------------");
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
  // First decompose the domain, so I know the dimensions of the local
  // grid to set up.  If nproc==1, then this sets the local domain to
  // be the full domain.
  //
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
  SimPM.levels[0].dx         = SimPM.Range[XX] / SimPM.NG[XX];
  SimPM.levels[0].simtime    = SimPM.simtime;
  SimPM.levels[0].dt         = 0.0;
  SimPM.levels[0].multiplier = 1;
  err                        = sub_domain->decomposeDomain(
      SimPM.ndim, SimPM.levels[0], SimPM.get_pbc_bools());
  if (err) {
    spdlog::error("{}: {}", "main: failed to decompose domain!", err);
    exit(1);
  }

  //
  // get a setup_grid class, and use it to set up the grid.
  //
  // CI.set_minimal_cell_data();
  class setup_fixed_grid *SimSetup = 0;
  SimSetup                         = new setup_fixed_grid_pllel();
  std::vector<class GridBaseClass *> cg;
  cg.resize(1);
  SimSetup->setup_grid(cg, SimPM);
  class GridBaseClass *grid = cg[0];
  if (!grid) {
    spdlog::error("{}: {}", "Grid setup failed", fmt::ptr(grid));
    exit(1);
  }
  spdlog::debug("\t\tg={}\tDX = {}", fmt::ptr(grid), grid->DX());
  err += SimSetup->setup_microphysics(SimPM);
  if (err) {
    spdlog::error("{}: {}", "Setup of microphysics and raytracing", err);
    exit(1);
  }
  class microphysics_base *MP = SimSetup->get_mp_ptr();

  //
  // This code needs 2d data to project...
  //
  if (SimPM.ndim != 2 || SimPM.coord_sys != COORD_CYL) {
    spdlog::error("projection needs 2D axisymmetric data: {}", SimPM.ndim);
    exit(1);
  }

  //
  // Setup X-ray emission tables.
  //
  class Xray_emission XR;

  spdlog::info("--------------- Finished Setting up Grid --------------");
  spdlog::info("-------------------------------------------------------");

  //
  // Output file: if multiple files, we will append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to
  // label output file.
  //
  class image_io imio;
  string filehandle("");
  string this_outfile("");

  // *****************************************************************
  // *****************************************************************
  // Images: we will write projected density, neutral density,
  // ionised density, emission measure, ...
  //
  double **img_array = 0;
  int n_images       = N_HD_SCALAR;
  size_t num_pix     = 0;
  int npix[2]        = {0, 0};
  double dx          = grid->DX();

  int n_extra = 0;
  if (angle_int != 90) {
    // size of image array must be larger than grid in Z-direction,
    // because we want to get glancing rays on both ends of the grid.
    // This means adding SimPM.Xmax[Rcyl]/tan(angle) to each end of the
    // grid.
    n_extra =
        static_cast<int>(fabs(SimPM.levels[0].Xmax[Rcyl] / tan(angle) / dx));
  }
  setup_image_array(
      &img_array, n_images, SimPM.Ncell, SimPM.NG, n_extra, npix, num_pix);

  //
  // array of image names for output files.
  //
  std::vector<string> im_name(n_images);
  for (int im = 0; im < n_images; im++) {
    switch (im) {
      case PROJ_D:
        im_name[im] = "Proj_SurfaceMass";
        break;
      case PROJ_NtD:
        im_name[im] = "Proj_NeutralDens";
        break;
      case PROJ_InD:
        im_name[im] = "Proj_IonizedDens";
        break;
      case PROJ_EM:
        im_name[im] = "Proj_EmissionMeasure";
        break;
      case PROJ_X00p1:
        im_name[im] = "Proj_XRAY_g00p1keV";
        break;
      case PROJ_X00p2:
        im_name[im] = "Proj_XRAY_g00p2keV";
        break;
      case PROJ_X00p3:
        im_name[im] = "Proj_XRAY_g00p3keV";
        break;
      case PROJ_X00p5:
        im_name[im] = "Proj_XRAY_g00p5keV";
        break;
      case PROJ_X01p0:
        im_name[im] = "Proj_XRAY_g01p0keV";
        break;
      case PROJ_X02p0:
        im_name[im] = "Proj_XRAY_g02p0keV";
        break;
      case PROJ_X05p0:
        im_name[im] = "Proj_XRAY_g05p0keV";
        break;
      case PROJ_X10p0:
        im_name[im] = "Proj_XRAY_g10p0keV";
        break;
      case PROJ_HA:
        im_name[im] = "Proj_Halpha";
        break;
      case PROJ_NII:
        im_name[im] = "Proj_NII_ll6584";
        break;
      case PROJ_BREMS6GHZ:
        im_name[im] = "Proj_BREMS6GHZ";
        break;
      default:
        spdlog::error("{}: {}", "Bad image count", im);
        exit(1);
        break;
    }
  }
  // *****************************************************************
  // *****************************************************************

  // *******************************************************************
  // loop over all files:
  // *******************************************************************

  spdlog::info("-------------------------------------------------------");
  spdlog::info("--------------- Starting Loop over all input files ----");
  spdlog::info("-------------------------------------------------------");
  unsigned int ifile = 0;
  clk.start_timer("total");
  double ttsf = 0.0;

  for (ifile = 0; ifile < static_cast<unsigned int>(nfiles);
       ifile += 1 + skip) {
    spdlog::debug(
        "--------- Starting Next Loop: ifile={}------\n========= reading file: {}",
        ifile, *ff);
    // cout <<"-------------------------------------------------------\n";
    // cout <<"--------------- Reading Simulation data to grid -------\n";

    //
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
    for (size_t q = 0; q <= skip; q++)
      ff++;

    //
    // Read header to get timestep info.
    //
    err = dataio.ReadHeader(infile, SimPM);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "(main) Didn't read header", 0, err);
      exit(1);
    }
    SimPM.grid_nlevels = 1;
    spdlog::debug("!! time = {}, step={}", SimPM.simtime, SimPM.timestep);

    //
    // Read data (this reader can read serial or parallel data.
    //
    spdlog::debug("reading from input file: {}", infile);
    err = dataio.ReadData(infile, cg, SimPM);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "(main) Failed to read data", 0, err);
      exit(1);
    }

    // cout <<"--------------- Finished Reading Data  ----------------\n";
    // cout <<"-------------------------------------------------------\n";
    // cout <<"--------------- Starting Data Analysis ----------------\n";
    clk.start_timer("analysis");
    double tsf = 0.0;
    //
    // ********************
    // * Analyse the Data *
    // ********************

    //
    // zero the image arrays.
    //
    for (size_t im = 0; im < static_cast<size_t>(n_images); im++)
      for (size_t ip = 0; ip < num_pix; ip++)
        img_array[im][ip] = 0.0;

    spdlog::debug("Xmin : {}", SimPM.Xmin);
    spdlog::debug("Xmax : {}", SimPM.Xmax);
    spdlog::debug("Angle w.r.t. symmetry axis={} degrees", angle_int);

    if (angle_int == 90) {
      err += generate_perpendicular_image(
          SimPM, MP, grid, XR, npix, num_pix, n_images, img_array);
    }
    else {
      err += generate_angle_image(
          SimPM, MP, grid, XR, angle, npix, num_pix, n_extra, n_images,
          img_array);
    }

    tsf = clk.stop_timer("analysis");
    spdlog::debug("\tFinished loop {}: loop time = {}", ifile, tsf);
    ttsf = clk.time_so_far("total");
    spdlog::debug(",\t total runtime={}\n", ttsf);
    spdlog::info("--------------- Finished Analysing this step ----------");
    spdlog::info("-------------------------------------------------------");
    spdlog::info("--------------- Writing image and getting next Im-file");

    //
    // **********************
    // * Write Data to file *
    // **********************
    //
    this_outfile = imio.get_output_filename(
        outfile, multi_opfiles, op_filetype, SimPM.timestep);
    err = imio.open_image_file(this_outfile, op_filetype, &filehandle);
    if (err) {
      spdlog::error("{}: {}", "failed to open output file", err);
      exit(1);
    }
    //
    // Write N images, here it is one for each variable.
    //
    double Xmin[2], pos[2];
    pixel_centre(SimPM.Xmin.data(), dx, n_extra, 0, 0, pos);
    Xmin[Zcyl]      = pos[Zcyl] * sin(angle);
    Xmin[Rcyl]      = grid->SIM_Xmin(Rcyl);
    double im_dx[2] = {grid->DX() * sin(angle), grid->DX()};
    // cout <<"saving image 0 of "<<n_images<<"\n";
    for (int outputs = 0; outputs < n_images; outputs++) {
      // cout <<"saving image "<<outputs<<" of "<<n_images<<"\n";
      // cout <<"@@ time = "<<SimPM.simtime<<", step="<<SimPM.timestep<<"\n";
      err = imio.write_image_to_file(
          filehandle, op_filetype, img_array[outputs], num_pix, 2, npix,
          im_name[outputs], Xmin, im_dx, SimPM.simtime, SimPM.timestep);
      if (err) {
        spdlog::error("{}: {}", "Failed to write image to file", err);
        exit(1);
      }
    }
    err = imio.close_image_file(filehandle);
    if (err) {
      spdlog::error("{}: {}", "failed to close output file", err);
      exit(1);
    }


  }  // Loop over all files.


  spdlog::info("-------------------------------------------------------");
  spdlog::info("-------------- Finished Analysing all Files -----------");
  spdlog::info("-------------------------------------------------------");
  spdlog::info("--------------- Clearing up and Exiting ---------------");

  //
  // free memory for images.
  //
  for (int v = 0; v < n_images; v++) {
    img_array[v] = mem.myfree(img_array[v]);
  }
  img_array = mem.myfree(img_array);

  if (grid != 0) {
    spdlog::info("\t Deleting Grid Data...");
    delete grid;
    grid = 0;
  }

  delete MP;
  MP = 0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------
