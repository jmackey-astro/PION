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
#include <iostream>
#include <silo.h>
#include <sstream>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "sim_params.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include "tools/timer.h"

#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"
#include "image_io.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_fixed_grid_MPI.h"

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
    double ***img_array,  ///< pointer to be initialised.
    size_t NIMG,          ///< number of images to make
    size_t Ncell,         ///< total number of grid cells
    int NG[],             ///< number of Grid Cells in each direction
    int n_extra,          ///< number of extra pixels in z-direction
    int npix[],           ///< OUTPUT: Number of pixels in each direction
    size_t &numpix        ///< OUTPUT: total number of pixels
)
{
  npix[0]    = NG[Zcyl] + 2 * n_extra;
  npix[1]    = NG[Rcyl];
  numpix     = npix[0] * npix[1];
  *img_array = mem.myalloc(*img_array, NIMG);
  for (int v = 0; v < NIMG; v++) {
    (*img_array)[v] =
        mem.myalloc((*img_array)[v], static_cast<long int>(numpix));
  }
  return;
}



// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  //
  // First initialise the comms class, since we need to define
  // PARALLEL to read a parallel file.
  //
  int err = COMM->init(&argc, &argv);
  //
  // Also initialise the MCMD class with myrank and nproc.
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);

  class MCMDcontrol *MCMD = &(SimPM.levels[0].MCMD);

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
    cout << "Use as follows:\n";
    cout
        << "executable-filename: <executable-filename> <input-path> <input-silo-file-base> <angle>\n";
    cout << "\t\t <output-file> <op-file-type>";
    cout << " <skip>\n";
    cout << "******************************************\n";
    cout << "input path:   path to input files.\n";
    cout << "input file:   base filename of sequence of files.\n";
    cout
        << "angle:        Angle with respect to symmetry-axis for projection (degrees, float)\n";
    cout << "output file:  filename for output file(s).\n";
    cout << "op-file-type: integer/string [1,fits,FITS], [3,vtk,VTK]\n";
    cout << "skip:         will skip this number of input files each loop. ";
    cout << "(0 means it will calculate every file)\n";
    rep.error("Bad number of args", argc);
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
    cerr << "bad angle: " << angle << ", must be in range (0,180)\n";
    exit(1);
  }
  int angle_int = static_cast<int>(angle * ONE_PLUS_EPS);
  angle *= M_PI / 180.0;

  //
  // outfile should contain the path as well (relative or absolute)
  //
  string outfile = argv[4];

  //
  // Redirect output to a text file if you want to:
  //
  //  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  ostringstream redir;
  redir.str("");
  redir << outfile << "_msg_" << myrank << "_";
  // rep.redirect(redir.str());


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
    cout << "\t\tsaving data to fits files.\n";
  }
  else if (optype == "2" || optype == "silo" || optype == "SILO") {
    op_filetype = OP_SILO;
    cout << "\t\tsaving data to silo files.\n";
    rep.error(
        "don't know how to output silo files yet... fix me please!", "sorry");
  }
  else if (optype == "3" || optype == "vtk" || optype == "VTK") {
    op_filetype = OP_VTK;
    cout << "\t\tsaving data to vtk files.\n";
  }
  else
    rep.error("What sort of output is this?", optype);

  // write one image file for each snapshot
  int multi_opfiles = 1.0;

  // how sparsely to sample the data files.
  size_t skip = static_cast<size_t>(atoi(argv[6]));


  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Getting List of Files to read ---------\n";

  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio(SimPM, "DOUBLE", MCMD);

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file, &files);
  if (err) rep.error("failed to get list of files", err);
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
  if (nfiles < 1) rep.error("Need at least one file, but got none", nfiles);

  cout << "--------------- Got list of Files ---------------------\n";
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Setting up Grid -----------------------\n";
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
  err = MCMD->decomposeDomain(SimPM.ndim, SimPM.levels[0]);
  if (err) rep.error("main: failed to decompose domain!", err);

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
  if (!grid) rep.error("Grid setup failed", grid);
  cout << "\t\tg=" << grid << "\tDX = " << grid->DX() << "\n";
  err += SimSetup->setup_microphysics(SimPM);
  if (err) rep.error("Setup of microphysics and raytracing", err);
  class microphysics_base *MP = SimSetup->get_mp_ptr();

  //
  // This code needs 2d data to project...
  //
  if (SimPM.ndim != 2 || SimPM.coord_sys != COORD_CYL) {
    rep.error("projection needs 2D axisymmetric data", SimPM.ndim);
  }

  //
  // Setup X-ray emission tables.
  //
  class Xray_emission XR;

  cout << "--------------- Finished Setting up Grid --------------\n";
  cout << "-------------------------------------------------------\n";

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
  string im_name[n_images];
  for (size_t im = 0; im < n_images; im++) {
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
        rep.error("Bad image count", im);
        break;
    }
  }
  // *****************************************************************
  // *****************************************************************

  // *******************************************************************
  // loop over all files:
  // *******************************************************************

  cout << "-------------------------------------------------------\n";
  cout << "--------------- Starting Loop over all input files ----\n";
  cout << "-------------------------------------------------------\n";
  unsigned int ifile = 0;
  clk.start_timer("total");
  double ttsf = 0.0;

  for (ifile = 0; ifile < static_cast<unsigned int>(nfiles);
       ifile += 1 + skip) {
    cout << "--------- Starting Next Loop: ifile=" << ifile << "------\n";
    cout << "========= reading file: " << *ff << "\n";
    cout.flush();
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
    rep.errorTest("(main) Didn't read header", 0, err);
    SimPM.grid_nlevels = 1;
    cout << "!! time = " << SimPM.simtime << ", step=" << SimPM.timestep
         << "\n";

    //
    // Read data (this reader can read serial or parallel data.
    //
    cout << "reading from input file: " << infile << "\n";
    err = dataio.ReadData(infile, cg, SimPM);
    rep.errorTest("(main) Failed to read data", 0, err);



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

    rep.printVec("Xmin", SimPM.Xmin, 2);
    rep.printVec("Xmax", SimPM.Xmax, 2);
    cout << "Angle w.r.t. symmetry axis=" << angle_int << " degrees\n";

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
    cout << "\tFinished loop " << ifile << ": loop time = " << tsf;
    ttsf = clk.time_so_far("total");
    cout << ",\t total runtime=" << ttsf << "\n";
    cout << "--------------- Finished Analysing this step ----------\n";
    cout << "-------------------------------------------------------\n";
    cout << "--------------- Writing image and getting next Im-file \n";

    //
    // **********************
    // * Write Data to file *
    // **********************
    //
    this_outfile = imio.get_output_filename(
        outfile, multi_opfiles, op_filetype, SimPM.timestep);
    err = imio.open_image_file(this_outfile, op_filetype, &filehandle);
    if (err) rep.error("failed to open output file", err);
    //
    // Write N images, here it is one for each variable.
    //
    double Xmin[2], pos[2];
    pixel_centre(SimPM.Xmin, dx, n_extra, 0, 0, pos);
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
      if (err) rep.error("Failed to write image to file", err);
    }
    err = imio.close_image_file(filehandle);
    if (err) rep.error("failed to close output file", err);


  }  // Loop over all files.


  cout << "-------------------------------------------------------\n";
  cout << "-------------- Finished Analysing all Files -----------\n";
  cout << "-------------------------------------------------------\n";
  cout << "--------------- Clearing up and Exiting ---------------\n";

  //
  // free memory for images.
  //
  for (int v = 0; v < n_images; v++) {
    img_array[v] = mem.myfree(img_array[v]);
  }
  img_array = mem.myfree(img_array);

  if (grid != 0) {
    cout << "\t Deleting Grid Data..." << endl;
    delete grid;
    grid = 0;
  }

  delete MP;
  MP = 0;
  delete COMM;
  COMM = 0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------
