///
/// file:    project_axi2Ddata.cc
/// author:  Jonathan Mackey
/// date:    2012-01-09
///
/// Description: Read in Silo data outputs for 2D axisymmetric (z,R)
/// simulations and generate projected 3D maps of the simulation in
/// various quantities.
///
/// Modifications:
///
/// - 2012-01-09 JM:  Adapted from pllel_analyse_data.cc
/// - 2012.01.10 JM: Debugged.  It now plots projected density.
/// - 2012.01.25 JM: Got rid of the GSL integration -- everything
///    could be integrated analytically for piecewise linear data, so
///    I just did that. Also implemented optically thin H-alpha
///    emission (i.e. no attenuation).
/// - 2012.11.23 JM: Updated the H-alpha emissivity to scale with
///    T^{-0.8} (Osterbrock, Spitzer textbooks).  I is now in units
///    of erg/cm3/s/steradian.
///    Set code to use MP->Temperature() to get gas temperature
///    (i.e. consistent with code).
///    Added in support for writing VTK files.
/// - 2013.10.08 JM: Moved from Betelgeuse project to RSG winds, and
///    changed the integration scheme.
/// - 2013.10.09 JM: Fixed bugs, corrected NII emission to account
///    for stellar wind being enriched in N and ISM having solar N.
/// - 2015.01.21 JM: Fixed bug in Halpha/NII lines where I was taking
///    steps that were too big in dTau, and absorption was sending
///    the intensity negative when passing through a dense shell!!!
///    This was corrected in the 1D projection code for W26 sims.
///    Also set it to skip the first file in the list (initial
///    conditions usually).
/// - 2015.01.21 JM: Added X-ray emission.
/// - 2015.01.21 JM: Switched off absorption and went back to the old
///    way of integrating.  When absorption coeff is very small, then
///    the source function becomes huge, and the method becomes
///    inaccurate.
/// - 2016.06.20 JM: Updated headers, new setup_grid structure.
///

#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
#include <cmath>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"
#include "tools/interpolate.h"

#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"
#include "dataIO/dataio_fits.h"
#include "image_io.h"

#include "MCMD_control.h"
#include "setup_fixed_grid_MPI.h"

#include "microphysics/microphysics_base.h"
#ifndef EXCLUDE_MPV1
#include "microphysics/microphysics.h"
#endif 
#ifndef EXCLUDE_HD_MODULE
#include "microphysics/microphysics_lowZ.h"
#endif 
#include "microphysics/mp_only_cooling.h"
#ifndef EXCLUDE_MPV2
#ifdef MP_V2_AIFA
#include "microphysics/mp_v2_aifa.h"
#endif
#endif 
#ifndef EXCLUDE_MPV3
#include "microphysics/mp_explicit_H.h"
#endif
#ifndef EXCLUDE_MPV4
#include "microphysics/mp_implicit_H.h"
#endif 

#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"






// ##################################################################
// ##################################################################
  size_t XNel=0;
  double *LT=0, *L1=0, *L2=0, *L3=0, *L4=0;
///
/// Function to set up the X-ray emissivity tables, from XSPEC
///
int setup_xray_tables(
      size_t *, ///< Number of elements read
      double **, ///< Log(T) array
      double **, ///< Log(L>0.1keV) array
      double **, ///< Log(L>0.5keV) array
      double **, ///< Log(L>1.0keV) array
      double **  ///< Log(L>5.0keV) array
      );

///
/// Function to free memory allocated for the X-ray emissivity tables
///
void free_xray_tables(
      double **, ///< Log(T) array
      double **, ///< Log(L>0.1keV) array
      double **, ///< Log(L>0.5keV) array
      double **, ///< Log(L>1.0keV) array
      double **  ///< Log(L>5.0keV) array
      );

///
/// Function to calculate the X-ray emissivity given an input gas
/// temperature.  Returns the emissivity in units of ergs.cm^3/sec,
/// and should be multiplied by n_e*n_H.
///
void get_xray_emissivity(
        const double, ///< Temperature (K)
        const size_t, ///< Number of elements in array.
        const double *, ///< Array of logT values.
        const double *, ///< Array of emissivities.
        const double *, ///< Array of emissivities.
        const double *, ///< Array of emissivities.
        const double *, ///< Array of emissivities.
        double *  ///< Results.
        );
// ##################################################################
// ##################################################################


#define OP_TEXT 0
#define OP_FITS 1
#define OP_SILO 2
#define OP_VTK  3

//
// data read from file is put in a 2D array with these values for
// the first index.
//
#define DATA_R   0 ///< radius
#define DATA_D   1 ///< density
#define DATA_P   2 ///< pressure
#define DATA_V   3 ///< velocity
#define DATA_TR0 4 ///< tracer 0
#define DATA_TR1 5 ///< tracer 1
#define DATA_T   6 ///< Temperature

//
// Variables to make projected values of
//
#define NIMG     10 ///< total number of images.
#define N_SCALAR 8  ///< number of scalar images.

#define PROJ_D   0  ///< projected density
#define PROJ_NtD 1  ///< projected neutral density
#define PROJ_InD 2  ///< projected ionised density
#define PROJ_EM  3  ///< projected emission measure.
#define PROJ_X01 4  ///< X-ray emission >0.1 keV
#define PROJ_X05 5  ///< X-ray emission >0.5 keV
#define PROJ_X10 6  ///< X-ray emission >1.0 keV
//#define PROJ_X20 7  ///< X-ray emission >2.0 keV
#define PROJ_X50 7  ///< X-ray emission >5.0 keV
#define PROJ_HA  8  ///< projected H-alpha emission.
#define PROJ_IML 9 ///< projected ionised metal-line.


//
// From the simulation data arrays, calculate various variables.
//
int get_emission_absorption_data(
      double const* const*, ///< raw data to get variable from
      const int,    ///< number of images to write
      const size_t, ///< Number of radial data elements
      double **,    ///< array for emission[img][rad] data.
      double **     ///< array for absorption[img][rad] data.
      );

//
// Project scalar quantities onto plane of the sky.
//
double calc_projection(
      const double *, ///< radius array
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius
      );

//
// Project quantities with emission and absorption onto sky.
//
double calc_projectionRT(
      const double *, ///< radius array.
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius.
      );


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
  class MCMDcontrol MCMD;
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  MCMD.set_myrank(myrank);
  MCMD.set_nproc(nproc);

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc<6) {
    cout << "Use as follows:\n";
    cout << "executable-filename: <executable-filename> <input-path> <input-silo-file-base>\n";
    cout << "\t\t <output-file> <op-file-type>";
    //cout << " <multi-opfiles>";
    cout << " <skip>\n";
    cout << "\t\t <ANY-EXTRA-STUFF> \n";
    cout <<"******************************************\n";
    cout <<"input path:   path to input files.\n";
    cout <<"input file:   base filename of sequence of filesn.\n";
    cout <<"output file:  filename for output file(s).\n";
    cout <<"op-file-type: integer/string [0,text,TEXT], [1,fits,FITS], [3,vtk,VTK]\n";
    //cout <<"muti-opfiles: integer. 1=one output file per step. (MUST CHOOSE 1))\n";
    cout <<"skip:         will skip this number of input files each loop.\n";
    rep.error("Bad number of args",argc);
  }

  //
  // The code will get a list of files matching 'input_file' in the
  // directory 'input_path' and do the same operation on all files in
  // this list.
  // 
  string input_path = argv[1];
  string input_file = argv[2];

  //
  // outfile should contain the path as well (relative or absolute)
  //
  string outfile    = argv[3];

  //
  // Redirect output to a text file if you want to:
  //
  //  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_"<<myrank<<"_";
  //rep.redirect(redir.str());


  //
  // What sort of output will depend on what sort of analysis we are
  // doing, here we are outputting one image per input file, so we should
  // choose multiple text or fits files.
  //
  string optype=argv[4];
  //
  // Set output file; if multiple files, append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to label output file.
  //
  int op_filetype;
  if      (optype=="0" || optype=="text" || optype=="TEXT") {
    op_filetype = OP_TEXT;
    cout <<"\t\toutputting data to text file.\n";
  }
  else if (optype=="1" || optype=="fits" || optype=="FITS") {
    op_filetype = OP_FITS;
    cout <<"\t\toutputting data to fits files.\n";
  }
  else if (optype=="2" || optype=="silo" || optype=="SILO") {
    op_filetype = OP_SILO;
    cout <<"\t\toutputting data to silo files.\n";
    rep.error("don't know how to output silo files yet... fix me please!","sorry");
  }
  else if (optype=="3" || optype=="vtk" || optype=="VTK") {
    op_filetype = OP_VTK;
    cout <<"\t\toutputting data to vtk files.\n";
  }
  else 
    rep.error("What sort of output is this?",optype);

  //
  // If we write a single output file for all steps or not.
  //
  int multi_opfiles = 1.0; //atoi(argv[5]);
  //switch (multi_opfiles) {
  //case 0:
  //  cerr <<"\t\tOutputting all timesteps in a single file.\n";
  //  rep.error("Can't do this!!!",multi_opfiles);
  //  break;
  //case 1:
  //  cout <<"\t\tOutputting timesteps in different files.\n";
  //  break;
  //default:
  //  rep.error("Bad multi-files value",multi_opfiles);
  //}

  //
  // how sparsely to sample the data files.
  //
  size_t skip = static_cast<size_t>(atoi(argv[5]));


  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Getting List of Files to read ---------\n";
  
  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio ("DOUBLE", &MCMD);

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file,  &files);
  if (err) rep.error("failed to get list of files",err);
  for (list<string>::iterator s=files.begin(); s!=files.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo")==string::npos) {
      cout <<"removing file "<<*s<<" from list.\n";
      files.erase(s);
    }
    else {
      cout <<"files: "<<*s<<endl;
    }
  }
  int nfiles = static_cast<int>(files.size());
  if (nfiles<1) rep.error("Need at least one file, but got none",nfiles);

  cout <<"--------------- Got list of Files ---------------------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Setting up Grid -----------------------\n";

  //
  // Set low-memory cells
  //
  CI.set_minimal_cell_data();

  //
  // Set up an iterator to run through all the files.
  //
  list<string>::iterator ff=files.begin();

  //
  // Open first file, read header, and setup grid
  //
  ostringstream temp; temp <<input_path<<"/"<<*ff;
  string first_file = temp.str();
  temp.str("");
  err = dataio.ReadHeader(first_file);
  if (err) rep.error("Didn't read header",err);

  //
  // First decompose the domain, so I know the dimensions of the local
  // grid to set up.  If nproc==1, then this sets the local domain to
  // be the full domain.
  //
  MCMD.decomposeDomain();
  if (err) rep.error("main: failed to decompose domain!",err);


  //
  // get a setup_grid class, and use it to set up the grid.
  //
  class setup_fixed_grid *SimSetup =0;
  SimSetup = new setup_fixed_grid_pllel();
  class GridBaseClass *grid = 0;
  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  SimSetup->setup_grid(&grid,&MCMD);
  if (!grid) rep.error("Grid setup failed",grid);
  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<"\n";
  double delr = grid->DX();

  //
  // This code needs 2d data to project...
  //
  if (SimPM.ndim!=2 || SimPM.coord_sys != COORD_CYL) {
    rep.error("projection needs 2D axisymmetric data",SimPM.ndim);
  }

  //
  // Now setup microphysics and raytracing classes
  //
  err += SimSetup->setup_microphysics();
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics and raytracing",err);

  //
  // Setup X-ray emission tables.
  //
  err += setup_xray_tables(&XNel, &LT, &L1, &L2, &L3, &L4);
  if (err) rep.error("Setup of xray-tables",err);
  //rep.printVec("LT",LT,XNel);
  //rep.printVec("L1",L1,XNel);
  //rep.printVec("L2",L2,XNel);
  //rep.printVec("L3",L3,XNel);
  //rep.printVec("L4",L4,XNel);

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Output file: if multiple files, we will append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to
  // label output file.
  //
  class image_io imio;
  string filehandle("");
  string this_outfile("");


  //*****************************************************************
  //*****************************************************************
  //
  // Set up data arrays for the radial data, and for the images.
  // Radial data has variables r, rho, p, v_r, tr0, tr1, temperature,
  // so 7 variables.
  //

  //
  // Images: we will write projected density, neutral density,
  // ionised density, emission measure, ...
  // 
  double **img_array=0;
  int n_images=NIMG;
  size_t num_pix = SimPM.Ncell; // Don't do reflected image yet.
  int N_R = SimPM.NG[Rcyl];
  int npix[NIMG] = {SimPM.NG[Zcyl],N_R};
  img_array = mem.myalloc(img_array,n_images);
  for (int v=0;v<n_images;v++)
    img_array[v] = mem.myalloc(img_array[v],
                               static_cast<long int>(num_pix));

  //
  // raw data
  //
#define NVAR 7
  double **raw_data=0;
  raw_data = mem.myalloc(raw_data,NVAR);
  for (int im=0; im<NVAR; im++) {
    raw_data[im] = mem.myalloc(raw_data[im], static_cast<long int>(N_R));
    for (int v=0; v<N_R; v++) raw_data[im][v] = 0.0;
  }

  //
  // columns of emission/absorption data for all output variables.
  //
  double **ems_data=0;
  double **abs_data=0;
  ems_data = mem.myalloc(ems_data,NIMG);
  for (int im=0; im<NIMG; im++) {
    ems_data[im] = mem.myalloc(ems_data[im],
                               static_cast<long int>(N_R));
    for (int v=0; v<N_R; v++) ems_data[im][v] = 0.0;
  }
  abs_data = mem.myalloc(abs_data,NIMG);
  for (int im=0; im<NIMG; im++) {
    abs_data[im] = mem.myalloc(abs_data[im],
                               static_cast<long int>(N_R));
    for (int v=0; v<N_R; v++) abs_data[im][v] = 0.0;
  }

  //
  // array of image names for output files.
  //
  string im_name[NIMG];
  for (size_t im=0; im<NIMG; im++) {
    switch (im) {
      case PROJ_D:   im_name[im] = "AA_SurfaceMass"; break;
      case PROJ_NtD: im_name[im] = "AA_ProjNeutralDens"; break;
      case PROJ_InD: im_name[im] = "AA_ProjIonisedDens"; break;
      case PROJ_EM:  im_name[im] = "AA_EmissionMeasure"; break;
      case PROJ_X01: im_name[im] = "AA_XRAY_g0p1keV"; break;
      case PROJ_X05: im_name[im] = "AA_XRAY_g0p5keV"; break;
      case PROJ_X10: im_name[im] = "AA_XRAY_g1p0keV"; break;
      //case PROJ_X20: im_name[im] = "XRAY_g2p0keV"; break;
      case PROJ_X50: im_name[im] = "AA_XRAY_g5p0keV"; break;
      case PROJ_HA:  im_name[im] = "AA_Halpha"; break;
      case PROJ_IML: im_name[im] = "AA_NII_ll6584"; break;
      default: rep.error("Bad image count",im); break;
    }
  }
  //*****************************************************************
  //*****************************************************************




  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";
  unsigned int ifile=0;
  clk.start_timer("total"); double ttsf=0.0;

  for (ifile=0; ifile<static_cast<unsigned int>(nfiles); ifile+=skip) {
    cout <<"--------- Starting Next Loop: ifile="<<ifile<<"------\n";
    cout <<"========= reading file: "<<*ff<<"\n";
    cout.flush();
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Reading Simulation data to grid -------\n";

    //
    // *******************
    // * Read Input Data *
    // *******************
    //
    // Get filename with path (ff is a list of string filenames)
    //
    temp.str("");
    temp <<input_path<<"/"<<*ff;
    string infile = temp.str();
    temp.str("");
    for (size_t q=0;q<skip;q++) ff++;

    //
    // skip first file in list.
    //
    if (ifile==0) {
      continue;
    }

    //
    // Read header to get timestep info.
    //
    err = dataio.ReadHeader(infile);
    if (err) rep.error("Didn't read header",err);

    //
    // Read data (this reader can read serial or parallel data.
    //
    err = dataio.parallel_read_any_data(infile, ///< file to read from
					grid    ///< pointer to data.
					);
    rep.errorTest("(main) Failed to read data",0,err);




    //cout <<"--------------- Finished Reading Data  ----------------\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Starting Data Analysis ----------------\n";
    clk.start_timer("analysis"); double tsf=0.0;
    //
    //********************
    //* Analyse the Data *
    //********************

    //
    // zero the image arrays.
    //
    for (size_t im=0; im< static_cast<size_t>(n_images); im++)
      for (size_t ip=0; ip<num_pix; ip++)
        img_array[im][ip] = 0.0;

    size_t iHp = SimPM.ftr;   // H+ fraction is 1st tracer.
    size_t iWf = SimPM.ftr+1; // Wind fraction is 2nd tracer.


    //
    // Loop over grid... Each z,R corresponds to a pixel.
    // Outer loop runs over z.
    //
    cell *cz = grid->FirstPt();
    cell *cy = 0;
    int iz=0;
    int iy=0;
    //double ori[2] = {0.0,0.0};
    double cpos[2];

    do {
      //cout <<"#+#+#+#+#+# New column, iz="<<iz<<"\n";
      //
      // First set the values of arrays in R for var,var_slope,R_cov.
      //
      cy = cz;
      iy = 0;

      //cout <<"\t\t assigning data from column to arrays.\n";
      do {
        CI.get_dpos(cy,cpos);

        raw_data[DATA_R][iy]   = cpos[Rcyl];
        raw_data[DATA_D][iy]   = cy->P[RO];
        raw_data[DATA_P][iy]   = cy->P[PG];
        raw_data[DATA_V][iy]   = cy->P[VY];
        raw_data[DATA_TR0][iy] = cy->P[iHp];
        raw_data[DATA_TR1][iy] = cy->P[iWf];
        raw_data[DATA_T][iy]   = MP->Temperature(cy->P,SimPM.gamma);

        // increment iy
        iy++;
      } while ((cy=grid->NextPt(cy,RPcyl))!=0 && cy->isgd);

      if (iy!=N_R)
        rep.error("Bad logic for radial grid size",iy-N_R);

      //
      // Now for this radial column of data, we set the value for
      // each output variable at the cell centre, including both an
      // emission and an absorption coefficient for some variables.
      //
      //cout <<"\t\t getting emission/absorption data.\n";
      err = get_emission_absorption_data(raw_data, n_images, N_R, ems_data, abs_data);


      //
      // Now go back to start of column in R, and for each pixel,
      // calculate the integral along a line of sight with impact
      // parameter b=R_i, for each variable.
      //
      double b = 0;

      for (int ivar=0; ivar<n_images; ivar++) {
        //cout <<"\t\t Calculating image im="<<ivar<<", name="<<im_name[ivar]<<"\n";
        for (int ipix=0; ipix<N_R; ipix++) {
          b = raw_data[DATA_R][ipix];
          if (ivar<N_SCALAR) {
            img_array[ivar][npix[Zcyl]*ipix +iz] = calc_projection(
                  raw_data[DATA_R],
                  ems_data[ivar], abs_data[ivar],
                  N_R, b, delr);
          }
          else {
            img_array[ivar][npix[Zcyl]*ipix +iz] = calc_projectionRT(
                  raw_data[DATA_R],
                  ems_data[ivar], abs_data[ivar],
                  N_R, b, delr);
          }
        }
      }

      //
      // increment iz and move to next radial column.
      //
      iz++;
    } while ((cz=grid->NextPt(cz,ZPcyl))!=0 && cz->isgd);

    //
    // Multiply X-ray intensity by 4*PI*(206165)^2*(3.086e18)^2 to
    // get units of luminosity per sq. parsec (erg/s/pc^2), so that
    // we can compare to Townsley et al. papers.
    // Otherwise it is erg/cm^2/s/sq.arcsec.
    //
    for (size_t ipix=0; ipix<num_pix; ipix++) {
      img_array[PROJ_X01][ipix] *= 5.09e48;
      img_array[PROJ_X05][ipix] *= 5.09e48;
      img_array[PROJ_X10][ipix] *= 5.09e48;
      img_array[PROJ_X50][ipix] *= 5.09e48;
    }

    tsf= clk.stop_timer("analysis");
    cout <<"\tFinished loop "<<ifile<<": loop time = "<<tsf;
    ttsf=clk.time_so_far("total");
    cout <<",\t total runtime="<<ttsf<<"\n";
    //cout <<"--------------- Finished Analysing this step ----------\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Writing image and getting next Im-file \n";

    //
    // **********************
    // * Write Data to file *
    // **********************
    //
    this_outfile = imio.get_output_filename(outfile, multi_opfiles, op_filetype, SimPM.timestep);
    err = imio.open_image_file(this_outfile, op_filetype, &filehandle);
    if (err) rep.error("failed to open output file",err);
    //
    // Write N images, here it is one for each variable.
    //
    //string im_name;
    //double *im=0;
    double Xmin[MAX_DIM];
    for (int v=0;v<SimPM.ndim;v++)
      Xmin[v] = grid->SIM_Xmin(static_cast<axes>(v));
    double im_dx[2] = {grid->DX(), grid->DX()};    
    for (int outputs=0;outputs<n_images;outputs++) {
      err = imio.write_image_to_file(filehandle, op_filetype,
                                    img_array[outputs], num_pix,
                                    2, npix, im_name[outputs],
                                    Xmin, im_dx, 
                                    SimPM.simtime, SimPM.timestep
                                    );
      if (err) rep.error("Failed to write image to file",err);
    }
    err = imio.close_image_file(filehandle);
    if (err) rep.error("failed to close output file",err);


  } // Loop over all files.    


  cout <<"-------------------------------------------------------\n";
  cout <<"-------------- Finished Analysing all Files -----------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";

  //
  // free memory for images.
  //
  for (int v=0;v<n_images;v++) {
    img_array[v] = mem.myfree(img_array[v]);
    ems_data[v] = mem.myfree(ems_data[v]);
    abs_data[v] = mem.myfree(abs_data[v]);
  }
  for (int v=0;v<NVAR;v++) {
    raw_data[v] = mem.myfree(raw_data[v]);
  }
  img_array = mem.myfree(img_array);
  ems_data  = mem.myfree(ems_data);
  abs_data  = mem.myfree(abs_data);
  raw_data  = mem.myfree(raw_data);


  if(grid!=0) {
    cout << "\t Deleting Grid Data..." << endl;
    delete grid; grid=0;
  }

  delete MP; MP=0;

  free_xray_tables(&LT,&L1,&L2,&L3,&L4);

  COMM->finalise();
  delete COMM; COMM=0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------




// ##################################################################
// ##################################################################



int get_emission_absorption_data(
      double const* const* data, ///< raw data to get variable from
      const int n_img,    ///< number of images to write
      const size_t Nr, ///< Number of radial data elements
      double **ems,    ///< array for emission[img][rad] data.
      double **abs     ///< array for absorption[img][rad] data.
      )
{
  double xr[5]; xr[0]=0.0,xr[1]=0.0,xr[2]=0.0,xr[3]=0.0,xr[4]=0.0;
  //
  // Loop over all image variables.
  //
  for (int im=0; im<n_img; im++) {

    //
    // Put data for variable ivar into ems[im][] array, and variable's
    // slope with radius in vsl[] array.
    // Use 1st order forward differencing for all slopes, except for
    // the last point which uses backwards differencing.
    //
    switch (im) {
      //
      // Projected density is easy.
      //
      case PROJ_D:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i];
        if (i==Nr-1) {
          abs[im][i] = abs[im][i-1];
        }
        else {
          abs[im][i] = (data[DATA_D][i+1]-data[DATA_D][i])/
                    (data[DATA_R][i+1]-data[DATA_R][i]);
        }
      }
      break;

      //
      // Projected Neutral density: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      //
      case PROJ_NtD:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i]*(1.0-data[DATA_TR0][i]);
      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/(data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // Projected ionised density: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      //
      case PROJ_InD:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i]*data[DATA_TR0][i];
      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/(data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // Projected emission measure: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      // Here we hardcode conversion units.   We assume n_e= 1.1n_p, 
      // appropriate if He ionization follows H.
      // We further assume X_H=0.715 (H mass fraction), similar to the
      // Asplund et al. (2009) value of 0.7154.
      // 
      case PROJ_EM:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i]*7.16e28;

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>0.1 keV
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X01:
      for (size_t i=0; i<Nr; i++) {
        get_xray_emissivity(data[DATA_T][i],XNel,LT,L1,L2,L3,L4,xr);
        ems[im][i] = 3.76e35*xr[0]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>0.5 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X05:
      for (size_t i=0; i<Nr; i++) {
        get_xray_emissivity(data[DATA_T][i],XNel,LT,L1,L2,L3,L4,xr);
        ems[im][i] = 3.76e35*xr[1]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>1.0 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X10:
      for (size_t i=0; i<Nr; i++) {
        get_xray_emissivity(data[DATA_T][i],XNel,LT,L1,L2,L3,L4,xr);
        ems[im][i] = 3.76e35*xr[2]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>5.0 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X50:
      for (size_t i=0; i<Nr; i++) {
        get_xray_emissivity(data[DATA_T][i],XNel,LT,L1,L2,L3,L4,xr);
        ems[im][i] = 3.76e35*xr[3]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      

      //
      // H-alpha emission:  here ems[im][] is the emissivity, and
      // abs[im][] the absorption coefficient.
      // Here we hardcode conversion units.  We assume n_e= 1.1n_p, 
      // appropriate if He ionization follows H.
      // We further assume X_H=0.715 (H mass fraction), similar to the
      // Asplund et al. (2009) value of 0.7154.
      //
      // Emissivity from Osterbrock j(Ha)=2.63e-33*n_e*n_p/T^0.9 in 
      // units of erg/cm3/s/sq.arcsec (adapted from table), and use
      //  n_e*n_p = rho^2 y^2 1.1(X_H/m_p)^2
      //
      // Absorption from Henney et al. 2009, where they assume 
      //  alpha = 5.0e-22 nH per cm (absorption by dust).
      //
      case PROJ_HA:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = 5.28e14 *data[DATA_D][i]*data[DATA_TR0][i]*
                          data[DATA_D][i]*data[DATA_TR0][i]*
                          pow(data[DATA_T][i], -0.9);
                          //exp(-0.9*log(data[DATA_T][i]));
        abs[im][i] = 0.0; //213.7 *data[DATA_D][i];
      }
      break;
      
      //
      // Ionised Metal-lines:  here ems[im][] is the emissivity, and
      // abs[im][] the absorption coefficient.
      // Here we hardcode conversion units.  We assume n_e= n_p, 
      // appropriate if He remains neutral.
      // We further assume H_H=0.715 (H mass fraction), similar to the
      // Asplund et al. (2009) value of 0.7154.
      //
      // Emissivity from Dopita (1973)
      //  j([NII] ll 6584) =
      //   6.82e-18 n_e*n_p*f(N)*exp(-chi/kT)/(4*pi*sqrt(T))
      // in units of erg/cm3/s/sr, and use
      //  n_e*n_p = rho^2 y^2 1.1(X_H/m_p)^2
      // and we convert to erg/cm2/s/sq.arcsec.
      // Also we have an exponential cutoff in emissivity at 10^5K.
      //
      // Furthermore, we assume N is strongly enhanced in the RSG
      // wind to f(N)=2.0e-4, or A(N)=8.3 (Brott+,2011), compared to
      // the ISM abundance of A(N)=7.85 or f(N)=7.08e-5, and we use
      // the second tracer to discriminate wind from ISM material.
      //
      // Absorption from Henney et al. 2009, where they assume 
      //  alpha = 5.0e-22 nH per cm (absorption by dust).
      //
      case PROJ_IML:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = 2.56e18 *data[DATA_D][i]*data[DATA_TR0][i]*
                          data[DATA_D][i]*data[DATA_TR0][i]*
                          exp(-2.1855e4/data[DATA_T][i])/
                          sqrt(data[DATA_T][i])
                          *exp(-(data[DATA_T][i]*data[DATA_T][i])/1.0e10);
        //ems[im][i] *= (1.0-data[DATA_TR1][i])*7.08e-5
        //              +data[DATA_TR1][i]*2.0e-4;
        ems[im][i] *= 7.08e-5;
        abs[im][i] = 0.0; //213.4 *data[DATA_D][i];
      }
      break;

      default:
      cerr <<"get_var_and_slope(): Don't know what to do for var ";
      cerr <<im<<"\n";
      return im;
      break;
    }
  }
  return 0;
}




// ##################################################################
// ##################################################################



double calc_projection(
      const double *r, ///< radius array
      const double *v, ///< array of values at each radius
      const double *s, ///< array of slopes at each radius
      const size_t Nr, ///< Size of arrays.
      double b,        ///< impact parameter of ray.
      const double dr  ///< spacing of points in radius
      )
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr-1]+0.5*dr;
  if (b > grid_max) {
    cout <<"calc_projection: Bad B value, b="<<b<<"\n";
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    //cout <<"b="<<b<<" is <r[0], so resetting b=r[0]+eps.\n";
    b = r[0]*1.00000001;
  }

  //
  // start at b, integrate outwards to rmax, and then multiply by 2.
  //
  size_t ir = 0;
  while ( (r[ir]+dr) < b) ir++;

  double result = 0.0;
  double Rmin=0.0, Rmax=0.0, Rmin2=0.0, Rmax2=0.0;
  double maxd=0.0, mind=0.0;
  double slope=0.0, offset=0.0, xdx=0.0, x2dx=0.0;

  do {
    //
    // Min/Max for this line segment.
    //
    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    Rmax2 = Rmax*Rmax;
    Rmin2 = Rmin*Rmin;
    maxd  = sqrt(Rmax2 - b*b);
    mind  = sqrt(Rmin2 - b*b);

    //
    // var(r) = slope*r+offset, where a=|dvar/dr|_i, b=var_i-|dvar/dr|_i*r_i
    //
    slope = s[ir];
    offset = v[ir] -s[ir]*Rmin;
    //cout <<"\ta="<<a<<", b="<<b;

    //
    // xdx is the integral xdx/sqrt(x^2-y^2).
    // x2dx is the integral x^2dx/sqrt(x^2-y^2).
    //
    xdx = maxd - mind;
    x2dx = 0.5*(Rmax*maxd - Rmin*mind +b*b*log((Rmax+maxd)/(Rmin+mind)));
    //cout <<"\txdx="<<xdx<<", x2dx="<<x2dx;

    //
    // now this line segment's integral is a*x2dx + xdx*b
    //
    xdx = slope*x2dx + offset*xdx;
    //cout <<"\tint="<<xdx<<"\n";
    
    //
    // Add to result, and increment counter
    //
    result += xdx;
    ir++;
  } while (ir<Nr);

  //
  // double result for the outward ray (by symmetry).
  //
  result *= 2.0;

  return result;
}



// ##################################################################
// ##################################################################



double calc_projectionRT(
      const double *r, ///< radius array
      const double *ve, ///< array of emission values at each radius.
      const double *va, ///< array of absorption values at each radius.
      const size_t Nr, ///< Size of arrays.
      double b,        ///< impact parameter of ray.
      const double dr  ///< spacing of points in radius
      )
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr-1]+0.5*dr;
  if (b > grid_max) {
    cout <<"calc_projectionRT: Bad B value, b="<<b<<"\n";
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    //cout <<"b="<<b<<" is <r[0], so resetting b=r[0].\n";
    b = r[0];
  }

  //
  // start at Rmax, integrate inwards to b, and then back out.
  // For projections with emission and absorption we use 1st order
  // integration with piecewise constant data, so the integral is
  // easier.
  //
  long int ir = Nr-1;

  double result = 0.0;
  double Rmin=0.0, Rmax=0.0;
  double maxd=0.0, mind=0.0;
  double dIds=0.0;

  while (r[ir]*ONE_PLUS_EPS>b) {
    //
    // Min/Max for this line segment.
    //
    //cout <<"IN  b="<<b<<", r["<<ir<<"] = "<<r[ir]<<"\n";

    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    maxd  = sqrt(Rmax*Rmax - b*b);
    mind  = sqrt(Rmin*Rmin - b*b);
    // 
    // dI/ds = j - alpha*I (emission minus absorption)
    //
    dIds = ve[ir] -va[ir]*result;
    //
    // the integral xdx/sqrt(x^2-y^2) = maxd-mind.
    // Add this segment to result, and increment counter
    //
    result += dIds*(maxd - mind);

    ir--;
    if (ir<0) {
      //cout <<"ir="<<ir<<"\n";
      break;
    }
  }

  //
  // Continue for the outward ray. We went one cell too far in the
  // inward loop, so increment before we begin.
  //
  ir++;
  while (ir<static_cast<long int>(Nr)) {
    //cout <<"OUT b="<<b<<", r["<<ir<<"] = "<<r[ir];

    //
    // Min/Max for this line segment.
    //
    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    maxd  = sqrt(Rmax*Rmax - b*b);
    mind  = sqrt(Rmin*Rmin - b*b);
    //
    // dI/ds = j - alpha*I (emission minus absorption)
    //
    dIds = ve[ir] -va[ir]*result;
    //
    // the integral xdx/sqrt(x^2-y^2) = maxd-mind.
    // Add this segment to result, and increment counter
    //
    result += dIds*(maxd - mind);


    ir++;
  }

  return result;
}


// ##################################################################
// ##################################################################




// ##################################################################
// ##################################################################


///
/// Function to set up the X-ray emissivity tables, from XSPEC
///
int setup_xray_tables(
      size_t *N,
      double **logT,
      double **logLum1,
      double **logLum2,
      double **logLum3,
      double **logLum4
      )
{
  //
  // Read data from xray-table.txt
  //
  //
  // read input parameters:
  //
  ifstream ff;
  ff.open("xray-table.txt");
  if (!ff.is_open()) {
    cerr<<"Error opening file.\n";
    return 1;
  }

  //
  // Read the first line, which is the column headers:
  //
  string line;
  getline(ff, line);

  //
  // Now read in the data line-by-line.  Cols. are:
  // log(T/K) T(K) E(keV) L(E>0.1keV) L(E>0.5keV) L(E>1keV) L(E>5keV)
  // 
  // We need to store logT, L1, L2, L3, L4
  //
  vector<double> LT, L1, L2, L3, L4;
  double lt=0, l1=0, l2=0, l3=0, l4=0;
  double temp=0;
  size_t Nel=0;
  while (!ff.eof()) {
    getline(ff, line);
    // If it's not a line containing a parameter, continue and read another line.
    if (line.empty()==true || (line.find("#")!=string::npos) ) { 
      continue;
    }
    else {
      // We have found a line of data, so read it into variables.
      istringstream iss(line);
      iss >> lt >> temp >> temp >> l1 >> l2 >> l3 >> l4;
#ifdef TESTING
      cout << lt<<"  "<<l1<<"  "<<l2<<"  "<<l3<<"  "<<l4<<"\n";
#endif
      //
      // Add all of them to the vectors:
      //
      LT.push_back(lt);
      L1.push_back(l1);
      L2.push_back(l2);
      L3.push_back(l3);
      L4.push_back(l4);
      Nel++;
    }
  }
  cout <<"Read "<<Nel<<" lines of data from input xray tables.\n";
  ff.close();


  //
  // allocate arrays
  //
  (*logT) = mem.myalloc((*logT),Nel);
  *logLum1 = mem.myalloc(*logLum1,Nel);
  *logLum2 = mem.myalloc(*logLum2,Nel);
  *logLum3 = mem.myalloc(*logLum3,Nel);
  *logLum4 = mem.myalloc(*logLum4,Nel);

  //cout <<"(*logT)="<<(*logT)<<" (logT)="<<(logT)<<"\n";
  
  for (size_t v=0; v<Nel; v++) {
    //cout <<"v="<<v<<" "<<LT[v];
    (*logT)[v] = LT[v];
    //cout <<" "<<L1[v];
    (*logLum1)[v] = log10(L1[v]);
    //cout <<" "<<L2[v];
    (*logLum2)[v] = log10(L2[v]);
    //cout <<" "<<L3[v];
    (*logLum3)[v] = log10(L3[v]);
    //cout <<" "<<L4[v];
    (*logLum4)[v] = log10(L4[v]);
    //cout <<"\n";
  }

  *N = Nel;

  return 0;
}

//  double t[17]= {5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2};
//  double l1[17]={5.8367e-27, 6.8811e-26, 2.5779e-24, 1.5294e-23, 1.5956e-23, 1.7678e-23, 1.7051e-23, 1.8891e-23, 2.0251e-23, 2.6735e-23, 2.8442e-23, 1.8402e-23, 1.5671e-23, 1.7602e-23, 2.1008e-23, 2.5335e-23, 3.0165e-23}; 
//  double l2[17]={, , , , , , , , , , , , , , , , , , , , , , 
//  double l3[17]={, , , , , , , , , , , , , , , , , , , , , , 
//  double l4[17]={, , , , , , , , , , , , , , , , , , , , , , 

void free_xray_tables(
      double **logT,
      double **logL1,
      double **logL2,
      double **logL3,
      double **logL4
      )
{
  *logT = mem.myfree(*logT);
  *logL1 = mem.myfree(*logL1);
  *logL2 = mem.myfree(*logL2);
  *logL3 = mem.myfree(*logL3);
  *logL4 = mem.myfree(*logL4);
  return;
}

// ##################################################################
// ##################################################################



/// Function to calculate the X-ray emissivity given an input gas
/// temperature.  Returns the emissivity in units of ergs.cm^3/sec,
/// and should be multiplied by n_e*n_H.
///
void get_xray_emissivity(
        const double T, ///< Temperature (K)
        const size_t Nel, ///< Number of elements in array.
        const double *LT, ///< Array of logT values.
        const double *L1, ///< Array of emissivities.
        const double *L2, ///< Array of emissivities.
        const double *L3, ///< Array of emissivities.
        const double *L4, ///< Array of emissivities.
        double *res  ///< Results.
        )
{
  double ln10 = 2.302585093;
  double lrate=0.0;
  double lt = log10(T);

  //
  // extrapolate linearly
  //
  if (lt<LT[0]) {

    //cout << "extrapolate, T="<<T<<" : ";
    //res[0] = L1[0]+(L1[1]-L1[0])*(lt-LT[0])/(LT[1]-LT[0]);
    //res[1] = L2[0]+(L2[1]-L2[0])*(lt-LT[0])/(LT[1]-LT[0]);
    //res[2] = L3[0]+(L3[1]-L3[0])*(lt-LT[0])/(LT[1]-LT[0]);
    //res[3] = L4[0]+(L4[1]-L4[0])*(lt-LT[0])/(LT[1]-LT[0]);
    //for (size_t v=0; v<4; v++) res[v] = exp(ln10*res[v]);
    //rep.printVec("Res",res,4);
    for (size_t v=0; v<4; v++) res[v] = 0.0;
  }
  else if (lt>LT[Nel-1]) {
    cout << "extrapolate, T="<<T<<" : ";
    res[0] = L1[Nel-1]+(L1[Nel-1]-L1[Nel-2])*(lt-LT[Nel-1])/(LT[Nel-1]-LT[Nel-2]);
    res[1] = L2[Nel-1]+(L2[Nel-1]-L2[Nel-2])*(lt-LT[Nel-1])/(LT[Nel-1]-LT[Nel-2]);
    res[2] = L3[Nel-1]+(L3[Nel-1]-L3[Nel-2])*(lt-LT[Nel-1])/(LT[Nel-1]-LT[Nel-2]);
    res[3] = L4[Nel-1]+(L4[Nel-1]-L4[Nel-2])*(lt-LT[Nel-1])/(LT[Nel-1]-LT[Nel-2]);
    for (size_t v=0; v<4; v++) res[v] = exp(ln10*res[v]);
    rep.printVec("Res",res,4);
  }
  else {
    //cout <<"interpolate, T="<<T<<" : ";
    interpolate.root_find_linear(LT, L1, Nel, lt, &lrate);
    res[0] = exp(ln10*lrate);
    interpolate.root_find_linear(LT, L2, Nel, lt, &lrate);
    res[1] = exp(ln10*lrate);
    interpolate.root_find_linear(LT, L3, Nel, lt, &lrate);
    res[2] = exp(ln10*lrate);
    interpolate.root_find_linear(LT, L4, Nel, lt, &lrate);
    res[3] = exp(ln10*lrate);
    //rep.printVec("Res",res,4);
  }

  return;
}


// ##################################################################
// ##################################################################





