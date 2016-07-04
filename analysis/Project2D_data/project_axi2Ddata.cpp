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


#include "emission_absorption.h"
#include "xray_emission.h"



#define OP_TEXT 0
#define OP_FITS 1
#define OP_SILO 2
#define OP_VTK  3

//
// Variables to make projected values of
//
#define NIMG     10 ///< total number of images.
#define N_SCALAR 8  ///< number of scalar images.


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
  class Xray_emission XR;

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

    size_t iHp = 0; // H+ fraction is 1st tracer variable (if exists).
    size_t iWf = 0; // Wind fraction is 2nd tracer variable.
    //
    // Figure out which tracer variable is H+,wind, (only if there is
    // a microphysics class):
    //
    if (MP) {
      int tr = MP->Tr("H1+___");
      if ( tr==DONT_CALL_ME || tr<0  || !isfinite(tr) ) {
        cout <<"No H+ tracer variable, assuming all gas is ionized\n";
        iHp = 0;
      }
      else {
        iHp = tr;
      }
      //
      // Now for wind fraction (colour tracer):
      //
      if (SimPM.ntracer>0 && iHp==0) {
        iWf = SimPM.ftr;
      }
      else if (SimPM.ntracer>1 && iHp>0) {
        iWf = SimPM.ftr+1;
      }
      else {
        iWf = 0;
      }
    }

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
        //
        // If there is no H+ tracer, assume all gas is ionized.
        //
        if (iHp>0) {
          raw_data[DATA_TR0][iy] = cy->P[iHp];
        }
        else {
          raw_data[DATA_TR0][iy] = 1.0;
        }
        //
        // If there is no Wind colour tracer, assume there is no wind.
        //
        if (iWf>0) {
          raw_data[DATA_TR1][iy] = cy->P[iWf];
        }
        else {
          raw_data[DATA_TR1][iy] = 0.0;
        }
        //
        // Get Temperature from microphysics.
        //
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
      err = get_emission_absorption_data(raw_data, n_images, N_R, XR, ems_data, abs_data);


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
    //for (size_t ipix=0; ipix<num_pix; ipix++) {
    //  img_array[PROJ_X01][ipix] *= 5.09e48;
    //  img_array[PROJ_X05][ipix] *= 5.09e48;
    //  img_array[PROJ_X10][ipix] *= 5.09e48;
    //  img_array[PROJ_X50][ipix] *= 5.09e48;
    //}

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

  //free_xray_tables(&LT,&L1,&L2,&L3,&L4);
  //free_xray_tables();

  COMM->finalise();
  delete COMM; COMM=0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------



