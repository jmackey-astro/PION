///
/// \file main_projection.cc
///
/// \author Jonathan Mackey
///
/// This file contains main() for calculating projections through a 3D sim.
/// It has the following program flow:\n
///  - get a directory listing of files matching the infile argument
///  - Assuming they are silo files, read the first one's header and set up grid.
///  - Set up image; add image positions to grid cells;
///  - Add integration points to image pixels; associate 4 cells with each point.
///  - Loop over input sequence of files:
///  -  - Read in data for the file.
///  -  - For each pixel, calculate the line-of-sight integral and add data to image array.
///  -  - Write image to output file.
///  - Clean up and exit.
///
/// File to analyse a sequence of files from a photo-evaporating random
/// clumps simulation.  First we get the directory listing, then for
/// each file we load it onto the grid, run some analysis on it, output
/// results to a file, and continue to the next file.
///
/// The main thing this code does it creates a 2D image of either
/// projected, mean, or integrated quantities along the line of sight in
/// each pixel. The lines of sight can be aligned with the computational
/// grid, in which case it is very easy; or else at an angle to the grid,
/// in which case we need to do a lot more work.  
///
/// Modifications:\n
/// - 2009-06-25 Created file, split from sim_projection.cc
/// - 2010-03-22 JM: Moved to MHD_ET_2010, modified to also get projected field.
///    Simulation domain is hard-coded to be a subset of the 384x256x256 sims.
/// - 2010-04-17 JM: minor mods, got it working on furfur (Turlough's machine).
/// - 2010-09-08 JM: Updated LOS velocity so that you can smooth with a user set width (e.g. 1km/s)
/// - 2010.12.13 JM: Added  NEW_STOKES_CALC ifdef to Makefile; the new
///    code in the ifdef does a different Stokes Q,U calculation and
///    replaces the projected Bx,By with values calculated from Q,U.
/// - 2012.12.05 JM: Added section to subtract off mean values from
///    projected total and neutral density.
/// - 2013.10.14 JM: Added microphysics support (for getting T, n_H,
///    n_e, etc.).  Added geometric grid.  Added [NII] 6584AA
///    emission calculation.  Added VTK output option.
/// - 2015.07.03 JM: Got rid of NEW_STOKES_CALC (because old code was
///    broken, so I just deleted the old code).
/// - 2015.07.03 JM: updated for pion_dev: uses MCMD, SimSetup,
///    constants.h
/// - 2015.07.13 JM: debugged and fixed a few things.


///
/// This calls the function RESET_DOMAIN() every timestep to make the
/// simulation domain smaller.  If not set, this does nothing.  If
/// set, you need to hardcode the new Xmin[],Xmax[] in the function.
///
//#define RESET_DOMAIN

///
/// If set, this subtracts the mean density from column density
/// images.
///
#define SUBTRACT_MEAN

#include <iostream>
#include <sstream>
#include <silo.h>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"

#include "dataIO/dataio.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"

#include "grid/uniform_grid.h"

#include "sim_projection.h"
#include "image_io.h"

#include "MCMD_control.h"
#include "setup_fixed_grid_MPI.h"


#include "raytracing/raytracer_SC.h"
// ----------- MICROPHYSICS --------------



#ifdef THREADS
#include "andys_threads/msvc_constants.h"
#if defined(_DEBUG) &&  defined(_MSC_VER) &&  defined(MSVC_DEBUG_NEW_TRACE_ON)
  #define CRTDBG_MAP_ALLOC
  #include <stdlib.h> 
  #include <crtdbg.h> 
  #define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif
#include "andys_threads/reefa_constants.h"
#include "andys_threads/logmessages.h"
#include "andys_threads/threadpool/threadpool.h"
//
// Global threading variables.
//
threadpool_t     tp; // main threadpool
int monsecs_gl=0;    // seconds since the start of the month

struct calc_pix_args {
  class image *IMG; ///< pointer to image class.
  struct pixel *px; ///< pointer to pixel
  int i;      ///< pixel id
  int what_to_integrate; ///< flag for what to integrate.
  double *im;       ///< array of pixel data.
  double *tot_mass; ///< general purpose counter for stuff.
  struct vel_prof_stuff *vps; ///< struct with info for velocity binning.
};

//
// void function for threading with Andy's threadpool library
//
void calculate_pixelW(void *arg)
{
  struct calc_pix_args *ta = reinterpret_cast<struct calc_pix_args *>(arg);
  ta->IMG->calculate_pixel(ta->px,
			   ta->vps,
			   ta->what_to_integrate,
			   ta->im,
			   ta->tot_mass
			   );
  delete ta; ta=0;
  return;
}
#endif //THREADS



// ----------- MICROPHYSICS --------------
///
/// Reset the radiation sources in the header to correspond to projected
/// quantities and not the sources used for the simulation.
///
void reset_radiation_sources(struct rad_sources *);

// ----------- MICROPHYSICS --------------


//
// ------------------------------------------------------------------------
// MODIFYING MIN/MAX SO THAT I ONLY READ IN A SUBDOMAIN OF THE FULL GRID
//
#ifdef RESET_DOMAIN
void reset_domain()
{
  rep.printVec("Old Xmin",SimPM.Xmin,SimPM.ndim);
  rep.printVec("Old Xmax",SimPM.Xmax,SimPM.ndim);
  for (int v=0;v<SimPM.ndim;v++) {
    SimPM.Xmin[v] = -66.16384e18;
    SimPM.NG[v] = static_cast<int>(ONE_PLUS_EPS*SimPM.NG[v]*(SimPM.Xmax[v]-SimPM.Xmin[v])/SimPM.Range[v]);
    SimPM.Range[v] = SimPM.Xmax[v] - SimPM.Xmin[v];
  }
  rep.printVec("New Xmin",SimPM.Xmin,SimPM.ndim);
  rep.printVec("New Xmax",SimPM.Xmax,SimPM.ndim);
  mpiPM.decomposeDomain();
  return;
}
#endif
//
// ------------------------------------------------------------------------
//

int main(int argc, char **argv)
{
  //
  // First initialise MPI, even though this is a single processor
  // piece of code.
  //
  int err = COMM->init(&argc, &argv);
  //
  // Also initialise the MCMD class with myrank and nproc.
  //
  class MCMDcontrol MCMD;
  int r=-1, np=-1;
  COMM->get_rank_nproc(&r,&np);
  MCMD.set_myrank(r);
  MCMD.set_nproc(np);

#ifdef THREADS
  tp_init(&tp,NUM_THREADS_MAIN,"Main Threadpool");
#endif //THREADS

  //cout <<"WARNING! SIMULATION DOMAIN HARD-CODED FOR MHD ET SIMULATIONS.\n";
  //cout <<"********************* USE WITH CAUTION! ********************\n\n";

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc<10) {
    cout << "Use as follows:\n";
    cout << "projection: <projection> <input-path> <input-silo-file-base>\n";
    cout << "\t\t <output-file> <op-file-type> <multi-opfiles> \n";
    cout << "\t\t <normalvec> <fixed_dir> <theta> <what2integrate> +[optional velocity args] \n";
    cout <<"******************************************\n";
    cout <<"input path:   path to input files.\n";
    cout <<"input file:   base filename of sequence of filesn.\n";
    cout <<"output file:  filename for output file(s).\n";
    cout <<"op-file-type: integer. 0=text, 1=fits, 3=vtk.\n";
    cout <<"muti-opfiles: integer. 0=only one output file. 1=one output file per step.\n";
    cout <<"normal-vec:   string direction for LOS viewing normal to calculate angle from: XN,XP,YN,YP,ZN,ZP.\n";
    cout <<"fixed-dir:    string direction for axis which we rotate view around (XN,XP mean same thing).\n";
    cout <<"theta:        angle (DEGREES, in [-89,89]) which LOS makes with normal-vec (staying perp. to fixed dir).\n"; 
    cout <<"what2integrate: integer: 0=density, 1=neutral num.density, 2=los velocity, 3=VX, 4=Halpha\n";
    cout <<"                         5=StokesQ,6=StokesU, 7=All-scalars, 8=|B|(LOS), 9=|B|(perp)";
    cout <<"OPTIONAL VELOCITY ARGS:\n";
    cout <<"Nbins:        integer number of bins in velocity profile.\n";
    cout <<"v_min:        minimum velocity to measure (in same units as output files!)\n";
    cout <<"v_max:        maximum velocity to measure (in same units as output files!)\n";
    cout <<"smooth:       integer =1 for constant broadening, =2 for doppler broadening by temperature.\n";
    cout <<"smoooth-val:  float for the velocity by which to smooth, if constant broadening (FWHM)\n";
    rep.error("Bad number of args",argc);
  }
  string input_path = argv[1];
  string input_file = argv[2];
  string outfile    = argv[3];

  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  rep.redirect(redir.str());

  //
  // start a timer, so I can see how long each step takes.
  //
  clk.start_timer("mainloop"); double mltsf=0.0;
  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* Starting code,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();


  int op_filetype = atoi(argv[4]);
  switch (op_filetype) {
  case 0:
    cout <<"\t\toutputting data to text file.\n";
    break;
  case 1:
    cout <<"\t\toutputting data to fits files.\n";
    break;
  case 3:
    cout <<"\t\toutputting data to VTK files.\n";
    break;
  default:
    rep.error("Bad outfile format",op_filetype);
  }

  int multi_opfiles = atoi(argv[5]);
  switch (multi_opfiles) {
  case 0:
    cout <<"\t\tOutputting all timesteps in a single file.\n";
    break;
  case 1:
    cout <<"\t\tOutputting timesteps in different files.\n";
    break;
  default:
    rep.error("Bad multi-files value",multi_opfiles);
  }

  //
  // Normal vector to trace rays along (and measure angle from (if non-zero)).
  //
  string n=argv[6]; enum direction normal=NO;
  if      (n=="XN") normal=XN;
  else if (n=="XP") normal=XP;
  else if (n=="YN") normal=YN;
  else if (n=="YP") normal=YP;
  else if (n=="ZN") normal=ZN;
  else if (n=="ZP") normal=ZP;
  else rep.error("Bad normal direction",n);
  //
  // Perpendicular direction, around which to rotate the view.
  //
  n=argv[7]; enum direction perpdir=NO;
  if      (n=="XN") perpdir=XN;
  else if (n=="XP") perpdir=XP;
  else if (n=="YN") perpdir=YN;
  else if (n=="YP") perpdir=YP;
  else if (n=="ZN") perpdir=ZN;
  else if (n=="ZP") perpdir=ZP;
  else rep.error("Bad perp direction",n);
  //
  // Angle to normal direction, and perp. to perp. direction.
  //
  int th=atoi(argv[8]);
  if (th<-89 || th>89 || isnan(th) || isinf(th))
    rep.error("Input angle is not on [-89,89]. please input a valid angle.",th);
  bool zero_angle;
  if (th==0) zero_angle=true;
  else      zero_angle=false;
  int angle = th;

  th = atoi(argv[9]);
  if (th<0 || th>9 || isnan(th) || isinf(th))
    rep.error("Input what-to-integrate is not in [0,1,2,3,4,7,8,9]",th);
  int what_to_integrate=th;
  // 0=density, 1=neutral density, 2=LOS velocity, 3=VX, 4=emission,
  // 7=all_scalars, [5,6]=stokesQU, [8,9]=|B|(LOS,PERP)

  //
  // Number of bins for each pixel value -- if getting velocity profiles we
  // will want this to be >1.  Should make this be dynamically allocateable.
  //
  int Nbins=1, smooth=0;
  double v_min=0.0, v_max=0.0, bin_size=0.0, broadening=0.0;
  struct vel_prof_stuff vps;
  vps.v_min = v_min;
  vps.v_max = v_max;
  vps.smooth = smooth;
  
  if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
    if (argc<13) rep.error("Need at least 13 args for velocity profiling",argc);
    Nbins = atoi(argv[10]);
    v_min = atof(argv[11]);
    v_max = atof(argv[12]);
    smooth= atoi(argv[13]);
    broadening=atof(argv[14]);
    bin_size = (v_max-v_min)/Nbins;
    cout <<"Velocity info: max="<<v_max<<", min="<<v_min<<", binsize="<<bin_size;
    cout <<", Nb="<<Nbins<<", smooth="<<smooth<<", broaden by "<<broadening<<" cm/s\n";
    vps.v_min = v_min;
    vps.v_max = v_max;
    vps.smooth = smooth;
    vps.broadening=broadening;
  }

  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Getting List of Files to read ---------\n";
  
  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio("DOUBLE",&MCMD);

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file,  &files);
  if (err) rep.error("failed to get list of files",err);
  //
  // Remove non-SILO files from list
  //
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
  //
  // Set up an iterator to run through all the files.
  //
  list<string>::iterator ff=files.begin();

  unsigned int nfiles = files.size();
  if (nfiles<1) rep.error("Need at least one file, but got none",nfiles);

  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* Files read,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();
  cout <<"--------------- Got list of Files ---------------------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Setting up Grid -----------------------\n";

  //
  // Set low-memory cells
  //
  CI.set_minimal_cell_data();

  //
  // Open first file, read header, and setup grid
  //
  ostringstream temp; temp <<input_path<<"/"<<*ff;
  string first_file = temp.str();
  temp.str("");
  err = dataio.ReadHeader(first_file);
  if (err) rep.error("Didn't read header",err);

  //
  // ------------------------------------------------------------------------
  // MODIFYING XMIN/XMAX SO THAT I ONLY READ IN A SUBDOMAIN OF THE FULL GRID
  // Full X-range, Yrange=[-0.75,0.75]pc, Zrange=Yrange.
  //
#ifdef RESET_DOMAIN
  reset_domain();
#endif
  //
  // ------------------------------------------------------------------------
  //

  //
  // First decompose the domain, so I know the dimensions of the local
  // grid to set up.  If nproc==1, then this sets the local domain to
  // be the full domain.
  //
  err  = MCMD.decomposeDomain();
  if (err) rep.error("main: failed to decompose domain!",err);
  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables (here just set it to zero).
  //
  // *****************************************************
  // Now delete all radiation "sources" from SimPM.RS, to avoid
  // allocating memory for column densities in the cell-data.
  // *****************************************************
  reset_radiation_sources(&(SimPM.RS));

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

  //
  // This code needs 3d data to project...
  //
  if (SimPM.ndim!=3) rep.error("projection needs 3D data to work on",SimPM.ndim);


  //
  // If doing MHD we may want to project the field components, but
  // definitely NOT if running hydro only!
  //
  int SIMeqns = SimPM.eqntype;
  if      (SIMeqns==EQEUL ||
	   SIMeqns==EQEUL_ISO) SIMeqns=1;
  else if (SIMeqns==EQMHD ||
	   SIMeqns==EQGLM ||
	   SIMeqns==EQFCD) SIMeqns=2;
  else rep.error("Bad equations",SIMeqns);

  //
  // Now setup microphysics and raytracing classes
  //
  err += SimSetup->setup_microphysics();
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics and raytracing",err);

  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* Grid setup,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Image Setup ------------------\n";


  //*******************************************************************
  // Set up output image dimensions, and set up pixel cell lists.
  //*******************************************************************
  //struct image_props image;
  //image.set_image_properties(map_dim, normal, perpdir, theta, zero_angle, Npix_max,
  //			     SimPM.Xmin, SimPM.Xmax, SimPM.Range, grid->DX());

  class image IMG (normal, angle, perpdir, grid);
  
  int npix[3]={-1,-1,-1}, num_pixels=0;
  IMG.get_npix(npix);
  num_pixels = npix[0]*npix[1];
  cout <<"npix = ["<<npix[0]<<", "<<npix[1]<<"] and total="<<num_pixels<<"\n"; 
  //
  // If getting velocity profiles, velocity is the third image dimension.
  //
  npix[2]=Nbins;
  for (int ii=0;ii<3;ii++) vps.npix[ii] = npix[ii];

  cout <<"<----- Setting cell positions in Image ----->\n"; cout.flush();
  IMG.set_cell_positions_in_image();
  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();
  cout <<"<----- Adding cells to pixels...       ----->\n"; cout.flush();
  IMG.add_cells_to_pixels();
  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();
  cout <<"<----- Adding Integration points to px ----->\n"; cout.flush();
  IMG.add_integration_pts_to_pixels();
  mltsf=clk.time_so_far("mainloop");
  cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
  cout.flush();
  cout <<"<----- Finished setting up pixels      ----->\n"; cout.flush();

  //
  // So now we have a list of pixels with their physical 3D box locations,
  // and each pixel has a list of cells that are within its box.
  //

  //
  // Set output file; if multiple files, append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to label output file.
  //
  class image_io imio;
  unsigned int ifile=0;
  string filehandle("");
  string this_outfile("");

  cout <<"--------------- Finished Image Setup/Coordinates ------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";

  //*******************************************************************
  // loop over all files:
  //*******************************************************************


  for (ifile=0; ifile< static_cast<unsigned int>(nfiles); ifile++) {
    cout <<"--------------- Starting Next Loop: ifile="<<ifile<<"------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Reading Simulation data to grid -------\n";
    //*******************
    //* Read Input Data *
    //*******************
    //
    // Get filename with path:
    //
    temp.str("");
    temp <<input_path<<"/"<<*ff;
    string infile = temp.str();
    temp.str("");
    ff++;
    //
    // Read header to get timestep info; 
    // also reset the domain to 1/2 the size in Y and Z (if needed).
    //
    err = dataio.ReadHeader(infile);
#ifdef RESET_DOMAIN
    reset_domain();
#endif
    if (err) rep.error("Didn't read header",err);
    if ( (err=MCMD.decomposeDomain()) !=0) 
      rep.error("Couldn't Decompose Domain!",err);

    cout <<"############ SIMULATION TIME: "<<SimPM.simtime/3.16e7;
    cout <<" yrs for step="<<ifile<<"   ############\n";
    cout.flush();

    //
    // Read data (this reader can read serial or parallel data.
    //
    err = dataio.parallel_read_any_data(infile, ///< file to read from
					grid    ///< pointer to data.
					);
    rep.errorTest("(main) Failed to read data",0,err);
    
    cout <<"--------------- Finished Reading Data  ----------------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Starting Data Analysis ----------------\n";
    //********************
    //* Analyse the Data *
    //********************
    //
    // Sample analysis, calculating mean density in pixel,
    // with hard-coded ion-fraction tracer number.
    //
    // I am going to assume:
    //      SimPM.ftr   = Hydrogen ion fraction
    //      SimPM.ftr+1 = Colour tracer [0=ambient, 1=clump]
    //    cout <<"SimPM.tracertype: "<<SimPM.trtype<<endl;
    //
    //

    //
    // If integrating velocity profiles, I need num_pixels*N_velocity_bins elements to my pixel array.
    // If getting velocity profiles, also want a 2D image flattened in the perp. direction
    //
    // im is a pointer to one of im1/2/3/4/5
    // 
    double *im=0, *im1=0, *im2=0, *im3=0, *im4=0, *im5=0, *im6=0, *im7=0, *im8=0;
    long int nels = num_pixels*Nbins; // Nbins=1 unless we want V_los or V_x

    int n_images=0;
    int *what2int=0;
    double **img_array=0;

    switch (what_to_integrate) {
    case I_DENSITY:
    case I_NEUTRAL_NH:
    case I_EMISSION:
      //
      // Only need one image:
      //
      n_images = 1;
      im1 = mem.myalloc(im1,nels);
      what2int  = mem.myalloc(what2int ,n_images);
      img_array = mem.myalloc(img_array,n_images);
      for (int v=0;v<nels;           v++) im1[v]=0.0;
      what2int[0] = what_to_integrate;
      img_array[0] = im1;
      break;
    case I_VEL_LOS:
    case I_VX:
      n_images = 1;
      im1 = mem.myalloc(im1,nels);
      im2 = mem.myalloc(im2,npix[0]*npix[2]);
      what2int  = mem.myalloc(what2int ,n_images);
      img_array = mem.myalloc(img_array,n_images);
      for (int v=0;v<nels;           v++) im1[v]=0.0;
      for (int v=0;v<npix[0]*npix[2];v++) im2[v]=0.0;
      what2int[0] = what_to_integrate;
      img_array[0] = im1;      
      break;

    case I_ALL_SCALARS:
      if (SIMeqns==1) 
	n_images = 4; // No field components
      else
	n_images = 8; // Project Stokes Q,U and BX,BT
      im1 = mem.myalloc(im1,nels);
      im2 = mem.myalloc(im2,nels);
      im3 = mem.myalloc(im3,nels);
      im4 = mem.myalloc(im4,nels);
      if (SIMeqns==2) { 
	im5 = mem.myalloc(im5,nels);
	im6 = mem.myalloc(im6,nels);
	im7 = mem.myalloc(im7,nels);
	im8 = mem.myalloc(im8,nels);
      }
      what2int  = mem.myalloc(what2int ,n_images);
      img_array = mem.myalloc(img_array,n_images);
      for (int v=0;v<nels; v++)
	im1[v] = im2[v] = im3[v] = im4[v] = 0.0;
      if (SIMeqns==2) { 
	for (int v=0;v<nels; v++)
	  im5[v] = im6[v] = im7[v] = im8[v] = 0.0;
      }
      what2int[0] = I_DENSITY;
      what2int[1] = I_NEUTRAL_NH;
      what2int[2] = I_EMISSION;
      what2int[3] = I_NII6584;
      if (SIMeqns==2) { 
	what2int[4] = I_B_STOKESQ;
	what2int[5] = I_B_STOKESU;
	what2int[6] = I_BXabs;
	what2int[7] = I_BYabs;
      }
      img_array[0] = im1;
      img_array[1] = im2;
      img_array[2] = im3;
      img_array[3] = im4;
      if (SIMeqns==2) { 
	img_array[4] = im5;
	img_array[5] = im6;
	img_array[6] = im7;
	img_array[7] = im8;
      }
      break;
    default:
      rep.error("bad what-to-integrate integer...",what_to_integrate);
    }


    double tot_mass = 0.0;
    struct pixel *px;
    int w2i=-1;

    for (int outputs=0;outputs<n_images;outputs++) {
#ifdef TESTING
      cout <<"starting image "<<outputs<<" calculation.\n";
      cout.flush();
#endif // TESTING
      im  = img_array[outputs];
      w2i =  what2int[outputs];
      tot_mass = 0.0;

      //
      // For each output image, loop over all pixels:
      // either multi-threaded or not...
      //
      clk.start_timer("makeimage"); double tsf=0.0;

#ifdef THREADS
      cout <<"Beginning analysis: NUMTHREADS="<<NUM_THREADS_MAIN<<"... ";
      cout <<"i="<<outputs<<", w2i = "<<w2i<<" ... ";
      //cout.flush();
#endif // THREADS

      for (int i=0;i<num_pixels;i++) {
	px = &(IMG.pix[i]);

#ifndef THREADS
	IMG.calculate_pixel(px,       ///< pointer to pixel
			    &vps,     ///< info for velocity profiling.
			    w2i,      ///< flag for what to integrate.
			    im,       ///< array of pixel data.
			    &tot_mass ///< general purpose counter for stuff.
			    );
#endif // not THREADS

#ifdef THREADS

	struct calc_pix_args *ta = new struct calc_pix_args;
	ta->px = px;
	ta->IMG = &IMG;
	ta->what_to_integrate = w2i;
	ta->vps = &vps;
	ta->im = im;
	ta->tot_mass = &tot_mass;
#ifdef TESTING
        //cout <<" - -- - adding pixel "<<i<<" to work-list.\n";
        //cout.flush();
#endif
	//calculate_pixelW(reinterpret_cast<void *>(ta));
	tp_addWork(&tp,calculate_pixelW,reinterpret_cast<void *>(ta),"main()");
#endif // THREADS
	
      }
#ifdef THREADS
#ifdef TESTING
      cout <<" - -- - waiting for "<<num_pixels<<" threads to finish.\n";
      cout.flush();
#endif
      //DbgMsg(" main(): waiting for %i threads...",num_pixels);
      tp_waitOnFinished(&tp,num_pixels);
      //DbgMsg(" main(): all threads finished.");
#ifdef TESTING
      cout <<" - -- - All threads are finished.\n";
      cout.flush();
#endif
#endif // THREADS
      tsf=clk.time_so_far("makeimage");
      cout <<"\t time = "<<tsf<<" secs."<<"\n";
      //cout.flush();
      clk.stop_timer("makeimage");
    } // loop over output images

    
    // ***************************************************************
    //
    // Replace projected |Bx|,|By| (images 6,7) with values calculated
    // from the Stokes Q and U values in images 4,5.
    //
    if (n_images==8) {
      double norm;
      for (int ix=0;ix<num_pixels;ix++) {
	norm = sqrt(img_array[4][ix]*img_array[4][ix]+
		    img_array[5][ix]*img_array[5][ix]);
	img_array[6][ix] =
	  norm*cos(0.5*atan2(img_array[5][ix],img_array[4][ix]));
	img_array[7][ix] =
	  norm*sin(0.5*atan2(img_array[5][ix],img_array[4][ix]));
      }
    }
    // ***************************************************************

    //
    // Now see if we got all the mass in the simulation domain:
    //  
    tot_mass *= grid->DA();
    //cout <<"\t\tANGLE, TOTAL MASS FROM PROJECTION, SUMMATION: "<<angle<<"\t"<<tot_mass;
    //tot_mass=0;
    //  //double posIMG[3], posSIM[3];
    //cell *c=grid->FirstPt();
    //do {
    //  tot_mass += c->P[RO];
    //  //IMG.get_image_Ipos(c->pos,posIMG);
    //  //IMG.get_sim_Dpos(posIMG, posSIM);
    //  //rep.printVec("CELL POS:",c->pos,3);
    //  //rep.printVec("IMG  POS:",posIMG,3);
    //  //rep.printVec("SIM  POS:",posSIM,3);
    //  //rep.printVec("IMG OOOO:",IMG.s_origin_img,2);
    //} while ( (c=grid->NextPt(c))!=0);
    //tot_mass *= grid->DV();
    //cout <<"\t"<<tot_mass<<endl;

#ifdef SUBTRACT_MEAN
    //
    // Here we want to subtract off the mean density and neutral
    // density from the images, to avoid linear gradients from the
    // cubic domain being projected at an angle.
    // We assume the "top" of the image is upstream undisturbed gas
    // and subtract the values in the top row from all rows below.
    //
    switch (what_to_integrate) {
    case I_DENSITY:
    case I_NEUTRAL_NH:
      //
      // Here we just do the first (and only) image.
      //
      for (int iy=0; iy<npix[1]; iy++)
        for (int ix=0; ix<npix[0]; ix++)
          img_array[0][npix[0]*iy+ix] -= img_array[0][npix[0]*(npix[1]-1)+ix];
      break;
    case I_ALL_SCALARS:
      //
      // Here we need to subtract from the first and second images.
      //
      for (int iy=0; iy<npix[1]; iy++) {
        for (int ix=0; ix<npix[0]; ix++) {
          img_array[0][npix[0]*iy+ix] -= img_array[0][npix[0]*(npix[1]-1)+ix];
          img_array[1][npix[0]*iy+ix] -= img_array[1][npix[0]*(npix[1]-1)+ix];
        }
      }
      break;
    default:
      // default action is to do nothing.
      break;
    }
#endif // SUBTRACT_MEAN

    //
    // If we got a P-V data-cube, also construct a 2D image projected
    // along the perpendicular direction.
    //
    if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
      int ix=0, iy=0, iz=0;
      for (long int v=0;v<nels; v++) {
	im2[npix[0]*iz+ix] += im1[v];
	ix++;
	if (ix>=npix[0]) {ix=0; iy++;}
	if (iy>=npix[1]) {iy=0; iz++;}
      }
    }

    cout <<"--------------- Finished Analysing this step ----------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Writing image and getting next Im-file \n";
    
    //double posIMG[3], posSIM[3];
    //IMG.get_image_Ipos(grid->FirstPt()->pos,posIMG);
    //IMG.get_sim_Dpos(posIMG, posSIM);
    //rep.printVec("CELL POS:",grid->FirstPt()->pos,3);
    //rep.printVec("IMG  POS:",posIMG,3);
    //rep.printVec("SIM  POS:",posSIM,3);
    double im_xmin[3];
    for (int v=0; v<3;v++) {
      im_xmin[v] = 0.0;  //posSIM[v] - (posIMG[v]+0.5)*grid->DX();
    }
    double im_dx[3] = {grid->DX(), grid->DX(), grid->DX()};
    if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
      im_xmin[2] = v_min;
      im_dx[2]   = bin_size;
    }
    rep.printVec("IMG XMIN:",im_xmin,3);
    rep.printVec("IMG DX:  ",im_dx,3);


    //**********************
    //* Write Data to file *
    //**********************
    this_outfile = imio.get_output_filename(outfile, multi_opfiles, op_filetype, ifile);
    err = imio.open_image_file(this_outfile, op_filetype, &filehandle);
    if (err) rep.error("failed to open output file",err);

    string *im_name = mem.myalloc(im_name, n_images);
    ostringstream t; t.fill('0');
    
    switch (what_to_integrate) {
    case I_DENSITY:
      t<<"Proj_Dens";
      im_name[0]=t.str();
      break;
    case I_NEUTRAL_NH:
      t<<"Proj_NH";
      im_name[0]=t.str();
      break;
    case I_VEL_LOS:
      t<<"Proj_LOS_V";
      im_name[0]=t.str();
      break;
    case I_VX:
      t<<"Proj_VX";
      im_name[0]=t.str();
      break;
    case I_EMISSION:
      t<<"Proj_Halpha";
      im_name[0]=t.str();
      break;
    case I_ALL_SCALARS:
      t<<"Proj_Dens";
      im_name[0]=t.str(); t.str("");
      t<<"Proj_NH";
      im_name[1]=t.str(); t.str("");
      t<<"Proj_Halpha";
      im_name[2]=t.str(); t.str("");
      t<<"Proj_NII6584";
      im_name[3]=t.str(); t.str("");
      if (SIMeqns==2) { 
	t<<"Proj_b_q";
	im_name[4]=t.str(); t.str("");
	t<<"Proj_b_u";
	im_name[5]=t.str(); t.str("");
	t<<"Proj_bxabs";
	im_name[6]=t.str(); t.str("");
	t<<"Proj_byabs";
	im_name[7]=t.str(); t.str("");
      }
      break;
    default:
      rep.error("bad what-to-integrate integer...",what_to_integrate);
    }

    //
    // Write N images, depending on what we were asked to output:
    //
    for (int outputs=0;outputs<n_images;outputs++) {
      im  = img_array[outputs];

      switch (what_to_integrate) {
      case I_DENSITY: case I_NEUTRAL_NH: case I_EMISSION: case I_ALL_SCALARS:
        err = imio.write_image_to_file(filehandle, op_filetype, im,
                                      num_pixels, 2, npix,
                                      im_name[outputs],
                                      im_xmin, im_dx,
                                      SimPM.simtime, SimPM.timestep
                                      );
	break;
      case I_VEL_LOS: case I_VX:
        err = imio.write_image_to_file(filehandle, op_filetype, im,
                                      num_pixels*Nbins, 3, npix,
                                      im_name[outputs],
                                      im_xmin, im_dx,
                                      SimPM.simtime, SimPM.timestep
                                      );
	break;
      default:
	rep.error("bad what-to-integrate integer...",what_to_integrate);
      }
      if (err) rep.error("Failed to write image to file",err);
    }

    //
    // Also write a 2D P-V diagram summed along the perpendicular
    // direction, if we have calculted LOS velocity or VX.
    //
    //if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
    //  t.str("");
    //  if (what_to_integrate==I_VEL_LOS)
    //    t<<"los_vel_proj_";
    //  else
    //    t<<"velx_proj";
    //  t.width(5); t<<SimPM.timestep;
    //  im_name[0] = t.str();
    //  int np[2]; np[0]=npix[0]; np[1]=npix[2];
    //  err = imio.write_image_to_file(filehandle, op_filetype, im2,
    //                                  npix[0]*npix[2], 2, np, im_name[0],
    //                                  im_xmin, im_dx,
    //                                  SimPM.simtime, SimPM.timestep
    //                                  );
    //  if (err) rep.error("Failed to write 2nd image to file",err);
    //}
    
    
    //*********************************************
    //* Close outfile if using multiple O/P files *
    //*********************************************
    if (multi_opfiles) {
      err = imio.close_image_file(filehandle);
      if (err) rep.error("failed to close output file",err);
    }
    im1 = mem.myfree(im1);
    im2 = mem.myfree(im2);
    im3 = mem.myfree(im3);
    im4 = mem.myfree(im4);
    im5 = mem.myfree(im5);
    im6 = mem.myfree(im6);
    im7 = mem.myfree(im7);
    im8 = mem.myfree(im8);
    im_name = mem.myfree(im_name);


    mltsf=clk.time_so_far("mainloop");
    cout <<"*-*-*-* Loop: "<< ifile <<",\t total time so far = "<<mltsf<<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
    cout.flush();

  } // Loop over all files.    

  //
  // Close file, if single file.
  //
  if (!multi_opfiles) {
    err = imio.close_image_file(filehandle);
    if (err) rep.error("failed to close output file",err);
  }

  cout <<"--------------- Finised Analysing all Files -----------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";
  //
  // Need to delete extra cell position before deleting grid.
  //
  IMG.delete_cell_positions();

  if(grid!=0) {
    cout << "\t Deleting Grid Data..." << "\n";
    delete grid; grid=0;
  }
  if (MP)     {delete MP; MP=0;}
  if (RT)     {delete RT; RT=0;}

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

///
/// Reset the radiation sources in the header to correspond to projected
/// quantities and not the sources used for the simulation.
///
void reset_radiation_sources(struct rad_sources *rs)
{
  //
  // struct rad_sources {
  //   int Nsources;
  //   std::vector<struct rad_src_info> sources;
  // };
  //

  //
  // delete any current sources
  //
  if (rs->Nsources!=0) {
    rs->sources.clear();
    rs->Nsources=0;
  }
  SimPM.EP.raytracing=0;

  return;
}



// ##################################################################
// ##################################################################





