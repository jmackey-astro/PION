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
/// - 2015.08.19 JM: Changed subtraction of mean to take the mean
///    from the first image and apply that to all subsequent images.
/// - 2015.08.20 JM: Changed image coordinates, so that the origin of
///    the simulation is projected onto the image origin.
/// - 2015.10.13 JM: added 20cm Bremsstrahlung and Emission measure


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
#include "tools/threads_AJL/msvc_constants.h"
#if defined(_DEBUG) &&  defined(_MSC_VER) &&  defined(MSVC_DEBUG_NEW_TRACE_ON)
  #define CRTDBG_MAP_ALLOC
  #include <stdlib.h> 
  #include <crtdbg.h> 
  #define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif
#include "tools/threads_AJL/reefa_constants.h"
#include "tools/threads_AJL/logmessages.h"
#include "tools/threads_AJL/threadpool/threadpool.h"
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
void reset_domain(class MCMDcontrol *MCMD)
{
  rep.printVec("Old Xmin",SimPM.Xmin,SimPM.ndim);
  rep.printVec("Old Xmax",SimPM.Xmax,SimPM.ndim);
  for (int v=0;v<SimPM.ndim;v++) {
    //SimPM.Xmin[v] = -309.686272e18;  // for BT3_v070_r2p5, to give 768^3 grid
    //SimPM.Xmin[v] = -278.085632e18;  // for BT3_v[020/030]_r2p5, to give 768^3 grid
    //SimPM.Xmin[v] = -227.524608e18;  // for BT3_v010_r2p5, to give 768^3 grid
    SimPM.NG[v] = static_cast<int>(ONE_PLUS_EPS*SimPM.NG[v]*(SimPM.Xmax[v]-SimPM.Xmin[v])/SimPM.Range[v]);
    SimPM.Range[v] = SimPM.Xmax[v] - SimPM.Xmin[v];
  }
  rep.printVec("New Xmin",SimPM.Xmin,SimPM.ndim);
  rep.printVec("New Xmax",SimPM.Xmax,SimPM.ndim);
  MCMD->decomposeDomain();
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
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  MCMD.set_myrank(myrank);
  MCMD.set_nproc(nproc);


#ifdef THREADS
  tp_init(&tp,NUM_THREADS_MAIN,"Main Threadpool");
#endif //THREADS

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
    cout <<"\t\twriting data to text file.\n";
    break;
  case 1:
    cout <<"\t\twriting data to fits files.\n";
    break;
  case 3:
    cout <<"\t\twriting data to VTK files.\n";
    break;
  default:
    rep.error("Bad outfile format",op_filetype);
  }

  int multi_opfiles = atoi(argv[5]);
  switch (multi_opfiles) {
  case 0:
    cout <<"\t\twriting all timesteps in a single file.\n";
    break;
  case 1:
    cout <<"\t\twriting timesteps in different files.\n";
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
  //
#ifdef RESET_DOMAIN
  reset_domain(&MCMD);
#endif
  //
  // ------------------------------------------------------------------------
  //

  //
  // First decompose the domain, so I know the dimensions of the local
  // grid to set up.  If nproc==1, then this sets the local domain to
  // be the full domain.
  // Get axis corresponding to perpdir, and decompose only along
  // this axis.
  //
  enum axes perpaxis = static_cast<axes>(static_cast<int>(perpdir)/2);
  cout <<"*** perpendicular axis = "<<perpaxis<<"\n";
  MCMD.decomposeDomain(perpaxis);
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
  // Setup arrays for the pixels in each image.
  //
  // If integrating velocity profiles, I need
  // num_pixels*N_velocity_bins elements to my pixel array.
  //
  // If getting velocity profiles, also want a 2D image flattened in
  // the perp. direction
  //
  // im is a pointer to one of im1/2/3/4/5
  // 
  double *im=0, *im1=0, *im2=0, *im3=0, *im4=0, *im5=0, *im6=0,
         *im7=0, *im8=0, *im9=0, *im10=0, *im11=0;
  long int nels = num_pixels*Nbins; // Nbins=1 unless we want V_los or V_x

  //
  // first processor needs a bigger array for the images if there is
  // more than one MPI process.
  //
  int rank0_npix[3]={-1,-1,-1}, rank0_num_pixels=0;
  if (myrank==0) {
    nels *= nproc;
    rank0_npix[0] = npix[0];
    rank0_npix[1] = npix[1]*nproc;
    rank0_npix[2] = npix[2];
    rank0_num_pixels = rank0_npix[0]*rank0_npix[1];
  }

  int n_images=0;
  int *what2int=0;
  double **img_array=0;
#ifdef SUBTRACT_MEAN
  double **mean_array=0;
#endif // SUBTRACT_MEAN

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
    what2int[0] = what_to_integrate;
    img_array[0] = im1;      
    break;

  case I_ALL_SCALARS:
    if (SIMeqns==1) 
      n_images = 6; // No B-field components (dens, NH0, HA, NII, EM, 20cm)
    else
      n_images = 11; // Project Stokes Q,U and BX,BT, RM
    im1 = mem.myalloc(im1,nels);
    im2 = mem.myalloc(im2,nels);
    im3 = mem.myalloc(im3,nels);
    im4 = mem.myalloc(im4,nels);
    im5 = mem.myalloc(im5,nels);
    im6 = mem.myalloc(im6,nels);
    if (SIMeqns==2) { 
      im7 = mem.myalloc(im7,nels);
      im8 = mem.myalloc(im8,nels);
      im9 = mem.myalloc(im9,nels);
      im10 = mem.myalloc(im10,nels);
      im11 = mem.myalloc(im11,nels);
    }
    what2int  = mem.myalloc(what2int ,n_images);
    img_array = mem.myalloc(img_array,n_images);
    what2int[0] = I_DENSITY;
    what2int[1] = I_NEUTRAL_NH;
    what2int[2] = I_EMISSION;
    what2int[3] = I_NII6584;
    what2int[4] = I_EM;
    what2int[5] = I_BREMS20CM;
    if (SIMeqns==2) { 
      what2int[6] = I_B_STOKESQ;
      what2int[7] = I_B_STOKESU;
      what2int[8] = I_BXabs;
      what2int[9] = I_BYabs;
      what2int[10]= I_RM;
    }
    img_array[0] = im1;
    img_array[1] = im2;
    img_array[2] = im3;
    img_array[3] = im4;
    img_array[4] = im5;
    img_array[5] = im6;
    if (SIMeqns==2) { 
      img_array[6] = im7;
      img_array[7] = im8;
      img_array[8] = im9;
      img_array[9] = im10;
      img_array[10]= im11;
    }
#ifdef SUBTRACT_MEAN
    //
    // allocate memory for top row of the first image, so that we
    // can subtract that from all subsequent images.  Only do this
    // for the two images of projected density and neutral density.
    //
    mean_array = mem.myalloc(mean_array,2);
    for (int v=0;v<2;v++) {
      mean_array[v] = mem.myalloc(mean_array[v],npix[0]);
    }
#endif // SUBTRACT_MEAN
    break;
  default:
    rep.error("bad what-to-integrate integer...",what_to_integrate);
  }

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
#ifdef TESTING
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Reading Simulation data to grid -------\n";
#endif
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
    reset_domain(&MCMD);
#endif
    if (err) rep.error("Didn't read header",err);
    if ( (err=MCMD.decomposeDomain(perpaxis)) !=0) 
      rep.error("Couldn't Decompose Domain!",err);

    cout <<"############ SIMULATION TIME: "<<SimPM.simtime/3.156e7;
    cout <<" yrs for step="<<ifile<<"   ############\n";
    cout.flush();

    //
    // Read data (this reader can read serial or parallel data.
    //
    err = dataio.parallel_read_any_data(infile, ///< file to read from
					grid    ///< pointer to data.
					);
    rep.errorTest("(main) Failed to read data",0,err);
    
#ifdef TESTING
    cout <<"--------------- Finished Reading Data  ----------------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Starting Data Analysis ----------------\n";
#endif
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
    // Initialize image arrays to zero
    //
    switch (what_to_integrate) {
    case I_DENSITY:
    case I_NEUTRAL_NH:
    case I_EMISSION:
      //
      // Only need one image:
      //
      for (int v=0;v<nels;           v++) im1[v]=0.0;
      break;
    case I_VEL_LOS:
    case I_VX:
      for (int v=0;v<nels;           v++) im1[v]=0.0;
      for (int v=0;v<npix[0]*npix[2];v++) im2[v]=0.0;
      break;

    case I_ALL_SCALARS:
      for (int v=0;v<nels; v++) im1[v] = 0.0;
      for (int v=0;v<nels; v++) im2[v] = 0.0;
      for (int v=0;v<nels; v++) im3[v] = 0.0;
      for (int v=0;v<nels; v++) im4[v] = 0.0;
      for (int v=0;v<nels; v++) im5[v] = 0.0;
      for (int v=0;v<nels; v++) im6[v] = 0.0;
      if (SIMeqns==2) { 
        for (int v=0;v<nels; v++) im7[v] = 0.0;
        for (int v=0;v<nels; v++) im8[v] = 0.0;
        for (int v=0;v<nels; v++) im9[v] = 0.0;
        for (int v=0;v<nels; v++) im10[v]= 0.0;
        for (int v=0;v<nels; v++) im11[v]= 0.0;
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

      // ------------------------------------------------------------
      // ------------------------------------------------------------
      // If nproc>1, then we need to gather all data on process 0.
      //
      if (nproc>1 && myrank==0) {
        //
        // allocate buffer to receive data.
        //
        //cout <<"RANK 0: RECEIVING DATA\n";
        //double *buf =0;
        long int ct = nels/nproc;
        //buf = mem.myalloc(buf,ct);
        //
        // loop over all the other processes to get data from them,
        // but not neccessarily in order.
        //
        for (int irank=1; irank<nproc; irank++) {
          string recv_id;
          int recv_tag=-1;
          int from_rank=-1;
          err = COMM->look_for_data_to_receive(
                       &from_rank, ///< rank of sender
                       recv_id,    ///< identifier for receive.
                       &recv_tag,  ///< comm_tag associated with data.
                       COMM_DOUBLEDATA ///< type of data we want.
                       );
          if (err) rep.error("look for cell data failed",err);

          //
          // Receive data into buffer.
          //
          cout <<"receiving from "<<from_rank<<"  "<<recv_id<<"  "<<recv_tag<<"\n";
          err = COMM->receive_double_data(
                  from_rank, ///< rank of process we are receiving from.
                  recv_tag,  ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
                  recv_id, ///< identifier for receive, for any book-keeping that might be needed.
                  ct, ///< number of doubles to receive
                  &(im[ct*from_rank])
                  //buf ///< Pointer to array to write to (must be already initialised).
                  );
          if (err) {
            cout <<from_rank <<"\t"<< recv_tag <<"\t"<< recv_id;
            cout <<"\t"<< ct <<"\t"<< irank<<"\n";
            rep.error("Receive image getdata failed",err);
          }

          //
          // put data into image array
          //
          //for (long int v=0; v<ct; v++) im[ct*from_rank +v] = buf[v];
        }
        //buf = mem.myfree(buf);
      } // if myrank==0
      else if (myrank!=0 && nproc>1) {
        //
        // send data to rank 0
        //
        cout <<"RANK "<<myrank<<": SENDING DATA\n";
        string id;
        cout <<"sending "<<nels<<" to rank 0.\n";
        cout.flush();
        err = COMM->send_double_data(
              0,       ///< rank to send to.
              nels,    ///< size of buffer, in number of doubles.
              im,      ///< pointer to double array.
              id,      ///< identifier for send, for tracking delivery later.
              BC_RTtag ///< comm_tag, to say what kind of send this is.
              );
        if (err) rep.error("Send image failed.",err);

        err = COMM->wait_for_send_to_finish(id);
        if (err) rep.error("wait for send to finish failed",err);
      } // if myrank !=0
      // ------------------------------------------------------------
      // ------------------------------------------------------------


    } // loop over output images

    
    // ***************************************************************
    //
    // Replace projected |Bx|,|By| (images 6,7) with values calculated
    // from the Stokes Q and U values in images 4,5.
    //
    if (n_images==9) {
      double norm;
      for (int ix=0;ix<num_pixels;ix++) {
	norm = sqrt(img_array[6][ix]*img_array[6][ix]+
		    img_array[7][ix]*img_array[7][ix]);
	img_array[8][ix] =
	  norm*cos(0.5*atan2(img_array[7][ix],img_array[6][ix]));
	img_array[9][ix] =
	  norm*sin(0.5*atan2(img_array[7][ix],img_array[6][ix]));
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
    // Again, only do this for rank 0 with the gathered image.
    //
    if (myrank==0) {
      switch (what_to_integrate) {
      case I_DENSITY:
      case I_NEUTRAL_NH:
        //
        // Here we just do the first (and only) image.
        //
        for (int iy=0; iy<rank0_npix[1]; iy++)
          for (int ix=0; ix<rank0_npix[0]; ix++)
            img_array[0][rank0_npix[0]*iy+ix] -= img_array[0][rank0_npix[0]*(rank0_npix[1]-1)+ix];
        break;
      case I_ALL_SCALARS:
        //
        // if we are on the first image, then we need to get the top
        // row of pixels and store them in mean_array[img][x-pix]
        //
        if (ifile==0) {
          for (int ix=0; ix<rank0_npix[0]; ix++) {
            mean_array[0][ix] = img_array[0][rank0_npix[0]*(rank0_npix[1]-1)+ix];
            mean_array[1][ix] = img_array[1][rank0_npix[0]*(rank0_npix[1]-1)+ix];
          }
          //rep.printVec("rho",mean_array[0],npix[0]);
          //rep.printVec("NH ",mean_array[1],npix[0]);
        }

        //
        // Here we need to subtract from the first and second images.
        //
        for (int iy=0; iy<rank0_npix[1]; iy++) {
          for (int ix=0; ix<rank0_npix[0]; ix++) {
            //if (!isfinite(img_array[0][npix[0]*iy+ix]) ||
            //    !isfinite(mean_array[0][ix])) {
            //  cout <<"not finite! "<<img_array[0][npix[0]*iy+ix];
            //  cout<<"  "<<mean_array[0][ix]<<"\n";
            //}
            img_array[0][rank0_npix[0]*iy+ix] -= mean_array[0][ix];
            img_array[1][rank0_npix[0]*iy+ix] -= mean_array[1][ix];
            //img_array[0][npix[0]*iy+ix] -= img_array[0][npix[0]*(npix[1]-1)+ix];
            //img_array[1][npix[0]*iy+ix] -= img_array[1][npix[0]*(npix[1]-1)+ix];
          }
        }
        break;
      default:
        // default action is to do nothing.
        break;
      }
    } // if myrank==0
#endif // SUBTRACT_MEAN

    //
    // If we got a P-V data-cube, also construct a 2D image projected
    // along the perpendicular direction.
    //
    if (myrank==0) {
      if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
        int ix=0, iy=0, iz=0;
        for (long int v=0;v<nels; v++) {
          im2[rank0_npix[0]*iz+ix] += im1[v];
          ix++;
          if (ix>=rank0_npix[0]) {ix=0; iy++;}
          if (iy>=rank0_npix[1]) {iy=0; iz++;}
        }
      }
    } // if myrank==0

#ifdef TESTING
    cout <<"--------------- Finished Analysing this step ----------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Writing image and getting next Im-file \n";
#endif
    
    //double posIMG[3], posSIM[3];
    //IMG.get_image_Ipos(grid->FirstPt()->pos,posIMG);
    //IMG.get_sim_Dpos(posIMG, posSIM);
    //rep.printVec("CELL POS:",grid->FirstPt()->pos,3);
    //rep.printVec("IMG  POS:",posIMG,3);
    //rep.printVec("SIM  POS:",posSIM,3);
    //
    // Want to set the image origin to project onto the simulation
    // origin.  This is a clunky way to find it, but it works...
    // Rank 0 always has local Xmin equal to global Xmin, so it works
    // for multiple cores.
    //
    double im_xmin[3], o2[3];
    pion_flt origin[3];
    for (int v=0; v<3;v++) {
      im_xmin[v] = 0.0;  //posSIM[v] - (posIMG[v]+0.5)*grid->DX();
      origin[v] = 0.0;
      o2[v] = 0.0;
    }
    CI.get_ipos_as_double(o2,o2);
    for (int v=0; v<3;v++) origin[v]=o2[v];
    IMG.get_image_Dpos(origin,origin);
    for (int v=0; v<3;v++) im_xmin[v] = -origin[v]*grid->DX();
    rep.printVec("sim origin in units of dx",origin,3);

    double im_dx[3] = {grid->DX(), grid->DX(), grid->DX()};
    if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
      im_xmin[2] = v_min;
      im_dx[2]   = bin_size;
    }
//#ifdef TESTING
    rep.printVec("IMG XMIN:",im_xmin,3);
    rep.printVec("IMG DX:  ",im_dx,3);
//#endif // TESTING


    //**********************
    //* Write Data to file *
    //**********************
    // only proc 0 writes data!!!

    if (myrank==0) {
      //this_outfile = imio.get_output_filename(outfile, multi_opfiles, op_filetype, ifile);
      this_outfile = imio.get_output_filename(outfile, multi_opfiles, op_filetype, SimPM.timestep);
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
        t<<"Proj_EM";
        im_name[4]=t.str(); t.str("");
        t<<"Proj_BREMS20CM";
        im_name[5]=t.str(); t.str("");
        if (SIMeqns==2) { 
          t<<"Proj_b_q";
          im_name[6]=t.str(); t.str("");
          t<<"Proj_b_u";
          im_name[7]=t.str(); t.str("");
          t<<"Proj_bxabs";
          im_name[8]=t.str(); t.str("");
          t<<"Proj_byabs";
          im_name[9]=t.str(); t.str("");
          t<<"Proj_RM";
          im_name[10]=t.str(); t.str("");
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
                                        rank0_num_pixels, 2, rank0_npix,
                                        im_name[outputs],
                                        im_xmin, im_dx,
                                        SimPM.simtime, SimPM.timestep
                                        );
          break;
        case I_VEL_LOS: case I_VX:
          err = imio.write_image_to_file(filehandle, op_filetype, im,
                                        rank0_num_pixels*Nbins, 3, rank0_npix,
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
      im_name = mem.myfree(im_name);
    } // if myrank==0

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

  im1 = mem.myfree(im1);
  im2 = mem.myfree(im2);
  im3 = mem.myfree(im3);
  im4 = mem.myfree(im4);
  im5 = mem.myfree(im5);
  im6 = mem.myfree(im6);
  im7 = mem.myfree(im7);
  im8 = mem.myfree(im8);
  im9 = mem.myfree(im9);
  im10= mem.myfree(im10);
  im11= mem.myfree(im11);
  img_array=mem.myfree(img_array);
  what2int=mem.myfree(what2int);
#ifdef SUBTRACT_MEAN
  mean_array[0]=mem.myfree(mean_array[0]);
  mean_array[1]=mem.myfree(mean_array[1]);
  mean_array=mem.myfree(mean_array);
#endif // SUBTRACT_MEAN
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





