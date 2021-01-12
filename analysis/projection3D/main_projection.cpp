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
/// - 2015.10.13 JM: added 6GHz Bremsstrahlung and Emission measure


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
//#define SUBTRACT_MEAN

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

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"

#include "grid/uniform_grid.h"

#include "sim_projection.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_grid_NG_MPI.h"


#include "raytracing/raytracer_SC.h"

// these files are in ../projection.
#include "../projection/image_io.h"
#include "../xray/xray_emission.h"
#include "../projection/projection_constants.h"


///
/// Reset the radiation sources in the header to correspond to projected
/// quantities and not the sources used for the simulation.
///
void reset_radiation_sources(
      class SimParams &SimPM
      );

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
  SimPM.grid_nlevels = 1;
  MCMD->decomposeDomain(SimPM,SimPM.levels[0]);
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
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  cout <<"Projection3D: myrank="<<myrank<<", nproc="<<nproc<<"\n";
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);


  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc<11) {
    cout << "Use as follows:\n";
    cout << "projection: <projection> <input-path> <input-silo-file-base>\n";
    cout << "\t\t <output-file> <op-file-type> <multi-opfiles> \n";
    cout << "\t\t <normalvec> <fixed_dir> <theta> <what2integrate> ";
    cout << " <skip>";
    cout <<" +[optional velocity args] \n";
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
    cout <<"                         5=StokesQ,6=StokesU, 7=All-scalars, 8=|B|(LOS), 9=|B|(perp)\n";
    cout <<"skip:         will skip this number of input files each loop. ";
    cout <<"(0 means it will calculate every file)\n";
    cout <<"OPTIONAL VELOCITY ARGS:\n";
    cout <<"Nbins:        integer number of bins in velocity profile.\n";
    cout <<"v_min:        minimum velocity to measure (in same units as output files!)\n";
    cout <<"v_max:        maximum velocity to measure (in same units as output files!)\n";
    cout <<"smooth:       integer =1 for constant broadening, =2 for doppler broadening by temperature.\n";
    cout <<"smoooth-val:  float for the velocity by which to smooth, if constant broadening (FWHM)\n";
    rep.error("Bad number of args",argc);
  }


  cout <<"argv[0] = command = "<< argv[0] <<"\n";
  cout <<"argv[1] = input path = "<< argv[1] <<"\n";
  cout <<"argv[2] = input file = "<< argv[2] <<"\n";
  cout <<"argv[3] = image file = "<< argv[3] <<"\n";
  cout <<"argv[4] = image filetype = "<< argv[4] <<"  (3=vtk)\n";
  cout <<"argv[5] = multi_opfiles = "<< argv[5] <<"\n";
  cout <<"argv[6] = Normal dir = "<< argv[6] <<"\n";
  cout <<"argv[7] = Perpendicular dir = "<< argv[7] <<"\n";
  cout <<"argv[8] = Angle to normal = "<< argv[8] <<"\n";
  cout <<"argv[9] = what-to-integrate = "<< argv[9] <<"\n";
  cout <<"argv[10] = skip = "<< argv[10] <<"\n";
/*
  cout <<"argv[] =  = "<< argv[] <<"\n";
*/

  string input_path = argv[1];
  string input_file = argv[2];
  string outfile    = argv[3];

  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  //rep.redirect(redir.str());
  rep.kill_stdout_from_other_procs(0);

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
  if (th<0 || isnan(th) || isinf(th))
    rep.error("Input what-to-integrate is not in [0,1,2,3,4,7,8,9]",th);
  int what_to_integrate=th;
  // 0=density, 1=neutral density, 2=LOS velocity, 3=VX, 4=emission,
  // 7=all_scalars, [5,6]=stokesQU, [8,9]=|B|(LOS,PERP)

  // how sparsely to sample the data files.
  size_t skip = static_cast<size_t>(atoi(argv[10]));

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
    Nbins = atoi(argv[11]);
    v_min = atof(argv[12]);
    v_max = atof(argv[13]);
    smooth= atoi(argv[14]);
    broadening=atof(argv[15]);
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
  class dataio_silo_utility 
                      dataio(SimPM,"DOUBLE",&(SimPM.levels[0].MCMD));

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
      s=files.begin();
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
  err = dataio.ReadHeader(first_file, SimPM);
  if (err) rep.error("Didn't read header",err);

  //
  // ------------------------------------------------------------------------
  // MODIFYING XMIN/XMAX SO THAT I ONLY READ IN A SUBDOMAIN OF THE FULL GRID
  //
#ifdef RESET_DOMAIN
  reset_domain(&(SimPM.levels[0].MCMD));
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
  class setup_grid_NG_MPI *SimSetup =0;
  SimSetup = new setup_grid_NG_MPI();
  SimSetup->setup_NG_grid_levels(SimPM);
  // have to re-do the domain decomposition because we only split
  // the domain on one axis.
  enum axes perpaxis = static_cast<axes>(static_cast<int>(perpdir)/2);
  cout <<"*** perpendicular axis = "<<perpaxis<<"\n";
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    SimPM.levels[l].MCMD.decomposeDomain(perpaxis,SimPM,SimPM.levels[l]);
  }
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    SimPM.levels[l].MCMD.set_NG_hierarchy(SimPM,l);
  }
  if (err) rep.error("main: failed to decompose domain!",err);
  // setup grids
  vector<class GridBaseClass *> G;
  G.resize(SimPM.grid_nlevels);
  SimSetup->setup_grid(G, SimPM);
  class GridBaseClass *grid = G[0];
  if (!grid) rep.error("Grid setup failed",grid);
  SimPM.dx = grid->DX();
  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;

  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables (here just set it to zero).
  //
  // *****************************************************
  // Now delete all radiation "sources" from SimPM.RS, to avoid
  // allocating memory for column densities in the cell-data.
  // *****************************************************
  reset_radiation_sources(SimPM);

  // This code needs 3d data to project...
  if (SimPM.ndim!=3) rep.error("projection needs 3D data",SimPM.ndim);

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
  err += SimSetup->setup_microphysics(SimPM);
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics and raytracing",err);
  
  // Assign boundary conditions to boundary points.
  //err = SimSetup->boundary_conditions(grid);
  //if (err) rep.error("boundary_conditions setup",err);

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

  class image *IMG[SimPM.grid_nlevels];
  for (size_t v=0;v<SimPM.grid_nlevels;v++) {
    IMG[v] = new class image (normal, angle, perpdir, G[v]);
  }

  int npix[3]={-1,-1,-1}, num_pixels=0;
  IMG[0]->get_npix(npix);    // all levels are topologically the same.
  num_pixels = npix[0]*npix[1];
  cout <<"npix = ["<<npix[0]<<", "<<npix[1]<<"] and total="<<num_pixels<<"\n"; 

  // If getting velocity profiles, velocity is the third image dim.
  npix[2]=Nbins;
  for (int ii=0;ii<3;ii++) vps.npix[ii] = npix[ii];

  // setup rays for each level
  for (size_t v=0;v<SimPM.grid_nlevels;v++) {
    cout <<"<----- LEVEL "<<v<<": setting up rays        ----->\n";
    cout <<"<----- Setting cell positions in Image ----->\n";
    cout.flush();
    IMG[v]->set_cell_positions_in_image();
    mltsf=clk.time_so_far("mainloop");
    cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf;
    cout <<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
    cout.flush();
    cout <<"<----- Adding cells to pixels...       ----->\n";
    cout.flush();
    IMG[v]->add_cells_to_pixels();
    mltsf=clk.time_so_far("mainloop");
    cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf;
    cout <<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
    cout.flush();
    cout <<"<----- Adding Integration points to px ----->\n";
    cout.flush();
    IMG[v]->add_integration_pts_to_pixels();
    mltsf=clk.time_so_far("mainloop");
    cout <<"*-*-*-* ... ,\t total time so far = "<<mltsf;
    cout <<" secs or "<<mltsf/3600.0<<" hours. *-*-*-*\n";
    cout.flush();
    cout <<"<----- Finished setting up pixels      ----->\n"; cout.flush();
  }

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
  double *im=0, *im1=0, *im2=0;
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
  case I_D:
  case I_NtD:
  case I_HA:
    // Only need one image:
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
      n_images = N_HD_SCALAR; // No B-field components (dens, NH0, HA, NII, EM, 6GHz, xray)
    else
      n_images = N_MHD_SCALAR; // Project Stokes Q,U and BX,BT, RM
    
    img_array = mem.myalloc(img_array,n_images);
    for (int v=0;v<n_images;v++)
      img_array[v] = mem.myalloc(img_array[v],nels);

    what2int  = mem.myalloc(what2int ,n_images);
    what2int[PROJ_D] = I_D;
    what2int[PROJ_NtD] = I_NtD;
    what2int[PROJ_InD] = I_InD;
    what2int[PROJ_HA] = I_HA;
    what2int[PROJ_NII] = I_NII6584;
    what2int[PROJ_EM] = I_EM;
    what2int[PROJ_BREMS6GHZ] = I_BREMS6GHZ;
    what2int[PROJ_X00p1] = I_X00p1;
    what2int[PROJ_X00p2] = I_X00p2;
    what2int[PROJ_X00p3] = I_X00p3;
    what2int[PROJ_X00p5] = I_X00p5;
    what2int[PROJ_X01p0] = I_X01p0;
    what2int[PROJ_X02p0] = I_X02p0;
    what2int[PROJ_X05p0] = I_X05p0;
    what2int[PROJ_X10p0] = I_X10p0;
    if (SIMeqns==2) { 
      what2int[PROJ_B_STOKESQ] = I_B_STOKESQ;
      what2int[PROJ_B_STOKESU] = I_B_STOKESU;
      what2int[PROJ_BXabs] = I_BXabs;
      what2int[PROJ_BYabs] = I_BYabs;
      what2int[PROJ_ROTNMEASURE]= I_ROTNMEASURE;
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

  // set up master image vectors for each level and each variable.
  cout <<"allocating imgmaster array... ";
  vector< vector< double *> > imgmaster;
  imgmaster.resize(n_images);
  for (int im=0;im<n_images;im++) {
    imgmaster[im].resize(SimPM.grid_nlevels);
    for (int lv=0;lv<SimPM.grid_nlevels;lv++) {
      imgmaster[im][lv] = mem.myalloc(imgmaster[im][lv],nels);
    }
  }
  cout <<"done!\n";
  

  cout <<"--------------- Finished Image Setup/Coordinates ------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";

  //*******************************************************************
  // loop over all files:
  //*******************************************************************


  for (ifile=0; ifile< static_cast<unsigned int>(nfiles); ifile+=1+skip) {
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
    for (size_t q=0;q<=skip;q++) ff++;

    // Read header to get timestep info; 
    // also reset the domain to 1/2 the size in Y and Z (if needed).
    err = dataio.ReadHeader(infile, SimPM);
#ifdef RESET_DOMAIN
    reset_domain(&(SimPM.levels[0].MCMD));
#endif
    if (err) rep.error("Didn't read header",err);

    cout <<"############ SIMULATION TIME: "<<SimPM.simtime/3.156e7;
    cout <<" yrs for step="<<ifile<<"   ############\n";
    cout.flush();

    // Read data (this reader can read serial or parallel data.
    err = dataio.ReadData(infile, G, SimPM);
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
    cout <<"setting images to zero.\n";
    for (int im=0;im<n_images;im++) {
      for (int lv=0;lv<SimPM.grid_nlevels;lv++) {
        for (int v=0;v<nels; v++) imgmaster[im][lv][v] = 0.0;
      }
    }
    
    // Loop over levels to save images for each level
    for (size_t lv=0; lv<SimPM.grid_nlevels; lv++) {
      cout <<"analysis for level "<<lv<<" starting \n";
      grid = G[lv];
      //cout <<"grid pointer="<<grid<<"\n";
  
      switch (what_to_integrate) {
      case I_D:
      case I_NtD:
      case I_HA:
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
        for (int j=0;j<n_images;j++) {
          for (int v=0;v<nels; v++) img_array[j][v] = 0.0;
        }
        break;

      default:
        rep.error("bad what-to-integrate integer...",what_to_integrate);
      }
      //cout <<"set image to zero...\n";
      
      double tot_mass = 0.0;
      struct pixel *px;
      int w2i=-1;
      cell *c = grid->FirstPt();
      // For Cartesian 3D grid, all cell faces have the same area:
      double cell_area = grid->CellInterface(c,XP,0.0);

      // loop over images.
      for (int outputs=0;outputs<n_images;outputs++) {
#ifdef TESTING
        cout <<"starting image "<<outputs<<" calculation.\n";
        cout.flush();
#endif // TESTING
        im  = img_array[outputs];
        w2i =  what2int[outputs];
        tot_mass = 0.0;

        // For each output image, loop over all pixels:
        // either multi-threaded or not...
        clk.start_timer("makeimage"); double tsf=0.0;

        for (int i=0;i<num_pixels;i++) {
          px = &(IMG[lv]->pix[i]);
          IMG[lv]->calculate_pixel(px,&vps,w2i,SimPM,im,&tot_mass);
        }
        tsf=clk.time_so_far("makeimage");
#ifdef TESTING
        cout <<"\t time = "<<tsf<<" secs."<<"\n";
        cout.flush();
#endif // TESTING
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
          // in any order.
          //
          for (int irank=1; irank<nproc; irank++) {
            string recv_id;
            int recv_tag=-1;
            int from_rank=-1;
            int comm_tag = irank;
            err = COMM->look_for_data_to_receive(
                         &from_rank, ///< rank of sender
                         recv_id,    ///< identifier for receive.
                         &recv_tag,  ///< comm_tag associated with data.
                         BC_RTtag,
                         COMM_DOUBLEDATA ///< type of data we want.
                         );
            if (err) rep.error("look for cell data failed",err);

            //
            // Receive data into buffer.
            //
            // cout <<"receiving from "<<from_rank<<"  ";
            // cout <<recv_id<<"  "<<recv_tag<<"\n";
            err = COMM->receive_double_data(
                    from_rank, recv_tag, recv_id, ct, &(im[ct*from_rank])
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
          //cout <<"sending "<<nels<<" to rank 0.\n";
          //cout.flush();
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
        COMM->barrier("outputs");

      } // loop over output images

      
      // ***************************************************************
      //
      // Replace projected |Bx|,|By| (images 6,7) with values calculated
      // from the Stokes Q and U values in images 4,5.
      //
      if (n_images==N_MHD_SCALAR) {
        double norm;
        for (int ix=0;ix<num_pixels;ix++) {
          norm = sqrt(img_array[PROJ_B_STOKESQ][ix]*img_array[PROJ_B_STOKESQ][ix]+
                      img_array[PROJ_B_STOKESU][ix]*img_array[PROJ_B_STOKESU][ix]);
          img_array[PROJ_BXabs][ix] =
            norm*cos(0.5*atan2(img_array[PROJ_B_STOKESU][ix],img_array[PROJ_B_STOKESQ][ix]));
          img_array[PROJ_BYabs][ix] =
            norm*sin(0.5*atan2(img_array[PROJ_B_STOKESU][ix],img_array[PROJ_B_STOKESQ][ix]));
        }
      }
      // ***************************************************************

      //
      // Now see if we got all the mass in the simulation domain:
      //  
      tot_mass *= cell_area;
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
        case I_D:
        case I_NtD:
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
              mean_array[IMG_DENSITY][ix] = img_array[IMG_DENSITY][rank0_npix[0]*(rank0_npix[1]-1)+ix];
              mean_array[IMG_NtD][ix] = img_array[IMG_NtD][rank0_npix[0]*(rank0_npix[1]-1)+ix];
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
              img_array[IMG_DENSITY][rank0_npix[0]*iy+ix] -= mean_array[IMG_DENSITY][ix];
              img_array[IMG_NtD][rank0_npix[0]*iy+ix] -= mean_array[IMG_NtD][ix];
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

      // save image data for this level in imgmaster.
      for (int im=0;im<n_images;im++) {
        for (int v=0;v<nels; v++) imgmaster[im][lv][v] = img_array[im][v];
      }

    } // loop over levels.

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
    grid = G[0];
    double im_xmin[3], o2[3];
    pion_flt origin[3];
    for (int v=0; v<3;v++) {
      im_xmin[v] = 0.0;  //posSIM[v] - (posIMG[v]+0.5)*grid->DX();
      origin[v] = 0.0;
      o2[v] = 0.0;
    }
    CI.get_ipos_as_double(o2,o2);
    for (int v=0; v<3;v++) origin[v]=o2[v];
    IMG[0]->get_image_Dpos(origin,origin);
    for (int v=0; v<3;v++) im_xmin[v] = -origin[v]*grid->DX();
#ifdef TESTING
    rep.printVec("sim origin in units of dx",origin,3);
#endif // TESTING

    grid = G[SimPM.grid_nlevels -1];
    double im_dx[3] = {grid->DX(), grid->DX(), grid->DX()};
    if (what_to_integrate==I_VEL_LOS || what_to_integrate==I_VX) {
      im_xmin[2] = v_min;
      im_dx[2]   = bin_size;
    }
#ifdef TESTING
    rep.printVec("IMG XMIN:",im_xmin,3);
    rep.printVec("IMG DX:  ",im_dx,3);
#endif // TESTING


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
      case I_D:
        t<<"Proj_SurfaceMass";
        im_name[0]=t.str();
        break;
      case I_NtD:
        t<<"Proj_NeutralDens";
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
      case I_HA:
        t<<"Proj_Halpha";
        im_name[0]=t.str();
        break;
      case I_ALL_SCALARS:
        for (size_t im=0; im<n_images; im++) {
          switch (im) {
            case PROJ_D:   im_name[im] = "Proj_SurfaceMass"; break;
            case PROJ_NtD: im_name[im] = "Proj_NeutralDens"; break;
            case PROJ_InD: im_name[im] = "Proj_IonizedDens"; break;
            case PROJ_EM:  im_name[im] = "Proj_EmissionMeasure"; break;
            case PROJ_X00p1: im_name[im] = "Proj_XRAY_g00p1keV"; break;
            case PROJ_X00p2: im_name[im] = "Proj_XRAY_g00p2keV"; break;
            case PROJ_X00p3: im_name[im] = "Proj_XRAY_g00p3keV"; break;
            case PROJ_X00p5: im_name[im] = "Proj_XRAY_g00p5keV"; break;
            case PROJ_X01p0: im_name[im] = "Proj_XRAY_g01p0keV"; break;
            case PROJ_X02p0: im_name[im] = "Proj_XRAY_g02p0keV"; break;
            case PROJ_X05p0: im_name[im] = "Proj_XRAY_g05p0keV"; break;
            case PROJ_X10p0: im_name[im] = "Proj_XRAY_g10p0keV"; break;
            case PROJ_HA:   im_name[im] = "Proj_Halpha"; break;
            case PROJ_NII:  im_name[im] = "Proj_NII_ll6584"; break;
            case PROJ_BREMS6GHZ: im_name[im] = "Proj_BREMS6GHZ"; break;
            case PROJ_B_STOKESQ: im_name[im] = "Proj_B_STOKESQ"; break;
            case PROJ_B_STOKESU: im_name[im] = "Proj_B_STOKESU"; break;
            case PROJ_BXabs:     im_name[im] = "Proj_BXabs"; break;
            case PROJ_BYabs:     im_name[im] = "Proj_BYabs"; break;
            case PROJ_ROTNMEASURE: im_name[im] = "Proj_ROTNMEASURE"; break;
            default: rep.error("Bad image count",im); break;
          }
        }
        break;

      default:
        rep.error("bad what-to-integrate integer...",what_to_integrate);
      }

      //
      // Full image at the resolution of the finest image.
      //
      int gnpix[3]; // number of pixels for full img
      size_t gnumpix=rank0_num_pixels; // total pix count full img
      int ipx=1; // num level-pixels per full-img pixel.
      for (int i=0;i<3;i++) gnpix[i] = rank0_npix[i];
      for (int lv=1;lv<SimPM.grid_nlevels;lv++) {
        gnumpix *= 4;
        ipx*=2;
        for (int i=0;i<3;i++) gnpix[i] *= 2;
      }
      // Add root grid.
      double *global_image = mem.myalloc(global_image,gnumpix);
      //
      // Loop over all images.
      //
      for (int outputs=0;outputs<n_images;outputs++) {
        // Make a global image by summing results from all levels.
        im = imgmaster[outputs][0];
#ifdef TESTING
        cout <<"ipx="<<ipx<<", npix= ["<<rank0_npix[0]<<", ";
        cout <<rank0_npix[1]<<"], \n";
        cout <<"big-img pix= ["<<gnpix[0]<<", "<<gnpix[1];
        cout <<"], tot="<<gnumpix<<"\n";
#endif // TESTING
        // level 0 first: populate the grid.
        for (int j=0; j<rank0_npix[1]; j++) {
          for (int i=0; i<rank0_npix[0]; i++) {
            for (int ky=0; ky<ipx; ky++) {
              for (int kx=0; kx<ipx; kx++) {
                global_image[rank0_npix[0]*ipx*(ipx*j+ky)+ipx*i+kx]
                    = im[rank0_npix[0]*j+i];
              }
            }
          }
        }
        int sz=ipx;
        for (int lv=1; lv<SimPM.grid_nlevels; lv++) {
#ifdef TESTING
          cout <<"populating level "<<lv<<" data onto image.\n";
#endif // TESTING
          // find lower-left corner of nested grid.
          grid = G[lv];
          sz /= 2;
          double lv_xmin[3];
          int corner[3];
          im = imgmaster[outputs][lv];
          for (int v=0; v<3;v++) origin[v]=o2[v];
          IMG[lv]->get_image_Dpos(origin,origin);
          for (int v=0; v<3;v++) lv_xmin[v] = -origin[v]*grid->DX();
#ifdef TESTING
          rep.printVec("sim origin in units of dx",origin,3);
          rep.printVec("level origin in units of dx",origin,3);
          rep.printVec("sim xmin",im_xmin,3);
          rep.printVec("level xmin",lv_xmin,3);
#endif // TESTING
          for (int v=0; v<3;v++) {
            corner[v] = static_cast<int>(round(
                              (lv_xmin[v]-im_xmin[v])*sz*ONE_PLUS_EPS/grid->DX()));
            //corner[v] = static_cast<int>(
            //                  (lv_xmin[v]-im_xmin[v])*sz*ONE_PLUS_EPS/grid->DX());
            //cout<<(lv_xmin[v]-im_xmin[v])*sz*ONE_PLUS_EPS/grid->DX() <<" , ";
          }
          //cout <<"\n";
#ifdef TESTING
          rep.printVec("corner for level",corner,3);
#endif // TESTING
          for (int j=0; j<rank0_npix[1]; j++) {
            for (int i=0; i<rank0_npix[0]; i++) {
              for (int ky=0; ky<sz; ky++) {
                for (int kx=0; kx<sz; kx++) {
                  global_image[rank0_npix[0]*ipx*(corner[1]+sz*j+ky)+corner[0]+sz*i+kx]
                      += im[rank0_npix[0]*j+i];
                  //cout <<ct<<", "<<j<<", "<<i<<", ";
                  //cout <<rank0_npix[0]*ipx*(ipx*j+ky)+ipx*i+kx;
                  //cout <<", "<< rank0_npix[0]*j+i<<"\n";
                  //ct++;
                }
              }
            }
          }
        } // loop over levels 1->n-1

        switch (what_to_integrate) {
        case I_D: case I_NtD: case I_HA: case I_ALL_SCALARS:
          err = imio.write_image_to_file(
              filehandle, op_filetype, global_image,
              static_cast<long int>(gnumpix), 2,
              gnpix,im_name[outputs],im_xmin, im_dx,
              SimPM.simtime, SimPM.timestep);
          break;
        case I_VEL_LOS: case I_VX:
          rep.error("Code no longer works for PPV datacubes",1);
          err = imio.write_image_to_file(
              filehandle, op_filetype, im, rank0_num_pixels*Nbins,
              3, rank0_npix, im_name[outputs], im_xmin, im_dx,
              SimPM.simtime, SimPM.timestep);
          break;
        default:
          rep.error("bad what-to-integrate integer...",what_to_integrate);
        }
        if (err) rep.error("Failed to write image to file",err);
      } // loop over images
      global_image = mem.myfree(global_image);

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

  if (n_images<2 && im)   im = mem.myfree(im);
  if (im1) im1 = mem.myfree(im1);
  if (im2) im2 = mem.myfree(im2);
  if (n_images>1) {
    for (int v=0;v<n_images;v++) {
      //cout <<"v="<<v<<" freeing image "<<img_array[v]<<"\n";
      if (img_array[v]) img_array[v] = mem.myfree(img_array[v]);
    }
  }
  img_array=mem.myfree(img_array);
  what2int=mem.myfree(what2int);
#ifdef SUBTRACT_MEAN
  mean_array[IMG_DENSITY]=mem.myfree(mean_array[IMG_DENSITY]);
  mean_array[1]=mem.myfree(mean_array[1]);
  mean_array=mem.myfree(mean_array);
#endif // SUBTRACT_MEAN

  // Need to delete extra cell position before deleting grids
  for (int i=0; i<SimPM.grid_nlevels; i++)
    IMG[i]->delete_cell_positions();

  //if(grid!=0) {
  //  cout << "\t Deleting Grid Data..." << "\n";
  //  delete grid; grid=0;
  //}
  if (MP)     {delete MP; MP=0;}

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
void reset_radiation_sources(
      class SimParams &SimPM
      )
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
  if (SimPM.RS.Nsources!=0) {
    SimPM.RS.sources.clear();
    SimPM.RS.Nsources=0;
  }
  SimPM.EP.raytracing=0;

  return;
}



// ##################################################################
// ##################################################################





