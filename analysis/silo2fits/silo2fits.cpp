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

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"

#include "MCMD_control.h"
#include "setup_fixed_grid_MPI.h"

#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"
#include "dataIO/dataio_fits.h"

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

#include "raytracing/raytracer_SC.h"


#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"

#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
using namespace std;



// ##################################################################
// ##################################################################


///
/// Reset the radiation sources in the header to correspond to projected
/// quantities and not the sources used for the simulation.
///
void reset_radiation_sources(struct rad_sources *);




// ##################################################################
// ##################################################################


///
/// Resize the domain, so we only read in and write out part of it.
///
void reset_domain(
  const double *xmin,
  const double *xmax,
  const int ndim,
  class MCMDcontrol &MCMD     ///< address of MCMD controller class.
  )
{

  rep.printVec("Old Xmin",SimPM.Xmin, ndim);
  rep.printVec("Old Xmax",SimPM.Xmax, ndim);

  SimPM.Ncell=1;

  for (int v=0; v<ndim; v++) {
    SimPM.Xmin[v]  = xmin[v];
    SimPM.Xmax[v]  = xmax[v];
    SimPM.NG[v]    = static_cast<int>(ONE_PLUS_EPS*SimPM.NG[v]*(SimPM.Xmax[v]-SimPM.Xmin[v])/SimPM.Range[v]);
    SimPM.Ncell *= SimPM.NG[v];
    SimPM.Range[v] = SimPM.Xmax[v] - SimPM.Xmin[v];
  }

  rep.printVec("New Xmin",SimPM.Xmin, ndim);
  rep.printVec("New Xmax",SimPM.Xmax, ndim);

  MCMD.decomposeDomain();
  return;
}



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
  class MCMDcontrol MCMD;
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  MCMD.set_myrank(myrank);
  MCMD.set_nproc(nproc);
  if (nproc>1)
    rep.error("This is serial code",nproc);

  //MP=0; RT=0; grid=0;

  //
  // Get input files and an output file.
  //
  if (argc<5) {
    cout << "Use as follows:\n";
    cout << "<executable-filename> <input-path>";
    cout << " <input-silo-file-base>";
    cout << " <output-path> <output-fits-file-base> [Nskip]";
    cout << " [<xmin> <xmax> [<ymin> <ymax> [<zmin> <zmax>]]]";
    cout << "\n";
    cout <<"******************************************\n";
    cout <<"input path:  path to input files.\n";
    cout <<"input file:  base filename of sequence of files including _0000 if parallel.\n";
    cout <<"output path: directory to write output files to.\n";
    cout <<"output file: filename for output FITS file(s).\n";
    cout <<"[Nskip]:     Only convert every Nskip file.\n";
    cout <<"xmin/xmax:   optional parameters to resize the domain.\n";
    rep.error("Bad number of args",argc);
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
  string op_path    = argv[3];
  string outfile    = argv[4];

  //
  // Only convert every Nskip file
  //
  size_t Nskip = 1;
  if (argc>5) {
    Nskip = static_cast<size_t>(atoi(argv[5]));
  }

  //
  // see about xmin/xmax resetting domain.
  //
  double xmin[MAX_DIM], xmax[MAX_DIM], old_xmin[MAX_DIM], old_xmax[MAX_DIM];
  size_t reset=0;
  if (argc>6)  {xmin[XX] = atof(argv[6]); reset+=1;}
  else          xmin[XX]=0.0;
  if (argc>7)  {xmax[XX] = atof(argv[7]); reset+=1;}
  else          xmax[XX]=0.0;

  if (argc>8)  {xmin[YY] = atof(argv[8]); reset+=1;}
  else          xmin[YY]=0.0;
  if (argc>9)  {xmax[YY] = atof(argv[9]); reset+=1;}
  else          xmax[YY]=0.0;

  if (argc>10)  {xmin[ZZ] = atof(argv[10]); reset+=1;}
  else          xmin[ZZ]=0.0;
  if (argc>11) {xmax[ZZ] = atof(argv[11]); reset+=1;}
  else          xmax[ZZ]=0.0;



  //
  // Redirect output to a text file if you want to:
  //
  ostringstream redir; redir.str("");
  redir<<op_path<<"/msg_"<<outfile<<"_rank"<<MCMD.get_myrank()<<"_";
  //rep.redirect(redir.str());

  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Getting List of Files to read ---------\n";
  
  //
  // set up dataio_utility class and fits-writer class.
  //
  class dataio_silo_utility dataio ("DOUBLE", &MCMD);
  class DataIOFits writer;

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file,  &files);
  if (err) rep.error("failed to get list of files",err);
  //
  // remove non-silo files
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
  size_t nfiles = files.size();
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
  // reset the domain if we need to.
  //
  int dir=0;
  if      (reset==0) dir=0;
  else if (reset==2) dir=1;
  else if (reset==4) dir=2;
  else if (reset==6) dir=3;
  else rep.error("Bad number of args in resetting domain",reset);
  //
  // Save old simulation extents for all directions, for checking
  // later.
  //
  for (size_t i=0;i<static_cast<size_t>(SimPM.ndim);i++) old_xmin[i]=SimPM.Xmin[i];
  for (size_t i=0;i<static_cast<size_t>(SimPM.ndim);i++) old_xmax[i]=SimPM.Xmax[i];
  //
  // leave domains that are not reset as the original extents:
  //
  for (int i=dir;i<SimPM.ndim;i++) xmin[i]=SimPM.Xmin[i];
  for (int i=dir;i<SimPM.ndim;i++) xmax[i]=SimPM.Xmax[i];
  //
  // Now reset domain.
  //
  reset_domain(xmin,xmax,SimPM.ndim,MCMD);

  //
  // write simulation xmin/xmax and radiation source position to a
  // text file.
  //
  redir.str("");
  redir<<op_path<<"/data_"<<outfile<<".txt";
  ofstream outf;
  if (outf.is_open())
    rep.error("Output text file is already open!",1);
  outf.open(redir.str().c_str());
  if (!outf.is_open())
    rep.error("Failed to open text file for writing",redir.str());
  outf.setf( ios_base::scientific );
  outf.precision(6);
  outf <<"Text data for simulation "<<input_file<<"\n\n";
  outf <<"## GRID PROPERTIES ##\n";
  outf <<"Xmin "<<SimPM.Xmin[XX]<<"  cm\n";
  outf <<"Xmax "<<SimPM.Xmax[XX]<<"  cm\n";
  outf <<"Ymin "<<SimPM.Xmin[YY]<<"  cm\n";
  outf <<"Ymax "<<SimPM.Xmax[YY]<<"  cm\n";
  outf <<"Zmin "<<SimPM.Xmin[ZZ]<<"  cm\n";
  outf <<"Zmax "<<SimPM.Xmax[ZZ]<<"  cm\n";
  outf <<"N_X "<<SimPM.NG[XX]<<"\n";
  outf <<"N_Y "<<SimPM.NG[YY]<<"\n";
  outf <<"N_Z "<<SimPM.NG[ZZ]<<"\n";
  outf <<"#\n";
  if (SimPM.RS.Nsources>0) {
    outf <<"## RADIATION SOURCE ##\n";
    outf <<"POS_X "<<SimPM.RS.sources[0].pos[XX]<<"  cm\n";
    outf <<"POS_Y "<<SimPM.RS.sources[0].pos[YY]<<"  cm\n";
    outf <<"POS_Z "<<SimPM.RS.sources[0].pos[ZZ]<<"  cm\n";
    outf <<"Strength "<<SimPM.RS.sources[0].strength<<" erg/s";
    outf <<  "  blackbody source\n";
    outf <<"T_star   "<<SimPM.RS.sources[0].Tstar<<" K. ";
    outf <<  "  Star's effective T.\n";
    outf <<"#\n";
  }
  outf.close();

  // *****************************************************
  // Now delete all radiation "sources" from SimPM.RS, to avoid allocating
  // memory for column densities in the cell-data.
  // *****************************************************
  reset_radiation_sources(&(SimPM.RS));

  for (size_t v=0; v<static_cast<size_t>(SimPM.ndim); v++) {
    cout <<"old_xmin="<<old_xmin[v]<<", SimPM.Xmin="<<SimPM.Xmin[v]<<", MCMD.LocalXmin="<<MCMD.LocalXmin[v]<<", xmin="<<xmin[v]<<"\n";
    cout <<"old_xmax="<<old_xmax[v]<<", SimPM.Xmax="<<SimPM.Xmax[v]<<", MCMD.LocalXmax="<<MCMD.LocalXmax[v]<<", xmax="<<xmax[v]<<"\n";
  }
  cout.flush();

  //
  // get a setup_grid class, and use it to set up the grid!
  //
  class setup_fixed_grid *SimSetup =0;
  SimSetup = new setup_fixed_grid_pllel();
  class GridBaseClass *grid = 0;
  err  = MCMD.decomposeDomain();
  if (err) rep.error("main: failed to decompose domain!",err);
  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  SimSetup->setup_grid(&grid,&MCMD);
  if (!grid) rep.error("Grid setup failed",grid);

  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;



  //
  // Now setup microphysics and raytracing classes
  //
  err += SimSetup->setup_microphysics();
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics and raytracing",err);

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";

  size_t ifile=0;

  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";
  clk.start_timer("analyse_data");

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile=0; ifile<nfiles; ifile += Nskip) {
    cout <<"------ Starting Next Loop: ifile="<<ifile<<", time so far=";
    cout <<clk.time_so_far("analyse_data")<<" ----\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Reading Simulation data to grid -------\n";
    cout.flush();

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
    for (size_t skip=0; skip<Nskip; skip++) ff++;

    //
    // delete any current radiation sources
    //
    reset_radiation_sources(&(SimPM.RS));

    //
    // Read header to get timestep info.
    //
    err = dataio.ReadHeader(infile);
    if (err) rep.error("Didn't read header",err);

    // *****************************************************
    // Delete radiation sources again
    // *****************************************************
    reset_radiation_sources(&(SimPM.RS));
    //cout.flush();
    //
    // Now reset domain.
    //
    reset_domain(xmin,xmax,SimPM.ndim,MCMD);
    //
    // Read data (this reader can read serial or parallel data.
    //
    err = dataio.parallel_read_any_data(infile, grid);
    rep.errorTest("(main) Failed to read data",0,err);
    
#ifdef TESTING
    for (size_t v=0; v<static_cast<size_t>(SimPM.ndim); v++) {
      cout <<"old_xmin="<<old_xmin[v]<<", SimPM.Xmin="<<SimPM.Xmin[v]<<", MCMD.LocalXmin="<<MCMD.LocalXmin[v]<<", xmin="<<xmin[v]<<", grid->SIM_Xmin="<<grid->SIM_iXmin(static_cast<axes>(v))<<"\n";
      cout <<"old_xmax="<<old_xmax[v]<<", SimPM.Xmax="<<SimPM.Xmax[v]<<", MCMD.LocalXmax="<<MCMD.LocalXmax[v]<<", xmax="<<xmax[v]<<", grid->SIM_Xmax="<<grid->SIM_iXmax(static_cast<axes>(v))<<"\n";
    }
    cout.flush();
#endif // TESTING

    //
    // If domain is bigger than original, then add more data.
    //
    cell *c=grid->FirstPt();
    for (size_t v=0; v<static_cast<size_t>(SimPM.ndim); v++) {
      double vals[SimPM.nvar];
      //
      // First check for extensions in the negative direction:
      //
      if ((old_xmin[v]-xmin[v]) < -0.5*grid->DX()) {
        cout <<"dimension "<<v<<" has enarged domain in -ve dir.\n";
        rep.printVec("old",old_xmin,SimPM.ndim);
        rep.printVec("new",xmin,SimPM.ndim);
        cout <<(old_xmin[v]-xmin[v]) << "  "<<-0.5*grid->DX()<<"\n";
        //
        // Navigate to first point on the original grid.
        //
        for (size_t d=0; d<static_cast<size_t>(SimPM.ndim); d++) {
          enum direction negdd=static_cast<direction>(2*d);
          enum direction posdd=grid->OppDir(negdd);
          while (CI.get_dpos(c,d) < old_xmin[d])
            c=grid->NextPt(c,posdd);
        }
        //
        // Copy the state vector
        //
        for (size_t d=0; d<static_cast<size_t>(SimPM.nvar); d++) {
          vals[d] = c->P[d];
        }
        //
        // Write all new cells with this state vector.
        //
        c=grid->FirstPt();
        do {
          if (CI.get_dpos(c,v) < old_xmin[v]) {
            for (size_t d=0; d<static_cast<size_t>(SimPM.nvar); d++) {
              c->P[d] = vals[d];
            }
          }
        } while ((c=grid->NextPt(c))!=0);
      } // if extension in negative direction
      //
      // Then extensions in the positive direction:
      //
      if ((xmax[v]-old_xmax[v]) > 0.5*grid->DX()) {
        cout <<"dimension "<<v<<" has enarged domain in +ve dir.\n";
        //
        // Navigate to last point on the old grid in this dir
        //
        c=grid->LastPt();
        for (size_t d=0; d<static_cast<size_t>(SimPM.ndim); d++) {
          enum direction negdd=static_cast<direction>(2*d);
          enum direction posdd=grid->OppDir(negdd);
          while ( CI.get_dpos(c,d) > old_xmax[d]) {
            c=grid->NextPt(c,negdd);
          }
        }
        //
        // Copy the state vector
        //
        for (size_t d=0; d<static_cast<size_t>(SimPM.nvar); d++) {
          vals[d] = c->P[d];
        }
        //
        // Write all new cells with this state vector.
        //
        c=grid->FirstPt();
        do {
          if (CI.get_dpos(c,v) > old_xmax[v]) {
            for (size_t d=0; d<static_cast<size_t>(SimPM.nvar); d++) {
              c->P[d] = vals[d];
            }
          }
        } while ((c=grid->NextPt(c))!=0);
      } // if extension in positive direction 
    }  // loop over directions.

    cout <<"--------------- Finished Reading Data  ----------------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Starting Writing Data  ----------------\n";

    temp.str("");
    temp <<op_path<<"/"<<outfile;
    writer.OutputData(temp.str(), grid, SimPM.timestep);
    cout <<"--------------- Finished Writing Data  ----------------\n";

  } // Loop over all files.    


  cout <<"-------------------------------------------------------\n";
  cout <<"---- Finised with all Files: time=";
  cout <<clk.stop_timer("analyse_data")<<"--------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";

  if(grid!=0) {
    cout << "\t Deleting Grid Data..." << endl;
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


