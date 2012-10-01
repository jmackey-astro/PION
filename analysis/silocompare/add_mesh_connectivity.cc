///
/// file:    add_mesh_connectivity.cc
/// author:  Jonathan Mackey
/// date:    2010-01-28
///
/// Description: This file reads in a silo file written with an arbitrary number
///  of MPI processes, and re-writes it with the number of processes this
///  program is run on.
///
///  It's original purpose was to add multi-mesh-adjacency info for VisIt to 
///  enable correct plotting of streamlines, but that is only needed now for
///  data generated with pre-2010 simulations.
///
///
/// Modifications:
///
///  - 2010-01-28 JM: Calculate a bunch of things for the Lim &
///   Mellema (2003) simulations.
///
///  - 2010-02-03 JM: Made it work for multi-core use, so can run with
///  8, 16, 32 cores etc. (and the number of cores doesn't have to
///  match the number of cores used to write the data i.e. independent
///  of the number of quadmeshes in the silo file).
///
///  - 2010-02-08 JM: Modified analyse_data.cc to make
///    add_mesh_connectivity.cc so I can plot streamlines in my MHD
///    sims.
///
///  - 2012.01.04 JM: Updated to work with new grid/source code.
///


#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
#include <cmath>
using namespace std;
#include "../../source/global.h"
#include "../../source/dataio_utility.h"
#include "../../source/uniformGrid.h"
#include "../../source/dataio_fits.h"


#define OP_TEXT 0
#define OP_FITS 1
#define OP_SILO 2



int main(int argc, char **argv)
{
  //
  // First initialise the comms class, since we need to define
  // PARALLEL to read a parallel file.  This is not ideal, but too
  // bad...
  //
  int err = COMM->init(&argc, &argv);

  double runtime = 10000.0; // about 3 hours.
  mpiPM.set_max_walltime(runtime);

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc!=4) {
    cout << "Use as follows:\n";
    cout << "executable-filename: <executable-filename> <input-path> <input-silo-file-base>\n";
    cout << "\t\t <output-file> \n";
    cout <<"******************************************\n";
    cout <<"input path:   path to input files.\n";
    cout <<"input file:   base filename of sequence of filesn.\n";
    cout <<"output file:  filename for output file(s) (with path).\n";
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
  string outfile    = argv[3];

  //
  // Redirect output to a text file if you want to:
  //
  //  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_"<<mpiPM.myrank<<"_";
  rep.redirect(redir.str());


  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Getting List of Files to read ---------\n";
  
  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio;

  //
  // Get list of files to read:
  //
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file,  &files);
  if (err) rep.error("failed to get list of files",err);
  for (list<string>::iterator s=files.begin(); s!=files.end(); s++)
  cout <<"files: "<<*s<<endl;
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
  if ( (err=mpiPM.decomposeDomain()) !=0) 
    rep.error("Couldn't Decompose Domain!",err);

  if (grid) rep.error("grid already setup, so bugging out",grid);
  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables (here just set it to zero).
  //
  CI.setup_extra_data(SimPM.RS,0,0);
  try {
    grid = new UniformGridParallel (SimPM.ndim, SimPM.nvar, SimPM.eqntype, mpiPM.LocalXmin, mpiPM.LocalXmax, mpiPM.LocalNG);
  }
  catch (std::bad_alloc) {
    rep.error("(trunks::setup_grid) Couldn't assign data!", grid);
  }
  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Output file: if multiple files, we will append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to
  // label output file.
  //
  //string filehandle("");
  // string this_outfile;
  // if (op_filetype==OP_TEXT) {
  //   this_outfile=outfile+".txt";
  // }
  // else {
  //   rep.error("bad op_filetype",op_filetype);
  // }
  unsigned int ifile=0;



  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";
  GS.start_timer("analyse_data");
  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile=0; ifile<nfiles; ifile++) {
    cout <<"------ Starting Next Loop: ifile="<<ifile<<", time so far="<<GS.time_so_far("analyse_data")<<" ----\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Reading Simulation data to grid -------\n";

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
    ff++;

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

    //********************
    //* Analyse the Data *
    //********************
    //

    //
    // TESTING: RESCALE BOX SIZE
    //
    // for (int v=0;v<SimPM.ndim; v++) {
    //   SimPM.Xmin[v] = -10.0;
    //   SimPM.Xmax[v] =  10.0;
    //   SimPM.Range[v] = 20.0;
    // }
    // SimPM.dx = SimPM.Range[XX]/SimPM.NG[XX];
    // mpiPM.decomposeDomain();
    //
    // TESTING: RESCALE BOX SIZE
    //

    //
    // Replace pressure with myrank, for testing the boundary smoothing.
    //
    // cell *c=grid->FirstPt();
    // do {
    //   c->P[PG] = static_cast<double>(mpiPM.myrank);
    // } while ( (c=grid->NextPt(c)) !=0);

    //cout <<"--------------- Writing image and getting next Im-file \n";

    // rep.printVec("Xmin  ",SimPM.Xmin,SimPM.ndim);
    // rep.printVec("Xmax  ",SimPM.Xmax,SimPM.ndim);
    // rep.printVec("Range ",SimPM.Range,SimPM.ndim);
    // rep.printVec("LXmin ",mpiPM.LocalXmin,SimPM.ndim);
    // rep.printVec("LXmax ",mpiPM.LocalXmax,SimPM.ndim);
    // rep.printVec("LRange",mpiPM.LocalRange,SimPM.ndim);

    //
    // TESTING TESTING
    //class dataio_silo_utility dataio2;
    //ostringstream testfile; testfile<<"./"<<outfile<<"_np"<<mpiPM.nproc;
    if (mpiPM.nproc>1)
      dataio.dataio_silo_pllel::OutputData(outfile, grid, SimPM.timestep);
    else 
      dataio.dataio_silo::OutputData(outfile, grid, SimPM.timestep);
    // TESTING TESTING
    //

  } // Loop over all files.    
  

  cout <<"-------------------------------------------------------\n";
  cout <<"---- Finised with all Files: time="<<GS.stop_timer("analyse_data")<<"--------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";

  if(grid!=0) {
    cout << "\t Deleting Grid Data..." << endl;
    delete grid; grid=0;
  }

  COMM->finalise();
  delete COMM; COMM=0;

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------


