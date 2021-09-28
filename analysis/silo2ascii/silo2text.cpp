/// \file silo2text.cpp
/// \author Jonathan Mackey
/// 
/// This file reads in silo files generated by serial and parallel code
/// and re-writes the data as a single text file. 
///
/// Mods:
/// - 2011.12.14 JM: converted from silocompare.cc
/// - 2013.09.05 JM: Fixed for pion; added microphysics so ascii
///    files can have gas temperature written out.
/// - 2013.10.04 JM: Added output frequency to skip silo files.
/// - 2014.06.07 JM: Added mpv8 -- heating/cooling microphysics class
/// - 2016.03.03 JM: renamed to silo2text.cpp and tried to get it working.

#ifndef PARALLEL
# error "define PARALLEL so this will work!"
#endif

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"

#include "decomposition/MCMD_control.h"
#include "grid/setup_fixed_grid_MPI.h"
#include "grid/uniform_grid.h"

#include <dirent.h>
#include <errno.h>
#include <list>
#include <iostream>
#include <sstream>
using namespace std;
#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#include <silo.h>

#include "microphysics/microphysics_base.h"

#ifdef CODE_EXT_HHE
#include "future/mpv9_HHe.h"
#endif

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
  //
  // get a setup_grid class, to set up the grid, and a grid pointer.
  //
  class setup_fixed_grid *SimSetup =0;
  SimSetup = new setup_fixed_grid_pllel();
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  cout <<"Projection3D: myrank="<<myrank<<", nproc="<<nproc<<"\n";
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);


  //
  // Get an input file and an output file.
  //
  if (argc!=5) {
    cerr << "Error: must call as follows...\n";
    cerr << "silo2text: <silo2text> <source-dir> <file-base>  <output-file> <op-freq>\n";
    cerr << "  op-freq: if this is e.g. 10, we only convert every";
    cerr << " 10th input file to a text file.\n";
    rep.error("Bad number of args",argc);
  }
  string fdir = argv[1];
  string firstfile= argv[2];
  string outfilebase =argv[3]; string outfile;
  size_t op_freq = atoi(argv[4]);

  string rts("msg_"); rts += outfilebase;
  //rep.redirect(rts);

  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio(SimPM,"DOUBLE",&(SimPM.levels[0].MCMD));
  //
  // Get list of first and second files to read, and make sure they match.
  //
  list<string> files;
  err += dataio.get_files_in_dir(fdir, firstfile,  &files);
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
  // Set up iterators to run through all the files.
  //
  list<string>::iterator ff=files.begin();
  unsigned int nfiles = files.size();
  vector<class GridBaseClass *> G;
  class GridBaseClass *grid = 0;
  
  //
  // loop over all files: open first and write a text file.
  //
  for (unsigned int fff=0; fff<nfiles; fff += op_freq) {
    ostringstream oo;
    oo.str(""); oo<<fdir<<"/"<<*ff;
    firstfile = oo.str(); 

    class file_status fstat;
    if (!fstat.file_exists(firstfile)) {
       cout <<"first file: "<<firstfile<<"\n";
       rep.error("First file doesn't exist",firstfile);
    }
    
    //
    // Read in first code header so i know how to setup grid.
    //
    err = dataio.ReadHeader(firstfile,SimPM);
    if (err) rep.error("Didn't read header",err);
    SimPM.grid_nlevels = 1;
    SimPM.levels[0].parent=0;
    SimPM.levels[0].child=0;
    SimPM.levels[0].Ncell = SimPM.Ncell;
    for (int v=0;v<MAX_DIM;v++) SimPM.levels[0].NG[v] = SimPM.NG[v];
    for (int v=0;v<MAX_DIM;v++) SimPM.levels[0].Range[v] = SimPM.Range[v];
    for (int v=0;v<MAX_DIM;v++) SimPM.levels[0].Xmin[v] = SimPM.Xmin[v];
    for (int v=0;v<MAX_DIM;v++) SimPM.levels[0].Xmax[v] = SimPM.Xmax[v];
    SimPM.levels[0].dx = SimPM.Range[XX]/SimPM.NG[XX];
    SimPM.levels[0].simtime = SimPM.simtime;
    SimPM.levels[0].dt = 0.0;
    SimPM.levels[0].multiplier = 1;
    err  = SimPM.levels[0].MCMD.decomposeDomain(SimPM.ndim,SimPM.levels[0]);
    if (err) rep.error("main: failed to decompose domain!",err);
  
    // *****************************************************
    //
    // Setup microphysics and grid
    //
    if (fff==0) {
      //
      // Now setup microphysics and raytracing classes
      //
      err += SimSetup->setup_microphysics(SimPM);
      //err += SimSetup->setup_raytracing(SimPM);
      err += SimSetup->set_equations(SimPM);
      if (err) rep.error("Setup of microphysics and raytracing",err);

      //
      // May need to setup extra data in each cell for ray-tracing optical
      // depths and/or viscosity variables (here just set it to zero).
      //
      CI.setup_extra_data(SimPM.RS,0,0,0);
      //
      // Setup grid.
      //
      G.resize(1);
      SimSetup->setup_grid(G, SimPM);
      grid = G[0];
      if (!grid) rep.error("Grid setup failed",grid);
      cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<"\n";
    } // first step, setup grid and microphysics.
    // *****************************************************

    
    //
    // Set output file name.
    //
    oo.str("");
    oo <<outfilebase<<".";
    oo.fill('0'); oo.width(8);
    oo<<SimPM.timestep<<".txt";
    outfile = oo.str();
    cout <<"\n**************************************************************************************\n";
    cout <<"fff="<<fff<<"\tinput file: "<<firstfile<<"\toutput file: "<<outfile<<"\n";
    
    //
    // Read data (this reader can read serial or parallel data.
    //
    err = dataio.ReadData(firstfile, G, SimPM);
    rep.errorTest("(main) Failed to read data",0,err);
    // ***************************************************************
    // ********* FINISHED FIRST FILE, MOVE ON TO OUTPUT FILE *********
    // ***************************************************************
    
    class dataio_text textio (SimPM);
    err += textio.OutputData(outfilebase, G, SimPM, SimPM.timestep);
    if (err) {cerr<<"\t Error writing data.\n"; return(1);}

    //
    // move onto next first and second files
    //
    for (size_t vv=0;vv<op_freq;vv++) ff++;
  } // move onto next file

  //
  // Finish up and quit.
  //
  delete COMM; COMM=0;
  return 0;
}

