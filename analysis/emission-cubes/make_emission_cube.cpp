///
/// file:    make_emission_cube.cpp
/// author:  Jonathan Mackey
/// date:    2020.07.16
///
/// Description:
/// - Generate FITS data cubes of various point emission quantities
///   in units of erg/cm^3/s (sometimes also per unit energy).
///
/// Modifications:
///
/// - 2020.07.16 JM: wrote code


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#include "constants.h"
#include "sim_params.h"

#include "decomposition/MCMD_control.h"
//#include "grid/setup_fixed_grid_MPI.h"
//#include "grid/uniform_grid.h"
#include "grid/setup_grid_NG_MPI.h"

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"

#include "microphysics/microphysics_base.h"
#include "microphysics/mp_only_cooling.h"
#ifndef EXCLUDE_MPV3
#include "microphysics/MPv3.h"
#endif
#ifndef EXCLUDE_MPV4
#include "microphysics/MPv4.h"
#endif 
#include "microphysics/MPv5.h"
#include "microphysics/MPv6.h"
#include "microphysics/MPv7.h"


#include "raytracing/raytracer_SC.h"
#include "../xray/xray_emission.h"

#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
using namespace std;


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
  int myrank=-1, nproc=-1;
  COMM->get_rank_nproc(&myrank,&nproc);
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(myrank);
  SimPM.levels[0].MCMD.set_nproc(nproc);
  if (nproc>1)
    rep.error("This is serial code",nproc);

  //
  // Get input files and an output file.
  //
  if (argc<5) {
    cout << "Use as follows:\n";
    cout << "./emission-cube <input-path>";
    cout << " <input-silo-file-base>";
    cout << " <output-path> <output-fits-file-base> <var>";
    cout << "\n";
    cout <<"******************************************\n";
    cout <<"input path:  path to input files.\n";
    cout <<"input file:  base filename of sequence of files including _0000 if parallel.\n";
    cout <<"output path: directory to write output files to.\n";
    cout <<"output file: filename for output FITS file(s).\n";
    cout <<"var:         what to plot (currently just Halpha)\n";
    rep.error("Bad number of args",argc);
  }

  string input_path = argv[1];
  string input_file = argv[2];
  string op_path    = argv[3];
  string outfile    = argv[4];

  string var = argv[5];
  cout <<"Ignoring var="<<var<<" directive while debugging code.\n";

  //
  // Redirect output to a text file if you want to:
  //
  ostringstream redir; redir.str("");
  redir<<op_path<<"/msg_"<<outfile<<"_rank"<<myrank<<"_";
  //rep.redirect(redir.str());

  //*******************************************************************
  // Get input files, read header, setup grid
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Getting List of Files to read ---------\n";
  class dataio_silo_utility dataio (SimPM, "DOUBLE", &(SimPM.levels[0].MCMD));
  list<string> files;
  err += dataio.get_files_in_dir(input_path, input_file,  &files);
  if (err) rep.error("failed to get list of files",err);
  for (list<string>::iterator s=files.begin(); s!=files.end(); s++) {
    if ((*s).find(".silo")==string::npos) {
      cout <<"removing file "<<*s<<" from list.\n";
      files.erase(s);
      s=files.begin();
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

  CI.set_minimal_cell_data();
  list<string>::iterator ff=files.begin();
  ostringstream temp; temp <<input_path<<"/"<<*ff;
  string first_file = temp.str();
  temp.str("");
  err = dataio.ReadHeader(first_file, SimPM);
  if (err) rep.error("Didn't read header",err);


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

  cout.flush();

  class setup_grid_NG_MPI *SimSetup =0;
  SimSetup = new setup_grid_NG_MPI();
  SimSetup->setup_NG_grid_levels(SimPM);
  vector<class GridBaseClass *> G;
  G.resize(SimPM.grid_nlevels);
  SimSetup->setup_grid(G, SimPM);
  class GridBaseClass *grid = G[0];
  if (!grid) rep.error("Grid setup failed",grid);
  SimPM.dx = grid->DX();
  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;

  err += SimSetup->setup_microphysics(SimPM);
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics",err);

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";

  size_t ifile=0;
  class DataIOFits_pllel writer (SimPM, &(SimPM.levels[0].MCMD));
  class Xray_emission xray ();

  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";
  clk.start_timer("analyse_data");

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile=0; ifile<nfiles; ifile++) {
    ifile += Nskip;
    cout <<"------ Starting Next Loop: ifile="<<ifile<<", time so far=";
    cout <<clk.time_so_far("analyse_data")<<" ----\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Reading Simulation data to grid -------\n";
    cout <<" reading file "<<*ff<<"\n";
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
    ff++;
    for (size_t skip=0; skip<Nskip; skip++) ff++;

    // Read header to get timestep info.
    err = dataio.ReadHeader(infile, SimPM);
    if (err) rep.error("Didn't read header",err);
    //cout.flush();

    // read data onto grid.
    err = dataio.ReadData(infile, G, SimPM);
    rep.errorTest("(main) Failed to read data",0,err);
    
    cout <<"--------------- Finished Reading Data  ----------------\n";
    cout <<"-------------------------------------------------------\n";

    cout <<"--------------- Calculating Emission   ----------------\n";
    double T=0.0, em=0.0, ne=0.0, np=0.0;
    for (int l=0; l<SimPM.grid_nlevels; l++) {
      grid = G[l];
      cell *c = grid->FirstPt();
      do {
        if (MP) {
          T  = MP->Temperature(c->Ph,SimPM.gamma);
          ne = MP->get_n_elec(c->Ph);
          np = MP->get_n_Hplus(c->Ph);
        }
        else {
          T = ne = np = 0.0;
          rep.error("need microphysics to get temperature",MP);
        }
        em = xray.Halpha_emissivity(T)*ne*np;
        c->P[PG] = em;
      } while ((c=grid->NextPt(c))!=0);
    }

    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Starting Writing Data  ----------------\n";

    temp.str("");
    temp <<op_path<<"/"<<outfile;
    writer.OutputData(temp.str(), G, SimPM, SimPM.timestep);
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

  COMM->finalise();
  delete COMM; COMM=0;

  return 0;
}


// ##################################################################
// ##################################################################


