///
/// file:    make_text_files.cpp
/// author:  Jonathan Mackey
/// date:    2013.03.10
///
/// Adapted from templates/HIIregion_KineticE.cpp
///
/// Description: This is a template file for analysing simulation 
/// data with N cores, but reading data which may have been written
/// from M!=N processors.  In its current form it can only write text
/// data to a single file.
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
///
/// - 2013.03.10 JM: Modified for the StarBench Spitzer test problem.
/// - 2013.06.27 JM: Modified for the StarBench Spitzer test problem.


#include <iostream>
#include <sstream>
#include <silo.h>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <global.h>
#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"

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


#ifdef NEW_METALLICITY
#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"
#endif // NEW_METALLICITY




#define OP_TEXT 0
#define OP_FITS 1
#define OP_SILO 2

// ##################################################################
// ##################################################################
///
/// set filename based on counter, outfile-base-name.
///
std::string set_filename(
      const std::string outfile,
      const long int    counter
      )
{
  ostringstream temp; temp.str("");
  temp << outfile <<".";
  //
  // only add timestep if counter >=0 (initial conditions get no timestep).
  //
  if (counter >= 0) {
    temp.width(8); temp.fill('0');
    temp << counter <<".";
  }
  temp <<"txt";
  return temp.str();
}


///
/// Reset the radiation sources in the header to correspond to projected
/// quantities and not the sources used for the simulation.
///
void reset_radiation_sources(struct rad_sources *);

///
/// Function to setup the microphysics class (just so I can calculate gas
/// temperature in a manner consistent with how it was done in the 
/// simulation).  Function copied from gridMethods.cc
///
int setup_microphysics();

// ##################################################################
// ##################################################################




int main(int argc, char **argv)
{
  //
  // First initialise the comms class
  //
  int err = COMM->init(&argc, &argv);

  if (mpiPM.nproc>1)
    rep.error("This is serial code so far",mpiPM.nproc);

  //MP=0; RT=0; grid=0;

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc!=5) {
    cout << "Use as follows:\n";
    cout << "<executable-filename> <input-path>";
    cout << " <input-silo-file-base>";
    cout << " <output-path> <output-file> \n";
    cout <<"******************************************\n";
    cout <<"input path:   path to input files.\n";
    cout <<"input file:   base filename of sequence of filesn.\n";
    cout <<"output path:  directory to write output files to.\n";
    cout <<"output file:  filename for output file(s).\n";
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
  // Redirect output to a text file if you want to:
  //
  //  ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  ostringstream redir; redir.str("");
  redir<<op_path<<"/"<<outfile<<"_msg_"<<mpiPM.myrank<<"_";
  rep.redirect(redir.str());

  int op_filetype = OP_TEXT;
  cout <<"\t\toutputting data to text file.\n";

  double runtime = 100000.0; // about 28 hours.
  mpiPM.set_max_walltime(runtime);

  //
  // If we write a single output file for all steps or not.
  //
  int multi_opfiles = 0;
  cout <<"\t\tOutputting all timesteps in a single file.\n";


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
  // First decompose the domain, so I know the dimensions of the local
  // grid to set up.  If nproc==1, then this sets the local domain to
  // be the full domain.
  //
  if ( (err=mpiPM.decomposeDomain()) !=0) 
    rep.error("Couldn't Decompose Domain!",err);

  // *****************************************************
  // Now delete all radiation "sources" from SimPM.RS, to avoid allocating
  // memory for column densities in the cell-data.
  // *****************************************************
  reset_radiation_sources(&(SimPM.RS));
  //
  // Setup extra data in each cell for ray-tracing optical depths and/or
  // multi-D viscosity variables (here just set it to zero).
  //
  CI.setup_extra_data(SimPM.RS,0,0);

  if (grid) rep.error("grid already setup, so bugging out",grid);

  if      (SimPM.coord_sys==COORD_CRT) {
    grid = new UniformGridParallel (SimPM.ndim, SimPM.nvar,
				    SimPM.eqntype,  mpiPM.LocalXmin,
				    mpiPM.LocalXmax, mpiPM.LocalNG);
  }
  else if (SimPM.coord_sys==COORD_CYL) {
    grid = new uniform_grid_cyl_parallel (SimPM.ndim, SimPM.nvar,
					  SimPM.eqntype,  mpiPM.LocalXmin,
					  mpiPM.LocalXmax, mpiPM.LocalNG);
  }
  else if (SimPM.coord_sys==COORD_SPH) {
    grid = new uniform_grid_sph_parallel (SimPM.ndim, SimPM.nvar,
					  SimPM.eqntype,  mpiPM.LocalXmin,
					  mpiPM.LocalXmax, mpiPM.LocalNG);
  }
  else {
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
  }

  if (!grid)
    rep.error("(setup_grid) Couldn't assign data!", grid);
  
  cout <<"\t\tg="<<grid<<"\tDX = "<<grid->DX()<<endl;


  //
  // Now setup microphysics and raytracing classes
  //
  err += setup_microphysics();
  //err += setup_raytracing();
  if (err) rep.error("Setup of microphysics and raytracing",err);

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Output file: if multiple files, we will append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to
  // label output file.
  //
  //string filehandle("");
  string this_outfile;
  this_outfile = op_path + "/" + outfile;
  size_t ifile=0;


  double pc  = 3.086e18;
  double Myr = 3.16e13;


  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";
  GS.start_timer("analyse_data");

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile=0; ifile<101; ifile+=5) {
    cout.flush();
    cout <<"------ Starting Next Loop: ifile="<<ifile<<", time so far=";
    cout <<GS.time_so_far("analyse_data")<<" ----\n";
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
    std::advance(ff,5);

    //
    // delete any current radiation sources
    //
    if (SimPM.RS.Nsources!=0) {
      SimPM.RS.sources.clear();
      SimPM.RS.Nsources=0;
    }

    //
    // Read header to get timestep info.
    //
    err = dataio.ReadHeader(infile);
    if (err) rep.error("Didn't read header",err);

    // *****************************************************
    // Delete radiation sources again
    // *****************************************************
    //reset_radiation_sources(&(SimPM.RS));
    //cout.flush();

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


    //cout <<"--------------- Finished Analysing this step ----------\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Writing data and getting next src-file \n";
    
    string tfile = set_filename(this_outfile, SimPM.timestep);
    ofstream outf(tfile.c_str());
    if(!outf.is_open()) {
      cerr << "Error opening file " << tfile << " for writing.  Quitting..." <<"\n";
      return(1);
    }
    
    outf << "# format: r density pressure temperature\n";
    outf << "# time = "<<SimPM.simtime/Myr<<" Myr.  timestep = "<<SimPM.timestep<<"\n";
    outf.setf( ios_base::scientific );
    outf.precision(6);

    // Go through every point, output one line per point.
    class cell *cpt=grid->FirstPt(); do {
      // First positions.
      outf << CI.get_dpos(cpt,0)/pc << "  ";
      // Next density pressure temperature
      outf <<cpt->P[RO]<<"  ";
      outf <<cpt->P[PG]<<"  ";
      if (MP) outf << MP->Temperature(cpt->P,SimPM.gamma);
      outf  <<"\n";
    } while ( (cpt=grid->NextPt(cpt))!=0);
    
    outf.setf(ios_base::fmtflags(0),ios_base::floatfield);
    outf.precision(6);
    outf.close();

  } // Loop over all files.    


  cout <<"-------------------------------------------------------\n";
  cout <<"---- Finised with all Files: time="<<GS.stop_timer("analyse_data")<<"--------\n";
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




// stolen from gridMethods.cc
int setup_microphysics()
{
  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  //
  // Setup Microphysics class, if needed.
  // First see if we want the only_cooling class (much simpler), and if
  // not then check for the one of the bigger microphysics classes.
  //
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout <<"\t******* Requested cooling but no chemistry... setting";
    cout <<" up mp_only_cooling() class, with timestep-limiting.\n";
    SimPM.EP.MP_timestep_limit = 1;
    MP = new mp_only_cooling(SimPM.nvar, &(SimPM.EP));
    if (!MP) rep.error("mp_only_cooling() init",MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    cout <<"TRTYPE: "<<SimPM.trtype<<"\n";
    string mptype;
    if (SimPM.trtype.size() >=6)
      mptype = SimPM.trtype.substr(0,6); // Get first 6 chars for type of MP.
    else mptype = "None";
    bool have_set_MP=false;


#ifndef EXCLUDE_MPV1
    if      (mptype=="ChAH__" || mptype=="onlyH_") {
      cout <<"\t******* setting up MP_Hydrogen microphysics module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MP_Hydrogen(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      cout <<"\t**---** WARNING, THIS MODULE HAS BEEN SUPERSEDED BY mp_implicit_H. **--**\n";
      have_set_MP=true;
    }
#endif // exclude MPv1


#ifndef EXCLUDE_HD_MODULE
    if (mptype=="lowZ__") {
      cout <<"\t******* setting up microphysics_lowz module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new microphysics_lowz(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }
#endif // exclude Harpreet's module

#ifndef EXCLUDE_MPV2
    if (mptype=="MPv2__") {
#ifdef MP_V2_AIFA
      cout <<"\t******* setting up mp_v2_aifa module *********\n";
      cout <<"\t******* N.B. Timestep limiting is enforced. **\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_v2_aifa(SimPM.nvar, SimPM.ntracer, SimPM.trtype);
      SimPM.EP.MP_timestep_limit = 1;
#else
      rep.error("Enable mp_v2_aifa as an ifdef if you really want to use it",2);
#endif
      have_set_MP=true;
    }
#endif // exclude MPv2


#ifndef EXCLUDE_MPV3
    if (mptype=="MPv3__") {
      cout <<"\t******* setting up mp_explicit_H module *********\n";
#if MPV3_DTLIMIT>=0 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. Timestep limiting is enforced by #def";
      cout <<" MPV3_DTLIMIT="<<MPV3_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif

      MP = new mp_explicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype
#ifdef NEW_METALLICITY
      , &(SimPM.EP)
#endif // NEW_METALLICITY
      );
      //if (SimPM.EP.MP_timestep_limit != 1)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv3


#ifndef EXCLUDE_MPV4
    if (mptype=="MPv4__") {
      cout <<"\t******* setting up mp_implicit_H module *********\n";
#if MPV4_DTLIMIT>=5 && MPV4_DTLIMIT<=12
      cout <<"\t******* N.B. dt05-12 Timestep limiting is enforced by #def";
      cout <<" DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =5;
#elif MPV4_DTLIMIT>=0 && MPV4_DTLIMIT<=4
      cout <<"\t******* N.B. dt00-04 Timestep limiting is enforced by #def";
      cout <<" MPV4_DTLIMIT="<<MPV4_DTLIMIT<<". **\n";
      SimPM.EP.MP_timestep_limit =4;
#else
#error "No timestep-limiting is defined in source/defines/functionality_flags.h"
#endif
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mp_implicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype);
      //SimPM.EP.MP_timestep_limit = 4;  // limit by recombination time only
      //if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv4


#ifdef NEW_METALLICITY
    if (mptype=="MPv5__") {
      cout <<"\t******* setting up mpv5_molecular module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv5_molecular(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv6__") {
      cout <<"\t******* setting up mpv6_PureH module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv6_PureH(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }

    if (mptype=="MPv7__") {
      cout <<"\t******* setting up mpv7_TwoTempIso module *********\n";
      SimPM.EP.MP_timestep_limit = 1;
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new mpv7_TwoTempIso(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      have_set_MP=true;
    }
#endif // NEW_METALLICITY




#ifndef EXCLUDE_MPV1
    //
    // Finally, if MP has not been set up yet, try to set up the v0
    // microphysics integrator, which is slow, but can model a number
    // of elements and ions.
    //
    if (!have_set_MP) {
      cout <<"\t******* setting up MicroPhysics (v0) module *********\n";
      if (have_set_MP) rep.error("MP already initialised",mptype);
      MP = new MicroPhysics(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
        rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv1/0

    if (!MP) rep.error("microphysics init",MP);
    if (!have_set_MP) rep.error("HUH? have_set_MP",have_set_MP);
  }
  else {
    cout <<"\t******** not doing microphysics.\n";
    MP=0;
  }

  cout <<"************************************************************\n";
  cout <<"***************** MICROPHYSICS SETUP ***********************\n";
  cout <<"************************************************************\n";
  return 0;
}


// ##################################################################
// ##################################################################




