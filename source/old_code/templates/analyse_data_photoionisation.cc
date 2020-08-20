///
/// file:    analyse_data.cc
/// author:  Jonathan Mackey
/// date:    2010-01-28
///
/// Description: This is a template file for analysing simulation data
/// using a single processor, but reading data which may have been
/// written from multiple processors.  In its current form it can only
/// write text data to a single file.
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


#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
#include <cmath>
using namespace std;
#include "../../../source/global.h"
#include "../../../source/dataio_utility.h"
#include "../../../source/uniformGrid.h"
#include "../../../source/dataio_fits.h"


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
  if (argc<6) {
    cout << "Use as follows:\n";
    cout << "executable-filename: <executable-filename> <input-path> <input-silo-file-base>\n";
    cout << "\t\t <output-file> <op-file-type> <multi-opfiles> \n";
    cout << "\t\t <ANY-EXTRA-STUFF> \n";
    cout <<"******************************************\n";
    cout <<"input path:   path to input files.\n";
    cout <<"input file:   base filename of sequence of filesn.\n";
    cout <<"output file:  filename for output file(s).\n";
    cout <<"op-file-type: integer/string [0,text,TEXT], [1,fits,FITS], [2,silo,SILO] (ONLY TEXT WORKING OUT OF THE BOX!!!)\n";
    cout <<"muti-opfiles: integer. 0=only one output file. 1=one output file per step. (ONLY SINGLE FILE WORKING!!!)\n";
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


  //
  // What sort of output will depend on what sort of analysis we are
  // doing, usually a single text file will be sufficient for global
  // properties of the simulation.
  //
  string optype=argv[4];
  int op_filetype;
  if      (optype=="0" || optype=="text" || optype=="TEXT") {
    op_filetype = OP_TEXT;
    cout <<"\t\toutputting data to text file.\n";
  }
  else if (optype=="1" || optype=="fits" || optype=="FITS") {
    op_filetype = OP_FITS;
    cout <<"\t\toutputting data to fits files.\n";
    rep.error("don't know how to output fits files yet... fix me please (see ET2009 projection code)!","sorry");
  }
  else if (optype=="2" || optype=="silo" || optype=="SILO") {
    op_filetype = OP_SILO;
    cout <<"\t\toutputting data to silo files.\n";
    rep.error("don't know how to output silo files yet... fix me please!","sorry");
  }
  else 
    rep.error("What sort of output is this?",optype);

  //
  // If we write a single output file for all steps or not.
  //
  int multi_opfiles = atoi(argv[5]);
  switch (multi_opfiles) {
  case 0:
    cout <<"\t\tOutputting all timesteps in a single file.\n";
    break;
  case 1:
    cout <<"\t\tOutputting timesteps in different files.\n";
    rep.error("haven't implemented this yet!!!",multi_opfiles);
    break;
  default:
    rep.error("Bad multi-files value",multi_opfiles);
  }


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
  string this_outfile;
  if (op_filetype==OP_TEXT) {
    this_outfile=outfile+".txt";
  }
  else {
    rep.error("bad op_filetype",op_filetype);
  }
  unsigned int ifile=0;


  //
  // Text File: proc 0 writes all the data.
  // 
  ofstream outf;
  if (mpiPM.myrank==0) {
    if (outf.is_open())
      rep.error("Output text file is open! Close file before opening a new one!",1);
    outf.open(this_outfile.c_str());
    if (!outf.is_open()) rep.error("Failed to open text file for writing",this_outfile);
    outf <<"#\n# writing: simtime timestep tot_mass n_mass cl_n_mass";
    outf <<"  vel_vwt nvel_vwt ivel_vwt  vel_mwt nvel_mwt ivel_mwt ";
    outf <<"  leading_edge <n_pos> <cl_n_pos> \n";
    outf <<"#  [time: 1 2]  [mass: 3 4 5] [vel_vol: 6 7 8] [vel_mass: 9 10 11] [posn: 12 13 14]\n";
    outf <<"#\n\n";
  }

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
    // Get tracer indices in state vector.
    //
    string tracers = SimPM.trtype.substr(6);
    // cout <<"tracers: "<<tracers<<endl;
    int len = (SimPM.trtype.length() +5)/6 -1; // should give number of tracers.
    int tr_Hp     = -1;     ///< tracer for ion fraction.
    int tr_colour = -1; ///< tracer for colour (clump material=1; ambient=0).
    string s;
    for (int v=0;v<len;v++) {
      s = tracers.substr(6*(v),6); // Get 'i'th tracer variable.
      if (s=="H1+___" || s=="HII__") {
	tr_Hp = SimPM.ftr +v;
      }
      if (s=="colour" || s=="Colour") {
	tr_colour = SimPM.ftr +v;
      }
    }
    if (tr_Hp<0 || tr_colour<0)
      rep.error("Need tracer values for colour and ion fraction!",tr_Hp);

    //
    // Define variables we want to track
    //
    double cell_vol = grid->DV();
    double nmass_norm=0.0, imass_norm=0.0, nvol_norm=0.0, ivol_norm=0.0, clnmass_norm=0.0;
    double nx_cut=0.9, ix_cut=0.05;
    //
    // mass measures
    //
    double total_mass   = 0.0;     ///< total gas mass.
    double neutral_mass = 0.0;     ///< total neutral gas mass.
    double clump_nmass = 0.0;      ///< total neutral clump material.
    //
    // velocities
    //
    double vel_volwt = 0.0;   ///< mean gas velocity (total, volume weighting).
    double vel_masswt = 0.0;  ///< mean gas velocity (total, mass weighting).
    double nvel_volwt = 0.0;   ///< mean gas velocity (neutral, volume weighting).
    double nvel_masswt = 0.0;  ///< mean gas velocity (neutral, mass weighting).
    double ivel_volwt = 0.0;   ///< mean gas velocity (ionised, volume weighting).
    double ivel_masswt = 0.0;  ///< mean gas velocity (ionised, mass weighting).
    //
    // x-positions.
    //
    double leading_edge = 1.e200; ///< position of first neutral cell (x<0.01).
    double npos       = 0.0; ///< mean neutral gas location (mass weighting).
    double clump_npos = 0.0; ///< mean neutral clump gas location (mass weighting).
    double xpos=0.0;

    cell *c=grid->FirstPt();
    do {
      xpos = CI.get_dpos(c,XX)-SimPM.Xmin[XX];

      //
      // mass
      //
      total_mass   += c->P[RO];
      neutral_mass += c->P[RO]*(1.0-c->P[tr_Hp]);
      clump_nmass  += c->P[RO]*(1.0-c->P[tr_Hp])*c->P[tr_colour];

      //
      // velocity
      //
      vel_volwt  += c->P[VX];
      vel_masswt += c->P[VX]*c->P[RO];
      //
      // position
      //
      if (c->P[tr_Hp]<ix_cut)
	leading_edge = std::min(leading_edge, xpos);

      //
      // Try this so that only un-ionised cells count in the neutrals
      // sum 
      //
      if (c->P[tr_Hp]<nx_cut) {
	nmass_norm += (1.0-c->P[tr_Hp])*c->P[RO];
	clnmass_norm += (1.0-c->P[tr_Hp])*c->P[RO]*c->P[tr_colour];
	nvol_norm  += 1.0;

	nvel_volwt  += c->P[VX]*(1.0-c->P[tr_Hp]);
	nvel_masswt += c->P[VX]*(1.0-c->P[tr_Hp])*c->P[RO];

	npos       += xpos*(1.0-c->P[tr_Hp])*c->P[RO];
	clump_npos += xpos*(1.0-c->P[tr_Hp])*c->P[RO]*c->P[tr_colour];
      }
      //
      // Try this so that only ionised cells count in the ions sum
      //
      if (c->P[tr_Hp]>ix_cut) {
	imass_norm += c->P[tr_Hp]*c->P[RO];
	ivol_norm  += 1.0;

	ivel_volwt  += c->P[VX]*c->P[tr_Hp];
	ivel_masswt += c->P[VX]*c->P[tr_Hp]*c->P[RO];
      }

    } while ((c=grid->NextPt(c))!=0);

    //
    // Now do global communication so that we get the global
    // max/min/mean and normalisations.
    //
    total_mass   = COMM->global_operation_double("SUM", total_mass);
    neutral_mass = COMM->global_operation_double("SUM", neutral_mass);
    clump_nmass  = COMM->global_operation_double("SUM", clump_nmass);
    vel_volwt    = COMM->global_operation_double("SUM", vel_volwt);
    vel_masswt   = COMM->global_operation_double("SUM", vel_masswt);
    nmass_norm   = COMM->global_operation_double("SUM", nmass_norm);
    clnmass_norm = COMM->global_operation_double("SUM", clnmass_norm);
    nvol_norm    = COMM->global_operation_double("SUM", nvol_norm);
    nvel_volwt   = COMM->global_operation_double("SUM", nvel_volwt);
    nvel_masswt  = COMM->global_operation_double("SUM", nvel_masswt);
    npos         = COMM->global_operation_double("SUM", npos);
    clump_npos   = COMM->global_operation_double("SUM", clump_npos);
    imass_norm   = COMM->global_operation_double("SUM", imass_norm);
    ivol_norm    = COMM->global_operation_double("SUM", ivol_norm);
    ivel_volwt   = COMM->global_operation_double("SUM", ivel_volwt);
    ivel_masswt  = COMM->global_operation_double("SUM", ivel_masswt);
    leading_edge = COMM->global_operation_double("MIN", leading_edge);

    vel_volwt /= SimPM.Ncell;
    nvel_volwt /= nvol_norm;
    ivel_volwt /= ivol_norm;

    //nvel_volwt /= SimPM.Ncell;
    //ivel_volwt /= SimPM.Ncell;
    //vel_masswt /= total_mass;
    //nvel_masswt /= neutral_mass;
    //ivel_masswt /= (total_mass-neutral_mass);
    //npos /= neutral_mass;
    //clump_npos /= clump_nmass;

    vel_masswt /= total_mass;
    nvel_masswt /= nmass_norm;
    ivel_masswt /= imass_norm;
    npos /= nmass_norm;
    clump_npos /= clnmass_norm;

    //
    // Turn mass variables from densities to masses
    //
    total_mass   *= cell_vol;
    neutral_mass *= cell_vol;
    clump_nmass  *= cell_vol;

    //
    // check position:
    //
    if (leading_edge>1.e100)
      leading_edge = SimPM.Range[XX];
    if (isnan(npos) || isinf(npos) || npos>1.e100)
      npos = SimPM.Range[XX];
    if (isnan(clump_npos) || isinf(clump_npos) || clump_npos>1.e100)
      clump_npos = SimPM.Range[XX];


    //cout <<"--------------- Finished Analysing this step ----------\n";
    //cout <<"-------------------------------------------------------\n";
    //cout <<"--------------- Writing image and getting next Im-file \n";

    if (mpiPM.myrank==0) {
      // **********************
      // * Write Data to file *
      // **********************
      //
      // In theory this can write to multiple output files, but I
      // haven't implemented that yet.  Also if it's not a text file or
      // a single file for all timesteps, code will have bugged out by
      // now.
      //
      outf <<SimPM.simtime<<"  "<<SimPM.timestep<<"  ";
      outf <<total_mass<<"  "<< neutral_mass <<"  "<< clump_nmass <<"  ";
      outf << vel_volwt  <<"  "<< nvel_volwt  <<"  "<< ivel_volwt  <<"  ";
      outf << vel_masswt <<"  "<< nvel_masswt <<"  "<< ivel_masswt <<"  ";
      outf << leading_edge <<"  "<< npos <<"  "<< clump_npos <<"  \n";
      outf.flush(); // make sure it gets to file in case we bug out!!!
      //
      // Close output file, if multiple files.
      //
      if (multi_opfiles) {
	outf.close();
      }

    } // if proc 0, write data
    //
    // TESTING TESTING
    //ostringstream testfile; testfile<<"tmp_data/tmpfile_np"<<mpiPM.nproc;
    //if (mpiPM.nproc==1)
    ///dataio.dataio_silo::OutputData(testfile.str(), grid);
    //else
    //dataio.dataio_silo_pllel::OutputData(testfile.str(), grid);
    // TESTING TESTING
    //

  } // Loop over all files.    
  

  //
  // Close file, if single file.
  //
  if (!multi_opfiles && mpiPM.myrank==0) {
    outf.close();
  }

  cout <<"-------------------------------------------------------\n";
  cout <<"---- Finised Analysing all Files: time="<<GS.stop_timer("analyse_data")<<"--------\n";
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


