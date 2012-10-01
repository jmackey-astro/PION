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


#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
#include <cmath>
using namespace std;
#include "../../../uniform_grid_code/trunk/source/global.h"
#include "../../../uniform_grid_code/trunk/source/dataio_utility.h"
#include "../../../uniform_grid_code/trunk/source/uniformGrid.h"
#include "../../../uniform_grid_code/trunk/source/dataio_fits.h"


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

  //*******************************************************************
  //*******************************************************************
  //
  // Get input files and an output file.
  //
  if (argc<7) {
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
  //ostringstream redir; redir.str(""); redir<<outfile<<"_msg_";
  //rep.redirect(redir.str());


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
  if (grid) rep.error("grid already setup, so bugging out",grid);
  try {
    grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
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
  string this_outfile("");
  this_outfile=outfile;
  unsigned int ifile=0;

  //
  // Text File
  //
  ofstream outf;
  if (outf.is_open())
    rep.error("Output text file is open! Close file before opening a new one!",1);
  outf.open(this_outfile.c_str());
  if (!outf.is_open()) rep.error("Failed to open text file for writing",this_outfile);
  outf <<"#\n# writing: simtime timestep WHATEVER OTHER STUFF I WANT TO WRITE PER LINE!!!\n#\n\n";

  
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";

  //*******************************************************************
  // loop over all files:
  //*******************************************************************

  for (ifile=0; ifile<nfiles; ifile++) {
    cout <<"--------------- Starting Next Loop: ifile="<<ifile<<"------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Reading Simulation data to grid -------\n";

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
    err = dataio.serial_read_any_data(infile, ///< file to read from
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
    double mass=0.0;
    double cell_vol = grid->DV();

    cell *c=grid->FirstPt();
    do {
      mass += c->P[RO];
    } while ((c=grid->NextPt(c))!=0);

    mass *= cell_vol;

    cout <<"--------------- Finished Analysing this step ----------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Writing image and getting next Im-file \n";

    //**********************
    //* Write Data to file *
    //**********************
    //
    // In theory this can write to multiple output files, but I
    // haven't implemented that yet.  Also if it's not a text file or
    // a single file for all timesteps, code will have bugged out by
    // now.
    //
    outf <<SimPM.simtime<<"\t"<<SimPM.timestep<<"\t"<<mass<<endl;


  } // Loop over all files.    

  //
  // Close file, if single file.
  //
  if (!multi_opfiles) {
    outf.close();
  }

  cout <<"--------------- Finised Analysing all Files -----------\n";
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


