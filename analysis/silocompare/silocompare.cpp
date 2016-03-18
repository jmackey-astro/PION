/// \file silocompare.cpp
/// \author Jonathan Mackey
/// 
/// This file reads in silo files generated by serial and parallel code
/// and calculates the differences between them if they have the same
/// domains.  This resulting output is outputted as it's own silo file.
///
/// Mods:
///  - 2009-12-25 JM: Added L1 and L2 error calculations.
///
///  - 2010-09-03 JM: Updated support for new dataio class interface.
///
/// - 2010.12.07 JM: Added setup_extra_data() call to cell_interface
///   before setting up grid.  Need to add Geomtric grid options!
/// - 2011.04.23 JM: new code for setup_extra_data().
/// - 2015.03.26 JM: updated for pion v0.2
/// - 2015.06.18 JM: updated to allow FLOAT or DOUBLE silo data.
/// - 2016.03.16 JM: renamed to silocompare.cpp
/// - 2016.03.18 JM: updated to work better with grid/double/float

#ifndef PARALLEL
# error "define PARALLEL so this will work!"
#endif

#include <list>
#include <iostream>
#include <sstream>
#include <silo.h>
#include <cmath>
using namespace std;

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "constants.h"

#include "sim_params.h"
#include "dataIO/dataio.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#include "grid/uniform_grid.h"

#include "MCMD_control.h"
#include "setup_fixed_grid_MPI.h"


// ##################################################################
// ##################################################################




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

  //
  // Get an input file and an output file.
  //
  if (argc!=9) {
    cerr << "Error: must call as follows...\n";
    cerr << "silocompare: <silocompare> <first-dir> <comp-dir> <first-file> <flt/dbl> <comp-file> <flt/dbl> <outfile> <fabs/plus-minus/L1/L2>\n";
    cerr << "\t 0: Diff image is relative error for rho/p_g (abs.val)\n";
    cerr << "\t 1: Diff image is relative error for rho/p_g (+/-)\n";
    cerr << "\t 2: Just calculate L1+L2 error, no difference image.\n";
    cerr << "<flt/dbl>: FLOAT or DOUBLE for whether you think the file is float or double precision.\n";
    rep.error("Bad number of args",argc);
  }
  string fdir = argv[1];
  string sdir = argv[2];
  string firstfile   =argv[3];
  string dtype1(argv[4]);
  string secondfile  =argv[5];
  string dtype2(argv[6]);
  string outfilebase =argv[7]; string outfile;
  int optype=atoi(argv[8]);
  if (optype<0 || optype>2)
    rep.error("Please set optype to 0 (abs val.) or 1 (+- val.) or 2 (L1/L2)",optype);
  cout <<"fdir="<<fdir<<"\tsdir="<<sdir<<endl;
  cout <<"first file: "<<firstfile<<"\tsecond file: "<<secondfile<<"\toutput file: "<<outfilebase<<"\n";

  //string rts("msg_"); rts += outfilebase;
  //rep.redirect(rts);

  //
  // set up dataio_utility class
  //
  class dataio_silo_utility dataio(dtype2,&MCMD);

  // ----------------------------------------------------------------
  // ----------------------------------------------------------------

  //
  // Get list of files to read:
  //
  list<string> ffiles, sfiles;
  err += dataio.get_files_in_dir(fdir, firstfile,  &ffiles);
  if (err) rep.error("failed to get first list of files",err);
  err += dataio.get_files_in_dir(sdir, secondfile, &sfiles);;
  if (err) rep.error("failed to get second list of files",err);

  for (list<string>::iterator s=ffiles.begin(); s!=ffiles.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo")==string::npos) {
      cout <<"removing file "<<*s<<" from list.\n";
      ffiles.erase(s);
    }
    else {
      cout <<"files: "<<*s<<endl;
    }
  }
  for (list<string>::iterator s=sfiles.begin(); s!=sfiles.end(); s++) {
    // If file is not a .silo file, then remove it from the list.
    if ((*s).find(".silo")==string::npos) {
      cout <<"removing file "<<*s<<" from list.\n";
      sfiles.erase(s);
    }
    else {
      cout <<"files: "<<*s<<endl;
    }
  }

  //
  // Set up iterators to run through all the files.
  //
  list<string>::iterator ff=ffiles.begin();
  list<string>::iterator ss=sfiles.begin();

  //
  // Set the number of files to be the length of the shortest list.
  //
  unsigned int nfiles = ffiles.size();
  nfiles = std::min(nfiles, static_cast<unsigned int>(sfiles.size()));
  if (nfiles<1) rep.error("Need at least one file, but got none",nfiles);
  
  // ----------------------------------------------------------------
  // ----------------------------------------------------------------

  //
  // Read the first file, and setup the grid based on its parameters.
  //
  ostringstream oo;
  oo.str(""); oo<<fdir<<"/"<<*ff; firstfile =oo.str();
  err = dataio.ReadHeader(firstfile);
  if (err) rep.error("Didn't read header",err);

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

  // ----------------------------------------------------------------
  // ----------------------------------------------------------------

  //
  // loop over all files: open first+second and write the difference.
  //
  for (unsigned int fff=0; fff<nfiles; fff++) {
    
    oo.str(""); oo<<fdir<<"/"<<*ff; firstfile  =oo.str(); 
    oo.str(""); oo<<sdir<<"/"<<*ss; secondfile =oo.str(); 
    oo.str(""); oo <<outfilebase<<"."; oo.fill('0'); oo.width(5); oo<<fff<<".silo";
    outfile = oo.str();
    cout <<"\n**************************************************************************************\n";
    cout <<"fff="<<fff<<"\tfirst file: "<<firstfile;
    cout <<"\tsecond file: "<<secondfile<<"\toutput file: "<<outfile<<"\n";

    class file_status fstat;
    if (!fstat.file_exists(firstfile) ||
        !fstat.file_exists(secondfile)) {
       cout <<"first file: "<<firstfile;
       cout <<"\tand second file: "<<secondfile<<endl;
       rep.error("First or second file doesn't exist",secondfile);
    }
    
    //
    // Setup output silo I/O classes.
    //
    class dataio_silo firstio("DOUBLE");
    
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------

    //
    // Read in first code header so i know how to setup grid.
    //
    err = dataio.ReadHeader(firstfile);
    if (err) rep.error("Didn't read header",err);
  
    // ----------------------------------------------------------------
    // ----------------------------------------------------------------

    //
    // Read data (this reader can read serial or parallel data).
    //
    err = dataio.parallel_read_any_data(firstfile, grid);
    rep.errorTest("(silocompare) Failed to read firstfile",0,err);

    // ----------------------------------------------------------------
    // ----------------------------------------------------------------
    
    //
    // Now put the data in c->Ph for first file:
    //
    cell *c = grid->FirstPt();
    do {
      for (int v=0; v<SimPM.nvar; v++)
	c->Ph[v] = c->P[v];
    } while ( (c=grid->NextPt(c))!=0);

    cout <<"FINISHED reading first file: "<<firstfile<<endl;

    // ***************************************************************
    // ********* FINISHED FIRST FILE, MOVE ON TO SECOND FILE *********
    // ***************************************************************
    
    //
    // Read data (this reader can read serial or parallel data).
    //
    err = dataio.parallel_read_any_data(secondfile, grid);
    rep.errorTest("(silocompare) Failed to read secondfile",0,err);

    cout <<"FINISHED reading second file: "<<secondfile<<endl;

    //c = grid->FirstPt();
    //do {
    //  if (c->Ph[RO] < 2.0e-30) {
    //    CI.print_cell(c);
    //  }
    //} while ( (c=grid->NextPt(c))!=0);
    c = grid->FirstPt();
    rep.printVec("P ",c->P ,SimPM.nvar);
    rep.printVec("Ph",c->Ph,SimPM.nvar);
    

    // *********************************************************************
    // ********* FINISHED SECOND FILE, MOVE ON TO WRITE COMPARISON *********
    // *********************************************************************
    //
    // write difference to output file.
    //
    c=grid->FirstPt(); 
    double maxdiff[SimPM.nvar]; double temp;
    double reldiff[SimPM.nvar];
    double absdiff[SimPM.nvar];
    for (int v=0;v<SimPM.nvar;v++) maxdiff[v] = 0.0;
    for (int v=0;v<SimPM.nvar;v++) reldiff[v] = 0.0;
    for (int v=0;v<SimPM.nvar;v++) absdiff[v] = 0.0;
    int ipos[SimPM.ndim];

    rep.printVec("P ",c->P ,SimPM.nvar);
    rep.printVec("Ph",c->Ph,SimPM.nvar);
    
    switch (optype) {
    case 0: 
      //
      // fabs(error)
      //
      do {
	CI.get_ipos(c,ipos);
	for (int v=0;v<SimPM.nvar;v++) {
	  temp = std::max(maxdiff[v], fabs(c->P[v]-c->Ph[v]));
	  if (temp>maxdiff[v]) {
	    //cout <<"new max, cell posn="; rep.printVec("pos",ipos,SimPM.ndim);
	  }
	  maxdiff[v] = temp;
	  //
	  // Relative difference:
	  //
	  temp = fabs(c->P[v]-c->Ph[v])/(fabs(c->P[v])+fabs(c->Ph[v]) +SMALLVALUE*SimPM.RefVec[v]);
	  reldiff[v] = std::max(reldiff[v],temp);
	  //
	  // Create difference image:
	  //
	  c->P[v] -= c->Ph[v];
	  if (v==RO) c->P[v] = temp; // relative differences in rho
	  if (v==PG) c->P[v] = temp; // relative differences in p_g
	  c->P[v] = fabs(c->P[v]);
	}
      } while ( (c=grid->NextPt(c))!=0);
      break;
    case 1:
      //
      // not abs(diff) in the image
      //
      do {
	CI.get_ipos(c,ipos);
	for (int v=0;v<SimPM.nvar;v++) {
	  temp = std::max(maxdiff[v], fabs(c->P[v]-c->Ph[v]));
	  if (temp>maxdiff[v]) {
	    //cout <<"new max, cell posn="; rep.printVec("pos",ipos,SimPM.ndim);
	  }
	  maxdiff[v] = temp;
	  //
	  // Relative difference:
	  //
	  temp = fabs(c->P[v]-c->Ph[v])/(fabs(c->P[v])+fabs(c->Ph[v]) +SMALLVALUE*SimPM.RefVec[v]);
	  reldiff[v] = std::max(reldiff[v],temp);
	  //
	  // Create difference image:
	  //
	  c->P[v] -= c->Ph[v];
	  if (v==RO) c->P[v] = temp; // relative differences in rho
	  if (v==PG) c->P[v] = temp; // relative differences in p_g
	}
      } while ( (c=grid->NextPt(c))!=0);
      break;
    case 2:
      //
      // get the L1 and L2 error. L1=absdiff; L2=reldiff;
      //
      for (int v=0;v<SimPM.nvar;v++) {
        absdiff[v]=0.0;
	reldiff[v]=0.0;
	maxdiff[v]=0.0;
      }
      do {
	//
	// L1 = (1/n)sum(fabs(p1[v]-p0[v]))
	// L2 = sqrt[(1/n)sum(fabs(p1[v]-p0[v])^2)]
	//
	for (int v=0;v<SimPM.nvar;v++) {
	  temp = fabs(c->P[v]-c->Ph[v]);
	  absdiff[v] += temp;
	  reldiff[v] += temp*temp;
	  maxdiff[v] = std::max(maxdiff[v],temp);
	}
      } while ( (c=grid->NextPt(c))!=0);
      //
      // Now divide by N and take sqrt for L2
      //
      for (int v=0;v<SimPM.nvar;v++) {
	absdiff[v] /= SimPM.Ncell;
	reldiff[v] /= SimPM.Ncell;
	reldiff[v] = sqrt(reldiff[v]);
      }
      break;
    default:
      rep.error("Input a valid optype!!! (should have caught this!)",optype);
      break;
    }


    //
    // Output differences to screen/file.
    // Write difference data to file.
    //
    ofstream outf;
    switch (optype) {
    case 0: case 1:
      cout <<"Max. diffs for each variable follow: i, fabs(diff)\n";
      for (int v=0;v<SimPM.nvar;v++) cout <<v<<"\t"<<maxdiff[v]<<endl;
      cout <<"\n\n";
      cout <<"Rel. diffs for each variable follow: i, fabs(diff)/fabs(sum)\n";
      for (int v=0;v<SimPM.nvar;v++) cout <<v<<"\t"<<reldiff[v]<<endl;
      cout <<"\n\n";
      firstio.OutputData(outfile, grid, fff);
      break;
    case 2:
      oo.str(""); oo <<outfilebase<<"_"<<r<<".txt"; outfile=oo.str(); oo.str("");
      outf.open(outfile.c_str(), fstream::app);
      cout <<"Writing data to Output File: "<<outfile<<"\n";
      outf <<"# file1="<<firstfile<<"\n# file2="<<secondfile<<endl;
      cout <<"\n\t***: var L1err (Mean-Absolute Error) L2err ";
      cout <<"max-err-abs ****\n";
      outf <<"\n***: var L1err L2err ";
      outf <<"(root-mean-squared error) max-err-abs ***\n";
      outf <<" L1err = (Mean-Absolute Error), ";
      outf <<" L2err = (root-mean-squared error).\n";
      for (int v=0;v<SimPM.nvar;v++) {
	cout <<v<<"\t"<<absdiff[v]<<"\t"<<reldiff[v]<<"\t"<<absdiff[v]<<endl;
	outf <<v<<"    "<<absdiff[v]<<"    "<<reldiff[v]<<"    "<<absdiff[v]<<endl;
      }
      cout <<"\t\t***********************\n";
      outf.close();
      break;
    default:
      rep.error("Input a valid optype!!! (should have caught this!)",optype);
      break;
    }

    //
    // move onto next first and second files
    //
    ff++;
    ss++;
  } // move onto next file

  //
  // Finish up and quit.
  //
  if (grid) {delete grid; grid=0;}
  COMM->finalise();
  delete COMM; COMM=0;
  //MPI_Finalize();
  return 0;
}



