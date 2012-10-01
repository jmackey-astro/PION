///
/// \file blastwave_multifile.cc
///
/// \author Jonathan Mackey
///
/// This code takes a input file-base, reads in a sequence of silo
/// files matching this base-filename, calculates the mean shock
/// radius as a function of time through the simulation, and also
/// outputs the radial profile of density, pressure and radial
/// velocity for the second (first with t>0) and last (t=t_f)
/// simulation outputs.
///

#include <iostream>
#include <sstream>
#include <silo.h>
#include <cmath>
using namespace std;
#include "../../source/global.h"
#include "../../source/dataio_utility.h"
#include "../../source/uniformGrid.h"


int main(int argc, char **argv)
{
  //
  // First initialise MPI via the comms class.
  //
  int err = COMM->init(&argc, &argv);

  //
  // Get input directory and file-base, output directory and two
  // output file-names from cmd-line args.
  //
  if (argc<6 || argc>7) {
    cout <<"Error: must call with at least 6 arguments...\n";
    cout <<"blastwave_multifile.cc: <executable> <input-dir> ";
    cout <<"<infile-base> <output-dir> <r(t)-outfile> ";
    cout <<"<P(R)-outfile-base> [redirect-stdout?]\n";
    cout <<"e.g.: ./bw_multifile ./ sim20_0000 ./output sim20_rt.txt sim20_profile msg\n";
    rep.error("Bad number of Args",argc);
  }
  string input_dir   = argv[1];
  string infile_base = argv[2];
  string output_dir  = argv[3];
  string outfile_rad = argv[4];
  string outfile_pro = argv[5];

  //
  // redirect stdout if requested
  //
  if (argc==7) {
    ostringstream redir; redir.str("");
    redir<<output_dir<<"/"<<argv[6];
    rep.redirect(redir.str());
  }
    

  cout <<"reading from file "<<infile_base<<" in dir: "<<input_dir<<"\n";
  cout <<"Writing shock radius vs. time to file "<<outfile_rad;
  cout <<" and radial profiles to files "<<outfile_pro<<" in dir: ";
  cout <<output_dir<<"\n";
  cout <<"**********************************************\n";

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
  err += dataio.get_files_in_dir(input_dir, infile_base,  &files);
  if (err) rep.error("failed to get list of files",err);
  for (list<string>::iterator s=files.begin(); s!=files.end(); s++)
  cout <<"files: "<<*s<<endl;
  int nfiles = static_cast<int>(files.size());
  if (nfiles<1) rep.error("Need at least one file, but got none",nfiles);

  cout <<"--------------- Got list of Files ---------------------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"------------- Reading first file ----------------------\n";
  cout <<"--------------- Setting up Grid -----------------------\n";

  //
  // Set up an iterator to run through all the files.
  //
  list<string>::iterator ff=files.begin();
  //
  // Open first file, read header, and setup grid
  //
  ostringstream temp; temp <<input_dir<<"/"<<*ff;
  string first_file = temp.str();
  temp.str("");
  err = dataio.ReadHeader(first_file);
  if (err) rep.error("Didn't read header",err);

#ifndef HCORR
#error "Need H-correction variables for this code!"
#endif
  //
  // Don't want too much extra data in cells, to minimise memory.
  //
  int rt_flag = 0;
  int hc_flag = 2; // set distance and velocity for each cell!
  CI.setup_extra_data(rt_flag, hc_flag);
  //
  // Set low-memory cells
  //
  CI.set_minimal_cell_data();

  //
  // Now we can setup the grid:
  //
  cout <<"(UniformFV::setup_grid) Setting up grid...\n";
  if (grid) rep.error("Grid already set up!",grid);

#ifdef GEOMETRIC_GRID
  if      (SimPM.coord_sys==COORD_CRT)
    grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_CYL)
    grid = new uniform_grid_cyl (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_SPH)
    grid = new uniform_grid_sph (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else 
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
#else  // GEOMETRIC_GRID
  grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype, SimPM.Xmin, SimPM.Xmax, SimPM.NG);
#endif // GEOMETRIC_GRID

  if (grid==0) rep.error("(setup_grid) Couldn't assign data!", grid);
  cout <<"(setup_grid) Done. g="<<grid;
  cout <<"\tDX = "<<grid->DX()<<"\n";
  cout <<"------------------------------------------------------\n\n";

  cout <<"--------------- Finished Setting up Grid --------------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Opening time series outfile -----------\n";

  string outfile = output_dir;
  outfile += '/';
  outfile += outfile_rad;
  if (dataio.file_exists(outfile)) {
    cout <<"WARNING:: file exists, I am overwriting a text file...\n";
  }
  cout <<"output file: "<<outfile<<"\n";
  ofstream outf(outfile.c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",outfile);
  cout <<"writing to file "<<outfile<<endl;
  outf.setf( ios_base::scientific );
  outf.precision(6);
  outf <<"# BW output.  file: "<<infile_base<<endl;
  outf <<"# Columns are time step shock_pos shock_inneredge shock_outerdege\n";
  outf <<"#\n#\n";

  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting Loop over all input files ----\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Assume we have the problem set up so that the centre of symmetry
  // is the origin:
  //
  double BW_origin [SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++)
    BW_origin[i]=0.0;
  
 

  //*******************************************************************
  // loop over all files:
  //*******************************************************************
  int ifile=0;
  for (ifile=0; ifile<nfiles; ifile++) {
    cout <<"--------------- Starting Next Loop: ifile="<<ifile<<"------\n";
    cout <<"-------------------------------------------------------\n";
    cout <<"--------------- Reading Simulation data to grid -------\n";

    string infile(input_dir);
    if (infile[infile.size()-1] != '/')
      infile += "/";
    infile += *ff; // this points to the current element in the list of input files.
    cout <<"---- file: "<<infile<<" ----\n";

    ff++; // increment file pointer for next loop.

    //
    // Read header to get timestep info; 
    //
    err = dataio.ReadHeader(infile);
    if (err) rep.error("Didn't read header",err);
    cout <<"############ SIMULATION TIME: "<<SimPM.simtime/3.16e7;
    cout <<" yrs for step="<<ifile<<"   ############\n";

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
    // We want to calculate the mean, min, max. radius of the shock
    // for all values.
    //
    // If we are on the second or last timestep, then also output N
    // values of position, velocity, density, pressure.
    //

    //
    // First establish criteria for where the shock is:
    // (1) If pressure is >2x the external pressure,
    // (2) If pressure is <0.7 of the central value,
    // then we are at the outside edge of the shock for an adiabatic
    // blastwave.
    //
    cell *c = grid->FirstPt();
    //
    // get to cell nearest to origin:
    //
    for (int id=0;id<SimPM.ndim;id++) {
      while (c!=0 && CI.get_dpos(c,id) <0.0) {
	c = grid->NextPt(c, static_cast<direction>(2*id+1));
      }
      if (!c)
	rep.error("Origin is not on the simulation domain!!!","HELP!");
    }
    double p_in = c->P[PG];
    c= grid->LastPt();
    double p_out = c->P[PG];
    double pratio = p_in/p_out;
    if (pratio<50.0)
      rep.error("low pressure ratio",pratio);
    
    double pmin, pmax;
    pmin = p_out*2.0;
    pmax = p_in*0.7;
    cout <<"P_in="<<p_in<<" P_out="<<p_out<<": pmin/max="<<pmin<<", "<<pmax<<"\n";

    double dist;
    vector<double> distance;
    c = grid->FirstPt();
    do {

      //
      // Get distance from origin to cell.
      //
#ifdef    GEOMETRIC_GRID
      //
      // Assume origin is at a cell vertex.
      //
      dist = grid->distance_vertex2cell(BW_origin,c);
#else  // GEOMETRIC_GRID
      //
      // This returns the geometric cell--centre
      //
      double d[SimPM.ndim];
      CI.get_dpos(c,d);
      dist = GS.distance(BW_origin,d,SimPM.ndim);
#endif // GEOMETRIC_GRID
      //
      // Set the XX value of H-correction data to be distance
      //
      CI.set_Hcorr(c,XX,dist); 

      if (c->P[PG]>pmin && c->P[PG]<pmax) {
	distance.push_back(dist);
	cout <<"d="<<dist<<" and p="<<c->P[PG]<<endl;
      }
    } while ((c=grid->NextPt(c)) !=0);

    //
    // Now calculate mean, min, max, shock-pos for this timestep, and
    // then write to the output_rad file.
    //
    double sum=0.0, lv=1.e99, uv=0.0;
    for (unsigned int i=0; i<distance.size(); i++) {
      sum += distance[i];
      lv = std::min(lv,distance[i]);
      uv = std::max(uv,distance[i]);
    }
    sum /= distance.size();
    cout <<"****** "<<SimPM.simtime<<"\t"<<sum<<"\t"<<lv<<"\t"<<uv<<"\n";
    cout <<"\n***************************************************\n";
    // outf <<"# Columns are time step shock_pos shock_inneredge shock_outerdege\n";
    outf << SimPM.simtime <<"\t"<< SimPM.timestep;
    outf <<"\t"<< sum <<"\t"<< lv <<"\t"<< uv <<"\n";

    //
    // If we are on the 2nd or last output, then write a radial profile:
    //
    if (ifile==1 || ifile==(nfiles-1)) {
      cout <<"-------------------------------------------------------\n";
      cout <<"------------- WRITING PROFILE FILE i="<<ifile<<" -----------------\n";
      //
      // open output file:
      //
      string outfile2(output_dir);
      if (outfile2[outfile2.size()-1] != '/')
	outfile2 += "/";
      outfile2 += outfile_pro;
      if (ifile==1) 
	outfile2 += "_early.txt";
      else 
	outfile2 += "_late.txt";

      if (dataio.file_exists(outfile2)) {
	cout <<"WARNING:: file exists, I am overwriting a text file...\n";
      }
      cout <<"output file: "<<outfile2<<"\n";
      ofstream outf2(outfile2.c_str());
      if(!outf2.is_open()) rep.error("couldn't open outfile2",outfile2);
      cout <<"writing to file "<<outfile2<<endl;
      outf2.setf( ios_base::scientific );
      outf2.precision(6);
      
      outf2 <<"# file: "<<outfile2<<"\n";
      outf2 <<"# input file: "<<infile<<"\n";
      outf2 <<"# simtime = "<<SimPM.simtime<<", Step = "<<SimPM.timestep<<"\n";
      outf2 <<"#\n# radius  density pressure velocity \n\n";

      long int Nsamples = std::min(SimPM.Ncell, static_cast<long int>(10000));
      int isample =static_cast<int>(SimPM.Ncell/Nsamples);
      long int ct=0;
      c = grid->FirstPt();
      do {
	if (ct%isample==0) {
	  //
	  // distance (already calculated)
	  //
	  outf2 << CI.get_Hcorr(c,XX);
	  outf2 <<"  ";
	  outf2 << c->P[RO];
	  outf2 <<"  ";
	  outf2 << c->P[PG];
	  outf2 <<"  ";
	  outf2 << sqrt(c->P[VX]*c->P[VX] +c->P[VY]*c->P[VY] +c->P[VZ]*c->P[VZ]);
	  outf2 <<"\n";
	}
	ct ++;
      } while ((c=grid->NextPt(c)) !=0);
      outf2.close();
      cout <<"------------- PROFILE FILE WRITTEN --------------------\n";
      cout <<"-------------------------------------------------------\n";
    } // if writing profile...

    distance.clear();
      
  } // Loop over all files.    
	
  outf.close();
      
  cout <<"--------------- Finised Analysing all Files -----------\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";
   
  delete grid; grid=0;

  COMM->finalise();
  delete COMM; COMM=0;
  return 0;
}

