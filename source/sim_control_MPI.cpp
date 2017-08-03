/// \file sim_control_MPI.cpp
/// 
/// \brief Parallel Grid Methods Class Member Function definitions.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for 
/// the "sim_control_fixedgrid_pllel" class, which is a modification
/// of the basic 1st/2nd order Finite Volume Solver according to the
/// method outlined in Falle, Komissarov, \& Joarder (1998),MNRAS,297,265.
/// 
/// Modifications:
/// - 2007-10-11 Started writing file.
/// - 2007-11-11 Basically working with uniform grid and all sorts of boundaries.
/// - 2010-01-05 JM: Put RT step limiting in an #ifdef so that it is
///    not set for test problems.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits.  This also
///    meant I don't need a new parallel output_data() function so I
///    got rid of it.  All dataio classes now choose the filename
///   themselves.
/// - 2010.07.23 JM: New RSP source position class interface.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///    extra_data.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.15 JM: Added CI.setup_extra_data() call in setup_grid()
///     function.
/// - 2011.01.12 JM: Set so only proc 0 displays RT timestep.
///     Commented out some stuff in calc_timestep which was used for
///     testing and made the code less efficient.
/// - 2011.02.17 JM: Raytracer header file moved to raytracing/ subdir.
///     More ray-tracing options for CI.setup_extra_data
/// - 2011.02.25 JM: removed HCORR ifdef around new code. Modified setup_extra_data
///     for ray-tracing again.  Now we have N column-density vars for N sources.
/// - 2011.03.01 JM: fixed timings bug in time_int()
/// - 2011.03.02 JM: Added parallelised raytracer-shielding class setup.
/// - 2011.03.21 JM: moved cell setup_extra_data() to its own function, to save 
///     on copy-paste for the parallel version.
///     Deleted commented out output_data() function, which is now the same as the
///     serial version.
/// - 2011.03.22 JM: Updated setup_raytracing() for improved raytracer functionality.
/// - 2011.10.24 JM: Updated setup_raytracing() again.  It's better now.
/// - 2012.01.16 JM: Added update_evolving_RT_source() function to time_int().
/// - 2012.05.14 JM: Added check in setup_raytracing() for sources at infinity
///    in different directions (must set up the parallelised raytracer then).
/// - 2012.05.16 JM: fixed bug from last change.
/// - 2012.08.06 JM: Added separators between functions for clarity.
/// - 2013.04.16 JM: Fixed FITS read functions for new filename convention.
/// - 2013.09.05 JM: changed RS position[] to pos[].
/// - 2013.10.13 JM: Tidied up a bit.
/// - 2015.[01.26-02.03] JM: CHANGED FILENAME TO SIM_CONTROL_MPI.CPP,
///    added ParallelParams class, and fixing code for non-global mpiPM.
/// - 2015.02.18 JM: moved setup functions to setup_fixed_grid_MPI
/// - 2016.03.14 JM: Worked on parallel Grid_v2 update (full
///    boundaries).  Changed default I/O to DOUBLE precision.
/// - 2017.08.03 JM: Changed silo dataio class to be the utility class. 

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "MCMD_control.h"
#include "sim_control_MPI.h"
#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_utility.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"
#endif // if FITS

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#ifdef PARALLEL




// ##################################################################
// ##################################################################


sim_control_fixedgrid_pllel::sim_control_fixedgrid_pllel()
  : sim_control_fixedgrid()
{
#ifdef TESTING
  cout <<"sim_control_fixedgrid_pllel constructor.\n";
#endif
}



// ##################################################################
// ##################################################################


sim_control_fixedgrid_pllel::~sim_control_fixedgrid_pllel()
{
#ifdef TESTING    
  cout <<"sim_control_fixedgrid_pllel destructor.\n";
#endif
}



// ##################################################################
// ##################################################################



int sim_control_fixedgrid_pllel::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      class GridBaseClass **grid
      )
{

  /// \userguide \section icfiles Parallel IC/Restart files.
  ///
  /// [FITS IS NOT COMPILED BY DEFAULT ANYMORE; USE SILO; SEE BELOW] 
  ///
  /// When running in parallel, outputs are named according to the rank as
  /// follows: \'outfile\_\<myrank\>.\<timestep\>.fits\'. Because of this, an
  /// initial condition file should never have the string \'_\<number\>.\' in
  /// it.  This will confuse the data I/O routines.  If reading from a single
  /// initial conditions (or restart) file, this string should not be in the
  /// filename, and all processes will read from the same file.  If restarting
  /// from multiple input files, however, the given input filename should be the
  /// restart file for processor 0 (i.e. restartfile\_0.\<timestep\>.fits), and
  /// the function sim_control_fixedgrid_pllel::init() will parse the input
  /// filename and replace the \"\_0.\" with \"\_\<myrank\>.\", followed by a
  /// call to the original serial function sim_control_fixedgrid::init() with the new
  /// filename, if it exists.
  ///
  /// For Silo data I/O the model is a little different.  The number of files
  /// can be smaller than the number of processors.  Serial files should still
  /// never have the substring "\_0" in them.  Parallel files are named in a
  /// similar way as for fits files, but each file can contain data from a
  /// number of processors, stored in different subdirectories within the file.
  /// Also, filenames are now stored as follows: \'ouftile\_0000.00000.silo\'
  /// where the first zeros are for the number of the file, and the second are
  /// for the timestep.  It made sense to have a fixed width number, for listing
  /// files as much as anything else.

#ifdef TESTING
  cout <<"(sim_control_fixedgrid_pllel::init) Initialising grid: infile = "<<infile<<"\n";
#endif
  int err=0;

  //
  // Setup the MCMDcontrol class with rank and nproc.
  //
  int myrank = -1, nproc = -1;
  COMM->get_rank_nproc(&myrank, &nproc);
  mpiPM.set_myrank(myrank);
  mpiPM.set_nproc(nproc);


  //
  // Now find out what kind of file we are to read, and then read it.
  //
  if (typeOfFile==1) {
      rep.error("Don't give me text file for parallel I/O! Crazy fool!",typeOfFile);
  }
  // -------------------------------

#ifdef FITS
  // ******** FITS FILE I/O ********
  else if (typeOfFile==2) {
    if (!dataio) dataio = new DataIOFits_pllel(&mpiPM);
    if (!dataio) rep.error("DataIOFits initialisation",dataio);
    if (dataio->file_exists(infile)) {
#ifdef TESTING
      cout <<"\t Reading from file : "<< infile<<"\n";
      cout <<"\t Assume if filename contains \'_0000.\', that it is the first of multiple files.\n";
#endif
      string::size_type pos =infile.find("_0000.");
      if (pos==string::npos) {
#ifdef TESTING
	cout <<"\t Couldn't find \'_0000.\' in file, so reading from single file.\n";
#endif
	mpiPM.ReadSingleFile = true;
      }
      else {
#ifdef TESTING
	cout <<"\t Found \'_0000.\' in file, so replacing that with myrank.\n";
#endif
	mpiPM.ReadSingleFile = false;
#ifdef TESTING
	cout <<"\t Old infile: "<<infile<<"\n";
#endif
	ostringstream t; t.str("");
        t<<"_"; t.width(4); t.fill('0'); t<<mpiPM.get_myrank()<<".";
        string t2=t.str();
	infile.replace(pos,6,t2);
#ifdef TESTING
	cout <<"\t New infile: "<<infile<<"\n";
#endif
	if (!dataio->file_exists(infile))
          rep.error("infile doesn't exist!",infile);
      }
    }
    else {
      cout <<"\tInfile doesn't exist: failing\n";
      return(1);
    }
  }
  // ******** FITS FILE I/O ********
#endif // FITS
  // -------------------------------

#ifdef SILO
  // -------------------------------
  // ******** SILO FILE I/0 ********
  else if (typeOfFile==5) {
    string::size_type pos =infile.find("_0");
    if (pos==string::npos) {
#ifdef TESTING
      cout <<"\t Couldn't find \'_0\' in file, so assume reading from single file.\n";
#endif
      //cout <<"\t But using parallel I/O, so need parallel multimesh objects... failing!\n";
      //return 99;
      mpiPM.ReadSingleFile = true;
    }
    else {
#ifdef TESTING
      cout <<"\t Found \'_0\' in file, so assume reading from multiple files.\n";
      cout <<"silo class knows how to get to the right file name when reading data.\n";
#endif
      mpiPM.ReadSingleFile = false;
    }

    if (!dataio) dataio = new dataio_silo_utility ("DOUBLE", &mpiPM);
    if (!dataio) rep.error("dataio_silo_pllel initialisation",dataio);
    if (!dataio->file_exists(infile)) {
      cout <<"\tInfile doesn't exist: failing\n";
      return(1);
    }
  }  
  // ******** SILO FILE I/0 ********
  // -------------------------------
#endif // if SILO

  else
    rep.error("Bad file type specifier for parallel grids (2=fits,\
               5=silo) IS IT COMPILED IN???",typeOfFile);
  
  //
  // We need to decompose the domain here, because setup_grid() needs
  // this, but this means we need to read the header to find out what
  // the grid dimensions are.  So the header is read twice, but this
  // should be ok because it only happens during initialisation.
  //
  err = dataio->ReadHeader(infile);
  if (err) rep.error("PLLEL Init(): failed to read header",err);
  err = mpiPM.decomposeDomain();
  if (err) rep.error("PLLEL Init():Couldn't Decompose Domain!",err);



#ifdef TESTING
  cout <<"(sim_control_fixedgrid_pllel::init) Calling serial code sim_control_fixedgrid::init() on infile."<<"\n";
#endif
  err = sim_control_fixedgrid::Init(infile, typeOfFile, narg, args, grid);
  if (err) rep.error("failed to do serial init",err);



  // If outfile-type is different to infile-type, we need to delete dataio and set it up again.
  // This is ifdeffed out in the serial code because parallel I/O classes need to be set up.
  if (SimPM.typeofip != SimPM.typeofop) {
#ifdef TESTING
    cout <<"(sim_control_fixedgrid_pllel::INIT) infile-type="<<SimPM.typeofip;
    cout <<" and outfile-type="<<SimPM.typeofop;
    cout <<", so deleting and renewing dataio.\n";
#endif

    if (dataio) {delete dataio; dataio=0;}
    switch (SimPM.typeofop) {
    case 1: // ascii, so do nothing.
      break;
#ifdef FITS
    case 2: // fits
    case 3: // fits
    case 4: // fits +ascii
      dataio = new DataIOFits_pllel(&mpiPM);
      break;
#endif
#ifdef SILO
    case 5: // silo
      dataio = new dataio_silo_utility ("DOUBLE", &mpiPM);
      break;
#endif // if SILO
    default:
      rep.error("PARALLEL INIT: Bad type of Output",SimPM.typeofop);
    }
    if (!dataio) rep.error("INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(eqn);
  if (SimPM.timestep==0) {
#ifdef TESTING
     cout << "(P\'LLEL INIT) Outputting initial data.\n";
#endif
     output_data(*grid);
  }
  cout <<"                                   ******************************\n";

  return(err);
}



// ##################################################################
// ##################################################################




/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int sim_control_fixedgrid_pllel::Time_Int(
        class GridBaseClass *grid
        )
{
  cout <<"                               **************************************\n";
  cout <<"(sim_control_fixedgrid_pllel::time_int) STARTING TIME INTEGRATION."<<"\n";
  int err=0;
  int log_freq=10;
  SimPM.maxtime=false;
  clk.start_timer("time_int"); double tsf=0.0;
  while (SimPM.maxtime==false) {
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources();
    if (err) {
      cerr <<"(TIME_INT::update_evolving_RT_sources()) something went wrong!\n";
      return err;
    }

    //clk.start_timer("advance_time");
    err+= advance_time(grid);
    //cout <<"advance_time took "<<clk.stop_timer("advance_time")<<" secs.\n";
    if (err!=0) {
      cerr<<"(TIME_INT::advance_time) err!=0 Something went bad"<<"\n";
      return(1);
      }

    if (mpiPM.get_myrank()==0 && (SimPM.timestep%log_freq)==0) {
      cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime;
      cout <<"\t timestep: "<<SimPM.timestep;
      tsf=clk.time_so_far("time_int");
      cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#ifdef TESTING
      cout.flush();
#endif // TESTING
    }
	
    //
    // check if we are at time limit yet.
    //
    tsf=clk.time_so_far("time_int");
    double maxt = COMM->global_operation_double("MAX", tsf);
    if (maxt > get_max_walltime()) {
      SimPM.maxtime=true;
      cout <<"RUNTIME>"<<get_max_walltime()<<" SECS.\n";
    }
	
    err+= output_data(grid);
    if (err!=0){
      cerr<<"(TIME_INT::output_data) err!=0 Something went bad"<<"\n";
      return(1);
    }

    err+= check_eosim();
    if (err!=0) {
      cerr<<"(TIME_INT::) err!=0 Something went bad"<<"\n";
      return(1);
    }
  }
  cout <<"(sim_control_fixedgrid_pllel::time_int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM."<<"\n";
  tsf=clk.time_so_far("time_int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep;
  cout <<", sim-time="<<SimPM.simtime;
  cout <<", wall-time=" <<tsf;
  cout <<", time/step="<<tsf/static_cast<double>(SimPM.timestep) <<"\n";
  if (RT!=0) {
    //
    // output raytracing timing info.  Have to start and stop timers to get 
    // the correct runtime (this is sort of a bug... there is no function to
    // get the elapsed time of a non-running timer.  I should add that.
    //
    string t1="totalRT", t2="waitingRT", t3="doingRT";
    double total=0.0, wait=0.0, run=0.0;
    clk.start_timer(t1); total = clk.pause_timer(t1);
    clk.start_timer(t2); wait  = clk.pause_timer(t2);
    clk.start_timer(t3); run   = clk.pause_timer(t3);
    cout <<"TOTALS RT#: active="<<run<<" idle="<<wait<<" total="<<total<<"\n";
  }
  cout <<"                               **************************************\n\n";
  return(0);
}



// ##################################################################
// ##################################################################


int sim_control_fixedgrid_pllel::calc_timestep(
        class GridBaseClass *grid
        )
{
  //
  // First get the local grid dynamics and microphysics timesteps.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt(grid);
  t_mp  = calc_microphysics_dt(grid);
  
  //
  // Now get global min over all grids for dynamics and microphysics timesteps.
  // We only need both if we are doing Dedner et al. 2002, mixed-GLM divergence
  // cleaning of the magnetic field, since there the dynamical dt is an important
  // quantity (see next block of code).
  //
  //SimPM.dt = t_dyn;
  t_dyn = COMM->global_operation_double("MIN", t_dyn);
  //cout <<"proc "<<mpiPM.myrank<<":\t my t_dyn="<<SimPM.dt<<" and global t_dyn="<<t_dyn<<"\n";
  //SimPM.dt = t_mp;
  t_mp = COMM->global_operation_double("MIN", t_mp);
  //cout <<"proc "<<mpiPM.myrank<<":\t my t_mp ="<<SimPM.dt<<" and global t_mp ="<<t_mp<<"\n";
  
  // Write step-limiting info every tenth timestep.
  if (t_mp<t_dyn && (SimPM.timestep%10)==0)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";

  //
  // if using MHD with GLM divB cleaning, the following sets the hyperbolic wavespeed.
  // If not, it does nothing.  By setting it here and using t_dyn, we ensure that the
  // hyperbolic wavespeed is equal to the maximum signal speed on the grid, and not
  // an artificially larger speed associated with a shortened timestep.
  //
  eqn->GotTimestep(t_dyn);

  //
  // Now the timestep is the min of the global microphysics and Dynamics timesteps.
  //
  SimPM.dt = min(t_dyn,t_mp);

#ifdef THERMAL_CONDUCTION
  //
  // In order to calculate the timestep limit imposed by thermal conduction,
  // we need to actually calcuate the multidimensional energy fluxes
  // associated with it.  So we store Edot in c->dU[ERG], to be multiplied
  // by the actual dt later (since at this stage we don't know dt).  This
  // later multiplication is done in eqn->preprocess_data()
  //
  double t_cond = calc_conduction_dt_and_Edot();
  t_cond = COMM->global_operation_double("MIN", t_cond);
  if (t_cond<SimPM.dt) {
    cout <<"PARALLEL CONDUCTION IS LIMITING TIMESTEP: t_c="<<t_cond<<", t_m="<<t_mp;
    cout <<", t_dyn="<<t_dyn<<"\n";
  }
  SimPM.dt = min(SimPM.dt, t_cond);
#endif // THERMAL CONDUCTION

  //
  // Check that the timestep doesn't increase too much between step, and that it 
  // won't bring us past the next output time or the end of the simulation.
  // This function operates on SimPM.dt, resetting it to a smaller value if needed.
  //
  timestep_checking_and_limiting();
  
  //
  // Tell the solver class what the resulting timestep is.
  //
  eqn->Setdt(SimPM.dt);
  
#ifdef TESTING
  //
  // Check (sanity) that if my processor has modified dt to get to either
  // an output time or finishtime, then all processors have done this too!
  // This is really paranoid testing, and involves communication, so only 
  // do it when testing.
  //
  t_dyn=SimPM.dt; t_mp = COMM->global_operation_double("MIN", t_dyn);
  if (!pconst.equalD(t_dyn,t_mp))
    rep.error("synchonisation trouble in timesteps!",t_dyn-t_mp);
#endif // TESTING

  return 0;
}


// ##################################################################
// ##################################################################



#endif // PARALLEL


