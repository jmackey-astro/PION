/** \file gridMethodsMPI.cc
 * 
 * \brief Parallel Grid Methods Class Member Function definitions.
 * 
 * \author Jonathan Mackey
 * 
 * This file contains the definitions of the member functions for ParallelIntUniformFV
 * class, which is a modification of the basic 1st/2nd order Finite Volume Solver according to
 * the method outlined in Falle, Komissarov, \& Joarder (1998), MNRAS, 297, 265.
 * 
 * Modifications:
 *  - 2007-10-11 Started writing file.
 *  - 2007-11-11 Basically working with uniform grid and all sorts of boundaries.
 *  - 2010-01-05 JM: Put RT step limiting in an #ifdef so that it is not set for test problems.
 * */
///
/// - 2010-04-21 JM: Changed filename setup so that i can write
///   checkpoint files with fname.999999.txt/silo/fits.  This also
///   meant I don't need a new parallel output_data() function so I
///   got rid of it.  All dataio classes now choose the filename
///   themselves.
///
/// - 2010.07.23 JM: New RSP source position class interface.
///
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data.
///
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///
/// - 2010.12.15 JM: Added CI.setup_extra_data() call in setup_grid()
///     function.
///
/// - 2011.01.12 JM: Set so only proc 0 displays RT timestep.
///     Commented out some stuff in calc_timestep which was used for
///     testing and made the code less efficient.
///
/// - 2011.02.17 JM: Raytracer header file moved to raytracing/ subdir.
///     More ray-tracing options for CI.setup_extra_data
///
/// - 2011.02.25 JM: removed HCORR ifdef around new code. Modified setup_extra_data
///     for ray-tracing again.  Now we have N column-density vars for N sources.
///
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

#include "grid.h"
#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#ifdef PARALLEL



// ##################################################################
// ##################################################################


ParallelIntUniformFV::ParallelIntUniformFV()
  : IntUniformFV()
{
#ifdef TESTING
  cout <<"ParallelIntUniformFV constructor.\n";
#endif
}



// ##################################################################
// ##################################################################


ParallelIntUniformFV::~ParallelIntUniformFV()
{
#ifdef TESTING    
  cout <<"ParallelIntUniformFV destructor.\n";
#endif
}



// ##################################################################
// ##################################################################



int ParallelIntUniformFV::setup_grid()
{
#ifdef TESTING
  cout <<"ParallelIntUniformFV: setting up parallel grid.\n";
#endif

  if (SimPM.gridType!=1) {
    rep.warning("gridType not set correctly: Only know Uniform finite volume grid, so resetting to 1!",1,SimPM.gridType);
    SimPM.gridType=1;
  }
  if (SimPM.ndim <1 || SimPM.ndim>3)  rep.error("Only know 1D,2D,3D methods!",SimPM.ndim);
  
  // First decompose the domain, so I know the dimensions of the local grid to set up.
  int err=0;
  if((err=mpiPM.decomposeDomain()))
    rep.error("Couldn't Decompose Domain!",err);
  
  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.  Cells cannot be created unless this
  // the number of such extra variables has been set.
  //
  setup_cell_extra_data();

  //
  // Now set up the parallel uniform grid.
  //
#ifdef TESTING
  cout <<"(ParallelIntUniformFV::setup_grid) Setting up grid...\n";
#endif

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


  if (grid==0) rep.error("(ParallelIntUniformFV::setup_grid) Couldn't assign data!", grid);

#ifdef TESTING
  cout <<"(ParallelIntUniformFV::setup_grid) Done. grid="<<grid;//<<"\n";
  cout <<"\t DX = "<<grid->DX()<<"\n";
#endif

  return(0);
}



// ##################################################################
// ##################################################################



int ParallelIntUniformFV::Init(string infile, int typeOfFile, int narg, string *args)
{
    /** \userguide
     * \section icfiles Parallel IC/Restart files.
     * [FITS IS NOT COMPILED BY DEFAULT ANYMORE; USE SILO; SEE BELOW] 
     * When running in parallel, outputs are named according to the rank as follows:
     * \'outfile\_\<myrank\>.\<timestep\>.fits\'. Because of this, an initial condition
     * file should never have the string \'_\<number\>.\' in it.  This will confuse the
     * data I/O routines.  If reading from a single initial conditions (or restart) file,
     * this string should not be in the filename, and all processes will read from the
     * same file.  If restarting from multiple input files, however, the given input filename
     * should be the restart file for processor 0 (i.e. restartfile\_0.\<timestep\>.fits),
     * and the function ParallelIntUniformFV::init() will parse the input filename and
     * replace the \"\_0.\" with \"\_\<myrank\>.\", followed by a call to the original serial
     * function IntUniformFV::init() with the new filename, if it exists.
     *
     * For Silo data I/O the model is a little different.  The number of files can be 
     * smaller than the number of processors.  Serial files should still never have the
     * substring "\_0" in them.  Parallel files are named in a similar way as for fits files,
     * but each file can contain data from a number of processors, stored in different
     * subdirectories within the file.  Also, filenames are now stored as follows: 
     * \'ouftile\_0000.00000.silo\' where the first zeros are for the number of the file, and
     * the second are for the timestep.  It made sense to have a fixed width number, for listing
     * files as much as anything else.
     */

#ifdef TESTING
  cout <<"(ParallelIntUniformFV::init) Initialising grid: infile = "<<infile<<"\n";
#endif
  int err=0;

  if (typeOfFile==1) {
      rep.error("Don't give me text file for parallel I/O! Crazy fool!",typeOfFile);
  }
  // -------------------------------

#ifdef FITS
  // ******** FITS FILE I/O ********
  else if (typeOfFile==2) {
    if (!dataio) dataio = new DataIOFits();
    if (!dataio) rep.error("DataIOFits initialisation",dataio);
    if (dataio->file_exists(infile)) {
#ifdef TESTING
      cout <<"\t Reading from file : "<< infile<<"\n";
      cout <<"\t Assume if filename contains \'_0.\', that it is the first of multiple files.\n";
#endif
      string::size_type pos =infile.find("_0.");
      if (pos==string::npos) {
#ifdef TESTING
	cout <<"\t Couldn't find \'_0.\' in file, so reading from single file.\n";
#endif
	mpiPM.ReadSingleFile = true;
      }
      else {
#ifdef TESTING
	cout <<"\t Found \'_0.\' in file, so replacing that with myrank.\n";
#endif
	mpiPM.ReadSingleFile = false;
#ifdef TESTING
	cout <<"\t Old infile: "<<infile<<"\n";
#endif
	ostringstream t; t.str(""); t<<"_"<<mpiPM.myrank<<"."; string t2=t.str();
	infile.replace(pos,3,t2);
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

    if (!dataio) dataio = new dataio_silo_pllel ();
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
    rep.error("Bad file type specifier for parallel grids (2=fits,5=silo) IS IT COMPILED IN???",typeOfFile);

#ifdef TESTING
  cout <<"(ParallelIntUniformFV::init) Calling serial code IntUniformFV::init() on infile."<<"\n";
#endif
  err=IntUniformFV::Init(infile,typeOfFile,narg,args);
  if (err) rep.error("failed to do serial init",err);



  // If outfile-type is different to infile-type, we need to delete dataio and set it up again.
  // This is ifdeffed out in the serial code because parallel I/O classes need to be set up.
  if (SimPM.typeofip != SimPM.typeofop) {
#ifdef TESTING
    cout <<"(ParallelIntUniformFV::INIT) infile-type="<<SimPM.typeofip;
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
      dataio = new DataIOFits();
      break;
#endif
#ifdef SILO
    case 5: // silo
      dataio = new dataio_silo_pllel ();
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
     output_data();
  }
  cout <<"                                   ******************************\n";

  return(err);
}



// ##################################################################
// ##################################################################


int ParallelIntUniformFV::setup_raytracing()
{
  //
  // This function is basically identical to the serial setup function, except
  // that it sets up parallelised versions of the raytracers.
  //
  // If not doing raytracing, return immediately.
  //
  if (!SimPM.EP.raytracing) {
    return 0;
  }

  //
  // Now we are doing raytracing, so set up a raytracer and add sources to it.
  //
  if (!MP) rep.error("can't do raytracing without microphysics",MP);
  cout <<"\n***************** RAYTRACER SETUP STARTING ***********************\n";
  RT=0;
  //
  // If the ionising source is at infinity then set up the simpler parallel
  // rays tracer.  Otherwise the more complicated one is required.
  //
  bool parallel_rays=true;
  int dir=-1;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (!SimPM.RS.sources[isrc].at_infinity) parallel_rays=false;
    //
    // source is at infinity, so make sure all sources at infinity have rays
    // travelling in the same direction (by checking direction to source).
    //
    else {
      for (int i=0;i<SimPM.ndim;i++) {
        if (fabs(SimPM.RS.sources[isrc].position[i])>1.e99) {
          if (dir==-1) dir=i;
          else if (dir!=i) parallel_rays=false;
        }
      }
    }
  }     // loop over sources.
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  parallel_rays=false;
  // HACK -- DISABLE PARALLEL RAYS APPROX ALWAYS SO I CAN DO NORMAL
  // DOMAIN DECOMPOSITION.
  if (parallel_rays) {
    //
    // set up single source at infinity tracer, if appropriate
    //
    RT = new raytracer_USC_infinity(grid,MP);
    if (!RT) rep.error("init pllel-rays raytracer error",RT);
  }
  else {
    //
    // set up regular tracer if simple one not already set up.
    //
    RT = new raytracer_USC_pllel(grid,MP);
    if (!RT) rep.error("init raytracer error 2",RT);
  }

  //
  // Now add the sources to the tracer.  Note that both the implicit and explicit
  // integrators can still only handle a single ionising source, so we do a check
  // for this and bug out if there is more than one.
  //
  int ion_count=0, uv_count=0, dif_count=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) {
      //
      // single sources have a flux (if at infinity) or a luminosity (if point
      // sources.
      //
      cout <<"Adding IONISING or UV single-source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MONO ||
          SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI)
        ion_count++;
      else 
        uv_count++;
    } // if ionising source
    else {
      // note that diffuse radiation must be at infinity, and the strength is assumed to
      // be an intensity not a flux, so it is multiplied by a solid angle appropriate
      // to its location in order to get a flux.
      cout <<"Adding DIFFUSE radiation source with id: ";
      cout << RT->Add_Source(&(SimPM.RS.sources[isrc])) <<"\n";
      uv_count++;
      dif_count++;
    } // if diffuse source
  } // loop over sources
  if (ion_count>1) {
    rep.error("Can only have one ionising source for currently implemented method",ion_count);
  }
  cout <<"Added "<<ion_count<<" ionising and "<<uv_count<<" non-ionising";
  cout <<" radiation sources, of which "<<dif_count<<" are diffuse radiation.\n";
  RT->Print_SourceList();

  //
  // Now that we have added all of the sources, we query the raytracer to get
  // all of the source properties into structs for the microphysics calls.
  // NOTE THAT IF THE NUMBER OF SOURCES OR THEIR PROPERTIES CHANGE OVER TIME,
  // I WILL HAVE TO WRITE NEW CODE TO UPDATE THIS!
  //
  FVI_nheat = RT->N_heating_sources();
  FVI_nion  = RT->N_ionising_sources();
  FVI_heating_srcs.resize(FVI_nheat);
  FVI_ionising_srcs.resize(FVI_nion);
  RT->populate_UVheating_src_list(FVI_heating_srcs);
  RT->populate_ionising_src_list( FVI_ionising_srcs);

  //
  // See if we need column densities for the timestep calculation
  //
  if (RT->type_of_RT_integration()==RT_UPDATE_EXPLICIT) {
    FVI_need_column_densities_4dt = true;
  }
  else if (RT && RT->type_of_RT_integration()==RT_UPDATE_IMPLICIT
            && SimPM.EP.MP_timestep_limit==5) {
    // For implicit updates to limit by xdot and/or edot
    // Here the raytracing has not already been done, so we call it here.
    FVI_need_column_densities_4dt = true;
  }
  else {
    FVI_need_column_densities_4dt = false;
  }
  
  cout <<"***************** RAYTRACER SETUP ***********************\n";
  return 0;
}



// ##################################################################
// ##################################################################




/*****************************************************************/
/*********************** TIME INTEGRATION ************************/
/*****************************************************************/
int ParallelIntUniformFV::Time_Int()
{
    cout <<"                               **************************************\n";
    cout <<"(ParallelIntUniformFV::time_int) STARTING TIME INTEGRATION."<<"\n";
    int err=0;
    SimPM.maxtime=false;
    GS.start_timer("time_int"); double tsf=0.0;
    while (SimPM.maxtime==false) {
        //
        // Update RT sources.
        //
        err = update_evolving_RT_sources();
        if (err) {
          cerr <<"(TIME_INT::update_evolving_RT_sources()) something went wrong!\n";
          return err;
        }

	//GS.start_timer("advance_time");
	err+= advance_time();
	//cout <<"advance_time took "<<GS.stop_timer("advance_time")<<" secs.\n";
	if (err!=0){cerr<<"(TIME_INT::advance_time) err!=0 Something went bad"<<"\n";return(1);}

	if (mpiPM.myrank==0 && (SimPM.timestep%10)==0) {
	  cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
	  tsf=GS.time_so_far("time_int");
	  cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
	}
	
        // check if we are at time limit yet.
	double maxt = COMM->global_operation_double("MAX", tsf);

	if (maxt > mpiPM.get_max_walltime()) {
	    SimPM.maxtime=true;
	    cout <<"RUNTIME>"<<mpiPM.get_max_walltime()<<" SECS.\n";
	}
	
	err+= output_data();
	if (err!=0){cerr<<"(TIME_INT::output_data) err!=0 Something went bad"<<"\n";return(1);}
	err+= check_eosim();
	if (err!=0){cerr<<"(TIME_INT::) err!=0 Something went bad"<<"\n";return(1);}
    }
    cout <<"(ParallelIntUniformFV::time_int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM."<<"\n";
    tsf=GS.time_so_far("time_int");
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
      GS.start_timer(t1); total = GS.pause_timer(t1);
      GS.start_timer(t2); wait  = GS.pause_timer(t2);
      GS.start_timer(t3); run   = GS.pause_timer(t3);
      cout <<"TOTALS RT#: active="<<run<<" idle="<<wait<<" total="<<total<<"\n";
    }
    cout <<"                               **************************************\n\n";
    return(0);
}



// ##################################################################
// ##################################################################


int ParallelIntUniformFV::calc_timestep()
{
  //
  // First get the local grid dynamics and microphysics timesteps.
  //
  double t_dyn=0.0, t_mp=0.0;
  t_dyn = calc_dynamics_dt();
  t_mp  = calc_microphysics_dt();
  // output step-limiting info every tenth timestep.
  if (t_mp<t_dyn && (SimPM.timestep%10)==0)
    cout <<"Limiting timestep by MP: mp_t="<<t_mp<<"\thydro_t="<<t_dyn<<"\n";
  
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
  if (!GS.equalD(t_dyn,t_mp))
    rep.error("synchonisation trouble in timesteps!",t_dyn-t_mp);
#endif // TESTING

  return 0;
}


// ##################################################################
// ##################################################################



#endif // PARALLEL


