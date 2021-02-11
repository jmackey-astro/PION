/// \file icgen_parallel.cpp
/// \author Jonathan Mackey
/// \date 2018.05.04
///
/// Parallel version of the initial conditions generator.
/// Moved code from serial version.
///
/// Modifications:
/// - 2018.05.04 JM: moved from icgen.cpp


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen_base.h"
#include "ics/icgen.h"
#include "ics/get_sim_info.h"

#include "dataIO/dataio_base.h"
#ifdef FITS
#include "dataIO/dataio_fits_MPI.h"
#endif // if FITS
#ifdef SILO
#include "dataIO/dataio_silo_MPI.h"
#endif // if SILO

#include "grid/grid_base_class.h"
#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "grid/setup_fixed_grid_MPI.h"

#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;


// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  int err = COMM->init(&argc, &argv);
  if (err) rep.error("comms init error",err);

  int r=-1, np=-1;
  COMM->get_rank_nproc(&r,&np);
  class SimParams SimPM;
  SimPM.levels.clear();
  SimPM.levels.resize(1);
  SimPM.levels[0].MCMD.set_myrank(r);
  SimPM.levels[0].MCMD.set_nproc(np);

  if (argc<2) {
    cerr<<"Error, please give a filename to read IC parameters from.\n";
    cerr<<"Usage <icgen> <paramfile> [ic-filetype]\n";
    return(1);
  }

  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) args[i] = argv[i];
  // Redirect stdout/stderr if required.
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      ostringstream path;
      path << outpath <<"_"<<SimPM.levels[0].MCMD.get_myrank()<<"_";
      outpath = path.str();
      if (SimPM.levels[0].MCMD.get_myrank()==0) {
        cout <<"Redirecting stdout to "<<outpath<<"info.txt"<<"\n";
      }
      rep.redirect(outpath); // Redirects cout and cerr to text file.
    }
  }
#ifndef TESTING
  rep.kill_stdout_from_other_procs(0);
#endif
  //cout << "rank: " << MCMD.get_myrank();
  //cout << " nproc: " << MCMD.get_nproc() << "\n";

  class DataIOBase   *dataio=0;
  class get_sim_info *siminfo=0;
  class ICsetup_base *ic=0;
  class ReadParams   *rp=0;
  MP=0;  // global microphysics class pointer.

  string pfile = argv[1];
  string icftype;
  if (argc>2) icftype=argv[2];
  else icftype="fits"; // This is the default for now.
  
  siminfo=0; siminfo = new class get_sim_info ();
  if (!siminfo) rep.error("Sim Info class init error",siminfo);
  err = 0;
  err += siminfo->read_gridparams(pfile, SimPM);
  if (err) rep.error("Read Grid Params Error",err);
  delete siminfo; siminfo=0;

  // have to do something with SimPM.levels[0] because this
  // is used to set the local domain size in decomposeDomain
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


  class setup_fixed_grid *SimSetup =0;
  SimSetup = new setup_fixed_grid_pllel();

  //
  // Set up the Xmin/Xmax/Range/dx of each level in the NG grid
  //
  vector<class GridBaseClass *> grid;
  // have to do something with SimPM.levels[0] because this
  // is used to set the local domain size in decomposeDomain
  err = SimPM.levels[0].MCMD.decomposeDomain(SimPM, SimPM.levels[0]);
  rep.errorTest("Couldn't Decompose Domain!",0,err);

  class MCMDcontrol *MCMD = &(SimPM.levels[0].MCMD);

  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  grid.resize(1);
  // Now set up the grid structure.
  err = SimSetup->setup_grid(grid,SimPM);
  SimPM.dx = grid[0]->DX();
  if (!grid[0]) rep.error("Grid setup failed",grid[0]);
  
  //
  // read in what kind of ICs we are setting up.
  //
  rp = new ReadParams;
  if (!rp) rep.error("icgen:: initialising RP",rp);
  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);

  //
  // Assign data to the grid, so we call whichever
  // function is requested.
  //
  string seek="ics";
  string ics = rp->find_parameter(seek);
  setup_ics_type(ics,&ic);
  ic->set_SimPM(&SimPM);
  ic->set_MCMD_pointer(MCMD);


  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    // don't need to set up the class, because it just does cooling and
    // there is no need to equilibrate anything.
  }
  else if (SimPM.ntracer>0 && (SimPM.EP.cooling || SimPM.EP.chemistry)) {
    SimSetup->setup_microphysics(SimPM);
    if (!MP) rep.error("microphysics init",MP);
  }
  // ----------------------------------------------------------------

  // have to setup jet simulation before setting up boundary
  // conditions because jet boundary needs some grid data values.
  if (ics=="Jet" || ics=="JET" || ics=="jet") {
    ic->setup_data(rp,grid[0]);
  }

  //
  // Set up the boundary conditions, since internal boundary data
  // should really be already set to its correct value in the initial
  // conditions file.
  //
  SimSetup->boundary_conditions(SimPM,grid);
  if (err) rep.error("icgen Couldn't set up boundaries.",err);
  err += SimSetup->setup_raytracing(SimPM, grid[0]);
  err += SimSetup->setup_evolving_RT_sources(SimPM);
  err += SimSetup->update_evolving_RT_sources(SimPM,SimPM.simtime,
                                                      grid[0]->RT);
  if (err) rep.error("icgen: Failed to setup raytracer and/or microphysics",err);

  // ----------------------------------------------------------------
  // call "setup" to set up the data on the computational grid.
  err += ic->setup_data(rp,grid[0]);
  if (err) rep.error("Initial conditions setup failed.",err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // if data initialised ok, maybe we need to equilibrate the 
  // chemistry...
  //
  if (SimPM.ntracer>0 && (SimPM.EP.chemistry)) {
    cout <<"MAIN: equilibrating the chemical species.\n";
    if (!MP) rep.error("microphysics init",MP);

    // first avoid cooling the gas in getting to equilbrium, by
    // setting update_erg to false.
    bool uerg = SimPM.EP.update_erg;
    SimPM.EP.update_erg = false;

    err = ic->equilibrate_MP(grid[0],MP,rp,SimPM);
    if (err)
      rep.error("setting chemical states to equilibrium failed",err);

    SimPM.EP.update_erg = uerg;
    cout <<"MAIN: finished equilibrating the chemical species.\n";
  }
  // ----------------------------------------------------------------


  cout <<"IC file-type is "<<icftype<<"\n";
  seek="OutputFile";
  string icfile = rp->find_parameter(seek);
  if (icfile=="") {
    cout <<"WARNING: no filename for ic file.  writing to IC_temp.****\n";
    icfile = "IC_temp";
  }
  dataio = 0; // zero the class pointer.

#ifdef FITS    
  if (icftype=="fits") {
    cout <<"WRITING FITS FILE: ";
    cout << icfile << "\n";
    dataio = 0; dataio = new DataIOFits_pllel (SimPM, MCMD);
  }
#endif // if fits.

#ifdef SILO
  if (icftype=="silo") {
    cout <<"WRITING SILO FILE: ";
    cout <<icfile <<"\n";
    dataio=0; dataio=new dataio_silo_pllel (SimPM, "DOUBLE",MCMD);
  }
#endif // if SILO defined.

  if (!dataio) rep.error("IO class initialisation: ",icftype);
  err = dataio->OutputData(icfile,grid, SimPM, 0);
  if (err) rep.error("File write error",err);
  delete dataio; dataio=0;
  cout <<icftype<<" FILE WRITTEN in";

  // delete everything and return
  if (MP)   {delete MP; MP=0;}
  if (rp)   {delete rp; rp=0;} // Delete the read_parameters class.
  if (ic)   {delete ic; ic=0;}
  if (SimSetup) {delete SimSetup; SimSetup =0;}
  delete grid[0];

  //
  // Also delete any dynamic memory in the stellarwind_list in
  // global.h, since we may have added wind sources in the
  // read-params functions.
  //
  while (SWP.params.size()>0) {
    int i=static_cast<int>(SWP.params.size())-1;
    struct stellarwind_params *temp = SWP.params[i];
    SWP.params.pop_back();   // remove struct from list.
    temp = mem.myfree(temp); // delete struct.
  }

  cout << "rank: " << MCMD->get_myrank();
  cout << " nproc: " << MCMD->get_nproc() << "\n";
  COMM->finalise();
  delete COMM; COMM=0;
  delete [] args; args=0;
  return err;
}


// ##################################################################
// ##################################################################





