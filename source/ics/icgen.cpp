/// \file icgen.cpp
/// \brief Program to generate Initial Conditions for a uniform grid
/// \author Jonathan Mackey
/// 
/// Modifications:
///  - 2007-06-28 Started to write an ND shocktube function.  Need to test it.
///  - 2007-07-17 ND MHD shock-cloud interaction problem added to list.  Works well.
///  - 2007-07-24 Added passive tracer variable support.
///  - 2007-08 to 2008-02 Added more functions.
///  - 2009-12-18 Added in Laser Ablation check.
/// - 2010.12.04 JM: Added geometry-dependent grids, in a
///   GEOMETRIC_GRID ifdef. (only serial so far).
/// - 2010.01.06 JM: New stellar wind interface.
/// - 2011.04.14 JM: Added mp_v2_aifa microphysics integration class.
/// - 2011.04.26 JM: removed endl statements.
/// - 2011.04.29 JM: in AddNoise2data() and equilibrate_MP() functions, check
///    for edge cells and internal boundary data, and don't add noise to these
///    cells.  If we have e.g. inflow boundaries, which are set by the edge
///    cell value, we don't that to be a random number!
///
/// - 2011.05.02 JM: AddNoise2data() -- don't add noise to MPI boundaries.
///
/// - 2011.10.13 JM: Added mp_[ex,im]plicit_H classes. (updated 2011.10.22)
/// - 2012.01.23 JM: Added ifdefs for microphysics classes.
/// - 2012.02.07 JM: Added class for Harpreet's 1D to 2D mapping.
/// - 2012.07.25 JM: Fixed bug where noise was added to edge cells in the
///    YZ-face for parallel grids.
/// - 2013.01.10 JM: Added setup lines for StarBench Tests.
/// - 2013.02.15 JM: Added NEW_METALLICITY flag for testing the new
///    microphysics classes.
/// - 2013.02.27 JM: Added IC_read_BBurkhart_data() class for
///    turbulent simulations.
/// - 2013.02.27 JM: Added class for Harpreet's 1D to 3D mapping.
/// - 2013.03.23 JM: Added setup lines for StarBench Tests.
/// - 2013.03.24 JM: Added another StarBench test.
/// - 2013.06.13 JM: Added StarBench_TremblinCooling test.
/// - 2013.08.23 JM: Added new mpv9_HHe module code.
/// - 2014.07.11 JM: Added isothermal noise perturbation option.
/// - 2015.01.(15-26) JM: Added new include statements for new PION version.
/// - 2015.02.03 JM: changed to use IC_base class MCMD pointer.
/// - 2016.05.02 JM: Changed order of code so that MP is initialised
///    before the "setup" function is called.
/// - 2017.08.03 JM: Don't write IC_filename.silo, just name it like
///    a normal snapshot.
/// - 2018.05.04 JM: moved parallel code to icgen_parallel.cpp

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
#include "dataIO/dataio_fits.h"
#endif // if FITS
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO

#include "grid/uniform_grid.h"
#include "grid/setup_fixed_grid.h"
#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;



// ##################################################################
// ##################################################################



int main(int argc, char **argv)
{
  class MCMDcontrol MCMD;

  if (argc<2) {
    cerr<<"Error, please give a filename to read IC parameters from.\n";
    cerr<<"Usage <icgen> <paramfile> [ic-filetype]\n";
    return(1);
  }

  int err=0;
  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) args[i] = argv[i];
  // Redirect stdout/stderr if required.
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      cout <<"Redirecting stdout to "<<outpath<<"info.txt"<<"\n";
      rep.redirect(outpath); // Redirects cout and cerr to text files in the directory specified.
    }
  }

  class DataIOBase   *dataio=0;
  class get_sim_info *siminfo=0;
  class ICsetup_base *ic=0;
  class ReadParams   *rp=0;
  class MCMDcontrol ppar; // unused for serial code.
  class SimParams SimPM;
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


  class setup_fixed_grid *SimSetup =0;
  SimSetup = new setup_fixed_grid();

  //
  // Set up the Xmin/Xmax/Range/dx of each level in the nested grid
  //
  //SimSetup->setup_nested_grid_levels(SimPM);

  vector<class GridBaseClass *> grid;

  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  int l=0;
  grid.resize(1);
  // Now set up the grid structure.
  cout <<"Init: &grid="<< &(grid[l])<<", and grid="<< grid[l] <<"\n";
  err = SimSetup->setup_grid(&(grid[l]),SimPM,&MCMD);
  cout <<"Init: &grid="<< &(grid[l])<<", and grid="<< grid[l] <<"\n";
  SimPM.dx = grid[0]->DX();
  if (!grid[0]) rep.error("Grid setup failed",grid[0]);
  
  //
  // read in what kind of ICs we are setting up.
  //
  rp = new ReadParams;
  if (!rp) rep.error("icgen:: initialising RP",rp);
  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);
  string seek="ics";
  string ics = rp->find_parameter(seek);
  setup_ics_type(ics,&ic);
  ic->set_SimPM(&SimPM);


  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    // don't need to set up the class, because it just does cooling and
    // there is no need to equilibrate anything.
  }
  else if (SimPM.ntracer>0 && (SimPM.EP.cooling || SimPM.EP.chemistry)) {
    cout <<"MAIN: setting up microphysics module\n";
    SimSetup->setup_microphysics(SimPM);
    if (!MP) rep.error("microphysics init",MP);
  }
  // ----------------------------------------------------------------

  //
  // Set up the boundary conditions and data
  //
  SimSetup->boundary_conditions(SimPM,ppar,grid[0]);
  if (err) rep.error("icgen Couldn't set up boundaries.",err);
  err += SimSetup->setup_raytracing(SimPM, grid[0]);
  err += SimSetup->setup_evolving_RT_sources(SimPM);
  err += SimSetup->update_evolving_RT_sources(SimPM,grid[0]->RT);
  if (err) rep.error("icgen: Failed to setup raytracer and/or microphysics",err);

  // ----------------------------------------------------------------
  // call "setup" to set up the data on the computational grid.
  err += ic->setup_data(rp,grid[0]);
  if (err) rep.error("Initial conditions setup failed.",err);
  // ----------------------------------------------------------------

  err = SimSetup->assign_boundary_data(SimPM,ppar,grid[0]);
  rep.errorTest("icgen::assign_boundary_data",0,err);

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
  string icfile   = rp->find_parameter(seek);
  if (icfile=="") {
    cout <<"WARNING: no filename for ic file.  writing to IC_temp.****\n";
    icfile = "IC_temp";
  }

  dataio = 0; // zero the class pointer.
  if (icftype=="text") {
    cout <<"WRITING ASCII TEXT FILE: ";
    cout << icfile << "\n";
  }

#ifdef FITS    
  if (icftype=="fits") {
    cout <<"WRITING FITS FILE: ";
    cout << icfile << "\n";
    dataio = 0; dataio = new DataIOFits (SimPM);
  }
#endif // if fits.

#ifdef SILO
  if (icftype=="silo") {
    cout <<"WRITING SILO FILE: ";
    //    icfile = icfile+".silo";
    cout <<icfile <<"\n";
    dataio=0; dataio=new dataio_silo (SimPM, "DOUBLE");
  }
#endif // if SILO defined.
  if (!dataio) rep.error("IO class initialisation: ",icftype);
  err = dataio->OutputData(icfile,grid, SimPM, 0);
  if (err) rep.error("File write error",err);
  delete dataio; dataio=0;
  cout <<icftype<<" FILE WRITTEN... exiting\n";

  // delete everything and return
  if (MP)   {delete MP; MP=0;}
  if (rp)   {delete rp; rp=0;} // Delete the read_parameters class.
  //if (grid) {??grid; grid=0;}
  if (ic)   {delete ic; ic=0;}
  if (SimSetup) {delete SimSetup; SimSetup =0;}

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

  delete [] args; args=0;
  return err;
}


// ##################################################################
// ##################################################################





