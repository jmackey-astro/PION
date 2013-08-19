/// \file icgen.cc
/// \brief Program to generate Initial Conditions for my Uniform Grid code.
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

#include "ics/icgen.h"
#include "ics/get_sim_info.h"
#include "dataIO/dataio.h"
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#include "grid/uniform_grid.h"


#include "microphysics/microphysics_base.h"

#ifndef EXCLUDE_MPV1
#include "microphysics/microphysics.h"
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

#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"

#ifdef HARPREETS_CODE_EXT
#ifndef EXCLUDE_HD_MODULE
#include "microphysics/microphysics_lowZ.h"
#include "contrib/HD_MetalFree.h"
#endif // EXCLUDE_HD_MODULE
#endif // HARPREETS_CODE_EXT



#include "dataIO/readparams.h"

#include <sstream>
using namespace std;

int equilibrate_MP(GridBaseClass *, MicroPhysicsBase *, ReadParams *);

int main(int argc, char **argv)
{
  /** \todo This routine needs to be parallelised: decomposing the 
   * domain, setting up a parallel grid, adding noise equally on
   * subdomains, writing the data.  This ignores the various setup
   * classes, who should check if they need alterations for parallel
   * grids on a case-by-case basis.
   */
#ifdef PARALLEL
  int err = COMM->init(&argc, &argv);
  if (err) rep.error("comms init error",err);
#endif

  if (argc<2) {
    cerr<<"Error, please give a filename to read IC parameters from.\n";
    cerr<<"Usage <icgen> <paramfile> [ic-filetype]\n";
    return(1);
  }

#ifdef PARALLEL
  string *args=0;
  args = new string [argc];
  for (int i=0;i<argc;i++) args[i] = argv[i];
  // Redirect stdout/stderr if required.
  for (int i=0;i<argc; i++) {
    if (args[i].find("redirect=") != string::npos) {
      string outpath = (args[i].substr(9));
      ostringstream path; path << outpath <<"_"<<mpiPM.myrank<<"_";
      outpath = path.str();
      if (mpiPM.myrank==0) {
        cout <<"Redirecting stdout to "<<outpath<<"info.txt"<<"\n";
      }
      rep.redirect(outpath); // Redirects cout and cerr to text files in the directory specified.
    }
  }
  //cout << "rank: " << mpiPM.myrank << " nproc: " << mpiPM.nproc << "\n";
#endif //PARALLEL
#ifdef SERIAL
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
#endif // SERIAL

  class DataIOBase   *dataio=0;
  class get_sim_info *siminfo=0;
  class ICsetup_base *ic=0;
  class ReadParams   *rp=0;

  string pfile = argv[1];
  string icftype;
  if (argc>2) icftype=argv[2];
  else icftype="fits"; // This is the default for now.
  
  siminfo=0; siminfo = new class get_sim_info ();
  if (!siminfo) rep.error("Sim Info class init error",siminfo);
  err = 0;
  err += siminfo->read_gridparams(pfile);
  if (err) rep.error("Read Grid Params Error",err);
  delete siminfo; siminfo=0;
 
  //
  // May need to setup extra data in each cell for ray-tracing optical
  // depths and/or viscosity variables.
  //
  int hc_flag = 0, dv_flag=0;
  if (SimPM.artviscosity==AV_LAPIDUS ||
      SimPM.eqntype==EQEUL_EINT) {
    // Need one var. for Div(v)
    dv_flag = 1;
  }
  if (SimPM.artviscosity==AV_HCORRECTION ||
      SimPM.artviscosity==AV_HCORR_FKJ98) {
    // need one var for each dimension here.
    hc_flag = SimPM.ndim;
  }
  CI.setup_extra_data(SimPM.RS, hc_flag, dv_flag);
 
#ifdef PARALLEL
  err  = mpiPM.decomposeDomain();
  if (err) rep.error("main: failed to decompose domain!",err);
  grid =0;
  //
  // Now set up the parallel uniform grid.
  //
  cout <<"(icgen::setup_grid) Setting up grid...\n";

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
#endif // if PARALLEL

#ifdef SERIAL
  // Now we have read in parameters from the file, so set up a grid
  grid = 0; // global grid pointer.
#ifdef GEOMETRIC_GRID
  if      (SimPM.coord_sys==COORD_CRT)
    grid = new UniformGrid (SimPM.ndim, SimPM.nvar, SimPM.eqntype,
          SimPM.Xmin, SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_CYL)
    grid = new uniform_grid_cyl (SimPM.ndim, SimPM.nvar,
         SimPM.eqntype, SimPM.Xmin,
         SimPM.Xmax, SimPM.NG);
  else if (SimPM.coord_sys==COORD_SPH)
    grid = new uniform_grid_sph (SimPM.ndim, SimPM.nvar,
         SimPM.eqntype, SimPM.Xmin,
         SimPM.Xmax, SimPM.NG);
  else 
    rep.error("Bad Geometry in setup_grid()",SimPM.coord_sys);
#else  // GEOMETRIC_GRID
  grid = new UniformGrid(SimPM.ndim,SimPM.nvar,SimPM.eqntype,
       SimPM.Xmin,SimPM.Xmax,SimPM.NG);
#endif // GEOMETRIC_GRID
#endif // SERIAL

  if (!grid) rep.error("Grid setup failed",grid);
  

  // read in what kind of ICs we are setting up.
  rp = new ReadParams;
  if (!rp) rep.error("icgen:: initialising RP",rp);
  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);
  // Now we want to assign data to the grid, so we call whichever function is requested.
  string seek="ics";
  string ics = rp->find_parameter(seek);
  
  // invoke an appropriate class for whatever 'ics' is.
  if      (ics=="ShockTube")
    ic = new IC_shocktube ();
  // some basic tests...
  else if (ics=="OrszagTang" ||
     ics=="Uniform" ||
     ics=="Advection" || ics=="AdvectSineWave" ||
     ics=="KelvinHelmholz" || ics=="KelvinHelmholzStone" ||
     ics=="FieldLoop" || ics=="FieldLoopVz" || ics=="FieldLoopStatic"
     )
    ic = new IC_basic_tests();

  else if (ics=="Jet" || ics=="JET" || ics=="jet")
    ic = new IC_jet();

  else if (ics=="DoubleMachRef")
    ic = new IC_basic_tests();

  else if (ics=="RadiativeShock" ||
     ics=="RadiativeShockOutflow")
    ic = new IC_radiative_shock ();

  else if (ics=="LaserAblationAxi" ||
     ics=="LaserAblation3D")
    ic = new IC_laser_ablation ();

  //  else if (ics=="Jet")        ic = new IC_jet ();

  else if (ics=="ShockCloud")
    ic = new IC_shock_cloud ();

  else if (ics=="BlastWave" ||
           ics=="BlastWave_File")
    ic = new IC_blastwave ();

  else if (ics=="PhotoEvaporatingClump" ||
     ics=="PhotoEvaporatingClump2" ||
     ics=="PhotoEvap_radial" ||
     ics=="PhotoEvap_powerlaw" ||
     ics=="PhotoEvap_paralleltest" ||
     ics=="PhotoEvap_CloudClump")
    ic = new IC_photoevaporatingclump();

  else if (ics=="PhotEvap_RandomClumps" ||
     ics=="PERC" || ics=="PERC2" ||
     ics=="PhotEvap_RandomClumps2")
    ic = new IC_photevap_random_clumps();

  else if (ics=="PhotEvap_MultiClumps_FixNum" ||
     ics=="PE_MC_FN" || ics=="PE_MC_FM" ||
     ics=="PhotEvap_MultiClumps_FixMass")
    ic = new IC_photevap_multi_clumps();

  else if (ics=="Clump_Spherical" ||
           ics=="Clump_Axisymmetric")
    ic = new IC_spherical_clump();

  else if (ics=="StarBench_ContactDiscontinuity1" ||
           ics=="StarBench_ContactDiscontinuity2" ||
           ics=="StarBench_ContactDiscontinuity3" ||
           ics=="StarBench_ContactDiscontinuity4") {
    ic = new IC_StarBench_Tests();
  }

  else if (ics=="StarBench_IFI_testA" ||
           ics=="StarBench_IFI_testB" ||
           ics=="StarBench_IFI_testC") {
    ic = new IC_StarBench_Tests();
  }
  else if (ics=="StarBench_IrrCloud_Uniform" ||
           ics=="StarBench_IrrCloud_IsoSph") {
    ic = new IC_StarBench_Tests();
  }
  else if (ics=="StarBench_TremblinCooling") {
    ic = new IC_StarBench_Tests();
  }

#ifdef HARPREETS_CODE_EXT
#ifndef EXCLUDE_HD_MODULE
  else if (ics=="HD_2D_ShockCloud")
    ic = new IC_HD_2D_ShockCloud();
//  else if (ics=="HD_3D_ShockCloud")
//    ic = new IC_HD_3D_ShockCloud();
#endif
#endif // HARPREETS_CODE_EXT

#ifdef BBTURBULENCE_CODE_EXT
  else if (ics=="ReadBBTurbulence") {
    ic = new IC_read_BBurkhart_data();
  }
#endif // BBTURBULENCE_CODE_EXT

  else rep.error("BAD IC identifier",ics);
  if (!ic) rep.error("failed to init",ics);


  // call setup on the class just setup.
  err += ic->setup_data(rp,grid);
  if (err) rep.error("Initial conditions setup failed.",err);

  // if data initialised ok, see if we need to init microphysics variables,
  // and give them an equilibrium value.
  MP=0;  // global microphysics class pointer.
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout <<"\t******* Requested cooling but no chemistry... using";
    cout <<" up mp_only_cooling() class, with timestep-limiting.\n";
    SimPM.EP.MP_timestep_limit = 1;
    // don't need to set up the class, because it just does cooling and
    // there is no need to equilibrate anything.
  }
  else if (SimPM.ntracer>0 && (SimPM.EP.cooling || SimPM.EP.chemistry)) {
    cout <<"MAIN: setting up microphysics module\n";
    cout <<"TRTYPE: "<<SimPM.trtype<<"\n";
    // first avoid cooling the gas in getting to equilbrium, by
    // setting update_erg to false.
    bool uerg = SimPM.EP.update_erg;
    SimPM.EP.update_erg = false;
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
      MP = new mp_explicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP)
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
      MP = new mp_implicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
      //SimPM.EP.MP_timestep_limit = 4;  // limit by recombination time only
      //if (SimPM.EP.MP_timestep_limit <0 || SimPM.EP.MP_timestep_limit >5)
      //  rep.error("BAD dt LIMIT",SimPM.EP.MP_timestep_limit);
      have_set_MP=true;
    }
#endif // exclude MPv4


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


    err = equilibrate_MP(grid,MP,rp);
    if (err) rep.error("setting ionisation states to equilibrium failed",err);
    SimPM.EP.update_erg = uerg;
  }
  else {
    cout <<"MAIN: not doing ray-tracing or microphysics.\n";
  }


  //
  // Set up the boundary conditions, since internal boundary data
  // should really be already set to its correct value in the initial
  // conditions file.
  //
  grid->SetupBCs(2, SimPM.typeofbc);
  

  GS.start_timer("io"); double tsf=0;
  // MPI: WRITE PARALLEL FILES HERE
  // write data to file.
  cout <<"IC file-type is "<<icftype<<"\n";
  seek="ICfilename";
  string icfile   = rp->find_parameter(seek);
  string outfile;
  if (icfile=="") {
    cout <<"WARNING: no filename for ic file.  writing to IC_temp.****\n";
    icfile = "IC_temp";
  }

  dataio = 0; // zero the class pointer.

  if (icftype=="text") {
    cout <<"WRITING ASCII TEXT FILE: ";
    outfile=icfile;
    cout << outfile << "\n";
  }

#ifdef FITS    
  else if (icftype=="fits") {
    cout <<"WRITING FITS FILE: ";
    outfile=icfile;
    cout << outfile << "\n";
    dataio = 0; dataio = new DataIOFits ();
  }
#endif // if fits.

#ifdef SILO
  else if (icftype=="silo") {
    cout <<"WRITING SILO FILE: ";
    outfile=icfile;
#ifdef SERIAL
    //    outfile = icfile+".silo";
    cout <<outfile <<"\n";
    dataio=0; dataio=new dataio_silo ();
#endif
#ifdef PARALLEL
    cout <<outfile <<"\n";
    dataio=0; dataio=new dataio_silo_pllel ();
#endif
  }
#endif // if SILO defined.

  else rep.error("Don't recognise I/O type (text/fits/silo)",icftype);
  if (!dataio) rep.error("IO class initialisation: ",icftype);
  
  err = dataio->OutputData(outfile,grid, -1);
  if (err) rep.error("File write error",err);
  delete dataio; dataio=0;
  cout <<icftype<<" FILE WRITTEN in";
  tsf=GS.stop_timer("io");
  cout <<" "<<tsf<<" seconds.\n";

  // delete everything and return
  if (MP)   {delete MP; MP=0;}
  if (rp)   {delete rp; rp=0;} // Delete the read_parameters class.
  if (grid) {delete grid; grid=0;}
  if (ic)   {delete ic; ic=0;}

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

#ifdef PARALLEL
  cout << "rank: " << mpiPM.myrank << " nproc: " << mpiPM.nproc << "\n";
  COMM->finalise();
  delete COMM; COMM=0;
#endif
  delete [] args; args=0;
  return err;
}

int equilibrate_MP(
  class GridBaseClass *gg,
  class MicroPhysicsBase *mp,
  class ReadParams *rp
  )
{

  if (!mp || !gg || !rp) rep.error("microphysics or grid not initialised.",mp);
  rep.printVec("Init left  vec",gg->FirstPt()->P,SimPM.nvar);
  rep.printVec("Init right vec", gg->LastPt()->P,SimPM.nvar);

  string seek="InitIons";
  string s=rp->find_parameter(seek);
  if (s=="" || s=="YES" || s=="Y" || s=="y") {
    // integrate ion fractions to equilibrium
    cell *c = gg->FirstPt();
    do {
      //
      // Check if cell is boundary data or not (can only be an internal boundary, such as
      // a stellar wind, since we are looping over cells which are grid data).  If it is
      // boundary data, then we don't want to update anything, so we skip it
      //
      if (c->isbd) {
#ifdef TESTING
        cout <<"skipping cell "<<c->id<<" in equilibrate_MP() c->isbd.\n";
#endif
      }
      else {
        mp->Init_ionfractions(c->P,SimPM.gamma,-1);
      }
    } while ((c=gg->NextPt(c)) !=0);
    rep.printVec("init left  vec",gg->FirstPt()->P,SimPM.nvar);
    rep.printVec("init right vec", gg->LastPt()->P,SimPM.nvar);
    
    // now do a long time integration to get to equilibrium.
    c = gg->FirstPt(); double *p = c->P;
    double tint = sqrt(SimPM.gamma*p[PG]/p[RO]);
    tint = 50.*gg->DX()/tint; // gives us 50 times the dynamical time for a cell.
    cout <<"time to step by="<<tint<<"\n";
    double tt=0.;
    c = gg->FirstPt();
    do {
      if (!c->isbd) {
        for (int i=0;i<50;i++){
          mp->TimeUpdateMP(c->P, c->P, tint, SimPM.gamma, 0, &tt);
          //rep.printVec("Final left  vec",gg->FirstPt()->P,SimPM.nvar);
          cout <<"t="<<tt<<"\n";
        }
      }
    } while ((c=gg->NextPt(c)) !=0);
    rep.printVec("Final left  vec",gg->FirstPt()->P,SimPM.nvar);
    rep.printVec("Final right vec", gg->LastPt()->P,SimPM.nvar);
    c = gg->FirstPt();
    do {
      if (!c->isbd) {
        for (int i=0;i<50;i++){
          mp->TimeUpdateMP(c->P, c->P, tint, SimPM.gamma, 0, &tt);
        }
      }
    } while ((c=gg->NextPt(c)) !=0);
    rep.printVec("Final left  vec",gg->FirstPt()->P,SimPM.nvar);
    rep.printVec("Final right vec", gg->LastPt()->P,SimPM.nvar);
  }
  else if (s=="NO" || s=="N" || s=="n" || s=="no") {
    // initial values should be read from paramfile.
    string vb="Tracer"; cell *c=0;
    for (int i=0; i<SimPM.ntracer; i++) {
      ostringstream t; t<<vb<<i;
      string var=t.str(); cout <<"var: "<<var<<"\n";
      s=rp->find_parameter(var);
      if (s=="") rep.error("Don't know what to set tracer to.",s);
      else {
        double tr=atof(s.c_str());
        c = gg->FirstPt();
        do {c->P[SimPM.ftr+i] = tr;} while ((c=gg->NextPt(c)) !=0);
      }
    }
    rep.printVec("Final left  vec",gg->FirstPt()->P,SimPM.nvar);
    rep.printVec("Final right vec", gg->LastPt()->P,SimPM.nvar);
  }
  else if (s=="LEAVE") {
    // do nothing! hopefully a subroutine has set them already.
  }
  else rep.error("Bad InitIons specifier:",s);

  return 0;
}



int ICsetup_base::AddNoise2Data(int n, double frac)
{
  int seed= 975;
#ifdef PARALLEL
  seed += mpiPM.myrank;
  bool true_edge=false;
#endif
  srand(seed);
  class cell *cpt;  double avg=0.; long int ct=0;
  switch (n) {
    case 1:
    cout <<"\tAdding random noise to pressure at fractional level of "<<frac<<"\n";
    cpt=grid->FirstPt();
    do {
      if (!cpt->isbd) {
        avg += cpt->P[PG];
        ct++;
      }
    } while ( (cpt=grid->NextPt(cpt)) !=0);
#ifdef PARALLEL
    avg = COMM->global_operation_double("SUM",avg);
#endif
    cout <<"avg = "<<avg<< "\t ct= "<<ct;
    avg /= static_cast<double>(SimPM.Ncell);
    avg *= frac;  // avg is now a fraction 'frac' of the mean pressure on the grid.
    cout <<"avg = "<<avg<<"\n";
    cpt=grid->FirstPt();
    do {
#ifdef SERIAL
       if (!cpt->isedge && !cpt->isbd)
#endif 
#ifdef PARALLEL
      //
      // We want to exclude edge cells, but only those at the edge of the
      // full domain, not the local domain.
      //
      true_edge=false;
      if (cpt->isedge) {
        //
        // find out which direction the edge is (may be more than one!).
        // If any edge has no neighbour process, then it must be a full-domain
        // boundary.  If the neigbour doesn't exist, ngbproc(dir)<0.
        //
        //cout <<"Edge?: cpt="<<cpt->id<<", isedge="<<cpt->isedge;
        //cout <<"  true_edge="<<true_edge<<"\n";
        // x-dir
        if      (cpt->isedge %3 ==1) { // XN boundary
          //cout <<"got XN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
          if (mpiPM.ngbprocs[XN] <0) true_edge=true;
          //cout <<" ngb="<<mpiPM.ngbprocs[XN]<<", true_edge="<<true_edge<<"\n";
        }
        else if (cpt->isedge %3 ==2) { // XP boundary
          if (mpiPM.ngbprocs[XP] <0) true_edge=true;
        }
        // y-dir
        if (SimPM.ndim>1) {
          if      ((cpt->isedge%9)/3 ==1) { // YN boundary
            //cout <<"got YN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
            if (mpiPM.ngbprocs[YN] <0) true_edge=true;
            //cout <<" ngb="<<mpiPM.ngbprocs[YN]<<", true_edge="<<true_edge<<"\n";
          }
          else if ((cpt->isedge%9)/3 ==2) { // YP boundary
            if (mpiPM.ngbprocs[YP] <0) true_edge=true;
          }
        }
        // z-dir
        if (SimPM.ndim>2) {
          if      (cpt->isedge/9 ==1) { // ZN boundary
            if (mpiPM.ngbprocs[ZN] <0) true_edge=true;
          }
          else if (cpt->isedge/9 ==2) { // ZP boundary
            if (mpiPM.ngbprocs[ZP] <0) true_edge=true;
          }
        }
        //cout <<"true_edge="<<true_edge<<"\n";
      }
      if (true_edge==false && !cpt->isbd)
#endif
       {    // Don't want to alter any edge cells.
        //cout <<"PG before : "<<cpt->P[PG];
        cpt->P[PG] += avg*(static_cast<double>(rand())/RAND_MAX -0.5);
        //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
      }
    } while ( (cpt=grid->NextPt(cpt)) !=0);
    break;

    case 2:
    cout <<"Adding adiabatic random noise to cells at fractional level of "<<frac<<"\n";
    cpt=grid->FirstPt();
    do {
#ifdef SERIAL
       if (!cpt->isedge && !cpt->isbd)
#endif 
#ifdef PARALLEL
      //
      // We want to exclude edge cells, but only those at the edge of the
      // full domain, not the local domain.
      //
      true_edge=false;
      if (cpt->isedge) {
        //
        // find out which direction the edge is (may be more than one!).
        // If any edge has no neighbour process, then it must be a full-domain
        // boundary.  If the neigbour doesn't exist, ngbproc(dir)<0.
        //
        //cout <<"Edge?: cpt="<<cpt->id<<", isedge="<<cpt->isedge;
        //cout <<"  true_edge="<<true_edge<<"\n";
        // x-dir
        if      (cpt->isedge %3 ==1) { // XN boundary
          //cout <<"got XN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
          if (mpiPM.ngbprocs[XN] <0) true_edge=true;
          //cout <<" ngb="<<mpiPM.ngbprocs[XN]<<", true_edge="<<true_edge<<"\n";
        }
        else if (cpt->isedge %3 ==2) { // XP boundary
          if (mpiPM.ngbprocs[XP] <0) true_edge=true;
        }
        // y-dir
        if (SimPM.ndim>1) {
          if      ((cpt->isedge%9)/3 ==1) { // YN boundary
            //cout <<"got YN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
            if (mpiPM.ngbprocs[YN] <0) true_edge=true;
            //cout <<" ngb="<<mpiPM.ngbprocs[YN]<<", true_edge="<<true_edge<<"\n";
          }
          else if ((cpt->isedge%9)/3 ==2) { // YP boundary
            if (mpiPM.ngbprocs[YP] <0) true_edge=true;
          }
        }
        // z-dir
        if (SimPM.ndim>2) {
          if      (cpt->isedge/9 ==1) { // ZN boundary
            if (mpiPM.ngbprocs[ZN] <0) true_edge=true;
          }
          else if (cpt->isedge/9 ==2) { // ZP boundary
            if (mpiPM.ngbprocs[ZP] <0) true_edge=true;
          }
        }
        //cout <<"true_edge="<<true_edge<<"\n";
      }
      if (true_edge==false && !cpt->isbd)
#endif
       {    // Don't want to alter any edge cells.

         double temp;
         //if(cpt->isedge==0) {    // Don't want to alter any edge cells.
         //cout <<"PG before : "<<cpt->P[PG];
         //
         // final pressure = initial pressure * (1 +/- frac)
         // rho = pressure^(1/gamma)
         //
         temp = 2.*frac*(static_cast<double>(rand())/RAND_MAX -0.5);
         cpt->P[PG] *= 1+temp;
         cpt->P[RO] *= exp(log(1+temp)/SimPM.gamma);
         //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
         //}
      }
    } while ( (cpt=grid->NextPt(cpt)) !=0);
    break;
     
    case 3:
    cout <<"Adding adiabatic wave to cells at fractional level of "<<frac<<"\n";
#ifdef PARALLEL
    rep.error("adiabatic wave noise not implemented in parallel",3);
#endif
    // First get average pressure value.
    cpt=grid->FirstPt();
    do {if (!cpt->isedge && !cpt->isbd) avg += cpt->P[PG]; ct++;}
    while ( (cpt=grid->NextPt(cpt)) !=0);
    cout <<"avg = "<<avg<< "\t ct= "<<ct<<"\n";
    avg /= static_cast<double>(ct);
    // Now add wave to preshock state.
    cpt=grid->FirstPt();
    do {
       double temp;
       if(cpt->isedge==0 && !cpt->isbd && cpt->P[PG]<avg) {    // Don't want to alter any edge cells.
        //cout <<"PG before : "<<cpt->P[PG];
        temp = frac*sin(2.*M_PI*(CI.get_dpos(cpt,YY)/SimPM.Range[YY]) *(SimPM.NG[YY]/50));
         cpt->P[PG] *= 1+temp;
         cpt->P[RO] *= exp(log(1+temp)/SimPM.gamma);
         //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
       }
     } while ( (cpt=grid->NextPt(cpt)) !=0);
    break;
     
    default:
    cerr <<"\t Error, don't know what type of noise corresponds to "<<n<<"...\n";
    return(1);
  }
  return(0);
} // AddNoise2Data()

int ICsetup_base::SmoothData(int smooth) {return(0);}


