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
/// - 2013.08.23 JM: Added new mpv9_HHe module code.
/// - 2014.07.11 JM: Added isothermal noise perturbation option.
/// - 2015.01.(15-26) JM: Added new include statements for new PION version.
/// - 2015.02.03 JM: changed to use IC_base class MCMD pointer.
/// - 2016.05.02 JM: Changed order of code so that MP is initialised
///    before the "setup" function is called.
/// - 2017.08.03 JM: Don't write IC_filename.silo, just name it like
///    a normal snapshot.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen.h"
#include "ics/get_sim_info.h"

#include "dataIO/dataio.h"
#include "dataIO/readparams.h"
#ifdef FITS
#include "dataIO/dataio_fits.h"
#include "dataIO/dataio_fits_MPI.h"
#endif // if FITS
#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO

#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"

#include "setup_fixed_grid.h"
#ifdef PARALLEL
#include "setup_fixed_grid_MPI.h"
#endif

#include "microphysics/microphysics_base.h"

#include <sstream>
using namespace std;

int equilibrate_MP(GridBaseClass *, microphysics_base *, ReadParams *, SimParams &);

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
#endif // PARALLEL

  class MCMDcontrol MCMD;
#ifdef PARALLEL
  int r=-1, np=-1;
  COMM->get_rank_nproc(&r,&np);
  MCMD.set_myrank(r);
  MCMD.set_nproc(np);
#endif // PARALLEL

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
      ostringstream path; path << outpath <<"_"<<MCMD.get_myrank()<<"_";
      outpath = path.str();
      if (MCMD.get_myrank()==0) {
        cout <<"Redirecting stdout to "<<outpath<<"info.txt"<<"\n";
      }
      rep.redirect(outpath); // Redirects cout and cerr to text files in the directory specified.
    }
  }
#ifndef TESTING
  rep.kill_stdout_from_other_procs(0);
#endif
  //cout << "rank: " << MCMD.get_myrank() << " nproc: " << MCMD.get_nproc() << "\n";
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
  class SimParams SimPM;

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
#if   defined (PARALLEL)
  SimSetup = new setup_fixed_grid_pllel();
#elif defined (SERIAL)
  SimSetup = new setup_fixed_grid();
#else
#error "Define SERIAL or PARALLEL"
#endif


  class GridBaseClass *grid = 0;
#ifdef PARALLEL
  err  = MCMD.decomposeDomain(SimPM);
  if (err) rep.error("main: failed to decompose domain!",err);
#endif // if PARALLEL
  //
  // Now we have read in parameters from the file, so set up a grid.
  //
  SimSetup->setup_grid(&grid,SimPM,&MCMD);
  if (!grid) rep.error("Grid setup failed",grid);
  
  //
  // read in what kind of ICs we are setting up.
  //
  rp = new ReadParams;
  if (!rp) rep.error("icgen:: initialising RP",rp);
  err += rp->read_paramfile(pfile);
  if (err) rep.error("Error reading parameterfile", pfile);
  //
  // Now we want to assign data to the grid, so we call whichever
  // function is requested.
  //
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

  else if (ics=="Jet" || ics=="JET" || ics=="jet") {
    ic = new IC_jet();
    ic->set_SimPM(&SimPM);
    err += ic->setup_data(rp,grid);
  }

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

#ifdef CODE_EXT_SBII
  else if (ics=="StarBench_ContactDiscontinuity1" ||
           ics=="StarBench_ContactDiscontinuity2" ||
           ics=="StarBench_ContactDiscontinuity3" ||
           ics=="StarBench_ContactDiscontinuity4" ||
           ics=="StarBench_IFI_testA"             ||
           ics=="StarBench_IFI_testB"             ||
           ics=="StarBench_IFI_testC"             ||
           ics=="StarBench_IFI_V2"                ||
           ics=="StarBench_IrrCloud_Uniform"      ||
           ics=="StarBench_IrrCloud_IsoSph"       ||
           ics=="StarBench_TremblinCooling"       ||
           ics=="StarBench_Cone") {
    ic = new IC_StarBench_Tests();
  }
#endif // CODE_EXT_SBII

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

  ic->set_SimPM(&SimPM);
#ifdef PARALLEL
  ic->set_MCMD_pointer(&MCMD);
#endif // PARALLEL

  // ----------------------------------------------------------------
  // We need to init microphysics class for some of the setups.
  //
  MP=0;  // global microphysics class pointer.

  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    // don't need to set up the class, because it just does cooling and
    // there is no need to equilibrate anything.
  }
  else if (SimPM.ntracer>0 && (SimPM.EP.cooling || SimPM.EP.chemistry)) {
    cout <<"MAIN: setting up microphysics module\n";

    cout <<"TRTYPE: \n";
    for (int i=0;i<SimPM.ntracer;i++) {
      cout <<"\t"<<i<<"\t"<<SimPM.tracers[i]<<"\n";
    }
    SimSetup->setup_microphysics(SimPM);
    if (!MP) rep.error("microphysics init",MP);
  }
  else {
    cout <<"MAIN: not doing ray-tracing or microphysics.\n";
  }
  // ----------------------------------------------------------------

  //
  // Set up the boundary conditions, since internal boundary data
  // should really be already set to its correct value in the initial
  // conditions file.
  //
  grid->SetupBCs(SimPM);
  if (err) rep.error("icgen Couldn't set up boundaries.",err);


  err += SimSetup->setup_raytracing(SimPM, grid);
  err += SimSetup->setup_evolving_RT_sources(SimPM);
  if (err) rep.error("icgen: Failed to setup raytracer and/or microphysics",err);

  // ----------------------------------------------------------------
  // call "setup" to set up the data on the computational grid.
  err += ic->setup_data(rp,grid);
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

    err = equilibrate_MP(grid,MP,rp,SimPM);
    if (err)
      rep.error("setting chemical states to equilibrium failed",err);

    SimPM.EP.update_erg = uerg;
    cout <<"MAIN: finished equilibrating the chemical species.\n";
  }
  // ----------------------------------------------------------------


  clk.start_timer("io"); double tsf=0;
  // MPI: WRITE PARALLEL FILES HERE
  // write data to file.
  cout <<"IC file-type is "<<icftype<<"\n";
  seek="OutputFile";
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
#ifdef SERIAL
    dataio = 0; dataio = new DataIOFits (SimPM);
#endif // SERIAL
#ifdef PARALLEL
    dataio = 0; dataio = new DataIOFits_pllel (SimPM, &MCMD);
#endif // PARALLEL
  }
#endif // if fits.

#ifdef SILO
  else if (icftype=="silo") {
    cout <<"WRITING SILO FILE: ";
    outfile=icfile;
#ifdef SERIAL
    //    outfile = icfile+".silo";
    cout <<outfile <<"\n";
    dataio=0; dataio=new dataio_silo (SimPM, "DOUBLE");
#endif
#ifdef PARALLEL
    cout <<outfile <<"\n";
    dataio=0; dataio=new dataio_silo_pllel (SimPM, "DOUBLE",&MCMD);
#endif
  }
#endif // if SILO defined.

  else rep.error("Don't recognise I/O type (text/fits/silo)",icftype);
  if (!dataio) rep.error("IO class initialisation: ",icftype);
  
  err = dataio->OutputData(outfile,grid, SimPM, 0);
  if (err) rep.error("File write error",err);
  delete dataio; dataio=0;
  cout <<icftype<<" FILE WRITTEN in";
  tsf=clk.stop_timer("io");
  cout <<" "<<tsf<<" seconds.\n";

  // delete everything and return
  if (MP)   {delete MP; MP=0;}
  if (rp)   {delete rp; rp=0;} // Delete the read_parameters class.
  if (grid) {delete grid; grid=0;}
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

#ifdef PARALLEL
  cout << "rank: " << MCMD.get_myrank() << " nproc: " << MCMD.get_nproc() << "\n";
  COMM->finalise();
  delete COMM; COMM=0;
#endif
  delete [] args; args=0;
  return err;
}





// ##################################################################
// ##################################################################




int equilibrate_MP(
      class GridBaseClass *gg,
      class microphysics_base *mp,
      class ReadParams *rp,
      class SimParams &SimPM
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
    c = gg->FirstPt(); pion_flt *p = c->P;
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




// ##################################################################
// ##################################################################




int ICsetup_base::AddNoise2Data(
      class GridBaseClass *grid,
      class SimParams &SimPM,
      int n,
      double frac
      )
{
  int seed= 975;
#ifdef PARALLEL
  seed += MCMD->get_myrank();
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
          if (MCMD->ngbprocs[XN] <0) true_edge=true;
          //cout <<" ngb="<<MCMD->ngbprocs[XN]<<", true_edge="<<true_edge<<"\n";
        }
        else if (cpt->isedge %3 ==2) { // XP boundary
          if (MCMD->ngbprocs[XP] <0) true_edge=true;
        }
        // y-dir
        if (SimPM.ndim>1) {
          if      ((cpt->isedge%9)/3 ==1) { // YN boundary
            //cout <<"got YN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
            if (MCMD->ngbprocs[YN] <0) true_edge=true;
            //cout <<" ngb="<<MCMD->ngbprocs[YN]<<", true_edge="<<true_edge<<"\n";
          }
          else if ((cpt->isedge%9)/3 ==2) { // YP boundary
            if (MCMD->ngbprocs[YP] <0) true_edge=true;
          }
        }
        // z-dir
        if (SimPM.ndim>2) {
          if      (cpt->isedge/9 ==1) { // ZN boundary
            if (MCMD->ngbprocs[ZN] <0) true_edge=true;
          }
          else if (cpt->isedge/9 ==2) { // ZP boundary
            if (MCMD->ngbprocs[ZP] <0) true_edge=true;
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
          if (MCMD->ngbprocs[XN] <0) true_edge=true;
          //cout <<" ngb="<<MCMD->ngbprocs[XN]<<", true_edge="<<true_edge<<"\n";
        }
        else if (cpt->isedge %3 ==2) { // XP boundary
          if (MCMD->ngbprocs[XP] <0) true_edge=true;
        }
        // y-dir
        if (SimPM.ndim>1) {
          if      ((cpt->isedge%9)/3 ==1) { // YN boundary
            //cout <<"got YN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
            if (MCMD->ngbprocs[YN] <0) true_edge=true;
            //cout <<" ngb="<<MCMD->ngbprocs[YN]<<", true_edge="<<true_edge<<"\n";
          }
          else if ((cpt->isedge%9)/3 ==2) { // YP boundary
            if (MCMD->ngbprocs[YP] <0) true_edge=true;
          }
        }
        // z-dir
        if (SimPM.ndim>2) {
          if      (cpt->isedge/9 ==1) { // ZN boundary
            if (MCMD->ngbprocs[ZN] <0) true_edge=true;
          }
          else if (cpt->isedge/9 ==2) { // ZP boundary
            if (MCMD->ngbprocs[ZP] <0) true_edge=true;
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

    //
    // Isothermal perturbation.
    //
   case 4:
    cout <<"Adding isothermal random noise to cells at fractional level of "<<frac<<"\n";
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
          if (MCMD->ngbprocs[XN] <0) true_edge=true;
          //cout <<" ngb="<<MCMD->ngbprocs[XN]<<", true_edge="<<true_edge<<"\n";
        }
        else if (cpt->isedge %3 ==2) { // XP boundary
          if (MCMD->ngbprocs[XP] <0) true_edge=true;
        }
        // y-dir
        if (SimPM.ndim>1) {
          if      ((cpt->isedge%9)/3 ==1) { // YN boundary
            //cout <<"got YN true boundary: cpt="<<cpt->id<<", isedge="<<cpt->isedge<<" ";
            if (MCMD->ngbprocs[YN] <0) true_edge=true;
            //cout <<" ngb="<<MCMD->ngbprocs[YN]<<", true_edge="<<true_edge<<"\n";
          }
          else if ((cpt->isedge%9)/3 ==2) { // YP boundary
            if (MCMD->ngbprocs[YP] <0) true_edge=true;
          }
        }
        // z-dir
        if (SimPM.ndim>2) {
          if      (cpt->isedge/9 ==1) { // ZN boundary
            if (MCMD->ngbprocs[ZN] <0) true_edge=true;
          }
          else if (cpt->isedge/9 ==2) { // ZP boundary
            if (MCMD->ngbprocs[ZP] <0) true_edge=true;
          }
        }
        //cout <<"true_edge="<<true_edge<<"\n";
      }
      if (true_edge==false && !cpt->isbd)
#endif
      {
        // Don't want to alter any edge cells.
        double temp;
        //
        // final pressure = initial pressure * (1 +/- frac)
        // final density  = initial density  * (1 +/- frac)
        //
        temp = 2.0*frac*(static_cast<double>(rand())/RAND_MAX -0.5);
        cpt->P[PG] *= 1.0+temp;
        cpt->P[RO] *= 1.0+temp;
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


