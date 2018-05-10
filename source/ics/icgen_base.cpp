/// \file icgen_base.cpp
/// \author Jonathan Mackey
/// \date 2018-05-04
///
/// Description:\n
/// This implements a set of routines that are common to serial,
/// parallel, and nested grid code for generating initial conditions.
///
/// Modifications:\n
/// 2018.05.04 JM: moved code from icgen.cpp/.h
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/timer.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "icgen_base.h"
#include "icgen.h"
#include "microphysics/microphysics_base.h"

#include <sstream>
using namespace std;




// ##################################################################
// ##################################################################


void setup_ics_type(
      string ics, ///< string giving type of ICs
      class ICsetup_base **ic ///< pointer to address of IC class
      )
{
  // invoke an appropriate class for whatever 'ics' is.
  if      (ics=="ShockTube")
    *ic = new IC_shocktube ();
  // some basic tests...
  else if (ics=="OrszagTang" ||
     ics=="Uniform" ||
     ics=="Advection" || ics=="AdvectSineWave" ||
     ics=="KelvinHelmholz" || ics=="KelvinHelmholzStone" ||
     ics=="FieldLoop" || ics=="FieldLoopVz" || ics=="FieldLoopStatic"
     )
    *ic = new IC_basic_tests();

  else if (ics=="Jet" || ics=="JET" || ics=="jet") {
    *ic = new IC_jet();
  }

  else if (ics=="DoubleMachRef")
    *ic = new IC_basic_tests();

  else if (ics=="RadiativeShock" ||
     ics=="RadiativeShockOutflow")
    *ic = new IC_radiative_shock ();

  else if (ics=="LaserAblationAxi" ||
     ics=="LaserAblation3D")
    *ic = new IC_laser_ablation ();

  else if (ics=="ShockCloud")
    *ic = new IC_shock_cloud ();

  else if (ics=="BlastWave" ||
           ics=="BlastWave_File")
    *ic = new IC_blastwave ();

  else if (ics=="PhotoEvaporatingClump" ||
     ics=="PhotoEvaporatingClump2" ||
     ics=="PhotoEvap_radial" ||
     ics=="PhotoEvap_powerlaw" ||
     ics=="PhotoEvap_paralleltest" ||
     ics=="PhotoEvap_CloudClump")
    *ic = new IC_photoevaporatingclump();

  else if (ics=="PhotEvap_RandomClumps" ||
     ics=="PERC" || ics=="PERC2" ||
     ics=="PhotEvap_RandomClumps2")
    *ic = new IC_photevap_random_clumps();

  else if (ics=="PhotEvap_MultiClumps_FixNum" ||
     ics=="PE_MC_FN" || ics=="PE_MC_FM" ||
     ics=="PhotEvap_MultiClumps_FixMass")
    *ic = new IC_photevap_multi_clumps();

  else if (ics=="Clump_Spherical" ||
           ics=="Clump_Axisymmetric")
    *ic = new IC_spherical_clump();

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
    *ic = new IC_StarBench_Tests();
  }
#endif // CODE_EXT_SBII

#ifdef HARPREETS_CODE_EXT
#ifndef EXCLUDE_HD_MODULE
  else if (ics=="HD_2D_ShockCloud")
    *ic = new IC_HD_2D_ShockCloud();
//  else if (ics=="HD_3D_ShockCloud")
//    *ic = new IC_HD_3D_ShockCloud();
#endif
#endif // HARPREETS_CODE_EXT

#ifdef BBTURBULENCE_CODE_EXT
  else if (ics=="ReadBBTurbulence") {
    *ic = new IC_read_BBurkhart_data();
  }
#endif // BBTURBULENCE_CODE_EXT

  else rep.error("BAD IC identifier",ics);
  if (!*ic) rep.error("failed to init",ics);

  return;
}



// ##################################################################
// ##################################################################



int ICsetup_base::equilibrate_MP(
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




// ##################################################################
// ##################################################################



int ICsetup_base::SmoothData(int smooth) {return(0);}


