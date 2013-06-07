///
/// file:    make_cooling_tables.cpp
/// author:  Jonathan Mackey
/// date:    2013.05.10
///
/// Adapted from templates/HIIregion_KineticE.cpp
///
/// Description: This is a template file for analysing simulation 
/// data, modified to just set up a microphysics class and write
/// heating and cooling rates to a text file.
///
/// Modifications:
///
///  - 2010-01-28 JM: Calculate a bunch of things for the Lim &
///   Mellema (2003) simulations.
/// - 2010-02-03 JM: Made it work for multi-core use, so can run with
///  8, 16, 32 cores etc. (and the number of cores doesn't have to
///  match the number of cores used to write the data i.e. independent
///  of the number of quadmeshes in the silo file).
/// - 2012.05.18 JM: Adapted template for calculating kinetic energy
///    of gas in simulations of HII regions around runaway stars.
/// - 2012.10.17 JM: Added more statistics.
/// - 2012.11.06 JM: Added yet more statistics, and output
///    mass-function files to their own directory.
///
/// - 2013.03.10 JM: Modified for the StarBench Spitzer test problem.
/// - 2013.03.15 JM: Modified for the Off-centre HII region test.
/// - 2013.03.21 JM: Modified for the TTYB79 test.
/// - 2013.05.10 JM: Modified for microphysics testing.


#include <iostream>
#include <sstream>
#include <silo.h>
#include <fitsio.h>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <global.h>

#include "microphysics/microphysics_base.h"

#ifndef EXCLUDE_MPV1
#include "microphysics/microphysics.h"
#endif 

#ifndef EXCLUDE_HD_MODULE
#include "microphysics/microphysics_lowZ.h"
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

#include "raytracing/raytracer_SC.h"


#include "microphysics/mpv5_molecular.h"
#include "microphysics/mpv6_PureH.h"
#include "microphysics/mpv7_TwoTempIso.h"




///
/// Function to setup the microphysics class (just so I can calculate gas
/// temperature in a manner consistent with how it was done in the 
/// simulation).  Function copied from gridMethods.cc
///
int setup_microphysics();

int setup_raytracing(); // hardcoded for zeta-oph for now!

// ##################################################################
// ##################################################################




int main(int argc, char **argv)
{

  //*******************************************************************
  //*******************************************************************
  //
  // Get command-line arguments
  //
  if (argc!=4) {
    cout << "Use as follows:\n";
    cout << "MP-TEST <microphysics-string>";
//    cout << " <ionising-src> <FUV-src>";
//    cout << " <Tstar> <Lstar> <LFUV> ";
    cout << " \n";
    cout <<"******************************************\n";
    cout <<"microphysics-string:  MPv[3/4/5/6/7] which microphysics to set up.\n";
//    cout <<"ionising-src: 0=no point source, 1=point source.\n";
//    cout <<"FUV-src:      0=no point source, 1=point source.\n";
//    cout <<"Tstar:        Effective temperature (K) of ionising source.\n";
//    cout <<"Lstar:        Luminosity (erg/s) of ionsing source.\n";
//    cout <<"LFUV:         FUV photon luminosity (/s) of point source.\n";
    cout <<"nH:  number density of Hydrogen in gas (cm^{-3}).\n";
    cout <<"r:   distance from radiation source (parsecs).\n";
    rep.error("Bad number of args",argc);
  }

  string MPtype = argv[1];
  double nH = atof(argv[2]);
  double r  = atof(argv[3])*GS.parsec();
//  int ion_src = atoi(argv[2]);
//  int fuv_src = atoi(argv[3]);
//  double Tstar = atof(argv[4]);
//  double Lstar = atof(argv[5]);
//  double Lfuv  = atof(argv[6]);


  // L = 4pi R^2 sigma*T^4
//  double Rstar = sqrt(Lstar/(4.0*M_PI*GS.StefanBoltzmannConst()))/(Tstar*Tstar);
//  cout <<"L*="<<Lstar<<", T*="<<Tstar<<", R*="<<Rstar<<"\n";

  //*******************************************************************
  //*******************************************************************


  //*******************************************************************
  //*******************************************************************
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Setting up microphysics ---------------\n";

  SimPM.nvar    = 6;
  SimPM.ftr     = 5;
  SimPM.ntracer = 1;
  SimPM.trtype = MPtype+"__H1+___";
  cout <<"tracer-type = "<<SimPM.trtype<<"\n";

  SimPM.EP.dynamics = 1;
  SimPM.EP.raytracing = 1;
  SimPM.EP.cooling = 1;
  SimPM.EP.chemistry = 1;
  SimPM.EP.coll_ionisation = 1;
  SimPM.EP.rad_recombination = 1;
  SimPM.EP.phot_ionisation = 1;
  SimPM.EP.update_erg = 1;
  SimPM.EP.MP_timestep_limit = 1;
  SimPM.EP.MinTemperature = 1.0e1;
  SimPM.EP.MaxTemperature = 1.0e9;
#ifdef THERMAL_CONDUCTION
  SimPM.EP.thermal_conduction = 0;
#endif // THERMAL CONDUCTION
  SimPM.EP.Helium_MassFrac = 0.2857;
  SimPM.EP.Metal_MassFrac = 0.0142;

  SimPM.gamma = 5.0/3.0;
  SimPM.coord_sys = COORD_CRT;
  SimPM.ndim = 3;

  //
  // Now setup microphysics and raytracing classes
  //
  int err = 0;
  err += setup_raytracing();
  err += setup_microphysics();
  if (err) rep.error("Setup of microphysics and raytracing",err);

  cout <<"--------------- Finished Setting up MP   --------------\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Output file: if multiple files, we will append _xxx to the name.
  // Initialise file handle to empty string.  Will use it later to
  // label output file.
  //
  //string filehandle("");
  ostringstream to;
  to <<MPtype<<"_heating_cooling_nH"<<nH<<"_r"<<r/GS.parsec()<<".txt";
  string this_outfile = to.str();

  //
  // Text File: proc 0 writes all the data.
  // 
  ofstream outf;
  if (outf.is_open())
    rep.error("Output text file is already open!",1);
  outf.open(this_outfile.c_str());
  if (!outf.is_open())
    rep.error("Failed to open text file for writing",this_outfile);
  outf <<"#\n# writing: T/K  G-L(no-src,y=0.002)  (no-src,y=0.1)  (no-src,y=0.995)";
  outf <<"  (src-10pc,y=0.002)  (src-10pc,y=0.1)  (src-10pc,y=0.995)\n";
  outf <<"#\n";
  outf.setf( ios_base::scientific );
  outf.precision(6);




  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Starting integration ----\n";
  cout <<"-------------------------------------------------------\n";

  //
  // Here we want to set gas density to nH=1.0/cm3, and for a large
  // range of temperatures calculate the rate E-dot.
  // We can do this by setting up a state vector with the appropriate
  // pressure for a given temperature and ionisation state.  Then we
  // can integrate forward in time by 1 second, and difference the
  // final and initial pressures, divide by (gamma-1), and this gives
  // E-dot in erg/cm3/s.
  //
  double p_in[SimPM.nvar], p_out[SimPM.nvar];
  p_in[RO] = 2.338e-24*nH;    // n(H)=1/c.c.
  p_in[PG] = 1e-12;        // arbitrary;
  p_in[VX] = p_in[VY] = p_in[VZ] = 0.0;
  p_in[SimPM.ftr] = 1.0e-3;

  double delt = 1.0e5;   // time to integrate for.
  int nheat=1, nion=1; // number of FUV-heating, and ionising srcs.
  
  double sf=1.0; // scaling factor for column densities
  //double r = 3.086e19; // distance to source.
  double dr= 0.01*r;   // thickness of shell that is being ionised.

  //
  // FUV heating from star.
  //
  struct rt_source_data hs;
  hs.id = 0;
  hs.type = RT_SRC_SINGLE;
  hs.strength = 3.0e47;  // What i used for simulations.
  hs.Vshell = 4.0*M_PI*((r+dr)*(r+dr)*(r+dr)-r*r*r)/3.0; // 3.697e56
  hs.dS     = dr;   // 0.01 pc thickness of shell.
  hs.Column = 7.215e-5*sf;  // 2.338e-24 g/cm3 for 10pc, total gas column.
  hs.DelCol = 7.215e-8*sf;  // 0.01pc, total gas coumn.
  std::vector<struct rt_source_data> Vhs;
  Vhs.clear();
  Vhs.push_back(hs);
  
  //
  // ionising radiation from star
  //
  struct rt_source_data is;
  is.id = 0;
  is.type = RT_SRC_SINGLE;
  is.strength = 4.42134e37;  // less than zeta Oph b/c BB approx. is v. bad!
  is.Vshell = 4.0*M_PI*((r+dr)*(r+dr)*(r+dr)-r*r*r)/3.0;
  is.dS     = dr;
  is.Column = 7.215e-10*sf;  // 2.338e-24 g/cm3 for 10pc, with 1e-5 neutral-frac.
  is.DelCol = 7.215e-13*sf;  // 0.01pc, at 1e-5 neutral fraction.
  std::vector<struct rt_source_data> Vis;
  Vis.clear();
  Vis.push_back(is);

  //
  // Also set up dummy sources with enormous column densities,
  // so that there is no photoionisation or photoheating from
  // the star
  //
  struct rt_source_data Dhs;
  Dhs.id = 0;
  Dhs.type = RT_SRC_SINGLE;
  Dhs.strength = 3.0e47;  // What i used for simulations.
  Dhs.Vshell = 4.0*M_PI*((r+dr)*(r+dr)*(r+dr)-r*r*r)/3.0;
  Dhs.dS     = dr;
  Dhs.Column = 7.215e+2/sf;  // 1e7*2.338e-24 g/cm3 for 10pc, total gas column.
  Dhs.DelCol = 7.215e-1/sf;  // 0.01pc, total gas coumn.
  std::vector<struct rt_source_data> VDhs;
  VDhs.clear();
  VDhs.push_back(Dhs);

  struct rt_source_data Dis;
  Dis.id = 0;
  Dis.type = RT_SRC_SINGLE;
  Dis.strength = 4.42134e37;  // less than zeta Oph b/c BB approx. is v. bad!
  Dis.Vshell = 4.0*M_PI*((r+dr)*(r+dr)*(r+dr)-r*r*r)/3.0;
  Dis.dS     = dr;
  Dis.Column = 7.215e+2/sf;  // 1e7*2.338e-24 g/cm3 for 10pc, with 1.0 neutral-frac.
  Dis.DelCol = 7.215e-1/sf;  // 0.01pc, at 1.0 neutral fraction.
  std::vector<struct rt_source_data> VDis;
  VDis.clear();
  VDis.push_back(Dis);

  //
  // Find equilibrium ion fractions by integrating for a long time...
  //
  double y_eq_sh=0.02;
  double y_eq_ns=0.99;
  double dy=0.0;

  //
  // First for the shielded case.
  //
  p_in[SimPM.ftr] = y_eq_sh;
  double T=1000.0, tt=0.0;
  size_t count=0;
  MP->Set_Temp(p_in,T,SimPM.gamma);
  do {
    delt = 0.01*MP->timescales_RT(p_in, nheat, VDhs, nion, VDis, SimPM.gamma);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, VDhs, nion, VDis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    //cout <<"dt="<<delt<<"\n";
    dy = fabs(p_out[SimPM.ftr]-p_in[SimPM.ftr])/p_in[SimPM.ftr];
    dy = std::max(dy,fabs(p_out[PG]-p_in[PG])/p_in[PG]);
    for (int v=0;v<SimPM.nvar;v++) p_in[v]=p_out[v];
    count ++;
  } while (dy>1.0e-8 && count <200);
  cout <<"y(H0)="<<1.0-p_out[SimPM.ftr]<<";  ";
  rep.printVec("p_in", p_in, SimPM.nvar);
  y_eq_sh = p_in[SimPM.ftr];

  //
  // Then for the unshielded case.
  //
  delt=1.0e5;
  p_in[SimPM.ftr] = y_eq_ns;
  T=5000.0; tt=0.0;
  count =0;
  MP->Set_Temp(p_in,T,SimPM.gamma);
  do {
    delt = 0.01*MP->timescales_RT(p_in, nheat, Vhs, nion, Vis, SimPM.gamma);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, Vhs, nion, Vis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    //cout <<"dt="<<delt<<"\n";
    dy = fabs(p_out[SimPM.ftr]-p_in[SimPM.ftr])/p_in[SimPM.ftr];
    dy = std::max(dy,fabs(p_out[PG]-p_in[PG])/p_in[PG]);
    for (int v=0;v<SimPM.nvar;v++) p_in[v]=p_out[v];
    count ++;
  } while (dy>1.0e-8 && count <200);
  rep.printVec("p_in", p_in, SimPM.nvar);
  y_eq_ns = p_in[SimPM.ftr];




  T=10.0; tt=0.0;
  double dt_max=3.16e5;
  do {
    outf <<T<<"    ";
    //
    // first with dummy radiation sources (giant extinction).
    //
    // y(H+)=equilibrium value
    //
    p_in[SimPM.ftr] = y_eq_sh;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, VDhs, nion, VDis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, VDhs, nion, VDis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";

    //
    // y(H+)=0.1
    //
    p_in[SimPM.ftr] = 0.1;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, VDhs, nion, VDis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, VDhs, nion, VDis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";

    //
    // y(H+)=0.995
    //
    p_in[SimPM.ftr] = 0.995;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, VDhs, nion, VDis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, VDhs, nion, VDis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";

    //
    // first with real radiation sources (little extinction).
    //
    // y(H+)=0.002
    //
    p_in[SimPM.ftr] = 0.002;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, Vhs, nion, Vis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, Vhs, nion, Vis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";
    //cout <<"LOW: T="<<T<<"; dp/p="<<fabs(p_in[PG]-p_out[PG])/p_in[PG];
    //cout <<",  dy/y="<<fabs(p_in[5]-p_out[5])/p_in[5] <<"\n";
    //rep.printVec("p_in ",p_in ,6);
    //rep.printVec("p_out",p_out,6);

    //
    // y(H+)=0.1
    //
    p_in[SimPM.ftr] = 0.1;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, Vhs, nion, Vis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, Vhs, nion, Vis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";
    //cout <<"T="<<T<<"; ";
    //cout <<"MED: T="<<T<<"; dp/p="<<fabs(p_in[PG]-p_out[PG])/p_in[PG];
    //cout <<",  dy/y="<<fabs(p_in[5]-p_out[5])/p_in[5] <<"\n";
    //rep.printVec("p_in ",p_in ,6);
    //rep.printVec("           p_out",p_out,6);

    //
    // y(H+)=equilibrium value
    //
    p_in[SimPM.ftr] = y_eq_ns;
    MP->Set_Temp(p_in,T,SimPM.gamma);
    delt = 0.01*MP->timescales_RT(p_in, nheat, Vhs, nion, Vis, SimPM.gamma);
    delt = std::min(delt,dt_max);
    err += MP->TimeUpdateMP_RTnew(p_in, nheat, Vhs, nion, Vis,
                           p_out, delt, SimPM.gamma, 0, &tt);
    outf <<(p_out[PG]-p_in[PG])/(SimPM.gamma-1.0)/delt<<"  ";
    //cout <<"T="<<T<<"; ";
    //cout <<"HIH: T="<<T<<"; dp/p="<<fabs(p_in[PG]-p_out[PG])/p_in[PG];
    //cout <<",  dy/y="<<fabs(p_in[5]-p_out[5])/p_in[5] <<"\n";
    //rep.printVec("p_in ",p_in ,6);
    //rep.printVec("           p_out",p_out,6);


    outf <<"\n";
    T *= 1.05;
  } while (T<5.0e4);

 


  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Finished integration ----\n";
  cout <<"-------------------------------------------------------\n";
  cout <<"--------------- Clearing up and Exiting ---------------\n";

  outf.close();

  if(grid!=0) {
    cout << "\t Deleting Grid Data..." << endl;
    delete grid; grid=0;
  }
  if (MP)     {delete MP; MP=0;}
  if (RT)     {delete RT; RT=0;}

  return 0;
}
// -------------------------------------------------------------
// *************************************************************
// **************** END MAIN MAIN MAIN MAIN ********************
// *************************************************************
// -------------------------------------------------------------


// ##################################################################
// ##################################################################




// stolen from gridMethods.cc
int setup_microphysics()
{
  cout <<"------------------------------------------------------------\n";
  cout <<"----------------- MICROPHYSICS SETUP -----------------------\n";
  cout <<"------------------------------------------------------------\n";
  //
  // Setup Microphysics class, if needed.
  // First see if we want the only_cooling class (much simpler), and if
  // not then check for the one of the bigger microphysics classes.
  //
  if (SimPM.EP.cooling && !SimPM.EP.chemistry) {
    cout <<"\t******* Requested cooling but no chemistry... setting";
    cout <<" up mp_only_cooling() class, with timestep-limiting.\n";
    SimPM.EP.MP_timestep_limit = 1;
    MP = new mp_only_cooling(SimPM.nvar, &(SimPM.EP));
    if (!MP) rep.error("mp_only_cooling() init",MP);
  }
  else if (SimPM.EP.chemistry) {
    //    MP = 0;
    cout <<"TRTYPE: "<<SimPM.trtype<<"\n";
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

      MP = new mp_explicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype, &(SimPM.EP));
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
      MP = new mp_implicit_H(SimPM.nvar, SimPM.ntracer, SimPM.trtype);
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
  }
  else {
    cout <<"\t******** not doing microphysics.\n";
    MP=0;
  }

  cout <<"************************************************************\n";
  cout <<"***************** MICROPHYSICS SETUP ***********************\n";
  cout <<"************************************************************\n";
  return 0;
}


// ##################################################################
// ##################################################################

int setup_raytracing()
{

  SimPM.RS.Nsources=2;

  //
  // first the ionising source.
  //
  struct rad_src_info rs_temp;
  for (int v=0;v<3;v++)
    rs_temp.position[v] = 0.0;
  rs_temp.strength = 4.42134e37;
  rs_temp.Rstar = 3.8478;  // solar radii
  rs_temp.Tstar = 3.05e4;
  rs_temp.id = 0;
  rs_temp.type = RT_SRC_SINGLE;
  rs_temp.update = RT_UPDATE_EXPLICIT;
  rs_temp.at_infinity = 0;
  rs_temp.effect = RT_EFFECT_PION_MULTI;
  rs_temp.opacity_src = RT_OPACITY_MINUS;
  rs_temp.opacity_var = SimPM.ftr;
  rs_temp.EvoFile = "NOFILE";
  SimPM.RS.sources.push_back(rs_temp);

  //
  // Then the FUV source.
  //
  struct rad_src_info rs_t2;
  for (int v=0;v<3;v++)
    rs_t2.position[v] = 0.0;
  rs_t2.strength = 3.0e47;
  rs_t2.Rstar = 1.0;  // solar radii
  rs_t2.Tstar = 1.0;
  rs_t2.id = 1;
  rs_t2.type = RT_SRC_SINGLE;
  rs_t2.update = RT_UPDATE_EXPLICIT;
  rs_t2.at_infinity = 0;
  rs_t2.effect = RT_EFFECT_UV_HEATING;
  rs_t2.opacity_src = RT_OPACITY_TOTAL;
  rs_t2.opacity_var = SimPM.ftr;
  rs_t2.EvoFile = "NOFILE";
  SimPM.RS.sources.push_back(rs_t2);

  return 0;
}

// ##################################################################
// ##################################################################


