/// \file laser_ablation.cc
/// File for setting up laser ablation problems.
/// Jonathan Mackey
/// 2009-12-18
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen_base.h"
#include "icgen.h"


// ##################################################################
// ##################################################################




IC_laser_ablation::IC_laser_ablation() 
{
  gg = 0; rp = 0;
  vel0 = rho0 = Pressure0 = gam = BT0 = BX0 = 0.0;
  return;
}


// ##################################################################
// ##################################################################




IC_laser_ablation::~IC_laser_ablation()
{
  return;
}



// ##################################################################
// ##################################################################




int IC_laser_ablation::setup_data(class ReadParams *rrp,    ///< pointer to parameter list.
				   class GridBaseClass *ggg ///< pointer to grid
				   )
{
  int err=0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  string seek, str;

  IC_laser_ablation::eqns = SimPM->eqntype;
  if      (eqns==EQEUL) eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD) eqns=2;
  else rep.error("Bad equations",eqns);

  // Find initial velocity (cm/s)
  seek = "LASERAB_vel0";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_laser_ablation::vel0 = atof(str.c_str());

  // Find initial gas density(g/cm^3)
  seek = "LASERAB_rho0";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_laser_ablation::rho0 = atof(str.c_str());

  // inside/outside density ratio (>>1)
  seek = "LASERAB_Dratio";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_laser_ablation::Dratio = atof(str.c_str());

  // Find initial gas temperature (K)
  seek = "LASERAB_Pressure0";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_laser_ablation::Pressure0 = atof(str.c_str());

  // inside/outside pressure ratio (>>1)
  seek = "LASERAB_Pratio";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_laser_ablation::Pratio = atof(str.c_str());

  // Find initial longitudinal B-field (Gauss)
  seek = "LASERAB_BX0";
  str = rp->find_parameter(seek);
  if (str=="") IC_laser_ablation::BX0 = 0.0;
  IC_laser_ablation::BX0 = atof(str.c_str());

  // Find initial transverse B-field (Gauss)
  seek = "LASERAB_BT0";
  str = rp->find_parameter(seek);
  if (str=="") IC_laser_ablation::BT0 = 0.0;
  IC_laser_ablation::BT0 = atof(str.c_str());

#ifdef NEW_B_NORM
  // convert from CGS to internal units (no factors of 4pi)
  BX0 /= sqrt(4.0*M_PI);
  BT0 /= sqrt(4.0*M_PI);
#endif


  IC_laser_ablation::gam = SimPM->gamma;


  // now make sure we are to do a radiative shock sim...
  string ics = rp->find_parameter("ics");

  if (ics=="") rep.error("didn't get any ics to set up.",ics);
  else if (ics=="LaserAblationAxi") {
    cout <<"\t\tsetting up test problem: "<<ics<<endl;
    if (SimPM->ndim!=2 || SimPM->coord_sys!=COORD_CYL) 
      rep.error("Asked for 3D problem but not 3D grid!",SimPM->ndim);
    err += setup_LaserAblationAxi();
  }
  else if (ics=="LaserAblation3D") {
    cout <<"\t\tsetting up test problem: "<<ics<<endl;
    if (SimPM->ndim!=3) rep.error("Asked for 3D problem but not 3D grid!",SimPM->ndim);
    err += setup_LaserAblation3D();
  }
  else rep.error("Don't know what Initial Condition is!",ics);

  // Add noise to data?  Smooth data?
  double noise=0.0; int smooth=0;
  ics = rp->find_parameter("noise");
  if (ics!="") noise = atof(ics.c_str());
  else noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number",noise);
  if (noise>0) err+= AddNoise2Data(gg, *SimPM, 2,noise);

  ics = rp->find_parameter("smooth");
  if (ics!="") smooth = atoi(ics.c_str());
  else smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number",smooth);
  if (smooth>0)err+= SmoothData(smooth);

  return err;
}



// ##################################################################
// ##################################################################





int IC_laser_ablation::setup_LaserAblationAxi()
{

  cout <<"\t\tSetting up laser ablation problem with v0="<<vel0;
  cout <<", rho="<<rho0<<", p0="<<Pressure0<<" ...\n";

  if (!pconst.equalD(BT0,0.0))
    rep.error("can't have transverse field in axisymmetric model!",BT0);

  double p0 = Pressure0;
  double p1 = Pressure0/Pratio;
  double r0 = rho0;
  double r1 = rho0/Dratio;  


  class cell *c = gg->FirstPt();
  double pos[SimPM->ndim];
  do {
    CI.get_dpos(c,pos);
    if (pos[XX]<0.0025 && pos[YY]<0.04) {
      c->P[RO] = r0;
      c->P[PG] = p0;
      c->P[VX] = vel0;
      c->P[VY] = c->P[VZ] = 0.0;
      if (eqns==2) {
	c->P[BY] = BT0; c->P[BX] = BX0; c->P[BZ] = 0.0;
      }
      // tracers (shouldn't be any except passive tracers.)
      for (int i=0;i<SimPM->ntracer;i++) c->P[SimPM->ftr+i] = 1.0;
    }
    else if (pos[XX]<0.0025 && pos[YY]<0.06) {
      c->P[RO] = r0 +50.0*(r1-r0)*(pos[YY]-0.04);
      c->P[PG] = p0 +50.0*(p1-p0)*(pos[YY]-0.04);
      c->P[VX] = vel0;
      c->P[VY] = c->P[VZ] = 0.0;
      if (eqns==2) {
	c->P[BY] = BT0; c->P[BX] = BX0; c->P[BZ] = 0.0;
      }
      // tracers (shouldn't be any except passive tracers.)
      for (int i=0;i<SimPM->ntracer;i++) c->P[SimPM->ftr+i] = 1.0;
    }
    else {
      c->P[RO] = r1;
      c->P[PG] = p1;
      c->P[VX] = vel0;
      c->P[VY] = c->P[VZ] = 0.0;
      if (eqns==2) {
	c->P[BY] = BT0; c->P[BX] = BX0; c->P[BZ] = 0.0;
      }
      // tracers (shouldn't be any except passive tracers.)
      for (int i=0;i<SimPM->ntracer;i++) c->P[SimPM->ftr+i] = 0.0;
    }
     // done.
  } while ((c=gg->NextPt(c))!=0);

  //c = gg->FirstPt();
  //for (int v=0;v<40;v++) {
  //  cout <<"rho="<<c->P[RO]<<"\tp_g="<<c->P[PG]<<"\ty="<<c->pos[YY]<<endl;
  //  c=gg->NextPt(c,YP);
  //}

  return 0;
}


// ##################################################################
// ##################################################################




int IC_laser_ablation::setup_LaserAblation3D()
{
  return 0;
}



// ##################################################################
// ##################################################################





