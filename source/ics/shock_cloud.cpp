/** \file shock_cloud.cc
 * 
 * File for setting up shock-cloud interaction problems.
 * */
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen.h"
#include <sstream>
using namespace std;

IC_shock_cloud::IC_shock_cloud() 
{
  preshock = postshock = 0;
  gg = 0; rp = 0;
  ndim = coords = eqns = -1;
  gam = dratio = pratio = Bratio = 0.0;
  clrad = cltr0 = 0.0;
  return;
}

IC_shock_cloud::~IC_shock_cloud()
{
  if (preshock)  {delete [] preshock;  preshock=0; }
  if (postshock) {delete [] postshock; postshock=0;}
  return;
}

int IC_shock_cloud::setup_data(class ReadParams *rrp,    ///< pointer to parameter list.
			     class GridBaseClass *ggg ///< pointer to grid
			     )
{
  int err=0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  IC_shock_cloud::ndim = SimPM->ndim;
  if (ndim!=2 && ndim!=3) rep.error("Shock-Cloud problem must be 2d or 3d",ndim);
  IC_shock_cloud::coords = SimPM->coord_sys;
  if (coords!=COORD_CRT && coords!=COORD_CYL) rep.error("Bad coord sys",coords);
  IC_shock_cloud::eqns = SimPM->eqntype;
  if      (eqns==EQEUL) eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD) eqns=2;
  else rep.error("Bad equations",eqns);

  // initialise pre and post shock vectors to zero.
  IC_shock_cloud::preshock = 0; IC_shock_cloud::postshock=0;
  preshock = new double [SimPM->nvar]; postshock = new double [SimPM->nvar];
  if (!preshock || !postshock) rep.error("malloc pre/post shock vecs",preshock);
  for (int v=0;v<SimPM->nvar;v++) preshock[v] = postshock[v] = 0.0;

  // set pre-/post- shock and cloud states.
  string seek, str;
  // cloud radius
  seek = "SCcloudradius";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_shock_cloud::clrad = atof(str.c_str());
  clrad *= SimPM->Range[YY];  // radius is given in units of y-dir range.

  seek = "SCcloudTR0";
  str = rp->find_parameter(seek);
  if (str=="") IC_shock_cloud::cltr0 = 0.0;
  else IC_shock_cloud::cltr0 = atof(str.c_str());
  seek = "SCcloudTR01";
  str = rp->find_parameter(seek);
  if (str=="") IC_shock_cloud::cltr1 = 0.0;
  else IC_shock_cloud::cltr1 = atof(str.c_str());

  // density ratio
  seek = "SCdratio";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_shock_cloud::dratio = atof(str.c_str());

  // pressure ratio
  seek = "SCpratio";
  str = rp->find_parameter(seek);
  if (str=="") IC_shock_cloud::pratio = 1.0; // default to pressure equilibrium
  else         IC_shock_cloud::pratio = atof(str.c_str());

  // B-field ratio
  seek = "SCBratio";
  str = rp->find_parameter(seek);
  if (str=="") IC_shock_cloud::Bratio = 1.0; // default to constant initial field
  else         IC_shock_cloud::Bratio = atof(str.c_str());

  // pre- and post-shock state vectors
  ostringstream temp;
  string v;
  v="RO";
  temp.str("");
  temp << "SCprevec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::preshock[RO] = atof(str.c_str());
  else         IC_shock_cloud::preshock[RO] = -1.0e99;
  temp.str("");
  temp << "SCpostvec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::postshock[RO] = atof(str.c_str());
  else         IC_shock_cloud::postshock[RO] = -1.0e99;

  v="PG";
  temp.str("");
  temp << "SCprevec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::preshock[PG] = atof(str.c_str());
  else         IC_shock_cloud::preshock[PG] = -1.0e99;
  temp.str("");
  temp << "SCpostvec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::postshock[PG] = atof(str.c_str());
  else         IC_shock_cloud::postshock[PG] = -1.0e99;

  v="VX";
  temp.str("");
  temp << "SCprevec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::preshock[VX] = atof(str.c_str());
  else         IC_shock_cloud::preshock[VX] = -1.0e99;
  temp.str("");
  temp << "SCpostvec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::postshock[VX] = atof(str.c_str());
  else         IC_shock_cloud::postshock[VX] = -1.0e99;

  v="VY";
  temp.str("");
  temp << "SCprevec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::preshock[VY] = atof(str.c_str());
  else         IC_shock_cloud::preshock[VY] = -1.0e99;
  temp.str("");
  temp << "SCpostvec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::postshock[VY] = atof(str.c_str());
  else         IC_shock_cloud::postshock[VY] = -1.0e99;

  v="VZ";
  temp.str("");
  temp << "SCprevec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::preshock[VZ] = atof(str.c_str());
  else         IC_shock_cloud::preshock[VZ] = -1.0e99;
  temp.str("");
  temp << "SCpostvec" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_shock_cloud::postshock[VZ] = atof(str.c_str());
  else         IC_shock_cloud::postshock[VZ] = -1.0e99;

  if (eqns ==2) { // mhd sim
    v="BX";
    temp.str("");
    temp << "SCprevec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::preshock[BX] = atof(str.c_str());
    else         IC_shock_cloud::preshock[BX] = -1.0e99;
    temp.str("");
    temp << "SCpostvec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::postshock[BX] = atof(str.c_str());
    else         IC_shock_cloud::postshock[BX] = -1.0e99;
    
    v="BY";
    temp.str("");
    temp << "SCprevec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::preshock[BY] = atof(str.c_str());
    else         IC_shock_cloud::preshock[BY] = -1.0e99;
    temp.str("");
    temp << "SCpostvec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::postshock[BY] = atof(str.c_str());
    else         IC_shock_cloud::postshock[BY] = -1.0e99;

    v="BZ";
    temp.str("");
    temp << "SCprevec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::preshock[BZ] = atof(str.c_str());
    else         IC_shock_cloud::preshock[BZ] = -1.0e99;
    temp.str("");
    temp << "SCpostvec" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::postshock[BZ] = atof(str.c_str());
    else         IC_shock_cloud::postshock[BZ] = -1.0e99;

    if (SimPM->eqntype==EQGLM) preshock[SI] = postshock[SI] = 0.;
  } // if mhd vars

  // tracer variables
  for (int t=0; t<SimPM->ntracer; t++) {
    temp.str("");
    temp << "SCprevecTR" << t;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::preshock[t+SimPM->ftr] = atof(str.c_str());
    else         IC_shock_cloud::preshock[t+SimPM->ftr] = -1.0e99;
    temp.str("");
    temp << "SCpostvecTR" << t;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_shock_cloud::postshock[t+SimPM->ftr] = atof(str.c_str());
    else         IC_shock_cloud::postshock[t+SimPM->ftr] = -1.0e99;
  }

  IC_shock_cloud::gam = SimPM->gamma;

  // now make sure we are to do a shock-cloud sim.
  string ics = rp->find_parameter("ics");
  if (ics=="") rep.error("didn't get any ics to set up.",ics);
  else if (ics=="ShockCloud" || ics=="shockcloud" || ics=="sc") {
    cout <<"\t\tSetting up Shock-Cloud problem.\n";
  }
  else rep.error("Don't know what Initial Condition is!",ics);

  // set cloud centre
  IC_shock_cloud::shockpos = SimPM->Xmin[XX]+SimPM->Range[XX]/8.0;
  IC_shock_cloud::cloudcentre = 0;
  cloudcentre = new double [ndim];
  if (!cloudcentre) rep.error("cloud centre malloc",cloudcentre);
  cloudcentre[XX] = shockpos +2.0*clrad;;
  if (coords==COORD_CRT) {
    // cartesian coords, so centre cloud accordingly.
    cloudcentre[YY] = SimPM->Xmin[YY]+SimPM->Range[YY]/2.0;
    if (ndim>2) cloudcentre[ZZ] = SimPM->Xmin[ZZ]+SimPM->Range[ZZ]/2.0;
  }
  else if (coords==COORD_CYL) {
    // cylindrical (Axial symmetry), so centre cloud on axis.
    cloudcentre[YY] = 0.0;
    if (ndim>2) rep.error("3d cylindrical coords not done yet!",ndim);
  }

  // setup the shock cloud.
  err += setup_shockcloud();

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

int IC_shock_cloud::setup_shockcloud()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with a shock";
  cout <<" hitting a circular cloud with overdensity of "<<dratio<<".\n";
  rep.printVec("preshock ",preshock, SimPM->nvar);
  rep.printVec("postshock",postshock,SimPM->nvar);
  int nsub;
  if (ndim==2) nsub=100; else nsub=32;
  class inside_sphere stest(cloudcentre,clrad,SimPM->dx,nsub,ndim);
  
  cout <<"\t\tAssigning primitive vectors.\n";
  double vfrac=0.0;
  double dpos[ndim];
  class cell *cpt = gg->FirstPt();
  do {
    CI.get_dpos(cpt,dpos);
    // Set values of primitive variables.
    if(dpos[XX] < shockpos)
      for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = postshock[v];
    else
      for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = preshock[v];
     
    // This is where I set the state inside the blast radius.
    if( (vfrac=stest.volumeFraction(cpt)) >0) {
      //cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
      cpt->P[RO] = vfrac*(dratio*cpt->P[RO]) + (1.-vfrac)*cpt->P[RO];
      cpt->P[PG] = vfrac*(pratio*cpt->P[PG]) + (1.-vfrac)*cpt->P[PG];
      if (eqns==2)
	cpt->P[BX] = vfrac*(Bratio*cpt->P[BX]) + (1.-vfrac)*cpt->P[BX];
      if (SimPM->ntracer>0)
	cpt->P[SimPM->ftr] = cltr0; // may be zero.
      if (SimPM->ntracer>1)
	cpt->P[SimPM->ftr+1] = cltr1; // may be zero.
      if (SimPM->ntracer>2) {
	//rep.error("Write proper cloud tracer read-in routine!",SimPM->ntracer);
        //
        // must be using Harpreet's microphysics which has lots of tracers, and 
        // will set values itself.
        //
        for (int v=0;v<SimPM->ntracer;v++) cpt->P[SimPM->ftr+v] = -2.0;
      }
    }
  } while ( (cpt=gg->NextPt(cpt))!=NULL);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.  
  
  return(0);
}
  

