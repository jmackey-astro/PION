/// \file jet.cc
/// 
///  \author Jonathan Mackey
///
/// File for setting up jet simulations.
///
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen_base.h"
#include "ics/icgen.h"
#include <sstream>

IC_jet::IC_jet() 
{
  ambient = 0;
  gg = 0; rp = 0;
  ndim = coords = eqns = -1;
  jetrad=0;
  jdens = jpres = jvel = jtr0 = 0.0;
  JP.jetic = 1;
  return;
}

IC_jet::~IC_jet()
{
  if (ambient)  {delete [] ambient; ambient=0; }
  return;
}

int IC_jet::setup_data(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg ///< pointer to grid
      )
{
  int err=0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  IC_jet::ndim = SimPM->ndim;
  if (ndim!=2 && ndim!=3) rep.error("Shock-Cloud problem must be 2d or 3d",ndim);
  IC_jet::coords = SimPM->coord_sys;
  if (coords!=COORD_CRT && coords!=COORD_CYL) rep.error("Bad coord sys",coords);
  IC_jet::eqns = SimPM->eqntype;
  if      (eqns==EQEUL) eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD) eqns=2;
  else rep.error("Bad equations",eqns);

  //
  // initialise jet and ambient vectors to zero.
  //
  IC_jet::ambient = 0;
  ambient = new double [SimPM->nvar];
  if (!ambient) rep.error("malloc pre/post shock vecs",ambient);
  for (int v=0;v<SimPM->nvar;v++) ambient[v] = 0.0;


  string seek, str;

  //------------------------------------------------------------------
  //
  // First get the jet properties:
  //
  // radius
  //
  seek = "JETradius";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::jetrad = atoi(str.c_str());

  //
  // density in jet
  //
  seek = "JETdensity";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::jdens = atof(str.c_str());

  //
  // pressure in jet
  //
  seek = "JETpressure";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::jpres = atof(str.c_str());

  //
  // Jet velocity
  //
  seek = "JETvelocity";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::jvel  = atof(str.c_str());

  //
  // Jet Axial B-field
  //
  seek = "JET_Bax";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::j_bax  = atof(str.c_str());

  //
  // Jet Toroidal B-field
  //
  seek = "JET_Btor";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_jet::j_btor  = atof(str.c_str());

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    j_btor /= sqrt(4.0*M_PI);
    j_bax /= sqrt(4.0*M_PI);
#endif


  //------------------------------------------------------------------

  //------------------------------------------------------------------
  //
  // Ambient state vector
  //
  ostringstream temp;
  string v;
  v="RO";
  temp.str("");
  temp << "JETamb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_jet::ambient[RO] = atof(str.c_str());
  else         IC_jet::ambient[RO] = -1.0e99;

  v="PG";
  temp.str("");
  temp << "JETamb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_jet::ambient[PG] = atof(str.c_str());
  else         IC_jet::ambient[PG] = -1.0e99;

  v="VX";
  temp.str("");
  temp << "JETamb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_jet::ambient[VX] = atof(str.c_str());
  else         IC_jet::ambient[VX] = -1.0e99;

  v="VY";
  temp.str("");
  temp << "JETamb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_jet::ambient[VY] = atof(str.c_str());
  else         IC_jet::ambient[VY] = -1.0e99;

  v="VZ";
  temp.str("");
  temp << "JETamb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_jet::ambient[VZ] = atof(str.c_str());
  else         IC_jet::ambient[VZ] = -1.0e99;

  if (eqns ==2) { // mhd sim
    v="BX";
    temp.str("");
    temp << "JETamb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_jet::ambient[BX] = atof(str.c_str());
    else         IC_jet::ambient[BX] = -1.0e99;
    
    v="BY";
    temp.str("");
    temp << "JETamb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_jet::ambient[BY] = atof(str.c_str());
    else         IC_jet::ambient[BY] = -1.0e99;

    v="BZ";
    temp.str("");
    temp << "JETamb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_jet::ambient[BZ] = atof(str.c_str());
    else         IC_jet::ambient[BZ] = -1.0e99;

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    ambient[BX] /= sqrt(4.0*M_PI);
    ambient[BY] /= sqrt(4.0*M_PI);
    ambient[BZ] /= sqrt(4.0*M_PI);
#endif

    if (SimPM->eqntype==EQGLM) ambient[SI] = 0.;
  } // if mhd vars

  // tracer variables
  for (int t=0; t<SimPM->ntracer; t++) {
    temp.str("");
    temp << "JETambTR" << t;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_jet::ambient[t+SimPM->ftr] = atof(str.c_str());
    else         IC_jet::ambient[t+SimPM->ftr] = -1.0e99;
  }
  //------------------------------------------------------------------
  rep.printVec("JetState",JP.jetstate,SimPM->nvar);

  //
  // now make sure we are to do a jet sim.
  //
  string ics = rp->find_parameter("ics");
  if (ics=="") rep.error("didn't get any ics to set up.",ics);
  else if (ics=="jet" || ics=="Jet" || ics=="JET") {
    cout <<"\t\tSetting up Jet simulation.\n";
  }
  else rep.error("Don't know what Initial Condition is!",ics);

  //
  // Set all data to ambient state vector.
  // All of the jet parameters are in the global JP class/struct, and will
  // be used by the boundary conditions setup when the sim is starting.
  //
  class cell *c = gg->FirstPt();
  do {
     for (int v=0;v<SimPM->nvar;v++) c->P[v] = ambient[v];
  } while ( (c=gg->NextPt(c))!=0);

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
