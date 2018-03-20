/// \file photoevaporating_clump.cc
/// 
/// File for setting up photo-evaporating clump problems.
///
///  - 2010-01-12 JM: Modified setup_radialprofile() so that the
///     profile is the same as that in Iliev et al. 2009, MNRAS, 400,
///     1283.  This means r<r0 rho=rho0; r>r0 rho=rho0*(r0/r)^2 Also
///     set the profile centre to be at cloudcentre[] not source_pos[]
///     to allow off-centre sources.
///
///  - 2010-01-15 JM: Modified setup_radialprofile() again, so that
///     pressure is also changed with density to give a constant
///     temperature profile.
///
///  - 2010-01-26 JM: modified setup_pec2() so that the outer 5 per
///     cent of the clump has a softened density profile (with linear
///     interpolation).  Effectively the function has been re-written.
///
///  - 2010.07.23 JM: New RSP source position class interface.
///
/// - 2011.02.28 JM: Got rid of RSP class.
/// - 2011.05.06 JM: Added geometric grid distance calculation to
///    photoevap-radial.
/// - 2011.06.23 JM: Added PhotoEvap_powerlaw function.
/// - 2012.09.16 JM: Added problem with 1/r^2 cloud profile and a
///    top-hat clump.
/// - 2013.08.23 JM: Change rad_src_info.position[] to .pos[]
/// - 2014.07.11 JM: Changed noise from adiabatic to isothermal
///    perturbation.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.05 JM: Added pion_flt datatype.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "icgen.h"
#include <sstream>



// ##################################################################
// ##################################################################



IC_photoevaporatingclump::IC_photoevaporatingclump() 
{
  ambient = 0;
  gg = 0; rp = 0;
  ndim = coords = eqns = -1;
  gam = dratio = pratio = Bratio = 0.0;
  clrad = 0.0;
  return;
}



// ##################################################################
// ##################################################################



IC_photoevaporatingclump::~IC_photoevaporatingclump()
{
  if (ambient) {delete [] ambient; ambient=0;}
  return;
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_data(class ReadParams *rrp, ///< pointer to parameter list.
					 class GridBaseClass *ggg ///< pointer to grid
					 )
{
  int err=0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  IC_photoevaporatingclump::ndim = SimPM->ndim;
  IC_photoevaporatingclump::coords = SimPM->coord_sys;
  IC_photoevaporatingclump::eqns = SimPM->eqntype;
  if      (eqns==EQEUL) eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD) eqns=2;
  else rep.error("Bad equations",eqns);

  // initialise pre and post shock vectors to zero.
  IC_photoevaporatingclump::ambient = 0;
  ambient = new double [SimPM->nvar];
  if (!ambient) rep.error("malloc pre/post shock vecs",ambient);
  for (int v=0;v<SimPM->nvar;v++) ambient[v] = 0.0;

  // set ambient and cloud states.
  string seek, str;
  // cloud radius
  seek = "PECcloudradius";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photoevaporatingclump::clrad = atof(str.c_str());
  if (ndim>1) clrad *= SimPM->Range[YY];  // radius is given in units of y-dir range.
  else        clrad *= SimPM->Range[XX];
  cout <<"Cloud Radius in cm is "<<clrad<<endl;

  cltr = mem.myalloc(cltr,SimPM->ntracer);
  ostringstream temp;

  for (int v=0;v<SimPM->ntracer;v++) {
    temp.str(""); temp<<"PECcloudTR"<<v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str=="") cltr[v] = 0.0;
    else         cltr[v] = atof(str.c_str());
  }

  // density ratio
  seek = "PECdratio";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photoevaporatingclump::dratio = atof(str.c_str());

  // pressure ratio
  seek = "PECpratio";
  str = rp->find_parameter(seek);
  if (str=="") IC_photoevaporatingclump::pratio = 1.0; // default to pressure equilibrium
  else         IC_photoevaporatingclump::pratio = atof(str.c_str());

  // B-field ratio
  seek = "PECBratio";
  str = rp->find_parameter(seek);
  if (str=="") IC_photoevaporatingclump::Bratio = 1.0; // default to constant initial field
  else         IC_photoevaporatingclump::Bratio = atof(str.c_str());

  // ambient medium state vector
  string v;
  v="RO";
  temp.str("");
  temp << "PEC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::ambient[RO] = atof(str.c_str());
  else         IC_photoevaporatingclump::ambient[RO] = -1.0e99;

  v="PG";
  temp.str("");
  temp << "PEC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::ambient[PG] = atof(str.c_str());
  else         IC_photoevaporatingclump::ambient[PG] = -1.0e99;

  v="VX";
  temp.str("");
  temp << "PEC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::ambient[VX] = atof(str.c_str());
  else         IC_photoevaporatingclump::ambient[VX] = -1.0e99;

  v="VY";
  temp.str("");
  temp << "PEC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::ambient[VY] = atof(str.c_str());
  else         IC_photoevaporatingclump::ambient[VY] = -1.0e99;

  v="VZ";
  temp.str("");
  temp << "PEC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::ambient[VZ] = atof(str.c_str());
  else         IC_photoevaporatingclump::ambient[VZ] = -1.0e99;

  if (eqns ==2) { // mhd sim
    v="BX";
    temp.str("");
    temp << "PEC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photoevaporatingclump::ambient[BX] = atof(str.c_str());
    else         IC_photoevaporatingclump::ambient[BX] = -1.0e99;
    
    v="BY";
    temp.str("");
    temp << "PEC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photoevaporatingclump::ambient[BY] = atof(str.c_str());
    else         IC_photoevaporatingclump::ambient[BY] = -1.0e99;

    v="BZ";
    temp.str("");
    temp << "PEC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photoevaporatingclump::ambient[BZ] = atof(str.c_str());
    else         IC_photoevaporatingclump::ambient[BZ] = -1.0e99;

    if (SimPM->eqntype==EQGLM) ambient[SI] = 0.;
  } // if mhd vars

  // tracer variables
  for (int t=0; t<SimPM->ntracer; t++) {
    temp.str("");
    temp << "PEC_ambTR" << t;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photoevaporatingclump::ambient[t+SimPM->ftr] = atof(str.c_str());
    else         IC_photoevaporatingclump::ambient[t+SimPM->ftr] = -1.0e99;
  }

  IC_photoevaporatingclump::gam = SimPM->gamma;

  // now make sure we are to do a photo-evaporation sim.
  string ics = rp->find_parameter("ics");
  if (ics=="") rep.error("didn't get any ics to set up.",ics);
  else if (ics=="PhotoEvaporatingClump" || ics=="PEC") {
    ics="PE_CLUMP";
    cout <<"\t\tSetting up PhotoEvaporating Clump problem.\n";
  }
  else if (ics=="PhotoEvaporatingClump2" || ics=="PEC2") {
    ics="PE_CLUMP2";
    cout <<"\t\tSetting up PhotoEvaporating Clump problem.\n";
  }
  else if (ics=="PhotoEvap_radial") {
    ics="PE_SLOPE";
    cout <<"\t\tSetting up Photoevaporation into radial density profile.\n";
  }
  else if (ics=="PhotoEvap_powerlaw") {
    ics="PE_POWERLAW";
    cout <<"\t\tSetting up Photoevaporation into power law density profile in Zcyl (XX).\n";
  }
  else if (ics=="PhotoEvap_paralleltest") {
    ics="PE_PARALLELTEST";
    cout <<"\t\tSetting up Photoevaporation of parallel rays with varying density.\n";
  }
  else if (ics=="PhotoEvap_CloudClump") {
    ics="PE_CLOUD_CLUMP";
    cout <<"\t\tSetting up Photoevaporation into cloud with clump.\n";
  }
  else rep.error("Don't know what Initial Condition is!",ics);

  // set cloud centre
  IC_photoevaporatingclump::cloudcentre = 0;
  cloudcentre = new double [ndim];
  if (!cloudcentre) rep.error("cloud centre malloc",cloudcentre);
   
  seek = "PEC_xpos";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  cloudcentre[XX] = atof(str.c_str());
   if (ndim>1) {
     seek = "PEC_ypos";
     str = rp->find_parameter(seek);
     if (str=="") rep.error("didn't find parameter",seek);
     cloudcentre[YY] = atof(str.c_str());
   }
   if (ndim>2) {
     seek = "PEC_zpos";
     str = rp->find_parameter(seek);
     if (str=="") rep.error("didn't find parameter",seek);
     cloudcentre[ZZ] = atof(str.c_str());
   }
   
  // set cloud to be near front of grid.
//  cloudcentre[XX] = SimPM->Xmin[XX] +SimPM->Range[XX]/4.0;
//  if (coords==COORD_CRT) {
//    // cartesian coords, so centre cloud accordingly.
//    cloudcentre[YY] = SimPM->Xmin[YY]+SimPM->Range[YY]/2.0;
//    if (ndim>2) cloudcentre[ZZ] = SimPM->Xmin[ZZ]+SimPM->Range[ZZ]/2.0;
//  }
//  else if (coords==COORD_CYL) {
//    // cylindrical (Axial symmetry), so centre cloud on axis.
//    cloudcentre[YY] = 0.0;
//    if (ndim>2) rep.error("3d cylindrical coords not done yet!",ndim);
//  }

  //
  // see if the ambient medium has a radial profile and/or core-radius in it.
  //
  IC_photoevaporatingclump::radial_slope=0.0;
  temp.str("");
  temp << "PEC_radialslope";
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::radial_slope = atof(str.c_str());
  else         IC_photoevaporatingclump::radial_slope = 0.0;
  cout <<"radial_slope="<<radial_slope<<"  *****************************\n";
  
  temp.str("");
  temp << "PEC_core_radius";
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photoevaporatingclump::core_radius = atof(str.c_str());
  else         IC_photoevaporatingclump::core_radius = 1.0;
  if (ndim>1) core_radius *= SimPM->Range[YY];  // r is fraction of y-range.
  else        core_radius *= SimPM->Range[XX];
  cout <<"core_radius="<<core_radius<<" cm *********************\n";
  
  //
  // setup the problem.
  //
  if      (ics=="PE_CLUMP")  err += setup_pec();
  else if (ics=="PE_CLUMP2") err += setup_pec2();
  else if (ics=="PE_SLOPE")  err += setup_radialprofile();
  else if (ics=="PE_POWERLAW")  err += setup_powerlaw_density();
  else if (ics=="PE_PARALLELTEST") err += setup_paralleltest();
  else if (ics=="PE_CLOUD_CLUMP") err += setup_cloud_clump();
  else rep.error("Bad ics",ics);

  //
  // Add noise to data?  Smooth data?
  //
  int smooth=0; double noise=0.0;
  ics = rp->find_parameter("noise");
  if (ics!="") noise = atof(ics.c_str());
  else noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number",noise);
  if (noise>0) {
    cout <<"\t\tNOISE!!! Adding random adiabatic noise";
    cout <<" at fractional level = "<<noise<<endl;
    err+= AddNoise2Data(gg, *SimPM, 4,noise);
  }
  ics = rp->find_parameter("smooth");
  if (ics!="") smooth = atoi(ics.c_str());
  else smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number",smooth);
  if (smooth>0)err+= SmoothData(smooth);

  delete [] cloudcentre; cloudcentre=0;
  cltr = mem.myfree(cltr);
  ambient = mem.myfree(ambient);
  return err;
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_pec2()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" hitting TWO circular clouds with overdensity of "<<dratio<<".\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);
  // ******* OLD CODE FOR HARD-EDGED CLUMPS *******
  //  int nsub;
  //  if (ndim==2) nsub=100; else nsub=32;
  // cloudcentre[YY] -= 0.5*clrad;
  // cloudcentre[XX] -= clrad;
  // class inside_sphere stest1(cloudcentre,clrad,SimPM->dx,nsub,ndim);
  
  // cloudcentre[YY] += clrad;
  // cloudcentre[XX] += 2.0*clrad;
  // class inside_sphere stest2(cloudcentre,clrad,SimPM->dx,nsub,ndim);
  // ******* OLD CODE FOR HARD-EDGED CLUMPS *******
  
  double cloud1[ndim], cloud2[ndim];
  for (int v=0;v<ndim;v++) cloud1[v]=cloud2[v]=cloudcentre[v];
  //
  // Hard-coded for model 3 in Lim & Mellema (2003).
  //
  cloud1[XX] -= 1.0*clrad;
  cloud1[YY] -= 0.5*clrad;
  
  cloud2[XX] += 1.0*clrad;
  cloud2[YY] += 0.5*clrad;
   

  cout <<"\t\tAssigning primitive vectors.\n";
 // double vfrac=0.0;
  double cell_pos[ndim], d;
  double rho1 = ambient[RO]*dratio;
  double pg1  = ambient[PG]*pratio;
  //  double BX1  = ambient[BX]*pratio;
  double edge=0.95*clrad;

  class cell *cpt = gg->FirstPt();
  do {
    //
    // Set values of primitive variables.
    //
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];

    //
    // Get distance of point from centre of cloud 1.  Then if we are
    // in the inner 95 per cent of the cloud, set the state to be the
    // cloud state.  If in the outer 5 per cent, then assign a linear
    // interpolation.
    //
    CI.get_dpos(cpt,cell_pos);
    d = gg->distance(cell_pos,cloud1);
    if (d < edge) {
      cpt->P[RO] = rho1;
      cpt->P[PG] = pg1;
      //if (eqns==2)
      //cpt->P[BX] = Bratio*cpt->P[BX];
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = cltr[v];
    }
    else if (d>=edge && d<clrad) {
      cpt->P[RO] = rho1 + (ambient[RO]-rho1)*(d-edge)/(clrad-edge);
      cpt->P[PG] =  pg1 + (ambient[PG]- pg1)*(d-edge)/(clrad-edge);
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = cltr[v] + 
	  (ambient[SimPM->ftr+v]-cltr[v])*(d-edge)/(clrad-edge);
    }

    //
    // Get distance of point from centre of cloud 2.  Then if we are
    // in the inner 95 per cent of the cloud, set the state to be the
    // cloud state.  If in the outer 5 per cent, then assign a linear
    // interpolation.
    //
    CI.get_dpos(cpt,cell_pos);
    d = gg->distance(cell_pos,cloud2);
    if (d < edge) {
      cpt->P[RO] = rho1;
      cpt->P[PG] = pg1;
      //if (eqns==2)
      //cpt->P[BX] = Bratio*cpt->P[BX];
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = cltr[v];
    }
    else if (d>=edge && d<clrad) {
      cpt->P[RO] = rho1 + (ambient[RO]-rho1)*(d-edge)/(clrad-edge);
      cpt->P[PG] =  pg1 + (ambient[PG]- pg1)*(d-edge)/(clrad-edge);
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = cltr[v] + 
	  (ambient[SimPM->ftr+v]-cltr[v])*(d-edge)/(clrad-edge);
    }

    // ******* OLD CODE FOR HARD-EDGED CLUMPS *******
    // // This is where I set the state inside the cloud.
    // if( (vfrac=stest1.volumeFraction(cpt)) >0) {
    //   //cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
    //   cpt->P[RO] = vfrac*(dratio*cpt->P[RO]) + (1.-vfrac)*cpt->P[RO];
    //   cpt->P[PG] = vfrac*(pratio*cpt->P[PG]) + (1.-vfrac)*cpt->P[PG];
    //   if (eqns==2)
    // 	cpt->P[BX] = vfrac*(Bratio*cpt->P[BX]) + (1.-vfrac)*cpt->P[BX];
    //   for (int v=0; v<SimPM->ntracer; v++)
    // 	cpt->P[SimPM->ftr+v] = vfrac*cltr[v] + (1.-vfrac)*cpt->P[SimPM->ftr+v];
    // }
    // if( (vfrac=stest2.volumeFraction(cpt)) >0) {
    //   //cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
    //   cpt->P[RO] = vfrac*(dratio*cpt->P[RO]) + (1.-vfrac)*cpt->P[RO];
    //   cpt->P[PG] = vfrac*(pratio*cpt->P[PG]) + (1.-vfrac)*cpt->P[PG];
    //   if (eqns==2)
    // 	cpt->P[BX] = vfrac*(Bratio*cpt->P[BX]) + (1.-vfrac)*cpt->P[BX];
    //   for (int v=0; v<SimPM->ntracer; v++)
    // 	cpt->P[SimPM->ftr+v] = vfrac*cltr[v] + (1.-vfrac)*cpt->P[SimPM->ftr+v];
    // }
    // ******* OLD CODE FOR HARD-EDGED CLUMPS *******

  } while ( (cpt=gg->NextPt(cpt))!=0);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.
  return(0);
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_pec()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" hitting a circular cloud with overdensity of "<<dratio<<".\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);
  int nsub;
  if (ndim==2) nsub=100; else nsub=32;
  class inside_sphere stest(cloudcentre,clrad,SimPM->dx,nsub,ndim);
  
  cout <<"\t\tAssigning primitive vectors.\n";
  double vfrac=0.0;
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
     
    // This is where I set the state inside the cloud.
    if( (vfrac=stest.volumeFraction(cpt)) >0) {
      //cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
      cpt->P[RO] = vfrac*(dratio*cpt->P[RO]) + (1.-vfrac)*cpt->P[RO];
      cpt->P[PG] = vfrac*(pratio*cpt->P[PG]) + (1.-vfrac)*cpt->P[PG];
      if (eqns==2)
	cpt->P[BX] = vfrac*(Bratio*cpt->P[BX]) + (1.-vfrac)*cpt->P[BX];
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = vfrac*cltr[v] + (1.-vfrac)*cpt->P[SimPM->ftr+v];
    }
  } while ( (cpt=gg->NextPt(cpt))!=NULL);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.
  return(0);
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_powerlaw_density()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" propagating into an ambient medium with 1/r^3 profile in Zcyl.\n";
  cout <<"WARNING: THIS FUNCTION HAS BEEN TAKEN OVER FOR A SPECIFIC PROBLEM!!\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);
  
  cout <<"\t\tAssigning primitive vectors.\n";
  //double dist=0.0;

  // *************  HARDCODED FOR NEW PROBLEM *****************
  // *************  HARDCODED FOR NEW PROBLEM *****************
  //
  // source position:
  //
  double srcpos[SimPM->ndim];
  for (int v=0;v<SimPM->ndim;v++) {
    srcpos[v] = SimPM->RS.sources[0].pos[v];
    if (isnan(srcpos[v])) rep.error("Bad source position",srcpos[v]);
  }

  //
  // Implementing a power law density profile going as x^3, assuming the 
  // source is at the origin, and rho=0 at x=-4pc, and nH=40 per cc at x=-3pc.
  //
  double xoffset = 12.344e18;
  double x0 = 3.086e18;
  double rho0 = 9.352e-23;
  double dpos[ndim];
  class cell *cpt = gg->FirstPt();
  do {
    CI.get_dpos(cpt,dpos);
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
    cpt->P[RO] = rho0*exp(3.0*log((dpos[XX]+xoffset)/x0));
  } while ( (cpt=gg->NextPt(cpt))!=0);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.
  return(0);
  // *************  HARDCODED FOR NEW PROBLEM *****************
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_cloud_clump()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" propagating into an ambient medium with 1/r^"<<radial_slope<<" profile.\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);

  //
  // Set up the clump-locator class.
  //
  int nsub;
  if (ndim==2) nsub=100; else nsub=32;
  class inside_sphere stest(cloudcentre,clrad,SimPM->dx,nsub,ndim);
  double vfrac=0.0;
  
  cout <<"\t\tAssigning primitive vectors.\n";
  double dist=0.0;

  //
  // source position:
  //
  double srcpos[SimPM->ndim];
  //if (SimPM->RS.Nsources!=1) {
  //  rep.error("Bad number of sources",SimPM->RS.Nsources);
  //}
  for (int v=0;v<SimPM->ndim;v++) {
    srcpos[v] = SimPM->RS.sources[0].pos[v];
    if (isnan(srcpos[v])) rep.error("Bad source position",srcpos[v]);
  }
  double ISM_centre[SimPM->ndim];
  for (int v=0;v<SimPM->ndim;v++) {
    ISM_centre[v]=0.0;
  }
  rep.printVec("source position ", srcpos,      SimPM->ndim);
  rep.printVec("clump centre    ", cloudcentre, SimPM->ndim);
  rep.printVec("ISM Cloud centre", ISM_centre,  SimPM->ndim);

  double dpos[ndim];
  class cell *cpt = gg->FirstPt();
  do {
    CI.get_dpos(cpt,dpos);
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
    //
    // Now if we have a radial profile in the slope, we need to adjust rho
    // Note the radial profile is measured from cloudcentre, not neccessarily
    // centred on the source (srcpos).  cloudcentre is set by (PEC_xpos,PEC_ypos,PEC_zpos).
    //
    if (!pconst.equalD(radial_slope,0.0)) {
      dist = gg->distance_vertex2cell(ISM_centre,cpt);

      //
      // We use rho=rho0/(1+(r/r0)^n)
      // We also leave the pressure constant.
      //
      cpt->P[RO] /= (1.0 + pow(dist/core_radius,radial_slope));
    }

    //
    // Now add in the clump.
    // Here dratio is the actual clump density.
    //
    if( (vfrac=stest.volumeFraction(cpt)) >0) {
      //cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
      cpt->P[RO] = std::max(static_cast<pion_flt>(vfrac*(dratio)),cpt->P[RO]);
      cpt->P[PG] = vfrac*(pratio*cpt->P[PG]) + (1.-vfrac)*cpt->P[PG];
      if (eqns==2)
	cpt->P[BX] = vfrac*(Bratio*cpt->P[BX]) + (1.-vfrac)*cpt->P[BX];
      for (int v=0; v<SimPM->ntracer; v++)
	cpt->P[SimPM->ftr+v] = vfrac*cltr[v] + (1.-vfrac)*cpt->P[SimPM->ftr+v];
    }
    //cout <<dpos[0]<<"  "<<cpt->P[RO]<<"  "<<cpt->P[PG]<<"\n";
  } while ( (cpt=gg->NextPt(cpt))!=0);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.


  return 0;
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_radialprofile()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" propagating into an ambient medium with 1/r^"<<radial_slope<<" profile.\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);
  
  cout <<"\t\tAssigning primitive vectors.\n";
  double dist=0.0;

  //
  // source position:
  //
  double srcpos[SimPM->ndim];
  //if (SimPM->RS.Nsources!=1) {
  //  rep.error("Bad number of sources",SimPM->RS.Nsources);
  //}
  for (int v=0;v<SimPM->ndim;v++) {
    srcpos[v] = SimPM->RS.sources[0].pos[v];
    if (isnan(srcpos[v])) rep.error("Bad source position",srcpos[v]);
  }
  rep.printVec("source position",srcpos,      SimPM->ndim);
  rep.printVec("cloud centre   ",cloudcentre, SimPM->ndim);
  //
  // set scale radius of core. (TAKE CLUMP RADIUS AS THE SCALE RADIUS OF CORE!!!)
  //
  double r0 = clrad; 

  double dpos[ndim];
  class cell *cpt = gg->FirstPt();
  do {
    CI.get_dpos(cpt,dpos);
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
    //
    // Now if we have a radial profile in the slope, we need to adjust rho
    // Note the radial profile is measured from cloudcentre, not neccessarily
    // centred on the source (srcpos).  cloudcentre is set by (PEC_xpos,PEC_ypos,PEC_zpos).
    //
    if (!pconst.equalD(radial_slope,0.0)) {
      dist = gg->distance_vertex2cell(cloudcentre,cpt);
      //
      // Following the Iliev et al 2009 test 6, we use rho=rho0(r0/r)^n if r>r0
      // We also change the pressure so there is a constant temperature state.
      //
      if (dist>r0) {
	cpt->P[RO] *= exp(radial_slope*log(r0/dist));
	cpt->P[PG] *= exp(radial_slope*log(r0/dist));
      }
    }
    //#define SINEWAVE_TESTING
#ifdef SINEWAVE_TESTING
    cpt->P[RO] *= (1.0 + 
		   pow(cos(16*M_PI*dpos[XX]/SimPM->Range[XX]),2.0)*
		   pow(cos(16*M_PI*dpos[YY]/SimPM->Range[YY]),2.0));
#endif
  } while ( (cpt=gg->NextPt(cpt))!=0);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  ( (cpt=nextPt(cpt))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.
  return(0);
}



// ##################################################################
// ##################################################################



int IC_photoevaporatingclump::setup_paralleltest()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" propagating into an ambient medium with parallel rays\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);
  
  cout <<"\t\tAssigning primitive vectors.\n";
  double srcpos[SimPM->ndim];
  for (int v=0;v<SimPM->ndim;v++) {
    srcpos[v] = SimPM->RS.sources[0].pos[v];
    if (isnan(srcpos[v])) rep.error("Bad source position",srcpos[v]);
  }
  rep.printVec("srcpos",srcpos,SimPM->ndim);

  class cell *c = gg->FirstPt(), *tmp=0;
  do {
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) c->P[v] = ambient[v];
    // Now we have a radial profile in the slope, so we need to adjust rho
    if (ndim>1 && SimPM->NG[YY]>2 && (tmp=gg->NextPt(c,YN))!=0) {
      c->P[RO] = 1.1*tmp->P[RO];
      c->P[PG] = 1.1*tmp->P[PG];
    }
  } while ( (c=gg->NextPt(c))!=0);
  //  c = firstPt();
  //  do {cout <<"c.rho = "<<c->P[RO]<<endl;} while  ( (c=nextPt(c))!=NULL);
  cout <<"\t\tGot through data successfully.\n";
  // Data done.
  return(0);
}

