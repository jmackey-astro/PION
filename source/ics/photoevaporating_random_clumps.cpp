/// \file photoevaporating_random_clumps.cc
/// 
/// File for setting up photo-evaporation of many random clumps.
///
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.02.03 JM: changed to use IC_base class MCMD pointer.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen_base.h"
#include "ics/icgen.h"
#include <sstream>



// ##################################################################
// ##################################################################



IC_photevap_random_clumps::IC_photevap_random_clumps() 
{
  ambient = 0;
  gg = 0; rp = 0;
  ndim = coords = eqns = Nclumps = profile = -1;
  gam = min_overdensity = max_overdensity = min_size = max_size = 0.0;
  radial_slope=0.0;
  cl = 0;
  return;
}



// ##################################################################
// ##################################################################


IC_photevap_random_clumps::~IC_photevap_random_clumps()
{
  if (ambient) {delete [] ambient; ambient=0;}
  if (cl) {delete [] cl; cl=0;}
  return;
}



// ##################################################################
// ##################################################################


int IC_photevap_random_clumps::setup_data(class ReadParams *rrp, ///< pointer to parameter list.
					 class GridBaseClass *ggg ///< pointer to grid
					 )
{
  int err=0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  IC_photevap_random_clumps::ndim = SimPM->ndim;
  if (ndim!=1 && ndim!=2 && ndim!=3) rep.error("Photoevaporation problem must be 1-3D",ndim);
  IC_photevap_random_clumps::coords = SimPM->coord_sys;
  if (coords!=COORD_CRT) rep.error("Bad coord sys",coords);
  IC_photevap_random_clumps::eqns = SimPM->eqntype;
  if      (eqns==EQEUL) eqns=1;
  else if (eqns==EQMHD ||
	   eqns==EQGLM ||
	   eqns==EQFCD) eqns=2;
  else rep.error("Bad equations",eqns);

  // initialise ambient vectors to zero.
  IC_photevap_random_clumps::ambient = 0;
  ambient = new double [SimPM->nvar];
  if (!ambient) rep.error("malloc ambient vecs",ambient);
  for (int v=0;v<SimPM->nvar;v++) ambient[v] = 0.0;

  // set ambient state, and clump parameters.
  string seek, str;
  // Number of Clumps
  seek = "PERCnumclumps";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::Nclumps = atoi(str.c_str());

  // profile of clumps
  seek = "PERCclump_profile";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::profile = atoi(str.c_str());

  cltr = mem.myalloc(cltr,SimPM->ntracer);
  ostringstream temp;

  for (int v=0;v<SimPM->ntracer;v++) {
    temp.str(""); temp<<"PERCcloudTR"<<v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str=="") cltr[v] = 0.0;
    else         cltr[v] = atof(str.c_str());
    cout <<"***********tracer["<<v<<"] = "<<cltr[v]<<endl;
  }

  seek = "PERCmin_size";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::min_size = atof(str.c_str());

  seek = "PERCmax_size";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::max_size = atof(str.c_str());


  // ambient medium state vector
  string v;
  v="RO";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::ambient[RO] = atof(str.c_str());
  else         IC_photevap_random_clumps::ambient[RO] = -1.0e99;

  // set this for now to the ambient density, and will take mass out of it later if needed.
  IC_photevap_random_clumps::ambdens = ambient[RO];

  v="PG";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::ambient[PG] = atof(str.c_str());
  else         IC_photevap_random_clumps::ambient[PG] = -1.0e99;

  v="VX";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::ambient[VX] = atof(str.c_str());
  else         IC_photevap_random_clumps::ambient[VX] = -1.0e99;

  v="VY";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::ambient[VY] = atof(str.c_str());
  else         IC_photevap_random_clumps::ambient[VY] = -1.0e99;

  v="VZ";
  temp.str("");
  temp << "PERC_amb" << v;
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::ambient[VZ] = atof(str.c_str());
  else         IC_photevap_random_clumps::ambient[VZ] = -1.0e99;

  if (eqns ==2) { // mhd sim
    v="BX";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photevap_random_clumps::ambient[BX] = atof(str.c_str());
    else         IC_photevap_random_clumps::ambient[BX] = -1.0e99;
    
    v="BY";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photevap_random_clumps::ambient[BY] = atof(str.c_str());
    else         IC_photevap_random_clumps::ambient[BY] = -1.0e99;

    v="BZ";
    temp.str("");
    temp << "PERC_amb" << v;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photevap_random_clumps::ambient[BZ] = atof(str.c_str());
    else         IC_photevap_random_clumps::ambient[BZ] = -1.0e99;

    if (SimPM->eqntype==EQGLM) ambient[SI] = 0.;
  } // if mhd vars

  // tracer variables
  for (int t=0; t<SimPM->ntracer; t++) {
    temp.str("");
    temp << "PERC_ambTR" << t;
    seek = temp.str();
    str = rp->find_parameter(seek);
    if (str!="") IC_photevap_random_clumps::ambient[t+SimPM->ftr] = atof(str.c_str());
    else         IC_photevap_random_clumps::ambient[t+SimPM->ftr] = -1.0e99;
  }

  IC_photevap_random_clumps::gam = SimPM->gamma;

  // now make sure we are to do a photo-evaporation sim.
  seek="ics";
  string ics = rp->find_parameter(seek);
  if (ics=="") rep.error("didn't get any ics to set up.",ics);
  else if (ics=="PhotEvap_RandomClumps" || ics=="PERC") {
    ics="PERC";
    cout <<"\t\tSetting up PhotoEvaporation problem with Random Clumps.\n";
  }
  else if (ics=="PhotEvap_RandomClumps2" || ics=="PERC2") {
    ics="PERC2";
    cout <<"\t\tSetting up PhotoEvaporation problem with Random Clumps (FIXED TOTAL MASS!).\n";
  }
  else rep.error("Don't know what Initial Condition is!",ics);

  // see if the ambient medium has a radial profile in it.
  IC_photevap_random_clumps::radial_slope=0.0;
  temp.str("");
  temp << "PERCradialslope";
  seek = temp.str();
  str = rp->find_parameter(seek);
  if (str!="") IC_photevap_random_clumps::radial_slope = atof(str.c_str());
  else         IC_photevap_random_clumps::radial_slope = 0.0;
  cout <<"radial_slope="<<radial_slope<<"  *****************************\n";
  
  // setup the problem.
  if      (ics=="PERC")  err += setup_perc();
  else if (ics=="PERC2") err += setup_perc_fixedmass();
  else rep.error("Bad ics",ics);

  //
  // Set tracer values!
  //
  //cout <<"ntracer = "<<SimPM->ntracer<<"\t amb[tr1]="<<ambient[SimPM->ftr+1]<<"\tcl[tr1]="<<cltr[1];
  cell *c=gg->FirstPt();
  cout <<"\tamb-density="<<ambdens<<" clump_mass="<<clump_mass<<endl;
  do {
    for (int v=SimPM->ftr; v<SimPM->nvar;v++) {
      //cout <<"trval before="<<c->P[v];
      c->P[v] = (ambdens/c->P[RO])*ambient[v] +(1.0-ambdens/c->P[RO])*cltr[v-SimPM->ftr];
      //cout <<"\tafter="<<c->P[v]<<endl;
    }
  } while ( (c=gg->NextPt(c)) !=0);

  // Add noise to data?  Smooth data?
  int smooth=0; double noise=0.0;
  seek = "noise";
  ics = rp->find_parameter(seek);
  if (ics!="") noise = atof(ics.c_str());
  else noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number",noise);
  if (noise>0.0) err+= AddNoise2Data(gg, *SimPM, 2,noise);

  seek = "smooth";
  ics = rp->find_parameter(seek);
  if (ics!="") smooth = atoi(ics.c_str());
  else smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number",smooth);
  if (smooth>0)err+= SmoothData(smooth);

  cltr = mem.myfree(cltr);
  ambient = mem.myfree(ambient);
  return err;
}


// ##################################################################
// ##################################################################



int IC_photevap_random_clumps::setup_perc_fixedmass()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" hitting "<<Nclumps<<" random clumps with total mass of ";

  string seek, str;
  seek = "PERCclump_mass";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::clump_mass = atof(str.c_str());

  seek = "PERCrandomseed";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  srand(atoi(str.c_str()));

  // clump_mass is as a fraction of the ambient mass.
  class cell *cpt = gg->FirstPt();
  double vol = gg->CellVolume(cpt);
  cout <<clump_mass*(ambient[RO]*vol*SimPM->Ncell)<<"\n\t\t(="<<clump_mass;
  cout <<" times ambient mass) and radial sizes of ";
  cout <<min_size<<" to "<<max_size<<" in units of y-range.\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);

  cout <<"\t\tAssigning primitive vectors.\n";
  double m1=0.0;
  do {
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
    cpt->P[RO] *= (1.0-clump_mass);
    m1 += cpt->P[RO] *vol;
  } while ( (cpt=gg->NextPt(cpt))!=0);
#ifdef PARALLEL
  double m2 = COMM->global_operation_double("SUM", m1);
  m1 = m2;
#endif // PARALLEL

  //
  // set ambient density according to how much mass we took out for clumps.
  //
  IC_photevap_random_clumps::ambdens *= (1.0-clump_mass);
  //
  // Convert clump mass to an actual total mass in clumps.
  //
  clump_mass *= ambient[RO]*vol*SimPM->Ncell;

  cout <<"setting up clump properties...\n";
#ifdef SERIAL
  int err=clumps_random_setup_fixedmass();
#endif //SERIAL
#ifdef PARALLEL
  int err=clumps_random_setup_pllel_fixedmass();
#endif // PARALLEL
  cout <<"adding clumps to grid...\n";
  cpt = gg->FirstPt();
  double mass=0.0;
  do {
    err += clumps_random_set_dens(cpt);
    mass += cpt->P[RO] *vol;
  } while ( (cpt=gg->NextPt(cpt))!=0);
#ifdef PARALLEL
  double mtot = COMM->global_operation_double("SUM", mass);
  mass = mtot;
#endif // PARALLEL
  cout <<"***** TOTAL MASS:"<<mass<<" CLUMP   MASS IS "<<clump_mass/mass<<" OF THIS.*****\n";
  cout <<"***** TOTAL MASS:"<<mass<<" AMBIENT MASS IS "<<        m1/mass<<" OF THIS.*****\n";
  cout <<"***** AMB+CLUMP:"<<clump_mass+m1<<"  *****\n";
  cout <<"***** if these don't match up, it's because some clump mass is off grid in XN/XP directions.\n";
  return(err);
}


// ##################################################################
// ##################################################################



int IC_photevap_random_clumps::setup_perc()
{
  cout <<"\t\tSetting up a "<<ndim<<"-D simulation with an I-front";
  cout <<" hitting "<<Nclumps<<" random clumps with overdensities of ";
  
  string seek, str;
  seek = "PERCmin_overdensity";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::min_overdensity = atof(str.c_str());

  seek = "PERCmax_overdensity";
  str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  IC_photevap_random_clumps::max_overdensity = atof(str.c_str());

  cout <<min_overdensity<<" to "<<max_overdensity<<" and radial sizes of";
  cout <<min_size<<" to "<<max_size<<" in units of y-range.\n";
  rep.printVec("ambient ",ambient, SimPM->nvar);

  cout <<"\t\tAssigning primitive vectors.\n";
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    for (int v=0;v<SimPM->nvar;v++) cpt->P[v] = ambient[v];
  } while ( (cpt=gg->NextPt(cpt))!=NULL);

  cout <<"setting up clump properties...\n";
#ifdef SERIAL
  int err=clumps_random_setup();
#endif //SERIAL
#ifdef PARALLEL
  int err=clumps_random_setup_pllel();
#endif // PARALLEL
  cout <<"adding clumps to grid...\n";
  cpt = gg->FirstPt();
  do {
    err += clumps_random_set_dens(cpt);
  } while ( (cpt=gg->NextPt(cpt))!=NULL);

  return(err);
}


// ##################################################################
// ##################################################################



double IC_photevap_random_clumps::random_frac()
{
  //
  // Returns a random value between 0 and 1.
  //
  return ((double) rand())/((double) RAND_MAX);
}


// ##################################################################
// ##################################################################



int IC_photevap_random_clumps::clumps_random_setup_fixedmass()
{
  int err=0;
  double xmax=SimPM->Range[XX];
  double ymax=SimPM->Range[YY];
  double zmax=SimPM->Range[ZZ];

  //
  // Assign random values within some sensible limits. 
  //
  double x3=max_size*ymax; // max. radius of clump (physical)
  double empty_xn=0.0, empty_xp=0.0;
  double empty_yn=0.0, empty_yp=0.0;
  double empty_zn=0.0, empty_zp=0.0;
  cout <<"x3 = "<<x3<<" and Xrange="<<SimPM->Range[XX]<<endl;

  string seek, str;
  seek = "PERCempty_xn";
  str = rp->find_parameter(seek);
  if (str=="") empty_xn = x3;
  else         empty_xn = atof(str.c_str())*xmax;
  seek = "PERCempty_xp";
  str = rp->find_parameter(seek);
  if (str=="") empty_xp = x3;
  else         empty_xp = atof(str.c_str())*xmax;
  seek = "PERCempty_yn";
  str = rp->find_parameter(seek);
  if (str=="") empty_yn = x3;
  else         empty_yn = atof(str.c_str())*ymax;
  seek = "PERCempty_yp";
  str = rp->find_parameter(seek);
  if (str=="") empty_yp = x3;
  else         empty_yp = atof(str.c_str())*ymax;
  seek = "PERCempty_zn";
  str = rp->find_parameter(seek);
  if (str=="") empty_zn = x3;
  else         empty_zn = atof(str.c_str())*zmax;
  seek = "PERCempty_zp";
  str = rp->find_parameter(seek);
  if (str=="") empty_zp = x3;
  else         empty_zp = atof(str.c_str())*zmax;

#define FIXED_MASS_RANGE_CLUMPS
  //#define FIXED_NUM_CLUMPS

#ifdef FIXED_MASS_RANGE_CLUMPS
  double min_mass=0.0, max_mass=0.0;
  // Get min and max clump masses, as fraction of total mass.
  seek = "PERCmin_mass";
  str = rp->find_parameter(seek);
  if (str=="") min_mass = 0.0;
  else         min_mass = atof(str.c_str())*clump_mass;
  seek = "PERCmax_mass";
  str = rp->find_parameter(seek);
  if (str=="") max_mass = 0.0;
  else         max_mass = atof(str.c_str())*clump_mass;
  cout <<"min/max masses: "<<min_mass<<", "<<max_mass<<",  and total mass: "<<clump_mass<<endl;

  std::vector<double> masses;
  double mass_remaining = clump_mass, this_mass=0.0;
  int nnn=0;
  do {
    this_mass = min_mass +random_frac()*(max_mass-min_mass);
    this_mass = std::min(this_mass, mass_remaining);
    masses.push_back(this_mass);
    mass_remaining -= this_mass;
    nnn++;
  } while (mass_remaining > (clump_mass*SMALLVALUE));
  cout <<"Used up all the mass in "<<nnn<<" clumps (Nclumps was "<<Nclumps<<")  Mass left over="<<mass_remaining<<"\n";
  Nclumps = nnn;
  // set up clumps
  cl = new struct clump [Nclumps];
  if (!cl) rep.error("Initialising clumps",cl);
  for (int j=0;j<Nclumps;j++) {
    cl[j].mass = masses[j];
  }
#endif // FIXED_MASS_RANGE_CLUMPS

#ifdef FIXED_NUM_CLUMPS
  // set up clumps
  cl = new struct clump [Nclumps];
  if (!cl) rep.error("Initialising clumps",cl);

  // choose random values on the interval [0,M_cl], and use them for clump masses
  std::vector<double> masses;
  masses.push_back(0.0);
  for (int j=1;j<Nclumps;j++) {
    masses.push_back(random_frac()*clump_mass);
  }
  masses.push_back(clump_mass);
  sort(masses.begin(),masses.end());
  for (int j=0;j<Nclumps+1;j++) cout <<"masses: "<<masses[j]<<endl;
  for (int j=0;j<Nclumps;j++) {
    cl[j].mass = masses[j+1]-masses[j];
  }
#endif // FIXED_NUM_CLUMPS

  for (int j=0;j<Nclumps;j++) {
    //
    // This restricts clump centres to be at least 0.1*Range from the edge.
    //
    //if (SimPM->typeofbc.find("YNper") != std::string::npos) {
    if (SimPM->BC_YN == "periodic") {
      //cout <<"find = "<<SimPM->typeofbc.find("YNper")<<endl;
      // for periodic BCs we don't care how near the y or z boundary a clump is, but
      // we do want it to be away from the incident radiation boundary.
      // max_size is a fraction of the y-range, so we have to modify it for x, to get
      // no clumps within a physical max_size of the XP boundary (so we don't lose 
      // lots of mass of the XP boundary).
      cl[j].centre[XX] = SimPM->Xmin[XX]+ empty_xn +(xmax-empty_xp-empty_xn)*random_frac();
      cl[j].centre[YY] = SimPM->Xmin[YY]+ ymax*random_frac();
      cl[j].centre[ZZ] = SimPM->Xmin[ZZ]+ zmax*random_frac();
    }
    else {
      //cout <<"not periodic bcs!\n";
      cl[j].centre[XX] = SimPM->Xmin[XX]+ empty_xn +(xmax-empty_xp-empty_xn)*random_frac();
      cl[j].centre[YY] = SimPM->Xmin[YY]+ empty_yn +(ymax-empty_yp-empty_yn)*random_frac();
      cl[j].centre[ZZ] = SimPM->Xmin[ZZ]+ empty_zn +(zmax-empty_zp-empty_zn)*random_frac();
    }
    // radius of clumps
    cl[j].size[XX] = min_size*ymax +(max_size-min_size)*ymax*random_frac();
    cl[j].size[YY] = min_size*ymax +(max_size-min_size)*ymax*random_frac();
    cl[j].size[ZZ] = min_size*ymax +(max_size-min_size)*ymax*random_frac();
    // orientation from grid axes to principal axes.
    cl[j].ang[XX] = random_frac()*M_PI;
    cl[j].ang[YY] = random_frac()*M_PI;
    cl[j].ang[ZZ] = random_frac()*M_PI;
    //
    // zero the z-dir if 2d
    if (SimPM->ndim==2) {
      cl[j].ang[ZZ] = cl[j].ang[YY] = cl[j].size[ZZ] = cl[j].centre[ZZ] = 0.0;
    }

    // set overdensity
    //
    cl[j].overdensity = cl[j].mass/ambient[RO];
    for (int v=0;v<SimPM->ndim;v++) cl[j].overdensity /= sqrt(M_PI)*cl[j].size[v];

    //
    // print clump...
    print_clump(&cl[j]);
    //
    // now set up the rotation matrix.

    cl[j].rm[XX][XX] =  cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[XX][YY] =  sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[XX][ZZ] = -sin(cl[j].ang[YY])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][XX] = -cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][YY] = -sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][ZZ] =  sin(cl[j].ang[YY])*sin(cl[j].ang[ZZ]);
    cl[j].rm[ZZ][XX] =  cos(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[ZZ][YY] =  sin(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[ZZ][ZZ] =  cos(cl[j].ang[YY]);
    // all done, move on to next clump.
  }
  return err;
}
 

// ##################################################################
// ##################################################################


   
int IC_photevap_random_clumps::clumps_random_setup()
{
  int err=0;

  string seek = "PERCrandomseed";
  string str = rp->find_parameter(seek);
  if (str=="") rep.error("didn't find parameter",seek);
  srand(atoi(str.c_str()));
  //  srand(12);
  //
  double xmax=SimPM->Range[XX];
  double ymax=SimPM->Range[YY];
  double zmax=SimPM->Range[ZZ];

  // set up clumps
  cl = new struct clump [Nclumps];
  if (!cl) rep.error("Initialising clumps",cl);

  //
  // Assign random values within some sensible limits. 
  //
  for (int j=0;j<Nclumps;j++) {
    cl[j].overdensity = min_overdensity +random_frac()*(max_overdensity-min_overdensity);
    //
    // This restricts clump centres to be at least 0.1*Range from the edge.
    //
    //if (SimPM->typeofbc.find("YNper") != std::string::npos) {
    if (SimPM->BC_YN == "periodic") {
      //cout <<"find = "<<SimPM->typeofbc.find("YNper")<<endl;
      // for periodic BCs we don't care how near the y or z boundary a clump is, but
      // we do want it to be away from the incident radiation boundary.
      // max_size is a fraction of the y-range, so we have to modify it for x, to get
      // no clumps within a physical max_size of the XP boundary (so we don't lose 
      // lots of mass of the XP boundary).
      cl[j].centre[XX] = 0.2*xmax +(0.8-(2.0*max_size*xmax/ymax))*xmax*random_frac();
      cl[j].centre[YY] = ymax*random_frac();
      cl[j].centre[ZZ] = zmax*random_frac();
    }
    else {
      //cout <<"not periodic bcs!\n";
      cl[j].centre[XX] = 0.1*xmax +0.8*xmax*random_frac();
      cl[j].centre[YY] = max_size*ymax +(1.0-2.0*max_size)*ymax*random_frac();
      cl[j].centre[ZZ] = max_size*zmax +(1.0-2.0*max_size)*zmax*random_frac();
    }
    // radius of clumps
    cl[j].size[XX] = min_size*ymax +max_size*ymax*random_frac();
    cl[j].size[YY] = min_size*ymax +max_size*ymax*random_frac();
    cl[j].size[ZZ] = min_size*ymax +max_size*ymax*random_frac();
    // orientation from grid axes to principal axes.
    cl[j].ang[XX] = random_frac()*M_PI;
    cl[j].ang[YY] = random_frac()*M_PI;
    cl[j].ang[ZZ] = random_frac()*M_PI;
    //
    // zero the z-dir if 2d
    if (SimPM->ndim==2) {
      cl[j].ang[ZZ] = cl[j].ang[YY] = cl[j].size[ZZ] = cl[j].centre[ZZ] = 0.0;
    }
    //
    // print clump...
    print_clump(&cl[j]);
    //
    // now set up the rotation matrix.

    cl[j].rm[XX][XX] =  cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[XX][YY] =  sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[XX][ZZ] = -sin(cl[j].ang[YY])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][XX] = -cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][YY] = -sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][ZZ] =  sin(cl[j].ang[YY])*sin(cl[j].ang[ZZ]);
    cl[j].rm[ZZ][XX] =  cos(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[ZZ][YY] =  sin(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[ZZ][ZZ] =  cos(cl[j].ang[YY]);

    //
    /*
    cl[j].rm[XX][XX] =  cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[XX][YY] =  cos(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) -sin(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[XX][ZZ] =  cos(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[YY][XX] =  sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*cos(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*sin(cl[j].ang[ZZ]);
    cl[j].rm[YY][YY] = -sin(cl[j].ang[XX])*cos(cl[j].ang[YY])*sin(cl[j].ang[ZZ]) +cos(cl[j].ang[XX])*cos(cl[j].ang[ZZ]);
    cl[j].rm[YY][ZZ] =  sin(cl[j].ang[XX])*sin(cl[j].ang[YY]);
    cl[j].rm[ZZ][XX] = -sin(cl[j].ang[YY])*cos(cl[j].ang[ZZ]);
    cl[j].rm[ZZ][YY] =  sin(cl[j].ang[YY])*sin(cl[j].ang[ZZ]);
    cl[j].rm[ZZ][ZZ] =  cos(cl[j].ang[YY]);
    */
    // all done, move on to next clump.
  }
  return err;
}


// ##################################################################
// ##################################################################



#ifdef PARALLEL
int IC_photevap_random_clumps::clumps_random_setup_pllel_fixedmass()
{
  int err=0;
  if (MCMD->get_myrank()==0) {
    err += IC_photevap_random_clumps::clumps_random_setup_fixedmass();
    err += COMM->broadcast_data(0,"INT",1,&Nclumps);
    //err += MPI_Bcast(&Nclumps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else {
    err += COMM->broadcast_data(0,"INT",1,&Nclumps);
    //err += MPI_Bcast(&Nclumps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // set up clumps
    if (cl) rep.error("I'm not root, but cl already initialised!",cl);
    cl = new struct clump [Nclumps];
    if (!cl) rep.error("Initialising clumps",cl);
  }
  if (err) {cerr<<"didn't set up clumps (pllel)\n"; return err;}

  //
  // Now broadcast the clump list from root to all other processors.
  //int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
  //
  double d[Nclumps];
  // overdensity
  if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].overdensity;
  err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
  //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].overdensity = d[i];

  // mass
  if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].mass;
  err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
  //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].mass = d[i];

  // centre
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].centre[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].centre[nd] = d[i];
  }

  // size
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].size[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].size[nd] = d[i];
  }

  // angle
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].ang[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].ang[nd] = d[i];
  }

  // rotation matrix
  for (int x=0; x<MAX_DIM; x++) {
    for (int y=0; y<MAX_DIM; y++) {
      if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].rm[x][y];
      err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
      //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].rm[x][y] = d[i];
    }
  } 

  return err;
}


// ##################################################################
// ##################################################################



int IC_photevap_random_clumps::clumps_random_setup_pllel()
{
  int err=0;
  if (MCMD->get_myrank()==0) {
    err += IC_photevap_random_clumps::clumps_random_setup();
  }
  else {
    // set up clumps
    if (cl) rep.error("I'm not root, but cl already initialised!",cl);
    cl = new struct clump [Nclumps];
    if (!cl) rep.error("Initialising clumps",cl);
  }
  if (err) {cerr<<"didn't set up clumps (pllel)\n"; return err;}

  // Now broadcast the clump list from root to all other processors.
  //int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
  double d[Nclumps];
  // overdensity
  if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].overdensity;
  err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
  //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].overdensity = d[i];

  // centre
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].centre[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].centre[nd] = d[i];
  }

  // size
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].size[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].size[nd] = d[i];
  }

  // angle
  for (int nd=0; nd<MAX_DIM; nd++) {
    if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].ang[nd];
    err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
    //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].ang[nd] = d[i];
  }

  // rotation matrix
  for (int x=0; x<MAX_DIM; x++) {
    for (int y=0; y<MAX_DIM; y++) {
      if (MCMD->get_myrank()==0) for (int i=0;i<Nclumps;i++) d[i] = cl[i].rm[x][y];
      err += COMM->broadcast_data(0,"DOUBLE",Nclumps,d);
      //err += MPI_Bcast(d, Nclumps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (MCMD->get_myrank()!=0) for (int i=0;i<Nclumps;i++) cl[i].rm[x][y] = d[i];
    }
  } 
  return err;
}
#endif //PARALLEL


// ##################################################################
// ##################################################################



void IC_photevap_random_clumps::print_clump(struct clump *rc)
{
  cout <<"--clump overdensity:"<<rc->overdensity<<"  mass:"<<rc->mass<<endl;
  rep.printVec("Centre",rc->centre,SimPM->ndim);
  rep.printVec("Radius",rc->size,SimPM->ndim);
  rep.printVec("Angles",rc->ang,SimPM->ndim);
  cout <<"-------------------\n";
  return;
}


// ##################################################################
// ##################################################################



int IC_photevap_random_clumps::clumps_random_set_dens(class cell *c
						      )
{
  int err=0;
  double x0[ndim],x1[ndim], dpos[ndim];
  CI.get_dpos(c,dpos);
  for (int j=0;j<Nclumps;j++) {
    //
    // get distance from centre to clump.
    // If we have periodic BCs, need to do a bit more work...
    //
    //if (SimPM->typeofbc.find("YNper") != std::string::npos) {
    if (SimPM->BC_YN == "periodic") {
      x0[XX] = dpos[XX]-cl[j].centre[XX];
      double temp;
      for (int i=1;i<ndim;i++) {
	temp = dpos[i]-cl[j].centre[i];
	if (temp>0.0) {
	  if (temp >  SimPM->Range[i]/2.) x0[i] = temp -SimPM->Range[i];
	  else                           x0[i] = temp;
	}
	else {
	  if (temp < -SimPM->Range[i]/2.) x0[i] = temp +SimPM->Range[i];
	  else                           x0[i] = temp;
	}
      }
    }
    else {
      // not periodic bcs.
      for (int i=0;i<ndim;i++) x0[i] = dpos[i]-cl[j].centre[i];
    }

    // rotate clump.
    for (int u=0;u<ndim;u++) {
      x1[u] = 0.0;
      for (int v=0;v<ndim;v++)
	x1[u] += x0[v] *cl[j].rm[u][v];
    }
    
    // add density from clump.
    if (profile==0) {
      // top-hat profile
      /** \todo Need to make this trace ellipses, not postage stamps!!! */
      bool inside=true;
      for (int k=0;k<ndim;k++) if (fabs(x1[k])> cl[j].size[k]) inside=false;
      if (inside) {
	c->P[RO] += ambient[RO]*cl[j].overdensity;
      }
    }
    else if (profile==1) {
      //
      // Gaussian Profile: Set exponential factor based on how far from centre we are.
      // Then set density based on that.
      // Set tracers according to fraction of mass that is ambient and clump:
      //
      double ef=0.0;
      for (int v=0;v<ndim;v++) ef += x1[v]*x1[v]/cl[j].size[v]/cl[j].size[v];
      c->P[RO] += ambient[RO]*cl[j].overdensity*exp(-ef);
    }
    else rep.error("Bad profile id in parameter-file",profile);
  }
  return err;
}


// ##################################################################
// ##################################################################




