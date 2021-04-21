/// \file basic_tests.cc
/// \author Jonathan Mackey
///
/// File for setting up some basic test problems which won't be ever used for
/// actual physics sims, just for testing the code.
///
/// 2009-12-07 JM: added switch in calling FieldLoop test to allow it to
/// have a non-zero Vz, which provides a more stringent test.
///
///
/// 2010.10.05 JM: Added ambient medium parameters to be read from paramter file
///   for the "uniform" initial conditions.
/// - 2013.10.17 JM: Added optional core/radial slope isothermal
///    sphere for the UNIFORM medium parameters (with optional
///    expansion velocity, for modelling winds).
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.05 JM: Added pion_flt datatype.
/// - 2016.02.22 JM: Updated Field Loop test to use new grid.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "coord_sys/VectorOps.h"
#include "ics/icgen.h"
#include "ics/icgen_base.h"
using namespace std;
#include <sstream>

// ##################################################################
// ##################################################################

IC_basic_tests::IC_basic_tests() {}

// ##################################################################
// ##################################################################

IC_basic_tests::~IC_basic_tests() {}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_data(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  int err = 0;

  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!", ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams", rp);

  IC_basic_tests::eqns = SimPM->eqntype;
  if (eqns == EQEUL)
    eqns = 1;
  else if (eqns == EQMHD || eqns == EQGLM || eqns == EQFCD)
    eqns = 2;
  else
    rep.error("Bad equations", eqns);

  int ndim   = SimPM->ndim;
  string ics = rp->find_parameter("ics");

  if (ics == "")
    rep.error("didn't get any ics to set up.", ics);
  else if (ics == "Uniform" || ics == "uniform") {
    // cout <<"\t\tSetting up Uniform ICs.\n";
    err += setup_uniformgrid(rrp, ggg);
  }
  else if (ics == "AdvectSineWave") {
    // cout <<"\t\tSetting up HD/MHD advection of clump with SINE WAVE
    // variation in VY.\n";
    err += setup_sinewave_velocity();
  }
  else if (ics == "Advection") {
    // cout <<"\t\tSetting up HD/MHD advection of overdense clump.\n";
    if (ndim != 2 && ndim != 3)
      rep.error("only know 2d/3d advected grids", ndim);
    err += setup_advection();
  }
  else if (ics == "DivBPeak") {
    // cout <<"\t\tSetting up divB peak problem.\n";
    err += setup_divBpeak();
  }
  else if (
      ics == "FieldLoop" || ics == "FieldLoopVz" || ics == "FieldLoopStatic") {
    // cout <<"\t\tSetting up Advection of Magnetic Field Loop problem.\n";
    //
    // We pass in the z-velocity to the fieldloop problem:
    //
    if (ics == "FieldLoop")
      err += setup_FieldLoop(0.0);
    else if (ics == "FieldLoopVz")
      err += setup_FieldLoop(1.0);
    else
      err += setup_FieldLoop(-1.0);
  }
  else if (ics == "OrszagTang") {
    // cout <<"\t\tSetting up Orszag-Tang test problem.\n";
    err += setup_OrszagTang();
  }
  else if (ics == "DoubleMachRef") {
    // cout <<"\t\tSetting up Double Mach Reflection test problem.\n";
    err += setup_DoubleMachRef();
  }
  else if (ics == "KelvinHelmholz") {
    // cout <<"\t\tSetting up Kelvin-Helmholz Instability problems.\n";
    err += setup_KelvinHelmholz();
  }
  else if (ics == "KelvinHelmholzStone") {
    // cout <<"\t\tSetting up Stone's Kelvin-Helmholz Instability
    // problem.\n";
    err += setup_KelvinHelmholz_Stone();
  }
  else if (ics == "LiskaWendroffImplosion") {
    cout << "Setting up Liska & Wendroff (2003) implosion test\n";
    err += setup_LWImplosion();
  }
  else
    rep.error("Don't know what Initial Condition is!", ics);

  // Add noise to data?  Smooth data?
  double noise = 0;
  int smooth   = 0;
  ics          = rp->find_parameter("noise");
  if (ics != "")
    noise = atof(ics.c_str());
  else
    noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number", noise);
  if (noise > 0) err += AddNoise2Data(gg, *SimPM, 2, noise);

  ics = rp->find_parameter("smooth");
  if (ics != "")
    smooth = atoi(ics.c_str());
  else
    smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number", smooth);
  if (smooth > 0) err += SmoothData(smooth);

  return err;
}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_uniformgrid(
    class ReadParams *rrp,    ///< pointer to parameter list.
    class GridBaseClass *ggg  ///< pointer to grid
)
{
  // SimPM->typeofop=1; // text output
  //
  // Get ambient medium properties:
  //
  double Amb[SimPM->nvar];

  string seek, str;
  seek = "UNIFORM_ambRO";
  str  = rrp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  Amb[RO] = atof(str.c_str());

  seek = "UNIFORM_ambPG";
  str  = rrp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  Amb[PG] = atof(str.c_str());

  seek = "UNIFORM_ambVX";
  str  = rrp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  Amb[VX] = atof(str.c_str());

  seek = "UNIFORM_ambVY";
  str  = rrp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  Amb[VY] = atof(str.c_str());

  seek = "UNIFORM_ambVZ";
  str  = rrp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  Amb[VZ] = atof(str.c_str());

  if (SimPM->eqntype == EQMHD || SimPM->eqntype == EQGLM
      || SimPM->eqntype == EQFCD) {
    seek = "UNIFORM_ambBX";
    str  = rrp->find_parameter(seek);
    if (str == "") rep.error("didn't find parameter", seek);
    Amb[BX] = atof(str.c_str());

    seek = "UNIFORM_ambBY";
    str  = rrp->find_parameter(seek);
    if (str == "") rep.error("didn't find parameter", seek);
    Amb[BY] = atof(str.c_str());

    seek = "UNIFORM_ambBZ";
    str  = rrp->find_parameter(seek);
    if (str == "") rep.error("didn't find parameter", seek);
    Amb[BZ] = atof(str.c_str());

    if (SimPM->eqntype == EQGLM) Amb[SI] = 0.0;

#ifdef NEW_B_NORM
    // convert from CGS to internal units (no factors of 4pi)
    Amb[BX] /= sqrt(4.0 * M_PI);
    Amb[BY] /= sqrt(4.0 * M_PI);
    Amb[BZ] /= sqrt(4.0 * M_PI);
#endif
  }

  // tracer variables
  ostringstream temp;
  for (int t = 0; t < SimPM->ntracer; t++) {
    temp.str("");
    temp << "UNIFORM_ambTR" << t;
    seek = temp.str();
    str  = rrp->find_parameter(seek);
    if (str != "")
      Amb[t + SimPM->ftr] = atof(str.c_str());
    else
      Amb[t + SimPM->ftr] = -1.0e99;
  }

  //
  // optional radial slope in density from a central point.
  //
  double centre[MAX_DIM];
  double core_radius = 0.0, radial_slope = 0.0, radial_vel = 0.0;
  bool use_core = true;

  seek = "UNIFORM_radial_slope";
  str  = rrp->find_parameter(seek);
  if (str == "") {
    radial_slope = 0.0;
    use_core     = false;
  }
  else {
    radial_slope = atof(str.c_str());
    // cout <<"\t***\t*** slope = "<<radial_slope<<"\n";
  }

  seek = "UNIFORM_radial_velocity";
  str  = rrp->find_parameter(seek);
  if (str == "") {
    radial_vel = 0.0;
  }
  else {
    radial_vel = atof(str.c_str());
  }

  seek = "UNIFORM_core_radius";
  str  = rrp->find_parameter(seek);
  if (str == "") {
    core_radius = 0.0;
    use_core    = false;
  }
  else {
    core_radius = atof(str.c_str());
    // cout <<"\t***\t*** core = "<< core_radius <<"\n";
  }

  seek = "UNIFORM_core_centre_XX";
  str  = rrp->find_parameter(seek);
  if (str == "")
    centre[XX] = 0.0;
  else
    centre[XX] = atof(str.c_str());

  seek = "UNIFORM_core_centre_YY";
  str  = rrp->find_parameter(seek);
  if (str == "")
    centre[YY] = 0.0;
  else
    centre[YY] = atof(str.c_str());

  seek = "UNIFORM_core_centre_ZZ";
  str  = rrp->find_parameter(seek);
  if (str == "")
    centre[ZZ] = 0.0;
  else
    centre[ZZ] = atof(str.c_str());

  //
  // Check that the radial slope is non-zero.  If it is zero, then
  // switch off the core.
  //
  if (use_core && pconst.equalD(radial_slope, 0.0)) use_core = false;

  // cout <<"\t\tAssigning values to data.\n";
  class cell *cpt = ggg->FirstPt();
  double distance, dpos[MAX_DIM];
  do {
    for (int v = 0; v < SimPM->nvar; v++)
      cpt->P[v] = cpt->Ph[v] = Amb[v];
    //
    // If we are using a core/slope for the density, then set up
    // an isothermal sphere with a constant density core.
    // rho = rho_0 /(1+ (r_core/r)^slope)
    //
    CI.get_dpos(cpt, dpos);
    distance = ggg->distance(centre, dpos);
    if (use_core) {
      cpt->P[RO] /= (1.0 + exp(radial_slope * log(core_radius / distance)));
      cpt->P[PG] /= (1.0 + exp(radial_slope * log(core_radius / distance)));
      cpt->P[VX] = radial_vel * (dpos[XX] - centre[XX]) / distance;
      cpt->P[VY] = radial_vel * (dpos[YY] - centre[YY]) / distance;
      cpt->P[VZ] = radial_vel * (dpos[ZZ] - centre[ZZ]) / distance;
      // cout <<"core density = "<<cpt->P[RO]<<", d="<<distance<<"\n";
    }
  } while ((cpt = ggg->NextPt(cpt)) != NULL);
  //  cpt = ggg->FirstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;}
  // while  ( (cpt=ggg->NextPt(cpt))!=NULL);
  // cout <<"\t\tGot through data successfully.\n";
  return (0);
}

// ##################################################################
// ##################################################################

// ************************************************
// ******** SINE WAVE IN VELOCITY, ADVECTION ******
// ************************************************

int IC_basic_tests::setup_sinewave_velocity()
{
  //   string seek, str;
  //   seek = "NDadv_thetaXY";
  //   str = rp->find_parameter(seek);
  //   if (str=="") rep.error("didn't find parameter",seek);
  //   double thetaXY = M_PI/180.* atof(str.c_str());
  //   seek = "NDadv_thetaXZ";
  //   str = rp->find_parameter(seek);
  //   if (str=="") rep.error("didn't find parameter",seek);
  //   double thetaXZ = M_PI/180.* atof(str.c_str());

  //   cout <<"Using angle to X-Y grid of "<<thetaXY<<" radians, X-Z of
  //   "<<thetaXZ<<" radians.\n"; cout <<"The X-Z angle is from the z-axis to
  //   the vector. (polar angle)\n"; cout <<"The X-Y angle is from the x-axis
  //   to the projection of vector onto xy plane (azimuthal angle)\n";

  int ndim = gg->Ndim();  // int nvar=gg->Nvar();
  if (ndim != 2 && ndim != 3)
    rep.error("Bad ndim in setup_sinewave_velocity()", ndim);

  //   if (ndim==2 && !pconst.equalD(thetaXZ,M_PI/2.)) {
  //     rep.warning("Given 3d angle, but 2d sim.  setting 3d angle to
  //     zero.",0,thetaXZ); thetaXZ=M_PI/2.;
  //   }

  //
  // initial conditions:  Total velocity =1.
  // vx = sin(xz)cos(xy), vy=sin(xz)sin(xy), vz=cos(xz)
  // xz = angle from z-axis to vec(v), xy = angle from x-axis to v projected
  // onto xy plane
  //
  double vx = 1.0;  // sin(thetaXZ)*cos(thetaXY);
  double vy = 1.0;  // sin(thetaXZ)*sin(thetaXY);
  double vz = 1.0;  // cos(thetaXZ);
  // cout <<"\t\tvx,vy,vz = "<<vx<<", "<<vy<<", "<<vz<<endl;
  double pout;
  pout = 1.0;
  double rhoin, rhoout;
  rhoin  = 10.0;
  rhoout = 1.0;

  if (vx > vy)
    SimPM->finishtime = 3. * SimPM->Range[0] / vx;
  else
    SimPM->finishtime = 5. * SimPM->Range[0] / vy;

  // Circle setup
  double centre[ndim];
  for (int i = 0; i < ndim; i++) {
    centre[i] = (SimPM->Xmax[i] - SimPM->Xmin[i]) / 2.;
  }
  double radius = (SimPM->Xmax[0] - SimPM->Xmin[0])
                  / 10.;  // radius is 1/5 of box diameter in x-dir.
  // Set up the inside_sphere class, with 100 subpoints per cell.
  int nsub;
  if (ndim == 2)
    nsub = 100;
  else
    nsub = 32;
  class inside_sphere stest(centre, radius, SimPM->dx, nsub, ndim);
  double vfrac;

  // data
  // cout <<"\t\tAssigning primitive vectors.\n";
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    cpt->P[RO] = rhoout;
    cpt->P[PG] = pout;
    cpt->P[VX] = vx;
    cpt->P[VY] = vy * sin(2.0 * M_PI * CI.get_dpos(cpt, YY) / SimPM->Range[YY]);
    cpt->P[VZ] = vz;
    if (SimPM->eqntype == EQMHD || SimPM->eqntype == EQGLM
        || SimPM->eqntype == EQFCD) {
      cpt->P[BX] = vx;
      cpt->P[BY] = vy;
      cpt->P[BZ] = vz;  // field aligned with flow
    }
    for (int i = 0; i < SimPM->ntracer; i++)
      cpt->P[SimPM->ftr + i] = 1.;

    // This is where I set the state inside the overdense clump.
    if ((vfrac = stest.volumeFraction(cpt)) > 0) {
      cpt->P[RO] = vfrac * (rhoin) + (1. - vfrac) * rhoout;
      for (int i = 0; i < SimPM->ntracer; i++)
        cpt->P[SimPM->ftr + i] = vfrac * (-1.) + (1. - vfrac) * 1.;
      //       cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
    }
  } while ((cpt = gg->NextPt(cpt)) != NULL);
  // cout <<"\t\tGot through data successfully.\n";
  // Data done.

  return 0;
}  // setup_sinewave()

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_advection()
{
  string seek, str;
  seek = "NDadv_thetaXY";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double thetaXY = M_PI / 180. * atof(str.c_str());
  seek           = "NDadv_thetaXZ";
  str            = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double thetaXZ = M_PI / 180. * atof(str.c_str());

  // cout <<"Using angle to X-Y grid of "<<thetaXY<<" radians, X-Z of
  // "<<thetaXZ<<" radians.\n"; cout <<"The X-Z angle is from the z-axis to
  // the vector. (polar angle)\n"; cout <<"The X-Y angle is from the x-axis to
  // the projection of vector onto xy plane (azimuthal angle)\n";

  int ndim = gg->Ndim();  // int nvar=gg->Nvar();
  if (ndim != 2 && ndim != 3) rep.error("Bad ndim in setupNDadvectionHD", ndim);

  if (ndim == 2 && !pconst.equalD(thetaXZ, M_PI / 2.)) {
    rep.warning(
        "Given 3d angle, but 2d sim.  setting 3d angle to zero.", 0, thetaXZ);
    thetaXZ = M_PI / 2.;
  }

  //
  // initial conditions:  Total velocity =1.
  // vx = sin(xz)cos(xy), vy=sin(xz)sin(xy), vz=cos(xz)
  // xz = angle from z-axis to vec(v), xy = angle from x-axis to v projected
  // onto xy plane
  //
  double vx = sin(thetaXZ) * cos(thetaXY);
  double vy = sin(thetaXZ) * sin(thetaXY);
  double vz = cos(thetaXZ);
  cout << "\t\tvx,vy,vz = " << vx << ", " << vy << ", " << vz << endl;
  double pout;
  pout = 1.0;
  double rhoin, rhoout;
  rhoin  = 10.0;
  rhoout = 1.0;

  if (vx > vy)
    SimPM->finishtime = 3. * SimPM->Range[0] / vx;
  else
    SimPM->finishtime = 5. * SimPM->Range[0] / vy;

  // Circle setup
  double centre[ndim];
  for (int i = 0; i < ndim; i++)
    centre[i] = (SimPM->Xmax[i] - SimPM->Xmin[i]) / 2.;
  double radius = (SimPM->Xmax[0] - SimPM->Xmin[0])
                  / 10.;  // radius is 1/5 of box diameter in x-dir.
  // Set up the inside_sphere class, with 100 subpoints per cell.
  int nsub;
  if (ndim == 2)
    nsub = 100;
  else
    nsub = 32;
  class inside_sphere stest(centre, radius, SimPM->dx, nsub, ndim);
  double vfrac;

  // data
  // cout <<"\t\tAssigning primitive vectors.\n";
  class cell *cpt = gg->FirstPt();
  do {
    // Set values of primitive variables.
    cpt->P[RO] = rhoout;
    cpt->P[PG] = pout;
    cpt->P[VX] = vx;
    cpt->P[VY] = vy;
    cpt->P[VZ] = vz;
    if (SimPM->eqntype == EQMHD || SimPM->eqntype == EQGLM
        || SimPM->eqntype == EQFCD) {
      cpt->P[BX] = vx;
      cpt->P[BY] = vy;
      cpt->P[BZ] = vz;  // field aligned with flow
    }
    for (int i = 0; i < SimPM->ntracer; i++)
      cpt->P[SimPM->ftr + i] = 1.;

    // This is where I set the state inside the blast radius.
    if ((vfrac = stest.volumeFraction(cpt)) > 0) {
      //       cpt->P[PG] = vfrac*(pin) + (1.-vfrac)*cpt->P[PG];
      cpt->P[RO] = vfrac * (rhoin) + (1. - vfrac) * rhoout;
      for (int i = 0; i < SimPM->ntracer; i++)
        cpt->P[SimPM->ftr + i] = vfrac * (-1.) + (1. - vfrac) * 1.;
      //       cout <<"Setting cell "<<cpt->id<<" to internal value.\n";
    }
  } while ((cpt = gg->NextPt(cpt)) != NULL);
  //  cpt = firstPt();
  //  do {cout <<"cpt.rho = "<<cpt->P[RO]<<endl;} while  (
  //  (cpt=nextPt(cpt))!=NULL);
  // cout <<"\t\tGot through data successfully.\n";
  // Data done.

  return (0);
}  // setup_advection

// ##################################################################
// ##################################################################

/**************************************************************************/
// Div B peak.
/**************************************************************************/
int IC_basic_tests::setup_divBpeak()
{
  // See Dedner et al. (2002) JCP, 175, 645 for details of the problem and
  // reference results.
  int ndim = gg->Ndim();
  if (ndim != 2) rep.error("divBpeak only works in 2D", ndim);
  int nvar  = gg->Nvar();
  double *s = new double[nvar];
  s[RO]     = 1.0;
  s[VX] = s[VY] = 1.;
  s[VZ]         = 0.0;
  s[PG]         = 6.0;
  s[BY]         = 0.0;
  s[BZ]         = 1.;
  s[BX]         = 0.;

  // SimPM->BC_XN = "periodic";
  // SimPM->BC_XP = "periodic";
  // SimPM->BC_YN = "periodic";
  // SimPM->BC_YP = "periodic";
  // SimPM->BC_Nint = 0;
  SimPM->gamma = 5. / 3.;
  if (fabs(SimPM->Xmin[XX] + 0.5) > 1.e-6
      || fabs(SimPM->Xmax[XX] - 1.5) > 1.e-6) {
    cout << fabs(SimPM->Xmin[XX] + 0.5) << "\t" << fabs(SimPM->Xmax[XX] - 1.5)
         << endl;
    cout << SimPM->Xmin[XX] << "\t" << SimPM->Xmax[XX] << endl;
    rep.error(
        "Set bounds properly for divBpeak!!! x=[-.5,1.5] y=[-.5,1.5]",
        SimPM->Xmin[XX]);
  }
  class cell *c = gg->FirstPt();
  double r2     = 0;
  double dpos[MAX_DIM];
  CI.get_dpos(c, dpos);
  // rep.printVec("position of first cell",dpos,ndim);
  do {
    CI.get_dpos(c, dpos);
    for (int v = 0; v < nvar; v++)
      c->P[v] = s[v];
    r2 = dpos[XX] * dpos[XX] + dpos[YY] * dpos[YY];
    //     c->P[BX] = 1.3*exp(-r2/.01);
    if (r2 < 1. / 8.)
      r2 = 4096. * r2 * r2 * r2 * r2 - 128. * r2 * r2 + 1.;
    else
      r2 = 0.;
    if (r2 < 0.) r2 = 0.;
    c->P[BX] = r2;
  } while ((c = gg->NextPt(c)) != 0);
  delete[] s;
  return (0);
}
/**************************************************************************/
// Div B peak.
/**************************************************************************/

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_FieldLoop(double vz  ///< Z-velocity of fluid
)
{
  cout << "\tSetting up Field Loop test, from Stone's code test page:\n";
  cout << "\t\thttp://www.astro.princeton.edu/~jstone/tests/field-loop/"
          "Field-loop.html\n";
  cout << "\tAlso see http://www.dias.ie/~fdc/MHDCodeComp/adv.html\n";
  int ndim = SimPM->ndim;
  if (ndim != 2) {
    rep.error("Bad ndim in setup_FieldLoop", ndim);
  }
  if (SimPM->eqntype != EQMHD && SimPM->eqntype != EQGLM
      && SimPM->eqntype != EQFCD) {
    rep.error(
        "Advection of Field Loop must be mhd! bad eqntype", SimPM->eqntype);
  }
  SimPM->gamma = 5. / 3.;  // just to make sure.
  // SimPM->BC_XN = "periodic";
  // SimPM->BC_XP = "periodic";
  // SimPM->BC_YN = "periodic";
  // SimPM->BC_YP = "periodic";
  // cout <<"\tMake sure x=[-1,1] and y=[-0.5,0.5]\n";

  double A_max = 0.001;  // Peak in vector potential.
  double vel = 2.0, rho = 1.0, p_g = 1.0;
  // double flow_angle=60.0*M_PI/180.0;
  double radius = 0.3;
  double dist   = 0.0;
  double centre[ndim];
  for (int v = 0; v < ndim; v++)
    centre[v] = 0.0;
  double dpos[ndim];

  if (vz < 0) {
    // cout <<"\tWARNING - negative vz received, so setting up static
    // problem!\n";
    vel = 0.0;
    vz  = 0.0;
  }

  // cout <<"Assigning primitive vectors.\n";
  class cell *c = gg->FirstPt();
  do {
    c->P[RO] = rho;
    c->P[PG] = p_g;
    c->P[VX] = vel;        // vel*sin(flow_angle);
    c->P[VY] = vel / 2.0;  // vel*cos(flow_angle);
    c->P[VZ] = vz;         // If vz!=0, this tests if B_z gets contaminated.
    CI.get_dpos(c, dpos);
    dist = gg->distance(centre, dpos);

    // cout <<"cell id="<<c->id<<",  dist="<<dist;
    // cout <<",  radius="<< radius<<",  ";
    // rep.printVec("dpos",dpos,ndim);

    //
    // poor man's b-field (has divB errors)
    //
    if (dist < radius) {
      c->P[BX] = A_max * dpos[YY] / dist;
      c->P[BY] = -A_max * dpos[XX] / dist;
    }
    else {
      c->P[BX] = 0.0;
      c->P[BY] = 0.0;
    }
    c->P[BZ] = 0.0;
    if (SimPM->eqntype == EQGLM) c->P[SI] = 0.0;
    for (int i = 0; i < SimPM->ntracer; i++) {
      c->P[SimPM->ftr + i] = 1.;
    }
    //
    // vector potential...
    //
    if (dist < radius) {
      c->Ph[BZ] = A_max * (radius - dist);
      // cout <<"setting A(z)="<<c->Ph[BZ]<<"\n";
    }
    else {
      c->Ph[BZ] = 0.0;
    }
    c->Ph[BX] = c->Ph[BY] = 0.0;
  } while ((c = gg->NextPt(c)) != 0);
  // cout <<"Got through data successfully.\n";
  // Data done.

  //
  // Take curl of vector...
  //
  int els[3] = {BX, BY, BZ};
  pion_flt ans[3];
  class VectorOps_Cart *vec = new VectorOps_Cart(ndim);
  c                         = gg->FirstPt();
  do {

    // if (!pconst.equalD(c->Ph[BZ],0.0)) {
    //  cout <<"A(z) = "<<c->Ph[BZ]<<"  ";
    //  rep.printVec("ans",ans,3);
    //  CI.print_cell(c);
    //}

    if (!c->isedge) {
      vec->Curl(c, 1, els, gg, ans);
      c->P[BX] = ans[0];
      c->P[BY] = ans[1];
      c->P[BZ] = ans[2];
    }
    else {
      c->P[BX] = 0;
      c->P[BY] = 0;
      c->P[BZ] = 0;
    }
  } while ((c = gg->NextPt(c)) != 0);
  delete vec;
  vec = 0;

  return (0);
}

// ##################################################################
// ##################################################################

///
/// Set up Orszag-Tang Vortex problem (from Dai & Woodward 1998,APJ,494,317)
/// This assumes the grid is unit size in both directions, so it automatically
/// works in serial and in parallel.
///
int IC_basic_tests::setup_OrszagTang()
{
  // set plasma beta parameter (ratio of gas to magnetic pressure)
  // set mach number of gas (hydro mach no.).
  string seek, str;
  seek = "OTVbeta";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double otvbeta = atof(str.c_str());

  seek = "OTVmach";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double otvmach = atof(str.c_str());

  int ndim = gg->Ndim();
  cout << "\tSetting up Orszag-Tang vortex problem with plasma beta = "
       << otvbeta;
  cout << " and flow mach no. = " << otvmach << endl;
  if (ndim != 2) rep.error("Bad ndim in setup_OrszagTang", ndim);
  if (SimPM->eqntype != EQMHD && SimPM->eqntype != EQGLM
      && SimPM->eqntype != EQFCD)
    rep.error("O-T vortex must be mhd! bad eqntype", SimPM->eqntype);

  SimPM->gamma = 5. / 3.;                                // just to make sure.
  double p0    = otvbeta / 2.;                           // constant pressure.
  double d0    = SimPM->gamma * otvmach * otvmach * p0;  // constant density.
  // cout <<"Assigning primitive vectors.\n";
  double dpos[ndim];
  class cell *c = gg->FirstPt();
  do {
    CI.get_dpos(c, dpos);
    // Set values of primitive variables.
    c->P[RO] = d0;
    c->P[PG] = p0;
    c->P[VX] = -sin(2. * M_PI * dpos[YY]);
    c->P[VY] = sin(2. * M_PI * dpos[XX]);
    c->P[VZ] = 0.;
    c->P[BX] = -sin(2. * M_PI * dpos[YY]);
    c->P[BY] = sin(4. * M_PI * dpos[XX]);
    c->P[BZ] = 0.0;
    if (SimPM->eqntype == EQGLM) c->P[SI] = 0.0;
    for (int i = 0; i < SimPM->ntracer; i++) {
      c->P[SimPM->ftr + i] = 1.;
    }
  } while ((c = gg->NextPt(c)) != 0);
  // cout <<"Got through data successfully.\n";
  // Data done.
  return (0);
}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_DoubleMachRef()
{
  // Set up Woodward & Colella (1984) Double Mach Reflection.
  // dmrmach is the mach number of the incident shock.
  string seek, str;
  seek = "DMRmach";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double dmrmach = atof(str.c_str());

  seek = "DMRtheta";
  str  = rp->find_parameter(seek);
  if (str == "") rep.error("didn't find parameter", seek);
  double dmrtheta = atof(str.c_str());

  int ndim = gg->Ndim();
  cout << "Setting up Double Mach Reflection problem...  ";
  if (dmrmach <= 1) rep.error("Mach number must be >1", dmrmach);
  if (dmrtheta < 0.1 || dmrtheta > 89.99)
    rep.error("Angle must be between 0 and 90 degrees", dmrtheta);
  cout << "with mach no. = " << dmrmach << " and angle " << dmrtheta
       << " degrees to x-axis.\n";
  if (ndim != 2) rep.error("Bad ndim in setup_DoubleMachRef", ndim);
  if (SimPM->eqntype != EQEUL)
    rep.error("DMR must be euler equations!", SimPM->eqntype);
  // SimPM->BC_XN = "inflow";
  // SimPM->BC_XP = "outflow";
  // SimPM->BC_YN = "reflecting";
  // SimPM->BC_YP = "DMR";
  // SimPM->BC_Nint = 1;
  SimPM->gamma = 1.4;
  // cout <<"*NB*: Assuming grid dimensions are {[0,4],[0,1]}; if not things
  // may/will go wrong!\n";

  dmrtheta *= M_PI / 180.0;

  // override mach no and theta:
  dmrmach  = 10.0;
  dmrtheta = M_PI / 3.0;
  // cout <<"Override: hardwired to mach no. = "<<dmrmach;
  // cout <<" and angle "<<dmrtheta*180./M_PI<<" degrees to x-axis.\n";
  // cout <<"If this is a problem, fix the boundary conditions in the
  // code.!\n";

  double x0  = 1. / 6.;  // initial location of shock.
  double ro0 = 1.4;
  double pg0 = 1.0;
  double vx0 = 0.0, vy0 = 0.0, vz0 = 0.0;
  double pg1 = pg0 * (2. * SimPM->gamma * dmrmach * dmrmach - SimPM->gamma + 1.)
               / (SimPM->gamma + 1.);
  double alpha = (SimPM->gamma + 1.) / (SimPM->gamma - 1.);
  double ro1   = ro0 * (1. + alpha * pg1 / pg0) / (alpha + pg1 / pg0);
  double vx1   = vx0
               + sin(dmrtheta) * (pg1 / pg0 - 1.)
                     * sqrt(SimPM->gamma * pg0 / ro0) / SimPM->gamma / dmrmach;
  double vy1 = vy0
               - cos(dmrtheta) * (pg1 / pg0 - 1.)
                     * sqrt(SimPM->gamma * pg0 / ro0) / SimPM->gamma / dmrmach;
  double vz1 = vz0;
  // cout <<"postshock state: ro="<<ro1<<", pg="<<pg1<<", vx="<<vx1<<",
  // vy="<<vy1<<endl;

  double xs     = 0.0;
  class cell *c = gg->FirstPt();
  double dpos[ndim];
  do {
    CI.get_dpos(c, dpos);
    xs = x0 + dpos[YY] / tan(dmrtheta);
    if (dpos[XX] <= xs) {
      c->P[RO] = ro1;
      c->P[PG] = pg1;
      c->P[VX] = vx1;
      c->P[VY] = vy1;
      c->P[VZ] = vz1;
    }
    else {
      c->P[RO] = ro0;
      c->P[PG] = pg0;
      c->P[VX] = vx0;
      c->P[VY] = vy0;
      c->P[VZ] = vz0;
    }
  } while ((c = gg->NextPt(c)) != 0);
  // cout <<"Got through data successfully.\n";
  // Data done.

  return (0);
}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_KelvinHelmholz_Stone()
{
  cout << "KH Instability: assuming x=[-0.5,0.5], y=[-0.5,0.5]\n";
  int err  = 0;
  int ndim = gg->Ndim();
  if (ndim != 2) rep.error("KH needs 2D problem domain", ndim);

  // The following is for Jim Stone's test at:
  // http://www.astro.princeton.edu/~jstone/tests/kh/kh.html
  // SimPM->typeofbc = "XNper_XPper_YNper_YPper_";
  cout << "KH Instability: using periodic BCs everywhere.\n";
  SimPM->gamma    = 1.4;
  double pressure = 2.5;
  double vx_high  = -0.5;
  double vx_low   = 0.5;
  double rho_high = 1.0;
  double rho_low  = 2.0;
  double Bx =
      0.5 / sqrt(4. * M_PI);  // think this is right, but not sure about 4Pi
  int seed = 975;
#ifdef PARALLEL
  seed += MCMD->get_myrank();
#endif
  srand(seed);
  double noise_amp = 0.01;  // absolute amplitude of noise.

  class cell *c = gg->FirstPt();
  double dpos[ndim];
  do {
    CI.get_dpos(c, dpos);
    c->P[VY] = 0.0;
    c->P[PG] = pressure;
    if (eqns == 2) {
      c->P[BX] = Bx;
      c->P[BY] = c->P[BZ] = 0.0;
    }
    if (fabs(dpos[YY]) > 0.25) {
      c->P[RO] = rho_high;
      c->P[VX] = vx_high;
      for (int i = 0; i < SimPM->ntracer; i++)
        c->P[SimPM->ftr + i] = -1.;
    }
    else {
      c->P[RO] = rho_low;
      c->P[VX] = vx_low;
      for (int i = 0; i < SimPM->ntracer; i++)
        c->P[SimPM->ftr + i] = 1.;
    }
    c->P[VX] += noise_amp * (static_cast<double>(rand()) / RAND_MAX - 0.5);
    c->P[VY] += noise_amp * (static_cast<double>(rand()) / RAND_MAX - 0.5);
  } while ((c = gg->NextPt(c)) != 0);
  // cout <<"Got through data successfully.\n";
  // Data done.

  return err;
}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_KelvinHelmholz()
{
  cout << "KH Instability: assuming x=[-0.5,0.5], y=[-0.5,0.5]\n";
  int err  = 0;
  int ndim = gg->Ndim();
  if (ndim != 2) rep.error("KH needs 2D problem domain", ndim);

  // The following is for Frank, Jones, Ryu, \& Gaalaas, 1996, ApJ, 460, 777.
  // SimPM->typeofbc = "XNper_XPper_YNref_YPref_";
  SimPM->gamma    = 5. / 3.;
  double pressure = 0.6;
  double rho      = 1.0;
  double Ux       = 1.0;
  double Bx       = 0.2;   // weak field is 0.2, strong is 0.4
  double Uy_amp   = 0.01;  // Amplitude of velocity perturbation
  double a        = SimPM->Range[YY] / 25.0;  // Thickness of shear layer.

  class cell *c = gg->FirstPt();
  double dpos[ndim];
  do {
    CI.get_dpos(c, dpos);
    c->P[RO] = rho;
    c->P[PG] = pressure;
    if (eqns == 2) {
      c->P[BX] = Bx;
      c->P[BY] = c->P[BZ] = 0.0;
    }
    c->P[VX] = -0.5 * Ux * tanh((dpos[YY]) / a);  //-0.5*SimPM->Range[YY])/a);
    // This makes the two densities different...
    // c->P[RO] += 0.1*rho*(c->P[VX]/Ux);
    // This adds two sin wave perturbations with 1,2 and 15 periods in the
    // box, exponentially damped in the y-direction.
    c->P[VY] = 0.3333 * Uy_amp
               * (sin(30.0 * M_PI * dpos[XX] / SimPM->Range[XX])
                  + sin(2.0 * M_PI * dpos[XX] / SimPM->Range[XX])
                  + sin(4.0 * M_PI * dpos[XX] / SimPM->Range[XX]))
               * exp(-(dpos[YY] * dpos[YY]) / 4.0 / a / a);
    c->P[VZ] = 0.0;
  } while ((c = gg->NextPt(c)) != 0);
  // cout <<"Got through data successfully.\n";
  // Data done.

  return err;
}

// ##################################################################
// ##################################################################

int IC_basic_tests::setup_LWImplosion()
{
  // See Liska & Wendroff (2003,SIAM Journal on Scientific Computing,
  // 2003, Vol. 25, No. 3 : pp. 995-1017
  // Preprint at:
  // http://kfe.fjfi.cvut.cz/~liska/ps.gz/compare_euler_SISC_03.pdf
  //
  int err  = 0;
  int ndim = gg->Ndim();
  if (ndim != 2) rep.error("LWI needs 2D problem domain", ndim);

  SimPM->gamma    = 1.4;
  double pressure = 1.0;
  double rho      = 1.0;

  class cell *c = gg->FirstPt();
  double dpos[ndim];
  do {
    CI.get_dpos(c, dpos);
    c->P[RO] = rho;
    c->P[PG] = pressure;
    c->P[VX] = 0.0;
    c->P[VY] = 0.0;
    c->P[VZ] = 0.0;
    if (dpos[XX] < 0.15 && dpos[YY] < (0.15 - dpos[XX])) {
      c->P[RO] = 0.125;
      c->P[PG] = 0.140;
    }
  } while ((c = gg->NextPt(c)) != 0);
  return err;
}

// ##################################################################
// ##################################################################
