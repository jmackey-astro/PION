/// \file StarBench_test.cpp
/// \author Jonathan Mackey
///
/// File for setting up initial conditions for test problems to be
/// run as part of the StarBench code-testing workshop.
///
/// - 2013.01.10 JM: Started on file.
///     Added Contact Discontinuity tests.
/// - 2013.03.23-24 JM: Added Irradiated cloud and I-front
///    instability test problems.
/// - 2013.06.13 JM: Added StarBench cooling/shadowing test from
///    Pascal Tremblin.
/// - 2013.06.17 JM: Changed cooling/shadowing boundary condition
///    so that it is off-grid and I don't need to worry about it
///    in this file.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.05 JM: Added pion_flt datatype.
/// - 2016.05.02 JM: Added Cone-IF test, and planar IF test.
/// - 2016.05.15 JM: Added perturbations to planar IF test.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "ics/icgen.h"
#include "microphysics/microphysics_base.h"
#include "coord_sys/VectorOps.h"
#include "dataIO/dataio.h"
#include <sstream>


// ##################################################################
// ##################################################################



IC_StarBench_Tests::IC_StarBench_Tests()
{
  return;
}


// ##################################################################
// ##################################################################



IC_StarBench_Tests::~IC_StarBench_Tests()
{}


// ##################################################################
// ##################################################################




int IC_StarBench_Tests::setup_data(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg ///< pointer to grid
      )
{
  //
  // Make sure we have valid pointers to a computational grid and
  // to the list of parameters.
  //
  ICsetup_base::gg = ggg;
  if (!gg) rep.error("null pointer to grid!",ggg);

  ICsetup_base::rp = rrp;
  if (!rp) rep.error("null pointer to ReadParams",rp);

  if (SimPM->eqntype != EQEUL)
    rep.error("Bad equations",SimPM->eqntype);

  int err=0;
  string ics = rp->find_parameter("ics");
  if      (ics=="") {
    rep.error("didn't get any ics to set up.",ics);
  }
  else if (ics=="StarBench_ContactDiscontinuity1" ||
           ics=="StarBench_ContactDiscontinuity2" ||
           ics=="StarBench_ContactDiscontinuity3" ||
           ics=="StarBench_ContactDiscontinuity4") {
    cout <<"\t\tSetting up "<<ics<<" test.\n";
    err += setup_ContactDiscontinuity(rrp,ggg,ics);
  }
  else if (ics=="StarBench_IFI_testA" ||
           ics=="StarBench_IFI_testB" ||
           ics=="StarBench_IFI_testC") {
    cout <<"\t\tSetting up StarBench Planar ionisation front A.\n";
    err += setup_StarBench_IFI(rrp,ggg,ics);
  }
  else if (ics=="StarBench_IrrCloud_Uniform" ||
           ics=="StarBench_IrrCloud_IsoSph") {
    cout <<"\t\tSetting up StarBench Irradiated Cloud test.\n";
    err += setup_StarBench_IrrCl(rrp,ggg,ics);
  }
  else if (ics=="StarBench_TremblinCooling") {
    cout <<"\t\tSetting up StarBench Shadowing/Mixing/Cooling test.\n";
    err += setup_StarBench_TremblinCooling(rrp,ggg,ics);
  }
  else if (ics=="StarBench_Cone") {
    cout <<"\t\tSetting up StarBench Cone test.\n";
    err += setup_StarBench_Cone(rrp,ggg,ics);
  }
  else if (ics=="StarBench_IFI_V2") {
    cout <<"\t\tSetting up planar ionization front.\n";
    err += setup_StarBench_planarIF(rrp,ggg,ics);
  }


  //else if (ics=="") {
  //  cout <<"\t\tSetting up .\n";
  //  err += setup_(rrp,ggg);
  //}
  else rep.error("Don't know what Initial Condition is!",ics);
  
  if (err) rep.error("Test setup returned error",err);

  //
  // Add noise to data?  Smooth data?
  //
  int smooth=0; double noise=0.0;
  ics = rp->find_parameter("noise");
  if (ics!="") noise = atof(ics.c_str());
  else noise = -1;
  if (isnan(noise)) rep.error("noise parameter is not a number",noise);
  if (noise>0) {
    cout <<"\t\tNOISE: Adding random adiabatic noise at fractional level = "<<noise<<"\n";
    err+= AddNoise2Data(gg, *SimPM, 2,noise);
  }
  ics = rp->find_parameter("smooth");
  if (ics!="") smooth = atoi(ics.c_str());
  else smooth = -1;
  if (isnan(smooth)) rep.error("Smooth parameter not a number",smooth);
  if (smooth>0)err+= SmoothData(smooth);



  return 0;
}



// ##################################################################
// ##################################################################



int IC_StarBench_Tests::setup_ContactDiscontinuity(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  //
  // First see what test we are doing.
  //
  int test_id=-1;
  if      (test=="StarBench_ContactDiscontinuity1") test_id=1;
  else if (test=="StarBench_ContactDiscontinuity2") test_id=2;
  else if (test=="StarBench_ContactDiscontinuity3") test_id=3;
  else if (test=="StarBench_ContactDiscontinuity4") test_id=4;

  //
  // Check that we got something sensibe, that ndim is correct for
  // the test, and that we have a passive tracer.
  //
  if (test_id <0) {
    cerr <<"bad test! "<<test<<"\n";
    return 1;
  }
  if (test_id <=2  && SimPM->ndim != 1) {
    rep.error("ContactDiscontinuity1/2 is 1D test",SimPM->ndim);
  }
  if (test_id >=3 && SimPM->ndim != 2) {
    rep.error("ContactDiscontinuity3/4 is 2D test",SimPM->ndim);
  }

  if (SimPM->ntracer != 1) {
    rep.error("Need exactly one colour tracer for CD tests",
              SimPM->ntracer);
  }
  
  //
  // Get advection velocity from the parameterfile.
  //
  double Vel[SimPM->ndim];
  cout << "** It is assumed all parameters are correctly set **\n";

  string seek = rrp->find_parameter("StarBench_ContDisc_VX");
  if (seek=="") rep.error("Need parameter StarBench_ContDisc_VX",1);
  else Vel[XX] = atof(seek.c_str());

  if (SimPM->ndim>1) {
    string seek = rrp->find_parameter("StarBench_ContDisc_VY");
    if (seek=="") rep.error("Need parameter StarBench_ContDisc_VY",2);
    else Vel[YY] = atof(seek.c_str());
  }

  cell *c=ggg->FirstPt();
  double pos[SimPM->ndim];

  if      (test_id==1) {
    do {
      CI.get_dpos(c,pos);
      if (pos[XX]<0.5) {
        c->P[RO] = 1.0;
        c->P[SimPM->ftr] = 0.0;
      }
      else {
        c->P[RO] = 10.0;
        c->P[SimPM->ftr] = 1.0;
      }

      c->P[PG] = 10.0;

      c->P[VX] = Vel[XX];
      c->P[VY] = Vel[YY];
      c->P[VZ] = 0.0;

    } while ( (c=ggg->NextPt(c)) !=0);
  }

  else if (test_id==2) {
    do {
      CI.get_dpos(c,pos);
      if (pos[XX]<0.5) {
        c->P[RO] = 1.0;
        c->P[SimPM->ftr] = 0.0;
      }
      else {
        c->P[RO] = 1000.0;
        c->P[SimPM->ftr] = 1.0;
      }

      c->P[PG] = 1000.0;

      c->P[VX] = Vel[XX];
      c->P[VY] = Vel[YY];
      c->P[VZ] = 0.0;

    } while ( (c=ggg->NextPt(c)) !=0);
  }


  else if (test_id==3 || test_id==4) {
    //
    // pre-calculate trigonometric quantities
    //
    double theta=1.0;
    double
      tt   = tan(theta),
      itt  = 1.0/tt,
      ifst = 1.0/(4.0*sin(theta));
    //
    // set density/pressure inside square
    //
    double rho_in=0.0, pg=0.0;
    if (test_id==3) {
      rho_in = 10.0;
      pg  = 10.0;
    }
    else {
      rho_in = 1000.0;
      pg  = 1000.0;
    }
    //
    // Now go through grid setting inside/outside values.
    //
    do {
      CI.get_dpos(c,pos);
      //
      // Inside the square is non-trivial to calculate.  We use the
      // equation of each boundary line and make sure we are the 
      // correct side of it, and if all inequalities hold then we are
      // inside the square.
      //
      bool inside = true;
      // edge 1, (+,+) direction from 1,1
      if (pos[YY] > 1.0 + itt + ifst -pos[XX]*itt) inside=false;
      // edge 3, (-,-) direction from 1,1
      if (pos[YY] < 1.0 + itt - ifst -pos[XX]*itt) inside=false;
      // edge 2, (-,+) direction from 1,1
      if (pos[YY] > tt*(pos[XX] - (1.0 - itt - ifst))) inside=false;
      // edge 4, (+,-) direction from 1,1
      if (pos[YY] < tt*(pos[XX] - (1.0 - itt + ifst))) inside=false;
      //
      // now inside==true only if we are inside the square.
      //
      if (inside) {
        c->P[RO] = rho_in;
        c->P[SimPM->ftr] = 1.0;
      }
      else {
        c->P[RO] = 1.0;
        c->P[SimPM->ftr] = 0.0;
      }
      c->P[PG] = pg;
      c->P[VX] = Vel[XX];
      c->P[VY] = Vel[YY];
      c->P[VZ] = 0.0;
    } while ( (c=ggg->NextPt(c)) !=0);
  }

  return 0;
}





// ##################################################################
// ##################################################################


int IC_StarBench_Tests::setup_StarBench_planarIF(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  //
  // ----------------------------------------------------------------
  // First get input parameters:
  // v_up, rho_up, T_neutral, T_ionized, v_down
  //
  double rho0=0.0, vel0=0.0, vel2=0.0;

  string seek = rrp->find_parameter("StarBench_IFI_rho0");
  if (seek=="") rep.error("Need parameter StarBench_IFI_rho0",1);
  else rho0 = atof(seek.c_str());

  seek = rrp->find_parameter("StarBench_IFI_vel0");
  if (seek=="") rep.error("Need parameter StarBench_IFI_vel0",1);
  else vel0 = atof(seek.c_str());

  seek = rrp->find_parameter("StarBench_IFI_vel2");
  if (seek=="") rep.error("Need parameter StarBench_IFI_vel2",1);
  else vel2 = atof(seek.c_str());

  // v_x is the velocity flowing into the shock, in the shock frame.
  double v_x = vel0;
  //
  // X_up are upstream parameters, in grid reference frame.
  double d_up = rho0, v_up = 0.0;
  //
  // X_sh are shocked neutral shell parameters, in grid ref. frame.
  double d_sh = 0.0, v_sh = 0.0;
  //
  // X_dn are downstream ionized gas parameters, in grid ref. frame.
  double d_dn = 0.0, v_dn = vel2;
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // c_n and c_i are the neutral and ionized gas sound speeds.
  // Set them based on MinTemperature and MaxTemperature.
  //
  double c_n = 0.0, c_i = 0.0;

  cell *c=ggg->FirstPt();
  c->P[RO] = d_up;
  c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;

  //
  // Neutral gas sound speed: MP->Set_Temp() sets the gas pressure to
  // correspond to the given temperature, then we can get sound speed
  // from c^2 = p/rho
  //
  c->P[PG] = 1.0e-10;
  for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 0.0;
  MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
  c_n = sqrt(c->P[PG]/c->P[RO]);

  //
  // Ionized gas
  // 
  for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
  MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
  c_i = sqrt(c->P[PG]/c->P[RO]);
  
  cout <<"c_n = "<< c_n <<", c_i = "<<c_i<<"\n";
  //
  // Still need to calculate v_up, v_sh, d_sh, d_dn
  // ----------------------------------------------------------------


  // ----------------------------------------------------------------
  // Get shell density from shock jump conditions:
  // rho_shell = rho_0 *M^2
  //
  d_sh = d_up * pow(v_x/c_n, 2);
  //
  // Still need to calculate v_up, v_sh, d_dn
  // ----------------------------------------------------------------
  
  // ----------------------------------------------------------------
  // Get v_shell from downstream exhaust velocity and the jump 
  // conditions.
  // First calculate the term in square root (b^2-4ac)
  //
  v_sh = pow(v_dn,4) + 2.0*pow(c_i*v_dn,2) +pow(c_i,4)
         - 4.0*pow(c_n*v_dn,2);
  //
  // Now calculate v_sh from the quadratic root formula.
  //
  v_sh = (0.5/v_dn) * (pow(v_dn,2) + pow(c_i,2) - sqrt(v_sh));
  //
  // Still need to calculate v_up, d_dn
  // ----------------------------------------------------------------
  
  // ----------------------------------------------------------------
  // Get v_up from shock velocity and v_x
  // v_shock = c_n^2/v_x - v_sh
  //
  v_up = pow(c_n,2)/v_x - v_sh;
  cout  <<"shock vel = "<<v_up;
  v_up = v_x - v_up;
  cout <<" : upstream vel = "<<v_up<<"\n";
  //
  // Still need to calculate d_dn
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // Get d_dn from continuity
  d_dn = d_sh*v_sh/v_dn;
  // Now we have the parameters for the three regions.
  cout <<"\t v_0 = "<<v_up<<", rho_0 = "<<d_up<<"\n";
  cout <<"\t v_1 = "<<v_sh<<", rho_1 = "<<d_sh<<"\n";
  cout <<"\t v_2 = "<<v_dn<<", rho_2 = "<<d_dn<<"\n";
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // Get the position of the I-front, assuming the flow is from right
  // to left, and the ionizing source is shining from the left.
  //
  c->P[RO] = d_dn;
  c->P[VX] = v_dn;
  c->P[PG] = 0.0;
  MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
  for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
  
  //
  // This is the local recombination rate (number per cm^3 per sec).
  //
  double IF_pos = MP->get_recombination_rate(0, c->P, SimPM->gamma);
  //
  // Length to absorb all photons is the source strength * IF_pos
  //
  if (SimPM->RS.Nsources<1) {
    rep.error("Need Radiation Source for test!",SimPM->RS.Nsources);
  }
  cout <<"Recombination rate = "<<IF_pos<<" :  ";
  IF_pos = SimPM->RS.sources[0].strength / IF_pos;

  cout <<"source ionizing flux = "<< SimPM->RS.sources[0].strength;
  cout <<" :  "<<"length to absorb photons = "<<IF_pos<<"\n";
  //
  // Add in a fudge factor of 0.65 because the downstream region is
  // not constant density but increases towards the I-front.  So we
  // get more recombinations than for constant density.
  //IF_pos = 0.65*IF_pos + ggg->SIM_Xmin(XX);
  IF_pos = sqrt(v_x/c_i)*IF_pos + ggg->SIM_Xmin(XX);

  //
  // Use given parameter, if present, for initial position of IF.
  //
  double tmp = IF_pos;
  seek = rrp->find_parameter("StarBench_IFI_xIF");
  if (seek!="") {
    IF_pos = atof(seek.c_str());
    cout <<"Using given IF_pos="<<IF_pos<<" cm.  Calculated value="<<tmp<<" cm\n";
  }
  else {
    cout <<"No IF_pos specified.  Using calculated IF_pos="<<IF_pos<<" cm.\n";
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // Now set the positions of the two discontinuities.
  // Shell Thickness is input in number of cells.
  double shock_pos = 0.0;
  seek = rrp->find_parameter("StarBench_IFI_shell_thickness");
  if (seek=="") rep.error("Need parameter StarBench_IFI_shell_thickness",1);
  else shock_pos = atof(seek.c_str())*ggg->DX() +IF_pos;
  cout <<"IF_pos = "<<IF_pos<<", shock_pos="<<shock_pos<<"\n";
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  double pos[SimPM->ndim];
  do {
    CI.get_dpos(c,pos);

    if (pos[XX]<= IF_pos) {
      //
      // Set to downstream properties, ionized.
      //
      c->P[RO] = d_dn;
      c->P[PG] = 1.0e-10;
      c->P[VX] = -v_dn;
      c->P[VY] = c->P[VZ] = 0.0;
      for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
      MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
    }
    else if (pos[XX]<= shock_pos) {
      //
      // Shell properties, neutral.
      //
      c->P[RO] = d_sh;
      c->P[PG] = 1.0e-10;
      c->P[VX] = -v_sh;
      c->P[VY] = c->P[VZ] = 0.0;
      for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0e-12;
      MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
    }
    else {
      //
      // Upstream unshocked properties, neutral.
      //
      c->P[RO] = d_up;
      c->P[PG] = 1.0e-10;
      c->P[VX] = -v_up;
      c->P[VY] = c->P[VZ] = 0.0;
      for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0e-12;
      MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
    }

  } while ( (c=ggg->NextPt(c)) !=0);
  // ----------------------------------------------------------------


  //
  // ----------------------------------------------------------------
  // Add a perturbation upstream:
  //
  int ptype = -1;
  seek = rrp->find_parameter("StarBench_IFI_perturbation");
  if      (seek=="1" || seek=="velocity") {
    ptype = 1;
  }
  else if (seek=="2" || seek=="deformation") {
    ptype = 2;
  }
  else if (seek=="3" || seek=="def_small") {
    ptype = 3;
  }
  else if (seek=="4" || seek=="density") {
    ptype = 4;
  }
  else if (seek=="0" || seek=="none") {
    ptype = 0;
  }
  else {
   rep.error("Need parameter StarBench_IFI_perturbation",seek);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  if (ptype==1) {
    //
    // upstream velocity perturbation, 0.75 of the neutral gas sound
    // speed.
    //
    double lambda = 0.125*SimPM->Range[YY];
    double A = 0.75*c_n;
    double x0 = shock_pos + 0.2*SimPM->Range[XX];
    double sig= 0.05*SimPM->Range[XX];
    c=ggg->FirstPt();
    do {
      CI.get_dpos(c,pos);
      c->P[VY] = A*sin(2.0*M_PI*(pos[YY]+0.5*SimPM->Range[YY])/lambda)
                      *exp(-0.5*pow((pos[XX]-x0)/sig,2.0));
    } while ( (c=ggg->NextPt(c)) !=0);
  } //  ptype==1
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  else if (ptype==2 || ptype==3) {
    //
    // Overwrite data, with curved shock.
    //
    double lambda = SimPM->Range[YY];
    if      (ptype==2) lambda *= 0.25;  // 4 wavelengths on the domain
    else if (ptype==3) lambda *= 1.00;  // just one wavelength

    double A = 0.0;
    if      (ptype==2) A = lambda / 8.0;  // 1/32 of the y-domain, 1/8 of l.
    else if (ptype==3) A = lambda /128.0;  // 1/128 of l and the y-domain (1 cell at r1)
    else               A = 0.0;
    double deflection=0.0;
    c=ggg->FirstPt();
    do {
      CI.get_dpos(c,pos);
      deflection = A*sin(2.0*M_PI*(pos[YY]+0.5*SimPM->Range[YY])/lambda);
      //if (pos[XX]<= IF_pos+deflection) {
      // This makes the IF smooth, but the SF is still distorted.
      if (pos[XX]<= IF_pos) {
        //
        // Set to downstream properties, ionized.
        //
        c->P[RO] = d_dn;
        c->P[PG] = 1.0e-10;
        c->P[VX] = -v_dn;
        c->P[VY] = c->P[VZ] = 0.0;
        for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
        MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
      }
      else if (pos[XX]<= shock_pos+deflection) {
        //
        // Shell properties, neutral.
        //
        c->P[RO] = d_sh;
        c->P[PG] = 1.0e-10;
        c->P[VX] = -v_sh;
        c->P[VY] = c->P[VZ] = 0.0;
        for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0e-12;
        MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
      }
      else {
        //
        // Upstream unshocked properties, neutral.
        //
        c->P[RO] = d_up;
        c->P[PG] = 1.0e-10;
        c->P[VX] = -v_up;
        c->P[VY] = c->P[VZ] = 0.0;
        for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0e-12;
        MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
      }

    } while ( (c=ggg->NextPt(c)) !=0);
  } //  ptype==2 || 3
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  else if (ptype==4) {
    c=ggg->FirstPt();
    c->P[RO] = d_dn;
    c->P[VX] = v_dn;
    c->P[PG] = 0.0;
    MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
    for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
    IF_pos = MP->get_recombination_rate(0, c->P, SimPM->gamma);
    IF_pos = 0.65*SimPM->RS.sources[0].strength / IF_pos + ggg->SIM_Xmin(XX);
    //
    // Density sinusoidal perturbation, gaussian smoothed in x.
    //
    double centre = ggg->SIM_Xmin(XX)+ 0.75*ggg->SIM_Range(XX);
    double sigma  = 0.05*ggg->SIM_Range(XX);
    double lambda = ggg->SIM_Range(YY);  // perturbation wavelength
    double amp    = 0.1;  // relative amplitude of density perturbation
    double deltarho = 0.0; // overdensity (calculated at each grid cell)
    //
    // Put data on grid:
    //
    do {
      CI.get_dpos(c,pos);
      if (pos[XX]<= IF_pos) {
        //
        // Set to downstream properties, ionized.
        //
        c->P[RO] = d_dn;
        c->P[PG] = 1.0e-10;
        c->P[VX] = -v_dn;
        c->P[VY] = c->P[VZ] = 0.0;
        for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0;
        MP->Set_Temp(c->P, SimPM->EP.MaxTemperature, SimPM->gamma);
      }
      else {
        deltarho = amp*sin(2.0*M_PI*(pos[YY]+0.5*ggg->SIM_Range(YY))/lambda);
        deltarho *= exp(- 0.5*pow((pos[XX] - centre)/sigma, 2.0));
        //
        // Upstream unshocked properties, neutral.
        //
        c->P[RO] = d_up * (1.0 + deltarho);
        c->P[PG] = 1.0e-10; // reset later with Set_Temp() call.
        c->P[VX] = -v_up;
        c->P[VY] = c->P[VZ] = 0.0;
        for (int v=0; v<SimPM->ntracer; v++) c->P[SimPM->ftr+v] = 1.0e-12;
        MP->Set_Temp(c->P, SimPM->EP.MinTemperature, SimPM->gamma);
      }
    } while ( (c=ggg->NextPt(c)) !=0);

  } //  ptype==4
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  else if (ptype==0) {
  } //  ptype==0
  // ----------------------------------------------------------------



  return 0;
}




// ##################################################################
// ##################################################################


int IC_StarBench_Tests::setup_StarBench_IFI(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  
  int id=0;
  if      (test=="StarBench_IFI_testA") id=1;
  else if (test=="StarBench_IFI_testB") id=2;
  else if (test=="StarBench_IFI_testC") id=3;
  else rep.error("Bad test name",test);

  cell *c=ggg->FirstPt();
  double pos[SimPM->ndim];
  do {
    //CI.get_dpos(c,pos);

    c->P[RO] = 44.0*pconst.m_p();      // Pure H with n=44/cm3.
    c->P[PG] = 44.0*pconst.kB()*10.0; // 10K neutral gas (pure H).
    c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;
    for (int v=0; v<SimPM->ntracer; v++)
      c->P[SimPM->ftr+v] = 0.0;

  } while ( (c=ggg->NextPt(c)) !=0);

  if (id==3) {
    double lambda = 0.125*SimPM->Range[YY];
    double A = 0.75 *sqrt(pconst.kB()*1.0e4/pconst.m_p());
    double x0 = SimPM->Xmin[XX] +0.12*SimPM->Range[XX];
    double sig= 0.05*SimPM->Range[XX];
    c=ggg->FirstPt();
    do {
      CI.get_dpos(c,pos);
      c->P[VY] = A*sin(2.0*M_PI*(pos[YY]+0.5*SimPM->Range[YY])/lambda)
                  *exp(-0.5*pow((pos[XX]-x0)/sig,2.0));
    } while ( (c=ggg->NextPt(c)) !=0);
  }

  return 0;
}


// ##################################################################
// ##################################################################

int IC_StarBench_Tests::setup_StarBench_IrrCl(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  
  int id=0;
  if      (test=="StarBench_IrrCloud_Uniform") id=1;
  else if (test=="StarBench_IrrCloud_IsoSph")  id=2;
  else rep.error("Bad test name",test);

  cell *c=ggg->FirstPt();
  double pos[SimPM->ndim];
  //
  // Run through grid and set uniform initial conditions.
  //
  do {
    //CI.get_dpos(c,pos);

    c->P[RO] = 50.0*pconst.m_p();      // From the test document.
    c->P[PG] = 50.0*pconst.kB()*1000.0; // Just pick some temperature...
    c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;
    for (int v=0; v<SimPM->ntracer; v++)
      c->P[SimPM->ftr+v] = 0.0;
  } while ( (c=ggg->NextPt(c)) !=0);

  //
  // Now run through again and add the cloud.
  //
  if (id==1) {
    cout <<"\t\tAdding uniform-density cloud.\n";
    //
    // add a uniform density cloud with radius 1pc, centred
    // at x=1.94 pc, y=z=0.  rho_cloud=1000*m_p
    //
    double radius = 3.086e18;
    double dist=0.0;
    double rho_cl = 1000.0*pconst.m_p();
    double cl_centre[SimPM->ndim];
    cl_centre[XX] = 1.92*3.086e18;
    for (int v=1;v<SimPM->ndim;v++) cl_centre[v]=0.0;

    c=ggg->FirstPt();
    do {
      CI.get_dpos(c,pos);
      dist = gg->distance(cl_centre,pos);
      if (dist<radius)
        c->P[RO] = rho_cl;
    } while ( (c=ggg->NextPt(c)) !=0);
  }

  else if (id==2) {
    cout <<"\t\tAdding cutoff-isothermal-sphere cloud.\n";
    //
    // Add a cutoff isothermal sphere cloud with central density of
    // rho_cloud=1000*m_p, and core radius r_c=0.5pc.
    //
    double r_core = 0.5*3.086e18;
    double rho_cl = 1000.0*pconst.m_p();
    double rho_cell=0.0;
    double cl_centre[SimPM->ndim];
    cl_centre[XX] = 1.92*3.086e18;
    for (int v=1;v<SimPM->ndim;v++) cl_centre[v]=0.0;
    double dist=0.0;

    c=ggg->FirstPt();
    do {
      //
      // calculate core density, and take max. of this value and ISM.
      //
      CI.get_dpos(c,pos);
      dist = gg->distance(cl_centre,pos);
      rho_cell = rho_cl*r_core*r_core/(r_core*r_core +dist*dist);
      c->P[RO] = std::max(c->P[RO],static_cast<pion_flt>(rho_cell));
    } while ( (c=ggg->NextPt(c)) !=0);
  }
  else rep.error("Bad test name",test);

  return 0;
}


// ##################################################################
// ##################################################################


int IC_StarBench_Tests::setup_StarBench_TremblinCooling(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  //
  // Set density based on parameter (giving n(H) in cm^{-3})
  //
  double density=0.0;
  string seek = rrp->find_parameter("StarBench_TremblinCooling_Rho");
  if (seek=="") rep.error("Need parameter StarBench_TremblinCooling_Rho",1);
  else density = atof(seek.c_str());


  double dpos[SimPM->ndim];
  cell *c=ggg->FirstPt();
  do {
    CI.get_dpos(c,dpos);
    
    //
    // set general ISM density, pressure, and ion fraction.
    //
    c->P[RO] = density*pconst.m_p();      // Pure H with n=0.5/cm3.
    c->P[PG] = 2.0*c->P[RO]*pconst.kB()*1.0e4/pconst.m_p(); // 10000K ionised gas (pure H).
    c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;
    for (int v=0; v<SimPM->ntracer; v++)
      c->P[SimPM->ftr+v] = 1.0; // fully ionised


  } while ( (c=ggg->NextPt(c)) !=0);

  return 0;
}


// ##################################################################
// ##################################################################



// ##################################################################
// ##################################################################


int IC_StarBench_Tests::setup_StarBench_Cone(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  
  double srcpos[SimPM->ndim];
  for (int v=0;v<SimPM->ndim;v++) {
    srcpos[v] = SimPM->RS.sources[0].pos[v];
    if (isnan(srcpos[v])) rep.error("Bad source position",srcpos[v]);
  }
  rep.printVec("source position",srcpos,      SimPM->ndim);
  
  //double density=0.0;
  //string seek = rrp->find_parameter("StarBench_Cone_Rho");
  //if (seek=="") rep.error("Need parameter StarBench_TremblinCooling_Rho",1);
  //else density = atof(seek.c_str());
  
  cell *c=ggg->FirstPt();
  double pos[SimPM->ndim];
  double dist = 0.0;
  double r0 = 3.086e17;    // core radius of cloud.
  double radial_slope=2.0; // power law slope in rho for r>r0
  double theta=0.0;
  do {
    CI.get_dpos(c,pos);
    theta = atan2(pos[Rcyl],pos[Zcyl]);

    c->P[RO] = 1.0e4*pconst.m_p(); // Pure H with n=10^4/cm3.
    c->P[PG] = 1.518e-10; // 100K neutral gas (pure H).
    c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;
    for (int v=0; v<SimPM->ntracer; v++)
      c->P[SimPM->ftr+v] = 1.0e-12;
    
    dist = ggg->distance_vertex2cell(srcpos,c);
    //
    // Following the Iliev et al 2009 test 6, we use rho=rho0(r0/r)^n if r>r0
    // We also change the pressure so there is a constant temperature state.
    //
    if (dist>r0) {
      c->P[RO] *= exp(radial_slope*log(r0/dist))*(1.0-0.25*cos(theta));
      c->P[PG] *= exp(radial_slope*log(r0/dist))*(1.0-0.25*cos(theta));
    }
  

  } while ( (c=ggg->NextPt(c)) !=0);


  return 0;
}

// ##################################################################
// ##################################################################



