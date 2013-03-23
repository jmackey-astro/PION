/// \file StarBench_test.cpp
/// \author Jonathan Mackey
///
/// File for setting up initial conditions for test problems to be
/// run as part of the StarBench code-testing workshop.
///
/// - 2013.01.10 JM: Started on file.
///     Added Contact Discontinuity tests.

#include "ics/icgen.h"
#include "coord_sys/VectorOps.h"
#include "dataIO/dataio.h"
//#include "dataIO/dataio_fits.h"
#include <sstream>

IC_StarBench_Tests::IC_StarBench_Tests()
{
  return;
}

IC_StarBench_Tests::~IC_StarBench_Tests()
{}


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

  if (SimPM.eqntype != EQEUL)
    rep.error("Bad equations",SimPM.eqntype);

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
  else if (ics=="StarBench_IFI_TestA" ||
           ics=="StarBench_IFI_TestB" ||
           ics=="StarBench_IFI_TestC") {
    cout <<"\t\tSetting up StarBench Planar ionisation front A.\n";
    err += setup_StarBench_IFI(rrp,ggg,ics);
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
    cout <<"\t\tNOISE: Adding random adiabatic noise at fractional level = "<<noise<<endl;
    err+= AddNoise2Data(2,noise);
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
  if (test_id <=2  && SimPM.ndim != 1) {
    rep.error("ContactDiscontinuity1/2 is 1D test",SimPM.ndim);
  }
  if (test_id >=3 && SimPM.ndim != 2) {
    rep.error("ContactDiscontinuity3/4 is 2D test",SimPM.ndim);
  }

  if (SimPM.ntracer != 1) {
    rep.error("Need exactly one colour tracer for CD tests",
              SimPM.ntracer);
  }
  
  //
  // Get advection velocity from the parameterfile.
  //
  double Vel[SimPM.ndim];
  cout << "** It is assumed all parameters are correctly set **\n";

  string seek = rrp->find_parameter("StarBench_ContDisc_VX");
  if (seek=="") rep.error("Need parameter StarBench_ContDisc_VX",1);
  else Vel[XX] = atof(seek.c_str());

  if (SimPM.ndim>1) {
    string seek = rrp->find_parameter("StarBench_ContDisc_VY");
    if (seek=="") rep.error("Need parameter StarBench_ContDisc_VY",2);
    else Vel[YY] = atof(seek.c_str());
  }

  cell *c=ggg->FirstPt();
  double pos[SimPM.ndim];

  if      (test_id==1) {
    do {
      CI.get_dpos(c,pos);
      if (pos[XX]<0.5) {
        c->P[RO] = 1.0;
        c->P[SimPM.ftr] = 0.0;
      }
      else {
        c->P[RO] = 10.0;
        c->P[SimPM.ftr] = 1.0;
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
        c->P[SimPM.ftr] = 0.0;
      }
      else {
        c->P[RO] = 1000.0;
        c->P[SimPM.ftr] = 1.0;
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
        c->P[SimPM.ftr] = 1.0;
      }
      else {
        c->P[RO] = 1.0;
        c->P[SimPM.ftr] = 0.0;
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


int IC_StarBench_Tests::setup_StarBench_IFI(
      class ReadParams *rrp,    ///< pointer to parameter list.
      class GridBaseClass *ggg, ///< pointer to grid
      string &test
      )
{
  
  int id=0;
  if      (test=="StarBench_IFI_TestA") id=1;
  else if (test=="StarBench_IFI_TestB") id=2;
  else if (test=="StarBench_IFI_TestC") id=3;
  else rep.error("Bad test name",test);

  cell *c=ggg->FirstPt();
  double pos[SimPM.ndim];
  do {
    //CI.get_dpos(c,pos);

    c->P[RO] = 44.0*GS.m_p();      // Pure H with n=44/cm3.
    c->P[PG] = 44.0*GS.kB()*10.0; // 10K neutral gas (pure H).
    c->P[VX] = c->P[VY] = c->P[VZ] = 0.0;
    for (int v=0; v<SimPM.ntracer; v++)
      c->P[SimPM.ftr+v] = 0.0;

  } while ( (c=ggg->NextPt(c)) !=0);

  if (id==3) {
    double lambda = 0.125*SimPM.Range[YY];
    double A = 0.75 *sqrt(GS.kB()*1.0e4/GS.m_p());
    c=ggg->FirstPt();
    do {
      CI.get_dpos(c,pos);
      c->P[VY] = A*sin(2.0*M_PI*(pos[YY]+0.5*SimPM.Range[YY])/lambda)
                  *exp(-4.0*pow((pos[XX]-SimPM.Xmin[XX])/SimPM.Range[XX],2.0));
    } while ( (c=ggg->NextPt(c)) !=0);
  }

  return 0;
}


// ##################################################################
// ##################################################################






