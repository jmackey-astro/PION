///
/// \file integrator.cc
/// \author Jonathan Mackey
///
/// This is a general purpose integrator with a few different methods.
/// The most useful are brute force Euler integration and a 5th order
/// Runge Kutta method with adaptive stepsize (from Numerical Recipes).
///
/// Mods:
///  - Pre-2009-12 changes are undocumented except in CVS/SVN logs.
///  - 2009-12-16/17 JM: Tightened up some error checking in the RK5 methods.
///
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2013.01.11 JM: changed stepper tolerance in adaptive substepper
///    from 20 to 50 iterations.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "microphysics/integrator.h"
#include <iostream>
using namespace std;

// ##################################################################
// ##################################################################

Integrator_Base::Integrator_Base()
{
  // parameters are from Numerical Recipes in C (1992), p717, sec.16.2.
  a2 = 0.2;
  a3 = 0.3;
  a4 = 0.6;
  a5 = 1.0;
  a6 = 0.875;

  b21 = 0.2;

  b31 = 3. / 40.;
  b32 = 9. / 40.;

  b41 = 0.3;
  b42 = -0.9;
  b43 = 1.2;

  b51 = -11. / 54.;
  b52 = 2.5;
  b53 = -70. / 27.;
  b54 = 35. / 27.;

  b61 = 1631. / 55296.;
  b62 = 175. / 512.;
  b63 = 575. / 13824.;
  b64 = 44275. / 110592.;
  b65 = 253. / 4096.;

  c1 = 37. / 378.;
  c3 = 250. / 621.;
  c4 = 125. / 594.;
  c6 = 512. / 1771.;  ///< c2=c5=0. these are the 5th order weights.

  dc1 = c1 - 2825. / 27648.;
  dc3 = c3 - 18575. / 48384.;
  dc4 = c4 - 13525. / 55296.;
  dc5 = -277. / 14336.;
  dc6 =
      c6 - 0.25;  ///< dc2=0.  These are the 5th order minus 4th order weights.

  k1 = k2 = k3 = k4 = k5 = k6 = 0;
  int_nvar                    = 0;
}

// ##################################################################
// ##################################################################

Integrator_Base::~Integrator_Base()
{
  k1 = mem.myfree(k1);
  k2 = mem.myfree(k2);
  k3 = mem.myfree(k3);
  k4 = mem.myfree(k4);
  k5 = mem.myfree(k5);
  k6 = mem.myfree(k6);
}

// ##################################################################
// ##################################################################

int Integrator_Base::Set_Nvar(int nv)
{
  Integrator_Base::int_nvar = nv;
  k1                        = mem.myfree(k1);
  k2                        = mem.myfree(k2);
  k3                        = mem.myfree(k3);
  k4                        = mem.myfree(k4);
  k5                        = mem.myfree(k5);
  k6                        = mem.myfree(k6);
  k1                        = mem.myalloc(k1, int_nvar);
  k2                        = mem.myalloc(k2, int_nvar);
  k3                        = mem.myalloc(k3, int_nvar);
  k4                        = mem.myalloc(k4, int_nvar);
  k5                        = mem.myalloc(k5, int_nvar);
  k6                        = mem.myalloc(k6, int_nvar);
  return 0;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Int_Euler(
    const int nv,      ///< number of elements in P array.
    const double* p0,  ///< value of P at initial value of t.
    const double t0,   ///< initial value of t.
    const double dt,   ///< Total step dt to take.
    const int nn,      ///< Number of integration steps to take.
    double* pf         ///< pointer to final P value.
)
{
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }
  double h = dt / (static_cast<double>(nn));
  double t = t0;
  double ptemp[int_nvar];
  int err = 0;
  for (int v = 0; v < int_nvar; v++)
    pf[v] = p0[v];

  // simplest 1st order Euler integration p += dpdt*dt
  for (int i = 0; i < nn; i++) {
    err += dPdt(nv, pf, ptemp);
    for (int v = 0; v < int_nvar; v++)
      pf[v] += h * ptemp[v];
    t += h;
  }
  if (!pconst.equalD(t, t0 + dt)) {
    cout.setf(ios_base::scientific, ios_base::floatfield);
    cout.precision(12);
    cout << "t: " << t << " dt: " << dt << " t0: " << t0 << " h: " << h
         << "\t eps: " << (t - t0 - dt) / (t + t0 + dt) << "\n";
    rep.error("Int_Euler coding error, h too small??", t - t0 + dt);
  }

  return err;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Int_DumbAdaptive_Euler(
    const int nv,         ///< number of elements in P array.
    const double* p0,     ///< value of P at initial value of t.
    const double t0,      ///< initial value of t.
    const double dt,      ///< Total step dt to take.
    const double errtol,  ///< error tolerance per step.
    double* pf,           ///< pointer to final P value.
    double* tf            ///< pointer to final t value.
)
{
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    Set_Nvar(nv);
  }
  if (errtol < 0)
    rep.error("Int_DumbAdaptive_Euler() ErrTol is negative!", errtol);
  if (errtol < MACHINEACCURACY)
    rep.error(
        "errtol beyond machine accuracy.  use more lenient value!\n", errtol);
  // Dumb adaptive integrator.  Start with 16 subpoints, keep doubling until
  // values differ by less than errtol
  int nsub     = 16;
  int max_iter = 20;  // should get us to 2^24=1.6e7 points.

  double p1[int_nvar], p2[int_nvar];

  double t = t0;
  *tf      = t + dt;
  double maxerr, tmperr;
  for (int v = 0; v < int_nvar; v++)
    if (fabs(p0[v]) < 1.e-100)
      cerr << "Int_DumbAdaptive_Euler() WARNING: tiny values, so error "
              "unreliable!\n";
  int err  = Int_Euler(int_nvar, p0, t0, dt, nsub, p1);
  int iter = 0;

  do {
    maxerr = 0.0;
    nsub *= 2;
    err += Int_Euler(int_nvar, p0, t0, dt, nsub, p2);
    for (int v = 0; v < int_nvar; v++) {
      // if we get back NAN, make sure we try a shorter step!
      if (isnan(p2[v]) || isnan(p1[v]))
        maxerr = max(maxerr, 1000.0);
      else {
        tmperr = fabs(p2[v] - p1[v]) / (fabs(p2[v]) + fabs(p1[v]) + 1.e-100);
        //	cout <<"tmperr ="<<tmperr<<"\n";
        maxerr = max(maxerr, tmperr);
      }
      p1[v] = p2[v];
    }
    iter++;
  } while ((iter <= max_iter) && (maxerr > errtol));

  for (int v = 0; v < int_nvar; v++)
    pf[v] = p1[v];
  if (iter > max_iter)
    cerr << "Int_DumbAdaptive_Euler() WARNING: not converged!\n";
  //  cout <<"iter="<<iter<<" nsub="<<nsub<<" maxerr="<<maxerr<<"\n";

  return err;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Step_RK4(
    const int nv,      ///< number of elements in P array.
    const double* p0,  ///< value of P at initial value of t.
    const double t0,   ///< initial value of t.
    const double dt,   ///< Total step dt to take.
    double* pf         ///< pointer to final P value.
)
{
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }
  // Fourth Order Runge-Kutta Method.
  // See NR p.711, eqn:16.1.3
  double ptemp[int_nvar];
  int err = 0;

  err += dPdt(int_nvar, p0, k1);
  for (int v = 0; v < int_nvar; v++) {
    k1[v] *= dt;
    ptemp[v] = p0[v] + 0.5 * k1[v];  // temp array as point to get slope at.
  }
  //  rep.printVec("k1",k1,int_nvar);
  // rep.printVec("ptemp",ptemp,int_nvar);
  err += dPdt(int_nvar, ptemp, k2);
  for (int v = 0; v < int_nvar; v++) {
    k2[v] *= dt;
    ptemp[v] = p0[v] + 0.5 * k2[v];  // temp array as point to get slope at.
  }
  err += dPdt(int_nvar, ptemp, k3);
  for (int v = 0; v < int_nvar; v++) {
    k3[v] *= dt;
    ptemp[v] = p0[v] + k3[v];  // temp array as point to get slope at.
  }
  err += dPdt(int_nvar, ptemp, k4);
  for (int v = 0; v < int_nvar; v++)
    k4[v] *= dt;

  for (int v = 0; v < int_nvar; v++)
    pf[v] = p0[v] + (k1[v] + 2.0 * k2[v] + 2.0 * k3[v] + k4[v]) / 6.0;

  return 0;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Step_RK5CK(
    const int nv,      ///< number of elements in P array.
    const double* p0,  ///< value of y at initial value of t
    const double t0,   ///< initial value of t.
    const double dt,   ///< stepsize, dt.
    double* pf,        ///< final value of P
    double* dp         ///< error estimate.
)
{
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }
  /* single variable version...
  double k1,k2,k3,k4,k5,k6;
  k1 = dt *dPdt(nv, P0);
  k2 = dt *dPdt(nv, P0 +b21*k1);
  k3 = dt *dPdt(nv, P0 +b31*k1 +b32*k2);
  k4 = dt *dPdt(nv, P0 +b41*k1 +b42*k2 +b43*k3);
  k5 = dt *dPdt(nv, P0 +b51*k1 +b52*k2 +b53*k3 +b54*k4);
  k6 = dt *dPdt(nv, P0 +b61*k1 +b62*k2 +b63*k3 +b64*k4 +b65*k5);

  // *ans is the 5th order result.
  *ans = P0 + c1*k1         + c3*k3 + c4*k4         + c6*k6;
  // *err is the difference between the 5th and 4th order result.
  *err =     dc1*k1         +dc3*k3 +dc4*k4 +dc5*k5 +dc6*k6;
  */

  double ptemp[int_nvar];

  int err = 0;
  err += dPdt(int_nvar, p0, k1);

  //
  // First see if dPdt is small enough to make it 1st order.
  // calculate sum((dP/dt)*dt/p0) = dlnP/dlnt, and if it's small,
  // don't bother with fifth order, just return 1st order.
  // need fabs(k1[v]) because rates can be negative, but all the
  // actual values are positive definite so no fabs(p0[v]) needed.
  //
  ptemp[0] = 0.0;
  // static int ntimes=0, nnot=0, fll=0;
  //  cout <<"pt0:"<<ptemp[0]<<" dt="<<dt<<"\n";
  //  rep.printVec("p0",p0,int_nvar);
  //  rep.printVec("dp",k1,int_nvar);
  for (int v = 0; v < int_nvar; v++)
    ptemp[0] += fabs(k1[v]) * dt / (p0[v] + 1.0e-100);
  if (ptemp[0] < 1.e-6) {
    // ntimes++;
    for (int v = 0; v < int_nvar; v++) {
      pf[v] = p0[v] + k1[v] * dt;
      dp[v] = k1[v] * dt;  // this is max absolute error
    }
    return err;
  }
  //  else nnot++;

  // now if we get this far, we are doing the 5th order Cash-Karp RK step.
  for (int v = 0; v < int_nvar; v++) {
    ptemp[v] = 0.0;
    k1[v] *= dt;
    ptemp[v] = p0[v] + b21 * k1[v];  // temp array as point to get slope at.
    //    cout <<"ptemp address: index "<<v<<" addr: "<<&ptemp[v]<<"\n";
  }
  //  rep.printVec("k1",k1,int_nvar);
  err += dPdt(int_nvar, ptemp, k2);
  for (int v = 0; v < int_nvar; v++) {
    k2[v] *= dt;
    ptemp[v] = p0[v] + b31 * k1[v] + b32 * k2[v];
  }
  //  rep.printVec("k2",k2,int_nvar);
  err += dPdt(int_nvar, ptemp, k3);
  for (int v = 0; v < int_nvar; v++) {
    k3[v] *= dt;
    ptemp[v] = p0[v] + b41 * k1[v] + b42 * k2[v] + b43 * k3[v];
  }
  //  rep.printVec("k3",k3,int_nvar);
  err += dPdt(int_nvar, ptemp, k4);
  for (int v = 0; v < int_nvar; v++) {
    k4[v] *= dt;
    ptemp[v] = p0[v] + b51 * k1[v] + b52 * k2[v] + b53 * k3[v] + b54 * k4[v];
  }
  //  rep.printVec("k4",k4,int_nvar);
  err += dPdt(int_nvar, ptemp, k5);
  for (int v = 0; v < int_nvar; v++) {
    k5[v] *= dt;
    ptemp[v] = p0[v] + b61 * k1[v] + b62 * k2[v] + b63 * k3[v] + b64 * k4[v]
               + b65 * k5[v];
  }
  //  rep.printVec("k5",k5,int_nvar);
  err += dPdt(int_nvar, ptemp, k6);
  for (int v = 0; v < int_nvar; v++) {
    k6[v] *= dt;
  }
  //  rep.printVec("k6",k6,int_nvar);

  for (int v = 0; v < int_nvar; v++) {
    // pf is the 5th order result.
    pf[v] = p0[v] + c1 * k1[v] + c3 * k3[v] + c4 * k4[v] + c6 * k6[v];
    // err is the difference between the 5th and 4th order result.
    dp[v] = dc1 * k1[v] + dc3 * k3[v] + dc4 * k4[v] + dc5 * k5[v] + dc6 * k6[v];
  }
  //  rep.printVec("pf", pf, int_nvar);
  //  rep.printVec("dp", dp, int_nvar);

  return err;
}

// ##################################################################
// ##################################################################

#define BISECTION_STEPPER
int Integrator_Base::Stepper_RKCK(
    const int nv,  ///< number of elements in P array.
    const double* p0,
    const double t0,
    const double htry,
    const double errtol,
    double* p1,
    double* hdid,
    double* hnext)
{
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }
  int rval = 0;
  double h = htry;
  if (h < 0) {
    cout << "Integrator_Base::Stepper_RKCK() positive stepsize please\n";
    return 1;
  }
  double tnew = t0;
  double eps  = 1.e-100;  // tiny number
#ifndef BISECTION_STEPPER
  double htemp;
#endif
  double maxerr;
  int ct = 0;

  double err[int_nvar], ptemp[int_nvar];

  for (int v = 0; v < int_nvar; v++) {
    err[v] = 0.0;
    // if (fabs(p0[v]/eps) <10.0) {
    //  cerr <<"WARNING: stepper_RKCK() encountered tiny values, accuracy is
    //  unreliable.\n";
    //      return 1;
    //}
  }

  do {
    // cout <<"\tSTEPPER:\t iter="<<ct<<"\th="<<h<<" and dt = "<<htry<<"\n";
    // cout <<"\tSTEPPER:\t"; rep.printVec("stepper p0",p0,nv);
    // cout <<"DDDD: "<<p0[0]<<"\n";
    // cout <<"DDDD: "<<p0[1]<<"\n";
    // cout <<"DDDD: "<<p0[2]<<"\n";
    // cout <<"DDDD: "<<p0[3]<<"\n";
    // cout <<"DDDD: "<<p0[4]<<"\n";
    rval += Step_RK5CK(nv, p0, t0, h, ptemp, err);  // returns absolute error
    // cout <<"\tSTEPPER:\t"; rep.printVec("stepper pt",ptemp,nv);
    maxerr = 0;
    for (int v = 0; v < int_nvar; v++) {
      // if we get back NAN or negative values, try a shorter step.
      if (!isfinite(err[v]) || !isfinite(ptemp[v]) || ptemp[v] < 0.0) {
        maxerr = max(maxerr, 1000.0);
        // cout <<"err["<<v<<"]="<<err[v]<<" and
        // ptemp["<<v<<"]="<<ptemp[v]<<"\n";
      }
      else {
        err[v] /= fabs(ptemp[v]) + eps;  // converts to relative error.
        err[v] = fabs(err[v] / errtol);
        maxerr = max(maxerr, err[v]);
      }
    }

    //
    // Shrink Step for next try
    //
#ifdef BISECTION_STEPPER
    if (maxerr > 1.) {
      h /= 2.0;  // Shrink step in half.
    }
#else
    if (maxerr > 1.) {
      htemp = 0.9 * h * exp(-0.25 * log(maxerr));  // NR stepsize shrinker.
      // cout <<"\tSTEPPER:\t shrinking stepsize: h_old = "<<h<<" and
      // h_new =
      // "<<htemp<<"\n";
      h = max(htemp, 0.1 * h);  // max factor of 10 reduction in stepsize.
    }
#endif

    tnew = t0 + h;
    if (tnew == t0) {
      cout << "Stepsize too small in stepper_RKCK()\n";
      rep.printVec("Rel.err.", err, int_nvar);
      rep.printVec("p0      ", p0, int_nvar);
      rep.printVec("ptemp   ", ptemp, int_nvar);
      return -2;
    }
    ct++;
  } while (maxerr > 1.0 && ct < 50);

  if (maxerr > 1.0) {
    cout << "stepper_RKCK() has large error estimate: ct=" << ct
         << " out of max. 50 iterations. rel.err.=" << maxerr * errtol << "\n";
    rval += ct + static_cast<int>(fabs(maxerr));
  }

//#define BISECTION_STEPPER
#ifdef BISECTION_STEPPER
  *hnext = h * 2.0;  // Try double the step for the next one.
#else
  htemp = 0.9 * h * exp(-0.2 * log(maxerr));
  *hnext = min(htemp, 5.0 * h);
#endif

  //
  // Copy result vector into output vector, and check is has real numbers in
  // it!
  //
  *hdid = h;
  for (int v = 0; v < int_nvar; v++) {
    p1[v] = ptemp[v];
    if (isnan(p1[v]) || isinf(p1[v])) {
#ifdef TESTING
      commandline.console("Error >");
#endif
      cerr << "\tSTEPPER:\t NANs encountered!\n";
      rep.printVec("Rel.err.", err, int_nvar);
      rep.printVec("p0      ", p0, int_nvar);
      rep.printVec("ptemp   ", ptemp, int_nvar);
      p1[v] = -1.e100;
      rval++;
    }
  }

#ifdef TESTING
  if (dp.c->id == 2773) {
    cout << "*********RKCK err estimate = " << maxerr * errtol << " after "
         << ct << " RKCK calls. rval=" << rval << " and dt=" << h << "\n";
    rep.printVec("p0      ", p0, int_nvar);
    rep.printVec("ptemp   ", ptemp, int_nvar);
  }
#endif

  return rval;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Int_Adaptive_RKCK(
    const int nv,         ///< number of elements in P array.
    const double* p0,     ///< initial state vector.
    const double t0,      ///< initial time
    const double dt,      ///< timestep to advance by.
    const double errtol,  ///< error tolerance per step.
    double* pf,           ///< final state vector
    double* tf            ///< pointer to final time.
)
{
  // cout <<"\t\t\tIntegrator: explicit step!\n";
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    Set_Nvar(nv);
  }
  if (errtol < 0) rep.error("Int_Adaptive_RKCK() ErrTol is negative!", errtol);
  if (errtol < MACHINEACCURACY)
    rep.error(
        "errtol beyond machine accuracy.  use more lenient value!\n", errtol);
  if (errtol > 2)
    rep.error(
        "errtol is too large: relative accuracy required is >1!\n", errtol);
  double t = t0;

  double p1[int_nvar];
  double p2[int_nvar];

  for (int v = 0; v < int_nvar; v++)
    p1[v] = p0[v];

  *tf      = t0 + dt;
  double h = dt, hdid = 0.0, hnext = 0.0;
  int err = 0;
  int ct = 0, ctmax = 25;

  do {
    // rep.printVec("adaptive p1",p1,nv);
    err += Stepper_RKCK(nv, p1, t, h, errtol, p2, &hdid, &hnext);
    t += hdid;
    h = min(hnext, *tf - t);
    ct++;
    for (int v = 0; v < int_nvar; v++)
      p1[v] = p2[v];
    // cout <<"ADAPTIVE INT:\t t_old="<<t-hdid<<" t_new="<<t<<"
    // iter="<<ct<<"\n"; cout <<"ADAPTIVE INT:\t ";
    // rep.printVec("pnew",p1,int_nvar);
  } while (t < *tf && (err == 0) && (ct < ctmax));

  if (err || ct > ctmax) {
    if (err) {
      cout << "Int_Adaptive_RKCK() errors encountered. nstep=" << ct
           << " and err=" << err << "\n";
      rep.printVec("p1", p1, int_nvar);
    }
    else {
      err = ct;
      cout << "integration took too many steps!!! ct=err=" << ct << "\t" << err
           << "\n";
      // if (pconst.equalD(t,*tf)) {
      //        cout<<"took too many steps, but got to end of step, so
      //        returning normally!\n";
      //	err=0;
      //}
    }
  }
  // if (ct>10) cout <<"IntBase: ADAPTIVE INT: took "<<ct<<" steps!\n";
  for (int v = 0; v < int_nvar; v++)
    pf[v] = p1[v];
  *tf = t;

  return err;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Int_Iterative_FC(
    const int nv,      ///< number of elements in P array.
    const double* p0,  ///< value of P at initial value of t.
    const double t0,   ///< initial value of t.
    const double dt,   ///< Total step dt to take.
    double* pf         ///< pointer to final P value.
)
{
  /** \section Coefficients
   * This function is explicitly for solving the rate equations
   * pdot=A*(1-p)-B*p We linearise the equation by fixing A,B=constants, with
   * parameters given by some combination of p(0) and p(n-1), and put these
   * into the analytic solution p(n) = A/(A+B) +
   * (p(0)-A/(A+B))exp(-(A+B)(t-t0)) And iterate until p(n)=p(n-1) to some
   * tolerance. We can't just use dPdt() for this, but need the actual
   * parameters A,B, given by the functions C_rate() = A and D_rate()=A+B.  So
   * to use this function, the derived class must have these defined right.
   */
  //  cout <<"iterative_fixedcoeff() starting.\n";
  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }
  int err = 0;
  double errmax;
  double errtol = 1.e-6;
  double eps    = 1.e-12;
  int Nmax      = 100;
  int nc        = 0;

  double p1[int_nvar], p2[int_nvar], peq[int_nvar], ti[int_nvar];

  for (int v = 0; v < int_nvar; v++)
    p2[v] = p0[v];
  do {
    for (int v = 0; v < int_nvar; v++)
      p1[v] = p2[v];
    err += C_rate(nv, p1, peq);
    err += D_rate(nv, p1, ti);
    for (int v = 0; v < int_nvar; v++) {
      ti[v] = 1.0 / ti[v];
      peq[v] *= ti[v];
      p2[v] = peq[v] + (p0[v] - peq[v]) * exp(-dt / ti[v]);
    }
    nc++;
    errmax = 0.0;
    for (int v = 0; v < int_nvar; v++)
      errmax =
          max(errmax, fabs(p1[v] - p2[v]) / (fabs(p1[v]) + fabs(p2[v]) + eps));
  } while ((errmax > errtol) && nc <= Nmax);

  if (nc > Nmax)
    cerr << "\t!!! failed to converge after " << Nmax << " steps. ";
  //  cout <<"\tAccuracy of solution = "<<errmax<<"\n";
  for (int v = 0; v < int_nvar; v++)
    pf[v] = p2[v];

  if (nc > Nmax)
    return nc;
  else
    return err;
}

// ##################################################################
// ##################################################################

int Integrator_Base::Int_Subcycle_BDF(
    const int nv,  ///< number of elements in P array.
    double* p0,    ///< value of P at initial value of t.
    double t0,     ///< initial value of t.
    double dt,     ///< Total step dt to take.
    double dptol,  ///< subcycling parameter: dt_sub = dPtol*P/dPdt
    double* pf     ///< pointer to final P value.
)
{
  int err = 0;

  int Nmax = 100;
  int nc   = 0;

  if (int_nvar != nv) {
    cerr << "Integrator_Base() nvar not equal to state vector length.\n";
    return 1;
  }

  double rate[int_nvar], pnow[int_nvar], pguess[int_nvar];

  double subt;
  double t = 0.0;
  for (int v = 0; v < int_nvar; v++)
    pnow[v] = p0[v];

  do {
    subt = 1.e100;
    err += dPdt(nv, pnow, rate);
    for (int v = 0; v < int_nvar; v++)
      subt = min(subt, pnow[v] / rate[v]);  // this is p/(dpdt)
    subt = min(dt - t, dptol * subt);

    // First get estimate of new value.
    err += C_rate(nv, pnow, rate);  // creation rate [per sec]
    for (int v = 0; v < int_nvar; v++)
      pguess[v] = pnow[v] + rate[v] * subt;
    err += D_rate(nv, pnow, rate);  // destruction rate [per sec]
    for (int v = 0; v < int_nvar; v++)
      pguess[v] /= 1.0 + rate[v] * subt;

    // Now use estimate to get new value via backwards differencing.
    err += C_rate(nv, pguess, rate);  // creation rate [per sec]
    for (int v = 0; v < int_nvar; v++)
      pnow[v] += rate[v] * subt;
    err += D_rate(nv, pguess, rate);  // destruction rate [per sec]
    for (int v = 0; v < int_nvar; v++)
      pnow[v] /= 1.0 + rate[v] * subt;

    // advance a subtimestep and repeat until we get to the end.
    t += subt;
    nc++;
  } while (t < dt && nc <= Nmax);

  //  cout <<"nsteps = "<<nc<<"\n";
  if (nc > Nmax)
    cerr << "BDF_subcycling() took more than " << Nmax << " steps. ";
  for (int v = 0; v < int_nvar; v++)
    pf[v] = pnow[v];

  if (err) cerr << "Int_Subcycle_BDF encountered errors!\n";

  if (nc > Nmax)
    return nc;
  else
    return err;
}

// ##################################################################
// ##################################################################
