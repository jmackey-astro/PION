/// \file microphysics.cc
/// \brief Class Definitions for MicroPhysics Routines.
/// \author Jonathan Mackey
/// 
/// Modifications:
///  - 2007-12-17 Started File with CoolingFn, UpdateMicroPhysics classes.
///  - 2007-12-21 Changed tracers to fractional abundances.
///  - 2008-01-21 Moved Cooling to its own file.  updated interface.
///  - 2008-02-26 Wrote OnlyCooling() functions for TimeUpdateMP() and dPdt()
///  -- stopping putting in mods -- see the CVS log, or my code log.
///
/// - JM 2009-12-16 Fixed microphysics for "only cooling" -- it was
///   broken.  Added in an ifdef in MP_Hydrogen for if the integration
///   (not Photoionisation) returns an ifrac>2, which would indicate
///   something went very wrong.  In that case I try splitting the
///   integral into 10 subintegrals and retrying.
///   If that fails, bug out.
/// - 2010-01-05 JM: Added in (public) function which returns
///    timescales for heating/cooling and recombination/ionisation
///    etc.
/// - 2010-01-15 JM: Changed criteria for setting incoming flux to
///    zero in the RT update for efficiency.  It was failing for large
///    domains, so I tried to test for the value of a
///    scale-independent quantity: photons_in*ds
/// - 2010-01-19 JM: propagated same change from MP_Hydrogen:: into
///    Microphysics:: removed bug i introduced over the weekend,
///    forcing recomb rate to be 2.59e-13
/// - 2010-01-21 JM: Changed ISOTHERMAL_MP stuff to have no reference
///    to (1-gamma).  corrected isothermal temperature calculation.
///    Changed some heap arrays to stack arrays in MP_H::Temperature
///    and prim2local(),local2prim().
/// - 2010-02-01 JM: if parallel, told only proc 0 to write
///    hummer_recomb.txt file
/// - 2010-02-09 JM: Took abs.value of rates in timescales() function
///    (so heating doesn't give negative time!)
/// - 2010-04-10 JM: some changes to MicroPhysics() class -- allowed
///    double-counting of recombination cooling; put in a note to get
///    a better recombination cooling calculation.
/// - 2010-08-18 JM: Added cooling time calculation for MicroPhysics
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.10.11 JM: Moved "MicroPhysics" class to microphysics_v1.cc
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.02.17 JM: Check for negative pressure in Temperature() and
///   Set_Temp() functions before converting variables, so that code doesn't
///   bug out.
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2011.02.28 JM: New comments in an ifdef in RTsinglesource update.
/// - 2011.03.02 JM: It now ignores any tracers not listed in trtype, so we 
///     can have as many passive tracers as we like!
/// - 2011.08.17 JM: timescales() limits for RT_TEST_PROBS added.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.01.26 JM: Got rid of mpiPM. call.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/interpolate.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifndef EXCLUDE_MPV1

#include "global.h"
#include "microphysics/microphysics.h"
using namespace std;



/***************************************************************/
/*************** MP_Hydrogen Functions *************************/
/***************************************************************/


// ##################################################################
// ##################################################################




MP_Hydrogen::MP_Hydrogen(const int nv,
			 const int ntracer,
			 const std::string &trtype,
			 struct which_physics *ephys
			 )
  :
  kB(pconst.kB()),
#ifdef RT_TEST_PROBS
  m_p(2.338e-24), // this is for comparison to Krumholz,stone,gardiner(2007)
#else
  m_p(pconst.m_p()),
#endif
  nv_prim(nv)
{
  cout <<"Welcome to MP_Hydrogen: an easier place to work than MicroPhysics!\n";

#ifdef ISOTHERMAL_MP
  cout <<"***** USING ISOTHERMAL EOS WITH T=T_0(1-X) +T_1(X) !!!!!!!!!!!!!! ******\n";
#endif // ISOTHERMAL_MP

  min_elecf = 1.e-12;  // minimum electron fraction (to seed reactions!).
  // set flags for what processes we are and aren't doing.
  ep.dynamics          = ephys->dynamics;
  ep.cooling           = ephys->cooling;
  ep.chemistry         = ephys->chemistry;
  ep.coll_ionisation   = ephys->coll_ionisation;
  ep.rad_recombination = ephys->rad_recombination;
  ep.phot_ionisation   = ephys->phot_ionisation;
  ep.raytracing        = ephys->raytracing;
  ep.update_erg        = ephys->update_erg;
  if (!ep.chemistry) {ep.coll_ionisation = ep.rad_recombination = ep.phot_ionisation = 0;}
  //  cout <<"\t\tExtra Physics flags set.\n";

  if (!ep.cooling) {
    cout <<"\t\t cooling not needed.\n";
    cool=0;
  }
  else {
      cool=0;
      cool = new CoolingFn(ep.cooling);
      if (!cool) rep.error("CoolingFn init",cool);
  }
  //  cout <<"\t\tCooling Function set up.\n";

  // Set up tracer variables.
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last "<<ntracer<<" variables in state vec.\n";
  int ftr = nv_prim -ntracer; // first tracer variable.
  string s;
  int len = (trtype.length() +5)/6 -1; // first 6 chars are the type, then list of tracers, each 6 chars long.
    cout <<"\t\ttrtype = "<<trtype<<"\n";
    cout <<"\t\tlen="<<len<<", ntr="<<ntracer<<"\n";
  if (len!=ntracer) {
    cout <<"warning: string doesn't match ntracer.  make sure this looks ok: "<<trtype<<"\n";
    //rep.error("string doesn't match ntracer",ntracer-len);
  }

  MP_Hydrogen::lv_nh   = 0;
  MP_Hydrogen::lv_eint = 1;
  MP_Hydrogen::lv_Hp = 2;
  MP_Hydrogen::nvl = 3;
  if (ep.phot_ionisation) {
    //    lvar["dtau"] = 3;
    MP_Hydrogen::lv_dtau = 3;
    nvl += 1;
  }
  Integrator_Base::Set_Nvar(nvl);
  cout <<"nvl = "<<nvl<<" and nv_prim="<<nv_prim<<"\n";

  // Find ionisation fraction in tracer variable list.
  MP_Hydrogen::pv_Hp=-1;
  for (int i=0;i<len;i++) {
    s = trtype.substr(6*(i+1),6); // Get 'i'th tracer variable.
    if (s=="H1+___" || s=="HII__") {
      lv_Hp = 2;
      pv_Hp = ftr+i;
    }
  }
  if (pv_Hp<0)
    rep.error("No H ionisation fraction found in tracer list",trtype);
  
  ion_pot = 13.59844*1.602e-12;

  // error tolerance for integration of microphysics rate equation.
  // 0.01 is not great, and 0.0001 takes a long time, so 0.001 is a good
  // middle ground -- it is very accurate as far as i can tell.
  MP_Hydrogen::errtol = 1.e-3;


#ifdef HUMMER_RECOMB
  //
  // Use Hummer (1994) recombination rates and recombination energy loss rates.  These
  // are given as tables, so we need to fit cubic splines to the data.
  //
  hr_t=hr_alpha=hr_alpha2=hr_beta=hr_beta2=0;
  hr_nspl = 31;

  hr_t      = mem.myalloc(hr_t      , hr_nspl);
  hr_alpha  = mem.myalloc(hr_alpha  , hr_nspl);
  hr_alpha2 = mem.myalloc(hr_alpha2 , hr_nspl);
  hr_beta   = mem.myalloc(hr_beta   , hr_nspl);
  hr_beta2  = mem.myalloc(hr_beta2  , hr_nspl);

  double caseb[31] = {9.283e-11, 8.823e-11, 8.361e-11, 7.898e-11, 7.435e-11,
		      6.973e-11, 6.512e-11, 6.054e-11, 5.599e-11, 5.147e-11,
		      4.700e-11, 4.258e-11, 3.823e-11, 3.397e-11, 2.983e-11,
		      2.584e-11, 2.204e-11, 1.847e-11, 1.520e-11, 1.226e-11,
		      9.696e-12, 7.514e-12, 5.710e-12, 4.257e-12, 3.117e-12,
		      2.244e-12, 1.590e-12, 1.110e-12, 7.642e-13, 5.199e-13,
		      3.498e-13};
  double coolb[31] = {8.287e-11, 7.821e-11, 7.356e-11, 6.892e-11, 6.430e-11,
		      5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, 4.170e-11,
		      3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11,
		      1.736e-11, 1.402e-11, 1.103e-11, 8.442e-12, 6.279e-12,
		      4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13,
		      6.023e-13, 3.738e-13, 2.268e-13, 1.348e-13, 7.859e-14,
		      4.499e-14};
  // this is to do the spline in linear/linear space... interpolation would
  // be better in log space, but then i'd have to take log(T) and exp(Rate)
  // for every call, and they are expensive...
  for (int i=0; i<hr_nspl; i++) {
    hr_t[i] = exp(log(10.0)*(1.0 +0.2*static_cast<double>(i)));
    hr_alpha[i] = caseb[i]/sqrt(hr_t[i]);
    hr_beta[i]  = coolb[i]/sqrt(hr_t[i]);
    hr_alpha2[i]= 0.0;
    hr_beta2[i] = 0.0;
  }
  interpolate.spline(hr_t, hr_alpha, hr_nspl, 1.e99, 1.e99, hr_alpha2);
  interpolate.spline(hr_t, hr_beta,  hr_nspl, 1.e99, 1.e99, hr_beta2 );


#ifdef TESTING
  //
  // can't have all procs fighting over file in parallel, so just
  // don't write if we are running parallel code.
  //
#ifndef PARALLEL
  ofstream outf("hummer_recomb.txt");
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"Hummer Recombination and Cooling Curve Data: Temperature(K) Rate(cm^3/s) Cool(erg.cm^3/s)\n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double t=10.0;
  do {
    outf << t <<"\t"<< rad_recomb_rate(t)*sqrt(t) <<"\t"<< rad_recomb_energy(t)*sqrt(t);
    outf <<"\t"<<rad_recomb_rate(t)<<"\t"<<rad_recomb_energy(t)<<"\t";
    outf <<3.41202e-10*exp(-0.782991*log(t))<<"\t"<< 3.41202e-10*exp(-0.782991*log(t)) *kB*t/(2./3.) << "\n";
    t *=1.03;
  } while (t<1.e7);
  outf.close();
#endif // PARALLEL
#endif // TESTING
#endif // HUMMER_RECOMB

  //  cout <<"\t\tinit done.\n";
  return;  
}


// ##################################################################
// ##################################################################



MP_Hydrogen::~MP_Hydrogen()
{
  if (cool) {delete cool; cool=0;}
#ifdef HUMMER_RECOMB
  hr_t      = mem.myfree(hr_t);
  hr_alpha  = mem.myfree(hr_alpha);
  hr_alpha2 = mem.myfree(hr_alpha2);
  hr_beta   = mem.myfree(hr_beta);
  hr_beta2  = mem.myfree(hr_beta2);
#endif // HUMMER_RECOMB
}


// ##################################################################
// ##################################################################




int MP_Hydrogen::Tr(string t) 
{
  if (t=="H1+" || t=="HII" || t=="Hp") return pv_Hp;
  else return -1;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::Set_Temp(double *p,  ///< primitive vector.
			  const double T, ///< temperature.
			  const double g  ///< eos gamma.
			  )
{
  gamma = g;
  double P[nvl];
  //
  // Check for negative pressure.  If density<0 then we should bug out because
  // there is no way to set a temperature, but if p<0 we can just overwrite it.
  //
  if (p[PG]<=0.0) {
    //cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p[PG] = 1.0e-12;  // It doesn't matter what this is as long as p>0
  } 
  int err = convert_prim2local(p,P,g);
  P[lv_eint] = (1.0+P[lv_Hp])*P[lv_nh]*kB*T/(gamma-1.);
  err += convert_local2prim(P, p, p, g);

  return err;
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::Temperature(const double *pv, ///< primitive vector
				const double g   ///< eos gamma
				)
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    cout <<"MP_Hydrogen::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }

  gamma = g;
  double P[nvl];
	
  convert_prim2local(pv,P,g);
  double t=(gamma-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.
#ifdef ISOTHERMAL_MP
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
  t = P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh];
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
#endif // ISOTHERMAL_MP

  return t;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::Init_ionfractions(double *p_prim,  ///< Primitive vector to be updated.
				   const double gam, ///< eos gamma.
				   const double temp ///< optional gas temperature to end up at. (negative means use pressure)
				   )
{
  gamma = gam;
  double T=temp;
  double p_local[nvl];

  int err = convert_prim2local(p_prim, p_local, gamma);
  if (T>0) p_local[lv_eint] = (1.0+p_local[lv_Hp])*p_local[lv_nh]*kB*T/(gamma-1.);
  else T = (gamma-1.)*p_local[lv_eint]/kB/(1.0+p_local[lv_Hp])/p_local[lv_nh]; // Temperature.
  //  cout <<"Temp = "<<T<<"\n";

  // Now set ionisation fraction to 0.5 to start with.
  p_local[lv_Hp] = 0.5;
  err += convert_local2prim(p_local, p_prim, p_prim, gamma);

  return err;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::convert_prim2local(const double *p_in,
				    double *p_local,
				    const double gam
				    )
{
  p_local[lv_nh]   = p_in[RO]/m_p;
  p_local[lv_eint] = p_in[PG]/(gam-1.);
  if (p_in[PG]<=0.) {
#ifdef TESTING
    commandline.console("Bummer... negative pressure input to MP! >");
#endif // testing
    cout <<"neg.pres. input to MP: e="<<p_local[lv_eint]<<"\n";
    rep.error("Negative pressure input to RT solver!",p_in[PG]);  
  }
#ifdef TESTING
  if (p_local[lv_eint]>1.e-5) {
    //cout <<"cell with crazy temperature:"; CI.print_cell(dp.c);
    //rep.error("goodbye",0);
  }
#endif // testing
  p_local[lv_Hp]   = p_in[pv_Hp];
  if (p_local[lv_Hp]<0.0 || (p_local[lv_Hp]-1.0>SMALLVALUE)) {
#ifdef TESTING
    //commandline.console("Bummer... bad ion frac. input to MP! >");
#endif // testing
    //    rep.warning("bad ion frac. input to MP_Hydrogen()",0.5,p_local[lv_Hp]);
    //rep.error("bad ion frac. input to MP_Hydrogen()",p_local[lv_Hp]);
    if (p_local[lv_Hp]<0.0) p_local[lv_Hp]=min_elecf;
    else if (p_local[lv_Hp]>=1.0) p_local[lv_Hp]=1.0-SMALLVALUE;
    else cout <<"bad ion frac, but not bad enough!!! : "<<p_local[lv_Hp]<<"\n";
  }
  if (ep.phot_ionisation) {
    p_local[lv_dtau] = 0.0;
    tau_cell = 0.0;
  }
  //cout <<"nh="<<p_local[lv_nh]<<"\n";

  //
  // Check for NAN/INF
  //
  for (int v=0;v<nvl;v++) {
    if (isnan(p_local[v]) || isinf(p_local[v])) rep.error("INF/NAN input to microphysics",p_local[v]);
  }

#ifdef ISOTHERMAL_MP
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
  p_local[lv_eint] = (1.0+p_local[lv_Hp])*p_local[lv_nh]*kB*(100.0+9900.0*p_local[lv_Hp]);
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
#endif // ISOTHERMAL_MP

  return 0;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::convert_local2prim(const double *p_local,
				     const double *p_in,
				     double *p_out,
				     const double gam
				     )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];
  p_out[PG] = p_local[lv_eint]*(gam-1.0);

#ifdef ISOTHERMAL_MP
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
  double t1 = p_local[lv_eint]/kB/(1.0+p_local[lv_Hp])/p_local[lv_nh];
  double t2 = 100.0+9900.0*p_local[lv_Hp];
  if (fabs((t1-t2)/(t1+t2))>0.1 || t1<100.0 || t2<100.0) {
    cout <<"integration obtained T="<<t1;
    cout <<"\t Temperature should be T="<<t2<<"; resetting to what it should be.\n";
  }
  p_out[PG] = p_local[lv_nh]*(1.0+p_local[lv_Hp])*kB*(100.0+9900.0*p_local[lv_Hp]);
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
#endif // ISOTHERMAL_MP

  if (p_out[PG]<0.) {
#ifdef TESTING
    commandline.console("Bummer... negative pressure! >");
#endif // testing
    cout <<"neg.pres. e="<<p_local[lv_eint]<<"\n";
    rep.error("Negative pressure output from RT solver!",p_out[PG]);
  }
  p_out[pv_Hp] = max(min_elecf, p_local[lv_Hp]);
  p_out[pv_Hp] = min(1.0, p_out[pv_Hp]);

  if (ep.phot_ionisation) {
    tau_cell = p_local[lv_dtau]; // this should be int(exp(-tau),dt)
  }

  //
  // Check for nans/infs (easier to do this from calling function who has cell pointer).
  //
  //for (int v=0;v<nv_prim;v++) {
  //  if (isnan(p_out[v]) || isinf(p_out[v])) rep.error("INF/NAN",p_out[v]);
  //}

  return 0;
}



// ##################################################################
// ##################################################################



int MP_Hydrogen::TimeUpdateMP(const double *p_in,
			       double *p_out,
			       const double dt,
			       const double g,
			       const int sw_int,
			       double *ttt
			       )
{
  int err = 0;
  gamma = g;
  //cout <<"gamma: "<<gamma<<"\n";

  double P[nvl];  // local state vector for microphysics = [n_h,Eint,x_Hp,d_tau]

  // put data from p_in into local vector.  
  err += convert_prim2local(p_in,P,gamma);

  //rep.printVec("p_prim",p_in,nv_prim);
  //rep.printVec("p_old",P,nvl);
  //  double errtol = 1.0e-5;
  double tout =0.0;
  if (sw_int==0) {
    //errtol = 1.0e-8;
    err = Integrator_Base::Int_Adaptive_RKCK(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int==1) {
    //errtol = 1.0e-3;
    err = Integrator_Base::Int_DumbAdaptive_Euler(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int==2) {
    // DANGEROUS!!!
    cout <<"Requesting single step RK4 integration for MP_Hydrogen!!!?\n";
    err = Integrator_Base::Step_RK4(nvl, P, 0.0, dt, P);
  }
  else rep.error("this integration method not known.",sw_int);
  if (err) rep.error("integration failed.",err);
  //rep.printVec("p_new",P,nvl);

#define CHECK_IFRAC
#ifdef CHECK_IFRAC
  if (P[lv_Hp]>2.0) {
    cout <<"H+ has i-frac="<<P[lv_Hp]<<"  something went wrong, attempting to split the integral.\n";
    rep.printVec("Perror",P,nvl);
    int nstep=10;
    err += convert_prim2local(p_in,P,gamma); // reset P to the initial values.
    rep.printVec("Pinit ",P,nvl);
    for (int el=0;el<nstep;el++) {
      err += Integrator_Base::Int_Adaptive_RKCK(nvl, P, 0.0, dt/static_cast<double>(nstep), errtol, P, &tout);
      //      rep.printVec("P",P,nvl);
    }
    rep.printVec("Pend  ",P,nvl);
    if (err!=0 || P[lv_Hp]>2.0) {
      cout <<"Tried to split integral, but still got errors: err="<<err<<", i-frac="<<P[lv_Hp]<<"\n";
      rep.error("integration failed!",err);
    }
  }
#endif // CHECK_IFRAC

  if (P[lv_Hp]>1.00001) {
    cout <<"H+ has i-frac="<<P[lv_Hp]<<"  ...setting to 1.\n";
    P[lv_Hp] = 1.0;
  }
  // put updated state vector into p_out.
  *ttt = (gamma-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.
  //if (*ttt >1.e3) 
  //  cout <<"\tT="<<*ttt<<" and n_atoms= "<<P[lv_nh]<<"\n";
  err += convert_local2prim(P, p_in, p_out, gamma);

  return err;
}


// ##################################################################
// ##################################################################




int MP_Hydrogen::TimeUpdate_RTsinglesrc(
      const double *p_in,   ///< Primitive Vector to be updated.
      double *p_out,        ///< Destination Vector for updated values.
      const double dt,      ///< Time Step to advance by.
      const double g,       ///< EOS gamma.
      const int sw_int,     ///< Switch for what type of integration to use.
                            ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      const double phot_in, ///< flux in per unit length along ray (F/ds or L/dV)
      const double ds,      ///< path length ds through cell.
      const double tau2cell, ///< Optical depth to entry point of ray into cell.
      double *deltau        ///< return optical depth through cell in this variable.
      )
{
#ifdef RT_TESTING
  //
  // display inputs:
  //
  cout <<"gamma="<<g<<", sw_int="<<sw_int<<", phot_in="<<phot_in;
  cout <<", ds="<<ds<<", tau2cell="<<tau2cell<<"\n";
#endif // RT_TESTING

  if (!ep.phot_ionisation) rep.error("RT requested, but phot_ionisation not set!",ep.phot_ionisation);
  MP_Hydrogen::tau_cell=0.0;
  MP_Hydrogen::photons_in = phot_in; // units are photons/cm^3/s (/Hz if using frequency info).

  //
  // with new interface, the flux is passed in unattenuated
  //
  photons_in *= exp(-tau2cell);
  //
  // Now we calculate the FUV flux and FUV optical depth.
  // This is based on equation A3 in Henney, Arthur, et al., 2009.
  //
  double f_UV = 0.5; // number of UV photons (6-13eV) per ionising photon.
  FUV_unattenuated_flux = phot_in*f_UV*ds/1.2e7; // Flux in units of "Habing Flux" according to Henney et al. 2009.
  FUV_extinction        = 1.086*5.0e-22*tau2cell/phot_xsection(1.0);  // Extinction from Henney et al. 2009.
  FUV_attenuated_flux   = FUV_unattenuated_flux*exp(-1.9*FUV_extinction);

  MP_Hydrogen::path_length = ds;
  double temp=0.0;
  srand(124356);
  if (!isfinite(photons_in)) rep.error("photons_in is NAN!",photons_in);

  //
  // SETTING FLUX TO ZERO IF IT IS *VERY* WEAK, TO SAVE COMPUTATION.
  //
  // This needs to be treated carefully -- if I have large volumes
  // photons_in can be a very small number...  We can interpret it as
  // flux per unit length, so the relevant quantity is then
  // (photons_in*ds), which should be very small.  How small is very?
  // If we have 1 photon per cm2 per sec, then the PI rate (at the
  // front edge of the grid is ~6.3e-18 b/c that's the cross-section.
  // This gives a PI time of 5e9 years, so at this stage we have
  // reached seriously cautious limits.
  //
  if (photons_in*path_length<1.0) {
    photons_in=0.0; //ep.phot_ionisation=0;
  }
  //else cout <<"\tphot_in: "<<photons_in<<"  path: "<<path_length<<" \n";

  int err = 0;
  MP_Hydrogen::gamma = g;
  //  cout <<"gamma: "<<gamma<<"\n";

  //
  // local state vectors for microphysics = [n_h,Eint,x_Hp,d_tau]
  //
  double *P=0, *P2=0, *R=0;
  P  = mem.myalloc(P, nvl);
  P2 = mem.myalloc(P2,nvl);
  R  = mem.myalloc(R, nvl);

  //
  // put data from p_in into local vector.  
  //
  err += convert_prim2local(p_in,P,gamma);
  for (int i=0;i<nvl;i++) P2[i]=R[i]=0.0;
  //cout <<"nh: "<<P[lv_nh]<<"\n";  cout <<"e_i:"<<P[lv_eint]<<"\n";  cout <<"xh: "<<P[lv_Hp]<<"\n";
  //cout <<"tau:"<<P[lv_dtau]<<"\n";
  //cout <<"nh: "<<P[0]<<"\n";  cout <<"e_i:"<<P[1]<<" and lv_Hp="<<lv_Hp<<"\n";
  //cout <<"xh: "<<P[2]<<"\n";  cout <<"tau:"<<P[3]<<"\n";


#ifdef COUNT_ENERGETICS
  have_counted_ergs=false;
  double *rrr=0;
  rrr = mem.myalloc(rrr,nvl);

  for (int i=0;i<nvl;i++) rrr[i]=0.0;
  dPdt(nvl,P,rrr);
  have_counted_ergs=true;
#endif

  // First calculate A_phot/(C_rec*n_e) to see if the reaction is strongly
  // driven to high ionisation or not.
  double T = (gamma-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.
  double irate = P[lv_nh]*(1.0-P[lv_Hp])*phot_xsection(T)*path_length; // Tau through cell

  if (irate<0.01)
    irate = photons_in*phot_xsection(T)*path_length;
  else
    irate = photons_in*(1.0 -exp(-irate))/P[lv_nh]/(1.0-P[lv_Hp]);
  // now we have A_phot, so divide by C_rec*n_e
  irate /= rad_recomb_rate(T)*P[lv_Hp]*P[lv_nh];

  //
  // irate is a measure of how strong the photo-ionisation is.  If it is weak,
  // we should be able to do a simple RK5 adaptive integration, and if not we
  // will have to try the explicit/implicit method.  I set the dividing line at
  // irate=20.
  //
  double t_now=0.0, t_final=dt, hh=dt;
  
  MP_Hydrogen::i_crit=0.95;
  if (irate<20.0) {
    // ionisation rate is not super strong, so do regular explicit integration.
    double t_out=0.0;
    if (sw_int==0) {
      err += Integrator_Base::Int_Adaptive_RKCK(nvl, P, 0.0, hh, errtol, P2, &t_out);
    }
    else rep.error("Only allowed to Adaptive_RKCK Integration here!",sw_int);
    //if (err) rep.error("integration failed.",err);
    //if (irate>30) {cout <<"irate = "<<irate<<"  ";rep.printVec("p_new",P2,nvl);}

    //
    // Do some careful error checking here to try to catch and recover from
    // failed integration attempts:
    //

    //
    // First if there was an error, just try to move on and do the full 
    // integration with the MP_Hydrogen integrator:
    //
    if (err) {
      cout <<"base integrator failed, so continuing on and trying MPH integrator.\n";
    }
    // if it overshot, reject step and go on to more detailed method below.
    if (P2[lv_Hp]>1.0) {
      // Don't need to do anything; keep P[] and t_now as they are.
      //cout <<"H+ has i-frac="<<P2[lv_Hp]<<" and t_out/hh="<<t_out/hh;
      //cout <<"  ...rejecting step. integrator gave code: err="<<err<<"\n";
      //P2[lv_Hp] = 1.0;
      //cout <<"irate = "<<irate<<"\n";
#ifdef TESTING
      commandline.console("Bummer>");
#endif
    }
    else if (!pconst.equalD(t_out,hh)) {
      t_now += t_out;
      for (int i=0;i<nvl;i++) P[i]=P2[i];
      //cout <<"Integration didn't go for full timestep.  Too much work!\n";
      //cout <<"Will try more careful method to get the rest of the way.\n";
      //cout <<"integrator gave code: err="<<err<<"; t="<<t_now<<" and dt="<<t_final<<"\n";
      //cout <<"irate = "<<irate<<"\n";
    }
    else {  // step worked.
      if (err) {
	cout <<"integrator gave code: err="<<err<<"; t="<<t_out<<" and dt="<<t_final<<"\n";
	rep.error("Integration failed, but I didn't catch why! Fix me.",err);
      }
      for (int i=0;i<nvl;i++) P[i]=P2[i];
      t_now += hh;
    }
    if (P[lv_eint]<0.0 || isnan(P[lv_eint])) rep.error("explicit integral gives negative energy!",P[lv_eint]);
  } // end of weak ionisation step.
  
  //
  // Now we know we have strong ionising radiation (or previous step failed),
  // so do the careful integration method.
  //
  err=0;
  while (P[lv_Hp]<i_crit && !pconst.equalD(t_now,t_final)) {
    //cout <<"Doing Explicit Step:\n";
    dPdt(nvl,P,R);
    temp = max(0.1,P[lv_Hp]);
    if (temp>0.5) temp=min(0.05,1.0-P[lv_Hp]);
    //    hh = min(dt,P[lv_Hp]/R[lv_Hp]);
    hh = min(dt,temp/R[lv_Hp]);
    //cout <<"strongly driven: irate="<<irate<<" and x/xdot="<<hh<<" with dt="<<dt<<"\n";
    //rep.printVec("P",P,nvl);
    //rep.printVec("R",R,nvl);
    
    hh = min(hh,t_final-t_now);
    if (hh<0.0) {
      //cout <<"hh="<<hh<<", tf-tnow="<<t_final-t_now;
      //cout <<", resetting hh to dt="<<t_final-t_now<<"\n";
      hh=t_final-t_now;
    }
    //cout <<"\tDoing Explicit Step: x="<<P[lv_Hp]<<", t_now="<<t_now<<" and hh="<<hh<<" and t_final="<<t_final<<"\n";
    double t_out=0.0;
    if (sw_int==0) {
      err += MP_Hydrogen::Int_Adaptive_RKCK(nvl, P, 0.0, hh, errtol, P, &t_out);
    }
    else rep.error("Only allowed to Adaptive_RKCK Integration here!",sw_int);
    if (err) rep.error("integration failed.",err);
    //rep.printVec("p_new",P,nvl);
    if (P[lv_Hp]>1.00001) {
      cout <<"H+ has i-frac="<<P[lv_Hp]<<"  ...setting to 1.\n";
      //      P[lv_Hp] = 1.0;
#ifdef TESTING
      commandline.console("Bummer>");
#endif
      rep.error("FAILURE of method on every level! FIX ME!",P[lv_Hp]);
    }
    if (!pconst.equalD(t_out,hh)) {
      //cout <<"integration overshot, so cut short! req: "<<hh<<" did: "<<t_out<<"\n";
      //rep.printVec("p_new",P,nvl);
      //commandline.console("interr Bummer>");
      //rep.error("integration didn't go for specified time!",(t_out-hh)/(t_out+hh));
      hh = t_out;
    }
    t_now += hh;
    if (P[lv_eint]<0.0 || isnan(P[lv_eint])) rep.error("explicit high-Irate integral gives negative energy!",P[lv_eint]);
  }

  //
  // now if t_now<t_0+dt, set delt=remaining time, set oldval[i]=10^6
  // do {
  //   do implicit timestep forwards by delt
  //      for energy, x_i, dtau
  //   err = sum(oldval[i]-newval[i]/(old+new+1.e-100)
  //   half the timestep.
  // } while (err>0.01);
  //
  if (t_now<t_final && !pconst.equalD(t_now,t_final)) {
    //rep.printVec("begin implicit P",P,nvl);
    hh = t_final-t_now;
    err += implicit_step(nvl,hh,P,P,errtol);
    t_now+=hh;
    if (P[lv_eint]<0.0 || isnan(P[lv_eint])) {
      rep.printVec("P",P,nvl);
      rep.error("Implicit integral gives negative energy!",P[lv_eint]);
    }
  }
  
  if (!pconst.equalD(t_now,t_final)) {
    cout <<"irate="<<irate<<" and integration didn't go for right timestep...\n";
    cout <<"dt="<<dt<<"\tt_now="<<t_now<<"\tt_final="<<t_final<<"\n";
    rep.error("integration didn't go for right timestep!",t_now/t_final);
  }

  // put updated state vector into p_out.
  T = (gamma-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.
  //if (T >1.e3) 
  //  cout <<"\tT="<<T<<" and n_atoms= "<<P[lv_nh]<<"\n";
  //rep.printVec("P",P,nvl);
  err += convert_local2prim(P, p_in, p_out, gamma);

  //
  // Convert integrated exp(-tau) into a time-averaged optical depth to return to the code.
  // Have to check that it is finite and positive, since log is sensitive to roundoff errors.
  //
  *deltau = -log(tau_cell/dt);
  if (isinf(*deltau) || isnan(*deltau)) {
    // this can happen if the ray is very attenuated, and we are taking -log(0.0)
    //cout <<" <exp(-tau)> = "<<tau_cell<<" dt="<<*deltau<<"\n";
    *deltau = 1000.0; // This should be big enough -- exp(-1000)=10^(-440)!!!
  }
  if (*deltau<0.0) {
    if (*deltau<-1.e-14)
      cout <<"MP_Hydrogen: <tau> significantly negative: "<<*deltau<<" setting to zero. tau_cell/dt-1="<<tau_cell/dt-1.0<<"\n";
    *deltau=0.0;
  }

  P  = mem.myfree(P);
  P2 = mem.myfree(P2);
  R  = mem.myfree(R);
#ifdef COUNT_ENERGETICS
  rrr= mem.myfree(rrr);
#endif
  return err;
}




// ##################################################################
// ##################################################################



int MP_Hydrogen::implicit_step(
      const int nv,    ///< number of variables we are expecting.
      const double dt, ///< timestep to use.
      const double *P, ///< Current state vector.
      double *p_out,   ///< Final State Vector.
      const double etol ///< error tolerance
      )
{
  //
  // Do some paranoid checking first -- right number of variables and no NAN/INF in 
  // input vector!!!
  //
  //cout <<"\t\t\timplicit step!\n";
  //rep.printVec("\t\tP implicit:",P,nvl);

  if (nv!=nvl) rep.error("variables wrong!",nv-nvl);
  for (int i=0;i<nvl;i++) {
    if (isnan(P[i]) || isinf(P[i])) {
      rep.printVec("\t\tP implicit:",P,nvl);
      rep.error("\t\tnan passed to implicit step.",i);
    }
  }

  //
  // Allocate memory for temp arrays.
  //
// #ifdef USE_MM
//   double *p_old=0, *p_now=0;
// #ifdef TESTING
//   p_old = mem.myalloc(p_old, nvl, "MP_H:Implicit_step: p_old");
//   p_now = mem.myalloc(p_now, nvl, "MP_H:Implicit_step: p_now");
// #else
//   p_old = mem.myalloc(p_old, nvl);
//   p_now = mem.myalloc(p_now, nvl);
// #endif //TESTING
// #else
//   double p_old[nvl],p_now[nvl];
// #endif
  double p_old[nvl],p_now[nvl];
  for (int i=0;i<nvl;i++) {
    p_old[i]=p_now[i]=1.e100;
  }
  double maxerr=1.e100;
  int ct=0, nsub=1, max_ct=30;
  double hh=dt/nsub;

  double x_inf,t_ion, A=0.0,B=0.0,C=0.0,LL=0,T=0.0, tau, tau_int, e_int;
  double e_phot, e_ci, e_rr;
  //  static long int count=0;

  //
  // Run a do-loop, where we double the number of sub-steps in each iteration, starting with
  // a single step over the whole timestep.
  // The loop runs at least twice so we can estimate the relative error.
  //
  do {
    //
    // prepare state vectors:
    //   p_old is the final result of the previous loop.
    //   p_now is what we start the current loop with, and increment.
    //   p_old could have NAN/INF, so set it to a huge value in that case.
    //
    for (int i=0;i<nvl;i++) {
      p_old[i]=p_now[i];
      p_now[i]=P[i];
      if (isnan(p_old[i]) || isinf(p_old[i])) p_old[i] = -HUGEVALUE*ct;
    }

    //
    // If temperature is <1K then something has gone wrong, so set it to 1K
    // 
//    T = (gamma-1.)*p_now[lv_eint]/kB/(1.0+p_now[lv_Hp])/p_now[lv_nh]; // Temperature.
//    if (T<1.0) {
//      cout <<"Implicit_step: fixing very low Temperature="<<T<<" and setting to 1K.\n";
//      p_now[lv_eint]=(1.0+p_now[lv_Hp])*p_now[lv_nh]*kB*1.0/(gamma-1.);
//    }
    
    
    //
    // Now loop over each substep, doing the analytic integration I
    // have derived in the notes. (Tau is integrated numerically, but
    // x and e are integrated analytically based on the rates).
    //
    //cout <<"count = "<<ct<<" and nsub="<<nsub<<"\n";
    //rep.printVec("\tp_start",p_now,nvl);

    for (int i=0;i<nsub;i++) {
      //edot=0.0;
      //
      // Set the temperature based on p_now[].  If the temperature is not valid,
      // there is no point continuing the integration, so set p_now[] to have crazy
      // values, and go on to the next iteration.
      T = (gamma-1.)*p_now[lv_eint]/kB/(1.0+p_now[lv_Hp])/p_now[lv_nh]; // Temperature.
#ifdef TESTING
      if (dp.c->id==17042) {
	cout <<"\t\t\ttemperature="<<T<<" gamma="<<gamma<<"\n";
	rep.printVec("state:",p_now,4);
	rep.printVec("  old:",p_old,4);
      }
#endif
	
      if (T<0.0 || isnan(T) || isinf(T)) {
#ifdef TESTING
	cout <<"\t\t\ttemperature="<<T<<" gamma="<<gamma<<"\n";
	rep.printVec("state:",p_now,4);
	//CI.print_cell(dp.c);
#endif
	for (int v=0;v<nvl;v++) {
	  p_now[v] = -HUGEVALUE*(i+1.1);
	}
      }
      else {
	//
	// I have functions for all these energy gains and losses so that the 
	// explicit and implicit integration methods boths use the same bits
	// of code!
	//
	e_phot = phot_ion_energy(T); // positive value.
	e_ci   = coll_ion_energy(T); // positive value.
#ifdef HUMMER_RECOMB
	e_rr   = rad_recomb_energy(T)*p_now[lv_Hp]*p_now[lv_nh]; // erg/s = lambda*n_e (positive value)
#endif // HUMMER_RECOMB
#ifndef HUMMER_RECOMB
	e_rr   = rad_recomb_energy(T); // erg (positive value)
#endif // not HUMMER_RECOMB
	
	if (ep.phot_ionisation) {
	  tau = p_now[lv_nh]*(1.0-p_now[lv_Hp])*phot_xsection(T)*path_length;
	  if (tau<0.01)
	    A = photons_in*phot_xsection(T)*path_length;
	  else
	    A = photons_in*(1.0 -exp(-tau))/p_now[lv_nh]/(1.0-p_now[lv_Hp]);
	}
	else rep.error("Why implicit if no phot_ionisation?",10);
	
	if (ep.coll_ionisation) {
	  B = coll_ion_rate(T)*p_now[lv_Hp]*p_now[lv_nh]; // This is A_ci*n_e
	}
	else B=0.0;
	
	if (ep.rad_recombination) {
	  C = rad_recomb_rate(T)*p_now[lv_Hp]*p_now[lv_nh]; // This is alpha_rr*n_e
	}
	else C=0.0;
	
	if (ep.cooling) {
	  LL = cool->CoolingRate(T,p_now[lv_Hp],p_now[lv_nh],FUV_unattenuated_flux,FUV_extinction)/p_now[lv_Hp]/p_now[lv_nh]; // This is Lambda*n_e
#ifdef NO_DOUBLECOUNTING
	e_rr = 0.0; // Don't want to double count recombination cooling, if it
	            // is already in the cooling function!
#endif // NO_DOUBLECOUNTING
	}
	else LL=0.0;
	
	// cout <<"A="<<A<<" B="<<B<<" C="<<C<<"\n";
	t_ion = 1.0/(A+B+C);
	x_inf = (A+B)*t_ion;
#ifdef MP_DEBUG
	//if (isnan(t_ion) || isnan(x_inf)) rep.error("xinf,tion failed",t_ion);
#endif //MP_DEBUG
	
	// Add to int(edot,dt)
	// First get int(x,dt), and store in temporary variable e_int.
	e_int = hh*x_inf - t_ion*(x_inf-p_now[lv_Hp])*(1.0-exp(-hh/t_ion));

	// now int(edot,dt) = (A_phot*e_phot -B_ci*E_ci)hh -(A_phot*E_phot -B_ci*E_ci +C_rr*E_rr +Lambda*n_e)int(x,dt)
#ifdef HUMMER_RECOMB
	// here e_rr is a cooling rate, just like LL
	e_int = hh*(A*e_phot -B*e_ci) -e_int*(A*e_phot -B*e_ci +e_rr   +LL);
#endif // HUMMER_RECOMB
#ifndef HUMMER_RECOMB
	// here e_rr is just an energy per recombination, which is not quite right since
	// it is supposed to be inside the integration over the boltzmann distribution in
	// the rate calculation.  hummer_recomb is better!
	e_int = hh*(A*e_phot -B*e_ci) -e_int*(A*e_phot -B*e_ci +C*e_rr +LL);
#endif // not HUMMER_RECOMB
      
#ifdef TESTING
	//cout <<"ct="<<ct<<" i="<<i<<" out of "<<nsub<<"; hh="<<hh<<" T="<<T<<", pnow=";rep.printVec("pnow",p_now,nvl);
	//cout <<"\t\tLL="<<LL<<" e_rr="<<e_rr<<" e_int="<<e_int<<" x_inf="<<x_inf<<" photons_in="<<photons_in<<" A="<<A<<" B="<<B<<" C="<<C<<"\n";
	//cout <<"*** 1st="<<hh*(A*e_phot -B*e_ci)<<" 2nd=("<<hh*x_inf - t_ion*(x_inf-p_now[lv_Hp])*(1.0-exp(-hh/t_ion))<<")*(";
	//cout <<(A*e_phot -B*e_ci +e_rr   +LL)<<")\n";
#endif
      
	if (ep.update_erg)
	  p_now[lv_eint] += e_int*p_now[lv_nh];

	if (p_now[lv_eint]<0.) {
	  //cout <<"e_old= "<<p_now[lv_eint]-e_int*p_now[lv_nh]<<"\tde = "<<e_int*p_now[lv_nh]<<"\n";
	  p_now[lv_eint] -= HUGEVALUE*nsub;  // make sure integration fails if we get negative energy!!!
	}
      
#ifdef ISOTHERMAL_MP
        // **************************************************
	// ************** BIG ISOTHERMAL HACK!!!! ***********
	// **************************************************
	p_now[lv_eint] = p_now[lv_nh]*(1.0+p_now[lv_Hp])*pconst.kB()*(100.0+9900.0*p_now[lv_Hp]);
        // **************************************************
	// ************** BIG ISOTHERMAL HACK!!!! ***********
	// **************************************************
#endif // ISOTHERMAL_MP


	//p_now[lv_nh] stays constant.

#define SIMPLE_TAU
#ifdef SIMPLE_TAU      
	// simple tau_int
	// This is 1st order integration (not used anymore)
	//tau = path_length*phot_xsection(T)*p_now[lv_nh]*(1.-p_now[lv_Hp]);
	//if (tau <0.01) tau_int = hh*(1.0-tau);
	//else tau_int = hh*exp(-tau);
	//p_now[lv_Hp] = x_inf+ (p_now[lv_Hp]-x_inf)*exp(-hh/t_ion);
	//p_now[lv_dtau] += tau_int;
	
	// This is 2nd order integration int=0.5h(x(0)+x(1))
	// first get lower limit
	tau = path_length*phot_xsection(T)*p_now[lv_nh]*(1.-p_now[lv_Hp]);
	if (tau <0.01) tau_int = (1.0-tau);
	else tau_int = exp(-tau);
	// update ion fraction
	p_now[lv_Hp] = x_inf+ (p_now[lv_Hp]-x_inf)*exp(-hh/t_ion);
	// get upper limit and make integral from 0.5*hh*(lower+upper)
	tau = path_length*phot_xsection(T)*p_now[lv_nh]*(1.-p_now[lv_Hp]);
	if (tau <0.01) tau_int = 0.5*hh*((1.0-tau) +tau_int);
	else tau_int = 0.5*hh*(exp(-tau) +tau_int);
	// add integral to sum.
	p_now[lv_dtau] += tau_int;
#endif
#ifndef SIMPLE_TAU
	// new tau int... more accurate, but takes much longer to run with all
	// the exp() calls, so it is better to use the simple guesstimate above.
	double xx[5],tau0,tau1; xx[0]=0.0;
	if (t_ion<hh/10.0) {
	  xx[1] = t_ion/100.;
	  xx[2] = t_ion/10.;
	  xx[3] = t_ion;
	  xx[4] = hh;
	}
	else {
	  xx[1] = hh/1000.;
	  xx[2] = hh/100.;
	  xx[3] = hh/10.;
	  xx[4] = hh;
	}
	tau_int=0.0;
	tau0 = exp(-path_length*phot_xsection(T)*p_now[lv_nh]*(1.-p_now[lv_Hp]));
	for (int ii=1;ii<5;ii++) {
	  p_now[lv_Hp] = x_inf+ (p_now[lv_Hp]-x_inf)*exp(-xx[ii]/t_ion);
	  tau1 = exp(-path_length*phot_xsection(T)*p_now[lv_nh]*(1.-p_now[lv_Hp]));
	  //tau_int += 0.5*(xx[ii]-xx[ii-1])*(tau1+tau0);
	  tau_int += (xx[ii]-xx[ii-1])*(tau1);
	  tau0=tau1;
	}
	p_now[lv_dtau] += tau_int;
#endif
      }  // else (we have real temperature, so do the integrations
      
    }    // end of single substep

    //rep.printVec("\tp_now",p_now,nvl);
    //rep.printVec("\tp_old",p_old,nvl);


    //
    // Now we have done all the substeps in this iteration, so check if the old and new
    // values are similar, and if they satisfy the relative error tolerance.
    //
    maxerr = TINYVALUE;
    for (int v=0;v<nvl;v++) {
      maxerr = max(maxerr,fabs(p_now[v]-p_old[v])/(fabs(p_now[v])+fabs(p_old[v])+TINYVALUE));
      if (isnan(p_now[v]) || isinf(p_now[v])) maxerr = 1.e100;
      if (isnan(maxerr) || isinf(maxerr)) maxerr = 2.e100;
    }
    if (p_now[lv_eint]<=0.0) maxerr = 1.e100;

    //
    // Increment counter, double no. of substeps, halve the timestep in preparation
    // for the next iteration (if needed).
    //
    ct++;
    nsub *= 2;
    hh   /= 2.0;
  } while (maxerr>etol && ct<max_ct);

  //  count += ct;   

  if (ct>=max_ct) {
    cout <<"\timplicit_step() failed to converge to tol="<<etol<<" in "<<ct<<" steps.\n";
    cout <<"\timplicit_step() 1-x="<<1.0-P[lv_Hp]<<" ";rep.printVec("P",P,nvl);
    cout <<"\timplicit_step() 1-x="<<1.0-p_now[lv_Hp]<<" ";rep.printVec("Pnow",p_now,nvl);
    cout <<"\timplicit_step()";rep.printVec("Pold",p_old,nvl); cout <<"\n\n";
    rep.error("GG",3);
  }
  else ct=0;
  for (int v=0;v<nvl;v++) p_out[v]=p_now[v];

// #ifdef USE_MM
// #ifdef TESTING
//   p_old = mem.myfree(p_old, "MP_H:Implicit_step: p_old");
//   p_now = mem.myfree(p_now, "MP_H:Implicit_step: p_now");
// #else
//   p_old = mem.myfree(p_old);
//   p_now = mem.myfree(p_now);
// #endif //TESTING
// #endif

  return ct;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::Int_Adaptive_RKCK(const int nv,   ///< number of elements in P array.
				   const double *p0, ///< initial state vector.
				   const double t0,  ///< initial time
				   const double dt,  ///< timestep to advance by.
				   const double etol, ///< error tolerance per step.
				   double *pf,  ///< final state vector
				   double *tf   ///< pointer to final time.
				   )
{
  //cout <<"\t\t\tMP_Hydrogen explicit step!\n";
  if (nvl!=nv) {
    cerr <<"Integrator_Base() nvar not equal to state vector length.\n";
    Set_Nvar(nv);
  }
  if (etol<0) rep.error("Int_Adaptive_RKCK() ErrTol is negative!",etol);
  if (etol<MACHINEACCURACY)
    rep.error("errtol beyond machine accuracy.  use more lenient value!\n",etol);
  if (etol>2)
    rep.error("errtol is too large: relative accuracy required is >1!\n",etol);
  double t = t0;

  double p1[nvl];
  double p2[nvl];

  for (int v=0;v<nvl;v++) p1[v] = p0[v];
  
  *tf = t0+dt;
  double h=dt, hdid=0.0, hnext=0.0;
  int err=0;
  int ct=0, ctmax=200;
  
  //
  // Modified this step in two ways:
  //  1) if stepper went too far, reject the step, halve dt and go again.
  //  2) if H ionised fraction >i_crit, quit loop.
  // Note that the stepper now checks for NAN/INF in the temp and output vectors, so should
  // always return a valid vector if it returns without errors.
  //
  do {
    //rep.printVec("adaptive p1",p1,nv);
    err += Stepper_RKCK(nv, p1, t, h, etol, p2, &hdid, &hnext);
    //cout <<"XX: h="<<h<<" hdid="<<hdid<<" hnext="<<hnext<<" told="<<t<<" tnew="<<t+hdid<<" to-go="<<*tf-t-hdid<<"\n";
    if (p2[lv_Hp]>=0.999) {
      // reject step in this case.
      cout <<"rejecting step...\n";
      h = hdid/2.0;
    }
    else {
      // accept the step.
#ifdef TESTING
      if (dp.c->id==2773) {
	cout <<"\t\t*** energy before = "<<p1[lv_eint]<<" and after step = "<<p2[lv_eint]<<"\n";
      }
#endif
      //cout <<"before="<<t;
      t += hdid;
      //cout <<" and after="<<t<<"\n";
      //t = t+hdid;
      h = min(hnext, *tf-t);
      ct++;
      for (int v=0;v<nvl;v++) p1[v] = p2[v];
      //cout <<"ADAPTIVE INT:\t t_old="<<t-hdid<<" t_new="<<t<<" iter="<<ct<<"\n";
      //cout <<"ADAPTIVE INT:\t "; rep.printVec("pnew",p1,nvl);
      //cout <<"ADAPTIVE INT:\t\thdid="<<hdid<<" tf-t="<<*tf-t<<"\n";
      //cout <<"s: "<<sizeof(h)<<"\ttf-t-hdid="<<*tf-t-hdid<<"\th="<<h<<"\n";
      //cout <<"\t\terr="<<err<<"\tct="<<ct<<"\tp1[lv_Hp]-i_crit="<<p1[lv_Hp]-i_crit<<"\n";
    }
  } while (!pconst.equalD(*tf,t) && (err==0) && (ct<ctmax) && p1[lv_Hp]<i_crit);
  
  if (err || ct>=ctmax) {
    cerr<<"MP_HYDROGEN::Int_Adaptive_RKCK() errors encountered. nstep="<<ct<<"\n";
    rep.printVec("p1",p1,nvl);
    if (!err) err=ct;
#ifdef TESTING
      commandline.console("bad luck! >");
#endif
    
  }
  if (ct>0.75*ctmax) cout <<"ADAPTIVE INT: took "<<ct<<" steps!\n";
  for (int v=0;v<nvl;v++) pf[v] = p1[v];
  *tf = t;

  return err;
}


// ##################################################################
// ##################################################################



int MP_Hydrogen::dPdt(const int nv,    ///< number of variables we are expecting.
		      const double *P, ///< Current state vector.
		      double *R        ///< Rate Vector to write to.
		      )
{
  if (nv!=nvl) rep.error("variables wrong!",nv-nvl);

  for (int i=0;i<nvl; i++) R[i]=0.0;
  double temp=0.0;
  double T=0.0;

  T = (gamma-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.
  //cout <<"T="<<T;

  //
  // Check for Bad input state: 
  // (if nan/inf, set T to be negative to trigger large random values in rate)
  //
  for (int v=0;v<nvl; v++) 
    if (isnan(P[v]) || isinf(P[v])) T=-10.0;
  if (T<0.0 || P[lv_Hp]<0.0 || P[lv_Hp]>1.0) {
    //
    // want to reject this step, so choose random numbers to make sure
    // the integration doesn't converge.
    //
    //cout <<"Bad input to MP_Hydrogen::dPdt(): "<<T<<"\t"<<P[lv_Hp]<<"\n";
    for (int i=0;i<nvl; i++)
      R[i]=((double) rand())/((double) RAND_MAX);
    return 0;
  }


  if (ep.coll_ionisation) {
    temp  = coll_ion_rate(T)*(1.0-P[lv_Hp])*P[lv_Hp]*P[lv_nh];
    R[lv_Hp]   += temp; // coll.ion. to H+ adds to R[i]
    R[lv_eint] -= temp*coll_ion_energy(T); // reduces energy by the amount it took to ionise ion (i-1)
#ifdef COUNT_ENERGETICS
    if (!have_counted_ergs) {
      GLOBAL_CE->ci_rate = coll_ion_rate(T)*(1.0-P[lv_Hp])*P[lv_Hp]*P[lv_nh]*P[lv_nh];
      GLOBAL_CE->ci_cooling = GLOBAL_CE->ci_rate *coll_ion_energy(T);
    }
#endif //COUNT_ENERGETICS
  }
  //cout <<"  ci: rate="<<R[lv_Hp];

  if (ep.rad_recombination) {
    temp  = rad_recomb_rate(T)*P[lv_Hp]*P[lv_Hp]*P[lv_nh]; // rate [1/s]
    R[lv_Hp] -= temp; // recomb from i to i-1 reduces fraction.
#ifdef NO_DOUBLECOUNTING
    if (!ep.cooling) // avoid double counting!
#endif // NO_DOUBLECOUNTING
#ifdef HUMMER_RECOMB
      R[lv_eint] -= rad_recomb_energy(T) *P[lv_Hp]*P[lv_Hp]*P[lv_nh]; // rate [erg/s]
#ifdef TESTING
      if (dp.c->id==2773) {
	cout <<"\t\t*** energy rate after  recomb   ="<<R[lv_eint]<<"\n";
      }
#endif
#endif // HUMMER_RECOMB
#ifndef HUMMER_RECOMB
      R[lv_eint] -= temp *kB*T/(gamma-1.0);
    // also takes energy of e- out of gas (assumed radiated away)
    // if we are using a cooling curve, then we don't want to double count this!
    // This is an overestimate, maybe by a factor of 2 or more.
#endif // not HUMMER_RECOMB
#ifdef COUNT_ENERGETICS
      if (!have_counted_ergs) {
	GLOBAL_CE->rr_rate    = rad_recomb_rate(T)   *P[lv_Hp]*P[lv_Hp]*P[lv_nh]*P[lv_nh];
	GLOBAL_CE->rr_cooling = rad_recomb_energy(T) *P[lv_Hp]*P[lv_Hp]*P[lv_nh]*P[lv_nh];
      }
#endif //COUNT_ENERGETICS
  }
  //cout <<"  rr: rate="<<R[lv_Hp];
  //rep.printVec("\t\trr rate",R,nvl);

  if (ep.phot_ionisation) {
    R[lv_dtau] = P[lv_nh]*(1.0-P[lv_Hp])*phot_xsection(T)*path_length; // this is optical depth through cell.
    if (R[lv_dtau]>0.01)
      temp = photons_in*(1.0 -exp(-R[lv_dtau]))/P[lv_nh]/(1.0-P[lv_Hp]);
    else 
      temp = photons_in*phot_xsection(T)*path_length; // approx for 1-x<<1
    R[lv_Hp] += temp*(1.0-P[lv_Hp]);
    R[lv_eint] += temp*(1.0-P[lv_Hp])*phot_ion_energy(T); // this adds in X.XeV of energy per photo-ionisation.
#ifdef TESTING
    //if (dp.c->id==649090) {
    //  cout <<"photons_in="<<photons_in<<"\ttau="<<R[lv_dtau];
    //  cout <<"\tdxdt="<<R[lv_Hp]<<"\tdedt="<<R[lv_eint]<<"\n";
    //}
#endif
    R[lv_dtau] = exp(-R[lv_dtau]);
#ifdef COUNT_ENERGETICS
    if (!have_counted_ergs) {
      GLOBAL_CE->pi_rate    = temp*(1.0-P[lv_Hp])*P[lv_nh];
      GLOBAL_CE->pi_heating = temp*(1.0-P[lv_Hp])*P[lv_nh]*phot_ion_energy(T);
    }
#endif //COUNT_ENERGETICS
  }
  //  cout <<"  pi: rate="<<R[lv_Hp];
  //rep.printVec("\t\tpi rate",R,nvl);

  //if (P[lv_Hp]>1.0) {
    //cout<<"P = "<<P[lv_Hp]<<" and R = "<<R[lv_Hp]<<"\n";
    //R[lv_Hp]=min(0.0,R[lv_Hp]);
  //}
  //cout <<"  total: rate="<<R[lv_Hp]<<"\n";
  
  //#ifdef TESTING
  //  if (dp.c->id==10446) cout <<"pre- cooling R[eint]="<<R[lv_eint]<<"\n";
  //#endif
  if (ep.cooling) {
    temp = cool->CoolingRate(T,P[lv_Hp],P[lv_nh],FUV_unattenuated_flux,FUV_extinction) /P[lv_nh];

#ifdef TESTING
    if (dp.c->id==2773) {
      cout <<"\t\t*** energy gain/loss to cooling = "<<-temp<<"\n";
    }
#endif
    R[lv_eint] -= temp;
#ifdef COUNT_ENERGETICS
    if (!have_counted_ergs) {
      GLOBAL_CE->fn_cooling = temp*P[lv_nh];
    }
#endif //COUNT_ENERGETICS
  }
  //#ifdef TESTING
  //  if (dp.c->id==10446) cout <<"post-cooling R[eint]="<<R[lv_eint]<<"\n";
  //#endif

  R[lv_eint] *= P[lv_nh]; // To convert to energy per unit volume per unit time.
  
#ifdef COUNT_ENERGETICS
  if (!have_counted_ergs) {
    GLOBAL_CE->tot_heating = GLOBAL_CE->pi_heating;
    GLOBAL_CE->tot_cooling = GLOBAL_CE->rr_cooling +GLOBAL_CE->ci_cooling +GLOBAL_CE->fn_cooling;
    GLOBAL_CE->net_heating = GLOBAL_CE->tot_heating -GLOBAL_CE->tot_cooling;
    //cout <<"calculated heating="<<GLOBAL_CE->net_heating<<" and from dpdt()="<<R[lv_eint]<<"\n";
    GLOBAL_CE->cooling_time = P[lv_eint]/GLOBAL_CE->tot_cooling /pconst.s_per_yr();
    GLOBAL_CE->recomb_time  = P[lv_Hp]/GLOBAL_CE->rr_rate       /pconst.s_per_yr();
  }
#endif //COUNT_ENERGETICS

  if (!ep.update_erg) {
    //    cout <<"not updating energy.\n";
    R[lv_eint]=0.0;
  }


#ifdef ISOTHERMAL_MP
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
  //
  // Here I assume Tdot = (T1-T0)xdot
  // and use edot = nH*kB/(gamma-1) d/dt((1+x)(T0+(T1-T0)x))
  //   edot = nH*kB/(gamma-1)(T1*xdot+2(T1-T0)x*xdot)
  // This is appropriate for a gas temperature set by the ion fraction
  // according to T(x) = T0*(1-x) +T1*x   (Here T1=1.0e4, T0=1.0e2)
  // This is commonly used in Astrophysics, e.g. Williams, Ward-Thompson, Whitworth (2001)
  // 
  // Since gamma=1 I'm leaving out the (g-1) factor everywhere in the internal energy calcs.
  //
  //cout <<"isothermal rate!!!\n";
  R[lv_eint] = (P[lv_nh]*pconst.kB())*(1.0e4+1.98e4*P[lv_Hp])*R[lv_Hp];
  // **************************************************
  // ************** BIG ISOTHERMAL HACK!!!! ***********
  // **************************************************
#endif // ISOTHERMAL_MP

#ifdef TESTING
  if (dp.c->id==2773) {
    rep.printVec("\t\ttt rate",R,nvl);
  }
#endif

  return 0;
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::phot_xsection(double T)
{
  return 6.3e-18; // in cm^2 -- this is crude!
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::phot_ion_energy(double T)
{
  //
  // HACK HACK! ILIEV ET AL. 10^5K BLACKBODY WITH PHOTON ENERGY 16EV
  //
  //return 2.566e-11; // This is 16.0184eV for the Iliev et al test.
  //
  // HACK HACK! ILIEV ET AL. 10^5K BLACKBODY WITH PHOTON ENERGY 16EV
  //

  return 1.602e-12*5.0;  // 5.0 eV per ionisations.
  //  return 1.602e-12*2.4; // 2.4 eV per ionisation.
}



// ##################################################################
// ##################################################################




double MP_Hydrogen::coll_ion_rate(double T   ///< Precalculated Temperature.
				  )
{
  // This uses fitting formulae from Voronov (1997) ADANDT, 65, 1.
  // (Atomic Data And Nuclear Data Tables)
  // rate returned in cm^3/s
  if (T < 5.0e3) return 0.0;
  double A,X,K; int PP;
  PP=0; A=2.91e-8; X=0.232; K=0.39;
  //  cout <<"T="<<T<<" K,";
  T = ion_pot/kB/T;
  //  cout <<"\tPsi/kT="<<T<<" and CIR="<<A*(1.+PP*sqrt(T))*exp(K*log(T) -T)/(X+T)<<" cm^3/s\n";
  return A*(1.+PP*sqrt(T))*exp(K*log(T) -T)/(X+T);
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::coll_ion_energy(double T    ///< Precalculated Temperature.
				    )
{
  /** This is very simple -- each collisional ionisation removes the
   * ionisation energy from the system.  I don't know if it would ever
   * be more complicated, but it's in a function just in case.
   */
  return ion_pot; // This is already in ergs.
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::rad_recomb_rate(double T   ///< Precalculated Temperature.
				    )
{
#ifdef RT_TEST_PROBS
  return 2.59e-13; // set to const so can compare to analytic result
#endif //RT_TEST_PROBS

  double r;
#ifndef HUMMER_RECOMB
  // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
  // For case B recombination coefficient. (power law fit)
  r = 3.41202e-10*exp(-0.782991*log(T));
#endif // not HUMMER_RECOMB
#ifdef HUMMER_RECOMB
  if (T>5.0e6) T=5.0e6; // To avoid runoff at the end of spline fitting.
  if (T<10.0) T=10.0;
  interpolate.splint(hr_t, hr_alpha, hr_alpha2, hr_nspl, T, &r);
#endif // HUMMER_RECOMB
  return r;
}


// ##################################################################
// ##################################################################



double MP_Hydrogen::rad_recomb_energy(double T   ///< Precalculated Temperature.
				      )
{
#ifdef HUMMER_RECOMB
  /** This is the Case B energy loss rate in erg*cm^3/s. Fitting function from
   * Hummer, 1994, MNRAS, 268, 109.
   */
  double rate=0.0;
  if (T>5.0e6) T=5.0e6; // To avoid runoff at the end of spline fitting.
  if (T<10.0) T=10.0;
  interpolate.splint(hr_t, hr_beta, hr_beta2, hr_nspl, T, &rate);
  //  if (T>2.e4) cout <<"T="<<T<<" rate="<<rate*kB*T<<"\n";
  return rate*kB*T;
#endif // HUMMER_RECOMB
#ifndef HUMMER_RECOMB
  // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
  // For case B recombination coefficient. (power law fit), multiplied
  // by the mean electron energy (kT/(g-1)).  This is wrong, because
  // the energy should be integrated with the energy dependence of the 
  // rate, but it is wrong by less than a factor of 2.
  return 3.41202e-10*exp(-0.782991*log(T)) *kB*T/(gamma-1.0);
#endif // not HUMMER_RECOMB
}


// ##################################################################
// ##################################################################



///
/// This returns the minimum timescale of the times flagged in the
/// arguments.  Time is returned in seconds.
///
double MP_Hydrogen::timescales(const double *p_in,  ///< Current cell primitive vector.
			       const double gam,    ///< EOS gamma.
			       const bool f_cool,   ///< set to true if including cooling time.
			       const bool f_recomb, ///< set to true if including recombination time.
			       const bool f_photoion ///< set to true if including photo-ionsation time.
			       )
{
  //
  // First convert to MP variables and get temperature.
  //
  int err = 0;
  //  cout <<"gamma: "<<gamma<<"\n";
  double P[nvl];
  err += convert_prim2local(p_in,P,gam);
  double T = (gam-1.)*P[lv_eint]/kB/(1.0+P[lv_Hp])/P[lv_nh]; // Temperature.

  //
  // Now get the timescale for each requested process.
  // Note this is does not include radiative transfer for now, so the fluxes
  // are set to zero, and DON'T REQUEST THE PHOTO-IONISATION TIME.  That is 
  // a wish-list feature for the future.
  //
  double mintime=1.0e99;
  double rate=0.0;

  //
  // we skip the cooling time if the temperature is already very low;
  // a reasonable temperature could be 10K; there is unlikely to be
  // significant cooling below this temperature, and even if there is
  // it won't be by a large fraction.
  //
  if ((ep.cooling) && (f_cool) && (T>=10.0) ) {
    //
    // This gets the cooling rate per unit volume, first from the
    // cooling function, and then adding on the recombination cooling.
    //
    rate = cool->CoolingRate(T,P[lv_Hp],P[lv_nh], 0.0,0.0);
    if (ep.rad_recombination)
      rate += rad_recomb_energy(T) *P[lv_Hp]*P[lv_Hp]*P[lv_nh]*P[lv_nh];
    //
    // Cooling time is then the energy per unit volume divided by the rate.
    // Bear in mind we could have net heating, so take fabs(rate)
    //
    mintime = min(mintime, P[lv_eint]/fabs(rate));
  }

  //
  // We skip the recombination time if gas is neutral; the cutoff
  // point is arbitrary, but setting it to 0.03 seems quite
  // reasonable, and sufficiently conservative.
  //
  if ((ep.rad_recombination) && (f_recomb) && (P[lv_Hp]>=0.03) ) {
    //
    // Get the radiative recombination rate per ION.
    //
    rate = rad_recomb_rate(T) *P[lv_Hp]*P[lv_nh];
    //
    // The recombination time is then the reciprocal of the rate.
    // The timestep will then be multiplied by 0.3 in gridMethods,
    // to give a total MP step of 0.1*t_rec (0.33333*0.3/RRR)
    //
    mintime = min(mintime, 0.33333/fabs(rate));
  }

#ifdef RT_TEST_PROBS
  if ((f_recomb)) {
    //
    // Get the radiative recombination rate per ATOM here (may not be ions).
    // Temperature is irrelevant for test probs b/c rrr is set to 2.59e-13.
    //
    rate = rad_recomb_rate(1.0e4) *P[lv_nh];
    //
    // The recombination time is then the reciprocal of the rate.
    //
    mintime = min(mintime, 0.33333/fabs(rate));
    //cout <<"mintime="<<mintime<<"\n";
  }
#endif

  if ((ep.phot_ionisation) && (f_photoion)) {
    rep.error("Haven't coded for photo-ionisation time!!!",f_photoion);
  }

  return mintime;
}


// ##################################################################
// ##################################################################




#endif // if not excluding MPv1

