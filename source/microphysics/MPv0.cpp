/// \file microphysics_v1.cc
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
/// - JM 2009-12-16 Fixed microphysics for "only cooling" -- it was broken.
///   Added in an ifdef in MP_Hydrogen for if the integration (not
///   Photoionisation) returns an ifrac>2, which would indicate something went
///   very wrong.  In that case I try splitting the integral into 10
///   subintegrals and retrying.
//   If that fails, bug out.
/// - 2010-01-05 JM: Added in (public) function which returns timescales for
/// heating/cooling
///    and recombination/ionisation etc.
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
/// - 2010.10.11 JM: Split microphysics.cc into two -- this file has
///    the original class, which doesn't work for radiative transfer.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2012.01.20-26 JM: wrapped code in ifndef so that its compilation can
///    be disabled, since I never use it now.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2015.08.05 JM: tidied up code; added pion_flt datatype.
/// - 2018.03.20 JM: Renamed class to MPv0.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#ifdef LEGACY_CODE

#include "microphysics/MPv0.h"
using namespace std;

// ##################################################################
// ##################################################################

MPv0::MPv0(
    const int nv,
    const int ntr,
    const std::string chem_code,  ///< type of chemistry we are running.
    const std::string *tracers,   ///< List of what the tracer variables mean.
    struct which_physics *ephys,  ///< pointer to extra physics flags.
    struct rad_sources *rsrcs     ///< radiation sources.
    ) :
    microphysics_base(nv, ntr, tracers, ephys, rsrcs),
    kB(pconst.kB()),
    m_p(pconst.m_p())
{
  cout << "\t\tMPv0 constructor.\n";
  cout << "WARNING: THIS CODE HAS NEVER BEEN USED IN PUBLICATIONS.\n";
  cout << "IT SHOULD BE CONSIDERED UNTESTED AND UNRELIABLE.\n";
  set_atomic_data();  // sets atomic data arrays to right values.
  //  cout <<"\t\tAtomic data set.\n";
  min_elecf = 1.e-5;  // minimum electron fraction (to seed reactions!).

  // set flags for what processes we are and aren't doing.
  ep.dynamics          = ephys->dynamics;
  ep.cooling           = ephys->cooling;
  ep.chemistry         = ephys->chemistry;
  ep.coll_ionisation   = ephys->coll_ionisation;
  ep.rad_recombination = ephys->rad_recombination;
  ep.phot_ionisation   = ephys->phot_ionisation;
  ep.raytracing        = ephys->raytracing;
  ep.update_erg        = ephys->update_erg;
  if (!ep.chemistry) {
    ep.coll_ionisation = ep.rad_recombination = ep.phot_ionisation = 0;
  }
  //  cout <<"\t\tExtra Physics flags set.\n";

  if (!ep.cooling) {
    cout << "\t\t cooling not needed.\n";
    cool = 0;
  }
  else {
    cool = 0;
    cool = new CoolingFn(ep.cooling);
    if (!cool) rep.error("CoolingFn init", cool);
  }
  //  cout <<"\t\tCooling Function set up.\n";

  // Set up tracer variables.
  cout << "\t\tSetting up Tracer Variables.  Assuming tracers are last " << ntr
       << " variables in state vec.\n";
  int ftr = nv_prim - ntr;  // first tracer variable.
  string s;

  int len = ntr;

  MPv0::lvar["n_h"] =
      0;  // 1st element of local vector is hydrogen number density.
  MPv0::lvar["Eint"] =
      1;  // Second element of local state vector is internal energy/vol.
  MPv0::lv_nh   = 0;
  MPv0::lv_eint = 1;
  MPv0::lv_elec = -1;
  MPv0::pv_elec = -1;
  int firstion = 2;  // location in p-vec of first ion/electron tracer variable.
  // if we are doing photo-ionisation, we want to track the change in optical
  // depth over a substep.
  if (ep.phot_ionisation) {
    lvar["dtau"]  = 2;
    MPv0::lv_dtau = 2;
    firstion += 1;
  }
  int ct = 0, colour = 0;
  MPv0::nions = 0;
  MPv0::nels  = 0;

  for (int i = 0; i < len; i++) {
    // Now pick out the chemistry tracers and pass to microphysics
    // constructor
    s = tracers[i];  // Get 'i'th tracer variable.

    if (s == "e-____") {
      s       = "e-";
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      pv_elec = pvar[s];
      lv_elec = lvar[s];
      ct++;
    }
    else if (s == "trace_" || s == "colour" || s == "Colour") {
      s = "tr";
      colour++;  // do nothing with non-chemistry tracers.
    }
    else if (s == "H1+___" || s == "HII__") {
      s       = "H1+";
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      lv_hp = lvar[s];
      ions.push_back(s);
      nions++;
      els.push_back("H");
      nels++;
    }
    else if (s == "He1+__" || s == "HeII_") {
      s       = "He1+";
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      ions.push_back(s);
      nions++;
      els.push_back("He");
      nels++;
    }
    else if (s == "He2+__" || s == "HEIII") {
      s       = "He2+";
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      ions.push_back(s);
      nions++;
      bool fff = true;
      for (int j = 0; j < nels; j++)
        if (els[j] == "He") fff = false;
      if (fff) {
        els.push_back("He");
        nels++;
      }  // unlikely, but you never know...
    }
    else if (s.substr(0, 1) == "C") {
      bool fff = true;
      for (int j = 0; j < nels; j++)
        if (els[j] == "C") fff = false;
      if (fff) {
        els.push_back("C");
        nels++;
      }
      if (s.substr(0, 2) == "C0")
        s = "C0";
      else
        s = s.substr(0, 3);
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      ions.push_back(s);
      nions++;
    }
    else if (s.substr(0, 1) == "N") {
      bool fff = true;
      for (int j = 0; j < nels; j++)
        if (els[j] == "N") fff = false;
      if (fff) {
        els.push_back("N");
        nels++;
      }
      if (s.substr(0, 2) == "N0")
        s = "N0";
      else
        s = s.substr(0, 3);
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      ions.push_back(s);
      nions++;
    }
    else if (s.substr(0, 1) == "O") {
      bool fff = true;
      for (int j = 0; j < nels; j++)
        if (els[j] == "O") fff = false;
      if (fff) {
        els.push_back("O");
        nels++;
      }
      if (s.substr(0, 2) == "O0")
        s = "O0";
      else
        s = s.substr(0, 3);
      pvar[s] = ftr + i;
      lvar[s] = firstion + ct;
      ct++;
      ions.push_back(s);
      nions++;
    }
    else
      rep.error("Bad Tracer type!", s);
  }
  if ((ct + colour) != len)
    rep.error("didn't get all tracers", len - ct - colour);
  if (pvar.size() != static_cast<unsigned int>(ct))
    rep.error(
        "tracers init error", pvar.size() - static_cast<unsigned int>(ct));
  //  cout <<"\t\tset up tracers\n";

  MPv0::nvl = lvar.size();
  Integrator_Base::Set_Nvar(nvl);
  //  cout <<"\t\tset int_nvar\n";

  // Now initialize chemistry class.
  MPv0::chemtype = chem_code;

  if (chemtype == "color" || chemtype == "colour") chemtype = "None";
  // this has tracers, but none relating to chemistry
  else if (chemtype == "ChAH__") {  // Legacy Code...
    //    ions.push_back("H1+"); // only ion is H 1+
    //    els.push_back("H"); // only element is H
    if (nions != 1) {
      cerr << "Nions setup wrong.\n";
      nions = 1;
    }
    if (nels != 1) {
      nels = 1;
    }
  }
  else if (chemtype == "ChALL") {
    // do nothing.
  }
  else if (chemtype == "None" || chemtype == "NONE") {
    // no chemistry, so better be doing cooling, otherwise no point
    // existing...
    if (!cool)
      rep.error("set up MP but no microphysics and not doing cooling!", nions);
    if (nions != 0) {
      cout << "nions = " << nions << ", setting to zero.\n";
      nions = 0;
    }
    ii   = 0;
    ee   = 0;
    nels = 0;
    return;
  }
  else
    rep.error("No known chemistry defined.", chemtype);

  ii = 0;
  ee = 0;
  ii = new ion_struct *[nions];
  ee = new element_struct *[nels];
  if (!ii || !ee) rep.error("MPv0 init, ii/ee malloc", ii);

  for (int i = 0; i < nions; i++) {
    ii[i]        = &ion_props[ions[i]];
    ii[i]->index = lvar[ions[i]];  // set ion's index in the local state vector.
    ii[i]->pv_index = pvar[ions[i]];
  }
  for (int i = 0; i < nels; i++) {
    ee[i] = &el_props[els[i]];
    // now set the indices of the element's ions in the local state vector.
    ct = 0;
    // count number of ions we are using which are associated with the
    // current element.
    for (unsigned int v = 0; v < ee[i]->ions.size(); v++) {
      for (int vv = 0; vv < nions; vv++) {
        //	cout <<"ions: "<<ions[vv]<<" and els:
        //"<<ee[i]->ions[v]<<"\n";
        if (ions[vv] == ee[i]->ions[v]) ct++;
      }
    }
    // now add them if needed, and zero indices otherwise.
    if (ct > 0) {
      ee[i]->nspecies = ct;
      if (ee[i]->ion_indices) {
        delete[] ee[i]->ion_indices;
        ee[i]->ion_indices = 0;
      }
      ee[i]->ion_indices = new int[ct];
      if (!ee[i]->ion_indices) rep.error("meminit MPv0()", 0);
      ct = 0;
      for (unsigned int v = 0; v < ee[i]->ions.size(); v++) {
        for (int vv = 0; vv < nions; vv++)
          if (ions[vv] == ee[i]->ions[v]) {
            cout << "\t\tFor element " << ee[i]->el << " adding ion "
                 << ions[vv] << " with index " << ii[vv]->index << "\n";
            ee[i]->ion_indices[ct] = ii[vv]->index;
            ct++;
          }
      }
    }
    else {
      cout << "!!!!!!!!!!!!-- No ions found for element " << ee[i]->el << "\n";
      if (ee[i]->ion_indices) {
        delete[] ee[i]->ion_indices;
        ee[i]->ion_indices = 0;
      }
      ee[i]->nspecies = 0;
    }
  }

  // finally, we may not be including all ionisation stages of an element, so
  // we need to maybe set ip1 and im1 to zero...
  // At the moment I assume that im1 should always be tracked, as I don't
  // explicitly include the neutral stage as a variable usually, but I do want
  // to track it.
  cout << "setting higher/lower stages.\n";
  for (int i = 0; i < nions; i++) {
    bool higherstage = false, lowerstage = false;
    for (int j = 0; j < nions; j++) {
      if (ii[i]->ip1 == ii[j]->ion) higherstage = true;
      if (ii[i]->im1 == ii[j]->ion) lowerstage = true;
    }
    if (!higherstage) {
      ii[i]->ip1  = "";
      ii[i]->iip1 = 0;
    }
    else {
      ii[i]->iip1 = &ion_props[ii[i]->ip1];
    }
    if (!lowerstage) {
      // ii[i]->im1 = ""; ii[i]->iim1 = 0;
      if (ii[i]->im1 != "")
        ii[i]->iim1 = &ion_props[ii[i]->im1];
      else
        ii[i]->iim1 = 0;
    }
    else {
      ii[i]->iim1 = &ion_props[ii[i]->im1];
    }
  }
  cout << "setting higher/lower stages.\n";

  for (int i = 0; i < nions; i++)
    show_ion_struct(ii[i]);
  /*
  for (int i=0;i<ion.size(); i++) cout <<" ion: "<<ion[i];
  cout<<"\n";
  for (int i=0;i<el.size(); i++) cout <<" el: "<<el[i];
  cout<<"\n";
  for (int i=0;i<ion.size(); i++) cout <<" i2el: "<<i2el[ion[i]];
  cout<<"\n";
  for (int i=0;i<el.size(); i++) cout <<" numfrac: "<<numfrac[el[i]];
  cout<<"\n";
  */
  //  cout <<"\t\tinit done.\n";
  return;
}

// ##################################################################
// ##################################################################

void MPv0::copy_ion_struct(const ion_struct src, ion_struct *dest)
{
  rep.error("unused function. test me!", -99);
  dest->ion     = src.ion;
  dest->ip1     = src.ip1;
  dest->im1     = src.im1;
  dest->el      = src.el;
  dest->charge  = src.charge;
  dest->ion_pot = src.ion_pot;
  dest->i       = src.i;
  return;
}

// ##################################################################
// ##################################################################

void MPv0::copy_element_struct(
    const struct element_struct src, struct element_struct *dest)
{
  rep.error("unused function. test me!", -99);
  dest->el       = src.el;
  dest->nspecies = src.nspecies;
  dest->mass     = src.mass;
  dest->numfrac  = src.numfrac;
  if (!dest->ions.empty()) dest->ions.clear();
  for (unsigned int i = 0; i < src.ions.size(); i++)
    dest->ions.push_back(src.ions[i]);
  if (dest->ion_indices) {
    delete[] dest->ion_indices;
    dest->ion_indices = 0;
  }
  if (dest->nspecies > 0) dest->ion_indices = new int[dest->nspecies];
  if (!dest->ion_indices)
    rep.error("meminit copy_element_struct()", dest->ion_indices);
  for (int i = 0; i < dest->nspecies; i++)
    dest->ion_indices[i] = src.ion_indices[i];
  return;
}

// ##################################################################
// ##################################################################

void MPv0::show_ion_struct(const ion_struct *i)
{
  cout << "\t\tion: " << i->ion << "\tel: " << i->el << " index:" << i->index
       << " pv_index:" << i->pv_index;
  cout << "\tpot:" << i->ion_pot << "  charge:" << i->charge << "\n";
  cout << "\t\t\tip1: " << i->ip1 << " im1: " << i->im1 << " el:" << i->el
       << "\n";
  cout << "\t\t\tiip1: " << i->iip1 << "  iim1: " << i->iim1 << "\n";
  return;
}

// ##################################################################
// ##################################################################

MPv0::~MPv0()
{
  if (cool) {
    delete cool;
    cool = 0;
  }
  if (ii) {
    delete[] ii;
    ii = 0;
  }
  if (ee) {
    for (int i = 0; i < nels; i++) {
      cout << "deleting " << i << "\n";
      delete[] ee[i]->ion_indices;
    }
    delete[] ee;
    ee = 0;
  }
  ion_list.clear();
  el_list.clear();
}

// ##################################################################
// ##################################################################

double MPv0::Get_nH(const double rho  ///< gas density.
)
{
  double nh = 0.0;
  for (int i = 0; i < nels; i++)
    nh += ee[i]->numfrac * ee[i]->mass;
  nh = rho / m_p / nh;
  return nh;
}

// ##################################################################
// ##################################################################

double MPv0::Get_Temp(const double *P  ///< state vector
)
{
  //  rep.printVec("\t\t\tPPPPPPPPPPPP",P,nvl);
  //  cout <<"Get_nTot(P): "<<Get_nTot(P)<<"\n";
  return (gamma - 1.) * P[lv_eint] / kB / Get_nTot(P);
}

// ##################################################################
// ##################################################################

int MPv0::Set_Eint(
    double *P,      ///< state vector
    const double T  ///< Temperature we want to set to.
)
{
  if (T <= 0) {
    cout << "negative temperature in Set_Eint().\n";
    return 1;
  }
  P[lv_eint] = Get_nTot(P) * kB * T / (gamma - 1.);
  return 0;
}

// ##################################################################
// ##################################################################

double MPv0::Get_nTot(const double *P  ///< state vector
)
{
  // This assumes n_h has already been set, as that is the first thing
  // TimeUpdateMP does.
  double nf = P[lv_elec];
  for (int i = 0; i < nels; i++)
    nf += ee[i]->numfrac;
  return P[lv_nh] * nf;
}

// ##################################################################
// ##################################################################

double MPv0::Get_nIons(const double *P  ///< state vector
)
{
  double nf = 0.0;
  for (int i = 0; i < nions; i++)
    nf += ii[i]->e->numfrac * P[ii[i]->index];
  return P[lv_nh] * nf;
}

// ##################################################################
// ##################################################################

double MPv0::neutral_fraction(const double *P, const ion_struct *i)
{
  //  rep.printVec("nf vec",P,nvl);
  double nf         = 1.0;
  element_struct *e = i->e;
  //  cout <<"nions: "<<e->ions.size()<<" and first el="<<e->ions[1]<<"\n";
  for (int j = 0; j < e->nspecies; j++) {
    // cout <<"j="<<j<<"\te->ion_indices[j] = "<< e->ion_indices[j]<<"\t:
    // "<<P[e->ion_indices[j]]<<"\n";
    nf -= P[e->ion_indices[j]];
  }
  //  for (int j=0; j<e->nspecies; j++) nf -= P[lvar[e->ions[j+1]]];
  if (nf < 0.) {
    // rep.error("neutral_fraction() too small",nf);
    //    cout <<"neutral fraction negative: "<<nf<<"   correcting."<<"\n";
    nf = 0.0;
  }
  else if (nf > 1.0) {
    //    rep.printVec("P",P,nvl);
    //    rep.error("neutral_fraction() too large",nf);
    nf = 1.0;
  }
  return nf;
}

// ##################################################################
// ##################################################################

int MPv0::Tr(string t)
{
  if (pvar.find(t) == pvar.end())
    return -1;
  else
    return pvar[t];
}

// ##################################################################
// ##################################################################

int MPv0::Set_Temp(
    pion_flt *p,     ///< primitive vector.
    const double T,  ///< temperature.
    const double g   ///< eos gamma.
)
{
  gamma = g;
  double P[nvl];
  int err = convert_prim2local(p, P, g);
  err += Set_Eint(P, T);
  err += convert_local2prim(P, p, p, g);
  return err;
}

// ##################################################################
// ##################################################################

double MPv0::Temperature(
    const pion_flt *p,  ///< primitive vector
    const double g      ///< eos gamma
)
{
  if (nions > 0) {
    gamma = g;
    double P[nvl];
    convert_prim2local(p, P, g);
    return Get_Temp(P);
  }
  else
    return p[PG] * m_p * 0.61 / kB
           / p[RO];  // assumes fully ionised gas with mu=1.22
}

// ##################################################################
// ##################################################################

int MPv0::convert_prim2local(
    const pion_flt *p_in, double *p_local, const double gam)
{
  /** \todo Speed up this function by removing string maps!!! e.g. set up
   * a map<int,int> for relating primitive indices to local indices, and vice
   * versa. this would be much faster.
   */
  int xi         = 0;  // index for current ion.
  p_local[lv_nh] = Get_nH(p_in[RO]);
  // cout <<"nh = "<<p_local[lv_nh]<<"\n";
  // This is generic for all ideal gas EOSs, but it would be nice to have a
  // pointer to the equations.
  p_local[lv_eint] = p_in[PG] / (gam - 1.);

  if (nions > 0) {
    //    cout <<"setting electrion fraction\n";
    p_local[lv_elec] = p_in[pv_elec];
    p_local[lv_elec] = max(p_local[lv_elec], min_elecf);
    for (int i = 0; i < nions; i++) {
      // loop over all ions.
      // cout <<"setting ion fraction.\n";
      xi = ii[i]->index;  // current ion index in local array.
      if (xi >= nvl) rep.error("more ions than variables", xi - nvl);
      p_local[xi] = p_in[ii[i]->pv_index];

      if (p_local[xi] > 1.01) {
        rep.printVec("P", p_local, nvl);
        cerr << "ion fraction >1! for ion " << ii[i]->ion << " at value "
             << p_local[xi] - 1. << "\n";
        cout << "neutral fraction: " << neutral_fraction(p_local, ii[i])
             << "\n";
        ;
        //	return 1;
        // p_local[xi] = 1.0;
      }
      p_local[xi] = std::max(p_local[xi], 0.0);
      p_local[xi] = std::min(p_local[xi], 1.0);
    }
  }
  else if (nions == 0) {
    // no chemical species, not even electron fraction.
  }
  else {
    rep.error("ni<0 in MPv0::TimeUpdateMP", nions);
  }

  if (ep.phot_ionisation) {
    p_local[lv_dtau] =
        0.0;  // p_local[lv_nh]*(1.0-p_local[lv_hp])*phot_xsection(ii[i])*path_length;
    tau_cell = 0.0;
  }

  return 0;
}

// ##################################################################
// ##################################################################

int MPv0::convert_local2prim(
    const double *p_local,
    const pion_flt *p_in,
    pion_flt *p_out,
    const double gam)
{
  // Now return updated values, note density, velocity, mag.field don't
  // change. only variables affected are pressure and tracers.
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];
  p_out[PG]     = p_local[lv_eint] * (gam - 1.0);
  double min_if = 1.0e-12;
  if (nions > 0) {
    p_out[pv_elec] = max(p_local[lv_elec], min_elecf);
    for (int i = 0; i < nions; i++)
      p_out[ii[i]->pv_index] = max(p_local[ii[i]->index], min_if);
  }
  if (ep.phot_ionisation) {
    tau_cell = p_local[lv_dtau];  // this should be int(tau,dt)
  }
  return 0;
}

// ##################################################################
// ##################################################################

int MPv0::TimeUpdate_OnlyCooling(
    const pion_flt *p_in,
    pion_flt *p_out,
    const double dt,
    const double g,
    const int sw_int,
    double *ttt)
{
  //  cout <<"only cooling: nvl = "<<nvl<<"\n";
  if (nvl != 2) rep.error("nvl wrong", nvl);
  int err = 0;
  gamma   = g;

  // set up local state vector.
  double P[nvl];
  P[lv_nh]   = p_in[RO] / m_p / 1.22;  // roughly cosmic abundances give mu=1.22
  P[lv_eint] = p_in[PG] / (gamma - 1.);

  /*
  if(p_in[PG]>1.e-7) {
    rep.printVec("start P",P,nvl);
    rep.printVec("p in   ",p_in,nv_prim);
  }
  */

  double errtol = 1.0e-3;
  double tout   = 0.0;
  if (sw_int == 0) {
    errtol = 1.0e-3;
    err    = Int_Adaptive_RKCK(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int == 1) {
    errtol = 1.0e-3;
    err    = Int_DumbAdaptive_Euler(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int == 2) {
    // DANGEROUS!!!
    // cout <<"Requesting single step RK4 integration for MPv0! Are you sure
    // about this?\n";
    err = Step_RK4(nvl, P, 0.0, dt, P);
  }
  else
    rep.error("this integration method not known.", sw_int);
  if (err) rep.error("integration failed.", err);
  //  rep.printVec("p_new",P,nvl);

  *ttt = P[lv_eint] * (gamma - 1.0) / kB
         / (2.0 * P[lv_nh]);  // fully ionised hydrogen
  // set output state vector -- only thing that changes is the pressure due to
  // cooling.
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];
  p_out[PG] = P[lv_eint] * (gamma - 1.0);

  /*
  if(p_out[PG]>1.e-8) {
    rep.printVec("large P",P,nvl);
    rep.printVec("p in   ",p_in,nv_prim);
    rep.printVec("p out  ",p_out,nv_prim);
  }
  */
  return 0;
}

// ##################################################################
// ##################################################################

int MPv0::TimeUpdateMP(
    const pion_flt *p_in,
    pion_flt *p_out,
    const double dt,
    const double g,
    const int sw_int,
    double *ttt)
{
  // First if we are just doing cooling with no microphysics, call different
  // function
  if (nions == 0 && ep.cooling)
    return TimeUpdate_OnlyCooling(p_in, p_out, dt, g, sw_int, ttt);

  // else do the regular update:
  int err = 0;
  gamma   = g;
  double P[nvl];  // local state vector for microphysics = [n_h,Eint,x_e,(x_i)]
  // put data from p_in into local vector.
  err += convert_prim2local(p_in, P, gamma);

  // cout <<"start  T= "<<Get_Temp(P) <<" and n_atoms=
  // "<<Get_nTot(P)-P[lv_nh]*P[lv_elec]<<"\n";
  //  rep.printVec("p_old",P,nvl);
  double errtol = 1.0e-5;
  double tout   = 0.0;
  //  double perr[nvl];//,P2[nvl];
  if (sw_int == 0) {
    errtol = 1.0e-8;
    err    = Int_Adaptive_RKCK(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int == 1) {
    errtol = 1.0e-3;
    err    = Int_DumbAdaptive_Euler(nvl, P, 0.0, dt, errtol, P, &tout);
  }
  else if (sw_int == 2) {
    // DANGEROUS!!!
    cout << "Requesting single step RK4 integration for MPv0! Are you sure "
            "about this?\n";
    err = Step_RK4(nvl, P, 0.0, dt, P);
  }
  else
    rep.error("this integration method not known.", sw_int);
  if (err) rep.error("integration failed.", err);
  //  rep.printVec("p_new",P,nvl);

  int xi = 0;
  for (int i = 0; i < nions; i++) {
    xi = ii[i]->index;
    if (P[xi] >= 1.0001) {
      cout << "ion: " << ii[i]->ion << " has i-frac=" << P[xi]
           << "  ...setting to 1.\n";
      P[xi] = 1.0;
    }
  }
  if (P[lv_elec] > 1.0001) {
    //    cout <<"elec: "<<P[lv_elec];
    P[lv_elec] = 0.0;
    for (int i = 0; i < nions; i++) {
      P[lv_elec] += P[ii[i]->index] * ii[i]->charge * ii[i]->e->numfrac;
    }
    P[lv_elec] = max(min_elecf, P[lv_elec]);
    //    cout <<"  new: "<<P[lv_elec]<<"\n";
  }
  //  cout <<"finish T= "<<Get_Temp(P) <<" and n_atoms=
  //  "<<Get_nTot(P)-P[lv_nh]*P[lv_elec]<<"\n";
  // put updated state vector into p_out.
  *ttt = Get_Temp(P);
  // if (*ttt >1.e4) cout <<"\tT="<<*ttt<<" and n_atoms=
  // "<<Get_nTot(P)-P[lv_nh]*P[lv_elec]<<"\n";
  err += convert_local2prim(P, p_in, p_out, gamma);
  return err;
}

// ##################################################################
// ##################################################################

int MPv0::TimeUpdate_RTsinglesrc(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    pion_flt *p_out,       ///< Destination Vector for updated values.
    const double dt,       ///< Time Step to advance by.
    const double g,        ///< EOS gamma.
    const int sw_int,      ///< Switch for what type of integration to use.
                       ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
    const double phot_in,  ///< flux in per unit length along ray (F/ds or L/dV)
    const double ds,       ///< path length ds through cell.
    const double tau2cell,  ///< Optical depth to entry point of ray into cell.
    double *deltau  ///< return optical depth through cell in this variable.
)
{
  if (!ep.phot_ionisation)
    rep.error("RT requested, but phot_ionisation not set!", ep.phot_ionisation);
  MPv0::tau_cell = 0.0;
  MPv0::photons_in =
      phot_in;  // units are photons/cm^3/s (/Hz if using frequency info).

  //
  // with new interface, the flux is passed in unattenuated
  //
  photons_in *= exp(-tau2cell);
  //
  // Now we calculate the FUV flux and FUV optical depth.
  // This is based on equation A3 in Henney, Arthur, et al., 2009.
  //
  double f_UV = 0.5;  // number of UV photons (6-13eV) per ionising photon.
  FUV_unattenuated_flux =
      phot_in * f_UV * ds / 1.2e7;  // Flux in units of "Habing Flux"
                                    // according to Henney et al. 2009.
  FUV_extinction = 1.086 * 5.0e-22 * tau2cell
                   / 6.3e-18;  // Extinction from Henney et al. 2009.
  FUV_attenuated_flux = FUV_unattenuated_flux * exp(-1.9 * FUV_extinction);
  rep.error("New interface is only for MP_Hydrogen so far", -1234567);

  MPv0::path_length = ds;
  double temp       = 0.0;
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
  if (photons_in * path_length < 1.0) {
    photons_in = 0.0;  // ep.phot_ionisation=0;
  }
  //  else cout <<"\tphot_in: "<<photons_in<<"  path: "<<path_length<<" \n";
  int err = TimeUpdateMP(p_in, p_out, dt, g, sw_int, &temp);
  *deltau = tau_cell / dt;  // tau_cell is the integral of optical depth
                            // through cell over timestep.
  // deltau is now time-averaged optical depth through cell.
  ep.phot_ionisation = 1;
  return err;
}

// ##################################################################
// ##################################################################

int MPv0::dPdt_OnlyCooling(
    const int nv,     ///< number of variables we are expecting.
    const double *P,  ///< Current state vector.
    double *R         ///< Rate Vector to write to.
)
{
  R[lv_nh] = 0.0;
  //  R[lv_eint] = ne*ni*Lambda(T)
  double temp = (gamma - 1.) * P[lv_eint] / kB
                / (2.0 * P[lv_nh]);  // for fully ionised hydrogen.
  R[lv_eint] = -cool->CoolingRate(
      temp, 1.0, P[lv_nh], FUV_unattenuated_flux, FUV_extinction);
  if (R[lv_eint] > 0) {
    cout << "rate positive!: " << R[lv_eint] << "\n";
    R[lv_eint] = 0;
  }
  // R[lv_eint] /= 10.0; ///??? What's this for????
  return 0;
}

// ##################################################################
// ##################################################################

int MPv0::dPdt(
    const int nv,     ///< number of variables we are expecting.
    const double *P,  ///< Current state vector.
    double *R         ///< Rate Vector to write to.
)
{
  if (nv != nvl) rep.error("variables wrong!", nv - nvl);

  // first if just cooling, call a different function
  if (nions == 0) return dPdt_OnlyCooling(nv, P, R);

  // else do the microphysics routine:
  int xi = 0;  // index for current ion.
  for (int i = 0; i < nvl; i++)
    R[i] = 0.0;
  double temp        = 0.0;
  double temperature = 0.0;

  if (nions < 0) rep.error("nions<0 in MPv0::dPdt", nions);

  struct ion_struct *iip1 = 0, *iim1 = 0;
  temperature = Get_Temp(P);
  for (int i = 0; i < nions; i++) {
    iip1 = ii[i]->iip1;
    iim1 = ii[i]->iim1;

    xi   = ii[i]->index;
    temp = 0.0;
    //        cout <<"ion :"<<ions[i]<<" temp="<<temperature<<" "<<lv_nh<<"
    //        "<<lv_eint<<" "<<lv_elec<<" "<<xi<<"\n"; cout <<"
    //        temp="<<temperature<<"\n";
    //  rep.printVec("\t\tP",P,nvl);

    if (ep.coll_ionisation) {
      if (iim1) {
        // ionisation from lower stage.
        temp = Coll_Ion_rate(temperature, iim1) * P[lv_elec] * P[lv_nh];
        if (iim1->index < 0)
          temp *= neutral_fraction(P, iim1);
        else
          temp *= P[iim1->index];
        R[xi] += temp;  // coll.ion. to ion i adds to R[i]
        R[lv_eint] -= temp * iim1->ion_pot
                      * ii[i]->e->numfrac;  // reduces energy by the amount it
                                            // took to ionise ion (i-1)
      }
      if (iip1) {
        // ionisation to higher stage, takes ions out.
        temp = Coll_Ion_rate(temperature, ii[i]) * P[lv_elec] * P[lv_nh];
        temp *= P[xi];
        R[xi] -= temp;  // lose ions at stage i to stage i+1
                        // energy loss already counted by lower stage
      }
    }
    // rep.printVec("\t\tci rate",R,nvl);

    /*
    // charge exchange
    if (iim1 && iim1->ch_ex) {
      // ch.ex. ionisation of iim1 adds to the rate R(ch-ex)*n_(H0)*x(iim1)
    [per sec] temp =
    charge_exchange_rate(temperature,iim1,1)*(1.0-P[lv_hp])*P[iim1->index]*P[lv_nh];
      R[xi] += temp;
    }
    if (ii[i]->ch_ex) {
      // can have recombination or ionisation exchange of current ion with
    H0 or H+ R[xi] +=
    charge_exchange_rate(temperature,ii[i],1)*(1.0-P[lv_hp])*P[xi]*P[lv_nh];
      R[xi] -= charge_exchange_rate(temperature,ii[i],0)*     P[lv_hp]
    *P[xi]*P[lv_nh];
    }
    if (iip1 && iip1->ch_ex) {
      // ch.ex. recombination of iip1 adds to the rate
    R(ch-ex)*n_(H+)*x(iip1) [per sec] R[xi] -=
    charge_exchange_rate(temperature,iip1,0)*      P[lv_hp]
    *P[iip1->index]*P[lv_nh];
    }
    */

    if (ep.rad_recombination) {
      if (iim1) {
        // must have lower stage, but checked anyway.
        temp = Rad_Recomb_rate(temperature, ii[i]) * P[lv_elec]
               * P[lv_nh];      // rate [1/s]
        R[xi] -= temp * P[xi];  // recomb from i to i-1 reduces fraction.
        // if (!ep.cooling || temperature<1.0e4) // avoid double
        // counting!
        //
        // This is a bad approximation!!! UPDATE ME!!!
        //
        R[lv_eint] -=
            temp * P[xi] * ii[i]->e->numfrac * kB * temperature / (gamma - 1.0);
        // also takes energy of e- out of gas (assumed radiated away)
        // if we are using a cooling curve, then we don't want to double
        // count this!
      }
      //      else rep.error("no lower stage to recombine to!",iim1);
      if (iip1) {
        // may or may not have higher stage.
        temp = Rad_Recomb_rate(temperature, iip1) * P[lv_elec]
               * P[lv_nh];  // rate [1/s]
        temp *= P[iip1->index];
        R[xi] += temp;  // recomb to ion i adds to fraction.
        // Count the energy loss of the electron when calculating R[i+1]
      }
    }
    // rep.printVec("\t\trr rate",R,nvl);

    if (ep.phot_ionisation) {
      if (ii[i]->i == H_1p) {  // only know Hydrogen photo-ionisation for now.
        // if (P[xi]>1.) cout <<"P[xi]="<<P[xi]<<"\n"; //1.0-1.e-12;
        R[lv_dtau] = P[lv_nh] * (1.0 - P[xi]) * phot_xsection(ii[i])
                     * path_length;  // this is optical depth through cell.
        if (R[lv_dtau] > 0.01)
          temp = photons_in * (1.0 - exp(-R[lv_dtau]))
                 / P[lv_nh];  // this is per sec.
        else
          temp = photons_in * phot_xsection(ii[i]) * path_length
                 * (1.0 - P[xi]);  // approx for 1-x<<1
        //  cout <<"\t\tphot rate: "<<temp<<"  and dtau_cell:
        //  "<<R[lv_dtau]<<" x="<<P[xi]<<"\n";
        // if (temp>1. && P[xi]>0.99) {
        //  cout <<"\t\tphot rate: "<<temp<<"  and dtau_cell:
        //  "<<R[lv_dtau]<<" x="<<P[xi]<<"\n"; R[xi] = 1.0;
        //}
        R[xi] += temp;
        R[lv_eint] += temp * 1.602e-12 * 2.4;  // this adds in 2.4eV of energy
                                               // per photo-ionisation.
        // have calculated the full rate for ionised hydrogen now, so
        // the rate for dtau is simple:
        R[lv_dtau] = exp(-R[lv_dtau]);
      }
    }
    // rep.printVec("\t\tpi rate",R,nvl);

    if (P[xi] >= 1.0 - MACHINEACCURACY) {
      // cout<<"P = "<<P[xi]<<" and R = "<<R[xi]<<"\n";
      R[xi] = min(0.0, R[xi]);
    }

    // add current ion production rate to electron production rate.
    R[lv_elec] += R[xi] * ii[i]->charge * ii[i]->e->numfrac;
    // rep.printVec("\t\txe rate",R,nvl);
  }

  if (ep.cooling) {
    temp = cool->CoolingRate(
        temperature, P[lv_hp], P[lv_nh], FUV_unattenuated_flux, FUV_extinction);
    R[lv_eint] -=
        temp / P[lv_nh];  // (because we need to multiply by nH below...
  }

  R[lv_eint] *=
      P[lv_nh];  // To convert to energy per unit volume per unit time.

  if (!ep.update_erg) {
    //    cout <<"not updating energy.\n";
    R[lv_eint] = 0.0;
  }
  // rep.printVec("\t\ttt rate",R,nvl);
  return 0;
}

// ##################################################################
// ##################################################################

int MPv0::C_rate(
    const int nv,     ///< number of variables we are expecting.
    const double *P,  ///< Current state vector.
    double *R         ///< Rate Vector to write to.
)
{
  return 1;
}

// ##################################################################
// ##################################################################

int MPv0::D_rate(
    const int nv,     ///< number of variables we are expecting.
    const double *P,  ///< Current state vector.
    double *R         ///< Rate Vector to write to.
)
{
  return 1;
}

// ##################################################################
// ##################################################################

double MPv0::phot_xsection(const struct ion_struct *ci)
{
  if (ci->i == H_1p)
    return 6.3e-18;  // in cm^2
  else {
    cout << "unknown photo-ionisation x-section\n";
    return 0.0;
  }
}

// ##################################################################
// ##################################################################

double MPv0::Coll_Ion_rate(
    double temperature,          ///< Precalculated Temperature.
    const struct ion_struct *ci  ///< current ion.
)
{
  // This uses fitting formulae from Voronov (1997) ADANDT, 65, 1.
  // (Atomic Data And Nuclear Data Tables)
  // rate returned in cm^3/s
  double A = 0.0, X = 0.0, K = 0.0;
  int PP = 0;
  if (ci->i == H_0) {
    if (temperature < 5.0e3) return 0.0;
    PP = 0;
    A  = 2.91e-8;
    X  = 0.232;
    K  = 0.39;
  }
  else if (ci->i == He0) {
    if (temperature < 5.0e3) return 0.0;
    PP = 0;
    A  = 1.75e-8;
    X  = 0.180;
    K  = 0.35;
  }
  else if (ci->i == He1p) {
    if (temperature < 1.0e4) return 0.0;
    PP = 1;
    A  = 2.05e-9;
    X  = 0.265;
    K  = 0.25;
  }

  // carbon
  else if (ci->i == C0) {
    if (temperature < 3.0e2) return 0.0;
    PP = 0;
    A  = 0.685e-7;
    X  = 0.193;
    K  = 0.25;
  }
  else if (ci->i == C1p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.186e-7;
    X  = 0.286;
    K  = 0.24;
  }
  else if (ci->i == C2p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.635e-8;
    X  = 0.427;
    K  = 0.21;
  }
  else if (ci->i == C3p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.150e-8;
    X  = 0.416;
    K  = 0.13;
  }
  else if (ci->i == C4p) {
    if (temperature < 5.0e4) return 0.0;
    PP = 1;
    A  = 0.299e-9;
    X  = 0.666;
    K  = 0.02;
  }
  else if (ci->i == C5p) {
    if (temperature < 5.0e4) return 0.0;
    PP = 1;
    A  = 0.123e-9;
    X  = 0.620;
    K  = 0.16;
  }

  // nitrogen
  else if (ci->i == N0) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.482e-7;
    X  = .0652;
    K  = 0.42;
  }
  else if (ci->i == N1p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.298e-7;
    X  = 0.310;
    K  = 0.30;
  }
  else if (ci->i == N2p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.810e-8;
    X  = 0.350;
    K  = 0.24;
  }
  else if (ci->i == N3p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.371e-8;
    X  = 0.549;
    K  = 0.18;
  }
  else if (ci->i == N4p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.151e-8;
    X  = .0167;
    K  = 0.74;
  }
  else if (ci->i == N5p) {
    if (temperature < 1.0e5) return 0.0;
    PP = 0;
    A  = 0.371e-9;
    X  = 0.546;
    K  = 0.29;
  }
  else if (ci->i == N6p) {
    if (temperature < 1.0e5) return 0.0;
    PP = 1;
    A  = 0.777e-10;
    X  = 0.624;
    K  = 0.16;
  }

  // oxygen
  else if (ci->i == O0) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.359e-7;
    X  = 0.073;
    K  = 0.34;
  }
  else if (ci->i == O1p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.139e-7;
    X  = 0.212;
    K  = 0.22;
  }
  else if (ci->i == O2p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.931e-8;
    X  = 0.270;
    K  = 0.27;
  }
  else if (ci->i == O3p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.102e-7;
    X  = 0.614;
    K  = 0.27;
  }
  else if (ci->i == O4p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 1;
    A  = 0.219e-8;
    X  = 0.630;
    K  = 0.17;
  }
  else if (ci->i == O5p) {
    if (temperature < 1.0e3) return 0.0;
    PP = 0;
    A  = 0.195e-8;
    X  = 0.360;
    K  = 0.54;
  }
  else if (ci->i == O6p) {
    if (temperature < 1.0e5) return 0.0;
    PP = 0;
    A  = 0.212e-9;
    X  = 0.396;
    K  = 0.35;
  }
  else if (ci->i == O7p) {
    if (temperature < 1.0e5) return 0.0;
    PP = 1;
    A  = 0.521e-10;
    X  = 0.629;
    K  = 0.16;
  }
  else
    rep.error("Bad ion in Coll_Ion_rate()", ii);

  //  if (ci->i == C3p) {
  //  cout <<"temp="<<ci->ion_pot/kB/temperature<<" and ci
  //  rate="<<A*(1.+PP*sqrt(temperature))*exp(K*log(temperature)
  //  -temperature)/(X+temperature)<<"  exp: "<<-temperature<<"\n";
  //}
  /*
  if    (ci->i ==C1p)
    {return 3.93e-11*exp(0.5*log(temperature) -2.83e5/temperature);}
  else if (ci->i ==C2p)
    {return 2.04e-11*exp(0.5*log(temperature) -5.556e5/temperature);}
  else if (ci->i ==C3p)
    {return 0.0;}
  else {
    temperature = ci->ion_pot/kB/temperature;
    return A*(1.+PP*sqrt(temperature))*exp(K*log(temperature)
  -temperature)/(X+temperature);
  }
  */
  temperature = ci->ion_pot / kB / temperature;
  return A * (1. + PP * sqrt(temperature))
         * exp(K * log(temperature) - temperature) / (X + temperature);
}

// ##################################################################
// ##################################################################

double MPv0::Rad_Recomb_rate(
    double temperature,          ///< Precalculated Temperature.
    const struct ion_struct *ci  ///< current ion.
)
{
#ifdef NEW_RATES
  double rate = 0.0;
  rate += rad_recomb(temperature, ci->i);
  rate += dielec_recomb(temperature, ci->i);
  return rate;
#endif  // NEW_RATES
#ifdef RAGA_RATES
  // rate is for recombination from current ion *ci to one stage lower.
  double r, a1, a2, a3, a4, a5, a6;
  if (ci->i == H_1p) {
    // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
    // For case B recombination coefficient.
    r = 3.41202e-10 * exp(-0.782991 * log(temperature));
  }
  else if (ci->i == He1p) {
    // recombination from He(1+) (from Fabio... Raga, deColle, et al., 2007,
    // A&A, 465, 879).
    // r  = 4.3e-13*exp(-0.672*log(temperature/1.e4));
    // r += 0.0019*exp(-1.5*log(temperature)
    // -4.7e5/temperature)*(1.0+0.3*exp(-9.4e4/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1 = 9.356e-10;
    a2 = 0.7892;
    a3 = 4.266e-2;
    a4 = 4.677e6;
    r  = a1
        / (sqrt(temperature / a3)
           * exp((1. - a2) * log(1. + sqrt(temperature / a3))
                 + (1. + a2) * log(1. + sqrt(temperature / a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature, ci->i);
  }
  else if (ci->i == He2p) {
    // recombination from He(2+) (from Fabio... Raga, deColle, et al., 2007,
    // A&A, 465, 879).
    // r  = 2.21e-9*exp(-0.79*log(temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1 = 1.891e-10;
    a2 = 0.7524;
    a3 = 9.370;
    a4 = 2.774e6;
    r  = a1
        / (sqrt(temperature / a3)
           * exp((1. - a2) * log(1. + sqrt(temperature / a3))
                 + (1. + a2) * log(1. + sqrt(temperature / a4))));
  }
  // recombination for Carbon (from Raga, deColle, et al., 2007, A&A, 465,
  // 879).
  else if (ci->i == C1p) {
    a1 = 4.7e-13;
    a2 = 0.624;
    a3 = 6.9e-4;
    a4 = 1.1e5;
    a5 = 3.0;
    a6 = 4.9e4;
    // r = 1.0e-22; //a1*exp(a2*log(temperature) +a3/temperature);
    r = a1 * exp(-a2 * log(temperature / 1.e4))
        + a3 * exp(-1.5 * log(temperature) - a4 / temperature)
              * (1. + a5 * exp(-a6 / temperature));
  }
  else if (ci->i == C2p) {
    a1 = 2.3e-12;
    a2 = 0.645;
    a3 = 0.007;
    a4 = 1.5e5;
    a5 = 0.5;
    a6 = 2.3e5;
    //    r = a1*exp(a2*log(temperature) +a3/temperature);
    r = a1 * exp(-a2 * log(temperature / 1.e4))
        + a3 * exp(-1.5 * log(temperature) - a4 / temperature)
              * (1. + a5 * exp(-a6 / temperature));
  }
  else if (ci->i == C3p) {
    a1 = 3.2e-12;
    a2 = 0.770;
    a3 = 3.8e-3;
    a4 = 9.1e4;
    a5 = 2.0;
    a6 = 3.7e5;
    //    r = a1*exp(a2*log(temperature) +a3/temperature);
    r = a1 * exp(-a2 * log(temperature / 1.e4))
        + a3 * exp(-1.5 * log(temperature) - a4 / temperature)
              * (1. + a5 * exp(-a6 / temperature));
  }
  else if (ci->i == C4p) {
    // a1=; a2=; a3=; a4=; a5=; a6=;
    // r = a1*exp(a2*log(temperature) +a3/temperature);
    // r = a1*exp(-a2*log(temperature/1.e4))
    //  +a3*exp(-1.5*log(temperature)
    //  -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1 = 8.540e-11;
    a2 = 0.5247;
    a3 = 5.014e2;
    a4 = 1.479e7;  // a5=; a6=;
    r  = a1
        / (sqrt(temperature / a3)
           * exp((1. - a2) * log(1. + sqrt(temperature / a3))
                 + (1. + a2) * log(1. + sqrt(temperature / a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature, ci->i);
  }
  else if (ci->i == C5p) {
    // a1=; a2=; a3=; a4=; a5=; a6=;
    // r = a1*exp(a2*log(temperature) +a3/temperature);
    // r = a1*exp(-a2*log(temperature/1.e4))
    //  +a3*exp(-1.5*log(temperature)
    //  -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1 = 2.765e-10;
    a2 = 0.6858;
    a3 = 1.535e2;
    a4 = 2.556e7;  // a5=; a6=;
    r  = a1
        / (sqrt(temperature / a3)
           * exp((1. - a2) * log(1. + sqrt(temperature / a3))
                 + (1. + a2) * log(1. + sqrt(temperature / a4))));
    // ******* Mazzotta et al (1998) dielectronic rate.
    r += dielec_recomb(temperature, ci->i);
  }
  else if (ci->i == C6p) {
    // a1=; a2=; a3=; a4=; a5=; a6=;
    // r = a1*exp(a2*log(temperature) +a3/temperature);
    // r = a1*exp(-a2*log(temperature/1.e4))
    //  +a3*exp(-1.5*log(temperature)
    //  -a4/temperature)*(1.+a5*exp(-a6/temperature));
    // ******* Verner & Ferland (1996) ApJS, 103, 467. Four parameter fit.
    a1 = 6.556e-10;
    a2 = 0.7567;
    a3 = 6.523e1;
    a4 = 2.446e7;  // a5=; a6=;
    r  = a1
        / (sqrt(temperature / a3)
           * exp((1. - a2) * log(1. + sqrt(temperature / a3))
                 + (1. + a2) * log(1. + sqrt(temperature / a4))));
  }
  else
    rep.error("unknown ion in Rad_Recomb_rate()", ii);
  return r;
#endif  // RAGA_RATES
}

// ##################################################################
// ##################################################################

double MPv0::rad_recomb(double T, enum species i)
{
  // rates are from Verner & Ferland (1996) ApJS, 103, 467 where possible,
  // and from Pequignot, Petitjean, & Boisson, (1991) A&A, 251, 680 otherwise,
  // refitted with VF96 formula.
  //
  // Fortran subroutine with list of rates at:
  // http://hea-www.harvard.edu/~mazzotta/sub_mazz2/rrfit.f
  // ftp://gradj.pa.uky.edu//dima//rec//rrfit.f via
  // http://www.pa.uky.edu/~verner/fortran.html
  double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, r = 0.0;
  if (i == H_1p) {
    // For hydrogen, we calculate the VF96 case A recomb. coeff., and
    // then use the minimum of this and the Storey & Hummer (1995) coeff.
    // a1=7.982e-11; a2=0.7480; a3=3.148; a4=7.036e5;
    // r =
    // a1/(sqrt(T/a3)*exp((1.-a2)*log(1.+sqrt(T/a3))+(1.+a2)*log(1.+sqrt(T/a4))));
    //    if (T>1.e3)
    // cout <<"T="<<T<<"\t r(VF)="<<r<<"\t
    // r(SH)="<<3.41202e-10*exp(-0.782991*log(T))<<"\n"; r =
    // min(r,3.41202e-10*exp(-0.782991*log(T)));

    // This is fit to data in Storey & Hummer (1995), MNRAS, 272, 41.
    // For case B recombination coefficient.
    r = 3.41202e-10 * exp(-0.782991 * log(T));
    // r = 2.59e-13; // basic approx.
  }
  else if (i == He1p) {
    a1 = 9.356e-10;
    a2 = 0.7892;
    a3 = 4.266e-2;
    a4 = 4.677e6;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == He2p) {
    a1 = 1.891e-10;
    a2 = 0.7524;
    a3 = 9.370;
    a4 = 2.774e6;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C1p) {
    // 7.651E-09,0.8027,1.193E-03,9.334E+12
    a1 = 7.651E-09;
    a2 = 0.8027;
    a3 = 1.193E-03;
    a4 = 9.334E+12;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C2p) {
    // 8.577E-10,0.7837,7.286E-01,1.140E+07
    a1 = 8.577E-10;
    a2 = 0.7837;
    a3 = 7.286E-01;
    a4 = 1.140E+07;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C3p) {
    // 2.020E-09,0.7798,6.690E-01,2.425E+06
    a1 = 2.020E-09;
    a2 = 0.7798;
    a3 = 6.690E-01;
    a4 = 2.425E+06;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C4p) {
    // 8.540e-11,0.5247,5.014e+02,1.479e+07
    a1 = 8.540e-11;
    a2 = 0.5247;
    a3 = 5.014e2;
    a4 = 1.479e7;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C5p) {
    // 2.765e-10,0.6858,1.535e+02,2.556e+07
    a1 = 2.765e-10;
    a2 = 0.6858;
    a3 = 1.535e2;
    a4 = 2.556e7;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else if (i == C6p) {
    // 6.556e-10,0.7567,6.523e+01,2.446e+07
    a1 = 6.556e-10;
    a2 = 0.7567;
    a3 = 6.523e1;
    a4 = 2.446e7;
    r  = a1
        / (sqrt(T / a3)
           * exp((1. - a2) * log(1. + sqrt(T / a3))
                 + (1. + a2) * log(1. + sqrt(T / a4))));
  }
  else
    rep.error("unknown ion in Rad_Recomb_rate()", ii);
  return r;
}

// ##################################################################
// ##################################################################

double MPv0::dielec_recomb(double T, enum species i)
{
  // i is the ion we are recombining *from*
  // rates are from Mazzotta et al. (1998)
  double c[4];
  c[0] = c[1] = c[2] = c[3] = 0.0;
  double E[4];
  E[0] = E[1] = E[2] = E[3] = 0.0;
  if (i == H_1p) {
    return 0.0;
  }
  else if (i == He1p) {
    c[0] = 1.12e-9;
    E[0] = 39.70;
  }
  else if (i == He2p) {
    return 0.0;
  }
  else if (i == C1p) {
    c[0] = 1.0422e-09;
    c[1] = 5.8484e-10;
    c[2] = 5.6306e-11;
    E[0] = 12.57;
    E[1] = 162.90;
    E[2] = 6.30;
  }
  else if (i == C2p) {
    c[0] = 4.6178e-12;
    c[1] = 2.8234e-09;
    E[0] = 0.49;
    E[1] = 11.78;
  }
  else if (i == C3p) {
    c[0] = 5.3858e-10;
    c[1] = 1.5056e-09;
    c[2] = 6.0332e-10;
    E[0] = 185.96;
    E[1] = 17.99;
    E[2] = 2.41;
  }
  else if (i == C4p) {
    c[0] = 1.4008e-08;
    E[0] = 287.34;
  }
  else if (i == C5p) {
    c[0] = 3.3558e-08;
    E[0] = 356.46;
  }
  else if (i == C6p) {
    return 0.0;
  }
  // else if (i=="") {
  //}
  else {
    cout << "unknown ion: " << i << "\n";
    return -1.0;
  }

  T /= 1.16e4;  // convert from K to eV.
  //  if (T<=1.e-100) {cout <<"tiny temperature\n";return -1.0;}

  double rate = 0.0;
  for (int i = 0; i < 4; i++)
    rate += c[i] * exp(-E[i] / T);
  rate *= exp(-1.5 * log(T));
  return rate;
}

// ##################################################################
// ##################################################################

double MPv0::charge_exchange_rate(
    double T,                     ///< Precalculated Temperature.
    const struct ion_struct *ci,  ///< current ion.
    const int sw  ///< 0 for recombination from ion, 1 for ionisation of ion.
)
{
  // Charge Exchange rates are from Kingdon & Ferland, 1996, 106, 205.
  double r = 0.0;
  return r;
}

// ##################################################################
// ##################################################################

int MPv0::Init_ionfractions(
    pion_flt *p_prim,  ///< Primitive vector to be updated.
    const double gam,  ///< eos gamma.
    const double temp  ///< optional gas temperature to end up at. (negative
                       ///< means use pressure)
)
{
  gamma    = gam;
  double T = temp, p_local[nvl];
  int err  = convert_prim2local(p_prim, p_local, gamma);

  //  rep.printVec("p_local start",p_local,nvl);
  if (T > 0)
    err += Set_Eint(p_local, T);
  else
    T = Get_Temp(p_local);
  //  cout <<"Temp = "<<T<<"\n";
  if (nions <= 0) {
    cout << "no ions to init.\n";
    err += convert_local2prim(p_local, p_prim, p_prim, gamma);
    return err;
  }

  // cout <<"start Temp = "<<T<<"\t n_h="<<p_local[lv_nh]<<"\n";
  // This sets the first tracked ion of each element to have abundance=1
  // and all the others to have abundance 0.  This is the best way to
  // initialise, and we can run a long time update to get to equilibrium.
  for (int i = 0; i < nels; i++) {
    for (int j = 0; j < ee[i]->nspecies; j++) {
      if (j == 0)
        p_local[ee[i]->ion_indices[j]] = 1.0;
      else
        p_local[ee[i]->ion_indices[j]] = 0.0;
    }
  }
  p_local[lv_elec] = 0.0;
  for (int i = 0; i < nions; i++) {
    p_local[lv_elec] +=
        p_local[ii[i]->index] * ii[i]->charge * ii[i]->e->numfrac;
  }
  p_local[lv_elec] = max(min_elecf, p_local[lv_elec]);
  //  cout <<"end Temp = "<<T<<"\t n_h="<<p_local[lv_nh]<<"\n";
  err += convert_local2prim(p_local, p_prim, p_prim, gamma);
  return err;

  /*
  // we have ions, and safe to assume hydrogen is one of them, so use H
  // to set electron fraction using Saha Eqn.
  double saha_prefactor = 2.41e15;
  ion_struct *H0    = &ion_props["H0"];
  //  ion_struct *Hplus = &ion_props["H1+"];
  int hp = lvar["H1+"];
  double B = saha_prefactor*exp(1.5*log(T) -H0->ion_pot/kB/T) /p_local[lv_nh];
  if (B>1.e12) p_local[hp] = p_local[lv_elec] = 1.0;
  else if (B<1.e-50) p_local[hp] = min_elecf;
  else p_local[hp] = -B/2.0 *(1.0-sqrt(1.0+4.0/B));
  p_local[lv_elec] = max(min_elecf, p_local[hp]);
  if (isnan(p_local[hp])) {
    cout <<"nans!!!!!!!!! B="<<B<<" and p_local[H1+]="<<p_local[hp]<<"\n";
    cout <<"T="<<T<<" and p_local[lv_nh]="<<p_local[lv_nh]<<"\n";
  }

  for (int i=0; i<nels; i++) {
    if (ee[i]->el !="H") {
      int ni = ee[i]->ions.size();
      double ifracs[ni]; for (int e=0;e<ni;e++) ifracs[e] = 0.0;
      for (int e=1;e<ni;e++) {
        ion_struct *iii = &ion_props[ee[i]->ions[e]];
        //cout <<"ion: "<<iii->ion<<"  iii->im1: "<<iii->im1<<"\n";
        ion_struct *im1 = &ion_props[iii->im1];
        ifracs[e] = saha_prefactor*2.0*exp(1.5*log(T) -im1->ion_pot/kB/T);
        ifracs[e] *= (iii->g_stat/im1->g_stat)
  /p_local[lv_elec]/p_local[lv_nh]; if (e>1) for (int eee=e; eee<ni; eee++)
  ifracs[eee] *= ifracs[eee-1];
      }
      double nfrac = 1.0;
      for (int e=1;e<ni;e++) nfrac += ifracs[e];
      nfrac = 1.0/nfrac;
      for (int e=1;e<ni;e++) ifracs[e] *= nfrac;
      //     cout <<"ee = "<<ee[i]->el<<" and nfrac="<<nfrac<<" with ifracs:
  ";
      //      for (int e=1;e<ni;e++) cout <<ifracs[e]<<", "; cout <<"\n";
      // now set local state vector to have calculated ion fractions
      for (int e=1;e<ni;e++) p_local[ee[i]->ion_indices[e-1]] = ifracs[e];
    }
  }
  p_local[lv_elec] = 0.0;
  for (int i=0;i<nions; i++) {
    p_local[lv_elec] += p_local[ii[i]->index] *ii[i]->charge
  *ii[i]->e->numfrac;
  }
  p_local[lv_elec] = max(min_elecf,p_local[lv_elec]);
  //  cout <<"elec fraction = "<<p_local[lv_elec]<<"\n";
  //  err += Set_Eint(p_local,T);
  //  rep.printVec("p_local finish",p_local,nvl);
  */
}

// ##################################################################
// ##################################################################

/// \page userguide
/// \section atomicdata Atomic Data
/// Abundance by number (N) From Kaye & Laby:
/// http://www.kayelaby.npl.co.uk/chemistry/3_1/3_1_3.html \n Abundance A(EI)
/// from Lodders, 2003, ApJ, 591, 1220-1247, table 2, which refers to solar
/// system abundances.\n Masses from Wikipedia, in AMU.  The abundances A(EI)
/// are used in the code.
///
/// <TABLE>
/// <tr> <td>element</td>	  <td>mass</td>	  <td>abundance(N)</td>
/// <td>abundance(A(EI)=log10(n(EI)/n(H))+12)</td>  </tr> <tr> <td>H</td>
/// <td>1.00</td>   <td>2.8e10</td>      <td>12</td>    </tr> <tr> <td>He</td>
/// <td>4.00</td>   <td>2.7e9</td>       <td>10.98</td> </tr> <tr> <td>Li</td>
/// <td>6.94</td>   <td>57.0</td>	       <td>3.35</td>  </tr> <tr>
/// <td>Be</td>	  <td>9.01</td>   <td>0.7</td>	       <td>1.48</td>  </tr> <tr>
/// <td>B</td>	  <td>10.81</td>  <td>21.0</td>	       <td>2.85</td>  </tr> <tr>
/// <td>C</td>	  <td>12.01</td>  <td>1.0e7</td>       <td>8.46</td>  </tr> <tr>
/// <td>N</td>	  <td>14.01</td>  <td>3.1e6</td>       <td>7.90</td>  </tr> <tr>
/// <td>O</td>	  <td>16.00</td>  <td>2.4e7</td>       <td>8.76</td>  </tr> <tr>
/// <td>F</td>	  <td>19.00</td>  <td>8.5e2</td>       <td>4.53</td>  </tr> <tr>
/// <td>Ne</td>	  <td>20.18</td>  <td>3.0e6</td>       <td>7.95</td>  </tr> <tr>
/// <td>Na</td>	  <td>22.99</td>  <td>5.7e4</td>       <td>6.37</td>  </tr> <tr>
/// <td>Mg</td>	  <td>24.31</td>  <td>1.1e6</td>       <td>7.62</td>  </tr> <tr>
/// <td>Al</td>	  <td>26.98</td>  <td>8.5e4</td>       <td>6.54</td>  </tr> <tr>
/// <td>Si</td>	  <td>28.09</td>  <td>1.0e6</td>       <td>7.61</td>  </tr> <tr>
/// <td>P</td>	  <td>30.97</td>  <td>1.0e4</td>       <td>5.54</td>  </tr>
/// </table>

// ##################################################################
// ##################################################################

void MPv0::set_atomic_data()
{
  // list of all possible ions to treat in the code.
  ion_list.push_back("H1+");
  ion_list.push_back("He1+");
  ion_list.push_back("He2+");

  // list of all possible elements to treat in the code.
  el_list.push_back("H");
  el_list.push_back("He");
  el_list.push_back("Li");
  //  els.push_back = "";

  double eV_per_Erg = 1.602e-12;      // number of eV in 1erg
  double amu_to_mp  = 1.661 / 1.673;  // number of amu in one proton mass.
  struct ion_struct i;
  struct element_struct e;

  /** \section numfrac Number Fractions
   * Number fractions are from Lodders,K., 2003, ApJ, 591, 1220-1247
   *  \section Masses
   * Masses are from Wikipedia, in Atomic Mass Units (amu), converted to units
   * of proton mass.
   * */
  e.ion_indices = 0;

  e.el = "H";
  e.ions.clear();
  e.ions.push_back("H0");
  e.ions.push_back("H1+");
  e.nspecies    = 1;
  e.mass        = 1.01 * amu_to_mp;
  e.numfrac     = 1.00;
  el_props["H"] = e;

  e.el = "He";
  e.ions.clear();
  e.ions.push_back("He0");
  e.ions.push_back("He1+");
  e.ions.push_back("He2+");
  e.nspecies     = 2;
  e.mass         = 4.00 * amu_to_mp;
  e.numfrac      = pow(10.0, -1.02);  // 0.0955
  el_props["He"] = e;

  e.el = "Li";
  e.ions.clear();
  e.ions.push_back("Li0");
  e.ions.push_back("Li1+");
  e.ions.push_back("Li2+");
  e.ions.push_back("Li3+");
  e.nspecies     = 3;
  e.mass         = 6.94 * amu_to_mp;
  e.numfrac      = pow(10.0, 3.35 - 12.0);
  el_props["Li"] = e;

  // C 	12.01 	1.0e7 	8.46
  e.el = "C";
  e.ions.clear();
  e.ions.push_back("C0");
  e.ions.push_back("C1+");
  e.ions.push_back("C2+");
  e.ions.push_back("C3+");
  e.ions.push_back("C4+");
  e.ions.push_back("C5+");
  e.ions.push_back("C6+");
  e.nspecies    = 6;
  e.mass        = 12.01 * amu_to_mp;
  e.numfrac     = pow(10.0, 8.46 - 12.0);  // 2.884e-4
  el_props["C"] = e;

  // N 	14.01 	3.1e6 	7.90
  e.el = "N";
  e.ions.clear();
  e.ions.push_back("N0");
  e.ions.push_back("N1+");
  e.ions.push_back("N2+");
  e.ions.push_back("N3+");
  e.ions.push_back("N4+");
  e.ions.push_back("N5+");
  e.ions.push_back("N6+");
  e.ions.push_back("N7+");
  e.nspecies    = 7;
  e.mass        = 14.01 * amu_to_mp;
  e.numfrac     = pow(10.0, 7.90 - 12.0);
  el_props["N"] = e;

  // O 	16.00 	2.4e7 	8.76
  e.el = "O";
  e.ions.clear();
  e.ions.push_back("O0");
  e.ions.push_back("O1+");
  e.ions.push_back("O2+");
  e.ions.push_back("O3+");
  e.ions.push_back("O4+");
  e.ions.push_back("O5+");
  e.ions.push_back("O6+");
  e.ions.push_back("O7+");
  e.ions.push_back("O8+");
  e.nspecies    = 8;
  e.mass        = 16.00 * amu_to_mp;
  e.numfrac     = pow(10.0, 8.76 - 12.0);
  el_props["O"] = e;

  // Si 	28.09 	1.0e6 	7.61
  e.el = "Si";
  e.ions.clear();
  e.ions.push_back("Si0");
  e.ions.push_back("Si1+");
  e.ions.push_back("Si2+");
  e.ions.push_back("Si3+");
  e.ions.push_back("Si4+");
  e.ions.push_back("Si5+");
  e.ions.push_back("Si6+");
  e.ions.push_back("Si7+");
  e.ions.push_back("Si8+");
  //  e.ions.push_back("Si14+"); // can include up to 14 ionisation stages!
  e.nspecies     = 8;
  e.mass         = 28.09 * amu_to_mp;
  e.numfrac      = pow(10.0, 7.61 - 12.0);
  el_props["Si"] = e;

  // S 32.07 7.26pm0.04
  e.el = "S";
  e.ions.clear();
  e.ions.push_back("S0");
  e.ions.push_back("S1+");
  e.ions.push_back("S2+");
  //  e.ions.push_back("S1+");
  //  e.ions.push_back("+"); // can include up to  ionisation stages!
  e.nspecies    = 2;
  e.mass        = 32.07 * amu_to_mp;
  e.numfrac     = pow(10.0, 7.26 - 12.0);
  el_props["S"] = e;
  /*
  //
  e.el = "";
  e.ions.clear();
  e.ions.push_back("0" );
  e.ions.push_back("1+");
  //  e.ions.push_back("+"); // can include up to  ionisation stages!
  e.nspecies = ;
  e.mass = *amu_to_mp;
  e.numfrac = pow(10.0,-12.0);
  el_props["" ] = e;
  */
  // etc. for more elements.

  i.ion           = "H0";
  i.ip1           = "H1+";
  i.im1           = "";
  i.el            = "H";
  i.charge        = 0;
  i.ion_pot       = eV_per_Erg * 13.59844;
  i.i             = H_0;
  i.index         = -1;
  i.pv_index      = -1;
  i.e             = &(el_props["H"]);
  i.g_stat        = 2;
  i.ch_ex         = false;
  ion_props["H0"] = i;

  i.ion = "H1+";
  i.ip1 = "";
  i.im1 = "H0";
  //  i.ip1 = NONE;
  //  i.im1 = H_0;
  i.el             = "H";
  i.charge         = 1;
  i.ion_pot        = -1.0e99;
  i.i              = H_1p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["H"]);
  i.g_stat         = 1;
  i.ch_ex          = false;
  ion_props["H1+"] = i;

  i.ion = "He0";
  //  i.ip1 = He1p;
  //  i.im1 = NONE;
  i.ip1            = "He1+";
  i.im1            = "";
  i.el             = "He";
  i.charge         = 0;
  i.ion_pot        = eV_per_Erg * 24.58741;
  i.i              = He0;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["He"]);
  i.g_stat         = 1;
  i.ch_ex          = false;
  ion_props["He0"] = i;

  i.ion = "He1+";
  //  i.ip1 = He2p;
  //  i.im1 = He0;
  i.ip1             = "He2+";
  i.im1             = "He0";
  i.el              = "He";
  i.charge          = 1;
  i.ion_pot         = eV_per_Erg * 54.41778;
  i.i               = He1p;
  i.index           = -1;
  i.pv_index        = -1;
  i.e               = &(el_props["He"]);
  i.g_stat          = 2;
  i.ch_ex           = false;
  ion_props["He1+"] = i;

  i.ion = "He2+";
  //  i.ip1 = NONE;
  //  i.im1 = He1p;
  i.ip1             = "";
  i.im1             = "He1+";
  i.el              = "He";
  i.charge          = 2;
  i.ion_pot         = -1.e99;
  i.i               = He2p;
  i.index           = -1;
  i.pv_index        = -1;
  i.e               = &(el_props["He"]);
  i.g_stat          = 1;
  i.ch_ex           = false;
  ion_props["He2+"] = i;

  // carbon: C, C1+, C2+, C3+, C4+, C5+, C6+
  i.ion           = "C0";
  i.ip1           = "C1+";
  i.im1           = "";
  i.el            = "C";
  i.charge        = 0;
  i.ion_pot       = eV_per_Erg * 11.3;
  i.i             = C0;
  i.index         = -1;
  i.pv_index      = -1;
  i.e             = &(el_props["C"]);
  i.g_stat        = 1;
  ion_props["C0"] = i;
  i.ch_ex         = true;

  i.ion            = "C1+";
  i.ip1            = "C2+";
  i.im1            = "C0";
  i.el             = "C";
  i.charge         = 1;
  i.ion_pot        = eV_per_Erg * 24.4;
  i.i              = C1p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C1+"] = i;
  i.ch_ex          = true;

  i.ion            = "C2+";
  i.ip1            = "C3+";
  i.im1            = "C1+";
  i.el             = "C";
  i.charge         = 2;
  i.ion_pot        = eV_per_Erg * 47.9;
  i.i              = C2p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C2+"] = i;
  i.ch_ex          = true;

  i.ion            = "C3+";
  i.ip1            = "C4+";
  i.im1            = "C2+";
  i.el             = "C";
  i.charge         = 3;
  i.ion_pot        = eV_per_Erg * 64.5;
  i.i              = C3p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C3+"] = i;
  i.ch_ex          = true;

  i.ion            = "C4+";
  i.ip1            = "C5+";
  i.im1            = "C3+";
  i.el             = "C";
  i.charge         = 4;
  i.ion_pot        = eV_per_Erg * 392.1;
  i.i              = C4p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C4+"] = i;
  i.ch_ex          = true;

  i.ion            = "C5+";
  i.ip1            = "C6+";
  i.im1            = "C4+";
  i.el             = "C";
  i.charge         = 5;
  i.ion_pot        = eV_per_Erg * 490.0;
  i.i              = C5p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C5+"] = i;
  i.ch_ex          = false;

  i.ion            = "C6+";
  i.ip1            = "";
  i.im1            = "C5+";
  i.el             = "C";
  i.charge         = 6;
  i.ion_pot        = -1.e99;
  i.i              = C6p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["C"]);
  i.g_stat         = 1;
  ion_props["C6+"] = i;
  i.ch_ex          = false;

  // nitrogen: N, N1+, N2+, N3+, N4+, N5+, N6+, N7+
  i.ion           = "N0";
  i.ip1           = "N1+";
  i.im1           = "";
  i.el            = "N";
  i.charge        = 0;
  i.ion_pot       = eV_per_Erg * 14.5;
  i.i             = N0;
  i.index         = -1;
  i.pv_index      = -1;
  i.e             = &(el_props["N"]);
  i.g_stat        = 1;
  ion_props["N0"] = i;
  i.ch_ex         = true;

  i.ion            = "N1+";
  i.ip1            = "N2+";
  i.im1            = "N0";
  i.el             = "N";
  i.charge         = 1;
  i.ion_pot        = eV_per_Erg * 29.6;
  i.i              = N1p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N1+"] = i;
  i.ch_ex          = true;

  i.ion            = "N2+";
  i.ip1            = "N3+";
  i.im1            = "N1+";
  i.el             = "N";
  i.charge         = 2;
  i.ion_pot        = eV_per_Erg * 47.5;
  i.i              = N2p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N2+"] = i;
  i.ch_ex          = false;

  i.ion            = "N3+";
  i.ip1            = "N4+";
  i.im1            = "N2+";
  i.el             = "N";
  i.charge         = 3;
  i.ion_pot        = eV_per_Erg * 77.5;
  i.i              = N3p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N3+"] = i;
  i.ch_ex          = false;

  i.ion            = "N4+";
  i.ip1            = "N5+";
  i.im1            = "N3+";
  i.el             = "N";
  i.charge         = 4;
  i.ion_pot        = eV_per_Erg * 97.9;
  i.i              = N4p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N4+"] = i;
  i.ch_ex          = false;

  i.ion            = "N5+";
  i.ip1            = "N6+";
  i.im1            = "N+4";
  i.el             = "N";
  i.charge         = 5;
  i.ion_pot        = eV_per_Erg * 552.1;
  i.i              = N5p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N5+"] = i;
  i.ch_ex          = false;

  i.ion            = "N6+";
  i.ip1            = "N7+";
  i.im1            = "N5+";
  i.el             = "N";
  i.charge         = 6;
  i.ion_pot        = eV_per_Erg * 667.0;
  i.i              = N6p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N6+"] = i;
  i.ch_ex          = false;

  i.ion            = "N7+";
  i.ip1            = "";
  i.im1            = "N6+";
  i.el             = "N";
  i.charge         = 7;
  i.ion_pot        = -1.e99;
  i.i              = N7p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["N"]);
  i.g_stat         = 1;
  ion_props["N7+"] = i;
  i.ch_ex          = false;

  // oxygen: O, O1+, O2+, O3+, O4+, O5+, O6+, O7+, O8+,
  i.ion           = "O0";
  i.ip1           = "O1+";
  i.im1           = "";
  i.el            = "O";
  i.charge        = 0;
  i.ion_pot       = eV_per_Erg * 13.6;
  i.i             = O0;
  i.index         = -1;
  i.pv_index      = -1;
  i.e             = &(el_props["O"]);
  i.g_stat        = 1;
  ion_props["O0"] = i;
  i.ch_ex         = true;

  i.ion            = "O1+";
  i.ip1            = "O2+";
  i.im1            = "O0";
  i.el             = "O";
  i.charge         = 1;
  i.ion_pot        = eV_per_Erg * 35.1;
  i.i              = O1p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O1+"] = i;
  i.ch_ex          = true;

  i.ion            = "O2+";
  i.ip1            = "O3+";
  i.im1            = "O1+";
  i.el             = "O";
  i.charge         = 2;
  i.ion_pot        = eV_per_Erg * 54.9;
  i.i              = O2p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O2+"] = i;
  i.ch_ex          = false;

  i.ion            = "O3+";
  i.ip1            = "O4+";
  i.im1            = "O2+";
  i.el             = "O";
  i.charge         = 3;
  i.ion_pot        = eV_per_Erg * 77.4;
  i.i              = O3p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O3+"] = i;
  i.ch_ex          = false;

  i.ion            = "O4+";
  i.ip1            = "O5+";
  i.im1            = "O3+";
  i.el             = "O";
  i.charge         = 4;
  i.ion_pot        = eV_per_Erg * 113.9;
  i.i              = O4p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O4+"] = i;
  i.ch_ex          = false;

  i.ion            = "O5+";
  i.ip1            = "O6+";
  i.im1            = "O4+";
  i.el             = "O";
  i.charge         = 5;
  i.ion_pot        = eV_per_Erg * 138.1;
  i.i              = O5p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O5+"] = i;
  i.ch_ex          = false;

  i.ion            = "O6+";
  i.ip1            = "O7+";
  i.im1            = "O5+";
  i.el             = "O";
  i.charge         = 6;
  i.ion_pot        = eV_per_Erg * 739.3;
  i.i              = O6p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O6+"] = i;
  i.ch_ex          = false;

  i.ion            = "O7+";
  i.ip1            = "O8+";
  i.im1            = "O6+";
  i.el             = "O";
  i.charge         = 7;
  i.ion_pot        = eV_per_Erg * 871.4;
  i.i              = O7p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O7+"] = i;
  i.ch_ex          = false;

  i.ion            = "O8+";
  i.ip1            = "";
  i.im1            = "O7+";
  i.el             = "O";
  i.charge         = 8;
  i.ion_pot        = -1.e99;
  i.i              = O8p;
  i.index          = -1;
  i.pv_index       = -1;
  i.e              = &(el_props["O"]);
  i.g_stat         = 1;
  ion_props["O8+"] = i;
  i.ch_ex          = false;

  /*
    //ion
  i.ion = ""; i.ip1 = ""; i.im1 = ""; i.el  = "";
  i.charge = ; i.ion_pot = eV_per_Erg *;
  i.i = ; i.index = -1; i.pv_index = -1;
  i.e = &(el_props[""]); i.g_stat = 1; ion_props["" ] = i;
  */
  // etc. for other ions.

  return;
}

// ##################################################################
// ##################################################################

///
/// This returns the minimum timescale of the times flagged in the
/// arguments.  Time is returned in seconds.
///
double MPv0::timescales(
    const pion_flt *p_in,  ///< Current cell primitive vector.
    const double gam,      ///< EOS gamma.
    const bool f_cool,     ///< set to true if including cooling time.
    const bool f_recomb,   ///< set to true if including recombination time.
    const bool f_photoion  ///< set to true if including photo-ionsation time.
)
{
  //
  // This only works for "only cooling" so far.
  //
  if (nions == 0 && ep.cooling) {
    //
    // First convert to MP variables and get temperature.
    //
    //  cout <<"gamma: "<<gamma<<"\n";
    //  cout <<"only cooling: nvl = "<<nvl<<"\n";
    if (nvl != 2) rep.error("nvl wrong", nvl);
    gamma = gam;
    // set up local state vector.
    double P[nvl];
    P[lv_nh]   = p_in[RO] / m_p / 1.0;  // mu=1.0
    P[lv_eint] = p_in[PG] / (gam - 1.);
    double T   = (gam - 1.) * P[lv_eint] / kB / 2.0 / P[lv_nh];  // Temperature.

    //
    // Now get the timescale for each requested process.
    // Note this is does not include radiative transfer for now, so the
    // fluxes are set to zero, and DON'T REQUEST THE PHOTO-IONISATION TIME.
    // That is a wish-list feature for the future.
    //
    double mintime = 1.0e99;
    double rate    = 0.0;

    //
    // we skip the cooling time if the temperature is already very low;
    // a reasonable temperature could be 10K; there is unlikely to be
    // significant cooling below this temperature, and even if there is
    // it won't be by a large fraction.
    //
    if ((ep.cooling) && (f_cool) && (T >= 10.0)) {
      //
      // This gets the cooling rate per unit volume, first from the
      // cooling function, and then adding on the recombination cooling.
      //
      rate = cool->CoolingRate(T, 1.0, P[lv_nh], 0.0, 0.0);
      //
      // Cooling time is then the energy per unit volume divided by the
      // rate. Bear in mind we could have net heating, so take fabs(rate)
      //
      mintime = min(mintime, P[lv_eint] / fabs(rate));
    }
    return mintime;
  }
  else {
    rep.error("Can't limit timestep in Microphysics unless only cooling", 99);
  }
  return VERY_LARGE_VALUE;
}

// ##################################################################
// ##################################################################

#endif  // LEGACY_CODE
