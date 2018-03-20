///
/// \file mpv8_TwoTempIso.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running calculations with a simple 
/// two-temperature isothermal equation of state, where T is T_low
/// when gas is neutral, and T_high when ionised, and linearly
/// interpolated for partial ionisation.
///
/// Even if H is molecular, we assume that there is a PDR ahead of
/// the ionisation front, so the ionising photons always hit atomic
/// hydrogen.
///
///
/// - THE CODE IS NOT TESTED AGAINST STANDARD SOLUTIONS.
/// - THE CODE IS NOT TESTED IN MULTI-D, SO THE INTERPOLATION MAY BE
///    SKETCHY.
/// - THE CODE IS VERY UNSTABLE AT THE I-FRONT, DRIVING A FORCED
///    OSCILLATION INTO THE HII REGION INTERIOR!
///
/// Modifications:
/// - getting it written: mods up until 2013.08.13

#include "microphysics/mpv8_TwoTempIso.h"



using namespace std;



// ##################################################################
// ##################################################################


mpv8_TwoTempIso::mpv8_TwoTempIso(
          const int nv,              ///< Total number of variables in state vector
          const int ntracer,         ///< Number of tracer variables in state vector.
          const std::string &trtype,  ///< List of what the tracer variables mean.
          struct which_physics *ephys  ///< extra physics stuff.
	  )
  :
  nv_prim(nv)
//,  MPv3(nv,ntracer,trtype,ephys)
{
#ifdef TESTING
  cout <<"mpv8_TwoTempIso constructor setting up.\n";
#endif

  k_B = GS.kB();  // Boltzmann constant.
  m_p = GS.m_p(); // Proton mass.
  gamma   = SimPM.gamma;   // Gas has a constant ratio of specific heats.
  gamma_minus_one = SimPM.gamma -1.0;
  EP = ephys;

  nvl     = 3;    // three local variables.
  lv_nH   = 0;    // local H number density.
  lv_H0   = 1;    // x(H0).
  lv_eint = 2;    // E_{int}.

  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  //
  double X = EP->H_MassFrac;
  mean_mass_per_H = m_p/X; // this is mass per H nucleon.
  //
  // Number of ionised particles per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NION  = 1.0; // +0.25*EP->Helium_MassFrac/X;

  //
  // Number of electrons per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NELEC = 1.0; // +0.25*EP->Helium_MassFrac/X;

  //
  // The metal MassFrac is unused, so we use that to decide if the
  // neutral medium is neutral or molecular
  // If Metal_MassFrac >0.5 then neutral H is molecular (larger mu).
  // Otherwise, neutral H is atomic (smaller mu).
  //
  if (EP->Metal_MassFrac > 0.5) {
    TTI_Mol = 0.5;
  }
  else {
    TTI_Mol = 1.0;
  }
  TTI_Nnt = TTI_Mol +0.25*EP->Helium_MassFrac/X;

  TTI_Thi = EP->MaxTemperature;
  TTI_Tlo = EP->MinTemperature;
  Min_NeutralFrac     = 1.0e-12;
  Max_NeutralFrac     = 1.0-Min_NeutralFrac;

  //
  // Assume H+ fraction is the first tracer.
  //
  pv_Hp = nv-ntracer;
  if (pv_Hp<5) rep.error("mpv8_TwoTempIso::C bug: pv_Hp",pv_Hp);

#ifdef TESTING
  cout <<"mpv8_TwoTempIso: Y="<< EP->Helium_MassFrac;
  cout <<", Z="<< EP->Metal_MassFrac <<", mmpH="<<mean_mass_per_H;
  cout <<", NION="<< JM_NION <<", NELEC="<< JM_NELEC<<"\n";
#endif // TESTING
  return;
}

// ##################################################################
// ##################################################################

mpv8_TwoTempIso::~mpv8_TwoTempIso()
{
#ifdef TESTING
  cout <<"mpv8_TwoTempIso destructor.\n";
#endif
  return;
}


// ##################################################################
// ##################################################################


int mpv8_TwoTempIso::convert_prim2local(
          const double *p_in, ///< primitive vector from grid cell (length nv_prim)
          double *p_local
          )
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  p_local[lv_nH] = p_in[RO]/mean_mass_per_H;

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0]   = 1.0-p_in[pv_Hp];
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  //
  // internal energy is not used, so this is not really needed.
  //
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);


#ifdef TESTING
  //
  // Check for NAN/INF
  //
  for (int v=0;v<2;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics",p_local[v]);
  }
  if (p_local[lv_nH]<0.0 || !isfinite(p_local[lv_nH]))
    rep.error("Bad density input to mpv8_TwoTempIso::convert_prim2local",p_local[lv_nH]);
#endif // TESTING
  
  return 0;
}


// ##################################################################
// ##################################################################


int mpv8_TwoTempIso::convert_local2prim(
            const double *p_local,
            const double *p_in, ///< input primitive vector from grid cell (length nv_prim)
            double *p_out       ///< updated primitive vector for grid cell (length nv_prim)
            )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];

  //
  // Set output H+ fraction
  //
  p_out[pv_Hp] = 1.0-p_local[lv_H0];
  p_out[pv_Hp] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_out[pv_Hp]));

  // 
  // Output pressure is determined by ion fraction and number dens.:
  // - H+ fraction is y,
  // - Total number density, n_tot = nH*((1-y)*Nnt +y*(Nelec+Nion)),
  // - Temperature is T=Tlo + y*(Thi-Tlo) (from function),
  // so p = n_tot*k*T.
  //
  p_out[PG] = get_ntot(p_local[lv_nH],p_out[pv_Hp])*k_B
              *get_temperature(p_local[lv_nH],0.0,p_out[pv_Hp]);

#ifdef TESTING
  //cout <<"nH="<< p_local[lv_nH] <<", xp="<< p_out[pv_Hp] <<", ntot=";
  //cout <<get_ntot(p_local[lv_nH],p_out[pv_Hp])<<"\n";
#endif

#ifdef TESTING
  if (p_out[pv_Hp]<0.0 || p_out[pv_Hp]>1.0*(1.0+JM_RELTOL) || !isfinite(p_out[pv_Hp]))
    rep.error("Bad output H+ value in mpv8_TwoTempIso::convert_local2prim",p_out[pv_Hp]-1.0);
  if (p_out[PG]<0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in mpv8_TwoTempIso::convert_local2prim",p_out[PG]);
#endif // TESTING

  return 0;
}


// ##################################################################
// ##################################################################


double mpv8_TwoTempIso::get_temperature(
    const double,    ///< nH (per c.c.): UNUSED
    const double,    ///< E_int (per unit volume): UNUSED
    const double xp  ///< x(H+)
    )
{
  //
  // returns gas temperature according to T=(2yT_hi+(1-y)*T_lo)/(1+y)
  // which is RJR William's estimate of the correct temperature to
  // get the mixed cell physics right.
  //
  //double frac=0.1;
  //if (xp>frac) return TTI_Thi;
  //else return (xp/frac*(2.0*TTI_Thi-TTI_Tlo) +TTI_Tlo)/(1.0+xp/frac);

  // Dodgy logarithmic interpolation.
  //if (xp<1.0e-4) return TTI_Tlo;
  //else return TTI_Tlo +(TTI_Thi-TTI_Tlo)*(log10(xp)+4.0)/4.0;

  //
  // This is my estimate for gas temperature to try to account for
  // the mixed-cell physics in unresolved ionisation fronts.
  // It smoothes out the pressure, but in the end provides a worse
  // solution than other forms.
  //
  //return TTI_Tlo*TTI_Thi/(TTI_Thi-xp*(TTI_Thi-TTI_Tlo));

  //
  // returns gas temperature according to T=Tlo + y*(Thi-Tlo),
  //
  return (TTI_Tlo +xp*(TTI_Thi-TTI_Tlo));
}



// ##################################################################
// ##################################################################



double mpv8_TwoTempIso::Temperature(
            const double *pv, ///< primitive vector
            const double      ///< eos gamma
            )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    cout <<"MPv3::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];
  //convert_prim2local(pv,P);
  return get_temperature(0, 0, pv[pv_Hp]);
}



// ##################################################################
// ##################################################################




int mpv8_TwoTempIso::Set_Temp(
          double *p_pv,   ///< primitive vector.
          const double T, ///< temperature
          const double g  ///< eos gamma.
          )
{
  //
  // Check for T out of range
  //
  if      (T>TTI_Thi) {
    return (Set_Temp(p_pv,TTI_Thi,g));  
  }
  else if (T<TTI_Tlo) {
    return (Set_Temp(p_pv,TTI_Tlo,g));  
  }
  //
  // set x(H+) according to x=(T-T_lo)/(T_hi-T_lo).
  //
  p_pv[pv_Hp] = (T-TTI_Tlo)/(TTI_Thi-TTI_Tlo);
  //
  // now convert to local variables and back to primitive, and in the
  // 2nd conversion pressure will be set according to
  // get_temperature().
  //
  double P[nvl];
  int err = convert_prim2local(p_pv,P);
  err += convert_local2prim(P, p_pv, p_pv);
  return err;
}



// ##################################################################
// ##################################################################



double mpv8_TwoTempIso::get_ntot(
    const double nH, ///< nH
    const double xp  ///< x(H+) N.B. This is ion fraction, not neutral fraction.
    )
{
  //
  // This allows for molecular H neutral gas, with TTI_Mol, which is
  // 0.5 if molecular.  This is the (H0/H2) + (He) + (elect.+ions).
  //
  return
    ((1.0-xp)*TTI_Mol +(TTI_Nnt-TTI_Mol) +xp*(JM_NELEC+JM_NION))*nH;
}



// ##################################################################
// ##################################################################



double mpv8_TwoTempIso::get_recombination_rate(
          const int,          ///< ion index in tracer array (optional).
          const double *p_in, ///< input state vector (primitive).
          const double g      ///< EOS gamma (optional)
          )
{
#ifdef FUNCTION_ID
  cout <<"mpv8_TwoTempIso::get_recombination_rate()\n";
#endif // FUNCTION_ID
  double rate=0.0;
  double P[nvl];
  //
  // First convert to local variables.
  //
  convert_prim2local(p_in,P);
  //
  // Now get rate
  //
  rate = 2.7e-13*P[lv_nH]*P[lv_nH]*(1.0-P[lv_H0])*(1.0-P[lv_H0])*JM_NELEC;
  //cout <<"rate="<<rate<<"\n";

#ifdef FUNCTION_ID
  cout <<"mpv8_TwoTempIso::get_recombination_rate()\n";
#endif // FUNCTION_ID
  return rate;
}



// ##################################################################
// ##################################################################


int mpv8_TwoTempIso::TimeUpdateMP_RTnew(
          const double *p_in, ///< Primitive Vector to be updated.
          const int N_heat,   ///< Number of UV heating sources.
          const std::vector<struct rt_source_data> &heat_src,
            ///< list of UV-heating column densities and source properties.
          const int N_ion,    ///< number of ionising radiation sources.
          const std::vector<struct rt_source_data> &ion_src,
            ///< list of ionising src column densities and source properties.
          double *p_out,  ///< Destination Vector for updated values
                          ///< (can be same as first Vector).
          const double dt, ///< Time Step to advance by.
          const double,    ///< EOS gamma.
          const int, ///< Switch for what type of integration to use.
                     ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
          double *random_stuff ///< final temperature (not strictly needed).
          )
{
#ifdef FUNCTION_ID
  cout <<"mpv8_TwoTempIso::TimeUpdateMP_RTnew()\n";
#endif // FUNCTION_ID

  if (N_ion!=1) rep.error("No ionising sources!",N_ion);

  //
  // First set local variables for state vector and radiation source
  // properties.
  //
  int err=0;
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to MPv3::TimeUpdateMP_RTnew()",err);
  }

  //
  // The photoionisation equilibrium "radiation source" that
  // contains the fractional and cumulative attenuation for the cell
  // should be the first ionising source.
  //
  // so R = col, dR = delcol.
  //
  double R  = ion_src[0].Column;
  double dR = ion_src[0].DelCol;
  //if (SimPM.timestep>3000) cout <<"R="<<R<<", dR="<<dR;//<<"\n";
  //if      (R+dR >= 2.0) P[lv_H0] = 1.0;
  //else if (R+dR <= 1.0) P[lv_H0] = 0.0;
  //else                  P[lv_H0] = R+dR-1.0;

  if      (R >= 1.0) {
    //
    // all photons absorbed before ray reaches cell, so cell is fully
    // neutral.
    // 
    P[lv_H0] = 1.0;
  }
  else if (R+dR <= 1.0) {
    //
    // some photons pass through cell, so cell is fully ionised.
    //
    P[lv_H0] = 0.0;
  }
  else {
    //
    // I-front is within cell, so do a linear interpolation assuming
    // the I-front is perfectly sharp.  R<1, R+dR>1.
    //
    P[lv_H0] = (dR-1.0+R)/dR;
  }
  //if (SimPM.timestep>3000) cout <<", Nt.frac.="<<P[lv_H0]<<"\n";
  //
  // That's it!  The conversion function sets the output pressure
  // based on the neutral fraction.
  //
  err = convert_local2prim(P,p_in,p_out);

#ifdef FUNCTION_ID
  cout <<"mpv8_TwoTempIso::TimeUpdateMP_RTnew()\n";
#endif // FUNCTION_ID
  return err;
}


// ##################################################################
// ##################################################################



