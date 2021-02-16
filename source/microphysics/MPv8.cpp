///
/// \file MPv8.h
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running calculations with a simple
/// heating and cooling prescription for photoionization calculations
/// in the StarBench Workshop test problems.
///
/// Modifications:
/// - getting it written: mods up until 2014.06.XX
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.07 JM: New trtype array structure in constructor.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifdef LEGACY_CODE

#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "microphysics/MPv8.h"

using namespace std;

// ##################################################################
// ##################################################################

MPv8::MPv8(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr,   ///< Number of tracer variables in state vector.
    const std::string* tracers,   ///< List of what the tracer variables mean.
    struct which_physics* ephys,  ///< pointer to extra physics flags.
    struct rad_sources* rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    MPv3(nd, csys, nv, ntr, tracers, ephys, rsrcs, g)
{
#ifdef TESTING
  cout << "MPv8 constructor setting up.\n";
#endif
  gamma_minus_one = eos_gamma - 1.0;

  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  //
  double X        = 1.0 - EP->Helium_MassFrac;
  mean_mass_per_H = m_p / X;  // this is mass per H nucleon.
  //
  // Number of ionised particles per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NION = 1.0;  // +0.25*EP->Helium_MassFrac/X;

  //
  // Number of electrons per ionised H nucleon.
  // (assume He is inert, so no photons can ionise it).
  //
  JM_NELEC = 1.0;  // +0.25*EP->Helium_MassFrac/X;

  //
  // The metal MassFrac is unused, so we use that to decide if the
  // neutral medium is neutral or molecular
  // If Metal_MassFrac >0.5 then neutral H is molecular (larger mu).
  // Otherwise, neutral H is atomic (smaller mu).
  //
  if (EP->Metal_MassFrac > 0.5) {
    SBHC_Mol = 0.5;
  }
  else {
    SBHC_Mol = 1.0;
  }
  SBHC_Nnt = SBHC_Mol + 0.25 * EP->Helium_MassFrac / X;

  //
  // Set the high and low temperature variables to equal the cooling
  // rate at these temperatures, so that we can impose two equilibrium
  // temperatures.
  //
  double T   = EP->MaxTemperature;
  SBHC_EEqHi = 2.0e-19 * exp(-1.184e5 / (T + 1.0e3))
               + 2.8e-28 * sqrt(T) * exp(-92.0 / T);
  T          = EP->MinTemperature;
  SBHC_EEqLo = 2.0e-19 * exp(-1.184e5 / (T + 1.0e3))
               + 2.8e-28 * sqrt(T) * exp(-92.0 / T);
  cout << "\tPI-heating=" << SBHC_EEqHi << ", NT-heating=" << SBHC_EEqLo
       << "\n";

#ifdef TESTING
  cout << "MPv8: Y=" << EP->Helium_MassFrac;
  cout << ", Z=" << EP->Metal_MassFrac << ", mmpH=" << mean_mass_per_H;
  cout << ", NION=" << JM_NION << ", NELEC=" << JM_NELEC << "\n";
#endif  // TESTING
  return;
}

// ##################################################################
// ##################################################################

MPv8::~MPv8()
{
#ifdef TESTING
  cout << "MPv8 destructor.\n";
#endif
  return;
}

// ##################################################################
// ##################################################################

int MPv8::convert_prim2local(
    const double* p_in,  ///< primitive vector from grid cell (length nv_prim)
    double* p_local)
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  mpv_nH = p_in[RO] / mean_mass_per_H;

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0] = 1.0 - p_in[pv_Hp];
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  //
  // internal energy per unit volume
  //
  p_local[lv_eint] = p_in[PG] / (gamma_minus_one);

#ifdef TESTING
  //
  // Check for NAN/INF
  //
  for (int v = 0; v < 2; v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics", p_local[v]);
  }
  if (mpv_nH < 0.0 || !isfinite(mpv_nH))
    rep.error("Bad density input to MPv8::convert_prim2local", mpv_nH);
#endif  // TESTING

  return 0;
}

// ##################################################################
// ##################################################################

int MPv8::convert_local2prim(
    const double* p_local,
    const double*
        p_in,      ///< input primitive vector from grid cell (length nv_prim)
    double* p_out  ///< updated primitive vector for grid cell (length nv_prim)
)
{
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];

  //
  // Set output H+ fraction
  //
  p_out[pv_Hp] = 1.0 - p_local[lv_H0];
  p_out[pv_Hp] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_out[pv_Hp]));

  //
  // Output pressure is determined by ion fraction and number dens.:
  // - H+ fraction is x,
  // - Total number density, n_tot = nH*((1-x)*Nnt +x*(Nelec+Nion)),
  // - Pressure is E_int*(gamma-1)
  //
  p_out[PG] = p_local[lv_eint] * (gamma_minus_one);

#ifdef TESTING
  if (p_out[pv_Hp] < 0.0 || p_out[pv_Hp] > 1.0 * (1.0 + JM_RELTOL)
      || !isfinite(p_out[pv_Hp]))
    rep.error(
        "Bad output H+ value in MPv8::convert_local2prim", p_out[pv_Hp] - 1.0);
  if (p_out[PG] < 0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in MPv8::convert_local2prim", p_out[PG]);
#endif  // TESTING

  return 0;
}

// ##################################################################
// ##################################################################

double MPv8::get_temperature(
    const double nH,  ///< nH (per c.c.)
    const double E,   ///< E_int (per unit volume)
    const double xp   ///< x(H+)
)
{
  //
  // simple ideal gas EOS
  //
  return gamma_minus_one * E / (k_B * get_ntot(nH, xp));
}

// ##################################################################
// ##################################################################

double MPv8::get_ntot(
    const double nH,  ///< nH
    const double xp   ///< x(H+)
)
{
  return ((1.0 - xp) * SBHC_Mol + (SBHC_Nnt - SBHC_Mol)
          + xp * (JM_NELEC + JM_NION))
         * nH;
}

// ##################################################################
// ##################################################################

int MPv8::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double*          ///< extra user-data vector (UNUSED)
)
{

#ifdef TESTING
  // cout <<"MPv8::ydot(): Y="<< EP->Helium_MassFrac;
  // cout <<", Z="<< EP->Metal_MassFrac <<", mmpH="<<mean_mass_per_H;
  // cout <<", NION="<< JM_NION <<", NELEC="<< JM_NELEC;
  // cout <<", Nnt="<<SBHC_Nnt<<"\n";
#endif  // TESTING

  //
  // Set initial values
  //
  double OneMinusX = max(NV_Ith_S(y_now, lv_H0), Min_NeutralFrac);
  double x_in      = 1.0 - OneMinusX;
  double ne        = JM_NELEC * x_in * mpv_nH;
  double T         = get_temperature(mpv_nH, NV_Ith_S(y_now, lv_eint), x_in);

  double temp1 = 0.0;
  // double temp2=0.0;

  //
  // oneminusx_dot is in units of 1/s
  //
  double oneminusx_dot = 0.0;
  //
  // Edot is irrelevant here, so just leave it at zero so it never
  // affects the timesteps.
  //
  double Edot = 0.0;

  //
  // photo-ionisation of H:
  // photoionisation rate uses equation 18 in Mellema et al. 2006 (C2-ray
  // paper), noting that their Gamma is the rate per neutral H, so we multiply
  // by 1-x, as in their equation 11.
  //
  if (N_ion_srcs) {
    //
    // set current cell dTau0
    //

    temp1 = mpv_nH * mpv_delta_S * OneMinusX
            * Hi_monochromatic_photo_ion_xsection(JUST_IONISED);

    switch (ion_src_type) {
      case RT_EFFECT_MFION:
        //
        // Rather than divide the discretised rate by n(H0) and then
        // multiply by (1-x) to get oneminusx_dot, we simply divide by
        // n(H) since this is more numerically stable.  To do this, n(H)
        // is passed to the rate function instead of n(H0).
        //
        // Also, instead of using mpv_dTau0 for the cell optical depth,
        // we use the current optical depth (nH*mpv_delta_S*OneMinusX)
        // so that we allow the photoionisation rate to decrease as the
        // number of neutral atoms decreases during the timestep.
        //
        temp1 = Hi_discrete_multifreq_photoion_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_delta_S, mpv_Vshell);
        oneminusx_dot -= temp1;
        //
        // Instead of doing the multifrequency heating, we heat the gas
        // to get a specific ionized gas temperature, as for
        // monochromatic radiation.
        //
        // Edot += Hi_discrete_multifreq_photoheating_rate(mpv_Tau0,
        // temp1, mpv_nH,
        //                                    mpv_delta_S, mpv_Vshell);
        Edot += temp1 * SBHC_EEqHi / 2.7e-13;
        break;

      case RT_EFFECT_PION_MONO:
        //
        // hardcoded for a hv-13.6eV = 5.0eV monochromatic source.
        //
#define PHOTON_ENERGY 2.98e-11  // 5 eV
#define EXCESS_ENERGY 8.01e-12
        //#define PHOTON_ENERGY 2.24e-11
        //#define EXCESS_ENERGY 0.64e-12
        //#define PHOTON_ENERGY 2.499e-11  // 2 eV
        //#define EXCESS_ENERGY 3.204e-12
        temp1 = Hi_discrete_mono_photoion_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_NIdot, PHOTON_ENERGY, mpv_delta_S,
            mpv_Vshell);
        oneminusx_dot -= temp1;
        //
        // Excess energy is calculated by using heating/cooling balance
        // at T=EP->MaxTemperature in a fully-ionized plasma.
        //
        // Edot += temp1*EXCESS_ENERGY;
        Edot += temp1 * SBHC_EEqHi / 2.7e-13;
        // cout <<"PI-HR="<<temp1*EXCESS_ENERGY<<"\n";
        break;

      default:
        rep.error("Bad ion_src_type in dYdt()", ion_src_type);
        break;
    }  // switch
  }

  //
  // radiative recombination of H+, with hardcoded Recomb. rate.
  //
  // oneminusx_dot += Hii_rad_recomb_rate(T) *x_in*ne;
  oneminusx_dot += 2.7e-13 * x_in * ne;

  //
  // Cooling: Koyama and Inutsuka: 2002, ApJL, 564, 97 (KI02)
  // KI02 cooling has a cooling rate propto n^2 and a heating rate propto n.
  // It seems like this rate is independent of the ion fraction, so they must
  // assume a transition temperature between ionised and neutral, and this
  // must be implicitly included in their cooling rate function.
  //
  // First the cooling rate:  (KI02, eq.4) multiplied by n^2 to get volume
  // rate (see eq.3), and corrected by V-S etal. 2007, ApJ, 657, 870.
  // For n<10^6 per c.c.,  heating dominates for T>5K, so MinTemp is set to
  // 5K. If T<MinTemp there is no point calculating the cooling, so we skip
  // it.
  //
  Edot -= mpv_nH
          * (2.0e-19 * exp(-1.184e5 / (T + 1.0e3))
             + 2.8e-28 * sqrt(T) * exp(-92.0 / T));

  //
  // Artificial heating function for cold gas, scaling with nH*nH, to
  // get a certain equilibrium temperature independent of density.
  //
  Edot += mpv_nH * SBHC_EEqLo * EP->MinTemperature / T;

  //
  // now multiply Edot by nH to get units of energy loss/gain per unit volume
  // per second.
  //
  Edot *= mpv_nH;

  NV_Ith_S(y_dot, lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot, lv_eint) = Edot;
  return 0;
}

// ##################################################################
// ##################################################################

double MPv8::get_recombination_rate(
    const int,           ///< ion index in tracer array (optional).
    const double* p_in,  ///< input state vector (primitive).
    const double g       ///< EOS gamma (optional)
)
{
#ifdef FUNCTION_ID
  cout << "MPv8::get_recombination_rate()\n";
#endif  // FUNCTION_ID
  double rate = 0.0;
  double P[nvl];
  //
  // First convert to local variables.
  //
  convert_prim2local(p_in, P);
  //
  // Now get rate
  //
  rate = 2.7e-13 * mpv_nH * mpv_nH * (1.0 - P[lv_H0]) * (1.0 - P[lv_H0])
         * JM_NELEC;

#ifdef FUNCTION_ID
  cout << "MPv8::get_recombination_rate()\n";
#endif  // FUNCTION_ID
  return rate;
}

// ##################################################################
// ##################################################################

#endif  // LEGACY_CODE
