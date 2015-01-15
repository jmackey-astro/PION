///
/// \file mp_v2_aifa.cc
/// \author Jonathan Mackey
/// \date 2011.03.15
///
/// Description:
/// This class is an update on the microphysics class used for my thesis.  
///
/// - It uses multi-frequency photoionisation including spectral hardening with
///   optical depth, using the method outlined in Frank & Mellema
///   (1994,A&A,289,937), and updated by Mellema et al. (2006,NewA,11,374).
///
/// - It uses the Hummer (1994,MN,268,109) rates for radiative recombination and
///   its associated cooling, and Bremsstrahlung cooling.
///
/// - For collisional ionisation the function of Voronov (1997,ADNDT,65,1) is used.
///
/// - Collisional excitation of neutral H, table from Raga, Mellema, & Lundqvist
///   (1997,ApJS,109,517), using data from Aggarwal (1993).
///
/// - For cooling due to heavy elements, which are not explicitly included, we use
/// the CIE cooling curve of Sutherland & Dopita (1993,ApJS,88,253), but this time I
/// subtract the zero-metals curve from the solar-metallicity curve so that I only
/// have cooling due to metals.  This eliminates the potential for double-counting
/// which was in my previous cooling function.
///
/// - Then I use various formulae from Henney et al. (2009,MN,398,157) for cooling due
/// to collisional excitation of photoionised O,C,N (eq.A9), collisional excitation of
/// neutral metals (eq.A10), and the Wiersma+2009 CIE metals-only curve.  I take the max of
/// the WSS09 curve and the Henney et al. functions, to avoid double counting.  For 
/// neutral gas I can use the cooling curve of Henney et al 2009 eq.A14.
///
/// - Photoheating from ionisation is discussed above.  Cosmic ray heating will use a
/// standard value, X-ray heating is ignored.  UV heating due to the interstellar 
/// radiation field (ISRF) is according to the optical depth from the edge of the 
/// domain to the cell in question, using e.g. HEA09 eq.A3, if requested.  UV heating
/// from the star uses the same equation, but with the optical depth from the 
/// source (using a total H-nucleon column density).
///
///
/// The integration method uses the CVODES solver from the SUNDIALS package by
/// (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in Physics, 10, 138) available from 
/// https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton iteration.
///
/// Electrons, ions, ion fraction:  This modules makes the (crude) approximation that
/// He is an identical atom to H -- that it is only ever singly ionised and its ionisation
/// and recombination rates are identical to H.  But because most photons are below the
/// ionisation potential of He0 the opacity of He0 atoms is not counted in the optical 
/// depth of the medium to ionising photons.  This is not quite self-consistent, but it is
/// a better approximation than assuming He absorbs 14eV photons.
///
/// Modifications:
/// - 2011.03.29 JM: Wrote class mp_rates_ExpH_ImpMetals() functions.
/// - 2011.03.31 JM: Added solid angle vector for diffuse radiation.
///                  Finished coding, fixed a lot of bugs, need to test it now.
/// - 2011.04.12 JM: Added some 'todo's.
/// - 2011.04.14 JM: Fixed bugs; I am testing the non-RT part now.  Looks good so far.
/// - 2011.04.15 JM: Bugfixes
/// - 2011.04.17 JM: Added ifdefs for RT_TESTING to omit various processes in dYdt().
///     Debugging.  Added ifdef to integrate neutral fraction instead of ion fraction.
/// - 2011.04.18 JM: Bugfixes.
/// - 2011.04.22 JM: Debugged UV/IR heating.  Disabled X-ray heating.
/// - 2011.05.02 JM: Updated convert_prim2local to correct negative pressure/ion fraction
///    inputs rather than bugging out.
/// - 2011.05.02 JM: Added set_multifreq_source_properties() function
/// - 2011.05.04 JM: Added a discretised multifreq photoionisation rate with
///    an approximation for dtau<<1.  Fixed bugs, simplified code.
/// - 2011.05.06 JM: Set cooling limiting as T approaches SimPM.EP.MinTemperature so we
///    don't cool to too low a temperature.
/// - 2011.05.10 JM: Output cooling rates only if myrank==0 for parallel (so processes
///    don't fight over the file and slow down the code (by a lot!)).
/// - 2011.05.25 JM: Fixed a bad bug -- UV heating was not kicking in unless there was a 
///    diffuse sources, so the pt-src UV flux was useless up to now.
///    Got rid of "int-ion-frac" dydt function.  Too confusing to have two functions.
///    If for some reason I want to integrate the ion fraction in future I can go back to
///    an older revision.
/// - 2011.06.21 JM: Fixed bug in UV heating where Tau(FUV) was incorrectly calculated
///    Simplified timestep-limit calculation for xdot, since 2nd order integration is more 
///    accurate and so requires less tuning and allows larger timesteps.
/// - 2011.07.03 JM: Changed how Helium is treated appproximately (it is not ignored in
///    the photoionisation optical depths, which is a better approx. than including it).
/// - 2011.07.12 JM: Fixed the timescales function to have a less restrictive timestep
///    limitation based only on dt=0.15/|xdot|.
/// - 2011.09.21 JM: Set the cooling to be C2 (cooling=15 in cooling.cc) for RT_TEST_PROBS
/// - 2011.10.13 JM: Added switch in TimeUpdate_RTnew() so that if nothing is changing
///    much over a timestep then just do an Euler integration.
/// - 2015.01.15 JM: Added new include statements for new PION version.
///
/// NOTE: Oxygen abundance is set to 5.81e-4 from Lodders (2003,ApJ,591,1220)
///       which is the 'proto-solar nebula' value. The photospheric value is lower
///       4.9e-4, and that is used by Wiersma et al. (2009,MN,393,99).
///

#ifdef MP_V2_AIFA

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#ifndef EXCLUDE_MPV2

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/mp_v2_aifa.h"
#include "global.h"
using namespace std;

#define JM_NELEC 1.1 ///< ionised gas has 1.1 electrons per hydrogen nucleon (H+,He+, no He2+!)
#define JM_NION  1.1 ///< There are 1.1 nucleons per H nucleon (10% He by number).
#define JM_RELTOL 1.0e-2    ///< relative-error tolerance is set to 0.01 of this value.
#define JM_MINFRAC 1.0e-12  ///< min ion fraction i care about.

//#define MPV2_DEBUG

// ----------------------------------------------------------------------------
// -----------------  START OF mp_rates_ExpH_ImpMetals  -----------------------
// ----------------------------------------------------------------------------

mp_rates_ExpH_ImpMetals::mp_rates_ExpH_ImpMetals()
{
  //
  // This class has only two variables to integrate, Y(H+) and E_int.
  //
  n_eq = 2;
  n_xd = 0;

  k_B = 1.381e-16; // Boltzmann constant.
  m_p = 1.67e-24; // mass of proton.
  //
  // We want to set up the CIE cooling function for metals-only from WSS09
  // (i.e. with H+He cooling subtracted out).
  //
#ifdef RT_TEST_PROBS
  setup_SD93_cie();
#else
  setup_WSS09_CIE_OnlyMetals();
#endif

  //
  // initialise all the radiation variables to values that limit their heating
  // and cooling abilities.
  //
  nH = 1.0;
  Vshell = 1.0e54;
  NH0 = 1.0e24;
  dNH0 = 1.0e18;
  G0_UV = 0.0;
  G0_IR = 0.0;
  delta_S = 0.0;
  NIdot = 0.0;

  ion_src_type=-1;
  ion = diff = -1;
  return; 
}

mp_rates_ExpH_ImpMetals::~mp_rates_ExpH_ImpMetals()
{
  return;
}

void mp_rates_ExpH_ImpMetals::set_gamma_and_srcs(
          const double g,
          const int d, // No diffuse sources  --> 0, otherwise --> 1
          const int i  // No ionising sources --> 0, otherwise --> 1
          )
{
  gamma = SimPM.gamma; // This is a global constant in any simulation (for now!)
  gamma_minus_one = gamma-1.0;

  diff = d;
  ion  = i;

  //
  // If we have an ionising source we need to get the luminosity.
  //
  int got=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MONO) {
      NIdot = SimPM.RS.sources[isrc].strength;
      ion_src_type = RT_EFFECT_PION_MONO;
      got++;
      cout <<"mp_rates_ExpH_ImpMetals::set_gamma_and_srcs: Nisrc="<<ion<<", got="<<got<<": ";
      cout <<" found mono-p-ion src id="<<isrc<<" with NIdot="<<NIdot<<"\n";
    }
    else if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI) {
      ion_src_type = RT_EFFECT_PION_MULTI;
      NIdot = SimPM.RS.sources[isrc].strength;
      //cout <<"Please code for multi-frequency photoionisation!\n";
      //rep.error("multifreq not implemented","mp_rates_ExpH_ImpMetals::set_gamma_and_srcs");
      got++;
      cout <<"mp_rates_ExpH_ImpMetals::set_gamma_and_srcs: Nisrc="<<ion<<", got="<<got<<": ";
      cout <<" found multifreq-p-ion src id="<<isrc<<" with L="<<NIdot<<"\n";
    }
  }
  if (ion>0) {
    if (got!=1) rep.error("mismatch between ion and got",ion-got);
  }

  return;
}



void mp_rates_ExpH_ImpMetals::get_error_tolerances(
                double &reltol, ///< relative error tolerance.
                std::vector<double> &atol ///< absolute error tolerances
                )
{
  reltol = 0.01*JM_RELTOL;
  atol[0] = 0.1*JM_MINFRAC; ///< minimum ion fraction I care about.
  atol[1] = 1.0e-17; ///< for n=1.0, T=1.0e4, ==> E=2.07e-12, so say 1e-17?
  return;
}



void mp_rates_ExpH_ImpMetals::set_parameters_for_current_step(
          const std::vector<double> &P ///< list of parameters in an array.
          )
{
  //
  // We need number density of H nucleons, n(H).
  //
  nH = P[0];
  //
  // Next if we have a point radiation source we need Vshell.
  //
  Vshell = P[1];
  //
  // If the point source is ionising, we need N(H0), dN(H0).  The 
  // luminosity is already set in the MP class constructor by a call to
  // Setup_photoionisation_rate_table(Tstar,Rstar,etc).
  // Most H-ionising photons can't ionise He, so best that it doesn't contribute to the opacity.
  // BUT KEEP IN MIND THIS IS AN APPROXIMATION!
  //
  NH0 = P[2];
  dNH0 = P[3];
  //
  // If there is a UV heating point source, we need the UV heating flux for
  // direct absorption of UV radiation (G0_UV), and for absorption of
  // UV radiation which has been re-radiated in the far-IR (G0_IR).
  //
  G0_UV = P[4];
  G0_IR = P[5];

  //
  // Path length through cell.
  //
  delta_S = P[6];

  //
  // If there is a diffuse radiation source, we need A_v (and maybe
  // source intensity at some stage), due to dust extinction and hence
  // proportional to total column density to cell from all sided of the domain.
  //
  //Av_diffuse = P[6];
  //
  // EOS Gamma:
  //
  //gamma  = P[7];
  //gamma_minus_one = gamma-1.0;

//#ifdef RT_TESTING
#ifdef MPV2_DEBUG
  cout <<"MPR set params.  nH="<<nH<<", Vsh="<<Vshell<<", NH0="<<NH0<<", dNH0=";
  cout <<dNH0<<", G0_UV="<<G0_UV<<", G0_IR="<<G0_IR<<", delta_S="<<delta_S<<", Nion="<<ion<<"\n";
#endif
//#endif
  
  return;
}




void mp_rates_ExpH_ImpMetals::get_problem_size(
                  int *ne, ///< number of equations
                  int *np  ///< number of parameters in user_data vector.
                  )
{
  *ne = n_eq;
  *np = 0;
  return;
}




//
// returns gas temperature according to E=nkT/(g-1) with n=nH*(1.1+1.1*x_in),
// appropriate for a gas with 10% Helium by number, and if He is singly ionised
// whenever H is ionised (n_i=1.1n_H, n_e=1.1n_H).
//
double mp_rates_ExpH_ImpMetals::get_temperature(
    const double n, ///< nH
    const double E, ///< E_int
    const double x  ///< x(H+)
    )
{
  double ntotk= k_B*n*(JM_NION+JM_NELEC*x);
  //cout <<"mp_rates_ExpH_ImpMetals::get_temperature(): n.K="<<T<<", T="<<gamma_minus_one*E/T<<"\n";

  return gamma_minus_one*E/ntotk;
}




// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
int mp_rates_ExpH_ImpMetals::dYdt(
          const double OneMinusX, // neutral fraction!
          const double E_in,
          double *xdot, // change in neutral fraction!
          double *Edot
          )
{
  //
  //  1-x
  //
  //double OneMinusX = max(1.0-x_in, 1.0e-14); // fixes min-neutral-fraction to 1e-14
  double x_in = 1.0-OneMinusX;
  double ne  = JM_NELEC*x_in*nH;
  //double ni  = JM_NION*x_in*nH;
  //double nH0 = OneMinusX*nH;
  //double nHp = x_in*nH;

  //
  // First get the temperature.  We assume the total particle number density
  // is given by 1.1*nH*(1+x_in), appropriate for a gas with 10% Helium by 
  // number, and if He is singly ionised whenever H is.
  //
  double T = get_temperature(nH, E_in, x_in);


  double temp1=0.0, temp2=0.0;
  *xdot=0.0; // xdot is in units of 1/s
  *Edot=0.0; // Edot is calculated in units of erg/s per H nucleon, multiplied by nH at the end.

#ifdef RT_TEST_PROBS
  if (SimPM.EP.coll_ionisation) {
#endif 

    //
    // collisional ionisation of H, with its associated cooling.
    // scales with n_e*nH0
    //
    Hi_coll_ion_rates(T, &temp1, &temp2);
    *xdot -= temp1*ne*OneMinusX; // the nH is divided out on both sides.
    *Edot -= temp2*ne*OneMinusX;
#ifdef MPV2_DEBUG
    cout <<"T="<<T<<" x_in="<<x_in<<" E_in="<<E_in<<"    CI="<<-temp2*ne*OneMinusX<<"\n";
#endif // MPV2_DEBUG

#ifdef RT_TEST_PROBS
  } // if coll_ionisation
#endif
  
  //
  // photo-ionisation of H:
  // photoionisation rate uses equation 18 in Mellema et al. 2006 (C2-ray paper),
  // noting that their Gamma is the rate per neutral H, so we multiply by 1-x, as
  // in their equation 11.
  //
  if (ion) {
    switch (ion_src_type) {
      case RT_EFFECT_PION_MULTI:
      //
      // Rather than divide the discretised rate by n(H0) and then multiply by (1-x) to
      // get xdot, we simply divide by n(H) since this is more numerically stable.  To
      // do this, n(H) is passed to the rate function instead of n(H0).
      //
      // Also, instead of using dNH0 for the cell optical depth, we use the current
      // optical depth (nH*delta_S*OneMinusX) so that we allow the photoionisation
      // rate to decrease as the number of neutral atoms decreases during the timestep.
      // This is more stable.
      //
      *xdot -= Hi_discrete_multifreq_photoion_rate(    NH0, nH*delta_S*OneMinusX,
                                                   nH, delta_S, Vshell);
      *Edot += Hi_discrete_multifreq_photoheating_rate(NH0, nH*delta_S*OneMinusX,
                                                   nH, delta_S, Vshell);
      //cout <<"multifreq!   PI-rate="<<(temp1-temp2)<<"\n";
#ifdef MPV2_DEBUG
      cout <<"  PI="<<(temp1-temp2);
#endif // MPV2_DEBUG
      break;

      case RT_EFFECT_PION_MONO:
      //
      // hardcoded for a hv-13.6eV = 5.0eV monochromatic source.
      //
#define PHOTON_ENERGY 2.98e-11
#define EXCESS_ENERGY 8.01e-12
      temp1 = Hi_discrete_mono_photoion_rate(NH0, nH*delta_S*OneMinusX, nH*OneMinusX, NIdot, 
                                             PHOTON_ENERGY, delta_S, Vshell)*OneMinusX;
      *xdot -= temp1;
      *Edot += temp1*EXCESS_ENERGY;
#ifdef MPV2_DEBUG
      cout <<" PIR: NIdot="<<NIdot<<" nH="<<nH<<" NH0="<<NH0;
      cout <<" dNH0="<<dNH0<<" ds="<<delta_S<<" Vsh="<<Vshell;
      cout <<" 1-x="<<OneMinusX<<"... RATE="<<temp1<<"\n";
#endif
#ifdef MPV2_DEBUG
      cout <<"  xdot_PI="<<*xdot<<"\n";
#endif
      break;

      default:
      rep.error("Bad ion_src_type in dYdt()",ion_src_type);
      break;
    } // switch
  }

#ifdef RT_TEST_PROBS
  if (SimPM.EP.rad_recombination) {
    *xdot += Hii_rad_recomb_rate(T) *x_in*x_in*nH;
    //*Edot -= Hii_total_cooling(T) *x_in*x_in*nH;
    *Edot -= Hii_rad_recomb_cooling(T) *x_in*x_in*nH;
#else
    //
    // radiative recombination of H+
    //
    *xdot += Hii_rad_recomb_rate(T) *x_in*ne;
    //
    // Total H+ cooling: recombination plus free-free
    //
    *Edot -= Hii_total_cooling(T) *x_in*ne;
#endif
#ifdef MPV2_DEBUG
    cout <<"  RR+FF="<<-Hii_total_cooling(T)<<" xin="<<x_in<<" ne="<<ne<<"\t\t";
#endif // MPV2_DEBUG

#ifdef RT_TEST_PROBS
  } // if rad_recombination
  //else {
  //  // need free-free cooling at least
  //  *Edot -= 1.68e-27*sqrt(T) *x_in*ne;
  //}
#endif

#ifndef RT_TEST_PROBS // ************ HENNEY+2009 COOLING FOR REAL SIMS *************
  //
  // Add Helium free-free (Z^2*n(He)/n(H) = 0.4 of the H+ free-free rate)
  // The normalisation is scaled so that I multiply by ne*nHp to get the 
  // correct cooling rate (i.e. the abundance of He is included in the prefactor).
  //
  *Edot -= 6.72e-28*sqrt(T) *x_in*ne;
#ifdef MPV2_DEBUG
  cout <<"  He-FF="<<-6.72e-28*sqrt(T) *x_in*ne;
#endif // MPV2_DEBUG

  //
  // collisional excitation cooling of H0 Aggarwal (1993)
  //
  *Edot -= Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne;
#ifdef MPV2_DEBUG
  cout <<"  CExH0="<<-Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne;
#endif // MPV2_DEBUG
  //
  // --------- END OF HYDROGEN COOLING, MOVING TO METAL COOLING --------
  //
  // Now we get to the sketchier elements of the cooling/heating function.
  // First we do the heating:
  //

  if (diff) {
    //
    // UV heating due to both diffuse radiation and point source radiation.
    // The quantity G0_UV is as defined in Henney et al. (2009) Appendix A1, Eq.A3,
    // and is set in set_parameters_for_current step()
    //
    //cout <<"adding diffuse heating!\n";
    *Edot += 1.9e-26*G0_UV/(1.0+6.4*(G0_UV/nH));
#ifdef MPV2_DEBUG
    cout <<"  DfUV="<<1.9e-26*G0_UV/(1.0+6.4*(G0_UV/nH));
#endif // MPV2_DEBUG

    //
    // IR heating (HAdCM09 eq.A6) from point source and/or diffuse radiation.
    // There is a different G0 parameter because the attenuation is according to 
    // exp(-0.05Av) rather than before where the coefficient was 1.9.
    //
    *Edot += 7.7e-32*G0_IR/pow(1.0+3.0e4/nH,2.0);
#ifdef MPV2_DEBUG
    cout <<"  DfIR="<<7.7e-32*G0_IR/pow(1.0+3.0e4/nH,2.0);
#endif // MPV2_DEBUG
  }

  //
  // X-ray heating (HAdCM09 eq.A5)
  // Massive stars have x-ray luminosities of ~1.0e32 erg/s, so use this.
  // THIS IS A HACK JUST TO GET IT GOING! FIX IT LATER.
  //
  //*Edot += 6.0e9*delta_S/Vshell;
#ifdef MPV2_DEBUG
  //cout <<"  XRay="<<6.0e9*delta_S/Vshell;
#endif // MPV2_DEBUG

  //
  // Cosmic ray heating (HAdCM09 eq.A7)
  //
  *Edot += 5.0e-28;
#ifdef MPV2_DEBUG
  cout <<"  CR="<< 5.0e-28;
#endif // MPV2_DEBUG
  
  //
  // Now COOLING:
  // First forbidden line cooling of e.g. OII,OIII, dominant in HII regions.
  // This is collisionally excited lines of photoionised metals. (HAdCM09 eq.A9)
  // I have exponentially damped this at high temperatures! This was important!
  // Oxygen abundance set to 5.81e-4 from Lodders et al. (2003,ApJ,591,1220,Tab.2).
  //
  temp1 = 1.69e-22 *exp(-33610.0/T -(2180.0*2180.0/T/T)) *x_in*ne *exp(-T*T/5.0e10);
#ifdef MPV2_DEBUG
  cout <<"  CExMi="<<-1.69e-22 *exp(-33610.0/T -(2180.0*2180.0/T/T)) *x_in*ne *exp(-T*T/5.0e10);
#endif // MPV2_DEBUG
  //
  // Collisionally excited lines of neutral metals: (HAdCM09 eq.A10).
  // Assumes the neutral metal fraction is the same as neutral H fraction.
  // Oxygen abundance set to 5.81e-4 from Lodders et al. (2003,ApJ,591,1220,Tab.2).
  //
  temp1+= 2.60e-23 *exp(-28390.0/T -(1780.0*1780.0/T/T)) *ne*OneMinusX;
#ifdef MPV2_DEBUG
  cout <<"  CExMn="<<-2.60e-23 *exp(-28390.0/T -(1780.0*1780.0/T/T)) *ne*OneMinusX;
#endif // MPV2_DEBUG


  //
  // Now the Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve.
  // We take the actual cooling rate to be the max of SD93-CIE and the
  // previous two terms.
  //
  temp2 = cooling_rate_SD93CIE(T) *x_in*x_in*nH;
#ifdef MPV2_DEBUG
  cout <<"  CIEc="<<-cooling_rate_SD93CIE(T) *x_in*x_in*nH<<"\n";
  cout <<"WSS09 ccoling rate="<<temp2<<", metal line cooling rate="<<temp1<<"\n";
#endif // MPV2_DEBUG
  *Edot -= max(temp1,temp2);

  //
  // Finally "molecular" or PDR cooling from (HAdCM09 eq.A14), scaling with rho^{1.6}.
  // TODO: Query whether this should be applied for all ion fractions, or only for
  //       very low ion fractions.
  //
  if (T<1.0e4 && x_in<0.2) {
    temp1 = 70.0 +220.0*pow(nH/1.0e6, 0.2);
    temp2 = 3.981e-27*pow(nH,0.6)*sqrt(T)*exp(-temp1/T);
#ifdef MPV2_DEBUG
    cout <<"Without PDR, Edot="<<(*Edot)*nH<<", PDR cooling rate="<<temp2*nH<<"\n";
#endif // MPV2_DEBUG
    *Edot -= temp2;
    //cout <<", PDR cooling rate="<<temp2*nH<<", nH="<<nH<<"\n";
  }

#else // RT_TEST_PROBS ******************** SIMPLE C2 COOLING FOR TESTS *************
  //
  // SD93-CIE cooling
  //
  //cout <<"new C2 cooling.\n";
  *Edot -= cooling_rate_SD93CIE(T) *x_in*x_in*nH;
  //
  // Forbidden line cooling
  //
  if (T<2.0e4 && T>100.0) 
    *Edot -= 2.0e-24*(T/8000.0)*x_in*x_in*nH;
  //
  // Toy model corresponding to case 15 in cooling.cc
  //
  temp1 = max(0.0, (1.0+x_in)*(1.0-x_in)*(1.0-x_in)*k_B*(T-100.0)/(SimPM.gamma-1.0)/3.16e11 *exp(-T/1.e4) );
  *Edot -= temp1;
  
#endif // NOT RT_TEST_PROBS ******************** SIMPLE C2 COOLING FOR TESTS *************



  //
  // now multiply Edot by nH to get units of energy loss/gain per unit volume per second.
  //
  *Edot *= nH;

  //
  // We want to limit cooling as we approach the minimum temperature, so we scale
  // the rate to linearly approach zero as we reach Tmin.
  //
  if (*Edot<0.0 && T<2.0*SimPM.EP.MinTemperature) {

#ifdef MPV2_DEBUG
    cout <<"limiting cooling: Edot="<<*Edot<<", T="<<T;
#endif // MPV2_DEBUG

    *Edot = min(0.0, (*Edot)*(T-SimPM.EP.MinTemperature)/SimPM.EP.MinTemperature);

#ifdef MPV2_DEBUG
    cout <<"... resetting Edot to "<<*Edot<<"\n";
#endif // MPV2_DEBUG

  }
#ifdef RT_TEST_PROBS
  if (!SimPM.EP.update_erg) {
    *Edot = 0.0;
  }
#endif

  return 0;
}
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------

// ----------------------------------------------------------------------------
// -----------------  END OF mp_rates_ExpH_ImpMetals  -------------------------
// ----------------------------------------------------------------------------

//
// This has to be defined somewhere, so may as well make it as local as possible.
//
class mp_rates_ExpH_ImpMetals MPR;





// ----------------------------------------------------------------------------
// -----------------  START OF mp_v2_aifa WRAPPER CLASS -----------------------
// ----------------------------------------------------------------------------

mp_v2_aifa::mp_v2_aifa(
          const int nv,              ///< Total number of variables in state vector
          const int ntracer,         ///< Number of tracer variables in state vector.
    const std::string &trtype  ///< List of what the tracer variables mean.
    )
:
  nv_prim(nv)
{
  cout <<"mp_v2_aifa: new microphysics class.\n";

  //
  // Set up tracer variables (i.e. just find which one is H+).
  //
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  int ftr = nv_prim -ntracer; // first tracer variable.
  string s;
  int len = (trtype.length() +5)/6 -1; // first 6 chars are the type, then list of tracers, each 6 chars long.
    cout <<"\t\ttrtype = "<<trtype<<"\n";
    cout <<"\t\tlen="<<len<<", ntr="<<ntracer<<"\n";
  if (len!=ntracer) {
    cout <<"warning: string doesn't match ntracer.  make sure this looks ok: "<<trtype<<"\n";
    //rep.error("string doesn't match ntracer",ntracer-len);
  }

  //
  // Find ionisation fraction in tracer variable list.
  //
  pv_Hp=-1;
  for (int i=0;i<len;i++) {
    s = trtype.substr(6*(i+1),6); // Get 'i'th tracer variable.
    if (s=="H1+___" || s=="HII__") {
      pv_Hp = ftr+i;
    }
  }
  if (pv_Hp<0)
    rep.error("No H ionisation fraction found in tracer list",trtype);
  
  //
  // We only have two local variables: ion fraction and internal energy density.
  //
  nvl     = 2;
  lv_Hp   = 0;
  lv_eint = 1;
  lv_nH   = 0.0;
  gamma = SimPM.gamma;
  gamma_minus_one = SimPM.gamma -1.0;

  //
  // Also set gamma for MPR, and flags for whether we have diffuse and ionising
  // radiation sources.
  //
  int ion=0, diff=0;
  for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
    if ((SimPM.RS.sources[isrc].type==RT_SRC_DIFFUSE) ||
        (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE && SimPM.RS.sources[isrc].effect==RT_EFFECT_UV_HEATING)
        ) diff++;
    if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE &&
        (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MONO ||
         SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI))  ion ++;
  }
  cout <<"\t\tmp_v2_aifa:: found "<<diff<<" diffuse and "<<ion<<" ionising sources.\n";
  MPR.set_gamma_and_srcs(gamma,diff,ion);

  k_B = GS.kB();
  m_p = GS.m_p();
  mean_mass_per_H = 1.40*m_p;
  Min_IonFrac=JM_MINFRAC;
  Max_IonFrac=1.0-JM_MINFRAC;
  
  // ----------------------- DIFFUSE RADIATION -------------------------------
  // Set up solid angles for diffuse radiation, whether or not they are needed.
  //
  diff_angle.clear();
  if      (SimPM.coord_sys==COORD_CRT && SimPM.ndim==3) {
    diff_angle.resize(6);
    for (int v=0;v<6;v++) diff_angle[v] = 4.0*M_PI/6.0;
  }
  else if (SimPM.coord_sys==COORD_CRT && SimPM.ndim==2) {
    diff_angle.resize(4);
    for (int v=0;v<4;v++) diff_angle[v] = 2.0*M_PI/4.0;
  }
  else if (SimPM.coord_sys==COORD_CRT && SimPM.ndim==1) {
    diff_angle.resize(2);
    for (int v=0;v<2;v++) diff_angle[v] = 1.0;
  }
  else if (SimPM.coord_sys==COORD_CYL && SimPM.ndim==2) {
    diff_angle.resize(3);
    //
    // for each source in turn, we get its direction and set the angle accordingly.
    //
    int count=0;
    int dir=-1;
    for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
      if (SimPM.RS.sources[isrc].type==RT_SRC_DIFFUSE) {
        for (int v=0;v<SimPM.ndim;v++) {
          if (SimPM.RS.sources[isrc].position[v] > 1.0e99) dir = 2*v+1;
          if (SimPM.RS.sources[isrc].position[v] <-1.0e99) dir = 2*v;
        }
        if (dir<0) rep.error("Diffuse source not at infinity!",isrc);
        //
        // if direction is in Z then angle is as for 3D, and if R+ then 4x3D values.
        //
        if (dir==ZNcyl || dir==ZPcyl) diff_angle[count] = 4.0*M_PI/6.0;
        else if (dir==RPcyl)          diff_angle[count] = 16.0*M_PI/6.0;
        else rep.error("Bad source direction",dir);
        count++;
      }
    }
    cout <<"Angles for diffuse sources: ["<<diff_angle[0]<<", ";
    cout <<diff_angle[1]<<", "<<diff_angle[2]<<"]\n";
  }
  else if (SimPM.coord_sys==COORD_SPH && SimPM.ndim==1) {
    //
    // for spherical symmetry the only diffuse source is at infinity
    //
    diff_angle.resize(1);
    diff_angle[0] = 4.0*M_PI;
  }
  else {
    rep.error("Unhandled coord-sys/ndim combination in mp_v2_aifa::mp_v2_aifa",SimPM.ndim);
  }
  // ----------------------- DIFFUSE RADIATION -------------------------------


  // ------------------------- IONISING SOURCE ----------------------
  if (ion) {
    if (ion>1) rep.error("too many ionising source in mp_v2_aifa()",ion);
    //
    // Need to set up the multifrequency tables if needed.
    //
    for (int isrc=0; isrc<SimPM.RS.Nsources; isrc++) {
      if (SimPM.RS.sources[isrc].type==RT_SRC_SINGLE) {
        if (SimPM.RS.sources[isrc].effect==RT_EFFECT_PION_MULTI) {
          int err=set_multifreq_source_properties(&SimPM.RS.sources[isrc]);
          if (err)
            rep.error("multifreq photoionisation setup failed in mp_v2_aifa const.",err);
        }
      }
    }
    
  } 
  // ------------------------- IONISING SOURCE ----------------------


  // --------------------------- CVODES ----------------------
  // Initialise the CVODES solver memory etc.
  //
  have_setup_cvodes=false;
  setup_cvodes_solver_without_Jacobian();
  // --------------------------- CVODES ----------------------

#ifdef MPV2_DEBUG
  //
  // ------------- output cooling rates for various temperatures ---------
  //
#ifdef PARALLEL
  if (mpiPM.myrank==0) {
#endif 
    double p[nv_prim];
    p[RO]=2.34e-24; p[PG]=1.0e-12;
    p[pv_Hp] = 0.99;

    string opfile("cooling_mpv2_aifa.txt");
    ofstream outf(opfile.c_str());
    if(!outf.is_open()) rep.error("couldn't open outfile",1);
    outf <<"Cooling Curve Data: Temperature(K) Rates(erg/cm^3/s) x=0.99999, x=0.00001, x=0.5 (n=1 per cc)\n";
    outf.setf( ios_base::scientific );
    outf.precision(6);
    double t=SimPM.EP.MinTemperature, Edi=0.0,Edn=0.0,Edpi=0.0,junk=0.0;
    do {
      p[pv_Hp] = 0.99999;
      Set_Temp(p,t,junk);
      MPR.dYdt(1.0-p[pv_Hp],p[PG]/gamma_minus_one,&junk,&Edi);

      p[pv_Hp] = 0.5;
      Set_Temp(p,t,junk);
      MPR.dYdt(1.0-p[pv_Hp],p[PG]/gamma_minus_one,&junk,&Edpi);

      p[pv_Hp] = 0.00001;
      Set_Temp(p,t,junk);
      MPR.dYdt(1.0-p[pv_Hp],p[PG]/gamma_minus_one,&junk,&Edn);
      outf << t <<"\t"<< Edi <<"  "<< Edn <<"  "<< Edpi <<"\n";
      t *=1.05;
    } while (T<min(1.0e9,SimPM.EP.MaxTemperature));
    outf.close();  
#ifdef PARALLEL
  }
#endif 
  //
  // ------------- output cooling rates for various temperatures ---------
  //
#endif // MPV2_DEBUG

  cout <<"mp_v2_aifa: Constructor finished and returning.\n";
  return;
}




mp_v2_aifa::~mp_v2_aifa()
{
  diff_angle.clear();
  //
  // Free vector memory
  //
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
  N_VDestroy_Serial(abstol);
  //
  // Free integrator memory
  //
  CVodeFree(&cvode_mem);
}




int mp_v2_aifa::Tr(const string t)
{
  return lv_Hp;
}




//
// Set the properties of a multifrequency ionising radiation source.
//
int mp_v2_aifa::set_multifreq_source_properties(
              const struct rad_src_info *rsi
              )
{
  //
  // Some sanity checks:
  // - make sure source is multi-freq and ionising
  // - make sure Rstar and Tstar are positive and finite
  //
  if (rsi->effect!=RT_EFFECT_PION_MULTI)
    rep.error("Source is not multi-frequency!", rsi->id);
  if (rsi->Rstar<0 || !isfinite(rsi->Rstar))
    rep.error("Source has bad Rstar parameter", rsi->Rstar);
  if (rsi->Tstar<0 || !isfinite(rsi->Tstar))
    rep.error("Source has bad Tstar parameter", rsi->Tstar);

  double mincol=SimPM.dx*1.0e-11, maxcol=1.0e24, Emax=1000.0*1.602e-12;
  int Nspl=150, Nsub=800;
  cout <<"#################### mp_v2_aifa::set_multifreq_source_properties() MinCol="<<mincol<<"\n";

  MPR.Setup_photoionisation_rate_table(
              rsi->Tstar,
              rsi->Rstar*6.96e10,
              rsi->strength,
              mincol,  maxcol, Emax,  Nsub, Nspl);

  return 0;
}




int mp_v2_aifa::TimeUpdateMP(
        const double *p_in,   ///< Primitive Vector to be updated.
        double *p_out,        ///< Destination Vector for updated values.
        const double dt,      ///< Time Step to advance by.
        const double,         ///< EOS gamma.
        const int,            ///< Switch for what type of integration to use.
        double *random_stuff  ///< Vector of extra data (column densities, etc.).
        )
{
  std::vector<struct rt_source_data> temp;
  int err = TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}




int mp_v2_aifa::TimeUpdate_RTsinglesrc(
        const double *, ///< Primitive Vector to be updated.
        double *,       ///< Destination Vector for updated values.
        const double,   ///< Time Step to advance by.
        const double,   ///< EOS gamma.
        const int,      ///< Switch for what type of integration to use.
        const double,   ///< flux in per unit length along ray (F/ds or L/dV)
        const double,   ///< path length ds through cell.
        const double,   ///< Optical depth to entry point of ray into cell.
        double *        ///< return optical depth through cell in this variable.
        )
{
  cout <<"mp_v2_aifa::TimeUpdate_RTsinglesrc() is not implemented!\n";
  return 1;
}





int mp_v2_aifa::TimeUpdateMP_RTnew(
                   const double *p_in, ///< Primitive Vector to be updated.
              const int N_heat,   ///< Number of UV heating sources.
                   const std::vector<struct rt_source_data> &heat_src,
                   ///< list of UV-heating column densities and source properties.
                   const int N_ion,    ///< number of ionising radiation sources.
                   const std::vector<struct rt_source_data> &ion_src,
                   ///< list of ionising src column densities and source properties.
       double *p_out,  ///< Destination Vector for updated values
                       ///< (can be same as first Vector.
       const double dt,   ///< Time Step to advance by.
       const double,   ///< EOS gamma.
       const int, ///< Switch for what type of integration to use.
                  ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
       double *random_stuff ///< final temperature (not strictly needed).
       )
{
  int err=0;
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to mp_v2_aifa::TimeUpdateMP_RTnew()",err);
  }

  setup_radiation_source_parameters(p_in,P, N_heat, heat_src, N_ion, ion_src);

  //
  // update radiation source properties, if needed (re-calculate multi-frequency
  // photoionisation rates if the source properties have changed).
  // TODO: CODE THIS SOMEWHERE!! BUT MAYBE put this somewhere else -- 
  //       update the source properties when they change,
  //       and through a different interface function!!
  //

  vector<double> Y0(nvl), Yf(nvl,0.0);
  Y0[lv_Hp]   = 1.0-P[lv_Hp];
  Y0[lv_eint] = P[lv_eint];

  double maxdelta=0.0;
  err = MPR.dYdt(Y0[lv_Hp],Y0[lv_eint], &(Yf[lv_Hp]), &(Yf[lv_eint]));
  if (err) 
    rep.error("dYdt() returned an error in mp_v2_aifa::TimeUpdateMP_RTnew()",err);
  for (int v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(Yf[v]*dt/Y0[v]));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < 0.05) {
    for (int v=0;v<nvl;v++) {
      Yf[v] = Y0[v] + dt*Yf[v];
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    err = integrate_cvodes_step(Y0, 0, 0.0, dt, Yf);
    if (err) {
      rep.error("integration failed: mp_v2_aifa::TimeUpdateMP_RTnew()",err);
    }
  }

  P[lv_Hp]   = 1.0-Yf[lv_Hp];
  P[lv_eint] = Yf[lv_eint];

  err = convert_local2prim(P,p_in,p_out);
  
  Y0.clear(); Yf.clear();
  return err;
}




double mp_v2_aifa::Temperature(
            const double *pv, ///< primitive vector
            const double      ///< eos gamma
            )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    cout <<"mp_v2_aifa::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  double P[nvl];
  convert_prim2local(pv,P);

#ifdef MP_DEBUG
  //double T1 = MPR.get_temperature(lv_nH, P[lv_eint], P[lv_Hp]);
  //cout <<"TEMPERATURE: nH="<<lv_nH<<", eint="<< P[lv_eint]<<", x="<<P[lv_Hp];
  //double T2=(gamma_minus_one)*P[lv_eint]/k_B/(JM_NION+JM_NELEC*P[lv_Hp])/lv_nH; // Temperature.
  //cout <<"\tT(MPR)="<<T1<<", T(local)="<<T2<<"\n";
#endif // MP_DEBUG

  return MPR.get_temperature(lv_nH, P[lv_eint], P[lv_Hp]);
}




int mp_v2_aifa::Set_Temp(
          double *p,     ///< primitive vector.
          const double T, ///< temperature
          const double  ///< eos gamma.
          )
{
  double P[nvl];  
  //
  // Check for negative pressure.  If density<0 then we should bug out because
  // there is no way to set a temperature, but if p<0 we can just overwrite it.
  //
  if (p[PG]<=0.0) {
    //cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p[PG] = 1.0e-12;  // It doesn't matter what this is as long as p>0
  }
  int err = convert_prim2local(p,P);
  P[lv_eint] = (JM_NION+JM_NELEC*P[lv_Hp])*lv_nH*k_B*T/(gamma_minus_one);
  err += convert_local2prim(P, p, p);

  return err;

}




int mp_v2_aifa::Init_ionfractions(
        double *, ///< Primitive vector to be updated.
        const double, ///< eos gamma.
        const double  ///< optional gas temperature to end up at. (negative means use pressure)
        )
{
  cout <<"mp_v2_aifa::Init_ionfractions() is not implemented! Write me!\n";
  return 1;
}




int mp_v2_aifa::convert_prim2local(
          const double *p_in, ///< primitive vector from grid cell (length nv_prim)
          double *p_local
          )
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);
  p_local[lv_Hp]   = p_in[pv_Hp];
  lv_nH = p_in[RO]/mean_mass_per_H;

  //
  // Check for negative ion fraction, and set to a minimum value if found.
  //
  if (p_local[lv_Hp]<0.0) {
    cout <<"mp_v2_aifa::convert_prim2local: negative ion fraction input: ";
    cout <<p_local[lv_Hp] <<", setting to 1.0e-12.\n";
    p_local[lv_Hp] = Min_IonFrac;
  }

  //
  // Check for negative pressure (note this shouldn't happen, so we output a
  // warning) and set to 10K if we find it.
  //
  if (p_local[lv_eint]<=0.0) {
    cout <<"mp_v2_aifa::convert_prim2local: negative pressure input: p=";
    cout <<p_local[lv_eint]<<", setting to 10K.\n";
    p_local[lv_eint] = (JM_NION+JM_NELEC*p_local[lv_Hp])*lv_nH*k_B*10.0/(gamma_minus_one);
  }

  //
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  p_local[lv_Hp] = max(Min_IonFrac, min(Max_IonFrac, p_local[lv_Hp]));

#ifdef MP_DEBUG
  //
  // Check for NAN/INF
  //
  for (int v=0;v<nvl;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics",p_local[v]);
  }
  if (lv_nH<0.0 || !isfinite(lv_nH))
    rep.error("Bad density input to mp_v2_aifa::convert_prim2local",lv_nH);
#endif // MP_DEBUG
  
  return 0;
}




int mp_v2_aifa::convert_local2prim(
            const double *p_local,
            const double *p_in, ///< input primitive vector from grid cell (length nv_prim)
            double *p_out      ///< updated primitive vector for grid cell (length nv_prim)
            )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];

  p_out[PG]    = p_local[lv_eint]*(gamma_minus_one);
  p_out[pv_Hp] = p_local[lv_Hp];

#ifdef MP_DEBUG
  if (p_out[pv_Hp]<0.0 || p_out[pv_Hp]>1.0*(1.0+JM_RELTOL) || !isfinite(p_out[pv_Hp]))
    rep.error("Bad output H+ value in mp_v2_aifa::convert_local2prim",p_out[pv_Hp]-1.0);
  if (p_out[PG]<0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in mp_v2_aifa::convert_local2prim",p_out[PG]);
#endif // MP_DEBUG

  //
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  p_out[pv_Hp] = max(Min_IonFrac, min(Max_IonFrac, p_out[pv_Hp]));

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected xHp from p_out[]).
  //
  double T = MPR.get_temperature(lv_nH, p_local[lv_eint], p_out[pv_Hp]);
  if (T>1.01*SimPM.EP.MaxTemperature) {
    //cout <<"mp_v2_aifa::convert_local2prim() HIGH temperature encountered. ";
    //cout <<"T="<<T<<", obtained from nH="<<lv_nH<<", eint="<<p_local[lv_eint]<<", x="<<p_out[pv_Hp];
    //cout <<"...  limiting to T="<<SimPM.EP.MaxTemperature<<"\n";
    Set_Temp(p_out,SimPM.EP.MaxTemperature,0);
  }
  if (T<0.99*SimPM.EP.MinTemperature) {
    //cout <<"mp_v2_aifa::convert_local2prim() LOW  temperature encountered. ";
    //cout <<"T="<<T<<", obtained from nH="<<lv_nH<<", eint="<<p_local[lv_eint]<<", x="<<p_out[pv_Hp];
    //cout <<"...  limiting to T="<<SimPM.EP.MinTemperature<<"\n";
    Set_Temp(p_out,SimPM.EP.MinTemperature,0);
  }

  return 0;
}



double mp_v2_aifa::timescales(
          const double *p_in, ///< Current cell state vector.
          const double,   ///< EOS gamma.
          const bool, ///< set to 'true' if including cooling time.
          const bool, ///< set to 'true' if including recombination time.
          const bool  ///< set to 'true' if including photo-ionsation time.
          )
{
  //cout <<"mp_v2_aifa::timescales() is not implemented! use new timescales fn.\n";
  //return 1.0e200;
#ifdef MP_DEBUG
  if (SimPM.RS.Nsources!=0) {
    cout <<"WARNING: mp_v2_aifa::timescales() using non-RT version!\n";
  }
#endif // MP_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin= timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}




///
/// This returns the minimum timescale of all microphysical processes, including
/// reaction times for each species and the total heating/cooling time for the gas.
/// It requires the radiation field as an input, so it has substantially greater
/// capability than the other timescales function.
///
double mp_v2_aifa::timescales_RT(
                    const double *p_in, ///< Current cell state vector.
                    const int N_heat,      ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &heat_src,
                    ///< list of UV-heating column densities and source properties.
                    const int N_ion,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &ion_src,
                    ///< list of ionising src column densities and source properties.
                    const double   ///< EOS gamma.
                    )
{
  int err=0;
  //
  // First convert to local variables.
  //
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to mp_v2_aifa::timescales_RT()",err);
  }

  //
  // Next give MPR the radiation properties of the current cell.
  //
  setup_radiation_source_parameters(p_in,P, N_heat, heat_src, N_ion, ion_src);

  //
  // Now calculate y-dot[]...
  //
  double xdot=0.0, Edot=0.0;
  err = MPR.dYdt(1.0-P[lv_Hp], P[lv_eint], &xdot, &Edot);
  if (err) {
    rep.error("dYdt() returned an error in mp_v2_aifa::timescales_RT()",err);
  }

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t=HUGEVALUE;
  //
  // First get the ionisation timescale.  Since the timestep obtained will eventually
  // be multiplied by 0.3, this is effectively dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the energy change 
  // (heating timescale) is not required for accurately tracking ionisation fronts 
  // (although it may be needed for cooling!).
  //
  t = min(t,0.8333/(fabs(xdot)+TINYVALUE)); // hard-code that x can't change by more than 0.25
#ifdef MP_DEBUG
  cout <<"MP timescales: xdot="<<xdot<<", Edot="<<Edot<<" t_x="<<t;
#endif // MP_DEBUG
  //
  // Now cooling/heating time (this will then be multiplied by 0.3 in
  // gridmethods.cc:IntUniformFV::calc_microphysics_dt() to limit to 30% change
  // in energy).  THIS IS NO LONGER NEEDED.
  //t = min(t,P[lv_eint]/(fabs(Edot)+TINYVALUE));
#ifdef MP_DEBUG
  cout <<" and min(t_x,t_e)="<<t<<",  "; rep.printVec("P[x,E]",P,nvl);
#endif // MP_DEBUG

  return t;
}





void mp_v2_aifa::setup_radiation_source_parameters(
                    const double *p_in, ///< primitive input state vector.
                    double *P,  ///< local input state vector (x_in,E_int)
                    const int N_heat, ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &heat_src,
                    ///< list of UV-heating column densities and source properties.
                    const int N_ion,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &ion_src
                    ///< list of ionising src column densities and source properties.
                    )
{
  //-------------------- RADIATION SOURCE INFO -----------------------
  //
  // First deal with the column densities and source strengths to get the UV heating
  // and EUV ionisation+heating rates.
  //
#ifdef MP_DEBUG
  if (heat_src.size() != static_cast<unsigned int>(N_heat)) {
    rep.error("Timescales: N_heating_srcs doesn't match vector size in Harpreet's MP integrator",
              heat_src.size());
  }
  if (ion_src.size() != static_cast<unsigned int>(N_ion)) {
    rep.error("Timescales: N_ionising_srcs doesn't match vector size in Harpreet's MP integrator",
              ion_src.size());
  }
#endif // MP_DEBUG

  //
  // first if we have an ionising source, we get the path length ds, from (x,nH), and
  // source.DelCol = rho*ds*(1-x)
  //
  double delta_s = 0.0;
  //
  // We set Vshell to huge value so that if there is no source, 
  // it will divide the effect to zero (if it gets used).
  //
  double Vshell  = 1.0e200; 
  bool single_src = false;
  if (N_ion>0) {
    if (N_ion>1)
      rep.error("Fix mp_v2_aifa to deal with more than one ionising point source",N_ion);
    delta_s = ion_src[0].dS;
    Vshell  = ion_src[0].Vshell;
    single_src = true;
  }
#ifndef MP_DEBUG
  // if debugging, we want to run this loop regardless of whether ds is already set.
  if (N_heat>0 && fabs(delta_s)<SMALLVALUE) {
#endif // MP_DEBUG
    //
    // if no ionising sources, there may still be a point UV-heating source.
    // here DelCol = rho*ds.
    //
    for (int v=0; v<N_heat; v++) {
      if (heat_src[v].type == RT_SRC_SINGLE) {
#ifdef MP_DEBUG
      if (P[lv_Hp]<0.99) {
        cout <<"setup_rad_src_params: heating:  ds="<<heat_src[0].DelCol/p_in[RO];
        cout <<", Vshell="<<heat_src[0].Vshell<<"\n";
      }
#endif // MP_DEBUG
        delta_s = heat_src[0].dS;
        Vshell  = heat_src[0].Vshell;
        single_src=true;
      }
    }
#ifndef MP_DEBUG
  }
#endif // MP_DEBUG


  //
  // Hard-coded to assume UV radiation opacity comes from dust with sigma=5e-22 cm2,
  // and that the input is \int \rho ds, integrated to the front edge of the cell.
  // Diffuse radiation strength is an intensity (i.e. per solid angle), so the total flux
  // going into the solver will be roughly: sum(I(Omega)*delta-Omega)exp(-1.9Av).
  // This comes from Henney et al. (2009) eq. A3
  // TODO: Check value of 1.9 for extinction, and the cross section of 5e-22
  //
  double UV_diffuse_flux=0.0, IR_diffuse_flux=0.0, temp=0.0;
  int i_diff=0;
  double Av_UV = 1.9*1.086*5.0e-22/mean_mass_per_H;
  double Av_IR = Av_UV*0.05/1.9;

  for (int v=0; v<N_heat; v++) {
    if (heat_src[v].type == RT_SRC_DIFFUSE) {
      temp = heat_src[v].strength *diff_angle[i_diff];
#ifdef MP_DEBUG
      cout <<"setup_rad_src_params:\tdiffuse src: id="<<heat_src[v].id<<" 1.9Av="<<Av_UV*heat_src[v].Column;
      cout <<", strength="<<heat_src[v].strength<<", angle="<<diff_angle[i_diff];
      cout <<": attenuated flux="<<temp*exp(-Av_UV*heat_src[v].Column)<<"\n";
#endif // MP_DEBUG
      UV_diffuse_flux += temp*exp(-Av_UV*heat_src[v].Column);
      IR_diffuse_flux += temp*exp(-Av_IR*heat_src[v].Column);
      i_diff++;
      //cout <<"UV_diff_flux="<<temp*exp(-Av_UV*heat_src[v].Column);
      //cout <<" Col="<<heat_src[v].Column<<" Av="<<Av_UV*heat_src[v].Column/1.9<<"\n";
    }
    else {
      //
      // This source must be a point source of UV heating. In this case the strength is 
      // the photon luminosity, so flux = L*ds*exp(-1.9Av)/Vshell
      //
      temp = heat_src[v].strength*delta_s/heat_src[v].Vshell;
#ifdef MP_DEBUG
      cout <<"setup_rad_src_params:\tpoint   src: id="<<heat_src[v].id<<" 1.9Av="<<Av_UV*heat_src[v].Column;
      cout <<", strength="<<heat_src[v].strength<<", ds="<<delta_s;
      cout <<", Vshell="<<heat_src[v].Vshell;
      cout <<": attenuated flux="<<temp*exp(-Av_UV*heat_src[v].Column)<<"\n";
#endif // MP_DEBUG
      UV_diffuse_flux += temp*exp(-Av_UV*heat_src[v].Column);
      IR_diffuse_flux += temp*exp(-Av_IR*heat_src[v].Column);
      //cout <<"UV_ptsc_flux="<<temp*exp(-Av_UV*heat_src[v].Column)<<"\n";
      //cout <<"UV_ptsc_flux="<<temp*exp(-Av_UV*heat_src[v].Column);
      //cout <<" Col="<<heat_src[v].Column<<" Av="<<Av_UV*heat_src[v].Column/1.9<<"\n";
    }
  } // loop over heating sources.

  //
  // now divide by 1.2e7 to get it normalised correctly for Will's equation A3.
  //
  UV_diffuse_flux /= 1.2e7;
  IR_diffuse_flux /= 1.2e7;

#ifdef MP_DEBUG
  cout <<"\tTotal UV attenuated flux = "<<UV_diffuse_flux<<" in units of 1.2e7 phot/cm2/s\n";
#endif // MP_DEBUG

  //
  // Now ionising sources (we should only have one!).
  //
  if (N_ion>1) rep.error("Code for more than one ionising source!",N_ion);

  //
  // Now put the RT data into a vector and pass to MPR rate-calculation class.
  //
  std::vector<double> params(7);
  params[0] = lv_nH;
  params[1] = Vshell; // this is a huge value if no point sources
  if (N_ion>0) {
    //
    // columns are rho*ds*(1-x), and we want nH*ds*(1-x), so divide by mu_H
    //
    params[2] = ion_src[0].Column/mean_mass_per_H;
    params[3] = ion_src[0].DelCol/mean_mass_per_H;
  }
  else {
    params[2] = 0.0;
    params[3] = 0.0;
  }
  //
  // Now we need diffuse flux (these are zero if no heating sources)
  //
  params[4] = UV_diffuse_flux;
  params[5] = IR_diffuse_flux;
  //
  // And dS, which is needed for something...
  //
  params[6] = delta_s;
  //
  // now send all these to the MPR class so it knows the parameters when it
  // calculates dY/dt
  //
#ifdef RT_TESTING
  if (P[lv_Hp]<0.99) {
      cout <<"mpv2aifa: ionising: ds="<<delta_s<<", Vshell="<<Vshell;
      cout <<": NH0="<<params[2]<<", dNH0="<<params[3]<<", nH="<<lv_nH<<"\n";
  }
#endif
  MPR.set_parameters_for_current_step(params);
  //
  //---------------- END RADIATION SOURCE INFO -----------------------
  //
  params.clear();
  return;
}




int mp_v2_aifa::setup_cvodes_solver_without_Jacobian()
{
  if (have_setup_cvodes) {
    cout <<"Error! Trying to setup CVODES solver twice!\n";
    return 1;
  }
  int err = 0;
  //
  // Get the number of equations to solve, and number of extra-data elements needed.
  // Should be [2,0].
  //
  n_eq=0; n_xd=0;
  MPR.get_problem_size(&n_eq, &n_xd);
  //cout <<"\t\tn_eq="<<n_eq<<"\n"; 
  //
  // Allocate memory for CVodes y-vectors, and err-tol vector.
  //
  y_in  = N_VNew_Serial(n_eq);
  y_out = N_VNew_Serial(n_eq);
  abstol = N_VNew_Serial(n_eq);
  

  //
  // Call CVodeCreate to create the solver memory and specify the 
  // Backward Differentiation Formula and the use of a Newton iteration.
  //
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (!cvode_mem) {
    cerr <<"setup_cvodes_solver() error: cvode_mem="<<cvode_mem<<"\n";
    return 2;
  }

  //
  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the inital time T0, and
  // the initial dependent variable vector y.
  //
  double t=0.0;
  err = CVodeInit(cvode_mem, Ydot_for_cvodes, t, y_in);
  if (err!= CV_SUCCESS) {
    cerr <<"setup_cvodes_solver() CVodeInit error: "<<err<<"\n";
    return 3;
  }

  //
  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  //
  double reltol=0.0;
  vector<double> atol(n_eq);
  MPR.get_error_tolerances(reltol,atol); // both args passed by reference.
  for (int v=0;v<n_eq;v++) {
    NV_Ith_S(abstol,v) = atol[v];
  }
  //cout <<"reltol="<<reltol<<", atol=["<<atol[0]<<", "<<atol[1]<<", "<<atol[2]<<"]\n";
  err = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (err!= CV_SUCCESS) {
    cerr <<"setup_cvodes_solver() CVodeSVtolerances: err="<<err<<"\n";
    return 4;
  }

  //
  // Call CVDense to specify the CVDENSE dense linear solver
  //
  err = CVDense(cvode_mem, n_eq);
  if (err != CVDLS_SUCCESS) {
    cerr <<"setup_cvodes_solver() CVDense(): err="<<err<<"\n";
    return 5;
  }

  //
  // Set the Jacobian routine to Jac (user-supplied) WE SKIP THIS STEP, B/C
  // THE JACOBIAN IS TOO COMPLICATED TO CALCULATE
  //
  //err = CVDlsSetDenseJacFn(cvode_mem, Jacobian_for_cvodes);
  //if (err != CVDLS_SUCCESS) {
  //  cerr <<"setup_cvodes_solver() CVDlsSetDenseJacFn: err="<<err<<"\n";
  //  return 6;
  //}

  // Should be all done now, so return.
  have_setup_cvodes=true;
  return 0;


}

int mp_v2_aifa::integrate_cvodes_step(
              const vector<double> &Y0, ///< input vector
              double *params, ///< parameters for user_data
              double t_now,   ///< start time.
              double dt,      ///< time-step.
              vector<double> &Yf ///< output vector.
              )
{
#ifdef MP_DEBUG
  if (!have_setup_cvodes) {
    cout <<"Please setup cvodes solver before integrating cvodes step!\n";
    return 1;
  }
#endif // MP_DEBUG
  int err = 0;
  if (Y0.size() != static_cast<unsigned int>(n_eq)) {
    cerr <<"integrate_cvodes_step() Y0 has wrong size: "<<Y0.size()<<"\n";
    return 1;
  }

  //
  // copy Y0 to cvodes struct y_in
  //
  for (int v=0; v<n_eq; v++) NV_Ith_S(y_in,v) = Y0[v];

  //
  // user data are parameters which won't change during the step, but 
  // may be passed in to integrator in variable params.
  //
  if (n_xd >0) {
    err = CVodeSetUserData(cvode_mem, reinterpret_cast<void *>(params));
    if (err != CV_SUCCESS) {
      cerr <<"integrate_cvodes_step() CVodeSetUserData: err="<<err<<"\n";
      return 3;
    }
  }

  //
  // Now a do-while loop to get to the end of the step.  We may need to do it
  // in steps, which is why it is in a loop.
  //
  double t_temp = t_now;
  double tf = t_now +dt;
  int fail_ct=0, step_ct=0;
  do {
    //
    // Re-initialise cvodes to the current Y-value and start-time.
    //
    err = CVodeReInit(cvode_mem, t_now, y_in);
    if (err != CV_SUCCESS) {
      cerr <<"integrate_cvodes_step() CVodeReInit(): err="<<err<<"\n";
      return 4;
    }
    //
    // integrate one timestep, returning answer to temp array y_out, and new time to 
    // temporary variable t_temp.
    //
    err = CVode(cvode_mem, t_now+dt, y_out, &t_temp, CV_NORMAL);
    if (err!=CV_SUCCESS) {
      cout <<"error with integrate: err="<<err<<", trying smaller step.\n";
      cout <<" old-t="<<t_now<<", returned t="<<t_temp<<", dt was = "<<dt<<"\n";
      //cout <<"FAILED loop ct="<<fail_ct<<": yout = [";
      //for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(y_out,v)<<", ";
      //cout << NV_Ith_S(y_out,n_eq-1) <<"]\n";
      //cout <<"FAILED loop ct="<<fail_ct<<":    y = [";
      //for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(y_in,v)<<", ";
      //cout << NV_Ith_S(y_in,n_eq-1) <<"]\n";
      //
      // step failed, so don't copy yout to y, do shrink the timestep,
      // and don't update t.
      // 
      dt /=2.0;
      t_temp = t_now;
      fail_ct++;
    }
    else {
      //
      // Must have succeeded, copy y_out to y, make the next timestep
      // larger (if needed), and make sure dt is not too large, 
      // or running past the end of the integration.
      //
      //cout <<"EEEEEE old-t="<<t_now<<", returned t="<<t_temp<<", dt was = "<<dt<<"\n";
      //cout <<"EEEEEE loop ct="<<fail_ct<<": yout = [";
      //for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(y_out,v)<<", ";
      //cout << NV_Ith_S(y_out,n_eq-1) <<"]\n";
      //cout <<"EEEEEE loop ct="<<fail_ct<<":    y = [";
      //for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(y_in,v)<<", ";
      //cout << NV_Ith_S(y_in,n_eq-1) <<"]\n";
      //
      for (int v=0;v<n_eq;v++) NV_Ith_S(y_in,v) = NV_Ith_S(y_out,v);
      t_now = t_temp;
      dt *= 1.5;
      dt = min(dt,tf-t_now);
      //
      // step succeeded, so reset fail_ct to zero and increment step_ct
      //
      fail_ct=0;
      step_ct++;
    }
    //dt = min(dt,tf-t);
    //dt = max(dt,1.0);
  } while (t_now<tf && fail_ct<=MAX_TRYS && step_ct<=MAX_STEPS);


  if (fail_ct>=MAX_TRYS) {
    cout <<"Integration failed after bisecting the step "<<MAX_TRYS<<"times.\n";
    return fail_ct;
  }
  if (step_ct>=MAX_STEPS) {
    cout <<"Integration took "<<MAX_STEPS<<" and still didn't get to end of step.  giving up.\n";
    return step_ct;
  }

#ifdef MP_DEBUG
  cout <<"Final version step_ct="<<step_ct<<": t="<<t_now<<", y = [";
  for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(y_out,v)<<", ";
  cout << NV_Ith_S(y_out,n_eq-1) <<"]\n";
#endif // MP_DEBUG

  for (int v=0;v<n_eq;v++) Yf[v] = NV_Ith_S(y_out,v);
  return 0;
}




// ----------------------------------------------------------------------------
// -----------------    END OF mp_v2_aifa WRAPPER CLASS -----------------------
// ----------------------------------------------------------------------------

int Ydot_for_cvodes(
          double t, ///< current time
          N_Vector y,   ///< current Y-value
          N_Vector yd, ///< vector for Y-dot values
          void *data    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          )
{
  //
  // Call the class member ydot(x_in, E_in, &x_out, &E_out) function:
  //
  int err = MPR.dYdt(NV_Ith_S(y,0), NV_Ith_S(y,1), &(NV_Ith_S(yd,0)), &(NV_Ith_S(yd,1)));

#ifdef MP_DEBUG
  //
  // cout to make sure it worked.
  //
  for (int v=0; v<2; v++) {
    cout <<"\t\tYDOT: y["<<v<<"] = "<< NV_Ith_S(y,v);
    cout <<", yd["<<v<<"] = "<< NV_Ith_S(yd,v)<<"\n";
  }
#endif // MP_DEBUG
  
  return err;
}

#endif // if not excluding MPv2
#endif // if  MP_V2_AIFA is defined

