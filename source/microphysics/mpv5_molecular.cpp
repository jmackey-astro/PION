///
/// \file mpv5_molecular.cpp
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for modelling photoevaporation of dense molecular
/// clouds.
/// The file inherits from mpv3, and instead of the Wolfire et al.
/// (2003) neutral gas heating/cooling rates, it uses the Henney et
/// al. (2009) molecular cooling rate and heating rates.  So the only
/// significant difference is in the Ydot() function.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15
/// - 2013.03.21 JM: Fixed Helium free-free to use X(He).
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.

#include "microphysics/mpv5_molecular.h"
#include "global.h"

using namespace std;



// ##################################################################
// ##################################################################


mpv5_molecular::mpv5_molecular(
          const int nv,              ///< Total number of variables in state vector
          const int ntracer,         ///< Number of tracer variables in state vector.
          const std::string &trtype,  ///< List of what the tracer variables mean.
          struct which_physics *ephys  ///< extra physics stuff.
	  )
:
  mp_explicit_H(nv,ntracer,trtype,ephys)
{
#ifdef TESTING
  cout <<"mpv5_molecular constructor setting up.\n";
#endif
  return;
}

// ##################################################################
// ##################################################################

mpv5_molecular::~mpv5_molecular()
{
#ifdef TESTING
  cout <<"mpv5_molecular destructor.\n";
#endif
  return;
}


// ##################################################################
// ##################################################################


int mpv5_molecular::ydot(
          double,               ///< current time (UNUSED)
          const N_Vector y_now, ///< current Y-value
          N_Vector y_dot,       ///< vector for Y-dot values
          const double *        ///< extra user-data vector (UNUSED)
          )
{
  //
  // fixes min-neutral-fraction to Min_NeutralFrac
  //
  double OneMinusX = max(NV_Ith_S(y_now,lv_H0),Min_NeutralFrac);
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double x_in      = 1.0-OneMinusX;
  double ne        = JM_NELEC*x_in*mpv_nH;

  //
  // First get the temperature.  We assume the total particle number density
  // is given by 1.1*nH*(1+x_in), appropriate for a gas with 10% Helium by 
  // number, and if He is singly ionised whenever H is.
  //
  double T = get_temperature(mpv_nH, E_in, x_in);


  double temp1=0.0, temp2=0.0;
  double oneminusx_dot=0.0; // oneminusx_dot is in units of 1/s
  double Edot=0.0; // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH at the end.

  //
  // We set a minimum electron density based on the idea that Carbon is singly
  // ionised in low density gas.  y(C)=1.5e-4 in the gas phase (by number)
  // (Sofia,1997), approximately, so I add this to the electron number density
  // with an exponential cutoff at high densities.
  //
  ne += mpv_nH*1.5e-4*METALLICITY*exp(-mpv_nH/1.0e4);

  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0
  //
  Hi_coll_ion_rates(T, &temp1, &temp2);
  oneminusx_dot -= temp1*ne*OneMinusX; // the nH is divided out on both sides.
  Edot -= temp2*ne*OneMinusX;

  //
  // photo-ionisation of H:
  // photoionisation rate uses equation 18 in Mellema et al. 2006 (C2-ray paper),
  // noting that their Gamma is the rate per neutral H, so we multiply by 1-x, as
  // in their equation 11.
  //
  if (N_ion_srcs) {
    //
    // set current cell dTau0
    //
    temp1 = mpv_nH*mpv_delta_S*OneMinusX*
            Hi_monochromatic_photo_ion_xsection(JUST_IONISED);

    switch (ion_src_type) {
      case RT_EFFECT_PION_MULTI:
      //
      // Rather than divide the discretised rate by n(H0) and then multiply by (1-x) to
      // get oneminusx_dot, we simply divide by n(H) since this is more numerically stable.  To
      // do this, n(H) is passed to the rate function instead of n(H0).
      //
      // Also, instead of using mpv_dTau0 for the cell optical depth, we use the current
      // optical depth (nH*mpv_delta_S*OneMinusX) so that we allow the photoionisation
      // rate to decrease as the number of neutral atoms decreases during the timestep.
      // This is more stable.
      //
      oneminusx_dot -= Hi_discrete_multifreq_photoion_rate(mpv_Tau0, temp1,
                                  mpv_nH, mpv_delta_S, mpv_Vshell);
      Edot += Hi_discrete_multifreq_photoheating_rate(mpv_Tau0, temp1,
                                  mpv_nH, mpv_delta_S, mpv_Vshell);
      break;

      case RT_EFFECT_PION_MONO:
      //
      // hardcoded for a hv-13.6eV = 5.0eV monochromatic source.
      //
#define PHOTON_ENERGY 2.98e-11
#define EXCESS_ENERGY 8.01e-12
//#define PHOTON_ENERGY 2.24e-11
//#define EXCESS_ENERGY 0.64e-12
      temp1 = Hi_discrete_mono_photoion_rate(mpv_Tau0, temp1, mpv_nH*OneMinusX, mpv_NIdot, 
                                             PHOTON_ENERGY, mpv_delta_S, mpv_Vshell)*OneMinusX;
      oneminusx_dot -= temp1;
      Edot += temp1*EXCESS_ENERGY;
      break;

      default:
      rep.error("Bad ion_src_type in dYdt()",ion_src_type);
      break;
    } // switch
  }

  //
  // radiative recombination of H+
  //
  oneminusx_dot += Hii_rad_recomb_rate(T) *x_in*ne;
  //
  // Total H+ cooling: recombination plus free-free
  //
  Edot -= Hii_total_cooling(T) *x_in*ne;

  // ************ HENNEY+2009 COOLING FOR REAL SIMS *************
  //
  // Add Helium free-free (Z^2*n(He)/n(H) = X(He)/X(H) of the H+ free-free rate)
  // The normalisation is scaled so that I multiply by ne*nHp to get the 
  // correct cooling rate (i.e. the abundance of He is included in the prefactor).
  //
  Edot -= 1.68e-27*EP->Helium_MassFrac/(1.0-EP->Helium_MassFrac)*sqrt(T)*x_in*ne;

  //
  // collisional excitation cooling of H0 Aggarwal (1983) and Raga+(1997,ApJS).
  //
  Edot -= Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10);

  //
  // --------- END OF HYDROGEN COOLING, MOVING TO METAL COOLING --------
  //
  //
  // Now we get to the less certain elements of the cooling/heating function.
  // First we do the heating:
  //
  if (N_diff_srcs) {
    //
    // UV heating due to both diffuse radiation and point source radiation.
    // The quantity mpv_G0_UV is as defined in Henney et al. (2009) Appendix A1, Eq.A3,
    // and is set in set_parameters_for_current step()
    //
    //cout <<"adding diffuse heating! ";
    Edot += 1.9e-26*METALLICITY*mpv_G0_UV/(1.0+6.4*(mpv_G0_UV/mpv_nH));
    //cout <<"  DfUV="<<1.9e-26*METALLICITY*mpv_G0_UV/(1.0+6.4*(mpv_G0_UV/mpv_nH));

    //
    // IR heating (HAdCM09 eq.A6) from point source and/or diffuse radiation.
    // There is a different G0 parameter because the attenuation is according to 
    // exp(-0.05Av) rather than before where the coefficient was 1.9.
    //
    Edot += 7.7e-32*METALLICITY*mpv_G0_IR/pow(1.0+3.0e4/mpv_nH,2.0);
    //cout <<"  DfIR="<<7.7e-32*METALLICITY*mpv_G0_IR/pow(1.0+3.0e4/mpv_nH,2.0)<<"\n";
  }

  //
  // X-ray heating (HAdCM09 eq.A5)
  // Massive stars have x-ray luminosities of ~1.0e32 erg/s, so use this.
  // NOT IMPLEMENTED YET... NEEDS SOMETHING MORE CAREFUL.
  //
  ////Edot += 6.0e9*mpv_delta_S/mpv_Vshell;

  //
  // Cosmic ray heating (HAdCM09 eq.A7).
  //
  Edot += 5.0e-28*OneMinusX;

  //
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar neighbourhood.
  //
  oneminusx_dot -= 1.8e-17*OneMinusX;
  
  //
  // Now COOLING:
  // First forbidden line cooling of e.g. OII,OIII, dominant in HII regions.
  // This is collisionally excited lines of photoionised metals. (HAdCM09 eq.A9)
  // I have exponentially damped this at high temperatures! This was important!
  // Oxygen abundance set to 4.90e-4 from Asplund+(2009,ARAA).
  //
  temp1 = 1.42e-22*METALLICITY *exp(-33610.0/T -(2180.0*2180.0/T/T)) *x_in*ne*exp(-T*T/5.0e10);

  //
  // Collisionally excited lines of neutral metals: (HAdCM09 eq.A10).
  // Assumes the neutral metal fraction is the same as neutral H fraction.
  // Oxygen abundance set to 4.90e-4 from Asplund+(2009,ARAA).
  //
  temp1+= 2.19e-23*METALLICITY *exp(-28390.0/T -(1780.0*1780.0/T/T)) *ne*OneMinusX;

  //
  // Now the Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve.
  // We take the actual cooling rate to be the max of SD93-CIE and the
  // previous two terms.
  //
  temp2 = cooling_rate_SD93CIE(T) *x_in*x_in*mpv_nH*METALLICITY;
  Edot -= max(temp1,temp2);


  //
  // Finally "molecular" or PDR cooling from (HAdCM09 eq.A14), scaling with rho^{1.6}.
  // I have modified the Henney equation by multiplying by the neutral fraction
  // squared and also exponentially cutting off the function for T>10^5K.
  //
  temp1 = 70.0 +220.0*pow(mpv_nH/1.0e6, 0.2);
  temp2 = 3.981e-27*METALLICITY*pow(mpv_nH,0.6)*sqrt(T)*exp(-temp1/T);
  //
  // ********************* HACK FOR LOW DENSITY COOLING *******************
  //
  // This makes the cooling rate scale with nH^2 at low densities, and matches
  // the high density nH^1.6 scaling at nH=100 per c.c., at which point the 
  // scaling changes discontinuously (but the function is obviously continuous!).
  // 0.1585 is 100^(-0.4).
  //
  if (mpv_nH<100.0) temp2 *= 0.1585*pow(mpv_nH,0.4);
  // ********************* HACK FOR LOW DENSITY COOLING *******************
  //
  //
  // This should only apply to neutral gas, so multiply it by neutral-frac
  // squared, and also exponentially cut off the cooling above 10^5K.
  //
  temp2 *= OneMinusX*OneMinusX*exp(-T*T/1e10);
  Edot -= temp2;

  //
  // now multiply Edot by nH to get units of energy loss/gain per unit volume per second.
  //
  Edot *= mpv_nH;
#ifdef HIGHDENS_CUTOFF
  if (Edot<0.0) Edot *= exp(-mpv_nH*mpv_nH/1.0e6);
#endif //HIGHDENS_CUTOFF 

  //
  // We want to limit cooling as we approach the minimum temperature, so we scale
  // the rate to linearly approach zero as we reach Tmin.
  //
  if (Edot<0.0 && T<2.0*EP->MinTemperature) {
    Edot = min(0.0, (Edot)*(T-EP->MinTemperature)/SimPM.EP.MinTemperature);
  }

  NV_Ith_S(y_dot,lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot,lv_eint) = Edot;
  return 0;
}




