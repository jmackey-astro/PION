///
/// \file mpv6_PureH.cpp
/// \author Jonathan Mackey
/// \date 2013.02.15
///
/// This class is for running the cosmological radiative transfer
/// test problems I (Iliev et al.,2006,MNRAS,371,1057) and II
/// (Iliev et al.,2009,MNRAS,400,1283).
///
/// The file inherits from mpv3, and instead of the Wolfire et al.
/// (2003) and Henney et al. (2009) heating/cooling rates, it uses
/// pure atomic hydrogen heating/cooling, and also sets the He and
/// metal abundances to zero (hardcoded!) regardless of the settings
/// in the parameterfile.
///
/// Modifications:
/// - getting it written: mods up until 2013.02.15
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/mpv6_PureH.h"



using namespace std;



// ##################################################################
// ##################################################################


mpv6_PureH::mpv6_PureH(
          const int nv,              ///< Total number of variables in state vector
          const int ntracer,         ///< Number of tracer variables in state vector.
          const std::string &trtype,  ///< List of what the tracer variables mean.
          struct which_physics *ephys  ///< extra physics stuff.
	  )
:
  mp_explicit_H(nv,ntracer,trtype,ephys)
{
#ifdef TESTING
  cout <<"mpv6_PureH constructor setting up.\n";
#endif
  //
  // Here we set JM_NELEC and JM_NION to 1.0 because there is only H.
  //
  EP->Helium_MassFrac = 0.0;
  EP->Metal_MassFrac = 0.0;
  JM_NION  = 1.0;
  JM_NELEC = 1.0;
  METALLICITY = 0.0;
  mean_mass_per_H = m_p;
  
#ifdef TESTING
  cout <<"mpv6_PureH: Y="<< EP->Helium_MassFrac;
  cout <<", Z="<< EP->Metal_MassFrac <<", mmpH="<<mean_mass_per_H;
  cout <<", NION="<< JM_NION <<", NELEC="<< JM_NELEC<<"\n";
#endif // TESTING
  return;
}

// ##################################################################
// ##################################################################

mpv6_PureH::~mpv6_PureH()
{
#ifdef TESTING
  cout <<"mpv6_PureH destructor.\n";
#endif
  return;
}


// ##################################################################
// ##################################################################


int mpv6_PureH::ydot(
          double,               ///< current time (UNUSED)
          const N_Vector y_now, ///< current Y-value
          N_Vector y_dot,       ///< vector for Y-dot values
          const double *        ///< extra user-data vector (UNUSED)
          )
{

#ifdef TESTING
  //cout <<"mpv6_PureH::ydot(): Y="<< EP->Helium_MassFrac;
  //cout <<", Z="<< EP->Metal_MassFrac <<", mmpH="<<mean_mass_per_H;
  //cout <<", NION="<< JM_NION <<", NELEC="<< JM_NELEC<<"\n";
#endif // TESTING

  //
  // fixes min-neutral-fraction to Min_NeutralFrac
  //
  double OneMinusX = max(NV_Ith_S(y_now,lv_H0),Min_NeutralFrac);
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double x_in      = 1.0-OneMinusX;
  double ne        = x_in*mpv_nH;

  //
  // First get the temperature.  We assume the total particle number density
  // is given by 1.1*nH*(1+x_in), appropriate for a gas with 10% Helium by 
  // number, and if He is singly ionised whenever H is.
  //
  double T = get_temperature(mpv_nH, E_in, x_in);


  double temp1=0.0, temp2=0.0;
  //
  // oneminusx_dot is in units of 1/s
  //
  double oneminusx_dot=0.0;
  //
  // Edot is calculated in units of erg/s per H nucleon, multiplied by
  // mpv_nH at the end.
  //
  double Edot=0.0;

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
  //
  // collisional excitation cooling of H0 Aggarwal (1983) and Raga+(1997,ApJS).
  //
  Edot -= Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10);

  //
  // Cosmic ray heating (HAdCM09 eq.A7).
  //
  //Edot += 5.0e-28*OneMinusX;

  //
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar neighbourhood.
  //
  //oneminusx_dot -= 1.8e-17*OneMinusX;
  
  //
  // now multiply Edot by nH to get units of energy loss/gain per unit volume per second.
  //
  Edot *= mpv_nH;

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





