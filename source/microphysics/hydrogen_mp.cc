/// \file hydrogen_mp.cc
/// \author Jonathan Mackey
/// \date 08.03.2011
///
/// This contains collected classes and functions for hydrogen microphysics.
///
/// Modifications:
/// - 2011.03.08 JM: written. 
/// - 2011.03.14 JM: Added photoionisation inheritance.  removed redundant PI x-section.
/// - 2011.03.29 JM: Added coll_ion_rates() function to return both ionisation- and
///    heating rates.
/// - 2011.04.14 JM: Fixed bugs.
/// - 2011.06.20 JM: Got rid of non-ANSI-C exp10 functions
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/hydrogen_mp.h"
#include "microphysics/hydrogen_recomb_Hummer94.h"
#include "global.h"
using namespace std;


Hydrogen_chem::Hydrogen_chem()
 : Hummer94_Hrecomb(), hydrogen_photoion()
{
  //cout <<"Setting up Hydrogen_chem().\n";

  cx_T = cx_rate = cx_rate2 = 0;
  cx_Nspl = 26;
  cx_T     = mem.myalloc(cx_T     , cx_Nspl);
  cx_rate  = mem.myalloc(cx_rate  , cx_Nspl);
  cx_rate2 = mem.myalloc(cx_rate2 , cx_Nspl);
  setup_Hi_coll_excitation_rate();
  return;
}

Hydrogen_chem::~Hydrogen_chem()
{
  cout <<"Deleting Hydrogen_chem() class.\n";
  cx_T     = mem.myfree(cx_T);
  cx_rate  = mem.myfree(cx_rate);
  cx_rate2 = mem.myfree(cx_rate2);
  return;
}

void Hydrogen_chem::setup_Hi_coll_excitation_rate()
{
  // -------------------------------------------------------------------
  // We need to set up a spline fit for the Aggarwal (1993) collisional excitation
  // of HI data.
  //
  //
  // Tabulated values for the cooling rate from Raga, Mellema, \& Lundqvist,
  // 1997, ApJS, 109, 517.
  // The rates are for collisionally excited cooling from neutral hydrogen (H I), 
  // and the rates are per electron, per H I atom, in units erg.cm3/s
  //
  double T[26] = {3162.2776602, 3981.0717055, 5011.8723363, 6309.5734448, 
                  7943.2823472, 10000.0000000, 12589.2541179, 15848.9319246,
                  19952.6231497, 25118.8643151, 31622.7766017, 39810.7170553, 
                  50118.7233627, 63095.7344480, 79432.8234724, 100000.0000000, 
                  125892.5411794, 158489.3192461, 199526.2314969, 251188.6431510, 
                  316227.7660168, 398107.1705535, 501187.2336273, 630957.3444802, 
                  794328.2347243, 1000000.0000000};
  double R[26] = {1.150800e-34, 2.312065e-31, 9.571941e-29, 1.132400e-26, 4.954502e-25, 
                  9.794900e-24, 1.035142e-22, 6.652732e-22, 2.870781e-21, 9.036495e-21, 
                  2.218196e-20, 4.456562e-20, 7.655966e-20, 1.158777e-19, 1.588547e-19, 
                  2.013724e-19, 2.393316e-19, 2.710192e-19, 2.944422e-19, 3.104560e-19, 
                  3.191538e-19, 3.213661e-19, 3.191538e-19, 3.126079e-19, 3.033891e-19, 
                  2.917427e-19};

  for (int i=0; i<cx_Nspl; i++) {
    cx_T[i] = log10(T[i]);
    cx_rate[i] = log10(R[i]);
  }
  GS.spline(cx_T, cx_rate, cx_Nspl, 1.e99, 1.e99, cx_rate2);
  //
  // Logarithmic slopes at either end of the domain.
  //
  cx_minT = cx_T[0];
  cx_maxT = cx_T[cx_Nspl-1];
  cx_MinSlope = (cx_rate[1]-cx_rate[0])/(cx_T[1]-cx_T[0]);
  cx_MaxSlope = (cx_rate[cx_Nspl-1]-cx_rate[cx_Nspl-2])/
                (cx_T[cx_Nspl-1]-cx_T[cx_Nspl-2]);
  //cout << "\t\tAlpha min-slope="<<cx_MinSlope<<" max-slope="<<cx_MaxSlope<<"\n";
  //
  // End of HI collisional excitation fit.
  // --------------------------------------------------------------------
  return;
}

double Hydrogen_chem::Hi_coll_excitation_cooling_rate(double T)
{
  //
  // Spline is fit in log-log space, and the slopes off the end of the fit
  // are also logarithmic, so we take the log of T, get log10(rate), and then
  // return exp10() of the rate.
  //
  double rate = 0.0;
  T = log10(T);

  if (T<cx_minT) {
    // y = y0 +m*(x-x0)
    rate = cx_rate[0]         + cx_MinSlope*(T - cx_minT);
  }
  else if (T>cx_maxT) {
    // y = y0 +m*(x-x0)
    rate = cx_rate[cx_Nspl-1] + cx_MaxSlope*(T - cx_maxT);
  }
  else {
    GS.splint(cx_T, cx_rate, cx_rate2, cx_Nspl, T, &rate);
  }
  return exp(2.302585093*rate);
}




double Hydrogen_chem::Hi_coll_ion_rate(double T)
{
  //
  // Voronov (1997) fit to collisional ionisation rate (cm3/s).
  //
  if (T < 5.0e3) {
    return 0.0;
  }
  else {
    //double A=2.91e-8, X=0.232, K=0.39, IP=13.6*1.602e-12; // ionisation potential.
    //T = ion_pot/kB/T;
    //return A*exp(K*log(T) -T)/(X+T);
    //
    // Here we avoid allocating new variables, and just hard-code the equation 
    // with its parameters.
    //
    T=1.578e5/T;  // Inverse of temperature in units of the ionisation potential.
    return 2.91e-8*exp(0.39*log(T)-T)/(0.232+T);
  }
}

double Hydrogen_chem::Hi_coll_ion_cooling_rate(double T)
{
  //
  // This is just the collisional ionisation rate times the ionisation potential,
  // which for hydrogen is 13.6*1.602e-12 ergs.  Answer is returned in erg.cm3/s
  //
  T=1.578e5/T;
  return 6.34e-19*exp(0.39*log(T)-T)/(0.232+T);
}

void Hydrogen_chem::Hi_coll_ion_rates(double T, double *cir, double *cicr)
{
  //
  // see Hi_coll_ion[_cooling]_rate(double T) for more details.
  // Answer is returned in erg.cm3/s
  //
  T=1.578e5/T; // Inverse of temperature in units of the ionisation potential.
  *cir  = 2.91e-8*exp(0.39*log(T)-T)/(0.232+T);
  *cicr = 2.18e-11*(*cir);
  return;
}



double Hydrogen_chem::Hi_monochromatic_photo_ion_heating(const double E)
{
  //
  // For monochromatic photoionisation, the heating is just (E-E0) for every
  // ionisation.  For H, E0 = 13.6*1.602e-12 ergs = 2.18e-11 ergs.
  //
  return (E - 2.18e-11);
}


