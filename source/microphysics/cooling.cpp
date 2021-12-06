/// \file cooling.cc
///
/// \brief Class Definitions for Cooling Function Routines.
/// \author Jonathan Mackey
///
/// Modifications:
///  - 2009-01-21 Moved Cooling out from microphysics.cc to its own file.
///  - 2009-04-28 Modified interface so function returns volume density cooling
///  rate.
///  - 2010-01-05 JM: Added in low-ionisation-only cooling hack IFDEF. (deffed
///  in global.h)
///  - 2010-01-12 JM: modified cooling==5 so that it is appropriate for the
///  Iliev et al. tests.
///  - 2010-01-15 JM: comments added.
///  - 2010-06-24 JM: Corrected KI02 function (with typos fixed from Vaz.-Sem.
///  et al 2007.
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"


#include <spdlog/spdlog.h>

#ifndef NDEBUG
#include "grid/cell_interface.h"
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "microphysics/cooling.h"

using namespace std;

#define GammaMinusOne 0.667

// ##################################################################
// ##################################################################

CoolingFn::CoolingFn(int flag)
{
  kB     = pconst.kB();
  Temp   = 0;
  Lamb   = 0;
  int id = 0;

#ifdef FUNCTION_ID
  string fname = "CoolingFn::CoolingFn";
  spdlog::info("\t\tCooling Function Constructor)"
#endif  // FUNCTION_ID
  if (flag == 1) {
    WhichFunction = 1;
    spdlog::debug(
        "\t\tFlag = {} corresponding to Sutherland & Dopita (1993, ApJS, 88, 253) cooling function\n\t\tfrom file pk6ff75.neq from http://www.mso.anu.edu.au/~ralph/data/cool/\n",
        flag);
    Temp = 0;
    Lamb = 0;
    Nspl = 69;

    Temp = mem.myalloc(Temp, Nspl);
    Lamb = mem.myalloc(Lamb, Nspl);

    double temp1[69] = {
        7.500, 7.445, 7.396, 7.347, 7.297, 7.246, 7.195, 7.143, 7.091, 7.039,
        6.987, 6.934, 6.882, 6.829, 6.776, 6.723, 6.670, 6.617, 6.563, 6.510,
        6.457, 6.403, 6.350, 6.296, 6.243, 6.189, 6.136, 6.082, 6.028, 5.975,
        5.921, 5.868, 5.814, 5.760, 5.707, 5.654, 5.600, 5.547, 5.493, 5.440,
        5.389, 5.336, 5.283, 5.230, 5.176, 5.122, 5.069, 5.015, 4.961, 4.907,
        4.854, 4.800, 4.746, 4.694, 4.647, 4.604, 4.546, 4.495, 4.450, 4.409,
        4.361, 4.318, 4.269, 4.225, 4.185, 4.141, 4.099, 4.053, 4.007};
    double temp2[69] = {
        -22.73, -22.74, -22.74, -22.74, -22.73, -22.72, -22.70, -22.67, -22.63,
        -22.60, -22.58, -22.59, -22.60, -22.62, -22.62, -22.61, -22.58, -22.56,
        -22.53, -22.51, -22.46, -22.39, -22.27, -22.12, -22.00, -21.96, -21.96,
        -21.98, -21.99, -21.99, -21.99, -21.96, -21.91, -21.86, -21.82, -21.79,
        -21.76, -21.70, -21.61, -21.53, -21.47, -21.43, -21.41, -21.40, -21.41,
        -21.43, -21.46, -21.50, -21.54, -21.59, -21.65, -21.70, -21.77, -21.83,
        -21.89, -21.94, -22.00, -22.04, -22.09, -22.14, -22.21, -22.29, -22.40,
        -22.51, -22.62, -22.74, -22.86, -22.98, -23.08};
    for (int i = 0; i < Nspl; i++) {
      Temp[i] = exp(log(10.0) * temp1[Nspl - 1 - i]);
      Lamb[i] = exp(log(10.0) * temp2[Nspl - 1 - i]);
    }
    MinTemp = Temp[0];
    MaxTemp = Temp[Nspl - 1];
#ifdef LINSLP
    MinSlope = (Lamb[1] - Lamb[0]) / (Temp[1] - Temp[0]);
    spdlog::debug("\t\tMinSlope (linear) = {}", MinSlope);
#endif
#ifdef LOGSLP
    MinSlope = (temp2[67] - temp2[68]) / (temp1[67] - temp1[68]);
    //    MinSlope *=2.0;
    MinSlope = 4.0;
    spdlog::debug("\t\tMinSlope (logarithmic) = {}", MinSlope);
#endif
    // The two large values passed to spline tell it to use zero second
    // derivative boundary conditions for extrapolation beyond the range of
    // the data.  It is dangerous to go beyond the range, but this boundary
    // condition means that extrapolation has a chance of being reasonable.
    spline(Temp, Lamb, Nspl, 0.0, 0.0, id);
  }  // SD93 function (NEQ)

  else if (flag == 2) {
    WhichFunction = 2;
    Temp          = 0;
    Lamb          = 0;
    ;
    spdlog::debug(
        "\t\tFlag = {} corresponding to Koyami & Inutsuka (2002, ApJL, 564, L97) cooling function\n\t\tThis is a double exponential fitting function, using their equations 4 and 5.\n\t\tN.B. The KI02 equation had two typos: (1.148e5 instead of 1.184e5, and 14 instead of 1.4e-2)\n\t\tThese were corrected by Vasquez-Semadeni et al. 2007, ApJ, 657, 870.\n",
        flag);
    MinTemp = 5.0;
    MaxTemp = 5.0e8;
  }  // KI02 function

  else if (flag == 3) {
    WhichFunction = 3;
    Temp          = 0;
    Lamb          = 0;
    spdlog::debug(
        "\t\tFlag = {} corresponding to Dalgarno & McCray (1972, ARAA, 10,375) cooling function\n",
        flag);
    Nspl = 31;

    Temp = mem.myalloc(Temp, Nspl);
    Lamb = mem.myalloc(Lamb, Nspl);

    double temp1[31] = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
                        3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2,
                        5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0};
    double temp2[31] = {-27.5,  -27.3,  -27.1, -26.9, -26.7, -26.5,  -26.3,
                        -26.1,  -25.95, -25.8, -25.7, -25.6, -25.55, -25.5,
                        -25.45, -25.4,  -21.7, -21.5, -21.4, -21.3,  -21.2,
                        -21.3,  -21.4,  -21.5, -21.7, -22.3, -22.6,  -22.7,
                        -22.7,  -22.7,  -22.7};
    for (int i = 0; i < Nspl; i++) {
      Temp[i] = temp1[i];  // exp(log(10.0)*temp1[i]);
      Lamb[i] = temp2[i];  // exp(log(10.0)*temp2[i]);
    }
    MinTemp = exp(log(10.0) * Temp[0]);
    MaxTemp = exp(log(10.0) * Temp[Nspl - 1]);
#ifdef LINSLP
    MinSlope = (Lamb[1] - Lamb[0]) / (Temp[1] - Temp[0]);
    spdlog::debug("\t\tMinSlope (linear) = {}", MinSlope);
#endif
#ifdef LOGSLP
    // MinSlope = (temp2[67]-temp2[68])/(temp1[67]-temp1[68]);
    //    MinSlope *=2.0;
    MinSlope = 1.0;
    spdlog::debug("\t\tMinSlope (logarithmic) = {}", MinSlope);
#endif
    // The two large values passed to spline tell it to use zero second
    // derivative boundary conditions for extrapolation beyond the range of
    // the data.  It is dangerous to go beyond the range, but this boundary
    // condition means that extrapolation has a chance of being reasonable.
    spline(Temp, Lamb, Nspl, 0.0, 0.0, id);
    spdlog::error(
        "{}: {}",
        "Dalgarno and McCray cooling function is not usably coded -- get "
        "a better function",
        999);
  }  // DMcC72 function

  else if (flag == 4) {
    WhichFunction = 4;
    spdlog::debug(
        "\t\tFlag = {} corresponding to Sutherland & Dopita (1993, ApJS, 88, 253) CIE cooling function\n\t\tfrom file m-00.cie from http://www.mso.anu.edu.au/~ralph/data/cool/\n",
        flag);
    Temp = 0;
    Lamb = 0;
    Nspl = 91;

    Temp = mem.myalloc(Temp, Nspl);
    Lamb = mem.myalloc(Lamb, Nspl);

    double temp1[91] = {
        4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30, 4.35, 4.40, 4.45, 4.50, 4.55,
        4.60, 4.65, 4.70, 4.75, 4.80, 4.85, 4.90, 4.95, 5.00, 5.05, 5.10, 5.15,
        5.20, 5.25, 5.30, 5.35, 5.40, 5.45, 5.50, 5.55, 5.60, 5.65, 5.70, 5.75,
        5.80, 5.85, 5.90, 5.95, 6.00, 6.05, 6.10, 6.15, 6.20, 6.25, 6.30, 6.35,
        6.40, 6.45, 6.50, 6.55, 6.60, 6.65, 6.70, 6.75, 6.80, 6.85, 6.90, 6.95,
        7.00, 7.05, 7.10, 7.15, 7.20, 7.25, 7.30, 7.35, 7.40, 7.45, 7.50, 7.55,
        7.60, 7.65, 7.70, 7.75, 7.80, 7.85, 7.90, 7.95, 8.00, 8.05, 8.10, 8.15,
        8.20, 8.25, 8.30, 8.35, 8.40, 8.45, 8.50};
    double temp2[91] = {
        -23.06, -22.46, -22.17, -21.92, -21.79, -21.80, -21.86, -21.90, -21.88,
        -21.82, -21.73, -21.63, -21.53, -21.42, -21.32, -21.22, -21.14, -21.07,
        -21.01, -20.98, -20.99, -21.02, -21.03, -21.01, -20.98, -20.97, -20.96,
        -20.96, -20.99, -21.13, -21.35, -21.55, -21.66, -21.71, -21.71, -21.71,
        -21.76, -21.86, -21.93, -21.95, -21.96, -21.96, -21.96, -21.95, -21.94,
        -21.97, -22.07, -22.20, -22.31, -22.39, -22.44, -22.48, -22.50, -22.53,
        -22.56, -22.59, -22.60, -22.60, -22.59, -22.57, -22.57, -22.59, -22.62,
        -22.65, -22.68, -22.70, -22.72, -22.73, -22.73, -22.73, -22.73, -22.72,
        -22.71, -22.70, -22.68, -22.67, -22.65, -22.64, -22.62, -22.60, -22.58,
        -22.56, -22.54, -22.53, -22.51, -22.49, -22.47, -22.45, -22.43, -22.40,
        -22.38};
    for (int i = 0; i < Nspl; i++) {
      Temp[i] = exp(log(10.0) * temp1[i]);
      Lamb[i] = exp(log(10.0) * temp2[i]);
    }
    MinTemp = Temp[0];
    MaxTemp = Temp[Nspl - 1];
#ifdef LINSLP
    MinSlope = (Lamb[1] - Lamb[0]) / (Temp[1] - Temp[0]);
    spdlog::debug("\t\tMinSlope (linear) = {}", MinSlope);
#endif
#ifdef LOGSLP
    MinSlope = (temp2[1] - temp2[0]) / (temp1[1] - temp1[0]);
    spdlog::debug("\t\tMinSlope (logarithmic) = {}", MinSlope);
#endif
    spline(Temp, Lamb, Nspl, 0.0, 0.0, id);
  }  // SD93 function (CIE)

  else if (flag == 5) {
    // C2Ray -- free-free cooling and collisional excitation cooling due to
    // HI
    spdlog::info(
        "\t\tC2RAY-like cooling due to free-free and HI coll.excitation.\t\t the recombination and ionisation cooling are in the microphysics routine.\n\t\t I approx the ff cooling with 1.34e-11*sqrt(T_e)*ne*ni, good for 10^3<T<10^4 or so.\n\t\t I ignore inverse compton cooling off the cmb, since i'm at z=0.\n\t\t I use 7.5e-19*exp(-118348/T)/(1+sqrt(T/1.e5))*n_e*n_H0 as the c.ex. cooling rate...\n");
    MinTemp       = 0.000001;
    MaxTemp       = 1.e8;
    WhichFunction = 5;
    Temp          = 0;
    Lamb          = 0;
  }  // c2ray

  else if (flag == 6 || flag == 7 || flag == 8 || flag == 9 || flag == 10) {
    // cout <<"\t\tPower law cooling function";
    MinTemp       = 0.000001;
    MaxTemp       = 1.e8;
    WhichFunction = flag;

    switch (WhichFunction) {
      case 6:
        MinSlope = 0.5;
        break;
      case 7:
        MinSlope = 1.0;
        break;
      case 8:
        MinSlope = 2.5;
        break;
      case 9:
        MinSlope = 5.0;
        break;
      case 10:
        MinSlope = 10.0;
        break;
      default:
        spdlog::error("{}: {}", "bad cooling flag", WhichFunction);
    }
    // cout <<" with log slope="<<MinSlope<<"\n\t\t and pivot at 8000K";
    // cout <<" with rate 1e-24 at 8000K.\n";
  }

  else if (
      flag == 11 || flag == 12 || flag == 13 || flag == 14 || flag == 15
      || flag == 16) {
    spdlog::debug(
        "\t\tSD93-CIE function with power law contribution to account for\t\tforbidden line cooling in ionised gas.\n");
    if (flag == 12)
      spdlog::debug(
          "\t\tflag=={}: Plus Exponential Cooling in neutral gas (TOY MODEL n_H*t_c=1e6yrs)\n",
          flag);
    if (flag == 13)
      spdlog::debug(
          "\t\tflag=={}: Plus Exponential Cooling in neutral gas (TOY MODEL n_H*t_c=1e7yrs)\n",
          flag);
    if (flag == 14)
      spdlog::debug(
          "\t\tflag=={}: Plus Exponential Cooling in neutral gas (TOY MODEL t_c=1e3yrs)\n",
          flag);
    if (flag == 15)
      spdlog::debug(
          "\t\tflag=={}: Plus Exponential Cooling in neutral gas (TOY MODEL t_c=1e4yrs)\n",
          flag);
    if (flag == 16)
      spdlog::debug(
          "\t\tflag=={}: Plus Henney et al (2009) cooling terms.\n", flag);
    WhichFunction = flag;
    Temp          = 0;
    Lamb          = 0;
    Nspl          = 91;

    Temp = mem.myalloc(Temp, Nspl);
    Lamb = mem.myalloc(Lamb, Nspl);

    double temp1[91] = {
        4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30, 4.35, 4.40, 4.45, 4.50, 4.55,
        4.60, 4.65, 4.70, 4.75, 4.80, 4.85, 4.90, 4.95, 5.00, 5.05, 5.10, 5.15,
        5.20, 5.25, 5.30, 5.35, 5.40, 5.45, 5.50, 5.55, 5.60, 5.65, 5.70, 5.75,
        5.80, 5.85, 5.90, 5.95, 6.00, 6.05, 6.10, 6.15, 6.20, 6.25, 6.30, 6.35,
        6.40, 6.45, 6.50, 6.55, 6.60, 6.65, 6.70, 6.75, 6.80, 6.85, 6.90, 6.95,
        7.00, 7.05, 7.10, 7.15, 7.20, 7.25, 7.30, 7.35, 7.40, 7.45, 7.50, 7.55,
        7.60, 7.65, 7.70, 7.75, 7.80, 7.85, 7.90, 7.95, 8.00, 8.05, 8.10, 8.15,
        8.20, 8.25, 8.30, 8.35, 8.40, 8.45, 8.50};
    double temp2[91] = {
        -23.06, -22.46, -22.17, -21.92, -21.79, -21.80, -21.86, -21.90, -21.88,
        -21.82, -21.73, -21.63, -21.53, -21.42, -21.32, -21.22, -21.14, -21.07,
        -21.01, -20.98, -20.99, -21.02, -21.03, -21.01, -20.98, -20.97, -20.96,
        -20.96, -20.99, -21.13, -21.35, -21.55, -21.66, -21.71, -21.71, -21.71,
        -21.76, -21.86, -21.93, -21.95, -21.96, -21.96, -21.96, -21.95, -21.94,
        -21.97, -22.07, -22.20, -22.31, -22.39, -22.44, -22.48, -22.50, -22.53,
        -22.56, -22.59, -22.60, -22.60, -22.59, -22.57, -22.57, -22.59, -22.62,
        -22.65, -22.68, -22.70, -22.72, -22.73, -22.73, -22.73, -22.73, -22.72,
        -22.71, -22.70, -22.68, -22.67, -22.65, -22.64, -22.62, -22.60, -22.58,
        -22.56, -22.54, -22.53, -22.51, -22.49, -22.47, -22.45, -22.43, -22.40,
        -22.38};
    for (int i = 0; i < Nspl; i++) {
      Temp[i] = exp(log(10.0) * temp1[i]);
      Lamb[i] = exp(log(10.0) * temp2[i]);
    }
    MinTemp  = Temp[0];
    MaxTemp  = Temp[Nspl - 1];
    MinSlope = (temp2[1] - temp2[0]) / (temp1[1] - temp1[0]);
    spdlog::debug("\t\tMinSlope (logarithmic) = {}", MinSlope);
    spline(Temp, Lamb, Nspl, 0.0, 0.0, id);
  }  // SD93-CIE-ForbiddenLine

  else
    spdlog::error("{}: {}", "Bad flag in CoolingFn Constructor", flag);

  CoolingFn::spline_id = id;

#ifdef COOL_TESTING
  ofstream outf("coolingcurve.txt");
  if (!outf.is_open()) spdlog::error("{}: {}", "couldn't open outfile", 1);
  outf << "Cooling Curve Data: Temperature(K) Rate(erg/cm^3/s) (n=1 per cc) "
          "ifracs: 1e-6, 1e-3, 0.1, 0.5, 0.99, 0.999999, and then the same for "
          "n=1e6/cc\n";
  outf.setf(ios_base::scientific);
  outf.precision(6);
  double t = 1.e0;
  do {
    outf << t << "\t" << fabs(CoolingRate(t, 1.e-6, 1.0, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.001, 1.0, 0, 0));
    outf << "\t" << fabs(CoolingRate(t, 0.1, 1.0, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.5, 1.0, 0, 0));
    outf << "\t" << fabs(CoolingRate(t, 0.99, 1.0, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.999999, 1.0, 0, 0));  //<<"\n";

    outf << "\t" << fabs(CoolingRate(t, 1.e-6, 1.e6, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.001, 1.e6, 0, 0));
    outf << "\t" << fabs(CoolingRate(t, 0.1, 1.e6, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.5, 1.e6, 0, 0));
    outf << "\t" << fabs(CoolingRate(t, 0.99, 1.e6, 0, 0)) << "\t"
         << fabs(CoolingRate(t, 0.999999, 1.e6, 0, 0)) << "\n";
    t *= 1.05;
  } while (t < 1.e7);
  outf.close();
#endif  // COOL_TESTING
  // spdlog::error("{}: {}", "quitting",100);
  return;
}

// ##################################################################
// ##################################################################

CoolingFn::~CoolingFn()
{
  if (Temp) Temp = mem.myfree(Temp);
  if (Lamb) Lamb = mem.myfree(Lamb);
}

// ##################################################################
// ##################################################################

//
// Returns Cooling rate in erg/cm3/s i.e. volumetric rate of energy loss.
//
double CoolingFn::CoolingRate(
    const double T,
    const double xHp,
    const double nH,
    const double FUV_flux,
    const double FUV_extinction)
{
  if (xHp > 1.0 || T <= 0.0 || !isfinite(T)) {
    // cout <<"input ion fraction >1 in CoolingRate()... x="<<xHp<<"\n";
    return 0.0;
  }

  double rate = 0.0;
  //
  // All these rates were originally coded for returning erg.cm^3/s, so I need
  // to Multiply them by n_e.n_i = (xHp.nH)^2 to get the volume density rate.
  //
  if (WhichFunction == 1 || WhichFunction == 4) {
    if (T > MaxTemp) {
      spdlog::debug(
          "Temp out of range!! Too large: T={} and MAX.T={}", T, MaxTemp);
#ifndef NDEBUG
      CI.print_cell(dp.c);
#endif
      spdlog::debug("Returning Lambda(MaxTemp) = Lambda({}", MaxTemp);
      splint(Temp, Lamb, spline_id, Nspl, MaxTemp, &rate);
    }
    else if (T <= 0.)
      rate = 0.0;
    else if (T < MinTemp) {
      //      cout <<"Temp too small... returning smooth value.\n";
#ifdef LINSLP
      //      cout <<"T="<<T<<", MinTemp="<<MinTemp<<", rate = "<<max(0.0,
      //      Lamb[0]+MinSlope*(T-MinTemp))<<"\n";
      rate = max(0.0, Lamb[0] + MinSlope * (T - MinTemp));
#endif
#ifdef LOGSLP
      //      cout <<"T="<<T<<", MinTemp="<<MinTemp<<", fraction =
      //      "<<exp(MinSlope*(log(T)-log(MinTemp)))<<"\n";
      rate = max(1.e-50, Lamb[0] * exp(MinSlope * (log(T) - log(MinTemp))));
#endif
    }
    else {
      splint(Temp, Lamb, spline_id, Nspl, T, &rate);
    }
    //
    // Multiply by n_e.n_i
    //
    rate *= xHp * xHp * nH * nH;
  }  // SD93 funtion (NEQ and CIE)

  else if (WhichFunction == 2) {
    if (T > MaxTemp) {
      spdlog::debug(
          "Warning: very large temperature!: T={} and MAX.T={}", T, MaxTemp);
    }
    //
    // Koyama and Inutsuka: 2002, ApJL, 564, 97 (KI02)
    // KI02 cooling has a cooling rate propto n^2 and a heating rate propto
    // n. It seems like this rate is independent of the ion fraction, so
    // they must assume a transition temperature between ionised and
    // neutral, and this must be implicitly included in their cooling rate
    // function.
    //
    // First the cooling rate:  (KI02, eq.4) multiplied by n^2 to get volume
    // rate (see eq.3), and corrected by V-S etal. 2007, ApJ, 657, 870.
    // For n<10^6 per c.c.,  heating dominates for T>5K, so MinTemp is set
    // to 5K. If T<MinTemp there is no point calculating the cooling, so we
    // skip it.
    //
    if (T > MinTemp)
      rate += nH * nH
              * (2.0e-19 * exp(-1.184e5 / (T + 1.0e3))
                 + 2.8e-28 * sqrt(T) * exp(-92.0 / T));
    //
    // Now the heating rate:
    //
    rate -= nH * 2.0e-26;
#ifndef NDEBUG
    spdlog::debug("KI02 cooling: T={}, n={}, rate={}", T, nH, rate);
#endif
  }  // KI02 function.

  else if (WhichFunction == 3) {
    spdlog::info(
        "DALGARNO AND MCCRAY FUNCTION IS VERY DODGY -- BAD INTERPOLATION");
    if (T > MaxTemp) {
      spdlog::debug(
          "Temp out of range!! Too large: T={} and MAX.T={}", T, MaxTemp);
#ifndef NDEBUG
      CI.print_cell(dp.c);
#endif
      spdlog::debug("Returning Lambda(MaxTemp) = Lambda({}", MaxTemp);
      splint(Temp, Lamb, spline_id, Nspl, log10(MaxTemp), &rate);
      rate = exp(log(10.0) * rate);
    }
    else if (T <= 0.)
      rate = 0.0;
    else if (T < MinTemp) {
      //      cout <<"Temp too small... returning smooth value.\n";
      rate = 0.0;
#ifdef LINSLP
      rate = max(0.0, Lamb[0] + MinSlope * (T - MinTemp));
#endif
#ifdef LOGSLP
      rate = max(1.e-50, Lamb[0] * exp(MinSlope * (log(T) - log(MinTemp))));
#endif
    }
    else {
      splint(Temp, Lamb, spline_id, Nspl, log10(T), &rate);
      rate = exp(log(10.0) * rate);
    }
    //
    // Multiply by n_e.n_i
    //
    rate *= xHp * xHp * nH * nH;
  }  // DMcC72 funtion

  else if (WhichFunction == 5) {
    //
    //  THIS IS NOW PURE HYDROGEN COOLING, FOR ILIEV ET AL TESTS!!!!
    //
    // I can add Osterbrock's (1989) function for the free-free H+
    // component here, as Hummer says it is good to 30% in the worst
    // case.
    //
    rate += 1.34e-11 * kB * sqrt(T) * nH * nH * xHp * xHp;

    //
    // collisionally excited neutral H cooling from Raga, Mellema,
    // Lundqvist, 1997, ApJS, 109, 517, appendix. I fit it with gnuplot, BY
    // EYE!!!, SO THIS IS NOT FOR PRODUCTION!!! log10(cool) =
    // -18.5*(1+exp(-(log10(T)-3.44)/0.4)) in erg.cm3/s
    //
    //    rate +=
    //    exp(-18.5*(1.0+exp(-2.5*(log10(T)-3.44)))*log(10.0))*nH*nH*xHp*(1.0-xHp);
    //
    // This is from somewhere else, and is almost the same...
    //
    rate += 7.5e-19 * exp(-118348.0 / T) / (1.0 + sqrt(T / 1.0e5)) * nH * nH
            * xHp * (1.0 - xHp);
    //
    // And the recombination cooling is already taken care of in
    // microphysics
    //
    return rate;

    //
    // OLD C2RAY GUESSES.
    //
    // I think this is C2Ray-like cooling, appropriate for pure
    // hydrogen/helium gas.
    //
    // first free-free for hydrogen only.
    //    rate += kB*sqrt(T)*1.25e-11; // rate erg.cm^3/s is
    //    kT*1.25e-11/sqrt(T)
    // cout <<"\t\tT="<<T<<", ff rate="<<rate;
    // now excitation cooling for hydrogen only.
    //    if (xHp<0.1*SMALLVALUE || xHp>(1.0+SMALLVALUE)){
    //      spdlog::error("{}: {}", "ion fraction too small in
    //      cool->CoolingRate()",xHp);
    // THIS MUST BE A CRAZY INTEGRATION ATTEMPT, SO RETURN A MASSIVE COOLING
    // RATE...
    //      rate += 1.0/TINYVALUE;
    //    }
    // double temp
    // = 7.5e-19*exp(-118348.0/T)/(1.0+sqrt(T/1.0e5))*(1.0-xHp)/xHp; cout
    // <<", c.ex. rate="<<temp; rate += temp;
    //    rate
    //    += 7.5e-19*exp(-118348.0/T)/(1.0+sqrt(T/1.0e5))*(1.0-xHp)/xHp;
    // cout <<", ::: total="<<rate<<"\n";

    //
    // Multiply by n_e.n_i
    //
    //    rate *= xHp*xHp*nH*nH;
  }

  else if (
      WhichFunction == 6 || WhichFunction == 7 || WhichFunction == 8
      || WhichFunction == 9 || WhichFunction == 10) {
    // power law cooling.
    rate = 1.0e-24 * exp(MinSlope * log(T / 8000.0));
    //
    // Multiply by n_e.n_i
    //
    rate *= xHp * xHp * nH * nH;
  }  // power law cooling.

  else if (WhichFunction == 11) {
    //
    // This corresponds to SD93 CIE cooling function, with an extra
    // component for cooling due to collisionally excited ions of C,O in
    // photo-ionised plasma, from Osterbrock, or Spitzer.
    //
    if (T > MaxTemp) {
      // cout <<"Temp out of range!! Too large: T="<<T<<" and
      // MAX.T="<<MaxTemp<<"\n";
#ifndef NDEBUG
      // CI.print_cell(dp.c);
#endif
      // cout <<"Returning Lambda(MaxTemp) = Lambda("<<MaxTemp<<")\n";
      splint(Temp, Lamb, spline_id, Nspl, MaxTemp, &rate);
    }
    else if (T <= 0.)
      rate = 0.0;
    else if (T < MinTemp) {
      // cout <<"T="<<T<<", MinTemp="<<MinTemp<<", fraction =
      // "<<exp(MinSlope*(log(T)-log(MinTemp)))<<"\n";
      rate = max(1.e-50, Lamb[0] * exp(25.0 * (log(T) - log(MinTemp))));
    }
    else {
      splint(Temp, Lamb, spline_id, Nspl, T, &rate);
    }
    if (T < 2.0e4)
      rate += xHp * 2.0e-24 * T / 8000.0;  // exp(2.0*log(T/8000.0));

    //
    // Multiply by n_e.n_i
    //
    rate *= xHp * xHp * nH * nH;

    //
    // If the dynamics gives a negative pressure, we may need to ensure some
    // heating happens... (BUT THIS IS DANGEROUS -- CAN MAKE CLUMPS EXPLODE,
    // SO USE WITH EXTREME CAUTION!!!!!) double t_c=3.16e11, T_inf=50.0; if
    // (T<T_inf) {
    //  //cout <<"\t\tcooling_rate(): artificially heating the
    //  gas...T="<<T<<"K, and we want T>=100K.\n";
    //  //cout <<"\t\tcooling_rate(): converting rate from "<<rate<<" to ";
    //  rate += nH*nH*kB*(T-T_inf)/(GammaMinusOne)/t_c;
    //  //cout <<rate<<"\n";
    //}
  }  // SD93-CIE + Forbidden line

  else if (WhichFunction >= 12 && WhichFunction <= 15) {
    //
    // This corresponds to SD93 CIE cooling function, with an extra
    // component for cooling due to collisionally excited ions of C,O in
    // photo-ionised plasma, from Osterbrock, or Spitzer. Plus a toy model
    // for cooling neutral dense gas to 100K.
    //

    //
    // First the SD CIE Model:
    //
    if (T > MaxTemp) {
      // cout <<"Temp out of range!! Too large: T="<<T<<" and
      // MAX.T="<<MaxTemp<<"\n";
#ifndef NDEBUG
      // CI.print_cell(dp.c);
#endif
      // cout <<"Returning Lambda(MaxTemp) = Lambda("<<MaxTemp<<")\n";
      splint(Temp, Lamb, spline_id, Nspl, MaxTemp, &rate);
    }
    else if (T <= 0.)
      rate = 0.0;
    else if (T < MinTemp) {
      // cout <<"T="<<T<<", MinTemp="<<MinTemp<<", fraction =
      // "<<exp(25.0*(log(T)-log(MinTemp)))<<"\n";
      rate = max(1.e-50, Lamb[0] * exp(25.0 * (log(T) - log(MinTemp))));
    }
    else {
      splint(Temp, Lamb, spline_id, Nspl, T, &rate);
    }

    //
    // Now add the Line-Cooling term:
    // Add in the minimum Temperature so that dodgy calls won't get strange
    // results.
    //
    if (T < 2.0e4 && T > 100.0)
      rate += 2.0e-24 * T / 8000.0;  // exp(2.0*log(T/8000.0));

    //
    // Multiply by n_e.n_i
    //
    rate *= xHp * xHp * nH * nH;

    //
    // Toy model -- Newton's law of cooling, with time constant given by
    // t_c = XX*10^6yrs/((1-x)^2 nH)
    // (nH in the denominator of t_c gives time constant shorter with
    // density) (Change t_c to get a different normalisation.) Have an
    // exponential cutoff at 1e4K so that hot gas doesn't have this cooling.
    //
    double t_c = 0.0, t_h = 0.0, T_inf = 100.0, T_min = 5.0;
    switch (WhichFunction) {
      case 12:
        t_c = 3.16e13 / nH;
        t_h = 3.16e11 / nH;
        break;  // 1e6 yrs/nH
      case 13:
        t_c = 3.16e14 / nH;
        t_h = 3.16e12 / nH;
        break;  // 1e7 yrs/nH
      case 14:
        t_c = 1.0e11;
        t_h = 3.16e9;
        break;  // ~3e3 yrs: Takes away dependency on density!
      case 15:
        t_c = 3.16e11;
        t_h = 3.16e9;
        break;  // 1e4 yrs: Takes away dependency on density!
      default:
        spdlog::error(
            "{}: {}", "Bad flag in toy model cooling!!!", WhichFunction);
    }

    //
    // THIS IS A BIG HACK -- IF I GET A NEGATIVE PRESSURE, I START WITH A
    // VERY LOW TEMPERATURE, AND THE FOLLOWING LINE MEANS LOW TEMPERATURE
    // CELLS HAVE VERY SHORT ENERGY GAIN TIMESCALES, SO THERE IS **VERY**
    // STRONG HEATING IN THE MICROPHYSICS
    //
    if (T < T_min) {
      //
      // means large energy gain for very low temperature cells
      // (counteract neg.press.).
      //
      rate += nH * kB * (T - T_min) / (GammaMinusOne) / t_h;
    }
    else {
      //
      // "Standard" exponential cooling, where only cool neutral gas is
      // affected. The max(0,R) means that we don't get heating if
      // T_min<T<T_inf The exponential part of the function has no effect
      // in this case.
      //
      rate +=
          max(0.0, (1.0 + xHp) * (1.0 - xHp) * (1.0 - xHp) * nH * kB
                       * (T - T_inf) / (GammaMinusOne) / t_c * exp(-T / 1.e4));
    }
  }  // SD93-CIE + Forbidden line [+toy-model Molecular cooling.]

  else if (WhichFunction == 16) {
    //
    // This function uses the cooling components from the appendix of
    // Henney et al. (2009) for cooling by excited neutral metals,
    // cooling by ionised excited metals, and cooling by neutral and
    // molecular dense gas.
    //

    double z0 = 5.0e-4;  // ROUGH APPROXIMATION TO ISM OXYGEN ABUNDANCE
    //
    // Now the ionised metal line cooling (Henney et al. 2009, eq. A9)
    //
    rate += 2.905e-19 * z0 * xHp * xHp * nH * nH
            * exp(-33610.0 / T - (2180.0 * 2180.0 / T / T));

    //
    // Now the neutral metal line cooling (Henney et al. 2009, eq. A10)
    //
    rate += 4.477e-20 * z0 * xHp * (1.0 - xHp) * nH * nH
            * exp(-28390.0 / T - (1780.0 * 1780.0 / T / T));

    //
    // The H+ recombination (Henney et al. 2009, eq. A11) is taken
    // care of in microphysics (but maybe it should be done here??? I
    // think that is a better idea for the future)
    //
    // I can add Osterbrock's (1989) function for the free-free H+
    // component here, as Hummer says it is good to 30% in the worst
    // case.
    //
    rate += 1.34e-11 * kB * sqrt(T) * nH * nH * xHp * xHp;

    //
    // I can't find the collisionally excited neutral H cooling
    // (Henney et al. 2009, eq. A12) anywhere (it's certainly not in
    // Hummer (1994)), but I think it is less important than that of
    // the ions, so I'll ignore it and hope for the best.
    //

    //
    // Now the CIE component (Henney et al. 2009, eq. A13)
    // I have fit a linear smoothing to this, spread over 20,000K,
    // so that there isn't such a sharp spike in the cooling rate.
    //
    if (T > 5.0e4) {
      double temp = 3.485e-15 * z0 * exp(-0.63 * log(T))
                    * (1.0 - exp(-pow(1.0e-5 * T, 1.63)));
      temp *= min(1.0, (T - 5.0e4) / (2.0e4));
      temp *= nH * nH * xHp * xHp;
      rate += temp;
    }

    //
    // Now the molecular cooling (Henney et al. 2009, eq. A14)
    //
    rate += (1.0 - xHp) * (1.0 - xHp) * 3.981e-27 * exp(1.6 * log(nH)) * sqrt(T)
            * exp(-(70.0 + 220.0 * exp(0.2 * log(nH * nH / 1.e12))) / T);
    //
    // HACK ALERT!!! I AM REDUCING THE DENSITY SCALING COMPARED TO WHAT
    // HENNEY USES!
    //
    // rate += (1.0-xHp)*(1.0-xHp)*3.981e-27*exp(1.3*log(nH))*sqrt(T)
    //  *exp(-(70.0+220.0*exp(0.2*log(nH*nH/1.e12)))/T);

    //
    // FUV heating (Henney et al. 2009, eq. A3)
    //
    rate -= 1.9e-26 * nH * FUV_flux * exp(-1.9 * FUV_extinction)
            / (1.0 + 6.4 * FUV_flux * exp(-1.9 * FUV_extinction) / nH);
    //
    // X-ray heating (Henney et al. 2009, eq. A5)
    // This makes no difference at all as far as I can tell...
    //
    //    rate -= 6.0e-23*nH*FUV_flux*8.0e-16;
    //
    // IR Heating (Henney et al. 2009, eq. A6)
    //
    rate -= 7.7e-32 * nH * FUV_flux * exp(-0.05 * FUV_extinction)
            * exp(-2.0 * log(1.0 + 3.0e4 / nH));
    //
    // Cosmic Ray Heating (Henney et al. 2009, eq. A7)
    // HACK ALERT!!! I AM INCREASING THIS BY 10X TO COMPENSATE FOR NO XRAY
    // HEATING!!!
    //
    rate -= 5.0e-27 * nH;

    //
    // We want to have a "soft landing" to the equilibrium neutral
    // temperature, so we could multiply the cooling rate by (T-Tinf)
    // Note if T<T_inf, this means that we are changing the sign of
    // the cooling rate, to be a net heating rate, so we set the rate
    // to zero in this case.
    // NB Cooling is positive, hence max() and not min().  The max() is
    // to avoid artificial heating in cold gas, which can cause explosions
    // if the initial conditions are too cold...
    //
    double T_inf  = 200.0;
    double T_soft = 300.0;
    if (T < T_soft && rate > 0.0) {
      rate = max(0.0, rate * (T - T_inf) / (T_soft - T_inf));
    }

  }  // SD93-CIE + Henney et al. (2009) model. (C16)

  else
    spdlog::error(
        "{}: {}", "Which Function? CoolingFn::CoolingRate", WhichFunction);

  // Now whichever function we use, it should have calculated 'rate', so we
  // can just return it.
  return rate;
}

// ##################################################################
// ##################################################################
