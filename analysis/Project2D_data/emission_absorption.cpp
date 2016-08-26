/// \file emission_absorption.cpp
/// \author Jonathan Mackey
/// \date 2016.07.04
///
/// Purpose:
/// - various functions for calculating emission/absorption
///   coefficients, and calculating projected quantities.
///


#include <iostream>
//#include <sstream>
#include <cmath>
//#include <vector>

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
//#include "tools/reporting.h"
//#include "tools/mem_manage.h"
//#include "tools/interpolate.h"
#include "constants.h"
#include "sim_params.h"


#include "emission_absorption.h"
#include "xray_emission.h"
using namespace std;



// ##################################################################
// ##################################################################



int get_emission_absorption_data(
      double const* const* data, ///< raw data to get variable from
      const int n_img,    ///< number of images to write
      const size_t Nr, ///< Number of radial data elements
      class Xray_emission &XR,  ///< pointer to X-ray emission class.
      double **ems,    ///< array for emission[img][rad] data.
      double **abs     ///< array for absorption[img][rad] data.
      )
{
  double xr[5]; xr[0]=0.0,xr[1]=0.0,xr[2]=0.0,xr[3]=0.0,xr[4]=0.0;
  //
  // 1/Mean mass per H atom:  n(H) = immph*rho
  //
  double immpH = SimPM.EP.H_MassFrac/pconst.m_p();
  //cout <<"immpH="<<immpH<<", mmpH="<<1.0/immpH<<"\n";
  //
  // square arcsec in radians:
  //
  double sq_arcsec = 1.0/(4.0*M_PI*206265*206265);
  //
  // variable for pre-factor
  //
  double prefactor = 0.0;

  //
  // Loop over all image variables.
  //
  for (int im=0; im<n_img; im++) {

    //
    // Put data for variable ivar into ems[im][] array, and variable's
    // slope with radius in vsl[] array.
    // Use 1st order forward differencing for all slopes, except for
    // the last point which uses backwards differencing.
    //
    switch (im) {
      //
      // Projected density is easy.
      //
      case PROJ_D:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i];
        if (i==Nr-1) {
          abs[im][i] = abs[im][i-1];
        }
        else {
          abs[im][i] = (data[DATA_D][i+1]-data[DATA_D][i])/
                    (data[DATA_R][i+1]-data[DATA_R][i]);
        }
      }
      break;

      //
      // Projected Neutral density: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      //
      case PROJ_NtD:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i]*(1.0-data[DATA_TR0][i]);
      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/(data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // Projected ionised density: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      //
      case PROJ_InD:
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = data[DATA_D][i]*data[DATA_TR0][i];
      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/(data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // Projected emission measure: here do two loops, one to set the
      // variable and one to get the slope (b/c the var is non-linear).
      // Here we hardcode conversion units.   We assume n_e= 1.1n_p, 
      // appropriate if He ionization follows H.
      // 
      case PROJ_EM:
      prefactor = immpH*immpH*1.21;
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = prefactor*
                      data[DATA_D][i]*data[DATA_TR0][i]*
                      data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;

      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>0.1 keV
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X01:
      prefactor = 1.1*immpH*immpH*sq_arcsec;
      //cout <<"X01 prefactor="<<prefactor<<"\n";
      for (size_t i=0; i<Nr; i++) {
        XR.get_xray_emissivity(data[DATA_T][i],xr);
        ems[im][i] = prefactor*xr[0]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>0.5 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X05:
      prefactor = 1.1*immpH*immpH*sq_arcsec;
      for (size_t i=0; i<Nr; i++) {
        XR.get_xray_emissivity(data[DATA_T][i],xr);
        ems[im][i] = prefactor*xr[1]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>1.0 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X10:
      prefactor = 1.1*immpH*immpH*sq_arcsec;
      for (size_t i=0; i<Nr; i++) {
        XR.get_xray_emissivity(data[DATA_T][i],xr);
        ems[im][i] = prefactor*xr[2]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      
      //
      // X-ray emission: ignore absorption, just calculate emmisivity
      // from the table, for E>5.0 keV.
      // n_e *n_H = 1.1 nH^2 = 1.1*(X_H*rho/m_p)^2 = 2.01e47*rho^2
      // emissivity per unit solid angle... divide by 4*pi*206265^2
      // so that we have intensity in erg/cm3/s/sq.arcsec
      //
      case PROJ_X50:
      prefactor = 1.1*immpH*immpH*sq_arcsec;
      for (size_t i=0; i<Nr; i++) {
        XR.get_xray_emissivity(data[DATA_T][i],xr);
        ems[im][i] = prefactor*xr[3]*
                     data[DATA_D][i]*data[DATA_TR0][i]*
                     data[DATA_D][i]*data[DATA_TR0][i];

      }
      for (size_t i=0; i<Nr-1; i++) {
        abs[im][i] = (ems[im][i+1]-ems[im][i])/
                     (data[DATA_R][i+1]-data[DATA_R][i]);
      }
      abs[im][Nr-1] = abs[im][Nr-2];  // backward diff is same as forward@(i-1)
      break;
      

      //
      // H-alpha emission:  here ems[im][] is the emissivity, and
      // abs[im][] the absorption coefficient.
      // Here we hardcode conversion units.  We assume n_e= 1.1n_p, 
      // appropriate if He ionization follows H.
      // We further assume X_H=0.715 (H mass fraction), similar to the
      // Asplund et al. (2009) value of 0.7154.
      //
      // Emissivity from Osterbrock j(Ha)=2.63e-33*n_e*n_p/T^0.9 in 
      // units of erg/cm3/s/sq.arcsec (adapted from table), and use
      //  n_e*n_p = rho^2 y^2 1.1(X_H/m_p)^2
      //
      // 5.28e14
      //
      // Absorption from Henney et al. 2009, where they assume 
      //  alpha = 5.0e-22 nH per cm (absorption by dust). ~213.7
      //
      case PROJ_HA:
      prefactor = immpH*immpH*1.1*2.63e-33;
      //cout <<"Halpha prefactor="<<prefactor<<", expecting ~5.28e14";
      //cout <<".  Abs="<<immpH*5.0e-22<<", expecting 213.7\n";
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = prefactor*data[DATA_D][i]*data[DATA_TR0][i]*
                          data[DATA_D][i]*data[DATA_TR0][i]*
                          pow(data[DATA_T][i], -0.9);
                          //exp(-0.9*log(data[DATA_T][i]));
      //if (!isfinite(ems[im][i])) cout <<"ems="<< ems[im][i]<<"\n";
#ifdef ABSORPTION
        abs[im][i] = immpH*5.0e-22 *data[DATA_D][i];
#else
        abs[im][i] = 0.0;
#endif
      }
      break;
      
      //
      // Ionised Metal-lines:  here ems[im][] is the emissivity, and
      // abs[im][] the absorption coefficient.
      // Here we hardcode conversion units.  We assume n_e= n_p, 
      // appropriate if He remains neutral.
      //
      // Emissivity from Dopita (1973)
      //  j([NII] ll 6584) =
      //   6.82e-18 n_e*n_p*f(N)*exp(-chi/kT)/(4*pi*sqrt(T))
      // in units of erg/cm3/s/sr, and use
      //  n_e*n_p = rho^2 y^2 1.1(X_H/m_p)^2
      // and we convert to erg/cm2/s/sq.arcsec.
      // Also we have an exponential cutoff in emissivity at 10^5K.
      //
      // If -DNII set, we assume N is strongly enhanced in the RSG
      // wind to f(N)=2.0e-4, or A(N)=8.3 (Brott+,2011), compared to
      // the ISM abundance of A(N)=7.85 or f(N)=7.08e-5, and we use
      // the second tracer to discriminate wind from ISM material.
      // Otherwise, use f(N)=7.08e-5 everywhere.
      //
      // Absorption from Henney et al. 2009, where they assume 
      //  alpha = 5.0e-22 nH per cm (absorption by dust). ~213.4
      //
      case PROJ_IML:
      prefactor = immpH*immpH*1.1*6.82e-18*sq_arcsec;
      //cout <<"[NII] prefactor="<<prefactor<<", expecting ~2.56e18";
      //cout <<".  Abs="<<immpH*5.0e-22<<", expecting 213.4\n";
      for (size_t i=0; i<Nr; i++) {
        ems[im][i] = prefactor *data[DATA_D][i]*data[DATA_TR0][i]*
                          data[DATA_D][i]*data[DATA_TR0][i]*
                          exp(-2.1855e4/data[DATA_T][i])/
                          sqrt(data[DATA_T][i])
                          *exp(-(data[DATA_T][i]*data[DATA_T][i])/1.0e10);
#ifdef NII
        ems[im][i] *= (1.0-data[DATA_TR1][i])*7.08e-5
                      +data[DATA_TR1][i]*2.0e-4;
#else
        ems[im][i] *= 7.08e-5;
#endif // NII

#ifdef ABSORPTION
        abs[im][i] = immpH*5.0e-22 *data[DATA_D][i];
#else
        abs[im][i] = 0.0;
#endif
      }
      break;

      default:
      cerr <<"get_var_and_slope(): Don't know what to do for var ";
      cerr <<im<<"\n";
      return im;
      break;
    }
  }
  return 0;
}




// ##################################################################
// ##################################################################



double calc_projection(
      const double *r, ///< radius array
      const double *v, ///< array of values at each radius
      const double *s, ///< array of slopes at each radius
      const size_t Nr, ///< Size of arrays.
      double b,        ///< impact parameter of ray.
      const double dr  ///< spacing of points in radius
      )
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr-1]+0.5*dr;
  if (b > grid_max) {
    cout <<"calc_projection: Bad B value, b="<<b<<"\n";
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    //cout <<"b="<<b<<" is <r[0], so resetting b=r[0]+eps.\n";
    b = r[0]*1.00000001;
  }

  //
  // start at b, integrate outwards to rmax, and then multiply by 2.
  //
  size_t ir = 0;
  while ( (r[ir]+dr) < b) ir++;

  double result = 0.0;
  double Rmin=0.0, Rmax=0.0, Rmin2=0.0, Rmax2=0.0;
  double maxd=0.0, mind=0.0;
  double slope=0.0, offset=0.0, xdx=0.0, x2dx=0.0;

  do {
    //
    // Min/Max for this line segment.
    //
    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    Rmax2 = Rmax*Rmax;
    Rmin2 = Rmin*Rmin;
    maxd  = sqrt(Rmax2 - b*b);
    mind  = sqrt(Rmin2 - b*b);

    //
    // var(r) = slope*r+offset, where a=|dvar/dr|_i, b=var_i-|dvar/dr|_i*r_i
    //
    slope = s[ir];
    offset = v[ir] -s[ir]*Rmin;
    //cout <<"\ta="<<a<<", b="<<b;

    //
    // xdx is the integral xdx/sqrt(x^2-y^2).
    // x2dx is the integral x^2dx/sqrt(x^2-y^2).
    //
    xdx = maxd - mind;
    x2dx = 0.5*(Rmax*maxd - Rmin*mind +b*b*log((Rmax+maxd)/(Rmin+mind)));
    //cout <<"\txdx="<<xdx<<", x2dx="<<x2dx;

    //
    // now this line segment's integral is a*x2dx + xdx*b
    //
    xdx = slope*x2dx + offset*xdx;
    //cout <<"\tint="<<xdx<<"\n";
    
    //
    // Add to result, and increment counter
    //
    result += xdx;
    ir++;
  } while (ir<Nr);

  //
  // double result for the outward ray (by symmetry).
  //
  result *= 2.0;

  return result;
}



// ##################################################################
// ##################################################################



double calc_projectionRT(
      const double *r, ///< radius array
      const double *ve, ///< array of emission values at each radius.
      const double *va, ///< array of absorption values at each radius.
      const size_t Nr, ///< Size of arrays.
      double b,        ///< impact parameter of ray.
      const double dr  ///< spacing of points in radius
      )
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr-1]+0.5*dr;
  if (b > grid_max) {
    cout <<"calc_projectionRT: Bad B value, b="<<b<<"\n";
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    //cout <<"b="<<b<<" is <r[0], so resetting b=r[0].\n";
    b = r[0];
  }

  //
  // start at Rmax, integrate inwards to b, and then back out.
  // For projections with emission and absorption we use the solution
  // to the equation of radiative transfer for a constant source
  // function through each line segment, with dTau modified to use
  // delta-ell (real path length) instead of delta-r (integrating
  // variable).  Rybicki & Lightman (1978), eq.1.30
  //
  long int ir = Nr-1;

  double result = 0.0;
  double Rmin=0.0, Rmax=0.0;
  //double maxd=0.0, mind=0.0;
  //double dIds=0.0;
  double dl=0.0, this_r=0.0;

  while (r[ir]*ONE_PLUS_EPS>b) {
    //
    // Min/Max for this line segment.
    //
    //cout <<"IN  b="<<b<<", r["<<ir<<"] = "<<r[ir]<<"\n";

    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    this_r = 0.5*(Rmin+Rmax);
    //
    // at the innermost interval, we evaluate dl exactly (using
    // Pythagoras), but otherwise use the approximate expression.
    //
    if (pconst.equalD(Rmin,b)) {
      // dr*sqrt(1+2b/dr)
      dl = (Rmax-Rmin)*sqrt(1.0+2.0*b/(Rmax-Rmin));
    }
    else {
      // r*dr/sqrt(r^2-b^2)
      dl = this_r*(Rmax-Rmin)/sqrt(this_r*this_r-b*b);
    }
    //
    // Then integrate along the line segment, and increment counter
    // source function is (ve[ir]/va[ir])
    //
    result = (ve[ir]/va[ir]) +exp(-va[ir]*dl)*(result- (ve[ir]/va[ir]));
    ir--;
    if (ir<0) {
      //cout <<"ir="<<ir<<"\n";
      break;
    }
  }

  //
  // Continue for the outward ray. We went one cell too far in the
  // inward loop, so increment before we begin.
  //
  ir++;
  while (ir<static_cast<long int>(Nr)) {
    //cout <<"OUT b="<<b<<", r["<<ir<<"] = "<<r[ir];

    //
    // Min/Max for this line segment.
    //
    Rmin = std::max(b    , r[ir]);
    Rmax = std::min(grid_max, r[ir]+dr);
    this_r = 0.5*(Rmin+Rmax);
    //
    // at the innermost interval, we evaluate dl exactly (using
    // Pythagoras), but otherwise use the approximate expression.
    //
    if (pconst.equalD(Rmin,b)) {
      // dr*sqrt(1+2b/dr)
      dl = (Rmax-Rmin)*sqrt(1.0+2.0*b/(Rmax-Rmin));
    }
    else {
      // r*dr/sqrt(r^2-b^2)
      dl = this_r*(Rmax-Rmin)/sqrt(this_r*this_r-b*b);
    }

    //
    // Then integrate along the line segment, and increment counter
    //
    result = (ve[ir]/va[ir]) +exp(-va[ir]*dl)*(result- (ve[ir]/va[ir]));


    ir++;
  }

  return result;
}


// ##################################################################
// ##################################################################




