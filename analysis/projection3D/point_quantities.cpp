///
/// \file point_quantities.cpp
/// \author Jonathan Mackey
/// \date 2015.08.06
///
/// This file defines a number or functions to calculate physical
/// properties of a point in space, by bilinear interpolation from
/// four other coplanar points.
///
/// Modifications:
/// - 2015.08.06 JM: Trying to avoid code duplication in
///  point_velocity and image classes, so both can derive from this
///  class.
/// - 2015.08.19 JM: Added get_point_RotationMeasure()
/// - 2015.08.21 JM: MP->Temperature() is NOT threadsafe, so had to
///    hardcode the temperature estimate if using threads.
/// - 2015.10.13 JM: added 6GHz Bremsstrahlung and Emission measure
/// - 2018.01.25 JM: added functions to request n(H+),n(H0),n(e-)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"

#include "constants.h"
#include "sim_params.h"
#include "grid/cell_interface.h"
#include "microphysics/microphysics_base.h"

#include "point_quantities.h"
#include "../xray/xray_emission.h"
#include "../projection/projection_constants.h"
using namespace std;



// ##################################################################
// ##################################################################



double point_quantities::get_point_density(const struct point_4cellavg *pt)
{
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf)
      val += pt->wt[v] *pt->ngb[v]->P[RO];
  }
  return val;
}


// ##################################################################
// ##################################################################



double point_quantities::get_point_electron_numberdensity(
        const struct point_4cellavg *pt
        )
{
  double val=0.0;
  //
  // Use the microphysics function to get electron number density.
  //
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v] *MP->get_n_elec(pt->ngb[v]->P);
    }
  }
  return val;
}


// ##################################################################
// ##################################################################


double point_quantities::get_point_ionizedH_numberdensity(
        const struct point_4cellavg *pt
        )
{
  double val=0.0;
  //
  // Use the microphysics function to get ionized H number density.
  //
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v] *MP->get_n_Hplus(pt->ngb[v]->P);
    }
  }
  return val;
}


// ##################################################################
// ##################################################################


double point_quantities::get_point_neutralH_numberdensity(
        const struct point_4cellavg *pt
        )
{
  double val=0.0;
  //
  // Use the microphysics function to get neutral H number density.
  //
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v] *MP->get_n_Hneutral(pt->ngb[v]->P);
    }
  }
  return val;
}




// ##################################################################
// ##################################################################



double point_quantities::get_point_temperature(
        const struct point_4cellavg *pt,
        const double gamma   ///< EOS gamma
        )
{
  double val = 0.0;
  //
  // If microphysics is set up, then use MP->Temperature() to get the
  // temperature.  Otherwise assume a pure Hydrogen gas.
  //
  if (MP) {
    for (int v=0;v<4;v++) {
      if (pt->ngb[v] && pt->ngb[v]->isleaf) {
        val += pt->wt[v] *MP->Temperature(pt->ngb[v]->P,gamma);
        //if (!isfinite(MP->Temperature(pt->ngb[v]->P,gamma)) ||
        //    MP->Temperature(pt->ngb[v]->P,gamma)==0.0) {
        //  cout <<"Invalid Temperature in loop="<<val<<"  "<<gamma<<"  "<<endl;
        //  rep.printVec("pv", pt->ngb[v]->P, 5);
        //  rep.printVec("pos",pt->ngb[v]->pos, 3);
        //  //rep.printVec("img",pt->ngb[v]->Ph, SimPM.ndim);
        //}
      }
    }
  }
  else {
    rep.error("get_point_temperature(): no microphysics class",1);
  }
  //if (!isfinite(val))
  //cout <<"Invalid Temperature="<<val<<endl;
  return val;
}
    



// ##################################################################
// ##################################################################


//#define SQRT_DENSITY
//#define NO_DENSITY_WT
#define LINMAX_DENSITY
#define MAXDENS 25000.0

double point_quantities::get_point_StokesQ(
      struct point_4cellavg *pt,
      const int ,
      const int bx, const int by, const int bz,
      const int signx, const int , const int signz,
      const double sintht, const double costht
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //

  //
  // New Q calculated from Bx,By according to
  // Q = sum_i wt[i]*|f(n_H)*(Bx^2-By^2)/sqrt(Bx^2+By^2)|_i
  //
  double val=0.0, bx2=0.0, by2=0.0, btot=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      //
      // Get Bx along line of sight:
      //
      bx2 = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      bx2 *= bx2;
      //
      // By is easier, and so then is btot.
      //
      by2 = pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      btot = sqrt(bx2+by2);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p())*(bx2-by2)/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *(bx2-by2)/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *(bx2-by2)/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  //
  // End of new calc
  //

  return val;
}



// ##################################################################
// ##################################################################



double point_quantities::get_point_StokesU(
      struct point_4cellavg *pt,
      const int ,
      const int bx, const int by, const int bz,
      const int signx, const int signy, const int signz,
      const double sintht, const double costht
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //

  //
  // New Q calculated from Bx,By according to
  // Q = sum_i wt[i]*|f(n_H)*(2*Bx*By)/sqrt(Bx^2+By^2)|_i
  //
  double val=0.0, btot=0.0, bxy=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      //
      // Calculate Bx and add Bx^2 to btot
      //
      bxy = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      btot = bxy*bxy;
      //
      // Multipy Bx by By, add By^2 to btot
      //
      bxy *= signy*pt->ngb[v]->P[by];
      btot += pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      //
      // Btot = sqrt(Bx^2+By^2) now.  bxy=Bx*By.
      //
      btot = sqrt(btot);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p())*2.0*bxy/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *2.0*bxy/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *2.0*bxy/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  //
  // End of new calc
  //

  return val;
}



// ##################################################################
// ##################################################################



double point_quantities::get_point_BXabs(
      struct point_4cellavg *pt,
      const int,
      const int bx, const int by, const int bz,
      const int signx, const int , const int signz,
      const double sintht, const double costht
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0, btot=0.0, bx2=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      bx2 = signx*pt->ngb[v]->P[bx]*costht -signz*pt->ngb[v]->P[bz]*sintht;
      bx2 *= bx2;
      btot = sqrt(pt->ngb[v]->P[bx]*pt->ngb[v]->P[bx] +
		  pt->ngb[v]->P[by]*pt->ngb[v]->P[by] +
		  pt->ngb[v]->P[bz]*pt->ngb[v]->P[bz]);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p()) *bx2/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *bx2/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *bx2/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  return val;
}



// ##################################################################
// ##################################################################



double point_quantities::get_point_BYabs(
      struct point_4cellavg *pt,
      const int ,
      const int bx, const int by, const int bz,
      const int , const int , const int ,
      const double , const double
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and neighbouring
  // cells.
  //
  double val=0.0, btot=0.0, by2=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      by2 = pt->ngb[v]->P[by]*pt->ngb[v]->P[by];
      btot = sqrt(pt->ngb[v]->P[bx]*pt->ngb[v]->P[bx] +
		  by2 +
		  pt->ngb[v]->P[bz]*pt->ngb[v]->P[bz]);
#if defined (SQRT_DENSITY)
      val += pt->wt[v] *sqrt(pt->ngb[v]->P[RO]/pconst.m_p()) *by2/btot;
#elif defined (NO_DENSITY_WT)
      val += pt->wt[v] *by2/btot;
#elif defined (LINMAX_DENSITY)
      val += pt->wt[v] *std::min(pt->ngb[v]->P[RO]/pconst.m_p(),MAXDENS) *by2/btot;
#else
#error "try to define some sort of weighting for magnetic field!!!"
#endif
    }
  }
  return val;
}



// ##################################################################
// ##################################################################


double point_quantities::get_point_RotationMeasure(
      struct point_4cellavg *pt, ///< pt
      const int bx,    ///< bx index (image coords)
      const int bz,    ///< bz index (image coords)
      const int sx,    ///< sign(xx)
      const int sz,    ///< sign(zz)
      const double st, ///< sin(theta)
      const double ct  ///< cos(theta)
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and
  // neighbouring cells.
  // The point value is RM = 0.81*n_e(cm^{-3})*B_los/sqrt(4pi)
  // This gets multiplied at the end by the path length through each
  // element of the integral (hh), by sqrt(4pi), and divided by 1 pc.
  //
  // the sqrt(4pi) factor is because the code uses units for B that
  // are Gauss/sqrt(4pi), so that I don't have to do any multiplying
  // or dividing by these factors in the code.
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    // If point exists, add its contribution, with weight.
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v] *
          (sx*pt->ngb[v]->P[bx]*st +sz*pt->ngb[v]->P[bz]*ct) *
          MP->get_n_elec(pt->ngb[v]->P);
    }
  }
  return val;
};


// ##################################################################
// ##################################################################


///
/// Get the absorption and emission coefficients for H-alpha
/// recombination radiation, with a fit to Ostebrock (1989)'s tables.
///
void point_quantities::get_point_Halpha_params(
      const struct point_4cellavg *pt, ///< point in question.
      const int ifrac, ///< index of Prim.Vector with Ion. fraction.
      const double gamma,   ///< EOS gamma
      double *alpha,   ///< absorption coefficient (photons/cm)
      double *j        ///< emission coeff (phot/cm^3/s/ster)
      )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn, ne;
  T  = get_point_temperature(pt,gamma);
  nn = get_point_neutralH_numberdensity(pt);
  ni = get_point_ionizedH_numberdensity(pt);
  ne = get_point_electron_numberdensity(pt);
  //
  // First absorption, from Henney et al. (2009) assuming the opacity
  // is from dust, so that neutrals and ions both count.
  //
  *alpha = (ni+nn) *5.0e-22; // cgs units hardcoded.
  //*alpha = 0.0; // zero absorption (cuts out simulation edges).
  //
  // Emissivity in H-alpha is
  //   j = 1.12e-22*n_e*n_p/pow(T,0.9) erg/cm3/s/sr,
  // from Osterbrock (book, edition 2006, table 4.4).
  // Converted to per square arcsec this is
  //   j = 2.63e-33*n_e*n_p/pow(T,0.9) erg/cm3/s/sq.arcsec.
  // Assume n_e=n_p (i.e. ignore electrons from Helium).
  //
  if (T<1.0) {
    // can get zero temperature if point is off-grid.
    *j = 0.0;
  }
  else {
    *j = ni * ne * Halpha_emissivity(T); // from Xray class
  }
  return;
}



// ##################################################################
// ##################################################################



///
/// Get the absorption and emission coefficients for [N II] 6584AA
/// forbidden line emission, according to the formula in Dopita 
/// (1973,A&A,29,387)
/// 
void point_quantities::get_point_NII6584_params(
      const struct point_4cellavg *pt, ///< point in question.
      const int ifrac, ///< index of Prim.Vector with Ion. fraction.
      const double gamma,   ///< EOS gamma
      double *alpha,   ///< absorption coefficient (/cm)
      double *j        ///< emission coeff (erg/cm^3/s/sq.arcsec)
      )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn, ne;
  T  = get_point_temperature(pt,gamma);
  nn = get_point_neutralH_numberdensity(pt);
  ni = get_point_ionizedH_numberdensity(pt);
  ne = get_point_electron_numberdensity(pt);
  //
  // First absorption, from Henney et al. (2009) assuming the opacity
  // is from dust, so that neutrals and ions both count.
  //
  *alpha = (ni+nn) *5.0e-22; // cgs units hardcoded.
  //*alpha = 0.0; // zero absorption (cuts out simulation edges).
  //
  // Emissivity in [N II] 6584AA is
  //  j([NII] ll 6584) =
  //   6.82e-18 n_e*n_p*f(N)*exp(-chi/kT)/(4*pi*sqrt(T))
  // in units of erg/cm3/s/sr, and we convert to erg/cm2/s/sq.arcsec.
  // Furthermore, we assume N has its solar ISM abundance of
  // A(N)=7.85 or f(N)=7.08e-5.
  //
  if (T<1.0) {
    // can get zero temperature if point is off-grid!!!
    *j = 0.0;
  }
  else {
    //double n_N1p;
    //double fNp = 7.08e-5; // ISM abundance of Nitrogen, by number.
    //n_N1p = fNp * ni;
    *j = 7.08e-5 * ni * ne * NII6584_emissivity(T);
  }
  return;
}



// ##################################################################
// ##################################################################


double point_quantities::get_point_EmissionMeasure(
      const struct point_4cellavg *pt
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and
  // neighbouring cells.
  // The point value is EM = n_e^2(cm^{-3})
  // This gets multiplied at the end by the path length through each
  // element of the integral (hh) and divided by 1 pc.
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    // If point exists, add its contribution, with weight.
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v] 
              *pow(MP->get_n_elec(pt->ngb[v]->P),2.0);
    }
  }
  return val;
};



// ##################################################################
// ##################################################################


double point_quantities::get_point_Bremsstrahlung6GHz(
      const struct point_4cellavg *pt,
      const double gamma   ///< EOS gamma
      )
{
  //
  // Bilinear interpolation with pre-calculated weights and
  // neighbouring cells.
  // The point value is j = 3.27e-23 n_e^2 T^-0.35 nu^-0.1 (MJy/sr/cm)
  // This gets multiplied at the end by the path length through each
  // element of the integral (hh).
  //
  double val=0.0;
  for (int v=0;v<4;v++) {
    // If point exists, add its contribution, with weight.
    if (pt->ngb[v] && pt->ngb[v]->isleaf) {
      val += pt->wt[v]
              * MP->get_n_elec(pt->ngb[v]->P)
              * MP->get_n_elec(pt->ngb[v]->P)
              * Brems6GHz_emissivity(MP->Temperature(pt->ngb[v]->P,gamma));
    }
  }
  return val;
};



// ##################################################################
// ##################################################################

///
/// Get the emission coefficients for X-ray emissivity at E>0.1kev.
///
void point_quantities::get_point_Xray_params(
        const struct point_4cellavg *pt, ///< point in question.
        const int ifrac, ///< index of Prim.Vector with Ion. fraction.
        const int index, ///< which X-ray emissivity to use (index in array).
        const double gamma,   ///< EOS gamma
        double *alpha,   ///< absorption coefficient (/cm)
        double *j        ///< emission coefficient (erg/cm^3/s/sq.arcsec)
        )
{
  // Need the electron number density and temperature.
  double xr[8], T, ne, ni;
  double per_angle = 1.0/(4.0*pconst.pi()*pconst.sqasec_per_sr());
  T  = get_point_temperature(pt,gamma);
  ne = get_point_electron_numberdensity(pt);
  ni = get_point_ionizedH_numberdensity(pt);
 
  if (T<1.0) {
    // can get zero temperature if point is off-grid.
    *j = 0.0;
    *alpha = 0.0;
  }
  else {
    get_xray_emissivity(T,xr); // volume emissivity per e- per H+
    // prefactor converts to intensity in erg/cm2/s/square-arcsec,
    // once multiplied by a path length.
    *j = xr[index] * ne * ni * per_angle;
    *alpha = 0.0;
  }
  return;
}


// ##################################################################
// ##################################################################


