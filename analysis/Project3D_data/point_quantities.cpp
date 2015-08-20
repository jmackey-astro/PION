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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"

#include "constants.h"
#include "sim_params.h"
#include "grid/cell_interface.h"
#include "microphysics/microphysics_base.h"

#include "point_quantities.h"

using namespace std;



// ##################################################################
// ##################################################################



double point_quantities::get_point_density(const struct point_4cellavg *pt)
{
  double val=0.0;
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) val += pt->wt[v] *pt->ngb[v]->P[RO];
  }
  return val;
}


// ##################################################################
// ##################################################################



double point_quantities::get_point_neutralH_numberdensity(
        const struct point_4cellavg *pt,
        const int ifrac
        )
{
  double val=0.0;
  //
  // First get neutral mass density.
  //
  for (int v=0;v<4;v++) {
    if (pt->ngb[v]) {
      val += pt->wt[v] *(pt->ngb[v]->P[RO]*(1.0-pt->ngb[v]->P[ifrac]));
    }
  }
  //
  // Then scale by H mass fraction and divide by mass of H.
  //
  val *= SimPM.EP.H_MassFrac/pconst.m_p();
  return val;
}




// ##################################################################
// ##################################################################



double point_quantities::get_point_temperature(
        const struct point_4cellavg *pt,
        const int ifrac
        )
{
  double val = 0.0;
  //
  // If microphysics is set up, then use MP->Temperature() to get the
  // temperature.  Otherwise assume a pure Hydrogen gas.
  //
  if (MP) {
    for (int v=0;v<4;v++) {
      if (pt->ngb[v]) {
        val += pt->wt[v] *MP->Temperature(pt->ngb[v]->P,SimPM.gamma);
      }
    }
  }
  else {
    rep.error("get_point_temperature(): no microphysics class",1);
    //  
    // First get the mean of p/rho/(1+x)
    //
    for (int v=0;v<4;v++) {
      if (pt->ngb[v]) {
        //cout <<"p,ro,if = "<<pt->ngb[v]->P[PG]<<", "<< pt->ngb[v]->P[RO] <<", "<<pt->ngb[v]->P[ifrac]<<endl;
        val += pt->wt[v] *(pt->ngb[v]->P[PG]/pt->ngb[v]->P[RO]/(1.0+pt->ngb[v]->P[ifrac]));
      }
    }
    //
    // multiply by m_p/k_B = 1.67e-24/1.38e-16 = 1.21e-8
    //
    val *= 1.21e-8;
    //  cout <<"Temperature="<<val<<endl;
  }
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
    if (pt->ngb[v]) {
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
    if (pt->ngb[v]) {
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
    if (pt->ngb[v]) {
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
    if (pt->ngb[v]) {
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
      const int ifrac, ///< ifrac
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
    //
    // If point exists, add its contribution, with weight.
    //
    if (pt->ngb[v]) {
      val += pt->wt[v] *
          (sx*pt->ngb[v]->P[bx]*st +sz*pt->ngb[v]->P[bz]*ct) *
          (pt->ngb[v]->P[RO]*pt->ngb[v]->P[ifrac]);
    }
  }
  //
  // convert to cm^{-3} and return:
  //
  return val*SimPM.EP.H_MassFrac/pconst.m_p();
};


// ##################################################################
// ##################################################################


///
/// Get the absorption and emission coefficients for H-alpha
/// recombination radiation, according to Hummer94 and Henney et al.
/// (2005)'s generic formulae.
/// 
/// I found emissivities according to Storey & Hummer
/// (1995,MNRAS,272,41), and they drop linearly with temperature.
/// Not sure why that is, but it is for Ha,Hb,Hg, etc.
///
/// UPDATE: Replaced with a fit to Ostebrock (1989)'s tables, which
/// seem to be more reliable.
///
void point_quantities::get_point_Halpha_params(
        const struct point_4cellavg *pt, ///< point in question.
        const int ifrac, ///< index of Prim.Vector with Ion. fraction.
        double *alpha,   ///< absorption coefficient (photons/cm)
        double *j        ///< emission coeff (phot/cm^3/s/ster)
        )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn;
  T  = get_point_temperature(pt,ifrac);
  //
  // get_point_density() is mass density, so we divide by m_p().
  //
  nn = get_point_neutralH_numberdensity(pt,ifrac);
  ni = get_point_density(pt)*SimPM.EP.H_MassFrac/pconst.m_p();
  ni = std::max(0.0, ni-nn);
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
    // can get zero temperature if point is off-grid!!!
    *j = 0.0;
  }
  else {
    *j = 2.63e-33*ni*ni*exp(-0.9*log(T));
    // OLD VALUES *j = 2.056e-14*ni*ni/T;
    // Emissivity scales more steeply than recombination rate
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
        double *alpha,   ///< absorption coefficient (/cm)
        double *j        ///< emission coeff (erg/cm^3/s/sq.arcsec)
        )
{
  //
  // I need the H+ number density, neutral number density, and
  // temperature. ion density is calculated from (density minus
  // neutral_density).
  // 
  double T, ni, nn;
  T  = get_point_temperature(pt,ifrac);
  //
  // get_point_density() is mass density, so we divide by m_p().
  //
  nn = get_point_neutralH_numberdensity(pt,ifrac);
  ni = get_point_density(pt)*SimPM.EP.H_MassFrac/pconst.m_p();
  ni = std::max(0.0, ni-nn);
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
  // in units of erg/cm3/s/sr, and use
  //  n_e*n_p = rho^2 y^2 (X_H/m_p)^2
  // and we convert to erg/cm2/s/sq.arcsec.
  //
  // Furthermore, we assume N has its solar ISM abundance of
  // A(N)=7.85 or f(N)=7.08e-5.
  //
  if (T<1.0) {
    // can get zero temperature if point is off-grid!!!
    *j = 0.0;
  }
  else {
    *j = 9.03e-34*ni*ni*exp(-2.1855e4/T)/sqrt(T);
  }
  return;
}



// ##################################################################
// ##################################################################



