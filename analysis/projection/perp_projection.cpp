
///
/// file:    perp_projection.cpp
/// author:  Jonathan Mackey
/// date:    2019-11-26
///
/// Description: Routines for calculating quantities along rays that
/// are not perpendicular to the grid in 2D simulations.


#include <cmath>

#include <sstream>
using namespace std;
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "../xray/xray_emission.h"
#include "constants.h"
#include "sim_params.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#include <spdlog/spdlog.h>

#include "perp_projection.h"
#include "projection_constants.h"


// ##################################################################
// ##################################################################



int generate_perpendicular_image(
    class SimParams &SimPM,       ///< simulation parameters
    class microphysics_base *MP,  ///< microphysics class
    class GridBaseClass *grid,    ///< computational grid
    class Xray_emission &XR,      ///< pointer to class.
    int npix[],                   ///< Number of pixels in each direction
    size_t num_pix,               ///< total number of pixels
    size_t NIMG,                  ///< number of images to make
    double **img_array            ///< pointer to the image arrays.
)
{
  //
  // zero the image arrays.
  //
  for (size_t im = 0; im < static_cast<size_t>(NIMG); im++)
    for (size_t ip = 0; ip < num_pix; ip++)
      img_array[im][ip] = 0.0;


  //
  // Loop over grid... Each z,R corresponds to a pixel.
  // Outer loop runs over z.
  //
  cell *cz = grid->FirstPt();
  int iz   = 0;
  int err  = 0;

#ifdef PROJ_OMP
#pragma omp parallel for
#endif
  for (iz = 0; iz < npix[0]; iz++) {
    // cout <<"#+#+#+#+#+# New column, iz="<<iz<<"\n";
    err = calculate_column(
        SimPM, MP, cz, XR, grid, iz, npix[1], NIMG, npix, grid->DX(),
        img_array);
    if (0 != err)
      spdlog::error("{}: Expected {} but got {}", "calculate_column", 0, err);
    cz = grid->NextPt(cz, ZPcyl);
  }

#ifndef NDEBUG
  spdlog::debug(" - -- - waiting for {} threads to finish", iz);
#endif
#ifndef NDEBUG
  spdlog::debug(" - -- - All threads are finished.\n");
#endif

  //
  // Multiply X-ray intensity by 4*PI*(206165)^2*(3.086e18)^2 to
  // get units of luminosity per sq. parsec (erg/s/pc^2), so that
  // we can compare to Townsley et al. papers.
  // Otherwise it is erg/cm^2/s/sq.arcsec.
  //
  // for (size_t ipix=0; ipix<num_pix; ipix++) {
  //  img_array[PROJ_X01][ipix] *= 5.09e48;
  //  img_array[PROJ_X05][ipix] *= 5.09e48;
  //  img_array[PROJ_X10][ipix] *= 5.09e48;
  //  img_array[PROJ_X50][ipix] *= 5.09e48;
  //}

  return 0;
}



// ##################################################################
// ##################################################################



//
// calculate projection for a column of cells in the R-direction
//
int calculate_column(
    class SimParams &SimPM,       ///< simulation parameters
    class microphysics_base *MP,  ///< microphysics class
    class cell *cz,               ///< current cell at start of column.
    class Xray_emission &XR,      ///< Xray emission class.
    class GridBaseClass *grid,    ///< pointer to grid.
    int iz,                       ///< counter of number of columns completed.
    int N_R,                      ///< Number of cells in radial direction
    int n_images,                 ///< Number of images.
    int *npix,                    ///< Number of pixels
    double delr,                  ///< Cell diameter
    double **img_array            ///< image array.
)
{
  //
  // Allocate memory for
  // columns of emission/absorption data for all output variables.
  //
  int err           = 0;
  double **ems_data = 0;
  double **abs_data = 0;
  double **raw_data = 0;
  ems_data          = mem.myalloc(ems_data, n_images);
  for (int im = 0; im < n_images; im++) {
    ems_data[im] = mem.myalloc(ems_data[im], N_R);
    for (int v = 0; v < N_R; v++)
      ems_data[im][v] = 0.0;
  }
  abs_data = mem.myalloc(abs_data, n_images);
  for (int im = 0; im < n_images; im++) {
    abs_data[im] = mem.myalloc(abs_data[im], N_R);
    for (int v = 0; v < N_R; v++)
      abs_data[im][v] = 0.0;
  }
#define NVAR 7
  raw_data = mem.myalloc(raw_data, NVAR);
  for (int im = 0; im < NVAR; im++) {
    raw_data[im] = mem.myalloc(raw_data[im], N_R);
    for (int v = 0; v < N_R; v++)
      raw_data[im][v] = 0.0;
  }

  // check for H+ tracer and Wind tracer
  size_t iHp = 0, iWf = 0;
  if (!MP)
    iHp = iWf = 0;
  else {
    int tr = MP->Tr("H1+");
    if (tr == DONT_CALL_ME || tr < 0 || !isfinite(tr)) {
#ifndef NDEBUG
      spdlog::warn("No H+ tracer variable, assuming all gas is ionized");
#endif
      iHp = 0;
    }
    else {
      iHp = tr;
    }
    tr = MP->Tr("WIND");
    if (tr == DONT_CALL_ME || tr < 0 || !isfinite(tr)) {
#ifndef NDEBUG
      spdlog::warn("No WIND tracer variable, assuming all gas is ISM");
#endif
      iWf = 0;
    }
    else {
      iWf = tr;
    }
  }

  // Starting cell.
  cell *cy = cz;
  int iy   = 0;
  std::array<double, MAX_DIM> cpos;

  // cout <<"\t\t assigning data from column to arrays.\n";
  do {
    CI.get_dpos(cy, cpos);
    raw_data[DATA_R][iy] = cpos[Rcyl];
    raw_data[DATA_D][iy] = cy->P[RO];
    raw_data[DATA_P][iy] = cy->P[PG];
    raw_data[DATA_V][iy] = cy->P[VY];
    //
    // If there is no H+ tracer, assume all gas is ionized.
    //
    if (iHp > 0) {
      raw_data[DATA_TR0][iy] = cy->P[iHp];
    }
    else {
      raw_data[DATA_TR0][iy] = 1.0;
    }
    //
    // If there is no Wind colour tracer, assume there is no wind.
    //
    if (iWf > 0) {
      raw_data[DATA_TR1][iy] = cy->P[iWf];
    }
    else {
      raw_data[DATA_TR1][iy] = 0.0;
    }
    //
    // Get Temperature from microphysics, or set it to sound speed.
    //
    if (MP)
      raw_data[DATA_T][iy] = MP->Temperature(cy->P, SimPM.gamma);
    else
      raw_data[DATA_T][iy] = sqrt(SimPM.gamma * cy->P[PG] / cy->P[RO]);
    // increment iy
    iy++;
  } while ((cy = grid->NextPt(cy, RPcyl)) != 0 && cy->isgd);
  if (iy != N_R)
    spdlog::error("{}: {}", "Bad logic for radial grid size", iy - N_R);

  //
  // Now for this radial column of data, we set the value for
  // each output variable at the cell centre, including both an
  // emission and an absorption coefficient for some variables.
  //
  // cout <<"\t\t getting emission/absorption data.\n";
  err += get_emission_absorption_data(
      SimPM, MP, grid, cz, raw_data, n_images, N_R, XR, ems_data, abs_data);

  //
  // Now go back to start of column in R, and for each pixel,
  // calculate the integral along a line of sight with impact
  // parameter b=R_i, for each variable.
  //
  double b = 0;

  for (int ivar = 0; ivar < n_images; ivar++) {
    // cout <<"\t\t Calculating image im="<<ivar<<",
    // name="<<im_name[ivar]<<"\n";
    for (int ipix = 0; ipix < N_R; ipix++) {
      b = raw_data[DATA_R][ipix];
      if (ivar < N_HD_SCALAR) {
#ifndef NDEBUG
        spdlog::debug(
            "ivar={}, ipix = {}: itotal={}, numpix={}", ivar, ipix,
            npix[Zcyl] * ipix + iz, npix[0] * npix[1]);
#endif
        img_array[ivar][npix[Zcyl] * ipix + iz] = calc_projection_column(
            raw_data[DATA_R], ems_data[ivar], abs_data[ivar], N_R, b, delr);
      }
      else {
#ifdef ABSORPTION
        img_array[ivar][npix[Zcyl] * ipix + iz] = calc_projectionRT_column(
            raw_data[DATA_R], ems_data[ivar], abs_data[ivar], N_R, b, delr);
#else
        img_array[ivar][npix[Zcyl] * ipix + iz] = calc_projection_column(
            raw_data[DATA_R], ems_data[ivar], abs_data[ivar], N_R, b, delr);
#endif
      }
    }
  }

  //
  // free memory.
  //
  for (int v = 0; v < n_images; v++) {
    ems_data[v] = mem.myfree(ems_data[v]);
    abs_data[v] = mem.myfree(abs_data[v]);
  }
  ems_data = mem.myfree(ems_data);
  abs_data = mem.myfree(abs_data);
  for (int v = 0; v < NVAR; v++) {
    raw_data[v] = mem.myfree(raw_data[v]);
  }
  raw_data = mem.myfree(raw_data);

  return err;
}



// ##################################################################
// ##################################################################



int get_emission_absorption_data(
    class SimParams &SimPM,       ///< simulation parameters
    class microphysics_base *MP,  ///< microphysics class
    class GridBaseClass *grid,    ///< pointer to grid.
    class cell *cz,               ///< cell at start of column.
    double const *const *data,    ///< raw data to get variable from
    const int n_img,              ///< number of images to write
    const size_t Nr,              ///< Number of radial data elements
    class Xray_emission &XR,      ///< pointer to X-ray emission class.
    double **ems,                 ///< array for emission[img][rad] data.
    double **abs                  ///< array for absorption[img][rad] data.
)
{
  double xr[8];
  for (int v = 0; v < 8; v++)
    xr[v] = 0.0;
  double n_e = 0.0, n_Hp = 0.0, n_H0 = 0.0, T = 0.0;
  double n_N1p;
  double fNp       = 7.08e-5;  // ISM abundance of Nitrogen, by number.
  double per_angle = 1.0 / (4.0 * pconst.pi() * pconst.sqasec_per_sr());


  // loop over all cells in column and set cell-centre quantities
  cell *c = cz;
  for (size_t i = 0; i < Nr; i++) {

    // calculate data for each cell:
    n_e  = MP->get_n_elec(c->P);
    n_Hp = MP->get_n_Hplus(c->P);
    n_H0 = MP->get_n_Hneutral(c->P);
    // T = MP->Temperature(c->P,SimPM.gamma);
    T = data[DATA_T][i];
    XR.get_xray_emissivity(T, xr);
    if (MP->Tr("N1p") > 0)
      n_N1p = MP->get_n_ion("N1p", c->P);
    else
      n_N1p = fNp * n_Hp;

    ems[PROJ_D][i]         = c->P[RO];
    ems[PROJ_NtD][i]       = n_H0;
    ems[PROJ_InD][i]       = n_Hp;
    ems[PROJ_EM][i]        = n_e * n_e / pconst.parsec();
    ems[PROJ_X00p1][i]     = n_e * n_Hp * xr[0] * per_angle;
    ems[PROJ_X00p2][i]     = n_e * n_Hp * xr[1] * per_angle;
    ems[PROJ_X00p3][i]     = n_e * n_Hp * xr[2] * per_angle;
    ems[PROJ_X00p5][i]     = n_e * n_Hp * xr[3] * per_angle;
    ems[PROJ_X01p0][i]     = n_e * n_Hp * xr[4] * per_angle;
    ems[PROJ_X02p0][i]     = n_e * n_Hp * xr[5] * per_angle;
    ems[PROJ_X05p0][i]     = n_e * n_Hp * xr[6] * per_angle;
    ems[PROJ_X10p0][i]     = n_e * n_Hp * xr[7] * per_angle;
    ems[PROJ_HA][i]        = n_e * n_Hp * XR.Halpha_emissivity(T);
    ems[PROJ_NII][i]       = n_e * n_N1p * XR.NII6584_emissivity(T);
    ems[PROJ_BREMS6GHZ][i] = n_e * n_Hp * XR.Brems6GHz_emissivity(T);

    c = grid->NextPt(c, RPcyl);
  }

  // 2nd loop over all cells in column to set slope/abs quantities
  c         = cz;
  double dr = 0.0;
  for (size_t i = 0; i < Nr; i++) {

    // calculate data for each cell:
    // n_e = MP->get_n_elec(c->P);
    n_Hp = MP->get_n_Hplus(c->P);
    n_H0 = MP->get_n_Hneutral(c->P);

    if (i < Nr - 1) {
      dr                 = data[DATA_R][i + 1] - data[DATA_R][i];
      abs[PROJ_D][i]     = (ems[PROJ_D][i + 1] - ems[PROJ_D][i]) / dr;
      abs[PROJ_NtD][i]   = (ems[PROJ_NtD][i + 1] - ems[PROJ_NtD][i]) / dr;
      abs[PROJ_InD][i]   = (ems[PROJ_InD][i + 1] - ems[PROJ_InD][i]) / dr;
      abs[PROJ_EM][i]    = (ems[PROJ_EM][i + 1] - ems[PROJ_EM][i]) / dr;
      abs[PROJ_X00p1][i] = (ems[PROJ_X00p1][i + 1] - ems[PROJ_X00p1][i]) / dr;
      abs[PROJ_X00p2][i] = (ems[PROJ_X00p2][i + 1] - ems[PROJ_X00p2][i]) / dr;
      abs[PROJ_X00p3][i] = (ems[PROJ_X00p3][i + 1] - ems[PROJ_X00p3][i]) / dr;
      abs[PROJ_X00p5][i] = (ems[PROJ_X00p5][i + 1] - ems[PROJ_X00p5][i]) / dr;
      abs[PROJ_X01p0][i] = (ems[PROJ_X01p0][i + 1] - ems[PROJ_X01p0][i]) / dr;
      abs[PROJ_X02p0][i] = (ems[PROJ_X02p0][i + 1] - ems[PROJ_X02p0][i]) / dr;
      abs[PROJ_X05p0][i] = (ems[PROJ_X05p0][i + 1] - ems[PROJ_X05p0][i]) / dr;
      abs[PROJ_X10p0][i] = (ems[PROJ_X10p0][i + 1] - ems[PROJ_X10p0][i]) / dr;
#ifdef ABSORPTION
      abs[PROJ_HA][i]  = 5.0e-22 * (n_Hp + n_H0);
      abs[PROJ_NII][i] = 5.0e-22 * (n_Hp + n_H0);
#else
      abs[PROJ_HA][i] = 0.0;
      abs[PROJ_NII][i] = 0.0;
#endif
      abs[PROJ_BREMS6GHZ][i] =
          (ems[PROJ_BREMS6GHZ][i + 1] - ems[PROJ_BREMS6GHZ][i]) / dr;
    }
    else {
      // last point, use backwards differencing.
      abs[PROJ_D][i]     = abs[PROJ_D][i - 1];
      abs[PROJ_NtD][i]   = abs[PROJ_NtD][i - 1];
      abs[PROJ_InD][i]   = abs[PROJ_InD][i - 1];
      abs[PROJ_EM][i]    = abs[PROJ_EM][i - 1];
      abs[PROJ_X00p1][i] = abs[PROJ_X00p1][i - 1];
      abs[PROJ_X00p2][i] = abs[PROJ_X00p2][i - 1];
      abs[PROJ_X00p3][i] = abs[PROJ_X00p3][i - 1];
      abs[PROJ_X00p5][i] = abs[PROJ_X00p5][i - 1];
      abs[PROJ_X01p0][i] = abs[PROJ_X01p0][i - 1];
      abs[PROJ_X02p0][i] = abs[PROJ_X02p0][i - 1];
      abs[PROJ_X05p0][i] = abs[PROJ_X05p0][i - 1];
      abs[PROJ_X10p0][i] = abs[PROJ_X10p0][i - 1];
#ifdef ABSORPTION
      abs[PROJ_HA][i]  = 5.0e-22 * (n_Hp + n_H0);
      abs[PROJ_NII][i] = 5.0e-22 * (n_Hp + n_H0);
#else
      abs[PROJ_HA][i] = 0.0;
      abs[PROJ_NII][i] = 0.0;
#endif
      abs[PROJ_BREMS6GHZ][i] = abs[PROJ_BREMS6GHZ][i - 1];
    }

    c = grid->NextPt(c, RPcyl);
  }  // 2nd loop

  return 0;
}



// ##################################################################
// ##################################################################


double calc_projection_column(
    const double *r,  ///< radius array
    const double *v,  ///< array of values at each radius
    const double *s,  ///< array of slopes at each radius
    const size_t Nr,  ///< Size of arrays.
    double b,         ///< impact parameter of ray.
    const double dr   ///< spacing of points in radius
)
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr - 1] + 0.5 * dr;
  if (b > grid_max) {
    spdlog::debug("calc_projection_column: Bad B value, b={}", b);
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    // cout <<"b="<<b<<" is <r[0], so resetting b=r[0]+eps.\n";
    b = r[0] * 1.00000001;
  }

  //
  // start at b, integrate outwards to rmax, and then multiply by 2.
  //
  size_t ir = 0;
  while ((r[ir] + dr) < b)
    ir++;

  double result = 0.0;
  double Rmin = 0.0, Rmax = 0.0, Rmin2 = 0.0, Rmax2 = 0.0;
  double maxd = 0.0, mind = 0.0;
  double slope = 0.0, offset = 0.0, xdx = 0.0, x2dx = 0.0;

  do {
    //
    // Min/Max for this line segment.
    //
    Rmin  = std::max(b, r[ir]);
    Rmax  = std::min(grid_max, r[ir] + dr);
    Rmax2 = Rmax * Rmax;
    Rmin2 = Rmin * Rmin;
    maxd  = sqrt(Rmax2 - b * b);
    mind  = sqrt(Rmin2 - b * b);

    //
    // var(r) = slope*r+offset, where a=|dvar/dr|_i, b=var_i-|dvar/dr|_i*r_i
    //
    slope  = s[ir];
    offset = v[ir] - s[ir] * Rmin;
    // cout <<"\ta="<<a<<", b="<<b;

    //
    // xdx is the integral xdx/sqrt(x^2-y^2).
    // x2dx is the integral x^2dx/sqrt(x^2-y^2).
    //
    xdx  = maxd - mind;
    x2dx = 0.5
           * (Rmax * maxd - Rmin * mind
              + b * b * log((Rmax + maxd) / (Rmin + mind)));
    // cout <<"\txdx="<<xdx<<", x2dx="<<x2dx;

    //
    // now this line segment's integral is a*x2dx + xdx*b
    //
    xdx = slope * x2dx + offset * xdx;
    // cout <<"\tint="<<xdx<<"\n";

    //
    // Add to result, and increment counter
    //
    result += xdx;
    ir++;
  } while (ir < Nr);

  //
  // double result for the outward ray (by symmetry).
  //
  result *= 2.0;

  return result;
}



// ##################################################################
// ##################################################################



double calc_projectionRT_column(
    const double *r,   ///< radius array
    const double *ve,  ///< array of emission values at each radius.
    const double *va,  ///< array of absorption values at each radius.
    const size_t Nr,   ///< Size of arrays.
    double b,          ///< impact parameter of ray.
    const double dr    ///< spacing of points in radius
)
{
  //
  // N.B. if b > Nr*dr then return zero!
  //
  double grid_max = r[Nr - 1] + 0.5 * dr;
  if (b > grid_max) {
    spdlog::debug("calc_projectionRT_column: Bad B value, b={}", b);
    return 0.0;
  }

  //
  // If simulation doesn't start at the origin, we need to fill in
  // the empty gap at the centre by setting b = r[0]
  //
  if (b < r[0]) {
    // cout <<"b="<<b<<" is <r[0], so resetting b=r[0].\n";
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
  long int ir = Nr - 1;

  double result = 0.0;
  double Rmin = 0.0, Rmax = 0.0;
  // double maxd=0.0, mind=0.0;
  // double dIds=0.0;
  double dl = 0.0, this_r = 0.0;

  while (r[ir] * ONE_PLUS_EPS > b) {
    //
    // Min/Max for this line segment.
    //
    // cout <<"IN  b="<<b<<", r["<<ir<<"] = "<<r[ir]<<"\n";

    Rmin   = std::max(b, r[ir]);
    Rmax   = std::min(grid_max, r[ir] + dr);
    this_r = 0.5 * (Rmin + Rmax);
    //
    // at the innermost interval, we evaluate dl exactly (using
    // Pythagoras), but otherwise use the approximate expression.
    //
    if (pconst.equalD(Rmin, b)) {
      // dr*sqrt(1+2b/dr)
      dl = (Rmax - Rmin) * sqrt(1.0 + 2.0 * b / (Rmax - Rmin));
    }
    else {
      // r*dr/sqrt(r^2-b^2)
      dl = this_r * (Rmax - Rmin) / sqrt(this_r * this_r - b * b);
    }
    //
    // Then integrate along the line segment, and increment counter
    // source function is (ve[ir]/va[ir])
    //
    result =
        (ve[ir] / va[ir]) + exp(-va[ir] * dl) * (result - (ve[ir] / va[ir]));
    ir--;
    if (ir < 0) {
      // cout <<"ir="<<ir<<"\n";
      break;
    }
  }

  //
  // Continue for the outward ray. We went one cell too far in the
  // inward loop, so increment before we begin.
  //
  ir++;
  while (ir < static_cast<long int>(Nr)) {
    // cout <<"OUT b="<<b<<", r["<<ir<<"] = "<<r[ir];

    //
    // Min/Max for this line segment.
    //
    Rmin   = std::max(b, r[ir]);
    Rmax   = std::min(grid_max, r[ir] + dr);
    this_r = 0.5 * (Rmin + Rmax);
    //
    // at the innermost interval, we evaluate dl exactly (using
    // Pythagoras), but otherwise use the approximate expression.
    //
    if (pconst.equalD(Rmin, b)) {
      // dr*sqrt(1+2b/dr)
      dl = (Rmax - Rmin) * sqrt(1.0 + 2.0 * b / (Rmax - Rmin));
    }
    else {
      // r*dr/sqrt(r^2-b^2)
      dl = this_r * (Rmax - Rmin) / sqrt(this_r * this_r - b * b);
    }

    //
    // Then integrate along the line segment, and increment counter
    //
    result =
        (ve[ir] / va[ir]) + exp(-va[ir] * dl) * (result - (ve[ir] / va[ir]));


    ir++;
  }

  return result;
}



// ##################################################################
// ##################################################################
