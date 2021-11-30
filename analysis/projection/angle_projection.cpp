
///
/// file:    angle_projection.cpp
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
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include "angle_projection.h"
#include "projection_constants.h"

#ifdef PROJ_OMP
#include <omp.h>
#endif


// ##################################################################
// ##################################################################

//
// pixel coords: i=x-axis (Z), j=y-axis (R)
// Xmin = SimPM.Xmin + n_extra*DX.
// Xmax = SimPM.Xmax + n_extra*DX.
//
void pixel_centre(
    double *Xmin,   // SimPM.Xmin
    double dx,      // grid size, dx
    int n_extra,    // n_extra: number of extra cells in Z-dir off end of grid
    int i,          // i, pixel index
    int j,          // j, pixel index
    double *centre  // pixel centre [OUTPUT]
)
{
  double zmin  = Xmin[Zcyl] - n_extra * dx;
  centre[Zcyl] = zmin + (i + 0.5) * dx;
  centre[Rcyl] = Xmin[Rcyl] + (j + 0.5) * dx;
  return;
}


// ##################################################################
// ##################################################################


//
// Get Z-coordinate of ray, given pixel coordinates and an R-coord.
//
double get_Z_from_R(
    double *pos,   ///< position vector of pixel/cell
    double R,      ///< R coordinate of point sought.
    double angle,  ///< angle of LOS with respect to +ve z-axis (radians)
    int inout      ///< =1 for incoming ray, -1 for outgoing.
)
{
  if (!isfinite(tan(angle))) {
    spdlog::error("get_Z_from_R() Request Z coord when theta=pi/2");
    return pos[Zcyl];
  }
  else {
    return pos[Zcyl] + inout * sqrt(R * R - pos[Rcyl] * pos[Rcyl]) / tan(angle);
  }
}


// ##################################################################
// ##################################################################


//
// Get R-coordinate of ray, given pixel coordinates and a Z-coord.
//
double get_R_from_Z(
    double *pos,  ///< position vector of pixel/cell
    double Z,     ///< Z coordinate of point sought.
    double angle  ///< angle of LOS with respect to +ve z-axis (radians)
)
{
  if (!isfinite(tan(angle))) {
    spdlog::error("get_R_from_Z() all R have same Z when theta=pi/2");
    return -1.0e200;
  }
  return sqrt(
      pos[Rcyl] * pos[Rcyl]
      + (Z - pos[Zcyl]) * (Z - pos[Zcyl]) * tan(angle) * tan(angle));
}



// ##################################################################
// ##################################################################



//
// Get cell containing start of ray (most distant from observer),
// given pixel coordinates and an angle.
//
cell *get_start_of_ray(
    class GridBaseClass *grid,  ///< pointer to grid
    double *pos,                ///< position vector of pixel/cell
    double angle  ///< angle of ray with respect to radial direction
)
{
  //
  // Try assuming ray enters grid from R-pos boundary
  //
  int inout = 0;
  if (angle < pconst.halfpi())
    inout = 1;
  else
    inout = -1;

  double r_ini = grid->SIM_Xmax(Rcyl);
  double z_ini = get_Z_from_R(pos, r_ini, angle, inout);

  // Set it to z_ini or zmin depending on inout.
  if (z_ini > grid->SIM_Xmax(Zcyl)) {
    //
    // Assumption above was wrong --> ray enters grid from Z-pos boundary.
    //
    z_ini = grid->SIM_Xmax(Zcyl);
    r_ini = get_R_from_Z(pos, z_ini, angle);
  }
  else if (z_ini < grid->SIM_Xmin(Zcyl)) {
    //
    // Assumption above was wrong --> ray enters grid from Z-pos boundary.
    //
    z_ini = grid->SIM_Xmin(Zcyl);
    r_ini = get_R_from_Z(pos, z_ini, angle);
  }


  //
  // Get cell associated with this ray end-point
  //
  cell *c = grid->FirstPt();
  while ((CI.get_dpos(c, Rcyl) + 0.5 * grid->DX()) < r_ini
         && grid->NextPt(c, RPcyl) && grid->NextPt(c, RPcyl)->isgd) {
    c = grid->NextPt(c, RPcyl);
  }
  while ((CI.get_dpos(c, Zcyl) + 0.5 * grid->DX()) < z_ini
         && grid->NextPt(c, ZPcyl) && grid->NextPt(c, ZPcyl)->isgd) {
    c = grid->NextPt(c, ZPcyl);
  }

  if (!c) spdlog::error("{}: {}", "cell error", fmt::ptr(c));

#ifdef DEBUG
  if (grid->SIM_Xmax(Rcyl) - pos[Rcyl] < grid->DX()) {
    spdlog::debug(
        "{}, {}, {}, {}\n{}, {}, {}, {}", z_ini, grid->SIM_Xmin(Zcyl),
        grid->SIM_Xmax(Zcyl), CI.get_dpos(c, Zcyl), r_ini, grid->SIM_Xmin(Rcyl),
        grid->SIM_Xmax(Rcyl), CI.get_dpos(c, Rcyl));
  }
#endif

  return c;
}



// ##################################################################
// ##################################################################



//
// Calculate the positions where a (curved) ray enters and exits a
// cell in the (z,R) plane for the given input parameters.
//
void get_entry_exit_points(
    class GridBaseClass *grid,  ///< pointer to grid
    double *pos,                ///< position vector of pixel.
    double angle,           ///< angle of ray with respect to radial direction
    int inout,              ///< =1 for incoming ray, -1 for outgoing.
    cell *ray,              ///< cell we are at.
    vector<double> &entry,  ///< entry point of ray into cell  (output)
    vector<double> &exit,   ///< exit point of ray from cell  (output)
    cell **next             ///< cell the ray exits into.              (output)
)
{
  std::array<double, MAX_DIM> cpos;
  CI.get_dpos(ray, cpos);
  double dxo2  = 0.5 * grid->DX();
  double dx    = grid->DX();
  bool flag    = false;
  int dir_flag = inout;

  //
  // If theta<pi/2, then the ray starts at z=infinity and goes to
  // z=-infinity (by convention of how theta is defined).
  // Otherwise, ray goes from z=-infty to z=+infty.
  // Set entry z-value to be max/min possible: z +/- dx/2
  //
  if (angle < pconst.halfpi()) {
    entry[Zcyl] = cpos[Zcyl] + dxo2;
    exit[Zcyl]  = cpos[Zcyl] - dxo2;
  }
  else {
    entry[Zcyl] = cpos[Zcyl] - dxo2;
    exit[Zcyl]  = cpos[Zcyl] + dxo2;
  }

  // R coordinate of entry point is unique for given z.
  entry[Rcyl] = get_R_from_Z(pos, entry[Zcyl], angle);
  exit[Rcyl]  = get_R_from_Z(pos, exit[Zcyl], angle);
  // cout <<"1st R: "<< entry[Rcyl] <<"  "<< exit[Rcyl]<<"  "<<cpos[Rcyl] <<"
  // "<< pos[Rcyl]  <<endl;
  if (pconst.equalD(exit[Rcyl], entry[Rcyl])
      && pconst.equalD(pos[Zcyl], cpos[Zcyl])) {
    // maybe at start cell, or else at steep angle, so try and see
    // if there is another cell at RN/RP that is on the ray.
    if (!pconst.equalD(pos[Rcyl], cpos[Rcyl])) {
      flag = true;
    }
  }

  // if entry, exit, and midpt (below) have R coordinate outside of
  // the cell, then bug out (set outputs to impossible values)
  //
  double midpt = get_R_from_Z(pos, cpos[Zcyl], angle);
#ifdef DEBUG
  spdlog::debug("pos={}, cpos={}", pos[Zcyl], cpos[Zcyl]);
#endif
  double mn = (cpos[Rcyl] - dxo2) * ONE_MINUS_EPS;
  double mx = (cpos[Rcyl] + dxo2) * ONE_PLUS_EPS;
#ifdef DEBUG
  spdlog::debug(
      "Check reasonable: {}  {}  {}  {}  {}", entry[Rcyl], exit[Rcyl], midpt,
      mn, mx);
#endif
  if ((entry[Rcyl] < mn && exit[Rcyl] < mn && midpt < mn)
      || (entry[Rcyl] > mx && exit[Rcyl] > mx && midpt > mx)) {
#ifdef DEBUG
    spdlog::warn("impossible values");
#endif
    *next       = 0;
    entry[Rcyl] = 1.0e99;
    exit[Rcyl]  = 1.0e99;
    entry[Zcyl] = 1.0e99;
    exit[Zcyl]  = 1.0e99;
    return;
  }


  // ENTRY:
  // If R is outside the cell, then the ray must enter the cell
  // from R +/- dx/2, not from z +/- dx/2, so we reset entry points.
  if (entry[Rcyl] < cpos[Rcyl] - dxo2) {
    // entry point is below cell edge.
    entry[Rcyl] = cpos[Rcyl] - dxo2;
    if (inout == 0) dir_flag = 1;  // entry is on incoming ray
    entry[Zcyl] = get_Z_from_R(pos, entry[Rcyl], angle, dir_flag);
  }
  else if (entry[Rcyl] > cpos[Rcyl] + dxo2) {
    // entry point is above cell edge.
    entry[Rcyl] = cpos[Rcyl] + dxo2;
    if (inout == 0) dir_flag = 1;  // entry is on incoming ray
    entry[Zcyl] = get_Z_from_R(pos, entry[Rcyl], angle, dir_flag);
  }

  // EXIT:
  // If R is outside the cell, then the ray must exit the cell
  // from R +/- dx/2, not from z +/- dx/2, so we reset exit points.
  //
  if (exit[Rcyl] < cpos[Rcyl] - dxo2) {
    // exit point is below cell edge.
    exit[Rcyl] = cpos[Rcyl] - dxo2;
    if (inout == 0) dir_flag = -1;  // exit is on outgoing ray
    exit[Zcyl] = get_Z_from_R(pos, exit[Rcyl], angle, dir_flag);
  }
  else if (exit[Rcyl] > cpos[Rcyl] + dxo2) {
    // exit point is above cell edge.
    exit[Rcyl] = cpos[Rcyl] + dxo2;
    if (inout == 0) dir_flag = -1;  // exit is on outgoing ray
    exit[Zcyl] = get_Z_from_R(pos, exit[Rcyl], angle, dir_flag);
  }

  //
  // Now set the next cell in the ray, taking account of where the
  // ray leaves the cell.  If it leaves at a corner then we move
  // both vertically and horizontally, but have to be careful of
  // grazing incidence on the corner of a cell.
  //
  *next = ray;
  // use dir_flag for exiting ray, in case we are at turning point
  // and have changed the 'inout' flag above when calculating exit
  // point.
  if (fabs(exit[Rcyl] - (cpos[Rcyl] + dxo2)) / dx < 1.0e-8 && dir_flag == -1) {
    *next = grid->NextPt(ray, RPcyl);
  }
  else if (
      fabs(exit[Rcyl] - (cpos[Rcyl] - dxo2)) / dx < 1.0e-8 && dir_flag == 1) {
    *next = grid->NextPt(ray, RNcyl);
  }

  if (fabs(exit[Zcyl] - (cpos[Zcyl] + dxo2)) / dx < 1.0e-8
      && angle > pconst.halfpi()) {
    *next = grid->NextPt(*next, ZPcyl);
  }
  else if (
      fabs(exit[Zcyl] - (cpos[Zcyl] - dxo2)) / dx < 1.0e-8
      && angle < pconst.halfpi()) {
    *next = grid->NextPt(*next, ZNcyl);
  }

  // If the flag was set, then we must treat cells specially:
  // inward  flag, force exit to be in RN direction
  // outward flag, force entry to be in RP direction
  if (flag == true) {
    if (inout == 1) {
      exit[Rcyl] = cpos[Rcyl] - dxo2;
      exit[Zcyl] = get_Z_from_R(pos, exit[Rcyl], angle, inout);
      *next      = grid->NextPt(ray, RNcyl);
    }
    else if (inout == -1) {
      entry[Rcyl] = cpos[Rcyl] - dxo2;
      entry[Zcyl] = get_Z_from_R(pos, entry[Rcyl], angle, inout);
      if (pconst.equalD(exit[Rcyl], (cpos[Rcyl] + dxo2))) {
        *next = grid->NextPt(ray, RPcyl);
      }
    }
    return;
  }

  if (!(*next) || (*next)->isbd || !(*next)->isgd) *next = 0;
  return;
}



// ##################################################################
// ##################################################################


//
// Add emission (or whatever projected quantity is) from current cell
// assuming data within the cell is constant, and that the geometric
// part of the integral is pre-calculated in "integral" variable.
// So the integral reduces to the sum:
//   SUM_{cells} integral * emissivity(P)
//
void add_cell_emission_to_ray(
    class SimParams &SimPM,       ///< simulation parameters
    class microphysics_base *MP,  ///< microphysics class
    double integral,          ///< the path-length/geometry part of the integral
    double *P,                ///< State vector of cell (primitive vars)
    class Xray_emission &XR,  ///< pointer to class for getting X-ray em
    vector<double> &ans)
{

  double xr[8];
  for (int v = 0; v < 8; v++)
    xr[v] = 0.0;
  double n_e       = MP->get_n_elec(P);
  double n_Hp      = MP->get_n_Hplus(P);
  double n_H0      = MP->get_n_Hneutral(P);
  double T         = MP->Temperature(P, SimPM.gamma);
  double per_angle = 1.0 / (4.0 * pconst.pi() * pconst.sqasec_per_sr());
  // cout <<n_e<<", "<< n_H0 <<", "<< n_Hp <<", "<< T <<"\n";
  if (!isfinite(T)) {
    spdlog::error("{}: {}", "bugging out, T is infinite", T);
  }

  ans[PROJ_D] += integral * P[RO];

  ans[PROJ_NtD] += integral * n_H0;
  ans[PROJ_InD] += integral * n_Hp;
  ans[PROJ_EM] += integral * n_e * n_e;

  XR.get_xray_emissivity(T, xr);
  ans[PROJ_X00p1] += integral * n_e * n_Hp * xr[0] * per_angle;
  ans[PROJ_X00p2] += integral * n_e * n_Hp * xr[1] * per_angle;
  ans[PROJ_X00p3] += integral * n_e * n_Hp * xr[2] * per_angle;
  ans[PROJ_X00p5] += integral * n_e * n_Hp * xr[3] * per_angle;
  ans[PROJ_X01p0] += integral * n_e * n_Hp * xr[4] * per_angle;
  ans[PROJ_X02p0] += integral * n_e * n_Hp * xr[5] * per_angle;
  ans[PROJ_X05p0] += integral * n_e * n_Hp * xr[6] * per_angle;
  ans[PROJ_X10p0] += integral * n_e * n_Hp * xr[7] * per_angle;

  ans[PROJ_HA] += integral * n_e * n_Hp * XR.Halpha_emissivity(T);
  double n_N1p;
  double fNp = 7.08e-5;  // ISM abundance of Nitrogen, by number.
  if (MP->Tr("N1p") > 0)
    n_N1p = MP->get_n_ion("N1p", P);
  else
    n_N1p = fNp * n_Hp;
  ans[PROJ_NII] += integral * n_e * n_N1p * XR.NII6584_emissivity(T);
  ans[PROJ_BREMS6GHZ] += integral * n_e * n_Hp * XR.Brems6GHz_emissivity(T);

  return;
}



// ##################################################################
// ##################################################################



//
// Integrate all the lines of sight for data on the grid.
//
//#define DEBUG
int generate_angle_image(
    class SimParams &SimPM,       ///< simulation parameters
    class microphysics_base *MP,  ///< microphysics class
    class GridBaseClass *grid,    ///< computational grid
    class Xray_emission &XR,      ///< pointer to class.
    double angle,       ///< angle between LOS and symmetry axis [1,89]
    int npix[],         ///< Number of pixels in each direction
    size_t num_pix,     ///< total number of pixels
    int n_extra,        ///< number of extra pixels w.r.t. cells.
    size_t NIMG,        ///< number of images to make
    double **img_array  ///< pointer to the image arrays.
)
{

  double dx = grid->DX();

  //
  // Loop over all pixels
  //
  spdlog::info("looping over pixels");
#ifdef PROJ_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < npix[0]; i++) {
    //
    // Loop over pixels... Each z,R corresponds to a pixel.
    // Outer loop runs over z.  The pixel grid is a superset of the
    // computational grid, but each pixel should have a ray that
    // crosses at least one grid cell.  Each pixel is the size of a
    // grid cell, but they get narrower in the Z-direction as the
    // viewing angle deviates from 90 degrees.
    //
    size_t ipix = 0;
    std::array<double, MAX_DIM> pos, startpos;
    vector<double> entry(2), exit(2);
    double invst    = 1.0 / sin(angle);
    double integral = 0.0;

#ifdef DEBUG
    spdlog::debug("column = {}/{}", i, npix[0]);
#endif

    //
    // Loop over all pixels in column (R direction, for given Z)
    //
    for (int j = 0; j < npix[1]; j++) {
      int ncell = 0;
      pixel_centre(SimPM.Xmin.data(), dx, n_extra, i, j, &pos[0]);
#ifdef DEBUG
      //        double temp[3];
      spdlog::debug("Pixel: [{}, {}] of [{},{}]   ", i, j, npix[0], npix[1]);
      spdlog::debug("pos : {}", pos);
#endif

      cell *c         = grid->FirstPt();
      cell *ray       = get_start_of_ray(grid, &pos[0], angle);
      cell *next_cell = 0;
      CI.get_dpos(ray, startpos);
#ifdef DEBUG
      spdlog::debug("ray={}", ray);
      spdlog::debug("starting cell : {}", startpos);
#endif

      // if pixel is on the grid, then find the cell that this
      // corresponds to.  Important to know this because then the ray
      // has a turning point at this cell.  If not then we only have
      // one half of the ray (incoming or outgoing).
      //
      bool pix_on_grid = false;
      c                = 0;
      if (pos[Zcyl] < SimPM.Xmax[Zcyl] && pos[Zcyl] > SimPM.Xmin[Zcyl]) {
        pix_on_grid = true;
        c           = grid->FirstPt();
        while (!pconst.equalD(pos[Zcyl], CI.get_dpos(c, Zcyl)))
          c = grid->NextPt(c, ZPcyl);
        while (!pconst.equalD(pos[Rcyl], CI.get_dpos(c, Rcyl)))
          c = grid->NextPt(c, RPcyl);
      }
#ifdef DEBUG
      spdlog::debug("pix on grid? c={}", c);
#endif

      // set up array for answer
      vector<double> ans(NIMG);
      for (size_t v = 0; v < NIMG; v++)
        ans[v] = 0.0;
      int inout;  // to show whether on inward/outward trajectory.

      //
      // See if we have an incoming ray.  If so, then integrate it.
      //
      bool at_source = false;
      if (pos[Zcyl] < SimPM.Xmax[Zcyl] - dx) {
#ifdef DEBUG
        spdlog::debug("incoming ray starting");
#endif
        // We have an incoming ray.  Trace ray until we get to pos[]
        // or the edge of the grid.
        inout = 1;  // to show that we are on the inward trajectory.
        // first check that the ray passes through cell *ray:
        get_entry_exit_points(
            grid, &pos[0], angle, inout, ray, entry, exit, &next_cell);
        if (entry[Rcyl] > 1.0e90) {
#ifdef DEBUG
          spdlog::error("INCOMING RAY ERROR: doesn't intersect grid, skipping");
#endif
        }
        else {

          do {
            //
            // Get entry and exit points of ray, and pointer to next cell.
            //
            get_entry_exit_points(
                grid, &pos[0], angle, inout, ray, entry, exit, &next_cell);

            //
            // Integrate through cell along path (constant gas density)
            // int(k.dl) = int(k/sin(theta) R.dR/sqrt(R^2-y^2))
            //           = (k/sin(theta)).[sqrt(R+^2-y^2)-sqrt(R-^2-y^2)]
            //
            integral =
                invst
                * abs(sqrt(entry[Rcyl] * entry[Rcyl] - pos[Rcyl] * pos[Rcyl])
                      - sqrt(exit[Rcyl] * exit[Rcyl] - pos[Rcyl] * pos[Rcyl]));
            // For each emission variable, add to ans[] vector
            add_cell_emission_to_ray(SimPM, MP, integral, ray->P, XR, ans);
            if (ans[PROJ_NtD] < 0.0) {
              spdlog::debug(
                  "Ray state vec : {}",
                  std::vector<double>(ray->P, ray->P + SimPM.nvar));
              spdlog::debug(
                  "ray  {},  next_cell={}\n{}, {}", fmt::ptr(ray),
                  fmt::ptr(next_cell), ray->isgd, ray->isbd);
              spdlog::error("{}: {}", "error", 99);
            }

            //
            // move to next cell.
            //
            ray = next_cell;
            ncell++;
            if (pix_on_grid) {
              if (ray == c)
                at_source = true;
              else
                at_source = false;
            }
            else {
              if (!(ray) || !(ray->isgd)
                  || pconst.equalD(CI.get_dpos(ray, Zcyl), pos[Zcyl]))
                at_source = true;
              else
                at_source = false;
            }
          } while (!at_source);
#ifdef DEBUG
          spdlog::info("incoming ray finished");
#endif
        }
      }

      //
      // Now the turning point cell, if it exists.
      //
      if (pix_on_grid) {
#ifdef DEBUG
        spdlog::debug("pix on grid starting: ray={}", ray);
        std::array<double, MAX_DIM> temp;
        CI.get_dpos(ray, temp);
        spdlog::debug("cell pos on ray : {}", temp);
#endif
        //
        // The pixel/cell itself is a special case because we only go
        // to the cell-centre and back out again, and because the
        // slope of the path changes sign.
        //
        inout = 0;
        // spdlog::error("{}: {}", "got to cell",1);
        get_entry_exit_points(
            grid, &pos[0], angle, inout, ray, entry, exit, &next_cell);
        // Path length is twice the inward part.
        integral = invst * 2.0
                   * sqrt(entry[Rcyl] * entry[Rcyl] - pos[Rcyl] * pos[Rcyl]);
        // For each emission variable, add to ans[] vector
        add_cell_emission_to_ray(SimPM, MP, integral, ray->P, XR, ans);
        ray = next_cell;
        ncell++;
#ifdef DEBUG
        spdlog::debug("pix on grid finishing: ray=");
        if (ray) {
          spdlog::debug("{}", ray);
          CI.get_dpos(ray, temp);
          spdlog::debug("cell pos on ray : {}", temp);
        }
        else {
          spdlog::warn("0, doesn't exist!");
        }
#endif
      }
#ifdef DEBUG
      else
        spdlog::warn("  PIXEL NOT ON GRID  ");
#endif

      //
      // See if we have an outgoing ray.  If so, then integrate it.
      //
      if (pos[Zcyl] > SimPM.Xmin[Zcyl] + dx && ray != 0) {
#ifdef DEBUG
        spdlog::debug(
            "outgoing ray starting, ncell={}\npos[Zcyl] = {}, Xmin={}, Xmin+dx={}",
            ncell, pos[Zcyl], SimPM.Xmin[Zcyl], SimPM.Xmin[Zcyl] + dx);
#endif
        // We have an outgoing ray.  Trace ray until we get to the
        // edge of the grid.
        inout = -1;
        // first check that the ray passes through cell *ray:
        get_entry_exit_points(
            grid, &pos[0], angle, inout, ray, entry, exit, &next_cell);
        if (entry[Rcyl] > 1.0e90) {
#ifdef DEBUG
          spdlog::error(
              "  OUTGOING RAY ERROR\n ray={}, next_cell={} pos=", ray,
              next_cell);
          spdlog::debug(" : {}", pos);
          spdlog::error("{}: {}", "rays", 1);
#endif
        }
        else {
#ifdef DEBUG
          spdlog::debug("entering loop for outgoing ray, ncell={}", ncell);
#endif
          do {
            //
            // Get entry and exit points of ray, and pointer to next cell.
            //
#ifdef DEBUG
            spdlog::debug("getting entry/exit points of ray, ncell={}", ncell);
#endif
            get_entry_exit_points(
                grid, &pos[0], angle, inout, ray, entry, exit, &next_cell);
#ifdef DEBUG
            spdlog::debug("entry/exit next_cell={}", next_cell);
#endif
            //
            // Integrate through cell along path (constant gas density)
            // int(k.dl) = int(k/sin(theta) R.dR/sqrt(R^2-y^2))
            //           = (k/sin(theta)).[sqrt(R+^2-y^2)-sqrt(R-^2-y^2)]
            //
            integral =
                invst
                * abs(sqrt(entry[Rcyl] * entry[Rcyl] - pos[Rcyl] * pos[Rcyl])
                      - sqrt(exit[Rcyl] * exit[Rcyl] - pos[Rcyl] * pos[Rcyl]));
            // For each emission variable, add to ans[] vector
            add_cell_emission_to_ray(SimPM, MP, integral, ray->P, XR, ans);
            //
            // move to next cell.
            //
            ray = next_cell;
            ncell++;
          } while (ray != 0 && ray->isgd);
#ifdef DEBUG
          spdlog::debug("outgoing ray finishing");
#endif
        }
      }

      // Now we've traced the ray for pixel [i,j], so tidy up and
      // move on to the next pixel

      // convert emission measure to units cm^{-6}.pc
      ans[PROJ_EM] /= pconst.parsec();
      //
      // save results into pixel ipix
      //
      ipix = npix[0] * j + i;
      for (size_t v = 0; v < NIMG; v++)
        img_array[v][ipix] = ans[v];
#ifdef DEBUG
      spdlog::debug("ans : {}", ans);
      spdlog::debug("i,j={}, {} ... Ray crossed {} cells.", i, j, ncell);
#endif
    }
  }
  return 0;
}



// ##################################################################
// ##################################################################
