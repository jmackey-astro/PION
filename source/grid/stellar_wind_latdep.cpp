///
/// \file stellar_wind_latdep.cpp
/// \author Jonathan Mackey
/// \date 2021.02.03
///
/// Modifications:
/// - 2021.02.03 JM: getting it working.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "grid/grid_base_class.h"
#include "grid/stellar_wind_latdep.h"
#include "microphysics/microphysics_base.h"
#include <sstream>


#include "constants.h"
#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;


// ##################################################################
// ##################################################################


// Constructor
stellar_wind_latdep::stellar_wind_latdep(
    const int nd,           ///< ndim
    const int nv,           ///< nvar
    const int nt,           ///< ntracer
    const int ft,           ///< ftr
    const std::string *tr,  ///< List of tracer variable names.
    const int cs,           ///< coord_sys
    const int eq,           ///< eqn_type
    const double mt,        ///< Minimum temperature allowed on grid
    const double ss,        ///< Simulation start time.
    const double sf,        ///< Simulation finish time.
    const double xi         ///< exponent of wind equatorial enhancement
    ) :
    stellar_wind(nd, nv, nt, ft, tr, cs, eq, mt),
    stellar_wind_evolution(nd, nv, nt, ft, tr, cs, eq, mt, ss, sf)

{
  // Constants for wind functions
  stellar_wind_latdep::c_gamma = 0.35;
  stellar_wind_latdep::c_xi    = xi;

  // Number of points in theta, Omega and Teff vectors
  stellar_wind_latdep::Ntheta = 25;
  stellar_wind_latdep::NOmega = 25;

  // Generate interpolating tables
  setup_tables();
}


// ##################################################################
// ##################################################################


// Destructor
stellar_wind_latdep::~stellar_wind_latdep() {}


// ##################################################################
// ##################################################################



// Generate interpolating tables for wind density function
void stellar_wind_latdep::setup_tables()
{
  // Min, mid and max values of theta - 5pts between min and mid,
  // 20pts between mid and max
  double t_min = 0.00;
  double t_mid = 60.0;
  double t_max = 90.0;

  theta_vec.resize(Ntheta);
  for (int k = 0; k < Ntheta; k++) {
    if (k <= 4)
      theta_vec[k] =
          (t_min + k * ((t_mid - t_min) / 4.0)) * pconst.pi() / 180.0;
    else
      theta_vec[k] = (t_mid + (k - 4) * ((t_max - t_mid) / (Ntheta - 5.0)))
                     * pconst.pi() / 180.0;
  }

  // Set up Omega vector
  log_mu_vec.resize(NOmega);
  Omega_vec.resize(NOmega);
  // Iterate log(mu) values - evenly spaced
  for (int i = 0; i < NOmega; i++)
    log_mu_vec[NOmega - i - 1] = -4.0 + i * (4.0 / (NOmega - 1));
  // Iterate Omega values - log spacing - spacing decreases as
  // Omega approaches 1
  for (int i = 0; i < NOmega; i++)
    Omega_vec[i] = 1.0 - pconst.pow_fast(10.0, log_mu_vec[i]);

  // set up interpolating table for (1-Omega*sin(theta))^0.35
  // to get vinf(theta,Omega) via 2D interpolation.
  vinf_vec.resize(Ntheta);
  for (int i = 0; i < Ntheta; i++)
    vinf_vec[i].resize(NOmega);
  for (int i = 0; i < Ntheta; i++) {
    for (int j = 0; j < NOmega; j++) {
      // function normalised without the vinf.
      vinf_vec[i][j] = fn_v_inf(theta_vec[i], Omega_vec[j], 1.0);
    }
  }
  vinf_vsize.resize(2);
  vinf_vsize[0] = Ntheta;
  vinf_vsize[1] = NOmega;

  // set up interpolating table for wind density lat-dependence,
  // f(theta,Omega) = sin(theta)*(1-Omega*sin(theta))^xi
  // to get rho(r,theta,Omega) via 2D interpolation.
  f_vec.resize(Ntheta);
  for (int i = 0; i < Ntheta; i++)
    f_vec[i].resize(NOmega);
  for (int i = 0; i < Ntheta; i++) {
    for (int j = 0; j < NOmega; j++) {
      // function normalised without the vinf.
      f_vec[i][j] = f(theta_vec[i], Omega_vec[j]);
    }
  }
  f_vsize.resize(2);
  f_vsize[0] = Ntheta;
  f_vsize[1] = NOmega;

  // set up interpolating table for wind density normalisation,
  // int_0^pi/2 f(theta,Omega) sin(theta) dtheta
  // to get scaling factor A via interpolation.
  norm_vec.resize(NOmega);
  for (int j = 0; j < NOmega; j++) {
    norm_vec[j] = integrate_Simpson(0.0, 0.5 * pconst.pi(), 1000, Omega_vec[j]);
  }

  return;
}



// ##################################################################
// ##################################################################



// Integrand for integral in delta function
double stellar_wind_latdep::f(
    double theta,  // Co-latitude angle (radians)
    double Omega   // Omega (v_rot/v_crit)
)
{
  return sin(theta) * pconst.pow_fast(1.0 - Omega * sin(theta), c_xi);
}



// ##################################################################
// ##################################################################



double stellar_wind_latdep::integrate_Simpson(
    const double min,    ///< lower limit
    const double max,    ///< upper limit
    const long int npt,  ///< Number of points (must be even)
    const double Omega   ///<  Omega
)
{
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  double hh  = (max - min) / npt;
  double ans = 0.0;
  // f(0) lower limit
  ans += f(min, Omega) * sin(min);
  // f(N) upper limit
  ans += f(max, Omega) * sin(max);
  // Intermediate points.
  int wt   = 4;
  double x = 0.0;
  for (long int i = 1; i < npt; i++) {
    x = min + i * hh;
    ans += wt * f(x, Omega) * sin(x);
    wt = 6 - wt;
  }
  // finally multiply by hh/3.0
  ans *= hh / 3.0;
  return ans;
}



// ##################################################################
// ##################################################################



// Terminal wind velocity function
double stellar_wind_latdep::fn_v_inf(
    double theta,  // Co-latitude angle (radians)
    double Omega,  // Omega (v_rot/v_crit)
    double v_inf   // terminal velocity at pole (cm/s)
)
{
  return v_inf * pconst.pow_fast(1.0 - Omega * sin(theta), c_gamma);
}



// ##################################################################
// ##################################################################



// Terminal wind velocity function
double stellar_wind_latdep::interp_v_inf(
    double theta,  // Co-latitude angle (radians)
    double Omega,  // Omega (v_rot/v_crit)
    double v_inf   // terminal velocity at pole (cm/s)
)
{
  vector<double> input(2);
  input[0] = theta;
  input[1] = Omega;
  return v_inf
         * interpolate.root_find_bilinear_vec(
               theta_vec, Omega_vec, vinf_vec, vinf_vsize, input);
  ;
}



// ##################################################################
// ##################################################################



double stellar_wind_latdep::interp_density(
    double r,       // Radius (cm)
    double theta,   // Co-latitude angle (radians)
    double vtheta,  // terminal velocity at co-latitude theta (cm/s)
    double Omega,   // Omega (v_rot/v_crit)
    double mdot,    // Mass loss rate (g/s)
    double md0      // mass-loss rate for non-rotating star (g/s)
)
{
  // Vector for delta input (Omega, Teff)
  vector<double> input(2);
  input[0]   = theta;
  input[1]   = Omega;
  double lat = interpolate.root_find_bilinear_vec(
      theta_vec, Omega_vec, f_vec, f_vsize, input);

  // get normalisation factor A
  double A = 0.0;
  interpolate.root_find_linear_vec(Omega_vec, norm_vec, Omega, A);
  A = (mdot / md0 - 1.0) / A;

  return md0 * (1.0 + A * lat) / (4.0 * pconst.pi() * r * r * vtheta);
}



// ##################################################################
// ##################################################################



void stellar_wind_latdep::set_wind_cell_reference_state(
    class GridBaseClass *grid,
    struct wind_cell *wc,
    const struct wind_source *WS,
    const double eos_gamma  ///< EOS gamma
)
{
  struct stellarwind_params *WP = WS->pars;
  //
  // In this function we set the density, pressure, velocity, and
  // tracer values for the reference state of the cell.  Every
  // timestep the cell-values will be reset to this reference state.
  //
  bool set_rho = true;
  if (wc->dist < 0.75 * WP->radius && ndim > 1) {
    wc->p[RO] = 1.0e-31;
    wc->p[PG] = 1.0e-31;
    set_rho   = false;
  }

  std::array<double, MAX_DIM> pp;
  CI.get_dpos(wc->c, pp);
  double Omega = std::min(0.999, WP->Vrot / WP->Vcrit);

  // calculate terminal wind velocity
  double Vinf = interp_v_inf(wc->theta, Omega, WP->Vinf);

  //
  // 3D geometry: either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
  //
  if (set_rho) {
    wc->p[RO] =
        interp_density(wc->dist, wc->theta, Vinf, Omega, WP->Mdot, WS->Md0);
    //
    // Set pressure based on wind density/temperature at the stellar radius,
    // assuming adiabatic expansion outside Rstar, and that we don't care what
    // the temperature is inside Rstar (because this function will make it
    // hotter than Teff, which is not realistic):
    //
    // rho_star = Mdot/(4.pi.R_star^2.v_inf),
    //   p_star = rho_star.k.T_star/(mu.m_p)
    // So then p(r) = p_star (rho(r)/rho_star)^gamma
    //
    wc->p[PG] = WP->Tstar * pconst.kB() / pconst.m_p();  // taking mu = 1
    wc->p[PG] *= pconst.pow_fast(
        wc->p[RO] * pconst.pow_fast(wc->dist / WP->Rstar, 2), 1.0 - eos_gamma);
    wc->p[PG] *= pconst.pow_fast(wc->p[RO], eos_gamma);
  }

  cell *c  = wc->c;
  double x = 0.0, y = 0.0, z = 0.0, xf = 0.0, yf = 0.0;
  switch (ndim) {
    case 1:
      x = grid->difference_vertex2cell(WP->dpos, c, XX);
      y = 0.0;
      z = 0.0;
      break;
    case 2:
      x = grid->difference_vertex2cell(WP->dpos, c, XX);
      y = grid->difference_vertex2cell(WP->dpos, c, YY);
      z = 0.0;
      break;
    case 3:
      x = grid->difference_vertex2cell(WP->dpos, c, XX);
      y = grid->difference_vertex2cell(WP->dpos, c, YY);
      z = grid->difference_vertex2cell(WP->dpos, c, ZZ);
      break;
    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }

  // Velocities: cell-average values, i.e. values at the
  // centre-of-volume.  Same as for non-latitude-dependent wind, but
  // with Vinf that varies with latitude, calculated above.
  //
  // TODO: for general J vector, what is rotational component.
  // TODO: Add axi-symmetric BC so that VZ,BZ not reflected at
  //       symmetry axis.  Otherwise 2D with rotation won't work.
  switch (ndim) {
    case 1:
      wc->p[VX] = Vinf * x / wc->dist;
      wc->p[VY] = 0.0;
      wc->p[VZ] = 0.0;
      break;

    case 2:
      wc->p[VX] = Vinf * x / wc->dist;
      wc->p[VY] = Vinf * y / wc->dist;
      wc->p[VZ] = 0.0;
      break;

    case 3:
      wc->p[VX] = Vinf * x / wc->dist;
      wc->p[VY] = Vinf * y / wc->dist;
      wc->p[VZ] = Vinf * z / wc->dist;

      // add non-radial component to x/y-dir from rotation.
      // J is hardcoded to be parallel to z-axis
      xf = -WP->Vrot * WP->Rstar * y / pconst.pow_fast(wc->dist, 2);
      yf = WP->Vrot * WP->Rstar * x / pconst.pow_fast(wc->dist, 2);
      wc->p[VX] += xf;
      wc->p[VY] += yf;
      xf /= Vinf * x / wc->dist;  // fraction of x-vel in non-radial dir.
      yf /= Vinf * y / wc->dist;
      break;

    default:
      spdlog::error(
          "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
      break;
  }
  // include stellar space velocity if appropriate
  if (WP->moving_star) {
    for (int v = 0; v < ndim; v++)
      wc->p[VX + v] += WP->velocity[v];
  }

  // Magnetic field: cell-average values, i.e. values at the
  // centre-of-volume.
  // Use a split monopole plus a rotational term adding toroidal
  // component.
  // TODO: for general J vector, what is rotational component.
  if (eqntype == EQMHD || eqntype == EQGLM) {
    double B_s = WP->Bstar / sqrt(4.0 * M_PI);  // code units for B_surf
    double D_s = WP->Rstar / wc->dist;          // 1/d in stellar radii
    double D_2 = D_s * D_s;                     // 1/d^2 in stellar radii
    // this multiplies the toroidal component:
    double beta_B_sint = (WP->Vrot / Vinf) * B_s * D_s;

    switch (ndim) {
      case 1:
        spdlog::error("{}: {}", "1D spherical but MHD?", ndim);
        break;
      case 2:
        // split monopole
        wc->p[BX] = B_s * D_2 * fabs(x) / wc->dist;
        wc->p[BY] = B_s * D_2 / wc->dist;
        wc->p[BY] = (x > 0.0) ? y * wc->p[BY] : -y * wc->p[BY];
        // toroidal component
        beta_B_sint = beta_B_sint * y / wc->dist;
        wc->p[BZ]   = (x > 0.0) ? -beta_B_sint : beta_B_sint;
        break;

      case 3:
        // split monopole along z-axis, parallel to J
        wc->p[BX] = B_s * D_2 / wc->dist;
        wc->p[BX] = (z > 0.0) ? x * wc->p[BX] : -x * wc->p[BX];

        wc->p[BY] = B_s * D_2 / wc->dist;
        wc->p[BY] = (z > 0.0) ? y * wc->p[BY] : -y * wc->p[BY];

        wc->p[BZ] = B_s * D_2 * fabs(z) / wc->dist;

        // toroidal component in x-y plane from rotation, such that
        // we have a Parker spiral, inward winding for z<0 and
        // outwards for z>0.
        beta_B_sint *= sqrt(x * x + y * y) / wc->dist;
        beta_B_sint = (z > 0.0) ? -beta_B_sint : beta_B_sint;

        // NDEBUG
        // modulate strength near the equator by linearly reducing
        // torodial component for |theta|<1 degree from equator
        // See Pogorelov et al (2006,ApJ,644,1299).  This is for
        // testing the code.
        // t = fabs(z)/wc->dist * 180.0 / M_PI; // angle in degrees.
        // if (t < 2.0) beta_B_sint *= 0.5*t;
        // NDEBUG

        wc->p[BX] += -beta_B_sint * y / wc->dist;
        wc->p[BY] += beta_B_sint * x / wc->dist;
        break;

      default:
        spdlog::error(
            "{}: {}", "bad ndim in set_wind_cell_reference_state", ndim);
        break;
    }
  }
  if (eqntype == EQGLM) {
    wc->p[SI] = 0.0;
  }

  //
  // Set the H+ ion fraction so that it goes from y=1 at T>tp to
  // y=1.0e-7 at T<tm, with linear interpolation.  //
  double tm = 1.0e4, tp = 1.5e4;
  if (WS->Hplus >= 0) {
    if (WP->Tstar > tp)
      WP->tr[WS->iHplus] = 1.0;
    else if (WP->Tstar < tm)
      WP->tr[WS->iHplus] = 1.0e-10;
    else
      WP->tr[WS->iHplus] =
          1.0e-10 + (WP->Tstar - tm) * (1.0 - 1.0e-10) / (tp - tm);
  }
  // update tracers
  for (int v = 0; v < ntracer; v++)
    wc->p[ftr + v] = WP->tr[v];


#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // Set the minimum temperature to be Tstar in the wind...
  //
  if (MP) {
    if (MP->Temperature(wc->p, eos_gamma) < Tmin) {
      MP->Set_Temp(wc->p, Tmin, eos_gamma);
    }
  }
  else {
    // appropriate for a neutral medium, He+M mass fraction 0.285.
    wc->p[PG] =
        max(static_cast<double>(wc->p[PG]),
            Tmin * wc->p[RO] * pconst.kB() * 0.78625 / pconst.m_p());
  }
#endif  // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

#ifndef NDEBUG
  spdlog::debug("wc->p : {}", std::vector<double>(wc->p, wc->p + nvar));
#endif
#ifdef TEST_INF
  for (int v = 0; v < nvar; v++) {
    if (!isfinite(wc->p[v])) {
      spdlog::error("NAN in wind source {} {}", v, wc->p[v]);
      spdlog::error("{}: {}", "NAN in wind source", wc->p[v]);
    }
  }
#endif

  return;
}



// ##################################################################
// ##################################################################



int stellar_wind_latdep::add_evolving_source(
    const double t_now,              ///< current time.
    struct stellarwind_params *pars  ///< pointer to wind parameters struct
)
{
  if (pars->type != WINDTYPE_LATDEP) {
    spdlog::error(
        "{}: {}", "Bad wind type for evolving stellar wind (lat-dep)!",
        pars->type);
  }
  //
  // First we will read the file, and see when the source should
  // switch on in the simulation (it may not be needed for a while).
  //
  spdlog::debug(
      "\t\tsw-evo: adding source from file {}", pars->evolving_wind_file);


  //
  // Some properties of the wind source are specific to this module,
  // such as the number of points in
  // the array, and the timings (offset, update frequency).
  // They are stored in local data.
  //
  struct evolving_wind_data *temp = 0;
  double om                       = 0.0;
  temp                            = mem.myalloc(temp, 1);
  int err = read_evolution_file(pars->evolving_wind_file, temp);
  if (err)
    spdlog::error(
        "{}: {}", "couldn't read wind evolution file",
        pars->evolving_wind_file);

  // assign data for v_esc from v_crit.
  for (int i = 0; i < temp->Npt; i++)
    temp->vesc_evo.push_back(temp->vcrt_evo[i] * pconst.sqrt2());

  // optionally enhance Mdot artificially
  // Mdot = Mdot_0 * (1 + {omega-0.5}/0.5*enhance) for omega>0.5
  if (pars->enhance_mdot) {
    spdlog::debug(
        "Enhancing Mdot by factor of (1+{}) at Omega=1", pars->enhance_mdot);
    for (int i = 0; i < temp->Npt; i++) {
      om = temp->vrot_evo[i] / temp->vcrt_evo[i];
      if (om > 0.5)
        temp->Mdot_evo[i] *= 1.0 + 2.0 * pars->enhance_mdot * (om - 0.5);
    }
  }

  // set offsets and scaling for evolutionary time (all in seconds)
  temp->offset        = pars->time_offset / pars->t_scalefactor;
  temp->tstart        = temp->time_evo[0];
  temp->tfinish       = temp->time_evo[temp->Npt - 1];
  temp->update_freq   = pars->update_freq / pars->t_scalefactor;
  temp->t_next_update = max(temp->tstart, t_now);

  spdlog::debug(
      "\t\t tstart={}, next update={}, and tfinish={}", temp->tstart,
      temp->t_next_update, temp->tfinish);

  //
  // Optional time offset between simulation time and evolutionary
  // time.  Also optional scaling.
  //
  for (int i = 0; i < temp->Npt; i++) {
    temp->time_evo[i] += pars->time_offset;
    temp->time_evo[i] /= pars->t_scalefactor;
  }
  //
  // Decide if the wind src is active yet.  If it is, then
  // set up a constant wind source for updating its properties.
  //
  double mdot = 0.0, vinf = 0.0, Twind = 0.0, vrot = 0.0, rstar = 0.0,
         vcrt  = 0.0;
  double Lstar = 0.0, Mstar = 0.0;
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  if (((t_now + temp->update_freq) > temp->tstart
       || pconst.equalD(temp->tstart, t_now))
      && t_now < temp->tfinish) {
    temp->is_active = true;
    // Get the current values for wind properties.
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Teff_evo, t_now, Twind);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Mdot_evo, t_now, mdot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vinf_evo, t_now, vinf);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vrot_evo, t_now, vrot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vcrt_evo, t_now, vcrt);
    interpolate.root_find_linear_vec(temp->time_evo, temp->R_evo, t_now, rstar);
    interpolate.root_find_linear_vec(temp->time_evo, temp->L_evo, t_now, Lstar);
    interpolate.root_find_linear_vec(temp->time_evo, temp->M_evo, t_now, Mstar);

    // get tracer values for elements.
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_H_evo, t_now, xh);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->X_He_evo, t_now, xhe);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_C_evo, t_now, xc);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_N_evo, t_now, xn);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_O_evo, t_now, xo);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_Z_evo, t_now, xz);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_D_evo, t_now, xd);

    spdlog::info("Source is Active");
  }
  else {
    spdlog::debug(
        "WARNING: Source is not yet active: tnow={}, tstart={}. Setting wind source to INACTIVE",
        t_now, temp->tstart);

    temp->is_active = false;
    mdot            = -100.0;
    vinf            = -100.0;
    vrot            = -100.0;
    Twind           = -100.0;
  }
  pars->Tstar = Twind;
  pars->Mdot  = mdot;
  pars->Vinf  = vinf;
  pars->Vrot  = vrot;
  pars->Vcrit = vcrt;
  pars->Rstar = rstar;
  pars->Mass  = Mstar;
  // set tracer values for elements
  set_element_indices(temp);
  if (temp->i_XH >= 0) pars->tr[temp->i_XH] = xh;
  if (temp->i_XHe >= 0) pars->tr[temp->i_XHe] = xhe;
  if (temp->i_XC >= 0) pars->tr[temp->i_XC] = xc;
  if (temp->i_XN >= 0) pars->tr[temp->i_XN] = xn;
  if (temp->i_XO >= 0) pars->tr[temp->i_XO] = xo;
  if (temp->i_XZ >= 0) pars->tr[temp->i_XZ] = xz;
  if (temp->i_XD >= 0) pars->tr[temp->i_XD] = xd;

  // Set B-field of star
  // TODO: Decide how to set this better!  For now pick B=10G at
  //       radius 10 R_sun, and scale with R^-2 for constant flux.
  // One solution: just set it to a constant throughout evolution.
  // double Bstar= 10.0*pconst.pow_fast(10.0*pconst.Rsun()/rstar,2.0);
  // double Bstar= 1.0; // fix to 1 Gauss at the surface, regardless of radius

  // get the mass-loss rate of an equivalent non-rotating star from
  // Brott et al. (2011) recipe
  double Z   = (xc + xn + xo + xz + xd) / 0.014;
  double md0 = Mdot_Brott(
      Lstar / pconst.Lsun(), Mstar / pconst.Msun(), Twind,
      rstar / pconst.Rsun(), Z);

  // Now add source using rotating star version.
  add_rotating_source(pars);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_latdep::add_rotating_source(
    struct stellarwind_params *wp  ///< pointer to wind parameters struct
)
{
  struct wind_source *ws = 0;
  ws                     = mem.myalloc(ws, 1);
  ws->pars               = wp;
  ws->ncell              = 0;
  switch (wp->type) {
    case WINDTYPE_LATDEP:
      spdlog::debug("\tAdding latitude-dependent wind source as id={}", wp->id);
      break;
    default:
      spdlog::error(
          "{}: {}", "What type of source is this?  add a new type?", wp->type);
      break;
  }

  // if using microphysics, find H+ tracer variable, if it exists.
  ws->Hplus  = -1;
  ws->iHplus = -1;
  int hplus  = -1;
  if (MP) {
    hplus = MP->Tr("H1+");
  }
  ws->Hplus = hplus;
  if (hplus >= 0) ws->iHplus = hplus - nvar + ntracer;

  ws->cells_added = false;
  if (!ws->wcells.empty())
    spdlog::error(
        "{}: {}", "wind_source: wcells not empty!", ws->wcells.size());

  //
  // Make sure the source position is compatible with the geometry:
  //
  if (coordsys == COORD_SPH) {
    if (!pconst.equalD(wp->dpos[Rsph], 0.0))
      spdlog::error(
          "{}: {}", "Spherical symmetry but source not at origin!",
          wp->dpos[Rsph]);
  }
  if (coordsys == COORD_CYL && ndim == 2) {
    //
    // Axisymmetry
    //
    if (!pconst.equalD(wp->dpos[Rcyl], 0.0))
      spdlog::error(
          "{}: {}", "Axisymmetry but source not at R=0!", wp->dpos[Rcyl]);
  }

  wlist.push_back(ws);
  nsrc++;

  spdlog::debug(
      "\tAdded wind source id={} to list of {} elements", nsrc - 1, nsrc);

  return wp->id;
}



// ##################################################################
// ##################################################################



void stellar_wind_latdep::update_source(
    class GridBaseClass *grid,
    struct evolving_wind_data *wd,
    const double t_now,
    const double eos_gamma)
{
  //
  // We have a source that needs updating.  If it is not active, and
  // needs activating then we set that.
  //
  struct stellarwind_params *wp = wd->ws->pars;
  if (!wd->is_active) {
    array<double, MAX_DIM> wpos;
    for (int v = 0; v < MAX_DIM; v++)
      wpos[v] = wp->dpos[v];
    spdlog::debug(
        "stellar_wind_latdep::update_source() activating source id={} at Simulation time t=",
        wp->id, t_now);
    // spdlog::debug("Source position : {}", wpos);
    wd->is_active = true;
  }

  if (t_now < wd->tstart) {
    spdlog::error(
        "{}: {}", "Requested updating inactive source", wd->tstart - t_now);
  }

  wd->t_next_update = t_now;  // (We update every timestep now)
  wd->t_next_update = min(wd->t_next_update, wd->tfinish);

  //
  // Now we update Mdot, Vinf, Teff by linear interpolation.
  //
  double mdot = 0.0, vinf = 0.0, Twind = 0.0, vrot = 0.0, rstar = 0.0,
         vcrit = 0.0;
  double Lstar = 0.0, Mstar = 0.0;
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  interpolate.root_find_linear_vec(wd->time_evo, wd->Teff_evo, t_now, Twind);
  interpolate.root_find_linear_vec(wd->time_evo, wd->Mdot_evo, t_now, mdot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vinf_evo, t_now, vinf);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vrot_evo, t_now, vrot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vcrt_evo, t_now, vcrit);
  interpolate.root_find_linear_vec(wd->time_evo, wd->R_evo, t_now, rstar);
  interpolate.root_find_linear_vec(wd->time_evo, wd->L_evo, t_now, Lstar);
  interpolate.root_find_linear_vec(wd->time_evo, wd->M_evo, t_now, Mstar);

  // all in cgs units already.
  wp->Mdot  = mdot;
  wp->Vinf  = vinf;
  wp->Vrot  = vrot;
  wp->Vcrit = vcrit;
  wp->Tstar = Twind;
  wp->Rstar = rstar;

  // get tracer values for elements.
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_H_evo, t_now, xh);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_He_evo, t_now, xhe);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_C_evo, t_now, xc);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_N_evo, t_now, xn);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_O_evo, t_now, xo);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_Z_evo, t_now, xz);
  interpolate.root_find_linear_vec(wd->time_evo, wd->X_D_evo, t_now, xd);

  if (wd->i_XH >= 0) wp->tr[wd->i_XH] = xh;
  if (wd->i_XHe >= 0) wp->tr[wd->i_XHe] = xhe;
  if (wd->i_XC >= 0) wp->tr[wd->i_XC] = xc;
  if (wd->i_XN >= 0) wp->tr[wd->i_XN] = xn;
  if (wd->i_XO >= 0) wp->tr[wd->i_XO] = xo;
  if (wd->i_XZ >= 0) wp->tr[wd->i_XZ] = xz;
  if (wd->i_XD >= 0) wp->tr[wd->i_XD] = xd;

  // Set B-field of star
  // TODO: Decide how to set this better!  For now pick B=10G at
  //       radius 10 R_sun, and scale with R^-2 for constant flux.
  // One solution: leave it constant
  // wd->ws->Bstar= 10.0*pconst.pow_fast(10.0*pconst.Rsun()/rstar,2.0);

  // get the mass-loss rate of an equivalent non-rotating star from
  // Brott et al. (2011) recipe
  double Z    = (xc + xn + xo + xz + xd) / 0.014;
  wd->ws->Md0 = Mdot_Brott(
      Lstar / pconst.Lsun(), Mstar / pconst.Msun(), Twind,
      rstar / pconst.Rsun(), Z);
  wd->ws->Md0 = std::min(wd->ws->Md0, wp->Mdot);  // just in case

  //
  // Now re-assign state vector of each wind-boundary-cell with
  // updated values.
  //
  for (int i = 0; i < wd->ws->ncell; i++) {
    set_wind_cell_reference_state(grid, wd->ws->wcells[i], wd->ws, eos_gamma);
  }

  return;
}



// ##################################################################
// ##################################################################
