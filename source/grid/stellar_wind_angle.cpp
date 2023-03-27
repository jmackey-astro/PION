///
/// \file stellar_wind_angle.cpp
/// \author Jonathan Mackey
/// \date 2017.07.01
///
/// Modifications:
/// - 2017.07.[20-31] RK/JM: Getting it working.
/// - 2017.08.[15-?] RK: reworked to include trilinear interpolation
/// of density, to accommodate for addition of Teff dimension to
/// setup_tables
/// - 2017-08-30 JM: added flag for enhancing wind

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
#endif  // NDEBUG

#include "grid/grid_base_class.h"
#include "grid/stellar_wind_angle.h"
#include "microphysics/microphysics_base.h"
#include <sstream>

#include <cstdio>
#include <vector>
using namespace std;



// ##################################################################
// ##################################################################



// Constructor
stellar_wind_angle::stellar_wind_angle(
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
  stellar_wind_angle::c_gamma = 0.35;
  stellar_wind_angle::c_xi    = xi;
  stellar_wind_angle::c_beta  = -1.0;

  // Number of points in theta, omega and Teff vectors
  stellar_wind_angle::npts_theta = 25;
  stellar_wind_angle::npts_omega = 25;
  stellar_wind_angle::npts_Teff  = 22;

  // Generate interpolating tables
  setup_tables();
}



// ##################################################################
// ##################################################################



// Destructor
stellar_wind_angle::~stellar_wind_angle() {}



// ##################################################################
// ##################################################################



// Generate interpolating tables for wind density function
void stellar_wind_angle::setup_tables()
{
  // Min, mid and max values of theta - 5pts between min and mid, 20pts
  // between mid and max
  double theta_min = 0.1;
  double theta_mid = 60.0;
  double theta_max = 89.9;

  theta_vec.resize(npts_theta);

  for (int k = 0; k < npts_theta; k++) {
    if (k <= 4)
      theta_vec[k] = (theta_min + k * ((theta_mid - theta_min) / 4.0))
                     * (pconst.pi() / 180.0);
    else
      theta_vec[k] =
          (theta_mid + (k - 4) * ((theta_max - theta_mid) / (npts_theta - 5)))
          * (pconst.pi() / 180.0);
  }

  //
  // Set up omega vector
  //

  log_mu_vec.resize(npts_omega);

  // Iterate log(mu) values - evenly spaced
  for (int i = 0; i < npts_omega; i++)
    log_mu_vec[npts_omega - i - 1] = -4.0 + i * (4.0 / (npts_omega - 1));

  omega_vec.resize(npts_omega);

  // Iterate omega values - log spacing - spacing decreases as omega
  // approaches
  // 1
  for (int j = 0; j < npts_omega; j++)
    omega_vec[j] = 1 - pconst.pow_fast(10, log_mu_vec[j]);

  //
  // Set up Teff vector
  //

  // Temperature ranges from Eldridge et al. (2006, MN, 367, 186) (K) + upper
  // and lower Teff limits
  double T0 = 1000.0, T1 = 3600.0, T2 = 6000.0, T3 = 8000.0, T4 = 10000.0,
         T5 = 20000.0, T6 = 22000.0, T7 = 150000.0;

  Teff_vec.resize(npts_Teff);

  // Set values - based on plots for alpha and delta vs. Teff
  for (int i = 0; i < npts_Teff; i++) {
    if (i == 0) Teff_vec[i] = T0;
    if (i == 1) Teff_vec[i] = T1;
    if (2 <= i && i <= 6) Teff_vec[i] = T1 + i * ((T2 - T1) / 6);
    if (i == 7) Teff_vec[i] = T2;
    if (8 <= i && i <= 10) Teff_vec[i] = T2 + (i - 6) * ((T3 - T2) / 4);
    if (i == 11) Teff_vec[i] = T3;
    if (12 <= i && i <= 14) Teff_vec[i] = T3 + (i - 10) * ((T4 - T3) / 4);
    if (i == 15) Teff_vec[i] = T4;
    if (i == 16) Teff_vec[i] = T5;
    if (17 <= i && i <= 19) Teff_vec[i] = T5 + (i - 15) * ((T6 - T5) / 4);
    if (i == 20) Teff_vec[i] = T6;
    if (i == 21) Teff_vec[i] = T7;
  }

  //
  // Write delta table
  //

  delta_vec.resize(npts_omega);
  for (int i = 0; i < npts_omega; i++)
    delta_vec[i].resize(npts_Teff);

  // Set values - delta(omega, Teff)
  for (int i = 0; i < npts_omega; i++) {
    for (int j = 0; j < npts_Teff; j++) {
      delta_vec[i][j] = fn_delta(omega_vec[i], Teff_vec[j]);
#ifndef NDEBUG
      if (!isfinite(delta_vec[i][j])) {
        spdlog::debug(
            "infinite delta!!! {} {} {} {} {}", i, j, omega_vec[i], Teff_vec[j],
            delta_vec[i][j]);
        spdlog::error("{}: {}", "bug", delta_vec[i][j]);
      }
#endif
    }
  }

  //
  // Write alpha table
  //

  alpha_vec.resize(npts_omega);

  for (int i = 0; i < npts_omega; i++) {
    alpha_vec[i].resize(npts_theta);
    for (int j = 0; j < npts_theta; j++) {
      alpha_vec[i][j].resize(npts_Teff);
    }
  }

  // Set values - alpha(omega, theta, Teff)
  for (int i = 0; i < npts_omega; i++) {
    for (int j = 0; j < npts_theta; j++) {
      for (int k = 0; k < npts_Teff; k++) {
        alpha_vec[i][j][k] = fn_alpha(omega_vec[i], theta_vec[j], Teff_vec[k]);
      }
    }
  }

  return;
}



// ##################################################################
// ##################################################################



// Integrand for integral in delta function
double stellar_wind_angle::integrand(
    double theta,  // Co-latitude angle (radians)
    double omega,  // Omega (v_rot/v_crit)
    double Teff    // Teff (K)
)
{
  return fn_alpha(omega, theta, Teff)
         * pconst.pow_fast(1.0 - omega * sin(theta), c_xi) * sin(theta);
}



// ##################################################################
// ##################################################################



double stellar_wind_angle::integrate_Simpson(
    const double min,    ///< lower limit
    const double max,    ///< upper limit
    const long int npt,  ///< Number of points (must be even)
    const double omega,  ///<  omega
    const double Teff    ///< Teff (K)
)
{

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh  = (max - min) / npt;
  double ans = 0.0;
  //
  // f(0) lower limit
  //
  ans += integrand(min, omega, Teff);
  //
  // f(N) upper limit
  //
  ans += integrand(max, omega, Teff);
  //
  // Intermediate points.
  //
  int wt   = 4;
  double x = 0.0;
  for (long int i = 1; i < npt; i++) {
    x = min + i * hh;
    ans += wt * integrand(x, omega, Teff);
    wt = 6 - wt;
  }
  //
  // finally multiply by hh/3.0
  //
  ans *= hh / 3.0;
  return ans;
}



// ##################################################################
// ##################################################################



// Phi' function
double stellar_wind_angle::fn_phi(
    double omega,  // omega (v_rot/v_crit)
    double theta,  // Co-latitude angle (radians)
    double Teff    // Teff (K)
)
{
  double ans = (omega / (22.0 * pconst.sqrt2() * beta(Teff))) * sin(theta)
               * pconst.pow_fast(1.0 - omega * sin(theta), -c_gamma);
  return std::min(ans, 0.5 * pconst.pi() * ONE_MINUS_EPS);
}



// ##################################################################
// ##################################################################



// Alpha function
double stellar_wind_angle::fn_alpha(
    double omega,  // Omega (v_rot/v_crit)
    double theta,  // Co-latitude angle (radians)
    double Teff    // Teff (K)
)
{
  return pconst.pow_fast(
      cos(fn_phi(omega, theta, Teff))
          + pconst.pow_fast(tan(theta), -2.0)
                * (1.0
                   + c_gamma
                         * (omega * sin(theta) / (1.0 - omega * sin(theta))))
                * fn_phi(omega, theta, Teff) * sin(fn_phi(omega, theta, Teff)),
      -1.0);
}  // the cotan term will diverge here if theta = 0.0



// ##################################################################
// ##################################################################



// Delta function
double stellar_wind_angle::fn_delta(
    double omega,  // Omega (v_rot/v_crit)
    double Teff    // Teff (K)
)
{
  return 2.0
         * pconst.pow_fast(
             integrate_Simpson(0.001, pconst.pi() / 2.0, 230, omega, Teff),
             -1.0);
}  // 230 points determined to give sufficient accuracy



// ##################################################################
// ##################################################################



// Terminal wind velocity function
double stellar_wind_angle::fn_v_inf(
    double omega,  // Omega (v_rot/v_crit)
    double v_inf,  // terminal velocity at pole (cm/s)
    double theta   // Co-latitude angle (radians)
)
{

  omega = std::min(omega, 0.999);
  return std::max(
      0.5e5, v_inf * pconst.pow_fast(1.0 - omega * sin(theta), c_gamma));
}



// ##################################################################
// ##################################################################



// Analytic wind density function
double stellar_wind_angle::fn_density(
    double omega,   // Omega (v_rot/v_crit)
    double v_inf,   // terminal velocity at pole (cm/s)
    double mdot,    // Mass loss rate (g/s)
    double radius,  // Radius (cm)
    double theta,   // Co-latitude angle (radians)
    double Teff     // Teff (K)
)
{
  return (mdot * fn_alpha(omega, theta, Teff) * fn_delta(omega, Teff)
          * pconst.pow_fast(1.0 - omega * sin(theta), c_xi))
         / (8.0 * pconst.pi() * pconst.pow_fast(radius, 2.0)
            * fn_v_inf(omega, v_inf, theta));
}



// ##################################################################
// ##################################################################



// Interpolated wind density function
double stellar_wind_angle::fn_density_interp(
    double omega,   // Omega (v_rot/v_crit)
    double v_inf,   // terminal velocity at pole (cm/s)
    double mdot,    // Mass loss rate (g/s)
    double radius,  // Radius (cm)
    double theta,   // Co-latitude angle (radians)
    double Teff     // Teff (K)
)
{
  omega = std::min(omega, 0.999);
  //
  // Use tables to interpolate the value of delta
  //
  double delta_interp;

  // Vector for delta interpolation vector sizes
  vector<size_t> delta_vec_size(2);
  delta_vec_size[0] = npts_omega;
  delta_vec_size[1] = npts_Teff;

  // Vector for delta input (omega, Teff)
  vector<double> delta_input(2);
  delta_input[0] = omega;
  delta_input[1] = Teff;

  delta_interp = root_find_bilinear_vec(
      omega_vec, Teff_vec, delta_vec, delta_vec_size, delta_input);

  // Use tables to interpolate the value of alpha
  double alpha_interp;

  // Vector for delta interpolation vector sizes
  vector<size_t> alpha_vec_size(3);
  alpha_vec_size[0] = npts_omega;
  alpha_vec_size[1] = npts_theta;
  alpha_vec_size[2] = npts_Teff;

  // Vector for delta input (omega, Teff)
  vector<double> alpha_input(3);
  alpha_input[0] = omega;
  alpha_input[1] = theta;
  alpha_input[2] = Teff;

  alpha_interp = root_find_trilinear_vec(
      omega_vec, theta_vec, Teff_vec, alpha_vec, alpha_vec_size, alpha_input);

  // Return interpolated density
  double result =
      (mdot * alpha_interp * delta_interp
       * pconst.pow_fast(1.0 - omega * sin(theta), c_xi));
  result /=
      (8.0 * pconst.pi() * pconst.pow_fast(radius, 2.0)
       * fn_v_inf(omega, v_inf, theta));

#ifndef NDEBUG
  if (!isfinite(result)) {
    spdlog::debug(
        "{} {} {} {} {}\n{} {} {} {} {} {} {}", delta_interp, alpha_interp,
        result, fn_v_inf(omega, v_inf, theta), omega, mdot, alpha_interp,
        delta_interp, pconst.pow_fast(1.0 - omega * sin(theta), c_xi),
        pconst.pow_fast(radius, 2.0), fn_v_inf(omega, v_inf, theta), radius);
  }
#endif
  return result;
}



// ##################################################################
// ##################################################################



void stellar_wind_angle::set_wind_cell_reference_state(
    class GridBaseClass &grid,
    struct wind_cell &wc,
    const struct wind_source &WS,
    const double eos_gamma  ///< EOS gamma
)
{
  struct stellarwind_params *WP = WS.pars;
  //
  // In this function we set the density, pressure, velocity, and tracer
  // values for the reference state of the cell.  Every timestep the
  // cell-values will be reset to this reference state.
  //
  if (wc.dist < WP->Rstar && ndim > 1) {
    set_stellar_interior_values(wc, *WP);
    return;
  }

  spdlog::debug(
      "id {}: vinf {:12.6e}, vrot {:12.6e}, vcrit {:12.6e}", WS.pars->id,
      WP->Vinf, WP->Vrot, WP->Vcrit);

  //
  // 3D geometry: either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
  //
  wc.p[RO] = fn_density_interp(
      std::min(0.9999, WP->Vrot / WP->Vcrit), WP->Vinf, WP->Mdot, wc.dist,
      wc.theta, WP->Tstar);
  //
  // Set pressure based on wind density/temperature at the stellar radius,
  // assuming adiabatic expansion outside Rstar
  //
  wc.p[PG] = WP->Tstar * pconst.kB() / pconst.m_p();  // taking mu = 1
  wc.p[PG] *= pconst.pow_fast(
      fn_density_interp(
          std::min(0.9999, WP->Vrot / WP->Vcrit), WP->Vinf, WP->Mdot, WP->Rstar,
          wc.theta, WP->Tstar),
      1.0 - eos_gamma);
  wc.p[PG] *= pconst.pow_fast(wc.p[RO], eos_gamma);

  // calculate terminal wind velocity
  double Vinf =
      fn_v_inf(std::min(0.9999, WP->Vrot / WP->Vcrit), WP->Vinf, wc.theta);

  cell *c = wc.c;
  std::array<double, MAX_DIM> r;
  set_wind_cell_offset(grid, wc, WS, r);
  set_wind_cell_velocity_components_vinf(r, wc, *WP, Vinf, WP->Vrot);
  if (eqntype == EQMHD || eqntype == EQGLM) {
    set_wind_cell_B_components_vinf(r, wc, *WP, Vinf, WP->Vrot);
  }
  set_wind_cell_tracers(wc, WS);

#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // Set the minimum temperature to be Tstar in the wind...
  //
  if (MP) {
    if (MP->Temperature(wc.p, eos_gamma) < Tmin) {
      MP->Set_Temp(wc.p, Tmin, eos_gamma);
    }
  }
  else {
    // appropriate for a neutral medium, He+M mass fraction 0.285.
    wc.p[PG] =
        max(static_cast<double>(wc.p[PG]),
            Tmin * wc.p[RO] * pconst.kB() * 0.78625 / pconst.m_p());
  }
#endif  // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

#ifndef NDEBUG
  spdlog::debug("wc.p : {}", std::vector<double>(wc.p, wc.p + ndim));
#endif
#ifdef TEST_INF
  for (int v = 0; v < nvar; v++) {
    if (!isfinite(wc.p[v])) {
      spdlog::error("NAN in wind source {} {}", v, wc.p[v]);
      spdlog::error("{}: {}", "NAN in wind source", wc.p[v]);
    }
  }
#endif

  return;
}



// ##################################################################
// ##################################################################



int stellar_wind_angle::add_evolving_source(
    const double t_now,  ///< current sim time, to see if src is active.
    struct stellarwind_params *pars  ///< pointer to wind parameters struct
)
{
  if (pars->type != WINDTYPE_ANGLE) {
    spdlog::error(
        "{}: {}", "Bad wind type for evolving stellar wind (rotating star)!",
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
  temp                            = mem.myalloc(temp, 1);
  int err = read_evolution_file(pars->evolving_wind_file, temp);
  if (err) {
    spdlog::error(
        "couldn't read wind evolution file: {}", pars->evolving_wind_file);
  }

  // assign data for v_esc from v_crit.
  for (int i = 0; i < temp->Npt; i++)
    temp->vesc_evo.push_back(temp->vcrt_evo[i] * pconst.sqrt2());

  // set offsets and scaling for evolutionary time (all in seconds)
  temp->offset        = pars->time_offset / pars->t_scalefactor;
  temp->tstart        = temp->time_evo[0];
  temp->tfinish       = temp->time_evo[temp->Npt - 1];
  temp->update_freq   = pars->update_freq / pars->t_scalefactor;
  temp->t_next_update = max(temp->tstart, t_now);

  spdlog::debug(
      "\t\t tstart={} "
      ", next update={}"
      ", and tfinish={}",
      temp->tstart, temp->t_next_update, temp->tfinish);

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
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  if (((t_now + temp->update_freq) > temp->tstart
       || pconst.equalD(temp->tstart, t_now))
      && t_now < temp->tfinish) {
    temp->is_active = true;
    // Get the current values for wind properties.
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Teff_evo, t_now, pars->Tstar);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->Mdot_evo, t_now, pars->Mdot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vinf_evo, t_now, pars->Vinf);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vrot_evo, t_now, pars->Vrot);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->vcrt_evo, t_now, pars->Vcrit);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->R_evo, t_now, pars->Rstar);

    // get tracer values for elements.
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_H_evo, t_now, xh);
    interpolate.root_find_linear_vec(
        temp->time_evo, temp->X_He_evo, t_now, xhe);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_C_evo, t_now, xc);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_N_evo, t_now, xn);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_O_evo, t_now, xo);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_Z_evo, t_now, xz);
    interpolate.root_find_linear_vec(temp->time_evo, temp->X_D_evo, t_now, xd);

    spdlog::info("Source is Activ");
  }
  else {
    spdlog::warn(
        "Source is not yet active: tnow={}"
        ", tstart={}"
        ". Setting wind source to INACTIVE",
        t_now, temp->tstart);

    temp->is_active = false;
  }

  // set tracer values for elements
  set_element_indices(temp);
  if (temp->i_XH >= 0) pars->tr[temp->i_XH] = xh;
  if (temp->i_XHe >= 0) pars->tr[temp->i_XHe] = xhe;
  if (temp->i_XC >= 0) pars->tr[temp->i_XC] = xc;
  if (temp->i_XN >= 0) pars->tr[temp->i_XN] = xn;
  if (temp->i_XO >= 0) pars->tr[temp->i_XO] = xo;
  if (temp->i_XZ >= 0) pars->tr[temp->i_XZ] = xz;
  if (temp->i_XD >= 0) pars->tr[temp->i_XD] = xd;

  //
  // Now add source using rotating star version.
  //
  add_rotating_source(pars);
  temp->ws = wlist.back();
  wdata_evol.push_back(temp);
  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_angle::add_rotating_source(
    struct stellarwind_params *wp  ///< pointer to wind parameters struct
)
{
  struct wind_source *ws = 0;
  ws                     = mem.myalloc(ws, 1);
  ws->pars               = wp;
  ws->ncell              = 0;
  switch (wp->type) {
    case WINDTYPE_ANGLE:
      spdlog::debug("\tAdding rotating wind source as id={}", wp->id);
      break;
    default:
      spdlog::error("What type of source is this? {}", wp->type);
      exit(1);
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
  if (!ws->wcells.empty()) {
    spdlog::error("wind_source: wcells not empty :( {}", ws->wcells.size());
    exit(1);
  }

  //
  // Make sure the source position is compatible with the geometry:
  //
  if (coordsys == COORD_SPH) {
    if (!pconst.equalD(wp->dpos[Rsph], 0.0)) {
      spdlog::error(
          "Spherical symmetry but source not at origin {}", wp->dpos[Rsph]);
      exit(1);
    }
  }
  if (coordsys == COORD_CYL && ndim == 2) {
    //
    // Axisymmetry
    //
    if (!pconst.equalD(wp->dpos[Rcyl], 0.0)) {
      spdlog::error("Axisymmetry but source not at R=0 : {}", wp->dpos[Rcyl]);
      exit(1);
    }
  }

  wlist.push_back(ws);
  nsrc++;

  spdlog::debug(
      "Added wind source id={} to list of {} elements", nsrc - 1, nsrc);

  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_angle::update_source(
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
        "stellar_wind_angle::update_source() activating source id={}"
        " at Simulation time t={}",
        wp->id, t_now);
    // spdlog::debug("Source position : {}", wd->ws->dpos);
    wd->is_active = true;
  }

  if (t_now < wd->tstart) {
    spdlog::warn(
        "{}: {}", "Updating source, not yet active!", wd->tstart - t_now);
    return 0;
  }
  else if (t_now >= wd->tfinish) {
    spdlog::warn(
        "{}: {}", "Updating source: source no longer active!",
        wd->tstart - t_now);
    return 1;
  }

  wd->t_next_update = t_now;  // (We update every timestep now)
  wd->t_next_update = min(wd->t_next_update, wd->tfinish);

  //
  // Now we update Mdot, Vinf, Teff by linear interpolation.
  //
  double mdot = 0.0, vinf = 0.0, Twind = 0.0, vrot = 0.0, rstar = 0.0,
         vcrit = 0.0;
  double xh = 0.0, xhe = 0.0, xc = 0.0, xn = 0.0, xo = 0.0, xz = 0.0, xd = 0.0;

  interpolate.root_find_linear_vec(wd->time_evo, wd->Teff_evo, t_now, Twind);
  interpolate.root_find_linear_vec(wd->time_evo, wd->Mdot_evo, t_now, mdot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vinf_evo, t_now, vinf);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vrot_evo, t_now, vrot);
  interpolate.root_find_linear_vec(wd->time_evo, wd->vcrt_evo, t_now, vcrit);
  interpolate.root_find_linear_vec(wd->time_evo, wd->R_evo, t_now, rstar);

  // all in cgs units already.
  wp->Mdot  = mdot;
  wp->Vinf  = vinf;
  wp->Vrot  = vrot;
  wp->Vcrit = vcrit;
  wp->Tstar = std::min(Twind, Teff_vec.back());
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

  //
  // Now re-assign state vector of each wind-boundary-cell with
  // updated values.
  //
  for (int i = 0; i < wd->ws->ncell; i++) {
    set_wind_cell_reference_state(
        *grid, *(wd->ws->wcells[i]), *(wd->ws), eos_gamma);
  }

  return 0;
}



// ##################################################################
// ##################################################################



int stellar_wind_angle::set_cell_values(
    class GridBaseClass *grid,
    const int id,       ///< src id
    const double t_now  ///< simulation time
)
{
  int err = 0;
  if (id < 0 || id >= nsrc) {
    spdlog::error("{}: {}", "bad src id", id);
  }

  // check for evolving source
  struct wind_source *ws = wlist[id];
  if (ws->pars->evolving_wind_file != "NOFILE") {
    struct evolving_wind_data *wd = wdata_evol[id];
    // Check if the source needs to be updated, and if so, update it
    if (t_now >= wd->t_next_update) {
      // NB: This function has hardcoded EOS Gamma to 5/3
      err = update_source(grid, wd, t_now, 5. / 3.);
      if (err) {
        spdlog::error("update_source() returned error {}", err);
        return err;
      }
    }
    if (!wd->is_active) {
      return 0;
    }
  }
  else {
    // if source is not evolving, set the cell state.
    // spdlog::info("setting wind cell reference state src {}",id);
    // need to set Vcrit because this is not a parameter.
    ws->pars->Vcrit = sqrt(pconst.G() * ws->pars->Mass / ws->pars->Rstar);
    spdlog::debug(
        "Vcrt {:12.6e}, G {:12.6e}, M {:12.6e}, R {:12.6e}", ws->pars->Vcrit,
        pconst.G(), ws->pars->Mass, ws->pars->Rstar);

    for (int i = 0; i < ws->ncell; i++) {
      set_wind_cell_reference_state(*grid, *(ws->wcells[i]), *ws, 5. / 3.);
    }
  }

  // Update the wind cell values using the constant wind functions.
  // spdlog::info("setting cell values src {}",id);
  err += stellar_wind::set_cell_values(grid, id, t_now);

  return err;
}



// ##################################################################
// ##################################################################
