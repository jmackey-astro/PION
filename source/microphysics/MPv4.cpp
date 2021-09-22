///
/// \file MPv4.cpp
/// \author Jonathan Mackey
/// \date 2011.10.12
///
/// Description:
/// THIS IS LEGACY CODE USED ONLY FOR TESTING (Mackey,2012,A&A).
/// This class is an update on the microphysics class used for JMs
/// thesis (MPV1).  It uses the same implicit raytracing scheme, so
/// photon transmission through the cell is integrated with the rate
/// and energy equations to obtain a time-averaged value.
/// The main advantage here is that we can use multi-frequency
/// photoionisation rates which depend on the optical depth.
/// It was used in Mackey (2012,A&A) but found to be less efficient
/// than MPv3.
///
/// Modifications:
/// - 2011.10.13 JM: Debugging.
/// - 2011.10.17 JM: Debugging. (2011.10.22 also).
/// - 2013.02.14 JM: Tidied up file.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifdef LEGACY_CODE

//#define MPV3_DEBUG

#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "microphysics/MPv4.h"

using namespace std;

//
// The timestep-limiting is set by ifdef in
//  source/defines/functionality_flags.h
// DT05 is recommended for MPv4 (see Mackey,2012,
// A&A,539,A147), but DT02 is fine for monochromatic radiation.
//
#if MPV4_DTLIMIT == 0
#define DTFRAC 1.00
#elif MPV4_DTLIMIT == 1
#define DTFRAC 0.30
#elif MPV4_DTLIMIT == 2
#define DTFRAC 0.10
#elif MPV4_DTLIMIT == 3
#define DTFRAC 0.03
#elif MPV4_DTLIMIT == 4
#define DTFRAC 0.003
#elif MPV4_DTLIMIT == 5
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.5
#elif MPV4_DTLIMIT == 6
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.25
#elif MPV4_DTLIMIT == 7
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.125
#elif MPV4_DTLIMIT == 8
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.0625
#elif MPV4_DTLIMIT == 9
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.5
#elif MPV4_DTLIMIT == 10
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.25
#elif MPV4_DTLIMIT == 11
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.125
#elif MPV4_DTLIMIT == 12
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.0625
#else
#error "Must define a DTXX timestep limit in mp_implicit"
#endif

// ##################################################################
// ##################################################################

void MPv4::get_error_tolerances(
    double *reltol,  ///< relative error tolerance.
    double atol[]    ///< absolute error tolerances
)
{
  MPv3::get_error_tolerances(reltol, atol);
  //
  // Now the absolute tolerance for the integrated optical depth
  // This is exp(-dTau)*delta-t, so it can be tiny.  Basically
  // for every addition of 2 in optical depth, this gets smaller
  // by an order of magnitude.  So try 1e-40 to start with.
  //
  atol[lv_dtau] = 1.0e-30;
  return;
}

// ##################################################################
// ##################################################################

void MPv4::get_problem_size(
    int *ne,  ///< number of equations
    int *np   ///< number of parameters in user_data vector.
)
{
  *ne = N_equations;
  *np = N_extradata;
  return;
}

// ##################################################################
// ##################################################################

void MPv4::setup_local_vectors()
{
  nvl         = 3;  // three local variables to integrate
  N_extradata = 0;
  N_equations = 3;

  //
  // MPv3 has already set these up with 2 variables each, so we
  // delete and re-init them.
  //
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
  y_in    = N_VNew_Serial(N_equations);
  y_out   = N_VNew_Serial(N_equations);
  lv_H0   = 0;  // x(H0) is the first element in the array
  lv_eint = 1;  // E_{int} is the second element.
  lv_dtau = 2;  // integrated photon transmission fraction.
  return;
}

// ##################################################################
// ##################################################################

int MPv4::convert_prim2local(
    const double *p_in,  ///< primitive vector from grid cell (length nv_prim)
    double *p_local)
{
  //
  // First set x_Hp and E_int using the explicit solver
  //
  int err = MPv3::convert_prim2local(p_in, p_local);

  //
  // Set the integrated optical depth through the cell to zero.
  //
  p_local[lv_dtau] = 0.0;

  return err;
}

// ##################################################################
// ##################################################################

double MPv4::get_timeaveraged_rhodx(
    const double edtau,  ///< this is int(exp(-dtau)dt)
    const double dt      ///< this is dt
)
{
  //
  // Convert the integrated optical depth to a projected mass density
  // of ionised gas. (code taken mostly from microphysics.cc)
  // Convert integrated exp(-tau) into a time-averaged optical depth
  // to return to the code.
  // Have to check that it is finite and positive, since log is
  // sensitive to roundoff errors.
  //
  double deltau = -log(edtau / dt);
#ifdef RT_TESTING
  cout << "**** MPv4:  INT(exp(-dtau),dt) = " << edtau << " <deltau>=" << deltau
       << ", dt=" << dt << "\n";
#endif
  if (!isfinite(deltau)) {
    // this can happen if the ray is very attenuated, and we are taking
    // -log(0.0)
#ifdef MPV3_DEBUG
    cout << " INT(exp(-dtau),dt) = " << edtau << " <deltau>=" << deltau
         << ", resetting\n";
#endif                // MPV3_DEBUG
    deltau = 1000.0;  // This should be big enough -- exp(-1000)=10^(-440)!!
  }
  if (deltau < 0.0) {
    if (deltau < -JM_RELTOL) {
      cout << "==EEEE== MPv4: <tau> significantly negative: " << deltau;
      cout << " setting to zero. edtau/dt-1=" << edtau / dt - 1.0 << "\n";
      rep.error("negative Tau", deltau);
    }
    deltau = 0.0;
  }
  //
  // Now convert it to a projected neutral mass density Sigma(H0) by
  // multiplying by mu*m_p/crosssection()
  //
  deltau *= mean_mass_per_H / Hi_monochromatic_photo_ion_xsection(2.178721e-11);
#ifdef RT_TESTING
  cout << "**** MPv4:  INT(exp(-dtau),dt) = " << edtau << " <deltau>=" << deltau
       << ", dt=" << dt << "\n";
#endif
  return deltau;
}

// ##################################################################
// ##################################################################

int MPv4::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double *         ///< extra user-data vector (UNUSED)
)
{
  //
  // First get x_Hp_dot and E_int_dot
  //
  int err = MPv3::ydot(0, y_now, y_dot, 0);

  //
  // exp(-delta_Tau) = exp(-nH*(1-x)*sigma*ds), where we evaluate the
  // cross section at 13.6 eV = 2.1787e-11 ergs, the ionisation
  // energy of H.
  //
  NV_Ith_S(y_dot, lv_dtau) =
      exp(-mpv_nH * NV_Ith_S(y_now, lv_H0) * mpv_delta_S
          * Hi_monochromatic_photo_ion_xsection(2.178721e-11));

  return err;
}

// ##################################################################
// ##################################################################

MPv4::MPv4(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr    ///< Number of tracer variables in state vector.
    const std::string *tracers,   ///< List of what the tracer variables mean.
    struct which_physics *ephys,  ///< extra physics stuff.
    struct rad_sources *rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    MPv3(nd, csys, nv, ntr, tracers, ephys, rsrcs, g)
{
  //
  // All of the setup is in the explicit solver; the only changes
  // here are that I have re-implemented some of the functions which
  // the explicit constructor calls.
  //
  cout << "MPv4::MPv4() constructor.\n";
  setup_local_vectors();
  //
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  setup_cvode_solver_without_Jacobian();
  return;
}

// ##################################################################
// ##################################################################

MPv4::~MPv4()
{
  cout << "MPv4::MPv4() destructor.\n";
  return;
}

// ##################################################################
// ##################################################################

int MPv4::TimeUpdateMP_RTnew(
    const double *p_in,  ///< Primitive Vector to be updated.
    const int N_heat,    ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    double *p_out,    ///< Destination Vector for updated values
                      ///< (can be same as first Vector.
    const double dt,  ///< Time Step to advance by.
    const double,     ///< EOS gamma.
    const int,        ///< Switch for what type of integration to use.
                      ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
    double *dsigma    ///< Time-averaged projected neutral mass density, passed
                      ///< back to raytracer.
)
{
  //
  // This does the integration, and leaves int(exp(-dtau)dt) in the
  // variable NV_Ith_S(y_out,lv_dtau).
  //
  int err = MPv3::TimeUpdateMP_RTnew(
      p_in, N_heat, heat_src, N_ion, ion_src, p_out, dt, 0, 0, dsigma);
  //
  // Now we get the time-averaged column density (in g/cm2) to return
  // to the raytracer
  //
  *dsigma = get_timeaveraged_rhodx(NV_Ith_S(y_out, lv_dtau), dt);
  return err;
}

// ##################################################################
// ##################################################################

double MPv4::timescales(
    const double *p_in,  ///< Current cell state vector.
    const double g,      ///< EOS gamma.
    const bool tc,       ///< set to 'true' if including cooling time.
    const bool tr,       ///< set to 'true' if including recombination time.
    const bool ti        ///< set to 'true' if including photo-ionsation time.
)
{
#ifdef MPV3_DEBUG
  if (RS->Nsources != 0) {
    cout << "WARNING: MPv3::timescales() using non-RT version!\n";
  }
#endif  // MPV3_DEBUG
  int err = 0;
  //
  // First convert to local variables.
  //
  double P[nvl];
  err = convert_prim2local(p_in, P);
  if (err) {
    rep.error("Bad input state to MPv3::timescales()", err);
  }
  NV_Ith_S(y_in, lv_H0)   = P[lv_H0];
  NV_Ith_S(y_in, lv_eint) = P[lv_eint];

  //
  // Next set the (dummy) radiation properties of the current cell.
  //
  std::vector<struct rt_source_data> temp;
  setup_radiation_source_parameters(p_in, P, 0, temp, 0, temp);
  temp.clear();

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in MPv3::timescales()", err);
  }
  double T = get_temperature(
      mpv_nH, NV_Ith_S(y_in, lv_eint), 1.0 - NV_Ith_S(y_in, lv_H0));

  //
  // And finally get the smallest timescale over which things are
  // varying.
  //
  double tmin = HUGEVALUE;

  if (tc && T > EP->MinTemperature) {
#ifdef RT_TESTING
    cout << "timestep limiting: cooling time=";
    cout << NV_Ith_S(y_in, lv_eint)
                / (fabs(NV_Ith_S(y_out, lv_eint)) + TINYVALUE)
         << "\n";
#endif
    tmin =
        min(tmin, DTFRAC * NV_Ith_S(y_in, lv_eint)
                      / (fabs(NV_Ith_S(y_out, lv_eint)) + TINYVALUE));
  }

  if (tr) {
    //
    // The recombination time is then the reciprocal of the rate.
    // to give a total MP step of 0.1*t_rec (0.1/RRR)
    //
#ifdef RT_TEST_PROBS
    double rate = Hii_rad_recomb_rate(1.0e4) * mpv_nH;
#else
    double rate = Hii_rad_recomb_rate(T) * (1.0 - NV_Ith_S(y_in, lv_H0))
                  * mpv_nH * JM_NELEC;
#endif
    tmin = min(tmin, DTFRAC / (fabs(rate) + TINYVALUE));

#ifdef RT_TESTING
    cout << "timestep limiting: recomb time=";
    cout << 1.0 / (fabs(rate) + TINYVALUE) << "\n";
#endif
  }

  if (ti) {
    rep.error("ionisation time not possible in timescales()", ti);
  }

  return tmin;
}

// ##################################################################
// ##################################################################

///
/// This returns the minimum timescale of all microphysical
/// processes, including reaction times for each species and the
/// total heating/cooling time for the gas.  It requires the
/// radiation field as an input, so it has substantially greater
/// capability than the other timescales function.
/// Default setting is DT05, which limits by 1.0*E/Edot and 1.0/ydot
///
double MPv4::timescales_RT(
    const double *p_in,  ///< Current cell state vector.
    const int N_heat,    ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    const double  ///< EOS gamma.
)
{
  int err = 0;
  //
  // First convert to local variables.
  //
  double P[nvl];
  err = convert_prim2local(p_in, P);
  if (err) {
    rep.error("Bad input state to MPv4::timescales_RT()", err);
  }
  NV_Ith_S(y_in, lv_H0)   = P[lv_H0];
  NV_Ith_S(y_in, lv_eint) = P[lv_eint];
  NV_Ith_S(y_in, lv_dtau) = P[lv_dtau];

  //
  // Next set the radiation properties of the current cell.
  //
  setup_radiation_source_parameters(p_in, P, N_heat, heat_src, N_ion, ion_src);

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in MPv3::timescales_RT()", err);
  }

  //
  // And finally get the smallest timescale over which things are
  // varying.
  //
  double t = HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on
  // the energy change (heating timescale) is not required for
  // accurately tracking ionisation fronts (although it may be needed
  // for cooling!).
  // For testing purposes there are ifdefs to allow the code to use a
  // relative change in neutral fraction and/or the relative change
  // in energy as the timestep criterion, rather than the default of
  // absolute change in neutral fraction.
  //
#ifdef USE_RELATIVE_NEUFRAC_DTLIMIT
  t =
      min(t, DTFRAC * max(5.0e-2, NV_Ith_S(y_in, lv_H0))
                 / (fabs(NV_Ith_S(y_out, lv_H0)) + TINYVALUE));
  // cout <<"using neutral fraction.\n";
#else
  t = min(t, DTFRAC / (fabs(NV_Ith_S(y_out, lv_H0)) + TINYVALUE));
  // cout <<"limit by dx: dt="<<DTFRAC/(fabs(NV_Ith_S(y_out,
  // lv_H0))+TINYVALUE)<<"\n";
#endif

#ifdef MPV3_DEBUG
  if (t < 3.16e9) {
    cout << "MP timescales: xdot=" << NV_Ith_S(y_out, lv_H0);
    cout << ", Edot=" << NV_Ith_S(y_out, lv_eint) << " t_x=" << t;
  }
#endif  // MPV3_DEBUG

#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  //
  // Now cooling/heating time to limit to X% change in energy).
  //
  t = min(
      t, DTFRAC * P[lv_eint] / (fabs(NV_Ith_S(y_out, lv_eint)) + TINYVALUE));
  // cout <<"using fractional energy change.\n";
  // cout <<"limit by dE: dt="<<DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out,
  // lv_eint))+TINYVALUE)<<"\n";
#endif

#ifdef MPV3_DEBUG
  if (t < 3.16e9) {
    cout << " and min(t_x,t_e)=" << t << ",  ";
    rep.printVec("P[1-x,E]", P, nvl);
  }
#endif  // MPV3_DEBUG

  //
  // Change in Tau should be <N (here <30) per step.
  //
  t =
      min(t, 30.0
                 / (fabs(
                        NV_Ith_S(y_out, lv_H0) * mpv_nH * mpv_delta_S
                        * Hi_monochromatic_photo_ion_xsection(2.178721e-11))
                    + TINYVALUE));

  return t;
}

// ##################################################################
// ##################################################################

#endif  // LEGACY_CODE
