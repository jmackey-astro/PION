///
/// \file MPv3.cpp
/// \author Jonathan Mackey
/// \date 2011.10.06
///
/// Description:
/// This class is an update on the microphysics class used for JM's
/// thesis.  It uses an explicit raytracing scheme, so the photon
/// transmission through the cell does not have to be integrated with
/// the rate and energy equations.
///
/// - It uses multi-frequency photoionisation including spectral
///   hardening with optical depth, using the method outlined in
///   Frank & Mellema (1994,A&A,289,937), and updated by Mellema et
///   al. (2006,NewA,11,374).
/// - It uses the Hummer (1994,MN,268,109) rates for radiative
///   recombination and its associated cooling, and Bremsstrahlung
///   cooling.
/// - For collisional ionisation the function of Voronov
///   (1997,ADNDT,65,1) is used.
/// - Collisional excitation of neutral H, table from Raga, Mellema,
///   Lundqvist (1997,ApJS,109,517), using data from Aggarwal (1993).
/// - Then I use various formulae from Henney et al. (2009,MN,398,
///  157) for cooling due to collisional excitation of photoionised
///  O,C,N (eq.A9), collisional excitation of neutral metals
///  (eq.A10), and the Wiersma+2009 CIE metals-only curve.  I take
///  the max of the WSS09 curve and the Henney et al. functions, to
///  avoid double counting.  For neutral gas I use heating and
///  cooling rates from Wolfire et al. (2003).
/// - Photoheating from ionisation is discussed above.  Cosmic ray
///  heating will use a standard value, X-ray heating is ignored.
///  UV heating due to the interstellar radiation field (ISRF) is
///  according to the optical depth from the edge of the domain to
///  the cell in question, using e.g. HEA09 eq.A3, if requested.  UV
///  heating from the star uses the same equation, but with the
///  optical depth from the source (using a total H-nucleon column
///  density).
///
/// The integration method uses the CVODES solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from
///  https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing (i.e. implicit) with Newton
/// iteration.
///
/// Electrons, ions, ion fraction:  This modules makes the (crude)
/// approximation that He is an identical atom to H -- that it is
/// only ever singly ionised and its ionisation and recombination
/// rates are identical to H.  But because most photons are below the
/// ionisation potential of He0 the opacity of He0 atoms is not
/// counted in the optical depth of the medium to ionising photons.
/// This is not quite self-consistent, but it is a better
/// approximation than assuming He absorbs 14eV photons.
///
/// Modifications:
/// - 2011.03.29 JM: Wrote class MPv3() functions.
/// - 2011.03.31 JM: Added solid angle vector for diffuse radiation.
///    Finished coding, fixed a lot of bugs, need to test it now.
/// - 2011.04.12 JM: Added some 'todo's.
/// - 2011.04.14 JM: Fixed bugs; I am testing the non-RT part now.
///    Looks good so far.
/// - 2011.04.15 JM: Bugfixes
/// - 2011.04.17 JM: Added ifdefs for RT_TESTING to omit various
///    processes in dYdt(). Debugging.  Added ifdef to integrate
///    neutral fraction instead of ion fraction.
/// - 2011.04.18 JM: Bugfixes.
/// - 2011.04.22 JM: Debugged UV/IR heating.  Disabled X-ray heating.
/// - 2011.05.02 JM: Updated convert_prim2local to correct negative
///    pressure/ion fraction inputs rather than bugging out.
/// - 2011.05.02 JM: Added set_multifreq_source_properties() function
/// - 2011.05.04 JM: Added a discretised multifreq photoionisation
///    rate with an approximation for dtau<<1.  Fixed bugs,
///    simplified code.
/// - 2011.05.06 JM: Set cooling limiting as T approaches
///    SimPM.EP.MinTemperature so we don't cool to too low a
///    temperature.
/// - 2011.05.10 JM: Output cooling rates only if myrank==0 for
///    parallel (so processes don't fight over the file and slow
///    down the code (by a lot!)).
/// - 2011.05.25 JM: Fixed a bad bug -- UV heating was not kicking in
///    unless there was a diffuse sources, so the pt-src UV flux was
///    useless up to now. Got rid of "int-ion-frac" dydt function.
///    Too confusing to have two functions.  If for some reason I
///    want to integrate the ion fraction in future I can go back to
///    an older revision.
/// - 2011.06.21 JM: Fixed bug in UV heating where Tau(FUV) was
///    incorrectly calculated.  Simplified timestep-limit calculation
///    for xdot, since 2nd order integration is more accurate and so
///    requires less tuning and allows larger timesteps.
/// - 2011.07.03 JM: Changed how Helium is treated appproximately (it
///    is not ignored in the photoionisation optical depths, which is
///    a better approx. than including it).
/// - 2011.07.12 JM: Fixed the timescales function to have a less
///    restrictive timestep limitation based only on dt=0.15/|xdot|.
/// - 2011.09.21 JM: Set the cooling to be C2 (cooling=15 in
///    cooling.cc) for RT_TEST_PROBS.
///
/// Modifications:
/// - 2011.10.06 JM: NEW FILE STARTED, FROM mp_v2_aifa.cc.  Uses
///    different integrator.
/// - 2011.10.13/14 JM: Debugging.  Fixed logic bug in evaluating
///    G0_UV/G0_IR
/// - 2011.10.17 JM: Debugging. (2011.10.22 also, and 2011.11.17).
/// - 2012.01.26 JM: Minor mods.  Added "clean" version of ydot()
///    because it was unreadable with all the ifdefs and if
///    statements.
/// - 2012.02.27 JM: Debugging comments added/removed.
/// - 2012.04.19 JM: Added "PUREHYDROGEN" ifdef for the Iliev et al.
///    2009 tests.
/// - 2012.06.22 JM: Added Wolfire+,2003,ApJ,587,278) corrections to
///    Henney+ scheme.
/// - 2012.06.25 JM: Added more heating/cooling from Wolfire+(2003).
/// - 2012.10.03 JM: Added METALLICITY multiplier to all heating and
///    cooling rates that depend on dust or metals, and also to the
///    dust opacity of the ISM to FUV radiation.  It's a crude hack,
///    because low-Z galaxies are usually overabundant in CNO, but
///    it's a start.
/// - 2013.02.14 JM: Started modifying this to include He/Metal mass
///    fractions as EP parameters, to make metallicity and mu into
///    variables that can be set from the parameterfile.
/// - 2013.02.15 JM: Moved much ifdef stuff into new classes.
///    Changed Oxygen abundance from Lodders2003 to Asplund+2009.
/// - 2013.03.21 JM: Removed redundant ifdeffed stuff.
/// - 2013.08.12 JM: added get_recombination_rate() public function.
/// - 2013.09.03 JM: fixed minor bug in He-Bremsstrahlung function.
///    Added HE_INERT ifdef for case where He is always neutral.
/// - 2013.09.05 JM: Added critical density scaling for CII cooling
///    from electron impact excitation (Wolfire+,03,eq.C2).
///    Changed forbidden-line cooling to fit to Raga,Mellema, &
///    Lundqvist (1997) tables for singly ionised C,N,O.
/// - 2014.03.27 JM: fixed bug in discrete monochromatic PI rate.
/// - 2014.09.22 JM: Added  total_cooling_rate() function to get the
///    cooling rates per cell for postprocessing.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.07 JM: New trtype array structure in constructor.
/// - 2015.07.16 JM: added pion_flt datatype (double or float).
/// - 2015.08.05 JM: more pion_flt datatype changes.
/// - 2016.06.21 JM: Temperature() threadsafe.
/// - 2018.03.20 JM: Renamed file.
///
/// NOTE: Oxygen abundance is set to 5.81e-4 from Lodders (2003,ApJ,
///       591,1220) which is the 'proto-solar nebula' value. The
///       photospheric value is lower 4.9e-4, and that is used by
///       Wiersma et al. (2009,MN,393,99).
/// UPDATE: changed to 5.37e-4 to match Asplund+(2009)
/// UPDATE: changed to 5.37e-4*0.77 = 4.1e-4 because Lodders (2003)
///       says that 23 per cent of O is in solid phase.
///

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ================================================================
// ================================================================

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

// ##################################################################
// ##################################################################

#include "constants.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG

#include "microphysics/MPv3.h"

using namespace std;

//#define HIGHDENS_CUTOFF ///< decreases CIE cooling exponentially with
// exp(-(nH/1000)^2)
#define HE_INERT
//#define MPV3_DEBUG

//#define DUSTCOOL

//
// The timestep-limiting is set by ifdef in
// source/defines/functionality_flags.h
// DT02 is recommended for MPv3 (see Mackey,2012,A&A,539,A147)
//
#if MPV3_DTLIMIT == 0
#define DTFRAC 1.0
#elif MPV3_DTLIMIT == 1
#define DTFRAC 0.5
#elif MPV3_DTLIMIT == 2
#define DTFRAC 0.25
#elif MPV3_DTLIMIT == 3
#define DTFRAC 0.125
#elif MPV3_DTLIMIT == 4
#define DTFRAC 0.0625
#elif MPV3_DTLIMIT == 5
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.5
#elif MPV3_DTLIMIT == 6
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.25
#elif MPV3_DTLIMIT == 7
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.125
#elif MPV3_DTLIMIT == 8
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.0625
#elif MPV3_DTLIMIT == 9
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.5
#elif MPV3_DTLIMIT == 10
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.25
#elif MPV3_DTLIMIT == 11
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.125
#elif MPV3_DTLIMIT == 12
#define USE_RELATIVE_NEUFRAC_DTLIMIT
#define ENERGY_CHANGE_TIMESTEP_LIMIT
#define DTFRAC 0.0625
#else
#error "Must define a DTXX timestep limit in mp_explicit"
#endif

// ##################################################################
// ##################################################################

void MPv3::get_error_tolerances(
    double *reltol,  ///< relative error tolerance.
    double atol[]    ///< absolute error tolerances
)
{
  *reltol       = JM_RELTOL;
  atol[lv_H0]   = JM_MINNEU;  ///< minimum neutral fraction I care about.
  atol[lv_eint] = JM_MINERG;  ///< E_{int}: for n=1.0, T=1.0e4, ==>
                              ///< E=2.07e-12, so say 1e-17?
  return;
}

// ##################################################################
// ##################################################################

void MPv3::get_problem_size(
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

MPv3::MPv3(
    const int nd,    ///< grid dimensions
    const int csys,  ///< Coordinate System flag
    const int nv,    ///< Total number of variables in state vector
    const int ntr,   ///< Number of tracer variables in state vector.
    const std::string *tracers,   ///< List of what the tracer variables mean.
    struct which_physics *ephys,  ///< pointer to extra physics flags.
    struct rad_sources *rsrcs,    ///< radiation sources.
    const double g                ///< EOS Gamma
    ) :
    microphysics_base(nv, ntr, tracers, ephys, rsrcs),
    ndim(nd), eos_gamma(g), coord_sys(csys)
{
  spdlog::debug("MPv3: a microphysics class");

  // ----------------------------------------------------------------
  // ------- Set up tracer variables (find which one is H+). --------
  // ----------------------------------------------------------------
  spdlog::debug(
      "\t\tSetting up Tracer Variables.  Assuming tracers are last {} variables in state vec",
      ntracer);
  int len = ntracer;

  //
  // Find H+ fraction in tracer variable list.
  //
  int ftr = nv_prim - ntracer;  // first tracer variable.
  string s;
  pv_Hp   = -1;
  pv_WIND = -1;

  for (int i = 0; i < len; i++) {
    s = tracers[i];  // Get 'i'th tracer variable.
    if (s == "H1+___" || s == "HII__" || s == "H1+" || s == "HII") {
      pv_Hp = ftr + i;
      spdlog::debug(
          "\t\tGot H+ as the {}th element of P[] (zero offset)", pv_Hp);
    }
    if (s == "WIND") {
      pv_WIND = ftr + i;
      spdlog::debug(
          "\t\tGot WIND tracer as the {}th element of P[] (zero offset)",
          pv_WIND);
    }
  }
  if (pv_Hp < 0)
    spdlog::error(
        "{}: {}", "No H ionisation fraction found in tracer list", tracers[0]);

  // ================================================================
  // ================================================================
  spdlog::debug("MPv3:: EP and RS: {}\t{}", fmt::ptr(EP), fmt::ptr(RS));

  // ----------------------------------------------------------------
  // --- Set up local variables: ion fraction and internal energy density.
  // ----------------------------------------------------------------
  k_B = pconst.kB();   // Boltzmann constant.
  m_p = pconst.m_p();  // Proton mass.

  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  // Assume metal content is low enough to ignore it.
  //
  double X        = 1.0 - EP->Helium_MassFrac;
  mean_mass_per_H = m_p / X;
  //
  // Assume He is singly ionised whenever H is, and no metals exist
  // as far as number density of electrons/ions are concerned.
  //
  JM_NION  = 1.0 + 0.25 * EP->Helium_MassFrac / X;
  JM_NELEC = 1.0 + 0.25 * EP->Helium_MassFrac / X;
#ifdef HE_INERT
  // JM_NION is the number of neutral atoms per H atom
  JM_NION  = 1.0 + 0.25 * EP->Helium_MassFrac / X;
  JM_NELEC = 1.0;                             // if He is always neutral.
#endif                                        // HE_INERT
  METALLICITY = EP->Metal_MassFrac / 0.0142;  // in units of solar.
  spdlog::debug("Metallicity = {} of solar (0.0142)", METALLICITY);

  setup_local_vectors();
  gamma_minus_one = eos_gamma - 1.0;
  Min_NeutralFrac = JM_MINNEU;
  Max_NeutralFrac = 1.0 - JM_MINNEU;

  //
  // initialise all the radiation variables to values that limit their heating
  // and cooling abilities.
  //
  mpv_nH       = 1.0;  // Hydrogen number density (density of H+ and H)
  mpv_Vshell   = 1.0e54;
  mpv_Tau0     = 1.0e6;
  mpv_dTau0    = 1.0;
  mpv_G0_UV    = 0.0;
  mpv_G0_IR    = 0.0;
  mpv_delta_S  = 0.0;
  mpv_NIdot    = 0.0;
  ion_src_type = -1;
  N_ion_srcs   = 0;
  N_diff_srcs  = 0;
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // Set flags for whether we have diffuse and ionising radiation sources.
  // ----------------------------------------------------------------
  for (int isrc = 0; isrc < RS->Nsources; isrc++) {
    // diffuse source at infinity, or point source for UV heating.
    if ((RS->sources[isrc].type == RT_SRC_DIFFUSE)
        || (RS->sources[isrc].type == RT_SRC_SINGLE
            && RS->sources[isrc].effect == RT_EFFECT_UV_HEATING))
      N_diff_srcs++;
    // point source with monochromatic or multi-frequency ionising photons.
    if (RS->sources[isrc].type == RT_SRC_SINGLE
        && (RS->sources[isrc].effect == RT_EFFECT_PION_MONO
            || RS->sources[isrc].effect == RT_EFFECT_MFION)) {
      N_ion_srcs++;
      //
      // Set the source type and Luminosity for the ionising source:
      //
      mpv_NIdot = RS->sources[isrc].strength;
      if (RS->sources[isrc].effect == RT_EFFECT_PION_MONO)
        ion_src_type = RT_EFFECT_PION_MONO;
      else
        ion_src_type = RT_EFFECT_MFION;
    }
  }
  spdlog::debug(
      "\t\tMPv3: got {} diffuse and {} ionising sources", N_diff_srcs,
      N_ion_srcs);
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // ----------------------- DIFFUSE RADIATION ----------------------
  // Set up solid angles for diffuse radiation, if needed.
  // ----------------------------------------------------------------
  if (N_diff_srcs > 0) setup_diffuse_RT_angle();
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // ----------------------------------------------------------------
  // ------------------------- IONISING SOURCE ----------------------
  if (N_ion_srcs) {
    if (N_ion_srcs > 1)
      spdlog::error("{}: {}", "too many ionising source in MPv3()", N_ion_srcs);
    //
    // Need to set up the multifrequency tables if needed.
    // Monochromatic sources don't need any setup.
    //
    for (int isrc = 0; isrc < RS->Nsources; isrc++) {
      if ((RS->sources[isrc].type == RT_SRC_SINGLE)
          && (RS->sources[isrc].effect == RT_EFFECT_MFION)
          && (RS->sources[isrc].EvoFile == "NOFILE")) {
        int err =
            set_multifreq_source_properties(&RS->sources[isrc], &mpv_NIdot);
        if (err)
          spdlog::error(
              "{}: {}", "multifreq photoionisation setup failed in MPv3 const.",
              err);
      }
    }
  }
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // We want to set up the CIE cooling function for metals-only from WSS09
  // (i.e. with H+He cooling subtracted out).
  // ----------------------------------------------------------------
  setup_WSS09_CIE_OnlyMetals();
  // ================================================================
  // ================================================================
  gen_mpv3_lookup_tables();

  // ----------------------------------------------------------------
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  int err = setup_cvode_solver();
  if (err) {
    spdlog::debug("Failed to setup cvode solve: {}", err);
    exit(1);
  }
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // ---------- output cooling rates for various temperatures -------
  // ----------------------------------------------------------------
#ifdef MPV3_DEBUG
  double p[nv_prim], P[nvl], temp1 = 0.0, ne = 0.0;
  vector<rt_source_data> rt;
  N_ion_srcs = 0;
  p[RO]      = 2.338e-24;
  p[PG]      = 1.0e-12;
  p[pv_Hp]   = 0.99;
  if (pv_WIND >= 0) p[pv_WIND] = 0.1;
  mpv_nH = 1.0e0;

  string opfile("cooling_MPv3.txt");
  ofstream outf(opfile.c_str());
  if (!outf.is_open()) spdlog::error("{}: {}", "couldn't open outfile", 1);
  outf << "Cooling Curve Data: Temperature(K) Rates(erg/cm^3/s) x=0.99999,";
  outf << " x=0.00001, x=0.5 (n=1 per cc)\n";
  outf.setf(ios_base::scientific);
  outf.precision(6);
  double T = EP->MinTemperature, Edi = 0.0, Edn = 0.0, Edpi = 0.0, junk = 0.0;
  setup_radiation_source_parameters(p, P, 0, rt, 0, rt);
  do {
    // Get T-vector index
    int iT     = 0;
    size_t ihi = lt.NT - 1, ilo = 0, imid = 0;
    do {
      imid = ilo + floor((ihi - ilo) / 2.0);
      if (lt.T[imid] < T)
        ilo = imid;
      else
        ihi = imid;
    } while (ihi - ilo > 1);
    iT        = ilo;
    double dT = T - lt.T[iT];

    // Get ne-vector index
    int ie = 0;
    ne     = 0.5;
    ihi = lt.NT - 1, ilo = 0, imid = 0;
    do {
      imid = ilo + floor((ihi - ilo) / 2.0);
      if (lt.ne[imid] < ne)
        ilo = imid;
      else
        ihi = imid;
    } while (ihi - ilo > 1);
    ie         = ilo;
    double dne = ne - lt.ne[ie];

    p[pv_Hp] = 0.99999;
    Set_Temp(p, T, junk);
    NV_Ith_S(y_in, lv_H0)   = 1.0 - p[pv_Hp];
    NV_Ith_S(y_in, lv_eint) = p[PG] / gamma_minus_one;
    ydot(0, y_in, y_out, 0);
    Edi = NV_Ith_S(y_out, lv_eint);

    p[pv_Hp] = 0.05;
    Set_Temp(p, T, junk);
    NV_Ith_S(y_in, lv_H0)   = 1.0 - p[pv_Hp];
    NV_Ith_S(y_in, lv_eint) = p[PG] / gamma_minus_one;
    ydot(0, y_in, y_out, 0);
    Edpi = NV_Ith_S(y_out, lv_eint);

    p[pv_Hp] = 0.00001;
    Set_Temp(p, T, junk);
    NV_Ith_S(y_in, lv_H0)   = 1.0 - p[pv_Hp];
    NV_Ith_S(y_in, lv_eint) = p[PG] / gamma_minus_one;
    ydot(0, y_in, y_out, 0);
    Edn = NV_Ith_S(y_out, lv_eint);

    outf << T << "\t" << Edi / mpv_nH / mpv_nH;
    outf << "  " << Edn / mpv_nH / mpv_nH;
    outf << "  " << Edpi / mpv_nH / mpv_nH;

    outf << "  " << lt.cirh[iT] + dT * lt.s_cirh[iT];
    outf << "  " << lt.C_cih0[iT] + dT * lt.s_C_cih0[iT];
    outf << "  " << (lt.rrhp[iT] + dT * lt.s_rrhp[iT]);
    outf << "  " << (lt.C_rrh[iT] + dT * lt.s_C_rrh[iT]);
    outf << "  " << (lt.C_ffhe[iT] + dT * lt.s_C_ffhe[iT]);
    outf << "  " << (lt.C_cxh0[iT] + dT * lt.s_C_cxh0[iT]);
    outf << "  "
         << lt.H_pah[iT][ie] + dT * lt.st_H_pah[iT][ie]
                + dne * lt.se_H_pah[iT][ie];
    outf << "  " << (lt.C_fbdn[iT] + dT * lt.s_C_fbdn[iT]);
    outf << "  " << (lt.C_cie[iT] + dT * lt.s_C_cie[iT]);
    outf << "  " << (lt.C_cxch[iT] + dT * lt.s_C_cxch[iT]);
    outf << "  " << (lt.C_cxo[iT] + dT * lt.s_C_cxo[iT]);
    outf << "  " << (lt.C_cxce[iT] + dT * lt.s_C_cxce[iT]);
    outf << "  "
         << lt.C_pah[iT][ie] + dT * lt.st_C_pah[iT][ie]
                + dne * lt.se_C_pah[iT][ie];
    outf << "  " << (lt.C_dust[iT] + dT * lt.s_C_dust[iT]);
    // outf <<"  "<< ;
    outf << "\n";
    T *= 1.02;
  } while (T < min(1.0e5, EP->MaxTemperature));
  outf.close();
#endif  // Debug
        // ================================================================
        // ================================================================

  spdlog::debug("MPv3: Constructor finished and returning.");
}

// ##################################################################
// ##################################################################

void MPv3::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  nvl         = 2;  // two local variables to integrate
  N_extradata = 0;
  N_equations = 2;
#if defined(CVODE6)
  y_in  = N_VNew_Serial(N_equations, sunctx);
  y_out = N_VNew_Serial(N_equations, sunctx);
#else
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
#endif
  lv_H0   = 0;  // x(H0) is the first element in the array
  lv_eint = 1;  // E_{int} is the second element.
  // cout<<"!!!!!!!!!!!!!!!!!! nvl="<<nvl<<"\n";
  return;
}

// ##################################################################
// ##################################################################

void MPv3::setup_diffuse_RT_angle()
{
  //
  // Set up solid angles for diffuse radiation.
  //
  diff_angle.clear();
  if (coord_sys == COORD_CRT && ndim == 3) {
    diff_angle.resize(6);
    for (int v = 0; v < 6; v++)
      diff_angle[v] = 4.0 * M_PI / 6.0;
  }
  else if (coord_sys == COORD_CRT && ndim == 2) {
    diff_angle.resize(4);
    for (int v = 0; v < 4; v++)
      diff_angle[v] = 2.0 * M_PI / 4.0;
  }
  else if (coord_sys == COORD_CRT && ndim == 1) {
    diff_angle.resize(2);
    for (int v = 0; v < 2; v++)
      diff_angle[v] = 1.0;
  }
  else if (coord_sys == COORD_CYL && ndim == 2) {
    diff_angle.resize(3);
    //
    // for each source in turn, we get its direction and set the angle
    // accordingly.
    //
    int count = 0;
    int dir   = -1;
    for (int isrc = 0; isrc < RS->Nsources; isrc++) {
      if (RS->sources[isrc].type == RT_SRC_DIFFUSE) {
        for (int v = 0; v < ndim; v++) {
          if (RS->sources[isrc].pos[v] > 1.0e99) dir = 2 * v + 1;
          if (RS->sources[isrc].pos[v] < -1.0e99) dir = 2 * v;
        }
        if (dir < 0)
          spdlog::error("{}: {}", "Diffuse source not at infinity!", isrc);
        //
        // if direction is in Z then angle is as for 3D, and if R+ then
        // 4x3D values.
        //
        if (dir == ZNcyl || dir == ZPcyl)
          diff_angle[count] = 4.0 * M_PI / 6.0;
        else if (dir == RPcyl)
          diff_angle[count] = 16.0 * M_PI / 6.0;
        else
          spdlog::error("{}: {}", "Bad source direction", dir);
        count++;
      }
    }
    spdlog::debug(
        "Angles for diffuse sources: [{}, {}, {}]", diff_angle[0],
        diff_angle[1], diff_angle[2]);
  }
  else if (coord_sys == COORD_SPH && ndim == 1) {
    //
    // for spherical symmetry the only diffuse source is at infinity
    //
    diff_angle.resize(1);
    diff_angle[0] = 4.0 * M_PI;
  }
  else {
    spdlog::error("{}: {}", "Unhandled coord-sys/ndim combination", ndim);
  }
  return;
}

// ##################################################################
// ##################################################################

MPv3::~MPv3()
{
  diff_angle.clear();
  //
  // Free vector memory
  //
  N_VDestroy_Serial(y_in);
  N_VDestroy_Serial(y_out);
}

// ##################################################################
// ##################################################################

int MPv3::Tr(const string s)
{
  if (s == "H1+___" || s == "HII__" || s == "H1+" || s == "HII") {
    return pv_Hp;
  }
  if (s == "WIND") {
    return pv_WIND;
  }
  else {
    return -1;
  }
}

// ##################################################################
// ##################################################################

int MPv3::set_multifreq_source_properties(
    const struct rad_src_info *rsi,  ///< source data
    double *str  ///< O/P source strength in different energy bins.
)
{
  // cout <<"MPv3: updating radiation source properties\n";
  //
  // Some sanity checks:
  // - make sure source is multi-freq and ionising
  // - make sure Rstar and Tstar are positive and finite
  //
  if (rsi->effect != RT_EFFECT_MFION)
    spdlog::error("{}: {}", "Source is not multi-frequency!", rsi->id);
  if (rsi->Rstar < 0 || !isfinite(rsi->Rstar))
    spdlog::error("{}: {}", "Source has bad Rstar parameter", rsi->Rstar);
  if (rsi->Tstar < 0 || !isfinite(rsi->Tstar))
    spdlog::error("{}: {}", "Source has bad Tstar parameter", rsi->Tstar);

  //
  // mincol is the minimum Tau we care about, maxcol is the max Tau
  // we care about, Emax is the max. energy (in ergs) to integrate
  // to in the fit.  I think 100eV is OK
  //
  double mincol = 1.0e-4, maxcol = 1.0e6, Emax = 1000.0 * 1.602e-12;
  int Nspl = 75, Nsub = 800;

#ifdef HE_INERT
  //
  // If He is always neutral, then we can't have any He-ionising
  // photons, so just integrate up to its first ionisation energy.
  //
  Emax = 24.59 * 1.602e-12;
#else
  Emax  = 54.41778 * 1.602e-12;  // assume nothing doubly-ionized He.
  // Emax = 100.0*1.602e-12;  // This is better for cosmology RT tests.
#endif  // HE_INERT

  //
  // JMackey (2019):
  // For ionisation rate, we need to rescale the Blackbody luminosity so
  // that it is much smaller for T<30000K, b/c the actual ionising photon
  // luminosity of these stars is much less than indicated by BB curve.
  // I took data from Table 1 of Diaz-Miller, Franco, & Shore,
  // (1998,ApJ,501,192), compared them to the ionising photon luminosity
  // of a BB with the same radius and Teff, and got the following scaling
  // factor, using file conversion.py in code_misc/testing/planck_fn/
  // The radius is fitted vs. ionizing photon luminosity, the scaling
  // is then applied to the luminosity of the star.
  //
  double Lcorr = 1.0, Rcorr = 1.0;
  if (rsi->Tstar < 33979.25687) {
    // Luminosity is in erg/s (rsi->strength); the fit is logspace
    // in solar units.  Converted to linear units, the relation is:
    //   L(corrected) = L * 10**(-2*21.09342323) * T**(2.0*4.65513741)
    // or
    //   L(corrected) = 6.50359577123e-43 * L * T**9.31027482
    // for T<33979.25687 K.
    //
    Lcorr = 6.50359577123e-43 * pow(rsi->Tstar, 9.31027482);
    Rcorr = sqrt(Lcorr);
    // Rcorr = sqrt(rsi->strength*Lcorr/
    //  (4.0*pconst.pi()*pconst.StefanBoltzmannConst()*pow(rsi->Tstar, 4.0))
    //  ) /pconst.Rsun();
    // cout <<"Rstar="<<Rcorr<<", rsi->Rstar="<<rsi->Rstar<<",
    // sqrt(Lcorr)="<<sqrt(Lcorr)<<"\t"; Rcorr /= rsi->Rstar; cout
    // <<"Lcorr="<<Lcorr<<", Rcorr="<<Rcorr<<"\n";
    //
    // Now need to multiply rsi->strength and rsi->Rstar by the two
    // correction factors when setting up the tables.
  }
  // cout <<"MPv3: Updating Source: T="<<rsi->Tstar<<", L=";
  // cout <<rsi->strength*Lcorr<<", R=";
  // cout <<rsi->Rstar*pconst.Rsun()*Rcorr<<"\n";

  //
  // Call the function in hydrogen_photoion.
  //
  Setup_photoionisation_rate_table(
      rsi->Tstar, rsi->Rstar * pconst.Rsun() * Rcorr, rsi->strength * Lcorr,
      mincol, maxcol, Emax, Nsub, Nspl);

  *str = rsi->strength;  // this doesn't do much.
  return 0;
}

// ##################################################################
// ##################################################################

double MPv3::get_temperature(
    const double nH,  ///< nH (per c.c.)
    const double E,   ///< E_int (per unit volume)
    const double xp   ///< x(H+)
)
{
  //
  // returns gas temperature according to E=nkT/(g-1) with n=nH*(1.1+1.1*x),
  // appropriate for a gas with 10% Helium by number, and if He is singly
  // ionised whenever H is ionised (n_i=1.1n_H, n_e=1.1n_H).
  //
  return gamma_minus_one * E / (k_B * get_ntot(nH, xp));
}

// ##################################################################
// ##################################################################

double MPv3::get_ntot(
    const double nH,  ///< nH
    const double xp   ///< x(H+)
)
{
  return (JM_NION + JM_NELEC * xp) * nH;
}

// ##################################################################
// ##################################################################

double MPv3::get_n_elec(const pion_flt *pv  ///< primitive state vector.
)
{
  double P[nvl];
  convert_prim2local(pv, P);
  // cout <<"nH="<<mpv_nH<<" y(H+)="<<(1.0-P[lv_H0]);
  // cout <<", scale="<<JM_NELEC<<"\n";
  return mpv_nH * (1.0 - P[lv_H0]) * JM_NELEC;
}

// ##################################################################
// ##################################################################

double MPv3::get_n_Hplus(const pion_flt *pv  ///< primitive state vector.
)
{
  double P[nvl];
  convert_prim2local(pv, P);
  return mpv_nH * (1.0 - P[lv_H0]);
}

// ##################################################################
// ##################################################################

double MPv3::get_n_ion(
    std::string ion,    ///< ion name
    const pion_flt *pv  ///< primitive state vector.
)
{
  if (ion != "H1+") {
    spdlog::debug("Bad ion! {}", ion);
    return -1.0e99;
  }
  double P[nvl];
  convert_prim2local(pv, P);
  return mpv_nH * (1.0 - P[lv_H0]) * JM_NION;
}

// ##################################################################
// ##################################################################

double MPv3::get_n_Hneutral(const pion_flt *pv  ///< primitive state vector.
)
{
  double P[nvl];
  convert_prim2local(pv, P);
  return mpv_nH * P[lv_H0];
}

// ##################################################################
// ##################################################################

int MPv3::convert_prim2local(
    const pion_flt *p_in,  ///< primitive vector (length nv_prim)
    double *p_local)
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  p_local[lv_eint] = p_in[PG] / (gamma_minus_one);
  p_local[lv_H0]   = 1.0 - p_in[pv_Hp];
#ifdef DUSTCOOL
  if (pv_WIND > 0)
    f_dust = max(0.0, min(1.0, 1 - 0 - p_in[pv_WIND]));
  else
    f_dust = 0.0;
#endif
  mpv_nH = p_in[RO] / mean_mass_per_H;

#ifdef MPV3_DEBUG
  //
  // Check for negative ion fraction, and set to a minimum value if found.
  //
  if (p_local[lv_H0] > 1.0) {
    spdlog::debug(
        "MPv3::convert_prim2local: negative ion fraction input: x(H0)={}, resetting to {}",
        p_local[lv_H0], Max_NeutralFrac);
    p_local[lv_H0] = Max_NeutralFrac;
  }
  //
  // Check for bad values:
  //
  if (p_local[lv_H0] > 1.01 || p_local[lv_H0] < -0.01) {
    spdlog::warn(
        "MPv3::convert_prim2local: bad ion fraction: x(H0)={}, resetting to [0,1]"
        << p_local[lv_H0]);
  }
#endif

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  // Check for temperature too low.
  p_local[lv_eint] =
      max(p_local[lv_eint], get_ntot(mpv_nH, p_in[pv_Hp]) * k_B
                                * EP->MinTemperature / (gamma_minus_one));

#ifdef MPV3_DEBUG
  //
  // Check for NAN/INF
  //
  for (int v = 0; v < 2; v++) {
    if (!isfinite(p_local[v]))
      spdlog::error("{}: {}", "INF/NAN input to microphysics", p_local[v]);
  }
  if (mpv_nH < 0.0 || !isfinite(mpv_nH))
    spdlog::error(
        "{}: {}", "Bad density input to MPv3::convert_prim2local", mpv_nH);
#endif  // MPV3_DEBUG

  return 0;
}

// ##################################################################
// ##################################################################

int MPv3::convert_local2prim(
    const double *p_local,
    const pion_flt *p_in,  ///< input primitive vector
    pion_flt *p_out        ///< updated primitive vector
)
{
  for (int v = 0; v < nv_prim; v++)
    p_out[v] = p_in[v];

  p_out[PG]    = p_local[lv_eint] * (gamma_minus_one);
  p_out[pv_Hp] = 1.0 - p_local[lv_H0];

#ifdef MPV3_DEBUG
  if (p_out[pv_Hp] < -10.0 * JM_RELTOL || p_out[pv_Hp] > 1.0 * (1.0 + JM_RELTOL)
      || !isfinite(p_out[pv_Hp])) {
    spdlog::debug("p_in : {}", std::vector<double>(p_in, p_in + nv_prim));
    spdlog::debug("p_out : {}", std::vector<double>(p_out, p_out + nv_prim));
    spdlog::debug(
        "p_local : {}", std::vector<double>(p_local, p_local + nv_prim));
    spdlog::error(
        "{}: {}", "Bad output H+ value in MPv3::convert_local2prim",
        p_out[pv_Hp] - 1.0);
  }
  if (p_out[PG] < 0.0 || !isfinite(p_out[PG]))
    spdlog::error(
        "{}: {}", "Bad output pressure in MPv3::convert_local2prim", p_out[PG]);
#endif  // MPV3_DEBUG

  //
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  p_out[pv_Hp] = max(
      Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[pv_Hp])));

  //
  // Set output pressure to be within required temperature range (use the
  // possibly corrected x(H+) from p_out[]).
  //
  double T = get_temperature(mpv_nH, p_local[lv_eint], p_out[pv_Hp]);
  if (T > 1.001 * EP->MaxTemperature) {
    Set_Temp(p_out, EP->MaxTemperature, 0);
    spdlog::debug(
        "MPv3::convert_local2prim() HIGH temperature encountered. T={}, obtained from nH={}, eint={}, x={}...  limiting to T={}",
        T, mpv_nH, p_local[lv_eint], p_out[pv_Hp], EP->MaxTemperature);
  }
  if (T < 0.999 * EP->MinTemperature) {
    Set_Temp(p_out, EP->MinTemperature, 0);
    spdlog::debug(
        "MPv3::convert_local2prim() LOW  temperature encountered. T={}, obtained from nH={}, eint={}, x={}...  limiting to T={}",
        T, mpv_nH, p_local[lv_eint], p_out[pv_Hp], EP->MinTemperature);
  }

  return 0;
}

// ##################################################################
// ##################################################################

double MPv3::Temperature(
    const pion_flt *pv,  ///< primitive vector
    const double         ///< eos gamma
)
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO] <= 0.0 || pv[PG] <= 0.0) {
    // cout <<"MPv3::Temperature() negative rho="<<pv[RO]<<" or
    // p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
#ifdef TEST_INF
  for (int v = 0; v < nv_prim; v++) {
    if (!isfinite(pv[v])) {
      spdlog::debug("NAN in MPv3::Temperature(): {} {}", v, pv[v]);
      spdlog::debug("prim vec : {}", pv);
      return -1.0e99;
    }
  }
#endif
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];
  convert_prim2local(pv, P);
  return get_temperature(pv[RO] / mean_mass_per_H, P[lv_eint], 1.0 - P[lv_H0]);
}

// ##################################################################
// ##################################################################

int MPv3::Set_Temp(
    pion_flt *p_pv,  ///< primitive vector.
    const double T,  ///< temperature
    const double     ///< eos gamma.
)
{
  //
  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.
  //
  if (p_pv[PG] <= 0.0) {
    // cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p_pv[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }
  double P[nvl];
  int err    = convert_prim2local(p_pv, P);
  P[lv_eint] = get_ntot(mpv_nH, p_pv[pv_Hp]) * k_B * T / (gamma_minus_one);
  err += convert_local2prim(P, p_pv, p_pv);
  return err;
}

// ##################################################################
// ##################################################################

void MPv3::get_dtau(
    const rad_source *s,   ///< pointer to radiation source struct
    const pion_flt ds,     ///< ds, thickness of the cell
    const pion_flt *p_in,  ///< input primitive vector
    pion_flt *dtau         ///< output dtau vector
)
{
  // dTau for MPv3 is n(H0)*sigma(H0)*ds at 13.6eV for the ionizing
  // source, and n(H)*sigma(dust)*ds in UV for UV heating source.
  double yh0 = 0.0;
  switch (s->s->effect) {
    case RT_EFFECT_PION_MONO:
    case RT_EFFECT_MFION:
      yh0   = max(Min_NeutralFrac, min(Max_NeutralFrac, 1.0 - p_in[pv_Hp]));
      *dtau = p_in[RO] * yh0 / mean_mass_per_H
              * Hi_monochromatic_photo_ion_xsection(JUST_IONISED) * ds;
      break;

    case RT_EFFECT_UV_HEATING:
      // cout <<"uv heating: "<< s->s->id <<"  ";
      *dtau = p_in[RO] * 5.348e-22 * METALLICITY / mean_mass_per_H * ds;
      // cout << *dtau <<"\n";
      break;

    default:
      spdlog::debug("id={}, effect={}", s->s->id, s->s->effect);
      spdlog::error(
          "{}: {}", "Bad source effect in MPv3::get_dtau()", s->s->effect);
      break;
  }
  return;
}

// ##################################################################
// ##################################################################

int MPv3::TimeUpdateMP(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    pion_flt *p_out,       ///< Destination Vector for updated values.
    const double dt,       ///< Time Step to advance by.
    const double,          ///< EOS gamma.
    const int,             ///< Switch for what type of integration to use.
    double *random_stuff   ///< Vector of extra data (column densities, etc.).
)
{
  //
  // Call the new update function, but with zero radiation sources.
  //
  std::vector<struct rt_source_data> temp;
  int err =
      TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}

// ##################################################################
// ##################################################################

int MPv3::TimeUpdateMP_RTnew(
    const pion_flt *p_in,  ///< Primitive Vector to be updated.
    const int N_heat,      ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    pion_flt *p_out,  ///< Destination Vector for updated values
                      ///< (can be same as first Vector.
    const double dt,  ///< Time Step to advance by.
    const double,     ///< EOS gamma.
    const int,        ///< Switch for what type of integration to use.
                      ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
    double *random_stuff  ///< final temperature (not strictly needed).
)
{
  // First set local variables for state vector and radiation source
  // properties.
  int err = 0;
  std::vector<double> P(nvl);
  err = convert_prim2local(p_in, &P[0]);
  if (err) {
    spdlog::error(
        "{}: {}", "Bad input state to MPv3::TimeUpdateMP_RTnew()", err);
  }
  setup_radiation_source_parameters(
      p_in, &P[0], N_heat, heat_src, N_ion, ion_src);

  for (int v = 0; v < nvl; v++)
    NV_Ith_S(y_in, v) = P[v];
  // Calculate y-dot[] to see if anything is changing significantly over dt
  double maxdelta = 0.0;
  err             = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in MPv3::TimeUpdateMP_RTnew()",
        err);
    exit(1);
  }
  for (int v = 0; v < nvl; v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out, v) * dt / NV_Ith_S(y_in, v)));
  }

  // Now if nothing is changing much, just to a forward Euler integration.
  if (maxdelta < EULER_CUTOFF) {
    for (int v = 0; v < nvl; v++) {
      NV_Ith_S(y_out, v) = NV_Ith_S(y_in, v) + dt * NV_Ith_S(y_out, v);
    }
  }
  // Otherwise do the implicit CVODE integration
  else {
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      spdlog::debug("Plocal : {}", P);
      spdlog::debug("Pprim  : {}", std::vector<double>(p_in, p_in + nv_prim));
      for (int v = 0; v < nvl; v++)
        P[v] = NV_Ith_S(y_out, v);
      spdlog::debug("Pfinal : {}", P);
      spdlog::error(
          "ds={}, vs={}, tau0={}, dtau={}", mpv_delta_S, mpv_Vshell, mpv_Tau0,
          mpv_dTau0);
      spdlog::error(
          "{}: {}", "integration failed again: MPv3::TimeUpdateMP_RTnew()",
          err);
      exit(1);
    }
  }

  // Now put the result into p_out[] and return.
  for (int v = 0; v < nvl; v++)
    P[v] = NV_Ith_S(y_out, v);
  err = convert_local2prim(&P[0], p_in, p_out);
#ifdef TEST_INF
  for (int v = 0; v < nv_prim; v++) {
    if (!isfinite(p_in[v]) || !isfinite(p_out[v])) {
      spdlog::error("NAN in MPv3 update: {}\n", v);
      spdlog::debug("Pin  : {}", std::vector<double>(p_in, p_in + nv_prim));
      spdlog::debug("Pout : {}", std::vector<double>(p_out, p_out + nv_prim));
      spdlog::debug("Ploc  : {}", P);
      // spdlog::error("{}: {}", "NAN in MPv3",P[2]);
      return 1;
    }
  }
#endif

  return err;
}



// ##################################################################
// ##################################################################



double MPv3::timescales(
    const pion_flt *p_in,  ///< Current cell state vector.
    const double,          ///< EOS gamma.
    const bool,            ///< set to 'true' if including cooling time.
    const bool,            ///< set to 'true' if including recombination time.
    const bool             ///< set to 'true' if including photo-ionsation time.
)
{
  if (RS->Nsources != 0) {
    spdlog::warn("WARNING: MPv3::timescales() using non-RT version!");
  }
  std::vector<struct rt_source_data> temp;
  double tmin = timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}

// ##################################################################
// ##################################################################

///
/// This returns the minimum timescale of all microphysical processes, including
/// reaction times for each species and the total heating/cooling time for the
/// gas. It requires the radiation field as an input, so it has substantially
/// greater capability than the other timescales function. Default setting is
/// DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double MPv3::timescales_RT(
    const pion_flt *p_in,  ///< Current cell state vector.
    const int N_heat,      ///< Number of UV heating sources.
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
    spdlog::error("{}: {}", "Bad input state to MPv3::timescales_RT()", err);
  }
  NV_Ith_S(y_in, lv_H0)   = P[lv_H0];
  NV_Ith_S(y_in, lv_eint) = P[lv_eint];

  //
  // Next set the radiation properties of the current cell.
  //
  setup_radiation_source_parameters(p_in, P, N_heat, heat_src, N_ion, ion_src);

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt() returned an error in MPv3::timescales_RT()", err);
  }

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t = HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the
  // energy change (heating timescale) is not required for accurately tracking
  // ionisation fronts (although it may be needed for cooling!). For testing
  // purposes there are ifdefs to allow the code to use a relative change in
  // neutral fraction and/or the relative change in energy as the timestep
  // criterion, rather than the default of absolute change in neutral
  // fraction.
  //
#ifdef USE_RELATIVE_NEUFRAC_DTLIMIT
  t =
      min(t, DTFRAC * max(5.0e-2, NV_Ith_S(y_in, lv_H0))
                 / (fabs(NV_Ith_S(y_out, lv_H0)) + TINYVALUE));
  // cout <<"using neutral fraction.\n";
#else
  t = min(t, DTFRAC / (fabs(NV_Ith_S(y_out, lv_H0)) + TINYVALUE));
  // if (t<1.0e8) {
  //  cout <<"limit by dx: dt=";
  //  cout <<DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\t";
  //  cout <<"... lvH0 = "<<NV_Ith_S(y_in, lv_H0)<<"... ydot =
  //  "<<NV_Ith_S(y_out, lv_H0)<<" \n"; rep.printVec("p_in",p_in,nv_prim);
  //}
#endif

#ifdef MPV3_DEBUG
  if (t < 3.16e9) {
    spdlog::debug(
        "MP timescales: xdot={}, Edot={} t_x={}", NV_Ith_S(y_out, lv_H0),
        NV_Ith_S(y_out, lv_eint), t);
  }
#endif  // MPV3_DEBUG

#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  //
  // Now cooling/heating time to limit to X% change in energy).
  //
  t = min(
      t, DTFRAC * P[lv_eint] / (fabs(NV_Ith_S(y_out, lv_eint)) + TINYVALUE));
#endif

#ifdef MPV3_DEBUG
  if (t < 3.16e9) {
    spdlog::debug(" and min(t_x,t_e)={}", t);
    spdlog::debug("P[1-x,E] : {}", P);
  }
#endif  // MPV3_DEBUG
  return t;
}

// ##################################################################
// ##################################################################

double MPv3::total_cooling_rate(
    const pion_flt *p_in,  ///< primitive input state vector.
    const int N_heat,      ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src,
    ///< list of ionising src column densities and source properties.
    const double  ///< EOS gamma.
)
{
  double P[nvl];
  int err = convert_prim2local(p_in, P);
  if (err) {
    spdlog::error(
        "{}: {}", "Bad input state to MPv3::total_cooling_rate()", err);
  }
  NV_Ith_S(y_in, lv_H0)   = P[lv_H0];
  NV_Ith_S(y_in, lv_eint) = P[lv_eint];

  // set the radiation properties of the current cell.
  setup_radiation_source_parameters(p_in, P, N_heat, heat_src, N_ion, ion_src);
  // Now calculate y-dot[]...
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    spdlog::error(
        "{}: {}", "dYdt(): error in MPv3::total_cooling_rate()()", err);
  }
  return NV_Ith_S(y_out, lv_eint);
}

// ##################################################################
// ##################################################################

double MPv3::get_recombination_rate(
    const int,             ///< ion index in tracer array (optional).
    const pion_flt *p_in,  ///< input state vector (primitive).
    const double g         ///< EOS gamma (optional)
)
{
  spdlog::debug("MPv3::get_recombination_rate()");
  double rate = 0.0;
  double P[nvl];
  //
  // First convert to local variables.
  //
  convert_prim2local(p_in, P);
  //
  // Now get rate
  //
  rate =
      Hii_rad_recomb_rate(get_temperature(mpv_nH, P[lv_eint], 1.0 - P[lv_H0]))
      * mpv_nH * mpv_nH * (1.0 - P[lv_H0]) * (1.0 - P[lv_H0]) * JM_NELEC;

  spdlog::debug("MPv3::get_recombination_rate()");
  return rate;
}

// ##################################################################
// ##################################################################

void MPv3::setup_radiation_source_parameters(
    const pion_flt *p_in,  ///< primitive input state vector.
    double *P,             ///< local input state vector (x_in,E_int)
    const int N_heat,      ///< Number of UV heating sources.
    const std::vector<struct rt_source_data> &heat_src,
    ///< list of UV-heating column densities and source properties.
    const int N_ion,  ///< number of ionising radiation sources.
    const std::vector<struct rt_source_data> &ion_src
    ///< list of ionising src column densities and source properties.
)
{
  //-------------------- RADIATION SOURCE INFO -----------------------
  //
  // First deal with the column densities and source strengths to get the UV
  // heating and EUV ionisation+heating rates.
  //
#ifdef MPV3_DEBUG
  if (heat_src.size() != static_cast<unsigned int>(N_heat)) {
    spdlog::error(
        "{}: {}", "Timescales: N_heating_srcs doesn't match vector size in MP3",
        heat_src.size());
  }
  if (ion_src.size() != static_cast<unsigned int>(N_ion)) {
    spdlog::error(
        "{}: {}",
        "Timescales: N_ionising_srcs doesn't match vector size in MP3",
        ion_src.size());
  }
#endif  // MPV3_DEBUG

  //
  // We need:
  // - nH  (already obtained and set in convert_prim2local())
  // - mpv_Vshell
  // - mpv_Tau0
  // - mpv_dTau0
  // - mpv_G0_UV
  // - mpv_G0_IR
  // - mpv_delta_S
  //

  //-------------------- VSHELL AND DELTA_S -----------------------
  //
  // We set mpv_Vshell to huge value so that if there is no source,
  // it will divide the effect to zero (if it gets used).
  //
  mpv_delta_S = 0.0;
  mpv_Vshell  = 1.0e200;
  if (N_ion > 0) {
    mpv_delta_S = ion_src[0].dS;
    mpv_Vshell  = ion_src[0].Vshell;
  }
  if (N_ion == 0 && N_heat > 0) {
    //
    // if no ionising sources, there may still be a point UV-heating source
    // which we need mpv_Vshell for.  here DelCol = rho*ds.
    //
    for (int v = 0; v < N_heat; v++) {
      if (heat_src[v].type == RT_SRC_SINGLE) {
#ifdef MPV3_DEBUG
        if (P[lv_H0] > 0.01) {
          spdlog::debug(
              "setup_rad_src_params: heating:  ds={}, mpv_Vshell={}",
              heat_src[v].DelCol[0] / p_in[RO], heat_src[v].Vshell);
        }
#endif  // MPV3_DEBUG
        mpv_delta_S = heat_src[v].dS;
        mpv_Vshell  = heat_src[v].Vshell;
      }
    }
  }
  //-------------------- VSHELL AND DELTA_S -----------------------

  //-------------------- mpv_Tau0 and mpv_dTau0 -----------------------
  //
  // If the point source is ionising, we need N(H0), dN(H0).  The
  // luminosity is already set in the MP class constructor by a call to
  // Setup_photoionisation_rate_table(Tstar,Rstar,etc).
  // Most H-ionising photons can't ionise He, so best that it doesn't
  // contribute to the opacity.
  // Column and DelCol are Tau and delta-Tau, to the cell and through
  // the cell, respectively.
  //
  //
  if (N_ion > 0) {
    mpv_Tau0  = ion_src[0].Column[0];
    mpv_dTau0 = ion_src[0].DelCol[0];
  }
  else {
    mpv_Tau0  = 0.0;
    mpv_dTau0 = 0.0;
  }
  //-------------------- mpv_Tau0 and mpv_dTau0 -----------------------

  //-------------------- mpv_G0_UV and mpv_G0_IR -----------------------
  //
  // Hard-coded to assume UV radiation opacity comes from dust with
  // sigma=5e-22 cm2, and that the input is \int \rho ds, integrated
  // to the front edge of the cell. Diffuse radiation strength is an
  // intensity (i.e. per solid angle), so the total flux going into
  // the solver will be roughly: sum(I(Omega)*delta-Omega)exp(-1.9Av)
  // This comes from Henney et al. (2009) eq. A3.
  //
  // Column = Av to cell from source, DelCol is delta-Av through cell
  // TODO: Check value of 1.9 for extinction.
  //
  mpv_G0_UV = 0.0;
  mpv_G0_IR = 0.0;
  if (N_heat > 0) {
    double temp  = 0.0;
    int i_diff   = 0;
    double Av_UV = 1.90;
    double Av_IR = 0.05;

    for (int v = 0; v < N_heat; v++) {
      if (heat_src[v].type == RT_SRC_DIFFUSE) {
        temp = heat_src[v].strength[0] * diff_angle[i_diff];
        spdlog::debug(
            "setup_rad_src_params:\tdiffuse src: id={} 1.9Av={}, strength={}, angle={}: attenuated flux={}",
            heat_src[v].id, Av_UV * heat_src[v].Column[0],
            heat_src[v].strength[0], diff_angle[i_diff],
            temp * exp(-Av_UV * heat_src[v].Column[0]));
        mpv_G0_UV += temp * exp(-Av_UV * heat_src[v].Column[0]);
        mpv_G0_IR += temp * exp(-Av_IR * heat_src[v].Column[0]);
        i_diff++;
        spdlog::debug(
            "UV_diff_flux={} Col={} Av={}",
            temp * exp(-Av_UV * heat_src[v].Column[0]), heat_src[v].Column[0],
            Av_UV * heat_src[v].Column[0] / 1.9);
      }
      else {
        //
        // This source must be a point source of UV heating.
        // In this case the strength is
        // the photon luminosity, so flux = L*ds*exp(-1.9Av)/mpv_Vshell
        //
        // cout <<"heat_src[v].strength "<<heat_src[v].strength[0];
        temp = heat_src[v].strength[0] * mpv_delta_S / heat_src[v].Vshell;
        spdlog::debug(
            "setup_rad_src_params:\tpoint   src: id={} 1.9Av={}, strength={}, ds={}, mpv_Vshell={}: attenuated flux={}",
            heat_src[v].id, Av_UV * heat_src[v].Column[0],
            heat_src[v].strength[0], mpv_delta_S, heat_src[v].Vshell,
            temp * exp(-Av_UV * heat_src[v].Column[0]));
        mpv_G0_UV += temp * exp(-Av_UV * heat_src[v].Column[0]);
        mpv_G0_IR += temp * exp(-Av_IR * heat_src[v].Column[0]);
        spdlog::debug(
            "UV_ptsc_flux={}\nUV_ptsc_flux={} Col={} Av={}",
            temp * exp(-Av_UV * heat_src[v].Column[0]),
            temp * exp(-Av_UV * heat_src[v].Column[0]), heat_src[v].Column[0],
            Av_UV * heat_src[v].Column[0] / 1.9);
      }
    }  // loop over heating sources.
    //
    // now divide by 1.2e7 to get it normalised correctly for Will's
    // equation A3.
    //
    mpv_G0_UV /= 1.2e7;
    mpv_G0_IR /= 1.2e7;
    // cout <<"  "<<mpv_G0_UV<<"  "<<mpv_G0_IR<<"\n";
#ifdef MPV3_DEBUG
    if (mpv_G0_UV > 1.0)
      spdlog::debug(
          "\tTotal UV attenuated flux = {} in units of 1.2e7 phot/cm2/s\n",
          mpv_G0_UV);
#endif  // MPV3_DEBUG
  }     // If there are UV heating sources
  else {
    mpv_G0_UV = 0.0;
    mpv_G0_IR = 0.0;
  }
  //-------------------- mpv_G0_UV and mpv_G0_IR -----------------------

  spdlog::debug(
      "MPv3: ionising: ds={}, mpv_Vshell={}: mpv_Tau0={}, mpv_dTau0={}, nH={}",
      mpv_delta_S, mpv_Vshell, mpv_Tau0, mpv_dTau0, mpv_nH);
  //}
}



// ##################################################################
// ##################################################################



int MPv3::ydot(
    double,                ///< current time (UNUSED)
    const N_Vector y_now,  ///< current Y-value
    N_Vector y_dot,        ///< vector for Y-dot values
    const double *         ///< extra user-data vector (UNUSED)
)
{
  //
  // fixes min-neutral-fraction to Min_NeutralFrac
  //
  double OneMinusX = max(NV_Ith_S(y_now, lv_H0), Min_NeutralFrac);
  double E_in      = NV_Ith_S(y_now, lv_eint);
  double x_in      = 1.0 - OneMinusX;

  // First get the temperature.
  double T = get_temperature(mpv_nH, E_in, x_in);

  // now get electron density
  double ne    = JM_NELEC * x_in * mpv_nH;
  double expnh = exp(-mpv_nH / 1.0e4);
  // We set a minimum electron density based on the idea that Carbon is singly
  // ionised in low density gas.  y(C)=1.5e-4 in the gas phase (by number)
  // (Sofia,1997), with an exponential cutoff at high densities.
  ne += mpv_nH * 1.5e-4 * METALLICITY * expnh;

  //  int iT = -1;
  //  do {iT++;} while (T >= lt.T[iT]);
  //  if (iT>0) iT--;
  //  double dT = T - lt.T[iT];
  //  //cout <<"T="<<T<<", T[it]="<<lt.T[iT]<<", T[it+1]="<<lt.T[iT+1]<<"\n";
  //
  //  int ie = -1;
  //  do {ie++;} while (ne >= lt.ne[ie]);
  //  if (ie>0) ie--;
  //  double dne = ne - lt.ne[ie];

  // Get T-vector index
  int iT     = 0;
  size_t ihi = lt.NT - 1, ilo = 0, imid = 0;
  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (lt.T[imid] < T)
      ilo = imid;
    else
      ihi = imid;
  } while (ihi - ilo > 1);
  iT        = ilo;
  double dT = T - lt.T[iT];
  // cout <<"T="<<T<<", T[it]="<<lt.T[iT]<<", T[it+1]="<<lt.T[iT+1]<<"\n";

  // Get ne-vector index
  int ie = 0;
  ihi = lt.NT - 1, ilo = 0, imid = 0;
  do {
    imid = ilo + floor((ihi - ilo) / 2.0);
    if (lt.ne[imid] < ne)
      ilo = imid;
    else
      ihi = imid;
  } while (ihi - ilo > 1);
  ie         = ilo;
  double dne = ne - lt.ne[ie];
  // cout <<"ne="<<ne<<", ne[ie]="<<lt.ne[ie]<<",
  // ne[ie+1]="<<lt.ne[ie+1]<<"\n";

#ifdef TEST_INF
  if (!isfinite(mpv_nH) || !isfinite(mpv_delta_S) || !isfinite(mpv_Tau0)
      || !isfinite(mpv_Vshell) || !isfinite(E_in) || !isfinite(OneMinusX)) {
    spdlog::error(
        "NAN in ydot: {}  {}  {}  {}  {}  {}", mpv_nH, mpv_delta_S, mpv_Tau0,
        mpv_Vshell, E_in, OneMinusX);
  }
#endif

  double temp1 = 0.0, temp2 = 0.0;
  double oneminusx_dot = 0.0;  // oneminusx_dot is in units of 1/s
  double Edot          = 0.0;
  // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH
  // at the end.

  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0
  //
  temp1 = lt.cirh[iT] + dT * lt.s_cirh[iT];
  oneminusx_dot -=
      temp1 * ne * OneMinusX;  // the nH is divided out on both sides.

  temp2 = lt.C_cih0[iT] + dT * lt.s_C_cih0[iT];
  Edot -= temp2 * ne * OneMinusX;
  // cout <<"CI-CR="<< temp2*ne*OneMinusX<<"\n";
  // if (x_in>0.98 && T>2000.0) cout <<T<<"  "<< temp2*ne*OneMinusX <<"  ";

  //
  // photo-ionisation of H: photoionisation rate uses equation 18 in Mellema
  // et al. 2006 (C2-ray paper), noting that their Gamma is the rate per
  // neutral H, so we multiply by 1-x, as in their equation 11.
  //
  if (N_ion_srcs) {
    //
    // set current cell dTau0
    //
    temp1 = mpv_nH * mpv_delta_S * OneMinusX
            * Hi_monochromatic_photo_ion_xsection(JUST_IONISED);

    switch (ion_src_type) {
      case RT_EFFECT_MFION:
        //
        // Rather than divide the discretised rate by n(H0) and then
        // multiply by (1-x) to get oneminusx_dot, we simply divide by
        // n(H) since this is more numerically stable.  To do this, n(H)
        // is passed to the rate function instead of n(H0).
        //
        // Also, instead of using mpv_dTau0 for the cell optical depth,
        // we use the current optical depth (nH*mpv_delta_S*OneMinusX)
        // so that we allow the photoionisation rate to decrease as the
        // number of neutral atoms decreases during the timestep.  This
        // is more stable.
        //
        oneminusx_dot -= Hi_discrete_multifreq_photoion_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_delta_S, mpv_Vshell);
        Edot += Hi_discrete_multifreq_photoheating_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_delta_S, mpv_Vshell);
        // if (x_in>0.98 && T>2000.0) cout <<"PI-heat "<<
        // Hi_discrete_multifreq_photoheating_rate(mpv_Tau0, temp1,
        // mpv_nH, mpv_delta_S, mpv_Vshell)<<"\n";
        break;

      case RT_EFFECT_PION_MONO:
        //
        // hardcoded for a hv-13.6eV = 5.0eV monochromatic source.
        //
#define PHOTON_ENERGY 2.98e-11  // 5 eV
#define EXCESS_ENERGY 8.01e-12
        //#define PHOTON_ENERGY 2.24e-11
        //#define EXCESS_ENERGY 0.64e-12
        //#define PHOTON_ENERGY 2.499e-11  // 2 eV
        //#define EXCESS_ENERGY 3.204e-12
        temp1 = Hi_discrete_mono_photoion_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_NIdot, PHOTON_ENERGY, mpv_delta_S,
            mpv_Vshell);
        oneminusx_dot -= temp1;
        Edot += temp1 * EXCESS_ENERGY;
        // cout <<"PI-HR="<<temp1*EXCESS_ENERGY<<"\n";
        break;

      default:
        spdlog::error("{}: {}", "Bad ion_src_type in dYdt()", ion_src_type);
        break;
    }  // switch
  }

  // radiative recombination of H+
  oneminusx_dot += (lt.rrhp[iT] + dT * lt.s_rrhp[iT]) * x_in * ne;
  // Total H+ cooling: recombination plus free-free
  Edot -= (lt.C_rrh[iT] + dT * lt.s_C_rrh[iT]) * x_in * ne;
  // if (x_in>0.98 && T>2000.0) cout <<"HII-TC="<<(lt.C_rrh[iT] +
  // dT*lt.s_C_rrh[iT]) *x_in*ne <<"\n";

  // Add Helium free-free
#ifndef HE_INERT
  // Only if He is ionised, otherwise it has no free-free.
  Edot -= (lt.C_ffhe[iT] + dT * lt.s_C_ffhe[iT]) * x_in * ne;
// if (x_in>0.98 && T>2000.0) cout <<"He FF "<< (lt.C_ffhe[iT] +
// dT*lt.s_C_ffhe[iT]) *x_in*ne<<"\n";
#endif  // HE_INERT

  // collisional excitation cooling of H0 Aggarwal (1983) and
  // Raga+(1997,ApJS).
  Edot -= (lt.C_cxh0[iT] + dT * lt.s_C_cxh0[iT]) * OneMinusX * ne;
  // if (x_in>0.98 && T>2000.0) cout <<"CE-CR="<<(lt.C_cxh0[iT] +
  // dT*lt.s_C_cxh0[iT])*OneMinusX*ne<<"\n";
  //
  // --------- END OF HYDROGEN COOLING, MOVING TO METAL COOLING --------
  //
  //
  // Now we get to the less certain elements of the cooling/heating function.
  // First we do the heating:
  //
  if (N_diff_srcs) {
    //
    // UV heating due to both diffuse radiation and point source radiation.
    // The quantity mpv_G0_UV is as defined in Henney et al. (2009) Appendix
    // A1, Eq.A3, and is set in set_parameters_for_current step()
    //
    // cout <<"adding diffuse heating! ";
    Edot +=
        1.9e-26 * METALLICITY * mpv_G0_UV / (1.0 + 6.4 * (mpv_G0_UV / mpv_nH));
    // if (x_in>0.98 && T>2000.0)cout
    // <<"DfUV="<<1.9e-26*METALLICITY*mpv_G0_UV/(1.0+6.4*(mpv_G0_UV/mpv_nH));

    //
    // IR heating (HAdCM09 eq.A6) from point source and/or diffuse
    // radiation. There is a different G0 parameter because the attenuation
    // is according to exp(-0.05Av) rather than before where the coefficient
    // was 1.9.
    //
    temp1 = 1.0 + 3.0e4 / mpv_nH;
    temp1 = temp1 * temp1;
    Edot += 7.7e-32 * METALLICITY * mpv_G0_IR / temp1;
    // if (x_in>0.98 && T>2000.0)cout
    // <<"DfIR="<<7.7e-32*METALLICITY*mpv_G0_IR/temp1<<"\n";
  }

  // Cosmic ray heating (HAdCM09 eq.A7).
  Edot += 5.0e-28 * OneMinusX;
  // if (x_in>0.98 && T>2000.0)cout <<"CR-HR="<<5.0e-28*OneMinusX<<"\n";
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar neighbourhood.
  oneminusx_dot -= 1.8e-17 * OneMinusX;

  // Diffuse UV Heating rate (Wolfire+,2003,eq.20,21, Fig.10,b).
  // Using the first term of eq.19, with
  // phi_{PAH}=0.5 and G_0=1.7.  multiply by the neutral fraction OneMinusX
  // because this heating term is only calculated for warm neutral medium.
  temp1 =
      lt.H_pah[iT][ie] + dT * lt.st_H_pah[iT][ie] + dne * lt.se_H_pah[iT][ie];
  Edot += OneMinusX * temp1;
  // if (x_in>0.98 && T>2000.0) cout <<"PAH-HR="<< OneMinusX * temp1<<"\n";

  //
  // COOLING: First forbidden line cooling of e.g. OII,OIII, dominant in
  // HII regions.  This is collisionally excited lines of photoionised metals.
  // (HAdCM09 eq.A9) I have exponentially damped this at high temperatures
  // Oxygen abundance set to 5.37e-4 from
  // Asplund+(2009,ARAA), times 0.77 to account for 23% of O in solid dust.
  //
  temp1 = (lt.C_fbdn[iT] + dT * lt.s_C_fbdn[iT]) * x_in * ne;

  //
  // Now the Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve. We
  // take the actual cooling rate to be the max of SD93-CIE and the previous
  // two terms.
  //
  temp2 = (lt.C_cie[iT] + dT * lt.s_C_cie[iT]) * x_in * x_in * mpv_nH;
  // This is the CII cooling by electron collisions.
  // This rate has a very low critical density (Goldsmith, Langer et al.,
  // 2012ApJS..203...13G), at n_c=20 cm^{-3} at 1000K, so we use their
  // temperature scaling and divide by density according to
  // rate = rate/(1.0 + 0.05*nH*(T/2000K)^(-0.37))
  temp2 +=
      (lt.C_cxce[iT][ie] + dT * lt.st_C_cxce[iT][ie]
       + dne * lt.se_C_cxce[iT][ie]);
  Edot -= max(temp1, temp2);
  // if (x_in>0.98 && T>2000.0) cout <<"High-T"<< max(temp1,temp2)<<"\n";

  // Instead of the PDR cooling from Henney, use Wolfire's eq.C1,C3 for
  // collisional cooling of CII and OI by neutral H atoms.  In eq.C3 I have
  // absorbed the (100K)^{-0.4} into the prefactor, and used x^a=exp(a*ln(x)).
  // I have cut off equation C1 at high densities to be consistent with the
  // ion fraction of C that I assumed above for the electron density.
  //
  Edot -= (lt.C_cxch[iT] + dT * lt.s_C_cxch[iT]) * mpv_nH * OneMinusX * expnh;
  Edot -= (lt.C_cxo[iT] + dT * lt.s_C_cxo[iT]) * mpv_nH * OneMinusX;
  // if (x_in>0.98 && T>2000.0) cout <<"CII-HI "<< (lt.C_cxch[iT] +
  // dT*lt.s_C_cxch[iT])*mpv_nH*OneMinusX*expnh <<"\n"; if (x_in>0.98 &&
  // T>2000.0) cout <<"OI-HI  "<< (lt.C_cxo[iT]  + dT*lt.s_C_cxo[iT]
  // )*mpv_nH*OneMinusX <<"\n";

  //
  // PAH cooling: eq. 21 in Wolfire+,2003.  I think they should have
  // multiplied their equation by 1.3 for the increased PAH abundance...
  //
  Edot -=
      (lt.C_pah[iT][ie] + dT * lt.st_C_pah[iT][ie] + dne * lt.se_C_pah[iT][ie]);
  // if (x_in>0.98 && T>2000.0) cout <<"PAH-C "<<(lt.C_pah[iT][ie] +
  // dT*lt.st_C_pah[iT][ie] + dne*lt.se_C_pah[iT][ie]) <<"\n";

#ifdef DUSTCOOL
  // Dust cooling in hot gas (following Everett & Churchwell (2010)
  // figure 9, which is calculated from CLOUDY).
  if (pv_WIND > 0) {
    Edot -= ne * x_in * f_dust * (lt.C_dust[iT] + dT * lt.s_C_dust[iT]);
    // if (x_in>0.98 && T>2000.0) cout <<"DUST-C "<<
    // ne*x_in*f_dust*(lt.C_dust[iT] + dT*lt.s_C_dust[iT])<<"\n";
  }
#endif

  //
  // Multiply Edot by nH to get units of energy loss/gain per unit
  // volume per second.
  //
  Edot *= mpv_nH;
#ifdef HIGHDENS_CUTOFF
  if (Edot < 0.0) Edot *= exp(-mpv_nH * mpv_nH / 1.0e6);
#endif  // HIGHDENS_CUTOFF

  //
  // We want to limit cooling as we approach the minimum temperature, so we
  // scale the rate to linearly approach zero as we reach Tmin.
  //
  if (Edot < 0.0 && T < 2.0 * EP->MinTemperature) {
    Edot = min(0.0, (Edot) * (T - EP->MinTemperature) / EP->MinTemperature);
  }

  NV_Ith_S(y_dot, lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot, lv_eint) = Edot;

#ifdef TEST_INF
  if (!isfinite(Edot) || !isfinite(oneminusx_dot)) {
    spdlog::error(
        "NAN returns from ydot: {}  {} T={}\nydot INPUTS: {}  {}  {}  {}  {}  {}",
        oneminusx_dot, Edot, T, mpv_nH, mpv_delta_S, mpv_Tau0, mpv_Vshell, E_in,
        OneMinusX);
  }
  if (mpv_Vshell < 1.0e45) {
    spdlog::warn(
        "ydot dodgy Vshell INPUTS: {}  {}  {}  {}", mpv_nH, mpv_delta_S,
        mpv_Tau0, mpv_Vshell);
  }
#endif

  if (1 == 0 && x_in > 0.98 && T > 40.0) {
    spdlog::debug(
        "T, rates: {}  {}  {}  {}  {} \n {}  {}  {}", T, OneMinusX,
        oneminusx_dot, E_in, Edot, ne, ne * ne,
        (lt.rrhp[iT] + dT * lt.s_rrhp[iT]) * x_in * ne * 1.5 * 1.602e-12);
    temp1 = mpv_nH * mpv_delta_S * OneMinusX
            * Hi_monochromatic_photo_ion_xsection(JUST_IONISED);
    spdlog::debug(
        "\t{} {}", temp1,
        Hi_discrete_multifreq_photoheating_rate(
            mpv_Tau0, temp1, mpv_nH, mpv_delta_S, mpv_Vshell));
  }

  return 0;
}

// ##################################################################
// ##################################################################

void MPv3::gen_mpv3_lookup_tables()
{
  //  Start by generating the logarithmic temperature scale:
  lt.NT        = 200;
  double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
  lt.dlogT =
      (log10(EP->MaxTemperature) - log10(EP->MinTemperature)) / (lt.NT - 1);
  lt.dlogne = 12.0 / (lt.NT - 1);  // ne from 1e-6 to 1e6.
  lt.T.resize(lt.NT);
  lt.ne.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++) {
    lt.T[i]  = pow(10.0, log10(EP->MinTemperature) + i * lt.dlogT);
    lt.ne[i] = pow(10.0, -6.0 + i * lt.dlogne);
  }

  // collisional ionization of H
  lt.cirh.resize(lt.NT);
  lt.C_cih0.resize(lt.NT);
  lt.rrhp.resize(lt.NT);
  lt.C_rrh.resize(lt.NT);
  lt.C_ffhe.resize(lt.NT);
  lt.C_cxh0.resize(lt.NT);
  lt.C_fbdn.resize(lt.NT);

  lt.C_cie.resize(lt.NT);
  lt.C_cxch.resize(lt.NT);
  lt.C_cxo.resize(lt.NT);
  // lt.C_cxce.resize(lt.NT);
  lt.C_dust.resize(lt.NT);

  for (size_t i = 0; i < lt.NT; i++) {
    Hi_coll_ion_rates(lt.T[i], &temp1, &temp2);
    lt.cirh[i]   = temp1;
    lt.C_cih0[i] = temp2;
    lt.rrhp[i]   = Hii_rad_recomb_rate(lt.T[i]);
    lt.C_rrh[i]  = Hii_total_cooling(lt.T[i]);
    lt.C_ffhe[i] = 1.68e-27 * (JM_NION - 1.0) * sqrt(lt.T[i]);
    lt.C_cxh0[i] = Hi_coll_excitation_cooling_rate(lt.T[i])
                   * exp(-lt.T[i] * lt.T[i] / 5.0e10);
    lt.C_fbdn[i] =
        1.20e-22 * METALLICITY * 1.0
        * exp(-33610.0 / lt.T[i] - (2180.0 * 2180.0 / lt.T[i] / lt.T[i]))
        * exp(-lt.T[i] * lt.T[i] / 5.0e10);
    lt.C_cie[i]  = METALLICITY * cooling_rate_SD93CIE(lt.T[i]);
    lt.C_cxch[i] = 3.15e-27 * METALLICITY * exp(-92.0 / lt.T[i]);
    lt.C_cxo[i] =
        3.96e-28 * METALLICITY * exp(0.4 * log(lt.T[i]) - 228.0 / lt.T[i]);
    // lt.C_cxce[i]
    // = 1.4e-23*METALLICITY*exp(-0.5*log(lt.T[i])-92.0/lt.T[i]);
    lt.C_dust[i] = 1.0e-17 * exp(1.5 * log(lt.T[i] / 2.5e8));
  }

  // PAH cooling/heating and CII-e cooling depend non-linearly on ne,
  // so 2D lookup table
  lt.H_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.H_pah[i].resize(lt.NT);
  lt.C_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.C_pah[i].resize(lt.NT);
  lt.C_cxce.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.C_cxce[i].resize(lt.NT);

  for (size_t i = 0; i < lt.NT; i++) {
    for (size_t j = 0; j < lt.NT; j++) {
      lt.H_pah[i][j] = 1.083e-25 * METALLICITY
                       / (1.0 + 9.77e-3 * pow(sqrt(lt.T[i]) / lt.ne[j], 0.73));

      lt.C_pah[i][j] = 3.02e-30 * METALLICITY
                       * exp(0.94 * log(lt.T[i])
                             + 0.74 * pow(lt.T[i], -0.068)
                                   * log(3.4 * sqrt(lt.T[i]) / lt.ne[j]))
                       * lt.ne[j];

      lt.C_cxce[i][j] =
          1.4e-23 * METALLICITY * exp(-0.5 * log(lt.T[i]) - 92.0 / lt.T[i])
          * lt.ne[j] / (1.0 + 0.05 * lt.ne[j] * pow(lt.T[i] / 2000.0, -0.37));
    }
  }

  lt.s_cirh.resize(lt.NT);
  lt.s_C_cih0.resize(lt.NT);
  lt.s_rrhp.resize(lt.NT);
  lt.s_C_rrh.resize(lt.NT);
  lt.s_C_ffhe.resize(lt.NT);
  lt.s_C_cxh0.resize(lt.NT);
  lt.s_C_fbdn.resize(lt.NT);
  lt.s_C_cie.resize(lt.NT);
  lt.s_C_cxch.resize(lt.NT);
  lt.s_C_cxo.resize(lt.NT);
  // lt.s_C_cxce.resize(lt.NT);
  lt.s_C_dust.resize(lt.NT);
  for (size_t i = 0; i < lt.NT - 1; i++) {
    lt.s_cirh[i] = (lt.cirh[i + 1] - lt.cirh[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cih0[i] =
        (lt.C_cih0[i + 1] - lt.C_cih0[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_rrhp[i]  = (lt.rrhp[i + 1] - lt.rrhp[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_rrh[i] = (lt.C_rrh[i + 1] - lt.C_rrh[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_ffhe[i] =
        (lt.C_ffhe[i + 1] - lt.C_ffhe[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cxh0[i] =
        (lt.C_cxh0[i + 1] - lt.C_cxh0[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_fbdn[i] =
        (lt.C_fbdn[i + 1] - lt.C_fbdn[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cie[i] = (lt.C_cie[i + 1] - lt.C_cie[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cxch[i] =
        (lt.C_cxch[i + 1] - lt.C_cxch[i]) / (lt.T[i + 1] - lt.T[i]);
    lt.s_C_cxo[i] = (lt.C_cxo[i + 1] - lt.C_cxo[i]) / (lt.T[i + 1] - lt.T[i]);
    // lt.s_C_cxce[i] = (lt.C_cxce[i+1]-lt.C_cxce[i]) / (lt.T[i+1]-lt.T[i]);
    lt.s_C_dust[i] =
        (lt.C_dust[i + 1] - lt.C_dust[i]) / (lt.T[i + 1] - lt.T[i]);
  }
  lt.s_cirh[lt.NT - 1]   = 0.0;
  lt.s_C_cih0[lt.NT - 1] = 0.0;
  lt.s_rrhp[lt.NT - 1]   = 0.0;
  lt.s_C_rrh[lt.NT - 1]  = 0.0;
  lt.s_C_ffhe[lt.NT - 1] = 0.0;
  lt.s_C_cxh0[lt.NT - 1] = 0.0;
  lt.s_C_fbdn[lt.NT - 1] = 0.0;
  lt.s_C_cie[lt.NT - 1]  = 0.0;
  lt.s_C_cxch[lt.NT - 1] = 0.0;
  lt.s_C_cxo[lt.NT - 1]  = 0.0;
  // lt.s_C_cxce[lt.NT-1] = 0.0;
  lt.s_C_dust[lt.NT - 1] = 0.0;

  // slopes for bilinear interpolation
  lt.st_H_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.st_H_pah[i].resize(lt.NT);
  lt.se_H_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.se_H_pah[i].resize(lt.NT);
  lt.st_C_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.st_C_pah[i].resize(lt.NT);
  lt.se_C_pah.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.se_C_pah[i].resize(lt.NT);
  lt.st_C_cxce.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.st_C_cxce[i].resize(lt.NT);
  lt.se_C_cxce.resize(lt.NT);
  for (size_t i = 0; i < lt.NT; i++)
    lt.se_C_cxce[i].resize(lt.NT);
  for (size_t i = 0; i < lt.NT - 1; i++) {
    for (size_t j = 0; j < lt.NT - 1; j++) {
      lt.st_H_pah[i][j] =
          (lt.H_pah[i + 1][j] - lt.H_pah[i][j]) / (lt.T[i + 1] - lt.T[i]);
      lt.se_H_pah[i][j] =
          (lt.H_pah[i][j + 1] - lt.H_pah[i][j]) / (lt.ne[j + 1] - lt.ne[j]);
      lt.st_C_pah[i][j] =
          (lt.C_pah[i + 1][j] - lt.C_pah[i][j]) / (lt.T[i + 1] - lt.T[i]);
      lt.se_C_pah[i][j] =
          (lt.C_pah[i][j + 1] - lt.C_pah[i][j]) / (lt.ne[j + 1] - lt.ne[j]);
      lt.st_C_cxce[i][j] =
          (lt.C_cxce[i + 1][j] - lt.C_cxce[i][j]) / (lt.T[i + 1] - lt.T[i]);
      lt.se_C_cxce[i][j] =
          (lt.C_cxce[i][j + 1] - lt.C_cxce[i][j]) / (lt.ne[j + 1] - lt.ne[j]);
    }
  }
  for (size_t i = 0; i < lt.NT; i++) {
    lt.st_H_pah[i][lt.NT - 1]  = 0.0;
    lt.se_H_pah[i][lt.NT - 1]  = 0.0;
    lt.st_C_pah[i][lt.NT - 1]  = 0.0;
    lt.se_C_pah[i][lt.NT - 1]  = 0.0;
    lt.st_C_cxce[i][lt.NT - 1] = 0.0;
    lt.se_C_cxce[i][lt.NT - 1] = 0.0;
  }
  for (size_t j = 0; j < lt.NT; j++) {
    lt.st_H_pah[lt.NT - 1][j]  = 0.0;
    lt.se_H_pah[lt.NT - 1][j]  = 0.0;
    lt.st_C_pah[lt.NT - 1][j]  = 0.0;
    lt.se_C_pah[lt.NT - 1][j]  = 0.0;
    lt.st_C_cxce[lt.NT - 1][j] = 0.0;
    lt.se_C_cxce[lt.NT - 1][j] = 0.0;
  }

  return;
}

// ##################################################################
// ##################################################################
