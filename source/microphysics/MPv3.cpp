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
/// - For cooling due to heavy elements, which are not explicitly
///  included, we use the CIE cooling curve of Sutherland & Dopita
///  (1993,ApJS,88,253), but this time I subtract the zero-metals
///  curve from the solar-metallicity curve so that I only have
///  cooling due to metals.  This eliminates the potential for 
///  double-counting which was in my previous cooling function.
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


#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "microphysics/MPv3.h"

using namespace std;

//#define HIGHDENS_CUTOFF ///< decreases CIE cooling exponentially with exp(-(nH/1000)^2)
#define HE_INERT
//#define BETELGEUSE
//#define MPV3_DEBUG


//
// The timestep-limiting is set by ifdef in
// source/defines/functionality_flags.h
// DT02 is recommended for MPv3 (see Mackey,2012,A&A,539,A147)
//
#if   MPV3_DTLIMIT==0
  #define DTFRAC 1.0
#elif MPV3_DTLIMIT==1
  #define DTFRAC 0.5
#elif MPV3_DTLIMIT==2
  #define DTFRAC 0.25
#elif MPV3_DTLIMIT==3
  #define DTFRAC 0.125
#elif MPV3_DTLIMIT==4
  #define DTFRAC 0.0625
#elif MPV3_DTLIMIT==5
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.5
#elif MPV3_DTLIMIT==6
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.25
#elif MPV3_DTLIMIT==7
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.125
#elif MPV3_DTLIMIT==8
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.0625
#elif MPV3_DTLIMIT==9
  #define USE_RELATIVE_NEUFRAC_DTLIMIT
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.5
#elif MPV3_DTLIMIT==10
  #define USE_RELATIVE_NEUFRAC_DTLIMIT
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.25
#elif MPV3_DTLIMIT==11
  #define USE_RELATIVE_NEUFRAC_DTLIMIT
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.125
#elif MPV3_DTLIMIT==12
  #define USE_RELATIVE_NEUFRAC_DTLIMIT
  #define ENERGY_CHANGE_TIMESTEP_LIMIT
  #define DTFRAC 0.0625
#else
#error "Must define a DTXX timestep limit in mp_explicit"
#endif


// ##################################################################
// ##################################################################


void MPv3::get_error_tolerances(
      double *reltol, ///< relative error tolerance.
      double atol[] ///< absolute error tolerances
      )
{
  *reltol = JM_RELTOL;
  atol[lv_H0]   = JM_MINNEU; ///< minimum neutral fraction I care about.
  atol[lv_eint] = JM_MINERG; ///< E_{int}: for n=1.0, T=1.0e4, ==> E=2.07e-12, so say 1e-17?
  return;
}


// ##################################################################
// ##################################################################


void MPv3::get_problem_size(
      int *ne, ///< number of equations
      int *np  ///< number of parameters in user_data vector.
      )
{
  *ne = N_equations;
  *np = N_extradata;
  return;
}

// ##################################################################
// ##################################################################


MPv3::MPv3(
      const int nd,   ///< grid dimensions
      const int csys,   ///< Coordinate System flag
      const int nv,   ///< Total number of variables in state vector
      const int ntracer,         ///< Number of tracer variables in state vector.
      const std::string *tracers,  ///< List of what the tracer variables mean.
      struct which_physics *ephys, ///< pointer to extra physics flags.
      struct rad_sources *rsrcs,   ///< radiation sources.
      const double g  ///< EOS Gamma
      )
: microphysics_base(ephys,rsrcs),
  ndim(nd), nv_prim(nv), eos_gamma(g), coord_sys(csys)
{
  cout <<"\n---------------------------------------------------------------------\n";
  cout <<"MPv3: a microphysics class.\n";

  // ----------------------------------------------------------------
  // ------- Set up tracer variables (find which one is H+). --------
  // ----------------------------------------------------------------
  cout <<"\t\tSetting up Tracer Variables.  Assuming tracers are last ";
  cout <<ntracer<<" variables in state vec.\n";
  int len = ntracer;

  //
  // Find H+ fraction in tracer variable list.
  //
  int ftr = nv_prim -ntracer; // first tracer variable.
  string s; pv_Hp=-1;

  for (int i=0;i<len;i++) {
    s = tracers[i]; // Get 'i'th tracer variable.
    if (s=="H1+___" || s=="HII__" || s=="H1+" || s=="HII") {
      pv_Hp = ftr+i;
      cout <<"\t\tGot H+ as the "<<pv_Hp<<"th element of P[] (zero offset).\n";
    }
  }
  if (pv_Hp<0)
   rep.error("No H ionisation fraction found in tracer list",tracers[0]);

  // ================================================================
  // ================================================================
#ifdef TESTING
  cout <<"MPv3:: EP and RS: "<<EP<<"\t"<<RS<<endl;
#endif


  // ----------------------------------------------------------------
  // --- Set up local variables: ion fraction and internal energy density.
  // ----------------------------------------------------------------
  k_B = pconst.kB();  // Boltzmann constant.
  m_p = pconst.m_p(); // Proton mass.

  //
  // Get the mean mass per H atom from the He and Z mass fractions.
  // Assume metal content is low enough to ignore it.
  //
  double X = 1.0-EP->Helium_MassFrac;
  mean_mass_per_H = m_p/X;
  //
  // Assume He is singly ionised whenever H is, and no metals exist
  // as far as number density of electrons/ions are concerned.
  //
  JM_NION  = 1.0 +0.25*EP->Helium_MassFrac/X;
  JM_NELEC = 1.0 +0.25*EP->Helium_MassFrac/X;
#ifdef HE_INERT
  // JM_NION is the number of neutral atoms per H atom
  JM_NION  = 1.0 +0.25*EP->Helium_MassFrac/X;
  JM_NELEC = 1.0; // if He is always neutral.
#endif // HE_INERT
  METALLICITY = EP->Metal_MassFrac/0.0142; // in units of solar.
  cout <<"Metallicity = "<<METALLICITY<<" of solar (0.0142)\n";

  setup_local_vectors();
  gamma_minus_one = eos_gamma -1.0;
  Min_NeutralFrac     = JM_MINNEU;
  Max_NeutralFrac     = 1.0-JM_MINNEU;

  //
  // initialise all the radiation variables to values that limit their heating
  // and cooling abilities.
  //
  mpv_nH   = 1.0;  // Hydrogen number density (density of H+ and H)
  mpv_Vshell  = 1.0e54;
  mpv_Tau0     = 1.0e6;
  mpv_dTau0    = 1.0;
  mpv_G0_UV   = 0.0;
  mpv_G0_IR   = 0.0;
  mpv_delta_S = 0.0;
  mpv_NIdot   = 0.0;
  ion_src_type = -1;
  N_ion_srcs  = 0;
  N_diff_srcs = 0;
  // ================================================================
  // ================================================================
  

  // ----------------------------------------------------------------
  // Set flags for whether we have diffuse and ionising radiation sources.
  // ----------------------------------------------------------------
  for (int isrc=0; isrc<RS->Nsources; isrc++) {
    // diffuse source at infinity, or point source for UV heating.
    if ((RS->sources[isrc].type==RT_SRC_DIFFUSE) ||
        (RS->sources[isrc].type==RT_SRC_SINGLE &&
         RS->sources[isrc].effect==RT_EFFECT_UV_HEATING)
        )
      N_diff_srcs++;
    // point source with monochromatic or multi-frequency ionising photons.
    if (RS->sources[isrc].type==RT_SRC_SINGLE &&
        (RS->sources[isrc].effect==RT_EFFECT_PION_MONO ||
         RS->sources[isrc].effect==RT_EFFECT_PION_MULTI)) {
      N_ion_srcs++;
      //
      // Set the source type and Luminosity for the ionising source:
      //
      mpv_NIdot = RS->sources[isrc].strength;
      if (RS->sources[isrc].effect==RT_EFFECT_PION_MONO)
        ion_src_type=RT_EFFECT_PION_MONO;
      else
        ion_src_type=RT_EFFECT_PION_MULTI;
    }
  }
  cout <<"\t\tMPv3: got "<<N_diff_srcs<<" diffuse and ";
  cout <<N_ion_srcs<<" ionising sources.\n";
  // ================================================================
  // ================================================================


  // ----------------------------------------------------------------
  // ----------------------- DIFFUSE RADIATION ----------------------
  // Set up solid angles for diffuse radiation, if needed.
  // ----------------------------------------------------------------
  if (N_diff_srcs>0) setup_diffuse_RT_angle();
  // ================================================================
  // ================================================================


  // ----------------------------------------------------------------
  // ----------------------------------------------------------------
  // ------------------------- IONISING SOURCE ----------------------
  if (N_ion_srcs) {
    if (N_ion_srcs>1)
      rep.error("too many ionising source in MPv3()",N_ion_srcs);
    //
    // Need to set up the multifrequency tables if needed.
    // Monochromatic sources don't need any setup.
    //
    for (int isrc=0; isrc<RS->Nsources; isrc++) {
      if (  (RS->sources[isrc].type  ==RT_SRC_SINGLE) &&
            (RS->sources[isrc].effect==RT_EFFECT_PION_MULTI) ) {
          int err=set_multifreq_source_properties(&RS->sources[isrc]);
          if (err)
            rep.error("multifreq photoionisation setup failed in MPv3 const.",err);
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


  // ----------------------------------------------------------------
  // --------------------------- CVODES ----------------------------
  // Initialise the CVODES solver memory etc.
  // ----------------------------------------------------------------
  setup_cvode_solver_without_Jacobian();
  // ================================================================
  // ================================================================

  // ----------------------------------------------------------------
  // ---------- output cooling rates for various temperatures -------
  // ----------------------------------------------------------------
#ifdef MPV3_DEBUG
  double p[nv_prim];
  p[RO]=2.338e-24; p[PG]=1.0e-12;
  p[pv_Hp] = 0.99;
  mpv_nH=1.0e0;

  string opfile("cooling_MPv3.txt");
  ofstream outf(opfile.c_str());
  if(!outf.is_open()) rep.error("couldn't open outfile",1);
  outf <<"Cooling Curve Data: Temperature(K) Rates(erg/cm^3/s) x=0.99999,";
  outf <<" x=0.00001, x=0.5 (n=1 per cc)\n";
  outf.setf( ios_base::scientific );
  outf.precision(6);
  double T=EP->MinTemperature, Edi=0.0,Edn=0.0,Edpi=0.0,junk=0.0;
  do {
    p[pv_Hp] = 0.99999;
    Set_Temp(p,T,junk);
    NV_Ith_S(y_in,lv_H0)   = 1.0-p[pv_Hp];
    NV_Ith_S(y_in,lv_eint) = p[PG]/gamma_minus_one;
    ydot(0,y_in,y_out,0);
    Edi = NV_Ith_S(y_out,lv_eint);

    p[pv_Hp] = 0.05;
    Set_Temp(p,T,junk);
    NV_Ith_S(y_in,lv_H0)   = 1.0-p[pv_Hp];
    NV_Ith_S(y_in,lv_eint) = p[PG]/gamma_minus_one;
    ydot(0,y_in,y_out,0);
    Edpi = NV_Ith_S(y_out,lv_eint);

    p[pv_Hp] = 0.00001;
    Set_Temp(p,T,junk);
    NV_Ith_S(y_in,lv_H0)   = 1.0-p[pv_Hp];
    NV_Ith_S(y_in,lv_eint) = p[PG]/gamma_minus_one;
    ydot(0,y_in,y_out,0);
    Edn = NV_Ith_S(y_out,lv_eint);

    outf << T <<"\t"<< Edi/mpv_nH/mpv_nH;
    outf <<"  "<< Edn/mpv_nH/mpv_nH;
    outf <<"  "<< Edpi/mpv_nH/mpv_nH <<"\n";
    T *=1.05;
  } while (T<min(1.0e9,EP->MaxTemperature));
  outf.close();  
#endif //Debug
  // ================================================================
  // ================================================================

  cout <<"MPv3: Constructor finished and returning.\n";
  cout <<"---------------------------------------------------------------------\n\n";
  return;
}


// ##################################################################
// ##################################################################


void MPv3::setup_local_vectors()
{
  //
  // This is in a function so it can be replaced in an inherited class.
  //
  nvl     = 2;    // two local variables to integrate
  N_extradata = 0;
  N_equations = 2;
  y_in  = N_VNew_Serial(N_equations);
  y_out = N_VNew_Serial(N_equations);
  lv_H0   = 0;    // x(H0) is the first element in the array
  lv_eint = 1;    // E_{int} is the second element.
  //cout<<"!!!!!!!!!!!!!!!!!! nvl="<<nvl<<"\n";
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
  if      (coord_sys==COORD_CRT && ndim==3) {
    diff_angle.resize(6);
    for (int v=0;v<6;v++) diff_angle[v] = 4.0*M_PI/6.0;
  }
  else if (coord_sys==COORD_CRT && ndim==2) {
    diff_angle.resize(4);
    for (int v=0;v<4;v++) diff_angle[v] = 2.0*M_PI/4.0;
  }
  else if (coord_sys==COORD_CRT && ndim==1) {
    diff_angle.resize(2);
    for (int v=0;v<2;v++) diff_angle[v] = 1.0;
  }
  else if (coord_sys==COORD_CYL && ndim==2) {
    diff_angle.resize(3);
    //
    // for each source in turn, we get its direction and set the angle accordingly.
    //
    int count=0;
    int dir=-1;
    for (int isrc=0; isrc<RS->Nsources; isrc++) {
      if (RS->sources[isrc].type==RT_SRC_DIFFUSE) {
        for (int v=0;v<ndim;v++) {
          if (RS->sources[isrc].pos[v] > 1.0e99) dir = 2*v+1;
          if (RS->sources[isrc].pos[v] <-1.0e99) dir = 2*v;
        }
        if (dir<0) rep.error("Diffuse source not at infinity!",isrc);
        //
        // if direction is in Z then angle is as for 3D, and if R+ then 4x3D values.
        //
        if (dir==ZNcyl || dir==ZPcyl) diff_angle[count] = 4.0*M_PI/6.0;
        else if (dir==RPcyl)          diff_angle[count] = 16.0*M_PI/6.0;
        else rep.error("Bad source direction",dir);
        count++;
      }
    }
    cout <<"Angles for diffuse sources: ["<<diff_angle[0]<<", ";
    cout <<diff_angle[1]<<", "<<diff_angle[2]<<"]\n";
  }
  else if (coord_sys==COORD_SPH && ndim==1) {
    //
    // for spherical symmetry the only diffuse source is at infinity
    //
    diff_angle.resize(1);
    diff_angle[0] = 4.0*M_PI;
  }
  else {
    rep.error("Unhandled coord-sys/ndim combination",ndim);
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
  if (s=="H1+___" || s=="HII__" || s=="H1+" || s=="HII") {
    return pv_Hp;
  }
  else {
    return -1;
  }
}


// ##################################################################
// ##################################################################


//
// Set the properties of a multifrequency ionising radiation source.
//
int MPv3::set_multifreq_source_properties(
      const struct rad_src_info *rsi
      )
{
  //
  // Some sanity checks:
  // - make sure source is multi-freq and ionising
  // - make sure Rstar and Tstar are positive and finite
  //
  if (rsi->effect!=RT_EFFECT_PION_MULTI)
    rep.error("Source is not multi-frequency!", rsi->id);
  if (rsi->Rstar<0 || !isfinite(rsi->Rstar))
    rep.error("Source has bad Rstar parameter", rsi->Rstar);
  if (rsi->Tstar<0 || !isfinite(rsi->Tstar))
    rep.error("Source has bad Tstar parameter", rsi->Tstar);

  //
  // mincol is the minimum Tau we care about, maxcol is the max Tau
  // we care about, Emax is the max. energy (in ergs) to integrate
  // to in the fit.  I think 1000eV is a bit excessive...
  //
  double mincol=1.0e-3, maxcol=1.0e6, Emax=1000.0*1.602e-12;
  int Nspl=50, Nsub=800;

#ifdef HE_INERT
  //
  // If He is always neutral, then we can't have any He-ionising
  // photons, so just integrate up to its first ionisation energy.
  //
  Emax = 24.59*1.602e-12;
#else
  //Emax = 54.41778*1.602e-12;  // assume nothing doubly-ionized He.
  Emax = 100.0*1.602e-12;  // This is better for cosmology RT tests.
#endif // HE_INERT
  //
  // Call the function in hydrogen_photoion.
  //
  Setup_photoionisation_rate_table(
              rsi->Tstar,
              rsi->Rstar*6.96e10,
              rsi->strength,
              mincol,  maxcol, Emax,  Nsub, Nspl);


  return 0;
}


// ##################################################################
// ##################################################################


double MPv3::get_temperature(
      const double nH, ///< nH (per c.c.)
      const double E, ///< E_int (per unit volume)
      const double xp  ///< x(H+)
      )
{
  //
  // returns gas temperature according to E=nkT/(g-1) with n=nH*(1.1+1.1*x),
  // appropriate for a gas with 10% Helium by number, and if He is singly ionised
  // whenever H is ionised (n_i=1.1n_H, n_e=1.1n_H).
  //
  return gamma_minus_one*E/(k_B*get_ntot(nH,xp));
}



// ##################################################################
// ##################################################################


double MPv3::get_ntot(
      const double nH, ///< nH
      const double xp  ///< x(H+) 
      )
{
  return
    (JM_NION+JM_NELEC*xp)*nH;
}


// ##################################################################
// ##################################################################


int MPv3::convert_prim2local(
      const pion_flt *p_in, ///< primitive vector from grid cell (length nv_prim)
      double *p_local
      )
{
  //
  // Set internal energy density, H+ fraction, and number density of H.
  //
  p_local[lv_eint] = p_in[PG]/(gamma_minus_one);
  p_local[lv_H0]   = 1.0-p_in[pv_Hp];
  mpv_nH = p_in[RO]/mean_mass_per_H;


#ifdef MPV3_DEBUG
  //
  // Check for negative ion fraction, and set to a minimum value if found.
  //
  if (p_local[lv_H0]>1.0) {
    cout <<"MPv3::convert_prim2local: negative ion fraction input: ";
    cout <<"x(H0)="<<p_local[lv_H0] <<", resetting to "<<Max_NeutralFrac<<"\n";
    p_local[lv_H0] = Max_NeutralFrac;
  }
  //
  // Check for bad values:
  //
  if (p_local[lv_H0]>1.01 || p_local[lv_H0]<-0.01) {
    cout <<"MPv3::convert_prim2local: bad ion fraction: ";
    cout <<"x(H0)="<<p_local[lv_H0] <<", resetting to [0,1]\n";
  }
#endif

  //
  // Set x(H0) to be within the required range (not too close to zero or 1).
  //
  p_local[lv_H0] = max(Min_NeutralFrac, min(Max_NeutralFrac, p_local[lv_H0]));

  //
  // Check for negative pressure (note this shouldn't happen, so we output a
  // warning) and set to 10K if we find it.
  //
  if (p_local[lv_eint]<=0.0) {
#ifdef MPV3_DEBUG
    cout <<"MPv3::convert_prim2local: negative pressure input: p=";
    cout <<p_local[lv_eint]<<", setting to "<<EP->MinTemperature<<"K.\n";
#endif
    p_local[lv_eint] = get_ntot(mpv_nH,p_in[pv_Hp])*k_B*EP->MinTemperature/(gamma_minus_one);
  }


#ifdef MPV3_DEBUG
  //
  // Check for NAN/INF
  //
  for (int v=0;v<2;v++) {
    if (!isfinite(p_local[v]))
      rep.error("INF/NAN input to microphysics",p_local[v]);
  }
  if (mpv_nH<0.0 || !isfinite(mpv_nH))
    rep.error("Bad density input to MPv3::convert_prim2local",mpv_nH);
#endif // MPV3_DEBUG
  
  return 0;
}


// ##################################################################
// ##################################################################


int MPv3::convert_local2prim(
      const double *p_local,
      const pion_flt *p_in, ///< input primitive vector from grid cell (length nv_prim)
      pion_flt *p_out       ///< updated primitive vector for grid cell (length nv_prim)
      )
{
  for (int v=0;v<nv_prim;v++) p_out[v] = p_in[v];

  p_out[PG]    = p_local[lv_eint]*(gamma_minus_one);
  p_out[pv_Hp] = 1.0-p_local[lv_H0];

#ifdef MPV3_DEBUG
  if (p_out[pv_Hp]<-10.0*JM_RELTOL || p_out[pv_Hp]>1.0*(1.0+JM_RELTOL) || !isfinite(p_out[pv_Hp])) {
    rep.printVec("p_in",p_in, nv_prim);
    rep.printVec("p_out",p_out, nv_prim);
    rep.printVec("p_local",p_local, nvl);
    rep.error("Bad output H+ value in MPv3::convert_local2prim",p_out[pv_Hp]-1.0);
  }
  if (p_out[PG]<0.0 || !isfinite(p_out[PG]))
    rep.error("Bad output pressure in MPv3::convert_local2prim",p_out[PG]);
#endif // MPV3_DEBUG

  //
  // Set xHp to be within the required range (not too close to zero or 1).
  //
  p_out[pv_Hp] = max(Min_NeutralFrac, min(Max_NeutralFrac, static_cast<double>(p_out[pv_Hp])));

  //
  // Set output pressure to be within required temperature range (use the 
  // possibly corrected x(H+) from p_out[]).
  //
  double T = get_temperature(mpv_nH, p_local[lv_eint], p_out[pv_Hp]);
  if (T>1.01*EP->MaxTemperature) {
    Set_Temp(p_out,EP->MaxTemperature,0);
#ifdef MPV3_DEBUG
    cout <<"MPv3::convert_local2prim() HIGH temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_Hp]<<"...  limiting to T="<<EP->MaxTemperature<<"\n";
#endif // MPV3_DEBUG
  }
  if (T<0.99*EP->MinTemperature) {
    Set_Temp(p_out,EP->MinTemperature,0);
#ifdef MPV3_DEBUG
    cout <<"MPv3::convert_local2prim() LOW  temperature encountered. ";
    cout <<"T="<<T<<", obtained from nH="<<mpv_nH<<", eint="<<p_local[lv_eint];
    cout <<", x="<<p_out[pv_Hp]<<"...  limiting to T="<<EP->MinTemperature<<"\n";
#endif // MPV3_DEBUG
  }

  return 0;
}


// ##################################################################
// ##################################################################

double MPv3::Temperature(
      const pion_flt *pv, ///< primitive vector
      const double      ///< eos gamma
      )
{
  //
  // Check for negative pressure/density!  If either is found, return -1.0e99.
  //
  if (pv[RO]<=0.0 || pv[PG]<=0.0) {
    //cout <<"MPv3::Temperature() negative rho="<<pv[RO]<<" or p="<<pv[PG]<<"\n";
    return -1.0e99;
  }
  //
  // generate vector of (nH,y(H0),Eint), and get Temperature from it.
  //
  double P[nvl];
  convert_prim2local(pv,P);
  return get_temperature(pv[RO]/mean_mass_per_H, P[lv_eint], 1.0-P[lv_H0]);
}


// ##################################################################
// ##################################################################


int MPv3::Set_Temp(
      pion_flt *p_pv,   ///< primitive vector.
      const double T, ///< temperature
      const double    ///< eos gamma.
      )
{
  //
  // Check for negative pressure.  If density<0 then we should bug
  // out because there is no way to set a temperature, but if p<0 we
  // can just overwrite it.
  //
  if (p_pv[PG]<=0.0) {
    //cout <<"MP_Hydrogen::Set_Temp() correcting negative pressure.\n";
    p_pv[PG] = 1.0e-12;  // Need p>0 for prim-to-local conversion.
  }
  double P[nvl];
  int err = convert_prim2local(p_pv,P);
  P[lv_eint] = get_ntot(mpv_nH, p_pv[pv_Hp])*k_B*T/(gamma_minus_one);
  err += convert_local2prim(P, p_pv, p_pv);
  return err;
}


// ##################################################################
// ##################################################################


int MPv3::TimeUpdateMP(
      const pion_flt *p_in,   ///< Primitive Vector to be updated.
      pion_flt *p_out,        ///< Destination Vector for updated values.
      const double dt,      ///< Time Step to advance by.
      const double,         ///< EOS gamma.
      const int,            ///< Switch for what type of integration to use.
      double *random_stuff  ///< Vector of extra data (column densities, etc.).
      )
{
  //
  // Call the new update function, but with zero radiation sources.
  //
  std::vector<struct rt_source_data> temp;
  int err = TimeUpdateMP_RTnew(p_in, 0, temp, 0, temp, p_out, dt, 0, 0, random_stuff);
  return err;
}


// ##################################################################
// ##################################################################


int MPv3::TimeUpdateMP_RTnew(
      const pion_flt *p_in, ///< Primitive Vector to be updated.
      const int N_heat,   ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &heat_src,
        ///< list of UV-heating column densities and source properties.
      const int N_ion,    ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &ion_src,
        ///< list of ionising src column densities and source properties.
      pion_flt *p_out,  ///< Destination Vector for updated values
                      ///< (can be same as first Vector.
      const double dt,   ///< Time Step to advance by.
      const double,   ///< EOS gamma.
      const int, ///< Switch for what type of integration to use.
                 ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
      double *random_stuff ///< final temperature (not strictly needed).
      )
{
  //
  // First set local variables for state vector and radiation source
  // properties.
  //
  int err=0;
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to MPv3::TimeUpdateMP_RTnew()",err);
  }
  setup_radiation_source_parameters(p_in, P, N_heat, heat_src, N_ion, ion_src);

  //
  // update radiation source properties, if needed (re-calculate
  // multi-frequency photoionisation rates if the source properties
  // have changed).
  // TODO: CODE THIS SOMEWHERE!! BUT MAYBE put this somewhere else -- 
  //       update the source properties when they change,
  //       and through a different interface function!!
  //

  for (int v=0;v<nvl;v++) NV_Ith_S(y_in,v) = P[v];
  //
  // Calculate y-dot[] to see if anything is changing significantly over dt
  //
  double maxdelta=0.0;
  err = ydot(0, y_in, y_out, 0);
  if (err) 
    rep.error("dYdt() returned an error in MPv3::TimeUpdateMP_RTnew()",err);
  for (int v=0;v<nvl;v++) {
    maxdelta = max(maxdelta, fabs(NV_Ith_S(y_out,v)*dt/NV_Ith_S(y_in,v)));
  }

  //
  // Now if nothing is changing much, just to a forward Euler integration.
  //
  if (maxdelta < EULER_CUTOFF) {
    for (int v=0;v<nvl;v++) {
      NV_Ith_S(y_out,v) = NV_Ith_S(y_in,v) + dt*NV_Ith_S(y_out,v);
    }
  }
  //
  // Otherwise do the implicit CVODE integration
  //
  else {
    err = integrate_cvode_step(y_in, 0, 0.0, dt, y_out);
    if (err) {
      rep.error("integration failed: MPv3::TimeUpdateMP_RTnew()",err);
    }
  }

  //
  // Now put the result into p_out[] and return.
  //
  for (int v=0;v<nvl;v++) P[v] = NV_Ith_S(y_out,v);
  err = convert_local2prim(P,p_in,p_out);
  return err;
}


// ##################################################################
// ##################################################################


double MPv3::timescales(
      const pion_flt *p_in, ///< Current cell state vector.
      const double,   ///< EOS gamma.
      const bool, ///< set to 'true' if including cooling time.
      const bool, ///< set to 'true' if including recombination time.
      const bool  ///< set to 'true' if including photo-ionsation time.
      )
{
#ifdef MPV3_DEBUG
  if (RS->Nsources!=0) {
    cout <<"WARNING: MPv3::timescales() using non-RT version!\n";
  }
#endif // MPV3_DEBUG
  std::vector<struct rt_source_data> temp;
  double tmin= timescales_RT(p_in, 0, temp, 0, temp, 0.0);
  temp.clear();
  return tmin;
}


// ##################################################################
// ##################################################################


///
/// This returns the minimum timescale of all microphysical processes, including
/// reaction times for each species and the total heating/cooling time for the gas.
/// It requires the radiation field as an input, so it has substantially greater
/// capability than the other timescales function.
/// Default setting is DT02, which limits by 0.25/ydot (and not by E/Edot)
///
double MPv3::timescales_RT(
      const pion_flt *p_in, ///< Current cell state vector.
      const int N_heat,      ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &heat_src,
      ///< list of UV-heating column densities and source properties.
      const int N_ion,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &ion_src,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      )
{
  int err=0;
  //
  // First convert to local variables.
  //
  double P[nvl];
  err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to MPv3::timescales_RT()",err);
  }
  NV_Ith_S(y_in,lv_H0  ) = P[lv_H0];
  NV_Ith_S(y_in,lv_eint) = P[lv_eint];

  //
  // Next set the radiation properties of the current cell.
  //
  setup_radiation_source_parameters(p_in,P, N_heat, heat_src, N_ion, ion_src);

  //
  // Now calculate y-dot[]...
  //
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt() returned an error in MPv3::timescales_RT()",err);
  }

  //
  // And finally get the smallest timescale over which things are varying.
  //
  double t=HUGEVALUE;
  //
  // First get the ionisation timescale, limited to dt = 0.25/|xdot|.
  // Tests have shown this is good enough, and that a restriction on the energy change 
  // (heating timescale) is not required for accurately tracking ionisation fronts 
  // (although it may be needed for cooling!).
  // For testing purposes there are ifdefs to allow the code to use a relative 
  // change in neutral fraction and/or the relative change in energy as the 
  // timestep criterion, rather than the default of absolute change in neutral
  // fraction.
  //
#ifdef USE_RELATIVE_NEUFRAC_DTLIMIT
  t = min(t,DTFRAC*max(5.0e-2,NV_Ith_S(y_in, lv_H0))/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"using neutral fraction.\n";
#else
  t = min(t,DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE));
  //cout <<"limit by dx: dt="<<DTFRAC/(fabs(NV_Ith_S(y_out, lv_H0))+TINYVALUE)<<"\n";
#endif


#ifdef MPV3_DEBUG
  if (t<3.16e9) {
  cout <<"MP timescales: xdot="<<NV_Ith_S(y_out, lv_H0);
  cout <<", Edot="<<NV_Ith_S(y_out, lv_eint)<<" t_x="<<t;
  }
#endif // MPV3_DEBUG

#ifdef ENERGY_CHANGE_TIMESTEP_LIMIT
  //
  // Now cooling/heating time to limit to X% change in energy).
  //
  t = min(t,DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE));
  //cout <<"using fractional energy change.\n";
  //cout <<"limit by dE: dt="<<DTFRAC*P[lv_eint]/(fabs(NV_Ith_S(y_out, lv_eint))+TINYVALUE)<<"\n";
#endif

#ifdef MPV3_DEBUG
  if (t<3.16e9) {
  cout <<" and min(t_x,t_e)="<<t<<",  "; rep.printVec("P[1-x,E]",P,nvl);
  }
#endif // MPV3_DEBUG
  return t;
}



// ##################################################################
// ##################################################################




double MPv3::total_cooling_rate(
      const pion_flt *p_in, ///< primitive input state vector.
      const int N_heat, ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &heat_src,
      ///< list of UV-heating column densities and source properties.
      const int N_ion,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &ion_src,
      ///< list of ionising src column densities and source properties.
      const double   ///< EOS gamma.
      )
{
  double P[nvl];
  int err = convert_prim2local(p_in,P);
  if (err) {
    rep.error("Bad input state to MPv3::total_cooling_rate()",err);
  }
  NV_Ith_S(y_in,lv_H0  ) = P[lv_H0];
  NV_Ith_S(y_in,lv_eint) = P[lv_eint];

  // set the radiation properties of the current cell.
  setup_radiation_source_parameters(p_in,P, N_heat, heat_src, N_ion, ion_src);
  // Now calculate y-dot[]...
  err = ydot(0, y_in, y_out, 0);
  if (err) {
    rep.error("dYdt(): error in MPv3::total_cooling_rate()()",err);
  }
  return NV_Ith_S(y_out, lv_eint);
}


// ##################################################################
// ##################################################################



double MPv3::get_recombination_rate(
      const int,          ///< ion index in tracer array (optional).
      const pion_flt *p_in, ///< input state vector (primitive).
      const double g      ///< EOS gamma (optional)
      )
{
#ifdef FUNCTION_ID
  cout <<"MPv3::get_recombination_rate()\n";
#endif // FUNCTION_ID
  double rate=0.0;
  double P[nvl];
  //
  // First convert to local variables.
  //
  convert_prim2local(p_in,P);
  //
  // Now get rate
  //
  rate = Hii_rad_recomb_rate(get_temperature(mpv_nH, P[lv_eint], 1.0-P[lv_H0]))
          *mpv_nH*mpv_nH *(1.0-P[lv_H0])*(1.0-P[lv_H0]) *JM_NELEC;

#ifdef FUNCTION_ID
  cout <<"MPv3::get_recombination_rate()\n";
#endif // FUNCTION_ID
  return rate;
}



// ##################################################################
// ##################################################################



void MPv3::setup_radiation_source_parameters(
      const pion_flt *p_in, ///< primitive input state vector.
      double *P,  ///< local input state vector (x_in,E_int)
      const int N_heat, ///< Number of UV heating sources.
      const std::vector<struct rt_source_data> &heat_src,
      ///< list of UV-heating column densities and source properties.
      const int N_ion,      ///< number of ionising radiation sources.
      const std::vector<struct rt_source_data> &ion_src
      ///< list of ionising src column densities and source properties.
      )
{
  //-------------------- RADIATION SOURCE INFO -----------------------
  //
  // First deal with the column densities and source strengths to get the UV heating
  // and EUV ionisation+heating rates.
  //
#ifdef MPV3_DEBUG
  if (heat_src.size() != static_cast<unsigned int>(N_heat)) {
    rep.error("Timescales: N_heating_srcs doesn't match vector size in MP3",
              heat_src.size());
  }
  if (ion_src.size() != static_cast<unsigned int>(N_ion)) {
    rep.error("Timescales: N_ionising_srcs doesn't match vector size in MP3",
              ion_src.size());
  }
#endif // MPV3_DEBUG

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
  mpv_delta_S  = 0.0;
  mpv_Vshell   = 1.0e200; 
  if (N_ion>0) {
    mpv_delta_S = ion_src[0].dS;
    mpv_Vshell  = ion_src[0].Vshell;
  }
  if (N_ion==0 && N_heat>0) {
    //
    // if no ionising sources, there may still be a point UV-heating source
    // which we need mpv_Vshell for.  here DelCol = rho*ds.
    //
    for (int v=0; v<N_heat; v++) {
      if (heat_src[v].type == RT_SRC_SINGLE) {
#ifdef MPV3_DEBUG
        if (P[lv_H0]>0.01) {
          cout <<"setup_rad_src_params: heating:  ds="<<heat_src[v].DelCol[0]/p_in[RO];
          cout <<", mpv_Vshell="<<heat_src[v].Vshell<<"\n";
        }
#endif // MPV3_DEBUG
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
  // Most H-ionising photons can't ionise He, so best that it doesn't contribute to the opacity.
  // BUT KEEP IN MIND THIS IS AN APPROXIMATION!
  // Column and DelCol are rho*ds*(1-x), and we want nH*ds*(1-x), so divide by mu_H
  //
  //
  if (N_ion>0) {
    mpv_Tau0  = ion_src[0].Column[0]*Hi_monochromatic_photo_ion_xsection(JUST_IONISED)/mean_mass_per_H;
    mpv_dTau0 = ion_src[0].DelCol[0]*Hi_monochromatic_photo_ion_xsection(JUST_IONISED)/mean_mass_per_H;
  }
  else {
    mpv_Tau0  = 0.0;
    mpv_dTau0 = 0.0;
  }
  //-------------------- mpv_Tau0 and mpv_dTau0 -----------------------


  //-------------------- mpv_G0_UV and mpv_G0_IR -----------------------
  //
  // Hard-coded to assume UV radiation opacity comes from dust with sigma=5e-22 cm2,
  // and that the input is \int \rho ds, integrated to the front edge of the cell.
  // Diffuse radiation strength is an intensity (i.e. per solid angle), so the total flux
  // going into the solver will be roughly: sum(I(Omega)*delta-Omega)exp(-1.9Av).
  // This comes from Henney et al. (2009) eq. A3
  // TODO: Check value of 1.9 for extinction, and the cross section of 5e-22
  //
  mpv_G0_UV = 0.0;
  mpv_G0_IR = 0.0;
  if (N_heat>0) {
    double temp=0.0;
    int i_diff=0;
    double Av_UV = 1.9*1.086*5.0e-22*METALLICITY/mean_mass_per_H;
    double Av_IR = Av_UV*0.05/1.9;

    for (int v=0; v<N_heat; v++) {
      if (heat_src[v].type == RT_SRC_DIFFUSE) {
        temp = heat_src[v].strength *diff_angle[i_diff];
#ifdef MPV3_DEBUG
        cout <<"setup_rad_src_params:\tdiffuse src: id="<<heat_src[v].id<<" 1.9Av="<<Av_UV*heat_src[v].Column[0];
        cout <<", strength="<<heat_src[v].strength<<", angle="<<diff_angle[i_diff];
        cout <<": attenuated flux="<<temp*exp(-Av_UV*heat_src[v].Column[0])<<"\n";
#endif // MPV3_DEBUG
        mpv_G0_UV += temp*exp(-Av_UV*heat_src[v].Column[0]);
        mpv_G0_IR += temp*exp(-Av_IR*heat_src[v].Column[0]);
        i_diff++;
#ifdef MPV3_DEBUG
        cout <<"UV_diff_flux="<<temp*exp(-Av_UV*heat_src[v].Column[0]);
        cout <<" Col="<<heat_src[v].Column[0]<<" Av="<<Av_UV*heat_src[v].Column[0]/1.9<<"\n";
#endif // MPV3_DEBUG
      }
      else {
        //
        // This source must be a point source of UV heating. In this case the strength is 
        // the photon luminosity, so flux = L*ds*exp(-1.9Av)/mpv_Vshell
        //
        //cout <<"heat_src[v].strength "<<heat_src[v].strength;
        temp = heat_src[v].strength*mpv_delta_S/heat_src[v].Vshell;
#ifdef MPV3_DEBUG
        cout <<"setup_rad_src_params:\tpoint   src: id="<<heat_src[v].id<<" 1.9Av="<<Av_UV*heat_src[v].Column[0];
        cout <<", strength="<<heat_src[v].strength<<", ds="<<mpv_delta_S;
        cout <<", mpv_Vshell="<<heat_src[v].Vshell;
        cout <<": attenuated flux="<<temp*exp(-Av_UV*heat_src[v].Column[0])<<"\n";
#endif // MPV3_DEBUG
        mpv_G0_UV += temp*exp(-Av_UV*heat_src[v].Column[0]);
        mpv_G0_IR += temp*exp(-Av_IR*heat_src[v].Column[0]);
#ifdef MPV3_DEBUG
        cout <<"UV_ptsc_flux="<<temp*exp(-Av_UV*heat_src[v].Column[0])<<"\n";
        cout <<"UV_ptsc_flux="<<temp*exp(-Av_UV*heat_src[v].Column[0]);
        cout <<" Col="<<heat_src[v].Column[0]<<" Av="<<Av_UV*heat_src[v].Column[0]/1.9<<"\n";
#endif // MPV3_DEBUG
      }
    } // loop over heating sources.
    //
    // now divide by 1.2e7 to get it normalised correctly for Will's equation A3.
    //
    mpv_G0_UV /= 1.2e7;
    mpv_G0_IR /= 1.2e7;
    //cout <<"  "<<mpv_G0_UV<<"  "<<mpv_G0_IR<<"\n";
#ifdef MPV3_DEBUG
    if (mpv_G0_UV>1.0) {
      cout <<"\tTotal UV attenuated flux = "<<mpv_G0_UV<<" in units of 1.2e7 phot/cm2/s\n";
    }
#endif // MPV3_DEBUG
  } // If there are UV heating sources
  else {
    mpv_G0_UV = 0.0;
    mpv_G0_IR = 0.0;
  }
  //-------------------- mpv_G0_UV and mpv_G0_IR -----------------------

#ifdef RT_TESTING
  //if (P[lv_H0]>0.01) {
      cout <<"MPv3: ionising: ds="<<mpv_delta_S<<", mpv_Vshell="<<mpv_Vshell;
      cout <<": mpv_Tau0="<<mpv_Tau0<<", mpv_dTau0="<<mpv_dTau0<<", nH="<<mpv_nH<<"\n";
  //}
#endif
  return;
}


// ##################################################################
// ##################################################################


int MPv3::ydot(
      double,               ///< current time (UNUSED)
      const N_Vector y_now, ///< current Y-value
      N_Vector y_dot,       ///< vector for Y-dot values
      const double *        ///< extra user-data vector (UNUSED)
      )
{
  //
  // fixes min-neutral-fraction to Min_NeutralFrac
  //
  double OneMinusX = max(NV_Ith_S(y_now,lv_H0),Min_NeutralFrac);
  double E_in      = NV_Ith_S(y_now,lv_eint);
  double x_in      = 1.0-OneMinusX;

  double ne        = JM_NELEC*x_in*mpv_nH;

  //
  // First get the temperature.  We assume the total particle number density
  // is given by 1.1*nH*(1+x_in), appropriate for a gas with 10% Helium by 
  // number, and if He is singly ionised whenever H is.
  //
  double T = get_temperature(mpv_nH, E_in, x_in);


  double temp1=0.0, temp2=0.0;
  double oneminusx_dot=0.0; // oneminusx_dot is in units of 1/s
  double Edot=0.0;
  // Edot is calculated in units of erg/s per H nucleon, multiplied by mpv_nH
  // at the end.

  //
  // We set a minimum electron density based on the idea that Carbon is singly
  // ionised in low density gas.  y(C)=1.5e-4 in the gas phase (by number)
  // (Sofia,1997), approximately, so I add this to the electron number density
  // with an exponential cutoff at high densities.
  //
  ne += mpv_nH*1.5e-4*METALLICITY*exp(-mpv_nH/1.0e4);


  //
  // collisional ionisation of H, with its associated cooling.
  // scales with n_e*nH0
  //
  Hi_coll_ion_rates(T, &temp1, &temp2);
  oneminusx_dot -= temp1*ne*OneMinusX; // the nH is divided out on both sides.
  Edot -= temp2*ne*OneMinusX;
  //cout <<"CI-CR="<< temp2*ne*OneMinusX<<"\n";

  //
  // photo-ionisation of H: photoionisation rate uses equation 18 in Mellema et
  // al. 2006 (C2-ray paper), noting that their Gamma is the rate per neutral
  // H, so we multiply by 1-x, as in their equation 11.
  //
  if (N_ion_srcs) {
    //
    // set current cell dTau0
    //
    temp1 = mpv_nH*mpv_delta_S*OneMinusX*
            Hi_monochromatic_photo_ion_xsection(JUST_IONISED);

    switch (ion_src_type) {
      case RT_EFFECT_PION_MULTI:
      //
      // Rather than divide the discretised rate by n(H0) and then multiply by
      // (1-x) to get oneminusx_dot, we simply divide by n(H) since this is
      // more numerically stable.  To do this, n(H) is passed to the rate
      // function instead of n(H0).
      //
      // Also, instead of using mpv_dTau0 for the cell optical depth, we use
      // the current optical depth (nH*mpv_delta_S*OneMinusX) so that we allow
      // the photoionisation rate to decrease as the number of neutral atoms
      // decreases during the timestep.  This is more stable.
      //
      oneminusx_dot -= Hi_discrete_multifreq_photoion_rate(mpv_Tau0, temp1,
                                          mpv_nH, mpv_delta_S, mpv_Vshell);
      Edot += Hi_discrete_multifreq_photoheating_rate(mpv_Tau0, temp1, mpv_nH,
                                          mpv_delta_S, mpv_Vshell);
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
      temp1 = Hi_discrete_mono_photoion_rate(mpv_Tau0, temp1, mpv_nH, mpv_NIdot, 
                                             PHOTON_ENERGY, mpv_delta_S, mpv_Vshell);
      oneminusx_dot -= temp1;
      Edot += temp1*EXCESS_ENERGY;
      //cout <<"PI-HR="<<temp1*EXCESS_ENERGY<<"\n";
      break;

      default:
      rep.error("Bad ion_src_type in dYdt()",ion_src_type);
      break;
    } // switch
  }

  //
  // radiative recombination of H+
  //
  oneminusx_dot += Hii_rad_recomb_rate(T) *x_in*ne;
  //
  // Total H+ cooling: recombination plus free-free
  //
  Edot -= Hii_total_cooling(T) *x_in*ne;
  //cout <<"HII-TC="<<Hii_total_cooling(T) *x_in*ne<<"\n";

  //
  // Add Helium free-free (Z^2*n(He)/n(H) = 0.25*X(He)/X(H) of the H+ free-free
  // rate) The normalisation is scaled so that I multiply by ne*nHp to get the
  // correct cooling rate (i.e. the abundance of He is included in the
  // prefactor).
  //
#ifndef HE_INERT
  // Only if He is ionised, otherwise it has no free-free.
  //
  Edot -= 1.68e-27*(JM_NION-1.0)*sqrt(T)*x_in*ne;
#endif // HE_INERT

  //
  // collisional excitation cooling of H0 Aggarwal (1983) and Raga+(1997,ApJS).
  //
  Edot -= Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10);
  //cout <<"CE-CR="<<Hi_coll_excitation_cooling_rate(T)*OneMinusX*ne *exp(-T*T/5.0e10)<<"\n";

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
    //cout <<"adding diffuse heating! ";
    Edot += 1.9e-26*METALLICITY*mpv_G0_UV/(1.0+6.4*(mpv_G0_UV/mpv_nH));
    //cout <<"
    //DfUV="<<1.9e-26*METALLICITY*mpv_G0_UV/(1.0+6.4*(mpv_G0_UV/mpv_nH));

    //
    // IR heating (HAdCM09 eq.A6) from point source and/or diffuse radiation.
    // There is a different G0 parameter because the attenuation is according
    // to exp(-0.05Av) rather than before where the coefficient was 1.9.
    //
    Edot += 7.7e-32*METALLICITY*mpv_G0_IR/pow(1.0+3.0e4/mpv_nH,2.0);
    //cout <<"DfIR="<<7.7e-32*METALLICITY*mpv_G0_IR/pow(1.0+3.0e4/mpv_nH,2.0)<<"\n";
  }

  //
  // X-ray heating (HAdCM09 eq.A5) Massive stars have x-ray luminosities of
  // ~1.0e32 erg/s, so use this.  NOT IMPLEMENTED YET... NEEDS SOMETHING MORE
  // CAREFUL.
  //
  ////Edot += 6.0e9*mpv_delta_S/mpv_Vshell;

  //
  // Cosmic ray heating (HAdCM09 eq.A7).
  //
  Edot += 5.0e-28*OneMinusX;
  //cout <<"CR-HR="<<5.0e-28*OneMinusX<<"\n";

  //
  // Cosmic Ray ionisation rate (Wolfire+,2003,eq.16) in solar neighbourhood.
  //
  oneminusx_dot -= 1.8e-17*OneMinusX;

  //
  // Diffuse UV Heating rate (Wolfire+,2003,eq.20,21, Fig.10,b).  This is a fit
  // to the curve in the top-right panel of Fig.10.
  //Edot += 1.66e-26*pow(mpv_nH,0.2602)*OneMinusX;
  //
  // This is a better function, using the first term of eq.19, with
  // phi_{PAH}=0.5 and G_0=1.7.  I multiply by the neutral fraction OneMinusX
  // because this heating term is only calculated for warm neutral medium.
  //
  // HACK!!!
#ifndef BETELGEUSE
  Edot += 1.083e-25*METALLICITY*OneMinusX/(1.0+9.77e-3*pow(sqrt(T)/ne,0.73));
  //cout<<"FUV-HR="<<1.083e-25*METALLICITY*OneMinusX
  //                 /(1.0+9.77e-3*pow(sqrt(T)/ne,0.73))<<"\n";
#endif // BETELGEUSE
  
//#ifdef BETELGEUSE
//  Edot -= METALLICITY*mpv_nH* (2.0e-19*exp(-1.184e5/(T+1.0e3)) +2.8e-28*sqrt(T)*exp(-92.0/T));
//#endif // BETELGEUSE

//#ifndef BETELGEUSE
  //
  // Now COOLING: First forbidden line cooling of e.g. OII,OIII, dominant in
  // HII regions.  This is collisionally excited lines of photoionised metals.
  // (HAdCM09 eq.A9) I have exponentially damped this at high temperatures!
  // This was important!  Oxygen abundance set to 5.37e-4 from
  // Asplund+(2009,ARAA), times 0.77 to account for 23% of O in solid dust.
  //
  temp1 = 1.20e-22*METALLICITY *exp(-33610.0/T -(2180.0*2180.0/T/T))
                               *x_in*ne*exp(-T*T/5.0e10);
  //
  // Fit to Raga, Mellema, Lundqvist (1997) rates for CNO if all are only
  // singly ionised, and for gas phase abundances of CNO of - n(C)/nH =
  // 2.95e-4*0.508 = 1.5e-4 (Sofia+,1997) - n(N)/nH = 7.41e-5 - n(O)/nH =
  // 5.37e-4*0.77 (0.23 goes in solids) These are taken from Asplund et al.
  // 2009.
  //
  //temp1 = 3.0e-22*METALLICITY*exp(-pow(1.4e5/T,0.6))
  //        *exp(-sqrt(mpv_nH/3.0e4)) *x_in*ne*exp(-T*T/5.0e10);


  //
  // Now the Wiersma et al (2009,MN393,99) (metals-only) CIE cooling curve.  We
  // take the actual cooling rate to be the max of SD93-CIE and the previous
  // two terms.
  //
  temp2 = cooling_rate_SD93CIE(T) *x_in*x_in*mpv_nH*METALLICITY;
  Edot -= max(temp1,temp2);


  //
  // Instead of the PDR cooling from Henney, use Wolfire's eq.C1,C3 for
  // collisional cooling of CII and OI by neutral H atoms.  In eq.C3 I have
  // absorbed the (100K)^{-0.4} into the prefactor, and used x^a=exp(a*ln(x)).
  // I have cut off equation C1 at high densities to be consistent with the ion
  // fraction of C that I assumed above for the electron density.
  //
  Edot -= 3.15e-27*METALLICITY*exp(-92.0/T)*mpv_nH*OneMinusX*exp(-mpv_nH/1.0e4);
  Edot -= 3.96e-28*METALLICITY*exp(0.4*log(T)-228.0/T)*mpv_nH*OneMinusX;

  //
  // This is the CII cooling by electron collisions, cutoff at high density
  // again, for consistency, again with sqrt(100K) absorbed into the prefactor.
  // This rate has a very low critical density (Goldsmith, Langer et al.,
  // 2012ApJS..203...13G), at n_c=20 cm^{-3} at 1000K, so we use their
  // temperature scaling and divide by density according to rate = rate/(1.0 +
  // 0.05*nH*(T/2000K)^(-0.37))
  //
  Edot -= 1.4e-23*METALLICITY*exp(-0.5*log(T)-92.0/T)*ne
          *exp(-mpv_nH/1.0e4)
#ifndef BETELGEUSE
       //   /(1.0 + 0.05*mpv_nH*pow(T/2000.0,-0.37))
#endif // BETELGEUSE
          ;

  //
  // PAH cooling: eq. 21 in Wolfire+,2003.  I think they should have multiplied
  // their equation by 1.3 for the increased PAH abundance...
  //
#ifndef BETELGEUSE
  //Edot -= 2.325e-30*METALLICITY*exp(0.94*log(T) +0.74*pow(T,-0.068)*log(3.4*sqrt(T)/ne))*ne;
  Edot -= 3.02e-30*METALLICITY*exp(0.94*log(T) +0.74*pow(T,-0.068)*log(3.4*sqrt(T)/ne))*ne;
#endif // BETELGEUSE
//#endif // BETELGEUSE

  //
  // now multiply Edot by nH to get units of energy loss/gain per unit volume
  // per second.
  //
  Edot *= mpv_nH;
#ifdef HIGHDENS_CUTOFF
  if (Edot<0.0) Edot *= exp(-mpv_nH*mpv_nH/1.0e6);
#endif //HIGHDENS_CUTOFF 

  //
  // We want to limit cooling as we approach the minimum temperature, so we
  // scale the rate to linearly approach zero as we reach Tmin.
  //
  if (Edot<0.0 && T<2.0*EP->MinTemperature) {
    Edot = min(0.0, (Edot)*(T-EP->MinTemperature)/EP->MinTemperature);
  }

  NV_Ith_S(y_dot,lv_H0)   = oneminusx_dot;
  NV_Ith_S(y_dot,lv_eint) = Edot;
  return 0;
}


