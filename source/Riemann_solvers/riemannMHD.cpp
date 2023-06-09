/// \file riemannMHD.cc
/// \brief Class Member function definitions for linear MHD Riemann Solver.
///
///  \author Jonathan Mackey
///
/// Written 2006-11-20
/// Modified :
///  - 2007-02-01 Documented, added new eigenvector calculation.
///  - 2007-07-13 New Class members for setdirection and setavgstate.
///  - 2007-07-16 fixed problem in solver2codevars() introduced by
///  SetDirection() on 13/7.
///  - 2008-03-03 Fixed a few minor bugs.
///  - 2009-10-24 New class structure so that it inherits from the MHD equations
///  class.
///  - 2010-02-10 JM: Replaced exit() calls with spdlog::error("{}: {}", ) calls
///  so parallel code exits cleanly.
///
/// This is a linear MHD Riemann Solver for the ideal MHD equations.
/// It is a "Roe-type" solver (but using an arithmetic aveage instead
/// of a Roe-average, following Falle, Kommisarov and Joarder (1998)),
/// linearising the Jacobian matrix with an average state and jumping
/// across waves from the left and/or right to get to the state at x=0.
///
/// - 2010.10.01 JM: Cut out testing myalloc/myfree
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.23 JM: Moved to Riemann_solvers/ directory.
///   Made RS_left/right/pstar/meanp private local data to avoid
///   confusing inheritance issues.  Tidied up code.  Put some tests
///   for finiteness in RS_TESTING ifdef.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "riemannMHD.h"
using namespace std;
//#define RS_TESTING ///< Comment this out if not testing the solver; will speed
// things up.

// ##################################################################
// ##################################################################

/// operators for the enums.
waves &operator++(waves &d)
{
  return (d = waves(d + 1));
}
conserved &operator++(conserved &d)
{
  return (d = conserved(d + 1));
}
rsvars &operator++(rsvars &d)
{
  return (d = rsvars(d + 1));
}

// ##################################################################
// ##################################################################

//
// Constructor of the riemann_MHD class.
//
riemann_MHD::riemann_MHD(
    const int nv,           ///< length of state vectors.
    const pion_flt *state,  ///< reference vector.
    const double g          ///< equation of state gamma.
    ) :
    eqns_base(nv),
    eqns_mhd_ideal(nv)
{
#ifdef RS_TESTING
  spdlog::debug("(riemann_MHD::riemann_MHD) Initialising Riemann Solver Class");
  if (eq_nvar < 8) {
    spdlog::error(
        "\tProblem with MHD Riemann Solver... Nelements!=8.  {}", eq_nvar);
    exit(1);
  }
#endif

  RS_nvar = 7;  // Local state vector has length 7, because BX is not included.
  eq_gamma =
      g;  // gamma is passed to solve() function, which assigns it each time.

  RS_evalue.resize(RS_nvar);
  RS_pdiff.resize(RS_nvar);
  RS_strength.resize(RS_nvar);

  RS_left.resize(eq_nvar);
  RS_right.resize(eq_nvar);
  RS_meanp.resize(eq_nvar);
  RS_pstar.resize(eq_nvar);

  //
  // zero all the arrays
  //
  for (int v = 0; v < RS_nvar; v++) {
    RS_pdiff[v] = RS_evalue[v] = RS_strength[v] = 0.;
    for (int w = 0; w < RS_nvar; w++)
      RS_leftevec[v][w] = RS_rightevec[v][w] = 0.;
  }
  for (int v = 0; v < eq_nvar; v++)
    eq_refvec[v] = RS_pstar[v] = RS_left[v] = RS_right[v] = RS_meanp[v] = 0.0;
  spdlog::debug(" Done");

  onaxis = offaxis = 0;
  samestate        = 0;
  totalsolve       = 0;
  smallB = MACHINEACCURACY;  // This is a small number to use when roundoff
                             // errors happen.
  tinyB = smallB * smallB * smallB;  // Set it to be really insignificant.

  SetDirection(XX);  // Initially assume we are looking along the x-axis.

  SetAvgState(state, eq_gamma);

  spdlog::debug("(riemann_MHD::riemann_MHD) Finished Setup");
}

// ##################################################################
// ##################################################################

// Destructor.
riemann_MHD::~riemann_MHD()
{
  spdlog::debug("(riemann_MHD::riemann_MHD) Class Destructing!!!");
  spdlog::debug("\t Failed Tests: OnAxis: {} and OffAxis: {}", onaxis, offaxis);
}

// ##################################################################
// ##################################################################

//
// Public interface to the riemann_MHD solver -- it takes in a left
// and right state, and then sets the solver going, which returns the
// state P* in 'result'.  err_code solve(left, right, result );
//
int riemann_MHD::JMs_riemann_solve(
    const pion_flt *l,
    const pion_flt *r,
    pion_flt *ans,
    const int mode,  ///< Solve Type (1=LinearRS,2=ExactRS,3=HybridRS, 4=RoeRS)
    const double g)
{
  int err = 0;

  //
  // First check that we know how to do the solver which
  // was requested:
  //
  RS_mode = mode;
  switch (RS_mode) {
    case (FLUX_RSlinear):
      // cout <<"\tMODE 1: All Linear Solves.\n";
      break;
    default:
      spdlog::error("\tMODE i: Don't know what to do {}", RS_mode);
      exit(2);
      break;
  }
  totalsolve++;

  eq_gamma = g;

  //
  // Assign class variables for left and rights states from values passed
  //
  assign_data(l, r);

  //
  // Set the mean state vector:
  //
  get_average_state();

  //
  // BX: B_x is a parameter in the riemann solver. I set it to be the
  // mean of the left and right states.  The mixed-GLM solver sets
  // both left and right values of BX to be the resolved state from
  // the GLM method solution.
  //
  RS_pstar[RBX] = ansBX = RS_meanp[RBX];
  // cout <<"bx(in)="<<RS_pstar[RBX];

#ifdef RS_TESTING
  // cout <<"meanp[RBX]="<<RS_meanp[RBX];
  // cout <<"  and 0.5(l+r)="<<0.5*(l[RBX]+r[RBX])<<"\n";
  for (int i = 0; i < eq_nvar; i++) {
    if (!isfinite(RS_left[i]) || !isfinite(RS_right[i])) {
      spdlog::debug(
          "MHD Left  : {}", std::vector<double>(RS_left, RS_left + eq_nvar));
      spdlog::debug(
          "MHD Right : {}", std::vector<double>(RS_right, RS_right + eq_nvar));
      return (-100);
    }
  }
#endif  // RS_TESTING

  //
  // First test if the left and right states are the same, and if they
  // are then return the mean state (left[]+right[])/2
  //
  double diff = 0.;
  for (int i = 0; i < RS_nvar; i++)
    diff += fabs(RS_right[i] - RS_left[i]) / (fabs(eq_refvec[i]) + TINYVALUE);

  if (diff < 1.e-6) {
    samestate++;  // cout <<"same states: "<<samestate<<"\n";
    //
    // Assign pstar values to be the mean of left and right states.
    //
    for (int v = 0; v < RS_nvar; v++)
      RS_pstar[v] = RS_meanp[v];

#ifdef RS_TESTING
    //
    // Check for negative pressure/density!  This should never happen
    // if we are just taking the mean of left and right states!
    //
    if (RS_pstar[RPG] < 0.) {
      RS_pstar[RPG] = eq_refvec[RPG] * BASEPG;  // BASEPG is 1.e-8
      spdlog::error(
          "(reimannMHD::solve) Negative pressure in INPUT STATES. pstar[] = [");
      for (int v = 0; v < (eq_nvar - 1); v++)
        spdlog::error("{}, ", RS_pstar[v]);
      spdlog::error("{}]", RS_pstar[eq_nvar - 1]);
      spdlog::error("Negative pressure in INPUT STATES {}", RS_pstar[RPG]);
      exit(1);
    }
    if (RS_pstar[RRO] < 0.) {
      RS_pstar[RRO] = eq_refvec[RRO] * BASEPG;  // BASEPG is 1.e-8
      spdlog::error(
          "(reimannMHD::solve) Negative density in INPUT STATES!. pstar[] = [");
      for (int v = 0; v < (eq_nvar - 1); v++)
        spdlog::error("{}, ", RS_pstar[v]);
      spdlog::error("]", RS_pstar[eq_nvar - 1]);
      spdlog::error("Negative density in INPUT STATES {}", RS_pstar[RRO]);
      exit(3);
    }
#endif  // RS_TESTING

    //
    // Re-order velocity and field componenets for returning to code
    //
    // solver2codevars(RS_left); // don't care about this -- private data.
    // solver2codevars(RS_right);// don't care about this -- private data.
    solver2codevars(&RS_pstar[0]);
    for (int v = 0; v < eq_nvar; v++)
      ans[v] = RS_pstar[v];
    // cout <<"  SS  bx(out)="<<ans[eqBX]<<"\t";
    return 0;
  }  // same state finished.

  //*******************************************************************
  // At this point, we need to solve a non-trivial Riemann problem.
  // Now we do the main part of the solution.
  // Follows S.2.1 in Falle, Komissarov & Joarder, 1998, MNRAS, 297, 265.
  // ******************************************************************
  switch (mode) {
    case (1):  // Should push this out into a function called LinearSolve()
               // if I ever add a second mode.

      err += get_sound_speeds();
      if (err != 0) {
        spdlog::error(
            "riemann_MHD::(get_sound_speeds) returned with error", err);
        exit(4);
      }
      //    cout <<"Speeds: (ch2,ca2,ct2,cs2,cf2): ("<<ch*ch
      //      << " ,"<<ca*ca<<" ,"<< bt*bt<< " ,"
      //      <<cs*cs<< " ,"<<cf*cf<<" )\n";

      get_eigenvalues();

      err += RoeBalsara_evectors();
      if (err != 0) {
        spdlog::error(
            "riemann_MHD::(RoeBalsara_evectors) returned with error", err);
        exit(5);
      }

#ifdef RS_TESTING
      //
      // This makes sure the eigenvectors are normalised and have the
      // correct orthogonality properties.
      //
      err += check_evectors();
      if (err != 0) {
        spdlog::error(
            "riemann_MHD::(check_evectors): ERROR code: {} Failing...\nSpeeds: (ch2,ca2,ct2,cs2,cf2): ({} ,{} ,{} ,{} ,{} )\n",
            err, ch * ch, ca * ca, bt * bt, cs * cs, cf * cf);
        spdlog::error(
            "Left state  : {}",
            std::vector<double>(RS_left, RS_left + eq_nvar));
        spdlog::error(
            "Right state : {}",
            std::vector<double>(RS_right, RS_right + eq_nvar));
        spdlog::error("Pstar : {}", RS_pstar);
        spdlog::error("{}: {}", "riemann_MHD::(check_evectors)", err);
        exit(6);
      }
#endif

      calculate_wave_strengths();

      err += get_pstar();
      if (err != 0) {
        spdlog::error("riemann_MHD::(get_pstar)  returned with error", err);
        exit(err);
      }
      //  cout << "(riemann_MHD::solve) Got P*..." << "\n";
      //  cout << "P* (Got Solution) rho: " << RS_pstar[RO];
      //  cout << "  v: " <<  RS_pstar[VX];
      //  cout << "   p: " << RS_pstar[PG] << "\n";

      //
      // Now we have gone through all these, so if we have a solution it
      // should be the right one!  We could have negative pressure from
      // very strong rarefactions though.
      //
      if (RS_pstar[RPG] < 0.) {
        // cerr << "(reimannMHD::solve) Negative pressure... Returning
        // "; cerr <<"vacuum conditions, leaving velocity unchanged.\n";
        // rep.printVec("MHD Left ",RS_left,eq_nvar);
        // rep.printVec("MHD Right",RS_right,eq_nvar);
        // rep.printVec("MHD Pstar",RS_pstar,eq_nvar);
        RS_pstar[RPG] = eq_refvec[RPG] * BASEPG;
      }
      if (RS_pstar[RRO] < 0.) {
        RS_pstar[RRO] = eq_refvec[RRO] * BASEPG;
        // cerr << "(reimannMHD::solve) Negative density... Returning ";
        // cerr <<"vacuum conditions, leaving velocity unchanged.\n";
        // rep.printVec("MHD Left ",RS_left,eq_nvar);
        // rep.printVec("MHD Right",RS_right,eq_nvar);
        // rep.printVec("MHD Pstar",RS_pstar,eq_nvar);
        // rep.printVec("Evalues: ",RS_evalue,RS_nvar);
      }
      //  cout << "(riemann_MHD::solve) Success!" << "\n";
      break;  // END CASE MODE=1

    default:
      spdlog::debug("MODE not known.  Only know 1.  Please enter a valid mode");
      spdlog::error("{}: {}", "Bad solve mode in riemann_MHD", mode);
      exit(7);
      break;
  }

#ifdef RS_TESTING
  if ((!pconst.equalD(RS_pstar[RBZ], 0.0) && eq_dir == XX)
      || (!pconst.equalD(RS_pstar[RBY], 0.0) && eq_dir == YY)) {
    spdlog::debug("*************** eqBBZ = {}  BBZ={}", eqBBZ, BBZ);
    spdlog::debug(
        "left :  : {}", std::vector<double>(RS_left, RS_left + eq_nvar));
    spdlog::debug(
        "right:  : {}", std::vector<double>(RS_right, RS_right + eq_nvar));
    spdlog::debug(
        "meanp:  : {}", std::vector<double>(RS_meanp, RS_meanp + eq_nvar));
    spdlog::debug(
        "pstar:  : {}", std::vector<double>(RS_pstar, RS_pstar + eq_nvar));
    // rep.printVec("flux : ",FS_flux ,eq_nvar);
    spdlog::debug(
        "strength:  : {}",
        std::vector<double>(RS_strength, RS_strength + eq_nvar));
    spdlog::debug(
        "evalues :  : {}", std::vector<double>(RS_evalue, RS_evalue + eq_nvar));
    spdlog::debug(
        "pdiff :  : {}", std::vector<double>(RS_pdiff, RS_pdiff + eq_nvar));
    spdlog::debug("Bx={}\t alphaf,s={} {}", ansBX, alphaf, alphas);
    check_evectors();
  }
#endif  // RS_TESTING

  //
  // Need to get pstar back into primitive var. form before returning
  // from solve.
  //
  solver2codevars(&RS_pstar[0]);
  for (int v = 0; v < eq_nvar; v++)
    ans[v] = RS_pstar[v];
    // cout <<"  bx(out)="<<ans[eqBX]<<"\t";

#ifdef RS_TESTING
  //
  // Finally, check for NAN/INF
  //
  for (int i = 0; i < eq_nvar; i++) {
    if (!isfinite(RS_left[i]) || !isfinite(RS_right[i])
        || !isfinite(RS_pstar[i])) {
      spdlog::debug(
          "MHD Left  : {}", std::vector<double>(RS_left, RS_left + eq_nvar));
      spdlog::debug(
          "MHD Right : {}", std::vector<double>(RS_right, RS_right + eq_nvar));
      spdlog::debug(
          "MHD Pstar : {}", std::vector<double>(RS_pstar, RS_pstar + eq_nvar));
      return (-100);
    }
  }
#endif

  return (0);
}

// ##################################################################
// ##################################################################

void riemann_MHD::assign_data(const pion_flt *l, const pion_flt *r)
{
  //
  // The Riemann problem has different state vector ordering, so we
  // first copy the vectors to local arrays, and then re-order the
  // local copies.
  //

#ifdef RS_TESTING
  if (l[eqRO] < TINYVALUE || l[eqPG] < TINYVALUE || r[eqRO] < TINYVALUE
      || r[eqPG] < TINYVALUE) {
    spdlog::error(
        "riemann_MHD::assign_data() Density/Pressure too small {} {}", r[eqRO],
        r[eqPG]);
    exit(1);
  }
#endif  // RS_TESTING

  //
  // Copy each element in turn, and then re-order using a temp-array
  // in code2solvervars()
  //
  for (int v = 0; v < eq_nvar; v++) {
    RS_left[v]  = l[v];
    RS_right[v] = r[v];
  }
  code2solvervars(&RS_left[0]);
  code2solvervars(&RS_right[0]);
  return;
}

// ##################################################################
// ##################################################################

void riemann_MHD::code2solvervars(pion_flt *statevec)
{
  //
  // Need to re-order velocity and B-field elements.
  // The Riemann Solver should never access elements with
  // v>7, so we can ignore them.

  pion_flt temp[8];
  temp[0]       = statevec[eqRO];
  temp[1]       = statevec[eqPG];
  temp[2]       = statevec[eqVX];
  temp[3]       = statevec[eqVY];
  temp[4]       = statevec[eqVZ];
  temp[5]       = statevec[eqBX];
  temp[6]       = statevec[eqBY];
  temp[7]       = statevec[eqBZ];
  statevec[RRO] = temp[0];
  statevec[RPG] = temp[1];
  statevec[RVX] = temp[2];
  statevec[RVY] = temp[3];
  statevec[RVZ] = temp[4];
  statevec[RBX] = temp[5];
  statevec[RBY] = temp[6];
  statevec[RBZ] = temp[7];

  return;
}

// ##################################################################
// ##################################################################

void riemann_MHD::solver2codevars(pion_flt *statevec)
{
  pion_flt temp[8];
  // Density, pressure.
  // Velocities
  // B-Field
  temp[0] = statevec[RRO];
  temp[1] = statevec[RPG];
  temp[2] = statevec[RVX];
  temp[3] = statevec[RVY];
  temp[4] = statevec[RVZ];
  temp[5] = statevec[RBX];
  temp[6] = statevec[RBY];
  temp[7] = statevec[RBZ];

  statevec[eqRO] = temp[0];
  statevec[eqPG] = temp[1];
  statevec[eqVX] = temp[2];
  statevec[eqVY] = temp[3];
  statevec[eqVZ] = temp[4];
  statevec[eqBX] = temp[5];
  statevec[eqBY] = temp[6];
  statevec[eqBZ] = temp[7];

  //
  // Tracers and any variables with v>7 are unaffected by the MHD solver!
  //

  return;
}

// ##################################################################
// ##################################################################

void riemann_MHD::failerror(int err, string text)
{
  if (err != 0) {
    spdlog::error("riemann_MHD::({}): ERROR code: {} Exiting...\n", text, err);
    spdlog::error("{}: {}", text, err);
    exit(err);
  }
}

// ##################################################################
// ##################################################################

void riemann_MHD::get_average_state()
{
  //
  // We can construct the average state any way we want.  Falle et
  // al. (1998) advocate just the arithmetic mean of the left and
  // right states.
  //

  // ARITHMETIC MEAN:
  for (int v = 0; v < eq_nvar; v++)
    RS_meanp[v] = 0.5 * (RS_left[v] + RS_right[v]);

  return;
}

// ##################################################################
// ##################################################################

int riemann_MHD::get_sound_speeds()
{
  //
  // The equations for these are on p.2 of Falle, Komissarov, & Joarder, 1998.
  // The Average State sound speeds are calculated based on the average state,
  // not on the left and right sound speeds.  I think this is sensible, and it
  // is what Andy does.
  // Hydrodynamic sound speeds
  //

  riemann_MHD::ch = sqrt(eq_gamma * RS_meanp[RPG] / RS_meanp[RRO]);
  // MHD speeds:
  riemann_MHD::bx = ansBX / sqrt(RS_meanp[RRO]);
  riemann_MHD::ca = fabs(bx);

  /** \section stability Numerical Stability
   *
   * The Fast and Slow speeds are subject to roundoff error.
   *   - First the term \f$ [(a^2+b^2)^2 - 4a^2b_x^2] \f$ can be very small if
   *     \f$b_t\f$ is very small and \f$b_x^2 \simeq a^2\f$.  So we need to
   * check that this is positive.
   *   - Secondly for the slow speed, the term \f$[(a^2+b^2)-\sqrt{\ldots}]\f$
   * can be very small.  This happens when \f$a \ll 1\f$ or when \f$b_x
   * \ll1\f$. So we also need to check that this is positive. The question is
   * what to do when these turn out to be negative.  This means that the
   * computer can't calculate it accurately.  I can give it an approximation
   * to calculate it accurately, or I can arbitrarily set it to be a value
   * slightly larger than the machine precision.  I decided the simplest thing
   * is the latter.
   *
   * In two cases where the slow speed goes to zero (namely \f$a\rightarrow
   * 0\f$ and/or \f$b_x\rightarrow 0\f$), it is well approximated by \f[ c_s^2
   * \simeq \frac{a^2b_x^2}{a^2+b^2} + \mbox{h.o.t.} \f] so I can use this.
   * There is also the hydrodynamic limit, where \f$b\rightarrow 0\f$. In this
   * case \f$ c_s^2 \simeq b_x^2
   * +\mbox{h.o.t.}\f$.
   *
   * In cases where \f$ [(a^2+b^2)^2 - 4a^2b_x^2]\rightarrow 0 \f$, the fast
   * and slow speed become the same.
   *
   * So at this point I will determine whether we have any numerical
   * instability present, and will get rid of it.
   *
   * Three quantities determine the stability: \f$ \{b_x,b_t,a\}\f$.
   * - if \f$a\rightarrow 0\f$ Then we are in trouble.  In this case MHD
   * probably isn't applicable, because the magnetic field energy is \f$
   * 10^{15} \f$ times stronger than the gas energy.  This is not likely to
   * happen in the problems we are considering.  So I think if this is
   * encountered we should bug out of the code, with an explanation of why, of
   * course.  So for the rest of the cases we assume that \f$a\f$ is of order
   * unity.
   * - if \f$b_x\rightarrow 0\f$ but \f$b_t\simeq 1\f$, then the Alfv\'en and
   * slow waves disappear into the contact discontinuity.  The fast and slow
   * speeds are given by \f$c_f^2 \simeq a^2+b_t^2\f$ and \f$c_s^2 \simeq
   * \frac{a^2b_x^2}{a^2+b^2}\f$ respectively.  The normalisation quantities
   * \f$\alpha_{f,s}\f$ will be calculated without difficulty.  So in this
   * case I may want to approximate \f$c_s\f$.  Or, what would be better would
   * be to have a pseudo-hydro solver, as Falle points out that this is ok.
   * Certainly simpler than solving the MHD problem.
   * - if \f$b_x\rightarrow 0 \f$ and \f$b_t\rightarrow 0\f$ we get the
   * hydrodynamic limit, where the magnetic field is irrelevant.  I can just
   * call the hydro solver in this case, and leave the magnetic field
   * unchanged.
   * - The remaining case is where \f$b_t\rightarrow 0\f$ and \f$b_x\simeq
   * 1\f$. This is really three subcases.  Obviously the Alfv\'en waves need
   * attention regardless.  Falle suggests setting \f$\beta_{f,s} =
   * 1/\sqrt{2}\f$, saying that it really doesn't matter what you choose.  I
   * wonder is it better to give them random signs, so you don't have any
   * chance to add up error???  For the slow and fast waves here are the three
   * cases:
   *   - When \f$b_x^2<a^2\f$, there is no problem as long as \f$a^2-b_x^2 >
   * \epsilon\f$ where \f$\epsilon\f$ is a small number close to the machine
   * precision.
   *   - When \f$b_x^2>a^2\f$, there is no problem as long as \f$b_x^2-a^2 >
   * \epsilon\f$ where \f$\epsilon\f$ is a small number close to the machine
   * precision.
   *   - When \f$b_x^2 \simeq a^2\f$, things become difficult.  The slow and
   * fast speeds become the same.  The fast and slow normalisations can vary
   * from zero to one very sensitively, based on the small difference in the
   * sound speeds.  Falle, and Roe and Balsara both claim that it doesn't
   * matter what you choose for this case, that the flux you get out at the
   * end is insensitive to what you choose.  I should test this.
   *
   * So, basically, if the slow or fast speeds are the same as the sound
   * speed, then I set them to be chydro+-smallB.  If the tangential field is
   * less than smallB^3, I set it to be tinyB.  If alpha_f,s cannot be
   * determined accurately, it means that the speeds are very similar, so I
   * set them to both to be equal at 1/root2.
   * */

  /** \subsection bt Tangential Field
   * The tangential field absolute value is calculated as
   * \f$b_t = \sqrt{(B_y^2+B_z^2)/\rho}\f$.  The two normalised components of
   * the field are given by \f$\beta_{y,z} = B_{y,z}/(\sqrt{\rho}b_t)\f$ if
   * \f$b_t > \mbox{\tt tinyB}\f$, and \f$\beta_{y,z} = 1/\sqrt{2}\f$
   * otherwise.
   * */
  riemann_MHD::bt = sqrt(
      (RS_meanp[RBY] * RS_meanp[RBY] + RS_meanp[RBZ] * RS_meanp[RBZ])
      / RS_meanp[RRO]);
  if (bt > tinyB) {
    betay = RS_meanp[RBY] / sqrt(RS_meanp[RRO]) / bt;
    betaz = RS_meanp[RBZ] / sqrt(RS_meanp[RRO]) / bt;
  }
  else {
    // Set betayz to be equal, and the tangential field is at 45deg.
    betay = 1. / sqrt(2.);
    betaz = 1. / sqrt(2.);
  }

  /** \subsection checka Check that 'a' is not tiny.
   * We check that the hydro sound speed is not insignificantly small compared
   * to the magnetic field strength.  If this is the case, then MHD is
   * probably no good anyway, and nothing will be calculated accurately.\n The
   * condition is that if \f$a/\max{(|b_x|,b_t)} < \sqrt{\mbox{machine
   * precision}}\f$ then we bug out and complain.
   * */
  if ((ch / max(ca, bt)) < sqrt(smallB)) {
    spdlog::error(
        "(riemann_MHD::get_sound_speeds) hydro sound speed insignificantly small compared to magnetic speeds. Bugging out...\n\t ch = {} and c_a = {}, c_t = {}",
        ch, ca, bt);
    return (1);
  }

  /** \subsection cfcs Fast and Slow Speeds.
   * At the moment, if I have roundoff errors in the calculation in the two
   * subtractions that I need to do, I just set the difference to be
   * smallB squared (as it always gets square-rooted in the next line).  So
   * this way I am deliberately adding an error equal to a number slightly
   * larger than the machine precision.\n
   * Also, I check that the fast/slow speed is larger/smaller than the
   * hydrodynamic sound speed, and if it's not I set it to be larger by
   * smallB.
   * */
  double temp1 = ch * ch + bx * bx + bt * bt;
  double temp2 = 4. * ch * ch * bx * bx;
  if ((temp2 = temp1 * temp1 - temp2) < MACHINEACCURACY)
    temp2 = MACHINEACCURACY;  // This is as good as the computer can get.
  riemann_MHD::cf = sqrt((temp1 + sqrt(temp2)) / 2.);
  if ((temp2 = temp1 - sqrt(temp2)) < MACHINEACCURACY)
    temp2 = MACHINEACCURACY;  // Again, as good as the machine can get.
  riemann_MHD::cs = sqrt(temp2 / 2.);

  if (cs > ch) {
    // cout<<"cs>ch: "<<cs/ch-1.<<"!!!\n";
    cs = ch - smallB;
  }
  if (ch > cf) {
    // cout<<"ch>cf: "<<ch/cf-1.<<"!!!\n";
    cf = ch + smallB;
  }

  if (cs > ca) cs = ca - smallB;
  if (cs <= 0. || cs > ca) cs = ca / 2.;
  if (ca > cf) cf = ca + smallB;

  /** \subsection alphafs Roe+Balsara Normalisation constants
   * Calculate the eigenvector normalisation constants, as follows:
   * \f[ \alpha_f = \frac{a^2-c^2_s}{c^2_f-c^2_s} \;,\qquad
   *     \alpha_s = \frac{c^2_f-a^2}{c^2_f-c^2_s} \f]
   * If I am clever, I have already ensured that the sound speeds are all
   * finite, and as a byproduct, they are all greater than zero.  So I
   * shouldn't have any trouble with machine accuracy here.  I check anyway.\n
   * Where I might have trouble is that I have set cf=ch+smallB, cs=ch-smallB,
   * so this enforces alphaf^2=alphax^2=1/2.  But this is only enforced
   * sometimes.
   * */
  double cf2diff;
  if ((cf2diff = cf * cf - cs * cs) > smallB) {
    if ((alphaf = ch * ch - cs * cs) <= smallB) alphaf = 0.;
    if ((alphas = cf * cf - ch * ch) <= smallB) alphas = 0.;
    if ((alphaf = sqrt(alphaf / cf2diff)) > 1.) {
      // cout<<"alpha_f = "<<alphaf<<" !!!\n";
      alphaf = 1.;
    }
    if ((alphas = sqrt(alphas / cf2diff)) > 1.) {
      // cout<<"alpha_s>1!!!\n";
      alphas = 1.;
    }
  }
  else {
    //
    // We are near the triple degeneracy point, and haven't fixed things
    //  already (hopefully won't happen)..
    //
    spdlog::warn("Near Triple degeneracy point! (and didn't realise it!)"
                 " THINGS ARE PROBABLY GOING BAD...");
    alphaf = alphas = 1. / sqrt(2.);
    spdlog::error("left  state : {}", RS_left);
    spdlog::error("right state : {}", RS_right);
    spdlog::error("Avg.  state : {}", RS_meanp);
    spdlog::error("{}: {}", "Bugging out for now...", 99);
    exit(8);
  }

  //
  // Check sound speeds are positive definite!
  //
  if ((cf <= 0.) || (cs < 0.) || (ca < 0.) || (ch <= 0.)) {
    spdlog::error(
        "(riemann_MHD::get_sound_speeds) Error... sound speeds are negative.\n");
    spdlog::debug(
        "Speeds: (ch2,ca2,ct2,cs2,cf2): ({}, {}, {}, {}, {})", ch * ch, ca * ca,
        bt * bt, cs * cs, cf * cf);
    return (1);
  }

  //
  // PHEW!!! If we get here then we got sensible wave-speeds.  As you
  // can tell, this is the function where things start to go wrong
  // most of the time, hence the extra work to check everything is
  // sane.
  //
  return (0);
}

// ##################################################################
// ##################################################################

void riemann_MHD::get_eigenvalues()
{
  // This constructs the eigenvalues of the average matrix
  // \bar{A}. The 'average' sound speeds are already calculated.
  RS_evalue[FN] = RS_meanp[RVX] - cf;
  RS_evalue[FP] = RS_meanp[RVX] + cf;
  RS_evalue[AN] = RS_meanp[RVX] - ca;
  RS_evalue[AP] = RS_meanp[RVX] + ca;
  RS_evalue[SN] = RS_meanp[RVX] - cs;
  RS_evalue[SP] = RS_meanp[RVX] + cs;
  RS_evalue[CT] = RS_meanp[RVX];
  //  for (waves i=FN; i<=FP; ++i) {
  //    cout  << "(riemann_MHD::get_eigenvalues) evalue[" << i << "]: ";
  //    cout  << RS_evalue[i] << "\n";
  //  }

  return;
}

// ##################################################################
// ##################################################################

inline double riemann_MHD::dot_product(pion_flt *v1, pion_flt *v2, int nd)
{
  double temp = 0.0;
  for (int i = 0; i < nd; i++) {
    temp += v1[i] * v2[i];
  }
  return (temp);
}

// ##################################################################
// ##################################################################

void riemann_MHD::calculate_wave_strengths()
{
  riemann_MHD::getPdiff();
  for (waves i = FN; i <= FP; ++i) {
    RS_strength[i] =
        dot_product(RS_leftevec[i].data(), RS_pdiff.data(), RS_nvar);
    // cout  << "wavestrength[" << i << "]: " << RS_strength[i] << "\n";
  }
  //  rep.printVec("strenth",RS_strength,RS_nvar);

  return;
}

// ##################################################################
// ##################################################################

void riemann_MHD::getPdiff()
{
  RS_pdiff[RRO] = RS_right[RRO] - RS_left[RRO];
  RS_pdiff[RPG] = RS_right[RPG] - RS_left[RPG];
  RS_pdiff[RVX] = RS_right[RVX] - RS_left[RVX];
  RS_pdiff[RVY] = RS_right[RVY] - RS_left[RVY];
  RS_pdiff[RVZ] = RS_right[RVZ] - RS_left[RVZ];
  RS_pdiff[RBY] = RS_right[RBY] - RS_left[RBY];
  RS_pdiff[RBZ] = RS_right[RBZ] - RS_left[RBZ];
  return;
}

// ##################################################################
// ##################################################################

int riemann_MHD::get_pstar()
{
  /** \section arrays Arrays
   * pdiff[] is not used for the rest of the calculation,
   * and I need an extra array in this function, so I use pdiff[] to
   * store the result obtained by starting at the right state and
   * crossing waves to the left until I reach x=0.  pstar[] holds the
   * result of starting at the left and crossing waves to the right
   * until I reach x=0.
   * */

  //
  // Now go across waves from the left, assign result to pstar[]
  //
  int i = 0;
  // If the eigenvalue is this close to zero, do averaging.
  double evalacc = 1.e-4;

  for (int j = 0; j < RS_nvar; j++)
    RS_pstar[j] = RS_left[j];

  while ((i < RS_nvar) && (RS_evalue[i] < 0.)) {
    //    cout << "P* from left:  wave["<<i<<"]"<< "\n";
    for (int j = 0; j < RS_nvar; j++) {
      RS_pstar[j] += RS_strength[i] * RS_rightevec[i][j];
    }
    i++;
  }

  /** \section contact Contact Discontinuity
   * If pave[RVX]=0, then I need to do something to decide which value of
   * density to choose.  the most sensible thing is to look at pstar[RVX],
   * and decide on the basis of that.
   * */
  //  cout <<"pave rvx = "<<pave[RVX]<<"\n";

  //
  // if Pstar[Vx] is small, then set Pstar to be the mean of
  // left/right starred states.
  //
  if (fabs(RS_meanp[RVX]) < (evalacc * ch)) {
    // cout <<"contact discontinuity is stationary. ch="<<ch;
    // cout <<"\t pstar[RVX]="<<RS_pstar[RVX]<<"\n";

    //
    // Now go across waves from the right, assign result to pdiff[]
    //
    i = RS_nvar - 1;
    for (int j = 0; j < RS_nvar; j++)
      RS_pdiff[j] = RS_right[j];
    //  while ((RS_evalue[i]>(evalacc*ch)) && (i>=0)) {
    while ((i >= 0) && (RS_evalue[i] > 0.)) {
      // cout << "P* from right: wave["<<i<<"]"<< "\n";
      for (int j = 0; j < RS_nvar; j++) {
        RS_pdiff[j] -= RS_strength[i] * RS_rightevec[i][j];
      }
      i--;
    }

    //
    // Now pstar is the mean of the values obtained by going from left
    // to right, and then right to left.
    //
    for (int v = 0; v < RS_nvar; v++) {
      RS_pstar[v] = 0.5 * (RS_pstar[v] + RS_pdiff[v]);
      RS_pdiff[v] = RS_pstar[v];
    }
  }

#ifdef RS_TESTING
  for (i = 0; i < RS_nvar; i++) {
    if (fabs(RS_pstar[i] - RS_pdiff[i])
            / (fabs(RS_pstar[i]) + fabs(RS_pdiff[i]) + SMALLVALUE)
        > 5.e-2) {
      spdlog::debug(
          "(riemann_MHD::get_pstar) left-ans[{0}] - right-ans[{0}] = {}\nEvalues[] = [",
          i, i, RS_pstar[i] - RS_pdiff[i]);
      for (int j = 0; j < (RS_nvar - 1); j++)
        spdlog::debug("{}, ", RS_evalue[j]);
      spdlog::debug("{}]\nleft[]    = [", RS_evalue[RS_nvar - 1]);
      for (int j = 0; j < (RS_nvar - 1); j++)
        spdlog::debug("{}, ", RS_left[j]);
      spdlog::debug("{}]\nright[]   = [", RS_left[RS_nvar - 1]);
      for (int j = 0; j < (RS_nvar - 1); j++)
        spdlog::debug("{}, ", RS_right[j]);
      spdlog::debug("{}]\npstar_l[] = [", RS_right[RS_nvar - 1]);
      for (int j = 0; j < (RS_nvar - 1); j++)
        spdlog::debug("{}, ", RS_pstar[j]);
      spdlog::debug("{}]\npstar_r[] = [", RS_pstar[RS_nvar - 1]);
      for (int j = 0; j < (RS_nvar - 1); j++)
        spdlog::debug("{}, ", RS_pdiff[j]);
      spdlog::debug(
          "{}]\nBx = {}\t alphaf,s = {} {}", RS_pdiff[RS_nvar - 1], ansBX,
          alphaf, alphas);
      return (1);
    }
  }
#endif

  return (0);
}

// ##################################################################
// ##################################################################

int riemann_MHD::RoeBalsara_evectors()
{
  //
  // This calculates the Roe and Balsara (1996) eigenvectors.  They
  // are more numerically robust than the evectors proposed by
  // Falle, Komissarov and Joarder (1998) for the fast and slow waves.
  // I still use FKJ98's Alfven wave normalisation.
  //

  double r2 = sqrt(2.);
  int sBx;  // Sign of B_x
  if (ansBX < 0.)
    sBx = -1;
  else if (ansBX > 0.)
    sBx = 1;
  else {
    // cout << "Bx = 0, assigning sign(Bx) = +1... make sure this is ok!!!"
    // <<
    // "\n";
    sBx = 1;
  }
  // Have got alpha[f,s], beta[y,z] from get_sound_speeds().
  // Eigenvectors: first index refers to which wave it corresponds to,
  // and the second to which element of the vector it is.
  // ****LEFT EIGENVECTORS****
  // Left Fast Magnetosonic Wave
  RS_leftevec[FN][RRO] = 0.0;
  RS_leftevec[FN][RVX] = -alphaf * cf;
  RS_leftevec[FN][RVY] = alphas * cs * sBx * betay;
  RS_leftevec[FN][RVZ] = alphas * cs * sBx * betaz;
  RS_leftevec[FN][RPG] = alphaf / RS_meanp[RRO];
  RS_leftevec[FN][RBY] = alphas * ch * betay / sqrt(RS_meanp[RRO]);
  RS_leftevec[FN][RBZ] = alphas * ch * betaz / sqrt(RS_meanp[RRO]);
  // RS_Left Alfven Wave
  RS_leftevec[AN][RRO] = 0.;
  RS_leftevec[AN][RVX] = 0.;
  RS_leftevec[AN][RVY] = sBx * betaz / r2;
  RS_leftevec[AN][RVZ] = -sBx * betay / r2;
  RS_leftevec[AN][RPG] = 0.;
  RS_leftevec[AN][RBY] = betaz / sqrt(RS_meanp[RRO]) / r2;
  RS_leftevec[AN][RBZ] = -betay / sqrt(RS_meanp[RRO]) / r2;
  // Left Slow Magnetosonic Wave
  RS_leftevec[SN][RRO] = 0.0;
  RS_leftevec[SN][RVX] = -alphas * cs;
  RS_leftevec[SN][RVY] = -alphaf * cf * sBx * betay;
  RS_leftevec[SN][RVZ] = -alphaf * cf * sBx * betaz;
  RS_leftevec[SN][RPG] = alphas / RS_meanp[RRO];
  RS_leftevec[SN][RBY] = -alphaf * ch * betay / sqrt(RS_meanp[RRO]);
  RS_leftevec[SN][RBZ] = -alphaf * ch * betaz / sqrt(RS_meanp[RRO]);
  // Contact Discontinuity
  RS_leftevec[CT][RRO] = 1.;
  RS_leftevec[CT][RVX] = 0.;
  RS_leftevec[CT][RVY] = 0.;
  RS_leftevec[CT][RVZ] = 0.;
  RS_leftevec[CT][RPG] = -1 / ch / ch;
  RS_leftevec[CT][RBY] = 0.;
  RS_leftevec[CT][RBZ] = 0.;
  // Right Slow Magnetosonic Wave
  RS_leftevec[SP][RRO] = 0.0;
  RS_leftevec[SP][RVX] = -RS_leftevec[SN][RVX];
  RS_leftevec[SP][RVY] = -RS_leftevec[SN][RVY];
  RS_leftevec[SP][RVZ] = -RS_leftevec[SN][RVZ];
  RS_leftevec[SP][RPG] = RS_leftevec[SN][RPG];
  RS_leftevec[SP][RBY] = RS_leftevec[SN][RBY];
  RS_leftevec[SP][RBZ] = RS_leftevec[SN][RBZ];
  // Right Alfven Wave
  RS_leftevec[AP][RRO] = 0.;
  RS_leftevec[AP][RVX] = 0.;
  RS_leftevec[AP][RVY] = RS_leftevec[AN][RVY];
  RS_leftevec[AP][RVZ] = RS_leftevec[AN][RVZ];
  RS_leftevec[AP][RPG] = 0.;
  RS_leftevec[AP][RBY] = -RS_leftevec[AN][RBY];
  RS_leftevec[AP][RBZ] = -RS_leftevec[AN][RBZ];
  // Right Fast Magnetosonic wave
  RS_leftevec[FP][RRO] = 0.0;
  RS_leftevec[FP][RVX] = -RS_leftevec[FN][RVX];
  RS_leftevec[FP][RVY] = -RS_leftevec[FN][RVY];
  RS_leftevec[FP][RVZ] = -RS_leftevec[FN][RVZ];
  RS_leftevec[FP][RPG] = RS_leftevec[FN][RPG];
  RS_leftevec[FP][RBY] = RS_leftevec[FN][RBY];
  RS_leftevec[FP][RBZ] = RS_leftevec[FN][RBZ];
  // ****RIGHT EIGENVECTORS****
  // RS_Left Fast Magnetosonic Wave
  RS_rightevec[FN][RRO] = alphaf * RS_meanp[RRO];
  RS_rightevec[FN][RVX] = RS_leftevec[FN][RVX];
  RS_rightevec[FN][RVY] = RS_leftevec[FN][RVY];
  RS_rightevec[FN][RVZ] = RS_leftevec[FN][RVZ];
  RS_rightevec[FN][RPG] = alphaf * RS_meanp[RRO] * ch * ch;
  RS_rightevec[FN][RBY] = RS_leftevec[FN][RBY] * RS_meanp[RRO];
  RS_rightevec[FN][RBZ] = RS_leftevec[FN][RBZ] * RS_meanp[RRO];
  // RS_Left Alfven Wave
  RS_rightevec[AN][RRO] = 0.;
  RS_rightevec[AN][RVX] = 0.;
  RS_rightevec[AN][RVY] = RS_leftevec[AN][RVY];
  RS_rightevec[AN][RVZ] = RS_leftevec[AN][RVZ];
  RS_rightevec[AN][RPG] = 0.;
  RS_rightevec[AN][RBY] = RS_leftevec[AN][RBY] * RS_meanp[RRO];
  RS_rightevec[AN][RBZ] = RS_leftevec[AN][RBZ] * RS_meanp[RRO];
  // Left Slow Magnetosonic Wave
  RS_rightevec[SN][RRO] = alphas * RS_meanp[RRO];
  RS_rightevec[SN][RVX] = RS_leftevec[SN][RVX];
  RS_rightevec[SN][RVY] = RS_leftevec[SN][RVY];
  RS_rightevec[SN][RVZ] = RS_leftevec[SN][RVZ];
  RS_rightevec[SN][RPG] = alphas * RS_meanp[RRO] * ch * ch;
  RS_rightevec[SN][RBY] = RS_leftevec[SN][RBY] * RS_meanp[RRO];
  RS_rightevec[SN][RBZ] = RS_leftevec[SN][RBZ] * RS_meanp[RRO];
  // Contact Discontinuity
  RS_rightevec[CT][RRO] = 1.0;
  RS_rightevec[CT][RVX] = 0.;
  RS_rightevec[CT][RVY] = 0.;
  RS_rightevec[CT][RVZ] = 0.;
  RS_rightevec[CT][RPG] = 0.;
  RS_rightevec[CT][RBY] = 0.;
  RS_rightevec[CT][RBZ] = 0.;
  // Right Slow Magnetosonic Wave
  RS_rightevec[SP][RRO] = RS_rightevec[SN][RRO];
  RS_rightevec[SP][RVX] = -RS_rightevec[SN][RVX];
  RS_rightevec[SP][RVY] = -RS_rightevec[SN][RVY];
  RS_rightevec[SP][RVZ] = -RS_rightevec[SN][RVZ];
  RS_rightevec[SP][RPG] = RS_rightevec[SN][RPG];
  RS_rightevec[SP][RBY] = RS_rightevec[SN][RBY];
  RS_rightevec[SP][RBZ] = RS_rightevec[SN][RBZ];
  // Right Alfven Wave
  RS_rightevec[AP][RRO] = 0.;
  RS_rightevec[AP][RVX] = 0.;
  RS_rightevec[AP][RVY] = RS_rightevec[AN][RVY];
  RS_rightevec[AP][RVZ] = RS_rightevec[AN][RVZ];
  RS_rightevec[AP][RPG] = 0.;
  RS_rightevec[AP][RBY] = -RS_rightevec[AN][RBY];
  RS_rightevec[AP][RBZ] = -RS_rightevec[AN][RBZ];
  // Right Fast Magnetosonic wave
  RS_rightevec[FP][RRO] = RS_rightevec[FN][RRO];
  RS_rightevec[FP][RVX] = -RS_rightevec[FN][RVX];
  RS_rightevec[FP][RVY] = -RS_rightevec[FN][RVY];
  RS_rightevec[FP][RVZ] = -RS_rightevec[FN][RVZ];
  RS_rightevec[FP][RPG] = RS_rightevec[FN][RPG];
  RS_rightevec[FP][RBY] = RS_rightevec[FN][RBY];
  RS_rightevec[FP][RBZ] = RS_rightevec[FN][RBZ];

  // Now multiply left fast and slow evectors by 1/(2a^2) to make them
  // normalised.
  double a22 = 1. / (2. * ch * ch);
  for (int i = 0; i < RS_nvar; i++) {
    RS_leftevec[FN][i] *= a22;
    RS_leftevec[SN][i] *= a22;
    RS_leftevec[SP][i] *= a22;
    RS_leftevec[FP][i] *= a22;
  }
  return (0);
}

// ##################################################################
// ##################################################################

int riemann_MHD::check_evectors()
{
  // Printing the eigenvectors to screen!
  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug(" leftevec[{}] = [", static_cast<int>(i));
    for (int j = 0; j < RS_nvar - 1; j++) {
      spdlog::debug("{}, ", RS_leftevec[i][j]);
    }
    spdlog::debug("{}]", RS_leftevec[i][RS_nvar - 1]);
  }
  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug("rightevec[] = ", RS_rightevec[i]);
    // spdlog::debug("rightevec[{}] = [", i);
    // for (int j = 0; j < RS_nvar - 1; j++) {
    //  spdlog::debug("{}, ", RS_rightevec[i][j]);
    //}
    // spdlog::debug("{}]", RS_rightevec[i][RS_nvar - 1]);
  }

  // Print the dot products
  // On-axis dot products
  spdlog::debug("left[i].right[i] = [");
  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug(
        "{}, ",
        dot_product(RS_leftevec[i].data(), RS_rightevec[i].data(), RS_nvar));
  }
  spdlog::debug("]");
  // off-axis dot products.
  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug("left[{}].right[j] = ", static_cast<int>(i));
    for (waves j = FN; j <= FP; ++j) {
      spdlog::debug(
          "\t{}",
          dot_product(RS_leftevec[i].data(), RS_rightevec[j].data(), RS_nvar));
    }
    spdlog::debug("]");
  }

  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug("left[{}].left[j] = [ ", static_cast<int>(i));
    for (waves j = FN; j <= FP; ++j) {
      spdlog::debug(
          "\t{}",
          dot_product(RS_leftevec[i].data(), RS_leftevec[j].data(), RS_nvar));
    }
    spdlog::debug("]");
  }

  for (waves i = FN; i <= FP; ++i) {
    spdlog::debug("right[{}].right[j] = [ ", static_cast<int>(i));
    for (waves j = FN; j <= FP; ++j) {
      spdlog::debug(
          "\t{}",
          dot_product(RS_rightevec[i].data(), RS_rightevec[j].data(), RS_nvar));
    }
    spdlog::debug("]");
  }

  //
  // Check the dot products are what they should be
  // First the i=j products; should be 1.
  // Setting the tolerance: 1.e-12 is roughly truncation error.
  //
  double test;
  double errtol = 1.e-8;
  for (waves i = FN; i <= FP; ++i) {
    if ((test = fabs(
             dot_product(RS_leftevec[i].data(), RS_rightevec[i].data(), RS_nvar)
             - 1.))
        > errtol) {
      spdlog::debug(
          "evectors not properly normalised: (l[{0}].r[{0}] -1) = {1}",
          static_cast<int>(i), test);
      onaxis++;
      return (1);
    }
  }
  for (waves i = FN; i <= FP; ++i) {
    for (waves j = FN; j <= FP; ++j) {
      if (j != i) {
        if ((test = fabs(dot_product(
                 RS_leftevec[i].data(), RS_rightevec[j].data(), RS_nvar)))
            >= errtol) {
          spdlog::debug(
              "evectors not properly normalised: left[{}].right[{}] = {}",
              static_cast<int>(i), static_cast<int>(i), test);
          offaxis++;
          return (1);
        }
      }
    }
  }
  return (0);
}

// ##################################################################
// ##################################################################
