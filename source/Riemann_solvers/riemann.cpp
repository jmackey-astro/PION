///
/// \file riemann.cc
/// \brief Class definition for Hydrodynamics Riemann Solver
/// 
///  \author Jonathan Mackey
/// 
/// This is a Hydrodynamics Riemann Solver for the Euler Equations.
/// It works by first calling a linear solver, and if the results
/// are within certain tolerances it returns the answer, but if not
/// it calls an exact solver and returns the result from that.
/// 
/// References:
/// -- Hirsch, C.,  p.206 for the shock equations.
///    "Numerical Computation of Internal and External Flows: Volume 2" (1990, Wiley)
/// -- Toro, E.F.,  for structure of approximate/linear riemann solvers.
///    "Riemann Solvers and Numerical Methods for Fluid Dynamics" (1999, Springer)
/// 
/// Written 2006-11-09-Thursday
/// Modified :
/// - 2006-12-20 Testing for linear solver added.
/// - 2007-01-12 Lots more testing added; linearOK() changed.
/// - 2007-07-10 Added Index variables as class members, and a set-direction function.
/// - 2009-10-20 Changed class heirarchy, changed functions, and some of the interface.
/// - 2009-10-23 made Riemann_Euler inherit from findroot, and
///    redefine the function to get the root of!
/// - 2010-07-27 JM: cleaned up a bit, updated comments.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.23 JM: Got rid of last rsvar enum references.  Got rid
///   of RS_* arrays and made new private data rs_*[].  This will
///   avoid confusion in derived classes.  Inherited data arrays is 
///   a bad idea.
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer 
///    fluxes are dealt with by the flux-solver classes.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.14 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"

#include "Riemann_solvers/riemann.h"
#include <iostream>
using namespace std;


// ##################################################################
// ##################################################################


// This version is specific to solving the equation in the Exact Riemann Solver.
// Other versions could be specified, taking more or fewer parameters.
int riemann_Euler::FR_find_root(
      pion_flt *ans, ///< pointer to result
      const pion_flt p1,  ///< p_g left
      const pion_flt p2,  ///< p_g right
      const pion_flt p3,  ///< eq_gamma
      const pion_flt p4,  ///< cs_left
      const pion_flt p5   ///< cs_right
      )
{

  FR_param1 = p1;  ///< e.g. left state pressure
  FR_param2 = p2;  ///< e.g. right state pressure
  FR_param3 = p3;  ///< e.g. EOS gamma
  FR_param4 = p4;  ///< e.g. left state sound speed
  FR_param5 = p5;  ///< e.g. right state sound speed

  //
  // Set initial values for x1, x2
  // My initial guesses are 1/3 of and 3 times the arithmetic mean 
  // of the left and right pressures.
  //
  pion_flt x1 = (FR_param1+ FR_param2)/6.0;
  pion_flt x2 = x1*9.0;

  // 
  // Call the common solver, now that parameters are set properly, to
  // find the root in [x1,x2] given that the root is positive, and put
  // the solution in 'ans'
  //
  int err = findroot::solve_pos(x1, x2, ans);
  if (err!=0) {
    cerr << "(riemann_Euler::findroot::solve_riemann) solve_pos exited abnormally" << "\n";
    return(1);
  }
  //  cout << "(riemann_Euler::findroot::solve_riemann) Success: ans = " << *ans << "\n";
  return(0);
}


// ##################################################################
// ##################################################################


pion_flt riemann_Euler::FR_root_function(
      pion_flt pp
      )
{
  //
  // Wrapper function to solve across both waves given p*=pp, which returns
  // the velocity difference (u*(R) - u*(L)) = *uu
  // This takes in the test value of p*, and calculates the resulting
  // values of u* for the left and right waves, and returns the 
  // difference.  We want to find where the difference is zero, so 
  // basically we want to find the root of this 'equation'.
  // 
  int err=0;
  pion_flt ustarL, ustarR;
  err = HydroWave(XN, pp,  rs_left, &ustarL, eq_gamma);
  err+= HydroWave(XP, pp, rs_right, &ustarR, eq_gamma);
  if (err!=0) {
    cerr << "(riemann_Euler::FR_root_function) hydro_wave calls returned abnormally. " << "\n";
    return(1.e99);
  }
  return(ustarR - ustarL);
}


// ##################################################################
// ##################################################################


//
// Constructor of the riemann class: set up variables, allocate memory, etc.
//
riemann_Euler::riemann_Euler(
      const int nv, ///< number of vars.
      const pion_flt *state, ///< reference vec.
      const double g ///< gas EOS.
      )
  : eqns_base(nv),
    eqns_Euler(nv),
    findroot(), rs_nvar(5)
{
#ifdef FUNCTION_ID
  cout <<"riemann_Euler::riemann_Euler ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  cout <<"riemann_Euler::riemann_Euler: eqnvar="<<eq_nvar<<"\n";
#endif
  if(eq_nvar<5) {
    rep.error("Problem with HD Riemann Solver... eq_nvar<5. Please use [rho,p_g,vx,vy,vz]. Quitting!!!",eq_nvar);
  }

  eq_gamma = g;    // gamma is the gas constant, e.g. 5/3 or 7/5
  //cout <<"!!!!!!!!gamma="<<eq_gamma<<"\n";

#ifdef RSTESTING
  linearpstar = mem.myalloc(linearpstar, rs_nvar);
#endif //RSTESTING  // **************

  //
  // set counters to zero
  //
  linct = exact = total = 0;
#ifdef RSTESTING
  fails = failplr = failpst = failust = failrare = failcomp = faildens = 0;
#endif //RSTESTING  // **************
  
  SetDirection(XX); // Initially assume we are looking along the x-axis.
  
  //
  // Set the reference state vector (so we know what order of
  // magnitude data to expect).
  //
  SetAvgState(state,eq_gamma);
  
  samestate=0;

  //
  // Allocate memory for arrays.
  //
  rs_left  = mem.myalloc(rs_left ,rs_nvar);
  rs_right = mem.myalloc(rs_right,rs_nvar);
  rs_meanp = mem.myalloc(rs_meanp,rs_nvar);
  rs_pstar = mem.myalloc(rs_pstar,rs_nvar);

#ifdef FUNCTION_ID
  cout <<"riemann_Euler::riemann_Euler ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


//
// Destructor: delete memory etc.
//
riemann_Euler::~riemann_Euler()
{
#ifdef FUNCTION_ID
  cout <<"riemann_Euler::~riemann_Euler ...starting.\n";
#endif //FUNCTION_ID

#ifdef RSTESTING
  riemann_Euler::testing();
#endif //RSTESTING
  rs_left   = mem.myfree(rs_left);
  rs_right  = mem.myfree(rs_right);
  rs_meanp  = mem.myfree(rs_meanp);
  rs_pstar  = mem.myfree(rs_pstar);
#ifdef RSTESTING
  linearpstar = mem.myfree(linearpstar);
#endif //RSTESTING

#ifdef FUNCTION_ID
  cout <<"riemann_Euler::~riemann_Euler ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################


//
// Output diagnostic info -- how many of each type of solution we did.
//
void riemann_Euler::testing()
{
  cout <<"(riemann_Euler::testing) No. Solves: linear="<<linct<<" exact="<<exact<<" total="<<total<<"."<<"\n";
  cout <<"\tFraction linear/total: "<< static_cast<double>(linct)/(linct+exact) <<"\n";
  cout <<"Number of solves with same-state = "<<samestate<< " out of total numbe of solves above.\n";
#ifdef RSTESTING
  cout <<"\tNumber of fails: "<<fails<<"  out of "<<linct<<" linear solves corresponding to ";
  cout <<static_cast<double>(fails*100)/linct<< "%."<<"\n";
  cout<<"\tWhy not linear solver: pratio:"<<failplr<<" pstar:"<<failpst<<" ustar<<1:"<<failust;
  cout <<" rarefaction:"<< failrare <<" compression:"<<failcomp <<" dratio:"<<faildens<<"\n";
#endif //RSTESTING
}


// ##################################################################
// ##################################################################



//
// Public interface to the riemann solver -- it takes in a left and right state, 
// and then sets the solver going, which returns the state P* in 'result'.
// err_code  solve(left,     right,    result  );
//
int riemann_Euler::JMs_riemann_solve(
      const pion_flt *l,
      const pion_flt *r,
      pion_flt *ans,
      const int mode,
      ///< Solve Type (1=LinearRS,2=ExactRS,3=HybridRS)
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"riemann_Euler::JMs_riemann_solve ...starting.\n";
#endif //FUNCTION_ID

  total++;
  int err = 0;
#ifdef RSTESTING
  for (int v=0;v<rs_nvar;v++) {
    if (isnan(l[v]) || isnan(r[v])) {
      cout <<"NAN's detected!!!\n";
      rep.printVec("left ",l,rs_nvar);
      rep.printVec("right",r,rs_nvar);
      return 99;
    }
  }
#endif

  /// \section Gamma
  /// I assign gamma to a class variable every time.  This is so that all
  /// functions in the class have access to it without having to pass it 
  /// around everywhere.  It will never change in the solver, unless the 
  /// left and right states have different gammas.  That's a big upgrade in my
  /// degree of accuracy though, because then p* will also have its own
  /// gamma, and the solver will need a rewrite if I get to that stage.
  /// 
  eq_gamma = g;
  //cout << "(riemann_Euler::solve) Assigning data, gamma: " << riemann_Euler::eq_gamma << "\n";

  //
  // Assign class variables from values passed 
  //
  if (l[eqRO]<TINYVALUE || l[eqPG]<TINYVALUE ||
      r[eqRO]<TINYVALUE || r[eqPG]<TINYVALUE) {
    rep.printVec("left ",l,rs_nvar);
    rep.printVec("right",r,rs_nvar);
    rep.error("riemann::solve() Density/Pressure too small",l[eqRO]);
  }
  for (int v=0; v<rs_nvar; v++)
    rs_left[v] = l[v];
  for (int v=0; v<rs_nvar; v++)
    rs_right[v] = r[v];

  //  
  // First test if the left and right states are the same, and if they
  // are then return the simple answer.
  //
  double diff=0.;
  for(int i=0;i<rs_nvar; i++)
    diff += fabs(rs_right[i]-rs_left[i])/(fabs(eq_refvec[i])+TINYVALUE);
  if (diff <1.e-6) {
    samestate++;
    for(int i=0;i<rs_nvar; i++) {
      rs_pstar[i] = (l[i]+r[i])/2.;
      ans[i] = rs_pstar[i];
    }
    return(0);
  }

  //
  // Sound speeds.
  //
  cl = chydro(rs_left,  eq_gamma);
  cr = chydro(rs_right, eq_gamma);
  //  cout << "(riemann_Euler::solve) Assigned data, here it is: " << eq_gamma << "\n";

  //
  // Analytically we can show that:
  // (1) We get cavitation if     (u_R-u_L) >= 2(c_R-c_L)/(g-1)
  // (2) We get 2 rarefactions if (u_R-u_L) >= 2/(g-1)*[c_L+ sqrt((g-1)/2g)c_R]
  //
  // If two rarefactions or cavitation are definitely present, then I
  // can write down the solution.  Since these are likely to be rare,
  // we first check to see if we need to do the standard solution.
  //
  // First we test to see if the 2 rarefactions condition *doesn't*
  // hold, in which case we have to do a normal solve (either
  // approximate or iterative):
  // 
  if( (rs_right[eqVX]-rs_left[eqVX])<=
      2.*(cl+sqrt((eq_gamma-1.)/2./eq_gamma)*cr)/(eq_gamma-1.) ) {
    //
    //Do the normal solver
    //
    // *******************************************************************
    // Now we do the main part of the solution.
    // ******************************************************************
    switch (mode) {
    case (FLUX_RSlinear): // ALL LINEAR SOLVES.
      err = linear_solver();
      if (err != 0) {
	cerr << "(riemann_Euler::solve) linear_solver() returned abnormally!" << "\n";
	rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar[eqVX] = TINYVALUE;
	for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
	return(1);
      }
      linct++;
      break;
    case (FLUX_RSexact): // ALL EXACT SOLVES
      err = exact_solver();
      if (err != 0) {
	cerr << "(riemann_Euler::solve) exact_solver() returned abnormally - negative pressure!" << "\n";
	rs_pstar[eqPG] = rs_pstar[eqRO] = TINYVALUE;
	for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
	return(1);
      }
      exact++;
      break;
    case (FLUX_RShybrid):
      //
      // LINEAR SOLVES WHERE POSSIBLE, AND EXACT IN THE HARD CASES.
      // First call the linear solver and if it fails or, if the
      // states differ by too much, then call the exact solver.
      //
      err = linear_solver();
      linct++;
      if (err != 0) {
	cerr << "(riemann_Euler::solve) linear_solver() returned abnormally!" << "\n";
	rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar[eqVX] = TINYVALUE;
	//for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
	//return(1);
      }
#ifdef RSTESTING
      // 
      // For testing we do both the linear and exact solver and compare the results.
      // Assign the linear result to a temp variable, to compare with the exact one.
      //
      for(int i=0; i<rs_nvar;i++) {
	linearpstar[i] = rs_pstar[i];
      }
#endif //RSTESTING
      //
      // If linear solver failed, or if the left and right states are
      // too different, then call the exact solver.
      //
      if(err!=0 || linearOK()!=0) {
	err = exact_solver();
	if (err != 0) {
	  cerr << "(riemann_Euler::solve) exact_solver() returned abnormally - negative pressure!" << "\n";
	  rs_pstar[eqPG] = rs_pstar[eqRO] = TINYVALUE;
	  for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
	  return(err);
	}
	exact++;linct--;
      }
#ifdef RSTESTING
      // ***TESTING***
      else {
	err = exact_solver();
	if (err != 0) {
	  cerr << "(riemann_Euler::solve) exact_solver() returned abnormally!" << "\n";
	  rs_pstar[eqPG] = rs_pstar[eqRO] = TINYVALUE;
	  for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
	  return(1);
	}
	// Check the exact pstar compared to the linear one.
	//Verbose output
	//    cout <<"Linear Solver obtained: ";
	// 	   cout << "[rho,v,p] = ["<< linearpstar[eqRO] <<", "<< linearpstar[eqVX] <<", "<< linearpstar[eqPG] << " ]"<< "\n";
	// 	   cout <<"Difference in liner/exact solver: delta";
	// 	   cout << "[rho,v,p] = ["<< (rs_pstar[eqRO]-linearpstar[eqRO])/rs_pstar[eqRO];
	// 	   cout <<", "<< (rs_pstar[eqVX]-linearpstar[eqVX]) <<", ";
	// 	   cout << (rs_pstar[eqPG]-linearpstar[eqPG])/rs_pstar[eqPG] << " ]"<<"\n";
	// 	   cout <<"Differences are fractional error in rho,p, and absolute error in vx."<<"\n";
	lineartol = 0.1;
	if( (fabs(rs_pstar[eqRO]-linearpstar[eqRO])/(rs_pstar[eqRO]+linearpstar[eqRO])>lineartol) ||
	    (fabs(rs_pstar[eqVX]-linearpstar[eqVX])/(eq_refvec[eqVX]+fabs(rs_pstar[eqVX]))>lineartol) ||
	    (fabs(rs_pstar[eqPG]-linearpstar[eqPG])/(rs_pstar[eqPG]+linearpstar[eqPG])>lineartol) ) {
	  rep.printVec("refvec",eq_refvec,rs_nvar);
	  cout <<"!!! "<<fabs(rs_pstar[eqRO]-linearpstar[eqRO])/(rs_pstar[eqRO]+linearpstar[eqRO]);
	  cout <<"    "<<fabs(rs_pstar[eqVX]-linearpstar[eqVX])/(eq_refvec[eqVX]+fabs(rs_pstar[eqVX]));
	  cout <<"    "<<fabs(rs_pstar[eqPG]-linearpstar[eqPG])/(rs_pstar[eqPG]+linearpstar[eqPG])<<"\n";
	  fails++;
	  cout << "Left " << "                    [rho,v,p] =";
	  cout << " ["<< rs_left[eqRO]  <<", "<< rs_left[eqVX]  <<", "<< rs_left[eqPG]  << " ]\n";
	  cout << "Right" << "                    [rho,v,p] =";
	  cout << " ["<< rs_right[eqRO] <<", "<< rs_right[eqVX] <<", "<< rs_right[eqPG] << " ]\n";
	  cout <<"Linear Solver obtained: ";
	  cout << "[rho,v,p] = ["<< linearpstar[eqRO] <<", "<< linearpstar[eqVX] <<", "<< linearpstar[eqPG] << " ]"<< "\n";
	  cout <<" Exact Solver obtained: ";
	  cout << "[rho,v,p] = ["<< rs_pstar[eqRO] <<", "<< rs_pstar[eqVX] <<", "<< rs_pstar[eqPG] << " ]"<< "\n";
	  cout <<"Difference :       delta";
	  cout << "[rho,v,p] = ["<< (rs_pstar[eqRO]-linearpstar[eqRO])/rs_pstar[eqRO];
	  cout <<", "<< (rs_pstar[eqVX]-linearpstar[eqVX]) <<", ";
	  cout << (rs_pstar[eqPG]-linearpstar[eqPG])/rs_pstar[eqPG] << " ]";
	  cout <<"c*: "<< sqrt(eq_gamma*rs_pstar[eqPG]/rs_pstar[eqRO]) << " u*: "<< rs_pstar[eqVX]<<"\n";
	}
	// Set pstar back to linear result, as that is what the code will calculate when not in testing mode.
	rs_pstar[eqRO] = linearpstar[eqRO];
	rs_pstar[eqPG] = linearpstar[eqPG];
	rs_pstar[eqVX] = linearpstar[eqVX];
      }
#endif //RSTESTING
      break;
    default:
      rep.error("Error - switch not known in RiemannEul::solve()",mode);
      break;
    } // End of switch (mode) which determines what sort of solver to use.
  } //End of IF we don't have two rarefactions/cavitation, then do the normal solution.
  

  //
  // This checks to see if we *don't* have cavitation, in which case
  // we must have two rarefactions, since we already determined that
  // the normal solution doesn't apply.
  //
  else if ((rs_right[eqVX]-rs_left[eqVX])<=2.*(cl+cr)/(eq_gamma-1.) ) {
    // Two rarefactions, so write down solution.
    //    cout << "(reimann::solve) Definitely two rarefactions!" <<"\n";
    //      cout << "Left " << " [rho,v,p] = ["<< rs_left[eqRO]  <<", "<< rs_left[eqVX]  <<", "<< rs_left[eqPG]  << " ]\n";
    //      cout << "Right" << " [rho,v,p] = ["<< rs_right[eqRO] <<", "<< rs_right[eqVX] <<", "<< rs_right[eqPG] << " ]\n";
    //      cout << " (ur-ul) - 2(cl+sqrt((g-1)/2g)cr)/(g-1) = ";
    //      cout << (rs_right[eqVX]-rs_left[eqVX]) -2.*(cl+sqrt((eq_gamma-1.)/2./eq_gamma)*cr)/(eq_gamma-1.) << "\n";
    //      cout << " (ur-ul) - 2(cl+cr)/(g-1) = ";
    //      cout << (rs_right[eqVX]-rs_left[eqVX]) -2.*(cl+cr)/(eq_gamma-1.);
    //      cout << " ... cl,cr="<<cl<<", "<<cr << "\n";
    // Solve for two rarefactions (checking of wave locations is done inside solve_rarerare()).
    err = solve_rarerare();
    if (err != 0) {
      cerr << "(riemann_Euler::solve) solve_rarerare() didn't work! " << "\n";
      rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar [eqVX] = -1.9e99;
      for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
      return(err);
    }
  }

  //
  // If neither of the above cases were satisfied, then we must have
  // CAVITATION!  This usually means something has gone wrong, but
  // often the code can recover in the next step.
  //
  else {
    //    cout << "(reimann::solve) Cavitation!"<<"\n";
    //      cout << "Left " << " [rho,v,p] = ["<< rs_left[eqRO]  <<", "<< rs_left[eqVX]  <<", "<< rs_left[eqPG]  << " ]\n";
    //      cout << "Right" << " [rho,v,p] = ["<< rs_right[eqRO] <<", "<< rs_right[eqVX] <<", "<< rs_right[eqPG] << " ]\n";
    //      cout << " (ur-ul) - 2(cl+cr)/(g-1) = ";
    //      cout << (rs_right[eqVX]-rs_left[eqVX]) -2.*(cl+cr)/(eq_gamma-1.);
    //      cout << " ... cl,cr="<<cl<<", "<<cr << "\n";
    // Solve for cavitation solution.
    err = solve_cavitation();
    if (err) {
      cerr << "(riemann_Euler::solve) solve_cavitation() didn't work! " << "\n";
      rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar [eqVX] = -1.9e100;
      for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
      return(err);
    }
    // Can't call check_wave_locations() here because the wave structure is different.
    // Have checked the locations in the solve_cavitation() function itself.
  }

  //
  // Need to assign values to v_y,v_z if present.
  // They only change across the contact discontinuity.
  //
  if (rs_pstar[eqVX]>0) {
    rs_pstar[eqVY] = rs_left[eqVY];		    
    rs_pstar[eqVZ] = rs_left[eqVZ];
  }
  else {
    rs_pstar[eqVY] = rs_right[eqVY];		    
    rs_pstar[eqVZ] = rs_right[eqVZ];
  }
  //  cout << "(riemann_Euler::solve) Success!" << "\n";

  //
  //  if(rs_pstar[eqPG]<0.) cerr<<"error in riemann solver! pg="<<rs_pstar[eqPG]<<"\n";
  // TINYVALUE is 1.0e-50
  //
  if (rs_pstar[eqPG]<=TINYVALUE) { 
    rs_pstar[eqPG] = BASEPG*eq_refvec[eqPG];
  }
  if (rs_pstar[eqRO]<=TINYVALUE) {
    rs_pstar[eqRO] = BASEPG*eq_refvec[eqRO];
  }
  
  //
  // Now copy pstar into the ans[] pointer.
  //
  for(int i=0;i<rs_nvar; i++) ans[i] = rs_pstar[i];
  //
  // [DEBUGGING] Check that we haven't altered rs_left or rs_right!
  //
  // for (int v=0; v<rs_nvar; v++) {
  //   if (!pconst.equalD(l[v],rs_left[v] )) rep.error("changed  left state vector!",v);
  //   if (!pconst.equalD(r[v],rs_right[v])) rep.error("changed right state vector!",v);
  // }

#ifdef FUNCTION_ID
  cout <<"riemann_Euler::JMs_riemann_solve ...returning.\n";
#endif //FUNCTION_ID
  return 0;
}


// ##################################################################
// ##################################################################


int riemann_Euler::check_wave_locations()
{
  //
  // Checks that the starred state is valid at the origin.
  // For rarefactions, the wave can straddle the origin, so we need to check that it doesn't!!!
  // It could also be swept downstream by a supersonic flow, so we need to check for that too.
  // For each wave, we first check if it's a rarefaction; 
  // if true then if it could be swept past the origin by supersonic flow (if true assign/return);
  // if not, then if the origin could be within the rarefaction fan (if true assign/return);
  // if not, then move on to the right wave and start again.
  //

  //
  // First the left wave.
  //
  if (rs_pstar[eqPG] < rs_left[eqPG]) { // Left wave is rarefaction.
    if (rs_left[eqVX] >= cl) {
      // In this case everything is swept to the right and the result is the left statE
      rs_pstar[eqPG] = rs_left[eqPG];
      rs_pstar[eqRO] = rs_left[eqRO];
      rs_pstar[eqVX] = rs_left[eqVX];
      //cout << "Shock-Rarefaction... because (v_l > c_l), the resolved state IS the left state. cl= " 
      //     << cl << "\n";
      return(0);
    }
    else if (rs_pstar[eqVX] > 0.){
      //
      // In this case we need to check if the rarefaction straddles x=0 or not.
      //
      double cstar = chydro(rs_pstar, eq_gamma);
      if (rs_pstar[eqVX] > cstar) {
	//
	// in this case the rarefaction straddles the origin.
	//
	//	cout << "Xi_B >0 :x=0 inside left rarefaction: u*-c* = " << rs_pstar[eqVX]-cstar << "\n";
	// Have to find where x=0, and solve for p,u,rho there.
	// Equation 15,16,17 from my rarefaction.ps notes, with \xi=0
	rs_pstar[eqVX] = (2.*cl + rs_left[eqVX]*(eq_gamma -1.))/(eq_gamma +1.);
	rs_pstar[eqRO] = rs_left[eqRO] *exp(2./(eq_gamma-1.)*log(rs_pstar[eqVX]/cl));
	rs_pstar[eqPG] = exp(eq_gamma *log(rs_pstar[eqRO]/rs_left[eqRO]))*rs_left[eqPG];
	//	cout << "Left rarefaction straddles origin." << "\n";
	return(0);
      }
    }
  }
  //
  // Now the right wave.
  //
  if (rs_pstar[eqPG] < rs_right[eqPG]) { // Right wave is rarefaction.
    if (rs_right[eqVX] <= -cr) {
      //
      // In this case everything is swept to the left and the result is the right state
      //
      rs_pstar[eqPG] = rs_right[eqPG];
      rs_pstar[eqRO] = rs_right[eqRO];
      rs_pstar[eqVX] = rs_right[eqVX];
      //cout << "Shock-Rarefaction... because (v_r < -c_r), the resolved state IS the right state. cr= " 
      //	   << cr << "\n";
      return(0);
    }
    else if (rs_pstar[eqVX] < 0.) {
      //
      // In this case we need to check if the rarefaction straddles x=0 or not.
      //
      double cstar = chydro(rs_pstar, eq_gamma);
      if (rs_pstar[eqVX] < -cstar) {
	//
	// in this case the rarefaction straddles the origin.
	//
	//	cout << "Xi_B <0 :x=0 inside right rarefaction: u*-c* = " << rs_pstar[eqVX]-cstar << "\n";
	//	cout << "Starred rho: " << rs_pstar[eqRO] << "  v: " <<  rs_pstar[eqVX] << "   p: " << rs_pstar[eqPG] << " cstar: " << cstar << "\n";
	// Have to find where x=0, and solve for p,u,rho there.
	// Equation 15,16,17 from my rarefaction.ps notes, with \xi=0
	rs_pstar[eqVX] = (-2.*cr + rs_right[eqVX]*(eq_gamma -1.))/(eq_gamma +1.);
	rs_pstar[eqRO] = rs_right[eqRO] *exp(2./(eq_gamma-1.)*log(-rs_pstar[eqVX]/cr));
	rs_pstar[eqPG] = exp(eq_gamma *log(rs_pstar[eqRO]/rs_right[eqRO]))*rs_right[eqPG];
	//	cout << "Right rarefaction straddles origin." << "\n";
	return(0);
      }
    }
  }
  //
  // Now for shocks...  can a shock wave get swept down stream?  I suppose so.
  // If it does happen, we need to assign the preshock state to the result.
  // Shocks are thin, so they never straddle the origin, so we only need one check
  // per wave.
  //
  if (rs_pstar[eqPG]>1.0000001*rs_right[eqPG]) {
    //
    // right wave is a shock.
    // Check if the shock could get swept to the left so x=0 corresponds to right state
    // This equation from Hirsch p.206 (vol.2) where my vsh is denoted 'C' in his equations.
    //
    double vsh = rs_right[eqVX] +
      (rs_pstar[eqPG]/rs_right[eqPG] -1.)*cr*cr/eq_gamma/(rs_pstar[eqVX]-rs_right[eqVX]);
    if (vsh <0.) {
       rs_pstar[eqPG] = rs_right[eqPG];
       rs_pstar[eqRO] = rs_right[eqRO];
       rs_pstar[eqVX] = rs_right[eqVX];
       return(0);
    }
    //    cout << "right shock, vsh= " << vsh << "\n";
  }
  if (rs_pstar[eqPG]> 1.0000001*rs_left[eqPG]) {
    //
    // left wave is a shock.
    // Check if the shock could get swept to the right so x=0 corresponds
    // to left state.
    //
    double vsh = rs_left[eqVX] + (rs_pstar[eqPG]/rs_left[eqPG] -1.)*cl*cl/eq_gamma/(rs_pstar[eqVX]-rs_left[eqVX]);
    if (vsh >0.) {
      rs_pstar[eqPG] = rs_left[eqPG];
      rs_pstar[eqRO] = rs_left[eqRO];
      rs_pstar[eqVX] = rs_left[eqVX];
      //      cout << "riemann_Euler::solve  Left shock swept to the right!!!" << "\n";
      //      cout << "left shock, vsh= " << vsh << "\n";
      return(0);
    }
    //    cout << "left shock, vsh= " << vsh << "\n";
  }
  //  cout << "No problems with waves overlapping the origin or being swept up. " << "\n";
  return(0);
}

// ##################################################################
// ##################################################################


int riemann_Euler::linearOK()
{
  //
  // Return 0 if linear result is ok, and 1 if it's not.
  //
  /// \section ok Linear or Exact.
  /// Here we have to decide whether to use the linear result or not.
  /// Obviously we'd like to because it's cheaper numerically.\n
  /// I am using the mean matrix constructed from (P_L+P_R)/2. For this averaging,
  /// - contact discontinuities are always treated exactly,
  ///   so the density jump is always correct.
  /// - Pressure: I get 2% error if I use the linear one in the range P_R/P_L in [1/1.5, 1.5].
  /// - Velocity: Error gets very large if (v_L-v_R) is much different from zero, b/c that means
  ///   we have either a converging or diverging flow, and that's harder to handle.  Also when
  ///   v_L/R are close to any of the sound speeds we have trouble because the linear solver 
  ///   isn't good at telling exactly what the wavespeeds are, so it can put you on the wrong
  ///   side of a wave.  This is the largest source of error in the solver.
   
  //
  // First if the left and right states are within some percentage of each other, then we are ok.
  if ( (max(rs_left[eqPG],rs_right[eqPG])/min(rs_left[eqPG],rs_right[eqPG])<1.4) &&
       (max(rs_left[eqRO],rs_right[eqRO])/min(rs_left[eqRO],rs_right[eqRO])<1.4) &&
       ( fabs(rs_right[eqVX]-rs_left[eqVX])/min(cl,cr)<0.03) ) return(0);
  //
  // Next check for failing criteria: in the end I decided it was more
  // trouble than it was worth, so i just say everything that doesn't
  // pass the previous test is going to fail.  This is the simplest
  // thing to do and it works well, although it is very conservative.
  //
  else return(1);
  
#ifdef RSTESTING
  // PRESSURE VARIATIONS
  double phi = max(rs_left[eqPG],rs_right[eqPG]);
  double plo = min(rs_left[eqPG],rs_right[eqPG]);
  if (phi/plo>2.0) {
    //    cout <<"P_L/P_R either too low or too high.\n";
    failplr++;
    return(1);
  }
  if (rs_pstar[eqPG]<plo || rs_pstar[eqPG]>phi) {
    //    cout <<"P* either too low or too high.\n";
    failpst++;
    return(1);
  }
  // VELOCITY VARIATIONS
  // If u* is very small in the linear solution, then it's sign could be wrong...
  if(fabs(rs_pstar[eqVX])/cl < 1.e-2) {
    //     cout <<"u* too small, don't trust linear solution.\n";
    failust++;
    return(1);
  }
  // If we have a rarefaction of any significance then do the full solver:
  if((rs_right[eqVX]-rs_left[eqVX])/min(cl,cr) > 0.05) {
    //    cout <<"Rarefaction\n";
    failrare++;
    return(1);
  }
  // If we have a strong compression, then do the full solver.
  //  cout <<"ccccccccc: "<<left[eqVX]-rs_right[eqVX] <<" ";
  //  cout <<min(cl,cr)<<" "<<cl<<" "<<cr <<"\n";
  if( (rs_left[eqVX]-rs_right[eqVX])/min(cl,cr) > 0.1) {
    //    cout <<"Compression\n";
    failcomp++;
    return(1);
  }
  // DENSITY VARIATIONS
  if(max(rs_left[eqRO]/rs_right[eqRO],rs_right[eqRO]/rs_left[eqRO])>2.) {
    //    cout <<"Density variation too large"<<"\n";
    faildens++;
    return(1);
  }
  // If we pass all the tests, then use the linear one.
  //  cout << "Linear should be ok!"<<"\n";  
  return(0); // Return zero if we don't want to run the exact solver after this.
#endif //RSTESTING
}

// ##################################################################
// ##################################################################


int riemann_Euler::linear_solver()
{
  //
  // First assign the mean values for the primitive variable.
  //
  // Can do a straight average...
  //
  for (int i=0; i<rs_nvar; i++) rs_meanp[i] = (rs_left[i]+rs_right[i])/2.;
  //
  // Or can try the Roe average (NOT CORRECTLY WRITTEN..., FIX THIS!!!
  // IT MIGHT BE BETTER FOR STELLAR WINDS) just need to put in the
  // enthalpy average and get the pressure average from this.
  //
  // double rr = sqrt(rs_right[eqRO]/left[eqRO]);
  // rs_meanp[eqRO] = rr*left[eqRO];
  // rs_meanp[eqVX] = (rr*rs_right[eqVX]+left[eqVX])/(rr+1.);
  // rs_meanp[eqPG] = ???
  
  double mcs = chydro(rs_meanp,eq_gamma);
  //
  // TESTER FOR THE EIGENVALUE SIGNS... (u_av-c_av, u_av, u_av+c_av)
  //
  if (rs_meanp[eqVX]-mcs >=0.) {
    //
    // if(u_av-c_av>0) pstar = left state
    //
    for (int i=0; i<rs_nvar; i++) rs_pstar[i] = rs_left[i];
    return(0);
  }
  else if (rs_meanp[eqVX]+mcs <=0.) {
    //
    // else if(u_av+c_av<=0) pstar = right state
    //
    for (int i=0; i<rs_nvar; i++) rs_pstar[i] = rs_right[i];
    return(0);
  }
  else {
    //
    // else we are in the starred region, so p_star, u_star give the
    // solution.
    //
    rs_pstar[eqPG] = 0.5*(rs_left[eqPG]+rs_right[eqPG] - rs_meanp[eqRO]*mcs*(rs_right[eqVX]-rs_left[eqVX]));
    rs_pstar[eqVX] = 0.5*(rs_left[eqVX]+rs_right[eqVX] - (rs_right[eqPG]-rs_left[eqPG])/rs_meanp[eqRO]/mcs);
    //
    // Now just need to determine rho_star
    // if (u* near zero) rho_star is mean of left and right starred state 
    //
    if (fabs(rs_pstar[eqVX]/mcs) <= 1.e-6) {
      //
      //     cout <<"u*<<1, taking average"<<"\n";
      //
      rs_pstar[eqRO] = rs_meanp[eqRO]*(2.+(rs_left[eqVX]-rs_right[eqVX])/mcs)/2.;
    }     
    else if (rs_pstar[eqVX] >0) {
      //
      // else if(u*>0) rs_pstar = left starred state
      //
      rs_pstar[eqRO] = rs_left[eqRO]+ rs_meanp[eqRO]*(rs_left[eqVX]-rs_pstar[eqVX])/mcs; // rho_(*L)
    }
    else if (rs_pstar[eqVX] <0) {
      //
      // else if(u*<0) pstar = right starred state
      //
      rs_pstar[eqRO] = rs_right[eqRO]+ rs_meanp[eqRO]*(rs_pstar[eqVX]-rs_right[eqVX])/mcs; // rho_(*R)
    }
    else {
      //
      // else bug out
      //
      cerr <<"(Linear Solver) Error in logic!!! Grave."<<"\n";
      return(1);
    }
  }
  return(0);
}

// ##################################################################
// ##################################################################


int riemann_Euler::exact_solver()
{  
  int err=0;
  // ****************************************
  // Exact Solver -- solution is starred state [rho*[L/R],u*,p*]
  // ****************************************
  //
  // Rootfinder. This finds the value of p* which matches u*(L) and u*(R)
  // across the contact discontinuity.
  //
  err += FR_find_root(&(rs_pstar[eqPG]), rs_left[eqPG], rs_right[eqPG], eq_gamma, cl, cr);
  if (err) {
    cout <<"Riemann--exact_solver::FR_find_root: error: " <<rs_left[eqPG];
    cout <<"  "<< rs_right[eqPG] <<"  "<< eq_gamma;
    cout <<"  "<< cl <<"  "<< cr <<"  "<< rs_pstar[eqPG] <<"\n";
  }

  // 
  // ************* Now get u* and rho*left ****************************
  // This gets rs_pstar[eqVX]
  // 
  err += HydroWaveFull(XN, rs_pstar[eqPG], rs_left, &(rs_pstar[eqVX]), &(rs_pstar[eqRO]), eq_gamma);

  //
  // Now get rho*, using either rho*L or rho*R or their arithmetic mean if u* is close to zero.
  // 
  pion_flt rhostar, temp;
  if ( (rs_pstar[eqVX] > 0) && (fabs(rs_pstar[eqVX]/cr)>1.e-6) ) {
    err += HydroWaveFull(XN, rs_pstar[eqPG],  rs_left, &temp, &rhostar, eq_gamma);
  }
  else if ( (rs_pstar[eqVX] < 0) && (fabs(rs_pstar[eqVX]/cr)>1.e-6) ) {
    err += HydroWaveFull(XP, rs_pstar[eqPG], rs_right, &temp, &rhostar, eq_gamma);
  }
  else if (fabs(rs_pstar[eqVX]/cr)<=1.e-6) {
    err += HydroWaveFull(XN, rs_pstar[eqPG],  rs_left, &temp, &rhostar, eq_gamma); 
    err += HydroWaveFull(XP, rs_pstar[eqPG], rs_right, &temp, &(rs_pstar[eqRO]), eq_gamma);
    rhostar = (rhostar +rs_pstar[eqRO])/2.0;
  }
  else {
    // don't want to get here!
    rs_pstar[eqRO] = -1.0;
    cerr << "(riemann_Euler::exact_solver) something went very wrong!!!" << "\n";
    cerr << "P*  rho: " << rs_pstar[eqRO] << "  v: " <<  rs_pstar[eqVX] << "   p: " << rs_pstar[eqPG] << "\n";
    rep.printVec("rs_left", rs_left, rs_nvar); 
    rep.printVec("rs_right", rs_right, rs_nvar); 
    return(1);
  }
  rs_pstar[eqRO] = rhostar;

  //
  // Make sure everything went ok.
  // 
  if (err != 0) {
    cerr << "(riemann_Euler::solve) err!=0, something went wrong.! " << "\n";
    rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar [eqVX] = -1.9;
    return(1);
  }


  //
  // Now we check wave locations to get the state at x=0.
  // 
  err = check_wave_locations();
  if (err != 0) {
    cerr << "(riemann_Euler::solve) checking rarefactions/shocks didn't work! " << "\n";
    rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar [eqVX] = -1.9e98;
    return(1);
  }
  return(0);
}

// ##################################################################
// ##################################################################


int riemann_Euler::solve_rarerare()
{
  ///
  /// This function is called only if I know for sure there are two rarefactions,
  /// and no cavitation.  In this case the solution can be found in closed form
  /// analytically (see e.g. my notes on exact riemann solvers).
  /// 
  /// The condition for definitely having two rarefactions it that 
  /// \f[ \frac{2}{\gamma -1} \left( c_L + \sqrt{\frac{\gamma -1}{2\gamma}}c_R \right)
  ///     \leq u_R-u_L \leq \frac{2}{\gamma -1}(c_R+c_L) \f]
  ///
  /// \section p_equation Pressure equation 
  /// The resolved state pressure is given analytically as 
  /// \f[ p_{*} = \left[ \left(c_L+c_R -\frac{\gamma-1}{2}(u_R-u_L)\right)
  /// \left(\frac{c_L}{p_L^{(\gamma-1)/2\gamma}} +\frac{c_R}{p_R^{(\gamma-1)/2\gamma}}\right)^{-1}\right]^{2\gamma/(\gamma-1)}
  /// \f]
  ///

  rs_pstar[eqPG] = pow(  (cl+cr -(eq_gamma-1.)/2.*(rs_right[eqVX]-rs_left[eqVX])) 
		    /( (cl*exp(-(eq_gamma-1.)/2./eq_gamma*log(rs_left[eqPG]))) 
		       +(cr*exp(-(eq_gamma-1.)/2./eq_gamma*log(rs_right[eqPG]))) ), 2.*eq_gamma/(eq_gamma-1.));  
  //
  // Velocity Equation: v*-vl = 2cl/(g-1) *[1 -(ppl)^((g-1)/2g)]
  //
  rs_pstar[eqVX] = rs_left[eqVX] + 2.*cl/(eq_gamma-1.) *(1. -exp((eq_gamma-1.)/2./eq_gamma*log(rs_pstar[eqPG]/rs_left[eqPG])));
  //
  // Density equation: use rho*L/R depending on sign of v*
  // Solve by assuming adiabatic gas, and get answer from p*,rho_0,p_0
  //
  if ( (rs_pstar[eqVX]>0) && (fabs(rs_pstar[eqVX]/cr) > 1.e-6) ) {
    rs_pstar[eqRO] = rs_left[eqRO]*exp(log(rs_pstar[eqPG]/rs_left[eqPG])/eq_gamma);
  }
  else if ( (rs_pstar[eqVX]<0) && (fabs(rs_pstar[eqVX]/cr) > 1.e-6) ) {
    rs_pstar[eqRO]  = rs_right[eqRO]*exp(log(rs_pstar[eqPG]/rs_right[eqPG])/eq_gamma);
  }
  else if (fabs(rs_pstar[eqVX]/cr) <= 1.e-6){
    rs_pstar[eqRO]  = ( (rs_right[eqRO]*exp(log(rs_pstar[eqPG]/rs_right[eqPG])/eq_gamma))
		   +( rs_left[eqRO]*exp(log(rs_pstar[eqPG]/rs_left[eqPG])/eq_gamma)) )/2.0;
  }
  else{
    //don't want to get here!
    rs_pstar[eqRO] = -1.0;
    cerr << "(riemann_Euler::rare_rare) something went very wrong!!!" << "\n";
    return(1);
  }
  //
  // Now check that the starred state is valid at the origin.
  // For rarefactions, the wave can straddle the origin, so we need to check that it doesn't!!!
  //
  int err = check_wave_locations();
  if (err != 0) {
    cerr << "(riemann_Euler::solve) checking rarefactions/shocks didn't work! " << "\n";
    rs_pstar[eqPG] = rs_pstar[eqRO] = rs_pstar [eqVX] = -1.9;
    return(1);
  }
  return(0);
}


// ##################################################################
// ##################################################################


int riemann_Euler::solve_cavitation()
{
  /// This only gets called if I know I have cavitation.
  /// The solution for this is in my exact riemann solver notes.
  /// I go through the cases one by one, for where each wave is.  I know
  /// the locations of all the boundaries b/c they only depend on the left 
  /// and right states.\n
  /// If I always check the right boundary of the region, I only need one
  /// check per region.
  /// 
  /// The condition for cavitation is that 
  /// \f[ (u_R-u_L) >= \frac{2}{\gamma -1}(c_R+c_L) \f]
   
  if( (rs_left[eqVX]-cl)>=0. ) {
    //
    // Everything swept to the right.
    //
    //cout <<" Everything swept to the right..."<<"\n";
    for (int i=0;i<rs_nvar;i++)
      rs_pstar[i] = rs_left[i];
    return(0);
  }
  double temp = 2./(eq_gamma-1.);
  if ( (rs_left[eqVX]+temp*cl)>=0. ) {
    //
    // Inside the left rarefaction (See my rarefaction notes).
    //
    //cout <<"Inside left rarefaction..."<<"\n";
    rs_pstar[eqVX] = (2.*cl + rs_left[eqVX]*(eq_gamma -1.))/(eq_gamma +1.);
    rs_pstar[eqRO] = rs_left[eqRO] *exp(2./(eq_gamma-1.)*log(rs_pstar[eqVX]/cl));
    rs_pstar[eqPG] = exp(eq_gamma *log(rs_pstar[eqRO]/rs_left[eqRO]))*rs_left[eqPG];
    return(0);
  }
  if ( (rs_right[eqVX]-temp*cr)>=0. ) {
    //
    // Inside the cavity.
    //
    //cout <<"Inside the cavity..."<<"\n";
    rs_pstar[eqRO] = eq_refvec[eqRO]*BASEPG;
    rs_pstar[eqPG] = eq_refvec[eqPG]*BASEPG;
    rs_pstar[eqVX] = eq_refvec[eqVX]*BASEPG;  // Not sure if this is right or not!  but it doesn't matter.
    return(0);
  }
  if ( (rs_right[eqVX]+cr)>0. ) {
    //
    // Inside the right rarefaction (See my rarefaction notes).
    //
    //cout <<"Inside the right rarefaction..."<<"\n";
    rs_pstar[eqVX] = (-2.*cr + rs_right[eqVX]*(eq_gamma -1.))/(eq_gamma +1.);
    rs_pstar[eqRO] = rs_right[eqRO] *exp(2./(eq_gamma-1.)*log(-rs_pstar[eqVX]/cr));
    rs_pstar[eqPG] = exp(eq_gamma *log(rs_pstar[eqRO]/rs_right[eqRO]))*rs_right[eqPG];
    return(0);
  }
  if ( (rs_right[eqVX]+cr)<=0. ) {
    //
    // Everything swept to the left.
    //
    //cout <<"Everything swept to the left..."<<"\n";
    for (int i=0;i<rs_nvar;i++)
      rs_pstar[i] = rs_right[i];
    return(0);
  }
  //
  // Don't want to get here!
  //
  cerr << "(riemann_Euler::solve_cavitation) something went very wrong!!!" << "\n";
  rep.printVec("rs_pstar",rs_pstar,rs_nvar);
  cout <<"sound speed = "<<cl<<", "<<cr<<"\n";
  return(1);
}


// ##################################################################
// ##################################################################



