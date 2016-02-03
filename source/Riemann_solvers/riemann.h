/// \file riemann.h
/// \brief Class definition for Hydrodynamics Riemann Solver
/// 
///  \author Jonathan Mackey
/// 
/// Written 2006-11-09-Thurssday
/// Modified :2006-12-20 Testing for linear solver added.
///           2007-01-12 Lots more testing added; linearOK() changed.
/// 
/// This is a Hydrodynamics Riemann Solver for the Euler Equations.
/// It works by first calling a linear solver, and if the results
/// are within certain tolerances it returns the answer, but if not
/// it calls an exact solver and returns the result from that.
/// 
/// Modifications: 
///  - 2006-06-27 Changed counters to long ints, so there'e no wrap-around.
///  - 2007-07-13 New Functions: SetDirection() and SetAvgState()
///  - 2009-10-20 Split hydro and MHD solver, and put them in a new class heirarchy.
/// - 2010.12.23 JM: Moved to Riemann_solver/ directory.  Moved MHD
///   enums to riemannMHD.h.  Got rid of riemann_base class from 
///   eqns_bsae.h so that these don't derive from it.  Got rid of the
///   inherited data arrays RS_pstar/left/right since they have caused
///   bugs.
/// - 2011.03.03 JM: Added rs_nvar=5 for local state vectors.  New code versions
///    can handle up to 70 tracers, so it would hugely slow down the code if the
///    Riemann solver used all that memory when it only needs 5 vars.  Tracer 
///    fluxes are dealt with by the flux-solver classes.
///

#ifndef RIEMANN_H
#define RIEMANN_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "findroot.h"
#include "equations/eqns_hydro_adiabatic.h"
//#define RSTESTING ///< If testing the Riemann Solver.



/// \brief Hydrodynamics Riemann Solver.
/// 
/// This class solves the Riemann Problem for the Euler Equations in 
/// 1, 2, or 3D.  (It is always a 1D problem, but the other dimensions
/// can be added in and can vary across the contact discontinuity).  
/// 
/// The linear solver is based on Toro (1993) 
/// and Falle, Komissarov \& Joarder (1998), and is simple and fast.  I think
/// it would be classes as a Roe-type solver, but it doesn't use the Roe
/// averaging.
/// 
/// The exact solver  is base on an idea Andy had (which I think he uses 
/// for his exact MHD solver).  You have the left and right states, and make
/// a guess for the resolved state pressure.  With this you calculate the 
/// resolved state velocity from the left and right waves respectively, and
/// iterate until these velocities match.  Then you have the resolved state
/// and can calculate the wave locations to get the state at x=0.
/// 
/// References:\n
/// Toro, E.F., 1993, in 'Numerical Methods for Fluid Dynamics 4' 
/// (Oxford, U.K.: Clarendon Press), Eds. M.J. Baines & K.W. Morton, p.375.\n
/// Falle, S.A.E.G., Komissarov, S.S., & Joarder, P., 1998, M.N.R.A.S., 297, 265.\n
/// 
/// \note
/// This class is designed for single-threaded programming.  It is set
/// up at the start of a simulation, and then called many times in sequence
/// to solve riemann problems between grid cells.  It has no checking to make
/// sure it is not already doing a solve when a new solve is called.
/// 
/// \section Mode
/// I have a switch in the constructor to say
/// if you want it to always use a linear solver, or always an exact solver, or
/// a combination solver (linear for the easy problems, and exact for the hard
/// ones).  So:
///  - Mode 1 : Always Linear Solves.
///  - Mode 2 : Always Exact Solves.
///  - Mode 3 : Always a combination.
/// 
/// As far as I can tell, the difference in shock tube problems between all 
/// these modes is small, of the order of one percent.
///
class riemann_Euler : virtual public eqns_Euler,
                      virtual public findroot
{
  public:
  /// Constructor; assigns memory for arrays.
  /// 
  /// Assumes an ideal gas equation of state.  The eos gamma
  /// is passed in to the public solve function.
  ///
  riemann_Euler(
      const int,      ///< Length of State Vectors, nvar
      const pion_flt *, ///< Mean values of primitive variables on grid [vector, length nvar]
      const double    ///< Gamma for state vector.
      );

  /// \brief Destructor: deletes dynamically allocated data.
  ~riemann_Euler();

  ///
  /// Gets the solution to the Riemann Problem, either a linearised,
  /// exact, or hybrid solution depending on the 4th argument.
  ///  
  int JMs_riemann_solve(
      const pion_flt *, ///< Left Primitive var. state vector.
      const pion_flt *, ///< Right Primitive var. state vector.
      pion_flt *,       ///< Result Primitive var. state vector.
      const int,      ///< Solve Type (1=LinearRS,2=ExactRS,3=HybridRS)
      const double    ///< Gas constant gamma.
      );

   /// \brief Prints out info on what the solver has done, and how well it's done it.
   void testing();
   
 protected:
 private:
   const int rs_nvar; ///< length of state vectors required for the Riemann solver, i.e.5.
   pion_flt *rs_left;  ///< Local copy of the left state.
   pion_flt *rs_right; ///< Local copy of the right state.
   pion_flt *rs_pstar; ///< Local copy of result state.
   pion_flt *rs_meanp; ///< The average state.
   long int linct; ///< Counter for how many linear solves we have done.
   long int exact; ///< Counter for how many exact solves we have done.
   long int total; ///< counter for how many solves the instance has done.
   long int samestate; ///< Debugging counter, for counting how many solves had the same left and right state.
#ifdef RSTESTING
   pion_flt *linearpstar; ///< (TESTING) Testing variable, for the result state from the linear solver.
   int fails; ///< Counter for how many linear solves weren't good enough when linearOK() thought they should have been.
   int failplr, failpst, failust, failrare, failcomp, faildens;
   double lineartol; ///< (TESTING) How closely we require the linear solution to match the exact one.
#endif   // RSTESTING
   double cl;  ///< The hydro sound speed in the left state.
   double cr;  ///< The hydro sound speed in the right state.
   
  ///
  /// Re-definition of the public root-finding function, for solving the jump conditions
  /// across hydrodynamic discontinuities.
  ///
  virtual int FR_find_root(
      pion_flt *, ///< pointer to result 
      const pion_flt,  ///< parameter 1
      const pion_flt,  ///< parameter 2
      const pion_flt,  ///< parameter 3
      const pion_flt,  ///< parameter 4
      const pion_flt   ///< parameter 5
      );

  ///
  /// Re-definition of the funtion to find the root of.
  /// Input an x value and this returns f(x).  The aim of this class
  /// is to find the x value for which f(x)=0, so we are trying to
  /// minimize the return value of this function.
  ///
  virtual pion_flt FR_root_function(
      const pion_flt  ///< x-value
      );
   

  /// \brief The hydrodynamics 1D linear riemann solver.
  /// 
  /// Simple hydro linear riemann solver.
  /// Currently uses straight average for the average state.
  /// \retval 0 success
  /// \retval 1 failure
  /// 
  int linear_solver();
   
   /// \brief Decides whether or not the linear solver result is ok.
   /// 
   /// This is an important part of the solver, as I want the check to 
   /// be fast, but also to catch problem cases well, and also to make sure
   /// that the overwhelming majority of solves are linear and not exact
   /// (for speed).
   /// 
   /// \note I am currently testing a number of ideas.
   /// \retval 0 success
   /// \retval 1 failure
   /// 
   int linearOK();

   /// \brief The hydrodynamics 1D exact riemann solver.
   /// 
   /// The best iterative exact solver I can come up with.  It is 
   /// very robust, but doesn't handle cavitation.  The public solve
   /// function is expected to catch cavitation and treat it separately
   /// by calling solve_cavitation().
   /// 
   /// This solver returns the resolved state, and chooses the value of 
   /// rho_star depending on whether the contact moves to the left or 
   /// to the right.\n
   /// It then calls check_wave_locations(), which checks where the resolved
   /// state is w.r.t. the origin, and gets the correct state at x=0 if it is
   /// not the resolved state (e.g. if x=0 is within the rarefaction).
   /// \retval 0 success
   /// \retval 1 failure
   /// 
   int exact_solver();
 
   /// \brief Given resolved state (from exact solver) and left/right states, checks where x=0 is.
   /// 
   /// This function takes p* and the left and right states, and checks where x=0
   /// is in relation to the waves.  For example if a shock is swept downstream 
   /// past the origin by a strong flow, then the solution will not be the starred
   /// state but the left or right state.  This function finds where x=0 is in 
   /// relation to the waves, and assigns pstar[] to the appropriate value.  It can
   /// open up rarefaction fans and put in the exact values.
   /// \retval 0 success
   /// \retval 1 failure
   /// 
   int check_wave_locations();

   /// \brief Solves for when I know we have cavitation.
   /// 
   /// There is a simple test for cavitation, and if it occurs, then you can 
   /// just write down the solution without doing an iterative solve.  So that's
   /// what this function does.  If you know you have cavitation, call this function
   /// and it will give the correct solution at x=0.
   /// \retval 0 success
   /// \retval 1 failure
   /// 
   int solve_cavitation();

   /// \brief Solves for when I know we have a double rarefaction.
   /// 
   /// There is a condition for the left and right states which tells us that we
   /// definitely have a double rarefaction.  In this case there is an exact closed-
   /// form solution which we can write down.  So if you evaluate this condition and
   /// determine that you have double rarefaction (and not cavitation, which is a 
   /// more restictive condition), then call this function and it gives back the 
   /// correct solution at x=0.
   /// \retval 0 success
   /// \retval 1 failure
   /// 
   int solve_rarerare();
};

#endif

