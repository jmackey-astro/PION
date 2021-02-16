/// \file findroot.h
///
/// \brief Class definition for findroot class
///
/// \author Jonathan Mackey
///
/// Declarations of the equation classes and the rootfinding
/// class.  This class was written for finding roots to
/// equations encountered in solving the exact Riemann Problem
/// in hydrodynamics, but should be general enough to be applied
/// to find the roots of any one-dimensional equation.
///
/// Modifications:
/// - 2009-10-23 made Riemann_Euler inherit from findroot, and
///   redefine the function to get the root of!
///   So now findroot can't be set up as is, a derived class needs to
///   define the function to get the root of.
/// - 2010.12.23 JM: Moved to Riemann_solver/ directory.
/// - 2015.08.03 JM: Added pion_flt for pion_flt* arrays (allow floats)

#ifndef FINDROOT_H
#define FINDROOT_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

/// \brief The findroot class: Use this to find roots of equations by various
/// methods.
///
/// This class can find the roots of equations by either bisection or by Brent's
/// method. Bisection is slow but sure; Brent's method is from Numerical Recipes
/// (p.361 of ansi-c version).  I am using their zbrent() subroutine.
///
/// The class is designed so that *eqn is a pointer to a base equation class,
/// and you can define any derived equation class you want, and then when you
/// call a solve just make eqn point to your equation and call solve() and you
/// should get the answer out.
///
/// Any new solve_(problem) functions must set *eqn to point to a proper
/// function to find the roots of. It must also take pion_flt *ans as an
/// argument, and any data parameters it needs, which should be put in the
/// declaration of FR_base_eqn. Then it calls the private function
/// findroot::solve(ans); to get the answer.
///
/// The idea is that all the machinery of solving equations is in private
/// functions, and you just define a new public interface function which assigns
/// *eqn to your equation and calls either solve_pm() or solve_pos() to get the
/// solution.  The rest of the code should not ever have to be modified, unless
/// you want to add faster root-finding methods.
///
class findroot {
public:
  /// \brief Constructor
  ///
  /// The constructor just sets the required accuracy to be \f$10^{-8}\f$.
  ///
  findroot();
  /// Trivial destructor (no dynamically allocated memory)
  ~findroot();

  /// \brief My wrapper function to solve the hydro Riemann Problem pressure
  /// equation.
  ///
  /// This function finds a positive root to the function FR_root_function(),
  /// which must be defined in a derived class.
  ///
  virtual int FR_find_root(
      pion_flt*,       ///< pointer to result
      const pion_flt,  ///< parameter 1
      const pion_flt,  ///< parameter 2
      const pion_flt,  ///< parameter 3
      const pion_flt,  ///< parameter 4
      const pion_flt   ///< parameter 5
  );
  /// \brief This function doesn't actually do anything useful
  ///
  /// It is there for testing only.  Has a simple equation which I can find
  /// the roots of, and a variable to assign data into.
  ///
  int solve_test(
      pion_flt*,      ///< pointer to result.
      const pion_flt  ///< param1.
  );

protected:
  /// \brief The general (private) solver function to get the roots of your
  /// equation
  ///         in the range \f$ [-\infty, \infty]\f$
  ///
  /// Once the public solve has been called, the parameters will be set, and
  /// *eqn will point to the right equation, and solve just gets the roots,
  /// given initial guesses for bracketing the root [x1,x2] It first brackets
  /// the root (if possible), and then calls a findroot to find the root.
  /// \retval 0 Success
  /// \retval 1 Failure (with error message hopefully).
  ///
  int solve_pm(
      pion_flt,  ///< Initial bracketing value guess x1
      pion_flt,  ///< Initial bracketing value guess x2
      pion_flt*  ///< Pointer to ans.
  );

  /// \brief Same as solve_pm() but int the range \f$ [0, \infty]\f$.
  ///
  /// \retval 0 Success.
  /// \retval 1 Failure.
  ///
  int solve_pos(pion_flt, pion_flt, pion_flt*);

  ///
  /// Input an x value and this returns f(x).  The aim of this class is to
  /// find the x value for which f(x)=0, so we are trying to minimize the
  /// return value of this function.
  ///
  virtual pion_flt FR_root_function(const pion_flt  ///< x-value
                                    ) = 0;
  ///
  /// Various parameters we can use for this class.
  ///
  pion_flt FR_param1, FR_param2, FR_param3, FR_param4, FR_param5;
  pion_flt FR_var1, FR_var2, FR_var3;

private:
  pion_flt errtol;  ///<  The error tolerance of the solution, set to 1e-8 by
                    ///<  default.
  ///
  /// Test function for testing the root-finding algorithm.
  ///
  pion_flt FR_test_function(const pion_flt  ///< x-value
  );

  /// \brief Brackets a root in the range \f$ [0, \infty]\f$
  ///
  /// This function returns brackets around a root.
  /// If it cannot find a root in 50 expansions it tries taking the lower
  /// bound to be exactly zero, and if this still doesn't work it returns
  /// a failing value.
  /// \retval 0 Success.
  /// \retval 1 Failure.
  ///
  int bracket_root_pos(
      pion_flt*,  ///< The lower bracketing value x1
      pion_flt*   ///< The upper bracketing value x2
  );

  /// \brief Brackets a root in the range \f$ [-\infty, \infty]\f$
  ///
  /// This function returns brackets around a root, given initial guesses.
  /// It will perform up to 50 geometric expansions of the range by
  /// a given factor, and if it can't bracket the root it will return
  /// a failing value.
  /// \retval 0 Success.
  /// \retval 1 Failure.
  ///
  int bracket_root_pm(
      pion_flt*,  ///< The lower bracketing value x1
      pion_flt*   ///< The upper bracketing value x2
  );

  /// \brief Given brackets around a root, finds the root by the bisection
  /// method.
  ///
  /// This method is good in that it cannot fail if the root is initially
  /// bracketed. It is slow, however, and almost always takes the same number
  /// of iterations to get a solution.  Brent's method is about 20 to 30\%
  /// faster.
  ///
  /// \retval 0 success.
  /// \retval 1 fail.
  ///
  int find_root_bisection(
      pion_flt*,  ///< lower bound x1
      pion_flt*,  ///< upper bound x2
      pion_flt,   ///< required accuracy
      pion_flt*   ///< pointer to answer
  );

  /// \brief Given brackets around a root, finds the root by the Brent's
  /// method.
  ///
  /// This method is fast and robust, and I haven't found it to fail yet.
  /// Got it from Numerical Recipes, p.361 subroutine zbrent()
  ///
  /// \retval 0 success.
  /// \retval 1 fail.
  ///
  int find_root_zbrent(
      pion_flt x1,   ///< lower bound x1
      pion_flt x2,   ///< upper bound x2
      pion_flt err,  ///< required accuracy
      pion_flt* ans  ///< pointer to answer
  );
};

#endif
