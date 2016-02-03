/// - 2011.01.14 JM: moved to microphysics/ sub-dir.

#ifndef INTEGRATOR_H
#define INTEGRATOR_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


/** \brief Integration class, with a bunch of different integration methods in it.
 * 
 * List includes:\n
 *  - N-point Euler Integration.
 *  - Adaptive RK5 method.
 *
 * This is an abstract base class, with the function dPdt() not defined, 
 * and it is intended that any class that wants to use it will define the
 * rate function.  e.g. the MicroPhysics class will have dPdt() defined.
 * */
class Integrator_Base {
 public:
  Integrator_Base();
  virtual ~Integrator_Base();
  int Set_Nvar(int ///< length of state vectors.
	       );
  /** \brief Rate function defining the equation we are integrating.
   * This is not defined in Integrator_Base, but should be in a derived
   * class wanting to integrate stuff.
   *
   * This version is for when dP/dt depends explicitly on P only.
   */
  virtual int dPdt(const int,      ///< length of state vector (for checking).
		   const double *, ///< current state vector P.
		   double *        ///< Rate vector to write to, R=dPdt(P)
		   )=0;
  /** \brief This is for if we are solving the rate equation, and returns the
   * creation rate of some quantity. xdot=A*(1-x)-B*x, so this returns A(x).
   */
  virtual int C_rate(const int,      ///< length of state vector.
		     const double *, ///< current state vector P.
		     double *        ///< Creation rate vector to write to.
		     )=0;
  /** \brief This is for if we are solving the rate equation, and returns the
   * destruction rate of some quantity. xdot=A*(1-x)-B*x, so this returns A(x)+B(x).
   */
  virtual int D_rate(const int,      ///< length of state vector.
		     const double *, ///< current state vector P.
		     double *        ///< Destruction rate vector to write to.
		     )=0;

  /** \brief Do an n-point 1st order Euler integration of P(t) to P(t+dt). */
  int Int_Euler(const int,      ///< number of elements in P array.
		const double *, ///< value of P at initial value of t.
		const double,   ///< initial value of t.
		const double,   ///< Total step dt to take.
		const int,      ///< Number of integration steps to take.
		double *        ///< pointer to final P value.
		);
  /** \brief Do Euler integrations, doubling number of points until we get to
   * the required relative accuracy. 
   *
   * This does some error checking: if any element in the returned state is NAN,
   * then it doubles the number of points and tries again.  If it has to try more than
   * 2^24 points, it returns the best guess and fails, and the best guess may contain 
   * NANs.
   **/
  int Int_DumbAdaptive_Euler(const int ,     ///< number of elements in P array.
			     const double *, ///< value of P at initial value of t.
			     const double,   ///< initial value of t.
			     const double,   ///< Total step dt to take.
			     const double,   ///< error tolerance per step.
			     double *,       ///< pointer to final P value.
			     double *        ///< pointer to final t value.
			     );
  /** \brief Take a single 4th order Runge-Kutta step. 
   * 
   * Pointer to final state can be same as pointer to input state.
   * This is based on the algorithm description in Numerical Recipes in C (1992), ch16.1.
   * It is assumed that the Chemistry class contains a function to calculate the 
   * rate dPdt().
   */
  int Step_RK4(const int,      ///< number of elements in P array.
	       const double *, ///< value of P at initial value of t.
	       const double,   ///< initial value of t.
	       const double,   ///< Total step dt to take.
	       double *       ///< pointer to final P value.
	       );
  /** \brief take a single 5th order Cash-Karp Runge-Kutta step.  
   * 
   * Pointer to final state can be same as pointer to input state.
   * This is based on the algorithm description in Numerical Recipes in C (1992), ch16.1.
   * It is assumed that the Chemistry class contains a function to calculate the 
   * rate dPdt().
   */
  int Step_RK5CK(const int,      ///< number of elements in P array.
		 const double *, ///< value of P at initial value of t.
		 const double,   ///< initial value of t.
		 const double,   ///< stepsize, dt.
		 double *,       ///< final value of P
		 double *       ///< error estimate.
		 );
  /** \brief function to take a single RKCK step to a required accuracy.  
   * 
   * Pointer to final state can be same as pointer to input state.
   * This is based on the algorithm description in Numerical Recipes in C (1992), ch16.2.
   * It is assumed that the Chemistry class contains a function to calculate the 
   * rate dPdt().
   */
  int Stepper_RKCK(const int,      ///< number of elements in P array.
		   const double *, ///< value of P at initial value of t.
		   const double,   ///< initial value of t.
		   const double,   ///< stepsize to try.
		   const double,   ///< allowed error tolerance.
		   double *, ///< final value of P(t)
		   double *, ///< actual step taken.
		   double *  ///< estimate of next step to try.
		   );
  /** \brief Do an adaptive 5th order Cash-Karp integration to a given relative accuracy. 
   * 
   * Pointer to final state can be same as pointer to input state.
   * This is based on the algorithm description in Numerical Recipes in C (1992), ch16.2.
   * It is assumed that the Chemistry class contains a function to calculate the 
   * rate dPdt().
   */
  virtual int Int_Adaptive_RKCK(const int,      ///< number of elements in P array.
			const double *, ///< value of P at initial value of t.
			const double,   ///< initial value of t.
			const double,   ///< Total step dt to take.
			const double, ///< Required fractional accuracy.
			double *,  ///< pointer to final P value.
			double *  ///< pointer to final t value.
			);
  /** \brief Do an iterative integration of a system of equations, treating
   * coefficients as constants.
   *
   * I'm not sure how to make this work with an array of variables, one of which is 
   * the internal energy.  So don't use this for now.
   */
  int Int_Iterative_FC(const int,      ///< number of elements in P array.
		       const double *, ///< value of P at initial value of t.
		       const double,   ///< initial value of t.
		       const double,   ///< Total step dt to take.
		       double * ///< pointer to final P value.
		       );
  /** \brief Do a (sort of) Backwards Differencing subcycling integration.
   *
   * I'm not sure how to make this work with an array of variables, one of which is 
   * the internal energy.  So don't use this for now.
   **/
  int Int_Subcycle_BDF(const int,///< number of elements in P array.
		       double *, ///< value of P at initial value of t.
		       double,   ///< initial value of t.
		       double,   ///< Total step dt to take.
		       double,   ///< subcycling parameter: dt_sub = dPtol*P/dPdt
		       double * ///< pointer to final P value.
		       );
 private:
   double a2,a3,a4,a5,a6;
   double b21, b31,b32, b41,b42,b43, b51,b52,b53,b54, b61,b62,b63,b64,b65;
   double c1,c3,c4,c6;
   double dc1,dc3,dc4,dc5,dc6;
   int int_nvar;
   double *k1,*k2,*k3,*k4,*k5,*k6;
};

#endif // INTEGRATOR_H
