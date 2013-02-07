///
/// \file cvode_integrator.h
/// \author Jonathan Mackey
///
/// This file contains a class which is intended as a general
/// integrator of microphysics equations, which can be inherited by
/// specific implementations with a number of species and different
/// heating/cooling rates, and even integration of the optical depth
/// for the implicit C2-ray-type method.
///
/// The integration method uses the CVODE solver from the SUNDIALS
/// package by (Cohen, S. D., & Hindmarsh, A. C. 1996, Computers in
/// Physics, 10, 138) available from 
///   https://computation.llnl.gov/casc/sundials/main.html
/// The method is backwards differencing with Newton iteration.
///
/// Modifications:
/// - 2011.10.06 JM: Wrote file, based on old code in mp_v2_aifa.h and a test
///   integrator in active/code_misc/sundials/test_prog3/
/// - 2013.02.07 JM: Changed int to long int in Jacobian function for
///    compatibility with sundials 2.5.0.
///

#ifndef CVODE_INTEGRATOR_H
#define CVODE_INTEGRATOR_H


//#define CVODE_DEBUG

//
// The native cvode solver seems to be faster than the Lapack one.  But for
// future applications with larger networks that may no longer be true, so 
// there is an option here to use the LAPACK dense solver by enabling the 
// following #define.
//
//#define LAPACK

//
// Header files with a description of contents used
//
#include <cvode/cvode.h>             // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h>  // serial N_Vector types, fcts., and macros

#ifdef LAPACK
// LAPACK SOLVER
#include <cvode/cvode_lapack.h>      /* prototype for CVLapackDense */
#else
// NON-LAPACK SOLVER
#include <cvode/cvode_dense.h>       // prototype for CVDense
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
#endif


//
// ydot needs to be a stand-alone function for cvode, so this function will
// just provide the interface to call solver::ydot();
//
int Ydot_for_cvode(
          double,   ///< current time
          N_Vector, ///< current Y-value
          N_Vector, ///< vector for Y-dot values
          void *    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          );

int Jacobian_for_cvode(
          long int, ///< N (not sure what this is for! Must be internal)
          double,   ///< time, t
          N_Vector, ///< y
          N_Vector, ///< ydot
          DlsMat,   ///< Jacobian matrix
          void *,   ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          N_Vector, ///< temp vector, must be for internal use
          N_Vector, ///< temp vector, must be for internal use
          N_Vector  ///< temp vector, must be for internal use
          );


class cvode_solver {

  public:
  
  friend int Ydot_for_cvode(
          double,   ///< current time
          N_Vector, ///< current Y-value
          N_Vector, ///< vector for Y-dot values
          void *    ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          );
  
  friend int Jacobian_for_cvode(
          int,    ///< N (not sure what this is for! Must be internal)
          double, ///< time, t
          N_Vector, ///< y
          N_Vector, ///< ydot
          DlsMat,   ///< Jacobian matrix
          void *,   ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          N_Vector, ///< temp vector, must be for internal use
          N_Vector, ///< temp vector, must be for internal use
          N_Vector  ///< temp vector, must be for internal use
          );

  cvode_solver();
  ~cvode_solver();

  /// This initialises the solver with a user-supplied Jacobian function.
  ///
  int setup_cvode_solver();

  ///
  /// This initialises the solver without using a Jacobian function.
  /// This is the function I have used so far.
  ///
  int setup_cvode_solver_without_Jacobian();

  ///
  /// This takes a step dt, returning a non-zero error code if the error fails.
  ///
  int integrate_cvode_step(
              N_Vector, ///< input vector (may be overwritten during integration!)
              void *,  ///< parameters for user_data (pointer to instance of this class!)
              double,  ///< start time.
              double,  ///< time-step.
              N_Vector ///< output vector.
              );

  private:
  N_Vector abstol; ///< vector of absolute error tolerances for y-elements.
  void *cvode_mem; ///< pointer to memory allocation for the solver.
  int n_eq; ///< number of equations to solve.
  int n_xd; ///< number of elements in user-data array.
  bool have_setup_cvodes; ///< flag to make sure we only set up CVODES once.

  //---------------------------------------------------------------------------
  //------------ STUFF TO BE DEFINED IN DERVIED CLASS FOLLOWS -----------------
  //---------------------------------------------------------------------------
  public:
  ///
  /// calculate dy/dt for the vector of y-values (NOT IMPLEMENTED HERE).
  ///
  virtual int ydot(
          double,         ///< current time (probably not needed for rate equations)
          const N_Vector, ///< current Y-value
          N_Vector,       ///< vector for Y-dot values
          const double *  ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          )=0;

  ///
  /// Calculate the Jacobian matrix d(dy_i/dt)/dy_j for a vector of y-values.
  /// (NOT IMPLEMENTED HERE).
  ///
  virtual int Jacobian(
          int,            ///< N (not sure what this is for! Must be internal)
          double,         ///< time, t
          const N_Vector, ///< current Y-value
          const N_Vector, ///< vector for Y-dot values
          const double *, ///< extra user-data vector, P, for evaluating ydot(y,t,p)
          DlsMat          ///< Jacobian matrix
          )=0;

  ///
  /// Get the number of extra parameters and the number of equations.
  /// (NOT IMPLEMENTED HERE).
  ///
  virtual void get_problem_size(int *, ///< number of equations
                        int *  ///< number of parameters in user_data vector.
                        )=0;

  protected:
  ///
  /// set the relative and absolute error tolerances
  ///
  virtual void get_error_tolerances(
          double *, ///< relative error tolerance (single value)
          double []  ///< absolute error tolerance (array)
          )=0;

  //---------------------------------------------------------------------------
  //------------ END OF STUFF TO BE DEFINED IN DERVIED CLASS ------------------
  //---------------------------------------------------------------------------
};


#endif // CVODE_INTEGRATOR_H



