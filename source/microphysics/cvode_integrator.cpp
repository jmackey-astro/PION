///
/// \file cvode_integrator.cpp
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
/// - 2011.10.06 JM: Wrote file, based on old code in mp_v2_aifa.cc
///    and a test integrator in active/code_misc/sundials/test_prog3
/// - 2011.10.17 JM: Debugging.
/// - 2012.09.28 JM: Added separators between functions for clarity.
/// - 2013.02.07 JM: Changed int to long int in Jacobian function for
///    compatibility with sundials 2.5.0.

#include "microphysics/cvode_integrator.h"
#include <cmath>
#include <iostream>
using namespace std;

#define MAX_TRYS 50
#define MAX_STEPS 100

// ##################################################################
// ##################################################################

cvode_solver::cvode_solver()
{
  abstol    = 0;
  cvode_mem = 0;
  n_eq      = -1;  // number of equations i.e. number of species in Y.
  n_xd = -1;  // number of extra paramaters needed in ydot/jacobian functions
              // (maybe none).
  have_setup_cvodes = false;
}

// ##################################################################
// ##################################################################

cvode_solver::~cvode_solver()
{
#if defined CVODE2
  N_VDestroy_Serial(abstol);
#else
  N_VDestroy(abstol);
#endif
  CVodeFree(&cvode_mem);
}

// ##################################################################
// ##################################################################

int cvode_solver::setup_cvode_solver()
{
  if (have_setup_cvodes) {
    cout << ">>>---- Warning! Setting up CVODES solver twice! will delete";
    cout << " and re-init. ----<<<\n";
    if (abstol) {
#if defined CVODE2
      N_VDestroy_Serial(abstol);
#else
      N_VDestroy(abstol);
#endif
    }
    if (cvode_mem) {
      CVodeFree(&cvode_mem);
    }
  }

  int err = setup_cvode_solver_without_Jacobian();
  if (err) {
    cout << "Failed to setup solver without jacobian. err=" << err << "\n";
    return err;
  }
  //
  // Set the Jacobian routine to Jac (user-supplied)
  //
#if defined CVODE2
  err = CVDlsSetDenseJacFn(cvode_mem, Jacobian_for_cvode);
  if (err != CVDLS_SUCCESS) {
    cerr << "setup_cvode_solver() CVDlsSetDenseJacFn: err=" << err << "\n";
    return 6;
  }
#elif defined CVODE3
  err = CVDlsSetJacFn(cvode_mem, Jacobian_for_cvode);
  if (err != CVDLS_SUCCESS) {
    cerr << "setup_cvode_solver() CVDlsSetDenseJacFn: err=" << err << "\n";
    return 6;
  }
#else
  err = CVodeSetJacFn(cvode_mem, Jacobian_for_cvode);
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVDlsSetJacFn: err=" << err << "\n";
    return 6;
  }
#endif
  // All done now, so return.
  have_setup_cvodes = true;
  return 0;
}

// ##################################################################
// ##################################################################

int cvode_solver::setup_cvode_solver_without_Jacobian()
{
  if (have_setup_cvodes) {
    cout << ">>>---- Warning! Setting up CVODES solver twice! will delete";
    cout << " and re-init. ----<<<\n";
    if (abstol) {
      N_VDestroy_Serial(abstol);
      abstol = 0;
    }
    if (cvode_mem) {
      CVodeFree(&cvode_mem);
    }
  }

  int err = 0;

  get_problem_size(&n_eq, &n_xd);
#ifdef CVODE_DEBUG
  cout << "\t\tN_equations=" << n_eq << ", N_extra_data=" << n_xd << "\n";
#endif
  if (n_eq < 0 || n_xd < 0) {
    cout << "Error in values for problem size.\n";
    return 1;
  }

  //
  // initialise vector for abstol.
  //
  if (abstol) {
    cout << "Error! setup_cvode_solver(), vectors already initialised!\n";
    return 1;
  }

  abstol = N_VNew_Serial(n_eq);

  //
  // Call CVodeCreate to create the solver memory and specify the
  // Backward Differentiation Formula and the use of a Newton iteration.
  //
#if defined(CVODE2) || defined(CVODE3)
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#else
  cvode_mem = CVodeCreate(CV_BDF);
#endif
  if (!cvode_mem) {
    cerr << "setup_cvode_solver() error: cvode_mem=" << cvode_mem << "\n";
    return 2;
  }

  //
  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the inital time T0, and
  // the initial dependent variable vector y.
  //
  double t = 0.0;
  err      = CVodeInit(cvode_mem, Ydot_for_cvode, t, abstol);
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVodeInit error: " << err << "\n";
    return 3;
  }

  //
  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  //
  double reltol = 0.0;
  double *atol  = NV_DATA_S(abstol);
  get_error_tolerances(&reltol, atol);  // both args passed by reference.

#ifdef CVODE_DEBUG
  cout << "\t\treltol=" << reltol << ", atol=[" << NV_Ith_S(abstol, 0);
  for (int v = 1; v < n_eq; v++)
    cout << ", " << NV_Ith_S(abstol, v);
  cout << "]\n";
#endif

  err = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVodeSVtolerances: err=" << err << "\n";
    return 4;
  }

#if defined(CVODE2)

  // Call CVDense to specify the CVDENSE dense linear solver
#ifdef LAPACK
  err = CVLapackDense(cvode_mem, n_eq);
#else
  err = CVDense(cvode_mem, n_eq);
#endif
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVDense(): err=" << err << "\n";
    return 5;
  }

#elif defined(CVODE3)

  msetup = SUNDenseMatrix(n_eq, n_eq);
  vsetup = N_VNew_Serial(n_eq);
  LS     = SUNDenseLinearSolver(vsetup, msetup);
  err    = CVDlsSetLinearSolver(cvode_mem, LS, msetup);
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVDlsSetLinearSolver(): err=" << err << "\n";
    return 5;
  }

#else

  msetup = SUNDenseMatrix(n_eq, n_eq);
  vsetup = N_VNew_Serial(n_eq);
  LS     = SUNLinSol_Dense(vsetup, msetup);
  err    = CVodeSetLinearSolver(cvode_mem, LS, msetup);
  if (err != CV_SUCCESS) {
    cerr << "setup_cvode_solver() CVDlsSetLinearSolver(): err=" << err << "\n";
    return 5;
  }

#endif

  // Should be all done now, so return.
  have_setup_cvodes = true;
  return 0;
}

// ##################################################################
// ##################################################################

int cvode_solver::integrate_cvode_step(
    N_Vector Y_Input,
    void *user_data,  ///<  user_data
    double t_now,
    double dt,
    N_Vector Y_Output)
{
  int err = 0;
  //
  // user_data is a pointer to this class, but I can just use the 'this'
  // pointer instead.
  //
  err = CVodeSetUserData(cvode_mem, static_cast<void *>(this));
  if (err != CV_SUCCESS) {
    cerr << "integrate_cvode_step() CVodeSetUserData: err=" << err << "\n";
    return 3;
  }

  //
  // First try to do the integral in a single step and, if that fails, go into
  // a do-loop where we split the integral into substeps.
  //
  double t_temp = t_now;
  double tf     = t_now + dt;
  //
  // Re-initialise cvode to the current Y-value and start-time.
  //
  err = CVodeReInit(cvode_mem, t_now, Y_Input);
  if (err != CV_SUCCESS) {
    cerr << "integrate_cvode_step() CVodeReInit(): err=" << err << "\n";
    return 4;
  }
  //
  // integrate one timestep, returning answer to temp array Y_Output, and new
  // time to temporary variable t_temp.
  //
  err = CVode(cvode_mem, t_now + dt, Y_Output, &t_temp, CV_NORMAL);
  if (err == CV_SUCCESS) {
#ifdef CVODE_DEBUG
    t_now = t_temp;
    cout << "Success on first try: t=" << t_now << ", y = [";
    for (int v = 0; v < n_eq - 1; v++)
      cout << NV_Ith_S(Y_Output, v) << ", ";
    cout << NV_Ith_S(Y_Output, n_eq - 1) << "]\n";
#endif  // CVODE_DEBUG
    return 0;
  }

  //
  // If we get to here then the loop failed, so we reset the time and shorten
  // dt. Run a do-while loop to get to the end of the step.  We may need to do
  // it in steps, which is why it is in a loop.
  //
#ifdef CVODE_DEBUG
  cout << "First try failed\n";
#endif  // CVODE_DEBUG
  t_temp = t_now;
  dt *= 0.5;
  int fail_ct = 0, step_ct = 0;
  do {
    //
    // Re-initialise cvode to the current Y-value and start-time.
    //
    err = CVodeReInit(cvode_mem, t_now, Y_Input);
    if (err != CV_SUCCESS) {
      cerr << "integrate_cvode_step() CVodeReInit(): err=" << err << "\n";
      return 4;
    }
    //
    // integrate one timestep, returning answer to temp array Y_Output,
    // and new time to temporary variable t_temp.
    //
    err = CVode(cvode_mem, t_now + dt, Y_Output, &t_temp, CV_NORMAL);
    if (err != CV_SUCCESS) {
#ifdef CVODE_DEBUG
      cout << "FAILED STEP: err=" << err << ".";
      cout << "\tOld-t=" << t_now << ", returned t=" << t_temp
           << ", dt was = " << dt;
      cout << ".\tFAILED loop_ct=" << fail_ct << ": Y_Output = [";
      for (int v = 0; v < n_eq - 1; v++)
        cout << NV_Ith_S(Y_Output, v) << ", ";
      cout << NV_Ith_S(Y_Output, n_eq - 1) << "]\n";
      // cout <<"FAILED loop ct="<<fail_ct<<":    y = [";
      // for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(Y_Input,v)<<", ";
      // cout << NV_Ith_S(Y_Input,n_eq-1) <<"]\n";
#endif  // CVODE_DEBUG
        //
        // step failed, so don't copy yout to y, do shrink the timestep,
        // and don't update t.
        //
      dt *= 0.5;
      t_temp = t_now;
      fail_ct++;
    }
    else {
      //
      // Must have succeeded, copy Y_Output to y, make the next timestep
      // larger (if needed), and make sure dt is not too large,
      // or running past the end of the integration.
      //
#ifdef CVODE_DEBUG
      cout << "STEPPER: old-t=" << t_now << ", returned t=" << t_temp
           << ", dt was = " << dt;
      cout << ".\t This_loop_fail_ct=" << fail_ct << ", step_ct=" << step_ct
           << ": yout = [";
      for (int v = 0; v < n_eq - 1; v++)
        cout << NV_Ith_S(Y_Output, v) << ", ";
      cout << NV_Ith_S(Y_Output, n_eq - 1) << "]\n";
      // cout <<"EEEEEE loop ct="<<fail_ct<<":    y = [";
      // for (int v=0;v<n_eq-1;v++) cout << NV_Ith_S(Y_Input,v)<<", ";
      // cout << NV_Ith_S(Y_Input,n_eq-1) <<"]\n";
#endif  // CVODE_DEBUG

      for (int v = 0; v < n_eq; v++)
        NV_Ith_S(Y_Input, v) = NV_Ith_S(Y_Output, v);
      t_now = t_temp;
      dt *= 1.5;
      dt = min(dt, tf - t_now);
      //
      // step succeeded, so reset fail_ct to zero and increment step_ct
      //
      fail_ct = 0;
      step_ct++;
    }
  } while (t_now < tf && fail_ct <= MAX_TRYS && step_ct <= MAX_STEPS);

  if (fail_ct > MAX_TRYS) {
    cout << "integrate_cvode_step(): Integration failed after bisecting the "
            "step "
         << MAX_TRYS << "times.\n";
    return fail_ct;
  }
  if (step_ct >= MAX_STEPS) {
    cout << "integrate_cvode_step(): Integration took " << MAX_STEPS
         << " and not at end of step.  giving up.\n";
    return step_ct;
  }

#ifdef CVODE_DEBUG
  cout << "Final version step_ct=" << step_ct << ": t=" << t_now << ", y = [";
  for (int v = 0; v < n_eq - 1; v++)
    cout << NV_Ith_S(Y_Output, v) << ", ";
  cout << NV_Ith_S(Y_Output, n_eq - 1) << "]\n";
#endif  // CVODE_DEBUG

  // for (int v=0;v<n_eq;v++) NV_Ith_S(Yf,v) = NV_Ith_S(Y_Output,v);
  return 0;
}

// ##################################################################
// ##################################################################

// ------- NOT CLASS FUNCTIONS -- THESE ARE THE FUNCTIONS THAT CVODE
// ------- USES AS FUNCTION POINTERS.

// ##################################################################
// ##################################################################

int Ydot_for_cvode(
    double t,     ///< current time
    N_Vector y,   ///< current Y-value
    N_Vector yd,  ///< vector for Y-dot values
    void *data    ///< extra user-data pointer to the solver class
)
{
  //
  // Now call the class member ydot() function:
  //
  class cvode_solver *S = static_cast<cvode_solver *>(data);
  return S->ydot(t, y, yd, 0);
}

// ##################################################################
// ##################################################################

int Jacobian_for_cvode(
#if defined CVODE2
    long int N,  ///< N (not sure what this is for!)
#endif
    double t,       ///< time, t
    N_Vector y,     ///< y
    N_Vector yd,    ///< ydot
    CVMatrix J,     ///< Jacobian matrix
    void *data,     ///< extra user-data pointer to the solver class
    N_Vector tmp1,  ///< temp vector, must be for internal use
    N_Vector tmp2,  ///< temp vector, must be for internal use
    N_Vector tmp3   ///< temp vector, must be for internal use
)
{
  //
  // Now call class member function to get the Jacobian.
  // This is stored in J as a 'CVMatrix' stuct, so we don't
  // to post-process anything afterwards.
  //
  class cvode_solver *S = static_cast<cvode_solver *>(data);
#if defined CVODE2
  return S->Jacobian(N, t, y, yd, 0, J);
#else
  return S->Jacobian(t, y, yd, 0, J);
#endif
}

// ##################################################################
// ##################################################################
