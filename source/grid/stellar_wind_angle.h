#ifndef STELLAR_WIND_ANGLE_H
#define STELLAR_WIND_ANGLE_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"

#include "sim_constants.h"
#include "sim_params.h"
#include "tools/reporting.h"

#include <stellar_wind_BC.h>


//
// Stellar wind class for angle dependent winds, developed for BSGs etc.
// Added by Robert Kavanagh (21/7/17)
//
class stellar_wind_angle : virtual public stellar_wind {
	public:
	///
	/// Constructor: 
	///
	stellar_wind_angle();

	///
	/// Destructor (Delete spline arrays, etc.)
	///
	~stellar_wind_angle();


  //
  // setup tables for interpolation.
  //
  void setup_tables();

  //
  // Wind terminal velocity function
  //
  double fn_v_inf(
	double, ///< omega
	double, ///< escape velocity (cm/s)
	double ///< theta
	);

  //
  // Wind density function
  //
  double fn_density(
	double, ///< omega
	double, ///< escape velocity (cm/s)
	double, ///< mass loss rate (g/s)
	double, ///< radius
	double ///< theta
  );
    
  //
  // Interpolated wind density function
  //
  double fn_density_interp(
	double, ///< omega
	double, ///< escape velocity (cm/s)
	double, ///< mass loss rate (g/s)
	double, ///< radius
	double ///< theta
	);
  
  private:
  //
  // Simpson's rule integration of function f
  //
  double integrate_Simpson(
    const double, ///< lower limit
    const double, ///< upper limit
    const long int, ///< Number of points (must be even)
    const double ///< omega
    );

  //
  // Phi' function
  //
  double fn_phi(
	double, ///< omega
	double ///< theta
	);

  //
  // Alpha function
  //
  double fn_alpha(
	double, ///< omega
	double ///< theta
	);

  //
  // Integrand to be integrated to get delta
  //
  double integrand(
    double, ///< theta
    double  ///< omega
    );
  
  //
  // Delta function
  //
  double fn_delta(
	double ///< omega
	);

  double c_gamma; ///< exponent in velocity formula
  double c_xi;    ///< exponent in density formula
  double c_zeta;  ///< scaling for wind velocity

  int npts; ///< number of points in each table vector
  vector<double> theta_vec; ///< theta table
  vector<double> log_mu_vec; ///< log(mu) = log(1 - omega) table
  vector<double> omega_vec; ///< omega table
  vector<double> delta_vec; ///< delta table
  vector< vector<double> > alpha_vec; ///< alpha table
  
};


#endif
