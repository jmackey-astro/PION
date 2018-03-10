#ifndef STELLAR_WIND_ANGLE_H
#define STELLAR_WIND_ANGLE_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "sim_constants.h"
#include "sim_params.h"
#include "tools/reporting.h"
#include "grid/stellar_wind_BC.h"


//
// Stellar wind class for angle dependent winds, developed for BSGs etc.
// Added by Robert Kavanagh (21/7/17)
//
class stellar_wind_angle
: virtual public stellar_wind_evolution,
  virtual public interpolate_arrays  {
  public:
  ///
  /// Constructor: 
  ///
  stellar_wind_angle(
      const int, ///< ndim
      const int, ///< nvar
      const int, ///< ntracer
      const int, ///< ftr
      const int, ///< coord_sys
      const int  ///< eqn_type
      );

  ///
  /// Destructor (Delete spline arrays, etc.)
  ///
  ~stellar_wind_angle();

  ///
  /// Add an evolving source for rotating star: takes in the filename of the source data, and 
  /// a time offset between the start of the simulation and the time in the 
  /// stellar model (may need to be <0 so that wind feedback starts immediately).
  ///
  virtual int add_evolving_source(
      const double *, ///< position (physical units).
      const double,   ///< radius (physical units).
      const int,      ///< type (must be WINDTYPE_ANGLE).
      const double,   ///< Radius at which to get gas pressure from Teff
      const pion_flt *, ///< Any (constant) wind tracer values.
      const string,   ///< file name to read data from.
      const int,      ///< enhance mdot based on rotation (0=no,1=yes).
      const double,   ///< time offset = [t(sim)-t(wind_file)]
      const double,   ///< current time.
      const double,   ///< frequency with which to update wind properties.
      const double    ///< scale factor for time (t(sim)=[t(evo_file)-offset]/scalefactor
      );
  
  ///
  /// Add a wind source, returns source id (count from zero).
  /// Note the temperature is in Kelvin if we have a pure neutral atomic
  /// hydrogen gas, otherwise it will be modified accordingly.
  ///
  virtual int add_rotating_source(
      const double *, ///< position (cm from grid origin)
      const double,   ///< radius (cm)
      const int,      ///< type (2=lat-dep.)
      const double,   ///< Mdot (g/s)
      const double,   ///< Vesc (cm/s)
      const double,   ///< Vrot (cm/s)
      const double,   ///< Wind Temperature (p_g.m_p/(rho.k_b))
      const double,   ///< Radius where T=Twind (to get gas pressure)
      const pion_flt *  ///< Tracer values of wind (if any)
      );

  // Function to replace pow(a, b) - exp(b*log(a)) is twice as fast
  double pow_fast(
	double, ///< a
	double ///< b
	);

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
	double, ///< theta
	double ///< Teff(K)
	);

  //
  // Analytic wind density function
  //
  double fn_density(
	double, ///< omega
	double, ///< escape velocity (cm/s)
	double, ///< mass loss rate (g/s)
	double, ///< radius
	double, ///< theta
	double ///< Teff (K)
  );
    
  //
  // Interpolated wind density function
  //
  double fn_density_interp(
	double, ///< omega
	double, ///< escape velocity (cm/s)
	double, ///< mass loss rate (g/s)
	double, ///< radius
	double, ///< theta
	double ///< Teff
	);
  
  private:
  //
  // Eldridge et al. (2006, MN, 367, 186).
  // v_inf = sqrt(beta)*v_esc
  //
  double beta(const double); ///< Teff
  
  //
  // Simpson's rule integration of function f
  //
  double integrate_Simpson(
    const double, ///< lower limit
    const double, ///< upper limit
    const long int, ///< Number of points (must be even)
    const double, ///< omega
	const double ///< Teff (K)
    );

  //
  // Phi' function
  //
  double fn_phi(
	double, ///< omega
	double, ///< theta
	double ///< Teff (K)
	);

  //
  // Alpha function
  //
  double fn_alpha(
	double, ///< omega
	double, ///< theta
	double ///< Teff (K)
	);

  //
  // Integrand to be integrated to get delta
  //
  double integrand(
    double, ///< theta
    double,  ///< omega
	const double ///< Teff
    );
  
  //
  // Delta function
  //
  double fn_delta(
	double, ///< omega
	double ///< Teff (K)
	);

  double c_gamma; ///< exponent in velocity formula
  double c_xi;    ///< exponent in density formula
  double c_beta; ///< exponent in factor added to mass loss - Mdot = Mdot*(1 - omega)^c_beta 

  int npts_theta; ///< number of points in theta vector
  int npts_omega; ///< number of points in omega vector
  int npts_Teff;  ///< number of points in Teff vector
  
  vector<double> theta_vec;  ///< theta vector
  vector<double> log_mu_vec; ///< log(mu) = log(1 - omega) vector
  vector<double> omega_vec;  ///< omega vector
  vector<double> Teff_vec;   ///< Teff vector
      
  vector<vector<double> > delta_vec;          ///< delta table
  vector<vector<vector<double> > > alpha_vec; ///< alpha table
  
  vector<double> time_evo;
  vector<double> M_evo;
  vector<double> L_evo;
  vector<double> R_evo;
  vector<double> Teff_evo;
  vector<double> Mdot_evo;
  vector<double> vrot_evo;
  vector<double> vesc_evo;

  protected:

  ///
  /// Set values of wind_cell reference state based on Wind-Source properties
  /// and the cell-to-source distance.
  ///
  void set_wind_cell_reference_state(
      class GridBaseClass *,
      struct wind_cell *,
      const struct wind_source *
      );

  ///
  /// Update both the wind properties and the state vectors of all of the
  /// wind cells, for a rotating star with latitude-dependent wind.
  ///
  virtual void update_source(
      class GridBaseClass *,
      struct evolving_wind_data *, ///< source to update.
      const double ///< current simulation time.
      );
 
};


#endif
