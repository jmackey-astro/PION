#ifndef STELLAR_WIND_ANGLE_H
#define STELLAR_WIND_ANGLE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid/stellar_wind_BC.h"
#include "sim_constants.h"
#include "sim_params.h"
#include "tools/interpolate.h"


//
// Stellar wind class for angle dependent winds, developed for BSGs etc.
// Added by Robert Kavanagh (21/7/17)
//
class stellar_wind_angle :
    virtual public stellar_wind_evolution,
    virtual public interpolate_arrays {
public:
  ///
  /// Constructor:
  ///
  stellar_wind_angle(
      const int,            ///< ndim
      const int,            ///< nvar
      const int,            ///< ntracer
      const int,            ///< ftr
      const std::string *,  ///< List of tracer variable names.
      const int,            ///< coord_sys
      const int,            ///< eqn_type
      const double,         ///< minimum temperature allowed
      const double,         ///< Simulation start time.
      const double,         ///< Simulation finish time.
      const double          ///< xi value for equatorial flattening of wind
  );

  ///
  /// Destructor (Delete spline arrays, etc.)
  ///
  ~stellar_wind_angle();

  ///
  /// Add an evolving source for rotating star: takes in the filename of the
  /// source data, and a time offset between the start of the simulation and
  /// the time in the stellar model (may need to be <0 so that wind feedback
  /// starts immediately).
  ///
  virtual int add_evolving_source(
      const double,                ///< current time.
      struct stellarwind_params *  ///< pointer to wind parameters struct
  );

  ///
  /// Add a wind source, returns source id (count from zero).
  /// Note the temperature is in Kelvin if we have a pure neutral atomic
  /// hydrogen gas, otherwise it will be modified accordingly.
  /// - Note the v_crit should be input from the evolution file or
  ///   pre-calculated, and does not have to be consistent with v_esc.
  ///   The reason for this is that for rotating stars these
  ///   quantities are averaged over the surface (or surface layers!)
  ///
  virtual int add_rotating_source(
      struct stellarwind_params *  ///< pointer to wind parameters struct
  );

  /// setup tables for interpolation.
  void setup_tables();

  /// Wind terminal velocity function
  double fn_v_inf(
      double,  ///< omega (v_rot/v_crit)
      double,  ///< terminal velocity at pole (cm/s)
      double   ///< theta (latitude, measured from pole)
  );

  /// Analytic wind density function
  double fn_density(
      double,  ///< omega (v_rot/v_crit)
      double,  ///< terminal velocity at pole (cm/s)
      double,  ///< mass loss rate (g/s)
      double,  ///< radius (cm)
      double,  ///< theta (latitude, measured from pole)
      double   ///< Teff (K)
  );

  /// Interpolated wind density function
  double fn_density_interp(
      double,  ///< omega (v_rot/v_crit)
      double,  ///< terminal velocity at pole (cm/s)
      double,  ///< mass loss rate (g/s)
      double,  ///< radius
      double,  ///< theta
      double   ///< Teff
  );

private:
  /// Simpson's rule integration of function f
  double integrate_Simpson(
      const double,    ///< lower limit
      const double,    ///< upper limit
      const long int,  ///< Number of points (must be even)
      const double,    ///< omega
      const double     ///< Teff (K)
  );

  /// Phi' function
  double fn_phi(
      double,  ///< omega
      double,  ///< theta
      double   ///< Teff (K)
  );

  /// Alpha function
  double fn_alpha(
      double,  ///< omega
      double,  ///< theta
      double   ///< Teff (K)
  );

  /// Integrand to be integrated to get delta
  double integrand(
      double,       ///< theta
      double,       ///< omega
      const double  ///< Teff
  );

  /// Delta function
  double fn_delta(
      double,  ///< omega
      double   ///< Teff (K)
  );

  double c_gamma;  ///< exponent in velocity formula
  double c_xi;     ///< exponent in density formula
  double c_beta;   ///< exponent in factor added to mass loss - Mdot = Mdot*(1
                   ///< - omega)^c_beta

  int npts_theta;  ///< number of points in theta vector
  int npts_omega;  ///< number of points in omega vector
  int npts_Teff;   ///< number of points in Teff vector

  vector<double> theta_vec;   ///< theta vector
  vector<double> log_mu_vec;  ///< log(mu) = log(1 - omega) vector
  vector<double> omega_vec;   ///< omega vector
  vector<double> Teff_vec;    ///< Teff vector

  vector<vector<double> > delta_vec;           ///< delta table
  vector<vector<vector<double> > > alpha_vec;  ///< alpha table

protected:
  ///
  /// Set values of wind_cell reference state based on Wind-Source
  /// properties and the cell-to-source distance.
  ///
  void set_wind_cell_reference_state(
      class GridBaseClass &,       ///< grid class
      struct wind_cell &,          ///< cell to calculate for
      const struct wind_source &,  ///< wind source struct
      const double                 ///< EOS gamma
  );

  ///
  /// Update both the wind properties and the state vectors of all of
  /// the wind cells, for rotating star with latitude-dependent wind.
  ///
  virtual int update_source(
      class GridBaseClass *,
      struct evolving_wind_data *,  ///< source to update.
      const double,                 ///< current simulation time.
      const double                  ///< EOS Gamma
  );

  /// set values in cells.  This function decides if the wind is evolving or
  /// not and, if so, updates the source before setting cell values.
  int set_cell_values(
      class GridBaseClass *,
      const int,    ///< src id
      const double  ///< simulation time
  );
};
#endif
