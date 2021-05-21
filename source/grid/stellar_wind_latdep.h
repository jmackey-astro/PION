#ifndef STELLAR_WIND_LATDEP_H
#define STELLAR_WIND_LATDEP_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid/stellar_wind_BC.h"
#include "sim_constants.h"
#include "sim_params.h"
#include "tools/interpolate.h"
#include "tools/reporting.h"


///
/// Stellar wind class for angle dependent winds, developed for BSGs etc.
/// Modified version of the Langer et al. (1999) prescription that
/// ensures that the polar wind is not enhanced by rotation.
///
class stellar_wind_latdep :
    virtual public stellar_wind_evolution,
    virtual public interpolate_arrays {
public:
  ///
  /// Constructor:
  ///
  stellar_wind_latdep(
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
      const double          ///< exponent of wind equatorial enhancement
  );

  ///
  /// Destructor (Delete spline arrays, etc.)
  ///
  ~stellar_wind_latdep();

  ///
  /// Add an evolving source for rotating star: takes in the filename
  /// of the source data, and a time offset between the start of the
  /// simulation and the time in the stellar model.
  ///
  virtual int add_evolving_source(
      const double *,  ///< position (physical units).
      const double,    ///< radius (physical units).
      const int,       ///< type (must be WINDTYPE_LATDEP).
      pion_flt *,      ///< Any (constant) wind tracer values.
      const string,    ///< file name to read data from.
      const int,       ///< enhance mdot based on rotation (0=no,1=yes).
      const double,    ///< Surface B field (G)
      const double,    ///< time offset = [t(sim)-t(wind_file)]
      const double,    ///< current time.
      const double,    ///< frequency with which to update wind properties.
      const double     ///< scale factor for time
                       ///< (t(sim)=[t(evo_file)-offset]/scalefactor
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
      const double *,  ///< position (cm from grid origin)
      const double,    ///< radius (cm)
      const int,       ///< type (2=lat-dep.)
      const double,    ///< Mdot (g/s)
      const double,    ///< Md0, equiv. non-rotating star (g/s)
      const double,    ///< Vinf (cm/s)
      const double,    ///< Vrot (cm/s)
      const double,    ///< Vcrit (cm/s)
      const double,    ///< Wind Temperature (p_g.m_p/(rho.k_b))
      const double,    ///< Radius of star (cm)
      const double,    ///< Surface B-field (Gauss)
      pion_flt *       ///< Tracer values of wind (if any)
  );

  /// setup tables for interpolation.
  void setup_tables();

  /// Wind terminal velocity function
  double fn_v_inf(
      double,  ///< theta (latitude, measured from pole)
      double,  ///< Omega (v_rot/v_crit)
      double   ///< terminal velocity at pole (cm/s)
  );

  /// Wind terminal velocity function (interpolated from table)
  double interp_v_inf(
      double,  ///< theta (latitude, measured from pole)
      double,  ///< Omega (v_rot/v_crit)
      double   ///< terminal velocity at pole (cm/s)
  );


  /// Interpolated wind density function
  double interp_density(
      double,  ///< radius (cm)
      double,  ///< theta (latitude, measured from pole)
      double,  ///< terminal velocity at co-latitude theta (cm/s)
      double,  ///< Omega (v_rot/v_crit)
      double,  ///< actual mass loss rate (g/s)
      double   ///< mass-loss rate for non-rotating star (g/s).
  );

private:
  /// Simpson's rule integration of function f(t,O)
  double integrate_Simpson(
      const double,    ///< lower limit
      const double,    ///< upper limit
      const long int,  ///< Number of points (must be even)
      const double     ///< Omega
  );

  /// Integrand to be integrated to get normalisation A
  double f(
      double,  ///< theta
      double   ///< Omega
  );

  double c_gamma;  ///< exponent in velocity formula
  double c_xi;     ///< exponent in density formula

  int Ntheta;  ///< number of points in theta vector
  int NOmega;  ///< number of points in Omega vector

  vector<double> theta_vec;   ///< theta vector
  vector<double> log_mu_vec;  ///< log(mu) = log(1 - Omega) vector
  vector<double> Omega_vec;   ///< Omega vector

  vector<vector<double> > vinf_vec;  ///< vinf interpolation table
  vector<size_t> vinf_vsize;
  vector<vector<double> > f_vec;  ///< density interpolation table
  vector<size_t> f_vsize;
  vector<double> norm_vec;  ///< normalising func interpolation table

protected:
  ///
  /// Set values of wind_cell reference state based on Wind-Source
  /// properties and the cell-to-source distance.
  ///
  void set_wind_cell_reference_state(
      class GridBaseClass *,
      struct wind_cell *,
      const struct wind_source *,
      const double  ///< EOS gamma
  );

  ///
  /// Update both the wind properties and the state vectors of all of
  /// the wind cells, for rotating star with latitude-dependent wind.
  ///
  virtual void update_source(
      class GridBaseClass *,
      struct evolving_wind_data *,  ///< source to update.
      const double,                 ///< current simulation time.
      const double                  ///< EOS Gamma
  );
};


#endif
