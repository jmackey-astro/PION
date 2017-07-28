#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "tools/interpolate.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "grid/grid_base_class.h"
#include "grid/stellar_wind_angle.h"
#include "microphysics/microphysics_base.h"
#include <sstream>


#include "wind.h"
#include "constants.h"
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


// ##################################################################
// ##################################################################


// Constructor
stellar_wind_angle::stellar_wind_angle()
{
	// Constants for wind functions
	stellar_wind_angle::c_gamma = 0.35;
	stellar_wind_angle::c_zeta  = 1.0;
	stellar_wind_angle::c_xi    = -0.43;
	stellar_wind_angle::npts 	= 25; // change depending on tests

	setup_tables();
}


// ##################################################################
// ##################################################################


// Destructor
stellar_wind_angle::~stellar_wind_angle()
{
  // *****************************
  // Need to delete wind cell struct?
  // *****************************
}


// ##################################################################
// ##################################################################


// Generate interpolating tables for wind density function
void stellar_wind_angle::setup_tables()
{   
    //
	// Set up theta array
    //
    
    // Min, mid and max values of theta - 5pts between min and mid, 20pts between mid and max
    double theta_min = 0.1;
    double theta_mid = 60.0;
    double theta_max = 89.9;
    
	theta_vec.resize(npts);
    
    for (int k = 0; k < npts; k++)
    {
        if (k <= 4) theta_vec[k] = (theta_min + k*((theta_mid - theta_min)/4.0))*(pconst.pi()/180.0);
        else theta_vec[k] = (theta_mid + (k - 4)*((theta_max - theta_mid)/(npts - 5)))*(pconst.pi()/180.0);
	}
    
	//
    // Set up omega array
    //
    
	log_mu_vec.resize(npts);    

    // Iterate log(mu) values - evenly spaced (fixed s.t. smallest element of omega_vec is first)
    for (int i = 0; i < npts; i++) log_mu_vec[npts - i - 1] = -4.0 + i*(4.0/npts);

	omega_vec.resize(npts);
    
    // Iterate omega values - log spacing - spacing decreases as omega approaches 1
    for (int j = 0; j < npts; j++) omega_vec[j] = 1 - pow(10, log_mu_vec[j]);

    //
    // Write delta table
    //

	delta_vec.resize(npts);

	for (int i = 0; i < npts; i++) delta_vec[i] = fn_delta(omega_vec[i]);

    //
    // Write alpha table
    //
	
	alpha_vec.resize(npts);
	for (int i = 0; i < npts; i++) alpha_vec[i].resize(npts);

	for (int x = 0; x < npts; x++){
		for (int y = 0; y < npts; y++){
			alpha_vec[y][x] = fn_alpha(omega_vec[x], theta_vec[y]);
		}
	}
  return;
}


// ##################################################################
// ##################################################################


// Integrand for integral in delta function
double stellar_wind_angle::integrand(
	double theta, // Co-latitude angle (radians)
	double omega // Omega (v_rot/v_esc)
    )
{
	return fn_alpha(omega, theta) * pow(1.0 - omega*sin(theta), c_xi) * sin(theta);
} 


// ##################################################################
// ##################################################################


// Simpson's rule integration scheme
double  stellar_wind_angle::integrate_Simpson(
    const double min, ///< lower limit
    const double max, ///< upper limit
    const long int npt,    ///< Number of points (must be even)
    const double omega ///<  omega
    )
{

  //
  // Simpson's rule  int=(f(0)+4f(1)+2f(2)+...+2f(n-2)+4f(n-1)+f(n))*h/3
  // where n is even.
  //
  double hh = (max-min)/npt;
  double ans = 0.0;
  //
  // f(0) lower limit
  //
  ans += integrand(min, omega);
  //
  // f(N) upper limit
  //
  ans += integrand(max, omega);
  //
  // Intermediate points.
  //
  int wt = 4; double x=0.0;
  for (long int i=1; i<npt; i++) {
    x = min + i*hh;
    ans += wt*integrand(x, omega);
    wt = 6-wt;
  }
  //
  // finally multiply by hh/3.0
  //
  ans *= hh/3.0;
  return ans;
} 


// ##################################################################
// ##################################################################


// Phi' function
double stellar_wind_angle::fn_phi(
	double omega, // omega (v_rot/v_esc)
	double theta // Co-latitude angle (radians)
	)
{
	return (omega/(22.0*pconst.sqrt2()*c_zeta)) * sin(theta) * pow(1.0 - omega*sin(theta), -c_gamma);
}


// ##################################################################
// ##################################################################


// Alpha function
double stellar_wind_angle::fn_alpha(
	double omega, // Omega (v_rot/v_esc)
	double theta // Co-latitude angle (radians)
	)
{
	return pow(cos(fn_phi(omega, theta)) + pow(tan(theta),-2.0) * 
		   (1.0 + c_gamma*( omega*sin(theta) / (1.0 - omega*sin(theta)) )) * fn_phi(omega, theta) *
		   sin(fn_phi(omega, theta)), -1.0);
} // the cotan term will diverge here if theta = 0.0


// ##################################################################
// ##################################################################


// Delta function
double stellar_wind_angle::fn_delta(
	double omega // Omega (v_rot/v_esc)
	)
{
	return 2.0*pow(integrate_Simpson(0.001, pconst.pi()/2.0, 230, omega), -1.0);
}


// ##################################################################
// ##################################################################


// Terminal wind velocity function
double stellar_wind_angle::fn_v_inf(
	double omega, // Omega (v_rot/v_esc)
	double v_esc, // Escape velocity (cm/s)
	double theta // Co-latitude angle (radians)
	)
{
	return c_zeta * v_esc * pow(1.0 - omega*sin(theta), c_gamma);
}


// ##################################################################
// ##################################################################


// Wind density function (analytic - ~68 times slower than interpolated wind density function)
double stellar_wind_angle::fn_density(
	double omega, // Omega (v_rot/v_esc)
	double v_esc, // Escape velocity (cm/s)
	double mdot, // Mass loss rate (g/s)
	double radius, // Radius (cm)
	double theta // Co-latitude angle (radians)
    )
{
    return (mdot * fn_alpha(omega, theta) * fn_delta(omega) * pow(1.0 - omega*sin(theta), c_xi)) /
           (8.0 * pconst.pi() * pow(radius, 2.0) * fn_v_inf(omega, v_esc, theta));
}


// ##################################################################
// ##################################################################


// Interpolated wind density function
double stellar_wind_angle::fn_density_interp(
	double omega, // Omega (v_rot/v_esc)
	double v_esc, // Escape velocity (cm/s)
	double mdot, // Mass loss rate (g/s)
	double radius, // Radius (cm)
	double theta // Co-latitude angle (radians)
    )
{
    //
    // Use tables to interpolate the value of delta
    //

    double delta_interp;
    vector<double> d2_omega;
    d2_omega.resize(npts);
    for (int i = 0; i < npts; i++) d2_omega[i] = 0.0;	

    // Calculate second derivatives of omega wrt delta - stored in d2_omega
    spline_vec(omega_vec, delta_vec, npts, 1e32, 1e32, d2_omega);

    // Returns interpolated value of delta from input omega - stored in delta_interp
    splint_vec(omega_vec, delta_vec, d2_omega, npts, omega, delta_interp);

    //
    // Use tables to interpolate the value of alpha
    //

    double alpha_interp;

    vector<double> seek (2); // (omega, theta) input
    seek[0] = omega;
    seek[1] = theta;
    
    vector<size_t> nxy (2);
    nxy[0] = npts;
    nxy[1] = npts; // should we change this to accomodate for different sizes?

    alpha_interp = root_find_bilinear_vec(omega_vec, theta_vec, alpha_vec, nxy, seek);

    return (mdot * alpha_interp * delta_interp * pow(1.0 - omega*sin(theta), c_xi)) /
    (8.0 * pconst.pi() * pow(radius, 2.0) * fn_v_inf(omega, v_esc, theta));
}


// ##################################################################
// ##################################################################


void stellar_wind_angle::set_wind_cell_reference_state(
        class GridBaseClass *grid,
        struct wind_cell *wc,
        const struct wind_source *WS
        )
{
  //
  // In this function we set the density, pressure, velocity, and tracer values
  // for the reference state of the cell.  Every timestep the cell-values will
  // be reset to this reference state.
  //
  double pp[SimPM.ndim];
  CI.get_dpos(wc->c,pp);
  //rep.printVec("cell pos", pp, SimPM.ndim);
  //cout <<"dist="<<wc->dist<<"\n";

    //
    // 3D geometry, so either 3D-cartesian, 2D-axisymmetry, or 1D-spherical.
    //
    wc->p[RO] = fn_density_interp(pconst.sqrt2()*WS->v_rot/WS->v_esc, WS->v_esc, WS->Mdot, WS->radius, wc->theta);

    //
    // Set pressure based on wind density/temperature at the stellar radius,
    // assuming adiabatic expansion outside Rstar, and that we don't care what
    // the temperature is inside Rstar (because this function will make it 
    // hotter than Teff, which is not realistic):
    //
    // rho_star = Mdot/(4.pi.R_star^2.v_inf),
    //   p_star = rho_star.k.T_star/(mu.m_p)
    // So then p(r) = p_star (rho(r)/rho_star)^gamma
    // ******************************************************************************
    // *********** WARNING MU=1 HERE, PROBABLY SHOULD BE O.6 (IONISED) 1.3 (NEUTRAL).
    // ******************************************************************************
    //

	// **********************************
    // Is this valid for angle dependent wind?
    // **********************************

    wc->p[PG] = pconst.kB()*WS->Tw/pconst.m_p();
    wc->p[PG]*= exp((SimPM.gamma-1.0)*log(4.0*pconst.pi()*WS->Rstar*WS->Rstar*WS->Vinf/WS->Mdot));
    wc->p[PG]*= exp((SimPM.gamma)*log(wc->p[RO]));

  //
  // NOW VELOCITIES: These should be cell-average values, which are
  // the values at the centre-of-volume, so we call the geometry-aware
  // grid functions.
  //

  // calculate terminal wind velocity
  double Vinf = fn_v_inf(pconst.sqrt2()*WS->v_rot/WS->v_esc, WS->v_esc, wc->theta)

  cell *c = wc->c;

  // Velocity x component
  wc->p[VX] = Vinf*grid->difference_vertex2cell(WS->dpos,c,XX)/wc->dist;

  // Velocity y component (2D or 3D)
  if (SimPM.ndim>1)
    wc->p[VY] = Vinf*grid->difference_vertex2cell(WS->dpos,c,YY)/wc->dist;
  else
    wc->p[VY] = 0.0;

  // Velocity z component (3D)
  if (SimPM.ndim>2)
    wc->p[VZ] = Vinf*grid->difference_vertex2cell(WS->dpos,c,ZZ)/wc->dist;
  else
    wc->p[VZ] = 0.0;

  if (SimPM.eqntype!=EQEUL && SimPM.eqntype!=EQEUL_EINT)
    rep.error("Need to code B into winds model!",SimPM.eqntype);


  // update tracers
  for (int v=0;v<SimPM.ntracer;v++)
    wc->p[SimPM.ftr+v] = WS->tracers[v];
  //
  // HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK! HACK!
  // Assume the first tracer variable is the H+ ion fraction, and set it so
  // that it goes from y=1 at T>12500K to y=1.0e-7 at T<10000, with linear
  // interpolation.
  //
#ifdef HACK_WARNING
#error "REMOVE HACK setting ion fraction of H+ in winds"
#endif
  if (SimPM.ntracer>0) {
    if      (WS->Tw >1.25e4)
      wc->p[SimPM.ftr] = 1.0;
    else if (WS->Tw <1.00e4)
      wc->p[SimPM.ftr] = 1.0e-7;
    else
      wc->p[SimPM.ftr] = std::max((WS->Tw-1.0e4)*4e-4,1.0e-7);
  }

//#define TESTING
#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  //
  // Set the minimum temperature to be 10K in the wind...
  //
  if (MP) {
    if (MP->Temperature(wc->p,SimPM.gamma) <SimPM.EP.MinTemperature) {
      MP->Set_Temp(wc->p,SimPM.EP.MinTemperature,SimPM.gamma);
    }
  }
  else {
    // appropriate for a neutral medium.
    wc->p[PG] = max(static_cast<double>(wc->p[PG]), 
      SimPM.EP.MinTemperature*wc->p[RO]*pconst.kB()*(1.0-0.75*SimPM.EP.Helium_MassFrac)/pconst.m_p());
  }
#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
#ifdef TESTING
  cout << "\n";
#endif

  return;
}

