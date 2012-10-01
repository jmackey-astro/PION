/** \file cooling.h
 * 
 * \brief Class Declaration for cooling functions
 * \author Jonathan Mackey
 * 
 * Modifications:
 *  - 2000-01-21 Moved Cooling out from microphysics.h to its own file.
 * */
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)

#ifndef COOLING_H
#define COOLING_H


//#define COOL_TESTING ///< output cooling rate to a file.
#define LOGSLP
//#define LINSLP ///< Use linear extrapolation of cooling function for low temperatures.

/** \brief Class to tabulate and interpolate various cooling functions, with a public 
 * function to return the rate for any input temperature.
 * */
class CoolingFn {
  public:
   /** \brief Constructor.  Integer flag determines what cooling curve(s) to use.
    * Currently I only have the Sutherland \& Dopita (1993) Function (flag==1).
    * */
   CoolingFn(int ///< integer flag, 1=DopitaSutherland
	     );
   ~CoolingFn(); ///< destructor.
   /** \brief Returns the cooling rate for a given input State Vector.
    * 
    * Input Temperature should be in Kelvin, number density in per c.c., and ion
    * fraction between 0 and 1.  The cooling rate is returned in cgs units,
    * (rate = erg cm^-3 s^-1).  This should be the actual cooling rate, so no need to
    * multiply it by number density or electron/ion densities.
    */
   double CoolingRate(const double, ///< Input Temperature.
		      const double, ///< Hydrogen ion fraction [0,1]
		      const double, ///< nucleon number density, n_H
		      const double, ///< FUV flux G0 from eq.A3 Henney et al.2009.
		      const double  ///< FUV_extinction Av
		      );
  private:
   int WhichFunction; ///< Which Cooling Function we are using. (1=DS93)
   double MinTemp;    ///< Minimum Temperature we have cooling rate for.
   double MaxTemp;    ///< Maximum Temperature we have cooling rate for.
   double MinSlope;   ///< Slope of the Cooling Rate at the minimum temperature.
   double *Temp;      ///< Temperature data.
   double *Lamb;      ///< Normalised Cooling rate (erg cm^3 s^-1) = Lambda(T)/rho^2
   double *Lam2;      ///< Second derivative of the Cooling Rate for Cubic Spline Interpolation.
   int Nspl; ///< Number of points to use for the spline/splint data functions.
   double kB; ///< boltzmann constant.
};

#endif // COOLING_H
