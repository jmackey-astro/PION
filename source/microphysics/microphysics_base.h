///
/// \file microphysics_base.h
/// \author Jonathan Mackey
/// \date 2011.01.14
///
/// Microphysics base class, moved here from global.h.  This file is 
/// currently sourced by global.h.  I need to sort that out
/// eventually.
///
/// - 2011.02.25 JM: removed NEW_RT_MP_INTERFACE ifdef (it is assumed now)
/// - 2011.03.16 JM: Added TimeUpdateMP_RTnew() function to integrate the rate equations
///     NOT in the C2Ray way.  i.e. I assume the timestep is short enough that the
///     instantaneous optical depth is a good enough approximation to the time-averaged
///     optical depth.
/// - 2011.03.21 JM: Updated  RTnew() interface for more sources.  It is now simpler.
/// - 2011.05.02 JM: Added set_multifreq_source_properties() function
/// - 2012.01.16 JM: Tabbing.

#ifndef MICROPHYSICS_BASE_H
#define MICROPHYSICS_BASE_H

#include <vector>
#include "../grid/cell_interface.h"
//#define MP_DEBUG

///
/// Radiation Source data struct, used for passing info to microphysics classes.
///
struct rt_source_data {
  int id;   ///< source id.
  int type; ///< diffuse-radiation or a real source.
  double strength; ///< Luminosity (or flux if source at infinity).
  double Vshell;   ///< Shell volume for discrete photo-ionisation/-heating rates.
  double dS;       ///< Path length through cell.
  double Column;   ///< integral of density along line of sight to far edge of cell.
  double DelCol;   ///< integral of density along line of sight through cell.
};




/** \brief pure virtual base/interface class for in-cell microphysics update.*/
class MicroPhysicsBase {
  public :
   virtual ~MicroPhysicsBase() {} ///< non-virtual destructor.

   /** \brief Non-RT microphysics update, so cooling and heating and chemistry.
    * 
    * This uses various integration methods to update the elements and the
    * internal energy.
    * */
   virtual int TimeUpdateMP(const double *, ///< Primitive Vector to be updated.
			    double *,       ///< Destination Vector for updated values.
			    const double,   ///< Time Step to advance by.
			    const double,   ///< EOS gamma.
			    const int,      ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
			    double *        ///< final temperature.
			    )=0;

   ///
   /// If doing ray-tracing, the tracer can call this function to
   /// integrate the microphysics variables forward one timestep
   /// given an external radiation flux.
   ///
   virtual int TimeUpdate_RTsinglesrc(
              const double *, ///< Primitive Vector to be updated.
              double *,       ///< Destination Vector for updated values.
              const double,   ///< Time Step to advance by.
              const double,   ///< EOS gamma.
              const int,      ///< Switch for what type of integration to use. (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
              const double,   ///< flux in per unit length along ray (F/ds or L/dV)
              const double,   ///< path length ds through cell.
              const double,   ///< Optical depth to entry point of ray into cell.
              double *        ///< return optical depth through cell in this variable.
              )=0;

   /** \brief Initialise microphysics ionisation fractions to an equilibrium value. */
   virtual int Init_ionfractions(double *, ///< Primitive vector to be updated.
				 const double, ///< eos gamma.
				 const double  ///< optional gas temperature to end up at. (negative means use pressure)
				 )=0;

  ///
  /// Given a cell, calculated converged time-averaged rates of 
  /// change for species abundances (x(E)) and heating/cooling, including
  /// radiative effects.
  /// 
  /// The rates are to be stored in the Ph vector of the cell.  This
  /// means that this function must not be called on a half-timestep, and
  /// is supposed to be completely operator split from the dynamics.  It
  /// should only be called by the raytracing class.
  /// 
  /// The idea is that for a infinitesimally weak radiation source, this
  /// will produce an identical update to the non-radiative TimeUpdateMP
  /// function, once the time-averaged rates are added to the primitive 
  /// vector to get the updated value.
  /// THIS IS AN UNIMPLEMENTED FEATURE, WHICH WOULD ALLOW ME TO CALCULATE
  /// THE C2RAY METHOD UPDATES FOR MANY SOURCES.
  ///
  virtual int Get_MP_Rates(
                      cell *,       ///< current cell.
                      const double, ///< timestep
                      const double, ///< column to cell
                      const double, ///< Flux entering cell?
                      double *      ///< initial (and final) column density through cell?
                      ) {return 0;}

  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative transfer,
  /// where the column densities for diffuse and direct radiation are included as 
  /// parameters.  The input list of column densities is ordered in a list of UV
  /// heating sources and a list of ionising sources.
  ///
  virtual int TimeUpdateMP_RTnew(
                    const double *, ///< Primitive Vector to be updated.
                    const int,      ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of UV-heating column densities and source properties.
                    const int,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of ionising src column densities and source properties.
                    double *,       ///< Destination Vector for updated values
                    ///< (can be same as first Vector.
                    const double,   ///< Time Step to advance by.
                    const double,   ///< EOS gamma.
                    const int, ///< Switch for what type of integration to use.
                    ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
                    double *    ///< final temperature (not strictly needed).
                    ) {return 1;}

   /** \brief Given a string, return the index of the tracer that it refers to.
    * 
    * \retval int>0 success: int=tracer value
    * \retval int<0 failure: string doesn't correspond to any tracer.
    * */
   virtual int Tr(const string ///< name of tracer we are looking for.
		  )=0;

   /** \brief Sets the pressure in p-vec to be such that temperature is
    * what you want it to be. */
   virtual int Set_Temp(double *, ///< primitive vector.
			const double, ///< temperature.
			const double ///< eos gamma.
			)=0;

   /** \brief Returns the gas temperature (not very optimized though) 
    *
    * Assumes primitive vector is in cgs units.
    */
   virtual double Temperature(const double *, ///< primitive vector
			      const double    ///< eos gamma
			      )=0;

   ///
   /// This returns the minimum timescale of the times flagged in the
   /// arguments.  Time is returned in seconds.
   ///
   virtual double timescales(const double *, ///< Current cell.
			     const double,   ///< EOS gamma.
			     const bool, ///< set to true if including cooling time.
			     const bool, ///< set to true if including recombination time.
			     const bool  ///< set to true if including photo-ionsation time.
			     )=0;

  ///
  /// This returns the minimum timescale of all microphysical processes, including
  /// reaction times for each species and the total heating/cooling time for the gas.
  /// It requires the radiation field as an input, so it has substantially greater
  /// capability than the other timescales function.
  ///
  virtual double timescales_RT(
                    const double *, ///< Current cell.
                    const int,      ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of UV-heating column densities and source properties.
                    const int,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of ionising src column densities and source properties.
                    const double   ///< EOS gamma.
                    )=0;

  ///
  /// Set the properties of a multifrequency ionising radiation source.
  ///
  virtual int set_multifreq_source_properties(
                const struct rad_src_info *
                ) {return -999;}

};

#endif // MICROPHYSICS_BASE_H
