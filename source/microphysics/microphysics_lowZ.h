///
/// \file microphysics_lowZ.h
/// \author 
/// \date 2010.10.12
///
/// Microphysics class for low metallicity gas at high redshift.
///
/// Modifications:
///
/// - 2010.10.12 JM: Created template class.
///
/// - 2011.01.14 JM: moved to microphysics/ sub-dir.
/// - 2011.03.16 JM: updated to get ready for HD's code to be included.
/// - 2011.03.21 JM: Updated  RTnew() interface for more sources.  It is now simpler.
/// - 2011.11.09 JM: Added Set_Temp() function for correcting negative pressure
///                  in dynamics solver.
/// - 2012.01.19 JM: made timescales() function work.
/// - 2012.03.12 JM: Added function to deal with multiple column densities for
///    ionising, UV-heating, dissociating radiation, and diffuse UV background.



#ifndef MICROPHYSICS_LOWZ
#define MICROPHYSICS_LOWZ

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#ifndef EXCLUDE_HD_MODULE


#include "microphysics/microphysics_base.h"

#include "microphysics/cooling.h"


//#ifdef INCLUDE_HARPREETS_MODULE

#include <vector>
#include "contrib/HD_MetalFree.h"

//#define MP_DEBUG  ///< Use this flag to wrap debugging info and checks


//-------------------------------------
// This class inherits from one base-class as follows:
// MicrophysicsBase is the base interface class in microphysics_base.h, all
// functions in that class must at least be defined here, even if they
// do nothing.
//
// Terminology:
// "primitive variables" are [density, pressure, vx,vy,vz, [Bx,By,Bz], tracers 0,1,...]
// "local variables" are internal to microphysics and can be whatever you want.
// 
// Functions declared here must be defined either here or in microphysics_lowZ.cc,
// but must have only one definition.
//------------------------------------- 

///
/// Microphysics class for low and/or zero metallicity gas, including
/// molecular chemistry and maybe dust [EXPAND ME WHEN DECIDED!]
///
class microphysics_lowz
:
  public MicroPhysicsBase,
  public solver
{
  //
  // Public functions are the interface to the class.
  //
 public:
  ///
  /// Constructor -- sets up data, allocates memory if needed.
  /// (Feel free to add more parameters; these are a minimal set).
  ///
  microphysics_lowz(const int,          ///< Total number of variables in state vector
		    const int,          ///< Number of tracer variables in state vector.
		    const std::string &, ///< List of what the tracer variables mean.
		    struct which_physics * ///< pointer to "which-physics" flags.
		    );
  ///
  /// Destructor -- frees any memory that was dynamically allocated.
  ///
  ~microphysics_lowz();

  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics (no Radiative transfer!).
  ///
  int TimeUpdateMP(const double *, ///< Primitive Vector to be updated.
		   double *,       ///< Destination Vector for updated values
		                   ///< (can be same as first Vector.
		   const double,   ///< Time Step to advance by.
		   const double,   ///< EOS gamma.
		   const int, ///< Switch for what type of integration to use.
		              ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
		   double *    ///< final temperature (not strictly needed).
		   );

  ///
  /// This takes a copy of the primitive vector and advances it in time over
  /// the step requested, and at the end copies the updated vector into the
  /// destination vector.  For fully local microphysics but WITH radiative transfer,
  /// where the column densities for diffuse and direct radiation are included as 
  /// parameters.  The input list of column densities is ordered by the number of 
  /// sources in each category in the vector of integers.
  ///
   ///
  int TimeUpdateMP_RTnew(
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
		   );


  ///
  /// Not used b/c we have no photoionisation
  ///
  int TimeUpdate_RTsinglesrc(const double *, ///< Primitive Vector to be updated.
			     double *,       ///< Destination Vector for updated values.
			     const double,   ///< Time Step to advance by.
			     const double,   ///< EOS gamma.
			     const int, ///< Switch for what type of integration to use.
			                ///< (0=adaptive RK5, 1=adaptive Euler,2=onestep o4-RK)
			     const double,   ///< flux in per unit length along ray (F/ds or L/dV)
			     const double,   ///< path length ds through cell.
			     const double,   ///< Optical depth to entry point of ray into cell.
			     double *        ///< return optical depth through cell here.
			     );

  ///
  /// Returns the gas temperature.  This is only needed for data output, so
  /// there is no need to make it highly optimized.
  ///
  double Temperature(const pion_flt *, ///< primitive vector
		     const double    ///< eos gamma
		     );

  ///
  /// This returns the minimum timescale of the times flagged in the
  /// arguments.  Time is returned in seconds.  Sometimes needed if the
  /// cooling/heating times are much shorter than the gas-dynamics 
  /// timestep and we need to limit the global timestep (this is a runtime
  /// option which is not really well integrated into the code yet).
  ///
  double timescales(const double *, ///< Current cell.
		    const double,   ///< EOS gamma.
		    const bool, ///< set to 'true' if including cooling time.
		    const bool, ///< set to 'true' if including recombination time.
		    const bool  ///< set to 'true' if including photo-ionsation time.
		    );

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
                    );

  ///
  /// Returns element number of named tracer variable in state vector.
  /// This is not important, so I define it here.  If you want you can
  /// redefine it in microphysics_lowz.cc
  ///
  int Tr(string ///< tracer we want to get index for;
	 ) {return -1;}

  ///
  /// Initialise microphysics ionisation fractions to an equilibrium value.
  /// Only needed if you want this feature in the initial condition generator.
  /// This will call Harpreet's Yinit() function to set initial values for tracers.
  ///
  int Init_ionfractions(double *, ///< Primitive vector to be updated.
			const double, ///< eos gamma.
			const double  ///< optional gas temperature to end up at. (-ve means use pressure)
			);

  ///
  /// Set the gas temperature to a specified value.
  ///
  int Set_Temp(double *,     ///< primitive vector.
	       const double, ///< temperature
	       const double  ///< eos gamma.
	       );

  //
  // Temporary function
  //
  //void Interface_with_JMs_code(
  //            std::vector<double>, ///< Yinput
  //            std::vector<double>, ///< Youtput
  //            double, ///< density (g/cm3)
  //            double, ///< input internal energy (erg/cm3)
   //           double &, ///< output internal energy (erg/cm3)
   //           double, ///< Column density (g/cm2)
   //           double  ///< timestep to integrate
   //           ) {return;}

 protected:
  ///
  /// convert state vector from grid cell into local microphysics vector.
  /// You can do this however you like.  Primitive vector contains mass
  /// density, pressure, and tracers for various species which must be
  /// dimensionless fractional abundances.  They don't have to sum up to
  /// unity or anything, just they must be dimensionsless numbers.
  ///
  int convert_prim2local(const double *, ///< primitive vector from grid cell (length nv_prim)
			 const double    ///< eos gamma
			 );

  ///
  /// Convert local microphysics vector into state vector for grid cell.
  /// This is the inverse of convert_prim2local.
  ///
  int convert_local2prim(
            const double *, ///< input primitive vector from grid cell (length nv_prim)
            double *,       ///< updated primitive vector for grid cell (length nv_prim)
            const double    ///< eos gamma.
            );

  ///
  /// This function organises a 5-element array of data for Harpreet's module:
  ///  - the neutral H column density to point source,
  ///  - the total column density to point sourc, 
  ///  - and the H2 column density to point sourc,
  ///  - distance to point source
  ///  - minimum column density to diffuse radiation.
  ///
  void get_column_densities(
                    const int,      ///< Number of UV heating sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of UV-heating column densities and source properties.
                    const int,      ///< number of ionising radiation sources.
                    const std::vector<struct rt_source_data> &,
                    ///< list of ionising src column densities and source properties.
                    std::vector<double> & ///< array for column densities
                    );

 private:
  //
  // Any data which can be calculated/set at the start of the simulation
  // can be defined here.
  //
  const double kB;  ///<  Boltzmanns constant.
  const double m_p; ///< Mass of proton.
  double gamma;     ///< EOS gamma for ideal gas.
  struct which_physics ep; ///< struct with flags for which extra physics we are or aren't doing.

  const int nv_prim; ///< Number of variables in state vector.
  int  Nspecies;     ///< number of chemistry variables in Yi,Yf vectors.
  int Yvector_length; ///< Length of Yi,Yf vectors, equals Nspecies+1 (for internal energy).
  std::vector<double> Yi, Yf; ///< Vector for the species, input and output.
  double
    density, ///< mass density within current cell to update. (g/cm3)
    Eint,    ///< internal energy density (ergs/cm3)
    Column_density; ///< Minimum projected mass density in the 2N coordinate directions.

  bool have_pt_src; ///< false if no point sources.
  bool have_diff_r; ///< false if no diffuse radiation field.
  bool col_data_set; ///< false initially, set to true once data for srcs set.
  int 
    ion_src_index,  ///< index of ionising source in ionising_srcs[v]
    uvh_src_index,  ///< index of dissociating source in heating_srcs[v]
    dis_src_index;  ///< index of dissociating source in heating_srcs[v]
  std::vector<int> ds_index; ///< indices of diffuse srcs in heating_srcs[]
};



#endif // if not excluding Harpreet's module.
#endif // MICROPHYSICS_LOWZ
