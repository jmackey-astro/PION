///
/// \file sim_params.h
/// \author Jonathan Mackey
///
/// Structure with simulation parameters.
///
/// Modifications:
/// - 2015.01.09 JM: Wrote file, moved stuff from global.h

#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

#include <string>
#include <vector>
#include "sim_constants.h"

// *******************************************************************
///
/// struct with lots of flags for what extra physics we are using.
/// 
/// All flags set to 0 (not used) by default, and the initial
/// conditions file can set whichever ones it wants, (or command line
/// overrides).  Some of these must be set if others are set...
/// i.e. they are not all independent.
///
struct which_physics {
  int dynamics;
    ///< if we are doing dynamics or not (default is yes).
  int raytracing;
    ///< if we are doing ray-tracing or not.
  int cooling;
    ///< radiative cooling.
  int chemistry;
    ///< if we have multiple chemical species and solve rate equations.
  int coll_ionisation;
    ///< collisional ionisation.
  int rad_recombination;
    ///< radiative recombination.
  int phot_ionisation;
    ///< if we have photo-ionisation (with or w/o raytracing).
  int update_erg;
    ///< flag for whether to update the energy with the rates.
  ///
  /// Set this if we want to limit timesteps by cooling/recomb times.
  /// - limit=0 : only dynamical timestep
  /// - limit=1 : limit by dynamical and cooling/heating times
  /// - limit=2 : limit by dynamical + cooling/heating + recombination times.
  /// - limit=3 : limit by dyn. +cool/heat +recomb +ionisation times.
  /// - limit=4 : limit by dyn. +recomb times (NOT cooling/heating).
  ///
  int MP_timestep_limit;
  double MinTemperature; ///< Minimum temperature to allow in the simulation.
  double MaxTemperature; ///< Maximum temperature to allow in the simulation.
#ifdef THERMAL_CONDUCTION
  int thermal_conduction; ///< 0 if no conductivity, 1 if using it.
#endif // THERMAL CONDUCTION

  ///
  /// Mass fraction of H, X, used for calculating mean mass per
  /// particle by assuming the rest is He.
  ///
  double H_MassFrac;
  ///
  /// Mass fraction of He, Y, used for calculation electron/ion
  /// densities as a function of H number density, and for setting
  /// the mean mass per particle, mu.
  ///
  double Helium_MassFrac;
  ///
  /// Mass fraction of metals, Z, used for heating/cooling in
  /// microphysics (but doesn't contribute to mean mass per particle,
  /// mu).
  ///
  double Metal_MassFrac;
};
// *******************************************************************


// *******************************************************************
///
/// Star struct, for storing data from a stellar evolution file.
/// I need to get rid of Rstar, because it is not used, and is in the
/// wrong units.
///
struct star {
  std::string file_name;
    ///< Name of stellar evolution file.
  size_t
    Nlines,    ///< Number of lines in file.
    last_line;
      ///< line in file corresponding to most recent time t<=t_sim
  int
    src_id; ///< ID of source in SimPM.RS
    //fuv_src_id;      ///< ID of (possible) Far-UV heating source in SimPM.RS
  std::vector<double>
    time,  ///< Array of time, in seconds.
    Log_L, ///< log10 of Luminosity (Lsun)
    Log_T, ///< Log10 of Teff (K)
    Log_R, ///< Log10 of Radius (solar radii)
    Log_V; ///< Log10 of Wind (cm/s)
  double
    t_now, ///< CurrentTime in sec.
    Lnow,  ///< Current Luminosity in erg/s
    Tnow,  ///< Current Temperature (K)
    Rnow,  ///< Current Radius (Rsun)
    Vnow;  ///< Current wind velocity (cm/s).
};


///
/// Radiation source struct.
///
struct rad_src_info {
  double pos[MAX_DIM]; ///< src position (physical units).
  double strength; ///< src strength (photons/sec, or ergs/sec for multifreq.)
  double Rstar; ///< stellar radius in solar radii (for multifreq. photoionisation).
  double Tstar; ///< stellar effective temperature (for multifreq. photoionisation).
  int id;   ///< src identifier
  int type; ///< src type: either RT_SRC_DIFFUSE or RT_SRC_SINGLE.
  int update; ///< how the source is updated: RT_UPDATE_IMPLICIT=1, RT_UPDATE_EXPLICIT=2
  int at_infinity; ///< set to true if source is at infinity.
  ///
  /// "effect" is what the source does, and this defines many of its
  /// properties implicitly.  Options are:
  /// - RT_EFFECT_UV_HEATING,
  /// - RT_EFFECT_PION_MONO,
  /// - RT_EFFECT_PION_MULTI,
  /// - RT_EFFECT_PHOTODISS, (UNUSED)
  /// - RT_EFFECT_PION_EQM, (UNUSED--for photoionisation equilibrium)
  /// - RT_EFFECT_HHE_MFQ
  ///
  int effect;
  ///
  /// "NTau" sets the number of quantities traced from the source.
  ///
  int NTau;

  int opacity_src; ///< What provides the opacity: RT_OPACITY_TOTAL, RT_OPACITY_MINUS, RT_OPACITY_TRACER.
  int opacity_var; ///< optional tracer variable index in state vector, for opacity calculation.
  std::string EvoFile;  ///< Optional text file with output from stellar evolution code for time-varying source.
};

///
/// Basic list of all radiation sources.  Used in SimParams class.
///
struct rad_sources {
  int Nsources;
  std::vector<struct rad_src_info> sources;
};
// *******************************************************************



// *******************************************************************
///
/// This struct contains global data for a stellar wind source (for
/// reading and writing to file).  The grid setup functions should
/// read this struct and set up the appropriate wind boundary.
///
struct stellarwind_params {
  int id;    ///< if we have multiple sources, this identifies them.
  int type;  ///< what type of stellar wind?  see stellar_wind.h
  double dpos[MAX_DIM]; ///< position of source (physical coords).
  double Mdot;  ///< mass loss rate in solar masses per year.
  double Vinf;  ///< wind terminal velocity in km/s.
  double Rstar; ///< radius at which to set pressure based on Tstar.
  double Tstar; ///< stellar temperature (sets pressure at r=Rstar).
  double tr[MAX_NVAR]; ///< tracer values in wind at Rstar.
  double radius;
  ///< Radius out to which boundary condition is imposed (physical).
  std::string evolving_wind_file; ///< name of file containing evolving wind data.
  double time_offset;   ///< time offset between wind-data-file and sim-time (YEARS!).
  double update_freq;   ///< how often to update wind-data from file info (YEARS!).
  /// Evolution of wind happens this factor faster than normal (for factor>1)
  /// Should probably only be used for Main Sequence evolution!
  double t_scalefactor;
};

struct stellarwind_list {
  std::vector<struct stellarwind_params *> params; ///< list of params.
  int Nsources; ///< number of sources.
};

extern struct stellarwind_list SWP;
// *******************************************************************





// *******************************************************************
///
/// \brief Simulation Parameters Class
/// 
/// This is a holder for various flags and settings for the simulation.
///
class SimParams {
  public:
   SimParams();
   ~SimParams();
   int gridType;   ///< Uniform, Finite Volume, cubic cell grid: 1 = only option.
   int eqntype;    ///< Euler=1, Ideal-MHD=2, GLM-MHD=3, FCD-MHD=4, IsoHydro=5.
   int coord_sys;  ///< Cartesian=1, Cylindrical=2, Spherical=3
   int solverType; ///< 0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS.
   int eqnNDim;    ///< Dimensionality of equations, set to 3 for now.
   int ndim;       ///< Dimensionality of grid (can be one of [1,2,3]).
   int nvar;       ///< Length of State Vectors (number of variables).
   int ntracer;    ///< Number of tracer variables.
   int ftr;        ///< Position of first tracer variable in state vector.
   std::string trtype;  ///< String saying what type of tracer we are using.
   // Timing
   double simtime;    ///< current time in simulation. 
   double starttime;  ///< initial time to start simulation at. 
   double finishtime; ///< Time at which to finish the simulation.
   bool maxtime;      ///< False if simulation is to continue, true if time to stop.
   int timestep;      ///< Integer count of the number of timesteps since the start. 
   double dt;         ///< timestep size for simulation (all cells have same step).
   double last_dt;    ///< Remember the last timestep, and don't increase the next step by more than 30%.
   double min_timestep; ///< Minimum value timestep can be.  If it falls below this then bug out!

   // Grid Point data
   int NG[MAX_DIM];       ///< Number of 'real' grid-points in each direction (Total for sim)
   long int Ncell;            ///< Total number of 'real' grid points (within the range) (Total for sim).
   double Range[MAX_DIM]; ///< Size of domain in x,y,z-direction.
   double Xmin[MAX_DIM];  ///< Min value of x,y,z in domain.
   double Xmax[MAX_DIM];  ///< Max value of x,y,z in domain.
   double dx;            ///< Linear side length of (uniform, cubic, cartesian) grid cells.
   // Boundary cell data.
   std::string typeofbc; ///< Type of boundary condition(s).
   int Nbc;         ///< Depth of boundary/ghost cells from edge of grid.
   // Integration accuracy
   int spOOA;  ///< Spatial Order of Accuracy in the code.
   int tmOOA;  ///< Time Order of Accuracy in the code.

   // Physics
   double gamma;            ///< Ideal gas constant.
   double CFL;              ///< Courant factor, must be less than one, should be about 0.5.
   int artviscosity;        ///< Integer flag. 0=No Artificial Viscosity; 1=Falle's version;
                            ///< 2=Lapidus; 3=H-correction.
   double etav;             ///< Artificial viscosity coefficient, should be between 0.01 and 0.3
   struct which_physics EP; ///< flags for what extra physics we are going to use.
   struct rad_sources   RS; ///< list of radiation sources.
   std::vector<struct star>  STAR; ///< Data from stellar evolution file(s).

   // File I/O
   int typeofop;       ///< Output FileType: Integer flag. 1=ascii, 2=fits, 3=fitstable, 4=FITSandASCII, 5=Silo etc.
   int typeofip;       ///< Input FileType:  Integer flag. 1=ascii, 2=fits, 3=fitstable, 4=FITSandASCII, 5=Silo etc.

   std::string outFileBase; ///< Filename, with path, to write data to.
   int opfreq;         ///< Output file every 'opfreq'th timestep.
   int op_criterion;   ///< =0 for per n-steps, =1 for per n-years.
   double next_optime; ///< if op_criterion=1, then this is the next time an output will happen.
   double opfreq_time; ///< if op_criterion=1, then this is the output frequency in seconds.
   int checkpoint_freq; ///< how often to output a checkpoint file (# timesteps).

   // Initial Conditions.
   double addnoise; ///< Whether to add noise or not, 0=no, 1-3 are for different types of noise.
   double RefVec[MAX_NVAR];  ///< Reference state vector for simulation.
};

extern class SimParams SimPM;
// *******************************************************************




#endif // SIM_PARAMS_H

