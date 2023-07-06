///
/// \file sim_params.h
/// \author Jonathan Mackey
///
/// Structure with simulation parameters.
///
/// Modifications:
/// - 2015.01.09 JM: Wrote file, moved stuff from global.h
/// - 2015.07.03 JM: Started to change tracer setup in files.
///    Put units and JetParams classes into this file from global.h.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
/// - 2015.10.19 JM: Fixed wind-tracer to always use pion_flt.
/// - 2017.11.07-22 JM: updating boundary setup.

#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

#include <array>
#include <string>
#include <vector>

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "constants.h"
#include "dataIO/readparams.h"
#include "raytracing/rad_src_data.h"
#include "sim_constants.h"

#ifdef PARALLEL
#include "sub_domain/sub_domain.h"
#endif


// ##################################################################
// ##################################################################


///
/// struct with lots of flags for what extra physics we are using.
///
/// All flags set to 0 (not used) by default, and the initial
/// conditions file can set whichever ones it wants, (or command line
/// overrides).  Some of these must be set if others are set...
/// i.e. they are not all independent.
///
struct which_physics {
  /// if we are doing dynamics or not (1=yes,0=no).
  int dynamics;
  /// if we are doing ray-tracing or not (1=yes,0=no).
  int raytracing;
  /// radiative cooling (1-8=yes,0=no).
  int cooling;
  /// if we have multiple chemical species and solve rate equations
  /// (1=yes,0=no).
  int chemistry;
  /// collisional ionisation (1=yes,0=no).
  int coll_ionisation;
  /// radiative recombination (1=yes,0=no).
  int rad_recombination;
  /// if we have photo-ionisation (with or w/o raytracing) (1=yes,0=no).
  int phot_ionisation;

  /// if we are using thermal conduction (1=yes, 0=no)
  int sat_thermal_cond;
  /// the value of the TC scaling (1.0=0.33*Spitzer)
  double tc_strength;

  /// flag for whether to update the energy with the rates (1=yes,0=no).
  int update_erg;
  ///
  /// Set this if we want to limit timesteps by cooling/recomb times.
  /// - limit=0 : only dynamical timestep
  /// - limit=1 : limit by dynamical and cooling/heating times
  /// - limit=2 : limit by dynamical + cooling/heating + recombination times.
  /// - limit=3 : limit by dyn. +cool/heat +recomb +ionisation times.
  /// - limit=4 : limit by dyn. +recomb times (NOT cooling/heating).
  ///
  int MP_timestep_limit;
  /// Minimum temperature to allow in the simulation.
  double MinTemperature;
  /// Maximum temperature to allow in the simulation.
  double MaxTemperature;

  /// Mass fraction of H, X, used for calculating mean mass per
  /// particle by assuming the rest is He.
  double H_MassFrac;
  /// Mass fraction of He, Y, used for calculation electron/ion
  /// densities as a function of H number density, and for setting
  /// the mean mass per particle, mu.
  double Helium_MassFrac;
  /// microphysics (but doesn't contribute to mean mass per particle,
  /// mu).
  double Metal_MassFrac;

  /// If including Compton cooling
  int compton_cool;
  /// If including wind acceleration from stars
  int wind_acceleration;
};



// ##################################################################
// ##################################################################


/// Struct to hold parameters for the GLM-MHD divergence-cleaning
/// method
struct glm_mhd {
  /// value in [0,1] to scale hyperbolic wave speed by (default=1)
  double glm_chyp_multiplier;
  /// max value of exponent for parabolic dissipation of divB (default=0.3)
  double glm_par_limiter;
};


///
/// Star struct, for storing data from a stellar evolution file.
///
struct star {
  std::string file_name;
  ///< Name of stellar evolution file.
  size_t Nlines,  ///< Number of lines in file.
      last_line;
  ///< line in file corresponding to most recent time t<=t_sim
  int src_id;  ///< ID of source in SimPM.RS
  // fuv_src_id;      ///< ID of (possible) Far-UV heating source in SimPM.RS
  std::vector<double> time,  ///< Array of time, in seconds.
      Log_L,                 ///< log10 of Luminosity (Lsun)
      Log_T,                 ///< Log10 of Teff (K)
      Log_R,                 ///< Log10 of Radius (solar radii)
      Log_V;                 ///< Log10 of Wind (cm/s)
  double t_now,              ///< CurrentTime in sec.
      Lnow,                  ///< Current Luminosity in erg/s
      Tnow,                  ///< Current Temperature (K)
      Rnow,                  ///< Current Radius (Rsun)
      Vnow;                  ///< Current wind velocity (cm/s).
};



// ##################################################################
// ##################################################################



/// This struct contains global data for a stellar wind source (for
/// reading and writing to file).  The grid setup functions should
/// read this struct and set up the appropriate wind boundary.
///
struct stellarwind_params {
  /// time-dep position of source (input/stored cm).
  std::array<double, MAX_DIM> dpos;
  /// previous assigned position of source (input/stored cm).
  std::array<double, MAX_DIM> lastpos;
  /// current assigned position of source (input/stored cm).
  std::array<double, MAX_DIM> thispos;
  /// velocity vector of star (calculated)
  std::array<double, MAX_DIM> velocity;
  double Mdot;    ///< mass loss rate (input Msun/yr, stored g/s).
  double Mass;    ///< Mass of star (input Msun, stored g).
  double Vinf;    ///< wind terminal velocity (input km/s, stored cm/s)
  double Vrot;    ///< wind rotational velocity (input km/s, stored cm/s)
  double Vcrit;   ///< wind critical velocity (calculated cm/s)
  double Rstar;   ///< Radius of star (input/stored cm).
  double Bstar;   ///< Surface split-monopole magnetic field strength (G)
  double Tstar;   ///< stellar temperature (sets pressure at r=Rstar, K).
  double radius;  ///< Radius of boundary region (input cm, stored cm).
  double current_radius;  ///< current (modified) radius of boundary region
  double time_offset;     ///< time offset between wind-data-file and sim-time
                          ///< (input yr, stored s).
  double update_freq;     ///< how often to update wind-data from file info
                          ///< (input yr, stored s).
  /// Evolution of wind happens this factor faster than normal (for factor>1)
  /// Should probably only be used for Main Sequence evolution!
  double t_scalefactor;
  /// for rotating stars, this sets how strongly the mass-loss is
  /// focussed onto the equator.  It is a power law index on a term
  /// (1-Omega*sin(theta))^xi.  Default is xi = -0.43.
  /// See e.g. Langer, Garcia-Segura & Mac Low (1999,ApJ,520,L49).
  double xi;

#ifdef ANALYTIC_ORBITS
  double eccentricity;    ///< eccentricity of the orbit
  double OrbPeriod;       ///< orbital period (input in years, stored in s)
  double PeriastronX;     ///< (x_com - x_star) at t=0 at periastron (cm)
  double PeriastronY;     ///< (y_com - y_star) at t=0 at periastron (cm)
#endif                    // ANALYTIC_ORBITS
  pion_flt tr[MAX_NVAR];  ///< tracer values in wind at Rstar.
  int id;            ///< if we have multiple sources, this identifies them.
  int type;          ///< what type of stellar wind?  see stellar_wind.h
  int moving_star;   ///< is the star moving? 1=yes, 0=no
  int acc;           ///< 0 = launch at vinf; 1 = accelerate with beta=1
  int enhance_mdot;  ///< 0=no, 1=yes, for rapidly rotating stars.
  std::string
      evolving_wind_file;  ///< name of file containing evolving wind data.
};

struct stellarwind_list {
  std::vector<struct stellarwind_params *> params;  ///< list of params.
  int Nsources;                                     ///< number of sources.
};

extern struct stellarwind_list SWP;  ///< global scope list of wind sources


// ##################################################################
// ##################################################################



///
/// Grid data for each level of the NG grid.
///
class level {
public:
  class GridBaseClass *grid;    ///< grid at this level
  class GridBaseClass *parent;  ///< pointer to parent grid.
  class GridBaseClass *child;   ///< pointer to child grid.
  long int Ncell;  ///< Total number of 'real' grid zones within the domain
                   ///< (Total for level).
  long int step;   ///< number of steps taken on this level.
  std::array<double, MAX_DIM>
      Range;  ///< Size of domain in x,y,z-direction on this level
  std::array<double, MAX_DIM>
      Xmin;  ///< Min value of x,y,z in domain on this level
  std::array<double, MAX_DIM>
      Xmax;        ///< Max value of x,y,z in domain on this level
  double dx;       ///< Linear side length of grid cells on this level
  double simtime;  ///< simulation time at this level.
  double dt;       ///< current timestep at this level.
  double last_dt;  ///< previous timestep at this level.
  std::array<int, MAX_DIM> NG;  ///< Number of 'real' grid zones in each
                                ///< direction (Total for level).
  int multiplier;               ///< 2^l, l=0=coarsest.
#ifdef PARALLEL
  class Sub_domain sub_domain;  ///< domain decomposition on this level
#endif
};



// ##################################################################
// ##################################################################



///
/// \brief Simulation Parameters Class
///
/// This is a holder for various flags and settings for the simulation.
///
class SimParams {
public:
  SimParams();
  SimParams(const std::string);
  ~SimParams();
  int gridType;    ///< Uniform, Finite Volume, cubic cell grid: 1 = only
                   ///< option.
  int eqntype;     ///< Euler=1, Ideal-MHD=2, GLM-MHD=3, FCD-MHD=4, IsoHydro=5.
  int coord_sys;   ///< Cartesian=1, Cylindrical=2, Spherical=3
  int solverType;  ///< 0=Lax-Friedrichs,1=LinearRS,2=ExactRS,3=HybridRS.
  int eqnNDim;     ///< Dimensionality of equations, set to 3 for now.
  int ndim;        ///< Dimensionality of grid (can be one of [1,2,3]).
  int nvar;        ///< Length of State Vectors (number of variables).

  // Tracer variables:
  int ntracer;  ///< Number of tracer variables.
  int ftr;      ///< Position of first tracer variable in state vector.
  std::string chem_code;  ///< Code for what kind of chemistry we are running.
  std::string
      TRTYPE;  ///< LEGACY CODE for what kind of chemistry we are running.
  std::string *tracers;  ///< array of strings for the tracer type

  // Timing
  double simtime;     ///< current time in simulation.
  double starttime;   ///< initial time to start simulation at.
  double finishtime;  ///< Time at which to finish the simulation.
  bool maxtime;  ///< False if simulation is to continue, true if time to stop.
  int timestep;  ///< Integer count of the number of timesteps since the
                 ///< start.
  double dt;     ///< timestep size for simulation (all cells have same step).
  double last_dt;       ///< Remember the last timestep.
  double min_timestep;  ///< Minimum value timestep can be; if
                        ///< dt<min_timestep, but out.

  // Grid data
  std::array<int, MAX_DIM> NG;  ///< Number of 'real' grid zones in each
                                ///< direction (Total for sim) (top level).
  long int Ncell;  ///< Total number of 'real' grid zones (within the range)
                   ///< (Total for sim) (top level).
  std::array<double, MAX_DIM>
      Range;  ///< Size of domain in x,y,z-direction (top level).
  std::array<double, MAX_DIM>
      Xmin;  ///< Min value of x,y,z in domain (top level).
  std::array<double, MAX_DIM>
      Xmax;   ///< Max value of x,y,z in domain (top level).
  double dx;  ///< Linear side length of (uniform, cubic, cartesian) grid
              ///< cells (top level).

  // NG grid data
  int grid_nlevels;
  int grid_aspect_ratio[MAX_DIM];
  double NG_centre[MAX_DIM];
  int NG_refine[MAX_DIM];
  std::vector<class level> levels;

  // Boundary cell data.
  std::string BC_XN;      ///< Type of boundary condition.
  std::string BC_XP;      ///< Type of boundary condition.
  std::string BC_YN;      ///< Type of boundary condition.
  std::string BC_YP;      ///< Type of boundary condition.
  std::string BC_ZN;      ///< Type of boundary condition.
  std::string BC_ZP;      ///< Type of boundary condition.
  int BC_Nint;            ///< Number of internal boundary regions
  std::string *BC_INT;    ///< List of internal boundary regions.
  std::string BC_STRING;  ///< For reading pre-2018 data-files.

  int Nbc;     ///< Depth of boundary/ghost cells from edge of grid.
  int Nbc_DD;  ///< Depth of boundary cells for MPI boundaries.
  // Integration accuracy
  int spOOA;  ///< Spatial Order of Accuracy in the code.
  int tmOOA;  ///< Time Order of Accuracy in the code.

  // Physics
  double gamma;  ///< Ideal gas constant.
  double CFL;  ///< Courant factor, must be less than one, should be about 0.5.
  int artviscosity;  ///< Integer flag. 0=No Artificial Viscosity; 1=Falle's
                     ///< version; 2=Lapidus; 3=H-correction.
  double etav;  ///< Artificial viscosity coefficient, should be between 0.01
                ///< and 0.3
  struct which_physics
      EP;                 ///< flags for what extra physics we are going to use.
  struct rad_sources RS;  ///< list of radiation sources.
  std::vector<struct star> STAR;  ///< Data from stellar evolution file(s).
  struct glm_mhd glm_data;        ///< parameter for the Dedner mixed-GLM method

  // File I/O
  int typeofop;  ///< Output FileType: Integer flag. 1=ascii, 2=fits,
                 ///< 3=fitstable, 4=FITSandASCII, 5=Silo etc.
  int typeofip;  ///< Input FileType:  Integer flag. 1=ascii, 2=fits,
                 ///< 3=fitstable, 4=FITSandASCII, 5=Silo etc.

  std::string outFileBase;  ///< Filename, with path, to write data to.
  int opfreq;               ///< Output file every 'opfreq'th timestep.
  int op_criterion;         ///< =0 for per n-steps, =1 for per n-years.
  double next_optime;   ///< if op_criterion=1, then this is the next time an
                        ///< output will happen.
  double opfreq_time;   ///< if op_criterion=1, then this is the output
                        ///< frequency in seconds.
  int checkpoint_freq;  ///< how often to output a checkpoint file (#
                        ///< timesteps).

  // Initial Conditions.
  double addnoise;  ///< Whether to add noise or not, 0=no, 1-3 are for
                    ///< different types of noise.
  pion_flt RefVec[MAX_NVAR];  ///< Reference state vector for simulation.

  std::vector<int> get_pbc_bools()
      const;  ///< Return a boolean vector for periodic boundary conditions

private:
  int read_gridparams(const std::string  ///< parameter file.
  );

  /** \brief Reads in extra physics params from the text file, if they are
   * there. */
  int read_extra_physics();

  /** \brief if we are doing raytracing, read in the source list. */
  int read_radsources();

  /** \brief Reads in units params from the text file, if they are there. */
  int read_units();

  ///
  /// If we got one or more stellar wind sources, this function reads
  /// their properties
  ///
  int read_wind_sources();

  ///
  /// Read in parameters for a stellar jet, if requested by N_jet !=0
  ///
  int read_jet_params(
      class JetParams &jpar  ///< pointer to jet parameters class.
  );

  class ReadParams rp;
};

// *******************************************************************
// *******************************************************************
// extern class SimParams SimPM;
// *******************************************************************
// *******************************************************************

// *******************************************************************
///
/// Conversion factors between code units and physical units.  Currently
/// not used, as the code can handle cgs units without much trouble.
///
class Units {
public:
  std::string unitsys;  ///< Name of System of Units: MKS, CGS
  std::string density,  ///< Reference Density Units.
      length,           ///< Reference Length Units.
      velocity,         ///< Reference Velocity Units.
      bfield;           ///< Reference B-Field Units.
  double rhoVal,        ///< One code unit is this number of reference units.
      lenVal,           ///< One code unit is this number of reference units.
      velVal,           ///< One code unit is this number of reference units.
      magVal;           ///< One code unit is this number of reference units.
};

// *******************************************************************
// *******************************************************************
extern class Units
    uc;  ///< Unit conversion quantities from code to physical units.
// *******************************************************************
// *******************************************************************

// *******************************************************************
///
/// If we are running a jet simulations, the physical properties of the
/// jet are stored in this function.
///
class JetParams {
public:
  JetParams();
  ~JetParams();
  int jetic;      ///< 0=not a jet sim, 1=is a jet sim
  int jetradius;  ///< Radius of jet in units of cellsize. Jet always centred
                  ///< at origin of XN boundary.
  pion_flt *jetstate;  ///< State vector for jet inflow.
};
// *******************************************************************
// *******************************************************************
extern class JetParams JP;
// *******************************************************************
// *******************************************************************

#endif  // SIM_PARAMS_H
