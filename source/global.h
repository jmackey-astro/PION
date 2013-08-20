/// \file global.h
/// \author Jonathan Mackey
///
/// Contains some global variables and classes.
///
/// Enums used to access elements of the state vectors are defined here.
///
/// Modifications:
///
/// - 2007-06-27 Worked on Equations Class.  Changed method of
///   handling directional operations.
/// - 2007-07-05 Equations class, added in GLMEqn class derived from
///   IdealMHDEqn.
/// - 2007-07-10 Moved 'rsvars' and 'waves' enums to riemann.h; bit of
///   documentation of GLM Equations.
/// - 2007-07-13 got new structure working.
/// - 2007-07-13 new structure was too complex, so i simplified it,
///    and got it working properly.
/// - 2007-10-11 Added interactive command line interface class.
/// - 2007-10-16 Moved equations classes into equations.h
///  - 2010-01-05 JM: Added SimPM.EP.MP_timestep_limit flag: true=use
///  cooling/recomb times to limit the timestep.  Also added a public
///  microphysics function timescales() which returns the minimum of a
///  selection of timescales.
/// - 2010-01-19 JM: Added new #define for VERY_TINY_VALUE=1.0e-200
/// - 2010-01-22 JM: Moved #defs for (NON_)CELL_CENTRED_SRC from raytracer_SC.h
///     to here.
/// - 2010-01-26 JM: Added get_iXmin/iXmax/iRange() functions to gridBaseClass
///     to get the integer positions.
/// - 2010-02-05 JM: Added offsets to parallelparams.
/// - 2010-02-05 JM: Added function returning all abutting domains.
/// - 2010-04-11 JM: Tidied up comments to have shorter lines.
/// - 2010-04-25 JM: Added min_timestep parameter to SimPM to bug out
///   if dt gets too small.
/// - 2010-07-20 JM: Made OA1,OA2 into #defs, set spOOA,tmOOA to be integers.
/// - 2010-07-21 JM: rep.error() needed a mod for serial code if the
///    error was in the output routine, since this led to infinite loop!
/// - 2010-07-24 JM: Added stellar wind class.
/// - 2010.10.01 JM: Added spherical coordinate system.
/// - 2010.10.04 JM: Added flags for extra screen output on Blastwave and
///    FieldLoop test problems (moved flags from gridMethods.cc).
/// - 2010.10.05 JM: Moved stellar winds to their own file.
/// - 2010.10.13 JM: Added a function to display command-line options.
/// - 2010.11.12-15 JM: Moved cell interface to its own files; minor
///   changes to defines and documentation.
/// - 2010.12.04 JM: Added distance functions to base Grid class which
///  change with the geometry of the grid.
///  This is in an ifdef for now (GEOMETRIC_GRID)
/// - 2010.12.27 JM: added call to defines/functionality_flags.h.
///   This file defines which parts of the code get compiled and
///   which parts are left out.
/// - 2010.12.28 JM: Added EQEUL_EINT=9 hash-def.
///
/// - 2011.01.06 JM: Moved all flags which affect what gets compiled
///   and what doesn't to external files, so that not every piece of
///   code needs to call global.h
/// - 2011.01.06 JM: New stellar wind interface.
/// - 2011.01.17 JM: Added checkpoint frequency to SimPM.
/// - 2011.01.31 JM: Added min/max temperature values to SimPM.EP
/// - 2011.02.15 JM: Added extra parameters for time-varying stellar winds.
/// - 2011.02.17 JM: Added hash-defines for types of RT.
/// - 2011.02.24 JM: Added Struct SimPM.RS to handle radiation sources more simply.
///     Also added a new function to RSP so that I can query it more easily.
/// - 2011.02.25 JM: Updated interface to ray-tracing (need to move it out of this file).
/// - 2011.02.28 JM: Got rid of ray-tracing parameters class (RSP pointer).
/// - 2011.03.02 JM: Got rid of MAX_NTR. Increased MAX_NVAR to 70.
/// - 2011.03.21 JM: Updated RT interface; added extra params to SimPM.RS struct, and
///    extra #defines for RT source properties. (and 22.03)
/// - 2011.04.06 JM: Added idifference_vertex2cell() for conductivity
///    calculation.
/// - 2011.04.15 JM: small changes to RT flags (and 04.17).
/// - 2011.06.02 JM: RefVec is now a stack array rather than dynamically allocated.
/// - 2011.10.14 JM: Commented out RT_DIFF class, added new RT interface functions.
/// - 2011.12.01 JM: Added GeneralStuff::root_find_linear()
///
/// - 2012.01.14 JM: Added 'star' struct, for data from stellar evolution file.
/// - 2012.01.16 JM: Added extra constants to GeneralStuff
/// - 2012.05.15 JM: Added global ixmin/ixmax/irange functions to grid base class
/// - 2013.02.14 JM: Added He/Metal mass fractions as EP parameters,
///    to make metallicity and mu into parameterfile settings.
/// - 2013.04.18 JM: Removed NEW_METALLICITY flag.
/// - 2013.08.12 JM: Added RT_EFFECT_PION_EQM flag for radiation
///    sources where I assume photoionisation equilibrium.
/// - 2013.08.19/20 JM: Changed GeneralStuff constants (will remove
///    them eventually).  Added new RT_EFFECT_HHE_MFQ flag for H-He
///    chem.  Added NTau to rad_src_info struct.

#ifndef GLOBAL_H
#define GLOBAL_H

//
// These tell code what to compile and what to leave out.
//
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


///
/// TO DO:\n
/// - Move the grid base class into its own header file and include it from here.
///

//#define LOW_IONISATION_COOLING_HACK ///< Sets exponential cooling c12-c15 only in cell with x<0.01
 

//
// If running a jet simulation, this slows the jet linearly for the
// outer 25 per cent of the radius.
//
//#define SOFTJET

//
// If this is set, we correct the total energy for the new B-field
// update.  Can also set a flag to do temporal averaging for
// v-cross-B.
//
// NOTE: THE FIELD-CD METHOD DIDN'T TURN OUT TOO WELL AND MAY NOT EVEN
// WORK ANYMORE. (IT DOESN'T -- NOT IMPLEMENTED IN UPDATED SOLVERS).
//
//#define DO_ENERGY_FIX_FIELDCD
//#define FIELDCD_TEMPORAL_AVERAGING




#define MAX_DIM 3 ///< Max. number of spatial dimensions for a simulation.

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

//
// It's a good idea to include cmath here, since nearly every file needs it!
//
#ifndef INTEL
#include <cmath>     // Header file from gcc
#else
#include <mathimf.h> // Header file from Intel Compiler
#endif


#ifdef TESTING
#include <readline/readline.h>
#include <readline/history.h>
#endif

using namespace std;
//-------------------------------------------------------
// END OF (MOST OF) HEADER INCLUDES
//--------------------------------------------------------

///
/// Function to output commandline options for code.  It is here to
/// avoid code duplication in main.cc and mainMPI.cc.
///
void print_command_line_options(int, char **);


//-------------------------------------------------------
// Defines for integer constants used to describe things.
//--------------------------------------------------------

#define EQEUL 1  ///< Equation type: 1=Euler equations.
#define EQMHD 2  ///< Equation type: 2=Ideal MHD equations (Falle, Komissarov, Joarder, 1998). (ONLY GOOD FOR 1D!)
#define EQGLM 3  ///< Equation type: 3=MHD equations with GLM method for divergence cleaning (Dedner etal, 2002).
#define EQFCD 4  ///< Equation type: 4=MHD equations with Field-CD method for ensuring divB=0 (Toth,2000). (DON'T USE!!!)
#define EQEUL_ISO 5 ///< Equation type: 5=Isothermal Hydrodynamics
#define EQMHD_ISO 6 ///< Equation type: 6=Isothermal MHD
#define EQGLM_ISO 7 ///< Equation type: 7=Isothermal MHD with Dedner et al. GLM method  for Divergence Cleaning
#define EQFCD_ISO 8 ///< Equation type: 8=Isothermal MHD with Field-CD method for ensuring divB=0 (Toth,2000).
/// Equation type: 9=Adiabatic Hydro, integrating INTERNAL and not
/// total energy.
#define EQEUL_EINT 9

#define COORD_CRT 1 ///< Coordinates: Cartesian coordinate system.
#define COORD_CYL 2 ///< Coordinates: Cylindrical coordinate system.
#define COORD_SPH 3 ///< Coordinates: Spherical coordinate system.

#define MAX_NVAR 70 ///< Maximum length of state vectors.
//#define MAX_NTR   5 ///< Maximum number of tracer variables (feel free to expand later!).

//
// Defines for stuff related to machine accuracy and truncation error.
//
#define MACHINEACCURACY 5.e-16
///< This should be found out before running the code!

///
/// The pressure floor for the riemann solver and cell updates. 
/// (P_G(floor)=refvec[PG]*BASEPG)
/// This may need to be changed -- it should be larger than the largest
/// pressure ratio in the simulation you want to run.  Maybe it should 
/// be a runtime parameter...
///
#define BASEPG 1.e-5

/// For negative density, set rho=BASE_RHO*refvec[RO]
#define BASE_RHO 1.0e-5

#define SMALLVALUE 1.0e-12
///< A small value, for comparing with something we expect to be of
///< order unity.

#define TINYVALUE 1.0e-100
///< If any conserved or primitive variables are smaller than this,
///< change it!!!

#define VERY_TINY_VALUE 1.0e-200
///< If any conserved or primitive variables are smaller than this,
///< change it!!!

#define VERY_LARGE_VALUE 1.0e100
///< not quite HUGEVALUE, but still bigger than a quantity should be
///< in most situations.

#define HUGEVALUE        1.0e+200
///< A huge number that none of my physical variables should ever
///< approach.

#define ONE_PLUS_EPS (1.0+SMALLVALUE)
///< A number just greater than unity, for rounding doubles down to
///< ints.

#define ONE_MINUS_EPS (1.0-SMALLVALUE)
///< A number just less than unity, for rounding doubles up to ints.

/// Error code for a function which should not be called!
#define DONT_CALL_ME 99

///
/// Definitions for cylindrical axes and directions w.r.t. cartesian
/// ones.
///
#define Zcyl  XX
#define Rcyl  YY
#define Tcyl  ZZ
#define ZNcyl  XN
#define ZPcyl  XP
#define RNcyl  YN
#define RPcyl  YP
#define TNcyl  ZN
#define TPcyl  ZP

///
/// Definitions for Spherical coordinate axes and directions w.r.t. cartesian
/// ones. [R,T,P] stands for [r,theta,phi]
///
#define Rsph  XX
#define Tsph  YY
#define Psph  ZZ
#define RNsph  XN
#define RPsph  XP
#define TNsph  YN
#define TPsph  YP
#define PNsph  ZN
#define PPsph  ZP

///
/// Type of Solver:
///
#define FLUX_LF       0
#define FLUX_RSlinear 1
#define FLUX_RSexact  2
#define FLUX_RShybrid 3
#define FLUX_RSroe    4
#define FLUX_RSroe_pv 5
#define FLUX_FVS      6

/** \brief Primitive Variables Enum 
 * The code requires that the velocities follow each other in order,
 * and the same for the Magnetic Field.  I'm not sure if it will work 
 * without them being in that order, but I'm writing the code assuming 
 * they are in order.
 */
enum primitive {
    RO=0, ///< Density
    PG=1, ///< Gas pressure
    VX=2, ///< x-velocity
    VY=3, ///< y-velocity
    VZ=4, ///< z-velocity
    BX=5, ///< x-magnetic field
    BY=6, ///< y-magnetic field
    BZ=7, ///< z-magnetic field
    SI=8  ///< The GLM scalar field, psi.
};
/** \brief Conserved Variables Enum */
enum conserved {
    RHO=0,  ///< Density
    ERG=1,  ///< Total Energy
    MMX=2,  ///< x-momentum
    MMY=3,  ///< y-momentum
    MMZ=4,  ///< z-momentum
    BBX=5,  ///< x-magnetic field
    BBY=6,  ///< y-magnetic field
    BBZ=7,  ///< z-magnetic field
    PSI=8   ///< The GLM scalar field, psi.
};

/** \brief Enum for each Spatial Dimension. */
enum axes {
    XX=0, ///< X-Direction.
    YY=1, ///< Y-Direction.
    ZZ=2  ///< Z-Direction.
};

/** \brief Direction Vectors enum. 
 * This gets used for pointing to elements of arrays,
 * so the six real directions must comprise the numbers
 * zero to five.  They must also be arranged so that e.g.
 * XP = 2*axes.XX +1, XN = 2*axes.XX, and similarly for 
 * Y and Z directions.
 */
enum direction{
    NO=-1,   ///< No Direction.
    XN=0,  ///< Negative x-direction x-
    XP=1,  ///< Positive x-direction x+
    YN=2,  ///< Negative y-direction y-
    YP=3,  ///< Positive y-direction y+
    ZN=4,  ///< Negative z-direction z-
    ZP=5   ///< Positive z-direction z+
}; // {x-,x+,y-,y+,z-,z+}

//
// Space and Time Integration Order of Accuracy
//
#define OA1 1 ///< First Order Space.
#define OA2 2 ///< Second Order Space.

//
// Viscosity identifier flags
//
#define AV_NONE        0 ///< No artificial viscosity.
#define AV_FKJ98_1D    1 ///< 1D linear viscosity of Falle, Komissarov, Joarder (1998).
#define AV_LAPIDUS     2 ///< Lapidus-type viscosity (BROKEN).
#define AV_HCORRECTION 3 ///< H-correction of Sanders et al. (1998).
#define AV_HCORR_FKJ98 4 ///< combination of the H-corr + FKJ98 (not needed really).
#define AV_VonNeuRicht 5 ///< von Neumann & Richtmeyer (1950) multi-D viscosity.

/** \brief struct with lots of flags for what extra physics we are using.
 * 
 * All flags set to 0 (not used) by default, and the initial conditions 
 * file can set whichever ones it wants, (or command line overrides).
 * Some of these must be set if others are set... i.e. they are not all
 * independent.
 * */
struct which_physics {
  int dynamics;   ///< if we are doing dynamics or not (default is yes).
  int raytracing; ///< if we are doing ray-tracing or not.
  int cooling;    ///< radiative cooling.
  int chemistry;  ///< if we have multiple chemical species and solve rate equations.
  int coll_ionisation;   ///< collisional ionisation.
  int rad_recombination; ///< radiative recombination.
  int phot_ionisation;   ///< if we have photo-ionisation (with or w/o raytracing).
  int update_erg;    ///< flag for whether to update the energy with the rates.
  ///
  /// Set this if we want to limit timesteps by cooling/recomb times.
  /// - limit=0 : only dynamical timestep
  /// - limit=1 : limit by dynamical and cooling/heating times
  /// - limit=2 : limit by dynamical + cooling/heating + recombination times.
  /// - limit=3 : limit by dyn. +cool/heat +recomb +ionisation times.
  /// - limit=4 : limit by dyn. +recomb times (NOT cooling/heating).
  ///
  int MP_timestep_limit;
//#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  double MinTemperature; ///< Minimum temperature to allow in the simulation.
  double MaxTemperature; ///< Maximum temperature to allow in the simulation.
//#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
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


//
// Method of calculating opacity and updating microphysics:
//
#define RT_UPDATE_IMPLICIT   1 ///< C2-ray update, where time-integration happens during ray-trace.
#define RT_UPDATE_EXPLICIT  2 ///< Only instantaneous column densities calculated during ray-trace.

//
// Types of source -- either a single source or one of 2*Ndim diffuse sources.
//
#define RT_SRC_DIFFUSE 2
#define RT_SRC_SINGLE  1

//
// Effect of source:
//
#define RT_EFFECT_UV_HEATING 1   ///< UV heating source.
#define RT_EFFECT_PION_MONO  2   ///< monochromatic photoion source.
#define RT_EFFECT_PION_MULTI 3   ///< multifrequency photoion source.
#define RT_EFFECT_PHOTODISS  4   ///< photodissociation of molecules.
#define RT_EFFECT_PION_EQM   5   ///< Assume photoion. equilibrium.
#define RT_EFFECT_HHE_MFQ    6   ///< Frank&Mellema H+He p-ion scheme

//
// Source of opacity
//
#define RT_OPACITY_TOTAL  1 ///< opacity is proportional to total mass density, so just integrate rho*dr
#define RT_OPACITY_MINUS  2 ///< opacity provided by tracer i: integrate rho*(1-y_i)*dr
#define RT_OPACITY_TRACER 3 ///< opacity provided by tracer i: integrate rho*y_i*dr
#define RT_OPACITY_VSHELL 4 ///< should never be used.  RT uses it internally to set Vshell in each cell.
#define RT_OPACITY_HALPHA 5 ///< Used only for analysis to calculate projected H-alpha emission.
#define RT_OPACITY_NII_FL 6 ///< Used only for analysis to calculate projected [NII] forbidden line emission.
#define RT_OPACITY_RR     7 ///< Recombination rate (for ph-ion.eqm.)
#define RT_OPACITY_HHE    9 ///< H0,He0,He+,Dust.

///
/// Star struct, for storing data from a stellar evolution file.
/// I need to get rid of Rstar, because it is not used, and is in the wrong
/// units.
///
struct star {
  string file_name; ///< Name of stellar evolution file.
  size_t
    Nlines,    ///< Number of lines in file.
    last_line; ///< line in file corresponding to most recent time t<=t_sim
  int
    src_id; ///< ID of source in SimPM.RS
    //fuv_src_id;      ///< ID of (possible) Far-UV heating source in SimPM.RS
  vector<double>
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
  string EvoFile;  ///< Optional text file with output from stellar evolution code for time-varying source.
};

///
/// Basic list of all radiation sources.  Used in SimParams class.
///
struct rad_sources {
  int Nsources;
  std::vector<struct rad_src_info> sources;
};



/** \brief Simulation Parameters Class
 * 
 * This is a holder for various flags and settings for the simulation.
 * */
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
   string trtype;  ///< String saying what type of tracer we are using.
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
   string typeofbc; ///< Type of boundary condition(s).
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
   vector<struct star>  STAR; ///< Data from stellar evolution file(s).

   // File I/O
   int typeofop;       ///< Output FileType: Integer flag. 1=ascii, 2=fits, 3=fitstable, 4=FITSandASCII, 5=Silo etc.
   int typeofip;       ///< Input FileType:  Integer flag. 1=ascii, 2=fits, 3=fitstable, 4=FITSandASCII, 5=Silo etc.

   string outFileBase; ///< Filename, with path, to write data to.
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



/** \brief Conversion factors between code units and physical units.  Currently
 * not used, as the code can handle physical units without much trouble.
 * */
class Units {
  public:
   string unitsys; ///< Name of System of Units: code,mks,cgs,esu,etc.
   string density, ///< Reference Density Units.
     length,       ///< Reference Length Units.
     velocity,     ///< Reference Velocity Units.
     bfield;       ///< Reference B-Field Units.
   double rhoVal,  ///< One code unit is this number of reference units.
     lenVal,       ///< One code unit is this number of reference units.
     velVal,       ///< One code unit is this number of reference units.
     magVal;       ///< One code unit is this number of reference units.
};

extern class Units uc; ///< Unit conversion quantities from code to physical units.



/** \brief If we are running a jet simulations, the physical properties of the
 * jet are stored in this function.
 * */
class JetParams {
  public:
   JetParams();
   ~JetParams();
   int jetic;        ///< 0=not a jet sim, 1=is a jet sim
   int jetradius;    ///< Radius of jet in units of cellsize. Jet always centred at origin of XN boundary.
   double *jetstate; ///< State vector for jet inflow.
};
extern class JetParams JP;




// *******************************************************************
// *** Stellar winds have a rather complicated boundary condition ****
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
  string evolving_wind_file; ///< name of file containing evolving wind data.
  double time_offset;   ///< time offset between wind-data-file and sim-time (YEARS!).
  double update_freq;   ///< how often to update wind-data from file info (YEARS!).
  /// Evolution of wind happens this factor faster than normal (for factor>1)
  /// Should probably only be used for Main Sequence evolution!
  double t_scalefactor;
};

struct stellarwind_list {
  vector<struct stellarwind_params *> params; ///< list of params.
  int Nsources; ///< number of sources.
};

extern struct stellarwind_list SWP;

//#include "stellar_wind_BC.h"
//extern class stellar_wind SW;

// *******************************************************************
// *******************************************************************


#ifdef PARALLEL

//
// integer flags for MPI communication labels.
//
#define BC_ANYtag 0 ///< works for either sort of communication.
#define BC_MPItag 1 ///< This is an integer tag on send/receive operations, to label that this communicates MPI boundary data.
#define BC_PERtag 2 ///< Integer tag to say it is for periodic BC.
#define BC_RTtag  3 ///< Integer tag to say we are transferring a radiative transfer column density tag.

/** \brief Class to hold parameters only relevant to parallel code.
 * 
 * For parallel processing, the domain is split between processors, so the 
 * local domain is smaller than the full domain.  That information is stored 
 * here.
 * */
class ParallelParams {
 public:
  ParallelParams();
  ~ParallelParams();
  int nproc;  ///< Number of processors.
  int myrank; ///< Which processor am I?
  int LocalNG[MAX_DIM];  ///< Number of 'real' grid-zones in each direction, on this processor's domain.
  int LocalNcell;  ///< Total number of 'real' grid zones in this processor's domain.
  int offsets[MAX_DIM]; ///< number of zones from this proc's 1st zone to the sim's 1st zone.
  double LocalXmin[MAX_DIM];  ///< Min value of x,y,z in Processor's domain.
  double LocalXmax[MAX_DIM];  ///< Max value of x,y,z in Processor's domain.
  double LocalRange[MAX_DIM]; ///< Size of Processor's domain in x,y,z-direction.
  int *ngbprocs;  ///< list with processor rank of neighbours in each direction.
  bool ReadSingleFile; ///< If the ICs are in a single file, set this to true.
  bool WriteSingleFile; ///< If you want all the processors to write to one file, set this (BUGGY!)
  bool WriteFullImage;  ///< If you want multiple fits files, but each one is the full domain size (bad!)
  /** \brief Decompose the domain into blocks for each processor, and set up
   * a structure which contains the domain of each processor.
   * */
  int decomposeDomain();
  ///
  /// Get a list of all abutting domains, including corner/edge intersections.
  ///
  void get_abutting_domains(std::vector<int> & ///< write list to this vector.
          );
  /** \brief Set the maximum runtime to a new value. Should be set in main(). */
  void set_max_walltime(double ///< New Max. runtime in seconds.
      );
  /** \brief Get the maximum runtime in seconds. */
  double get_max_walltime();
  //
  // Returns the ix array for any requested rank.
  //
  void get_domain_ix(const int, ///< rank ix requested for.
         int *      ///< array to put ix into.
         );
 protected:
  int ix[MAX_DIM];  ///< this proc's position in the block of domains (zero offset)
  int nx[MAX_DIM];  ///< size of block of domains in each direction. (one block = one unit).
  std::vector<int> full_ngb_list; ///< list of abutting domains.
  double max_walltime; ///< Max. walltime to run for, in seconds, after which we output data and finish.
  /** \brief Celled by decomposeDomain() to set neighbouring processor ids.
   * 
   * This works if processors are ranked from zero to n-1 in integer steps.
   * If a subdomain is on the full domain boundary, the processor neighbour
   * in that direction is set to -999.
   * 
   * The rank of processor i is defined by
   * \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
   * where the 'n's refer to number of processors that span the domain
   * in each direction, and 'i's refer to the location of the current
   * processor along that direction.
   * */
  int pointToNeighbours();
};
extern class ParallelParams mpiPM; ///< Holds data about current processor's subdomain.
#endif //PARALLEL




/** \brief General Purpose class to hold generally useful functions and data.
 * 
 * e.g. contains physical constants.  Also a function to determine if two
 * doubles are the same to within some tolerance.  Also a function to calculate
 * the distance between two vector positions.
 * */
class GeneralStuff {
  public:
   GeneralStuff();
   ~GeneralStuff();
   /** \brief tests if two doubles are equal to an accuracy of 1e-8. */
   bool equalD(const double, 
         const double
         );
   /** \brief Sets up Cubic Spline interpolation (from Martin White's Code, from NR) */
   void spline(const double *, ///< Array of x values.
         const double *, ///< Array of y values.
         const int ,     ///< Length of arrays.
         double ,  ///< First Derivative of interpolating function at x[1] (unit offset array) (>1.e30 for natural spline)
         double ,  ///< First Derivative of interpolating function at x[n] (unit offset array) (>1.e30 for natural spline)
         double *  ///< Empty array to store d2y/dx2 in.
         );
   /** \brief Performs cubic spline interpolation to get y(x) (from Martin White's Code, from NR) */
   void splint(const double *, ///< Array of x values.
         const double *, ///< Array of y values.
         const double *, ///< Array of d2y/dx2 values.
         const int ,     ///< nspl
         const double ,  ///< x we are searching for.
         double *  ///< pointer to result.
         );

  ///
  /// Given a vector of x-values, and corresponding y-values, and an input
  /// x value, find the corresponding y-value by bisection and then linear
  /// interopolation.
  ///
  void root_find_linear(
        const double *, ///< Array of x values.
        const double *, ///< Array of y values.
        const size_t,   ///< Array sizes
        const double ,  ///< x we are searching for.
        double *  ///< pointer to result.
        );

   /** \brief calculate distance between two points in n-dimensional space.
    * Assumes cartesian geometry.
    */
   double distance(const double *, ///< position 1
       const double *, ///< position 2
       const int       ///< number of spatial dimensions.
       );
   /** \brief calculate distance between two points in n-dimensional space.
    * Assumes cartesian geometry.
    */
   double idistance(const int *, ///< position 1
        const int *, ///< position 2
        const int    ///< number of spatial dimensions.
        );
   /** \brief start a timer, identified by a string. If timer exists already, this
    * function assumes it has been paused and sets it counting again. 
    * */
   void start_timer(string);
   /** \brief Pause a timer, identified by a string, and return time in seconds. 
    * Sets time value to number of seconds it has been running. */
   double pause_timer(string);
   /** \brief Return total time since timer started, and keep running. */
   double time_so_far(string);
   /** \brief stop a timer, identified by a string, and return time in seconds. 
    * Deletes the timer. */
   double stop_timer(string);

   inline double kB()  {return 1.381e-16;}
   inline double m_p() {return 1.673e-24;}
   inline double Ith_H() {return 13.59844;}
   inline double K_per_eV() {return 1.16045e4;}
   inline double s_per_yr() {return 3.15576e7;}
   inline double Msun() {return 1.989e33;}
   inline double parsec() {return 3.086e18;}
   inline double Lsun() {return 1.989e33;}
   inline double Rsun() {return 6.960e10;}
   inline double StefanBoltzmannConst() {return 5.670e-5;}
   inline double ln10() {return 2.302585093;}
   inline double c() {return 2.9979e10;}
    private:
   map<string,double> timers; ///< maps a string identifier for a timer to an index in the start vector.
};
extern class GeneralStuff GS; ///< physical constants, utility functions.



//
// This contains the definition of a grid cell, and the access
// functions for its data, position, neighbours, etc.
//
#include "grid/cell_interface.h"
#ifdef COUNT_ENERGETICS
extern struct energetics *GLOBAL_CE;
 ///< for tracking rates in microphysics/raytracing.
#endif
extern class cell_interface CI;
  ///< Global Instance of cell interface class.


#ifndef GRIDV2

///
/// Abstract Base Class to define interface for the Grid data and
/// member functions.
///
class GridBaseClass {
  public:
   virtual ~GridBaseClass() {}
   /** \brief copies contents of one cell into another. */
   virtual void CopyCell(const cell *, ///< Pointer to original.
       cell *       ///< Pointer to copy.
       )=0;
   /** \brief Prints contents of cell. */
   virtual void PrintCell(const cell *)=0;
   /**\brief returns a pointer to a neighbouring cell in the given direction.
    * */
   virtual cell * NextPt(const cell *,         ///< Current cell.
       const enum direction ///< direction of neighbour.
       )=0;
   /**\brief returns a pointer to the next cell.*/
   virtual cell * NextPt(const cell *)=0;
   /** \brief Returns the opposite direction to the one passed in to the function. */
   virtual enum direction OppDir(enum direction )=0;
   /** \brief Like nextPt(cell,dir), but in reverse direction. */
   virtual class cell* PrevPt(const class cell*, ///< Current Point.
            enum direction     ///< Direction to go in.
            )=0;
   /** \brief Runs through ghost boundary cells and does the appropriate time update on them.*/
   virtual int TimeUpdateExternalBCs(const int, ///< Current step number in the timestep.
             const int  ///< Maximum step number in timestep.
             )=0;
   /** \brief Runs through boundary cells which are grid cells and does the appropriate time update on them.*/
   virtual int TimeUpdateInternalBCs(const int, ///< Current step number in the timestep.
             const int  ///< Maximum step number in timestep.
             )=0;
   virtual double DX() const =0; ///< Returns dx.
   virtual double DA() const =0; ///< Returns dA.
   virtual double DV() const =0; ///< Returns dV.
   virtual int Ndim() const =0; ///< Returns dimensionality of grid.
   virtual int Nvar() const =0; ///< Returns length of state vectors.
   virtual double  Xmin(enum axes) const =0; ///< Returns x,y,z lower bounds.
   virtual double  Xmax(enum axes) const =0; ///< Returns x,y,z upper bounds.
   virtual double Range(enum axes) const =0; ///< Returns x,y,z range.
   virtual int iXmin(enum axes) const =0; ///< Returns x,y,z lower bounds in cell integer coordinates.
   virtual int iXmax(enum axes) const =0; ///< Returns x,y,z upper bounds in cell integer coordinates.
   virtual int iRange(enum axes) const =0; ///< Returns x,y,z range in cell integer coordinates.
#ifdef GEOMETRIC_GRID
   virtual int SIM_iXmin(enum axes)  const =0; ///< Returns GLOBAL x,y,z lower bounds in cell integer coordinates.
   virtual int SIM_iXmax(enum axes)  const =0; ///< Returns GLOBAL x,y,z upper bounds in cell integer coordinates.
   virtual int SIM_iRange(enum axes) const =0; ///< Returns GLOBAL x,y,z range in cell integer coordinates.
#endif // GEOMETRIC_GRID
   /** \brief Returns pointer to the first grid point. */
   virtual cell * FirstPt()=0;
   /** \brief Returns pointer to the last grid point.  */
   virtual cell * LastPt()=0;
   virtual int SetupBCs(int,   ///< Depth of Boundary cells, 1,2,etc.
      string ///< string containing info on types of BCs on all sides.
      )=0;
#ifdef PLLEL_RT
   /** \brief Setup lists of processors to receive data from and send data to, 
    * and setup extra boundaries at corners. */
   virtual int Setup_RT_Boundaries(const int   ///< id of source.
           )=0;
   /** \brief Receive all optical depths for boundaries closer to source. */
   virtual int Receive_RT_Boundaries(const int ///< source id
             )=0;
   /** \bried Send all optical depths for boundaries to domains further from source. */
   virtual int Send_RT_Boundaries(const int ///< source id
          )=0;
#endif // PLLEL_RT

#ifdef GEOMETRIC_GRID
   ///
   /// Calculate distance between two points, where the two position
   /// are interpreted in the appropriate geometry.
   /// This function takes input in physical units, and outputs in 
   /// physical units.
   ///
   virtual double distance(const double *, ///< position 1 (physical)
         const double *  ///< position 2 (physical)
         )=0;

   ///
   /// Calculate distance between two points, where the two position
   /// are interpreted in the appropriate geometry.
   /// This function takes input in code integer units, and outputs in
   /// integer units (but obviously the answer is not an integer).
   ///
   virtual double idistance(const int *, ///< position 1 (integer)
          const int *  ///< position 2 (integer)
          )=0;
   
   ///
   /// Calculate distance between two cell--centres (will be between
   /// centre-of-volume of cells if non-cartesian geometry).
   /// Result returned in physical units (e.g. centimetres).
   ///
   virtual double distance_cell2cell(const cell *, ///< cell 1
             const cell *  ///< cell 2
             )=0;

   ///
   /// Calculate distance between two cell--centres (will be between
   /// centre-of-volume of cells if non-cartesian geometry).
   /// Result returned in grid--integer units (one cell has a diameter
   /// two units).
   ///
   virtual double idistance_cell2cell(const cell *, ///< cell 1
              const cell *  ///< cell 2
              )=0;
   
   ///
   /// Calculate distance between a cell-vertex and a cell--centres
   /// (will be between centre-of-volume of cells if non-cartesian
   /// geometry).  Here both input and output are physical units.
   ///
   virtual double distance_vertex2cell(const double *, ///< vertex (physical)
               const cell *    ///< cell
               )=0;

   ///
   /// Calculate distance between a cell-vertex and a cell--centres
   /// (will be between centre-of-volume of cells if non-cartesian
   /// geometry).  Here both input and output are code-integer units.
   ///
   virtual double idistance_vertex2cell(const int *, ///< vertex (integer)
          const cell * ///< cell
          )=0;

   ///
   /// As idistance_vertex2cell(int,cell) but for a single component
   /// of the position vector, and not the absolute value.  It returns
   /// the *cell* coordinate minus the *vertex* coordinate.
   ///
   virtual double idifference_vertex2cell(const int *,  ///< vertex (integer)
            const cell *, ///< cell
            const axes    ///< Axis to calculate.
            )=0;
   ///
   /// As idifference_vertex2cell(int,cell,axis) but for the coordinate
   /// difference between two cell positions along a given axis.
   /// It returns *cell2* coordinate minus *cell1* coordinate.
   ///
   virtual double idifference_cell2cell(
              const cell *, ///< cell 1
              const cell *, ///< cell 2
              const axes    ///< Axis.
              )=0;
#endif // GEOMETRIC_GRID
  protected:
};

/** \brief Global Pointer to the Computational Grid.
 * 
 * This is a pointer to an abstract base class, or interface class,
 * which should be appropriate for any sort of finite volume grid,
 * providing enough access functions to do whatever I need to.  It
 * can be initialised to any derived classes, as long as they implement
 * all the functions of the base class.  For example my UniformGrid class
 * has been modified to be derived from the base class, and the 
 * parallel implementation of the uniform grid will be further derived
 * from this.
 * */
extern class GridBaseClass *grid;
/*****************************************************************/
#endif // GRIDV2


/**********************************************************************/
/************************ MULTI-PROCESS COMMS *************************/
//
// I want the communication to be independent of MPI or any other parallelisation library,
// so I setup an abstract base class in comms.h, and implement it in either comm_mpi or
// comm_files, and it should abstract away all the details of how the communication happens.
//
#ifdef PARALLEL
#include "comms/comms.h"

#if   defined USE_MPI
#include "comms/comm_mpi.h"
#elif defined USE_FILE_COMMS
#include "comms/comm_files.h"
#else
#error "MUST DEFINE EITHER USE_MPI or USE_FILE_COMMS"
#endif

extern class comms_base *COMM;
#endif // PARALLEL
/************************ MULTI-PROCESS COMMS *************************/
/**********************************************************************/



/************************* GRID METHODS **************************/
/** \brief Abstract base class for finite volume grids.
 * */
class IntegratorBaseFV
{
  public:
   virtual ~IntegratorBaseFV() {}
   virtual int Init(string, ///< Name of input file.
        int,    ///< Type of File (1=ASCII, 2=FITS, 3=fitstable, 4=fits and ascii, ...).
        int,    ///< Number of command-line arguments.
        string * ///< Pointer to array of command-line arguments.
        ) =0;
   virtual int Time_Int() =0;
   virtual int Finalise() =0;
};
extern class IntegratorBaseFV *integrator;
/************************* GRID METHODS **************************/


/************************ EQN SOLVER ***************************/
/*****************************************************************/



/************************ MICROPHYSICS ***************************/
#include "microphysics/microphysics_base.h"
extern class MicroPhysicsBase *MP;

/************************ MICROPHYSICS ***************************/


/*****************************************************************/
/************************** RAY TRACER ***************************/
/** \brief Pure virtual ray-tracer base class. 
 *
 * This provides the interface to the raytracer, regardless of the 
 * implementation.  So far these are the functions I'll need for 
 * a uniform grid with a single source, so I may need to add more
 * later.
 * */
class RayTracingBase {
 public:
  virtual ~RayTracingBase() {} ///< non-virtual destructor.

  ///
  /// \brief Adds a source to the list of sources to trace.
  /// 
  /// Returns the sources id, which starts at zero and increases 
  /// linearly.  So if we have 10 sources, and add another source,
  /// its id will be 10.  Note the ID is contained in the rad_src_info
  /// struct though, so each particular class may not have consecutively 
  /// numbered sources.
  ///
  virtual int Add_Source(struct rad_src_info * ///< ptr to source info.
       )=0;
  ///
  /// Processes a single source's effect on the grid over a timestep.
  ///
  virtual int RayTrace_SingleSource(const int, ///< Source id
                                    const double, ///< Timestep
                                    const double  ///< EOS Gamma.
                                    )=0;
  ///
  /// Just calculate the column densities required for RT.
  ///
  virtual int RayTrace_Column_Density(
                const int,    ///< Source id
                const double, ///< Timestep
                const double  ///< EOS Gamma.
                )=0;

  /** \brief Prints list of sources with id, location, strength. */
  virtual void Print_SourceList()=0;
  /** \brief Returns the number of sources to track. */
  virtual int NSources()=0;

    ///
  /// Returns whether we are doing an implicit (==0) or an explicit (==1)
  /// integration of the raytracing/microphysics.
  ///
  virtual int type_of_RT_integration()=0;

  ///
  /// This sets the number of ionising and UV heating sources of radiation,
  /// and makes sure the rt_source_data structs are populated correctly.
  /// It can be used to change the radiation sources if e.g. the luminosity
  /// changes over time.
  ///
  virtual void update_RT_source_properties(
                  const struct rad_src_info * ///< ptr to source info.
                  )=0;

  /// Returns the number of ionising sources
  virtual int N_ionising_sources()=0;

  /// Returns the number of UV heating sources
  virtual int N_heating_sources()=0;

  ///
  /// This function copies the ionising source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_ionising_src_list(
                std::vector<struct rt_source_data> & ///< list of data for ionising sources
                )=0;

  ///
  /// This function copies the UV-heating source data into two
  /// vectors of structs which are returned by reference.
  /// If rt-testing flags are set it will check that the input vector matches
  /// the number of sources to add to the list, but otherwise there is no
  /// checking.
  ///
  virtual int populate_UVheating_src_list(
                std::vector<struct rt_source_data> & ///< list of data for UV-heating sources
                )=0;

};

extern class RayTracingBase *RT; ///< Raytracer for all radiation sources
/************************** RAY TRACER ***************************/
/*****************************************************************/



/************************* ERROR REPORTING ***********************/
/** \brief Global Class For Writing Error messages and Reports.
 * This determines where to write different messages to.
 * */
class reporting {
  public:
   reporting(); ///< Default Constructor, write to std[out/err].
   ~reporting(); ///< Destructor.
   
   /** \brief Redirects stdout/stderr to files in the path specified.*/
   int redirect(const string & ///< Location of files to write reporting to.
    );
   
   /** \brief This exits the code cleanly if I know I have an error.
    *
    * This exits from the code, printing an error message and the  
    * offending value.
    * I would like to have this in global.cc, but it is a template
    * function, and the export keyword is not implemented in any
    * compilers, so I have to have the definitions in the header.
    */
  template<class T> inline void error(string msg, ///< Error Message.
               T err       ///< Value which is wrong.
               )
  {
    cerr <<msg<<"\t error code: "<<err<<" ...exiting.\n";
    cerr <<"\t timestep="<<SimPM.timestep<<" and simtime="<<SimPM.simtime<<endl;
    flush(cout);
    // set maxtime and then Finalise will output data in the current state.
    SimPM.maxtime=true;
    if (integrator) {
#ifdef SERIAL
      //
      // would like to output data and quit, but if the error is in
      // the output routine, then this puts us in an infinite loop!
      //
      //if (grid) integrator->Finalise();
#endif // SERIAL
      delete integrator; integrator=0;
    }
    if (grid) {delete grid; grid=0;}
    if (MP) {delete MP; MP=0;}
    if (RT) {delete RT; RT=0;}
#ifdef PARALLEL
    if (COMM) {
      COMM->abort();
    }
#endif // PARALLEL
    exit(1);
  }
   
   /** \brief Prints a warning message.
    * If an Value is not what was expected, but I don't want
    * to bug out or call it an 'error', then this function displays a
    * message, and the expected and received values.
    * */
   template<class T1, class T2> inline void warning(string msg,  ///< Message.
                T1 expected,///< Expected Value
                T2 found    ///< Received Value
                )
     {
       cout<<"WARNING: "<<msg<<"\t Expected: "<<expected<<" but got "<<found<<endl;
       return;
     }
   
   /** \brief Error if actual value is less than Expected Value
    * T1 and T2 must be compatible types for testing 'less-than' operation.
    * */
   template<class T1, class T2>inline void errorLT(string msg, ///< Message.
               T1 exptd,  ///< Expected Value
               T2 recvd   ///< Actual Value
               )
     {
       if(recvd<exptd) {
   cout<<"ERROR DETECTED: "<<msg<<"\t expected less than "<<exptd;
   cout <<" but got "<<recvd<<endl; 
   error(msg,recvd);
       }
       return;
     }
   
   /** \brief Tests if two values are the same, and error if not.
    * If the expected value is not equal to the actual value, then print
    * an error message and exit. */
   inline void errorTest(string msg, ///< Error Message.
       int exptd, ///< Expected Value.
       int recvd  ///< Actual Value.
       )
     {
       if(exptd!=recvd) {
   cerr<<"ERROR DETECTED: "<<msg<<"\t expected "<<exptd<<" but got "<<recvd<<endl; 
   error(msg,recvd);
       }
     }
   
   /** \brief Print out a vector. */
   template<class T> inline void printVec(string msg, ///< Name of Vector
            T* vec,   ///<pointer to vector
            int nd    ///< length of vector.
            )
     {
       cout <<"Vector "<<msg<<" : [";
       for (int i=0;i<nd-1;i++) cout <<vec[i]<<", ";
       cout <<vec[nd-1]<<" ]\n";
       return;
     }
   
 private:
   ofstream errmsg, infomsg;//, iomsg;
   streambuf *saved_buffer_cout, *saved_buffer_cerr;
};
extern class reporting rep;
/************************* ERROR REPORTING ***********************/


/************************* MEMORY MANAGEMENT ***********************/
/** \brief This class contains generic functions for dynamically allocating
 * and freeing memory.  Hopefully everywhere in the code will use this and
 * not have to write it all out every time.
 */
class memory_management {
public:
  memory_management() {}
  ~memory_management() {}
  /** \brief Allocates memory for an array of any valid type.
   *
   * If pointer passed in is not null (=0), a warning is displayed and
   * the pointer is returned unchanged.  If it is null, then the 
   * function attempts to initialise it, and if it fails, it exits the
   * code, since memory allocation failures are generally non-recoverable
   * (in the sense that if I ask for memory, I need it for something...).
   */
  ///
  /// Turns out this isn't particularly useful since I've read in various 
  /// places that allocation requests NEVER fail unless they eat up
  /// all the RAM and SWAP, and maybe even more.  OS's are designed to
  /// at least try to deal with any memory requests.
  ///
  template<class T> T * myalloc(T *ptr, ///< uninitialised null pointer.
        const long int n_el ///< number elements to initialise.
#ifdef TESTING
        //, const string msg ///< message to display if allocation fails.
#endif
        )
    {
#ifndef CHECK_NEW_EXCEP_ON
#error "Only use memory_management with try/catch exceptions."
#endif
      //if (ptr) {
      //#ifdef TESTING
      //  cerr <<"mem_alloc() from: "<<msg<<endl;
      //#endif
      //  cerr <<"mem_alloc() pointer is initialised! returning unchanged.\n";
      //}
      //else {
  try {ptr = new T [n_el];}
  catch (std::bad_alloc) {
#ifdef TESTING
    //    cerr <<"mem_alloc() from: "<<msg<<endl;
#endif
    cerr <<"mem_alloc() pointer initialisation failed.\n";
    rep.error("mem_alloc() pointer initialisation failed.",ptr);
  }
  //}
      return ptr;
    }

  /** \brief Delete a pointer to an array.
   *
   * If the pointer is already null (=0), print a warning and don't
   * try to delete it again!  Otherwise delete it and set it to zero.
   */
  template<class T> T * myfree(T *ptr
#ifdef TESTING
             //, const string msg ///< message to display if allocation fails.
#endif
             )
    {
      //
      // if ptr is already null, can't free anything, so just return.
      //
      if (!ptr) {
#ifdef TESTING
  //cerr <<"mem_myfree() from: "<<msg<<endl;
  //cerr <<"mem_myfree() pointer is already 0!\n";
#endif
      }
      else {
  delete [] ptr; ptr=0;
      }
      return 0;
    }
  /** \brief Delete a pointer to single object (no [] in delete).
   *
   * If the pointer is already null (=0), print a warning and don't
   * try to delete it again!  Otherwise delete it and set it to zero.
   */
  template<class T> T * myfree_single(T *ptr
#ifdef TESTING
              //, const string msg ///< message to display if allocation fails.
#endif
              )
    {
      //
      // if ptr is already null, can't free anything, so just return.
      //
      if (!ptr) {
#ifdef TESTING
  //cerr <<"mem_myfree_single() from: "<<msg<<endl;
  //cerr <<"mem_myfree_single() pointer is already 0!\n";
#endif
      }
      else {
#ifdef TESTING
  //cout <<"mem_myfree_single() from: "<<msg<<endl;
  //cout <<"\tptr="<<ptr<<endl;
#endif
  delete ptr; ptr=0;
      }
      return ptr;
    }
};

extern class memory_management mem;
/************************* MEMORY MANAGEMENT ***********************/
/*******************************************************************/



#ifdef TESTING
/** \brief class for debugging the code, tracking energy/momentum to make
 * sure it is conserved, etc.
 * 
 * Also contains a cell pointer, so we know what cell we are working on at
 * all times.
 * The grid and MP pointers are globally accessible pointers to the 
 * grid data and microphysics classes, so I don't need them here, but I 
 * do need pointers to the solver and the integrator, so they are 
 * included, and are set whenever the respective classes are initialised.
 * */
class DebugParams {
  public:
   DebugParams();
   ~DebugParams();
   class cell *c; ///< Pointer that can be used to track which cell we are at.
   //class BaseFVSolver *solver; ///< pointer to solver.
   double initERG;   ///< Initial total energy on the grid (not including ghost cells).
   double initMMX;   ///< Initial X-momentum on the grid (not including ghost cells).
   double initMMY;   ///< Initial Y-momentum on the grid (not including ghost cells).
   double initMMZ;   ///< Initial Z-momentum on the grid (not including ghost cells).
   double ergTotChange; ///< For tracking energy entering/leaving domain.
   double mmxTotChange; ///< For tracking x-momentum entering/leaving domain.
   double mmyTotChange; ///< For tracking y-momentum entering/leaving domain.
   double mmzTotChange; ///< For tracking z-momentum entering/leaving domain.
   double vec1[MAX_NVAR];
   double vec2[MAX_NVAR];
   double vec3[MAX_NVAR];
   double vec4[MAX_NVAR];
};
extern class DebugParams dp;

/** \brief Command Line Interface for interactive debugging of code.
 * 
 * Adapted from base code Andy gave me.
 */
class CommandLineInterface {
  public:
   CommandLineInterface(); ///< Do-nothing Constructor.
   ~CommandLineInterface(); ///< Do-nothing Destructor.
   /** \brief This function calls up the CLI for interactive debugging. */
   void console(char * ///< Optional text for prompt.
    );
   void auto_console(char * ///< Optional text for prompt.
         );
  private:
   /** \brief Reads a string, and return a pointer to it.  Returns NULL on EOF. */
   char *rl_gets(char *,  ///< pointer to line of text.
                 char *   ///< text inputted to prompt.
                 ); 
   /** \brief Executes any commands the CLI recognises in the string passed in. */
   int execute(char *);
   void cmd1(); ///< An example command, with no arguments.
   void cmd2(char *); ///< another example command, with arguments.
   void bigcmd(); ///< A further example command.
   void print_cell(); ///< prints cell which dp.c points to.(DebugParams)
   void next_point(const string); ///< moves to cell in direction given (XN,XP,YN,YP,ZN,ZP)
   void fpt(); ///< moves dp.c to the first point.
   void lpt(); ///< moved dp.c to the last point.
   void end_of_col(const string); ///< moves dp.c to the edge of the grid in direction specified.
   void print_flux(); ///< prints flux data from intercell flux function.
   enum direction parse_dir(const string); ///< given a string of XN,YN,etc, return direction.
};
extern class CommandLineInterface commandline;
#endif //TESTING

#endif // GLOBAL_H
   

/** \page todo Stuff I Need to Work On
 * 
 * \section Introduction
 * This page contains a list of things I need to work on in the code.  They 
 * are not bugs as such.  Every comment in the code with a "\todo" tag  will
 * show up on this page, as well as anything else I put on the page.
 * 
 * */

/** \mainpage 
 * \author Jonathan Mackey
 * \section Introduction
 * This is my code.
 * */

/** \page userguide Users Guide
 * \section Introduction
 * This is intended as a guide to using the code.
 * 
 * \section compilation Compilation
 * There are two makefiles... <tt>Makefile.serial</tt> and <tt>Makefile.parallel</tt>
 * which should be modified according to the compiler and architecture
 * you are using.  <tt>Makefile.serial</tt> compiles the serial (single processor)
 * version of the code, and produces an executable called 
 * <tt>main_serial</tt>, whereas <tt>Makefile.parallel</tt> comiles the 
 * parallel version, using some distribution of the message passing
 * interface (MPI), producing an executable called <tt>main_parallel</tt>.
 * 
 * \section tracers Tracer Variables.
 * Tracer Variables are any variable that obey a simple linear advection
 * equation in the hydrodynamics.  For example the fractional abundance of 
 * electrons is passively transported by the hydrodynamics, and then 
 * updated in the microphysics routines.  The class SimParams has three
 * members relating to the tracers:
 *  - SimParams::ntracer  This is the integer number of tracer variables
 *  - SimParams::ftr      This is the integer index of the first tracer variable.  The
 * tracers are the last <tt>ntracer</tt> elements of the primitive vector, so they are
 * indexed as <tt>[ftr, ftr+1, ftr+2, ..., ftr+ntracer-1]</tt>.
 *  - SimParams::trtype   This is a string containing information on what the
 * tracers are for.  It is split into 5 character units.  The first unit is
 * the type of chemistry we are using (ex. atomic hydrogen only), and the
 * remaining ntracer units describe each tracer.  The code parses this string
 * during setup to determine what to do with the tracer variables.  Note that 
 * the numbers must match for the code to run -- there must be exactly the 
 * right number of units in trtype to give a type to each and every tracer
 * variable.
 * 
 * The first unit of SimParams::trtype is the type of tracer we are using,
 * so if we have chemistry we need to know what kind of chemistry, and if 
 * they are just colour tracers we need to know that too.  So here is the list
 * of known values for the first string unit:
 *  - <tt>"color"</tt> corresponds to only passive colour tracers;
 *  - <tt>"chAH_"</tt> corresponds to atomic hydrogen chemistry only;
 *  - <tt>"chAMH"</tt> means atomic and molecular hydrogen chemistry;
 *  - Others will be added to this list in due course.
 * 
 * The remaining units of SimParams::trtype describe each tracer individually,
 * and here is the list so far of what they do:
 *  - <tt>"trace"</tt> non-chemistry tracer, so just advected as a colour tracer.
 *  - <tt>"e-___"</tt> electron fractional abundance.
 *  - <tt>"HI___"</tt> neutral hydrogen fractional abundance.
 *  - <tt>"HII__"</tt> ionized hydrogen fractional abundance.
 *  - <tt>"H2___"</tt> molecular hydrogen fractional abundance.
 *  - <tt>"HEI__"</tt> neutral helium fractional abundance.
 *  - <tt>"HEII_"</tt> singly ionized helium fractional abundance.
 *  - <tt>"HEIII"</tt> doubly ionized helium fractional abundance.
 *  - <tt>any others?</tt>
 * Note that tracers should be dimensionless quantities that scale with the
 * density.
 *
 * \section mpint MicroPhysics Internal Variables
 * Microphysics has its own state vector, without the extra stuff that it
 * doesn't need, elements as follows:
 *  - 0: "n_h" : hydrogen number density, obtained from p_in[RO].
 *  - 1: "Eint": Internal Energy per unit volume 
 *  (obtained from p_in[RO]/(gamma-1) ).
 *  - 2: "e-"  : electron fraction n(e-)/n(H).  From here on all variables are
 *  optional.
 *  - 3-: "ION": number fraction of ion X relative to total number density of
 *  that element.
 * */

/** \page algorithms Code Algorithms
 * 
 * \section glm GLM Method
 * The GLM method of divergence cleaning is described in Dedner, Kemm, 
 * Kroner, et al., 2002, JCP, 175, 645, and in Dedner's thesis which I 
 * downloaded from his webpage (section 8).  I am using the mixed-glm method
 * with some modifications.
 * 
 * \subsection chyp Hyperbolic Wavespeed
 * Dedner advocates using the hyperbolic wavespeed
 * as \f$c_{\mbox{hyp}} = \mbox{CFL}\times\delta x/\delta t\f$, but in 3d I
 * found that this doesn't work.  I get checkerboard instabilities developing
 * in the field (for CFL>0.35).
 * The hyperbolic wavespeed shouldn't be a problem, because if
 * the CFL number is set to give a stable scheme then it should remain stable
 * when I add in a wave with the above speed.
 * Anyway, if I modify the hyperbolic speed to be 
 * \f$c_{\mbox{hyp}} = \mbox{CFL}\times\delta x/\delta t/N_{\mbox{dim}}\f$, 
 * then it recovers most of its stability range again.  This is a hack, and
 * I'm not at all sure why this works, or even why it is a good thing to do.
 * In fact I have taken out the ndim division because I don't like it.  I suspect the
 * reduced range of stability is a problem with the method and not of my
 * implementation of it.
 * 
 * \subsection cpar Parabolic Speed
 * The parabolic component of the mixed glm method has a parabolic speed,
 * which is basically a diffusion constant, \f$c_p\f$, with units of distance
 * over square root of time.  Dedner suggests using a different constant
 * instead of setting \f$c_p\f$, where the new constant is 
 * \f$c_r \equiv c_p^2/c_h\f$, but this has units of length, so I multiply
 * it by the grid spacing \f$\delta x\f$ (this is what Andy does).  In his
 * paper he says using a value of 0.18 for \f$c_r\f$ is good, but he uses the
 * reciprocal in his thesis, and I found that the reciprocal reproduces the
 * figures in the paper whereas 0.18 doesn't.  I experimented a bit and 
 * settled on \f$c_r = 4.0\delta x\f$ as a good value.
 * In his thesis (p.121 eq.8.22) he suggests using cr as a length related 
 * to the box-size.  I tried his formula (given that it is the reciprocal of
 * what I am using), and found it gave a larger value of cr for my test problem,
 * but I don't think it scales the way it should.
 * 
 * \subsection wrap Conclusions
 *  - It looks like this Dedner method is not quite as good as claimed, but works
 * well for low CFL numbers (<0.5 for 2d and <0.35 for 3d).  This is 
 * restrictive, but maybe not prohibitively so.  
 *  - When there is an instability,
 * it starts out with a checkerboard pattern in Psi, which then appears in
 * divB, and finally appears in the field components, or else crazy pressure
 * and density values.
 *  - It would be nice to have a better scheme, but I'm not sure there is one.
 * The constrained Transport (CT) schemes seem to be
 * better, although they are very complex to look at. But they also have problems,
 * in that you are again post-processing the magnetic field after the flux updates,
 * purely in order to give a divergence free field.  You keep divB=0, but you
 * lose in other ways, as far as I can tell.
 *  - I can reduce the Psi speed to get better stability.  Not sure why I 
 * would need to do this, but I'd say it's more to do with the operator
 * splitting and coupling between the parabolic and hyperbolic components
 * than it is to do with the actual wavespeed being too fast.
 *  - Dedner says in his thesis (8.5.1, p.124) that increasing cr is better, 
 * but leads to an unstable solution, whereas he says earlier that things 
 * should be unconditionally stable for any cr.  So consequently for me, 
 * decreasing cr should decrease stability.  Things are clearly not quite as
 * good as he makes out in their paper.
 * */

