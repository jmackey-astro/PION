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
/// - 2013.09.06 JM: Added difference_vertex2cell() functions to grid
///    base class.
/// - 2015.01.08 JM: moved grid base class to grid/grid_base_class.h.
/// - 2015.01.12 JM: moved reporting to tools/reporting.h.
///    Moved commandline interface to tools/command_line_interface.h.
///    Moved memory management to tools/mem_manage.h.
/// - 2015.01.26 JM: removed IntegratorBaseFV (no point to this base
///    class, as far as I can see (famous last words...))
///    Removed ParallelParams class (to sim_control_MPI.h), and COMM.

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


#include "sim_constants.h"
#include "sim_params.h"

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

//#define MAX_NTR   5 ///< Maximum number of tracer variables (feel free to expand later!).

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
#define RT_EFFECT_HHE_MFQ    9   ///< Frank&Mellema H+He p-ion scheme

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

