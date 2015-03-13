///
/// \file constants.h
/// \author Jonathan Mackey
/// \date 2013.08.19
///
/// Class with inline functions for many physics and astronomy
/// constants.  All constants are in cgs units unless stated
/// otherwise.
///
///
/// Modifications:
/// - getting it written: mods up until 2013.08.20
/// - 2014.01.20 JM: Added more significant digits
/// - 2015.03.10 JM: added equalD() function.
/// - 2015.03.11 JM: added more constants from global.h


#ifndef CONSTANTS_H
#define CONSTANTS_H


//
// It's a good idea to include cmath here, since nearly every file needs it!
//
#ifndef INTEL
#include <cmath>     // Header file from gcc
#else
#include <mathimf.h> // Header file from Intel Compiler
#endif

class constants {
  public:
  constants();
  ~constants();
  //
  // Mathematics
  //
  inline double ln10() {return 2.3025850929940459;}
  inline double pi()     {return 3.14159265358979324;}
  //
  // Physics
  //
  inline double c()      {return 2.99792458e+10;}
  inline double kB()     {return 1.3806488e-16;}
  inline double h()      {return 6.62606957e-27;}
  inline double StefanBoltzmannConst() {return 5.670373e-5;}
  inline double m_p()    {return 1.672621777e-24;}
  //
  // Atomic data
  //
  inline double Eth_H()  {return 13.59844*eV();}
  inline double NuTh_H() {return 3.288e15;}
  //
  // Unit conversion
  //
  inline double K_per_eV() {return 1.1604519e4;}
  inline double year()     {return 3.1558150e7;}
  inline double eV()       {return 1.602176565e-12;}
  //
  // Astronomy
  //
  inline double Msun()   {return 1.9891e33;}
  inline double Lsun()   {return 3.839e33;}
  inline double Rsun()   {return 6.955e10;}
  inline double AU()     {return 1.49597870700e+13;}
  inline double parsec() {return 3.0856776e18;}
  
  ///
  /// tests if two doubles are equal to an accuracy of 1e-8.
  ///
  bool equalD(
        const double, 
        const double
        );

};

extern class constants pconst;

#define ION_DUST -1
#define ION_H_N  0
#define ION_HE_N 1
#define ION_HE_P 2

#define EL_H  0
#define EL_HE 1



//-------------------------------------------------------
// Defines for different sized numbers.
//-------------------------------------------------------
#define SMALLVALUE       1.0e-12
#define MACHINEACCURACY  5.e-16
#define TINYVALUE        1.0e-100
#define VERY_TINY_VALUE  1.0e-200
#define VERY_LARGE_VALUE 1.0e100
#define HUGEVALUE        1.0e+200
#define ONE_PLUS_EPS     (1.0+SMALLVALUE)
#define ONE_MINUS_EPS    (1.0-SMALLVALUE)




//-------------------------------------------------------
// Equations to solve
//--------------------------------------------------------
/// Equation type: 1=Euler equations.
#define EQEUL 1
/// Equation type: 2=Ideal MHD equations (Falle,
/// Komissarov, Joarder, 1998). (ONLY GOOD FOR 1D!)
#define EQMHD 2
/// Equation type: 3=MHD equations with GLM method for
/// divergence cleaning (Dedner etal, 2002).
#define EQGLM 3
/// Equation type: 4=MHD equations with Field-CD method
/// for ensuring divB=0 (Toth,2000). (DON'T USE!!!)
#define EQFCD 4  
/// Equation type: 5=Isothermal Hydrodynamics (UNTESTED)
#define EQEUL_ISO 5
/// Equation type: 6=Isothermal MHD (UNTESTED)
#define EQMHD_ISO 6
/// Equation type: 7=Isothermal MHD with Dedner et al. GLM
// method  for Divergence Cleaning (UNTESTED)
#define EQGLM_ISO 7
/// Equation type: 8=Isothermal MHD with Field-CD method
/// for ensuring divB=0 (Toth,2000). (UNTESTED)
#define EQFCD_ISO 8
/// Equation type: 9=Adiabatic Hydro, integrating INTERNAL and not
/// total energy.  (UNTESTED)
#define EQEUL_EINT 9


//--------------------------------------------------------
// Coordinate systems: Cartesian, cylindrical, spherical.
//--------------------------------------------------------
/// Coordinates: Cartesian coordinate system.
#define COORD_CRT 1 
/// Coordinates: Cylindrical coordinate system.
#define COORD_CYL 2
/// Coordinates: Spherical coordinate system.
#define COORD_SPH 3

//--------------------------------------------------------
///
/// Definitions for cylindrical axes and directions w.r.t. cartesian
/// ones.
///
//--------------------------------------------------------
#define Zcyl  XX
#define Rcyl  YY
#define Tcyl  ZZ
#define ZNcyl  XN
#define ZPcyl  XP
#define RNcyl  YN
#define RPcyl  YP
#define TNcyl  ZN
#define TPcyl  ZP

//--------------------------------------------------------
///
/// Definitions for Spherical coordinate axes and directions w.r.t. cartesian
/// ones. [R,T,P] stands for [r,theta,phi]
///
//--------------------------------------------------------
#define Rsph  XX
#define Tsph  YY
#define Psph  ZZ
#define RNsph  XN
#define RPsph  XP
#define TNsph  YN
#define TPsph  YP
#define PNsph  ZN
#define PPsph  ZP

//--------------------------------------------------------
///
/// Type of Solver:
///
//--------------------------------------------------------
#define FLUX_LF       0
#define FLUX_RSlinear 1
#define FLUX_RSexact  2
#define FLUX_RShybrid 3
#define FLUX_RSroe    4
#define FLUX_RSroe_pv 5
#define FLUX_FVS      6


//--------------------------------------------------------
/// \brief Primitive Variables Enum 
/// The code requires that the velocities follow each other in order,
/// and the same for the Magnetic Field.  I'm not sure if it will work 
/// without them being in that order, but I'm writing the code assuming 
/// they are in order.
///
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
//--------------------------------------------------------

//--------------------------------------------------------
/// \brief Conserved Variables Enum
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
//--------------------------------------------------------

//--------------------------------------------------------
/// \brief Enum for each Spatial Dimension.
enum axes {
    XX=0, ///< X-Direction.
    YY=1, ///< Y-Direction.
    ZZ=2  ///< Z-Direction.
};
//--------------------------------------------------------


///
/// \brief Direction Vectors enum. 
/// This gets used for pointing to elements of arrays,
/// so the six real directions must comprise the numbers
/// zero to five.  They must also be arranged so that e.g.
/// XP = 2*axes.XX +1, XN = 2*axes.XX, and similarly for 
/// Y and Z directions.
///
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




#endif // CONSTANTS_H


