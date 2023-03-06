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
#include <cmath>  // Header file from gcc
#else
#include <mathimf.h>  // Header file from Intel Compiler
#endif

///
/// Physical constants are from CODATA 2014, converted to CGS.
/// http://dx.doi.org/10.5281/zenodo.22826
/// Masses of elements from http://www.ciaaw.org/atomic-weights.htm
///
class constants {
public:
  constants();
  ~constants();

  //
  // Mathematics
  //
  inline double ln10() { return 2.3025850929940459; }
  inline double pi() { return 3.14159265358979324; }
  inline double halfpi() { return 0.5 * pi(); }
  inline double root4pi() { return 3.5449077018110318; }
  inline double sqrt2() { return 1.4142135623730950; }
  //
  // Physics
  //
  inline double c() { return 2.99792458e+10; }
  inline double kB() { return 1.38064852e-16; }
  inline double h() { return 6.62606957e-27; }
  inline double StefanBoltzmannConst() { return 5.670367e-5; }
  inline double eV() { return 1.6021766208e-12; }
  inline double ThomsonXSection() { return 6.6524587321e-25; }
  ///
  /// Unified atomic mass unit.
  ///
  inline double uamu() { return 1.660539040e-24; }
  ///
  /// Electron mass (CODATA 2018) in grams
  ///
  inline double m_e() { return 9.1093837015e-28; }
  ///
  /// proton mass (CODATA 2014) in grams
  ///
  inline double m_p() { return 1.672621898e-24; }
  ///
  /// Hydrogen mass (IUPAC value is 1.0080+/-0.0002 * uamu).  in grams
  ///
  inline double m_H() { return 1.6738e-24; }
  ///
  /// Helium mass (CIAWW online table) in grams
  ///
  inline double m_He() { return 6.6464768e-24; }

  ///
  /// Carbon mass (CIAWW online table) in grams
  /// 0.5*(12.0096 + 12.0116) uamu
  ///
  inline double m_C() { return 1.994374e-23; }

  ///
  /// Nitrogen mass 0.5*(14.00643 + 14.00728) uamu in grams
  /// (CIAWW online table)
  ///
  inline double m_N() { return 2.325892e-23; }

  ///
  /// Oxygen mass (CIAWW online table) in grams
  /// 0.5*(15.99903 + 15.99977) uamu
  ///
  inline double m_O() { return 2.6567628e-23; }

  /// Gravitational constant (CODATA 2014)
  ///
  inline double G() { return 6.67408e-8; }
  ///
  /// Ionization potential of H (eV)
  ///
  inline double Eth_H() { return 2.178710264e-11; }
  ///
  /// Ionization potential of H, in Hertz.
  ///
  inline double NuTh_H() { return 3.28808819e+15; }
  //
  // Unit conversion
  //
  inline double K_per_eV() { return 1.16045221e4; }
  inline double year() { return 3.1558150e7; }
  inline double asec_per_rad() { return 206264.806; }
  inline double sqasec_per_sr() { return 4.2545250225e10; }
  //
  // Astronomy
  //
  inline double Msun() { return 1.9891e33; }
  inline double Lsun() { return 3.839e33; }
  inline double Rsun() { return 6.955e10; }
  inline double AU() { return 1.49597870700e+13; }
  inline double parsec() { return 3.0856776e18; }

  ///
  /// tests if two doubles are equal to an accuracy of 1e-8.
  ///
  bool equalD(const double, const double);
  /// Function to replace pow(a, b) - exp(b*log(a)) is twice as fast
  double pow_fast(
      double,  ///< a
      double   ///< b
  );
};

extern class constants pconst;

#define ION_DUST -1
#define ION_H_N 0
#define ION_HE_N 1
#define ION_HE_P 2

#define EL_H 0
#define EL_HE 1

//-------------------------------------------------------
// Defines for different sized numbers.
//-------------------------------------------------------
#define SMALLVALUE 1.0e-12
#define MACHINEACCURACY 5.e-16
#define TINYVALUE 1.0e-100
#define VERY_TINY_VALUE 1.0e-200
#define VERY_LARGE_VALUE 1.0e100
#define HUGEVALUE 1.0e+200
#define ONE_PLUS_EPS (1.0 + SMALLVALUE)
#define ONE_MINUS_EPS (1.0 - SMALLVALUE)

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
#define Zcyl XX
#define Rcyl YY
#define Tcyl ZZ
#define ZNcyl XN
#define ZPcyl XP
#define RNcyl YN
#define RPcyl YP
#define TNcyl ZN
#define TPcyl ZP

//--------------------------------------------------------
///
/// Definitions for Spherical coordinate axes and directions w.r.t. cartesian
/// ones. [R,T,P] stands for [r,theta,phi]
///
//--------------------------------------------------------
#define Rsph XX
#define Tsph YY
#define Psph ZZ
#define RNsph XN
#define RPsph XP
#define TNsph YN
#define TPsph YP
#define PNsph ZN
#define PPsph ZP

//--------------------------------------------------------
///
/// Type of Solver:
///
//--------------------------------------------------------
#define FLUX_LF 0
#define FLUX_RSlinear 1  ///< Hydro/MHD: not very robust, not recommended
#define FLUX_RSexact                                                           \
  2  ///< Hydro only: exact solver. very slow, not recommended
#define FLUX_RShybrid 3  ///< Hydro combination of linear/exact solver. good.
#define FLUX_RSroe 4     ///< Hydro/MHD: Roe conserved vars solver
#define FLUX_RSroe_pv                                                          \
  5                 ///< Hydro: Roe solver in primitive vars (not recommended)
#define FLUX_FVS 6  ///< Hydro: flux vector splitting solver
#define FLUX_RS_HLLD                                                           \
  7  ///< MHD: The HLLD solver with switches to HLL for stability
#define FLUX_RS_HLL 8
///< Hydro/MHD the diffusive but robust HLL solver
#define FLUX_RCV_HLL 9
///< Hydro Roe-CV solver + HLL for large density contrast

//--------------------------------------------------------
/// \brief Primitive Variables Enum
/// The code requires that the velocities follow each other in order,
/// and the same for the Magnetic Field.  I'm not sure if it will work
/// without them being in that order, but I'm writing the code assuming
/// they are in order.
///
enum primitive {
  RO = 0,  ///< Density
  PG = 1,  ///< Gas pressure
  VX = 2,  ///< x-velocity
  VY = 3,  ///< y-velocity
  VZ = 4,  ///< z-velocity
  BX = 5,  ///< x-magnetic field
  BY = 6,  ///< y-magnetic field
  BZ = 7,  ///< z-magnetic field
  SI = 8   ///< The GLM scalar field, psi.
};
//--------------------------------------------------------

//--------------------------------------------------------
/// \brief Conserved Variables Enum
enum conserved {
  RHO = 0,  ///< Density
  ERG = 1,  ///< Total Energy
  MMX = 2,  ///< x-momentum
  MMY = 3,  ///< y-momentum
  MMZ = 4,  ///< z-momentum
  BBX = 5,  ///< x-magnetic field
  BBY = 6,  ///< y-magnetic field
  BBZ = 7,  ///< z-magnetic field
  PSI = 8   ///< The GLM scalar field, psi.
};
//--------------------------------------------------------

//--------------------------------------------------------
/// \brief Enum for each Spatial Dimension.
enum axes {
  XX = 0,  ///< X-Direction.
  YY = 1,  ///< Y-Direction.
  ZZ = 2   ///< Z-Direction.
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
enum direction {
  NO = -1,  ///< No Direction.
  XN = 0,   ///< Negative x-direction x-
  XP = 1,   ///< Positive x-direction x+
  YN = 2,   ///< Negative y-direction y-
  YP = 3,   ///< Positive y-direction y+
  ZN = 4,   ///< Negative z-direction z-
  ZP = 5    ///< Positive z-direction z+
};          // {x-,x+,y-,y+,z-,z+}

//
// Space and Time Integration Order of Accuracy
//
#define OA1 1  ///< First Order Space.
#define OA2 2  ///< Second Order Space.

//
// Viscosity identifier flags
//
#define AV_NONE 0  ///< No artificial viscosity.
#define AV_FKJ98_1D                                                            \
  1  ///< 1D linear viscosity of Falle, Komissarov, Joarder (1998).
#define AV_LAPIDUS 2      ///< Lapidus-type viscosity (BROKEN).
#define AV_HCORRECTION 3  ///< H-correction of Sanders et al. (1998).
#define AV_HCORR_FKJ98                                                         \
  4  ///< combination of the H-corr + FKJ98 (not needed really).
#define AV_VonNeuRicht                                                         \
  5  ///< von Neumann & Richtmeyer (1950) multi-D viscosity.

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

// *******************************************************************
// **************  Constants for Radiative Trasnfer   ****************
// *******************************************************************
//
// Method of calculating opacity and updating microphysics:
//
#define RT_UPDATE_IMPLICIT                                                     \
  1  ///< C2-ray update, where time-integration happens during ray-trace.
#define RT_UPDATE_EXPLICIT                                                     \
  2  ///< Only instantaneous column densities calculated during ray-trace.

//
// Types of source -- either a single source or one of 2*Ndim diffuse sources.
//
#define RT_SRC_DIFFUSE 2
#define RT_SRC_SINGLE 1

//
// Effect of source:
//
#define RT_EFFECT_UV_HEATING 1  ///< UV heating source.
#define RT_EFFECT_PION_MONO 2   ///< monochromatic photoion source.
#define RT_EFFECT_MFION 3       ///< multifrequency photoion source.
#define RT_EFFECT_PHOTODISS 4   ///< photodissociation of molecules.
#define RT_EFFECT_PION_EQM 5    ///< Assume photoion. equilibrium.
#define RT_EFFECT_HHE_MFQ                                                      \
  9  ///< Frank&Mellema H+He p-ion scheme (not fully implemented)

//
// Source of opacity
//
#define RT_OPACITY_TOTAL                                                       \
  1  ///< opacity is proportional to total mass density, so just integrate
     ///< rho*dr
#define RT_OPACITY_MINUS                                                       \
  2  ///< opacity provided by tracer i: integrate rho*(1-y_i)*dr
#define RT_OPACITY_TRACER                                                      \
  3  ///< opacity provided by tracer i: integrate rho*y_i*dr
#define RT_OPACITY_VSHELL                                                      \
  4  ///< should never be used.  RT uses it internally to set Vshell in each
     ///< cell.
#define RT_OPACITY_HALPHA                                                      \
  5  ///< Used only for analysis to calculate projected H-alpha emission.
#define RT_OPACITY_NII_FL                                                      \
  6  ///< Used only for analysis to calculate projected [NII] forbidden line
     ///< emission.
#define RT_OPACITY_RR 7   ///< Recombination rate (for ph-ion.eqm.)
#define RT_OPACITY_HHE 9  ///< H0,He0,He+,Dust.
#define RT_OPACITY_MP 10  ///< Get opacity from microphysics.
// *******************************************************************
// **************  Constants for Radiative Trasnfer   ****************
// *******************************************************************

#endif  // CONSTANTS_H
