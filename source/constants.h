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


#define SMALLVALUE       1.0e-12
#define MACHINEACCURACY  5.e-16
#define TINYVALUE        1.0e-100
#define VERY_TINY_VALUE  1.0e-200
#define VERY_LARGE_VALUE 1.0e100
#define HUGEVALUE        1.0e+200
#define ONE_PLUS_EPS     (1.0+SMALLVALUE)
#define ONE_MINUS_EPS    (1.0-SMALLVALUE)


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


#endif // CONSTANTS_H


