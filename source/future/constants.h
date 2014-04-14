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


#ifndef CONSTANTS_H
#define CONSTANTS_H

class constants {
  public:
  //
  // Mathematics
  //
  inline double ln10() {return 2.3025850929940459;}
  //
  // Physics
  //
  inline double pi()     {return 3.14159265358979324;}
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
};


#define ION_DUST -1
#define ION_H_N  0
#define ION_HE_N 1
#define ION_HE_P 2

#define EL_H  0
#define EL_HE 1


#endif // CONSTANTS_H


