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
/// - getting it written: mods up until 2013.08.19


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
  inline double c()      {return 2.9979e10;}
  inline double kB()     {return 1.381e-16;}
  inline double h()      {return 6.626e-27;}
  inline double StefanBoltzmannConst() {return 5.670e-5;}
  inline double m_p()    {return 1.673e-24;}
  //
  // Atomic data
  //
  inline double Eth_H()  {return 13.59844*eV();}
  inline double NuTh_H() {return 3.288e15;}
  //
  // Unit conversion
  //
  inline double K_per_eV() {return 1.16045e4;}
  inline double year()     {return 3.15576e7;}
  inline double eV()       {return 1.602e-12;}
  //
  // Astronomy
  //
  inline double Msun()   {return 1.989e33;}
  inline double parsec() {return 3.086e18;}
  inline double Lsun()   {return 1.989e33;}
  inline double Rsun()   {return 6.960e10;}
};


#endif // CONSTANTS_H


