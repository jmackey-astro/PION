/// \file icgen_base.cpp
/// \author Jonathan Mackey
/// \date 2018-05-04
///
/// Description:\n
/// This implements a set of routines that are common to serial,
/// parallel, and NG grid code for generating initial conditions.
///
/// Modifications:\n
/// 2018.05.04 JM: moved code from icgen.cpp/.h
///

#ifndef ICGEN_BASE_H
#define ICGEN_BASE_H

#include <string>

class ICsetup_base {
public:
  virtual ~ICsetup_base() {}

  virtual int setup_data(
      class ReadParams *,    ///< pointer to parameter list.
      class GridBaseClass *  ///< pointer to grid
      ) = 0;

  void set_SimPM(class SimParams *sp)
  {
    SimPM = sp;
    return;
  }

  /// allows you to set a pointer to a microphysics class
  void set_mp_pointer(
      class microphysics_base *ptr  ///< pointer to microphysics class
  )
  {
    mp = ptr;
    return;
  }

#ifdef PARALLEL
  void set_sub_domain_pointer(class Sub_domain *m)
  {
    sub_domain = m;
    return;
  }
#endif  // PARALLEL

  int equilibrate_MP(
      class GridBaseClass *,
      class microphysics_base *,
      class ReadParams *,
      class SimParams &);

protected:
  int AddNoise2Data(
      class GridBaseClass *,
      class SimParams &,
      int,    ///< type of noise (1=pressure,2=adiabatic,3=adiabatic wave)
      double  ///< Noise level (in pressure) in fractional level.
  );

  int SmoothData(int  ///< Number of cell diameters to smooth over.
  );

  class GridBaseClass *gg;  ///< pointer to grid.
  class ReadParams *rp;     ///< pointer to readparams.
  class SimParams *SimPM;   ///< pointer to simulation parameters
#ifdef PARALLEL
  class Sub_domain *sub_domain;
#endif                          // PARALLEL
  class microphysics_base *mp;  ///< pointer to microphysics class
};

///
/// Assigns a type of initial conditions based on input string, and
/// intialised the correct clas and returns a pointed to it.
///
void setup_ics_type(
    std::string,           ///< string giving type of ICs
    class ICsetup_base **  ///< pointer to address of IC class
);

#endif  // ICGEN_BASE_H
