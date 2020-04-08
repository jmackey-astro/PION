/// \file emission_absorption.h
/// \author Jonathan Mackey
/// \date 2016.07.04
///
/// Purpose:
/// - various functions for calculating emission/absorption
///   coefficients, and calculating projected quantities.
///

#ifndef EMISSION_ABSORPTION_H
#define EMISSION_ABSORPTION_H

//
// data read from file is put in a 2D array with these values for
// the first index.
//
#define DATA_R   0 ///< radius
#define DATA_D   1 ///< density
#define DATA_P   2 ///< pressure
#define DATA_V   3 ///< velocity
#define DATA_TR0 4 ///< tracer 0
#define DATA_TR1 5 ///< tracer 1
#define DATA_T   6 ///< Temperature

#define PROJ_D   0  ///< projected density
#define PROJ_NtD 1  ///< projected neutral density
#define PROJ_InD 2  ///< projected ionised density
#define PROJ_EM  3  ///< projected emission measure.
#define PROJ_X01 4  ///< X-ray emission >0.1 keV
#define PROJ_X05 5  ///< X-ray emission >0.5 keV
#define PROJ_X10 6  ///< X-ray emission >1.0 keV
//#define PROJ_X20 7  ///< X-ray emission >2.0 keV
#define PROJ_X50 7  ///< X-ray emission >5.0 keV
#define PROJ_HA  8  ///< projected H-alpha emission.
#define PROJ_IML 9 ///< projected ionised metal-line.




//
// From the simulation data arrays, calculate various variables.
//
int get_emission_absorption_data(
      double const* const*, ///< raw data to get variable from
      const int,    ///< number of images to write
      const size_t, ///< Number of radial data elements
      class Xray_emission &,  ///< pointer to X-ray emission class.
      double **,    ///< array for emission[img][rad] data.
      double **     ///< array for absorption[img][rad] data.
      );

//
// Project scalar quantities onto plane of the sky.
//
double calc_projection(
      const double *, ///< radius array
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius
      );

//
// Project quantities with emission and absorption onto sky.
//
double calc_projectionRT(
      const double *, ///< radius array.
      const double *, ///< array of emission vals at each radius
      const double *, ///< array of absorption vals at each radius
      const size_t ,  ///< Size of arrays.
      const double ,  ///< impact parameter of ray.
      const double    ///< spacing of points in radius.
      );
#endif // EMISSION_ABSORPTION_H

