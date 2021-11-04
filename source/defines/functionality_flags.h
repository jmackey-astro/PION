///
/// \file functionality_flags.h
/// \author Jonathan Mackey
/// \date 2010.12.27
///
/// This file sets defines which can enable or disable large sections
/// of code for compilation.  The idea is to allow very efficient
/// code to be compiled for production runs.  It can also be used to
/// switch on/off code which is currently under development and may
/// not even compile without errors.
///
/// Modifications:
///
/// - 2010.12.27 JM: Created.
/// - 2012.01.26 JM: Added flags for disabling complilation of microphysics
///    modules.
/// - 2012.04.19 JM: Added PUREHYDROGEN class (this should get its own
///    microphysics derived class soon).
///

#ifndef FUNCTIONALITY_FLAGS_H
#define FUNCTIONALITY_FLAGS_H

// -------------------------------------------------------------------
// Functionality:
// -------------------------------------------------------------------

//
// Whether to use single or double precision
//
#define PION_DATATYPE_DOUBLE
#define pion_flt double
//#define PION_DATATYPE_FLOAT
//#define pion_flt float

//
// This flag prevents compilation if any parts of the code are flagged as HACKs
//
//#define HACK_WARNING

// this flag enables input B-fields in CGS units (Gauss)
#define NEW_B_NORM

//
// These flags disable microphysics modules
//
#define EXCLUDE_HD_MODULE  // Harpreet's module.
#define EXCLUDE_MPV9       // mpv9_HHe module.

//
// Comment out to disable PDR cooling in mp_explicit/implicit_H.cc
// (only do this if there is no molecular gas in the simulation!)
//
#define USE_PDR_COOLING
//#define PUREHYDROGEN ///< in MPv3, disables all He+metal terms.

//
// These flags set the timestep-limiting option for MPv3 and
// MPv4.  Only modify them if you know what you are doing.
// See Mackey (2012,A&A,539,A147) for details.
//
#define MPV3_DTLIMIT 2
//#define MPV3_DTLIMIT 3
#define MPV4_DTLIMIT 5

///
/// With this enabled, we reset negative pressures such that the
/// temperature is 10K, rather than a fixed pressure floor value.
/// N.B. This must be unset for sims with no MP-> module, because
/// it calls the MP->Temperature() function!
///
#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

///
/// This prevents the timestep increasing by more than 30% per step,
/// which is reputed to help with stability (can't remember refs).
///
#define TIMESTEP_LIMITING

// -------------------------------------------------------------------
// Flags for parts of code under development.
// -------------------------------------------------------------------

///
/// Flag for whether to do thermal conduction or not...
/// This is not yet working well -- the solution is very sketchy.
#define THERMAL_CONDUCTION

//
// Stuff for the internal energy solver (which only works in cartesian
// coordinates, and is *VERY* diffusive, so not recommended).
//
//#define INCLUDE_EINT_ADI_HYDRO
//#define EINT_ETOT_PARALLEL
//#define VNR_VISCOSITY

//#define ISOTHERMAL_SOLVERS_ENABLED

// -------------------------------------------------------------------
// Flags which should be removed/changed because the code is working
// and tested.
// -------------------------------------------------------------------

#endif  // FUNCTIONALITY_FLAGS_H
