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


//
// These flags disable JM's microphysics modules and (in the last case)
// Harpreet's microphysics module.
//
//#define EXCLUDE_MPV1       // both microphysics_v1 and H-only MP.
#define EXCLUDE_MPV2       // mp_v2_aifa
//#define EXCLUDE_MPV3       // mp_explicit_H 
#define EXCLUDE_MPV4       // mp_implicit_H
#define EXCLUDE_HD_MODULE  // Harpreet's module.
//#define EXCLUDE_MPV9       // mpv9_HHe module.

//
// Comment out to disable PDR cooling in mp_explicit/implicit_H.cc
// (only do this if there is no molecular gas in the simulation!)
//
#define USE_PDR_COOLING
//#define PUREHYDROGEN ///< in MPv3, disables all He+metal terms.

///
/// This enables the mp_v2_aifa code, which is an old version of mp_explicit_H
///
//#define MP_V2_AIFA


//
// These flags set the timestep-limiting option for mp_explicit_H and
// mp_implicit_H.  Only modify them if you know what you are doing.
// See Mackey (2011,A&A,539,A147) for details.
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
/// For the new time-integration scheme
///
#define NEW_TIME_UPDATE

///
/// Flag for whether to do thermal conduction or not...
/// This is not yet working well -- the solution is very sketchy.
//#define THERMAL_CONDUCTION

//
// Stuff for the internal energy solver (which only works in cartesian
// coordinates, and is *VERY* diffusive, so not recommended).
//
//#define INCLUDE_EINT_ADI_HYDRO
//#define EINT_ETOT_PARALLEL
//#define VNR_VISCOSITY

//#define ISOTHERMAL_SOLVERS_ENABLED

//#define LAPIDUS_VISCOSITY_ENABLED
//#define TEST_LAPIDUS

// -------------------------------------------------------------------
// Flags which should be removed/changed because the code is working
// and tested.
// -------------------------------------------------------------------


///
/// Use try{}catch(){} for bad_alloc in dynamic memory allocation.
///
#define CHECK_NEW_EXCEP_ON
#ifdef CHECK_NEW_EXCEP_ON
#include <new>
#endif

/// Use the memory management class for dynamic memory allocations and frees.
#define USE_MM
#ifndef USE_MM
#error "MUST USE MEMORY MANAGEMENT NOW -- ALL NEW CODE ASSUMES IT IS PRESENT!!!!"
#endif

#define PLLEL_RT ///< For point sources with multiple processes (ALWAYS SET NOW)

///
/// Use an updated interface to the microphysics. (MUST USE THIS NOW!)
/// (removed 2011.02.25 JM)
//#define NEW_RT_MP_INTERFACE

//
// Must define either CELL_CENTRED_SRC or NON_CELL_CENTRED_SRC
//
// This enables code for ray-tracing point sources where the src is
//  at a cell centre.
//
//#define CELL_CENTRED_SRC 
//
// This enables code for non-cell-centred sources 
//
#define NON_CELL_CENTRED_SRC

#endif //FUNCTIONALITY_FLAGS_H
