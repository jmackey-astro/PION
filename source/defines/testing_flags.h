///
/// \file testing_flags.h
/// \author Jonathan Mackey
/// \date 2010.12.27
///
/// This file sets defines which can enable or disable large sections
/// of code for compilation.  Flags defined here are just for testing
/// the code, and most (all?) should be disabled for production runs.
///
/// Modifications:
///
/// - 2011.01.06 JM: Created.
///

#ifndef TESTING_FLAGS_H
#define TESTING_FLAGS_H

//
// Various defines for testing bits of the code.
//
//#define TESTING ///< Comment this out if not testing the code.

//
// This counts the energy loss/gain in microphysics (SLOW!!!)
//
//#define COUNT_ENERGETICS

//
// prints a line naming a function every time it is called.
//
//#define FUNCTION_ID

#define RT_TESTING ///< Enable this for debugging info on raytracing.

//#define RSTESTING ///< If testing the Riemann Solvers.

//#define AVFALLE_TESTING ///< This is for testing the FKJ98 Art.Viscosity Function.



///
/// Sets fixed x-sections for ionisation and recombination,
/// and disables timestep limiting in first 3 steps.
///
//#define RT_TEST_PROBS
#define RT_TESTING_OUTPUTCOL ///< output column density info to files.

/// for the field loop advection test problem only.
//#define CHECK_MAGP

/// for a 1D spherical blast wave test problem only.
//#define BLAST_WAVE_CHECK


//#define GRIDV2

#endif // TESTING_FLAGS_H
