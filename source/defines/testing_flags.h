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

#define REPORT_RANK0  // only reports stdout from rank 0.
#define C2F_FULLSTEP  // Coarse-to-fine boundaries updated every 2 steps

#define NEWGRIDDATA
//#define SKIP_C2F_BC
//#define SKIP_F2C_BC
//#define SKIP_BC89_FLUX
//#define TEST_BC89FLUX        // debugging info for BC89 flux correction
//#define TEST_COMMS           // debugging info for MPI
//#define TEST_CONSERVATION    // test mass/mom/energy conservation
//#define TEST_INF             // test for infinity/nan
//#define TEST_SYMMETRY // dangerous for real simulations (assumes units for
// rho,v) #define RT_TESTING ///< Enable this for debugging info on raytracing.

//
// Various defines for testing bits of the code (lots of sdtout).
//
// ///< Comment this out if not testing the code.

// displays debugging information for time integration.
//#define TEST_INT

// This counts the energy loss/gain in microphysics (SLOW!!!)
//
//#define COUNT_ENERGETICS

//
// prints a line naming a function every time it is called.
//
//#define FUNCTION_ID

//#define RSTESTING ///< If testing the Riemann Solvers.

//#define AVFALLE_TESTING ///< This is for testing the FKJ98 Art.Viscosity
// Function.

///
/// Sets fixed x-sections for ionisation and recombination,
/// and disables timestep limiting in first 3 steps.
///
//#define RT_TEST_PROBS
#define RT_TESTING_OUTPUTCOL  ///< output column density info to files.

/// for the field loop advection test problem only.
//#define CHECK_MAGP

/// for a 1D spherical blast wave test problem only.
//#define BLAST_WAVE_CHECK

//#define GRIDV2

#endif  // NDEBUG_FLAGS_H
