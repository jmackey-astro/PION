///
/// \file sim_constants.h
/// \author Jonathan Mackey
///
/// Constants that are used in many places in the code.
///
/// Modifications:
/// - 2015.01.10 JM: Wrote file, moved constants from global.h

#ifndef SIM_CONSTANTS_H
#define SIM_CONSTANTS_H

///
/// Max. number of spatial dimensions for a simulation.
///
#define MAX_DIM 3
///
/// Maximum length of state vectors.
///
#define MAX_NVAR 70

/// Minimum number of cells radius for the stellar-wind boundary
#define MIN_WIND_RAD 7

#endif  // SIM_CONSTANTS_H
