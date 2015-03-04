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


#define SMALLVALUE       1.0e-12
#define MACHINEACCURACY  5.e-16
#define TINYVALUE        1.0e-100
#define VERY_TINY_VALUE  1.0e-200
#define VERY_LARGE_VALUE 1.0e100
#define HUGEVALUE        1.0e+200
#define ONE_PLUS_EPS     (1.0+SMALLVALUE)
#define ONE_MINUS_EPS    (1.0-SMALLVALUE)


#endif // SIM_CONSTANTS_H
