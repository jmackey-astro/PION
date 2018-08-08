/// \file update_boundaries.cpp
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class definitions for routines to update grid boundaries with
/// different boundary conditions.
///
/// Modifications:\n
/// - 2018.05.10 JM: moved code from uniform_grid.cpp
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "sim_control/update_boundaries.h"


// ##################################################################
// ##################################################################

update_boundaries::update_boundaries()
{
  spatial_solver=0;
  return;
}

update_boundaries::~update_boundaries()
{
  if (spatial_solver) {delete spatial_solver; spatial_solver=0;}
  return;
}


// ##################################################################
// ##################################################################


