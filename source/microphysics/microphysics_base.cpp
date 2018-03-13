/// \file microphysics_base.cpp
/// \author Jonathan Mackey
/// \date 2015.03.26


#include "microphysics_base.h"

///
/// Global pointer to microphysics class.
///
class MicroPhysicsBase *MP=0;


MicroPhysicsBase::MicroPhysicsBase(
      struct which_physics *ephys, ///< which physics to calculate.
      struct rad_sources *rsrcs    ///< radiation sources.
      )
{
  EP = ephys;
  RS = rsrcs;
  return;
}


