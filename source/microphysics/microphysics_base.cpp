/// \file microphysics_base.cpp
/// \author Jonathan Mackey
/// \date 2015.03.26


#include "microphysics_base.h"
#include <iostream>
using namespace std;

///
/// Global pointer to microphysics class.
///
class microphysics_base *MP=0;


microphysics_base::microphysics_base(
      struct which_physics *ephys, ///< which physics to calculate.
      struct rad_sources *rsrcs    ///< radiation sources.
      )
{
#ifdef TESTING
  cout <<"microphysics_base:: setting up EP and RS: ";
#endif

  EP = ephys;
  RS = rsrcs;

#ifdef TESTING
  cout <<EP<<"\t"<<RS<<"\n";
#endif

  return;
}


