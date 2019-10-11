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
      const int nv,   ///< Total number of variables in state vector
      const int ntr,  ///< Number of tracer variables in state vector
      const std::string *tr,  ///< List of tracer variable names.
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

  // set up a map of tracer names to indices in state vec.
  int offset = nv-ntr;
  for (int i=0;i<ntr;i++) {
    tr_map[tr[i]] = i+offset;
    // second map just for elements.
    if (tr[i].substr(0,2) == "X_") {
      el_map[tr[i]] = i+offset;
    }
  }
  //for (int v=0;v<nels;v++}
  for (auto elem : el_map) {
    cout << "Elements Map: " << elem.first << " " << elem.second << "\n";
  }
  for (auto elem : tr_map) {
    cout << "Tracers Map: " << elem.first << " " << elem.second << "\n";
  }
  
  return;
}


