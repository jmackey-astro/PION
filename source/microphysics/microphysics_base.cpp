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



// ##################################################################
// ##################################################################



microphysics_base::microphysics_base(
      const int nv,   ///< Total number of variables in state vector
      const int ntr,  ///< Number of tracer variables in state vector
      const std::string *tr,  ///< List of tracer variable names.
      struct which_physics *ephys, ///< which physics to calculate.
      struct rad_sources *rsrcs    ///< radiation sources.
      )
 : nv_prim(nv), ntracer(ntr)
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
  n_el = 0;
  for (int i=0;i<ntr;i++) {
    tr_map[tr[i]] = i+offset;
    // second map just for elements.
    if (tr[i].substr(0,2) == "X_") {
      el_map[tr[i]] = i+offset;
      el_index.push_back(i+offset);
      n_el++;
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



// ##################################################################
// ##################################################################



void microphysics_base::sCMA(
    pion_flt *corrector, ///< input corrector vector
    const pion_flt *p_in ///< input primitive vector (nv_prim)
    )
{
  //  Re-initialise corrector every step
  for (int i=0;i<nv_prim;i++) corrector[i] = 1.0;
  int print_flagg = 0;
  double corr = 0;
  
  // sum over all elements and add up mass fraction.
  for (int v=0;v<n_el;v++) {
    corr += p_in[el_index[v]];
  }
  corr = 1.0 / corr; // correction factor.
  // set correction factor for elements.
  for (int v=0;v<n_el;v++) {
    corrector[el_index[v]] = corr;
  }
  
  return;
}



// ##################################################################
// ##################################################################



int microphysics_base::Tr(const string s)
{
  if (tr_map.find(s)==tr_map.end())
    return -1;
  else
    return tr_map[s]; 
}




// ##################################################################
// ##################################################################



