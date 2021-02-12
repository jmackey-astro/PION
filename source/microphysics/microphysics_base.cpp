/// \file microphysics_base.cpp
/// \author Jonathan Mackey
/// \date 2015.03.26

#include "microphysics_base.h"
#include "tools/reporting.h"
#include <iostream>
using namespace std;

///
/// Global pointer to microphysics class.
///
class microphysics_base* MP = 0;

// ##################################################################
// ##################################################################

microphysics_base::microphysics_base(
    const int nv,           ///< Total number of variables in state vector
    const int ntr,          ///< Number of tracer variables in state vector
    const std::string* tr,  ///< List of tracer variable names.
    struct which_physics* ephys,  ///< which physics to calculate.
    struct rad_sources* rsrcs     ///< radiation sources.
    ) :
    nv_prim(nv),
    ntracer(ntr)
{
#ifdef TESTING
    cout << "microphysics_base:: setting up EP and RS: ";
#endif

    EP = ephys;
    RS = rsrcs;

#ifdef TESTING
    cout << EP << "\t" << RS << "\n";
#endif

    // set up a map of tracer names to indices in state vec.
    int offset = nv - ntr;
    n_el       = 0;
    for (int i = 0; i < ntr; i++) {
        tr_map[tr[i]] = i + offset;
        tr_index.push_back(i + offset);
        // second map just for elements.
        if (tr[i].substr(0, 2) == "X_") {
            el_map[tr[i]] = i + offset;
            el_index.push_back(i + offset);
            n_el++;
        }
    }

#ifdef DEBUG_MP
    // for (auto elem : el_map) { // c++11 only (not working with intel compiler
    // 2018r4)
    for (std::map<std::string, int>::iterator it = el_map.begin();
         it != el_map.end(); ++it) {
        cout << "Elements Map: " << it->first << " " << it->second << "\n";
    }
    // for (auto elem : tr_map) { // c++11 only (not working with intel compiler
    // 2018r4)
    for (std::map<std::string, int>::iterator it = tr_map.begin();
         it != tr_map.end(); ++it) {
        cout << "Tracers Map: " << it->first << " " << it->second << "\n";
    }
#endif

    return;
}

// ##################################################################
// ##################################################################

//#define DEBUG_SCMA

void microphysics_base::sCMA(
    pion_flt* corrector,  ///< input corrector vector
    const pion_flt* p_in  ///< input primitive vector (nv_prim)
)
{
#ifdef DEBUG_SCMA
    int print_flagg = 0;
    double val      = 0.0;
#endif  // DEBUG_SCMA
    double corr = 0.0, ptemp[nv_prim];
    for (int i = 0; i < nv_prim; i++)
        corrector[i] = 1.0;
    // enforce all tracers to be within [0,1]
    for (int v = 0; v < ntracer; v++) {
        ptemp[tr_index[v]] = min(1.0, max(0.0, p_in[tr_index[v]]));
    }

    // sum over all elements and add up mass fraction, correcting
    // state vector to be within allowed limits for mass fraction.
    for (int v = 0; v < n_el; v++) {
        corr += ptemp[el_index[v]];
    }
    corr = 1.0 / corr;  // correction factor for elements

    // set multiplier to ensure that corrected tracer value
    // is within [0,1] for all tracers, if current value is out of
    // range
    for (int v = 0; v < ntracer; v++) {
#ifdef DEBUG_SCMA
        val = ptemp[tr_index[v]] / p_in[tr_index[v]];
#endif  // DEBUG_SCMA
        // if p_in[v] is out of range, then set corrector to fix this
        corrector[tr_index[v]] = (p_in[tr_index[v]] < 0.0) ? 0.0 : 1.0;
        corrector[tr_index[v]] =
            (p_in[tr_index[v]] > 1.0) ? 1.0 / p_in[tr_index[v]] : 1.0;
    }

    // renormalise element values to sum to 1.
    for (int v = 0; v < n_el; v++) {
        corrector[el_index[v]] *= corr;
    }

#ifdef DEBUG_SCMA
    if (fabs(val - 1.0) > 1.0e-1) {
        cout << "------------  sCMA  ---------------\n";
        rep.printVec("p_in     : ", p_in, nv_prim);
        rep.printVec("ptemp    : ", ptemp, nv_prim);
        rep.printVec("Corrector: ", corrector, nv_prim);
    }
#endif  // DEBUG_SCMA

    return;
}

// ##################################################################
// ##################################################################

int microphysics_base::Tr(const string s)
{
    if (tr_map.find(s) == tr_map.end())
        return -1;
    else
        return tr_map[s];
}

// ##################################################################
// ##################################################################
