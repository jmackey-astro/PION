/// \file xray-emission.cpp
/// \author Jonathan Mackey
/// \date 2016.07.04
///
/// Purpose:
/// - setup table for calculating X-ray emission.
///

#include <cmath>

#include <sstream>
#include <vector>

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/interpolate.h"
#include "tools/mem_manage.h"



#include "xray_emission.h"
using namespace std;


// ##################################################################
// ##################################################################


Xray_emission::Xray_emission()
{
  XNel    = 0;
  LT      = 0;
  L1      = 0;
  L2      = 0;
  L3      = 0;
  L4      = 0;
  int err = setup_xray_tables_priv(&XNel, &LT, &L1, &L2, &L3, &L4);
  if (err) {
    cerr << "Failed to setup xray tables correctly: " << err << "\n";
    exit(err);
  }
  return;
}


// ##################################################################
// ##################################################################



///
/// Function to set up the X-ray emissivity tables, from XSPEC
///
int Xray_emission::setup_xray_tables_priv(
    size_t *N,
    double **logT,
    double **logLum1,
    double **logLum2,
    double **logLum3,
    double **logLum4)
{
  //
  // Read data from xray-table.txt
  //
  //
  // read input parameters:
  //
  ifstream ff;
  ff.open("xray-table.txt");
  if (!ff.is_open()) {
    cerr << "Error opening file.\n";
    return 1;
  }

  //
  // Read the first line, which is the column headers:
  //
  string line;
  getline(ff, line);

  //
  // Now read in the data line-by-line.  Cols. are:
  // log(T/K) T(K) E(keV) L(E>0.1keV) L(E>0.5keV) L(E>1keV) L(E>5keV)
  //
  // We need to store logT, L1, L2, L3, L4
  //
  vector<double> lLT, lL1, lL2, lL3, lL4;
  double lt = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
  double temp = 0;
  size_t Nel  = 0;
  while (!ff.eof()) {
    getline(ff, line);
    // If it's not a line containing a parameter, continue and read another
    // line.
    if (line.empty() == true || (line.find("#") != string::npos)) {
      continue;
    }
    else {
      // We have found a line of data, so read it into variables.
      istringstream iss(line);
      iss >> lt >> temp >> temp >> l1 >> l2 >> l3 >> l4;
#ifndef NDEBUG
      cout << lt << "  " << l1 << "  " << l2 << "  " << l3 << "  " << l4
           << "\n";
#endif
      //
      // Add all of them to the vectors:
      //
      lLT.push_back(lt);
      lL1.push_back(l1);
      lL2.push_back(l2);
      lL3.push_back(l3);
      lL4.push_back(l4);
      Nel++;
    }
  }
  cout << "Read " << Nel << " lines of data from input xray tables.\n";
  ff.close();


  //
  // allocate arrays
  //
  (*logT)  = mem.myalloc((*logT), Nel);
  *logLum1 = mem.myalloc(*logLum1, Nel);
  *logLum2 = mem.myalloc(*logLum2, Nel);
  *logLum3 = mem.myalloc(*logLum3, Nel);
  *logLum4 = mem.myalloc(*logLum4, Nel);

  for (size_t v = 0; v < Nel; v++) {
    (*logT)[v]    = lLT[v];
    (*logLum1)[v] = log10(lL1[v]);
    (*logLum2)[v] = log10(lL2[v]);
    (*logLum3)[v] = log10(lL3[v]);
    (*logLum4)[v] = log10(lL4[v]);
  }

  *N = Nel;

  return 0;
}


// ##################################################################
// ##################################################################



void Xray_emission::free_xray_tables_priv(
    double **logT,
    double **logL1,
    double **logL2,
    double **logL3,
    double **logL4)
{
  *logT  = mem.myfree(*logT);
  *logL1 = mem.myfree(*logL1);
  *logL2 = mem.myfree(*logL2);
  *logL3 = mem.myfree(*logL3);
  *logL4 = mem.myfree(*logL4);
  return;
}

// ##################################################################
// ##################################################################



/// Function to calculate the X-ray emissivity given an input gas
/// temperature.  Returns the emissivity in units of ergs.cm^3/sec,
/// and should be multiplied by n_e*n_H.
///
void Xray_emission::get_xray_emissivity(
    const double T,  ///< Temperature (K)
    double *res      ///< Results.
)
{
  double ln10  = 2.302585093;
  double lrate = 0.0;
  double lt    = log10(T);

  //
  // extrapolate linearly
  //
  if (lt < LT[0]) {
    for (size_t v = 0; v < 4; v++)
      res[v] = 0.0;
  }
  else if (lt > LT[XNel - 1]) {
    cout << "extrapolate, T=" << T << " : ";
    res[0] = L1[XNel - 1]
             + (L1[XNel - 1] - L1[XNel - 2]) * (lt - LT[XNel - 1])
                   / (LT[XNel - 1] - LT[XNel - 2]);
    res[1] = L2[XNel - 1]
             + (L2[XNel - 1] - L2[XNel - 2]) * (lt - LT[XNel - 1])
                   / (LT[XNel - 1] - LT[XNel - 2]);
    res[2] = L3[XNel - 1]
             + (L3[XNel - 1] - L3[XNel - 2]) * (lt - LT[XNel - 1])
                   / (LT[XNel - 1] - LT[XNel - 2]);
    res[3] = L4[XNel - 1]
             + (L4[XNel - 1] - L4[XNel - 2]) * (lt - LT[XNel - 1])
                   / (LT[XNel - 1] - LT[XNel - 2]);
    for (size_t v = 0; v < 4; v++)
      res[v] = exp(ln10 * res[v]);
    rep.printVec("Res", res, 4);
  }
  else {
    interpolate.root_find_linear(LT, L1, XNel, lt, &lrate);
    res[0] = exp(ln10 * lrate);
    interpolate.root_find_linear(LT, L2, XNel, lt, &lrate);
    res[1] = exp(ln10 * lrate);
    interpolate.root_find_linear(LT, L3, XNel, lt, &lrate);
    res[2] = exp(ln10 * lrate);
    interpolate.root_find_linear(LT, L4, XNel, lt, &lrate);
    res[3] = exp(ln10 * lrate);
  }

  return;
}


// ##################################################################
// ##################################################################
