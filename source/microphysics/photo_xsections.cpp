///
/// \file photo_xsections.cpp
/// \author Maggie Celeste Goulden & Jonathan Mackey
/// \date 2019.08.14
///
/// This class calculates the photoionisation cross-sections of
/// ions using VernerFerlandEA1996_ApJ_465_187.
/// References:
/// - VF96: Verner and Ferland (1996,ApJ,465,187) Table 1 & Eqn 1
/// - DYL99: Dalgarno, Yan, Liu (1999,ApJS,125,237) Table 7
/// - MS05:  Meijerink & Spaans (2005,A&A,436,397)  App. B
/// - PCP12: Panoglou,Cabrit,+  (2012,A&A,538,A2)
/// - Yan (1997) PhD thesis:
///     http://adsabs.harvard.edu/abs/1997PhDT.........4Y
/// - Meijerink (2006) PhD thesis:

#include "photo_xsections.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <map>
#include <stdlib.h>

#include "tools/mem_manage.h"

using namespace std;

#ifndef INTEL
#include <cmath>     // Header file from gcc
#else
#include <mathimf.h> // Header file from Intel Compiler
#endif


// ##################################################################
// ##################################################################
///
/// constructor: inputs are Emin, Emax, nbins.
/// Emin must be >= 0.03 keV  
/// Emax must be <= 20   keV  ???
/// 

//TODO: DELETE MAIN() ONCE YOU'VE ASKED JONATHAN WHERE THIS FUNCTION IS ACTUALLY CALLED FROM...
/*int main()
{double emax[15] = {13.59844, 14.5, 24.4, 24.58741, 29.6, 47.5, 47.9, 54.41778, 64.5, 77.5, 97.9, 392.1, 490.0, 552.1, 667.0};//bin edges correspond to ionisation energies
  double emin[15] = {11.3, 13.59844, 14.5, 24.4, 24.58741, 29.6, 47.5, 47.9, 54.41778, 64.5, 77.5, 97.9, 392.1, 490.0, 552.1};
  //convert to kev:
  for (int i=0; i<15; i++){
      emax[i] *= 1e-3;
      emin[i] *= 1e-3;
  }
  photo_xsections *qs = new photo_xsections();//&emin[0], &emax[0], 15);
}*/
    
photo_xsections::photo_xsections(
    /*const double *emin, ///< Emin  (eV)
    const double *emax, ///< Emax  (eV)
    const int    nbins   ///< nbins*/
    )
{
  // initialise ions list.
  n_ions = 78;
  ions = mem.myalloc(ions,78);
  // list of ions in the order they appear in VF96.  Note this is
  // photoionisation rate INTOOO these ions.
  string temp_ion[78] = 
    {"H1+", "He1+","He2+", "Li1+","Li2+","Li3+", "Be1+","Be2+","Be3+","Be4+",
     "B1+","B2+","B3+","B4+","B5+", "C1+","C2+","C3+","C4+","C5+","C6+",
     "N1+","N2+","N3+","N4+","N5+","N6+","N7+", 
     "O1+","O2+","O3+","O4+","O5+","O6+","O7+","O8+",
     "F1+","F2+","F3+","F4+","F5+","F6+","F7+","F8+","F9+",
     "Ne1+","Ne2+","Ne3+","Ne4+","Ne5+","Ne6+","Ne7+","Ne8+","Ne9+","Ne10+",
     "Na1+","Na2+","Na3+","Na4+","Na5+","Na6+","Na7+","Na8+","Na9+","Na10+","Na11+",
     "Mg1+","Mg2+","Mg3+","Mg4+","Mg5+","Mg6+","Mg7+","Mg8+","Mg9+","Mg10+","Mg11+","Mg12s+"
    };
  for (int v=0;v<n_ions;v++) ions[v] = temp_ion[v];

  double emax[12] = {11.2, 13.6, 16.3, 21.56, 24.6, 29.6, 35.12, 40.96, 47.89, 54.4, 64.5, 77.0};//bin edges correspond to ionisation energies
  double emin[12] = {7.64, 11.2, 13.6, 16.3, 21.56, 24.6, 29.6, 35.12, 40.96, 47.89, 54.4, 64.5};
  int nbins = 12;
  // load in parameter file from VF1996
  // TODO: ADD IN OTHER IONS BEYOND H0, He0, He1+
  FILE * pf = fopen("./source/microphysics/VernerFerlandEA1996_ApJ_465_187_photoionisation_crosssection.dat", "r");
  for (int i=0; i<n_ions; i++){
    struct params_struct p;
    int ignore = fscanf(pf, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
           &p.Z, &p.N, &p.E_th, &p.E_max, &p.E_0, &p.sigma_0, &p.y_a, &p.P, &p.y_w, &p.y_0, &p.y_1);
    if (ignore != 11) cerr << "Not all parameters read in! photo_xsections \n";
    parameters[ ions[i] ] = p;
  }
  //
  // set up photo_xsections arrays for the energy bins.
  //
  
  photo_xsections::Nbin = nbins;
  photo_xsections::Emin = emin[0];
  photo_xsections::Emax = emax[nbins-1];

  // energies
  photo_xsections::EBmin = new double [nbins];
  photo_xsections::EBmax = new double [nbins];
  photo_xsections::E0    = new double [nbins];
  photo_xsections::dE    = new double [nbins];

  
  //
  // Energy of each bin is the midpoint of the bin (lin-space).
  //
  for (unsigned int v=0;v<Nbin;v++) {
    EBmin[v] = emin[v];
    EBmax[v] = emax[v];
    E0[v] = 0.5*(EBmin[v]+EBmax[v]);
    dE[v] = EBmax[v]-EBmin[v];
  }
  
  //
  // Setting the cross-section is not trivial...
  // Iterate over N ions to get N cross-section arrays
  // units of cm^2
  for (int N=0; N<n_ions; N++){
    // cross section at each energy
    photo_xsections::xs = new double [nbins];
    //Iterate over each energy bin
    for (unsigned int v=0;v<Nbin;v++) {
        xs[v] = set_bin_mean_xsection(EBmin[v],EBmax[v], ions[N]);
        cout << "mean_xsection=" << xs[v] << " for " << ions[N] << " at E<" << EBmax[v] << "eV\n";
    //cout <<v <<"  "<< EBmin[v] <<"  "<< EBmax[v] <<"  ";
    //cout << E0[v] <<"  "<< dE[v] <<"  "<< xs[v] <<"\n";
    }
    parameters[ions[N]].xsections = xs;
    /*cout << xs[8] << "\n";
    cout << parameters[ions[N]].xsections[8] << "\n";*/
  }
  cout << parameters["H1+"].xsections[5] << "\n";
  //cout <<"E min/max="<<Emin<<", "<<Emax<<",  Nbins="<<Nbin<<"\n";
  return;
}


// ##################################################################
// ##################################################################



///
/// Destructor: cleans up memory.
///
photo_xsections::~photo_xsections()
{
  delete [] E0; E0=0;
  delete [] xs; xs=0;
  delete [] EBmin; EBmin=0;
  delete [] EBmax; EBmax=0;
  delete [] dE; dE=0;

#ifdef NEON
  delete [] Ne_xs; Ne_xs=0;
#endif // NEON

  return;
}


// ##################################################################
// ##################################################################


///
/// Calculate cross-section of the gas, per ion nucleus, as a 
/// function of energy
///
double photo_xsections::ion_xsection(
  const double E, ///< energy of photons, in eV
  string ion ///< name of ion
  )
{
  params_struct p = parameters[ion];
  if ( E < p.E_th) return 0;
  
  double x = E/p.E_0 - p.y_0;
  double y = sqrt(x*x + p.y_1*p.y_1);
  double F_y = ((x-1.)*(x-1.) + p.y_w*p.y_w)* pow(y, 0.5*p.P - 5.5) * pow(1 + sqrt(y/p.y_a), -p.P);
  double sigma_E = p.sigma_0*F_y;
  //cout << "Sigma=" << sigma_E << " for " << ion << " at E=" << E << "eV\n";
  return sigma_E*1e-18;
}


// ##################################################################
// ##################################################################


///
/// Calculates the mean cross-section of an energy bin between
/// two limits, using a mean attenuation.
///
double photo_xsections::set_bin_mean_xsection(
  const double em, ///< Minimum energy in bin (eV)
  const double ep,  ///< Maximum energy in bin (eV)
  string ion_name ///< name of ion
  )
{
  double ea = 0.5*(em+ep);
  unsigned int N=100; // number of intervals for simpsons rule
  double h = (ep-em)/N; 
  double NH = 1.0/ion_xsection(ea, ion_name);

  double ans = 0.0;
  //
  // Simpson's rule, first add the lower and upper limits:
  //
  ans += exp(-NH*ion_xsection(em, ion_name));
  ans += exp(-NH*ion_xsection(ep, ion_name));
  //
  // now intermediate points
  //
  unsigned int weight=4;
  for (unsigned int i=1;i<N;i++) {
    ea = em+h*i;
    ans += weight*exp(-NH*ion_xsection(ea, ion_name));
    weight = 6-weight;
  }
  ans *= h/3.0;  // this is the integral of exp(-tau) over em->ep.
  //
  // Now 
  // - divide out the energy range to get a mean exp(-tau),
  // - take log to get -<tau>
  // - divide by NH to get <sigma>
  //
  ans = -1.0*log(ans/(ep-em))/NH;
  if (!isfinite(ans)) return 0;
  return ans;
}

