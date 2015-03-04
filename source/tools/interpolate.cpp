///
/// \file interpolate.cpp
///  
/// \brief General purpose class for interpolating tables.
/// \author Jonathan Mackey
/// 
/// Modifications:
/// - 2015.03.04 JM: moved from global.h GeneralStuff class.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/interpolate.h"
#include "sim_constants.h"
#include "tools/reporting.h"

#ifndef INTEL
#include <cmath>     // Header file from gcc
#else
#include <mathimf.h> // Header file from Intel Compiler
#endif

///
/// global instance of class, defined in interpolate.cpp.
///
class interpolate_arrays interpolate;

// ##################################################################
// ##################################################################

interpolate_arrays::interpolate_arrays()
{
}

// ##################################################################
// ##################################################################

interpolate_arrays::~interpolate_arrays() 
{
  //cout <<"timers.size: "<<timers.size()<<"\n";
}



// ##################################################################
// ##################################################################



void interpolate_arrays::spline(const double *x,
			  const double *y,
			  const int n,
			  double yp1,
			  double ypn,
			  double *y2
			  )
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double [n];
  
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  //cout <<"x[0]="<<x[0]<<" y[0]="<<y[0]<<" y2[0]="<<y2[0]<<"\n";
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    //cout <<"x[i]="<<x[i]<<" y[i]="<<y[i]<<" y2[i]="<<y2[i]<<"\n";
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  //rep.printVec("y2",y2,50);
  delete [] u;
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::splint(const double xa[],
			  const double ya[],
			  const double y2a[],
			  const int n,
			  const double x,
			  double *y
			  )
{
  int klo,khi,k;
  double h,b,a;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
  }
  //cout <<"khi="<<khi<<" klo="<<klo;
  //cout <<"\t\txhi="<<xa[khi]<<" xlo="<<xa[klo];
  //cout <<"\t\tyhi="<<ya[khi]<<" ylo="<<ya[klo];
  //cout <<"\t\ty2hi="<<y2a[khi]<<" y2lo="<<y2a[klo]<<"\n";
  h=xa[khi]-xa[klo];
  //if (h == 0.0) { rep.error("Bad xa input to routine splint",h); }
  if (h < VERY_TINY_VALUE) { rep.error("Bad xa input to routine splint",h); }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}


// ##################################################################
// ##################################################################



void interpolate_arrays::root_find_linear(
        const double *xarr, ///< Array of x values.
        const double *yarr, ///< Array of y values.
        const size_t len,     ///< Array sizes
        const double xreq,  ///< x we are searching for.
        double *yreq  ///< pointer to result.
        )
{
  //
  // Given a vector of x-values, and corresponding y-values, and an input
  // x value, find the corresponding y-value by bisection and then linear
  // interopolation.
  //
  // First we find the two x-points in the array which bracket the requested
  // x-value, with bisection.
  //
  size_t
    ihi = len,  // upper bracketing value
    ilo = 0,    // lower bracketing value
    imid= 0;    // midpoint
  //int count=0;
  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (xarr[imid] < xreq) ilo = imid;
    else                   ihi = imid;
    //count ++;
  } while (ihi-ilo >1);
  //cout.precision(12);
  //cout <<"count="<<count<<", ihi="<<ihi<<" and ilo="<<ilo<<", x[ihi]=";
  //cout <<xarr[ihi]<<" and x[ilo]="<<xarr[ilo]<<" and xreq="<<xreq<<"\n";
  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  *yreq = yarr[ilo] + (yarr[ihi]-yarr[ilo])*(xreq-xarr[ilo])/(xarr[ihi]-xarr[ilo]);
  return;
}


// ##################################################################
// ##################################################################



void interpolate_arrays::spline_vec(
        const std::vector<double> &x,
        const std::vector<double> &y,
        const int n,
        double yp1,
        double ypn,
        std::vector<double> &y2
        )
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double [n];
  
  //
  // Large values in the 4th,5th args tell it to use natural boundary conditions,
  // which means set the second derivative to zero at the endpoints.
  // A small value (<1.0e30) indicates that this is the actual value of the first
  // derivative at the boundary values (4th is lower limit, 5th is upper limit).
  //

  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  //cout <<"x[0]="<<x[0]<<" y[0]="<<y[0]<<" y2[0]="<<y2[0]<<"\n";
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    //cout <<"x[i]="<<x[i]<<" y[i]="<<y[i]<<" y2[i]="<<y2[i]<<"\n";
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    //cout <<"constant gradient\n";
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  //rep.printVec("y2",y2,50);
  delete [] u;
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::splint_vec(
        const std::vector<double> &xa,
        const std::vector<double> &ya,
        const std::vector<double> &y2a,
        const int n,
        const double x,
        double *y
        )
{
  int klo,khi,k;
  double h,b,a;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
  }
  //cout <<"khi="<<khi<<" klo="<<klo;
  //cout <<"\t\txhi="<<xa[khi]<<" xlo="<<xa[klo];
  //cout <<"\t\tyhi="<<ya[khi]<<" ylo="<<ya[klo];
  //cout <<"\t\ty2hi="<<y2a[khi]<<" y2lo="<<y2a[klo]<<"\n";
  h=xa[khi]-xa[klo];
  if (h < 1.0e-150) { rep.error("Bad xa input to routine splint",h); }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}




