///
/// \file interpolate.cpp
///  
/// \brief General purpose class for interpolating tables.
/// \author Jonathan Mackey
/// 
/// Modifications:
/// - 2015.03.04 JM: moved from global.h GeneralStuff class.
/// - 2015.03.23 JM: added bilinear interolation, added bounds
///   checking for linear/bilinear int., fixed bug in linear int.
/// - 2017.07.26 JM: added root_find_bilinear_vec() function.



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/interpolate.h"
#include "constants.h"
#include <iostream>

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



void interpolate_arrays::spline(
      const double *x,
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



void interpolate_arrays::splint(
      const double xa[],
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
  h=xa[khi]-xa[klo];
  if (h < VERY_TINY_VALUE) {
    cerr << "Bad xa input to routine splint: h="<< h<<"\n";
    *y = VERY_LARGE_VALUE;
    return;
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}



// ##################################################################
// ##################################################################



void interpolate_arrays::root_find_linear_vec(
        const vector<double> &xarr, ///< Array of x values.
        const vector<double> &yarr, ///< Array of y values.
        const double xreq,  ///< x we are searching for.
        double &yreq        ///< pointer to result.
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
  size_t len = xarr.size();
  size_t
    ihi = len-1,  // upper bracketing value
    ilo = 0,    // lower bracketing value
    imid= 0;    // midpoint
  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (xarr[imid] < xreq) ilo = imid;
    else                   ihi = imid;
  } while (ihi-ilo >1);

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  double xval=0.0;
  if      (xreq>xarr[ihi]) xval = xarr[ihi];
  else if (xreq<xarr[ilo]) xval = xarr[ilo];
  else                     xval = xreq;

  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  yreq = yarr[ilo] + (yarr[ihi]-yarr[ilo])*(xval-xarr[ilo])/(xarr[ihi]-xarr[ilo]);
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
    ihi = len-1,  // upper bracketing value
    ilo = 0,    // lower bracketing value
    imid= 0;    // midpoint
  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (xarr[imid] < xreq) ilo = imid;
    else                   ihi = imid;
  } while (ihi-ilo >1);

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  double xval=0.0;
  if      (xreq>xarr[ihi]) xval = xarr[ihi];
  else if (xreq<xarr[ilo]) xval = xarr[ilo];
  else                     xval = xreq;

  //
  // Now we linearly interpolate the y value between the two adjacent
  // bracketing points.
  //
  *yreq = yarr[ilo] + (yarr[ihi]-yarr[ilo])*(xval-xarr[ilo])/(xarr[ihi]-xarr[ilo]);
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
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
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
        double &y
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
  h=xa[khi]-xa[klo];
  if (h < VERY_TINY_VALUE) {
    cerr << "Bad xa input to routine splint: h="<< h<<"\n";
    y = VERY_LARGE_VALUE;
    return;
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}


// ##################################################################
// ##################################################################


void interpolate_arrays::root_find_bilinear(
        const double *x,     ///< Array of x values.
        const double *y,     ///< Array of y values.
        double **f,  ///< Array of function values
        const size_t *len,  ///< Array sizes
        const double *xreq, ///< (x,y) we are searching for.
        double *res          ///< pointer to result.
        )
{
  //
  // Given a vector of (x,y)-values, and corresponding f-values, 
  // and an input (x,y) value, find the corresponding f-value by
  // bisection and then bi-linear interopolation.
  //
  // First we find the two x/y-points in the array which bracket the requested
  // x/y-value, with bisection.
  //
  size_t
    ihi = len[0]-1, jhi = len[1]-1,  // upper bracketing value
    ilo = 0,      jlo=0,    // lower bracketing value
    imid= 0,      jmid=0;    // midpoint
  int count=0;
  double xval=0.0, yval=0.0;

  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (x[imid] < xreq[0]) ilo = imid;
    else                   ihi = imid;
    count ++;
  } while (ihi-ilo >1);
  cout.precision(6);

  count=0;
  do {
    jmid = jlo + floor((jhi-jlo)/2.0);
    if (y[jmid] < xreq[1]) jlo = jmid;
    else                   jhi = jmid;
    count ++;
  } while (jhi-jlo >1);

  if (ihi-ilo != 1) {
    cerr <<"root_find_bilinear: Couldn't bracket root i: ";
    cerr <<ihi<<", "<<ilo<<"\n";
  }

  if (jhi-jlo != 1) {
    cerr <<"root_find_bilinear: Couldn't bracket root j: ";
    cerr <<jhi<<", "<<jlo<<"\n";
  }

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  if      (xreq[0]>x[ihi]) xval = x[ihi];
  else if (xreq[0]<x[ilo]) xval = x[ilo];
  else                     xval = xreq[0];
  if      (xreq[1]>y[jhi]) yval = y[jhi];
  else if (xreq[1]<y[jlo]) yval = y[jlo];
  else                     yval = xreq[1];

  //
  // Now we use bilinear interpolation to get the result
  // f(x,y)=( f(lo,lo)(xhi-x)(yhi-y)+
  //          f(hi,lo)(x-xlo)(yhi-y)+
  //          f(lo,hi)(xhi-x)(y-ylo)+
  //          f(hi,hi)(x-xlo)(y-ylo) )/(dx*dy)
  //
  *res = (f[ilo][jlo] * (x[ihi]-xval) * (y[jhi]-yval) +
          f[ihi][jlo] * (xval-x[ilo]) * (y[jhi]-yval) +
          f[ilo][jhi] * (x[ihi]-xval) * (yval-y[jlo]) +
          f[ihi][jhi] * (xval-x[ilo]) * (yval-y[jlo]));
  *res /= ( (x[ihi]-x[ilo]) * (y[jhi]-y[jlo]) );
  return;
}


// ##################################################################
// ##################################################################


double interpolate_arrays::root_find_bilinear_vec(
        const vector<double> &x,     ///< Array of x values.
        const vector<double> &y,     ///< Array of y values.
        const vector<vector<double> > &f,  ///< Array of function values
        const vector<size_t> &len,  ///< Array sizes
        const vector<double> &xreq ///< (x,y) we are searching for.
        )
{
  //
  // Given a vector of (x,y)-values, and corresponding f-values, 
  // and an input (x,y) value, find the corresponding f-value by
  // bisection and then bi-linear interopolation.
  //
  // First we find the two x/y-points in the array which bracket the requested
  // x/y-value, with bisection.
  //
  size_t
    ihi = len[0]-1, jhi = len[1]-1,  // upper bracketing value
    ilo = 0,      jlo=0,    // lower bracketing value
    imid= 0,      jmid=0;    // midpoint
  int count=0;
  double xval=0.0, yval=0.0, result=0.0;

  do {
    imid = ilo + floor((ihi-ilo)/2.0);
    if (x[imid] < xreq[0]) ilo = imid;
    else                   ihi = imid;
    count ++;
  } while (ihi-ilo >1);
  cout.precision(6);
  //cout <<"count="<<count<<", ihi="<<ihi<<" and ilo="<<ilo<<", x[ihi]=";
  //cout <<x[ihi]<<" and x[ilo]="<<x[ilo]<<" and xreq="<<xreq[0];

  count=0;
  do {
    jmid = jlo + floor((jhi-jlo)/2.0);
    if (y[jmid] < xreq[1]) jlo = jmid;
    else                   jhi = jmid;
    count ++;
  } while (jhi-jlo >1);
  //cout <<"count="<<count<<", jhi="<<jhi<<" and jlo="<<jlo<<", y[jhi]=";
  //cout <<y[jhi]<<" and y[jlo]="<<y[jlo]<<" and yreq="<<xreq[1]<<"\n";

  if (ihi-ilo != 1) {
    cerr <<"root_find_bilinear: Couldn't bracket root i: ";
    cerr <<ihi<<", "<<ilo<<"\n";
  }

  if (jhi-jlo != 1) {
    cerr <<"root_find_bilinear: Couldn't bracket root j: ";
    cerr <<jhi<<", "<<jlo<<"\n";
  }

  //
  // Array bounds checking: if we are extrapolating do it
  // with zero slope (take the edge value).
  //
  if      (xreq[0]>x[ihi]) xval = x[ihi];
  else if (xreq[0]<x[ilo]) xval = x[ilo];
  else                     xval = xreq[0];
  //cout <<", xval="<<xval<<"\n";
  if      (xreq[1]>y[jhi]) yval = y[jhi];
  else if (xreq[1]<y[jlo]) yval = y[jlo];
  else                     yval = xreq[1];

  //
  // Now we use bilinear interpolation to get the result
  // f(x,y)=( f(lo,lo)(xhi-x)(yhi-y)+
  //          f(hi,lo)(x-xlo)(yhi-y)+
  //          f(lo,hi)(xhi-x)(y-ylo)+
  //          f(hi,hi)(x-xlo)(y-ylo) )/(dx*dy)
  //
  result = (f[ilo][jlo] * (x[ihi]-xval) * (y[jhi]-yval) +
          f[ihi][jlo] * (xval-x[ilo]) * (y[jhi]-yval) +
          f[ilo][jhi] * (x[ihi]-xval) * (yval-y[jlo]) +
          f[ihi][jhi] * (xval-x[ilo]) * (yval-y[jlo]));
  result /= ( (x[ihi]-x[ilo]) * (y[jhi]-y[jlo]) );
  return result;
}


// ##################################################################
// ##################################################################


double interpolate_arrays::root_find_trilinear_vec(
    const vector<double> &x_vec,     ///< Array of x values.
    const vector<double> &y_vec,     ///< Array of y values.
    const vector<double> &z_vec,     ///< Array of z values.
    const vector<vector<vector<double> > > &f,  ///< Array of function values f(x,y,z)
    const vector<size_t> &len,  ///< Array sizes
    const vector<double> &xreq ///< (x,y,z) we are searching for.
    )
{
  //
  // Determine (x0,y0,z0) and (x1,y1,z1) - nearest neighbours to input (x,y,z)
  //

  // Initialise nearest neighbour variables and x, y, z
  double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
  double x = xreq[0], y = xreq[1], z = xreq[2];
  
  // Loops to determine the value of the nearest neighbours of x, y and z
  int x_index = 0, y_index = 0, z_index = 0;
  
  while(x > x_vec[x_index]) x_index++;
  while(y > y_vec[y_index]) y_index++;
  while(z > z_vec[z_index]) z_index++;
  
  // Set nearest neighbours
  x0 = x_vec[x_index - 1], x1 = x_vec[x_index];
  y0 = y_vec[y_index - 1], y1 = y_vec[y_index];
  z0 = z_vec[z_index - 1], z1 = z_vec[z_index];
  
  //
  // Calculate delta x, delta y and delta z terms for trilinear interpolation
  //
  double dx = (x - x0)/(x1 - x0);
  double dy = (y - y0)/(y1 - y0);
  double dz = (z - z0)/(z1 - z0);
  
  //
  // Calculate f000, f001 etc. terms for coefficients
  //
  double f000 = f[x_index - 1][y_index - 1][z_index - 1]; // f(x0, y0, z0)
  
  double f001 = f[x_index - 1][y_index - 1][z_index]; // f(x0, y0, z1)
  double f010 = f[x_index - 1][y_index][z_index - 1]; // f(x0, y1, z0)
  double f100 = f[x_index][y_index - 1][z_index - 1]; // f(x1, y0, z0)
  
  double f110 = f[x_index][y_index][z_index - 1]; // f(x1, y1, z0)
  double f011 = f[x_index - 1][y_index][z_index]; // f(x0, y1, z1)
  double f101 = f[x_index][y_index - 1][z_index]; // f(x1, y0, z1)
  
  double f111 = f[x_index][y_index][z_index]; // f(x1, y1, z1)
  
  //
  // Calculate c coefficients for trilinear interpolation
  //
  double c0 = f000;
  double c1 = f100 - f000;
  double c2 = f010 - f000;
  double c3 = f001 - f000;
  double c4 = f110 - f010 - f100 + f000;
  double c5 = f011 - f001 - f010 + f000;
  double c6 = f101 - f001 - f100 + f000;
  double c7 = f111 - f011 - f101 - f110 + f100 + f001 + f010 - f000;
  
  //
  // Return f(x,y,z)
  //
  return c0 + c1*dx + c2*dy + c3*dz + c4*dx*dy + c5*dy*dz + c6*dz*dx + c7*dx*dy*dz;
}


// ##################################################################
// ##################################################################





