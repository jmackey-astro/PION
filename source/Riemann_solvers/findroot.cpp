/// \file findroot.cc
/// 
/// \brief Rootfinding routines for Riemann Solver
/// 
/// \author Jonathan Mackey
/// 
/// Two solve options are available, solve_pos() for when you know the root
/// is positive, i.e. in \f$ [0,\infty]\f$, and solve_pm() when the root could
/// be anywhere.
/// 
/// The rootfinding algorithms implemented so far are:
///    Bisection:  Slow(ish) but very reliable.
///    Brent's Method: from NR, I find it fast and reliable.
/// 
/// - 2009-10-23 made Riemann_Euler inherit from findroot, and
///  redefine the function to get the root of!  So now findroot
///  can't be set up as is, a derived class needs to define the
///  function to get the root of.
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2010.12.23 JM: Moved to Riemann_solver/ directory.
/// - 2015.03.10 JM: Tidied up a bit.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)
///

#include "findroot.h"
#include "constants.h"
#include <iostream>
using namespace std;



// ##################################################################
// ##################################################################




findroot::findroot()
{
  findroot::errtol = 1.0e-8;
  // This is the fractional accuracy we want to find the root to.
}


// ##################################################################
// ##################################################################



findroot::~findroot()
{
}


// ##################################################################
// ##################################################################



// This version is specific to solving the equation in the Exact Riemann Solver.
// Other versions could be specified, taking more or fewer parameters.
int findroot::FR_find_root(
      pion_flt *ans, ///< pointer to result
      const pion_flt p1,  ///< parameter 1
      const pion_flt p2,  ///< parameter 2
      const pion_flt p3,  ///< parameter 3
      const pion_flt p4,  ///< parameter 4
      const pion_flt p5   ///< parameter 5
      )
{

  FR_param1 = p1;  ///< e.g. left state pressure
  FR_param2 = p2;  ///< e.g. right state pressure
  FR_param3 = p3;  ///< e.g. left state sound speed
  FR_param4 = p4;  ///< e.g. right state sound speed
  FR_param5 = p5;  ///< e.g. EOS gamma

  // Set initial values for x1, x2
  // Modify this for your equations
  //  pion_flt x1 = 1.;
  //  pion_flt x2 = 2.;
  // My initial guesses are 1/3 of and 3 times the arithmetic mean 
  // of the left and right pressures.
  pion_flt x1 = (FR_param1+ FR_param2)/6.0;
  pion_flt x2 = x1*9.0;
  // This guess is just [0,1], which is not a bad starting point.
  //  pion_flt x1 = 0.;
  //  pion_flt x2 = 1.;
 
  // Call the common solver, now that parameters are set properly.
  int err = findroot::solve_pos(x1, x2, ans);
  if (err!=0) {
    cerr << "(findroot::solve_riemann) solve_pos exited abnormally" << "\n";
    return(1);
  }
  //  cout << "(findroot::solve_riemann) Success: ans = " << *ans << "\n";
  return(0);
}


// ##################################################################
// ##################################################################



// This version is for testing
int findroot::solve_test(
        pion_flt *ans, ///< pointer to result
	const pion_flt p1   ///< pointer to parameter data.
        )
{
  //
  // assign parameters which the test function uses (here only one function).
  // 
  FR_param1 = p1;
  
  //
  // Set initial values for x1, x2
  // 
  pion_flt x1 = 0.001;
  pion_flt x2 = 1.;
  
  //
  // Call the common solver, now that parameters are set properly.
  //
  int err = findroot::solve_pos(x1, x2, ans);
  if (err!=0) {
    cerr << "(findroot::solve_test) exited abnormally, must be a bug!" << "\n";
    return(1);
  }
  cout << "(findroot::solve_test) Success: ans = " << *ans << "\n";
  return(0);
}


// ##################################################################
// ##################################################################



pion_flt findroot::FR_test_function(const pion_flt x)
{
  // f(x) = x^2-gamma, where gamma is a parameter passed into the solver.
  //   return(x*x-(*gamma));
  
  // use this for solving a 2D stromgen sphere (no recombs) with 1/r profile.
  return x-log(1.+x)-(FR_param1);
}



// ##################################################################
// ##################################################################



int findroot::solve_pos(
      pion_flt x1,
      pion_flt x2,
      pion_flt *ans
      )
{
  //  cout << "(fr::solve) ans = " << *ans << "\n";
  int err = bracket_root_pos(&x1,&x2);
  if (err!=0) {
    cerr << "(findroot::solve_pos) bracket exited abnormally" << "\n";
    *ans = -1.0;
    return(1);
  }
  //  cout <<"(fr::solve) Bracketed Root: x1="<<x1<<" and x2="<<x2<<"\n";
  //err = find_root_bisection(&x1,&x2,errtol,ans);
  err = find_root_zbrent(x1,x2,errtol,ans);
  if (err!=0) {
     cerr << "(findroot::solve_pos) couldn't find root in range [0, 1e10]" << "\n";
    *ans = -1.0;
    return(1);
  }
  //  cout << "(fr::solve) Got answer: p* = " << *ans << "\n";
  return(0);
}


// ##################################################################
// ##################################################################



int findroot::solve_pm(
      pion_flt x1,
      pion_flt x2,
      pion_flt *ans
      )
{
   //  cout << "(fr::solve) ans = " << *ans << "\n";
   //  int err = bracket_root_pos(&x1,&x2);
  int err = bracket_root_pm(&x1,&x2);
  if (err!=0) {
     //cerr << "(findroot::solve_pos) bracket exited abnormally" << "\n";
    *ans = -1.0;
    return(1);
  }
  err = find_root_bisection(&x1,&x2,errtol,ans);
  //err = find_root_zbrent(x1,x2,errtol,ans);
  if (err!=0) {
     //cerr << "(findroot::solve_pos) couldn't find root in range [0, 1e10]" << "\n";
    *ans = -1.0;
     return(1);
  }
  //  cout << "(fr::solve) ans = " << *ans << "\n";
  return(0);
}
 

// ##################################################################
// ##################################################################


 

/// Brackets a root in the range \f$ [-\infty, \infty]\f$
/// 
/// This function returns brackets around a root.
///
int findroot::bracket_root_pm(
      pion_flt *x1,
      pion_flt *x2
      )
{
   // This is based on the NR root bracketing procedure.
   // It works for functions of x=[-infty,infty].
   float factor=0.2;
   if (*x1 == *x2) {
      cerr << "(bracket) error -- x1,x2 are the same." << "\n";
      return(1);
   }
   if (*x1 > *x2) {
      //      cout << "Reordering... setting x1<x2" << "\n";
      pion_flt temp = *x1;
      *x1 = *x2;
      *x2 = temp;
   }
   pion_flt f1 = FR_root_function(*x1);
   pion_flt f2 = FR_root_function(*x2);
   //   cout << "F1,F2 " << f1 << ", " << f2 << "\n";
   for (int j=0; j<50; j++) {
      //cout << "(bracket) have done " << j << " expansions. (x1,x2)= " << *x1 << ", " << *x2 << "\n";
      //cout << "F1,F2 " << f1 << ", " << f2 << "\n";
      if (f1*f2 < 0) {
	 //cout << "(bracket) Root bracketed after " << j << " expansions. (x1,x2)= " << *x1 << ", " << *x2 << "\n";
	 //cout << "F1,F2 " << f1 << ", " << f2 << "\n";
	 return(0);
      }
      if (fabs(f1) < fabs(f2))
	 f1 = FR_root_function(*x1 = *x1-factor*(*x2-*x1));
      else
	 f2 = FR_root_function(*x2 = *x2+factor*(*x2-*x1));
   }
   //cerr << "(bracket) Error -- couldn't bracket root after 50 iterations.  Failing" << "\n";
   *x1=*x2=0.;
   return(1);
}


// ##################################################################
// ##################################################################



int findroot::bracket_root_pos(pion_flt *x1, pion_flt *x2)
{
  // This is based on the NR root bracketing procedure.
  // It only works for fucntion of x=[0,infty], and would need
  // to be modified to handle negative values of x.
  float factor=1.6;
  if (*x1 == *x2) {
    cerr << "(bracket) error -- x1,x2 are the same." << "\n";
    return(1);
  }
  if (*x1 > *x2) {
     //    cout << "Reordering... setting x1<x2" << "\n";
    pion_flt temp = *x1;
    *x1 = *x2;
    *x2 = temp;
  }
  pion_flt f1 = FR_root_function(*x1);
  pion_flt f2 = FR_root_function(*x2);
  //  cout << "F1,F2 " << f1 << ", " << f2 << "\n";
  for (int j=0; j<50; j++) {
    if (f1*f2 < 0) {
      //      cout << "(bracket) Root bracketed after " << j << " expansions. (x1,x2)= " << *x1 << ", " << *x2 << "\n";
      //      cout << "F1,F2 " << f1 << ", " << f2 << "\n";
      return(0);
    }
    if (fabs(f1) < fabs(f2))
      f1 = FR_root_function(*x1 *= 1./factor);
    else
      f2 = FR_root_function(*x2 *= factor);
  }
  // If it didn't work, try x1=0.  (note this routine only works for f(x>=0).
  f1 = FR_root_function(*x1 =0.);
  if (f1*f2 < 0) {
    //    cout << "(bracket) Root bracketed at x1=0. (x1,x2)= " << *x1 << ", " << *x2 << "\n";
    return(0);
  }
  cerr << "(bracket) Error -- couldn't bracket root after 50 iterations.  Failing (x1,x2)= " << *x1 << ", " << *x2  << "\n";
  *x1=*x2=0.;
  return(1);
}
  

// ##################################################################
// ##################################################################




int findroot::find_root_bisection(pion_flt *x1, pion_flt *x2, pion_flt err, pion_flt *ans)
{
  // This is a simple bisection routine to find the root to an accuracy
  // of 'err'.  It is based on the NR bisection routine.
  if (*x1>*x2) {
    cerr << "(findroot) Error, must have x1<x2, ordered brackets!" << "\n";
    return(1);
  }
  pion_flt f1 =  FR_root_function(*x1);
  pion_flt f2 =  FR_root_function(*x2);
  pion_flt xmid = (*x2+*x1)/2.0;
  pion_flt fmid=0.0;
  for (int j=0; j<100; j++) {
    if (fabs((*x1 - *x2)/ (*x2)) < err) {
      //      cout << "(find_root) Root found after " << j << " bisections to accuracy " << err << "\n";
      *ans = xmid;
      return(0);
    }
    fmid =  FR_root_function(xmid);
    if (f1<f2)
      (fmid >0) ? (*x2 = xmid) : (*x1 = xmid);
    else
      (fmid >0) ? (*x1 = xmid) : (*x2 = xmid);
    xmid = (*x2+*x1)/2.0;
  }
  cerr << "(findroot) Couldn't find root after 50 bisections, failing" << "\n";
  return(1);
}


// ##################################################################
// ##################################################################



int findroot::find_root_zbrent(
      pion_flt x1,
      pion_flt x2,
      pion_flt tol,
      pion_flt *ans)
{
   // This is the numerical recipes routine 'zbrent.c'
   // It uses bisection and a higher order method when it can.
   // See NR book, p.361 of ansi-c version.
   int ITMAX=100;
   pion_flt EPS=MACHINEACCURACY;
   int iter;
   pion_flt a=x1,b=x2,c=x2,d,e,min1,min2;
   pion_flt fa=FR_root_function(a),fb=FR_root_function(b),fc,p,q,r,s,tol1,xm;
  d=0.;e=0.;   
   if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
      cerr<<"Root must be bracketed in zbrent"<<"\n"; return(1);
   }
   fc=fb;
   for (iter=1;iter<=ITMAX;iter++) {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
	 c=a;
	 fc=fa;
	 e=d=b-a;
      }
      if (fabs(fc) < fabs(fb)) {
	 a=b;
	 b=c;
	 c=a;
	 fa=fb;
	 fb=fc;
	 fc=fa;
      }
     //     cout <<"fabs(b)="<<fabs(b)<<"\t"<<fabs(c)<<"\t"<<fabs(a)<<"\n";
     // This tolerance is first a relative accuracy, where EPS is roughly
     // the machine precision, and fabs(b) is the value of the pressure
     // at the lower or upper bound.
     // The second term was initially an absolute tolerance, but I don't 
     // want to have this, as I know my function is always positive, but
     // can sometimes be very small -- e.g. PG=1.e-11 is typical for cgs
     // units.  So I changed it to a relative tolerance too.  'tol' is 
     // usually larger than EPS, set to something like 1.e-8 as this is
     // still pretty good accuracy.
      tol1=2.0*EPS*fabs(b) +0.5*tol*fabs(b);
      xm=0.5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0.0) {
	*ans= b;
	//	cout <<"ITER="<<iter<<"\n";
	return(0);
      }
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
	 s=fb/fa;
	 if (a == c) {
	    p=2.0*xm*s;
	    q=1.0-s;
	 } else {
	    q=fa/fc;
	    r=fb/fc;
	    p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	    q=(q-1.0)*(r-1.0)*(s-1.0);
	 }
	 if (p > 0.0) q = -q;
	 p=fabs(p);
	 min1=3.0*xm*q-fabs(tol1*q);
	 min2=fabs(e*q);
	 if (2.0*p < (min1 < min2 ? min1 : min2)) {
	    e=d;
	    d=p/q;
	 } else {
	    d=xm;
	    e=d;
	 }
      } else {
	 d=xm;
	 e=d;
      }
      a=b;
      fa=fb;
      if (fabs(d) > tol1)
	 b += d;
      else
	 b += ((xm) >= 0.0 ? fabs(tol1) : -fabs(tol1));
      fb=FR_root_function(b);
   }
   cerr<<"Maximum number of iterations exceeded in zbrent"<<"\n";
   return(1);
}


// ##################################################################
// ##################################################################



//  For testing the findroot class
// int main (int argc, char **argv)
// {
//    if(argc!=2) {
//       cerr<< "Main: pleas pass 'pion_flt g' as argument."<< "\n";
//       return(1);
//    }
//    findroot fr;
//    int err=0;
//    pion_flt gamma = 5./3.;
//    pion_flt left[3]  = {1.0, 0., 10.};
//    pion_flt right[3] = {0.125, 0., 1.};
//    pion_flt cl = sqrt(gamma*left[2]/left[0]);
//    pion_flt cr = sqrt(gamma*right[2]/right[0]);
//    pion_flt ans[3] = {0.,0.,0.};
//    fr.solve_riemann(ans, left, right, &gamma, &cl, &cr);
//    cout << "ans: p=" << ans[0] << " " << ans[1] << " " << ans[2] << "\n";
//    cout <<" ****************************************************************** "<<"\n";
//    left[1] = -4.2;
//    fr.solve_riemann(ans, left, right, &gamma, &cl, &cr);
//    cout << "ans: p=" << ans[0] << " " << ans[1] << " " << ans[2] << "\n";
//    cout <<" ****************************************************************** "<<"\n";
//    left[1] = -0.2;
//    left[2] = 2500000;
//    fr.solve_riemann(ans, left, right, &gamma, &cl, &cr);
//    cout << "ans: p=" << ans[0] << " " << ans[1] << " " << ans[2] << "\n";
//    cout <<" ****************************************************************** "<<"\n";
//    cout <<" ****************************************************************** "<<"\n";
//    pion_flt vmin=-23.;
//    pion_flt vmax= 25.;
//    int Ng = 50;
//    left[2] = 10.;
//    for (int k=0;k<Ng;k++) { 
//       left[1] = (vmin+ k*(vmax-vmin)/(Ng-1));
//       //      cout <<" left: {d,u,p} " << left[0] << " " << left[1] << " " << left[2] << "\n";
//       //      cout <<"right: {d,u,p} " <<right[0] << " " <<right[1] << " " <<right[2] << "\n";
//       err = fr.solve_riemann(ans, left, right, &gamma, &cl, &cr);
//       //      cout <<err<<"\n";
//       if (err!=0) {
// 	 cout <<"Solve failed: left: {d,u,p} " << left[0] << " " << left[1] << " " << left[2] << "\n";
// 	 cout <<"Solve failed:right: {d,u,p} " <<right[0] << " " <<right[1] << " " <<right[2] << "\n";
//       }
//       cout << "ans: p=" << ans[0] << " " << ans[1] << " " << ans[2] << "\n";
//       //      cout <<" ****************************************************************** "<<"\n";
//    } 
//    cout <<" ****************************************************************** "<<"\n";
//    cout <<" ****************************************************************** "<<"\n";
//    gamma = atof(argv[1]);
//    fr.solve_test(ans,&gamma);
//    cout << "Test g=" <<gamma<< " Ans=" << ans[0] << " " << ans[1] << " " << ans[2] << "\n";
//    cout <<" ****************************************************************** "<<"\n";
//    cout <<" ****************************************************************** "<<"\n";
//    time_t t1,t2;
//    t1 = time(0);
//    int i=0;
//    for (i=0;i<1000000;i++) {
//       gamma = (static_cast<pion_flt>(i)+0.1)/1.e4;
//       err = fr.solve_test(ans,&gamma);
//       gamma = 5./3.;
//       err+= fr.solve_riemann(ans, left, right, &gamma, &cl, &cr);
//       if (err!=0) {
// 	 cout <<"Solve failed: gamma=" << gamma << "  ans="<<ans[0]<<"\n";
//       }
//    }
//    t2 = time(0);
//    pion_flt d;
//    d =difftime(t2,t1);
//    cout << "Timing completed" << "\n";
//    cout << "Nt: " << i << "\t solver took " << d << " seconds. "<< "\n";
//    return(0);
// }

