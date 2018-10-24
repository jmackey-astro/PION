///
/// \file VectorOps.cc
/// 
/// \author Jonathan Mackey
/// 
/// Function definitions for class members of the various VectorOps classes.
/// These work for a grid in which the cells are all the same size, and either 
/// square or cube-shaped for 2 and 3D respectively.
/// 
/// Modified:\n
///  - 2007-08-01 File Created
///  - 2007-08-07 Axi-symmetry working (at for Euler and MHD equations).
///  - 2007-11-30 Curl function working in cartesian coords.
///  - 2009-09    Added MinMod slope limiter option as an ifdef (TEMPORARY)
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.12.04 JM: Added constructor with only one argument.  Also
///   a set_dx() function.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.13 JM: Modified for new code structure; added the grid
///    pointer everywhere.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "constants.h"

#include "VectorOps.h"

using namespace std;

BaseVectorOps::~BaseVectorOps() {}

#define AVG_FALLE
//#define AVG_MINMOD

// ##################################################################
// ##################################################################

double BaseVectorOps::AvgFalle(
        const double a,
        const double b
        )
{
  // If the slopes have different signs, then set the slope to zero.
  // If the slopes are both very small, then set slope to zero.
  if (a*b <= VERY_TINY_VALUE) return(0.0);
  // If the slopes are both very small, then set slope to zero.
  //  else if ( temp <= 2*MACHINEACCURACY) return(0.);
#if   defined AVG_FALLE
  // If we get to here, then they are both the same sign, and finite, so we do 
  // the proper averaging.
  else return(a*b*(a+b)/(a*a+b*b));
#elif defined AVG_MINMOD
  double r=a/b;
  return (r>0.0) ? min(r,1.0)*b : 0.0; 
#else
#error "MUST DEFINE EITHER AVG_FALLE OR AVG_MINMOD"
#endif
}


// ##################################################################
// ##################################################################


double BaseVectorOps::DotProduct(
        const pion_flt *v1, // Vector 1.
        const pion_flt *v2, // Vector 2.
        const int n       // length of vectors.
        )
{
  double temp=0.;
  for (int i=0;i<n;i++) temp += v1[i]*v2[i];
  return temp;
}


// ##################################################################
// ##################################################################


int BaseVectorOps::CrossProduct(
        const pion_flt *a, // Vector 1.
        const pion_flt *b, // Vector 2.
        const int n,      // length of vectors.
        pion_flt *ans       // Result vector
        )
{
  if      (n==3) {
    ans[0] = a[1]*b[2] - a[2]*b[1];
    ans[1] = a[2]*b[0] - a[0]*b[2];
    ans[2] = a[0]*b[1] - a[1]*b[0];
  }
  else if (n==2) {
    *ans   = a[0]*b[1] - a[1]*b[0];
  }
  else {
    cerr <<"BaseVectorOps::CrossProduct() not 2 or 3d vector!\n";
    return(1);
  }
  return 0; 
}
  
// ##################################################################
// ##################################################################


/// *****************************************************************
/// *************  CARTESIAN COORDINATES X,Y,Z **********************
/// *****************************************************************


// ##################################################################
// ##################################################################


VectorOps_Cart::VectorOps_Cart(int n)
  : VOnd(n)
{
#ifdef TESTING
  cout <<"Setting up VectorOpsCart with ndim="<<VOnd<<"\n";
#endif
  if (VOnd>3) rep.error("Can't do more than 3D simulations!",VOnd);
}



// ##################################################################
// ##################################################################



VectorOps_Cart::~VectorOps_Cart() {}



// ##################################################################
// ##################################################################



double VectorOps_Cart::CellVolume(
      const cell *c,
      const double dR
      )
{
  return dR*dR*dR;
}



// ##################################################################
// ##################################################################



double VectorOps_Cart::CellInterface(
      const cell *,
      const direction,
      const double dR ///< cell diameter
      )
{
  return dR*dR;
}



// ##################################################################
// ##################################################################



double VectorOps_Cart::max_grad_abs(
        const cell *cpt,
        const int sv,
        const int var,
        class GridBaseClass *grid
        )
{
  double VOdx = grid->DX();
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(cpt,static_cast<direction>(i)))
      rep.error("VectorOps_Cart::max_grad_abs: Some neighbour cells don't exist",i);
#endif //TESTING
  
  double grad=0, temp=0;
  switch (sv) {
   case 0: // Use vector c->P
    for (int i=0;i<2*VOnd;i++) {
      temp = fabs( (grid->NextPt(cpt,static_cast<direction>(i)))->P[var] - cpt->P[var])/VOdx;
      if(temp>grad) grad = temp;
    }
    break;
   case 1: // Use Vector c-Ph
    for (int i=0;i<2*VOnd;i++) {
      temp = fabs( (grid->NextPt(cpt,static_cast<direction>(i)))->Ph[var] - cpt->Ph[var])/VOdx;
      if(temp>grad) grad = temp;
    }
    break;
   default:
    rep.error("Don't know what state vector to use for calculating gradient.",sv);
  }

  return(grad);
} // max_grad_abs


// ##################################################################
// ##################################################################


void VectorOps_Cart::Gradient(
        const cell *c,
        const int sv,
        const int var,
        class GridBaseClass *grid,
        pion_flt *grad
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cart::Grad: Some neighbour cells don't exist",i);
#endif //TESTING
  double VOdx = grid->DX();
  
  switch (sv) {
   case 0: // Use vector c->P
               grad[0] = (grid->NextPt(c,XP)->P[var] - grid->NextPt(c,XN)->P[var])/(2.*VOdx);
    if(VOnd>1) grad[1] = (grid->NextPt(c,YP)->P[var] - grid->NextPt(c,YN)->P[var])/(2.*VOdx);
    if(VOnd>2) grad[2] = (grid->NextPt(c,ZP)->P[var] - grid->NextPt(c,ZN)->P[var])/(2.*VOdx);
    break;
   case 1: // Use Vector c-Ph
               grad[0] = (grid->NextPt(c,XP)->Ph[var] - grid->NextPt(c,XN)->Ph[var])/(2.*VOdx);
    if(VOnd>1) grad[1] = (grid->NextPt(c,YP)->Ph[var] - grid->NextPt(c,YN)->Ph[var])/(2.*VOdx);
    if(VOnd>2) grad[2] = (grid->NextPt(c,ZP)->Ph[var] - grid->NextPt(c,ZN)->Ph[var])/(2.*VOdx);
    break;
   default:
    rep.error("Don't know what state vector to use for calculating gradiend.",sv);
  }
  return;
}



// ##################################################################
// ##################################################################



double VectorOps_Cart::CentralDiff(
      class GridBaseClass *grid,  ///< pointer to computational grid.
      class cell *c, ///< pointer to cell
      const int ax,    ///< axis along which to take difference
      const int sv,    ///< Which vector to take values from (P=0,Ph=1,dU=2)
      const int ii    ///< index in state vector of variable
      )
{
  cell *cn, *cp;
  enum direction ndir = static_cast<direction>(2*ax);
  enum direction pdir = static_cast<direction>(2*ax+1);

  cn = (grid->NextPt(c,ndir)) ? grid->NextPt(c,ndir) : c;
  cp = (grid->NextPt(c,pdir)) ? grid->NextPt(c,pdir) : c;

  double ans=0.0;
  switch (sv) {

    case 0:
    ans = (cp->P[ii]-cn->P[ii]);
    break;

    case 1:
    ans = (cp->Ph[ii]-cn->Ph[ii]);
    break;

    default:
    rep.error("state vector for calculating CentralDiff.",sv);
  }

  return ans;
}



// ##################################################################
// ##################################################################


double VectorOps_Cart::Divergence(
        cell *c, 
        const int sv, 
        const int *var,
        class GridBaseClass *grid
        )
{ // get divergence of vector quantity.
  
  cell *ngb[6];
  double dx[3];
  double VOdx = grid->DX();
  enum direction ndir;
  enum direction pdir;

  // Assign pointers to neighbouring cells if they exist.
  // If not, then use the current cell and do a one-sided calc.
  // If we are at the edge of the grid, then we don't really care
  // what value we obtain because this is a outer boundary zone.
  for (int v=0;v<VOnd;v++) {
    ndir = static_cast<direction>(2*v);
    pdir = static_cast<direction>(2*v+1);
    ngb[ndir] = (grid->NextPt(c,ndir)) ? grid->NextPt(c,ndir) : c;
    ngb[pdir] = (grid->NextPt(c,pdir)) ? grid->NextPt(c,pdir) : c;
    if (ngb[ndir]==c || ngb[pdir]==c) dx[v]=VOdx;
    else                              dx[v]=2.0*VOdx;
  }
  //cout <<"c="<<c<<" : ";
  //rep.printVec("ngb",ngb,2*VOnd);
  //rep.printVec("dx",dx,VOnd);

  double divv=0.0;
  switch (sv) {
    case 0: // Use vector c->P
    for (int v=0;v<VOnd;v++) {
      ndir = static_cast<direction>(2*v);
      pdir = static_cast<direction>(2*v+1);
      divv += (ngb[pdir]->P[var[v]] - ngb[ndir]->P[var[v]]) / dx[v];
    }
    break;

    case 1: // Use Vector c-Ph
    for (int v=0;v<VOnd;v++) {
      ndir = static_cast<direction>(2*v);
      pdir = static_cast<direction>(2*v+1);
      divv += (ngb[pdir]->Ph[var[v]] - ngb[ndir]->Ph[var[v]]) /dx[v];
#ifdef TESTING
      if (!isfinite(divv)) {
        cout <<"divv["<<v<<"] = ";
        cout <<(ngb[pdir]->Ph[var[v]] - ngb[ndir]->Ph[var[v]]) /dx[v];
        cout <<": dx="<<dx[v]<<" cn="<<ngb[ndir]<<", cp="<<ngb[pdir];
        cout <<": cp-ph="<<ngb[pdir]->Ph[var[v]];
        cout <<", cn-ph="<<ngb[ndir]->Ph[var[v]];
        cout <<": var="<<var[v]<<"\n";
      }
#endif
    }
    break;

    default:
    rep.error("state vectorfor calculating divergence.",sv);
  }
  return(divv);
} // Div


// ##################################################################
// ##################################################################


void VectorOps_Cart::Curl(
        const cell *c, 
        const int vec, 
        const int *var, 
        class GridBaseClass *grid,
        pion_flt *ans
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cart::Curl: Some neighbour cells don't exist",i);
#endif //TESTING
  if (!c->isgd) {
    //cout <<"curl of non-grid-cell, returning 0";
    ans[0] = ans[1] = ans[2] =0.0;
    return;
    //rep.error("Not Grid Cell! can't get curl. id:",c->id);
  }
  double VOdx = grid->DX();

  int vx=var[0], vy=var[1], vz=var[2];
  pion_flt *vxp=0, *vxn=0, *vyp=0, *vyn=0, *vzp=0, *vzn=0;
  switch (vec) {
   case 0: // c->P
//    cout <<"using P.\n";
    vxp = grid->NextPt(c,XP)->P; vxn = grid->NextPt(c,XN)->P;
    if (VOnd>1) {
      vyp = grid->NextPt(c,YP)->P;
      vyn = grid->NextPt(c,YN)->P;
    }
    if (VOnd>2) {
      vzp = grid->NextPt(c,ZP)->P;
      vzn = grid->NextPt(c,ZN)->P;
    }
    break;
   case 1: // c->Ph
//    cout <<"using Ph = ["<<c->Ph[vx]<<", "<<c->Ph[vy]<<", "<<c->Ph[vz]<<" ]\n";
    vxp = grid->NextPt(c,XP)->Ph; vxn = grid->NextPt(c,XN)->Ph;
    if (VOnd>1) {
      vyp = grid->NextPt(c,YP)->Ph;
      vyn = grid->NextPt(c,YN)->Ph;
    }
    if (VOnd>2) {
      vzp = grid->NextPt(c,ZP)->Ph;
      vzn = grid->NextPt(c,ZN)->Ph;
    }
    break;
   case 2: // c->dU
//    cout <<"using dU.\n";
    vxp = grid->NextPt(c,XP)->dU; vxn = grid->NextPt(c,XN)->dU;
    if (VOnd>1) {
      vyp = grid->NextPt(c,YP)->dU;
      vyn = grid->NextPt(c,YN)->dU;
    }
    if (VOnd>2) {
      vzp = grid->NextPt(c,ZP)->dU;
      vzn = grid->NextPt(c,ZN)->dU;
    }
    break;
   default:
    rep.error("Which vector to calculate on?  (VecCurl), don't know!",vec);
  }

  ans[0] =   0.0;
  ans[1] = -(vxp[vz] - vxn[vz]);
  ans[2] =  (vxp[vy] - vxn[vy]);
  if (VOnd>1) {
    ans[0] +=  (vyp[vz] - vyn[vz]);
    ans[1] +=   0.0;
    ans[2] += -(vyp[vx] - vyn[vx]);
  }
  if (VOnd>2) {
    ans[0] += -(vzp[vy] - vzn[vy]);
    ans[1] +=  (vzp[vx] - vzn[vx]);
    ans[2] +=  0.0;
  }
  // All differences have the same denominator in cartesian cubic cells:
  ans[0] /= (2.*VOdx);
  ans[1] /= (2.*VOdx);
  ans[2] /= (2.*VOdx);
  return;
}
  

// ##################################################################
// ##################################################################


int VectorOps_Cart::SetEdgeState(
        const cell *c,      ///< Current Cell.
        const direction d,  ///< Add or subtract the slope depending on direction.
        const int nv,       ///< length of state vectors.
        const pion_flt *dpdx, ///< Slope vector.
        pion_flt *edge,       ///< vector for edge state. 
        const int OA,       ///< Order of spatial Accuracy.
        class GridBaseClass *grid
        )
{
  double VOdx = grid->DX();
  switch (OA) {
    
   case OA1: // First Order Spatial Accuracy.
    for (int v=0;v<nv;v++) edge[v] = c->Ph[v];
    break; // O1
    
   case OA2: // Second Order Spatial Accuracy.
    switch (d) {
     case XP: case YP: case ZP:
      for (int v=0;v<nv;v++) edge[v] = c->Ph[v] + dpdx[v]*VOdx/2.;
      break; //XP
     case XN: case YN: case ZN:
      for (int v=0;v<nv;v++) edge[v] = c->Ph[v] - dpdx[v]*VOdx/2.;
      break; //XN
     default:
      cerr<<"\t(SetEdgeState) wrong direction!\n"; 
      return(1);
    } // Which Direction.
    break; // O2
   default:
    cerr<<"\t(SetEdgeState) Only know 1st and 2nd order spatial accuracy!!!\n"; 
    cerr<<"\tOA="<<OA<<"\n";
    return(1);
  } // Order of accuracy.
  return(0);
} // SetEdgeState
  

// ##################################################################
// ##################################################################


int VectorOps_Cart::SetSlope(
        const cell *c, ///< Current Cell.
        const axes d,  ///< Which direction to calculate slope in.
        const int nv,  ///< length of state vectors.
        pion_flt *dpdx,  ///< Slope vector to be written to.
        const int  OA, ///< Order of spatial Accuracy.
        class GridBaseClass *grid
        )
{
  if (OA==OA1){ // first order accurate so zero slope.
    for(int v=0;v<nv;v++) dpdx[v]=0.;
  } // 1st order.
  
  else if (OA==OA2) { // second order spatial accuracy.
    pion_flt slpn[nv], slpp[nv];
    cell *cn=0,*cp=0;
    enum direction dp=NO,dn=NO;
    double dx=grid->DX();
    switch (d) {
     case XX: dp=XP; dn=XN; break;
     case YY: dp=YP; dn=YN; break;
     case ZZ: dp=ZP; dn=ZN; break;
     default: rep.error("Bad direction in SetSlope",d);
    }
    cp = grid->NextPt(c,dp); cn=grid->NextPt(c,dn);
//#ifdef TESTING
    if (cp==0 || cn==0) rep.error("No left or right cell in SetSlope",cp);
//#endif //TESTING
    for (int v=0;v<nv;v++) {
      slpn[v] = (c->Ph[v] -cn->Ph[v])/dx;
      slpp[v] = (cp->Ph[v]- c->Ph[v])/dx;
      dpdx[v] = AvgFalle(slpn[v],slpp[v]);
#ifdef FIX_TRACER_SLOPES_TO_DENSITY
      if (v>=SimPM.ftr) {
	dpdx[v] = 0.0;
      }
#endif //FIX_TRACER_SLOPES_TO_DENSITY
    }
  } // 2nd order accurate
  else {
    cerr <<"Error: Only know how to do 1st or 2nd order slope calculation.\n";
    exit(1);
  }
  return(0);
} // SetSlope
  

// ##################################################################
// ##################################################################


int VectorOps_Cart::DivStateVectorComponent(
        const cell *c,    ///< current cell.
        class GridBaseClass *grid,
        const axes d,     ///< current coordinate axis we are looking along.
        const int nv,     ///< length of state vectors.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        pion_flt *dudt      ///< Vector to assign divergence component to.
        )
{
  //enum direction dp;
  //switch (d) {
  // case XX: dp=XP; break;
  // case YY: dp=YP; break;
  // case ZZ: dp=ZP; break;
  // default: rep.error("Bad direction in DivStateVectorComponent",d);
  //}
  
  /** \section Sign
   * Note that this function returns the negative of the i-th component of 
   * the divergence.  This is b/c it is used in the finite volume time update
   * where what is needed is the negative of div(F).
   * */
  double dx=grid->DX();
  for (int v=0;v<nv;v++) {
    dudt[v] = (fn[v]-fp[v])/dx;
  }
  return(0);
} // DivStateVectorComponent


// ##################################################################
// ##################################################################



///******************************************************************
///************  CYLINDRICAL COORDINATES Z,R,THETA ******************
///******************************************************************

// ##################################################################
// ##################################################################



VectorOps_Cyl::VectorOps_Cyl(int n)
  : VectorOps_Cart(n)
{
#ifdef TESTING
  cout <<"Setting up VectorOps_Cyl with ndim="<<VOnd<<"\n";
#endif
  //if (VOnd!=2 && VOnd!=3) rep.error("Why use cylindrical coords in not 2 or 3D?",VOnd);
  if (VOnd>2) rep.warning("VectorOps_Cyl NOT TESTED IN 3D YET",2,VOnd);
}



// ##################################################################
// ##################################################################



VectorOps_Cyl::~VectorOps_Cyl() {}



// ##################################################################
// ##################################################################



double VectorOps_Cyl::CellVolume(
      const cell *c,
      const double dR
      )
{
  /// cell vol = pi * (R+^2 - R-^2) * dz, where dz=dR
  double r = CI.get_dpos(c,Rcyl);
  r = (r+0.5*dR)*(r+0.5*dR) - (r-0.5*dR)*(r-0.5*dR);
  return(M_PI*r*dR);
}




// ##################################################################
// ##################################################################

double VectorOps_Cyl::CellInterface(
      const cell *c, 
      const direction dir,
      const double dR ///< cell diameter
      )
{
  /** \brief Calculation
   * For cylindrical coordinates, different faces have different 
   * surface areas.  
   *  - In the \f$\theta\f$ direction, the area is just
   * \f$\delta A = \delta R \delta z\f$.
   *  - In the \f$z\f$ direction, the area is 
   * \f$\delta A = \pi ((R_i+\frac{\delta R}{2})^2-(R_i-\frac{\delta R}{2})^2)\f$.
   *  - In the positive/negative radial direction, the area is 
   * \f$\delta A = 2\pi \delta z (R_i \pm \frac{\delta R}{2})\f$.
   * */
  //cout <<"CellInterface:"<<dir<<", "<<dR<<"\n";
  double dZ = dR;
  double pos[MAX_DIM];
  CI.get_dpos(c,pos);
  //rep.printVec("pos",pos,2);

  switch (dir) {
   case ZNcyl: case ZPcyl:
    return M_PI*( (pos[Rcyl]+0.5*dR)*(pos[Rcyl]+0.5*dR) - 
                  (pos[Rcyl]-0.5*dR)*(pos[Rcyl]-0.5*dR) );
    break;
   case TNcyl: case TPcyl:
    rep.error("3D cylindrical not implemented in CellInterface",dir);
    break;
   case RNcyl:
    if (pos[Rcyl]>0 && pos[Rcyl]<dR) return 0.0;
    else return 2.0*M_PI*dZ*(pos[Rcyl]-0.5*dR);
    break;
   case RPcyl:
    //cout <<"area: "; rep.printVec("pos",pos,2);
    if (pos[Rcyl]<0 && pos[Rcyl]>-dR) {
      return 0.0;
    }
    else {
      //cout <<"area: "<<dZ<<", "<<pos[Rcyl]<<", "<<dR<<"\n";
      return 2.0*M_PI*dZ*(pos[Rcyl]+0.5*dR);
    }
    break;
   default:
    rep.error("Bad direction to VectorOps_Cyl::CellInterface",dir);
  }
  return(-1.0);
}



// ##################################################################
// ##################################################################


double VectorOps_Cyl::max_grad_abs(
        const cell *c, 
        const int sv, 
        const int var,
        class GridBaseClass *grid
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cyl::max_grad_abs: Some neighbour cells don't exist",i);
#endif //TESTING

  double dT=0.; // physical length R*dTheta
  double dR=grid->DX();
  double dZ=dR;

  double grad=0, temp=0;
  switch (sv) {
   case 0: // Use vector c->P
    // Z-dir
    temp = fabs( grid->NextPt(c,ZPcyl)->P[var] - c->P[var])/dZ;
    if(temp>grad) grad = temp;
    temp = fabs( grid->NextPt(c,ZNcyl)->P[var] - c->P[var])/dZ;
    if(temp>grad) grad = temp;
    if (VOnd>1) { // R-dir
      temp = fabs( grid->NextPt(c,RPcyl)->P[var] - c->P[var])/ (R_com(grid->NextPt(c,RPcyl),dR) - R_com(c,dR));
      if(temp>grad) grad = temp;
      temp = fabs( grid->NextPt(c,RNcyl)->P[var] - c->P[var])/ (R_com(c,dR) - R_com(grid->NextPt(c,RNcyl),dR));
      if(temp>grad) grad = temp;
    }
    if (VOnd>2) { // Theta-dir, need scale factor 1/R
      dT = (CI.get_dpos(grid->NextPt(c,TPcyl),Tcyl) - CI.get_dpos(c,Tcyl))*R_com(c,dR);
      temp = fabs( grid->NextPt(c,TPcyl)->P[var] - c->P[var])/dT;
      if(temp>grad) grad = temp;
      dT = (CI.get_dpos(c,Tcyl) -  CI.get_dpos(grid->NextPt(c,TNcyl),Tcyl))*R_com(c,dR);
      temp = fabs( grid->NextPt(c,TNcyl)->P[var] - c->P[var])/dT;
      if(temp>grad) grad = temp;
    }
    break;
   case 1: // Use Vector c-Ph
    // Z-dir
    temp = fabs( grid->NextPt(c,ZPcyl)->Ph[var] - c->Ph[var])/dZ;
    if(temp>grad) grad = temp;
    temp = fabs( grid->NextPt(c,ZNcyl)->Ph[var] - c->Ph[var])/dZ;
    if(temp>grad) grad = temp;
    if (VOnd>1) { // R-dir
      temp = fabs( grid->NextPt(c,RPcyl)->Ph[var] - c->Ph[var])/ (R_com(grid->NextPt(c,RPcyl),dR) - R_com(c,dR));
      if(temp>grad) grad = temp;
      temp = fabs( grid->NextPt(c,RNcyl)->Ph[var] - c->Ph[var])/ (R_com(c,dR) - R_com(grid->NextPt(c,RNcyl),dR));
      if(temp>grad) grad = temp;
    }
    if (VOnd>2) { // Theta-dir, need scale factor 1/R
      dT = ( CI.get_dpos(grid->NextPt(c,TPcyl),Tcyl) - CI.get_dpos(c,Tcyl))*R_com(c,dR);
      temp = fabs( grid->NextPt(c,TPcyl)->Ph[var] - c->Ph[var])/dT;
      if(temp>grad) grad = temp;
      dT =  (CI.get_dpos(c,Tcyl) - CI.get_dpos(grid->NextPt(c,TNcyl),Tcyl))*R_com(c,dR);
      temp = fabs( grid->NextPt(c,TNcyl)->Ph[var] - c->Ph[var])/dT;
      if(temp>grad) grad = temp;
    }
    break;
   default:
    rep.error("Don't know what state vector to use for calculating gradient.",sv);
  }
  return(grad);
} // max_grad_abs


// ##################################################################
// ##################################################################


void VectorOps_Cyl::Gradient(
        const cell *c, 
        const int sv, 
        const int var, 
        class GridBaseClass *grid,
        pion_flt *grad
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cart::Grad: Some neighbour cells don't exist",i);
#endif //TESTING
  
  cell *cn,*cp;
  double dR=grid->DX();
  double dZ=dR;
  
  switch (sv) {
   case 0: // Use vector c->P
      cn = grid->NextPt(c,ZNcyl); cp=grid->NextPt(c,ZPcyl);
      grad[0] = (cp->P[var] - cn->P[var])/(2.*dZ);
    if(VOnd>1) {
      cn = grid->NextPt(c,RNcyl); cp=grid->NextPt(c,RPcyl);
      grad[1] = (cp->P[var] - cn->P[var])/(R_com(cp,dR) - R_com(cn,dR));
    }
    if(VOnd>2) {
      cn = grid->NextPt(c,TNcyl); cp=grid->NextPt(c,TPcyl);
      grad[2] = (cp->P[var] - cn->P[var])/(CI.get_dpos(cp,Tcyl) - CI.get_dpos(cn,Tcyl))/R_com(c,dR);
    }
    break;
   case 1: // Use Vector c-Ph
      cn = grid->NextPt(c,ZNcyl); cp=grid->NextPt(c,ZPcyl);
      grad[0] = (cp->Ph[var] - cn->Ph[var])/(2.*dZ);
    if(VOnd>1) {
      cn = grid->NextPt(c,RNcyl); cp=grid->NextPt(c,RPcyl);
      grad[1] = (cp->Ph[var] - cn->Ph[var])/(R_com(cp,dR) - R_com(cn,dR));
    }
    if(VOnd>2) {
      cn = grid->NextPt(c,TNcyl); cp=grid->NextPt(c,TPcyl);
      grad[2] = (cp->Ph[var] - cn->Ph[var])/(CI.get_dpos(cp,Tcyl) - CI.get_dpos(cn,Tcyl))/R_com(c,dR);
    }
    break;
   default:
    rep.error("Don't know what state vector to use for calculating gradiend.",sv);
  }
  return;
}


// ##################################################################
// ##################################################################

double VectorOps_Cyl::Divergence(
        cell *c, 
        const int sv, 
        const int *var,
        class GridBaseClass *grid
        )
{
  // get divergence of vector quantity.


  double dT=0.;
  double dR=grid->DX();
  double dZ=dR;
  double dx[3];
  double divv=0.;
  cell *cn,*cp;
  cell *ngb[6];
  enum direction ndir;
  enum direction pdir;
  

  for (int v=0;v<VOnd;v++) {
    ndir = static_cast<direction>(2*v);
    pdir = static_cast<direction>(2*v+1);
    ngb[ndir] = (grid->NextPt(c,ndir)) ? grid->NextPt(c,ndir) : c;
    ngb[pdir] = (grid->NextPt(c,pdir)) ? grid->NextPt(c,pdir) : c;
    if (ngb[ndir]==c || ngb[pdir]==c) dx[v]=dZ;
    else                              dx[v]=2.0*dZ;
  }

  switch (sv) {
    case 0: // Use vector c->P
    // d(V_z)/dz
    cn = ngb[ZNcyl];
    cp = ngb[ZPcyl];
    divv  = (cp->P[var[0]] - cn->P[var[0]])/dx[0];
    if (VOnd>1) { // d(R* V_R)/(R*dR)
      cn = ngb[RNcyl];
      cp = ngb[RPcyl];
      double rn = R_com(cn,dR);
      double rp = R_com(cp,dR);
      divv  += 2.0*(rp*cp->P[var[1]] - rn*cn->P[var[1]])/(rp*rp-rn*rn);
    }
    if (VOnd>2) { // d(V_theta)/(R*dtheta)
      cn = ngb[TNcyl];
      cp = ngb[TPcyl];
      dT = (CI.get_dpos(cp,Tcyl) - CI.get_dpos(cn,Tcyl))*R_com(c,dR);
      divv  += (cp->P[var[2]] - cn->P[var[2]])/dT;
    }
    break;

    case 1: // Use Vector c-Ph
    // d(V_z)/dz
    cn = ngb[ZNcyl];
    cp = ngb[ZPcyl];
    divv  = (cp->Ph[var[0]] - cn->Ph[var[0]])/dx[0];
    if (VOnd>1) { // d(R* V_R)/(R*dR)
      cn = ngb[RNcyl];
      cp = ngb[RPcyl];
      double rn = R_com(cn,dR);
      double rp = R_com(cp,dR);
      divv  += 2.0*(rp*cp->Ph[var[1]] - rn*cn->Ph[var[1]])/(rp*rp-rn*rn);
    }
    if (VOnd>2) { // d(V_theta)/(R*dtheta)
      cn = ngb[TNcyl];
      cp = ngb[TPcyl];
      dT = (CI.get_dpos(cp,Tcyl) - CI.get_dpos(cn,Tcyl))*R_com(c,dR);
      divv  += (cp->Ph[var[2]] - cn->Ph[var[2]])/dT;
    }
    break;
   default:
    rep.error("CYL: state vector for calculating divergence.",sv);
  }
  return(divv);
} // Div


// ##################################################################
// ##################################################################

void VectorOps_Cyl::Curl(
        const cell *c, 
        const int vec, 
        const int *var, 
        class GridBaseClass *grid,
        pion_flt *ans
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cyl::Curl: Some neighbour cells don't exist",i);
#endif //TESTING
  if (!c->isgd) rep.error("Not Grid Cell! can't calculate curl. id follows",c->id);
  cout <<"Cyl_Curl() is not tested!!! make sure it works!!!\n";
  // variables z,R,theta = z,r,t
  int vz=var[0], vr=var[1], vt=var[2];
  pion_flt *vzp=0, *vzn=0, *vrp=0, *vrn=0, *vtp=0, *vtn=0;
  double rp=0., rn=0.;
  double dR=grid->DX();
  double dZ=dR;
  
  switch (vec) {
   case 0: // c->P
    cout <<"using P.\n";
    vzp = grid->NextPt(c,ZPcyl)->P; vzn = grid->NextPt(c,ZNcyl)->P;
    vrp = grid->NextPt(c,RPcyl)->P; vrn = grid->NextPt(c,RNcyl)->P;
    vtp = grid->NextPt(c,TPcyl)->P; vtn = grid->NextPt(c,TNcyl)->P;
    break;
   case 1: // c->Ph
    cout <<"using Ph = ["<<c->Ph[vz]<<", "<<c->Ph[vr]<<", "<<c->Ph[vt]<<" ]\n";
    vzp = grid->NextPt(c,ZPcyl)->Ph; vzn = grid->NextPt(c,ZNcyl)->Ph;
    vrp = grid->NextPt(c,RPcyl)->Ph; vrn = grid->NextPt(c,RNcyl)->Ph;
    vtp = grid->NextPt(c,TPcyl)->Ph; vtn = grid->NextPt(c,TNcyl)->Ph;
    break;
   case 2: // c->dU
    cout <<"using dU.\n";
    vzp = grid->NextPt(c,ZPcyl)->dU; vzn = grid->NextPt(c,ZNcyl)->dU;
    vrp = grid->NextPt(c,RPcyl)->dU; vrn = grid->NextPt(c,RNcyl)->dU;
    vtp = grid->NextPt(c,TPcyl)->dU; vtn = grid->NextPt(c,TNcyl)->dU;
    break;
   default:
    rep.error("Which vector to calculate on?  (VecCurl), don't know!",vec);
  }

  // First dz derivatives
  rp = R_com(grid->NextPt(c,ZPcyl),dR);
  rn = R_com(grid->NextPt(c,ZNcyl),dR);
  ans[Zcyl] =   0.0;
  ans[Rcyl] = -(rp*vzp[vt] - rn*vzn[vt]) /(2.*dZ);
  ans[Tcyl] =  (vzp[vr] - vzn[vr]) /(2.*dZ);
  if (VOnd>1) { // now R derivatives, if present
    rp = R_com(grid->NextPt(c,RPcyl),dR);
    rn = R_com(grid->NextPt(c,RNcyl),dR);
    ans[Zcyl] +=  (rp*vrp[vt] - rn*vrn[vt])/(rp-rn);
    ans[Rcyl] +=   0.0;
    ans[Tcyl] += -(vrp[vz] - vrn[vz])/(rp-rn);
  }
  if (VOnd>2) { // now Theta derivatives if present.
    rp = CI.get_dpos(grid->NextPt(c,TPcyl),Tcyl); rn = CI.get_dpos(grid->NextPt(c,TNcyl),Tcyl);
    ans[Zcyl] += -(vtp[vr] - vtn[vr])/(rp-rn);
    ans[Rcyl] +=  (vtp[vz] - vtn[vz])/(rp-rn);
    ans[Tcyl] +=  0.0;
  }
  ans[Zcyl] /= R_com(c,dR);
  ans[Rcyl] /= R_com(c,dR);
  // Theta component needs no further modification.

  return;
}// VecCurl


// ##################################################################
// ##################################################################


int VectorOps_Cyl::SetEdgeState(
        const cell *c,       ///< Current Cell.
        const direction dir, ///< Add or subtract the slope depending on direction.
        const int nv,        ///< length of state vectors.
        const pion_flt *dpdx,  ///< Slope vector.
        pion_flt *edge,        ///< vector for edge state. 
        const int OA,        ///< Order of spatial Accuracy.
        class GridBaseClass *grid
        )
{
  if (OA==OA1) { // 1st order, constant data.
    for (int v=0;v<nv;v++) edge[v] = c->Ph[v];
  }
  
  else if (OA==OA2) {
    double del=0.;
    double dR=grid->DX();
    double dZ=dR;
    switch (dir) {
     case ZPcyl: del = dZ/2.; break;
     case ZNcyl: del =-dZ/2.; break;
     case RPcyl:
      del = CI.get_dpos(c,Rcyl)+dR/2. - R_com(c,dR);
      break;
     case RNcyl:
      del = CI.get_dpos(c,Rcyl)-dR/2.  - R_com(c,dR);
      break;
     case TPcyl:
       del = R_com(c,dR)*(CI.get_dpos(grid->NextPt(c,TPcyl),Tcyl) - CI.get_dpos(c,Tcyl));
      break;
     case TNcyl:
       del = R_com(c,dR)*(CI.get_dpos(grid->NextPt(c,TNcyl),Tcyl) - CI.get_dpos(c,Tcyl));
      break;
     default:
      rep.error("Bad direction in SetEdgeState",dir);
    } // setting del based on direction.
    
    for (int v=0;v<nv;v++) edge[v] = c->Ph[v] + dpdx[v]*del;
  } // OA2
  else rep.error("Bad order of accuracy in SetEdgeState -- only know 1st and 2nd order space",OA);
  
  return(0);
} // SetEdgeState
  

// ##################################################################
// ##################################################################


int VectorOps_Cyl::SetSlope(
        const cell *c, ///< Current Cell.
        const axes d,  ///< Which direction to calculate slope in.
        const int nv,  ///< length of state vectors.
        pion_flt *dpdx,  ///< Slope vector to be written to.
        const int  OA, ///< Order of spatial Accuracy.
        class GridBaseClass *grid
        )
{
  double dR = grid->DX();
  double dZ = dR;
  if (OA==OA1){ // first order accurate so zero slope.
    for(int v=0;v<nv;v++) dpdx[v]=0.;
  } // 1st order.
  
  else if (OA==OA2) { // second order spatial accuracy.
    pion_flt slpn[nv], slpp[nv];
    cell *cn=0,*cp=0;
    enum direction dp=NO,dn=NO;
    switch (d) {
     case Zcyl: dp=ZPcyl; dn=ZNcyl; break;
     case Rcyl: dp=RPcyl; dn=RNcyl; break;
     case Tcyl: dp=TPcyl; dn=TNcyl; break;
     default: rep.error("Bad direction in SetSlope",d);
    }
    cp = grid->NextPt(c,dp); cn=grid->NextPt(c,dn);
#ifdef TESTING
    if (cp==0 || cn==0) rep.error("No left or right cell in SetSlope",cp);
#endif //TESTING
    switch (d) {
     case Zcyl: 
      for (int v=0;v<nv;v++) {
	slpn[v] = (c->Ph[v] -cn->Ph[v])/dZ;
	slpp[v] = (cp->Ph[v]- c->Ph[v])/dZ;
	dpdx[v] = AvgFalle(slpn[v],slpp[v]);
      }
#ifdef TESTING
      rep.printVec("Z slpn",slpn,nv);
      rep.printVec("Z slpp",slpp,nv);
      rep.printVec("Z dpdx",dpdx,nv);
      cout <<"R_com(c) = "<<R_com(c,dR)<<" vodz="<<dR<<"\n";
#endif //TESTING
      break;
     case Rcyl:
      for (int v=0;v<nv;v++) {
	slpn[v] = (c->Ph[v] -cn->Ph[v])/ (R_com(c,dR) - R_com(cn,dR));
	slpp[v] = (cp->Ph[v]- c->Ph[v])/ (R_com(cp,dR)- R_com(c,dR) );
	dpdx[v] = AvgFalle(slpn[v],slpp[v]);
      }
#ifdef TESTING
      rep.printVec("R slpn",slpn,nv);
      rep.printVec("R slpp",slpp,nv);
      rep.printVec("R dpdx",dpdx,nv);
      cout <<"R_com(c) = "<<R_com(c,dR)<<"\n";
#endif //TESTING
      break;
     case Tcyl: // Need extra scale factor in denominator to get R*dtheta
      for (int v=0;v<nv;v++) {
	slpn[v] = (c->Ph[v] -cn->Ph[v])/ R_com(c,dR)/(CI.get_dpos(c,Tcyl)  - CI.get_dpos(cn,Tcyl));
	slpp[v] = (cp->Ph[v]- c->Ph[v])/ R_com(c,dR)/(CI.get_dpos(cp,Tcyl) - CI.get_dpos(c,Tcyl));
	dpdx[v] = AvgFalle(slpn[v],slpp[v]);
      }
      break;
     default:
      rep.error("Bad axis in SetSlope",d);
    } // calculate slope in direction d.
#ifdef FIX_TRACER_SLOPES_TO_DENSITY
    if (SimPM.ntracer>0) {
      for (int v=SimPM.ftr; v<nv; v++) {
	//if (!pconst.equalD(dpdx[v],0.0)) 
        //  cout <<"setting dpdx["<<v<<"] from "<<dpdx[v]<<" to zero.\n";
	dpdx[v] =0.0;
      }
    }
#endif //FIX_TRACER_SLOPES_TO_DENSITY
  } // 2nd order accurate
  else {
    cerr <<"Error: Only know how to do 1st or 2nd order slope calculation.\n";
    exit(1);
  }
  return(0);
} // SetSlope


// ##################################################################
// ##################################################################


int VectorOps_Cyl::DivStateVectorComponent(
        const cell *c,    ///< current cell.
        class GridBaseClass *grid, 
        const axes d,     ///< current coordinate axis we are looking along.
        const int nv,     ///< length of state vectors.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        pion_flt *dudt      ///< Vector to assign divergence component to.
        )
{
  //enum direction dp;
  //switch (d) {
  // case Zcyl: dp=ZPcyl; break;
  // case Rcyl: dp=RPcyl; break;
  // case Tcyl: dp=TPcyl; break;
  // default: rep.error("Bad direction in DivStateVectorComponent",d);
  //}

  /** \section Sign
   * Note that this function returns the negative of the i-th component of 
   * the divergence.  This is b/c it is used in the finite volume time update
   * where what is needed is the negative of div(F).
   * */

  double dZ = grid->DX();
  double dR = dZ;
  
  if (d==Zcyl) for (int v=0;v<nv;v++) dudt[v] = (fn[v]-fp[v])/dZ;
  
  else if (d==Rcyl) {
    double rp = CI.get_dpos(c,Rcyl)+dR/2.; double rn=rp-dR;
    for (int v=0;v<nv;v++)
      dudt[v] = 2.0*(rn*fn[v]-rp*fp[v])/(rp*rp-rn*rn);
  }
  else if (d==Tcyl) {
    double dth = R_com(c,dR)*(CI.get_dpos(grid->NextPt(c,TPcyl),Tcyl) -
			   CI.get_dpos(c,Tcyl));
    for (int v=0;v<nv;v++) dudt[v] = (fn[v]-fp[v])/dth;
  }
  else {
    rep.error("Bad axis in DivStateVectorComponent",d);
  }
  return(0);
} // DivStateVectorComponent


// ##################################################################
// ##################################################################


