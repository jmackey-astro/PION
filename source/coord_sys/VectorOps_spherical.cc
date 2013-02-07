///
/// \file VectorOps_spherical.cc
///
/// \author Jonathan Mackey
///
/// Class for spherical coordinates.
///
/// Created 2010.10.01
///
/// Modifications:
/// - 2010.10.04 JM: Fixed bugs, so now it works.
/// - 2010.11.03 JM: Changed finite diff. for divergence.
/// - 2010.12.04 JM: Added constructor with only one argument.  Also
///   a set_dx() function.
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
///


/// ***************************************
/// ******** SPHERICAL COORDINATES ********
/// ***************************************

#include "VectorOps_spherical.h"
using namespace std;


// ##################################################################
// ##################################################################

VectorOps_Sph::VectorOps_Sph(int n, double del)
  : VectorOps_Cart(n,del), VectorOps_Cyl(n,del)
{
#ifdef TESTING
  cout <<"Setting up 1D spherical coordinates with ndim="<<VOnd;
  cout <<" and dR="<<VOdR<<"\n";
#endif
  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",
                          VOnd);
  
  VOdV = 4.0*M_PI*VOdR;
  VOdA = 4.0*M_PI;
  return;
}


// ##################################################################
// ##################################################################

VectorOps_Sph::VectorOps_Sph(int n)
  : VectorOps_Cart(n), VectorOps_Cyl(n)
{
  cout <<"Setting up 1D spherical coordinates with ndim="<<VOnd<<"\n";
  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);
  have_set_dx=false;
  return;
}


// ##################################################################
// ##################################################################

void VectorOps_Sph::set_dx(const double x)
{
  //
  // This will set VOdx, VOdA, VOdV for Cartesian and Cylindrical.
  //
  VectorOps_Cyl::set_dx(x);

  //
  // Now modify Volume and Area for Spherical coords.
  //
  VOdV = 4.0*M_PI*VOdR;
  VOdA = 4.0*M_PI;
  return;
}


// ##################################################################
// ##################################################################

VectorOps_Sph::~VectorOps_Sph()
{}


// ##################################################################
// ##################################################################

double VectorOps_Sph::CellVolume(const cell *c)
{
  ///
  /// The volume of a cell is 
  /// \f$\delta V =4\pi R_i^2 \delta R \f$, where \f$R_i\f$ is the
  /// cell centre in the radial direction.
  ///
  double temp = CI.get_dpos(c,Rsph);
  return VOdV*temp*temp;
}


// ##################################################################
// ##################################################################

double VectorOps_Sph::CellInterface(
        const cell *c,
        const direction dir
        )
{
  double temp=0.0;

  if (VOnd==1) {
    switch (dir) {
    case RNsph:
      temp =  CI.get_dpos(c,Rsph) - 0.5*VOdR;
      return 4.0*M_PI*temp*temp;
      break;
    case RPsph:
      temp =  CI.get_dpos(c,Rsph) + 0.5*VOdR;
      return 4.0*M_PI*temp*temp;
      break;
    default:
      rep.error("Bad direction in 1D spherical coords!",dir);
      break;
    }
  }

  rep.error("VectorOps_Sph::CellInterface","fix me");
  return -1.0;
}
    

// ##################################################################
// ##################################################################

double VectorOps_Sph::maxGradAbs(
        const cell *c,
        const int sv,
        const int var
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Sph::maxGradAbs: Some neighbour cells don't exist",i);
#endif //TESTING

  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);

  //
  // 1D gradient in radial direction is just df/dr
  //
  double grad=0, temp=0;
  cell *cn;
  switch (sv) {
  case 0: // Use vector c->P
    cn = grid->NextPt(c,RPsph);
    temp = fabs( cn->P[var] - c->P[var])/(R_com(cn)-R_com(c));
    if(temp>grad) grad = temp;
    cn = grid->NextPt(c,RNsph);
    temp = fabs( cn->P[var] - c->P[var])/(R_com(c)-R_com(cn));
    if(temp>grad) grad = temp;
    break;
  case 1: // Use Vector c-Ph
    cn = grid->NextPt(c,RPsph);
    temp = fabs( cn->Ph[var] - c->Ph[var])/(R_com(cn)-R_com(c));
    if(temp>grad) grad = temp;
    cn = grid->NextPt(c,RNsph);
    temp = fabs( cn->Ph[var] - c->Ph[var])/(R_com(c)-R_com(cn));
    if(temp>grad) grad = temp;
    break;
   default:
    rep.error("Don't know what state vector to use for calculating gradient.",sv);
  }
  return(grad);
} // maxGradAbs


// ##################################################################
// ##################################################################


void VectorOps_Sph::Grad(
        const cell *c,
        const int sv,
        const int var,
        double *grad
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Cart::maxGradAbs: Some neighbour cells don't exist",i);
#endif //TESTING
  
  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);
  //
  // 1D gradient in radial direction is just df/dr.  The sign is
  // automatically correct for r<0, so we don't need to make any
  // modifications.
  //
  cell *cn,*cp;
  double rn=0.0, rp=0.0;
  switch (sv) {
  case 0: // Use vector c->P
    cn = grid->NextPt(c,RNsph); cp=grid->NextPt(c,RPsph);
    rn = R_com(cn); rp = R_com(cp);
    grad[0] = (cp->P[var] - cn->P[var])/(rp-rn);
    break;
  case 1: // Use Vector c-Ph
    cn = grid->NextPt(c,RNsph); cp=grid->NextPt(c,RPsph);
    rn = R_com(cn); rp = R_com(cp);
    grad[0] = (cp->Ph[var] - cn->Ph[var])/(rp-rn);
    break;
  default:
    rep.error("Don't know what state vector to use for calculating gradient.",sv);
  }
  return;
}


// ##################################################################
// ##################################################################

// get divergence of vector quantity.
double VectorOps_Sph::Div(
        const cell *c,
        const int sv,
        const int *var
        )
{

#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Sph::Div: Some neighbour cells don't exist",i);
#endif //TESTING

  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);
  //
  // 1D divergence in radial direction is (1/r^2)d/dr(r^2*V_r)
  //
  double divv=0.0, rn=0.0, rp=0.0;
  cell *cn,*cp;
  
  switch (sv) {
  case 0: // Use vector c->P
    // r^{-2}d(r^2 V_r)/dr or (2/r)V_r +d(V_r)/dr
    cn = grid->NextPt(c,RNsph); cp=grid->NextPt(c,RPsph);
    rn = R_com(cn); rp = R_com(cp);
    //divv = (rp*rp*cp->P[var[0]] - rn*rn*cn->P[var[0]])*3.0/(pow(rp,3.0)-pow(rn,3.0));
    divv = 2.0*c->P[var[0]]/R_com(c) +(cp->P[var[0]]-cn->P[var[0]])/(rp-rn);
    break;
  case 1: // Use Vector c-Ph
    // r^{-2}d(r^2 V_r)/dr or (2/r)V_r+dV_r/dr
    cn = grid->NextPt(c,RNsph); cp=grid->NextPt(c,RPsph);
    rn = R_com(cn); rp = R_com(cp);
    //divv = (rp*rp*cp->P[var[0]] - rn*rn*cn->P[var[0]])*3.0/(pow(rp,3.0)-pow(rn,3.0));
    divv = 2.0*c->P[var[0]]/R_com(c) +(cp->P[var[0]]-cn->P[var[0]])/(rp-rn);
    break;
  default:
    rep.error("Don't know what state vector to use for calculating divergence.",sv);
  }
  return(divv);
} // Div


// ##################################################################
// ##################################################################


void VectorOps_Sph::Curl(
        const cell *c,
        const int vec,
        const int *var,
        double *ans
        )
{
#ifdef TESTING
  for (int i=0;i<2*VOnd; i++)
    if (!grid->NextPt(c,static_cast<direction>(i)))
      rep.error("VectorOps_Sph::Curl: Some neighbour cells don't exist",i);
#endif //TESTING
  if (!c->isgd) rep.error("Not Grid Cell! can't calculate curl. id follows",c->id);
  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);

  //
  // 1D curl is always zero
  //
  for (int v=0;v<MAX_DIM;v++) ans[v]=0.0;
  return;
} // VecCurl


// ##################################################################
// ##################################################################


int VectorOps_Sph::SetEdgeState(
        const cell *c,       ///< Current Cell.
        const direction dir, ///< Add or subtract the slope depending on direction.
        const int nv,        ///< length of state vectors.
        const double *dpdx,  ///< Slope vector.
        double *edge,        ///< vector for edge state. 
        const int OA         ///< Order of spatial Accuracy.
        )
{
  
  if (VOnd!=1) rep.error("Spherical coordinates only work in 1D!",VOnd);

  //
  // 1st order, constant data.
  //
  if (OA==OA1) {
    for (int v=0;v<nv;v++) edge[v] = c->Ph[v];
  }
  
  //
  // 2nd order, linear data, pivot point is centre of gravity in R.
  //
  else if (OA==OA2) {
    double del=0.;
    switch (dir) {
    case RPsph:
      del = CI.get_dpos(c,Rsph)+0.5*VOdR - R_com(c);
      break;
    case RNsph:
      del = CI.get_dpos(c,Rsph)-0.5*VOdR - R_com(c);
      break;
     default:
      rep.error("Bad direction in SetEdgeState",dir);
    } // setting del based on direction
    
    for (int v=0;v<nv;v++) edge[v] = c->Ph[v] + dpdx[v]*del;
  } // OA2

  else rep.error("Bad order of accuracy in SetEdgeState \
                  -- only know 1st and 2nd order space",OA);
  
  return 0;
} // SetEdgeState


// ##################################################################
// ##################################################################


int VectorOps_Sph::SetSlope(
        const cell *c, ///< Current Cell.
        const axes d,  ///< Which direction to calculate slope in.
        const int nv,  ///< length of state vectors.
        double *dpdx,  ///< Slope vector to be written to.
        const int  OA  ///< Order of spatial Accuracy.
        )
{
  //
  // first order accurate so zero slope.
  //
  if (OA==OA1){
    for(int v=0;v<nv;v++) dpdx[v]=0.;
  } // 1st order.
  
  //
  // second order spatial accuracy.
  // 
  else if (OA==OA2) {
    double slpn[nv], slpp[nv];
    cell *cn=0,*cp=0;
    enum direction dp=NO,dn=NO;

    //
    // Set directions and neigbour cells
    //
    switch (d) {
     case Rsph: dp=RPsph; dn=RNsph; break;
     default: rep.error("Bad direction in SetSlope",d);
    }
    cp = grid->NextPt(c,dp); cn=grid->NextPt(c,dn);
#ifdef TESTING
    if (cp==0 || cn==0) rep.error("No left or right cell in SetSlope",cp);
#endif //TESTING

    //
    // Calculate slope based on axes we are looking along.
    //
    switch (d) {
    case Rsph:
      for (int v=0;v<nv;v++) {
	slpn[v] = (c->Ph[v] -cn->Ph[v])/ (R_com(c) - R_com(cn));
	slpp[v] = (cp->Ph[v]- c->Ph[v])/ (R_com(cp)- R_com(c) );
	dpdx[v] = AvgFalle(slpn[v],slpp[v]);
      }
      break;
    default:
      rep.error("Bad axis in SetSlope",d);
    } // calculate slope in direction d.

#ifdef FIX_TRACER_SLOPES_TO_DENSITY
    if (SimPM.ntracer>0) {
      for (int v=SimPM.ftr; v<nv; v++) {
	//if (!GS.equalD(dpdx[v],0.0)) 
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

  return 0;
} // SetSlope


// ##################################################################
// ##################################################################


int VectorOps_Sph::DivStateVectorComponent(
        const cell *c,    ///< current cell.
        const axes d,     ///< current coordinate axis we are looking along.
        const int nv,     ///< length of state vectors.
        const double *fn, ///< Negative direction flux.
        const double *fp, ///< Positive direction flux.
        double *dudt      ///< Vector to assign divergence component to.
        )
{
  ///
  /// \section Sign
  /// Note that this function returns the negative of the i-th component of 
  /// the divergence.  This is b/c it is used in the finite volume time update
  /// where what is needed is the negative of div(F).
  ///
  /// The finite difference is from Boss & Myhill (1992, ApJS, 83, 311), eq. 15.
  ///
  /// For boundary cells with r<0, the cubic nature of the denominator
  /// ensures the signs cancel and the answer is correct.
  ///

  if      (d==Rsph) {
    //    cout <<"updating flux divergence R\n";
    double rc = CI.get_dpos(c,Rsph);
    double rp = rc+0.5*VOdR; double rn=rp-VOdR;
    rc = (pow(rp,3.0) -pow(rn,3.0))/3.0;
    for (int v=0;v<nv;v++) dudt[v] = (rn*rn*fn[v]-rp*rp*fp[v])/rc;
  }
  else {
    rep.error("Bad axis in DivStateVectorComponent",d);
  }

  return 0;
} // DivStateVectorComponent


// ##################################################################
// ##################################################################

