///
/// \file solver_eqn_hydro_adi.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 written, based on old solver.cc classes.
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via either Lax-Friedrichs
/// or Riemann solver (linear and/or exact).  Adds viscosity if asked for, and tracks flux
/// of N passive tracers.
///
/// - 2009-12-18 JM: Added Axisymmetric Class (cyl_FV_solver_Hydro_Euler)
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.10.01 JM: Added spherical coordinate system.
/// - 2010.11.03 JM: Fixed source terms for spherical coordinates
///   (although results don't change much).
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.22 JM: Added new Riemann solver classes for Hydro.
/// - 2010.12.23 JM: Removed references to riemann_base class.
///    added extra variable to inviscid_flux() function.
///    Moved UtoP() etc. from solver to flux-solver.
/// - 2010.12.28 JM: removed some debugging comments.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2013.08.19 JM: tested a bunch of approximations, but nothing
///    was an improvement so I left it the way it was.

#include "solver_eqn_hydro_adi.h"
using namespace std;

// *********************************
// ***** FV SOLVER HYDRO EULER *****
// *********************************


// ##################################################################
// ##################################################################

FV_solver_Hydro_Euler::FV_solver_Hydro_Euler(
        const int nv, ///< number of variables in state vector.
        const int nd, ///< number of space dimensions in grid.
        const double cflno,   ///< CFL number
        const double cellsize,    ///< dx, cell size.
        const double gam,     ///< gas eos gamma.
        double *state,     ///< State vector of mean values for simulation.
        const double avcoeff, ///< Artificial Viscosity Parameter etav.
        const int ntr         ///< Number of tracer variables.
        )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv), riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize)
{
#ifdef TESTING
  cout <<"FV_solver_Hydro_Euler::FV_solver_Hydro_Euler() constructor.\n";
  cout <<"FV_solver_Hydro_Euler::FV_solver_Hydro_Euler() gamma = "<<eq_gamma<<"\n";
#endif
  return;
}


// ##################################################################
// ##################################################################

FV_solver_Hydro_Euler::~FV_solver_Hydro_Euler()
{
#ifdef TESTING
  cout <<"FV_solver_Hydro_Euler::~FV_solver_Hydro_Euler() destructor.\n";
#endif
  return;
}



// ##################################################################
// ##################################################################


///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_Hydro_Euler::dU_Cell(cell *c,          // Current cell.
				   const axes d,     // Which axis we are looking along.
				   const double *fn, // Negative direction flux.
				   const double *fp, // Positive direction flux.
				   const double *,   // slope vector for cell c.
				   const int,        // spatial order of accuracy.
				   const double,     // cell length dx.
				   const double      // cell TimeStep, dt.
				   )
{
  double u1[eq_nvar];
  //
  // This calculates -dF/dx
  //
  //if (d!=eq_dir) rep.error("direction problem!!!!!!!!",d);
  int err = DivStateVectorComponent(c,d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  return(err);
}


// ##################################################################
// ##################################################################

///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_Hydro_Euler::CellAdvanceTime(class cell *c,
					   const double *Pin, // Initial State Vector.
					   double *dU, // Update vector dU
					   double *Pf, // Final state vector (can be same as initial vec.).
					   double *dE, // TESTING Tracks change of energy for negative pressure correction.
					   const double, // gas EOS gamma.
					   const double  // Cell timestep dt.
					   )
{
  double u1[eq_nvar];
#ifdef TESTING
  double u2[eq_nvar];
#endif //TESTING

  //
  // First convert from Primitive to Conserved Variables
  //
  PtoU(Pin,u1, eq_gamma);

  //
  // Now add dU[] to U[], and change back to primitive variables.
  // This can give negative pressures, so check for that and fix it if needed.
  //
  for (int v=0;v<eq_nvar;v++) {
    u1[v] += dU[v];   // Update conserved variables
    dU[v] = 0.;       // Reset the dU array for the next timestep.
  }
  if(UtoP(u1,Pf, eq_gamma)!=0) {
    cout<<"(FV_solver_Hydro_Euler::CellAdvanceTime) UtoP complained \
           (maybe about negative pressure...) fixing\n";
#ifdef TESTING
    //grid->PrintCell(dp.c);
    //grid->PrintCell(dp.c->ngb[XP]);
    //rep.printVec("pin",Pin,SimPM.nvar);
    //PtoU(Pin,u1, eq_gamma);
    //rep.printVec("Uin",u1,SimPM.nvar);
    //rep.printVec("dU ",dU, SimPM.nvar);
    PtoU(Pf, u2, eq_gamma);
    *dE += (u2[ERG]-u1[ERG]);
    UtoP(u2,Pf, eq_gamma);
#endif //TESTING
  }

//#ifdef TEST_MP7
//  double Ta = Pf[PG]*GS.m_p()/((1.0+Pf[SimPM.ftr])*GS.kB()*Pf[RO]);
//  double Tt = MP->Temperature(Pf,SimPM.gamma);
//  if (fabs(Ta-Tt)/Tt >1.0e-2) {
//    cout <<"T_now="<<Ta<<", but it should be "<<Tt<<"\n";
//  }
//  MP->Set_Temp(Pf,Tt,SimPM.gamma);
//#endif // TEST_MP7

  //  for (int v=0;v<eq_nvar;v++) dU[v] = 0.; // Reset the dU array for the next timestep.

  return 0;
}


// ##################################################################
// ##################################################################

///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_Hydro_Euler::CellTimeStep(const cell *c, ///< pointer to cell
					   const double, ///< gas EOS gamma.
					   const double  ///< Cell size dx.
					   )
{
  /** \section Algorithm
   * First Get the maximum fluid velocity in each of the three directions.
   * Add the sound speed to the abs. value of this to get the max signal 
   * propagation speed.  Then dt=dx/max.vel.
   * Next test if the pressure gradient can induce a larger speed in the
   * timestep, and if it can, divide the timestep by the appropriate factor.
   * Finally multiply by the CFl no. and return.
   * */

  //
  // Get Max velocity along a grid direction.
  //
  double temp = fabs(c->P[eqVX]);
  if (FV_gndim>1) temp = max(temp,fabs(c->P[eqVY]));
  if (FV_gndim>2) temp = max(temp,fabs(c->P[eqVZ]));
  
  //
  // Add the sound speed to this, and it is the max wavespeed.
  //
  temp += chydro(c->P,eq_gamma);
  FV_dt = FV_dx/temp;

  
  //
  // Check the gradient of pressure with neighbouring cells, since this can 
  // dramatically shorten the timestep.  (CAN'T REMEMBER WHY!!!)
  //
  // int pg = static_cast<int>(eqPG);
  // if (c->isgd) {
  //   double grad = maxGradAbs(c,0,pg)/c->P[RO];
  //   if( (temp = grad*FV_dt/temp) >1.) {
  //     FV_dt /= temp;
  //     //      if (temp>3) cout <<"gradient adjusting timestep\n";
  //   }
  // }
  // else if (grid->NextPt(c,XP)) {
  //   //  cout <<"grid="<<grid<<"\n";
  //   double grad = fabs(c->P[PG]-grid->NextPt(c,XP)->P[PG])/FV_dx/c->P[RO];
  //   if( (temp = grad*FV_dt/temp) >1.) {
  //     FV_dt /= temp;
  //   }
  // }
  // else rep.error("No neighbour to gradient test",c);    
  
  //
  // Now scale the max. allowed timestep by the CFL number we are using (<1).
  //
  FV_dt *= FV_cfl;
  return FV_dt;
}


// ##################################################################
// ##################################################################

////////////////////////////////////////////////////////////////////////
/// CYLINDRICAL (AXISYMMETRIC) COORDINATES
////////////////////////////////////////////////////////////////////////


// ##################################################################
// ##################################################################

cyl_FV_solver_Hydro_Euler::cyl_FV_solver_Hydro_Euler(
        const int nv, ///< number of variables in state vector.
        const int nd, ///< number of space dimensions in grid.
        const double cflno, ///< CFL number
        const double cellsize, ///< dx, cell size.
        const double gam,     ///< gas eos gamma.
        double *state, ///< State vector of mean values for simulation.
        const double avcoeff, ///< Artificial Viscosity Parameter etav.
        const int ntr         ///< Number of tracer variables.
        )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv), riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_Hydro_Euler(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
#ifdef TESTING
  cout <<"cyl_FV_solver_Hydro_Euler::cyl_FV_solver_Hydro_Euler() constructor.\n";
#endif
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for \
                        2d axial symmetry so far.  Sort it out!",nd);
  return;
}


// ##################################################################
// ##################################################################

cyl_FV_solver_Hydro_Euler::~cyl_FV_solver_Hydro_Euler()
{
#ifdef TESTING
  cout <<"cyl_FV_solver_Hydro_Euler DESTRUCTOR.\n";
#endif
}


// ##################################################################
// ##################################################################

int cyl_FV_solver_Hydro_Euler::dU_Cell(cell *c, ///< Current cell.
				       const axes d, ///< Which axis we are looking along.
				       const double *fn, ///< Negative direction flux.
				       const double *fp, ///< Positive direction flux.
				       const double *dpdx, ///< slope vector for cell c.
				       const int OA,      ///< spatial order of accuracy.
				       const double, ///< cell length dx.
				       const double  ///< cell TimeStep, dt.
				       )
{
  double u1[eq_nvar];
  //
  // This calculates -dF/dx
  //
  //if (d!=eq_dir) rep.error("direction problem!!!!!!!!",d);
  int err = DivStateVectorComponent(c,d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  //
  // Add source term for the radial direction.
  //
  if (d==Rcyl) {
    switch (OA) {
     case OA1:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]/CI.get_dpos(c,Rcyl));
      //u1[0] = CI.get_dpos(c,Rsph)-0.5*VOdR;
      //u1[1] = u1[0]+VOdR;
      //u1[2] = 2.0*(u1[1]*u1[3] - u1[0]*u1[2])/(u1[1]*u1[1]-u1[0]*u1[0]);
      //cout <<"O1: rel.diff. Src Term: "<< (c->Ph[eqPG]/CI.get_dpos(c,Rcyl) -u1[2])/u1[2];
      //cout <<"  old = "<< FV_dt*c->Ph[eqPG]/CI.get_dpos(c,Rcyl)<<",  new = "<<FV_dt*u1[2]<<",  p="<<c->Ph[eqPG]<<"\n";
      //c->dU[eqMMX] += FV_dt*u1[2];
      break;
     case OA2:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG] + dpdx[eqPG]*(CI.get_dpos(c,Rcyl)-R_com(c)))/CI.get_dpos(c,Rcyl);
      //u1[0] = CI.get_dpos(c,Rsph)-0.5*VOdR;              // r-
      //u1[1] = u1[0]+VOdR;                                // r+
      //u1[2] = c->Ph[eqPG] +dpdx[eqPG]*(u1[0]-R_com(c));  // P-
      //u1[3] = c->Ph[eqPG] +dpdx[eqPG]*(u1[1]-R_com(c));  // P+
      //u1[4] = 2.0*(u1[1]*u1[3] - u1[0]*u1[2])/(u1[1]*u1[1]-u1[0]*u1[0]) -
      //            (u1[3]-u1[2])/(u1[1]-u1[0]);
      //c->dU[eqMMX] += FV_dt*u1[4];
      //cout <<"O2: r-="<<u1[0]<<", r+="<<u1[1]<<", p-="<<u1[2]<<", p+="<<u1[3]<<", <p/r>="<<u1[4]<<", old p/r="<<(c->Ph[eqPG] + dpdx[eqPG]*(CI.get_dpos(c,Rcyl)-R_com(c)))/CI.get_dpos(c,Rcyl)<<"\n";
      break;
     default:
      rep.error("Bad OOA in cyl_FV_solver_Hydro_Euler::dU, only know 1st,2nd",OA);
    }
  }

  return err;
}


// ##################################################################
// ##################################################################


////////////////////////////////////////////////////////////////////////
/// SPHERICAL (POLAR) COORDINATES
////////////////////////////////////////////////////////////////////////


// ##################################################################
// ##################################################################

sph_FV_solver_Hydro_Euler::sph_FV_solver_Hydro_Euler(
        const int nv, ///< number of variables in state vector.
        const int nd, ///< number of space dimensions in grid.
        const double cflno, ///< CFL number
        const double cellsize, ///< dx, cell size.
        const double gam,     ///< gas eos gamma.
        double *state, ///< State vector of mean values for simulation.
        const double avcoeff, ///< Artificial Viscosity Parameter etav.
        const int ntr         ///< Number of tracer variables.
        )
  : eqns_base(nv),
    //riemann_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv), riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_Hydro_Euler(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize),
    VectorOps_Sph(nd,cellsize)
{
#ifdef TESTING
  cout <<"sph_FV_solver_Hydro_Euler::sph_FV_solver_Hydro_Euler() constructor.\n";
#endif
  if (nd!=1) rep.error("Spherical coordinates only implemented for 1D \
                        spherical symmetry so far.  Sort it out!",nd);
  return;
}


// ##################################################################
// ##################################################################

sph_FV_solver_Hydro_Euler::~sph_FV_solver_Hydro_Euler()
{
  //  cout <<"sph_FV_solver_Hydro_Euler DESTRUCTOR.\n";
}


// ##################################################################
// ##################################################################

int sph_FV_solver_Hydro_Euler::dU_Cell(cell *c, ///< Current cell.
				       const axes d, ///< Which axis we are looking along.
				       const double *fn, ///< Negative direction flux.
				       const double *fp, ///< Positive direction flux.
				       const double *dpdx, ///< slope vector for cell c.
				       const int OA,      ///< spatial order of accuracy.
				       const double, ///< cell length dx.
				       const double  ///< cell TimeStep, dt.
				       )
{
  double u1[eq_nvar];
  //
  // This calculates the negative of the ith component of divergence
  //
  int err = DivStateVectorComponent(c,d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  //
  // Add source term for the radial direction (2p_g/r).
  //
  // This is more complicated than axisymmetry, but the equations
  // are easy enough to derive.
  //
  if (d==Rsph) {
    switch (OA) {
      case OA1:
      c->dU[eqMMX] += FV_dt*2.0*c->Ph[eqPG]/R3(c);
      //u1[0] = CI.get_dpos(c,Rsph)-0.5*VOdR;
      //u1[1] = u1[0]+VOdR;
      //u1[2] = 3.0*(u1[1]*u1[1]*c->Ph[eqPG] - u1[0]*u1[0]*c->Ph[eqPG])/(u1[1]*u1[1]*u1[1]-u1[0]*u1[0]*u1[0]);
      //if (fabs(u1[2]-2.0*c->Ph[eqPG]/R3(c))/u1[2] >1.0e-15) {
      //  cout <<"O1: rel.diff. Src Term: "<< (2.0*c->Ph[eqPG]/R3(c) -u1[2])/u1[2];
      //  cout <<"  old = "<< 2.0*c->Ph[eqPG]/R3(c)<<",  new = "<<u1[2]<<",  p="<<c->Ph[eqPG]<<"\n";
      //}
      //c->dU[eqMMX] += FV_dt*u1[2];
      break;

      case OA2:
      c->dU[eqMMX] += FV_dt*2.0*
       ( (c->Ph[eqPG]-dpdx[eqPG]*R_com(c))/R3(c) +dpdx[eqPG] );
      //u1[0] = CI.get_dpos(c,Rsph)-0.5*VOdR;              // r-
      //u1[1] = u1[0]+VOdR;                                // r+
      //u1[2] = c->Ph[eqPG] +dpdx[eqPG]*(u1[0]-R_com(c));  // P-
      //u1[3] = c->Ph[eqPG] +dpdx[eqPG]*(u1[1]-R_com(c));  // P+
      //u1[4] = 3.0*(u1[1]*u1[1]*u1[3] - u1[0]*u1[0]*u1[2])/(pow(u1[1],3.0)-pow(u1[0],3.0)) -
      //            (u1[3]-u1[2])/(u1[1]-u1[0]);
      //c->dU[eqMMX] += FV_dt*u1[4];
      //if (fabs(u1[4]-2.0*((c->Ph[eqPG]-dpdx[eqPG]*R_com(c))/R3(c) +dpdx[eqPG]))/u1[4] >1.0e-15)
      //  cout <<"O2: r-="<<u1[0]<<", r+="<<u1[1]<<", p-="<<u1[2]<<", p+="<<u1[3]<<", new <2p/r>="<<u1[4]<<", old <2p/r>="<<2.0*( (c->Ph[eqPG]-dpdx[eqPG]*R_com(c))/R3(c) +dpdx[eqPG] )<<"\n";
      break;
     default:
      rep.error("Bad OOA in sph_FV_solver_Hydro_Euler::dU, only know 1st,2nd",OA);
    }
  }


  //#define GRAVITY
#ifdef GRAVITY
  //
  // Here we impose an external gravitational acceleration, with a
  // radial power law and a normalisation K=G.M, so that: a=-K/r^alpha
  //
#define GRAV_NORM 1.33e29  // 1.989e37*6.67e-8 = 10^4 Msun*G
#define GRAV_ALPHA 2.0     // point mass has 1/r potential, 1/r^2 acc.

  //
  // Mean value of k.r^{-a} in a cell: for a!=3, cell centred at r=ri,
  // and cell diameter delta, this is
  // k.[(ri+delta/2)^{3-alpha}-(ri-delta/2)^{3-alpha}]/[(3-alpha)*delta*ri*R3]
  //
  // Note this is the same for 1st and 2nd order because we are
  // calculating it from the exact potential integrated through the
  // cell.
  //
  if (d==Rsph) {
    //cout <<"FVdx="<<FV_dx;
    double rc = CI.get_dpos(c,Rsph);
    // double temp = GRAV_NORM*(pow(rc+0.5*FV_dx,3.0-GRAV_ALPHA) -
    // 			     pow(rc-0.5*FV_dx,3.0-GRAV_ALPHA))
    //   /( (3.0-GRAV_ALPHA)*rc*FV_dx*R3(c) );
    double temp = c->Ph[eqRO]*GRAV_NORM/(rc*R3(c)); // GM/r^2
    if (rc<0.0) temp *= -1.0; // so that we always point towards the origin! 
    //cout <<"  r="<<rc<<"  GM/r^2="<<temp;
    //cout <<" d(MMX)/d(MMX,grav) = "<<c->dU[eqMMX]<<" / "<< -FV_dt*temp;
    c->dU[eqMMX] -= FV_dt*temp;
    //cout <<" d(ERG)/d(ERG,grav) = "<<c->dU[eqERG]<<" / "<< -FV_dt*temp*c->Ph[eqVX]<<"\n";
    c->dU[eqERG] -= FV_dt*c->Ph[eqVX]*temp;
  }

#endif // GRAVITY

  return err;
}


// ##################################################################
// ##################################################################


