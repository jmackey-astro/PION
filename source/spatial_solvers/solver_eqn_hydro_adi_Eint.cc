///
/// \file solver_eqn_hydro_adi_Eint.cc
/// \author Jonathan Mackey
///
/// Solver for the adiabatic Euler Equations which integrates the
/// INTERNAL ENERGY and NOT the total energy.  Calculates flux via
/// either Lax-Friedrichs or Riemann solver (linear and/or exact).
/// Also tracks flux of N passive tracers.
///
/// Modifications:\n
/// - 2010.12.28 JM: Started.
/// - 2010.12.30 JM: Can now integrate both internal and total energy.
/// - 2010.12.31 JM: updated logic in correcting negative pressure,
///    and in when to use the internal energy.  Also got spherical
///    and cylindrical correction terms for p.div(v) implemented
///    correctly (I think -- it will need further testing).
///
/// - 2011.01.03 JM: Added PostProcess_dU() for artificial viscosity.
///   Added functions to calculate the Q-viscosity.  Added
///   preprocess_dU function to set the Q,div(v) values (moved from
///   gridMethods.cc).
///
/// - 2011.01.04 JM: Added viscous terms for cylindrical/spherical
///   coordinates, but not working yet.  Needs testing/fixing.

#include "../defines/functionality_flags.h"
#ifdef INCLUDE_EINT_ADI_HYDRO

//
// I don't really need to include all of these, but this tells me
// which files the class depends on...
//

#include "../global.h"
//#include "../grid/cell_interface.h"
//#include "../equations/eqns_base.h"
//#include "../equations/eqns_hydro_adiabatic.h"
#include "eqns_hydro_adi_Eint.h"

#include "../Riemann_solvers/riemann.h"
#include "../Riemann_solvers/Roe_Hydro_PrimitiveVar_solver.h"

//#include "../flux_calc/flux_base.h"
//#include "../flux_calc/flux_hydro_adiabatic.h"
#include "../flux_calc/flux_hydro_adi_Eint.h"

//#include "solver_eqn_base.h"
//#include "solver_eqn_hydro_adi.h"
#include "solver_eqn_hydro_adi_Eint.h"

#include "../coord_sys/VectorOps.h"
#include "../coord_sys/VectorOps_spherical.h"

using namespace std;



FV_solver_Hydro_Euler_Eint::FV_solver_Hydro_Euler_Eint(const int nv,
 ///< number of variables in state vector.
					     const int nd,
 ///< number of space dimensions in grid.
					     const double cflno,
 ///< CFL number
					     const double cellsize,
 ///< dx, cell size.
					     const double gam,
 ///< gas eos gamma.
					     double *state,
 ///< State vector of mean values for simulation.
					     const double avcoeff,
 ///< Artificial Viscosity Parameter etav.
					     const int ntr
 ///< Number of tracer variables.
					     )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv),
    riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_Hydro_Euler(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_Euler_Eint(nv),
    flux_solver_hydro_adi_Eint(nv,state,avcoeff,gam,ntr)
    
{
  cout <<"FV_solver_Hydro_Euler_Eint::FV_solver_Hydro_Euler_Eint() constructor\n";
  Qvisc_eta = 4.0;
  return;
}

FV_solver_Hydro_Euler_Eint::~FV_solver_Hydro_Euler_Eint()
{
  cout <<"FV_solver_Hydro_Euler_Eint::~FV_solver_Hydro_Euler_Eint() destructor.\n";
  return;
}

///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_Hydro_Euler_Eint::dU_Cell(cell *c,          // Current cell.
				   const axes d,     // Which axis we are looking along.
				   const double *fn, // Negative direction flux.
				   const double *fp, // Positive direction flux.
				   const double *dpdx, // slope vector for cell c.
				   const int OA,      // spatial order of accuracy.
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

#ifdef        EINT_ETOT_PARALLEL
///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_Hydro_Euler_Eint::CellAdvanceTime(class cell *c,
 ///< current cell
					   const double *Pin,
 // Initial State Vector.
					   double *dU,
 // Update vector dU
					   double *Pf,
 // Final state vector (can be same as initial vec.).
					   double *dE,
 // TESTING Tracks change of energy for negative pressure correction.
					   const double,
 // gas EOS gamma.
					   const double
  // Cell timestep dt.
					   )
{
  double u1[eq_nvar];
#ifdef TESTING
  double u2[eq_nvar];
#endif //TESTING
  static int counter=0;

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
    if (counter<1000){
      cout<<"(FV_solver_Hydro_Euler_Eint::CellAdvanceTime) UtoP complained \
           (maybe about negative pressure...) fixing\n";
    }
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

#ifdef        EINT_ETOT_PARALLEL
  double divv = CI.get_DivV(c);
  //cout <<"divv="<<divv<<"\n";
  if (divv>0.0) {
    //if (Pf[eqEINT]<0.0) {
    //  cout << "rarefaction, so using Eint. dv="<<divv;
    //  cout << " P="<<Pf[eqPG]<<" from eint, P="<<Pf[eqEINT] <<"\n";
    //}
    if (Pf[eqPG] <=0.5*Pf[eqEINT])
      Pf[eqPG] = Pf[eqEINT];
    //Pf[eqPG] = max(Pf[eqPG],Pf[eqEINT]);
  } 
#endif //     EINT_ETOT_PARALLEL


  //  for (int v=0;v<eq_nvar;v++) dU[v] = 0.; // Reset the dU array for the next timestep.

  return 0;
}
#endif //     EINT_ETOT_PARALLEL

double FV_solver_Hydro_Euler_Eint::CellTimeStep(const cell *c, ///< pointer to cell
						const double, ///< gas EOS gamma.
						const double  ///< Cell size dx.
						)
{
  //
  // First we call the basic hydro timestep calculation.
  //
  double deltat = FV_solver_Hydro_Euler::CellTimeStep(c,0.0,0.0);

  if (!(c->isgd)) {
    return FV_dt; // should be set in Euler version.
  }
  //
  // Next we calculate the local velocity divergence and a timestep
  // limit imposed by viscous diffusion, eqn.78 in Stone & Norman
  // (1992,ApJS,80,753).
  //
  int indices[MAX_DIM];
  indices[0]=eqVX; indices[1]=eqVY; indices[2]=eqVZ;
  double divv=fabs(Div(c,1,indices));
  // if (FV_cfl/(4.0*Qvisc_eta*divv) < deltat) {
  //   cout <<" Hydro step="<<deltat<<"  Visc.step=";
  //   cout <<FV_cfl/(4.0*Qvisc_eta*divv)<<"\n";
  // }
  deltat = min(deltat,FV_cfl/(4.0*Qvisc_eta*divv));
  
  //
  // Set the internal data value to deltat (AAARRGH! THESE ARE
  // DANGEROUS!  I SHOULD GET RID OF THIS.
  //
  FV_dt = deltat;

  return FV_dt;
}


//
// Multi-dimensional calculations to be performed on every cell before
// the flux calculations.
//
int FV_solver_Hydro_Euler_Eint::preprocess_data(const int csp, const int ctm)
{
  //    cout <<"HELLO!!!\n";
  //
  // First call the base solver (which sets div(v) if Lapidus
  // viscosity or the H-correction eta values in each direction if
  // using the H-correction.
  //
  int err = FV_solver_base::preprocess_data(csp,ctm);


#ifdef  EINT_ETOT_PARALLEL
#ifndef VNR_VISCOSITY
  //
  // If not using Q-viscosity then we still need div(v) for the energy update.
  //cout <<"\t\t Eint: Calculating div(v).\n";
  class cell* c = grid->FirstPt();
  do {
    set_div_v(c);
  } while ( (c =grid->NextPt(c)) !=0);
  //
#endif // not VNR_VISCOSITY
#endif //EINT_ETOT_PARALLEL

#ifdef VNR_VISCOSITY
  //
  // von Neumann & Richtmeyer multi-D viscosity needs div(v) and Ndim
  // Q values (the diagonal elements of a traceless rank 2 tensor).
  //
  // WARNING!! THIS REALLY NEEDS TO BE CALCULATED FOR THE FIRST
  // BOUNDARY CELLS TOO IN ORDER TO WORK FOR PARALLEL GRIDS AND
  // PERIODIC DOMAINS.  NEED TO FIX THE GRID FOR THIS!
  //
  //if (SimPM.artviscosity==AV_VonNeuRicht) {
  class cell* c = grid->FirstPt();
  //    cout <<"HELLO!!!\n";
  SetDirection(XX); // just to make sure!
  double divv=0.0;
  if (ctm!=SimPM.tmOOA) {
    do {
      set_div_v(c);
    } while ( (c =grid->NextPt(c)) !=0);
  }
  else {
    do {
      set_div_v(c);
      divv = CI.get_DivV(c);
      set_Qvisc(c,XX,divv);
      if (FV_gndim>1)
	set_Qvisc(c,YY,divv);
      if (FV_gndim>2)
	set_Qvisc(c,ZZ,divv);
    } while ( (c =grid->NextPt(c)) !=0);
  }
  //} // if AV_VonNeuRicht
#endif // VNR_VISCOSITY
  
  return err;
} // preprocess_data()



void FV_solver_Hydro_Euler_Eint::set_Qvisc(cell *c,
					   ///< cell to operate on.
					   const axes ax,
					   ///< axis we are looking along (Q1/2/3).
					   const double divv
					   ///< value of divv (should be pre-calculated).
					   )
{
  //
  // Now if div(v)>0 then Q=0, so we check that here:
  //
  if (divv>=0.0) {
    CI.set_Hcorr(c,ax,0.0);
    return;
  }

  //
  // This is for Cartesian geometries, where there are no geometric
  // terms!  First add -1/3 div(v)
  //
  double value = -divv/3.0;

  //cout <<"divv="<<divv;
  //
  // Now we need to get a gradient.
  //
  cell *cn, *cp; int v=-1;
  switch (ax) {
  case XX:
    cn = grid->NextPt(c,XN);
    cp = grid->NextPt(c,XP);
    v = eqVX;
    break;
  case YY:
    cn = grid->NextPt(c,YN);
    cp = grid->NextPt(c,YP);
    v = eqVY;
    break;
  case ZZ:
    cn = grid->NextPt(c,ZN);
    cp = grid->NextPt(c,ZP);
    v = eqVZ;
    break;
  default:
    rep.error("Bad axes",ax);
  }
  //
  // WON'T NEED TO CHECK IF THEY EXIST ONCE I UPDATE THE GRID TO HAVE
  // CORNER CELLS.  CAN REMOVE THIS THEN.
  //
  if      (cn && cp) {
    //
    // Add velocity gradient along requested direction.
    //
    value += 0.5*(cp->Ph[v]-cn->Ph[v])/FV_dx;
    //cout <<" grad(v)="<<0.5*(cp->Ph[v]-cn->Ph[v])/FV_dx;
  }
  else if (!cn && !cp) {
    CI.print_cell(c);
    rep.error("isolated cell -- can't get gradient",c->id);
  }
  else if (!cn) {
    value += (cp->Ph[v]-c->Ph[v])/FV_dx;
  }
  else { // if (!cp)
    value += (c->Ph[v]-cn->Ph[v])/FV_dx;
  }

  //
  // Now multiply the lot by the prefactor.
  //
  value *=Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO]*divv;
  //cout <<" prefactor="<<Qvisc_eta*Qvisc_eta*c->Ph[eqRO]*divv;
  //cout <<" ... total="<<value<<"\n";
  //
  // Finally add this to the Cell's extra data using the H-correction
  // interface.
  //
  CI.set_Hcorr(c,ax,value);
  return;
}

double FV_solver_Hydro_Euler_Eint::get_Qvisc(const cell *c,
					     ///< cell to get Q from.
					     const axes ax
					     ///< get Q for this axis.
					     )
{
  return CI.get_Hcorr(c,ax);
}


int FV_solver_Hydro_Euler_Eint::PostProcess_dU(const int csp,
					       ///< Space order of acc for this call.
					       const int ctm
					       ///< Time order of acc for this call.
					       )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_Hydro_Euler_Eint::PostProcess_dU ...starting.\n";
#endif //FUNCTION_ID
  int err=0;

#ifdef VNR_VISCOSITY
  //
  // Loop through every real grid cell and add viscous contributions.
  // Refs:
  // Tscharnuter & Winkler (1979) Computer Physics Communications, 18,171.
  // Stone & Norman (1992) ApJS, 80, 753.
  //
  // Qii = l^2 rho div(v)*(grad(v)-div(v)/3)
  // Momentum corrections are divergence of each of these components.
  // Energy correction is (l^2 rho div(v)/3)[velocity gradients].
  // 
  // l is the viscosity coefficient in units of the grid spacing dx.
  //
  SetDirection(XX);
  class cell* c = grid->FirstPt();
  double prefactor, QP, QN, QZ;
  cell *cp, *cn;
  do {
    //
    // Subtract off P.div(v) term.
    //
#ifdef EINT_ETOT_PARALLEL
    c->dU[eqEINT] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#else
    c->dU[eqERG] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#endif // EINT_ETOT_PARALLEL

    //
    // div(v), Qii should be stored and pre-calculated.  Only add
    // viscosity for the full update.  The Q-viscosity is non-zero
    // only for compression, so we first check div(v):
    //
    if (ctm==SimPM.tmOOA && CI.get_DivV(c) <0.0) {
      //if (CI.get_DivV(c) <0.0) {
      //
      // First the x-direction: take 2nd order central differences for
      // divergence and gradient.
      // The momenta get modified by the ith component of div(Q)
      //
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);
      if (!cn || !cp) rep.error("No boundary data!",cn);
      QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,XX) : 0.0;
      QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,XX) : 0.0;
      //cout <<"QX: dU="<< c->dU[eqMMX] <<" corr="<< FV_dt*0.5*(QP-QN)/FV_dx<<"\n";
      c->dU[eqMMX] -= FV_dt*0.5*(QP-QN)/FV_dx;
      //
      // Y-direction, if relevant.
      //
      if (FV_gndim>1) {
	cp = grid->NextPt(c,YP);
	cn = grid->NextPt(c,YN);
	if (!cn || !cp) rep.error("No boundary data!",cn);
	QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,YY) : 0.0;
	QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,YY) : 0.0;
	c->dU[eqMMY] -= FV_dt*0.5*(QP-QN)/FV_dx;
      }
      //
      // Z-direction, if relevant.
      //
      if (FV_gndim>2) {
	cp = grid->NextPt(c,ZP);
	cn = grid->NextPt(c,ZN);
	if (!cn || !cp) rep.error("No boundary data!",cn);
	QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,ZZ) : 0.0;
	QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,ZZ) : 0.0;
	c->dU[eqMMZ] -= FV_dt*0.5*(QP-QN)/FV_dx;
      }
      //
      // Now the energy correction:
      // We need partial derivatives d(vx)/dx, d(vy)/dy, d(vz)/dz.
      // Use QN, QP, QZ for the three derivatives
      // We have already bugged out if the neighbour cells don't 
      // exist, so no need to check again.
      //
      prefactor = Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO]*CI.get_DivV(c)/3.0;
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);    
      QN = 0.5*(cp->Ph[eqVX]-cn->Ph[eqVX])/FV_dx;
      if (FV_gndim>1) {
	cp = grid->NextPt(c,YP);
	cn = grid->NextPt(c,YN);    
	QP = 0.5*(cp->Ph[eqVY]-cn->Ph[eqVY])/FV_dx;
      }
      else {
	QP = 0.0;
      }
      if (FV_gndim>2) {
	cp = grid->NextPt(c,ZP);
	cn = grid->NextPt(c,ZN);    
	QZ = 0.5*(cp->Ph[eqVZ]-cn->Ph[eqVZ])/FV_dx;
      }
      else {
	QZ = 0.0;
      }
#ifdef EINT_ETOT_PARALLEL
      c->dU[eqEINT] -= FV_dt*prefactor*
	( (QN-QP)*(QN-QP) +(QN-QZ)*(QN-QZ) +(QZ-QP)*(QZ-QP) );
#else
      //cout <<"QE: dU="<< c->dU[eqERG] <<" corr="<< FV_dt*prefactor*
      //	( (QN-QP)*(QN-QP) +(QN-QZ)*(QN-QZ) +(QZ-QP)*(QZ-QP) ) <<"\n";
      c->dU[eqERG]  -= FV_dt*prefactor*
	( (QN-QP)*(QN-QP) +(QN-QZ)*(QN-QZ) +(QZ-QP)*(QZ-QP) );
#endif // EINT_ETOT_PARALLEL
    } // if we have compression.


  } while ( (c =grid->NextPt(c)) !=0);
  
#endif // VNR_VISCOSITY

#ifdef FUNCTION_ID
  cout <<"FV_solver_Hydro_Euler_Eint::PostProcess_dU ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


////////////////////////////////////////////////////////////////////////
/// CYLINDRICAL (AXISYMMETRIC) COORDINATES
////////////////////////////////////////////////////////////////////////


cyl_FV_solver_Hydro_Euler_Eint::cyl_FV_solver_Hydro_Euler_Eint(const int nv,
 ///< number of variables in state vector.
					     const int nd,
 ///< number of space dimensions in grid.
					     const double cflno,
 ///< CFL number
					     const double cellsize,
 ///< dx, cell size.
					     const double gam,
 ///< gas eos gamma.
					     double *state,
 ///< State vector of mean values for simulation.
					     const double avcoeff,
 ///< Artificial Viscosity Parameter etav.
					     const int ntr
 ///< Number of tracer variables.
					     )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv),
    riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_Hydro_Euler(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_Euler_Eint(nv),
    flux_solver_hydro_adi_Eint(nv,state,avcoeff,gam,ntr),
    FV_solver_Hydro_Euler_Eint(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
  cout <<"cyl_FV_solver_Hydro_Euler_Eint::cyl_FV_solver_Hydro_Euler_Eint() constructor\n";
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for 2d axial symmetry \
                        so far.  Sort it out!",nd);
  cout <<"\n";
  cout <<"\t\t************************\n";
  cout <<"\tINTERNAL ENERGY SOLVER NOT YET WORKING IN AXISYMMETRY.\n";
  cout <<"\tTHE CODE IS THERE, BUT THE VISCOSITY IS NOT RIGHT.\n";
  cout <<"\tFEEL FREE TO TRY AND GET IT WORKING BUT BE WARNED: THE\n";
  cout <<"\tSOLVER IS MUCH MORE DIFFUSIVE THAN THE TOTAL ENERGY ONE.\n";
  cout <<"\t\t************************\n";
  rep.error("solver is broken.",-9);
  return;
}

cyl_FV_solver_Hydro_Euler_Eint::~cyl_FV_solver_Hydro_Euler_Eint()
{
  cout <<"cyl_FV_solver_Hydro_Euler_Eint::~cyl_FV_solver_Hydro_Euler_Eint() destructor.\n";
  return;
}

///
/// Adds the contribution from flux in the current direction to dU.
///
int cyl_FV_solver_Hydro_Euler_Eint::dU_Cell(cell *c,          // Current cell.
				   const axes d,     // Which axis we are looking along.
				   const double *fn, // Negative direction flux.
				   const double *fp, // Positive direction flux.
				   const double *dpdx, // slope vector for cell c.
				   const int OA,      // spatial order of accuracy.
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
  
  //
  // Add geometric source term for the radial momentum:
  //
  if (d==Rcyl) {
    switch (OA) {
     case OA1:
       c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]/CI.get_dpos(c,Rcyl));
      break;
     case OA2:
       c->dU[eqMMX] += FV_dt*(c->Ph[eqPG] + dpdx[eqPG]*(CI.get_dpos(c,Rcyl)-R_com(c)))/CI.get_dpos(c,Rcyl);
      break;
     default:
      rep.error("Bad OOA in cyl_FV_solver_Hydro_Euler_Eint::dU, only know 1st,2nd",OA);
    }
  }

  //
  // The p.div(v) source term in the energy equation used to be here,
  // but I moved it to PostProcess_dU() so I can use the
  // pre-calculated div(v).
  //

  return(err);
}
 
 
void cyl_FV_solver_Hydro_Euler_Eint::set_Qvisc(cell *c,
					   ///< cell to operate on.
					   const axes ax,
					   ///< axis we are looking along (Q1/2/3).
					   const double divv
					   ///< value of divv (should be pre-calculated).
					   )
{
  //cout <<"hello!";
  //
  // This is based on the cartesian version, but radial gradients must
  // be handled differently.
  //
  // double value=0.0, dx=0.0;
  // cell *cn, *cp; int v=-1;
  // switch (ax) {
  // case Zcyl:
  //   cn = grid->NextPt(c,XN);
  //   cp = grid->NextPt(c,XP);
  //   v = eqVX;
  //   break;
  // case Rcyl:
  //   cn = grid->NextPt(c,YN);
  //   cp = grid->NextPt(c,YP);
  //   v = eqVY;
  //   break;
  // case ZZ:
  //   rep.error("Bad axis in cyl_FV_solver_Hydro_Euler_Eint::set_Qvisc()",ax);
  //   cn = grid->NextPt(c,ZN);
  //   cp = grid->NextPt(c,ZP);
  //   v = eqVZ;
  //   break;
  // default:
  //   rep.error("Bad axes",ax);
  // }
  // //
  // // Now make sure we have cells to use to calculate divergence.
  // //
  // //dx = 2.0*FV_dx;
  // if (!cn && !cp) {
  //   CI.print_cell(c);
  //   rep.error("isolated cell -- can't get gradient",c->id);
  // }
  // else if (!cn) {
  //   cn = c;
  //   //dx *= 0.5;
  // }
  // else if (!cp) {
  //   cp = c;
  //   //dx *= 0.5;
  // }

  // if (ax==Zcyl) {
  //   value += (cp->Ph[v]-cn->Ph[v]); // /dx;
  // }
  // else {
  //   value += (cp->Ph[v]-cn->Ph[v]); // /(R_com(cp)-R_com(cn));
  // }
  
  // //
  // // So now we have gradient of vel. along direction.
  // // We need to square it and multiply by l^2 rho.
  // //
  // if (value >= 0.0) {
  //   CI.set_Hcorr(c,ax,0.0);
  // }
  // else {
  //   value *= value;
  //   //value *= Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO];
  //   value *= Qvisc_eta*Qvisc_eta*c->Ph[eqRO];
  //   CI.set_Hcorr(c,ax,value);
  // }
  // return;




  //
  // Now if div(v)>0 then Q=0, so we check that here:
  //
  if (divv>=0.0) {
    CI.set_Hcorr(c,ax,0.0);
    return;
  }

  //
  //  First add -1/3 div(v)
  //
  double value = -divv/3.0;

  //cout <<"divv="<<divv;
  //
  // Now we need to get a gradient.
  //
  cell *cn, *cp; int v=-1;
  switch (ax) {
  case Zcyl:
    cn = grid->NextPt(c,XN);
    cp = grid->NextPt(c,XP);
    v = eqVX;
    break;
  case Rcyl:
    cn = grid->NextPt(c,YN);
    cp = grid->NextPt(c,YP);
    v = eqVY;
    break;
  case ZZ:
    rep.error("Bad axis in cyl_FV_solver_Hydro_Euler_Eint::set_Qvisc()",ax);
    cn = grid->NextPt(c,ZN);
    cp = grid->NextPt(c,ZP);
    v = eqVZ;
    break;
  default:
    rep.error("Bad axes",ax);
  }
  //
  // WON'T NEED TO CHECK IF THEY EXIST ONCE I UPDATE THE GRID TO HAVE
  // CORNER CELLS.  CAN REMOVE THIS THEN.
  //
  if      (cn && cp) {
    //
    // Add velocity gradient along requested direction.
    //
    if (ax==Zcyl) {
      value += 0.5*(cp->Ph[v]-cn->Ph[v])/FV_dx;
    }
    else {
      value += (cp->Ph[v]-cn->Ph[v])/(R_com(cp)-R_com(cn));
    }
    //cout <<" grad(v)="<<0.5*(cp->Ph[v]-cn->Ph[v])/FV_dx;
  }
  else if (!cn && !cp) {
    CI.print_cell(c);
    rep.error("isolated cell -- can't get gradient",c->id);
  }
  else if (!cn) {
    if (ax==Zcyl) {
      value += (cp->Ph[v]-c->Ph[v])/FV_dx;
    }
    else {
      value += (cp->Ph[v]-c->Ph[v])/(R_com(cp)-R_com(c));
    }
  }
  else { // if (!cp)
    if (ax==Zcyl) {
      value += (c->Ph[v]-cn->Ph[v])/FV_dx;
    }
    else {
      value += (c->Ph[v]-cn->Ph[v])/(R_com(c)-R_com(cn));
    }
  }

  //
  // Now multiply the lot by the prefactor.
  //
  value *=Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO]*divv;
  //cout <<" prefactor="<<Qvisc_eta*Qvisc_eta*c->Ph[eqRO]*divv;
  //cout <<" ... total="<<value<<"\n";
  //
  // Finally add this to the Cell's extra data using the H-correction
  // interface.
  //
  CI.set_Hcorr(c,ax,value);
  return;
}

int cyl_FV_solver_Hydro_Euler_Eint::PostProcess_dU(const int csp,
					       ///< Space order of acc for this call.
					       const int ctm
					       ///< Time order of acc for this call.
					       )
{
#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_Hydro_Euler_Eint::PostProcess_dU ...starting.\n";
#endif //FUNCTION_ID
  int err=0;

#ifdef VNR_VISCOSITY
  //
  // Loop through every real grid cell and add viscous contributions.
  // Refs:
  // Tscharnuter & Winkler (1979) Computer Physics Communications, 18,171.
  // Stone & Norman (1992) ApJS, 80, 753.
  //
  // Qii = l^2 rho div(v)*(grad(v)-div(v)/3)
  // Momentum corrections are divergence of each of these components.
  // Energy correction is (l^2 rho div(v)/3)[velocity gradients].
  // 
  // l is the viscosity coefficient in units of the grid spacing dx.
  //
  SetDirection(XX);
  class cell* c = grid->FirstPt();
  double prefactor, QP, QN, Q3;
  cell *cp, *cn; // cells to the right and left of cell c.
  double rp, rn; // cell centre-of-volume in radial dirn. for cn,cp.
  do {
    //
    // Subtract off P.div(v) term.
    //
#ifdef EINT_ETOT_PARALLEL
    c->dU[eqEINT] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#else
    c->dU[eqERG] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#endif // EINT_ETOT_PARALLEL

    // ******************************************************
    // Let's try something simpler -- the vN-R 1D viscosity.
    // ******************************************************
    if (1==2 && 
	ctm==SimPM.tmOOA && 
	c->isgd 
	//&& CI.get_DivV(c) <0.0
	) {
      //
      // Z-momentum:
      //
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);
      if (!cn || !cp) rep.error("No boundary data!",cn);
      QP = get_Qvisc(cp,XX);
      QN = get_Qvisc(cn,XX);
      c->dU[eqMMX] -= FV_dt*0.5*(QP-QN)/FV_dx;
      //
      // R-momentum
      //
      if (FV_gndim>1) {
	cp = grid->NextPt(c,YP);
	cn = grid->NextPt(c,YN);
	if (!cn || !cp) rep.error("No boundary data!",cn);
	//
	// CYLINDRICAL: GRADIENT NEEDS CENTRE-OF-VOL. POSITIONS.
	//
	rp = R_com(cp);
	rn = R_com(cn);
	QP = get_Qvisc(cp,YY);
	QN = get_Qvisc(cn,YY);
	c->dU[eqMMY] -= 
	  //  FV_dt*( 2.0*(rp*QP-rn*QN)/(rp*rp-rn*rn));
	  FV_dt*( (QP-QN)/(rp-rn) );
      }
      //
      // Energy: - Qz*dvz/dz -QR*dvR/dR
      //
#ifdef EINT_ETOT_PARALLEL
      c->dU[eqEINT] -= get_Qvisc(c,Zcyl)*0.5*
	(grid->NextPt(c,XP)->Ph[eqVX]-grid->NextPt(c,XN)->Ph[eqVX])/FV_dx;
      if (FV_gndim>1) {
	c->dU[eqEINT] -= get_Qvisc(c,Rcyl)*
	  (cp->Ph[eqVY]-cn->Ph[eqVY])/(rp-rn);
      }
#else
      c->dU[eqERG] -= get_Qvisc(c,Zcyl)*0.5*
	(grid->NextPt(c,XP)->Ph[eqVX]-grid->NextPt(c,XN)->Ph[eqVX])/FV_dx;
      if (FV_gndim>1) {
	c->dU[eqERG] -= get_Qvisc(c,Rcyl)*
	  (cp->Ph[eqVY]-cn->Ph[eqVY])/(rp-rn);
      }
#endif // EINT_ETOT_PARALLEL
    }
    // ******************************************************
    // something simpler -- the vN-R 1D viscosity... not really working.
    // ******************************************************

    // ******************************************************
    // TENSOR Q-VISCOSITY -- ALSO NOT REALLY WORKING!
    // ******************************************************
    //
    // div(v), Qii should be stored and pre-calculated.
    //
    // - Only add viscosity for the full 2nd order time update.
    // - Only do the calculation for grid-cells (boundary data
    // irrelevant).
    // - The Q-viscosity is non-zero only for compression, so we first
    // check div(v):
    //
    if (//1==2 && 
	//ctm==SimPM.tmOOA && 
	c->isgd && CI.get_DivV(c) <0.0) {
      //if (CI.get_DivV(c) <0.0) {
      //
      // First the x-direction: take 2nd order central differences for
      // divergence and gradient.
      // The momenta get modified by the ith component of div(Q)
      //
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);
      if (!cn || !cp) rep.error("No boundary data!",cn);
      QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,XX) : 0.0;
      QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,XX) : 0.0;
      //cout <<"QX: dU="<< c->dU[eqMMX] <<" corr="<< FV_dt*0.5*(QP-QN)/FV_dx<<"\n";
      c->dU[eqMMX] -= FV_dt*0.5*(QP-QN)/FV_dx;

      //
      // R-direction, if relevant.
      //
      if (FV_gndim>1) {
	cp = grid->NextPt(c,YP);
	cn = grid->NextPt(c,YN);
	if (!cn || !cp) rep.error("No boundary data!",cn);
	//
	// CYLINDRICAL: GRADIENT NEEDS CENTRE-OF-VOL. POSITIONS.
	//
	rp = R_com(cp);
	rn = R_com(cn);
	QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,YY) : 0.0;
	QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,YY) : 0.0;
	//
	// Here we are taking a divergence, so we need the radial 
	// divergence equation, and there is a source term too.
	// So we want d(QRR)/dR + (2(QRR)+(QZZ))/R
	// I am going to ignore the extra correction to the second
	// term which would come from the gradient of QRR over the
	// cell.
	// ************ DIFFERENCE FROM CARTESIAN! *************
	//
	// CYLINDRICAL: GRADIENT NEEDS CENTRE-OF-VOL. POSITIONS.
	//
	c->dU[eqMMY] -= 
	  //FV_dt*( (QP-QN)/(rp-rn) 
	  //	  +(2.0*get_Qvisc(c,Rcyl) +get_Qvisc(c,Zcyl))/R_com(c)
	  //	  );
	  FV_dt*( 3.0*(rp*rp*QP-rn*rn*QN)/(pow(rp,3.0)-pow(rn,3.0)) 
		  +get_Qvisc(c,Zcyl)/R_com(c)
		  );
      }
      //
      // Z-direction is not implemented.
      //

      //
      // Now the energy correction: We need partial derivatives
      // d(vx)/dx, d(vy)/dy.  Use QN, QP for the two derivatives.
      // Because they are covariant derivatives, the 3rd component has
      // a 'source term' of V_R/R even though there can be no theta
      // derivatives.
      // We have already bugged out if the neighbour cells don't 
      // exist, so no need to check again.
      //
      prefactor = Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO]*CI.get_DivV(c)/3.0;
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);    
      QN = 0.5*(cp->Ph[eqVX]-cn->Ph[eqVX])/FV_dx;
      if (FV_gndim>1) {
	cp = grid->NextPt(c,YP);
	cn = grid->NextPt(c,YN);
	//
	// CYLINDRICAL: GRADIENT NEEDS CENTRE-OF-VOL. POSITIONS.
	//
	QP = (cp->Ph[eqVY]-cn->Ph[eqVY])/(rp-rn);
	//
	// CYLINDRICAL: ALSO NEED SOURCE TERM FOR 3RD COMPONENT.
	//
	Q3 = c->Ph[eqVY]/R_com(c);
      }
      else {
	QP = 0.0;
	Q3 = 0.0;
      }
#ifdef EINT_ETOT_PARALLEL
      c->dU[eqEINT] -= FV_dt*prefactor*
	( (QN-QP)*(QN-QP) +(QN-Q3)*(QN-Q3) +(Q3-QP)*(Q3-QP) );
#else
      //cout <<"QE: dU="<< c->dU[eqERG] <<" corr="<< FV_dt*prefactor*
      // 	( (QN-QP)*(QN-QP) +(QN*QN) +(QP*QP) ) <<"\n";
      c->dU[eqERG]  -= FV_dt*prefactor*
	//( (QN-QP)*(QN-QP) +(QN*QN) +(QP*QP) );
	( (QN-QP)*(QN-QP) +(QN-Q3)*(QN-Q3) +(Q3-QP)*(Q3-QP) );
#endif // EINT_ETOT_PARALLEL
    // ******************************************************
    // TENSOR Q-VISCOSITY -- ALSO NOT REALLY WORKING!
    // ******************************************************

    } // if we have compression.


  } while ( (c =grid->NextPt(c)) !=0);
  
#endif // VNR_VISCOSITY

#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_Hydro_Euler_Eint::PostProcess_dU ...returning.\n";
#endif //FUNCTION_ID
  return err;
}

////////////////////////////////////////////////////////////////////////
/// SPHERICAL (POLAR) COORDINATES
////////////////////////////////////////////////////////////////////////

sph_FV_solver_Hydro_Euler_Eint::sph_FV_solver_Hydro_Euler_Eint(const int nv,
 ///< number of variables in state vector.
					     const int nd,
 ///< number of space dimensions in grid.
					     const double cflno,
 ///< CFL number
					     const double cellsize,
 ///< dx, cell size.
					     const double gam,
 ///< gas eos gamma.
					     double *state,
 ///< State vector of mean values for simulation.
					     const double avcoeff,
 ///< Artificial Viscosity Parameter etav.
					     const int ntr
 ///< Number of tracer variables.
					     )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_Euler(nv),
    riemann_Euler(nv,state,gam),
    Riemann_FVS_Euler(nv,gam),
    Riemann_Roe_Hydro_PV(nv,gam),
    Riemann_Roe_Hydro_CV(nv,gam),
    flux_solver_hydro_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_Hydro_Euler(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_Euler_Eint(nv),
    flux_solver_hydro_adi_Eint(nv,state,avcoeff,gam,ntr),
    FV_solver_Hydro_Euler_Eint(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize),
    VectorOps_Sph(nd,cellsize)
{
  cout <<"sph_FV_solver_Hydro_Euler_Eint::sph_FV_solver_Hydro_Euler_Eint() constructor\n";
  if (nd!=1) rep.error("Spherical coordinates only implemented for 1D spherical symmetry \
                        so far.  Sort it out!",nd);
  rep.error("SOLVER WITH INTERNAL ENERGY NOT REALLY WORKING FOR SPHERICAL.",-8);
  return;
}

sph_FV_solver_Hydro_Euler_Eint::~sph_FV_solver_Hydro_Euler_Eint()
{
  cout <<"sph_FV_solver_Hydro_Euler_Eint::~sph_FV_solver_Hydro_Euler_Eint() destructor.\n";
  return;
}

///
/// Adds the contribution from flux in the current direction to dU.
///
int sph_FV_solver_Hydro_Euler_Eint::dU_Cell(cell *c,          // Current cell.
				   const axes d,     // Which axis we are looking along.
				   const double *fn, // Negative direction flux.
				   const double *fp, // Positive direction flux.
				   const double *dpdx, // slope vector for cell c.
				   const int OA,      // spatial order of accuracy.
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
  
  //
  // Add geometric source term for the radial momentum:
  //
  // This is more complicated than axisymmetry, but the equations
  // are easy enough to derive.
  //
  if (d==Rsph) {
    switch (OA) {
     case OA1:
       c->dU[eqMMX] += FV_dt*2.0*c->Ph[eqPG]/R3(c);
      break;
     case OA2:
       c->dU[eqMMX] += FV_dt*2.0*
	 ( (c->Ph[eqPG]-dpdx[eqPG]*R_com(c))/R3(c) +dpdx[eqPG] );
       break;
     default:
      rep.error("Bad OOA in sph_FV_solver_Hydro_Euler_Eint::dU, only know 1st,2nd",OA);
    }
  }

  //
  // The p.div(v) source term in the energy equation used to be here,
  // but I moved it to PostProcess_dU() so I can use the
  // pre-calculated div(v).
  //

  return(err);
}

int sph_FV_solver_Hydro_Euler_Eint::PostProcess_dU(const int csp,
					       ///< Space order of acc for this call.
					       const int ctm
					       ///< Time order of acc for this call.
					       )
{
#ifdef FUNCTION_ID
  cout <<"sph_FV_solver_Hydro_Euler_Eint::PostProcess_dU ...starting.\n";
#endif //FUNCTION_ID
  int err=0;

#ifdef VNR_VISCOSITY
  //
  // Loop through every real grid cell and add viscous contributions.
  // Refs:
  // Tscharnuter & Winkler (1979) Computer Physics Communications, 18,171.
  // Stone & Norman (1992) ApJS, 80, 753.
  //
  // Qii = l^2 rho div(v)*(grad(v)-div(v)/3)
  // Momentum corrections are divergence of each of these components.
  // Energy correction is (l^2 rho div(v)/3)[velocity gradients].
  // 
  // l is the viscosity coefficient in units of the grid spacing dx.
  //
  SetDirection(XX);
  class cell* c = grid->FirstPt();
  double prefactor, QP, QN;
  cell *cp, *cn;
  double rp, rn;
  do {
    //
    // Subtract off P.div(v) term.
    //
#ifdef EINT_ETOT_PARALLEL
    c->dU[eqEINT] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#else
    c->dU[eqERG] -= FV_dt*c->Ph[eqPG]*CI.get_DivV(c);
#endif // EINT_ETOT_PARALLEL

    //
    // div(v), Qii should be stored and pre-calculated.
    //
    // - Only add viscosity for the full 2nd order time update.
    // - Only do the calculation for grid-cells (boundary data
    // irrelevant).
    // - The Q-viscosity is non-zero only for compression, so we first
    // check div(v):
    //
    if (ctm==SimPM.tmOOA && c->isgd && CI.get_DivV(c) <0.0) {
      //if (CI.get_DivV(c) <0.0) {
      //
      // First the r-direction: take 2nd order central differences for
      // divergence and gradient.
      // The momenta get modified by the ith component of div(Q)
      //
      cp = grid->NextPt(c,XP);
      cn = grid->NextPt(c,XN);
      if (!cn || !cp) rep.error("No boundary data!",cn);
      //
      // SPHERICAL: NEED POSITIONS OF CELL-CENTRES.
      //
      rp = R_com(cp);
      rn = R_com(cn);
      QP = (CI.get_DivV(cp) <0.0) ? get_Qvisc(cp,XX) : 0.0;
      QN = (CI.get_DivV(cn) <0.0) ? get_Qvisc(cn,XX) : 0.0;
      //cout <<"QX: dU="<< c->dU[eqMMX] <<" corr="<< FV_dt*0.5*(QP-QN)/FV_dx<<"\n";
      //
      // Here we are taking a divergence, so we need the radial 
      // divergence equation, and there is a source term too.
      // Tscharnuter & Winkler (1979) eqn.47
      // SPHERICAL: *** DIFFERENCE FROM CARTESIAN! *************
      //
      c->dU[eqMMX] -= 
	//FV_dt*( (QP*pow(rp,3.0)-QN*pow(rn,3.0))*3.0/
	//	((pow(rp,3.0)-pow(rn,3.0))*R_com(c))
	//	);
	FV_dt*( (QP*pow(rp,3.0)-QN*pow(rn,3.0))/
		((rp-rn)*pow(R_com(c),3.0))
		);
      //
      // 2nd and 3rd directions are not implemented.
      //

      //
      // Now the energy correction:
      // We need partial derivative d(vx)/dx.
      // Use QN, QP for the two derivatives
      // We have already bugged out if the neighbour cells don't 
      // exist, so no need to check again.
      //
      prefactor = get_Qvisc(cp,Rsph);
      //cp = grid->NextPt(c,XP);
      //cn = grid->NextPt(c,XN);    
      QN = (cp->Ph[eqVX]-cn->Ph[eqVX])/(rp-rn);
#ifdef EINT_ETOT_PARALLEL
      c->dU[eqEINT] -= FV_dt*prefactor*QN;
#else
      //cout <<"QE: dU="<< c->dU[eqERG] <<" corr="<< FV_dt*prefactor*
      // 	(QN-CI.get_DivV(c)/3.0)*(QN-CI.get_DivV(c)/3.0) <<"\n";
      c->dU[eqERG]  -= FV_dt*prefactor*QN;
#endif // EINT_ETOT_PARALLEL
    } // if we have compression.


  } while ( (c =grid->NextPt(c)) !=0);
  
#endif // VNR_VISCOSITY

#ifdef FUNCTION_ID
  cout <<"sph_FV_solver_Hydro_Euler_Eint::PostProcess_dU ...returning.\n";
#endif //FUNCTION_ID
  return err;
}

void sph_FV_solver_Hydro_Euler_Eint::set_Qvisc(cell *c,
					   ///< cell to operate on.
					   const axes ax,
					   ///< axis we are looking along (Q1/2/3).
					   const double divv
					   ///< value of divv (should be pre-calculated).
					   )
{
  //
  // This is based on the cartesian version, but radial gradients must
  // be handled differently.
  //

  //
  // If div(v)>0 then Q=0, so we check that here:
  //
  if (divv>=0.0) {
    CI.set_Hcorr(c,ax,0.0);
    return;
  }

  //
  //  First add -1/3 div(v)
  //
  double value = -divv/3.0;

  //cout <<"divv="<<divv;
  //
  // Now we need to get a gradient.
  //
  cell *cn, *cp; int v=-1;
  switch (ax) {
  case Rsph:
    cn = grid->NextPt(c,XN);
    cp = grid->NextPt(c,XP);
    v = eqVX;
    break;
  default:
    rep.error("Bad axes",ax);
  }
  if (!cn && !cp)
    rep.error("no neighbours!",cp);
  if (!cn) cn=c;
  if (!cp) cp=c;

  //
  // Add velocity gradient along requested direction.
  //
  value += (cp->Ph[v]-cn->Ph[v])/(R_com(cp)-R_com(cn));

  //
  // Now multiply the lot by the prefactor.
  //
  value *=Qvisc_eta*Qvisc_eta*FV_dx*FV_dx*c->Ph[eqRO]*divv;
  //cout <<" prefactor="<<Qvisc_eta*Qvisc_eta*c->Ph[eqRO]*divv;
  //cout <<" ... total="<<value<<"\n";
  //
  // Finally add this to the Cell's extra data using the H-correction
  // interface.
  //
  CI.set_Hcorr(c,ax,value);
  return;
}


#endif // if INCLUDE_EINT_ADI_HYDRO
