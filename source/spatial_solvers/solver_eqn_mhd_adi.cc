///
/// \file solver_eqn_mhd_adi.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 written, based on old solver.cc classes.
///
/// Solver for the adiabatic Euler Equations.  Calculates flux via either Lax-Friedrichs
/// or Riemann solver (linear and/or exact).  Adds viscosity if asked for, and tracks flux
/// of N passive tracers.
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
/// - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
/// - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
/// - 2010.12.27 JM: Removed riemann_base references.  Added
///   function_id identifiers.
/// - 2010.12.30 JM: Added cell pointer to dU_cell()
/// - 2011.04.15 JM: Change in UtoP() for tracers (to try to correct for negative density!).
/// - 2013.02.07 JM: Tidied up for pion v.0.1 release.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.08.03 JM: Added pion_flt for double* arrays (allow floats)

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"

#include "solver_eqn_mhd_adi.h"
using namespace std;

// *****************************************
// ***** FV SOLVER MHD Ideal Adiabatic *****
// *****************************************


// ##################################################################
// ##################################################################

FV_solver_mhd_ideal_adi::FV_solver_mhd_ideal_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno, ///< CFL number
      const double cellsize, ///< dx, cell size.
      const double gam, ///< gas eos gamma.
      pion_flt *state, ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr ///< Number of tracer variables.
      )
  : eqns_base(nv),
    //riemann_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    flux_solver_mhd_ideal_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize)
{
#ifdef FUNCTION_ID
  cout <<"::FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  cout <<"::FV_solver_mhd_ideal_adi() constructor.\n";
  //cout <<"::FV_solver_mhd_ideal_adi() gamma = "<<eq_gamma<<"\n";
#endif

#ifdef FUNCTION_ID
  cout <<"::FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################

FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTING
  cout <<"FV_solver_mhd_ideal_adi::~FV_solver_mhd_ideal_adi() destructor.\n";
#endif

#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::PtoU(
      const pion_flt *p,
      pion_flt *u,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::PtoU ...starting.\n";
#endif //FUNCTION_ID

  eqns_mhd_ideal::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


int FV_solver_mhd_ideal_adi::UtoP(
      const pion_flt *u,
      pion_flt *p,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::UtoP ...starting.\n";
#endif //FUNCTION_ID

  int err=eqns_mhd_ideal::UtoP(u,p,g);
  // we use u[eqRO] because if there was a negative density, then usually
  // the tracer sign will follow, and this way we get a positive primitive
  // variable tracer back.  usually negative densities are unrecoverable, 
  // but you never know...
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/u[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::UtoP ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::PUtoFlux(
      const pion_flt *p,
      const pion_flt *u,
      pion_flt *f
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::PUtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_mhd_ideal::PUtoFlux(p,u,f);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = p[eqTR[t]]*f[eqRHO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::PUtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


void FV_solver_mhd_ideal_adi::UtoFlux(
      const pion_flt *u,
      pion_flt *f,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::UtoFlux ...starting.\n";
#endif //FUNCTION_ID

  eqns_mhd_ideal::UtoFlux(u,f,g);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = u[eqTR[t]]*f[eqRHO]/u[eqRHO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::UtoFlux ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################



///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_mhd_ideal_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c,          // Current cell.
        const axes d,     // Which axis we are looking along.
        const pion_flt *fn, // Negative direction flux.
        const pion_flt *fp, // Positive direction flux.
        const pion_flt *,   // slope vector for cell c.
        const int,        // spatial order of accuracy.
        const double,     // cell length dx.
        const double      // cell TimeStep, dt.
        )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::dU_Cell ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar];
  //
  // This calculates -dF/dx
  //
  //if (d!=eq_dir) rep.error("direction problem!!!!!!!!",d);
  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::dU_Cell ...returning.\n";
#endif //FUNCTION_ID
  return(err);
}


// ##################################################################
// ##################################################################

///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_mhd_ideal_adi::CellAdvanceTime(
      class cell *c, // cell to update.
      const pion_flt *Pin, // Initial State Vector.
      pion_flt *dU, // Update vector dU
      pion_flt *Pf, // Final state vector (can be same as initial vec.).
      pion_flt *dE, // Tracks change of energy if I have to correct for negative pressure
      const double, // gas EOS gamma.
      const double  // Cell timestep dt.
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellAdvanceTime ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar], u2[eq_nvar];

  //
  // First convert from Primitive to Conserved Variables
  //
  PtoU(Pin,u1, eq_gamma);

  //
  // Now add dU[] to U[], and change back to primitive variables.
  // This can give negative pressures, so check for that and fix it if needed.
  //
  for (int v=0;v<eq_nvar;v++) {
    //cout <<"du["<<v<<"] = "<<dU[v]<<"\n";
    u1[v] += dU[v];   // Update conserved variables
    dU[v] = 0.;       // Reset the dU array for the next timestep.
  }
  if(UtoP(u1,Pf, eq_gamma)!=0) {
    cout<<"(FV_solver_mhd_ideal_adi::CellAdvanceTime) UtoP complained (maybe about negative pressure...) fixing\n";
#ifdef TESTING
    //grid->PrintCell(dp.c);
    //grid->PrintCell(dp.c->ngb[XP]);
#endif //TESTING
    //rep.printVec("pin",Pin,SimPM.nvar);
    //PtoU(Pin,u1, eq_gamma);
    //rep.printVec("Uin",u1,SimPM.nvar);
    //rep.printVec("dU ",dU, SimPM.nvar);
    PtoU(Pf, u2, eq_gamma);
    *dE += (u2[ERG]-u1[ERG]);
    UtoP(u2,Pf, eq_gamma);
  }

#ifdef TESTING
  //  else if (dp.c->id==5887) {
  //    grid->PrintCell(dp.c->ngb[XN]);
  //    grid->PrintCell(dp.c);
  //    grid->PrintCell(dp.c->ngb[XP]);
  //    grid->PrintCell(dp.c->ngb[XP]->ngb[XP]);
  //  }
#endif //TESTING
 
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellAdvanceTime ...returning.\n";
#endif //FUNCTION_ID
  return 0;
}



// ##################################################################
// ##################################################################


///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_mhd_ideal_adi::CellTimeStep(
        const cell *c, ///< pointer to cell
        const double, ///< gas EOS gamma.
        const double  ///< Cell size dx.
        )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellTimeStep ...starting.\n";
#endif //FUNCTION_ID

  //
  // Get Max velocity along a grid direction.
  //
  pion_flt u1[eq_nvar];
  pion_flt temp = fabs(c->P[eqVX]);
  if (FV_gndim>1) temp = max(temp,fabs(c->P[eqVY]));
  if (FV_gndim>2) temp = max(temp,fabs(c->P[eqVZ]));
  
  /*  if (fabs(c->P[VX])>fabs(c->P[VY])) {temp = fabs(c->P[VX]);}
   else {temp = fabs(c->P[VY]);}
   if (fabs(c->P[VZ])>temp2) temp2=c->P[VZ];
   */
  
  //
  // First rotate the state vector to the fastest directions,
  // then add the fast speed to the max-velocity, and this is
  // the max wavespeed.
  //
  enum axes newdir;
  if (FV_gndim==1) temp += cfast(c->P,eq_gamma);
  else { // We may have to rotate the state vector to find the fastest fast speed.
    newdir = XX;
    if (fabs(c->P[BY])<fabs(c->P[BX])) {
      newdir =YY;
      if (fabs(c->P[BZ])<fabs(c->P[BY])) newdir=ZZ;
    }
    else if (fabs(c->P[BZ])<fabs(c->P[BX])) newdir=ZZ;

    if (newdir !=XX) {
      for (int v=0;v<eq_nvar;v++) u1[v] = c->P[v];
      rotate(u1,XX,newdir);
      temp += cfast(u1, eq_gamma);
    }
    else temp += cfast(c->P, eq_gamma);
  }
  FV_dt = FV_dx/temp;

  //
  // Check the gradient of pressure with neighbouring cells, since this can 
  // dramatically shorten the timestep.  (CAN'T REMEMBER WHY!!!)
  //
  //if (c->isgd) {
  //  int pg = static_cast<int>(eqPG);
  //  double grad = max_grad_abs(c,0,pg, grid)/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) FV_dt /= temp;
  //}
  //else if (grid->NextPt(c,XP)) {
  //  double grad = fabs(c->P[PG]-grid->NextPt(c,XP)->P[PG])/FV_dx/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) {
  //    FV_dt /= temp;
  //  }
  //}
  //else rep.error("No neighbour to gradient test",c);    

  
  //
  // Now scale the max. allowed timestep by the CFL number we are using (<1).
  //
  FV_dt *= FV_cfl;

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_ideal_adi::CellTimeStep ...returning.\n";
#endif //FUNCTION_ID
  return FV_dt;
}



// ##################################################################
// ##################################################################



// ********************************************************************************
// FV_solver_mhd_mixedGLM_adi class, for the Dedner-GLM divergence cleaning method.
// ********************************************************************************


// ##################################################################
// ##################################################################

FV_solver_mhd_mixedGLM_adi::FV_solver_mhd_mixedGLM_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno,   ///< CFL number
      const double cellsize,    ///< dx, cell size.
      const double gam,     ///< gas eos gamma.
      pion_flt *state,     ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr         ///< Number of tracer variables.
      )
  : eqns_base(nv),
    //riemann_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    flux_solver_mhd_ideal_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_mhd_mixedGLM(nv),
    flux_solver_mhd_mixedGLM_adi(nv,state,avcoeff,gam,ntr)
{
#ifdef FUNCTION_ID
  cout <<"::FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTNG
  cout <<"FV_solver_mhd_mixedGLM_adi::FV_solver_mhd_mixedGLM_adi() constructorg.\n";
  //cout <<"FV_solver_mhd_mixedGLM_adi::FV_solver_mhd_mixedGLM_adi() gamma = "<<eq_gamma<<"\n";
#endif

#ifdef FUNCTION_ID
  cout <<"::FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################

FV_solver_mhd_mixedGLM_adi::~FV_solver_mhd_mixedGLM_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef TESTNG
  cout <<"FV_solver_mhd_mixedGLM_adi::~FV_solver_mhd_mixedGLM_adi() destructor.\n";
#endif

#ifdef FUNCTION_ID
  cout <<"::~FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################

///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_mhd_mixedGLM_adi::CellAdvanceTime(
      class cell *c,
      const pion_flt *Pin, // Initial State Vector.
      pion_flt *dU, // Update vector dU
      pion_flt *Pf, // Final state vector (can be same as initial vec.).
      pion_flt *dE, // Tracks change of energy if I have to correct for negative pressure
      const double, // gas EOS gamma.
      const double // Cell timestep dt.
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...starting.\n";
#endif //FUNCTION_ID

  int err=FV_solver_mhd_ideal_adi::CellAdvanceTime(c,Pin,dU,Pf,dE,0,0);
  GLMsource(&(Pf[eqSI]),FV_dt);

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::CellAdvanceTime ...returning.\n";
#endif //FUNCTION_ID
  return(err);
}


// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::PtoU(
      const pion_flt *p,
      pion_flt *u,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::PtoU ...starting.\n";
#endif //FUNCTION_ID

  eqns_mhd_mixedGLM::PtoU(p,u,g);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::PtoU ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################

int FV_solver_mhd_mixedGLM_adi::UtoP(
      const pion_flt *u,
      pion_flt *p,
      const double g
      )
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::UtoP ...starting.\n";
#endif //FUNCTION_ID

  int err=eqns_mhd_mixedGLM::UtoP(u,p,g);
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/p[eqRO];

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::UtoP ...returning.\n";
#endif //FUNCTION_ID
  return err;
}


// ##################################################################
// ##################################################################

void FV_solver_mhd_mixedGLM_adi::GotTimestep(
        const double delt ///< timestep dt.
	)
{
#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep ...starting.\n";
#endif //FUNCTION_ID

  //     cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep() setting wave speeds.\n";
  GLMsetPsiSpeed(FV_cfl,FV_dx,delt);

#ifdef FUNCTION_ID
  cout <<"FV_solver_mhd_mixedGLM_adi::GotTimestep ...returning.\n";
#endif //FUNCTION_ID
}



// ##################################################################
// ##################################################################




/// ---------------------------------------------------------------------
/// -------------------  AXI-SYMMETRIC EQUATIONS ------------------------
/// ---------------------------------------------------------------------




// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_ideal_adi::cyl_FV_solver_mhd_ideal_adi(
      const int nv, ///< number of variables in state vector.
      const int nd, ///< number of space dimensions in grid.
      const double cflno, ///< CFL number
      const double cellsize, ///< dx, cell size.
      const double gam, ///< gas eos gamma.
      pion_flt *state,   ///< State vector of mean values for simulation.
      const double avcoeff, ///< Artificial Viscosity Parameter etav.
      const int ntr    ///< Number of tracer variables.
      )
  : eqns_base(nv),
    // riemann_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam),
    Riemann_Roe_MHD_CV(nv,gam),
    flux_solver_mhd_ideal_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

//  cout <<"cyl_FV_solver_mhd_ideal_adi CONSTRUCTOR\n";
//  cout <<"IdealMHD Equations; Riemann Solver Method; Cylindrical Coordinates.\n";
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for 2d axial symmetry \
			 so far.  Sort it out!",nd);

#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}


// ##################################################################
// ##################################################################


cyl_FV_solver_mhd_ideal_adi::~cyl_FV_solver_mhd_ideal_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_ideal_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_ideal_adi ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################

int cyl_FV_solver_mhd_ideal_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c, ///< Current cell.
        const axes d, ///< Which axis we are looking along.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        const pion_flt *dpdx, ///< slope vector for cell c.
        const int OA,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        )
{
#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_ideal_adi::dU_Cell ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar];

  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  if (d==Rcyl) {
    double pm = (c->Ph[eqBX]*c->Ph[eqBX] +c->Ph[eqBY]*c->Ph[eqBY] +c->Ph[eqBZ]*c->Ph[eqBZ])/2.;
    switch (OA) {
     case OA1:
       c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm)/CI.get_dpos(c,Rcyl);
      break;
     case OA2:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm + 
			  (CI.get_dpos(c,Rcyl)-R_com(c))*
		   (dpdx[eqPG] +c->Ph[eqBX]*dpdx[eqBX] +c->Ph[eqBY]*dpdx[eqBY] +c->Ph[eqBZ]*dpdx[eqBZ]))
                   /CI.get_dpos(c,Rcyl);
      break;
     default:
      rep.error("Bad OOA in cyl_IdealMHD_RS::dU, only know 1st,2nd",OA);
    }
  }

#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_ideal_adi::dU_Cell ...returning.\n";
#endif //FUNCTION_ID
  return(err);
} 


// ##################################################################
// ##################################################################

//------------------------------------------------------------------//
cyl_FV_solver_mhd_mixedGLM_adi::cyl_FV_solver_mhd_mixedGLM_adi(
        const int nv, ///< number of variables in state vector.
        const int nd, ///< number of space dimensions in grid.
        const double cflno, ///< CFL number
        const double cellsize, ///< dx, cell size.
        const double gam, ///< gas eos gamma.
        pion_flt *state,   ///< State vector of mean values for simulation.
        const double avcoeff, ///< Artificial Viscosity Parameter etav.
        const int ntr    ///< Number of tracer variables.
        )
  : eqns_base(nv),
    flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,gam,avcoeff,ntr),
    eqns_mhd_ideal(nv),
    riemann_MHD(nv,state,gam), 
    Riemann_Roe_MHD_CV(nv,gam),
    flux_solver_mhd_ideal_adi(nv,state,avcoeff,gam,ntr),
    VectorOps_Cart(nd,cellsize),
    FV_solver_mhd_ideal_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    eqns_mhd_mixedGLM(nv),
    flux_solver_mhd_mixedGLM_adi(nv,state,avcoeff,gam,ntr),
    FV_solver_mhd_mixedGLM_adi(nv,nd,cflno,cellsize,gam,state,avcoeff,ntr),
    VectorOps_Cyl(nd,cellsize)
{
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

//  cout <<"cyl_FV_solver_mhd_mixedGLM_adi CONSTRUCTOR\n";
//  cout <<"glmMHD Equations; Riemann Solver Method; Cylindrical Coordinates.\n";
  if (nd!=2) rep.error("Cylindrical coordinates only implemented for \
                        2d axial symmetry so far.  Sort it out!",nd);
#ifdef FUNCTION_ID
  cout <<"::cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
  return;
}



// ##################################################################
// ##################################################################

cyl_FV_solver_mhd_mixedGLM_adi::~cyl_FV_solver_mhd_mixedGLM_adi()
{
#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_mixedGLM_adi ...starting.\n";
#endif //FUNCTION_ID

#ifdef FUNCTION_ID
  cout <<"::~cyl_FV_solver_mhd_mixedGLM_adi ...returning.\n";
#endif //FUNCTION_ID
}


// ##################################################################
// ##################################################################

int cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell(
        class GridBaseClass *grid,
        cell *c, ///< Current cell.
        const axes d, ///< Which axis we are looking along.
        const pion_flt *fn, ///< Negative direction flux.
        const pion_flt *fp, ///< Positive direction flux.
        const pion_flt *dpdx, ///< slope vector for cell c.
        const int OA,      ///< spatial order of accuracy.
        const double, ///< cell length dx.
        const double  ///< cell TimeStep, dt.
        )
{
#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell ...starting.\n";
#endif //FUNCTION_ID

  pion_flt u1[eq_nvar];
  int err = DivStateVectorComponent(c, grid, d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  if (d==Rcyl) {
    double pm = (c->Ph[eqBX]*c->Ph[eqBX] +c->Ph[eqBY]*c->Ph[eqBY] +c->Ph[eqBZ]*c->Ph[eqBZ])/2.;
    switch (OA) {
     case OA1:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm)/CI.get_dpos(c,Rcyl);
      c->dU[eqBBX] += FV_dt*c->Ph[eqSI]/CI.get_dpos(c,Rcyl);
      break;
     case OA2:
      c->dU[eqMMX] += FV_dt*(c->Ph[eqPG]+pm + 
			 (CI.get_dpos(c,Rcyl)-R_com(c))*
			 (dpdx[eqPG] +c->Ph[eqBX]*dpdx[eqBX] +c->Ph[eqBY]*dpdx[eqBY] +c->Ph[eqBZ]*dpdx[eqBZ]))
	                /CI.get_dpos(c,Rcyl);
      c->dU[eqBBX] += FV_dt*(c->Ph[eqSI] +(CI.get_dpos(c,Rcyl)-R_com(c))*dpdx[eqSI])/CI.get_dpos(c,Rcyl);
      break;
     default:
      rep.error("Bad OOA in cyl_glmMHD_RS::dU, only know 1st,2nd",OA);
    }
  }

#ifdef FUNCTION_ID
  cout <<"cyl_FV_solver_mhd_mixedGLM_adi::dU_Cell ...returning.\n";
#endif //FUNCTION_ID
  return(err);
} 


// ##################################################################
// ##################################################################

