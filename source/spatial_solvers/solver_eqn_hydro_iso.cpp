///
/// \file solver_eqn_hydro_iso.cc
/// \author Jonathan Mackey
/// History: 2009-10-21 Created file, partly from bits of solver.cc classes.
///
/// Solver for the isothermal Euler Equations.  Calculates flux via either Lax-Friedrichs
/// or Riemann solver (linear and/or exact).  Adds viscosity if asked for, and tracks flux
/// of N passive tracers.
///
/// - 2010-02-10 JM: Forced the flux calculation to use the Roe solver, with a
///   warning when the flag is set to a different solver.  This is because the
///   other solvers are all useless at the moment.
///
/// - 2010-07-20 JM: changed order of accuracy variables to integers.
///
///  - 2010.09.30 JM: Worked on Lapidus AV (added Cl,Cr pointers to flux functions).
///
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///   Made InterCellFlux general for all classes (moved to FV_solver_base)
///
/// - 2010.12.27 JM: Put all isothermal dynamics in an ifdef b/c I
///   updated the code structure which has broken everything and I
///   don't have time to fix isothermal stuff now...
///

#ifdef ISOTHERMAL_SOLVERS_ENABLED

#include "solver_eqn_hydro_iso.h"
using namespace std;

// **************************************
// ***** FV SOLVER HYDRO ISOTHERMAL *****
// **************************************
FV_solver_Hydro_iso::FV_solver_Hydro_iso(const int nv, ///< number of variables in state vector.
					 const int nd, ///< number of space dimensions in grid.
					 const double cflno,   ///< CFL number
					 const double cellsize,    ///< dx, cell size.
					 const double,     ///< gas eos gamma.
					 double *state,     ///< State vector of mean values for simulation.
					 const double avcoeff, ///< Artificial Viscosity Parameter etav.
					 const int ntr         ///< Number of tracer variables.
					 )
  : eqns_base(nv), riemann_base(nv), flux_solver_base(nv,avcoeff,ntr),
    FV_solver_base(nv,nd,cflno,cellsize,0,avcoeff,ntr),
    eqns_IsoEuler(nv), riemann_hydro_iso(nv,state), flux_solver_hydro_iso(nv,state,avcoeff,ntr),
    VectorOps_Cart(nd,cellsize)
{
  cout <<"FV_solver_Hydro_iso::FV_solver_Hydro_iso() constructor which does nothing.\n";
  return;
}

FV_solver_Hydro_iso::~FV_solver_Hydro_iso()
{
  cout <<"FV_solver_Hydro_iso::~FV_solver_Hydro_iso() destructor which does nothing.\n";
  return;
}


void FV_solver_Hydro_iso::PtoU(const double* p, double* u, const double)
{
  eqns_IsoEuler::PtoU(p,u,0);
  for (int t=0;t<FS_ntr;t++) u[eqTR[t]] = p[eqTR[t]]*p[eqRO];
  return;
}

int FV_solver_Hydro_iso::UtoP(const double *u, double *p, const double)
{
  int err=eqns_IsoEuler::UtoP(u,p,0);
  for (int t=0;t<FS_ntr;t++) p[eqTR[t]] = u[eqTR[t]]/p[eqRO];
  return err;
}

void FV_solver_Hydro_iso::PUtoFlux(const double *p, const double *u, double *f)
{
  eqns_IsoEuler::PUtoFlux(p,u,f);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = p[eqTR[t]]*f[eqRHO];
  return;
}

void FV_solver_Hydro_iso::UtoFlux(const double *u, double *f, const double)
{
  eqns_IsoEuler::UtoFlux(u,f,0);
  for (int t=0;t<FS_ntr;t++) f[eqTR[t]] = u[eqTR[t]]*f[eqRHO]/u[eqRHO];
  return;
}


///
/// Adds the contribution from flux in the current direction to dU.
///
int FV_solver_Hydro_iso::dU_Cell(cell *c,          // Current cell.
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
  int err = DivStateVectorComponent(c,d,eq_nvar,fn,fp,u1);
  for (int v=0;v<eq_nvar;v++) c->dU[v] += FV_dt*u1[v];
  
  return(err);
}

///
/// General Finite volume scheme for updating a cell's
/// primitive state vector, for homogeneous equations.
///
int FV_solver_Hydro_iso::CellAdvanceTime(const double *Pin, // Initial State Vector.
					 double *dU, // Update vector dU
					 double *Pf, // Final state vector (can be same as initial vec.).
					 double *dE, // Tracks change of energy if I have to correct for negative pressure
					 const double, // gas EOS gamma.
					 const double  // Cell timestep dt.
					 )
{
  double u1[eq_nvar], u2[eq_nvar];

  //
  // First convert from Primitive to Conserved Variables
  //
  PtoU(Pin,u1,0);
#ifdef TESTING
//  if (dp.c->id==984) {
//    cout <<"____________________DU()DU()DU()DU()____id="<<dp.c->id<<"________________\n";
//    CI.print_cell(dp.c);
//    rep.printVec("U0",u1,eq_nvar);
//    cout <<"____________________DU()DU()DU()DU()____id="<<dp.c->id<<"________________\n";
//  }
#endif

  //
  // Now add dU[] to U[], and change back to primitive variables.
  // This can give negative pressures, so check for that and fix it if needed.
  //
  for (int v=0;v<eq_nvar;v++) {
    u1[v] += dU[v];   // Update conserved variables
    dU[v] = 0.;       // Reset the dU array for the next timestep.
  }
  if(UtoP(u1,Pf, 0)!=0) {
    cout<<"(LF_FVSolver::CellAdvanceTime) UtoP complained (maybe about negative pressure...) fixing\n";
#ifdef TESTING
    //grid->PrintCell(dp.c);
    //grid->PrintCell(dp.c->ngb[XP]);
#endif //TESTING
    //rep.printVec("pin",Pin,SimPM.nvar);
    //PtoU(Pin,u1,0);
    //rep.printVec("Uin",u1,SimPM.nvar);
    //rep.printVec("dU ",dU, SimPM.nvar);
    PtoU(Pf, u2, 0);
    *dE += (u2[ERG]-u1[ERG]);
    UtoP(u2,Pf, 0);
  }

#ifdef TESTING
  //  else if (dp.c->id==5887) {
  //    grid->PrintCell(dp.c->ngb[XN]);
  //    grid->PrintCell(dp.c);
  //    grid->PrintCell(dp.c->ngb[XP]);
  //    grid->PrintCell(dp.c->ngb[XP]->ngb[XP]);
  //  }
#endif //TESTING

  return 0;
}

///
/// Given a cell, calculate the hydrodynamic timestep.
///
double FV_solver_Hydro_iso::CellTimeStep(const cell *c, ///< pointer to cell
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
  //cout <<"hello jonathan from isothermal hydro solver\n";

  //
  // Get Max velocity along a grid direction.
  //
  double temp = fabs(c->P[eqVX]);
  if (FV_gndim>1) temp = max(temp,fabs(c->P[eqVY]));
  if (FV_gndim>2) temp = max(temp,fabs(c->P[eqVZ]));
  
  //
  // Add the sound speed to this to get the max wavespeed.
  //
  temp += c->P[eqAA];
  FV_dt = FV_dx/temp;

//   if (c->id==1062) {
//     cout <<"FV_dt="<<FV_dt<<"\tc_hydro="<<temp<<"\n";
//     CI.print_cell(c);
//     CI.print_cell(grid->NextPt(c,XP));
//     CI.print_cell(grid->NextPt(c,XN));
//     CI.print_cell(grid->NextPt(c,YP));
//     CI.print_cell(grid->NextPt(c,YN));
//   }
  
  //
  // Check the gradient of pressure with neighbouring cells, since this can 
  // dramatically shorten the timestep.
  //
  //int aa = static_cast<int>(eqAA);
  //if (c->isgd) {
  //  double grad = maxGradAbs(c,0,aa)/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) {
  //    FV_dt /= temp;
  //    //      if (temp>3) cout <<"gradient adjusting timestep\n";
  //  }
  //}
  //else if (grid->NextPt(c,XP)) {
  //  //  cout <<"grid="<<grid<<"\n";
  //  double grad = fabs(c->P[eqAA]-grid->NextPt(c,XP)->P[eqAA])/FV_dx/c->P[RO];
  //  if( (temp = grad*FV_dt/temp) >1.) {
  //    FV_dt /= temp;
  //  }
  //}
  //  else rep.error("No neighbour to gradient test",c);    
  
  //
  // Now scale the max. allowed timestep by the CFL number we are using (<1).
  FV_dt *= FV_cfl;
  return FV_dt;
}


#endif // ISOTHERMAL_SOLVERS_ENABLED
