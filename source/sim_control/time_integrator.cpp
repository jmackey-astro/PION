/// \file time_integrator.cpp
/// \brief time integration routines for sim_control class.
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for sim_control 
/// class, which is a 1st/2nd order Finite Volume Solver following
/// Falle, Komissarov, \& Joarder (1998), MNRAS, 297, 265.
/// 
/// Modifications:
/// - 2018.01.24 JM: worked on making SimPM non-global
/// - 2018.05.10 JM: moved calc_timestep function to its own class.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"

#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "sim_control/time_integrator.h"
#include "dataIO/dataio_base.h"

#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include "spatial_solvers/solver_eqn_base.h"


#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <time.h>
using namespace std;

// ##################################################################
// ##################################################################



time_integrator::time_integrator()
{
  return;
}



// ##################################################################
// ##################################################################



time_integrator::~time_integrator()
{
  return;
}



// ##################################################################
// ##################################################################

double time_integrator::advance_time(
      const int level,          ///< level in grid hierarchy
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // Check order-of-accuracy (OOA) requested, and perform the
  // appropriate update.
  //
  if      (SimPM.tmOOA==OA1 && SimPM.spOOA==OA1 ) {
    //
    // Send in full timestep, and the order-of-accuracy, so the 
    // function knows whether to update P[i] as well as Ph[i]
    //
    //cout <<"First order update\n";
    err += first_order_update(SimPM.dt, SimPM.tmOOA, grid);
    if (err)
      rep.error("first_order_update() returned error",err);

    // Update boundary data.
    err += TimeUpdateInternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA1, OA1);
    err += TimeUpdateExternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA1, OA1);
    if (err) 
      rep.error("second_order_update: error from bounday update",err);
  }

  else if (SimPM.tmOOA==OA2 && SimPM.spOOA==OA2) {
    //cout <<"Second order update\n";
    err += first_order_update( 0.5*SimPM.dt, OA2, grid);
    if (err)
      rep.error("1st order time-update returned error",err);
    // Update boundary data.
    err += TimeUpdateInternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA1, OA2);
    err += TimeUpdateExternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA1, OA2);
    if (err) 
      rep.error("second_order_update: error from bounday update",err);

    err += second_order_update(SimPM.dt,     OA2, grid);
    if (err)
      rep.error("Second order time-update returned error",err);

    // Update boundary data.
    err += TimeUpdateInternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA2, OA2);
    err += TimeUpdateExternalBCs(SimPM, level, grid, spatial_solver,
                                        SimPM.simtime,   OA2, OA2);
    if (err) 
      rep.error("second_order_update: error from bounday update",err);
  }
  //
  // Add in 3rd order PPM at some stage???
  //
  else {
    rep.error("Bad OOA requests; choose (1,1) or (2,2)",SimPM.tmOOA);
  }

  //
  // Update timestepping variables to new state.
  //
  //  cout <<"now dt = "<<SimPM.dt<<"\n";
  SimPM.simtime +=SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;

  return SimPM.dt;
}



// ##################################################################
// ##################################################################



int time_integrator::first_order_update(
      const double dt,
      const int   ooa,
      class GridBaseClass *grid ///< Computational grid.
      )
{
  // NB Only used for uniform grid.  update for NG grid is in
  // sim_control_NG.cpp
  int err=0;
  //
  // Set dt for equations class
  //
  spatial_solver->Setdt(dt);

  //
  // May need to do raytracing, if it wasn't needed for calculating
  // the timestep.
  //
  if (!FVI_need_column_densities_4dt && grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid,0);
    if (err) 
      rep.error("first_order_update: error from first calc_rt_cols()",err);
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(dt, grid);
  err += calc_dynamics_dU(dt,TIMESTEP_FIRST_PART, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt,TIMESTEP_FIRST_PART, grid);
#endif // THERMAL_CONDUCTION
  if (err) 
    rep.error("first_order_update: error from calc_*_dU",err);
  
  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(dt,TIMESTEP_FIRST_PART,ooa, grid);
  if (err) 
    rep.error("first_order_update: error from state-vec update",err);

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::second_order_update(
      const double dt,
      const int   ooa,
      class GridBaseClass *grid ///< Computational grid.
      )
{
  // NB Only used for uniform grid.  update for NG grid is in
  // sim_control_NG.cpp
  int err=0;
  //
  // Set dt for equations class
  // MULTITHREADING RISK!
  //
  spatial_solver->Setdt(dt);

  //
  // Raytracing, to get column densities for microphysics update.
  //
  if (grid->RT) {
    err += calculate_raytracing_column_densities(SimPM,grid,0);
    if (err) {
      rep.error("second_order_update: error from first calc_rt_cols()",err);
    }
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(      dt,      grid);
  err += calc_dynamics_dU(          dt, TIMESTEP_FULL, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt, TIMESTEP_FULL, grid);
#endif // THERMAL_CONDUCTION
  if (err) 
    rep.error("second_order_update: error from calc_*_dU",err);
  
  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(  dt, TIMESTEP_FULL, ooa, grid);
  if (err) 
    rep.error("second_order_update: error from state-vec update",err);

  return 0;
}




// ##################################################################
// ##################################################################



int time_integrator::calc_microphysics_dU(
      const double delt, ///< timestep to integrate MP eqns.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //cout <<"\tcalc_microphysics_dU starting.\n";

  //
  // If we are not doing microphysics or raytracing, the modules will
  // not be initialised, so just return zero.  (RT requires MP)
  //
  if (!MP) return 0;
  
#ifdef TESTING
  cout <<"calc_microphysics_dU() Updating MicroPhysics. ";
  cout <<"  RT-Nsrc="<<SimPM.RS.Nsources<<"\n";
#endif // TESTING
  int err = 0;

  if (SimPM.RS.Nsources==0) {
    //
    // If no radiative transfer, then just do a simple MP update.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling calc_noRT_microphysics_dU()\n";
#endif // RT_TESTING
    err += calc_noRT_microphysics_dU(delt, grid);
  }

  else {
    //
    // must have at least one radiation source, so we call a newer
    // update function with a bit more overhead.
    //
#ifdef RT_TESTING
    cout <<"\t\t--- calling calc_RT_microphysics_dU()\n";
#endif // RT_TESTING
    err += calc_RT_microphysics_dU(delt, grid);
  }
    
  //cout <<"\tcalc_microphysics_dU finished.\n";
  return err;
}


// ##################################################################
// ##################################################################



int time_integrator::calc_RT_microphysics_dU(
      const double delt, // timestep to integrate
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef RT_TESTING
  if (!grid->RT)
    rep.error("Logic error: must have RT unless i'm an idiot",
              "GENERAL-RT");
  if (SimPM.RS.Nsources<1)
    rep.error("Need at least one source for diffuse-RT update",
              SimPM.RS.Nsources);
#endif // RT_TESTING

  int err=0;
  //
  // Do MP update.  Only new MP classes have the  XX_RTnew() defined,
  // so if we try to call this with old code then it should return with 
  // a non-zero error code.
  //
  cell *c = grid->FirstPt();
  pion_flt p[SimPM.nvar]; // temporary state vector for output state.
  pion_flt ui[SimPM.nvar], uf[SimPM.nvar]; // conserved variable states.

  double tt=0.; // temperature returned at end of microphysics step.
  do {
    //
    // Check if cell is internal boundary data or not.  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (!c->isdomain) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in calc_RT_microphysics_dU() c->isdomain.\n";
#endif
    }
    else {
      //
      // Get column densities and Vshell in struct for each source.
      //
      for (int v=0; v<FVI_nheat; v++) {
        FVI_heating_srcs[v].Vshell  = CI.get_cell_Vshell(c, FVI_heating_srcs[v].id);
        FVI_heating_srcs[v].dS      = CI.get_cell_deltaS(c, FVI_heating_srcs[v].id);
        CI.get_cell_col(c, FVI_heating_srcs[v].id, FVI_heating_srcs[v].DelCol);
        CI.get_col(     c, FVI_heating_srcs[v].id, FVI_heating_srcs[v].Column);
        for (short unsigned int iC=0; iC<FVI_heating_srcs[v].NTau; iC++)
          FVI_heating_srcs[v].Column[iC] -= FVI_heating_srcs[v].DelCol[iC];
      }
      for (int v=0; v<FVI_nion; v++) {
        FVI_ionising_srcs[v].Vshell = CI.get_cell_Vshell(c, FVI_ionising_srcs[v].id);
        FVI_ionising_srcs[v].dS     = CI.get_cell_deltaS(c, FVI_ionising_srcs[v].id);
        CI.get_cell_col(c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].DelCol);
        CI.get_col(     c, FVI_ionising_srcs[v].id, FVI_ionising_srcs[v].Column);
        for (short unsigned int iC=0; iC<FVI_ionising_srcs[v].NTau; iC++)
          FVI_ionising_srcs[v].Column[iC] -= FVI_ionising_srcs[v].DelCol[iC];
      }
      //
      // 4th and 5th args are for ionising sources.
      //
      err += MP->TimeUpdateMP_RTnew(c->P, FVI_nheat,FVI_heating_srcs,
                                    FVI_nion, FVI_ionising_srcs,
                                    p, delt, SimPM.gamma, 0, &tt);

//#define NO_COOLING_ON_AXIS
#ifdef NO_COOLING_ON_AXIS
      //cout <<"hello\n";
//#error "Fix HACK in time_integrator.cpp"
      if (SimPM.coord_sys==COORD_CYL && 
          //!grid->NextPt(c,YN)->isgd &&
          p[RO] > 1.4e-20 &&    // density more than twice background density
          //p[RO] > 0.77e-20 &&    // density more than 1.1x background density
          c->pos[Rcyl] < 5 &&   // get the first three radial cells (R=0,2,4)
          c->pos[Zcyl] > 0 &&   // only consider cells with z>0 (upstream)
          p[PG] < c->P[PG] &&   // only consider cells that were cooled (not heated)
          p[SimPM.ftr] < 0.5    // only consider mostly neutral gas
          ) {
        //tt = MP->Temperature(p,SimPM.gamma);
        //if (tt < 1.0e3) {	
        //  MP->Set_Temp(p,1.0e3,SimPM.gamma);
        //}
        // Just set the pressure equal to what it was before cooling,
        // so that the gas is adiabatic.
        p[PG] = c->P[PG];
      }
#endif // NO_COOLING_ON_AXIS

      //
      // New state is p[], old state is c->P[].  Get dU from these.
      //
      spatial_solver->PtoU(c->P,ui,SimPM.gamma);
      spatial_solver->PtoU(p,   uf,SimPM.gamma);
      for (int v=0;v<SimPM.nvar;v++) c->dU[v] += uf[v]-ui[v];

    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"calc_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
} // RT microphysics update.



// ##################################################################
// ##################################################################



int time_integrator::calc_noRT_microphysics_dU(
      const double delt, ///< timestep to integrate
      class GridBaseClass *grid ///< Computational grid.
      )
{
#ifdef TESTING
  cout <<"calc_noRT_microphysics_dU starting.\n";
#endif
  //
  // No radiation sources and no diffuse radiation optical depths,
  // so call a simple microphysics update.
  //
  cell *c = grid->FirstPt();
  pion_flt p[SimPM.nvar]; // temporary state vector for output state.
  pion_flt ui[SimPM.nvar], uf[SimPM.nvar]; // conserved variable states.
  double tt=0.; // temperature returned at end of microphysics step.
  int err=0;
  do {
    //
    // Check if cell is internal boundary data or not.  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (!c->isdomain) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in calc_noRT_microphysics_dU() c->isdomain.\n";
#endif
    }
    else {
      //
      // integer 5th argument determines type of integration substepping:
      // 0 = adaptive RK5 Cash-Karp method.
      // 1 = adaptive euler integration.
      // 2 = single step RK4 method (at your own risk!)
      //
      err += MP->TimeUpdateMP(c->P, p, delt, SimPM.gamma, 0, &tt);
      //rep.printVec("Original vector P",c->P,nvar);
      //rep.printVec("Updated  vector p",p   ,nvar);
      if (err)
        rep.error("calc_noRT_microphysics_dU returned error: cell id",c->id);

      //
      // New state is p[], old state is c->P[].  Get dU from these.
      //
      spatial_solver->PtoU(c->P,ui,SimPM.gamma);
      spatial_solver->PtoU(p,   uf,SimPM.gamma);
      for (int v=0;v<SimPM.nvar;v++) c->dU[v] += uf[v]-ui[v];

    } // if not boundary data.
  } while ( (c=grid->NextPt(c)) !=0);
  //    cout <<"calc_noRT_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
} 



// ##################################################################
// ##################################################################


  
int time_integrator::calc_dynamics_dU(
      const double dt, ///< timestep to integrate
      const int step, ///< whether TIMESTEP_FIRST_PART or TIMESTEP_FULL.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  //cout <<"\tcalc_dynamics_dU starting.\n";
  //
  // first check if we are doing dynamics, and return if not.
  //
  if (!SimPM.EP.dynamics) return 0;
  int err=0;

#ifdef TESTING
  if (step==TIMESTEP_FIRST_PART)
    cout <<"*****Updating dynamics: OA1\n";
  else if (step==TIMESTEP_FULL)
    cout <<"*****Updating dynamics: OA2\n";
  else rep.error("Bad ooa",step);
#endif //TESTING

  //
  // First we pre-process the cells, if needed.  This is required for
  // multi-dimensional viscosity such as the H-Correction.
  //
  err = spatial_solver->preprocess_data(step, SimPM, grid);

  //
  // Now calculate the directionally-unsplit time update for the
  // conserved variables:
  //
  err = set_dynamics_dU(dt, step, grid); //,time_ooa);
  rep.errorTest("calc_dynamics_dU() set_dynamics_dU returned error.", 0,err);

  //
  // This function is used for refined grids, to make the flux across
  // grid boundaries be consistent across all levels.
  //
  // Other potential uses of post-processing include if we are doing
  // something like Constrained Transport, where we have to change
  // the B-field update.
  //
  err = spatial_solver->PostProcess_dU(dt, step, SimPM, grid);;
  rep.errorTest("calc_dynamics_dU() spatial_solver->PostProcess_dU()",0,err);

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::set_dynamics_dU(
      const double dt,     ///< timestep for this calculation
      const int step, ///< whether half-step or full-step
      class GridBaseClass *grid ///< Computational grid.
      )
{
  // 
  // Allocate arrays
  //
  int return_value=0;
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP; posdirs[1] = YP; posdirs[2] = ZP;
  negdirs[0] = XN; negdirs[1] = YN; negdirs[2] = ZN;
  axis[0] = XX; axis[1] = YY; axis[2] = ZZ;

  int space_ooa;
  if (step == TIMESTEP_FIRST_PART)
    space_ooa=OA1;
  else
    space_ooa=OA2;

#ifdef TEST_INT
  // get current level of grid in hierarchy.
  int level=0;
  if (SimPM.grid_nlevels >1) {
    for (int v=0;v<SimPM.grid_nlevels;v++) {
      if (grid == SimPM.levels[v].grid) level = v;
    }
  }
  cout <<"*** Calculating DU dynamics on level "<<level<<".\n";
#endif

  //
  // Loop over all directions, and in each direction, calculate fluxes
  // in all columns of cells in that direction (it does work!).
  // This function depends on cells being labelled as on-grid or as
  // boundary cells, and also on dynamics_dU_column returning 0 on successful
  // completion, and -1 if the column ends up at the last cell in the domain.
  // Any other return value will signal an error in this function, stopping the code.
  //
  // 2011.04.29 JM: changed logic here so we check for cells which are not grid
  // cells.  Checking for boundary data is not correct, since we can have 
  // internal boundaries which are also grid data.
  //
  for (int i=0;i<SimPM.ndim;i++) {
    spatial_solver->SetDirection(axis[i]);
    class cell *cpt    = grid->FirstPt_All();
    class cell *marker = cpt;

#ifdef TEST_INT
    cout <<"Direction="<<axis[i]<<", i="<<i<<"\n";
    rep.printVec("cpt",cpt->pos,SimPM.ndim);
#endif
    
    //
    // loop over the number of cells in the line/plane of starting
    // cells.
    //
    enum direction d1 = posdirs[(i+1)%3];
    enum direction d2 = posdirs[(i+2)%3];
    enum axes x1 = axis[(i+1)%3];
    enum axes x2 = axis[(i+2)%3];

    //
    // loop over the two perpendicular axes, to trace out a plane of
    // starting cells for calculating fluxes along columns along this
    // axis.  Note that the NG_All() array is initialised so that
    // unused dimensions have NG=1, so the plane can be a single cell
    // (in 1D) or a line (in 2D) or a plane (in 3D).
    //
    for (int ax2=0; ax2<grid->NG_All(x2); ax2++) {
      for (int ax1=0; ax1<grid->NG_All(x1); ax1++) {
#ifdef TEST_INT
        cout <<"ax1="<<ax1<<", ax2="<<ax2<<", i="<<i;
        cout<<", cpt="<<cpt<<":  ";
        //CI.print_cell(cpt);
        cout <<"\n";
#endif
        return_value = dynamics_dU_column(cpt,posdirs[i],negdirs[i], dt,
                                        space_ooa, grid);
        rep.errorTest("set_dynamics_dU: column",0,return_value);
        cpt = grid->NextPt(cpt,d1);
      }
      marker = grid->NextPt(marker,d2);
      cpt = marker;
    } // loop over all columns

  } // Loop over three directions.
  spatial_solver->SetDirection(axis[0]); // Reset fluxes to x-dir, (just to be safe!).
  return 0;
}   // set_dynamics_dU()



// ##################################################################
// ##################################################################

  

int time_integrator::dynamics_dU_column(
      class cell *startingPt, ///< starting point of column.
      const enum direction posdir, ///< direction to trace column.
      const enum direction negdir, ///< reverse direction
      const double dt, ///< timestep we are advancing by.
      const int csp,  ///< spatial order-of-accuracy for this step.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  if ( (SimPM.spOOA>2) || (SimPM.tmOOA>2) || (csp>2)  ) {
    cerr <<"(dynamics_dU_column) Error, ooa="<< SimPM.spOOA <<", ";
    cerr << SimPM.tmOOA <<".\n";
    return(1);
  }
  //cout <<"dynamics_dU_column: d+="<<posdir<<", d-="<<negdir;
  //cout <<", csp="<<csp<<", OOA="<<SimPM.spOOA<<"\n";
  int err = 0;
  enum axes axis = spatial_solver->GetDirection();
  double dx = grid->DX();
#ifdef TEST_CONSERVATION 
  double dE=0.0, dMX=0.0, dMY=0.0, dMZ=0.0, dM=0.0;
#endif

  // Calculate Flux at positive (right) boundary of cell for the
  // current cell (Fr_this) and the negative neighbour (Fr_prev).
  pion_flt *Fr_this=0, *Fr_prev=0, *temp=0, *slope_cpt=0,
           *slope_npt=0, *edgeR=0, *edgeL=0, *pstar=0;
  Fr_prev   = mem.myalloc(Fr_prev,   SimPM.nvar);
  Fr_this   = mem.myalloc(Fr_this,   SimPM.nvar);
  slope_cpt = mem.myalloc(slope_cpt, SimPM.nvar);
  slope_npt = mem.myalloc(slope_npt, SimPM.nvar);
  edgeL     = mem.myalloc(edgeL,     SimPM.nvar);
  edgeR     = mem.myalloc(edgeR,     SimPM.nvar);
  pstar     = mem.myalloc(pstar,     SimPM.nvar);

  //
  // Set starting point, and next two points in the column.
  //
  cell *cpt = startingPt; 
  if(cpt==0) {
    cerr<<"(dynamics_dU_column) error finding left boundary cell.\n";
    return(1);
  }
  cell *npt  = grid->NextPt(cpt,posdir);
  cell *n2pt = grid->NextPt(npt,posdir);
#ifdef TEST_INT
  cout <<"Column: "<<cpt<<", "<<npt<<", "<<n2pt<<"\n";
#endif
  if (npt==0 || n2pt==0) rep.error("Couldn't find two real cells in column",0);
  
  //
  // Left Ghost Cell (doesn't get updated)
  //
  for (int v=0;v<SimPM.nvar;v++) { 
    slope_cpt[v] = 0.;
    slope_npt[v] = 0.;
    Fr_prev[v]   = 0.;
    Fr_this[v]   = 0.;
    edgeL[v]     = 0.;
    edgeR[v]     = 0.;
  }
  
  //
  // Now go through all cells in the column and calculate fluxes and add to dU vector.
  //
  do {
#ifdef TESTING
    dp.c = cpt;
#endif
#ifdef TEST_INT
    cout<<"First Cell:"; CI.print_cell(cpt);
    cout<<"Next Cell: "; CI.print_cell(npt);
#endif
    // Get the flux from left and right states, adding artificial
    // viscosity if needed.
    err += spatial_solver->SetEdgeState(cpt, posdir, SimPM.nvar,
                                slope_cpt, edgeL, csp, grid);
    err += spatial_solver->SetSlope(npt, axis, SimPM.nvar, slope_npt,
                                csp, grid);
    err += spatial_solver->SetEdgeState(npt, negdir, SimPM.nvar,
                                slope_npt, edgeR, csp, grid);
    err += spatial_solver->InterCellFlux(grid, cpt, npt, edgeL,
                                edgeR, Fr_this, SimPM.solverType,
                                SimPM.artviscosity, SimPM.gamma, dx);
    err += spatial_solver->dU_Cell(grid, cpt, axis, Fr_prev, Fr_this,
                                slope_cpt, csp, dx, dt);

    // record flux entering and leaving domain
    if (cpt->isbd_ref[negdir]) {
      for (int v=0;v<SimPM.nvar;v++) cpt->F[v] = Fr_this[v];
    }
    if (npt->isbd_ref[posdir]) {
      for (int v=0;v<SimPM.nvar;v++) npt->F[v] = Fr_this[v];
    }

#ifdef TEST_INT
    for (int v=0;v<SimPM.nvar;v++) {
      if(!isfinite(cpt->dU[v])) {
        rep.printVec("Fl",Fr_prev,SimPM.nvar);
        rep.printVec("Fr",Fr_this,SimPM.nvar);
        rep.printVec("El",edgeL,SimPM.nvar);
        rep.printVec("Er",edgeR,SimPM.nvar);
        cout <<"dt:"<<dt<<"\tdx="<<dx<<"\n";
        //  rep.printVec("dU",&cpt->dU[v],1);
        CI.print_cell(cpt);
        CI.print_cell(npt);
        rep.error("nans!!!",2);
      }
    }
#endif //TESTING

#ifdef TEST_CONSERVATION 
    // Track energy, momentum entering/leaving domain, if outside
    // boundary
    double dA=spatial_solver->CellInterface(cpt,posdir,dx);
    if(csp==SimPM.tmOOA &&
       pconst.equalD(grid->Xmin(axis),SimPM.Xmin[axis]) &&
       !(cpt->isdomain) && npt->isgd && npt->isleaf) {
//#ifdef TESTING
//      if (!pconst.equalD(Fr_this[MMX],0.0) && !pconst.equalD(dA,0.0)) {
//        cout <<"entering domain: F="<<Fr_this[MMX]<<". ";
//        cout <<", dA="<<dA<<"\n";
//        rep.printVec("pos",cpt->pos,SimPM.ndim);
//        rep.printVec("Ph",cpt->Ph,SimPM.nvar);
//        CI.print_cell(npt);
//      }
//#endif
      dM += Fr_this[RHO]*dt*dA;
      dE += Fr_this[ERG]*dt*dA;
      dMX += Fr_this[MMX]*dt*dA;
      dMY += Fr_this[MMY]*dt*dA;
      dMZ += Fr_this[MMZ]*dt*dA;
    }
    else if (csp==SimPM.tmOOA &&
       pconst.equalD(grid->Xmax(axis),SimPM.Xmax[axis]) &&
       cpt->isgd && cpt->isleaf && !npt->isdomain) {
      //cout <<"leaving domain\n";
      dM -= Fr_this[RHO]*dt*dA;
      dE -= Fr_this[ERG]*dt*dA;
      dMX -= Fr_this[MMX]*dt*dA;
      dMY -= Fr_this[MMY]*dt*dA;
      dMZ -= Fr_this[MMZ]*dt*dA;
    }
#endif // TEST_CONSERVATION
    //
    // Now move temp arrays for moving on to the next cell
    //
    temp=Fr_prev;
    Fr_prev = Fr_this;
    Fr_this = temp;
    temp = slope_cpt;
    slope_cpt = slope_npt;
    slope_npt = temp;
    cpt = npt; npt = n2pt;
  }  while ( (n2pt=grid->NextPt(n2pt,posdir)) );

  //
  // Now n2pt=null. npt = bd-data, cpt= (gd/bd)-data.
  // So have to do something different.
  //
#ifdef TESTING
  dp.c = cpt;
#endif
  // last cell 1st order.
  err += spatial_solver->SetEdgeState(
          cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
  for (int v=0;v<SimPM.nvar;v++) slope_npt[v] = 0.;
  err += spatial_solver->SetEdgeState(
          npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
  err += spatial_solver->InterCellFlux(
          grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType,
          SimPM.artviscosity, SimPM.gamma, dx);
  err += spatial_solver->dU_Cell(
          grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, dx, dt);
  // record flux entering and leaving domain
  if (cpt->isbd_ref[negdir]) {
    for (int v=0;v<SimPM.nvar;v++) cpt->F[v] = Fr_this[v];
  }
  if (npt->isbd_ref[posdir]) {
    for (int v=0;v<SimPM.nvar;v++) npt->F[v] = Fr_this[v];
  }

#ifdef TEST_CONSERVATION 
  // Track energy, momentum entering/leaving domain, if outside
  // boundary
  double dA=spatial_solver->CellInterface(cpt,posdir,dx);
  if (csp==SimPM.tmOOA &&
      pconst.equalD(grid->Xmax(axis),SimPM.Xmax[axis]) &&
      cpt->isgd && cpt->isleaf && !npt->isdomain) {
    //cout <<"leaving domain 2\n";
    dM -= Fr_this[RHO]*dt*dA;
    dE -= Fr_this[ERG]*dt*dA;
    dMX -= Fr_this[MMX]*dt*dA;
    dMY -= Fr_this[MMY]*dt*dA;
    dMZ -= Fr_this[MMZ]*dt*dA;
  }
#endif // TEST_CONSERVATION

  //
  // Right Ghost Cell-- have to calculate it's left interface differently,
  //
  cpt=npt;
#ifdef TESTING
  dp.c = cpt;
#endif
  for (int v=0;v<SimPM.nvar;v++) cpt->dU[v] +=0.; // nothing to calculate for it.
  rep.errorTest("(dU_Column) encountered an error",0,err);
  
  Fr_this   = mem.myfree(Fr_this);
  Fr_prev   = mem.myfree(Fr_prev);
  slope_cpt = mem.myfree(slope_cpt);
  slope_npt = mem.myfree(slope_npt);
  edgeL     = mem.myfree(edgeL);
  edgeR     = mem.myfree(edgeR);
  pstar     = mem.myfree(pstar);

#ifdef TEST_CONSERVATION
#ifdef PARALLEL
  dM = COMM->global_operation_double("SUM",dM);
  dE = COMM->global_operation_double("SUM",dE);
  dMX = COMM->global_operation_double("SUM",dMX);
  dMY = COMM->global_operation_double("SUM",dMY);
  dMZ = COMM->global_operation_double("SUM",dMZ);
  //cout <<"d="<<dM<<", "<<dE<<", "<<dMX<<", "<<dMY<<"\n";
#endif
  initMASS += dM;
  initERG  += dE;
  initMMX  += dMX;
  initMMY  += dMY;
  initMMZ  += dMZ;
#endif // TEST_CONSERVATION

  return 0;
}



// ##################################################################
// ##################################################################

  
int time_integrator::grid_update_state_vector(
      const double dt,  ///< timestep
      const int step, ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa,   ///< Full order of accuracy of simulation
      class GridBaseClass *grid ///< Computational grid.
      )
{
  int err=0;
  //
  // temp variable to handle change of energy when correcting for negative pressure.
  //
  pion_flt temperg =0.0;

  //
  // Loop through grid, updating Ph[] with CellAdvanceTime function.
  //
  class cell* c = grid->FirstPt_All();
  do {

#ifdef TESTING
    double dx = grid->DX();
    dp.ergTotChange = 0.;temperg =0.;
    dp.c = c;
#endif
    if (!c->isdomain) {
      // skip cell if is has been cut out of the domain.
      for (int v=0;v<SimPM.nvar;v++) c->dU[v]=0.0;
    }
    else {
      err += spatial_solver->CellAdvanceTime(c, c->P, c->dU, c->Ph,
                &temperg, SimPM.gamma, SimPM.EP.MinTemperature, dt);
    }

#ifdef TESTING
    if (err) {
      cout <<"______ Error in Cell-advance-time: ";
      CI.print_cell(c);
      CI.print_cell(c->npt);
      err=0;
    }
#else
    // ignore negative pressures and try to continue
    if (err) err=0;
#endif // TESTING
    
    //
    // If the current step is the full update, then also set
    // the updated base-state-vector to updated value.
    //
    if (step==ooa) {
      for (int v=0;v<SimPM.nvar;v++) c->P[v] = c->Ph[v];

#ifdef TESTING
      //
      // Update Total Energy from fixing negative pressures. Reset
      // update variables.
      //
      dp.ergTotChange = temperg;
      dp.initERG += dp.ergTotChange*spatial_solver->CellVolume(c,dx);
#endif // TESTING

    }

  } while ( (c =grid->NextPt_All(c)) !=0);

#ifdef TESTING
  cout <<"\tcalc_dynamics_dU done. error="<<err<<"\n";
#endif // TESTING
  return err; 
}


// ##################################################################
// ##################################################################








