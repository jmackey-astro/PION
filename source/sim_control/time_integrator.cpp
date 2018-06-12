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

int time_integrator::advance_time(
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
  }

  else if (SimPM.tmOOA==OA2 && SimPM.spOOA==OA2) {
    //cout <<"Second order update\n";
    err += first_order_update( 0.5*SimPM.dt, OA2, grid);
    if (err)
      rep.error("1st order time-update returned error",err);
    // Update boundary data.
    err += TimeUpdateInternalBCs(SimPM, grid, SimPM.simtime,   OA1, OA2);
    err += TimeUpdateExternalBCs(SimPM, grid, SimPM.simtime,   OA1, OA2);
    if (err) 
      rep.error("second_order_update: error from bounday update",err);

    err += second_order_update(SimPM.dt,     OA2, grid);
    if (err)
      rep.error("Second order time-update returned error",err);
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

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::first_order_update(
      const double dt,
      const int   ooa,
      class GridBaseClass *grid ///< Computational grid.
      )
{
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
    err += calculate_raytracing_column_densities(SimPM,grid->RT);
    if (err) 
      rep.error("first_order_update: error from first calc_rt_cols()",err);
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(dt, grid);
  err += calc_dynamics_dU(dt,OA1, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt,OA1, grid);
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
    err += calculate_raytracing_column_densities(SimPM,grid->RT);
    if (err) {
      rep.error("second_order_update: error from first calc_rt_cols()",err);
    }
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_microphysics_dU(      dt,      grid);
  err += calc_dynamics_dU(          dt, OA2, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt, OA2, grid);
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
  if (!RT)
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
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in calc_RT_microphysics_dU() c->isbd.\n";
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
#ifdef TESTING
        //cout <<"HEAT: Vs="<<FVI_heating_srcs[v].Vshell<<", dS="<<FVI_heating_srcs[v].dS<<", dC="<<FVI_heating_srcs[v].DelCol<<", Col="<<FVI_heating_srcs[v].Column<<"\n";
#endif
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
    // Check if cell is boundary data or not (can only be an internal boundary, such as
    // a stellar wind, since we are looping over cells which are grid data).  If it is
    // boundary data, then we don't want to update anything, so we skip it
    //
    if (c->isbd) {
#ifdef TESTING
      cout <<"skipping cell "<<c->id<<" in calc_noRT_microphysics_dU() c->isbd.\n";
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
      const int space_ooa, ///< spatial order of accuracy for update.
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
  if (space_ooa==OA1)
    cout <<"*****Updating dynamics: OA1\n";
  else if (space_ooa==OA2)
    cout <<"*****Updating dynamics: OA2\n";
  else rep.error("Bad ooa",space_ooa);
#endif //TESTING

  //
  // First we pre-process the cells, if needed.  This is required for
  // genuinely multi-dimensional viscosity such as Lapidus-like AV or
  // the H-Correction.
  //
  err = spatial_solver->preprocess_data(space_ooa, SimPM, grid);

  //
  // Now calculate the directionally-unsplit time update for the
  // conserved variables:
  //
  err = set_dynamics_dU(dt, space_ooa, grid); //,time_ooa);
  rep.errorTest("calc_dynamics_dU() set_dynamics_dU returned error.", 0,err);

  //
  // Post-processing is for if we are doing something like Constrained
  // Transport, where we have to change the B-field update.  At the
  // moment there is *NO* solver which does anything here since I
  // found the Dedner et al. (2002) divergence cleaning to be more
  // robust than e.g. Toth (2000) Field-CT method.  (well the internal
  // energy solver uses it, but it's not really worth using).
  //
  err = spatial_solver->PostProcess_dU(dt, space_ooa, grid); //,time_ooa);
  rep.errorTest("calc_dynamics_dU() spatial_solver->PostProcess_dU()",0,err);

  return 0;
}


// ##################################################################
// ##################################################################



int time_integrator::set_dynamics_dU(
      const double dt,     ///< timestep for this calculation
      const int space_ooa, ///< space OOA for this calculation
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
    class cell *cpt    = grid->FirstPt();
    class cell *marker = grid->FirstPt();
    
#ifdef TESTING
    rep.printVec("cpt",cpt->pos,SimPM.ndim);
    rep.printVec("+XX",(grid->NextPt(cpt,XP))->pos,SimPM.ndim);
    if (SimPM.ndim>1) rep.printVec("+YY",(grid->NextPt(cpt,YP))->pos,SimPM.ndim);
    if (SimPM.ndim>2) rep.printVec("+ZZ",(grid->NextPt(cpt,ZP))->pos,SimPM.ndim);
#endif
    
    while (
      (return_value = dynamics_dU_column(cpt,posdirs[i],negdirs[i], dt,
#ifdef TESTING
      // this is a hack, assuming spatial o-o-a is the same as the
      // time o-o-a.  But this is only needed for checking energy
      // and momentum conservation (to know if we are on the full
      // or half step) so it is not too important.
                                        space_ooa,
#endif
                                        space_ooa, grid)) ==0) {
      if ( !(cpt=grid->NextPt(cpt,posdirs[(i+1)%SimPM.ndim]))->isgd ) {
        if ( !(cpt=grid->NextPt(marker,posdirs[(i+2)%SimPM.ndim]))->isgd ) {
          CI.print_cell(cpt);
          rep.error("set_dynamics_dU: Got to edge of box before last point!",cpt);
        }
        marker = cpt;
      } // if null pointer.
    } // Loop over columns.
    if (return_value!=-1) rep.error("dUdtColumn returned abnormally.",return_value);
  } // Loop over three directions.
  spatial_solver->SetDirection(axis[0]); // Reset fluxes to x-dir, (just to be safe!).
  return 0;
}   // set_dynamics_dU()



// ##################################################################
// ##################################################################

  

int time_integrator::dynamics_dU_column
      (
      const class cell *startingPt, ///< sterting point of column.
      const enum direction posdir, ///< direction to trace column.
      const enum direction negdir, ///< reverse direction
      const double dt, ///< timestep we are advancing by.
#ifdef TESTING
      const int ctm, ///< time order-of-accuracy (for conservation)
#endif
      const int csp,  ///< spatial order-of-accuracy for this step.
      class GridBaseClass *grid ///< Computational grid.
      )
{
  if ( (SimPM.spOOA>2) || (SimPM.tmOOA>2) || (csp>2)  ) {
    cerr<<"(RSMethod::calc_dUdt) Error, only know 1st and 2nd order accurate methods.\n";
    return(1);
  }
  int err = 0;
#ifdef TESTING
  int ct=0;
#endif
  enum axes axis = spatial_solver->GetDirection();
  double dx = grid->DX();

  // Calculate Flux at positive (right) boundary of cell and store in temporary arrays
  // for the current cell (Fr_this) and the negative neighbour (Fr_prev)
  // Have to do it this way b/c ISO C++ forbids re-assignment of arrays.
  pion_flt *Fr_this=0, *Fr_prev=0, *temp=0, *slope_cpt=0, *slope_npt=0, *edgeR=0, *edgeL=0, *pstar=0;
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
  cell *cpt = grid->NextPt(startingPt,posdir); //=grid->NextPt(startingPt,negdir);
  while (grid->NextPt(cpt,negdir)) {cpt = grid->NextPt(cpt,negdir);} // CI.print_cell(cpt);}
  if(cpt==0) {cerr<<"(RSMethod::calc_dUdt) error finding left boundary cell.\n";return(1);}
  cell *npt  = grid->NextPt(cpt,posdir);
  cell *n2pt = grid->NextPt(npt,posdir);
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
    cout<<"First Cell:"; CI.print_cell(cpt);
    cout<<"Next Cell: "; CI.print_cell(npt);
#endif
    // Get the flux from left and right states, adding artificial viscosity if needed.
    err += spatial_solver->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
    err += spatial_solver->SetSlope(npt, axis, SimPM.nvar, slope_npt, csp, grid);
    err += spatial_solver->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
    err += spatial_solver->InterCellFlux(grid, cpt, npt, edgeL, edgeR, Fr_this,
                              SimPM.solverType, SimPM.artviscosity, SimPM.gamma, dx);
    err += spatial_solver->dU_Cell(grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, dx, dt);

#ifdef TESTING
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
    // Track energy, momentum entering domain.
    if(ctm==SimPM.tmOOA && !(cpt->isgd) && npt->isgd) {
      ct++; if (ct>1) rep.error("Entering domain more than once! (dUcolumn)",ct);
      //      cout <<"Entering Domain Dir = "<<posdir<<" and interface area = "<<spatial_solver->CellInterface(cpt,posdir)<<"\n";
      dp.initERG += Fr_this[ERG]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMX += Fr_this[MMX]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMY += Fr_this[MMY]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMZ += Fr_this[MMZ]*dt*spatial_solver->CellInterface(cpt,posdir);
      //      if (posdir==YP && fabs(Fr_this[MMY])>2*MACHINEACCURACY) {
      //  cout <<"R-momentum flux entering domain from R=0(!) = "<<Fr_this[MMY]<<"\n";
      //  cout <<"v_R in first cell = "<<npt->Ph[VY]<<", "<<npt->P[VY]<<"\n";
      //      }
    }
    else if (ctm==SimPM.tmOOA && !(npt->isgd) && cpt->isgd) {
      ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
      //      cout <<"Leaving Domain Dir = "<<posdir<<" and interface area = "<<spatial_solver->CellInterface(cpt,posdir)<<"\n";
      dp.initERG -= Fr_this[ERG]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMX -= Fr_this[MMX]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMY -= Fr_this[MMY]*dt*spatial_solver->CellInterface(cpt,posdir);
      dp.initMMZ -= Fr_this[MMZ]*dt*spatial_solver->CellInterface(cpt,posdir);
    }
#endif //TESTING


    //
    // Now move temp arrays to prepare for moving on to the next cell.
    //
    temp=Fr_prev;
    Fr_prev = Fr_this;
    Fr_this = temp; // just point to the free memory to be overwritten next step.
    temp = slope_cpt;
    slope_cpt = slope_npt;
    slope_npt = temp;
    cpt = npt; npt = n2pt;
  }  while ( (n2pt=grid->NextPt(n2pt,posdir)) );
  
  
  //
  // Now n2pt=null. npt = bd-data, cpt= (gd/bd)-data. So have to do something different.
  //
#ifdef TESTING
  dp.c = cpt;
#endif

  err += spatial_solver->SetEdgeState(cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
  for (int v=0;v<SimPM.nvar;v++) slope_npt[v] = 0.; // last cell must be 1st order.
  err += spatial_solver->SetEdgeState(npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
  err += spatial_solver->InterCellFlux(grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.solverType, SimPM.artviscosity, SimPM.gamma, dx);
  err += spatial_solver->dU_Cell(grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, dx, dt);

#ifdef TESTING
  if (ctm==SimPM.tmOOA && cpt->isgd && !(npt->isgd)) {
    ct++; if (ct>2) rep.error("Leaving domain more than once! (dUcolumn)",ct);
    dp.initERG -= Fr_this[ERG]*dt*spatial_solver->CellInterface(cpt,posdir);
    dp.initMMX -= Fr_this[MMX]*dt*spatial_solver->CellInterface(cpt,posdir);
    dp.initMMY -= Fr_this[MMY]*dt*spatial_solver->CellInterface(cpt,posdir);
    dp.initMMZ -= Fr_this[MMZ]*dt*spatial_solver->CellInterface(cpt,posdir);
  }
#endif //TESTING
 
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

  //
  // Check if this is the last column or not. (first track back from bd to grid).
  // If it is the last column, return -1 instead of 0 to indicate this. (A positive
  // return value from errors is picked up in the calling function).
  //
  do{} while( !(cpt=grid->NextPt(cpt,negdir))->isgd );
  if (cpt->id == grid->LastPt()->id) return(-1);
  else return(0);
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
  class cell* c = grid->FirstPt();
  do {

#ifdef TESTING
    dp.ergTotChange = 0.;temperg =0.;
    dp.c = c;
#endif
    
    err += spatial_solver->CellAdvanceTime(c, c->P, c->dU, c->Ph,
                &temperg, SimPM.gamma, SimPM.EP.MinTemperature, dt);

#ifdef TESTING
    if (err) {
      cout <<"______ Error in Cell-advance-time: ";
      CI.print_cell(c);
      err=0;
    }
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
      dp.initERG += dp.ergTotChange*spatial_solver->CellVolume(c);
#endif // TESTING

    }

  } while ( (c =grid->NextPt(c)) !=0);

#ifdef TESTING
  cout <<"\tcalc_dynamics_dU done. error="<<err<<"\n";
#endif // TESTING
  return err; 
}


// ##################################################################
// ##################################################################








