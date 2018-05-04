/// \file sim_control_nested.cpp
/// 
/// \brief Simulation Control Class for Nested Grids.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for
/// the nested-grid simulation control class.  This is built on top
/// of the control class for uniform grids, and so doesn't add too
/// much, just the moving up and down between levels.
/// 
/// Modifications:
/// - 2018.05.03 JM: Started on neste grid simulation control.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "nested_grid/sim_control_nested.h"

//#include "microphysics/microphysics_base.h"
//#include "raytracing/raytracer_SC.h"
//#include "dataIO/dataio_base.h"
//#ifdef SILO
//#include "dataIO/dataio_silo.h"
//#endif // if SILO
//#ifdef FITS
//#include "dataIO/dataio_fits.h"
//#endif // if FITS

//#include "spatial_solvers/solver_eqn_hydro_adi.h"
//#include "spatial_solvers/solver_eqn_mhd_adi.h"

//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <sys/time.h>
//#include <time.h>
//#include <climits>
using namespace std;




// ##################################################################
// ##################################################################


sim_control_nestedgrid::sim_control_nestedgrid()
{
}



// ##################################################################
// ##################################################################


sim_control_nestedgrid::~sim_control_nestedgrid()
{
#ifdef TESTING
  cout << "(sim_control_nestedgrid::Destructor) starting" <<"\n";
#endif
#ifdef TESTING
  cout << "(sim_control_nestedgrid::Destructor) Done." <<"\n";
#endif
}



// ##################################################################
// ##################################################################


int sim_control_nestedgrid::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
#ifdef TESTING
  cout <<"(sim_control_nestedgrid::Init) Initialising grid"<<"\n";
#endif
  int err=0;
  
  SimPM.typeofip=typeOfFile;
  setup_dataio_class(typeOfFile);
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("(INIT::get_parameters) err!=0 Something went wrong",0,err);
  
  // Now see if any commandline args override the Parameters from the file.
  err = override_params(narg, args);
  rep.errorTest("(INIT::override_params) err!=0 Something went wrong",0,err);
  
  //
  // Set up the Xmin/Xmax/Range/dx of each level in the nested grid
  //
  setup_nested_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
  RT.resize(SimPM.grid_nlevels);

  err = setup_grid(grid,SimPM,&mpiPM);
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(INIT::setup_grid) Something went wrong",0,err);

  err = set_equations();
  rep.errorTest("(INIT::set_equations) err!=0 Fix me!",0,err);
  spatial_solver->set_dx(SimPM.dx);
  spatial_solver->SetEOS(SimPM.gamma);

  err = setup_microphysics(SimPM);
  rep.errorTest("(INIT::setup_microphysics) err!=0",0,err);
  
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(INIT::assign_initial_data) err!=0 Something went wrong",0,err);

  //
  // For each grid in the nested grid, set Ph[] = P[],
  // and then implement the boundary conditions on the grid and ghost cells.
  //
  for (int l=0; l<SimPM.grid_nlevels; l++) {

    spatial_solver->set_dx(SimPM.nest_levels[l].dx);
    CI.set_dx(SimPM.nest_levels[l].dx);

    //
    // Set Ph=P in every cell.
    //
    cell *c = grid[l]->FirstPt();
    do {
      for(int v=0;v<SimPM.nvar;v++) c->Ph[v]=c->P[v];
    } while ((c=grid[l]->NextPt(c))!=0);
    //
    // If I'm using the GLM method, make sure Psi variable is initialised.
    //
    if (SimPM.eqntype==EQGLM && SimPM.timestep==0) {
#ifdef TESTING
      cout <<"Initial state, zero-ing glm variable.\n";
#endif
      for (int l=0; l<SimPM.grid_nlevels; l++) {
        c = grid[l]->FirstPt(); do {
          c->P[SI] = c->Ph[SI] = 0.;//grid->divB(c);
        } while ( (c=grid[l]->NextPt(c)) !=0);
      }
    }

  }

  err = boundary_conditions(grid, SimPM);
  rep.errorTest("(INIT::boundary_conditions) error",0,err);
  err += setup_raytracing(SimPM, grid, RT);
  rep.errorTest("(INIT::setup_raytracing) error",0,err);


  //
  // If using opfreq_time, set the next output time correctly.
  //
  if (SimPM.op_criterion==1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!",SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
    double tmp = 
      ((SimPM.simtime/SimPM.opfreq_time)-
       floor(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
    SimPM.next_optime-= tmp;
  }

  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    dataio->SetSolver(spatial_solver);
    if (!dataio) rep.error("INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

  if (SimPM.timestep==0) {
    cout << "(INIT) Saving initial data.\n";
    err=output_data(grid);
    if (err)
      rep.error("Failed to write file!","maybe dir does not exist?");
  }
  cout <<"------------------------------------------------------------\n";
  
#ifdef TESTING
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    cell *c = (grid[l])->FirstPt_All();
    do {
      if (pconst.equalD(c->P[RO],0.0)) {
        cout <<"zero data in cell: ";
        CI.print_cell(c);
      }
    } while ( (c=(grid[l])->NextPt_All(c)) !=0 );
  }
#endif // TESTING
  
  return(0);
}



// ##################################################################
// ##################################################################




int sim_control_nestedgrid::Time_Int(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control_nestedgrid::Time_Int) STARTING TIME INTEGRATION."<<"\n";
  cout <<"------------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  clk.start_timer("Time_Int"); double tsf=0;
  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid[0]);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid);
#endif
    //
    // Update RT sources.
    //
    err = update_evolving_RT_sources(SimPM,RT[0]);
    rep.errorTest("(TIME_INT::update_evolving_RT_sources()) error",0,err);

    //clk.start_timer("advance_time");
    err+= advance_time(grid,RT);
    //cout <<"advance_time took "<<clk.stop_timer("advance_time")<<" secs.\n";
    rep.errorTest("TIME_INT::advance_time()",0,err);

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.dt<<"\tNew time: "<<SimPM.simtime<<"\t timestep: "<<SimPM.timestep;
    tsf=clk.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err+= output_data(grid);
    rep.errorTest("TIME_INT::output_data()",0,err);
    err+= check_eosim();
    rep.errorTest("TIME_INT::check_eosim()",0,err);
  }

  cout <<"(sim_control_nestedgrid::Time_Int) TIME_INT FINISHED.  MOVING ON TO FINALISE SIM.\n";

  tsf=clk.time_so_far("Time_Int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep<<" wall-time=";
  cout <<tsf<<" time/step="<<tsf/static_cast<double>(SimPM.timestep)<<"\n";
  cout <<"STEPS "<<SimPM.timestep;
  cout.setf( ios_base::scientific );
  cout.precision(6);
  cout <<"\t"<<tsf<<"\t"<<tsf/static_cast<double>(SimPM.timestep);
  cout <<"\t"<<static_cast<double>(SimPM.timestep*SimPM.Ncell)/tsf<<"\n";
  cout <<"------------------------------------------------------------\n";

  return(0);
}



// ##################################################################
// ##################################################################



#ifdef CHECK_MAGP
///
/// This is only for a test problem -- it checks the magnetic
/// pressure on the full domain and outputs it to screen
///
void sim_control_nestedgrid::calculate_magnetic_pressure(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the total magnetic pressure on the domain, normalised to the
  // initial value.
  //
  double magp=0.0, cellvol=0.0;
  static double init_magp=-1.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    spatial_solver->set_dx(SimPM.nest_levels[l].dx);
    CI.set_dx(SimPM.nest_levels[l].dx);
    
    cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd) 
        magp += (spatial_solver->Ptot(c->P,0.0) - c->P[PG]) * spatial_solver->CellVolume(c);
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }
  if (init_magp<0) init_magp = magp;
  cout <<SimPM.simtime<<"\t"<<magp/init_magp<<"\t"<<magp<<"\n";
  return;
}
#endif // CHECK_MAGP



// ##################################################################
// ##################################################################



#ifdef BLAST_WAVE_CHECK
///
/// If running a 1D spherical blast wave, calculate the shock position
/// and output to screen.
///
void sim_control_nestedgrid::calculate_blastwave_radius(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  //
  // Calculate the blast wave outer shock position.
  // If a nested grid, start on the finest grid and work outwards
  //
  double shockpos=0.0;
  static double old_pos=0.0;
  bool shock_found = false;
  //  static double last_dt=0.0;
  for (int l=SimPM.grid_nlevels-1; l>=0; l++) {
    spatial_solver->set_dx(SimPM.nest_levels[l].dx);
    CI.set_dx(SimPM.nest_levels[l].dx);

    if (shock_found) continue;
    cell *c=grid->LastPt();
    if (fabs(c->P[VX])>=1.0e4) {
      cout<<"level "<<l<<" does not contain shock.\n";
    }
    else {
      do {
        c = grid->NextPt(c,RNsph);
        //cout <<c->id<<", vx="<<c->P[VX]<<"\n";
      } while ( c!=0 && fabs(c->P[VX])<1.0e4);
      if (c && (c->P[VX] >= 1.0e4)) {
        shockpos = CI.get_dpos(c,Rsph);
        shock_found=true;
      }
    }
  }
  if (pconst.equalD(old_pos,0.0))
    old_pos = shockpos;
  cout <<SimPM.simtime<<"\t"<<shockpos;
  //cout <<"\t"<<(shockpos-old_pos)/(SimPM.dt+TINYVALUE);
  cout <<"\n";
  old_pos=shockpos;
  return;
}
#endif // BLAST_WAVE_CHECK



// ##################################################################
// ##################################################################



int sim_control_nestedgrid::advance_time(
      vector<class GridBaseClass *> &g,  ///< grid pointer
      vector<class RayTracingBase *> &r  ///< raytracer for this grid.
      )
{
  return 0;
}


// ##################################################################
// ##################################################################




