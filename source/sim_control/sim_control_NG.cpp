/// \file sim_control_NG.cpp
/// 
/// \brief Simulation Control Class for Nested Grids.
/// 
/// \author Jonathan Mackey
/// 
/// This file contains the definitions of the member functions for
/// the NG-grid simulation control class.  This is built on top
/// of the control class for uniform grids, and so doesn't add too
/// much, just the moving up and down between levels.
/// 
/// Modifications:
/// - 2018.05.03 JM: Started on NG grid simulation control.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/command_line_interface.h"
#include "tools/reporting.h"
#include "tools/timer.h"
#include "constants.h"

#include "sim_control/sim_control_NG.h"

using namespace std;


//#define TEST_BC89FLUX

// ##################################################################
// ##################################################################



sim_control_NG::sim_control_NG()
{
#ifdef TESTING
  cout << "(sim_control_NG::Constructor)\n";
#endif
}



// ##################################################################
// ##################################################################



sim_control_NG::~sim_control_NG()
{
#ifdef TESTING
  cout << "(sim_control_NG::Destructor)\n";
#endif
}



// ##################################################################
// ##################################################################



int sim_control_NG::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
#ifdef TESTING
  cout <<"(sim_control_NG::Init) Initialising grid"<<"\n";
#endif
  int err=0;
  
  SimPM.typeofip=typeOfFile;
  setup_dataio_class(SimPM, typeOfFile);
  err = dataio->ReadHeader(infile, SimPM);
  rep.errorTest("(NG_INIT::get_parameters) error",0,err);

  // Check for commandline args that override the file parameters.
  err = override_params(narg, args);
  rep.errorTest("(NG_INIT::override_params) error",0,err);
  
  //
  // Set up the Xmin/Xmax/Range/dx of each level in the NG grid
  //
  setup_NG_grid_levels(SimPM);
  grid.resize(SimPM.grid_nlevels);
  err = setup_grid(grid,SimPM);
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(INIT::setup_grid) Something went wrong",0,err);

  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  // ----------------------------------------------------------------
  err = set_equations(SimPM);
  rep.errorTest("(NG_INIT::set_equations)",0,err);
  spatial_solver->SetEOS(SimPM.gamma);

  // ----------------------------------------------------------------
  err = setup_microphysics(SimPM);
  rep.errorTest("(NG_INIT::setup_microphysics)",0,err);
  
  // assign data to the grid from snapshot file.
  // ----------------------------------------------------------------
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(NG_INIT::assign_initial_data)",0,err);

  // For each grid in the NG grid, set Ph[] = P[],
  // and then implement the boundary conditions on the grid and
  // ghost cells.
  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    // Set Ph=P in every cell.
    cell *c = grid[l]->FirstPt();
    do {
      for(int v=0;v<SimPM.nvar;v++) c->Ph[v]=c->P[v];
    } while ((c=grid[l]->NextPt(c))!=0);

    if (SimPM.eqntype==EQGLM && SimPM.timestep==0) {
#ifdef TESTING
      cout <<"Initial state, zero-ing glm variable.\n";
#endif
      c = grid[l]->FirstPt(); do {
        c->P[SI] = c->Ph[SI] = 0.;//grid->divB(c);
      } while ( (c=grid[l]->NextPt(c)) !=0);
    }
  } // loop over levels

  // Assign boundary conditions to boundary points.
  // ----------------------------------------------------------------
  err = boundary_conditions(SimPM, grid);
  rep.errorTest("(NG_INIT::boundary_conditions) err!=0",0,err);

  // Setup Raytracing on each grid, if needed.
  // ----------------------------------------------------------------
  err += setup_raytracing(SimPM, grid);
  rep.errorTest("Failed to setup raytracer",0,err);

  // ----------------------------------------------------------------
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err = assign_boundary_data(SimPM, l, grid[l]);
    rep.errorTest("NG_INIT::assign_boundary_data",0,err);
  }
  // ----------------------------------------------------------------



  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"updating external boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateExternalBCs(SimPM,l,grid[l], spatial_solver,
                            SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("NG_INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TESTING
    cout <<"updating internal boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateInternalBCs(SimPM,l,grid[l], spatial_solver,
                            SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("NG_INIT: error from bounday update",0,err);
  // ----------------------------------------------------------------

  //
  // If testing the code, this calculates the momentum and energy on
  // the domain.
  //
  initial_conserved_quantities(grid);

  //
  // If using opfreq_time, set the next output time correctly.
  //
  if (SimPM.op_criterion==1) {
    if (SimPM.opfreq_time < TINYVALUE)
      rep.error("opfreq_time not set right and is needed!",
                SimPM.opfreq_time);
    SimPM.next_optime = SimPM.simtime+SimPM.opfreq_time;
    double tmp = 
      ((SimPM.simtime/SimPM.opfreq_time)-
       floor(SimPM.simtime/SimPM.opfreq_time))*SimPM.opfreq_time;
    SimPM.next_optime-= tmp;
  }

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    setup_dataio_class(SimPM,SimPM.typeofop);
    if (!dataio)
      rep.error("NG_INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

  
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



int sim_control_NG::initial_conserved_quantities(
      vector<class GridBaseClass *> &grid
      )
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION 
  pion_flt u[SimPM.nvar];
  //dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange = 0.0;
  //  cout <<"initERG: "<<dp.initERG<<"\n";
  initERG = 0.;  initMMX = initMMY = initMMZ = 0.; initMASS = 0.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    double dx = SimPM.levels[l].dx;
    double dv = 0.0;
    class cell *c=grid[l]->FirstPt();
    do {
      if (c->isdomain && c->isleaf) {
        dv = spatial_solver->CellVolume(c,dx);
        spatial_solver->PtoU(c->P,u,SimPM.gamma);
        initERG += u[ERG]*dv;
        initMMX += u[MMX]*dv;
        initMMY += u[MMY]*dv;
        initMMZ += u[MMZ]*dv;
        initMASS += u[RHO]*dv;
      }
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }

  cout <<"(conserved quantities) ["<< initERG <<", ";
  cout << initMMX <<", ";
  cout << initMMY <<", ";
  cout << initMMZ <<", ";
  cout << initMASS <<"]\n";

#endif // TEST_CONSERVATION
  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_NG::Time_Int(
      vector<class GridBaseClass *> &grid  ///< vector of grids
      )
{
  cout <<"-------------------------------------------------------\n";
  cout <<"(sim_control_NG::Time_Int) STARTING TIME INTEGRATION\n";
  cout <<"-------------------------------------------------------\n";
  int err=0;
  SimPM.maxtime=false;
  bool first_step=true;
  bool restart=true;
  if (SimPM.timestep!=0) first_step=false;
  if (SimPM.timestep==0) restart=false;
  clk.start_timer("Time_Int"); double tsf=0;

  // make sure all levels start at the same time.
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    SimPM.levels[l].dt = 0.0;
    SimPM.levels[l].simtime = SimPM.simtime;
  }

  // Do raytracing on all levels, and update internal and external
  // boundaries to populate the column densities correctly.
  // Even if there is not RT, this updates the boundaries.
  err = RT_all_sources_levels(SimPM);
  rep.errorTest("sim_control_NG: RT_all_sources_levels",0,err);
  if (SimPM.timestep==0) {
    cout << "(step=0) Writing initial data.\n";
    err=output_data(grid);
    rep.errorTest("Failed to write file... path?",0,err);
  }

  
  cout <<"-------------------------------------------------------\n";

  while (SimPM.maxtime==false) {

#if defined (CHECK_MAGP)
    calculate_magnetic_pressure(grid[0]);
#elif defined (BLAST_WAVE_CHECK)
    calculate_blastwave_radius(grid);
#endif

    // Get timestep on each level
    int scale = 1;
    double mindt = 1.0e99;
#ifdef DERIGS
    spatial_solver->set_max_speed(0.0);
#endif
    //err = RT_all_sources_levels(SimPM);
    //rep.errorTest("sim_control_NG: RT_all_sources_levels",0,err);
    
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
#ifdef TEST_INT
      cout <<"Calculate timestep, level "<<l<<", dx=";
      cout <<SimPM.levels[l].dx<<"\n";
#endif
      if (!restart && !first_step) {
        SimPM.last_dt = SimPM.levels[l].last_dt;
      }
      else {
        SimPM.levels[l].last_dt = SimPM.last_dt/
                                  SimPM.levels[l].multiplier;
      }

      err += calculate_timestep(SimPM, grid[l],spatial_solver,l);
      rep.errorTest("TIME_INT::calc_timestep()",0,err);
      
      mindt = std::min(mindt, SimPM.dt/scale);
#ifdef TEST_INT
      cout <<"level "<<l<<" got dt="<<SimPM.dt<<" and ";
      cout <<SimPM.dt/scale <<"\n";
#endif
      SimPM.levels[l].dt = SimPM.dt;
      scale *= 2;
    }
    // make sure all levels use same step (scaled by factors of 2).
    scale = 1;
    for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
      SimPM.levels[l].dt = mindt*scale;
      scale *= 2;
#ifdef TEST_INT
      cout <<"new dt="<<SimPM.levels[l].dt<<", t=";
      cout <<SimPM.levels[l].simtime<<"\n";
#endif
    }
    if (first_step) {
      // take a ~3x smaller timestep for the first timestep.
      for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
        //cout <<"level "<<l<<", orig dt="<<SimPM.levels[l].dt;
        SimPM.levels[l].dt *=0.3;
      }
      first_step=false;
    }
    if (restart) restart=false;
    SimPM.last_dt = SimPM.levels[0].last_dt;

    //
    // If using MHD with GLM divB cleaning, the following sets the
    // hyperbolic wavespeed.  If not, it does nothing.  By setting it
    // here and using t_dyn, we ensure that the hyperbolic wavespeed is
    // equal to the maximum signal speed on the grid, and not an
    // artificially larger speed associated with a shortened timestep.
    //
#ifdef DERIGS
    double cr=0.0;
    for (int d=0;d<SimPM.ndim;d++)
      cr += 1.0/(SimPM.levels[0].Range[d]*SimPM.levels[0].Range[d]);
    cr = M_PI*sqrt(cr);
    spatial_solver->Set_GLM_Speeds(SimPM.levels[0].dt,
                                   SimPM.levels[0].dx, cr);
#endif

    //clk.start_timer("advance_time");
    //
    // Use a recursive algorithm to update the coarsest level.  This
    // function also updates the next level twice, by calling itself
    // for the finer level, and so on.
    //
    advance_time(0,SimPM.levels[0].grid);
    SimPM.simtime = SimPM.levels[0].simtime;

#if ! defined (CHECK_MAGP)
#if ! defined (BLAST_WAVE_CHECK)
    cout <<"dt="<<SimPM.levels[0].dt<<"\tNew time: "<<SimPM.simtime;
    cout <<"\t timestep: "<<SimPM.timestep;
    tsf=clk.time_so_far("Time_Int");
    cout <<"\t runtime so far = "<<tsf<<" secs."<<"\n";
#endif
#endif

    err += check_energy_cons(grid);

    err+= output_data(grid);
    rep.errorTest("TIME_INT::output_data()",0,err);
    err+= check_eosim();
    rep.errorTest("TIME_INT::check_eosim()",0,err);
  }

  cout <<"(sim_control_NG::Time_Int) TIME_INT FINISHED.";
  cout <<" MOVING ON TO FINALISE SIM.\n";

  tsf=clk.time_so_far("Time_Int");
  cout <<"TOTALS ###: Nsteps="<<SimPM.timestep<<" wall-time=";
  cout <<tsf<<" time/step=";
  cout <<tsf/static_cast<double>(SimPM.timestep)<<"\n";
  cout <<"STEPS "<<SimPM.timestep;
  cout.setf( ios_base::scientific );
  cout.precision(6);
  cout <<"\t"<<tsf<<"\t"<<tsf/static_cast<double>(SimPM.timestep);
  cout <<"\t"<<static_cast<double>(SimPM.timestep*SimPM.Ncell)/tsf;
  cout <<"\n";
  cout <<"-------------------------------------------------------\n";

  return(0);
}



// ##################################################################
// ##################################################################





#ifdef CHECK_MAGP
///
/// This is only for a test problem -- it checks the magnetic
/// pressure on the full domain and outputs it to screen
///
void sim_control_NG::calculate_magnetic_pressure(
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      )
{
  //
  // Calculate the total magnetic pressure on the domain, normalised
  // to the initial value.
  //
  double magp=0.0, cellvol=0.0;
  static double init_magp=-1.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    
    cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd) 
        magp += (spatial_solver->Ptot(c->P,0.0) - c->P[PG])
                * spatial_solver->CellVolume(c);
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
void sim_control_NG::calculate_blastwave_radius(
      vector<class GridBaseClass *> &grid  ///< grid pointers.
      )
{
  //
  // Calculate the blast wave outer shock position.
  // If a NG grid, start on the finest grid and work outwards
  //
  double shockpos=0.0;
  static double old_pos=0.0;
  bool shock_found = false;
  //  static double last_dt=0.0;
  for (int l=SimPM.grid_nlevels-1; l>=0; l++) {

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


int sim_control_NG::Finalise(
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
  int err=0;
  cout <<"------------------------------------------------------------\n";
  cout <<"(sim_control::Finalise) FINALISING SIMULATION."<<"\n";
  err += check_energy_cons(grid);
  err+= output_data(grid);
  rep.errorTest("(FINALISE::output_data) Something went wrong",0,err);
  cout <<"\tSimTime = "<<SimPM.simtime<<"   #timesteps = "<<SimPM.timestep<<"\n";
#ifdef TESTING
  cout <<"(sim_control::Finalise) DONE.\n";
#endif
  cout <<"------------------------------------------------------------\n";
  return(0);
}


// ##################################################################
// ##################################################################



double sim_control_NG::advance_time(
      const int l,       ///< level to advance.
      class GridBaseClass *grid ///< grid pointer
      )
{
#ifdef TESTING
  cout <<"advance_time, level="<<l<<", starting.\n";
#endif

  double step=0.0;
  if (SimPM.tmOOA==1) {
    step = advance_step_OA1(l);
  }
  else if (SimPM.tmOOA==2) {
#ifdef TEST_INT
    //cout <<"Calling advance_step_OA2: level "<<l<<"\n";
#endif
    step = advance_step_OA2(l);
  }
  return step;
}




// ##################################################################
// ##################################################################



double sim_control_NG::advance_step_OA1(
      const int l       ///< level to advance.
      )
{
#ifdef TESTING
  cout <<"advance_step_OA1, level="<<l<<", starting.\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  class MCMDcontrol ppar; // unused for serial code.
  class GridBaseClass *grid = SimPM.levels[l].grid;

  err = update_evolving_RT_sources(
            SimPM,SimPM.levels[l].simtime,grid->RT);
  rep.errorTest("NG TIME_INT::update_RT_sources error",0,err);

  // take the first finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA1(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;

  // now calculate dU, change in conserved variables on this grid,
  // for this step.
  spatial_solver->Setdt(SimPM.levels[l].dt);
  // May need to do raytracing, if not already completed for getting
  // the timestep.
  if (!FVI_need_column_densities_4dt ||
    (SimPM.levels[l].step%SimPM.levels[l].multiplier !=0) ) {
    err += do_ongrid_raytracing(SimPM,grid,l);
    rep.errorTest("scn::advance_step_OA1: calc_rt_cols()",0,err);
  }
  err += calc_microphysics_dU(SimPM.levels[l].dt, grid);
  err += calc_dynamics_dU(SimPM.levels[l].dt,TIMESTEP_FIRST_PART,
                                                            grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(SimPM.levels[l].dt,
                                        TIMESTEP_FIRST_PART, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA1: calc_x_dU",0,err);

  // take the second finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA1(l+1);
  }

  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  // First correct fluxes
  //
  if (l < SimPM.grid_nlevels-1) {
    err += recv_BC89_fluxes_F2C(l,SimPM.levels[l].dt,OA1,OA1);
    rep.errorTest("scn::advance_step_OA1: recv_BC89_flux",0,err);
  }
  err += grid_update_state_vector(SimPM.levels[l].dt,OA1,OA1, grid);
  rep.errorTest("scn::advance_step_OA1: state-vec update",0,err);

  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l==0) SimPM.last_dt = SimPM.levels[l].dt;

  //
  // update internal and external boundaries.
  //
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                                      SimPM.simtime, OA1, OA1);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                                      SimPM.simtime, OA1, OA1);

#ifdef TESTING
  cout <<"advance_step_OA1, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step="<<SimPM.levels[l].step;
  cout <<", next dt="<<SimPM.levels[l].dt<<" next time=";
  cout << SimPM.levels[l].simtime + SimPM.levels[l].dt <<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}




// ##################################################################
// ##################################################################




double sim_control_NG::advance_step_OA2(
      const int l       ///< level to advance.
      )
{
#ifdef TEST_INT
  cout <<"advance_step_OA2, level="<<l<<", starting. ";
  cout <<SimPM.levels[l].simtime<<", step=";
  cout <<SimPM.levels[l].step<<"\n";
#endif
  int err=0;
  double dt2_fine=0.0; // timestep for two finer level steps.
  double dt2_this=0.0; // two timesteps for this level.
  class GridBaseClass *grid = SimPM.levels[l].grid;
  err = update_evolving_RT_sources(
            SimPM,SimPM.levels[l].simtime,grid->RT);
  rep.errorTest("NG TIME_INT::update_RT_sources error",0,err);

  // take the first finer grid step, if there is a finer grid.
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA2(l+1);
    
    // timestep for this level is equal to two steps of finer level,
    // where we take the sum of the fine step just taken and the next
    // step (not yet taken).
    SimPM.levels[l].dt = dt2_fine;
  }
  dt2_this = SimPM.levels[l].dt;

  // Predictor step: use 0.5*dt to get to time-centred state
  double dt_now = dt2_this*0.5;             // half of the timestep
  double ctime = SimPM.levels[l].simtime; // current time
  spatial_solver->Setdt(dt_now);
  // May need to do raytracing
  if (!FVI_need_column_densities_4dt ||
    (SimPM.levels[l].step%SimPM.levels[l].multiplier !=0) ) {
    err += do_ongrid_raytracing(SimPM,grid,l);
    rep.errorTest("scn::advance_step_OA2: calc_rt_cols()",0,err);
  }

  //
  // now calculate dU, the change in conserved variables on this grid
  // for the half step of the 2nd order step.
  //
#ifdef TEST_INT
  cout <<"l="<<l<<" half step, start calc_microphysics_dU\n";
#endif
  err += calc_microphysics_dU(dt_now, grid);
#ifdef TEST_INT
  cout <<"l="<<l<<" half step, start calc_dynamics_dU\n";
#endif
  err += calc_dynamics_dU(dt_now,TIMESTEP_FIRST_PART, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now,TIMESTEP_FIRST_PART, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA1",0,err);
#ifdef TEST_INT
  cout <<"l="<<l<<" half step, grid_update_state_vector\n";
#endif
  err += grid_update_state_vector(dt_now,
                                    TIMESTEP_FIRST_PART,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA1",0,err);  
  // Update boundary data.
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                                    ctime+dt_now, OA1, OA2);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                                    ctime+dt_now, OA1, OA2);
  rep.errorTest("scn::advance_step_OA2: bounday update OA1",0,err);

  //
  // Now calculate dU for the full step (OA2)
  //
  dt_now = dt2_this;  // full step
  spatial_solver->Setdt(dt_now);
  err += do_ongrid_raytracing(SimPM,grid,l);
  rep.errorTest("scn::advance_time: calc_rt_cols() OA2",0,err);
#ifdef TEST_INT
  cout <<"l="<<l<<" full step, start calc_microphysics_dU\n";
#endif
  err += calc_microphysics_dU(dt_now, grid);
#ifdef TEST_INT
  cout <<"l="<<l<<" full step, start calc_dynamics_dU\n";
#endif
  err += calc_dynamics_dU(dt_now, TIMESTEP_FULL, grid);
#ifdef THERMAL_CONDUCTION
  err += calc_thermal_conduction_dU(dt_now,TIMESTEP_FULL, grid);
#endif // THERMAL_CONDUCTION
  rep.errorTest("scn::advance_step_OA2: calc_x_dU OA2",0,err);

  // take the second finer grid step, if there is a finer grid.
#ifdef TEST_INT
  cout <<"l="<<l<<" full step, call 2nd l+1 update\n";
#endif
  if (l<SimPM.grid_nlevels-1) {
    dt2_fine = advance_step_OA2(l+1);
  }

  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  // First correct fluxes
  //
#ifdef TEST_INT
  cout <<"l="<<l<<" full step, grid_update_state_vector\n";
#endif
  spatial_solver->Setdt(dt2_this);
  if (l < SimPM.grid_nlevels-1) {
    err += recv_BC89_fluxes_F2C(l,dt_now,TIMESTEP_FULL,OA2);
    rep.errorTest("scn::advance_step_OA1: recv_BC89_flux",0,err);
  }
  err += grid_update_state_vector(dt_now,TIMESTEP_FULL,OA2, grid);
  rep.errorTest("scn::advance_step_OA2: state-vec update OA2",0,err);  

  // increment time and timestep for this level
  SimPM.levels[l].simtime += SimPM.levels[l].dt;
  SimPM.levels[l].step ++;
  if (l==SimPM.grid_nlevels-1) {
    SimPM.timestep ++;
  }
  SimPM.levels[l].last_dt = SimPM.levels[l].dt;
  if (l==0) SimPM.last_dt = SimPM.levels[l].dt;

  //
  // update internal and external boundaries.
  //
  err += TimeUpdateInternalBCs(SimPM, l, grid, spatial_solver,
                                SimPM.levels[l].simtime, OA2, OA2);
  err += TimeUpdateExternalBCs(SimPM, l, grid, spatial_solver,
                                SimPM.levels[l].simtime, OA2, OA2);


#ifdef TEST_INT
  cout <<"advance_step_OA2, level="<<l<<", returning. t=";
  cout <<SimPM.levels[l].simtime<<", step=";
  cout <<SimPM.levels[l].step<<"\n";
#endif
  return dt2_this + SimPM.levels[l].dt;
}



// ##################################################################
// ##################################################################



int sim_control_NG::check_energy_cons(
      vector<class GridBaseClass *> &grid
      )
{
  // Energy, and Linear Momentum in x-direction.
#ifdef TEST_CONSERVATION 
  pion_flt u[SimPM.nvar];
  nowERG=0.;
  nowMMX = 0.;
  nowMMY = 0.;
  nowMMZ = 0.;
  nowMASS = 0.0;
  double totmom=0.0;
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    double dx = SimPM.levels[l].dx;
    double dv = 0.0;
    class cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd && c->isgd) {
        dv = spatial_solver->CellVolume(c,dx);
        spatial_solver->PtoU(c->P,u,SimPM.gamma);
        nowERG += u[ERG]*dv;
        nowMMX += u[MMX]*dv;
        nowMMY += u[MMY]*dv;
        nowMMZ += u[MMZ]*dv;
        nowMASS += u[RHO]*dv;
        totmom += sqrt(u[MMX]*u[MMX] + u[MMY]*u[MMY] + u[MMZ]*u[MMZ])
                   *dv;
      }
    } while ( (c =grid[l]->NextPt(c)) !=0);
  }

  cout <<"(conserved quantities) ["<< nowERG <<", ";
  cout << nowMMX <<", ";
  cout << nowMMY <<", ";
  cout << nowMMZ <<", ";
  cout << nowMASS <<"]\n";
  cout <<"(relative error      ) [";
  cout << (nowERG-initERG)/(initERG) <<", ";
  cout << (nowMMX-initMMX)/(totmom) <<", ";
  cout << (nowMMY-initMMY)/(totmom) <<", ";
  cout << (nowMMZ-initMMZ)/(totmom) <<", ";
  cout << (nowMASS-initMASS)/initMASS <<"]\n";

#endif // TEST_CONSERVATION
  return(0);
}



// ##################################################################
// ##################################################################



int sim_control_NG::recv_BC89_fluxes_F2C(
      const int level,      ///< My level in grid hierarchy.
      const double dt,  ///< timestep
      const int step,   ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int ooa     ///< Full order of accuracy of simulation
      )
{
  //return 0;
  if (step != ooa) {
    cout <<"don't receive fluxes on half step\n";
    return 1;
  }
  if (level==SimPM.grid_nlevels-1) {
    rep.error("finest level trying to receive data from l+1",level);
  }
#ifdef TEST_BC89FLUX
  cout <<"level "<<level<<": correcting fluxes from finer grid\n";
#endif

  int err=0;
  class GridBaseClass *grid = SimPM.levels[level].grid;
  class GridBaseClass *fine = SimPM.levels[level].child;
  enum axes ax;
  
  // loop over boundaries, and, for each direction that has a
  // non-zero boundary, correct the coarse fluxes.
  for (unsigned int d=0; d<grid->flux_update_recv.size(); d++) {
#ifdef TEST_BC89FLUX
    cout <<"flux update: "<<d<<", looping through faces.\n";
#endif
    ax = static_cast<enum axes>(d/2);

    if (grid->flux_update_recv[d].fi[0] ==0) {
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      cout <<"flux update: "<<d<<", receiving.\n";
#endif
      err += recv_BC89_flux_boundary(grid, dt,
                      fine->flux_update_send[d],
                      grid->flux_update_recv[d], d, ax);
    }

#ifdef TEST_BC89FLUX
    cout <<"Direction: "<<d<<", finished.\n";
#endif
  } // loop over directions
  return err;
}



// ##################################################################
// ##################################################################



int sim_control_NG::recv_BC89_flux_boundary(
      class GridBaseClass *grid, ///< pointer to coarse grid
      const double dt,  ///< timestep
      struct flux_update &send, ///< data for fine grid
      struct flux_update &recv, ///< data for coarse grid
      const unsigned int d,   ///< direction of outward normal
      const axes ax ///< axis of normal direction.
      )
{ 
  struct flux_interface *fc=0;
  struct flux_interface *ff=0;
  double ftmp[SimPM.nvar],utmp[SimPM.nvar];
  for (int v=0;v<SimPM.nvar;v++) ftmp[v]=0.0;
  for (int v=0;v<SimPM.nvar;v++) utmp[v]=0.0;

  if (send.fi.size() != recv.fi.size()) {
    cout <<"send="<<send.fi.size()<<", recv="<<recv.fi.size()<<"\n";
    rep.error("fine and parent face arrays r different size",2);
  }

#ifdef TEST_BC89FLUX
  cout <<"SERIAL BC89 RECV: send-d="<<send.dir<<", recv-d=";
  cout <<recv.dir<<", d="<<d<<", ax="<<ax<<endl;
#endif
  for (unsigned int f=0; f<recv.fi.size(); f++) {
    fc = recv.fi[f];
    ff = send.fi[f];

#ifdef TEST_BC89FLUX
    cout <<"f="<<f<<":  grid="<<fc<<", flux =  ";
    rep.printVec("fc->flux",fc->flux,SimPM.nvar);
    cout <<"f="<<f<<":  fine="<<ff<<", flux =  ";
    rep.printVec("ff->flux",ff->flux,SimPM.nvar);
    for (int v=0;v<SimPM.nvar;v++) {
      fc->flux[v] /= fc->area[0];
    }
    for (int v=0;v<SimPM.nvar;v++) {
      ff->flux[v] /= fc->area[0];
    }
    cout <<"f="<<f<<":  grid="<<fc<<", flux =  ";
    rep.printVec("fc->flux",fc->flux,SimPM.nvar);
    cout <<"f="<<f<<":  fine="<<ff<<", flux =  ";
    rep.printVec("ff->flux",ff->flux,SimPM.nvar);
#endif
    
    for (int v=0;v<SimPM.nvar;v++) fc->flux[v] += ff->flux[v];

#ifdef TEST_BC89FLUX
    rep.printVec("     dU          ",fc->c[0]->dU,SimPM.nvar);
    rep.printVec("flux",fc->flux,SimPM.nvar);
    if (d==3) {
      rep.printVec("c state",fc->c[0]->Ph,SimPM.nvar);
      rep.printVec("+ state",grid->NextPt(fc->c[0],YP)->Ph,SimPM.nvar);
      rep.printVec("- state",grid->NextPt(fc->c[0],YN)->Ph,SimPM.nvar);
    }
#endif

    for (int v=0;v<SimPM.nvar;v++) fc->flux[v] /= fc->area[0];
    for (int v=0;v<SimPM.nvar;v++) ftmp[v]=0.0;
    
#ifdef TEST_BC89FLUX
    for (int v=0;v<SimPM.nvar;v++) fc->flux[v] *= fc->area[0];
#endif
    //
    // fc->flux is now the error in dU made for both coarse cells.
    // Correct dU in outer cell only, because inner cell is on 
    // top of fine grid and gets overwritten anyway.
    //
    // If we are at a negative boundary, then it is the positive
    // face of the cell that we correct, and vice versa for a
    // positive boundary.  So we call the DivStateVectorComponent
    // function with the flux in different order, depending on 
    // what direction the boundary is in.
    // This function is aware of geometry, so it calculates the
    // divergence correctly for curvilinear grids (important!!!).
    // The other face of the cell is set to zero flux.
    if (d%2 ==0) {
      spatial_solver->DivStateVectorComponent(fc->c[0], grid, ax,
                                  SimPM.nvar,ftmp,fc->flux,utmp);
    }
    else {
      spatial_solver->DivStateVectorComponent(fc->c[0], grid, ax,
                                  SimPM.nvar,fc->flux,ftmp,utmp);
    }
#ifdef TEST_BC89FLUX
    rep.printVec("**********  Error",utmp, SimPM.nvar);
    cout <<"Flux rho: "<<fc->flux[0]<<": "<<fc->c[0]->dU[0]<<", "<<utmp[0]<<"\n";
#endif

    for (int v=0;v<SimPM.nvar;v++) fc->c[0]->dU[v] += utmp[v];

#ifdef TEST_BC89FLUX
    spatial_solver->PtoU(fc->c[0]->P,utmp,SimPM.gamma);
    rep.printVec(" U",utmp,SimPM.nvar);
    rep.printVec("dU",fc->c[0]->dU,SimPM.nvar);
#endif
  } // loop over interfaces.
  return 0;
}



// ##################################################################
// ##################################################################



int sim_control_NG::do_ongrid_raytracing(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid, ///< Computational grid.
      const int l                ///< level in NG
      )
{
  if (!grid->RT) return 0;
  int err=0;
  //
  // If we have raytracing, we call the ray-tracing routines 
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    if (!SimPM.RS.sources[isrc].ongrid) continue;
#ifdef RT_TESTING
    cout <<"calc_raytracing_col_dens: SRC-ID: "<<isrc<<"\n";
#endif
    err += grid->RT->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      cout <<"isrc="<<isrc<<"\t"; 
      rep.error("ongrid RT: RT return",err);
    } // if error
  } // loop over sources
  return err;
}



// ##################################################################
// ##################################################################



int sim_control_NG::do_offgrid_raytracing(
      class SimParams &par,      ///< pointer to simulation parameters
      class GridBaseClass *grid,       ///< Computational grid.
      const int
      )
{
  if (!grid->RT) return 0;
  int err=0;
  //
  // If we have raytracing, we call the ray-tracing routines 
  // to get Tau0, dTau, Vshell in cell->extra_data[].
  //
  for (int isrc=0; isrc<par.RS.Nsources; isrc++) {
    if (SimPM.RS.sources[isrc].ongrid) {
      //cout <<"skipping source "<<isrc<<" b/c ongrid\n";
      continue;
    }
#ifdef RT_TESTING
    cout <<"calc_raytracing_col_dens: SRC-ID: "<<isrc<<"\n";
#endif
  err += grid->RT->RayTrace_Column_Density(isrc, 0.0, par.gamma);
    if (err) {
      cout <<"isrc="<<isrc<<"\t"; 
      rep.error("offgrid RT: RT return",err);
    } // if error
  } // loop over sources
  return err;
}



// ##################################################################
// ##################################################################



int sim_control_NG::RT_all_sources_levels(
      class SimParams &par  ///< simulation parameters
      )
{
  /// Do this in 3 passes: 1st we go from coarse to fine, tracing the
  /// off-grid sources.  This gets 1/2 of those rays right.
  /// Then go from fine to coarse, tracing all sources and 
  /// updating column densities as we go.
  /// Finally go from coarse to fine again, updating boundary data.
  int err=0;
  class GridBaseClass *grid = 0;

  // --------------------------------------------------------------
  // Update off-grid sources and external boundaries.
  for (int l=0; l<par.grid_nlevels; l++) {
#ifdef TEST_INT
    cout <<"updating external boundaries for level "<<l<<"\n";
#endif
    grid = par.levels[l].grid;
    err = TimeUpdateExternalBCs(par, l, grid,
              spatial_solver, par.simtime,par.tmOOA,par.tmOOA);
    rep.errorTest("NG RT_all_sources_levels: pass 1 BC-ext",0,err);
    err = do_offgrid_raytracing(par,grid,l);
    rep.errorTest("NG RT_all_sources_levels: pass 1 RT",0,err);
  }
  // --------------------------------------------------------------

  // --------------------------------------------------------------
  // update internal boundaries and then all sources (fine first)
  for (int l=par.grid_nlevels-1; l>=0; l--) {
#ifdef TEST_INT
    cout <<"updating internal boundaries for level "<<l<<"\n";
#endif
    grid = par.levels[l].grid;
    // F2C gets data from child grid onto this grid.
    err = TimeUpdateInternalBCs(par, l, grid, spatial_solver,
                            par.simtime,par.tmOOA,par.tmOOA);
    rep.errorTest("NG RT_all_sources_levels: pass 2 BC-int",0,err);
    err = do_ongrid_raytracing(par,grid,l);
    err = do_offgrid_raytracing(par,grid,l);
    rep.errorTest("NG RT_all_sources_levels: pass 2 RT",0,err);
  }
  rep.errorTest("sim_control_NG: internal boundary update",0,err);
  // --------------------------------------------------------------
  
  // --------------------------------------------------------------
  // Update external boundaries.
  for (int l=0; l<par.grid_nlevels; l++) {
#ifdef TEST_INT
    cout <<"Pass 3: external boundaries for level "<<l<<"\n";
#endif
    grid = par.levels[l].grid;
    // C2F gets data from parent grid onto this grid.
    err = TimeUpdateExternalBCs(par, l, grid,
              spatial_solver, par.simtime,par.tmOOA,par.tmOOA);
    rep.errorTest("NG RT_all_sources_levels: pass 3 BC-ext",0,err);
  }
  // --------------------------------------------------------------

  return err;
}



// ##################################################################
// ##################################################################




