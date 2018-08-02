/// \file sim_init_nested.cpp
/// \author Jonathan Mackey
/// \date 2018.05.10
///
/// Description:\n
/// Class declaration for sim_init_nested, which sets up a PION simulation
/// and gets everything ready to run.
///
/// Modifications:\n
/// - 2018.05.11 JM: moved code from sim_control.cpp
///

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "tools/reporting.h"
#include "tools/command_line_interface.h"
#include "sim_control/sim_init_nested.h"

#include "raytracing/raytracer_SC.h"
#include "microphysics/microphysics_base.h"
#include "spatial_solvers/solver_eqn_hydro_adi.h"
#include "spatial_solvers/solver_eqn_mhd_adi.h"


#include "dataIO/dataio_base.h"
#include "dataIO/dataio_text.h"
#ifdef SILO
#include "dataIO/dataio_silo_nestedgrid.h"
#endif // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif // if FITS

#include <climits>
using namespace std;



// ##################################################################
// ##################################################################


sim_init_nested::sim_init_nested()
{
#ifdef TESTING
  cout << "(sim_init_nested::Constructor)\n";
#endif
  return;
}



// ##################################################################
// ##################################################################


sim_init_nested::~sim_init_nested()
{
#ifdef TESTING
  cout << "(sim_init_nested::Destructor)\n";
#endif
  return;
}

// ##################################################################
// ##################################################################


int sim_init_nested::Init(
      string infile,
      int typeOfFile,
      int narg,
      string *args,
      vector<class GridBaseClass *> &grid  ///< address of vector of grid pointers.
      )
{
#ifdef TESTING
  cout <<"(sim_init_nested::Init) Initialising grid"<<"\n";
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
  err = setup_grid(grid,SimPM,&mpiPM);
  SimPM.dx = grid[0]->DX();
  rep.errorTest("(INIT::setup_grid) Something went wrong",0,err);

  //
  // All grid parameters are now set, so set up the appropriate
  // equations/solver class.
  //
  err = set_equations();
  rep.errorTest("(INIT::set_equations) err!=0 Fix me!",0,err);
  spatial_solver->set_dx(SimPM.dx);
  spatial_solver->SetEOS(SimPM.gamma);

  err = setup_microphysics(SimPM);
  rep.errorTest("(INIT::setup_microphysics) err!=0",0,err);
  
  //
  // Now assign data to the grid, either from file, or via some function.
  //
  err = dataio->ReadData(infile, grid, SimPM);
  rep.errorTest("(INIT::assign_initial_data) err!=0 Something went wrong",0,err);

  //
  // For each grid in the nested grid, set Ph[] = P[],
  // and then implement the boundary conditions on the grid and ghost cells.
  //
  for (int l=0; l<SimPM.grid_nlevels; l++) {
    spatial_solver->set_dx(SimPM.levels[l].dx);
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

  //
  // Assign boundary conditions to boundary points.
  //
  err = boundary_conditions(SimPM, grid);
  rep.errorTest("(INIT::boundary_conditions) err!=0",0,err);
  //
  // Setup Raytracing on each grid, if needed.
  //
  err += setup_raytracing(SimPM, grid);
  rep.errorTest("Failed to setup raytracer",0,err);

  // ----------------------------------------------------------------
  for (int l=0;l<SimPM.grid_nlevels;l++) {
    err = assign_boundary_data(SimPM, grid[l], SimPM.levels[l].parent, SimPM.levels[l].child);
    rep.errorTest("icgen_nest::assign_boundary_data",0,err);
  }
  // ----------------------------------------------------------------



  // ----------------------------------------------------------------
  for (int l=0; l<SimPM.grid_nlevels; l++) {
#ifdef TESTING
    cout <<"updating external boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateExternalBCs(SimPM, grid[l], l, SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("sim_init_nested: error from bounday update",0,err);
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  for (int l=SimPM.grid_nlevels-1; l>=0; l--) {
    spatial_solver->set_dx(SimPM.levels[l].dx);
#ifdef TESTING
    cout <<"updating internal boundaries for level "<<l<<"\n";
#endif
    err += TimeUpdateInternalBCs(SimPM, grid[l], l, SimPM.simtime,SimPM.tmOOA,SimPM.tmOOA);
  }
  rep.errorTest("sim_init_nested: error from bounday update",0,err);
  // ----------------------------------------------------------------

  //
  // If testing the code, this calculates the momentum and energy on the domain.
  //
  initial_conserved_quantities(grid);

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

  //
  // If outfile-type is different to infile-type, we need to delete
  // dataio and set it up again.
  //
  if (SimPM.typeofip != SimPM.typeofop) {
    if (dataio) {delete dataio; dataio=0;}
    if (textio) {delete textio; textio=0;}
    setup_dataio_class(SimPM.typeofop);
    if (!dataio) rep.error("INIT:: dataio initialisation",SimPM.typeofop);
  }
  dataio->SetSolver(spatial_solver);
  if (textio) textio->SetSolver(spatial_solver);

#ifdef SERIAL
  if (SimPM.timestep==0) {
    cout << "(INIT) Writing initial data.\n";
    err=output_data(grid);
    if (err)
      rep.error("Failed to write file!","maybe dir does not exist?");
  }
  cout <<"------------------------------------------------------------\n";
#endif // SERIAL
  
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



void sim_init_nested::setup_dataio_class(
      const int typeOfFile ///< type of I/O: 1=text,2=fits,5=silo
      )
{
  //
  // set up the right kind of data I/O class depending on the input.
  //
  switch (typeOfFile) {

#ifdef SILO
  case 5: // Start from Silo snapshot.
    dataio = new dataio_nested_silo (SimPM, "DOUBLE");
    break; 
#endif // if SILO

  default:
    rep.error("sim_init_nested::Init unhandled filetype",typeOfFile);
  }
  return;
}



// ##################################################################
// ##################################################################



int sim_init_nested::initial_conserved_quantities(
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
    spatial_solver->set_dx(SimPM.levels[l].dx);
    double dx = SimPM.levels[l].dx;
    double dv = 0.0;
    class cell *c=grid[l]->FirstPt();
    do {
      if (!c->isbd && c->isgd) {
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




