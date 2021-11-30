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
/// - 2021.07.13 JM: updated advance_time() so that boundary data is updated at
/// the correct time.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#include "dataIO/dataio_base.h"
#include "sim_control/time_integrator.h"

#include "microphysics/microphysics_base.h"
#include "raytracing/raytracer_SC.h"

#ifdef SILO
#include "dataIO/dataio_silo.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif  // if FITS

#include "spatial_solvers/solver_eqn_base.h"

#include <fstream>

#include <sstream>
#include <sys/time.h>
#include <time.h>

#ifdef PION_OMP
#include <omp.h>
#endif

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
    const int level,           ///< level in grid hierarchy
    class GridBaseClass *grid  ///< Computational grid.
)
{
  int err = 0;
  //
  // Check order-of-accuracy (OOA) requested, and perform the
  // appropriate update.
  //
  if (SimPM.tmOOA == OA1 && SimPM.spOOA == OA1) {
    //
    // Send in full timestep, and the order-of-accuracy, so the
    // function knows whether to update P[i] as well as Ph[i]
    //
    // cout <<"First order update\n";
    err += first_order_update(SimPM.dt, SimPM.tmOOA, grid);
    if (err)
      spdlog::error("{}: {}", "first_order_update() returned error", err);

    // Update boundary data to new state
    err += TimeUpdateInternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + SimPM.dt, OA1, OA1);
    err += TimeUpdateExternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + SimPM.dt, OA1, OA1);
    if (err)
      spdlog::error(
          "{}: {}", "second_order_update: error from bounday update", err);
  }

  else if (SimPM.tmOOA == OA2 && SimPM.spOOA == OA2) {
    // cout <<"Second order update\n";
    // take the 1st-order half-step to predict the time-centred state.
    err += first_order_update(0.5 * SimPM.dt, OA2, grid);
    if (err)
      spdlog::error("{}: {}", "1st order time-update returned error", err);

    // Update boundary data to intermediate state
    err += TimeUpdateInternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + 0.5 * SimPM.dt, OA1,
        OA2);
    err += TimeUpdateExternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + 0.5 * SimPM.dt, OA1,
        OA2);
    if (err)
      spdlog::error(
          "{}: {}", "second_order_update: error from bounday update", err);

    // take the time-centred 2nd-order step.
    err += second_order_update(SimPM.dt, OA2, grid);
    if (err)
      spdlog::error("{}: {}", "Second order time-update returned error", err);

    // Update boundary data to new state.
    err += TimeUpdateInternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + 1.0 * SimPM.dt, OA2,
        OA2);
    err += TimeUpdateExternalBCs(
        SimPM, level, grid, spatial_solver, SimPM.simtime + 1.0 * SimPM.dt, OA2,
        OA2);
    if (err)
      spdlog::error(
          "{}: {}", "second_order_update: error from bounday update", err);
  }
  //
  // Add in 3rd order PPM at some stage???
  //
  else {
    spdlog::error(
        "{}: {}", "Bad OOA requests; choose (1,1) or (2,2)", SimPM.tmOOA);
  }

  // Update time variables to new values.
  SimPM.simtime += SimPM.dt;
  SimPM.last_dt = SimPM.dt;
  SimPM.timestep++;


  return SimPM.dt;
}

// ##################################################################
// ##################################################################

int time_integrator::first_order_update(
    const double dt,
    const int ooa,
    class GridBaseClass *grid  ///< Computational grid.
)
{
  // NB Only used for uniform grid.  update for NG grid is in
  // sim_control_NG.cpp
  int err = 0;
  //
  // Set dt for equations class
  //
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt);
#ifdef PION_OMP
  }
#endif

  //
  // May need to do raytracing, if it wasn't needed for calculating
  // the timestep.
  //
  if (!FVI_need_column_densities_4dt && grid->RT) {
#ifdef RT_TESTING
    spdlog::info(" 1st order step doing RT ");
#endif
    err += RT_all_sources(SimPM, grid, 0);
    if (err)
      spdlog::error("{}: {}", "first_order_update: first calc_rt_cols()", err);
  }

  //
  // Calculate updates for each physics module
  //
  err += calc_conduction_dU(dt, OA1, grid);
  err += calc_microphysics_dU(dt, grid);
  err += calc_dynamics_dU(dt, OA1, grid);
  if (err)
    spdlog::error("{}: {}", "first_order_update: error from calc_*_dU", err);

  //
  // Now update Ph[i] to new values (and P[i] also if full step).
  //
  err += grid_update_state_vector(dt, OA1, ooa, grid);
  if (err) spdlog::error("{}: {}", "first_order_update: state-vec update", err);

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::second_order_update(
    const double dt,
    const int ooa,
    class GridBaseClass *grid  ///< Computational grid.
)
{
  // NB Only used for uniform grid.  update for NG grid is in
  // sim_control_NG.cpp
  int err = 0;
  // Set dt for equations class
#ifdef PION_OMP
  #pragma omp parallel
  {
#endif
    spatial_solver->Setdt(dt);
#ifdef PION_OMP
  }
#endif

  // Raytracing, to get column densities for microphysics update.
  err += RT_all_sources(SimPM, grid, 0);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "second_order_update: RT", 0, err);

  // Calculate updates for each physics module
  err += calc_conduction_dU(dt, OA2, grid);
  err += calc_microphysics_dU(dt, grid);
  err += calc_dynamics_dU(dt, OA2, grid);
  if (err)
    spdlog::error("{}: {}", "second_order_update: error from calc_*_dU", err);

  // Now update Ph[i] to new values (and P[i] also if full step).
  err += grid_update_state_vector(dt, OA2, ooa, grid);
  if (err)
    spdlog::error("{}: {}", "second_order_update: state-vec update", err);

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::calc_conduction_dU(
    const double delt,         ///< timestep to integrate
    const int step,            ///< whether OA1 or OA2.
    class GridBaseClass *grid  ///< Computational grid.
)
{
  int err = 0;
  // ----------------------------------------------------------------
  // cout <<"\ttime_integrator::calc_conduction_dU setting TC Edot\n";
  if (step != OA1) err = set_thermal_conduction_Edot(SimPM, step, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "calc_conduction_dU: set_thermal_conduction_Edot()", 0, err);
  // We need to multiply dU[ERG] by dt to convert from Edot to Delta-E
  class cell *c = grid->FirstPt();
  do {
    c->dU[ERG] *= delt;
  } while ((c = grid->NextPt(c)) != 0);
  // ----------------------------------------------------------------
  return err;
}



// ##################################################################
// ##################################################################



int time_integrator::calc_microphysics_dU(
    const double delt,         ///< timestep to integrate MP eqns.
    class GridBaseClass *grid  ///< Computational grid.
)
{
  // cout <<"\tcalc_microphysics_dU starting.\n";

  //
  // If we are not doing microphysics or raytracing, the modules will
  // not be initialised, so just return zero.  (RT requires MP)
  //
  if (!MP) return 0;

  spdlog::info("calc_microphysics_dU() Updating MicroPhysics");
#ifndef NDEBUG
  spdlog::debug("  RT-Nsrc={}", SimPM.RS.Nsources);
#endif  // NDEBUG
  int err = 0;

  if (SimPM.RS.Nsources == 0 || !SimPM.EP.raytracing) {
    //
    // If no radiative transfer, then just do a simple MP update.
    //
#ifdef RT_TESTING
    spdlog::info("\t\t--- calling calc_noRT_microphysics_dU()");
#endif  // RT_TESTING
    err += calc_noRT_microphysics_dU(delt, grid);
  }

  else {
    //
    // must have at least one radiation source, so we call a newer
    // update function with a bit more overhead.
    //
#ifdef RT_TESTING
    spdlog::info("\t\t--- calling calc_RT_microphysics_dU()");
#endif  // RT_TESTING
    err += calc_RT_microphysics_dU(delt, grid);
  }

  // cout <<"\tcalc_microphysics_dU finished.\n";
  if (0 != err)
    spdlog::error("{}: Expected {} but got {}", "calc_microphysics_dU", 0, err);
  return err;
}



// ##################################################################
// ##################################################################



int time_integrator::calc_RT_microphysics_dU(
    const double delt,         // timestep to integrate
    class GridBaseClass *grid  ///< Computational grid.
)
{
#ifdef RT_TESTING
  if (!grid->RT)
    spdlog::error(
        "{}: {}", "Logic error: must have RT unless i'm an idiot",
        "GENERAL-RT");
  if (SimPM.RS.Nsources < 1)
    spdlog::error(
        "{}: {}", "Need at least one source for diffuse-RT update",
        SimPM.RS.Nsources);
#endif  // RT_TESTING

  int err = 0;
  cell *c = grid->FirstPt_All();
  std::vector<pion_flt> p(
      SimPM.nvar);  // temporary state vector for output state.
  std::vector<pion_flt> ui(SimPM.nvar),
      uf(SimPM.nvar);  // conserved variable states.
  double tt = 0.;      // temperature returned at end of microphysics step.
  int index[3];
  // copies of the source data (one copy for each thread)
  // seems to be memory issues when using class member data.
  vector<struct rt_source_data> heating;
  vector<struct rt_source_data> ionize;

#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,index,p,ui,uf,tt,heating,ionize,err)
#endif
    // loop through cells in 1st y-z plane and calculate the MP update for the
    // x-column of cells associated with each starting cell.
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        heating  = FVI_heating_srcs;
        ionize   = FVI_ionising_srcs;
        do {

          //
          // Check if cell is internal boundary data or not.  If it is
          // boundary data, then we don't want to update anything, so we skip it
          //
          if (!c->isdomain) {
#ifndef NDEBUG
            spdlog::debug(
                "skipping cell {} in calc_RT_microphysics_dU() c->isdomain.\n",
                c->id);
#endif
          }
          else {
            //
            // Get column densities and Vshell in struct for each source.
            //
            for (int v = 0; v < FVI_nheat; v++) {
              heating[v].Vshell = CI.get_cell_Vshell(c, heating[v].id);
              heating[v].dS     = CI.get_cell_deltaS(c, heating[v].id);
              CI.get_cell_col(c, heating[v].id, heating[v].DelCol);
              CI.get_col(c, heating[v].id, heating[v].Column);
              for (short unsigned int iC = 0; iC < heating[v].NTau; iC++)
                heating[v].Column[iC] -= heating[v].DelCol[iC];
            }
            for (int v = 0; v < FVI_nion; v++) {
              ionize[v].Vshell = CI.get_cell_Vshell(c, ionize[v].id);
              ionize[v].dS     = CI.get_cell_deltaS(c, ionize[v].id);
              CI.get_cell_col(c, ionize[v].id, ionize[v].DelCol);
              CI.get_col(c, ionize[v].id, ionize[v].Column);
              for (short unsigned int iC = 0; iC < ionize[v].NTau; iC++)
                ionize[v].Column[iC] -= ionize[v].DelCol[iC];
              if (ionize[v].Column[0] < 0.0) {
#ifdef RT_TESTING
                spdlog::debug("dx={}", grid->DX());
                CI.print_cell(c);
                spdlog::error(
                    "{}: {}", "time_int:calc_RT_microphysics_dU tau<0", 1);
#endif
                for (short unsigned int iC = 0; iC < ionize[v].NTau; iC++)
                  ionize[v].Column[iC] = max(0.0, ionize[v].Column[iC]);
              }
            }
#ifdef RT_TESTING
            for (int v = 3; v < FVI_nion; v++) {
              spdlog::debug(
                  "{}  {}  {}  {}  {}  {}  {}  {}  {}", v, ionize[v].Vshell,
                  ionize[v].dS, ionize[v].strength[0], ionize[v].id,
                  ionize[v].type, ionize[v].NTau, ionize[v].Column[0],
                  ionize[v].DelCol[0]);
            }
#endif

            // 4th and 5th args are for ionising sources.
            err += MP->TimeUpdateMP_RTnew(
                c->P, FVI_nheat, heating, FVI_nion, ionize, p.data(), delt,
                SimPM.gamma, 0, &tt);

            // New state is p[], old state is c->P[].  Get dU from these.
            spatial_solver->PtoU(c->P, ui.data(), SimPM.gamma);
            spatial_solver->PtoU(p.data(), uf.data(), SimPM.gamma);
            for (int v = 0; v < SimPM.nvar; v++)
              c->dU[v] += uf[v] - ui[v];
#ifdef TEST_INF
            for (int v = 0; v < SimPM.nvar; v++) {
              if (!isfinite(c->P[v]) || !isfinite(p[v])
                  || !isfinite(c->dU[v])) {
                spdlog::error("NAN/INF in calc_RT_microphysics_dU() ");
                for (int s = 0; s < FVI_nion; s++) {
                  for (short unsigned int iC = 0;
                       iC < FVI_ionising_srcs[s].NTau; iC++) {
                    spdlog::debug(
                        "ion: Tau = {}, dTau = {}",
                        FVI_ionising_srcs[s].Column[iC],
                        FVI_ionising_srcs[s].DelCol[iC]);
                  }
                }
                spdlog::debug("Pin  : {}", c->P);
                spdlog::debug("Pout : {}", p);
                spdlog::debug("dU : {}", c->dU);
                CI.print_cell(c);
                spdlog::error(
                    "{}: {}", "NAN/INF in calc_RT_microphysics_dU()", v);
              }
            }
#endif
          }                                        // if not boundary data.
        } while ((c = grid->NextPt(c, XP)) != 0);  // loop over x-column cells
      }                                            // ax2
    }                                              // ax3
#ifdef PION_OMP
  }
#endif
  //    cout <<"calc_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
}  // RT microphysics update.



// ##################################################################
// ##################################################################



int time_integrator::calc_noRT_microphysics_dU(
    const double delt,         ///< timestep to integrate
    class GridBaseClass *grid  ///< Computational grid.
)
{
#ifndef NDEBUG
  spdlog::info("calc_noRT_microphysics_dU starting");
#endif
  //
  // No radiation sources and no diffuse radiation optical depths,
  // so call a simple microphysics update.
  //
  cell *c = grid->FirstPt_All();
  std::vector<pion_flt> p(
      SimPM.nvar);  // temporary state vector for output state.
  std::vector<pion_flt> ui(SimPM.nvar),
      uf(SimPM.nvar);  // conserved variable states.
  double tt = 0.;      // temperature returned at end of microphysics step.
  int err   = 0;
  int index[3];
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,index,p,ui,uf,tt,err)
#endif
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        do {
          // Check if cell is internal boundary data or not.  If it is
          // boundary data, then we don't want to update anything, so we skip it
          if (!c->isdomain) {
#ifndef NDEBUG
            spdlog::debug(
                "skipping cell {} in calc_noRT_microphysics_dU() c->isdomain.",
                c->id);
#endif
          }
          else {
            // integer 5th argument determines type of integration substepping:
            // 0 = adaptive RK5 Cash-Karp method.
            // 1 = adaptive euler integration.
            // 2 = single step RK4 method (at your own risk!)
            err += MP->TimeUpdateMP(c->P, p.data(), delt, SimPM.gamma, 0, &tt);
            if (err)
              spdlog::error(
                  "{}: {}", "calc_noRT_microphysics_dU returned error: cell id",
                  c->id);
            // New state is p[], old state is c->P[].  Get dU from these.
            spatial_solver->PtoU(c->P, ui.data(), SimPM.gamma);
            spatial_solver->PtoU(p.data(), uf.data(), SimPM.gamma);
            for (int v = 0; v < SimPM.nvar; v++)
              c->dU[v] += uf[v] - ui[v];
          }                                        // if not boundary data.
        } while ((c = grid->NextPt(c, XP)) != 0);  // loop over x-column cells
      }                                            // ax2
    }                                              // ax3
#ifdef PION_OMP
  }
#endif
  //    cout <<"calc_noRT_microphysics_dU() Updating MicroPhysics Done!\n";
  return err;
}



// ##################################################################
// ##################################################################



int time_integrator::calc_dynamics_dU(
    const double dt,           ///< timestep to integrate
    const int step,            ///< whether OA1 or OA2.
    class GridBaseClass *grid  ///< Computational grid.
)
{
  // cout <<"\tcalc_dynamics_dU starting.\n";
  //
  // first check if we are doing dynamics, and return if not.
  //
  if (!SimPM.EP.dynamics) return 0;
  int err = 0;

#ifndef NDEBUG
  if (step == OA1)
    spdlog::info("*****Updating dynamics: OA1");
  else if (step == OA2)
    spdlog::info("*****Updating dynamics: OA2");
  else
    spdlog::error("{}: {}", "Bad ooa", step);
#endif  // NDEBUG

  //
  // First we pre-process the cells, if needed.  This is required for
  // multi-dimensional viscosity such as the H-Correction.
  //
  err = preprocess_data(step, SimPM, grid);

  //
  // Now calculate the directionally-unsplit time update for the
  // conserved variables:
  //
  err = set_dynamics_dU(dt, step, grid);  //,time_ooa);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "calc_dynamics_dU() set_dynamics_dU returned error.", 0, err);

  //
  // This function is used for refined grids, to make the flux across
  // grid boundaries be consistent across all levels.
  //
  // Other potential uses of post-processing include if we are doing
  // something like Constrained Transport, where we have to change
  // the B-field update.
  //
  err = spatial_solver->PostProcess_dU(dt, step, SimPM, grid);
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}",
        "calc_dynamics_dU() spatial_solver->PostProcess_dU()", 0, err);

  return 0;
}

// ##################################################################
// ##################################################################

//
// Multi-dimensional calculations to be performed on every cell before
// the flux calculations.
//
int time_integrator::preprocess_data(
    const int csp,             ///< order of accuracy required.
    class SimParams &SimPM,    ///< pointer to simulation parameters
    class GridBaseClass *grid  ///< pointer to grid.
)
{
  //  cout <<"\t\t\tpreprocess_data(): Starting: ndim = "<<SimPM.ndim<<"\n";
  int err = 0;

  // ----------------------------------------------------------------
  // For the H-correction we need a maximum speed in each direction
  // for every cell.
  if (SimPM.artviscosity == AV_HCORRECTION
      || SimPM.artviscosity == AV_HCORR_FKJ98) {
    err += calc_Hcorrection(csp, SimPM, grid);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // HLLD has a switch based on velocity divergence, where it can
  // reduce to HLL near strong shocks.  So set divV here.
  if (SimPM.solverType == FLUX_RS_HLLD) {
    class cell *c = grid->FirstPt_All();
    double gradp  = 0.0;
    int indices[MAX_DIM];
    indices[0] = VX;
    indices[1] = VY;
    indices[2] = VZ;

    int index[3];
#ifdef PION_OMP
    #pragma omp parallel
    {
      #pragma omp for collapse(2) private(c,index)
#endif
      for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
        for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
          index[0] = 0;
          index[1] = ax2;
          index[2] = ax3;
          c        = grid->get_cell_all(index[0], index[1], index[2]);
          do {
            CI.set_DivV(c, spatial_solver->Divergence(c, 1, indices, grid));
            gradp = 0.0;
            for (int i = 0; i < SimPM.ndim; i++)
              gradp += spatial_solver->GradZone(grid, c, i, 1, PG);
            CI.set_MagGradP(c, gradp);
          } while ((c = grid->NextPt(c, XP)) != 0);
        }  // ax2
      }    // ax3
#ifdef PION_OMP
    }
#endif
  }
  // ----------------------------------------------------------------

  return err;
}

// ##################################################################
// ##################################################################

int time_integrator::calc_Hcorrection(
    const int csp,
    class SimParams &SimPM,  ///< pointer to simulation parameters
    class GridBaseClass *grid)
{
#ifndef NDEBUG
  spdlog::debug("\t\t\tcalc_Hcorrection() ndim = {}", SimPM.ndim);
#endif  // NDEBUG

  // Allocate arrays for direction values.
  int err = 0;
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP;
  posdirs[1] = YP;
  posdirs[2] = ZP;
  negdirs[0] = XN;
  negdirs[1] = YN;
  negdirs[2] = ZN;
  axis[0]    = XX;
  axis[1]    = YY;
  axis[2]    = ZZ;

  // Loop through each dimension.
  for (int idim = 0; idim < SimPM.ndim; idim++) {
#ifndef NDEBUG
    spdlog::debug("\t\t\tidim={}", idim);
#endif  // NDEBUG

#ifdef PION_OMP
    #pragma omp parallel
    {
      spatial_solver->SetDirection(axis[idim]);
    }
#else
    spatial_solver->SetDirection(axis[idim]);
#endif  // PION_OMP

    class cell *cpt = grid->FirstPt_All();

#ifdef TEST_INT
    spdlog::debug("Direction={}, i={}", axis[idim], idim);
    spdlog::debug("cpt : {}", cpt->pos);
#endif

    // loop over the number of cells in the line/plane of starting
    // cells.
    enum axes x1 = axis[(idim + 1) % 3];
    enum axes x2 = axis[(idim + 2) % 3];
    enum axes x3 = axis[idim];

    // loop over the two perpendicular axes, to trace out a plane of
    // starting cells for calculating fluxes along columns along this
    // axis.  The plane can be a single cell
    // (in 1D) or a line (in 2D) or a plane (in 3D).
    int index[3];
#ifdef PION_OMP
    #pragma omp parallel
    {
      #pragma omp for collapse(2) private(index,cpt)
#endif
      for (int ax2 = 0; ax2 < grid->NG_All(x2); ax2++) {
        for (int ax1 = 0; ax1 < grid->NG_All(x1); ax1++) {
          index[x1] = static_cast<int>(ax1);
          index[x2] = static_cast<int>(ax2);
          index[x3] = 0;
          cpt       = grid->get_cell_all(index[0], index[1], index[2]);
          // Slope and edge state temporary arrays
          pion_flt *slope_cpt = 0, *slope_npt = 0, *temp = 0;
          slope_cpt = mem.myalloc(slope_cpt, SimPM.nvar);
          slope_npt = mem.myalloc(slope_npt, SimPM.nvar);
          std::vector<pion_flt> edgeR(SimPM.nvar), edgeL(SimPM.nvar);

          // Set three cell pointers (2nd order slopes have a 3-point
          // stencil).
          cell *npt  = grid->NextPt(cpt, posdirs[idim]);
          cell *n2pt = grid->NextPt(npt, posdirs[idim]);
          if (npt == 0 || n2pt == 0)
            spdlog::error("{}: {}", "Couldn't find three cells in column", 0);
          // Need to get slopes and edge states if 2nd order (csp==OA2).
          for (int v = 0; v < SimPM.nvar; v++)
            slope_cpt[v] = 0.;

          // --------------------------------------------------------
          // Run through column, calculating slopes, edge-states, and
          // eta[] values as we go.
          // --------------------------------------------------------
          do {
            err += spatial_solver->SetEdgeState(
                cpt, posdirs[idim], SimPM.nvar, slope_cpt, edgeL.data(), csp,
                grid);
            err += spatial_solver->SetSlope(
                npt, axis[idim], SimPM.nvar, slope_npt, csp, grid);
            err += spatial_solver->SetEdgeState(
                npt, negdirs[idim], SimPM.nvar, slope_npt, edgeR.data(), csp,
                grid);
            spatial_solver->set_Hcorrection(
                cpt, axis[idim], edgeL.data(), edgeR.data(), SimPM.gamma);

            cpt       = npt;
            npt       = n2pt;
            temp      = slope_cpt;
            slope_cpt = slope_npt;
            slope_npt = temp;
          } while ((n2pt = grid->NextPt(n2pt, posdirs[idim])) != 0);

          // last cell must be 1st order.
          err += spatial_solver->SetEdgeState(
              cpt, posdirs[idim], SimPM.nvar, slope_cpt, edgeL.data(), csp,
              grid);
          for (int v = 0; v < SimPM.nvar; v++)
            slope_npt[v] = 0.;
          err += spatial_solver->SetEdgeState(
              npt, negdirs[idim], SimPM.nvar, slope_npt, edgeR.data(), csp,
              grid);
          spatial_solver->set_Hcorrection(
              cpt, axis[idim], edgeL.data(), edgeR.data(), SimPM.gamma);

          // --------------------------------------------------------
          // Finished H-correction calculation for the column.
          // --------------------------------------------------------
          slope_cpt = mem.myfree(slope_cpt);
          slope_npt = mem.myfree(slope_npt);
        }
      }
#ifdef PION_OMP
    }
#endif  // PION_OMP
  }     // Loop over three directions.

#ifdef PION_OMP
  #pragma omp parallel
  {
    spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
  }
#else
  spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
#endif  // PION_OMP

  return err;
}  // calc_Hcorrection()



// ##################################################################
// ##################################################################



int time_integrator::set_dynamics_dU(
    const double dt,           ///< timestep for this calculation
    const int step,            ///< whether half-step or full-step
    class GridBaseClass *grid  ///< Computational grid.
)
{
  //
  // Allocate arrays
  //
  int return_value = 0;
  enum direction posdirs[MAX_DIM], negdirs[MAX_DIM];
  enum axes axis[MAX_DIM];
  posdirs[0] = XP;
  posdirs[1] = YP;
  posdirs[2] = ZP;
  negdirs[0] = XN;
  negdirs[1] = YN;
  negdirs[2] = ZN;
  axis[0]    = XX;
  axis[1]    = YY;
  axis[2]    = ZZ;

  int space_ooa;
  if (step == OA1)
    space_ooa = OA1;
  else
    space_ooa = OA2;

#ifdef TEST_INT
  // get current level of grid in hierarchy.
  int level = 0;
  if (SimPM.grid_nlevels > 1) {
    for (int v = 0; v < SimPM.grid_nlevels; v++) {
      if (grid == SimPM.levels[v].grid) level = v;
    }
  }
  spdlog::debug("*** Calculating DU dynamics on level {}", level);
#endif

  //
  // Loop over all directions, and in each direction, calculate fluxes
  // in all columns of cells in that direction.
  //
  for (int i = 0; i < SimPM.ndim; i++) {
#ifdef PION_OMP
    #pragma omp parallel
    {
      spatial_solver->SetDirection(axis[i]);
    }
#else
    spatial_solver->SetDirection(axis[i]);
#endif  // PION_OMP
    class cell *cpt = grid->FirstPt_All();

#ifdef TEST_INT
    spdlog::debug("Direction={}, i={}", axis[i], i);
    spdlog::debug("cpt : {}", cpt->pos);
#endif

    //
    // loop over the number of cells in the line/plane of starting
    // cells.
    //
    enum axes x1 = axis[(i + 1) % 3];
    enum axes x2 = axis[(i + 2) % 3];
    enum axes x3 = axis[i];

    //
    // loop over the two perpendicular axes, to trace out a plane of
    // starting cells for calculating fluxes along columns along this
    // axis.  Note that the NG_All() array is initialised so that
    // unused dimensions have NG=1, so the plane can be a single cell
    // (in 1D) or a line (in 2D) or a plane (in 3D).
    //
    int index[3];
#ifdef PION_OMP
    #pragma omp parallel
    {
      #pragma omp for collapse(2) private(index,cpt)
#endif
      for (int ax2 = 0; ax2 < grid->NG_All(x2); ax2++) {
        for (int ax1 = 0; ax1 < grid->NG_All(x1); ax1++) {
          index[x1]    = static_cast<int>(ax1);
          index[x2]    = static_cast<int>(ax2);
          index[x3]    = 0;
          cpt          = grid->get_cell_all(index[0], index[1], index[2]);
          return_value = dynamics_dU_column(
              cpt, posdirs[i], negdirs[i], dt, space_ooa, grid);
          if (0 != return_value)
            spdlog::error(
                "{}: Expected {} but got {}", "set_dynamics_dU: column", 0,
                return_value);
        }
      }
#ifdef PION_OMP
    }
#endif  // PION_OMP

  }  // Loop over three directions.

#ifdef PION_OMP
  #pragma omp parallel
  {
    spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
  }
#else
  spatial_solver->SetDirection(axis[0]);  // Reset fluxes to x-dir.
#endif  // PION_OMP
  return 0;
}  // set_dynamics_dU()



// ##################################################################
// ##################################################################



int time_integrator::dynamics_dU_column(
    class cell *startingPt,       ///< starting point of column.
    const enum direction posdir,  ///< direction to trace column.
    const enum direction negdir,  ///< reverse direction
    const double dt,              ///< timestep we are advancing by.
    const int csp,                ///< spatial order-of-accuracy for this step.
    class GridBaseClass *grid     ///< Computational grid.
)
{
  if ((SimPM.spOOA > 2) || (SimPM.tmOOA > 2) || (csp > 2)) {
    spdlog::error(
        "(dynamics_dU_column) Error, ooa={}, {}", SimPM.spOOA, SimPM.tmOOA);
    return (1);
  }
  // cout <<"dynamics_dU_column: d+="<<posdir<<", d-="<<negdir;
  // cout <<", csp="<<csp<<", OOA="<<SimPM.spOOA<<"\n";
  int err        = 0;
  enum axes axis = spatial_solver->GetDirection();
  double dx      = grid->DX();
#ifdef TEST_CONSERVATION
  double dE = 0.0, dMX = 0.0, dMY = 0.0, dMZ = 0.0, dM = 0.0;
#endif

  // Calculate Flux at positive (right) boundary of cell for the
  // current cell (Fr_this) and the negative neighbour (Fr_prev).
  pion_flt *Fr_this = 0, *Fr_prev = 0, *temp = 0, *slope_cpt = 0,
           *slope_npt = 0, *edgeR = 0, *edgeL = 0, *pstar = 0;
  Fr_prev   = mem.myalloc(Fr_prev, SimPM.nvar);
  Fr_this   = mem.myalloc(Fr_this, SimPM.nvar);
  slope_cpt = mem.myalloc(slope_cpt, SimPM.nvar);
  slope_npt = mem.myalloc(slope_npt, SimPM.nvar);
  edgeL     = mem.myalloc(edgeL, SimPM.nvar);
  edgeR     = mem.myalloc(edgeR, SimPM.nvar);
  pstar     = mem.myalloc(pstar, SimPM.nvar);

  //
  // Set starting point, and next two points in the column.
  //
  cell *cpt = startingPt;
  if (cpt == 0) {
    spdlog::error("(dynamics_dU_column) error finding left boundary cell");
    return (1);
  }
  cell *npt  = grid->NextPt(cpt, posdir);
  cell *n2pt = grid->NextPt(npt, posdir);
#ifdef TEST_INT
  // cout <<"Column: "<<cpt<<", "<<npt<<", "<<n2pt<<"\n";
#endif
  if (npt == 0 || n2pt == 0)
    spdlog::error("{}: {}", "Couldn't find two real cells in column", 0);

  //
  // Left Ghost Cell (doesn't get updated)
  //
  for (int v = 0; v < SimPM.nvar; v++) {
    slope_cpt[v] = 0.;
    slope_npt[v] = 0.;
    Fr_prev[v]   = 0.;
    Fr_this[v]   = 0.;
    edgeL[v]     = 0.;
    edgeR[v]     = 0.;
  }

  //
  // Now go through all cells in the column and calculate fluxes and add to dU
  // vector.
  //
  do {
#ifndef NDEBUG
    dp.c = cpt;
#endif
#ifdef TEST_INT
    cout << "First Cell:";
    CI.print_cell(cpt);
    cout << "Next Cell: ";
    CI.print_cell(npt);
#endif
    // Get the flux from left and right states, adding artificial
    // viscosity if needed.
    err += spatial_solver->SetEdgeState(
        cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
    err +=
        spatial_solver->SetSlope(npt, axis, SimPM.nvar, slope_npt, csp, grid);
#ifdef ZERO_SLOPE_TRACERS
    // not recommended b/c it is diffusive, but it does make the code
    // more symmetric.
    for (int v = SimPM.ftr; v < SimPM.nvar; v++)
      slope_npt[v] = 0.0;
#endif
    err += spatial_solver->SetEdgeState(
        npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
    err += spatial_solver->InterCellFlux(
        SimPM, grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.gamma, dx);
    err += spatial_solver->MHDsource(
        grid, cpt, npt, edgeL, edgeR, axis, posdir, negdir, dt);
    err += spatial_solver->dU_Cell(
        grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, dx, dt);

    // record flux entering and leaving domain
    if (cpt->isbd_ref[negdir]) {
      for (int v = 0; v < SimPM.nvar; v++)
        cpt->F[axis][v] = Fr_this[v];
    }
    if (npt->isbd_ref[posdir]) {
      for (int v = 0; v < SimPM.nvar; v++)
        npt->F[axis][v] = Fr_this[v];
    }

#ifdef TEST_INT
    for (int v = 0; v < SimPM.nvar; v++) {
      if (!isfinite(cpt->dU[v])) {
        spdlog::debug("Fl : {}", Fr_prev);
        spdlog::debug("Fr : {}", Fr_this);
        spdlog::debug("El : {}", edgeL);
        spdlog::debug("Er : {}", edgeR);
        spdlog::debug("dt:{}\tdx={}", dt, dx);
        //  rep.printVec("dU",&cpt->dU[v],1);
        CI.print_cell(cpt);
        CI.print_cell(npt);
        spdlog::error("{}: {}", "nans!!!", 2);
      }
    }
#endif

#ifdef TEST_CONSERVATION
    // Track energy, momentum entering/leaving domain, if outside
    // boundary
    double dA = spatial_solver->CellInterface(cpt, posdir, dx);
    if (csp == SimPM.tmOOA && pconst.equalD(grid->Xmin(axis), SimPM.Xmin[axis])
        && !(cpt->isdomain) && npt->isgd && npt->isleaf) {
      dM += Fr_this[RHO] * dt * dA;
      dE += Fr_this[ERG] * dt * dA;
      dMX += Fr_this[MMX] * dt * dA;
      dMY += Fr_this[MMY] * dt * dA;
      dMZ += Fr_this[MMZ] * dt * dA;
    }
    else if (
        csp == SimPM.tmOOA && pconst.equalD(grid->Xmax(axis), SimPM.Xmax[axis])
        && cpt->isgd && cpt->isleaf && !npt->isdomain) {
      // cout <<"leaving domain\n";
      dM -= Fr_this[RHO] * dt * dA;
      dE -= Fr_this[ERG] * dt * dA;
      dMX -= Fr_this[MMX] * dt * dA;
      dMY -= Fr_this[MMY] * dt * dA;
      dMZ -= Fr_this[MMZ] * dt * dA;
    }
#endif  // TEST_CONSERVATION
        //
        // Now move temp arrays for moving on to the next cell
        //
    temp      = Fr_prev;
    Fr_prev   = Fr_this;
    Fr_this   = temp;
    temp      = slope_cpt;
    slope_cpt = slope_npt;
    slope_npt = temp;
    cpt       = npt;
    npt       = n2pt;
  } while ((n2pt = grid->NextPt(n2pt, posdir)));

  //
  // Now n2pt=null. npt = bd-data, cpt= (gd/bd)-data.
  // So have to do something different.
  //
#ifndef NDEBUG
  dp.c = cpt;
#endif
  // last cell 1st order.
  err += spatial_solver->SetEdgeState(
      cpt, posdir, SimPM.nvar, slope_cpt, edgeL, csp, grid);
  for (int v = 0; v < SimPM.nvar; v++)
    slope_npt[v] = 0.;
  err += spatial_solver->SetEdgeState(
      npt, negdir, SimPM.nvar, slope_npt, edgeR, csp, grid);
  err += spatial_solver->InterCellFlux(
      SimPM, grid, cpt, npt, edgeL, edgeR, Fr_this, SimPM.gamma, dx);
  err += spatial_solver->MHDsource(
      grid, cpt, npt, edgeL, edgeR, axis, posdir, negdir, dt);
  err += spatial_solver->dU_Cell(
      grid, cpt, axis, Fr_prev, Fr_this, slope_cpt, csp, dx, dt);
  // record flux entering and leaving domain
  if (cpt->isbd_ref[negdir]) {
    for (int v = 0; v < SimPM.nvar; v++)
      cpt->F[axis][v] = Fr_this[v];
  }
  if (npt->isbd_ref[posdir]) {
    for (int v = 0; v < SimPM.nvar; v++)
      npt->F[axis][v] = Fr_this[v];
  }

#ifdef TEST_CONSERVATION
  // Track energy, momentum entering/leaving domain, if outside
  // boundary
  double dA = spatial_solver->CellInterface(cpt, posdir, dx);
  if (csp == SimPM.tmOOA && pconst.equalD(grid->Xmax(axis), SimPM.Xmax[axis])
      && cpt->isgd && cpt->isleaf && !npt->isdomain) {
    dM -= Fr_this[RHO] * dt * dA;
    dE -= Fr_this[ERG] * dt * dA;
    dMX -= Fr_this[MMX] * dt * dA;
    dMY -= Fr_this[MMY] * dt * dA;
    dMZ -= Fr_this[MMZ] * dt * dA;
  }
#endif  // TEST_CONSERVATION

  //
  // Right Ghost Cell-- have to calculate it's left interface differently,
  //
  cpt = npt;
#ifndef NDEBUG
  dp.c = cpt;
#endif
  for (int v = 0; v < SimPM.nvar; v++)
    cpt->dU[v] += 0.;  // nothing to calculate for it.
  if (0 != err)
    spdlog::error(
        "{}: Expected {} but got {}", "(dU_Column) encountered an error", 0,
        err);

  Fr_this   = mem.myfree(Fr_this);
  Fr_prev   = mem.myfree(Fr_prev);
  slope_cpt = mem.myfree(slope_cpt);
  slope_npt = mem.myfree(slope_npt);
  edgeL     = mem.myfree(edgeL);
  edgeR     = mem.myfree(edgeR);
  pstar     = mem.myfree(pstar);

#ifdef TEST_CONSERVATION
#ifdef PARALLEL
  dM  = SimPM.levels[0].sub_domain.global_operation_double("SUM", dM);
  dE  = SimPM.levels[0].sub_domain.global_operation_double("SUM", dE);
  dMX = SimPM.levels[0].sub_domain.global_operation_double("SUM", dMX);
  dMY = SimPM.levels[0].sub_domain.global_operation_double("SUM", dMY);
  dMZ = SimPM.levels[0].sub_domain.global_operation_double("SUM", dMZ);
  // cout <<"d="<<dM<<", "<<dE<<", "<<dMX<<", "<<dMY<<"\n";
#endif
  initMASS += dM;
  initERG += dE;
  initMMX += dMX;
  initMMY += dMY;
  initMMZ += dMZ;
#endif  // TEST_CONSERVATION

  return 0;
}



// ##################################################################
// ##################################################################



int time_integrator::grid_update_state_vector(
    const double dt,           ///< timestep
    const int step,            ///< OA1 or OA2
    const int ooa,             ///< Full order of accuracy of simulation
    class GridBaseClass *grid  ///< Computational grid.
)
{
  int err = 0;
  // temp variable to handle change of energy when correcting for negative
  // pressure.
  pion_flt temperg = 0.0;

  // Loop through grid, updating Ph[] with CellAdvanceTime function.
  class cell *c = grid->FirstPt_All();
  int index[3];
  double T = 0.0;
#ifdef PION_OMP
  #pragma omp parallel
  {
    #pragma omp for collapse(2) private(c,index,err,temperg,T)
#endif
    // loop through cells in 1st y-z plane and calculate the MP update for the
    // x-column of cells associated with each starting cell.
    for (int ax3 = 0; ax3 < grid->NG_All(ZZ); ax3++) {
      for (int ax2 = 0; ax2 < grid->NG_All(YY); ax2++) {
        index[0] = 0;
        index[1] = ax2;
        index[2] = ax3;
        c        = grid->get_cell_all(index[0], index[1], index[2]);
        do {

#ifndef NDEBUG
          double dx       = grid->DX();
          dp.ergTotChange = 0.;
          temperg         = 0.;
          dp.c            = c;
#endif
          if (!c->isdomain || !c->isleaf) {
            // skip cell if it has been cut out of the domain.
            for (int v = 0; v < SimPM.nvar; v++)
              c->dU[v] = 0.0;
          }
          else {
            err += spatial_solver->CellAdvanceTime(
                c, c->P, c->dU, c->Ph, &temperg, SimPM.gamma,
                SimPM.EP.MinTemperature, dt);
          }

#ifndef NDEBUG
          if (err) {
            spdlog::error("______ Error in Cell-advance-time: ");
            CI.print_cell(c);
            CI.print_cell(c->npt);
            err = 0;
          }
#else
        // ignore negative pressures and try to continue
        if (err) err = 0;
#endif  // NDEBUG

          if (MP) {
            T = MP->Temperature(c->Ph, SimPM.gamma);
            if (T > SimPM.EP.MaxTemperature) {
              // cout <<"warning, temperature too large: "<<T<<"\n";
              MP->Set_Temp(c->Ph, SimPM.EP.MaxTemperature, SimPM.gamma);
            }
          }

          // If the current step is the full update, then also set
          // the updated base-state-vector to updated value.
          if (step == ooa) {
            for (int v = 0; v < SimPM.nvar; v++)
              c->P[v] = c->Ph[v];

#ifndef NDEBUG
            // Update Total Energy from fixing negative pressures. Reset
            // update variables.
            dp.ergTotChange = temperg;
            dp.initERG += dp.ergTotChange * spatial_solver->CellVolume(c, dx);
#endif  // NDEBUG
          }
        } while ((c = grid->NextPt(c, XP)) != 0);
      }  // ax2
    }    // ax3
#ifdef PION_OMP
  }
#endif


#ifndef NDEBUG
  spdlog::debug("\tgrid_update_state_vector done. error={}", err);
#endif  // NDEBUG
  return err;
}



// ##################################################################
// ##################################################################
