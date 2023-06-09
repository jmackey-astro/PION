/// \file stellar_wind_boundaries.cpp
/// \brief Class definitions for stellar_wind boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include "boundaries/stellar_wind_boundaries.h"
#include "grid/stellar_wind_BC.h"
#include "grid/stellar_wind_angle.h"
#include "grid/stellar_wind_latdep.h"
#include "sim_params.h"
#include "tools/mem_manage.h"

#include <array>
#include <spdlog/spdlog.h>
#include <sstream>
#include <vector>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

using namespace std;

// ##################################################################
// ##################################################################
struct starpos {
  int id;
  double mass;
  double pos[MAX_DIM];
  double vel[MAX_DIM];
  double acc[MAX_DIM];
};


stellar_wind_bc::stellar_wind_bc()
{
  return;
}

stellar_wind_bc::~stellar_wind_bc()
{
  if (outf.is_open()) {
    outf.flush();
    outf.close();
  }
  return;
}

// ##################################################################
// ##################################################################

//
// Add internal stellar wind boundaries -- these are (possibly
// time-varying) winds defined by a mass-loss-rate and a terminal
// velocity.  A region within the domain is given fixed values
// corresponding to a freely expanding wind from a
// cell-vertex-located source.
//
int stellar_wind_bc::BC_assign_STWIND(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    boundary_data *b,
    class microphysics_base *mp  ///< pointer to microphysics
)
{
  //
  // Check that we have an internal boundary struct, and that we have
  // a stellar wind source to set up.
  //
  if (b->dir != NO) {
    spdlog::error("STWIND not external boundary {}", static_cast<int>(b->dir));
    exit(1);
  }
#ifndef NDEBUG
  spdlog::debug("Assigning data to STWIND boundary. Nsrc={}", SWP.Nsources);
#endif
  if (SWP.Nsources < 1) {
    spdlog::error("BC_assign_STWIND() No Sources {}", SWP.Nsources);
    exit(1);
  }

  //
  // Setup reference state vector and initialise to zero.
  //
  if (b->refval) {
    spdlog::error("Initialised STWIND boundary refval {}", fmt::ptr(b->refval));
    exit(1);
  }
  b->refval = mem.myalloc(b->refval, par.nvar);
  if (!b->data.empty()) {
    spdlog::error("BC_assign_STWIND: Not empty boundary data {}", b->itype);
    exit(1);
  }
  for (int v = 0; v < par.nvar; v++)
    b->refval[v] = 0.0;

  //
  // New structure: we need to initialise the stellar wind class with
  // all of the wind sources in the global parameter list (this was
  // formerly done in DataIOBase::read_simulation_parameters()).
  //
  // The type of class we set up is determined first.
  // stellar_wind_evolution is derived from stellar_wind, and
  // stellar_wind_angle is derived from stellar_wind_evolution.
  //
  int err = 0, wtype = 0;
  int Ns = SWP.Nsources;
  for (int isw = 0; isw < Ns; isw++) {
    if (SWP.params[isw]->type == WINDTYPE_EVOLVING) wtype = 1;
  }
  for (int isw = 0; isw < Ns; isw++) {
    if (SWP.params[isw]->type == WINDTYPE_ANGLE) wtype = 2;
  }
  for (int isw = 0; isw < Ns; isw++) {
    if (SWP.params[isw]->type == WINDTYPE_LATDEP) wtype = 3;
  }

  //
  // check values of xi.  At the moment we assume it is the same for
  // all wind sources, so it is not source-dependent.
  //
  double xi = 0.0;
  for (int isw = 0; isw < Ns; isw++) {
    if (isw == 0)
      xi = SWP.params[isw]->xi;
    else {
      if (xi != SWP.params[isw]->xi) {
        spdlog::error(
            "wind xi values don't match {} {}", xi, SWP.params[isw]->xi);
        exit(1);
      }
    }
  }

  if (Ns > 0) {
    // cout <<"\n----------- SETTING UP STELLAR WIND CLASS ----------\n";
    if (wtype == 0) {
#ifndef NDEBUG
      spdlog::info("Setting up stellar_wind class");
#endif
      grid->Wind = new stellar_wind(
          par.ndim, par.nvar, par.ntracer, par.ftr, par.tracers, par.coord_sys,
          par.eqntype, par.EP.MinTemperature);
    }
    else if (wtype == 1) {
#ifndef NDEBUG
      spdlog::info("Setting up stellar_wind_evolution class");
#endif
      grid->Wind = new stellar_wind_evolution(
          par.ndim, par.nvar, par.ntracer, par.ftr, par.tracers, par.coord_sys,
          par.eqntype, par.EP.MinTemperature, par.starttime, par.finishtime);
      err = 0;
    }
    else if (wtype == 2) {
#ifndef NDEBUG
      spdlog::info("Setting up stellar_wind_angle class");
#endif
      grid->Wind = new stellar_wind_angle(
          par.ndim, par.nvar, par.ntracer, par.ftr, par.tracers, par.coord_sys,
          par.eqntype, par.EP.MinTemperature, par.starttime, par.finishtime,
          xi);
    }
    else if (wtype == 3) {
#ifndef NDEBUG
      spdlog::info("Setting up stellar_wind_angle class");
#endif
      grid->Wind = new stellar_wind_latdep(
          par.ndim, par.nvar, par.ntracer, par.ftr, par.tracers, par.coord_sys,
          par.eqntype, par.EP.MinTemperature, par.starttime, par.finishtime,
          xi);
    }
  }

  grid->Wind->SetMicrophysics(mp);

  // initialise "current_radius" var to value in parameter file.
  for (int isw = 0; isw < Ns; isw++) {
    SWP.params[isw]->current_radius = SWP.params[isw]->radius;
    for (int v = par.ndim; v < MAX_DIM; v++)
      SWP.params[isw]->dpos[v] = 0.0;
    for (int v = 0; v < MAX_DIM; v++)
      SWP.params[isw]->thispos[v] = SWP.params[isw]->dpos[v];
    spdlog::debug(
        "source {}, thispos {} dpos {}", isw, SWP.params[isw]->dpos,
        SWP.params[isw]->thispos);
  }

  // check acceleration flag.
  for (int isw = 0; isw < Ns; isw++) {
    spdlog::debug(
        "source {}, acc {}, wind-acc {}", isw, SWP.params[isw]->acc,
        par.EP.wind_acceleration);
    if (SWP.params[isw]->acc && !par.EP.wind_acceleration) {
      spdlog::error("requesting wind acceleration but par not set");
      exit(1);
    }
  }

  //
  // Run through sources and add sources.
  //
  for (int isw = 0; isw < Ns; isw++) {
#ifndef NDEBUG
    spdlog::debug("BC_assign_STWIND: Adding source {}", isw);
#endif
    if (SWP.params[isw]->type == WINDTYPE_CONSTANT) {
      //
      // This is for spherically symmetric winds that are constant
      // in time.
      //
      err = grid->Wind->add_source(SWP.params[isw]);
#ifndef NDEBUG
#ifdef ANALYTIC_ORBITS
      spdlog::debug(
          "add_source call: {},{}", SWP.params[isw]->PeriastronX,
          SWP.params[isw]->PeriastronY);
#endif
#endif
    }
    else if (SWP.params[isw]->evolving_wind_file == "NOFILE") {
      // star is not evolving, but also not type==0, so must be rotating
      spdlog::info("Adding rotating source {}", isw);
      err = grid->Wind->add_rotating_source(SWP.params[isw]);
    }
    else {
      //
      // This works for spherically symmetric winds and for
      // latitude-dependent winds that evolve over time.
      //
      // cout <<"Adding evolving source "<<isw<<" with filename ";
      // cout <<SWP.params[isw]->evolving_wind_file<<"\n";
      err = grid->Wind->add_evolving_source(par.simtime, SWP.params[isw]);
    }
    if (err) {
      spdlog::error("Error adding wind source {}", isw);
      exit(1);
    }
  }

  //
  // loop over sources, adding cells to boundary data list in order.
  //
  for (int isw = 0; isw < Ns; isw++) {
#ifndef NDEBUG
    spdlog::debug("BC_assign_STWIND: Adding cells to source {}", isw);
#endif
    for (int v = 0; v < MAX_DIM; v++)
      SWP.params[isw]->thispos[v] = SWP.params[isw]->dpos[v];
    grid->Wind->set_this_source_position(isw, SWP.params[isw]->dpos);

#ifndef NDEBUG
    spdlog::info(
        "setting wind radius: {:9.3e}  {:9.3e}",
        SWP.params[isw]->current_radius, SWP.params[isw]->radius);
#endif

    BC_set_wind_radius(par, grid, isw);

#ifndef NDEBUG
    spdlog::info(
        "setting wind radius: {:9.3e}  {:9.3e}",
        SWP.params[isw]->current_radius, SWP.params[isw]->radius);
#endif

    BC_assign_STWIND_add_cells2src(par, grid, isw);
  }

  // if any star is moving, then open a file to follow its trajectory
  bool moving = false;
  for (int isw = 0; isw < Ns; isw++) {
    if (SWP.params[isw]->moving_star == 1) moving = true;
  }
  if (moving && grid->level() == par.grid_nlevels - 1
#ifdef PARALLEL
      && par.levels[0].sub_domain.get_myrank() == 0
#endif
  ) {
    // spdlog::info("OPENING OUTPUT FILE trajectory.txt");
    outf.open("trajectory.txt", std::ios_base::app);
    outf.setf(ios_base::scientific);
    outf.precision(6);
    if (!outf.is_open()) {
      spdlog::error("Couldn't open file {}", "trajectory.txt");
      exit(1);
    }
    outf
        << "\n time         timestep       star 0 x       star 0 y       star 0 vx      star 0 vy      star 1 x       star 1 y        star 1 vx      star1 vy\n";
    outf.flush();
  }


  //
  // Now we should have set everything up, so we assign the boundary
  // cells with their boundary values.
  //
  err += BC_update_STWIND(par, 0, grid, par.simtime, 0.0, b, 0, 0);
#ifndef NDEBUG
  spdlog::info(
      "Finished setting up wind parameters\n------ DONE SETTING UP STELLAR WIND CLASS ----------\n");
#endif
  return err;
}



// ##################################################################
// ##################################################################



int stellar_wind_bc::BC_assign_STWIND_add_cells2src(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int id                ///< source id
)
{
  //
  // We run through each cell, and if it is within the
  // source's radius of influence, then we add it to the lists.
  //
  int err       = 0;
  int ncell     = 0;
  double srcrad = SWP.params[id]->current_radius;
  array<double, MAX_DIM> srcpos;
  grid->Wind->get_this_source_position(id, srcpos);
  // grid->Wind->get_src_posn(id, srcpos);
  // for (int v=0;v<MAX_DIM;v++) srcpos[v] = SWP.params[id]->thispos[v];

  spdlog::debug(
      "star id {}, lev {}, radius {:8.2e}, dx {:8.2e}", id, grid->level(),
      srcrad, grid->DX());

  if (grid->Wind->get_num_cells(id) != 0) {
    spdlog::error(
        "adding cells to source that already has cells: {}",
        grid->Wind->get_num_cells(id));
    exit(1);
  }

  // find min/max of coordinates of cube that contains spherical wind source
  double dx = grid->DX(), aneg = 0.0, apos = 0.0, xmin = 0.0, xmax = 0.0;
  array<int, MAX_DIM> ineg = {0, 0, 0}, ipos = {0, 0, 0};
  for (int v = 0; v < par.ndim; v++) {
    xmin    = grid->Xmin_all(static_cast<axes>(v));
    xmax    = grid->Xmax_all(static_cast<axes>(v));
    aneg    = srcpos[v] - srcrad - dx;  // add extra safety factor dx/2
    aneg    = max(xmin, aneg);
    aneg    = min(xmax, aneg);
    apos    = srcpos[v] + srcrad + dx;  // add extra safety factor dx/2
    apos    = min(xmax, apos);
    apos    = max(xmin, apos);
    ineg[v] = (aneg - xmin) / dx;
    ipos[v] = (apos - xmin) / dx;
  }
  for (int v = par.ndim; v < MAX_DIM; v++) {
    ineg[v] = 0;
    ipos[v] = 1;
  }

#ifndef NDEBUG
  spdlog::info("*** srcrad={}, ineg {}, ipos {}", srcrad, ineg, ipos);
#endif

  // loop over cube of data and add only cells with position within the
  // source radius
  enum axes x1 = XX;
  enum axes x2 = YY;
  enum axes x3 = ZZ;
  //#ifdef PION_OMP
  //  #pragma omp parallel for collapse(2)
  //#endif
  for (int ax3 = ineg[x3]; ax3 < ipos[x3]; ax3++) {
    for (int ax2 = ineg[x2]; ax2 < ipos[x2]; ax2++) {
      array<double, MAX_DIM> cpos;
      int index[3];
      index[x1] = ineg[x1];
      index[x2] = static_cast<int>(ax2);
      index[x3] = static_cast<int>(ax3);
      cell *c   = grid->get_cell_all(index[0], index[1], index[2]);
      for (int ax1 = index[x1]; ax1 < ipos[x1]; ax1++) {
        CI.get_dpos(*c, cpos);
        if (grid->distance(srcpos, cpos) <= srcrad) {
          ncell++;
          err += grid->Wind->add_cell(grid, id, *c);
          // spdlog::info("src id {}, add cell {}",id,c->id);
        }
        c = grid->NextPt(*c, XP);
      }
    }
  }

#ifndef NDEBUG
  spdlog::debug(
      "BC_assign_STWIND_add_cells2src: Added {} cells to wind boundary for WS {}",
      ncell, id);
#endif
  return err;
}



// ##################################################################
// ##################################################################



int stellar_wind_bc::BC_update_STWIND(
    class SimParams &par,       ///< pointer to simulation parameters
    const int l,                ///< level in grid hierarchy
    class GridBaseClass *grid,  ///< pointer to grid.
    const double simtime,       ///< current simulation time
    const double dt,            ///< timestep
    boundary_data *,            ///< Boundary to update (unused)
    const int,                  ///< current fractional step being taken
    const int                   ///< final step
)
{
  int err = 0;
//#define ANALYTIC_ORBITS
#ifdef ANALYTIC_ORBITS
  // Moving wind stuff
  // loop over wind-sources
  for (int id = 0; id < grid->Wind->Nsources(); id++) {
    array<double, MAX_DIM> srcpos;
    double srcrad;
    double eccentricity_fac = 0.0, periastron_vec[2], orbital_period = 0.0,
           initial_position[MAX_DIM];  // orbital parameters
    double newpos[MAX_DIM];  // new source-position that will be assigned
    // elipse parameters
    double cos_a = 0.0, sin_a = 0.0, a = 0.0, b = 0.0, e = 0.0, sin_t = 0.0,
           cos_t = 0.0;
    int moving   = 0;
    // Get the orbital parameters for the source
    grid->Wind->get_src_orbit(
        id, &moving, &eccentricity_fac, &periastron_vec[0], &periastron_vec[1],
        &orbital_period, initial_position);
    // execute only if orbit exists == orbital_period!=0
    if (orbital_period > 1.0e-8) {
#ifndef NDEBUG
      spdlog::debug(
          "   >>> orbital period <<< {}, simtime={}", orbital_period, simtime);
#endif
      grid->Wind->get_src_posn(id, srcpos);
      srcrad = SWP.params[id]->radius;
#ifndef NDEBUG
      cout << "current pos, size = [" << srcpos[0] << "," << srcpos[1] << "]  "
           << srcrad << "\n";
#endif
      // delete wind-cell list
      err += grid->Wind->remove_cells(id);
      // Update positions
      // Elipse needs to be rotated -> get rotation matrix entries
      if (periastron_vec[0] != 0) {
        cos_a = -1.0 * periastron_vec[0] / abs(periastron_vec[0])
                * cos(atan(periastron_vec[1] / periastron_vec[0]));
      }
      else {
        cos_a = -1.0;  // both omitted terms are zero is vec[0]==0
      }
      if (periastron_vec[1] != 0) {
        sin_a =
            sin(-1.0 * periastron_vec[1] / abs(periastron_vec[1])
                * acos(cos_a));  // sqrt(1-cos_a*cos_a);
      }
      else {
        sin_a = sin(-1.0 * acos(cos_a));
      }
#ifndef NDEBUG
      cout << "(cos_a,sin_a) = " << cos_a << ", " << sin_a;
      cout << " ;  peri-vec= (" << periastron_vec[0] << ", "
           << periastron_vec[1] << ")\n";
#endif
      // Get elipse parameters
      a = sqrt(
              periastron_vec[0] * periastron_vec[0]
              + periastron_vec[1] * periastron_vec[1])
          / (1.0 - eccentricity_fac);
      b = a * sqrt(1 - eccentricity_fac * eccentricity_fac);
#ifndef NDEBUG
      cout << "(a,b)=(" << a << "," << b << "), e-fac= " << eccentricity_fac
           << "\n";
#endif
      sin_t = sin(2.0 * pconst.pi() * simtime / (orbital_period));
      cos_t = cos(2.0 * pconst.pi() * simtime / (orbital_period));
      // Set new position from orbital parameters
      newpos[0] = initial_position[0] - a * cos_a + cos_a * a * cos_t
                  - sin_a * b * sin_t;
      newpos[1] = initial_position[1] - a * sin_a + sin_a * a * cos_t
                  + cos_a * b * sin_t;
      // Orbit has to be in the x-y-plane -> z-component unchanged
      if (MAX_DIM > 2) newpos[2] = initial_position[2];
#ifndef NDEBUG
      rep.printVec("original pos", initial_position, 3);
      rep.printVec("updated  pos", newpos, 3);
#endif
      // Set new source position
      grid->Wind->set_src_posn(id, newpos);
      // Asign new Boundary cells
      grid->Wind->get_src_posn(id, srcpos);
      srcrad = SWP.params[id]->radius;
      BC_assign_STWIND_add_cells2src(par, grid, id);
#ifndef NDEBUG
      cout << "new pos, size = [" << srcpos[0] << "," << srcpos[1] << "]  "
           << srcrad << "\n";
#endif
      // End of source loop
    }
  }
#else  // NOT ANALYTIC_ORBITS

  // loop over wind-sources, and calculate the orbital evolution using
  // leapfrog integrator, drift-kick-drift method.
  // https://en.wikipedia.org/wiki/Leapfrog_integration
  // First get list of stars that are moving according to gravity
  vector<struct starpos> stars;
  int ndim = grid->Ndim();
  for (int i = 0; i < SWP.Nsources; i++) {
    if (SWP.params[i]->moving_star > 0) {
      struct starpos s;
      s.id   = i;
      s.mass = SWP.params[i]->Mass;
      for (int v = 0; v < ndim; v++)
        s.pos[v] = SWP.params[i]->dpos[v];
      for (int v = 0; v < ndim; v++)
        s.vel[v] = SWP.params[i]->velocity[v];
      for (int v = 0; v < ndim; v++)
        s.acc[v] = 0.0;
      stars.push_back(s);
    }
  }

  // only move source if we are on the finest level:
  if (par.grid_nlevels == l + 1) {
    // with masses, positions and velocities, calculate acceleration,
    // brute-force n^2 method
    double sep[MAX_DIM], dist, acc;
    for (unsigned long i = 0; i < stars.size(); i++) {
      // cout<<"init "<<i<<", "; rep.printVec("star pos",stars[i].pos,3);
      // cout<<"     "<<i<<", "; rep.printVec("star vel",stars[i].vel,3);
      // first get acceleration
      for (unsigned long j = 0; j < i; j++) {
        for (int v = 0; v < ndim; v++)
          sep[v] = stars[j].pos[v] - stars[i].pos[v];
        dist = 0.0;
        for (int v = 0; v < ndim; v++)
          dist += sep[v] * sep[v];
        dist = sqrt(dist);
        acc  = pconst.G() / (dist * dist * dist);
        for (int v = 0; v < ndim; v++) {
          stars[i].acc[v] += sep[v] * stars[j].mass * acc;
          stars[j].acc[v] -= sep[v] * stars[i].mass * acc;
        }
      }
    }

    // step forward in time
    for (unsigned long i = 0; i < stars.size(); i++) {
      for (int v = 0; v < ndim; v++) {
        stars[i].vel[v] += stars[i].acc[v] * 0.5 * dt;
        stars[i].pos[v] += stars[i].vel[v] * dt;
      }
      // if (i==1) {
      //  cout <<i<<" dt= "<<dt<<"  "; rep.printVec("vel",stars[i].vel,2);
      //  cout <<i<<" pos    "; rep.printVec("pos",stars[i].pos,2);
      //}
    }

    // need acceleration at t+dt to get v at t+dt:
    for (unsigned long i = 0; i < stars.size(); i++) {
      for (int v = 0; v < ndim; v++)
        stars[i].acc[v] = 0.0;
      // first get acceleration
      for (unsigned long j = 0; j < i; j++) {
        for (int v = 0; v < ndim; v++)
          sep[v] = stars[j].pos[v] - stars[i].pos[v];
        dist = 0.0;
        for (int v = 0; v < ndim; v++)
          dist += sep[v] * sep[v];
        dist = sqrt(dist);
        acc  = pconst.G() / (dist * dist * dist);
        for (int v = 0; v < ndim; v++) {
          stars[i].acc[v] += sep[v] * stars[j].mass * acc;
          stars[j].acc[v] -= sep[v] * stars[i].mass * acc;
        }
      }
    }
    for (unsigned long i = 0; i < stars.size(); i++) {
      for (int v = 0; v < ndim; v++) {
        stars[i].vel[v] += stars[i].acc[v] * 0.5 * dt;
      }
    }

    // update star parameters
    for (unsigned long i = 0; i < stars.size(); i++) {
      for (int v = 0; v < ndim; v++)
        SWP.params[stars[i].id]->dpos[v] = stars[i].pos[v];
      for (int v = 0; v < ndim; v++)
        SWP.params[stars[i].id]->velocity[v] = stars[i].vel[v];
    }
  }

  //
  // Decide if we need to update star position based on whether its new
  // position is more than 0.1*dx away from the old position
  //
  double dx_fine = par.levels[par.grid_nlevels - 1].dx;
  static std::array<double, MAX_DIM> my_pos;
  bool moved = false;

  for (unsigned long i = 0; i < stars.size(); i++) {
    bool move = false;
    // if star has moved > 0.1dx
    grid->Wind->get_this_source_position(stars[i].id, my_pos);

    if (par.grid_nlevels == l + 1) {
      if ((SWP.params[i]->moving_star > 0)
          && (grid->distance(SWP.params[i]->dpos, my_pos) > 0.1 * dx_fine)) {
        move = true;
        if (!moved) moved = true;
      }
    }
    else {
      if (!pconst.equalD(0.0, grid->distance(SWP.params[i]->thispos, my_pos))) {
        move = true;
      }
    }

    // if star has moved > 0.1dx
    if (move || par.timestep == 0) {
      if (par.grid_nlevels == l + 1) {
        for (int v = 0; v < MAX_DIM; v++)
          SWP.params[i]->thispos[v] = SWP.params[i]->dpos[v];
      }

      // delete wind-cell list (we have to move the star, so everything
      // gets rewritten)
      err += grid->Wind->remove_cells(stars[i].id);

      grid->Wind->set_this_source_position(stars[i].id, SWP.params[i]->thispos);
      BC_set_wind_radius(par, grid, stars[i].id);
      spdlog::debug(
          "l {}, source moved >0.1 * {}, updating pos from {} to {}", l,
          dx_fine, my_pos, SWP.params[i]->thispos);

      // add cells back to wind source, based on new position
      BC_assign_STWIND_add_cells2src(par, grid, stars[i].id);
    }  // if move
  }    // loop over stars to move positions

  static long unsigned int count = 0;

  // write trajectory to file.
  if (moved && (par.grid_nlevels == l + 1)
#ifdef PARALLEL
      && par.levels[0].sub_domain.get_myrank() == 0
#endif
  ) {
    spdlog::debug("stars have moved, ct={}", count);
    if ((count % 10) == 0) {
      spdlog::debug("writing to trajectory file, ct={}", count);
      // if (!outf.is_open())
      //  spdlog::error("output file is not open!");
      outf << simtime << "  " << dt << "  ";
      for (unsigned long i = 0; i < stars.size(); i++) {
        for (int v = 0; v < ndim; v++)
          outf << SWP.params[i]->dpos[v] << "  ";
        for (int v = 0; v < ndim; v++)
          outf << SWP.params[i]->velocity[v] << "  ";
      }
      outf << "\n";
      // outf.flush();
    }
    count++;
  }  // print to file

#endif  // ANALYTIC_ORBITS

  //
  // The stellar_wind class already has a list of cells to update
  // for each source, together with pre-calculated state vectors,
  // so we just call the set_cell_values() function.
  // The stellar evolution (if present) is updated in
  // Wind->set_cell_values()
  //
#ifndef NDEBUG
  spdlog::info("stellar_wind_bc: updating wind boundary");
#endif
  // int err=0;
  for (int id = 0; id < grid->Wind->Nsources(); id++) {
#ifndef NDEBUG
    spdlog::debug(
        "stellar_wind_bc: updating wind boundary for id={} of {}", id,
        grid->Wind->Nsources());
#endif
    err += grid->Wind->set_cell_values(grid, id, simtime);
    if (err == -1) {
      spdlog::info("Star reached end of life: set_cell_values {}", err);
      par.maxtime = true;
      err         = 0;
    }
    else if (err > 0) {
      spdlog::error(
          "error in setting cell values for wind src {}: {}", id, err);
      exit(err);
    }
    // set wind acceleration in each cell
    BC_set_windacc_radflux(par, grid, id);
  }
#ifndef NDEBUG
  spdlog::debug("stellar_wind_bc: finished updating wind boundaries");
#endif

  return err;
}



// ##################################################################
// ##################################################################



void stellar_wind_bc::BC_set_wind_radius(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int id                ///< source id
)
{
  int lev = par.grid_nlevels - 1;
  //
  // Determine the radius we need for the wind region based on star location,
  // imposing a minimum of MIN_WIND_RAD on any grid that the boundary
  // intersects with.
  double rad = SWP.params[id]->radius;
  if (par.ndim == 1) {
    if (lev == 0)
      rad = max(rad, 2.0 * par.dx);
    else
      rad = max(rad, 2.0 * par.levels[lev].dx);
    SWP.params[id]->current_radius = rad;
    spdlog::info(
        "1D sim: resetting wind radius {:9.3e} to {:9.3e}",
        SWP.params[id]->radius, SWP.params[id]->current_radius);
    return;
  }
  // spdlog::info("BC_set_wind_radius 1: {:9.3e}  {:9.3e}", rad,
  // SWP.params[id]->current_radius);

  // impose that wind radius must be at least 2 cells larger than Rstar
  // if on finest level
  rad = max(rad, grid->Wind->get_min_wind_radius(id, par.levels[lev].dx));
  // spdlog::info("BC_set_wind_radius 2: {:9.3e}  {:9.3e}  {:9.3e}", rad,
  // SWP.params[id]->current_radius,par.levels[lev].dx);

  bool fin = false, on = true;
  do {
    // check that radius satisfies minimum radius criterion:
    rad = max(rad, MIN_WIND_RAD * par.levels[lev].dx);
    // spdlog::info("BC_set_wind_radius 3: {:9.3e}  {:9.3e}", rad,
    //              SWP.params[id]->current_radius);

    // see if any part of wind boundary is outside level. If not: break out,
    // if so: continue to next coarser level.
    for (int v = 0; v < par.ndim; v++) {
      if (par.ndim == 2 && v == 1 && par.coord_sys == COORD_CYL) {
        continue;
      }
      if ((SWP.params[id]->thispos[v] - rad < par.levels[lev].Xmin[v])
          || (SWP.params[id]->thispos[v] + rad > par.levels[lev].Xmax[v]))
        on = false;
    }
    if (on || lev == 0) {
      rad = max(rad, grid->Wind->get_min_wind_radius(id, par.levels[lev].dx));
      // spdlog::info("BC_set_wind_radius 4: {:9.3e}  {:9.3e}", rad,
      //            SWP.params[id]->current_radius);
      SWP.params[id]->current_radius = rad;
      fin                            = true;
    }
    else {
      lev--;
      on = true;
    }
  } while (lev >= 0 && !fin);

  // Now wind radius should be at least MIN_WIND_RAD cells radius on the
  // coarsest level containing the source, and at least 2 cells larger than
  // the stellar radius.
  // spdlog::info(
  //    "BC_set_wind_radius: orig: {:9.3e} now: {:9.3e}",
  //    SWP.params[id]->radius, SWP.params[id]->current_radius);
  return;
}



// ##################################################################
// ##################################################################



void stellar_wind_bc::BC_set_windacc_radflux(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to grid.
    const int id                ///< source id
)
{
  //
  // For each cell, store radiation energy density if compton cooling is needed
  // For cells at <100 stellar radii, calculate the wind acceleration,
  // if requested
  //
  if (!par.EP.compton_cool && !par.EP.wind_acceleration) {
    return;
  }

  int ncell    = 0;
  double r_acc = 100.0 * SWP.params[id]->Rstar;
  // u_prefactor is a prefactor to get the radiation energy density
  double u_prefactor = 2.0 * pconst.StefanBoltzmannConst()
                       * pow(SWP.params[id]->Tstar, 4) / pconst.c();
  array<double, MAX_DIM> srcpos;
  grid->Wind->get_this_source_position(id, srcpos);

#ifndef NDEBUG
  spdlog::debug(
      "star id {}, lev {}, radius {:8.2e}, dx {:8.2e}", id, grid->level(),
      r_acc, grid->DX());
#endif

  // find min/max of coordinates of cube that contains spherical wind source
  double dx = grid->DX(), aneg = 0.0, apos = 0.0, xmin = 0.0, xmax = 0.0;
  array<int, MAX_DIM> ineg = {0, 0, 0}, ipos = {0, 0, 0};
  for (int v = 0; v < par.ndim; v++) {
    xmin    = grid->Xmin_all(static_cast<axes>(v));
    xmax    = grid->Xmax_all(static_cast<axes>(v));
    aneg    = srcpos[v] - r_acc - 0.5 * dx;
    aneg    = max(xmin, aneg);
    aneg    = min(xmax, aneg);
    apos    = srcpos[v] + r_acc + 0.5 * dx;
    apos    = min(xmax, apos);
    apos    = max(xmin, apos);
    ineg[v] = (aneg - xmin) / dx;
    ipos[v] = (apos - xmin) / dx;
  }
  for (int v = par.ndim; v < MAX_DIM; v++) {
    ineg[v] = 0;
    ipos[v] = 1;
  }

#ifndef NDEBUG
  spdlog::info("*** r_acc={}, ineg {}, ipos {}", r_acc, ineg, ipos);
#endif

  // loop over cube of data and add only cells with position within the
  // source radius
  enum axes x1 = XX;
  enum axes x2 = YY;
  enum axes x3 = ZZ;
  double vinf  = SWP.params[id]->Vinf;
  // set launching speed to be 4x the sound speed at the star's surface.
  double v0    = 4.0 * sqrt(pconst.kB() * SWP.params[id]->Tstar / pconst.m_p());
  double rstar = SWP.params[id]->Rstar;


  // we shut off acceleration if T>10^6K

  if (!par.EP.compton_cool) {
    //#ifdef PION_OMP
    //  #pragma omp parallel for collapse(2)
    //#endif
    for (int ax3 = ineg[x3]; ax3 < ipos[x3]; ax3++) {
      for (int ax2 = ineg[x2]; ax2 < ipos[x2]; ax2++) {
        std::array<double, MAX_DIM> cpos;
        int index[3];
        index[x1] = ineg[x1];
        index[x2] = static_cast<int>(ax2);
        index[x3] = static_cast<int>(ax3);
        double v = 0.0, dist = 0.0, T = 0.0;
        std::vector<double> acc;
        acc.resize(par.ndim);
        cell *c = grid->get_cell_all(index[0], index[1], index[2]);
        for (int ax1 = index[x1]; ax1 < ipos[x1]; ax1++) {
          CI.get_dpos(*c, cpos);
          if ((dist = grid->distance(srcpos, cpos)) <= r_acc) {
            // approximate T by assuming mu=1
            T = c->Ph[PG] * pconst.m_p() / (pconst.kB() * c->Ph[RO]);
            if (T > 2.0e6) {
              // For high temperatures switch off acc b/c no ions left.
              for (int d = 0; d < par.ndim; d++)
                acc[d] = 0.0;
            }
            else {
              // v = v0 + (vinf-v0)* exp(beta * log(1.0-rstar/dist));
              // v = v * beta * (vinf-v0) *  exp((beta-1.0) *
              // log(1.0-rstar/dist))
              // * rstar / pow(dist,3);
              // first assume beta==1 for simplicity
              v = v0 + (vinf - v0) * (1.0 - rstar / dist);
              v = v * (vinf - v0) * rstar / (dist * dist * dist);
              // calculate each component.
              for (int d = 0; d < par.ndim; d++) {
                acc[d] = v * (cpos[d] - srcpos[d]);
              }
              if (T > 1.0e6) {
                // linearly decrease acc to zero in range 1e6-2e6 K
                for (int d = 0; d < par.ndim; d++)
                  acc[d] *= (2.0e6 - T) / 1.0e6;
              }
            }
            CI.set_wind_acceleration(*c, id, acc);
            // spdlog::info("src id {}, acceleration {}",id,r_acc);
          }
          else {
            for (int d = 0; d < par.ndim; d++) {
              acc[d] = 0.0;
            }
            CI.set_wind_acceleration(*c, id, acc);
          }
          c = grid->NextPt(*c, XP);
          // spdlog::info("src id {} index {}",index[x1]);
        }
      }
    }
  }
  else {
    // do compton cooling in every cell, and wind acceleration too.
    enum axes x1 = XX;
    enum axes x2 = YY;
    enum axes x3 = ZZ;
    int nx2      = grid->NG_All(x2);
    int nx3      = grid->NG_All(x3);
#ifdef PION_OMP
    #pragma omp parallel for collapse(2)
#endif
    for (int ax3 = 0; ax3 < nx3; ax3++) {
      for (int ax2 = 0; ax2 < nx2; ax2++) {
        int index[3];
        index[x1] = 0;
        index[x2] = static_cast<int>(ax2);
        index[x3] = static_cast<int>(ax3);
        double v = 0.0, dist = 0.0, T = 0.0;
        std::vector<double> acc;
        acc.resize(par.ndim);
        std::array<double, MAX_DIM> cpos;
        cell *c = grid->get_cell_all(index[0], index[1], index[2]);
        do {
          CI.get_dpos(*c, cpos);
          dist = grid->distance(srcpos, cpos);
          // set wind acceleration
          if (par.EP.wind_acceleration) {
            if (dist <= r_acc) {
              // approximate T by assuming mu=1
              T = c->Ph[PG] * pconst.m_p() / (pconst.kB() * c->Ph[RO]);
              if (T > 2.0e6) {
                // For high temperatures switch off acc b/c no ions left.
                for (int d = 0; d < par.ndim; d++)
                  acc[d] = 0.0;
              }
              else {
                // v = v0 + (vinf-v0)* exp(beta * log(1.0-rstar/dist));
                // v = v * beta * (vinf-v0) *  exp((beta-1.0) *
                // log(1.0-rstar/dist))
                // * rstar / pow(dist,3);
                // first assume beta==1 for simplicity
                v = v0 + (vinf - v0) * (1.0 - rstar / dist);
                v = v * (vinf - v0) * rstar / (dist * dist * dist);
                // calculate each component.
                for (int d = 0; d < par.ndim; d++) {
                  acc[d] = v * (cpos[d] - srcpos[d]);
                }
                if (T > 1.0e6) {
                  // linearly decrease acc to zero in range 1e6-2e6 K
                  for (int d = 0; d < par.ndim; d++)
                    acc[d] *= (2.0e6 - T) / 1.0e6;
                }
              }
              CI.set_wind_acceleration(*c, id, acc);
              // spdlog::info("src id {}, acceleration {}",id,r_acc);
            }
            else {
              for (int d = 0; d < par.ndim; d++) {
                acc[d] = 0.0;
              }
              CI.set_wind_acceleration(*c, id, acc);
            }
          }
          // set compton cooling info (radiation density from star, erg/cm3)
          // U = 2 sigma_SB T^4 /c * (1 - sqrt(1-(R/r)^2))
          // which includes the dilution factor accounting for finite size of
          // the star
          // only do this for leaf cells because coarse-grid cells can have
          // very outlandish temperatures because of averaging.
          if (c->isleaf && c->isdomain) {
            v = SWP.params[id]->Rstar / dist;
            v *= v;
            v = u_prefactor * (1.0 - sqrt(1.0 - v));
          }
          else {
            v = 0.0;
          }
          CI.set_compton_urad(*c, id, v);
        } while ((c = grid->NextPt(*c, XP)) != 0);
      }  // axis 2
    }    // axis 3
  }      // if compton cooling

#ifndef NDEBUG
  spdlog::debug(
      "BC_assign_STWIND_add_cells2src: Added {} cells to wind boundary for WS {}",
      ncell, id);
#endif
  return;
}


// ##################################################################
// ##################################################################
