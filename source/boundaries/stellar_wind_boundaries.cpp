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
  outf.open("trajectory.txt");
  outf.setf(ios_base::scientific);
  outf.precision(6);
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
    spdlog::error("{}: {}", "STWIND not external boundary!", b->dir);
    exit(1);
  }
#ifndef NDEBUG
  spdlog::debug("Assigning data to STWIND boundary. Nsrc={}", SWP.Nsources);
#endif
  if (SWP.Nsources < 1) {
    spdlog::error("{}: {}", "BC_assign_STWIND() No Sources!", SWP.Nsources);
    exit(1);
  }

  //
  // Setup reference state vector and initialise to zero.
  //
  if (b->refval) {
    spdlog::error(
        "{}: {}", "Initialised STWIND boundary refval", fmt::ptr(b->refval));
    exit(1);
  }
  b->refval = mem.myalloc(b->refval, par.nvar);
  if (!b->data.empty()) {
    spdlog::error(
        "{}: {}", "BC_assign_STWIND: Not empty boundary data", b->itype);
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
            "{}: Expected {} but got {}", "wind xi values don't match", xi,
            SWP.params[isw]->xi);
        exit(1);
      }
    }
  }
  // initialise "current_radius" var to value in parameter file.
  for (int isw = 0; isw < Ns; isw++) {
    SWP.params[isw]->current_radius = SWP.params[isw]->radius;
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

  //
  // Run through sources and add sources.
  //
  for (int isw = 0; isw < Ns; isw++) {
#ifndef NDEBUG
    spdlog::debug("\tBC_assign_STWIND: Adding source {}", isw);
#endif
    if (SWP.params[isw]->type == WINDTYPE_CONSTANT) {
      //
      // This is for spherically symmetric winds that are constant
      // in time.
      //
      err = grid->Wind->add_source(SWP.params[isw]);
#ifndef NDEBUG
      spdlog::debug(
          "add_source call: {},{}", SWP.params[isw]->PeriastronX,
          SWP.params[isw]->PeriastronY);
#endif
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
      spdlog::error("{}: {}", "Error adding wind source", isw);
      exit(1);
    }
    // if star is moving, then set initial velocity to values at periastron
    // if (SWP.params[isw]->moving_star) {
    //  double pre = 2.0 * pconst.pi() * sqrt(1.0 -
    //}
  }

  //
  // loop over sources, adding cells to boundary data list in order.
  //
  // Perform only if t=0

  for (int id = 0; id < Ns; id++) {
#ifndef NDEBUG
    spdlog::debug("\tBC_assign_STWIND: Adding cells to source {}", id);
#endif
    BC_assign_STWIND_add_cells2src(par, grid, id);
  }

  //
  // Now we should have set everything up, so we assign the boundary
  // cells with their boundary values.
  //
  err += BC_update_STWIND(par, 0, grid, par.simtime, 0.0, b, 0, 0);
#ifndef NDEBUG
  spdlog::info(
      "\tFinished setting up wind parameters\n------ DONE SETTING UP STELLAR WIND CLASS ----------\n");
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
  grid->Wind->get_src_posn(id, srcpos);


  if (grid->Wind->get_num_cells(id) != 0) {
    spdlog::error(
        "{}: {}", "adding cells to source that already has cells",
        grid->Wind->get_num_cells(id));
    exit(1);
  }

  // find min/max of coordinates of cube that contains spherical wind source
  double dx = grid->DX(), aneg = 0.0, apos = 0.0, xmin = 0.0, xmax = 0.0;
  array<int, MAX_DIM> ineg = {0, 0, 0}, ipos = {0, 0, 0};
  for (int v = 0; v < par.ndim; v++) {
    xmin    = grid->Xmin_all(static_cast<axes>(v));
    xmax    = grid->Xmax_all(static_cast<axes>(v));
    aneg    = srcpos[v] - srcrad - 0.5 * dx;  // add extra safety factor dx/2
    aneg    = max(xmin, aneg);
    aneg    = min(xmax, aneg);
    apos    = srcpos[v] + srcrad + 0.5 * dx;  // add extra safety factor dx/2
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
          // spdlog::info("src id {}, add cell {}",id,index[x1]);
        }
        c = grid->NextPt(*c, XP);
        // spdlog::info("src id {} index {}",index[x1]);
      }
    }
  }

  err += grid->Wind->set_num_cells(id, ncell);

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
    boundary_data *b,           ///< Boundary to update.
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
      err += grid->Wind->remove_cells(grid, id);
      err += grid->Wind->set_num_cells(id, 0);
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
#else   // NOT ANALYTIC_ORBITS

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
      // delete wind-cell list (we have to movet the star, so everything
      // gets rewritten)
      err += grid->Wind->remove_cells(grid, i);
      err += grid->Wind->set_num_cells(i, 0);
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
    // Set new source position
    if (stars.size() >= 1) {
      outf << simtime << "  " << dt << "  ";
      for (unsigned long i = 0; i < stars.size(); i++) {
        for (int v = 0; v < ndim; v++)
          SWP.params[stars[i].id]->dpos[v] = stars[i].pos[v];
        for (int v = 0; v < ndim; v++)
          SWP.params[stars[i].id]->velocity[v] = stars[i].vel[v];
        for (int v = 0; v < ndim; v++)
          outf << SWP.params[i]->dpos[v] << "  ";
        for (int v = 0; v < ndim; v++)
          outf << SWP.params[i]->velocity[v] << "  ";
      }
      outf << "\n";
    }
    // Determine the radius we need for the wind region based on star location,
    // imposing a minimum of MIN_WIND_RAD on any grid that the boundary
    // intersects with.
    static int count = 0;
    for (unsigned long i = 0; i < stars.size(); i++) {
      int lev    = l;
      double rad = SWP.params[stars[i].id]->radius;
      bool fin = false, on = true;
      do {
        // check that radius satisfies minimum radius criterion:
        rad = max(rad, MIN_WIND_RAD * par.levels[lev].dx);

        // see if any part of wind boundary is outside level. If not: break out,
        // if so: continue to next coarser level.
        for (int v = 0; v < par.ndim; v++) {
          if ((stars[i].pos[v] - rad < par.levels[lev].Xmin[v])
              || (stars[i].pos[v] + rad > par.levels[lev].Xmax[v]))
            on = false;
        }
        if (on) {
          SWP.params[stars[i].id]->current_radius = rad;
          fin                                     = true;
        }
        else {
          lev--;
          on = true;
        }
      } while (lev >= 0 && !fin);
      if (count % 512 == 0)
        spdlog::info("star {} on level {}, radius = {:12.6e}", i, lev, rad);
    }
    count++;
  }

  for (int i = 0; i < SWP.Nsources; i++) {
    if (SWP.params[i]->moving_star > 0) {
      // Asign new Boundary cells
      BC_assign_STWIND_add_cells2src(par, grid, SWP.params[i]->id);
      // cout<<"end  "<<i<<", "; rep.printVec("star pos",stars[i].pos,3);
      // cout<<"     "<<i<<", "; rep.printVec("star vel",stars[i].vel,3);
    }
  }
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
    if (err) {
      spdlog::info("Star reached end of life? set_cell_values {}", err);
      par.maxtime = true;
    }
  }
#ifndef NDEBUG
  spdlog::debug("stellar_wind_bc: finished updating wind boundaries");
#endif

  return err;
}

// ##################################################################
// ##################################################################
