///
/// \file NG_MPI_BC89flux.cpp
///
/// \author Jonathan Mackey
///
/// Function definitions for fluxes leaving/entering
/// levels on a nested grid, to make sure that fluxes between
/// adjacent levels are consistent with each other.  This file
/// includes extra code needed for parallel execution.
///
/// \author Jonathan Mackey
///
/// Modifications:n
/// - 2019.12.04 JM: moved functions from uniform_grid and
///   sim_control_NG classes to gather everything in one file.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_MPI_BC89flux.h"
#include "constants.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#include <iostream>
using namespace std;

// ##################################################################
// ##################################################################

int NG_MPI_BC89flux::setup_flux_recv(
    class SimParams &par,       ///< simulation params (including BCs)
    class GridBaseClass *grid,  // pointer to coarse grid.
    const int lp1               ///< fine level to receive from
)
{
#ifdef TEST_BC89FLUX
  cout << "NG_MPI_BC89flux::setup_flux_recv() recv from level=";
  cout << lp1 << "\n";
#endif
  int l     = lp1 - 1;  // my level
  size_t nc = 1;        // number of cells in each interface
  int ixmin[MAX_DIM], ixmax[MAX_DIM], ncell[MAX_DIM];  // interface region
  int f_lxmin[MAX_DIM], f_lxmax[MAX_DIM];              // finer grid
  int d_xmin[MAX_DIM], d_xmax[MAX_DIM];                // full domain
  struct flux_interface *fi = 0;
  class MCMDcontrol *MCMD   = &(par.levels[l].MCMD);
  vector<struct cgrid> cg;
  MCMD->get_child_grid_info(cg);
  int nchild = cg.size();

  // coarse grid quantities
  int cg_ixmin[MAX_DIM], cg_ixmax[MAX_DIM];
  int idx = grid->idx();

  // if only one MPI process, then no send/recv necessary and we
  // call the serial version:
  if (par.levels[0].MCMD.get_nproc() == 1) {
#ifdef TEST_BC89FLUX
    cout << "setup_flux_recv(): nproc=1, calling serial version.\n";
#endif
    int err = NG_BC89flux::setup_flux_recv(par, grid, lp1);
    rep.errorTest("UniformGrid flux setup", 0, err);
    err = flux_update_recv[l].size();
    for (int i = 0; i < err; i++) {
      flux_update_recv[l][i].rank.push_back(0);
      flux_update_recv[l][i].dir = i;
    }
    return 0;
  }

  //
  // Domain sizes on different levels, and local grid properties
  //
  CI.get_ipos_vec(par.levels[lp1].Xmin, f_lxmin);
  CI.get_ipos_vec(par.levels[lp1].Xmax, f_lxmax);
  CI.get_ipos_vec(par.levels[0].Xmin, d_xmin);
  CI.get_ipos_vec(par.levels[0].Xmax, d_xmax);
  for (int v = 0; v < par.ndim; v++)
    cg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    cg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

#ifdef TEST_BC89FLUX
  cout << "NG_MPI_BC89flux::setup_flux_recv: " << nchild << " child grids\n";
#endif

  // Two cases: if there are children, then part or all of grid is
  // within l+1 level.  If there are no children, then at most one
  // face of my grid could be an outer face of the l+1 grid.
  if (nchild > 0) {
    bool recv[nchild * 2 * par.ndim];   // whether to get data in this dir
    size_t nel[nchild * 2 * par.ndim];  // number of interfaces in each dir

    for (int ic = 0; ic < nchild; ic++) {
      int off = ic * 2 * par.ndim;
      CI.get_ipos_vec(cg[ic].Xmin, ixmin);
      CI.get_ipos_vec(cg[ic].Xmax, ixmax);

      // define interface region of fine and coarse grids, and if
      // each direction is to be included or not.
      for (int ax = 0; ax < par.ndim; ax++) {
        if (f_lxmin[ax] == ixmin[ax] && f_lxmin[ax] != d_xmin[ax]
            && ixmin[ax] == cg_ixmin[ax])
          recv[off + 2 * ax] = false;
        else if (f_lxmin[ax] == ixmin[ax] && f_lxmin[ax] != d_xmin[ax])
          recv[off + 2 * ax] = true;
        else
          recv[off + 2 * ax] = false;

        if (f_lxmax[ax] == ixmax[ax] && f_lxmax[ax] != d_xmax[ax]
            && ixmax[ax] == cg_ixmax[ax])
          recv[off + 2 * ax + 1] = false;
        else if (f_lxmax[ax] == ixmax[ax] && f_lxmax[ax] != d_xmax[ax])
          recv[off + 2 * ax + 1] = true;
        else
          recv[off + 2 * ax + 1] = false;

        ncell[ax] = (ixmax[ax] - ixmin[ax]) / idx;
        if ((ixmax[ax] - ixmin[ax]) % idx != 0) {
          rep.error("interface region not divisible!", ixmax[ax] - ixmin[ax]);
        }
      }  // all dimensions
      for (int d = 0; d < 2 * par.ndim; d++)
        nel[off + d] = 0;
    }  // all child grids.

    for (int ic = 0; ic < nchild; ic++) {
      int off = ic * 2 * par.ndim;
      // different number of interfaces depending on dimensionality.
      switch (par.ndim) {
        case 1:
          if (recv[off + XN]) nel[off + XN] = 1;
          if (recv[off + XP]) nel[off + XP] = 1;
          break;
        case 2:
          if (recv[off + XN]) nel[off + XN] = ncell[YY];
          if (recv[off + XP]) nel[off + XP] = ncell[YY];
          if (recv[off + YN]) nel[off + YN] = ncell[XX];
          if (recv[off + YP]) nel[off + YP] = ncell[XX];
          break;
        case 3:
          if (recv[off + XN]) nel[off + XN] = ncell[YY] * ncell[ZZ];
          if (recv[off + XP]) nel[off + XP] = ncell[YY] * ncell[ZZ];
          if (recv[off + YN]) nel[off + YN] = ncell[XX] * ncell[ZZ];
          if (recv[off + YP]) nel[off + YP] = ncell[XX] * ncell[ZZ];
          if (recv[off + ZN]) nel[off + ZN] = ncell[XX] * ncell[YY];
          if (recv[off + ZP]) nel[off + ZP] = ncell[XX] * ncell[YY];
          break;
        default:
          rep.error("bad ndim in setup_flux_recv", par.ndim);
          break;
      }  // dims
    }    // all child grids

    // initialize arrays
    flux_update_recv[l].resize(nchild * 2 * par.ndim);
    for (int ic = 0; ic < nchild; ic++) {
#ifdef TEST_BC89FLUX
      cout << "NG_MPI_BC89flux::setup_flux_recv: ";
      cout << ic << " has " << 2 * par.ndim << "boundaries,\n";
#endif
      int off = ic * 2 * par.ndim;
      for (int d = 0; d < 2 * par.ndim; d++) {
        int el = off + d;
#ifdef TEST_BC89FLUX
        cout << "NG_MPI_BC89flux::setup_flux_recv: ";
        cout << ic << " d=" << d << ", recv[el]=" << recv[el] << ", nel=";
        cout << nel[el] << "\n";
#endif
        if (recv[el] == true) {
          flux_update_recv[l][el].fi.resize(nel[el]);
          flux_update_recv[l][el].dir    = d;
          flux_update_recv[l][el].ax     = d / 2;
          flux_update_recv[l][el].Ncells = nc;
          flux_update_recv[l][el].rank.push_back(cg[ic].rank);
          for (size_t i = 0; i < nel[el]; i++) {
            flux_update_recv[l][el].fi[i] =
                mem.myalloc(flux_update_recv[l][el].fi[i], 1);
            fi = flux_update_recv[l][el].fi[i];
            fi->c.resize(nc);
            fi->area.resize(nc);
            fi->flux = mem.myalloc(fi->flux, par.nvar);
            for (int v = 0; v < par.nvar; v++)
              fi->flux[v] = 0.0;
          }
        }
        else {
          flux_update_recv[l][el].fi.resize(1);
          flux_update_recv[l][el].fi[0]  = 0;
          flux_update_recv[l][el].dir    = d;
          flux_update_recv[l][el].ax     = d / 2;
          flux_update_recv[l][el].Ncells = nc;
        }
      }  // loop over dims
      CI.get_ipos_vec(cg[ic].Xmin, ixmin);
      CI.get_ipos_vec(cg[ic].Xmax, ixmax);

      // For each interface, find the cell that is off the fine grid
      // and that includes the interface.
      for (int d = 0; d < 2 * par.ndim; d++) {
        if (recv[off + d]) {
#ifdef TEST_BC89FLUX
          cout << "parallel setup_flux_recv(): adding cells for child grids.\n";
#endif
          add_cells_to_face(
              par, grid, static_cast<enum direction>(d), ixmin, ixmax, ncell, 1,
              flux_update_recv[l][off + d]);
        }
      }  // loop over dims.
    }    // loop over child grids.
  }      // if there are child grids.

  else {
#ifdef TEST_BC89FLUX
    cout << "parallel setup_flux_recv(): no child grids, looking for "
            "neighbours.\n";
#endif
    // There are no children, but my grid might have a boundary in
    // common with the l+1 level's outer boundary.  There is at most
    // one of these.  Use MCMD functions to test this.
    int edge = -1, axis = -1, perp[2];
    size_t nsub = 0;
    size_t nel  = 0;
    std::vector<std::vector<struct cgrid> > cgngb;

    MCMD->get_level_lp1_ngb_info(cgngb);
    if (cgngb.size() != static_cast<size_t>(2 * par.ndim)) {
      rep.error("l+1 neigbouring grids vector not set up right", cgngb.size());
    }

    // try to find a direction (there is at most one).  If we find
    // one then populate flux_update_recv[d] with values.
    // default is that there are no neighbours:
    flux_update_recv[l].resize(1);
    flux_update_recv[l][0].fi.resize(1);
    flux_update_recv[l][0].fi[0] = 0;
    for (int d = 0; d < 2 * par.ndim; d++) {
      nsub = cgngb[d].size();
#ifdef TEST_BC89FLUX
      cout << "pllel setup_flux_recv(): dir=" << d << ", nsub=" << nsub << "\n";
#endif
      if (nsub > 0) {
        flux_update_recv[l].resize(nsub);
        edge    = d;
        axis    = d / 2;
        perp[0] = (axis + 1) % par.ndim;
        perp[1] = (axis + 2) % par.ndim;

        for (size_t ic = 0; ic < nsub; ic++) {
          for (int i = 0; i < MAX_DIM; i++)
            ncell[i] = 0;

          flux_update_recv[l][ic].rank.push_back(cgngb[d][ic].rank);
          flux_update_recv[l][ic].dir =
              grid->OppDir(static_cast<enum direction>(edge));
          flux_update_recv[l][ic].ax = axis;

          // figure out how many cells in the interface region:
          // "nel" is the number of cells in the interface region on
          // the coarse grid.
          CI.get_ipos_vec(cgngb[d][ic].Xmin, ixmin);  // l+1 grid
          CI.get_ipos_vec(cgngb[d][ic].Xmax, ixmax);  // l+1 grid

          switch (par.ndim) {

            case 1:
              nel         = 1;
              ncell[axis] = 1;
              break;

            case 2:
              ixmin[perp[0]] = std::max(ixmin[perp[0]], cg_ixmin[perp[0]]);
              ixmax[perp[0]] = std::min(ixmax[perp[0]], cg_ixmax[perp[0]]);
              ncell[axis]    = 1;
              ncell[perp[0]] = (ixmax[perp[0]] - ixmin[perp[0]]) / idx;
              nel            = ncell[perp[0]];
#ifdef TEST_BC89FLUX
              cout << "BC89 recv 2D from l+1 ngb: xmax=" << ixmax[perp[0]];
              cout << ", xmin=" << ixmin[perp[0]] << ", nel=" << nel << "\n";
#endif
              break;

            case 3:
              ixmin[perp[0]] = std::max(ixmin[perp[0]], cg_ixmin[perp[0]]);
              ixmax[perp[0]] = std::min(ixmax[perp[0]], cg_ixmax[perp[0]]);
              ixmin[perp[1]] = std::max(ixmin[perp[1]], cg_ixmin[perp[1]]);
              ixmax[perp[1]] = std::min(ixmax[perp[1]], cg_ixmax[perp[1]]);
              ncell[axis]    = 1;
              ncell[perp[0]] = (ixmax[perp[0]] - ixmin[perp[0]]) / idx;
              ncell[perp[1]] = (ixmax[perp[1]] - ixmin[perp[1]]) / idx;
              nel            = 1;
              for (int p = 0; p < 2; p++)
                nel *= ncell[perp[p]];
#ifdef TEST_BC89FLUX
              cout << "BC89 recv 3D from l+1 ngb: nel=" << nel << "\n";
              rep.printVec("interface xmin", ixmin, 3);
              rep.printVec("interface xmax", ixmax, 3);
#endif
              break;
          }

          flux_update_recv[l][ic].fi.resize(nel);
          for (size_t i = 0; i < nel; i++) {
            flux_update_recv[l][ic].fi[i] =
                mem.myalloc(flux_update_recv[l][ic].fi[i], 1);
            fi = flux_update_recv[l][ic].fi[i];
            fi->c.resize(nc);
            fi->area.resize(nc);
            fi->flux = mem.myalloc(fi->flux, par.nvar);
            for (int v = 0; v < par.nvar; v++)
              fi->flux[v] = 0.0;
          }
#ifdef TEST_BC89FLUX
          cout << "ic=" << ic << "parallel setup_flux_recv(): ";
          cout << "adding cells for l+1 neigbouring grids.\n";
#endif
          add_cells_to_face(
              par, grid, grid->OppDir(static_cast<enum direction>(edge)), ixmin,
              ixmax, ncell, 1, flux_update_recv[l][ic]);

        }  // loop over sub-domains in direction d
      }    // if sub-domains exist
    }      // loop over directions
  }        // child grids?

#ifdef TEST_BC89FLUX
  for (int l = 0; l < par.grid_nlevels; l++) {
    for (unsigned int d = 0; d < flux_update_send[l].size(); d++) {
      struct flux_update *fup = &(flux_update_send[l][d]);
      cout << "l=" << l << ", FUP_SEND: d=" << d << " info: \n";
      cout << "\t ranks: ";
      for (unsigned int r = 0; r < fup->rank.size(); r++)
        cout << fup->rank[r] << "  ";
      cout << "\n";
      cout << "\t nel = " << fup->fi.size() << " dir=" << fup->dir
           << ", ax=" << fup->ax << "\n";
    }
    for (unsigned int d = 0; d < flux_update_recv[l].size(); d++) {
      struct flux_update *fup = &(flux_update_recv[l][d]);
      cout << "l=" << l << ", FUP_RECV: d=" << d << " info: \n";
      cout << "\t ranks: ";
      for (unsigned int r = 0; r < fup->rank.size(); r++)
        cout << fup->rank[r] << "  ";
      cout << "\n";
      cout << "\t nel = " << fup->fi.size() << " dir=" << fup->dir
           << ", ax=" << fup->ax << "\n";
    }
  }
#endif

  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_BC89flux::setup_flux_send(
    class SimParams &par,       ///< simulation params (including BCs)
    class GridBaseClass *grid,  // pointer to finer grid.
    const int lm1               ///< level to send to
)
{
#ifdef TEST_BC89FLUX
  cout << " NG_MPI_BC89flux::setup_flux_send() send to level=";
  cout << lm1 << " from MY LEVEL l=" << lm1 + 1 << "\n";
#endif

  int err = NG_BC89flux::setup_flux_send(par, grid, lm1);
  rep.errorTest("NG_BC89flux::setup_flux_send", 0, err);

  // Add ranks for each send, based on parent rank.
  int l = lm1 + 1;  // my level
  int fg_ixmin[MAX_DIM], fg_ixmax[MAX_DIM];
  for (int v = 0; v < par.ndim; v++)
    fg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    fg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

  // Size of flux_update_send[] is ALWAYS 2*ndim
  int ns = flux_update_send[l].size();
  if (ns != 2 * par.ndim) rep.error("bad flux send size", ns);

  if (par.levels[l].MCMD.get_nproc() == 1) {
    for (int d = 0; d < ns; d++) {
      flux_update_send[l][d].rank.push_back(0);
      flux_update_send[l][d].dir = d;
    }
  }
  else {
    // Usually send to parent, sometimes send to neighbour of parent
    // if boundary between parent and neighbour sits on the edge of
    // my level.
    class MCMDcontrol *MCMD = &(par.levels[l].MCMD);

    // get data on parent grid, and neighbours of parent grid.
    struct cgrid pg;
    std::vector<struct cgrid> pgngb;
    MCMD->get_parent_grid_info(&pg);
    MCMD->get_parent_ngb_grid_info(pgngb);
    int pg_ixmin[MAX_DIM], pg_ixmax[MAX_DIM];
    CI.get_ipos_vec(pg.Xmin, pg_ixmin);
    CI.get_ipos_vec(pg.Xmax, pg_ixmax);

    for (int ax = 0; ax < par.ndim; ax++) {
#ifdef TEST_BC89FLUX
      cout << "setup_flux_send: axis=" << ax << endl;
#endif
      int d = -1, pp = -1;
      // negative direction
      d = 2 * ax;
      if (flux_update_send[l][d].fi[0] != 0) {
        // if cells outside interface are not on parent, then send to
        // neighbour of parent.
        if (fg_ixmin[ax] == pg_ixmin[ax]) {
          pp = pgngb[d].rank;
#ifdef TEST_BC89FLUX
          cout << "setup_flux_send: negative direction, ngb of parent " << pp
               << endl;
#endif
        }
        // else it is parent.
        else {
          pp = pg.rank;
#ifdef TEST_BC89FLUX
          cout << "negative direction, send to parent " << pp << endl;
#endif
        }
        flux_update_send[l][d].dir = d;
        flux_update_send[l][d].ax  = ax;
        flux_update_send[l][d].rank.push_back(pp);
#ifdef TEST_BC89FLUX
        cout << "setup_flux_send: axis=" << ax << ", d=" << d
             << ", send rank = " << pp << endl;
#endif
      }
      // positive direction
      d = 2 * ax + 1;
      if (flux_update_send[l][d].fi[0] != 0) {
        // if cells outside interface are not on parent, then send to
        // neighbour of parent.
        if (fg_ixmax[ax] == pg_ixmax[ax]) {
          pp = pgngb[d].rank;
#ifdef TEST_BC89FLUX
          cout << "setup_flux_send: positive direction, ngb of parent " << pp
               << endl;
#endif
        }
        // else it is parent.
        else {
          pp = pg.rank;
#ifdef TEST_BC89FLUX
          cout << "setup_flux_send: positive direction, send to parent " << pp
               << endl;
#endif
        }
        flux_update_send[l][d].dir = d;
        flux_update_send[l][d].ax  = ax;
        flux_update_send[l][d].rank.push_back(pp);
#ifdef TEST_BC89FLUX
        cout << "axis=" << ax << ", d=" << d << ", send rank = " << pp << endl;
#endif
      }
    }  // loop over axes
  }    // if nproc>1

  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_BC89flux::clear_sends_BC89_fluxes()
{
  for (unsigned int ib = 0; ib < BC89_flux_send_list.size(); ib++) {
#ifdef TEST_BC89FLUX
    cout << "clear_sends_BC89_fluxes: waiting for send " << ib;
    cout << " of " << BC89_flux_send_list.size() << endl;
#endif
    COMM->wait_for_send_to_finish(BC89_flux_send_list[ib]);
  }
  BC89_flux_send_list.clear();
  return;
}

// ##################################################################
// ##################################################################

int NG_MPI_BC89flux::send_BC89_fluxes_F2C(
    class SimParams &par,  ///< simulation params (including BCs)
    const int l,           ///< My level in grid hierarchy.
    const int step,        ///< OA1 or OA2
    const int ooa          ///< Full order of accuracy of simulation
)
{
#ifdef SKIP_BC89_FLUX
  return 0;
#endif
  if (step != ooa) {
    cout << "l=" << l << ": don't send fluxes on half step\n";
    return 1;
  }
  if (par.levels[l].step % 2 != 0) {
    rep.error("Don't call BC89 Flux-SEND on odd steps", 1);
  }
  if (l == 0) {
    rep.error("Coarsest level trying to send flux data", l);
  }
  int err = 0;

  // else we have to send data to at least one other MPI process:
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  int n_send              = flux_update_send[l].size();
  for (int isend = 0; isend < n_send; isend++) {
    // loop over boundaries (some send to more than one process)
    struct flux_update *fup = &(flux_update_send[l][isend]);
    // some sends may be null, so we skip them:
    if (fup->fi[0] == 0) {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": BC89_FLUX send " << isend << " of " << n_send;
      cout << " is null, continuing..." << endl;
#endif
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": BC89_FLUX send " << isend;
      cout << " is not null: sending data." << endl;
#endif
    }
    size_t n_el    = fup->fi.size();
    size_t n_data  = n_el * par.nvar;
    pion_flt *data = 0;  // buffer for sending by MPI
    data           = mem.myalloc(data, n_data);
    size_t iel     = 0;
    for (size_t ii = 0; ii < n_el; ii++) {
      for (int v = 0; v < par.nvar; v++) {
        data[iel] = fup->fi[ii]->flux[v];
        iel++;
      }
    }
    //
    // Send data using a non-blocking MPI send
    //
    string id;

    // loop over procs on level l-1 to send to.
    for (unsigned int ii = 0; ii < fup->rank.size(); ii++) {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": isend=" << isend << ", send " << ii << " of ";
      cout << fup->rank.size() << " to rank " << fup->rank[ii] << endl;
#endif
      if (fup->rank[ii] == MCMD->get_myrank()) {
#ifdef TEST_BC89FLUX
        cout << "l=" << l << ": ignoring BC89 for isend=" << isend
             << " (to myrank)." << endl;
#endif
        continue;
      }
      // unique as long as isend<30, l<100.
      int comm_tag = BC_MPI_FLUX_tag + 100 * isend + l;
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": BC89 FLUX: Sending " << n_data;
      cout << " doubles from proc " << MCMD->get_myrank();
      cout << " to parent proc " << fup->rank[ii] << endl;
#endif
      err += COMM->send_double_data(fup->rank[ii], n_data, data, id, comm_tag);
      if (err) rep.error("FLUX_F2C send_data failed.", err);
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": BC89 FLUX: returned with id=" << id;
      cout << endl;
#endif
      // store ID to clear the send later
      BC89_flux_send_list.push_back(id);
    }  // loop over ranks
    data = mem.myfree(data);
  }  // loop over send boundaries.
#ifdef TEST_BC89FLUX
  cout << "BC89 MPI flux send finished" << endl;
#endif
  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_BC89flux::recv_BC89_fluxes_F2C(
    class FV_solver_base *spatial_solver,  ///< solver, for gradients
    class SimParams &par,  ///< simulation params (including BCs)
    const int l,           ///< My level in grid hierarchy.
    const int step,        ///< OA1 or OA2
    const int ooa          ///< Full order of accuracy of simulation
)
{
#ifdef SKIP_BC89_FLUX
  return 0;
#endif
  if (step != ooa) {
    cout << "don't receive fluxes on half step\n";
    rep.error("trying to receive BC89 flux on half step", l);
  }
  if (l == par.grid_nlevels - 1) {
    rep.error("finest level trying to receive data from l+1", l);
  }
  int err = 0;

  // loop over all boundaries that we need to receive
  double ftmp[par.nvar], utmp[par.nvar];
  for (int v = 0; v < par.nvar; v++)
    ftmp[v] = 0.0;
  class GridBaseClass *grid = par.levels[l].grid;
  class MCMDcontrol *MCMD   = &(par.levels[l].MCMD);
  int n_bd                  = flux_update_recv[l].size();
  double dt                 = par.levels[l].dt;

  for (int irecv = 0; irecv < n_bd; irecv++) {
    struct flux_update *fup = &(flux_update_recv[l][irecv]);
    // some recvs may be null, so we skip them:
    if (fup->fi[0] == 0) {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": recv " << irecv << " is null, continuing...\n";
#endif
      continue;
    }
    else if (fup->rank[0] == MCMD->get_myrank()) {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": calling serial BC89 for irecv=" << irecv;
      cout << " (same rank). dir=" << fup->dir << endl;
#endif
      int d                     = fup->dir;
      struct flux_update *fsend = &(flux_update_send[l + 1][d]);
      err                       = recv_BC89_flux_boundary(
          spatial_solver, par, grid, dt, *fsend, *fup, d,
          static_cast<axes>(fup->ax));
      rep.errorTest("serial BC89 call", 0, err);
      continue;
    }
    else {
#ifdef TEST_BC89FLUX
      cout << "l=" << l << ": recv " << irecv << " is not null: ";
      cout << "recving data from " << fup->rank[0] << endl;
#endif
    }

    // Normally we have nchild*2*ndim boundaries, so the direction
    // associated with a given "irecv" is irecv%(2*ndim).
    // But if only the grid boundary is in common with level l+1,
    // then there are only 2^(ndim-1) elements, one for each l+1
    // grid that is adjacent to the level l grid.
    int dir = -1;
    // if (n_bd >= 2*par.ndim) dir = irecv % 2*par.ndim;
    // else                      dir = irecv;
    dir = fup->dir;

    struct flux_interface *fi = 0;
    size_t n_el               = fup->fi.size();
    size_t n_data             = n_el * par.nvar;
    pion_flt *buf             = 0;
    buf                       = mem.myalloc(buf, n_data);

    // receive data: comm_tag ensures that we select the data
    // relating to this value of "irecv".
    string recv_id;
    int recv_tag  = -1;
    int from_rank = fup->rank[0];
    // unique as long as isend<30, l<10, rank<10000.
    int comm_tag = BC_MPI_FLUX_tag + 100 * dir + (l + 1);
#ifdef TEST_BC89FLUX
    cout << "looking for data with tag: " << comm_tag << endl;
#endif
    err = COMM->look_for_data_to_receive(
        &from_rank,  // rank of sender (output)
        recv_id,     // identifier for receive (output).
        &recv_tag,   // comm_tag associated with data (output)
        comm_tag,
        COMM_DOUBLEDATA  // type of data we want.
    );
    if (err) rep.error("FLUX look for double data failed", err);
#ifdef TEST_BC89FLUX
    cout << "l=" << l << ": BC89 FLUX: found data from rank ";
    cout << from_rank << ", and ID " << recv_id << endl;
#endif

    //
    // Receive data into buffer.  Data stored for each fine cell:
    // flux[nv]
    //
    err = COMM->receive_double_data(from_rank, recv_tag, recv_id, n_data, buf);
    if (err) rep.error("(flux BC89) getdata failed", err);

    size_t iel = 0;
    for (size_t ii = 0; ii < n_el; ii++) {
      fi = fup->fi[ii];
#ifdef TEST_BC89FLUX
      cout << "f=" << ii << ":coarse=" << fi << ", flux =  ";
      rep.printVec("fc->flux", fi->flux, par.nvar);
      cout << "f=" << ii << ":  fine=" << fi << ", flux =  ";
      rep.printVec("ff->flux", &(buf[iel]), par.nvar);
#endif

      // construct error in flux:
      for (int v = 0; v < par.nvar; v++) {
        // make flux[v] be the difference of coarse and fine.
        fi->flux[v] += buf[iel + v];
#ifdef TEST_BC89FLUX
        if (!isfinite(buf[iel])) {
          cout << "l=" << l << ": element " << ii << " of FLUX C2F RECV: ";
          cout << " var " << iel << " is " << buf[iel] << endl;
          rep.error("infinite", buf[ii]);
        }
#endif
      }
      iel += par.nvar;

      // divide by face area so that it is a flux.
      for (int v = 0; v < par.nvar; v++)
        fi->flux[v] /= fi->area[0];
      for (int v = 0; v < par.nvar; v++)
        ftmp[v] = 0.0;
      // re-calculate dU based on error in flux.
      if (fup->dir % 2 == 0) {
        spatial_solver->DivStateVectorComponent(
            fi->c[0], grid, static_cast<axes>(fup->ax), par.nvar, ftmp,
            fi->flux, utmp);
      }
      else {
        spatial_solver->DivStateVectorComponent(
            fi->c[0], grid, static_cast<axes>(fup->ax), par.nvar, fi->flux,
            ftmp, utmp);
      }
#ifdef TEST_BC89FLUX
      rep.printVec("**********  Error", utmp, par.nvar);
      int q = 1;
      cout << "Flux Energy: " << fi->flux[q] << ": " << fi->c[0]->dU[q];
      cout << ", " << utmp[q] << "\n";
#endif
      // correct dU so that coarse level is consistent with fine.
      for (int v = 0; v < par.nvar; v++)
        fi->c[0]->dU[v] += utmp[v];

    }  // loop over elements
    if (iel != n_data) rep.error("ndata", iel - n_data);
    buf = mem.myfree(buf);
#ifdef TEST_BC89FLUX
    cout << "l=" << l << ": BC89 FLUX: finished with recv ID ";
    cout << recv_id << endl;
#endif

  }  // loop over faces in this direction.

#ifdef TEST_BC89FLUX
  cout << "l=" << l << ": BC89 FLUX: finished with recv" << endl;
#endif

  return 0;
}

// ##################################################################
// ##################################################################
