/// \file NG_MPI_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_MPI_coarse_to_fine boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.09.06 JM: started writing code.

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include "boundaries/NG_MPI_coarse_to_fine_boundaries.h"
#include "tools/mem_manage.h"
#include <sstream>
using namespace std;

//#define TEST_C2F
//#define TEST_MPI_NG

// ##################################################################
// ##################################################################

int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_SEND(
    class SimParams &par,  ///< simulation parameters
    const int l,           ///< level of this grid.
    boundary_data *b       ///< boundary data
)
{

  class GridBaseClass *grid = par.levels[l].grid;
  // see how many child grids I have
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
  vector<struct cgrid> fg;
  sub_domain->get_child_grid_info(fg);
  int nchild = fg.size();

#ifdef TEST_C2F
  if (nchild == 0) {
    spdlog::debug("COARSE_TO_FINE_SEND: no children");
  }
  else {
    spdlog::debug("COARSE_TO_FINE_SEND: {} child grids", nchild);
  }
#endif

  // Two cases:
  // (1) if level l+1 has a boundary within my domain, then
  // send data to child grids.
  // (2) if level l+1 has a boundary coincident with my boundary,
  // then I need to send data if none of the l+1 domain intersects
  // my grid.  Otherwise another process can do it.

  std::array<int, MAX_DIM> cl_ixmin, cl_ixmax;
  std::array<int, MAX_DIM> fl_ixmin, fl_ixmax;
  std::array<int, MAX_DIM> cg_ixmin, cg_ixmax;
  std::array<int, MAX_DIM> fg_ixmin, fg_ixmax;

  CI.get_ipos_vec(par.levels[l].Xmin, cl_ixmin);
  CI.get_ipos_vec(par.levels[l].Xmax, cl_ixmax);
  CI.get_ipos_vec(par.levels[l + 1].Xmin, fl_ixmin);
  CI.get_ipos_vec(par.levels[l + 1].Xmax, fl_ixmax);
  for (int v = 0; v < par.ndim; v++)
    cg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    cg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

  // loop over child grids
  for (int i = 0; i < nchild; i++) {

    if (sub_domain->get_myrank() == fg[i].rank) {
      // if child is on my process, do nothing because child grid
      // can grab the data directly using serial C2F code.
#ifdef TEST_C2F
      spdlog::debug(
          "C2F_SEND: child {}, my rank ({}) == child rank ({}), no need to set up COARSE_TO_FINE_SEND\n",
          i, sub_domain->get_myrank(), fg[i].rank);
#endif
    }
    else {
      // if child is on another process, make a list of cells that
      // we need to send.  One list for each external boundary that
      // we need to update.  Only external boundaries of the whole
      // domain are included.
#ifdef TEST_C2F
      spdlog::debug(
          "C2F_SEND: child {}, my rank != child rank ({}, {}) running parallel COARSE_TO_FINE_SEND\n",
          i, sub_domain->get_myrank(), fg[i].rank);
#endif
      // get dimensions of child grid from struct
      CI.get_ipos_vec(fg[i].Xmin, fg_ixmin);
      CI.get_ipos_vec(fg[i].Xmax, fg_ixmax);

      // loop over dimensions
      for (int d = 0; d < par.ndim; d++) {
        // if child xmin == its level xmin, but > my level xmin,
        // then we need to send data, so set up a list.
        if ((fg_ixmin[d] == fl_ixmin[d]) && (fg_ixmin[d] > cl_ixmin[d])
            && (fg_ixmin[d] != cg_ixmin[d])) {
#ifdef TEST_C2F
          spdlog::debug("C2F_SEND: child {}, dim {} NEG DIR", i, d);
          spdlog::debug("localxmin : {}", sub_domain->get_Xmin());
          spdlog::debug("Childxmin : {}", fg[i].Xmin);
#endif
          struct c2f *bdata = new struct c2f;
          bdata->rank       = fg[i].rank;
          bdata->dir        = 2 * d;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(
              par, grid, bdata, fg_ixmin, fg_ixmax, l + 1, fl_ixmin, fl_ixmax);
          b->NGsendC2F.push_back(bdata);
#ifdef TEST_C2F
          spdlog::debug(
              "added {} cells to C2F send el {}", bdata->c.size(),
              b->NGsendC2F.size() - 1);
#endif
        }
        // if child xmax == its level xmax, but < my level xmax,
        // then we need to send data, so set up a list.
        if ((fg_ixmax[d] == fl_ixmax[d]) && (fg_ixmax[d] < cl_ixmax[d])
            && (fg_ixmax[d] != cg_ixmax[d])) {
#ifdef TEST_C2F
          spdlog::debug("C2F_SEND: child {}, dim {} POS DIR", i, d);
#endif
          struct c2f *bdata = new struct c2f;
          bdata->rank       = fg[i].rank;
          bdata->dir        = 2 * d + 1;
          bdata->c.clear();

          // find cells along this boundary.
          add_cells_to_C2F_send_list(
              par, grid, bdata, fg_ixmin, fg_ixmax, l + 1, fl_ixmin, fl_ixmax);
          b->NGsendC2F.push_back(bdata);
#ifdef TEST_C2F
          spdlog::debug(
              "added {} cells to C2F send el {}", bdata->c.size(),
              b->NGsendC2F.size() - 1);
#endif
        }
      }  // loop over dimensions
    }    // if child is not on my process
  }      // loop over child grids

  // We've now dealt with C2F boundaries that are within my domain,
  // so we have to also consider boundaries coincident with my
  // domain boundary, where level l+1 only intersects at the boundary
  // this only applies if nproc>1.
  if (sub_domain->get_nproc() > 1) {
    std::vector<std::vector<struct cgrid> > fgngb;
    sub_domain->get_level_lp1_ngb_info(fgngb);
    if (fgngb.size() != static_cast<size_t>(2 * par.ndim)) {
      spdlog::error(
          "{}: {}", "l+1 neigbouring grids vector not set up right",
          fgngb.size());
    }
    for (int d = 0; d < par.ndim; d++) {
      // negative direction
      int dir = 2 * d;
      if (cg_ixmin[d] == fl_ixmax[d]) {
        if (fgngb[dir].size() != 0) {
          // there are neighbouring grids on l+1, so add all of them
          // to the send list.
          for (size_t f = 0; f < fgngb[dir].size(); f++) {
            if (fgngb[dir][f].rank == sub_domain->get_myrank()) {
              // if child is on my process, do nothing because child
              // grid can grab the data directly using serial C2F
              // code.
#ifdef TEST_C2F
              spdlog::debug(
                  "C2F_SEND: ngb {}, dir={}, my rank ({}) == child-ngb rank ({}), no need to set up COARSE_TO_FINE_SEND\n",
                  f, dir, sub_domain->get_myrank(), fgngb[dir][f].rank);
#endif
            }
            else {
              CI.get_ipos_vec(fgngb[dir][f].Xmin, fg_ixmin);
              CI.get_ipos_vec(fgngb[dir][f].Xmax, fg_ixmax);
              struct c2f *bdata = new struct c2f;
              bdata->rank       = fgngb[dir][f].rank;
              bdata->dir        = 2 * d + 1;  // outward normal of child is +ve
              bdata->c.clear();
              add_cells_to_C2F_send_list(
                  par, grid, bdata, fg_ixmin, fg_ixmax, l + 1, fl_ixmin,
                  fl_ixmax);
              b->NGsendC2F.push_back(bdata);
            }
          }
        }
      }
      // positive direction
      dir = 2 * d + 1;
      if (cg_ixmax[d] == fl_ixmin[d]) {
        if (fgngb[dir].size() != 0) {
          // there are neighbouring grids on l+1, so add all of them
          // to the send list.
          for (size_t f = 0; f < fgngb[dir].size(); f++) {
            if (fgngb[dir][f].rank == sub_domain->get_myrank()) {
              // if child is on my process, do nothing because child
              // grid can grab the data directly using serial C2F
              // code.
#ifdef TEST_C2F
              spdlog::debug(
                  "C2F_SEND: ngb {}, dir={}, my rank ({}) == child-ngb rank ({}), no need to set up COARSE_TO_FINE_SEND\n",
                  f, dir, sub_domain->get_myrank(), fgngb[dir][f].rank);
#endif
            }
            else {
              CI.get_ipos_vec(fgngb[dir][f].Xmin, fg_ixmin);
              CI.get_ipos_vec(fgngb[dir][f].Xmax, fg_ixmax);
              struct c2f *bdata = new struct c2f;
              bdata->rank       = fgngb[dir][f].rank;
              bdata->dir        = 2 * d;  // outward normal of child is -ve dir
              bdata->c.clear();
              add_cells_to_C2F_send_list(
                  par, grid, bdata, fg_ixmin, fg_ixmax, l + 1, fl_ixmin,
                  fl_ixmax);
              b->NGsendC2F.push_back(bdata);
            }
          }
        }
      }  // positive direction
    }    // loop over dimensions
  }      // if nproc>1

  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_SEND(
    class SimParams &par,          ///< simulation parameters
    class GridBaseClass *grid,     ///< pointer to coarse-level grid
    class FV_solver_base *solver,  ///< pointer to equations
    const int l,                   ///< level in the NG hierarchy
    struct boundary_data *b,       ///< pointer to boundary struct
    const int cstep,               ///< fractional step
    const int maxstep              ///< number of fractional steps
)
{
#ifdef TEST_C2F
  spdlog::debug(
      "MPI C2F SEND: starting from level {} to level {}... ", l, l + 1);
  if (b->NGsendC2F.size() == 0) {
    spdlog::warn("empty send list, so just returning now.");
    return 0;
  }
  else
    spdlog::debug(" send-list size = {}", b->NGsendC2F.size());
#endif
#ifdef TEST_C2F
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
#endif
  int err = 0;

  //
  // if on an odd-numbered step, then need to update the data on the
  // coarse grid to half way through a coarse step.  Assume dU has
  // been calculated for the coarse grid, but not updated to Ph[]
  //
  if (cstep != maxstep) {
#ifdef C2F_FULLSTEP
    return 0;
#endif
#ifdef TEST_C2F
    spdlog::debug("MPI C2F SEND: odd step, interpolating data in time");
#endif
    std::vector<double> U(par.nvar);
    for (unsigned int ib = 0; ib < b->NGsendC2F.size(); ib++) {
      for (unsigned int c_iter = 0; c_iter < b->NGsendC2F[ib]->c.size();
           c_iter++) {
        cell *c = b->NGsendC2F[ib]->c[c_iter];
        solver->PtoU(c->P.data(), U.data(), par.gamma);
        for (int v = 0; v < par.nvar; v++)
          U[v] += 0.5 * c->dU[v];
        solver->UtoP(U.data(), c->Ph, par.EP.MinTemperature, par.gamma);
#ifdef TEST_INF
        for (int v = 0; v < par.nvar; v++) {
          if (!isfinite(c->Ph[v])) {
            spdlog::debug("NAN c->P  : {}", c->P);
            spdlog::debug("NAN c->Ph : {}", c->Ph);
          }
        }
#endif
      }  // loop over cells in boundary
    }    // loop over send boundaries
  }      // if not at a full step update

#ifdef TEST_C2F
  spdlog::debug("C2F SEND: {} num send={}", l, b->NGsendC2F.size());
#endif

  // loop over send boundaries, pack and send the data.
  for (unsigned int ib = 0; ib < b->NGsendC2F.size(); ib++) {
    //
    // if 1st order accuracy, then just need Ph[]+cell-vol.
    // if 2nd order accuracy, also a slope vector for each dimension
    //
    size_t n_cell = b->NGsendC2F[ib]->c.size();

#ifdef TEST_C2F
    spdlog::debug(
        "C2F SEND: {}, sending {} elements to process: {}",
        sub_domain->get_myrank(), n_cell, b->NGsendC2F[ib]->rank);
#endif

    size_t n_el = 0;
    if (par.spOOA == OA1)
      n_el = n_cell * (par.nvar + 1 + par.ndim);
    else if (par.spOOA == OA2)
      n_el = n_cell * ((1 + par.ndim) * par.nvar + 1 + par.ndim);
    else
      spdlog::error("{}: {}", "bad spOOA in MPI C2F", par.spOOA);
    vector<pion_flt> buf(n_el);
    std::vector<double> slope(par.nvar);
    std::array<double, MAX_DIM> cpos;

    // loop over cells, add Ph[], cell-vol, slopes to send buffer
    size_t ibuf = 0;
    for (unsigned int c_iter = 0; c_iter < b->NGsendC2F[ib]->c.size();
         c_iter++) {
      cell *c = b->NGsendC2F[ib]->c[c_iter];
      for (int v = 0; v < par.nvar; v++)
        buf[ibuf + v] = c->Ph[v];
      ibuf += par.nvar;
      buf[ibuf] = grid->CellVolume(*c, 0);
      ibuf++;
      CI.get_dpos(*c, cpos);
      for (int v = 0; v < par.ndim; v++)
        buf[ibuf + v] = cpos[v];
      ibuf += par.ndim;

      if (par.spOOA == OA2) {
        for (int idim = 0; idim < par.ndim; idim++) {
          enum axes a = static_cast<axes>(idim);
          // cout <<"BC_update_COARSE_TO_FINE_SEND: el="<<c_iter;
          // cout <<" idim="<<idim<<" calling setslope on cell ";
          // cout <<c->id<<", isbd="<<c->isbd<<",
          // isgd="<<c->isgd<<"\n";
          solver->SetSlope(*c, a, par.nvar, &slope[0], OA2, grid);
          for (int v = 0; v < par.nvar; v++)
            buf[ibuf + v] = slope[v];
          ibuf += par.nvar;
        }  // idim
      }    // if 2nd order

    }  // loop over cells for this send boundary

    if (ibuf != n_el) spdlog::error("{}: {}", "C2F MPI SEND counting", ibuf);

    //
    // Send data using a non-blocking MPI send
    //
    // ostringstream tmp;
    // tmp <<"C2F_"<<sub_domain->get_myrank()<<"_to_"<<b->NGsendC2F[ib]->rank;
    string id;  // = tmp.str();
    // Need to add direction to comm-tag because we might be sending
    // more than one boundary to the same process.  Also add level,
    // because it can happen that more than one level sends the same
    // boundary to the same proc.
    // So the tag is BC_MPI_NGC2F_tag + 100*dir + level+1.
    //
    int comm_tag = BC_MPI_NGC2F_tag + 100 * b->NGsendC2F[ib]->dir + l + 1;
#ifdef TEST_C2F
    spdlog::debug(
        "BC_update_COARSE_TO_FINE_SEND: l={}, Sending {} doubles from proc {} to child proc {}, id={}",
        l, n_el, sub_domain->get_myrank(), b->NGsendC2F[ib]->rank, comm_tag);
#endif
    err += par.levels[l].sub_domain.send_double_data(
        b->NGsendC2F[ib]->rank, n_el, buf, id, comm_tag);
    if (err) spdlog::error("{}: {}", "Send_C2F send_data failed.", err);
#ifdef TEST_C2F
    spdlog::debug("BC_update_COARSE_TO_FINE_SEND: returned with id={}", id);
#endif
    // store ID to clear the send later (and delete the MPI temp data)
    par.levels[l].sub_domain.NG_C2F_send_list.push_back(id);
  }  // loop over send boundaries
  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::BC_COARSE_TO_FINE_SEND_clear_sends(
    class Sub_domain &sub_domain)
{
  if ((sub_domain.NG_C2F_send_list).empty()) {
    return;
  }
  for (unsigned int ib = 0; ib < sub_domain.NG_C2F_send_list.size(); ib++) {
#ifdef TEST_C2F
    spdlog::debug(
        "C2F_send: clearing send # {} of {}, id={}...", ib + 1,
        sub_domain.NG_C2F_send_list.size(), sub_domain.NG_C2F_send_list[ib]);
#endif
    sub_domain.wait_for_send_to_finish(sub_domain.NG_C2F_send_list[ib]);
#ifdef TEST_C2F
    spdlog::debug(" ... done!");
#endif
  }
  sub_domain.NG_C2F_send_list.clear();
  return;
}

// ##################################################################
// ##################################################################

int NG_MPI_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE_RECV(
    class SimParams &par,  ///< simulation parameters
    const int l,           ///< level of this grid.
    boundary_data *b       ///< boundary data
)
{
  // The boundary is already a regular external boundary, so just
  // need to decide if parent is on the same MPI process (and call
  // the serial version) or on a different MPI process (and set up
  // some data structures to receive the coarse-grid data for
  // interpolation).
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
  class GridBaseClass *grid    = par.levels[l].grid;

  std::array<int, MAX_DIM> cl_ixmin, cl_ixmax;
  std::array<int, MAX_DIM> fl_ixmin, fl_ixmax;
  std::array<int, MAX_DIM> cg_ixmin, cg_ixmax;
  std::array<int, MAX_DIM> fg_ixmin, fg_ixmax;

  CI.get_ipos_vec(par.levels[l - 1].Xmin, cl_ixmin);
  CI.get_ipos_vec(par.levels[l - 1].Xmax, cl_ixmax);
  CI.get_ipos_vec(par.levels[l].Xmin, fl_ixmin);
  CI.get_ipos_vec(par.levels[l].Xmax, fl_ixmax);
  for (int v = 0; v < par.ndim; v++)
    fg_ixmin[v] = grid->iXmin(static_cast<axes>(v));
  for (int v = 0; v < par.ndim; v++)
    fg_ixmax[v] = grid->iXmax(static_cast<axes>(v));

  // get data on parent grid, and neighbours of parent grid.
  struct cgrid pg;
  std::vector<struct cgrid> pgngb;
  sub_domain->get_parent_grid_info(&pg);
  sub_domain->get_parent_ngb_grid_info(pgngb);
  CI.get_ipos_vec(pg.Xmin, cg_ixmin);
  CI.get_ipos_vec(pg.Xmax, cg_ixmax);

  bool send2pg = true;
  switch (b->dir) {
    case XN:
      if (cg_ixmin[XX] == fg_ixmin[XX]) send2pg = false;
      break;
    case XP:
      if (cg_ixmax[XX] == fg_ixmax[XX]) send2pg = false;
      break;
    case YN:
      if (cg_ixmin[YY] == fg_ixmin[YY]) send2pg = false;
      break;
    case YP:
      if (cg_ixmax[YY] == fg_ixmax[YY]) send2pg = false;
      break;
    case ZN:
      if (cg_ixmin[ZZ] == fg_ixmin[ZZ]) send2pg = false;
      break;
    case ZP:
      if (cg_ixmax[ZZ] == fg_ixmax[ZZ]) send2pg = false;
      break;
    case NO:
      spdlog::error("{}: {}", "bad direction", b->dir);
      break;
    default:
      spdlog::error("{}: {}", "bad direction", b->dir);
      break;
  }
  if (!send2pg) {
    // must need to receive from neighbour of parent grid
    if (pgngb[b->dir].rank < 0)
      spdlog::error("{}: {}", "ngb of parent is null", b->dir);
    CI.get_ipos_vec(pgngb[b->dir].Xmin, cg_ixmin);
    CI.get_ipos_vec(pgngb[b->dir].Xmax, cg_ixmax);
    b->NGrecvC2F_parent = pgngb[b->dir].rank;
#ifdef TEST_C2F
    spdlog::debug(
        "C2F MPI Recv, recv from rank{}, in outward normal direction {}",
        b->NGrecvC2F_parent, b->dir);
#endif
  }
  else {
    b->NGrecvC2F_parent = pg.rank;
  }

  if (sub_domain->get_myrank() == b->NGrecvC2F_parent) {
#ifdef TEST_C2F
    spdlog::debug(
        "my rank == parent rank, setting up serial COARSE_TO_FINE_RECV");
#endif
    int err = NG_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
        par, grid, b, par.levels[l].parent);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "serial C2F BC setup", 0, err);
  }
  else {
#ifdef TEST_C2F
    spdlog::debug(
        "my rank != parent rank, so will need MPI for COARSE_TO_FINE_RECVme={}, parent={}",
        sub_domain->get_myrank(), b->NGrecvC2F_parent);
#endif

    //
    // group cells into lists that are child cells of a coarse cell
    //
    size_t n_cell = b->data.size();
    for (int idim = 0; idim < par.ndim; idim++)
      n_cell /= 2;
    b->NGrecvC2F.resize(n_cell);
    size_t ic                     = 0;
    list<cell *>::iterator f_iter = b->data.begin();

    if (par.ndim == 1) {
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        cell *c = (*f_iter);
        b->NGrecvC2F[ic].push_back(c);
        f_iter++;
        c = (*f_iter);
        b->NGrecvC2F[ic].push_back(c);
        ic++;
      }  // loop over cells
    }    // if 1D

    else if (par.ndim == 2) {
      // 4 fine cells per coarse cell
      f_iter  = b->data.begin();
      cell *c = (*f_iter), *temp = 0;
      int row_y = c->pos[YY];
      int idx   = grid->idx();
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        c = (*f_iter);
        if (c->pos[YY] - row_y == 2 * idx) {
          // move to next row of cells
          row_y = c->pos[YY];
        }
        if (c->pos[YY] == row_y) {
          // on same row of cells so continue adding cells
          b->NGrecvC2F[ic].push_back(c);
          f_iter++;
          temp = (*f_iter);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(*c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(*temp, YP));
          ic++;
        }
        else if (c->pos[YY] - row_y == idx) {
          // odd-numbered row, already dealt with.
          continue;
        }
        else
          spdlog::error(
              "{}: {}", "error in 2d logic C2FR_setup", c->pos[YY] - row_y);
      }  // loop over cells
    }    // if 2D

    else if (par.ndim == 3) {
      // 8 fine cells per coarse cell
      f_iter  = b->data.begin();
      cell *c = (*f_iter), *temp = 0;
      int row_y = c->pos[YY];
      int row_z = c->pos[ZZ];
      int idx   = grid->idx();
      for (f_iter = b->data.begin(); f_iter != b->data.end(); ++f_iter) {
        c = (*f_iter);
        if (c->pos[ZZ] - row_z == 2 * idx) {
          // move to next plane of cells
          row_z = c->pos[ZZ];
        }
        if (c->pos[YY] - row_y == 2 * idx) {
          // move to next row of cells in z-plane
          row_y = c->pos[YY];
        }
        if (row_y > c->pos[YY]) {
          row_y = c->pos[YY];
        }

        if (c->pos[YY] == row_y && c->pos[ZZ] == row_z) {
          // on same row of cells so continue adding cells
          f_iter++;
          temp = (*f_iter);
          b->NGrecvC2F[ic].push_back(c);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(*c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(*temp, YP));
          c    = grid->NextPt(*c, ZP);
          temp = grid->NextPt(*temp, ZP);
          b->NGrecvC2F[ic].push_back(c);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(*c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(*temp, YP));

          ic++;
        }
        else if (c->pos[YY] - row_y == idx || c->pos[ZZ] - row_z == idx) {
          // odd-numbered row/plane, already added cells.
          continue;
        }
        else
          spdlog::error(
              "{}: {}", "error in 3d logic C2FR_setup", c->pos[YY] - row_y);
      }  // loop over cells
    }    // if 3D
  }      // if parent is on other MPI process

  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_coarse_to_fine_bc::BC_update_COARSE_TO_FINE_RECV(
    class SimParams &par,          ///< simulation parameters
    class FV_solver_base *solver,  ///< pointer to equations
    const int l,                   ///< level in the NG grid structure
    struct boundary_data *b,       ///< boundary to update
    const int step                 ///< timestep on this (fine) grid
)
{
#ifdef C2F_FULLSTEP
  if ((step + 2) % 2 != 0) return 0;
#endif
#ifdef TEST_C2F
  spdlog::debug(
      "C2F_MPI: receiving boundary data to level {}, updating boundary dir = {}",
      l, b->dir);
#endif
  int err                      = 0;
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
  class GridBaseClass *grid    = par.levels[l].grid;

  if (sub_domain->get_myrank() == b->NGrecvC2F_parent) {
#ifdef TEST_C2F
    spdlog::debug("my rank == parent rank, calling serial COARSE_TO_FINE");
#endif
    NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(par, solver, l, b, step);
  }
  else {
#ifdef TEST_C2F
    spdlog::debug(
        "my rank != parent rank, so updating COARSE_TO_FINE_RECV: me={}, parent={}",
        sub_domain->get_myrank(), b->NGrecvC2F_parent);
#endif

    // receive data.
    string recv_id;
    int recv_tag  = -1;
    int from_rank = -1;
    int comm_tag  = BC_MPI_NGC2F_tag + 100 * b->dir + l;
#ifdef TEST_C2F
    spdlog::debug(
        "BC_update_COARSE_TO_FINE_RECV: looking for data tag {}", comm_tag);
#endif
    err = par.levels[l].sub_domain.look_for_data_to_receive(
        &from_rank, recv_id, &recv_tag, comm_tag, COMM_DOUBLEDATA);
    if (err) spdlog::error("{}: {}", "look for double data failed", err);
#ifdef TEST_C2F
    spdlog::debug(
        "BC_update_COARSE_TO_FINE_RECV: found data from rank {}, with tag {} and id {}.  Looked for comm_tag={}",
        from_rank, recv_tag, recv_id, comm_tag);
#endif

    // receive the data: nel is the number of coarse grid cells
    size_t n_cell = b->data.size();
    for (int idim = 0; idim < par.ndim; idim++)
      n_cell /= 2;
    size_t n_el = 0;
    if (par.spOOA == OA1)
      n_el = n_cell * (par.nvar + 1 + par.ndim);
    else if (par.spOOA == OA2)
      n_el = n_cell * ((1 + par.ndim) * par.nvar + 1 + par.ndim);
    else
      spdlog::error("{}: {}", "bad spOOA in MPI C2F", par.spOOA);
    vector<pion_flt> buf(n_el);
#ifdef TEST_C2F
    if (1 == 1) {
      spdlog::debug(
          "BC_update_COARSE_TO_FINE_RECV: get {} cells, and {} doubles", n_cell,
          n_el);
    }
#endif
    //
    // Receive data into buffer.  Data stored for each coarse cell:
    // Ph[nv],cellvol,cellpos[nd],slopeX[nv],slopeY[nv],slopeZ[nv]
    //
    err = par.levels[l].sub_domain.receive_double_data(
        from_rank, recv_tag, recv_id, n_el, buf);
    if (err)
      spdlog::error("{}: {}", "(BC_update_C2F_RECV) getdata failed", err);

    //
    // Get Ph and slopes from coarse cells, and interpolate into
    // groups of fine cells that are children of each coarse cell.
    //
    size_t ibuf = 0;
    if (par.spOOA == OA1) {
      // 1st order, so no averaging.  Fine cells have exactly the
      // same state as the coarse one.
      std::vector<double> Ph(par.nvar);
      std::vector<double> cpos(par.ndim);
#ifdef TEST_C2F
      std::vector<double> off(par.ndim);
#endif
      // double c_vol=0.0;
      cell *c = 0;
      for (unsigned int ic = 0; ic < b->NGrecvC2F.size(); ic++) {
        // read data for this coarse cell into arrays
        for (int v = 0; v < par.nvar; v++)
          Ph[v] = buf[ibuf + v];
        ibuf += par.nvar;
        // c_vol = buf[ibuf];
        ibuf++;
        for (int v = 0; v < par.ndim; v++)
          cpos[v] = buf[ibuf + v];
        ibuf += par.ndim;

        list<cell *>::iterator f_iter = b->NGrecvC2F[ic].begin();
        for (f_iter = b->NGrecvC2F[ic].begin();
             f_iter != b->NGrecvC2F[ic].end(); ++f_iter) {
          c = (*f_iter);
#ifdef TEST_C2F
          // for (int v=0;v<par.ndim;v++) off[v] = cpos[v]- c->pos[v];
          // cout <<"ic="<<ic<<", cell is "<<c<<"  ";
          // rep.printVec("offset is:",off,par.ndim);
#endif
          for (int v = 0; v < par.nvar; v++)
            c->Ph[v] = Ph[v];
          for (int v = 0; v < par.nvar; v++)
            c->P[v] = Ph[v];
        }  // loop over fine cells
      }    // loop over coarse cells
    }      // if 1st order accurate

    else if (par.spOOA == OA2) {
      // 2nd order needs a slope for each dimension, so split the
      // put different ndim in if/else statements.
      if (par.ndim == 1) {
        std::vector<double> Ph(par.nvar);
        std::vector<double> cpos(par.ndim);
        double c_vol = 0.0;
        std::vector<double> sx(par.nvar);
        cell *f[2];
        for (unsigned int ic = 0; ic < b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v = 0; v < par.nvar; v++)
            Ph[v] = buf[ibuf + v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v = 0; v < par.ndim; v++)
            cpos[v] = buf[ibuf + v];
          ibuf += par.ndim;
          for (int v = 0; v < par.nvar; v++)
            sx[v] = buf[ibuf + v];
          ibuf += par.nvar;
          list<cell *>::iterator f_iter = b->NGrecvC2F[ic].begin();
          for (int v = 0; v < 2; v++) {
            f[v] = *f_iter;
            f_iter++;
          }
          interpolate_coarse2fine1D(
              par, grid, solver, Ph.data(), c_vol, sx.data(), *f[0], *f[1]);
        }  // loop over coarse cells
      }    // if 1D
      else if (par.ndim == 2) {
        std::vector<double> Ph(par.nvar);
        std::array<double, MAX_DIM> cpos;
        double c_vol = 0.0;
        std::array<int, MAX_DIM> ipos;
        ipos[0] = 0;
        ipos[1] = 0;
        std::vector<double> sx(par.nvar), sy(par.nvar);
        cell *f[4];
        for (unsigned int ic = 0; ic < b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v = 0; v < par.nvar; v++)
            Ph[v] = buf[ibuf + v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v = 0; v < par.ndim; v++)
            cpos[v] = buf[ibuf + v];
          ibuf += par.ndim;
          for (int v = 0; v < par.nvar; v++)
            sx[v] = buf[ibuf + v];
          ibuf += par.nvar;
          for (int v = 0; v < par.nvar; v++)
            sy[v] = buf[ibuf + v];
          ibuf += par.nvar;
          list<cell *>::iterator f_iter = b->NGrecvC2F[ic].begin();
          for (int v = 0; v < 4; v++) {
            f[v] = *f_iter;
            f_iter++;
          }
          // rep.printVec("cpos",cpos,par.ndim);
          // rep.printVec("Ph",Ph,par.nvar);
          CI.get_ipos_vec(cpos, ipos);
          interpolate_coarse2fine2D(
              par, grid, solver, Ph.data(), ipos.data(), c_vol, sx.data(),
              sy.data(), *f[0], *f[1], *f[2], *f[3]);

        }  // loop over coarse cells
      }    // if 2D
      else {
        std::vector<double> Ph(par.nvar);
        std::array<double, MAX_DIM> cpos;
        double c_vol = 0.0;
#ifdef TEST_C2F
        std::array<double, MAX_DIM> cpos2;
#endif
        std::array<int, MAX_DIM> ipos;
        std::vector<double> sx(par.nvar), sy(par.nvar), sz(par.nvar);
        cell *fch[8];
        for (unsigned int ic = 0; ic < b->NGrecvC2F.size(); ic++) {
          // read data for this coarse cell into arrays
          for (int v = 0; v < par.nvar; v++)
            Ph[v] = buf[ibuf + v];
          ibuf += par.nvar;
          c_vol = buf[ibuf];
          ibuf++;
          for (int v = 0; v < par.ndim; v++)
            cpos[v] = buf[ibuf + v];
          ibuf += par.ndim;
          for (int v = 0; v < par.nvar; v++)
            sx[v] = buf[ibuf + v];
          ibuf += par.nvar;
          for (int v = 0; v < par.nvar; v++)
            sy[v] = buf[ibuf + v];
          ibuf += par.nvar;
          for (int v = 0; v < par.nvar; v++)
            sz[v] = buf[ibuf + v];
          ibuf += par.nvar;
#ifdef TEST_C2F
          // for (int v=0;v<par.ndim;v++) cpos2[v] = cpos[v]/3.086e18;
          // rep.printVec("*********** cpos",cpos2,par.ndim);
#endif
          list<cell *>::iterator f_iter = b->NGrecvC2F[ic].begin();
          for (int v = 0; v < 8; v++) {
            fch[v] = *f_iter;
#ifdef TEST_C2F
            // CI.get_dpos_vec(fch[v]->pos, cpos2);
            // for (int v=0;v<par.ndim;v++) cpos2[v] /= 3.086e18;
            // rep.printVec("fpos", cpos2, par.ndim);
#endif
            f_iter++;
          }
          CI.get_ipos_vec(cpos, ipos);
          interpolate_coarse2fine3D(
              par, grid, solver, Ph.data(), ipos.data(), c_vol, sx.data(),
              sy.data(), sz.data(), fch);
        }  // loop over coarse cells
      }    // if 3D
    }      // if 2nd order accurate
  }        // if parent proc is different to my proc.

#ifdef TEST_C2F
  spdlog::info("NG_MPI_C2F_bc::BC_update_COARSE_TO_FINE_RECV() done");
#endif
  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list(
    class SimParams &par,             ///< pointer to simulation parameters
    class GridBaseClass *grid,        ///< pointer to coarse-level grid
    struct c2f *bdata,                ///< pointer to list of cells
    std::array<int, MAX_DIM> &ixmin,  ///< child grid xmin (integer)
    std::array<int, MAX_DIM> &ixmax,  ///< child grid xmax (integer)
    const int lf,                     ///< level of fine grid.
    const std::array<int, MAX_DIM> &fl_xmin,  ///< level xmin of fine grid.
    const std::array<int, MAX_DIM> &fl_xmax   ///< level xmax of fine grid.
)
{
  // In XN,XP direction we add cells with faces that touch the fine-
  // level grid from the outside.
  // In YN,YP direction we also add corner data
  // In ZN,ZP direction we also add edge data
  //
  // easier to have a function for different grid dimensions.
#ifdef TEST_C2F
  spdlog::debug("Adding cells to C2F Send list: ");
#endif
  if (par.ndim == 1) {
#ifdef TEST_C2F
    spdlog::debug("1D");
#endif
    add_cells_to_C2F_send_list_1D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  else if (par.ndim == 2) {
#ifdef TEST_C2F
    spdlog::debug("2D");
#endif
    add_cells_to_C2F_send_list_2D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  else {
#ifdef TEST_C2F
    spdlog::debug("3D");
#endif
    add_cells_to_C2F_send_list_3D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_1D(
    class SimParams &par,             ///< pointer to simulation parameters
    class GridBaseClass *grid,        ///< pointer to coarse-level grid
    struct c2f *bdata,                ///< pointer to list of cells
    std::array<int, MAX_DIM> &ixmin,  ///< child grid xmin (integer)
    std::array<int, MAX_DIM> &ixmax,  ///< child grid xmax (integer)
    const int lf,                     ///< level of fine grid.
    const std::array<int, MAX_DIM> &fl_xmin,  ///< level xmin of fine grid.
    const std::array<int, MAX_DIM> &fl_xmax   ///< level xmax of fine grid.
)
{
  int bsize = grid->idx() * par.Nbc / 2;  // idx is >=2, Nbc is >=1.

  // define domain of boundary region
  int xn = 0, xp = 0;
  switch (bdata->dir) {
    case XN:
      xn = ixmin[XX] - bsize;
      xp = ixmin[XX];
      break;

    case XP:
      xn = ixmax[XX];
      xp = ixmax[XX] + bsize;
      break;

    default:
      spdlog::error("{}: {}", "bad direction in 1D C2F", bdata->dir);
  }

  cell *c = grid->FirstPt_All();
  do {
    if (c->pos[XX] > xn && c->pos[XX] < xp) bdata->c.push_back(c);
  } while ((c = grid->NextPt_All(*c)) != 0);

  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_2D(
    class SimParams &par,             ///< pointer to simulation parameters
    class GridBaseClass *grid,        ///< pointer to coarse-level grid
    struct c2f *bdata,                ///< pointer to list of cells
    std::array<int, MAX_DIM> &ixmin,  ///< child grid xmin (integer)
    std::array<int, MAX_DIM> &ixmax,  ///< child grid xmax (integer)
    const int lf,                     ///< level of fine grid.
    const std::array<int, MAX_DIM> &fl_xmin,  ///< level xmin of fine grid.
    const std::array<int, MAX_DIM> &fl_xmax   ///< level xmax of fine grid.
)
{
  // depth of boundary region, in integer coordinates.
  int bsize = grid->idx() * par.Nbc / 2;  // idx is >=2, Nbc is >=1.
  int n     = 0;
  // cout <<"bsize="<<bsize<<", and idx="<<grid->idx()<<", Nbc=";
  // cout <<par.Nbc<<"\n";
  // rep.printVec("ixmin",ixmin,par.ndim);
  // rep.printVec("ixmax",ixmax,par.ndim);

  // define domain of boundary region
  int xn = 0, xp = 0, yn = 0, yp = 0;
  switch (bdata->dir) {
    case XN:
      xn = ixmin[XX] - bsize;
      xp = ixmin[XX];
      yn = ixmin[YY];
      yp = ixmax[YY];
      break;

    case XP:
      xn = ixmax[XX];
      xp = ixmax[XX] + bsize;
      yn = ixmin[YY];
      yp = ixmax[YY];
      break;

    case YN:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      yn = ixmin[YY] - bsize;
      yp = ixmin[YY];
      break;

    case YP:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      yn = ixmax[YY];
      yp = ixmax[YY] + bsize;
      break;

    default:
      spdlog::error("{}: {}", "bad direction in 2D C2F", bdata->dir);
  }

  // cout <<"boundary: x in ["<<xn<<","<<xp<<"], y in["<<yn<<","<<yp<<"]\n";
  size_t ct = 0;
  cell *c   = grid->FirstPt_All();
  do {
    // rep.printVec("cpos",c->pos,par.ndim);
    if (c->pos[XX] > xn && c->pos[XX] < xp && c->pos[YY] > yn
        && c->pos[YY] < yp) {
      bdata->c.push_back(c);
      ct++;
    }
  } while ((c = grid->NextPt_All(*c)) != 0);
#ifdef TEST_C2F
  spdlog::debug("add_cells_to_C2F_send_list_2D: added {} cells", ct);
#endif

  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_3D(
    class SimParams &par,             ///< pointer to simulation parameters
    class GridBaseClass *grid,        ///< pointer to coarse-level grid
    struct c2f *bdata,                ///< pointer to list of cells
    std::array<int, MAX_DIM> &ixmin,  ///< child grid xmin (integer)
    std::array<int, MAX_DIM> &ixmax,  ///< child grid xmax (integer)
    const int lf,                     ///< level of fine grid.
    const std::array<int, MAX_DIM> &fl_xmin,  ///< level xmin of fine grid.
    const std::array<int, MAX_DIM> &fl_xmax   ///< level xmax of fine grid.
)
{
#ifdef TEST_C2F
  spdlog::debug("C2F Setup Send: ixmin : {}", ixmin);
  spdlog::debug("C2F Setup Send: ixmax : {}", ixmax);
#endif

  int bsize = grid->idx() * par.Nbc / 2;  // idx is >=2, Nbc is >=1.
  int n     = 0;

  // define domain of boundary region
  int xn = 0, xp = 0, yn = 0, yp = 0, zn = 0, zp = 0;
  switch (bdata->dir) {
    case XN:
      xn = ixmin[XX] - bsize;
      xp = ixmin[XX];
      yn = ixmin[YY];
      yp = ixmax[YY];
      zn = ixmin[ZZ];
      zp = ixmax[ZZ];
      break;

    case XP:
      xn = ixmax[XX];
      xp = ixmax[XX] + bsize;
      yn = ixmin[YY];
      yp = ixmax[YY];
      zn = ixmin[ZZ];
      zp = ixmax[ZZ];
      break;

    case YN:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      xp = ixmax[XX] + n;
      yn = ixmin[YY] - bsize;
      yp = ixmin[YY];
      zn = ixmin[ZZ];
      zp = ixmax[ZZ];
      break;

    case YP:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      yn = ixmax[YY];
      yp = ixmax[YY] + bsize;
      zn = ixmin[ZZ];
      zp = ixmax[ZZ];
      break;

    case ZN:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      n  = (ixmin[YY] == fl_xmin[YY]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      yn = ixmin[YY] - n;
      n  = (ixmax[YY] == fl_xmax[YY]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      yp = ixmax[YY] + n;
      zn = ixmin[ZZ] - bsize;
      zp = ixmin[ZZ];
      break;

    case ZP:
      n = (ixmin[XX] == fl_xmin[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xn = ixmin[XX] - n;
      n  = (ixmax[XX] == fl_xmax[XX]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      xp = ixmax[XX] + n;
      n  = (ixmin[YY] == fl_xmin[YY]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      yn = ixmin[YY] - n;
      n  = (ixmax[YY] == fl_xmax[YY]) ? grid->idx() * par.Nbc / 2 :
                                       grid->idx() * par.Nbc_DD / 2;
      yp = ixmax[YY] + n;
      zn = ixmax[ZZ];
      zp = ixmax[ZZ] + bsize;
      break;

    default:
      spdlog::error("{}: {}", "bad direction in 3D C2F", bdata->dir);
  }

#ifdef TEST_C2F
  spdlog::debug(
      "xn={}, xp={}\nyn={}, yp={}\nzn={}, zp={}", xn, xp, yn, yp, zn, zp);
#endif

  int ct  = 0;
  cell *c = grid->FirstPt_All();
  do {
    if ((c->pos[XX] > xn && c->pos[XX] < xp)
        && (c->pos[YY] > yn && c->pos[YY] < yp)
        && (c->pos[ZZ] > zn && c->pos[ZZ] < zp)) {
      bdata->c.push_back(c);
      // rep.printVec("c",c->pos,3);
      // cout <<"adding cell "<<ct<<" to list.\n";
      ct++;
    }
  } while ((c = grid->NextPt_All(*c)) != 0);
#ifdef TEST_C2F
  spdlog::debug("Added {} cells to c2f list", ct);
#endif
}

// ##################################################################
// ##################################################################
