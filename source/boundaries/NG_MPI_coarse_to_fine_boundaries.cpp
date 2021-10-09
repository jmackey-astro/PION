/// \file NG_MPI_coarse_to_fine_boundaries.cpp
/// \brief Class definitions for NG_MPI_coarse_to_fine boundaries
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.09.06 JM: started writing code.

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
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
  vector<struct cgrid> fg;
  MCMD->get_child_grid_info(fg);
  int nchild = fg.size();

#ifdef TEST_C2F
  if (nchild == 0) {
    cout << "COARSE_TO_FINE_SEND: no children.\n";
  }
  else {
    cout << "COARSE_TO_FINE_SEND: " << nchild << " child grids.\n";
  }
#endif

  // Two cases:
  // (1) if level l+1 has a boundary within my domain, then
  // send data to child grids.
  // (2) if level l+1 has a boundary coincident with my boundary,
  // then I need to send data if none of the l+1 domain intersects
  // my grid.  Otherwise another process can do it.

  int cl_ixmin[MAX_DIM], cl_ixmax[MAX_DIM];
  int fl_ixmin[MAX_DIM], fl_ixmax[MAX_DIM];
  int cg_ixmin[MAX_DIM], cg_ixmax[MAX_DIM];
  int fg_ixmin[MAX_DIM], fg_ixmax[MAX_DIM];

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

    if (MCMD->get_myrank() == fg[i].rank) {
      // if child is on my process, do nothing because child grid
      // can grab the data directly using serial C2F code.
#ifdef TEST_C2F
      cout << "C2F_SEND: child " << i << ", ";
      cout << "my rank (" << MCMD->get_myrank() << ") == child rank (";
      cout << fg[i].rank << "), no need to set up ";
      cout << "COARSE_TO_FINE_SEND\n";
#endif
    }
    else {
      // if child is on another process, make a list of cells that
      // we need to send.  One list for each external boundary that
      // we need to update.  Only external boundaries of the whole
      // domain are included.
#ifdef TEST_C2F
      cout << "C2F_SEND: child " << i << ", ";
      cout << "my rank != child rank (" << MCMD->get_myrank();
      cout << ", " << fg[i].rank << ") running parallel ";
      cout << "COARSE_TO_FINE_SEND\n";
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
          cout << "C2F_SEND: child " << i << ", dim " << d << " NEG DIR\n";
          rep.printVec("localxmin", MCMD->LocalXmin, 3);
          rep.printVec("Childxmin", fg[i].Xmin, 3);
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
          cout << "added " << bdata->c.size() << " cells to C2F send el ";
          cout << b->NGsendC2F.size() - 1 << "\n";
#endif
        }
        // if child xmax == its level xmax, but < my level xmax,
        // then we need to send data, so set up a list.
        if ((fg_ixmax[d] == fl_ixmax[d]) && (fg_ixmax[d] < cl_ixmax[d])
            && (fg_ixmax[d] != cg_ixmax[d])) {
#ifdef TEST_C2F
          cout << "C2F_SEND: child " << i << ", dim " << d << " POS DIR\n";
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
          cout << "added " << bdata->c.size() << " cells to C2F send el ";
          cout << b->NGsendC2F.size() - 1 << "\n";
#endif
        }
      }  // loop over dimensions
    }    // if child is not on my process
  }      // loop over child grids

  // We've now dealt with C2F boundaries that are within my domain,
  // so we have to also consider boundaries coincident with my
  // domain boundary, where level l+1 only intersects at the boundary
  // this only applies if nproc>1.
  if (MCMD->get_nproc() > 1) {
    std::vector<std::vector<struct cgrid> > fgngb;
    MCMD->get_level_lp1_ngb_info(fgngb);
    if (fgngb.size() != static_cast<size_t>(2 * par.ndim)) {
      rep.error("l+1 neigbouring grids vector not set up right", fgngb.size());
    }
    for (int d = 0; d < par.ndim; d++) {
      // negative direction
      int dir = 2 * d;
      if (cg_ixmin[d] == fl_ixmax[d]) {
        if (fgngb[dir].size() != 0) {
          // there are neighbouring grids on l+1, so add all of them
          // to the send list.
          for (size_t f = 0; f < fgngb[dir].size(); f++) {
            if (fgngb[dir][f].rank == MCMD->get_myrank()) {
              // if child is on my process, do nothing because child
              // grid can grab the data directly using serial C2F
              // code.
#ifdef TEST_C2F
              cout << "C2F_SEND: ngb " << f << ", dir=" << dir << ", ";
              cout << "my rank (" << MCMD->get_myrank();
              cout << ") == child-ngb rank (";
              cout << fgngb[dir][f].rank << "), no need to set up ";
              cout << "COARSE_TO_FINE_SEND\n";
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
            if (fgngb[dir][f].rank == MCMD->get_myrank()) {
              // if child is on my process, do nothing because child
              // grid can grab the data directly using serial C2F
              // code.
#ifdef TEST_C2F
              cout << "C2F_SEND: ngb " << f << ", dir=" << dir << ", ";
              cout << "my rank (" << MCMD->get_myrank();
              cout << ") == child-ngb rank (";
              cout << fgngb[dir][f].rank << "), no need to set up ";
              cout << "COARSE_TO_FINE_SEND\n";
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
  cout << "MPI C2F SEND: starting from level " << l << " to level " << l + 1
       << "... ";
  if (b->NGsendC2F.size() == 0) {
    cout << "empty send list, so just returning now.\n";
    return 0;
  }
  else
    cout << " send-list size = " << b->NGsendC2F.size();
  cout << endl;
#endif
#ifdef TEST_C2F
  class MCMDcontrol *MCMD = &(par.levels[l].MCMD);
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
    cout << "MPI C2F SEND: odd step, interpolating data in time.\n";
#endif
    double U[par.nvar];
    for (unsigned int ib = 0; ib < b->NGsendC2F.size(); ib++) {
      for (unsigned int c_iter = 0; c_iter < b->NGsendC2F[ib]->c.size();
           c_iter++) {
        cell *c = b->NGsendC2F[ib]->c[c_iter];
        solver->PtoU(c->P, U, par.gamma);
        for (int v = 0; v < par.nvar; v++)
          U[v] += 0.5 * c->dU[v];
        solver->UtoP(U, c->Ph, par.EP.MinTemperature, par.gamma);
#ifdef TEST_INF
        for (int v = 0; v < par.nvar; v++) {
          if (!isfinite(c->Ph[v])) {
            rep.printVec("NAN c->P ", c->P, par.nvar);
            rep.printVec("NAN c->Ph", c->Ph, par.nvar);
          }
        }
#endif
      }  // loop over cells in boundary
    }    // loop over send boundaries
  }      // if not at a full step update

#ifdef TEST_C2F
  cout << "C2F SEND: " << l << " num send=" << b->NGsendC2F.size() << "\n";
#endif

  // loop over send boundaries, pack and send the data.
  for (unsigned int ib = 0; ib < b->NGsendC2F.size(); ib++) {
    //
    // if 1st order accuracy, then just need Ph[]+cell-vol.
    // if 2nd order accuracy, also a slope vector for each dimension
    //
    size_t n_cell = b->NGsendC2F[ib]->c.size();

#ifdef TEST_C2F
    cout << "C2F SEND: " << MCMD->get_myrank() << ", sending ";
    cout << n_cell << " elements to process: ";
    cout << b->NGsendC2F[ib]->rank << "\n";
    cout.flush();
#endif

    size_t n_el = 0;
    if (par.spOOA == OA1)
      n_el = n_cell * (par.nvar + 1 + par.ndim);
    else if (par.spOOA == OA2)
      n_el = n_cell * ((1 + par.ndim) * par.nvar + 1 + par.ndim);
    else
      rep.error("bad spOOA in MPI C2F", par.spOOA);
    pion_flt *buf = new pion_flt[n_el];
    double slope[par.nvar];
    double cpos[par.ndim];

    // loop over cells, add Ph[], cell-vol, slopes to send buffer
    size_t ibuf = 0;
    for (unsigned int c_iter = 0; c_iter < b->NGsendC2F[ib]->c.size();
         c_iter++) {
      cell *c = b->NGsendC2F[ib]->c[c_iter];
      for (int v = 0; v < par.nvar; v++)
        buf[ibuf + v] = c->Ph[v];
      ibuf += par.nvar;
      buf[ibuf] = grid->CellVolume(c, 0);
      ibuf++;
      CI.get_dpos(c, cpos);
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
          solver->SetSlope(c, a, par.nvar, slope, OA2, grid);
          for (int v = 0; v < par.nvar; v++)
            buf[ibuf + v] = slope[v];
          ibuf += par.nvar;
        }  // idim
      }    // if 2nd order

    }  // loop over cells for this send boundary

    if (ibuf != n_el) rep.error("C2F MPI SEND counting", ibuf);

    //
    // Send data using a non-blocking MPI send
    //
    // ostringstream tmp;
    // tmp <<"C2F_"<<MCMD->get_myrank()<<"_to_"<<b->NGsendC2F[ib]->rank;
    string id;  // = tmp.str();
    // Need to add direction to comm-tag because we might be sending
    // more than one boundary to the same process.  Also add level,
    // because it can happen that more than one level sends the same
    // boundary to the same proc.
    // So the tag is BC_MPI_NGC2F_tag + 100*dir + level+1.
    //
    int comm_tag = BC_MPI_NGC2F_tag + 100 * b->NGsendC2F[ib]->dir + l + 1;
#ifdef TEST_C2F
    cout << "BC_update_COARSE_TO_FINE_SEND: l=" << l << ", Sending " << n_el;
    cout << " doubles from proc " << MCMD->get_myrank();
    cout << " to child proc " << b->NGsendC2F[ib]->rank << endl;
#endif
    err +=
        COMM->send_double_data(b->NGsendC2F[ib]->rank, n_el, buf, id, comm_tag);
    if (err) rep.error("Send_C2F send_data failed.", err);
#ifdef TEST_C2F
    cout << "BC_update_COARSE_TO_FINE_SEND: returned with id=" << id;
    cout << endl;
#endif
    // store ID to clear the send later (and delete the MPI temp data)
    NG_C2F_send_list.push_back(id);
    buf = mem.myfree(buf);
  }  // loop over send boundaries
  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::BC_COARSE_TO_FINE_SEND_clear_sends()
{
  for (unsigned int ib = 0; ib < NG_C2F_send_list.size(); ib++) {
#ifdef TEST_C2F
    cout << "C2F_send: clearing send # " << ib + 1 << " of ";
    cout << NG_C2F_send_list.size() << ", id=";
    cout << NG_C2F_send_list[ib] << "...";
    cout.flush();
#endif
    COMM->wait_for_send_to_finish(NG_C2F_send_list[ib]);
#ifdef TEST_C2F
    cout << " ... done!\n";
    cout.flush();
#endif
  }
  NG_C2F_send_list.clear();
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
  class MCMDcontrol *MCMD   = &(par.levels[l].MCMD);
  class GridBaseClass *grid = par.levels[l].grid;

  int cl_ixmin[MAX_DIM], cl_ixmax[MAX_DIM];
  int fl_ixmin[MAX_DIM], fl_ixmax[MAX_DIM];
  int cg_ixmin[MAX_DIM], cg_ixmax[MAX_DIM];
  int fg_ixmin[MAX_DIM], fg_ixmax[MAX_DIM];

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
  MCMD->get_parent_grid_info(&pg);
  MCMD->get_parent_ngb_grid_info(pgngb);
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
    default:
      rep.error("bad direction", b->dir);
  }
  if (!send2pg) {
    // must need to receive from neighbour of parent grid
    if (pgngb[b->dir].rank < 0) rep.error("ngb of parent is null", b->dir);
    CI.get_ipos_vec(pgngb[b->dir].Xmin, cg_ixmin);
    CI.get_ipos_vec(pgngb[b->dir].Xmax, cg_ixmax);
    b->NGrecvC2F_parent = pgngb[b->dir].rank;
#ifdef TEST_C2F
    cout << "C2F MPI Recv, recv from rank" << b->NGrecvC2F_parent;
    cout << ", in outward normal direction " << b->dir << "\n";
#endif
  }
  else {
    b->NGrecvC2F_parent = pg.rank;
  }

  if (MCMD->get_myrank() == b->NGrecvC2F_parent) {
#ifdef TEST_C2F
    cout << "my rank == parent rank, setting up serial ";
    cout << "COARSE_TO_FINE_RECV\n";
#endif
    int err = NG_coarse_to_fine_bc::BC_assign_COARSE_TO_FINE(
        par, grid, b, par.levels[l].parent);
    rep.errorTest("serial C2F BC setup", 0, err);
  }
  else {
#ifdef TEST_C2F
    cout << "my rank != parent rank, so will need MPI for ";
    cout << "COARSE_TO_FINE_RECV";
    cout << "me=" << MCMD->get_myrank();
    cout << ", parent=" << b->NGrecvC2F_parent << "\n";
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
          b->NGrecvC2F[ic].push_back(grid->NextPt(c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp, YP));
          ic++;
        }
        else if (c->pos[YY] - row_y == idx) {
          // odd-numbered row, already dealt with.
          continue;
        }
        else
          rep.error("error in 2d logic C2FR_setup", c->pos[YY] - row_y);
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
          b->NGrecvC2F[ic].push_back(grid->NextPt(c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp, YP));
          c    = grid->NextPt(c, ZP);
          temp = grid->NextPt(temp, ZP);
          b->NGrecvC2F[ic].push_back(c);
          b->NGrecvC2F[ic].push_back(temp);
          b->NGrecvC2F[ic].push_back(grid->NextPt(c, YP));
          b->NGrecvC2F[ic].push_back(grid->NextPt(temp, YP));

          ic++;
        }
        else if (c->pos[YY] - row_y == idx || c->pos[ZZ] - row_z == idx) {
          // odd-numbered row/plane, already added cells.
          continue;
        }
        else
          rep.error("error in 3d logic C2FR_setup", c->pos[YY] - row_y);
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
  cout << "C2F_MPI: receiving boundary data to level ";
  cout << l << ", updating boundary dir = " << b->dir << endl;
#endif
  int err                   = 0;
  class MCMDcontrol *MCMD   = &(par.levels[l].MCMD);
  class GridBaseClass *grid = par.levels[l].grid;

  if (MCMD->get_myrank() == b->NGrecvC2F_parent) {
#ifdef TEST_C2F
    cout << "my rank == parent rank, calling serial ";
    cout << "COARSE_TO_FINE" << endl;
#endif
    NG_coarse_to_fine_bc::BC_update_COARSE_TO_FINE(par, solver, l, b, step);
  }
  else {
#ifdef TEST_C2F
    cout << "my rank != parent rank, so updating ";
    cout << "COARSE_TO_FINE_RECV: ";
    cout << "me=" << MCMD->get_myrank();
    cout << ", parent=" << b->NGrecvC2F_parent << endl;
#endif

    // receive data.
    string recv_id;
    int recv_tag  = -1;
    int from_rank = -1;
    int comm_tag  = BC_MPI_NGC2F_tag + 100 * b->dir + l;
#ifdef TEST_C2F
    if (1 == 1) {
      cout << "BC_update_COARSE_TO_FINE_RECV: looking for data tag ";
      cout << comm_tag << endl;
    }
#endif
    err = COMM->look_for_data_to_receive(
        &from_rank, recv_id, &recv_tag, comm_tag, COMM_DOUBLEDATA);
    if (err) rep.error("look for double data failed", err);
#ifdef TEST_C2F
    if (1 == 1) {
      cout << "BC_update_COARSE_TO_FINE_RECV: found data from rank ";
      cout << from_rank << ", with tag " << recv_tag << " and id ";
      cout << recv_id << ".  Looked for comm_tag=" << comm_tag << endl;
    }
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
      rep.error("bad spOOA in MPI C2F", par.spOOA);
    pion_flt *buf = 0;
    buf           = mem.myalloc(buf, n_el);
#ifdef TEST_C2F
    if (1 == 1) {
      cout << "BC_update_COARSE_TO_FINE_RECV: get " << n_cell;
      cout << " cells, and " << n_el << " doubles.\n";
    }
#endif
    //
    // Receive data into buffer.  Data stored for each coarse cell:
    // Ph[nv],cellvol,cellpos[nd],slopeX[nv],slopeY[nv],slopeZ[nv]
    //
    err = COMM->receive_double_data(from_rank, recv_tag, recv_id, n_el, buf);
    if (err) rep.error("(BC_update_C2F_RECV) getdata failed", err);

    //
    // Get Ph and slopes from coarse cells, and interpolate into
    // groups of fine cells that are children of each coarse cell.
    //
    size_t ibuf = 0;
    if (par.spOOA == OA1) {
      // 1st order, so no averaging.  Fine cells have exactly the
      // same state as the coarse one.
      double Ph[par.nvar];
      double cpos[par.ndim];
#ifdef TEST_C2F
      double off[par.ndim];
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
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol = 0.0;
        double sx[par.nvar];
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
              par, grid, solver, Ph, c_vol, sx, f[0], f[1]);
        }  // loop over coarse cells
      }    // if 1D
      else if (par.ndim == 2) {
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol = 0.0;
        int ipos[par.ndim];
        ipos[0] = 0;
        ipos[1] = 0;
        double sx[par.nvar], sy[par.nvar];
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
              par, grid, solver, Ph, ipos, c_vol, sx, sy, f[0], f[1], f[2],
              f[3]);

        }  // loop over coarse cells
      }    // if 2D
      else {
        double Ph[par.nvar];
        double cpos[par.ndim], c_vol = 0.0;
#ifdef TEST_C2F
        double cpos2[par.ndim];
#endif
        int ipos[par.ndim];
        double sx[par.nvar], sy[par.nvar], sz[par.nvar];
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
              par, grid, solver, Ph, ipos, c_vol, sx, sy, sz, fch);
        }  // loop over coarse cells
      }    // if 3D
    }      // if 2nd order accurate

    buf = mem.myfree(buf);
    buf = 0;
  }  // if parent proc is different to my proc.

#ifdef TEST_C2F
  cout << "NG_MPI_C2F_bc::BC_update_COARSE_TO_FINE_RECV() done\n";
  cout.flush();
#endif
  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to coarse-level grid
    struct c2f *bdata,          ///< pointer to list of cells
    int *ixmin,                 ///< child grid xmin (integer)
    int *ixmax,                 ///< child grid xmax (integer)
    const int lf,               ///< level of fine grid.
    const int *fl_xmin,         ///< level xmin of fine grid.
    const int *fl_xmax          ///< level xmax of fine grid.
)
{
  // In XN,XP direction we add cells with faces that touch the fine-
  // level grid from the outside.
  // In YN,YP direction we also add corner data
  // In ZN,ZP direction we also add edge data
  //
  // easier to have a function for different grid dimensions.
#ifdef TEST_C2F
  cout << "Adding cells to C2F Send list: ";
#endif
  if (par.ndim == 1) {
#ifdef TEST_C2F
    cout << "1D \n";
#endif
    add_cells_to_C2F_send_list_1D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  else if (par.ndim == 2) {
#ifdef TEST_C2F
    cout << "2D \n";
#endif
    add_cells_to_C2F_send_list_2D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  else {
#ifdef TEST_C2F
    cout << "3D \n";
#endif
    add_cells_to_C2F_send_list_3D(
        par, grid, bdata, ixmin, ixmax, lf, fl_xmin, fl_xmax);
  }
  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_1D(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to coarse-level grid
    struct c2f *bdata,          ///< pointer to list of cells
    int *ixmin,                 ///< child grid xmin (integer)
    int *ixmax,                 ///< child grid xmax (integer)
    const int lf,               ///< level of fine grid.
    const int *fl_xmin,         ///< level xmin of fine grid.
    const int *fl_xmax          ///< level xmax of fine grid.
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
      rep.error("bad direction in 1D C2F", bdata->dir);
  }

  cell *c = grid->FirstPt_All();
  do {
    if (c->pos[XX] > xn && c->pos[XX] < xp) bdata->c.push_back(c);
  } while ((c = grid->NextPt_All(c)) != 0);

  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_2D(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to coarse-level grid
    struct c2f *bdata,          ///< pointer to list of cells
    int *ixmin,                 ///< child grid xmin (integer)
    int *ixmax,                 ///< child grid xmax (integer)
    const int lf,               ///< level of fine grid.
    const int *fl_xmin,         ///< level xmin of fine grid.
    const int *fl_xmax          ///< level xmax of fine grid.
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
      rep.error("bad direction in 2D C2F", bdata->dir);
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
  } while ((c = grid->NextPt_All(c)) != 0);
#ifdef TEST_C2F
  cout << "add_cells_to_C2F_send_list_2D: added " << ct << " cells.\n";
#endif

  return;
}

// ##################################################################
// ##################################################################

void NG_MPI_coarse_to_fine_bc::add_cells_to_C2F_send_list_3D(
    class SimParams &par,       ///< pointer to simulation parameters
    class GridBaseClass *grid,  ///< pointer to coarse-level grid
    struct c2f *bdata,          ///< pointer to list of cells
    int *ixmin,                 ///< child grid xmin (integer)
    int *ixmax,                 ///< child grid xmax (integer)
    const int lf,               ///< level of fine grid.
    const int *fl_xmin,         ///< level xmin of fine grid.
    const int *fl_xmax          ///< level xmax of fine grid.
)
{
#ifdef TEST_C2F
  rep.printVec("C2F Setup Send: ixmin", ixmin, 3);
  rep.printVec("C2F Setup Send: ixmax", ixmax, 3);
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
      rep.error("bad direction in 3D C2F", bdata->dir);
  }

#ifdef TEST_C2F
  cout << "xn=" << xn << ", xp=" << xp << "\n";
  cout << "yn=" << yn << ", yp=" << yp << "\n";
  cout << "zn=" << zn << ", zp=" << zp << "\n";
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
  } while ((c = grid->NextPt_All(c)) != 0);
#ifdef TEST_C2F
  cout << "Added " << ct << " cells to c2f list\n";
#endif

  return;
}

// ##################################################################
// ##################################################################
