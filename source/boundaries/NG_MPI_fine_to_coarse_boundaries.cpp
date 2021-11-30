/// \file NG_MPI_fine_to_coarse_boundaries.cpp
/// \brief Class definitions for NG_MPI_fine_to_coarse boundaries
/// \author Jonathan Mackey
///
/// \section Description
/// Coarse-to-Fine boundary updates send grid data from a coarse grid
/// to update the external boundaries of a fine grid contained within
/// the coarse grid.  This version is parallelised with MPI.
///
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include "boundaries/NG_MPI_fine_to_coarse_boundaries.h"
#include "tools/mem_manage.h"
using namespace std;

//#define TEST_MPI_NG_F2C
//#define NG_F2C_POS

// ##################################################################
// ##################################################################

int NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_SEND(
    class SimParams &par,  ///< pointer to simulation parameters
    const int l,           ///< level of this grid.
    boundary_data *b       ///< boundary data
)
{
#ifdef TEST_MPI_NG_F2C
  spdlog::debug("***F2C F2C F2C*** Adding F2C SEND for l={}", l);
#endif
  // Check if parent grid is on my MPI process
  int pproc                    = par.levels[l].sub_domain.get_parent_proc();
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);

  //
  // If we are doing raytracing, then also send the column densities
  // from fine to coarse grids.  Here we set up the number of them.
  //
  struct rad_src_info *s;
  F2C_Nxd = 0;
  F2C_tauoff.resize(par.RS.Nsources);
  for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
    s                = &(par.RS.sources[isrc]);
    F2C_tauoff[isrc] = F2C_Nxd;
    F2C_Nxd += 2 * s->NTau;  // Need col2cell and cell_col for each Tau.
#ifdef TEST_MPI_NG_F2C
    spdlog::debug(
        "F2C_MPI_BC: RT Source {}: adding {} cols, running total = {}", isrc,
        2 * s->NTau, F2C_Nxd);
#endif
  }

  // If so, do nothing because the parent grid will call the
  // serial NG setup function, which just grabs the data from
  // this grid.
  if (sub_domain->get_myrank() == pproc) {
#ifdef TEST_MPI_NG_F2C
    spdlog::debug(
        "my rank == parent rank, not setting up FINE_TO_COARSE_SEND\n");
#endif
  }
  else {
    // Else, make a vector of structs with a list of cells contained
    // in each coarse cell (2,4,or 8) of the parent grid.  This grid
    // is (by design) entirely contained within the parent.
    class GridBaseClass *grid = par.levels[l].grid;
    int nc                    = 1;  // number of fine cells per coarse cell
    for (int i = 0; i < par.ndim; i++)
      nc *= 2;
    int nel = grid->Ncell() / nc;  // number of coarse cells
    // initialise the "avg" vector to hold the average state.
    b->avg.resize(nel);
    for (int v = 0; v < nel; v++) {
      b->avg[v].c.resize(nc);
      b->avg[v].avg_state =
          mem.myalloc(b->avg[v].avg_state, par.nvar + F2C_Nxd);
    }
    // add cells from this grid to the avg vectors.
    add_cells_to_avg(par.ndim, grid, nel, b->avg);
#ifdef TEST_MPI_NG_F2C
    spdlog::debug("***F2C added {} avg cells to this grid\n", nel);
#endif

  }  // if parent is not on a different MPI process

  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_fine_to_coarse_bc::BC_update_FINE_TO_COARSE_SEND(
    class SimParams &par,          ///< pointer to simulation parameters
    class FV_solver_base *solver,  ///< pointer to equations
    const int l,                   ///< level in the NG grid structure
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{

  // Check if parent grid is on my MPI process
  int pproc                    = par.levels[l].sub_domain.get_parent_proc();
  int err                      = 0;
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);

  // If so, do nothing because the parent grid will call the
  // serial NG setup function, which just grabs the data from
  // this grid.
  if (sub_domain->get_myrank() == pproc) {
#ifdef TEST_MPI_NG_F2C
    spdlog::debug(
        "my rank == parent rank, no need to send anything FINE_TO_COARSE_SEND\n");
#endif
    return 0;
  }

  // else we have to send data to another MPI process:
  class GridBaseClass *grid = par.levels[l].grid;
  int nc                    = b->avg[0].c.size();
  int nel                   = b->avg.size();

  // data to send will be ordered as position,conserved-var,X-data
  // for each averaged cell.  Position only needed for testing.
  pion_flt *data = 0;
#ifdef NG_F2C_POS
  data = mem.myalloc(data, nel * (par.nvar + F2C_Nxd + par.ndim));
#else
  data = mem.myalloc(data, nel * (par.nvar + F2C_Nxd));
#endif

  // loop through avg vector and add cells and positions to
  // data array.
  int v  = 0;
  int nv = par.nvar + F2C_Nxd;  // conserved vector and optical depths
  std::vector<double> cd(nv);
  size_t ct = 0;
  for (v = 0; v < nel; v++) {
    for (int j = 0; j < nv; j++)
      cd[j] = 0.0;
    average_cells(
        par, solver, grid, nc, b->avg[v].c, b->avg[v].cpos, cd.data());
#ifdef NG_F2C_POS
    for (int i = 0; i < par.ndim; i++)
      data[ct + i] = b->avg[v].cpos[i];
    ct += par.ndim;
#endif
    for (int i = 0; i < nv; i++)
      data[ct + i] = cd[i];
    ct += nv;

  }  // go to next avg element.

  //
  // Send data using a non-blocking MPI send
  //
  string id;
  int comm_tag = BC_MPI_NGF2C_tag + l;
  // id <<"F2C_"<<sub_domain->get_myrank()<<"_to_"<<pproc;
#ifdef TEST_MPI_NG_F2C
  spdlog::debug(
      "BC_update_FINE_TO_COARSE_SEND: Sending {} doubles from proc {} to parent proc {}, tag={}",
      ct, sub_domain->get_myrank(), pproc, comm_tag);
#endif
  err += sub_domain->send_double_data(pproc, ct, data, id, comm_tag);
  if (err) spdlog::error("{}: {}", "Send_F2C send_data failed.", err);
#ifdef TEST_MPI_NG_F2C
  spdlog::debug("BC_update_FINE_TO_COARSE_SEND: returned with id={}", id);
#endif

  // store ID to clear the send later (and delete the MPI temp data)
  sub_domain->NG_F2C_send_list.push_back(id);
#ifdef TEST_MPI_NG_F2C
  spdlog::debug(
      "F2C_Send: id=[ {} ]  size={}", id, sub_domain->NG_F2C_send_list.size());
#endif
  data = mem.myfree(data);

  return 0;
}

// ##################################################################
// ##################################################################

void NG_MPI_fine_to_coarse_bc::BC_FINE_TO_COARSE_SEND_clear_sends(
    class Sub_domain &sub_domain)
{
  for (unsigned int i = 0; i < sub_domain.NG_F2C_send_list.size(); i++) {
#ifdef TEST_MPI_NG_F2C
    spdlog::debug(
        "F2C_send: clearing send # {} of {}, id={}...", i + 1,
        sub_domain.NG_F2C_send_list.size(), sub_domain.NG_F2C_send_list[i]);
#endif
    // cout <<"waiting for send : "<<sub_domain.NG_F2C_send_list[i]<<"\n";
    sub_domain.wait_for_send_to_finish(sub_domain.NG_F2C_send_list[i]);
#ifdef TEST_MPI_NG_F2C
    spdlog::debug(" ... done!");
#endif
    // cout <<"cleared F2C send "<<i<<"\n";
  }
  sub_domain.NG_F2C_send_list.clear();
  // cout <<"clear F2C sends returning\n";
  return;
}

// ##################################################################
// ##################################################################

int NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_RECV(
    class SimParams &par,  ///< pointer to simulation parameters
    const int l,           ///< level of this grid.
    boundary_data *b       ///< boundary data
)
{
#ifdef TEST_MPI_NG_F2C
  spdlog::info(
      "NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_RECV(): starting... ");
#endif

  // Check if child grids exist or are on my MPI process
  int err                      = 0;
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
  vector<struct cgrid> cg;
  sub_domain->get_child_grid_info(cg);
  int nchild = cg.size();
#ifdef TEST_MPI_NG_F2C
  spdlog::debug("level = {}, nchild={}", l, nchild);
#endif
  b->NGrecvF2C.resize(nchild);
  b->NGrecvF2C_ranks.resize(nchild);
  for (int i = 0; i < nchild; i++) {
    b->NGrecvF2C_ranks[i] = cg[i].rank;
  }

  // loop over children:
  for (int i = 0; i < nchild; i++) {
#ifdef TEST_MPI_NG_F2C
    spdlog::debug("child = {}, processing.", i);
#endif

    // get dimensions of child grid from struct
    std::array<int, MAX_DIM> ixmin, ixmax;
    CI.get_ipos_vec(cg[i].Xmin, ixmin);
    CI.get_ipos_vec(cg[i].Xmax, ixmax);

#ifdef TEST_MPI_NG_F2C
    spdlog::debug("F2C MPI: child xmin : {}", ixmin);
    spdlog::debug("F2C MPI: child xmax : {}", ixmax);
#endif

    class GridBaseClass *grid = par.levels[l].grid;
    cell *c                   = grid->FirstPt();
    size_t ct                 = 0;
    do {
      bool ongrid = true;
      for (int j = 0; j < par.ndim; j++) {
        if (c->pos[j] < ixmin[j] || c->pos[j] > ixmax[j]) ongrid = false;
      }
      if (ongrid) {
        // c->isbd = true;
        // c->isleaf = false;
#ifdef TEST_MPI_NG_F2C
        spdlog::debug("cell over fine grid : {}", c->pos);
#endif
        b->NGrecvF2C[i].push_back(c);
        ct++;
      }
    } while ((c = grid->NextPt(c)) != 0);

#ifdef TEST_MPI_NG_F2C
    spdlog::debug(
        "i={} added {} cells to boundary, Ncell/2^ndim={}", i, ct,
        grid->Ncell() / pow(2, par.ndim));
#endif

    if (sub_domain->get_myrank() == b->NGrecvF2C_ranks[i]) {
      // If child is on my grid call serial version that grabs data
      // directly from the child grid, and puts them onto the local
      // grid.
#ifdef TEST_MPI_NG_F2C
      spdlog::debug(
          "my rank == child rank, calling serial FINE_TO_COARSE\nlevel {}: my grid={}, child grid={}",
          l, par.levels[l].grid, par.levels[l].child);
#endif
      err = NG_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE(
          par, par.levels[l].grid, b, par.levels[l].child, i);
      if (0 != err)
        spdlog::error(
            "{}: Expected {} but got {}",
            "BC_assign_FINE_TO_COARSE_RECV serial", 0, err);
    }
    else {
      //
      // else we have to create a list of cells for each
      // child, and add them to b->NGrecv[i]
      // This is just done above.
      //
#ifdef TEST_MPI_NG_F2C
      spdlog::info(
          "my rank != child rank, running parallel FINE_TO_COARSE_RECV");
#endif
    }  // if child is on a different MPI process

  }  // loop over child grids

  return 0;
}

// ##################################################################
// ##################################################################

int NG_MPI_fine_to_coarse_bc::BC_update_FINE_TO_COARSE_RECV(
    class SimParams &par,          ///< pointer to simulation parameters
    class FV_solver_base *solver,  ///< pointer to equations
    const int l,                   ///< level in the NG grid structure
    struct boundary_data *b,
    const int cstep,
    const int maxstep)
{
  // Check if child grids exist or are on my MPI process
  int err                      = 0;
  class Sub_domain *sub_domain = &(par.levels[l].sub_domain);
  vector<struct cgrid> cg;
  sub_domain->get_child_grid_info(cg);
  int nchild = b->NGrecvF2C.size();
  if (static_cast<int>(cg.size()) != nchild)
    spdlog::error(
        "{}: {}", "BC_update_FINE_TO_COARSE_RECV: bad size",
        nchild - cg.size());
  int count = 0;

  // loop over children twice, once for child grids that are on my
  // MPI process, and once for grids that on other processes
  for (int i = 0; i < nchild; i++) {
#ifdef TEST_MPI_NG_F2C
    spdlog::debug("F2C_RECV: child {}, receiving...", i);
#endif

    if (sub_domain->get_myrank() == cg[i].rank) {
      // If child is on my grid call serial version and grab data
      // directly from the child grid.
#ifdef TEST_MPI_NG_F2C
      spdlog::info(
          "my rank == child rank, calling serial FINE_TO_COARSE update");
#endif
      err = NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
          par, solver, l, b, i, cstep, maxstep);
      if (0 != err)
        spdlog::error(
            "{}: Expected {} but got {}",
            "BC_update_FINE_TO_COARSE_RECV serial", 0, err);
      count++;
    }

    else {
#ifdef TEST_MPI_NG_F2C
      spdlog::info(
          "my rank != child rank\nF2C_RECV: child {}, receiving via MPI...\nmy rank != child rank, running parallel FINE_TO_COARSE_RECV\n",
          count);
#endif
      //
      // receive data.
      //
      string recv_id;
      int recv_tag  = -1;
      int from_rank = cg[i].rank;
      int comm_tag  = BC_MPI_NGF2C_tag + l + 1;
      err           = sub_domain->look_for_data_to_receive(
          &from_rank,  // rank of sender (output)
          recv_id,     // identifier for receive (output).
          &recv_tag,   // comm_tag associated with data (output)
          comm_tag,
          COMM_DOUBLEDATA  // type of data we want.
      );
      if (err) spdlog::error("{}: {}", "look for double data failed", err);
#ifdef TEST_MPI_NG_F2C
      spdlog::debug(
          "BC_update_FINE_TO_COARSE_RECV: found data from rank {}", from_rank);
#endif

      // check association of data with the child grid:
      if (b->NGrecvF2C_ranks[i] == from_rank) {
#ifdef TEST_MPI_NG_F2C
        spdlog::debug(
            "getting data from rank {} associated with child grid {}",
            from_rank, i);
#endif
      }
      else {
        spdlog::error("{}: {}", "data from wrong child", from_rank);
      }
#ifdef TEST_MPI_NG_F2C
      spdlog::debug(
          "Found data to receive from {} associated with child grid {}",
          from_rank, i);
#endif

      // receive the data
      size_t nel = b->NGrecvF2C[i].size();

#ifdef NG_F2C_POS
      size_t ct = nel * (par.nvar + F2C_Nxd + par.ndim);
#else
      size_t ct = nel * (par.nvar + F2C_Nxd);
#endif

      pion_flt *buf = 0;
      buf           = mem.myalloc(buf, ct);
#ifdef TEST_MPI_NG_F2C
      spdlog::debug("BC_update_FINE_TO_COARSE_RECV: get {} cells.\n", nel);
#endif
      //
      // Receive data into buffer.
      //
      err = sub_domain->receive_double_data(
          from_rank, recv_tag, recv_id, ct, buf);
      if (err)
        spdlog::error("{}: {}", "(BC_update_F2C_RECV) getdata failed", err);

      // Go through list of cells in child i, and put the received
      // data onto these cells.
      list<cell *>::iterator c_iter = b->NGrecvF2C[i].begin();
      cell *c                       = 0;
#ifdef NG_F2C_POS
      std::array<pion_flt, MAX_DIM> pos;
#endif
      std::vector<pion_flt> prim(par.nvar);
      size_t i_el = 0;
      for (c_iter = b->NGrecvF2C[i].begin(); c_iter != b->NGrecvF2C[i].end();
           ++c_iter) {
        c = (*c_iter);
#ifdef NG_F2C_POS
        CI.get_dpos(c, pos);
        bool errors = false;
        for (int v = 0; v < par.ndim; v++) {
          if (!pconst.equalD(pos[v], buf[i_el + v])) errors = true;
        }
        if (errors) {
          spdlog::error("ERROR in positions! ");
          spdlog::debug("cell pos : {}", pos);
          spdlog::debug("recv pos : {}", &(buf[i_el]));
        }
        i_el += par.ndim;
#endif
        solver->UtoP(
            &(buf[i_el]), prim.data(), par.EP.MinTemperature, par.gamma);
        i_el += par.nvar;

#ifdef TEST_INF
        for (int v = 0; v < par.nvar; v++)
          if (!isfinite(prim[v])) spdlog::error("{}: {}", "NAN F2C recv", v);
#endif

        for (int v = 0; v < par.nvar; v++)
          c->Ph[v] = prim[v];
        if (cstep == maxstep) {
          for (int v = 0; v < par.nvar; v++)
            c->P[v] = c->Ph[v];
        }

        // set coarse cell optical depths for any radiation sources by
        // taking values received (see get_F2C_Tau() for ordering)
        struct rad_src_info *s;
        int off;
        for (int isrc = 0; isrc < par.RS.Nsources; isrc++) {
          s   = &(par.RS.sources[isrc]);
          off = i_el + F2C_tauoff[isrc];
          CI.set_col(c, s->id, &(buf[off]));
          CI.set_cell_col(c, s->id, &(buf[off + s->NTau]));
        }
        i_el += F2C_Nxd;

      }  // loop over cells

#ifdef TEST_MPI_NG_F2C
      spdlog::debug(
          "(BC_update_F2C_RECV) i_el={} of {} total elements.\n", i_el, ct);
#endif
      buf = mem.myfree(buf);
      count++;

    }  // else child is not on my proc
  }    // loop until we got through all child grids.

  return 0;
}

// ##################################################################
// ##################################################################
