

#include "sim_params.h"
#include "sub_domain.h"
#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

#include <sstream>

using namespace std;

//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------

// ##################################################################
// ##################################################################


void Sub_domain::get_domain_coordinates(const int r, int *arr) const
{
  MPI_Cart_coords(cart_comm, r, m_ndim, arr);
}

// ##################################################################
// ##################################################################


int Sub_domain::get_rank_from_grid_location(
    const class SimParams &par,  ///< simulation parameters
    const vector<double> &loc,   ///< location sought
    const int l                  ///< grid level
    ) const
{
  int grid_coordinates[par.ndim];
  for (int i = 0; i < par.ndim; i++) {
    grid_coordinates[i] = static_cast<int>(floor(
        num_subdomains[i] * (loc[i] - par.levels[l].Xmin[i])
        / par.levels[l].Range[i]));
  }
  /* check validity (0 <= x <= num_subdomains[] - 1) */
  for (int i = 0; i < par.ndim; i++) {
    if ((grid_coordinates[i] * (grid_coordinates[i] + 1 - num_subdomains[i])
         > 0)) {
      return -1;
    }
  }
  /* MPI command to get rank from index in Cartesian grid */
  int proc;
  MPI_Cart_rank(cart_comm, grid_coordinates, &proc);
  return proc;
}


// ##################################################################
// ##################################################################


void Sub_domain::get_num_subdomains(int *tmp) const
{
  for (int i = 0; i < MAX_DIM; i++)
    tmp[i] = num_subdomains[i];
  return;
}


// ##################################################################
// ##################################################################


void Sub_domain::get_parent_grid_info(struct cgrid *cg) const
{
  cg->rank = pgrid.rank;
  for (int i = 0; i < MAX_DIM; i++)
    cg->Xmin[i] = pgrid.Xmin[i];
  for (int i = 0; i < MAX_DIM; i++)
    cg->Xmax[i] = pgrid.Xmax[i];
}

// ##################################################################
// ##################################################################


void Sub_domain::get_parent_ngb_grid_info(vector<struct cgrid> &pgngb) const
{
  pgngb.resize(pgrid_ngb.size());
  for (size_t iter = 0; iter < pgrid_ngb.size(); iter++) {
    pgngb[iter].rank = pgrid_ngb[iter].rank;
    for (int i = 0; i < MAX_DIM; i++)
      pgngb[iter].Xmin[i] = pgrid_ngb[iter].Xmin[i];
    for (int i = 0; i < MAX_DIM; i++)
      pgngb[iter].Xmax[i] = pgrid_ngb[iter].Xmax[i];
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::get_child_grid_info(vector<struct cgrid> &cg) const
{
  cg.resize(child_procs.size());
  for (size_t iter = 0; iter < child_procs.size(); iter++) {
    cg[iter].rank = child_procs[iter].rank;
    for (int i = 0; i < MAX_DIM; i++)
      cg[iter].Xmin[i] = child_procs[iter].Xmin[i];
    for (int i = 0; i < MAX_DIM; i++)
      cg[iter].Xmax[i] = child_procs[iter].Xmax[i];
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::get_level_lp1_ngb_info(
    vector<vector<struct cgrid> > &cgngb) const
{
  cgngb.resize(cgrid_ngb.size());
  for (size_t p = 0; p < cgrid_ngb.size(); p++) {
    cgngb[p].resize(cgrid_ngb[p].size());
    for (size_t q = 0; q < cgrid_ngb[p].size(); q++) {
      cgngb[p][q].rank = cgrid_ngb[p][q].rank;
      for (int i = 0; i < MAX_DIM; i++)
        cgngb[p][q].Xmin[i] = cgrid_ngb[p][q].Xmin[i];
      for (int i = 0; i < MAX_DIM; i++)
        cgngb[p][q].Xmax[i] = cgrid_ngb[p][q].Xmax[i];
    }
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::gather_ncells(int *recv_buffer, const int &root) const
{
  MPI_Gather(&Ncell, 1, MPI_INT, recv_buffer, 1, MPI_INT, root, cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::allgather_ncells(vector<int> &ncells_list) const
{
  ncells_list.resize(nproc);
  MPI_Allgather(&Ncell, 1, MPI_INT, &ncells_list[0], 1, MPI_INT, cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::gather_extents(double *recv_buffer, const int &root) const
{
  double extents[2 * m_ndim];
  for (int i = 0; i < m_ndim; i++) {
    extents[i]          = Xmin[i];
    extents[m_ndim + i] = Xmax[i];
  }
  MPI_Gather(
      &extents, 2 * m_ndim, MPI_DOUBLE, recv_buffer, 2 * m_ndim, MPI_DOUBLE,
      root, cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::gather_offsets(
    vector<int> &offsets_list, const int &root) const
{
  offsets_list.resize(nproc * m_ndim);
  MPI_Gather(
      &offsets[0], m_ndim, MPI_INT, &offsets_list[0], m_ndim, MPI_INT, root,
      cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::gather_directional_Ncells(
    vector<int> &directional_Ncells_list, const int &root) const
{
  directional_Ncells_list.resize(nproc * m_ndim);
  MPI_Gather(
      &directional_Ncells[0], m_ndim, MPI_INT, &directional_Ncells_list[0],
      m_ndim, MPI_INT, root, cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::gather_abutting_domains(
    vector<int> &abutting_domains_list,
    int *rank_displacements_in_abutting_domains_list,
    int *num_abutting_domains_by_rank,
    const int &root)
{
  create_abutting_domains_list();
  int num_abutting_domains = abutting_domains.size();
  MPI_Gather(
      &num_abutting_domains, 1, MPI_INT, &num_abutting_domains_by_rank[0], 1,
      MPI_INT, root, cart_comm);

  /*
   * calculate running sum of abutting domains by rank and displacements for
   * storing local abutting domains lists in a global list
   */
  if (myrank == 0) {
    rank_displacements_in_abutting_domains_list[0] = 0;
    int i;
    for (i = 0; i < nproc - 1; i++) {
      rank_displacements_in_abutting_domains_list[i + 1] =
          rank_displacements_in_abutting_domains_list[i]
          + num_abutting_domains_by_rank[i];
    }
    abutting_domains_list.resize(
        rank_displacements_in_abutting_domains_list[i]
        + num_abutting_domains_by_rank[i]);
  }
  MPI_Gatherv(
      &abutting_domains[0], abutting_domains.size(), MPI_INT,
      &abutting_domains_list[0], &num_abutting_domains_by_rank[0],
      rank_displacements_in_abutting_domains_list, MPI_INT, root, cart_comm);
}

// ##################################################################
// ##################################################################


void Sub_domain::print_grid(
    std::vector<int> &coords, int current_dimension) const
{
  stringstream grid;
  for (coords[current_dimension] = 0;
       coords[current_dimension] < num_subdomains[current_dimension];
       coords[current_dimension]++) {
    if (current_dimension == m_ndim - 1) {
      int proc;
      MPI_Cart_rank(cart_comm, &coords[0], &proc);
      grid << proc << " ";
    }
    else {
      print_grid(coords, current_dimension + 1);
    }
  }
  spdlog::debug("{}", grid.str());
}
