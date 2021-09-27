/// \file MCMD_control.cpp
///
/// \brief Defines class for controlling Multi-Core-Multi-Domain
///        simulations.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.27 JM: moved from sim_control_MPI.cpp
/// - 2016.02.02 JM: Added option to decompose only along one axis.
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "MCMD_control.h"
#include <algorithm>
#include <tools/command_line_interface.h>
#include <tools/mem_manage.h>
#include <tools/reporting.h>

using namespace std;

//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------

// ##################################################################
// ##################################################################


MCMDcontrol::MCMDcontrol()
{
  nproc      = -1;
  myrank     = -1;
  LocalNcell = -1;
  periodic   = vector<int>(MAX_DIM, 0);
  for (int i = 0; i < MAX_DIM; i++) {
    LocalNG[i] = offsets[i] = ix[i] = -1;
    nx[i]                           = 0;
    LocalXmin[i] = LocalXmax[i] = LocalRange[i] = -1.e99;
  }
  pgrid.rank = -1;

  ReadSingleFile =
      true;  ///< If the ICs are in a single file, set this to true.
  WriteSingleFile = false;  ///< If all processes to write to one file, set this
  WriteFullImage =
      false;  ///< If multiple fits files, each is the full domain size.
}

// ##################################################################
// ##################################################################


MCMDcontrol::~MCMDcontrol()
{
  /* TODO: can't free because COMM object is deleted before this one */
  // if (cart_comm) {
  //    MPI_Comm_free(&cart_comm);
  //}
}

// ##################################################################
// ##################################################################


void MCMDcontrol::calculate_process_topology(vector<float> &&ratios)
{
  /* account for hardcoded decomposition (i.e, nx[i] already non-zero) */
  int free_processes = nproc;
  for (int i = 0; i < m_ndim; i++) {
    if (nx[i] > 0) {
      free_processes /= nx[i];
      ratios[i] = 0;
    }
    else {
      nx[i] = 1;
    }
  }

  for (; free_processes > 1; free_processes /= 2) {
    const int max_index =
        max_element(ratios.begin(), ratios.end()) - ratios.begin();
    ratios[max_index] /= 2;
    nx[max_index] *= 2;
  }
}

// ##################################################################
// ##################################################################

void MCMDcontrol::print_grid(
    std::vector<int> &coords, int current_dimension) const
{
  for (coords[current_dimension] = 0;
       coords[current_dimension] < nx[current_dimension];
       coords[current_dimension]++) {
    if (current_dimension == m_ndim - 1) {
      int proc;
      MPI_Cart_rank(cart_comm, &coords[0], &proc);
      cout << proc << " ";
    }
    else {
      print_grid(coords, current_dimension + 1);
    }
  }
  cout << "\n";
}

// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
    const int &ndim,          ///< number of dimensions
    const class level &level  ///< parameters for NG grid level
)
{
#ifndef NDEBUG
  cout << "---MCMDcontrol::decomposeDomain() decomposing domain. ";
  cout << " Nproc=" << nproc << ", myrank=" << myrank << "\n";
#endif

  m_ndim = ndim;

  /*
   * the following sets nx to application topology aware values
   * dimensions of nx not set to 0 are restricted to the previously held value
   */
  vector<float> ratios(begin(level.Range), begin(level.Range) + m_ndim);
  calculate_process_topology(move(ratios));

  /* TODO: reordering of ranks is temporarily disabled so old communicator works
   */
  MPI_Cart_create(MPI_COMM_WORLD, m_ndim, nx, &periodic[0], 0, &cart_comm);

  /* find my rank in the Cartesian communicator */
  MPI_Comm_rank(cart_comm, &myrank);

  MPI_Cart_coords(cart_comm, myrank, m_ndim, ix);

  LocalNcell = 1;
  ngbprocs   = vector<int>(2 * m_ndim);
  for (int i = 0; i < m_ndim; i++) {
    LocalRange[i] = level.Range[i] / nx[i];
    LocalXmin[i]  = level.Xmin[i] + ix[i] * LocalRange[i];
    LocalXmax[i]  = level.Xmin[i] + (ix[i] + 1) * LocalRange[i];
    LocalNG[i]    = level.NG[i] / nx[i];
    offsets[i]    = LocalNG[i] * ix[i];
    LocalNcell *= LocalNG[i];

    /* determine neighbours in ith dimension */
    MPI_Cart_shift(cart_comm, i, 1, &ngbprocs[2 * i], &ngbprocs[2 * i + 1]);
  }

#ifndef NDEBUG
  if (myrank == 0) {
    cout << "Global Ncell = " << level.Ncell << "\n";
    cout << "Dim\t Range\t Xmin\t Xmax\t Grid size\n";
    for (int i = 0; i < m_ndim; i++) {
      cout << i << "\t " << level.Range[i] << "\t " << level.Xmin[i] << "\t "
           << level.Xmax[i] << "\t " << nx[i] << "\n";
    }
    std::vector<int> coords(m_ndim, 0);
    cout << "Process topology:\n";
    print_grid(coords, 0);
  }

  cout << "Proc " << myrank << ":\n";
  cout << "Ncell = " << LocalNcell << "\n";
  cout << "Dim\t Range\t Xmin\t Xmax\t -ngbr\t +ngbr\t coords\n";
  for (int i = 0; i < m_ndim; i++) {
    cout << i << "\t " << LocalRange[i] << "\t " << LocalXmin[i] << "\t "
         << LocalXmax[i] << "\t " << ngbprocs[2 * i] << "\t "
         << ngbprocs[2 * i + 1] << "\t " << ix[i] << "\n";
  }
  cout << "---MCMDcontrol::decomposeDomain() Domain decomposition done.\n"
       << endl;
#endif
  return 0;
}

// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
    const enum axes daxis,  ///< Axis to decompose domain along.
    const int &ndim,        ///< pointer to simulation parameters
    const class level
        &level  ///< pointer to domain parameters for NG grid level
)
{
  for (int i = 0; i < ndim; i++) {
    /* don't decompose domain by default */
    nx[i] = 1;
  }
  /* do decompose in daxis dimension */
  nx[daxis] = 0;

  return decomposeDomain(ndim, level);
}

// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
    const int &ndim,           ///< number of dimensions
    const class level &level,  ///< parameters for NG grid level
    vector<int> &&pbc          ///< boolean array of whether each face has pbc
)
{
  periodic = move(pbc);
#ifndef NDEBUG
  rep.printVec("periodic", &periodic[0], ndim);
#endif
  return decomposeDomain(ndim, level);
}

// ##################################################################
// ##################################################################


int MCMDcontrol::decomposeDomain(
    const enum axes daxis,  ///< Axis to decompose domain along.
    const int &ndim,        ///< pointer to simulation parameters
    const class level
        &level,        ///< pointer to domain parameters for NG grid level
    vector<int> &&pbc  ///< boolean array of whether each face has pbc
)
{
  periodic = move(pbc);
#ifndef NDEBUG
  rep.printVec("periodic", &periodic[0], ndim);
#endif
  return decomposeDomain(daxis, ndim, level);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::create_abutting_domains_list(
    int current_dimension,  ///< current dimension to traverse
    int cursor[],           ///< position to traverse relative to
    bool explore)
{
  if (current_dimension != 0) {
    create_abutting_domains_list(current_dimension - 1, cursor, explore);
  }

  if (!explore) {
    /* have neighbours in the negative direction */
    if (--cursor[current_dimension] > -1) {
      int neighbour_rank;
      MPI_Cart_rank(cart_comm, cursor, &neighbour_rank);
      abutting_domains.push_back(neighbour_rank);
    }
    /* have neighbours in the positive direction */
    if ((cursor[current_dimension] += 2) < nx[current_dimension]) {
      int neighbour_rank;
      MPI_Cart_rank(cart_comm, cursor, &neighbour_rank);
      abutting_domains.push_back(neighbour_rank);
    }
    cursor[current_dimension]--;
  }
  else if (current_dimension != 0) {
    /* have neighbours in the negative direction */
    if (--cursor[current_dimension] > -1) {
      create_abutting_domains_list(current_dimension - 1, cursor, false);
      create_abutting_domains_list(current_dimension - 1, cursor, true);
    }
    /* have neighbours in the positive direction */
    if ((cursor[current_dimension] += 2) < nx[current_dimension]) {
      create_abutting_domains_list(current_dimension - 1, cursor, false);
      create_abutting_domains_list(current_dimension - 1, cursor, true);
    }
    cursor[current_dimension]--;
  }
}

// ##################################################################
// ##################################################################


void MCMDcontrol::create_abutting_domains_list()
{
  int cursor[m_ndim];
  for (int i = 0; i < m_ndim; i++) {
    cursor[i] = ix[i];
  }
  if (abutting_domains.empty()) {
    create_abutting_domains_list(m_ndim - 1, cursor, false);
    create_abutting_domains_list(m_ndim - 1, cursor, true);
  }
}

// ##################################################################
// ##################################################################


void MCMDcontrol::get_domain_coordinates(const int r, int *arr) const
{
  MPI_Cart_coords(cart_comm, r, m_ndim, arr);
}

// ##################################################################
// ##################################################################


int MCMDcontrol::get_rank_from_grid_location(
    const class SimParams &par,  ///< simulation parameters
    const vector<double> &loc,   ///< location sought
    const int l                  ///< grid level
    ) const
{
  int grid_coordinates[par.ndim];
  for (int i = 0; i < par.ndim; i++) {
    grid_coordinates[i] = static_cast<int>(floor(
        nx[i] * (loc[i] - par.levels[l].Xmin[i]) / par.levels[l].Range[i]));
  }
  /* check validity (0 <= x <= nx[] - 1) */
  for (int i = 0; i < par.ndim; i++) {
    if ((grid_coordinates[i] * (grid_coordinates[i] + 1 - nx[i]) > 0)) {
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


void MCMDcontrol::determine_parent_processes(
    const class SimParams &par, const int level)
{
  vector<double> centre(par.ndim);
  for (int i = 0; i < par.ndim; i++)
    centre[i] = 0.5 * (LocalXmin[i] + LocalXmax[i]);

  /* pgrid contains information about the parent process' domain */
  pgrid.rank = get_rank_from_grid_location(par, centre, level - 1);
  /* parent grid's neighbouring ranks */
  pgrid_ngb.resize(2 * par.ndim);

  int domain_coordinates[par.ndim];
  get_domain_coordinates(pgrid.rank, domain_coordinates);
  for (int i = 0; i < par.ndim; i++) {
    /* parent level has range 2x my range in each direction */
    pgrid.Xmin[i] = par.levels[level - 1].Xmin[i]
                    + domain_coordinates[i] * 2.0 * LocalRange[i];
    pgrid.Xmax[i] = pgrid.Xmin[i] + 2.0 * LocalRange[i];

    /* parent's neighbours in this dimension */
    /* negative x-dir */
    if (--domain_coordinates[i] < 0)
      pgrid_ngb[2 * i].rank = -999;
    else {
      MPI_Cart_rank(cart_comm, domain_coordinates, &(pgrid_ngb[2 * i].rank));
      for (int j = 0; j < par.ndim; j++) {
        pgrid_ngb[2 * i].Xmin[j] =
            par.levels[level - 1].Xmin[j]
            + domain_coordinates[j] * 2.0 * LocalRange[j];
        pgrid_ngb[2 * i].Xmax[j] =
            pgrid_ngb[2 * i].Xmin[j] + 2.0 * LocalRange[j];
      }
    }
    /* positive x-dir */
    if ((domain_coordinates[i] += 2) == nx[i])
      pgrid_ngb[2 * i + 1].rank = -999;
    else {
      MPI_Cart_rank(
          cart_comm, domain_coordinates, &(pgrid_ngb[2 * i + 1].rank));
      for (int j = 0; j < par.ndim; j++) {
        pgrid_ngb[2 * i + 1].Xmin[j] =
            par.levels[level - 1].Xmin[j]
            + domain_coordinates[j] * 2.0 * LocalRange[j];
        pgrid_ngb[2 * i + 1].Xmax[j] =
            pgrid_ngb[2 * i + 1].Xmin[j] + 2.0 * LocalRange[j];
      }
    }
    domain_coordinates[i] -= 1;
  }

#ifndef NDEBUG
  cout << "level " << level << ", parent proc = " << pgrid.rank
       << ", parent grid xmin/xmax:\n";
  rep.printVec("Xmin", pgrid.Xmin, par.ndim);
  rep.printVec("Xmax", pgrid.Xmax, par.ndim);
  for (int d = 0; d < 2 * par.ndim; d++) {
    if (pgrid_ngb[d].rank > -1) {
      cout << "pproc ngb in direction " << d << " has neighbour proc ";
      cout << pgrid_ngb[d].rank << "\n";
      rep.printVec("Xmin", pgrid_ngb[d].Xmin, par.ndim);
      rep.printVec("Xmax", pgrid_ngb[d].Xmax, par.ndim);
    }
  }
  cout.flush();
#endif
}

// ##################################################################
// ##################################################################


void MCMDcontrol::determine_child_ranks(
    const class SimParams &par,
    const int level,
    vector<int> &children,
    vector<double> &centre,
    const int current_dimension) const
{
  /* split domain in 4 in each dimension */
  for (int i = 0; i < 4; i++) {
    centre[current_dimension] =
        LocalXmin[current_dimension] + 0.25 * i * LocalRange[current_dimension];
    if (current_dimension == m_ndim - 1) {
      int child_rank = get_rank_from_grid_location(par, centre, level + 1);
      if (child_rank > -1) children.push_back(child_rank);
    }
    else {
      determine_child_ranks(
          par, level, children, centre, current_dimension + 1);
    }
  }
}

// ##################################################################
// ##################################################################


void MCMDcontrol::determine_child_processes(
    const class SimParams &par, const int level)
{
  child_procs.clear();
  if (nproc == 1) /* special case */ {
    struct cgrid child;
    child.rank = myrank;
    for (int v = 0; v < par.ndim; v++) {
      child.Xmin[v] = par.levels[level + 1].Xmin[v];
      child.Xmax[v] = par.levels[level + 1].Xmax[v];
    }
    child_procs.push_back(child);
    return;
  }

  vector<int> child_ranks;
  vector<double> centre(par.ndim);
  determine_child_ranks(par, level, child_ranks, centre, 0);

  if (!child_ranks.empty()) {
    sort(child_ranks.begin(), child_ranks.end());
    child_ranks.erase(
        unique(child_ranks.begin(), child_ranks.end()), child_ranks.end());

    std::vector<struct cgrid> children;
    for (auto c : child_ranks) {
      struct cgrid child;
      child.rank = c;
      int domain_coordinates[par.ndim];
      get_domain_coordinates(child.rank, domain_coordinates);
      for (int i = 0; i < par.ndim; i++) {
        child.Xmin[i] = par.levels[level + 1].Xmin[i]
                        + domain_coordinates[i] * 0.5 * LocalRange[i];
        child.Xmax[i] = child.Xmin[i] + 0.5 * LocalRange[i];
      }
      child_procs.push_back(child);
    }
  }
#ifndef NDEBUG
  cout << "child procs on level " << level + 1 << "\n";
  for (auto c : child_procs) {
    cout << "child rank " << c.rank << ":\n";
    rep.printVec("Xmin", c.Xmin, par.ndim);
    rep.printVec("Xmax", c.Xmax, par.ndim);
    cout.flush();
  }
#endif
}

// ##################################################################
// ##################################################################


void MCMDcontrol::determine_child_neighbour_ranks(
    const class SimParams &par,
    const int level,
    vector<double> &centre,
    vector<vector<int> > &neighbours,
    const int current_dimension,
    const int neighbour_dimension) const
{
  const int next_dimension = (current_dimension + 1) % par.ndim;
  for (int i = 0; i < 4; i++) {
    centre[current_dimension] =
        LocalXmin[current_dimension] + 0.25 * i * LocalRange[current_dimension];
    if (next_dimension != neighbour_dimension) {
      determine_child_neighbour_ranks(
          par, level, centre, neighbours, next_dimension, neighbour_dimension);
    }
    else {
      /* negative direction */
      centre[neighbour_dimension] =
          LocalXmin[neighbour_dimension]
          - LocalRange[neighbour_dimension] / LocalNG[neighbour_dimension];
      int child_rank = get_rank_from_grid_location(par, centre, level + 1);
      if (child_rank > -1) neighbours[0].push_back(child_rank);

      /* positive direction */
      centre[neighbour_dimension] =
          LocalXmax[neighbour_dimension]
          + LocalRange[neighbour_dimension] / LocalNG[neighbour_dimension];
      child_rank = get_rank_from_grid_location(par, centre, level + 1);
      if (child_rank > -1) neighbours[1].push_back(child_rank);
    }
  }
}

// ##################################################################
// ##################################################################


void MCMDcontrol::determine_child_neighbours(
    const class SimParams &par, const int level)
{
  cgrid_ngb.resize(2 * par.ndim);
  for (auto c : cgrid_ngb) {
    c.clear();
  }

  /* ranks of neighbours in each direction */
  for (int i = 0; i < par.ndim; i++) {
    vector<double> centre(par.ndim);
    vector<vector<int> > neighbour_ranks(2);
    determine_child_neighbour_ranks(
        par, level, centre, neighbour_ranks, (i + 1) % par.ndim, i);

    for (int j = 0; j < 2; j++) {
      if (!neighbour_ranks[j].empty()) {
        sort(neighbour_ranks[j].begin(), neighbour_ranks[j].end());
        neighbour_ranks[j].erase(
            unique(neighbour_ranks[j].begin(), neighbour_ranks[j].end()),
            neighbour_ranks[j].end());

        for (auto c : neighbour_ranks[j]) {
          struct cgrid child;
          child.rank = c;
          int domain_coordinates[par.ndim];
          get_domain_coordinates(child.rank, domain_coordinates);
          for (int i = 0; i < par.ndim; i++) {
            child.Xmin[i] = par.levels[level + 1].Xmin[i]
                            + domain_coordinates[i] * 0.5 * LocalRange[i];
            child.Xmax[i] = child.Xmin[i] + 0.5 * LocalRange[i];
          }
          cgrid_ngb[2 * i + j].push_back(child);
        }
      }
    }
  }
#ifndef NDEBUG
  cout << "Child process' neighbours:\n";
  for (int d = 0; d < 2 * par.ndim; d++) {
    cout << "dir " << d << ":\n";
    for (auto ngb : cgrid_ngb[d]) {
      cout << "rank = " << ngb.rank << "\n";
      rep.printVec("Xmin", ngb.Xmin, par.ndim);
      rep.printVec("Xmax", ngb.Xmax, par.ndim);
      cout.flush();
    }
  }
#endif
}

// ##################################################################
// ##################################################################


void MCMDcontrol::set_NG_hierarchy(
    class SimParams &par,  ///< simulation parameters
    const int l            ///< level to work on
)
{
#ifndef NDEBUG
  cout << "Setting up NG hierarchy (MCMDcontrol) on level " << l << endl;
  ;
#endif
  // set rank of parent for each grid except root level 0
  // Every child grid must have a parent because the nested grid is
  // entirely within the coarser level.
  if (l > 0) {
    determine_parent_processes(par, l);
  }

  // set rank of child grids, if they exist.
  // Domain is intersection of child full grid and this local grid.
  // Must be split in half/quadrant/octant of this grid, unless
  // nproc==1, for which there is only one child so it is trivial.
  if (l < par.grid_nlevels - 1) {
    determine_child_processes(par, l);
    if (nproc != 1) {
      determine_child_neighbours(par, l);
    }
  }
}


// ##################################################################
// ##################################################################


void MCMDcontrol::get_nx_subdomains(int *tmp) const
{
  for (int i = 0; i < MAX_DIM; i++)
    tmp[i] = nx[i];
  return;
}


// ##################################################################
// ##################################################################


void MCMDcontrol::get_parent_grid_info(struct cgrid *cg) const
{
  cg->rank = pgrid.rank;
  for (int i = 0; i < MAX_DIM; i++)
    cg->Xmin[i] = pgrid.Xmin[i];
  for (int i = 0; i < MAX_DIM; i++)
    cg->Xmax[i] = pgrid.Xmax[i];
}

// ##################################################################
// ##################################################################


void MCMDcontrol::get_parent_ngb_grid_info(vector<struct cgrid> &pgngb) const
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


void MCMDcontrol::get_child_grid_info(vector<struct cgrid> &cg) const
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


void MCMDcontrol::get_level_lp1_ngb_info(
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


void MCMDcontrol::gather_ncells(int *recv_buffer, const int &root) const
{
  MPI_Gather(&LocalNcell, 1, MPI_INT, recv_buffer, 1, MPI_INT, root, cart_comm);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::allgather_ncells(vector<int> &ncells_list) const
{
  ncells_list.resize(nproc);
  MPI_Allgather(
      &LocalNcell, 1, MPI_INT, &ncells_list[0], 1, MPI_INT, cart_comm);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::gather_extents(double *recv_buffer, const int &root) const
{
  double extents[2 * m_ndim];
  for (int i = 0; i < m_ndim; i++) {
    extents[i]          = LocalXmin[i];
    extents[m_ndim + i] = LocalXmax[i];
  }
  MPI_Gather(
      &extents, 2 * m_ndim, MPI_DOUBLE, recv_buffer, 2 * m_ndim, MPI_DOUBLE,
      root, cart_comm);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::gather_offsets(
    vector<int> &offsets_list, const int &root) const
{
  offsets_list.resize(nproc * m_ndim);
  MPI_Gather(
      &offsets[0], m_ndim, MPI_INT, &offsets_list[0], m_ndim, MPI_INT, root,
      cart_comm);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::gather_localNG(
    vector<int> &localNG_list, const int &root) const
{
  localNG_list.resize(nproc * m_ndim);
  MPI_Gather(
      &LocalNG[0], m_ndim, MPI_INT, &localNG_list[0], m_ndim, MPI_INT, root,
      cart_comm);
}

// ##################################################################
// ##################################################################


void MCMDcontrol::gather_abutting_domains(
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

//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------
