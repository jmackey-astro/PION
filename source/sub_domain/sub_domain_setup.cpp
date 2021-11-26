/// \file sub_domain.cpp
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

#include <algorithm>
#include <tools/command_line_interface.h>
#include <tools/mem_manage.h>

#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <spdlog/fmt/bundled/ranges.h>

#include "sub_domain.h"

using namespace std;

//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------


unsigned int Sub_domain::m_count;

// ##################################################################
// ##################################################################


Sub_domain::Sub_domain()
{
  Ncell = pgrid.rank = -1;
  periodic           = vector<int>(MAX_DIM, 0);
  for (int i = 0; i < MAX_DIM; i++) {
    directional_Ncells[i] = offsets[i] = coordinates[i] = -1;
    num_subdomains[i]                                   = 0;
    Xmin[i] = Xmax[i] = range[i] = -1.e99;
  }

  ReadSingleFile =
      true;  ///< If the ICs are in a single file, set this to true.
  WriteSingleFile = false;  ///< If all processes to write to one file, set this
  WriteFullImage =
      false;  ///< If multiple fits files, each is the full domain size.

  /* multiple Sub_domains may be initialised by each process, we need to ensure
   * only one of these takes responsibility for initialising and finialising MPI
   */
  if (m_count == 0) {
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  ++m_count;

  NG_F2C_send_list.clear();
  NG_C2F_send_list.clear();
  NG_C2F_send_list.resize(0);
  BC89_flux_send_list.clear();
}

// ##################################################################
// ##################################################################


Sub_domain::~Sub_domain()
{
  if (cart_comm) {
    MPI_Comm_free(&cart_comm);
  }
  if (m_count == 1) {
    MPI_Finalize();
  }

  --m_count;

#ifndef NDEBUG
  spdlog::info("Sub_domain Destructor: done");
#endif
}

// ##################################################################
// ##################################################################


void Sub_domain::calculate_process_topology(vector<float> &&ratios)
{
  /* account for hardcoded decomposition (i.e, num_subdomains[i] already
   * non-zero) */
  int free_processes = nproc;
  for (int i = 0; i < m_ndim; i++) {
    if (num_subdomains[i] > 0) {
      free_processes /= num_subdomains[i];
      ratios[i] = 0;
    }
    else {
      num_subdomains[i] = 1;
    }
  }

  for (; free_processes > 1; free_processes /= 2) {
    const int max_index =
        max_element(ratios.begin(), ratios.end()) - ratios.begin();
    ratios[max_index] /= 2;
    num_subdomains[max_index] *= 2;
  }
}

// ##################################################################
// ##################################################################


int Sub_domain::decomposeDomain(
    const int &ndim,          ///< number of dimensions
    const class level &level  ///< parameters for NG grid level
)
{
#ifndef NDEBUG
  spdlog::info(
      "---Sub_domain::decomposeDomain() decomposing domain.  Nproc={}, myrank={}",
      nproc, myrank);
#endif

  m_ndim = ndim;

  /*
   * the following sets num_subdomains to application topology aware values
   * dimensions of num_subdomains not set to 0 are restricted to the previously
   * held value
   */
  vector<float> ratios(begin(level.Range), begin(level.Range) + m_ndim);
  calculate_process_topology(move(ratios));

  /* TODO: reordering of ranks is temporarily disabled so old communicator works
   */
  MPI_Cart_create(
      MPI_COMM_WORLD, m_ndim, num_subdomains.data(), &periodic[0], 0,
      &cart_comm);

  /* find my rank in the Cartesian communicator */
  MPI_Comm_rank(cart_comm, &myrank);

  MPI_Cart_coords(cart_comm, myrank, m_ndim, coordinates.data());

  Ncell           = 1;
  neighbour_ranks = vector<int>(2 * m_ndim);

  for (int i = 0; i < m_ndim; i++) {
    range[i]              = level.Range[i] / num_subdomains[i];
    Xmin[i]               = level.Xmin[i] + coordinates[i] * range[i];
    Xmax[i]               = level.Xmin[i] + (coordinates[i] + 1) * range[i];
    directional_Ncells[i] = level.NG[i] / num_subdomains[i];
    offsets[i]            = directional_Ncells[i] * coordinates[i];
    Ncell *= directional_Ncells[i];

    /* determine neighbours in ith dimension */
    MPI_Cart_shift(
        cart_comm, i, 1, &neighbour_ranks[2 * i], &neighbour_ranks[2 * i + 1]);
  }

#ifndef NDEBUG
  if (myrank == 0) {
    spdlog::debug("Global Ncell = {}", level.Ncell);
    spdlog::debug("Dim\t Range\t Xmin\t Xmax\t Grid size");
    for (int i = 0; i < m_ndim; i++) {
      spdlog::debug(
          "{:^3}\t {:^5}\t {:^4}\t {:^4}\t {:^9}", i, level.Range[i],
          level.Xmin[i], level.Xmax[i], num_subdomains[i]);
    }
    std::vector<int> coords(m_ndim, 0);
    spdlog::debug("Process topology:");
    print_grid(coords, 0);
  }

  spdlog::debug("Proc {}:\tNcell = {}", myrank, Ncell);
  spdlog::debug("Dim\t Range\t Xmin\t Xmax\t -ngbr\t +ngbr\t coords");
  for (int i = 0; i < m_ndim; i++) {
    spdlog::debug(
        "{:^3}\t {:^5}\t {:^4}\t {:^4}\t {:^5}\t {:^5}\t {:^6}", i, range[i],
        Xmin[i], Xmax[i], neighbour_ranks[2 * i], neighbour_ranks[2 * i + 1],
        coordinates[i]);
  }
  spdlog::info("---Sub_domain::decomposeDomain() Domain decomposition done");
#endif
  return 0;
}

// ##################################################################
// ##################################################################


int Sub_domain::decomposeDomain(
    const enum axes daxis,  ///< Axis to decompose domain along.
    const int &ndim,        ///< pointer to simulation parameters
    const class level
        &level  ///< pointer to domain parameters for NG grid level
)
{
  for (int i = 0; i < ndim; i++) {
    /* don't decompose domain by default */
    num_subdomains[i] = 1;
  }
  /* do decompose in daxis dimension */
  num_subdomains[daxis] = 0;

  return decomposeDomain(ndim, level);
}

// ##################################################################
// ##################################################################


int Sub_domain::decomposeDomain(
    const int &ndim,           ///< number of dimensions
    const class level &level,  ///< parameters for NG grid level
    vector<int> &&pbc          ///< boolean array of whether each face has pbc
)
{
  periodic = move(pbc);
#ifndef NDEBUG
  spdlog::debug("periodic : {}", periodic);
#endif
  return decomposeDomain(ndim, level);
}

// ##################################################################
// ##################################################################


int Sub_domain::decomposeDomain(
    const enum axes daxis,  ///< Axis to decompose domain along.
    const int &ndim,        ///< pointer to simulation parameters
    const class level
        &level,        ///< pointer to domain parameters for NG grid level
    vector<int> &&pbc  ///< boolean array of whether each face has pbc
)
{
  periodic = move(pbc);
#ifndef NDEBUG
  spdlog::debug("periodic : {}", periodic);
#endif
  return decomposeDomain(daxis, ndim, level);
}

// ##################################################################
// ##################################################################


void Sub_domain::create_abutting_domains_list(
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
    if ((cursor[current_dimension] += 2) < num_subdomains[current_dimension]) {
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
    if ((cursor[current_dimension] += 2) < num_subdomains[current_dimension]) {
      create_abutting_domains_list(current_dimension - 1, cursor, false);
      create_abutting_domains_list(current_dimension - 1, cursor, true);
    }
    cursor[current_dimension]--;
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::create_abutting_domains_list()
{
  int cursor[m_ndim];
  for (int i = 0; i < m_ndim; i++) {
    cursor[i] = coordinates[i];
  }
  if (abutting_domains.empty()) {
    create_abutting_domains_list(m_ndim - 1, cursor, false);
    create_abutting_domains_list(m_ndim - 1, cursor, true);
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::determine_parent_processes(
    const class SimParams &par, const int level)
{
  vector<double> centre(par.ndim);
  for (int i = 0; i < par.ndim; i++)
    centre[i] = 0.5 * (Xmin[i] + Xmax[i]);

  /* pgrid contains information about the parent process' domain */
  pgrid.rank = get_rank_from_grid_location(par, centre, level - 1);
  /* parent grid's neighbouring ranks */
  pgrid_ngb.resize(2 * par.ndim);

  int domain_coordinates[par.ndim];
  get_domain_coordinates(pgrid.rank, domain_coordinates);
  for (int i = 0; i < par.ndim; i++) {
    /* parent level has range 2x my range in each direction */
    pgrid.Xmin[i] =
        par.levels[level - 1].Xmin[i] + domain_coordinates[i] * 2.0 * range[i];
    pgrid.Xmax[i] = pgrid.Xmin[i] + 2.0 * range[i];

    /* parent's neighbours in this dimension */
    /* negative x-dir */
    if (--domain_coordinates[i] < 0)
      pgrid_ngb[2 * i].rank = -999;
    else {
      MPI_Cart_rank(cart_comm, domain_coordinates, &(pgrid_ngb[2 * i].rank));
      for (int j = 0; j < par.ndim; j++) {
        pgrid_ngb[2 * i].Xmin[j] = par.levels[level - 1].Xmin[j]
                                   + domain_coordinates[j] * 2.0 * range[j];
        pgrid_ngb[2 * i].Xmax[j] = pgrid_ngb[2 * i].Xmin[j] + 2.0 * range[j];
      }
    }
    /* positive x-dir */
    if ((domain_coordinates[i] += 2) == num_subdomains[i])
      pgrid_ngb[2 * i + 1].rank = -999;
    else {
      MPI_Cart_rank(
          cart_comm, domain_coordinates, &(pgrid_ngb[2 * i + 1].rank));
      for (int j = 0; j < par.ndim; j++) {
        pgrid_ngb[2 * i + 1].Xmin[j] = par.levels[level - 1].Xmin[j]
                                       + domain_coordinates[j] * 2.0 * range[j];
        pgrid_ngb[2 * i + 1].Xmax[j] =
            pgrid_ngb[2 * i + 1].Xmin[j] + 2.0 * range[j];
      }
    }
    domain_coordinates[i] -= 1;
  }

#ifndef NDEBUG
  spdlog::debug(
      "level {}, parent proc = {}, parent grid xmin/xmax:", level, pgrid.rank);
  spdlog::debug("Xmin : {}", pgrid.Xmin);
  spdlog::debug("Xmax : {}", pgrid.Xmax);
  for (int d = 0; d < 2 * par.ndim; d++) {
    if (pgrid_ngb[d].rank > -1) {
      spdlog::debug("pproc ngb in direction {} has neighbour proc ", d);
      spdlog::debug("{}", pgrid_ngb[d].rank);
      spdlog::debug("Xmin : {}", pgrid_ngb[d].Xmin);
      spdlog::debug("Xmax : {}", pgrid_ngb[d].Xmax);
    }
  }
#endif
}

// ##################################################################
// ##################################################################


void Sub_domain::determine_child_ranks(
    const class SimParams &par,
    const int level,
    vector<int> &children,
    vector<double> &centre,
    const int current_dimension) const
{
  /* split domain in 4 in each dimension */
  for (int i = 0; i < 4; i++) {
    centre[current_dimension] =
        Xmin[current_dimension] + 0.25 * i * range[current_dimension];
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


void Sub_domain::determine_child_processes(
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
                        + domain_coordinates[i] * 0.5 * range[i];
        child.Xmax[i] = child.Xmin[i] + 0.5 * range[i];
      }
      child_procs.push_back(child);
    }
  }
#ifndef NDEBUG
  spdlog::debug("child procs on level {}", level + 1);
  for (auto c : child_procs) {
    spdlog::debug("child rank {}", c.rank);
    spdlog::debug("Xmin : {}", c.Xmin);
    spdlog::debug("Xmax : {}", c.Xmax);
  }
#endif
}

// ##################################################################
// ##################################################################


void Sub_domain::determine_child_neighbour_ranks(
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
        Xmin[current_dimension] + 0.25 * i * range[current_dimension];
    if (next_dimension != neighbour_dimension) {
      determine_child_neighbour_ranks(
          par, level, centre, neighbours, next_dimension, neighbour_dimension);
    }
    else {
      /* negative direction */
      centre[neighbour_dimension] =
          Xmin[neighbour_dimension]
          - range[neighbour_dimension]
                / directional_Ncells[neighbour_dimension];
      int child_rank = get_rank_from_grid_location(par, centre, level + 1);
      if (child_rank > -1) neighbours[0].push_back(child_rank);

      /* positive direction */
      centre[neighbour_dimension] =
          Xmax[neighbour_dimension]
          + range[neighbour_dimension]
                / directional_Ncells[neighbour_dimension];
      child_rank = get_rank_from_grid_location(par, centre, level + 1);
      if (child_rank > -1) neighbours[1].push_back(child_rank);
    }
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::determine_child_neighbours(
    const class SimParams &par, const int level)
{
  cgrid_ngb.resize(2 * par.ndim);
  for (auto c : cgrid_ngb) {
    c.clear();
  }
  if (par.ndim == 1) return;

  /* ranks of neighbours in each direction */
  for (int i = 0; i < par.ndim; ++i) {
    vector<double> centre(par.ndim);
    vector<vector<int> > neighbours(2);
    determine_child_neighbour_ranks(
        par, level, centre, neighbours, (i + 1) % par.ndim, i);

    for (int j = 0; j < 2; ++j) {
      if (!neighbours[j].empty()) {
        sort(neighbours[j].begin(), neighbours[j].end());
        neighbours[j].erase(
            unique(neighbours[j].begin(), neighbours[j].end()),
            neighbours[j].end());

        for (auto c : neighbours[j]) {
          struct cgrid child;
          child.rank = c;
          int domain_coordinates[par.ndim];
          get_domain_coordinates(child.rank, domain_coordinates);
          for (int k = 0; k < par.ndim; ++k) {
            child.Xmin[k] = par.levels[level + 1].Xmin[k]
                            + domain_coordinates[k] * 0.5 * range[k];
            child.Xmax[k] = child.Xmin[k] + 0.5 * range[k];
          }
          cgrid_ngb[2 * i + j].push_back(child);
        }
      }
    }
  }
#ifndef NDEBUG
  spdlog::debug("Child process' neighbours:");
  for (int d = 0; d < 2 * par.ndim; d++) {
    spdlog::debug("dir {}:", d);
    for (auto ngb : cgrid_ngb[d]) {
      spdlog::debug("rank = {}", ngb.rank);
      spdlog::debug("Xmin : {}", ngb.Xmin);
      spdlog::debug("Xmax : {}", ngb.Xmax);
    }
  }
#endif
}

// ##################################################################
// ##################################################################


void Sub_domain::set_NG_hierarchy(
    class SimParams &par,  ///< simulation parameters
    const int l            ///< level to work on
)
{
#ifndef NDEBUG
  spdlog::info("Setting up NG hierarchy (Sub_domain) on level {}", l);
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
    determine_child_neighbours(par, l);
  }
}


//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------
