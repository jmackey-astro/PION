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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_params.h"
#include "sub_domain.h"
#include <tools/mem_manage.h>

#include <algorithm>

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

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
  if (m_count++ == 0) {
#ifndef NDEBUG
    spdlog::debug("Subdomain::Subdomain() initialising MPI");
#endif
    MPI_Init(nullptr, nullptr);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  NG_F2C_send_list.clear();
  NG_C2F_send_list.clear();
  NG_C2F_send_list.resize(0);
  BC89_flux_send_list.clear();
}

// ##################################################################
// ##################################################################


Sub_domain::~Sub_domain()
{
  NG_F2C_send_list.clear();
  NG_C2F_send_list.clear();
  NG_C2F_send_list.resize(0);
  BC89_flux_send_list.clear();
  send_list.clear();
  recv_list.clear();

  if (--m_count == 0) {
#ifndef NDEBUG
    spdlog::debug("Subdomain::Subdomain() finalising MPI");
#endif
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
  }

#ifndef NDEBUG
  spdlog::debug("Sub_domain Destructor: done");
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
  spdlog::info("num_subdomains: {}", num_subdomains);
#endif

  Ncell           = 1;
  neighbour_ranks = vector<int>(2 * m_ndim, -1);

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
  spdlog::debug("---Sub_domain::decomposeDomain() Domain decomposition done");
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
  num_subdomains[daxis] = nproc;

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
  periodic = std::move(pbc);
#ifndef NDEBUG
  spdlog::debug("periodic : {}", periodic);
#endif
  m_ndim = ndim;

  /*
   * the following sets num_subdomains to application topology aware values
   * dimensions of num_subdomains not set to 0 are restricted to the previously
   * held value
   */
  vector<float> ratios(begin(level.Range), begin(level.Range) + m_ndim);
  calculate_process_topology(std::move(ratios));

  /* TODO: reordering of ranks is temporarily disabled so old communicator works
   */
  MPI_Cart_create(
      MPI_COMM_WORLD, m_ndim, num_subdomains.data(), &periodic[0], 0,
      &cart_comm);

  m_is_cartcomm = true;

  /* find my rank in the Cartesian communicator */
  MPI_Comm_rank(cart_comm, &myrank);

  MPI_Cart_coords(cart_comm, myrank, m_ndim, coordinates.data());

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
  periodic = std::move(pbc);
#ifndef NDEBUG
  spdlog::debug("periodic : {}", periodic);
#endif
  m_ndim = ndim;

  for (int i = 0; i < ndim; i++) {
    /* don't decompose domain by default */
    num_subdomains[i] = 1;
  }
  /* do decompose in daxis dimension */
  num_subdomains[daxis] = nproc;

  /*
   * the following sets num_subdomains to application topology aware values
   * dimensions of num_subdomains not set to 0 are restricted to the previously
   * held value
   */
  vector<float> ratios(begin(level.Range), begin(level.Range) + m_ndim);
  calculate_process_topology(std::move(ratios));

  /* TODO: reordering of ranks is temporarily disabled so old communicator works
   */
  MPI_Cart_create(
      MPI_COMM_WORLD, m_ndim, num_subdomains.data(), &periodic[0], 0,
      &cart_comm);

  m_is_cartcomm = true;

  /* find my rank in the Cartesian communicator */
  MPI_Comm_rank(cart_comm, &myrank);

  MPI_Cart_coords(cart_comm, myrank, m_ndim, coordinates.data());

  return decomposeDomain(daxis, ndim, level);
}

// ##################################################################
// ##################################################################


/* The abutting domains vector is required by Silo to be in the following order;
 *      XN,XP,YN,YP,ZN,ZP
 *      YNXN, YNXP, YPXN, YPXP,
 *      ZNXN, ZNXP, ZNYN, ZNYP,
 *      ZNYNXN, ZNYNXP, ZNYPXN, ZNYPXP
 *      ZPXN, ZPXP, ZPYN, ZPYP,
 *      ZPYNXN, ZPYNXP, ZPYPXN, ZPYPXP
 *
 * These can be found in the correct order by first adding (in order of
 * dimension) each basis vector which corresponds to a valid neighbour
 * (XN, XP, ..., ZP). Then, in a queue-like fashion, checking for combinations
 * of a basis vector (eg. YN) with known valid neighbours in the queue (eg. XN,
 * XP).
 */
void Sub_domain::create_abutting_domains_list()
{
  if (!abutting_domains.empty()) return;

  int max_neighbours = pow(3, m_ndim) - 1;

  vector<vector<int> > neighbours;
  neighbours.reserve(max_neighbours);
  /* add all single direction neighbours, all other abutting domains are sums
   * of these bases */
  vector<int> neighbour(coordinates.begin(), coordinates.end());
  for (int i = 0; i < m_ndim; i++) {
    /* add negative direction neighbour in ith dimension */
    if (coordinates[i] > 0) {
      --neighbour[i];
      neighbours.push_back(neighbour);
      ++neighbour[i];
    }
    /* add positive direction neighbour in ith dimension */
    if (coordinates[i] < num_subdomains[i] - 1) {
      ++neighbour[i];
      neighbours.push_back(neighbour);
      --neighbour[i];
    }
  }

  /* each dimension checks if it can create a new valid neighbour given the
   * valid neighbours formed from lower-order dimensions */
  for (int i = 1; i < m_ndim; i++) {
    /* loop over those neighbours which have only traversed lower-order
     * dimensions. e.g for i = 1, we combine YN with XN, XP (to make YNXN, YNXP)
     * but not YN, YP, ZN, or ZP. Valid combinations are added to the list */
    auto it = neighbours.begin() - 1;
    do {
      it = find_if(it + 1, neighbours.end(), [i, this](vector<int> &neighbour) {
        return equal(
            neighbour.begin() + i, neighbour.end(), coordinates.begin() + i);
      });
      if (it == neighbours.end()) break; /* no more valid combinations in i */
      /* it is a neighbour we can create new neighbours from by moving it in
       * negative and positive directions of the ith dimension */
      if (coordinates[i] > 0) {
        --(*it)[i];
        neighbours.push_back(*it);
        ++(*it)[i];
      }
      if (coordinates[i] < num_subdomains[i] - 1) {
        ++(*it)[i];
        neighbours.push_back(*it);
        --(*it)[i];
      }
    } while (true);
  }

  /* now that we have the coordinates of each abutting domain, we are left to
   * determine their ranks */
  abutting_domains.reserve(neighbours.size());
#ifndef NDEBUG
  spdlog::debug("abutting_domains:\n {}", neighbours);
#endif
  for (auto const &n : neighbours) {
    int neighbour_rank;
    MPI_Cart_rank(cart_comm, n.data(), &neighbour_rank);
    abutting_domains.push_back(neighbour_rank);
  }
}

// ##################################################################
// ##################################################################


void Sub_domain::evaluate_parent_process(
    const class SimParams &par, const int level)
{
  vector<double> centre(par.ndim, 0.0);
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

/*
 * recursively work through the dimensions to move cursor and try to find
 * child domains.
 * Begining at the lowest-order dimension, cursor is set to a position on the
 * domain and the function is called recursively for the next-order dimension.
 * When current_dimension is the highest-order dimension (therefore cursor has
 * a value in each dimension), it is determined whether cursor corresponds to a
 * valid location on the domain and, if so, the rank of the processor who owns
 * that location on the next grid level is added to the unordered_set, children.
 * When the recursion returns to a given dimension, the cursor is moved in that
 * dimension and the recursion is invoked again for the higher-order dimensions.
 * In this way, the recursion tree explores all possible locations for children.
 */
void Sub_domain::determine_child_ranks(
    const class SimParams &par,
    const int level,
    unordered_set<int> &children,
    vector<double> &cursor,
    const int current_dimension) const
{
  cursor[current_dimension] =
      Xmin[current_dimension] + 0.25 * range[current_dimension];

  if (current_dimension == m_ndim - 1) { /* check if cursor found child */
    cursor[current_dimension] =
        Xmin[current_dimension] + 0.25 * range[current_dimension];
    int child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) children.insert(child_rank);

    cursor[current_dimension] += 0.5 * range[current_dimension];
    child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) children.insert(child_rank);
  }
  else {
    determine_child_ranks(par, level, children, cursor, current_dimension + 1);

    cursor[current_dimension] += 0.5 * range[current_dimension];

    determine_child_ranks(par, level, children, cursor, current_dimension + 1);
  }
}

// ##################################################################
// ##################################################################

/*
 * populate the child_procs vector with cgrid structs. Each cgrid has the rank
 * of a child and the range of it's domain
 */
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

  unordered_set<int> child_ranks;
  vector<double> cursor(par.ndim, 0.0);
  determine_child_ranks(par, level, child_ranks, cursor, 0);

  if (!child_ranks.empty()) {
    child_procs.reserve(child_ranks.size());
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
  for (auto const &c : child_procs) {
    spdlog::debug("child rank {}", c.rank);
    spdlog::debug("Xmin : {}", c.Xmin);
    spdlog::debug("Xmax : {}", c.Xmax);
  }
#endif
}

// ##################################################################
// ##################################################################

/*
 * Recurse dimensions to get child process' neighbouring ranks in
 * neighbour_dimension. cursor is used to explore the space, set in each other
 * dimension first (two locations per dimension), before exploring for valid
 * neighbours in neighbour_dimension. Similar logic to calculate_child_processes
 *
 * Note: this function is called from evaluate_child_neighbours for each
 * neighbour_dimension
 */
void Sub_domain::determine_child_neighbour_ranks(
    const class SimParams &par,
    const int level,
    vector<double> &cursor,
    array<unordered_set<int>, 2> &neighbours,
    const int current_dimension,
    const int neighbour_dimension) const
{
  const int next_dimension = (current_dimension + 1) % par.ndim;
  cursor[current_dimension] =
      Xmin[current_dimension] + 0.25 * range[current_dimension];

  if (next_dimension != neighbour_dimension) {
    determine_child_neighbour_ranks(
        par, level, cursor, neighbours, next_dimension, neighbour_dimension);

    cursor[current_dimension] += 0.5 * range[current_dimension];

    determine_child_neighbour_ranks(
        par, level, cursor, neighbours, next_dimension, neighbour_dimension);
  }
  else {
    /* negative direction */
    cursor[neighbour_dimension] =
        Xmin[neighbour_dimension]
        - range[neighbour_dimension] / directional_Ncells[neighbour_dimension];
    int child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) neighbours[0].insert(child_rank);

    cursor[current_dimension] += 0.5 * range[current_dimension];

    child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) neighbours[0].insert(child_rank);

    /* positive direction */
    cursor[neighbour_dimension] =
        Xmax[neighbour_dimension]
        + range[neighbour_dimension] / directional_Ncells[neighbour_dimension];
    child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) neighbours[1].insert(child_rank);

    cursor[current_dimension] -= 0.5 * range[current_dimension];

    child_rank = get_rank_from_grid_location(par, cursor, level + 1);
    if (child_rank > -1) neighbours[1].insert(child_rank);
  }
}

// ##################################################################
// ##################################################################

/*
 * cgrid_ngb is filled with cgrid structs corresponding to the neighbours of
 * children in order such that cgrid_ngb[2 * dim + dir] contains a vector of
 * cgrids that are neighbours in dimension, dim, in direction dir (neg = 0, pos
 * = 1).
 */
void Sub_domain::evaluate_child_neighbours(
    const class SimParams &par, const int level)
{
  cgrid_ngb.resize(2 * par.ndim);
  for (auto &c : cgrid_ngb) {
    c.clear();
  }
  if (par.ndim == 1) return;

  for (int i = 0; i < par.ndim; ++i) { /* per dimension */
    vector<double> cursor(par.ndim, 0.0);
    /* ranks of neighbours in dimension i in negative and positive directions,
     * neighbours[0] is negative, neighbours[1] is positive */
    array<unordered_set<int>, 2> neighbours;
    determine_child_neighbour_ranks(
        par, level, cursor, neighbours, (i + 1) % par.ndim, i);

    for (int j = 0; j < 2; ++j) { /* per direction (negative/positive) */
      if (!neighbours[j].empty()) {
        cgrid_ngb[2 * i + j].reserve(neighbours[j].size());
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
    for (auto const &ngb : cgrid_ngb[d]) {
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
  spdlog::debug("Setting up NG hierarchy (Sub_domain) on level {}", l);
  ;
#endif
  // set rank of parent for each grid except root level 0
  // Every child grid must have a parent because the nested grid is
  // entirely within the coarser level.
  if (l > 0) {
    evaluate_parent_process(par, l);
  }

  // set rank of child grids, if they exist.
  // Domain is intersection of child full grid and this local grid.
  // Must be split in half/quadrant/octant of this grid, unless
  // nproc==1, for which there is only one child so it is trivial.
  if (l < par.grid_nlevels - 1) {
    determine_child_processes(par, l);
    evaluate_child_neighbours(par, l);
  }
}


//------------------------------------------------
//-------------- MPI PARAMETERS ------------------
//------------------------------------------------
