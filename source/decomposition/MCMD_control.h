/// \file decomposition/MCMD_control.h
///
/// \brief Declares class for controlling Multi-Core-Multi-Domain
///        simulations.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.27 JM: moved from sim_control_MPI.h
/// - 2016.02.02 JM: Added option to decompose only along one axis.
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef MCMD_CONTROL_H
#define MCMD_CONTROL_H

#include <constants.h>
#include <defines/functionality_flags.h>
#include <defines/testing_flags.h>
#include <iostream>
#include <mpi.h>
#include <sim_constants.h>
#include <vector>

// ##################################################################
// ##################################################################

/// struct to hold rank of child grids and Xmin/Xmax arrays.
struct cgrid {
  double Xmin[MAX_DIM];  ///< Xmin of the child grid.
  double Xmax[MAX_DIM];  ///< Xmax of the child grid.
  int rank;              ///< rank of process that holds child grid.
};

// ##################################################################
// ##################################################################

///
/// Class to hold information about Multi-Core-Multi-Domain
/// simulations, including the domain decomposition, the rank and
/// nproc, and the ranks of abutting domains.
///
/// For parallel processing, the domain is split between processors, so the
/// local domain is smaller than the full domain.  That information is stored
/// here.
///
class MCMDcontrol {
public:
  MCMDcontrol();

  ~MCMDcontrol();

  /// Number of real grid-zones in each direction, on this
  /// processor's domain.
  int LocalNG[MAX_DIM];

  /// Total number of real grid zones in this processor's domain.
  int LocalNcell;

  /// number of zones from this proc's 1st zone to the sim's 1st zone
  int offsets[MAX_DIM];
  double LocalXmin[MAX_DIM];   ///< Min value of x,y,z in Processor's domain.
  double LocalXmax[MAX_DIM];   ///< Max value of x,y,z in Processor's domain.
  double LocalRange[MAX_DIM];  ///< Size of Processor's domain in
                               ///< x,y,z-direction.

  std::vector<int> ngbprocs;  ///< list with processor rank of neighbours in
                              ///< each direction.

  bool ReadSingleFile;   ///< If the ICs are in a single file, set this to true.
  bool WriteSingleFile;  ///< If you want all the processors to write to one
                         ///< file, set this
  bool WriteFullImage;   ///< If you want multiple fits files, but each one is
                         ///< the full domain size

  /// choose the topology of subdomains to minimise surface area
  void calculate_process_topology(
      std::vector<float> &&  ///< ranges of grid data
  );
  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  ///
  int decomposeDomain(
      const int &,         ///< number of dimensions
      const class level &  ///< pointer to domain parameters for NG grid level
  );

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  ///
  int decomposeDomain(
      const int &,          ///< number of dimensions
      const class level &,  ///< pointer to domain parameters for NG grid level
      std::vector<int> &&pbc  ///< boolean array of whether each face has pbc
  );

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  /// *** This version decomposes only along one axis ***
  ///
  int decomposeDomain(
      const enum axes,
      const int &,         ///< number of dimensions
      const class level &  ///< domain parameters for NG grid level
  );

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  /// *** This version decomposes only along one axis ***
  ///
  int decomposeDomain(
      const enum axes,
      const int &,            ///< number of dimensions
      const class level &,    ///< domain parameters for NG grid level
      std::vector<int> &&pbc  ///< boolean array of whether each face has pbc
  );

  ///
  /// Create a vector of child grids of this processes on level + 1
  ///
  void determine_parent_processes(
      const class SimParams &par,  ///< simulation parameters
      const int level              ///< level in grid hierarchy
  );

  ///
  /// Determine the neighbours of this processes' children on level + 1
  /// in neighbour_dimension recursively
  ///
  void determine_child_neighbour_ranks(
      const class SimParams &par,   ///< simulation parameters
      const int level,              ///< level in grid hierarchy
      std::vector<double> &centre,  ///< centre of current local grid
      std::vector<std::vector<int> >
          &neighbours,  ///< resultant vector of child neighbour ranks
      const int current_dimension,  ///< current dimension to traverse
      const int
          neighbour_dimension  ///< the dimension to search for neighbours in
      ) const;

  ///
  /// Determine the neighbours of this processes' children on level + 1
  ///
  void determine_child_neighbours(
      const class SimParams &par,  ///< simulation parameters
      const int level              ///< level in grid hierarchy
  );

  ///
  /// Determine the ranks of this processes children on level + 1
  ///
  void determine_child_ranks(
      const class SimParams &par,   ///< simulation parameters
      const int level,              ///< level in grid hierarchy
      std::vector<int> &children,   ///< resultant vector of child process ranks
      std::vector<double> &centre,  ///< centre of current local grid
      const int current_dimension   ///< current dimension to traverse
      ) const;

  ///
  /// Create a vector of child grids of this processes on level + 1
  ///
  void determine_child_processes(
      const class SimParams &par,  ///< simulation parameters
      const int level              ///< level in grid hierarchy
  );

  ///
  /// For nested grid, set process ranks of parent and child grid(s).
  ///
  void set_NG_hierarchy(
      class SimParams &,  ///< simulation parameters
      const int           ///< level in grid hierarchy.
  );

  ///
  /// Get a list of all abutting domains, including corner/edge
  /// intersections.
  /// Helper function to recurse traversal of dimensions
  ///
  void create_abutting_domains_list(
      int,    ///< current dimension to traverse
      int[],  ///< position to traverse relative to
      bool);

  ///
  /// Get a list of all abutting domains, including corner/edge
  /// intersections.
  ///
  void create_abutting_domains_list();

  ///
  /// Returns the rank of the parent process
  ///
  int get_parent_proc() { return pgrid.rank; }

  ///
  /// Returns the ix array for any requested rank.
  ///
  void get_domain_coordinates(
      const int,  ///< rank ix requested for.
      int *       ///< array to put ix into.
      ) const;

  ///
  /// Return rank of process that has grid that contains the
  /// requested position on a given refinement level.
  ///
  int get_rank_from_grid_location(
      const class SimParams &,      ///< simulation parameters
      const std::vector<double> &,  ///< location sought
      const int                     ///< grid level
      ) const;

  /// get number of active dimensions
  int get_ndim() const { return m_ndim; }

  /// get my process rank
  int get_myrank() const { return myrank; }
  /// set my process rank
  void set_myrank(const int r)
  {
    myrank = r;
    return;
  }
  /// get number of MPI processes
  int get_nproc() const { return nproc; }
  /// set number of MPI processes
  void set_nproc(const int n)
  {
    nproc = n;
    return;
  }

  /// get nx[] array of subdomains
  void get_nx_subdomains(int *) const;

  /// get data on parent grid, if it exists
  void get_parent_grid_info(struct cgrid *) const;

  /// get data on neighbouring grids to parent grid, if they exist
  void get_parent_ngb_grid_info(std::vector<struct cgrid> &) const;

  /// get data on child grids, if they exist
  void get_child_grid_info(std::vector<struct cgrid> &) const;

  /// get data on neighbouring grids on level l+1, if they exist.
  void get_level_lp1_ngb_info(std::vector<std::vector<struct cgrid> > &) const;

  /// get the Cartesian communicator
  MPI_Comm get_communicator() const { return cart_comm; }

  /// gather each LocalNcell variable at the root process
  void gather_ncells(int *, const int &) const;

  /// gather each LocalNcell variable at each process
  void allgather_ncells(std::vector<int> &) const;

  /// gather LocalXmax and LocalXmin arrays at the root process
  void gather_extents(double *recv_buffer, const int &root) const;

  /// gather offsets arrays at the root process
  void gather_offsets(std::vector<int> &offsets_list, const int &root) const;

  /// gather localNG arrays at the root process
  void gather_localNG(std::vector<int> &localNG_list, const int &root) const;

  /// gather rank ordered list of number of local abutting domains and their
  /// displacements in the global list
  void gather_abutting_domains(
      std::vector<int> &abutting_domains_list,
      int *rank_displacements_in_abutting_domains_list,
      int *num_abutting_domains_by_rank,
      const int &root);

  void print_grid(std::vector<int> &coords, int current_dimension) const;

protected:
  int nproc;   ///< Number of processors.
  int myrank;  ///< Which processor am I?

  /// this proc's position in the block of domains (zero offset)
  int ix[MAX_DIM];

  ///< size of block of domains in each direction. (one block = one unit).
  int nx[MAX_DIM];

  /// booleans for whether a dimension has pbc
  std::vector<int> periodic;

  /// number of active dimensions
  int m_ndim;

  ///< list of abutting domains.
  std::vector<int> abutting_domains;

  /// data for the parent grid (xmin,xmax,etc)
  struct cgrid pgrid;

  /// data for neigbouring grids of parent process, if NG and l>0
  /// Up to 6 in 3D, 4 in 2D, 2 in 1D
  std::vector<struct cgrid> pgrid_ngb;

  /// list of level l+1 grids that share an external face with, but
  /// don't intersect with, my grid on level l.  The outer index of
  /// the 2D array is the outward normal direction of the grid face,
  /// and the inner one is the list of grids (up to 2**(ND-1)).
  std::vector<std::vector<struct cgrid> > cgrid_ngb;

  /// a process can have up to 2**NDIM child grids.
  std::vector<struct cgrid> child_procs;

  /// communicator created by MPI_Cart_create with Cartesian domain
  /// decomposition
  MPI_Comm cart_comm;
};

#endif  // MCMD_CONTROL_H
