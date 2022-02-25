/// \file sub_domain/sub_domain.h
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

#ifndef sub_domain_CONTROL_H
#define sub_domain_CONTROL_H

#include <constants.h>
#include <defines/functionality_flags.h>
#include <defines/testing_flags.h>

#include <list>
#include <mpi.h>
#include <sim_constants.h>
#include <string>
#include <unordered_set>
#include <vector>

#include "grid/cell_interface.h"

#ifdef SILO
#include <silo.h>
/* prevent clang format from reordering these headers, undefined references
 * otherwise */
#include <pmpio.h>
#endif

// ##################################################################
// ##################################################################

/// struct to hold rank of child grids and Xmin/Xmax arrays.
struct cgrid {
  std::array<double, MAX_DIM> Xmin;  ///< Xmin of the child grid.
  std::array<double, MAX_DIM> Xmax;  ///< Xmax of the child grid.
  int rank;                          ///< rank of process that holds child grid.
};

// ##################################################################
// ##################################################################

#define COMM_CELLDATA 1
#define COMM_DOUBLEDATA 2

struct Send_info {
  MPI_Request request;  ///< MPI handle for the send.
  std::string id;       ///< string that code uses as handle for send.
  int comm_tag;         ///< the tag used to describe the send.
  int from_rank;        ///< sender.
  int to_rank;          ///< recipient
  int type;             ///< type of data: COMM_CELLDATA=char array;
  std::vector<char> send_buff;
  std::vector<double> send_buff_double;
};

struct Recv_info {
  MPI_Status status;  ///< MPI handle for the receive.
  std::string id;     ///< string that code uses as handle for send.
  int comm_tag;       ///< the tag used to describe the send.
  int from_rank;      ///< sender.
  int to_rank;        ///< recipient
  void *data;         ///< pointer to data.
  int type;           ///< type of data: COMM_CELLDATA=char array;
};

// ##################################################################
// ##################################################################

enum mpi_op { MAX, MIN, SUM };
enum mpi_type { INT, FLOAT, DOUBLE };

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
class Sub_domain {
public:
  Sub_domain();
  ~Sub_domain();

  Sub_domain(const Sub_domain &sub_domain) : Sub_domain(){};
  Sub_domain(Sub_domain &sub_domain) : Sub_domain(){};

  /* Setup functions */

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
  void evaluate_parent_process(
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
      std::vector<double> &cursor,  ///< cursor for exploring domain
      std::array<std::unordered_set<int>, 2>
          &neighbours,  ///< resultant vector of child neighbour ranks
      const int current_dimension,  ///< current dimension to traverse
      const int
          neighbour_dimension  ///< the dimension to search for neighbours in
      ) const;

  ///
  /// Determine the neighbours of this processes' children on level + 1
  ///
  void evaluate_child_neighbours(
      const class SimParams &par,  ///< simulation parameters
      const int level              ///< level in grid hierarchy
  );

  ///
  /// Determine the ranks of this processes children on level + 1
  ///
  void determine_child_ranks(
      const class SimParams &par,  ///< simulation parameters
      const int level,             ///< level in grid hierarchy
      std::unordered_set<int>
          &children,                ///< resultant set of child process ranks
      std::vector<double> &cursor,  ///< cursor for exploring domain
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
  ///
  void create_abutting_domains_list();


  /* Complex accessors */

  ///
  /// Returns the coordinates array for any requested rank.
  ///
  void get_domain_coordinates(
      const int,  ///< rank coordinates requested for.
      int *       ///< array to put coordinates into.
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

  /// get num_subdomains[] array of subdomains
  void get_num_subdomains(int *) const;

  /// get data on parent grid, if it exists
  void get_parent_grid_info(struct cgrid *) const;

  /// get data on neighbouring grids to parent grid, if they exist
  void get_parent_ngb_grid_info(std::vector<struct cgrid> &) const;

  /// get data on child grids, if they exist
  void get_child_grid_info(std::vector<struct cgrid> &) const;

  /// get data on neighbouring grids on level l+1, if they exist.
  void get_level_lp1_ngb_info(std::vector<std::vector<struct cgrid> > &) const;

  /// gather each Ncell variable at the root process
  void gather_ncells(int *, const int &) const;

  /// gather each Ncell variable at each process
  void allgather_ncells(std::vector<int> &) const;

  /// gather Xmax and Xmin arrays at the root process
  void gather_extents(double *recv_buffer, const int &root) const;

  /// gather offsets arrays at the root process
  void gather_offsets(std::vector<int> &offsets_list, const int &root) const;

  /// gather directional_Ncells arrays at the root process
  void gather_directional_Ncells(
      std::vector<int> &directional_Ncells_list, const int &root) const;

  /// gather rank ordered list of number of local abutting domains and their
  /// displacements in the global list
  void gather_abutting_domains(
      std::vector<int> &abutting_domains_list,
      int *rank_displacements_in_abutting_domains_list,
      int *num_abutting_domains_by_rank,
      const int &root);

  /* Member accessors */

  ///
  /// Returns the rank of the parent process
  ///
  int get_parent_proc() const { return pgrid.rank; }

  /// get number of active dimensions
  int get_ndim() const { return m_ndim; }

  /// get my process rank
  int get_myrank() const { return myrank; }

  /// get number of MPI processes
  int get_nproc() const { return nproc; }

  /// returns the number of real grid-zones in each direction, on this
  /// processor's domain.
  const std::array<int, MAX_DIM> get_directional_Ncells() const
  {
    return directional_Ncells;
  }

  /// returns the number of real grid-zones in a given direction, on this
  /// processor's domain.
  int get_directional_Ncells(const int i) const
  {
    return directional_Ncells[i];
  }

  /// get the total number of real grid zones in this processor's domain.
  int get_Ncell() const { return Ncell; }

  ///< get min value of x,y,z in Processor's domain.
  const std::array<double, MAX_DIM> get_Xmin() const { return Xmin; }

  ///< get min value of x,y, or z in Processor's domain.
  double get_Xmin(const int i) const { return Xmin[i]; }

  ///< get max value of x,y,z in Processor's domain.
  const std::array<double, MAX_DIM> get_Xmax() const { return Xmax; }

  ///< get max value of x,y, or z in Processor's domain.
  double get_Xmax(const int i) const { return Xmax[i]; }

  /// get list with processor rank of neighbours in
  /// each direction.
  std::vector<int> get_neighbour_ranks() const { return neighbour_ranks; }

  /// get processor rank of neighbour in a given direction
  int get_neighbour_rank(const int i) const { return neighbour_ranks[i]; }

  /// Are the ICs are in a single file?
  bool get_ReadSingleFile() const { return ReadSingleFile; }

  /// Should all the processors write to one file?
  bool get_WriteSingleFile() const { return WriteSingleFile; }

  /// Should there be multiple fits files, with each being the full domain size?
  bool get_WriteFullImage() const { return WriteFullImage; }

  /* Member modifiers */

  /// set my process rank
  void set_myrank(const int r) { myrank = r; }

  /// set number of MPI processes
  void set_nproc(const int n) { nproc = n; }

  /// Set the number of real grid-zones in a given direction, on this
  /// processor's domain.
  void set_directional_Ncells(const int i, const int n)
  {
    directional_Ncells[i] = n;
  }

  /// Set whether the ICs are in a single file
  void set_ReadSingleFile(const bool b) { ReadSingleFile = b; }

  /// Set whether all the processors write to one file
  void set_WriteSingleFile(const bool b) { WriteSingleFile = b; }

  /// Set whether there are multiple fits files, with each being the full domain
  /// size
  void set_WriteFullImage(const bool b) { WriteFullImage = b; }

  /* Debug */
  void print_grid(std::vector<int> &coords, int current_dimension) const;

private:
  int nproc;   ///< Number of processors.
  int myrank;  ///< Which processor am I?

  /// this proc's position in the block of domains (zero offset)
  std::array<int, MAX_DIM> coordinates;

  ///< size of block of domains in each direction. (one block = one unit).
  std::array<int, MAX_DIM> num_subdomains;

  /// booleans for whether a dimension has pbc
  std::vector<int> periodic;

  /// number of active dimensions
  int m_ndim;

  /// Number of real grid-zones in each direction, on this
  /// processor's domain.
  std::array<int, MAX_DIM> directional_Ncells;

  /// Total number of real grid zones in this processor's domain.
  int Ncell;

  /// number of zones from this proc's 1st zone to the sim's 1st zone
  std::array<int, MAX_DIM> offsets;
  std::array<double, MAX_DIM>
      Xmin;  ///< Min value of x,y,z in Processor's domain.
  std::array<double, MAX_DIM>
      Xmax;  ///< Max value of x,y,z in Processor's domain.
  std::array<double, MAX_DIM> range;  ///< Size of Processor's domain in
                                      ///< x,y,z-direction.

  std::vector<int> neighbour_ranks;  ///< list with processor rank of neighbours
                                     ///< in each direction.

  bool ReadSingleFile;   ///< If the ICs are in a single file, set this to true.
  bool WriteSingleFile;  ///< If you want all the processors to write to one
                         ///< file, set this
  bool WriteFullImage;   ///< If you want multiple fits files, but each one is
                         ///< the full domain size

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
  MPI_Comm cart_comm = MPI_COMM_WORLD;
  /// whether the cartesion communicator has been created
  bool m_is_cartcomm = false;

  /* MPI communication */
public:
  /// List of IDs for fine-to-coarse sends from this level, should be cleared at
  /// the beginning of each timestep.
  std::vector<std::string> NG_F2C_send_list;

  /// List of IDs for coarse-to-fine sends from this level, should be cleared at
  /// the beginning of each timestep.
  std::vector<std::string> NG_C2F_send_list;

  /// List of IDs for MPI sends related to the BC89 Flux
  /// correction algorithms.  Should be cleared at the beginning
  /// of each timestep.
  std::vector<std::string> BC89_flux_send_list;

  /// Tell other processes to abort!
  int abort() const;

  /// Set up a barrier, and return when all processes have also
  /// set up their own barrier.
  int barrier() const;

  /// Do a global operation on local data, options are MAX,MIN,SUM.
  /// Returns the global value when local value is passed in.
  ///
  double global_operation_double(
      const mpi_op,  ///< MAX,MIN,SUM
      const double   ///< this process's local value.
      ) const;

  ///
  /// Do a global operation on local data, options are MAX,MIN,SUM.
  /// Receives in local value, performs global operation on it,
  /// and returns global value.  This function works with
  /// N-element arrays, and should work "in-place" i.e. send and recv
  /// buffers are the same.
  ///
  void global_op_double_array(
      const mpi_op,  ///< MAX,MIN,SUM
      const size_t,  ///< Number of elements in array.
      double *       ///< pointer to this process's data array.
      ) const;

  /// Broadcast data from one process to all others.
  int broadcast_data(
      const int,       ///< rank of sender.
      const mpi_type,  ///< Type of data INT,DOUBLE,etc.
      const int,       ///< number of elements
      void *           ///< pointer to data.
      ) const;

  /// Send cell data to another processor, return immediately, but
  /// keep a record of the send for later to know the send has been
  /// received.
  ///
  /// receives a list of cells, puts data from each cell into a
  /// buffer to send, sends it to another processor, which should
  /// know how to unpack the data from it's matching receive call.
  ///
  /// Ccopies the data into its own buffer, so the user is
  /// free to delete the list of cells even if the send is not
  /// complete.
  ///
  int send_cell_data(
      const int,            ///< rank to send to.
      std::list<cell *> &,  ///< list of cells to get data from.
      const int,            ///< nvar
      std::string &,        ///< identifier for send, delivery tracking
      const int             ///< comm_tag, kind of send this is.
  );

  ///
  /// Send an array of n doubles to another processor, return but
  /// keep a record for later when the send has been received.
  ///
  /// copies the data into its own buffer, so the user is
  /// free to delete the data on return, before the send is complete.
  ///
  int send_double_data(
      const int,              ///< rank to send to.
      const long int,         ///< size of buffer, in number of doubles.
      std::vector<double> &,  ///< vector of data
      std::string &,          ///< identifier for send, for tracking delivery
      const int               ///< comm_tag, for what kind of send this is.
  );

  /// Called when we need to make sure a send has been received.
  /// It returns once the receiver confirms that it has got the data.
  ///
  int wait_for_send_to_finish(
      std::string &  ///< identifier for the send we are waiting on.
  );

  /// Look for some data that is being sent to us.  Does not return
  /// until it finds some, so it is up to you to prevent deadlock!
  ///
  int look_for_data_to_receive(
      int *,          ///< rank of sender
      std::string &,  ///< identifier for receive.
      int *,          ///< comm_tag associated with data.
      const int,      ///< comm_tag requested: (PER,MPI,F2C,C2F)
      const int       ///< type of data to look for (e.g. COMM_DOUBLEDATA)
  );

  /// Receive Cell data from a specific process rank.
  ///
  /// It is up to the caller to make sure that data was packed in the
  /// right way, and that the cells were sent in the same order that
  /// they are received.  No checking of this is performed.
  ///
  int receive_cell_data(
      const int,            ///< rank of process we are receiving from.
      std::list<cell *> &,  ///< list of cells to get data for.
      const int,            ///< nvar
      const int,            ///< comm_tag: (PER,MPI,F2C,C2F)
      const std::string &   ///< identifier for receive.
  );

  /// Receive array of doubles from a specific process rank.
  int receive_double_data(
      const int,             ///< rank of process we are receiving from.
      const int,             ///< comm_tag: (PER,MPI,F2C,C2F)
      const std::string &,   ///< identifier for receive.
      const long int,        ///< number of doubles to receive
      std::vector<double> &  ///< Pointer to array to write to (initialised).
  );

#ifdef SILO
  int silo_pllel_init(
      const int,          ///< number of files to write
      const std::string,  ///< READ or WRITE
      const std::string,  ///< identifier for this read/write.
      int *,              ///< group rank.
      int *               ///< rank in group
  );

  int silo_pllel_wait_for_file(
      const std::string,  ///< identifier for this read/write.
      const std::string,  ///< File Name
      const std::string,  ///< Directory to open in file.
      DBfile **           ///< pointer that file gets returned in.
  );

  int silo_pllel_finish_with_file(
      const std::string,  ///< identifier for this read/write.
      DBfile **           ///< pointer to file we have been working on.
  );

  void silo_pllel_get_ranks(
      std::string id,  ///< identifier for this read/write.
      const int,       ///< processor rank
      int *,           ///< rank of group processor is in.
      int *            ///< rank of processor within group.
      ) const;

#endif

private:
  std::list<struct Send_info> send_list;
  std::list<struct Recv_info> recv_list;

  long int m_num_send_cells = 0;

#ifdef SILO
  PMPIO_baton_t *m_baton;
  std::string m_silo_id;
#endif

  static unsigned int m_count;
};

#endif  // sub_domain_CONTROL_H
