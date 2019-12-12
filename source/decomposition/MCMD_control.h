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

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_constants.h"
#include "constants.h"
#include <vector>



// ##################################################################
// ##################################################################



/// struct to hold rank of child grids and Xmin/Xmax arrays.
struct cgrid {
  double Xmin[MAX_DIM]; ///< Xmin of the child grid.
  double Xmax[MAX_DIM]; ///< Xmax of the child grid.
  int rank; ///< rank of process that holds child grid.
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
  double LocalXmin[MAX_DIM];  ///< Min value of x,y,z in Processor's domain.
  double LocalXmax[MAX_DIM];  ///< Max value of x,y,z in Processor's domain.
  double LocalRange[MAX_DIM]; ///< Size of Processor's domain in x,y,z-direction.

  int *ngbprocs;  ///< list with processor rank of neighbours in each direction.

  int parent_proc; ///< process of the parent grid, if NG and l>0
  /// a process can have up to 2**NDIM child grids.
  std::vector<struct cgrid> child_procs;

  bool ReadSingleFile; ///< If the ICs are in a single file, set this to true.
  bool WriteSingleFile; ///< If you want all the processors to write to one file, set this (BUGGY!)
  bool WriteFullImage;  ///< If you want multiple fits files, but each one is the full domain size (bad!)

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  ///
  int decomposeDomain(
      class SimParams &,  ///< pointer to simulation parameters
      class level &     ///< pointer to domain parameters for NG grid level
      );

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  /// *** This version decomposes only along one axis ***
  ///
  int decomposeDomain(
      const enum axes,
      class SimParams &,  ///< simulation parameters
      class level &     ///< domain parameters for NG grid level
      );

  ///
  /// For nested grid, set process ranks of parent and child grid(s).
  ///
  void set_NG_hierarchy(
      class SimParams &, ///< simulation parameters
      const int  ///< level in grid hierarchy.
      );

  ///
  /// Get a list of all abutting domains, including corner/edge
  /// intersections.
  ///
  void get_abutting_domains(
      const int,  ///< grid dimensions
      std::vector<int> & ///< write list to this vector.
      );

  ///
  /// Returns the ix array for any requested rank.
  ///
  void get_domain_ix(
      const int,  ///< grid dimensions
      const int, ///< rank ix requested for.
      int *      ///< array to put ix into.
      );

  ///
  /// Return rank of process that has grid that contains the
  /// requested position on a given refinement level.
  ///
  int get_grid_rank(
      class SimParams &,  ///< simulation parameters
      const double *, ///< location sought
      const int     ///< grid level
      );

  /// get my process rank
  int get_myrank() {return myrank;}
  /// set my process rank
  void set_myrank(const int r) {myrank=r; return;}
  /// get number of MPI processes
  int get_nproc()  {return nproc;}
  /// set number of MPI processes
  void set_nproc(const int n)  {nproc=n; return;}

  /// get data on parent grid, if it exists
  void get_parent_grid_info(struct cgrid *cg) {cg = &pgrid;return;}

 protected:
  int nproc;  ///< Number of processors.
  int myrank; ///< Which processor am I?

  /// this proc's position in the block of domains (zero offset)
  int ix[MAX_DIM];

  ///< size of block of domains in each direction. (one block = one unit).
  int nx[MAX_DIM];

  ///< list of abutting domains.
  std::vector<int> full_ngb_list;

  /// data for the parent grid (xmin,xmax,etc)
  struct cgrid  pgrid;

  /// data for neigbouring grids of parent process, if NG and l>0
  /// Up to 6 in 3D, 4 in 2D, 2 in 1D 
  std::vector<struct cgrid>  pgrid_ngb;

  /// list of level l+1 grids that share an external face with, but 
  /// don't intersect with, my grid on level l.  The outer index of
  /// the 2D array is the outward normal direction of the grid face,
  /// and the inner one is the list of grids (up to 2**(ND-1)).
  std::vector< std::vector<struct cgrid> >  cgrid_ngb;


  ///
  /// Called by decomposeDomain() to set neighbouring processor ids.
  /// 
  /// This works if processors are ranked from zero to n-1 in integer steps.
  /// If a subdomain is on the full domain boundary, the processor neighbour
  /// in that direction is set to -999.
  /// 
  /// The rank of processor i is defined by
  /// \f[ \mbox{myrank} = n_x*n_y*i_z + n_x*i_y + i_x \f]
  /// where the 'n's refer to number of processors that span the domain
  /// in each direction, and 'i's refer to the location of the current
  /// processor along that direction.
  ///
  int pointToNeighbours(
      class SimParams &,  ///< pointer to simulation parameters
      class level &     ///< pointer to domain parameters for NG grid level
      );
};

#endif // MCMD_CONTROL_H

