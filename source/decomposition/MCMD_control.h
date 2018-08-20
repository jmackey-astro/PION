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

  bool ReadSingleFile; ///< If the ICs are in a single file, set this to true.
  bool WriteSingleFile; ///< If you want all the processors to write to one file, set this (BUGGY!)
  bool WriteFullImage;  ///< If you want multiple fits files, but each one is the full domain size (bad!)

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  ///
  int decomposeDomain(
      class SimParams &,  ///< pointer to simulation parameters
      struct level &     ///< pointer to domain parameters for SMR grid level
      );

  ///
  /// Decompose the domain into blocks for each processor, and set up
  /// a structure which contains the domain of each processor.
  /// *** This version decomposes only along one axis ***
  ///
  int decomposeDomain(
      const enum axes,
      class SimParams &,  ///< pointer to simulation parameters
      struct level &     ///< pointer to domain parameters for SMR grid level
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

  /// get my process rank
  int get_myrank() {return myrank;}
  /// set my process rank
  void set_myrank(const int r) {myrank=r; return;}
  /// get number of MPI processes
  int get_nproc()  {return nproc;}
  /// set number of MPI processes
  void set_nproc(const int n)  {nproc=n; return;}


 protected:
  int nproc;  ///< Number of processors.
  int myrank; ///< Which processor am I?

  /// this proc's position in the block of domains (zero offset)
  int ix[MAX_DIM];

  ///< size of block of domains in each direction. (one block = one unit).
  int nx[MAX_DIM];

  ///< list of abutting domains.
  std::vector<int> full_ngb_list;

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
      struct level &     ///< pointer to domain parameters for SMR grid level
      );
};

#endif // MCMD_CONTROL_H

