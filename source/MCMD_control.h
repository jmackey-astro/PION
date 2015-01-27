/// \file MCMD_control.h
/// 
/// \brief Declares class for controlling Multi-Core-Multi-Domain
///        simulations.
/// 
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.01.27 JM: moved from sim_control_MPI.h

#ifndef MCMD_CONTROL_H
#define MCMD_CONTROL_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "sim_constants.h"
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

  int nproc;  ///< Number of processors.
  int myrank; ///< Which processor am I?

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
  int decomposeDomain();

  ///
  /// Get a list of all abutting domains, including corner/edge
  /// intersections.
  ///
  void get_abutting_domains(std::vector<int> & ///< write list to this vector.
          );

  ///
  /// Returns the ix array for any requested rank.
  ///
  void get_domain_ix(const int, ///< rank ix requested for.
         int *      ///< array to put ix into.
         );

 protected:
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
  int pointToNeighbours();
};

#endif // MCMD_CONTROL_H

