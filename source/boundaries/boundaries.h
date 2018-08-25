/// \file boundaries.h
/// \brief constants and structures for boundary data
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2018.08.08 JM: moved code.

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include <string>
#include <list>


// ##################################################################
// ##################################################################


///
/// enum for the types of boundary condition.
///
enum BoundaryTypes {
    PERIODIC   = 1, ///< periodic bcs.
    OUTFLOW    = 2, ///< outflow or absorbing bcs.
    INFLOW     = 3, ///< inflow bcs., value on boundary doesn't change.
    REFLECTING = 4, ///< reflecting bcs., hard wall.
    FIXED      = 5, ///< fixed bcs, means every point on boundary has same unchanging value.
    JETBC      = 6, ///< A jet boundary, internal boundary.
    JETREFLECT = 7, ///< Sort-of reflection for bi-directional jet, 
                    ///< where normal B field passes through, but tangential is reversed.
    DMACH      = 8, ///< Outflow boundary for double mach reflection test problem only.
    DMACH2     = 9, ///< Fixed boundary on y=0, x in [0,1/6] fixed to postshock state.
    BCMPI      =10,
    ///< boundary between processor domains in parallel grid (on one
    ///< level)
    RADSHOCK   =11,
    ///< Boundary condition adjacent to cold wall for radiative shock
    ///< test problem.
    RADSH2     =12, ///< Outflow augmented with fixed outflow speed.
    ONEWAY_OUT =13,
    ///< One-way valve -- allows outflow but not inflow (zero
    ///< gradient OR reflecting).
    STWIND     =14,
    ///< Stellar wind sources exist, so apply internal boundaries.
    STARBENCH1 =15, ///< StarBench test for mixing with a solid wall.
    FINE_TO_COARSE  =16,
    ///< data in this cell should be obtained from finer-scale data
    ///< on a NG grid.
    COARSE_TO_FINE   =17,
    ///< data should be obtained from coarser level in a NG grid.
    FINE_TO_COARSE_SEND = 18,
    ///< boundary for sending data from this grid to coarser-level
    ///< grid on another MPI process.
    FINE_TO_COARSE_RECV = 19,
    ///< boundary for receiving data onto this grid from finer-level
    ///< grid on another MPI process.
    COARSE_TO_FINE_SEND = 20,
    ///< boundary for sending data from this grid to external
    ///< boundary of finer-level grid on another MPI process.
    COARSE_TO_FINE_RECV = 21
    ///< boundary for receiving data onto the external boundary of
    ///< this grid from a coarser-level grid on another MPI process.
};


// ##################################################################
// ##################################################################


//
// integer flags for MPI communication labels.
//
#define BC_ANYtag 0 ///< works for either sort of communication.
#define BC_MPItag 1 ///< Integer tag on MPI send/receive operations, to label that this communicates MPI boundary data.
#define BC_PERtag 2 ///< Integer tag to say it is for periodic BC.
#define BC_RTtag  3 ///< Integer tag to say we are transferring a radiative transfer column density tag.
#define BC_MPI_NG_tag 4 ///< Integer tag on MPI send/receive ops that transfer data between levels on a nested grid.


// ##################################################################
// ##################################################################

///
/// Struct to hold data on 2^ndim cells on a fine grid that should
/// be averaged to be sent to a coarser grid.
struct averaging {
  std::vector<cell *> c; ///< list of cells to be averaged
  pion_flt cpos[MAX_DIM];
  pion_flt *avg_state;
};


// ##################################################################
// ##################################################################


///
/// Struct to contain all the information for a grid boundary.
///
struct boundary_data {
  enum direction dir; ///< Outward Normal direction of boundary (NO dir if internal).
  enum direction ondir; ///< direction back onto grid.
  std::string type; ///< What type of boundary it is (Periodic, Absorbing, Fixed, Reflective, etc.).
  int itype;         ///< Integer flag for boundary type.
  int bloc;          ///< boundary location, e.g. x=0
  bool bpos;         ///< whether boundary is in +ve direction?
  enum axes baxis;   ///< index in position vector relating to bpos.
  std::list<cell *> data; ///< STL linked list for boundary data cells.
  ///
  /// STL linked list for grid cells to send to neighbouring
  /// processor (parallel only; for serial code this is unused).
  ///
  std::list<cell *> send_data;
  ///
  /// STL linked list for grid cells in a parent/child grid needed
  /// for the external boundaries of a child grid (unused for uniform
  /// grid) or the non-leaf data of a parent grid.
  ///
  std::list<cell*> NG;
  pion_flt *refval;  ///< Optional reference state vector.

  /// vector of data to be sent to coarser level grid (MPI-NG only)
  std::vector<struct averaging *> avg;
};

///
/// Struct to hold data on 2^ndim cells on a fine grid that should
/// be averaged to be sent to a coarser grid.
struct averaging {
  std::vector<cell *> c; ///< list of cells to be averaged
  pion_flt cpos[MAX_DIM];
  pion_flt *state;  ///< averaged state vector.

#endif // BOUNDARIES_H
