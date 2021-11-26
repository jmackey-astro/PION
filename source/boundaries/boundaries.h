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
#include "grid/cell_interface.h"
#include <list>
#include <string>
#include <vector>

// ##################################################################
// ##################################################################

//#define GLM_ZERO_BOUNDARY // Set this flag to make Psi=0 on boundary cells.
#define GLM_NEGATIVE_BOUNDARY  // Set this flag for Psi[boundary cell]=-Psi[edge
                               // cell]

//#define TEST_MPI_NG
//#define TEST_C2F
//#define TEST_MPI_NG_F2C
//#define TEST_MPI_BC

///
/// enum for the types of boundary condition.
///
enum BoundaryTypes {
  PERIODIC   = 1,  ///< periodic bcs.
  OUTFLOW    = 2,  ///< outflow or absorbing bcs.
  INFLOW     = 3,  ///< inflow bcs., value on boundary doesn't change.
  REFLECTING = 4,  ///< reflecting bcs., hard wall.
  FIXED      = 5,  ///< fixed bcs, means every point on boundary has same
                   ///< unchanging value.
  JETBC      = 6,  ///< A jet boundary, internal boundary.
  JETREFLECT = 7,  ///< Sort-of reflection for bi-directional jet,
                   ///< where normal B field passes through, but tangential is
                   ///< reversed.
  DMACH  = 8,      ///< Outflow boundary for double mach reflection.
  DMACH2 = 9,      ///< Fixed boundary on y=0, x in [0,1/6].
  BCMPI  = 10,
  ///< boundary between processor domains in parallel grid (on one
  ///< level)
  RADSHOCK = 11,
  ///< Boundary condition adjacent to cold wall for radiative shock
  ///< test problem.
  RADSH2     = 12,  ///< Outflow augmented with fixed outflow speed.
  ONEWAY_OUT = 13,
  ///< One-way valve -- allows outflow but not inflow (zero
  ///< gradient OR reflecting).
  STWIND = 14,
  ///< Stellar wind sources exist, so apply internal boundaries.
  STARBENCH1     = 15,  ///< StarBench test for mixing with a solid wall.
  FINE_TO_COARSE = 16,
  ///< data in this cell should be obtained from finer-scale data
  ///< on a NG grid.
  COARSE_TO_FINE = 17,
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
  COARSE_TO_FINE_RECV = 21,
  ///< boundary for receiving data onto the external boundary of
  ///< this grid from a coarser-level grid on another MPI process.
  AXISYMMETRIC = 22
  ///< reflecting about R=0 symmetry axis, but assumes rotational
  ///< symmetry for theta-component of V and B.
};

// ##################################################################
// ##################################################################

//
// integer flags for MPI communication labels.
//
#define BC_ANYtag 0  ///< works for any sort of communication.
#define BC_MPItag                                                              \
  1  ///< Integer tag on MPI send/receive operations, to label that this
     ///< communicates MPI boundary data.
#define BC_PERtag 2  ///< Integer tag to say it is for periodic BC.
#define BC_RTtag                                                               \
  3  ///< Integer tag to say we are transferring a radiative transfer column
     ///< density tag.
#define BC_MPI_NGF2C_tag 10    ///< MPI send/receive from fine to coarse grid.
#define BC_MPI_NGC2F_tag 1000  ///< MPI send/receive from coarse to fine grid.
// this gives a unique identifier as long as n_levels<10, nproc<10000.
#define BC_MPI_FLUX_tag 10000  ///< Berger & Colella (1989) flux correction.

// ##################################################################
// ##################################################################

///
/// Struct to hold data on 2^ndim cells on a fine grid that should
/// be averaged to be sent to a coarser grid.
struct averaging {
  std::vector<cell *> c;  ///< list of cells to be averaged
  std::array<pion_flt, MAX_DIM>
      cpos;  ///< position of cell-centre of coarse cell
  pion_flt *avg_state;
};

// ##################################################################
// ##################################################################

/// struct to hold cells that should be sent to a finer-level grid
struct c2f {
  std::vector<cell *> c;  ///< list of cells to be sent/received
  int rank;               ///< rank of process to send to/receive from.
  int dir;                ///< location of boundary on finer-level grid.
};

// ##################################################################
// ##################################################################

///
/// Struct to contain all the information for a grid boundary.
///
struct boundary_data {
  enum direction
      dir;  ///< Outward Normal direction of boundary (NO dir if internal).
  enum direction ondir;  ///< direction back onto grid.
  std::string type;      ///< What type of boundary it is (Periodic, Absorbing,
                         ///< Fixed, Reflective, etc.).
  int itype;             ///< Integer flag for boundary type.
  int bloc;              ///< boundary location, e.g. x=0
  int depth;             ///< how many cells deep the boundary is.
  bool bpos;             ///< whether boundary is in +ve direction?
  enum axes baxis;       ///< index in position vector relating to bpos.
  std::list<cell *> data;  ///< STL linked list for boundary data cells.
  ///
  /// STL linked list for grid cells to send to neighbouring
  /// MPI process (parallel only), on the same grid level.
  ///
  std::list<cell *> send_data;
  /// vector of data to be sent to coarser level grid (MPI-NG only)
  std::vector<struct averaging> avg;
  /// Vector of lists of cells, for a coarse grid that
  /// receives data from a number of child grids to replace the
  /// on-grid data.  Vector length is the number of children.
  /// - Used by serial NG code:
  ///    - setup_NG_grid::setup_boundary_structs()
  ///    - NG_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE()
  ///    - NG_fine_to_coarse_bc::BC_update_FINE_TO_COARSE(
  ///    - sim_control_NG::do_ongrid_raytracing()
  /// - Used by parallel NG code:
  ///    - NG_MPI_fine_to_coarse_bc::BC_assign_FINE_TO_COARSE_RECV()
  ///    - NG_MPI_fine_to_coarse_bc::BC_update_FINE_TO_COARSE_RECV()
  std::vector<std::list<cell *> > NGrecvF2C;
  std::vector<int> NGrecvF2C_ranks;
  /// as NGrecvF2C, but C2F
  std::vector<std::list<cell *> > NGrecvC2F;
  int NGrecvC2F_parent;  ///< process containing coarse level data.
  /// (MPI-NG only) list of lists of cells, for a coarse grid that
  /// sends data to a number of child grids for their external
  /// boundaries.  Vector length is the number of child boundaries
  /// to update.
  std::vector<struct c2f *> NGsendC2F;
  pion_flt *refval;  ///< Optional reference state vector.
};

#endif  // BOUNDARIES_H
