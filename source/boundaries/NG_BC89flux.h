/// \file NG_BC89flux.h
///
/// \brief Declares class that calculates fluxes leaving/entering
/// levels on a nested grid, to make sure that fluxes between
/// adjacent levels are consistent with each other.
///
/// \author Jonathan Mackey
///
/// Modifications:n
/// - 2019.12.04 JM: moved functions from uniform_grid and
///   sim_control_NG classes to gather everything in one file.

#ifndef NG_BC89FLUX_H
#define NG_BC89FLUX_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/uniform_grid.h"
#include "sim_constants.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

// ##################################################################
// ##################################################################

///
/// Struct for adding up the flux crossing a grid-boundary interface,
/// for ensuring conservation between different refinement levels in
/// a NG/AMR simulation.
///
struct flux_interface {
  std::vector<class cell *> c;  ///< list of cells with faces on this interface
  std::vector<double> area;  ///< area of each face for cells c, ordered as is c
  pion_flt *flux;            ///< flux through interface.
};

// ##################################################################
// ##################################################################

///
/// struct to hold all the interface fluxes, for ensuring
/// conservation between different refinement levels in
/// a NG/AMR simulation.  (Berger & Colella, 1989).
///
struct flux_update {
  vector<struct flux_interface *> fi;
  int Ncells;             ///< number of cells contributing.
  std::vector<int> rank;  ///< list of grid ranks to send to/recv from
  int dir;                ///< direction of face (outward for send)
  int ax;                 ///< axis of normal vector to face.
};

// ##################################################################
// ##################################################################

///
/// Class to set up flux arrays at level boundaries in a nested grid
/// structure, to save the fluxes after a calculation of dU(), and to
/// make the fluxes between levels consistent after each timestep.
/// The implementation follows the paper of Berger & Colella
/// (1989,JCP,82,64).
///
class NG_BC89flux {
public:
  NG_BC89flux();

  virtual ~NG_BC89flux();

  void setup_flux_vectors(const int  ///< number of levels in nested grid.
  );

  ///
  /// Setup the flux struct flux_update_recv with list of interfaces
  /// that need to be updated with fluxes from a finer level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual int setup_flux_recv(
      class SimParams &,      ///< simulation params (including BCs)
      class GridBaseClass *,  ///< pointer to coarse grid.
      const int               ///< finer level to receive from
  );

  ///
  /// Setup the flux struct flux_update_send with list of interfaces
  /// that need to be sent to a coarser level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual int setup_flux_send(
      class SimParams &,      ///< simulation params (including BCs)
      class GridBaseClass *,  ///< pointer to finer grid.
      const int               ///< level to send to
  );

  ///
  /// Add fluxes from boundary cells to grid structures, for cells
  /// at the grid boundary, to be sent to the coarser level grid.
  /// Assumes that fluxes have already been calculated and saved in
  /// array "F" of all relevant cells.
  ///
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).
  ///
  virtual void save_fine_fluxes(
      class SimParams &,  ///< simulation parameters
      const int           ///< level to save.
  );

  ///
  /// Add fluxes from internal cells to grid structures, for cells
  /// that sit above a grid boundary at a finer level.
  /// Assumes that fluxes have already been calculated and saved in
  /// array "F" of all relevant cells.
  ///
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels.
  ///
  virtual void save_coarse_fluxes(
      class SimParams &,  ///< simulation parameters
      const int           ///< level to save.
  );

  ///
  /// Receive fine-level fluxes at level boundary onto coarser parent
  /// grid(s) for static mesh refinement.  This version is for when the
  /// child and parent grid are on the same MPI process.
  ///
  virtual int recv_BC89_fluxes_F2C(
      class FV_solver_base *,  ///< spatial solver, for gradients
      class SimParams &,       ///< simulation parameters
      const int,               ///< My level in grid hierarchy.
      const int,               ///< TIMESTEP_FULL or TIMESTEP_FIRST_PART
      const int                ///< Full order of accuracy of simulation
  );

  /// For a given boundary, implement the BC89 flux-correction from
  /// fine to coarse grid, so that conserved quantities are conserved
  int recv_BC89_flux_boundary(
      class FV_solver_base *,  ///< spatial solver, for gradients
      class SimParams &,       ///< simulation parameters
      class GridBaseClass *,   ///< pointer to coarse grid
      const double,            ///< timestep
      struct flux_update &,    ///< data for fine grid
      struct flux_update &,    ///< data for coarse grid
      const unsigned int,      ///< direction of outward normal
      const axes               ///< axis of normal direction.
  );

protected:
  ///
  /// Add cells to the flux_update_recv and flux_update_send lists,
  /// for keeping track of flux entering/leaving a grid on one level.
  /// This is book-keeping for NG/AMR to ensure fluxes are
  /// consistent across all levels.
  ///
  int add_cells_to_face(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      enum direction,         ///< which direction we're facing
      int *,                  ///< xmin of interface region (integer units)
      int *,                  ///< xmax of interface region (integer units)
      int *,                  ///< number of elements in interface region
      const int,              ///< number of cells per face, per dim.
      struct flux_update &    ///< struct with list to populate
  );

  /// array of interfaces for a finer grid within this grid.
  /// Each element contains a single cell that needs to be updated.
  ///
  /// First dimension is for levels, second for directions on that
  /// level.  Directions correspond to outward normal of fine level.
  ///
  /// The list has each direction as an element of the outer array
  /// index, and each interface in that direction as the inner index.
  /// Note that it is not a matix because each direction doesn't
  /// necessarily have any elements.
  ///
  /// These fluxes are used to correct the fluxes on the coarse grid.
  ///
  std::vector<std::vector<struct flux_update> > flux_update_recv;

  /// Array of interfaces for a coarser grid encompassing this grid.
  /// Each element contains 2^(ndim-1) cells.
  ///
  /// First dimension is for levels, second for directions on that
  /// level.  Directions correspond to outward normal of fine level.
  ///
  /// The list has each direction as an element of the outer array
  /// index, and each interface in that direction as the inner index.
  /// Note that it is not a matix because each direction doesn't
  /// necessarily have any elements.
  ///
  /// These fluxes are used to correct the fluxes on the coarse grid.
  ///
  std::vector<std::vector<struct flux_update> > flux_update_send;
};

#endif  // NG_BC89FLUX_H
