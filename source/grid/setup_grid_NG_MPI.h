/// \file setup_grid_NG_MPI.h
///
/// \brief Declares a class for setting up NG-MPI grids.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.09.05 JM: worked on code

#ifndef SETUP_GRID_NG_MPI_H
#define SETUP_GRID_NG_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/NG_MPI_BC89flux.h"
#include "boundaries/assign_update_bcs_NG_MPI.h"
#include "decomposition/MCMD_control.h"
#include "grid/grid_base_class.h"
#include "grid/setup_NG_grid.h"
#include "grid/setup_fixed_grid_MPI.h"
#include "grid/uniform_grid.h"
#include "setup_fixed_grid.h"
#include "spatial_solvers/solver_eqn_base.h"

#include "dataIO/dataio_base.h"
#ifdef SILO
#include "dataIO/dataio_silo_utility.h"
#endif  // if SILO
#ifdef FITS
#include "dataIO/dataio_fits.h"
#endif  // if FITS

///
/// Set up a static NG grid structure.  Serial code, so each
/// level of the NG has a single grid.
///
class setup_grid_NG_MPI :
    virtual public assign_update_bcs_NG_MPI,
    virtual public setup_NG_grid,
    virtual public setup_fixed_grid_pllel,
    virtual public NG_MPI_BC89flux {
public:
  setup_grid_NG_MPI();
  ~setup_grid_NG_MPI();

  ///
  /// Populate the array SimPM.levels with Xmin,Xmax,Range,dx,etc.
  ///
  void setup_NG_grid_levels(
      class SimParams&  ///< pointer to simulation parameters
  );

  ///
  /// Sets up a NG grid.
  ///
  int setup_grid(
      vector<class GridBaseClass*>&,  ///< vector of grids.
      class SimParams&                ///< pointer to simulation parameters
  );

  ///
  /// Decide if I need to setup RT class and, if so, set up a
  /// raytracer associated with each grid.
  ///
  int setup_raytracing(
      class SimParams&,              ///< pointer to simulation parameters
      vector<class GridBaseClass*>&  ///< address of vector of grid.
  );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  int boundary_conditions(
      class SimParams&,              ///< pointer to simulation parameters
      vector<class GridBaseClass*>&  ///< address of vector of grid.
  );

  //---------------------------------------
protected:
  //---------------------------------------

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams&,      ///< reference to SimParams list.
      class GridBaseClass*,  ///< pointer to grid.
      const int              ///< level of grid in NG
  );

  /// function to setup parallel data-I/O class.
  void setup_dataio_class(
      class SimParams&,  ///< pointer to simulation parameters
      const int          ///< type of I/O: 2=fits,5=silo
  );

};  // setup_grid_NG_MPI

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif  // if not SETUP_GRID_NG_MPI_H
