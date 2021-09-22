/// \file setup_fixed_grid.h
///
/// \brief Declares a class for setting up grids.
///
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2015.02.09 JM: Split sim_control class into a setup class and
///   a derived class for running simulations.
/// - 2017.08.24 JM: moved evolving_RT_sources functions to setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef SETUP_FIXED_GRID_H
#define SETUP_FIXED_GRID_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/assign_update_bcs.h"
#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "spatial_solvers/solver_eqn_base.h"

#define TIMESTEP_FULL 2
#define TIMESTEP_FIRST_PART 1

///
/// The simplest finite volume grid - a uniform grid with cells that
/// are cubes in the chosen coordinates.  This sets up the grid and
/// other things to get a simulation ready to run, so it is useful
/// for simulation analysis.  PION itself uses a derived class to
/// setup and run simulations.
///
class setup_fixed_grid : virtual public assign_update_bcs {
public:
  setup_fixed_grid();
  virtual ~setup_fixed_grid();

  ///
  /// Setup cell extra data through the cell_interface class CI.
  ///
  void setup_cell_extra_data(
      class SimParams &  ///< pointer to simulation parameters
  );

  /// \brief initialise the grid class with appropriate parameters.
  ///
  /// This function sets up the appropriate grid; so far I only have a
  /// UniformGrid class -- uniform, cartesian, finite volume grid.
  ///
  virtual int setup_grid(
      vector<class GridBaseClass *> &,  ///< grid pointers.
                                        // class GridBaseClass **, ///< address
                                        // of pointer to grid.
      class SimParams &                 ///< pointer to simulation parameters
  );

  ///
  /// Decide if I need to setup MP class, and do it if i need to.
  ///
  virtual int setup_microphysics(
      class SimParams &  ///< pointer to simulation parameters
  );

  ///
  /// Decide if I need to setup RT class, and do it if i need to.
  ///
  virtual int setup_raytracing(
      class SimParams &,     ///< pointer to simulation parameters
      class GridBaseClass *  ///< pointer to computational grid
  );

  ///
  /// Check for any time-evolving radiation sources, and read
  /// the evolution file if there are any.  Data is stored in
  /// global struct SimPM.STAR[v]
  ///
  virtual int setup_evolving_RT_sources(
      class SimParams &  ///< pointer to simulation parameters
  );

  ///
  /// Check for any time-evolving radiation sources, and update source
  /// properties from global struct SimPM.STAR[v] if needed
  ///
  virtual int update_evolving_RT_sources(
      class SimParams &,      ///< pointer to simulation parameters
      const double,           ///< current time
      class RayTracingBase *  ///< pointer to raytracing class [OUTPUT]
  );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  virtual int boundary_conditions(
      class SimParams &,               ///< pointer to simulation parameters
      vector<class GridBaseClass *> &  ///< grid pointers.
      // class GridBaseClass *  ///< pointer to grid.
  );

  ///
  /// Initialise the correct Equations to solve, based on paramters.
  ///
  int set_equations(class SimParams &  ///< simulation parameters
  );

  ///
  /// Get pointer to equations class.
  ///
  class FV_solver_base *get_solver_ptr() { return spatial_solver; }

  //---------------------------------------
protected:
  //---------------------------------------

  /// flag: true if calc_timestep needs raytracing column densities
  bool FVI_need_column_densities_4dt;

  int FVI_nheat;  ///< number of RT heating sources

  int FVI_nion;  ///< number of ionising sources

  /// vector of RT heating sources, of size FVI_nheat
  std::vector<struct rt_source_data> FVI_heating_srcs;

  /// vector of RT ionising sources, of size FVI_nion
  std::vector<struct rt_source_data> FVI_ionising_srcs;

  ///
  /// Pointer to equations to solve, and routines for calculating
  /// fluxes on the grid.
  ///
  class FV_solver_base *spatial_solver;

  ///
  /// pointer to class for reading/writing data.
  ///
  class DataIOBase *dataio;

  ///
  /// pointer to class for reading/writing ascii-text-data.
  ///
  class DataIOBase *textio;

  /// function to setup data-I/O class.
  virtual void setup_dataio_class(
      class SimParams &,  ///< pointer to simulation parameters
      const int           ///< type of I/O: 1=text,2=fits,5=silo
  );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,      ///< reference to SimParams list.
      class GridBaseClass *,  ///< pointer to grid.
      const int               ///< level
  );

};  // setup_fixed_grid

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif  // if not SETUP_FIXED_GRID_H
