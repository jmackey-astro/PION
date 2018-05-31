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

#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "grid/uniform_grid.h"
#include "decomposition/MCMD_control.h"

#define TIMESTEP_FULL 1
#define TIMESTEP_FIRST_PART 2

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
    BCMPI      =10, ///< boundary between processor domains in parallel grid.
    RADSHOCK   =11, ///< Boundary condition adjacent to cold wall for radiative shock test problem.
    RADSH2     =12, ///< Outflow augmented with fixed outflow speed.
    ONEWAY_OUT =13, ///< One-way valve -- allows outflow but not inflow (zero gradient OR reflecting).
    STWIND     =14, ///< Stellar wind sources exist, so apply internal boundaries.
    STARBENCH1 =15, ///< StarBench test for mixing with a solid wall.
    FINE_TO_COARSE  =16, ///< data in this cell should be obtained from finer-scale data on a nested grid.
    COARSE_TO_FINE   =17  ///< data should be obtained from coarser level in a nested grid.
};


///
/// The simplest finite volume grid - a uniform grid with cells that
/// are cubes in the chosen coordinates.  This class sets up the grid and
/// other things to get a simulation ready to run, so it is useful
/// for simulation analysis.  PION itself uses a derived class to
/// setup and run simulations.
///
class setup_fixed_grid
{
  public:
  setup_fixed_grid();  ///< Simple constructor, initialises value.
  virtual ~setup_fixed_grid(); ///< Deletes any dynamic memory, if not already done.


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
      class GridBaseClass **, ///< address of pointer to grid.
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol *     ///< address of MCMD controller class.
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
      class SimParams &,  ///< pointer to simulation parameters
      class GridBaseClass * ///< pointer to computational grid
      );

  ///
  /// Check for any time-evolving radiation sources, and read the evolution
  /// file if there are any.  Data is stored in global struct SimPM.STAR[v]
  ///
  virtual int setup_evolving_RT_sources(
      class SimParams &,  ///< pointer to simulation parameters
      class RayTracingBase * ///< pointer to raytracing class
      );

  ///
  /// Determines what kind of boundary conditions are needed and
  /// creates the boundary data structures.  Asks the grid to create
  /// grid cells for the external boundaries, and label internal
  /// boundary cells as such.
  ///
  int boundary_conditions(
      class SimParams &,  ///< pointer to simulation parameters
      class GridBaseClass *  ///< pointer to grid.
      );

  ///
  /// Assigns data to each boundary.
  ///
  virtual int assign_boundary_data(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *  ///< pointer to grid.
      );

  //---------------------------------------
  protected:
  //---------------------------------------

  /// flag: true if timestep limit needs raytracing column densities
  bool FVI_need_column_densities_4dt;

  int FVI_nheat; ///< number of RT heating sources

  int FVI_nion;  ///< number of ionising sources

  /// vector of RT heating sources, of size FVI_nheat
  std::vector<struct rt_source_data> FVI_heating_srcs;

  /// vector of RT ionising sources, of size FVI_nion
  std::vector<struct rt_source_data> FVI_ionising_srcs;

  ///
  /// Check for any time-evolving radiation sources, and update source
  /// properties from global struct SimPM.STAR[v] if needed
  ///
  virtual int update_evolving_RT_sources(
      class SimParams &,  ///< pointer to simulation parameters
      class RayTracingBase * ///< pointer to raytracing class [OUTPUT]
      );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int setup_boundary_structs(
      class SimParams &,  ///< reference to SimParams list.
      class GridBaseClass *  ///< pointer to grid.
      );

  /// Assigns data to a periodic boundary.
  virtual int BC_assign_PERIODIC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on an outflow (zero gradient) boundary.
  int BC_assign_OUTFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  ///
  /// Assigns data on a one-way outflow (zero gradient) boundary.
  /// If the flow is off-domain, then I use zero-gradient, but if flow
  /// is onto domain then I set the boundary cells to have zero normal 
  /// velocity.
  ///
  int BC_assign_ONEWAY_OUT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on an inflow (fixed) boundary.
  int BC_assign_INFLOW(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on a reflecting boundary.
  int BC_assign_REFLECTING(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on a fixed boundary.
  int BC_assign_FIXED(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  ///
  /// Sets some boundary cells to be fixed to the Jet inflow
  /// condition.
  ///
  virtual int BC_assign_JETBC(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  ///
  /// Assigns data on a JetReflect boundary, which is the same
  /// as a reflecting boundary, except the B-field is reflected differently.
  /// It is the base of a jet, so there is assumed to be another jet
  /// in the opposite direction, so the normal B-field is unchanged across
  /// the boundary and the tangential field is reversed.
  ///
  int BC_assign_JETREFLECT(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on a boundary for the Double Mach Reflection Problem.
  int BC_assign_DMACH(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  /// Assigns data on The other DMR test problem boundary
  virtual int BC_assign_DMACH2(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  ///
  /// Add internal stellar wind boundaries -- these are (possibly time-varying)
  /// winds defined by a mass-loss-rate and a terminal velocity.  A region within
  /// the domain is given fixed values corresponding to a freely expanding wind
  /// from a cell-vertex-located source.
  ///
  int BC_assign_STWIND(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      boundary_data *
      );

  ///
  /// Add cells to both the Wind class, and to the boundary data list
  /// of cells.  This is re-defined for cylindrical and spherical
  /// coords below.
  ///
  virtual int BC_assign_STWIND_add_cells2src(
      class SimParams &,      ///< pointer to simulation parameters
      class GridBaseClass *,  ///< pointer to grid.
      //boundary_data *, ///< boundary ptr.
      const int    ///< source id
      );

}; // setup_fixed_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_FIXED_GRID_H
