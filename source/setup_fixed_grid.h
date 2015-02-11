/// \file setup_fixed_grid.h
/// 
/// \brief Declares a class for setting up grids.
/// 
/// \author Jonathan Mackey
/// 
/// Modifications :\n
/// - 2015.02.09 JM: Split sim_control class into a setup class and
///   a derived class for running simulations.

#ifndef SETUP_FIXED_GRID_H
#define SETUP_FIXED_GRID_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "global.h"
#include "grid/grid_base_class.h"
#include "spatial_solvers/solver_eqn_base.h"
#include "grid/uniform_grid.h"
#include "MCMD_control.h"

///
/// The simplest finite volume grid - a uniform grid with cubic cells
/// in the chosen coordinates.  This class can set up the grid and
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
  void setup_cell_extra_data();

  /// \brief initialise the grid class with appropriate parameters.
  ///
  /// This function sets up the appropriate grid; so far I only have a 
  /// UniformGrid class -- uniform, cartesian, finite volume grid.
  ///
  virtual int setup_grid(
        class GridBaseClass **, ///< address of pointer to computational grid.
        class MCMDcontrol *     ///< address of MCMD controller class.
        );

  ///
  /// Determines what kind of boundary conditions are needed.
  /// Sets gp.Nbc to the appropriate value for the order of accuracy used.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int boundary_conditions(
        class GridBaseClass * 
        );   

  ///
  /// Decide if I need to setup MP class, and do it if i need to.
  ///
  virtual int setup_microphysics();

  ///
  /// Decide if I need to setup RT class, and do it if i need to.
  ///
  virtual int setup_raytracing(
        class GridBaseClass * 
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


}; // setup_fixed_grid
   
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#endif // if not SETUP_FIXED_GRID_H
