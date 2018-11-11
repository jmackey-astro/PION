
#ifndef RAYTRACER_SC_PLLEL_H
#define RAYTRACER_SC_PLLEL_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"
#include "boundaries/RT_MPI_boundaries.h"
#include "sim_params.h"
#include "decomposition/MCMD_control.h"

#include "raytracing/raytracer_base.h"
#include "raytracing/raytracer_SC.h"

#ifdef PARALLEL

///
/// Distributed Memory Raytracer, for a source which can be on grid, off grid,
///or at infinity in one coordinate direction.
///
class raytracer_USC_pllel : public raytracer_USC, public RT_MPI_bc {
 public:
  /// Constructor
  raytracer_USC_pllel(
      class GridBaseClass *, ///< Pointer to grid
      class microphysics_base *, ///< icroPhysics Class.
      class SimParams *,     ///< simulation parameters
      class MCMDcontrol *,   ///< domain decomposition info
      int,  ///< grid dimensionality
      int,  ///< coordinate system
      int,  ///< number of variables in state vector
      int,  ///< index of first tracer variable in state vector
      int   ///< Number of radiation sources
      );
  
  ~raytracer_USC_pllel(); ///< destructor.

  /// Adds a source to the list of sources to trace.
  ///
  /// Differs from serial version in that it also calls 
  /// gridptr->setup_RT_boundaries() so that the grid can decide which (if any)
  /// extra corner boundaries it needs to set up for sending and receiving data.
  ///
  virtual int Add_Source(
      struct rad_src_info * ///< ptr to source info.
      );

  /// Processes a source's effect on the grid over a timestep.
  /// This is a generic algorithm for 1D,2D,3D.  It calls different 
  /// functions depending on the dimensionality of the simulation.
  ///
  /// This parallel version receives boundary column densities, calls
  /// the serial raytracer, and then sends outgoing boundary column densities.
    
  virtual int RayTrace_SingleSource(
      const int,    ///< Source id
      const double, ///< Timestep
      const double  ///< eos gamma.
      );


  protected:
  /// Traces a 1D column from a starting cell, in a direction,
  /// to the edge of the grid. 
  int trace_column(
      const rad_source *,  ///< source we are tracing from.
      cell *,              ///< cell to start from.
      const enum direction ///< direction we are looking.
      );


  ///
  /// Short Characteristic Method of getting column density to cell. 
  /// The parallel version assumes boundary cells exist, so doesn't check that we
  /// are on the grid.
  ///
  virtual void col2cell_2d(
      const rad_source *,     ///< source we are working on.
      const cell *,           ///< cell to get column to.
      const enum direction,   ///< face ray enters cell through.
      const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
      const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
      double []               ///< Column densities.
      );

  ///\brief Short Characteristic Method of getting column density to cell.
  ///The parallel version assumes boundary cells exist, so doesn't check that we
  ///are on the grid.
  ///
  virtual void col2cell_3d(
      const rad_source *,     ///< source we are working on.
      const cell *,           ///< cell to get column to.
      const enum direction,   ///< face ray enters cell through.
      const enum direction *, ///< perp direction(s) towards source. (1 el array in 2d)
      const double *,         ///< fabs tan theta (angle(s) between 0 and 45deg) (1 el array in 2d)
      double []               ///< Column densities.
      );

  class SimParams *par; ///< pointer to simulation parameters
  class MCMDcontrol *MCMD; ///< pointer to domain decomposition

};
#endif //PARALLEL

// ##################################################################
// ##################################################################

#endif // RAYTRACER_SC_PLLEL_H

