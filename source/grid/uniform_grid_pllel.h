/// \file uniform_grid_pllel.h
///
/// \brief Declares MPI-parallelised version of the uniform grid.
///
/// \author Jonathan Mackey
///


#ifndef UNIFORM_GRID_PLLEL_H
#define UNIFORM_GRID_PLLEL_H



#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <list>
using namespace std;

#include "sim_constants.h"
#include "sim_params.h"


#include "grid/grid_base_class.h"
#include "grid/uniform_grid.h"
#include "grid/stellar_wind_BC.h"
#include "grid/stellar_wind_angle.h"
#include "coord_sys/VectorOps.h"
#include "coord_sys/VectorOps_spherical.h"
#include "decomposition/MCMD_control.h"


#ifdef PARALLEL


///
/// Parallel implementation of the serial uniform grid.
/// 
/// This differs mostly in that it has to treat the boundaries differently.
/// There are new internal boundaries between processes, and periodic 
/// boundaries may or may not need to get data from a different process to 
/// update themselves.
///
class UniformGridParallel
: virtual public UniformGrid
{

  public:
  /// 
  /// Constructor. Sets up a grid in the same way as the serial grid.
  /// 
  UniformGridParallel(
      int,         ///< ndim
      int,         ///< nvar
      int,         ///< equation type
      int, ///< number of boundary cells to use.
      double *,    ///< local xmin
      double *,    ///< local xmax
      int *,       ///< local number of grid zones
      double *, ///< array of min. x/y/z for full simulation.
      double *  ///< array of max. x/y/z for full simulation.
      );

  /// 
  /// Deletes the grid.
  /// 
  ~UniformGridParallel() {return;}

  /// Returns Simulation xyz lower bounds (code units)
  virtual double SIM_Xmin(enum axes a) const
  {return(Sim_xmin[a] );}

  /// Returns Simulation xyz upper bounds (code units)
  virtual double SIM_Xmax(enum axes a) const
  {return(Sim_xmax[a] );}

  /// Returns Simulation range (code units)
  virtual double SIM_Range(enum axes a) const
  {return(Sim_range[a]);}

  /// Returns Simulation xyz lower bounds (integer units, 1cell=2units)
  virtual int  Sim_iXmin(enum axes a) const
  {return(Sim_ixmin[a] );}

  /// Returns Simulation xyz upper bounds (integer units, 1cell=2units)
  virtual int  Sim_iXmax(enum axes a) const
  {return(Sim_ixmax[a] );}

  /// Returns Simulation xyz range (integer units, 1cell=2units)
  virtual int Sim_iRange(enum axes a) const
  {return(Sim_irange[a]);}

  ///
  /// Setup the flux struct flux_update_recv with list of interfaces
  /// that need to be updated with fluxes from a finer level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).  This is an MPI-parallelised
  /// version which can deal with grids not on this process.
  ///
  int setup_flux_recv(
      class SimParams &,  ///< simulation params (including BCs)
      const int           ///< level to receive from
      );

  ///
  /// Setup the flux struct flux_update_send with list of interfaces
  /// that need to be sent to a coarser level grid.
  /// These fluxes are used to correct the fluxes on the coarse grid,
  /// to ensure that they are consistent across all levels, following
  /// Berger & Colella (1989,JCP,82,64).  This is an MPI-parallelised
  /// version which can deal with grids not on this process.
  ///
  int setup_flux_send(
      class SimParams &,  ///< simulation params (including BCs)
      const int           ///< level to receive from
      );


  protected:

};

///
/// Uniform Grid in cylindrical coordinates, for parallel simulations
/// (i.e. each parallel grid is a part of the overall simulation
/// domain).  The cylindrical version takes account of the fact that
/// the cell radial coordinate is not the midpoint of the cell because
/// of radial divergence.  It is the centre-of-volume, which is at a
/// larger radius than the midpoint (although the difference becomes
/// negligible at large radii it has a significant effect at small
/// radii).
///
class uniform_grid_cyl_parallel
: virtual public UniformGridParallel,
  virtual public uniform_grid_cyl
{
 public:
  ///
  /// The constructor won't do very much:
  ///
  uniform_grid_cyl_parallel(
      int, ///< ndim, length of position vector.
      int, ///< nvar, length of state vectors.
      int, ///< eqntype, which equations we are using (needed by BCs).
      int, ///< number of boundary cells to use.
      double *, ///< array of minimum values of x,y,z.
      double *, ///< array of maximum values of x,y,z.
      int *, ///< array of number of cells in x,y,z directions.
      double *, ///< array of min. x/y/z for full simulation.
      double *  ///< array of max. x/y/z for full simulation.
      );

  ///
  /// Nor will the destructor
  ///
  ~uniform_grid_cyl_parallel();

  ///
  /// Returns the centre of volume of a cell (in the radial
  /// direction) in the dimensionless integer coordinate system.
  /// It is redefined here because we need the radius calculated from
  /// the global simulation Xmin[Rcyl], not the grid Xmin.
  ///
  virtual double iR_cov(const cell *);
};

///
/// Uniform Grid in spherical coordinates, for parallel simulations
/// (i.e. each parallel grid is a part of the overall simulation
/// domain).  The spherical version takes account of the fact that
/// the cell radial coordinate is not the midpoint of the cell because
/// of radial divergence.  It is the centre-of-volume, which is at a
/// larger radius than the midpoint (although the difference becomes
/// negligible at large radii it has a significant effect at small
/// radii).
///
class uniform_grid_sph_parallel
: virtual public UniformGridParallel,
  virtual public uniform_grid_sph
{
 public:
  ///
  /// The constructor won't do very much:
  ///
  uniform_grid_sph_parallel(
      int, ///< ndim, length of position vector.
      int, ///< nvar, length of state vectors.
      int, ///< eqntype, which equations we are using (needed by BCs).
      int, ///< number of boundary cells to use.
      double *, ///< array of minimum values of x,y,z.
      double *, ///< array of maximum values of x,y,z.
      int *, ///< array of number of cells in x,y,z directions.
      double *, ///< array of min. x/y/z for full simulation.
      double *  ///< array of max. x/y/z for full simulation.
      );

  ///
  /// Nor will the destructor
  ///
  ~uniform_grid_sph_parallel();

  ///
  /// Returns the centre of volume of a cell (in the radial
  /// direction) in the dimensionless integer coordinate system.
  /// It is re-defined here because we need the radius calculated from
  /// the global simulation Xmin[Rcyl], not the grid Xmin.
  ///
  virtual double iR_cov(const cell *);
};


#endif // PARALLEL

#endif // UNIFORM_GRID_PLLEL_H
