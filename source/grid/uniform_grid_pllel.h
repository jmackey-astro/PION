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


#ifdef PARALLEL

//
// integer flags for MPI communication labels.
//
#define BC_ANYtag 0 ///< works for either sort of communication.
#define BC_MPItag 1 ///< This is an integer tag on send/receive operations, to label that this communicates MPI boundary data.
#define BC_PERtag 2 ///< Integer tag to say it is for periodic BC.
#define BC_RTtag  3 ///< Integer tag to say we are transferring a radiative transfer column density tag.



enum rt_dirs {
  dir_NO     =-1,
  dir_XN     =0,
  dir_XP     =1,
  dir_YN     =2,
  dir_YP     =3,
  dir_ZN     =4,
  dir_ZP     =5,
};

///
/// Struct containing the information required for a boundary communication in
/// the ray-tracing algorithm.  RT_bd contains a list of cells which participate
/// in the communication (either the receiving or sending cells).
///
struct RT_boundary_list_element {
  int rank;
  int dir;
  struct boundary_data *RT_bd;
};

///
/// This is a boundary data struct.  Each boundary needs a source ID, a list of
/// boundaries to receive (which may be empty), and a list of boundaries to send
/// (which may also be empty).  Rather than have three separate lists which have
/// the same indexing, we just make a struct.
///
struct RT_source_comms_info {
  int source_id; ///< id of the source.
  std::vector<struct RT_boundary_list_element>
    RT_recv_list, ///< list of processors to receive data from, for each source.
    RT_send_list; ///< list of processors to send data to, for each source.
};


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
  protected:

  ///
  /// Assigns data to each boundary.  Called by SetupBCs().
  ///
  virtual int assign_boundary_data(
      const double,   ///< current simulation time (for DMACH)
      const double, ///< Simulation start time.
      const double,  ///< Simulation finish time.
      const double ///< minimum temperature allowed
      );
  
  ///
  /// Assigns data to a periodic boundary, getting data from another
  /// process if necessary.
  ///
  virtual int BC_assign_PERIODIC(  boundary_data *);
  
  ///
  /// Updates data on a periodic boundary, getting data from another
  /// process if necessary. 
  virtual int BC_update_PERIODIC(
        boundary_data *, ///< Boundary to update.
        const int,  ///< current fractional step being taken.
        const int   ///< final step.
        );
   
  ///
  /// Get boundary data from other processor.
  /// Int tag is to distinguish between periodic and internal boundaries, 
  /// as the domain can be split between two processors, so that proc 0 is
  /// getting a periodic and an internal boundary from proc 1, and vice
  /// versa.
  ///
  virtual int BC_assign_BCMPI(
        boundary_data *, ///< pointer to boundary we are assigning.
        int              ///< tag, either BC_MPItag or BC_PERtag.
        );

  ///
  /// Updates data on an inter-process communicating boundary.
  ///
  virtual int BC_update_BCMPI(
        boundary_data *, ///< Boundary to update.
        const int,  ///< current fractional step being taken.
        const int,   ///< final step.
        int              ///< tag, either BC_MPItag or BC_PERtag.
        );

  ///
  /// Set the boundary conditions string and initialise BC_bd
  ///
  virtual int BC_setBCtypes(
        class SimParams &  ///< reference to SimParams list.
        );

  int BC_select_data2send(
      list<cell *> *,  ///< list of cells (returned by this func.)
      int *,         ///< number of cells in list.
      boundary_data *  ///< pointer to boundary data.
      );

  ///
  /// This is the list where element i corresponds to source i, and the struct
  /// contains the list of boundaries we need to send and receive, which in
  /// turn contain the (ordered) list of cells which participate in the given
  /// boundary communication.
  ///
  std::vector<struct RT_source_comms_info> RT_source_list;
	
  ///
  /// If the source is for the purpose of calculating optical depth to diffuse
  /// radiation, then there is at most one send and one receive boundary element,
  /// and this function finds them and sets them up.
  ///
  int setup_RT_infinite_src_BD(
      const int, ///< Source id.
      struct rad_src_info &,
      std::vector<struct RT_boundary_list_element>  &, ///< RECV list for this source.
      std::vector<struct RT_boundary_list_element>  &  ///< SEND list for this source.
      );

  ///
  /// If we have a source at infinity, this function returns the direction
  /// from the grid to the source.
  ///
  enum direction RT_src_at_infty_direction(
      const int, ///< source id.
      struct rad_src_info &
      );

  ///
  /// If the source is a monochromatic point source (not at infinity), then this
  /// function finds all the send and recv abutting domains and adds them, and a
  /// list of their constituent cells (in the correct order!) to the send and 
  /// recv lists.
  ///
  int setup_RT_finite_ptsrc_BD(
      const int, ///< Source id.
      struct rad_src_info &, ///< SimParams list of radiation sources.
      std::vector<struct RT_boundary_list_element>  &, ///< RECV list for this source.
      std::vector<struct RT_boundary_list_element>  &  ///< SEND list for this source.
      );


  ///
  /// allocate memory for new cells and attach them to the grid.
  ///
  int setup_RT_recv_boundary(
      struct RT_boundary_list_element & ///< pointer to boundary info.
      );

  ///
  /// find cells needed for send boundary, and add them to the list.
  ///
  int setup_RT_send_boundary(
      struct RT_boundary_list_element & ///< pointer to boundary info.
      );

  ///
  /// Add cells to the receive boundary list, so we know what to
  /// expect. 
  ///
  int RT_populate_recv_boundary(
      struct boundary_data *, ///< pointer to RT boundary data.
      const struct boundary_data *, ///< pointer to BC boundary data.
      const enum direction ///< face direction
      );

  ///
  /// multi-core decomposition class.
  ///
  class MCMDcontrol *mpiPM;

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
      double *,  ///< array of max. x/y/z for full simulation.
      class MCMDcontrol * ///< pointer to MPI domain decomposition
      );

  /// 
  /// Deletes the grid.
  /// 
  ~UniformGridParallel();


  /// 
  /// Runs through ghost boundary cells and does the appropriate time update on them.
  /// 
  /// This is different from the serial version, as some boundaries need to get
  /// data from other processors.
  ///
  virtual int TimeUpdateExternalBCs(
      const double,   ///< current simulation time
      const int, ///< Current step number in tsimtime,mestep.
      const int  ///< Maximum step number in timestep.
      );

  ///
  /// Setup lists of processors to receive data from and send data to, 
  /// and setup extra boundaries at corners.
  ///
  virtual int Setup_RT_Boundaries(
      const int,  ///< source id
      struct rad_src_info &
      );

  ///
  /// Receive all optical depths for boundaries closer to source.
  ///
  virtual int Receive_RT_Boundaries(
      const int,  ///< source id
      struct rad_src_info &
      );

  ///
  /// Send all optical depths for boundaries to domains further from
  /// source.
  ///
  virtual int Send_RT_Boundaries(
      const int,  ///< source id
      struct rad_src_info &
      );

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
      double *,  ///< array of max. x/y/z for full simulation.
      class MCMDcontrol * ///< pointer to MPI domain decomposition
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
      double *,  ///< array of max. x/y/z for full simulation.
      class MCMDcontrol * ///< pointer to MPI domain decomposition
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
