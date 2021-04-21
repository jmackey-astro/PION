/// \file RT_MPI_boundaries.h
/// \brief Class declarations for raytracing across a parallel
///   decomposed domain, by communicating with neighbouring MPI
///   processes
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.09 JM: moved code from parallel uniform grid.

#ifndef RT_MPI_BOUNDARIES_H
#define RT_MPI_BOUNDARIES_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/boundaries.h"
#include "decomposition/MCMD_control.h"
#include "grid/grid_base_class.h"
#include "sim_params.h"

#include <vector>

// ##################################################################
// ##################################################################

enum rt_dirs {
  dir_NO = -1,
  dir_XN = 0,
  dir_XP = 1,
  dir_YN = 2,
  dir_YP = 3,
  dir_ZN = 4,
  dir_ZP = 5,
};

// ##################################################################
// ##################################################################

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

// ##################################################################
// ##################################################################

///
/// This is a boundary data struct.  Each boundary needs a source ID, a list of
/// boundaries to receive (which may be empty), and a list of boundaries to send
/// (which may also be empty).  Rather than have three separate lists which have
/// the same indexing, we just make a struct.
///
struct RT_source_comms_info {
  int source_id;  ///< id of the source.
  std::vector<struct RT_boundary_list_element>
      RT_recv_list,  ///< list of processors to receive data from, for each
                     ///< source.
      RT_send_list;  ///< list of processors to send data to, for each source.
};

// ##################################################################
// ##################################################################

///
/// Exchanging of column densities / optical depths across boundaries
/// between MPI processes, for radiative transfer.
///
class RT_MPI_bc {

public:
  ~RT_MPI_bc();

  ///
  /// Setup lists of processors to receive data from and send data to,
  /// and setup extra boundaries at corners.
  ///
  virtual int Setup_RT_Boundaries(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &,    ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      const int,              ///< source id
      struct rad_src_info &);

  ///
  /// Receive all optical depths for boundaries closer to source.
  ///
  virtual int Receive_RT_Boundaries(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &,    ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      const int,              ///< source id
      struct rad_src_info &);

  ///
  /// Send all optical depths for boundaries to domains further from
  /// source.
  ///
  virtual int Send_RT_Boundaries(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &,    ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      const int,              ///< source id
      struct rad_src_info &);

protected:
  ///
  /// This is the list where element i corresponds to source i, and the struct
  /// contains the list of boundaries we need to send and receive, which in
  /// turn contain the (ordered) list of cells which participate in the given
  /// boundary communication.
  ///
  std::vector<struct RT_source_comms_info> RT_source_list;

  ///
  /// If the source is for the purpose of calculating optical depth to diffuse
  /// radiation, then there is at most one send and one receive boundary
  /// element, and this function finds them and sets them up.
  ///
  int setup_RT_infinite_src_BD(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &,    ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      const int,              ///< Source id.
      struct rad_src_info &,
      std::vector<struct RT_boundary_list_element> &,  ///< RECV list for this
                                                       ///< source.
      std::vector<struct RT_boundary_list_element> &   ///< SEND list for this
                                                       ///< source.
  );

  ///
  /// If we have a source at infinity, this function returns the direction
  /// from the grid to the source.
  ///
  enum direction RT_src_at_infty_direction(
      class SimParams &,  ///< pointer to simulation parameters
      const int,          ///< source id.
      struct rad_src_info &);

  ///
  /// If the source is a monochromatic point source (not at infinity), then
  /// this function finds all the send and recv abutting domains and adds
  /// them, and a list of their constituent cells (in the correct order!) to
  /// the send and recv lists.
  ///
  int setup_RT_finite_ptsrc_BD(
      class SimParams &,      ///< pointer to simulation parameters
      class MCMDcontrol &,    ///< domain decomposition info
      class GridBaseClass *,  ///< pointer to grid.
      const int,              ///< Source id.
      struct rad_src_info &,  ///< SimParams list of radiation sources.
      std::vector<struct RT_boundary_list_element> &,  ///< RECV list for this
                                                       ///< source.
      std::vector<struct RT_boundary_list_element> &   ///< SEND list for this
                                                       ///< source.
  );

  ///
  /// allocate memory for new cells and attach them to the grid.
  ///
  int setup_RT_recv_boundary(
      class GridBaseClass *,             ///< pointer to grid.
      struct RT_boundary_list_element &  ///< pointer to boundary info.
  );

  ///
  /// find cells needed for send boundary, and add them to the list.
  ///
  int setup_RT_send_boundary(
      class GridBaseClass *,             ///< pointer to grid.
      struct RT_boundary_list_element &  ///< pointer to boundary info.
  );

  ///
  /// Add cells to the receive boundary list, so we know what to
  /// expect.
  ///
  int RT_populate_recv_boundary(
      struct boundary_data *,        ///< pointer to RT boundary data.
      const struct boundary_data *,  ///< pointer to BC boundary data.
      const enum direction           ///< face direction
  );
};

#endif  // RT_MPI_BOUNDARIES_H
