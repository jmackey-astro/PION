/// \file dataio_silo_nestedgrid.h
/// \author Jonathan Mackey
/// \date   2018-05-04
/// 
///
/// Modifications:\n
/// - 2018.05.04 JM: new read/write functions for a nested grid.
/// 

#ifndef DATAIO_SILO_NESTEDGRID_H
#define DATAIO_SILO_NESTEDGRID_H


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#ifdef SILO

#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"
#include <silo.h>
#include <vector>




///
/// Class for reading and writing uniform grid data to silo files. 
///
class dataio_nested_silo : public dataio_silo {
  public:
  ///
  /// Constructor. 
  ///
  dataio_nested_silo(
      class SimParams &,  ///< pointer to simulation parameters
      std::string ///< read/write either FLOAT or DOUBLE to/from file
      );

  ///
  /// destructor
  ///
  virtual ~dataio_nested_silo();

  ///
  /// This writes the header for the simulation parameters,
  /// and then the data.
  /// 
  /// If the solver pointer is not null, it also writes data for the
  /// internal energy/Temperature (and the Magnetic Field divergence
  /// and total Pressure if the solver is an MHD solver).
  ///
  virtual int OutputData(
      const string, ///< File to write to
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &,  ///< pointer to simulation parameters
      const long int ///< number to stamp file with (e.g. timestep)
      );

  ///
  /// This reads the data variables in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the header, which should be read first.
  ///
  virtual int ReadData(
      string, ///< file to read from
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &  ///< pointer to simulation parameters
      );

  protected:

  ///
  /// Call once to setup arrays with the properties of the grid, for
  /// writing to the Silo Quadmesh.
  ///
  virtual int setup_grid_properties(
      class GridBaseClass *, ///< pointer to data.
      class SimParams &  ///< pointer to simulation parameters
      );
  
};



#endif // if SILO

#endif // if not DATAIO_SILO_NESTEDGRID_H

