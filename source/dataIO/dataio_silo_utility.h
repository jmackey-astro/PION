/// \file dataio_silo_utility.h
/// \author Jonathan Mackey
///
/// Utility class for reading silo files.
///
/// - 2010-02-02 JM: Added support for N procs to read data written
///     by M procs, where N not equal to M.
/// - 2010-02-03 JM: Fixed all the bugs in yesterday's work (as far
///     as i could find).
/// - 2013.02.19 JM: Moved dataio_utility's directory-listing
///    functions into file_status, so I deleted dataio_utility.
///    Renamed file from dataio_utility.h to dataio_silo_utility.h
/// - 2015.03.26 JM: updated for pion v0.2
/// - 2015.06.13/18 JM: tidied up code.

#ifndef DATAIO_UTILITY_H
#define DATAIO_UTILITY_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <dirent.h>
#include <errno.h>
#include <list>
using namespace std;

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_silo.h"
#include "dataIO/dataio_silo_MPI.h"
#include "dataIO/file_status.h"
#include "sub_domain/sub_domain.h"

#ifndef SILO
#error "Define SILO so this code will work!"
#endif

#ifndef PARALLEL
#error "Define PARALLEL so this code will work!"
#endif

class dataio_silo_utility : public dataio_silo_pllel {
public:
  dataio_silo_utility(
      class SimParams &,  ///< pointer to simulation parameters
      std::string,        ///< FLOAT or DOUBLE for files.
      class Sub_domain *);

  ~dataio_silo_utility() {}

  ///
  /// Parallel code to read any data.
  ///
  int ReadData(
      string,                           ///< file to read from
      vector<class GridBaseClass *> &,  ///< address of vector of grid pointers.
      class SimParams &                 ///< pointer to simulation parameters
  );

protected:
  ///
  /// Parallel code to read any data onto one level of a nested grid.
  ///
  int ReadLevelData(
      string,                 ///< file to read from
      class GridBaseClass *,  ///< address of vector of grid pointers.
      class SimParams &,      ///< pointer to simulation parameters
      const int               ///< level in grid hierarchy
  );

  ///
  /// Read data from a parallel file using a single processor.
  ///
  int serial_read_pllel_silodata(
      const string,           ///< filename
      class GridBaseClass *,  ///< pointer to data.
      const int,              ///< number of files
      const int,              ///< number of groups
      class SimParams &,      ///< pointer to simulation parameters
      class Sub_domain *      ///< pointer to class with nproc.
  );

  ///
  /// Given a Variable and a Silo File, read the data onto the grid.
  /// This can read a file that was written either in serial or parallel, with
  /// any number of processors.
  ///
  int SRAD_read_var2grid(
      DBfile *,               ///< pointer to silo file.
      class GridBaseClass *,  ///< pointer to data.
      const string,           ///< variable name to read.
      const long int,         ///< number of cells expected.
      class SimParams &,      ///< pointer to simulation parameters
      class Sub_domain *      ///< pointer to class with nproc.
  );

  ///
  /// get nproc and numfiles from the header, and returns an error
  /// if the Silo file is a serial file (i.e. if it fails to find numfiles).
  ///
  int SRAD_get_nproc_numfiles(
      const string,  ///< name of file to open
      int *,         ///< number of processes.
      int *          ///< number of files per timestep
  );

  ///
  /// Given a point, determine if it is on my local domain or not.
  /// Returns true if it is.
  ///
  bool SRAD_point_on_my_domain(
      const cell &,       ///< pointer to cell
      class SimParams &,  ///< pointer to simulation parameters
      class Sub_domain *  ///< pointer to class with nproc.
  );

  ///
  /// Parse filename and replace current filenumber
  /// (e.g. fname_0001.X.silo) with new one (e.g. fname_0002.X.silo).
  /// The integer argument to the function is the new filenumber.
  ///
  void set_pllel_filename(
      std::string &,  ///< filename
      const int       ///< file number
  );

  ///
  /// Get Xmin[],Xmax[] for a quadmesh, and return them.
  ///
  void get_quadmesh_extents(
      DBfile *,          ///< pointer to silo file.
      const string,      ///< directory of mesh
      const string,      ///< name of mesh
      double *,          ///< integer Xmin for mesh (output)
      double *,          ///< integer Xmax for mesh (output)
      class SimParams &  ///< pointer to simulation parameters
  );

  ///
  /// Get Xmin[],Xmax[] for a quadmesh, convert to integer positions,
  /// and return them.  Note that a grid must be set up for this to
  /// work!
  ///
  void get_quadmesh_integer_extents(
      DBfile *,                    ///< pointer to silo file.
      class GridBaseClass *,       ///< pointer to data.
      class SimParams &,           ///< pointer to simulation parameters
      const string,                ///< directory of mesh
      const string,                ///< name of mesh
      std::array<int, MAX_DIM> &,  ///< integer Xmin for mesh (output)
      std::array<int, MAX_DIM> &   ///< integer Xmax for mesh (output)
  );

  ///
  /// This function allows multi-core jobs to read data written by a
  /// single-core simulation.
  ///
  int parallel_read_serial_silodata(
      string,                 ///< file to read from
      class GridBaseClass *,  ///< pointer to data.
      class SimParams &       ///< pointer to simulation parameters
  );

  ///
  /// This function allows N-core jobs to read data written by M-core
  /// jobs, where M and N are not equal.
  ///
  int parallel_read_parallel_silodata(
      string,                   ///< file to read from
      class GridBaseClass *,    ///< pointer to data.
      class SimParams &,        ///< pointer to simulation parameters
      const int,                ///< number of files
      const int,                ///< number of groups
      const int,                ///< number of processes used to write file.
      const int,                ///< level in grid hierarchy
      const std::vector<int> &  ///< number of cells in each file
  );

  ///
  /// Given a Variable and a Silo File, read the data onto the grid.
  /// This can read a file that was written either in serial or parallel, with
  /// any number of processors.
  ///
  int PP_read_var2grid(
      DBfile *,               ///< pointer to silo file.
      class GridBaseClass *,  ///< pointer to data.
      class SimParams &,      ///< pointer to simulation parameters
      const string,           ///< variable name to read.
      const long int,         ///< number of cells expected (not needed)
      const int *,            ///< integer Xmin for mesh
      const int *             ///< integer Xmax for mesh
  );
};

#endif  // DATAIO_UTILITY_H
