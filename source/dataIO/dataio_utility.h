/// \file dataio_utility.h
/// \author Jonathan Mackey
/// 
/// Utility class for reading silo files.
///
///  - 2010-02-02 JM: Added support for N procs to read data written
///     by M procs, where N not equal to M.
///
///  - 2010-02-03 JM: Fixed all the bugs in yesterday's work (as far
///     as i could find).
///

#include <dirent.h>
#include <errno.h>
#include <list>
using namespace std;
#include "../global.h"
#include "dataio.h"
#include "dataio_silo.h"

#ifndef DATAIO_UTILITY_H
#define DATAIO_UTILITY_H

#ifndef SILO
#error "Define SILO so this code will work!"
#endif

#ifndef PARALLEL
#error "Define PARALLEL so this code will work!"
#endif

class dataio_utility {
 public:
  dataio_utility() {}
  virtual ~dataio_utility() {}
  //
  // Given a directory and a string to match files to (may be blank),
  // return sorted list of files.
  //
  int get_files_in_dir(const string,  ///< directory to list.
		       const string,  ///< string that files start with
		       list<string> * ///< list to put filenames in.
		       );
 protected:
  //
  // Return list of files in a given directory.
  //
  int get_dir_listing (const string,   ///< directory to list.
		       list<string> *  ///< list to put filenames in.
		       );
};

class dataio_silo_utility : public dataio_silo_pllel, public dataio_utility  {
 public:
  dataio_silo_utility() {}
  ~dataio_silo_utility() {}
  int serial_read_any_data(string, ///< file to read from
			   class GridBaseClass * ///< pointer to data.
			   );
  int parallel_read_any_data(string, ///< file to read from
			     class GridBaseClass * ///< pointer to data.
			     );
 protected:
  /** \brief Read data from a parallel file using a single processor. */
  int serial_read_pllel_silodata(const string, ///< filename
				 class GridBaseClass *, ///< pointer to data.
				 const int, ///< number of files
				 const int, ///< number of groups
				 class ParallelParams * ///< pointer to class with nproc.
				 );

  /** \brief Given a Variable and a Silo File, read the data onto the grid.
   * This can read a file that was written either in serial or parallel, with
   * any number of processors. */
  int SRAD_read_var2grid(DBfile *, ///< pointer to silo file.
			 class GridBaseClass *, ///< pointer to data.
			 const string,   ///< variable name to read.
			 const long int,  ///< number of cells expected.
			 class ParallelParams * ///< pointer to class with nproc.
			 );
  /** \brief get nproc and numfiles from the header, and returns an error
   * if the Silo file is a serial file (i.e. if it fails to find numfiles).
   */
  int SRAD_get_nproc_numfiles(const string, ///< name of file to open
			      int *,   ///< number of processes.
			      int *    ///< number of files per timestep
			      );
  /** \brief Given a point, determine if it is on my local domain or not.
   * Returns true if it is.
   */
  bool SRAD_point_on_my_domain(const cell *, ///< pointer to cell
			       class ParallelParams * ///< pointer to class with nproc.
			       );

  ///
  /// Parse filename and replace current filenumber
  /// (e.g. fname_0001.X.silo) with new one (e.g. fname_0002.X.silo).
  /// The integer argument to the function is the new filenumber.
  ///
  void set_pllel_filename(std::string &,  ///< filename
			  const int       ///< file number
			  );
  ///
  /// get the name of a mesh, currently "unigridXXXX" where XXXX=rank
  ///
  void mesh_name(const int, ///< rank
		 string &   ///< mesh name returned in this string.
		 );
  ///
  /// Given myrank and mygroup, create directory string where data is
  /// held.  This is determined as rank_[global_rank]_domain_[rank_within_group]
  ///
  void set_dir_in_file(std::string &,  ///< directory name.
		       const int,       ///< myrank (global).
		       const int        ///< myrank in group.
		       );
  ///
  /// Get Xmin[],Xmax[] for a quadmesh, and return them.
  ///
  void get_quadmesh_extents(DBfile *,     ///< pointer to silo file.
			    const string, ///< directory of mesh
			    const string, ///< name of mesh
			    double *, ///< integer Xmin for mesh (output)
			    double *  ///< integer Xmax for mesh (output)
			    );
  ///
  /// Get Xmin[],Xmax[] for a quadmesh, convert to integer positions,
  /// and return them.  Note that a grid must be set up for this to
  /// work!
  ///
  void get_quadmesh_integer_extents(DBfile *,        ///< pointer to silo file.
				    const string, ///< directory of mesh
				    const string, ///< name of mesh
				    int *, ///< integer Xmin for mesh (output)
				    int *  ///< integer Xmax for mesh (output)
				    );
  ///
  /// This function allows multi-core jobs to read data written by a
  /// single-core simulation.
  ///
  int parallel_read_serial_silodata(string,      ///< file to read from
				    class GridBaseClass * ///< pointer to data.
				    );
  ///
  /// This function allows N-core jobs to read data written by M-core
  /// jobs, where M and N are not equal.
  ///
  int parallel_read_parallel_silodata(string,    ///< file to read from
				      class GridBaseClass *, ///< pointer to data.
				      const int, ///< number of files
				      const int, ///< number of groups
				      const int  ///< number of processes used to write file.
				      );
  ///
  /// Given a Variable and a Silo File, read the data onto the grid.
  /// This can read a file that was written either in serial or parallel, with
  /// any number of processors.
  ///
  int PP_read_var2grid(DBfile *, ///< pointer to silo file.
		       class GridBaseClass *, ///< pointer to data.
		       const string,   ///< variable name to read.
		       const long int,  ///< number of cells expected (not needed)
		       const int *, ///< integer Xmin for mesh
		       const int *  ///< integer Xmax for mesh
		       );
};


#endif // DATAIO_UTILITY_H
