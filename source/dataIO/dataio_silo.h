/// \file silo.h
/// This file contains class declarations for the serial and parallel
/// versions of the Silo dataio Class.  Info on Silo can be found at
/// https://wci.llnl.gov/codes/silo/ on LLNL website.  It's advantage
/// is that it has good parallel I/O interface that lets VisIt read in 
/// all the files as though the data were just in one file.
/// 
/// - 2010-02-04 JM: Added multimesh-adjacency object write so that I
///    can get streamlines to cross boundaries when plotting with Visit.
/// - 2010-04-11 JM: parallel class gets its own
/// setup_write_variables() class so that it can save disk space by
/// only writing primitive variables. (tidied up comments too).
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010-04-25 JM: renamed parallel choose_filename to choose_pllel_filename()
/// - 2010-07-20/22 JM: Work on new dataio structure with a list of parameters
///    to read and write.
/// - 2010.07.23 JM: removed obselete read_header(),
///    write_header() functions.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2011.03.22 JM: Removed parallel setup_write_variables() function.  There
///    was too much duplicate code with the serial version.
/// - 2011.06.02 JM: Added WriteHeader() function so I can over-write header
///    parameters and restart a simulation with e.g. different microphysics.
/// - 2015.01.28 JM: Tidied up a lot, added new mpiPM pointer.

#ifndef DATAIO_SILO_H
#define DATAIO_SILO_H

#ifdef SILO

#include "dataIO/dataio.h"
#include "MCMD_control.h"
#include <silo.h>
#include <vector>

//
// These are here because I am initially couldn't write double
// precision data to the silo files, but I figured it out.
//
//#define FAKE_DOUBLE float //#define FAKE_DATATYPE DB_FLOAT
#define FAKE_DOUBLE double
#define FAKE_DATATYPE DB_DOUBLE

///
/// write vectors (velocity,B) in component form as scalar variables
///
#define SILO_SCALARS
///
/// write vectors (velocity,B) as vectors (needs more memory! and no
/// advantage)
///
//#define SILO_VECTORS

// one of the following must be defined -- which driver to use (HDF or Silo's PDB)
//#define SILO_FILETYPE DB_HDF5
#define SILO_FILETYPE DB_PDB

///
/// Class for reading and writing uniform grid data to silo files. 
///
class dataio_silo :public DataIOBase {
  public:
  ///
  /// Constructor. 
  ///
  dataio_silo();

  ///
  /// sets equation pointer to zero.
  ///
  virtual ~dataio_silo();

  ///
  /// Class can run with or without a solver, but this function
  /// allows you to set a pointer to the solver.#
  ///
  void SetSolver(FV_solver_base * ///< Pointer to the solver (for Eint,divB,Ptot)
		  );
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
        class GridBaseClass *, ///< pointer to data.
        const long int ///< number to stamp file with (e.g. timestep)
        );
  ///
  /// Reads the header from the file specified, and puts the
  /// simulation parameters in the global data class declared in global.h.
  /// This should be called first if we want to set up a uniform grid with
  /// the right dimensions to take the data from the file.
  ///
   virtual int ReadHeader(string ///< file to read from
			  );


  ///
  /// Write simulation header info to file (file must already exist!)
  ///
  virtual int WriteHeader(
          const string ///< file to write to (full, exact filename).
          );


  ///
  /// This reads the data variables in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the header, which should be read first.
  ///
   virtual int ReadData(string, ///< file to read from
			class GridBaseClass * ///< pointer to data.
			);
  protected:
   class GridBaseClass *gp; ///< pointer to computational grid.
   class FV_solver_base *eqn; ///< pointer to the solver, which knows the equations we are solving.
   string silofile;         ///< filename to write data to.
   int silo_filetype;       ///< What sort of file to write, can be DB_PDB or DB_HDF5, among others.
   int ndim;                ///< Dimensionality of Grid.
   int vec_length;          ///< Length of Vectors (always 3d for now).
   unsigned int strlength;  ///< length of character arrays in class.
   FAKE_DOUBLE **node_coords; ///< Ndim array of coordinates of nodes of uniform grid (zone corners!).
   FAKE_DOUBLE 
     *nodex,        ///< array of x node positions (addressed through node_coords).
     *nodey,        ///< array of y node positions (addressed through node_coords).
     *nodez;        ///< array of z node positions (addressed through node_coords).
   int *nodedims;   ///< ndim element array with number of nodes along each direction.
   int *zonedims;   ///< ndim element array with number of zones/cells along each direction.
   bool have_setup_gridinfo;     ///< false initially, set to true on first file write.
   bool have_setup_writevars;    ///< false initially, set to true on first file write.
   std::vector<string> varnames; ///< list of names of variables to write to file.
   std::vector<string> readvars; ///< list of variable names to read from file when starting.
   FAKE_DOUBLE 
     *data0,     ///< array for grid data to write to the mesh.
     *data1,     ///< array for grid data to write to the mesh.
     *data2,     ///< array for grid data to write to the mesh.
     **vec_data; ///< array of pointers to data for vector data.

   DBoptlist 
     *GridOpts; ///< List of options relating to mesh.

   DBfile
     **db_ptr; ///< pointer to Silo file pointer.

  ///
  /// Choose a FileName to write to.
  ///
   virtual int choose_filename(const string, ///< filebase passed in from main code.
			       const int     ///< file counter to use (e.g. timestep).
			       );
  ///
  /// Call once to setup arrays with the properties of the grid, for
  /// writing to the Silo Quadmesh.
  ///
   virtual int setup_grid_properties(
        class GridBaseClass * ///< pointer to data.
        );
  ///
  /// Choose what data to write to file, based on equations being solved.*/
  ///
   virtual int setup_write_variables();

  ///
  /// Generate Quadmesh and write it to silo file.
  ///
   int generate_quadmesh(DBfile *, ///< pointer to silo file to write to.
			 string    ///< name of mesh to write.
			 );
  ///
  /// Given a file, meshname, variable, write that variable to the mesh.*/
  ///
   int write_variable2mesh(DBfile *, ///< pointer to silo file.
			   string,   ///< name of mesh to write to.
			   string    ///< variable name to write.
			   );
  ///
  /// allocate memory for arrays used to write data to files.*/
  ///
   virtual void create_data_arrays();
  ///
  /// deallocate memory for arrays used to write data to files.*/
  ///
   void delete_data_arrays();
  ///
  /// Get Scalar variable into array ready to write to file.
  ///
   int get_scalar_data_array(string, ///< variable name to get.
			     FAKE_DOUBLE * ///< array to write to.
			     );
  ///
  /// Get Vector variable into 2d array ready to write to file.
  ///
   int get_vector_data_array(string, ///< variable name to get.
			     FAKE_DOUBLE ** ///< array to write to.
			     );
  ///
  /// Writes a scalar array to the specified mesh and file.
  ///
   int write_scalar2mesh(DBfile *, ///< silo file pointer.
			 string,   ///< mesh name
			 string,   ///< variable name
			 FAKE_DOUBLE *  ///< pointer to data array.
			 );
  ///
  /// Writes a vector array to the specified mesh and file.
  ///
   int write_vector2mesh(DBfile *, ///< silo file pointer.
			 string,   ///< mesh name
			 string,   ///< variable name
			 FAKE_DOUBLE **  ///< pointer to data array.
			 );

   ///
   /// Function which defines how to get the data from a silo file.
   ///
   int read_header_param(class pm_base *);
   ///
   /// Function which defines how to write a parameter to a silo file.
   ///
   int write_header_param(class pm_base *);

  ///
  /// Set list of variable names to read, based on equations to solve.
  ///
   int set_readvars(int ///< SimPM.eqntype -- the type of equations we are solving.
		    );
  ///
  /// Given a Variable and a Silo File, read the data onto the grid.
  ///
   int read_variable2grid(DBfile *, ///< pointer to silo file.
			  string,   ///< name of mesh to read from.
			  string,   ///< variable name to read.
			  long int  ///< number of cells expected.
			  );
};



#ifdef PARALLEL
///
/// parallel I/O using the PMPIO interface, based on serial class.
///
class dataio_silo_pllel : public dataio_silo {
 public:
  ///
  /// Constructor.
  ///
  dataio_silo_pllel(
      class MCMDcontrol *  ///< address of MCMD controller class.
  );

  /// Destructor (doensn't have much to do).
  virtual ~dataio_silo_pllel();

  ///
  /// This writes the header and data for the simulation parameters.
  /// 
  /// If the solver pointer is not null, it also writes some derived variables
  /// such as Temperature, Div(B), etc.
  ///
  int OutputData(const string, ///< File to write to
		 class GridBaseClass *, ///< pointer to data.
		 const long int ///< number to stamp file with (e.g. timestep)
		 );

  ///
  /// Reads the header from the file specified, and puts the
  /// simulation parameters in the global data class declared in global.h.
  /// This should be called first if we want to set up a uniform grid with
  /// the right dimensions to take the data from the file later.  This
  /// parallel version uses PMPIO to let the processes take turns to read 
  /// the header from the root file, regardless of the number of files 
  /// present.
  ///
  int ReadHeader(string ///< file to read from
		 );

  ///
  /// This reads the data variables in turn, and puts the data into
  /// the grid points, assuming the grid has been set up with paramters
  /// from the header, which should be read first.  It also assumes the
  /// domain decomposition has been done already, and will bug out if it
  /// hasn't.  The routine uses the PMPIO interface.
  ///
  int ReadData(string, ///< file to read from
	       class GridBaseClass * ///< pointer to data.
	       );

 protected:
  class MCMDcontrol *mpiPM;
  ///
  /// Choose filename based on outfile, group_rank, and counter.
  ///
  virtual int choose_pllel_filename(
        const string, ///< filebase passed in from main code.
        const int,    ///< group_rank (i.e. which file I write to)
        const int,    ///< file counter to use (e.g. timestep).
        string &      ///< string to return filename in.
        );

  ///
  /// Call once to setup arrays with the properties of the grid, for
  /// writing to the Silo Quadmesh.  This is for the local domain of the current
  /// processor.
  ///
  int setup_grid_properties();

  ///
  /// allocate memory for arrays used to write data to files.
  /// Arrays are different sizes for parallel grid -- local Ncell instead
  /// of global Ncell.
  ///
  void create_data_arrays();
  
  int numfiles; ///< number of files to split the data into.

  ///
  /// Write a mulitmesh adjacency object (obsolete and doesn't work!)
  /// 
  int write_multimeshadj(
        DBfile *, ///< pointer to silo file.
        class GridBaseClass *, ///< pointer to data.
        string, ///< multimesh name
        string  ///< multimesh adjacency name.
        );

  ///
  /// Write an MRG tree object (replaces multimesh adjacency)
  /// 
  int write_MRGtree(
        DBfile *, ///< pointer to silo file.
        class GridBaseClass *, ///< pointer to data.
        string, ///< multimesh name
        string  ///< MRG tree name.
        );

  ///
  /// Write a mulitmesh object
  /// 
  int write_multimesh(
        DBfile *, ///< pointer to silo file.
        class GridBaseClass *, ///< pointer to data.
        string, ///< multimesh name
        string  ///< multimesh adjacency name.
        );
};

#endif // if PARALLEL

#endif // if SILO

#endif // if not DATAIO_SILO_H

