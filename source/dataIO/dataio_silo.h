/// \file dataio_silo.h
///
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
/// - 2015.06.13 JM: Changed datatype (FLOAT/DOUBLE) to a runtime
///    parameter, set in the constructor.
/// - 2016.03.18 JM: better treatment of float/double coordinates.
/// - 2018.01.24 JM: worked on making SimPM non-global

#ifndef DATAIO_SILO_H
#define DATAIO_SILO_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SILO

#include "dataIO/dataio_base.h"
#include "decomposition/MCMD_control.h"
#include <silo.h>
#include <vector>

//#define WRITE_GHOST_ZONES
//#define TEST_SILO_IO

///
/// write vectors (velocity,B) in component form as scalar variables
///
#define SILO_SCALARS
///
/// write vectors (velocity,B) as vectors (needs more memory! and no
/// advantage)
///
//#define SILO_VECTORS

// one of the following must be defined -- which driver to use (HDF or Silo's
// PDB)
//#define SILO_FILETYPE DB_HDF5
#define SILO_FILETYPE DB_PDB

///
/// Class for reading and writing uniform grid data to silo files.
///
class dataio_silo : public DataIOBase {
  public:
    ///
    /// Constructor.
    ///
    dataio_silo(
        class SimParams&,  ///< pointer to simulation parameters
        std::string        ///< read/write either FLOAT or DOUBLE to/from file
    );

    ///
    /// sets equation pointer to zero.
    ///
    virtual ~dataio_silo();

    ///
    /// Class can run with or without a solver, but this function
    /// allows you to set a pointer to the solver.#
    ///
    void SetSolver(
        FV_solver_base*  ///< Pointer to the solver (for Eint,divB,Ptot)
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
        const string,                   ///< File to write to
        vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
        class SimParams&,               ///< pointer to simulation parameters
        const long int  ///< number to stamp file with (e.g. timestep)
    );

    ///
    /// Reads the header from the file specified, and puts the
    /// simulation parameters in the global data class declared in global.h.
    /// This should be called first if we want to set up a uniform grid with
    /// the right dimensions to take the data from the file.
    ///
    virtual int ReadHeader(
        string,           ///< file to read from
        class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// Write simulation header info to file (file must already exist!)
    ///
    virtual int WriteHeader(
        const string,     ///< file to write to (full, exact filename).
        class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// This reads the data variables in turn, and puts the data into
    /// the grid points, assuming the grid has been set up with paramters
    /// from the header, which should be read first.
    ///
    virtual int ReadData(
        string,                         ///< file to read from
        vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
        class SimParams&                ///< pointer to simulation parameters
    );

  protected:
    class GridBaseClass* gp;  ///< pointer to computational grid.

    class FV_solver_base* eqn;  ///< pointer to the solver, which knows the
                                ///< equations we are solving.

    ///
    /// Whether to write float or double data (DB_FLOAT or DB_DOUBLE)
    ///
    int silo_dtype;
    ///
    /// what the coordinate system is DB_CARTESIAN, etc.
    ///
    int silo_coordsys;

    string silofile;         ///< filename to write data to.
    int silo_filetype;       ///< What sort of file to write, can be DB_PDB or
                             ///< DB_HDF5, among others.
    int ndim;                ///< Dimensionality of Grid.
    int vec_length;          ///< Length of Vectors (always 3d for now).
    unsigned int strlength;  ///< length of character arrays in class.
    int* nodedims;  ///< ndim element array with number of nodes along each
                    ///< direction.
    int* zonedims;  ///< ndim element array with number of zones/cells along
                    ///< each direction.
    bool have_setup_gridinfo;   ///< false initially, set to true on first file
                                ///< write.
    bool have_setup_writevars;  ///< false initially, set to true on first file
                                ///< write.
    std::vector<string>
        varnames;  ///< list of names of variables to write to file.
    std::vector<string>
        readvars;  ///< list of variable names to read from file when starting.

    void *data0,     ///< array for grid data to write to the mesh.
        *data1,      ///< array for grid data to write to the mesh.
        *data2,      ///< array for grid data to write to the mesh.
        *mask,       ///< array for mask for nested grid.
        **vec_data;  ///< array of pointers to data for vector data.

    ///
    /// Ndim array of coordinates of nodes of uniform grid (zone corners!).
    ///
    void** node_coords;
    void
        *nodex,  ///< array of x node positions (addressed through node_coords).
        *nodey,  ///< array of y node positions (addressed through node_coords).
        *nodez;  ///< array of z node positions (addressed through node_coords).

    DBoptlist* GridOpts;  ///< List of options relating to mesh.

    DBfile** db_ptr;  ///< pointer to Silo file pointer.

    ///
    /// Choose a FileName to write to.
    ///
    virtual int choose_filename(
        const string,  ///< filebase passed in from main code.
        const int      ///< file counter to use (e.g. timestep).
    );

    ///
    /// Call once to setup arrays with the properties of the grid, for
    /// writing to the Silo Quadmesh.
    ///
    virtual int setup_grid_properties(
        class GridBaseClass*,  ///< pointer to data.
        class SimParams&       ///< pointer to simulation parameters
    );

    ///
    /// Choose what data to write to file, based on equations being solved.*/
    ///
    virtual int setup_write_variables(
        class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// Generate Quadmesh and write it to silo file.
    ///
    int generate_quadmesh(
        DBfile*,          ///< pointer to silo file to write to.
        string,           ///< name of mesh to write.
        class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// Given a file, meshname, variable, write that variable to the mesh.*/
    ///
    int write_variable2mesh(
        class SimParams& SimPM,  ///< pointer to simulation parameters
        DBfile*,                 ///< pointer to silo file.
        string,                  ///< name of mesh to write to.
        string                   ///< variable name to write.
    );

    ///
    /// allocate memory for arrays used to write data to files.*/
    ///
    virtual void create_data_arrays(
        class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// deallocate memory for arrays used to write data to files.*/
    ///
    void delete_data_arrays();

    ///
    /// Get Scalar variable into integer array ready to write to file.
    ///
    int get_int_scalar_data_array(
        string,            ///< variable name to get.
        class SimParams&,  ///< pointer to simulation parameters
        void*              ///< array to write to
    );

    ///
    /// Get Scalar variable into array ready to write to file.
    ///
    int get_scalar_data_array(
        string,            ///< variable name to get.
        class SimParams&,  ///< pointer to simulation parameters
        void*              ///< array to write to.
    );

    ///
    /// Get Vector variable into 2d array ready to write to file.
    ///
    int get_vector_data_array(
        string,            ///< variable name to get.
        class SimParams&,  ///< pointer to simulation parameters
        void**             ///< array to write to.
    );

    ///
    /// Writes a scalar array to the specified mesh and file.
    ///
    int write_scalar2mesh(
        DBfile*,  ///< silo file pointer.
        string,   ///< mesh name
        string,   ///< variable name
        void*     ///< pointer to data array.
    );

    ///
    /// Writes a vector array to the specified mesh and file.
    ///
    int write_vector2mesh(
        DBfile*,  ///< silo file pointer.
        string,   ///< mesh name
        string,   ///< variable name
        void**    ///< pointer to data array.
    );

    ///
    /// Function which defines how to get the data from a silo file.
    ///
    int read_header_param(class pm_base*);

    ///
    /// Function which defines how to write a parameter to a silo file.
    ///
    int write_header_param(class pm_base*);

    ///
    /// Set list of variable names to read, based on equations to solve.
    ///
    int set_readvars(class SimParams&  ///< pointer to simulation parameters
    );

    ///
    /// Given a Variable and a Silo File, read the data onto the grid.
    ///
    int read_variable2grid(
        class SimParams&,  ///< pointer to simulation parameters
        DBfile*,           ///< pointer to silo file.
        string,            ///< name of mesh to read from.
        string,            ///< variable name to read.
        long int           ///< number of cells expected.
    );
};

#endif  // if SILO

#endif  // if not DATAIO_SILO_H
