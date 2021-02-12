/// \file dataio_silo_MPI.h
/// \author Jonathan Mackey
///
/// Description:\n
/// Routines for writing data when PION is run in parallel with
/// MPI-based domain decomposition

#ifndef DATAIO_SILO_MPI_H
#define DATAIO_SILO_MPI_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SILO

#include "dataIO/dataio_base.h"
#include "dataIO/dataio_silo.h"
#include "decomposition/MCMD_control.h"
#include <silo.h>
#include <vector>

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
        class SimParams&,   ///< pointer to simulation parameters
        std::string,        ///< FLOAT or DOUBLE for files.
        class MCMDcontrol*  ///< address of MCMD controller class.
    );

    /// Destructor (doensn't have much to do).
    virtual ~dataio_silo_pllel();

    ///
    /// Loops through levels in the nested-grid structure, and calls
    /// SaveLevelData() on each level.
    ///
    int OutputData(
        const string,                   ///< File to write to
        vector<class GridBaseClass*>&,  ///< grid pointers.
        class SimParams&,               ///< pointer to simulation parameters
        const long int  ///< number to stamp file with (e.g. timestep)
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
    int ReadHeader(
        string,           ///< file to read from
        class SimParams&  ///< simulation parameters
    );

    ///
    /// This reads the data variables in turn, and puts the data into
    /// the grid points, assuming the grid has been set up with paramters
    /// from the header, which should be read first.  It also assumes the
    /// domain decomposition has been done already, and will bug out if it
    /// hasn't.  The routine uses the PMPIO interface.
    ///
    virtual int ReadData(
        string,                         ///< file to read from
        vector<class GridBaseClass*>&,  ///< grid pointers.
        class SimParams&                ///< simulation parameters
    );

  protected:
    class MCMDcontrol* mpiPM;
    ///
    /// Choose filename based on outfile, group_rank, and counter.
    ///
    virtual int choose_pllel_filename(
        const string,  ///< filebase passed in from main code.
        const int,     ///< group_rank (i.e. which file I write to)
        const int,     ///< file counter to use (e.g. timestep).
        string&        ///< string to return filename in.
    );

    ///
    /// get the name of a mesh, currently "unigridXXXX" where XXXX=rank
    ///
    void mesh_name(
        const int,  ///< rank
        string&     ///< mesh name returned in this string.
    );

    ///
    /// Given myrank and mygroup, create directory string where data is
    /// held.  This is determined as
    /// rank_[global_rank]_domain_[rank_within_group]
    ///
    void set_dir_in_file(
        std::string&,  ///< directory name.
        const int,     ///< myrank (global).
        const int,     ///< myrank in group.
        const int      ///< level of grid in hierarchy
    );

    ///
    /// This writes the header and data for the simulation parameters, for
    /// a given level in the nested-grid structure.
    ///
    /// If the solver pointer is not null, it also writes some derived
    /// variables such as Temperature, Div(B), etc.
    ///
    int SaveLevelData(
        const string,          ///< File-base to write to
        class GridBaseClass*,  ///< grid pointer.
        class SimParams&,      ///< simulation parameters
        const long int         ///< timestep
    );

    ///
    /// Call once to setup arrays with the properties of the grid, for
    /// writing to the Silo Quadmesh.  This is for the local domain of the
    /// current processor.
    ///
    int setup_grid_properties(
        class GridBaseClass*,  ///< pointer to data.
        class SimParams&       ///< pointer to simulation parameters
    );

    ///
    /// allocate memory for arrays used to write data to files.
    /// Arrays are different sizes for parallel grid -- local Ncell instead
    /// of global Ncell.
    ///
    void create_data_arrays(
        class SimParams&  ///< pointer to simulation parameters
    );

    int numfiles;  ///< number of files to split the data into.

    ///
    /// Write a mulitmesh adjacency object (obsolete and doesn't work!)
    ///
    int write_multimeshadj(
        class SimParams&,      ///< pointer to simulation parameters
        DBfile*,               ///< pointer to silo file.
        class GridBaseClass*,  ///< pointer to data.
        string,                ///< multimesh name
        string                 ///< multimesh adjacency name.
    );

    ///
    /// Write an MRG tree object (replaces multimesh adjacency)
    ///
    int write_MRGtree(
        class SimParams&,      ///< pointer to simulation parameters
        DBfile*,               ///< pointer to silo file.
        class GridBaseClass*,  ///< pointer to data.
        string,                ///< multimesh name
        string                 ///< MRG tree name.
    );
};

#endif  // if PARALLEL

#endif  // if SILO

#endif  // if not DATAIO_SILO_MPI_H
