///
/// \file dataio_silo_MPI.cc
///
/// \author Jonathan Mackey
///
/// This file contains the member functions of the parallel I/O class
/// for writing silo files.  It uses the PMPIO interface.  Note the
/// #def variable SILO must be set in the Makefile.
///
/// Modifications:
/// - 2010-02-03 JM: changed a few things to fix compiler warnings;
///   esp. tests for equality.
///  - 2010-02-04 JM: Added multimesh-adjacency object write so that I
///    can get streamlines to cross boundaries when plotting with Visit.
///  - 2010-02-05 JM: multimesh-adjacency and MRG tree connectivity
///    doesn't work with VisIt.  Very annoying.
///  - 2010-02-06 JM: Found a way to get multimesh-adjacency
///    connectivity working with VisIt.  Still no streamlines across
///    domains, but countours match (in 2D, 3D has a bug...)
///  - 2010-02-17 JM: Set numfiles to make files with max. size of
///    about 1GB.  For MHD there are about 120bytes per cell in the
///    file.
///  - 2010-04-11 JM: parallel class gets its own
///  setup_write_variables() class so that it can save disk space by
///  only writing primitive variables. (tidied up comments too).
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010-04-25 JM: renamed parallel choose_filename to choose_pllel_filename()
/// - 2010-07-20/21 JM: Work on new dataio structure: replaced dbfile with
/// dp_ptr
///    where appropriate.  Need a class pointer for the generic I/O interface.
/// - 2010.07.23 JM: removed obselete read_header(),
///    write_header() functions.
/// - 2010.10.01 JM: Spherical coordinates added.
/// - 2010.11.03 JM: Removed MM/Testing ifdefs for memory management
///    Also added Ndigits for width of counter in filename.
///    Also removed endl for c-style end of lines to avoid flushing o/p.
/// - 2010.11.19 JM: Got rid of testing myalloc() myfree() functions.
/// - 2011.03.02 JM: Added ability to write multiple column density data.
///                  Improved tracer variable handling (MAX_NVAR possible now).
/// - 2011.03.21 JM: Updated column-density variables for new cell interface
/// functions.
/// - 2011.03.22 JM: Removed setup_write_variables() function -- now use serial
/// version.
/// - 2011.10.24 JM: wrapped most of the info-reporting with ifdef-testing.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.01.28 JM: Changes for new code structure.
/// - 2015.06.13/15 JM: Changed datatype (FLOAT/DOUBLE) to a runtime
///    parameter, set in the constructor. (more 2015.06.18)

#ifdef SILO

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"


#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#include "tools/command_line_interface.h"
#endif  // NDEBUG

#ifndef PARALLEL
#error "PARALLEL not defined!  don't compile dataio_silo_MPI.cc without it!"
#endif
#include "dataio_silo_MPI.h"
#include <cstring>
#include <mpi.h>
#include <sstream>

// ##################################################################
// ##################################################################

dataio_silo_pllel::dataio_silo_pllel(
    class SimParams &SimPM,  ///< pointer to simulation parameters
    std::string dtype,       ///< FLOAT or DOUBLE for files.
    class Sub_domain *p) :
    dataio_silo(SimPM, dtype),
    mpiPM(p)
{
#ifndef NDEBUG
  spdlog::info("Setting up parallel Silo I/O class.");
#endif
  numfiles = -1;
  return;
}

// ##################################################################
// ##################################################################

dataio_silo_pllel::~dataio_silo_pllel()
{
#ifndef NDEBUG
  spdlog::info("Deleting parallel Silo I/O class.");
#endif
  return;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::ReadHeader(
    string infile,          ///< file to read from
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  int err  = 0;
  silofile = infile;
#ifndef NDEBUG
  spdlog::debug(
      "Rank: {}\tReading Header from file: {}", mpiPM->get_myrank(), silofile);
#endif

  if (!file_exists(silofile)) {
    spdlog::error(
        "{}: {}", "dataio_silo_pllel::ReadHeader() File not found, \
               myrank follows",
        mpiPM->get_myrank());
  }

  int group_rank = 0, myrank_group = 0;
  string file_id = "read_header";
  int num_files  = 1;
  err            = mpiPM->silo_pllel_init(
      num_files, "READ", file_id, &group_rank, &myrank_group);
  if (err) {
    spdlog::error("{}: {}", "mpiPM->silo_pllel_init() returned err", err);
  }

  //
  // Now wait for baton, and open into /header directory in file.
  //
  *db_ptr      = 0;
  string mydir = "/header";
  err = mpiPM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr)) {
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_wait_for_file() returned err", err);
  }

  //
  // Now read the header, and also NUM_FILES, which tells me how many files
  // there are for when I need to read in the data later on.
  //
  err               = read_simulation_parameters(SimPM);
  dataio_silo::ndim = SimPM.ndim;
  if (err) {
    spdlog::error(
        "{}: {}", "dataio_silo_MPI::ReadHeader() error reading header \
               from silo file",
        err);
  }

  err += DBReadVar(*db_ptr, "NUM_FILES", &numfiles);
  if (err) {
    numfiles = 1;
    err      = 0;
#ifndef NDEBUG
    spdlog::warn("Warning didn't read NUM_FILES from silo file.");
#endif
  }

  //
  // Finished Local work; hand off baton and wait!
  //
  err = mpiPM->silo_pllel_finish_with_file(file_id, db_ptr);
  if (err)
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_finish_with_file() returned err", err);

  SimPM.levels.resize(SimPM.grid_nlevels);
  mpiPM = &SimPM.levels[0].sub_domain;

#ifndef NDEBUG
  spdlog::debug(
      "Rank: {}\tFINISHED reading Header from file: {}", mpiPM->get_myrank(),
      silofile);
#endif
  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::ReadData(
    string infile,
    vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
    class SimParams &SimPM              ///< pointer to simulation parameters
)
{
  dataio_silo::gp = cg[0];
  if (!gp)
    spdlog::error(
        "{}: {}", "dataio_silo_pllel::ReadData() null pointer to grid!",
        fmt::ptr(gp));

  int err = 0;
#ifndef NDEBUG
  spdlog::debug(
      "\n----Rank: {}\tReading Data from files: {}", mpiPM->get_myrank(),
      infile);
#endif

  // set grid properties for quadmesh
  if (!have_setup_gridinfo) {
    err = setup_grid_properties(gp, SimPM);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}",
          "dataio_silo_pllel::ReadData() error setting up grid_props", 0, err);
  }

  // check numfiles was read from header
  if (numfiles < 0)
    spdlog::error(
        "{}: {}", "Failed to get numfiles from header! but error not detected!",
        numfiles);
  //  cout <<"\tnumfiles="<<numfiles<<" and nproc="<<mpiPM->nproc<<"\n";

  //
  // Get ready to open the file:
  //
  int group_rank = 0, myrank_ingroup = 0;
  string file_id = "read_data";
  err            = mpiPM->silo_pllel_init(
      numfiles, "READ", file_id, &group_rank, &myrank_ingroup);
  if (err)
    spdlog::error("{}: {}", "mpiPM->silo_pllel_init() returned err", err);

  // Choose correct filename based on numfiles.
  // Also choose correct directory name, assuming that nproc is the same as
  // for when the file was generated.
  //
  // replace 0000 in filename with my group rank, based on numfiles.
  //
  ostringstream temp;
  temp.fill('0');
  temp.str("");
  string::size_type pos = infile.find("0000");
  if (pos == string::npos) {
    spdlog::debug(
        "didn't fine 0000 in file, so assuming reading a single file");
    // spdlog::error("{}: {}", "not a parallel i/o filename!",infile);
  }
  else {
    temp.str("");
    temp.width(4);
    temp << group_rank;
    infile.replace(pos, 4, temp.str());
#ifndef NDEBUG
    spdlog::debug("\tNew infile: {}", infile);
#endif
    temp.str("");
  }
  silofile.clear();
  silofile = infile;

  //
  // choose my directory name (this function hardcoded for reading
  // only level 0 of grid!)
  //
  string mydir;
  set_dir_in_file(mydir, mpiPM->get_myrank(), myrank_ingroup, 0);
  if (SimPM.grid_nlevels > 0)
    spdlog::error("{}: {}", "dataio silo MPI levels", SimPM.grid_nlevels);
  // temp.str("");
  // temp << "level_"; temp.width(4); temp << mpiPM->get_myrank() <<
  // "_domain_"; temp.width(4); temp<<myrank_ingroup; mydir = temp.str();
  // temp.str("");

  *db_ptr = 0;
  err     = mpiPM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr))
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_wait_for_file() returned err", err);

  if (silo_filetype == DB_HDF5) {
    // char *compress = DBGetCompression();
    //    cout <<"compression="<<compress<<"\n";
    // DBSetCompression(compress);
    int friendly = DBGetFriendlyHDF5Names();
    DBSetFriendlyHDF5Names(friendly);
  }

  //
  // Check extents of quadmesh are what i have determined my domain to be...
  //  err = check_subdomain_extents(db_ptr);
  //
  temp.str("");
  temp << "unigrid";
  temp.width(4);
  temp << mpiPM->get_myrank();
  string meshname = temp.str();
  DBquadmesh *qm  = DBGetQuadmesh(*db_ptr, meshname.c_str());
  if (!qm)
    spdlog::error(
        "{}: {}", "failed to find quadmesh named as follows", meshname);
  int n;
  //
  // first check dimensions.
  //
  if ((n = qm->ndims) != ndim)
    spdlog::error("{}: {}", "bad dimensionality in qmesh", n);
  //
  // check number of cells.
  //
  int nn = 1;
  for (int i = 0; i < ndim; i++)
    nn *= mpiPM->get_directional_Ncells(i) + 1;
  if ((n = qm->nnodes) != nn) {
    spdlog::error("nnodes={} and ncell={}", n, mpiPM->get_Ncell());
    spdlog::error("{}: {}", "bad number of nodes", n - mpiPM->get_Ncell());
  }
  //
  // ---------- check grid nodes are where we expect ----------
  //
  // Check that datatype is what we are expecting!  If not, then
  // delete data arrays, reset datatype, and re-create data arrays.
  //
  if (qm->datatype != silo_dtype) {
    spdlog::debug(
        "HMMM: file has datatype {} but I am trying to read datatype {}\n    DB_INT=16, DB_SHORT=17, DB_LONG=18, DB_FLOAT=19, DB_DOUBLE=20, DB_CHAR=21, DB_LONG_LONG=22, DB_NOTYPE=25\nReadData() quadmesh has type={} but expecting type={}\n\t... resetting datatype for this file",
        qm->datatype, silo_dtype, qm->datatype, silo_dtype);
    delete_data_arrays();
    silo_dtype = qm->datatype;
    create_data_arrays(SimPM);
  }
  //
  // qm->coords is a <void **> pointer, so we have to reinterpret
  // it as either float or double.
  //
  if (silo_dtype == DB_FLOAT) {
    float **meshcoords = reinterpret_cast<float **>(qm->coords);
    //
    // check origin of subgrid is at the right place.
    //
    float *nx = reinterpret_cast<float *>(nodex);
    if (!pconst.equalD(nx[0], meshcoords[XX][0])) {
      spdlog::debug(
          "I think x[0]={}, silo file says x[0]={}", nx[0], meshcoords[XX][0]);
      spdlog::error(
          "{}: {}", "XX mesh not at right place...", nx[0] - meshcoords[XX][0]);
    }
    if (ndim > 1) {
      nx = reinterpret_cast<float *>(nodey);
      if (!pconst.equalD(nx[0], meshcoords[YY][0])) {
        spdlog::debug(
            "nodey = {} and coords[y]= {}\n", nx[0], meshcoords[YY][0]);
        spdlog::error(
            "{}: {}", "YY mesh not at right place...",
            nx[0] - meshcoords[YY][0]);
      }
    }
    if (ndim > 2) {
      nx = reinterpret_cast<float *>(nodez);
      if (!pconst.equalD(nx[0], meshcoords[ZZ][0])) {
        spdlog::debug(
            "nodez = {} and coords[z]= {}\n", nx[0], meshcoords[ZZ][0]);
        spdlog::error(
            "{}: {}", "ZZ mesh not at right place...",
            nx[0] - meshcoords[ZZ][0]);
      }
    }
  }
  else {
    //
    // do the same thing for double precision.
    //
    double **meshcoords = reinterpret_cast<double **>(qm->coords);
    //
    // check origin of subgrid is at the right place.
    //
    double *nx = reinterpret_cast<double *>(nodex);
    if (!pconst.equalD(nx[0], meshcoords[XX][0])) {
      spdlog::error(
          "{}: {}", "XX mesh not at right place...", nx[0] - meshcoords[XX][0]);
    }
    if (ndim > 1) {
      nx = reinterpret_cast<double *>(nodey);
      if (!pconst.equalD(nx[0], meshcoords[YY][0])) {
        spdlog::debug(
            "nodey = {} and coords[y]= {}\n", nx[0], meshcoords[YY][0]);
        spdlog::error(
            "{}: {}", "YY mesh not at right place...",
            nx[0] - meshcoords[YY][0]);
      }
    }
    if (ndim > 2) {
      nx = reinterpret_cast<double *>(nodez);
      if (!pconst.equalD(nx[0], meshcoords[ZZ][0])) {
        spdlog::debug(
            "nodez = {} and coords[z]= {}\n", nx[0], meshcoords[ZZ][0]);
        spdlog::error(
            "{}: {}", "ZZ mesh not at right place...",
            nx[0] - meshcoords[ZZ][0]);
      }
    }
  }
  DBFreeQuadmesh(qm);  // qm=0;

  //
  // now read each variable in turn from the mesh
  //
  err = set_readvars(SimPM);
  if (err) spdlog::error("{}: {}", "failed to set readvars in ReadData", err);
  for (std::vector<string>::iterator i = readvars.begin(); i != readvars.end();
       ++i) {
    err =
        read_variable2grid(SimPM, *db_ptr, meshname, (*i), mpiPM->get_Ncell());
    if (err)
      spdlog::error(
          "{}: {}", "dataio_silo_pllel::ReadData() error reading variable",
          (*i));
  }
  // Now assign Ph to be equal to P for each cell.
  cell *cpt = gp->FirstPt();
  do {
    for (int v = 0; v < SimPM.nvar; v++)
      cpt->Ph[v] = cpt->P[v];
  } while ((cpt = gp->NextPt(cpt)) != 0);

  //
  // Finished Local work; hand off baton and wait!
  //
  err = mpiPM->silo_pllel_finish_with_file(file_id, db_ptr);
  if (err)
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_finish_with_file() returned err", err);
  mpiPM->barrier("dataio_silo_pllel__ReadData");

#ifndef NDEBUG
  spdlog::debug(
      "----Rank: {}\tFINISHED reading Data from file: {}", mpiPM->get_myrank(),
      silofile);
#endif
  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::OutputData(
    const string outfilebase,
    vector<class GridBaseClass *> &cg,  ///< grid pointers.
    class SimParams &SimPM,             ///< simulation parameters
    const long int file_counter         ///< timestep
)
{
  int err = 0;
  // loop over grid refinement levels, save one file per level
  for (int l = 0; l < SimPM.grid_nlevels; l++) {

    if (!cg[l])
      spdlog::error(
          "{}: {}", "dataio_silo::OutputData() null pointer!", fmt::ptr(cg[l]));
    dataio_silo::gp = cg[l];
    mpiPM           = &(SimPM.levels[l].sub_domain);

    // write a different file for each level in the NG grid.
    string fbase;
    if (SimPM.grid_nlevels > 1) {
      ostringstream temp;
      temp << outfilebase << "_level";
      temp.width(2);
      temp.fill('0');
      temp << l;
      fbase = temp.str();
    }
    else {
      fbase = outfilebase;
    }

    err = SaveLevelData(fbase, cg[l], SimPM, file_counter);
    if (0 != err)
      spdlog::error("{}: Expected {} but got {}", "saveleveldata", 0, err);
  }
  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::SaveLevelData(
    const string outfilebase,    ///< filename
    class GridBaseClass *cg,     ///< grid pointers.
    class SimParams &SimPM,      ///< simulation parameters
    const long int file_counter  ///< timestep
)
{
  int err         = 0;
  dataio_silo::gp = cg;
  if (!gp)
    spdlog::error(
        "{}: {}", "dataio_silo_pllel::SaveLevelData() null grid!",
        fmt::ptr(gp));

  int level = 0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (gp == SimPM.levels[l].grid) level = l;
    // cout <<"saving level "<<level<<"\n";
  }

#ifndef NDEBUG
  spdlog::info(
      "----dataio_silo_pllel::SaveLevelData() Writing data to filebase: {}",
      outfilebase);
#endif

  // First initialise the I/O.
  if (silo_filetype == DB_HDF5) {
    DBSetCompression("METHOD=GZIP LEVEL=1");
    DBSetFriendlyHDF5Names(1);
#ifndef NDEBUG
    spdlog::info(" *** setting compression.");
#endif
  }

  ///
  /// numfiles is hardcoded for now, but could set it as a runtime
  /// parameter at some stage.  It produces files of up to about
  /// 1.2GB at the moment, and then split into 2 files.
  ///
  int threshold = static_cast<int>(SimPM.Ncell / 3.0e7) + 1;
  if (threshold < 1) threshold = 1;
  //
  // Have to make sure that we write at most one file per process!
  //
  if (mpiPM->get_nproc() > threshold)
    numfiles = threshold;
  else
    numfiles = mpiPM->get_nproc();

#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() running pmpio_init");
#endif

  int group_rank = 0, myrank_ingroup = 0;
  string file_id = "write_data";
  err            = mpiPM->silo_pllel_init(
      numfiles, "WRITE", file_id, &group_rank, &myrank_ingroup);
  if (err)
    spdlog::error("{}: {}", "mpiPM->silo_pllel_init() returned err", err);

#ifndef NDEBUG
  spdlog::debug(
      "myrank: {}\tnumfiles: {}\tmy_group: {}\tmy_index_in_group: {}",
      mpiPM->get_myrank(), numfiles, group_rank, myrank_ingroup);
#endif  // NDEBUG

  //
  // Choose output filename:
  //
#ifndef NDEBUG
  spdlog::debug("setting strings... outfilebase={}", outfilebase);
#endif  // NDEBUG

  choose_pllel_filename(outfilebase, group_rank, file_counter, silofile);
#ifndef NDEBUG
  spdlog::info("string for outfile set...");
#endif  // NDEBUG

  //
  // Choose directory within silo file.
  //
  string mydir;
  set_dir_in_file(mydir, mpiPM->get_myrank(), myrank_ingroup, level);
  ostringstream temp;
  temp.fill('0');
  temp.str("");
  // temp << "level_";
  // temp.width(4);
  // temp << lev << "_domain_";
  // temp.width(4);
  // temp << myrank_ingroup;
  // mydir = temp.str();
  // temp.str("");
  if (mydir.size() + 1 >= strlength)
    spdlog::error("{}: {}", "string too large", mydir);

  err = 0;
  // set grid properties for quadmesh
#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() grid properties");
#endif
  err = setup_grid_properties(gp, SimPM);
  if (err)
    spdlog::error(
        "{}: {}", "dataio_silo_pllel::SaveLevelData() grid_props", err);
#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() grid props done");
#endif

  if (!have_setup_writevars) {
    // set what data to write to the mesh.
#ifndef NDEBUG
    spdlog::info("----dataio_silo_pllel::SaveLevelData() write variables");
#endif
    err = setup_write_variables(SimPM);
    if (err)
      spdlog::error(
          "{}: {}", "dataio_silo_pllel::SaveLevelData() write-vars", err);
#ifndef NDEBUG
    spdlog::info("----dataio_silo_pllel::SaveLevelData() write vars done");
#endif
  }

  //
  // performance options:
  //
  /* must be communicated before serial region */
  int ext_size = 2 * ndim;
  int nmesh    = mpiPM->get_nproc();
  double extents[ext_size * nmesh];
  int zonecounts[nmesh];

  SimPM.levels[level].sub_domain.gather_ncells(zonecounts, 0);
  SimPM.levels[level].sub_domain.gather_extents(extents, 0);

  //
  // Wait for my turn to write to the file.
  //
#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() waiting for baton");
#endif
  /*
   * TODO: non-root processes call this to pass required data to rank 0 when it
   * calls it later on. This is horrific and must be changed. It also relies on
   * rank 0 writing to file first
   */
  if (mpiPM->get_myrank() > 0) {
    write_multimeshadj(SimPM, nullptr, gp, "", "");
  }

  *db_ptr = 0;
  err     = mpiPM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr))
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_wait_for_file() returned err", err);

    //
    // Have got the baton, now, so the file is mine to write to.
    // local work here... each proc write their part of the grid.
    //
    // Generate Mesh in file.
    //
#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() generating quadmesh");
#endif
  string meshname;
  mesh_name(mpiPM->get_myrank(), meshname);
  err = generate_quadmesh(*db_ptr, meshname, SimPM);
  if (err)
    spdlog::error("{}: {}", "dataio_silo_pllel::SaveLevelData() quadmesh", err);
#ifndef NDEBUG
  spdlog::info("----dataio_silo_pllel::SaveLevelData() quadmesh generated");
#endif

  //
  // now write each variable in turn to the mesh
  //
#ifndef NDEBUG
  spdlog::debug("----dataio_silo_pllel::SaveLevelData() create data arrays");
#endif

  create_data_arrays(SimPM);

#ifndef NDEBUG
  spdlog::debug("----dataio_silo_pllel::SaveLevelData() arrays created");
#endif

  for (std::vector<string>::iterator i = varnames.begin(); i != varnames.end();
       ++i) {
#ifndef NDEBUG
    spdlog::debug("\twriting variable {} to file {}", (*i), silofile);
#endif

    err = write_variable2mesh(SimPM, *db_ptr, meshname, (*i));

#ifndef NDEBUG
    spdlog::debug("\t\tvariable {} written.", (*i));
#endif

    if (err)
      spdlog::error(
          "{}: {}", "dataio_silo_pllel::SaveLevelData() writing variable",
          (*i));
  }

#ifndef NDEBUG
  spdlog::debug("----dataio_silo_pllel::SaveLevelData() deleting arrays");
#endif

  delete_data_arrays();

#ifndef NDEBUG
  spdlog::debug("----dataio_silo_pllel::SaveLevelData() arrays deleted");
#endif

  // ONLY DO THIS IF I AM ROOT PROCESSOR *IN EACH GROUP*
  // GET CURRENT DIR, MOVE OUT, MAKE HEADER,
  // MOVE BACK TO CURRENT DIR
  if (myrank_ingroup == 0) {
    //    cout <<"myrank_ingroup="<<myrank_ingroup<<" so writing header to
    //    file.\n";
    // Am keeping file open, so no need to re-open it anymore.
    //    *db_ptr =0;
    //    *db_ptr = DBOpen(fname,DB_UNKNOWN,DB_APPEND);

    char datadir[strlength];
    err += DBGetDir(*db_ptr, datadir);
    DBSetDir(*db_ptr, "/");

    // WRITE HEADER
    DBMkDir(*db_ptr, "header");
    DBSetDir(*db_ptr, "/header");
    //
    // Write Number of files to header directory of file.
    //
    int dim1[1];
    dim1[0]   = 1;
    int nproc = mpiPM->get_nproc();
    err += DBWrite(*db_ptr, "NUM_FILES", &numfiles, dim1, 1, DB_INT);
    err += DBWrite(*db_ptr, "MPI_nproc", &nproc, dim1, 1, DB_INT);
    err += DBWrite(*db_ptr, "grid_level", &level, dim1, 1, DB_INT);
    dim1[0] = 3;
    err += DBWrite(
        *db_ptr, "level_xmin", &(SimPM.levels[level].Xmin), dim1, 1, DB_DOUBLE);
    err += DBWrite(
        *db_ptr, "level_xmax", &(SimPM.levels[level].Xmax), dim1, 1, DB_DOUBLE);
    //    err = write_header(*db_ptr);
    err = write_simulation_parameters(SimPM);
    if (err)
      spdlog::error("{}: {}", "dataio_silo_pllel::SaveLevelData() header", err);

    DBSetDir(*db_ptr, "/");
    DBSetDir(*db_ptr, datadir);
    //    cout <<"Finished writing header to file.\n";
  }  // if root processor of *group*

  if (mpiPM->get_myrank() == 0) {
    //
    // WRITE MULTIMESH DATA NOW TO FIRST PARTIAL FILE IF ROOT PROC.
    //
    string mm_name   = "MultiMesh";
    string mma_name  = "Domain_Decomposition";
    string mrgt_name = "";  // MRGtree";

    //    cout <<"Writing multimesh data to file...\n";
    // err = write_multimesh(*db_ptr, gp, mm_name, mma_name);
    // if (err) spdlog::error("{}: {}", "Failed to write multimesh Object",err);
    char datadir[strlength];
    int err = DBGetDir(*db_ptr, datadir);
    DBSetDir(*db_ptr, "/");

    //
    // Strip the path from outfilebase to get just the filename; do this by
    // reverse finding '/' from outfilebase and taking the substring
    // following its position.  As far as I can tell this is foolproof...
    //
    string fname;
    string::size_type p = outfilebase.rfind("/");
    if (p != string::npos) {
      fname = outfilebase.substr(p + 1);
    }
    else
      fname = outfilebase;
    //    cout <<"outfilebase="<<outfilebase<<" and fname="<<fname<<"\n";

    //  string mm_name="mesh";
    char **mm_names = 0;
    int *meshtypes  = 0;
    int *groups     = 0,  // lists which group each process is in.
        *ranks      = 0;  // lists which rank each process has in its group.

    meshtypes = mem.myalloc(meshtypes, nmesh);
    mm_names  = mem.myalloc(mm_names, nmesh);
    for (int i = 0; i < nmesh; i++)
      mm_names[i] = mem.myalloc(mm_names[i], strlength);
    groups = mem.myalloc(groups, nmesh);
    ranks  = mem.myalloc(ranks, nmesh);

    //
    // Assign groups and rank_in_groups to each process:
    //
    for (int i = 0; i < nmesh; i++) {
      mpiPM->silo_pllel_get_ranks(file_id, i, &(groups[i]), &(ranks[i]));
    }

    int csys = -1;
    // csys = DB_QUADMESH;
    if (SimPM.coord_sys == COORD_CRT)
      csys = DB_QUAD_RECT;
    else if (SimPM.coord_sys == COORD_CYL)
      csys = DB_QUAD_RECT;
    else if (SimPM.coord_sys == COORD_SPH)
      csys = DB_QUAD_RECT;
    else
      spdlog::error("{}: {}", "bad coord system", SimPM.coord_sys);
    for (int i = 0; i < nmesh; i++) {
      meshtypes[i] = csys;
      temp.str("");
      string tfn;
      if (groups[i] != group_rank) {
        //
        // the root processor's group file doesn't need to be in the
        // path, so only other groups' files have to get listed.
        //
        // temp << fname <<"_"; temp.width(4); temp <<groups[i]<<".";
        // temp.width(6); temp <<SimPM.timestep<<".silo:";
        choose_pllel_filename(fname, groups[i], file_counter, tfn);
        temp << tfn.c_str() << ":";
      }
      tfn.clear();
      set_dir_in_file(tfn, i, ranks[i], level);
      // temp << "/level_"; temp.width(4); temp << lev << "_domain_";
      // temp.width(4); temp<<ranks[i];
      temp << tfn.c_str();
      tfn.clear();
      mesh_name(i, tfn);
      temp << "/" << tfn.c_str();
      // temp <<"/unigrid"; temp.width(4); temp<<i;
      if (temp.str().length() > strlength - 1)
        spdlog::error(
            "{}: {}", "directory name too long!!!", temp.str().length());
      strcpy(mm_names[i], temp.str().c_str());
    }

    DBoptlist *mm_opts = DBMakeOptlist(7);
    DBAddOption(mm_opts, DBOPT_DTIME, &SimPM.simtime);
    DBAddOption(mm_opts, DBOPT_CYCLE, &SimPM.timestep);
    //  DBAddOption(mm_opts,DBOPT_ADJACENCY_NAME,mma_name.c_str()); //
    //  doesn't exist!!!
    // char mrgt[256]; strcpy(mrgt,mrgt_name.c_str());
    // DBAddOption(mm_opts,DBOPT_MRGTREE_NAME,mrgt);
    int blockorigin = 0;
    DBAddOption(mm_opts, DBOPT_BLOCKORIGIN, &blockorigin);

    int externalzones[nmesh];
    for (int v = 0; v < nmesh; v++) {
      externalzones[v] = 0;  // set to 1 if domain has zones outside multimesh.
    }

    DBAddOption(mm_opts, DBOPT_EXTENTS_SIZE, &ext_size);
    DBAddOption(mm_opts, DBOPT_EXTENTS, extents);
    DBAddOption(mm_opts, DBOPT_ZONECOUNTS, zonecounts);
    DBAddOption(mm_opts, DBOPT_HAS_EXTERNAL_ZONES, externalzones);

    // write the multimesh
    err = DBPutMultimesh(
        *db_ptr, mm_name.c_str(), nmesh, mm_names, meshtypes, mm_opts);
    if (err)
      spdlog::error(
          "{}: {}", "dataio_silo_pllel::SaveLevelData()multimesh", err);
#ifndef NDEBUG
    spdlog::debug("----dataio_silo_pllel::SaveLevelData() quadmesh written");
#endif
    DBClearOptlist(mm_opts);

    // now write the multivars
    // err = DBPutMultivar(*db_ptr, mv_name, nmvar, mv_names, vartypes,
    // mv_opts);
    // I'm going to reuse all the multimesh vars to the multivar object...
    for (std::vector<string>::iterator i = varnames.begin();
         i != varnames.end(); ++i) {
      // multivar name
      string vname;
      vname.clear();
      vname = (*i);
      //
      // multivar names:
      //
      for (int ii = 0; ii < nmesh; ii++) {
        meshtypes[ii] = DB_QUADVAR;
        temp.str("");
        string tfn;
        if (groups[ii] != group_rank) {
          //
          // the root processor's group file doesn't need to be in the
          // path, so only other groups' files have to get listed.
          //
          // temp << fname <<"_"; temp.width(4); temp
          // <<groups[ii]<<"."; temp.width(6); temp
          // <<SimPM.timestep<<".silo:";
          choose_pllel_filename(fname, groups[ii], file_counter, tfn);
          temp << tfn.c_str() << ":";
        }
        tfn.clear();
        set_dir_in_file(tfn, ii, ranks[ii], level);
        // temp << "/level_"; temp.width(4); temp << lev << "_domain_";
        // temp.width(4); temp<<ranks[ii];
        temp << tfn.c_str();
        temp << "/" << vname;
        if (temp.str().length() > strlength)
          spdlog::error(
              "{}: {}", "multivar name too long!!!", temp.str().length());
        strcpy(mm_names[ii], temp.str().c_str());
      }

      //
      // Now the performance options for the multivars
      //
      DBAddOption(mm_opts, DBOPT_DTIME, &SimPM.simtime);
      DBAddOption(mm_opts, DBOPT_CYCLE, &SimPM.timestep);
      DBAddOption(mm_opts, DBOPT_BLOCKORIGIN, &blockorigin);
      // ext_size=2;
      // extents[] is the max/min values of the variable -- no time to
      // calculate that.
      char mn[256];
      strcpy(mn, mm_name.c_str());
      DBAddOption(mm_opts, DBOPT_MMESH_NAME, mn);

      err = DBPutMultivar(
          *db_ptr, vname.c_str(), nmesh, mm_names, meshtypes, mm_opts);
      if (err)
        spdlog::error(
            "{}: {}", "dataio_silo_pllel::SaveLevelData() variable", (*i));
      DBClearOptlist(mm_opts);
    }

    //
    // Free memory
    //
    DBFreeOptlist(mm_opts);
    for (int i = 0; i < nmesh; i++)
      mm_names[i] = mem.myfree(mm_names[i]);
    mm_names  = mem.myfree(mm_names);
    meshtypes = mem.myfree(meshtypes);
    groups    = mem.myfree(groups);
    ranks     = mem.myfree(ranks);

    DBSetDir(*db_ptr, "/");
    DBSetDir(*db_ptr, datadir);

    //    cout <<"Finished writing multimesh data to file...\n";

    //
    // Now write a mulitmeshadj object (but only if nproc>1).
    //
    if (mpiPM->get_nproc() > 1) {
      err = write_multimeshadj(SimPM, *db_ptr, gp, mm_name, mma_name);
      if (err)
        spdlog::error(
            "{}: {}", "Failed to write multimesh Adjacency Object", err);
    }
  }  // if root processor of whole simulation

  //
  // Finished Local work; hand off baton and wait!
  //
  err = mpiPM->silo_pllel_finish_with_file(file_id, db_ptr);
  if (err)
    spdlog::error(
        "{}: {}", "mpiPM->silo_pllel_finish_with_file() returned err", err);
  mpiPM->barrier("dataio_silo_pllel__SaveLevelData");

  //   cout <<"Got past barrier... finished outputting silo data.\n\n";

#ifndef NDEBUG
  spdlog::debug(
      "----dataio_silo_pllel::SaveLevelData() Finished writing data to file: {}",
      silofile);
#endif
  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::choose_pllel_filename(
    const string fbase,      ///< filebase passed in from main code.
    const int igroup,        ///< group_rank (i.e. which file I write to)
    const int file_counter,  ///< file counter to use (e.g. timestep).
    string &outfile          ///< write filename to this string.
)
{
  ostringstream temp;
  temp.str("");
  temp.fill('0');
  temp << fbase << "_";
  temp.width(4);
  temp << igroup << ".";
  if (file_counter >= 0) {
    temp.width(Ndigits);
    temp << file_counter << ".";
  }
  temp << "silo";
  outfile.clear();
  outfile = temp.str();
  if (outfile.size() + 1 >= strlength)
    spdlog::error("{}: {}", "string too large", outfile);
  temp.str("");
  return 0;
}

// ##################################################################
// ##################################################################

int dataio_silo_pllel::setup_grid_properties(
    class GridBaseClass *grid,  ///< pointer to data.
    class SimParams &SimPM      ///< pointer to simulation parameters
)
{
  // set grid parameters
  // This version is for the local domain of the current processor.
  if (!grid)
    spdlog::error(
        "{}: {}", "dataio_silo::setup_grid_properties() null ptr",
        fmt::ptr(grid));
  double dx               = grid->DX();
  dataio_silo::ndim       = SimPM.ndim;
  dataio_silo::vec_length = SimPM.eqnNDim;

  if (!nodedims) nodedims = mem.myalloc(nodedims, ndim);
  if (!zonedims) zonedims = mem.myalloc(zonedims, ndim);

    //
    // node_coords is a void pointer, so if we are writing data in
    // single or double precision then we need different allocation
    // calls.  Same for nodex, nodey, nodez.
    //
    // setup arrays with locations of nodes in coordinate directions.
    // This differs from the serial code in that we use directional_Ncells, not
    // the global number of points NG.
    //
#ifdef WRITE_GHOST_ZONES
  int nx = mpiPM->get_directional_Ncells(XX) + 2 * SimPM.Nbc
           + 1;  // N cells, have N+1 nodes.
  int ny = mpiPM->get_directional_Ncells(YY) + 2 * SimPM.Nbc
           + 1;  // N cells, have N+1 nodes.
  int nz = mpiPM->get_directional_Ncells(ZZ) + 2 * SimPM.Nbc
           + 1;  // N cells, have N+1 nodes.
#else
  int nx =
      mpiPM->get_directional_Ncells(XX) + 1;  // for N cells, have N+1 nodes.
  int ny =
      mpiPM->get_directional_Ncells(YY) + 1;  // for N cells, have N+1 nodes.
  int nz =
      mpiPM->get_directional_Ncells(ZZ) + 1;  // for N cells, have N+1 nodes.
#endif

  if (silo_dtype == DB_FLOAT) {
    //
    // Allocate memory for node_coords, and set pointers.
    //
    float **d   = 0;
    float *posx = 0, *posy = 0, *posz = 0;

    if (node_coords) {
      posx = reinterpret_cast<float *>(nodex);
      posy = reinterpret_cast<float *>(nodey);
      posz = reinterpret_cast<float *>(nodez);
    }
    else {
      d           = mem.myalloc(d, ndim);
      node_coords = reinterpret_cast<void **>(d);
      posx        = mem.myalloc(posx, nx);
      if (ndim > 1) posy = mem.myalloc(posy, ny);
      if (ndim > 2) posz = mem.myalloc(posy, nz);
    }

    //
    // Assign data for nodex, nodey, nodez for this grid
    //
    for (int i = 0; i < nx; i++)
      posx[i] = static_cast<float>(mpiPM->get_Xmin(XX) + i * dx);
    nodex           = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;

    if (ndim > 1) {
      for (int i = 0; i < ny; i++)
        posy[i] = static_cast<float>(mpiPM->get_Xmin(YY) + i * dx);
      nodey           = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim > 2) {
      posz = mem.myalloc(posz, nz);
      for (int i = 0; i < nz; i++)
        posz[i] = static_cast<float>(mpiPM->get_Xmin(ZZ) + i * dx);
      nodez           = reinterpret_cast<void *>(posz);
      node_coords[ZZ] = nodez;
    }
  }
  else {
    //
    // Allocate double-precision memory for node_coords, and set
    // pointers.
    //
    double **d   = 0;
    double *posx = 0, *posy = 0, *posz = 0;

    if (node_coords) {
      posx = reinterpret_cast<double *>(nodex);
      posy = reinterpret_cast<double *>(nodey);
      posz = reinterpret_cast<double *>(nodez);
    }
    else {
      d           = mem.myalloc(d, ndim);
      node_coords = reinterpret_cast<void **>(d);
      posx        = mem.myalloc(posx, nx);
      if (ndim > 1) posy = mem.myalloc(posy, ny);
      if (ndim > 2) posz = mem.myalloc(posz, nz);
    }

    //
    // Assign data for nodex, nodey, nodez for this grid
    //
    for (int i = 0; i < nx; i++) {
      posx[i] = static_cast<double>(mpiPM->get_Xmin(XX) + i * dx);
    }
    nodex           = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;

    if (ndim > 1) {
      for (int i = 0; i < ny; i++) {
        posy[i] = static_cast<double>(mpiPM->get_Xmin(YY) + i * dx);
      }
      nodey           = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim > 2) {
      posz = mem.myalloc(posz, nz);
      for (int i = 0; i < nz; i++) {
        posz[i] = static_cast<double>(mpiPM->get_Xmin(ZZ) + i * dx);
      }
      nodez           = reinterpret_cast<void *>(posz);
      node_coords[ZZ] = nodez;
    }
  }

  nodedims[0] = nx;
  zonedims[0] = nx - 1;

  if (ndim > 1) {
    nodedims[1] = ny;
    zonedims[1] = ny - 1;
  }
  if (ndim > 2) {
    nodedims[2] = nz;
    zonedims[2] = nz - 1;
  }

  int nopts             = 4;
  dataio_silo::GridOpts = DBMakeOptlist(nopts);
  if (SimPM.coord_sys == COORD_CRT)
    silo_coordsys = DB_CARTESIAN;
  else if (SimPM.coord_sys == COORD_CYL)
    silo_coordsys = DB_CYLINDRICAL;
  else if (SimPM.coord_sys == COORD_SPH)
    silo_coordsys = DB_SPHERICAL;
  else
    spdlog::error("{}: {}", "bad coord system", SimPM.coord_sys);
  DBAddOption(GridOpts, DBOPT_COORDSYS, &silo_coordsys);
  DBAddOption(GridOpts, DBOPT_DTIME, &SimPM.simtime);
  DBAddOption(GridOpts, DBOPT_CYCLE, &SimPM.timestep);
  DBAddOption(GridOpts, DBOPT_NSPACE, &SimPM.ndim);
  // char temp[strlength];
  // strcpy(temp,uc.length.c_str());
  // DBAddOption(GridOpts,DBOPT_XUNITS,temp);

  have_setup_gridinfo = true;
  return 0;
}

// ##################################################################
// ##################################################################

void dataio_silo_pllel::create_data_arrays(
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{
  // first check if we have the data arrays set up yet.
  // We need at least one array for a scalar, and two more for a
  // vector.
  // This differs from serial code in that we use Ncell instead
  // of the global Ncell to allocate memory.
  //
  // cout <<"local Ncell="<<mpiPM->get_Ncell();
  // cout <<" and global Ncell is "<<SimPM.Ncell<<"\n";
  //
  // data0 is a void pointer, so we need to do different things for
  // float and double data (same for data1,data2).
  //
  if (!data0) {
    if (silo_dtype == DB_FLOAT) {
      float *d = 0;
      d        = mem.myalloc(d, mpiPM->get_Ncell());
      data0    = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d         = mem.myalloc(d, mpiPM->get_Ncell());
      data0     = reinterpret_cast<void *>(d);
    }
  }

  // set up array for mask variable for nested grid.
  if (!mask) {
    int *m = 0;
    m      = mem.myalloc(m, mpiPM->get_Ncell());
    mask   = reinterpret_cast<void *>(m);
  }

  //
  // If we are only writing scalar data, we don't need data1,data2,vec_data
  // so by setting vec_length=0 they don't get initialised.
  // This saves some memory.
  //
#if defined(SILO_SCALARS)
  vec_length = 0;
#endif

  //
  // If we need data1 and data2, and vec_data, create them too.
  //
  if ((vec_length > 1) && (!data1)) {
    if (silo_dtype == DB_FLOAT) {
      float *d = 0;
      d        = mem.myalloc(d, mpiPM->get_Ncell());
      data1    = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d         = mem.myalloc(d, mpiPM->get_Ncell());
      data1     = reinterpret_cast<void *>(d);
    }
  }
  if ((vec_length > 2) && (!data2)) {
    if (silo_dtype == DB_FLOAT) {
      float *d = 0;
      d        = mem.myalloc(d, mpiPM->get_Ncell());
      data2    = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d         = mem.myalloc(d, mpiPM->get_Ncell());
      data2     = reinterpret_cast<void *>(d);
    }
  }

  if ((vec_length > 1) && (!vec_data)) {
    if (silo_dtype == DB_FLOAT) {
      float **d = 0;
      d         = mem.myalloc(d, vec_length);
      vec_data  = reinterpret_cast<void **>(d);
    }
    else {
      double **d = 0;
      d          = mem.myalloc(d, vec_length);
      vec_data   = reinterpret_cast<void **>(d);
    }
    vec_data[0] = data0;
    if (vec_length > 1) vec_data[1] = data1;
    if (vec_length > 2) vec_data[2] = data2;
  }

  return;
}

// ##################################################################
// ##################################################################

//
// Write a mulitmesh adjacency object
//
int dataio_silo_pllel::write_multimeshadj(
    class SimParams &SimPM,    ///< pointer to simulation parameters
    DBfile *dbfile,            ///< pointer to silo file.
    class GridBaseClass *ggg,  ///< pointer to data.
    string mm_name,            ///< multimesh  name
    string mma_name            ///< multimeshadj name.
)
{
  int err         = 0;
  const int nmesh = mpiPM->get_nproc();

  int level = 0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (gp == SimPM.levels[l].grid) level = l;
    // cout <<"saving level "<<level<<"\n";
  }

  int meshtypes[nmesh], Nngb[nmesh], Sk[nmesh];
  for (int v = 0; v < nmesh; v++) {
    // meshtypes[v] = DB_QUAD_RECT;
    meshtypes[v] = DB_QUADMESH;
    Nngb[v] = Sk[v] = 0;
  }

  //
  // get number of neighbours in each mesh, and also the running sum of the
  // total number of neighbours, because we need some arrays of this length.
  //
  // Note neighbouring domains are not just in the six directions, but
  // also corner-abutting and edge-abutting domains.
  //
  // Stot is the total number of neighbours over all meshes.  Sk[v] is
  // the number of meshes in all neighbours up to and including v-1.
  //

  std::vector<int> ngb_list;
  mpiPM->gather_abutting_domains(ngb_list, Sk, Nngb, 0);

  int Stot = ngb_list.size();

  std::vector<int> offsets_list, directional_Ncells_list, ix_list;
  mpiPM->gather_offsets(offsets_list, 0);
  mpiPM->gather_directional_Ncells(directional_Ncells_list, 0);

  /* non-root ranks were only here to share data with root rank */
  if (mpiPM->get_myrank() > 0) {
    return 0;
  }

#ifndef NDEBUG
  spdlog::info("Writing multimesh adjacency object into Silo file.");
#endif

  //
  // Y ngb      = list of neighbour ids.
  // Y back     = reverse list of neighbour ids.
  // N nzones   = number of zones adjacent to each interface. (UNUSED)
  // N zonelist = index ranges of each adjacent zone. (UNUSED)
  // Y nnodes   = number of elements in nodelist == 15 for structured grid!
  // * nodelist = index ranges of each list of nodes..
  // ct       = counter for where we are in neighbour lists.
  //
  int ngb[Stot], back[Stot], nnodes[Stot], *nodelist[Stot];
  //
  // Each nodelist element has 15 elements regardless of ndim, and we
  // set some to Q=-1 if ndim<3
  //
  int Q      = -1;  // for unused values.
  int nl_len = 15;  // size of elements.
  for (int s = 0; s < Stot; s++) {
    nnodes[s]   = nl_len;
    ngb[s]      = ngb_list[s];
    nodelist[s] = mem.myalloc(nodelist[s], nl_len);
    for (int v = 0; v < nl_len; v++) {
      nodelist[s][v] = Q;
    }
  }

  //
  // loop over meshes and populate the neighbour lists.
  //
  int *offsets, *directional_Ncells;
  for (int v = 0; v < nmesh; v++) {
    offsets            = &offsets_list[v * mpiPM->get_ndim()];
    directional_Ncells = &directional_Ncells_list[v * mpiPM->get_ndim()];
    long int off1      = Sk[v];

    //
    // Assign reverse neighbour's id for each of myrank's neighbours.
    // If myrank=v, ngb[off1+i]=dom2, then back[off1+i] is v's index
    // in dom2's list of neigbours.
    //
    int dom2;
    long int off2;
    for (int i = 0; i < Nngb[v]; i++) {
      dom2 = ngb[off1 + i];
      off2 = Sk[dom2];
      if (off2 > Stot)
        spdlog::error(
            "{}: {}", "Counting error in loop to find back array", off2 - Stot);

      for (int s = 0; s < Nngb[dom2]; s++) {
        if (ngb[off2 + s] == v) back[off1 + i] = s;
      }  // loop over dom2's neighbours
    }    // loop over dom1's neighbours

    //
    // Now the nodelists[]: [0,1] are xmin,xmax [2,3] ymin,ymax [4,5] are
    // zmin,zmax
    // https://email.ornl.gov/pipermail/visit-users/attachments/20081017/a72e7c92/multimesh_cyrus.cc
    //
    // This page says that the ordering should be:
    // local[x0,x1,y0,y1,z0,z1]boundary[x0,x1,y0,y1,z0,z1]orientation[1,2,3]
    // It also requires domains which only share a corner, and it
    // requires nodes to have a global indexing system.
    //
    // We loop over mesh v's neighbours and assign nodelist[Sk[v]+i]
    //
    for (int i = 0; i < Nngb[v]; i++) {
      off1 = Sk[v] + i;
      if (off1 > Stot)
        spdlog::error(
            "{}: {}", "Counting error in loop over neighbours for mesh v",
            off1 - Stot);
      //
      // Local nodelist is the same for all of v's nodelists.
      //
      nodelist[off1][0] = offsets[XX] + 0;
      nodelist[off1][1] = offsets[XX] + directional_Ncells[XX];
      if (ndim > 1) {
        nodelist[off1][2] = offsets[YY] + 0;
        nodelist[off1][3] = offsets[YY] + directional_Ncells[YY];
      }
      if (ndim > 2) {
        nodelist[off1][4] = offsets[ZZ] + 0;
        nodelist[off1][5] = offsets[ZZ] + directional_Ncells[ZZ];
      }

      //
      // For neighbours nodelists, get their relative position in the
      // block structure.
      //
      int my_ix[MAX_DIM], ngb_ix[MAX_DIM], nx[MAX_DIM];
      // std::array<int,MAX_DIM> my_ix, ngb_ix, nx;
      for (int ii = 0; ii < MAX_DIM; ii++)
        my_ix[ii] = ngb_ix[ii] = nx[ii] = -1;
      mpiPM->get_domain_coordinates(v, my_ix);
      mpiPM->get_domain_coordinates(ngb[off1], ngb_ix);
      mpiPM->get_num_subdomains(nx);

      //
      // X-dir first.
      //
      if ((my_ix[XX] - ngb_ix[XX]) == 1) {
        nodelist[off1][6] = offsets[XX];
        nodelist[off1][7] = offsets[XX];
      }
      else if ((my_ix[XX] - ngb_ix[XX]) == -1) {
        nodelist[off1][6] = offsets[XX] + directional_Ncells[XX];
        nodelist[off1][7] = offsets[XX] + directional_Ncells[XX];
      }
      else if (my_ix[XX] == ngb_ix[XX]) {
        nodelist[off1][6] = offsets[XX];
        nodelist[off1][7] = offsets[XX] + directional_Ncells[XX];
      }
      else {
        spdlog::debug(
            "i= {} v={}  ix {}  {}  {}", i, v, my_ix[XX], ngb_ix[XX], nx[XX]);
        // spdlog::debug("my_ix : {}", my_ix);
        // spdlog::debug("ngb_ix : {}", ngb_ix);
        spdlog::error(
            "{}: {}", "domains don't touch (X-dir)!", my_ix[XX] - ngb_ix[XX]);
      }

      //
      // Now Y-dir
      //
      if (ndim > 1) {
        if ((my_ix[YY] - ngb_ix[YY]) == 1) {  // neighbour below us
          nodelist[off1][8] = offsets[YY];
          nodelist[off1][9] = offsets[YY];
        }
        else if ((my_ix[YY] - ngb_ix[YY]) == -1) {  // neighbour above
                                                    // us
          nodelist[off1][8] = offsets[YY] + directional_Ncells[YY];
          nodelist[off1][9] = offsets[YY] + directional_Ncells[YY];
        }
        else if (my_ix[YY] == ngb_ix[YY]) {  // neighbour level with us.
          nodelist[off1][8] = offsets[YY];
          nodelist[off1][9] = offsets[YY] + directional_Ncells[YY];
        }
        else
          spdlog::error(
              "{}: {}", "domains don't touch! (Y-dir)", my_ix[YY] - ngb_ix[YY]);
      }  // at least 2D

      //
      // Now Z-dir
      //
      if (ndim > 2) {
        if ((my_ix[ZZ] - ngb_ix[ZZ]) == 1) {  // neighbour below us
          nodelist[off1][10] = offsets[ZZ];
          nodelist[off1][11] = offsets[ZZ];
        }
        else if ((my_ix[ZZ] - ngb_ix[ZZ]) == -1) {  // neighbour above
                                                    // us
          nodelist[off1][10] = offsets[ZZ] + directional_Ncells[ZZ];
          nodelist[off1][11] = offsets[ZZ] + directional_Ncells[ZZ];
        }
        else if (my_ix[ZZ] == ngb_ix[ZZ]) {  // neighbour level with us.
          nodelist[off1][10] = offsets[ZZ];
          nodelist[off1][11] = offsets[ZZ] + directional_Ncells[ZZ];
        }
        else
          spdlog::error(
              "{}: {}", "domains don't touch! (Z-dir)", my_ix[ZZ] - ngb_ix[ZZ]);
      }  // 3D

      //
      // The last 3 elements are the same for grids with same
      // orientations.
      //
      nodelist[off1][12] = 1;
      nodelist[off1][13] = 2;
      nodelist[off1][14] = 3;
    }  // loop over neighbours for this mesh
  }    // loop over meshes

  // for (int s=0;s<Stot;s++) {
  //  for (int v=0; v<nl_len; v++) {
  // cout <<"zonelist["<<s<<"]["<<v<<"] = "<<zonelist[s][v]<<"\n";
  // cout <<"nodelist["<<s<<"]["<<v<<"] = "<<nodelist[s][v]<<"\n";
  //  }
  //}

  //
  // So now we have set nmesh, meshtypes[], Nngb[], ngb[S], back[S],
  // nnodes[S], nodelist[S][], and that is all we need.
  //
  char datadir[strlength];
  err += DBGetDir(dbfile, datadir);
  DBSetDir(dbfile, "/");
  DBMkDir(dbfile, "Decomposition");
  DBSetDir(dbfile, "/Decomposition");

  DBoptlist *mma_opts = DBMakeOptlist(2);
  DBAddOption(mma_opts, DBOPT_DTIME, &SimPM.simtime);
  DBAddOption(mma_opts, DBOPT_CYCLE, &SimPM.timestep);

  //
  // write the multimeshadj
  //
  err = DBPutMultimeshadj(
      dbfile, mma_name.c_str(), nmesh, meshtypes, Nngb, ngb, back, nnodes,
      nodelist, 0, 0, mma_opts);
  if (err)
    spdlog::error(
        "{}: {}", "dataio_silo_pllel::OutputData() multimesh info", err);

  DBClearOptlist(mma_opts);
  DBFreeOptlist(mma_opts);

  //
  // Need to write a variable called "NumDomains"
  //
  int dims[ndim];
  for (int i = 0; i < ndim; i++)
    dims[i] = 1;
  DBWrite(dbfile, "NumDomains", &nmesh, dims, 1, DB_INT);

  DBSetDir(dbfile, "/");
  DBSetDir(dbfile, datadir);

  for (int s = 0; s < Stot; s++) {
    nodelist[s] = mem.myfree(nodelist[s]);
  }
#ifndef NDEBUG
  spdlog::info("Finished writing mulitmesh adjacency info.");
#endif

  return 0;
}

// ##################################################################
// ##################################################################

void dataio_silo_pllel::set_dir_in_file(
    std::string &mydir,       ///< directory name.
    const int my_rank,        ///< myrank (global).
    const int my_group_rank,  ///< myrank in group.
    const int level           ///< level in grid heirarchy
)
{
  ostringstream temp;
  temp.fill('0');
  temp.str("");
  temp << "/rank_";
  temp.width(4);
  temp << my_rank;
  // temp << "_l"; temp.width(2); temp << level;
  temp << "_domain_";
  temp.width(4);
  temp << my_group_rank;
  mydir = temp.str();
  temp.str("");
  // cout <<"\t\tdomain: "<<mydir<<"\n";
  return;
}

// ##################################################################
// ##################################################################

void dataio_silo_pllel::mesh_name(
    const int rank,  ///< rank
    string &mesh_name)
{
  //
  // Get mesh_name from rank
  //
  ostringstream temp;
  temp.str("");
  temp.fill('0');
  temp << "unigrid";
  temp.width(4);
  temp << rank;
  mesh_name = temp.str();
  return;
}

// ##################################################################
// ##################################################################

#endif  // if SILO
