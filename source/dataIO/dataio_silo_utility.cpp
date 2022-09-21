/// \file dataio_silo_utility.cpp
/// \author Jonathan Mackey
///
/// This is code for analysing silo data files in serial mode;
/// written so that a single function will determine if the file is
/// serial or parallel and read in the data regardless.  It gets the
/// parameters for the grid from the header.
///
///
///  - 2010-02-02 JM: Added support for N procs to read data written
///     by M procs, where N not equal to M.
///  - 2010-02-03 JM: Fixed all the bugs in yesterday's work (as far
///     as i could find).
///  - 2010-04-27 JM: renamed 'ngroups' to 'groupsize', and updated
///    logic so the last file can have fewer domains.
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2011.03.02 JM: Better handling of tracer variables (up to
///    MAX_NVAR now).
/// - 2012.05.17 JM: Fixed bug in how pllel_read_any_data() dealt
///    with silo databases where files don't all have the same number
///    of meshes.
/// - 2013.02.19 JM: Got rid of dataio_utility class, and moved its
///    functions into file_status class, which now has its own file.
///    Renamed file to dataio_silo_utility.cpp from dataio_utility.cc
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.03.26 JM: updated for pion v0.2
/// - 2015.06.13 JM: started updating to work with void* pointers for
///    coordinates and data variables (they can be float or double).
/// - 2016.03.18 JM: updated to work better with large grids and with
///    FLOAT and DOUBLE data (no buggy integer positions anymore...).
/// - 2016.04.04 JM: fixed bug in get_quadmesh_integer_extents() for
///    xmin/xmax<0.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"

#include "tools/timer.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#ifndef NDEBUG
#endif  // NDEBUG

#include "dataIO/dataio_silo_utility.h"

#include <sstream>
using namespace std;

//#define TEST_SILO_IO

/********************************************************/
/*************** dataio_silo_utility ********************/
/********************************************************/

dataio_silo_utility::dataio_silo_utility(
    class SimParams &SimPM,  ///< pointer to simulation parameters
    std::string dtype,       ///< FLOAT or DOUBLE for files.
    class Sub_domain *p) :
    dataio_silo_pllel(SimPM, dtype, p)
{
#ifdef TEST_SILO_IO
  spdlog::info("Setting up utility Silo I/O class");
#endif
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::SRAD_get_nproc_numfiles(string fname, int *np, int *nf)
{
  int err = 0;

  //
  // open file
  //
#ifdef TEST_SILO_IO
  spdlog::debug("opening file: {}", fname);
#endif
  DBfile *dbfile = 0;
  dbfile         = DBOpen(fname.c_str(), DB_UNKNOWN, DB_READ);
  if (!dbfile)
    spdlog::error("{}: {}", "open silo file failed.", fmt::ptr(dbfile));

  //
  // read nproc, numfiles from header
  //
  DBSetDir(dbfile, "/header");
  int nproc = 0, numfiles = 0;  //, groupsize=0;
  err += DBReadVar(dbfile, "MPI_nproc", &nproc);
  err += DBReadVar(dbfile, "NUM_FILES", &numfiles);
  if (err) {
#ifdef TEST_SILO_IO
    spdlog::error(
        "must be serial file -- failed to find NUM_FILES and MPI_nproc")
        spdlog::error("continuing assuming serial file....");
#endif
    // spdlog::error("{}: {}", "error reading params from file",fname);
    nproc    = 1;
    numfiles = 1;  // groupsize=1;
  }
  else {
#ifdef TEST_SILO_IO
    spdlog::debug("\tRead nproc={}\tand numfiles={}", nproc, numfiles);
#endif
    //    groupsize = nproc/numfiles;
    // nproc = nproc;
  }
  DBClose(dbfile);
  dbfile = 0;

  string::size_type pos = fname.find("_0000");
  if (pos == string::npos) {
    spdlog::warn("didn't find _0000 in file, so we are reading serial file");
    err = 1;
  }

  *np = nproc;
  *nf = numfiles;
#ifdef TEST_SILO_IO
  spdlog::debug(
      "dataio_silo_utility::SRAD_get_nproc_numfiles returning {}", err);
#endif
  return err;
}

// ##################################################################
// ##################################################################

bool dataio_silo_utility::SRAD_point_on_my_domain(
    const cell &c,            ///< pointer to cell
    class SimParams &SimPM,   ///< pointer to simulation parameters
    class Sub_domain *filePM  ///< pointer to class with nproc.
)
{
  //
  // Assume point is on domain, and set to false if found to be off.
  //
  bool on = true;
  std::array<double, MAX_DIM> dpos;
  CI.get_dpos(c, dpos);
  for (int i = 0; i < SimPM.ndim; i++) {
    if (dpos[i] < filePM->get_Xmin(i)) on = false;
    if (dpos[i] > filePM->get_Xmax(i)) on = false;
  }
  return on;
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::SRAD_read_var2grid(
    DBfile *dbfile,            ///< pointer to silo file.
    class GridBaseClass *ggg,  ///< pointer to data.
    const string variable,     ///< variable name to read.
    const long int npt,        ///< number of cells expected.
    class SimParams &SimPM,    ///< pointer to simulation parameters
    class Sub_domain *filePM   ///< pointer to class with nproc.
)
{
  //
  // The dbfile pointer should already be in the directory containing
  // the named variable to read, so it's ok to bug out if we don't
  // find it.
  //
  DBquadvar *silodata = 0;
  silodata            = DBGetQuadvar(dbfile, variable.c_str());
  if (!silodata) {
    spdlog::error(
        "{}: {}", "dataio_silo::read_variable2grid() failed to read variable",
        variable);
  }
  if (silodata->nels != npt) {
    spdlog::error(
        "{}: {}", "dataio_silo::read_variable2grid() wrong number of cells",
        silodata->nels - SimPM.Ncell);
  }

  //
  // Check that datatype is what we are expecting!  If not, then
  // delete data arrays, reset datatype, and re-create data arrays.
  //
  if (silodata->datatype != silo_dtype) {
#ifdef TEST_SILO_IO
    spdlog::debug(
        "\n\tSRAD_read_var2grid() quadvar has type={} but expecting type={}\n\t... resetting datatype for this file.\n    DB_INT=16, DB_SHORT=17, DB_LONG=18, DB_FLOAT=19, DB_DOUBLE=20, DB_CHAR=21, DB_LONG_LONG=22, DB_NOTYPE=25",
        silodata->datatype, silo_dtype);
#endif
    delete_data_arrays();
    silo_dtype = silodata->datatype;
    create_data_arrays(SimPM);
  }

  //
  // Create a pointer to the data in the silo stuct DBquadvar.  This
  // is a void pointer, so I have to reinterpret it to get data that
  // PION can understand.
  //
  void **data    = silodata->vals;
  float **fdata  = 0;
  double **ddata = 0;
  if (silo_dtype == DB_FLOAT) {
    fdata = reinterpret_cast<float **>(data);
  }
  else {
    ddata = reinterpret_cast<double **>(data);
  }

  if (variable == "Velocity" || variable == "MagneticField") {
    int v1, v2, v3;
    bool B      = false;
    double norm = 1.0 / sqrt(4.0 * M_PI);
    if (variable == "Velocity") {
      v1 = VX;
      v2 = VY;
      v3 = VZ;
    }
    else {
      v1 = BX;
      v2 = BY;
      v3 = BZ;
      B  = true;
    }
    //    cout <<"name: "<<silodata->name<<"\tnels="<<silodata->nels<<"\n";
    //    cout <<"ndims: "<<silodata->ndims<<"\tnvals:
    //    "<<silodata->nvals<<"\n";
    // cout <<"reading variable "<<variable<<" into element "<<v1<<" of
    // state vec.\n";
    cell *c     = ggg->FirstPt();
    long int ct = 0;
    do {
      if (SRAD_point_on_my_domain(*c, SimPM, filePM)) {
        //      cout <<"ct="<<ct<<"\t and ncell="<<npt<<"\n";
        if (silo_dtype == DB_FLOAT) {
          c->P[v1] = fdata[0][ct];
          c->P[v2] = fdata[1][ct];
          c->P[v3] = fdata[2][ct];
          ct++;
        }
        else {
          c->P[v1] = ddata[0][ct];
          c->P[v2] = ddata[1][ct];
          c->P[v3] = ddata[2][ct];
          ct++;
        }

#ifdef NEW_B_NORM
        if (B) {
          // scale values from CGS to code units.
          c->P[v1] *= norm;
          c->P[v2] *= norm;
          c->P[v3] *= norm;
        }
#endif
      }
    } while ((c = ggg->NextPt(*c)) != 0);
    if (ct != npt)
      spdlog::error(
          "{}: {}", "wrong number of points read for vector variable",
          ct - npt);
  }  // vector variable

  else {
    int v1 = 0;
    bool B = false;
    if (variable == "Density")
      v1 = RO;
    else if (variable == "Pressure")
      v1 = PG;
    else if (variable == "VelocityX")
      v1 = VX;
    else if (variable == "VelocityY")
      v1 = VY;
    else if (variable == "VelocityZ")
      v1 = VZ;
    else if (variable == "MagneticFieldX") {
      v1 = BX;
      B  = true;
    }
    else if (variable == "MagneticFieldY") {
      v1 = BY;
      B  = true;
    }
    else if (variable == "MagneticFieldZ") {
      v1 = BZ;
      B  = true;
    }
    else if (variable == "glmPSI")
      v1 = SI;
    //
    // Now loop over up to MAX_NVAR tracers...
    //
    else if (variable.substr(0, 2) == "Tr") {
      int itr = atoi(variable.substr(2, 3).c_str());
      if (!isfinite(itr) || itr < 0 || itr >= MAX_NVAR) {
        spdlog::error(
            "{}: {}", "Bad diffuse Column-density identifier.", variable);
      }
      v1 = SimPM.ftr + itr;
    }
    else
      spdlog::error("{}: {}", "what var to read???", variable);
    // cout <<"reading variable "<<variable<<" into element "<<v1<<" of
    // state vec.\n";

    //
    // First get to start cell in local domain:
    //
    enum direction posdir[3] = {XP, YP, ZP};
    cell *start              = ggg->FirstPt();
    long int ct              = 0;
    for (int i = 0; i < SimPM.ndim; i++) {
      while (CI.get_dpos(*start, i) < filePM->get_Xmin(i))
        start = ggg->NextPt(*start, posdir[i]);
    }
    //
    // Now use NG for-loops to only go through cells in the local domain.
    // We assume we go through x-dir column first, then along y, and finally
    // Z
    //
    class cell *cx = start, *cy = start, *cz = start;
    double norm = 1.0 / sqrt(4.0 * M_PI);
    if (SimPM.ndim < 3) filePM->set_directional_Ncells(ZZ, 1);
    if (SimPM.ndim < 2) filePM->set_directional_Ncells(YY, 1);
    for (int k = 0; k < filePM->get_directional_Ncells(ZZ); k++) {
      cy = cz;
      for (int j = 0; j < filePM->get_directional_Ncells(YY); j++) {
        cx = cy;
        for (int i = 0; i < filePM->get_directional_Ncells(XX); i++) {
          if (!SRAD_point_on_my_domain(*cx, SimPM, filePM))
            spdlog::error(
                "{}: {}", "FAST READ IS IN THE WRONG PLACE!!!", cx->pos[XX]);
          if (silo_dtype == DB_FLOAT) {
            cx->P[v1] = fdata[0][ct];
          }
          else {
            cx->P[v1] = ddata[0][ct];
          }
#ifdef NEW_B_NORM
          // scale values from CGS to code units.
          if (B) {
            // cout <<"SRAD_read_var2grid: scaling B var "<<v1;
            // cout <<" by "<<norm<<"\n";
            cx->P[v1] *= norm;
          }
#endif
          ct++;
          cx = ggg->NextPt(*cx, posdir[XX]);
        }
        if (SimPM.ndim > 1) cy = ggg->NextPt(*cy, posdir[YY]);
      }
      if (SimPM.ndim > 2) cz = ggg->NextPt(*cz, posdir[ZZ]);
    }
    if (ct != npt) {
      spdlog::error(
          "{}: {}", "wrong number of points read for scalar variable",
          ct - npt);
    }

  }  // scalar variable

  //  cout <<"Read variable "<<variable<<"\n";
  DBFreeQuadvar(silodata);  // silodata=0;
  fdata = 0;
  ddata = 0;
  return 0;
}

// ##################################################################
// ##################################################################

void dataio_silo_utility::set_pllel_filename(
    std::string &fname,  ///< filename
    const int ifile      ///< file number to replace name with.
)
{
  //
  // If we have a parallel file, Parse filename, and replace file number
  // with new ifile, store in name 'fname'.
  //
  ostringstream temp;
  temp.fill('0');
  string::size_type pos = fname.find("_0000");
  if (pos == string::npos) {
    spdlog::error(
        "didn't find _0000 in file, but claim we are reading pllel file!");
    spdlog::error("{}: {}", "not a parallel i/o filename!", fname);
  }
  else {
    temp.str("");
    temp << "_";
    temp.width(4);
    temp << ifile;
    fname.replace(pos, 5, temp.str());
    // cout <<"\tNew fname: "<<fname<<"\n";
    temp.str("");
  }
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::serial_read_pllel_silodata(
    const string firstfile,    ///< filename
    class GridBaseClass *ggg,  ///< pointer to data.
    const int numfiles,        ///< number of files
    const int groupsize,       ///< number of groups
    class SimParams &SimPM,    ///< pointer to simulation parameters
    class Sub_domain *filePM   ///< number of processes used to write file.
)
{
  int err = 0;

  int level = 0;
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    if (ggg == SimPM.levels[l].grid) level = l;
    // cout <<"saving level "<<level<<"\n";
  }

  //
  // First loop over all files:
  //
  for (int ifile = 0; ifile < numfiles; ifile++) {
    string infile = firstfile;
    //
    // Replace filename with new file corresponding to current 'ifile'
    // value.
    //
    set_pllel_filename(infile, ifile);

    DBfile *dbfile = DBOpen(infile.c_str(), DB_UNKNOWN, DB_READ);
    if (!dbfile)
      spdlog::error("{}: {}", "open first silo file failed.", fmt::ptr(dbfile));
    //
    // loop over domains within this file.  The last file may have fewer
    // domains, so we set ng to be the minimum of groupsize or all
    // remaining domains.
    //
    int ng = min(filePM->get_nproc() - ifile * groupsize, groupsize);
    for (int igroup = 0; igroup < ng; igroup++) {
      DBSetDir(dbfile, "/");
      //
      // choose myrank, and decompose domain accordingly.
      //
      filePM->set_myrank(ifile * groupsize + igroup);
      filePM->decomposeDomain(SimPM.ndim, SimPM.levels[level]);

      //
      // set directory in file.
      //
      string mydir;
      set_dir_in_file(mydir, filePM->get_myrank(), igroup, level);
      DBSetDir(dbfile, mydir.c_str());

      //
      // set variables to read: (from dataio_silo class)
      //
      set_readvars(SimPM);

      //
      // now read each variable in turn from the mesh
      //
      for (std::vector<string>::iterator i = readvars.begin();
           i != readvars.end(); ++i) {
        err = SRAD_read_var2grid(
            dbfile, ggg, (*i), filePM->get_Ncell(), SimPM, filePM);
        if (err) spdlog::error("{}: {}", "error reading variable", (*i));
      }
    }  // loop over domains within a file.

    //
    // Close this file
    //
    DBClose(dbfile);
    dbfile = 0;
  }  // loop over files

  // cout <<"read parallel data successfully.\n";
  return 0;
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::ReadData(
    string firstfile,                   ///< file to read from
    vector<class GridBaseClass *> &cg,  ///< grid pointers.
    class SimParams &SimPM              ///< pointer to simulation parameters
)
{
  spdlog::debug("(pion-ng mpi)  Reading data from file: {}", firstfile);
  silofile = firstfile;

  int err = 0;
  // Loop over grid levels, and read data for each level.
  for (int l = 0; l < SimPM.grid_nlevels; l++) {

    spdlog::debug(
        " reading data on level {}, nlevels={}", l, SimPM.grid_nlevels);
    // for now read a different file for each level in the NG grid.
    // If more than one level of grid, look for level in filename:
    string::size_type p;
    if ((p = silofile.find("_level")) == string::npos
        && SimPM.grid_nlevels > 1) {
      spdlog::error(
          "{}: {}", "dataio_silo_utility::ReadData() level", silofile);
    }
    else if (SimPM.grid_nlevels > 1) {
      ostringstream temp;
      temp.str("");
      temp.width(2);
      temp.fill('0');
      temp << l;
      silofile.replace(p + 6, 2, temp.str());
#ifdef TEST_SILO_IO
      spdlog::debug("p={}  string={}, silofile={}", p, temp.str(), silofile);
#endif
    }

    if (!cg[l])
      spdlog::error(
          "{}: {}", "dataio_silo_utility::ReadData() null grid!",
          fmt::ptr(cg[l]));
    dataio_silo::gp = cg[l];
    mpiPM           = &(SimPM.levels[l].sub_domain);

    //
    // set grid properties for quadmesh: each level of the NG
    // grid has different zone and node coordinates so we need to
    // call this each time.
    //
    // err = setup_grid_properties(gp, SimPM);
    // rep.errorTest("IO_silo_utility::ReadData() setup_grid_properties"
    //                                                        ,0, err);

    // now call the code that works for each level:
    err = ReadLevelData(silofile, gp, SimPM, l);
    if (0 != err) {
      spdlog::error(
          "{}: Expected {} but got {}", "IO_silo_utility:: ReadLevelData", 0,
          err);
      exit(1);
    }
#ifdef TEST_SILO_IO
    spdlog::debug(
        "Finished reading data on level {} of {}", l, SimPM.grid_nlevels);
#endif
  }

  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::ReadLevelData(
    string firstfile,         ///< file to read from
    class GridBaseClass *cg,  ///< grid pointer.
    class SimParams &SimPM,   ///< simulation parameters
    const int l               ///< level in grid hierarchy
)
{

  class GridBaseClass *ggg = cg;
  if (!ggg)
    spdlog::error(
        "{}: {}", "null pointer to computational grid!", fmt::ptr(ggg));
  int err = 0;

  //
  // The idea behind this is that a program running on N cores can
  // read data written by M cores, where N and M can be any positive
  // integers (possibly powers of 2 for the domain decomposition to
  // work).
  //
  // If we are a parallel program, then mpiPM should be already set,
  // and the domain decomposition into N sub-domains is already done,
  // and the grid is set up for this sub-domain.
  //
  int nproc = 0, numfiles = 0, groupsize = 0;
  err = SRAD_get_nproc_numfiles(firstfile, &nproc, &numfiles);
  if (err) {
    //
    // must be reading serial file, so we only want part of the domain
    // read onto the local grid, so we need a new function to read
    // this.
    //
    groupsize = 1;
    numfiles  = 1;
    err       = parallel_read_serial_silodata(firstfile, ggg, SimPM);
    if (0 != err)
      spdlog::error(
          "{}: Expected {} but got {}", "(silocompare) Failed to read data", 0,
          err);
  }
  else {
    //
    // must be reading parallel file, so want to read in every
    // subdomain onto grid.  Use local functions for this:
    //
    // groupsize is decided interestingly by PMPIO (silo): if numfiles
    // divided nproc evenly then groupsize is nproc/numfiles, but if not
    // then groupsize is
    // ((int) (nproc/numfiles))+1 for the first NX domains, and
    // (nproc/numfiles) for the remainder.  Here NX=(nproc%numfiles)
    //
#ifdef TEST_SILO_IO
    spdlog::debug("READING PLLEL DATA for level {}", l);
#endif
    groupsize = nproc / numfiles;
    if ((nproc % numfiles) != 0) groupsize++;

    //
    // We should take it in turns reading data.
    // Try allowing up to 16 simultaneous reads.
    //
    int max_reads = 16;
    int nloops    = 0;
    if (mpiPM->get_nproc() < max_reads)
      nloops = 1;
    else {
      nloops = mpiPM->get_nproc() / max_reads;
      if (mpiPM->get_nproc() % max_reads != 0) {
        nloops += 1;  // this shouldn't happen, but anyway...
        spdlog::debug("Nproc not a power of 2!  This will cause trouble.\n");
        spdlog::error(
            "{}: {}", "dataio_silo_utility::ReadLevelData()",
            mpiPM->get_nproc());
      }
    }

    clk.start_timer("readdata");
    double tsf = 0;
#ifdef TEST_SILO_IO
    spdlog::debug(
        "READING PLLEL DATA: {}  {}  {} {} {}", mpiPM->get_myrank(),
        mpiPM->get_nproc(), numfiles, groupsize, l);
#endif
    std::vector<int> cells_per_file;
    mpiPM->allgather_ncells(cells_per_file);
    for (int count = 0; count < nloops; count++) {
      if ((mpiPM->get_myrank() + nloops) % nloops == count) {
#ifdef TEST_SILO_IO
        spdlog::debug(
            "!READING DATA!!... myrank={}  i={}", mpiPM->get_myrank(), count);
#endif
        err = parallel_read_parallel_silodata(
            firstfile, ggg, SimPM, numfiles, groupsize, nproc, l,
            cells_per_file);
        if (0 != err) {
          spdlog::error(
              "{}: Expected {} but got {}", "Failed to read parallel data", 0,
              err);
          exit(1);
        }
      }
      else {
#ifdef TEST_SILO_IO
        spdlog::debug(
            "waiting my turn... myrank={}  i={}", mpiPM->get_myrank(), count);
#endif
      }
      mpiPM->barrier();
      tsf = clk.time_so_far("readdata");
#ifdef TEST_SILO_IO
      spdlog::debug("\t time = {} secs.", tsf);
#endif
    }
    clk.stop_timer("readdata");
  }
#ifdef TEST_SILO_IO
  spdlog::info("read data successfully");
#endif
  return 0;
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::parallel_read_serial_silodata(
    string firstfile,          ///< file to read from
    class GridBaseClass *ggg,  ///< pointer to data.
    class SimParams &SimPM     ///< pointer to simulation parameters
)
{
  int err = 0;
  //
  // This is quite simple -- the local domain reads in a subset of the
  // uniform grid in the file.
  //
  DBfile *dbfile = DBOpen(firstfile.c_str(), DB_UNKNOWN, DB_READ);
  if (!dbfile)
    spdlog::error("{}: {}", "open first silo file failed.", fmt::ptr(dbfile));

  int ftype = DBGetDriverType(dbfile);
  if (ftype == DB_HDF5) {
    // cout <<"READING HDF5 FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    //
    // Not sure if we need this for reading, but set it anyway.
    //
    int friendly = DBGetFriendlyHDF5Names();
    DBSetFriendlyHDF5Names(friendly);
  }

  //
  // Set variables to read based on what equations we are using (this
  // is read from the header previously)
  //
  err = set_readvars(SimPM);
  if (err) spdlog::error("{}: {}", "failed to set readvars in ReadData", err);

  string qm_dir("/rank_0000_domain_0000");
  DBSetDir(dbfile, qm_dir.c_str());
  string qm_name = "unigrid0000";

  //
  // Get max and min quadmesh positions (integers)
  //
  std::array<int, MAX_DIM> mesh_iXmin, mesh_iXmax;
  get_quadmesh_integer_extents(
      dbfile, ggg, SimPM, qm_dir, qm_name, mesh_iXmin, mesh_iXmax);
  // if (err) spdlog::error("{}: {}", "Failed to get quadmesh extents!",qm_dir);

  //
  // now read each variable in turn from the mesh, using the parallel
  // read function.
  //
  for (std::vector<string>::iterator i = readvars.begin(); i != readvars.end();
       ++i) {
    err = PP_read_var2grid(
        dbfile, ggg, SimPM, (*i), SimPM.Ncell, &mesh_iXmin[0], &mesh_iXmax[0]);
    if (err)
      spdlog::error(
          "{}: {}", "dataio_silo::ReadData() error reading variable", (*i));
  }

  DBSetDir(dbfile, "/");
  DBClose(dbfile);
  dbfile = 0;

  return err;
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::parallel_read_parallel_silodata(
    string firstfile,                       ///< file to read from
    class GridBaseClass *ggg,               ///< pointer to data.
    class SimParams &SimPM,                 ///< simulation parameters
    const int numfiles,                     ///< number of files
    const int groupsize,                    ///< number of groups
    const int nmesh,                        ///< number of domains in file.
    const int l,                            ///< level in grid hierarchy
    const std::vector<int> &cells_per_file  ///< number of cells in each file
)
{
  int err = 0;

  //
  // First loop over all files:
  //
  for (int ifile = 0; ifile < numfiles; ifile++) {
    string infile = firstfile;
    //
    // Replace filename with new file corresponding to current 'ifile'
    // value.
    //
    set_pllel_filename(infile, ifile);

    DBfile *dbfile = DBOpen(infile.c_str(), DB_UNKNOWN, DB_READ);
    if (!dbfile)
      spdlog::error("{}: {}", "open first silo file failed.", fmt::ptr(dbfile));
    //
    // loop over domains within this file. The number of domains per file is
    // not constant if (R=nproc%numfiles)!=0.  One extra domain is added to
    // the first R files.  If R==0, then we can ignore this, since all files
    // have ng=groupsize domains
    //
    int R = nmesh % numfiles;
    if (R < 0) R += numfiles;
    int ng = 0;
    if (R != 0) {
      if (ifile < R)
        ng = groupsize;
      else
        ng = groupsize - 1;
    }
    else {
      ng = groupsize;
    }
    // int ng=min(nmesh-ifile*groupsize, groupsize);

    for (int igroup = 0; igroup < ng; igroup++) {
#ifdef TEST_SILO_IO
      spdlog::info("igroup {}, ngroup {}", igroup, ng);
#endif
      DBSetDir(dbfile, "/");
      //
      // choose pseudo-rank for reading file.
      //
      int pseudo_rank;
      if (R != 0) {
        if (ifile < R)
          pseudo_rank = ifile * groupsize + igroup;
        else
          pseudo_rank =
              (R * groupsize) + (ifile - R) * (groupsize - 1) + igroup;
      }
      else {
        pseudo_rank = ifile * groupsize + igroup;
      }

      // mpiPM = &(filePM);
      // err = setup_grid_properties(gp, SimPM);
      // rep.errorTest("IO_silo_utility::p_r_pd() setup_grid_properties"
      //                                                       ,0,err);

      //
      // set directory in file.
      //
      string qm_dir;
      set_dir_in_file(qm_dir, pseudo_rank, igroup, l);
      // cout << "DBSetDir 2:" << qm_dir.c_str() << std::endl;
      DBSetDir(dbfile, qm_dir.c_str());

      //
      // Set mesh_name from rank. quadmesh is in the current directory
      // and called "unigridXXXX" where XXXX=filePM.myrank
      //
      string qm_name;
      mesh_name(pseudo_rank, qm_name);
#ifdef TEST_SILO_IO
      spdlog::debug("got mesh name= {} in mesh dir= {}", qm_name, qm_dir);
#endif

      //
      // Get max and min quadmesh positions (integers)
      //
      std::array<int, MAX_DIM> mesh_iXmin, mesh_iXmax;
      get_quadmesh_integer_extents(
          dbfile, ggg, SimPM, qm_dir, qm_name, mesh_iXmin, mesh_iXmax);
      if (err)
        spdlog::error("{}: {}", "Failed to get quadmesh extents!", qm_dir);

      //
      // Get max and min grid positions (integers)
      //
      std::vector<int> localmin(ndim), localmax(ndim);
      for (int v = 0; v < ndim; v++) {
        localmin[v] = ggg->iXmin(static_cast<axes>(v));
        localmax[v] = ggg->iXmax(static_cast<axes>(v));
      }

      //
      // See if quadmesh intersects local domain at all.
      //
      bool get_data = true;
      for (int v = 0; v < ndim; v++) {
        if ((mesh_iXmax[v] <= localmin[v]) || (mesh_iXmin[v] >= localmax[v]))
          get_data = false;
      }

      if (!get_data) {
#ifdef TEST_SILO_IO
        spdlog::info(
            "*** skipping mesh {} because not on local domain", qm_name);
#endif
      }
      else {
#ifdef TEST_SILO_IO
        spdlog::info(
            "**** reading mesh {} because it is on local domain", qm_name);
        spdlog::info("mesh_iXmin : {}", mesh_iXmin);
        spdlog::info("mesh_iXmax : {}", mesh_iXmax);
        spdlog::info("local iXmin : {}", localmin);
        spdlog::info("local iXmax : {}", localmax);
#endif
        //
        // set variables to read: (from dataio_silo class)
        //
        set_readvars(SimPM);

        //
        // now read each variable in turn from the mesh, using the
        // parallel-parallel read function
        //
        for (std::vector<string>::iterator i = readvars.begin();
             i != readvars.end(); ++i) {
#ifdef TEST_SILO_IO
          spdlog::info(
              "reading var {} to grid, pdeudorank {}, cells-per-file {}", *i,
              pseudo_rank, cells_per_file);
#endif
          err = PP_read_var2grid(
              dbfile, ggg, SimPM, (*i), cells_per_file[pseudo_rank],
              &mesh_iXmin[0], &mesh_iXmax[0]);
          if (err) spdlog::error("{}: {}", "error reading variable", (*i));
        }  // loop over variables.
      }
    }  // loop over domains within a file.

#ifdef TEST_SILO_IO
    spdlog::info("finished reading domains, parallel_read_parallel_silodata");
#endif
    DBClose(dbfile);
    dbfile = 0;
  }  // loop over files

  return 0;
}

// ##################################################################
// ##################################################################

void dataio_silo_utility::get_quadmesh_extents(
    DBfile *dbfile,         ///< pointer to silo file.
    const string mesh_dir,  ///< directory of mesh
    const string qm_name,   ///< name of mesh
    double *mesh_xmin,      ///< Xmin for mesh (output)
    double *mesh_xmax,      ///< Xmax for mesh (output)
    class SimParams &SimPM  ///< pointer to simulation parameters
)
{

  //
  // Make sure we're in the right dir
  //
  DBSetDir(dbfile, mesh_dir.c_str());

  //
  // Now get the mesh from the file.
  //
  DBquadmesh *qm = 0;
  qm             = DBGetQuadmesh(dbfile, qm_name.c_str());
  if (!qm) spdlog::error("{}: {}", "failed to get quadmesh", qm_name);

  //
  // Check that datatype is what we are expecting!  If not, then
  // delete data arrays, reset datatype, and re-create data arrays.
  //
  if (qm->datatype != silo_dtype) {
#ifdef TEST_SILO_IO
    spdlog::debug(
        "\n\tget_quadmesh_extents() quadvar has type={} but expecting type={}\n\t... resetting datatype for this file.\n    DB_INT=16, DB_SHORT=17, DB_LONG=18, DB_FLOAT=19, DB_DOUBLE=20, DB_CHAR=21, DB_LONG_LONG=22, DB_NOTYPE=25",
        qm->datatype, silo_dtype);
#endif
    delete_data_arrays();
    silo_dtype = qm->datatype;
    create_data_arrays(SimPM);
  }

  if (silo_dtype == DB_FLOAT) {
    float *fqmmin = qm->min_extents;
    float *fqmmax = qm->max_extents;
    for (int v = 0; v < ndim; v++) {
      mesh_xmin[v] = fqmmin[v];
      mesh_xmax[v] = fqmmax[v];
#ifndef NDEBUG
      spdlog::debug(
          "dir: {}\t min={} and max={}", v, mesh_xmin[v], mesh_xmax[v]);
#endif  // NDEBUG
    }
  }
  else {
    double *dqmmin = reinterpret_cast<double *>(qm->min_extents);
    double *dqmmax = reinterpret_cast<double *>(qm->max_extents);
    for (int v = 0; v < ndim; v++) {
      mesh_xmin[v] = dqmmin[v];
      mesh_xmax[v] = dqmmax[v];
#ifndef NDEBUG
      spdlog::debug(
          "dir: {}\t min={} and max={}", v, mesh_xmin[v], mesh_xmax[v]);
#endif  // NDEBUG
    }
  }
  DBFreeQuadmesh(qm);  // qm=0;
  return;
}

// ##################################################################
// ##################################################################

void dataio_silo_utility::get_quadmesh_integer_extents(
    DBfile *dbfile,                   ///< pointer to silo file.
    class GridBaseClass *ggg,         ///< pointer to data.
    class SimParams &SimPM,           ///< pointer to simulation parameters
    const string mesh_dir,            ///< directory of mesh
    const string qm_name,             ///< name of mesh
    std::array<int, MAX_DIM> &iXmin,  ///< integer Xmin for mesh (output)
    std::array<int, MAX_DIM> &iXmax   ///< integer Xmax for mesh (output)
)
{
  //
  // First get the double precision extents
  //
  std::array<double, MAX_DIM> mesh_xmin, mesh_xmax;
  get_quadmesh_extents(
      dbfile, mesh_dir, qm_name, &mesh_xmin[0], &mesh_xmax[0], SimPM);

  //
  // Now use the cell interface to get the integer extents (Note that
  // this will fail and bug out if the global grid class isn't set up,
  // since that defines the coordinate system).
  //
  if (silo_dtype == DB_FLOAT) {
    //
    // we need an extra buffer here to put the variable on the +ve
    // side of the cell border.  Have to do it a bit carefully,
    // because the xmin/xmax values can be > or < 0.
    //
    double buffer = 0.0;
    for (int v = 0; v < ndim; v++) {
      buffer = std::min(1.0e-5, 0.1 / ggg->SIM_iRange(static_cast<axes>(v)));

      if (mesh_xmin[v] < 0.0)
        mesh_xmin[v] *= (1.0 - buffer);
      else
        mesh_xmin[v] *= (1.0 + buffer);

      if (mesh_xmax[v] < 0.0)
        mesh_xmax[v] *= (1.0 - buffer);
      else
        mesh_xmax[v] *= (1.0 + buffer);
    }
  }
  for (int v = 0; v < MAX_DIM; v++) {
    iXmin[v] = 0;
    iXmax[v] = 0;
  }
  CI.get_ipos_vec(mesh_xmin, iXmin);
  CI.get_ipos_vec(mesh_xmax, iXmax);

#ifndef NDEBUG
  spdlog::debug("get_quadmesh_integer_extents: mesh_Xmin : {}", mesh_xmin);
  spdlog::debug("get_quadmesh_integer_extents: mesh_Xmax : {}", mesh_xmax);
  spdlog::debug("get_quadmesh_integer_extents: iXmin : {}", iXmin);
  spdlog::debug("get_quadmesh_integer_extents: iXmax : {}", iXmax);
#endif  // NDEBUG
}

// ##################################################################
// ##################################################################

int dataio_silo_utility::PP_read_var2grid(
    DBfile *dbfile,            ///< pointer to silo file.
    class GridBaseClass *ggg,  ///< pointer to data.
    class SimParams &SimPM,    ///< pointer to simulation parameters
    const string variable,     ///< variable name to read.
    const long int,            ///< number of cells expected (defunct!)
    const int *iXmin,          ///< integer Xmin for mesh
    const int *iXmax           ///< integer Xmax for mesh
)
{
  //
  // The dbfile pointer should already be in the directory containing
  // the named variable to read, so it's ok to bug out if we don't
  // find it.
  //
  DBquadvar *qv = 0;
  qv            = DBGetQuadvar(dbfile, variable.c_str());
  if (!qv)
    spdlog::error(
        "{}: {}", "dataio_silo::read_variable2grid() failed to read variable",
        variable);

  //
  // Check that datatype is what we are expecting!  If not, then
  // delete data arrays, reset datatype, and re-create data arrays.
  //
  if (qv->datatype != silo_dtype) {
#ifdef TEST_SILO_IO
    spdlog::debug(
        "\n\tPP_read_var2grid() quadvar has type={} but expecting type={}\n\t... resetting datatype for this file.\n    DB_INT=16, DB_SHORT=17, DB_LONG=18, DB_FLOAT=19, DB_DOUBLE=20, DB_CHAR=21, DB_LONG_LONG=22, DB_NOTYPE=25",
        qv->datatype, silo_dtype);
#endif
    delete_data_arrays();
    silo_dtype = qv->datatype;
    create_data_arrays(SimPM);
  }

  //
  // So now part of the quadmesh intersects the local domain, so we
  // run through the data and pick out the ones we want.  Silo stores
  // the data in a big 1D array with elements stored in the order
  // D[NY*NY*iz+NX*iy+ix], so we go along x-columns, then y, then z.
  //
  //
  // Set a pointer to the data in the quadmesh
  //
  void **data    = qv->vals;
  float **fdata  = 0;
  double **ddata = 0;
  if (silo_dtype == DB_FLOAT) {
    fdata = reinterpret_cast<float **>(data);
  }
  else {
    ddata = reinterpret_cast<double **>(data);
  }

  //
  // Set variables to read, first check for vector and then scalar
  // data
  //
  int v1 = -1, v2 = -1, v3 = -1;
  bool B = false;
  if (variable == "Velocity" || variable == "MagneticField") {
    if (variable == "Velocity") {
      v1 = VX;
      v2 = VY;
      v3 = VZ;
    }
    else {
      v1 = BX;
      v2 = BY;
      v3 = BZ;
      B  = true;
    }
  }
  else {
    if (variable == "Density")
      v1 = RO;
    else if (variable == "Pressure")
      v1 = PG;
    else if (variable == "VelocityX")
      v1 = VX;
    else if (variable == "VelocityY")
      v1 = VY;
    else if (variable == "VelocityZ")
      v1 = VZ;
    else if (variable == "MagneticFieldX") {
      v1 = BX;
      B  = true;
    }
    else if (variable == "MagneticFieldY") {
      v1 = BY;
      B  = true;
    }
    else if (variable == "MagneticFieldZ") {
      v1 = BZ;
      B  = true;
    }
    else if (variable == "glmPSI")
      v1 = SI;
    //
    // Now loop over up to MAX_NVAR tracers...
    //
    else if (variable.substr(0, 2) == "Tr") {
      int itr = atoi(variable.substr(2, 3).c_str());
      if (!isfinite(itr) || itr < 0 || itr >= MAX_NVAR) {
        spdlog::error(
            "{}: {}", "Bad diffuse Column-density identifier.", variable);
      }
      v1 = SimPM.ftr + itr;
    }
    else
      spdlog::error("{}: {}", "what var to read???", variable);
  }

  //
  // Get to first cell in the quadmesh/grid intersection region.
  //
  cell *c = ggg->FirstPt();

  for (int v = 0; v < ndim; v++) {
#ifndef NDEBUG
    spdlog::debug("c->pos[v]={} : iXmin[v]={}", c->pos[v], iXmin[v]);
#endif
    enum direction posdir = static_cast<direction>(2 * v + 1);

    while (c != 0 && c->pos[v] < iXmin[v]) {
      c = ggg->NextPt(*c, posdir);
    }
    if (!c) {
      spdlog::error(
          "{}: {}", "Went off end of grid looking for starting cell",
          iXmin[v] - ggg->FirstPt()->pos[v]);
    }
  }

  //
  // Get the starting indices for the quadmesh -- the starting cell is
  // at least dx/2 greater than iXmin in every direction, and the
  // difference is an odd number for cell-size=2, so this integer
  // division will always get the right answer, since all compilers
  // round down for integer division of positive numbers.
  //
  int dx = ggg->idx();
  int qm_start[ndim];
  int qm_ix[ndim], qm_NX[ndim];
  for (int v = 0; v < ndim; v++) {
    qm_start[v] = (c->pos[v] - iXmin[v]) / dx;
#ifdef TEST_SILO_IO
    spdlog::debug(
        "\t\tv={} start={} pos={}, xmin={} dims={}, var = {}", v, qm_start[v],
        c->pos[v], iXmin[v], qv->dims[v], variable);
#endif
    qm_ix[v] = qm_start[v];
    //
    // Get number of elements in each direction for this subdomain.
    // Can read it from the quadvar struct or else we could get it
    // from mpiPM->get_directional_Ncells() I suppose...  N.B. qv->dims is the
    // number of data entries in each direction (by contrast quadmesh has
    // qm->dims[] = num.nodes = qv->dims[]+1).
    //
    qm_NX[v] = qv->dims[v];
  }

  class cell *cx = c, *cy = c, *cz = c;
  long int ct       = 0;
  long int qv_index = 0;
  double norm       = 1.0 / sqrt(4.0 * M_PI);

  while (cz != 0) {
    //
    // Trace an x-y plane in the current z-plane.
    // Set index, and reset the y-index counter
    //
    qm_ix[YY] = qm_start[YY];
    cy        = cz;

    while (cy != 0) {
      //
      // Trace an x-column starting at the mesh start point.
      //
      cx = cy;
      //
      // Get to starting index [NX*NY*iz+NX*iy+ix]
      //
      qv_index = 0;
      if (ndim > 2) qv_index += qm_NX[XX] * qm_NX[YY] * qm_ix[ZZ];
      if (ndim > 1) qv_index += qm_NX[XX] * qm_ix[YY];
      qm_ix[XX] = qm_start[XX];
      qv_index += qm_ix[XX];

      while ((cx != 0) && cx->pos[XX] < iXmax[XX]) {
#ifdef TEST_SILO_IO
        // rep.printVec("cpos",cx->pos,ndim);
        // rep.printVec("P",cx->P,SimPM.nvar);
#endif
        //
        // Different pointers if float or double.
        //
        if (silo_dtype == DB_FLOAT) {
          cx->P[v1] = fdata[0][qv_index];
          if (v2 > 0) cx->P[v2] = fdata[1][qv_index];
          if (v3 > 0) cx->P[v3] = fdata[2][qv_index];
        }
        else {
          cx->P[v1] = ddata[0][qv_index];
          if (v2 > 0) cx->P[v2] = ddata[1][qv_index];
          if (v3 > 0) cx->P[v3] = ddata[2][qv_index];
        }
#ifdef NEW_B_NORM
        if (B) {
          // cout <<"PP_read_var2grid: scaling B var "<<v1;
          // cout <<" val="<< cx->P[v1]<<" by "<<norm<<"\n";
          cx->P[v1] *= norm;
          if (v2 > 0) cx->P[v2] *= norm;
          if (v2 > 0) cx->P[v3] *= norm;
        }
#endif

        cx = ggg->NextPt(*cx, XP);
        qv_index++;
        qm_ix[XX]++;
        ct++;
      }  // x-column

      if (ndim > 1) {
        //
        // move to next x-column in YP direction, if it exists, and if
        // it is on the mesh domain.  Also increment qm_ix[YY] to
        // indicate this.
        //
        cy = ggg->NextPt(*cy, YP);
        if (cy != 0 && cy->pos[YY] > iXmax[YY]) cy = 0;
        qm_ix[YY]++;
      }
      else {
        //
        // ndim==1, so we want to break out of the y-dir loop
        //
        cy = 0;
      }

    }  // y-loop

    if (ndim > 2) {
      //
      // move to next XY-plane in the ZP direction, if it exists and
      // if it is on the mesh domain.  Also increment the qm_ix[ZZ]
      // counter.
      //
      cz = ggg->NextPt(*cz, ZP);
      if (cz != 0 && cz->pos[ZZ] > iXmax[ZZ]) cz = 0;
      qm_ix[ZZ]++;
    }
    else {
      //
      // ndim<=2, so we want to break out of the z-dir loop
      //
      cz = 0;
    }

  }  // z-loop.

#ifdef TEST_SILO_IO
  spdlog::debug("Read variable {}, got {} cells", variable, ct);
#endif
  DBFreeQuadvar(qv);  // qv=0;
  data  = 0;
  fdata = 0;
  ddata = 0;
  return 0;
}

// ##################################################################
// ##################################################################

/********************************************************/
/*************** dataio_silo_utility ********************/
/********************************************************/
