/// \file dataio_fits_pllel.cc
/// \author Jonathan Mackey
///
/// This file contains the class definitions for the DataIOFits class, which
/// uses FITS, from
/// http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html at NASA.  I
/// have tried to make the functions atomic as much as possible, so bits of code
/// can be reused.
///
/// One thing to note is that writing to a single file in parallel
/// (from multiple processes) is broken -- there is no useful queueing
/// system to make processes take turns, so it can and probably will
/// seg.fault.  Workaround is to write multiple files and then stitch
/// them together later with grid/analysis/stitchfits.cc
///
/// modified:\n
/// - 2015.01.28 JM: parallel version of the serial functions.

#ifdef FITS

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
#endif  // NDEBUG

#include "dataIO/dataio_fits_MPI.h"
#ifdef RT_TESTING_OUTPUTCOL
#include "raytracing/raytracer_base.h"
#endif  // RT_TESTING_OUTPUTCOL
#include "microphysics/microphysics_base.h"

#include "fitsio.h"
#include <cstring>
#include <sstream>
#include <vector>
using namespace std;

// ##################################################################
// ##################################################################

//-------------------------------------------------
//-------------  FITS DATA I/O CLASS  -------------
//-------------------------------------------------

// ##################################################################
// ##################################################################

DataIOFits_pllel::DataIOFits_pllel(
    class SimParams &SimPM,  ///< pointer to simulation parameters
    class Sub_domain *p) :
    DataIOFits(SimPM)
{
  spdlog::info("Setting up DataIOFits_pllel class");
  DataIOFits_pllel::mpiPM = p;
}

// ##################################################################
// ##################################################################

DataIOFits_pllel::~DataIOFits_pllel()
{
  spdlog::info("Deleting DataIOFits_pllel class");
  DataIOFits_pllel::mpiPM = 0;
}

// ##################################################################
// ##################################################################

std::string DataIOFits_pllel::choose_filename(
    const std::string fbase,  ///< filebase passed in from main code.
    const int file_counter    ///< file counter to use (e.g. timestep).
)
{
  //
  // Choose filename based on the basename and the counter passed to
  // this function.
  //
  string outfile;
  ostringstream temp;
  temp.str("");
  temp << fbase;
  //
  // Add _RANK to filename if running in parallel and writing multiple
  // files.
  //
  if (!mpiPM->get_WriteSingleFile()) {
    temp << "_";
    temp.width(4);
    temp.fill('0');
    temp << mpiPM->get_myrank();
  }
  temp << ".";
  if (file_counter >= 0) {
    temp.width(Ndigits);
    temp.fill('0');
    temp << file_counter << ".";
  }
  temp << "fits";
  outfile = temp.str();
  temp.str("");
  return outfile;
}

// ##################################################################
// ##################################################################

int DataIOFits_pllel::OutputData(
    string outfilebase,                 ///< base filename
    vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
    class SimParams &SimPM,             ///< pointer to simulation parameters
    const long int file_counter  ///< number to stamp file with (e.g. timestep)
)
{
  int err = 0;
  // loop over grid refinement levels, save one file per level
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug("SAVING DATA FOR LEVEL {} IN FITS FORMAT", l);
    if (!cg[l])
      spdlog::error(
          "{}: {}", "dataio_silo::OutputData() null pointer!", fmt::ptr(cg[l]));
    DataIOFits_pllel::gp = cg[l];
    mpiPM                = &(SimPM.levels[l].sub_domain);

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

    err = SaveLevelData(fbase, l, cg[l], SimPM, file_counter);
    if (0 != err)
      spdlog::error("{}: Expected {} but got {}", "saveleveldata", 0, err);
  }
  return err;
}

// ##################################################################
// ##################################################################

int DataIOFits_pllel::SaveLevelData(
    string outfilebase,          ///< base filename
    const int l,                 ///< level to save
    class GridBaseClass *cg,     ///< address of vector of grid pointers.
    class SimParams &SimPM,      ///< pointer to simulation parameters
    const long int file_counter  ///< number to stamp file with (e.g. timestep)
)
{
  string fname = "DataIOFits_pllel::SaveLevelData";
  mpiPM        = &(SimPM.levels[l].sub_domain);

  if (!cg)
    spdlog::error(
        "{}: {}", "DataIOFits_pllel::OutputData() null pointer to grid!",
        fmt::ptr(cg));
  DataIOFits_pllel::gp = cg;

  if (DataIOFits_pllel::eqn == 0) {
    // cout <<"WARNING: DataIOFits_pllel::OutputData() Set up Equations
    // pointer before outputting data!\n"; cout <<"WARNING:
    // DataIOFits_pllel::OutputData() Not outputting Eint/divB/Ptot b/c no
    // way to calculate it!\n";
  }
  fitsfile *ff = 0;
  int status = 0, err = 0;

  // Add variables to list based on what equations we are solving.
  int nvar        = SimPM.nvar;
  string *extname = 0;
  if (SimPM.ntracer > 5)
    spdlog::error(
        "{}: {}",
        "OutputFitsData:: only handles up to 5 tracer variables! Add "
        "more if needed.",
        SimPM.ntracer);

#ifdef RT_TESTING_OUTPUTCOL
  // output column densities!
  if (gp->RT != 0 && SimPM.RS.Nsources > 0) {
    for (int si = 0; si < SimPM.RS.Nsources; si++) {
      nvar += SimPM.RS.sources[si].NTau;  // for column densities
    }
  }
#endif  // RT_TESTING_OUTPUTCOL

  if (SimPM.eqntype == EQEUL || SimPM.eqntype == EQEUL_ISO
      || SimPM.eqntype == EQEUL_EINT) {
    extname         = mem.myalloc(extname, nvar + 1);
    string pvar[10] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                       "TR0",     "TR1",     "TR2",   "TR3",   "TR4"};
    for (int i = 0; i < SimPM.nvar; i++)
      extname[i] = pvar[i];
    if (DataIOFits::eqn != 0 && mp == 0) {
      extname[nvar] = "Eint";
      nvar += 1;
    }
    else if (mp != 0) {
      extname[nvar] = "Temp";
      nvar += 1;
    }
  }
  else if (SimPM.eqntype == EQMHD || SimPM.eqntype == EQFCD) {
    extname         = mem.myalloc(extname, nvar + 3);
    string pvar[13] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                       "Bx",      "By",      "Bz",    "TR0",   "TR1",
                       "TR2",     "TR3",     "TR4"};
    for (int i = 0; i < SimPM.nvar; i++)
      extname[i] = pvar[i];
    if (DataIOFits::eqn != 0 && mp == 0) {
      extname[nvar]     = "Eint";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (DataIOFits::eqn != 0 && mp != 0) {
      extname[nvar]     = "Temp";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (mp != 0) {
      extname[nvar] = "Temp";
      nvar += 1;
    }
  }
  else if (SimPM.eqntype == EQGLM) {
    extname         = mem.myalloc(extname, nvar + 3);
    string pvar[14] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                       "Bx",      "By",      "Bz",    "psi",   "TR0",
                       "TR1",     "TR2",     "TR3",   "TR4"};
    for (int i = 0; i < SimPM.nvar; i++)
      extname[i] = pvar[i];
    spdlog::debug("EQN={}, MP={}", fmt::ptr(DataIOFits::eqn), fmt::ptr(mp));
    if (DataIOFits::eqn != 0 && mp == 0) {
      extname[nvar]     = "Eint";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (DataIOFits::eqn != 0 && mp != 0) {
      extname[nvar]     = "Temp";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (mp != 0) {
      extname[nvar] = "Temp";
      nvar += 1;
    }
  }
  else {
    extname = mem.myalloc(extname, 10);
    spdlog::error("{}: {}", "What equations?!", SimPM.eqntype);
  }

#ifdef RT_TESTING_OUTPUTCOL
  // output column densities!
  if (gp->RT != 0 && SimPM.RS.Nsources > 0) {
    if (extname[SimPM.nvar] != "")
      spdlog::error("{}: {}", "Tau not writeable!", extname[SimPM.nvar]);
    //
    // Loop over all sources, and all variables for each source:
    //
    ostringstream var;
    unsigned int ivar = SimPM.nvar;
    for (int v = 0; v < SimPM.RS.Nsources; v++) {
      for (int iT = 0; iT < SimPM.RS.sources[v].NTau; iT++) {
        var.str("");
        var << "Col_Src_" << v << "_T" << iT;
        extname[ivar] = var.str();
        ivar++;
      }  // loop over Tau variables for source v.
    }    // loop over Nsources
  }      // if RT
#endif   // RT_TESTING_OUTPUTCOL

  //
  // Choose filename based on the basename and the counter passed to
  // this function.
  //
  string outfile = choose_filename(outfilebase, file_counter);
  ostringstream temp;
  temp.str("");

  // -------------------------------------------------------
  // -------------------------------------------------------
  if (!mpiPM->get_WriteSingleFile()) {
    // This is the default -- each process writes its own file
    spdlog::debug("DataIOFits_pllel::OutputData() writing multiple files");
    temp.str("");
    spdlog::debug(
        "Proc {}:\t writing to file {}", mpiPM->get_myrank(), outfile);
    if (file_exists(outfile)) {
      spdlog::debug(
          "Proc {}:\t file exists... overwriting!", mpiPM->get_myrank());
      temp.str("");
      temp << "!" << outfile;
      outfile = temp.str();
    }
    //    if(acquire_lock(outfile)) spdlog::error("{}: {}", "Failed to lock
    //    file",err);

    // Create fits file.
    fits_create_file(&ff, outfile.c_str(), &status);
    if (status) {
      spdlog::error("Creating new file went bad");
      exit(1);
    }

    // write fits header
    //    err += write_fits_header(ff);
    //    if(err) spdlog::error("{}: {}", "DataIOFits_pllel::OutputData()
    //    couldn't write fits header",err);
    // --------------------------------------------------------
    //
    // create HDU for header
    //
    fits_create_img(ff, DOUBLE_IMG, 0, 0, &status);
    if (status) {
      fits_report_error(stderr, status);
    }
    //
    // set file pointer for the header writing function.
    //
    file_ptr = ff;
    err      = write_simulation_parameters(SimPM);
    if (err)
      spdlog::error(
          "{}: {}", "DataIOFits_pllel::OutputData() couldn't write fits header",
          err);
    ff = file_ptr;

    // err = fits_open_file(&ff, outfile.c_str(), READWRITE, &status);
    // if(status) {fits_report_error(stderr,status); return(err);}
    int num = -1;
    fits_get_hdu_num(ff, &num);
    if (num != 1) ffmahd(ff, 1, 0, &status);
    if (status) {
      fits_report_error(stderr, status);
      spdlog::error("{}: {}", "NG can't find fits header", status);
    }
    char key[128];
    int lev = l;
    strcpy(key, "grid_level");
    err += fits_update_key(ff, TINT, key, &lev, 0, &status);
    strcpy(key, "level_xmin0");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmin[0]), 0, &status);
    strcpy(key, "level_xmin1");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmin[1]), 0, &status);
    strcpy(key, "level_xmin2");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmin[2]), 0, &status);
    strcpy(key, "level_xmax0");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmax[0]), 0, &status);
    strcpy(key, "level_xmax1");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmax[1]), 0, &status);
    strcpy(key, "level_xmax2");
    err += fits_update_key(
        file_ptr, TDOUBLE, key, &(SimPM.levels[l].Xmax[2]), 0, &status);

    // --------------------------------------------------------

    //
    // for each image, create image and write my portion of it.
    //
    double *data = 0;
    for (int i = 0; i < nvar; i++) {
      if (mpiPM->get_WriteFullImage()) {  // write full image, with only local
                                          // part being non-zero.
        err += create_fits_image(ff, extname[i], SimPM.ndim, SimPM.NG.data());
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->get_Ncell(), &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->get_Xmin().data(),
            SimPM.levels[l].Xmin.data(), SimPM.levels[l].dx, SimPM.ndim,
            mpiPM->get_directional_Ncells().data(), mpiPM->get_Ncell(), data);
      }
      else {  // Write only part of image that is on local grid.
        err += create_fits_image(
            ff, extname[i], SimPM.ndim, mpiPM->get_directional_Ncells().data());
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->get_Ncell(), &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->get_Xmin().data(), mpiPM->get_Xmin().data(),
            SimPM.levels[l].dx, SimPM.ndim,
            mpiPM->get_directional_Ncells().data(), mpiPM->get_Ncell(), data);
      }
    }
    if (err)
      spdlog::error(
          "{}: {}", "DataIOFits_pllel::OutputData() Image Writing went bad",
          err);
    data = mem.myfree(data);

    // Close file
    err += fits_close_file(ff, &status);
    //    release_lock(outfile);
    spdlog::debug(
        "Proc {}: file created, written, and unlocked. err={}",
        mpiPM->get_myrank(), err);
  }

  else if (mpiPM->get_WriteSingleFile()) {
    spdlog::warn(
        "THIS IS NOT SYNCHRONOUS!  WILL FAIL RANDOMLY AND CAUSE CRASH");
    if (file_exists(outfile)) {
      spdlog::debug("Proc {}: file exists...", mpiPM->get_myrank());
      // If file exists, wait for access and then lock it.
      acquire_lock(outfile);
    }
    else {
      // If file doesn't exist, lock it and create it.
      spdlog::debug(
          "Proc {}: file doesn't exist, lock it and create it",
          mpiPM->get_myrank());
      acquire_lock(outfile);
      // Create fits file.
      fits_create_file(&ff, outfile.c_str(), &status);
      if (status) {
        fits_report_error(stderr, status);
        spdlog::error("Creating new file went bad");
        exit(1);
      }

      // write fits header
      //      err += write_fits_header(ff);
      //      if(err) spdlog::error("{}: {}", "DataIOFits_pllel::OutputData()
      //      couldn't write fits header",err);
      // --------------------------------------------------------
      //
      // create HDU for header
      //
      fits_create_img(ff, DOUBLE_IMG, 0, 0, &status);
      if (status) {
        fits_report_error(stderr, status);
      }
      //
      // set file pointer for the header writing function.
      //
      file_ptr = ff;
      err      = write_simulation_parameters(SimPM);
      if (err)
        spdlog::error(
            "{}: {}",
            "DataIOFits_pllel::OutputData() couldn't write fits header", err);
      ff = file_ptr;
      // --------------------------------------------------------

      //
      // for each image, create image and write my portion of it.
      //
      for (int i = 0; i < nvar; i++) {
        err += create_fits_image(ff, extname[i], SimPM.ndim, SimPM.NG.data());
        double *data = 0;
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->get_Ncell(), &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->get_Xmin().data(),
            SimPM.levels[l].Xmin.data(), SimPM.levels[l].dx, SimPM.ndim,
            mpiPM->get_directional_Ncells().data(), mpiPM->get_Ncell(), data);
        data = mem.myfree(data);
      }
      if (err)
        spdlog::error(
            "{}: {}",
            "DataIOFits_pllel::OutputData() SingleFile Image Writing went bad",
            err);
      // Close file
      err += fits_close_file(ff, &status);
      release_lock(outfile);
      spdlog::debug(
          "Proc {}: file created, written, and unlocked. err={}",
          mpiPM->get_myrank(), err);
      delete[] extname;
      extname = 0;
      return (err);
    }
    //--------------------------------------------------------
    // Now file exists, and we have access to it.
    err = fits_open_file(&ff, outfile.c_str(), READWRITE, &status);
    if (status) {
      fits_report_error(stderr, status);
      return (err);
    }
    // Move to the right hdu image and check it is the right size.
    char temp[256];
    for (int i = 0; i < nvar; i++) {
      strcpy(temp, extname[i].c_str());
      err += fits_movnam_hdu(ff, ANY_HDU, temp, 0, &status);
      if (status) {
        fits_report_error(stderr, status);
        spdlog::error("Couldn't find hdu {}", temp);
        return (err);
      }
      err = check_fits_image_dimensions(
          ff, extname[i], SimPM.ndim, SimPM.NG.data());
      if (err != 0)
        spdlog::error(
            "{}: {}",
            "DataIOFits_pllel::OutputData() SingleFile: image dimensions "
            "don't match",
            err);
      // write my portion of image.
      double *data = 0;
      err += put_variable_into_data_array(
          SimPM, extname[i], mpiPM->get_Ncell(), &data);
      err += write_fits_image(
          ff, extname[i], mpiPM->get_Xmin().data(), SimPM.levels[l].Xmin.data(),
          SimPM.levels[l].dx, SimPM.ndim,
          mpiPM->get_directional_Ncells().data(), mpiPM->get_Ncell(), data);
      data = mem.myfree(data);
    }
    if (err)
      spdlog::error(
          "{}: {}",
          "DataIOFits_pllel::OutputData() SingleFile: Error writing image",
          err);
    // Close file
    err += fits_close_file(ff, &status);
    if (status) {
      fits_report_error(stderr, status);
      spdlog::error(
          "DataIOFits_pllel::OutputData() SingleFile: Error writing to existing file");
      return (err);
    }
    release_lock(outfile);
  }  // If write to a single file.

  // -------------------------------------------------------
  // -------------------------------------------------------
  else
    spdlog::error(
        "{}: {}", "Logic Error in DataIOFits_pllel::OutputData()",
        mpiPM->get_WriteSingleFile());

  extname = mem.myfree(extname);

  if (status) {
    fits_report_error(stderr, status);
    fits_clear_errmsg();
    return (status);
  }
  return err;
}

// ##################################################################
// ##################################################################

int DataIOFits_pllel::ReadData(
    string infile,
    vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
    class SimParams &SimPM              ///< pointer to simulation parameters
)
{
  string fname = "DataIOFits_pllel::ReadData";
  int err      = 0;

  // -------------------------------------------------------
  int nvar    = SimPM.nvar;
  string *var = 0;
  if (SimPM.ntracer > 5) {
    spdlog::error(
        "{}: {}", "Fits_pllel::ReadData() too many tracer variables!",
        SimPM.ntracer);
  }
  if (SimPM.eqntype == EQEUL || SimPM.eqntype == EQEUL_ISO
      || SimPM.eqntype == EQEUL_EINT) {
    var = mem.myalloc(var, 10);
    if (SimPM.nvar > 10) {
      spdlog::error(
          "{}: {}", "DataIOFits_pllel::ReadData() too many tracers.",
          SimPM.nvar);
    }
    string t[10] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                    "TR0",     "TR1",     "TR2",   "TR3",   "TR4"};
    for (int i = 0; i < 10; i++)
      var[i] = t[i];
  }
  else if (
      SimPM.eqntype == EQMHD || SimPM.eqntype == EQGLM
      || SimPM.eqntype == EQFCD) {
    var = mem.myalloc(var, 14);
    if (SimPM.nvar > 14) {
      spdlog::error(
          "{}: {}", "DataIOFits_pllel::ReadData() too many tracers.",
          SimPM.nvar);
    }
    string t[14] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                    "Bx",      "By",      "Bz",    "psi",   "TR0",
                    "TR1",     "TR2",     "TR3",   "TR4"};
    for (int i = 0; i < 14; i++)
      var[i] = t[i];
  }
  else {
    var = mem.myalloc(var, 10);
    spdlog::error("{}: {}", "What equations?!", SimPM.eqntype);
  }
  // -------------------------------------------

  // -------------------------------------------
  // -------------------------------------------
  // loop over grid refinement levels, read one file per level per proc
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    spdlog::debug("READING DATA FOR LEVEL {} IN FITS FORMAT", l);
    if (!cg[l])
      spdlog::error(
          "{}: {}", "dataio_silo::ReadData() null pointer!", fmt::ptr(cg[l]));
    DataIOFits_pllel::gp = cg[l];
    mpiPM                = &(SimPM.levels[l].sub_domain);

    // read a different file for each level in the NG grid.
    // If more than one level of grid, look for level in filename:
    string::size_type p;
    if ((p = infile.find("_level")) == string::npos && SimPM.grid_nlevels > 1) {
      spdlog::error("{}: {}", "::ReadData() level", infile);
    }
    else if (SimPM.grid_nlevels > 1) {
      ostringstream temp;
      temp.str("");
      temp.width(2);
      temp.fill('0');
      temp << l;
      infile.replace(p + 6, 2, temp.str());
      spdlog::debug("p={}  string={}, silofile={}", p, temp.str(), infile);
    }

    int status = 0;
    fitsfile *ff;
    spdlog::debug(
        "DataIOFits_pllel::ReadData() opening fits file to read data...");
    err = fits_open_file(&ff, infile.c_str(), READONLY, &status);
    if (status) {
      fits_report_error(stderr, status);
      return (err);
    }
    // cout <<"done.\n";
    // Move to first data hdu; should be 2nd hdu;
    int num;
    fits_get_hdu_num(ff, &num);
    if (num != 1) err += ffmahd(ff, 1, 0, &status);
    if (status) {
      fits_report_error(stderr, status);
      return (err);
    }
    err += ffmrhd(ff, 1, 0, &status);
    fits_get_hdu_num(ff, &num);
    spdlog::debug("Current hdu: {}\t and err={}", num, err);

    // -------------------------------------------
    // Loop over all Variables and read from file.
    // -------------------------------------------
    char temp[32];
    for (int i = 0; i < nvar; i++) {
      strcpy(temp, var[i].c_str());
      int v = 0;
      if (var[i] == "GasDens")
        v = static_cast<int>(RO);
      else if (var[i] == "GasPres")
        v = static_cast<int>(PG);
      else if (var[i] == "GasVX")
        v = static_cast<int>(VX);
      else if (var[i] == "GasVY")
        v = static_cast<int>(VY);
      else if (var[i] == "GasVZ")
        v = static_cast<int>(VZ);
      else if (var[i] == "Bx")
        v = static_cast<int>(BX);
      else if (var[i] == "By")
        v = static_cast<int>(BY);
      else if (var[i] == "Bz")
        v = static_cast<int>(BZ);
      else if (var[i] == "psi")
        v = static_cast<int>(SI);
      else if (var[i] == "TR0") {
        v = SimPM.ftr;
        spdlog::debug("reading from first tracer var: {}", v);
      }
      else if (var[i] == "TR1")
        v = SimPM.ftr + 1;
      else if (var[i] == "TR2")
        v = SimPM.ftr + 2;
      else if (var[i] == "TR3")
        v = SimPM.ftr + 3;
      else if (var[i] == "TR4")
        v = SimPM.ftr + 4;
      else
        spdlog::error(
            "{}: {}", "Bad variable index in fits read routine", var[i]);
      err += fits_movnam_hdu(ff, ANY_HDU, temp, 0, &status);

      if (err != 0) {
        // If can't find variable, set them all to zero.
        cell *c = gp->FirstPt();
        do {
          c->P[v] = 0.;
        } while ((c = gp->NextPt(c)) != 0);
        if (status) {
          fits_report_error(stderr, status);
        }
        err = 0;
        fits_clear_errmsg();
        status = 0;
        spdlog::debug(
            "couldn't get data for variable {}; will set data to zero and hope for the best",
            temp);
      }

      else {
        // Variable found, check we're at the right hdu and read data.
        fits_get_hdu_num(ff, &num);
        spdlog::debug("Current hdu: {}\t i={} and var[i] = {}", num, i, var[i]);
        // -----------------------------------------------------------------
        // --- Now call read function differently depending on if infile
        // ---
        // --- is a single file or split already between processors. ---
        // -----------------------------------------------------------------
        if (!mpiPM->get_ReadSingleFile()) {
          // This is where each process reads from its own file.
          spdlog::debug(
              "DataIOFits_pllel::ReadData() Reading from multiple files.\nProc {}:\t reading from file {}\n, directional_Ncells=",
              mpiPM->get_myrank(), infile);
          spdlog::debug("  : {}", mpiPM->get_directional_Ncells());
          err += check_fits_image_dimensions(
              ff, var[i], SimPM.ndim, mpiPM->get_directional_Ncells().data());
          if (err) spdlog::error("{}: {}", "image wrong size.", err);
          err += read_fits_image(
              SimPM, ff, var[i], mpiPM->get_Xmin().data(),
              mpiPM->get_Xmin().data(), mpiPM->get_directional_Ncells().data(),
              mpiPM->get_Ncell());
          if (err) spdlog::error("{}: {}", "error reading image.", err);
        }
        else if (mpiPM->get_ReadSingleFile()) {
          // All processes read from a single ic/restart file.
          spdlog::debug(
              "DataIOFits_pllel::ReadData() Reading from single file.\nProc {}:\t reading from file {}\n",
              mpiPM->get_myrank(), infile);
          err += check_fits_image_dimensions(
              ff, var[i], SimPM.ndim, SimPM.NG.data());
          if (err) spdlog::error("{}: {}", "image wrong size.", err);
          err += read_fits_image(
              SimPM, ff, var[i], mpiPM->get_Xmin().data(), SimPM.Xmin.data(),
              mpiPM->get_directional_Ncells().data(), mpiPM->get_Ncell());
          if (err) spdlog::error("{}: {}", "error reading image.", err);
        }
        else
          spdlog::error(
              "{}: {}", "DataIOFits_pllel::ReadData() logic error",
              mpiPM->get_WriteSingleFile());
        //------------------------------------------------------------------
      }  // got real hdu and read data.

    }  // Loop over all Primitive Variables

    var = mem.myfree(var);
    //  cout <<"Closing fits file. err="<<err<<"\n";
    fits_close_file(ff, &status);
    if (status) {
      fits_report_error(stderr, status);
      fits_clear_errmsg();
    }
    //  cout <<"Closed fits file. err="<<err<<"\n";

    // Now assign Ph to be equal to P for each cell.
    cell *cpt = gp->FirstPt();
    do {
      for (int v = 0; v < nvar; v++)
        cpt->Ph[v] = cpt->P[v];
    } while ((cpt = gp->NextPt(cpt)) != 0);
  }
  return err;
}

// ##################################################################
// ##################################################################

#endif  // if FITS
