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
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

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
    class MCMDcontrol *p) :
    DataIOFits(SimPM)
{
  cout << "Setting up DataIOFits_pllel class.\n";
  DataIOFits_pllel::mpiPM = p;
}

// ##################################################################
// ##################################################################

DataIOFits_pllel::~DataIOFits_pllel()
{
  cout << "Deleting DataIOFits_pllel class.\n";
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
  if (!mpiPM->WriteSingleFile) {
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
    cout << "SAVING DATA FOR LEVEL " << l << " IN FITS FORMAT\n";
    if (!cg[l]) rep.error("dataio_silo::OutputData() null pointer!", cg[l]);
    DataIOFits_pllel::gp = cg[l];
    mpiPM                = &(SimPM.levels[l].MCMD);

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
    rep.errorTest("saveleveldata", 0, err);
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
  mpiPM        = &(SimPM.levels[l].MCMD);

  if (!cg)
    rep.error("DataIOFits_pllel::OutputData() null pointer to grid!", cg);
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
    rep.error(
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
    if (DataIOFits::eqn != 0 && MP == 0) {
      extname[nvar] = "Eint";
      nvar += 1;
    }
    else if (MP != 0) {
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
    if (DataIOFits::eqn != 0 && MP == 0) {
      extname[nvar]     = "Eint";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (DataIOFits::eqn != 0 && MP != 0) {
      extname[nvar]     = "Temp";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (MP != 0) {
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
    cout << "EQN=" << DataIOFits::eqn << ", MP=" << MP << "\n";
    if (DataIOFits::eqn != 0 && MP == 0) {
      extname[nvar]     = "Eint";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (DataIOFits::eqn != 0 && MP != 0) {
      extname[nvar]     = "Temp";
      extname[nvar + 1] = "divB";
      extname[nvar + 2] = "Ptot";
      nvar += 3;
    }
    else if (MP != 0) {
      extname[nvar] = "Temp";
      nvar += 1;
    }
  }
  else {
    extname = mem.myalloc(extname, 10);
    rep.error("What equations?!", SimPM.eqntype);
  }

#ifdef RT_TESTING_OUTPUTCOL
  // output column densities!
  if (gp->RT != 0 && SimPM.RS.Nsources > 0) {
    if (extname[SimPM.nvar] != "")
      rep.error("Tau not writeable!", extname[SimPM.nvar]);
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
  if (!mpiPM->WriteSingleFile) {
    // This is the default -- each process writes its own file
    cout << "DataIOFits_pllel::OutputData() writing multiple files.\n";
    temp.str("");
    cout << "Proc " << mpiPM->get_myrank() << ":\t writing to file " << outfile
         << "\n";
    if (file_exists(outfile)) {
      cout << "Proc " << mpiPM->get_myrank()
           << ":\t file exists... overwriting!\n";
      temp.str("");
      temp << "!" << outfile;
      outfile = temp.str();
    }
    //    if(acquire_lock(outfile)) rep.error("Failed to lock file",err);

    // Create fits file.
    fits_create_file(&ff, outfile.c_str(), &status);
    if (status) {
      cerr << "Creating new file went bad.\n";
      exit(1);
    }

    // write fits header
    //    err += write_fits_header(ff);
    //    if(err) rep.error("DataIOFits_pllel::OutputData() couldn't write
    //    fits header",err);
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
      rep.error(
          "DataIOFits_pllel::OutputData() couldn't write fits header", err);
    ff = file_ptr;

    // err = fits_open_file(&ff, outfile.c_str(), READWRITE, &status);
    // if(status) {fits_report_error(stderr,status); return(err);}
    int num = -1;
    fits_get_hdu_num(ff, &num);
    if (num != 1) ffmahd(ff, 1, 0, &status);
    if (status) {
      fits_report_error(stderr, status);
      rep.error("NG can't find fits header", status);
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
      if (mpiPM->WriteFullImage) {  // write full image, with only local
                                    // part being non-zero.
        err += create_fits_image(ff, extname[i], SimPM.ndim, SimPM.NG);
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->LocalNcell, &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->LocalXmin, SimPM.levels[l].Xmin,
            SimPM.levels[l].dx, SimPM.ndim, mpiPM->LocalNG, mpiPM->LocalNcell,
            data);
      }
      else {  // Write only part of image that is on local grid.
        err += create_fits_image(ff, extname[i], SimPM.ndim, mpiPM->LocalNG);
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->LocalNcell, &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->LocalXmin, mpiPM->LocalXmin,
            SimPM.levels[l].dx, SimPM.ndim, mpiPM->LocalNG, mpiPM->LocalNcell,
            data);
      }
    }
    if (err)
      rep.error("DataIOFits_pllel::OutputData() Image Writing went bad", err);
    data = mem.myfree(data);

    // Close file
    err += fits_close_file(ff, &status);
    //    release_lock(outfile);
    cout << "Proc " << mpiPM->get_myrank()
         << ": file created, written, and unlocked. err=" << err << "\n";
  }

  else if (mpiPM->WriteSingleFile) {
    cout << "WARNING! THIS IS NOT SYNCHRONOUS!  WILL FAIL RANDOMLY AND CAUSE "
            "CRASH.\n";
    if (file_exists(outfile)) {
      cout << "Proc " << mpiPM->get_myrank() << ": file exists...";
      // If file exists, wait for access and then lock it.
      acquire_lock(outfile);
    }
    else {
      // If file doesn't exist, lock it and create it.
      cout << "Proc " << mpiPM->get_myrank()
           << ": file doesn't exist, lock it and create it.\n";
      acquire_lock(outfile);
      // Create fits file.
      fits_create_file(&ff, outfile.c_str(), &status);
      if (status) {
        fits_report_error(stderr, status);
        cerr << "Creating new file went bad.\n";
        exit(1);
      }

      // write fits header
      //      err += write_fits_header(ff);
      //      if(err) rep.error("DataIOFits_pllel::OutputData() couldn't
      //      write fits header",err);
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
        rep.error(
            "DataIOFits_pllel::OutputData() couldn't write fits header", err);
      ff = file_ptr;
      // --------------------------------------------------------

      //
      // for each image, create image and write my portion of it.
      //
      for (int i = 0; i < nvar; i++) {
        err += create_fits_image(ff, extname[i], SimPM.ndim, SimPM.NG);
        double *data = 0;
        err += put_variable_into_data_array(
            SimPM, extname[i], mpiPM->LocalNcell, &data);
        err += write_fits_image(
            ff, extname[i], mpiPM->LocalXmin, SimPM.levels[l].Xmin,
            SimPM.levels[l].dx, SimPM.ndim, mpiPM->LocalNG, mpiPM->LocalNcell,
            data);
        data = mem.myfree(data);
      }
      if (err)
        rep.error(
            "DataIOFits_pllel::OutputData() SingleFile Image Writing went bad",
            err);
      // Close file
      err += fits_close_file(ff, &status);
      release_lock(outfile);
      cout << "Proc " << mpiPM->get_myrank()
           << ": file created, written, and unlocked. err=" << err << "\n";
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
        cerr << "Couldn't find hdu " << temp << ".\n";
        return (err);
      }
      err = check_fits_image_dimensions(ff, extname[i], SimPM.ndim, SimPM.NG);
      if (err != 0)
        rep.error(
            "DataIOFits_pllel::OutputData() SingleFile: image dimensions "
            "don't match",
            err);
      // write my portion of image.
      double *data = 0;
      err += put_variable_into_data_array(
          SimPM, extname[i], mpiPM->LocalNcell, &data);
      err += write_fits_image(
          ff, extname[i], mpiPM->LocalXmin, SimPM.levels[l].Xmin,
          SimPM.levels[l].dx, SimPM.ndim, mpiPM->LocalNG, mpiPM->LocalNcell,
          data);
      data = mem.myfree(data);
    }
    if (err)
      rep.error(
          "DataIOFits_pllel::OutputData() SingleFile: Error writing image",
          err);
    // Close file
    err += fits_close_file(ff, &status);
    if (status) {
      fits_report_error(stderr, status);
      cerr << "DataIOFits_pllel::OutputData() SingleFile: Error writing to "
              "existing file\n";
      return (err);
    }
    release_lock(outfile);
  }  // If write to a single file.

  // -------------------------------------------------------
  // -------------------------------------------------------
  else
    rep.error(
        "Logic Error in DataIOFits_pllel::OutputData()",
        mpiPM->WriteSingleFile);

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
    rep.error(
        "Fits_pllel::ReadData() too many tracer variables!", SimPM.ntracer);
  }
  if (SimPM.eqntype == EQEUL || SimPM.eqntype == EQEUL_ISO
      || SimPM.eqntype == EQEUL_EINT) {
    var = mem.myalloc(var, 10);
    if (SimPM.nvar > 10) {
      rep.error("DataIOFits_pllel::ReadData() too many tracers.", SimPM.nvar);
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
      rep.error("DataIOFits_pllel::ReadData() too many tracers.", SimPM.nvar);
    }
    string t[14] = {"GasDens", "GasPres", "GasVX", "GasVY", "GasVZ",
                    "Bx",      "By",      "Bz",    "psi",   "TR0",
                    "TR1",     "TR2",     "TR3",   "TR4"};
    for (int i = 0; i < 14; i++)
      var[i] = t[i];
  }
  else {
    var = mem.myalloc(var, 10);
    rep.error("What equations?!", SimPM.eqntype);
  }
  // -------------------------------------------

  // -------------------------------------------
  // -------------------------------------------
  // loop over grid refinement levels, read one file per level per proc
  for (int l = 0; l < SimPM.grid_nlevels; l++) {
    cout << "READING DATA FOR LEVEL " << l << " IN FITS FORMAT\n";
    if (!cg[l]) rep.error("dataio_silo::ReadData() null pointer!", cg[l]);
    DataIOFits_pllel::gp = cg[l];
    mpiPM                = &(SimPM.levels[l].MCMD);

    // read a different file for each level in the NG grid.
    // If more than one level of grid, look for level in filename:
    string::size_type p;
    if ((p = infile.find("_level")) == string::npos && SimPM.grid_nlevels > 1) {
      rep.error("::ReadData() level", infile);
    }
    else if (SimPM.grid_nlevels > 1) {
      ostringstream temp;
      temp.str("");
      temp.width(2);
      temp.fill('0');
      temp << l;
      infile.replace(p + 6, 2, temp.str());
      cout << "p=" << p << "  string=" << temp.str() << ", silofile=";
      cout << infile << "\n";
    }

    int status = 0;
    fitsfile *ff;
    cout << "DataIOFits_pllel::ReadData() opening fits file to read data...";
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
    cout << "Current hdu: " << num << "\t and err=" << err << "\n";

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
        cout << "reading from first tracer var: " << v << "\n";
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
        rep.error("Bad variable index in fits read routine", var[i]);
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
        cout << "couldn't get data for variable " << temp
             << "; will set data to zero and hope for the best.\n";
      }

      else {
        // Variable found, check we're at the right hdu and read data.
        fits_get_hdu_num(ff, &num);
        cout << "Current hdu: " << num << "\t i=" << i
             << " and var[i] = " << var[i] << "\n";
        // -----------------------------------------------------------------
        // --- Now call read function differently depending on if infile
        // ---
        // --- is a single file or split already between processors. ---
        // -----------------------------------------------------------------
        if (!mpiPM->ReadSingleFile) {
          // This is where each process reads from its own file.
          cout << "DataIOFits_pllel::ReadData() Reading from multiple files.\n";
          cout << "Proc " << mpiPM->get_myrank() << ":\t reading from file "
               << infile << "\n";
          cout << ", localNG=";
          rep.printVec(" ", mpiPM->LocalNG, 2);
          err += check_fits_image_dimensions(
              ff, var[i], SimPM.ndim, mpiPM->LocalNG);
          if (err) rep.error("image wrong size.", err);
          err += read_fits_image(
              SimPM, ff, var[i], mpiPM->LocalXmin, mpiPM->LocalXmin,
              mpiPM->LocalNG, mpiPM->LocalNcell);
          if (err) rep.error("error reading image.", err);
        }
        else if (mpiPM->ReadSingleFile) {
          // All processes read from a single ic/restart file.
          cout << "DataIOFits_pllel::ReadData() Reading from single file.\n";
          cout << "Proc " << mpiPM->get_myrank() << ":\t reading from file "
               << infile << "\n";
          err += check_fits_image_dimensions(ff, var[i], SimPM.ndim, SimPM.NG);
          if (err) rep.error("image wrong size.", err);
          err += read_fits_image(
              SimPM, ff, var[i], mpiPM->LocalXmin, SimPM.Xmin, mpiPM->LocalNG,
              mpiPM->LocalNcell);
          if (err) rep.error("error reading image.", err);
        }
        else
          rep.error(
              "DataIOFits_pllel::ReadData() logic error",
              mpiPM->WriteSingleFile);
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
