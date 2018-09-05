///
/// \file   dataio_NG_silo.cc
/// \author Jonathan Mackey
/// \date   2018-05-04
/// 
///
/// Modifications:\n
/// - 2018.05.04 JM: new read/write functions for a NG grid.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SILO

#include "dataIO/dataio_silo_NG.h"
#include "microphysics/microphysics_base.h"

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#include "raytracing/raytracer_base.h"
#endif // TESTING
#ifdef RT_TESTING_OUTPUTCOL
#include "raytracing/raytracer_base.h"
#endif // RT_TESTING_OUTPUTCOL

#include <cstring>
#include <sstream>
using namespace std;




// ##################################################################
// ##################################################################


dataio_NG_silo::dataio_NG_silo(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      std::string dtype // read/write either FLOAT or DOUBLE to/from file
      )
: dataio_silo(SimPM,dtype)
{
#ifdef TESTING
  cout <<"setting up dataio_NG_silo class.\n";
#endif
  return;
}




// ##################################################################
// ##################################################################


dataio_NG_silo::~dataio_NG_silo()
{
#ifdef TESTING
  cout <<"deleting dataio_NG_silo class.\n";
#endif
}



// ##################################################################
// ##################################################################



int dataio_NG_silo::OutputData(
      const string outfile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      const long int file_counter   ///< number to stamp file with (e.g. timestep)
      )
{
  int err=0;
  if (!have_setup_writevars) {
    // set what data to write to the mesh.
    err = setup_write_variables(SimPM);
    rep.errorTest("dataio_NG_silo::OutputData() setup_write_variables",0,err);
  }

  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the NG grid.
    ostringstream temp; temp << outfile << "_level";
    temp.width(2); temp.fill('0');
    temp << l;

    if (!cg[l])
      rep.error("dataio_NG_silo::OutputData() null pointer to grid!",cg[l]);
    dataio_silo::gp = cg[l];

    int err=0;

    err = choose_filename(temp.str(), file_counter);
    if (err) {
      cerr<<"dataio_NG_silo::OutputData() error choosing filename.\n";
      return err;
    }
    cout <<"\tWriting to file: "<<silofile<<"\n";
    ofstream fff;
    fff.open(silofile.c_str());
    if (!fff) rep.error("!!!**** Can't write to file.  Does directory exist???",silofile);
    fff.close();

    //
    // Create file
    //
    *db_ptr = 0;
    if (silo_filetype==DB_HDF5) {
      DBSetCompression("METHOD=GZIP LEVEL=1");
      DBSetFriendlyHDF5Names(1);
    }
    *db_ptr = DBCreate(silofile.c_str(), DB_CLOBBER, DB_LOCAL,
                       "Nested-Grid PION data", silo_filetype);
    if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);

    //
    // set grid properties for quadmesh: each level of the NG
    // grid has different zone and node coordinates so we need to
    // call this each time.
    //
    err = setup_grid_properties(gp, SimPM);
    rep.errorTest("dataio_NG_silo::OutputData() setup_grid_properties",0, err);

    //
    // now write the simulation parameters to the header part of the file.
    //
    DBSetDir(*db_ptr,"/");
    DBMkDir(*db_ptr,"header");
    DBSetDir(*db_ptr,"/header");
    err = write_simulation_parameters(SimPM);
    if (err)
      rep.error("dataio_NG_silo::OutputData() error writing header to silo file",err);

    //
    // Create data directory, generate the mesh in the file, and then write each
    // variable in turn to the mesh.
    //
    DBSetDir(*db_ptr,"/");
    string meshname="UniformGrid";
    err = generate_quadmesh(*db_ptr, meshname,SimPM);
    if (err)
      rep.error("dataio_NG_silo::OutputData() error writing quadmesh to silo file",err);

    create_data_arrays(SimPM);
    for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
    err = dataio_silo::write_variable2mesh(SimPM, *db_ptr, meshname, (*i));
    if (err)
      rep.error("dataio_NG_silo::OutputData() error writing variable",(*i));
    }
    delete_data_arrays();
    DBSetDir(*db_ptr,"/");
    DBClose(*db_ptr); //*db_ptr=0;

  } // loop over levels.

  return 0;
}



// ##################################################################
// ##################################################################



int dataio_NG_silo::ReadData(
      string infile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  dataio_silo::gp = cg[0];
  silofile=infile;
  int err=0;

  // Loop over grid levels, and read data for each level.
  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the NG grid.

    string::size_type p;
    if ((p=silofile.find("_level"))==string::npos) {
      rep.error("dataio_NG_silo::ReadData() level",silofile);
    }
    else {
      ostringstream temp; temp.str("");
      temp.width(2); temp.fill('0');
      temp << l;
      silofile.replace(p+6,2,temp.str());
      cout <<"p="<<p<<"  string="<<temp.str()<<", silofile="<<silofile<<"\n";
    }

    if (!cg[l])
      rep.error("dataio_NG_silo::OutputData() null pointer to grid!",cg[l]);
    dataio_silo::gp = cg[l];

    *db_ptr = DBOpen(silofile.c_str(), DB_UNKNOWN, DB_READ);
    if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);

    int ftype = DBGetDriverType(*db_ptr);
    if (ftype==DB_HDF5) {
      int friendly=DBGetFriendlyHDF5Names();
      DBSetFriendlyHDF5Names(friendly);
    }  

    err = set_readvars(SimPM);
    if (err) rep.error("failed to set readvars in ReadData",err);
    //
    // set grid properties for quadmesh: each level of the NG
    // grid has different zone and node coordinates so we need to
    // call this each time.
    //
    err = setup_grid_properties(gp, SimPM);
    rep.errorTest("dataio_NG_silo::OutputData() setup_grid_properties",0, err);

    DBSetDir(*db_ptr,"/");
    string meshname="UniformGrid";

    // now read each variable in turn from the mesh
    for (std::vector<string>::iterator i=readvars.begin(); i!=readvars.end(); ++i) {
#ifdef WRITE_GHOST_ZONES
      err = dataio_silo::read_variable2grid(SimPM, *db_ptr, meshname, (*i), gp->Ncell_all());
#else
      err = dataio_silo::read_variable2grid(SimPM, *db_ptr, meshname, (*i), gp->Ncell());
#endif
      if (err)
        rep.error("dataio_NG_silo::ReadData() error reading variable",(*i));
    }
    DBSetDir(*db_ptr,"/");
    DBClose(*db_ptr); //*db_ptr=0; 
  } // levels

  return err;
}


// ##################################################################
// ##################################################################


#endif // SILO

