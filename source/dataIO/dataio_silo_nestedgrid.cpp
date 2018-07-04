///
/// \file   dataio_nested_silo.cc
/// \author Jonathan Mackey
/// \date   2018-05-04
/// 
///
/// Modifications:\n
/// - 2018.05.04 JM: new read/write functions for a nested grid.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SILO

#include "dataIO/dataio_silo_nestedgrid.h"
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


dataio_nested_silo::dataio_nested_silo(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      std::string dtype // read/write either FLOAT or DOUBLE to/from file
      )
: dataio_silo(SimPM,dtype)
{
#ifdef TESTING
  cout <<"setting up dataio_nested_silo class.\n";
#endif
  return;
}




// ##################################################################
// ##################################################################


dataio_nested_silo::~dataio_nested_silo()
{
#ifdef TESTING
  cout <<"deleting dataio_nested_silo class.\n";
#endif
}



// ##################################################################
// ##################################################################



int dataio_nested_silo::OutputData(
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
    rep.errorTest("dataio_nested_silo::OutputData() setup_write_variables",0,err);
  }

  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the nested grid.
    ostringstream temp; temp << outfile << "_level";
    temp.width(2); temp.fill('0');
    temp << l;

    if (!cg[l])
      rep.error("dataio_nested_silo::OutputData() null pointer to grid!",cg[l]);
    dataio_silo::gp = cg[l];

    int err=0;

    err = choose_filename(temp.str(), file_counter);
    if (err) {
      cerr<<"dataio_nested_silo::OutputData() error choosing filename.\n";
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
    *db_ptr = DBCreate(silofile.c_str(), DB_CLOBBER, DB_LOCAL, "Nested-Grid PION data", silo_filetype);
    if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);
    //cout <<"\tdb_ptr="<<db_ptr<<"\n";
    //cout <<"\t*db_ptr="<<*db_ptr<<"\n";

    //
    // set grid properties for quadmesh: each level of the nested
    // grid has different zone and node coordinates so we need to
    // call this each time.
    //
    err = setup_grid_properties(gp, SimPM);
    rep.errorTest("dataio_nested_silo::OutputData() setup_grid_properties",0, err);

    //
    // now write the simulation parameters to the header part of the file.
    //
    DBSetDir(*db_ptr,"/");
    DBMkDir(*db_ptr,"header");
    DBSetDir(*db_ptr,"/header");
    err = write_simulation_parameters(SimPM);
    if (err)
      rep.error("dataio_nested_silo::OutputData() error writing header to silo file",err);

    //
    // Create data directory, generate the mesh in the file, and then write each
    // variable in turn to the mesh.
    //
    DBSetDir(*db_ptr,"/");
    string meshname="UniformGrid";
    err = generate_quadmesh(*db_ptr, meshname,SimPM);
    if (err)
      rep.error("dataio_nested_silo::OutputData() error writing quadmesh to silo file",err);

    create_data_arrays(SimPM);
    for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
    err = dataio_silo::write_variable2mesh(SimPM, *db_ptr, meshname, (*i));
    if (err)
      rep.error("dataio_nested_silo::OutputData() error writing variable",(*i));
    }
    delete_data_arrays();
    DBSetDir(*db_ptr,"/");
    DBClose(*db_ptr); //*db_ptr=0;

  } // loop over levels.

  return 0;
}



// ##################################################################
// ##################################################################



int dataio_nested_silo::ReadData(
      string infile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (!cg[0])
    rep.error("dataio_nested_silo::ReadData() null pointer to grid!",cg[0]);
  dataio_silo::gp = cg[0];
  silofile=infile;
  int err=0;

  // Loop over grid levels, and read data for each level.
  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the nested grid.

    string::size_type p;
    if ((p=silofile.find("_level"))==string::npos) {
      rep.error("dataio_nested_silo::ReadData() level",silofile);
    }
    else {
      ostringstream temp; temp.str("");
      temp.width(2); temp.fill('0');
      temp << l;
      silofile.replace(p+6,2,temp.str());
      cout <<"p="<<p<<"  string="<<temp.str()<<", silofile="<<silofile<<"\n";
    }

    if (!cg[l])
      rep.error("dataio_nested_silo::OutputData() null pointer to grid!",cg[l]);
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
    // set grid properties for quadmesh: each level of the nested
    // grid has different zone and node coordinates so we need to
    // call this each time.
    //
    err = setup_grid_properties(gp, SimPM);
    rep.errorTest("dataio_nested_silo::OutputData() setup_grid_properties",0, err);

    DBSetDir(*db_ptr,"/");
    string meshname="UniformGrid";

    // now read each variable in turn from the mesh
    for (std::vector<string>::iterator i=readvars.begin(); i!=readvars.end(); ++i) {
      err = dataio_silo::read_variable2grid(SimPM, *db_ptr, meshname, (*i), SimPM.Ncell);
      if (err)
        rep.error("dataio_nested_silo::ReadData() error reading variable",(*i));
    }
    DBSetDir(*db_ptr,"/");
    DBClose(*db_ptr); //*db_ptr=0; 
  } // levels

  return err;
}


// ##################################################################
// ##################################################################



int dataio_nested_silo::setup_grid_properties(
      class GridBaseClass *grid, 
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  // set grid parameters -- UNIFORM FIXED GRID
  if (!grid)
    rep.error("dataio_nested_silo::setup_grid_properties() null grid pointer!",grid);
  double dx=grid->DX();
  //if (node_coords || nodedims || zonedims ||
  //    nodex || nodey || nodez) {
  //  cerr<<"Have already setup variables for grid props! ";
  //  cerr<<"You shouldn't have called me! crazy fool...\n";
  //  return 1000000;
  //}

  dataio_silo::ndim = SimPM.ndim;
  dataio_silo::vec_length = SimPM.eqnNDim;

  if (!nodedims) nodedims = mem.myalloc(nodedims,ndim);
  if (!zonedims) zonedims = mem.myalloc(zonedims,ndim);

  //
  // node_coords is a void pointer, so if we are writing silo data in
  // single or double precision then we need different allocation
  // calls.  Same for nodex, nodey, nodez.
  //
  // We setup arrays with locations of nodes in coordinate directions.
  //
  int nx = SimPM.NG[XX]+1; // for N cells, have N+1 nodes.
  int ny = SimPM.NG[YY]+1; // for N cells, have N+1 nodes.
  int nz = SimPM.NG[ZZ]+1; // for N cells, have N+1 nodes.
  
  //node_coords = mem.myalloc(node_coords,ndim);
  if (silo_dtype==DB_FLOAT) {
    //
    // Allocate memory for node_coords, and set pointers.
    //
    float **d = 0;
    float *posx=0, *posy=0, *posz=0;

    if (node_coords) {
      //d = reinterpret_cast<float **>(node_coords);
      posx = reinterpret_cast<float *>(nodex);
      posy = reinterpret_cast<float *>(nodey);
      posz = reinterpret_cast<float *>(nodez);
    }
    else {
      d = mem.myalloc(d,ndim);
      node_coords = reinterpret_cast<void **>(d);
      posx = mem.myalloc(posx,nx);
      if (ndim>1) posy = mem.myalloc(posy,ny);
      if (ndim>2) posz = mem.myalloc(posz,nz);
    }

    //
    // Assign data for nodex, nodey, nodez for this grid
    //
    for (int i=0;i<nx;i++)
      posx[i] = static_cast<float>(grid->Xmin(XX)+i*dx);
    nodex = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;
    if (ndim>1) {
      for (int i=0;i<ny;i++)
        posy[i] = static_cast<float>(grid->Xmin(YY)+i*dx);
      nodey = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim>2) {
      for (int i=0;i<nz;i++)
        posz[i] = static_cast<float>(grid->Xmin(ZZ)+i*dx);
      nodez = reinterpret_cast<void *>(posz);
      node_coords[ZZ] = nodez;
    }
  }
  else {
    //
    // Allocate memory for node_coords, and set pointers.
    //
    double **d=0;
    double *posx=0, *posy=0, *posz=0;

    if (node_coords) {
      //d = reinterpret_cast<float **>(node_coords);
      posx = reinterpret_cast<double *>(nodex);
      posy = reinterpret_cast<double *>(nodey);
      posz = reinterpret_cast<double *>(nodez);
    }
    else {
      d = mem.myalloc(d,ndim);
      node_coords = reinterpret_cast<void **>(d);
      posx = mem.myalloc(posx,nx);
      if (ndim>1) posy = mem.myalloc(posy,ny);
      if (ndim>2) posz = mem.myalloc(posz,nz);
    }
    d = mem.myalloc(d,ndim);
    node_coords = reinterpret_cast<void **>(d);
    //
    // Assign data for nodex, nodey, nodez for this grid
    //
    for (int i=0;i<nx;i++)
      posx[i] = static_cast<double>(grid->Xmin(XX)+i*dx);
    nodex = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;

    if (ndim>1) {
      for (int i=0;i<ny;i++)
        posy[i] = static_cast<double>(grid->Xmin(YY)+i*dx);
      nodey = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim>2) {
      for (int i=0;i<nz;i++)
        posz[i] = static_cast<double>(grid->Xmin(ZZ)+i*dx);
      nodez = reinterpret_cast<void *>(posz);
      node_coords[ZZ] = nodez;
    }
  }

  nodedims[0] = nx;
  zonedims[0] = nx-1;

  if (ndim>1) {
    nodedims[1] = ny;
    zonedims[1] = ny-1;
  }
  if (ndim>2) {
    nodedims[2] = nz;
    zonedims[2] = nz-1;
  }


  int nopts=4;
  dataio_silo::GridOpts = DBMakeOptlist(nopts);
  if      (SimPM.coord_sys==COORD_CRT) silo_coordsys=DB_CARTESIAN;
  else if (SimPM.coord_sys==COORD_CYL) silo_coordsys=DB_CYLINDRICAL;
  else if (SimPM.coord_sys==COORD_SPH) silo_coordsys=DB_SPHERICAL;
  else rep.error("bad coord system",SimPM.coord_sys);
  DBAddOption(GridOpts,DBOPT_COORDSYS,&silo_coordsys);
  DBAddOption(GridOpts,DBOPT_DTIME,&SimPM.simtime);
  DBAddOption(GridOpts,DBOPT_CYCLE,&SimPM.timestep);
  DBAddOption(GridOpts,DBOPT_NSPACE,&SimPM.ndim);
  // labels don't seem to display right in VisIt...
  //char s[strlength];
  //strcpy(s,"XXXX"); DBAddOption(GridOpts,DBOPT_XLABEL,s);
  //strcpy(s,"YYYY"); DBAddOption(GridOpts,DBOPT_YLABEL,s);
  //strcpy(s,"ZZZZ"); DBAddOption(GridOpts,DBOPT_ZLABEL,s);
  //char temp[strlength];
  //strcpy(temp,uc.length.c_str());
  //DBAddOption(GridOpts,DBOPT_XUNITS,temp);
  //DBAddOption(GridOpts,DBOPT_YUNITS,temp);
  //DBAddOption(GridOpts,DBOPT_ZUNITS,temp);

  have_setup_gridinfo = true;
  return 0;
}



// ##################################################################
// ##################################################################


#endif // SILO

