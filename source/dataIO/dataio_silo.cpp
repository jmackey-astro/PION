///
/// \file dataio_silo.cc
/// This has the class definitions for the serial I/O class
/// dataio_silo. Note that the variable SILO must be set in
/// the Makefile for this code to be compiled.
///
///  - JM 2009-12-15: Added possibility to write curl(B) in 2D MHD sims.
///    This is only needed if using cgs since VisIt does single precision
///    calculations and seems to have trouble with small numbers.  Should
///    be possible to switch it off easily, but I haven't done that yet.
///  - 2010-02-03 JM: small changes to fix compiler warnings -- use
///     GS.equalD() in a few places instead of testing for equality.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010-07-20/22 JM: Work on new dataio structure with a list of
///    parameters to read and write.  read_simulation_parameters()
///    function replaces read_header() so I will be able to delete a
///    lot of code!  Same for write_header().
/// - 2010.07.23 JM: removed obselete read_header(),
///    write_header() functions.
/// - 2010.10.01 JM: Spherical coordinates added.
///    Got rid of testing myalloc/myfree commands.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
///    Also replaced endl with c-style line-break for JUROPA.
/// - 2010.11.12 JM: Changed ->col to use cell interface for
///   extra_data. (2010.11.15 JM fixed bug introduced here!)
/// - 2011.02.25 JM: removed HCORR ifdef around new code.
/// - 2011.03.01 JM: Added outputting of diffuse-RT column density when RT_TESTING.
/// - 2011.03.02 JM: Better support for tracer variables (with or without names).
///    Can handle an arbitrary number of tracers now.
/// - 2011.03.21 JM: Updated column-density variables for new cell interface
///    functions.
/// - 2011.03.22 JM: Parallel code now uses serial setup_write_vars().  I put
///    extra data (Ptot,divB,T) in an #ifdef SERIAL so it is not written out
///    for parallel code.
/// - 2011.06.02 JM: Added WriteHeader() function so I can over-write header
///    parameters and restart a simulation with e.g. different microphysics.
/// - 2011.10.14 JM: commented out RT_DIFF
/// - 2012.03.01 JM: Added spacing between functions.
/// - 2013.02.07 JM: Made code less verbose.
/// - 2013.08.20 JM: Modified cell_interface for optical depth vars.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.01.28 JM: Removed parallel code (new class).
/// - 2015.06.13 JM: Changed datatype (FLOAT/DOUBLE) to a runtime
///    parameter, set in the constructor.  (More 18.06)
/// - 2015.07.07 JM: New trtype array structure in constructor.


#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#ifdef SILO

#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#include "raytracing/raytracer_base.h"
#endif // TESTING

#include "microphysics/microphysics_base.h"

#ifdef RT_TESTING_OUTPUTCOL
#include "raytracing/raytracer_base.h"
#endif // RT_TESTING_OUTPUTCOL

#include "dataio_silo.h"
#include <cstring>
#include <sstream>
using namespace std;




// ##################################################################
// ##################################################################


dataio_silo::dataio_silo(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      std::string dtype // read/write either FLOAT or DOUBLE to/from file
      )
: DataIOBase(SimPM)
{
#ifdef TESTING
  cout <<"setting up dataio_silo class.\n";
#endif
  dataio_silo::eqn = 0;
  dataio_silo::gp  = 0;
  silofile.erase();
  ndim = -1;
  nodedims=0;
  zonedims=0;
  have_setup_gridinfo=false;
  have_setup_writevars=false;
  varnames.clear();
  readvars.clear();
  silo_filetype=SILO_FILETYPE;
  strlength=256;
  db_ptr = mem.myalloc(db_ptr,1);
  GridOpts=0;

  //
  // These are either double or float, depending on dtype value.
  // So here they are set to void pointers.
  //
  node_coords=0;
  nodex=nodey=nodez=0;
  data0=data1=data2=0;
  vec_data=0;

  //
  // choose what sort of data to read/write
  //
  if (dtype=="FLOAT") {
    silo_dtype  = DB_FLOAT;    ///< defined in <silo.h> (ext.lib.)
  }
  else if (dtype=="DOUBLE") {
    silo_dtype  = DB_DOUBLE;   ///< defined in <silo.h> (ext.lib.)
  }
  else {
    rep.error("Bad datatype for silo initialisation",dtype);
  }
  
  return;
}




// ##################################################################
// ##################################################################


dataio_silo::~dataio_silo()
{
#ifdef TESTING
  cout <<"deleting dataio_silo class.\n";
#endif
  dataio_silo::eqn = 0;
  dataio_silo::gp  = 0;
  silofile.erase();
  varnames.clear();
  readvars.clear();

  nodedims = mem.myfree(nodedims);
  zonedims = mem.myfree(zonedims);
  //
  // freeing memory for void arrays:
  //
  if (silo_dtype == DB_FLOAT) {
    mem.myfree(reinterpret_cast<float *>(nodex));
    mem.myfree(reinterpret_cast<float *>(nodey));
    mem.myfree(reinterpret_cast<float *>(nodez));
    mem.myfree(reinterpret_cast<float **>(node_coords));
  }
  if (silo_dtype == DB_DOUBLE) {
    mem.myfree(reinterpret_cast<double *>(nodex));
    mem.myfree(reinterpret_cast<double *>(nodey));
    mem.myfree(reinterpret_cast<double *>(nodez));
    mem.myfree(reinterpret_cast<double **>(node_coords));
  }
  nodex    = 0;
  nodey    = 0;
  nodez    = 0;
  node_coords = 0;
  delete_data_arrays();

  // Have to check if we used the grid options for writing data.
  if(GridOpts) {
    DBClearOptlist(GridOpts);
    DBFreeOptlist(GridOpts);
  }
  //  *db_ptr=0;
  db_ptr = mem.myfree(db_ptr);
}




// ##################################################################
// ##################################################################


void dataio_silo::SetSolver(FV_solver_base *solver)
{
#ifdef TESTING
  cout <<"dataio_silo::SetSolver() Setting solver pointer.\n";
#endif
  dataio_silo::eqn = solver;
}




// ##################################################################
// ##################################################################



int dataio_silo::WriteHeader(
      const string overwritefile, ///< file to write to (full, exact filename).
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // File must already exist, so we call Open() with the APPEND setting
  // so we can overwrite existing parameters.
  //
  *db_ptr = 0;
  *db_ptr = DBOpen(overwritefile.c_str(), DB_UNKNOWN, DB_APPEND);
  if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);

  //
  // Create header directory if it doesn't exist, move to it,
  // and write the sim parameters.
  //
  int err = DBSetDir(*db_ptr,"/header");
  if (err==-1) {
    DBSetDir(*db_ptr,"/");
    DBMkDir(*db_ptr,"header");
    DBSetDir(*db_ptr,"/header");
  }
  err = write_simulation_parameters(SimPM);
  if (err)
    rep.error("dataio_silo::WriteHeader() error writing header to silo file",err);

  DBClose(*db_ptr); //*db_ptr=0; 
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::OutputData(
      const string outfile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      const long int file_counter   ///< number to stamp file with (e.g. timestep)
      )
{
  int err=0;

  if (!cg[0])
    rep.error("dataio_silo::OutputData() null pointer to grid!",cg[0]);
  dataio_silo::gp = cg[0];

  if (!have_setup_writevars) {
    // set what data to write to the mesh.
    err = dataio_silo::setup_write_variables(SimPM);
    if (err)
      rep.error("dataio_silo::OutputData() write vars",err);
  }

  err = dataio_silo::choose_filename(outfile, file_counter);
  if (err) {
    cerr<<"dataio_silo::OutputData() error choosing filename.\n";
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
  *db_ptr = DBCreate(silofile.c_str(), DB_CLOBBER, DB_LOCAL, "PION data", silo_filetype);
  if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);
  //cout <<"\tdb_ptr="<<db_ptr<<"\n";
  //cout <<"\t*db_ptr="<<*db_ptr<<"\n";

  if (!have_setup_gridinfo) {
    // set grid properties for quadmesh
    err = dataio_silo::setup_grid_properties(gp, SimPM);
    if (err)
      rep.error("dataio_silo::OutputData() error setting up grid_props", err);
  }

  //
  // now write the simulation parameters to the header part of the file.
  //
  DBSetDir(*db_ptr,"/");
  DBMkDir(*db_ptr,"header");
  DBSetDir(*db_ptr,"/header");
  err = write_simulation_parameters(SimPM);
  if (err)
    rep.error("dataio_silo::OutputData() error writing header to silo file",err);

  //
  // Create data directory, generate the mesh in the file, and then write each
  // variable in turn to the mesh.
  //
  DBSetDir(*db_ptr,"/");
  string meshname="UniformGrid";
  err = dataio_silo::generate_quadmesh(*db_ptr, meshname,SimPM);
  if (err)
    rep.error("dataio_silo::OutputData() error writing quadmesh to silo file",err);

  dataio_silo::create_data_arrays(SimPM);
  for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
    err = dataio_silo::write_variable2mesh(SimPM, *db_ptr, meshname, (*i));
    if (err)
      rep.error("dataio_silo::OutputData() error writing variable",(*i));
  }
  dataio_silo::delete_data_arrays();
  DBSetDir(*db_ptr,"/");
  DBClose(*db_ptr); //*db_ptr=0;

  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::ReadHeader(
      string infile, ///< file to read from
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  int err=0;
  silofile=infile;
#ifdef TESTING
  cout <<"Reading Header from file: "<<silofile<<"\n";
#endif

  // Create file
  //*db_ptr=0;
  if (silofile.size()>=strlength-1) {
    rep.error("string too large",silofile);
  }
  char temp[strlength];
  strcpy(temp,silofile.c_str());
  *db_ptr = DBOpen(temp, DB_UNKNOWN, DB_READ);
  if (!(*db_ptr)) {
    rep.error("open silo file failed.",*db_ptr);
  }

  DBSetDir(*db_ptr,"/header");
  err = read_simulation_parameters(SimPM);
  if (err) {
    rep.error("dataio_silo::ReadHeader() error reading header \
               from silo file",err);
  }
  dataio_silo::ndim = SimPM.ndim;

  DBClose(*db_ptr); //*db_ptr=0;
#ifdef TESTING
  cout <<"FINISHED reading Header from file: "<<silofile<<"\n";
#endif
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::ReadData(
      string infile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (!cg[0])
    rep.error("dataio_silo::ReadData() null pointer to grid!",cg[0]);
  dataio_silo::gp = cg[0];
  silofile=infile;

  int err=0;
  if (!have_setup_gridinfo) {
    // set grid properties for quadmesh,
    // also check grid pointer is not null.
    err = dataio_silo::setup_grid_properties(gp, SimPM);
    if (err)
      rep.error("dataio_silo::ReadData() error setting up grid_props", err);
  }

  *db_ptr = DBOpen(silofile.c_str(), DB_UNKNOWN, DB_READ);
  if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);

  int ftype = DBGetDriverType(*db_ptr);
  if (ftype==DB_HDF5) {
    //cout <<"READING HDF5 FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    //
    // I think you don't need to get the compression for reading a file...
    // Either way this code is buggy and it works fine without it.
    //
    //char *ccc = DBGetCompression();
    //char compress[strlength];
    //strcpy(compress,ccc);
    //cout <<"compression="<<compress<<"\n";
    //DBSetCompression(compress);
    //
    // Not sure if we need this either for reading, but set it anyway.
    //
    int friendly=DBGetFriendlyHDF5Names();
    DBSetFriendlyHDF5Names(friendly);
  }

// This is another way of reading the data, but not really better.
//   DBtoc *toc=0;
//   DBSetDir(db_ptr,"/data");
//   toc = DBGetToc(db_ptr);
//   if (!toc) rep.error("didn't get table of contents",toc);
//   cout <<"READING TOC:  nvar="<<toc->nqvar<<"\n";
//   for(int i=0;i<toc->nqvar;i++){ 
//     if(i==0){
//       printf("%d Quad mesh variables \n",toc->nqvar);
//     }
//     printf(" %s \n",toc->qvar_names[i]);
//   }
//   if (toc->nqvar <SimPM.nvar)
//     rep.error("not enough variables in file!",SimPM.nvar-toc->nqvar);
//   // now read each variable in turn from the mesh
//   for (int i=0;i<toc->nqvar;i+) {
//     err = dataio_silo::read_variable2grid(db_ptr, meshname, toc->qvar_names[i]);
//     if (err)
//       rep.error("dataio_silo::ReadData() error reading variable",toc->qvar_names[i]);
//   }
//   DBSetDir(db_ptr,"/");
  

  err = set_readvars(SimPM);
  if (err) rep.error("failed to set readvars in ReadData",err);

  //DBSetDir(*db_ptr,"/data");
  DBSetDir(*db_ptr,"/");
  string meshname="UniformGrid";

  // now read each variable in turn from the mesh
  for (std::vector<string>::iterator i=readvars.begin(); i!=readvars.end(); ++i) {
    err = dataio_silo::read_variable2grid(SimPM, *db_ptr, meshname, (*i), SimPM.Ncell);
    if (err)
      rep.error("dataio_silo::ReadData() error reading variable",(*i));
  }
  DBSetDir(*db_ptr,"/");

  DBClose(*db_ptr); //*db_ptr=0; 

  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::choose_filename(
      const string codefile,
      const int counter
      )
{
  //
  // for serial files this is easy -- just codefile.cycle.silo. 
  //
  // If codefile already contains .silo we assume this is the full
  // filename i want to write.  Alternatively if counter<0 then it is
  // assumed we just want to append .silo to the codefile string.
  // Finally if counter>0, then we write to codefile.counter.silo
  //
  if (codefile.find(".silo")!=string::npos) {
    // we're hopefully writing initial conditions, so don't append cycle.
    silofile=codefile;
  }
  else if (counter <0) {
    ostringstream temp; temp.str("");
    temp << codefile.c_str() <<".silo";
    silofile=temp.str();
    temp.str("");
  }
  else {
    ostringstream temp; temp.str("");
    temp << codefile.c_str() <<".";
    temp.width(Ndigits); temp.fill('0');
    temp << counter << ".silo";
    silofile=temp.str();
    temp.str("");
  }
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::setup_grid_properties(
      class GridBaseClass *grid, 
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  // set grid parameters -- UNIFORM FIXED GRID
  if (!grid)
    rep.error("dataio_silo::setup_grid_properties() null grid pointer!",grid);
  double dx=grid->DX();
  if (node_coords || nodedims || zonedims ||
      nodex || nodey || nodez) {
    cerr<<"Have already setup variables for grid props! ";
    cerr<<"You shouldn't have called me! crazy fool...\n";
    return 1000000;
  }

  dataio_silo::ndim = SimPM.ndim;
  dataio_silo::vec_length = SimPM.eqnNDim;
  //  dataio_silo::vec_length = ndim;
  //  cout <<"************VEC_LENGTH="<<vec_length<<"\n";

  nodedims = mem.myalloc(nodedims,ndim);
  zonedims = mem.myalloc(zonedims,ndim);

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
    d = mem.myalloc(d,ndim);
    node_coords = reinterpret_cast<void **>(d);
    //
    // Allocate memory for nodex, nodey, nodez
    //
    float *posx=0, *posy=0, *posz=0;
    posx = mem.myalloc(posx,nx);
    for (int i=0;i<nx;i++)
      posx[i] = static_cast<float>(SimPM.Xmin[XX]+i*dx);
    nodex = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;
    if (ndim>1) {
      posy = mem.myalloc(posy,ny);
      for (int i=0;i<ny;i++)
        posy[i] = static_cast<float>(SimPM.Xmin[YY]+i*dx);
      nodey = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim>2) {
      posz = mem.myalloc(posz,nz);
      for (int i=0;i<nz;i++)
        posz[i] = static_cast<float>(SimPM.Xmin[ZZ]+i*dx);
      nodez = reinterpret_cast<void *>(posz);
      node_coords[ZZ] = nodez;
    }
  }
  else {
    //
    // Allocate memory for node_coords, and set pointers.
    //
    double **d=0;
    d = mem.myalloc(d,ndim);
    node_coords = reinterpret_cast<void **>(d);
    //
    // Allocate memory for nodex, nodey, nodez
    //
    double *posx=0, *posy=0, *posz=0;
    posx = mem.myalloc(posx,nx);
    for (int i=0;i<nx;i++)
      posx[i] = static_cast<double>(SimPM.Xmin[XX]+i*dx);
    nodex = reinterpret_cast<void *>(posx);
    node_coords[XX] = nodex;

    if (ndim>1) {
      posy = mem.myalloc(posy,ny);
      for (int i=0;i<ny;i++)
        posy[i] = static_cast<double>(SimPM.Xmin[YY]+i*dx);
      nodey = reinterpret_cast<void *>(posy);
      node_coords[YY] = nodey;
    }
    if (ndim>2) {
      posz = mem.myalloc(posz,nz);
      for (int i=0;i<nz;i++)
        posz[i] = static_cast<double>(SimPM.Xmin[ZZ]+i*dx);
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



int dataio_silo::setup_write_variables(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (!varnames.empty())
    rep.error("dataio_silo::setup_write_variables() variable list not empty!",varnames.size());

  // select variables based on what equations we are using.
  // All equations have density, pressure, and velocity.
  varnames.push_back("Density");
  varnames.push_back("Pressure");
#if defined (SILO_VECTORS)
  varnames.push_back("Velocity");
#elif defined (SILO_SCALARS)
  varnames.push_back("VelocityX");
  varnames.push_back("VelocityY");
  varnames.push_back("VelocityZ");
#else
#error "need to define SILO_SCALARS OR SILO_VECTORS"
#endif

  // MHD has B-field, and maybe Psi for glm-mhd
  if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQFCD || SimPM.eqntype==EQGLM) {
#if defined (SILO_VECTORS)
    varnames.push_back("MagneticField");
#elif defined (SILO_SCALARS)
    varnames.push_back("MagneticFieldX");
    varnames.push_back("MagneticFieldY");
    varnames.push_back("MagneticFieldZ");
#else
#error "need to define SILO_SCALARS OR SILO_VECTORS"
#endif
    if (SimPM.eqntype==EQGLM) 
      varnames.push_back("glmPSI");

#ifdef SERIAL
    //
    // if equations are set up, can output divB and Ptot
    // If doing a 2D sim, also output Curl(B) which is a 
    // scalar for a 2D field.
    // WE ONLY WANT TO DO THIS FOR SERIAL CODE (SAVE DISK SPACE FOR BIG SIMS).
    //
    if (dataio_silo::eqn!=0) {
      varnames.push_back("DivB");
      varnames.push_back("Ptot");
      if (ndim==2 && SimPM.coord_sys==COORD_CRT) varnames.push_back("CurlB");
    }
#endif // SERIAL
  }

  //#ifdef SERIAL
  // if equations are set up, can get temperature/internal energy.
  if (dataio_silo::eqn!=0) {
    if (MP) varnames.push_back("Temperature");
    else varnames.push_back("InternalEnergy");
  }
  //#endif // SERIAL

#ifdef COUNT_ENERGETICS
  //
  // output extra variables if doing ionising-RT and we want to look at energetics.
  //
  if (SimPM.RS.Nsources>0 && SimPM.EP.phot_ionisation) {
    varnames.push_back("ci_cooling");
    varnames.push_back("rr_cooling");
    varnames.push_back("fn_cooling");
    varnames.push_back("pi_heating");
    varnames.push_back("ci_rate");
    varnames.push_back("rr_rate");
    varnames.push_back("pi_rate");
    varnames.push_back("tot_heating");
    varnames.push_back("tot_cooling");
    varnames.push_back("net_heating");
    varnames.push_back("cooling_time");
    varnames.push_back("recomb_time");
  }
#endif

#ifdef RT_TESTING_OUTPUTCOL
  //
  // If testing diffuse/ionising RT, output extra column density data:
  //
  for (int v=0;v<SimPM.RS.Nsources;v++) {
    ostringstream var;
    for (int iT=0; iT<SimPM.RS.sources[v].NTau; iT++) {
      var.str("");
      switch (SimPM.RS.sources[v].type) {
        case RT_SRC_SINGLE:
        var << "Col_Src_" << v;
        var << "_T"<<iT;
        varnames.push_back(var.str());
        break;

        case RT_SRC_DIFFUSE:
        var << "ColDiff_" << v;
        var << "_T"<<iT;
        varnames.push_back(var.str());
        break;

        default:
        rep.error("Bad radiation source type",SimPM.RS.sources[v].type);
        break;
      } // switch
    } // loop over Tau vars for source
  } // loop over Nsources
#endif // RT_TESTING_OUTPUTCOL

  //
  // if there are any tracer variables, get their names from SimPM.tracers,
  // if it has the info.  All we really need is "TrXXX", where the Xs are 
  // the tracer number.
  //
  if (SimPM.ntracer>0) {
    string s;
    ostringstream temp;
    for (int i=0; i<SimPM.ntracer; i++) {
      s.erase();
      temp.str("");
      temp<< "Tr";
      temp.width(3); temp.fill('0'); temp << i;
      temp <<"_"<< SimPM.tracers[i];

      s=temp.str();
      // replace "+" with "p", and "-" with "m"
      string::size_type p=s.find("+");
      if (p!=string::npos) s.replace(p,1,"p");
      p=s.find("-");
      if (p!=string::npos) s.replace(p,1,"m");
      //      cout <<"tracer = "<<s<<"\n";
      varnames.push_back(s);
    }
  } //tracers

#ifdef TESTING
  cout <<"list of vars: ";
  for (unsigned int i=0; i<varnames.size(); i++)
    cout <<varnames[i]<<"  ";
  cout <<"\n";
#endif //TESTING
  have_setup_writevars=true;
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_header_param(class pm_base *p)
{
  int err=0;
  int i=p->type;
  if      (i==MY_INT) {
    int dim1=1;
    int *x = static_cast<int *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim1,1,DB_INT);
  }
  else if (i==MY_DOUBLE) {
    int dim1=1;
    double *x = static_cast<double *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim1,1,DB_DOUBLE);
  }
  else if (i==MY_FLOAT) {
    int dim1=1;
    float *x = static_cast<float *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim1,1,DB_FLOAT);
  }
  else if (i==MY_LONG) {
    int dim1=1;
    long int *x = static_cast<long int *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim1,1,DB_LONG);
  }
  else if (i==MY_STRING) {
    //
    // strings are harder -- need to get pointer and copy to char[]
    //
    int dim2=strlength; 
    string x(*(static_cast<string *>(p->get_ptr())));
    char temp[strlength]; strcpy(temp,x.c_str());
    err += DBWrite(*db_ptr, p->name.c_str(), temp, &dim2,1,DB_CHAR);
  }
  else if (i==MY_DDIMARR) {
    int dim3=MAX_DIM;
    double *x = static_cast<double *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim3,1,DB_DOUBLE);
  }
  else if (i==MY_IDIMARR) {
    int dim3=MAX_DIM;
    int *x = static_cast<int *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dim3,1,DB_INT);
  }
  else if (i==MY_DVARARR) {
    int dimN=MAX_NVAR;
    double *x = static_cast<double *>(p->get_ptr());
    err += DBWrite(*db_ptr, p->name.c_str(), x, &dimN,1,DB_DOUBLE);
  }

  if (err) {
    cout <<"\t"<<p->name<<":  ERROR WRITING VAR!\n";
  }
  // else {
  //   cout <<"\t"<<p->name<<":  "; p->show_val(); cout <<"\n";
  // }

  return err;
}



// ##################################################################
// ##################################################################




int dataio_silo::generate_quadmesh(
      DBfile *dbfile,
      string meshname,
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  DBClearOption(GridOpts,DBOPT_DTIME);
  DBAddOption(GridOpts,DBOPT_DTIME,&SimPM.simtime);
  DBClearOption(GridOpts,DBOPT_CYCLE);
  DBAddOption(GridOpts,DBOPT_CYCLE,&SimPM.timestep);
  
  //DBClearOption(GridOpts,DBOPT_COORDSYS);
  //int csys=0;
  //if      (SimPM.coord_sys==COORD_CRT) csys=DB_CARTESIAN;
  //else if (SimPM.coord_sys==COORD_CYL) csys=DB_CYLINDRICAL;
  //else if (SimPM.coord_sys==COORD_SPH) csys=DB_SPHERICAL;
  //else rep.error("bad coord system",SimPM.coord_sys);
  //DBAddOption(GridOpts,DBOPT_COORDSYS,&csys);
 
  //
  // set coordinate axis names.  This has to be char **, so I can't just
  // send in strings to the silo function.
  //
  char **coordnames=0;
  coordnames = mem.myalloc(coordnames,ndim);
  for (int i=0;i<ndim;i++) {
    coordnames[i] = mem.myalloc(coordnames[i],32);
  }

  std::vector<std::string> s;
  s.push_back("X"); s.push_back("Y"); s.push_back("Z");
  for (int i=0;i<ndim;i++) {
    strcpy(coordnames[i],s[i].c_str());
  }
  //cout <<"coords: "; for (int i=0;i<ndim;i++) cout << coordnames[i]<<"  "; cout <<"\n";

  //  cout <<"dbfile: "<<dbfile<<"\ttemp:"<<temp<<"\tcoordnames:"<<coordnames<<"\n";
  //  cout <<"nodecoords:"<<node_coords<<"\tnodedims:"<<nodedims<<"\tndim:"<<ndim<<"\n";
  //  cout <<"gridopts:"<<GridOpts<<"\n";

  //
  // DBPutQuadmesh requires the data to be (void **), with the actual
  // datatype in silo_dtype.  This is why node_coords is void **.
  //
  int err = DBPutQuadmesh(dbfile, meshname.c_str(), coordnames,
                          node_coords, nodedims, ndim, silo_dtype,
                          DB_COLLINEAR, GridOpts);

  for (int i=0;i<ndim;i++) coordnames[i] = mem.myfree(coordnames[i]);
  coordnames = mem.myfree(coordnames);

  return err;
}



// ##################################################################
// ##################################################################



void dataio_silo::create_data_arrays(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // first check if we have the data arrays set up yet.
  // We need at least one array for a scalar, plus 2 for a vector.
  //
  // data0 is a void pointer, so we need to do different things for
  // float and double data (same for data1,data2).
  //
  if (!data0) {
    if (silo_dtype==DB_FLOAT) {
      float *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data0 = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data0 = reinterpret_cast<void *>(d);
    }
  }

  //
  // If we are only writing scalar data, we don't need data1,data2,vec_data
  // so by setting vec_length=0 they don't get initialised.
  // This saves some memory.
  //
#if defined (SILO_SCALARS)
  vec_length=0;
#endif

  //
  // If we need data1 and data2, and vec_data, create them too.
  //
  if ((vec_length>1) && (!data1)) {
    if (silo_dtype==DB_FLOAT) {
      float *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data1 = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data1 = reinterpret_cast<void *>(d);
    }
    //data1 = mem.myalloc(data1, SimPM.Ncell);
  }
  if ((vec_length>2) && (!data2)) {
    if (silo_dtype==DB_FLOAT) {
      float *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data2 = reinterpret_cast<void *>(d);
    }
    else {
      double *d = 0;
      d = mem.myalloc(d,SimPM.Ncell);
      data2 = reinterpret_cast<void *>(d);
    }
    //data2 = mem.myalloc(data2, SimPM.Ncell);
  }

  if ((vec_length>1) && (!vec_data)) {
    if (silo_dtype==DB_FLOAT) {
      float **d = 0;
      d = mem.myalloc(d,vec_length);
      vec_data = reinterpret_cast<void **>(d);
    }
    else {
      double **d = 0;
      d = mem.myalloc(d,vec_length);
      vec_data = reinterpret_cast<void **>(d);
    }
    //vec_data = mem.myalloc(vec_data, vec_length);
    vec_data[0] = data0;
    if (vec_length>1) vec_data[1] = data1;
    if (vec_length>2) vec_data[2] = data2;
  }

  return;
}



// ##################################################################
// ##################################################################



void dataio_silo::delete_data_arrays()
{
  //
  // freeing memory for void arrays:
  //
  if (silo_dtype == DB_FLOAT) {
    mem.myfree(reinterpret_cast<float *>(data0));
    mem.myfree(reinterpret_cast<float *>(data1));
    mem.myfree(reinterpret_cast<float *>(data2));
    mem.myfree(reinterpret_cast<float **>(vec_data));
  }
  if (silo_dtype == DB_DOUBLE) {
    mem.myfree(reinterpret_cast<double *>(data0));
    mem.myfree(reinterpret_cast<double *>(data1));
    mem.myfree(reinterpret_cast<double *>(data2));
    mem.myfree(reinterpret_cast<double **>(vec_data));
  }
  data0    = 0;
  data1    = 0;
  data2    = 0;
  vec_data = 0;
  return;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_variable2mesh(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      DBfile *dbfile,  ///< pointer to silo file.
      string meshname, ///< name of mesh to write to.
      string variable  ///< variable name to write.
      )
{
  if (!data0) rep.error("allocate data arrays before trying to write data!",data0);
  int err=0;
  
  if (variable=="Velocity" || variable=="MagneticField") {
    // put data into vec array.
    err = get_vector_data_array(variable, SimPM, vec_data);
    if (err) rep.error("failed to get vector data for var.",variable);
    // write data to mesh.
    err = write_vector2mesh(dbfile, meshname, variable, vec_data);
    if (err) rep.error("failed to write vector data for var.",variable);
  } // vector variable

  else {
    // scalar variable, already have array, so just get data and write it.
    err = get_scalar_data_array(variable, SimPM,  data0);
    if (err) rep.error("failed to get scalar data for var.",variable);
    err = write_scalar2mesh(dbfile, meshname, variable, data0);
    if (err) rep.error("failed to write scalar data for var.",variable);
  } // scalar variable
  
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::get_scalar_data_array(
      string variable, ///< variable name to get.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      void *data_array     ///< array to write to.
      )
{
  int v=999;
  if      (variable=="Density")         {v=static_cast<int>(RO);}
  else if (variable=="Pressure")        {v=static_cast<int>(PG);}
  else if (variable=="VelocityX")       {v=static_cast<int>(VX);}
  else if (variable=="VelocityY")       {v=static_cast<int>(VY);}
  else if (variable=="VelocityZ")       {v=static_cast<int>(VZ);}
  else if (variable=="MagneticFieldX")  {v=static_cast<int>(BX);}
  else if (variable=="MagneticFieldY")  {v=static_cast<int>(BY);}
  else if (variable=="MagneticFieldZ")  {v=static_cast<int>(BZ);}
  else if (variable=="glmPSI")     {v=static_cast<int>(SI);}
  //
  // Now loop over up to MAX_NVAR tracers...
  //
  else if (variable.substr(0,2)=="Tr") {
    int itr = atoi(variable.substr(2,3).c_str());
    if (!isfinite(itr) || itr<0 || itr>=MAX_NVAR) {
      rep.error("Bad tracer variable identifier.",variable);
    }
    v = SimPM.ftr +itr;
  }

  else if (variable=="Temperature" ||
	   variable=="InternalEnergy")  {v = -1;}
  else if (variable=="DivB")            {v = -2;}
  else if (variable=="Ptot")            {v = -3;}
  else if (variable=="CurlB")           {v = -5;}
#ifdef COUNT_ENERGETICS
  else if (variable=="ci_cooling")          {v = -105;}
  else if (variable=="ci_rate")             {v = -106;}
  else if (variable=="pi_heating")          {v = -107;}
  else if (variable=="pi_rate")             {v = -108;}
  else if (variable=="rr_cooling")          {v = -109;}
  else if (variable=="rr_rate")             {v = -110;}
  else if (variable=="fn_cooling")          {v = -111;}
  else if (variable=="tot_heating")         {v = -112;}
  else if (variable=="tot_cooling")         {v = -113;}
  else if (variable=="net_heating")         {v = -114;}
  else if (variable=="cooling_time")        {v = -115;}
  else if (variable=="recomb_time")         {v = -116;}
#endif // COUNT_ENERGETICS

#ifdef RT_TESTING_OUTPUTCOL
  //
  // Also pick up optional diffuse-RT column density data (-10 >= v > -20)
  //
  else if (variable.find("ColDiff") !=string::npos) {
    int tdv = atoi(variable.substr(8).c_str());
    if (!isfinite(tdv) || tdv<0 || tdv>9) {
      rep.error("Bad diffuse Column-density identifier.",variable);
    }
    v=-10-tdv;
  }
  //
  // Ionising source Optical depth variable, with the id of the source:
  //
  else if (variable.find("Col_Src") !=string::npos) {
    int tdv = atoi(variable.substr(8).c_str());
    if (!isfinite(tdv) || tdv<0 || tdv>9) {
      rep.error("Bad Ionising source identifier.",variable);
    }
    v=-20-tdv;
  }
  else rep.error("Bad variable requested for dataio_silo::get_scalar_data_array()",variable);
#endif // RT_TESTING_OUTPUTCOL

  //
  // Now pick out the data requested cell by cell, and put it into
  // the 1D array.
  //
  // data_array is a void pointer, so we need a temporary data array
  // for floats and doubles to write the numbers to data_array.
  //
  float  *farr = 0;
  double *darr = 0;
  if (silo_dtype==DB_FLOAT) {
    farr = reinterpret_cast<float *>(data_array);
  }
  else {
    darr = reinterpret_cast<double *>(data_array);
  }

  cell *c = gp->FirstPt(); long int ct=0;
  if (v>=0) {
    if (silo_dtype==DB_FLOAT) {
      do {
        farr[ct] = static_cast<float>(c->P[v]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        darr[ct] = static_cast<double>(c->P[v]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

  else if (v==-1) {
    //
    // internal energy (or temperature if we have microphysics)
    //
    if (MP) {
      if (silo_dtype==DB_FLOAT) {
        do {
          farr[ct] = static_cast<float>(MP->Temperature(c->P,SimPM.gamma));
          ct++;
        } while ( (c=gp->NextPt(c))!=0 );
      }
      else {
        do {
          darr[ct] = static_cast<double>(MP->Temperature(c->P,SimPM.gamma));
          ct++;
        } while ( (c=gp->NextPt(c))!=0 );
      }
    }
    else {
      if (silo_dtype==DB_FLOAT) {
        do {
          farr[ct] = static_cast<float>(eqn->eint(c->P,SimPM.gamma));
          ct++;
        } while ( (c=gp->NextPt(c))!=0 );
      }
      else {
        do {
          darr[ct] = static_cast<double>(eqn->eint(c->P,SimPM.gamma));
          ct++;
        } while ( (c=gp->NextPt(c))!=0 );
      }
    }
  }

  else if (v==-2) { // divB
    int vars[3];
    vars[0] = static_cast<int>(BX);
    vars[1] = static_cast<int>(BY);
    vars[2] = static_cast<int>(BZ);
    if (silo_dtype==DB_FLOAT) {
      do {
        farr[ct] = static_cast<float>(eqn->Divergence(c,0,vars, gp));
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        darr[ct] = static_cast<double>(eqn->Divergence(c,0,vars, gp));
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

  else if (v==-5) { // CurlB (for 2D data only!)
    int vars[3];
    vars[0] = static_cast<int>(BX);
    vars[1] = static_cast<int>(BY);
    vars[2] = static_cast<int>(BZ);
    pion_flt crl[3]; for (int el=0;el<3;el++) crl[el]=0.0;
    if (silo_dtype==DB_FLOAT) {
      do {
        eqn->Curl(c,0,vars, gp, crl);
        farr[ct] = static_cast<float>(crl[2]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        eqn->Curl(c,0,vars, gp, crl);
        darr[ct] = static_cast<double>(crl[2]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

  else if (v==-3) { // total pressure.
    if (silo_dtype==DB_FLOAT) {
      do {
        farr[ct] = static_cast<float>(eqn->Ptot(c->P,SimPM.gamma));
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        darr[ct] = static_cast<double>(eqn->Ptot(c->P,SimPM.gamma));
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

#ifdef RT_TESTING_OUTPUTCOL
  else if (v<=-20 && v>-30) {
    // ionising-RT column density variable.
#ifdef TESTING
    cout <<"writing variable "<<v<<" corresponding to NH0 RT ";
    cout <<"variable "<<variable<<"\n";
#endif
    double Tau[MAX_TAU];
    int col_id=abs(v+20);
    // which Tau variable?  get from string.
    int iT = atoi(variable.substr(11).c_str());
    if (silo_dtype==DB_FLOAT) {
      do {
        CI.get_col(c, col_id, Tau);
        farr[ct] = static_cast<float>(Tau[iT]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        CI.get_col(c, col_id, Tau);
        darr[ct] = static_cast<double>(Tau[iT]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

  else if (v<=-10 && v>-20) {
    // diffuse-RT column density variable.
#ifdef TESTING
    cout <<"writing variable "<<v<<" corresponding to Ntot RT ";
    cout <<"variable "<<variable<<"\n";
#endif
    double Tau[MAX_TAU];
    int col_id=abs(v+10);
    // which Tau variable?  get from string.
    int iT = atoi(variable.substr(11).c_str());
    if (silo_dtype==DB_FLOAT) {
      do {
        CI.get_col(c, col_id, Tau);
        farr[ct] = static_cast<float>(Tau[iT]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        CI.get_col(c, col_id, Tau);
        darr[ct] = static_cast<double>(Tau[iT]);
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }
#endif // RT_TESTING_OUTPUTCOL


#ifdef COUNT_ENERGETICS
  else if (v==-105) {
    //
    // Hardcode this for double arrays (to save coding).
    //
    if (silo_dtype==DB_FLOAT)
      rep.error("(silo) Use double precision for debugging!",silo_dtype);
    do {darr[ct] = c->e.ci_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-106) {
    do {darr[ct] = c->e.ci_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-107) {
    do {darr[ct] = c->e.pi_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-108) {
    do {darr[ct] = c->e.pi_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-109) {
    do {darr[ct] = c->e.rr_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-110) {
    do {darr[ct] = c->e.rr_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-111) {
    do {darr[ct] = c->e.fn_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-112) {
    do {darr[ct] = c->e.tot_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-113) {
    do {darr[ct] = c->e.tot_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-114) {
    do {darr[ct] = c->e.net_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-115) {
    do {darr[ct] = c->e.cooling_time; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-116) {
    do {darr[ct] = c->e.recomb_time; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
    
#endif // COUNT_ENERGETICS

  else rep.error("Don't understand what variable to write.",v);

  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::get_vector_data_array(
      string variable, ///< variable name to get.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      void **buffer  ///< array to write to.
      )
{
  int err=0;
  if (variable=="Velocity") {
    err += get_scalar_data_array("VelocityX", SimPM, buffer[0]);
    if (vec_length>1) err += get_scalar_data_array("VelocityY", SimPM, buffer[1]);
    if (vec_length>2) err += get_scalar_data_array("VelocityZ", SimPM, buffer[2]);
  }
  else if (variable=="MagneticField") {
    err += get_scalar_data_array("MagneticFieldX", SimPM, buffer[0]);
    if (vec_length>1) err += get_scalar_data_array("MagneticFieldY", SimPM, buffer[1]);
    if (vec_length>2) err += get_scalar_data_array("MagneticFieldZ", SimPM, buffer[2]);
  }
  else {
    cerr <<"Don't know variable: "<<variable<<" as a vector array.\n";
    err=99;
  }
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_scalar2mesh(DBfile *dbfile,  ///< silo file pointer.
				   string meshname, ///< mesh name
				   string variable, ///< variable name
				   void *data     ///< pointer to data array.
				   )
{
  //cout <<"writing variable "<<variable<<" to mesh.\n";
  //
  // data has to be passed to the function as void ** in recent
  // versions of silo.  Datatype is specified with silo_dtype.
  //
  int err = DBPutQuadvar1(dbfile, variable.c_str(), meshname.c_str(),
                          data, zonedims, ndim, 0,0, silo_dtype,
                          DB_ZONECENT, 0);
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_vector2mesh(DBfile *dbfile,  ///< silo file pointer.
				   string meshname, ///< mesh name
				   string variable, ///< variable name
				   void **data    ///< pointer to data array.
				   )
{
  int err=0;
  if (variable!="Velocity" && variable!="MagneticField")
    rep.error("Don't know what vector to write!!!",variable);
  //cout <<"writing variable "<<variable<<" to mesh.\n";
  //
  // Set up array of vector element names.  Has to be done with (char **).
  //
  char **vnames =0;
  vnames = mem.myalloc(vnames, vec_length);
  for (int i=0; i<vec_length; i++) 
    vnames[i] = mem.myalloc(vnames[i], strlength);

  for (int i=0; i<vec_length; i++) {
    string temp=variable.c_str()+i;
    strcpy(vnames[i],temp.c_str());
  }

  //
  // data has to be passed to the function as void ** in recent
  // versions of silo.  Datatype is specified with silo_dtype.
  //
  err = DBPutQuadvar(dbfile, variable.c_str(), meshname.c_str(), vec_length, vnames, data, zonedims, ndim, 0,0, silo_dtype, DB_ZONECENT, 0);

  for (int i=0;i<ndim;i++)
    vnames[i] = mem.myfree(vnames[i]);
  vnames = mem.myfree(vnames);

  return err;
}




// ##################################################################
// ##################################################################




/************************************************************
 **************** READ FUNCTIONS ****************************
 ************************************************************/
int dataio_silo::read_header_param(class pm_base *p)
{
  int err=0;

  //
  // the pointer is unavailable explicitly, so we read to a temp var
  // and copy it over.
  //
  int i=p->type;
  if      (i==MY_INT) {
    int x;
    err += DBReadVar(*db_ptr, p->name.c_str(), &x);
    p->assign_val(&x);
  }
  else if (i==MY_DOUBLE) {
    double x;
    err += DBReadVar(*db_ptr, p->name.c_str(), &x);
    p->assign_val(&x);
  }
  else if (i==MY_FLOAT) {
    float x;
    err += DBReadVar(*db_ptr, p->name.c_str(), &x);
    p->assign_val(&x);
  }
  else if (i==MY_LONG) {
    long int x;
    err += DBReadVar(*db_ptr, p->name.c_str(), &x);
    p->assign_val(&x);
  }
  else if (i==MY_STRING) {
    char x[strlength];
    err += DBReadVar(*db_ptr, p->name.c_str(), x);
    string temp(x);
    p->assign_val(&temp);
  }
  else if (i==MY_DDIMARR) {
    double x[MAX_DIM];
    err += DBReadVar(*db_ptr, p->name.c_str(), x);
    p->assign_val(x);
  }
  else if (i==MY_IDIMARR) {
    int x[MAX_DIM];
    err += DBReadVar(*db_ptr, p->name.c_str(), x);
    p->assign_val(x);
  }
  else if (i==MY_DVARARR) {
    double x[MAX_NVAR];
    err += DBReadVar(*db_ptr, p->name.c_str(), x);
    p->assign_val(x);
  }

  if (err) {
    cout <<"\t"<<p->name<<":  ERROR READING VAR!\n";
  }
  // else {
  //   cout <<"\t"<<p->name<<":  "; p->show_val(); cout <<"\n";
  // }
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::set_readvars(
      class SimParams &SimPM
      )
{
  // select variables based on what equations we are using.
  // All equations have density, pressure, and velocity.
  if (!readvars.empty()) {
    //cout <<"dataio_silo::set_readvars() list not empty!! clearing it now.\n";
    readvars.clear();
  }
  readvars.push_back("Density");
  readvars.push_back("Pressure");

#if defined (SILO_SCALARS)
  readvars.push_back("VelocityX");
  readvars.push_back("VelocityY");
  readvars.push_back("VelocityZ");
#elif defined (SILO_VECTORS)
  readvars.push_back("Velocity");
#else
#error "Must have either scalar components or vector variables defined!"
#endif
  
  // MHD has B-field, and maybe Psi for glm-mhd
  if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQFCD || SimPM.eqntype==EQGLM) {
#if defined (SILO_SCALARS)
    readvars.push_back("MagneticFieldX");
    readvars.push_back("MagneticFieldY");
    readvars.push_back("MagneticFieldZ");
#elif defined (SILO_VECTORS)
    readvars.push_back("MagneticField");
#else
#error "Must have either scalar components or vector variables defined!"
#endif
    if (SimPM.eqntype==EQGLM) 
      readvars.push_back("glmPSI");
  }
  //
  // if there are any tracer variables, get their names from SimPM.tracers,
  // if it contains the details.  At least it's name will start with 
  // "TrXXX", which is all we really need.
  //
  if (SimPM.ntracer>0) {
    string s; ostringstream temp;
    for (int i=0; i<SimPM.ntracer; i++) {
      s.erase(); temp.str("");
      temp<< "Tr";
      temp.width(3); temp.fill('0'); temp << i;
      temp <<"_"<< SimPM.tracers[i];

      s=temp.str();
      // replace "+" with "p", and "-" with "m"
      string::size_type p=s.find("+");
      if (p!=string::npos) s.replace(p,1,"p");
      p=s.find("-");
      if (p!=string::npos) s.replace(p,1,"m");
      readvars.push_back(s);
      //      cout <<"tracer = "<<s<<"\n";
    }
  } //tracers
  return 0;
}



// ##################################################################
// ##################################################################




int dataio_silo::read_variable2grid(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      DBfile *dbfile,  ///< pointer to silo file.
      string, ///< name of mesh to read from (can use it for debugging)
      string variable, ///< variable name to read.
      long int npt     ///< number of points we are expecting.
      )
{
  //
  // get data from silo file: DBGetQuadvar will allocate memory and
  // return a pointer with all the data in it.
  //
  DBquadvar *silodata=0;
  silodata = DBGetQuadvar(dbfile,variable.c_str());
  if (!silodata) {
    rep.error("dataio_silo::read_variable2grid() failed to read variable",
              variable);
  }
  if (silodata->nels != npt) {
    rep.error("dataio_silo::read_variable2grid() wrong number of cells",
              silodata->nels-SimPM.Ncell);
  }
  //
  // Check that datatype is what we are expecting!  If not, then 
  // delete data arrays, reset datatype, and re-create data arrays.
  //
  if (silodata->datatype != silo_dtype) {
    cout <<"HMMM: file has datatype "<<silodata->datatype;
    cout <<" but I am trying to read datatype "<<silo_dtype<<"\n";
    cout <<"    DB_INT=16, DB_SHORT=17, DB_LONG=18, DB_FLOAT=19, ";
    cout <<"DB_DOUBLE=20, DB_CHAR=21, DB_LONG_LONG=22, DB_NOTYPE=25\n";
    cout <<"SRAD_read_var2grid() quadvar has type="<<silodata->datatype;
    cout <<" but expecting type="<<silo_dtype<<"\n";
    cout <<"\t... resetting datatype for this file.\n";
    delete_data_arrays();
    silo_dtype = silodata->datatype;
    create_data_arrays(SimPM);
  }

  //
  // Create a pointer to the data in the silo stuct DBquadvar.  This
  // is a void pointer, so I have to reinterpret it to get data that
  // PION can understand.
  //
  float **fdata=0;
  double **ddata=0;
  if (silo_dtype==DB_FLOAT) {
    fdata = reinterpret_cast<float **>(silodata->vals);
  }
  else {
    ddata = reinterpret_cast<double **>(silodata->vals);
  }
  
  //
  // first read in vector data, if it exists in the file.
  //
  if (variable=="Velocity" || variable=="MagneticField") {
    //
    // choose elements of primitive variables for Velocity or B-Field
    //
    int v1,v2,v3;
    if (variable=="Velocity") {
      v1=VX;
      v2=VY;
      v3=VZ;
    }
    else {
      v1=BX;
      v2=BY;
      v3=BZ;
    }
    //    cout <<"name: "<<silodata->name<<"\tnels="<<silodata->nels<<"\n";
    //    cout <<"ndims: "<<silodata->ndims<<"\tnvals: "<<silodata->nvals<<"\n";
    //cout <<"reading variable "<<variable<<" into element "<<v1<<" of state vec.\n";
    cell *c=gp->FirstPt();
    long int ct=0;

    if (silo_dtype==DB_FLOAT) {
      do {
        c->P[v1] = fdata[0][ct];
        c->P[v2] = fdata[1][ct];
        c->P[v3] = fdata[2][ct];
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        c->P[v1] = ddata[0][ct];
        c->P[v2] = ddata[1][ct];
        c->P[v3] = ddata[2][ct];
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    if (ct != npt) rep.error("wrong number of points read for vector variable",ct-npt);
  } // vector variable

  //
  // Read in scalar data from the file.
  //
  else {
    int v1=0;
    if      (variable=="Density")         v1=RO;
    else if (variable=="Pressure")        v1=PG;
    else if (variable=="VelocityX")       v1=VX;
    else if (variable=="VelocityY")       v1=VY;
    else if (variable=="VelocityZ")       v1=VZ;
    else if (variable=="MagneticFieldX")  v1=BX;
    else if (variable=="MagneticFieldY")  v1=BY;
    else if (variable=="MagneticFieldZ")  v1=BZ;
    else if (variable=="glmPSI")          v1=SI;
    //
    // Now loop over up to MAX_NVAR tracers...
    //
    else if (variable.substr(0,2)=="Tr") {
      int itr = atoi(variable.substr(2,3).c_str());
      if (!isfinite(itr) || itr<0 || itr>=MAX_NVAR) {
        rep.error("Bad tracer variable identifier.",variable);
      }
      v1 = SimPM.ftr +itr;
    }
    else rep.error("what var to read???",variable);

    //cout <<"reading variable "<<variable<<" into element "<<v1<<" of state vec.\n";
    cell *c=gp->FirstPt();
    long int ct=0;
    if (silo_dtype==DB_FLOAT) {
      do {
        c->P[v1] = fdata[0][ct];
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
        c->P[v1] = ddata[0][ct];
        ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
    if (ct != npt)
      rep.error("wrong number of points read for scalar variable",ct-npt);
  } // scalar variable

  //  cout <<"Read variable "<<variable<<"\n";
  DBFreeQuadvar(silodata); //silodata=0;
  fdata=0;
  ddata=0;
  return 0;
}



// ##################################################################
// ##################################################################




#endif // if SILO
