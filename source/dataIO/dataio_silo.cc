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
#ifdef SILO

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "dataio_silo.h"
#include <cstring>
#include <sstream>
using namespace std;




// ##################################################################
// ##################################################################


dataio_silo::dataio_silo()
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
  node_coords=0;
  nodex=nodey=nodez=0;
  have_setup_gridinfo=false;
  have_setup_writevars=false;
  varnames.clear();
  readvars.clear();
  data0=data1=data2=0;
  vec_data=0;
  silo_filetype=SILO_FILETYPE;
  strlength=256;
  db_ptr = mem.myalloc(db_ptr,1);
  GridOpts=0;
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
  nodex    = mem.myfree(nodex   );
  nodey    = mem.myfree(nodey   );
  nodez    = mem.myfree(nodez   );
  node_coords = mem.myfree(node_coords);
  data0    = mem.myfree(data0   );
  data1    = mem.myfree(data1   );
  data2    = mem.myfree(data2   );
  vec_data = mem.myfree(vec_data);

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
          const string overwritefile ///< file to write to (full, exact filename).
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
  err = write_simulation_parameters();
  if (err)
    rep.error("dataio_silo::OutputData() error writing header to silo file",err);

  DBClose(*db_ptr); //*db_ptr=0; 
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::OutputData(
        const string outfile,
        class GridBaseClass *cg,
        const long int file_counter   ///< number to stamp file with (e.g. timestep)
        )
{
  if (!cg)
    rep.error("dataio_silo::OutputData() null pointer to grid!",cg);
  dataio_silo::gp = cg;

  int err=0;

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
  *db_ptr = DBCreate(silofile.c_str(), DB_CLOBBER, DB_LOCAL, "JM's astro code data", silo_filetype);
  if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);
  //cout <<"\tdb_ptr="<<db_ptr<<"\n";
  //cout <<"\t*db_ptr="<<*db_ptr<<"\n";

  if (!have_setup_gridinfo) {
    // set grid properties for quadmesh
    err = dataio_silo::setup_grid_properties(gp);
    if (err)
      rep.error("dataio_silo::OutputData() error setting up grid_props", err);
  }
  if (!have_setup_writevars) {
    // set what data to write to the mesh.
    err = dataio_silo::setup_write_variables();
    if (err)
      rep.error("dataio_silo::OutputData() error settting up variables to write to",err);
  }

  //
  // now write the simulation parameters to the header part of the file.
  //
  DBSetDir(*db_ptr,"/");
  DBMkDir(*db_ptr,"header");
  DBSetDir(*db_ptr,"/header");
  err = write_simulation_parameters();
  if (err)
    rep.error("dataio_silo::OutputData() error writing header to silo file",err);

  //
  // Create data directory, generate the mesh in the file, and then write each
  // variable in turn to the mesh.
  //
  DBSetDir(*db_ptr,"/");
  string meshname="UniformGrid";
  err = dataio_silo::generate_quadmesh(*db_ptr, meshname);
  if (err)
    rep.error("dataio_silo::OutputData() error writing quadmesh to silo file",err);

  dataio_silo::create_data_arrays();
  for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
    err = dataio_silo::write_variable2mesh(*db_ptr, meshname, (*i));
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



int dataio_silo::ReadHeader(string infile ///< file to read from
			   )
{
  int err=0;
  silofile=infile;
#ifdef TESTING
  cout <<"Reading Header from file: "<<silofile<<"\n";
#endif

  // Create file
  //*db_ptr=0;
  if (silofile.size()>=strlength-1) rep.error("string too large",silofile);
  char temp[strlength]; strcpy(temp,silofile.c_str());
  *db_ptr = DBOpen(temp, DB_UNKNOWN, DB_READ);
  if (!(*db_ptr)) rep.error("open silo file failed.",*db_ptr);

  DBSetDir(*db_ptr,"/header");
  err = read_simulation_parameters();
  if (err)
    rep.error("dataio_silo::ReadHeader() error reading header from silo file",err);
  dataio_silo::ndim = SimPM.ndim;

  DBClose(*db_ptr); //*db_ptr=0;
#ifdef TESTING
  cout <<"FINISHED reading Header from file: "<<silofile<<"\n";
#endif
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::ReadData(string infile,
			  class GridBaseClass *cg
			  )
{
  if (!cg)
    rep.error("dataio_silo::ReadData() null pointer to grid!",cg);
  dataio_silo::gp = cg;
  silofile=infile;

  int err=0;
  if (!have_setup_gridinfo) {
    // set grid properties for quadmesh,
    // also check grid pointer is not null.
    err = dataio_silo::setup_grid_properties(gp);
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
  

  err = set_readvars(SimPM.eqntype);
  if (err) rep.error("failed to set readvars in ReadData",err);

  //DBSetDir(*db_ptr,"/data");
  DBSetDir(*db_ptr,"/");
  string meshname="UniformGrid";

  // now read each variable in turn from the mesh
  for (std::vector<string>::iterator i=readvars.begin(); i!=readvars.end(); ++i) {
    err = dataio_silo::read_variable2grid(*db_ptr, meshname, (*i), SimPM.Ncell);
    if (err)
      rep.error("dataio_silo::ReadData() error reading variable",(*i));
  }
  DBSetDir(*db_ptr,"/");

  DBClose(*db_ptr); //*db_ptr=0; 

  // Now assign Ph to be equal to P for each cell.
  //cell *cpt = gp->FirstPt();
  //do {for(int v=0;v<SimPM.nvar;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=gp->NextPt(cpt))!=0);

  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::choose_filename(const string codefile,
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
        class GridBaseClass *grid
        )
{
  // set grid parameters -- EXPLICITLY UNIFORM FIXED GRID
  if (!grid)
    rep.error("dataio_silo::setup_grid_properties() null grid pointer!",grid);
  //double dx=gp->DX();
  double dx=SimPM.dx;
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
  node_coords = mem.myalloc(node_coords,ndim);

  // now setup arrays with locations of nodes in coordinate directions.
  int nn = SimPM.NG[XX]+1; // for N cells, have N+1 nodes.
  nodex = mem.myalloc(nodex,nn);

  for (int i=0;i<nn;i++)
    nodex[i] = SimPM.Xmin[XX]+static_cast<FAKE_DOUBLE>(i)*dx;
  node_coords[0] = nodex;
  nodedims[0] = nn;
  zonedims[0] = nn-1;

  if (ndim>1) {
    nn = SimPM.NG[YY]+1;
    nodey = mem.myalloc(nodey,nn);

    for (int i=0;i<nn;i++)
      nodey[i] = SimPM.Xmin[YY]+static_cast<FAKE_DOUBLE>(i)*dx;
    node_coords[1] = nodey;
    nodedims[1] = nn;
    zonedims[1] = nn-1;
  }
  if (ndim>2) {
    nn = SimPM.NG[ZZ]+1;
    nodez = mem.myalloc(nodez,nn);

    for (int i=0;i<nn;i++)
      nodez[i] = SimPM.Xmin[ZZ]+static_cast<FAKE_DOUBLE>(i)*dx;
    node_coords[2] = nodez;
    nodedims[2] = nn;
    zonedims[2] = nn-1;
  }

  int nopts=4; int csys;
  dataio_silo::GridOpts = DBMakeOptlist(nopts);
  if      (SimPM.coord_sys==COORD_CRT) csys=DB_CARTESIAN;
  else if (SimPM.coord_sys==COORD_CYL) csys=DB_CYLINDRICAL;
  else if (SimPM.coord_sys==COORD_SPH) csys=DB_SPHERICAL;
  else rep.error("bad coord system",SimPM.coord_sys);
  DBAddOption(GridOpts,DBOPT_COORDSYS,&csys);
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



int dataio_silo::setup_write_variables()
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
  if (RT!=0 && SimPM.EP.phot_ionisation) {
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
  if (RT) {
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
  } // if RT
#endif // RT_TESTING_OUTPUTCOL

  //
  // if there are any tracer variables, get their names from SimPM.trtype,
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
      if (static_cast<int>(SimPM.trtype.size()) > 6*(i+1)) {
        temp<<"_"<< SimPM.trtype.substr(6*(i+1),6);
      }
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




int dataio_silo::generate_quadmesh(DBfile *dbfile, string meshname)
{
  DBClearOption(GridOpts,DBOPT_DTIME);
  DBAddOption(GridOpts,DBOPT_DTIME,&SimPM.simtime);
  DBClearOption(GridOpts,DBOPT_CYCLE);
  DBAddOption(GridOpts,DBOPT_CYCLE,&SimPM.timestep);
 
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
  // DBPutQuadmesh requires the data to be (float **), even though it will allow
  // doubles to be written and will store them correctly, so I have to reinterpret/cast
  // the data to float ** before writing it.
  //
  //  float **nc = reinterpret_cast<float **>(node_coords);
  void **nc = reinterpret_cast<void **>(node_coords);
  int err = DBPutQuadmesh(dbfile, meshname.c_str(), coordnames, nc, nodedims, ndim, FAKE_DATATYPE, DB_COLLINEAR, GridOpts);
  //char temp[256]; strcpy(temp,meshname.c_str());
  //int err = DBPutQuadmesh(dbfile, temp, coordnames, nc, nodedims, ndim, FAKE_DATATYPE, DB_COLLINEAR, GridOpts);
  

  for (int i=0;i<ndim;i++) coordnames[i] = mem.myfree(coordnames[i]);
  coordnames = mem.myfree(coordnames);

  return err;
}



// ##################################################################
// ##################################################################



void dataio_silo::create_data_arrays()
{
  //
  // first check if we have the data arrays set up yet.
  // We need at least one array for a scalar, and two more for a vector.
  //
  if (!data0) {
    data0 = mem.myalloc(data0, SimPM.Ncell);
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
    data1 = mem.myalloc(data1, SimPM.Ncell);
  }
  if ((vec_length>2) && (!data2)) {
    data2 = mem.myalloc(data2, SimPM.Ncell);
  }

  if ((vec_length>1) && (!vec_data)) {
    vec_data = mem.myalloc(vec_data, vec_length);
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
  data0    = mem.myfree(data0   );
  data1    = mem.myfree(data1   );
  data2    = mem.myfree(data2   );
  vec_data = mem.myfree(vec_data);
  return;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_variable2mesh(DBfile *dbfile,  ///< pointer to silo file.
				     string meshname, ///< name of mesh to write to.
				     string variable  ///< variable name to write.
				     )
{
  if (!data0) rep.error("allocate data arrays before trying to write data!",data0);
  int err=0;
  
  if (variable=="Velocity" || variable=="MagneticField") {
    // put data into vec array.
    err = get_vector_data_array(variable, vec_data);
    if (err) rep.error("failed to get vector data for var.",variable);
    // write data to mesh.
    err = write_vector2mesh(dbfile, meshname, variable, vec_data);
    if (err) rep.error("failed to write vector data for var.",variable);
  } // vector variable

  else {
    // scalar variable, already have array, so just get data and write it.
    err = get_scalar_data_array(variable, data0);
    if (err) rep.error("failed to get scalar data for var.",variable);
    err = write_scalar2mesh(dbfile, meshname, variable, data0);
    if (err) rep.error("failed to write scalar data for var.",variable);
  } // scalar variable
  
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::get_scalar_data_array(string variable, ///< variable name to get.
				       FAKE_DOUBLE *data     ///< array to write to.
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

  // Now pick out the data requested cell by cell, and put it into
  // the 1D array.
  cell *c = gp->FirstPt(); long int ct=0;
  if (v>=0) {
    //    cout <<"writing variable v="<<v<<" corresponding to "<<variable<<"\n";
    do {
      data[ct] = c->P[v];
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
  }

  else if (v==-1) { // internal energy (or temperature if we have microphysics)
    //    cout <<"writing variable v="<<v<<" corresponding to "<<variable<<"\n";
    if (MP) {
      do {
	data[ct] = MP->Temperature(c->P,SimPM.gamma);
	ct++;
	//cout <<"temp="<<data[ct-1]<<"\n";
      } while ( (c=gp->NextPt(c))!=0 );
    }
    else {
      do {
	data[ct] = eqn->eint(c->P,SimPM.gamma);
	//      cout <<"data ["<<ct<<"] = "<<data[ct] <<"\n";
	ct++;
      } while ( (c=gp->NextPt(c))!=0 );
    }
  }

  else if (v==-2) { // divB
    //    cout <<"writing variable v="<<v<<" corresponding to "<<variable<<"\n";
    int vars[3];
    vars[0] = static_cast<int>(BX);
    vars[1] = static_cast<int>(BY);
    vars[2] = static_cast<int>(BZ);
    do {data[ct] = eqn->Divergence(c,0,vars, gp); ct++;}
    while ( (c=gp->NextPt(c))!=0 );
  }

  else if (v==-5) { // CurlB (for 2D data only!)
    //    cout <<"writing variable v="<<v<<" corresponding to "<<variable<<"\n";
    int vars[3];
    vars[0] = static_cast<int>(BX);
    vars[1] = static_cast<int>(BY);
    vars[2] = static_cast<int>(BZ);
    double crl[3]; for (int el=0;el<3;el++) crl[el]=0.0;
    do {eqn->Curl(c,0,vars, gp, crl); data[ct] = crl[2]; ct++;}
    while ( (c=gp->NextPt(c))!=0 );
  }

  else if (v==-3) { // total pressure.
    //    cout <<"writing variable v="<<v<<" corresponding to "<<variable<<"\n";
    do {data[ct] = eqn->Ptot(c->P,SimPM.gamma); ct++;}
    while ( (c=gp->NextPt(c))!=0 );
  }

#ifdef RT_TESTING_OUTPUTCOL
  else if (v<=-20 && v>-30) {
    // ionising-RT column density variable.
#ifdef TESTING
    cout <<"writing variable "<<v<<" corresponding to NH0 RT variable "<<variable<<"\n";
#endif
    double Tau[MAX_TAU];
    int col_id=abs(v+20);
    // which Tau variable?  get from string.
    int iT = atoi(variable.substr(11).c_str());
    do {
      CI.get_col(c, col_id, Tau);
      data[ct] = Tau[iT];
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v<=-10 && v>-20) {
    // diffuse-RT column density variable.
#ifdef TESTING
    cout <<"writing variable "<<v<<" corresponding to Ntot RT variable "<<variable<<"\n";
#endif
    double Tau[MAX_TAU];
    int col_id=abs(v+10);
    // which Tau variable?  get from string.
    int iT = atoi(variable.substr(11).c_str());
    do {
      CI.get_col(c, col_id, Tau);
      data[ct] = Tau[iT];
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
  }
#endif // RT_TESTING_OUTPUTCOL

#ifdef COUNT_ENERGETICS
  else if (v==-105) {
    do {data[ct] = c->e.ci_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-106) {
    do {data[ct] = c->e.ci_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-107) {
    do {data[ct] = c->e.pi_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-108) {
    do {data[ct] = c->e.pi_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-109) {
    do {data[ct] = c->e.rr_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-110) {
    do {data[ct] = c->e.rr_rate; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-111) {
    do {data[ct] = c->e.fn_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-112) {
    do {data[ct] = c->e.tot_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-113) {
    do {data[ct] = c->e.tot_cooling; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-114) {
    do {data[ct] = c->e.net_heating; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-115) {
    do {data[ct] = c->e.cooling_time; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
  else if (v==-116) {
    do {data[ct] = c->e.recomb_time; ct++; 
    } while ( (c=gp->NextPt(c))!=0 );
  }
    
#endif // COUNT_ENERGETICS

  else rep.error("Don't understand what variable to write.",v);

  // This won't work for parallel grid... b/c need mpiPM.LocalNcell...
  //  if (ct!=SimPM.Ncell) rep.error("Counting cells error",ct-SimPM.Ncell); 
  return 0;
}



// ##################################################################
// ##################################################################



int dataio_silo::get_vector_data_array(string variable, ///< variable name to get.
				       FAKE_DOUBLE **buffer  ///< array to write to.
				       )
{
  int err=0;
  if (variable=="Velocity") {
    err += get_scalar_data_array("VelocityX",buffer[0]);
    if (vec_length>1) err += get_scalar_data_array("VelocityY",buffer[1]);
    if (vec_length>2) err += get_scalar_data_array("VelocityZ",buffer[2]);
  }
  else if (variable=="MagneticField") {
    err += get_scalar_data_array("MagneticFieldX",buffer[0]);
    if (vec_length>1) err += get_scalar_data_array("MagneticFieldY",buffer[1]);
    if (vec_length>2) err += get_scalar_data_array("MagneticFieldZ",buffer[2]);
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
				   FAKE_DOUBLE *data     ///< pointer to data array.
				   )
{
  //cout <<"writing variable "<<variable<<" to mesh.\n";
  //
  // data has to be passed to the function as float *, even though it will
  // store doubles correctly, so we reinterpret/cast the data first.
  //
  //  float *dd = reinterpret_cast<float *>(data);
  //  if (variable=="CurlB") cout <<"Writing curlB data!\n";
  void *dd = reinterpret_cast<void *>(data);
  int err = DBPutQuadvar1(dbfile, variable.c_str(), meshname.c_str(), dd, zonedims, ndim, 0,0, FAKE_DATATYPE, DB_ZONECENT, 0);

  //char var_str[strlength];  strcpy(var_str,variable.c_str());
  //char mesh_str[strlength]; strcpy(mesh_str,meshname.c_str());
  //int err = DBPutQuadvar1(dbfile, var_str, mesh_str, dd, zonedims, ndim, 0,0, FAKE_DATATYPE, DB_ZONECENT, 0);
  //  if (err) rep.error("couldn't write density to quadmesh",err);
  return err;
}



// ##################################################################
// ##################################################################



int dataio_silo::write_vector2mesh(DBfile *dbfile,  ///< silo file pointer.
				   string meshname, ///< mesh name
				   string variable, ///< variable name
				   FAKE_DOUBLE **data    ///< pointer to data array.
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
  // data has to be passed to the function as float **, even though it will
  // store doubles correctly, so we reinterpret/cast the data first.
  //
  //float **dd = reinterpret_cast<float **>(data);
  void **dd = reinterpret_cast<void **>(data);
  err = DBPutQuadvar(dbfile, variable.c_str(), meshname.c_str(), vec_length, vnames, dd, zonedims, ndim, 0,0, FAKE_DATATYPE, DB_ZONECENT, 0);

  //char var_str[strlength];  strcpy(var_str,variable.c_str());
  //char mesh_str[strlength]; strcpy(mesh_str,meshname.c_str());
  //err = DBPutQuadvar(dbfile, var_str, mesh_str, vec_length, vnames, dd, zonedims, ndim, 0,0, FAKE_DATATYPE, DB_ZONECENT, 0);
  //  if (err) rep.error("couldn't write velocity to quadmesh",err);

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



int dataio_silo::set_readvars(int eqns ///< equations we are solving.
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
    if (eqns==EQGLM) 
      readvars.push_back("glmPSI");
  }
  //
  // if there are any tracer variables, get their names from SimPM.trtype,
  // if it contains the details.  At least it's name will start with 
  // "TrXXX", which is all we really need.
  //
  if (SimPM.ntracer>0) {
    string s; ostringstream temp;
    for (int i=0; i<SimPM.ntracer; i++) {
      s.erase(); temp.str("");
      temp<< "Tr";
      temp.width(3); temp.fill('0'); temp << i;
      if (static_cast<int>(SimPM.trtype.size()) > 6*(i+1)) {
        temp<<"_"<< SimPM.trtype.substr(6*(i+1),6);
      }
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




int dataio_silo::read_variable2grid(DBfile *dbfile,  ///< pointer to silo file.
				    string, ///< name of mesh to read from (can use it for debugging)
				    string variable, ///< variable name to read.
				    long int npt     ///< number of points we are expecting.
				    )
{
  DBquadvar *silodata=0;
  silodata = DBGetQuadvar(dbfile,variable.c_str());
  if (!silodata)
    rep.error("dataio_silo::read_variable2grid() failed to read variable",variable);
  if (silodata->nels != npt)
    rep.error("dataio_silo::read_variable2grid() wrong number of cells",silodata->nels-SimPM.Ncell);

  //
  // Create a pointer to the data in the silo stuct DBquadvar.  This is hardcoded
  // in silo to be floats, but doubles can be written to it, so I have to 
  // re-interpret/cast it in order to read doubles correctly.
  //
  FAKE_DOUBLE **data = reinterpret_cast<FAKE_DOUBLE **>(silodata->vals);

  if (variable=="Velocity" || variable=="MagneticField") {
    int v1,v2,v3;
    if (variable=="Velocity") {v1=VX;v2=VY;v3=VZ;}
    else                      {v1=BX;v2=BY;v3=BZ;}
    //    cout <<"name: "<<silodata->name<<"\tnels="<<silodata->nels<<"\n";
    //    cout <<"ndims: "<<silodata->ndims<<"\tnvals: "<<silodata->nvals<<"\n";
    //cout <<"reading variable "<<variable<<" into element "<<v1<<" of state vec.\n";
    cell *c=gp->FirstPt(); long int ct=0;
    do {
      //      cout <<"ct="<<ct<<"\t and ncell="<<npt<<"\n";
      c->P[v1] = data[0][ct];
      c->P[v2] = data[1][ct];
      c->P[v3] = data[2][ct];
      //c->P[v1] = silodata->vals[0][ct];
      //c->P[v2] = silodata->vals[1][ct];
      //c->P[v3] = silodata->vals[2][ct];
      //cout <<"ct="<<ct<<"\t and ncell="<<npt<<"\n";
      ct++;
    } while ( (c=gp->NextPt(c))!=0 );
    if (ct != npt) rep.error("wrong number of points read for vector variable",ct-npt);
  } // vector variable

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
    //else if (variable.substr(0,3)=="Tr0") {v1=SimPM.ftr;}
    //else if (variable.substr(0,3)=="Tr1") {v1=SimPM.ftr+1;}
    //else if (variable.substr(0,3)=="Tr2") {v1=SimPM.ftr+2;}
    //else if (variable.substr(0,3)=="Tr3") {v1=SimPM.ftr+3;}
    //else if (variable.substr(0,3)=="Tr4") {v1=SimPM.ftr+4;}
    //else if (variable.substr(0,3)=="Tr5") {v1=SimPM.ftr+5;}
    //else if (variable.substr(0,3)=="Tr6") {v1=SimPM.ftr+6;}
    //else if (variable.substr(0,3)=="Tr7") {v1=SimPM.ftr+7;}
    //else if (variable.substr(0,3)=="Tr8") {v1=SimPM.ftr+8;}
    //else if (variable.substr(0,3)=="Tr9") {v1=SimPM.ftr+9;}
    else rep.error("what var to read???",variable);

    //cout <<"reading variable "<<variable<<" into element "<<v1<<" of state vec.\n";
    cell *c=gp->FirstPt();
    long int ct=0;
    do {
      //cout <<"val="<<silodata->vals[0][ct]<<" and data="<<data[0][ct]<<"\n";
      //c->P[v1] = silodata->vals[0][ct];
      c->P[v1] = data[0][ct];
      ct++;
      //      cout <<"ct="<<ct<<"\t and ncell="<<npt<<"\n";
    } while ( (c=gp->NextPt(c))!=0 );
    if (ct != npt) rep.error("wrong number of points read for scalar variable",ct-npt);
  } // scalar variable

  //  cout <<"Read variable "<<variable<<"\n";
  DBFreeQuadvar(silodata); //silodata=0;
  data=0;
  return 0;
}



// ##################################################################
// ##################################################################




#endif // if SILO
