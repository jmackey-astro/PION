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
/// - 2010-07-20/21 JM: Work on new dataio structure: replaced dbfile with dp_ptr
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
/// - 2011.03.21 JM: Updated column-density variables for new cell interface functions.
/// - 2011.03.22 JM: Removed setup_write_variables() function -- now use serial version.
/// - 2011.10.24 JM: wrapped most of the info-reporting with ifdef-testing.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.01.28 JM: Changes for new code structure.
///
#ifdef SILO

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#ifndef PARALLEL
#error "PARALLEL not defined!  don't compile dataio_silo_MPI.cc without it!"
#endif
#include "dataio_silo.h"
#include <cstring>
#include <sstream>

dataio_silo_pllel::dataio_silo_pllel(
      class MCMDcontrol *p
      )
{
#ifdef TESTING
  cout <<"Setting up parallel Silo I/O class.\n";
#endif
  numfiles=-1;
  dataio_silo_pllel::mpiPM = p;
  return;
}

dataio_silo_pllel::~dataio_silo_pllel()
{
#ifdef TESTING
  cout <<"Deleting parallel Silo I/O class.\n";
#endif
  return;
}

int dataio_silo_pllel::ReadHeader(string infile ///< file to read from
				  )
{
  int err=0;
  silofile=infile;
#ifdef TESTING
  cout <<"Rank: "<<mpiPM->myrank<<"\tReading Header from file: "<<silofile<<"\n";
#endif
  if (!file_exists(silofile))
    rep.error("dataio_silo_pllel::ReadHeader() File not found, myrank follows",mpiPM->myrank);

  int group_rank=0, myrank_group=0;
  string file_id="read_header";
  int num_files=1;
  err = COMM->silo_pllel_init(num_files,"READ", file_id, &group_rank, &myrank_group);
  if (err) rep.error("COMM->silo_pllel_init() returned err",err);

  //
  // Now wait for baton, and open into /header directory in file.
  //
  *db_ptr=0;
  string mydir = "/header";
  err = COMM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr)) rep.error("COMM->silo_pllel_wait_for_file() returned err",err);

  //
  // Now read the header, and also NUM_FILES, which tells me how many files there are for
  // when I need to read in the data later on.
  //
  //  err = dataio_silo::read_header(*db_ptr);
  err = read_simulation_parameters();
  dataio_silo::ndim = SimPM.ndim;
  if (err)
    rep.error("dataio_silo::ReadHeader() error reading header from silo file",err);
  err += DBReadVar(*db_ptr,"NUM_FILES",   &numfiles);
  if (err) {
    numfiles=1; err =0;
#ifdef TESTING
    cout <<"Warning didn't read NUM_FILES from silo file.\n";
#endif
    //    rep.error("dataio_silo::ReadHeader() error reading NUM_FILES from silo file",err);
  }

  //
  // Finished Local work; hand off baton and wait!
  //
  err = COMM->silo_pllel_finish_with_file(file_id,db_ptr);
  if (err) rep.error("COMM->silo_pllel_finish_with_file() returned err",err);

#ifdef TESTING
  cout <<"Rank: "<<mpiPM->myrank<<"\tFINISHED reading Header from file: "<<silofile<<"\n";
#endif
  return err;
}


int dataio_silo_pllel::ReadData(
        string infile,
	class GridBaseClass *cg
	)
{
  if (!cg)
    rep.error("dataio_silo_pllel::ReadData() null pointer to grid!",cg);
  dataio_silo::gp = cg;

  int err=0;
#ifdef TESTING
  cout <<"\n----Rank: "<<mpiPM->myrank<<"\tReading Data from files: "<<infile;
#endif

  // set grid properties for quadmesh 
  if (!have_setup_gridinfo) {
    err = setup_grid_properties(gp);
    rep.errorTest("dataio_silo_pllel::ReadData() error setting up grid_props", 0, err);
  }

  // check numfiles was read from header
  if (numfiles<0)
    rep.error("Failed to get numfiles from header! but error not detected!",numfiles);
  //  cout <<"\tnumfiles="<<numfiles<<" and nproc="<<mpiPM->nproc<<"\n";


  //
  // Get ready to open the file:
  //
  int group_rank=0, myrank_ingroup=0;
  string file_id="read_data";
  err = COMM->silo_pllel_init(numfiles,"READ", file_id, &group_rank, &myrank_ingroup);
  if (err) rep.error("COMM->silo_pllel_init() returned err",err);

  // get group rank and index from hard-coding...
  //int group_rank = mpiPM->myrank*numfiles/mpiPM->nproc;
  //cout <<"myrank: "<<mpiPM->myrank<<"\tnumfiles: "<<numfiles<<"\tmy_group: "<<group_rank<<"\n";
  //int myrank_ingroup = mpiPM->myrank % (mpiPM->nproc/numfiles);
  //cout <<"myrank: "<<mpiPM->myrank<<"\tmy_index_in_group: "<<myrank_ingroup<<"\n";

  // Choose correct filename based on numfiles.
  // Also choose correct directory name, assuming that nproc is the same as for
  // when the file was generated.
  // *** I would like to relax this requirement eventually, but it would take a lot of
  // *** coding for the two cases of reading more than one mesh to the local domain, and reading
  // *** parts of a larger mesh to multiple local domains.

  //
  // replace 0000 in filename with my group rank, based on numfiles.
  //
  ostringstream temp; 
  temp.fill('0'); temp.str("");
  string::size_type pos = infile.find("0000");
  if (pos==string::npos) {
    cout <<"didn't fine 0000 in file, so assuming reading a single file.\n";
    //rep.error("not a parallel i/o filename!",infile);
  }
  else {
    temp.str("");
    temp.width(4);
    temp<<group_rank;
    infile.replace(pos,4,temp.str());
#ifdef TESTING
    cout <<"\tNew infile: "<<infile<<"\n";
#endif
    temp.str("");
  }
  silofile.clear();
  silofile=infile;

  //
  // choose my directory name:
  //
  string mydir; temp.str("");
  temp << "rank_"; temp.width(4); temp << mpiPM->myrank << "_domain_";
  temp.width(4); temp<<myrank_ingroup; mydir = temp.str(); temp.str("");
  
  *db_ptr=0;
  err = COMM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr)) rep.error("COMM->silo_pllel_wait_for_file() returned err",err);

  if (silo_filetype==DB_HDF5) {
    char *compress = DBGetCompression();
    //    cout <<"compression="<<compress<<"\n";
    DBSetCompression(compress);
    int friendly=DBGetFriendlyHDF5Names();
    DBSetFriendlyHDF5Names(friendly);
  }

  //
  // Check extents of quadmesh are what i have determined my domain to be...
  //  err = check_subdomain_extents(db_ptr);
  //
  temp.str(""); temp <<"unigrid";
  temp.width(4);
  temp<<mpiPM->myrank;
  string meshname=temp.str();
  DBquadmesh *qm = DBGetQuadmesh(*db_ptr,meshname.c_str());
  if (!qm) rep.error("failed to find quadmesh named as follows",meshname);
  int n; 
  //
  // first check dimensions.
  //
  if ((n=qm->ndims)!=ndim) rep.error("bad dimensionality in qmesh",n);
  //
  // check number of cells.
  //
  int nn=1; for (int i=0;i<ndim;i++) nn*=mpiPM->LocalNG[i]+1;
  if ((n=qm->nnodes)!=nn) {
    cerr<<"nnodes="<<n<<" and ncell="<<mpiPM->LocalNcell<<"\n";
    rep.error("bad number of nodes",n-mpiPM->LocalNcell);
  }
  //
  // Have to reinterpret qm->coords, as they are doubles, but the interface thinks
  // they are floats.
  //
  FAKE_DOUBLE **c = reinterpret_cast<FAKE_DOUBLE **>(qm->coords);
  //
  // check origin of subgrid is at the right place.
  //
  if (!GS.equalD(nodex[0],c[XX][0])) rep.error("XX mesh not at right place...",nodex[0]-c[XX][0]);
  if (ndim>1)
    if (!GS.equalD(nodey[0], c[YY][0])) {
      cout <<"nodey = "<<nodey[0]<<" and coords[y]= "<<c[YY][0]<<"\n";
      rep.error("YY mesh not at right place...",nodey[0]-c[YY][0]);
    }
  if (ndim>2)
    if (!GS.equalD(nodez[0], c[ZZ][0]))
      rep.error("ZZ mesh not at right place...",nodez[0]-c[ZZ][0]);
  //
  // this is probably enough checking! -- if all these are right then it should all be right.
  //
  DBFreeQuadmesh(qm); //qm=0;

  //
  // now read each variable in turn from the mesh
  //
  err = set_readvars(SimPM.eqntype);
  if (err) rep.error("failed to set readvars in ReadData",err);
  for (std::vector<string>::iterator i=readvars.begin(); i!=readvars.end(); ++i) {
    err = read_variable2grid(*db_ptr, meshname, (*i), mpiPM->LocalNcell);
    if (err)
      rep.error("dataio_silo_pllel::ReadData() error reading variable",(*i));
  }
  // Now assign Ph to be equal to P for each cell.
  cell *cpt = gp->FirstPt();
  do {
    for(int v=0;v<SimPM.nvar;v++) cpt->Ph[v]=cpt->P[v];
  } while ((cpt=gp->NextPt(cpt))!=0);

  //
  // Finished Local work; hand off baton and wait!
  //
  err = COMM->silo_pllel_finish_with_file(file_id,db_ptr);
  if (err) rep.error("COMM->silo_pllel_finish_with_file() returned err",err);
  COMM->barrier("dataio_silo_pllel__ReadData");

#ifdef TESTING
  cout <<"----Rank: "<<mpiPM->myrank<<"\tFINISHED reading Data from file: "<<silofile<<"\n"<<"\n";
#endif
  return err;
}


int dataio_silo_pllel::OutputData(const string outfilebase,
				  class GridBaseClass *cg,
				  const long int file_counter   ///< number to stamp file with (e.g. timestep)
				  )
{
  int err=0;
  if (!cg)
    rep.error("dataio_silo_pllel::OutputData() null pointer to grid!",cg);
  dataio_silo::gp = cg;

  //  COMM->barrier("dataio_silo_pllel::OutputData() starting");
#ifdef TESTING
  cout <<"\n----dataio_silo_pllel::OutputData() Writing data to filebase: "<<outfilebase<<"\n";
#endif
  // First initialise the I/O.
  if (silo_filetype==DB_HDF5) {
    DBSetCompression("METHOD=GZIP LEVEL=1");
    DBSetFriendlyHDF5Names(1);
#ifdef TESTING
    cout <<" *** setting compression.\n";
#endif
  }

  ///
  /// \todo not too sure what to do about numfiles -- I have it hardcoded for
  /// now, but could set it as a runtime parameter at some stage. 
  /// I have set it to produce files of up to about 1.2GB and then split.
  ///
  int threshold=static_cast<int>(SimPM.Ncell/1.0e7)+1;
  if (threshold<1) threshold=1;
  //
  // Have to make sure that we write at most one file per process!
  //
  if (mpiPM->nproc>threshold)
    numfiles=threshold;
  else
    numfiles=mpiPM->nproc;

  //  cout <<"----dataio_silo_pllel::OutputData() running pmpio_init\n";
  int group_rank=0, myrank_ingroup=0;
  string file_id="write_data";
  err = COMM->silo_pllel_init(numfiles,"WRITE", file_id, &group_rank, &myrank_ingroup);
  if (err) rep.error("COMM->silo_pllel_init() returned err",err);

#ifdef TESTING
  cout <<"myrank: "<<mpiPM->myrank<<"\tnumfiles: "<<numfiles<<"\tmy_group: ";
  cout <<group_rank <<"\tmy_index_in_group: "<<myrank_ingroup<<"\n";
#endif // TESTING

  //
  // Choose output filename:
  //
  //  cout <<"setting strings... outfilebase="<<outfilebase<<"\n";
  choose_pllel_filename(outfilebase,group_rank,file_counter,silofile);
  //  cout <<"string for outfile set...\n";

  //
  // Choose directory within silo file.
  //
  string mydir;
  ostringstream temp; temp.fill('0'); temp.str("");
  temp << "rank_";
  temp.width(4);
  temp << mpiPM->myrank << "_domain_";
  temp.width(4);
  temp << myrank_ingroup;
  mydir = temp.str();
  temp.str("");
  if (mydir.size()+1>=strlength) rep.error("string too large",mydir);
  //  cout <<"string for my directory in file set...\n";


  //
  // Wait for my turn to write to the file.
  //
  //  cout <<"----dataio_silo_pllel::OutputData() waiting for baton\n";
  *db_ptr=0;
  err = COMM->silo_pllel_wait_for_file(file_id, silofile, mydir, db_ptr);
  if (err || !(*db_ptr)) rep.error("COMM->silo_pllel_wait_for_file() returned err",err);

  //
  // Have got the baton, now, so the file is mine to write to.
  // local work here... each proc write their part of the grid.
  //
  err=0;
  if (!have_setup_gridinfo) {
    // set grid properties for quadmesh 
    // CHANGED FOR MPI LOCAL GRID
    //cout <<"----dataio_silo_pllel::OutputData() setting up grid properties\n";
    err = setup_grid_properties(gp);
    if (err)
      rep.error("dataio_silo_pllel::OutputData() error setting up grid_props", err);
    //cout <<"----dataio_silo_pllel::OutputData() grid props setup done\n";
  }
  if (!have_setup_writevars) {
    // set what data to write to the mesh.
    // NOT THE SAME AS SERIAL VERSION -- ONLY PRIMITIVE VARS!!!
    //cout <<"----dataio_silo_pllel::OutputData() setting up write variables\n";
    err = setup_write_variables();
    if (err)
      rep.error("dataio_silo_pllel::OutputData() error settting up variables to write to",err);
    //cout <<"----dataio_silo_pllel::OutputData() write vars setup done\n";
  }


  //
  // Generate Mesh in file.
  //
  //  cout <<"----dataio_silo_pllel::OutputData() generating quadmesh\n";
  temp.str(""); temp <<"unigrid";
  temp.width(4);
  temp<<mpiPM->myrank;
  string meshname=temp.str();
  err = generate_quadmesh(*db_ptr, meshname);
  if (err)
    rep.error("dataio_silo_pllel::OutputData() error writing quadmesh to silo file",err);
  //  cout <<"----dataio_silo_pllel::OutputData() quadmesh generated\n";

  //
  // now write each variable in turn to the mesh
  //
  //  cout <<"----dataio_silo_pllel::OutputData() creating data arrays\n";
  create_data_arrays();
  //  cout <<"----dataio_silo_pllel::OutputData() arrays created\n";
  for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
    //    cout <<"\twriting variable "<<(*i)<<" to file "<<silofile<<"\n";
    err = write_variable2mesh(*db_ptr, meshname, (*i));
    //    cout <<"\t\tvariable "<<(*i)<<" written.\n";
    if (err) rep.error("dataio_silo_pllel::OutputData() error writing variable",(*i));
  }
  //  cout <<"----dataio_silo_pllel::OutputData() deleting arrays\n";
  delete_data_arrays();
  //  cout <<"----dataio_silo_pllel::OutputData() arrays deleted\n";

  // ONLY DO THIS IF I AM ROOT PROCESSOR *IN EACH GROUP*
  // GET CURRENT DIR, MOVE OUT, MAKE HEADER,
  // MOVE BACK TO CURRENT DIR
  if (myrank_ingroup==0) {
    //    cout <<"myrank_ingroup="<<myrank_ingroup<<" so writing header to file.\n";
    // Am keeping file open, so no need to re-open it anymore.
    //    *db_ptr =0;
    //    *db_ptr = DBOpen(fname,DB_UNKNOWN,DB_APPEND);
    
    char datadir[strlength];
    err += DBGetDir(*db_ptr,datadir);
    DBSetDir(*db_ptr,"/");

    // WRITE HEADER
    DBMkDir(*db_ptr,"header");
    DBSetDir(*db_ptr,"/header");
    //
    // Write Number of files to header directory of file.
    //
    int dim1[1]; dim1[0]=1;
    err += DBWrite(*db_ptr,"NUM_FILES",   &numfiles,    dim1,1,DB_INT);
    err += DBWrite(*db_ptr,"MPI_nproc",   &mpiPM->nproc, dim1,1,DB_INT);
    //    err = write_header(*db_ptr);
    err = write_simulation_parameters();
    if (err)
      rep.error("dataio_silo_pllel::OutputData() error writing header to silo file",err);

    DBSetDir(*db_ptr,"/");
    DBSetDir(*db_ptr,datadir);
    //DBClose(*db_ptr); *db_ptr=0;
    //    cout <<"Finished writing header to file.\n";
  } // if root processor of *group*

  if (mpiPM->myrank==0) {
    //
    // WRITE MULTIMESH DATA NOW TO FIRST PARTIAL FILE IF ROOT PROC.
    //
    string mm_name ="MultiMesh";
    string mma_name="Domain_Decomposition";
    string mrgt_name=""; //MRGtree";
    
    //    cout <<"Writing multimesh data to file...\n";
    //err = write_multimesh(*db_ptr, gp, mm_name, mma_name);
    //if (err) rep.error("Failed to write multimesh Object",err);
    char datadir[strlength];
    int err = DBGetDir(*db_ptr,datadir);
    DBSetDir(*db_ptr,"/");
    // Not putting the multimesh in its own directory anymore.
    //    DBMkDir(*db_ptr,"mmesh");
    //    DBSetDir(*db_ptr,"/mmesh");
    
    //
    // Strip the path from outfilebase to get just the filename; do this by 
    // reverse finding '/' from outfilebase and taking the substring following
    // its position.  As far as I can tell this is foolproof...
    //
    string fname;
    string::size_type p = outfilebase.rfind("/");
    if (p != string::npos) {
      fname = outfilebase.substr(p+1);
    }
    else fname = outfilebase;
    //    cout <<"outfilebase="<<outfilebase<<" and fname="<<fname<<"\n";
    
    
    //  string mm_name="mesh";
    int nmesh = mpiPM->nproc;
    char **mm_names=0;
    int *meshtypes=0;
    int 
      *groups=0, // lists which group each process is in.
      *ranks=0;  // lists which rank each process has in its group.

    meshtypes = mem.myalloc(meshtypes, nmesh);
    mm_names = mem.myalloc(mm_names, nmesh);
    for (int i=0;i<nmesh;i++)
      mm_names[i] = mem.myalloc(mm_names[i], strlength);
    groups = mem.myalloc(groups, nmesh);
    ranks  = mem.myalloc(ranks,  nmesh);
    
    //
    // Assign groups and rank_in_groups to each process:
    //
    for (int i=0;i<nmesh;i++) {
      COMM->silo_pllel_get_ranks(file_id, i, &(groups[i]), &(ranks[i]));
    }
    
    int csys=-1;
    //csys = DB_QUADMESH;
    if      (SimPM.coord_sys==COORD_CRT) csys=DB_QUAD_RECT;
    else if (SimPM.coord_sys==COORD_CYL) csys=DB_QUAD_RECT;
    else if (SimPM.coord_sys==COORD_SPH) csys=DB_QUAD_RECT;
    else rep.error("bad coord system",SimPM.coord_sys);
    for (int i=0;i<nmesh;i++) {
      meshtypes[i] = csys;
      temp.str("");
      if (groups[i] != group_rank) {
	//
	// the root processor's group file doesn't need to be in the path,
	// so only other groups' files have to get listed.
	//
	//temp << fname <<"_"; temp.width(4); temp <<groups[i]<<".";
	//temp.width(6); temp <<SimPM.timestep<<".silo:";
	string temp_filename;
	choose_pllel_filename(fname, groups[i], file_counter, temp_filename);
	temp << temp_filename.c_str() <<":";
      }
      temp << "/rank_"; temp.width(4); temp << i << "_domain_";
      temp.width(4); temp<<ranks[i];
      temp <<"/unigrid"; temp.width(4); temp<<i;
      if (temp.str().length() >strlength-1)
	rep.error("directory name too long!!!",temp.str().length());
      strcpy(mm_names[i],temp.str().c_str());
    }

    DBoptlist *mm_opts = DBMakeOptlist(7);
    DBAddOption(mm_opts,DBOPT_DTIME,&SimPM.simtime);
    DBAddOption(mm_opts,DBOPT_CYCLE,&SimPM.timestep);
    //  DBAddOption(mm_opts,DBOPT_ADJACENCY_NAME,mma_name.c_str()); // doesn't exist!!!
    //char mrgt[256]; strcpy(mrgt,mrgt_name.c_str());
    //DBAddOption(mm_opts,DBOPT_MRGTREE_NAME,mrgt);
    int blockorigin=0;
    DBAddOption(mm_opts,DBOPT_BLOCKORIGIN,&blockorigin);
    //
    // performance options:
    //
    int ext_size = 2*ndim;
    double extents[ext_size*nmesh];
    int zonecounts[nmesh];
    int externalzones[nmesh];
    class MCMDcontrol pp; pp.nproc = mpiPM->nproc;
    for (int v=0; v<nmesh; v++) {
      pp.myrank = v;
      pp.decomposeDomain();
      zonecounts[v] = pp.LocalNcell;
      externalzones[v] = 0; // set to 1 if domain has zones outside multimesh.
      for (int i=0; i<ndim;i++) {
	//extents[ext_size*v+2*i  ] = pp.LocalXmin[i];
	//extents[ext_size*v+2*i+1] = pp.LocalXmax[i];
	extents[ext_size*v+i     ] = SimPM.Xmin[i];
	extents[ext_size*v+ndim+i] = SimPM.Xmax[i];
      }
    } // loop over meshes
    DBAddOption(mm_opts,DBOPT_EXTENTS_SIZE,&ext_size);
    DBAddOption(mm_opts,DBOPT_EXTENTS,extents);
    DBAddOption(mm_opts,DBOPT_ZONECOUNTS,zonecounts);
    DBAddOption(mm_opts,DBOPT_HAS_EXTERNAL_ZONES,externalzones);


    // write the multimesh
    err = DBPutMultimesh(*db_ptr, mm_name.c_str(), nmesh, mm_names, meshtypes, mm_opts);
    if (err) rep.error("dataio_silo_pllel::OutputData() error writing multimesh info",err);
    DBClearOptlist(mm_opts);
   
    // now write the multivars
    //err = DBPutMultivar(*db_ptr, mv_name, nmvar, mv_names, vartypes, mv_opts);
    // I'm going to reuse all the multimesh vars to the multivar object...
    for (std::vector<string>::iterator i=varnames.begin(); i!=varnames.end(); ++i) {
      // multivar name
      string vname;
      vname.clear(); vname=(*i);
      //
      // multivar names:
      //
      for (int ii=0;ii<nmesh;ii++) {
	meshtypes[ii] = DB_QUADVAR;
	temp.str("");
	if (groups[ii] != group_rank) {
	  //
	  // the root processor's group file doesn't need to be in the path,
	  // so only other groups' files have to get listed.
	  //
	  //temp << fname <<"_"; temp.width(4); temp <<groups[ii]<<".";
	  //temp.width(6); temp <<SimPM.timestep<<".silo:";
	  string temp_filename;
	  choose_pllel_filename(fname, groups[ii], file_counter, temp_filename);
	  temp << temp_filename.c_str() <<":";
	}
	temp << "/rank_"; temp.width(4); temp << ii << "_domain_";
	temp.width(4); temp<<ranks[ii]; temp <<"/"<<vname;
	if (temp.str().length() >strlength)
	  rep.error("multivar name too long!!!",temp.str().length());
	strcpy(mm_names[ii],temp.str().c_str());
      }
      
      //
      // Now the performance options for the multivars
      //
      DBAddOption(mm_opts,DBOPT_DTIME,&SimPM.simtime);
      DBAddOption(mm_opts,DBOPT_CYCLE,&SimPM.timestep);
      DBAddOption(mm_opts,DBOPT_BLOCKORIGIN,&blockorigin);
      // ext_size=2;
      // extents[] is the max/min values of the variable -- no time to calculate that.
      char mn[256]; strcpy(mn, mm_name.c_str());
      DBAddOption(mm_opts,DBOPT_MMESH_NAME,mn);

      err = DBPutMultivar(*db_ptr, vname.c_str(), nmesh, mm_names, meshtypes, mm_opts);
      if (err)
	rep.error("dataio_silo_pllel::OutputData() error writing variable",(*i));
      DBClearOptlist(mm_opts);
    }

    //
    // Free memory
    //
    DBFreeOptlist(mm_opts);
    for (int i=0;i<nmesh;i++)
      mm_names[i] = mem.myfree(mm_names[i]);
    mm_names = mem.myfree(mm_names);
    meshtypes = mem.myfree(meshtypes);
    groups = mem.myfree(groups);
    ranks  = mem.myfree(ranks);
    
    DBSetDir(*db_ptr,"/");
    DBSetDir(*db_ptr,datadir);
  
    //    cout <<"Finished writing multimesh data to file...\n";
    
    // 
    // Now write a mulitmeshadj object (but only if nproc>1).
    // 
    if (mpiPM->nproc>1) {
      err = write_multimeshadj(*db_ptr, gp, mm_name, mma_name);
      if (err) rep.error("Failed to write multimesh Adjacency Object",err);
    }
    //
    // Write an MRG Tree object (NOT IMPLEMENTED IN SILO/VISIT YET!!)
    //
    //err = write_MRGtree(*db_ptr, gp, mm_name, mrgt_name);
    //if (err) rep.error("Failed to write MRG tree Object",err);
    
  } // if root processor of whole simulation

  //
  // Finished Local work; hand off baton and wait!
  //
  err = COMM->silo_pllel_finish_with_file(file_id,db_ptr);
  if (err) rep.error("COMM->silo_pllel_finish_with_file() returned err",err);
  COMM->barrier("dataio_silo_pllel__OutputData");

  //   cout <<"Got past barrier... finished outputting silo data.\n\n";
  
#ifdef TESTING
  cout <<"----dataio_silo_pllel::OutputData() Finished writing data to file: "<<silofile<<"\n"<<"\n";
#endif
  return err;
}

int dataio_silo_pllel::choose_pllel_filename(const string fbase, ///< filebase passed in from main code.
					     const int igroup,   ///< group_rank (i.e. which file I write to)
					     const int file_counter,  ///< file counter to use (e.g. timestep).
					     string &outfile     ///< write filename to this string.
					     )
{
  //  long int counter;
  //  if (file_counter<0) counter = 0;
  //  else counter = file_counter;

  ostringstream temp; temp.str(""); temp.fill('0');  
  temp << fbase   <<"_";
  temp.width(4);
  temp << igroup  <<".";
  if (file_counter >= 0) {
    temp.width(Ndigits);
    temp << file_counter <<".";
  }
  temp <<"silo";
  outfile.clear();
  outfile=temp.str();
  if (outfile.size()+1>=strlength) rep.error("string too large",outfile);
  temp.str("");
  return 0;
}

int dataio_silo_pllel::setup_grid_properties(
        class GridBaseClass *grid ///< pointer to data.
        )
{
  // set grid parameters -- EXPLICITLY UNIFORM FIXED GRID
  // This version is for the local domain of the current processor.
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

  nodedims    = mem.myalloc(nodedims,    ndim);
  zonedims    = mem.myalloc(zonedims,    ndim);
  node_coords = mem.myalloc(node_coords, ndim);

  // now setup arrays with locations of nodes in coordinate directions.
  int nn = mpiPM->LocalNG[XX]+1;
  nodex = mem.myalloc(nodex,nn);
  for (int i=0;i<nn;i++)
    nodex[i] = mpiPM->LocalXmin[XX]+static_cast<FAKE_DOUBLE>(i)*dx;
  node_coords[0] = nodex;
  nodedims[0] = nn;
  zonedims[0] = nn-1;

  if (ndim>1) {
    nn = mpiPM->LocalNG[YY]+1;
    nodey = mem.myalloc(nodey,nn);
    for (int i=0;i<nn;i++)
      nodey[i] = mpiPM->LocalXmin[YY]+static_cast<FAKE_DOUBLE>(i)*dx;
    node_coords[1] = nodey;
    nodedims[1] = nn;
    zonedims[1] = nn-1;
  }
  if (ndim>2) {
    nn = mpiPM->LocalNG[ZZ]+1;
    nodez = mem.myalloc(nodez,nn);
    for (int i=0;i<nn;i++)
      nodez[i] = mpiPM->LocalXmin[ZZ]+static_cast<FAKE_DOUBLE>(i)*dx;
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
  //char temp[strlength];
  //strcpy(temp,uc.length.c_str());
  //DBAddOption(GridOpts,DBOPT_XUNITS,temp);

  have_setup_gridinfo = true;
  return 0;
}
 
// int dataio_silo_pllel::write_header(DBfile *dbfile)
// {
//   int err = dataio_silo::write_header(dbfile);
//   if (err) rep.error("writing serial part of header failed",err);

//   // for writing a single variable.
//   int dim1[1]; dim1[0]=1;
//   // for writing a three element vector
//   int dim3[1]; dim3[0]=3;

//   // NOT SURE IF I REALLY NEED THIS!!!  I WANT TO WRITE THE FILES SO
//   // THAT I CAN RESTART WITH A DIFFERENT NUMBER OF PROCESSORS IF
//   // POSSIBLE, AT LEAST IF THE NUMBER DIFFERS BY A FACTOR OF 2.
//   err += DBWrite(dbfile,"MPI_rank",   &mpiPM->myrank,    dim1,1,DB_INT);
//   err += DBWrite(dbfile,"MPI_nproc",  &mpiPM->nproc,     dim1,1,DB_INT);
//   err += DBWrite(dbfile,"MPI_NCell",  &mpiPM->LocalNcell,dim1,1,DB_LONG);
//   err += DBWrite(dbfile,"MPI_NGrid",   mpiPM->LocalNG,   dim3,1,DB_INT);
//   err += DBWrite(dbfile,"MPI_Xmin",    mpiPM->LocalXmin, dim3,1,DB_DOUBLE);
//   err += DBWrite(dbfile,"MPI_Xmax",    mpiPM->LocalXmax, dim3,1,DB_DOUBLE);

//   return err;
// }

void dataio_silo_pllel::create_data_arrays()
{
  // first check if we have the data arrays set up yet.
  // We need at least one array for a scalar, and two more for a vector.
  //cout <<"local Ncell="<<mpiPM->LocalNcell<<" and global Ncell is "<<SimPM.Ncell<<"\n";
  if (!data0) {
    data0 = mem.myalloc(data0, mpiPM->LocalNcell);
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
    data1 = mem.myalloc(data1, mpiPM->LocalNcell);
  }
  if ((vec_length>2) && (!data2)) {
    data2 = mem.myalloc(data2, mpiPM->LocalNcell);
  }

  if ((vec_length>1) && (!vec_data)) {
    vec_data = mem.myalloc(vec_data, vec_length);
    vec_data[0] = data0;
    if (vec_length>1) vec_data[1] = data1;
    if (vec_length>2) vec_data[2] = data2;
  }

  return;
}

//
// Write a mulitmesh adjacency object
// 
int dataio_silo_pllel::write_multimeshadj(DBfile *dbfile, ///< pointer to silo file.
					  class GridBaseClass *ggg, ///< pointer to data.
					  string mm_name, ///< multimesh  name
					  string mma_name ///< multimeshadj name.
					  )
{
  int err=0;
#ifdef TESTING
  cout <<"Writing multimesh adjacency object into Silo file.\n";
#endif

  //
  // Create temporary params class for repeated domain decomposition.
  //
  class MCMDcontrol pp;
  int nmesh = mpiPM->nproc;
  pp.nproc = mpiPM->nproc;
  
  int meshtypes[nmesh], Nngb[nmesh], Sk[nmesh];
  for (int v=0;v<nmesh;v++) {
    //meshtypes[v] = DB_QUAD_RECT;
    meshtypes[v] = DB_QUADMESH;
    Nngb[v] = Sk[v] = 0;
  }
  
  //
  // Now go through all meshes, get number of neighbours in each of
  // them, and also the running sum of the total number of neighbours,
  // because we need some arrays of this length.
  //
  // Note neighbouring domains are not just in the six directions, but
  // also corner-abutting and edge-abutting domains.
  //
  int Stot=0;
  std::vector<int> ngb_list;
  for (int v=0;v<nmesh;v++) {
    pp.myrank = v;
    pp.decomposeDomain();
    //
    // Get list of abutting domains.
    //
    std::vector<int> dl;
    pp.get_abutting_domains(dl);
    Nngb[v] = dl.size();
    for (unsigned int i=0; i<static_cast<unsigned int>(Nngb[v]); i++)
      ngb_list.push_back(dl[i]);
    //for (vector<int>::iterator i=dl.begin(); i!=dl.end(); i++)
    //  ngb_list.push_back(*i);
    
    if (v<nmesh-1)
      Sk[v+1] = Sk[v]+ Nngb[v];
    else
      Stot = Sk[v]+ Nngb[v];
    //cout <<"Sk["<<v<<"]="<<Sk[v]<<"\n";
  }

  //
  // Stot is the total number of neighbours over all meshes.  Sk[v] is
  // the number of meshes in all neighbours up to and including v-1.
  //
  //cout <<"Stot="<<Stot<<"\n";
  if (Stot!=static_cast<int>(ngb_list.size()))
    rep.error("counting error!",Stot-ngb_list.size());

  //
  // Y ngb      = list of neighbour ids. 
  // Y back     = reverse list of neighbour ids.
  // N nzones   = number of zones adjacent to each interface. (UNUSED)
  // N zonelist = index ranges of each adjacent zone. (UNUSED)
  // Y nnodes   = number of elements in nodelist == 15 for structured grid!
  // * nodelist = index ranges of each list of nodes..
  // ct       = counter for where we are in neighbour lists.
  //
  int ngb[Stot], back[Stot], 
    nnodes[Stot], *nodelist[Stot];
  //
  // Each nodelist element has 15 elements regardless of ndim, and we
  // set some to Q=-1 if ndim<3
  //
  int Q      =-1;  // for unused values.
  int nl_len = 15; // size of elements.
  for (int s=0;s<Stot;s++) {
    nnodes[s] = nl_len;
    ngb[s]    = ngb_list[s];
    nodelist[s] = mem.myalloc(nodelist[s], nl_len);
    for (int v=0; v<nl_len; v++) {
      nodelist[s][v] = Q;
    }
  }

  //
  // loop over meshes again and populate the neighbour lists.
  //
  //long int ct=0;
  for (int v=0;v<nmesh;v++) {
    pp.myrank = v;
    pp.decomposeDomain();
    long int off1 = Sk[v]; // this should be the same as ct (maybe don't need ct then!)

    //
    // Assign reverse neighbour's id for each of myrank's neighbours.
    // If myrank=v, ngb[off1+i]=dom2, then back[off1+i] is v's index
    // in dom2's list of neigbours.
    // 
    int dom2;
    long int off2;
    for (int i=0; i<Nngb[v]; i++) {
      dom2 = ngb[off1+i];
      off2 = Sk[dom2];
      if (off2 > Stot) 
	rep.error("Counting error in loop to find back array",off2-Stot);

      for (int s=0; s<Nngb[dom2]; s++) {
	if (ngb[off2+s]==v) back[off1+i] = s;
      } // loop over dom2's neighbours
    } // loop over dom1's neighbours

    //
    // Now the nodelists[]: [0,1] are xmin,xmax [2,3] ymin,ymax [4,5] are zmin,zmax
    // https://email.ornl.gov/pipermail/visit-users/attachments/20081017/a72e7c92/multimesh_cyrus.cc
    //
    // This page says that the ordering should be:
    // local[x0,x1,y0,y1,z0,z1]boundary[x0,x1,y0,y1,z0,z1]orientation[1,2,3]
    // It also requires domains which only share a corner, and it
    // requires nodes to have a global indexing system.
    //
    // We loop over mesh v's neighbours and assign nodelist[Sk[v]+i]
    //
    for (int i=0; i<Nngb[v]; i++) {
      off1 = Sk[v]+i;
      if (off1 > Stot) 
	rep.error("Counting error in loop over neighbours for mesh v",off1-Stot);
      //
      // Local nodelist is the same for all of v's nodelists.
      //
      nodelist[off1][0] = pp.offsets[XX] + 0;
      nodelist[off1][1] = pp.offsets[XX] + pp.LocalNG[XX];
      if (ndim>1) {
	nodelist[off1][2] = pp.offsets[YY] + 0;
	nodelist[off1][3] = pp.offsets[YY] + pp.LocalNG[YY];
      }
      if (ndim>2) {
	nodelist[off1][4] = pp.offsets[ZZ] + 0;
	nodelist[off1][5] = pp.offsets[ZZ] + pp.LocalNG[ZZ];
      }

      //
      // For neighbours nodelists, get their relative position in the
      // block structure.
      //
      int my_ix[MAX_DIM], ngb_ix[MAX_DIM];
      for (int ii=0;ii<MAX_DIM;ii++) my_ix[ii] = ngb_ix[ii] = -1;
      pp.get_domain_ix(pp.myrank,my_ix);
      pp.get_domain_ix(ngb[off1],ngb_ix);
      
      //
      // X-dir first.
      //
      if      ((my_ix[XX]-ngb_ix[XX]) == 1) {
	nodelist[off1][6] = pp.offsets[XX];
	nodelist[off1][7] = pp.offsets[XX];
      }
      else if ((my_ix[XX]-ngb_ix[XX]) ==-1) {
	nodelist[off1][6] = pp.offsets[XX] + pp.LocalNG[XX];
	nodelist[off1][7] = pp.offsets[XX] + pp.LocalNG[XX];
      }
      else if (my_ix[XX]==ngb_ix[XX]) {
	nodelist[off1][6] = pp.offsets[XX];
	nodelist[off1][7] = pp.offsets[XX] + pp.LocalNG[XX];
      }
      else 
	rep.error("domains don't touch (X-dir)!", my_ix[XX]-ngb_ix[XX]);

      //
      // Now Y-dir
      //
      if (ndim>1) {
	if      ((my_ix[YY]-ngb_ix[YY]) == 1) { // neighbour below us
	  nodelist[off1][8] = pp.offsets[YY];
	  nodelist[off1][9] = pp.offsets[YY];
	}
	else if ((my_ix[YY]-ngb_ix[YY]) ==-1) { // neighbour above us
	  nodelist[off1][8] = pp.offsets[YY] + pp.LocalNG[YY];
	  nodelist[off1][9] = pp.offsets[YY] + pp.LocalNG[YY];
	}
	else if (my_ix[YY]==ngb_ix[YY])       { // neighbour level with us.
	  nodelist[off1][8] = pp.offsets[YY];
	  nodelist[off1][9] = pp.offsets[YY] + pp.LocalNG[YY];
	}
	else 
	  rep.error("domains don't touch! (Y-dir)", my_ix[YY]-ngb_ix[YY]);
      } // at least 2D

      //
      // Now Z-dir
      //
      if (ndim>2) {
	if      ((my_ix[ZZ]-ngb_ix[ZZ]) == 1) { // neighbour below us
	  nodelist[off1][10] = pp.offsets[ZZ];
	  nodelist[off1][11] = pp.offsets[ZZ];
	}
	else if ((my_ix[ZZ]-ngb_ix[ZZ]) ==-1) { // neighbour above us
	  nodelist[off1][10] = pp.offsets[ZZ] + pp.LocalNG[ZZ];
	  nodelist[off1][11] = pp.offsets[ZZ] + pp.LocalNG[ZZ];
	}
	else if (my_ix[ZZ]==ngb_ix[ZZ])       { // neighbour level with us.
	  nodelist[off1][10] = pp.offsets[ZZ];
	  nodelist[off1][11] = pp.offsets[ZZ] + pp.LocalNG[ZZ];
	}
	else 
	  rep.error("domains don't touch! (Z-dir)", my_ix[ZZ]-ngb_ix[ZZ]);
      } // 3D

      //
      // The last 3 elements are the same for grids with same orientations.
      //
      nodelist[off1][12] = 1;
      nodelist[off1][13] = 2;
      nodelist[off1][14] = 3;
    } // loop over neighbours for this mesh
  } // loop over meshes

  //for (int s=0;s<Stot;s++) {
  //  for (int v=0; v<nl_len; v++) {
  //cout <<"zonelist["<<s<<"]["<<v<<"] = "<<zonelist[s][v]<<"\n";
  //cout <<"nodelist["<<s<<"]["<<v<<"] = "<<nodelist[s][v]<<"\n";
  //  }
  //}
  

  //
  // So now we have set nmesh, meshtypes[], Nngb[], ngb[S], back[S],
  // nnodes[S], nodelist[S][], and that is all we need.
  // 
  char datadir[strlength];
  err += DBGetDir(dbfile,datadir);
  DBSetDir(dbfile,"/");
  DBMkDir(dbfile,"Decomposition");
  DBSetDir(dbfile,"/Decomposition");
  
  DBoptlist *mma_opts = DBMakeOptlist(2);
  DBAddOption(mma_opts,DBOPT_DTIME,&SimPM.simtime);
  DBAddOption(mma_opts,DBOPT_CYCLE,&SimPM.timestep);
  
  //
  // write the multimeshadj
  // 
  err = DBPutMultimeshadj(dbfile, mma_name.c_str(), nmesh, meshtypes, Nngb, ngb, back, nnodes, nodelist, 0, 0, mma_opts);
  if (err) rep.error("dataio_silo_pllel::OutputData() error writing multimesh info",err);
  
  DBClearOptlist(mma_opts);
  DBFreeOptlist(mma_opts);
  
  //
  // Need to write a variable called "NumDomains"
  //
  int dims[ndim]; for (int i=0;i<ndim;i++) dims[i]=1;
  DBWrite(dbfile, "NumDomains", &nmesh, dims, 1, DB_INT);

  DBSetDir(dbfile,"/");
  DBSetDir(dbfile,datadir);
  
  for (int s=0;s<Stot;s++) {
    nodelist[s] = mem.myfree(nodelist[s]);
  }
#ifdef TESTING
  cout <<"Finished writing mulitmesh adjacency info.\n";
#endif

  return 0;
}  

//
// Write a mulitmesh object
// 
int dataio_silo_pllel::write_multimesh(DBfile *dbfile, ///< pointer to silo file.
				       class GridBaseClass *ggg, ///< pointer to data.
				       string mm_name, ///< multimesh  name
				       string mma_name ///< multimeshadj name.
				       )
{
  //    cout <<"Finished writing multimesh data to file...\n";
  rep.error("don't call me!","write_multimesh");
  return 0;
}

//
// Write a MRG tree object
// 
int dataio_silo_pllel::write_MRGtree(DBfile *dbfile, ///< pointer to silo file.
				     class GridBaseClass *ggg, ///< pointer to data.
				     string mm_name, ///< multimesh  name
				     string mrgt_name ///< MRG tree name.
				     )
{
  char datadir[strlength];
  int err= DBGetDir(dbfile,datadir);
  DBSetDir(dbfile,"/");

  //
  // Create the MRGtree object
  //
  int mesh_type = DB_MULTIMESH;
  int max_children = 10;
  DBoptlist *mt_opts = DBMakeOptlist(2);
  DBAddOption(mt_opts,DBOPT_DTIME,&SimPM.simtime);
  DBAddOption(mt_opts,DBOPT_CYCLE,&SimPM.timestep);
  
  DBmrgtree *tree = DBMakeMrgtree(mesh_type, 0, max_children, mt_opts);
  if (!tree) rep.error("Failed to create mrgtree!",tree);

  //
  // make top level MRG tree node named 'neighbors' (American spelling!)
  //
  char root[256]; strcpy(root, "neighbors");
  err = DBAddRegion(tree, root, 0, max_children, 0,0,0,0,0,0);
  if (err<0) rep.error("Failed to add root neighbours node",err);
  err=0;     
  DBSetCwr(tree, root);

  //
  // Next create a 'groupel map' with Nmesh segments.  Segment i will
  // be of type DB_BLOCKCENT and enumerate the neighbouring blocks of
  // block i.
  // 
  // from silo.h:
  // typedef struct _DBgroupelmap {
  //     char *name;
  //     int num_segments;
  //     int *groupel_types;
  //     int *segment_lengths;
  //     int *segment_ids; // can be  zero if ordered 0...(n-1)
  //     int **segment_data;
  //     void **segment_fracs; // set to null if not needed (hope not!)
  //     int fracs_data_type;  // set to null if previous el is null.
  // } DBgroupelmap;
  //
  DBgroupelmap gmap;
  char gmap_block[256]; strcpy(gmap_block, "gmap_block");
  int nsegs = mpiPM->nproc;
  int seg_types[nsegs];
  int seg_ids[nsegs];
  int seg_lens[nsegs];   // number of neighbours for mesh i
  int *seg_data[nsegs]; // list of neighbours for mesh i
  for (int i=0; i<nsegs; i++) {
    seg_types[i] = DB_BLOCKCENT;
    seg_ids[i] = i;
    seg_lens[i] = 0;
  }

  //
  // Each seg_data[i] refers to a mesh, so we put its neighbour ids in
  // the array.
  //
  MCMDcontrol pp;
  pp.nproc = nsegs;
  int ct=0;
  for (int v=0; v<nsegs; v++) {
    pp.myrank = v;
    pp.decomposeDomain();
    //
    // Count how many neighbours we have
    //
    for (int i=0; i<2*ndim; i++) {
      if (pp.ngbprocs[i]>=0) seg_lens[v] ++;
    }
    //
    // Allocate an array for those neighbours.
    //
    seg_data[v] = mem.myalloc(seg_data[v], seg_lens[v]);
    //
    // populate the array.
    //
    ct=0;
    for (int i=0; i<2*ndim; i++) {
      if (pp.ngbprocs[i]>=0) {
	seg_data[v][ct] = pp.ngbprocs[i];
	ct++;
      }
    } // loop over dirs.
  } // loop over meshes.

  //
  // Group all this in a groupel map struct:
  //
  gmap.name = gmap_block;
  gmap.num_segments = nsegs;
  gmap.groupel_types = seg_types;
  gmap.segment_ids = seg_ids;
  gmap.segment_lengths = seg_lens;
  gmap.segment_data = seg_data;
  

  //
  // Define a child node of the root named 'neighborhoods'.  Under
  // this node definean array of regions, one for each block of the
  // multiblock mesh and associate the groupel map with this array of
  // regions.
  //
  char child[256]; strcpy(child, "neighborhoods");
  err += DBAddRegion(tree, child, 0, max_children, 0,0,0,0,0,0);
  if (err<0) rep.error("Failed to add neighborhoods region",err);
  err=0;
  DBSetCwr(tree, child);
  
  char *regn_names[1]; regn_names[0] = strdup("@blocklist_%03d@n");
  err += DBAddRegionArray(tree, 
			  nsegs,      // number of regions (# of meshes)
			  regn_names, // names of regions, from scheme above.
			  0,          // info_bits, unused
			  gmap.name,  // name of the groupel map (optional)
			  1,          // number of map segs per region
			  gmap.segment_ids, // list of ids (0...(n-1) here)
			  gmap.segment_lengths, // list of lengths for each neighbour list.
			  gmap.groupel_types,   // list of block types.
			  0                     // DBoptlist
			  );
  if (err<0) rep.error("Failed to add first region array",err);
  err=0; // addregion returns positive value on success!!!
  //
  // Write the map to file.
  //
  err += DBPutGroupelmap(dbfile, gmap.name, gmap.num_segments, gmap.groupel_types, 
			 gmap.segment_lengths, gmap.segment_ids, gmap.segment_data, 
			 0,0,0);
  if (err) rep.error("Failed to add first gmap",err);

  //
  // structured grid: define a 2nd groupel map with Nmesh segments.
  // Segment i will be of type DB_NODECENT and will list the slabs of
  // nodes the block i shares with each neighbour in the same order as
  // those are listed in the neighbor map.  segment will have length
  // Nngb[i]*6 (two 3-element lists with the min and max node ranges).
  // NOTE THIS IS NODES, NOT ZONES!!!
  //
  char gmap_nodes[256]; strcpy(gmap_nodes, "gmap_nodes");
  int ns_types[nsegs];
  int ns_ids[nsegs];
  int ns_lens[nsegs];   // number of neighbours for mesh i
  int *ns_data[nsegs]; // list of neighbours for mesh i
  for (int i=0; i<nsegs; i++) {
    ns_types[i] = DB_NODECENT;
    ns_ids[i] = i;
    ns_lens[i] = 0;
  }

  //
  // loop over domains again.
  //
  pp.nproc = nsegs;
  ct=0;
  for (int v=0; v<nsegs; v++) {
    pp.myrank = v;
    pp.decomposeDomain();

    //
    // Count how many neighbours we have; length is 6 ints per
    // neighbour. (# neighbours from prev. map).
    //
    ns_lens[v] = seg_lens[v]*6;

    //
    // Allocate an array for those neighbours.
    //
    ns_data[v] = mem.myalloc(ns_data[v], ns_lens[v]);

    //
    // populate the array.
    //
    ct=0;
    enum direction dir;
    for (int i=0; i<2*ndim; i++) {
      dir = static_cast<direction>(i);
      if (pp.ngbprocs[i]>=0) {
	//
	// xmin (nodes run from 0...N)
	//
	if (dir==XP) ns_data[v][ct  ] = pp.LocalNG[XX];
	else         ns_data[v][ct  ] = 0;
	//
	// ymin
	//
	if (dir==YP) ns_data[v][ct+1] = pp.LocalNG[YY];
	else         ns_data[v][ct+1] = 0;
	//
	// zmin
	//
	if (dir==ZP) ns_data[v][ct+2] = pp.LocalNG[ZZ];
	else         ns_data[v][ct+2] = 0;
	//
	// xmax
	//
	if (dir==XN) ns_data[v][ct+3] = 0;
	else         ns_data[v][ct+3] = pp.LocalNG[XX];
	//
	// ymax
	//
	if (dir==YN) ns_data[v][ct+4] = 0;
	else         ns_data[v][ct+4] = pp.LocalNG[YY];
	//
	// zmax
	//
	if (dir==ZN) ns_data[v][ct+5] = 0;
	else         ns_data[v][ct+5] = pp.LocalNG[ZZ];
	//
	// increment counter for next neighbour.
	//
	ct += 6;
      } // if neighbour exists
    } // loop over dirs.

    //
    // Check we got all the neighbours.
    //
    if (ct != ns_lens[v])
      rep.error("number of neighbours doesn't match!",ct-ns_lens[v]);

  } // loop over meshes.
  
  //
  // Group all this in a groupel map struct:
  //
  DBgroupelmap gmap2;
  gmap2.name = gmap_nodes;
  gmap2.num_segments = nsegs;
  gmap2.groupel_types = ns_types;
  gmap2.segment_ids = ns_ids;
  gmap2.segment_lengths = ns_lens;
  gmap2.segment_data = ns_data;
  //
  // Write this second groupel map:
  //
  err += DBPutGroupelmap(dbfile, gmap2.name, gmap2.num_segments, gmap2.groupel_types, 
			 gmap2.segment_lengths, gmap2.segment_ids, gmap2.segment_data, 
			 0,0,0);
  if (err)
    rep.error("Failed to add second gmap",err);
 
  //
  // Maybe add this as a second region???
  //
  char *regn_names2[1]; regn_names2[0] = strdup("@nodelist_%03d@n");
  err += DBAddRegionArray(tree, 
			  nsegs,       // number of regions (# of meshes)
			  regn_names2, // names of regions, from scheme above.
			  0,           // info_bits, unused
			  gmap2.name,   // name of the groupel map (optional)
			  1,           // number of map segs per region
			  gmap2.segment_ids,     // list of ids (0...(n-1) here)
			  gmap2.segment_lengths, // list of lengths for each neighbour list.
			  gmap2.groupel_types,   // list of block types.
			  0                      // DBoptlist
			  );
  if (err<0) rep.error("Failed to add second region array",err);
  err=0; // addregion returns positive value on success!!!

  //
  // I think this should be the full mrgtree now, so write it to file.
  //
  char tname[256]; strcpy(tname, mrgt_name.c_str());
  char mname[256]; strcpy(mname, mm_name.c_str());
  err += DBPutMrgtree(dbfile, tname, mname, tree, mt_opts);
  if (err<0) rep.error("Failed to write MRG tree to file",err);
  err = 0;

  //
  // Delete data:
  //
  DBClearOptlist(mt_opts);
  DBFreeOptlist(mt_opts);
  DBFreeMrgtree(tree);
  for (int v=0;v<nsegs;v++) {
    seg_data[v] = mem.myfree(seg_data[v]);
    ns_data[v]  = mem.myfree(ns_data[v]);
  }
  free(regn_names[0]);
  free(regn_names2[0]);

  DBSetDir(dbfile,"/");
  DBSetDir(dbfile,datadir);
  return err;
}
#endif // if SILO
