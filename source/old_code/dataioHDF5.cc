/** \file dataioHDF5.cc
 * \author Jonathan Mackey
 * 
 * This file is not used, but has some commented out code for reading
 * and writing hdf5 data and header info.  I'm keeping it around in 
 * case I resurrect it at some stage.
 * 
 * */

/*
void LFMethod::myh5error(int e, string msg)
{
  if (e<0) {
    cerr <<"Error in HDF5 command: "<<msg<<"\t...exiting.\n";
    exit(1);
  }
  // Following line for testing/diagnostics only.
  //else cout <<"message: "<<message<<"\t code: "<<code<<endl;
}


int LFMethod::getParametersHDF5(string infile)
{
  // Opens the HDF5 file and gets the parameters from it.
  cout <<"(LFMethod::get_parameters) from file "<<infile<<" starting.\n";
  herr_t myh5_err=0;
  
  // Set datatype for Grid Parameter Class.
  hid_t myh5t_intNDVec, myh5t_dblNDVec; 
  hsize_t NDArray[1] = {3};
  if ( (myh5t_intNDVec=H5Tarray_create(H5T_NATIVE_INT, 1, NDArray, 0)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  if ( (myh5t_dblNDVec=H5Tarray_create(H5T_NATIVE_DOUBLE, 1, NDArray, 0)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  int err =0;
  hid_t myh5t_String32  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String32,32); myh5error(err,"Set string datatype");
  hid_t myh5t_String128  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String128,128); myh5error(err,"Set string datatype");
  hid_t myh5t_GP;
  if ( (myh5t_GP=H5Tcreate (H5T_COMPOUND, sizeof(GridParams))) <0)
    { cerr <<"myHDF5: Couldn't create GridParams datatype.\n"; return(1); }
  err = H5Tinsert(myh5t_GP, "Grid Type",   HOFFSET(GridParams, gridType  ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. Type",   HOFFSET(GridParams, eqntype   ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Solver Type", HOFFSET(GridParams, solverType),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Grid NDIM", HOFFSET(GridParams, ndim),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NDIM", HOFFSET(GridParams, eqnNDim ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. SimPM.nvar", HOFFSET(GridParams, nvar ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Sim Time", HOFFSET(GridParams, simtime),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Start Time", HOFFSET(GridParams, starttime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Finish Time", HOFFSET(GridParams, finishtime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Time Step", HOFFSET(GridParams, timestep),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "NGrid[3]", HOFFSET(GridParams, SimNG   ), myh5t_intNDVec);
  err+= H5Tinsert(myh5t_GP, "Total NCell",  HOFFSET(GridParams, Ncell), H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Range[3]", HOFFSET(GridParams, SimRange), myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmin[3]",  HOFFSET(GridParams, SimXmin),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmax[3]",  HOFFSET(GridParams, SimXmax),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Cell diameter",  HOFFSET(GridParams, dx), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Face Area", HOFFSET(GridParams, dA), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Vol. Area", HOFFSET(GridParams, dV), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type Of B.C.", HOFFSET(GridParams, typeOfBC), myh5t_String32);
  err+= H5Tinsert(myh5t_GP, "Total Point No.", HOFFSET(GridParams, Nbc),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Space OOA",    HOFFSET(GridParams, spOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Time  OOA",    HOFFSET(GridParams, tmOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "ArtVisc Flag", HOFFSET(GridParams, artviscosity),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Gas Cnst Gamma", HOFFSET(GridParams, gamma), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "CFL parameter",  HOFFSET(GridParams, CFL),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Visc parameter", HOFFSET(GridParams, etav),  H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type of O/P",    HOFFSET(GridParams, typeofop),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Output File", HOFFSET(GridParams, outFileBase), myh5t_String128);
  err+= H5Tinsert(myh5t_GP, "O/P Frequency",  HOFFSET(GridParams, opfreq),  H5T_NATIVE_INT);
  //  err+= H5Tinsert(myh5t_GP, "", HOFFSET(GridParams, ),);
  myh5error(err,"Insert Parameter into ParamType");

  cout <<"Opening File...";
  hid_t myh5_InFile = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  myh5error(myh5_InFile,"opening infile");
  // Get Datatypes for Grid Parameters, Units.
  cout <<"\t getting datatypes...";
  hid_t myh5_Group = H5Gopen(myh5_InFile, "/datatypes");
  myh5error(myh5_Group,"opening grid dataset");
  // hid_t myh5t_GP = H5Topen(myh5_Group, "myh5t_Parameters");
  //  myh5error(myh5t_GP,"Getting parameter datatype.");
  hid_t myh5t_UC = H5Topen(myh5_Group, "myh5t_Units");
  myh5error(myh5t_GP,"Getting Units datatype.");
  myh5_err = H5Gclose(myh5_Group);
  // Open Parameter Dataset and read in params to grid::gp
  cout <<"\tReading in Parameters...";
  hid_t indataset = H5Dopen(myh5_InFile, "/parameters/params");  
  myh5_err = H5Dread(indataset, myh5t_GP, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gp);  
  myh5error(myh5_err,"reading parameters dataset to memory");
  myh5_err = H5Dclose(indataset);
  // Read units to grid::uc
  cout <<"\tReading in Units...";
  indataset = H5Dopen(myh5_InFile, "/units/units");
  myh5_err = H5Dread(indataset, myh5t_UC, H5S_ALL, H5S_ALL, H5P_DEFAULT, &uc);  
  myh5error(myh5_err,"reading units dataset to memory");
  myh5_err = H5Dclose(indataset);
  cout <<" Done.\n";
  myh5_err = H5Tclose(myh5t_GP);    myh5error(myh5_err,"Closing Param dtype");
  myh5_err = H5Tclose(myh5t_UC);    myh5error(myh5_err,"Closing Unit  dtype");
  myh5_err = H5Tclose(myh5t_String32);   myh5error(myh5_err,"Closing str32 dtype");
  myh5_err = H5Tclose(myh5t_String128);  myh5error(myh5_err,"Closing st128 dtype");
  myh5_err = H5Tclose(myh5t_intNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Tclose(myh5t_dblNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Fclose(myh5_InFile);
  
  // Some Diagnostics follow.
//  cout<<"Density units: "<<uc.density<<" "<<uc.rhoVal<<endl;
//  cout<<"Simtime: "<<SimPM.simtime<<endl;
//  cout<<"Starttime: "<<SimPM.starttime<<endl;
//  cout<<"BCS: "<<SimPM.typeOfBC<<endl;
//  cout<<"Outputfilebase: "<<SimPM.outFileBase<<endl;
  
  // Copy GridParams data to some code variables:
//  LFMethod::typeofic = "ICFILE"; // Just set this.
//  LFMethod::typeofbc = SimPM.typeOfBC;
//  LFMethod::outpath  = SimPM.outFileBase;
//  LFMethod::spOOA = static_cast<enum ooa>(SimPM.spOOA);
//  LFMethod::tmOOA = static_cast<enum ooa>(SimPM.tmOOA);
  SimPM.dt = 0.;
  SimPM.maxtime = false;
  
  // Check that all the defines are set to what they should be:
#if defined (UNIFORM_FV)
  if (SimPM.gridType!=1)
    rep.error("(getParametersHDF5) Wrong Gridtype for type defined in compilation.",SimPM.gridType);
#else
# error "UNIFORM_FV not defined in Makefile, this is required!\n"
#endif
#if defined (HD)
  if (SimPM.eqntype!=1) rep.warning("(getParametersHDF5) Equations not set to HD in restartfile.",SimPM.eqntype,1);
#elif defined (MHD)
  if (SimPM.eqntype!=2) rep.warning("(getParametersHDF5) Equations not set to MHD in restartfile.",SimPM.eqntype,2);
#else
# error "Equations not defined in Makefile, this is required!\n"
#endif
#if defined (LAX_FRIEDRICHS)
  if (SimPM.solverType!=1) rep.warning("(getParametersHDF5) Solver set to LAX_FRIEDRICHS in makefile, but not in restartfile.",1,SimPM.solverType);  
#elif defined (GODUNOV)
  if (SimPM.solverType!=2) rep.warning("(getParametersHDF5) Solver set to GODUNOV in makefile, but not in restartfile.",2,SimPM.solverType);
#else
# error "Solver not defined in Makefile, this is required!\n"
#endif
  if (SimPM.ndim  != SimPM.ndim) rep.error("(getParametersHDF5) Dimensions of Grid wrong in makefile; should be",SimPM.ndim);
  if (EQNDIM != SimPM.eqnNDim)  rep.error("(getParametersHDF5) Dimensions of Equations wrong in makefile; should be",SimPM.eqnNDim);
  if (SimPM.nvar   != SimPM.nvar)  rep.warning("(getParametersHDF5) Number of Variables in ic-file",SimPM.nvar,SimPM.nvar);
  cout <<"\t\tnvar: "<<SimPM.nvar<<endl;

  // Check that all parameters are set
  if (SimPM.simtime <0) {
    rep.warning("(getParametersHDF5) simtime not set, setting to zero. THIS MAY BE A BIG ERROR",1.,SimPM.simtime);
    SimPM.simtime = 0.;
  }
  if (SimPM.starttime <0) {
    rep.warning("(getParametersHDF5) startime <0, PROBABLY AN ERROR!",0.,SimPM.starttime);
    SimPM.starttime = 0.;
  }
  if (SimPM.finishtime <=0 ) {rep.warning("(getParametersHDF5) finishtime not set.",1.,SimPM.finishtime);SimPM.finishtime=-1.;}
  if (SimPM.timestep <0) {
    rep.warning("(getParametersHDF5) timestep <0, PROBABLY AN ERROR!",1,SimPM.timestep);
    SimPM.timestep = 0;
  }
  if (SimPM.SimNG[0]<0) rep.error("(getParametersHDF5) SimNG[0]<0 -- must not have read from file",SimPM.SimNG[0]);
  if (SimPM.ndim>1)
    if (SimPM.SimNG[1]<0) rep.error("(getParametersHDF5) SimNG[1]<0 -- must not have read from file",SimPM.SimNG[1]);
  if (SimPM.ndim>2)
    if (SimPM.SimNG[2]<0) rep.error("(getParametersHDF5) SimNG[2]<0 -- must not have read from file",SimPM.SimNG[2]);
  if (SimPM.SimNcell<0)   rep.error("(getParametersHDF5) Ncell<0 -- Error reading from file",SimPM.SimNcell);
  if (SimPM.SimRange[0] <0)  rep.error("(getParametersHDF5) SimRange[0]<0 -- Error reading from file",SimPM.SimRange[0]);
  if (SimPM.ndim>1)
    if (SimPM.SimRange[1]<0) rep.error("(getParametersHDF5) SimRange[1]<0 -- must not have read from file",SimPM.SimRange[1]);
  if (SimPM.ndim>2)
    if (SimPM.SimRange[2]<0) rep.error("(getParametersHDF5) SimRange[2]<0 -- must not have read from file",SimPM.SimRange[2]);
  if (SimPM.SimXmin[0] <0)  rep.error("(getParametersHDF5) SimXmin[0]<0 -- Error reading from file",SimPM.SimXmin[0]);
  if (SimPM.ndim>1)
    if (SimPM.SimXmin[1]<0) rep.error("(getParametersHDF5) SimXmin[1]<0 -- must not have read from file",SimPM.SimXmin[1]);
  if (SimPM.ndim>2)
    if (SimPM.SimXmin[2]<0) rep.error("(getParametersHDF5) SimXmin[2]<0 -- must not have read from file",SimPM.SimXmin[2]);
  if (SimPM.SimXmax[0] <0)  rep.error("(getParametersHDF5) SimXmax[0]<0 -- Error reading from file",SimPM.SimXmax[0]);
  if (SimPM.ndim>1)
    if (SimPM.SimXmax[1]<0) rep.error("(getParametersHDF5) SimXmax[1]<0 -- must not have read from file",SimPM.SimXmax[1]);
  if (SimPM.ndim>2)
    if (SimPM.SimXmax[2]<0) rep.error("(getParametersHDF5) SimXmax[2]<0 -- must not have read from file",SimPM.SimXmax[2]);
  if (SimPM.dx <0) {
    rep.warning("(getParametersHDF5) dx<0 May be an error, setting it correctly",1.,SimPM.dx);
    setCellSize();
  }
  if (SimPM.typeOfBC ==0) rep.error("(getParametersHDF5) Didn't get type of BC from file. ",1);
  if (SimPM.Nbc <0) {
    rep.warning("(getParametersHDF5) Nbc not obtained from file. Will infer it from BCs",100,SimPM.Nbc);
  }
  if (SimPM.spOOA<0 || SimPM.tmOOA<0) {
    rep.warning("(getParametersHDF5) Order of Accuracy not set from file. Default to second order",2,SimPM.spOOA);
    SimPM.spOOA=2; SimPM.tmOOA=2;
    LFMethod::spOOA = static_cast<enum ooa>(SimPM.spOOA);
    LFMethod::tmOOA = static_cast<enum ooa>(SimPM.tmOOA);
  }
#ifdef LAX_FRIEDRICHS
  SimPM.spOOA = SimPM.tmOOA = 1;
  LFMethod::spOOA = static_cast<enum ooa>(SimPM.spOOA);
  LFMethod::tmOOA = static_cast<enum ooa>(SimPM.tmOOA);
#endif  
  cout <<"(getParametersHDF5) spOOA: "<<SimPM.spOOA<<" SimPM.tmOOA: "<<SimPM.tmOOA<<endl;
  if (SimPM.gamma<0) rep.error("(getParametersHDF5) gamma<0 -- must not have read it from file",SimPM.gamma);
  if (SimPM.CFL<0)   rep.error("(getParametersHDF5) CFL<0   -- must not have read it from file",SimPM.CFL);
  if (SimPM.artviscosity<0) rep.error("(getParametersHDF5) artviscosity<0 -- must not have read from file",SimPM.artviscosity);
  if (SimPM.etav<0 && SimPM.artviscosity>0) rep.error("(getParametersHDF5) etav<0 -- must not have read from file",SimPM.etav);
  else if (SimPM.etav<0 && SimPM.artviscosity<0) {rep.warning("(getParametersHDF5) etav<0 but art.visc. not enabled so it's ok.",static_cast<double>(SimPM.artviscosity),SimPM.etav); SimPM.etav=0.15;}
  if (SimPM.typeofop <0) {
    rep.warning("(getParametersHDF5) type of output file not set, setting to HDF5",2,SimPM.typeofop);
    SimPM.typeofop = 2;
  }
  //SimPM.typeofop = 1; //TESTING ... FORCE OUTPUT TO TXT.
  if (SimPM.opfreq <0) {
    rep.warning("(getParametersHDF5) Output frequency not set.  Setting to 100",100,SimPM.opfreq);
    SimPM.opfreq = 100;
  }
  if (SimPM.outFileBase ==0) {
    rep.warning("(getParametersHDF5) output file base not set, setting to ./noname",1,1);
    strcpy(SimPM.outFileBase,"noname"); outpath = SimPM.outFileBase;
  }
  
  for (int i=0;i<SimPM.ndim;i++) SimPM.SimRange[i] = SimPM.SimXmax[i]-SimPM.SimXmin[i];
  setCellSize();
  // All Parameters Checked.  Can now continue.
  return(0);
}


int LFMethod::assignDataFromHDF5(string infile)
{
  // Opens the HDF5 file and gets the data from it.
  // The grid parameters must be already read, and the grid set up to hold the data.

  int nel = SimPM.SimNcell;
  double *dblarr;
  int    *intarr;
  dblarr = new double [nel]; if (!dblarr) {cerr<<"memory assignment dblarr.\n";exit(1);}
  intarr = new int [nel]; if (!intarr) {cerr<<"memory assignment intarr.\n";exit(1);}
  //  hsize_t dim[1];
  //  dim[0] = nel;
  hsize_t dim[SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++) {
    dim[i] = SimPM.SimNG[i];
  } // This is for N-D data.

  //cout <<"Opening File...";
  hid_t myh5_InFile = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  myh5error(myh5_InFile,"opening infile");
  //cout <<"\tOpening Data Group...";
  hid_t myh5_Group = H5Gopen(myh5_InFile, "/data");
  myh5error(myh5_Group,"opening grid dataset");
  //cout <<"\tCreating 1D Dataspace...";
  //  hid_t myh5_Dspace = H5Screate_simple(1, dim, 0);
  hid_t myh5_Dspace = H5Screate_simple(SimPM.ndim, dim, 0);
  myh5error(myh5_Dspace,"creating dataspace");
  //cout <<"\t Creating Dataset...";
  hid_t myh5_Dset = H5Dopen(myh5_Group, "X-Position");
  myh5error(myh5_Dset,"Open x-pos data-dataset");// cout <<"Done!\n";
  cout <<"\t Reading X-Position...\n";
  herr_t myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
  myh5error(myh5_err,"Reading data From myh5_Dset");// cout <<"Done!\n";
  cell* cpt=g->FirstPt(); int i=0;
  do {cpt->x[0] = dblarr[i];i++;} while ((cpt=g->NextPt(cpt))!=0);
  rep.errorTest("Correct Number of points read",nel,i);
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  if(SimPM.ndim>1) {
    cout <<"\tReading Y-Position.\n";
    myh5_Dset = H5Dopen(myh5_Group, "Y-Position");
    myh5error(myh5_Dset,"read y-pos data-dataset");
    cpt=g->FirstPt(); i=0;
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    //cout <<"y-pos, h5err = "<<myh5_err<<endl;
    myh5error(myh5_err,"reading data to dblarr");
    do {cpt->x[1] = dblarr[i];i++;} while ((cpt=g->NextPt(cpt))!=0);
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  if(SimPM.ndim>2) {  
    cout <<"\tReading Z-Position.\n";
    myh5_Dset = H5Dopen(myh5_Group, "Z-Position");
    myh5error(myh5_Dset,"read z-pos data-dataset");
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"reading data to myh5_Dset");
    cpt=g->FirstPt(); i=0;
    do {cpt->x[2] = dblarr[i];i++;} while ((cpt=g->NextPt(cpt))!=0);
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  // read ids.
  //cout <<"Opening Dataset for IDs.\n";
  myh5_Dset = H5Dopen(myh5_Group, "ID");
  myh5error(myh5_Dset,"create data-dataset");
  //cout<<" Reading IDs...";
  myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intarr);
  //cout <<"IDs, h5err = "<<myh5_err<<endl;
  myh5error(myh5_err,"writing data to myh5_Dset");
  //  cout <<"Read... assigning to grid...";
  cpt=g->FirstPt(); i=0;
  do {cpt->id = intarr[i]; i++;} while ((cpt=g->NextPt(cpt))!=0);
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");
  //cout <<"\tDone.\n";
  // read State Variables.
  string s[8];
  s[RO]="Density"; s[PG]="Gas-Pressure"; s[VX]="X-Velocity";
  s[VY]="Y-Velocity"; s[VZ]="Z-Velocity"; s[BX]="Bx"; s[BY]="By"; s[BZ]="Bz";

  int nv = SimPM.nvar; cout <<"\tnvar: "<<nv<<endl;
  if(nv>SimPM.nvar)
    rep.error("(getParametersHDF5) Code not compiled for MHD.",SimPM.eqntype);
  else if (nv<SimPM.nvar)
    rep.warning("(getParametersHDF5) No B-Field in restartfile, solving HD problem with MHD solver.",SimPM.eqntype,2);
  
  for (int vct=0;vct<nv;vct++) {// Loop 5 or 8 times, for HD and MHD.  
    cout <<"\tOpening Dataset for "<<s[vct]<<"... ";
    myh5_Dset = H5Dopen(myh5_Group, (s[vct]).c_str());
    myh5error(myh5_Dset,"create data-dataset");
    cout <<"Reading Data... ";
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"reading data to dblarr");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");
    cout <<"Read, assigning to grid... ";
    // Assign data to grid: Assumes grid is set up for same point ordering as was written.
    cpt=g->FirstPt(); i=0; enum primitive var = static_cast<primitive>(vct);
    do {
      cpt->P[var] = dblarr[i];
      i++;
    } while ((cpt=g->NextPt(cpt))!=0);
    //  cout <<"i = "<<i<<"\t Ncell = "<<SimPM.SimNcell<<endl;
    if (i!=SimPM.SimNcell) {cerr<<"Error: too many gridpoints in output!\n";exit(100);}
    cout <<"\tDone.\n";
  }
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  myh5_err = H5Fclose(myh5_InFile);  myh5error(myh5_err,"Closing myh5_InFile file");
  delete [] dblarr;
  delete [] intarr;
  
  if(nv<SimPM.nvar) {// initialise B-field to zero
    cpt = g->FirstPt();
    do {
      cpt->P[BY]=cpt->P[BY]=cpt->P[BZ]=0.;
    } while ((cpt=g->NextPt(cpt))!=0);
  }
  
  cpt = g->FirstPt();
  do {for(int v=0;v<SimPM.nvar;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=g->NextPt(cpt))!=0);
  //  AddNoise2Data(1);
  return(0);
}
*/
/*
int LFMethod::outputHDF5Data(string outfile)
{
  // Set datatype for Grid Parameter Class.
  hid_t myh5t_intNDVec, myh5t_dblNDVec; 
  hsize_t NDArray[1] = {3};
  if ( (myh5t_intNDVec=H5Tarray_create(H5T_NATIVE_INT, 1, NDArray, 0)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  if ( (myh5t_dblNDVec=H5Tarray_create(H5T_NATIVE_DOUBLE, 1, NDArray, 0)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  int err =0;
  hid_t myh5t_String32  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String32,32); myh5error(err,"Set string datatype");
  hid_t myh5t_String128  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String128,128); myh5error(err,"Set string datatype");
  hid_t myh5t_GP;
  if ( (myh5t_GP=H5Tcreate (H5T_COMPOUND, sizeof(GridParams))) <0)
    { cerr <<"myHDF5: Couldn't create GridParams datatype.\n"; return(1); }
  err = H5Tinsert(myh5t_GP, "Grid Type",   HOFFSET(GridParams, gridType  ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. Type",   HOFFSET(GridParams, eqntype   ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Solver Type", HOFFSET(GridParams, solverType),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Grid NDIM", HOFFSET(GridParams, ndim),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NDIM", HOFFSET(GridParams, eqnNDim ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. SimPM.nvar", HOFFSET(GridParams, nvar ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Sim Time", HOFFSET(GridParams, simtime),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Start Time", HOFFSET(GridParams, starttime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Finish Time", HOFFSET(GridParams, finishtime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Time Step", HOFFSET(GridParams, timestep),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "NGrid[3]", HOFFSET(GridParams, SimNG   ), myh5t_intNDVec);
  err+= H5Tinsert(myh5t_GP, "Total NCell",  HOFFSET(GridParams, Ncell), H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Range[3]", HOFFSET(GridParams, SimRange), myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmin[3]",  HOFFSET(GridParams, SimXmin),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmax[3]",  HOFFSET(GridParams, SimXmax),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Cell diameter",  HOFFSET(GridParams, dx), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Face Area", HOFFSET(GridParams, dA), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Vol. Area", HOFFSET(GridParams, dV), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type Of B.C.", HOFFSET(GridParams, typeOfBC), myh5t_String32);
  err+= H5Tinsert(myh5t_GP, "Total Point No.", HOFFSET(GridParams, Nbc),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Space OOA",    HOFFSET(GridParams, spOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Time  OOA",    HOFFSET(GridParams, tmOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "ArtVisc Flag", HOFFSET(GridParams, artviscosity),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Gas Cnst Gamma", HOFFSET(GridParams, gamma), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "CFL parameter",  HOFFSET(GridParams, CFL),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Visc parameter", HOFFSET(GridParams, etav),  H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type of O/P",    HOFFSET(GridParams, typeofop),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Output File", HOFFSET(GridParams, outFileBase), myh5t_String128);
  err+= H5Tinsert(myh5t_GP, "O/P Frequency",  HOFFSET(GridParams, opfreq),  H5T_NATIVE_INT);
  //  err+= H5Tinsert(myh5t_GP, "", HOFFSET(GridParams, ),);
  myh5error(err,"Insert Parameter into ParamType");
  
  // DataTypes for Grid-Data and File-Data.
  int nel = SimPM.SimNcell;
  double *dblarr;
  int    *intarr;
  dblarr = new double [nel]; if (!dblarr) {cerr<<"memory assignment dblarr.\n";exit(1);}
  intarr = new int [nel]; if (!intarr) {cerr<<"memory assignment intarr.\n";exit(1);}
//  hsize_t dim[1];
//  dim[0] = nel; // This is for writing out the data in a big 1D array.
  hsize_t dim[SimPM.ndim];
  for (int i=0;i<SimPM.ndim;i++) {
    dim[i] = SimPM.SimNG[i];
  } // This is for N-D data.
  herr_t myh5_err=0;
  
  //  cout<<"Creating dataspace.\n";
  //  hid_t myh5_Dspace = H5Screate_simple(1, dim, 0);
  hid_t myh5_Dspace = H5Screate_simple(SimPM.ndim, dim, 0);
  myh5error(myh5_Dspace,"create myh5_Dspace");
  //  cout <<"Creating File.\n";
  hid_t myh5_OutputFile = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  myh5error(myh5_OutputFile,"create myh5_OutputFile outfile");
  //  cout <<"Creating /data Group.\n";
  hid_t myh5_Group  = H5Gcreate(myh5_OutputFile, "/data", 0);
  myh5error(myh5_Group,"create datagroup");
  // Write positions:
  //  cout <<"Creating Dataset x-pos.\n";
  hid_t myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/X-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
  myh5error(myh5_Dset,"create data-dataset");
  cell* cpt=g->FirstPt(); int i=0;
  do {dblarr[i] = cpt->x[0];i++;} while ((cpt=g->NextPt(cpt))!=0);
  //  cout <<"Writing x-pos to file.\n";
  myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
  myh5error(myh5_err,"writing data to myh5_Dset");
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  if(SimPM.ndim>1) {
    //    cout <<"Creating Dataset y-pos.\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/Y-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=g->FirstPt(); i=0;
    do {dblarr[i] = cpt->x[1];i++;} while ((cpt=g->NextPt(cpt))!=0);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  if(SimPM.ndim>2) {  
    //cout <<"creating Dataset z-pos.\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/Z-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=g->FirstPt(); i=0;
    do {dblarr[i] = cpt->x[2];i++;} while ((cpt=g->NextPt(cpt))!=0);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  // Write ids.
  //  cout <<"Creating Dataset for IDs.\n";
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/ID", H5T_NATIVE_INT, myh5_Dspace, H5P_DEFAULT);
  myh5error(myh5_Dset,"create data-dataset");
  cpt=g->FirstPt(); i=0;
  do {intarr[i] = cpt->id;i++;} while ((cpt=g->NextPt(cpt))!=0);
  myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intarr);
  myh5error(myh5_err,"writing data to myh5_Dset");
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    

  // Write State Variables.
  string states[11] = {"Density","Gas-Pressure","X-Velocity","Y-Velocity",
      "Z-Velocity","Specific-Eint","Total-Energy","Total-Pressure"
      ,"Bx","By","Bz"};
  for (int vct=0;vct<SimPM.nvar+3;vct++) {// Loop 8 or 11 times, for HD and MHD.
    ostringstream sdir;
    sdir << "/data/"<<states[vct]; string t2 = sdir.str();
    //    cout <<"Creating Dataset for "<<t2<<"\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, t2.c_str(), H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=g->FirstPt(); i=0; string temp; double b2=0.; double Utemp[SimPM.nvar];
    do {
      if (i>SimPM.SimNcell) {cerr<<"Error: too many gridpoints in output!\n";exit(100);}
      // Run through all the cases... don't want to assume element [vct] is a certain quantity.
      temp = states[vct];
      if      (temp=="Density")         dblarr[i] = cpt->P[RO];
      else if (temp=="X-Velocity")      dblarr[i] = cpt->P[VX];
      else if (temp=="Y-Velocity")      dblarr[i] = cpt->P[VY];
      else if (temp=="Z-Velocity")      dblarr[i] = cpt->P[VZ];
      else if (temp=="Gas-Pressure")    dblarr[i] = cpt->P[PG];
      else if (temp=="Specific-Eint")   dblarr[i] = cpt->P[PG]/(SimPM.gamma-1.)/cpt->P[RO];
      else if (temp=="Total-Energy")    {
	eqn->PtoU(cpt->P,Utemp,SimPM.gamma);
	dblarr[i] = Utemp[ERG];
      }
      else if (temp=="Total-Pressure")  {
	dblarr[i] = cpt->P[PG];
	if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
	  b2 = cpt->P[BX]*cpt->P[BX] +cpt->P[BY]*cpt->P[BY] +cpt->P[BZ]*cpt->P[BZ];
	  dblarr[i] += b2/2.;
	}
      }
      else if (temp=="Bx")              dblarr[i] = cpt->P[BX];
      else if (temp=="By")              dblarr[i] = cpt->P[BY];
      else if (temp=="Bz")              dblarr[i] = cpt->P[BZ];
      else {cerr<<"Don't know what to write to file...\n"; exit(100);}
      // End of running through all the cases.
      i++;
    } while ((cpt=g->NextPt(cpt))!=0);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
    
  // Now Write parameters to file.
  myh5_Group = H5Gcreate(myh5_OutputFile, "/parameters", 0);
  myh5error(myh5_Group,"create paramgroup");
  hsize_t paramDim[1] = {1};
  myh5_Dspace = H5Screate_simple(1, paramDim, 0);
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/parameters/params", myh5t_GP, myh5_Dspace, H5P_DEFAULT);
  myh5_err = H5Dwrite(myh5_Dset, myh5t_GP, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gp);
  myh5error(myh5_err,"writing data to myh5_ParameterDset");
  myh5_err = H5Dclose(myh5_Dset); myh5error(myh5_err,"Closing myh5_Dset dataspace");
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");

  // Unit Datatype
  hid_t myh5t_UC;
  if ( (myh5t_UC=H5Tcreate (H5T_COMPOUND, sizeof(units))) <0)
    { cerr <<"myHDF5: Couldn't create Units datatype.\n"; return(1); }
  err = H5Tinsert(myh5t_UC, "Length",   HOFFSET(units, length),   myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Density",  HOFFSET(units, density),  myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Velocity", HOFFSET(units, velocity), myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "B-Field",  HOFFSET(units, bfield),   myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Density Conv.Factor",  HOFFSET(units, rhoVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "Length Conv.Factor",   HOFFSET(units, lenVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "Velocity Conv.Factor", HOFFSET(units, velVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "BField Conv.Factor",   HOFFSET(units, magVal), H5T_NATIVE_DOUBLE);
  myh5error(err,"Insert element of units to H5 type");
  // Write units to file.
  //  cout <<"Creating group for units.\n";
  myh5_Group  = H5Gcreate(myh5_OutputFile, "/units", 0);
  myh5error(myh5_Group,"create unitgroup");
  hsize_t unitDim[1] = {1};
  //  Creating dataspace for units.
  myh5_Dspace = H5Screate_simple(1, unitDim, 0);
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/units/units", myh5t_UC, myh5_Dspace, H5P_DEFAULT);
  //  Writing dataspace for units.
  myh5_err = H5Dwrite(myh5_Dset, myh5t_UC, H5S_ALL, H5S_ALL, H5P_DEFAULT, &uc);
  myh5error(myh5_err,"writing data to myh5_ParameterDset");
  myh5_err = H5Dclose(myh5_Dset); myh5error(myh5_err,"Closing myh5_Dset dataspace");
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  
  // Commit Datatypes to file.
  myh5_Group  = H5Gcreate(myh5_OutputFile, "/datatypes", 0);
  myh5error(myh5_Group,"create TypeGroup");
  myh5_err = H5Tcommit(myh5_Group,"myh5t_Parameters", myh5t_GP);
  myh5_err = H5Tcommit(myh5_Group,"myh5t_Units",  myh5t_UC);
  myh5error(myh5_err,"Committing Types");  
  
  // Close everything
  myh5_err = H5Tclose(myh5t_GP);    myh5error(myh5_err,"Closing Param dtype");
  myh5_err = H5Tclose(myh5t_UC);    myh5error(myh5_err,"Closing Unit  dtype");
  myh5_err = H5Tclose(myh5t_String32);   myh5error(myh5_err,"Closing str32 dtype");
  myh5_err = H5Tclose(myh5t_String128);  myh5error(myh5_err,"Closing st128 dtype");
  myh5_err = H5Tclose(myh5t_intNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Tclose(myh5t_dblNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  myh5_err = H5Fclose(myh5_OutputFile);  myh5error(myh5_err,"Closing myh5_OutputFile outfile");
  delete [] dblarr;
  delete [] intarr;
  
  
  return(0);
}
*/  
  



/*
   // from icgen.cc... may or may not be useful someday...
   //
   int outputHDF5Data(string outfile)
   {   
   int ndim=grid->Ndim(); int nvar=grid->Nvar();
   
   // Set datatype for Grid Parameter Class.
   hid_t myh5t_intNDVec, myh5t_dblNDVec; 
   hsize_t NDArray[1] = {3};
  if ( (myh5t_intNDVec=H5Tarray_create(H5T_NATIVE_INT, 1, NDArray, NULL)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  if ( (myh5t_dblNDVec=H5Tarray_create(H5T_NATIVE_DOUBLE, 1, NDArray, NULL)) <0)
  {cerr<<"myHDF5: Couldn't create state vector datatype.\n"; return(1); }
  int err =0;
  hid_t myh5t_String32  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String32,32); myh5error(err,"Set string datatype");
  hid_t myh5t_String128  = H5Tcopy (H5T_C_S1);
  err = H5Tset_size (myh5t_String128,128); myh5error(err,"Set string datatype");
  hid_t myh5t_GP;
  if ( (myh5t_GP=H5Tcreate (H5T_COMPOUND, sizeof(GridParams))) <0)
    { cerr <<"myHDF5: Couldn't create GridParams datatype.\n"; return(1); }
  err = H5Tinsert(myh5t_GP, "Grid Type",   HOFFSET(GridParams, gridType  ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. Type",   HOFFSET(GridParams, eqntype   ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Solver Type", HOFFSET(GridParams, solverType),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Grid NDIM", HOFFSET(GridParams, ndim),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NDIM", HOFFSET(GridParams, eqnNDim ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NVAR", HOFFSET(GridParams, eqnNVar ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Sim Time", HOFFSET(GridParams, simtime),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Start Time", HOFFSET(GridParams, starttime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Finish Time", HOFFSET(GridParams, finishtime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Time Step", HOFFSET(GridParams, timestep),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "NGrid[3]", HOFFSET(GridParams, NG   ), myh5t_intNDVec);
  err+= H5Tinsert(myh5t_GP, "Total NCell",  HOFFSET(GridParams, Ncell), H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Range[3]", HOFFSET(GridParams, range), myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmin[3]",  HOFFSET(GridParams, xmin),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmax[3]",  HOFFSET(GridParams, xmax),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Cell diameter",  HOFFSET(GridParams, dx), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Face Area", HOFFSET(GridParams, dA), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Vol. Area", HOFFSET(GridParams, dV), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type Of B.C.", HOFFSET(GridParams, typeofbc), myh5t_String32);
  err+= H5Tinsert(myh5t_GP, "Total Point No.", HOFFSET(GridParams, Nbc),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Space OOA",    HOFFSET(GridParams, spOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Time  OOA",    HOFFSET(GridParams, tmOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "ArtVisc Flag", HOFFSET(GridParams, artviscosity),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Gas Cnst Gamma", HOFFSET(GridParams, gamma), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "CFL parameter",  HOFFSET(GridParams, CFL),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Visc parameter", HOFFSET(GridParams, etav),  H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type of O/P",    HOFFSET(GridParams, typeofop),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Output File", HOFFSET(GridParams, outFileBase), myh5t_String128);
  err+= H5Tinsert(myh5t_GP, "O/P Frequency",  HOFFSET(GridParams, opfreq),  H5T_NATIVE_INT);
  //  err+= H5Tinsert(myh5t_GP, "", HOFFSET(GridParams, ),);
  myh5error(err,"Insert Parameter into ParamType");
  
  // DataTypes for Grid-Data and File-Data.
  //  cout <<"Setting up data arrays for writing.\n";
  int nel = SimPM.Ncell;
  double *dblarr;
  int    *intarr;
  dblarr = new double [nel]; if (!dblarr) {cerr<<"memory assignment dblarr.\n";exit(1);}
  intarr = new int [nel]; if (!intarr) {cerr<<"memory assignment intarr.\n";exit(1);}

  hsize_t dim[ndim];
  for (int i=0;i<ndim;i++) {
    dim[i] = SimPM.NG[i];
  } // This is for N-D data.
  herr_t myh5_err=0;
  
  //  cout<<"Creating dataspace.\n";
  //  hid_t myh5_Dspace = H5Screate_simple(1, dim, NULL);
  hid_t myh5_Dspace = H5Screate_simple(ndim, dim, NULL);
  myh5error(myh5_Dspace,"create myh5_Dspace");
  //  cout <<"Creating File.\n";
  hid_t myh5_OutputFile = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  myh5error(myh5_OutputFile,"create myh5_OutputFile outfile");
  //  cout <<"Creating /data Group.\n";
  hid_t myh5_Group  = H5Gcreate(myh5_OutputFile, "/data", 0);
  myh5error(myh5_Group,"create datagroup");

  // Write positions:
  //  cout <<"Creating Dataset x-pos.\n";
  hid_t myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/X-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
  myh5error(myh5_Dset,"create data-dataset");
  cell* cpt=grid->FirstPt(); int i=0;
  do {dblarr[i] = cpt->x[0];i++;} while ((cpt=grid->NextPt(cpt))!=NULL);
  //  cout <<"Writing x-pos to file.\n";
  myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
  myh5error(myh5_err,"writing data to myh5_Dset");
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  if(ndim>1) {
  //    cout <<"Creating Dataset y-pos.\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/Y-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=grid->FirstPt(); i=0;
    do {dblarr[i] = cpt->x[1];i++;} while ((cpt=grid->NextPt(cpt))!=NULL);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  if(ndim>2) {  
  //    cout <<"Creating Dataset z-pos.\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/Z-Position", H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=grid->FirstPt(); i=0;
    do {dblarr[i] = cpt->x[2];i++;} while ((cpt=grid->NextPt(cpt))!=NULL);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }

  // Write ids.
  //  cout <<"Creating Dataset for IDs.\n";
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/data/ID", H5T_NATIVE_INT, myh5_Dspace, H5P_DEFAULT);
  myh5error(myh5_Dset,"create data-dataset");
  cpt=grid->FirstPt(); i=0;
  do {intarr[i] = cpt->id;i++;} while ((cpt=grid->NextPt(cpt))!=NULL);
  myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intarr);
  myh5error(myh5_err,"writing data to myh5_Dset");
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    

  // Write State Variables.
  string states[8] = {"Density","Gas-Pressure","X-Velocity","Y-Velocity",
      "Z-Velocity","Bx","By","Bz"};
  for (int vct=0;vct<nvar;vct++) {// Loop 5 or 8 times, for HD and MHD.
    ostringstream sdir;
    sdir << "/data/"<<states[vct]; string t2 = sdir.str();
    //    cout <<"Creating Dataset for "<<t2<<"\n";
    myh5_Dset = H5Dcreate(myh5_OutputFile, t2.c_str(), H5T_NATIVE_DOUBLE, myh5_Dspace, H5P_DEFAULT);
    myh5error(myh5_Dset,"create data-dataset");
    cpt=grid->FirstPt(); i=0; string temp;
    do {
      if (i>SimPM.Ncell) {cerr<<"Error: too many gridpoints in output!\n";exit(100);}
      // Run through all the cases... don't want to assume element [vct] is a certain quantity.
      temp = states[vct];
      if      (temp=="Density")         dblarr[i] = cpt->P[RO];
      else if (temp=="X-Velocity")      dblarr[i] = cpt->P[VX];
      else if (temp=="Y-Velocity")      dblarr[i] = cpt->P[VY];
      else if (temp=="Z-Velocity")      dblarr[i] = cpt->P[VZ];
      else if (temp=="Gas-Pressure")    dblarr[i] = cpt->P[PG];
      else if (temp=="Bx")              dblarr[i] = cpt->P[BX];
      else if (temp=="By")              dblarr[i] = cpt->P[BY];
      else if (temp=="Bz")              dblarr[i] = cpt->P[BZ];
      else {cerr<<"Don't know what to write to file...\n"; exit(100);}
      // End of running through all the cases.
      i++;
    } while ((cpt=grid->NextPt(cpt))!=NULL);
    myh5_err = H5Dwrite(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"writing data to myh5_Dset");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
    
  // Now Write parameters to file.
  myh5_Group = H5Gcreate(myh5_OutputFile, "/parameters", 0);
  myh5error(myh5_Group,"create paramgroup");
  hsize_t paramDim[1] = {1};
  myh5_Dspace = H5Screate_simple(1, paramDim, NULL);
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/parameters/params", myh5t_GP, myh5_Dspace, H5P_DEFAULT);
  myh5_err = H5Dwrite(myh5_Dset, myh5t_GP, H5S_ALL, H5S_ALL, H5P_DEFAULT, gp);
  myh5error(myh5_err,"writing data to myh5_ParameterDset");
  myh5_err = H5Dclose(myh5_Dset); myh5error(myh5_err,"Closing myh5_Dset dataspace");
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");

  // Unit Datatype
  hid_t myh5t_UC;
  if ( (myh5t_UC=H5Tcreate (H5T_COMPOUND, sizeof(Units))) <0)
    { cerr <<"myHDF5: Couldn't create Units datatype.\n"; return(1); }
  err = H5Tinsert(myh5t_UC, "Length",   HOFFSET(Units, length),   myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Density",  HOFFSET(Units, density),  myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Velocity", HOFFSET(Units, velocity), myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "B-Field",  HOFFSET(Units, bfield),   myh5t_String32);
  err+= H5Tinsert(myh5t_UC, "Density Conv.Factor",  HOFFSET(Units, rhoVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "Length Conv.Factor",   HOFFSET(Units, lenVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "Velocity Conv.Factor", HOFFSET(Units, velVal), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_UC, "BField Conv.Factor",   HOFFSET(Units, magVal), H5T_NATIVE_DOUBLE);
  myh5error(err,"Insert element of Units to H5 type");
  // Write units to file.
  //  cout <<"Creating group for Units.\n";
  myh5_Group  = H5Gcreate(myh5_OutputFile, "/units", 0);
  myh5error(myh5_Group,"create unitgroup");
  hsize_t unitDim[1] = {1};
  //  Creating dataspace for units.
  myh5_Dspace = H5Screate_simple(1, unitDim, NULL);
  myh5_Dset = H5Dcreate(myh5_OutputFile, "/units/units", myh5t_UC, myh5_Dspace, H5P_DEFAULT);
  //  Writing dataspace for units.
  myh5_err = H5Dwrite(myh5_Dset, myh5t_UC, H5S_ALL, H5S_ALL, H5P_DEFAULT, uc);
  myh5error(myh5_err,"writing data to myh5_ParameterDset");
  myh5_err = H5Dclose(myh5_Dset); myh5error(myh5_err,"Closing myh5_Dset dataspace");
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  
  // Commit Datatypes to file.
  myh5_Group  = H5Gcreate(myh5_OutputFile, "/datatypes", 0);
  myh5error(myh5_Group,"create TypeGroup");
  myh5_err = H5Tcommit(myh5_Group,"myh5t_Parameters", myh5t_GP);
  myh5_err = H5Tcommit(myh5_Group,"myh5t_Units",  myh5t_UC);
  myh5error(myh5_err,"Committing Types");  
  
  // Close everything
  myh5_err = H5Tclose(myh5t_GP);    myh5error(myh5_err,"Closing Param dtype");
  myh5_err = H5Tclose(myh5t_UC);    myh5error(myh5_err,"Closing Unit  dtype");
  myh5_err = H5Tclose(myh5t_String32);   myh5error(myh5_err,"Closing str32 dtype");
  myh5_err = H5Tclose(myh5t_String128);  myh5error(myh5_err,"Closing st128 dtype");
  myh5_err = H5Tclose(myh5t_intNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Tclose(myh5t_dblNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  myh5_err = H5Fclose(myh5_OutputFile);  myh5error(myh5_err,"Closing myh5_OutputFile outfile");
  
  return(0);
  }
*/
