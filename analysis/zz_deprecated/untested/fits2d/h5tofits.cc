
#include <iostream>
#include <fstream>
#include <cmath>
#include "hdf5.h"
#include "fitsio.h"
using namespace std;

#define MACHINEACCURACY 5.e-16 ///< This should be found out before running the code!

#define GNDIM 2
#define EQNDIM 3

#if defined (HD)
#define NVAR (2+EQNDIM)
#elif defined (MHD)
#define NVAR (2+2*EQNDIM) // density, pressure, plus velocity and magnetic field.
#else
#error HD and MHD not defined! Please specify which equations are to be solved.
#define NVAR -1 /**< \brief Number of variables in the equations.*/
#endif


/** \brief Primitive Variables Enum 
 * The code requires that the velocities follow each other in order,
 * and the same for the Magnetic Field.  I'm not sure if it will work 
 * without them being in that order, but I'm writing the code assuming 
 * they are in order.
 */
enum primitive {
    RO=0, ///< Density
    PG=1, ///< Gas pressure
    VX=2, ///< x-velocity
    VY=3, ///< y-velocity
    VZ=4, ///< z-velocity
    BX=5, ///< x-magnetic field
    BY=6, ///< y-magnetic field
    BZ=7  ///< z-magnetic field
};
/** \brief Conserved Variables Enum */
enum conserved {
    RHO=0,  ///< Density
    ERG=1,  ///< Total Energy
    MMX=2,  ///< x-momentum
    MMY=3,  ///< y-momentum
    MMZ=4,  ///< z-momentum
    BBX=5,  ///< x-magnetic field
    BBY=6,  ///< y-magnetic field
    BBZ=7   ///< z-magnetic field
};
enum axes {
    XX=0, ///< X-Direction.
    YY=1, ///< Y-Direction.
    ZZ=2  ///< Z-Direction.
};
enum direction{
    NO=-1,   ///< No Direction.
    XN=0,  ///< Negative x-direction x-
    XP=1,  ///< Positive x-direction x+
    YN=2,  ///< Negative y-direction y-
    YP=3,  ///< Positive y-direction y+
    ZN=4,  ///< Negative z-direction z-
    ZP=5   ///< Positive z-direction z+
}; // {x-,x+,y-,y+,z-,z+}


class GridParams {
  public :
   // Overall grid properties (some from Makefile)
   int gridType; ///< Uniform, Finite Volume, Cartesian = 1 = only option.
   int eqnType;  ///< HD=1, MHD=2.  These are the only options for now.
   int solverType; ///< Lax-Friedrichs=1, Godunov-RS=2.
   int gridNDim; ///< Dimensionality of grid. Must match GNDIM.
   int eqnNDim;  ///< Dimensionality of equations, set to 3 for now.
   int eqnNVar;  ///< Number of variables in equations,  =5(HD) or 8(MHD).
   // Timing
   double simtime; /**< \brief current time in simulation. */
   double starttime; /**< \brief initial time to start simulation at. */
   double finishtime; ///< Time at which to finish the simulation.
   int timestep; /**< \brief Integer count of the number of timesteps since the start. */
   // Grid Point data
   int NG[GNDIM]; ///< Number of 'real' grid-points in each direction.
   int Ncell;  ///< Total number of 'real' grid points (within the range).
   double range[GNDIM]; ///< Size of domain in x,y,z-direction.
   double xmin[GNDIM];  ///< Min value of x,y,z in domain.
   double xmax[GNDIM];  ///< Max value of x,y,z in domain.
   double dx;  ///< Linear side length of (uniform, cubic, cartesian) grid cells.
   double dA;  ///< Area of one surface of the (uniform, cubic, cartesian) grid cells.
   double dV;  ///< Volume of one cubic, cartesian grid cell (same for all cells).
   // Boundary cell data.
   char typeOfBC[32];  /**< \brief Type of boundary condition.*/
   int Nbc;   /**< \brief Depth of boundary/ghost cells perpendicular to edge.*/
   // Integration accuracy
   int spOOA;   ///< Spatial Order of Accuracy in the code.
   int tmOOA;   ///< Time Order of Accuracy in the code.
   // Physics
   double gamma;     /**< \brief Ideal gas constant.*/
   double CFL;       /**< \brief Courant factor, must be less than one, should be about 0.5.*/
   int artviscosity; ///< Integer flag. 0=No Artificial Viscosity; 1=Falle's version; Others??
   double etav;      ///< Artificial viscosity coefficient, should be between 0.01 and 0.3
   // File I/O
   int typeOfOP;   ///< Integer flag. 1=ascii, 2=hdf5, etc.
   char outFileBase[128];  ///< Filename, with path, to write data to.
   int opfreq;         ///< Output file every 'opfreq'th timestep.
};


/** \brief Data Point for finite-volume grid. 
 * Stripped down grid point, with just position, primitve state vector,
 * and an id.
 * */
class cell{
 public:
   double x[GNDIM]; /**< \brief Position of centre of cell */
   double P[NVAR];  /**< \brief Primitive Variables. For 1D-HD they are \f$ [\rho, u, p]\f$. */
   int id;          /**< \brief Grid point's id. */
   cell *ngb[2*GNDIM]; ///< Pointers to cells Neigbours.
   cell *npt;  ///< Pointer to next point in list, or some other cell.
};


class tempgrid {
  private:
   const int nc, nx, ny;
   const double dx;
  public:
   tempgrid(int nt, int nxx, int nyy, double dxx)
     : nc(nt), nx(nxx), ny(nyy), dx(dxx)
     {
       data = new cell [nc];
       for(int i=0;i<nc;i++) data[i].id=i;
       if(nx*ny != nc) {cerr<<"bad nx,ny,nc\n";exit(1);}
     }
   ~tempgrid() {delete [] data;}
   
   cell *data;
   cell *pt; // for use as a general purpose pointer.
   
   cell * firstPt() {return (&data[0]);}
   
   cell * nextPt(cell *cpt) {
     if((cpt->id +1)<nc) return (&data[cpt->id +1]);
     else return(0);
   }
   
   int assignNGB() {
     cell *cpt = firstPt();
     do {
       int t;
       if( (t=cpt->id+nx) <nc) cpt->ngb[YP] = &data[t];
       if( (t=cpt->id-nx) >=0) cpt->ngb[YN] = &data[t];
       if( (t= (cpt->id+1)%nx)    > cpt->id%nx) cpt->ngb[XP] = &data[cpt->id+1];
       if( (t= (cpt->id-1+nx)%nx) < cpt->id%nx) cpt->ngb[XN] = &data[cpt->id-1];
//       cout <<"id="<<cpt->id<<" me:"<<cpt<<"\t ngb :";
//       for (int i=0; i<2*GNDIM; i++) cout <<"\t"<<cpt->ngb[i]; cout <<"\n";
     } while ( (cpt=nextPt(cpt)) !=0);
     return(0);
   }
   cell * nextPt(const cell *cpt, const direction dir) 
   {    return(cpt->ngb[dir]);   }
   
   int assignDataFromHDF5(string);
   
   double divB(const cell *cpt)
   { // get divergence of B
     double divb = nextPt(cpt,XP)->P[BX] - nextPt(cpt,XN)->P[BX];
     if (GNDIM>1) divb += nextPt(cpt,YP)->P[BY] - nextPt(cpt,YN)->P[BY];
     if (GNDIM>2) divb += nextPt(cpt,ZP)->P[BZ] - nextPt(cpt,ZN)->P[BZ];
     divb /= 2.*dx;
     return (divb);
   }
   
   void rotate(double *v, double theta)
   {// Rotate by theta in x-y plane.
     double ct=cos(theta); double st=sin(theta);
     
     double vx = v[VX]*ct - v[VY]*st;
     double vy = v[VX]*st + v[VY]*ct;
     v[VX] = vx; v[VY]=vy;
#ifdef MHD     
     vx = v[BX]*ct - v[BY]*st;
     vy = v[BX]*st + v[BY]*ct;
     v[BX] = vx; v[BY]=vy;
#endif
   }

};


double distance(double *p1, double *p2) {
  double d=0.;  for (int i=0;i<GNDIM;i++) d+=(p1[i]-p2[i])*(p1[i]-p2[i]);
  return(sqrt(d));
}

bool equalD(const double a, const double b)
{
  if(a==b) return(true);
  if(fabs(a+b)<1.e-100) {
    cout <<"tiny numbers in equalD(a,b); a,b <1.e-100... a="<<a<<", b="<<b<<"; returning true.\n";
    return(true);
  }
  if( (fabs(a-b)/fabs(a+b+1.e-100)) </*10.* MACHINEACCURACY*/1.e-6) return(true); // true is 1
  else return(false); // false is zero.
}


// Declarations
int getParametersHDF5(string);
int assignDataFromHDF5(string);

GridParams gp;


int main(int argc, char **argv)
{
  if (argc!=4) {
    cerr << "Main: must call with 3 arguments... <main> <restartfile.h5> <outfile> <param>"<< endl;
    exit(1);
  }
  string infile = argv[1];
  string opfile = argv[2];
  int param=atoi(argv[3]);

  
  int err=0;
  err =  getParametersHDF5(infile);
  if (err!=0) {cerr<<"failed to get parameters.\n"; return(1);}
  const int ncell = gp.Ncell;
  cout <<"\tNcell: "<<gp.Ncell<<"\t"<<gp.NG[XX]<<"\t"<<gp.NG[YY]<<"\t"<<endl;
  tempgrid g(ncell,gp.NG[XX], gp.NG[YY], gp.dx);
  g.assignNGB();
  g.assignDataFromHDF5(infile);
  
  // Now I have a grid of data in 'g', specified by parameters in 'gp',
  // exactly as I would have in my grid code.
  // 

  fitsfile *ff; ///< fits file pointer.
  int status =0, nkeys;
  err=0;
  cout <<"opening fits file...";
  // This creates a new file with name given.  status is non-zero if it fails,
  // which happens if the file exists.
  fits_create_file(&ff, opfile.c_str(), &status);
  
  // Opens a fits file, can use READWRITE or READONLY as options.  If rw, then
  // it will create a new file if it doesn't exist, if ro then it fails if the
  // file doesn't exist.
  //  fits_open_file(&ff, opfile.c_str(), READWRITE, &status);
  if(status) {cerr<<"file open went bad.\n";exit(1);}
  cout <<"done.\n";

  
  double var;
  long int pix[GNDIM];  for (int i=0;i<GNDIM;i++) pix[i]=gp.NG[i];
  fits_create_img(ff,DOUBLE_IMG,GNDIM,pix,&status);
  err += fits_write_key(ff, TINT, "gridType"  ,&gp.gridType, "UniformGrid=1", &status);
  err += fits_write_key(ff, TINT, "eqnType"   ,&gp.eqnType,  "HD=1,MHD=2", &status);
  err += fits_write_key(ff, TINT, "solverType",&gp.solverType, "LF=1,GD=2", &status);
  err += fits_write_key(ff, TINT, "gridNDim",  &gp.gridNDim, "Num. spatial dimensions of grid", &status);
  err += fits_write_key(ff, TINT, "eqnNDim",   &gp.eqnNDim,  "Num. spatial dimensions in vectors", &status);
  err += fits_write_key(ff, TINT, "eqnNVar",   &gp.eqnNVar,  "Num. variables in equations", &status);
  err += fits_write_key(ff, TDOUBLE, "simtime",&gp.simtime, "current time in simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "starttime",&gp.starttime, "initial time to start simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "finishtime",&gp.finishtime, "Time to finish the simulation", &status);
  err += fits_write_key(ff, TINT, "timestep",&gp.timestep, "TIMESTEP SINCE START OF SIM", &status);
  err += fits_write_key(ff, TINT, "NG[XX]",&gp.NG[XX], "Number of grid-points in X-direction", &status);
  if(GNDIM>1)
    err += fits_write_key(ff, TINT, "NG[YY]",&gp.NG[YY], "Number of grid-points in Y-direction", &status);
  if(GNDIM>2)
    err += fits_write_key(ff, TINT, "NG[ZZ]",&gp.NG[ZZ], "Number of grid-points in Z-direction", &status);
  err += fits_write_key(ff, TINT, "Ncell",&gp.Ncell, "Total Number of Grid Points", &status);
  g.pt = g.firstPt();
  for (int i=0;i<GNDIM;i++) pix[i]=1;
  do {
    pix[XX] = static_cast<int>( (g.pt->x[XX]+gp.dx/1.99)/gp.dx); //truncates to integer part
    pix[YY] = static_cast<int>( (g.pt->x[YY]+gp.dx/1.99)/gp.dx);
    //fits_write_col(*fptr, datatype, colnum, firstpix, firstelem, nelements, data, &status)
    if (param==-1) {
      var = g.pt->P[PG]/g.pt->P[RO]/(gp.gamma -1.);
      err += fits_write_pix(ff, TDOUBLE, pix, 1, &var, &status);
    }
    else if(param<NVAR)
      err += fits_write_pix(ff, TDOUBLE, pix, 1, &(g.pt->P[param]), &status);
    else {cerr<<"bad param\n";exit(1);}
  } while ( (g.pt=g.nextPt(g.pt))!=0 );
  if(status) fits_report_error(stderr,status);
  
  
  /*
  // Write a second image if you want to.
  fits_create_img(ff,DOUBLE_IMG,GNDIM,pix,&status);
  g.pt = g.firstPt();
  for (int i=0;i<GNDIM;i++) pix[i]=1;
  do {
    pix[XX] = static_cast<int>( (g.pt->x[XX]+gp.dx/1.99)/gp.dx);
    pix[YY] = static_cast<int>( (g.pt->x[YY]+gp.dx/1.99)/gp.dx);
    //fits_write_col(*fptr, datatype, colnum, firstpix, firstelem, nelements, data, &status)
    err += fits_write_pix(ff, TDOUBLE, pix, 1, &(g.pt->P[PG]), &status);
  } while ( (g.pt=g.nextPt(g.pt))!=0 );
  if(status) fits_report_error(stderr,status);
  */
  
/*  int  fits_create_tbl(fitsfile ptr.to.file, ///< fits file pointer.
		  int tbltype=bin/acs, ///< 
		  int num.reserved.rows,///< should be set to zero unless i have a good reason not to.
		  int num.cols, ///< number of columns of table.
		  char** col.name.arr, ///< array of c_strings with names of columns
		  char** col.datattype.arr, ///< array of c_strings with datatypes J=int,E=float,D=double
		  char** col.phys.units.arr, ///< array of phsyical units, can be null
		  cahr* table.name, ///< name of table, can be null.
		  int &status ///< error code for returning.
		  )
*/
  /*
  // Write a binary table if you want to.
  int nels=32;
  int ncols = 3+NVAR;
#ifdef MHD
  ncols +=1; // for divB
#endif
  string colnames[] = {"X","Y","Density","Gas-Pressure","X-Velocity","Y-Velocity","Z-Velocity","InternalEnergy","Bx","By","Bz","divB"};
  string coltypes[] = {"1D","1D","1D","1D","1D","1D","1D","1D","1D","1D","1D","1D"};
//  if(GNDIM==2) coltypes[0] = "2D"; if(GNDIM==3) coltypes[0] = "3D";
  string siunits[] = {"m","m","kg/m^3","kg/m/s^2","m/s","m/s","m/s","J","T","T","T","T/m"};
  cout <<"setup string arrays\n";
  cout <<"position vector: "<<coltypes[0]<<endl;
  char *ttype[ncols];
  char *tform[ncols];
  char *tunit[ncols];
  for (int i=0; i<ncols; i++) {
    ttype[i] = new char [nels]; tform[i] = new char [nels]; tunit[i] = new char [nels];
  }
  for (int i=0; i<ncols; i++) {
    strcpy(ttype[i],colnames[i].c_str());
    strcpy(tform[i],coltypes[i].c_str());
    strcpy(tunit[i],siunits[i].c_str() );
  }
  char tblname[nels]; strcpy(tblname,"mytable");
  err += fits_create_tbl(ff, BINARY_TBL, 0, ncols,ttype,tform,tunit,tblname,&status);
  //  err += fits_modify_vector_len(ff,1,GNDIM,&status); // change 1st column to n-dim vector.
  
  
  err += fits_write_key(ff, TINT, "gridType"  ,&gp.gridType, "UniformGrid=1", &status);
  err += fits_write_key(ff, TINT, "eqnType"   ,&gp.eqnType,  "HD=1,MHD=2", &status);
  err += fits_write_key(ff, TINT, "solverType",&gp.solverType, "LF=1,GD=2", &status);
  err += fits_write_key(ff, TINT, "gridNDim",  &gp.gridNDim, "Num. spatial dimensions of grid", &status);
  err += fits_write_key(ff, TINT, "eqnNDim",   &gp.eqnNDim,  "Num. spatial dimensions in vectors", &status);
  err += fits_write_key(ff, TINT, "eqnNVar",   &gp.eqnNVar,  "Num. variables in equations", &status);
  err += fits_write_key(ff, TDOUBLE, "simtime",&gp.simtime, "current time in simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "starttime",&gp.starttime, "initial time to start simulation", &status);
  err += fits_write_key(ff, TDOUBLE, "finishtime",&gp.finishtime, "Time to finish the simulation", &status);
  err += fits_write_key(ff, TINT, "timestep",&gp.timestep, "TIMESTEP SINCE START OF SIM", &status);
  err += fits_write_key(ff, TINT, "NG[XX]",&gp.NG[XX], "Number of grid-points in X-direction", &status);
  if(GNDIM>1)
    err += fits_write_key(ff, TINT, "NG[YY]",&gp.NG[YY], "Number of grid-points in Y-direction", &status);
  if(GNDIM>2)
    err += fits_write_key(ff, TINT, "NG[ZZ]",&gp.NG[ZZ], "Number of grid-points in Z-direction", &status);
  err += fits_write_key(ff, TINT, "Ncell",&gp.Ncell, "Total Number of Grid Points", &status);
  // Insert rows in the binary table.
  fits_insert_rows(ff,0,gp.Ncell,&status);
  g.pt = g.firstPt();
  int row=0;
  do {
    row++; // fits uses unit offset.
    //fits_write_col(*fptr, datatype, colnum, firstrow, firstelem, nelements, data, &status)
    err += fits_write_col(ff, TDOUBLE, 1, row, 1, 1, &g.pt->x[XX], &status);
    err += fits_write_col(ff, TDOUBLE, 2, row, 1, 1, &g.pt->x[YY], &status);
    err += fits_write_col(ff, TDOUBLE, 3, row, 1, 1, &(g.pt->P[RO]), &status);
  } while ( (g.pt=g.nextPt(g.pt))!=0 );
*/
    
    err += fits_close_file(ff,&status);
  
/*  ofstream opf(opfile.c_str());
  if(!opf.is_open()) {cerr<<"outfile not open\n";exit(1);}
  opf <<"# twod2oned outputfile for data from "<<infile<<endl;
  opf <<"# format: x,rho,pg,vx,vy,vz,e_int,[Bx,By,Bz,divB]\n";
  opf.setf( ios_base::scientific );
  opf.precision(6);

  double *origin = new double[GNDIM]; for (int i=0;i<GNDIM;i++) origin[i]=0.;
  
  g.pt = g.firstPt();
  do {
    if ((g.nextPt(g.pt,XN)!=0) && (g.nextPt(g.pt,YP)!=0) && !equalD(g.pt->P[PG],(g.nextPt(g.nextPt(g.pt,XN),YP))->P[PG])) {
      cout<<"Warning, ngb not the same! "<<g.pt->id;
      cout<<"\tcpt: "<<g.pt->P[PG]<<"\tnpt: "<<(g.nextPt(g.nextPt(g.pt,XN),YP))->P[PG]<<endl;
    }
    // rotate vectors so diagonal becomes x-axis.
    g.rotate(g.pt->P,-M_PI/4.);
    // output to file.
    opf << distance(g.pt->x,origin) << "  ";
    opf <<g.pt->P[RO] <<"  "<<g.pt->P[PG];
    opf <<"  "<< g.pt->P[VX] <<"  "<<  g.pt->P[VY] <<"  "<<  g.pt->P[VZ];
    opf <<"  "<< g.pt->P[PG]/g.pt->P[RO]/(gp.gamma-1.);
#ifdef MHD
    b2 = g.pt->P[BX]*g.pt->P[BX] +g.pt->P[BY]*g.pt->P[BY] +g.pt->P[BZ]*g.pt->P[BZ];
    opf <<"  "<< g.pt->P[BX] <<"  "<< g.pt->P[BY] <<"  "<< g.pt->P[BZ] <<"  "<< g.pt->P[PG]+b2/2.<<"  "<< g.divB(g.pt);
#endif //MHD
    opf  <<endl;
    //if (g.nextPt(g.pt,YP) !=0) {cout<<"getting next point: id="; g.pt=g.nextPt(g.pt,YP);cout <<g.pt->id<<endl;}
  } while ( (g.pt=g.nextPt(g.pt,XP))!=0 && (g.pt=g.nextPt(g.pt,YP))!=0);

  opf.setf(ios_base::fmtflags(0),ios_base::floatfield);
  opf.precision(6);
  opf.close();
 */
  return(0);
}


void myh5error(int e, string msg)
{
  if (e<0) {
    cerr <<"Error in HDF5 command: "<<msg<<"\t...exiting.\n";
    exit(1);
  }
  // Following line for testing/diagnostics only.
  //else cout <<"message: "<<message<<"\t code: "<<code<<endl;
}



int getParametersHDF5(string infile)
{
  // Opens the HDF5 file and gets the parameters from it.
  cout <<"(baseGrid::get_parameters) from file "<<infile<<" starting.\n";
  herr_t myh5_err=0;
  
  // Set datatype for Grid Parameter Class.
  hid_t myh5t_intNDVec, myh5t_dblNDVec; 
  hsize_t NDArray[1] = {GNDIM};
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
  err = H5Tinsert(myh5t_GP, "Grid Type",      HOFFSET(GridParams, gridType  ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. Type",      HOFFSET(GridParams, eqnType   ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Solver Type",    HOFFSET(GridParams, solverType),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Grid NDIM",      HOFFSET(GridParams, gridNDim),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NDIM",      HOFFSET(GridParams, eqnNDim ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Eqn. NVAR",      HOFFSET(GridParams, eqnNVar ),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Sim Time",       HOFFSET(GridParams, simtime),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Start Time",     HOFFSET(GridParams, starttime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Finish Time",    HOFFSET(GridParams, finishtime), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Time Step",      HOFFSET(GridParams, timestep),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "NGrid[GNDIM]",   HOFFSET(GridParams, NG   ), myh5t_intNDVec);
  err+= H5Tinsert(myh5t_GP, "Total NCell",    HOFFSET(GridParams, Ncell), H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Range[GNDIM]",   HOFFSET(GridParams, range), myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmin[GNDIM]",    HOFFSET(GridParams, xmin),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Xmax[GNDIM]",    HOFFSET(GridParams, xmax),  myh5t_dblNDVec);
  err+= H5Tinsert(myh5t_GP, "Cell diameter",  HOFFSET(GridParams, dx), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Face Area", HOFFSET(GridParams, dA), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Cell Vol. Area", HOFFSET(GridParams, dV), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type Of B.C.",   HOFFSET(GridParams, typeOfBC), myh5t_String32);
  err+= H5Tinsert(myh5t_GP, "Total Point No.",HOFFSET(GridParams, Nbc),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Space OOA",      HOFFSET(GridParams, spOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Time  OOA",      HOFFSET(GridParams, tmOOA),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "ArtVisc Flag",   HOFFSET(GridParams, artviscosity),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Gas Cnst Gamma", HOFFSET(GridParams, gamma), H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "CFL parameter",  HOFFSET(GridParams, CFL),   H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Visc parameter", HOFFSET(GridParams, etav),  H5T_NATIVE_DOUBLE);
  err+= H5Tinsert(myh5t_GP, "Type of O/P",    HOFFSET(GridParams, typeOfOP),  H5T_NATIVE_INT);
  err+= H5Tinsert(myh5t_GP, "Output File",    HOFFSET(GridParams, outFileBase), myh5t_String128);
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
  cout <<"\tReading in Parameters...\n";
  hid_t indataset = H5Dopen(myh5_InFile, "/parameters/params");  
  myh5_err = H5Dread(indataset, myh5t_GP, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gp);
  myh5error(myh5_err,"reading parameters dataset to memory");
  myh5_err = H5Dclose(indataset);
  // Read units to grid::uc
/*  cout <<"\tReading in Units...";
  indataset = H5Dopen(myh5_InFile, "/units/units");
  myh5_err = H5Dread(indataset, myh5t_UC, H5S_ALL, H5S_ALL, H5P_DEFAULT, &uc);
  myh5error(myh5_err,"reading units dataset to memory");
  myh5_err = H5Dclose(indataset);
  cout <<" Done.\n";
*/
  myh5_err = H5Tclose(myh5t_GP);    myh5error(myh5_err,"Closing Param dtype");
  myh5_err = H5Tclose(myh5t_UC);    myh5error(myh5_err,"Closing Unit  dtype");
  myh5_err = H5Tclose(myh5t_String32);   myh5error(myh5_err,"Closing str32 dtype");
  myh5_err = H5Tclose(myh5t_String128);  myh5error(myh5_err,"Closing st128 dtype");
  myh5_err = H5Tclose(myh5t_intNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Tclose(myh5t_dblNDVec);   myh5error(myh5_err,"Closing NDvec dtype");
  myh5_err = H5Fclose(myh5_InFile);
  
  // Some Diagnostics follow.
//  cout<<"Density units: "<<uc.density<<" "<<uc.rhoVal<<endl;
//  cout<<"Simtime: "<<gp.simtime<<endl;
//  cout<<"Starttime: "<<gp.starttime<<endl;
//  cout<<"BCS: "<<gp.typeOfBC<<endl;
//  cout<<"Outputfilebase: "<<gp.outFileBase<<endl;
  
  // Copy GridParams data to some code variables:
/*  baseGrid::typeofic = "ICFILE"; // Just set this.
  baseGrid::typeofbc = gp.typeOfBC;
  baseGrid::outpath  = gp.outFileBase;
  baseGrid::spOOA = static_cast<enum ooaSpace>(gp.spOOA);
  baseGrid::tmOOA = static_cast<enum ooaTime>(gp.tmOOA);
  baseGrid::dt = 0.;
  baseGrid::maxtime = false;
*/
  return(0);

}


int tempgrid::assignDataFromHDF5(string infile)
{
  // Opens the HDF5 file and gets the data from it.
  // The grid parameters must be already read, and the grid set up to hold the data.

  int nel = gp.Ncell;
  double *dblarr;
  int    *intarr;
  dblarr = new double [nel]; if (!dblarr) {cerr<<"memory assignment dblarr.\n";exit(1);}
  intarr = new int [nel]; if (!intarr) {cerr<<"memory assignment intarr.\n";exit(1);}
  //  hsize_t dim[1];
  //  dim[0] = nel;
  hsize_t dim[GNDIM];
  for (int i=0;i<GNDIM;i++) {
    dim[i] = gp.NG[i];
  } // This is for N-D data.

  //cout <<"Opening File...";
  hid_t myh5_InFile = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  myh5error(myh5_InFile,"opening infile");
  //cout <<"\tOpening Data Group...";
  hid_t myh5_Group = H5Gopen(myh5_InFile, "/data");
  myh5error(myh5_Group,"opening grid dataset");
  //cout <<"\tCreating 1D Dataspace...";
  //  hid_t myh5_Dspace = H5Screate_simple(1, dim, 0);
  hid_t myh5_Dspace = H5Screate_simple(GNDIM, dim, 0);
  myh5error(myh5_Dspace,"creating dataspace");
  //cout <<"\t Creating Dataset...";
  hid_t myh5_Dset = H5Dopen(myh5_Group, "X-Position");
  myh5error(myh5_Dset,"Open data-dataset");// cout <<"Done!\n";
  //cout <<"\t Reading X-position...";
  herr_t myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
  myh5error(myh5_err,"Reading data From myh5_Dset");// cout <<"Done!\n";
  cell* cpt=firstPt(); int i=0;
  do {cpt->x[0] = dblarr[i];i++;} while ((cpt=nextPt(cpt))!=0);
  if(nel!=i) {cerr<<"Not Correct Number of points read\n";exit(1);}
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  if(GNDIM>1) {
    //cout <<"Opening Dataset y-pos.\n";
    myh5_Dset = H5Dopen(myh5_Group, "Y-Position");
    myh5error(myh5_Dset,"create data-dataset");
    cpt=firstPt(); i=0;
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"reading data to dblarr");
    do {cpt->x[1] = dblarr[i];i++;} while ((cpt=nextPt(cpt))!=0);
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  if(GNDIM>2) {  
    //cout <<"Opening Dataset z-pos.\n";
    myh5_Dset = H5Dopen(myh5_Group, "Z-Position");
    myh5error(myh5_Dset,"create data-dataset");
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"reading data to myh5_Dset");
    cpt=firstPt(); i=0;
    do {cpt->x[2] = dblarr[i];i++;} while ((cpt=nextPt(cpt))!=0);
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");    
  }
  // read ids.
  //cout <<"Opening Dataset for IDs.\n";
  myh5_Dset = H5Dopen(myh5_Group, "ID");
  myh5error(myh5_Dset,"create data-dataset");
  //cout<<" Reading IDs...";
  myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intarr);
  myh5error(myh5_err,"writing data to myh5_Dset");
  //cout <<"Read... assigning to grid...";
  cpt=firstPt(); i=0;
  do {cpt->id = intarr[i]; i++;} while ((cpt=nextPt(cpt))!=0);
  myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");
  //cout <<"\tDone.\n";
  // read State Variables.
  string s[8];
  s[RO]="Density"; s[PG]="Gas-Pressure"; s[VX]="X-Velocity";
  s[VY]="Y-Velocity"; s[VZ]="Z-Velocity"; s[BX]="Bx"; s[BY]="By"; s[BZ]="Bz";
  
  for (int vct=0;vct<NVAR;vct++) {// Loop 5 or 8 times, for HD and MHD.  
    //cout <<"Opening Dataset for "<<s[vct]<<"... ";
    myh5_Dset = H5Dopen(myh5_Group, (s[vct]).c_str());
    myh5error(myh5_Dset,"create data-dataset");
    //cout <<"Reading Data... ";
    myh5_err = H5Dread(myh5_Dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dblarr);
    myh5error(myh5_err,"reading data to dblarr");
    myh5_err = H5Dclose(myh5_Dset);   myh5error(myh5_err,"Closing myh5_Dset");
    //cout <<"Read, assigning to grid... ";
    // Assign data to grid: Assumes grid is set up for same point ordering as was written.
    cpt=firstPt(); i=0; enum primitive var = static_cast<primitive>(vct);
    do {
      cpt->P[var] = dblarr[i];
      i++;
    } while ((cpt=nextPt(cpt))!=0);
    //    cout <<"i = "<<i<<"\t Ncell = "<<gp.Ncell<<endl;
    if (i!=gp.Ncell) {cerr<<"Error: too many gridpoints in output!\n";exit(100);}
    //cout <<"\tDone.\n";
  }
  myh5_err = H5Sclose(myh5_Dspace); myh5error(myh5_err,"Closing myh5_Dspace dataspace");
  myh5_err = H5Gclose(myh5_Group);  myh5error(myh5_err,"Closing myh5_Group");
  myh5_err = H5Fclose(myh5_InFile);  myh5error(myh5_err,"Closing myh5_InFile file");
  delete [] dblarr;
  delete [] intarr;
  
//  cpt = firstPt();
//  do {for(int v=0;v<NVAR;v++) cpt->Ph[v]=cpt->P[v];} while ((cpt=nextPt(cpt))!=0);
  return(0);
}
