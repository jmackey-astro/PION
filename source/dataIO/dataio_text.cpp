/// \file dataio_text.cpp
/// \author Jonathan Mackey
///
/// Class definitions for ASCII Text data I/O.
/// 
/// modified:\n
/// - 2018.05.01 JM: moved dataio_text from dataio.cpp.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "dataIO/dataio_text.h"
#include "dataIO/readparams.h"
#include "grid/stellar_wind_BC.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;



//----------------------------------------------
//------------------ TEXT I/O ------------------
//----------------------------------------------

dataio_text::dataio_text(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
: DataIOBase(SimPM)
{
  rp = 0;
  gp=0;
  eqn=0;
}



// ##################################################################
// ##################################################################


dataio_text::~dataio_text()
{
  if (rp) {
    delete rp; rp=0;
  }
  gp=0;
  //if (eqn) {
  //  delete eqn; eqn=0;
  //}
}



// ##################################################################
// ##################################################################


int dataio_text::ReadHeader(
      string pfile,           ///< Name of parameter file
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  // read info from parameterfile.
  cout <<"dataio_text::ReadHeader() Read simulation info from parameterfile.\n";
  if (!file_exists(pfile)) {
    cerr <<"dataio_text::ReadHeader() file doesn't exist.\n";
    return 99;
  }
  int err=get_parameters(pfile, SimPM);
  cout <<"dataio_text::ReadHeader() Header info all read in.\n";
  return err;
}



// ##################################################################
// ##################################################################


int dataio_text::ReadData(
      string,   ///< Name of file
      vector<class GridBaseClass *> &cg, ///< address of vector of grid pointers.
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (!cg[0])
    rep.error("dataio_text::ReadData() null pointer to grid!",cg[0]);
  dataio_text::gp = cg[0];

  cout <<"dataio_text::ReadData() Assigning initial data.\n";
  int err = assign_initial_data(SimPM);
  cout <<"dataio_text::ReadData() Assigned initial data.\n";

  // We have read in all the parameters we need, so delete it.
  //if(rp!=0) {
  //  cout <<"deleting readparameters class.\n";
  //  delete rp; rp=0;
  //} // Delete the ReadParams class.

  return err;
}



// ##################################################################
// ##################################################################


void dataio_text::SetSolver(FV_solver_base *solver)
{
  cout <<"dataio_text::SetSolver() Setting solver pointer.\n";
  dataio_text::eqn = solver;
  return;
}




// ##################################################################
// ##################################################################


int dataio_text::OutputData(
      const string outfile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      const long int counter   ///< number to stamp file with (e.g. timestep)
      )
{
  if (!cg[0])
    rep.error("dataio_text::output_ascii_data() null pointer to grid!",cg[0]);

  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the nested grid.
    CI.set_dx(SimPM.nest_levels[l].dx);
    ostringstream temp; temp << outfile << "_level";
    temp.width(2); temp.fill('0');
    temp << l;
    dataio_text::gp = cg[l];

    cout <<"dataio_text::OutputData() writing data.\n";
    string fname = set_filename(temp.str(), counter);
    int err = output_ascii_data(fname, SimPM);
    cout <<"dataio_text::OutputData() written data.\n";

  }
  return err;
}



// ##################################################################
// ##################################################################


///
/// set filename based on counter, outfile-base-name.
///
std::string dataio_text::set_filename(
      const std::string outfile,
      const long int    counter
      )
{
  ostringstream temp; temp.str("");
  temp << outfile <<".";
  //
  // only add timestep if counter >=0 (initial conditions get no timestep).
  //
  if (counter >= 0) {
    temp.width(Ndigits); temp.fill('0');
    temp << counter <<".";
  }
  temp <<"txt";
  return temp.str();
}



// ##################################################################
// ##################################################################



int dataio_text::get_parameters(
      string pfile,  ///< Name of parameter file
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  cout <<"(dataio_text::get_parameters) from file "<<pfile<<" starting.\n";

  rp = new ReadParams();
  if (!rp) {cerr<<"(init) Error initialising *rp...exiting."<<"\n";return(1);}

  // Parse the parameterfile.
  if(rp->read_paramfile(pfile) !=0)
  {cerr<<"(init) Error reading parameterfile...exiting."<<"\n";return(1);}
  //   rp->write_out_parameters();
  
  string ts;
  // Basic Grid Properties
  
  ts = rp->find_parameter("solver");
  if      (ts=="") {rep.warning("No solver specified; using LF.","LF",ts); SimPM.solverType = FLUX_LF;}
  else if (ts=="LF") {SimPM.solverType = FLUX_LF;}
  else if (ts=="RSlinear") {SimPM.solverType = FLUX_RSlinear;}
  else if (ts=="RSexact" ) {SimPM.solverType = FLUX_RSexact; }
  else if (ts=="RShybrid") {SimPM.solverType = FLUX_RShybrid;}
  else if (ts=="RSroe"   ) {SimPM.solverType = FLUX_RSroe;}
  else if (ts=="RSroe_pv") {SimPM.solverType = FLUX_RSroe_pv;}
  else if (ts=="RS_FVS"  ) {SimPM.solverType = FLUX_FVS;}
   else rep.error("No solver specified!",ts);
  
  /** \section eqnndim
   * If I ever use a set of equations that doesn't have 3D vectors like velocity and 
   * magnetic field, then I will need to set eqnndim more cleverly, but now I always 
   * set it to 3.
   * */
  ts = rp->find_parameter("eqnndim");
  if (ts=="") {
    rep.warning("Setting eqnndim=3, not found in param-file",-1,-1);
    SimPM.eqnNDim = 3;
  }
  else if (ts=="3") {SimPM.eqnNDim = 3;}
  else rep.error("eqnndim !=3 in paramfile.",ts);
  
  ts = rp->find_parameter("eqn");
  if      (ts=="") {
    rep.warning("No equations specified; using Euler.","euler",ts);
    SimPM.eqntype = 1; SimPM.nvar = 5;
  }
  else if (ts=="euler")   {SimPM.eqntype = EQEUL;        SimPM.nvar = 2+SimPM.eqnNDim;}
  else if (ts=="isohydro"){SimPM.eqntype = EQEUL_ISO;    SimPM.nvar = 2+SimPM.eqnNDim;}
  else if (ts=="i-mhd")   {SimPM.eqntype = EQMHD; SimPM.nvar = 2+2*SimPM.eqnNDim;}
  else if (ts=="glm-mhd") {SimPM.eqntype = EQGLM; SimPM.nvar = 3+2*SimPM.eqnNDim;}
  else if (ts=="fcd-mhd") {SimPM.eqntype = EQFCD; SimPM.nvar = 2+2*SimPM.eqnNDim;}
#ifdef INCLUDE_EINT_ADI_HYDRO
  else if (ts=="euler-Eint") {
    SimPM.eqntype = EQEUL_EINT;
    SimPM.nvar = 3+SimPM.eqnNDim;
  }
#endif // if INCLUDE_EINT_ADI_HYDRO
  else {
    cout <<" Please use \'euler\',\'i-mhd\',\'glm-mhd\',\'euler-Eint\' for eqn-type\n";
    rep.error("\'eqn\' not in param-file: Don't know what equations to set",ts);
  }

  ts=rp->find_parameter("gridtype"); 
  if(ts!="")
    SimPM.gridType =atoi(ts.c_str());
  else {rep.warning("No GridType set; using Uniform Finite-Volume.",1,SimPM.gridType); SimPM.gridType=1;}
  if (SimPM.gridType!=1) rep.error("GridType not set to Uniform Finite-Volume; need to set later.",SimPM.gridType);
  
  ts = rp->find_parameter("ndim");
  if (ts!="")  SimPM.ndim = atoi(ts.c_str());
  //  else rep.error("No grid dimensions specified in parameterfile",ts);
  else {
    cout <<"No ndim specified in parameterfile, assuming ndim=1!!!\n";
    SimPM.ndim=1;
  }
  
  ts = rp->find_parameter("ntracer");
  int tr=0;
  if (ts!="") tr = atoi(ts.c_str());
  else {
    cout<<"No tracer number specified in paramfile, using none.\n";
    tr = 0;
  }
  SimPM.ntracer = tr;
  SimPM.ftr = SimPM.nvar; // Set first tracer element in state vector to be after all primitive vars.
  
  int nv=0;
  ts = rp->find_parameter("nvar");
  if (ts!="") nv = atoi(ts.c_str());
  else {
    cout <<"nvar not specified in paramfile. Inferring from equations.\n";
    nv = SimPM.nvar;
  }
  if (nv < SimPM.nvar+SimPM.ntracer) rep.error("Bad nvar in paramfile (account for tracers).",nv);
  else if (nv>SimPM.nvar+SimPM.ntracer) rep.error("too many variables (which are tracers?)",nv);
  if (SimPM.nvar <= nv) {
    SimPM.nvar = nv;
  }
  else {
    rep.warning("nvar in paramfile less than that specified by eqn. type",SimPM.nvar,nv);
  }
  
  ts = rp->find_parameter("coordinates");
  if (ts=="") SimPM.coord_sys = COORD_CRT;
  else if (ts=="cartesian")   SimPM.coord_sys = COORD_CRT;
  else if (ts=="cylindrical") SimPM.coord_sys = COORD_CYL;
  else if (ts=="spherical")   SimPM.coord_sys = COORD_SPH;
  else rep.error("Don't recognise coordinate system in GetParameters",ts);
  
  // Now assign dataio_text variables with parameters from the file.
  // Get the base path/filename to write output to.
  string outpath = rp->find_parameter("OutputPath");
  string outfile = rp->find_parameter("OutputTextFile");
  ostringstream temp;
  temp << outpath << outfile;
  SimPM.outFileBase = temp.str();
  string oft = rp->find_parameter("OutputFileType");
  if      (oft=="TXT"  || oft=="txt" ||
     oft=="TEXT" || oft=="text") SimPM.typeofop =1;
#ifdef FITS
  else if (oft=="FITS" || oft=="fits") SimPM.typeofop =2;
  else if (oft=="BOTH" || oft=="both") SimPM.typeofop =4;
#endif // if FITS
#ifdef SILO
  else if (oft=="SILO" || oft=="silo") SimPM.typeofop =5;
#endif // if SILO
  else    {cerr<<"Error, bad type of outputfile in pfile.\n"; return(1);}
  SimPM.opfreq = atoi( (rp->find_parameter("OutputFrequency")).c_str() );
  cout <<"\tOutFile: "<<SimPM.outFileBase<< ".xxx\t Type="<<SimPM.typeofop;
  cout <<" every "<<SimPM.opfreq<<" timesteps."<<"\n";
  
  //
  // Boundary conditions.
  // First do the edges of the domain:
  //
  SimPM.BC_XN = rp->find_parameter("BC_XN");
  SimPM.BC_XP = rp->find_parameter("BC_XP");
  if (SimPM.ndim>1) {
    SimPM.BC_YN = rp->find_parameter("BC_YN");
    SimPM.BC_YP = rp->find_parameter("BC_YP");
  }
  else {
    SimPM.BC_YN = "NONE";
    SimPM.BC_YP = "NONE";
  }
  if (SimPM.ndim>2) {
    SimPM.BC_ZN = rp->find_parameter("BC_ZN");
    SimPM.BC_ZP = rp->find_parameter("BC_ZP");
  }
  else {
    SimPM.BC_ZN = "NONE";
    SimPM.BC_ZP = "NONE";
  }
  //
  // Now do the internal boundaries (if any).  Seek the string
  // for a given internal boundary, and add to a vector until
  // no more are found.
  //
  ts = rp->find_parameter("BC_Ninternal"); 
  if (ts=="")  SimPM.BC_Nint = 0;
  else         SimPM.BC_Nint = atoi(ts.c_str());
  if (SimPM.BC_Nint>0) {
    if (SimPM.BC_INT) rep.error("BC_INT already initialized",2);
    SimPM.BC_INT = mem.myalloc(SimPM.BC_INT,SimPM.BC_Nint);
  }
  //SimPM.BC_INT.clear();
  int v=0;
  do {
    ostringstream intbc; intbc.str("");
    intbc << "BC_INTERNAL_";
    intbc.width(3); intbc.fill('0');
    intbc << v;
    string temp = rp->find_parameter(intbc.str());
    if (temp != "") {
      SimPM.BC_INT[v] = temp;
      v++;
    }
    else
      rep.error("Didn't find internal BC",v);
  } while (v<SimPM.BC_Nint);
  //SimPM.BC_Nint = SimPM.BC_INT.size();

  // Number of boundary cells
  // (set automatically based on order of scheme)
  SimPM.Nbc = -1;

  
  // Physics:
  SimPM.gamma = atof( (rp->find_parameter("GAMMA")).c_str());
  SimPM.CFL   = atof( (rp->find_parameter("CFL")).c_str());
  if( (SimPM.artviscosity=atoi( (rp->find_parameter("ArtificialViscosity")).c_str())) ==0) {
    //    cout <<"\tNot using Artificial Viscosity.\n";
    SimPM.etav=0.;
  }
  else if ( SimPM.artviscosity >=0 && SimPM.artviscosity <= 4) {
    //    cout <<"\tUsing Artificial Viscosity number ";
    //cout <<SimPM.artviscosity<<" with eta = ";
    SimPM.etav = atof( (rp->find_parameter("EtaViscosity")).c_str() );
    //    cout <<SimPM.etav<<"\n";
  }
  else {cerr<<"\tUnknown viscosity requested... please update me.\n"; return(1);}
  
  // Overall Grid Properties
  SimPM.NG[0] = atoi( (rp->find_parameter("NGridX")).c_str());
  SimPM.Ncell = SimPM.NG[0];
  if (SimPM.ndim>1) {
    SimPM.NG[1] = atoi( (rp->find_parameter("NGridY")).c_str());
    SimPM.Ncell*= SimPM.NG[1];
  }
  if (SimPM.ndim>2) {
    SimPM.NG[2] = atoi( (rp->find_parameter("NGridZ")).c_str());
    SimPM.Ncell*= SimPM.NG[2];
  }

  // Grid Point properties.
//  cout <<"\t Getting (xmin,xmax) from Pfile, setting range.\n";
  SimPM.Xmin[0]   = atof( (rp->find_parameter("Xmin")).c_str());
  SimPM.Xmax[0]   = atof( (rp->find_parameter("Xmax")).c_str());
  SimPM.Range[0]  = SimPM.Xmax[0]-SimPM.Xmin[0];
  if (SimPM.ndim>1) {
    SimPM.Xmin[1]  = atof( (rp->find_parameter("Ymin")).c_str());
    SimPM.Xmax[1]  = atof( (rp->find_parameter("Ymax")).c_str());
    SimPM.Range[1] = SimPM.Xmax[1]-SimPM.Xmin[1];
  }
  if (SimPM.ndim>2) {
  SimPM.Xmin[2]  = atof( (rp->find_parameter("Zmin")).c_str());
  SimPM.Xmax[2]  = atof( (rp->find_parameter("Zmax")).c_str());
  SimPM.Range[2] = SimPM.Xmax[2]-SimPM.Xmin[2];
  }
  
  // Order of Accuracy of integration:
  SimPM.spOOA = atoi( (rp->find_parameter("OrderOfAccSpace")).c_str());
  SimPM.tmOOA = atoi( (rp->find_parameter("OrderOfAccTime")).c_str());
  if (SimPM.spOOA==OA1) {}// cout<<"\tFirst Order Spatial Accuracy.\n";
  else if (SimPM.spOOA==OA2) {}// cout<<"\tSecond Order Spatial Accuracy.\n";
  else {cout<<"\t Error: Spatial Accuracy requested is not available.\n";return(1);}
  if (SimPM.tmOOA==OA1) {}// cout<<"\tFirst Order Time Accuracy.\n";
  else if (SimPM.tmOOA==OA2) {} //cout<<"\tSecond Order Time Accuracy.\n";
  else {cout<<"\t Error: Time Accuracy requested is not available.\n";return(1);}

  // Timing Quantities.
  SimPM.starttime = atof( (rp->find_parameter("StartTime")).c_str());
  SimPM.simtime = SimPM.starttime;
  if ( (rp->find_parameter("FinishTime")) =="") SimPM.finishtime=-1.;
  else SimPM.finishtime = atof( (rp->find_parameter("FinishTime")).c_str());
  SimPM.dt = 0.;
  SimPM.timestep =0; // Counter for what timestep we are on.

#ifdef TESTING
  dp.initERG = dp.initMMX = dp.initMMY = dp.initMMZ = 0.;
  dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange = 0.;
#endif //TESTING
  SimPM.maxtime = false;

  // Ref state vec for riemann solver.
  for (int v=0;v<MAX_NVAR;v++) SimPM.RefVec[v] = 1.;
  
  // Code units
  string unit;
  if( (unit = rp->find_parameter("units")) == "") {
    cout <<"No units specified in parameterfile.  Using code units for input data.\n";
    uc.unitsys  = "code";
    uc.density  = "no units";
    uc.length   = "no units";
    uc.velocity = "no units";
    uc.bfield   = "no units";
    uc.rhoVal = uc.lenVal = uc.velVal = uc.magVal = 1.;
  }
  else if (unit=="MKS") {
    cout <<"\t Using MKS (SI) as reference units.\n";
    uc.unitsys  = "mks";
    uc.density  = "kg/m^3";
    uc.length   = "m";
    uc.velocity = "m/s";
    uc.bfield   = "Tesla";
    uc.rhoVal = atoi( (rp->find_parameter("rhoval")).c_str());
    uc.lenVal = atoi( (rp->find_parameter("lenval")).c_str());
    uc.velVal = atoi( (rp->find_parameter("velval")).c_str());
    uc.magVal = atoi( (rp->find_parameter("magval")).c_str());
  }
  else {rep.error("Don't recognise units system",unit);}
  
  cout <<"(dataio_text::get_parameters) Finished getting parameters.\n";
  return(0);
} // dataio_text::get_parameters




// ##################################################################
// ##################################################################




int dataio_text::output_ascii_data(
      string outfile,
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  ofstream outf(outfile.c_str());
  if(!outf.is_open()) 
  {
    cerr << "Error opening file " << outfile << " for writing.  Quitting..." <<"\n";
    return(1);
  }
  //  cout <<"(dataio_text::output_ascii_data) Writing data in format: x[Ndim], rho, p_g, v_x, v_y, v_z, e_int(erg/mass), [B_x, B_y, B_z, p_g+p_m].\n";
  double b2=0.; // magnetic field squared.
  double dx = gp->DX();
#ifdef RT_TESTING_OUTPUTCOL
  double Utemp[SimPM.nvar];
#endif // RT_TESTING_OUTPUTCOL
  //  outf.setf( ios_base::fixed,ios_base::floatfield );
  //  outf.precision(6);
  outf << "# format: x,[y,z,],rho,pg,vx,vy,vz,[Bx,By,Bz],[Tr0,Tr1,Tr2,..],T,[Tau0,Tau1,...]\n";
  outf << "# time = "<<SimPM.simtime<<"  timestep = "<<SimPM.timestep<<"\n";
  outf.setf( ios_base::scientific );
  outf.precision(14);

  // list of B-field vector elements for calculating div(B)
  int vars[3]; vars[0] = static_cast<int>(BX);
  vars[1] = static_cast<int>(BY); vars[2] = static_cast<int>(BZ);

  // Go through every point, output one line per point.
  class cell *cpt=gp->FirstPt(); do {
     if(CI.get_dpos(cpt,0)<gp->SIM_Xmin(XX)+dx) outf <<"\n"; // put in a blank line for gnuplot
     // First positions.
     outf << CI.get_dpos(cpt,0) << "  ";
     if (SimPM.ndim>1) outf << CI.get_dpos(cpt,1) << "  ";
     if (SimPM.ndim>2) outf << CI.get_dpos(cpt,2) << "\t";
     // Next all primitive variables.
     for (int v=0;v<SimPM.nvar;v++) outf <<cpt->P[v]<<"  ";

      // internal energy/ temperature.
    if      (MP) outf << MP->Temperature(cpt->P,SimPM.gamma);
    else if (eqn) {
      outf << eqn->eint(cpt->P,SimPM.gamma);
      // total energy, x-momentum
      //eqn->PtoU(cpt->P,Utemp,SimPM.gamma);
      //outf <<"  "<< Utemp[ERG] <<"  "<< Utemp[MMX];
    }
    // mhd vars.
    if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQGLM || SimPM.eqntype==EQFCD) {
      b2 = cpt->P[BX]*cpt->P[BX] +cpt->P[BY]*cpt->P[BY] +cpt->P[BZ]*cpt->P[BZ];
      //       outf <<"  "<< cpt->P[BX] <<"  "<< cpt->P[BY] <<"  "<< cpt->P[BZ] <<"  ";
      outf <<"  "<< cpt->P[PG]+b2/2.;
      outf <<"  "<< eqn->Divergence(cpt,0,vars,gp);
    }
#ifdef RT_TESTING_OUTPUTCOL
    for (int v=0;v<SimPM.RS.Nsources;v++) {
      //cout <<"hello";
      CI.get_col(cpt, v, Utemp);
      for (int iT=0; iT<SimPM.RS.sources[v].NTau; iT++)
        outf <<"  "<< Utemp[iT];
    }
#endif // RT_TESTING_OUTPUTCOL
    outf  <<"\n";
  } while ( (cpt=gp->NextPt(cpt))!=0);
  
  outf.setf(ios_base::fmtflags(0),ios_base::floatfield);
  outf.precision(6);
  outf.close();
  return(0);
}



// ##################################################################
// ##################################################################



int dataio_text::assign_initial_data(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  cout <<"(dataio_text::assign_initial_data) Assigning Data.\n";
  int err=0;
  double dx = gp->DX();

  // Get initial conditions, and assign them to the grid.
  string typeofic = rp->find_parameter("IC");
  if(typeofic=="SHOCKTUBE") {
    cout << "\t Using shock tube conditions.\n";
    int which_riemann = atoi((rp->find_parameter("RIEMANN")).c_str());
    double left[SimPM.nvar], right[SimPM.nvar];
    double interface=0.0;
    err += get_riemann_ics(SimPM, which_riemann, left, right, &interface);
    class cell *cpt=gp->FirstPt();
    if (SimPM.ndim==1) {
      //interface = 0.5;
       do {
   if (CI.get_dpos(cpt,0)<=interface) for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = left[v];
   else for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = right[v];
       } while( (cpt=gp->NextPt(cpt))!=0);
    } // if 1D

    else if (SimPM.ndim==2) {
      // Get Shock Angle.
      double theta;
      if ( (rp->find_parameter("ShockAngle"))=="" ) {
  rep.warning("dataio_text::assign_initial_data 2D sim, but no shock angle in paramfile.",1,1);
  theta=M_PI/4.;
      }
      else {
  theta = atof((rp->find_parameter("ShockAngle")).c_str()) *M_PI/180.;
      }
      cout <<"theta="<<theta<<"\n";
      // Assign position of dividing line: pivot point and slope.
      double xpivot,ypivot, slope;
      if (theta >0. && theta<M_PI/2.) {
  xpivot = gp->SIM_Xmin(XX)+ gp->SIM_Range(XX)/2.;
  ypivot = gp->SIM_Xmin(YY)+ gp->SIM_Range(YY)/2.;
  slope = tan(theta);
  //if (!eqn) rep.error("equations not set!",eqn);
  //eqn->rotateXY(left,-(M_PI/2.-theta)); eqn->rotateXY(right,-(M_PI/2.-theta));
  double ct=cos(-(M_PI/2.-theta)); double st=sin(-(M_PI/2.-theta));
  double vx, vy;
  vx = left[VX]*ct - left[VY]*st;
  vy = left[VX]*st + left[VY]*ct;
  left[VX] = vx; left[VY]=vy;
  vx = right[VX]*ct - right[VY]*st;
  vy = right[VX]*st + right[VY]*ct;
  right[VX] = vx; right[VY]=vy;
  if (SimPM.eqntype==EQMHD || SimPM.eqntype==EQGLM ||
      SimPM.eqntype==EQFCD) {
    vx = left[BX]*ct - left[BY]*st;
    vy = left[BX]*st + left[BY]*ct;
    left[BX] = vx; left[BY] = vy;
    vx = right[BX]*ct - right[BY]*st;
    vy = right[BX]*st + right[BY]*ct;
    right[BX] = vx; right[BY] = vy;
  }

      }
      else if (theta<=0.) {
  cout <<"using theta=0\n";
  xpivot = gp->SIM_Xmin(XX)+ gp->SIM_Range(XX)/2.;
  ypivot = gp->SIM_Xmin(YY)+ gp->SIM_Range(YY)/2.;
  slope = 0;
  theta = 0;
      }
      else {
  cout <<"please use an angle between 0 and 90 degrees.\n";
  exit(1);
      }
      
      if (fabs(theta)<1.e-6) { // If zero slope, then just run a vertical shock at x=interface.
  do {
    if (CI.get_dpos(cpt,XX) < interface) 
      for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = left[v];
    else for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = right[v];
  } while( (cpt=gp->NextPt(cpt))!=0);
      }
      else {
  cout <<"2D and non-zero slope, so setting angle to "<<theta<<" radians.\n";
  // Slope is non-zero, so set a shock at angle theta to x-axis,
  // passing through centre of domain.
  do {
    if      (CI.get_dpos(cpt,YY)-dx > ypivot +slope*(CI.get_dpos(cpt,XX)+dx/2.-xpivot))
      for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = left[v];
    else if (CI.get_dpos(cpt,YY)+dx <=ypivot +slope*(CI.get_dpos(cpt,XX)-dx/2.-xpivot))
      for (int v=0;v<SimPM.nvar;v++) cpt->P[v] = cpt->Ph[v] = right[v];
    else {
      //      cout <<"intermediate point. x="<<CI.get_dpos(cpt,XX)<<", y="<<CI.get_dpos(cpt,YY)<<",\t";
      // intermediate region, so we have to take an average.
      // Do this by splitting cell into 32x32 subcells.
      int nint=32; double dv = 1.;
      for (int v=0;v<SimPM.ndim;v++) {dv*= dx/nint;}
      double dxc = dx/nint;
      double startpt[SimPM.ndim]; double pos[SimPM.ndim]; double frac=0.;
      // Set first position.
      int ntot=1;
      for (int v=0;v<SimPM.ndim;v++) {
        startpt[v] = CI.get_dpos(cpt,v)-dx/2.+dxc/2.;
        ntot *=nint;
      }
      for (int i=0;i<ntot;i++) {
        pos[XX] = startpt[XX] + (i%nint)*dxc;
        if (SimPM.ndim>1) pos[YY] = startpt[YY] + ((i/nint)%nint)*dxc;
        if (SimPM.ndim>2) pos[ZZ] = startpt[ZZ] + ((i/nint/nint)%nint)*dxc;
        if (pos[YY] > ypivot+slope*(pos[XX]+dxc/2.-xpivot)) {
    frac +=dv;
        }
      }
      frac /= gp->DV();
      //cout <<"frac = "<<frac<<"\n";
      for (int v=0;v<SimPM.nvar;v++)
        cpt->P[v] = cpt->Ph[v] = frac*left[v] + (1.-frac)*right[v];
    }
  } while( (cpt=gp->NextPt(cpt))!=0);
      } // if slope is non-zero.
    } // if 2d
    if (SimPM.eqntype == EQGLM) {
      cpt=gp->FirstPt(); do {cpt->P[SI] = cpt->Ph[SI] = 0.;} while( (cpt=gp->NextPt(cpt))!=0);
    }
    if (err!=0) {cerr<<"\tError assigning data"<<"\n";return(1);}
  }
  else {
    cerr << "\tWhat sort of ICs??? only know 'SHOCKTUBE'"<<"\n";
    return(1);
  }
  
  // Initial Condition Smoothing?
  int smooth = atoi((rp->find_parameter("SmoothICS")).c_str());
  if (smooth>0) {
    cout <<"Can't smooth data... sorry!  Write a smoothing function.\n";
//    cout <<"\t Smoothing Data... ns="<<smooth;
//    err += gp->smoothData(smooth);
//    cout <<"  ...Done\n";
    if (err!=0) {cerr<<"\tSmoothing Data failed. exiting.\n"; return(1);}
  } // if smooth>0;
  
  // Initial Conditions:  Add noise?
  int noise = atoi((rp->find_parameter("NoisyICS")).c_str());
  if(noise>0) {
    cout <<"\t Adding noise to Data at 0.1% level to pressure, type="<<noise;
    err += add_noise2data(SimPM, noise, 0.001);
    cout <<"  ...Done\n";
    if (err!=0) {cerr<<"\tAdding noise to Data failed. exiting.\n"; return(1);}
  } // if noise>0;
    
  //  initial_conserved_quantities();
  cout <<"(dataio_text::assign_initial_data) Done.\n";
  return(0);
}



// ##################################################################
// ##################################################################



int dataio_text::get_riemann_ics(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      int sw, double *l, double *r, double *xm)
{
  // These are Toro's five tests on p.225 of his book.
  switch (sw) {
  case 1:
    /** case 1: Toro's test no.1 on p.225 of his book.\n*/
    l[RO]=1.;    l[PG]=1.;  l[VX]=0.75; l[VY]=l[VZ]=0.;
    //l[RO]=1.;    l[PG]=1.;  l[VX]=0.0; l[VY]=l[VZ]=0.;
    r[RO]=0.125; r[PG]=0.1; r[VX]=0.0;  r[VY]=r[VZ]=0.;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.3;
    SimPM.gamma=1.4;
    break;
   case 2:
    /** case 2: Toro's test no.2 on p.225 of his book.\n*/
     // ************************************************
     // Temporarily taken over for isothermal testing!!!
     // ************************************************
    l[RO]= 10.0;    l[PG]=1.0;  l[VX]= 0.0; l[VY]=l[VZ]=0.;
    r[RO]= 1.0;    r[PG]=1.0;  r[VX]= 0.0;  r[VY]=r[VZ]=0.;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.5;
    SimPM.gamma=1.0;
    //l[RO]=1.;    l[PG]=0.4;  l[VX]=-2.0; l[VY]=l[VZ]=0.;
    //r[RO]=1.;    r[PG]=0.4;  r[VX]=2.0;  r[VY]=r[VZ]=0.;
    //if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
    //  l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    //}
    //*xm = 0.5;
    //SimPM.gamma=1.4;
    break;
   case 3:
    /** case 3: Toro's test no.3 on p.225 of his book.\n*/
    l[RO]=1.; l[PG]=1000.;  l[VX]=0.0; l[VY]=l[VZ]=0.;
    r[RO]=1.; r[PG]=0.01;   r[VX]=0.0;  r[VY]=r[VZ]=0.;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM){
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.5;
    SimPM.gamma=1.4;
    break;
   case 4:
    /** case 4: Toro's test no.4 on p.225 of his book.\n*/
    l[RO]=5.99924;   l[PG]=460.894;  l[VX]=19.5975; l[VY]=l[VZ]=0.;
    r[RO]=5.99242;   r[PG]=46.0950;  r[VX]=-6.19633;  r[VY]=r[VZ]=0.;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.4;
    SimPM.gamma=1.4;
    break;
   case 5:
    /** case 5: Toro's test no.5 on p.225 of his book.\n*/
    l[RO]=1.;   l[PG]=1000.;   l[VX]=-19.59745; l[VY]=l[VZ]=0.;
    r[RO]=1.;   r[PG]=0.01;    r[VX]=-19.59745;  r[VY]=r[VZ]=0.;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.8;
    SimPM.gamma=1.4;
    break;
  // other wierd cases
  //Left                     [rho,v,p] = [0.604543, 1.876, 1.69426 ]
  //Right                    [rho,v,p] = [1, 2, 1 ]
  // This fools the linear solver...
   case 6:
    /** case 6: Slightly difficult case with an almost stationary rarefaction.*/
    l[RO]=0.604543;   l[PG]=1.69426;   l[VX]=1.876; l[VY]=l[VZ]=0.4;
    r[RO]=1;   r[PG]=1;    r[VX]=2;                 r[VY]=r[VZ]=0.5;
    if (SimPM.eqntype==2 || SimPM.eqntype==EQGLM) {
      l[BX]=l[BY]=l[BZ]=r[BX]=r[BY]=r[BZ]=0.;
    }
    *xm = 0.5;
    SimPM.gamma=1.4;
    break;
   // From here on I am using MHD test cases, so it won't work for hydro.
   case 7:
    /** case 7: Sam Falle's test 'BW', the Brio and Wu problem.\n */
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=1.; l[PG]=1.; l[VX]= l[VY]= l[VZ]=0.;
    l[BX]=0.75; l[BY]=1.; l[BZ]=0.;
    r[RO]=0.125; r[PG]=0.1; r[VX]= r[VY]= r[VZ]=0.;
    r[BX]=0.75; r[BY]=-1.; r[BZ]=0.;
    *xm = 0.5;
    SimPM.gamma=2.0;
    break;
   case 8:
    /** case 8: Sam Falle's test 'AW', an Alfven wave. (seems to have no difference in left and right states).\n */
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]= l[PG]=1.; l[VX]=0.; l[VY]= l[VZ]=1.;
    l[BX]= l[BY]=1.; l[BZ]=0.;
    r[RO]= r[PG]=1.; r[VX]=0.; r[VY]= r[VZ]=1.;
    r[BX]= r[BY]=1.; r[BZ]=0.;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    cout <<"THIS SHOCK TUBE TEST DOESN'T WORK...ALFVEN WAVE.\n"; return(1);
    break;
   case 9:
    /** case 9: Sam Falle's test 'FS', a fast shock.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=3.; l[PG]=16.333333; l[VX]=-0.732; l[VY]=-1.3333333; l[VZ]=0.;
    l[BX]=3.; l[BY]=2.309; l[BZ]=1.;
    r[RO]=1.; r[PG]=1.; r[VX]=-4.196; r[VY]=0.; r[VZ]=0.;
    r[BX]=3.; r[BY]=0.; r[BZ]=0.;
    //    l[VX]=-0.732+4.196; r[VX]=0.;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    break;
   case 10:
    /** case 10: Sam Falle's test 'SS', a slow shock.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=1.368;
    l[PG]=1.769;
    l[VX]=0.269; l[VY]=1.; l[VZ]=0.;
    l[BX]=1.; l[BY]= l[BZ]=0.;
    r[RO]=1.; 
    r[PG]=1.;
    r[VX]= r[VY]= r[VZ]=0.;
    r[BX]=1.; r[BY]=1.; r[BZ]=0.;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    break;
   case 11:
    /** case 11: Sam Falle's test 'FR', a fast rarefaction.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=1.;
    l[PG]=2.;
    l[VX]= l[VY]= l[VZ]=0.;
    l[BX]=1.; l[BY]=3.; l[BZ]=0.;
    r[RO]=0.2641; 
    r[PG]=0.2175;
    r[VX]=3.6; r[VY]=-2.551; r[VZ]=0.;
    r[BX]=1.; r[BY]= r[BZ]=0.;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    break;
   case 12:
    /** case 12: Sam Falle's test 'SR', a slow rarefaction.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=1.;
    l[PG]=2.;
    l[VX]= l[VY]= l[VZ]=0.;
    l[BX]=1.; l[BY]= l[BZ]=0.;
    r[RO]=0.2;
    r[PG]=0.1368; 
    r[VX]=1.186; r[VY]=2.967; r[VZ]=0.;
    r[BX]=1.; r[BY]=1.6405; r[BZ]=0.;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    break;
   case 13:
    /** case 13: Sam Falle's test 'OFS', an oblique fast shock.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    l[RO]=1.;
    l[PG]=1.;
    l[VX]=6.505; l[VY]=1.; l[VZ]=0.;
    l[BX]= l[BY]= l[BZ]=1.;
    r[RO]=3.;
    r[PG]=20.268;
    r[VX]=2.169; r[VY]=1.331; r[VZ]=0.331;
    r[BX]=1.; r[BY]=3.153; r[BZ]=3.153;
    *xm = 0.5;
    SimPM.gamma=5./3.;
    break;    
   case 14:
    /** case 14: Trivial case, only call it if you intend to add noise later.\n*/
    for (int v=0;v<SimPM.nvar;v++) {l[v]=r[v]=1.; *xm=0.5;}
    break;
    /** Ryu and Jones (1995) Shock Tube tests, 1a-5b follow. [Ryu \& Jones, 1995, ApJ,442,228].\n */
   case 15:
    /** case 15: Ryu and Jones test 1a.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=10.; l[VY]=l[VZ]=0.;
    l[BX]=l[BY]=5./sqrt(4*M_PI); l[BZ]=0.;
    l[PG]=20.;
    r[RO]=1.; r[VX]=-10.; r[VY]=r[VZ]=0.;
    r[BX]=r[BY]=5./sqrt(4*M_PI);r[BZ]=0.;
    r[PG]=1.;
    *xm = 0.5;
    break;
   case 16:
    /** case 16: Ryu and Jones test 1b.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=3./sqrt(4*M_PI); l[BY]=5./sqrt(4*M_PI); l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.1; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=3./sqrt(4*M_PI); r[BY]=2./sqrt(4*M_PI); r[BZ]=0.;
    r[PG]=10.;
    *xm = 0.5;
    break;
   case 17:
    /** case 17: Ryu and Jones test 2a.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.08; l[VX]=1.2; l[VY]=0.01; l[VZ]=0.5;
    l[BX]=2./sqrt(4.*M_PI); l[BY]=3.6/sqrt(4.*M_PI); l[BZ]=2./sqrt(4.*M_PI);
    l[PG]=0.95;
    r[RO]=1.; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=2./sqrt(4.*M_PI); r[BY]=4./sqrt(4.*M_PI); r[BZ]=2./sqrt(4.*M_PI);
    r[PG]=1.;
    *xm = 0.5;
    break;
   case 18:
    /** case 18: Ryu and Jones test 2b.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=3./sqrt(4.*M_PI); l[BY]=6./sqrt(4.*M_PI); l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.1; r[VX]=0.; r[VY]=2.; r[VZ]=1.;
    r[BX]=3./sqrt(4.*M_PI); r[BY]=1./sqrt(4.*M_PI); r[BZ]=0.;
    r[PG]=10.;
    *xm = 0.5;
    break;
   case 19:
    /** case 19: Ryu and Jones test 3a.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=0.1; l[VX]=50.; l[VY]=l[VZ]=0.;
    l[BX]=0.; l[BY]=-1./sqrt(4.*M_PI); l[BZ]=-2./sqrt(4.*M_PI);
    l[PG]=0.4;
    r[RO]=0.1; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=0.; r[BY]=1./sqrt(4.*M_PI); r[BZ]=2./sqrt(4.*M_PI);
    r[PG]=0.2;
    *xm = 0.5;
    break;
   case 20:
    /** case 20: Ryu and Jones test 3b.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=-1.; l[VY]=l[VZ]=0.;
    l[BX]=0.; l[BY]=1.; l[BZ]=0.;
    l[PG]=1.;
    r[RO]=1.; r[VX]=1.; r[VY]=r[VZ]=0.;
    r[BX]=0.; r[BY]=1.; r[BZ]=0.;
    r[PG]=1.;
    *xm = 0.5;
    break;
   case 21:
    /** case 21: Ryu and Jones test 4a.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=1.; l[BY]=1.; l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.2; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=1.;r[BY]=0.;r[BZ]=0.;
    r[PG]=0.1;
    *xm = 0.5;
    break;
   case 22:
    /** case 22: Ryu and Jones test 4b.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=0.4; l[VX]=-0.66991; l[VY]=0.98263; l[VZ]=0.;
    l[BX]=1.3; l[BY]=0.0025293; l[BZ]=0.;
    l[PG]=0.52467;
    r[RO]=1.; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=1.3; r[BY]=1.; r[BZ]=0.;
    r[PG]=1.;
    *xm = 0.5;
    break;
   case 23:
    /** case 23: Ryu and Jones test 4c.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=0.65; l[VX]=0.667; l[VY]=-0.257; l[VZ]=0.;
    l[BX]=0.75; l[BY]=0.55; l[BZ]=0.;
    l[PG]=0.5;
    r[RO]=1.; r[VX]=0.4; r[VY]=-0.94; r[VZ]=0.;
    r[BX]=0.75; r[BY]=0.; r[BZ]=0.;
    r[PG]=0.75;
    *xm = 0.5;
    break;
   case 24:
    /** case 24: Ryu and Jones test 4d.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=0.7; l[BY]=0.;l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.3; r[VX]=r[VY]=0.; r[VZ]=1.;
    r[BX]=0.7; r[BY]=1.; r[BZ]=0.;
    r[PG]=0.2;
    *xm = 0.5;
    break;
   case 25:
    /** case 25: Ryu and Jones test 5a.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=0.75; l[BY]=1.; l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.125; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=0.75; r[BY]=-1.; r[BZ]=0.;
    r[PG]=0.1;
    *xm = 0.5;
    break;
   case 26:
    /** case 26: Ryu and Jones test 5b.\n*/
    if(SimPM.eqntype!=EQMHD && SimPM.eqntype!=EQGLM) 
    {cerr<<"(dataio_text::get_riemann_ics) Not using MHD but asking for MHD test problem. Exiting.\n"; return(1);}
    SimPM.gamma = 5./3.; cout <<"\t Forcing gamma=5/3 for test problem.\n";
    l[RO]=1.; l[VX]=l[VY]=l[VZ]=0.;
    l[BX]=1.3; l[BY]=1.; l[BZ]=0.;
    l[PG]=1.;
    r[RO]=0.4; r[VX]=r[VY]=r[VZ]=0.;
    r[BX]=1.3; r[BY]=-1.; r[BZ]=0.;
    r[PG]=0.4;
    *xm = 0.5;
    break;
   default:
    cout <<"Error: only know 26 tests, but sw!={1,..,26}"<<"\n";
    return(1);
  }
  cout <<"(dataio_text::get_riemann_ics) Got test number: "<<sw<<"\n";
  return(0);
}



// ##################################################################
// ##################################################################


int dataio_text::add_noise2data(
      class SimParams &SimPM,  ///< pointer to simulation parameters
      int n,
      double frac
      )
{
  srand(9768975);
  class cell *cpt;  double avg=0.; long int ct=0;
  switch (n) {
  case 1:
    cout <<"\tAdding random noise to pressure at fractional level of "<<frac<<"\n";
    cpt=gp->FirstPt();
    do {avg += cpt->P[PG]; ct++;} while ( (cpt=gp->NextPt(cpt)) !=0);
    cout <<"avg = "<<avg<< "\t ct= "<<ct<<"\n";
    avg /= static_cast<double>(ct);
    avg *= frac;  // avg is now a fraction 'frac' of the mean pressure on the grid.
    cout <<"avg = "<<avg<<"\n";
    cpt=gp->FirstPt(); do {
      if(cpt->isedge==0) {    // Don't want to alter any edge cells.
  //cout <<"PG before : "<<cpt->P[PG];
  cpt->P[PG] += avg*(static_cast<double>(rand())/RAND_MAX -0.5);
  //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
      }
    } while ( (cpt=gp->NextPt(cpt)) !=0);
    break;
    
  case 2:
    cout <<"Adding adiabatic random noise to cells at fractional level of "<<frac<<"\n";
    cpt=gp->FirstPt(); do {
      double temp;
      if(cpt->isedge==0) {    // Don't want to alter any edge cells.
  //cout <<"PG before : "<<cpt->P[PG];
  temp = 2.*frac*(static_cast<double>(rand())/RAND_MAX -0.5);
  cpt->P[PG] *= 1+temp;
  cpt->P[RO] *= exp(log(1+temp)/SimPM.gamma);
  //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
      }
    } while ( (cpt=gp->NextPt(cpt)) !=0);
    break;
    
  case 3:
    cout <<"Adding adiabatic wave to cells at fractional level of "<<frac<<"\n";
    // First get average pressure value.
    cpt=gp->FirstPt();
    do {avg += cpt->P[PG]; ct++;} while ( (cpt=gp->NextPt(cpt)) !=0);
    cout <<"avg = "<<avg<< "\t ct= "<<ct<<"\n";
    avg /= static_cast<double>(ct);
    // Now add wave to preshock state.
    cpt=gp->FirstPt(); do {
      double temp;
      if(cpt->isedge==0 && cpt->P[PG]<avg) {    // Don't want to alter any edge cells.
  //cout <<"PG before : "<<cpt->P[PG];
  temp = frac*sin(2.*M_PI*(CI.get_dpos(cpt,YY)/SimPM.Range[YY]) *(SimPM.NG[YY]/50));
  cpt->P[PG] *= 1+temp;
  cpt->P[RO] *= exp(log(1+temp)/SimPM.gamma);
  //cout <<"\tPG after : "<<cpt->P[PG]<<"\n";
      }
    } while ( (cpt=gp->NextPt(cpt)) !=0);
    
    break;
    
  default:
    cerr <<"\t Error, don't know what type of noise corresponds to "<<n<<"...\n";
    return(1);
  }
  return(0);
} // add_noise2data()



// ##################################################################
// ##################################################################



