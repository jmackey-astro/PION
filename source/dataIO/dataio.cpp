/// \file dataio.cc
/// \author Jonathan Mackey
/// Class definitions for ASCII Text data I/O.
/// 
/// modified:\n
///  - 2007-10-25 got the parallel fits-io class to compile and write data.
///  - 2007-10-26 parallel fits-io class reads data from single and multiple
///    files. Put in serial ifdefs, so that it should work for serial code too.
///  - 2008-09-19 Moved fits I/O to dataio_fits.cc, and moved text dataio from
///     gridMethods.cc to here in a new class dataio_text.
///
///  - 2010-02-03 JM: removed unused variables from functions.
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2010-07-20/22 JM: Work on new dataio structure with a list of parameters
///    to read and write.  read/write_simulation_parameters() functions work now.
/// - 2010-09-03 JM: analysis software tried to add new RSP source
///    every read, so I changed the function to only add sources if
///    the source list is empty.
/// - 2010.10.01 JM: Spherical coordinates added.
///    Got rid of testing myalloc/myfree commands.
/// - 2010.10.05 JM: Added an extra parameter "Rstar" for stellar winds.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2010.11.03 JM: New Ndigits variable added.  Removed 'endl' statements.
/// - 2010.11.21 JM: Added more viscosity flags for dataio_text
/// - 2011.01.06 JM: New stellar wind interface.
/// - 2011.02.15 JM: Added new stellar wind parameters for evolving wind sources.
/// - 2011.02.24 JM: Added read/write for multiple radiation sources, with
///    additional parameters. Much simplified RT I/O by using new struct SimPM.RS.
/// - 2011.02.25 JM: got rid of HCORR ifdef wrapper.
/// - 2011.02.28 JM: got rid of RSP references. 03.01 few bug fixes for new RT.
/// - 2011.03.02 JM: Better support for tracer variables (with or without names).
///        Can handle an arbitrary number of tracers now.
/// - 2011.03.21 JM: Added extra radiation source parameters. (and 22.03)
/// - 2011.05.02 JM: New RT params.
/// - 2011.06.02 JM: WriteHeader() added; some minor mods to dataio_text.
/// - 2011.10.05 JM: Fixed problems in the wind-source reading (for analysis code).
/// - 2011.11.22 JM: Added Stellar wind t_scalefactor parameter.
/// - 2011.12.14 JM: Added radiation source EvoFile parameter (time-varying sources).
/// - 2012.03.23 JM: Removed warning about not finding Nbc as parameter (it is
///    obselete now).
/// - 2013.02.14 JM: Added He/Metal mass fractions as EP parameters,
///    to make metallicity and mu into parameterfile settings.
/// - 2013.02.19 JM: Moved file_status class definitions to new file.
/// - 2013.08.19 JM: Added Hydrogen MassFrac to EP parameter list.
/// - 2013.08.20 JM: Modified cell_interface for optical depth vars.
/// - 2013.09.05 JM: changed logic of writing T/Eint in ascii data.
/// - 2013.09.16 JM: Increased precision of ascii data to 14 digits.
/// - 2015.01.15 JM: Added new include statements for new PION version.
/// - 2015.07.0[6-8] JM: Started to change tracer setup in files.
/// - 2015.08.05 JM: tidied up code; added pion_flt datatype.
/// - 2015.10.19 JM: Fixed dvararr to always use pion_flt correctly.
/// - 2017.11.07-22 JM: updating boundary setup.
/// - 2018.01.24 JM: worked on making SimPM non-global

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#include "constants.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "dataIO/dataio.h"
#include "dataIO/readparams.h"
#include "grid/stellar_wind_BC.h"
#include "raytracing/raytracer_base.h"

#include <sstream>
using namespace std;


// ##################################################################
// ##################################################################




// -----------------------------------------------------
// ---------- PARAMETER CLASSES ------------------------

//
// default constructor
//
pm_int::pm_int()
 {type=MY_INT; defval=-1;}
pm_double::pm_double()
 {type=MY_DOUBLE; defval=-VERY_LARGE_VALUE;}
pm_float::pm_float()
 {type=MY_FLOAT; defval=-SMALLVALUE;}
pm_long::pm_long()
 {type=MY_LONG; defval=-1;}
pm_string::pm_string()
{type=MY_STRING; defval.assign("hello");}
pm_ddimarr::pm_ddimarr()
 {type=MY_DDIMARR; len=MAX_DIM; defval=0;}
pm_idimarr::pm_idimarr()
 {type=MY_IDIMARR; len=MAX_DIM; defval=0;}
pm_dvararr::pm_dvararr()
 {type=MY_DVARARR; len=MAX_NVAR; defval=0;}


// ##################################################################
// ##################################################################



//
// constructor with name only
//
pm_int::pm_int(const string s)
{
  type=MY_INT;    name.assign(s); defval=-1;
}
pm_double::pm_double(const string s)
{
  type=MY_DOUBLE; name.assign(s); defval=-VERY_LARGE_VALUE;
}
pm_float::pm_float(const string s)
{
  type=MY_FLOAT;  name.assign(s); defval=-SMALLVALUE;
}
pm_long::pm_long(const string s)
{
  type=MY_LONG;   name.assign(s); defval=-1;
}
pm_string::pm_string(const string s)
{
  type=MY_STRING; name.assign(s); defval.assign("HELLO");
}
pm_ddimarr::pm_ddimarr(const string s)
{
  type=MY_DDIMARR; len=MAX_DIM; name.assign(s); defval=0;
}
pm_idimarr::pm_idimarr(const string s)
{
  type=MY_IDIMARR; len=MAX_DIM; name.assign(s); defval=0;
}
pm_dvararr::pm_dvararr(const string s)
{
  type=MY_DVARARR; len=MAX_NVAR; name.assign(s); defval=0;
}


// ##################################################################
// ##################################################################



//
// constructor with name and pointer to data
//
pm_int::pm_int(const string s, int *p)
{
  type=MY_INT; name.assign(s); ptr=p; defval=-1;
}
pm_double::pm_double(const string s, double *p)
{
  type=MY_DOUBLE; name.assign(s); ptr=p; defval=-VERY_LARGE_VALUE;
}
pm_float::pm_float(const string s, float *p)
{
  type=MY_FLOAT; name.assign(s); ptr=p; defval=-SMALLVALUE;
}
pm_long::pm_long(const string s, long int *p)
{
  type=MY_LONG; name.assign(s); ptr=p; defval=-1;
}
pm_string::pm_string(const string s, string *p)
{
  type=MY_STRING; name.assign(s); ptr=p; defval.assign("HELLO");
}
pm_ddimarr::pm_ddimarr(const string s, double *p)
{
  type=MY_DDIMARR; len=MAX_DIM; name.assign(s);
  ptr=p; defval=0;
}
pm_idimarr::pm_idimarr(const string s, int *p)
{
  type=MY_IDIMARR; len=MAX_DIM; name.assign(s);
  ptr=p; defval=0;
}
pm_dvararr::pm_dvararr(const string s, pion_flt *p)
{
  type=MY_DVARARR; len=MAX_NVAR; name.assign(s);
  ptr=p; defval=0;
}


// ##################################################################
// ##################################################################



//
// constructor with name, pointer to data, default value.
//
pm_int::pm_int(const string s, int *p, const int def)
{
  type=MY_INT; name.assign(s); ptr=p; defval = def;
}
pm_double::pm_double(const string s, double *p, const double def)
{
  type=MY_DOUBLE; name.assign(s); ptr=p; defval=def;
}
pm_float::pm_float(const string s, float *p, const float def)
{
  type=MY_FLOAT; name.assign(s); ptr=p; defval=def;
}
pm_long::pm_long(const string s, long int *p, const long int def)
{
  type=MY_LONG; name.assign(s); ptr=p; defval=def;
}
pm_string::pm_string(const string s, string *p, const string def)
{
  type=MY_STRING; name.assign(s); ptr=p; defval.assign(def);
}
pm_ddimarr::pm_ddimarr(const string s, double *p, const double *def)
{
  type=MY_DDIMARR; len=MAX_DIM; name.assign(s); ptr=p;
  defval = mem.myalloc(defval,len);
  for (int v=0;v<len;v++) defval[v] = def[v];
}
pm_idimarr::pm_idimarr(const string s, int *p, const int *def)
{
  type=MY_IDIMARR; len=MAX_DIM; name.assign(s); ptr=p;
  defval = mem.myalloc(defval,len);
  for (int v=0;v<len;v++) defval[v] = def[v];
}
pm_dvararr::pm_dvararr(const string s, pion_flt *p, const pion_flt *def)
{
  type=MY_DVARARR; len=MAX_NVAR; name.assign(s); ptr=p;
  defval = mem.myalloc(defval,len);
  for (int v=0;v<len;v++) defval[v] = def[v];
}


// ##################################################################
// ##################################################################



//
// Some have destructors:
//
pm_ddimarr::~pm_ddimarr() {
  if (defval) defval = mem.myfree(defval);
  ptr=0;
}
pm_idimarr::~pm_idimarr() {
  if (defval) defval = mem.myfree(defval);
  ptr=0;
}
pm_dvararr::~pm_dvararr() {
  if (defval) defval = mem.myfree(defval);
  ptr=0;
}


// ##################################################################
// ##################################################################




void pm_int::assign_val(void *val)
 {*ptr=*(static_cast<int *>(val));}
void pm_double::assign_val(void *val)
 {*ptr=*(static_cast<double *>(val));}
void pm_float::assign_val(void *val)
 {*ptr=*(static_cast<float *>(val));}
void pm_long::assign_val(void *val)
 {*ptr = *(static_cast<long int *>(val));}
void pm_string::assign_val(void *val)
 {(*ptr).assign(*(static_cast<string *>(val)));}
void pm_ddimarr::assign_val(void *val)
{
  for (int i=0;i<len;i++) ptr[i]= (static_cast<double *>(val))[i];
}
void pm_idimarr::assign_val(void *val)
{
  for (int i=0;i<len;i++) ptr[i]= (static_cast<int *>(val))[i];
}
void pm_dvararr::assign_val(void *val)
{
  for (int i=0;i<len;i++) ptr[i]= (static_cast<pion_flt *>(val))[i];
}


// ##################################################################
// ##################################################################



void pm_int::set_ptr(void *p) {ptr=static_cast<int *>(p);}
void pm_double::set_ptr(void *p) {ptr=static_cast<double *>(p);}
void pm_float::set_ptr(void *p) {ptr=static_cast<float *>(p);}
void pm_long::set_ptr(void *p) {ptr=static_cast<long int *>(p);}
void pm_string::set_ptr(void *p) {ptr=static_cast<string *>(p);}
void pm_ddimarr::set_ptr(void *p) {ptr=static_cast<double *>(p);}
void pm_idimarr::set_ptr(void *p) {ptr=static_cast<int *>(p);}
void pm_dvararr::set_ptr(void *p) {ptr=static_cast<pion_flt *>(p);}


// ##################################################################
// ##################################################################




void pm_int::show_val() {cout<<*ptr;}
void pm_double::show_val() {cout<<*ptr;}
void pm_float::show_val() {cout<<*ptr;}
void pm_long::show_val() {cout<<*ptr;}
void pm_string::show_val() {cout<<*ptr<<", name="<<name;}
void pm_ddimarr::show_val()
{
  cout <<"[";
  for (int i=0;i<len-1;i++) cout<<ptr[i]<<", ";
  cout <<ptr[len-1]<<"]";
}
void pm_idimarr::show_val()
{
  cout <<"[";
  for (int i=0;i<len-1;i++) cout<<ptr[i]<<", ";
  cout <<ptr[len-1]<<"]";
}
void pm_dvararr::show_val()
{
  cout <<"[";
  for (int i=0;i<len-1;i++) cout<<ptr[i]<<", ";
  cout <<ptr[len-1]<<"]";
}



// ##################################################################
// ##################################################################



void pm_int::set_to_default() {*ptr=defval;}
void pm_double::set_to_default() {*ptr=defval;}
void pm_float::set_to_default() {*ptr=defval;}
void pm_long::set_to_default() {*ptr=defval;}
void pm_string::set_to_default() {(*ptr).assign(defval);}
void pm_ddimarr::set_to_default() {
  if (!defval) rep.error("No default value!","SET ME!");
  else for (int i=0;i<len;i++) ptr[i] = defval[i];
}
void pm_idimarr::set_to_default() {
  if (!defval) rep.error("No default value!","SET ME!");
  else for (int i=0;i<len;i++) ptr[i] = defval[i];
}

void pm_dvararr::set_to_default() {
  if (!defval) rep.error("No default value!","SET ME!");
  else for (int i=0;i<len;i++) ptr[i] = defval[i];
}


// ##################################################################
// ##################################################################



void * pm_int::get_ptr() {return static_cast<void *>(ptr);}
void * pm_double::get_ptr() {return static_cast<void *>(ptr);}
void * pm_float::get_ptr() {return static_cast<void *>(ptr);}
void * pm_long::get_ptr() {return static_cast<void *>(ptr);}
void * pm_string::get_ptr() {return static_cast<void *>(ptr);}
void * pm_ddimarr::get_ptr() {return static_cast<void *>(ptr);}
void * pm_idimarr::get_ptr() {return static_cast<void *>(ptr);}
void * pm_dvararr::get_ptr() {return static_cast<void *>(ptr);}



// ##################################################################
// ##################################################################




void     pm_int::set_default_val(void *v) {defval= *(static_cast<int *>(v));}
void  pm_double::set_default_val(void *v) {defval= *(static_cast<double *>(v));}
void   pm_float::set_default_val(void *v) {defval= *(static_cast<float *>(v));}
void    pm_long::set_default_val(void *v) {defval= *(static_cast<long *>(v));}
void  pm_string::set_default_val(void *v) {defval= *(static_cast<string *>(v));}
void pm_ddimarr::set_default_val(void *v) {
  if (!defval) defval=mem.myalloc(defval,len);
  else for (int i=0;i<len;i++) defval[i]=(static_cast<double *>(v))[i];
}
void pm_idimarr::set_default_val(void *v) {
  if (!defval) defval=mem.myalloc(defval,len);
  else for (int i=0;i<len;i++) defval[i]=(static_cast<double *>(v))[i];
}
void pm_dvararr::set_default_val(void *v) {
  if (!defval) defval=mem.myalloc(defval,len);
  else for (int i=0;i<len;i++) defval[i]=(static_cast<pion_flt *>(v))[i];
}



// void pm_int::
// void pm_double::
// void pm_float::
// void pm_long::
// void pm_string::
// void pm_ddimarr::
// void pm_idimarr::
// void pm_dvararr::
// ---------- PARAMETER CLASSES ------------------------
// -----------------------------------------------------

// -----------------------------------------------------
// --------- BASE DATA I/O CLASS DEFINITIONS -----------
// ------------------ NEW IN SVN-R202 ------------------
// -----------------------------------------------------


// ##################################################################
// ##################################################################


DataIOBase::DataIOBase(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
  : Ndigits(8)
{
  params.clear();
  jet_pm.clear();
  rt_src.clear();
  windsrc.clear();
  set_params(SimPM);
  return;
}



// ##################################################################
// ##################################################################


DataIOBase::~DataIOBase()
{
  clear_param_list(params);
  clear_param_list(jet_pm);
  clear_param_list(rt_src);
  clear_param_list(windsrc);

      
  return;
}



// ##################################################################
// ##################################################################


void DataIOBase::set_params(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // make sure list is empty:
  //
  if (!params.empty()) {
    cerr <<"WARNING! params list not empty. CLEARING IT, BUT THIS SHOULDN'T HAPPEN!\n";
    params.clear();
  }

  //
  // now create parameter structs and add them to the list.
  // CRITICAL parameters don't need a default value since the code will bug 
  // out if they are missing.  Other params should get one.
  //
  // START WITH GRID PROPERTIES
  //
  pm_base *p;

  pm_int     *p001 = new pm_int     
    ("gridtype",    &SimPM.gridType);
  p = p001; p->critical=true;  
  params.push_back(p);
  pm_int     *p002 = new pm_int     
    ("gridndim",    &SimPM.ndim);
  p = p002; p->critical=true;  
  params.push_back(p);
  pm_idimarr *p003 = new pm_idimarr 
    ("NGrid",       SimPM.NG);
  p = p003; p->critical=true;  
  params.push_back(p);
  pm_long    *p004 = new pm_long    
    ("Ncell",       &SimPM.Ncell);
  p = p004; p->critical=true;  
  params.push_back(p);
  pm_ddimarr *p005 = new pm_ddimarr 
    ("Xmin",        SimPM.Xmin);
  p = p005; p->critical=true;  
  params.push_back(p);
  pm_ddimarr *p006 = new pm_ddimarr 
    ("Xmax",        SimPM.Xmax);
  p = p006; p->critical=true;  
  params.push_back(p);

  //
  // Nested grid parameters
  //
  pm_int     *n001 = new pm_int     
    ("grid_nlevels",      &SimPM.grid_nlevels);
  p = n001; p->critical=false;  
  params.push_back(p);

  pm_idimarr *n002 = new pm_idimarr     
    ("grid_aspect_ratio", SimPM.grid_aspect_ratio);
  p = n002; p->critical=false;  
  params.push_back(p);

  pm_ddimarr *n003 = new pm_ddimarr 
    ("grid_nest_centre",  SimPM.grid_nest_centre);
  p = n003; p->critical=false;  
  params.push_back(p);



  //
  // Boundary conditions
  //
  // Number of internal boundaries
  pm_int     *p124 = new pm_int     
    ("BC_Ninternal",   &SimPM.BC_Nint);
  p = p124; p->critical=false;  
  params.push_back(p);
  have_setup_bc_pm=false; // so we know to populate list later.

  // LEGACY: look for a string parameter called BC, in case the
  // data file was written with a pre-2018 version of PION.
  pm_string     *p125 = new pm_string     
    ("typeofbc_str",   &SimPM.BC_STRING);
  p = p125; p->critical=false;  
  params.push_back(p);


  //
  // EQUATIONS
  //
  pm_int     *p008 = new pm_int     
    ("eqn_type",   &SimPM.eqntype);
  p = p008; p->critical=true;  
  params.push_back(p);
  pm_int     *p009 = new pm_int     
    ("eqn_ndim",   &SimPM.eqnNDim, 3);
  p = p009; p->critical=false;  
  params.push_back(p);
  pm_int     *p010 = new pm_int     
    ("eqn_nvar",   &SimPM.nvar);
  p = p010; p->critical=true;  
  params.push_back(p);

  //
  // Tracers
  //
  pm_int     *p011 = new pm_int     
    ("num_tracer", &SimPM.ntracer);
  p = p011; p->critical=true;  
  params.push_back(p);
  pm_string  *p012a = new pm_string  
    ("chem_code", &SimPM.chem_code, "CHEM_CODE");
  p = p012a; p->critical=false;  
  params.push_back(p);
  have_setup_tracers=false; // so we know to populate list later.
  // LEGACY: look for a string parameter called trtype, in case the
  // data file was written with a pre-2018 version of PION.
  pm_string     *p312 = new pm_string     
    ("tracer_str",   &SimPM.TRTYPE);
  p = p312; p->critical=false;  
  params.push_back(p);

  //
  // hydro solver parameters.
  //
  pm_int     *p013 = new pm_int     
    ("solver",     &SimPM.solverType);
  p = p013; p->critical=true;  
  params.push_back(p);
  pm_int     *p014 = new pm_int     
    ("coord_sys",  &SimPM.coord_sys);
  p = p014; p->critical=true;  
  params.push_back(p);
  pm_int     *p015 = new pm_int     
    ("Space_OOA", &SimPM.spOOA, 2);
  p = p015; p->critical=false;  
  params.push_back(p);
  pm_int     *p016 = new pm_int     
    ("Time_OOA", &SimPM.tmOOA, 2);
  p = p016; p->critical=false;  
  params.push_back(p);
  pm_double  *p017 = new pm_double  
    ("Gamma",      &SimPM.gamma);
  p = p017; p->critical=true;  
  params.push_back(p);
  pm_double  *p018 = new pm_double  
    ("CFL",        &SimPM.CFL);
  p = p018; p->critical=true;  
  params.push_back(p);
  pm_int     *p019 = new pm_int     
    ("art_visc",   &SimPM.artviscosity);
  p = p019; p->critical=true;  
  params.push_back(p);
  pm_double  *p020 = new pm_double  
    ("eta_visc",   &SimPM.etav);
  p = p020; p->critical=true;  
  params.push_back(p);

  pm_dvararr *p120 = new pm_dvararr 
    ("Ref_Vector",  SimPM.RefVec);
  p = p120; p->critical=true;  
  params.push_back(p);

  //
  // PHYSICS FLAGS
  //
  pm_int     *p021 = new pm_int     
    ("EP_dynamics",          &SimPM.EP.dynamics);
  p = p021; p->critical=true;  
  params.push_back(p);
  pm_int     *p022 = new pm_int     
    ("EP_raytracing",        &SimPM.EP.raytracing);
  p = p022; p->critical=true;  
  params.push_back(p);
  pm_int     *p023 = new pm_int     
    ("EP_cooling",           &SimPM.EP.cooling);
  p = p023; p->critical=true;  
  params.push_back(p);
  pm_int     *p024 = new pm_int     
    ("EP_chemistry",         &SimPM.EP.chemistry);
  p = p024; p->critical=true;  
  params.push_back(p);
  pm_int     *p025 = new pm_int     
    ("EP_coll_ionisation",   &SimPM.EP.coll_ionisation);
  p = p025; p->critical=true;  
  params.push_back(p);
  pm_int     *p026 = new pm_int     
    ("EP_phot_ionisation",   &SimPM.EP.phot_ionisation);
  p = p026; p->critical=true;  
  params.push_back(p);
  pm_int     *p027 = new pm_int     
    ("EP_rad_recombination", &SimPM.EP.rad_recombination);
  p = p027; p->critical=true;  
  params.push_back(p);
  pm_int     *p028 = new pm_int     
    ("EP_update_erg",        &SimPM.EP.update_erg);
  p = p028; p->critical=true;  
  params.push_back(p);
  pm_int     *p029 = new pm_int     
    ("EP_MP_timestep_limit", &SimPM.EP.MP_timestep_limit, 0);
  p = p029; p->critical=false;
  params.push_back(p);
//#ifdef SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE
  pm_double  *pMNT = new pm_double
    ("EP_Min_Temperature", &SimPM.EP.MinTemperature, 0.0);
  p = pMNT; p->critical=false;
  params.push_back(p);
  pm_double  *pMXT = new pm_double
    ("EP_Max_Temperature", &SimPM.EP.MaxTemperature, 1.0e100);
  p = pMXT; p->critical=false;
  params.push_back(p);
//#endif // SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE

  //
  // Hydrogen abundance (by mass) X.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double  *pXXX = new pm_double
    ("EP_Hydrogen_MassFrac", &SimPM.EP.H_MassFrac, 0.7154);
  p = pXXX; p->critical=false;
  params.push_back(p);

  //
  // Helium abundance (by mass) Y.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double  *pYYY = new pm_double
    ("EP_Helium_MassFrac", &SimPM.EP.Helium_MassFrac, 0.2703);
  p = pYYY; p->critical=false;
  params.push_back(p);
  //
  // Metal abundance (by mass) Z.
  // Default value is from Asplund et al. (2009,ARA&A,47,481)
  //
  pm_double  *pZZZ = new pm_double
    ("EP_Metal_MassFrac", &SimPM.EP.Metal_MassFrac, 0.0142);
  p = pZZZ; p->critical=false;
  params.push_back(p);

  //
  // TIMESTEPS
  //
  pm_double  *p030 = new pm_double  
    ("t_start",      &SimPM.starttime);
  p = p030; p->critical=true;  
  params.push_back(p);
  pm_double  *p031 = new pm_double  
    ("t_finish",     &SimPM.finishtime);
  p = p031; p->critical=true;  
  params.push_back(p);
  pm_int     *p032 = new pm_int     
    ("t_step",       &SimPM.timestep, 0);
  p = p032; p->critical=false;  
  params.push_back(p);
  pm_double  *p033 = new pm_double  
    ("t_sim",        &SimPM.simtime, 0.0);
  p = p033; p->critical=false;  
  params.push_back(p);
  pm_double  *p034 = new pm_double  
    ("min_timestep", &SimPM.min_timestep, TINYVALUE);
  p = p034; p->critical=false;  
  params.push_back(p);

  //
  // OUTPUT PARAMETERS
  //
  pm_int     *p040 = new pm_int     
    ("typeofop",   &SimPM.typeofop);
  p = p040; p->critical=true;  
  params.push_back(p);
  pm_string  *p041 = new pm_string  
    ("outfile",    &SimPM.outFileBase);
  p = p041; p->critical=true;  
  params.push_back(p);
  pm_int     *p042 = new pm_int     
    ("op_freq",    &SimPM.opfreq, 0);
  p = p042; p->critical=false;  
  params.push_back(p);
  pm_double  *p043 = new pm_double  
    ("opfreq_time",  &SimPM.opfreq_time, -1.0);
  p = p043; p->critical=false;  
  params.push_back(p);
  pm_int     *p044 = new pm_int     
    ("op_criterion", &SimPM.op_criterion, 0);
  p = p044; p->critical=false;  
  params.push_back(p);

  //
  // JET SIMULATION PARAMS (GET WIDTH/STAT LATER)
  //
  pm_int     *p050 = new pm_int     
    ("JetSim",     &JP.jetic, 0);
  p = p050; p->critical=false;  
  params.push_back(p);
  have_setup_jet_pm=false; // so we know to populate list later.

  //
  // RAY-TRACING PARAMS (GET STRENGTH/POSN LATER)
  //
  pm_int     *p060 = new pm_int     
    ("RT_Nsources",   &(SimPM.RS.Nsources), 0);
  p = p060; p->critical=false;  
  params.push_back(p);
  have_setup_rt_src=false; // so we know to populate list later.

  //
  // STELLAR WIND PARAMS (GET STRENGTH/POSN LATER)
  //
  pm_int     *w200 = new pm_int     
    ("WIND_Nsources",   &SWP.Nsources, 0);
  p = w200; p->critical=false;  
  params.push_back(p);
  have_setup_windsrc=false; // so we know to populate list later.

  //
  // UNITS (NOTE THESE ARE JUST FOR INFORMATION PURPOSES NOW!)
  //
  pm_string  *u001 = new pm_string  
    ("unitsys", &uc.unitsys, "CGS");
  p = u001; p->critical=false;  
  params.push_back(p);
  pm_string  *u002 = new pm_string  
    ("unitdens", &uc.density, "g.cm-3");
  p = u002; p->critical=false;  
  params.push_back(p);
  pm_string  *u003 = new pm_string  
    ("unitlen", &uc.length, "cm");
  p = u003; p->critical=false;  
  params.push_back(p);
  pm_string  *u004 = new pm_string  
    ("unitvel", &uc.velocity, "cm.s-1");
  p = u004; p->critical=false;  
  params.push_back(p);
  pm_string  *u005 = new pm_string  
    ("unitmagf", &uc.bfield, "Gauss/sqrt(4pi)");
  p = u005; p->critical=false;  
  params.push_back(p);
  //
  // For the values, I think the idea is that one code unit equals
  // this many units in the current unit system.  But I can't remember.
  //
  pm_double  *u006 = new pm_double  
    ("rhoval", &uc.rhoVal, 1.0);
  p = u006; p->critical=false;  
  params.push_back(p);
  pm_double  *u007 = new pm_double  
    ("lenval", &uc.lenVal, 1.0);
  p = u007; p->critical=false;  
  params.push_back(p);
  pm_double  *u008 = new pm_double  
    ("velval", &uc.velVal, 1.0);
  p = u008; p->critical=false;  
  params.push_back(p);
  pm_double  *u009 = new pm_double  
    ("magval", &uc.magVal, 1.0);
  p = u009; p->critical=false;  
  params.push_back(p);
  

  return;
}



// ##################################################################
// ##################################################################


void DataIOBase::clear_param_list(std::list<class pm_base *> &listptr)
{
  if (listptr.empty()) 
    return;
  else {
    do {
      list<pm_base *>::iterator i=listptr.begin();
      pm_base *p=(*i);
      p->set_ptr(0);
      delete p;
      listptr.erase(i);
    } while (!listptr.empty());
  }
  return;
}



// ##################################################################
// ##################################################################


int DataIOBase::read_simulation_parameters(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (params.empty())
    rep.error("Parameter list is empty -- make sure it populates itself!!",0);

  //
  // loop over all parameters.
  //
  pm_base *p=0;
  int err=0;
  for (list<pm_base *>::iterator iter=params.begin(); iter!=params.end(); ++iter) {
    p = (*iter);
    err = read_header_param(p);
    if (err) {
      if (p->critical) {
        rep.error("Error reading parameter",p->name);
      }
      else {
        //cout <<"parameter "<<p->name<<" not found. setting to default val.\n";
        p->set_to_default();
        err =0;
      }
    }
  }

  //
  // Read boundary conditions for each edge and internal BCs
  // Two options: lecacy code uses BC_STRING, new code uses individual
  // parameters.  Test for legacy code first.
  //
  if (SimPM.BC_STRING != "") {
    int i=0;
    string::size_type pos;
    string temp;

    // First get the external boundaries
    string d[6] = {"XN","XP","YN","YP","ZN","ZP"};
    string *par[6];
    par[0] = &(SimPM.BC_XN);
    par[1] = &(SimPM.BC_XP);
    par[2] = &(SimPM.BC_YN);
    par[3] = &(SimPM.BC_YP);
    par[4] = &(SimPM.BC_ZN);
    par[5] = &(SimPM.BC_ZP);
    //rep.printVec("par",par,6);
    for (i=0; i<2*SimPM.ndim; i++) {
      if ( (pos=SimPM.BC_STRING.find(d[i])) == string::npos)
        rep.error("Couldn't find boundary condition for ",d[i]);
      else {
        temp = SimPM.BC_STRING.substr(pos+2,3);
        if      (temp=="per") *par[i] = "periodic";
        else if (temp=="out") *par[i] = "outflow";
        else if (temp=="owo") *par[i] = "one-way-outflow";
        else if (temp=="inf") *par[i] = "inflow";
        else if (temp=="ref") *par[i] = "reflecting";
        else if (temp=="jrf") *par[i] = "equator-reflect";
        else if (temp=="fix") *par[i] = "fixed";
        else if (temp=="dmr") *par[i] = "DMR";
        else if (temp=="sb1") *par[i] = "SB1";
        else rep.error("Unrecognised BC type",SimPM.BC_STRING);
      }
    } // loop over external boundaries
    //cout <<SimPM.BC_XN <<"  "<<SimPM.BC_XP <<"  ";
    //cout <<SimPM.BC_YN <<"  "<<SimPM.BC_YP <<"  ";
    //cout <<SimPM.BC_ZN <<"  "<<SimPM.BC_ZP <<"\n";

    // Now look for any internal boundaries
    int len = SimPM.BC_STRING.length(); len = (len+5)/6;
    if (len > 2*SimPM.ndim) {
      SimPM.BC_Nint = len - 2*SimPM.ndim;
      if (!SimPM.BC_INT)
        SimPM.BC_INT = mem.myalloc(SimPM.BC_INT,SimPM.BC_Nint);
      //cout <<" got "<<SimPM.BC_Nint<<" internal boundaries\n";
    }
    for (i=2*SimPM.ndim; i<len; i++) {
      //cout <<"i="<<i<<", len="<<len<<"\n";
      if ( (pos=SimPM.BC_STRING.find("IN",i*6)) ==string::npos) {
        rep.error("internal boundary condition not found",SimPM.BC_STRING);
      }
      else {
        temp = SimPM.BC_STRING.substr(pos+2,3);
        if      (temp=="jet") SimPM.BC_INT[i-2*SimPM.ndim] = "jet";
        else if (temp=="dm2") SimPM.BC_INT[i-2*SimPM.ndim] = "DMR2";
        else if (temp=="rsh") SimPM.BC_INT[i-2*SimPM.ndim] = "RadShock";
        else if (temp=="rs2") SimPM.BC_INT[i-2*SimPM.ndim] = "RadShock2";
        else if (temp=="wnd") SimPM.BC_INT[i-2*SimPM.ndim] = "stellar-wind";
        else rep.error("Unrecognised INT BC type",SimPM.BC_STRING);
      }
    } // loop over internal boundaries

  } // if OLD LEGACY BC STRING

  else {
    // New boundary conditions format
    if (!have_setup_bc_pm) set_bc_pm_params(SimPM);
    if (bc_pm.empty()) rep.error("Boundary parameter list is empty!!",0);
    //
    // now read them:
    //
    int ct=0;
    for (list<pm_base *>::iterator iter=bc_pm.begin(); iter!=bc_pm.end(); ++iter) {
      p = (*iter);
      err = read_header_param(p);
      //cout <<"boundary list "<<ct<<", parameter "<<p->name<<"\n";
      if (err) rep.error("Error reading parameter",p->name);
      ct++;
    }
  }


  //
  // We now use num_tracer to set the position of the first tracer:
  //
  if (SimPM.ntracer>0) SimPM.ftr = SimPM.nvar-SimPM.ntracer;
  else                 SimPM.ftr = SimPM.nvar;
  //
  // Set up tracer parameters, based on ntracer and read them in
  //
  if (!have_setup_tracers) set_tracer_params(SimPM);
  
  //
  // Two options: first legacy code, which uses a string
  // called trtype to set all tracers, or new code, which
  // uses chem_code and then a list called TracerIJK.
  //
  if (SimPM.TRTYPE != "") {
    SimPM.chem_code = SimPM.TRTYPE.substr(0,6);
    int len = (SimPM.TRTYPE.length() +5)/6 -1;
    if (len!=SimPM.ntracer)
      rep.error("bad tracer string (LEGACY)",SimPM.TRTYPE);
    for (int i=0;i<len;i++) {
      SimPM.tracers[i] = SimPM.TRTYPE.substr(6*(i+1),6);
      //cout <<"tracer["<<i<<"] = "<<SimPM.tracers[i] <<"\n";
    }
      
  } // if LEGACY tracer variables used.

  else {
    for (list<pm_base *>::iterator iter=tr_pm.begin(); iter!=tr_pm.end(); ++iter) {
      p = (*iter);
      err = read_header_param(p);
      if (err) rep.error("Error reading parameter",p->name);
    }
  }


  //
  // Set SimPM.Range[] based on Xmin,Xmax
  //
  for (int i=0; i<SimPM.ndim; i++)
    SimPM.Range[i] = SimPM.Xmax[i]-SimPM.Xmin[i];

  //
  //
  // Read Jet parameters if doing a JET SIM
  //
  if (JP.jetic) {
    if (!have_setup_jet_pm) set_jet_pm_params();
    if (jet_pm.empty())
      rep.error("Jet parameter list is empty -- make sure it populates itself!!",0);
    //
    // now read them:
    //
    for (list<pm_base *>::iterator iter=jet_pm.begin(); iter!=jet_pm.end(); ++iter) {
      p = (*iter);
      err = read_header_param(p);
      if (err) rep.error("Error reading parameter",p->name);
    }
  }

  //
  // Now the radiation sources -- only read in the sources if we
  // haven't already done so!  This can happen when analysing multiple
  // timesteps.  First check if there are radiation sources -- SimPM.RS.Nsources
  // was already assigned by reading the header data from the input file.
  // THIS CODE ASSUMES THE NUMBER OF SOURCES DOES NOT CHANGE OVER TIME!!
  //
  //cout <<"RT-Nsrcs = "<<SimPM.RS.Nsources<<"\n";
  if (SimPM.RS.Nsources > 0) {
    //
    // First set up the list of parameters we need to read:
    //
    if (!have_setup_rt_src || SimPM.RS.sources.empty()) {
      have_setup_rt_src=false;
      clear_param_list(rt_src);
      set_rt_src_params(SimPM);
    }
    if (rt_src.empty()) {
      rep.error("RT-src parameter list is empty -- make sure it populates itself!!",0);
    }
    
    //
    // Now read properties for each source and populate SimPM.RS.sources vector.
    // Note the rt_src list has all parameters for Nsources radiation sources,
    // so Nsources*Nparams elements.
    //
    for (list<pm_base *>::iterator iter=rt_src.begin(); iter!=rt_src.end(); ++iter) {
      p = (*iter);
      err = read_header_param(p);
      if (err && p->critical) {
        rep.error("Error reading parameter",p->name);
      }
    }

    //
    // Now we have N sources in SimPM.RS.sources.  First assign an id
    // to each of them, and make sure their data got assigned.
    //
    for (int i=0; i<SimPM.RS.Nsources; i++) {
      SimPM.RS.sources.at(i).id = i;
      //
      // Set NTau based on source.effect.
      //
      if (SimPM.RS.sources.at(i).effect==RT_EFFECT_HHE_MFQ)
        SimPM.RS.sources.at(i).NTau = 4;
      else
        SimPM.RS.sources.at(i).NTau = 1;
      //
      // Check for sensible values:
      //
      if (SimPM.RS.sources.at(i).type <0) {
        rep.error("Failed to get source type for source id",i);
      }
      if (SimPM.RS.sources.at(i).at_infinity <0) {
        rep.error("Failed to get source_at_infty for source id",i);
      }
      if (SimPM.RS.sources.at(i).pos[XX] <-9.0e200) {
        rep.error("Failed to get source position for source id",i);
      }
      if (SimPM.RS.sources.at(i).EvoFile=="NOFILE") {
        //cout <<"\t|+|+|+|+|+| Non-evolving radiation source "<<i<<" found.\n";
      }
      else {
        //cout <<"\t|+|+|+|+|+| Evolving radiation source "<<i<<" detected = ";
        cout <<SimPM.RS.sources.at(i).EvoFile<<"\n";
      }
      if (SimPM.RS.sources.at(i).opacity_var+SimPM.ftr >=SimPM.nvar) {
        rep.error("Opacity var for source is off end of array (ftr offset!)",i);
      }
    }

  } // read RT data.

  //
  // Now read in wind parameters, if doing a wind simulation
  // THIS CODE ASSUMES THE NUMBER OF SOURCES DOES NOT CHANGE OVER TIME!!
  //
  if (SWP.Nsources > 0) {
    if (!have_setup_windsrc) set_windsrc_params();
    if (windsrc.empty())
      rep.error("wind-src parameter list is empty -- make sure it populates itself!!",0);
    //
    // Now read each property for each wind source and add the source
    // Need to read from disk to temp variables first.
    //
    list<pm_base *>::iterator iter=windsrc.begin();
    //
    // Problem: for data analysis, we read the wind source from each output
    // in turn, so we can't simply add it every time without deleting all the
    // elements beforehand.
    //
    struct stellarwind_params *temp_wind=0;
    while (SWP.params.size()>0) {
      temp_wind = SWP.params.back();
      SWP.params.pop_back();
      mem.myfree(temp_wind);
    }


    for (int isw=0; isw<SWP.Nsources; isw++) {
      // double Mdot, rad, posn[MAX_DIM], Vinf, Tw, trcr[MAX_NVAR], Rstar;
      // int type;
      ostringstream nm;
      struct stellarwind_params *wind=0;
      wind = mem.myalloc(wind,1);
      wind->id = isw;

      nm.str(""); nm << "WIND_"<<isw<<"_posn";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(posn));
      (*iter)->set_ptr(static_cast<void *>(wind->dpos));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_rad_";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&rad));
      (*iter)->set_ptr(static_cast<void *>(&wind->radius));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_Mdot";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&Mdot));
      (*iter)->set_ptr(static_cast<void *>(&wind->Mdot));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_Vinf";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&Vinf));
      (*iter)->set_ptr(static_cast<void *>(&wind->Vinf));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_Tw__";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&Tw));
      (*iter)->set_ptr(static_cast<void *>(&wind->Tstar));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_Rstr";
      if ( (*iter)->name.compare( nm.str() ) !=0)
        rep.error("Stellar wind parameters not ordered as expected!",
      (*iter)->name);
      //(*iter)->set_ptr(static_cast<void *>(&Rstar));
      (*iter)->set_ptr(static_cast<void *>(&wind->Rstar));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_type";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->type));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm.str()<<" = "<<wind->type<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_trcr";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //(*iter)->set_ptr(static_cast<void *>(trcr));
      (*iter)->set_ptr(static_cast<void *>(wind->tr));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      //
      // zero extra variables to avoid uninitialised data...
      // tracers are stored in the first ntracer values of array.
      //
      for (int v=SimPM.ntracer; v<MAX_NVAR; v++) wind->tr[v]=0.0;
      ++iter;
      //cout<<nm<<"\n";

      //
      // New stuff for evolving winds: data-file, time-offset, update-frequency.
      //
      nm.str(""); nm << "WIND_"<<isw<<"_evofile";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->evolving_wind_file));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm.str()<<" = "<<wind->evolving_wind_file<<"\n";

      //
      // Whether to enhance Mdot based on Omega (over and above what
      // the evolutionary code does).
      //
      nm.str(""); nm << "WIND_"<<isw<<"_enhance_mdot";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->enhance_mdot));
      err = read_header_param(*iter);
      if (err) {
        cout <<"failed to find WIND_"<<isw<<"_enhance_mdot parameter, setting to 0.\n";
        //rep.error("Error reading parameter",(*iter)->name);
      }
      ++iter;
      //cout<<nm.str()<<" = "<<wind->enhance_mdot<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_t_offset";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->time_offset));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm.str()<<" = "<<wind->time_offset<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_updatefreq";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->update_freq));
      err = read_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm.str()<<" = "<<wind->update_freq<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_t_scalefac";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&wind->t_scalefactor));
      err = read_header_param(*iter);
      if (err) {
        cout <<"Error reading parameter "<<(*iter)->name;
        cout <<" setting to default value of 1.0.\n";
        wind->t_scalefactor = 1.0;
        err=0;
      }
      ++iter;
      //cout<<nm.str()<<" = "<<wind->t_scalefactor<<"\n";

      //
      // Now we should have got all the sources, so add the source to 
      // the global list.
      //
      // NOTE THIS WIND STUCT DOESN'T GET FREED AT THE END OF THE SIMULATION,
      // SO IT IS TECHNICALLY A MEMORY LEAK.  I CAN'T THINK OF A GOOD WAY TO
      // GET RID OF IT THOUGH... IF I DELETE A DATAIO CLASS I DON'T WANT TO 
      // DELETE THE SWP DATA B/C SOMETIMES I'M NOT FINISHING THE SIMULATION.
      //
      SWP.params.push_back(wind);

    } // loop over sources
    
    //
    // Check we got all the sources:
    //
    // if (SWP.Nsources != SW.Nsources())
    //   rep.error("Got the wrong number of stellar wind sources!",
    //     SWP.Nsources-SW.Nsources());
    if (SWP.Nsources != static_cast<int>(SWP.params.size())) {
      cout <<"Num wind srcs="<<SWP.Nsources<<", SWP.params.size()="<<SWP.params.size()<<"\n";
      rep.error("Got the wrong number of stellar wind sources!",
                SWP.Nsources-static_cast<int>(SWP.params.size()));
    }
  } // read WIND data.


  //
  // Finally run some checks to make sure parameters are sane
  //
  err += check_header_parameters(SimPM);

  return err;
}  



// ##################################################################
// ##################################################################


void DataIOBase::set_tracer_params(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (have_setup_tracers) {
    rep.error("Trying to setup tracer parameters twice!",
              have_setup_tracers);
  }
  if (SimPM.ntracer==0) {
    //cout <<"\t *** No tracers to set up\n";
    have_setup_tracers = true;
    return;
  }
  //
  // So now we have tracers and we need to add them to the list.
  //
  if (!SimPM.tracers) {
    SimPM.tracers = mem.myalloc(SimPM.tracers,SimPM.ntracer);
  }

  class pm_base *p;

  for (int i=0;i<SimPM.ntracer;i++) {
    ostringstream temp; temp <<"Tracer";
    temp.width(3); temp.fill('0'); temp <<i;
    //cout <<"i="<<i<<", setting up tracer : "<<temp.str()<<"\n";

    pm_string  *ptemp = new pm_string  
      (temp.str(), &SimPM.tracers[i], "NEED_TRACER_VALUES!");
    p = ptemp; p->critical=true;  
    tr_pm.push_back(p);
  }
  have_setup_tracers = true;
  return;
}


// ##################################################################
// ##################################################################


void DataIOBase::set_windsrc_params()
{
  //
  // Sanity checks!
  //
  if (have_setup_windsrc || !windsrc.empty() || !SWP.Nsources) {
    cout <<"FLAG: "<< have_setup_windsrc;
    cout <<" EMPTY?: "<< windsrc.empty();
    cout <<" Wind-SIM?: "<< SWP.Nsources<<"\n";
    rep.error("WARNING! why set wind parameters?!",have_setup_windsrc);
  }
   
  //
  // Loop over sources and add strength+position for each
  //
  for (int n=0; n<SWP.Nsources; n++) {
    ostringstream temp1; temp1.str(""); temp1 << "WIND_"<<n<<"_posn";
    ostringstream temp2; temp2.str(""); temp2 << "WIND_"<<n<<"_rad_";
    ostringstream temp3; temp3.str(""); temp3 << "WIND_"<<n<<"_Mdot";
    ostringstream temp4; temp4.str(""); temp4 << "WIND_"<<n<<"_Vinf";
    ostringstream temp5; temp5.str(""); temp5 << "WIND_"<<n<<"_Tw__";
    ostringstream temp8; temp8.str(""); temp8 << "WIND_"<<n<<"_Rstr";
    ostringstream temp6; temp6.str(""); temp6 << "WIND_"<<n<<"_type";
    ostringstream temp7; temp7.str(""); temp7 << "WIND_"<<n<<"_trcr";
    //
    // New stuff for evolving winds:
    //
    ostringstream temp11; temp11.str(""); temp11<< "WIND_"<<n<<"_evofile";
    ostringstream temp13; temp13.str(""); temp13<< "WIND_"<<n<<"_enhance_mdot";
    ostringstream temp9;  temp9.str("");  temp9 << "WIND_"<<n<<"_t_offset";
    ostringstream temp10; temp10.str(""); temp10<< "WIND_"<<n<<"_updatefreq";
    ostringstream temp12; temp12.str(""); temp12<< "WIND_"<<n<<"_t_scalefac";

    
    pm_ddimarr *w001 = new pm_ddimarr (temp1.str()); // position of source.
    windsrc.push_back(w001);
    pm_double  *w002 = new pm_double  (temp2.str()); // radius of wind BC (cm)
    windsrc.push_back(w002);
    pm_double  *w003 = new pm_double  (temp3.str()); // Mdot (Msun/yr)
    windsrc.push_back(w003);
    pm_double  *w004 = new pm_double  (temp4.str()); // v_inf (cm/s)
    windsrc.push_back(w004);
    pm_double  *w005 = new pm_double  (temp5.str()); // Twind (K)
    windsrc.push_back(w005);
    pm_double  *w008 = new pm_double  (temp8.str()); // radius at which T=Twind
    windsrc.push_back(w008);
    pm_int     *w006 = new pm_int     (temp6.str()); // wind type flag
    windsrc.push_back(w006);
    pm_dvararr *w007 = new pm_dvararr (temp7.str()); // tracers
    windsrc.push_back(w007);
    //
    // New stuff for evolving winds:
    //
    // wind-evolution file.
    pm_string  *w011 = new pm_string  (temp11.str());
    w011->critical=false;
    windsrc.push_back(w011);

    // enhance mdot based on rotation?
    pm_int  *w013 = new pm_int        (temp13.str());
    w013->critical=false;
    windsrc.push_back(w013);

    // time offset
    pm_double  *w009 = new pm_double  (temp9.str());
    double dv=0.0;
    w009->critical=false; w009->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w009);

    // update frequency (in years)
    pm_double  *w010 = new pm_double  (temp10.str());
    w010->critical=false; dv=1000.0; w010->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w010);

    // scale factor (default must be 1, parameter must not be critical).
    pm_double  *w012 = new pm_double  (temp12.str());
    w012->critical=false; dv=1.0; w012->set_default_val(static_cast<void *>(&dv));
    windsrc.push_back(w012);
  }
  have_setup_windsrc=true;
  return;
}



// ##################################################################
// ##################################################################


void DataIOBase::set_rt_src_params(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //
  // Sanity checks!
  //
  if (have_setup_rt_src || !rt_src.empty() || !SimPM.RS.Nsources) {
    cout <<"FLAG: "<< have_setup_rt_src;
    cout <<" EMPTY?: "<< rt_src.empty();
    cout <<" RT-SIM?: "<< SimPM.RS.Nsources<<"\n";
    rep.error("WARNING! why set RT parameters?!",have_setup_rt_src);
  }
  
  //
  // First make sure that SimPM.RS.sources has the required number of elements.
  // If we are reading data, this should be empty, and we add Nsources elements.
  // If we are writing, it should have been set up before, so we shouldn't need
  // to add any elements.
  //
  if (SimPM.RS.sources.empty()) {
    cout <<" DataIOBase::set_rt_src_params() setting up SimPM.RS.sources: ";
    cout <<SimPM.RS.sources.size()<<"  "<<SimPM.RS.Nsources<<"\n";
    for (int n=0; n<SimPM.RS.Nsources; n++) {
      struct rad_src_info temp;

      for (int v=0;v<MAX_DIM; v++) temp.pos[v] = -1.0e200;
      temp.strength=-1.e200;
      temp.Rstar=-1.0e200;
      temp.Tstar=1.0e200;
      temp.id=n;
      temp.type=-1;
      temp.update=-1;
      temp.at_infinity=-1;
      temp.effect = -1;
      temp.NTau = 1;
      temp.opacity_src=-1;
      temp.opacity_var=-1;

      SimPM.RS.sources.push_back(temp);
      if (SimPM.RS.sources.size() != static_cast<unsigned int>(n+1)) {
        rep.error("Radiation source list size error, DataIOBase::set_rt_src_params",n);
      }
    }
  }
  if (SimPM.RS.sources.size() != static_cast<unsigned int>(SimPM.RS.Nsources)) {
    rep.error("wrong no. of srcs.",SimPM.RS.sources.size() );
  }

  //
  // Loop over sources and add strength, position, type, and (bool) location
  // for each
  //
#ifdef RT_TESTING
  cout <<"DataIOBase::set_rt_src_params() Nsrc="<<SimPM.RS.Nsources<<"\n";
#endif // RT_TESTING
  for (int n=0; n<SimPM.RS.Nsources; n++) {
    ostringstream temp2; temp2.str(""); temp2 << "RT_position_" <<n<<"_";
    ostringstream temp3; temp3.str(""); temp3 << "RT_strength_" <<n;
    ostringstream temp4; temp4.str(""); temp4 << "RT_src_type_" <<n;
    ostringstream temp5; temp5.str(""); temp5 << "RT_at_infty_" <<n;
    ostringstream temp6; temp6.str(""); temp6 << "RT_update___" <<n;
    ostringstream temp7; temp7.str(""); temp7 << "RT_Tau_src__" <<n;
    ostringstream temp8; temp8.str(""); temp8 << "RT_Tau_var__" <<n;
    ostringstream temp9; temp9.str(""); temp9 << "RT_effect___" <<n;
    ostringstream tmp10; tmp10.str(""); tmp10 << "RT_Rstar____" <<n;
    ostringstream tmp11; tmp11.str(""); tmp11 << "RT_Tstar____" <<n;
    ostringstream tmp12; tmp12.str(""); tmp12 << "RT_EVO_FILE_" <<n;

//ADD SOURCE_EFFECT VARIABLE!
    pm_ddimarr *rtpos = new pm_ddimarr (temp2.str(), (SimPM.RS.sources[n].pos));
    rt_src.push_back(rtpos);
    pm_double  *rtstr = new pm_double  (temp3.str(), &(SimPM.RS.sources[n].strength));
    rt_src.push_back(rtstr);
    pm_int     *rttyp = new pm_int     (temp4.str(), &(SimPM.RS.sources[n].type));
    rt_src.push_back(rttyp);
    pm_int     *rtinf = new pm_int     (temp5.str(), &(SimPM.RS.sources[n].at_infinity));
    rt_src.push_back(rtinf);
    pm_int     *rtupd = new pm_int     (temp6.str(), &(SimPM.RS.sources[n].update));
    rt_src.push_back(rtupd);
    pm_int     *rttsc = new pm_int     (temp7.str(), &(SimPM.RS.sources[n].opacity_src));
    rt_src.push_back(rttsc);
    pm_int     *rttvr = new pm_int     (temp8.str(), &(SimPM.RS.sources[n].opacity_var));
    rt_src.push_back(rttvr);
    pm_int     *rteff = new pm_int     (temp9.str(), &(SimPM.RS.sources[n].effect));
    rt_src.push_back(rteff);
    pm_double  *rtRst = new pm_double  (tmp10.str(), &(SimPM.RS.sources[n].Rstar), 0.0);
    rtRst->critical=false;
    rt_src.push_back(rtRst);
    pm_double  *rtTst = new pm_double  (tmp11.str(), &(SimPM.RS.sources[n].Tstar), 0.0);
    rtTst->critical=false;
    rt_src.push_back(rtTst);
    pm_string  *rtEvo = new pm_string  (tmp12.str(), &(SimPM.RS.sources[n].EvoFile), "NOFILE");
    rtEvo->critical=false;
    rt_src.push_back(rtEvo);
  }
  have_setup_rt_src=true;
  return;
}



// ##################################################################
// ##################################################################



void DataIOBase::set_jet_pm_params()
{
  if (have_setup_jet_pm || !jet_pm.empty() || !JP.jetic) {
    cout <<"FLAG: "<< have_setup_jet_pm;
    cout <<" EMPTY?: "<< jet_pm.empty();
    cout <<" JET-SIM?: "<< JP.jetic<<"\n";
    rep.error("WARNING! why set jet parameters?!",have_setup_jet_pm);
  }
  
  pm_int     *p051 = new pm_int     
    ("JetRadius",  &JP.jetradius);
  jet_pm.push_back(p051);
  pm_dvararr *p052 = new pm_dvararr 
    ("JetState",    JP.jetstate);
  jet_pm.push_back(p052);

  have_setup_jet_pm=true;
  return;
}


// ##################################################################
// ##################################################################


void DataIOBase::set_bc_pm_params(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  //cout <<"setting up BC parameters\n";
  if (have_setup_bc_pm || !bc_pm.empty()) {
    cout <<"BCs FLAG: "<< have_setup_bc_pm;
    cout <<" EMPTY?: "<< bc_pm.empty();
    cout <<"\n";
    rep.error("set bc parameters twice?!",have_setup_bc_pm);
  }
  
  pm_string  *p007 = new pm_string  
    ("BC_XN",&SimPM.BC_XN);
  p007->critical=true;  
  bc_pm.push_back(p007);
  pm_string  *p117 = new pm_string  
    ("BC_XP",&SimPM.BC_XP);
  p117->critical=true;  
  bc_pm.push_back(p117);
  pm_string  *p118 = new pm_string  
    ("BC_YN",&SimPM.BC_YN);
  //p118->critical=true;  
  bc_pm.push_back(p118);
  pm_string  *p119 = new pm_string  
    ("BC_YP",&SimPM.BC_YP);
  //p119->critical=true;  
  bc_pm.push_back(p119);
  pm_string  *p122 = new pm_string  
    ("BC_ZN",&SimPM.BC_ZN);
  //p122->critical=true;  
  bc_pm.push_back(p122);
  pm_string  *p123 = new pm_string  
    ("BC_ZP",&SimPM.BC_ZP);
  //p123->critical=true;  
  bc_pm.push_back(p123);

  //
  // Read internal boundary data
  //
  //SimPM.BC_INT.clear();
  //SimPM.BC_INT.resize(SimPM.BC_Nint);
  if (SimPM.BC_Nint>0) {
    if (!SimPM.BC_INT)
      SimPM.BC_INT = mem.myalloc(SimPM.BC_INT,SimPM.BC_Nint);
  }
  //cout <<"BC_Nint = "<<SimPM.BC_Nint<<"\n";
  for (int v=0; v<SimPM.BC_Nint; v++) {
    //cout <<"reading internal boundary "<<v<<"\n";
    ostringstream intbc; intbc.str("");
    intbc << "BC_INTERNAL_";
    intbc.width(3); intbc.fill('0');
    intbc << v;
    //cout <<"v="<<v<<", setting up Internal BC : "<<intbc.str()<<"\n";
    pm_string  *p124 = new pm_string (intbc.str(),&(SimPM.BC_INT[v]));
    //p124->critical=true;  
    bc_pm.push_back(p124);
  }
  have_setup_bc_pm=true;
  return;
}





// ##################################################################
// ##################################################################



int DataIOBase::write_simulation_parameters(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  if (params.empty())
    rep.error("Parameter list is empty -- make sure it populates itself!!",0);

  //
  // loop over all *normal* parameters, writing one-by-one.
  //
  pm_base *p=0;
  int err=0;
  for (list<pm_base *>::iterator iter=params.begin(); iter!=params.end(); ++iter) {
    p = (*iter);
    err = write_header_param(p);
    if (err) rep.error("Error writing parameter",(*iter)->name);
  }

  //
  // Write Boundary string parameters
  //
  if (!have_setup_bc_pm) set_bc_pm_params(SimPM);
  if (bc_pm.empty()) rep.error("bc_pm empty, Need boundaries!","Huh?");
  //
  // now write them:
  //
  for (list<pm_base *>::iterator iter=bc_pm.begin(); iter!=bc_pm.end(); ++iter) {
    p = (*iter);
    //cout <<p->name<<"    "; p->show_val(); cout <<"\n";
    err = write_header_param(p);
    if (err) rep.error("Error writing BC parameter",p->name);
  }

  //
  // Write tracer parameters
  //
  if (!have_setup_tracers) set_tracer_params(SimPM);
  //cout <<"Writing tracer names.\n";
  for (list<pm_base *>::iterator iter=tr_pm.begin(); iter!=tr_pm.end(); ++iter) {
    p = (*iter);
    //cout <<"tracer val: "; p->show_val(); cout <<"\n";
    err = write_header_param(p);
    if (err) rep.error("Error writing tracer parameter",p->name);
  }

  //
  // Write Jet parameters if doing a JET SIM
  //
  if (JP.jetic) {
    if (!have_setup_jet_pm) set_jet_pm_params();
    for (list<pm_base *>::iterator iter=jet_pm.begin(); iter!=jet_pm.end(); ++iter) {
      p = (*iter);
      err = write_header_param(p);
      if (err) rep.error("Error writing Jet parameter",p->name);
    }
  }
  
  //
  // Now the radiation sources
  //
  if (SimPM.RS.Nsources > 0) {
#ifdef TESTING
    cout <<"WRITING Radiation Source PARAMETERS, Nsrc="<<SimPM.RS.Nsources<<"\n";
#endif
    //
    // Check we have rt_src parameters
    //
    if (!have_setup_rt_src) set_rt_src_params(SimPM);
    if (rt_src.empty()) {
      rep.error("rt_src list empty, but running RT sim!?","HMMM");
    }
    //
    // Data should be already set up and ready to go.
    //
    //cout <<"rt_src.size = "<<rt_src.size()<<"\n";
    for (list<pm_base *>::iterator iter=rt_src.begin(); iter!=rt_src.end(); ++iter) {
      p = (*iter);
      err = write_header_param(p);
      if (err) rep.error("Error writing Radiation source parameter",p->name);
    }
  } // write RT data.
    
  // --------------------------------------
  // ---- Now the stellar wind sources ----
  // --------------------------------------
  if (SWP.Nsources > 0) {
    //
    // Check we have windsrc parameters
    //
    if (!have_setup_windsrc) set_windsrc_params();
    if (windsrc.empty())
      rep.error("windsrc list empty, but we are apparently ouputting \
                 a wind sim!?","HMMM");

    //
    // Now write each property from WS class.
    // Need to copy to a temp variable for the silo interface.
    // List elements are ordered src 0,1,2.
    //
    list<pm_base *>::iterator iter=windsrc.begin();

    for (int isw=0; isw<SWP.Nsources; isw++) {
      //double xd[MAX_NVAR];
      //int xi[MAX_NVAR];
      ostringstream nm;
      //
      // Make sure the pm name is right, then set the pointer and write the
      // data for each property of source number isw
      //
      nm.str(""); nm << "WIND_"<<isw<<"_posn";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_posn(isw,xd); 
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->dpos));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_rad_";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->radius));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_Mdot";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_Mdot(isw,xd); 
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Mdot));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_Vinf";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_Vinf(isw,xd); 
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Vinf));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_Tw__";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_Tw(isw,xd); 
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Tstar));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_Rstr";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_Rstar(isw,xd); 
      //(*iter)->set_ptr(static_cast<void *>(xd));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->Rstar));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_type";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      //SW.get_src_type(isw,xi); 
      //(*iter)->set_ptr(static_cast<void *>(xi));
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->type));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;

      nm.str(""); nm << "WIND_"<<isw<<"_trcr";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(SWP.params[isw]->tr));
      for (int v=SimPM.ntracer; v<MAX_NVAR; v++) {
        SWP.params[isw]->tr[v]=0.0;
      }
      err = write_header_param(*iter);
      if (err) rep.error("Error writing RT parameter",(*iter)->name);
      ++iter;
      //
      // New stuff for evolving winds: data-file, time-offset, update-frequency.
      //
      nm.str(""); nm << "WIND_"<<isw<<"_evofile";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->evolving_wind_file));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<" = "<<wind->evolving_wind_file<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_enhance_mdot";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->enhance_mdot));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<" = "<<wind->enhance_mdot<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_t_offset";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->time_offset));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<" = "<<wind->time_offset<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_updatefreq";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->update_freq));
      err = write_header_param(*iter);
      if (err) rep.error("Error reading parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<" = "<<wind->update_freq<<"\n";

      nm.str(""); nm << "WIND_"<<isw<<"_t_scalefac";
      if ( (*iter)->name.compare( nm.str() ) !=0) {
        rep.error("Stellar wind parameters not ordered as expected!",(*iter)->name);
      }
      (*iter)->set_ptr(static_cast<void *>(&SWP.params[isw]->t_scalefactor));
      err = write_header_param(*iter);
      if (err) rep.error("Error writing parameter",(*iter)->name);
      ++iter;
      //cout<<nm<<" = "<<wind->t_scalefactor<<"\n";

    } // loop over sources
  } // write Stellar wind data.
    
  return err;
}



// ##################################################################
// ##################################################################


int DataIOBase::check_header_parameters(
      class SimParams &SimPM  ///< pointer to simulation parameters
      )
{
  // This is where we check that nvar,ndim,ntracer,outfile, etc. are all set to
  // something sensible.
  if (SimPM.simtime <0) {
    rep.warning("(Dataio::check_header_parameters) simtime not set, setting to zero. THIS MAY BE A BIG ERROR",1.,SimPM.simtime);
    SimPM.simtime = 0.;
  }
  if (SimPM.starttime <0) {
    rep.warning("(Dataio::check_header_parameters) startime <0, PROBABLY AN ERROR!",0.,SimPM.starttime);
    SimPM.starttime = 0.;
  }
  if (SimPM.finishtime <=0 ) {rep.warning("(Dataio::check_header_parameters) finishtime not set.",1.,SimPM.finishtime);SimPM.finishtime=-1.;}
  if (SimPM.timestep <0) {
    rep.warning("(Dataio::check_header_parameters) timestep <0, PROBABLY AN ERROR!",1,SimPM.timestep);
    SimPM.timestep = 0;
  }
  if (SimPM.NG[0]<0) rep.error("(Dataio::check_header_parameters) NG[0]<0 -- must not have read from file",SimPM.NG[0]);
  if (SimPM.ndim>1) if (SimPM.NG[1]<0) rep.error("(Dataio::check_header_parameters) NG[1]<0 -- must not have read from file",SimPM.NG[1]);
  if (SimPM.ndim>2) if (SimPM.NG[2]<0) rep.error("(Dataio::check_header_parameters) NG[2]<0 -- must not have read from file",SimPM.NG[2]);

  if (SimPM.Ncell<0)   rep.error("(Dataio::check_header_parameters) Ncell<0 -- Error reading from file",SimPM.Ncell);
 
  if (SimPM.Xmin[0] <-1.e50)  rep.error("(Dataio::check_header_parameters) Xmin[0]<0 -- Error reading from file",SimPM.Xmin[0]);
  if (SimPM.ndim>1) if (SimPM.Xmin[1]<-1e50) rep.error("(Dataio::check_header_parameters) Xmin[1]<0 -- must not have read from file",SimPM.Xmin[1]);
  if (SimPM.ndim>2) if (SimPM.Xmin[2]<-1e50) rep.error("(Dataio::check_header_parameters) Xmin[2]<0 -- must not have read from file",SimPM.Xmin[2]);

  if (SimPM.Xmax[0] <=SimPM.Xmin[0])  rep.error("(Dataio::check_header_parameters) Xmax[0]<0 -- Error reading from file",SimPM.Xmax[0]);
  if (SimPM.ndim>1) if (SimPM.Xmax[1]<=SimPM.Xmin[1]) rep.error("(Dataio::check_header_parameters) Xmax[1]<0 -- must not have read from file",SimPM.Xmax[1]);
  if (SimPM.ndim>2) if (SimPM.Xmax[2]<=SimPM.Xmin[2]) rep.error("(Dataio::check_header_parameters) Xmax[2]<0 -- must not have read from file",SimPM.Xmax[2]);


  if (SimPM.spOOA!=OA1 && SimPM.spOOA!=OA2) {
    rep.warning("(Dataio::check_header_parameters) Order of Accuracy not set from file. Default to second order",2,SimPM.spOOA);
    SimPM.spOOA=OA2;
  }
  if (SimPM.tmOOA!=OA1 && SimPM.tmOOA!=OA2) {
    rep.warning("(Dataio::check_header_parameters) Order of Accuracy not set from file. Default to second order",2,SimPM.tmOOA);
    SimPM.tmOOA=OA2;
  }
  
  //  cout <<"(Dataio::check_header_parameters) spOOA: "<<SimPM.spOOA<<" tmOOA: "<<SimPM.tmOOA<<"\n";
  if (SimPM.gamma<0)
    rep.error("(Dataio::check_header_parameters) gamma<0 -- must not have read it from file",SimPM.gamma);
  if (SimPM.CFL<0)
    rep.error("(Dataio::check_header_parameters) CFL<0   -- must not have read it from file",SimPM.CFL);
  if (SimPM.artviscosity<0)
    rep.error("(Dataio::check_header_parameters) artviscosity<0 -- must not have read from file",SimPM.artviscosity);
  if (SimPM.etav<0 && SimPM.artviscosity>0)
    rep.error("(Dataio::check_header_parameters) etav<0 -- must not have read from file",SimPM.etav);
  else if (SimPM.etav<0 && SimPM.artviscosity<0) {
    rep.warning("(Dataio::check_header_parameters) etav<0 but art.visc. not enabled so it's ok.",
    static_cast<double>(SimPM.artviscosity),SimPM.etav);
    SimPM.etav=0.15;
  }
  if (SimPM.typeofop <0) {
    rep.warning("(Dataio::check_header_parameters) type of output file not set, setting to Fits",2,SimPM.typeofop);
    SimPM.typeofop = 2;
  }

  if (SimPM.opfreq <0) {
    rep.warning("(Dataio::check_header_parameters) Output frequency not set.  Setting to 50",100,SimPM.opfreq);
    SimPM.opfreq = 50;
  }
  if (SimPM.outFileBase =="") {
    rep.warning("(Dataio::check_header_parameters) output file base not set, setting to ./noname",1,1);
    SimPM.outFileBase = "noname";
  }

#ifdef TESTING
  dp.initERG = dp.initMMX = dp.initMMY = dp.initMMZ = 0.;
  dp.ergTotChange = dp.mmxTotChange = dp.mmyTotChange = dp.mmzTotChange = 0.;
#endif //TESTING
  SimPM.maxtime = false;
  SimPM.dt = 0.;

  if (SimPM.EP.dynamics!=1) {
    cout <<"WARNING: DYNAMICS SWITCHED OFF!!! *****************************\n";
  }

  //
  // Check that Helium_MassFrac and H_MassFrac are set.
  //
  if (SimPM.EP.H_MassFrac <0.001  || SimPM.EP.H_MassFrac>1.0) {
    cout <<"Warning: H_MassFrac = "<<SimPM.EP.H_MassFrac<<", ";
    cout <<"resetting to 0.7154\n";
    SimPM.EP.H_MassFrac = 0.7154;
  }
  if ((SimPM.EP.H_MassFrac+SimPM.EP.Helium_MassFrac)>1.1 ||
      (SimPM.EP.H_MassFrac+SimPM.EP.Helium_MassFrac)<0.9) {
    cout <<"Warning: H_MassFrac  = "<<SimPM.EP.H_MassFrac<<"\n";
    cout <<"Warning: He_MassFrac = "<<SimPM.EP.H_MassFrac<<"\n";
    cout <<"resetting to 0.7154 and 0.2846, respectively.\n";
    SimPM.EP.H_MassFrac = 0.7154;
    SimPM.EP.Helium_MassFrac = 0.2846;
  }
  if (SimPM.EP.Metal_MassFrac<0.0 || SimPM.EP.Metal_MassFrac>1.0) {
    cout <<"Warning: SimPM.EP.Metal_MassFrac = ";
    cout <<SimPM.EP.Metal_MassFrac<<" !! resetting to 0.0142.\n";
    SimPM.EP.Metal_MassFrac= 0.0142;
  }


  return 0;
}

// -----------------------------------------------------
// --------- BASE DATA I/O CLASS DEFINITIONS -----------
// -----------------------------------------------------


// ##################################################################
// ##################################################################

// ##################################################################
// ##################################################################



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
      class RayTracingBase *RT, ///< pointer to raytracing class
      const long int counter   ///< number to stamp file with (e.g. timestep)
      )
{
  if (!cg[0])
    rep.error("dataio_text::output_ascii_data() null pointer to grid!",cg[0]);
  dataio_text::gp = cg[0];

  cout <<"dataio_text::OutputData() writing data.\n";
  string fname = set_filename(outfile, counter);
  int err = output_ascii_data(fname, SimPM,RT);
  cout <<"dataio_text::OutputData() written data.\n";
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
      class SimParams &SimPM,  ///< pointer to simulation parameters
      class RayTracingBase *RT ///< pointer to raytracing class [OUTPUT]
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
    if (RT) {
      for (int v=0;v<SimPM.RS.Nsources;v++) {
        //cout <<"hello";
        CI.get_col(cpt, v, Utemp);
        for (int iT=0; iT<SimPM.RS.sources[v].NTau; iT++)
          outf <<"  "<< Utemp[iT];
      }
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



