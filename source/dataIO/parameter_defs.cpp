/// \file parameter_defs.cpp
/// \author Jonathan Mackey
///
/// Class definitions for simulation parameter types.
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
///    to read and write.  read/write_simulation_parameters() functions work
///    now.
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
/// - 2011.02.15 JM: Added new stellar wind parameters for evolving wind
/// sources.
/// - 2011.02.24 JM: Added read/write for multiple radiation sources, with
///    additional parameters. Much simplified RT I/O by using new struct
///    SimPM.RS.
/// - 2011.02.25 JM: got rid of HCORR ifdef wrapper.
/// - 2011.02.28 JM: got rid of RSP references. 03.01 few bug fixes for new RT.
/// - 2011.03.02 JM: Better support for tracer variables (with or without
/// names).
///        Can handle an arbitrary number of tracers now.
/// - 2011.03.21 JM: Added extra radiation source parameters. (and 22.03)
/// - 2011.05.02 JM: New RT params.
/// - 2011.06.02 JM: WriteHeader() added; some minor mods to dataio_text.
/// - 2011.10.05 JM: Fixed problems in the wind-source reading (for analysis
/// code).
/// - 2011.11.22 JM: Added Stellar wind t_scalefactor parameter.
/// - 2011.12.14 JM: Added radiation source EvoFile parameter (time-varying
/// sources).
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
/// - 2018.05.01 JM: moved DataIOBase and dataio_text to other files.

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/mem_manage.h"
#include "tools/reporting.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif  // TESTING

#include "dataIO/dataio_base.h"
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
{
    type   = MY_INT;
    defval = -1;
}
pm_double::pm_double()
{
    type   = MY_DOUBLE;
    defval = -VERY_LARGE_VALUE;
}
pm_float::pm_float()
{
    type   = MY_FLOAT;
    defval = -SMALLVALUE;
}
pm_long::pm_long()
{
    type   = MY_LONG;
    defval = -1;
}
pm_string::pm_string()
{
    type = MY_STRING;
    defval.assign("hello");
}
pm_ddimarr::pm_ddimarr()
{
    type   = MY_DDIMARR;
    len    = MAX_DIM;
    defval = 0;
}
pm_idimarr::pm_idimarr()
{
    type   = MY_IDIMARR;
    len    = MAX_DIM;
    defval = 0;
}
pm_dvararr::pm_dvararr()
{
    type   = MY_DVARARR;
    len    = MAX_NVAR;
    defval = 0;
}

// ##################################################################
// ##################################################################

//
// constructor with name only
//
pm_int::pm_int(const string s)
{
    type = MY_INT;
    name.assign(s);
    defval = -1;
}
pm_double::pm_double(const string s)
{
    type = MY_DOUBLE;
    name.assign(s);
    defval = -VERY_LARGE_VALUE;
}
pm_float::pm_float(const string s)
{
    type = MY_FLOAT;
    name.assign(s);
    defval = -SMALLVALUE;
}
pm_long::pm_long(const string s)
{
    type = MY_LONG;
    name.assign(s);
    defval = -1;
}
pm_string::pm_string(const string s)
{
    type = MY_STRING;
    name.assign(s);
    defval.assign("HELLO");
}
pm_ddimarr::pm_ddimarr(const string s)
{
    type = MY_DDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    defval = 0;
}
pm_idimarr::pm_idimarr(const string s)
{
    type = MY_IDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    defval = 0;
}
pm_dvararr::pm_dvararr(const string s)
{
    type = MY_DVARARR;
    len  = MAX_NVAR;
    name.assign(s);
    defval = 0;
}

// ##################################################################
// ##################################################################

//
// constructor with name and pointer to data
//
pm_int::pm_int(const string s, int* p)
{
    type = MY_INT;
    name.assign(s);
    ptr    = p;
    defval = -1;
}
pm_double::pm_double(const string s, double* p)
{
    type = MY_DOUBLE;
    name.assign(s);
    ptr    = p;
    defval = -VERY_LARGE_VALUE;
}
pm_float::pm_float(const string s, float* p)
{
    type = MY_FLOAT;
    name.assign(s);
    ptr    = p;
    defval = -SMALLVALUE;
}
pm_long::pm_long(const string s, long int* p)
{
    type = MY_LONG;
    name.assign(s);
    ptr    = p;
    defval = -1;
}
pm_string::pm_string(const string s, string* p)
{
    type = MY_STRING;
    name.assign(s);
    ptr = p;
    defval.assign("HELLO");
}
pm_ddimarr::pm_ddimarr(const string s, double* p)
{
    type = MY_DDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    ptr    = p;
    defval = 0;
}
pm_idimarr::pm_idimarr(const string s, int* p)
{
    type = MY_IDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    ptr    = p;
    defval = 0;
}
pm_dvararr::pm_dvararr(const string s, pion_flt* p)
{
    type = MY_DVARARR;
    len  = MAX_NVAR;
    name.assign(s);
    ptr    = p;
    defval = 0;
}

// ##################################################################
// ##################################################################

//
// constructor with name, pointer to data, default value.
//
pm_int::pm_int(const string s, int* p, const int def)
{
    type = MY_INT;
    name.assign(s);
    ptr    = p;
    defval = def;
    // cout <<"PM_INT: "<<s<<".. "<<ptr<<".. "<<defval<<"\n";
}
pm_double::pm_double(const string s, double* p, const double def)
{
    type = MY_DOUBLE;
    name.assign(s);
    ptr    = p;
    defval = def;
}
pm_float::pm_float(const string s, float* p, const float def)
{
    type = MY_FLOAT;
    name.assign(s);
    ptr    = p;
    defval = def;
}
pm_long::pm_long(const string s, long int* p, const long int def)
{
    type = MY_LONG;
    name.assign(s);
    ptr    = p;
    defval = def;
}
pm_string::pm_string(const string s, string* p, const string def)
{
    type = MY_STRING;
    name.assign(s);
    ptr = p;
    defval.assign(def);
}
pm_ddimarr::pm_ddimarr(const string s, double* p, const double* def)
{
    type = MY_DDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    ptr    = p;
    defval = mem.myalloc(defval, len);
    for (int v = 0; v < len; v++)
        defval[v] = def[v];
}
pm_idimarr::pm_idimarr(const string s, int* p, const int* def)
{
    type = MY_IDIMARR;
    len  = MAX_DIM;
    name.assign(s);
    ptr    = p;
    defval = mem.myalloc(defval, len);
    for (int v = 0; v < len; v++)
        defval[v] = def[v];
}
pm_dvararr::pm_dvararr(const string s, pion_flt* p, const pion_flt* def)
{
    type = MY_DVARARR;
    len  = MAX_NVAR;
    name.assign(s);
    ptr    = p;
    defval = mem.myalloc(defval, len);
    for (int v = 0; v < len; v++)
        defval[v] = def[v];
}

// ##################################################################
// ##################################################################

//
// Some have destructors:
//
pm_ddimarr::~pm_ddimarr()
{
    if (defval) defval = mem.myfree(defval);
    ptr = 0;
}
pm_idimarr::~pm_idimarr()
{
    if (defval) defval = mem.myfree(defval);
    ptr = 0;
}
pm_dvararr::~pm_dvararr()
{
    if (defval) defval = mem.myfree(defval);
    ptr = 0;
}

// ##################################################################
// ##################################################################

void pm_int::assign_val(void* val)
{
    *ptr = *(static_cast<int*>(val));
}
void pm_double::assign_val(void* val)
{
    *ptr = *(static_cast<double*>(val));
}
void pm_float::assign_val(void* val)
{
    *ptr = *(static_cast<float*>(val));
}
void pm_long::assign_val(void* val)
{
    *ptr = *(static_cast<long int*>(val));
}
void pm_string::assign_val(void* val)
{
    (*ptr).assign(*(static_cast<string*>(val)));
}
void pm_ddimarr::assign_val(void* val)
{
    for (int i = 0; i < len; i++)
        ptr[i] = (static_cast<double*>(val))[i];
}
void pm_idimarr::assign_val(void* val)
{
    for (int i = 0; i < len; i++)
        ptr[i] = (static_cast<int*>(val))[i];
}
void pm_dvararr::assign_val(void* val)
{
    for (int i = 0; i < len; i++)
        ptr[i] = (static_cast<pion_flt*>(val))[i];
}

// ##################################################################
// ##################################################################

void pm_int::set_ptr(void* p)
{
    ptr = static_cast<int*>(p);
}
void pm_double::set_ptr(void* p)
{
    ptr = static_cast<double*>(p);
}
void pm_float::set_ptr(void* p)
{
    ptr = static_cast<float*>(p);
}
void pm_long::set_ptr(void* p)
{
    ptr = static_cast<long int*>(p);
}
void pm_string::set_ptr(void* p)
{
    ptr = static_cast<string*>(p);
}
void pm_ddimarr::set_ptr(void* p)
{
    ptr = static_cast<double*>(p);
}
void pm_idimarr::set_ptr(void* p)
{
    ptr = static_cast<int*>(p);
}
void pm_dvararr::set_ptr(void* p)
{
    ptr = static_cast<pion_flt*>(p);
}

// ##################################################################
// ##################################################################

void pm_int::show_val()
{
    cout << *ptr;
}
void pm_double::show_val()
{
    cout << *ptr;
}
void pm_float::show_val()
{
    cout << *ptr;
}
void pm_long::show_val()
{
    cout << *ptr;
}
void pm_string::show_val()
{
    cout << *ptr << ", name=" << name;
}
void pm_ddimarr::show_val()
{
    cout << "[";
    for (int i = 0; i < len - 1; i++)
        cout << ptr[i] << ", ";
    cout << ptr[len - 1] << "]";
}
void pm_idimarr::show_val()
{
    cout << "[";
    for (int i = 0; i < len - 1; i++)
        cout << ptr[i] << ", ";
    cout << ptr[len - 1] << "]";
}
void pm_dvararr::show_val()
{
    cout << "[";
    for (int i = 0; i < len - 1; i++)
        cout << ptr[i] << ", ";
    cout << ptr[len - 1] << "]";
}

// ##################################################################
// ##################################################################

void pm_int::set_to_default()
{
    *ptr = defval;
}
void pm_double::set_to_default()
{
    *ptr = defval;
}
void pm_float::set_to_default()
{
    *ptr = defval;
}
void pm_long::set_to_default()
{
    *ptr = defval;
}
void pm_string::set_to_default()
{
    (*ptr).assign(defval);
}
void pm_ddimarr::set_to_default()
{
    if (!defval)
        rep.error("No default value!", "SET ME!");
    else
        for (int i = 0; i < len; i++)
            ptr[i] = defval[i];
}
void pm_idimarr::set_to_default()
{
    if (!defval)
        rep.error("No default value!", "SET ME!");
    else
        for (int i = 0; i < len; i++)
            ptr[i] = defval[i];
}

void pm_dvararr::set_to_default()
{
    if (!defval)
        rep.error("No default value!", "SET ME!");
    else
        for (int i = 0; i < len; i++)
            ptr[i] = defval[i];
}

// ##################################################################
// ##################################################################

void* pm_int::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_double::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_float::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_long::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_string::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_ddimarr::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_idimarr::get_ptr()
{
    return static_cast<void*>(ptr);
}
void* pm_dvararr::get_ptr()
{
    return static_cast<void*>(ptr);
}

// ##################################################################
// ##################################################################

void pm_int::set_default_val(void* v)
{
    defval = *(static_cast<int*>(v));
}
void pm_double::set_default_val(void* v)
{
    defval = *(static_cast<double*>(v));
}
void pm_float::set_default_val(void* v)
{
    defval = *(static_cast<float*>(v));
}
void pm_long::set_default_val(void* v)
{
    defval = *(static_cast<long*>(v));
}
void pm_string::set_default_val(void* v)
{
    defval = *(static_cast<string*>(v));
}
void pm_ddimarr::set_default_val(void* v)
{
    if (!defval)
        defval = mem.myalloc(defval, len);
    else
        for (int i = 0; i < len; i++)
            defval[i] = (static_cast<double*>(v))[i];
}
void pm_idimarr::set_default_val(void* v)
{
    if (!defval)
        defval = mem.myalloc(defval, len);
    else
        for (int i = 0; i < len; i++)
            defval[i] = (static_cast<double*>(v))[i];
}
void pm_dvararr::set_default_val(void* v)
{
    if (!defval)
        defval = mem.myalloc(defval, len);
    else
        for (int i = 0; i < len; i++)
            defval[i] = (static_cast<pion_flt*>(v))[i];
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

// ##################################################################
// ##################################################################
