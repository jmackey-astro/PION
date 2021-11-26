/// \file parameter_defs.h
/// \author Jonathan Mackey
///
/// modified:\n
/// - 2007-10-25 Got class structure that works for parallel output.
/// - 2007-10-26 parallel fits-io class reads data from single and multiple
///    files. Put in serial ifdefs, so that it should work for serial code too.
/// - 2008-01-16 Documented the classes.
/// - 2010-02-03 JM: made base class virtual
/// - 2010-04-21 JM: Changed filename setup so that i can write
///    checkpoint files with fname.999999.txt/silo/fits
/// - 2101-07-20/22 JM: Updated sim.file header I/O to be more general
///    and make it easier to add new parameters.
/// - 2010.10.13 JM: Removed NEW_SOLVER_STRUCT ifdefs.
/// - 2010.11.03 JM: Added Ndigits variable for naming files.
/// - 2011.06.02 JM: Added WriteHeader() function so I can over-write header
///    parameters and restart a simulation with e.g. different microphysics.
/// - 2011.11.22 JM: Added set_default_val() function to parameter classes.
/// - 2015.07.09 JM: updated tracer parameters (each tracer has its
///    own parameter name.
/// - 2015.08.05 JM: tidied up code; added pion_flt datatype.
/// - 2015.10.19 JM: Fixed dvararr to always use pion_flt correctly.
/// - 2017.08.03 JM: Added parallel_read_any_data() function to base class
/// - 2018.01.24 JM: worked on making SimPM non-global
/// - 2018.05.01 JM: moved DataIOBase and dataio_text to other files.

#ifndef PARAMETER_DEFS_H
#define PARAMETER_DEFS_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"


#include <cstring>
#include <list>
#include <string>

// -----------------------------------------------------------
// -----------------------------------------------------------

//
// datatypes for reading/writing parameters
//
#define MY_INT 1
#define MY_DOUBLE 2
#define MY_STRING 3
#define MY_FLOAT 4
#define MY_LONG 5
#define MY_DDIMARR 6  // double 3d array (MAX_NDIM)
#define MY_IDIMARR 7  // integer 3d array (MAX_NDIM)
#define MY_DVARARR 8  // double Nd array (MAX_NVAR)

///
/// base parameter class, each parameter gets an instance of a derived
/// class with the base class interface.
///
class pm_base {
public:
  virtual ~pm_base() {}
  /// sets value of data pointed to by ptr to value in argument.
  virtual void assign_val(void *) = 0;
  virtual void set_to_default() = 0;  ///< sets value of data to default value.
  virtual void set_ptr(void *)  = 0;  ///< changes pointer to new address.
  virtual void show_val()       = 0;  ///< returns value pointed to by ptr.
  virtual void *get_ptr()       = 0;  ///< returns address of ptr.
  virtual void set_default_val(
      void *) = 0;   ///< sets default value to fn. argument.
  std::string name;  ///< string identifier for parameter.
  int type;          ///< type of data int,long,dbl,array,string,etc.
  bool critical;     ///< true if parameter is required.
};

///
/// parameter class for Integer parameters
///
class pm_int : public pm_base {
public:
  pm_int();
  pm_int(const std::string);
  pm_int(const std::string, int *);
  pm_int(
      const std::string,  ///< name
      int *,              ///< data pointer
      const int           ///< default value
  );
  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  int *ptr;
  int defval;
};

///
/// parameter class for Double parameters
///
class pm_double : public pm_base {
public:
  pm_double();
  pm_double(const std::string);
  pm_double(const std::string, double *);
  pm_double(const std::string, double *, const double);
  void assign_val(void *val);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  double *ptr;
  double defval;
};

///
/// parameter class for FLOAT parameters
///
class pm_float : public pm_base {
public:
  pm_float();
  pm_float(const std::string);
  pm_float(const std::string, float *);
  pm_float(const std::string, float *, const float);
  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  float *ptr;
  float defval;
};

///
/// parameter class for LONG-INT parameters
///
class pm_long : public pm_base {
public:
  pm_long();
  pm_long(const std::string);
  pm_long(const std::string, long int *);
  pm_long(const std::string, long int *, long int);
  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  long int *ptr;
  long int defval;
};

///
/// parameter class for STRING parameters
///
class pm_string : public pm_base {
public:
  pm_string();
  pm_string(const std::string);
  pm_string(const std::string, std::string *);
  pm_string(const std::string, std::string *, const std::string);
  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  std::string *ptr;
  std::string defval;
};

///
/// parameter class for DOUBLE ARRAY [MAX_DIM] parameters
///
class pm_ddimarr : public pm_base {
public:
  pm_ddimarr();
  pm_ddimarr(const std::string);
  pm_ddimarr(const std::string, double *);
  pm_ddimarr(const std::string, double *, const double *);
  ~pm_ddimarr();

  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  double *ptr;
  int len;
  double *defval;
};

///
/// parameter class for INT ARRAY [MAX_DIM] parameters
///
class pm_idimarr : public pm_base {
public:
  pm_idimarr();
  pm_idimarr(const std::string);
  pm_idimarr(const std::string, int *);
  pm_idimarr(const std::string, int *, const int *);
  ~pm_idimarr();

  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  int *ptr;
  int len;
  int *defval;
};

///
/// parameter class for DOUBLE ARRAY [MAX_NVAR] parameters
/// This uses float or double, depending on the value of pion_flt
/// in functionality_flags.h
///
class pm_dvararr : public pm_base {
public:
  pm_dvararr();
  pm_dvararr(const std::string);
  pm_dvararr(const std::string, pion_flt *);
  pm_dvararr(const std::string, pion_flt *, const pion_flt *);
  ~pm_dvararr();

  void assign_val(void *);
  void set_ptr(void *);
  void *get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);

protected:
  pion_flt *ptr;
  int len;
  pion_flt *defval;
};
// -----------------------------------------------------------
// -----------------------------------------------------------

#endif  // PARAMETER_DEFS_H
