/// \file dataio.h
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

#ifndef DATAIO_H
#define DATAIO_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"

#include <cstring>
#include <list>

#include "spatial_solvers/solver_eqn_base.h"
#include "dataIO/file_status.h"

// -----------------------------------------------------------
// -----------------------------------------------------------

//
// datatypes for reading/writing parameters
//
#define MY_INT    1
#define MY_DOUBLE 2
#define MY_STRING 3
#define MY_FLOAT  4
#define MY_LONG   5
#define MY_DDIMARR 6 // double 3d array (MAX_NDIM)
#define MY_IDIMARR 7 // integer 3d array (MAX_NDIM)
#define MY_DVARARR 8 // double Nd array (MAX_NVAR)

///
/// base parameter class, each parameter gets an instance of a derived
/// class with the base class interface.
///
class pm_base {
public:
  virtual ~pm_base() {}
   /// sets value of data pointed to by ptr to value in argument.
  virtual void assign_val(void *)=0;
  virtual void set_to_default()=0;   ///< sets value of data to default value.
  virtual void set_ptr(void *)=0;    ///< changes pointer to new address.
  virtual void show_val()=0;   ///< returns value pointed to by ptr.
  virtual void * get_ptr()=0;  ///< returns address of ptr.
  virtual void set_default_val(void *)=0; ///< sets default value to fn. argument.
  string name;   ///< string identifier for parameter.
  int type;      ///< type of data int,long,dbl,array,string,etc.
  bool critical; ///< true if parameter is required.
};

///
/// parameter class for Integer parameters
///
class pm_int : public pm_base {
public:
  pm_int();
  pm_int(const string );
  pm_int(const string, int *);
  pm_int(const string, ///< name
	 int *,        ///< data pointer
	 const int     ///< default value
	 );
  void assign_val(void *);
  void set_ptr(void *);
  void * get_ptr();
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
  pm_double(const string );
  pm_double(const string, double *);
  pm_double(const string,
	    double *,
	    const double
	    );
  void assign_val(void *val);
  void set_ptr(void *); 
  void * get_ptr();
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
  pm_float(const string );
  pm_float(const string, float *);
  pm_float(const string,
	   float *,
	   const float
	   );
  void assign_val(void *);
  void set_ptr(void *); 
  void * get_ptr();
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
  pm_long(const string );
  pm_long(const string, long int *);
  pm_long(const string ,
	  long int *,
	  long int
	  );
  void assign_val(void *);
  void set_ptr(void *); 
  void * get_ptr();
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
  pm_string(const string );
  pm_string(const string, string *);
  pm_string(const string,
	    string *,
	    const string
	    );
  void assign_val(void *);
  void set_ptr(void *); 
  void * get_ptr();
  void show_val();
  void set_to_default();
  void set_default_val(void *);
protected:
  string *ptr;
  string defval;
};

///
/// parameter class for DOUBLE ARRAY [MAX_DIM] parameters
///
class pm_ddimarr : public pm_base {
public:
  pm_ddimarr();
  pm_ddimarr(const string);
  pm_ddimarr(const string, double *);
  pm_ddimarr(const string, double *, const double *);
  ~pm_ddimarr();

  void assign_val(void *);
  void set_ptr(void *); 
  void * get_ptr();
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
  pm_idimarr(const string);
  pm_idimarr(const string, int *);
  pm_idimarr(const string, int *, const int *);
  ~pm_idimarr();

  void assign_val(void *);
  void set_ptr(void *); 
  void * get_ptr();
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
  pm_dvararr(const string );
  pm_dvararr(const string, pion_flt *);
  pm_dvararr(const string, pion_flt *, const pion_flt *);
  ~pm_dvararr();

  void assign_val(void *);
  void set_ptr(void *);
  void * get_ptr();
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

///
/// Abstract base class for Data I/O in my finite-volume code.
///
/// This is a virtual (or interface) class, defining the allowed
/// public functions.
/// It also defines the list of parameters to read and write in 
/// file headers.
///
class DataIOBase : public file_status {
  public:
  DataIOBase();
  ~DataIOBase();
   /** \brief Set a pointer to the solver. */
   virtual void SetSolver(FV_solver_base *
      ///< Pointer to the solver (to get Eint,divB,Ptot)
			  )=0;
  
  ///
  /// Write simulation header info to file
  ///
  virtual int WriteHeader(
          const string ///< file to write to (full, exact filename).
          )=0;
  

  ///
  /// Write simulation data to file.
  ///
  virtual int OutputData(
        const string, ///< File to write to
        class GridBaseClass *, ///< pointer to data.
        const long int ///< number to stamp file with (e.g. timestep)
        )=0;

   /** \brief Reads simulation parameters from file. */
   virtual int ReadHeader(string ///< file to read from
			   )=0;
   /** Having set up a grid with parameters from the header, this reads
    * data from the file and puts it on the grid.
    * */
   virtual int ReadData(string, ///< file to read from
			class GridBaseClass * ///< pointer to data.
			)=0;
 protected:
   const int Ndigits; ///< number of digits for counter in filenames.
   ///
   /// list of parameters (with names, types, data-ptrs)
   ///
   std::list<class pm_base *> params;
   ///
   /// parameters for Jet simulations
   ///
   std::list<class pm_base *> jet_pm;
   bool have_setup_jet_pm; ///< we only want to add to the list once!
   void set_jet_pm_params(); ///< add jet parameters to jet_pm list.
   ///
   /// parameters for any radiation sources.
   ///
   std::list<class pm_base *> rt_src;
   bool have_setup_rt_src; ///< we only want to add to the list once!
   void set_rt_src_params(); ///< add source parameters to rt_src list.
  
#ifndef OLD_TRACER

  ///
  /// Parameters for tracers.
  ///
  std::list<class pm_base *> tr_pm;
  bool have_setup_tracers; ///< track whether we have set up tracers
  void set_tracer_params(); ///< add tracer parameters to tr_pm list.

#endif // OLD_TRACER

   ///
   /// parameters for any stellar wind sources.
   ///
   std::list<class pm_base *> windsrc;
   bool have_setup_windsrc; ///< we only want to add to the list once!
   void set_windsrc_params(); ///< add source parameters to windsrc list.
   ///
   /// this makes a list of the params to read and should be called
   /// by the class constructor.
   ///
   void set_params();
   ///
   /// deletes the data associated with all param lists.
   ///
   void clear_param_list(std::list<class pm_base *> &);
   ///
   /// general driver function to read parameters which all derived 
   /// classes can use.
   ///
   int read_simulation_parameters();
   ///
   /// general driver function to write parameters which all derived 
   /// classes can use.
   ///
   int write_simulation_parameters();
   ///
   /// Function which defines how to get the data from the file, so 
   /// each derived class will define this differently.
   ///
   virtual int read_header_param(class pm_base *)=0;
   ///
   /// Function which defines how to write the data to the file, so 
   /// each derived class will define this differently.
   ///
   virtual int write_header_param(class pm_base *)=0;
   ///
   /// Check that header parameters are sensible
   ///
   int check_header_parameters();

};

///
/// The Text dataio class.  Note this is only intended for 1d and 2d test problems
/// for debugging purposes.  Outputs do not contain simulation
/// parameters and so cannot be used for restarts.
///
class dataio_text : public DataIOBase {
 public:
  dataio_text();
  ~dataio_text();
   /** \brief Set a pointer to the solver. */
   virtual void SetSolver(FV_solver_base * ///< Pointer to the solver (to get Eint,divB,Ptot)
			  );
   /** \brief Output simulation data to file. */
   virtual int OutputData(const string, ///< File to write to
			  class GridBaseClass *, ///< pointer to data.
			  const long int ///< number to stamp file with (e.g. timestep)
			  );
   /** \brief Reads simulation parameters from file. */
   virtual int ReadHeader(string ///< file to read from
			  );

  ///
  /// Write simulation header info to file
  ///
  virtual int WriteHeader(
          const string ///< file to write to (full, exact filename).
          ) {cout<<"can't write text header!\n"; return 99;}

   /** Having set up a grid with parameters from the header, this reads
    * data from the file and puts it on the grid.
    * */
   virtual int ReadData(string, ///< file to read from
			class GridBaseClass * ///< pointer to data.
			);
 protected:
   class GridBaseClass *gp; ///< pointer to computational grid.
   class FV_solver_base *eqn; ///< pointer to the solver, which knows the equations we are solving.
   class ReadParams *rp;     ///< pointer to read_parameters class instance.

  ///
  /// set filename based on counter, outfile-base-name.
  ///
  std::string set_filename(
      const std::string, /// outfile
      const long int     ///< counter
      );


   /** \brief Reads parameters from a text file (almost vestigial at this stage).
    * 
    * The parameterfile is from a command-line argument.  This function is only used for 
    * test problems which are setup algorithmically.  Usually the code starts with an 
    * initial condition file which contains all the required parameters in its header.
    * \retval 0 success
    * \retval 1 failure
    * */
   int get_parameters(string /**< Name of parameterfile.*/);


   /**  \brief Get initial conditions and populate grid with them.
    * 
    * Get type of IC from parameterfile.\n
    * If shocktube, then get which_riemann and assign left and right states.\n
    * Assign data to grid.\n
    * 1D only so far.
    * \retval 0 success
    * \retval 1 failure
    * */
   int assign_initial_data();


   /**  \brief Gets Initial left and right states for a Riemann Problem to solve (1D).
    * 
    * You pass this an integer, and pointers to left and right states, and it
    * gives you the appropriate state.
    * \retval 0 success
    * \retval 1 failure
    * */
   int get_riemann_ics(int, ///< int value of which_riemann, to tell it which problem to solve
		       double *, ///< pointer to left state.
		       double *, ///< pointer to right state.
		       double *  ///< pointer to left/right interface value (as a fraction of the range).
		       );


   /** \brief Add a low level of pseudo-random noise to the data.
    * 
    * This adds random noise to the data.  The integer flag determines what 
    * kind of noise to add:
    *  - n=1: add random noise to pressure, at fixed fraction of average pressure on grid.
    *  - n=2: add random noise to each cell, at fixed fraction of specific cell value.
    * Noise is adiabatic, so density increases with pressure to keep p/rho^gamma constant.
    * 
    * \retval 0 success
    * \retval 1 failure
    * */
   int add_noise2data(int,   ///< flag saying what kind of noise to add.
		      double ///< fractional level of noise you want, wrt mean value on grid.
		      );
   /** \brief Output Data to ASCII file.
    * 
    * This is a text file with positions and various quantities.
    * Restarting is not possible from these files.
    * \retval 0 success
    * \retval 1 failure
    * */
   int output_ascii_data(string ///< File name to write to.
			 );
   int read_header_param(class pm_base *)
   {rep.error("WHY USE read_header_param FOR TEXTIO?",0); return 0;}
   int write_header_param(class pm_base *)
   {rep.error("WHY USE write_header_param FOR TEXTIO?",0); return 0;}
};


#endif // DATAIO_H
