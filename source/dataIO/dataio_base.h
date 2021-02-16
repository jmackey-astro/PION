/// \file dataio_base.h
/// \author Jonathan Mackey
///
/// modified:\n
/// - 2018.05.01 JM: moved DataIOBase from dataio.h

#ifndef DATAIO_BASE_H
#define DATAIO_BASE_H

#include "dataIO/file_status.h"
#include "parameter_defs.h"
#include "sim_params.h"
#include "spatial_solvers/solver_eqn_base.h"

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
  DataIOBase(class SimParams&  ///< pointer to simulation parameters
  );

  ~DataIOBase();

  /// Set a pointer to the solver.
  virtual void SetSolver(
      FV_solver_base*  ///< Pointer to solver (for Eint,divB,Ptot)
      ) = 0;

  ///
  /// Write simulation header info to file
  ///
  virtual int WriteHeader(
      const string,     ///< file to write to (full, exact filename).
      class SimParams&  ///< pointer to simulation parameters
      ) = 0;

  ///
  /// Write simulation data to file.
  ///
  virtual int OutputData(
      const string,                   ///< File to write to
      vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
      class SimParams&,               ///< pointer to simulation parameters
      const long int  ///< number to stamp file with (e.g. timestep)
      ) = 0;

  ///
  /// Reads simulation parameters from file.
  ///
  virtual int ReadHeader(
      string,           ///< file to read from
      class SimParams&  ///< pointer to simulation parameters
      ) = 0;

  ///
  /// Having set up a grid with parameters from the header, this reads
  /// data from the file and puts it on the grid.
  ///
  virtual int ReadData(
      string,                         ///< file to read from
      vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
      class SimParams&                ///< pointer to simulation parameters
      ) = 0;

protected:
  const int Ndigits;  ///< number of digits for counter in filenames.
  ///
  /// list of parameters (with names, types, data-ptrs)
  ///
  std::list<class pm_base*> params;
  ///
  /// parameters for Boundary Condition
  ///
  std::list<class pm_base*> bc_pm;
  bool have_setup_bc_pm;  ///< we only want to add to the list once!
  void set_bc_pm_params(class SimParams&  ///< pointer to simulation parameters
  );                                      ///< add BC parameters to bc_pm list.

  ///
  /// parameters for Jet simulations
  ///
  std::list<class pm_base*> jet_pm;
  bool have_setup_jet_pm;    ///< we only want to add to the list once!
  void set_jet_pm_params();  ///< add jet parameters to jet_pm list.

  ///
  /// parameters for any radiation sources.
  ///
  std::list<class pm_base*> rt_src;
  bool have_setup_rt_src;  ///< we only want to add to the list once!
  void set_rt_src_params(class SimParams&  ///< pointer to simulation parameters
  );  ///< add source parameters to rt_src list.

#ifndef OLD_TRACER

  ///
  /// Parameters for tracers.
  ///
  std::list<class pm_base*> tr_pm;
  bool have_setup_tracers;  ///< track whether we have set up tracers
  void set_tracer_params(class SimParams&  ///< pointer to simulation parameters
  );  ///< add tracer parameters to tr_pm list.

#endif  // OLD_TRACER

  ///
  /// parameters for any stellar wind sources.
  ///
  std::list<class pm_base*> windsrc;
  bool have_setup_windsrc;    ///< we only want to add to the list once!
  void set_windsrc_params();  ///< add source parameters to windsrc list.
  ///
  /// this makes a list of the params to read and should be called
  /// by the class constructor.
  ///
  void set_params(class SimParams&  ///< pointer to simulation parameters
  );
  ///
  /// deletes the data associated with all param lists.
  ///
  void clear_param_list(std::list<class pm_base*>&);
  ///
  /// general driver function to read parameters which all derived
  /// classes can use.
  ///
  int read_simulation_parameters(
      class SimParams&  ///< pointer to simulation parameters
  );

  ///
  /// general driver function to write parameters which all derived
  /// classes can use.
  ///
  int write_simulation_parameters(
      class SimParams&  ///< pointer to simulation parameters
  );

  ///
  /// Function which defines how to get the data from the file, so
  /// each derived class will define this differently.
  ///
  virtual int read_header_param(class pm_base*) = 0;
  ///
  /// Function which defines how to write the data to the file, so
  /// each derived class will define this differently.
  ///
  virtual int write_header_param(class pm_base*) = 0;

  ///
  /// Check that header parameters are sensible
  ///
  int check_header_parameters(
      class SimParams&  ///< pointer to simulation parameters
  );
};

#endif  // DATAIO_BASE_H
