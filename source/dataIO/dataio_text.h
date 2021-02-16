/// \file dataio_text.h
/// \author Jonathan Mackey
///
/// This file contains class declarations for ascii text data I/O.
///
/// modified:\n
/// - 2018.05.01 JM: moved dataio_text from dataio.cpp.

#ifndef DATAIO_TEXT_H
#define DATAIO_TEXT_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "dataio_base.h"

///
/// The Text dataio class.  Note this is only intended for 1d and 2d test
/// problems for debugging purposes.  Outputs do not contain simulation
/// parameters and so cannot be used for restarts.
///
class dataio_text : public DataIOBase {
public:
  dataio_text(class SimParams&  ///< pointer to simulation parameters
  );

  ~dataio_text();

  /// Set a pointer to the solver.
  virtual void SetSolver(
      FV_solver_base*  ///< Pointer to the solver (to get Eint,divB,Ptot)
  );

  /// Write simulation data to file.
  virtual int OutputData(
      const string,                   ///< File to write to
      vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
      class SimParams&,               ///< pointer to simulation parameters
      const long int  ///< number to stamp file with (e.g. timestep)
  );

  ///
  /// Reads simulation parameters from file.
  ///
  virtual int ReadHeader(
      string,                 ///< file to read from
      class SimParams& SimPM  ///< pointer to simulation parameters
  );

  ///
  /// Write simulation header info to file
  ///
  virtual int WriteHeader(
      const string,     ///< file to write to (full, exact filename).
      class SimParams&  ///< pointer to simulation parameters
  )
  {
    cout << "can't write text header!\n";
    return 99;
  }

  ///
  /// Read Simulation data form a file.
  ///
  virtual int ReadData(
      string,                         ///< file to read from
      vector<class GridBaseClass*>&,  ///< address of vector of grid pointers.
      class SimParams&                ///< pointer to simulation parameters
  );

protected:
  class GridBaseClass* gp;    ///< pointer to computational grid.
  class FV_solver_base* eqn;  ///< pointer to the solver, which knows the
                              ///< equations we are solving.
  class ReadParams* rp;       ///< pointer to read_parameters class instance.

  ///
  /// set filename based on counter, outfile-base-name.
  ///
  std::string set_filename(
      const std::string,  /// outfile
      const long int      ///< counter
  );

  /// Reads parameters from a text file (almost vestigial at this stage).
  ///
  /// The parameterfile is from a command-line argument.  This function is
  /// only used for test problems which are setup algorithmically.  Usually
  /// the code starts with an initial condition file which contains all the
  /// required parameters in its header. \retval 0 success \retval 1 failure
  ///
  int get_parameters(
      string,           ///< Name of parameterfile.
      class SimParams&  ///< pointer to simulation parameters
  );

  /// Get initial conditions and populate grid with them.
  ///
  /// Get type of IC from parameterfile.\n
  /// If shocktube, then get which_riemann and assign left and right states.\n
  /// Assign data to grid.\n
  /// 1D only so far.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int assign_initial_data(
      class SimParams&  ///< pointer to simulation parameters
  );

  ///  \brief Gets Initial left and right states for a Riemann Problem to
  ///  solve (1D).
  ///
  /// You pass this an integer, and pointers to left and right states, and it
  /// gives you the appropriate state.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int get_riemann_ics(
      class SimParams&,  ///< pointer to simulation parameters
      int,  ///< int value of which_riemann, to tell it which problem to solve
      double*,  ///< pointer to left state.
      double*,  ///< pointer to right state.
      double*   ///< pointer to left/right interface value (as a fraction of
                ///< the range).
  );

  /// Add a low level of pseudo-random noise to the data.
  ///
  /// This adds random noise to the data.  The integer flag determines what
  /// kind of noise to add:
  ///  - n=1: add random noise to pressure, at fixed fraction
  ///         of average pressure on grid.
  ///  - n=2: add random noise to each cell, at fixed fraction
  ///         of specific cell value.  Noise is adiabatic, so
  ///         density increases with pressure to keep p/rho^gamma
  ///         constant.
  ///
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int add_noise2data(
      class SimParams&,  ///< pointer to simulation parameters
      int,               ///< flag saying what kind of noise to add.
      double  ///< fractional level of noise you want, wrt mean value on grid.
  );

  /// Save Data to ASCII file.
  ///
  /// This is a text file with positions and various quantities.
  /// Restarting is not possible from these files.
  /// \retval 0 success
  /// \retval 1 failure
  ///
  int output_ascii_data(
      string,           ///< File name to write to.
      class SimParams&  ///< pointer to simulation parameters
  );

  int read_header_param(class pm_base*)
  {
    rep.error("WHY USE read_header_param FOR TEXTIO?", 0);
    return 0;
  }
  int write_header_param(class pm_base*)
  {
    rep.error("WHY USE write_header_param FOR TEXTIO?", 0);
    return 0;
  }
};

#endif  // DATAIO_TEXT_H
