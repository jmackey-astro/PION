///
/// \file command_line_interface.h
/// \author Andrew J. Lim, Jonathan Mackey
///
/// This file declares the command line interface, for debugging, and
/// the debugging class containing bookkeeping values and pointers.
///
/// Modifications:
/// - 2015.01.12 JM: created file, moved class from global.h

#ifndef COMMAND_LINE_INTERFACE_H
#define COMMAND_LINE_INTERFACE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "grid/grid_base_class.h"

#include <string>

#ifdef TESTING

#include <readline/history.h>
#include <readline/readline.h>

/// \brief class for debugging the code, tracking energy/momentum to make
/// sure it is conserved, etc.
///
/// Also contains a cell pointer, so we know what cell we are working on at
/// all times.
/// The grid and MP pointers are globally accessible pointers to the
/// grid data and microphysics classes, so I don't need them here, but I
/// do need pointers to the solver and the integrator, so they are
/// included, and are set whenever the respective classes are initialised.
///
class DebugParams {
public:
  DebugParams();
  ~DebugParams();
  class cell *c;  ///< Pointer that can be used to track which cell we are at.
  class GridBaseClass *grid;  ///< pointer to grid.
  // class BaseFVSolver *solver; ///< pointer to solver.
  double initERG;  ///< Initial total energy on the grid (not including ghost
                   ///< cells).
  double initMMX;  ///< Initial X-momentum on the grid (not including ghost
                   ///< cells).
  double initMMY;  ///< Initial Y-momentum on the grid (not including ghost
                   ///< cells).
  double initMMZ;  ///< Initial Z-momentum on the grid (not including ghost
                   ///< cells).
  double ergTotChange;  ///< For tracking energy entering/leaving domain.
  double mmxTotChange;  ///< For tracking x-momentum entering/leaving domain.
  double mmyTotChange;  ///< For tracking y-momentum entering/leaving domain.
  double mmzTotChange;  ///< For tracking z-momentum entering/leaving domain.
  double vec1[MAX_NVAR];
  double vec2[MAX_NVAR];
  double vec3[MAX_NVAR];
  double vec4[MAX_NVAR];
};

extern class DebugParams dp;

/// \brief Command Line Interface for interactive debugging of code.
///
/// Adapted from base code Andy gave me.
///
class CommandLineInterface {
public:
  CommandLineInterface();   ///< Do-nothing Constructor.
  ~CommandLineInterface();  ///< Do-nothing Destructor.
  /// This function calls up the CLI for interactive debugging.
  void console(std::string  ///< Optional text for prompt.
  );

  void auto_console(char *  ///< Optional text for prompt.
  );

private:
  /// Reads a string, and return a pointer to it.  Returns NULL on EOF.
  char *rl_gets(
      char *,  ///< pointer to line of text.
      char *   ///< text inputted to prompt.
  );
  /// Executes any commands the CLI recognises in the string passed in.
  int execute(char *);
  void cmd1();        ///< An example command, with no arguments.
  void cmd2(char *);  ///< another example command, with arguments.
  void bigcmd();      ///< A further example command.

  /// prints cell which dp.c points to.(DebugParams)
  void print_cell();

  /// moves to cell in direction given (XN,XP,YN,YP,ZN,ZP)
  void next_point(const std::string);

  /// moves dp.c to the first point.
  void fpt();

  /// moves dp.c to the last point.
  void lpt();

  /// moves dp.c to the edge of the grid in direction specified.
  void end_of_col(const std::string  ///< direction string
  );

  /// prints flux data from intercell flux function.
  void print_flux(const int  ///< length of state vectors
  );

  /// given a string of XN,YN,etc, return direction.
  enum direction parse_dir(const std::string);
};

extern class CommandLineInterface commandline;

#endif  // TESTING

#endif  // COMMAND_LINE_INTERFACE_H
