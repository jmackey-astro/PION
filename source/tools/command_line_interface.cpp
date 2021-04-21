///
/// \file command_line_interface.cpp
/// \author Andrew J. Lim, Jonathan Mackey
///
/// This file defines the command line interface, for debugging, and
/// the debugging class containing bookkeeping values and pointers.
///
/// Modifications:
/// - 2015.01.12 JM: created file, moved class from global.h

#include "constants.h"
#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "grid/cell_interface.h"
#include "tools/command_line_interface.h"
#include "tools/reporting.h"
using namespace std;

#ifdef TESTING
class DebugParams dp;
class CommandLineInterface commandline;

DebugParams::DebugParams()
{
  //   cout <<"initialising DebugParams...\n";
  initERG = initMMX = initMMY = initMMZ = 0.0;
  ergTotChange = mmxTotChange = mmyTotChange = mmzTotChange = 0.0;
  c                                                         = 0;
}

DebugParams::~DebugParams() {}

/************************ COMMAND LINE INTERFACE
 * ************************************/
CommandLineInterface::CommandLineInterface()
{
  //  cout <<"Setting up interactive debugging.\n";
}
CommandLineInterface::~CommandLineInterface()
{
  //  cout <<"Destructing interactive debugging.\n";
}

void CommandLineInterface::auto_console(char *prompt)
{
  cout << "\n Welcome to the AUTOPILOT command line\n";
  cout << "****************************************\n";
  cout << prompt << "\n";
  cout << "going to first point\n";
  fpt();
  print_cell();
  end_of_col("YP");
  print_cell();
  next_point("YP");
  print_cell();
  next_point("XN");
  print_cell();
  next_point("YN");
  print_cell();
  next_point("XP");
  print_cell();

  cout << "\n LEAVING the AUTOPILOT command line\n";
  cout << "****************************************\n";
  return;
}

void CommandLineInterface::console(std::string p)
{
  char prompt[512];
  strcpy(prompt, p.c_str());
  //
  // A command line interface which allows investigation of
  //  the grid during runtime.
  //
  fprintf(stderr, "\n Welcome to the command line - type \"help\"\n\n");
  read_history(".cmddemo_CLI.hist");
  //
  // Get a command.
  //
  char *comm = NULL;
  for (;;) {
    comm    = rl_gets(comm, prompt);
    char *c = comm;
    do {
      char comm1[512], *d = comm1;
      *d = 0;
      while (*c && ((*(c)) != ';')) {
        *d++     = *c++;
        *(d + 1) = 0;
      }  // <--- (A)
      if ((*c) == ';') c++;
      if (execute(comm1)) return;
    } while (*c);
  }
}

char *CommandLineInterface::rl_gets(char *line, char *prompt)
{
  //
  // Reads a string, and return a pointer to it.  Returns NULL on EOF.
  //
  if (line) {
    //
    // If the buffer has already been allocated, return the memory
    //      to the free pool.
    //
    free(line);
    line = (char *)NULL;
  }
  //
  // Get a line from the user.
  //
  line = readline(prompt);
  //
  // If the line has any text in it, save it on the history.
  //
  if (line && *line && line[0] != 'q') add_history(line);
  return (line);
}

int CommandLineInterface::execute(char *com)
{
  //
  // Executes the command in com[], Could be smarter than to hardcode the
  //  command lengths.
  //
  while ((*(com)) == ' ') {
    (*com)++;
  }  // strip off any leading white space
  if (!(*com))
    return 0;  // don't execute NULL commands!!!
  else if (!strncmp(com, "!", 1))
    system(com + 1);
  else if (!strncmp(com, "cmd1", 4))
    cmd1();
  else if (!strncmp(com, "cmd2", 4))
    cmd2(com + 4);
  else if (!strncmp(com, "bigcmd", 6))
    bigcmd();
  else if (!strncmp(com, "print_cell", 10))
    print_cell();
  else if (!strncmp(com, "print_flux", 10)) {
    string s = com;
    s        = s.substr(11, 2);
    cout << s << "\n";
    print_flux(atoi(s.c_str()));
  }
  else if (!strncmp(com, "next_point", 10)) {
    string s = com;
    s        = s.substr(11, 2);
    cout << s << "\n";
    next_point(s);
  }
  else if (!strncmp(com, "fpt", 3))
    fpt();
  else if (!strncmp(com, "lpt", 3))
    lpt();
  else if (!strncmp(com, "end_of_col", 10)) {
    string s = com;
    s        = s.substr(11, 2);
    cout << s << "\n";
    end_of_col(s);
  }
  else if (!strncmp(com, "help", 4)) {
    fprintf(stderr, "\n commands:\n");
    fprintf(stderr, "           cmd1 - demo cmd\n");
    fprintf(stderr, "           cmd2 - demo cmd with args\n");
    fprintf(stderr, "         bigcmd - no particular reason for this\n");
    cerr << "print_cell - prints info about the current cell.\n";
    cerr << "print_flux - If in dynamics solver, prints intercell fluxes.\n";
    cerr << "next_point - move from current cell to next one, in a direction "
            "(XY,XN,etc.)\n";
    cerr << "fpt() - moves dp.c to the first point.\n";
    cerr << "lpt() - moves dp.c to the last point.\n";
    cerr << "end_of_col(string) - moves dp.c to the edge of the grid in "
            "direction (XN,XP,ZN,ZP,YN,YP).\n";
    fprintf(stderr, "           help - prints this list\n");
    fprintf(stderr, "  q, quit, exit - end the program\n");
    fprintf(stderr, "\n leading \"!\" sends cmd to shell\n\n");
  }
  else if (!strncmp(com, "cont", 4)) {
    fprintf(stderr, " Continuing with code execution.\n");
    write_history(".cmddemo_CLI.hist");
    history_truncate_file(".cmddemo_CLI.hist", 50);
    return 1;
  }
  else if (
      (!strncmp(com, "q", 1)) || (!strncmp(com, "quit", 4))
      || (!strncmp(com, "exit", 4))) {
    write_history(".cmddemo_CLI.hist");
    history_truncate_file(".cmddemo_CLI.hist", 50);
    rep.error("Quitting at user request", com);
  }
  //   else if  (!strncmp(com,"shutdown",8)) system("halt -p"); Steph wouldn't
  //   like this...
  else {
    fprintf(stderr, " Unknown command: %s\n\n", com);
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                        Commands are below.                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////
void CommandLineInterface::cmd1()
{
  fprintf(stderr, " This is command 1\n");
  double x[MAX_DIM];
  CI.get_dpos(dp.c, x);
  rep.printVec("position", x, 2);
}

void CommandLineInterface::cmd2(char *args)
{
  fprintf(stderr, " This is command 2, args: \"%s\"\n", args);
  rep.printVec("P ", dp.c->P, 5);
  rep.printVec("Ph", dp.c->P, 5);
}

void CommandLineInterface::bigcmd()
{
  fprintf(stderr, " This is a big command (apparently).\n");
}

void CommandLineInterface::print_cell()
{
  CI.print_cell(dp.c);
  cout << "\n";
  return;
}

void CommandLineInterface::next_point(string s)
{
  enum direction dir = parse_dir(s);
  if (dir == NO) return;

  if (!(dp.c->ngb[dir])) {
    cout << "no neighbour in direction " << s << "; try again.\n";
    return;
  }
  else {
    cout << "moving from cell " << dp.c->id;
    dp.c = dp.c->ngb[dir];
    cout << " to cell " << dp.c->id << "\n";
  }
  return;
}

void CommandLineInterface::print_flux(
    const int nvar  ///< Length of state vectors
)
{
  rep.printVec("left ", dp.vec2, nvar);
  rep.printVec("right", dp.vec3, nvar);
  rep.printVec("pstar", dp.vec1, nvar);
  rep.printVec("flux ", dp.vec4, nvar);
  return;
}

void CommandLineInterface::fpt()
{
  dp.c = dp.grid->FirstPt();
  return;
}

void CommandLineInterface::lpt()
{
  dp.c = dp.grid->LastPt();
  return;
}

void CommandLineInterface::end_of_col(const string s)
{
  if (!(dp.c)) {
    cout << "dp.c is null pointer! can't get to end of column.\n";
    return;
  }
  enum direction dir = parse_dir(s);
  if (dir == NO) return;

  int ct = 0;
  while (dp.grid->NextPt(dp.c, dir)->isgd) {
    dp.c = dp.grid->NextPt(dp.c, dir);
    ct++;
  }
  cout << "moved " << ct << " cells in direction " << s << "\n";
  return;
}

enum direction CommandLineInterface::parse_dir(const string s)
{
  enum direction dir = NO;
  if (s == "XP")
    dir = XP;
  else if (s == "XN")
    dir = XN;
  else if (s == "YP")
    dir = YP;
  else if (s == "YN")
    dir = YN;
  else if (s == "ZP")
    dir = ZP;
  else if (s == "ZN")
    dir = ZN;
  else {
    cout << "bad direction string: " << s << ", try again.\n";
  }
  return dir;
}

#endif  // TESTING
