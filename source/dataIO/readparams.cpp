/** \file readparams.cc
 *
 * \brief Definition of member functions of read_parameters class.
 *
 * \author Jonathan Mackey
 *
 * See readparams.h for extensive documentation.
 * */
///
/// Modifications:\n
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
///

#include <fstream>

#include <sstream>
using namespace std;

#include "readparams.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */

// Reading parameters from parameterfile

ReadParams::ReadParams() {}

ReadParams::~ReadParams()
{
  //  cout <<"(ReadParams::~ReadParams) Deleting ReadParams class.\n";
  //  cout << " ... Done.\n";
}

int ReadParams::read_paramfile(const string &infile)
{
  ifstream paramfile;
  paramfile.open(infile.c_str());
  if (!paramfile.is_open()) {
    spdlog::error("Error opening file {} for writing.  Quitting...", infile);
    return (1);
  }
  string line;
  // initialise to empty.
  if (!params.empty()) params.clear();
  int iparam = 0;
  struct parameter p;
  while (!paramfile.eof()) {
    getline(paramfile, line);
    // If it's not a line containing a parameter, continue and read another
    // line.
    if (line.empty() == true || line.substr(0, 1) == "#") {
      continue;
    }
    // We have found a parameter, so read it into the array.
    istringstream fileline(line);
    p.name.clear();
    p.val.clear();
    fileline >> p.name >> p.val;
    params.push_back(p);
    // cout <<"param "<<iparam<<": name="<<params[iparam].name<<",
    // val="<<params[iparam].val<<"\n";
    iparam++;
  }
  paramfile.close();
  return 0;
}

void ReadParams::write_out_parameters()
{
  for (vector<struct parameter>::iterator i = params.begin(); i != params.end();
       ++i)
    spdlog::debug("{}\t{}", (*i).name, (*i).val);
  return;
}

// find and assign function
string ReadParams::find_parameter(const string &p)
{
  vector<struct parameter>::iterator i = params.begin();
  while (i != params.end() && (*i).name != p)
    ++i;
  if (i == params.end()) {
    spdlog::warn(
        "findparameter: couldn't find parameter: {} in file. Returning empty string.",
        p);
    return ("");
  }
  else
    return (*i).val;
}
