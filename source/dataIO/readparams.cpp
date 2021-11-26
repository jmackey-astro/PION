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

#include <spdlog/spdlog.h>

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
  // for (int i=0; i<arraylength; i++) {
  //   params[i][0] = "";
  //   params[i][1] = "";
  //}
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

    // cout << line<<"\t"<<iparam<<"\n";
    // if (iparam>=arraylength)
    // spdlog::error("{}: {}", "too many parameters",iparam);
    // fileline >> params[iparam][0] >> params[iparam][1];
    // iparam++;
  }
  paramfile.close();
  return 0;
}

void ReadParams::write_out_parameters()
{
  for (vector<struct parameter>::iterator i = params.begin(); i != params.end();
       ++i)
    spdlog::debug("{}\t{}", (*i).name, (*i).val);

  // for (int iparam=0; iparam<arraylength; iparam++)
  //    cout << params[iparam][0] << "\t" << params[iparam][1] << "\n";
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

  // int i = 0;
  // while(i<arraylength && params[i][0] != p) {
  //    i++;
  // }
  // if (i >= arraylength) {
  //    cerr << "Error: findparameter: couldn't find parameter: "<<p<<" in
  //    file. Returning empty string." << "\n"; return("");
  // }
  // return(params[i][1]);
}
