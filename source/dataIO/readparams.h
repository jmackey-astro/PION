/** \file readparams.h
 *
 * \brief Declaration of read_parameters class
 *
 * \author Jonathan Mackey
 *
 * This is a general class which will open a text file and parse it,
 * looking for parameter names and values, which it stores in an
 * array of strings (currently 50x2 elements).  The text file has to
 * be formatted in the proper way, with '#' for a comment, and parameters
 * given by "[parameter-name]  [parameter-value]".
 *
 * First read_paramfile(filename) must be called, and then string values
 * of parameters can be obtained by calling find_parameter(param-name).
 * */

#ifndef READPARAMS_H
#define READPARAMS_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <iostream>
#include <vector>
using namespace std;

/** \brief Class for parsing a text file to get parameter values.
 *
 *
 * */
struct parameter {
  std::string name;
  std::string val;
};

class ReadParams {
public:
  ReadParams();  /**< \brief Trivial Constructor */
  ~ReadParams(); /**< \brief Destructor; deletes parameter array. */
  /** \brief opens a file and parses it looking for parameters.
   *
   * This takes in a filename as an argument, and opens the file.
   * It then parses it, and stores any parameters it finds in
   * params[i][], with the param name in the first element and the
   * value in the second (as a string).  Closes the file before returning.
   *
   * \retval 0 success
   * \retval 1 failure
   * */
  int read_paramfile(const std::string& /**< Parameterfile to open. */);
  /** \brief prints all parameters and values to stdout.
   *
   * This function is just for testing really, as all it does
   * is print out all the parameters found in the file.
   * */
  void write_out_parameters();
  /** \brief Searches for a param-name, and returns its value
   *
   * This function is passed in a string, and searches the list
   * of parameters for a matching name, returning the value if
   * it finds a parameter, and a null string if it can't.
   * */
  std::string find_parameter(const std::string&);

private:
  std::vector<struct parameter> params;  ///< List of Parameters;
  /** \todo Make this a vector of strings! */
  // static const int arraylength=150;   /**< \brief size of array for
  // parameters (modify if needed).*/ std::string params[arraylength][2]; /**<
  // \brief array to store parameters in.*/
};

#endif
