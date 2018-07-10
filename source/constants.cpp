///
/// \file constants.cpp
/// \author Jonathan Mackey
/// \date 2015.03.10
///
/// Class with inline functions for many physics and astronomy
/// constants.  All constants are in cgs units unless stated
/// otherwise.
///
/// Modifications:
/// - 2015.03.10 JM: added equalD() function.


#include <iostream>
#include "constants.h"
using namespace std;


///
/// pion constants class.
///
class constants pconst;


// ##################################################################
// ##################################################################



constants::constants() {}


// ##################################################################
// ##################################################################



constants::~constants() {}


// ##################################################################
// ##################################################################




bool constants::equalD(
        const double a,
        const double b
        )
{
  if (a==b) {
    return(true);
  }
  if (fabs(a)+fabs(b)<TINYVALUE) {
#ifdef TESTING
    cout <<"tiny numbers in equalD(a,b); a,b <1.e-100... a=";
    cout <<a<<", b="<<b<<"; returning true.\n";
#endif
    return(true);
  }
  if ( (fabs(a-b)/(fabs(a)+fabs(b)+TINYVALUE)) <SMALLVALUE) {
    return(true);
  }
  else {
    return(false); // false is zero.
  }
}


// ##################################################################
// ##################################################################



