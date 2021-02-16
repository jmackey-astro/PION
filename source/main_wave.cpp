/** \file wave.cc
 * \author Jonathan Mackey
 *
 * File to calculate postshock state given pre-shock and a mach no.
 * */
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

int main(int argc, char** argv)
{

  if (argc < 8) {
    cerr << "Please input preshock state and mach no.\n";
    cerr << "Usage: <wave> lro lpg lvx lvy lvz gamma machno\n";
    exit(1);
  }

  double lro   = atof(argv[1]);
  double lpg   = atof(argv[2]);
  double lvx   = atof(argv[3]);
  double lvy   = atof(argv[4]);
  double lvz   = atof(argv[5]);
  double gamma = atof(argv[6]);
  double mach  = atof(argv[7]);

  cout << "Calculating Postshock state for a mach " << mach << " shock\n";
  cout << "with left state [" << lro << ", " << lpg << ", " << lvx << ", "
       << lvy << ", " << lvz << "]\n";
  cout << "with ideal gas eos and gamma = " << gamma << "\n";

  cout.setf(ios_base::scientific);
  cout.precision(10);
  double rpg   = lpg * (2. * gamma * mach * mach - gamma + 1.) / (gamma + 1.);
  double alpha = (gamma + 1.) / (gamma - 1.);
  double rro   = lro * (1. + alpha * rpg / lpg) / (alpha + rpg / lpg);
  double rvx = lvx + (rpg / lpg - 1.) * sqrt(gamma * lpg / lro) / gamma / mach;
  double rvy = lvy;
  double rvz = lvz;

  cout << "Pre-shock state: density, pressure, x-vel\n";
  cout << "[" << lro << ", " << lpg << ", " << lvx << "]" << endl;
  cout << "Postshock state: density, pressure, x-vel\n";
  cout << "[" << rro << ", " << rpg << ", " << rvx << "]" << endl;
  return (0);
}
