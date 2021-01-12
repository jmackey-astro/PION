
This code projects 2D axisymmetric simulations of the interstellar medium onto
the plane of the sky.  It assumes that PION is using CGS units.

It solves the equation of radiative transfer for line emission, and sums up
material for column densities and emission measure.  It produces emission maps of total intensity in lines (CGS units of intensity: erg/cm2/s/sq.arcsec).

Compile with (replacing ubuntu18 with the appropriate flag for your operating system)
$ ABS="" NII="" MAKE_UNAME=ubuntu18 make -j6

This compiles with flags for not including absorption of line
emission from gas in the simulation, and without enhancing the
nitrogen abundance above the Solar value, and for an Ubuntu 18.04
operating system.
Within the Makefile are directives for "debian9", "ubuntu18",
"ubuntu16", "OSX" and "SUPERMUC".  If you have a different OS then
either modify the Makefile yourself or contact the authors.

X-ray table was generated using Xspec 12 APEC model with Asplund
abundances, and is in units erg.cm^{3}.s^{-1}, and should be divided
by 4pi to get an emissivity per solid angle.

run with, e.g.,
./projection2D /path/to/data/ FILENAME_0000.000 30 IMGFILE 3 0



