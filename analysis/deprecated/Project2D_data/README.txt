
This code projects 2D axisymmetric simulations of the interstellar medium onto
the plane of the sky.  It assumes that PION is using CGS units, and that the
first tracer variable is the ionization fraction of Hydrogen.

It may have a few hacks here and there, so look in the source code to convince
yourself that it does what it says it does.

It solves the equation of radiative transfer for line emission, and sums up
material for column densities and emission measure.  It produces emission maps of total intensity in lines (CGS units of intensity: erg/cm2/s/sq.arcsec).


