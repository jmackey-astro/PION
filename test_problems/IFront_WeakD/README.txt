

* params_WeakD_nXXXX.txt
These are 1D Cartesian simulations of plane-parallel radiation hitting the left side of the domain, and an advecting ISM from right to left, leaving what should be a static D-type I-front on the domain.

* params_WeakD_Cartesian_nXXXX.txt
These are the same, but with a point source that has a 1/r^2 flux law.

* params_WeakD_Spherical_n0040.txt
This is in spherical coordinates, 1D, with a point source and an advecting domain.  It is not a physical setup, but just to see if the simulation has the velocity spike.  It does.

So I have found that if the I-front is static, then the velocity spike disappears.  It is only if I have a travelling I-front that it shows up.  And it must be to do with the mixed-cell physics.  It occurs in the densest cell that still has some ionisation, so deep in the I-front, where the density gradient is large, and the ionisation length is hard to resolve.

I can set up a travelling I-front, but it is difficult to resolve the ionisation length with the current setup.  I think I should increase the neutral gas temperature, so the sound speeds are not so different.  I could set T_n=1000K.  That would make the density discontinuity weaker, and so easier to resolve the ionisation length.


