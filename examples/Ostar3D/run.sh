#!/bin/bash

# generate the initial conditions from the parameter file
# This setup has 3 levels of refinement, and a magnetised wind from
# an O star that is rotating at 100 km/s.  It is very zoomed out from
# the star, so the wind has an almost purely toroidal component in
# the boundary region.
# The star is moving with 20 km/s relative to the interstellar medium
# and it is assumed that the star photoionizes its surroundings.
../../icgen_NG_parallel params_Ostar3D_n0128l3.txt silo

# run PION with 8 MPI processes using 8 cores, loading the Ostar3D
# initial conditions snapshot, and saving snapshots every 32
# timesteps into the current working directory.
mpirun -np 8 ../../pion_NG_parallel Ostar3D_n0128l3_level00_0000.00000000.silo opfreq=32


