#!/bin/bash

# generate the initial conditions from the parameter file
# This setup has 3 levels of refinement, and a magnetised wind from
# an O star that is rotating at 100 km/s.  It is very zoomed out from
# the star, so the wind has an almost purely toroidal component in
# the boundary region.
# The star is moving with 25 km/s relative to the interstellar medium
# and it is assumed that the star photoionizes its surroundings.
mpirun -np 16 ../../icgen-ng params_Ostar3D_B010_n0128l3.txt silo
mpirun -np 16 ../../icgen-ng params_Ostar3D_B100_n0128l3.txt silo

# run PION with 8 MPI processes using 8 cores, loading the Ostar3D
# initial conditions snapshot
mpirun -np 32 ../../pion-ng Ostar3D_B010_n0128l3_level00_0000.00000000.silo opfreq=512 \
redirect=/mnt/jmackey/jm/data/pion_codepaper/Ostar3D/log_Ostar3D_B010_n0128l3 \
outfile=/mnt/jmackey/jm/data/pion_codepaper/Ostar3D/Ostar3D_B010_n0128l3_node &

mpirun -np 32 ../../pion-ng Ostar3D_B100_n0128l3_level00_0000.00000000.silo opfreq=512 \
redirect=/mnt/jmackey/jm/data/pion_codepaper/Ostar3D/log_Ostar3D_B100_n0128l3 \
outfile=/mnt/jmackey/jm/data/pion_codepaper/Ostar3D/Ostar3D_B100_n0128l3_node &

wait




