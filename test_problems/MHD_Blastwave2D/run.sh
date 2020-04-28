#!/bin/bash

mpirun -np 1 ../../icgen-ng params_MHD_blastwave2D_l2_B010_n128.txt silo
mpirun -np 1 ../../icgen-ng params_MHD_blastwave2D_l2_B010_n256.txt silo
mpirun -np 1 ../../icgen-ng params_MHD_blastwave2D_l2_B010_n512.txt silo

mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B010_n128_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n128
mpirun -np 8 ../../pion-ng BW2d_StoneMHD_l2_B010_n256_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n256
mpirun -np 16 ../../pion-ng BW2d_StoneMHD_l2_B010_n512_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n512
exit

mpirun ../../icgen-ug params_MHD_blastwave2D_UG_B001_n256.txt silo
mpirun ../../icgen-ug params_MHD_blastwave2D_UG_B010_n256.txt silo
mpirun ../../icgen-ug params_MHD_blastwave2D_UG_B100_n256.txt silo

mpirun -np 4 ../../pion-ug BW2d_StoneMHD_UG_B100_n256_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=UG_B100_n256
mpirun -np 4 ../../pion-ug BW2d_StoneMHD_UG_B010_n256_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=UG_B010_n256
mpirun -np 4 ../../pion-ug BW2d_StoneMHD_UG_B001_n256_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=UG_B001_n256


mpirun ../../icgen-ng params_MHD_blastwave2D_l2_B000_n256.txt silo
mpirun ../../icgen-ng params_MHD_blastwave2D_l2_B001_n256.txt silo
mpirun ../../icgen-ng params_MHD_blastwave2D_l2_B010_n256.txt silo
mpirun ../../icgen-ng params_MHD_blastwave2D_l2_B100_n256.txt silo

mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B000_n256_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B000_n256
mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B001_n256_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B001_n256
mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B010_n256_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n256
mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B010_n256_level00_0000.00000000.silo solver=8 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n256_HLL
mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B010_n256_level00_0000.00000000.silo solver=7 ooa=1 cfl=0.24 opfreq_time=0.1 outfile=NG_B010_n256_HLLDOA1
mpirun -np 4 ../../pion-ng BW2d_StoneMHD_l2_B100_n256_level00_0000.00000000.silo solver=7 cfl=0.24 opfreq_time=0.1 outfile=NG_B100_n256

exit


