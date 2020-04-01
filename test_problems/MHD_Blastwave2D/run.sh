#!/bin/bash

../../icgen_NG_parallel params_MHD_blastwave2D_B001_n040.txt silo
../../icgen_NG_parallel params_MHD_blastwave2D_B001_n200.txt silo
../../pion_NG_parallel BW2d_StoneMHD_B001_n040_level00_0000.00000000.silo \
  outfile=HLLD_B001_n040_l2 solver=7
mpirun -np 4 ../../pion_NG_parallel \
  BW2d_StoneMHD_B001_n200_level00_0000.00000000.silo \
  outfile=HLLD_B001_n200_l2 solver=7

../../icgen_NG_parallel params_MHD_blastwave2D_B010_n040.txt silo
../../icgen_NG_parallel params_MHD_blastwave2D_B010_n200.txt silo
../../pion_NG_parallel BW2d_StoneMHD_B010_n040_level00_0000.00000000.silo \
  outfile=HLLD_B010_n040_l2 solver=7
mpirun -np 4 ../../pion_NG_parallel \
  BW2d_StoneMHD_B010_n200_level00_0000.00000000.silo \
  outfile=HLLD_B010_n200_l2 solver=7

../../icgen_NG_parallel params_MHD_blastwave2D_B100_n040.txt silo
../../icgen_NG_parallel params_MHD_blastwave2D_B100_n200.txt silo
../../pion_NG_parallel BW2d_StoneMHD_B100_n040_level00_0000.00000000.silo \
  outfile=HLLD_B100_n040_l2 solver=7
mpirun -np 4 ../../pion_NG_parallel \
  BW2d_StoneMHD_B100_n200_level00_0000.00000000.silo \
  outfile=HLLD_B100_n200_l2 solver=7


