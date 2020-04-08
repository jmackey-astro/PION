#!/bin/bash

../../icgen-ng params_ResStudy_HD_l3n0128.txt silo
../../icgen-ng params_ResStudy_HD_l3n0256.txt silo
../../icgen-ng params_ResStudy_HD_l3n0384.txt silo
../../icgen-ng params_ResStudy_HD_l3n0512.txt silo

echo "Running HD n128 simulation"
mpirun -np  4 ../../pion-ng ResStudy_HD_l3n0128_level00_0000.00000000.silo redirect=log_HD_l3n0128
echo "Running HD n256 simulation"
mpirun -np 16 ../../pion-ng ResStudy_HD_l3n0256_level00_0000.00000000.silo redirect=log_HD_l3n0256
echo "Running HD n384 simulation"
mpirun -np 64 ../../pion-ng ResStudy_HD_l3n0384_level00_0000.00000000.silo redirect=log_HD_l3n0384
echo "Running HD n512 simulation"
mpirun -np 64 ../../pion-ng ResStudy_HD_l3n0512_level00_0000.00000000.silo redirect=log_HD_l3n0512

../../icgen-ng params_ResStudy_MHD_l3n0128.txt silo
../../icgen-ng params_ResStudy_MHD_l3n0256.txt silo
../../icgen-ng params_ResStudy_MHD_l3n0384.txt silo
../../icgen-ng params_ResStudy_MHD_l3n0512.txt silo

echo "Running MHD n128 simulation"
mpirun -np  4 ../../pion-ng ResStudy_MHD_l3n0128_level00_0000.00000000.silo redirect=log_MHD_l3n0128
echo "Running MHD n256 simulation"
mpirun -np 16 ../../pion-ng ResStudy_MHD_l3n0256_level00_0000.00000000.silo redirect=log_MHD_l3n0256
echo "Running MHD n384 simulation"
mpirun -np 64 ../../pion-ng ResStudy_MHD_l3n0384_level00_0000.00000000.silo redirect=log_MHD_l3n0384
echo "Running MHD n512 simulation"
mpirun -np 64 ../../pion-ng ResStudy_MHD_l3n0512_level00_0000.00000000.silo redirect=log_MHD_l3n0512


