#!/bin/bash

mpirun -np  1 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np01_n128 redirect=log_np01_n128 &
mpirun -np  2 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np02_n128 redirect=log_np02_n128 &
mpirun -np  4 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np04_n128 redirect=log_np04_n128 &
mpirun -np  8 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np08_n128 redirect=log_np08_n128 &
mpirun -np 16 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np16_n128 redirect=log_np16_n128 &
mpirun -np 32 ../../pion-ng fastwind3d_n0128_l2_level00_0000.00000000.silo outfile=np32_n128 redirect=log_np32_n128 &

wait
exit

