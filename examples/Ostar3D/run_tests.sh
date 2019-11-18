#!/bin/bash


mpirun -np 4 ../../icgen_NG_parallel params_Ostar3D_n0128l1.txt silo
mpirun -np 4 ../../pion_NG_parallel Ostar3D_n0128l1_0000.00000000.silo outfile=Ostar3D_n0128l1_V2
exit

mpirun -np 4 ../../icgen_NG_parallel params_Ostar3D_n0192l1.txt silo
mpirun -np 4 ../../pion_NG_parallel Ostar3D_n0192l1_0000.00000000.silo



