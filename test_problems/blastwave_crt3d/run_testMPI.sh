#!/bin/bash

../../icgen_NG_parallel params_BW_NGcrt3D_NR032.txt silo

mpirun -np  1 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np01 opfreq=64 redirect=log_np01
mpirun -np  2 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np02 opfreq=64 redirect=log_np02
mpirun -np  4 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np04 opfreq=64 redirect=log_np04
mpirun -np  8 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np08 opfreq=64 redirect=log_np08
mpirun -np 16 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np16 opfreq=64 redirect=log_np16
mpirun -np 32 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np32 opfreq=64 redirect=log_np32
mpirun -np 64 ../../pion_NG_parallel BW_NGcrt3D_NR032_level00_0000.00000000.silo outfile=NG3D_NR032_np64 opfreq=64 redirect=log_np64

grep TOTALS log_np*.txt
exit

