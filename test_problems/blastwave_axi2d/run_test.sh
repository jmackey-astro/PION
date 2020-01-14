#!/bin/bash


mpirun -np 32 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np32 opfreq=1024 redirect=log_np32 &
mpirun -np 16 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np16 opfreq=1024 redirect=log_np16 &
mpirun -np  8 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np08 opfreq=1024 redirect=log_np08 &
mpirun -np  4 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04 opfreq=1024 redirect=log_np04 &
mpirun -np  2 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np02 opfreq=1024 redirect=log_np02 &
mpirun -np  1 ../../pion_NG_parallel BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np01 opfreq=1024 redirect=log_np01 &
wait

grep TOTALS log_np*.txt
exit

