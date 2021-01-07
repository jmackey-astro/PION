#!/bin/bash
# Run on machine with ~64 cores.

time mpirun -np 8 ../../pion-ng DTEHD_d3l3n0064_level00_0000.00024968.silo \
  redirect=log_DTEHD_d3l3n0064_v2 &
time mpirun -np 32 ../../pion-ng DTEHD_d3l2n0128_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l2n0128 &
wait
exit

mpirun -np 1 ../../icgen-ng params_DTEHD_d3l1n0032.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d3l2n0032.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d3l3n0032.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d3l1n0064.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d3l2n0064.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d3l3n0064.txt 
mpirun -np 8 ../../icgen-ng params_DTEHD_d3l1n0128.txt 
mpirun -np 8 ../../icgen-ng params_DTEHD_d3l2n0128.txt 
mpirun -np 8 ../../icgen-ng params_DTEHD_d3l3n0128.txt 

time mpirun -np 1 ../../pion-ng DTEHD_d3l1n0032_0000.00000000.silo \
  redirect=log_DTEHD_d3l1n0032 &
time mpirun -np 1 ../../pion-ng DTEHD_d3l2n0032_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l2n0032 &
time mpirun -np 1 ../../pion-ng DTEHD_d3l3n0032_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l3n0032 &
#exit

time mpirun -np 8 ../../pion-ng DTEHD_d3l1n0064_0000.00000000.silo \
  redirect=log_DTEHD_d3l1n0064 &
time mpirun -np 8 ../../pion-ng DTEHD_d3l2n0064_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l2n0064 &
time mpirun -np 8 ../../pion-ng DTEHD_d3l3n0064_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l3n0064 &
wait
#exit

time mpirun -np 32 ../../pion-ng DTEHD_d3l1n0128_0000.00000000.silo \
  redirect=log_DTEHD_d3l1n0128
time mpirun -np 32 ../../pion-ng DTEHD_d3l2n0128_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l2n0128
time mpirun -np 32 ../../pion-ng DTEHD_d3l3n0128_level00_0000.00000000.silo \
  redirect=log_DTEHD_d3l3n0128
exit




