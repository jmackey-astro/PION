#!/bin/bash

mpirun -np 1 ../../icgen-ng params_DTEHD_d2l1n0064.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l2n0064.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l3n0064.txt 
time mpirun -np 4 ../../pion-ng DTEHD_d2l1n0064_0000.00000000.silo \
  redirect=log_DTEHD_d2l1n0064
time mpirun -np 4 ../../pion-ng DTEHD_d2l2n0064_0000.00000000.silo \
  redirect=log_DTEHD_d2l2n0064
time mpirun -np 4 ../../pion-ng DTEHD_d2l3n0064_0000.00000000.silo \
  redirect=log_DTEHD_d2l3n0064
#exit

mpirun -np 1 ../../icgen-ng params_DTEHD_d2l1n0128.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l2n0128.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l2n0128.txt 
time mpirun -np 16 ../../pion-ng DTEHD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTEHD_d2l1n0128
time mpirun -np 16 ../../pion-ng DTEHD_d2l2n0128_0000.00000000.silo \
  redirect=log_DTEHD_d2l2n0128
time mpirun -np 16 ../../pion-ng DTEHD_d2l3n0128_0000.00000000.silo \
  redirect=log_DTEHD_d2l3n0128
exit

mpirun -np 1 ../../icgen-ng params_DTEHD_d2l1n0032.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l2n0032.txt 
mpirun -np 1 ../../icgen-ng params_DTEHD_d2l3n0032.txt 
time mpirun -np 1 ../../pion-ng DTEHD_d2l1n0032_0000.00000000.silo \
  redirect=log_DTEHD_d2l1n0032
time mpirun -np 1 ../../pion-ng DTEHD_d2l2n0032_level00_0000.00000000.silo \
  redirect=log_DTEHD_d2l2n0032
time mpirun -np 1 ../../pion-ng DTEHD_d2l3n0032_level00_0000.00000000.silo \
  redirect=log_DTEHD_d2l3n0032
exit


mpirun -np 1 ../../icgen-ng params_DTEHD_d2l1n0256.txt 
time mpirun -np 4 ../../pion-ng DTEHD_d2l1n0256_0000.00000000.silo \
  redirect=log_DTEHD_d2l1n0256
exit


