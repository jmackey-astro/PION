#!/bin/bash

mpirun -np 8 ../../icgen-ng params_DTE_HD_d3l1n0128.txt 
mpirun -np 8 ../../icgen-ng params_DTE_HD_d3l2n0128.txt 
mpirun -np 8 ../../icgen-ng params_DTE_HD_d3l3n0128.txt

mpirun -np 16 ../../pion-ng DTE_HD_d3l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d3l1n0128 opfreq=64   &
mpirun -np 16 ../../pion-ng DTE_HD_d3l2n0128_level00_0000.00000000.silo  \
  redirect=log_DTE_HD_d3l2n0128 opfreq=128  &
mpirun -np 16 ../../pion-ng DTE_HD_d3l3n0128_level00_0000.00000000.silo  \
  redirect=log_DTE_HD_d3l3n0128 opfreq=256  &
wait
exit


