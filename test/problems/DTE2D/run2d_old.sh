#!/bin/bash

mpirun -np 1 ../../icgen-ng params_DTE_HD_d2l1n0128.txt 
mpirun -np 1 ../../icgen-ng params_DTE_HD_d2l2n0128.txt 
mpirun -np 1 ../../icgen-ng params_DTE_HD_d2l3n0128.txt
mpirun -np 1 ../../icgen-ng params_DTE_HD_d2l4n0128.txt
mpirun -np 1 ../../icgen-ng params_DTE_HD_d2l5n0128.txt

mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=8 finishtime=1.0e12
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=16 finishtime=3.0e12
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=32 finishtime=1.0e13
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=64 finishtime=3.0e13
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=128 finishtime=1.0e14
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=256 finishtime=3.0e14
mpirun -np 8 ../../pion-ng DTE_HD_d2l1n0128_0000.00000000.silo \
  redirect=log_DTE_HD_d2l1n0128 opfreq=512 finishtime=1.1e15

mpirun -np 8 ../../pion-ng DTE_HD_d2l2n0128_level00_0000.00000000.silo \
  redirect=log_DTE_HD_d2l2n0128 opfreq=16 finishtime=1.0e12
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=32 finishtime=3.0e12
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=64 finishtime=1.0e13
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=128 finishtime=3.0e13
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=256 finishtime=1.0e14
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=512 finishtime=3.0e14
REST=`ls DTE_HD_d2l2n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l2n0128 opfreq=1024 finishtime=1.1e15

mpirun -np 8 ../../pion-ng DTE_HD_d2l3n0128_level00_0000.00000000.silo \
  redirect=log_DTE_HD_d2l3n0128 opfreq=32 finishtime=1.0e12
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=64 finishtime=3.0e12
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=128 finishtime=1.0e13
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=256 finishtime=3.0e13
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=512 finishtime=1.0e14
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=1024 finishtime=3.0e14
REST=`ls DTE_HD_d2l3n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l3n0128 opfreq=2048 finishtime=1.1e15

mpirun -np 8 ../../pion-ng DTE_HD_d2l4n0128_level00_0000.00000000.silo \
  redirect=log_DTE_HD_d2l4n0128 opfreq=64 finishtime=1.0e12
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=128 finishtime=3.0e12
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=256 finishtime=1.0e13
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=512 finishtime=3.0e13
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=1024 finishtime=1.0e14
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=2048 finishtime=3.0e14
REST=`ls DTE_HD_d2l4n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng  $REST \
  redirect=log_DTE_HD_d2l4n0128 opfreq=4096 finishtime=1.1e15

mpirun -np 8 ../../pion-ng DTE_HD_d2l5n0128_level00_0000.00000000.silo \
  redirect=log_DTE_HD_d2l5n0128 opfreq=128 finishtime=1.0e12
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=256 finishtime=3.0e12
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=512 finishtime=1.0e13
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=1024 finishtime=3.0e13
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=2048 finishtime=1.0e14
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=4096 finishtime=3.0e14
REST=`ls DTE_HD_d2l5n0128_level00_0000.* | tail -n1`
mpirun -np 8 ../../pion-ng $REST \
  redirect=log_DTE_HD_d2l5n0128 opfreq=8192 finishtime=1.1e15

exit


