#!/bin/bash

echo "RUNNING 3D RT-TESTING MODELS WITH RECOMBINATIONS AND NO DYNAMICS"

test_dir=$1
exe_dir=$2
data_dir=$3

#cd ${code_dir}
# test for main_serial and icgen ???

ICGEN="mpiexec -np 4 ${exe_dir}/icgen_parallel"
EXE="mpiexec -np 4 ${exe_dir}/gridcode_parallel"

#ICGEN=${exe_dir}/icgen_serial
#EXE=${exe_dir}/main_serial

##########################
## 3D w/ recombinations ##
##########################
$ICGEN ${test_dir}/pf_rtt3d_n32_nh1_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n32_nh2_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n64_nh1_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n64_nh2_rec.txt silo redirect=msg_temp_


$EXE ICTEST_rtt3D_n32_nh1_rec_0000.silo 5 1 cfl=3.2664 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh1_rec_dt010_ 
#$EXE ICTEST_rtt3D_n32_nh1_rec_0000.silo 5 1 cfl=0.32664 opfreq=2 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_rec_dt100 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh1_rec_dt100_ 
#$EXE ICTEST_rtt3D_n32_nh1_rec_0000.silo 5 1 cfl=0.065328 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_rec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh1_rec_dt500_ 

$EXE ICTEST_rtt3D_n32_nh2_rec_0000.silo 5 1 cfl=0.32664 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt010_ 
#$EXE ICTEST_rtt3D_n32_nh2_rec_0000.silo 5 1 cfl=0.032664 opfreq=2 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt100 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt100_ 
#$EXE ICTEST_rtt3D_n32_nh2_rec_0000.silo 5 1 cfl=0.0065328 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt500_ 

$EXE ICTEST_rtt3D_n64_nh1_rec_0000.silo 5 1 cfl=10.2075 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh1_rec_dt010_ 
#$EXE ICTEST_rtt3D_n64_nh1_rec_0000.silo 5 1 cfl=1.02075 opfreq=2 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_rec_dt100 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh1_rec_dt100_ 
#$EXE ICTEST_rtt3D_n64_nh1_rec_0000.silo 5 1 cfl=0.20415 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_rec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh1_rec_dt500_ 

$EXE ICTEST_rtt3D_n64_nh2_rec_0000.silo 5 1 cfl=1.02075 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh2_rec_dt010_ 
#$EXE ICTEST_rtt3D_n64_nh2_rec_0000.silo 5 1 cfl=0.102075 opfreq=2 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_rec_dt100 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh2_rec_dt100_ 
#$EXE ICTEST_rtt3D_n64_nh2_rec_0000.silo 5 1 cfl=0.020415 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_rec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh2_rec_dt500_ 

#################
## 3D Analysis ##
#################
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
#make -f Makefile.plotradius clean
make -f Makefile.plotradius

./plot_radius ${data_dir}/rtt3D_n32_nh1_rec_dt010  ${data_dir}/rtt3D_n32_nh1_rec_dt010_0000  0 1  5 silo
#./plot_radius ${data_dir}/rtt3D_n32_nh1_rec_dt100  ${data_dir}/rtt3D_n32_nh1_rec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt3D_n32_nh1_rec_dt500  ${data_dir}/rtt3D_n32_nh1_rec_dt500_0000  0 10 5 silo

./plot_radius ${data_dir}/rtt3D_n32_nh2_rec_dt010  ${data_dir}/rtt3D_n32_nh2_rec_dt010_0000  0 1  5 silo
#./plot_radius ${data_dir}/rtt3D_n32_nh2_rec_dt100  ${data_dir}/rtt3D_n32_nh2_rec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt3D_n32_nh2_rec_dt500  ${data_dir}/rtt3D_n32_nh2_rec_dt500_0000  0 10 5 silo

./plot_radius ${data_dir}/rtt3D_n64_nh1_rec_dt010 ${data_dir}/rtt3D_n64_nh1_rec_dt010_0000 0 1  5 silo
#./plot_radius ${data_dir}/rtt3D_n64_nh1_rec_dt100 ${data_dir}/rtt3D_n64_nh1_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt3D_n64_nh1_rec_dt500 ${data_dir}/rtt3D_n64_nh1_rec_dt500_0000 0 10 5 silo

./plot_radius ${data_dir}/rtt3D_n64_nh2_rec_dt010 ${data_dir}/rtt3D_n64_nh2_rec_dt010_0000 0 1  5 silo
#./plot_radius ${data_dir}/rtt3D_n64_nh2_rec_dt100 ${data_dir}/rtt3D_n64_nh2_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt3D_n64_nh2_rec_dt500 ${data_dir}/rtt3D_n64_nh2_rec_dt500_0000 0 10 5 silo

echo "FINISHED 3D W/RECOMBINATIONS SIMULATIONS AND ANALYSIS."
 
