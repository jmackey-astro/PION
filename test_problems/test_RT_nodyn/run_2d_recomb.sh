#!/bin/bash


test_dir=$1
exe_dir=$2
data_dir=$3

#cd ${code_dir}
# test for main_serial and icgen ???

#ICGEN=${exe_dir}/icgen_serial
#EXE=${exe_dir}/main_serial
ICGEN="mpiexec -np 4 ${exe_dir}/icgen_parallel"
EXE="mpiexec -np 4 ${exe_dir}/gridcode_parallel"


###############################################################################
## 2D runs for stromgen spheres with no dynamics, testing the raytracer.
## 23/11/2009
###############################################################################
##############################################

################################################
## 2D runs with no dynamics, with recombinations
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_rec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_rec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n256_nh1_rec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_rec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n256_nh2_rec.txt silo redirect=${data_dir}/msg_temp_


#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n32_nh1_rec_0000.silo 5 1 cfl=3.2664 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh1_rec_0000.silo 5 1 cfl=0.32664 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt100_ 
#$EXE ICTEST_rtt2D_n32_nh1_rec_0000.silo 5 1 cfl=0.065328 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt500_ 
#wait

$EXE ICTEST_rtt2D_n100_nh1_rec_0000.silo 5 1 cfl=10.2075 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh1_rec_0000.silo 5 1 cfl=1.02075 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt100_ 
#$EXE ICTEST_rtt2D_n100_nh1_rec_0000.silo 5 1 cfl=0.20415 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt500_ 
#wait

$EXE ICTEST_rtt2D_n256_nh1_rec_0000.silo 5 1 cfl=26.1312 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_rec_dt010_ 
$EXE ICTEST_rtt2D_n256_nh1_rec_0000.silo 5 1 cfl=2.61312 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_rec_dt100_ 
#$EXE ICTEST_rtt2D_n256_nh1_rec_0000.silo 5 1 cfl=0.522624 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh1_rec_dt500_ 
#wait

#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n32_nh2_rec_0000.silo 5 1 cfl=0.32664 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh2_rec_0000.silo 5 1 cfl=0.032664 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt100_ 
#$EXE ICTEST_rtt2D_n32_nh2_rec_0000.silo 5 1 cfl=0.0065328 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt500_ 
#wait

$EXE ICTEST_rtt2D_n100_nh2_rec_0000.silo 5 1 cfl=1.02075 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh2_rec_0000.silo 5 1 cfl=0.102075 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt100_ 
#$EXE ICTEST_rtt2D_n100_nh2_rec_0000.silo 5 1 cfl=0.020415 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt500_ 
#wait

$EXE ICTEST_rtt2D_n256_nh2_rec_0000.silo 5 1 cfl=2.61312 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_rec_dt010_ 
$EXE ICTEST_rtt2D_n256_nh2_rec_0000.silo 5 1 cfl=0.261312 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_rec_dt100_ 
#$EXE ICTEST_rtt2D_n256_nh2_rec_0000.silo 5 1 cfl=0.0522624 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh2_rec_dt500_ 
#wait

#
# nh=1000 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh3_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh3_rec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_rec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh3_rec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_rec_dt100_ 
#$EXE ICTEST_rtt2D_n32_nh3_rec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh3_rec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh3_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh3_rec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_rec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh3_rec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_rec_dt100_ 
#$EXE ICTEST_rtt2D_n100_nh3_rec_0000.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh3_rec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n256_nh3_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n256_nh3_rec_0000.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_rec_dt010_ 
$EXE ICTEST_rtt2D_n256_nh3_rec_0000.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_rec_dt100_ 
#$EXE ICTEST_rtt2D_n256_nh3_rec_0000.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_rec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh3_rec_dt500_ 

####################################

echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
#make -f Makefile.plotradius clean;
make -j4 -f Makefile.plotradius

#################
## 2D Analysis ##
#################
# nh=10 per cc
./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt010  ${data_dir}/rtt2D_n32_nh1_rec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt100  ${data_dir}/rtt2D_n32_nh1_rec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt500  ${data_dir}/rtt2D_n32_nh1_rec_dt500_0000  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt010 ${data_dir}/rtt2D_n100_nh1_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt100 ${data_dir}/rtt2D_n100_nh1_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt500 ${data_dir}/rtt2D_n100_nh1_rec_dt500_0000 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_rec_dt010 ${data_dir}/rtt2D_n256_nh1_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_rec_dt100 ${data_dir}/rtt2D_n256_nh1_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh1_rec_dt500 ${data_dir}/rtt2D_n256_nh1_rec_dt500_0000 0 10 5 silo

# nh=100 per cc
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt010  ${data_dir}/rtt2D_n32_nh2_rec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt100  ${data_dir}/rtt2D_n32_nh2_rec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt500  ${data_dir}/rtt2D_n32_nh2_rec_dt500_0000  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt010 ${data_dir}/rtt2D_n100_nh2_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt100 ${data_dir}/rtt2D_n100_nh2_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt500 ${data_dir}/rtt2D_n100_nh2_rec_dt500_0000 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh2_rec_dt010 ${data_dir}/rtt2D_n256_nh2_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh2_rec_dt100 ${data_dir}/rtt2D_n256_nh2_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh2_rec_dt500 ${data_dir}/rtt2D_n256_nh2_rec_dt500_0000 0 10 5 silo

# nh=1000 per cc
./plot_radius ${data_dir}/rtt2D_n32_nh3_rec_dt010  ${data_dir}/rtt2D_n32_nh3_rec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_rec_dt100  ${data_dir}/rtt2D_n32_nh3_rec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh3_rec_dt500  ${data_dir}/rtt2D_n32_nh3_rec_dt500_0000  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_rec_dt010 ${data_dir}/rtt2D_n100_nh3_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_rec_dt100 ${data_dir}/rtt2D_n100_nh3_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh3_rec_dt500 ${data_dir}/rtt2D_n100_nh3_rec_dt500_0000 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh3_rec_dt010 ${data_dir}/rtt2D_n256_nh3_rec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh3_rec_dt100 ${data_dir}/rtt2D_n256_nh3_rec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh3_rec_dt500 ${data_dir}/rtt2D_n256_nh3_rec_dt500_0000 0 10 5 silo

#
echo "FINISHED 2D W/RECOMBINATIONS SIMULATIONS AND ANALYSIS"

exit
