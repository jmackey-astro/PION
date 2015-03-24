#!/bin/bash


test_dir=$1
exe_dir=$2
data_dir=$3

#cd ${code_dir}
# test for main_serial and icgen ???

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

###############################################################################
## 3D runs for stromgen spheres with no dynamics, testing the raytracer.
## 23/11/2009
###############################################################################
##############################################

##########################
## 3D no recombinations ##
##########################
$ICGEN ${test_dir}/pf_rtt3d_n32_nh1_norec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n32_nh2_norec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n64_nh1_norec.txt silo redirect=${data_dir}/msg_temp_
$ICGEN ${test_dir}/pf_rtt3d_n64_nh2_norec.txt silo redirect=${data_dir}/msg_temp_

$EXE ICTEST_rtt3D_n32_nh1_norec.silo 5 1 cfl=0.32 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt010_ &
$EXE ICTEST_rtt3D_n32_nh1_norec.silo 5 1 cfl=0.032 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt100_ &
#$EXE ICTEST_rtt3D_n32_nh1_norec.silo 5 1 cfl=0.0064 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt500_ &

$EXE ICTEST_rtt3D_n32_nh2_norec.silo 5 1 cfl=0.32 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt010_ &
$EXE ICTEST_rtt3D_n32_nh2_norec.silo 5 1 cfl=0.032 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt100_ &
#$EXE ICTEST_rtt3D_n32_nh2_norec.silo 5 1 cfl=0.0064 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt500_ &
wait

$EXE ICTEST_rtt3D_n64_nh1_norec.silo 5 1 cfl=0.64 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt010_ &
$EXE ICTEST_rtt3D_n64_nh1_norec.silo 5 1 cfl=0.064 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt100_ &
#$EXE ICTEST_rtt3D_n64_nh1_norec.silo 5 1 cfl=0.0128 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt500_ &

$EXE ICTEST_rtt3D_n64_nh2_norec.silo 5 1 cfl=0.64 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt010_ &
$EXE ICTEST_rtt3D_n64_nh2_norec.silo 5 1 cfl=0.064 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt100_ &
#$EXE ICTEST_rtt3D_n64_nh2_norec.silo 5 1 cfl=0.0128 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt500_ &
wait


####################################

echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
#make -f Makefile.plotradius clean
make -j4 -f Makefile.plotradius

###############
## 3D Models ##
###############
./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt010  ${data_dir}/rtt3D_n32_nh1_norec_dt010  0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt100  ${data_dir}/rtt3D_n32_nh1_norec_dt100  0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt500  ${data_dir}/rtt3D_n32_nh1_norec_dt500  0 10 5 silo 

./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt010  ${data_dir}/rtt3D_n32_nh2_norec_dt010  0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt100  ${data_dir}/rtt3D_n32_nh2_norec_dt100  0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt500  ${data_dir}/rtt3D_n32_nh2_norec_dt500  0 10 5 silo 
#wait

./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt010 ${data_dir}/rtt3D_n64_nh1_norec_dt010 0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt100 ${data_dir}/rtt3D_n64_nh1_norec_dt100 0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt500 ${data_dir}/rtt3D_n64_nh1_norec_dt500 0 10 5 silo 

./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt010 ${data_dir}/rtt3D_n64_nh2_norec_dt010 0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt100 ${data_dir}/rtt3D_n64_nh2_norec_dt100 0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt500 ${data_dir}/rtt3D_n64_nh2_norec_dt500 0 10 5 silo 
#wait

./ssphere_norec_plots_3d.sh $test_dir $exe_dir $data_dir
echo "all done with 3D no recombinations!"
exit

