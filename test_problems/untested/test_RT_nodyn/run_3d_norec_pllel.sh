#!/bin/bash


test_dir=$1
exe_dir=$2
data_dir=$3
serialdir=$4

#cd ${code_dir}
# test for main_serial and icgen ???

ICGEN="mpiexec -np 4 ${exe_dir}/icgen_parallel"
EXE="mpiexec -np 4 ${exe_dir}/gridcode_parallel"

###############################################################################
## 3D runs for stromgen spheres with no dynamics, testing the raytracer.
## 23/11/2009
###############################################################################
##############################################

##########################
## 3D no recombinations ##
##########################
$ICGEN ${test_dir}/pf_rtt3d_n32_nh1_norec.txt silo redirect=${data_dir}/msg_temp
$ICGEN ${test_dir}/pf_rtt3d_n32_nh2_norec.txt silo redirect=${data_dir}/msg_temp
$ICGEN ${test_dir}/pf_rtt3d_n64_nh1_norec.txt silo redirect=${data_dir}/msg_temp
$ICGEN ${test_dir}/pf_rtt3d_n64_nh2_norec.txt silo redirect=${data_dir}/msg_temp

$EXE ICTEST_rtt3D_n32_nh1_norec_0000.silo 5 1 cfl=0.32 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt010
$EXE ICTEST_rtt3D_n32_nh1_norec_0000.silo 5 1 cfl=0.032 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt100
#$EXE ICTEST_rtt3D_n32_nh1_norec_0000.silo 5 1 cfl=0.0064 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh1_norec_dt500

$EXE ICTEST_rtt3D_n32_nh2_norec_0000.silo 5 1 cfl=0.32 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt010
$EXE ICTEST_rtt3D_n32_nh2_norec_0000.silo 5 1 cfl=0.032 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt100
#$EXE ICTEST_rtt3D_n32_nh2_norec_0000.silo 5 1 cfl=0.0064 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt500
#exit
$EXE ICTEST_rtt3D_n64_nh1_norec_0000.silo 5 1 cfl=0.64 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt010
$EXE ICTEST_rtt3D_n64_nh1_norec_0000.silo 5 1 cfl=0.064 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt100
#$EXE ICTEST_rtt3D_n64_nh1_norec_0000.silo 5 1 cfl=0.0128 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh1_norec_dt500

$EXE ICTEST_rtt3D_n64_nh2_norec_0000.silo 5 1 cfl=0.64 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt010
$EXE ICTEST_rtt3D_n64_nh2_norec_0000.silo 5 1 cfl=0.064 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt100
#$EXE ICTEST_rtt3D_n64_nh2_norec_0000.silo 5 1 cfl=0.0128 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n64_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt3D_n64_nh2_norec_dt500


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
./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt010  ${data_dir}/rtt3D_n32_nh1_norec_dt010_0000  0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt100  ${data_dir}/rtt3D_n32_nh1_norec_dt100_0000  0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n32_nh1_norec_dt500  ${data_dir}/rtt3D_n32_nh1_norec_dt500_0000  0 10 5 silo 
./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt010  ${data_dir}/rtt3D_n32_nh2_norec_dt010_0000  0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt100  ${data_dir}/rtt3D_n32_nh2_norec_dt100_0000  0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n32_nh2_norec_dt500  ${data_dir}/rtt3D_n32_nh2_norec_dt500_0000  0 10 5 silo 
#wait

./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt010 ${data_dir}/rtt3D_n64_nh1_norec_dt010_0000 0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt100 ${data_dir}/rtt3D_n64_nh1_norec_dt100_0000 0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n64_nh1_norec_dt500 ${data_dir}/rtt3D_n64_nh1_norec_dt500_0000 0 10 5 silo 
./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt010 ${data_dir}/rtt3D_n64_nh2_norec_dt010_0000 0 1  5 silo 
./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt100 ${data_dir}/rtt3D_n64_nh2_norec_dt100_0000 0 2  5 silo 
#./plot_radius ${data_dir}/rtt3D_n64_nh2_norec_dt500 ${data_dir}/rtt3D_n64_nh2_norec_dt500_0000 0 10 5 silo 
#wait

./ssphere_norec_plots_3d.sh $test_dir $exe_dir $data_dir
cp photoncons_norec3d* $data_dir

#############################
### silo file comparison ####
#############################
echo "############ Calculating DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"
cd ../../analysis/silocompare
make -j4 -f Makefile.silocompare
./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh1_norec_dt010 rtt3D_n32_nh1_norec_dt010 cmp01 2
./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh2_norec_dt010 rtt3D_n32_nh2_norec_dt010 cmp02 2
./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh1_norec_dt100 rtt3D_n32_nh1_norec_dt100 cmp03 2
./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh2_norec_dt100 rtt3D_n32_nh2_norec_dt100 cmp04 2
#./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh1_norec_dt500 rtt3D_n32_nh1_norec_dt500 cmp05 2
#./silocompare ${data_dir} ${serialdir} rtt3D_n32_nh2_norec_dt500 rtt3D_n32_nh2_norec_dt500 cmp06 2

./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh1_norec_dt010 rtt3D_n64_nh1_norec_dt010 cmp11 2
./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh2_norec_dt010 rtt3D_n64_nh2_norec_dt010 cmp12 2
./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh1_norec_dt100 rtt3D_n64_nh1_norec_dt100 cmp13 2
./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh2_norec_dt100 rtt3D_n64_nh2_norec_dt100 cmp14 2
#./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh1_norec_dt500 rtt3D_n64_nh1_norec_dt500 cmp15 2
#./silocompare ${data_dir} ${serialdir} rtt3D_n64_nh2_norec_dt500 rtt3D_n64_nh2_norec_dt500 cmp16 2

echo "############ PRINTING DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"
grep "^[0-5]" msg_cmp*fo.txt
rm msg_cmp*
cd -
echo "############ PRINTING DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"

echo "all done with 3D no recombinations!"
exit

