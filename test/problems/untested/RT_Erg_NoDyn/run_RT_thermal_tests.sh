#!/bin/bash
#
# 2011.04.21 JM: New file, copied from compare_RT_methods.sh in test_RT_nodyn/
# 2011.06.22 JM: Expanded to test the 2nd order convergence of the RT for different MP timestepping.
# 2011.07.13 JM: Modified to test different timestepping criteria for 2nd order convergence.

test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/RT_Erg_NoDyn
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp

#./ssphere_comparison.sh $test_dir $exe_dir $code_dir $data_dir
#exit

#cd ${code_dir}
#make -f Makefile.serial.code
#make -f Makefile.serial.icgenerator
#cd $test_dir

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

###############################################################################
## 2D and 3D runs for stromgen spheres with no dynamics, testing the raytracer
## and the energetics (heating and cooling due to H, He, metals).
##  2011.04.21
###############################################################################

fbase0=rtt2d_ERG_n032_nh2
fbase1=rtt2d_ERG_n100_nh1
fbase2=rtt2d_ERG_n100_nh2
fbase3=rtt2d_ERG_n100_nh3


##
## Generate initial conditions files (uniform neutral medium).
##
$ICGEN ${test_dir}/pf_${fbase0}.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_${fbase1}.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_${fbase2}.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_${fbase3}.txt silo redirect=msg_temp_

##
## Run test problems with recombinations, with different xdot timestep criteria
##
#EXE=${test_dir}/main_serial_dE10dX03
#$EXE ICTEST_${fbase0}.silo 5 1 cfl=0.3 opfreq=5 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_dE10dX03 \
#  redirect=${data_dir}/msg_${fbase0}_dE10dX03_ & 
#$EXE ICTEST_${fbase1}.silo 5 1 cfl=0.3 opfreq=5 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase1}_dE10dX03 \
#  redirect=${data_dir}/msg_${fbase1}_dE10dX03_ &
#$EXE ICTEST_${fbase2}.silo 5 1 cfl=0.3 opfreq=5 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase2}_dE10dX03 \
#  redirect=${data_dir}/msg_${fbase2}_dE10dX03_ &
#$EXE ICTEST_${fbase3}.silo 5 1 cfl=0.3 opfreq=5 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase3}_dE10dX03 \
#  redirect=${data_dir}/msg_${fbase3}_dE10dX03_ &
#wait
##
EXE=${test_dir}/main_serial_dX03
$EXE ICTEST_${fbase0}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_dX03 \
  redirect=${data_dir}/msg_${fbase0}_dX03_ & 
$EXE ICTEST_${fbase1}.silo 5 1 cfl=0.3 opfreq=5 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase1}_dX03 \
  redirect=${data_dir}/msg_${fbase1}_dX03_ &
$EXE ICTEST_${fbase2}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase2}_dX03 \
  redirect=${data_dir}/msg_${fbase2}_dX03_ &
$EXE ICTEST_${fbase3}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase3}_dX03 \
  redirect=${data_dir}/msg_${fbase3}_dX03_ &
wait
##
EXE=${test_dir}/main_serial_dX05
$EXE ICTEST_${fbase0}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_dX05 \
  redirect=${data_dir}/msg_${fbase0}_dX05_ & 
$EXE ICTEST_${fbase1}.silo 5 1 cfl=0.3 opfreq=5 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase1}_dX05 \
  redirect=${data_dir}/msg_${fbase1}_dX05_ &
$EXE ICTEST_${fbase2}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase2}_dX05 \
  redirect=${data_dir}/msg_${fbase2}_dX05_ &
$EXE ICTEST_${fbase3}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase3}_dX05 \
  redirect=${data_dir}/msg_${fbase3}_dX05_ &
wait
##
EXE=${test_dir}/main_serial_dX10
$EXE ICTEST_${fbase0}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_dX10 \
  redirect=${data_dir}/msg_${fbase0}_dX10_ & 
$EXE ICTEST_${fbase1}.silo 5 1 cfl=0.3 opfreq=5 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase1}_dX10 \
  redirect=${data_dir}/msg_${fbase1}_dX10_ &
$EXE ICTEST_${fbase2}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase2}_dX10 \
  redirect=${data_dir}/msg_${fbase2}_dX10_ &
$EXE ICTEST_${fbase3}.silo 5 1 cfl=0.3 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase3}_dX10 \
  redirect=${data_dir}/msg_${fbase3}_dX10_ &
wait
##
EXE=${test_dir}/main_serial_dX33
$EXE ICTEST_${fbase0}.silo 5 1 cfl=0.3 opfreq=1 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_dX33 \
  redirect=${data_dir}/msg_${fbase0}_dX33_ & 
$EXE ICTEST_${fbase1}.silo 5 1 cfl=0.3 opfreq=5 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase1}_dX33 \
  redirect=${data_dir}/msg_${fbase1}_dX33_ &
$EXE ICTEST_${fbase2}.silo 5 1 cfl=0.3 opfreq=1 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase2}_dX33 \
  redirect=${data_dir}/msg_${fbase2}_dX33_ &
$EXE ICTEST_${fbase3}.silo 5 1 cfl=0.3 opfreq=1 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase3}_dX33 \
  redirect=${data_dir}/msg_${fbase3}_dX33_ &
wait

########################################
######### DONE WITH SIMULATIONS ########
########################################
###----------------------------------------------------------------
########################################
echo moving on to analysis
cd ${test_dir}/../test_RT_nodyn
echo Now in directory: ${test_dir}/../test_RT_nodyn
pwd
make -f Makefile.plotradius clean; make -j4 -f Makefile.plotradius
#
./plot_radius ${data_dir}/${fbase0}_dE10dX03  ${data_dir}/${fbase0}_dE10dX03  0 5 5 silo
./plot_radius ${data_dir}/${fbase1}_dE10dX03  ${data_dir}/${fbase1}_dE10dX03  0 5 5 silo
./plot_radius ${data_dir}/${fbase2}_dE10dX03  ${data_dir}/${fbase2}_dE10dX03  0 5 5 silo
./plot_radius ${data_dir}/${fbase3}_dE10dX03  ${data_dir}/${fbase3}_dE10dX03  0 5 5 silo
#
./plot_radius ${data_dir}/${fbase0}_dX03  ${data_dir}/${fbase0}_dX03  0 2 5 silo
./plot_radius ${data_dir}/${fbase1}_dX03  ${data_dir}/${fbase1}_dX03  0 5 5 silo
./plot_radius ${data_dir}/${fbase2}_dX03  ${data_dir}/${fbase2}_dX03  0 2 5 silo
./plot_radius ${data_dir}/${fbase3}_dX03  ${data_dir}/${fbase3}_dX03  0 2 5 silo
#
./plot_radius ${data_dir}/${fbase0}_dX05  ${data_dir}/${fbase0}_dX05  0 2 5 silo
./plot_radius ${data_dir}/${fbase1}_dX05  ${data_dir}/${fbase1}_dX05  0 5 5 silo
./plot_radius ${data_dir}/${fbase2}_dX05  ${data_dir}/${fbase2}_dX05  0 2 5 silo
./plot_radius ${data_dir}/${fbase3}_dX05  ${data_dir}/${fbase3}_dX05  0 2 5 silo
#
./plot_radius ${data_dir}/${fbase0}_dX10  ${data_dir}/${fbase0}_dX10  0 2 5 silo
./plot_radius ${data_dir}/${fbase1}_dX10  ${data_dir}/${fbase1}_dX10  0 5 5 silo
./plot_radius ${data_dir}/${fbase2}_dX10  ${data_dir}/${fbase2}_dX10  0 2 5 silo
./plot_radius ${data_dir}/${fbase3}_dX10  ${data_dir}/${fbase3}_dX10  0 2 5 silo
#
./plot_radius ${data_dir}/${fbase0}_dX33  ${data_dir}/${fbase0}_dX33  0 1 5 silo
./plot_radius ${data_dir}/${fbase1}_dX33  ${data_dir}/${fbase1}_dX33  0 5 5 silo
./plot_radius ${data_dir}/${fbase2}_dX33  ${data_dir}/${fbase2}_dX33  0 1 5 silo
./plot_radius ${data_dir}/${fbase3}_dX33  ${data_dir}/${fbase3}_dX33  0 1 5 silo
#
#
echo "FINISHED 2D W/ ANALYSIS OF *NEW* 2D SIMS W/ AND W/O RECOMBINATIONS"
#
echo moving on to generate figs comparing different timestepping criteria.
cd ${test_dir}
./ssphere_comparison.sh $test_dir $exe_dir $code_dir $data_dir

exit

###----------------------------------------------------------------


##
## OLD C2RAY UPDATE RUNS 
##
## 2D runs with no dynamics, no recombinations
#
# nh=10 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt100_ 
$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.006534 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt100_ 
$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=0.02 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt500_ 

#
# nh=100 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt100_ 
$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.006534 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt100_ 
$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=0.02 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt500_ 

#
# nh=1000 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh3_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt100_ 
$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.006534 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh3_norec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt100_ 
$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=0.02 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt500_ 


################################################
## 2D runs with no dynamics, with recombinations
#
# nh=10 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=3.9734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=0.39734 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt100_ 
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=0.0794676 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=12.162162 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=1.2162162 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt100_ 
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=0.24324324 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt500_ 

#
# nh=100 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010_ 
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.039734 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt100_ 
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.00794676 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt500_ 

$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=1.2162162 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt010_ 
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=0.12162162 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt100_ 
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=0.024324324 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt500_ 

#############################################################################
######### DONE WITH SIMULATIONS ########
########################################
echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
make -f Makefile.plotradius clean; make -f Makefile.plotradius
#
# nh=10 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt010  ${data_dir}/rtt2D_n32_nh1_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt100  ${data_dir}/rtt2D_n32_nh1_norec_dt100  0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt500  ${data_dir}/rtt2D_n32_nh1_norec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt010 ${data_dir}/rtt2D_n100_nh1_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt100 ${data_dir}/rtt2D_n100_nh1_norec_dt100 0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt500 ${data_dir}/rtt2D_n100_nh1_norec_dt500 0 10 5 silo
#
# nh=100 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt010  ${data_dir}/rtt2D_n32_nh2_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt100  ${data_dir}/rtt2D_n32_nh2_norec_dt100  0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt500  ${data_dir}/rtt2D_n32_nh2_norec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt010 ${data_dir}/rtt2D_n100_nh2_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt100 ${data_dir}/rtt2D_n100_nh2_norec_dt100 0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt500 ${data_dir}/rtt2D_n100_nh2_norec_dt500 0 10 5 silo
#
# nh=1000 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt010  ${data_dir}/rtt2D_n32_nh3_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt100  ${data_dir}/rtt2D_n32_nh3_norec_dt100  0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt500  ${data_dir}/rtt2D_n32_nh3_norec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt010 ${data_dir}/rtt2D_n100_nh3_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt100 ${data_dir}/rtt2D_n100_nh3_norec_dt100 0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt500 ${data_dir}/rtt2D_n100_nh3_norec_dt500 0 10 5 silo

#################
## 2D Analysis ##
#################
# nh=10 per cc
./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt010  ${data_dir}/rtt2D_n32_nh1_rec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt100  ${data_dir}/rtt2D_n32_nh1_rec_dt100  0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_dt500  ${data_dir}/rtt2D_n32_nh1_rec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt010 ${data_dir}/rtt2D_n100_nh1_rec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt100 ${data_dir}/rtt2D_n100_nh1_rec_dt100 0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_dt500 ${data_dir}/rtt2D_n100_nh1_rec_dt500 0 10 5 silo
# nh=100 per cc
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt010  ${data_dir}/rtt2D_n32_nh2_rec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt100  ${data_dir}/rtt2D_n32_nh2_rec_dt100  0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_dt500  ${data_dir}/rtt2D_n32_nh2_rec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt010 ${data_dir}/rtt2D_n100_nh2_rec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt100 ${data_dir}/rtt2D_n100_nh2_rec_dt100 0 2  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_dt500 ${data_dir}/rtt2D_n100_nh2_rec_dt500 0 10 5 silo
#
echo "FINISHED 2D W/ ANALYSIS OF OLD 2D SIMS W/ AND W/O RECOMBINATIONS"

