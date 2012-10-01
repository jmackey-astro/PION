#!/bin/bash
#
# 2011.04.17 JM: New file, copied from run_2d_norec.sh
# 2011.04.18 JM: Expanded a lot, does recombinations runs also.
# 2011.06.22 JM: Expanded more; suffix for output files is now an input option.

SUFFIX=$1

test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/test_RT_nodyn
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp


cd ${code_dir}
#./clean.sh
#make -j 4 -f Makefile.serial.code
#make -j 4 -f Makefile.serial.icgenerator
cd $test_dir

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial
#EXE=./main_serial_$SUFFIX

###############################################################################
## 2D and 3D runs for stromgen spheres with no dynamics, testing the raytracer.
##  'NEW' refers to the new raytracer, and other runs use the C2Ray-type update.
##  2011.04.18
###############################################################################


##### ----------------------------------------------------------------------- #####
##### ----------------------------------------------------------------------- #####

#####################################
### NOW THE NEW RAYTRACER.  IT HAS AUTOMATIC TIMESTEP-LIMITING, 
### SO WE ONLY NEED TO RUN ONE OF EACH MODEL.
####################33
##
## First no recombinations:
##
$ICGEN ${test_dir}/pf_rtt2d_NEW_n32_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n32_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n32_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n100_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n100_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n100_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n32_nh1_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n100_nh1_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n32_nh2_rec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_NEW_n100_nh2_rec.txt silo redirect=msg_temp_


$EXE ICTEST_rtt2D_NEW_n32_nh1_norec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n32_nh2_norec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n32_nh3_norec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n100_nh1_norec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n100_nh2_norec.silo 5 1 cfl=0.3 opfreq=25 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n100_nh3_norec.silo 5 1 cfl=0.3 opfreq=50 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_${SUFFIX}_ &
wait

##
## Now with recombinations
##
$EXE ICTEST_rtt2D_NEW_n32_nh1_rec.silo 5 1 cfl=0.3 opfreq=05 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n100_nh1_rec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n32_nh2_rec.silo 5 1 cfl=0.3 opfreq=05 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_${SUFFIX}_ &

$EXE ICTEST_rtt2D_NEW_n100_nh2_rec.silo 5 1 cfl=0.3 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_${SUFFIX} \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_${SUFFIX}_ &
wait

########################################
######### DONE WITH SIMULATIONS ########
########################################
echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
make -f Makefile.plotradius clean; make -j 4 -f Makefile.plotradius
#
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_${SUFFIX}  ${data_dir}/rtt2D_n32_nh1_norec_${SUFFIX}  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_${SUFFIX} ${data_dir}/rtt2D_n100_nh1_norec_${SUFFIX} 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_${SUFFIX}  ${data_dir}/rtt2D_n32_nh2_norec_${SUFFIX}  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_${SUFFIX} ${data_dir}/rtt2D_n100_nh2_norec_${SUFFIX} 0 25 5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_${SUFFIX}  ${data_dir}/rtt2D_n32_nh3_norec_${SUFFIX}  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_${SUFFIX} ${data_dir}/rtt2D_n100_nh3_norec_${SUFFIX} 0 50 5 silo

./plot_radius ${data_dir}/rtt2D_n32_nh1_rec_${SUFFIX}  ${data_dir}/rtt2D_n32_nh1_rec_${SUFFIX}  0 05 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_rec_${SUFFIX} ${data_dir}/rtt2D_n100_nh1_rec_${SUFFIX} 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_rec_${SUFFIX}  ${data_dir}/rtt2D_n32_nh2_rec_${SUFFIX}  0 05 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_rec_${SUFFIX} ${data_dir}/rtt2D_n100_nh2_rec_${SUFFIX} 0 10 5 silo
#
echo "FINISHED 2D W/ ANALYSIS OF *NEW* 2D SIMS W/ AND W/O RECOMBINATIONS"
#

./ssphere_comparison.sh $test_dir $exe_dir $code_dir $data_dir $SUFFIX

exit

##### ----------------------------------------------------------------------- #####
##### ----------------------------------------------------------------------- #####

##
## OLD C2RAY UPDATE RUNS 
##
## 2D runs with no dynamics, no recombinations
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n32_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh3_norec.txt silo redirect=msg_temp_

#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt100_ &
$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt500_ &
#wait
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt100_ &
$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt500_ &
#wait
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt100_ &
$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt500_ &
wait

#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt100_ &
$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt500_ &
#wait
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt100_ &
$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt500_ &
#wait
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt100_ &
$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt500_ &
wait



################################################
## 2D runs with no dynamics, with recombinations
#
# nh=10 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=3.9734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=0.39734 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt100_ &
$EXE ICTEST_rtt2D_n32_nh1_rec.silo 5 1 cfl=0.0794676 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_rec_dt500_ &
wait

$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=12.162162 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=1.2162162 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt100_ &
$EXE ICTEST_rtt2D_n100_nh1_rec.silo 5 1 cfl=0.24324324 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_rec_dt500_ &
wait

#
# nh=100 per cc
#
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.039734 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt100_ &
$EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.00794676 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt500_ &
wait

$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_rec.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=1.2162162 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=0.12162162 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt100_ &
$EXE ICTEST_rtt2D_n100_nh2_rec.silo 5 1 cfl=0.024324324 opfreq=10 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_rec_dt500 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_rec_dt500_ &
wait

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

##### ----------------------------------------------------------------------- #####
##### ----------------------------------------------------------------------- #####



