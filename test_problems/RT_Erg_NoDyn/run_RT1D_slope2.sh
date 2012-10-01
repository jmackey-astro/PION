#!/bin/bash
#
# 2011.04.21 JM: New file, copied from compare_RT_methods.sh in test_RT_nodyn/
# 2011.06.22 JM: Expanded to test the 2nd order convergence of the RT for different MP timestepping.
# 2011.07.09 JM: Adapted file to work on 1D grid.
# 2011.07.12 JM: Adapted for dynamic model with 1/r^2 density profile.

test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_Erg_NoDyn
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/RT_Erg_NoDyn
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp

cd $test_dir

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

###############################################################################
## 1D runs for stromgren spheres with no dynamics, testing the raytracer
## and the energetics (heating and cooling due to H, He, metals).
##  2011.07.09
###############################################################################

fbase0=rtt_Sph1D_ERG_n128_nh3_slope2
#fbase1=rtt2d_ERG_n100_nh1
#fbase2=rtt2d_ERG_n100_nh2
#fbase3=rtt2d_ERG_n100_nh3


##
## Generate initial conditions files (uniform neutral medium).
##
$ICGEN ${test_dir}/pf_${fbase0}.txt silo redirect=msg_temp_

##
## Run test problems with recombinations, with different xdot timestep criteria
##
SUFFIX=dE10dX01
EXE=./main_serial_${SUFFIX}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

#SUFFIX=dE10dX03
#EXE=./main_serial_${SUFFIX}
#$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
#  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

#SUFFIX=dX03
#EXE=./main_serial_${SUFFIX}
#$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
#  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

SUFFIX=dX05
EXE=./main_serial_${SUFFIX}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

SUFFIX=dX10
EXE=./main_serial_${SUFFIX}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

SUFFIX=dX33
EXE=./main_serial_${SUFFIX}
$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

#SUFFIX=dX03HiAc
#EXE=./main_serial_${SUFFIX}
#$EXE IC_${fbase0}.silo 5 1 cfl=0.3 opfreq=10 optype=6 \
#  checkpt_freq=100000 outfile=${data_dir}/${fbase0}_${SUFFIX} \
#  redirect=${data_dir}/msg_${fbase0}_${SUFFIX}_ & 

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
#make -f Makefile.plotradius clean;
make -j4 -f Makefile.plotradius
#
SUFFIX=dE10dX01
./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
#SUFFIX=dE10dX03
#./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
#SUFFIX=dX03
#./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
#SUFFIX=dX03HiAc
#./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
SUFFIX=dX05
./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
SUFFIX=dX10
./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo
SUFFIX=dX33
./plot_radius ${data_dir}/${fbase0}_${SUFFIX}  ${data_dir}/${fbase0}_${SUFFIX}  0 10 5 silo


#
#
echo "FINISHED 1D W/ ANALYSIS OF *NEW* 1D SIMS W/ RECOMB/DYNAMICS"
#
#echo moving on to generate figs comparing different timestepping criteria.
#cd ${test_dir}
#./ssphere_comparison1D.sh $test_dir $exe_dir $code_dir $data_dir

exit

###----------------------------------------------------------------


