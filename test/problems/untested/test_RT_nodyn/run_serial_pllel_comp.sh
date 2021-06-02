#!/bin/bash
#
# 2011.03.01 JM: New file.  Runs some non-dynamics photoionisation tests in serial and
#                with 4 cores and looks for differences between them.

# call with ./run_serial_pllel_comp.sh $serial_dir $pllel_dir $exe_dir $test_dir $data_dir
#serial_dir=$1
# pllel_dir=$2
#   exe_dir=$3
#  test_dir=$4
#  data_dir=$5
serial_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
pllel_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_parallel
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin
test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_RT_nodyn
data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/RT_static_tests

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd $serial_dir
#
# need to define RT_TEST_PROBS
#
sed -i -e "s/\/\/#define RT_TEST_PROBS/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
echo "GREP"; grep "RT_TEST_PROBS" ../source/defines/testing_flags.h

#
# Now compile serial and parallel code versions.
#
./compile_code.sh
cd $pllel_dir
./compile_code.sh
echo "DONE WITH COMPILING"

S_ICGEN=${exe_dir}/icgen_serial
S_EXE=${exe_dir}/main_serial
P_ICGEN=${exe_dir}/icgen_parallel
P_EXE=${exe_dir}/gridcode_parallel


#
# We don't need to run too many models, or for many timesteps.  In principle a single
# should be sufficient to identify any errors.
#
# 2D, low-res, no-recombinations
#
$S_ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
$S_EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010_ 
mpirun -np 4 $P_ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
mpirun -np 4 $P_EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010
mpirun -np 8 $P_ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
mpirun -np 8 $P_EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010_np8 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010_np8


#
# 2D, low-res, w/recombinations
#
$S_ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=msg_temp_
$S_EXE ICTEST_rtt2D_n32_nh2_rec.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010_ 
mpirun -np 4 $P_ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=msg_temp
mpirun -np 4 $P_EXE ICTEST_rtt2D_n32_nh2_rec_0000.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010
mpirun -np 8 $P_ICGEN ${test_dir}/pf_rtt2d_n32_nh2_rec.txt silo redirect=msg_temp
mpirun -np 8 $P_EXE ICTEST_rtt2D_n32_nh2_rec_0000.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_rec_dt010_np8  \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_rec_dt010_np8 

#
# 3D, low-res, no-recombinations
#
$S_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_norec.txt silo redirect=msg_temp_
$S_EXE ICTEST_rtt3D_n32_nh2_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt01_
mpirun -np 4 $P_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_norec.txt silo redirect=msg_temp
mpirun -np 4 $P_EXE ICTEST_rtt3D_n32_nh2_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt010
mpirun -np 8 $P_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_norec.txt silo redirect=msg_temp
mpirun -np 8 $P_EXE ICTEST_rtt3D_n32_nh2_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_norec_dt010_np8  \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_norec_dt010_np8 

#
# 3D, low-res, w/recombinations
#
$S_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_rec.txt silo redirect=msg_temp_
$S_EXE ICTEST_rtt3D_n32_nh2_rec.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt010_ 
mpirun -np 4 $P_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_rec.txt silo redirect=msg_temp
mpirun -np 4 $P_EXE ICTEST_rtt3D_n32_nh2_rec_0000.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt010 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt010
mpirun -np 8 $P_ICGEN ${test_dir}/pf_rtt3d_n32_nh2_rec.txt silo redirect=msg_temp
mpirun -np 8 $P_EXE ICTEST_rtt3D_n32_nh2_rec_0000.silo 5 1 cfl=0.39734 \
  checkpt_freq=100000 outfile=${data_dir}/rtt3D_n32_nh2_rec_dt010_np8  \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_rec_dt010_np8 

## ---------------------------------------------------------------
## -- Compare serial vs. parallel np4 vs. parallel np8 results. --
## ---------------------------------------------------------------
cd ${serial_dir}/../analysis/silocompare
make -f Makefile.silocompare

FBASE=rtt2D_n32_nh1_norec_dt010
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_0000.0     cmp01 2 
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_np8_0000.0 cmp02 2 

FBASE=rtt2D_n32_nh2_rec_dt010
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_0000.0     cmp03 2 
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_np8_0000.0 cmp04 2 

FBASE=rtt3D_n32_nh2_norec_dt010
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_0000.0     cmp05 2 
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_np8_0000.0 cmp06 2 

FBASE=rtt3D_n32_nh2_rec_dt010
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_0000.0     cmp07 2 
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_np8_0000.0 cmp08 2 

#more cmp0*.txt


#
# Now switch off RT_TEST_PROBS flag
#
cd $serial_dir
sed -i -e "s/#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
echo "GREP"; grep "RT_TEST_PROBS" ../source/defines/testing_flags.h
#./compile_code.sh
#echo "MAKE SUCEEDED"
cd $test_dir


