#!/bin/bash
#
# 2011.03.01 JM: testing the diffuse-RT self-shielding algorithm.
# 2011.03.02 JM: comparing serial and parallel runs added.

# call with ./run_RTnodyn_tests.sh $test_dir $code_dir $data_dir
#test_dir=${1}/test_RT_nodyn
#code_dir=$2
#data_dir=$3
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/Diffuse_RT
data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests_201103/RT_diffuse_tests
exe_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
#
# need to define RT_TEST_PROBS
#
#sed -i -e "s/\/\/#define RT_TEST_PROBS/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
#sed -i -e "s/\/\/#define RT_TESTING/#define RT_TESTING/g" ../source/defines/testing_flags.h
echo "GREP";
grep "RT_TEST_PROBS" ../source/defines/testing_flags.h
grep "RT_TESTING"    ../source/defines/testing_flags.h
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"
cd $test_dir

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

echo "RUN_DIFFUSE-RT_TEST.SH:"
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_diffuse_XNXPYNYP.txt silo redirect=msg_temp_
$EXE ICTEST_rtt2D_n32_nh2_diffuse_XNXPYNYP.silo 5 1 \
  outfile=${data_dir}/rtt2D_n32_nh2_diffuse_XNXPYNYP \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_diffuse_XNXPYNYP_

$ICGEN ${test_dir}/pf_rtt3d_n32_nh2_diffuse_XYZ.txt silo redirect=msg_temp_
$EXE ICTEST_rtt3D_n32_nh2_diffuse_XYZ.silo 5 1 \
  outfile=${data_dir}/rtt3D_n32_nh2_diffuse_XYZ \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_diffuse_XYZ_

$ICGEN ${test_dir}/pf_rtt3d_n32_nh2_diffuseXYZ_Ionising.txt silo redirect=msg_temp_
$EXE ICTEST_rtt3D_n32_nh2_diffuseXYZ_Ionising.silo 5 1 cfl=0.39734 checkpt_freq=100000 \
  outfile=${data_dir}/rtt3D_n32_nh2_diffuseXYZ_Ionising \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_diffuseXYZ_Ionising_

cd ${code_dir}/../bin_parallel
make -f Makefile.pllel.code; make -f Makefile.pllel.icgenerator
ICGEN=${exe_dir}/icgen_parallel
EXE=${exe_dir}/gridcode_parallel
cd $test_dir
echo "RUNNING PARALLEL DIFFUSE-RT TESTS"
mpirun -np 4 $ICGEN ${test_dir}/pf_rtt3d_n32_nh2_diffuseXYZ_Ionising.txt \
  silo redirect=msg_temp
#echo mpirun -np 4 $EXE ICTEST_rtt3D_n32_nh2_diffuseXYZ_Ionising_0000.silo 5 1 
#echo    outfile=${data_dir}/rtt3D_n32_nh2_diffuseXYZ_Ionising redirect=msg_temp_np4
mpirun -np 4 $EXE ICTEST_rtt3D_n32_nh2_diffuseXYZ_Ionising_0000.silo 5 1 \
  cfl=0.39734 checkpt_freq=100000 \
  outfile=${data_dir}/rtt3D_n32_nh2_diffuseXYZ_Ionising \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_diffuseXYZ_Ionising_np4
#exit

mpirun -np 8 $ICGEN ${test_dir}/pf_rtt3d_n32_nh2_diffuseXYZ_Ionising.txt \
  silo redirect=msg_temp
mpirun -np 8 $EXE ICTEST_rtt3D_n32_nh2_diffuseXYZ_Ionising_0000.silo 5 1 \
  cfl=0.39734 checkpt_freq=100000 \
  outfile=${data_dir}/rtt3D_n32_nh2_diffuseXYZ_Ionising_np8 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_diffuseXYZ_Ionising_np8

mpirun -np 8 $ICGEN ${test_dir}/pf_rtt3d_n32_nh2_diffuse_XYZ.txt \
  silo redirect=msg_temp
mpirun -np 8 $EXE ICTEST_rtt3D_n32_nh2_diffuse_XYZ_0000.silo 5 1 \
  cfl=0.39734 checkpt_freq=100000 \
  outfile=${data_dir}/rtt3D_n32_nh2_diffuse_XYZ_np8 \
  redirect=${data_dir}/msg_rtt3D_n32_nh2_diffuse_XYZ_np8

## ---------------------------------------------------------------
## -- Compare serial vs. parallel np4 vs. parallel np8 results. --
## ---------------------------------------------------------------
cd ${code_dir}/../analysis/silocompare
make -f Makefile.silocompare

FBASE=rtt3D_n32_nh2_diffuseXYZ_Ionising
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_0000.0     cmp01 2 
./silocompare $data_dir $data_dir ${FBASE}.0 ${FBASE}_np8_0000.0 cmp02 2 




#
# Now switch off RT_TEST_PROBS flag
#
cd $code_dir
#sed -i -e "s/#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
#sed -i -e "s/#define RT_TESTING/\/\/#define RT_TESTING/g" ../source/defines/testing_flags.h
echo "GREP";
grep "RT_TEST_PROBS" ../source/defines/testing_flags.h
grep "RT_TESTING"    ../source/defines/testing_flags.h
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"
cd $test_dir

echo "RUN_RTNODYN_TESTS.SH: FINISHED."
