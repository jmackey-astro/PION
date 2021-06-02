#!/bin/bash
#
# 2011.02.28 JM: modified a bit for updated code and new directories.

# call with ./run_RTnodyn_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}/test_RT_nodyn
code_dir=$2
data_dir=${3}/RT_static_tests
plleldir=${data_dir}/../RT_pllel_static_tests
exe_dir=${code_dir}/../bin

#code_dir=/vol/aibn129/aibn129_1/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#test_dir=/vol/aibn129/aibn129_1/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_RT_nodyn
#data_dir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_tests
#plleldir=/vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_pllel_static_tests
#exe_dir=/vol/aibn129/aibn129_1/jmackey/active/projects/uniform_grid_code/trunk/bin

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir
mkdir $plleldir

cd ${code_dir}
echo "MAKE IN" $code_dir
#
# need to define RT_TEST_PROBS
#
sed -i -e "s/\/\/#define RT_TEST_PROBS/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
echo "GREP"; grep "RT_TEST_PROBS" ../source/defines/testing_flags.h
./compile_code.sh
cd ../bin_parallel
./compile_code.sh
echo "MAKE FINISHED"
cd $test_dir

echo "RUN_RTNODYN_TESTS.SH: 2D MODELS WITH NO RECOMBINATIONS"
./run_2d_norec.sh $test_dir $exe_dir $data_dir
./run_2d_norec_pllel.sh $test_dir $exe_dir $plleldir $data_dir
./ssphere_norec_plots_2d.sh $test_dir $exe_dir $data_dir

echo "RUN_RTNODYN_TESTS.SH: 3D MODELS WITH NO RECOMBINATIONS"
./run_3d_norec.sh $test_dir $exe_dir $data_dir
./run_3d_norec_pllel.sh $test_dir $exe_dir $plleldir $data_dir
./ssphere_norec_plots_3d.sh $test_dir $exe_dir $data_dir

echo "RUN_RTNODYN_TESTS.SH: 2D MODELS WITH RECOMBINATIONS"
./run_2d_recomb.sh $test_dir $exe_dir $data_dir
./ssphere_recomb_plots_2d.sh $test_dir $exe_dir $data_dir

echo "RUN_RTNODYN_TESTS.SH: 3D MODELS WITH RECOMBINATIONS"
./run_3d_recomb.sh $test_dir $exe_dir $data_dir
./ssphere_recomb_plots_3d.sh $test_dir $exe_dir $data_dir

#
# Now switch off RT_TEST_PROBS flag
#
cd $code_dir
sed -i -e "s/#define RT_TEST_PROBS/\/\/#define RT_TEST_PROBS/g" ../source/defines/testing_flags.h
echo "GREP"; grep "RT_TEST_PROBS" ../source/defines/testing_flags.h
#./compile_code.sh
#cd ../bin_parallel
#./compile_code.sh
#echo "MAKE FINISHED"
cd $test_dir

echo "RUN_RTNODYN_TESTS.SH: FINISHED."
