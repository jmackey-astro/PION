#!/bin/bash
#
# Compare blast wave test results between serial and parallel code.
#
# 2016.03.18 JM: program written.
#

# Test command-line arguments
if test "$#" -ne 4; then
    echo "UASGE: run_BW3D_parallel.sh test_dir code_dir data_dir resolution"
    exit
fi


# call with ./run_BW3D_parallel.sh $test_dir $cmp_dir $data_dir $resolution
test_dir=$1   # should be current directory
cmp_dir=$2
data_dir=$3   # should be sub-directory 'blastwave_axi2d' of the test-results directory.
resolution=$4

# In case it doesn't exist, create the destination directory.
mkdir -p $data_dir

cd ${cmp_dir}
echo "MAKE IN" ${cmp_dir}
MAKE_UNAME=standard  make -j8 clean
MAKE_UNAME=standard  make -j8
echo "MAKE FINISHED"
cp silocompare ${test_dir}
cd ${test_dir}

# <silocompare> <first-dir> <comp-dir> <first-file> <flt/dbl>
#    <comp-file> <flt/dbl> <outfile> <fabs/plus-minus/L1/L2>
for solver in "HYB_FKJav01" "RPV_FKJav01" "FVS_FKJav01" "RCV_FKJav01" "RCV_Hcorr"; do
  rm BWcrt3Dcomparison_NR${resolution}_${solver}*.txt
  ./silocompare ${data_dir} ${data_dir}  \
    BWcrt3Dser_Octant_NR${resolution}_${solver}.0 DOUBLE      \
    BWcrt3Dpll_Octant_NR${resolution}_${solver}_0000.0 DOUBLE \
    BWcrt3Dcomparison_NR${resolution}_${solver} 2 > tmp_NR${resolution}_${solver}.txt
  echo " *** Checking for L1, L2, Max, errors: ***"
  grep "^[0-9]" BWcrt3Dcomparison_NR${resolution}_${solver}_0.txt
done

