#!/bin/bash

##################################
## Overstable Radiative Shocks  ##
## 2009-12-16                   ##
##################################
#
# JM 2009-12-17: 
# This file runs a series of 2D radiative shock tests with a reflecting
# wall as one BC.  The results will then be plotted to get the shock
# position as a function of time to show oscillatory behaviour.
#
# 2011.04.14 JM: Modified a bit to update for new code.
#
# call with ./run_RSH_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}/test_OverstableShock
code_dir=$2
data_dir=$3/RSH
#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/test_OverstableShock
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/RSH

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir
#rm ${data_dir}/*

cd ${code_dir}
echo "MAKE IN" $code_dir
./compile_code.sh
echo "MAKE SUCEEDED"
#
./icgen ${test_dir}/pf_rsh_v100.txt silo
./icgen ${test_dir}/pf_rsh_v120.txt silo
./icgen ${test_dir}/pf_rsh_v130.txt silo
./icgen ${test_dir}/pf_rsh_v140.txt silo
./icgen ${test_dir}/pf_rsh_v150.txt silo
#
./icgen ${test_dir}/pf_rsh_chemV2_v100.txt silo
./icgen ${test_dir}/pf_rsh_chemV2_v120.txt silo
./icgen ${test_dir}/pf_rsh_chemV2_v130.txt silo
./icgen ${test_dir}/pf_rsh_chemV2_v140.txt silo
./icgen ${test_dir}/pf_rsh_chemV2_v150.txt silo
#
#
./main_serial IC_RSH2D_n128_v100_nochem.silo 5 1 redirect=${data_dir}/msg_RSH2D_n128_v100_nochem_ outfile=${data_dir}/RSH2D_n128_v100_nochem cooling=6 &
./main_serial IC_RSH2D_n128_v120_nochem.silo 5 1 redirect=${data_dir}/msg_RSH2D_n128_v120_nochem_ outfile=${data_dir}/RSH2D_n128_v120_nochem cooling=6 &
./main_serial IC_RSH2D_n128_v130_nochem.silo 5 1 redirect=${data_dir}/msg_RSH2D_n128_v130_nochem_ outfile=${data_dir}/RSH2D_n128_v130_nochem cooling=6 &
./main_serial IC_RSH2D_n128_v140_nochem.silo 5 1 redirect=${data_dir}/msg_RSH2D_n128_v140_nochem_ outfile=${data_dir}/RSH2D_n128_v140_nochem cooling=6 &
./main_serial IC_RSH2D_n128_v150_nochem.silo 5 1 redirect=${data_dir}/msg_RSH2D_n128_v150_nochem_ outfile=${data_dir}/RSH2D_n128_v150_nochem cooling=6 &
#wait

./main_serial IC_RSH2D_n256_v100_chemV2.silo   5 1 redirect=${data_dir}/msg_RSH2D_n256_v100_chemV2_ outfile=${data_dir}/RSH2D_n256_v100_chemV2 &
./main_serial IC_RSH2D_n256_v120_chemV2.silo   5 1 redirect=${data_dir}/msg_RSH2D_n256_v120_chemV2_ outfile=${data_dir}/RSH2D_n256_v120_chemV2 &
./main_serial IC_RSH2D_n256_v130_chemV2.silo   5 1 redirect=${data_dir}/msg_RSH2D_n256_v130_chemV2_ outfile=${data_dir}/RSH2D_n256_v130_chemV2 &
./main_serial IC_RSH2D_n256_v140_chemV2.silo   5 1 redirect=${data_dir}/msg_RSH2D_n256_v140_chemV2_ outfile=${data_dir}/RSH2D_n256_v140_chemV2 &
./main_serial IC_RSH2D_n256_v150_chemV2.silo   5 1 redirect=${data_dir}/msg_RSH2D_n256_v150_chemV2_ outfile=${data_dir}/RSH2D_n256_v150_chemV2 &
wait
#
#

#
# Now calculate the shock position as a function of time and plot it!
# Third argument is the base filename for the figures.
#
cd ${test_dir}
bash ./make_overstable_plots.sh ${test_dir} ${data_dir} ${test_dir}/OverstableShock

exit
