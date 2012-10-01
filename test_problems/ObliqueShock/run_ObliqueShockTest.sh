#!/bin/bash
#
# Run oblique shock test.
#
# 2010.12.07 JM: written.
#

# call with ./run_ObliqueShockTest.sh $test_dir $code_dir $data_dir
test_dir=${1}/ObliqueShock
code_dir=$2
data_dir=$3
#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/test_FieldLoop
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_results


# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"

./icgen ${test_dir}/pf_oblique_shock_M25.txt   silo redirect=tmp_
./icgen ${test_dir}/pf_oblique_shock_M40.txt   silo redirect=tmp_

./main_serial IC_ObliqueM40.silo 5 1 outfile=${data_dir}/ObliqueM40_Roe_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM40_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
./main_serial IC_ObliqueM40.silo 5 1 outfile=${data_dir}/ObliqueM40_Roe_Hcorr \
 redirect=${data_dir}/msg_ObliqueM40_Roe_Hcorr   cfl=0.2 AVtype=3 EtaVisc=0.1 solver=4 &

./main_serial IC_ObliqueM25.silo 5 1 outfile=${data_dir}/ObliqueM25_Roe_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM25_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
./main_serial IC_ObliqueM25.silo 5 1 outfile=${data_dir}/ObliqueM25_Roe_Hcorr \
 redirect=${data_dir}/msg_ObliqueM25_Roe_Hcorr   cfl=0.2 AVtype=3 EtaVisc=0.1 solver=4 &

wait
exit

