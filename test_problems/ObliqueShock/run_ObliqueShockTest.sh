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
ICGEN=$4
PION=$5

#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/test_FieldLoop
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_results


# Just in case it doesn't exist, create the destination directory.
mkdir -p $data_dir
echo $ICGEN
echo $PION
rm ${data_dir}/ObliqueM*.silo


mpirun -np 4 ${ICGEN} ${test_dir}/params_oblique_shock_M25.txt   silo redirect=tmp_


mpirun -np 4 ${PION} IC_ObliqueM25_0000.silo 5 1 outfile=${data_dir}/ObliqueM25_Roe_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM25_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4
mpirun -np 4 ${PION} IC_ObliqueM25_0000.silo 5 1 outfile=${data_dir}/ObliqueM25_Roe_Hcorr \
 redirect=${data_dir}/msg_ObliqueM25_Roe_Hcorr   cfl=0.4 AVtype=3 EtaVisc=0.1 solver=4
mpirun -np 4 ${PION} IC_ObliqueM25_0000.silo 5 1 outfile=${data_dir}/ObliqueM25_FVS_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM25_FVS_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=6

mpirun -np 4 ${ICGEN} ${test_dir}/params_oblique_shock_M40.txt   silo redirect=tmp_
mpirun -np 4 ${PION} IC_ObliqueM40_0000.silo 5 1 outfile=${data_dir}/ObliqueM40_Roe_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM40_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4
mpirun -np 4 ${PION} IC_ObliqueM40_0000.silo 5 1 outfile=${data_dir}/ObliqueM40_Roe_Hcorr \
 redirect=${data_dir}/msg_ObliqueM40_Roe_Hcorr   cfl=0.4 AVtype=3 EtaVisc=0.1 solver=4
mpirun -np 4 ${PION} IC_ObliqueM40_0000.silo 5 1 outfile=${data_dir}/ObliqueM40_FVS_FKJav01 \
 redirect=${data_dir}/msg_ObliqueM40_FVS_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=6


rm tmp_*

exit

