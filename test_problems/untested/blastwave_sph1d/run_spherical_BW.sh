#!/bin/bash
#
# Run Blast--Wave test problem in 1D spherical symmetry, testing all
# of the solvers: Hybrid, Roe-Conserved (FKJ98 and H-corr viscosity),
# Roe-Primitive, Flux-Vector-Splitting.
#
#
# 2010.12.08 JM: written.
#

# call with ./run_spherical_BW.sh $test_dir $code_dir $data_dir
test_dir=$1   # should be current directory
code_dir=$2
data_dir=$3   # should be sub-directory 'blastwave_sph1d' of the test-results directory.

#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/blastwave_sph1d
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_results/blastwave_sph1d


# In case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
./compile_code.sh
echo "MAKE SUCEEDED"

./icgen ${test_dir}/pf_sphBW_n128.txt   silo redirect=tmp_
./icgen ${test_dir}/pf_sphBW_n256.txt   silo redirect=tmp_
./icgen ${test_dir}/pf_sphBW_n512.txt   silo redirect=tmp_

#
# N=128
#
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${data_dir}/BW1D_phys_n128_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n128_Hyb_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 &
#
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${data_dir}/BW1D_phys_n128_RPV_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n128_RPV_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=5 &
#
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${data_dir}/BW1D_phys_n128_FVS_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n128_FVS_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=6 &
#
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${data_dir}/BW1D_phys_n128_Roe_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n128_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
#
./main_serial IC_BW1D_phys_n128.silo 5 1 outfile=${data_dir}/BW1D_phys_n128_Roe_Hcorr \
 redirect=${data_dir}/msg_BW1D_phys_n128_Roe_Hcorr   cfl=0.4 AVtype=3 EtaVisc=0.0 solver=4 &
#
#wait
#
# N=256
#
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${data_dir}/BW1D_phys_n256_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n256_Hyb_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 &
#
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${data_dir}/BW1D_phys_n256_RPV_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n256_RPV_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=5 &
#
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${data_dir}/BW1D_phys_n256_FVS_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n256_FVS_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=6 &
#
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${data_dir}/BW1D_phys_n256_Roe_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n256_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
#
./main_serial IC_BW1D_phys_n256.silo 5 1 outfile=${data_dir}/BW1D_phys_n256_Roe_Hcorr \
 redirect=${data_dir}/msg_BW1D_phys_n256_Roe_Hcorr   cfl=0.4 AVtype=3 EtaVisc=0.0 solver=4 &
#
#wait
#
# N=512
#
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${data_dir}/BW1D_phys_n512_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n512_Hyb_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=3 &
#
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${data_dir}/BW1D_phys_n512_RPV_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n512_RPV_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=5 &
#
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${data_dir}/BW1D_phys_n512_FVS_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n512_FVS_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=6 &
#
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${data_dir}/BW1D_phys_n512_Roe_FKJav01 \
 redirect=${data_dir}/msg_BW1D_phys_n512_Roe_FKJav01 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
#
./main_serial IC_BW1D_phys_n512.silo 5 1 outfile=${data_dir}/BW1D_phys_n512_Roe_Hcorr \
 redirect=${data_dir}/msg_BW1D_phys_n512_Roe_Hcorr   cfl=0.4 AVtype=3 EtaVisc=0.0 solver=4 &
#
wait



