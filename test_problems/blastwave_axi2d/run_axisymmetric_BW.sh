#!/bin/bash
#
# Run Blast--Wave test problem in 2D axi-symmetry, testing all
# of the solvers: Hybrid, Roe-Conserved (FKJ98 and H-corr viscosity),
# Roe-Primitive, Flux-Vector-Splitting.
#
#
# 2010.12.09 JM: written.
#

# call with ./run_axisymmetric_BW.sh $test_dir $code_dir $data_dir
test_dir=$1   # should be current directory
code_dir=$2
data_dir=$3   # should be sub-directory 'blastwave_axi2d' of the test-results directory.

#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems/blastwave_axi2d
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_results/blastwave_axi2d


# In case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
./compile_code.sh
echo "MAKE FINISHED"

./icgen ${test_dir}/pf_axi2dBW_HalfPlane_NR064.txt   silo redirect=tmp_
./icgen ${test_dir}/pf_axi2dBW_HalfPlane_NR128.txt   silo redirect=tmp_
./icgen ${test_dir}/pf_axi2dBW_HalfPlane_NR256.txt   silo redirect=tmp_
#./icgen ${test_dir}/pf_axi2dBW_HalfPlane_NR512.txt   silo redirect=tmp_

#
# First run some short simulations to make sure the early expansion is ok.
#
# N=064
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_RPV_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_FVS_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Roe_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Roe_Hcorr \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
#wait
#
# N=128
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_RPV_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_FVS_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Roe_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Roe_Hcorr \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
opfreq_time=1.0e20 finishtime=3.16e10 &
#
#wait
##
## N=256
##
#./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Hyb_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_RPV_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_FVS_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Roe_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Roe_Hcorr \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#wait
##
## N=512
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Hyb_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_RPV_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_FVS_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Roe_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Roe_Hcorr \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
#opfreq_time=1.0e20 finishtime=3.16e10 &
##
wait

#
# Then run simulations for the full 50kyr to make sure they all converge to roughly the same answer.
#
# N=064
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_RPV_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_FVS_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Roe_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n064.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n064_Roe_Hcorr \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n064_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
opfreq_time=1.0e20  &
#
#wait
#
# N=128
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_RPV_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_FVS_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Roe_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
opfreq_time=1.0e20  &
#
./main_serial IC_BWaxi2D_HalfPlane_n128.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n128_Roe_Hcorr \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n128_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
opfreq_time=1.0e20  &
#
#wait
##
## N=256
##
./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Hyb_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
opfreq_time=1.0e20  &
##
./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_RPV_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
opfreq_time=1.0e20  &
##
./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_FVS_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
opfreq_time=1.0e20  &
##
./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Roe_FKJav01 \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
opfreq_time=1.0e20  &
##
./main_serial IC_BWaxi2D_HalfPlane_n256.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n256_Roe_Hcorr \
 redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n256_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
opfreq_time=1.0e20  &
##
#wait
##
## N=512
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Hyb_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Hyb_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
#opfreq_time=1.0e20  &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_RPV_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_RPV_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
#opfreq_time=1.0e20  &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_FVS_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_FVS_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
#opfreq_time=1.0e20  &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Roe_FKJav01 \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Roe_FKJav01 cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
#opfreq_time=1.0e20  &
##
#./main_serial IC_BWaxi2D_HalfPlane_n512.silo 5 1 outfile=${data_dir}/BWaxi2D_HalfPlane_n512_Roe_Hcorr \
# redirect=${data_dir}/msg_BWaxi2D_HalfPlane_n512_Roe_Hcorr   cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
#opfreq_time=1.0e20  &
##
wait

# --------------------------------------------------------------------
# Now I should have a sequence for each simulation of:
# 0: initial conditions
# 1: Results at  1 kyr (3.16e10)
# 2: Results at 50 kyr (1.58e12)
# 3: checkpoint file 1 .9999998.silo
# 4: checkpoint file 2 .9999999.silo
#
# So let's get rid of the last two.
#
rm ${data_dir}/BWaxi2D_HalfPlane*.999999*.silo
#
# Now we can make visit images, and files 1 and 2 will be the ones of
# interest.  I can tile them for each resolution.
# --------------------------------------------------------------------
# 
