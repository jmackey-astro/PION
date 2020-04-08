#!/bin/bash
#
# Run Blast--Wave test problem in 2D axi-symmetry, testing all
# of the solvers: Hybrid, Roe-Conserved (FKJ98 and H-corr viscosity),
# Roe-Primitive, Flux-Vector-Splitting.
#
#
# 2010.12.09 JM: written.
# 2013.02.07 JM: updated to work with pion v.0.1.
#

# call with ./run_axisymmetric_BW.sh $test_dir $code_dir $data_dir
test_dir=$1   # should be current directory
code_dir=$2
data_dir=$3   # should be sub-directory 'blastwave_axi2d' of the test-results directory.


# In case it doesn't exist, create the destination directory.
mkdir -p $data_dir

cd ${code_dir}
echo "MAKE IN" ${code_dir}
#bash ./clean.sh
bash ./compile_code.sh
echo "MAKE FINISHED"
cp ../icgen_* ../pion_* ${test_dir}
cd ${test_dir}

ICGEN="./icgen-ugs"
PION="./pion-ugs"

${ICGEN} ${test_dir}/params_axi2dBW_HalfPlane_NR064.txt   silo redirect=tmp_NR064_
${ICGEN} ${test_dir}/params_axi2dBW_HalfPlane_NR128.txt   silo redirect=tmp_NR128_
${ICGEN} ${test_dir}/params_axi2dBW_HalfPlane_NR256.txt   silo redirect=tmp_NR256_
${ICGEN} ${test_dir}/params_axi2dBW_HalfPlane_NR512.txt   silo redirect=tmp_NR512_

#
#
for NNN in "n064" "n128" "n256" "n512"
  do

# 
# First run some short simulations to make sure the early expansion is ok.
#
  echo "Running"  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Hyb_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Hyb_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
    opfreq_time=1.0e20 finishtime=3.156e10 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_RPV_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_RPV_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
    opfreq_time=1.0e20 finishtime=3.156e10 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_FVS_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_FVS_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
    opfreq_time=1.0e20 finishtime=3.156e10 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Roe_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Roe_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
    opfreq_time=1.0e20 finishtime=3.156e10 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Roe_Hcorr \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Roe_Hcorr   \
    cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
    opfreq_time=1.0e20 finishtime=3.156e10 &
#
  wait
#exit
#
# Then run simulations for the full 50kyr to make sure they all converge to roughly the same answer.
#
# N=064
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Hyb_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Hyb_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
    opfreq_time=1.0e20 finishtime=1.578e12 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_RPV_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_RPV_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
    opfreq_time=1.0e20 finishtime=1.578e12 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_FVS_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_FVS_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
    opfreq_time=1.0e20 finishtime=1.578e12 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Roe_FKJav01 \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Roe_FKJav01 \
    cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
    opfreq_time=1.0e20 finishtime=1.578e12 &
#
  ${PION} BWaxi2D_HalfPlane_${NNN}.00000000.silo 5 1 \
    outfile=${data_dir}/BWaxi2D_HalfPlane_${NNN}_Roe_Hcorr \
    redirect=${data_dir}/msg_BWaxi2D_HalfPlane_${NNN}_Roe_Hcorr  \
    cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
    opfreq_time=1.0e20 finishtime=1.578e12 &
#
  wait
done
#exit


# --------------------------------------------------------------------
# Now I should have a sequence for each simulation of:
# 0: initial conditions
# 1: Results at  1 kyr (3.156e10)
# 2: Results at 50 kyr (1.578e12)
# Now we can make visit images, and files 1 and 2 will be the ones of
# interest.  I can tile them for each resolution.
# --------------------------------------------------------------------
# 
