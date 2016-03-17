#!/bin/bash
#
# Run Blast-Wave test problem in 3D, testing all of the solvers:
# Hybrid, Roe-Conserved (FKJ98 and H-corr viscosity),
# Roe-Primitive, Flux-Vector-Splitting.
#
# 2016.03.17 JM: program written.
#

# Test command-line arguments
if test "$#" -ne 4; then
    echo "UASGE: run_BW3D_parallel.sh test_dir code_dir data_dir resolution"
    exit
fi


# call with ./run_BW3D_parallel.sh $test_dir $code_dir $data_dir $resolution
test_dir=$1   # should be current directory
code_dir=$2
data_dir=$3   # should be sub-directory 'blastwave_axi2d' of the test-results directory.
resolution=$4

# In case it doesn't exist, create the destination directory.
mkdir -p $data_dir

cd ${code_dir}
echo "MAKE IN" ${code_dir}
#bash ./clean.sh
#bash ./compile_code.sh
echo "MAKE FINISHED"
cp ../icgen_parallel ../pion_parallel ${test_dir}
cd ${test_dir}

ICGEN="mpirun -np 16 ./icgen_parallel"
PION="mpirun -np 16 ./pion_parallel"

${ICGEN} ${test_dir}/params_BWcrt3D_Octant_NR${resolution}.txt   silo redirect=tmp_NR${resolution}_

# remove pre-existing simulation outputs from this test-run.
rm ${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_*.*.silo

# 
# First run some short simulations to make sure the early expansion is ok.
#
echo "Running"  ${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_HYB_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_HYB_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
  opfreq_time=1.578e11 finishtime=3.156e10
exit
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RPV_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RPV_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
  opfreq_time=1.578e11 finishtime=3.156e10
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_FVS_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_FVS_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
  opfreq_time=1.578e11 finishtime=3.156e10
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RCV_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RCV_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
  opfreq_time=1.578e11 finishtime=3.156e10
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RCV_Hcorr \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RCV_Hcorr   \
  cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
  opfreq_time=1.578e11 finishtime=3.156e10
#
#wait
#exit
#
# Then run simulations for the full 50kyr to make sure they all
# converge to roughly the same answer.
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_HYB_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_HYB_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=3 \
  opfreq_time=1.578e11 finishtime=1.578e12
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RPV_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RPV_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=5 \
  opfreq_time=1.578e11 finishtime=1.578e12
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_FVS_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_FVS_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=6 \
  opfreq_time=1.578e11 finishtime=1.578e12
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RCV_FKJav01 \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RCV_FKJav01 \
  cfl=0.3 AVtype=1 EtaVisc=0.1 solver=4 \
  opfreq_time=1.578e11 finishtime=1.578e12
#
${PION} IC_BWcrt3D_Octant_NR${resolution}_0000.silo 5 1 \
  outfile=${data_dir}/BWcrt3Dpll_Octant_NR${resolution}_RCV_Hcorr \
  redirect=${data_dir}/msg_BWcrt3Dpll_Octant_NR${resolution}_RCV_Hcorr  \
  cfl=0.3 AVtype=3 EtaVisc=0.0 solver=4 \
  opfreq_time=1.578e11 finishtime=1.578e12
#
#wait
#exit


# --------------------------------------------------------------------
# Now I should have a sequence for each simulation of:
# 0: initial conditions
# 1: Results at  1 kyr (3.156e10)
# 2: Results at 5 kyr
# 3: Results at 10 kyr
# ...
# 11: Results at 50 kyr (1.578e12)
#
# maybe: checkpoint file 1 .9999998.silo
# maybe: checkpoint file 2 .9999999.silo
#
# So let's get rid of the checkpoint files
#
rm ${data_dir}/BWcrt3Dpll_Octant*_00*.999999*.silo
#
# Now we can make visit images, and files 1 and 11 will be the ones of
# interest.  I can tile them for each resolution.
# --------------------------------------------------------------------
# 

