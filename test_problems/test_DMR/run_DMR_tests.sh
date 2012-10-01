#!/bin/bash
#
# 2010-09-21 JM: Added FVS solver tests.
#
# call with ./run_DMR_tests.sh $test_dir $code_dir $data_dir

test_dir=$1
code_dir=$2
data_dir=$3
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test_problems
#data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/temp_sims

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
./compile_code.sh
if [ ! -f ../bin/main_serial ] || [ ! -f ../bin/icgen_serial ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi

#
# Double Mach Reflection, same as Stone et al. 2008.
#
echo "DOUBLE MACH REFLECTION: GENERATE ICS"
./icgen ${test_dir}/test_DMR/pf_DMR_n260.txt silo
./icgen ${test_dir}/test_DMR/pf_DMR_n520.txt silo

#./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_EINT_av10 cfl=0.5 \
# redirect=${data_dir}/msg_DMRm10t60_n260_EINT_av10 solver=3 eqntype=9 AVtype=1 EtaVisc=0.1
#exit



#
echo "DOUBLE MACH REFLECTION: RUN LOW RES"
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_Hyb_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_Hyb_av00 solver=3 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_Hyb_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_Hyb_av10 solver=3 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_av00 solver=4 &
#wait
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_av10 solver=4 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RPV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_RPV_av00 solver=5 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RPV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_RPV_av10 solver=5 &
#wait
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_FVS_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_FVS_av00 solver=6 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_FVS_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_FVS_av10 solver=6 &
./main_serial IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_Hcor cfl=0.4 AVtype=3    \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_Hcor solver=4 &
wait

########## TEMP #########
#exit
########## TEMP #########

echo "DOUBLE MACH REFLECTION: RUN HIGH RES"
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_Hyb_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_Hyb_av00 solver=3 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_Hyb_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_Hyb_av10 solver=3 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_av00 solver=4 &
#wait
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_av10 solver=4 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RPV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_RPV_av00 solver=5 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RPV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_RPV_av10 solver=5 &
#wait
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_FVS_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_FVS_av00 solver=6 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_FVS_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_FVS_av10 solver=6 &
./main_serial IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_Hcor cfl=0.4 AVtype=3    \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_Hcor solver=4 &
wait

exit
