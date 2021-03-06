#!/bin/bash
#
# Author:  Jonathan Mackey
#
# 2010-09-21 JM: Added FVS solver tests.
# 2015.03.26 JM: Set up to run in parallel.
#
# call with ./run_DMR_tests.sh $test_dir $code_dir $data_dir

if [ "$1" = "" ] ||  [ "$2" = "" ] ||  [ "$3" = "" ]
then
 echo Usage: $0 test_dir code_dir data_dir
 exit
fi


test_dir=$1
code_dir=$2
data_dir=$3
resolution=$4

# Just in case it doesn't exist, create the destination directory.
mkdir -p $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
bash ./compile_code.sh
if [ ! -f ../pion-ugs ] || [ ! -f ../icgen-ugs ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
  cp ../pion-ugs ../icgen-ugs ${test_dir}
  cd ${test_dir}
fi

##############  TESTING  ##############
#echo "DOUBLE MACH REFLECTION: GENERATE ICS"
#./icgen-ugs ${test_dir}/params_DMR_n260.txt silo
#./icgen-ugs ${test_dir}/params_DMR_n520.txt silo
#echo "DOUBLE MACH REFLECTION: RUN LOW RES"
#./pion-ugs IC_DMRm10t60_n260.silo 5 1 \
# outfile=${data_dir}/DMRm10t60_n260_Hyb_av00 cfl=0.4 artvisc=0   \
# redirect=${data_dir}/msg_DMRm10t60_n260_Hyb_av00 solver=3 &
#./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_Hcor cfl=0.4 AVtype=3 \
# redirect=${data_dir}/msg_DMRm10t60_n260_RCV_Hcor solver=4 &
#./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_Hcor cfl=0.4 AVtype=3    \
# redirect=${data_dir}/msg_DMRm10t60_n520_RCV_Hcor solver=4 &
#wait
#exit
##############  TESTING  ##############



#
# Double Mach Reflection, same as Stone et al. 2008.
#
echo "DOUBLE MACH REFLECTION: GENERATE ICS"
./icgen-ugs ${test_dir}/params_DMR_n${resolution}.txt silo redirect=${data_dir}/iclog_n${resolution}_

echo "DOUBLE MACH REFLECTION: RUN LOW RES"
./pion-ugs DMRm10t60_n${resolution}.00000000.silo 5 1 outfile=${data_dir}/DMRm10t60_n${resolution}_Hyb_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n${resolution}_Hyb_av10 solver=3 &
./pion-ugs DMRm10t60_n${resolution}.00000000.silo 5 1 outfile=${data_dir}/DMRm10t60_n${resolution}_RCV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n${resolution}_RCV_av10 solver=4 &
./pion-ugs DMRm10t60_n${resolution}.00000000.silo 5 1 outfile=${data_dir}/DMRm10t60_n${resolution}_RPV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n${resolution}_RPV_av10 solver=5 &
./pion-ugs DMRm10t60_n${resolution}.00000000.silo 5 1 outfile=${data_dir}/DMRm10t60_n${resolution}_FVS_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n${resolution}_FVS_av10 solver=6 &
./pion-ugs DMRm10t60_n${resolution}.00000000.silo 5 1 outfile=${data_dir}/DMRm10t60_n${resolution}_RCV_Hcor cfl=0.4 AVtype=3    \
 redirect=${data_dir}/msg_DMRm10t60_n${resolution}_RCV_Hcor solver=4 &
wait
exit














#
echo "DOUBLE MACH REFLECTION: RUN LOW RES"
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_Hyb_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_Hyb_av00 solver=3 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_Hyb_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_Hyb_av10 solver=3 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_av00 solver=4 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_av10 solver=4 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RPV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_RPV_av00 solver=5 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RPV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_RPV_av10 solver=5 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_FVS_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n260_FVS_av00 solver=6 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_FVS_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n260_FVS_av10 solver=6 &
./pion-ugs IC_DMRm10t60_n260.silo 5 1 outfile=${data_dir}/DMRm10t60_n260_RCV_Hcor cfl=0.4 AVtype=3    \
 redirect=${data_dir}/msg_DMRm10t60_n260_RCV_Hcor solver=4 &
wait

########## TEMP #########
exit
########## TEMP #########

echo "DOUBLE MACH REFLECTION: RUN HIGH RES"
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_Hyb_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_Hyb_av00 solver=3 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_Hyb_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_Hyb_av10 solver=3 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_av00 solver=4 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_av10 solver=4 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RPV_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_RPV_av00 solver=5 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RPV_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_RPV_av10 solver=5 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_FVS_av00 cfl=0.4 artvisc=0   \
 redirect=${data_dir}/msg_DMRm10t60_n520_FVS_av00 solver=6 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_FVS_av10 cfl=0.4 artvisc=0.1 \
 redirect=${data_dir}/msg_DMRm10t60_n520_FVS_av10 solver=6 &
./pion-ugs IC_DMRm10t60_n520.silo 5 1 outfile=${data_dir}/DMRm10t60_n520_RCV_Hcor cfl=0.4 AVtype=3    \
 redirect=${data_dir}/msg_DMRm10t60_n520_RCV_Hcor solver=4 &
wait

exit
