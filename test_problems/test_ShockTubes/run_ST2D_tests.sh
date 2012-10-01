#!/bin/bash

#
# 2009-12-25 JM: Added MHD Alfven Waves stuff.
# 2010.12.07 JM: Grabs the correct file (without having to specify a
#  specific timestep) for each plot.
# 2010.12.28 JM: Added calculations for different solvers. (fixed bug 29.12)
# 2012.08.07 JM: Moved analysis to new file (ST2D_analyse_tests.sh).

#
# Usage:
# call with ./run_ST2D_tests.sh $test_dir $code_dir $data_dir
#
test_dir=${1}/test_ShockTubes/2D
code_dir=$2
data_dir=$3

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir} 
echo "MAKE IN" $code_dir
# switch off negative pressure correction based on temperature, since we are dimensionless.
DDD=`grep "^#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE" ../source/defines/functionality_flags.h`
if [ -n "$DDD" ]
then
  echo "resetting negative pressure flag."
  sed -i -e "s/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/" ../source/defines/functionality_flags.h
fi
# compile the code
./compile_code.sh
# Reset the negative pressure correction switch
sed -i -e "s/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/" ../source/defines/functionality_flags.h
if [ ! -f ../bin/main_serial ] || [ ! -f ../bin/icgen_serial ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi

##########################
# 2D shock tube tests    #
# 400x400 cells          #
##########################
# HYDRO
./icgen ${test_dir}/pf_st2D_toro1.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_toro2.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_toro3.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_toro4.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_toro5.txt silo redirect=${data_dir}/ic2d_
# MHD
./icgen ${test_dir}/pf_st2D_falle07.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_falle08.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_falle09.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_falle10.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_falle11.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_falle12.txt silo redirect=${data_dir}/ic2d_
#
./icgen ${test_dir}/pf_st2D_AWn16.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_AWn32.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_AWn64.txt silo redirect=${data_dir}/ic2d_
./icgen ${test_dir}/pf_st2D_AWn128.txt silo redirect=${data_dir}/ic2d_


#
# Now run the shock-tube tests
#
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_HYB redirect=${data_dir}/msg_Toro2D1_HYB artvisc=0.1 ooa=2 opfreq_time=5.0 solver=3 &
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_RCV redirect=${data_dir}/msg_Toro2D1_RCV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_RPV redirect=${data_dir}/msg_Toro2D1_RPV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=5 &
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_FVS redirect=${data_dir}/msg_Toro2D1_FVS artvisc=0   ooa=2 opfreq_time=5.0 solver=6 &
#wait
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_HYB redirect=${data_dir}/msg_Toro2D2_HYB artvisc=0.1 ooa=2 opfreq_time=5.0 solver=3 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_RCV redirect=${data_dir}/msg_Toro2D2_RCV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_RPV redirect=${data_dir}/msg_Toro2D2_RPV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=5 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_FVS redirect=${data_dir}/msg_Toro2D2_FVS artvisc=0   ooa=2 opfreq_time=5.0 solver=6 &
#wait
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_HYB redirect=${data_dir}/msg_Toro2D3_HYB artvisc=0.1 ooa=2 opfreq_time=5.0 solver=3 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_RCV redirect=${data_dir}/msg_Toro2D3_RCV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_RPV redirect=${data_dir}/msg_Toro2D3_RPV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=5 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_FVS redirect=${data_dir}/msg_Toro2D3_FVS artvisc=0   ooa=2 opfreq_time=5.0 solver=6 &
#wait
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_HYB redirect=${data_dir}/msg_Toro2D4_HYB artvisc=0.1 ooa=2 opfreq_time=5.0 solver=3 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_RCV redirect=${data_dir}/msg_Toro2D4_RCV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_RPV redirect=${data_dir}/msg_Toro2D4_RPV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=5 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_FVS redirect=${data_dir}/msg_Toro2D4_FVS artvisc=0   ooa=2 opfreq_time=5.0 solver=6 &
wait

./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_HYB redirect=${data_dir}/msg_Toro2D5_HYB artvisc=0.1 ooa=2 opfreq_time=5.0 solver=3 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_RCV redirect=${data_dir}/msg_Toro2D5_RCV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_RPV redirect=${data_dir}/msg_Toro2D5_RPV artvisc=0.1 ooa=2 opfreq_time=5.0 solver=5 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_FVS redirect=${data_dir}/msg_Toro2D5_FVS artvisc=0   ooa=2 opfreq_time=5.0 solver=6 &
#wait

#./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1 redirect=${data_dir}/msg_Toro2D1_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
#./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2 redirect=${data_dir}/msg_Toro2D2_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
#./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3 redirect=${data_dir}/msg_Toro2D3_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
#./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4 redirect=${data_dir}/msg_Toro2D4_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
#./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5 redirect=${data_dir}/msg_Toro2D5_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
#wait

./main_serial IC_BrioWu2D.silo   5 1 outfile=${data_dir}/BrioWu2D_FKJ   redirect=${data_dir}/msg_BrioWu2D_FKJ_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=1 & 
./main_serial IC_BrioWu2D.silo   5 1 outfile=${data_dir}/BrioWu2D_RCV   redirect=${data_dir}/msg_BrioWu2D_RCV_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 & 
./main_serial IC_BrioWu2D.silo   5 1 outfile=${data_dir}/BrioWu2D_HCR   redirect=${data_dir}/msg_BrioWu2D_HCR_ AVtype=3    ooa=2 opfreq_time=5.0 solver=4 & 
#wait
./main_serial IC_Falle2D_FS.silo  5 1 outfile=${data_dir}/Falle2D_FS_FKJ  redirect=${data_dir}/msg_Falle2D_FS_FKJ_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_Falle2D_FS.silo  5 1 outfile=${data_dir}/Falle2D_FS_RCV  redirect=${data_dir}/msg_Falle2D_FS_RCV_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Falle2D_FS.silo  5 1 outfile=${data_dir}/Falle2D_FS_HCR  redirect=${data_dir}/msg_Falle2D_FS_HCR_ AVtype=3    ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_Falle2D_SS.silo  5 1 outfile=${data_dir}/Falle2D_SS_FKJ  redirect=${data_dir}/msg_Falle2D_SS_FKJ_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_Falle2D_SS.silo  5 1 outfile=${data_dir}/Falle2D_SS_RCV  redirect=${data_dir}/msg_Falle2D_SS_RCV_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Falle2D_SS.silo  5 1 outfile=${data_dir}/Falle2D_SS_HCR  redirect=${data_dir}/msg_Falle2D_SS_HCR_ AVtype=3    ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_Falle2D_FR.silo  5 1 outfile=${data_dir}/Falle2D_FR_FKJ  redirect=${data_dir}/msg_Falle2D_FR_FKJ_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_Falle2D_FR.silo  5 1 outfile=${data_dir}/Falle2D_FR_RCV  redirect=${data_dir}/msg_Falle2D_FR_RCV_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Falle2D_FR.silo  5 1 outfile=${data_dir}/Falle2D_FR_HCR  redirect=${data_dir}/msg_Falle2D_FR_HCR_ AVtype=3    ooa=2 opfreq_time=5.0 solver=4 &
wait
./main_serial IC_Falle2D_SR.silo  5 1 outfile=${data_dir}/Falle2D_SR_FKJ  redirect=${data_dir}/msg_Falle2D_SR_FKJ_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_Falle2D_SR.silo  5 1 outfile=${data_dir}/Falle2D_SR_RCV  redirect=${data_dir}/msg_Falle2D_SR_RCV_ artvisc=0.1 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Falle2D_SR.silo  5 1 outfile=${data_dir}/Falle2D_SR_HCR  redirect=${data_dir}/msg_Falle2D_SR_HCR_ AVtype=3    ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_Falle2D_AW.silo  5 1 outfile=${data_dir}/Falle2D_AW_FKJ  redirect=${data_dir}/msg_Falle2D_AW_FKJ_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_Falle2D_AW.silo  5 1 outfile=${data_dir}/Falle2D_AW_RCV  redirect=${data_dir}/msg_Falle2D_AW_RCV_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_Falle2D_AW.silo  5 1 outfile=${data_dir}/Falle2D_AW_HCR  redirect=${data_dir}/msg_Falle2D_AW_HCR_ AVtype=3     ooa=2 opfreq_time=5.0 solver=4 &
#wait

./main_serial IC_AW2D_n128.silo  5 1 outfile=${data_dir}/AW2D_n128_FKJ  redirect=${data_dir}/msg_AW2D_n128_FKJ_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_AW2D_n128.silo  5 1 outfile=${data_dir}/AW2D_n128_RCV  redirect=${data_dir}/msg_AW2D_n128_RCV_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_AW2D_n128.silo  5 1 outfile=${data_dir}/AW2D_n128_HCR  redirect=${data_dir}/msg_AW2D_n128_HCR_ AVtype=3     ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_AW2D_n64.silo  5 1 outfile=${data_dir}/AW2D_n064_FKJ  redirect=${data_dir}/msg_AW2D_n064_FKJ_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_AW2D_n64.silo  5 1 outfile=${data_dir}/AW2D_n064_RCV  redirect=${data_dir}/msg_AW2D_n064_RCV_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_AW2D_n64.silo  5 1 outfile=${data_dir}/AW2D_n064_HCR  redirect=${data_dir}/msg_AW2D_n064_HCR_ AVtype=3     ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_AW2D_n32.silo  5 1 outfile=${data_dir}/AW2D_n032_FKJ  redirect=${data_dir}/msg_AW2D_n032_FKJ_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_AW2D_n32.silo  5 1 outfile=${data_dir}/AW2D_n032_RCV  redirect=${data_dir}/msg_AW2D_n032_RCV_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_AW2D_n32.silo  5 1 outfile=${data_dir}/AW2D_n032_HCR  redirect=${data_dir}/msg_AW2D_n032_HCR_ AVtype=3     ooa=2 opfreq_time=5.0 solver=4 &
#wait
./main_serial IC_AW2D_n16.silo  5 1 outfile=${data_dir}/AW2D_n016_FKJ  redirect=${data_dir}/msg_AW2D_n016_FKJ_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=1 &
./main_serial IC_AW2D_n16.silo  5 1 outfile=${data_dir}/AW2D_n016_RCV  redirect=${data_dir}/msg_AW2D_n016_RCV_ artvisc=0.05 ooa=2 opfreq_time=5.0 solver=4 &
./main_serial IC_AW2D_n16.silo  5 1 outfile=${data_dir}/AW2D_n016_HCR  redirect=${data_dir}/msg_AW2D_n016_HCR_ AVtype=3     ooa=2 opfreq_time=5.0 solver=4 &
wait

exit
