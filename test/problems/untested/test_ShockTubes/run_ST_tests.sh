#!/bin/bash

# call with ./run_ST_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}/test_ShockTubes
code_dir=$2
data_dir=$3
#test_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_ShockTubes
#code_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/bin_serial

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

bin_dir=${code_dir}/../bin
echo CDIR=$code_dir BINDIR=$bin_dir

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
if [ ! -f ../bin/main_serial ] || [ ! -f ../bin/icgen_serial ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi
# Reset the negative pressure correction switch
sed -i -e "s/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/" ../source/defines/functionality_flags.h
#
# copy the executables into the test directory.
cd $test_dir
rm icgen_serial main_serial
cp ${bin_dir}/icgen_serial $test_dir/icgen_serial
cp ${bin_dir}/main_serial $test_dir/main_serial

##########################
# 1D shock tube tests    #
# 200 cells              #
##########################
# HYDRO
./icgen_serial ${test_dir}/pf_st_toro1.txt silo redirect=${data_dir}/ic_Toro1_
./icgen_serial ${test_dir}/pf_st_toro2.txt silo redirect=${data_dir}/ic_Toro2_
./icgen_serial ${test_dir}/pf_st_toro3.txt silo redirect=${data_dir}/ic_Toro3_
./icgen_serial ${test_dir}/pf_st_toro4.txt silo redirect=${data_dir}/ic_Toro4_
./icgen_serial ${test_dir}/pf_st_toro5.txt silo redirect=${data_dir}/ic_Toro5_
# MHD
./icgen_serial ${test_dir}/pf_st_falle07.txt silo redirect=${data_dir}/ic_BrioWu_
./icgen_serial ${test_dir}/pf_st_falle08.txt silo redirect=${data_dir}/ic_Falle08_
./icgen_serial ${test_dir}/pf_st_falle09.txt silo redirect=${data_dir}/ic_Falle09_ 
./icgen_serial ${test_dir}/pf_st_falle10.txt silo redirect=${data_dir}/ic_Falle10_
./icgen_serial ${test_dir}/pf_st_falle11.txt silo redirect=${data_dir}/ic_Falle11_
./icgen_serial ${test_dir}/pf_st_falle12.txt silo redirect=${data_dir}/ic_Falle12_
./icgen_serial ${test_dir}/pf_st_falle13.txt silo redirect=${data_dir}/ic_Falle13_

#
# Now run the shock-tube tests
#
./main_serial IC_Toro1.silo 5 1 outfile=${data_dir}/Toro1 opfreq=0 redirect=${data_dir}/msg_Toro1_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro2.silo 5 1 outfile=${data_dir}/Toro2 opfreq=0 redirect=${data_dir}/msg_Toro2_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro3.silo 5 1 outfile=${data_dir}/Toro3 opfreq=0 redirect=${data_dir}/msg_Toro3_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro4.silo 5 1 outfile=${data_dir}/Toro4 opfreq=0 redirect=${data_dir}/msg_Toro4_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro5.silo 5 1 outfile=${data_dir}/Toro5 opfreq=0 redirect=${data_dir}/msg_Toro5_ artvisc=0.05 ooa=2 &
./main_serial IC_BrioWu.silo   5 1 outfile=${data_dir}/BrioWu   opfreq=0 redirect=${data_dir}/msg_BrioWu_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleFS.silo  5 1 outfile=${data_dir}/FalleFS  opfreq=0 redirect=${data_dir}/msg_FalleFS_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleSS.silo  5 1 outfile=${data_dir}/FalleSS  opfreq=0 redirect=${data_dir}/msg_FalleSS_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleFR.silo  5 1 outfile=${data_dir}/FalleFR  opfreq=0 redirect=${data_dir}/msg_FalleFR_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleSR.silo  5 1 outfile=${data_dir}/FalleSR  opfreq=0 redirect=${data_dir}/msg_FalleSR_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleOFS.silo 5 1 outfile=${data_dir}/FalleOFS opfreq=0 redirect=${data_dir}/msg_FalleOFS_ artvisc=0.05 ooa=2 &
./main_serial IC_FalleAW.silo  5 1 outfile=${data_dir}/FalleAW  opfreq=0 redirect=${data_dir}/msg_FalleAW_ artvisc=0.05 ooa=2 &
wait

exit
