# call with ./run_ST_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}/test_ShockTubes
code_dir=$2
data_dir=$3
#test_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_ShockTubes
#code_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/bin_serial

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
#
# switch off negative pressure correction based on temperature, since we are dimensionless.
#
DDD=`grep "^#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE" ../source/defines/functionality_flags.h`
if [ -n "$DDD" ]
then
  echo "resetting negative pressure flag."
  sed -i -e "s/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/" ../source/defines/functionality_flags.h
fi
# compile the code
./compile_code.sh
#
# Reset the negative pressure correction switch
#
sed -i -e "s/\/\/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/#define SET_NEGATIVE_PRESSURE_TO_FIXED_TEMPERATURE/" ../source/defines/functionality_flags.h
if [ ! -f ../bin/main_serial ] || [ ! -f ../bin/icgen_serial ]
then
  echo "Cannot compile code"
  exit
else
  echo "MAKE SUCEEDED"
fi

##########################
# 1D shock tube tests    #
# 200 cells              #
##########################
# HYDRO
./icgen ${test_dir}/pf_st_toro1_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_toro2_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_toro3_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_toro4_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_toro5_hires.txt silo redirect=${data_dir}/ichires_
# MHD
./icgen ${test_dir}/pf_st_falle07_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle08_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle09_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle10_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle11_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle12_hires.txt silo redirect=${data_dir}/ichires_
./icgen ${test_dir}/pf_st_falle13_hires.txt silo redirect=${data_dir}/ichires_

#
# Now run the shock-tube tests
#
rm ${data_dir}/Toro10k* ${data_dir}/BrioWu10k* ${data_dir}/Falle10k*
#
./main_serial IC_Toro10k1.silo 5 1 outfile=${data_dir}/Toro10k1 opfreq=0 redirect=${data_dir}/msg_Toro10k1_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro10k2.silo 5 1 outfile=${data_dir}/Toro10k2 opfreq=0 redirect=${data_dir}/msg_Toro10k2_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro10k3.silo 5 1 outfile=${data_dir}/Toro10k3 opfreq=0 redirect=${data_dir}/msg_Toro10k3_ artvisc=0.05 ooa=2 &
./main_serial IC_Toro10k4.silo 5 1 outfile=${data_dir}/Toro10k4 opfreq=0 redirect=${data_dir}/msg_Toro10k4_ artvisc=0.05 ooa=2 &
#wait
./main_serial IC_Toro10k5.silo 5 1 outfile=${data_dir}/Toro10k5 opfreq=0 redirect=${data_dir}/msg_Toro10k5_ artvisc=0.05 ooa=2 &
./main_serial IC_BrioWu10k.silo   5 1 outfile=${data_dir}/BrioWu10k   opfreq=0 redirect=${data_dir}/msg_BrioWu10k_ artvisc=0.05 ooa=2 &
./main_serial IC_Falle10kFS.silo  5 1 outfile=${data_dir}/Falle10kFS  opfreq=0 redirect=${data_dir}/msg_Falle10kFS_ artvisc=0.05 ooa=2 &
./main_serial IC_Falle10kSS.silo  5 1 outfile=${data_dir}/Falle10kSS  opfreq=0 redirect=${data_dir}/msg_Falle10kSS_ artvisc=0.05 ooa=2 &
#wait
./main_serial IC_Falle10kFR.silo  5 1 outfile=${data_dir}/Falle10kFR  opfreq=0 redirect=${data_dir}/msg_Falle10kFR_ artvisc=0.05 ooa=2 &
./main_serial IC_Falle10kSR.silo  5 1 outfile=${data_dir}/Falle10kSR  opfreq=0 redirect=${data_dir}/msg_Falle10kSR_ artvisc=0.05 ooa=2 &
./main_serial IC_Falle10kOFS.silo 5 1 outfile=${data_dir}/Falle10kOFS opfreq=0 redirect=${data_dir}/msg_Falle10kOFS_ artvisc=0.05 ooa=2 &
./main_serial IC_Falle10kAW.silo  5 1 outfile=${data_dir}/Falle10kAW  opfreq=0 redirect=${data_dir}/msg_Falle10kAW_ artvisc=0.05 ooa=2 &
wait

exit

