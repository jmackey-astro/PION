#!/bin/bash

# call with ./run_ST_tests.sh $test_dir $code_dir $data_dir
#test_dir=${1}/test_ShockTubes
#code_dir=$2
#data_dir=$3
#test_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_ShockTubes
#code_dir=/home/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/mnt/local/jm/temp_sims/code_test_dir
code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_ShockTubes
data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/temp_sims

# Just in case it doesn't exist, create the destination directory.
mkdir $data_dir

cd ${code_dir}
echo "MAKE IN" $code_dir
make -f Makefile.serial.code; make -f Makefile.serial.icgenerator
echo "MAKE SUCEEDED"

##########################
# 1D shock tube tests    #
# 200 cells              #
##########################
# HYDRO
./icgen ${test_dir}/pf_st_toro1.txt silo redirect=${data_dir}/ic_Toro1_
./icgen ${test_dir}/pf_st_toro2.txt silo redirect=${data_dir}/ic_Toro2_
./icgen ${test_dir}/pf_st_toro3.txt silo redirect=${data_dir}/ic_Toro3_
./icgen ${test_dir}/pf_st_toro4.txt silo redirect=${data_dir}/ic_Toro4_
./icgen ${test_dir}/pf_st_toro5.txt silo redirect=${data_dir}/ic_Toro5_


#
# Now run the shock-tube tests with the hybrid solver:
#
./main_serial IC_Toro1.silo 5 1 outfile=${data_dir}/Toro1 opfreq=0 redirect=${data_dir}/msg_Toro1_ artvisc=0.05 ooa=2 solver=3
./main_serial IC_Toro2.silo 5 1 outfile=${data_dir}/Toro2 opfreq=0 redirect=${data_dir}/msg_Toro2_ artvisc=0.05 ooa=2 solver=3
./main_serial IC_Toro3.silo 5 1 outfile=${data_dir}/Toro3 opfreq=0 redirect=${data_dir}/msg_Toro3_ artvisc=0.05 ooa=2 solver=3
./main_serial IC_Toro4.silo 5 1 outfile=${data_dir}/Toro4 opfreq=0 redirect=${data_dir}/msg_Toro4_ artvisc=0.05 ooa=2 solver=3
./main_serial IC_Toro5.silo 5 1 outfile=${data_dir}/Toro5 opfreq=0 redirect=${data_dir}/msg_Toro5_ artvisc=0.05 ooa=2 solver=3
#
# Now run the shock-tube tests with the FVS solver:
#
./main_serial IC_Toro1.silo 5 1 outfile=${data_dir}/Toro1_FVS opfreq=0 redirect=${data_dir}/msg_Toro1_FVS artvisc=0.05 ooa=2 solver=6
./main_serial IC_Toro2.silo 5 1 outfile=${data_dir}/Toro2_FVS opfreq=0 redirect=${data_dir}/msg_Toro2_FVS artvisc=0.05 ooa=2 solver=6
./main_serial IC_Toro3.silo 5 1 outfile=${data_dir}/Toro3_FVS opfreq=0 redirect=${data_dir}/msg_Toro3_FVS artvisc=0.05 ooa=2 solver=6
./main_serial IC_Toro4.silo 5 1 outfile=${data_dir}/Toro4_FVS opfreq=0 redirect=${data_dir}/msg_Toro4_FVS artvisc=0.05 ooa=2 solver=6
./main_serial IC_Toro5.silo 5 1 outfile=${data_dir}/Toro5_FVS opfreq=0 redirect=${data_dir}/msg_Toro5_FVS artvisc=0.05 ooa=2 solver=6

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
#
# Now run the shock-tube tests
#
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_FVSav0 solver=6 \
 redirect=${data_dir}/msg_Toro2D1_FVSav0_ artvisc=0 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_FVSav0 solver=6 \
 redirect=${data_dir}/msg_Toro2D2_FVSav0_ artvisc=0 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_FVSav0 solver=6 \
 redirect=${data_dir}/msg_Toro2D3_FVSav0_ artvisc=0 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_FVSav0 solver=6 \
 redirect=${data_dir}/msg_Toro2D4_FVSav0_ artvisc=0 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_FVSav0 solver=6 \
 redirect=${data_dir}/msg_Toro2D5_FVSav0_ artvisc=0 ooa=2 opfreq_time=5.0 &
wait
exit
#
# Now run the shock-tube tests
#
./main_serial IC_Toro2D1.silo 5 1 outfile=${data_dir}/Toro2D1_FVS solver=6 \
 redirect=${data_dir}/msg_Toro2D1_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D2.silo 5 1 outfile=${data_dir}/Toro2D2_FVS solver=6 \
 redirect=${data_dir}/msg_Toro2D2_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D3.silo 5 1 outfile=${data_dir}/Toro2D3_FVS solver=6 \
 redirect=${data_dir}/msg_Toro2D3_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D4.silo 5 1 outfile=${data_dir}/Toro2D4_FVS solver=6 \
 redirect=${data_dir}/msg_Toro2D4_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
./main_serial IC_Toro2D5.silo 5 1 outfile=${data_dir}/Toro2D5_FVS solver=6 \
 redirect=${data_dir}/msg_Toro2D5_ artvisc=0.15 ooa=2 opfreq_time=5.0 &
wait
exit
