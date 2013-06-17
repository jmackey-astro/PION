#!/bin/bash


#data_dir=/Users/jmackey/data/StarBench
data_dir=/vol/aibn129/aibn129_1/jmackey/scratch/StarBench
movie_dir=.
#rm ${movie_dir}/*.mp4
#rm ${data_dir}/*.silo

NP=4
ICG="mpirun -np $NP ../../../icgen_parallel"
PION="mpirun -np $NP ../../../pion_parallel"

$ICG params_TremblinMixing_D2_N064_NH0p5.txt silo
$ICG params_TremblinMixing_D2_N064_nH5.txt silo


$PION IC_TremblinMixing_D2_N064_nH5_0000.silo 5 1 \
 outfile=${data_dir}/TremblinMixing_D2_N064_nH5_Native \
 redirect=${data_dir}/log_TremblinMixing_D2_N064_nH5_Native_

$PION IC_TremblinMixing_D2_N064_NH0p5_0000.silo 5 1 \
 outfile=${data_dir}/TremblinMixing_D2_N064_nH0p5_Native \
 redirect=${data_dir}/log_TremblinMixing_D2_N064_nH0p5_Native_


./make_dens_Hplus_plots.sh \
 ${data_dir} TremblinMixing_D2_N064_nH5_Native_0000 \
 ${data_dir} fig_TremblinMixing_D2_N064_nH5_NativeP \
 "n(H)=5 test" 0.0 4.0 0.0 4.0  0.0 4.0  0 0

./make_dens_Hplus_plots.sh \
 ${data_dir} TremblinMixing_D2_N064_nH0p5_Native_0000 \
 ${data_dir} fig_TremblinMixing_D2_N064_nH0p5_NativeP \
 "n(H)=5 test" 0.0 4.0 0.0 4.0  -1.0 3.0  0 0

ffmpeg -f image2 -r 10.0 \
 -i ${data_dir}/fig_TremblinMixing_D2_N064_nH0p5_Native_%03d.png \
 -sameq -s 938x768 ${movie_dir}/TremblinMixing_D2_N064_nH0p5_Native.mp4

ffmpeg -f image2 -r 10.0 \
 -i ${data_dir}/fig_TremblinMixing_D2_N064_nH5_Native_%03d.png \
 -sameq -s 938x768 ${movie_dir}/TremblinMixing_D2_N064_nH5_Native.mp4


