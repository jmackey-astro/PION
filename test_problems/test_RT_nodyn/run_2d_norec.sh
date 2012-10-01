#!/bin/bash
#
# 2011.03.01 JM: changed from 32,256,100 cells to 32,256,100 cells.
#

test_dir=$1
exe_dir=$2
data_dir=$3

#cd ${code_dir}
# test for main_serial and icgen ???

ICGEN=${exe_dir}/icgen_serial
EXE=${exe_dir}/main_serial

###############################################################################
## 2D and 3D runs for stromgen spheres with no dynamics, testing the raytracer.
## 21/11/2009
###############################################################################
##############################################
## 2D runs with no dynamics, no recombinations
$ICGEN ${test_dir}/pf_rtt2d_n32_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n256_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n256_nh1_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n32_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh2_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n32_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n100_nh3_norec.txt silo redirect=msg_temp_
$ICGEN ${test_dir}/pf_rtt2d_n256_nh2_norec.txt silo redirect=msg_temp_

#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt100_ &
#$EXE ICTEST_rtt2D_n32_nh1_norec.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt500_ &
#wait
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt100_ &
#$EXE ICTEST_rtt2D_n32_nh2_norec.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt500_ &
#wait
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt010_ &
$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt100_ &
#$EXE ICTEST_rtt2D_n32_nh3_norec.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt500_ &
wait

#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt100_ &
#$EXE ICTEST_rtt2D_n100_nh1_norec.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt500_ &
#wait
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt100_ &
#$EXE ICTEST_rtt2D_n100_nh2_norec.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt500_ &
#wait
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt010_ &
$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt100_ &
#$EXE ICTEST_rtt2D_n100_nh3_norec.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt500_ &
wait


#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n256_nh1_norec.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt010_ &
$EXE ICTEST_rtt2D_n256_nh1_norec.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt100_ &
#$EXE ICTEST_rtt2D_n256_nh1_norec.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt500_ &
#wait
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n256_nh2_norec.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt010_ &
$EXE ICTEST_rtt2D_n256_nh2_norec.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt100_ &
#$EXE ICTEST_rtt2D_n256_nh2_norec.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt500_ &
#wait
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n256_nh3_norec.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt010_ &
$EXE ICTEST_rtt2D_n256_nh3_norec.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt100_ &
#$EXE ICTEST_rtt2D_n256_nh3_norec.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt500_ &
wait

####################################

echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
make -f Makefile.plotradius clean; make -j4 -f Makefile.plotradius
#
# nh=10 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt010  ${data_dir}/rtt2D_n32_nh1_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt100  ${data_dir}/rtt2D_n32_nh1_norec_dt100  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt500  ${data_dir}/rtt2D_n32_nh1_norec_dt500  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt010 ${data_dir}/rtt2D_n100_nh1_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt100 ${data_dir}/rtt2D_n100_nh1_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt500 ${data_dir}/rtt2D_n100_nh1_norec_dt500 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt010 ${data_dir}/rtt2D_n256_nh1_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt100 ${data_dir}/rtt2D_n256_nh1_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt500 ${data_dir}/rtt2D_n256_nh1_norec_dt500 0 10 5 silo

#
# nh=100 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt010  ${data_dir}/rtt2D_n32_nh2_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt100  ${data_dir}/rtt2D_n32_nh2_norec_dt100  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt500  ${data_dir}/rtt2D_n32_nh2_norec_dt500  0 10 5 silo
#
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt010 ${data_dir}/rtt2D_n100_nh2_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt100 ${data_dir}/rtt2D_n100_nh2_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt500 ${data_dir}/rtt2D_n100_nh2_norec_dt500 0 10 5 silo
#
./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt010 ${data_dir}/rtt2D_n256_nh2_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt100 ${data_dir}/rtt2D_n256_nh2_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt500 ${data_dir}/rtt2D_n256_nh2_norec_dt500 0 10 5 silo

#
# nh=1000 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt010  ${data_dir}/rtt2D_n32_nh3_norec_dt010  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt100  ${data_dir}/rtt2D_n32_nh3_norec_dt100  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt500  ${data_dir}/rtt2D_n32_nh3_norec_dt500  0 10 5 silo

./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt010 ${data_dir}/rtt2D_n100_nh3_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt100 ${data_dir}/rtt2D_n100_nh3_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt500 ${data_dir}/rtt2D_n100_nh3_norec_dt500 0 10 5 silo

./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt010 ${data_dir}/rtt2D_n256_nh3_norec_dt010 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt100 ${data_dir}/rtt2D_n256_nh3_norec_dt100 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt500 ${data_dir}/rtt2D_n256_nh3_norec_dt500 0 10 5 silo

./ssphere_norec_plots_2d.sh $test_dir $exe_dir $data_dir
cp photoncons* $data_dir
echo "all done with no-recombinations"
exit


