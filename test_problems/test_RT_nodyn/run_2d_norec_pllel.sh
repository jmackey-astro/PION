#!/bin/bash
#
# 2011.03.01 JM: changed from 32,256,100 cells to 32,256,100 cells.
# 2011.10.24 JM: New parallel RT test, to make sure results agree w/serial.

test_dir=$1
exe_dir=$2
data_dir=$3
serialdir=$4

#cd ${code_dir}

ICGEN="mpiexec -np 4 ${exe_dir}/icgen_parallel"
EXE="mpiexec -np 4 ${exe_dir}/gridcode_parallel"

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
$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt010 
$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt100 
#$EXE ICTEST_rtt2D_n32_nh1_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh1_norec_dt500 
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt010
$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt100 
#$EXE ICTEST_rtt2D_n32_nh2_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh2_norec_dt500 
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.3267 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt010 
$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.03267 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt100 
#$EXE ICTEST_rtt2D_n32_nh3_norec_0000.silo 5 1 cfl=0.006534 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n32_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n32_nh3_norec_dt500 


#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt010
$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt100
#$EXE ICTEST_rtt2D_n100_nh1_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh1_norec_dt500
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt010
$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt100
#$EXE ICTEST_rtt2D_n100_nh2_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh2_norec_dt500
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=1.0 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt010
$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=0.1 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt100
#$EXE ICTEST_rtt2D_n100_nh3_norec_0000.silo 5 1 cfl=0.02 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n100_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n100_nh3_norec_dt500



#
# nh=10 per cc
#
$EXE ICTEST_rtt2D_n256_nh1_norec_0000.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt010 
$EXE ICTEST_rtt2D_n256_nh1_norec_0000.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt100 
#$EXE ICTEST_rtt2D_n256_nh1_norec_0000.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh1_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh1_norec_dt500 
#
# nh=100 per cc
#
$EXE ICTEST_rtt2D_n256_nh2_norec_0000.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt010 
$EXE ICTEST_rtt2D_n256_nh2_norec_0000.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt100 
#$EXE ICTEST_rtt2D_n256_nh2_norec_0000.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh2_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh2_norec_dt500 
#
# nh=1000 per cc
#
$EXE ICTEST_rtt2D_n256_nh3_norec_0000.silo 5 1 cfl=2.54455 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt010 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt010
$EXE ICTEST_rtt2D_n256_nh3_norec_0000.silo 5 1 cfl=0.254455 opfreq=2 \
  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt100 \
  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt100 
#$EXE ICTEST_rtt2D_n256_nh3_norec_0000.silo 5 1 cfl=0.0508911 opfreq=10 \
#  checkpt_freq=100000 outfile=${data_dir}/rtt2D_n256_nh3_norec_dt500 \
#  redirect=${data_dir}/msg_rtt2D_n256_nh3_norec_dt500

####################################

echo moving on to analysis
cd ${test_dir}
echo Now in directory: ${test_dir}
pwd
make -f Makefile.plotradius clean; make -j4 -f Makefile.plotradius
#
# nh=10 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt010  ${data_dir}/rtt2D_n32_nh1_norec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt100  ${data_dir}/rtt2D_n32_nh1_norec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh1_norec_dt500  ${data_dir}/rtt2D_n32_nh1_norec_dt500_0000  0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt010 ${data_dir}/rtt2D_n100_nh1_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt100 ${data_dir}/rtt2D_n100_nh1_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh1_norec_dt500 ${data_dir}/rtt2D_n100_nh1_norec_dt500_0000 0 10 5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt010 ${data_dir}/rtt2D_n256_nh1_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt100 ${data_dir}/rtt2D_n256_nh1_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh1_norec_dt500 ${data_dir}/rtt2D_n256_nh1_norec_dt500_0000 0 10 5 silo

#
# nh=100 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt010  ${data_dir}/rtt2D_n32_nh2_norec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt100  ${data_dir}/rtt2D_n32_nh2_norec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh2_norec_dt500  ${data_dir}/rtt2D_n32_nh2_norec_dt500_0000  0 10 5 silo
#
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt010 ${data_dir}/rtt2D_n100_nh2_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt100 ${data_dir}/rtt2D_n100_nh2_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh2_norec_dt500 ${data_dir}/rtt2D_n100_nh2_norec_dt500_0000 0 10 5 silo
#
./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt010 ${data_dir}/rtt2D_n256_nh2_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt100 ${data_dir}/rtt2D_n256_nh2_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh2_norec_dt500 ${data_dir}/rtt2D_n256_nh2_norec_dt500_0000 0 10 5 silo

#
# nh=1000 per cc
#
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt010  ${data_dir}/rtt2D_n32_nh3_norec_dt010_0000  0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt100  ${data_dir}/rtt2D_n32_nh3_norec_dt100_0000  0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n32_nh3_norec_dt500  ${data_dir}/rtt2D_n32_nh3_norec_dt500_0000  0 10 5 silo

./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt010 ${data_dir}/rtt2D_n100_nh3_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt100 ${data_dir}/rtt2D_n100_nh3_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n100_nh3_norec_dt500 ${data_dir}/rtt2D_n100_nh3_norec_dt500_0000 0 10 5 silo

./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt010 ${data_dir}/rtt2D_n256_nh3_norec_dt010_0000 0 1  5 silo
./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt100 ${data_dir}/rtt2D_n256_nh3_norec_dt100_0000 0 2  5 silo
#./plot_radius ${data_dir}/rtt2D_n256_nh3_norec_dt500 ${data_dir}/rtt2D_n256_nh3_norec_dt500_0000 0 10 5 silo

./ssphere_norec_plots_2d.sh $test_dir $exe_dir $data_dir
cp photoncons* $data_dir
echo "------------ Finished with plot_radius, moving to silcocompare -------------"

#########################
### silo file compare ###
#########################
echo "############ Calculating DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"
cd ../../analysis/silocompare
make -j4 -f Makefile.silocompare
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh1_norec_dt010 rtt2D_n32_nh1_norec_dt010 cmp01 2
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh2_norec_dt010 rtt2D_n32_nh2_norec_dt010 cmp02 2
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh3_norec_dt010 rtt2D_n32_nh3_norec_dt010 cmp032 2
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh1_norec_dt100 rtt2D_n32_nh1_norec_dt100 cmp04 2
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh2_norec_dt100 rtt2D_n32_nh2_norec_dt100 cmp05 2
./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh3_norec_dt100 rtt2D_n32_nh3_norec_dt100 cmp06 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh1_norec_dt500 rtt2D_n32_nh1_norec_dt500 cmp07 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh2_norec_dt500 rtt2D_n32_nh2_norec_dt500 cmp08 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n32_nh3_norec_dt500 rtt2D_n32_nh3_norec_dt500 cmp09 2

./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh1_norec_dt010 rtt2D_n100_nh1_norec_dt010 cmp10 2
./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh2_norec_dt010 rtt2D_n100_nh2_norec_dt010 cmp11 2
./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh3_norec_dt010 rtt2D_n100_nh3_norec_dt010 cmp12 2
./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh1_norec_dt100 rtt2D_n100_nh1_norec_dt100 cmp13 2
./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh2_norec_dt100 rtt2D_n100_nh2_norec_dt100 cmp14 2
./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh3_norec_dt100 rtt2D_n100_nh3_norec_dt100 cmp15 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh1_norec_dt500 rtt2D_n100_nh1_norec_dt500 cmp16 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh2_norec_dt500 rtt2D_n100_nh2_norec_dt500 cmp17 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n100_nh3_norec_dt500 rtt2D_n100_nh3_norec_dt500 cmp18 2

./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh1_norec_dt010 rtt2D_n256_nh1_norec_dt010 cmp20 2
./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh2_norec_dt010 rtt2D_n256_nh2_norec_dt010 cmp21 2
./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh3_norec_dt010 rtt2D_n256_nh3_norec_dt010 cmp22 2
./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh1_norec_dt100 rtt2D_n256_nh1_norec_dt100 cmp23 2
./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh2_norec_dt100 rtt2D_n256_nh2_norec_dt100 cmp24 2
./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh3_norec_dt100 rtt2D_n256_nh3_norec_dt100 cmp25 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh1_norec_dt500 rtt2D_n256_nh1_norec_dt500 cmp26 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh2_norec_dt500 rtt2D_n256_nh2_norec_dt500 cmp27 2
#./silocompare ${data_dir} ${serialdir} rtt2D_n256_nh3_norec_dt500 rtt2D_n256_nh3_norec_dt500 cmp28 2

echo "############ PRINTING DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"
grep "^[0-5]" msg_cmp*fo.txt
rm msg_cmp*
echo "############ PRINTING DIFFS FROM ALL THE FILES -- SHOULD BE ALL ZEROS! ###########"

echo "all done with no-recombinations"
exit


