#!/bin/bash
#
# 2009-12-23 JM: makes 1D and 2D plots for shock tube problems.
# 2009-12-25 JM: Added variable number as last arg to plot2D_ST_MHDAW.sh
# 2010.12.07 JM: Grabs the correct file (without having to specify a
#  specific timestep) for each plot.
# 2010.12.28 JM: extra plots added for the various solvers.
#
# call with ./make_ST_plots.sh $test_dir $code_dir $data_dir
#
#
#test_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/test/problems/test_ShockTubes
#code_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/bin_serial
#data_dir=/users/jmackey/active/projects/uniform_grid_code/trunk/results/STMHD
test_dir=$1
code_dir=$2
data_dir=$3
ref_dir=$data_dir


cd $test_dir

#
# 1D plots
#
file1=`ls ${data_dir}/Toro1.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Toro10k1.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro1 "
./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro1

file1=`ls ${data_dir}/Toro2.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Toro10k2.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro2"
./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro2

file1=`ls ${data_dir}/Toro3.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Toro10k3.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro3"
./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro3

file1=`ls ${data_dir}/Toro4.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Toro10k4.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro4 "
./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro4

file1=`ls ${data_dir}/Toro5.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Toro10k5.0*[1-9]*.txt`
echo "./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro5 "
./plot_ST.sh $file1 $file2  ${test_dir}/STfig_Toro5

file1=`ls ${data_dir}/BrioWu.0*[1-9]*.txt`
file2=`ls ${ref_dir}/BrioWu10k.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_BrioWu"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_BrioWu

file1=`ls ${data_dir}/FalleFS.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Falle10kFS.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleFS"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleFS

file1=`ls ${data_dir}/FalleSS.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Falle10kSS.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleSS"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleSS

file1=`ls ${data_dir}/FalleFR.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Falle10kFR.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleFR"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleFR

file1=`ls ${data_dir}/FalleSR.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Falle10kSR.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleSR"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleSR

file1=`ls ${data_dir}/FalleOFS.0*[1-9]*.txt`
file2=`ls ${ref_dir}/Falle10kOFS.0*[1-9]*.txt`
echo "./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleOFS"
./plot_ST_MHD.sh $file1 $file2  ${test_dir}/STfig_FalleOFS

file1=`ls ${data_dir}/FalleAW.0*[1-9]*.txt`
file2=`ls ${data_dir}/Falle10kAW.0*[1-9]*.txt`
file3=`ls ${ref_dir}/Falle10kAW.0000*.txt`
echo "./plot_ST_MHDAW.sh $file1 $file3 $file2 ${test_dir}/STfig_FalleAW"
./plot_ST_MHDAW.sh $file1 $file3 $file2 ${test_dir}/STfig_FalleAW

#
# 2D plots
# Offsets are calculated from the interface at x0 via xoff = x0*(cos(t)-1) + 0.5*sin(t)
# Look in shock_tube.cc and the twoD2oneD.cc conversion program to see why.
#
file2=`ls ${data_dir}/Toro10k1.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D1_HYB_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T1_HYB 0.1919351
file2=`ls ${data_dir}/Toro10k2.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D2_HYB_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T2_HYB 0.1708204
file2=`ls ${data_dir}/Toro10k3.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D3_HYB_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T3_HYB 0.1708204
file2=`ls ${data_dir}/Toro10k4.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D4_HYB_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T4_HYB 0.1813777
file2=`ls ${data_dir}/Toro10k5.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D5_HYB_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T5_HYB 0.1391486
#
file2=`ls ${data_dir}/Toro10k1.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D1_RCV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T1_RCV 0.1919351
file2=`ls ${data_dir}/Toro10k2.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D2_RCV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T2_RCV 0.1708204
file2=`ls ${data_dir}/Toro10k3.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D3_RCV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T3_RCV 0.1708204
file2=`ls ${data_dir}/Toro10k4.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D4_RCV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T4_RCV 0.1813777
file2=`ls ${data_dir}/Toro10k5.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D5_RCV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T5_RCV 0.1391486
#
file2=`ls ${data_dir}/Toro10k1.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D1_RPV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T1_RPV 0.1919351
file2=`ls ${data_dir}/Toro10k2.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D2_RPV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T2_RPV 0.1708204
file2=`ls ${data_dir}/Toro10k3.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D3_RPV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T3_RPV 0.1708204
file2=`ls ${data_dir}/Toro10k4.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D4_RPV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T4_RPV 0.1813777
file2=`ls ${data_dir}/Toro10k5.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D5_RPV_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T5_RPV 0.1391486
#
file2=`ls ${data_dir}/Toro10k1.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D1_FVS_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T1_FVS 0.1919351
file2=`ls ${data_dir}/Toro10k2.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D2_FVS_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T2_FVS 0.1708204
file2=`ls ${data_dir}/Toro10k3.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D3_FVS_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T3_FVS 0.1708204
file2=`ls ${data_dir}/Toro10k4.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D4_FVS_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T4_FVS 0.1813777
file2=`ls ${data_dir}/Toro10k5.0*[1-9]*.txt`
./plot2D_ST_hydro.sh $file2  ${test_dir}/2D/Toro2D5_FVS_plotdata.txt ${test_dir}/2D/STfig_Toro2D_T5_FVS 0.1391486
#


file2=`ls ${data_dir}/BrioWu10k.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/BrioWu2D_FKJ_plotdata.txt   ${test_dir}/2D/STfig_BrioWu2D_FKJ 0.1708204
file2=`ls ${data_dir}/Falle10kFS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FS_FKJ_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FS_FKJ 0.1919351
file2=`ls ${data_dir}/Falle10kFR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FR_FKJ_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FR_FKJ 0.1708204
file2=`ls ${data_dir}/Falle10kSS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SS_FKJ_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SS_FKJ 0.1919351
file2=`ls ${data_dir}/Falle10kSR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SR_FKJ_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SR_FKJ 0.1708204
#
file2=`ls ${data_dir}/BrioWu10k.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/BrioWu2D_RCV_plotdata.txt   ${test_dir}/2D/STfig_BrioWu2D_RCV 0.1708204
file2=`ls ${data_dir}/Falle10kFS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FS_RCV_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FS_RCV 0.1919351
file2=`ls ${data_dir}/Falle10kFR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FR_RCV_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FR_RCV 0.1708204
file2=`ls ${data_dir}/Falle10kSS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SS_RCV_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SS_RCV 0.1919351
file2=`ls ${data_dir}/Falle10kSR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SR_RCV_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SR_RCV 0.1708204
#
file2=`ls ${data_dir}/BrioWu10k.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/BrioWu2D_HCR_plotdata.txt   ${test_dir}/2D/STfig_BrioWu2D_HCR 0.1708204
file2=`ls ${data_dir}/Falle10kFS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FS_HCR_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FS_HCR 0.1919351
file2=`ls ${data_dir}/Falle10kFR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_FR_HCR_plotdata.txt ${test_dir}/2D/STfig_Falle2D_FR_HCR 0.1708204
file2=`ls ${data_dir}/Falle10kSS.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SS_HCR_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SS_HCR 0.1919351
file2=`ls ${data_dir}/Falle10kSR.0*[1-9]*.txt`
./plot2D_ST_MHD.sh $file2 ${test_dir}/2D/Falle2D_SR_HCR_plotdata.txt ${test_dir}/2D/STfig_Falle2D_SR_HCR 0.1708204
#

./plot2D_ST_MHDAW.sh ${test_dir}/2D/AW2D_n016_FKJ_plotdata.txt ${test_dir}/2D/AW2D_n032_FKJ_plotdata.txt \
                     ${test_dir}/2D/AW2D_n064_FKJ_plotdata.txt ${test_dir}/2D/AW2D_n128_FKJ_plotdata.txt ${test_dir}/2D/STfig_AW2D_res_FKJ 9
./plot2D_ST_MHDAW.sh ${test_dir}/2D/AW2D_n016_RCV_plotdata.txt ${test_dir}/2D/AW2D_n032_RCV_plotdata.txt \
                     ${test_dir}/2D/AW2D_n064_RCV_plotdata.txt ${test_dir}/2D/AW2D_n128_RCV_plotdata.txt ${test_dir}/2D/STfig_AW2D_res_RCV 9
./plot2D_ST_MHDAW.sh ${test_dir}/2D/AW2D_n016_HCR_plotdata.txt ${test_dir}/2D/AW2D_n032_HCR_plotdata.txt \
                     ${test_dir}/2D/AW2D_n064_HCR_plotdata.txt ${test_dir}/2D/AW2D_n128_HCR_plotdata.txt ${test_dir}/2D/STfig_AW2D_res_HCR 9

#echo "FINISHED MAKING PLOTS.  OPENING A VIEWER NOW."
echo "COMPARE IMAGES TO http://www.astro.uni-bonn.de/~jmackey/jmac/node9.html"
#eog STfig_*.jpeg &
#eog 2D/STfig_*.jpeg &

mv STfig_* 2D/STfig_* ../  

exit


#./plot_ST.sh ${data_dir}/Toro1.000155.txt ${data_dir}/Toro10k1.007850.txt  ${test_dir}/Toro1 
#./plot_ST.sh ${data_dir}/Toro2.000118.txt ${data_dir}/Toro10k2.005890.txt  ${test_dir}/Toro2
#./plot_ST.sh ${data_dir}/Toro3.000186.txt ${data_dir}/Toro10k3.009148.txt  ${test_dir}/Toro3 
#./plot_ST.sh ${data_dir}/Toro4.000300.txt ${data_dir}/Toro10k4.014985.txt ${test_dir}/Toro4 
#./plot_ST.sh ${data_dir}/Toro5.000196.txt ${data_dir}/Toro10k5.009774.txt  ${test_dir}/Toro5 
#./plot_ST_MHD.sh ${data_dir}/BrioWu.000131.txt ${data_dir}/BrioWu10k.006521.txt ${test_dir}/BrioWu
#./plot_ST_MHD.sh ${data_dir}/FalleFS.000823.txt ${data_dir}/Falle10kFS.041121.txt ${test_dir}/FalleFS
#./plot_ST_MHD.sh ${data_dir}/FalleFR.000162.txt ${data_dir}/Falle10kFR.007985.txt ${test_dir}/FalleFR
#./plot_ST_MHD.sh ${data_dir}/FalleSS.000255.txt ${data_dir}/Falle10kSS.012667.txt ${test_dir}/FalleSS
#./plot_ST_MHD.sh ${data_dir}/FalleSR.000503.txt ${data_dir}/Falle10kSR.024477.txt ${test_dir}/FalleSR
#./plot_ST_MHD.sh ${data_dir}/FalleOFS.000368.txt ${data_dir}/Falle10kOFS.018371.txt ${test_dir}/FalleOFS
#./plot_ST_MHDAW.sh ${data_dir}/FalleAW.002539.txt ${data_dir}/Falle10kAW.000000.txt ${data_dir}/Falle10kAW.126473.txt ${test_dir}/FalleAW
#
##
## 2D plots
## Offsets are calculated from the interface at x0 via xoff = x0*(cos(t)-1) + 0.5*sin(t)
## Look in shock_tube.cc and the twoD2oneD.cc conversion program to see why.
##
#./plot2D_ST_hydro.sh ${data_dir}/Toro10k1.007850.txt  ${test_dir}/2D/Toro2D1_plotdata.txt ${test_dir}/2D/Toro2D_T1 0.1919351
#./plot2D_ST_hydro.sh ${data_dir}/Toro10k2.005890.txt  ${test_dir}/2D/Toro2D2_plotdata.txt ${test_dir}/2D/Toro2D_T2 0.1708204
#./plot2D_ST_hydro.sh ${data_dir}/Toro10k3.009148.txt  ${test_dir}/2D/Toro2D3_plotdata.txt ${test_dir}/2D/Toro2D_T3 0.1708204
#./plot2D_ST_hydro.sh ${data_dir}/Toro10k4.014985.txt  ${test_dir}/2D/Toro2D4_plotdata.txt ${test_dir}/2D/Toro2D_T4 0.1813777
#./plot2D_ST_hydro.sh ${data_dir}/Toro10k5.009774.txt  ${test_dir}/2D/Toro2D5_plotdata.txt ${test_dir}/2D/Toro2D_T5 0.1391486
##
#./plot2D_ST_MHD.sh ${data_dir}/BrioWu10k.006521.txt  ${test_dir}/2D/BrioWu2D_plotdata.txt ${test_dir}/2D/BrioWu2D 0.1708204
#./plot2D_ST_MHD.sh ${data_dir}/Falle10kFS.041121.txt ${test_dir}/2D/Falle2D_FS_plotdata.txt ${test_dir}/2D/Falle2D_FS 0.1919351
#./plot2D_ST_MHD.sh ${data_dir}/Falle10kFR.007985.txt  ${test_dir}/2D/Falle2D_FR_plotdata.txt ${test_dir}/2D/Falle2D_FR 0.1708204
#./plot2D_ST_MHD.sh ${data_dir}/Falle10kSS.012667.txt  ${test_dir}/2D/Falle2D_SS_plotdata.txt ${test_dir}/2D/Falle2D_SS 0.1919351
#./plot2D_ST_MHD.sh ${data_dir}/Falle10kSR.024477.txt  ${test_dir}/2D/Falle2D_SR_plotdata.txt ${test_dir}/2D/Falle2D_SR 0.1708204
##
#./plot2D_ST_MHDAW.sh ${test_dir}/2D/AW2D_n016_plotdata.txt ${test_dir}/2D/AW2D_n032_plotdata.txt \
#                     ${test_dir}/2D/AW2D_n064_plotdata.txt ${test_dir}/2D/AW2D_n128_plotdata.txt ${test_dir}/2D/AW2D_res 9
#
#echo "FINISHED MAKING PLOTS.  OPENING A VIEWER NOW."
#echo "COMPARE IMAGES TO http://www.astro.uni-bonn.de/~jmackey/jmac/node9.html"
#eog *.jpeg &
#
