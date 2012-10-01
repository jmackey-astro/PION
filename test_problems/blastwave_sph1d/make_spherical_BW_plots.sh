#!/bin/bash
#
# Make plots of the results for spherical blastwave radii as a function of time.
#


# call with ./run_spherical_BW.sh $test_dir $code_dir $data_dir
data_dir=$1   # should be sub-directory 'blastwave' of the test-results directory.
#data_dir=/vol/aibn129/aibn129_1/jmackey/data_etc/code_tests/blastwave

file1=`ls ${data_dir}/BW1D_phys_n128_Roe_Hcorr.0*   | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_Roe_Hcorr.0*   | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n512_Roe_Hcorr.0*   | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_Roe_HCorr \
$file1 $file2 $file3 \
"RCV, Hcorr, N=128" "RCV, Hcorr, N=256" "RCV, Hcorr, N=512"

file1=`ls ${data_dir}/BW1D_phys_n128_Roe_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_Roe_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n512_Roe_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_Roe_FKJ98 \
$file1 $file2 $file3 \
"RCV, FKJ98, N=128" "RCV, FKJ98, N=256" "RCV, FKJ98, N=512"

file1=`ls ${data_dir}/BW1D_phys_n128_RPV_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_RPV_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n512_RPV_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_RPV_FKJ98 \
$file1 $file2 $file3 \
"RPV, FKJ98, N=128" "RPV, FKJ98, N=256" "RPV, FKJ98, N=512"

file1=`ls ${data_dir}/BW1D_phys_n128_FVS_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_FVS_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n512_FVS_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_FVS_FKJ98 \
$file1 $file2 $file3 \
"FVS, FKJ98, N=128" "FVS, FKJ98, N=256" "FVS, FKJ98, N=512"

file1=`ls ${data_dir}/BW1D_phys_n128_Hyb_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_Hyb_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n512_Hyb_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_Hyb_FKJ98 \
$file1 $file2 $file3 \
"Hyb, FKJ98, N=128" "Hyb, FKJ98, N=256" "Hyb, FKJ98, N=512"

#
# Now compare different solvers at N=256
#
file1=`ls ${data_dir}/BW1D_phys_n256_Hyb_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_Roe_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n256_Roe_Hcorr.0*   | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_aaN256_HRc1Rc2 \
$file1 $file2 $file3 \
"Hyb, FKJ98, N=256" "RCV, FKJ98, N=256" "RCV, H-corr, N=256"

file1=`ls ${data_dir}/BW1D_phys_n256_Hyb_FKJav01.0* | tail -n1`
file2=`ls ${data_dir}/BW1D_phys_n256_RPV_FKJav01.0* | tail -n1`
file3=`ls ${data_dir}/BW1D_phys_n256_FVS_FKJav01.0* | tail -n1`
./plot_sphBW_dens_pres_profile.sh BW1D_Profile50kyr_aaN256_HRpFS \
$file1 $file2 $file3 \
"Hyb, FKJ98, N=256" "RPV, FKJ98, N=256" "FVS, FKJ98, N=256"


#exit

./plot_sphBW_radius.sh BW1D_radius_Roe_HCorr \
${data_dir}/msg_BW1D_phys_n128_Roe_Hcorrinfo.txt \
${data_dir}/msg_BW1D_phys_n256_Roe_Hcorrinfo.txt \
${data_dir}/msg_BW1D_phys_n512_Roe_Hcorrinfo.txt \
"RCV, Hcorr, N=128" "RCV, Hcorr, N=256" "RCV, Hcorr, N=512"

./plot_sphBW_radius.sh BW1D_radius_RCV_FKJ98 \
${data_dir}/msg_BW1D_phys_n128_Roe_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n256_Roe_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n512_Roe_FKJav01info.txt \
"RCV, FKJ98, N=128" "RCV, FKJ98, N=256" "RCV, FKJ98, N=512"

./plot_sphBW_radius.sh BW1D_radius_FVS_FKJ98 \
${data_dir}/msg_BW1D_phys_n128_FVS_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n256_FVS_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n512_FVS_FKJav01info.txt \
"FVS, FKJ98, N=128" "FVS, FKJ98, N=256" "FVS, FKJ98, N=512"

./plot_sphBW_radius.sh BW1D_radius_Hyb_FKJ98 \
${data_dir}/msg_BW1D_phys_n128_Hyb_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n256_Hyb_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n512_Hyb_FKJav01info.txt \
"Hyb, FKJ98, N=128" "Hyb, FKJ98, N=256" "Hyb, FKJ98, N=512"

./plot_sphBW_radius.sh BW1D_radius_RPV_FKJ98 \
${data_dir}/msg_BW1D_phys_n128_RPV_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n256_RPV_FKJav01info.txt \
${data_dir}/msg_BW1D_phys_n512_RPV_FKJav01info.txt \
"RPV, FKJ98, N=128" "RPV, FKJ98, N=256" "RPV, FKJ98, N=512"

exit
