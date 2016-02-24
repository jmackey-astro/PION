#!/bin/bash
#
# Field Loop Test
#
# -JM 2009-12-15 Added appending of figs.
# -JM 2009-12-17 Changed directories to be read in from args.
#


# call with ./run_FL_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}
code_dir=${2}
data_dir=${3}/FieldLoop
visit_cmd=${4}

# Just in case it doesn't exist, create the destination directory.
mkdir -p $data_dir
#rm ${data_dir}/*

cd ${code_dir}
echo "MAKE IN" $code_dir
./compile_code.sh
echo "MAKE SUCEEDED"
cd -

../../icgen_serial ${test_dir}/params_FieldLoop.txt       silo redirect=tmp_
../../icgen_serial ${test_dir}/params_FieldLoopVz.txt     silo redirect=tmp_
../../icgen_serial ${test_dir}/params_FieldLoopStatic.txt silo redirect=tmp_

../../pion_serial IC_FieldLoop200static.silo 5 1 outfile=${data_dir}/FieldLoop200static_Roe \
 redirect=${data_dir}/msg_FieldLoop200_Roe_static_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop200.silo 5 1 outfile=${data_dir}/FieldLoop200_Roe \
 redirect=${data_dir}/msg_FieldLoop200_Roe_planar_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
wait
../../pion_serial IC_FieldLoop200vz.silo 5 1 outfile=${data_dir}/FieldLoop200vz_Roe \
 redirect=${data_dir}/msg_FieldLoop200_Roe_VxVyVz_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop200vz.silo 5 1 outfile=${data_dir}/FieldLoop200vz_RoeAV0 \
 redirect=${data_dir}/msg_FieldLoop200_Roe_VxVyVz_av0 cfl=0.4 AVtype=0 EtaVisc=0.0 solver=4 &
wait

../../pion_serial IC_FieldLoop200static.silo 5 1 outfile=${data_dir}/FieldLoop200static_Lin \
 redirect=${data_dir}/msg_FieldLoop200_Lin_static_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop200.silo 5 1 outfile=${data_dir}/FieldLoop200_Lin \
 redirect=${data_dir}/msg_FieldLoop200_Lin_planar_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
wait
../../pion_serial IC_FieldLoop200vz.silo 5 1 outfile=${data_dir}/FieldLoop200vz_Lin \
 redirect=${data_dir}/msg_FieldLoop200_Lin_VxVyVz_av1 cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop200vz.silo 5 1 outfile=${data_dir}/FieldLoop200vz_LinAV0 \
 redirect=${data_dir}/msg_FieldLoop200_Lin_VxVyVz_av0 cfl=0.4 AVtype=0 EtaVisc=0.0 solver=4 &
wait

../../pion_serial IC_FieldLoop200static.silo 5 1 outfile=${data_dir}/FieldLoop200static_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop200_Roe_static_Hcorr cfl=0.4 AVtype=3 solver=4 &
../../pion_serial IC_FieldLoop200.silo 5 1 outfile=${data_dir}/FieldLoop200_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop200_Roe_planar_Hcorr cfl=0.4 AVtype=3 solver=4 &
../../pion_serial IC_FieldLoop200vz.silo 5 1 outfile=${data_dir}/FieldLoop200vz_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop200_Roe_VxVyVz_Hcorr cfl=0.4 AVtype=3 solver=4 &
wait

exit
echo "Sims finished, moving to FieldLoop dir to make figures."
cd ${test_dir}
./make_FL_img.sh ${data_dir} FieldLoop200static_Roe FieldLoop200static_Roe $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200_Roe FieldLoop200_Roe             $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200vz_Roe FieldLoop200vz_Roe         $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200vz_RoeAV0 FieldLoop200vz_RoeAV0   $visit_cmd

./make_FL_img.sh ${data_dir} FieldLoop200static_Lin FieldLoop200static_Lin $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200_Lin FieldLoop200_Lin             $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200vz_Lin FieldLoop200vz_Lin         $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200vz_LinAV0 FieldLoop200vz_LinAV0   $visit_cmd

./make_FL_img.sh ${data_dir} FieldLoop200static_Roe_Hcorr FieldLoop200static_Roe_Hcorr $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200_Roe_Hcorr       FieldLoop200_Roe_Hcorr       $visit_cmd
./make_FL_img.sh ${data_dir} FieldLoop200vz_Roe_Hcorr     FieldLoop200vz_Roe_Hcorr     $visit_cmd

convert ${data_dir}/FieldLoop200_Roe_MagP_00.tif ${data_dir}/FieldLoop200static_Roe_MagP_02.tif \
    ${data_dir}/FieldLoop200_Roe_MagP_02.tif +append ${test_dir}/FL200_Roe_MagP.tif
convert ${data_dir}/FieldLoop200_Lin_MagP_00.tif ${data_dir}/FieldLoop200static_Lin_MagP_02.tif \
    ${data_dir}/FieldLoop200_Lin_MagP_02.tif +append ${test_dir}/FL200_Lin_MagP.tif
convert ${data_dir}/FieldLoop200_Roe_Hcorr_MagP_00.tif ${data_dir}/FieldLoop200static_Roe_Hcorr_MagP_02.tif \
    ${data_dir}/FieldLoop200_Roe_Hcorr_MagP_02.tif +append ${test_dir}/FL200_Roe_Hcorr_MagP.tif

convert ${data_dir}/FieldLoop200_Roe_CurlB2D_00.tif ${data_dir}/FieldLoop200static_Roe_CurlB2D_02.tif \
    ${data_dir}/FieldLoop200_Roe_CurlB2D_02.tif +append ${test_dir}/FL200_Roe_CurlB2D.tif
convert ${data_dir}/FieldLoop200_Lin_CurlB2D_00.tif ${data_dir}/FieldLoop200static_Lin_CurlB2D_02.tif \
    ${data_dir}/FieldLoop200_Lin_CurlB2D_02.tif +append ${test_dir}/FL200_Lin_CurlB2D.tif
convert ${data_dir}/FieldLoop200_Roe_Hcorr_CurlB2D_00.tif ${data_dir}/FieldLoop200static_Roe_Hcorr_CurlB2D_02.tif \
    ${data_dir}/FieldLoop200_Roe_Hcorr_CurlB2D_02.tif +append ${test_dir}/FL200_Roe_Hcorr_CurlB2D.tif

convert ${test_dir}/FL200_Roe_MagP.tif    eps3:${test_dir}/FL200_Roe_MagP.eps
convert ${test_dir}/FL200_Roe_CurlB2D.tif eps3:${test_dir}/FL200_Roe_CurlB2D.eps
convert ${test_dir}/FL200_Lin_MagP.tif    eps3:${test_dir}/FL200_Lin_MagP.eps
convert ${test_dir}/FL200_Lin_CurlB2D.tif eps3:${test_dir}/FL200_Lin_CurlB2D.eps
convert ${test_dir}/FL200_Roe_Hcorr_MagP.tif    eps3:${test_dir}/FL200_Roe_Hcorr_MagP.eps
convert ${test_dir}/FL200_Roe_Hcorr_CurlB2D.tif eps3:${test_dir}/FL200_Roe_Hcorr_CurlB2D.eps

#convert ${test_dir}/FL200_Roe_MagP.tif ${test_dir}/FL200_Roe_MagP.png
#convert ${test_dir}/FL200_Roe_MagP.png eps2:${test_dir}/FL200_Roe_MagP.eps
#convert ${test_dir}/FL200_Roe_CurlB2D.tif ${test_dir}/FL200_Roe_CurlB2D.png
#convert ${test_dir}/FL200_Roe_CurlB2D.png eps2:${test_dir}/FL200_Roe_CurlB2D.eps
#convert ${test_dir}/FL200_Lin_MagP.tif ${test_dir}/FL200_Lin_MagP.png
#convert ${test_dir}/FL200_Lin_MagP.png eps2:${test_dir}/FL200_Lin_MagP.eps
#convert ${test_dir}/FL200_Lin_CurlB2D.tif ${test_dir}/FL200_Lin_CurlB2D.png
#convert ${test_dir}/FL200_Lin_CurlB2D.png eps2:${test_dir}/FL200_Lin_CurlB2D.eps

convert ${test_dir}/FL200_Roe_MagP.tif ${test_dir}/FL200_Roe_MagP.jpeg
convert ${test_dir}/FL200_Roe_CurlB2D.tif ${test_dir}/FL200_Roe_CurlB2D.jpeg
convert ${test_dir}/FL200_Lin_MagP.tif ${test_dir}/FL200_Lin_MagP.jpeg
convert ${test_dir}/FL200_Lin_CurlB2D.tif ${test_dir}/FL200_Lin_CurlB2D.jpeg
convert ${test_dir}/FL200_Roe_Hcorr_MagP.tif ${test_dir}/FL200_Roe_Hcorr_MagP.jpeg
convert ${test_dir}/FL200_Roe_Hcorr_CurlB2D.tif ${test_dir}/FL200_Roe_Hcorr_CurlB2D.jpeg

#DESTDIR=~/active/projects/MHD_ETs_2010/docs/rmhd_pillars/figs
#cp ${test_dir}/FL200_Lin_MagP.jpeg    $DESTDIR/FieldLoop_MagP.jpeg
#cp ${test_dir}/FL200_Lin_CurlB2D.jpeg $DESTDIR/FieldLoop_J.jpeg
#cp ${test_dir}/FL200_Lin_MagP.eps    $DESTDIR/FieldLoop_MagP.eps
#cp ${test_dir}/FL200_Lin_CurlB2D.eps $DESTDIR/FieldLoop_J.eps

./magp_decay.sh ${data_dir} FieldLoop2D_n200
exit

