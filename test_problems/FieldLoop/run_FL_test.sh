#!/bin/bash
#
# Field Loop Test
#
# -JM 2009-12-15 Added appending of figs.
# -JM 2009-12-17 Changed directories to be read in from args.
# - 2016.03.06 JM: updated to work better...
#


# call with ./run_FL_tests.sh $test_dir $code_dir $data_dir
test_dir=${1}
code_dir=${2}
data_dir=${3}/FieldLoop
visit_cmd=${4}
resolution=${5}

# Just in case it doesn't exist, create the destination directory.
mkdir -p $data_dir


cd ${code_dir}
echo "MAKE IN " $code_dir
./compile_code.sh
echo "MAKE SUCEEDED"
cd -

../../icgen_serial ${test_dir}/params_FieldLoop${resolution}.txt       silo redirect=tmp_
../../icgen_serial ${test_dir}/params_FieldLoop${resolution}vz.txt     silo redirect=tmp_
../../icgen_serial ${test_dir}/params_FieldLoop${resolution}Static.txt silo redirect=tmp_

../../pion_serial IC_FieldLoop${resolution}static.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}static_Roe \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_static_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop${resolution}.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}_Roe \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_planar_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
#wait
../../pion_serial IC_FieldLoop${resolution}vz.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}vz_Roe \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop${resolution}vz.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}vz_RoeAV0 \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_av0_ cfl=0.4 AVtype=0 EtaVisc=0.0 solver=4 &
wait

../../pion_serial IC_FieldLoop${resolution}static.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}static_Lin \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Lin_static_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop${resolution}.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}_Lin \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Lin_planar_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
#wait
../../pion_serial IC_FieldLoop${resolution}vz.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}vz_Lin \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Lin_VxVyVz_av1_ cfl=0.4 AVtype=1 EtaVisc=0.1 solver=4 &
../../pion_serial IC_FieldLoop${resolution}vz.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}vz_LinAV0 \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Lin_VxVyVz_av0_ cfl=0.4 AVtype=0 EtaVisc=0.0 solver=4 &
wait

../../pion_serial IC_FieldLoop${resolution}static.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}static_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_static_Hcorr_ cfl=0.4 AVtype=3 solver=4 &
../../pion_serial IC_FieldLoop${resolution}.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_planar_Hcorr_ cfl=0.4 AVtype=3 solver=4 &
../../pion_serial IC_FieldLoop${resolution}vz.silo 5 1 outfile=${data_dir}/FieldLoop${resolution}vz_Roe_Hcorr \
 redirect=${data_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_Hcorr_ cfl=0.4 AVtype=3 solver=4 &
wait


#exit

echo "Sims finished, moving to FieldLoop dir to make figures."
cd ${test_dir}
./make_FL_img.sh ${data_dir} FieldLoop${resolution}static_Roe FieldLoop${resolution}static_Roe $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}_Roe FieldLoop${resolution}_Roe             $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}vz_Roe FieldLoop${resolution}vz_Roe         $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}vz_RoeAV0 FieldLoop${resolution}vz_RoeAV0   $visit_cmd
#sleep 10

./make_FL_img.sh ${data_dir} FieldLoop${resolution}static_Lin FieldLoop${resolution}static_Lin $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}_Lin FieldLoop${resolution}_Lin             $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}vz_Lin FieldLoop${resolution}vz_Lin         $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}vz_LinAV0 FieldLoop${resolution}vz_LinAV0   $visit_cmd
#sleep 10

./make_FL_img.sh ${data_dir} FieldLoop${resolution}static_Roe_Hcorr FieldLoop${resolution}static_Roe_Hcorr $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}_Roe_Hcorr       FieldLoop${resolution}_Roe_Hcorr       $visit_cmd
#sleep 10
./make_FL_img.sh ${data_dir} FieldLoop${resolution}vz_Roe_Hcorr     FieldLoop${resolution}vz_Roe_Hcorr     $visit_cmd
#sleep 10

convert ${data_dir}/FieldLoop${resolution}_Roe_MagP_00.png ${data_dir}/FieldLoop${resolution}static_Roe_MagP_02.png \
    ${data_dir}/FieldLoop${resolution}_Roe_MagP_02.png +append ${test_dir}/FL${resolution}_Roe_MagP.png
convert ${data_dir}/FieldLoop${resolution}_Lin_MagP_00.png ${data_dir}/FieldLoop${resolution}static_Lin_MagP_02.png \
    ${data_dir}/FieldLoop${resolution}_Lin_MagP_02.png +append ${test_dir}/FL${resolution}_Lin_MagP.png
convert ${data_dir}/FieldLoop${resolution}_Roe_Hcorr_MagP_00.png ${data_dir}/FieldLoop${resolution}static_Roe_Hcorr_MagP_02.png \
    ${data_dir}/FieldLoop${resolution}_Roe_Hcorr_MagP_02.png +append ${test_dir}/FL${resolution}_Roe_Hcorr_MagP.png

convert ${data_dir}/FieldLoop${resolution}_Roe_CurlB2D_00.png ${data_dir}/FieldLoop${resolution}static_Roe_CurlB2D_02.png \
    ${data_dir}/FieldLoop${resolution}_Roe_CurlB2D_02.png +append ${test_dir}/FL${resolution}_Roe_CurlB2D.png
convert ${data_dir}/FieldLoop${resolution}_Lin_CurlB2D_00.png ${data_dir}/FieldLoop${resolution}static_Lin_CurlB2D_02.png \
    ${data_dir}/FieldLoop${resolution}_Lin_CurlB2D_02.png +append ${test_dir}/FL${resolution}_Lin_CurlB2D.png
convert ${data_dir}/FieldLoop${resolution}_Roe_Hcorr_CurlB2D_00.png ${data_dir}/FieldLoop${resolution}static_Roe_Hcorr_CurlB2D_02.png \
    ${data_dir}/FieldLoop${resolution}_Roe_Hcorr_CurlB2D_02.png +append ${test_dir}/FL${resolution}_Roe_Hcorr_CurlB2D.png


./magp_decay.sh ${data_dir} FieldLoop2D_n${resolution}
exit

