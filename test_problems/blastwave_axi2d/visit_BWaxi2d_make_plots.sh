#!/bin/bash

#
# Make plots for all the axisymmetric blast wave tests.
#
test_dir=$1
code_dir=$2
data_dir=$3
visit_cmd=$4


infile_base=BWaxi2D_HalfPlane_n128_FVS_FKJav01
outfile_base=BWaxi2D_HalfPlane_n128_FVS_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n128_Hyb_FKJav01
outfile_base=BWaxi2D_HalfPlane_n128_Hyb_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n128_Roe_FKJav01
outfile_base=BWaxi2D_HalfPlane_n128_Roe_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n128_Roe_Hcorr
outfile_base=BWaxi2D_HalfPlane_n128_Roe_Hcorr
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n128_RPV_FKJav01
outfile_base=BWaxi2D_HalfPlane_n128_RPV_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#

convert ${test_dir}/BWaxi2D_HalfPlane_n128*_01.jpeg +append ${test_dir}/N128_01.jpeg
convert ${test_dir}/BWaxi2D_HalfPlane_n128*_02.jpeg +append ${test_dir}/N128_02.jpeg
convert ${test_dir}/N128_01.jpeg ${test_dir}/N128_02.jpeg -append N128_EarlyLate.jpeg

#
#
#

infile_base=BWaxi2D_HalfPlane_n064_FVS_FKJav01
outfile_base=BWaxi2D_HalfPlane_n064_FVS_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n064_Hyb_FKJav01
outfile_base=BWaxi2D_HalfPlane_n064_Hyb_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n064_Roe_FKJav01
outfile_base=BWaxi2D_HalfPlane_n064_Roe_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n064_Roe_Hcorr
outfile_base=BWaxi2D_HalfPlane_n064_Roe_Hcorr
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#
infile_base=BWaxi2D_HalfPlane_n064_RPV_FKJav01
outfile_base=BWaxi2D_HalfPlane_n064_RPV_FKJav01
./make_BWaxi2D_image.sh $test_dir $code_dir $data_dir $visit_cmd $infile_base $outfile_base
#

convert ${test_dir}/BWaxi2D_HalfPlane_n064*_01.jpeg +append ${test_dir}/N064_01.jpeg
convert ${test_dir}/BWaxi2D_HalfPlane_n064*_02.jpeg +append ${test_dir}/N064_02.jpeg
convert ${test_dir}/N064_01.jpeg ${test_dir}/N064_02.jpeg -append N064_EarlyLate.jpeg

mv N064*.jpeg N128*.jpeg ${test_dir}/../
rm *.jpeg

