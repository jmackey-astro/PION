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

convert ${data_dir}/BWaxi2D_HalfPlane_n128*_01.png +append ${data_dir}/N128_01.png
convert ${data_dir}/BWaxi2D_HalfPlane_n128*_02.png +append ${data_dir}/N128_02.png
convert ${data_dir}/N128_01.png ${data_dir}/N128_02.png -append N128_EarlyLate.png

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

convert ${data_dir}/BWaxi2D_HalfPlane_n064*_01.png +append ${data_dir}/N064_01.png
convert ${data_dir}/BWaxi2D_HalfPlane_n064*_02.png +append ${data_dir}/N064_02.png
convert ${data_dir}/N064_01.png ${data_dir}/N064_02.png -append N064_EarlyLate.png

mv ${data_dir}/N064*.png ${data_dir}/N128*.png ${test_dir}/../
rm *.png

