#!/bin/bash

#
# This file should generate a gnu.plt file which gnuplot can use to read data and generate a figure.
#
# -JM 2009-12-16 File created.  currently not run by default for the test problems.
# -JM 2010.12.06 Changed so that it reads the test runs, and runs by default.
# - 2016.03.06 JM: updated to work for different resolutions.
#

input_dir=$1
outfile=$2
resolution=$3

cat << EOF  > gnu.plt_VA_FKJ
set terminal postscript enhanced color eps
set output "${input_dir}/${outfile}_VA_FKJ.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" offset 0.8,0.0
set xrange [0:2]
#set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop${resolution}_Lin_static_av1_info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Lin_planar_av1_info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Lin_VxVyVz_av1_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Lin_VxVyVz_av0_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF


cat << EOF  > gnu.plt_VA_ROE
set terminal postscript enhanced color eps
set output "${input_dir}/${outfile}_VA_ROE.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" offset 0.8,0.0
set xrange [0:2]
#set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop${resolution}_Roe_static_av1_info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_planar_av1_info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_av1_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_av0_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF


cat << EOF  > gnu.plt_VA_ROE_HCORR
set terminal postscript enhanced color eps
set output "${input_dir}/${outfile}_VA_ROE_Hcorr.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" offset 0.8,0.0
set xrange [0:2]
#set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop${resolution}_Roe_static_Hcorr_info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_planar_Hcorr_info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_Hcorr_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop${resolution}_Roe_VxVyVz_av0_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF


gnuplot gnu.plt_VA_FKJ
gnuplot gnu.plt_VA_ROE
gnuplot gnu.plt_VA_ROE_HCORR

convert -density 400 ${input_dir}/${outfile}_VA_FKJ.eps -background white -flatten ${input_dir}/${outfile}_VA_FKJ.png
convert -density 400 ${input_dir}/${outfile}_VA_ROE.eps -background white -flatten ${input_dir}/${outfile}_VA_ROE.png
convert -density 400 ${input_dir}/${outfile}_VA_ROE_Hcorr.eps -background white -flatten ${input_dir}/${outfile}_VA_ROE_Hcorr.png

exit


