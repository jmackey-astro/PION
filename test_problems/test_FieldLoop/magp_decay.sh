#!/bin/bash

#
# This file should generate a gnu.plt file which gnuplot can use to read data and generate a figure.
#
# -JM 2009-12-16 File created.  currently not run by default for the test problems.
# -JM 2010.12.06 Changed so that it reads the test runs, and runs by default.
#

#infile=./FL_LimVA_oa2av10_Lin_info.txt
input_dir=$1
outfile=$2

cat << EOF  > gnu.plt_VA_FKJ
set terminal postscript enhanced color eps
set output "${outfile}_VA_FKJ.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" offset 0.8,0.0
set xrange [0:2]
set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop200_Lin_static_av1info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Lin_planar_av1info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Lin_VxVyVz_av1info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Lin_VxVyVz_av0info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF

cat << EOF  > gnu.plt_VA_ROE
set terminal postscript enhanced color eps
set output "${outfile}_VA_ROE.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" offset 0.8,0.0
set xrange [0:2]
set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop200_Roe_static_av1info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_planar_av1info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_VxVyVz_av1info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_VxVyVz_av0info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF

#
# Once I code the H-correction for MHD I can also make this plot for
# H-corrected fluxes.  Although I doubt it will make much difference
# since the main effect seems to come from the slope limiting.
#
cat << EOF  > gnu.plt_VA_ROE_HCORR
set terminal postscript enhanced color eps
set output "${outfile}_VA_ROE_Hcorr.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" 0.8
set xrange [0:2]
set yrange [0.8:1.0]
set key top right
plot '${input_dir}/msg_FieldLoop200_Roe_static_Hcorr_info.txt' u 1:2 w l title 'Static {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_planar_Hcorr_info.txt' u 1:2 w l title 'Planar {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_VxVyVz_Hcorr_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.10', \
     '${input_dir}/msg_FieldLoop200_Roe_VxVyVz_av0_info.txt' u 1:2 w l title 'Vel.3D {/Symbol h}=0.00'
EOF


gnuplot gnu.plt_VA_FKJ
gnuplot gnu.plt_VA_ROE
gnuplot gnu.plt_VA_ROE_HCORR

convert -density 300 -quality 100 ${outfile}_VA_FKJ.eps ${outfile}_VA_FKJ.jpeg
convert -density 300 -quality 100 ${outfile}_VA_ROE.eps ${outfile}_VA_ROE.jpeg
convert -density 300 -quality 100 ${outfile}_VA_ROE_HCORR.eps ${outfile}_VA_ROE_HCORR.jpeg


exit


cat << EOF  > gnu.plt_VA_FKJ
set terminal postscript enhanced color eps
set output "${outfile}_VA_FKJ.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" 0.8
set xrange [0:2]
set yrange [0.8:1.0]
set key top right
plot './FL_LimVA_oa2av10_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.10', \
     './FL_LimVA_oa2av01_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.01', \
     './FL_LimVA_oa2av00_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.00'
EOF

cat << EOF  > gnu.plt_VA_ROE
set terminal postscript enhanced color eps
set output "${outfile}_VA_ROE.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" 0.8
set xrange [0:2]
set yrange [0.8:1.0]
set key top right
plot './FL_LimVA_oa2av10_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.10', \
     './FL_LimVA_oa2av01_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.01', \
     './FL_LimVA_oa2av00_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.00'
EOF

cat << EOF  > gnu.plt_MM_FKJ
set terminal postscript enhanced color eps
set output "${outfile}_MM_FKJ.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" 0.8
set xrange [0:2]
set yrange [0.65:1.0]
set key top right
plot './FL_LimMM_oa2av10_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.10', \
     './FL_LimMM_oa2av01_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.01', \
     './FL_LimMM_oa2av00_Lin_info.txt' u 1:2 w l title '{/Symbol h}=0.00'
EOF

cat << EOF  > gnu.plt_MM_ROE
set terminal postscript enhanced color eps
set output "${outfile}_MM_ROE.eps"
set size 0.5, 0.5
set xlabel "Time"
set ylabel "Magnetic Pressure" 0.8
set xrange [0:2]
set yrange [0.65:1.0]
set key top right
plot './FL_LimMM_oa2av10_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.10', \
     './FL_LimMM_oa2av01_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.01', \
     './FL_LimMM_oa2av00_Roe_info.txt' u 1:2 w l title '{/Symbol h}=0.00'
EOF

gnuplot gnu.plt_VA_FKJ
gnuplot gnu.plt_VA_ROE
#gnuplot gnu.plt_MM_FKJ
#gnuplot gnu.plt_MM_ROE

convert -density 300 -quality 100 ${outfile}_VA_FKJ.eps ${outfile}_VA_FKJ.jpeg
convert -density 300 -quality 100 ${outfile}_VA_ROE.eps ${outfile}_VA_ROE.jpeg
#convert -density 300 -quality 100 ${outfile}_MM_FKJ.eps ${outfile}_MM_FKJ.jpeg
#convert -density 300 -quality 100 ${outfile}_MM_ROE.eps ${outfile}_MM_ROE.jpeg

#cp FieldLoop_MM_FKJ.eps FieldLoop_VA_FKJ.eps   ~/active/jmthesis/trunk/figs/
#cp FieldLoop_MM_FKJ.jpeg FieldLoop_VA_FKJ.jpeg ~/active/jmthesis/trunk/figs/

exit
